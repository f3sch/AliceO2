// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DecayFitMath.h
/// \brief The arithmetic the decay vertex fits share: the helix from a point
///        to a plane in a uniform field, and the dense solves (Eigen).

#ifndef GTRK_STUDY_DECAYFITMATH_H_
#define GTRK_STUDY_DECAYFITMATH_H_

#include "GenericTracking/Study/Dual.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

#include <Eigen/Dense>

namespace o2::gtrk::study::decay_fit_detail
{

using o2::gtrk::study::dual::Dual;

/// The most iterations of a fit's loops, and the most halvings of one step:
/// double's mantissa bits. A cap, not a derivation: a step halved that often
/// no longer moves a double; Gauss-Newton under a fixed model converges
/// quadratically near its solution, far faster; a fit that reaches the cap
/// is reported as not converged.
inline constexpr int kMaxSteps = std::numeric_limits<double>::digits;

/// A tolerance, stated openly: how far the material model may still move a
/// converged solution, in its chi2 metric (step^T N step), and the fit count
/// as converged -- one standard deviation. The model is evaluated along the
/// solution from a material table that is piecewise constant along a
/// trajectory, so an exact fixed point need not exist; a solution the model's
/// ambiguity moves by less than its own error is taken, more is a failure.
inline constexpr double kModelTolerance = 1.;

/// A free state -- x, y, z, phi, theta, q/p -- with its derivatives by N
/// inputs.
template <int N>
using FreeState = std::array<Dual<N>, 6>;

/// \a from carried along its helix to \a plane (a plane along the beam, or
/// across it): the free state at the crossing nearest along the helix,
/// ahead or behind, and its loc0, loc1 there; false if the helix never
/// reaches the plane.
///
/// On a plane along the beam with normal angle phiL the crossing has
/// sin(phi_end - phiL) = sin(phi - phiL) - rho u, u the start's distance
/// along the normal; of the two solutions phi_end and every turn by 2 pi, the
/// smallest |rho s| = |phi_end - phi| is the nearest crossing. In double the
/// difference loses nothing that matters at small rho.
template <int N>
[[nodiscard]] bool crossToPlane(const FreeState<N>& from, const TrackingSurface& plane, double fieldScale, FreeState<N>& to, std::array<Dual<N>, 2>& local)
{
  using D = Dual<N>;
  const double nx = plane.normal[0];
  const double ny = plane.normal[1];
  const double nz = plane.normal[2];
  const D& phi = from[3];
  const D& theta = from[4];
  const D sinTheta = sin(theta);
  const D cotTheta = cos(theta) / sinTheta;
  const D rho = fieldScale * (from[5] / sinTheta);
  D s;
  if (nz == 0.) {
    const double phiL = std::atan2(ny, nx);
    const D u = (nx * (from[0] + (-static_cast<double>(plane.origin[0])))) + (ny * (from[1] + (-static_cast<double>(plane.origin[1]))));
    const D delta = phi + (-phiL);
    if (rho.v == 0.) {
      const D cosDelta = cos(delta);
      if (!(std::abs(cosDelta.v) > 0.)) {
        return false;
      }
      s = (-u) / cosDelta;
    } else {
      const D sinEnd = sin(delta) - (rho * u);
      if (!(std::abs(sinEnd.v) < 1.)) {
        return false; // missed, or grazing: no crossing to bind to
      }
      const D first = chain(sinEnd, std::asin(sinEnd.v), 1. / std::sqrt(1. - (sinEnd.v * sinEnd.v)));
      // The turn to each solution, brought within (-pi, pi]; the smaller wins.
      const auto turnTo = [&](const D& end) {
        D turn = end - delta;
        turn.v = std::remainder(turn.v, 2. * std::numbers::pi);
        return turn;
      };
      const D one = turnTo(first);
      const D two = turnTo((-first) + std::numbers::pi);
      s = (std::abs(one.v) <= std::abs(two.v) ? one : two) / rho;
    }
  } else if (nx == 0. && ny == 0.) {
    if (!(std::abs(cotTheta.v) > 0.)) {
      return false;
    }
    s = (((-1.) * from[2]) + static_cast<double>(plane.origin[2])) / cotTheta;
  } else {
    // A tilted plane: f(s) = n . (x(s) - origin) = 0 by Newton's method from
    // the straight line's crossing, while each step shrinks; s's derivatives
    // by the implicit function, -df / f'(s).
    const auto f = [&](const D& at) {
      const D half = 0.5 * (rho * at);
      const D chord = at * sinc(half);
      return (nx * ((from[0] + (chord * cos(phi + half))) + (-static_cast<double>(plane.origin[0])))) + (ny * ((from[1] + (chord * sin(phi + half))) + (-static_cast<double>(plane.origin[1])))) + (nz * ((from[2] + (at * cotTheta)) + (-static_cast<double>(plane.origin[2]))));
    };
    const auto slope = [&](double at) { return (nx * std::cos(phi.v + (rho.v * at))) + (ny * std::sin(phi.v + (rho.v * at))) + (nz * cotTheta.v); };
    const double ahead = slope(0.);
    if (!(std::abs(ahead) > 0.)) {
      return false;
    }
    double at = -f(D{0.}).v / ahead;
    double last = std::numeric_limits<double>::infinity();
    for (int step = 0; step < kMaxSteps; ++step) {
      const double fs = slope(at);
      if (!(std::abs(fs) > 0.)) {
        return false;
      }
      const double change = f(D{at}).v / fs;
      if (!(std::abs(change) < last)) {
        break; // no longer shrinking: at the root to rounding
      }
      last = std::abs(change);
      at -= change;
    }
    const D value = f(D{at});
    // Closed to the plane's own precision (its geometry is floats), or no
    // crossing: f's terms' size times float epsilon.
    {
      const double half = 0.5 * rho.v * at;
      const double chord = at * (half == 0. ? 1. : std::sin(half) / half);
      const std::array<double, 3> x{from[0].v + (chord * std::cos(phi.v + half)), from[1].v + (chord * std::sin(phi.v + half)), from[2].v + (at * cotTheta.v)};
      double size = 0.;
      for (int c = 0; c < 3; ++c) {
        size += std::abs(static_cast<double>(plane.normal[c])) * (std::abs(x[c]) + std::abs(static_cast<double>(plane.origin[c])));
      }
      if (!(std::abs(value.v) <= static_cast<double>(std::numeric_limits<Scalar>::epsilon()) * size)) {
        return false;
      }
    }
    s = D{at};
    const double fs = slope(at);
    for (int k = 0; k < N; ++k) {
      s.d[k] = -value.d[k] / fs;
    }
  }
  const D half = 0.5 * (rho * s);
  const D chord = s * sinc(half);
  to[0] = from[0] + (chord * cos(phi + half));
  to[1] = from[1] + (chord * sin(phi + half));
  to[2] = from[2] + (s * cotTheta);
  to[3] = phi + (rho * s);
  to[4] = theta;
  to[5] = from[5];
  const std::array<D, 3> d{to[0] + (-static_cast<double>(plane.origin[0])), to[1] + (-static_cast<double>(plane.origin[1])), to[2] + (-static_cast<double>(plane.origin[2]))};
  local[0] = (static_cast<double>(plane.axis0[0]) * d[0]) + (static_cast<double>(plane.axis0[1]) * d[1]) + (static_cast<double>(plane.axis0[2]) * d[2]);
  local[1] = (static_cast<double>(plane.axis1[0]) * d[0]) + (static_cast<double>(plane.axis1[1]) * d[1]) + (static_cast<double>(plane.axis1[2]) * d[2]);
  return true;
}

/// The plane's frame: its transverse origin and its normal's angle.
struct PerigeePlane {
  double x = 0.;
  double y = 0.;
  double phi = 0.;
};

[[nodiscard]] inline PerigeePlane perigeePlaneFrame(const TrackingSurface& plane)
{
  return {.x = plane.origin[0], .y = plane.origin[1], .phi = std::atan2(static_cast<double>(plane.normal[1]), static_cast<double>(plane.normal[0]))};
}

/// The helix from (vx, vy, vz) with direction phi, theta and q/p in a
/// uniform field (\a fieldScale = Bz B2C, the curvature per q/pT) to the
/// perigee plane \a plane: its loc0, loc1, phi, theta, q/p there; false if it
/// does not reach the plane on its nearer half turn.
template <int N>
[[nodiscard]] bool predictOnPlane(const Dual<N>& vx, const Dual<N>& vy, const Dual<N>& vz, const Dual<N>& phi, const Dual<N>& theta, const Dual<N>& qOverP, double fieldScale, const PerigeePlane& plane, std::array<Dual<N>, 5>& out)
{
  using T = Dual<N>;
  const double cosL = std::cos(plane.phi);
  const double sinL = std::sin(plane.phi);
  const T dx = vx + (-plane.x);
  const T dy = vy + (-plane.y);
  const T sinTheta = sin(theta);
  const T rho = fieldScale * (qOverP / sinTheta); // dphi / ds, s the transverse arc length
  const T delta = phi + (-plane.phi);
  const T sinDelta = sin(delta);
  const T cosDelta = cos(delta);
  const T u = (cosL * dx) + (sinL * dy); // the start's distance before the plane
  // sin(delta + rho s) = sin(delta) - rho u at the crossing; s = asin(w) / rho,
  // w = sin(rho s), written without the cancellation at small rho.
  const T sinEnd = sinDelta - (rho * u);
  if (!(std::abs(sinEnd.v) < 1.) || !(cosDelta.v > 0.)) {
    return false;
  }
  const T rootEnd = sqrt((-(sinEnd * sinEnd)) + 1.);
  // The turn rho s, from the end's direction on the near side (cos > 0):
  // asin recovers it only within a quarter turn, cos(rho s) > 0.
  if (!((rootEnd.v * cosDelta.v) + (sinEnd.v * sinDelta.v) > 0.)) {
    return false;
  }
  const T wOverRho = (-(u * cosDelta)) - ((sinDelta * u) * (sinEnd + sinDelta) / (cosDelta + rootEnd));
  const T w = rho * wOverRho;
  if (!(std::abs(w.v) < 1.)) {
    return false;
  }
  const T s = wOverRho * asinc(w);
  const T half = 0.5 * (rho * s);
  out[0] = ((-sinL) * dx) + (cosL * dy) + (s * (sin(delta + half) * sinc(half)));
  out[1] = vz + (s * (cos(theta) / sinTheta));
  out[2] = phi + (rho * s);
  out[3] = theta;
  out[4] = qOverP;
  return true;
}

/// The dense solves' matrices: row-major n x n in a std::vector, seen by Eigen.
using RowMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
[[nodiscard]] inline Eigen::Map<RowMatrix> matrixOf(std::vector<double>& a, int rows, int columns) { return {a.data(), rows, columns}; }
[[nodiscard]] inline Eigen::Map<const RowMatrix> matrixOf(const std::vector<double>& a, int rows, int columns) { return {a.data(), rows, columns}; }
[[nodiscard]] inline Eigen::Map<Eigen::VectorXd> vectorOf(std::vector<double>& a) { return {a.data(), static_cast<Eigen::Index>(a.size())}; }
[[nodiscard]] inline Eigen::Map<const Eigen::VectorXd> vectorOf(const std::vector<double>& a) { return {a.data(), static_cast<Eigen::Index>(a.size())}; }

/// A symmetric positive definite matrix's Cholesky factorisation; false if
/// it is not positive definite.
[[nodiscard]] inline bool factorise(Eigen::LLT<RowMatrix>& llt, const RowMatrix& a)
{
  llt.compute(a);
  return llt.info() == Eigen::Success;
}

/// ln det of a matrix from its Cholesky factorisation.
[[nodiscard]] inline double logDeterminant(const Eigen::LLT<RowMatrix>& llt)
{
  return 2. * llt.matrixLLT().diagonal().array().log().sum();
}

/// The inverse of the matrix factorised in \a llt (n x n, row-major).
[[nodiscard]] inline std::vector<double> inverseOf(const Eigen::LLT<RowMatrix>& llt, int n)
{
  std::vector<double> inverse(static_cast<std::size_t>(n * n));
  matrixOf(inverse, n, n) = llt.solve(RowMatrix::Identity(n, n));
  return inverse;
}

/// A track state's 5 x 5 covariance in double.
template <typename State>
[[nodiscard]] RowMatrix covarianceOf(const State& state)
{
  RowMatrix covariance(5, 5);
  for (int r = 0; r < 5; ++r) {
    for (int c = 0; c < 5; ++c) {
      covariance(r, c) = state.covAt(r, c);
    }
  }
  return covariance;
}

/// The weight (inverse covariance, 5 x 5, row-major) of the covariance
/// factorised in \a llt.
inline void weightOf(const Eigen::LLT<RowMatrix>& llt, std::array<double, 25>& weight)
{
  Eigen::Map<Eigen::Matrix<double, 5, 5, Eigen::RowMajor>>{weight.data()} = llt.solve(RowMatrix::Identity(5, 5));
}

/// A positive semi-definite matrix's factor L, L L^T = \a a (n x n,
/// row-major), from the pivoted LDL^T: L = P^T L' D^1/2, a pivot that is not
/// positive -- the matrix singular there, or below zero by rounding --
/// leaving its column zero. L is not triangular.
[[nodiscard]] inline std::vector<double> semidefiniteFactor(const std::vector<double>& a, int n)
{
  const Eigen::LDLT<RowMatrix> ldlt(matrixOf(a, n, n));
  const Eigen::VectorXd root = ldlt.vectorD().cwiseMax(0.).cwiseSqrt();
  RowMatrix l = RowMatrix(ldlt.matrixL()) * root.asDiagonal();
  l = ldlt.transpositionsP().transpose() * l;
  std::vector<double> factor(static_cast<std::size_t>(n * n));
  matrixOf(factor, n, n) = l;
  return factor;
}

} // namespace o2::gtrk::study::decay_fit_detail

#endif // GTRK_STUDY_DECAYFITMATH_H_
