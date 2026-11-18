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
///        to a plane in a uniform field, and the dense solves on
///        std::vector matrices (DenseSolve.h).

#ifndef GTRK_STUDY_DECAYFITMATH_H_
#define GTRK_STUDY_DECAYFITMATH_H_

#include "GenericTracking/Study/Dual.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/DenseSolve.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

namespace o2::gtrk::study::decay_fit_detail
{

using o2::gtrk::study::dual::Dual;

/// q/p's numerator for a particle of \a charge: |q|, 1 for a neutral one
/// (its q/p is 1/p).
[[nodiscard]] inline double numeratorOf(double charge) { return charge != 0. ? std::abs(charge) : 1.; }

/// A momentum: its vector and its magnitude [GeV].
template <typename T>
struct Momentum {
  std::array<T, 3> vector{};
  T magnitude{};
};

/// The momentum of direction \a phi, \a theta and \a qOverP with q/p's
/// \a numerator (`numeratorOf`), in doubles or dual numbers.
template <typename T>
[[nodiscard]] Momentum<T> momentumOf(const T& phi, const T& theta, const T& qOverP, double numerator)
{
  using dual::valueOf;
  using std::cos;
  using std::sin;
  const T p = T{valueOf(qOverP) > 0. ? numerator : -numerator} / qOverP;
  const T pt = p * sin(theta);
  return {{pt * cos(phi), pt * sin(phi), p * cos(theta)}, p};
}

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

/// \a from carried along its helix (\a fieldScale = Bz B2C, the curvature per
/// q/pT) to \a plane -- along the beam, across it, or tilted --: the free
/// state at the crossing, ahead or behind, and its loc0, loc1 there; false if
/// the helix does not reach the plane (along the beam: within a quarter turn).
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
    // Oriented so that the helix leaves through the plane's front (the
    // normal and the direction on the same side, cos(delta) > 0): the
    // crossing then has sin(delta + rho s) = sin(delta) - rho u, u the start's
    // distance before the plane, solved for w = sin(rho s) without the
    // cancellation at small rho (the line's -u / cos(delta) at rho = 0),
    // s = (w / rho) asinc(w). asin recovers rho s within a quarter turn only.
    const double phiL = std::atan2(ny, nx);
    const D delta = phi + (-phiL);
    const double side = cos(delta).v >= 0. ? 1. : -1.;
    const D u = side * ((nx * (from[0] + (-static_cast<double>(plane.origin[0])))) + (ny * (from[1] + (-static_cast<double>(plane.origin[1])))));
    const D sinDelta = side * sin(delta);
    const D cosDelta = side * cos(delta);
    const D sinEnd = sinDelta - (rho * u);
    if (!(std::abs(sinEnd.v) < 1.) || !(cosDelta.v > 0.)) {
      return false; // missed, grazing, or along the plane: no crossing to bind to
    }
    const D cosEnd = sqrt((-(sinEnd * sinEnd)) + 1.);
    if (!((cosEnd.v * cosDelta.v) + (sinEnd.v * sinDelta.v) > 0.)) {
      return false; // beyond a quarter turn
    }
    const D wOverRho = (-(u * cosDelta)) - ((sinDelta * u) * (sinEnd + sinDelta) / (cosDelta + cosEnd));
    const D w = rho * wOverRho;
    if (!(std::abs(w.v) < 1.)) {
      return false;
    }
    s = wOverRho * asinc(w);
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

/// A row-major matrix of \a columns columns in a std::vector, seen by the
/// dense solves (DenseSolve.h).
[[nodiscard]] inline dense::MatrixView<double> viewOf(std::vector<double>& a, int columns) { return {a.data(), columns}; }

/// A symmetric positive definite matrix's Cholesky factor.
struct Cholesky {
  int n = 0;
  std::vector<double> l;
  /// Factor \a a (\a size x \a size, row-major, its lower triangle read);
  /// false if it is not positive definite.
  [[nodiscard]] bool factor(const std::vector<double>& a, int size)
  {
    n = size;
    l = a;
    return dense::choleskyFactor(viewOf(l, n), n);
  }
  [[nodiscard]] double logDeterminant() { return dense::choleskyLogDeterminant(viewOf(l, n), n); }
  /// x <- a^-1 x.
  void solve(std::vector<double>& x) { dense::choleskySolve(viewOf(l, n), n, x.data()); }
  /// a^-1, n x n, row-major.
  [[nodiscard]] std::vector<double> inverse()
  {
    std::vector<double> out(static_cast<std::size_t>(n * n));
    dense::choleskyInverse(viewOf(l, n), n, viewOf(out, n));
    return out;
  }
};

/// A track state's 5 x 5 covariance in double, row-major.
template <typename State>
[[nodiscard]] std::vector<double> covarianceOf(const State& state)
{
  std::vector<double> covariance(25U);
  for (int r = 0; r < 5; ++r) {
    for (int c = 0; c < 5; ++c) {
      covariance[(r * 5) + c] = state.covAt(r, c);
    }
  }
  return covariance;
}

/// The weight (inverse covariance, 5 x 5, row-major) of the covariance
/// factored in \a cholesky.
inline void weightOf(Cholesky& cholesky, std::array<double, 25>& weight)
{
  dense::choleskyInverse(viewOf(cholesky.l, 5), 5, dense::MatrixView<double>{weight.data(), 5});
}

/// A positive semi-definite matrix's factor L, L L^T = \a a (n x n,
/// row-major): the pivoted L D L^T's P^T L D^1/2 (`dense::semidefiniteFactor`).
/// L is not triangular.
[[nodiscard]] inline std::vector<double> semidefiniteFactor(const std::vector<double>& a, int n)
{
  std::vector<double> work = a;
  std::vector<int> order(static_cast<std::size_t>(n));
  std::vector<double> factor(static_cast<std::size_t>(n * n));
  dense::semidefiniteFactor(viewOf(work, n), n, order.data(), viewOf(factor, n));
  return factor;
}

} // namespace o2::gtrk::study::decay_fit_detail

#endif // GTRK_STUDY_DECAYFITMATH_H_
