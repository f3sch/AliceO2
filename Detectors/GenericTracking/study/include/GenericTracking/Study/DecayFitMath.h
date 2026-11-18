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
/// \brief The arithmetic the decay vertex fits share: forward
///        differentiation (dual numbers), the helix from a point to a perigee
///        plane in a uniform field, and the dense symmetric solves.

#ifndef GTRK_STUDY_DECAYFITMATH_H_
#define GTRK_STUDY_DECAYFITMATH_H_

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

namespace o2::gtrk::study::decay_fit_detail
{

/// A value and its derivatives by N inputs: forward differentiation.
template <int N>
struct Dual {
  double v = 0.;
  std::array<double, N> d{};
};

template <int N>
Dual<N> operator+(const Dual<N>& a, const Dual<N>& b)
{
  Dual<N> r{a.v + b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a.d[i] + b.d[i];
  }
  return r;
}
template <int N>
Dual<N> operator-(const Dual<N>& a, const Dual<N>& b)
{
  Dual<N> r{a.v - b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a.d[i] - b.d[i];
  }
  return r;
}
template <int N>
Dual<N> operator-(const Dual<N>& a)
{
  Dual<N> r{-a.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = -a.d[i];
  }
  return r;
}
template <int N>
Dual<N> operator*(const Dual<N>& a, const Dual<N>& b)
{
  Dual<N> r{a.v * b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a.d[i] * b.v + a.v * b.d[i];
  }
  return r;
}
template <int N>
Dual<N> operator/(const Dual<N>& a, const Dual<N>& b)
{
  Dual<N> r{a.v / b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = (a.d[i] - r.v * b.d[i]) / b.v;
  }
  return r;
}
template <int N>
Dual<N> operator*(double a, const Dual<N>& b)
{
  Dual<N> r{a * b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a * b.d[i];
  }
  return r;
}
template <int N>
Dual<N> operator+(const Dual<N>& a, double b)
{
  Dual<N> r = a;
  r.v += b;
  return r;
}
/// f(a) with f'(a) given.
template <int N>
Dual<N> chain(const Dual<N>& a, double value, double derivative)
{
  Dual<N> r{value};
  for (int i = 0; i < N; ++i) {
    r.d[i] = derivative * a.d[i];
  }
  return r;
}
template <int N>
Dual<N> sin(const Dual<N>& a)
{
  return chain(a, std::sin(a.v), std::cos(a.v));
}
template <int N>
Dual<N> cos(const Dual<N>& a)
{
  return chain(a, std::cos(a.v), -std::sin(a.v));
}
template <int N>
Dual<N> sqrt(const Dual<N>& a)
{
  const double root = std::sqrt(a.v);
  return chain(a, root, 0.5 / root);
}
/// Below these arguments sinc and asinc take their series, through the x^6
/// term: there its first omitted term is below double's rounding, where the
/// closed forms' derivatives lose digits to cancellation. sinc's omitted term
/// is x^8 / 9!, asinc's 35 x^8 / 1152.
inline const double kSincSeries = std::pow(362880. * std::numeric_limits<double>::epsilon(), 1. / 8.);
inline const double kAsincSeries = std::pow((1152. / 35.) * std::numeric_limits<double>::epsilon(), 1. / 8.);

/// sin(x) / x, 1 at 0.
template <int N>
Dual<N> sinc(const Dual<N>& a)
{
  const double x = a.v;
  if (std::abs(x) < kSincSeries) {
    const double x2 = x * x;
    return chain(a, 1. - (x2 / 6.) + (x2 * x2 / 120.) - (x2 * x2 * x2 / 5040.), x * (-(1. / 3.) + (x2 / 30.) - (x2 * x2 / 840.)));
  }
  return chain(a, std::sin(x) / x, (std::cos(x) / x) - (std::sin(x) / (x * x)));
}
/// asin(x) / x, 1 at 0.
template <int N>
Dual<N> asinc(const Dual<N>& a)
{
  const double x = a.v;
  if (std::abs(x) < kAsincSeries) {
    const double x2 = x * x;
    return chain(a, 1. + (x2 / 6.) + (3. * x2 * x2 / 40.) + (5. * x2 * x2 * x2 / 112.), x * ((1. / 3.) + (3. * x2 / 10.) + (15. * x2 * x2 / 56.)));
  }
  const double root = std::sqrt(1. - (x * x));
  return chain(a, std::asin(x) / x, (1. / (x * root)) - (std::asin(x) / (x * x)));
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

/// A symmetric matrix's Cholesky factor in place (lower triangle, row-major
/// n x n); false if it is not positive definite.
[[nodiscard]] inline bool cholesky(std::vector<double>& a, int n)
{
  for (int j = 0; j < n; ++j) {
    double diagonal = a[j * n + j];
    for (int k = 0; k < j; ++k) {
      diagonal -= a[j * n + k] * a[j * n + k];
    }
    if (!(diagonal > 0.)) {
      return false;
    }
    const double root = std::sqrt(diagonal);
    a[j * n + j] = root;
    for (int i = j + 1; i < n; ++i) {
      double value = a[i * n + j];
      for (int k = 0; k < j; ++k) {
        value -= a[i * n + k] * a[j * n + k];
      }
      a[i * n + j] = value / root;
    }
  }
  return true;
}

/// Solve L L^T x = b with the factor of `cholesky`.
inline void choleskySolve(const std::vector<double>& factor, int n, std::span<double> b)
{
  for (int i = 0; i < n; ++i) {
    double value = b[i];
    for (int k = 0; k < i; ++k) {
      value -= factor[i * n + k] * b[k];
    }
    b[i] = value / factor[i * n + i];
  }
  for (int i = n - 1; i >= 0; --i) {
    double value = b[i];
    for (int k = i + 1; k < n; ++k) {
      value -= factor[k * n + i] * b[k];
    }
    b[i] = value / factor[i * n + i];
  }
}

/// The inverse of a symmetric positive definite matrix from its factor.
[[nodiscard]] inline std::vector<double> choleskyInverse(const std::vector<double>& factor, int n)
{
  std::vector<double> inverse(static_cast<std::size_t>(n * n), 0.);
  std::vector<double> column(static_cast<std::size_t>(n));
  for (int j = 0; j < n; ++j) {
    std::ranges::fill(column, 0.);
    column[j] = 1.;
    choleskySolve(factor, n, column);
    for (int i = 0; i < n; ++i) {
      inverse[i * n + j] = column[i];
    }
  }
  return inverse;
}

/// A positive semi-definite matrix's factor L, L L^T = \a a (n x n,
/// row-major; lower triangle): a pivot that is not positive -- the matrix
/// singular there, or below zero by rounding -- leaves its column zero.
[[nodiscard]] inline std::vector<double> semidefiniteFactor(const std::vector<double>& a, int n)
{
  std::vector<double> l(static_cast<std::size_t>(n * n), 0.);
  for (int j = 0; j < n; ++j) {
    double diagonal = a[j * n + j];
    for (int k = 0; k < j; ++k) {
      diagonal -= l[j * n + k] * l[j * n + k];
    }
    if (!(diagonal > 0.)) {
      continue;
    }
    const double root = std::sqrt(diagonal);
    l[j * n + j] = root;
    for (int i = j + 1; i < n; ++i) {
      double value = a[i * n + j];
      for (int k = 0; k < j; ++k) {
        value -= l[i * n + k] * l[j * n + k];
      }
      l[i * n + j] = value / root;
    }
  }
  return l;
}

} // namespace o2::gtrk::study::decay_fit_detail

#endif // GTRK_STUDY_DECAYFITMATH_H_
