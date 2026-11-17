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

#include "GenericTracking/Tracking/TrackingSurface.h"

#include <algorithm>
#include <array>
#include <cmath>
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
/// sin(x) / x, 1 at 0.
template <int N>
Dual<N> sinc(const Dual<N>& a)
{
  const double x = a.v;
  if (std::abs(x) < 1e-4) {
    const double x2 = x * x;
    return chain(a, 1. - (x2 / 6.), -x / 3.);
  }
  return chain(a, std::sin(x) / x, (std::cos(x) / x) - (std::sin(x) / (x * x)));
}
/// asin(x) / x, 1 at 0.
template <int N>
Dual<N> asinc(const Dual<N>& a)
{
  const double x = a.v;
  if (std::abs(x) < 1e-4) {
    const double x2 = x * x;
    return chain(a, 1. + (x2 / 6.), x / 3.);
  }
  const double root = std::sqrt(1. - (x * x));
  return chain(a, std::asin(x) / x, (1. / (x * root)) - (std::asin(x) / (x * x)));
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

} // namespace o2::gtrk::study::decay_fit_detail

#endif // GTRK_STUDY_DECAYFITMATH_H_
