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

/// \file DecayVertexFit.h
/// \brief One decay vertex fitted from its charged daughters' fitted tracks:
///        the vertex and every daughter's momentum at it.
///
/// The parameters are the vertex V and, per daughter, its direction and q/p
/// at V (phi, theta, q/p; Billoir's choice). A daughter enters as its fitted
/// state on its innermost measurement, which holds all its rows'
/// information, carried inward with its species' mass and the material
/// between (`carryInward`) to the perigee plane about the current vertex: the
/// plane through the vertex's transverse point along the beam, normal to the
/// daughter's direction at its closest approach (`perigeePlaneOf`). Its five
/// parameters there and their covariance are the measurement. The
/// prediction is the helix from V with the daughter's parameters to that
/// plane, exact in a uniform field: V and the plane are close, so no
/// material lies between. Its derivatives come from forward differentiation
/// of the same formula, exact to rounding.
///
/// Gauss-Newton on all parameters at once (one dense normal matrix, 3 + 3N
/// wide): every iteration carries the daughters anew to the plane about the
/// latest vertex, so the material and the transport follow the vertex. A
/// step that raises the chi2 is halved. The start is the closest approach of
/// the first two daughters' helices in the transverse plane.
///
/// The fit is in double: the states come in as floats.

#ifndef GTRK_STUDY_DECAYVERTEXFIT_H_
#define GTRK_STUDY_DECAYVERTEXFIT_H_

#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/HypothesisRefit.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/TrackFitterUtils.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <vector>

namespace o2::gtrk::study
{

/// How a decay vertex fit takes its daughters.
enum class DecayFitVariant : std::uint8_t {
  kFull,       ///< carried with their material to the plane about every iteration's vertex
  kCarryOnce,  ///< carried once, to the plane about the start; later iterations re-linearise only
  kNoMaterial, ///< carried every iteration, without material
  kCount
};
inline constexpr std::size_t kDecayFitVariants = static_cast<std::size_t>(DecayFitVariant::kCount);

/// Why a decay vertex fit stopped.
enum class DecayFitStatus : std::uint8_t {
  kOk,
  kNoStart,     ///< the first two daughters' helices give no start
  kCarryFailed, ///< a daughter could not be carried to its plane
  kNoCrossing,  ///< the prediction's helix does not reach a plane
  kSingular,    ///< the normal matrix or a measurement's covariance is not positive definite
  kNotConverged,
  kCount
};

/// A daughter as the fit takes it.
struct DecayDaughter {
  BoundTrackState innermost; ///< fitted, on `surface`
  TrackingSurface surface;   ///< its innermost measurement's plane
  Propagator propagator;     ///< the field, with the daughter's species' mass
  double mass = 0.;          ///< [GeV]
};

/// A decay vertex fit's result.
struct DecayVertexFitResult {
  DecayFitStatus status = DecayFitStatus::kNoStart;
  std::uint8_t iterations = 0U;
  std::array<double, 3> vertex{};
  /// Per daughter: phi, theta, q/p at the vertex.
  std::vector<std::array<double, 3>> momenta;
  /// The covariance of (vertex, the daughters' momenta in order), (3 + 3N)^2.
  std::vector<double> covariance;
  double chi2 = 0.;
  int ndf = 0;
  /// The daughters' last carried states' phi, theta, q/p: their momenta
  /// without the vertex's information.
  std::vector<std::array<double, 3>> carried;
  /// The covariance of each carried state's phi, theta, q/p, 3 x 3.
  std::vector<std::array<double, 9>> carriedCovariance;
  /// The material [x/X0] each daughter's last carry crossed.
  std::vector<double> carryX2X0;
};

namespace decay_fit_detail
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

/// The global position of \a state on \a surface [cm].
[[nodiscard]] inline std::array<double, 3> positionOf(const BoundTrackState& state, const TrackingSurface& surface)
{
  Scalar position[3]{};
  surfaceToGlobal(surface, state.loc0(), state.loc1(), position);
  return {position[0], position[1], position[2]};
}

/// The start: where the transverse circles of the first two daughters come
/// closest -- their crossing with the smaller z gap, or the point between
/// them on the line of their centres --, at the mean of their z there.
[[nodiscard]] inline std::optional<std::array<double, 3>> closestApproach(const DecayDaughter& a, const DecayDaughter& b)
{
  struct Circle {
    double cx, cy, radius, rho, px, py, pz, cotTheta;
  };
  const auto circleOf = [](const DecayDaughter& daughter) {
    const auto position = positionOf(daughter.innermost, daughter.surface);
    const double phi = daughter.innermost.phi();
    const double theta = daughter.innermost.theta();
    const double rho = daughter.propagator.curvature(daughter.innermost.qOverP(), std::sin(theta));
    return Circle{position[0] - (std::sin(phi) / rho), position[1] + (std::cos(phi) / rho), 1. / std::abs(rho), rho, position[0], position[1], position[2], std::cos(theta) / std::sin(theta)};
  };
  const auto one = circleOf(a);
  const auto two = circleOf(b);
  if (!std::isfinite(one.radius) || !std::isfinite(two.radius)) {
    return std::nullopt;
  }
  const double ux = two.cx - one.cx;
  const double uy = two.cy - one.cy;
  const double d = std::hypot(ux, uy);
  if (!(d > 0.)) {
    return std::nullopt;
  }
  const double ex = ux / d;
  const double ey = uy / d;
  std::vector<std::array<double, 2>> points;
  if (d > one.radius + two.radius) {
    const double along = one.radius + (0.5 * (d - one.radius - two.radius));
    points.push_back({one.cx + (along * ex), one.cy + (along * ey)});
  } else if (d < std::abs(one.radius - two.radius)) {
    const double sign = one.radius > two.radius ? 1. : -1.;
    const double along = 0.5 * ((sign * one.radius) + d + (sign * two.radius));
    points.push_back({one.cx + (along * ex), one.cy + (along * ey)});
  } else {
    const double along = ((one.radius * one.radius) - (two.radius * two.radius) + (d * d)) / (2. * d);
    const double across = std::sqrt(std::max(0., (one.radius * one.radius) - (along * along)));
    points.push_back({one.cx + (along * ex) - (across * ey), one.cy + (along * ey) + (across * ex)});
    points.push_back({one.cx + (along * ex) + (across * ey), one.cy + (along * ey) - (across * ex)});
  }
  // z on a circle's helix at a transverse point: the nearer way round from
  // the daughter's innermost position.
  const auto zAt = [](const Circle& circle, const std::array<double, 2>& point) {
    const double from = std::atan2(circle.py - circle.cy, circle.px - circle.cx);
    const double to = std::atan2(point[1] - circle.cy, point[0] - circle.cx);
    const double turn = std::remainder(to - from, 2. * std::numbers::pi);
    return circle.pz + ((turn / circle.rho) * circle.cotTheta);
  };
  std::optional<std::array<double, 3>> best;
  double bestGap = std::numeric_limits<double>::infinity();
  for (const auto& point : points) {
    const double z1 = zAt(one, point);
    const double z2 = zAt(two, point);
    if (std::abs(z1 - z2) < bestGap) {
      bestGap = std::abs(z1 - z2);
      best = std::array<double, 3>{point[0], point[1], 0.5 * (z1 + z2)};
    }
  }
  return best;
}

} // namespace decay_fit_detail

/// The vertex of \a daughters (at least two) under \a variant, with
/// \a material's material between each daughter's innermost measurement and
/// the vertex.
template <typename Material>
[[nodiscard]] DecayVertexFitResult fitDecayVertex(std::span<const DecayDaughter> daughters, const Material& material, DecayFitVariant variant, int maxIterations = 20)
{
  using namespace decay_fit_detail;
  constexpr int kDerivatives = 6; // vx, vy, vz, phi, theta, q/p
  DecayVertexFitResult result;
  const int count = static_cast<int>(daughters.size());
  const int n = 3 + (3 * count);
  if (count < 2) {
    return result;
  }
  const auto start = closestApproach(daughters[0], daughters[1]);
  if (!start) {
    return result;
  }

  struct Measurement {
    PerigeePlane plane;
    std::array<double, 5> value{};
    std::array<double, 25> weight{}; ///< the inverse covariance
  };
  std::vector<Measurement> measurements(static_cast<std::size_t>(count));
  result.carried.resize(static_cast<std::size_t>(count));
  result.carriedCovariance.resize(static_cast<std::size_t>(count));
  result.carryX2X0.resize(static_cast<std::size_t>(count));

  // Each daughter carried to the perigee plane about (x, y).
  const auto carry = [&](double x, double y) {
    for (int i = 0; i < count; ++i) {
      const auto& daughter = daughters[i];
      const auto plane = detail::hypothesis::perigeePlaneOf(daughter.propagator, daughter.innermost, daughter.surface, static_cast<Scalar>(x), static_cast<Scalar>(y));
      BoundTrackState state = daughter.innermost;
      auto factors = factorCovariance(state.cov);
      PropagationStatus status = PropagationStatus::kOk;
      if (variant == DecayFitVariant::kNoMaterial) {
        status = detail::hypothesis::carryInward(daughter.propagator, o2::gtrk::detail::GivenMaterial{}, daughter.surface, plane, state, factors);
        result.carryX2X0[i] = 0.;
      } else {
        result.carryX2X0[i] = detail::hypothesis::segmentRadiation(daughter.propagator, material, daughter.innermost, daughter.surface, plane);
        status = detail::hypothesis::carryInward(daughter.propagator, material, daughter.surface, plane, state, factors);
      }
      if (status != PropagationStatus::kOk) {
        return DecayFitStatus::kCarryFailed;
      }
      covarianceFromFactors(factors, state.cov);
      std::vector<double> covariance(25);
      for (int r = 0; r < 5; ++r) {
        for (int c = 0; c < 5; ++c) {
          covariance[r * 5 + c] = state.covAt(r, c);
        }
      }
      if (!cholesky(covariance, 5)) {
        return DecayFitStatus::kSingular;
      }
      const auto inverse = choleskyInverse(covariance, 5);
      auto& m = measurements[i];
      m.plane = perigeePlaneFrame(plane);
      for (int r = 0; r < 5; ++r) {
        m.value[r] = state.par[r];
      }
      std::ranges::copy(inverse, m.weight.begin());
      result.carried[i] = {state.phi(), state.theta(), state.qOverP()};
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          result.carriedCovariance[i][r * 3 + c] = state.covAt(r + 2, c + 2);
        }
      }
    }
    return DecayFitStatus::kOk;
  };

  // The parameters: the vertex, then each daughter's phi, theta, q/p.
  std::vector<double> parameters(static_cast<std::size_t>(n));
  std::copy(start->begin(), start->end(), parameters.begin());
  if (const auto status = carry((*start)[0], (*start)[1]); status != DecayFitStatus::kOk) {
    result.status = status;
    return result;
  }
  for (int i = 0; i < count; ++i) {
    for (int k = 0; k < 3; ++k) {
      parameters[3 + (3 * i) + k] = result.carried[i][k];
    }
  }

  // The chi2 at \a p and, with \a normal and \a gradient, the normal
  // equations' matrix and right-hand side; false if a helix misses its plane.
  const auto evaluate = [&](const std::vector<double>& p, double& chi2, std::vector<double>* normal, std::vector<double>* gradient) {
    chi2 = 0.;
    if (normal != nullptr) {
      std::ranges::fill(*normal, 0.);
      std::ranges::fill(*gradient, 0.);
    }
    for (int i = 0; i < count; ++i) {
      const auto& daughter = daughters[i];
      const double fieldScale = static_cast<double>(daughter.propagator.bzKiloGauss) * static_cast<double>(constants::math::B2C);
      using D = Dual<kDerivatives>;
      std::array<D, kDerivatives> in{};
      const int q = 3 + (3 * i);
      const std::array<int, kDerivatives> index{0, 1, 2, q, q + 1, q + 2};
      for (int k = 0; k < kDerivatives; ++k) {
        in[k].v = p[index[k]];
        in[k].d[k] = 1.;
      }
      std::array<D, 5> predicted{};
      if (!predictOnPlane(in[0], in[1], in[2], in[3], in[4], in[5], fieldScale, measurements[i].plane, predicted)) {
        return false;
      }
      std::array<double, 5> residual{};
      for (int r = 0; r < 5; ++r) {
        residual[r] = measurements[i].value[r] - predicted[r].v;
      }
      residual[2] = std::remainder(residual[2], 2. * std::numbers::pi);
      const auto& weight = measurements[i].weight;
      std::array<double, 5> weighted{};
      for (int r = 0; r < 5; ++r) {
        for (int c = 0; c < 5; ++c) {
          weighted[r] += weight[r * 5 + c] * residual[c];
        }
        chi2 += residual[r] * weighted[r];
      }
      if (normal == nullptr) {
        continue;
      }
      // H^T W H and H^T W r over this daughter's six parameters.
      for (int a = 0; a < kDerivatives; ++a) {
        double g = 0.;
        for (int r = 0; r < 5; ++r) {
          g += predicted[r].d[a] * weighted[r];
        }
        (*gradient)[index[a]] += g;
        for (int b = 0; b < kDerivatives; ++b) {
          double h = 0.;
          for (int r = 0; r < 5; ++r) {
            for (int c = 0; c < 5; ++c) {
              h += predicted[r].d[a] * weight[r * 5 + c] * predicted[c].d[b];
            }
          }
          (*normal)[index[a] * n + index[b]] += h;
        }
      }
    }
    return true;
  };

  std::vector<double> normal(static_cast<std::size_t>(n * n));
  std::vector<double> gradient(static_cast<std::size_t>(n));
  std::vector<double> factor;
  double chi2 = 0.;
  bool converged = false;
  for (int iteration = 0; iteration < maxIterations; ++iteration) {
    result.iterations = static_cast<std::uint8_t>(iteration + 1);
    if (variant != DecayFitVariant::kCarryOnce && iteration != 0) {
      if (const auto status = carry(parameters[0], parameters[1]); status != DecayFitStatus::kOk) {
        result.status = status;
        return result;
      }
    }
    if (!evaluate(parameters, chi2, &normal, &gradient)) {
      result.status = DecayFitStatus::kNoCrossing;
      return result;
    }
    factor = normal;
    if (!cholesky(factor, n)) {
      result.status = DecayFitStatus::kSingular;
      return result;
    }
    std::vector<double> step = gradient;
    choleskySolve(factor, n, step);
    // Halve a step that raises the chi2 (the measurements fixed).
    std::vector<double> trial(parameters.size());
    double trialChi2 = 0.;
    double scale = 1.;
    bool taken = false;
    for (int halving = 0; halving < 10; ++halving, scale *= 0.5) {
      for (int k = 0; k < n; ++k) {
        trial[k] = parameters[k] + (scale * step[k]);
      }
      if (evaluate(trial, trialChi2, nullptr, nullptr) && trialChi2 <= chi2 * (1. + 1e-12) + 1e-12) {
        taken = true;
        break;
      }
    }
    if (!taken) {
      converged = true; // no step lowers the chi2: the fit sits at its minimum
      break;
    }
    const auto inverse = choleskyInverse(factor, n);
    double largest = 0.;
    for (int k = 0; k < n; ++k) {
      largest = std::max(largest, std::abs(scale * step[k]) / std::sqrt(inverse[k * n + k]));
    }
    parameters = trial;
    chi2 = trialChi2;
    if (largest < 1e-3) {
      converged = true;
      break;
    }
  }
  result.status = converged ? DecayFitStatus::kOk : DecayFitStatus::kNotConverged;
  // The covariance at the solution, with the measurements of the last carry.
  if (!evaluate(parameters, chi2, &normal, &gradient)) {
    result.status = DecayFitStatus::kNoCrossing;
    return result;
  }
  factor = normal;
  if (!cholesky(factor, n)) {
    result.status = DecayFitStatus::kSingular;
    return result;
  }
  result.covariance = choleskyInverse(factor, n);
  std::copy_n(parameters.begin(), 3, result.vertex.begin());
  result.momenta.resize(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    for (int k = 0; k < 3; ++k) {
      result.momenta[i][k] = parameters[3 + (3 * i) + k];
    }
  }
  result.chi2 = chi2;
  result.ndf = (2 * count) - 3;
  return result;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYVERTEXFIT_H_
