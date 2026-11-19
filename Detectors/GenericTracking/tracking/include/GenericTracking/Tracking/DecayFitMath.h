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
///        to a plane in a uniform field, the dense solves on caller
///        storage (DenseSolve.h), fixed-capacity storage.

#ifndef GTRK_TRACKING_DECAYFITMATH_H_
#define GTRK_TRACKING_DECAYFITMATH_H_

#include "GenericTracking/Base/Dual.h"
#include "GenericTracking/Base/FixedVector.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Base/DenseSolve.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include <array>
#include <cstddef>
#include <limits>
#include <span>
#include <type_traits>

namespace o2::gtrk::decay_fit_detail
{

using o2::gtrk::dual::Dual;

/// Small arithmetic the fit's code shares, on the host and the device.
template <typename From, typename To>
GPUhdi() void copyN(const From* from, std::size_t n, To* to)
{
  for (std::size_t i = 0; i < n; ++i) {
    to[i] = from[i];
  }
}
template <typename T, typename V>
GPUhdi() void fillN(T* to, std::size_t n, const V& value)
{
  for (std::size_t i = 0; i < n; ++i) {
    to[i] = value;
  }
}

/// q/p's numerator for a particle of \a charge: |q|, 1 for a neutral one
/// (its q/p is 1/p).
[[nodiscard]] GPUhdi() double numeratorOf(double charge) { return charge != 0. ? o2::gpu::CAMath::Abs(charge) : 1.; }

/// A momentum: its vector and its magnitude [GeV].
template <typename T>
struct Momentum {
  std::array<T, 3> vector{};
  T magnitude{};
};

/// The momentum of direction \a phi, \a theta and \a qOverP with q/p's
/// \a numerator (`numeratorOf`), in doubles or dual numbers.
template <typename T>
[[nodiscard]] GPUhdi() Momentum<T> momentumOf(const T& phi, const T& theta, const T& qOverP, double numerator)
{
  using dual::valueOf;
  using Value = typename dual::ScalarOf<T>::type;
  const T p = T{static_cast<Value>(valueOf(qOverP) > 0. ? numerator : -numerator)} / qOverP;
  const T pt = p * sinOf(theta);
  return {{pt * cosOf(phi), pt * sinOf(phi), p * cosOf(theta)}, p};
}

/// The momentum of a daughter of \a charge whose phi, theta and q/p are the
/// inputs \a first, \a first + 1 and \a first + 2 of N, at \a phi, \a theta
/// and \a qOverP: its derivatives by them.
template <int N, typename M = double>
[[nodiscard]] GPUhdi() Momentum<Dual<N, M>> momentumOfInputs(double phi, double theta, double qOverP, int first, double charge)
{
  return momentumOf(dual::variable<N, M>(phi, first), dual::variable<N, M>(theta, first + 1), dual::variable<N, M>(qOverP, first + 2), numeratorOf(charge));
}

/// A particle's energy, sqrt(|p|^2 + m^2) with its \a mass [GeV].
template <typename T>
[[nodiscard]] GPUhdi() T energyOf(const Momentum<T>& momentum, double mass)
{
  return sqrtOf((momentum.magnitude * momentum.magnitude) + (mass * mass));
}

/// Momenta summed (a composite's, a resonance's), with their energies where
/// a mass is asked of the sum.
template <typename T>
struct MomentumSum {
  std::array<T, 3> vector{};
  T energy{};
  GPUhdi() void add(const std::array<T, 3>& momentum)
  {
    for (int c = 0; c < 3; ++c) {
      vector[c] = vector[c] + momentum[c];
    }
  }
  GPUhdi() void add(const std::array<T, 3>& momentum, const T& itsEnergy)
  {
    add(momentum);
    energy = energy + itsEnergy;
  }
  /// |P|^2.
  [[nodiscard]] GPUhdi() T squared() const { return (vector[0] * vector[0]) + (vector[1] * vector[1]) + (vector[2] * vector[2]); }
  /// The invariant mass squared, E^2 - |P|^2.
  [[nodiscard]] GPUhdi() T massSquared() const { return (energy * energy) - squared(); }
};

/// A momentum vector's direction and its |P|.
template <typename T>
struct DirectionOfMomentum {
  T phi{};
  T theta{};
  T p{};
};

template <typename T>
[[nodiscard]] GPUhdi() DirectionOfMomentum<T> directionOfMomentum(const std::array<T, 3>& momentum)
{
  const T pt = sqrtOf((momentum[0] * momentum[0]) + (momentum[1] * momentum[1]));
  return {.phi = atan2Of(momentum[1], momentum[0]), .theta = atan2Of(pt, momentum[2]), .p = sqrtOf((pt * pt) + (momentum[2] * momentum[2]))};
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
/// inputs, in the arithmetic T.
template <int N, typename T = double>
using FreeState = std::array<Dual<N, T>, 6>;

/// \a from carried along its helix (\a fieldScale = Bz B2C, the curvature per
/// q/pT) to \a plane -- along the beam, across it, or tilted --: the free
/// state at the crossing, ahead or behind, and its loc0, loc1 there; false if
/// the helix does not reach the plane (along the beam: within a quarter turn).
/// The point of the helix from \a start along \a phi, of curvature \a rho
/// and cot(theta) \a cotTheta, after the transverse path \a s.
template <typename V>
[[nodiscard]] GPUhdi() std::array<V, 3> helixPointOf(const std::array<V, 3>& start, const V& phi, const V& rho, const V& cotTheta, const V& s)
{
  using Value = typename dual::ScalarOf<V>::type;
  const V half = Value(0.5) * (rho * s);
  const V chord = s * sinc(half);
  return {start[0] + (chord * cosOf(phi + half)), start[1] + (chord * sinOf(phi + half)), start[2] + (s * cotTheta)};
}

template <int N, typename T>
[[nodiscard]] GPUhdi() bool crossToPlane(const FreeState<N, T>& from, const TrackingSurface& plane, double fieldScale, FreeState<N, T>& to, std::array<Dual<N, T>, 2>& local)
{
  using D = Dual<N, T>;
  const T nx = plane.normal[0];
  const T ny = plane.normal[1];
  const T nz = plane.normal[2];
  const D& phi = from[3];
  const D& theta = from[4];
  const D sinTheta = sinOf(theta);
  const D cotTheta = cosOf(theta) / sinTheta;
  const D rho = static_cast<T>(fieldScale) * (from[5] / sinTheta);
  D s;
  if (nz == 0.) {
    // Oriented so that the helix leaves through the plane's front (the
    // normal and the direction on the same side, cos(delta) > 0): the
    // crossing then has sin(delta + rho s) = sin(delta) - rho u, u the start's
    // distance before the plane, solved for w = sin(rho s) without the
    // cancellation at small rho (the line's -u / cos(delta) at rho = 0),
    // s = (w / rho) asin(w) / w. asin recovers rho s within a quarter turn only.
    const T phiL = atan2Of(ny, nx);
    const D delta = phi + (-phiL);
    const T side = cosOf(delta).v >= T(0.) ? T(1.) : T(-1.);
    const D u = side * ((nx * (from[0] + (-static_cast<T>(plane.origin[0])))) + (ny * (from[1] + (-static_cast<T>(plane.origin[1])))));
    const D sinDelta = side * sinOf(delta);
    const D cosDelta = side * cosOf(delta);
    const D sinEnd = sinDelta - (rho * u);
    if (!(o2::gpu::CAMath::Abs(sinEnd.v) < 1.) || !(cosDelta.v > 0.)) {
      return false; // missed, grazing, or along the plane: no crossing to bind to
    }
    const D cosEnd = sqrtOf((-(sinEnd * sinEnd)) + 1.);
    if (!((cosEnd.v * cosDelta.v) + (sinEnd.v * sinDelta.v) > 0.)) {
      return false; // beyond a quarter turn
    }
    const D wOverRho = (-(u * cosDelta)) - ((sinDelta * u) * (sinEnd + sinDelta) / (cosDelta + cosEnd));
    const D w = rho * wOverRho;
    if (!(o2::gpu::CAMath::Abs(w.v) < 1.)) {
      return false;
    }
    s = wOverRho * o2::gtrk::detail::asinOverArgument(w);
  } else if (nx == 0. && ny == 0.) {
    if (!(o2::gpu::CAMath::Abs(cotTheta.v) > 0.)) {
      return false;
    }
    s = (((-1.) * from[2]) + static_cast<T>(plane.origin[2])) / cotTheta;
  } else {
    // A tilted plane: f(s) = n . (x(s) - origin) = 0 by Newton's method from
    // the straight line's crossing, while each step shrinks; s's derivatives
    // by the implicit function, -df / f'(s).
    const std::array<D, 3> start{from[0], from[1], from[2]};
    const auto f = [&](const D& at) {
      const auto x = helixPointOf(start, phi, rho, cotTheta, at);
      return (nx * (x[0] + (-static_cast<T>(plane.origin[0])))) + (ny * (x[1] + (-static_cast<T>(plane.origin[1])))) + (nz * (x[2] + (-static_cast<T>(plane.origin[2]))));
    };
    const auto slope = [&](T at) { return (nx * cosOf(phi.v + (rho.v * at))) + (ny * sinOf(phi.v + (rho.v * at))) + (nz * cotTheta.v); };
    const T ahead = slope(0.);
    if (!(o2::gpu::CAMath::Abs(ahead) > 0.)) {
      return false;
    }
    T at = -f(D{0.}).v / ahead;
    T last = std::numeric_limits<T>::infinity();
    for (int step = 0; step < kMaxSteps; ++step) {
      const T fs = slope(at);
      if (!(o2::gpu::CAMath::Abs(fs) > 0.)) {
        return false;
      }
      const T change = f(D{at}).v / fs;
      if (!(o2::gpu::CAMath::Abs(change) < last)) {
        break; // no longer shrinking: at the root to rounding
      }
      last = o2::gpu::CAMath::Abs(change);
      at -= change;
    }
    const D value = f(D{at});
    // Closed to the plane's own precision (its geometry is floats), or no
    // crossing: f's terms' size times float epsilon.
    {
      const auto x = helixPointOf<T>({from[0].v, from[1].v, from[2].v}, phi.v, rho.v, cotTheta.v, at);
      T size = 0.;
      for (int c = 0; c < 3; ++c) {
        size += o2::gpu::CAMath::Abs(static_cast<T>(plane.normal[c])) * (o2::gpu::CAMath::Abs(x[c]) + o2::gpu::CAMath::Abs(static_cast<T>(plane.origin[c])));
      }
      if (!(o2::gpu::CAMath::Abs(value.v) <= static_cast<T>(std::numeric_limits<Scalar>::epsilon()) * size)) {
        return false;
      }
    }
    s = D{at};
    const T fs = slope(at);
    for (int k = 0; k < N; ++k) {
      s.d[k] = -value.d[k] / fs;
    }
  }
  const auto x = helixPointOf<D>({from[0], from[1], from[2]}, phi, rho, cotTheta, s);
  to[0] = x[0];
  to[1] = x[1];
  to[2] = x[2];
  to[3] = phi + (rho * s);
  to[4] = theta;
  to[5] = from[5];
  const std::array<D, 3> d{to[0] + (-static_cast<T>(plane.origin[0])), to[1] + (-static_cast<T>(plane.origin[1])), to[2] + (-static_cast<T>(plane.origin[2]))};
  local[0] = (static_cast<T>(plane.axis0[0]) * d[0]) + (static_cast<T>(plane.axis0[1]) * d[1]) + (static_cast<T>(plane.axis0[2]) * d[2]);
  local[1] = (static_cast<T>(plane.axis1[0]) * d[0]) + (static_cast<T>(plane.axis1[1]) * d[1]) + (static_cast<T>(plane.axis1[2]) * d[2]);
  return true;
}

/// A symmetric positive definite matrix's Cholesky factor in the arithmetic
/// S, packed, in storage the caller hands in (`l`, at least
/// `covSize(n)`).
template <typename S = double>
struct Cholesky {
  int n = 0;
  S* l = nullptr;
  /// Factor \a a (\a size x \a size, packed); false if it is not positive
  /// definite.
  [[nodiscard]] GPUhdi() bool factor(const S* a, int size)
  {
    n = size;
    copyN(a, covSize(n), l);
    return dense::choleskyFactor(dense::PackedView<S>{l}, n);
  }
  [[nodiscard]] GPUhdi() S logDeterminant() const { return dense::choleskyLogDeterminant(dense::PackedView<S>{l}, n); }
  /// x <- a^-1 x.
  GPUhdi() void solve(S* x) const { dense::choleskySolve(dense::PackedView<S>{l}, n, x); }
  /// a^-1 into \a inverse, packed.
  GPUhdi() void inverse(S* inverse) const { dense::choleskyInverse(dense::PackedView<S>{l}, n, dense::PackedView<S>{inverse}); }
};

/// A track state's covariance in S, packed (`covIndex`, the state's own).
template <typename S = double, typename State>
[[nodiscard]] GPUhdi() std::array<S, kBoundCovSize> covarianceOf(const State& state)
{
  std::array<S, kBoundCovSize> covariance{};
  copyN(state.cov, kBoundCovSize, covariance.data());
  return covariance;
}

/// A Gaussian mixture's collapse to one Gaussian (its mean and covariance)
/// and its weighted chi2, accumulated a component at a time or merged from
/// two partial collapses, without keeping the components: the largest
/// log-weight so far h, the weight sum W relative to exp(h), the weighted
/// mean mu, the scatter S = sum w (Sigma + (p - mu)(p - mu)^T) unnormalised,
/// the weighted chi2 sum. Merging A and B rescales both to the larger h and
/// combines means and scatters with their means' difference; adding one
/// component is a merge with a collapse of one. No raw moments
/// (sum w p p^T - mu mu^T), which cancel. Its mean and scatter live in
/// storage the caller hands in (n, and the scatter packed): the mean in double, as the
/// parameters are, the scatter in S, as covariances are; the weights in
/// double.
template <typename S = double>
struct MixtureCollapse {
  int n = 0;
  int count = 0; ///< the components taken
  double logScale = -std::numeric_limits<double>::infinity();
  double weight = 0.;
  double chi2 = 0.;
  std::span<double> mean;
  std::span<S> scatter; ///< packed (`covIndex`)

  GPUhdi() MixtureCollapse(std::span<double> meanStorage, std::span<S> scatterStorage) : n(static_cast<int>(meanStorage.size())), mean(meanStorage), scatter(scatterStorage)
  {
    fillN(mean.data(), mean.size(), 0.);
    fillN(scatter.data(), scatter.size(), S(0.));
  }

  /// Takes the component of log-weight \a logWeight, mean \a p, covariance
  /// \a sigma (packed) and chi2 \a componentChi2.
  GPUhdi() void add(double logWeight, std::span<const double> p, std::span<const S> sigma, double componentChi2)
  {
    combine(logWeight, 1., p, sigma, componentChi2, 1);
  }
  GPUhdi() void merge(const MixtureCollapse& other)
  {
    if (other.count > 0) {
      combine(other.logScale, other.weight, other.mean, other.scatter, other.chi2, other.count);
    }
  }
  /// The collapse's covariance, S / W, into \a covariance (packed).
  GPUhdi() void covarianceInto(S* covariance) const
  {
    for (std::size_t i = 0; i < scatter.size(); ++i) {
      covariance[i] = static_cast<S>(scatter[i] / weight);
    }
  }
  [[nodiscard]] GPUhdi() double meanChi2() const { return chi2 / weight; }

 private:
  GPUhdi() void combine(double logScaleB, double weightB, std::span<const double> meanB, std::span<const S> scatterB, double chi2B, int countB)
  {
    if (count == 0) {
      logScale = logScaleB;
      weight = weightB;
      copyN(meanB.data(), meanB.size(), mean.data());
      copyN(scatterB.data(), scatterB.size(), scatter.data());
      chi2 = chi2B;
      count = countB;
      return;
    }
    const double h = o2::gpu::CAMath::Max(logScale, logScaleB);
    const double a = o2::gpu::CAMath::Expd(logScale - h);
    const double b = o2::gpu::CAMath::Expd(logScaleB - h);
    const double u = a * weight;
    const double v = b * weightB;
    const double w = u + v;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j <= i; ++j) {
        const double di = meanB[i] - mean[i];
        const double dj = meanB[j] - mean[j];
        const int ij = covIndex(i, j);
        scatter[ij] = static_cast<S>((a * scatter[ij]) + (b * scatterB[ij]) + ((u * v / w) * di * dj));
      }
    }
    // the mean anchored on the larger contribution, the smaller added as a
    // correction: a contribution whose weight underflows leaves no trace
    for (int i = 0; i < n; ++i) {
      if (v > u) {
        mean[i] = meanB[i] + ((u / w) * (mean[i] - meanB[i]));
      } else {
        mean[i] += (v / w) * (meanB[i] - mean[i]);
      }
    }
    chi2 = (a * chi2) + (b * chi2B);
    weight = w;
    logScale = h;
    count += countB;
  }
};

} // namespace o2::gtrk::decay_fit_detail

#endif // GTRK_TRACKING_DECAYFITMATH_H_
