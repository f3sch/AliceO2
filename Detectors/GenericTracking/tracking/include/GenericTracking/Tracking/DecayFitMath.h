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
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Base/DenseSolve.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

namespace o2::gtrk::decay_fit_detail
{

using o2::gtrk::dual::Dual;

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
  using Value = typename dual::ScalarOf<T>::type;
  const T p = T{static_cast<Value>(valueOf(qOverP) > 0. ? numerator : -numerator)} / qOverP;
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
/// inputs, in the arithmetic T.
template <int N, typename T = double>
using FreeState = std::array<Dual<N, T>, 6>;

/// \a from carried along its helix (\a fieldScale = Bz B2C, the curvature per
/// q/pT) to \a plane -- along the beam, across it, or tilted --: the free
/// state at the crossing, ahead or behind, and its loc0, loc1 there; false if
/// the helix does not reach the plane (along the beam: within a quarter turn).
template <int N, typename T>
[[nodiscard]] bool crossToPlane(const FreeState<N, T>& from, const TrackingSurface& plane, double fieldScale, FreeState<N, T>& to, std::array<Dual<N, T>, 2>& local)
{
  using D = Dual<N, T>;
  const T nx = plane.normal[0];
  const T ny = plane.normal[1];
  const T nz = plane.normal[2];
  const D& phi = from[3];
  const D& theta = from[4];
  const D sinTheta = sin(theta);
  const D cotTheta = cos(theta) / sinTheta;
  const D rho = static_cast<T>(fieldScale) * (from[5] / sinTheta);
  D s;
  if (nz == 0.) {
    // Oriented so that the helix leaves through the plane's front (the
    // normal and the direction on the same side, cos(delta) > 0): the
    // crossing then has sin(delta + rho s) = sin(delta) - rho u, u the start's
    // distance before the plane, solved for w = sin(rho s) without the
    // cancellation at small rho (the line's -u / cos(delta) at rho = 0),
    // s = (w / rho) asin(w) / w. asin recovers rho s within a quarter turn only.
    const T phiL = std::atan2(ny, nx);
    const D delta = phi + (-phiL);
    const T side = cos(delta).v >= T(0.) ? T(1.) : T(-1.);
    const D u = side * ((nx * (from[0] + (-static_cast<T>(plane.origin[0])))) + (ny * (from[1] + (-static_cast<T>(plane.origin[1])))));
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
    s = wOverRho * o2::gtrk::detail::asinOverArgument(w);
  } else if (nx == 0. && ny == 0.) {
    if (!(std::abs(cotTheta.v) > 0.)) {
      return false;
    }
    s = (((-1.) * from[2]) + static_cast<T>(plane.origin[2])) / cotTheta;
  } else {
    // A tilted plane: f(s) = n . (x(s) - origin) = 0 by Newton's method from
    // the straight line's crossing, while each step shrinks; s's derivatives
    // by the implicit function, -df / f'(s).
    const auto f = [&](const D& at) {
      const D half = 0.5 * (rho * at);
      const D chord = at * sinc(half);
      return (nx * ((from[0] + (chord * cos(phi + half))) + (-static_cast<T>(plane.origin[0])))) + (ny * ((from[1] + (chord * sin(phi + half))) + (-static_cast<T>(plane.origin[1])))) + (nz * ((from[2] + (at * cotTheta)) + (-static_cast<T>(plane.origin[2]))));
    };
    const auto slope = [&](T at) { return (nx * std::cos(phi.v + (rho.v * at))) + (ny * std::sin(phi.v + (rho.v * at))) + (nz * cotTheta.v); };
    const T ahead = slope(0.);
    if (!(std::abs(ahead) > 0.)) {
      return false;
    }
    T at = -f(D{0.}).v / ahead;
    T last = std::numeric_limits<T>::infinity();
    for (int step = 0; step < kMaxSteps; ++step) {
      const T fs = slope(at);
      if (!(std::abs(fs) > 0.)) {
        return false;
      }
      const T change = f(D{at}).v / fs;
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
      const T half = T(0.5) * rho.v * at;
      const T chord = at * o2::gtrk::sinc(half);
      const std::array<T, 3> x{from[0].v + (chord * std::cos(phi.v + half)), from[1].v + (chord * std::sin(phi.v + half)), from[2].v + (at * cotTheta.v)};
      T size = 0.;
      for (int c = 0; c < 3; ++c) {
        size += std::abs(static_cast<T>(plane.normal[c])) * (std::abs(x[c]) + std::abs(static_cast<T>(plane.origin[c])));
      }
      if (!(std::abs(value.v) <= static_cast<T>(std::numeric_limits<Scalar>::epsilon()) * size)) {
        return false;
      }
    }
    s = D{at};
    const T fs = slope(at);
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
  const std::array<D, 3> d{to[0] + (-static_cast<T>(plane.origin[0])), to[1] + (-static_cast<T>(plane.origin[1])), to[2] + (-static_cast<T>(plane.origin[2]))};
  local[0] = (static_cast<T>(plane.axis0[0]) * d[0]) + (static_cast<T>(plane.axis0[1]) * d[1]) + (static_cast<T>(plane.axis0[2]) * d[2]);
  local[1] = (static_cast<T>(plane.axis1[0]) * d[0]) + (static_cast<T>(plane.axis1[1]) * d[1]) + (static_cast<T>(plane.axis1[2]) * d[2]);
  return true;
}

/// At most \a Capacity elements in place, a runtime count of them used: the
/// fit's storage of what its capacities bound, without allocation. Elements
/// a resize adds are value-initialised, as std::vector's.
template <typename T, std::size_t Capacity>
struct FixedVector {
  std::array<T, Capacity> items{};
  std::size_t count = 0;
  [[nodiscard]] static constexpr std::size_t capacity() { return Capacity; }
  [[nodiscard]] std::size_t size() const { return count; }
  [[nodiscard]] bool empty() const { return count == 0; }
  [[nodiscard]] T* data() { return items.data(); }
  [[nodiscard]] const T* data() const { return items.data(); }
  [[nodiscard]] T* begin() { return items.data(); }
  [[nodiscard]] T* end() { return items.data() + count; }
  [[nodiscard]] const T* begin() const { return items.data(); }
  [[nodiscard]] const T* end() const { return items.data() + count; }
  [[nodiscard]] T& operator[](std::size_t i) { return items[i]; }
  [[nodiscard]] const T& operator[](std::size_t i) const { return items[i]; }
  [[nodiscard]] T& back() { return items[count - 1]; }
  void clear() { count = 0; }
  /// The caller keeps within the capacity.
  void push_back(const T& value) { items[count++] = value; }
  void resize(std::size_t size)
  {
    for (std::size_t i = count; i < size; ++i) {
      items[i] = T{};
    }
    count = size;
  }
};

/// A row-major matrix of \a columns columns at \a a, seen by the dense
/// solves (DenseSolve.h).
template <typename S>
[[nodiscard]] dense::MatrixView<S> viewOf(S* a, int columns)
{
  return {a, columns};
}

/// A symmetric positive definite matrix's Cholesky factor in the arithmetic
/// S, in storage the caller hands in (`l`, at least n x n).
template <typename S = double>
struct Cholesky {
  int n = 0;
  S* l = nullptr;
  /// Factor \a a (\a size x \a size, row-major, its lower triangle read);
  /// false if it is not positive definite.
  [[nodiscard]] bool factor(const S* a, int size)
  {
    n = size;
    std::copy_n(a, n * n, l);
    return dense::choleskyFactor(viewOf(l, n), n);
  }
  [[nodiscard]] S logDeterminant() const { return dense::choleskyLogDeterminant(viewOf(l, n), n); }
  /// x <- a^-1 x.
  void solve(S* x) const { dense::choleskySolve(viewOf(l, n), n, x); }
  /// a^-1 into \a inverse, n x n, row-major.
  void inverse(S* inverse) const { dense::choleskyInverse(viewOf(l, n), n, viewOf(inverse, n)); }
};

/// A track state's 5 x 5 covariance in S, row-major.
template <typename S = double, typename State>
[[nodiscard]] std::array<S, 25> covarianceOf(const State& state)
{
  std::array<S, 25> covariance{};
  for (int r = 0; r < 5; ++r) {
    for (int c = 0; c < 5; ++c) {
      covariance[(r * 5) + c] = state.covAt(r, c);
    }
  }
  return covariance;
}

/// The weight (inverse covariance, 5 x 5, row-major) of the covariance
/// factored in \a cholesky.
template <typename S>
void weightOf(const Cholesky<S>& cholesky, std::array<S, 25>& weight)
{
  dense::choleskyInverse(viewOf(cholesky.l, 5), 5, dense::MatrixView<S>{weight.data(), 5});
}

/// A positive semi-definite matrix's factor L, L L^T = \a a (n x n,
/// row-major), into \a factor: the pivoted L D L^T's P^T L D^1/2
/// (`dense::semidefiniteFactor`; \a a overwritten, \a order n ints of
/// scratch). L is not triangular.
template <typename S>
void semidefiniteFactor(S* a, int n, int* order, S* factor)
{
  dense::semidefiniteFactor(viewOf(a, n), n, order, viewOf(factor, n));
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
/// storage the caller hands in (n and n x n): the mean in double, as the
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
  std::span<S> scatter; ///< n x n, row-major

  MixtureCollapse(std::span<double> meanStorage, std::span<S> scatterStorage) : n(static_cast<int>(meanStorage.size())), mean(meanStorage), scatter(scatterStorage)
  {
    std::ranges::fill(mean, 0.);
    std::ranges::fill(scatter, S(0.));
  }

  /// Takes the component of log-weight \a logWeight, mean \a p, covariance
  /// \a sigma (n x n) and chi2 \a componentChi2.
  void add(double logWeight, std::span<const double> p, std::span<const S> sigma, double componentChi2)
  {
    combine(logWeight, 1., p, sigma, componentChi2, 1);
  }
  void merge(const MixtureCollapse& other)
  {
    if (other.count > 0) {
      combine(other.logScale, other.weight, other.mean, other.scatter, other.chi2, other.count);
    }
  }
  /// The collapse's covariance, S / W, into \a covariance (n x n).
  void covarianceInto(S* covariance) const
  {
    for (std::size_t i = 0; i < scatter.size(); ++i) {
      covariance[i] = static_cast<S>(scatter[i] / weight);
    }
  }
  [[nodiscard]] double meanChi2() const { return chi2 / weight; }

 private:
  void combine(double logScaleB, double weightB, std::span<const double> meanB, std::span<const S> scatterB, double chi2B, int countB)
  {
    if (count == 0) {
      logScale = logScaleB;
      weight = weightB;
      std::ranges::copy(meanB, mean.begin());
      std::ranges::copy(scatterB, scatter.begin());
      chi2 = chi2B;
      count = countB;
      return;
    }
    const double h = std::max(logScale, logScaleB);
    const double a = std::exp(logScale - h);
    const double b = std::exp(logScaleB - h);
    const double u = a * weight;
    const double v = b * weightB;
    const double w = u + v;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        const double di = meanB[i] - mean[i];
        const double dj = meanB[j] - mean[j];
        scatter[(i * n) + j] = static_cast<S>((a * scatter[(i * n) + j]) + (b * scatterB[(i * n) + j]) + ((u * v / w) * di * dj));
      }
    }
    for (int i = 0; i < n; ++i) {
      mean[i] += (v / w) * (meanB[i] - mean[i]);
    }
    chi2 = (a * chi2) + (b * chi2B);
    weight = w;
    logScale = h;
    count += countB;
  }
};

} // namespace o2::gtrk::decay_fit_detail

#endif // GTRK_TRACKING_DECAYFITMATH_H_
