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

/// \file GSFTrackFitter.h
/// \brief A Gaussian-sum filter for electrons: outward and inward passes over a chain whose state is a mixture of Gaussians.
///
/// An electron's bremsstrahlung loss is far from Gaussian (BetheHeitler.h),
/// so a single Kalman state either misses the large losses or, given their
/// variance, the small ones. The Gaussian-sum filter carries a mixture: K
/// components, each a weight, a state and its covariance. At each plane
/// every component branches into the loss mixture's components, the mixture
/// is reduced back to K, and each component is predicted to the next plane on
/// its own trajectory and updated there as the Kalman filter does
/// (KFTrackFitter.h); its weight takes the measurement's likelihood under it.
///
/// Where the radiation is taken: the material between two measurements is cut
/// at the midpoint of their straight segment, and a plane's radiation is that
/// of the half-gaps on both sides of it, taken right after its update. The
/// half the pass is about to cross acts before the propagation it bends; the
/// half it just crossed is deferred past the update, which its loss does not
/// move (a loss at the end of a step moves no point on that step's plane).
/// Consecutive losses multiply into the loss of the summed thickness, so the
/// two halves are one branching. The ionisation stays the Kalman step's.

#ifndef GTRK_TRACKING_GSFTRACKFITTER_H_
#define GTRK_TRACKING_GSFTRACKFITTER_H_

#include "GenericTracking/Tracking/BetheHeitler.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/KFTrackFitter.h"
#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Tracking/MeasurementUpdate.h"
#include "GenericTracking/Tracking/Seed.h"

#include "GPUCommonMath.h"
#include "MathUtils/detail/trigonometric.h"

#include <cfloat>
#include <cstdint>
#include <numbers>

namespace o2::gtrk
{

/// The radiation model of the filter: Bethe and Heitler's loss as a mixture
/// (`betheHeitlerMixture`). A model gives at most `kComponents` components
/// for a thickness up to `kMaxX2X0`; the filter takes a thicker one in equal
/// parts.
struct BetheHeitlerLoss {
  static constexpr int kComponents = kBetheHeitlerComponents;
  static constexpr Scalar kMaxX2X0 = kBetheHeitlerMaxX2X0;
  GPUdi() static int mixture(Scalar x2X0, BetheHeitlerComponent (&components)[kComponents]) { return betheHeitlerMixture(x2X0, components); }
};

/// A mixture of at most \a Capacity Gaussian states, all on one plane: a
/// component's weight (its logarithm), its state and its covariance as
/// factors; `state.cov` is written only at a pass's end.
template <int Capacity>
struct GSFMixture {
  static constexpr int kCapacity = Capacity;
  int size = 0;
  Scalar logWeight[Capacity]{};
  BoundTrackState state[Capacity]{};
  CovarianceFactors factors[Capacity]{};
  Scalar chi2[Capacity]{}; ///< each component's summed chi2 increments
  /// Each component's reference q/p, the momentum its steps take the material
  /// at (`predictStep`): the pass's seed's, scaled with each loss the
  /// component took, carried by each accepted step.
  Scalar reference[Capacity]{};
};

/// One pass of the filter: its end mixture on the pass's last measurement,
/// with that plane's radiation taken (`state.cov` filled), and its summaries.
template <int K>
struct GSFPassResult {
  GSFMixture<K> mixture;
  /// The mixture as one Gaussian: its mean and covariance.
  BoundTrackState collapsed{};
  /// The modes of the q/p marginal and of the momentum's density [1/GeV, GeV]:
  /// the latter is not the inverse of the former, the density of p being the
  /// q/p one's times (q/p)^2.
  Scalar qOverPMode = 0.F;
  Scalar momentumMode = 0.F;
  /// -2 ln of the measurements' likelihood under the pass's predictions,
  /// the start's rows by the diffuse likelihood (`startMixture`), without the
  /// constant 2 ln 2 pi per coordinate pair, given that the components it
  /// kept propagated (a failed one's weight is dropped).
  Scalar minusTwoLogLikelihood = 0.F;
  PropagationStatus status = PropagationStatus::kOk;
  /// The measurement the failed step was propagating to.
  std::int16_t failedMeasurement = -1;
};

/// The outward and the inward pass, each with nothing known before its own
/// first measurements (as `fitMeasurements`): the inward pass is the one at
/// the vertex end. If the outward pass fails, the inward one does not run
/// and carries its status.
template <int K>
struct GSFFitResult {
  GSFPassResult<K> outward;
  GSFPassResult<K> inward;
};

namespace detail::gsf
{

/// A step's loss as one Gaussian: the mixture's mean and variance of z.
template <int M>
GPUhdi() BetheHeitlerComponent momentsOf(const BetheHeitlerComponent (&components)[M], int size)
{
  Scalar mean = 0.F;
  Scalar second = 0.F;
  for (int k = 0; k < size; ++k) {
    mean += components[k].weight * components[k].mean;
    second += components[k].weight * (components[k].variance + (components[k].mean * components[k].mean));
  }
  const Scalar variance = second - (mean * mean);
  return {.weight = 1.F, .mean = mean, .variance = variance > 0.F ? variance : 0.F};
}

/// The factor a loss component scales q/p by: 1 / z along the momentum, z
/// against it (`applyBetheHeitler`); a component's reference q/p takes it too.
GPUhdi() Scalar lossScale(const BetheHeitlerComponent& component, PropagationDirection direction)
{
  return direction == PropagationDirection::kForward ? 1.F / component.mean : component.mean;
}

/// A loss component on a diffuse covariance (`DiffuseFactors`): the
/// q/p scaling on both parts -- P_inf's root and P_star -- and the variance
/// on P_star, as `applyBetheHeitler` takes them on a full one.
GPUhdi() bool applyLossDiffuse(const BetheHeitlerComponent& component, PropagationDirection direction, BoundTrackState& state, DiffuseFactors& c)
{
  const bool along = direction == PropagationDirection::kForward;
  const Scalar scale = along ? 1.F / component.mean : component.mean;
  for (int k = 0; k < c.rank; ++k) {
    c.root[kQOverP][k] *= scale;
  }
  const Scalar qOverP = state.par[kQOverP] * scale;
  state.par[kQOverP] = qOverP;
  for (int row = 0; row < kQOverP; ++row) {
    c.star.upper[strictUpperIndex(row, kQOverP)] /= scale;
  }
  c.star.diagonal[kQOverP] *= scale * scale;
  const Scalar slope = along ? qOverP * scale : qOverP / scale;
  Scalar unit[kNBoundParams]{};
  unit[kQOverP] = 1.F;
  rankOneUpdateUD(c.star.upper, c.star.diagonal, slope * slope * component.variance, unit);
  return isFiniteScalar(qOverP);
}

/// The straight segment's midpoint between measurements \a a and \a b.
template <DetectorLayout Det, typename Chain>
GPUdi() void segmentMidpoint(const TrackingGeometry<Det>& geometry, const Chain& measurements, int a, int b, Scalar (&mid)[3])
{
  Scalar pa[3];
  Scalar pb[3];
  surfaceToGlobal(geometry.planeOf(measurements[a]), measurements[a].loc0, measurements[a].loc1, pa);
  surfaceToGlobal(geometry.planeOf(measurements[b]), measurements[b].loc0, measurements[b].loc1, pb);
  for (int i = 0; i < 3; ++i) {
    mid[i] = 0.5F * (pa[i] + pb[i]);
  }
}

/// Plane \a i's radiation [x/X0]: \a lut's material from the midpoint of
/// the segment to the measurement before to the midpoint of the segment to
/// the one after, along the two straight half-segments; nothing beyond the
/// chain's ends.
template <DetectorLayout Det, typename Material, typename Chain>
GPUdi() Scalar planeRadiation(const TrackingGeometry<Det>& geometry, const Material& lut, const Chain& measurements, int n, int i)
{
  Scalar point[3];
  surfaceToGlobal(geometry.planeOf(measurements[i]), measurements[i].loc0, measurements[i].loc1, point);
  Scalar x2X0 = 0.F;
  if (i > 0) {
    Scalar mid[3];
    segmentMidpoint(geometry, measurements, i - 1, i, mid);
    x2X0 += lut.getMatBudget(mid[0], mid[1], mid[2], point[0], point[1], point[2]).meanX2X0;
  }
  if (i + 1 < n) {
    Scalar mid[3];
    segmentMidpoint(geometry, measurements, i, i + 1, mid);
    x2X0 += lut.getMatBudget(point[0], point[1], point[2], mid[0], mid[1], mid[2]).meanX2X0;
  }
  return x2X0;
}

/// How many equal parts a plane's radiation \a x2X0 is taken in: the fewest
/// no thicker than the model's range.
template <typename Loss>
GPUhdi() int lossParts(Scalar x2X0)
{
  int parts = 1;
  if (!isFiniteScalar(x2X0)) {
    return parts; // the model refuses it
  }
  while (x2X0 > static_cast<Scalar>(parts) * Loss::kMaxX2X0) {
    ++parts;
  }
  return parts;
}

/// a's components, b's and the spread of their means as one Gaussian of
/// weight \a wa + \a wb (relative weights): the weighted mean, phi taken
/// relative to a's, and the covariance (wa Pa + wb Pb) / W + wa wb / W^2 d
/// d^T, d the means' difference, its factors by weighted Gram-Schmidt over
/// the factor columns of both and d (`factorColumns`): no covariance is
/// subtracted.
GPUdi() void mergeInto(Scalar wa, Scalar wb, BoundTrackState& a, CovarianceFactors& fa, const BoundTrackState& b, const CovarianceFactors& fb)
{
  const Scalar total = wa + wb;
  if (!(total > 0.F)) {
    return;
  }
  const Scalar sa = wa / total;
  const Scalar sb = wb / total;
  Scalar difference[kNBoundParams];
  for (int i = 0; i < kNBoundParams; ++i) {
    difference[i] = a.par[i] - b.par[i];
  }
  o2::math_utils::detail::bringToPMPi<Scalar>(difference[kPhi]);
  constexpr int kColumns = (2 * kNBoundParams) + 1;
  Scalar y[kNBoundParams][kColumns]{};
  Scalar w[kColumns]{};
  for (int k = 0; k < kNBoundParams; ++k) {
    for (int i = 0; i < k; ++i) {
      y[i][k] = fa.upper[strictUpperIndex(i, k)];
      y[i][kNBoundParams + k] = fb.upper[strictUpperIndex(i, k)];
    }
    y[k][k] = 1.F;
    y[k][kNBoundParams + k] = 1.F;
    w[k] = sa * fa.diagonal[k];
    w[kNBoundParams + k] = sb * fb.diagonal[k];
    y[k][2 * kNBoundParams] = difference[k];
  }
  w[2 * kNBoundParams] = sa * sb;
  factorColumns(y, w, fa);
  // a's mean moves towards b's by b's share of the difference.
  for (int i = 0; i < kNBoundParams; ++i) {
    a.par[i] -= sb * difference[i];
  }
  o2::math_utils::detail::bringToPMPi<Scalar>(a.par[kPhi]);
}

/// Runnalls' cost of merging two components, on the q/p marginal: the
/// discrimination the mixture loses, B = [W ln s_ab - wa ln s_a - wb ln s_b]
/// / 2 with s the q/p variances and s_ab the merged one, means' spread
/// included. A light component costs little, so negligible ones merge first.
GPUhdi() Scalar mergeCost(Scalar wa, Scalar ma, Scalar sa, Scalar logSa, Scalar wb, Scalar mb, Scalar sb, Scalar logSb)
{
  const Scalar total = wa + wb;
  if (!(total > 0.F)) {
    return 0.F;
  }
  const Scalar d = ma - mb;
  const Scalar merged = ((wa * sa) + (wb * sb)) / total + (wa * wb / (total * total)) * d * d;
  return 0.5F * ((total * o2::gpu::GPUCommonMath::Log(merged)) - (wa * logSa) - (wb * logSb));
}

/// \a from's components reduced to at most \a K into \a to: until K remain,
/// the pair of least `mergeCost` merges (`mergeInto`). The pair costs sit in
/// a triangular table over \a from's slots; a merge recomputes the merged
/// component's costs and retires the other's.
template <int K, int N>
GPUdi() void reduce(GSFMixture<N>& from, GSFMixture<K>& to)
{
  Scalar maxLog = -FLT_MAX;
  for (int k = 0; k < from.size; ++k) {
    maxLog = o2::gpu::GPUCommonMath::Max(maxLog, from.logWeight[k]);
  }
  Scalar weight[N];
  Scalar logVariance[N];
  bool alive[N];
  bool merged[N];
  for (int k = 0; k < from.size; ++k) {
    weight[k] = o2::gpu::GPUCommonMath::Exp(from.logWeight[k] - maxLog);
    logVariance[k] = o2::gpu::GPUCommonMath::Log(from.factors[k].diagonal[kQOverP]);
    alive[k] = true;
    merged[k] = false;
  }
  int remaining = from.size;
  if (remaining > K) {
    constexpr int kPairs = N * (N - 1) / 2;
    Scalar cost[kPairs > 0 ? kPairs : 1];
    const auto pairCost = [&](int a, int b) {
      return mergeCost(weight[a], from.state[a].qOverP(), from.factors[a].diagonal[kQOverP], logVariance[a],
                       weight[b], from.state[b].qOverP(), from.factors[b].diagonal[kQOverP], logVariance[b]);
    };
    for (int b = 1; b < from.size; ++b) {
      for (int a = 0; a < b; ++a) {
        cost[strictUpperIndex(a, b)] = pairCost(a, b);
      }
    }
    while (remaining > K) {
      // The least cost; a NaN never wins, the first alive pair stands in.
      int bestA = -1;
      int bestB = -1;
      Scalar best = FLT_MAX;
      for (int b = 1; b < from.size; ++b) {
        for (int a = 0; a < b && alive[b]; ++a) {
          if (alive[a] && (bestA < 0 || cost[strictUpperIndex(a, b)] < best)) {
            best = cost[strictUpperIndex(a, b)] < best ? cost[strictUpperIndex(a, b)] : best;
            bestA = a;
            bestB = b;
          }
        }
      }
      mergeInto(weight[bestA], weight[bestB], from.state[bestA], from.factors[bestA], from.state[bestB], from.factors[bestB]);
      from.chi2[bestA] = ((weight[bestA] * from.chi2[bestA]) + (weight[bestB] * from.chi2[bestB])) / (weight[bestA] + weight[bestB]);
      from.reference[bestA] = ((weight[bestA] * from.reference[bestA]) + (weight[bestB] * from.reference[bestB])) / (weight[bestA] + weight[bestB]);
      weight[bestA] += weight[bestB];
      merged[bestA] = true;
      logVariance[bestA] = o2::gpu::GPUCommonMath::Log(from.factors[bestA].diagonal[kQOverP]);
      alive[bestB] = false;
      --remaining;
      for (int k = 0; k < from.size; ++k) {
        if (alive[k] && k != bestA) {
          cost[k < bestA ? strictUpperIndex(k, bestA) : strictUpperIndex(bestA, k)] = k < bestA ? pairCost(k, bestA) : pairCost(bestA, k);
        }
      }
    }
  }
  to.size = 0;
  for (int k = 0; k < from.size; ++k) {
    if (!alive[k]) {
      continue;
    }
    const int slot = to.size++;
    // An unmerged component keeps its weight's logarithm as it was, however small.
    to.logWeight[slot] = merged[k] ? maxLog + o2::gpu::GPUCommonMath::Log(weight[k]) : from.logWeight[k];
    to.state[slot] = from.state[k];
    to.factors[slot] = from.factors[k];
    to.chi2[slot] = from.chi2[k];
    to.reference[slot] = from.reference[k];
  }
}

/// Every component of \a mixture with every loss component of the radiation
/// \a x2X0, reduced back to at most K (`reduce`); the radiation in equal
/// parts if it is thicker than the model's range. False if a component's
/// loss gave no finite state.
template <typename Loss, int K>
GPUdi() bool takeRadiation(Scalar x2X0, PropagationDirection direction, GSFMixture<K>& mixture)
{
  constexpr int M = Loss::kComponents;
  const int parts = lossParts<Loss>(x2X0);
  BetheHeitlerComponent loss[M];
  const int m = Loss::mixture(x2X0 / static_cast<Scalar>(parts), loss);
  if (m <= 0) {
    return false;
  }
  for (int part = 0; part < parts; ++part) {
    GSFMixture<K * M> branched;
    for (int k = 0; k < mixture.size; ++k) {
      for (int j = 0; j < m; ++j) {
        const int slot = branched.size++;
        branched.state[slot] = mixture.state[k];
        branched.factors[slot] = mixture.factors[k];
        branched.chi2[slot] = mixture.chi2[k];
        branched.reference[slot] = mixture.reference[k] * lossScale(loss[j], direction);
        branched.logWeight[slot] = mixture.logWeight[k] + o2::gpu::GPUCommonMath::Log(loss[j].weight);
        if (!applyBetheHeitler(loss[j], direction, branched.state[slot], branched.factors[slot])) {
          return false;
        }
      }
    }
    reduce(branched, mixture);
  }
  return true;
}

/// One component's step: the Kalman prediction on its own trajectory and the
/// update (`filterStep`'s arithmetic), with the update's ln det S for its
/// likelihood: \a logLikelihood = -(chi2 + ln det S) / 2. The component's
/// reference q/p \a referenceQOverP moves with an accepted step only.
template <CrossingRule Rule, DetectorLayout Det, typename Material, typename Prop>
GPUdi() PropagationStatus componentStep(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& state, CovarianceFactors& factors, Scalar& chi2, Scalar& referenceQOverP, Scalar& logLikelihood)
{
  SolenoidJacobian jacobian;
  StepNoise noise;
  Scalar nextReference = referenceQOverP;
  const auto predicted = detail::predictStep<false, Rule>(prop, geometry, material, measurement, direction, state, nextReference, jacobian, noise, 1.F);
  if (predicted != PropagationStatus::kOk) {
    return predicted;
  }
  if (!transportFactors(jacobian, noise, factors)) {
    return PropagationStatus::kNotFinite;
  }
  Scalar determinant = 0.F;
  const auto update = updateWithMeasurement<true>(state, factors, measurement, &determinant);
  if (!updateUsable(update, factors)) {
    return PropagationStatus::kNotFinite;
  }
  referenceQOverP = nextReference;
  chi2 += update.chi2;
  logLikelihood = -0.5F * (update.chi2 + o2::gpu::GPUCommonMath::Log(determinant));
  return PropagationStatus::kOk;
}

/// One branch's step in the diffuse start: the Kalman diffuse step
/// (`diffuseStep`'s arithmetic), with the update's terms of -2 ln of the
/// diffuse likelihood added to \a density (`updateDiffuse`); the branch's
/// reference q/p \a referenceQOverP moves with an accepted step only. Never
/// inlined on the device, as `diffuseStep`.
template <CrossingRule Rule, DetectorLayout Det, typename Material, typename Prop>
GPUdnii() PropagationStatus diffuseBranchStep(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& state, DiffuseFactors& c, Scalar& chi2, Scalar& referenceQOverP, Scalar& density)
{
  SolenoidJacobian jacobian;
  StepNoise noise;
  Scalar nextReference = referenceQOverP;
  const auto predicted = detail::predictStep<false, Rule>(prop, geometry, material, measurement, direction, state, nextReference, jacobian, noise, 1.F);
  if (predicted != PropagationStatus::kOk) {
    return predicted;
  }
  transportDiffuse(jacobian, noise, c);
  if (!updateDiffuse<true>(c, measurement, state, chi2, &density)) {
    return PropagationStatus::kNotFinite;
  }
  if (c.rank > 0) {
    o2::math_utils::detail::bringToPMPi<Scalar>(state.par[kPhi]);
    referenceQOverP = nextReference;
    return PropagationStatus::kOk;
  }
  if (!positiveDiagonal(c.star)) {
    return PropagationStatus::kNotFinite;
  }
  detail::canonicalizeDirection(state, c.star);
  referenceQOverP = nextReference;
  return PropagationStatus::kOk;
}

/// The weights normalised to sum to one; returns the logarithm of what
/// they summed to.
template <typename Mixture>
GPUhdi() Scalar normalise(Mixture& mixture)
{
  Scalar maxLog = -FLT_MAX;
  for (int k = 0; k < mixture.size; ++k) {
    maxLog = o2::gpu::GPUCommonMath::Max(maxLog, mixture.logWeight[k]);
  }
  Scalar sum = 0.F;
  for (int k = 0; k < mixture.size; ++k) {
    sum += o2::gpu::GPUCommonMath::Exp(mixture.logWeight[k] - maxLog);
  }
  const Scalar logSum = maxLog + o2::gpu::GPUCommonMath::Log(sum);
  for (int k = 0; k < mixture.size; ++k) {
    mixture.logWeight[k] -= logSum;
  }
  return logSum;
}

/// A mixture of diffuse states (`DiffuseFactors`), a pass's start.
template <int Capacity>
struct DiffuseMixture {
  int size = 0;
  Scalar logWeight[Capacity]{};
  BoundTrackState state[Capacity]{};
  DiffuseFactors factors[Capacity]{};
  Scalar chi2[Capacity]{};
  Scalar reference[Capacity]{}; ///< as `GSFMixture::reference`
};

GPUhdi() constexpr int power(int base, int exponent) { return exponent == 0 ? 1 : base * power(base, exponent - 1); }

/// The start's branches: every loss path over the planes a start of
/// `kPassStart` rows leaves.
template <typename Loss>
inline constexpr int kStartBranches = power(Loss::kComponents, kPassStart - 1);

/// The pass's start: the Kalman diffuse start (`startDiffuse`, `diffuseStep`)
/// over its first measurements until nothing is undetermined, each plane it
/// leaves branching on its radiation as later planes do. A loss on these
/// planes changes the curvature the start fits, so the branches end with
/// different states. Their weights take the diffuse likelihood (de Jong;
/// Koopman and Durbin): a coordinate that determines a diffuse direction
/// counts ln(h P_inf h), which the branches' different q/p scalings of P_inf
/// make differ, a tested one ln F + v^2 / F (`diffuseBranchStep`). A start longer than `kPassStart`
/// rows (a coordinate that determined nothing) takes the loss of the planes
/// beyond as one Gaussian with the mixture's moments. Ends reduced to at
/// most K components on the last plane it took, before that plane's
/// radiation; \a minusTwoLogLikelihood gets the start's -2 ln likelihood.
template <typename Loss, CrossingRule Rule, typename Prop, DetectorLayout Det, typename Material, typename Lut, typename Chain, int K>
GPUdi() PropagationStatus startMixture(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Lut& lut, const Chain& measurements, int n, PropagationDirection direction, BoundTrackState seed, GSFMixture<K>& mixture, int& taken, std::int16_t& failedMeasurement, Scalar& minusTwoLogLikelihood)
{
  constexpr int M = Loss::kComponents;
  constexpr int S = kStartBranches<Loss>;
  const int first = direction == PropagationDirection::kForward ? 0 : n - 1;
  const int increment = direction == PropagationDirection::kForward ? 1 : -1;
  DiffuseMixture<S> branches;
  branches.size = 1;
  branches.state[0] = seed;
  branches.reference[0] = seed.qOverP(); // the pass's reference, as `startPass`'s
  taken = 0;
  if (n < kPassStart || !startDiffuse(measurements[first], branches.state[0], branches.factors[0])) {
    failedMeasurement = static_cast<std::int16_t>(first);
    return PropagationStatus::kNotFinite;
  }
  const auto undetermined = [&]() {
    bool any = false;
    for (int b = 0; b < branches.size; ++b) {
      any = any || branches.factors[b].rank > 0;
    }
    return any;
  };
  for (taken = 1; taken < n && undetermined(); ++taken) {
    const int previous = first + ((taken - 1) * increment);
    const int index = previous + increment;
    const Scalar x2X0 = planeRadiation(geometry, lut, measurements, n, previous);
    const int parts = lossParts<Loss>(x2X0);
    BetheHeitlerComponent loss[M];
    const int m = Loss::mixture(x2X0 / static_cast<Scalar>(parts), loss);
    if (m <= 0) {
      failedMeasurement = static_cast<std::int16_t>(previous);
      return PropagationStatus::kNotFinite;
    }
    for (int part = 0; part < parts; ++part) {
      if (branches.size * m <= S) {
        DiffuseMixture<S> branched;
        for (int b = 0; b < branches.size; ++b) {
          for (int j = 0; j < m; ++j) {
            const int slot = branched.size++;
            branched.state[slot] = branches.state[b];
            branched.factors[slot] = branches.factors[b];
            branched.chi2[slot] = branches.chi2[b];
            branched.reference[slot] = branches.reference[b] * lossScale(loss[j], direction);
            branched.logWeight[slot] = branches.logWeight[b] + o2::gpu::GPUCommonMath::Log(loss[j].weight);
            if (!applyLossDiffuse(loss[j], direction, branched.state[slot], branched.factors[slot])) {
              failedMeasurement = static_cast<std::int16_t>(previous);
              return PropagationStatus::kNotFinite;
            }
          }
        }
        branches = branched;
      } else {
        const auto moments = momentsOf(loss, m);
        for (int b = 0; b < branches.size; ++b) {
          if (!applyLossDiffuse(moments, direction, branches.state[b], branches.factors[b])) {
            failedMeasurement = static_cast<std::int16_t>(previous);
            return PropagationStatus::kNotFinite;
          }
          branches.reference[b] *= lossScale(moments, direction);
        }
      }
    }
    const auto& stepProp = stepPropagator(prop, geometry, measurements, previous, index);
    int kept = 0;
    PropagationStatus failure = PropagationStatus::kOk;
    for (int b = 0; b < branches.size; ++b) {
      Scalar density = 0.F;
      const auto status = diffuseBranchStep<Rule>(stepProp, geometry, material, measurements[index], direction, branches.state[b], branches.factors[b], branches.chi2[b], branches.reference[b], density);
      if (status != PropagationStatus::kOk) {
        failure = status;
        continue;
      }
      branches.logWeight[kept] = branches.logWeight[b] - (0.5F * density);
      branches.state[kept] = branches.state[b];
      branches.factors[kept] = branches.factors[b];
      branches.chi2[kept] = branches.chi2[b];
      branches.reference[kept] = branches.reference[b];
      ++kept;
    }
    branches.size = kept;
    if (kept == 0) {
      failedMeasurement = static_cast<std::int16_t>(index);
      return failure;
    }
    minusTwoLogLikelihood -= 2.F * normalise(branches);
  }
  if (undetermined()) {
    failedMeasurement = static_cast<std::int16_t>(first + ((n - 1) * increment));
    return PropagationStatus::kNotFinite;
  }
  GSFMixture<S> determined;
  determined.size = branches.size;
  for (int b = 0; b < branches.size; ++b) {
    determined.logWeight[b] = branches.logWeight[b];
    determined.state[b] = branches.state[b];
    determined.factors[b] = branches.factors[b].star;
    determined.chi2[b] = branches.chi2[b];
    determined.reference[b] = branches.reference[b];
  }
  reduce(determined, mixture);
  return PropagationStatus::kOk;
}

/// The mixture as one Gaussian: the weighted mean (phi about the heaviest
/// component's, wrapped) and the covariance with the means' spread.
template <int K>
GPUhdi() BoundTrackState collapse(const GSFMixture<K>& mixture)
{
  int heaviest = 0;
  for (int k = 1; k < mixture.size; ++k) {
    if (mixture.logWeight[k] > mixture.logWeight[heaviest]) {
      heaviest = k;
    }
  }
  BoundTrackState out = mixture.state[heaviest];
  Scalar mean[kNBoundParams]{};
  Scalar weights[K];
  for (int k = 0; k < mixture.size; ++k) {
    weights[k] = o2::gpu::GPUCommonMath::Exp(mixture.logWeight[k]);
    for (int i = 0; i < kNBoundParams; ++i) {
      Scalar value = mixture.state[k].par[i] - out.par[i];
      if (i == kPhi) {
        o2::math_utils::detail::bringToPMPi<Scalar>(value);
      }
      mean[i] += weights[k] * value;
    }
  }
  Scalar cov[kBoundCovSize]{};
  for (int k = 0; k < mixture.size; ++k) {
    Scalar d[kNBoundParams];
    for (int i = 0; i < kNBoundParams; ++i) {
      d[i] = mixture.state[k].par[i] - out.par[i];
      if (i == kPhi) {
        o2::math_utils::detail::bringToPMPi<Scalar>(d[i]);
      }
      d[i] -= mean[i];
    }
    for (int i = 0; i < kNBoundParams; ++i) {
      for (int j = 0; j <= i; ++j) {
        cov[covIndex(i, j)] += weights[k] * (mixture.state[k].cov[covIndex(i, j)] + (d[i] * d[j]));
      }
    }
  }
  for (int i = 0; i < kNBoundParams; ++i) {
    out.par[i] += mean[i];
  }
  o2::math_utils::detail::bringToPMPi<Scalar>(out.par[kPhi]);
  for (int i = 0; i < kBoundCovSize; ++i) {
    out.cov[i] = cov[i];
  }
  return out;
}

/// The mode of the density sum_k w_k N(u; m_k, s_k) |u|^Power over u = q/p
/// (Power 0: the q/p marginal; 2: the momentum's density, in u): from each
/// component's mean, Newton's steps on the log density, a step that does not
/// raise the density halved, until the step no longer moves u; the highest
/// of the ends.
template <int Power, int K>
GPUhdi() Scalar marginalMode(const GSFMixture<K>& mixture)
{
  const auto logDensity = [&](Scalar u, Scalar& gradient, Scalar& curvature) {
    Scalar maxLog = -FLT_MAX;
    Scalar logTerm[K];
    for (int k = 0; k < mixture.size; ++k) {
      const Scalar s = mixture.state[k].cov[covIndex(kQOverP, kQOverP)];
      const Scalar d = u - mixture.state[k].qOverP();
      logTerm[k] = mixture.logWeight[k] - (0.5F * o2::gpu::GPUCommonMath::Log(s)) - (0.5F * d * d / s);
      maxLog = o2::gpu::GPUCommonMath::Max(maxLog, logTerm[k]);
    }
    Scalar f = 0.F;
    Scalar f1 = 0.F;
    Scalar f2 = 0.F;
    for (int k = 0; k < mixture.size; ++k) {
      const Scalar s = mixture.state[k].cov[covIndex(kQOverP, kQOverP)];
      const Scalar d = u - mixture.state[k].qOverP();
      const Scalar term = o2::gpu::GPUCommonMath::Exp(logTerm[k] - maxLog);
      f += term;
      f1 -= term * d / s;
      f2 += term * (((d * d) / (s * s)) - (1.F / s));
    }
    // ln g = ln f + Power ln |u|
    const Scalar g1 = (f1 / f) + (Power != 0 ? static_cast<Scalar>(Power) / u : 0.F);
    const Scalar g2 = (f2 / f) - ((f1 * f1) / (f * f)) - (Power != 0 ? static_cast<Scalar>(Power) / (u * u) : 0.F);
    gradient = g1;
    curvature = g2;
    return maxLog + o2::gpu::GPUCommonMath::Log(f) + (Power != 0 ? static_cast<Scalar>(Power) * o2::gpu::GPUCommonMath::Log(o2::gpu::GPUCommonMath::Abs(u)) : 0.F);
  };
  Scalar best = mixture.state[0].qOverP();
  Scalar bestValue = -FLT_MAX;
  for (int k = 0; k < mixture.size; ++k) {
    Scalar u = mixture.state[k].qOverP();
    Scalar gradient;
    Scalar curvature;
    Scalar value = logDensity(u, gradient, curvature);
    // Each accepted step raises the density; a step is at least halved until it does, and the walk ends when a
    // step leaves u unchanged in float, which a bounded number of halvings reaches.
    for (;;) {
      Scalar step = curvature < 0.F ? -gradient / curvature : gradient * mixture.state[k].cov[covIndex(kQOverP, kQOverP)];
      if (!isFiniteScalar(step)) {
        break;
      }
      Scalar next = u + step;
      Scalar nextGradient;
      Scalar nextCurvature;
      Scalar nextValue = logDensity(next, nextGradient, nextCurvature);
      while (!(nextValue > value) && next != u) {
        step *= 0.5F;
        next = u + step;
        nextValue = logDensity(next, nextGradient, nextCurvature);
      }
      if (next == u) {
        break;
      }
      u = next;
      value = nextValue;
      gradient = nextGradient;
      curvature = nextCurvature;
    }
    if (value > bestValue) {
      bestValue = value;
      best = u;
    }
  }
  return best;
}

/// One pass of the filter over the chain: the start (`startMixture`), then
/// on every plane its radiation (`takeRadiation`) and each component's step to
/// the next (`componentStep`) with its weight's likelihood; at the end the
/// last plane's radiation, the covariances written and the summaries.
template <int K, typename Loss, CrossingRule Rule, typename Prop, DetectorLayout Det, typename Material, typename Lut, typename Chain>
GPUdi() void runPass(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Lut& lut, const Chain& measurements, int n, PropagationDirection direction, const BoundTrackState& seed, GSFPassResult<K>& pass)
{
  const int first = direction == PropagationDirection::kForward ? 0 : n - 1;
  const int increment = direction == PropagationDirection::kForward ? 1 : -1;
  auto& mixture = pass.mixture;
  int taken = 0;
  pass.minusTwoLogLikelihood = 0.F;
  pass.status = startMixture<Loss, Rule>(prop, geometry, material, lut, measurements, n, direction, seed, mixture, taken, pass.failedMeasurement, pass.minusTwoLogLikelihood);
  if (pass.status != PropagationStatus::kOk) {
    return;
  }
  for (int step = taken; step <= n; ++step) {
    const int previous = first + ((step - 1) * increment);
    if (!takeRadiation<Loss>(planeRadiation(geometry, lut, measurements, n, previous), direction, mixture)) {
      pass.status = PropagationStatus::kNotFinite;
      pass.failedMeasurement = static_cast<std::int16_t>(previous);
      return;
    }
    if (step == n) {
      break;
    }
    const int index = previous + increment;
    const auto& stepProp = stepPropagator(prop, geometry, measurements, previous, index);
    int kept = 0;
    PropagationStatus failure = PropagationStatus::kOk;
    for (int k = 0; k < mixture.size; ++k) {
      Scalar logLikelihood = 0.F;
      const auto status = componentStep<Rule>(stepProp, geometry, material, measurements[index], direction, mixture.state[k], mixture.factors[k], mixture.chi2[k], mixture.reference[k], logLikelihood);
      if (status != PropagationStatus::kOk) {
        failure = status;
        continue;
      }
      mixture.logWeight[kept] = mixture.logWeight[k] + logLikelihood;
      mixture.state[kept] = mixture.state[k];
      mixture.factors[kept] = mixture.factors[k];
      mixture.chi2[kept] = mixture.chi2[k];
      mixture.reference[kept] = mixture.reference[k];
      ++kept;
    }
    mixture.size = kept;
    if (kept == 0) {
      pass.status = failure;
      pass.failedMeasurement = static_cast<std::int16_t>(index);
      return;
    }
    pass.minusTwoLogLikelihood -= 2.F * normalise(mixture);
  }
  normalise(mixture); // the radiation and the merges keep the sum, to rounding
  for (int k = 0; k < mixture.size; ++k) {
    covarianceFromFactors(mixture.factors[k], mixture.state[k].cov);
  }
  pass.collapsed = collapse(mixture);
  pass.qOverPMode = marginalMode<0>(mixture);
  pass.momentumMode = 1.F / o2::gpu::GPUCommonMath::Abs(marginalMode<2>(mixture));
}

} // namespace detail::gsf

/// The Gaussian-sum fit of an ordered chain, one measurement per logical
/// layer, under the propagator's mass hypothesis (an electron's): an outward
/// and an inward pass, each with nothing known but what its own first
/// measurements determine (`fitMeasurements`'s rule), the outward one
/// linearised from the chain's helix (`makeChainSeed`), the inward one from
/// the outward pass's collapsed end state, values only. Both take the
/// scattering of the chain's span as the Kalman fit does (`SpanMaterial`).
/// \a Loss the radiation model; at most \a K components after each
/// reduction, by default as many as one loss has (a plane's mixture keeps one
/// component per loss component; more changes nothing measurable).
template <typename Loss = BetheHeitlerLoss, int K = Loss::kComponents, CrossingRule Rule = CrossingRule::kFirstAhead, typename Prop, DetectorLayout Det, typename Material, typename Chain>
GPUdi() GSFFitResult<K> fitGaussianSum(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n)
{
  GSFFitResult<K> fit;
  const auto seed = makeChainSeed(prop, geometry, measurements, n);
  const auto span = detail::spanMaterial(geometry, material, measurements, n);
  detail::gsf::runPass<K, Loss, Rule>(prop, geometry, span, material, measurements, n, PropagationDirection::kForward, seed, fit.outward);
  if (fit.outward.status != PropagationStatus::kOk) {
    // The inward pass starts from the outward one's end: it does not run.
    fit.inward.status = fit.outward.status;
    fit.inward.failedMeasurement = fit.outward.failedMeasurement;
    return fit;
  }
  detail::gsf::runPass<K, Loss, Rule>(prop, geometry, span, material, measurements, n, PropagationDirection::kBackward, fit.outward.collapsed, fit.inward);
  return fit;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GSFTRACKFITTER_H_
