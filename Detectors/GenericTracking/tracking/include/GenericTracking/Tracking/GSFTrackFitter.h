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

#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Tracking/BetheHeitler.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Tracking/MeasurementUpdate.h"
#include "GenericTracking/Tracking/Seed.h"
#include "GenericTracking/Tracking/TrackFitterUtils.h"

#include "GPUCommonMath.h"
#include "MathUtils/detail/trigonometric.h"

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>
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

/// One Gaussian of a mixture: its weight (its logarithm), its state and its
/// covariance as \a Factors; `state.cov` is written only at a pass's end.
template <typename Factors>
struct GSFComponent {
  Scalar logWeight = 0.F;
  BoundTrackState state{};
  Factors factors{};
  Scalar chi2 = 0.F; ///< the summed chi2 increments
  /// The reference momentum the component's steps take the material at
  /// (`predictStep`), as the Kalman pass's (PassReference), scaled with each
  /// loss the component took.
  detail::PassReference reference{};
};

/// A mixture of at most \a Capacity components, all on one plane; a pass's
/// start holds diffuse ones (`DiffuseFactors`).
template <int Capacity, typename Factors = CovarianceFactors>
struct GSFMixture {
  static constexpr int kCapacity = Capacity;
  int size = 0;
  GSFComponent<Factors> components[Capacity]{};
};

/// The Gaussian-sum fit of a chain (`fitGaussianSum`): the outward and the
/// inward pass as the Kalman fit gives them (`fit`, TrackFitResult) -- the
/// mixtures' collapses as the states; as each pass's chi2 its mixture's
/// weighted mean of the components' chi2, the Kalman fit's without loss;
/// no kinks -- and what only a mixture has, outward ones prefixed as there:
/// the end mixtures (`state.cov` filled), the modes of the momentum's
/// density [GeV] (not the inverse of the q/p marginal's mode, the density
/// of p being the q/p one's times (q/p)^2), and each pass's -2 ln of the
/// measurements' likelihood under its predictions, the start's rows by the
/// diffuse likelihood (`startMixture`), without the constant 2 ln 2 pi per
/// coordinate pair, given that the components it kept propagated (a failed
/// one's weight is dropped).
template <int K>
struct GSFFitResult {
  TrackFitResult fit;
  GSFMixture<K> mixture;
  GSFMixture<K> outwardMixture;
  Scalar momentumMode = 0.F;
  Scalar outwardMomentumMode = 0.F;
  Scalar minusTwoLogLikelihood = 0.F;
  Scalar outwardMinusTwoLogLikelihood = 0.F;
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

/// A loss component on a diffuse covariance (`DiffuseFactors`): the q/p
/// scaling on both parts -- P_inf's root and P_star -- and the variance on
/// P_star, as `applyBetheHeitler` takes them on a full one.
GPUhdi() bool applyBetheHeitler(const BetheHeitlerComponent& component, PropagationDirection direction, BoundTrackState& state, DiffuseFactors& c)
{
  const Scalar scale = betheHeitlerScale(component, direction);
  for (int k = 0; k < c.rank; ++k) {
    c.root[kQOverP][k] *= scale;
  }
  return o2::gtrk::detail::applyBetheHeitlerFactors(component, direction, state, c.star);
}

/// A loss component on a mixture component: its state, its covariance and
/// its reference q/p scaled alike. False if the loss gave no finite state.
template <typename Factors>
GPUhdi() bool applyLoss(const BetheHeitlerComponent& loss, PropagationDirection direction, GSFComponent<Factors>& component)
{
  const Scalar scale = betheHeitlerScale(loss, direction);
  component.reference.qOverP *= scale;
  component.reference.lossScale *= scale;
  return applyBetheHeitler(loss, direction, component.state, component.factors);
}

/// Every component of \a from with each of the \a m loss components \a loss
/// into \a to: the weights multiplied (`applyLoss`). False if a loss gave no
/// finite state.
template <typename Factors, int N, int M>
GPUdi() bool branch(const GSFMixture<N, Factors>& from, const BetheHeitlerComponent* loss, int m, PropagationDirection direction, GSFMixture<M, Factors>& to)
{
  to.size = 0;
  for (int k = 0; k < from.size; ++k) {
    for (int j = 0; j < m; ++j) {
      auto& component = to.components[to.size++];
      component = from.components[k];
      component.logWeight += logOf(loss[j].weight);
      if (!applyLoss(loss[j], direction, component)) {
        return false;
      }
    }
  }
  return true;
}

/// \a step (a status of a component) on every component of \a mixture, the
/// ones it fails dropped and the rest kept in order: kOk if any is left,
/// else the last failure.
template <int N, typename Factors, typename Step>
GPUdi() PropagationStatus stepEach(GSFMixture<N, Factors>& mixture, const Step& step)
{
  int kept = 0;
  PropagationStatus failure = PropagationStatus::kOk;
  for (int k = 0; k < mixture.size; ++k) {
    const auto status = step(mixture.components[k]);
    if (status != PropagationStatus::kOk) {
      failure = status;
      continue;
    }
    mixture.components[kept++] = mixture.components[k];
  }
  mixture.size = kept;
  return kept > 0 ? PropagationStatus::kOk : failure;
}

/// Plane \a i's radiation [x/X0]: the link material's from the middle of
/// the step from the measurement before to the middle of the step to the one
/// after (LinkMaterial::step: the band reached and the second half of the
/// gap before, the band left and the first half of the gap after); nothing
/// beyond the chain's ends.
template <DetectorLayout Det, typename Chain>
GPUdi() Scalar planeRadiation(const TrackingGeometry<Det>& geometry, const Chain& measurements, int n, int i)
{
  Scalar point[3];
  const auto& plane = geometry.planeOf(measurements[i]);
  surfaceToGlobal(plane, measurements[i].loc0, measurements[i].loc1, point);
  Scalar x2X0 = 0.F;
  if (i > 0) {
    Scalar before[3];
    const auto& from = geometry.planeOf(measurements[i - 1]);
    surfaceToGlobal(from, measurements[i - 1].loc0, measurements[i - 1].loc1, before);
    const auto step = geometry.material().step(from.layer, plane.layer, before, point);
    x2X0 += step.after + (step.before - step.left - step.gapFirstHalf);
  }
  if (i + 1 < n) {
    Scalar after[3];
    const auto& to = geometry.planeOf(measurements[i + 1]);
    surfaceToGlobal(to, measurements[i + 1].loc0, measurements[i + 1].loc1, after);
    const auto step = geometry.material().step(plane.layer, to.layer, point, after);
    x2X0 += step.left + step.gapFirstHalf;
  }
  return x2X0;
}

/// How many equal parts a plane's radiation \a x2X0 is taken in: the fewest
/// no thicker than the model's range.
template <typename Loss>
GPUhdi() int lossParts(Scalar x2X0)
{
  int parts = 1;
  if (!isFiniteOf(x2X0)) {
    return parts; // the model refuses it
  }
  while (x2X0 > static_cast<Scalar>(parts) * Loss::kMaxX2X0) {
    ++parts;
  }
  return parts;
}

/// Components \a a and \a b as one Gaussian of weight \a wa + \a wb
/// (relative weights) into \a a: the weighted mean, phi taken relative to
/// a's, and the covariance (wa Pa + wb Pb) / W + wa wb / W^2 d d^T, d the
/// means' difference, its factors by weighted Gram-Schmidt over the factor
/// columns of both and d (`factorColumns`): no covariance is subtracted. The
/// chi2 and the reference momentum are the weighted means. Two weights that
/// are both zero (negligible beside the mixture's heaviest) leave \a a.
GPUdi() void mergeInto(Scalar wa, Scalar wb, GSFComponent<CovarianceFactors>& ca, const GSFComponent<CovarianceFactors>& cb)
{
  const Scalar total = wa + wb;
  if (!(total > 0.F)) {
    return;
  }
  auto& a = ca.state;
  auto& fa = ca.factors;
  const auto& b = cb.state;
  const auto& fb = cb.factors;
  const Scalar sa = wa / total;
  const Scalar sb = wb / total;
  Scalar difference[kNBoundParams];
  for (int i = 0; i < kNBoundParams; ++i) {
    difference[i] = a.par[i] - b.par[i];
  }
  o2::math_utils::detail::bringToPMPi<Scalar>(difference[kPhi]);
  constexpr int kDifference = 2 * kNBoundParams; // the column of d
  constexpr int kColumns = kDifference + 1;
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
    y[k][kDifference] = difference[k];
  }
  w[kDifference] = sa * sb;
  factorColumns(y, w, fa);
  // a's mean moves towards b's by b's share of the difference.
  for (int i = 0; i < kNBoundParams; ++i) {
    a.par[i] -= sb * difference[i];
  }
  o2::math_utils::detail::bringToPMPi<Scalar>(a.par[kPhi]);
  ca.chi2 = (sa * ca.chi2) + (sb * cb.chi2);
  ca.reference.qOverP = (sa * ca.reference.qOverP) + (sb * cb.reference.qOverP);
  ca.reference.lossScale = (sa * ca.reference.lossScale) + (sb * cb.reference.lossScale);
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
  return 0.5F * ((total * logOf(merged)) - (wa * logSa) - (wb * logSb));
}

/// \a from's components reduced to at most \a K into \a to: until K remain,
/// the pair of least `mergeCost` merges (`mergeInto`). The pair costs sit in
/// a triangular table over \a from's slots; a merge recomputes the merged
/// component's costs and retires the other's.
template <int K, int N>
GPUdi() void reduce(GSFMixture<N>& from, GSFMixture<K>& to)
{
  auto& c = from.components;
  Scalar maxLog = -FLT_MAX;
  for (int k = 0; k < from.size; ++k) {
    maxLog = o2::gpu::GPUCommonMath::Max(maxLog, c[k].logWeight);
  }
  Scalar weight[N];
  Scalar logVariance[N];
  bool alive[N];
  bool merged[N];
  for (int k = 0; k < from.size; ++k) {
    weight[k] = expOf(c[k].logWeight - maxLog);
    logVariance[k] = logOf(c[k].factors.diagonal[kQOverP]);
    alive[k] = true;
    merged[k] = false;
  }
  int remaining = from.size;
  if (remaining > K) {
    constexpr int kPairs = N * (N - 1) / 2;
    Scalar cost[kPairs > 0 ? kPairs : 1];
    const auto pairCost = [&](int a, int b) {
      return mergeCost(weight[a], c[a].state.qOverP(), c[a].factors.diagonal[kQOverP], logVariance[a],
                       weight[b], c[b].state.qOverP(), c[b].factors.diagonal[kQOverP], logVariance[b]);
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
      mergeInto(weight[bestA], weight[bestB], c[bestA], c[bestB]);
      weight[bestA] += weight[bestB];
      merged[bestA] = true;
      logVariance[bestA] = logOf(c[bestA].factors.diagonal[kQOverP]);
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
    auto& component = to.components[to.size++];
    component = c[k];
    // An unmerged component keeps its weight's logarithm as it was, however small.
    if (merged[k]) {
      component.logWeight = maxLog + logOf(weight[k]);
    }
  }
}

/// Every component of \a mixture with every loss component of the radiation
/// \a x2X0, reduced back to at most K (`reduce`); the radiation in equal
/// parts if it is thicker than the model's range. False if a component's
/// loss gave no finite state.
template <typename Loss, int K>
GPUdi() bool takeRadiation(Scalar x2X0, PropagationDirection direction, GSFMixture<K>& mixture)
{
  const int parts = lossParts<Loss>(x2X0);
  BetheHeitlerComponent loss[Loss::kComponents];
  const int m = Loss::mixture(x2X0 / static_cast<Scalar>(parts), loss);
  if (m <= 0) {
    return false;
  }
  for (int part = 0; part < parts; ++part) {
    GSFMixture<K * Loss::kComponents> branched;
    if (!branch(mixture, loss, m, direction, branched)) {
      return false;
    }
    reduce(branched, mixture);
  }
  return true;
}

/// The weights normalised to sum to one; returns the logarithm of what
/// they summed to.
template <typename Mixture>
GPUhdi() Scalar normalise(Mixture& mixture)
{
  Scalar maxLog = -FLT_MAX;
  for (int k = 0; k < mixture.size; ++k) {
    maxLog = o2::gpu::GPUCommonMath::Max(maxLog, mixture.components[k].logWeight);
  }
  Scalar sum = 0.F;
  for (int k = 0; k < mixture.size; ++k) {
    sum += expOf(mixture.components[k].logWeight - maxLog);
  }
  const Scalar logSum = maxLog + logOf(sum);
  for (int k = 0; k < mixture.size; ++k) {
    mixture.components[k].logWeight -= logSum;
  }
  return logSum;
}

GPUhdi() constexpr int power(int base, int exponent)
{
  int result = 1;
  for (int i = 0; i < exponent; ++i) {
    result *= base;
  }
  return result;
}

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
/// make differ, a tested one ln F + v^2 / F (`diffuseStepBody`). A start
/// longer than `kPassStart` rows (a coordinate that determined nothing)
/// takes the loss of the planes beyond as one Gaussian with the mixture's
/// moments. Ends reduced to at most K components on the last plane it took,
/// before that plane's radiation; \a minusTwoLogLikelihood gets the start's
/// -2 ln likelihood.
template <typename Loss, typename Prop, DetectorLayout Det, typename Material, typename Chain, int K>
GPUdi() PropagationStatus startMixture(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, PropagationDirection direction, BoundTrackState seed, GSFMixture<K>& mixture, int& taken, std::int16_t& failedMeasurement, Scalar& minusTwoLogLikelihood)
{
  constexpr int kBranches = kStartBranches<Loss>;
  using Branches = GSFMixture<kBranches, DiffuseFactors>;
  const int first = direction == PropagationDirection::kForward ? 0 : n - 1;
  const int increment = direction == PropagationDirection::kForward ? 1 : -1;
  Branches branches;
  branches.size = 1;
  auto& start = branches.components[0];
  start.state = seed;
  start.reference = PassReference{.qOverP = seed.qOverP()}; // as `startPass`'s
  taken = 0;
  if (n < kPassStart || !startDiffuse(measurements[first], start.state, start.factors)) {
    failedMeasurement = static_cast<std::int16_t>(first);
    return PropagationStatus::kNotFinite;
  }
  const auto undetermined = [&]() {
    bool any = false;
    for (int b = 0; b < branches.size; ++b) {
      any = any || branches.components[b].factors.rank > 0;
    }
    return any;
  };
  for (taken = 1; taken < n && undetermined(); ++taken) {
    const int previous = first + ((taken - 1) * increment);
    const int index = previous + increment;
    const Scalar x2X0 = planeRadiation(geometry, measurements, n, previous);
    const int parts = lossParts<Loss>(x2X0);
    BetheHeitlerComponent loss[Loss::kComponents];
    const int m = Loss::mixture(x2X0 / static_cast<Scalar>(parts), loss);
    if (m <= 0) {
      failedMeasurement = static_cast<std::int16_t>(previous);
      return PropagationStatus::kNotFinite;
    }
    for (int part = 0; part < parts; ++part) {
      bool finite = true;
      if (branches.size * m <= kBranches) {
        Branches branched;
        finite = branch(branches, loss, m, direction, branched);
        branches = branched;
      } else {
        const auto moments = momentsOf(loss, m);
        for (int b = 0; b < branches.size && finite; ++b) {
          finite = applyLoss(moments, direction, branches.components[b]);
        }
      }
      if (!finite) {
        failedMeasurement = static_cast<std::int16_t>(previous);
        return PropagationStatus::kNotFinite;
      }
    }
    const auto& stepProp = stepPropagator(prop, geometry, measurements, previous, index);
    const Scalar leftOut = leftOutQOverP(prop, geometry, measurements, n, index, direction);
    const auto status = stepEach(branches, [&](GSFComponent<DiffuseFactors>& component) {
      Scalar density = 0.F;
      const Scalar stepReference = referenceFor(component.reference, n, direction, [&] { return leftOut; });
      Scalar carried = stepReference;
      const auto stepped = o2::gtrk::detail::diffuseStepBody<true>(stepProp, geometry, material, measurements[index], direction, component.state, component.factors, component.chi2, carried, &density);
      advanceReference(component.reference, n, direction, stepReference, carried);
      component.logWeight -= 0.5F * density;
      return stepped;
    });
    if (status != PropagationStatus::kOk) {
      failedMeasurement = static_cast<std::int16_t>(index);
      return status;
    }
    minusTwoLogLikelihood -= 2.F * normalise(branches);
  }
  if (undetermined()) {
    failedMeasurement = static_cast<std::int16_t>(first + ((n - 1) * increment));
    return PropagationStatus::kNotFinite;
  }
  GSFMixture<kBranches> determined;
  determined.size = branches.size;
  for (int b = 0; b < branches.size; ++b) {
    const auto& from = branches.components[b];
    determined.components[b] = {.logWeight = from.logWeight, .state = from.state, .factors = from.factors.star, .chi2 = from.chi2, .reference = from.reference};
  }
  reduce(determined, mixture);
  return PropagationStatus::kOk;
}

/// The mixture as one Gaussian: the weighted mean (phi about the heaviest
/// component's, wrapped) and the covariance with the means' spread.
template <int K>
GPUhdi() BoundTrackState collapse(const GSFMixture<K>& mixture)
{
  const auto& c = mixture.components;
  int heaviest = 0;
  for (int k = 1; k < mixture.size; ++k) {
    if (c[k].logWeight > c[heaviest].logWeight) {
      heaviest = k;
    }
  }
  BoundTrackState out = c[heaviest].state;
  Scalar mean[kNBoundParams]{};
  Scalar weights[K];
  for (int k = 0; k < mixture.size; ++k) {
    weights[k] = expOf(c[k].logWeight);
    for (int i = 0; i < kNBoundParams; ++i) {
      Scalar value = c[k].state.par[i] - out.par[i];
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
      d[i] = c[k].state.par[i] - out.par[i];
      if (i == kPhi) {
        o2::math_utils::detail::bringToPMPi<Scalar>(d[i]);
      }
      d[i] -= mean[i];
    }
    for (int i = 0; i < kNBoundParams; ++i) {
      for (int j = 0; j <= i; ++j) {
        cov[covIndex(i, j)] += weights[k] * (c[k].state.cov[covIndex(i, j)] + (d[i] * d[j]));
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

/// The mode [GeV] of the momentum's density, in u = q/p sum_k w_k N(u; m_k,
/// s_k) u^2 (the q/p marginal's times the Jacobian (q/p)^2): from each
/// component's mean, Newton's steps on the log density, a step that does not
/// raise the density halved, until the step no longer moves u; the highest
/// of the ends.
template <int K>
GPUhdi() Scalar momentumModeOf(const GSFMixture<K>& mixture)
{
  const auto& c = mixture.components;
  const auto logDensity = [&](Scalar u, Scalar& gradient, Scalar& curvature) {
    Scalar maxLog = -FLT_MAX;
    Scalar logTerm[K];
    for (int k = 0; k < mixture.size; ++k) {
      const Scalar s = c[k].state.cov[covIndex(kQOverP, kQOverP)];
      const Scalar d = u - c[k].state.qOverP();
      logTerm[k] = c[k].logWeight - (0.5F * logOf(s)) - (0.5F * d * d / s);
      maxLog = o2::gpu::GPUCommonMath::Max(maxLog, logTerm[k]);
    }
    Scalar f = 0.F;
    Scalar f1 = 0.F;
    Scalar f2 = 0.F;
    for (int k = 0; k < mixture.size; ++k) {
      const Scalar s = c[k].state.cov[covIndex(kQOverP, kQOverP)];
      const Scalar d = u - c[k].state.qOverP();
      const Scalar term = expOf(logTerm[k] - maxLog);
      f += term;
      f1 -= term * d / s;
      f2 += term * (((d * d) / (s * s)) - (1.F / s));
    }
    // ln g = ln f + 2 ln |u|
    const Scalar g1 = (f1 / f) + (2.F / u);
    const Scalar g2 = (f2 / f) - ((f1 * f1) / (f * f)) - (2.F / (u * u));
    gradient = g1;
    curvature = g2;
    return maxLog + logOf(f) + (2.F * logOf(o2::gpu::GPUCommonMath::Abs(u)));
  };
  Scalar best = c[0].state.qOverP();
  Scalar bestValue = -FLT_MAX;
  for (int k = 0; k < mixture.size; ++k) {
    Scalar u = c[k].state.qOverP();
    Scalar gradient;
    Scalar curvature;
    Scalar value = logDensity(u, gradient, curvature);
    // Each accepted step raises the density; a step is at least halved until it does (a bounded number of
    // halvings leaves u unchanged in float), and the walk ends when a step leaves u unchanged. Newton's steps
    // double the digits they resolve near the mode and halved steps halve the distance to it, so float's
    // digits bound the steps that still move u; the walk is held to that many.
    for (int iteration = 0; iteration < std::numeric_limits<Scalar>::digits; ++iteration) {
      Scalar step = curvature < 0.F ? -gradient / curvature : gradient * c[k].state.cov[covIndex(kQOverP, kQOverP)];
      if (!isFiniteOf(step)) {
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
  return 1.F / o2::gpu::GPUCommonMath::Abs(best);
}

/// A mixture's end on a plane: the weights normalised (the radiation and the
/// merges keep their sum, to rounding), each component's `state.cov`
/// written, the mixture's collapse into \a collapsed.
template <int K>
GPUhdi() void collapseMixture(GSFMixture<K>& mixture, BoundTrackState& collapsed)
{
  normalise(mixture);
  for (int k = 0; k < mixture.size; ++k) {
    covarianceFromFactors(mixture.components[k].factors, mixture.components[k].state.cov);
  }
  collapsed = collapse(mixture);
}

/// `collapseMixture`, and the mode of the momentum's density [GeV] into
/// \a momentumMode.
template <int K>
GPUhdi() void finishMixture(GSFMixture<K>& mixture, BoundTrackState& collapsed, Scalar& momentumMode)
{
  collapseMixture(mixture, collapsed);
  momentumMode = momentumModeOf(mixture);
}

/// One pass of the filter over the chain from \a seed: the start
/// (`startMixture`), then on every plane its radiation (`takeRadiation`) and
/// each component's Kalman step to the next (`filterStepBody`) with its
/// weight's likelihood; at the end the last plane's radiation and the
/// summaries (`finishMixture`, `GSFFitResult`) into \a mixture, \a end,
/// \a chi2, \a momentumMode and \a minusTwoLogLikelihood. On a failure,
/// \a failedMeasurement is the measurement the failed step went to.
template <int K, typename Loss, typename Prop, DetectorLayout Det, typename Material, typename Chain>
GPUdi() PropagationStatus runPass(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, PropagationDirection direction, const BoundTrackState& seed, GSFMixture<K>& mixture, BoundTrackState& end, Scalar& chi2, Scalar& momentumMode, Scalar& minusTwoLogLikelihood, std::int16_t& failedMeasurement)
{
  const int first = direction == PropagationDirection::kForward ? 0 : n - 1;
  const int increment = direction == PropagationDirection::kForward ? 1 : -1;
  int taken = 0;
  minusTwoLogLikelihood = 0.F;
  const auto started = startMixture<Loss>(prop, geometry, material, measurements, n, direction, seed, mixture, taken, failedMeasurement, minusTwoLogLikelihood);
  if (started != PropagationStatus::kOk) {
    return started;
  }
  for (int step = taken; step <= n; ++step) {
    const int previous = first + ((step - 1) * increment);
    if (!takeRadiation<Loss>(planeRadiation(geometry, measurements, n, previous), direction, mixture)) {
      failedMeasurement = static_cast<std::int16_t>(previous);
      return PropagationStatus::kNotFinite;
    }
    if (step == n) {
      break;
    }
    const int index = previous + increment;
    const auto& stepProp = stepPropagator(prop, geometry, measurements, previous, index);
    const Scalar leftOut = leftOutQOverP(prop, geometry, measurements, n, index, direction);
    const auto status = stepEach(mixture, [&](GSFComponent<CovarianceFactors>& component) {
      Scalar logLikelihood = 0.F;
      const Scalar stepReference = referenceFor(component.reference, n, direction, [&] { return leftOut; });
      Scalar carried = stepReference;
      const auto stepped = o2::gtrk::detail::filterStepBody<false, true>(stepProp, geometry, material, measurements[index], direction, component.state, component.factors, component.chi2, carried, 1.F, &logLikelihood);
      advanceReference(component.reference, n, direction, stepReference, carried);
      component.logWeight += logLikelihood;
      return stepped;
    });
    if (status != PropagationStatus::kOk) {
      failedMeasurement = static_cast<std::int16_t>(index);
      return status;
    }
    minusTwoLogLikelihood -= 2.F * normalise(mixture);
  }
  finishMixture(mixture, end, momentumMode);
  chi2 = 0.F;
  for (int k = 0; k < mixture.size; ++k) {
    chi2 += expOf(mixture.components[k].logWeight) * mixture.components[k].chi2;
  }
  return PropagationStatus::kOk;
}

} // namespace detail::gsf

/// The Gaussian-sum fit of an ordered chain, one measurement per logical
/// layer, under the propagator's mass hypothesis (an electron's): an outward
/// and an inward pass, each with nothing known but what its own first
/// measurements determine (`fitMeasurements`'s rule), the outward one
/// linearised from the chain's helix (`makeChainSeed`), the inward one from
/// the outward pass's collapsed end state, values only. Both take the
/// scattering of the chain's span as the Kalman fit does (`SpanMaterial`).
/// If the outward pass fails, the inward one does not run (`failedPass`).
/// \a Loss the radiation model; at most \a K components after each
/// reduction, by default as many as one loss has (a plane's mixture keeps one
/// component per loss component; more changes nothing measurable).
template <typename Loss = BetheHeitlerLoss, int K = Loss::kComponents, typename Prop, DetectorLayout Det, typename Material, typename Chain>
GPUdi() GSFFitResult<K> fitGaussianSum(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n)
{
  GSFFitResult<K> result;
  auto& fit = result.fit;
  fit.ndf = static_cast<std::int16_t>((2 * n) - kNBoundParams);
  const auto seed = makeChainSeed(prop, geometry, measurements, n);
  const auto span = detail::spanMaterial(geometry, material, measurements, n);
  auto status = detail::gsf::runPass<K, Loss>(prop, geometry, span, measurements, n, PropagationDirection::kForward, seed, result.outwardMixture, fit.outwardState, fit.outwardChi2, result.outwardMomentumMode, result.outwardMinusTwoLogLikelihood, fit.failedMeasurement);
  if (status != PropagationStatus::kOk) {
    fit.status = status;
    fit.failedPass = FitPass::kOutward;
    return result;
  }
  status = detail::gsf::runPass<K, Loss>(prop, geometry, span, measurements, n, PropagationDirection::kBackward, fit.outwardState, result.mixture, fit.state, fit.chi2, result.momentumMode, result.minusTwoLogLikelihood, fit.failedMeasurement);
  if (status != PropagationStatus::kOk) {
    fit.status = status;
    fit.failedPass = FitPass::kInward;
  }
  return result;
}

#ifndef GPUCA_GPUCODE_DEVICE
/// The q/p marginal's distribution function at \a u, from the components'
/// `state.cov` (a pass's end): the probability integral transform of a true
/// q/p.
template <int K>
inline double marginalCdf(const GSFMixture<K>& mixture, double u)
{
  double sum = 0.;
  for (int k = 0; k < mixture.size; ++k) {
    const auto& c = mixture.components[k];
    const double sigma = std::sqrt(static_cast<double>(c.state.variance(kQOverP)));
    sum += std::exp(static_cast<double>(c.logWeight)) * 0.5 * std::erfc(-(u - static_cast<double>(c.state.qOverP())) / (sigma * std::numbers::sqrt2));
  }
  return sum;
}
#endif

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GSFTRACKFITTER_H_
