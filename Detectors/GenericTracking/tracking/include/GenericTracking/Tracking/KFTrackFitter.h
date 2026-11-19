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

/// \file KFTrackFitter.h
/// \brief Independent outward and inward Kalman filters over a chain, each step linearised about its filtered state and the pass's reference momentum.

#ifndef GTRK_TRACKING_KFTRACKFITTER_H_
#define GTRK_TRACKING_KFTRACKFITTER_H_

#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Seed.h"
#include "GenericTracking/Tracking/TrackFitterUtils.h"
#include "GenericTracking/Tracking/WorkCounters.h"

#include "GPUCommonMath.h"

#include <cstdint>

namespace o2::gtrk
{

/// How many steps of a filter pass may be kinks: a
/// step whose spatial increment is beyond \a gate, with a measurement after
/// it in the pass, is taken again from the surface before with the track's
/// direction freed (`freeDirection`) -- a break point, two degrees of
/// freedom. None by default.
/// The last step, with nothing after it to check a freed direction, never
/// is; the road fit's takes the kink allowance there instead
/// (`fitOutward`'s LastStepAllowance).
struct KinkBudget {
  int max = 0;
  Scalar gate = 0.F; ///< `chiSquareQuantile(2, nSigma)`
};

/// An iteration's kink budget, from its parameters' `maxKinks` and
/// `kinkChi2Min`: the refit's and the extension's.
template <typename Parameters>
[[nodiscard]] GPUhdi() GPUAlwaysInline() KinkBudget kinkBudget(const Parameters& parameters)
{
  return KinkBudget{.max = static_cast<int>(parameters.maxKinks), .gate = parameters.kinkChi2Min};
}

/// Whether an iteration's road fit takes its kink allowance on the last step
/// (`fitOutward`'s LastStepAllowance, against the road-end fakes): every iteration with an allowance. The one rule the host's and the
/// device's road-fit stages are instantiated by and `fitRoad` follows.
template <typename Parameters>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool takesLastStepAllowance(const Parameters& parameters)
{
  return parameters.kinkScale2 > 1.F;
}

/// A fit's per-measurement chi2 increments, for QA: called once per update
/// with the pass's direction, the measurement's index in the chain and its
/// increment (two dof). The default does nothing.
struct NoFitObserver {
  GPUhdi() GPUAlwaysInline() void operator()(PropagationDirection /*direction*/, int /*index*/, Scalar /*chi2*/) const {}
};

/// An observer that also takes each pass's filtered state on every
/// measurement once the pass has determined it -- from the start's last
/// measurement on (`startPass`) -- with its covariance's factors:
/// `observer.state(direction, index, state, factors)`. Only such an observer
/// compiles the call; `NoFitObserver` does not.
template <typename Observer>
concept ObservesFitStates = requires(const Observer& observer, PropagationDirection direction, int index, const BoundParameters& state, const CovarianceFactors& factors) {
  observer.state(direction, index, state, factors);
};

namespace detail
{

/// One step of a filter pass: \a state predicted on \a measurement's plane
/// (`predictStep`) and updated with it. The covariance is \a factors, the pass's
/// (CovarianceFactors): `state.cov` is neither read nor written. The pass's
/// reference momentum \a referenceQOverP (`predictStep`) moves to the plane
/// reached only with an accepted update. \a Scaled,
/// the step's scattering theta0 is multiplied by \a noiseScale, both halves:
/// the kink allowance's last step (`allowanceStep`); the nominal step is its
/// own instantiation. \a Rule is the propagation's choice of crossing
/// (`CrossingRule`), and with it every pass below. Never inlined on the
/// device: one copy for each instantiation, not one in every fit of every
/// stage variant.
template <bool Scaled = false, CrossingRule Rule = CrossingRule::kFirstAhead, typename Solve = ExactScalar, DetectorLayout Det, typename Material, typename Prop>
GPUdnii() PropagationStatus filterStep(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& GPUrestrict() state, CovarianceFactors& GPUrestrict() factors, Scalar& GPUrestrict() chi2, Scalar& GPUrestrict() referenceQOverP, [[maybe_unused]] Scalar noiseScale = 1.F)
{
  return filterStepBody<Scaled, Rule, false, Solve>(prop, geometry, material, measurement, direction, state, factors, chi2, referenceQOverP, noiseScale, nullptr);
}

/// A filter step (`filterStep`) on a diffuse covariance \a c: the prediction
/// the same, the transport `transportDiffuse`, the update `updateDiffuse`;
/// \a chi2 gets only what the update tests. Once P_inf is gone, \a c's
/// factors are the pass's. Never inlined on the device.
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Solve = ExactScalar, DetectorLayout Det, typename Material, typename Prop>
GPUdnii() PropagationStatus diffuseStep(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& GPUrestrict() state, DiffuseFactors& GPUrestrict() c, Scalar& GPUrestrict() chi2, Scalar& GPUrestrict() referenceQOverP)
{
  return diffuseStepBody<Rule, false, Solve>(prop, geometry, material, measurement, direction, state, c, chi2, referenceQOverP, nullptr);
}

/// \a budget for a pass over \a n measurements: no more kinks than leave the
/// fit a degree of freedom, `2 n - 5 - 2 k >= 1`, so `n - 3` -- one for four
/// measurements, two for five. Each kink frees the direction, two parameters;
/// beyond that the chain no longer constrains the track, and the chi2 cut by
/// degrees of freedom has no entry for it.
GPUhdi() GPUAlwaysInline() KinkBudget kinkBudgetFor(KinkBudget budget, int n)
{
  const int allowed = ((2 * n) - kNBoundParams - 1) / 2;
  if (budget.max > allowed) {
    budget.max = allowed > 0 ? allowed : 0;
  }
  return budget;
}

/// A kink's freed direction, as the exact diffuse filter (DiffuseFactors)
/// takes it: \a factors' covariance with phi and theta's rows and columns
/// zero as P_star -- position and q/p keep theirs -- and the two angles as
/// P_inf.
GPUdi() void freeDirection(const CovarianceFactors& factors, DiffuseFactors& c)
{
  Scalar y[kNBoundParams][kNBoundParams]{};
  Scalar w[kNBoundParams]{};
  factorColumnsOf(factors, y, w);
  for (int k = 0; k < kNBoundParams; ++k) {
    y[kPhi][k] = 0.F;
    y[kTheta][k] = 0.F;
  }
  factorColumns(y, w, c.star);
  c.root[kPhi][0] = 1.F;
  c.root[kTheta][1] = 1.F;
  c.rank = 2;
}

/// A filter step that may be a kink: taken as `filterStep`; if kinks are left,
/// \a followed -- a measurement comes after this one in the pass -- and its
/// increment is beyond the gate, taken again from the surface before with the
/// direction freed (`freeDirection`), and counted in \a kinks: the
/// measurement then determines the direction and adds nothing to the chi2.
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Solve = ExactScalar, DetectorLayout Det, typename Material, typename Prop>
GPUdi() PropagationStatus kinkStep(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& state, CovarianceFactors& factors, Scalar& chi2, Scalar& referenceQOverP, const KinkBudget& budget, bool followed, int& kinks)
{
  if (kinks >= budget.max || !followed) {
    return filterStep<false, Rule, Solve>(prop, geometry, material, measurement, direction, state, factors, chi2, referenceQOverP);
  }
  const BoundTrackState keptState = state;
  const CovarianceFactors keptFactors = factors;
  const Scalar keptChi2 = chi2;
  const Scalar keptReference = referenceQOverP;
  auto status = filterStep<false, Rule, Solve>(prop, geometry, material, measurement, direction, state, factors, chi2, referenceQOverP);
  if (status != PropagationStatus::kOk || !(chi2 - keptChi2 > budget.gate)) {
    return status;
  }
  state = keptState;
  chi2 = keptChi2;
  referenceQOverP = keptReference;
  DiffuseFactors freed;
  freeDirection(keptFactors, freed);
  status = diffuseStep<Rule, Solve>(prop, geometry, material, measurement, direction, state, freed, chi2, referenceQOverP);
  if (status == PropagationStatus::kOk && freed.rank != 0) {
    status = PropagationStatus::kNotFinite;
  }
  if (status == PropagationStatus::kOk) {
    factors = freed.star;
    ++kinks;
  }
  return status;
}

/// The road fit's last step with the kink allowance (`fitOutward`'s
/// LastStepAllowance, against the road-end fakes): taken as
/// `filterStep`; if kinks are left and its increment is beyond the gate,
/// taken again from the surface before with its scattering's theta0 scaled by
/// \a noiseScale, as the extension's kink step, and counted in \a kinks;
/// \a refused if the increment is still beyond the gate. It frees no
/// parameter.
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Solve = ExactScalar, DetectorLayout Det, typename Material>
GPUdi() PropagationStatus allowanceStep(const Propagator& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& state, CovarianceFactors& factors, Scalar& chi2, Scalar& referenceQOverP, const KinkBudget& budget, Scalar noiseScale, int& kinks, bool& refused)
{
  refused = false;
  if (kinks >= budget.max) {
    return filterStep<false, Rule, Solve>(prop, geometry, material, measurement, direction, state, factors, chi2, referenceQOverP);
  }
  const BoundTrackState keptState = state;
  const CovarianceFactors keptFactors = factors;
  const Scalar keptChi2 = chi2;
  const Scalar keptReference = referenceQOverP;
  auto status = filterStep<false, Rule, Solve>(prop, geometry, material, measurement, direction, state, factors, chi2, referenceQOverP);
  if (status != PropagationStatus::kOk || !(chi2 - keptChi2 > budget.gate)) {
    return status;
  }
  state = keptState;
  factors = keptFactors;
  chi2 = keptChi2;
  referenceQOverP = keptReference;
  status = filterStep<true, Rule, Solve>(prop, geometry, material, measurement, direction, state, factors, chi2, referenceQOverP, noiseScale);
  if (status == PropagationStatus::kOk) {
    ++kinks;
    refused = !(chi2 - keptChi2 <= budget.gate);
  }
  return status;
}

/// A pass's start, with nothing known but what its measurements determine
/// (`startDiffuse`, then `diffuseStep`) until nothing is undetermined: the
/// first fixes the position, the next two the direction and q/p, the third's
/// second coordinate the first the chi2 tests -- three measurements, unless a
/// coordinate carries no undetermined direction. Ends with \a state, its
/// \a factors and \a chi2 on the last measurement it took, their count in
/// \a taken, and the pass's \a reference momentum (PassReference) there,
/// from \a state's q/p as the pass begins. \a observe sees each
/// measurement's increment. A pass whose measurements leave a direction
/// undetermined determines no state.
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Solve = ExactScalar, typename Prop, DetectorLayout Det, typename Material, typename Chain, typename Observer>
GPUdi() PropagationStatus startPass(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, PropagationDirection direction, BoundTrackState& state, CovarianceFactors& factors, Scalar& chi2, PassReference& reference, int& taken, std::int16_t& failedMeasurement, const Observer& observe)
{
  const int first = direction == PropagationDirection::kForward ? 0 : n - 1;
  const int increment = direction == PropagationDirection::kForward ? 1 : -1;
  DiffuseFactors c;
  taken = 0;
  reference = PassReference{.qOverP = state.qOverP()};
  if (n < kPassStart || !startDiffuse(measurements[first], state, c)) {
    failedMeasurement = static_cast<std::int16_t>(first);
    return PropagationStatus::kNotFinite;
  }
  chi2 = 0.F;
  observe(direction, first, 0.F);
  for (taken = 1; taken < n && c.rank > 0; ++taken) {
    const int index = first + (taken * increment);
    const Scalar before = chi2;
    const auto& stepProp = stepPropagator(prop, geometry, measurements, index - increment, index);
    const Scalar stepReference = referenceFor(reference, n, direction, [&] { return leaveOneOutQOverP(prop, geometry, measurements, n, index); });
    Scalar carried = stepReference;
    const auto status = diffuseStep<Rule, Solve>(stepProp, geometry, material, measurements[index], direction, state, c, chi2, carried);
    if (status != PropagationStatus::kOk) {
      failedMeasurement = static_cast<std::int16_t>(index);
      return status;
    }
    advanceReference(reference, n, direction, stepReference, carried);
    observe(direction, index, chi2 - before);
  }
  if (c.rank != 0) {
    failedMeasurement = static_cast<std::int16_t>(first + ((n - 1) * increment));
    return PropagationStatus::kNotFinite;
  }
  factors = c.star;
  if constexpr (ObservesFitStates<Observer>) {
    observe.state(direction, first + ((taken - 1) * increment), state, factors);
  }
  return PropagationStatus::kOk;
}

/// A filter pass from \a seed over the chain: its start (`startPass`), then
/// the filter. Its covariance is carried as factors from the first update to
/// the last (CovarianceFactors) and written to \a result's `cov` once, at
/// the end. \a seed's values are the point its first steps linearise about,
/// its q/p the pass's reference momentum (`startPass`); a kink is not tested
/// before the start has determined the state.
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Solve = ExactScalar, typename Prop, DetectorLayout Det, typename Material, typename Chain, typename Observer>
GPUdi() PropagationStatus runPass(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, PropagationDirection direction, BoundTrackState seed, BoundTrackState& result, Scalar& chi2, std::int16_t& failedMeasurement, const KinkBudget& passBudget, int& kinks, std::int16_t& kinkMeasurement, const Observer& observe)
{
  const auto budget = kinkBudgetFor(passBudget, n);
  const int first = direction == PropagationDirection::kForward ? 0 : n - 1;
  const int increment = direction == PropagationDirection::kForward ? 1 : -1;

  BoundTrackState state = seed;
  CovarianceFactors factors;
  PassReference reference;
  int taken = 0;
  const auto started = startPass<Rule, Solve>(prop, geometry, material, measurements, n, direction, state, factors, chi2, reference, taken, failedMeasurement, observe);
  if (started != PropagationStatus::kOk) {
    return started;
  }
  for (int step = taken; step < n; ++step) {
    const int index = first + step * increment;
    const Scalar before = chi2;
    const int kinksBefore = kinks;
    const auto& stepProp = stepPropagator(prop, geometry, measurements, index - increment, index);
    const Scalar stepReference = referenceFor(reference, n, direction, [&] { return leaveOneOutQOverP(prop, geometry, measurements, n, index); });
    Scalar carried = stepReference;
    const auto status = kinkStep<Rule, Solve>(stepProp, geometry, material, measurements[index], direction, state, factors, chi2, carried, budget, step < n - 1, kinks);
    if (status != PropagationStatus::kOk) {
      failedMeasurement = static_cast<std::int16_t>(index);
      return status;
    }
    advanceReference(reference, n, direction, stepReference, carried);
    if (kinks != kinksBefore) {
      kinkMeasurement = static_cast<std::int16_t>(index);
    }
    observe(direction, index, chi2 - before);
    if constexpr (ObservesFitStates<Observer>) {
      observe.state(direction, index, state, factors);
    }
  }
  covarianceFromFactors(factors, state.cov);
  result = state;
  return PropagationStatus::kOk;
}

/// The inward pass of `fitMeasurements`: from the outward pass's end state
/// (the whole-chain helix carried there unfitted misses the outermost rows at
/// low momentum); its values only, so no measurement counts twice.
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Solve = ExactScalar, typename Prop, DetectorLayout Det, typename Material, typename Chain, typename Observer>
GPUdi() void inwardPass(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, const KinkBudget& kinkBudget, const Observer& observe, TrackFitResult& fit)
{
  int kinks = 0;
  const auto status = runPass<Rule, Solve>(prop, geometry, material, measurements, n, PropagationDirection::kBackward, fit.outwardState, fit.state, fit.chi2, fit.failedMeasurement, kinkBudget, kinks, fit.kinkMeasurement, observe);
  if (status != PropagationStatus::kOk) {
    fit.status = status;
    fit.failedPass = FitPass::kInward;
    return;
  }
  fit.kinks = static_cast<std::uint8_t>(kinks);
  fit.ndf = static_cast<std::int16_t>(fit.ndf - (2 * kinks));
}

} // namespace detail

/// Fit an ordered chain containing one measurement per logical layer.
///
/// Neither pass takes a prior: each starts with what its own first
/// measurements determine (`startPass`). The outward one starts from a
/// helix fitted through every measurement of the chain (`makeChainSeed`),
/// the inward one from the outward filter's end state; each step linearises
/// about the filtered state it starts from (`filterStep`), its material at the
/// pass's reference momentum, the start state's carried along
/// (`predictStep`). The inward filter takes no covariance from the outward
/// one, so no measurement's weight is counted twice in the returned state.
/// Both take the scattering of the chain's span
/// as a whole (`SpanMaterial`), as `fitInward` and `fitOutward` do.
///
/// The measurements' times are boxes, which the fit does not take: the
/// track's time is their intersection, set by the refit (TrackRefit.h).
///
/// \a measurements is anything indexed by `[i]` giving a `Measurement`, as
/// for `makeChainSeed`. Each pass may take \a kinkBudget's kinks. \a observe
/// sees every update's chi2 increments (`NoFitObserver`).
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Prop, DetectorLayout Det, typename Material, typename Chain, typename Observer = NoFitObserver>
GPUdi() TrackFitResult fitMeasurements(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, const KinkBudget& kinkBudget = {}, const Observer& observe = {})
{
  TrackFitResult fit{.ndf = static_cast<std::int16_t>((2 * n) - kNBoundParams)};
  const auto seed = makeChainSeed(prop, geometry, measurements, n);
  const auto span = detail::spanMaterial(geometry, material, measurements, n);
  int outwardKinks = 0;
  std::int16_t outwardKinkMeasurement = -1;
  auto status = detail::runPass<Rule>(prop, geometry, span, measurements, n, PropagationDirection::kForward, seed, fit.outwardState, fit.outwardChi2, fit.failedMeasurement, kinkBudget, outwardKinks, outwardKinkMeasurement, observe);
  fit.outwardKinks = static_cast<std::uint8_t>(outwardKinks);
  if (status != PropagationStatus::kOk) {
    fit.status = status;
    fit.failedPass = FitPass::kOutward;
    return fit;
  }
  detail::inwardPass<Rule>(prop, geometry, span, measurements, n, kinkBudget, observe, fit);
  return fit;
}

/// `fitMeasurements` given its outward pass's end state, \a outwardState
/// -- as the road fit keeps it for an accepted
/// road (`RoadFitResult::state`): only the inward pass runs. The outward
/// pass's chi2 and kinks are not known here and stay zero.
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Solve = ExactScalar, typename Prop, DetectorLayout Det, typename Material, typename Chain, typename Observer = NoFitObserver>
GPUdi() TrackFitResult fitInward(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, const BoundTrackState& outwardState, const KinkBudget& kinkBudget = {}, const Observer& observe = {})
{
  TrackFitResult fit{.ndf = static_cast<std::int16_t>((2 * n) - kNBoundParams)};
  fit.outwardState = outwardState;
  detail::inwardPass<Rule, Solve>(prop, geometry, detail::spanMaterial(geometry, material, measurements, n), measurements, n, kinkBudget, observe, fit);
  return fit;
}

/// The road fit's candidate pass: the outward filter
/// alone, as in `fitMeasurements`.
struct RoadFitResult {
  Scalar chi2 = 0.F;
  RoadFitStatus status = RoadFitStatus::kAccepted;
  std::uint8_t kinks = 0;
  /// The outward pass's state at the outermost measurement: the pT floor
  /// reads it, and `fitMeasurements`' outward pass ends on it, so the refit
  /// of the road takes it instead of running that pass again.
  BoundTrackState state;
};

/// Its chi2 is the whole road's; the inward pass only gives the state on the
/// innermost surface, which the refit of an accepted track provides. The cut
/// is \a chi2MaxByNdf read at the road's degrees of freedom, `2 N - 5` less
/// two per kink. A filter's chi2 only grows as measurements are added, and the
/// cut only falls as kinks are taken, so the pass stops as soon as the chi2
/// exceeds the cut so far: the decision the complete pass would reach.
///
/// With \a LastStepAllowance, for an iteration with a kink allowance, the
/// last step is `allowanceStep` with the allowance's theta0 factor
/// \a lastStepScale (against the road-end fakes): refused
/// beyond the gate, it frees no degree of freedom, and the result's kinks
/// count it. A compile-time choice: without it the last step is a `kinkStep`
/// like the others.
template <bool LastStepAllowance = false, CrossingRule Rule = CrossingRule::kFirstAhead, typename Solve = ExactScalar, typename Prop, DetectorLayout Det, typename Material, typename Chain>
GPUdi() RoadFitResult fitOutward(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, const Scalar* chi2MaxByNdf, const KinkBudget& kinkBudget = {}, [[maybe_unused]] Scalar lastStepScale = 1.F)
{
  const auto budget = detail::kinkBudgetFor(kinkBudget, n);
  BoundTrackState state = makeChainSeed(prop, geometry, measurements, n);
  const auto span = detail::spanMaterial(geometry, material, measurements, n);
  CovarianceFactors factors;
  RoadFitResult result;
  countWork(WorkCounter::kRoadStepsPossible, static_cast<std::uint64_t>(n - 1));
  // The start as `runPass` takes it, so that the kink gate, which reads the
  // increments, rounds as in the refit's outward pass.
  detail::PassReference reference;
  int taken = 0;
  std::int16_t failed = -1;
  if (detail::startPass<Rule, Solve>(prop, geometry, span, measurements, n, PropagationDirection::kForward, state, factors, result.chi2, reference, taken, failed, NoFitObserver{}) != PropagationStatus::kOk) {
    // The steps tried, the failing one included. Not a measurement of this
    // road (a negative or NaN chi2 would pass the cut below and rank
    // first): counted with the propagation failures.
    countWork(WorkCounter::kRoadSteps, static_cast<std::uint64_t>(failed > 0 ? failed : 0));
    countWork(WorkCounter::kRoadPropagation);
    result.status = RoadFitStatus::kPropagationFailed;
    return result;
  }
  countWork(WorkCounter::kRoadSteps, static_cast<std::uint64_t>(taken - 1));
  const int ndf = (2 * n) - kNBoundParams;
  int kinks = 0;
  // The allowance's last step, which frees no degree of freedom; counted by its
  // instantiation alone.
  [[maybe_unused]] int allowed = 0; // NOLINT(misc-const-correctness)
  // A road is at least a triplet, so `ndf` >= 1, and `kinkBudgetFor` keeps
  // `ndf - 2 kinks` there; the analyzer, not knowing n >= 3, takes n = 2.
  // NOLINTBEGIN(clang-analyzer-security.ArrayBound)
  for (int index = taken; index < n; ++index) {
    if (!(result.chi2 <= chi2MaxByNdf[ndf - (2 * kinks)])) {
      countWork(WorkCounter::kRoadChi2);
      result.status = RoadFitStatus::kRejected;
      return result;
    }
    countWork(WorkCounter::kRoadSteps);
    const auto& stepProp = stepPropagator(prop, geometry, measurements, index - 1, index);
    const Scalar stepReference = detail::referenceFor(reference, n, PropagationDirection::kForward, [&] { return detail::leaveOneOutQOverP(prop, geometry, measurements, n, index); });
    Scalar carried = stepReference;
    if constexpr (LastStepAllowance) {
      if (index == n - 1) {
        int kinksAfter = kinks;
        bool refused = false;
        const auto status = detail::allowanceStep<Rule, Solve>(stepProp, geometry, span, measurements[index], PropagationDirection::kForward, state, factors, result.chi2, carried, budget, lastStepScale, kinksAfter, refused);
        if (status != PropagationStatus::kOk) {
          countWork(WorkCounter::kRoadPropagation);
          result.status = RoadFitStatus::kPropagationFailed;
          return result;
        }
        if (refused) {
          countWork(WorkCounter::kRoadChi2);
          result.status = RoadFitStatus::kRejected;
          return result;
        }
        allowed = kinksAfter - kinks;
        continue;
      }
    }
    const auto status = detail::kinkStep<Rule, Solve>(stepProp, geometry, span, measurements[index], PropagationDirection::kForward, state, factors, result.chi2, carried, budget, index < n - 1, kinks);
    if (status != PropagationStatus::kOk) {
      countWork(WorkCounter::kRoadPropagation);
      result.status = RoadFitStatus::kPropagationFailed;
      return result;
    }
    detail::advanceReference(reference, n, PropagationDirection::kForward, stepReference, carried);
  }
  result.kinks = static_cast<std::uint8_t>(kinks + allowed);
  covarianceFromFactors(factors, state.cov);
  result.state = state;
  if (!(result.chi2 <= chi2MaxByNdf[ndf - (2 * kinks)])) {
    // NOLINTEND(clang-analyzer-security.ArrayBound)
    countWork(WorkCounter::kRoadChi2);
    result.status = RoadFitStatus::kRejected;
  }
  return result;
}

/// The outward road fit as an iteration of \a parameters takes it: its kink
/// budget, and with \a LastStepAllowance its kink allowance on the last step
/// (`fitOutward`'s), with the allowance's theta0 factor. The stages make the
/// choice at compile time (`RoadFitStage`).
template <bool LastStepAllowance, CrossingRule Rule = CrossingRule::kFirstAhead, typename Parameters, typename Prop, DetectorLayout Det, typename Material, typename Chain>
[[nodiscard]] GPUdi() RoadFitResult fitRoad(const Parameters& parameters, const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n)
{
  if constexpr (LastStepAllowance) {
    return fitOutward<true, Rule, Scalar>(prop, geometry, material, measurements, n, parameters.roadChi2MaxByNdf, kinkBudget(parameters), o2::gpu::CAMath::Sqrt(parameters.kinkScale2));
  } else {
    return fitOutward<false, Rule, Scalar>(prop, geometry, material, measurements, n, parameters.roadChi2MaxByNdf, kinkBudget(parameters));
  }
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_KFTRACKFITTER_H_
