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

/// \file HypothesisRefit.h
/// \brief A found track's rows refitted under a species hypothesis, its PID's, and carried to the production vertex.
///
/// The tracking fits every track as a pion. A consumer that knows a track's
/// species -- from its PID after the matching to the outer detectors -- refits
/// its rows under that hypothesis: an electron with the Gaussian-sum filter
/// (its bremsstrahlung), any other species with the Kalman fit at its mass.
/// The PID chooses; the consumer never picks a fitter. The result is the state
/// on the innermost measurement and the state on the perigee plane about a
/// transverse point (the beam spot or a vertex): the plane through that point
/// along the beam, normal to the track's direction at its closest approach to
/// the point (`beamLinePlane`, `perigeePhi`); loc0 the transverse
/// distance along the plane, the track's d0 where the track crosses the plane
/// normally, loc1 the z. An electron's mixture components cross this one
/// plane each at its own angle: its loc0 is the d0 to first order in the
/// components' spread.

#ifndef GTRK_TRACKING_HYPOTHESISREFIT_H_
#define GTRK_TRACKING_HYPOTHESISREFIT_H_

#include "GenericTracking/Tracking/BeamConstraint.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/GSFTrackFitter.h"
#include "GenericTracking/Tracking/KFTrackFitter.h"
#include "GenericTracking/Tracking/MeasurementUpdate.h"

#include "ReconstructionDataFormats/PID.h"

#include <cstdint>

namespace o2::gtrk
{

/// A refit under a species hypothesis (`refitWithHypothesis`). For an
/// electron the single states are the mixtures' collapses -- their mean and
/// covariance, a Gaussian a vertexer can take -- and the momentum estimate
/// apart: the mode of the momentum's density (`GSFFitResult::momentumMode`),
/// the estimate the mixture's posterior is skewed about.
template <int K = BetheHeitlerLoss::kComponents>
struct HypothesisFit {
  o2::track::PID pid;
  BoundTrackState innermost{}; ///< on the innermost measurement, the inward pass's end
  BoundTrackState perigee{};   ///< on `perigeePlane`
  TrackingSurface perigeePlane{};
  /// The momentum's estimate [GeV] on the innermost measurement and on the
  /// perigee plane: 1 / |q/p| of the Kalman fit's state; an electron's the
  /// mode of its mixture's momentum density.
  Scalar innermostMomentum = 0.F;
  Scalar perigeeMomentum = 0.F;
  /// An electron's mixtures on the innermost measurement and on the perigee
  /// plane; empty for another species.
  GSFMixture<K> innermostMixture{};
  GSFMixture<K> perigeeMixture{};
  /// The inward pass's chi2 and its degrees of freedom (TrackFitResult; an
  /// electron's the mixture's, `GSFFitResult`).
  Scalar chi2 = 0.F;
  std::int16_t ndf = 0;
  PropagationStatus status = PropagationStatus::kOk;
  /// The fit's pass that failed; kNone with a failed status: the
  /// extrapolation to the perigee plane failed.
  FitPass failedPass = FitPass::kNone;
  std::int16_t failedMeasurement = -1;
};

namespace detail::hypothesis
{

/// Whether \a to lies ahead of \a state on \a from, along its direction:
/// a vertex fit's start or iterate beyond the daughter's innermost
/// measurement, where it never was (the fit's solution can lie there).
GPUdi() bool planeAhead(const BoundTrackState& state, const TrackingSurface& from, const TrackingSurface& to)
{
  Scalar at[3]{};
  surfaceToGlobal(from, state.loc0(), state.loc1(), at);
  const Scalar sinTheta = sinOf(state.theta());
  const Scalar direction[3] = {sinTheta * cosOf(state.phi()), sinTheta * sinOf(state.phi()), cosOf(state.theta())};
  const Scalar offset[3] = {to.origin[0] - at[0], to.origin[1] - at[1], to.origin[2] - at[2]};
  return dotOf(to.normal, offset) * dotOf(to.normal, direction) > 0.F;
}

/// \a state, its covariance \a factors, carried from \a from to \a to against
/// the momentum, inward, as a Kalman step carries a state between two
/// measurements (`predictStep`, `transportStep`) with \a material's on the
/// straight segment between (its `segment`, LinkMaterial::segment: the band in front of
/// \a from's measurement and the share of its vertex gap down to the plane's
/// crossing), the ionisation over the segment undone. The perigee plane is no
/// sensor: its bounds are empty, and the propagation's out-of-bounds status
/// is the crossing (`predictBeamLine`); kNoMaterial if the table holds no
/// such segment. A plane ahead (`planeAhead`): the state carried forward
/// along its helix, no material on the way (none was crossed).
template <typename Material>
GPUdi() PropagationStatus carryInward(const Propagator& prop, const Material& material, const TrackingSurface& from, const TrackingSurface& to, BoundTrackState& state, CovarianceFactors& factors)
{
  if (planeAhead(state, from, to)) {
    SolenoidJacobian jacobian;
    const auto status = jacobianAnalytic(prop, state, from, to, PropagationDirection::kForward, jacobian, &state);
    if (status != PropagationStatus::kOk && status != PropagationStatus::kOutOfBounds) {
      return status;
    }
    state.surface = to.id;
    return transportFactors(jacobian, StepNoise{}, factors) ? PropagationStatus::kOk : PropagationStatus::kNotFinite;
  }
  Scalar start[3]{};
  surfaceToGlobal(from, state.loc0(), state.loc1(), start);
  const Scalar leaving = state.theta();
  SolenoidJacobian jacobian;
  const auto status = jacobianAnalytic(prop, state, from, to, PropagationDirection::kBackward, jacobian, &state);
  if (status != PropagationStatus::kOk && status != PropagationStatus::kOutOfBounds) {
    return status;
  }
  state.surface = to.id;
  Scalar end[3]{};
  surfaceToGlobal(to, state.loc0(), state.loc1(), end);
  const auto step = o2::gtrk::detail::stepMaterialOf(material.segment(from.layer, start, kInvalidLayerID, end));
  if (!step.found) {
    return PropagationStatus::kNoMaterial;
  }
  const Scalar momentum = 1.F / o2::gpu::GPUCommonMath::Abs(state.qOverP());
  const auto before = scatteringNoise(highlandTheta0(step.nearerX2X0, step.nearerLogX2X0, momentum, prop.massGeV), leaving);
  return o2::gtrk::detail::transportStep<false>(prop, step, momentum, before, PropagationDirection::kBackward, jacobian, state, factors, 1.F);
}

/// The perigee plane about (\a x, \a y) of \a state on \a surface: through
/// the point along the beam, normal to the direction at closest approach to it.
GPUdi() TrackingSurface perigeePlaneOf(const Propagator& prop, const BoundTrackState& state, const TrackingSurface& surface, Scalar x, Scalar y)
{
  return beamLinePlane(x, y, perigeePhi(prop, state, surface, x, y));
}

/// The radiation [x/X0] of the straight segment from \a state on \a surface
/// to where it crosses \a plane: the material inward of the innermost
/// measurement, which no plane of the chain holds (`planeRadiation`), as
/// `carryInward` takes it; -1 where it cannot be.
template <typename Material>
GPUdi() Scalar segmentRadiation(const Propagator& prop, const Material& material, const BoundTrackState& state, const TrackingSurface& surface, const TrackingSurface& plane)
{
  if (planeAhead(state, surface, plane)) {
    return 0.F; // carried forward, through nothing
  }
  BoundTrackState probe = state;
  SolenoidJacobian unused;
  const auto status = jacobianAnalytic(prop, state, surface, plane, PropagationDirection::kBackward, unused, &probe);
  if (status != PropagationStatus::kOk && status != PropagationStatus::kOutOfBounds) {
    return -1.F;
  }
  Scalar start[3]{};
  Scalar end[3]{};
  surfaceToGlobal(surface, state.loc0(), state.loc1(), start);
  surfaceToGlobal(plane, probe.loc0(), probe.loc1(), end);
  const auto segment = material.segment(surface.layer, start, kInvalidLayerID, end);
  return segment.found ? segment.before + segment.after : -1.F;
}

/// \a mixture's components carried inward from \a surface onto \a plane
/// (`carryInward`), the segment's radiation \a x2X0 (`segmentRadiation`)
/// taken first, branched and reduced as on a plane; kNotFinite if it cannot
/// be.
template <typename Mixture, typename Material>
GPUdi() PropagationStatus carryMixtureInward(const Propagator& prop, const Material& material, const TrackingSurface& surface, const TrackingSurface& plane, Scalar x2X0, Mixture& mixture)
{
  if (!(x2X0 >= 0.F) || !detail::gsf::takeRadiation<BetheHeitlerLoss>(x2X0, PropagationDirection::kBackward, mixture)) {
    return PropagationStatus::kNotFinite;
  }
  return detail::gsf::stepEach(mixture, [&](GSFComponent<CovarianceFactors>& component) {
    return carryInward(prop, material, surface, plane, component.state, component.factors);
  });
}

/// \a state with its covariance carried inward from \a surface onto
/// \a plane (`carryInward`).
template <typename Material>
GPUdi() PropagationStatus carryStateInward(const Propagator& prop, const Material& material, const TrackingSurface& surface, const TrackingSurface& plane, BoundTrackState& state)
{
  auto factors = factorCovariance(state.cov);
  const auto status = carryInward(prop, material, surface, plane, state, factors);
  if (status == PropagationStatus::kOk) {
    covarianceFromFactors(factors, state.cov);
  }
  return status;
}

} // namespace detail::hypothesis

/// \a n measurements of \a measurements, ordered from the innermost, refitted
/// under \a pid's hypothesis, up to the state on the innermost measurement.
/// The propagator's mass becomes the species', for the scattering and the
/// ionisation alike. An electron: the Gaussian-sum fit (`fitGaussianSum`),
/// its inward mixture kept. Another species: the Kalman fit
/// (`fitMeasurements`) at its mass. A species of charge other than one is
/// refused (kNotFinite, failedPass kNone): the energy loss takes a unit
/// charge. The perigee fields stay empty (`carryToPerigee`).
template <int K = BetheHeitlerLoss::kComponents, CrossingRule Rule = CrossingRule::kFirstAhead, DetectorLayout Det, typename Material, typename Chain>
GPUdi() HypothesisFit<K> refitRowsWithHypothesis(const Propagator& base, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, o2::track::PID pid)
{
  HypothesisFit<K> fit{.pid = pid};
  if (o2::gpu::GPUCommonMath::Abs(pid.getCharge()) != 1) {
    fit.status = PropagationStatus::kNotFinite;
    return fit;
  }
  Propagator prop = base;
  prop.massGeV = pid.getMass();
  // What both fits give alike (TrackFitResult); false if the fit failed.
  const auto take = [&](const TrackFitResult& result, Scalar momentum) {
    fit.status = result.status;
    fit.failedPass = result.failedPass;
    fit.failedMeasurement = result.failedMeasurement;
    if (result.status != PropagationStatus::kOk) {
      return false;
    }
    fit.innermost = result.state;
    fit.innermostMomentum = momentum;
    fit.chi2 = result.chi2;
    fit.ndf = result.ndf;
    return true;
  };
  if (pid == o2::track::PID::Electron) {
    const auto gsf = fitGaussianSum<BetheHeitlerLoss, K, Rule>(prop, geometry, material, measurements, n);
    if (take(gsf.fit, gsf.momentumMode)) {
      fit.innermostMixture = gsf.mixture;
    }
    return fit;
  }
  const auto kalman = fitMeasurements<Rule>(prop, geometry, material, measurements, n);
  take(kalman, 1.F / o2::gpu::GPUCommonMath::Abs(kalman.state.qOverP()));
  return fit;
}

/// \a fit, a rows' refit that succeeded (`refitRowsWithHypothesis`), its
/// innermost measurement on \a innermostSurface, carried to the perigee plane
/// about (\a x, \a y) with the material between: an electron's mixture, the
/// segment's radiation undone first (branched and reduced, as on a plane),
/// each component on its own trajectory onto the collapse's perigee plane;
/// another species' state. Its perigee fields and status set.
template <int K, DetectorLayout Det>
GPUdi() PropagationStatus carryToPerigee(const Propagator& base, const TrackingGeometry<Det>& geometry, const TrackingSurface& innermostSurface, HypothesisFit<K>& fit, Scalar x, Scalar y)
{
  Propagator prop = base;
  prop.massGeV = fit.pid.getMass();
  fit.perigeePlane = detail::hypothesis::perigeePlaneOf(prop, fit.innermost, innermostSurface, x, y);
  if (fit.pid == o2::track::PID::Electron) {
    GSFMixture<K> mixture = fit.innermostMixture;
    const Scalar x2X0 = detail::hypothesis::segmentRadiation(prop, geometry.material(), fit.innermost, innermostSurface, fit.perigeePlane);
    fit.status = detail::hypothesis::carryMixtureInward(prop, geometry.material(), innermostSurface, fit.perigeePlane, x2X0, mixture);
    if (fit.status == PropagationStatus::kOk) {
      detail::gsf::finishMixture(mixture, fit.perigee, fit.perigeeMomentum);
      fit.perigeeMixture = mixture;
    }
    return fit.status;
  }
  BoundTrackState state = fit.innermost;
  fit.status = detail::hypothesis::carryStateInward(prop, geometry.material(), innermostSurface, fit.perigeePlane, state);
  if (fit.status == PropagationStatus::kOk) {
    fit.perigee = state;
    fit.perigeeMomentum = 1.F / o2::gpu::GPUCommonMath::Abs(state.qOverP());
  }
  return fit.status;
}

/// The rows' refit (`refitRowsWithHypothesis`) carried to the perigee plane
/// about (\a x, \a y) (`carryToPerigee`); a failed carry leaves failedPass
/// kNone with the failed status.
template <int K = BetheHeitlerLoss::kComponents, CrossingRule Rule = CrossingRule::kFirstAhead, DetectorLayout Det, typename Material, typename Chain>
GPUdi() HypothesisFit<K> refitWithHypothesis(const Propagator& base, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, o2::track::PID pid, Scalar x, Scalar y)
{
  auto fit = refitRowsWithHypothesis<K, Rule>(base, geometry, material, measurements, n, pid);
  if (fit.status == PropagationStatus::kOk) {
    carryToPerigee(base, geometry, geometry.planeOf(measurements[0]), fit, x, y);
  }
  return fit;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_HYPOTHESISREFIT_H_
