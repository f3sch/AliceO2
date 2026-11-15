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
/// the point (`beamLinePlane`, `perigeePhiAbout`); loc0 the transverse
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
/// apart: the mode of the momentum's density (`GSFPassResult::momentumMode`,
/// not the inverse of the q/p marginal's mode: the density of p is the q/p
/// one's times (q/p)^2), the estimate the mixture's posterior is skewed about.
template <int K = BetheHeitlerLoss::kComponents>
struct HypothesisFit {
  o2::track::PID pid{};
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
  /// The inward pass's: the Kalman fit's chi2 and its degrees of freedom; for
  /// an electron the pass's -2 ln likelihood (`GSFPassResult`) and the same
  /// degrees of freedom.
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

/// \a state, its covariance \a factors, carried from \a from to \a to against
/// the momentum, inward, as a Kalman step carries a state between two
/// measurements (`predictStep`, `transportStep`) with \a lut's material on the
/// straight segment between: the half nearer \a from scatters before the
/// transport, the other half after, the ionisation over the segment undone.
/// Each half's kick sits at its end of the segment: a scatterer's lever arm
/// to the plane is its half's, as a fit's step takes it. The perigee plane is
/// no sensor: its bounds are empty, and the
/// propagation's out-of-bounds status is the crossing (`predictBeamLine`).
template <typename Material>
GPUdi() PropagationStatus carryInward(const Propagator& prop, const Material& lut, const TrackingSurface& from, const TrackingSurface& to, BoundTrackState& state, CovarianceFactors& factors)
{
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
  const auto step = o2::gtrk::detail::stepMaterial(lut, start, end);
  const Scalar momentum = 1.F / o2::gpu::GPUCommonMath::Abs(state.qOverP());
  const auto before = scatteringNoise(highlandTheta0(step.nearerX2X0, step.nearerLogX2X0, momentum, prop.massGeV), leaving);
  return o2::gtrk::detail::transportStep<false>(prop, step, momentum, before, PropagationDirection::kBackward, jacobian, state, factors, 1.F);
}

/// The direction [rad] of \a state's helix at its transverse closest
/// approach to the line along the beam through (\a x, \a y): `perigeePhi`'s
/// solution, the position taken from that point.
GPUdi() Scalar perigeePhiAbout(const Propagator& prop, const BoundTrackState& state, const TrackingSurface& surface, Scalar x, Scalar y)
{
  ExactScalar position[3]{};
  o2::gtrk::detail::boundToGlobal(surface, static_cast<ExactScalar>(state.loc0()), static_cast<ExactScalar>(state.loc1()), position);
  position[0] -= static_cast<ExactScalar>(x);
  position[1] -= static_cast<ExactScalar>(y);
  ExactScalar sinPhi = 0.;
  ExactScalar cosPhi = 0.;
  ExactScalar sinTheta = 0.;
  ExactScalar cosTheta = 0.;
  sincosOf(static_cast<ExactScalar>(state.phi()), sinPhi, cosPhi);
  sincosOf(static_cast<ExactScalar>(state.theta()), sinTheta, cosTheta);
  const ExactScalar curvature = pointPropagator(prop, surface, state).curvature(static_cast<ExactScalar>(state.qOverP()), sinTheta);
  const ExactScalar along = (position[0] * cosPhi) + (position[1] * sinPhi);
  const ExactScalar across = (-position[0] * sinPhi) + (position[1] * cosPhi);
  ExactScalar phi0 = static_cast<ExactScalar>(state.phi()) + atan2Of(-curvature * along, 1. + (curvature * across));
  math_utils::detail::bringToPMPi<ExactScalar>(phi0);
  return static_cast<Scalar>(phi0);
}

/// The perigee plane about (\a x, \a y) of \a state on \a surface: through
/// the point along the beam, normal to the direction at closest approach to it.
GPUdi() TrackingSurface perigeePlaneOf(const Propagator& prop, const BoundTrackState& state, const TrackingSurface& surface, Scalar x, Scalar y)
{
  return beamLinePlane(x, y, perigeePhiAbout(prop, state, surface, x, y));
}

/// The radiation [x/X0] of the straight segment from \a state on \a surface
/// to where it crosses \a plane: the material inward of the innermost
/// measurement, which no plane of the chain holds (`planeRadiation`).
template <typename Material>
GPUdi() Scalar segmentRadiation(const Propagator& prop, const Material& lut, const BoundTrackState& state, const TrackingSurface& surface, const TrackingSurface& plane)
{
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
  return lut.getMatBudget(start[0], start[1], start[2], end[0], end[1], end[2]).meanX2X0;
}

} // namespace detail::hypothesis

/// \a n measurements of \a measurements, ordered from the innermost, refitted
/// under \a pid's hypothesis and carried to the perigee plane about
/// (\a x, \a y). The propagator's mass becomes the species', for the
/// scattering and the ionisation alike. An electron: the Gaussian-sum fit
/// (`fitGaussianSum`), its inward mixture carried to the perigee plane, the
/// segment's radiation undone first (branched and reduced, as on a plane),
/// each component on its own trajectory onto the collapse's perigee plane.
/// Another species: the Kalman fit (`fitMeasurements`) at its mass, its
/// inward state carried so. A species of charge other than one is refused
/// (kNotFinite, failedPass kNone): the energy loss takes a unit charge.
template <int K = BetheHeitlerLoss::kComponents, CrossingRule Rule = CrossingRule::kFirstAhead, DetectorLayout Det, typename Material, typename Chain>
GPUdi() HypothesisFit<K> refitWithHypothesis(const Propagator& base, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, o2::track::PID pid, Scalar x, Scalar y)
{
  HypothesisFit<K> fit{.pid = pid};
  if (o2::gpu::GPUCommonMath::Abs(pid.getCharge()) != 1) {
    fit.status = PropagationStatus::kNotFinite;
    return fit;
  }
  Propagator prop = base;
  prop.massGeV = pid.getMass();
  const int innermostIndex = 0;
  const auto& innermostSurface = geometry.planeOf(measurements[innermostIndex]);
  if (pid == o2::track::PID::Electron) {
    const auto gsf = fitGaussianSum<BetheHeitlerLoss, K, Rule>(prop, geometry, material, measurements, n);
    const auto& inward = gsf.inward;
    if (inward.status != PropagationStatus::kOk) {
      fit.status = inward.status;
      fit.failedPass = gsf.outward.status != PropagationStatus::kOk ? FitPass::kOutward : FitPass::kInward;
      fit.failedMeasurement = inward.failedMeasurement;
      return fit;
    }
    fit.innermostMixture = inward.mixture;
    fit.innermost = inward.collapsed;
    fit.innermostMomentum = inward.momentumMode;
    fit.chi2 = inward.minusTwoLogLikelihood;
    fit.ndf = static_cast<std::int16_t>((2 * n) - kNBoundParams);
    fit.perigeePlane = detail::hypothesis::perigeePlaneOf(prop, inward.collapsed, innermostSurface, x, y);
    GSFMixture<K> mixture = inward.mixture;
    const Scalar x2X0 = detail::hypothesis::segmentRadiation(prop, material, inward.collapsed, innermostSurface, fit.perigeePlane);
    if (!(x2X0 >= 0.F) || !detail::gsf::takeRadiation<BetheHeitlerLoss>(x2X0, PropagationDirection::kBackward, mixture)) {
      fit.status = PropagationStatus::kNotFinite;
      return fit;
    }
    int kept = 0;
    PropagationStatus failure = PropagationStatus::kOk;
    for (int k = 0; k < mixture.size; ++k) {
      const auto status = detail::hypothesis::carryInward(prop, material, innermostSurface, fit.perigeePlane, mixture.state[k], mixture.factors[k]);
      if (status != PropagationStatus::kOk) {
        failure = status;
        continue;
      }
      mixture.logWeight[kept] = mixture.logWeight[k];
      mixture.state[kept] = mixture.state[k];
      mixture.factors[kept] = mixture.factors[k];
      mixture.chi2[kept] = mixture.chi2[k];
      ++kept;
    }
    mixture.size = kept;
    if (kept == 0) {
      fit.status = failure;
      return fit;
    }
    detail::gsf::normalise(mixture);
    for (int k = 0; k < mixture.size; ++k) {
      covarianceFromFactors(mixture.factors[k], mixture.state[k].cov);
    }
    fit.perigeeMixture = mixture;
    fit.perigee = detail::gsf::collapse(mixture);
    fit.perigeeMomentum = 1.F / o2::gpu::GPUCommonMath::Abs(detail::gsf::marginalMode<2>(mixture));
    return fit;
  }
  const auto kalman = fitMeasurements<Rule>(prop, geometry, material, measurements, n);
  if (kalman.status != PropagationStatus::kOk) {
    fit.status = kalman.status;
    fit.failedPass = kalman.failedPass;
    fit.failedMeasurement = kalman.failedMeasurement;
    return fit;
  }
  fit.innermost = kalman.state;
  fit.innermostMomentum = 1.F / o2::gpu::GPUCommonMath::Abs(kalman.state.qOverP());
  fit.chi2 = kalman.chi2;
  fit.ndf = kalman.ndf;
  fit.perigeePlane = detail::hypothesis::perigeePlaneOf(prop, kalman.state, innermostSurface, x, y);
  BoundTrackState state = kalman.state;
  auto factors = factorCovariance(state.cov);
  fit.status = detail::hypothesis::carryInward(prop, material, innermostSurface, fit.perigeePlane, state, factors);
  if (fit.status != PropagationStatus::kOk) {
    return fit;
  }
  covarianceFromFactors(factors, state.cov);
  fit.perigee = state;
  fit.perigeeMomentum = 1.F / o2::gpu::GPUCommonMath::Abs(state.qOverP());
  return fit;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_HYPOTHESISREFIT_H_
