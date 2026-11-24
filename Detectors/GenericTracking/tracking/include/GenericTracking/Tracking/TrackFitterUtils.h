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

/// \file TrackFitterUtils.h
/// \brief What the track fitters' steps share: a step's material, its prediction and update, and the exact diffuse filter that starts a pass.
///
/// The Kalman fit (KFTrackFitter.h) and the Gaussian-sum fit
/// (GSFTrackFitter.h) take the same step between two measurements and start
/// a pass alike; their passes are their own.

#ifndef GTRK_TRACKING_TRACKFITTERUTILS_H_
#define GTRK_TRACKING_TRACKFITTERUTILS_H_

#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/CovarianceTransport.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/MeasurementUpdate.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/Seed.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include "GPUCommonMath.h"
#include "MathUtils/detail/trigonometric.h"

#include <cfloat>
#include <cstdint>
#include <type_traits>

namespace o2::gtrk
{

/// Result of the two independent filter passes.
///
/// `state` is the inward result on the innermost measurement surface. The
/// outward state is the starting point for extension; its chi2 remains
/// available for QA.
struct TrackFitResult {
  BoundTrackState state{};
  BoundTrackState outwardState{};
  Scalar chi2 = 0.F;
  Scalar outwardChi2 = 0.F;
  /// Two coordinates per measurement against the five parameters, two fewer
  /// per kink.
  std::int16_t ndf = 0;
  /// kOk, or what the failed step's propagation gave (`failedPass`,
  /// `failedMeasurement` say where).
  PropagationStatus status = PropagationStatus::kOk;
  FitPass failedPass = FitPass::kNone;
  /// The measurement the failed step was propagating to.
  std::int16_t failedMeasurement = -1;
  /// The kinks each pass took. `ndf` is the inward pass's, two fewer per
  /// kink; the outward pass has `2 N - 5 - 2 outwardKinks`.
  std::uint8_t kinks = 0;
  std::uint8_t outwardKinks = 0;
  /// The measurement the inward pass's last kink updated with; -1 without.
  std::int16_t kinkMeasurement = -1;
};

/// The mass a fit step's energy loss takes: the propagator's hypothesis,
/// the one its scattering takes. A study's propagator type overloads it to
/// set the two apart (the study's TruthFit.h).
GPUhdi() GPUAlwaysInline() Scalar lossMassGeV(const Propagator& prop) { return prop.massGeV; }

/// The bremsstrahlung a fit step adds to its energy loss over \a x2X0
/// radiation lengths: none for the tracking's propagator. A study's
/// propagator type overloads it (the study's TruthFit.h) to change \a state's q/p
/// and add to the step's q/p variance \a varQOverP.
GPUhdi() GPUAlwaysInline() void applyRadiation(const Propagator& /*prop*/, Scalar /*x2X0*/, PropagationDirection /*direction*/, BoundParameters& /*state*/, Scalar& /*varQOverP*/) {}

namespace detail
{

/// A step's material in two halves about its midpoint, and its `rho * length`
/// for the energy loss.
///
/// The Highland noise of the half nearer the surface the step leaves goes into
/// the covariance *before* the transport, so that its kink becomes position
/// spread at the next surface through `F Q F^T`; the other half's goes in
/// after it, before the update. Added after the transport, all of it, every
/// prediction would miss `theta0^2 L^2` of the material near the start and the
/// outward state would carry one step's unconsumed `theta0^2` in its direction.
/// With a layer's material
/// attributed from half-gap to half-gap, the two halves are the stave just
/// left and the stave about to be reached.
///
/// Each half's Highland logarithm is taken at the scatterer it is a part of
/// (`highlandTheta0`): its own material for a step alone, as a candidate's
/// step in the follow; the fitted span's for a fit's step (`SpanMaterial`).
///
/// From the link material (LinkMaterial.h, `LinkStep`) the nearer part is the
/// band left and the gap, a kick at the gap's mean place -- \a lever, the
/// fraction of the step after it -- with the spread of its places about it,
/// \a spread; the farther part the band reached.
struct StepMaterial {
  Scalar nearerX2X0 = 0.F;     ///< x/X0 from the start to the midpoint (the link material's band left and gap)
  Scalar fartherX2X0 = 0.F;    ///< x/X0 from the midpoint to the end (the link material's band reached)
  Scalar nearerLogX2X0 = 0.F;  ///< x/X0 the nearer half's logarithm is taken at
  Scalar fartherLogX2X0 = 0.F; ///< x/X0 the farther half's logarithm is taken at
  Scalar xRho = 0.F;           ///< [g/cm^2] over the whole step
  Scalar lever = 1.F;          ///< 1 - <f> of the nearer part, f the fraction of the step
  Scalar spread = 0.F;         ///< <f^2> - <f>^2 of the nearer part
  bool found = true;           ///< false: the link material holds no such step (LinkMaterial::step)
};

/// What a step's material is: a policy the fits and the follow take by
/// type, every step's material the link material's (LinkMaterial.h, the
/// geometry's). The plain one: each step's own, Highland's logarithm at its
/// own x/X0.
struct LinkSteps {
};

/// A fit's material: \a Material's steps, with every step's Highland
/// logarithm at \a spanX2X0, the material of the fitted span (`spanX2X0`).
/// The scattering a fit takes over its span is the span's deflection, which
/// Highland's formula gives for the span as a whole; each step takes its
/// share (`highlandTheta0`). Both passes of a fit take the same span.
template <typename Material>
struct SpanMaterial {
  Scalar spanX2X0 = 0.F;
};

/// The step of a link-material step \a link (LinkMaterial::step,
/// `segment`): its `before` the nearer part, `after` the farther,
/// Highland's logarithm at the step's whole x/X0 for both.
GPUhdi() StepMaterial stepMaterialOf(const LinkStep& link)
{
  const Scalar span = link.before + link.after;
  return {.nearerX2X0 = link.before, .fartherX2X0 = link.after, .nearerLogX2X0 = span, .fartherLogX2X0 = span, .xRho = link.xRho, .lever = link.lever, .spread = link.spread, .found = link.found};
}

/// A step's material known before the step, not looked up: a follow step's
/// window's (`follow::predictedWindow`), for every candidate the window
/// holds. They leave the window's state for the same layer and end within
/// the window about its nominal crossing, so they cross what it crossed.
struct GivenMaterial {
  StepMaterial step;
};

GPUdi() StepMaterial stepMaterial(const GivenMaterial& material, const Scalar* /*start*/, const Scalar* /*end*/)
{
  return material.step;
}

/// A fit's material where its chain may hold two measurements on one layer
/// (overlap partners, OverlapPartners.h): the link material's between
/// layers, whose bands hold each layer's whole material, and none on a step
/// between two sensors of one layer, which the bands hold already.
struct LayerMaterial {
};

template <typename Material>
struct IsLayerMaterial : std::false_type {
};
template <>
struct IsLayerMaterial<LayerMaterial> : std::true_type {
};

/// A step's material from the link material (LinkMaterial.h) between the
/// planes \a from and \a to, `found` false where it holds none:
/// \a material's given step for a follow candidate (`GivenMaterial`), each
/// half's logarithm at the fit's span for a fit (`SpanMaterial`), at the
/// step's own otherwise.
template <typename Material>
struct IsSpanMaterial : std::false_type {
};
template <typename Material>
struct IsSpanMaterial<SpanMaterial<Material>> : std::true_type {
};
template <>
struct IsLayerMaterial<SpanMaterial<LayerMaterial>> : std::true_type {
};

template <DetectorLayout Det, typename Material>
GPUdi() StepMaterial stepMaterial(const TrackingGeometry<Det>& geometry, const TrackingSurface& from, const TrackingSurface& to, const Material& material, const Scalar start[3], const Scalar end[3])
{
  if constexpr (std::is_same_v<Material, GivenMaterial>) {
    return material.step;
  } else {
    if constexpr (IsLayerMaterial<Material>::value) {
      if (from.layer == to.layer) {
        StepMaterial none;
        if constexpr (IsSpanMaterial<Material>::value) {
          none.nearerLogX2X0 = material.spanX2X0;
          none.fartherLogX2X0 = material.spanX2X0;
        }
        return none;
      }
    }
    // From a point off the layers (a decay vertex's plane, no layer): the
    // segment's share of the material table's gaps.
    const auto link = from.layer == kInvalidLayerID ? geometry.material().segment(kInvalidLayerID, start, to.layer, end) : geometry.material().step(from.layer, to.layer, start, end);
    if (!link.found) {
      StepMaterial none;
      none.found = false;
      return none;
    }
    auto step = stepMaterialOf(link);
    if constexpr (IsSpanMaterial<Material>::value) {
      step.nearerLogX2X0 = material.spanX2X0;
      step.fartherLogX2X0 = material.spanX2X0;
    }
    return step;
  }
}

/// The mass hypothesis's energy loss over a step's \a xRho [g/cm^2], O2's
/// way, evaluated at \a lossQOverP (the pass's reference, or \a state's own
/// q/p): the loss's fraction of the momentum, its straggling's q/p variance
/// in \a varQOverP for the step's noise (StepNoise), and \a lossQOverP and
/// \a state both scaled by that fraction. Lost along the momentum, regained
/// against it. False if the particle stopped. The mass is
/// `lossMassGeV(prop)`, the propagator's hypothesis, as is the |q|.
template <typename Prop>
GPUhdi() GPUAlwaysInline() bool applyEnergyLoss(const Prop& prop, Scalar xRho, PropagationDirection direction, Scalar& lossQOverP, BoundParameters& state, Scalar& varQOverP)
{
  varQOverP = 0.F;
  const auto loss = energyLoss(lossQOverP, lossMassGeV(prop), prop.absCharge, direction == PropagationDirection::kForward ? -xRho : xRho);
  if (!loss.ok) {
    return false;
  }
  lossQOverP *= loss.momentumScale;
  state.par[kQOverP] *= loss.momentumScale;
  varQOverP = loss.varQOverP;
  return true;
}

/// \a theta0 times \a scale with \a Scaled, \a theta0 alone without: the
/// nominal step never multiplies.
template <bool Scaled>
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar scaledTheta0(Scalar theta0, [[maybe_unused]] Scalar scale)
{
  if constexpr (Scaled) {
    return scale * theta0;
  } else {
    return theta0;
  }
}

/// A step's noise, once \a state is on the plane the step reaches: the
/// farther half's scattering of \a step at the direction \a state arrives
/// with and at \a momentum, and the energy loss over \a step at
/// \a lossQOverP (`applyEnergyLoss`), applied to \a state and to
/// \a lossQOverP, its straggling's variance in \a noise; \a before is the
/// nearer half's (StepMaterial). With \a Scaled, \a noiseScale multiplies
/// the farther half's theta0 (`scaledTheta0`).
template <bool Scaled, typename Prop>
GPUdi() PropagationStatus stepNoise(const Prop& prop, const StepMaterial& step, Scalar momentum, const ScatteringNoise& before, PropagationDirection direction, BoundParameters& state, Scalar& lossQOverP, StepNoise& noise, Scalar noiseScale)
{
  const auto after = scatteringNoise(scaledTheta0<Scaled>(highlandTheta0(step.fartherX2X0, step.fartherLogX2X0, momentum, prop.massGeV, prop.absCharge), noiseScale), state.theta());
  noise = StepNoise{.phiBefore = before.phi, .thetaBefore = before.theta, .phiAfter = after.phi, .thetaAfter = after.theta, .lever = step.lever, .spread = step.spread};
  if (!applyEnergyLoss(prop, step.xRho, direction, lossQOverP, state, noise.qOverPAfter)) {
    return PropagationStatus::kStopped;
  }
  applyRadiation(prop, step.nearerX2X0 + step.fartherX2X0, direction, state, noise.qOverPAfter);
  return PropagationStatus::kOk;
}

/// A step's second half, once \a state is on the plane the step reaches and
/// \a jacobian carries the covariance there: its noise (`stepNoise`) and the
/// covariance \a factors transported with it. The fit's step (`filterStep`)
/// and the follow's (`transportCandidate`) both end so; each takes the nearer
/// half's noise before its propagation.
template <bool Scaled, typename Prop>
GPUdi() PropagationStatus transportStep(const Prop& prop, const StepMaterial& step, Scalar momentum, const ScatteringNoise& before, PropagationDirection direction, const SolenoidJacobian& jacobian, BoundParameters& state, CovarianceFactors& factors, Scalar noiseScale)
{
  StepNoise noise;
  Scalar lossQOverP = state.qOverP();
  const auto status = stepNoise<Scaled>(prop, step, momentum, before, direction, state, lossQOverP, noise, noiseScale);
  if (status != PropagationStatus::kOk) {
    return status;
  }
  if (!transportFactors(jacobian, noise, factors)) {
    return PropagationStatus::kNotFinite;
  }
  return PropagationStatus::kOk;
}

/// A filter step's prediction: \a state carried to \a measurement's plane,
/// linearised about the filtered state's position and direction and the
/// pass's reference momentum \a referenceQOverP (PassReference). That point is
/// propagated exactly, the Jacobian and the material -- the scattering and the
/// mean energy loss with its straggling -- are taken along it, and the
/// prediction adds the Jacobian's q/p column times the filtered q/p's
/// difference from the reference; the step's energy loss carries the
/// reference to the plane reached.
///
/// The position and direction are the filtered state's: a reference
/// trajectory that the filter never corrects lies off a scattered track by
/// the deflections taken since it started, and the transport's second-order
/// terms about it then bias the fitted angles of a low-momentum track, most
/// where long steps reach planes across the track (discs). The momentum is
/// the reference's: the filtered q/p is only as precise as the measurements
/// taken so far, and after a pass's start it rests on three, which on a weakly
/// bending track (planes along the beam) leave it far from the particle's
/// momentum. Propagated at that curvature, the prediction's second-order
/// error goes into the fitted angles, and scattering evaluated there
/// misweights every measurement that follows. \a jacobian and \a noise are
/// what carries a covariance along (`transportFactors`). \a Scaled, the
/// step's scattering theta0 is multiplied by \a noiseScale, both halves.
template <bool Scaled, CrossingRule Rule, typename Solve = ExactScalar, DetectorLayout Det, typename Material, typename Prop>
GPUdi() PropagationStatus predictStep(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& state, Scalar& referenceQOverP, SolenoidJacobian& jacobian, StepNoise& noise, Scalar noiseScale)
{
  const auto& from = geometry.planeOf(state);
  const auto& to = geometry.planeOf(measurement);

  Scalar start[3]{};
  surfaceToGlobal(from, state.loc0(), state.loc1(), start);
  const Scalar leaving = state.theta();
  const Scalar dQOverP = state.qOverP() - referenceQOverP;
  state.par[kQOverP] = referenceQOverP;
  // An intersection just outside the sensor's active area is still the point
  // to update at: the cluster is there, so the track was. Bounds are for
  // searches, not fits (`jacobianAnalytic` takes kOutOfBounds).
  const auto jacobianStatus = jacobianAnalytic<Scalar, Rule, Solve>(prop, state, from, to, direction, jacobian, &state);
  if (jacobianStatus != PropagationStatus::kOk) {
    return jacobianStatus;
  }
  // The filtered q/p again, so that the energy loss's fraction at the
  // reference (`stepNoise`) scales it too.
  state.par[kQOverP] += dQOverP;

  Scalar end[3]{};
  surfaceToGlobal(to, state.loc0(), state.loc1(), end);

  const auto step = stepMaterial(geometry, from, to, material, start, end);
  if (!step.found) {
    return PropagationStatus::kNoMaterial;
  }
  const Scalar momentum = prop.momentum(referenceQOverP);
  // The nearer half before the transport (StepMaterial), at the direction on
  // the surface the step leaves; the farther half after it.
  const auto before = scatteringNoise(scaledTheta0<Scaled>(highlandTheta0(step.nearerX2X0, step.nearerLogX2X0, momentum, prop.massGeV, prop.absCharge), noiseScale), leaving);
  const auto status = stepNoise<Scaled>(prop, step, momentum, before, direction, state, referenceQOverP, noise, noiseScale);
  if (status != PropagationStatus::kOk) {
    return status;
  }
  // The theta row of the Jacobian is the identity.
  state.par[kLoc0] += jacobian.loc0[kQOverP] * dQOverP;
  state.par[kLoc1] += jacobian.loc1[kQOverP] * dQOverP;
  state.par[kPhi] += jacobian.phi[kQOverP] * dQOverP;
  return PropagationStatus::kOk;
}

/// `filterStep`'s arithmetic. With \a Likelihood, \a logLikelihood gets
/// the measurement's log likelihood under the prediction, -(chi2 + ln det S)
/// / 2 without the constant: a Gaussian-sum component's weight
/// (GSFTrackFitter.h).
template <bool Scaled, CrossingRule Rule, bool Likelihood, typename Solve = ExactScalar, DetectorLayout Det, typename Material, typename Prop>
GPUdi() GPUAlwaysInline() PropagationStatus filterStepBody(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& state, CovarianceFactors& factors, Scalar& chi2, Scalar& referenceQOverP, [[maybe_unused]] Scalar noiseScale, [[maybe_unused]] Scalar* logLikelihood)
{
  SolenoidJacobian jacobian;
  StepNoise noise;
  Scalar nextReference = referenceQOverP;
  const auto predicted = predictStep<Scaled, Rule, Solve>(prop, geometry, material, measurement, direction, state, nextReference, jacobian, noise, noiseScale);
  if (predicted != PropagationStatus::kOk) {
    return predicted;
  }
  if (!transportFactors(jacobian, noise, factors)) {
    return PropagationStatus::kNotFinite;
  }
  [[maybe_unused]] Scalar determinant = 0.F;
  const auto update = updateWithMeasurement<Likelihood>(state, factors, measurement, Likelihood ? &determinant : nullptr);
  if (!updateUsable(update, factors)) {
    return PropagationStatus::kNotFinite;
  }
  referenceQOverP = nextReference;
  chi2 += update.chi2;
  if constexpr (Likelihood) {
    *logLikelihood = -0.5F * (update.chi2 + logOf(determinant));
  }
  return PropagationStatus::kOk;
}

/// The exact diffuse filter (de Jong, Ann. Stat. 19 (1991) 1073; Koopman,
/// JASA 92 (1997) 1630; Koopman and Durbin, J. Time Ser. Anal. 21 (2000)
/// 281): the covariance of a state some of whose parameters no measurement
/// has determined yet is P = kappa P_inf + P_star with kappa -> infinity,
/// P_inf spanning the undetermined directions. A filter pass starts so -- its
/// first measurement fixes the position, the angles and q/p are unknown --
/// and a kink frees the direction so. Each measured coordinate while P_inf
/// has rank determines one of its directions, with no chi2: nothing is there
/// yet to test it against. The limit is the pass without any prior: the
/// seed's values are the point the pass linearises about, never information,
/// and a pass's chi2 has exactly `2 N - 5` degrees of freedom, two fewer per
/// kink.
///
/// Both parts in square-root form: P_star = U D U^T as CovarianceFactors, D
/// zero along the undetermined directions, and P_inf = S S^T with `rank`
/// columns. Every operation is a weighted sum of outer products reduced by
/// modified weighted Gram-Schmidt (`transportFactors`'s method) or an
/// orthogonal rotation of S: no covariance is subtracted from another, and
/// float holds it.
struct DiffuseFactors {
  /// At most the angles and q/p: what a single measurement leaves undetermined.
  static constexpr int kMaxRank = 3;
  CovarianceFactors star;
  Scalar root[kNBoundParams][kMaxRank]{}; ///< S's columns
  int rank = 0;
};

/// The fewest measurements a pass's start takes (`startPass`): the first
/// fixes the position, the next two the direction and q/p.
inline constexpr int kPassStart = 3;

/// U D U^T = sum_k w_k y_k y_k^T over \a M columns of \a y (rows the five
/// parameters), into \a factors: modified weighted Gram-Schmidt from the last
/// row up (`transportFactors`), a zero pivot giving a zero column of U.
template <int M>
GPUdi() void factorColumns(Scalar (&y)[kNBoundParams][M], const Scalar (&w)[M], CovarianceFactors& factors)
{
  for (int j = kNBoundParams - 1; j >= 0; --j) {
    Scalar weighted[M];
    Scalar d = 0.F;
    for (int k = 0; k < M; ++k) {
      weighted[k] = w[k] * y[j][k];
      d += y[j][k] * weighted[k];
    }
    // A pivot below the smallest normal float is a zero pivot, as flushing
    // denormals to zero makes it: its reciprocal would overflow.
    if (d < FLT_MIN) {
      d = 0.F;
    }
    factors.diagonal[j] = d;
    const Scalar inverse = d > 0.F ? 1.F / d : 0.F;
    for (int i = 0; i < j; ++i) {
      Scalar dot = 0.F;
      for (int k = 0; k < M; ++k) {
        dot += y[i][k] * weighted[k];
      }
      const Scalar u = dot * inverse;
      factors.upper[strictUpperIndex(i, j)] = u;
      for (int k = 0; k < M; ++k) {
        y[i][k] -= u * y[j][k];
      }
    }
  }
}

/// U's columns into \a y's first five, D into \a w's.
template <int M>
GPUdi() void factorColumnsOf(const CovarianceFactors& factors, Scalar (&y)[kNBoundParams][M], Scalar (&w)[M])
{
  for (int k = 0; k < kNBoundParams; ++k) {
    for (int i = 0; i < k; ++i) {
      y[i][k] = factors.upper[strictUpperIndex(i, k)];
    }
    y[k][k] = 1.F;
    for (int i = k + 1; i < kNBoundParams; ++i) {
      y[i][k] = 0.F;
    }
    w[k] = factors.diagonal[k];
  }
}

/// x -> F x, F the step's Jacobian (SolenoidJacobian: theta and q/p kept).
GPUdi() void applyJacobian(const SolenoidJacobian& jac, Scalar (&x)[kNBoundParams])
{
  Scalar out[3]{};
  for (int j = 0; j < kNBoundParams; ++j) {
    out[0] += jac.loc0[j] * x[j];
    out[1] += jac.loc1[j] * x[j];
    out[2] += jac.phi[j] * x[j];
  }
  x[kLoc0] = out[0];
  x[kLoc1] = out[1];
  x[kPhi] = out[2];
}

/// The prediction of a diffuse covariance: P_star -> F (P_star + Q_before)
/// F^T + Q_after as `transportFactors` forms it, here with zero pivots
/// allowed, and P_inf -> F P_inf F^T, which no noise reaches.
GPUdi() void transportDiffuse(const SolenoidJacobian& jac, const StepNoise& noise, DiffuseFactors& c)
{
  constexpr int kColumns = kNBoundParams + 7;
  Scalar y[kNBoundParams][kColumns]{};
  Scalar w[kColumns]{};
  factorColumnsOf(c.star, y, w);
  for (int k = 0; k < kNBoundParams; ++k) {
    Scalar column[kNBoundParams];
    for (int i = 0; i < kNBoundParams; ++i) {
      column[i] = y[i][k];
    }
    applyJacobian(jac, column);
    for (int i = 0; i < kNBoundParams; ++i) {
      y[i][k] = column[i];
    }
  }
  // The nearer part's noise transported, F e_phi and F e_theta with the
  // positions' response at its lever (`StepNoise::lever`), its spread as
  // position noise alone; the farther half's and the straggling as they are.
  const Scalar spread = o2::gpu::GPUCommonMath::Sqrt(noise.spread);
  y[kLoc0][5] = jac.loc0[kPhi] * noise.lever;
  y[kLoc1][5] = jac.loc1[kPhi] * noise.lever;
  y[kPhi][5] = jac.phi[kPhi];
  w[5] = noise.phiBefore;
  y[kLoc0][6] = jac.loc0[kTheta] * noise.lever;
  y[kLoc1][6] = jac.loc1[kTheta] * noise.lever;
  y[kPhi][6] = jac.phi[kTheta];
  y[kTheta][6] = 1.F;
  w[6] = noise.thetaBefore;
  y[kPhi][7] = 1.F;
  w[7] = noise.phiAfter;
  y[kTheta][8] = 1.F;
  w[8] = noise.thetaAfter;
  y[kQOverP][9] = 1.F;
  w[9] = noise.qOverPAfter;
  y[kLoc0][10] = jac.loc0[kPhi] * spread;
  y[kLoc1][10] = jac.loc1[kPhi] * spread;
  w[10] = noise.phiBefore;
  y[kLoc0][11] = jac.loc0[kTheta] * spread;
  y[kLoc1][11] = jac.loc1[kTheta] * spread;
  w[11] = noise.thetaBefore;
  factorColumns(y, w, c.star);
  for (int k = 0; k < c.rank; ++k) {
    Scalar column[kNBoundParams];
    for (int i = 0; i < kNBoundParams; ++i) {
      column[i] = c.root[i][k];
    }
    applyJacobian(jac, column);
    for (int i = 0; i < kNBoundParams; ++i) {
      c.root[i][k] = column[i];
    }
  }
}

/// One coordinate h.x with variance \a r and residual \a v while P_inf has
/// rank: K = P_inf h / (h P_inf h), P_star -> (I - K h^T) P_star (I - K
/// h^T)^T + K K^T r (Koopman's diffuse update), and P_inf loses h's
/// direction: S is rotated (Householder) so that its first column alone
/// carries it, and that column is dropped. No chi2. False, nothing changed,
/// if h P_inf h is not positive: the coordinate carries no undetermined
/// direction. With \a Density, ln(h P_inf h) is added to \a density: the
/// coordinate's term of -2 ln of the diffuse likelihood (de Jong; Koopman
/// and Durbin), which compares filters that differ in their P_inf.
template <bool Density = false>
GPUdi() bool diffuseUpdate(DiffuseFactors& c, const Scalar (&h)[kNBoundParams], Scalar r, Scalar v, BoundTrackState& state, [[maybe_unused]] Scalar* density = nullptr)
{
  constexpr int kRank = DiffuseFactors::kMaxRank;
  // g = S^T h, K = S g / |g|^2.
  Scalar g[kRank]{};
  Scalar g2 = 0.F;
  for (int j = 0; j < c.rank; ++j) {
    for (int i = 0; i < kNBoundParams; ++i) {
      g[j] += h[i] * c.root[i][j];
    }
    g2 += g[j] * g[j];
  }
  if (!(g2 > 0.F)) {
    return false;
  }
  if constexpr (Density) {
    *density += logOf(g2);
  }
  Scalar gain[kNBoundParams]{};
  for (int i = 0; i < kNBoundParams; ++i) {
    for (int j = 0; j < c.rank; ++j) {
      gain[i] += c.root[i][j] * g[j];
    }
    gain[i] /= g2;
  }
  // H = I - 2 u u^T / |u|^2 with g^T H = (-+|g|, 0, ...): S H's first column
  // carries h, the others are orthogonal to it.
  Scalar u[kRank]{};
  for (int j = 0; j < c.rank; ++j) {
    u[j] = g[j];
  }
  const Scalar norm = o2::gpu::GPUCommonMath::Sqrt(g2);
  u[0] += g[0] >= 0.F ? norm : -norm;
  Scalar u2 = 0.F;
  for (int j = 0; j < c.rank; ++j) {
    u2 += u[j] * u[j];
  }
  for (auto& row : c.root) {
    Scalar su = 0.F;
    for (int j = 0; j < c.rank; ++j) {
      su += row[j] * u[j];
    }
    const Scalar scale = 2.F * su / u2;
    for (int j = 0; j < c.rank; ++j) {
      row[j] -= scale * u[j];
    }
    for (int j = 0; j + 1 < c.rank; ++j) {
      row[j] = row[j + 1];
    }
    row[c.rank - 1] = 0.F;
  }
  --c.rank;

  constexpr int kColumns = kNBoundParams + 1;
  Scalar y[kNBoundParams][kColumns]{};
  Scalar w[kColumns]{};
  factorColumnsOf(c.star, y, w);
  for (int k = 0; k < kNBoundParams; ++k) {
    Scalar hy = 0.F;
    for (int i = 0; i < kNBoundParams; ++i) {
      hy += h[i] * y[i][k];
    }
    for (int i = 0; i < kNBoundParams; ++i) {
      y[i][k] -= gain[i] * hy;
    }
  }
  for (int i = 0; i < kNBoundParams; ++i) {
    y[i][kNBoundParams] = gain[i];
  }
  w[kNBoundParams] = r;
  factorColumns(y, w, c.star);
  detail::applySpatialGain(state, gain, v);
  return true;
}

/// One coordinate h.x with variance \a r and residual \a v: the diffuse
/// update if P_inf projects onto it (`diffuseUpdate`); otherwise -- P_inf
/// gone, or orthogonal to h -- the filter's (Bierman, as
/// `updateWithMeasurement`) on P_star, P_inf kept, its chi2 added to
/// \a chi2. False for an innovation that is not positive. With \a Density,
/// the coordinate's term of -2 ln of the diffuse likelihood is added to
/// \a density: `diffuseUpdate`'s, or ln F + v^2 / F of the filter's
/// innovation F.
template <bool Density = false>
GPUdi() bool updateCoordinate(DiffuseFactors& c, const Scalar (&h)[kNBoundParams], Scalar r, Scalar v, BoundTrackState& state, Scalar& chi2, [[maybe_unused]] Scalar* density = nullptr)
{
  if (c.rank > 0 && diffuseUpdate<Density>(c, h, r, v, state, density)) {
    return true;
  }
  Scalar gain[kNBoundParams]{};
  const Scalar innovation = detail::updateUD(c.star.upper, c.star.diagonal, h, r, gain);
  if (!(innovation > 0.F)) {
    return false;
  }
  detail::applySpatialGain(state, gain, v);
  chi2 += v * v / innovation;
  if constexpr (Density) {
    *density += logOf(innovation) + (v * v / innovation);
  }
  return true;
}

/// A measurement's two coordinates, whitened as `updateWithMeasurement`
/// whitens them: the first raw, the second conditioned on it. While P_inf
/// has rank the first is the one P_inf projects more onto, so that it
/// carries a direction of P_inf; the joint update does not depend on the
/// order. With \a Density, each coordinate's term of -2 ln of the diffuse
/// likelihood is added to \a density (`updateCoordinate`).
template <bool Density = false>
GPUdi() bool updateDiffuse(DiffuseFactors& c, const Measurement& measurement, BoundTrackState& state, Scalar& chi2, [[maybe_unused]] Scalar* density = nullptr)
{
  int a = kLoc0;
  if (c.rank > 0) {
    Scalar projected0 = 0.F;
    Scalar projected1 = 0.F;
    for (int j = 0; j < c.rank; ++j) {
      projected0 += c.root[kLoc0][j] * c.root[kLoc0][j];
      projected1 += c.root[kLoc1][j] * c.root[kLoc1][j];
    }
    a = projected1 > projected0 ? kLoc1 : kLoc0;
  }
  const int b = a == kLoc0 ? kLoc1 : kLoc0;
  const Scalar value[2] = {measurement.loc0, measurement.loc1};
  const Scalar variance[2] = {measurement.varLoc0, measurement.varLoc1};
  Scalar h[kNBoundParams]{};
  h[a] = 1.F;
  if (!updateCoordinate<Density>(c, h, variance[a], value[a] - state.par[a], state, chi2, density)) {
    return false;
  }
  const Scalar correlation = measurement.covLoc01 / variance[a];
  h[a] = -correlation;
  h[b] = 1.F;
  const Scalar residual = (value[b] - correlation * value[a]) - (state.par[b] - correlation * state.par[a]);
  return updateCoordinate<Density>(c, h, variance[b] - correlation * measurement.covLoc01, residual, state, chi2, density);
}

/// `diffuseStep`'s arithmetic. With \a Density, the update's terms of -2 ln
/// of the diffuse likelihood are added to \a density (`updateDiffuse`): a
/// Gaussian-sum branch's weight (GSFTrackFitter.h).
template <CrossingRule Rule, bool Density, typename Solve = ExactScalar, DetectorLayout Det, typename Material, typename Prop>
GPUdi() GPUAlwaysInline() PropagationStatus diffuseStepBody(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& state, DiffuseFactors& c, Scalar& chi2, Scalar& referenceQOverP, [[maybe_unused]] Scalar* density)
{
  SolenoidJacobian jacobian;
  StepNoise noise;
  Scalar nextReference = referenceQOverP;
  const auto predicted = predictStep<false, Rule, Solve>(prop, geometry, material, measurement, direction, state, nextReference, jacobian, noise, 1.F);
  if (predicted != PropagationStatus::kOk) {
    return predicted;
  }
  transportDiffuse(jacobian, noise, c);
  if (!updateDiffuse<Density>(c, measurement, state, chi2, density)) {
    return PropagationStatus::kNotFinite;
  }
  if (c.rank > 0) {
    o2::math_utils::detail::bringToPMPi<Scalar>(state.par[kPhi]);
    referenceQOverP = nextReference;
    return PropagationStatus::kOk;
  }
  // Every variance a sum of positive terms once each direction is determined.
  if (!positiveDiagonal(c.star)) {
    return PropagationStatus::kNotFinite;
  }
  detail::canonicalizeDirection(state, c.star);
  referenceQOverP = nextReference;
  return PropagationStatus::kOk;
}

/// A pass's first measurement as the exact diffuse filter takes it
/// (DiffuseFactors): the position becomes the measurement's, with its
/// covariance as P_star; the angles and q/p are P_inf, undetermined, and keep
/// \a state's values only as the point the pass linearises about. False if
/// the measurement's covariance is not positive definite.
GPUdi() bool startDiffuse(const Measurement& measurement, BoundTrackState& state, DiffuseFactors& c)
{
  state.par[kLoc0] = measurement.loc0;
  state.par[kLoc1] = measurement.loc1;
  c = DiffuseFactors{};
  // As factors: loc1's variance, then loc0's given loc1.
  c.star.diagonal[kLoc1] = measurement.varLoc1;
  c.star.upper[strictUpperIndex(kLoc0, kLoc1)] = measurement.covLoc01 / measurement.varLoc1;
  c.star.diagonal[kLoc0] = measurement.varLoc0 - (measurement.covLoc01 * measurement.covLoc01 / measurement.varLoc1);
  // Any basis of the angles and q/p gives the same limit.
  c.root[kPhi][0] = 1.F;
  c.root[kTheta][1] = 1.F;
  c.root[kQOverP][2] = 1.F;
  c.rank = 3;
  return c.star.diagonal[kLoc1] > 0.F && c.star.diagonal[kLoc0] > 0.F;
}

/// \a measurements without the one at \a skip, for a helix through the
/// others (`leaveOneOutQOverP`).
template <typename Chain>
struct SkipChain {
  const Chain* chain;
  int skip;
  GPUdi() auto operator[](int i) const { return (*chain)[i < skip ? i : i + 1]; }
};

/// The q/p of the helix through every measurement of the chain but the one
/// at \a index (`makeChainSeed`), on the chain's first measurement.
template <typename Prop, DetectorLayout Det, typename Chain>
GPUdi() Scalar leaveOneOutQOverP(const Prop& prop, const TrackingGeometry<Det>& geometry, const Chain& measurements, int n, int index)
{
  return makeChainSeed(prop, geometry, SkipChain<Chain>{.chain = &measurements, .skip = index}, n - 1).qOverP();
}

/// A pass's reference momentum (`predictStep`) as the pass goes. Outward,
/// the step to a measurement takes the helix through every other measurement
/// of the chain (`leaveOneOutQOverP`), carried by \a lossScale, the energy
/// loss's fraction of the momentum so far: a measurement then never sets the
/// scattering its own step is tested with. A wrong measurement drawn into
/// the helix bends it, and the low momentum's scattering would open the gate
/// wide enough to take that measurement. A chain of the fewest measurements a
/// pass takes leaves no helix without one of them; it, and the inward pass,
/// take \a qOverP, the start state's q/p (the outward pass's end for the
/// inward one) carried through the energy loss.
struct PassReference {
  Scalar qOverP = 0.F;
  Scalar lossScale = 1.F;
};

/// Whether the pass over \a n measurements in \a direction takes the
/// leave-one-out reference (PassReference).
GPUdi() bool leavesOneOut(int n, PropagationDirection direction) { return direction == PropagationDirection::kForward && n > kPassStart; }

/// The helix q/p the step to the measurement at \a index takes if the pass
/// leaves one out (`leaveOneOutQOverP`), else 0: the same for every state a
/// pass carries (a Gaussian-sum fit's components).
template <typename Prop, DetectorLayout Det, typename Chain>
GPUdi() Scalar leftOutQOverP(const Prop& prop, const TrackingGeometry<Det>& geometry, const Chain& measurements, int n, int index, PropagationDirection direction)
{
  return leavesOneOut(n, direction) ? leaveOneOutQOverP(prop, geometry, measurements, n, index) : 0.F;
}

/// The reference momentum of a step (PassReference); \a leftOut() gives its
/// leave-one-out helix q/p (`leaveOneOutQOverP`), called only if the pass
/// takes it.
template <typename LeftOut>
GPUdi() Scalar referenceFor(const PassReference& reference, int n, PropagationDirection direction, const LeftOut& leftOut)
{
  return leavesOneOut(n, direction) ? leftOut() * reference.lossScale : reference.qOverP;
}

/// \a reference after a step that took \a before and carried it through its
/// energy loss to \a after.
GPUdi() void advanceReference(PassReference& reference, int n, PropagationDirection direction, Scalar before, Scalar after)
{
  if (!leavesOneOut(n, direction)) {
    reference.qOverP = after;
  } else if (before != 0.F) {
    reference.lossScale *= after / before;
  }
}

/// The material of a fit's span: the link material's steps between
/// consecutive measurements of \a measurements, the first \a n, which are
/// what the fit's steps read; the fit's material (`SpanMaterial`) for both
/// passes.
template <DetectorLayout Det, typename Material, typename Chain>
GPUdi() SpanMaterial<Material> spanMaterial(const TrackingGeometry<Det>& geometry, const Material& /*material*/, const Chain& measurements, int n)
{
  SpanMaterial<Material> span{};
  Scalar previous[3]{};
  LayerID previousLayer = kInvalidLayerID;
  for (int i = 0; i < n; ++i) {
    const auto measurement = measurements[i];
    Scalar point[3];
    const auto& plane = geometry.planeOf(measurement);
    surfaceToGlobal(plane, measurement.loc0, measurement.loc1, point);
    // A layer's two sensors: its material is in the steps to and from it.
    bool sameLayer = false; // NOLINT(misc-const-correctness): set under if constexpr, which the check does not see
    if constexpr (IsLayerMaterial<Material>::value) {
      sameLayer = previousLayer == plane.layer;
    }
    if (i > 0 && !sameLayer) {
      // A step the link material holds none of fails the fit (kNoMaterial).
      const auto step = geometry.material().step(previousLayer, plane.layer, previous, point);
      span.spanX2X0 += step.before + step.after;
    }
    previousLayer = plane.layer;
    previous[0] = point[0];
    previous[1] = point[1];
    previous[2] = point[2];
  }
  return span;
}

} // namespace detail

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRACKFITTERUTILS_H_
