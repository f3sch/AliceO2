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

/// \file TrackFitter.h
/// \brief Independent outward and inward Kalman filters over a chain.

#ifndef GTRK_TRACKING_TRACKFITTER_H_
#define GTRK_TRACKING_TRACKFITTER_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Tracking/CovarianceTransport.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/MeasurementUpdate.h"
#include "GenericTracking/Tracking/Seed.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
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

/// A fit's per-measurement chi2 increments, for QA: called once per update
/// with the pass's direction, the measurement's index in the chain and its
/// increment (two dof). The default does nothing.
struct NoFitObserver {
  GPUhdi() GPUAlwaysInline() void operator()(PropagationDirection /*direction*/, int /*index*/, Scalar /*chi2*/) const {}
};

/// The mass a fit step's energy loss takes: the propagator's hypothesis,
/// the one its scattering takes. A QA oracle's propagator type overloads it
/// to set the two apart (SpeciesOracle.h).
GPUhdi() GPUAlwaysInline() Scalar lossMassGeV(const Propagator& prop) { return prop.massGeV; }

/// The bremsstrahlung a fit step adds to its energy loss over \a x2X0
/// radiation lengths: none for the tracking's propagator. A QA oracle's
/// propagator type overloads it (SpeciesOracle.h) to change \a state's q/p
/// and add to the step's q/p variance \a varQOverP.
GPUhdi() GPUAlwaysInline() void applyRadiation(const Propagator& /*prop*/, Scalar /*x2X0*/, PropagationDirection /*direction*/, BoundTrackState& /*state*/, Scalar& /*varQOverP*/) {}

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
struct StepMaterial {
  Scalar nearerX2X0 = 0.F;  ///< x/X0 from the start to the midpoint
  Scalar fartherX2X0 = 0.F; ///< x/X0 from the midpoint to the end
  Scalar xRho = 0.F;        ///< [g/cm^2] over the whole step
};

template <typename Material>
GPUdi() StepMaterial stepMaterial(const Material& material, const Scalar start[3], const Scalar end[3])
{
  const Scalar mid[3] = {0.5F * (start[0] + end[0]), 0.5F * (start[1] + end[1]), 0.5F * (start[2] + end[2])};
  const auto nearer = material.getMatBudget(start[0], start[1], start[2], mid[0], mid[1], mid[2]);
  const auto farther = material.getMatBudget(mid[0], mid[1], mid[2], end[0], end[1], end[2]);
  return {.nearerX2X0 = nearer.meanX2X0, .fartherX2X0 = farther.meanX2X0, .xRho = (nearer.meanRho * nearer.length) + (farther.meanRho * farther.length)};
}

/// The mass hypothesis's energy loss over a step's \a xRho [g/cm^2], O2's
/// way: taken at \a state's momentum and applied to it, the straggling's q/p
/// variance returned in \a varQOverP for the step's noise (StepNoise). Lost
/// along the momentum, regained against it. False if the particle stopped.
/// The mass is `lossMassGeV(prop)`, the propagator's hypothesis.
template <typename Prop>
GPUhdi() GPUAlwaysInline() bool applyEnergyLoss(const Prop& prop, Scalar xRho, PropagationDirection direction, BoundTrackState& state, Scalar& varQOverP)
{
  varQOverP = 0.F;
  const auto loss = energyLoss(state.qOverP(), lossMassGeV(prop), direction == PropagationDirection::kForward ? -xRho : xRho);
  if (!loss.ok) {
    return false;
  }
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
/// with, and the energy loss over \a step at its momentum, applied to
/// \a state, its straggling's variance in \a noise; \a before is the nearer
/// half's (StepMaterial). With \a Scaled, \a noiseScale multiplies the
/// farther half's theta0 (`scaledTheta0`).
template <bool Scaled, typename Prop>
GPUdi() PropagationStatus stepNoise(const Prop& prop, const StepMaterial& step, Scalar betaP, const ScatteringNoise& before, PropagationDirection direction, BoundTrackState& state, StepNoise& noise, Scalar noiseScale)
{
  const auto after = scatteringNoise(scaledTheta0<Scaled>(highlandTheta0(step.fartherX2X0, betaP), noiseScale), state.theta());
  noise = StepNoise{.phiBefore = before.phi, .thetaBefore = before.theta, .phiAfter = after.phi, .thetaAfter = after.theta};
  if (!applyEnergyLoss(prop, step.xRho, direction, state, noise.qOverPAfter)) {
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
GPUdi() PropagationStatus transportStep(const Prop& prop, const StepMaterial& step, Scalar betaP, const ScatteringNoise& before, PropagationDirection direction, const SolenoidJacobian& jacobian, BoundTrackState& state, CovarianceFactors& factors, Scalar noiseScale)
{
  StepNoise noise;
  const auto status = stepNoise<Scaled>(prop, step, betaP, before, direction, state, noise, noiseScale);
  if (status != PropagationStatus::kOk) {
    return status;
  }
  if (!transportFactors(jacobian, noise, factors)) {
    return PropagationStatus::kNotFinite;
  }
  return PropagationStatus::kOk;
}

/// A filter step's prediction: \a state carried to \a measurement's plane,
/// linearised about \a reference -- the pass's reference trajectory on the
/// surface the step leaves, not the filtered state. The reference is
/// propagated exactly, with the mean energy loss, and the Jacobian, the
/// material and the scattering are taken along it; the state is predicted as
/// `f(reference) + F (state - reference)`, so the transport and its noise do
/// not depend on the measurements the pass has taken. The propagated
/// reference is \a nextReference, \a jacobian and \a noise what carries a
/// covariance along (`transportFactors`). \a Scaled, the step's scattering
/// theta0 is multiplied by \a noiseScale, both halves.
template <bool Scaled, CrossingRule Rule, DetectorLayout Det, typename Material, typename Prop>
GPUdi() PropagationStatus predictStep(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& state, const BoundTrackState& reference, BoundTrackState& nextReference, SolenoidJacobian& jacobian, StepNoise& noise, Scalar noiseScale)
{
  BoundTrackState prediction = reference;
  const auto& from = geometry.planeOf(prediction);
  const auto& to = geometry.planeOf(measurement);

  Scalar start[3]{};
  surfaceToGlobal(from, prediction.loc0(), prediction.loc1(), start);
  const Scalar leaving = prediction.theta();
  // An intersection just outside the sensor's active area is still the point
  // to update at: the cluster is there, so the track was. Bounds are for
  // searches, not fits (`jacobianAnalytic` takes kOutOfBounds).
  const auto jacobianStatus = jacobianAnalytic<Scalar, Rule>(prop, prediction, from, to, direction, jacobian, &prediction);
  if (jacobianStatus != PropagationStatus::kOk) {
    return jacobianStatus;
  }

  Scalar end[3]{};
  surfaceToGlobal(to, prediction.loc0(), prediction.loc1(), end);

  const auto step = stepMaterial(material, start, end);
  const Scalar momentum = 1.F / o2::gpu::GPUCommonMath::Abs(prediction.qOverP());
  const Scalar betaP = betaMomentum(momentum, prop.massGeV);
  // The nearer half before the transport (StepMaterial), at the direction on
  // the surface the step leaves; the farther half after it.
  const auto before = scatteringNoise(scaledTheta0<Scaled>(highlandTheta0(step.nearerX2X0, betaP), noiseScale), leaving);
  const auto noiseStatus = stepNoise<Scaled>(prop, step, betaP, before, direction, prediction, noise, noiseScale);
  if (noiseStatus != PropagationStatus::kOk) {
    return noiseStatus;
  }

  nextReference = prediction;
  Scalar delta[kNBoundParams];
  for (int i = 0; i < kNBoundParams; ++i) {
    delta[i] = state.par[i] - reference.par[i];
  }
  o2::math_utils::detail::bringToPMPi<Scalar>(delta[kPhi]);
  const Scalar* rows[3] = {jacobian.loc0, jacobian.loc1, jacobian.phi};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < kNBoundParams; ++j) {
      prediction.par[i] += rows[i][j] * delta[j];
    }
  }
  // The theta and q/p rows of the Jacobian are the identity.
  prediction.par[kTheta] += delta[kTheta];
  prediction.par[kQOverP] += delta[kQOverP];
  state = prediction;
  return PropagationStatus::kOk;
}

/// One step of a filter pass: \a state predicted on \a measurement's plane
/// (`predictStep`) and updated with it. The propagated reference, its
/// position moved onto the measurement after the update (as O2's ITS fit
/// shifts its reference to the cluster), is \a nextReference, the next
/// step's: its direction and momentum stay the propagated ones, never
/// updated by the filter. The covariance is \a factors, the pass's
/// (CovarianceFactors): `state.cov` is neither read nor written. \a Scaled,
/// the step's scattering theta0 is multiplied by \a noiseScale, both halves:
/// the kink allowance's last step (`allowanceStep`); the nominal step is its
/// own instantiation. \a Rule is the propagation's choice of crossing
/// (`CrossingRule`), and with it every pass below. Never inlined on the
/// device: one copy for each instantiation, not one in every fit of every
/// stage variant.
template <bool Scaled = false, CrossingRule Rule = CrossingRule::kFirstAhead, DetectorLayout Det, typename Material, typename Prop>
GPUdnii() PropagationStatus filterStep(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& state, CovarianceFactors& factors, Scalar& chi2, const BoundTrackState& reference, BoundTrackState& nextReference, [[maybe_unused]] Scalar noiseScale = 1.F)
{
  SolenoidJacobian jacobian;
  StepNoise noise;
  const auto predicted = predictStep<Scaled, Rule>(prop, geometry, material, measurement, direction, state, reference, nextReference, jacobian, noise, noiseScale);
  if (predicted != PropagationStatus::kOk) {
    return predicted;
  }
  if (!transportFactors(jacobian, noise, factors)) {
    return PropagationStatus::kNotFinite;
  }
  const auto update = updateWithMeasurement(state, factors, measurement);
  if (!updateUsable(update, factors)) {
    return PropagationStatus::kNotFinite;
  }
  nextReference.par[kLoc0] = measurement.loc0;
  nextReference.par[kLoc1] = measurement.loc1;
  chi2 += update.chi2;
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
  constexpr int M = kNBoundParams + 5;
  Scalar y[kNBoundParams][M]{};
  Scalar w[M]{};
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
  // The nearer half's noise transported, F e_phi and F e_theta; the farther
  // half's and the straggling as they are.
  y[kLoc0][5] = jac.loc0[kPhi];
  y[kLoc1][5] = jac.loc1[kPhi];
  y[kPhi][5] = jac.phi[kPhi];
  w[5] = noise.phiBefore;
  y[kLoc0][6] = jac.loc0[kTheta];
  y[kLoc1][6] = jac.loc1[kTheta];
  y[kPhi][6] = jac.phi[kTheta];
  y[kTheta][6] = 1.F;
  w[6] = noise.thetaBefore;
  y[kPhi][7] = 1.F;
  w[7] = noise.phiAfter;
  y[kTheta][8] = 1.F;
  w[8] = noise.thetaAfter;
  y[kQOverP][9] = 1.F;
  w[9] = noise.qOverPAfter;
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
/// direction.
GPUdi() bool diffuseUpdate(DiffuseFactors& c, const Scalar (&h)[kNBoundParams], Scalar r, Scalar v, BoundTrackState& state)
{
  constexpr int R = DiffuseFactors::kMaxRank;
  // g = S^T h, K = S g / |g|^2.
  Scalar g[R]{};
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
  Scalar gain[kNBoundParams]{};
  for (int i = 0; i < kNBoundParams; ++i) {
    for (int j = 0; j < c.rank; ++j) {
      gain[i] += c.root[i][j] * g[j];
    }
    gain[i] /= g2;
  }
  // H = I - 2 u u^T / |u|^2 with g^T H = (-+|g|, 0, ...): S H's first column
  // carries h, the others are orthogonal to it.
  Scalar u[R]{};
  for (int j = 0; j < c.rank; ++j) {
    u[j] = g[j];
  }
  const Scalar norm = o2::gpu::GPUCommonMath::Sqrt(g2);
  u[0] += g[0] >= 0.F ? norm : -norm;
  Scalar u2 = 0.F;
  for (int j = 0; j < c.rank; ++j) {
    u2 += u[j] * u[j];
  }
  for (int i = 0; i < kNBoundParams; ++i) {
    Scalar su = 0.F;
    for (int j = 0; j < c.rank; ++j) {
      su += c.root[i][j] * u[j];
    }
    const Scalar scale = 2.F * su / u2;
    for (int j = 0; j < c.rank; ++j) {
      c.root[i][j] -= scale * u[j];
    }
    for (int j = 0; j + 1 < c.rank; ++j) {
      c.root[i][j] = c.root[i][j + 1];
    }
    c.root[i][c.rank - 1] = 0.F;
  }
  --c.rank;

  constexpr int M = kNBoundParams + 1;
  Scalar y[kNBoundParams][M]{};
  Scalar w[M]{};
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
/// \a chi2. False for an innovation that is not positive.
GPUdi() bool updateCoordinate(DiffuseFactors& c, const Scalar (&h)[kNBoundParams], Scalar r, Scalar v, BoundTrackState& state, Scalar& chi2)
{
  if (c.rank > 0 && diffuseUpdate(c, h, r, v, state)) {
    return true;
  }
  Scalar gain[kNBoundParams]{};
  const Scalar innovation = detail::updateUD(c.star.upper, c.star.diagonal, h, r, gain);
  if (!(innovation > 0.F)) {
    return false;
  }
  detail::applySpatialGain(state, gain, v);
  chi2 += v * v / innovation;
  return true;
}

/// A measurement's two coordinates, whitened as `updateWithMeasurement`
/// whitens them: the first raw, the second conditioned on it. While P_inf
/// has rank the first is the one P_inf projects more onto, so that it
/// carries a direction of P_inf; the joint update does not depend on the
/// order.
GPUdi() bool updateDiffuse(DiffuseFactors& c, const Measurement& measurement, BoundTrackState& state, Scalar& chi2)
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
  if (!updateCoordinate(c, h, variance[a], value[a] - state.par[a], state, chi2)) {
    return false;
  }
  const Scalar correlation = measurement.covLoc01 / variance[a];
  h[a] = -correlation;
  h[b] = 1.F;
  const Scalar residual = (value[b] - correlation * value[a]) - (state.par[b] - correlation * state.par[a]);
  return updateCoordinate(c, h, variance[b] - correlation * measurement.covLoc01, residual, state, chi2);
}

/// A filter step (`filterStep`) on a diffuse covariance \a c: the prediction
/// the same, the transport `transportDiffuse`, the update `updateDiffuse`;
/// \a chi2 gets only what the update tests. Once P_inf is gone, \a c's
/// factors are the pass's. Never inlined on the device.
template <CrossingRule Rule = CrossingRule::kFirstAhead, DetectorLayout Det, typename Material, typename Prop>
GPUdnii() PropagationStatus diffuseStep(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& state, DiffuseFactors& c, Scalar& chi2, const BoundTrackState& reference, BoundTrackState& nextReference)
{
  SolenoidJacobian jacobian;
  StepNoise noise;
  const auto predicted = predictStep<false, Rule>(prop, geometry, material, measurement, direction, state, reference, nextReference, jacobian, noise, 1.F);
  if (predicted != PropagationStatus::kOk) {
    return predicted;
  }
  transportDiffuse(jacobian, noise, c);
  if (!updateDiffuse(c, measurement, state, chi2)) {
    return PropagationStatus::kNotFinite;
  }
  nextReference.par[kLoc0] = measurement.loc0;
  nextReference.par[kLoc1] = measurement.loc1;
  if (c.rank > 0) {
    o2::math_utils::detail::bringToPMPi<Scalar>(state.par[kPhi]);
    return PropagationStatus::kOk;
  }
  // Every variance a sum of positive terms once each direction is determined.
  if (!positiveDiagonal(c.star)) {
    return PropagationStatus::kNotFinite;
  }
  detail::canonicalizeDirection(state, c.star);
  return PropagationStatus::kOk;
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
/// The reference is the step's either way (`filterStep`).
template <CrossingRule Rule = CrossingRule::kFirstAhead, DetectorLayout Det, typename Material, typename Prop>
GPUdi() PropagationStatus kinkStep(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& state, CovarianceFactors& factors, Scalar& chi2, const KinkBudget& budget, bool followed, int& kinks, const BoundTrackState& reference, BoundTrackState& nextReference)
{
  if (kinks >= budget.max || !followed) {
    return filterStep<false, Rule>(prop, geometry, material, measurement, direction, state, factors, chi2, reference, nextReference);
  }
  const BoundTrackState keptState = state;
  const CovarianceFactors keptFactors = factors;
  const Scalar keptChi2 = chi2;
  auto status = filterStep<false, Rule>(prop, geometry, material, measurement, direction, state, factors, chi2, reference, nextReference);
  if (status != PropagationStatus::kOk || !(chi2 - keptChi2 > budget.gate)) {
    return status;
  }
  state = keptState;
  chi2 = keptChi2;
  DiffuseFactors freed;
  freeDirection(keptFactors, freed);
  status = diffuseStep<Rule>(prop, geometry, material, measurement, direction, state, freed, chi2, reference, nextReference);
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
template <CrossingRule Rule = CrossingRule::kFirstAhead, DetectorLayout Det, typename Material>
GPUdi() PropagationStatus allowanceStep(const Propagator& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, BoundTrackState& state, CovarianceFactors& factors, Scalar& chi2, const KinkBudget& budget, Scalar noiseScale, int& kinks, bool& refused, const BoundTrackState& reference, BoundTrackState& nextReference)
{
  refused = false;
  if (kinks >= budget.max) {
    return filterStep<false, Rule>(prop, geometry, material, measurement, direction, state, factors, chi2, reference, nextReference);
  }
  const BoundTrackState keptState = state;
  const CovarianceFactors keptFactors = factors;
  const Scalar keptChi2 = chi2;
  auto status = filterStep<false, Rule>(prop, geometry, material, measurement, direction, state, factors, chi2, reference, nextReference);
  if (status != PropagationStatus::kOk || !(chi2 - keptChi2 > budget.gate)) {
    return status;
  }
  state = keptState;
  factors = keptFactors;
  chi2 = keptChi2;
  status = filterStep<true, Rule>(prop, geometry, material, measurement, direction, state, factors, chi2, reference, nextReference, noiseScale);
  if (status == PropagationStatus::kOk) {
    ++kinks;
    refused = !(chi2 - keptChi2 <= budget.gate);
  }
  return status;
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

/// A pass's start, with nothing known but what its measurements determine
/// (`startDiffuse`, then `diffuseStep`) until nothing is undetermined: the
/// first fixes the position, the next two the direction and q/p, the third's
/// second coordinate the first the chi2 tests -- three measurements, unless a
/// coordinate carries no undetermined direction. Ends with \a state, its
/// \a factors, \a chi2 and the next step's \a reference on the last
/// measurement it took, their count in \a taken; \a observe sees each
/// measurement's increment. A pass whose measurements leave a direction
/// undetermined determines no state.
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Prop, DetectorLayout Det, typename Material, typename Chain, typename Observer>
GPUdi() PropagationStatus startPass(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, PropagationDirection direction, BoundTrackState& state, CovarianceFactors& factors, Scalar& chi2, BoundTrackState& reference, int& taken, std::int16_t& failedMeasurement, const Observer& observe)
{
  const int first = direction == PropagationDirection::kForward ? 0 : n - 1;
  const int increment = direction == PropagationDirection::kForward ? 1 : -1;
  DiffuseFactors c;
  taken = 0;
  if (n < kPassStart || !startDiffuse(measurements[first], state, c)) {
    failedMeasurement = static_cast<std::int16_t>(first);
    return PropagationStatus::kNotFinite;
  }
  chi2 = 0.F;
  reference = state;
  observe(direction, first, 0.F);
  for (taken = 1; taken < n && c.rank > 0; ++taken) {
    const int index = first + (taken * increment);
    const Scalar before = chi2;
    const auto& stepProp = stepPropagator(prop, geometry, measurements, index - increment, index);
    BoundTrackState nextReference;
    const auto status = diffuseStep<Rule>(stepProp, geometry, material, measurements[index], direction, state, c, chi2, reference, nextReference);
    if (status != PropagationStatus::kOk) {
      failedMeasurement = static_cast<std::int16_t>(index);
      return status;
    }
    reference = nextReference;
    observe(direction, index, chi2 - before);
  }
  if (c.rank != 0) {
    failedMeasurement = static_cast<std::int16_t>(first + ((n - 1) * increment));
    return PropagationStatus::kNotFinite;
  }
  factors = c.star;
  return PropagationStatus::kOk;
}

/// A filter pass from \a seed over the chain: its start (`startPass`), then
/// the filter. Its covariance is carried as factors from the first update to
/// the last (CovarianceFactors) and written to \a result's `cov` once, at
/// the end. Its reference trajectory starts from \a seed's values on the
/// first measurement (`filterStep`); a kink is not tested before the start
/// has determined the state.
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Prop, DetectorLayout Det, typename Material, typename Chain, typename Observer>
GPUdi() PropagationStatus runPass(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, PropagationDirection direction, BoundTrackState seed, BoundTrackState& result, Scalar& chi2, std::int16_t& failedMeasurement, const KinkBudget& passBudget, int& kinks, std::int16_t& kinkMeasurement, const Observer& observe)
{
  const auto budget = kinkBudgetFor(passBudget, n);
  const int first = direction == PropagationDirection::kForward ? 0 : n - 1;
  const int increment = direction == PropagationDirection::kForward ? 1 : -1;

  BoundTrackState state = seed;
  CovarianceFactors factors;
  BoundTrackState reference;
  int taken = 0;
  const auto started = startPass<Rule>(prop, geometry, material, measurements, n, direction, state, factors, chi2, reference, taken, failedMeasurement, observe);
  if (started != PropagationStatus::kOk) {
    return started;
  }
  for (int step = taken; step < n; ++step) {
    const int index = first + step * increment;
    const Scalar before = chi2;
    const int kinksBefore = kinks;
    const auto& stepProp = stepPropagator(prop, geometry, measurements, index - increment, index);
    BoundTrackState nextReference;
    const auto status = kinkStep<Rule>(stepProp, geometry, material, measurements[index], direction, state, factors, chi2, budget, step < n - 1, kinks, reference, nextReference);
    if (status != PropagationStatus::kOk) {
      failedMeasurement = static_cast<std::int16_t>(index);
      return status;
    }
    reference = nextReference;
    if (kinks != kinksBefore) {
      kinkMeasurement = static_cast<std::int16_t>(index);
    }
    observe(direction, index, chi2 - before);
  }
  covarianceFromFactors(factors, state.cov);
  result = state;
  return PropagationStatus::kOk;
}

/// The inward pass of `fitMeasurements`: from the outward pass's end state
/// as its reference, as ACTS's reverse filter starts from the forward one's
/// (the whole-chain helix carried there unfitted misses the outermost rows
/// at low momentum), independent of the outward pass's; its values only, so
/// no measurement counts twice.
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Prop, DetectorLayout Det, typename Material, typename Chain, typename Observer>
GPUdi() void inwardPass(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, const KinkBudget& kinkBudget, const Observer& observe, TrackFitResult& fit)
{
  int kinks = 0;
  const auto status = runPass<Rule>(prop, geometry, material, measurements, n, PropagationDirection::kBackward, fit.outwardState, fit.state, fit.chi2, fit.failedMeasurement, kinkBudget, kinks, fit.kinkMeasurement, observe);
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
/// measurements determine (`startPass`). The outward one linearises about a
/// helix fitted through every measurement of the chain (`makeChainSeed`),
/// the inward one about the outward filter's end state, as ACTS's reverse
/// filter; each about a reference trajectory of its own, started there and
/// never updated by its filter (`filterStep`). The inward filter takes no
/// covariance from the outward one, so no measurement's weight is counted
/// twice in the returned state.
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
  int outwardKinks = 0;
  std::int16_t outwardKinkMeasurement = -1;
  auto status = detail::runPass<Rule>(prop, geometry, material, measurements, n, PropagationDirection::kForward, seed, fit.outwardState, fit.outwardChi2, fit.failedMeasurement, kinkBudget, outwardKinks, outwardKinkMeasurement, observe);
  fit.outwardKinks = static_cast<std::uint8_t>(outwardKinks);
  if (status != PropagationStatus::kOk) {
    fit.status = status;
    fit.failedPass = FitPass::kOutward;
    return fit;
  }
  detail::inwardPass<Rule>(prop, geometry, material, measurements, n, kinkBudget, observe, fit);
  return fit;
}

/// `fitMeasurements` given its outward pass's end state, \a outwardState
/// -- as the road fit keeps it for an accepted
/// road (`RoadFitResult::state`): only the inward pass runs. The outward
/// pass's chi2 and kinks are not known here and stay zero.
template <CrossingRule Rule = CrossingRule::kFirstAhead, typename Prop, DetectorLayout Det, typename Material, typename Chain, typename Observer = NoFitObserver>
GPUdi() TrackFitResult fitInward(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, const BoundTrackState& outwardState, const KinkBudget& kinkBudget = {}, const Observer& observe = {})
{
  TrackFitResult fit{.ndf = static_cast<std::int16_t>((2 * n) - kNBoundParams)};
  fit.outwardState = outwardState;
  detail::inwardPass<Rule>(prop, geometry, material, measurements, n, kinkBudget, observe, fit);
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
template <bool LastStepAllowance = false, CrossingRule Rule = CrossingRule::kFirstAhead, typename Prop, DetectorLayout Det, typename Material, typename Chain>
GPUdi() RoadFitResult fitOutward(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& measurements, int n, const Scalar* chi2MaxByNdf, const KinkBudget& kinkBudget = {}, [[maybe_unused]] Scalar lastStepScale = 1.F)
{
  const auto budget = detail::kinkBudgetFor(kinkBudget, n);
  BoundTrackState state = makeChainSeed(prop, geometry, measurements, n);
  CovarianceFactors factors;
  RoadFitResult result;
  BoundTrackState reference;
  countWork(WorkCounter::kRoadStepsPossible, static_cast<std::uint64_t>(n - 1));
  // The start as `runPass` takes it, so that the kink gate, which reads the
  // increments, rounds as in the refit's outward pass.
  int taken = 0;
  std::int16_t failed = -1;
  if (detail::startPass<Rule>(prop, geometry, material, measurements, n, PropagationDirection::kForward, state, factors, result.chi2, reference, taken, failed, NoFitObserver{}) != PropagationStatus::kOk) {
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
    BoundTrackState nextReference;
    if constexpr (LastStepAllowance) {
      if (index == n - 1) {
        int taken = kinks;
        bool refused = false;
        const auto status = detail::allowanceStep<Rule>(stepProp, geometry, material, measurements[index], PropagationDirection::kForward, state, factors, result.chi2, budget, lastStepScale, taken, refused, reference, nextReference);
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
        allowed = taken - kinks;
        continue;
      }
    }
    const auto status = detail::kinkStep<Rule>(stepProp, geometry, material, measurements[index], PropagationDirection::kForward, state, factors, result.chi2, budget, index < n - 1, kinks, reference, nextReference);
    if (status != PropagationStatus::kOk) {
      countWork(WorkCounter::kRoadPropagation);
      result.status = RoadFitStatus::kPropagationFailed;
      return result;
    }
    reference = nextReference;
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
    return fitOutward<true, Rule>(prop, geometry, material, measurements, n, parameters.roadChi2MaxByNdf, kinkBudget(parameters), o2::gpu::CAMath::Sqrt(parameters.kinkScale2));
  } else {
    return fitOutward<false, Rule>(prop, geometry, material, measurements, n, parameters.roadChi2MaxByNdf, kinkBudget(parameters));
  }
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRACKFITTER_H_
