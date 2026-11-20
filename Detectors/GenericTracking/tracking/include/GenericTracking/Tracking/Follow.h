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

/// \file Follow.h
/// \brief A track followed beyond its rows, outward and then inward, as a CKF
///        follows it: the crossing and the update,
///        one step's search, the observer, the stage (`FollowStage`), the
///        passes after it, and the host's orchestration (`Follow`).
///
/// After an iteration's refit each of its tracks is followed, and, with
/// `followEarlier`, every earlier iteration's track first. The rows proposed
/// then pass the bidding between the tracks followed in their order, the
/// combination, the refit under the road fit's
/// chi2 cut, and the marks: the bid is over one iteration's tracks, and the
/// marks freeze the rows kept, so that later iterations do not build on them.

#ifndef GTRK_TRACKING_FOLLOW_H_
#define GTRK_TRACKING_FOLLOW_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/CovarianceTransport.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/KFTrackFitter.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"
#include "GenericTracking/Tracking/MeasurementUpdate.h"
#include "GenericTracking/Tracking/Profile.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/SensorCells.h"
#include "GenericTracking/Tracking/ShapeCut.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackRefit.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/WorkCounters.h"

#include "CommonConstants/MathConstants.h"
#include "GPUCommonMath.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <type_traits>

namespace o2::gtrk
{

/// A track's intersection with a layer's nominal cylinder or disc plane, in
/// the two coordinates the measurement grid indexes. The search is wider
/// than the difference between that nominal surface and an individual sensor;
/// the candidate is propagated to its actual sensor before it is chosen.
struct LayerCrossingPrediction {
  Scalar phi = 0.F;
  Scalar second = 0.F;
  PropagationStatus status = PropagationStatus::kNotFinite;
};

namespace follow
{

template <typename T>
GPUhdi() GPUAlwaysInline() T directedTurn(T turn, T curvature, PropagationDirection direction)
{
  const auto wanted = static_cast<T>(static_cast<int>(direction));
  if ((turn / curvature) * wanted < T{0}) {
    turn += wanted * (curvature > T{0} ? kTwoPi<T> : -kTwoPi<T>);
  }
  return turn;
}

/// Whether a crossing \a turn radians along the helix is past half a turn.
/// Such a crossing is a looper coming round again -- a track at eta ~ 0 reaches
/// a forward disc after about ten turns -- and not the next layer of the same
/// arc; its transported covariance is also past what float can factorise.
template <typename T>
GPUhdi() GPUAlwaysInline() bool beyondHalfTurn(T turn)
{
  return o2::gpu::CAMath::Abs(turn) > kPi<T>;
}

template <typename T>
GPUhdi() GPUAlwaysInline() bool nearerPath(T candidate, T chosen)
{
  return exactlyEqual(chosen, T{0}) || o2::gpu::CAMath::Abs(candidate) < o2::gpu::CAMath::Abs(chosen);
}

} // namespace follow

/// What every layer crossing of one state needs, taken once for the
/// state, in float: the crossing only names the layer's sensors to search
/// and orders the targets, the window and every state taken are the
/// double plane solve's. Per state: its global position, its direction's sines and cosines, its
/// curvature and, on a helix, the circle it turns on. \a status is not kOk
/// where the state runs along the beam axis.
struct LayerCrossingStart {
  Scalar position[3]{};
  Scalar start[3]{};
  Scalar sinPhi = 0.F;
  Scalar cosPhi = 0.F;
  Scalar sinTheta = 0.F;
  Scalar cosTheta = 0.F;
  Scalar curvature = 0.F;
  Scalar minAbsCosIncidence = 0.F;
  // The helix's circle in the transverse plane, where |curvature| >= 1e-12.
  Scalar inverse = 0.F;
  Scalar centreX = 0.F;
  Scalar centreY = 0.F;
  Scalar centreDistance = 0.F;
  Scalar circleRadius = 0.F;
  PropagationStatus status = PropagationStatus::kNotFinite;
};

/// \a state's `LayerCrossingStart`, \a from its surface.
template <typename Prop>
[[nodiscard]] GPUhdi() GPUAlwaysInline() LayerCrossingStart layerCrossingStart(const Prop& prop, const TrackingSurface& from, const BoundParameters& state)
{
  using Math = o2::gpu::CAMath;
  LayerCrossingStart s;
  const auto& propagator = pointPropagator(prop, from, state);
  surfaceToGlobal(from, state.loc0(), state.loc1(), s.position);
  for (int i = 0; i < 3; ++i) {
    s.start[i] = s.position[i];
  }
  sincosOf(state.phi(), s.sinPhi, s.cosPhi);
  sincosOf(state.theta(), s.sinTheta, s.cosTheta);
  if (Math::Abs(s.sinTheta) < 1.e-9F) {
    return s;
  }
  s.curvature = (state.qOverP() / s.sinTheta) * static_cast<Scalar>(propagator.bzKiloGauss) * static_cast<Scalar>(constants::math::B2C);
  s.minAbsCosIncidence = static_cast<Scalar>(propagator.minAbsCosIncidence);
  if (!(Math::Abs(s.curvature) < 1.e-12F)) {
    s.inverse = 1.F / s.curvature;
    s.centreX = s.start[0] - (s.sinPhi * s.inverse);
    s.centreY = s.start[1] + (s.cosPhi * s.inverse);
    s.centreDistance = (Math::Sqrt(((s.centreX * s.centreX) + (s.centreY * s.centreY))));
    s.circleRadius = Math::Abs(s.inverse);
  }
  s.status = PropagationStatus::kOk;
  return s;
}

/// Predict where the state of \a s crosses \a layer, without choosing a
/// sensor. The nominal radius or z is copied into IterationParameters on the
/// host, so the runtime-indexed path reads no detector constexpr table on a
/// device.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() LayerCrossingPrediction predictLayerCrossing(const IterationParameters<Det>& parameters, const LayerCrossingStart& s, LayerID layer, PropagationDirection direction)
{
  using Math = o2::gpu::CAMath;
  LayerCrossingPrediction result;
  if (s.status != PropagationStatus::kOk) {
    result.status = s.status;
    return result;
  }
  const Scalar* start = s.start;
  const Scalar sinPhi = s.sinPhi;
  const Scalar cosPhi = s.cosPhi;
  const Scalar sinTheta = s.sinTheta;
  const Scalar cosTheta = s.cosTheta;
  const Scalar curvature = s.curvature;
  const auto directionValue = static_cast<Scalar>(static_cast<int>(direction));
  Scalar path = 0.F;
  Scalar x = start[0];
  Scalar y = start[1];
  Scalar z = start[2];

  if (Detector<Det>::hasDiscs && !Detector<Det>::isBarrelLayer(layer)) {
    if (Math::Abs(cosTheta) < s.minAbsCosIncidence) {
      result.status = PropagationStatus::kGrazing;
      return result;
    }
    path = (static_cast<Scalar>(parameters.layerSpans[layer].position) - start[2]) / cosTheta;
    if (path * directionValue < 0.F) {
      result.status = PropagationStatus::kWrongDirection;
      return result;
    }
    const Scalar turn = curvature * path * sinTheta;
    if (follow::beyondHalfTurn(turn)) {
      result.status = PropagationStatus::kTooManyTurns;
      return result;
    }
    const Scalar transverse = path * sinTheta;
    Scalar sincTurn = 0.F;
    Scalar versineTurn = 0.F;
    sincVersine(turn, sincTurn, versineTurn);
    x += transverse * ((cosPhi * sincTurn) - (sinPhi * versineTurn));
    y += transverse * ((sinPhi * sincTurn) + (cosPhi * versineTurn));
    z = static_cast<Scalar>(parameters.layerSpans[layer].position);
  } else {
    const auto radius = static_cast<Scalar>(parameters.layerSpans[layer].position);
    const Scalar dot = (start[0] * cosPhi) + (start[1] * sinPhi);
    const Scalar radial = (start[0] * start[0]) + (start[1] * start[1]) - (radius * radius);
    if (Math::Abs(curvature) < 1.e-12F) {
      const Scalar discriminant = (dot * dot) - radial;
      if (discriminant < 0.F) {
        result.status = PropagationStatus::kNoIntersection;
        return result;
      }
      const auto root = (Math::Sqrt((discriminant)));
      const Scalar first = -dot - root;
      const Scalar second = -dot + root;
      Scalar transverse = 0.F;
      if (first * directionValue >= 0.F) {
        transverse = first;
      }
      if (second * directionValue >= 0.F && follow::nearerPath(second, transverse)) {
        transverse = second;
      }
      if (exactlyEqual(transverse, Scalar{0})) {
        result.status = PropagationStatus::kWrongDirection;
        return result;
      }
      path = transverse / sinTheta;
      x += transverse * cosPhi;
      y += transverse * sinPhi;
      z += path * cosTheta;
    } else {
      const Scalar centreX = s.centreX;
      const Scalar centreY = s.centreY;
      const Scalar centreDistance = s.centreDistance;
      const Scalar circleRadius = s.circleRadius;
      if (centreDistance > radius + circleRadius || centreDistance < Math::Abs(radius - circleRadius) || !(centreDistance > 0.F)) {
        result.status = PropagationStatus::kNoIntersection;
        return result;
      }
      const Scalar along = ((radius * radius) - (circleRadius * circleRadius) + (centreDistance * centreDistance)) / (2.F * centreDistance);
      const Scalar heightSquared = (radius * radius) - (along * along);
      if (heightSquared < 0.F) {
        result.status = PropagationStatus::kNoIntersection;
        return result;
      }
      const auto height = (Math::Sqrt((heightSquared)));
      const Scalar unitX = centreX / centreDistance;
      const Scalar unitY = centreY / centreDistance;
      Scalar chosenPath = 0.F;
      Scalar chosenX = 0.F;
      Scalar chosenY = 0.F;
      for (int side = -1; side <= 1; side += 2) {
        const auto sideValue = static_cast<Scalar>(side);
        const Scalar candidateX = (along * unitX) - (sideValue * height * unitY);
        const Scalar candidateY = (along * unitY) + (sideValue * height * unitX);
        const Scalar deltaX = candidateX - start[0];
        const Scalar deltaY = candidateY - start[1];
        const Scalar alongDirection = (deltaX * cosPhi) + (deltaY * sinPhi);
        const Scalar alongNormal = (-deltaX * sinPhi) + (deltaY * cosPhi);
        Scalar turn = atan2Of(curvature * alongDirection, 1.F - (curvature * alongNormal));
        turn = follow::directedTurn(turn, curvature, direction);
        const Scalar candidatePath = turn / (curvature * sinTheta);
        if (candidatePath * directionValue >= 0.F && follow::nearerPath(candidatePath, chosenPath)) {
          chosenPath = candidatePath;
          chosenX = candidateX;
          chosenY = candidateY;
        }
      }
      if (exactlyEqual(chosenPath, Scalar{0})) {
        result.status = PropagationStatus::kWrongDirection;
        return result;
      }
      // The nearer of the two circle crossings: if it is past half a turn, so
      // is the other.
      if (follow::beyondHalfTurn(curvature * chosenPath * sinTheta)) {
        result.status = PropagationStatus::kTooManyTurns;
        return result;
      }
      path = chosenPath;
      x = chosenX;
      y = chosenY;
      z += path * cosTheta;
    }
  }

  result.phi = azimuth(x, y);
  result.second = Detector<Det>::isBarrelLayer(layer) ? z : Math::Sqrt((x * x) + (y * y));
  result.status = isFiniteOf(result.phi) && isFiniteOf(result.second) ? PropagationStatus::kOk : PropagationStatus::kNotFinite;
  return result;
}

/// Predict where \a state crosses \a layer, \a from its surface: one
/// crossing from the state's `layerCrossingStart`.
template <typename Prop, DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() LayerCrossingPrediction predictLayerCrossing(const Prop& prop, const IterationParameters<Det>& parameters, const TrackingSurface& from, const BoundParameters& state, LayerID layer, PropagationDirection direction)
{
  return predictLayerCrossing(parameters, layerCrossingStart(prop, from, state), layer, direction);
}

/// Rows a track's follow would add in each direction. Outward rows are in
/// increasing road order; inward rows are in the order they are found,
/// outermost to innermost, and are reversed only when the combined RowChain
/// is read. With the end decided again (`followEnd`), \a replacedEnd is the
/// track's outermost row where the outward rows take its place; it goes only
/// with them, if outward rows are kept after the passes.
template <std::size_t K>
struct FollowProposal {
  std::uint32_t outward[K]{};
  std::uint32_t inward[K]{};
  std::uint32_t replacedEnd = std::numeric_limits<std::uint32_t>::max();
  /// A row on an inactive layer between two of the track's rows (LayerSpec.h),
  /// inserted before its row `interiorAt` (`FollowStage::interior`); none
  /// without.
  std::uint32_t interior = std::numeric_limits<std::uint32_t>::max();
  std::uint8_t interiorAt = 0U;
  std::uint8_t outwardCount = 0U;
  std::uint8_t inwardCount = 0U;
};

namespace follow
{

/// The follow's window carries the state to a layer's nominal surface with
/// it (`predictedWindow`), and `updateCandidate` is this plus the update.
/// \a traversed, where given, gets the step's material. The body, inlined
/// where it is called; the follow calls `transportCandidate`.
template <typename Solve = ExactScalar, typename Prop, DetectorLayout Det, typename Material>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool transportCandidateInline(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const TrackingSurface& to, PropagationDirection direction, const BoundParameters& GPUrestrict() input, const CovarianceFactors& GPUrestrict() inputFactors, BoundParameters& GPUrestrict() output, CovarianceFactors& GPUrestrict() outputFactors, detail::StepMaterial* GPUrestrict() traversed = nullptr)
{
  const auto& from = geometry.planeOf(input);
  const auto& propagator = pointPropagator(prop, from, input);
  Scalar start[3]{};
  surfaceToGlobal(from, input.loc0(), input.loc1(), start);
  SolenoidJacobian jacobian;
  PropagationResult propagated;
  if (jacobianAnalytic<Scalar, CrossingRule::kFirstAhead, Solve>(propagator, input, from, to, direction, jacobian, &output, &propagated) != PropagationStatus::kOk) {
    return false;
  }
  // The turn is kappa times the transverse path, (q/p) Bz B2C times the 3D
  // path: curvature() at sin(theta) = 1 against the 3D length.
  if constexpr (std::is_same_v<Solve, ExactScalar>) {
    if (beyondHalfTurn(propagator.curvature(static_cast<ExactScalar>(input.qOverP()), 1.) * static_cast<ExactScalar>(propagated.pathLength))) {
      return false;
    }
  } else {
    if (beyondHalfTurn(input.qOverP() * propagator.bzKiloGauss * static_cast<Scalar>(constants::math::B2C) * propagated.pathLength)) {
      return false;
    }
  }
  Scalar end[3]{};
  surfaceToGlobal(to, output.loc0(), output.loc1(), end);
  const auto step = detail::stepMaterial(geometry, from, to, material, start, end);
  if (traversed != nullptr) {
    *traversed = step;
  }
  const Scalar momentum = 1.F / o2::gpu::CAMath::Abs(input.qOverP());
  // As filterStep: the nearer half's noise before the transport, at the
  // direction the step leaves with.
  const auto before = scatteringNoise(highlandTheta0(step.nearerX2X0, step.nearerLogX2X0, momentum, propagator.massGeV), input.theta());
  outputFactors = inputFactors;
  return detail::transportStep<false>(propagator, step, momentum, before, direction, jacobian, output, outputFactors, 1.F) == PropagationStatus::kOk;
}

/// `transportCandidateInline`, never inlined on the device: the one copy
/// every search shares. Only
/// `FollowStage::passingOn` takes the body inline (`updateCandidate`'s
/// \a Inlined).
template <typename Solve = ExactScalar, typename Prop, DetectorLayout Det, typename Material>
[[nodiscard]] GPUhdnii() bool transportCandidate(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const TrackingSurface& to, PropagationDirection direction, const BoundParameters& GPUrestrict() input, const CovarianceFactors& GPUrestrict() inputFactors, BoundParameters& GPUrestrict() output, CovarianceFactors& GPUrestrict() outputFactors, detail::StepMaterial* GPUrestrict() traversed = nullptr)
{
  return transportCandidateInline<Solve>(prop, geometry, material, to, direction, input, inputFactors, output, outputFactors, traversed);
}

/// \a input, its covariance \a inputFactors (CovarianceFactors), carried to
/// \a measurement (`transportCandidate`) and updated with it:
/// \a output and \a outputFactors, and the update's \a chi2 -- and, with
/// \a Determinant, its innovation covariance's determinant into
/// \a determinant (the follow's EndDensity). False if the step or the
/// update is not usable. With
/// \a Inlined the transport is inlined (`transportCandidateInline`), else
/// called.
template <bool Determinant = false, bool Inlined = false, typename Solve = ExactScalar, typename Prop, DetectorLayout Det, typename Material>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool updateCandidate(const Prop& prop, const TrackingGeometry<Det>& geometry, const Material& material, const Measurement& measurement, PropagationDirection direction, const BoundParameters& GPUrestrict() input, const CovarianceFactors& GPUrestrict() inputFactors, BoundParameters& GPUrestrict() output, CovarianceFactors& GPUrestrict() outputFactors, Scalar& GPUrestrict() chi2, [[maybe_unused]] Scalar* GPUrestrict() determinant = nullptr)
{
  bool carried = false;
  if constexpr (Inlined) {
    carried = transportCandidateInline<Solve>(prop, geometry, material, geometry.planeOf(measurement), direction, input, inputFactors, output, outputFactors);
  } else {
    carried = transportCandidate<Solve>(prop, geometry, material, geometry.planeOf(measurement), direction, input, inputFactors, output, outputFactors);
  }
  if (!carried) {
    return false;
  }
  const auto update = updateWithMeasurement<Determinant>(output, outputFactors, measurement, determinant);
  chi2 = update.chi2;
  // A negative chi2 would otherwise win the best-candidate comparison: the
  // same test the fit passes make (`updateUsable`).
  return updateUsable(update, outputFactors);
}

/// The half-widths of the planes `predictedWindow` builds: no bound [cm].
inline constexpr Scalar kOpenBound = 1.e4F;

/// The window of one step of a follow, narrowing \a window about the nominal crossing \a nominal. The
/// state is
/// carried to the plane touching the layer there, with the material's noise
/// (`transportCandidate`), and the window's half-widths are:
///  - the gate's sigma of the prediction's spread, its covariance plus the
///    pair windows' resolution term, projected onto phi and the second
///    coordinate;
///  - plus the layer's envelope about its nominal position, its
///    measurements' radii on a barrel layer or z on a disc, times the track's
///    slope there: a sensor anywhere in it, tilted, overlapping or curved, is
///    crossed within that;
///  - at most \a window's as given.
/// Where the transport fails, \a window stands, and it returns false.
/// Which candidate is taken is still the update's chi2 on its own sensor.
/// \a reach, where given, gets the envelope's term alone, azimuth then the
/// second coordinate: how far from the nominal crossing the track itself
/// may cross a sensor of the layer (the follow's hole test). \a traversed,
/// where given, gets the material the state crossed to the plane: what the
/// window's candidates take (`detail::GivenMaterial`).
template <typename Prop, DetectorLayout Det, typename Material>
GPUhdi() bool predictedWindow(const Prop& prop, const IterationParameters<Det>& parameters, const TrackingGeometry<Det>& geometry, const Material& material, const BoundParameters& state, const CovarianceFactors& factors, LayerID target, PropagationDirection direction, const LayerCrossingPrediction& nominal, CandidateWindow& window, Scalar* reach = nullptr, detail::StepMaterial* traversed = nullptr)
{
  using Math = o2::gpu::CAMath;
  const bool barrel = ((parameters.barrelLayers >> target) & 1U) != 0U;
  // The plane touching the layer at the nominal crossing: on a barrel layer
  // loc0 along r phi and loc1 along z, on a disc loc0 along r and loc1
  // along r phi.
  const Scalar position = parameters.layerSpans[target].position;
  const Scalar radius = barrel ? position : nominal.second;
  if (!(radius > 0.F)) {
    return false;
  }
  Scalar sinPhi = 0.F;
  Scalar cosPhi = 0.F;
  sincosOf(nominal.phi, sinPhi, cosPhi);
  const Scalar origin[3] = {radius * cosPhi, radius * sinPhi, barrel ? nominal.second : position};
  TrackingSurface plane = barrel ? radialPlane(origin, sinPhi, cosPhi) : discPlane(origin, sinPhi, cosPhi);
  plane.bounds = {.halfX = kOpenBound, .halfY = kOpenBound};
  plane.layer = target;
  BoundParameters predicted;
  CovarianceFactors predictedFactors;
  // Navigation: the window's plane, centre and widths need no more than a
  // float solve (the states taken into the track are transported in double).
  if (!transportCandidate<Scalar>(prop, geometry, material, plane, direction, state, factors, predicted, predictedFactors, traversed)) {
    return false;
  }
  Scalar predictedCov[kBoundCovSize]{};
  covarianceFromFactors(predictedFactors, predictedCov);
  const Scalar variance0 = predictedCov[varIndex(kLoc0)] + parameters.followPositionVariance;
  const Scalar variance1 = predictedCov[varIndex(kLoc1)] + parameters.followPositionVariance;
  if (!(variance0 > 0.F) || !(variance1 > 0.F)) {
    return false;
  }
  // The window spans what the update's gate admits on one axis,
  // sqrt(clusterChi2Max), not the pair cut's nSigma: a row the gate takes is
  // never outside the window.
  const Scalar nSigma = Math::Sqrt(parameters.clusterChi2Max);
  const Scalar phiSpread = nSigma * Math::Sqrt(barrel ? variance0 : variance1) / radius;
  const Scalar secondSpread = nSigma * Math::Sqrt(barrel ? variance1 : variance0);
  // How far the crossing moves per unit of the layer's position, along the
  // track there: alpha is its transverse direction from the radial one. The
  // transverse path is 1 / cos(alpha) per unit of radius on a barrel layer,
  // tan(theta) per unit of z on a disc.
  Scalar sinAlpha = 0.F;
  Scalar cosAlpha = 0.F;
  sincosOf(predicted.phi() - nominal.phi, sinAlpha, cosAlpha);
  Scalar sinTheta = 0.F;
  Scalar cosTheta = 0.F;
  sincosOf(predicted.theta(), sinTheta, cosTheta);
  const Scalar transverse = barrel ? 1.F / Math::Abs(cosAlpha) : Math::Abs(sinTheta / cosTheta);
  const Scalar phiSlope = transverse * Math::Abs(sinAlpha) / radius;
  const Scalar secondSlope = transverse * (barrel ? Math::Abs(cosTheta / sinTheta) : Math::Abs(cosAlpha));
  const auto& span = parameters.layerSpans[target];
  const Scalar envelope = Math::Max(span.positionMax - position, position - span.positionMin);
  // The window given first: a NaN or infinite slope leaves it.
  window.phiHalfWidth = Math::Min(window.phiHalfWidth, phiSpread + (envelope * phiSlope));
  window.secondHalfWidth = Math::Min(window.secondHalfWidth, secondSpread + (envelope * secondSlope));
  if (reach != nullptr) {
    reach[0] = envelope * phiSlope;
    reach[1] = envelope * secondSlope;
  }
  return true;
}

/// The material a follow step from \a state (its covariance \a factors)
/// to \a layer in \a direction gives every row it updates there: its
/// window's (`predictedWindow`, the nominal step), into \a traversed. False
/// where the step makes no window. For the studies that replay a step
/// without its window.
template <typename Prop, DetectorLayout Det, typename Material>
[[nodiscard]] bool windowMaterial(const Prop& prop, const IterationParameters<Det>& parameters, const TrackingGeometry<Det>& geometry, const Material& material, const BoundParameters& state, const CovarianceFactors& factors, LayerID layer, PropagationDirection direction, detail::StepMaterial& traversed)
{
  const auto nominal = predictLayerCrossing(prop, parameters, geometry.planeOf(state), state, layer, direction);
  if (nominal.status != PropagationStatus::kOk) {
    return false;
  }
  CandidateWindow window{.phiHalfWidth = kPi<Scalar>, .secondHalfWidth = kOpenBound};
  return predictedWindow(prop, parameters, geometry, material, state, factors, layer, direction, nominal, window, nullptr, &traversed);
}

/// The nominal crossing \a prediction of \a layer as a point [cm], into
/// \a point: on a barrel layer its nominal radius at the crossing's azimuth
/// and z, on a disc its radius and azimuth at the disc's z.
template <DetectorLayout Det>
GPUhdi() GPUAlwaysInline() void crossingPoint(const IterationParameters<Det>& parameters, LayerID layer, const LayerCrossingPrediction& prediction, Scalar* point)
{
  const bool barrel = ((parameters.barrelLayers >> layer) & 1U) != 0U;
  const Scalar position = parameters.layerSpans[layer].position;
  const Scalar r = barrel ? position : prediction.second;
  Scalar sinPhi = 0.F;
  Scalar cosPhi = 0.F;
  sincosOf(prediction.phi, sinPhi, cosPhi);
  point[0] = r * cosPhi;
  point[1] = r * sinPhi;
  point[2] = barrel ? prediction.second : position;
}

/// The centre of a step's window on \a layer, in the rows' coordinates: the
/// nominal crossing \a prediction, which is the detector frame's, as its
/// point's azimuth and second coordinate about the beam line, as the
/// measurements have theirs (`gridPoint`).
template <DetectorLayout Det>
GPUhdi() GPUAlwaysInline() void windowCentre(const IterationParameters<Det>& parameters, LayerID layer, const LayerCrossingPrediction& prediction, Scalar& phi, Scalar& second)
{
  Scalar point[3]{};
  crossingPoint(parameters, layer, prediction, point);
  const Scalar x = point[0] - parameters.beamX;
  const Scalar y = point[1] - parameters.beamY;
  phi = azimuth(x, y);
  second = secondCoordinate<Det>(layer, o2::gpu::CAMath::Hypot(x, y), point[2]);
}

/// What a follow's observer sees (`FollowStage::follow`); the default sees
/// nothing and compiles to nothing:
///  - each direction's start, and each step: the layer it leaves, its
///    parameters, its covariance's factors and the rows the direction has
///    taken;
///  - each target the step searches, with its window;
///  - each row the search visits, with its outcome (`visit`) and the update's
///    chi2 (0 where the update was not made);
///  - a ranking by confirmation (the lookahead), with the rows passing;
///  - each hole test, with whether the track crosses a sensor;
///  - the row a step takes and its chi2;
///  - the last row dropped as nothing confirms it (the unconfirmed drop), refused at the
///    layout's end (EndDensity) or against the beam line (BeamEnd), and the
///    end.
struct NoObserver {
  GPUhdi() GPUAlwaysInline() void begin(PropagationDirection /*direction*/) const {}
  GPUhdi() GPUAlwaysInline() void step(LayerID /*layer*/, const BoundParameters& /*state*/, const CovarianceFactors& /*factors*/, std::uint32_t /*rows*/) const {}
  GPUhdi() GPUAlwaysInline() void target(LayerID /*layer*/, const CandidateWindow& /*window*/) const {}
  GPUhdi() GPUAlwaysInline() void candidate(LayerID /*layer*/, std::uint32_t /*row*/, StepOutcome /*outcome*/, Scalar /*chi2*/) const {}
  GPUhdi() GPUAlwaysInline() void lookahead(std::uint32_t /*passing*/) const {}
  GPUhdi() GPUAlwaysInline() void hole(LayerID /*layer*/, bool /*crossed*/) const {}
  /// the lookahead's confirmation of a nominal row (`confirm`): the row, its chi2, and
  /// the next step's row confirming it and that row's chi2 (kNoCluster: none).
  GPUhdi() GPUAlwaysInline() void lookaheadConfirmation(std::uint32_t /*row*/, Scalar /*chi2*/, std::uint32_t /*next*/, Scalar /*nextChi2*/) const {}
  GPUhdi() GPUAlwaysInline() void taken(std::uint32_t /*row*/, Scalar /*chi2*/) const {}
  GPUhdi() GPUAlwaysInline() void unconfirmed(std::uint32_t /*row*/) const {}
  /// The last row refused at the layout's end (EndDensity's layout-end refusal).
  GPUhdi() GPUAlwaysInline() void refusedEnd(std::uint32_t /*row*/) const {}
  /// The inward last row refused against the beam line (BeamEnd).
  GPUhdi() GPUAlwaysInline() void refusedBeamEnd(std::uint32_t /*row*/) const {}
  GPUhdi() GPUAlwaysInline() void end(FollowEnd /*why*/) const {}
};

/// What a step's search reads: the iteration, the frame and its grid, the
/// sensors, the material, the field and the used mask.
template <DetectorLayout Det, typename Material>
struct StepContext {
  const IterationParameters<Det>* parameters = nullptr;
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  const Material* material = nullptr;
  const Propagator* propagator = nullptr;
  const std::uint8_t* used = nullptr;
};

/// A step's row: the row, its update's chi2, and the state updated with it.
struct Child {
  std::uint32_t row = std::numeric_limits<std::uint32_t>::max();
  Scalar chi2 = 0.F;
  BoundParameters state;
  CovarianceFactors factors;
};

/// A step's row as EndDensity reads it: with its update's innovation covariance's determinant
/// [cm^4]. A type of its own, so that the follow without it keeps `Child`.
struct DensityChild : Child {
  Scalar determinant = 0.F;
};

/// The follow's row type with or without \a EndDensity.
template <bool EndDensity>
using ChildOf = std::conditional_t<EndDensity, DensityChild, Child>;

/// Where \a child's update writes its determinant (`updateCandidate`):
/// nowhere for a `Child`.
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar* determinantOf(Child& /*child*/) { return nullptr; }
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar* determinantOf(DensityChild& child) { return &child.determinant; }

/// One search of a step: each row of \a layer the grid holds in \a window, in
/// the grid's order, through the step's tests in the order `StepOutcome`
/// names them: free (with \a Own, \a own counts as free: the track's own end
/// decided again), inside the window's spatial bounds, its segment from
/// (\a r, \a z) -- where the track stands -- crossing none of \a holes (0: no
/// test, LinkHoles.h), and the update of \a state (its covariance \a factors)
/// with it under the gate (`updateCandidate`), every row's step taking
/// \a traversed, the window's material (`predictedWindow`). With
/// \a EndDensity a candidate holds its update's determinant.
/// \a report(row, outcome, chi2) sees every row; \a counter counts the
/// updates made.
///
/// The lanes of \a group walk the window in step
/// (`MeasurementGrid::visitCandidates`), each its share of the rows, and
/// every lane calls \a step(candidate) at every step of the walk:
/// \a candidate its row passing the gate, with the state updated with it and
/// its factors, or none there (its row `kNoCluster`). So \a step may call
/// the group's collectives. One lane has nothing to keep in step: it calls
/// \a step for each row passing the gate alone, right after its report, in
/// the grid's order. \a Inlined: the transport of each row inlined
/// (`updateCandidate`). With \a Shape, the iteration's `shapeTest`, a usable
/// update is refused where the row's slope interval does not hold the
/// updated state's direction (`slopeAdmits`, StepOutcome::kShape), before
/// the gate. The body, inlined where it is called; the follow calls `visit`.
template <bool Own = false, bool EndDensity = false, bool Inlined = false, bool Shape = false, DetectorLayout Det, typename Material, typename Group, typename Report, typename Step>
GPUhdi() GPUAlwaysInline() void visitInline(const StepContext<Det, Material>& context, const BoundParameters& state, const CovarianceFactors& factors, Scalar r, Scalar z, LayerID layer, PropagationDirection direction, const CandidateWindow& window, const detail::StepMaterial& traversed, std::uint64_t holes, WorkCounter counter, const Group& group, const Report& report, const Step& step, [[maybe_unused]] std::uint32_t own = std::numeric_limits<std::uint32_t>::max())
{
  const auto& parameters = *context.parameters;
  const auto& frame = *context.frame;
  const FollowGridScope followCounts;
  // The tests before the update, each reported where it refuses the row.
  const auto reaches = [&](std::uint32_t row) {
    if (context.used[row] != kFree && (!Own || row != own)) {
      report(row, StepOutcome::kUsed, 0.F);
      return false;
    }
    const auto measurement = frame.measurement(row);
    if (!spatiallyCompatible(window, measurement.phi, secondCoordinate<Det>(layer, measurement.r, measurement.z))) {
      report(row, StepOutcome::kOutsideWindow, 0.F);
      return false;
    }
    if (holes != 0U && segmentLeavesHole(holes, parameters.barrelLayers, parameters.layerSpans, r, z, measurement.r, measurement.z)) {
      report(row, StepOutcome::kHole, 0.F);
      return false;
    }
    return true;
  };
  // The update into \a candidate, its row set where it passes the gate.
  const auto update = [&](std::uint32_t row, ChildOf<EndDensity>& candidate) {
    countWork(counter);
    Scalar chi2 = 0.F;
    // The follow's states steer its own steps only (a proposal hands on
    // rows; the refit fits them again): a float solve.
    const bool usable = updateCandidate<EndDensity, Inlined, Scalar>(*context.propagator, *context.geometry, detail::GivenMaterial{traversed}, measurementAt(frame, row), direction, state, factors, candidate.state, candidate.factors, chi2, determinantOf(candidate));
    if (!usable) {
      report(row, StepOutcome::kNoUpdate, 0.F);
      return;
    }
    if constexpr (Shape) {
      const auto measurement = frame.measurement(row);
      const Scalar theta = candidate.state.theta();
      const Scalar phi = candidate.state.phi();
      const Scalar sinTheta = o2::gpu::CAMath::Sin(theta);
      if (!slopeAdmits(measurement.slopeLow, measurement.slopeHigh, measurement.x, measurement.y, measurement.r, o2::gpu::CAMath::Cos(phi) * sinTheta, o2::gpu::CAMath::Sin(phi) * sinTheta, o2::gpu::CAMath::Cos(theta))) {
        report(row, StepOutcome::kShape, chi2);
        return;
      }
    }
    const bool passes = chi2 <= parameters.clusterChi2Max;
    report(row, passes ? StepOutcome::kPasses : StepOutcome::kChi2, chi2);
    if (passes) {
      candidate.row = row;
      candidate.chi2 = chi2;
    }
  };
  if constexpr (Group::kLanes == 1U) {
    // One lane: the grid's walk of one thread, the candidate of the row
    // alone, made where it is updated.
    frame.grid().visitCandidates(layer, window, frame.grid().geometry().allBlocks(), [&](std::uint32_t row) {
      if (!reaches(row)) {
        return;
      }
      ChildOf<EndDensity> candidate;
      update(row, candidate);
      if (candidate.row != std::numeric_limits<std::uint32_t>::max()) {
        step(candidate);
      }
    });
  } else {
    // A group: one candidate a lane for the window, its row reset at every
    // step; the state and its factors are written by the update before they
    // are read.
    ChildOf<EndDensity> candidate;
    frame.grid().visitCandidates(
      layer, window, frame.grid().geometry().allBlocks(),
      [&](std::uint32_t row, bool has) {
        candidate.row = std::numeric_limits<std::uint32_t>::max();
        if (has && reaches(row)) {
          update(row, candidate);
        }
        step(candidate);
      },
      group);
  }
}

/// `visitInline`, never inlined on the device: one copy for each call site,
/// which the kernels then do not absorb.
template <bool Own = false, bool EndDensity = false, bool Shape = false, DetectorLayout Det, typename Material, typename Group, typename Report, typename Step>
GPUhdnii() void visit(const StepContext<Det, Material>& GPUrestrict() context, const BoundParameters& GPUrestrict() state, const CovarianceFactors& GPUrestrict() factors, Scalar r, Scalar z, LayerID layer, PropagationDirection direction, const CandidateWindow& GPUrestrict() window, const detail::StepMaterial& GPUrestrict() traversed, std::uint64_t holes, WorkCounter counter, const Group& group, const Report& report, const Step& step, [[maybe_unused]] std::uint32_t own = std::numeric_limits<std::uint32_t>::max())
{
  visitInline<Own, EndDensity, false, Shape>(context, state, factors, r, z, layer, direction, window, traversed, holes, counter, group, report, step, own);
}

/// \a candidate of lane \a from of \a group, on every lane: its row, chi2,
/// parameters, surface and factors (and with \a EndDensity its determinant).
template <typename Group, typename Row>
[[nodiscard]] GPUhdi() GPUAlwaysInline() decltype(auto) shared(const Group& group, const Row& candidate, std::uint32_t from)
{
  if constexpr (Group::kLanes == 1U) {
    // The candidate itself: one lane has it already.
    static_cast<void>(group);
    static_cast<void>(from);
    return (candidate);
  } else {
    struct Wire {
      std::uint32_t row;
      Scalar chi2;
      Scalar par[kNBoundParams];
      SensorID surface;
      Scalar surfacePhi;
      CovarianceFactors factors;
    };
    Wire wire{
      .row = candidate.row,
      .chi2 = candidate.chi2,
      .par = {},
      .surface = candidate.state.surface,
      .surfacePhi = candidate.state.surfacePhi,
      .factors = candidate.factors,
    };
    for (int i = 0; i < kNBoundParams; ++i) {
      wire.par[i] = candidate.state.par[i];
    }
    wire = group.broadcast(wire, from);
    Row out;
    if constexpr (std::is_same_v<Row, DensityChild>) {
      out.determinant = group.broadcast(candidate.determinant, from);
    }
    out.row = wire.row;
    out.chi2 = wire.chi2;
    for (int i = 0; i < kNBoundParams; ++i) {
      out.state.par[i] = wire.par[i];
    }
    out.state.surface = wire.surface;
    out.state.surfacePhi = wire.surfacePhi;
    out.factors = wire.factors;
    return out;
  }
}

/// What `FollowStage::passingOn` sees on a layer: the rows passing the nominal gate, whether the row asked
/// about is among them and its chi2, and the smallest chi2 of the others.
struct Passing {
  std::uint32_t count = 0U;
  bool own = false;
  Scalar ownChi2 = 0.F;
  Scalar otherChi2 = std::numeric_limits<Scalar>::max();
};

/// The witness of a confirmation: the row the next step takes, by the
/// nominal pass's choice -- the smallest chi2, the row breaking ties -- or
/// none (`kNoCluster`).
struct Witness {
  std::uint32_t row = std::numeric_limits<std::uint32_t>::max();
  Scalar chi2 = 0.F;
};

/// Whether \a row with \a chi2 beats \a best: the nominal pass's choice, as
/// the one-lane walk makes it row by row. Over the usable, gated chi2 (no
/// NaN) a strict total order on the rows, so the choice does not depend on
/// the order the rows come in.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool beats(std::uint32_t row, Scalar chi2, const Witness& best) { return best.row == std::numeric_limits<std::uint32_t>::max() || chi2 < best.chi2 || (!(best.chi2 < chi2) && row < best.row); }

/// The best of the lanes' \a own witnesses (`beats`), on every lane of
/// \a group.
template <typename Group>
[[nodiscard]] GPUhdi() GPUAlwaysInline() Witness bestOf(const Group& group, const Witness& own)
{
  if constexpr (Group::kLanes == 1U) {
    return own;
  } else {
    Witness best;
    for (auto lanes = group.ballot(own.row != std::numeric_limits<std::uint32_t>::max()); lanes != 0U; lanes &= lanes - 1U) {
      const auto other = group.broadcast(own, o2::gpu::CAMath::Ctz64(lanes));
      if (beats(other.row, other.chi2, best)) {
        best = other;
      }
    }
    return best;
  }
}

/// The most links leaving or entering one layer of \a Det: the targets a
/// step orders.
template <DetectorLayout Det>
inline constexpr std::size_t maxDegree = [] {
  std::size_t most = 1U;
  for (std::size_t layer = 0; layer < static_cast<std::size_t>(Detector<Det>::nLayers); ++layer) {
    most = std::max(most, outLinks<Det>.count(static_cast<LayerID>(layer)));
    most = std::max(most, inLinks<Det>.count(static_cast<LayerID>(layer)));
  }
  return most;
}();

/// One direction of a track's follow: where it stands, its covariance's
/// factors, the rows it took, the last one's chi2 (the unconfirmed drop removes it), the holes it
/// crossed, whether its latest step, one that took no row, crossed a
/// sensor without a row, and the track's own row its searches take as free
/// (the end decided again, `followEnd`; none otherwise).
template <std::size_t K>
struct Branch {
  BoundParameters state;
  CovarianceFactors factors;
  std::uint32_t rows[K]{};
  Scalar lastChi2 = 0.F;
  std::uint32_t count = 0U;
  std::uint32_t holes = 0U;
  std::uint32_t own = std::numeric_limits<std::uint32_t>::max();
  bool crossed = false;
};

/// A direction as EndDensity reads it: with its last row's update's
/// determinant [cm^4]. A type of its own, so that the follow without it
/// keeps `Branch`.
template <std::size_t K>
struct DensityBranch : Branch<K> {
  Scalar lastDeterminant = 0.F;
};

/// The follow's direction type with or without \a EndDensity.
template <std::size_t K, bool EndDensity>
using BranchOf = std::conditional_t<EndDensity, DensityBranch<K>, Branch<K>>;

/// What BeamEnd reads of an inward walk: the states its last two steps started from -- entry 0 the
/// latest -- their factors, and how many rows the walk took. A row the unconfirmed drop or the
/// layout's end removes makes the one before the last: entry 1.
struct BeamEndHistory {
  BoundParameters state[2];
  CovarianceFactors factors[2];
  std::uint32_t taken = 0U;

  /// A step from \a from (its factors \a fromFactors) taking a row.
  GPUhdi() GPUAlwaysInline() void remember(const BoundParameters& from, const CovarianceFactors& fromFactors)
  {
    state[1] = state[0];
    factors[1] = factors[0];
    state[0] = from;
    factors[0] = fromFactors;
    ++taken;
  }
};

/// The follow without BeamEnd keeps nothing.
struct NoBeamEndHistory {
};

/// The inward walk's history with or without \a BeamEnd.
template <bool BeamEnd>
using BeamEndHistoryOf = std::conditional_t<BeamEnd, BeamEndHistory, NoBeamEndHistory>;

/// \a rows with \a row inserted before their \a at-th, as the fit reads a
/// chain: a track with its interior row (`FollowStage::interior`).
template <DetectorLayout Det>
struct InsertedChain {
  const TimeFrame<Det>* frame = nullptr;
  const std::uint32_t* rows = nullptr;
  std::uint32_t at = 0U;
  std::uint32_t row = 0U;

  [[nodiscard]] GPUhdi() GPUAlwaysInline() Measurement operator[](int i) const
  {
    const auto k = static_cast<std::uint32_t>(i);
    return measurementAt(*frame, k < at ? rows[k] : (k == at ? row : rows[k - 1U]));
  }
};

/// A layer a step may reach: its nominal crossing and its distance from the
/// track, squared.
struct Target {
  Scalar distance2 = 0.F;
  LayerID layer = 0;
  LayerCrossingPrediction prediction;
};

} // namespace follow

/// One track's follow: outward from its outward
/// state, then inward from its inner state, into its proposal. It reads the
/// grid and the used mask, and writes only its own proposal.
///
/// A step, in either direction:
///  - **Navigation.** The layers every link of the layout reaches from the
///    one the track stands on (enters it, inward), skip links included, in
///    the order the track crosses them: by the distance of the nominal
///    crossing from the track. The first with a free row within the gate
///    gives the step's row. A layer the state cannot be carried to is
///    passed.
///  - **The window.** The prediction's own (`predictedWindow`) and the
///    track's time box, which every row taken narrows, outward then inward.
///  - **Holes.** A layer without a row is a hole only where the track crosses
///    one of its sensors within its bounds (SensorCells.h); a gap costs
///    nothing. Past the iteration's `followHoles` the direction ends; each
///    direction has its own count.
///  - **The lookahead.** Where more than one row passes, each is taken and the next step
///    made from there: the rows it confirms rank first, by their chi2 plus
///    the confirming row's; the others by their chi2; then the row.
///  - **The unconfirmed drop.** A direction whose last step crossed a sensor without a row
///    drops its last row, which nothing confirms -- unless its chi2 is below
///    the iteration's `followKeepChi2` (0: never), or, inward, unless that step
///    saw a row near: one within its window the update refused above the
///    gate (inward the sensor usually held the particle's row, refused;
///    outward the row near is a noise row too often). The row seen is not
///    taken.
///
/// Inward starts from the inner state updated with the beam line where the
/// iteration has the beam constraint (\a BeamLine, its own instantiation),
/// and does not run where the iteration shares innermost rows (\a inward).
/// The directions share one capacity, the longest chain less the track's
/// rows.
///
/// With \a DecidesEnd, the iteration's `followEnd`, the track's outermost row is decided again (`decideEnd`): from
/// the road fit's outward pass without it (`follow::EndFit`, a pass of its
/// own before the stage, \a endStarts), the walk outward, that row free to
/// it. Its first row the end: the walk's rows beyond it are the
/// outward follow. Another: the walk's rows replace the end. None: the end
/// stays and the follow is as without. An iteration without `follow` walks
/// one row, the end alone, and nothing inward.
///
/// With \a Greedy, the iteration's `followGreedy`, a step more than one row passes takes the smallest chi2, the
/// row breaking ties (`follow::beats`), instead of the lookahead's confirmation.
///
/// With \a EndDensity, where the iteration has `followEndMaxSqrtDet` or
/// `followKeepDensity`, a
/// direction's end is decided on its last row's own update, V = H P H^T + R
/// and the chi2:
///  - **inward at the layout's end** -- the layer it stands on has no
///    inward link -- the last row goes where sqrt|V| is above
///    `followEndMaxSqrtDet`: a prediction carried that far admits strangers;
///  - **outward, where the unconfirmed drop would remove it**, the last row stays where its
///    Gaussian density g = exp(-chi2 / 2) / (2 pi sqrt|V|) is above
///    `followKeepDensity`: the keep's chi2 with a penalty for a wide
///    prediction.
///
/// With \a BeamEnd, the iteration's `followBeamEnd` (only with \a BeamLine), the inward walk's last
/// row is tested against that step's
/// prediction updated with the beam line (`refusesBeamEnd`): refused above
/// the gate. Its own step saw only the rows outside; the beam line is the
/// evidence beyond the last row inward, where the strangers of the vertex
/// detector -- QED pairs in Pb-Pb -- pass the step alone.
template <DetectorLayout Det, typename Material, bool BeamLine = false, bool DecidesEnd = false, bool Greedy = false, bool EndDensity = false, bool BeamEnd = false, bool Shape = false>
class FollowStage
{
 public:
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  static_assert(BeamLine || !BeamEnd, "the beam-line end needs the iteration's beam constraint");
  /// Whether it decides the track's end again (`DecidesEnd`).
  static constexpr bool kDecidesEnd = DecidesEnd;
  using TrackView = typename BoundedColumns<TrackOf<kMaxClusters>::template RecordT>::ConstView;
  using Proposal = FollowProposal<kMaxClusters>;
  using Branch = follow::BranchOf<kMaxClusters, EndDensity>;
  using Child = follow::ChildOf<EndDensity>;

  FollowStage(const IterationParameters<Det>* parameters, const TimeFrame<Det>* frame, const TrackingGeometry<Det>* geometry, const Material* material, const Propagator& propagator, TrackView tracks, const std::uint8_t* used, const std::uint32_t* maxCandidateWidthPS, SensorCellsView cells, Proposal* proposals, bool inward, std::uint32_t trackOffset, const BoundTrackState* endStarts = nullptr)
    : mParameters(parameters),
      mFrame(frame),
      mGeometry(geometry),
      mMaterial(material),
      mPropagator(propagator),
      mTracks(tracks),
      mUsed(used),
      mMaxCandidateWidthPS(maxCandidateWidthPS),
      mCells(cells),
      mProposals(proposals),
      mInward(inward),
      mTrackOffset(trackOffset),
      mEndStarts(endStarts)
  {
  }

  GPUhdi() void operator()(std::uint32_t track) const { together(track, SoloGroup{}); }

  /// `operator()` for the lanes of \a group, which follow \a track together:
  /// every lane makes the same decisions, the
  /// first writes the proposal.
  template <typename Group>
  GPUhdi() void together(std::uint32_t track, const Group& group) const
  {
    const std::uint32_t input = mTrackOffset + track;
    Proposal proposal;
    for (std::size_t i = 0; i < kMaxClusters; ++i) {
      proposal.outward[i] = ClusterList<kMaxClusters>::kNoCluster;
      proposal.inward[i] = ClusterList<kMaxClusters>::kNoCluster;
    }
    if (mTracks.status[input] == PropagationStatus::kOk) {
      follow(mTracks.state[input], mTracks.outwardState[input], mTracks.clusters[input], mTracks.count[input], proposal, follow::NoObserver{}, mEndStarts != nullptr ? mEndStarts + track : nullptr, group);
    }
    if (group.lane() == 0U) {
      mProposals[track] = proposal;
    }
  }

  /// A track's follow from \a inner and \a outer, its inner and outward
  /// states, its \a count rows \a rows giving the time box, into
  /// \a proposal's rows (empty before),
  /// \a observe seeing every decision: what `operator()` does for a track,
  /// for the bench's replay too. With \a DecidesEnd, \a endStart is the
  /// track's state without its outermost row (`follow::EndFit`). The lanes of
  /// \a group follow it together.
  template <typename Observer, typename Group = SoloGroup>
  GPUhdi() void follow(const BoundTrackState& inner, const BoundTrackState& outer, const ClusterList<kMaxClusters>& rows, std::uint32_t count, Proposal& proposal, const Observer& observe, [[maybe_unused]] const BoundTrackState* endStart = nullptr, const Group& group = {}) const
  {
    // A track as long as a chain can be is not followed, nor one below the
    // follow's pT minimum: its windows, sized by the scattering at its
    // momentum, would take in whole layers.
    if (count >= kMaxClusters || !reachesPt(inner, mParameters->followMinPt)) {
      return;
    }
    if constexpr (DecidesEnd) {
      if (endStart != nullptr && endStart->surface != kInvalidSensorID && decideEnd(inner, *endStart, rows, count, proposal, observe, group)) {
        return;
      }
      if (!mParameters->follow) {
        return;
      }
    }
    std::uint32_t capacity = static_cast<std::uint32_t>(kMaxClusters) - count;
    auto box = trackTimeBox(mFrame->measurements(), rows, static_cast<std::uint8_t>(count));
    // A row on an inactive layer the track's road straddles (`interior`).
    if constexpr (Detector<Det>::inactiveLayers != 0U) {
      if (capacity > 1U && interior(outer, rows, count, box, proposal)) {
        --capacity;
      }
    }
    const std::uint64_t heldLayers = heldInactive(rows, count, proposal);
    proposal.outwardCount = walk(outer, PropagationDirection::kForward, proposal.outward, capacity, box, heldLayers, observe, std::numeric_limits<std::uint32_t>::max(), group);
    walkInward(inner, rows.index[0], capacity - proposal.outwardCount, box, heldLayers, proposal, observe, group);
  }

  /// The free rows of \a layer passing the nominal gate from \a state, its
  /// covariance \a factors, in \a direction within the time \a box: what a
  /// nominal step from there would see on that layer, \a own among them or
  /// not (`follow::TrimAmbiguous`). None where no link of the layout leads from
  /// the state's layer to \a layer, or the window cannot be made.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() follow::Passing passingOn(const BoundParameters& state, const CovarianceFactors& factors, LayerID layer, std::uint32_t own, const TimeEstPS& box, PropagationDirection direction) const
  {
    follow::Target targets[follow::maxDegree<Det>];
    const std::uint32_t nTargets = targetsOf(state, direction, targets);
    for (std::uint32_t t = 0; t < nTargets; ++t) {
      const auto& target = targets[t];
      if (target.layer != layer) {
        continue;
      }
      CandidateWindow window;
      Scalar reach[2]{};
      detail::StepMaterial traversed;
      if (!targetWindow(state, factors, box, target, direction, window, reach, traversed, SoloGroup{})) {
        return {};
      }
      follow::Passing passing;
      // The search inlined with its transports: called, it slows
      // `TrimAmbiguous`'s kernel on a device.
      follow::visitInline<false, false, true, Shape>(
        context(), state, factors, 0.F, 0.F, target.layer, direction, window, traversed, 0U, WorkCounter::kFollowTrimUpdates, SoloGroup{}, [](std::uint32_t /*row*/, StepOutcome /*outcome*/, Scalar /*chi2*/) {},
        [&](const follow::Child& candidate) {
          ++passing.count;
          if (candidate.row == own) {
            passing.own = true;
            passing.ownChi2 = candidate.chi2;
          } else {
            passing.otherChi2 = o2::gpu::CAMath::Min(passing.otherChi2, candidate.chi2);
          }
        });
      return passing;
    }
    return {};
  }

  /// The walk inward after the outward one, from \a inner (`inwardStart`)
  /// over \a capacity rows at most, where the stage follows inward and rows
  /// are left; \a heldLayers as `walk`'s.
  template <typename Observer, typename Group>
  GPUhdi() GPUAlwaysInline() void walkInward(const BoundTrackState& inner, std::uint32_t innermostRow, std::uint32_t capacity, TimeEstPS& box, std::uint64_t heldLayers, Proposal& proposal, const Observer& observe, const Group& group) const
  {
    if (mInward && capacity != 0U) {
      proposal.inwardCount = walk(inwardStart(inner, innermostRow), PropagationDirection::kBackward, proposal.inward, capacity, box, heldLayers, observe, std::numeric_limits<std::uint32_t>::max(), group);
    }
  }

  /// The track's outermost row decided again (\a DecidesEnd): from \a start,
  /// \a rows' first \a count - 1 fitted outward (`follow::EndFit`), the
  /// walk outward with the end free to it -- the follow's capacity, or one
  /// row without `follow`. Whether the walk took a row: then \a proposal
  /// holds the follow, inward too with `follow`; else it is untouched.
  template <typename Observer, typename Group>
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool decideEnd(const BoundTrackState& inner, const BoundTrackState& start, const ClusterList<kMaxClusters>& rows, std::uint32_t count, Proposal& proposal, const Observer& observe, const Group& group) const
  {
    const std::uint32_t end = rows.index[count - 1U];
    const std::uint32_t capacity = mParameters->follow ? static_cast<std::uint32_t>(kMaxClusters) - (count - 1U) : 1U;
    auto box = trackTimeBox(mFrame->measurements(), rows, static_cast<std::uint8_t>(count - 1U));
    std::uint32_t walked[kMaxClusters];
    const std::uint64_t heldLayers = heldInactive(rows, count - 1U, proposal);
    const std::uint8_t taken = walk(start, PropagationDirection::kForward, walked, capacity, box, heldLayers, observe, end, group);
    if (taken == 0U) {
      return false;
    }
    // The end taken again: the rows beyond it are the follow outward;
    // another row: the walk's rows replace the end.
    const bool replaced = walked[0] != end;
    const std::uint32_t first = replaced ? 0U : 1U;
    if (replaced) {
      countWork(WorkCounter::kFollowEndReplaced);
      proposal.replacedEnd = end;
    }
    for (std::uint32_t i = first; i < taken; ++i) {
      proposal.outward[i - first] = walked[i];
    }
    proposal.outwardCount = static_cast<std::uint8_t>(taken - first);
    const std::uint32_t held = count - (replaced ? 1U : 0U) + proposal.outwardCount;
    walkInward(inner, rows.index[0], mParameters->follow && held < kMaxClusters ? static_cast<std::uint32_t>(kMaxClusters) - held : 0U, box, heldLayers, proposal, observe, group);
    return true;
  }

  /// Where the inward follow starts: \a inner, with \a BeamLine updated with
  /// the beam line, \a innermost the track's innermost row -- the nominal
  /// update, or, where the iteration has a kink allowance and the nominal
  /// update is above the gate, the update with the scattering inside scaled
  /// by it, as `clearsBeamLine` admitted the road.
  [[nodiscard]] GPUhdi() BoundTrackState inwardStart(const BoundTrackState& inner, std::uint32_t innermost) const
  {
    if constexpr (BeamLine) {
      BoundTrackState constrained = inner;
      const auto prediction = predictBeamLine(mPropagator, *mParameters, *mMaterial, mGeometry->planeOf(inner), pairPointAt(mFrame->measurements(), innermost), inner);
      const auto update = updateWithBeamLine(*mParameters, prediction, constrained);
      if (update.ok && !(update.chi2 <= mParameters->beamChi2Max) && mParameters->kinkScale2 > 1.F) {
        constrained = inner;
        static_cast<void>(updateWithBeamLine(*mParameters, prediction, constrained, mParameters->kinkScale2));
      }
      return constrained;
    } else {
      static_cast<void>(innermost);
      return inner;
    }
  }

  /// One direction of a follow from \a start: at most \a capacity rows into
  /// \a rows, \a box the track's time box, \a own a row of the track its
  /// searches take as free (with
  /// \a DecidesEnd); how many rows. The lanes of \a group walk it together.
  template <typename Observer, typename Group = SoloGroup>
  GPUhdi() GPUAlwaysInline() std::uint8_t walk(const BoundTrackState& start, PropagationDirection direction, std::uint32_t* rows, std::uint32_t capacity, TimeEstPS& box, std::uint64_t heldLayers, const Observer& observe, std::uint32_t own = std::numeric_limits<std::uint32_t>::max(), const Group& group = {}) const
  {
    observe.begin(direction);
    Branch branch;
    branch.own = own;
    branch.state = static_cast<const BoundParameters&>(start); // the covariance travels as its factors
    branch.factors = factorCovariance(start.cov);
    auto end = FollowEnd::kCapacity;
    // Whether the latest step saw a row near where it took none (the unconfirmed drop, below).
    bool nearMiss = false;
    [[maybe_unused]] follow::BeamEndHistoryOf<BeamEnd> history; // NOLINT(misc-const-correctness): written in a dependent call the check does not see
    // The inactive layers' tail resumes from the state and time box before
    // the last row when the unconfirmed drop removes it.
    [[maybe_unused]] BoundParameters beforeLast = branch.state;
    [[maybe_unused]] CovarianceFactors beforeLastFactors = branch.factors;
    [[maybe_unused]] TimeEstPS beforeLastBox = box;
    while (branch.count < capacity) {
      const LayerID current = mGeometry->surface(branch.state.surface).layer;
      observe.step(current, branch.state, branch.factors, branch.count);
      follow::Target targets[follow::maxDegree<Det>];
      const std::uint32_t nTargets = targetsOf(branch.state, direction, targets);
      if (nTargets == 0U) {
        end = FollowEnd::kNoTarget;
        break;
      }
      branch.crossed = false;
      nearMiss = false;
      Child child;
      if (nominal(branch, box, direction, targets, nTargets, child, nearMiss, observe, group)) {
        observe.taken(child.row, child.chi2);
        if constexpr (BeamEnd) {
          if (direction == PropagationDirection::kBackward) {
            history.remember(branch.state, branch.factors);
          }
        }
        if constexpr (Detector<Det>::inactiveLayers != 0U) {
          beforeLast = branch.state;
          beforeLastFactors = branch.factors;
          beforeLastBox = box;
        }
        take(branch, box, child);
        continue;
      }
      end = branch.holes > mParameters->followHoles ? FollowEnd::kHoles : FollowEnd::kNoRow;
      break;
    }
    // The unconfirmed drop: the last step crossed a sensor
    // without a row; the last row, which nothing confirms, goes -- unless,
    // inward, the step saw a row near (inward that row is almost always the
    // particle's own, outward a noise row often enough to add fakes, so
    // the outward unconfirmed drop stands).
    // A last row under the iteration's followKeepChi2 stays (the outward unconfirmed drop
    // drops mostly the particle's own rows, the wrong ones at the higher
    // chi2).
    const bool evidence = (direction == PropagationDirection::kBackward && nearMiss) || branch.lastChi2 < mParameters->followKeepChi2;
    // With EndDensity also an
    // outward last row of a high density stays (`densityKept`).
    [[maybe_unused]] bool dropped = false;
    if (branch.crossed && !evidence && branch.count != 0U && !densityKept(branch, direction)) {
      --branch.count;
      dropped = true;
      countWork(WorkCounter::kFollowUnconfirmed);
      observe.unconfirmed(branch.rows[branch.count]);
    }
    if constexpr (EndDensity) {
      // At every exit: a row taken as the last the capacity allows leaves
      // the loop before the search for targets. Never where the unconfirmed drop removed a
      // row: a layer without an inward link has no targets, so the walk
      // crossed nothing there.
      if (refusesLayoutEnd(branch, direction)) {
        --branch.count;
        countWork(WorkCounter::kFollowEndRefused);
        observe.refusedEnd(branch.rows[branch.count]);
      }
    }
    if constexpr (BeamEnd) {
      // After the unconfirmed drop and the layout's end, on the row that stands last.
      if (refusesBeamEnd(branch, direction, history)) {
        --branch.count;
        countWork(WorkCounter::kFollowBeamEndRefused);
        observe.refusedBeamEnd(branch.rows[branch.count]);
      }
    }
    if constexpr (Detector<Det>::inactiveLayers != 0U) {
      // The inactive layers' tail: outward, where the walk ended before them
      // -- a hole, no row, or no link on -- one nominal step after another
      // from its last state onto the inactive layers it has no row on,
      // nearest crossing first (`inactiveLayerTail`), until a step takes none.
      // Where the unconfirmed drop removed the last row, from the state
      // before that row.
      if (dropped) {
        branch.state = beforeLast;
        branch.factors = beforeLastFactors;
        box = beforeLastBox;
      }
      if (direction == PropagationDirection::kForward && end != FollowEnd::kCapacity) {
        while (branch.count < capacity) {
          Child child;
          bool nearTail = false;
          if (!inactiveLayerTail(branch, heldLayers, box, direction, child, nearTail, observe, group)) {
            break;
          }
          countWork(WorkCounter::kFollowInactiveTail);
          observe.taken(child.row, child.chi2);
          take(branch, box, child);
        }
      }
    }
    observe.end(end);
    for (std::uint32_t i = 0; i < branch.count; ++i) {
      rows[i] = branch.rows[i];
    }
    return static_cast<std::uint8_t>(branch.count);
  }

 private:
  /// EndDensity's density keep, asked where the unconfirmed drop would remove \a branch's last row:
  /// whether, outward, it has its density g
  /// above the iteration's `followKeepDensity` -- in logs, chi2 + log|V|
  /// below -2 log(2 pi key). A determinant that is not positive never keeps.
  /// Without EndDensity, false.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool densityKept(const Branch& branch, PropagationDirection direction) const
  {
    if constexpr (EndDensity) {
      const Scalar key = mParameters->followKeepDensity;
      const bool kept = direction == PropagationDirection::kForward && key > 0.F && branch.lastDeterminant > 0.F && branch.lastChi2 + logOf(branch.lastDeterminant) < -2.F * logOf(kTwoPi<Scalar> * key);
      countWork(WorkCounter::kFollowDensityKept, kept ? 1U : 0U);
      return kept;
    } else {
      static_cast<void>(branch);
      static_cast<void>(direction);
      return false;
    }
  }

  /// EndDensity's layout-end refusal: whether \a branch, inward, stands at the layout's
  /// end -- no inward link from its layer, not a search whose every
  /// prediction failed -- on a row with sqrt|V| above the
  /// iteration's `followEndMaxSqrtDet`.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool refusesLayoutEnd(const Branch& branch, PropagationDirection direction) const
  {
    const Scalar key = mParameters->followEndMaxSqrtDet;
    if (direction != PropagationDirection::kBackward || !(key > 0.F) || branch.count == 0U) {
      return false;
    }
    const LayerID current = mGeometry->surface(branch.state.surface).layer;
    return mParameters->allInLinkOffsets[current] == mParameters->allInLinkOffsets[current + 1U] && branch.lastDeterminant > key * key;
  }

  /// BeamEnd's rule: whether \a branch, inward, ends on a row (taken from
  /// \a history's state) whose chi2 against that step's
  /// prediction updated with the beam line is above the gate. The
  /// prediction carried as the step carried it, the update at the row's
  /// sensor from the row's pair point (`updateWithBeamLine`), the residual's
  /// covariance in double as the update's. A prediction or an update that
  /// fails, or a determinant that is not positive, refuses nothing.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool refusesBeamEnd(const Branch& branch, PropagationDirection direction, const follow::BeamEndHistory& history) const
  {
    if (direction != PropagationDirection::kBackward || branch.count == 0U) {
      return false;
    }
    // The rows the unconfirmed drop and the layout's end dropped: at most one of the two acts.
    const std::uint32_t dropped = history.taken - branch.count;
    if (dropped > 1U) {
      return false;
    }
    const std::uint32_t row = branch.rows[branch.count - 1U];
    const auto measurement = measurementAt(*mFrame, row);
    const auto& plane = mGeometry->planeOf(measurement);
    BoundTrackState predicted;
    CovarianceFactors carried;
    if (!follow::transportCandidate(mPropagator, *mGeometry, *mMaterial, plane, direction, history.state[dropped], history.factors[dropped], predicted, carried)) {
      return false;
    }
    covarianceFromFactors(carried, predicted.cov);
    if (!updateWithBeamLine(mPropagator, *mParameters, *mMaterial, plane, pairPointAt(mFrame->measurements(), row), predicted).ok) {
      return false;
    }
    const auto v00 = static_cast<ExactScalar>(predicted.variance(kLoc0)) + static_cast<ExactScalar>(measurement.varLoc0);
    const auto v11 = static_cast<ExactScalar>(predicted.variance(kLoc1)) + static_cast<ExactScalar>(measurement.varLoc1);
    const auto v01 = static_cast<ExactScalar>(predicted.covAt(kLoc0, kLoc1)) + static_cast<ExactScalar>(measurement.covLoc01);
    const ExactScalar det = (v00 * v11) - (v01 * v01);
    const auto r0 = static_cast<ExactScalar>(measurement.loc0) - static_cast<ExactScalar>(predicted.loc0());
    const auto r1 = static_cast<ExactScalar>(measurement.loc1) - static_cast<ExactScalar>(predicted.loc1());
    return det > 0. && ((r0 * r0 * v11) - (2. * r0 * r1 * v01) + (r1 * r1 * v00)) / det > static_cast<ExactScalar>(mParameters->clusterChi2Max);
  }

  /// What a step's search reads (`follow::visit`).
  [[nodiscard]] GPUhdi() GPUAlwaysInline() follow::StepContext<Det, Material> context() const { return {.parameters = mParameters, .frame = mFrame, .geometry = mGeometry, .material = mMaterial, .propagator = &mPropagator, .used = mUsed}; }

  /// The window of a step from \a state (its covariance \a factors) to
  /// \a target: the track's time \a box, narrowed to the prediction
  /// (`follow::predictedWindow`), its \a reach and the material it crossed, \a traversed,
  /// beside. False where the prediction fails. Every lane of \a group makes
  /// the same window; they look its material up together
  /// (`detail::GroupMaterial`).
  template <typename Group>
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool targetWindow(const BoundParameters& state, const CovarianceFactors& factors, const TimeEstPS& box, const follow::Target& target, PropagationDirection direction, CandidateWindow& window, Scalar (&reach)[2], detail::StepMaterial& traversed, const Group& group) const
  {
    Scalar phiCentre = 0.F;
    Scalar secondCentre = 0.F;
    follow::windowCentre(*mParameters, target.layer, target.prediction, phiCentre, secondCentre);
    window = {
      .earliestPS = static_cast<std::int64_t>(box.lower()) - static_cast<std::int64_t>(mMaxCandidateWidthPS[target.layer]) + 1,
      .latestPS = static_cast<std::int64_t>(box.upper()) - 1,
      .afterPS = box.lower(),
      .phiCentre = phiCentre,
      .phiHalfWidth = kPi<Scalar>,
      .secondCentre = secondCentre,
      .secondHalfWidth = follow::kOpenBound,
    };
    if constexpr (Group::kLanes == 1U) {
      return follow::predictedWindow(mPropagator, *mParameters, *mGeometry, *mMaterial, state, factors, target.layer, direction, target.prediction, window, reach, &traversed);
    } else {
      const detail::GroupMaterial<Material, Group> material{.lut = mMaterial, .group = &group};
      return follow::predictedWindow(mPropagator, *mParameters, *mGeometry, material, state, factors, target.layer, direction, target.prediction, window, reach, &traversed);
    }
  }

  /// The layers a step from \a state may reach in \a direction, into
  /// \a targets, nearest crossing first; how many.
  [[nodiscard]] GPUhdi() std::uint32_t targetsOf(const BoundParameters& state, PropagationDirection direction, follow::Target* targets) const
  {
    const bool forward = direction == PropagationDirection::kForward;
    const LayerID current = mGeometry->surface(state.surface).layer;
    const auto* offsets = forward ? mParameters->allOutLinkOffsets : mParameters->allInLinkOffsets;
    const auto* indices = forward ? mParameters->allOutLinkIndices : mParameters->allInLinkIndices;
    // The state's part of every link's crossing, once.
    const auto start = layerCrossingStart(mPropagator, mGeometry->planeOf(state), state);
    const Scalar* here = start.position;
    std::uint32_t nTargets = 0U;
    for (auto at = offsets[current]; at < offsets[current + 1U]; ++at) {
      const auto link = indices[at];
      const LayerID layer = forward ? mParameters->linkTargets[link] : mParameters->linkSources[link];
      const auto prediction = predictLayerCrossing(*mParameters, start, layer, direction);
      if (prediction.status != PropagationStatus::kOk) {
        continue;
      }
      Scalar point[3]{};
      follow::crossingPoint(*mParameters, layer, prediction, point);
      const Scalar dx = point[0] - here[0];
      const Scalar dy = point[1] - here[1];
      const Scalar dz = point[2] - here[2];
      const follow::Target target{.distance2 = (dx * dx) + (dy * dy) + (dz * dz), .layer = layer, .prediction = prediction};
      // With one link there is nothing to sort; compiling the loop out also
      // spares GCC's -Warray-bounds a subscript it cannot prove unreachable.
      // The link tables are `outLinks<Det>`/`inLinks<Det>`'s (Iteration.h), so
      // a layer's links never outnumber `maxDegree<Det>`, the targets' size;
      // the analyzer, not seeing the tables' origin, takes them as any.
      // NOLINTBEGIN(clang-analyzer-security.ArrayBound)
      if constexpr (follow::maxDegree<Det> > 1U) {
        std::uint32_t slot = nTargets++;
        while (slot > 0U && (target.distance2 < targets[slot - 1U].distance2 || (!(targets[slot - 1U].distance2 < target.distance2) && target.layer < targets[slot - 1U].layer))) {
          targets[slot] = targets[slot - 1U];
          --slot;
        }
        targets[slot] = target;
      } else {
        targets[nTargets++] = target;
      }
      // NOLINTEND(clang-analyzer-security.ArrayBound)
    }
    return nTargets;
  }

  /// The nominal pass's walk over \a branch's \a nTargets \a targets: per
  /// target its window (`predictedWindow`; a target the state cannot be
  /// carried to is passed, the track taken not to cross it), and
  /// \a search(target, window, traversed) with the material the window
  /// crossed, whether it found a row there; a target
  /// without one crossed in a sensor is a hole, counted on the branch, past
  /// the iteration's holes the pass ends. Whether a search found a row. The
  /// lanes of \a group make each window together (`targetWindow`).
  template <typename Observer, typename Group, typename Search>
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool overTargets(Branch& branch, const TimeEstPS& box, PropagationDirection direction, const follow::Target* targets, std::uint32_t nTargets, const Observer& observe, const Group& group, const Search& search) const
  {
    for (std::uint32_t t = 0; t < nTargets; ++t) {
      const auto& target = targets[t]; // NOLINT(clang-analyzer-security.ArrayBound): nTargets <= maxDegree<Det> (`targetsOf`)
      CandidateWindow window;
      Scalar reach[2]{};
      detail::StepMaterial traversed;
      if (!targetWindow(branch.state, branch.factors, box, target, direction, window, reach, traversed, group)) {
        continue;
      }
      observe.target(target.layer, window);
      if (search(target, window, traversed)) {
        return true;
      }
      const bool crossed = crossesSensor(branch.state, target, reach, direction);
      observe.hole(target.layer, crossed);
      branch.crossed = branch.crossed || crossed;
      if (crossed && ++branch.holes > mParameters->followHoles) {
        return false;
      }
    }
    return false;
  }

  /// The inactive layers (\a Det's `inactiveLayers`) \a rows' first \a count
  /// and \a proposal's interior row lie on: a walk's tail adds no row there.
  [[nodiscard]] GPUhdi() std::uint64_t heldInactive(const ClusterList<kMaxClusters>& rows, std::uint32_t count, const Proposal& proposal) const
  {
    std::uint64_t held = 0U;
    if constexpr (Detector<Det>::inactiveLayers != 0U) {
      const auto& measurements = mFrame->measurements();
      const auto add = [&](std::uint32_t row) { held |= std::uint64_t{1} << mGeometry->surface(measurements.surface[row]).layer; };
      for (std::uint32_t i = 0; i < count; ++i) {
        add(rows.index[i]);
      }
      if (proposal.interior != ClusterList<kMaxClusters>::kNoCluster) {
        add(proposal.interior);
      }
      held &= Detector<Det>::inactiveLayers;
    }
    return held;
  }

  /// The interior row: the first inactive layer between two
  /// consecutive rows of the track, \a rows' i and i + 1 -- a link from the
  /// first's layer onto it, and one from it onto the second's -- which a road
  /// crosses without a row on it, the follow never visiting it. From the
  /// track's forward pass cut at row i (its first i + 1 rows, at least
  /// three), or else its inward pass from \a outer, its outward state, cut at
  /// row i + 1 (at least three rows), a nominal step onto it, with the
  /// track's time \a box: the free row with the best predicted chi2 in the
  /// gate, as the follow takes a row, goes into \a proposal before row i + 1
  /// and narrows \a box where the track's fit with it raises the chi2 by no
  /// more than the gate. Whether one did.
  [[nodiscard]] GPUhdi() bool interior(const BoundTrackState& outer, const ClusterList<kMaxClusters>& rows, std::uint32_t count, TimeEstPS& box, Proposal& proposal) const
  {
    const auto& measurements = mFrame->measurements();
    const auto linked = [&](LayerID from, LayerID to) {
      for (auto at = mParameters->allOutLinkOffsets[from]; at < mParameters->allOutLinkOffsets[from + 1U]; ++at) {
        if (mParameters->linkTargets[mParameters->allOutLinkIndices[at]] == to) {
          return true;
        }
      }
      return false;
    };
    for (std::uint32_t i = 0U; i + 1U < count; ++i) {
      const LayerID before = mGeometry->surface(measurements.surface[rows.index[i]]).layer;
      const LayerID after = mGeometry->surface(measurements.surface[rows.index[i + 1U]]).layer;
      for (std::uint64_t left = Detector<Det>::inactiveLayers; left != 0U; left &= left - 1U) {
        const auto layer = static_cast<LayerID>(o2::gpu::CAMath::Ctz64(left));
        if (layer == before || layer == after || !linked(before, layer) || !linked(layer, after)) {
          continue;
        }
        const bool forward = i + 1U >= 3U;
        if (!forward && count - i - 1U < 3U) {
          return false;
        }
        const auto direction = forward ? PropagationDirection::kForward : PropagationDirection::kBackward;
        const RowChain<Det> chain{.frame = mFrame, .rows = rows.index};
        BoundTrackState state;
        Scalar chi2 = 0.F;
        std::int16_t failed = -1;
        std::int16_t kinkAt = -1;
        int kinks = 0;
        // A pass takes its start's values alone (`runPass`): the whole
        // chain's helix forward, the outward state inward.
        const BoundTrackState start = forward ? makeChainSeed(mPropagator, *mGeometry, chain, static_cast<int>(count)) : outer;
        const RowChain<Det> side{.frame = mFrame, .rows = forward ? rows.index : rows.index + i + 1U};
        const auto sideCount = static_cast<int>(forward ? i + 1U : count - i - 1U);
        const auto span = detail::spanMaterial(*mGeometry, *mMaterial, side, sideCount);
        if (detail::runPass(mPropagator, *mGeometry, span, side, sideCount, direction, start, state, chi2, failed, KinkBudget{}, kinks, kinkAt, NoFitObserver{}) != PropagationStatus::kOk) {
          return false;
        }
        const auto prediction = predictLayerCrossing(mPropagator, *mParameters, mGeometry->planeOf(state), state, layer, direction);
        if (prediction.status != PropagationStatus::kOk) {
          return false;
        }
        const follow::Target target{.distance2 = 0.F, .layer = layer, .prediction = prediction};
        const auto factors = factorCovariance(state.cov);
        CandidateWindow window;
        Scalar reach[2]{};
        detail::StepMaterial traversed;
        if (!targetWindow(state, factors, box, target, direction, window, reach, traversed, SoloGroup{})) {
          return false;
        }
        countWork(WorkCounter::kFollowInteriorSearches);
        // The row the follow would take there: the best predicted chi2 in
        // the gate, kept where the track's fit with it raises its chi2 by no
        // more than the gate.
        std::uint32_t best = ClusterList<kMaxClusters>::kNoCluster;
        Scalar bestChi2 = 0.F;
        follow::visitInline<false, false, false, Shape>(
          context(), state, factors, 0.F, 0.F, layer, direction, window, traversed, 0U, WorkCounter::kFollowUpdates, SoloGroup{}, [](std::uint32_t /*row*/, StepOutcome /*outcome*/, Scalar /*chi2*/) {},
          [&](const auto& candidate) {
            if (best == ClusterList<kMaxClusters>::kNoCluster || candidate.chi2 < bestChi2) {
              best = candidate.row;
              bestChi2 = candidate.chi2;
            }
          });
        if (best == ClusterList<kMaxClusters>::kNoCluster) {
          countWork(WorkCounter::kFollowInteriorEmpty);
          return false;
        }
        const auto budget = kinkBudget(*mParameters);
        const auto plain = fitMeasurements(mPropagator, *mGeometry, *mMaterial, chain, static_cast<int>(count), budget);
        const follow::InsertedChain<Det> with{.frame = mFrame, .rows = rows.index, .at = i + 1U, .row = best};
        const auto fit = fitMeasurements(mPropagator, *mGeometry, *mMaterial, with, static_cast<int>(count + 1U), budget);
        if (plain.status != PropagationStatus::kOk || fit.status != PropagationStatus::kOk || fit.chi2 - plain.chi2 > mParameters->clusterChi2Max) {
          countWork(WorkCounter::kFollowInteriorRefused);
          return false;
        }
        proposal.interior = best;
        proposal.interiorAt = static_cast<std::uint8_t>(i + 1U);
        box += timeBoxAt(measurements, best);
        return true;
      }
    }
    return false;
  }

  /// The inactive layers' tail of \a branch's walk (\a Det's
  /// `inactiveLayers`): the inactive layers neither it nor the track
  /// (\a heldLayers, `heldInactive`) holds a row on and its state crosses
  /// ahead, nearest first, searched as a nominal step is (`nominal`), the
  /// holes not counted. Whether it found a row, into \a child.
  template <typename Observer, typename Group>
  [[nodiscard]] GPUhdi() bool inactiveLayerTail(Branch& branch, std::uint64_t heldLayers, const TimeEstPS& box, PropagationDirection direction, Child& child, bool& nearMiss, const Observer& observe, const Group& group) const
  {
    follow::Target targets[follow::maxDegree<Det>];
    std::uint32_t nTargets = 0U;
    const auto& from = mGeometry->planeOf(branch.state);
    Scalar here[3]{};
    surfaceToGlobal(from, branch.state.loc0(), branch.state.loc1(), here);
    for (std::uint64_t left = Detector<Det>::inactiveLayers; left != 0U && nTargets < follow::maxDegree<Det>; left &= left - 1U) {
      const auto layer = static_cast<LayerID>(o2::gpu::CAMath::Ctz64(left));
      bool held = ((heldLayers >> layer) & 1U) != 0U || mGeometry->surface(branch.state.surface).layer == layer;
      for (std::uint32_t i = 0; i < branch.count; ++i) {
        held = held || mGeometry->surface(mFrame->measurements().surface[branch.rows[i]]).layer == layer;
      }
      if (held) {
        continue;
      }
      const auto prediction = predictLayerCrossing(mPropagator, *mParameters, from, branch.state, layer, direction);
      if (prediction.status != PropagationStatus::kOk) {
        continue;
      }
      Scalar point[3]{};
      follow::crossingPoint(*mParameters, layer, prediction, point);
      const Scalar d2 = ((point[0] - here[0]) * (point[0] - here[0])) + ((point[1] - here[1]) * (point[1] - here[1])) + ((point[2] - here[2]) * (point[2] - here[2]));
      std::uint32_t slot = nTargets++;
      while (slot > 0U && d2 < targets[slot - 1U].distance2) {
        targets[slot] = targets[slot - 1U];
        --slot;
      }
      targets[slot] = {.distance2 = d2, .layer = layer, .prediction = prediction};
    }
    if (nTargets == 0U) {
      return false;
    }
    Branch tail = branch;
    tail.holes = 0U;
    return nominal(tail, box, direction, targets, nTargets, child, nearMiss, observe, group);
  }

  /// The nominal pass of a step of \a branch over its \a nTargets
  /// \a targets (`overTargets`): the first target with a row in the gate
  /// gives \a child -- the row, or, with more than one, the lookahead's best confirmed
  /// (`confirm`). Whether it found a row; its updates count as the follow's.
  /// \a nearMiss is set where a target held a row within its window that
  /// the update refused above the gate (the unconfirmed drop reads it).
  ///
  /// The lanes of \a group search a window in step: after each step of the
  /// walk the rows that passed, in lane order -- the grid's order -- are
  /// counted for the whole window; the first stays on its lane until a second
  /// comes, then both are confirmed, and every later one as it comes, each
  /// confirmation by the whole group from the row's state on every lane.
  template <typename Observer, typename Group>
  [[nodiscard]] GPUhdi() bool nominal(Branch& branch, const TimeEstPS& box, PropagationDirection direction, const follow::Target* targets, std::uint32_t nTargets, Child& child, bool& nearMiss, const Observer& observe, const Group& group) const
  {
    return overTargets(branch, box, direction, targets, nTargets, observe, group, [&](const follow::Target& target, const CandidateWindow& window, const detail::StepMaterial& traversed) {
      std::uint32_t passing = 0U;
      std::uint32_t firstLane = 0U; // \a child holds the first passing row on this lane
      Confirmed confirmed;
      bool gated = false; // this lane saw a row refused above the gate
      follow::visit<DecidesEnd, EndDensity, Shape>(
        context(), branch.state, branch.factors, 0.F, 0.F, target.layer, direction, window, traversed, 0U, WorkCounter::kFollowUpdates, group,
        [&](std::uint32_t row, StepOutcome outcome, Scalar chi2) {
          gated = gated || outcome == StepOutcome::kChi2;
          observe.candidate(target.layer, row, outcome, chi2);
        },
        [&](const Child& candidate) {
          for (auto lanes = group.ballot(candidate.row != ClusterList<kMaxClusters>::kNoCluster); lanes != 0U; lanes &= lanes - 1U) {
            const auto from = static_cast<std::uint32_t>(o2::gpu::CAMath::Ctz64(lanes));
            if constexpr (Greedy) {
              // Every lane holds the best so far, each passing row shared as
              // it comes.
              const auto one = follow::shared(group, candidate, from);
              if (++passing == 1U || follow::beats(one.row, one.chi2, follow::Witness{.row = child.row, .chi2 = child.chi2})) {
                child = one;
              }
              continue;
            }
            if (++passing == 1U) {
              firstLane = from;
              if (group.lane() == from) {
                child = candidate;
              }
              continue;
            }
            if (passing == 2U) {
              confirm(branch, box, direction, follow::shared(group, child, firstLane), confirmed, observe, group);
            }
            confirm(branch, box, direction, follow::shared(group, candidate, from), confirmed, observe, group);
          }
        },
        branch.own);
      // Any lane's: a collective, every lane calls it.
      nearMiss = nearMiss || group.ballot(gated) != 0U;
      if (passing == 0U) {
        return false;
      }
      if constexpr (Greedy) {
        return true;
      }
      if (passing == 1U) {
        if constexpr (Group::kLanes > 1U) {
          child = follow::shared(group, child, firstLane);
        }
      } else {
        countWork(WorkCounter::kFollowLookaheads);
        observe.lookahead(passing);
        child = confirmed.child;
      }
      return true;
    });
  }

  /// The nominal pass of a step of \a branch as a confirmation sees it
  /// (`overTargets`): the first target with a row in the gate gives
  /// \a witness, its smallest chi2, the row breaking ties. Whether it found
  /// one; its updates count on \a counter. The lanes of \a group search a
  /// window each its own rows, and agree on the best after it.
  template <typename Group>
  [[nodiscard]] GPUhdi() bool witnessOf(Branch& branch, const TimeEstPS& box, PropagationDirection direction, const follow::Target* targets, std::uint32_t nTargets, follow::Witness& witness, WorkCounter counter, const Group& group) const
  {
    return overTargets(branch, box, direction, targets, nTargets, follow::NoObserver{}, group, [&](const follow::Target& target, const CandidateWindow& window, const detail::StepMaterial& traversed) {
      follow::Witness own;
      follow::visit<DecidesEnd, false, Shape>(
        context(), branch.state, branch.factors, 0.F, 0.F, target.layer, direction, window, traversed, 0U, counter, group, [](std::uint32_t /*row*/, StepOutcome /*outcome*/, Scalar /*chi2*/) {},
        [&](const follow::Child& candidate) {
          if (candidate.row != ClusterList<kMaxClusters>::kNoCluster && follow::beats(candidate.row, candidate.chi2, own)) {
            own = {.row = candidate.row, .chi2 = candidate.chi2};
          }
        },
        branch.own);
      witness = follow::bestOf(group, own);
      return witness.row != ClusterList<kMaxClusters>::kNoCluster;
    });
  }

  /// The lookahead: of the rows of a window passing
  /// the gate, where more than one does, the one the next step confirms, by
  /// its chi2 plus the confirming row's; where none is confirmed, the
  /// smallest chi2; the row breaks ties. `nominal` hands it each passing row
  /// in the window's order from the second on, the first with it.
  struct Confirmed {
    Child child; // its row kNoCluster
    bool confirmed = false;
    Scalar key = 0.F;
  };

  /// \a taken, a row passing the gate from \a branch, against \a best; the
  /// lanes of \a group together. \a observe sees the confirmation.
  template <typename Observer, typename Group>
  GPUhdi() void confirm(const Branch& branch, const TimeEstPS& box, PropagationDirection direction, const Child& taken, Confirmed& best, const Observer& observe, const Group& group) const
  {
    std::uint32_t next = ClusterList<kMaxClusters>::kNoCluster;
    Scalar witness = 0.F;
    const bool confirmed = confirmation(branch, box, direction, taken, next, witness, WorkCounter::kFollowUpdates, group);
    const Scalar key = confirmed ? taken.chi2 + witness : taken.chi2;
    observe.lookaheadConfirmation(taken.row, taken.chi2, next, witness);
    const bool better = best.child.row == ClusterList<kMaxClusters>::kNoCluster || (confirmed != best.confirmed ? confirmed : follow::beats(taken.row, key, follow::Witness{.row = best.child.row, .chi2 = best.key}));
    if (better) {
      best.child = taken;
      best.confirmed = confirmed;
      best.key = key;
    }
  }

  /// What confirms \a taken, a row passing the gate from \a branch: the
  /// next step's row by the nominal pass from past it, into \a next and
  /// \a chi2, its updates counted on \a counter. Whether it does. The lanes
  /// of \a group together.
  template <typename Group>
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool confirmation(const Branch& branch, const TimeEstPS& box, PropagationDirection direction, const Child& taken, std::uint32_t& next, Scalar& chi2, WorkCounter counter, const Group& group) const
  {
    Branch probe = branch;
    auto probeBox = box;
    take(probe, probeBox, taken);
    follow::Target targets[follow::maxDegree<Det>];
    const std::uint32_t nTargets = targetsOf(probe.state, direction, targets);
    if (nTargets == 0U) {
      return false;
    }
    follow::Witness witness;
    if (!witnessOf(probe, probeBox, direction, targets, nTargets, witness, counter, group)) {
      return false;
    }
    next = witness.row;
    chi2 = witness.chi2;
    return true;
  }

  /// \a branch past \a child's row; \a box narrowed by it.
  GPUhdi() GPUAlwaysInline() void take(Branch& branch, TimeEstPS& box, const Child& child) const
  {
    branch.rows[branch.count++] = child.row;
    branch.state = child.state;
    branch.factors = child.factors;
    branch.lastChi2 = child.chi2;
    if constexpr (EndDensity) {
      branch.lastDeterminant = child.determinant;
    }
    branch.crossed = false;
    box += timeBoxAt(mFrame->measurements(), child.row);
  }

  /// Whether \a state crosses a sensor of \a target's layer within its
  /// bounds in \a direction: carried to each sensor listed near its nominal
  /// crossing, as far from it as the layer's envelope \a reach lets the
  /// crossing move; a curved sensor is tested on its tangent plane at the
  /// crossing's angle.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool crossesSensor(const BoundParameters& state, const follow::Target& target, const Scalar* reach, PropagationDirection direction) const
  {
    // An inactive layer is never a hole, nor a crossing the unconfirmed drop reads (LayerSpec.h).
    if constexpr (Detector<Det>::inactiveLayers != 0U) {
      if (Detector<Det>::isInactiveLayer(target.layer)) {
        return false;
      }
    }
    countWork(WorkCounter::kFollowHoleTests);
    const auto& from = mGeometry->planeOf(state);
    Scalar point[3]{};
    follow::crossingPoint(*mParameters, target.layer, target.prediction, point);
    const auto& propagator = pointPropagator(mPropagator, from, state);
    const Scalar curvature = state.qOverP() * propagator.bzKiloGauss * static_cast<Scalar>(constants::math::B2C);
    const bool crossed = mCells.anySensor(target.layer, target.prediction.phi, reach[0], target.prediction.second, reach[1], [&](std::uint32_t id) {
      const TrackingSurface plane = mGeometry->planeAt(static_cast<SensorID>(id), point);
      BoundParameters at = state;
      const auto result = propagate<CrossingRule::kFirstAhead, Scalar>(propagator, at, from, plane, direction);
      return result.status == PropagationStatus::kOk && !follow::beyondHalfTurn(curvature * result.pathLength);
    });
    countWork(WorkCounter::kFollowHoles, crossed ? 1U : 0U);
    return crossed;
  }

  const IterationParameters<Det>* mParameters = nullptr;
  const TimeFrame<Det>* mFrame = nullptr;
  const TrackingGeometry<Det>* mGeometry = nullptr;
  const Material* mMaterial = nullptr;
  Propagator mPropagator;
  TrackView mTracks;
  const std::uint8_t* mUsed = nullptr;
  const std::uint32_t* mMaxCandidateWidthPS = nullptr;
  SensorCellsView mCells;
  Proposal* mProposals = nullptr;
  bool mInward = true;
  std::uint32_t mTrackOffset = 0U;
  const BoundTrackState* mEndStarts = nullptr;
};

/// Whether the iteration drops tracks after the follow: it has a length
/// after the follow, or it follows and has a pT floor.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool dropsAfterFollow(const IterationParameters<Det>& parameters)
{
  bool floors = false;
  for (const auto floor : parameters.minPtByClusters) {
    floors = floors || floor > 0.F;
  }
  return parameters.minClustersAfterFollow != 0U || (parameters.follow && floors);
}

/// Whether a track of \a count clusters, \a state its refit as the follow
/// left it, stays after the follow: at least the iteration's length after
/// it, and, where the iteration follows, at or above its pT floor for that
/// count.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool keptAfterFollow(const IterationParameters<Det>& parameters, std::uint32_t count, const BoundTrackState& state)
{
  return count >= parameters.minClustersAfterFollow && (!parameters.follow || clearsPtFloor(parameters, count, state));
}

namespace follow
{

/// The road fit's outward pass of each track's rows less its outermost,
/// a pass of its own before the follow
/// decides the end again (`FollowStage`'s DecidesEnd): \a starts gets the
/// state, or one on no surface where the track has fewer than four rows --
/// three left, a degree of freedom -- or is as long as a chain can be, or the
/// fit refuses them. Not with a kink allowance: the road fit takes its last
/// step with it (`fitOutward`'s LastStepAllowance), this one does not.
template <DetectorLayout Det, typename Material>
struct EndFit {
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using TrackView = typename BoundedColumns<TrackOf<kMaxClusters>::template RecordT>::ConstView;

  const IterationParameters<Det>* parameters = nullptr;
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  const Material* material = nullptr;
  Propagator propagator;
  TrackView tracks;
  BoundTrackState* starts = nullptr;
  std::uint32_t trackOffset = 0U;

  GPUhdi() void operator()(std::uint32_t element) const
  {
    const std::uint32_t track = trackOffset + element;
    starts[element].surface = kInvalidSensorID;
    const std::uint32_t count = tracks.count[track];
    if (tracks.status[track] != PropagationStatus::kOk || count < 4U || count >= kMaxClusters) {
      return;
    }
    countWork(WorkCounter::kFollowEndFits);
    const auto fit = fitOutward(propagator, *geometry, *material, RowChain<Det>{.frame = frame, .rows = tracks.clusters[track].index}, static_cast<int>(count - 1U), parameters->roadChi2MaxByNdf, kinkBudget(*parameters));
    if (fit.status == RoadFitStatus::kAccepted) {
      starts[element] = fit.state;
    }
  }
};

template <std::size_t K>
struct Bid {
  const FollowProposal<K>* proposals = nullptr;
  std::uint32_t* slots = nullptr;
  std::uint32_t trackOffset = 0U;
  std::uint32_t trackCount = 0U;

  GPUhdi() void operator()(std::uint32_t track) const
  {
    const std::uint32_t priority = trackCount - (trackOffset + track);
    const auto& proposal = proposals[track];
    for (std::uint32_t i = 0; i < proposal.outwardCount; ++i) {
      o2::gpu::CAMath::AtomicMax(slots + proposal.outward[i], priority);
    }
    for (std::uint32_t i = 0; i < proposal.inwardCount; ++i) {
      o2::gpu::CAMath::AtomicMax(slots + proposal.inward[i], priority);
    }
    if (proposal.interior != ClusterList<K>::kNoCluster) {
      o2::gpu::CAMath::AtomicMax(slots + proposal.interior, priority);
    }
  }
};

/// Keep each direction only through the first row on which another track has
/// the higher output-order priority. A later row won by this track is still
/// dropped: one bidding round does not retry a truncated continuation.
/// \a earlier counts the earlier iterations' follow apart.
template <std::size_t K>
struct Resolve {
  FollowProposal<K>* proposals = nullptr;
  const std::uint32_t* slots = nullptr;
  std::uint32_t trackOffset = 0U;
  std::uint32_t trackCount = 0U;
  bool earlier = false;

  GPUhdi() void operator()(std::uint32_t track) const
  {
    const std::uint32_t priority = trackCount - (trackOffset + track);
    auto& proposal = proposals[track];
    const bool interior = proposal.interior != ClusterList<K>::kNoCluster;
    const std::uint32_t proposed = static_cast<std::uint32_t>(proposal.outwardCount) + proposal.inwardCount + (interior ? 1U : 0U);
    proposal.outwardCount = kept(proposal.outward, proposal.outwardCount, priority);
    proposal.inwardCount = kept(proposal.inward, proposal.inwardCount, priority);
    if (interior && slots[proposal.interior] != priority) {
      proposal.interior = ClusterList<K>::kNoCluster;
    }
    const std::uint32_t left = static_cast<std::uint32_t>(proposal.outwardCount) + proposal.inwardCount + (proposal.interior != ClusterList<K>::kNoCluster ? 1U : 0U);
    countWork(earlier ? WorkCounter::kFollowEarlierProposed : WorkCounter::kFollowProposed, proposed);
    countWork(earlier ? WorkCounter::kFollowEarlierOutbid : WorkCounter::kFollowOutbid, proposed - left);
  }

 private:
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint8_t kept(const std::uint32_t* rows, std::uint8_t count, std::uint32_t priority) const
  {
    for (std::uint8_t i = 0U; i < count; ++i) {
      if (slots[rows[i]] != priority) {
        return i;
      }
    }
    return count;
  }
};

template <std::size_t K>
struct Combine {
  const ClusterList<K>* original = nullptr;
  const std::uint8_t* originalCount = nullptr;
  const FollowProposal<K>* proposals = nullptr;
  ClusterList<K>* rows = nullptr;
  std::uint8_t* counts = nullptr;
  std::uint32_t trackOffset = 0U;

  GPUhdi() void operator()(std::uint32_t element) const
  {
    const std::uint32_t track = trackOffset + element;
    const auto& proposal = proposals[element];
    ClusterList<K> combined;
    for (std::size_t i = 0; i < K; ++i) {
      combined.index[i] = ClusterList<K>::kNoCluster;
    }
    for (std::uint32_t i = 0; i < proposal.inwardCount; ++i) {
      combined.index[i] = proposal.inward[proposal.inwardCount - i - 1U];
    }
    // A replaced end goes with the outward rows that replace it.
    const bool replaced = proposal.replacedEnd != ClusterList<K>::kNoCluster && proposal.outwardCount != 0U;
    std::uint32_t kept = originalCount[track] - (replaced ? 1U : 0U);
    // An interior row goes before the track's row it was found before.
    const bool interior = proposal.interior != ClusterList<K>::kNoCluster;
    for (std::uint32_t i = 0, at = proposal.inwardCount; i < kept; ++i) {
      if (interior && i == proposal.interiorAt) {
        combined.index[at++] = proposal.interior;
      }
      combined.index[at++] = original[track].index[i];
    }
    if (interior && proposal.interiorAt >= kept) {
      combined.index[static_cast<std::uint32_t>(proposal.inwardCount) + kept] = proposal.interior;
    }
    kept += interior ? 1U : 0U;
    for (std::uint32_t i = 0; i < proposal.outwardCount; ++i) {
      combined.index[static_cast<std::uint32_t>(proposal.inwardCount) + kept + i] = proposal.outward[i];
    }
    rows[element] = combined;
    counts[element] = static_cast<std::uint8_t>(kept + proposal.inwardCount + proposal.outwardCount);
  }
};

/// Whether \a fit, a track refitted after its follow, passes the road's
/// chi2 cut (`Refit`, `TrimAmbiguous`).
template <DetectorLayout Det, typename Fit>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool passesRoadChi2(const IterationParameters<Det>& parameters, const Fit& fit)
{
  return fit.status == PropagationStatus::kOk && fit.chi2 <= parameters.roadChi2MaxByNdf[fit.ndf];
}

/// The extended refit of each followed track, all or nothing. \a Rule is
/// the fit's choice of crossing (`CrossingRule`), the iteration's
/// `crossingRule`. With \a Prefix, the iteration's `followPrefix`, a refit the cut
/// refuses or that fails is tried again on the follow's rows' prefixes in
/// the walks' order -- the outward rows first, then the inward ones --
/// the longest first, and the first the refit accepts is kept.
template <DetectorLayout Det, typename Material, CrossingRule Rule = CrossingRule::kFirstAhead, bool Prefix = false>
struct Refit {
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using TrackView = typename BoundedColumns<TrackOf<kMaxClusters>::template RecordT>::View;

  const IterationParameters<Det>* parameters = nullptr;
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  const Material* material = nullptr;
  Propagator propagator;
  TrackView tracks;
  FollowProposal<kMaxClusters>* proposals = nullptr;
  const ClusterList<kMaxClusters>* rows = nullptr;
  const std::uint8_t* counts = nullptr;
  std::uint32_t trackOffset = 0U;
  bool earlier = false;

  GPUhdi() void operator()(std::uint32_t element) const
  {
    const std::uint32_t track = trackOffset + element;
    auto& proposal = proposals[element];
    const std::uint32_t added = static_cast<std::uint32_t>(proposal.outwardCount) + proposal.inwardCount + (proposal.interior != ClusterList<kMaxClusters>::kNoCluster ? 1U : 0U);
    if (added == 0U) {
      return;
    }
    const std::uint8_t count = counts[element];
    if (!accept(track, rows[element], count, proposal.outwardCount, proposal.inwardCount)) {
      if constexpr (Prefix) {
        if (keepPrefix(element, track, proposal, count)) {
          return;
        }
      }
      countWork(earlier ? WorkCounter::kFollowEarlierRefitRejected : WorkCounter::kFollowRefitRejected, added);
      proposal.outwardCount = 0U;
      proposal.inwardCount = 0U;
      proposal.interior = ClusterList<kMaxClusters>::kNoCluster;
    }
  }

 private:
  /// The track with \a list's rows (\a n of them, \a outward and \a inward
  /// of them the follow's) refitted: whether the fit passes the road's chi2
  /// cut; then it is written as the track, its added counts grown -- added
  /// to what an earlier follow added: an earlier iteration's track is
  /// followed again with `followEarlier`.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool accept(std::uint32_t track, const ClusterList<kMaxClusters>& list, std::uint8_t n, std::uint8_t outward, std::uint8_t inward) const
  {
    const auto fit = fitMeasurements<Rule>(propagator, *geometry, *material, RowChain<Det>{.frame = frame, .rows = list.index}, static_cast<int>(n), kinkBudget(*parameters));
    if (!passesRoadChi2(*parameters, fit)) {
      return false;
    }
    writeFit(tracks, track, fit, list, n, trackTimeBox(frame->measurements(), list, n));
    tracks.addedOutward[track] = static_cast<std::uint8_t>(tracks.addedOutward[track] + outward);
    tracks.addedInward[track] = static_cast<std::uint8_t>(tracks.addedInward[track] + inward);
    return true;
  }

  /// `Prefix`'s second chance: the combined rows are the inward ones
  /// innermost first, the track's, the outward ones; a prefix of \a outward
  /// outward and \a inward inward rows is the track's rows with the first
  /// of each. With an interior row (`FollowStage::interior`) the prefixes
  /// with it come first, then those without it, the whole walk's among them.
  /// Whether one was kept, its counts (and interior row) then the proposal's.
  GPUhdi() GPUAlwaysInline() bool keepPrefix(std::uint32_t element, std::uint32_t track, FollowProposal<kMaxClusters>& proposal, std::uint8_t count) const
  {
    const auto outwardAdded = static_cast<std::uint32_t>(proposal.outwardCount);
    const auto inwardAdded = static_cast<std::uint32_t>(proposal.inwardCount);
    const bool interior = proposal.interior != ClusterList<kMaxClusters>::kNoCluster;
    const std::uint32_t kept = count - outwardAdded - inwardAdded;
    const auto& combined = rows[element];
    for (int pass = 0; pass < (interior ? 2 : 1); ++pass) {
      const bool withInterior = pass == 0;
      // The interior row's place among the combined rows, skipped without it.
      const std::uint32_t skip = withInterior ? count : inwardAdded + proposal.interiorAt;
      for (std::uint32_t length = outwardAdded + inwardAdded + (withInterior ? 0U : 1U); length-- > 1U;) {
        const auto outward = length < outwardAdded ? length : outwardAdded;
        const auto inward = length - outward;
        ClusterList<kMaxClusters> prefix{};
        std::uint32_t n = 0U;
        for (std::uint32_t i = inward; i-- > 0U;) {
          prefix.index[n++] = proposal.inward[i];
        }
        for (std::uint32_t i = 0; i < kept; ++i) {
          if (inwardAdded + i != skip) {
            prefix.index[n++] = combined.index[inwardAdded + i];
          }
        }
        for (std::uint32_t i = 0; i < outward; ++i) {
          prefix.index[n++] = proposal.outward[i];
        }
        if (!accept(track, prefix, static_cast<std::uint8_t>(n), static_cast<std::uint8_t>(outward), static_cast<std::uint8_t>(inward))) {
          continue;
        }
        countWork(earlier ? WorkCounter::kFollowEarlierRefitRejected : WorkCounter::kFollowRefitRejected, outwardAdded + inwardAdded - length + (interior && !withInterior ? 1U : 0U));
        proposal.outwardCount = static_cast<std::uint8_t>(outward);
        proposal.inwardCount = static_cast<std::uint8_t>(inward);
        if (!withInterior) {
          proposal.interior = ClusterList<kMaxClusters>::kNoCluster;
        }
        return true;
      }
    }
    return false;
  }
};

/// The trim of an ambiguous last row, after the refit, where the iteration has a `followTrimMargin`: a
/// track whose last row is one its follow added outward is refitted without it. Where that
/// fit passes the road's chi2 cut, the track one row shorter keeps the cuts
/// after the follow (`keptAfterFollow`), and from the fit's outward state two
/// rows or more of that row's layer pass the nominal gate (`passingOn`: free
/// rows, the row itself among them -- nothing is marked before `Mark`) with
/// the row itself not among them, or the smallest other chi2 less the row's
/// below the margin (a margin above the gate takes every such
/// row), the row is dropped and the track written with that fit; else the
/// track stays as the refit left it. The row dropped is not marked: it is free for the
/// iterations after. Not with `followEnd`, whose replaced end goes with the
/// outward rows.
template <DetectorLayout Det, typename Material, CrossingRule Rule = CrossingRule::kFirstAhead, bool Shape = false>
struct TrimAmbiguous {
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using TrackView = typename BoundedColumns<TrackOf<kMaxClusters>::template RecordT>::View;

  /// The search of the step's rows, as a nominal step makes it (with \a Shape,
  /// the iteration's shape test).
  FollowStage<Det, Material, false, false, false, false, false, Shape> search;
  const IterationParameters<Det>* parameters = nullptr;
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  const Material* material = nullptr;
  Propagator propagator;
  TrackView tracks;
  FollowProposal<kMaxClusters>* proposals = nullptr;
  std::uint32_t trackOffset = 0U;

  GPUhdi() void operator()(std::uint32_t element) const
  {
    const std::uint32_t track = trackOffset + element;
    auto& proposal = proposals[element];
    const std::uint32_t outward = proposal.outwardCount;
    const std::uint32_t count = tracks.count[track];
    if (outward == 0U || tracks.status[track] != PropagationStatus::kOk) {
      return;
    }
    // The outward rows are the track's last (`Combine`).
    ClusterList<kMaxClusters> rows = tracks.clusters[track];
    const std::uint32_t row = rows.index[count - 1U];
    // Preserve the last row when it comes from an inactive layer.
    if constexpr (Detector<Det>::inactiveLayers != 0U) {
      if (Detector<Det>::isInactiveLayer(geometry->surface(frame->measurements().surface[row]).layer)) {
        return;
      }
    }
    rows.index[count - 1U] = ClusterList<kMaxClusters>::kNoCluster;
    const auto shorter = static_cast<std::uint8_t>(count - 1U);
    countWork(WorkCounter::kFollowTrimFits);
    const auto fit = fitMeasurements<Rule>(propagator, *geometry, *material, RowChain<Det>{.frame = frame, .rows = rows.index}, static_cast<int>(shorter), kinkBudget(*parameters));
    if (!passesRoadChi2(*parameters, fit) || !keptAfterFollow(*parameters, shorter, fit.state)) {
      return;
    }
    const auto box = trackTimeBox(frame->measurements(), rows, shorter);
    const LayerID layer = geometry->surface(frame->measurements().surface[row]).layer;
    const auto passing = search.passingOn(fit.outwardState, factorCovariance(fit.outwardState.cov), layer, row, box, PropagationDirection::kForward);
    if (passing.count < 2U || (passing.own && !(passing.otherChi2 - passing.ownChi2 < parameters->followTrimMargin))) {
      return;
    }
    writeFit(tracks, track, fit, rows, shorter, box);
    tracks.addedOutward[track] = static_cast<std::uint8_t>(tracks.addedOutward[track] - 1U);
    proposal.outwardCount = static_cast<std::uint8_t>(outward - 1U);
    countWork(WorkCounter::kFollowTrimmed);
  }
};

template <std::size_t K>
struct Mark {
  const FollowProposal<K>* proposals = nullptr;
  std::uint8_t* used = nullptr;

  GPUhdi() void operator()(std::uint32_t track) const
  {
    const auto& proposal = proposals[track];
    // A replaced end, its replacement kept through the refit, is free again:
    // no other track proposed it, so no other one marks it.
    if (proposal.replacedEnd != ClusterList<K>::kNoCluster && proposal.outwardCount != 0U) {
      used[proposal.replacedEnd] = kFree;
    }
    for (std::uint32_t i = 0; i < proposal.outwardCount; ++i) {
      used[proposal.outward[i]] = kExclusive;
    }
    for (std::uint32_t i = 0; i < proposal.inwardCount; ++i) {
      used[proposal.inward[i]] = kExclusive;
    }
    if (proposal.interior != ClusterList<K>::kNoCluster) {
      used[proposal.interior] = kExclusive;
    }
  }
};

} // namespace follow

template <DetectorLayout Det>
class Follow;

namespace follow
{

/// What the follow's passes read and write (`runSequence`): the iteration's
/// parameters as the passes read them and as the host decides on them (one
/// object on the host, the device's copy and the host's on a device), the
/// track model, the tracks from \a trackOffset to \a trackCount, the used
/// mask (\a rowCount rows) and the passes' storage, one element per track
/// followed -- the end starts only where it decides the end again.
template <DetectorLayout Det, typename Material>
struct Inputs {
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using TrackView = typename BoundedColumns<TrackOf<kMaxClusters>::template RecordT>::View;
  using ConstTrackView = typename BoundedColumns<TrackOf<kMaxClusters>::template RecordT>::ConstView;

  const IterationParameters<Det>* parameters = nullptr;
  const IterationParameters<Det>* choices = nullptr;
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  const Material* material = nullptr;
  Propagator propagator;
  TrackView tracks;
  ConstTrackView constTracks;
  std::uint32_t trackOffset = 0U;
  std::uint32_t trackCount = 0U;
  bool inward = false;
  bool earlier = false;
  std::uint8_t* used = nullptr;
  std::size_t rowCount = 0U;
  const std::uint32_t* maxCandidateWidthPS = nullptr;
  SensorCellsView cells;
  FollowProposal<kMaxClusters>* proposals = nullptr;
  std::uint32_t* bids = nullptr;
  ClusterList<kMaxClusters>* rows = nullptr;
  std::uint8_t* counts = nullptr;
  BoundTrackState* endStarts = nullptr;
};

/// The host's launches of the follow's passes: each in its profile zone,
/// none of its own for the earlier iterations' tracks (`Follow::run`).
template <typename ForEach>
struct HostLaunch {
  ForEach* forEach;
  bool earlier = false;

  [[nodiscard]] profile::Phase phase(profile::Phase own) const { return earlier ? profile::Phase::kNone : own; }

  template <typename Body>
  void pass(profile::Phase own, std::uint32_t count, const Body& body) const
  {
    const profile::Zone zone{phase(own)};
    (*forEach)(count, body);
  }

  template <typename Stage>
  void propose(profile::Phase own, std::uint32_t count, const Stage& stage) const
  {
    pass(own, count, stage);
  }

  /// The bids zeroed in parallel, as the selection's slots are.
  void clearBids(std::uint32_t* bids, std::size_t rows) const
  {
    const profile::Zone zone{phase(profile::Phase::kFollowPropose)};
    (*forEach)(static_cast<std::uint32_t>(rows), [bids](std::uint32_t row) { bids[row] = 0U; });
  }
};

/// The follow over \a in's tracks, the host's and a device's alike: the bids
/// cleared, the end fit where the iteration decides the end again (on its
/// own tracks only), the proposals, then bid, resolve, combine, the refit,
/// the trim of the ambiguous last rows
/// (`followTrimMargin`, not with the end decided again, on the iteration's
/// own tracks) and the mark. Each choice of the iteration is its own
/// instantiation. \a launch runs a pass (`pass(phase, count, body)`), the
/// proposals (`propose(phase, count, stage)`: a device's lanes) and
/// `clearBids(bids, rows)`.
template <DetectorLayout Det, typename Material, typename Launch>
void runSequence(const Launch& launch, const Inputs<Det, Material>& in)
{
  using profile::Phase;
  constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  const IterationParameters<Det>& choices = *in.choices;
  const std::uint32_t elements = in.trackCount - in.trackOffset;
  launch.clearBids(in.bids, in.rowCount);
  // The beam line starts the inward follow where the iteration has the
  // constraint, and the end is decided
  // again where it has followEnd, on its own tracks
  // only.
  const bool end = choices.followEnd && !in.earlier;
  if (end) {
    launch.pass(Phase::kFollowPropose, elements,
                EndFit<Det, Material>{
                  .parameters = in.parameters,
                  .frame = in.frame,
                  .geometry = in.geometry,
                  .material = in.material,
                  .propagator = in.propagator,
                  .tracks = in.constTracks,
                  .starts = in.endStarts,
                  .trackOffset = in.trackOffset,
                });
  }
  // The iteration's choice among the passing rows, its end decisions on the
  // last row's density and against the beam line (with the beam line only).
  withBools(
    [&]<bool BeamLine, bool DecidesEnd, bool Greedy, bool EndDensity, bool BeamEnd, bool Shape>(std::bool_constant<BeamLine>, std::bool_constant<DecidesEnd>, std::bool_constant<Greedy>, std::bool_constant<EndDensity>, std::bool_constant<BeamEnd>, std::bool_constant<Shape>) {
      if constexpr (BeamLine || !BeamEnd) {
        const FollowStage<Det, Material, BeamLine, DecidesEnd, Greedy, EndDensity, BeamEnd, Shape> stage{in.parameters, in.frame, in.geometry, in.material, in.propagator, in.constTracks, in.used, in.maxCandidateWidthPS, in.cells, in.proposals, in.inward, in.trackOffset, DecidesEnd ? in.endStarts : nullptr};
        launch.propose(Phase::kFollowPropose, elements, stage);
      }
    },
    choices.beamConstraint, end, choices.followGreedy, Follow<Det>::endDensity(choices), choices.followBeamEnd, choices.shapeTest);
  launch.pass(Phase::kFollowResolve, elements, Bid<kMaxClusters>{in.proposals, in.bids, in.trackOffset, in.trackCount});
  launch.pass(Phase::kFollowResolve, elements, Resolve<kMaxClusters>{in.proposals, in.bids, in.trackOffset, in.trackCount, in.earlier});
  launch.pass(Phase::kFollowResolve, elements, Combine<kMaxClusters>{in.constTracks.clusters, in.constTracks.count, in.proposals, in.rows, in.counts, in.trackOffset});
  // The iteration's choice of crossing and of
  // the prefix.
  withCrossingRule(choices.crossingRule, [&]<CrossingRule Rule>(std::integral_constant<CrossingRule, Rule>) { withBools([&]<bool Prefix>(std::bool_constant<Prefix>) { launch.pass(Phase::kFollowRefit, elements, Refit<Det, Material, Rule, Prefix>{in.parameters, in.frame, in.geometry, in.material, in.propagator, in.tracks, in.proposals, in.rows, in.counts, in.trackOffset, in.earlier}); }, choices.followPrefix); });
  // The ambiguous last rows trimmed, on the iteration's own
  // tracks, before the marks.
  if (choices.followTrimMargin > 0.F && !end && !in.earlier) {
    withBools([&]<bool Shape>(std::bool_constant<Shape>) {
      const FollowStage<Det, Material, false, false, false, false, false, Shape> search{in.parameters, in.frame, in.geometry, in.material, in.propagator, in.constTracks, in.used, in.maxCandidateWidthPS, in.cells, nullptr, false, in.trackOffset};
      withCrossingRule(choices.crossingRule, [&]<CrossingRule Rule>(std::integral_constant<CrossingRule, Rule>) {
        launch.pass(Phase::kFollowRefit, elements,
                    TrimAmbiguous<Det, Material, Rule, Shape>{
                      .search = search,
                      .parameters = in.parameters,
                      .frame = in.frame,
                      .geometry = in.geometry,
                      .material = in.material,
                      .propagator = in.propagator,
                      .tracks = in.tracks,
                      .proposals = in.proposals,
                      .trackOffset = in.trackOffset,
                    });
      });
    },
              choices.shapeTest);
  }
  launch.pass(Phase::kFollowMark, elements, Mark<kMaxClusters>{in.proposals, in.used});
}

} // namespace follow

/// Host orchestration and bounded storage of the follow: propose, then bid,
/// resolve, combine, refit, trim the
/// ambiguous last rows (`followTrimMargin`) and mark. The
/// tracker runs it on an iteration's own tracks and, with `followEarlier`,
/// on the earlier iterations' tracks first; a device runs the same sequence
/// (`follow::runSequence`) with its own launches.
template <DetectorLayout Det>
class Follow
{
 public:
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  static constexpr std::size_t kLayerCount = static_cast<std::size_t>(Detector<Det>::nLayers);
  using Proposal = FollowProposal<kMaxClusters>;
  using Tracks = typename TrackRefit<Det>::Tracks;

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool = std::move(pool);
    clearStorage();
  }

  /// The arrays back to the pool.
  void clearStorage() noexcept { deepVectorClear(mMemoryPool.get(), mProposals, mBids, mRows, mCounts, mEndStarts); }

  /// Follow the tracks from \a trackOffset on, inward too with \a inward,
  /// each layer's candidates bound by \a maxCandidateWidthPS, the frame's
  /// widest time box there, the sensors
  /// near a crossing from \a cells. \a earlier: the earlier iterations'
  /// tracks (`followEarlier`), counted and timed apart -- within the
  /// caller's `kFollowEarlier` zone, no phase of its own.
  template <typename Material, typename ForEach>
  [[nodiscard]] StageStatus run(ForEach forEach, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator, Tracks& tracks, std::uint32_t trackOffset, bool inward, bool earlier, std::span<std::uint8_t> used, std::span<const std::uint32_t, kLayerCount> maxCandidateWidthPS, SensorCellsView cells)
  {
    using profile::Phase;
    return guardOutOfMemory([&] {
      const auto trackCount = static_cast<std::uint32_t>(tracks.size());
      const auto elements = trackCount - trackOffset;
      const bool end = parameters.followEnd && !earlier;
      {
        const profile::Zone zone{earlier ? Phase::kNone : Phase::kFollowPropose};
        // Kept from iteration to iteration and filled by the passes, not
        // reallocated.
        resizeBoundedVector(mProposals, elements, mMemoryPool.get());
        resizeBoundedVector(mBids, used.size(), mMemoryPool.get());
        resizeBoundedVector(mRows, elements, mMemoryPool.get());
        resizeBoundedVector(mCounts, elements, mMemoryPool.get());
        if (end) {
          resizeBoundedVector(mEndStarts, elements, mMemoryPool.get());
        }
      }
      follow::runSequence(follow::HostLaunch<ForEach>{.forEach = &forEach, .earlier = earlier}, follow::Inputs<Det, Material>{
                                                                                                  .parameters = &parameters,
                                                                                                  .choices = &parameters,
                                                                                                  .frame = &frame,
                                                                                                  .geometry = &geometry,
                                                                                                  .material = &material,
                                                                                                  .propagator = propagator,
                                                                                                  .tracks = tracks.view(),
                                                                                                  .constTracks = static_cast<const Tracks&>(tracks).view(),
                                                                                                  .trackOffset = trackOffset,
                                                                                                  .trackCount = trackCount,
                                                                                                  .inward = inward,
                                                                                                  .earlier = earlier,
                                                                                                  .used = used.data(),
                                                                                                  .rowCount = used.size(),
                                                                                                  .maxCandidateWidthPS = maxCandidateWidthPS.data(),
                                                                                                  .cells = cells,
                                                                                                  .proposals = mProposals.data(),
                                                                                                  .bids = mBids.data(),
                                                                                                  .rows = mRows.data(),
                                                                                                  .counts = mCounts.data(),
                                                                                                  .endStarts = end ? mEndStarts.data() : nullptr,
                                                                                                });
      return StageStatus::kComplete;
    });
  }

  /// Each track's proposal after the passes: the rows it took.
  [[nodiscard]] std::span<const Proposal> proposals() const { return mProposals; }

  /// Whether an iteration's follow decides its ends on the last row's
  /// density (`FollowStage`'s EndDensity): either of its keys is set.
  [[nodiscard]] static bool endDensity(const IterationParameters<Det>& parameters) { return parameters.followEndMaxSqrtDet > 0.F || parameters.followKeepDensity > 0.F; }

  /// Proposal \a first + i taken from proposal kept[i], for every i, as the
  /// tracker drops short tracks in parallel: \a kept ascends
  /// and kept[i] >= first + i, so the moves go through a scratch.
  /// \a forEach(count, body) runs body(i) for i in [0, count).
  template <typename Index, typename ForEach>
  void keepProposals(std::span<const Index> kept, std::size_t first, const ForEach& forEach)
  {
    const auto count = static_cast<Index>(kept.size());
    const auto move = [&](auto& values) {
      using Value = typename std::remove_reference_t<decltype(values)>::value_type;
      BoundedVector<Value> scratch(std::pmr::polymorphic_allocator<Value>{mMemoryPool.get()});
      resizeBoundedVector(scratch, count);
      forEach(count, [&](Index at) { scratch[at] = values[kept[at]]; });
      forEach(count, [&](Index at) { values[first + at] = scratch[at]; });
    };
    move(mProposals);
  }
  void truncateProposals(std::size_t count) { mProposals.resize(count); }

 private:
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  BoundedVector<Proposal> mProposals;
  BoundedVector<std::uint32_t> mBids;
  BoundedVector<ClusterList<kMaxClusters>> mRows;
  BoundedVector<std::uint8_t> mCounts;
  BoundedVector<BoundTrackState> mEndStarts; ///< with followEnd: each track's state without its outermost row
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_FOLLOW_H_
