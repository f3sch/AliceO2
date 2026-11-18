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

/// \file ElectronChain.h
/// \brief Electrons through ToyBarrel8 as the Gaussian-sum filter models them, Bethe and Heitler's losses sampled.

#ifndef GTRK_TEST_ELECTRONCHAIN_H_
#define GTRK_TEST_ELECTRONCHAIN_H_

#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Tracking/BeamConstraint.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/TrackFitterUtils.h"
#include "ToyBarrel8.h"
#include "UniformMaterial.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <random>
#include <utility>
#include <vector>

namespace o2::gtrk::test
{

/// An electron through ToyBarrel8's eight layers (2-16 cm) on uniform
/// material, as the filter models it: on each step the scattering in two
/// kicks of half the chord's Highland theta0 (behind the sensor left, in
/// front of the sensor reached, the logarithm at the span's material), the
/// ionisation over the chord on reaching the next sensor (its mean and the
/// straggling the filter takes, Gaussian in q/p), and at each
/// plane the radiation of the half-gap before it, then of the half-gap after
/// it, each fraction z drawn from Bethe and Heitler's model itself (-ln z
/// Gamma-distributed with shape t / ln 2), not from the mixture. truth[i] is
/// the state on plane i after the half-gap before it: what the outward pass
/// ends with on the last plane and the inward one on the first. The chords
/// are those of the track without material. From a vertex (`generate`'s
/// \a fromVertex), the track starts at (`vertexX`, `vertexY`, 0), on the
/// perigee plane about it (d0 = z0 = 0), and reaches the first layer through the
/// segment's material as through a step's, the segment's radiation taken on
/// arriving there, before its measurement (the material inside the chain's
/// first plane, which the refit's extrapolation to the vertex undoes);
/// `vertex` holds that start.
struct ElectronChain {
  static constexpr std::size_t kLayers = 8;
  using Barrel = detectors::ToyBarrel8;
  Geometry<Barrel> geometry = toyGeometry<Barrel>();
  TrackingGeometry<Barrel> tracking{geometry};
  Propagator propagator;
  std::array<BoundTrackState, kLayers> truth;
  std::array<Measurement, kLayers> measurements;
  BoundTrackState vertex{};
  TrackingSurface vertexPlane{};
  Scalar vertexX = 0.F; ///< the vertex's transverse position, its z 0
  Scalar vertexY = 0.F;
  /// Each plane's retained fractions: of the half-gap before it, after it
  /// (before the first: the segment from the vertex).
  std::array<double, kLayers> before{};
  std::array<double, kLayers> after{};
  Scalar sigma = 5e-4F;

  ElectronChain() { propagator.massGeV = 0.000511F; }

  /// The share of the momentum radiated between rows \a first and \a last.
  [[nodiscard]] double radiatedBetween(std::size_t first, std::size_t last) const
  {
    double kept = 1.;
    for (std::size_t i = first; i < last; ++i) {
      kept *= after[i] * before[i + 1];
    }
    return 1. - kept;
  }

  [[nodiscard]] bool cross(BoundTrackState& state, LayerID layer) const
  {
    const auto& from = tracking.surface(state.surface);
    for (SensorID id = firstSensorOfLayer<Barrel>(layer); id < firstSensorOfLayer<Barrel>(static_cast<LayerID>(layer + 1)); ++id) {
      BoundTrackState probe = state;
      if (propagate(propagator, probe, from, tracking.surface(id), PropagationDirection::kForward).status == PropagationStatus::kOk) {
        state = probe;
        return true;
      }
    }
    return false;
  }

  static Scalar fraction(double x2X0, std::mt19937& rng)
  {
    if (!(x2X0 > 0.)) {
      return 1.F;
    }
    std::gamma_distribution<double> gamma(x2X0 / std::numbers::ln2, 1.);
    return static_cast<Scalar>(std::exp(-gamma(rng)));
  }

  /// A track of \a qOverP; false if it does not reach every layer.
  bool generate(Scalar qOverP, const UniformMaterial& material, bool radiate, std::mt19937& rng, bool fromVertex = false)
  {
    BoundTrackState start{};
    start.par[kPhi] = 0.05F;
    start.par[kTheta] = 1.55F;
    start.par[kQOverP] = qOverP;
    start.surface = firstSensorOfLayer<Barrel>(0);
    std::normal_distribution<double> gauss;
    before.fill(1.);
    after.fill(1.);
    if (fromVertex) {
      // The perigee plane about the origin, normal to the direction there; the
      // track through the origin; the first layer reached through the segment.
      vertex = BoundTrackState{};
      vertex.par[kPhi] = start.par[kPhi];
      vertex.par[kTheta] = start.par[kTheta];
      vertex.par[kQOverP] = qOverP;
      vertexPlane = beamLinePlane(vertexX, vertexY, vertex.phi());
      vertex.surface = vertexPlane.id;
      BoundTrackState probe = vertex;
      bool crossed = false;
      for (SensorID id = firstSensorOfLayer<Barrel>(0); id < firstSensorOfLayer<Barrel>(1) && !crossed; ++id) {
        probe = vertex;
        crossed = propagate(propagator, probe, vertexPlane, tracking.surface(id), PropagationDirection::kForward).status == PropagationStatus::kOk;
      }
      if (!crossed) {
        return false;
      }
      const Scalar a[3] = {vertexX, vertexY, 0.F};
      Scalar b[3]{};
      surfaceToGlobal(tracking.surface(probe.surface), probe.loc0(), probe.loc1(), b);
      const Scalar segment = std::hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
      const Scalar half = 0.5F * material.x2X0PerCm * segment;
      const Scalar theta0 = highlandTheta0(half, half, 1.F / std::abs(qOverP), propagator.massGeV);
      BoundTrackState leaving = vertex;
      leaving.par[kTheta] += static_cast<Scalar>(gauss(rng)) * theta0;
      leaving.par[kPhi] += static_cast<Scalar>(gauss(rng)) * theta0 / std::sin(leaving.theta());
      crossed = false;
      for (SensorID id = firstSensorOfLayer<Barrel>(0); id < firstSensorOfLayer<Barrel>(1) && !crossed; ++id) {
        start = leaving;
        crossed = propagate(propagator, start, vertexPlane, tracking.surface(id), PropagationDirection::kForward).status == PropagationStatus::kOk;
      }
      if (!crossed) {
        return false;
      }
      Scalar variance = 0.F;
      Scalar lossQOverP = start.qOverP();
      if (!detail::applyEnergyLoss(propagator, material.rho * segment, PropagationDirection::kForward, lossQOverP, start, variance)) {
        return false;
      }
      start.par[kQOverP] += static_cast<Scalar>(gauss(rng)) * std::sqrt(variance);
      start.par[kTheta] += static_cast<Scalar>(gauss(rng)) * theta0;
      start.par[kPhi] += static_cast<Scalar>(gauss(rng)) * theta0 / std::sin(start.theta());
      if (radiate) {
        before[0] = fraction(2. * static_cast<double>(half), rng);
        start.par[kQOverP] /= static_cast<Scalar>(before[0]);
      }
    }
    // The chords without material.
    std::array<Scalar, kLayers> chord{};
    BoundTrackState state = start;
    Scalar previous[3];
    surfaceToGlobal(tracking.surface(state.surface), state.loc0(), state.loc1(), previous);
    for (std::size_t layer = 1; layer < kLayers; ++layer) {
      if (!cross(state, static_cast<LayerID>(layer))) {
        return false;
      }
      Scalar point[3];
      surfaceToGlobal(tracking.surface(state.surface), state.loc0(), state.loc1(), point);
      chord[layer] = std::hypot(point[0] - previous[0], point[1] - previous[1], point[2] - previous[2]);
      std::copy(point, point + 3, previous);
    }
    Scalar span = 0.F;
    for (const Scalar c : chord) {
      span += material.x2X0PerCm * c;
    }
    const auto kick = [&](BoundTrackState& s, Scalar theta0) {
      const auto sinTheta = std::sin(s.theta());
      s.par[kTheta] += static_cast<Scalar>(gauss(rng)) * theta0;
      s.par[kPhi] += static_cast<Scalar>(gauss(rng)) * theta0 / sinTheta;
    };
    const auto halfGap = [&](std::size_t layer) { return 0.5 * static_cast<double>(material.x2X0PerCm * chord[layer]); };
    state = start;
    for (std::size_t layer = 0; layer < kLayers; ++layer) {
      if (layer > 0) {
        const Scalar theta0 = highlandTheta0(0.5F * material.x2X0PerCm * chord[layer], span, 1.F / std::abs(state.qOverP()), propagator.massGeV);
        kick(state, theta0);
        if (!cross(state, static_cast<LayerID>(layer))) {
          return false;
        }
        Scalar variance = 0.F;
        Scalar lossQOverP = state.qOverP();
        if (!detail::applyEnergyLoss(propagator, material.rho * chord[layer], PropagationDirection::kForward, lossQOverP, state, variance)) {
          return false;
        }
        state.par[kQOverP] += static_cast<Scalar>(gauss(rng)) * std::sqrt(variance); // the straggling the filter takes
        kick(state, theta0);
        if (radiate) {
          before[layer] = fraction(halfGap(layer), rng);
          state.par[kQOverP] /= static_cast<Scalar>(before[layer]);
        }
      }
      truth[layer] = state;
      auto& m = measurements[layer];
      m = Measurement{};
      m.surface = state.surface;
      m.loc0 = state.loc0() + static_cast<Scalar>(gauss(rng)) * sigma;
      m.loc1 = state.loc1() + static_cast<Scalar>(gauss(rng)) * sigma;
      m.varLoc0 = sigma * sigma;
      m.varLoc1 = sigma * sigma;
      m.clusterIndex = static_cast<std::uint32_t>(layer);
      if (radiate && layer + 1 < kLayers) {
        after[layer] = fraction(halfGap(layer + 1), rng);
        state.par[kQOverP] /= static_cast<Scalar>(after[layer]);
      }
      if (!isFiniteOf(state.qOverP()) || std::abs(state.qOverP()) > 50.F) {
        return false;
      }
    }
    return true;
  }
};

/// The median and half the 16-84 % width of \a values.
inline std::pair<double, double> medianWidth(std::vector<double> values)
{
  std::ranges::sort(values);
  const auto at = [&](double q) { return values[static_cast<std::size_t>(q * static_cast<double>(values.size() - 1))]; };
  return {at(0.5), 0.5 * (at(0.84) - at(0.16))};
}

} // namespace o2::gtrk::test

#endif // GTRK_TEST_ELECTRONCHAIN_H_
