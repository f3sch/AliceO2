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
#include "GenericTracking/Base/Statistics.h"
#include "GenericTracking/Tracking/BeamConstraint.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/TrackFitterUtils.h"

#include "BandMaterial.h"
#include "ScatteringKick.h"
#include "ToyBarrel8.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

namespace o2::gtrk::test
{

/// An electron chain's material: each layer's band holds \a x2X0 and
/// \a xRho [g/cm^2] at normal incidence (`bandMaterial`), nothing between.
struct ChainBands {
  Scalar x2X0 = 0.F;
  Scalar xRho = 0.F;
};

/// An electron through ToyBarrel8's eight layers (2-16 cm) on its bands
/// (`ChainBands`, the chain's link material), as the filter models it
/// (LinkMaterial::step): at each sensor the band's part behind the
/// measurement on leaving it and the part in front on arriving, the band
/// split at the measurement's depth, at the step chord's incidence; each part
/// a kick of its Highland theta0 (the logarithm at the span's material) and
/// its radiation, each fraction z drawn from Bethe and Heitler's model itself
/// (-ln z Gamma-distributed with shape t / ln 2), not from the mixture; the
/// ionisation of both parts on arriving (its mean and the straggling the
/// filter takes, Gaussian in q/p). truth[i] is the state on plane i after the
/// part in front of it: what the outward pass ends with on the last plane and
/// the inward one on the first. The chords are those of the track without
/// material. From a vertex (`generate`'s \a fromVertex), the track starts at
/// (`vertexX`, `vertexY`, 0), on the perigee plane about it (d0 = z0 = 0),
/// and reaches the first layer through nothing but its band's part in front,
/// taken on arriving there, before its measurement (what the refit's carry
/// to the vertex undoes); `vertex` holds that start.
struct ElectronChain {
  static constexpr std::size_t kLayers = 8;
  using Barrel = detectors::ToyBarrel8;
  Geometry<Barrel> geometry = toyGeometry<Barrel>();
  /// On the heap: with its link material it is megabytes. Its link material
  /// is the chain's bands' (`generate`).
  std::unique_ptr<TrackingGeometry<Barrel>> tracking = std::make_unique<TrackingGeometry<Barrel>>(geometry);
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
  ChainBands held{}; ///< the bands the link material holds

  ElectronChain()
  {
    propagator.massGeV = 0.000511F;
    tracking->material() = noMaterial<Barrel>();
  }

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
    const auto& from = tracking->surface(state.surface);
    for (SensorID id = firstSensorOfLayer<Barrel>(layer); id < firstSensorOfLayer<Barrel>(static_cast<LayerID>(layer + 1)); ++id) {
      BoundTrackState probe = state;
      if (propagate(propagator, probe, from, tracking->surface(id), PropagationDirection::kForward).status == PropagationStatus::kOk) {
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

  /// A track of \a qOverP through \a bands; false if it does not reach
  /// every layer.
  bool generate(Scalar qOverP, const ChainBands& bands, bool radiate, std::mt19937& rng, bool fromVertex = false)
  {
    if (bands.x2X0 != held.x2X0 || bands.xRho != held.xRho) {
      std::array<double, kLayers> x{};
      std::array<double, kLayers> rho{};
      x.fill(bands.x2X0);
      rho.fill(bands.xRho);
      tracking->material() = *bandMaterial<Barrel>(x, {}, {.low = -20., .high = 20.}, rho);
      held = bands;
    }
    BoundTrackState start{};
    start.par[kPhi] = 0.05F;
    start.par[kTheta] = 1.55F;
    start.par[kQOverP] = qOverP;
    start.surface = firstSensorOfLayer<Barrel>(0);
    std::normal_distribution<double> gauss;
    before.fill(1.);
    after.fill(1.);
    // Each band's part in front of (\a front) or behind a point on it, the
    // chord \a d's path there: x/X0 and rho x.
    const auto part = [&](LayerID layer, const Scalar* point, const Scalar* d, bool front) {
      const auto [share, path] = bandShare<Barrel>(layer, point, d);
      const Scalar k = (front ? share : 1.F - share) * path;
      return std::pair{bands.x2X0 * k, bands.xRho * k};
    };
    // The ionisation of \a xRho with the straggling the filter takes.
    const auto lose = [&](BoundTrackState& state, Scalar xRho) {
      Scalar variance = 0.F;
      Scalar lossQOverP = state.qOverP();
      if (!detail::applyEnergyLoss(propagator, xRho, PropagationDirection::kForward, lossQOverP, state, variance)) {
        return false;
      }
      state.par[kQOverP] += static_cast<Scalar>(gauss(rng)) * std::sqrt(variance);
      return true;
    };
    if (fromVertex) {
      // The perigee plane about the origin, normal to the direction there; the
      // track through the origin; the first layer reached, its band in front.
      vertex = BoundTrackState{};
      vertex.par[kPhi] = start.par[kPhi];
      vertex.par[kTheta] = start.par[kTheta];
      vertex.par[kQOverP] = qOverP;
      vertexPlane = beamLinePlane(vertexX, vertexY, vertex.phi());
      vertex.surface = vertexPlane.id;
      bool crossed = false;
      for (SensorID id = firstSensorOfLayer<Barrel>(0); id < firstSensorOfLayer<Barrel>(1) && !crossed; ++id) {
        start = vertex;
        crossed = propagate(propagator, start, vertexPlane, tracking->surface(id), PropagationDirection::kForward).status == PropagationStatus::kOk;
      }
      if (!crossed) {
        return false;
      }
      const Scalar a[3] = {vertexX, vertexY, 0.F};
      Scalar b[3]{};
      surfaceToGlobal(tracking->surface(start.surface), start.loc0(), start.loc1(), b);
      const Scalar d[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
      const auto [x, xRho] = part(0, b, d, true);
      if (!lose(start, xRho)) {
        return false;
      }
      scatteringKickKickedTheta(start, highlandTheta0(x, x, propagator.momentum(qOverP), propagator.massGeV, propagator.absCharge), gauss, rng);
      if (radiate) {
        before[0] = fraction(static_cast<double>(x), rng);
        start.par[kQOverP] /= static_cast<Scalar>(before[0]);
      }
    }
    // Each step's band parts, on the chords without material: what the fit's
    // span holds.
    std::array<Scalar, kLayers> leaving{};
    std::array<Scalar, kLayers> leavingRho{};
    std::array<Scalar, kLayers> reaching{};
    std::array<Scalar, kLayers> reachingRho{};
    BoundTrackState state = start;
    Scalar previous[3];
    surfaceToGlobal(tracking->surface(state.surface), state.loc0(), state.loc1(), previous);
    Scalar span = 0.F;
    for (std::size_t layer = 1; layer < kLayers; ++layer) {
      if (!cross(state, static_cast<LayerID>(layer))) {
        return false;
      }
      Scalar point[3];
      surfaceToGlobal(tracking->surface(state.surface), state.loc0(), state.loc1(), point);
      const Scalar d[3] = {point[0] - previous[0], point[1] - previous[1], point[2] - previous[2]};
      std::tie(leaving[layer - 1], leavingRho[layer - 1]) = part(static_cast<LayerID>(layer - 1), previous, d, false);
      std::tie(reaching[layer], reachingRho[layer]) = part(static_cast<LayerID>(layer), point, d, true);
      span += leaving[layer - 1] + reaching[layer];
      std::copy(point, point + 3, previous);
    }
    state = start;
    for (std::size_t layer = 0; layer < kLayers; ++layer) {
      if (layer > 0) {
        if (!cross(state, static_cast<LayerID>(layer))) {
          return false;
        }
        if (!lose(state, leavingRho[layer - 1] + reachingRho[layer])) {
          return false;
        }
        scatteringKick(state, highlandTheta0(reaching[layer], span, propagator.momentum(state.qOverP()), propagator.massGeV, propagator.absCharge), gauss, rng);
        if (radiate) {
          before[layer] = fraction(static_cast<double>(reaching[layer]), rng);
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
      if (layer + 1 < kLayers) {
        scatteringKick(state, highlandTheta0(leaving[layer], span, propagator.momentum(state.qOverP()), propagator.massGeV, propagator.absCharge), gauss, rng);
        if (radiate) {
          after[layer] = fraction(static_cast<double>(leaving[layer]), rng);
          state.par[kQOverP] /= static_cast<Scalar>(after[layer]);
        }
      }
      if (!isFiniteOf(state.qOverP()) || std::abs(state.qOverP()) > 50.F) {
        return false;
      }
    }
    return true;
  }
};

/// The median (quantile 0.5) and the half central 68 % width (halfWidth68)
/// of \a values.
inline std::pair<double, double> medianWidth(std::vector<double> values)
{
  const double median = quantile(values, 0.5);
  return {median, halfWidth68(values)};
}

} // namespace o2::gtrk::test

#endif // GTRK_TEST_ELECTRONCHAIN_H_
