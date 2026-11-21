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

/// \file testInactiveLayers.cxx
/// \brief Inactive layers (LayerSpec.h), as ALICE 3's TOF in TRKFT3TOF: never
///        a hole, found by the follow between two rows of a road that
///        straddles one (the interior row) and after the outward walk (the
///        tail).
///
/// ToyInactive8 is ToyBarrel8's eight layers with layers 3 and 7 inactive:
/// layer 3 sits between 2 and 4, which also link directly (2 -> 4); layer 7
/// is the last, as ALICE 3's oTOF.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK inactive layers
#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/KFTrackFitter.h"
#include "GenericTracking/Tracking/LinkHoles.h"
#include "GenericTracking/Tracking/MeasurementGridBuilder.h"
#include "GenericTracking/Tracking/SensorCells.h"
#include "GenericTracking/Tracking/SlicePlan.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include "AnalyticLayout.h"
#include "BandMaterial.h"
#include "GridTestSupport.h"
#include "SerialReferences.h"
#include "ToyBarrel8.h"
#include "UniformMaterial.h"

#include "CommonConstants/PhysicsConstants.h"

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <vector>

namespace o2::gtrk
{

namespace detectors
{

/// ToyBarrel8 with layers 3 and 7 inactive (see above).
struct ToyInactive8 {
  static constexpr const char* name = "ToyInactive8";
  static constexpr auto barrelSpecs = ToyBarrel8::barrelSpecs;
  static constexpr std::array<DiscLayerSpec, 0> discSpecs{};
  static constexpr std::uint64_t inactiveLayers = (std::uint64_t{1} << 3U) | (std::uint64_t{1} << 7U);
  static constexpr std::array<LayerLink, 8> layerLinks{{
    {.from = 0, .to = 1, .sourceMin = -1.50, .sourceMax = 1.50, .skipsLayer = false},
    {.from = 1, .to = 2, .sourceMin = -1.50, .sourceMax = 1.50, .skipsLayer = false},
    {.from = 2, .to = 3, .sourceMin = -1.50, .sourceMax = 1.50, .skipsLayer = false},
    {.from = 2, .to = 4, .sourceMin = -1.50, .sourceMax = 1.50, .skipsLayer = true},
    {.from = 3, .to = 4, .sourceMin = -1.50, .sourceMax = 1.50, .skipsLayer = false},
    {.from = 4, .to = 5, .sourceMin = -1.50, .sourceMax = 1.50, .skipsLayer = false},
    {.from = 5, .to = 6, .sourceMin = -1.50, .sourceMax = 1.50, .skipsLayer = false},
    {.from = 6, .to = 7, .sourceMin = -1.50, .sourceMax = 1.50, .skipsLayer = false},
  }};
  static constexpr auto barrelLayers = placedLayers(barrelSpecs);
  static constexpr auto discLayers = placedLayers(discSpecs);
  static constexpr int nSensors = countSensors(barrelLayers, discLayers);
};

inline const bool kToyInactive8Published = publishToy<ToyInactive8>();

} // namespace detectors

namespace
{

using Det = detectors::ToyInactive8;
using Frame = TimeFrame<Det>;
constexpr std::size_t kLayers = 8U;
constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
using Tracks = BoundedColumns<TrackOf<kMaxClusters>::template RecordT>;
using Proposal = FollowProposal<kMaxClusters>;

const TrackingGeometry<Det>& trackingGeometry()
{
  static const auto detector = toyGeometry<Det>();
  static const TrackingGeometry<Det> tracking{detector};
  return tracking;
}

const SensorCells<Det>& cells()
{
  static const SensorCells<Det> sensorCells(trackingGeometry());
  return sensorCells;
}

struct Hit {
  LayerID layer = kInvalidLayerID;
  SensorID surface = kInvalidSensorID;
  Scalar x = 0.F;
  Scalar y = 0.F;
  Scalar z = 0.F;
  Scalar loc0 = 0.F;
  Scalar loc1 = 0.F;
  TimeEstPS time = test::toyTimeBox(test::kToyTrackTimeNS);
};

/// A helix from the origin at 2 T on every layer: the first sensor whose
/// plane it meets inside its bounds (as ToyBarrel8Fixture.h's `sensorHelixHits`).
std::vector<Hit> helixHits(double phi0, double slope, double ptGeV)
{
  const auto& tracking = trackingGeometry();
  const double curvature = 1. / curvatureRadius(ptGeV, 20.);
  const auto at = [&](double arc, double position[3]) {
    const double direction = phi0 + (curvature * arc);
    position[0] = (std::sin(direction) - std::sin(phi0)) / curvature;
    position[1] = (std::cos(phi0) - std::cos(direction)) / curvature;
    position[2] = slope * arc;
  };
  std::vector<Hit> hits;
  for (std::size_t layer = 0; layer < kLayers; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    const double crossing = 2. * std::asin(0.5 * curvature * layerDistance<Det>(id)) / curvature;
    for (SensorID sensor = firstSensorOfLayer<Det>(id); sensor < firstSensorOfLayer<Det>(static_cast<LayerID>(layer + 1)); ++sensor) {
      const auto& surface = tracking.surface(sensor);
      const auto height = [&](double arc) {
        double position[3]{};
        at(arc, position);
        double result = 0.;
        for (int i = 0; i < 3; ++i) {
          result += (position[i] - static_cast<double>(surface.origin[i])) * static_cast<double>(surface.normal[i]);
        }
        return result;
      };
      double low = 0.8 * crossing;
      double high = 1.3 * crossing;
      if ((height(low) < 0.) == (height(high) < 0.)) {
        continue;
      }
      for (int i = 0; i < 100; ++i) {
        const double middle = 0.5 * (low + high);
        ((height(middle) < 0.) == (height(low) < 0.) ? low : high) = middle;
      }
      double position[3]{};
      at(0.5 * (low + high), position);
      const Scalar global[3]{static_cast<Scalar>(position[0]), static_cast<Scalar>(position[1]), static_cast<Scalar>(position[2])};
      Scalar local[3]{};
      surfaceToLocal(surface, global, local);
      if (!surface.bounds.isInside(local[0], local[1])) {
        continue;
      }
      hits.push_back({.layer = id, .surface = sensor, .x = global[0], .y = global[1], .z = global[2], .loc0 = local[0], .loc1 = local[1]});
      break;
    }
  }
  return hits;
}

SearchWindowParams windowParameters()
{
  SearchWindowParams result;
  result.minPtGeV = 1.;
  result.bzKiloGauss = 20.;
  result.materialX2X0 = 0.;
  result.positionSigma = 0.003;
  result.nSigma = 3.;
  return result;
}

/// One timeframe of \a hits, a measurement each: its cluster index its hit's.
struct Fixture {
  explicit Fixture(const std::vector<Hit>& hits)
    : geometry(test::cappedGrid<Det>(20., test::withToyFrames(GridParams{.blockWidth = 20'000}), windowParameters(), VertexConstraint{.z = 0., .sigmaZ = 0.02, .sigmaTransverse = 0.001}, 16, 8)),
      offsets(static_cast<std::size_t>(geometry.bucketCount()) + 1U),
      frame(test::frameOf(
        geometry, hits.size(),
        [&](std::size_t i) {
          const auto& hit = hits[i];
          double phi = std::atan2(hit.y, hit.x);
          if (phi < 0.) {
            phi += kTwoPi<double>;
          }
          return toybarrel8InputPoint(hit, static_cast<Scalar>(phi));
        },
        [&](std::size_t i) { return hits[i].time; },
        [](std::size_t /*i*/) { return test::PointVariances{.z = 1.e-6F, .r = 2.e-6F, .rPhi = 1.e-6F}; },
        columns, offsets,
        [&](std::size_t i, std::uint32_t row) {
          columns.surface[row] = hits[i].surface;
          columns.loc0[row] = hits[i].loc0;
          columns.loc1[row] = hits[i].loc1;
          columns.varLoc0[row] = 1.e-6F;
          columns.varLoc1[row] = 1.e-6F;
        }))
  {
  }

  struct InputPoint {
    LayerID layer = kInvalidLayerID;
    Scalar x = 0.F;
    Scalar y = 0.F;
    Scalar z = 0.F;
    Scalar r = 0.F;
    Scalar phi = 0.F;
  };
  static InputPoint toybarrel8InputPoint(const Hit& hit, Scalar phi) { return {.layer = hit.layer, .x = hit.x, .y = hit.y, .z = hit.z, .r = std::hypot(hit.x, hit.y), .phi = phi}; }

  [[nodiscard]] std::uint32_t row(std::size_t layer) const { return static_cast<std::uint32_t>(std::ranges::find(columns.clusterIndex, static_cast<std::uint32_t>(layer)) - columns.clusterIndex.begin()); }

  GridGeometry<Det> geometry;
  std::vector<std::uint32_t> offsets;
  MeasurementStore columns;
  std::unique_ptr<Frame> frame;
};

/// The track of \a fixture's rows on \a layers, innermost first, fitted.
Tracks fittedTrack(const Fixture& fixture, std::initializer_list<std::size_t> layers, const std::shared_ptr<BoundedMemoryResource>& pool)
{
  ClusterList<kMaxClusters> rows;
  std::ranges::fill(rows.index, ClusterList<kMaxClusters>::kNoCluster);
  std::uint8_t count = 0;
  for (const auto layer : layers) {
    rows.index[count++] = fixture.row(layer);
  }
  const auto fit = fitMeasurements(Propagator{}, trackingGeometry(), UniformMaterial{}, RowChain<Det>{.frame = fixture.frame.get(), .rows = rows.index}, count);
  BOOST_REQUIRE(fit.status == PropagationStatus::kOk);
  Tracks tracks;
  tracks.setMemoryResource(pool.get());
  tracks.resize(1U);
  auto view = tracks.view();
  view.state[0] = fit.state;
  view.outwardState[0] = fit.outwardState;
  view.chi2[0] = fit.chi2;
  view.ndf[0] = fit.ndf;
  view.status[0] = fit.status;
  view.clusters[0] = rows;
  view.count[0] = count;
  view.iteration[0] = 0U;
  view.addedOutward[0] = 0U;
  view.addedInward[0] = 0U;
  return tracks;
}

IterationParameters<Det> followParameters()
{
  IterationSettings settings{.windows = windowParameters(), .vertexSource = VertexSource::kBeamLine};
  settings.follow = true;
  return iterationParameters<Det>(settings, VertexConstraint{}, noMaterial<Det>());
}

Proposal followTrack(const Fixture& fixture, const Tracks& tracks, const std::vector<std::uint8_t>& used, const IterationParameters<Det>& parameters)
{
  const auto widths = frameBounds(*fixture.frame).widthPS;
  std::vector<Proposal> proposals(tracks.size());
  const UniformMaterial material;
  const FollowStage<Det, UniformMaterial, false, false, false, false, false> stage{&parameters, fixture.frame.get(), &trackingGeometry(), &material, Propagator{}, tracks.view(), used.data(), widths.data(), cells().view(), proposals.data(), true, 0U, nullptr};
  stage(0U);
  return proposals.front();
}

} // namespace

/// Both the interior row and the tail use the iteration's cluster gate.
BOOST_AUTO_TEST_CASE(InactiveLayersUseTheIterationsClusterGate)
{
  static_assert(Detector<Det>::isInactiveLayer(3) && Detector<Det>::isInactiveLayer(7) && !Detector<Det>::isInactiveLayer(2));
  static_assert(Detector<detectors::ToyBarrel8>::inactiveLayers == 0U);
  auto hits = helixHits(0.42, 0.025, 1.8);
  BOOST_REQUIRE(hits.size() == kLayers);
  for (const std::size_t layer : {3U, 7U}) {
    hits[layer].z += 0.002F;
    hits[layer].loc1 += 0.002F;
  }
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, {0, 1, 2, 4, 5}, pool);
  const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  auto parameters = followParameters();

  const auto accepted = followTrack(fixture, tracks, used, parameters);
  BOOST_TEST(accepted.interior == fixture.row(3));
  BOOST_REQUIRE_EQUAL(static_cast<int>(accepted.outwardCount), 2);
  BOOST_TEST(accepted.outward[1] == fixture.row(7));

  parameters.clusterChi2Max = 0.1F;
  const auto rejected = followTrack(fixture, tracks, used, parameters);
  BOOST_TEST(rejected.interior == ClusterList<kMaxClusters>::kNoCluster);
  BOOST_REQUIRE_EQUAL(static_cast<int>(rejected.outwardCount), 1);
  BOOST_TEST(rejected.outward[0] == fixture.row(6));
}

/// A pair crossing an inactive layer has no hole there: 2 -> 4 skips
/// layer 3 and crosses nothing; in ToyBarrel8's numbering a 2 -> 4 pair
/// would cross layer 3.
BOOST_AUTO_TEST_CASE(AnInactiveLayerIsNeverAHole)
{
  BOOST_TEST(linkHoleLayers<Det>(2, 4) == 0U);
  BOOST_TEST(linkHoleLayers<detectors::ToyBarrel8>(2, 4) == (std::uint64_t{1} << 3U));
  BOOST_TEST(linkHoleLayers<Det>(1, 4) == (std::uint64_t{1} << 2U));
}

/// A step over the inactive layer crosses its material: the follow's step on
/// the direct link 2 -> 4 holds what the steps 2 -> 3 and 3 -> 4 through it
/// hold, outward and inward, as ALICE 3's ML6 -> ML7 and ML6 -> iTOF -> ML7.
BOOST_AUTO_TEST_CASE(AStepOverTheInactiveLayerCrossesItsMaterial)
{
  static_cast<void>(trackingGeometry()); // the layout's envelopes
  std::array<double, kLayers> x2X0{};
  x2X0.fill(0.005);
  x2X0[3] = 0.02;
  const auto material = bandMaterial<Det>(x2X0);
  // A straight line from the beam line through the layers' nominal surfaces.
  constexpr double kPhi = 0.3;
  constexpr double kSlope = 0.1;
  const auto at = [](LayerID layer) {
    const double r = PlacedEnvelopes<Det>::get(layer).position;
    return std::array{static_cast<Scalar>(r * std::cos(kPhi)), static_cast<Scalar>(r * std::sin(kPhi)), static_cast<Scalar>(0.2 + (kSlope * r))};
  };
  const auto total = [&](LayerID from, LayerID to) {
    const auto step = material->step(from, to, at(from).data(), at(to).data());
    BOOST_TEST(step.found);
    return static_cast<double>(step.before + step.after);
  };
  const double direct = total(2, 4);
  const double through = total(2, 3) + total(3, 4);
  BOOST_TEST(direct == through, boost::test_tools::tolerance(1.e-2));
  BOOST_TEST(total(4, 2) == direct, boost::test_tools::tolerance(1.e-2));
  BOOST_TEST(total(4, 3) + total(3, 2) == through, boost::test_tools::tolerance(1.e-2));
  // The inactive layer's band is most of it.
  BOOST_TEST(direct > x2X0[3]);
}

/// A road on 0, 1, 2, 4, 5 straddles the inactive layer 3 through the direct
/// 2 -> 4 link: the follow proposes its row there, before row 4, and the
/// rows beyond the road outward, 6 and 7.
BOOST_AUTO_TEST_CASE(TheFollowTakesTheInteriorRowAndTheTail)
{
  const auto hits = helixHits(0.42, 0.025, 1.8);
  BOOST_REQUIRE(hits.size() == kLayers);
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, {0, 1, 2, 4, 5}, pool);
  const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  const auto proposal = followTrack(fixture, tracks, used, followParameters());
  BOOST_TEST(proposal.interior == fixture.row(3));
  BOOST_TEST(proposal.interiorAt == 3U);
  BOOST_REQUIRE_EQUAL(static_cast<int>(proposal.outwardCount), 2);
  BOOST_TEST(proposal.outward[0] == fixture.row(6));
  BOOST_TEST(proposal.outward[1] == fixture.row(7));
  BOOST_TEST(proposal.inwardCount == 0U);
}

/// Without their rows the inactive layers cost nothing: no interior row, no
/// hole on layer 7 with none allowed, and the track keeps layer 6's row.
BOOST_AUTO_TEST_CASE(WithoutTheirRowsTheInactiveLayersCostNothing)
{
  const auto hits = helixHits(0.42, 0.025, 1.8);
  BOOST_REQUIRE(hits.size() == kLayers);
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, {0, 1, 2, 4, 5}, pool);
  std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  used[fixture.row(3)] = kExclusive;
  used[fixture.row(7)] = kExclusive;
  const auto proposal = followTrack(fixture, tracks, used, followParameters());
  BOOST_TEST(proposal.interior == ClusterList<kMaxClusters>::kNoCluster);
  BOOST_REQUIRE_EQUAL(static_cast<int>(proposal.outwardCount), 1);
  BOOST_TEST(proposal.outward[0] == fixture.row(6));
}

/// The tail: where the outward walk ends before an inactive layer -- here a
/// hole on layer 6, whose row another track holds, with none allowed -- a
/// step from the walk's last state onto the inactive layers it has no row on
/// takes layer 7's row.
BOOST_AUTO_TEST_CASE(TheTailReachesTheInactiveLayerPastTheWalksEnd)
{
  const auto hits = helixHits(0.42, 0.025, 1.8);
  BOOST_REQUIRE(hits.size() == kLayers);
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, {0, 1, 2, 3, 4, 5}, pool);
  std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  used[fixture.row(6)] = kExclusive;
  const auto proposal = followTrack(fixture, tracks, used, followParameters());
  BOOST_REQUIRE_EQUAL(static_cast<int>(proposal.outwardCount), 1);
  BOOST_TEST(proposal.outward[0] == fixture.row(7));
  // Taken too, the tail finds nothing.
  used[fixture.row(7)] = kExclusive;
  BOOST_TEST(followTrack(fixture, tracks, used, followParameters()).outwardCount == 0U);
}

/// Layer 6's hole drops the row on layer 5. The tail can take layer 7
/// whose time overlaps the retained rows but is disjoint from layer 5's.
BOOST_AUTO_TEST_CASE(TheTailUsesTheTimeBoxBeforeTheUnconfirmedRow)
{
  auto hits = helixHits(0.42, 0.025, 1.8);
  BOOST_REQUIRE(hits.size() == kLayers);
  hits[5].time = {.start = 4'000U, .width = 2'000U};
  hits[7].time = {.start = 10'000U, .width = 2'000U};
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, {0, 1, 2, 3}, pool);
  std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  used[fixture.row(6)] = kExclusive;
  const auto proposal = followTrack(fixture, tracks, used, followParameters());
  BOOST_REQUIRE_EQUAL(static_cast<int>(proposal.outwardCount), 2);
  BOOST_TEST(proposal.outward[0] == fixture.row(4));
  BOOST_TEST(proposal.outward[1] == fixture.row(7));
}

/// A road already on the inactive layer has no interior row to find: the
/// follow adds what lies beyond it, as on any layer.
BOOST_AUTO_TEST_CASE(ARoadOnTheInactiveLayerHasNoInteriorRow)
{
  const auto hits = helixHits(1.1, -0.02, 2.2);
  BOOST_REQUIRE(hits.size() == kLayers);
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  const auto tracks = fittedTrack(fixture, {0, 1, 2, 3, 4, 5}, pool);
  const std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  const auto proposal = followTrack(fixture, tracks, used, followParameters());
  BOOST_TEST(proposal.interior == ClusterList<kMaxClusters>::kNoCluster);
  BOOST_REQUIRE_EQUAL(static_cast<int>(proposal.outwardCount), 2);
  BOOST_TEST(proposal.outward[1] == fixture.row(7));
}

/// The follow's passes (`Follow::run`): the track takes its interior row in
/// path order, between rows 2 and 4, and the outward rows, and owns them.
BOOST_AUTO_TEST_CASE(TheTrackHoldsItsInteriorRowInPathOrderAndOwnsIt)
{
  const auto hits = helixHits(0.42, 0.025, 1.8);
  BOOST_REQUIRE(hits.size() == kLayers);
  const Fixture fixture{hits};
  const auto pool = std::make_shared<BoundedMemoryResource>();
  auto tracks = fittedTrack(fixture, {0, 1, 2, 4, 5}, pool);
  std::vector<std::uint8_t> used(fixture.frame->size(), kFree);
  const auto parameters = followParameters();
  Follow<Det> follow;
  follow.setMemoryPool(pool);
  const auto status = follow.run(SerialForEach{}, parameters, *fixture.frame, trackingGeometry(), UniformMaterial{}, Propagator{}, tracks, 0U, true, false, used, frameBounds(*fixture.frame).widthPS, cells().view());
  BOOST_REQUIRE(status == StageStatus::kComplete);
  const auto view = tracks.view();
  BOOST_REQUIRE_EQUAL(static_cast<std::size_t>(view.count[0]), kLayers);
  for (std::size_t layer = 0; layer < kLayers; ++layer) {
    BOOST_TEST(view.clusters[0].index[layer] == fixture.row(layer));
    // The follow owns the rows it adds; the road's are the selection's.
    const bool added = layer == 3U || layer == 6U || layer == 7U;
    BOOST_TEST(used[fixture.row(layer)] == (added ? kExclusive : kFree));
  }
}

} // namespace o2::gtrk
