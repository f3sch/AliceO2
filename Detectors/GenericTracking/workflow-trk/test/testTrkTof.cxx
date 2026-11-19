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

/// \file testTrkTof.cxx
/// \brief TRKFT3TOF, ALICE 3's tracker with its TOF layers, against TRKFT3:
///        layer mapping, native settings and grid, frame times, measurement
///        depths and IOTOF cluster conversion (WorkflowTRK/TRKTOF.h).

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK TRKTOF
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/TimeEstPS.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Detectors/TRKFT3TOF.h"
#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Workflow/BeamLine.h"
#include "GenericTracking/Workflow/GridConfig.h"
#include "GenericTracking/Workflow/RecoParam.h"
#include "GenericTracking/WorkflowTRK/TRKTOF.h"

#include "CommonConstants/LHCConstants.h"
#include "DataFormatsIOTOF/Cluster.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

#include "UniformMaterial.h"

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <span>
#include <vector>

namespace
{
using namespace o2::gtrk;
using Base = detectors::TRKFT3;
using Tof = detectors::TRKFT3TOF;
constexpr int kBaseLayers = Detector<Base>::nLayers;

[[nodiscard]] const PlacedLayerSpec& layerSpec(auto layout, int layer)
{
  using L = decltype(layout);
  return layer < Detector<L>::nBarrelLayers ? L::barrelLayers[static_cast<std::size_t>(layer)] : L::discLayers[static_cast<std::size_t>(layer - Detector<L>::nBarrelLayers)];
}

} // namespace

/// TRKFT3's layers in their order with the two TOF layers among them: the
/// iTOF between ML 6 and 7, the oTOF after the OT, the discs after it; the
/// same sensors, overlaps and links, and the TOF layers inactive.
BOOST_AUTO_TEST_CASE(TRKFT3TOFIsTRKFT3WithTheTOFLayers)
{
  BOOST_TEST(Detector<Tof>::nLayers == kBaseLayers + 2);
  BOOST_TEST(Detector<Tof>::nBarrelLayers == Detector<Base>::nBarrelLayers + 2);
  BOOST_TEST(Detector<Tof>::inactiveLayers == (std::uint64_t{1} << Tof::kInnerTOFLayer | std::uint64_t{1} << Tof::kOuterTOFLayer));
  BOOST_TEST(Detector<Base>::inactiveLayers == 0U);
  // The TOF's clusters carry their own boxes; TRKFT3's are frame records.
  BOOST_TEST(Detector<Tof>::hitTimedLayers == Detector<Tof>::inactiveLayers);
  BOOST_TEST(Detector<Base>::hitTimedLayers == 0U);
  BOOST_TEST(Tof::kFirstTOFSensor == Base::nSensors);
  std::vector<bool> taken(static_cast<std::size_t>(Detector<Tof>::nLayers), false);
  for (int layer = 0; layer < kBaseLayers; ++layer) {
    const auto to = detectors::trkft3tof::fromTRKFT3(static_cast<LayerID>(layer));
    BOOST_TEST_CONTEXT("TRKFT3 layer " << layer << " -> " << static_cast<int>(to))
    {
      BOOST_TEST(!Detector<Tof>::isInactiveLayer(to));
      BOOST_TEST(!taken[to]);
      taken[to] = true;
      BOOST_TEST(Detector<Tof>::isBarrelLayer(to) == Detector<Base>::isBarrelLayer(static_cast<LayerID>(layer)));
      BOOST_TEST(layerSpec(Tof{}, to).sensorCount == layerSpec(Base{}, layer).sensorCount);
      BOOST_TEST(layerSpec(Tof{}, to).sensor.halfX == layerSpec(Base{}, layer).sensor.halfX);
      BOOST_TEST(((Tof::overlapLayers >> to) & 1U) == ((Base::overlapLayers >> layer) & 1U));
    }
  }
  BOOST_TEST(!taken[Tof::kInnerTOFLayer]);
  BOOST_TEST(!taken[Tof::kOuterTOFLayer]);
  BOOST_TEST(((Tof::overlapLayers >> Tof::kInnerTOFLayer) & 1U) == 0U);
  // Every TRKFT3 link, renumbered, with its source range.
  for (const auto& link : Base::layerLinks) {
    const auto found = std::ranges::find_if(Tof::layerLinks, [&](const LayerLink& other) { return other.from == detectors::trkft3tof::fromTRKFT3(link.from) && other.to == detectors::trkft3tof::fromTRKFT3(link.to); });
    BOOST_TEST_CONTEXT("link " << static_cast<int>(link.from) << " -> " << static_cast<int>(link.to))
    {
      BOOST_REQUIRE(found != Tof::layerLinks.end());
      BOOST_TEST(found->sourceMin == link.sourceMin);
      BOOST_TEST(found->sourceMax == link.sourceMax);
    }
  }
  // ML 6 -> ML 7 stays direct, beside ML 6 -> iTOF -> ML 7.
  const auto has = [](LayerID from, LayerID to) { return std::ranges::any_of(Tof::layerLinks, [&](const LayerLink& link) { return link.from == from && link.to == to; }); };
  BOOST_TEST(has(6, detectors::trkft3tof::fromTRKFT3(7)));
  BOOST_TEST(has(6, Tof::kInnerTOFLayer));
  BOOST_TEST(has(Tof::kInnerTOFLayer, detectors::trkft3tof::fromTRKFT3(7)));
}

/// A TOF configuration can be set independently of the TRK configuration.
BOOST_AUTO_TEST_CASE(TheTofScheduleUsesItsOwnKeys)
{
  const auto trk = recoParam<Base>().iterationMinPtGeV;
  o2::conf::ConfigurableParam::updateFromString("GTRKRecoTRKTOF.iterationMinPtGeV=[0.5];GTRKRecoTRKTOF.iterationLayerMask=[248]");
  const auto& tof = recoParam<Tof>();
  BOOST_TEST(tof.iterationMinPtGeV == (std::vector<double>{0.5}), boost::test_tools::per_element());
  BOOST_TEST(tof.iterationLayerMask == (std::vector<long>{248}), boost::test_tools::per_element());
  BOOST_TEST(recoParam<Base>().iterationMinPtGeV == trk, boost::test_tools::per_element());
}

/// The default mask admits iTOF links and keeps oTOF for the follower.
BOOST_AUTO_TEST_CASE(TheDefaultMaskKeepsOuterTofOutOfSeeding)
{
  RecoSettings tof;
  tof.iterationMinPtGeV = {0.1, 0.2};
  o2::base::Propagator::Instance(true);
  for (int iteration = 0; iteration < workflow::iterationCount(tof); ++iteration) {
    const auto links = takenLinks<Tof>(workflow::iterationSettings(tof, iteration));
    bool innerSeeded = false;
    bool outerSeeded = false;
    for (std::size_t i = 0; i < links.size(); ++i) {
      const auto& link = Tof::layerLinks[i];
      innerSeeded |= links[i] && (link.from == Tof::kInnerTOFLayer || link.to == Tof::kInnerTOFLayer);
      outerSeeded |= links[i] && (link.from == Tof::kOuterTOFLayer || link.to == Tof::kOuterTOFLayer);
    }
    BOOST_TEST(innerSeeded);
    BOOST_TEST(!outerSeeded);
  }
}

/// TRKFT3TOF's frames and depths are TRKFT3's on its layers; the TOF's a
/// bunch crossing from BC 0 at depth 0.
BOOST_AUTO_TEST_CASE(TheFramesAndDepthsAreTRKFT3sOnItsLayers)
{
  std::vector<FrameTiming> frames;
  std::vector<double> depths;
  for (int layer = 0; layer < kBaseLayers; ++layer) {
    frames.push_back({.lengthBC = 18 + layer, .firstStartBC = layer % 3});
    depths.push_back(-0.001 * (layer + 1));
  }
  const auto tofFrames = alice::trktofFrameTiming(frames);
  const auto tofDepths = alice::trktofMeasurementDepthsCM(depths);
  BOOST_REQUIRE(tofFrames.size() == static_cast<std::size_t>(Detector<Tof>::nLayers));
  BOOST_REQUIRE(tofDepths.size() == static_cast<std::size_t>(Detector<Tof>::nLayers));
  for (int layer = 0; layer < kBaseLayers; ++layer) {
    const auto to = detectors::trkft3tof::fromTRKFT3(static_cast<LayerID>(layer));
    BOOST_TEST(tofFrames[to].lengthBC == frames[static_cast<std::size_t>(layer)].lengthBC);
    BOOST_TEST(tofFrames[to].firstStartBC == frames[static_cast<std::size_t>(layer)].firstStartBC);
    BOOST_TEST(tofDepths[to] == depths[static_cast<std::size_t>(layer)]);
  }
  for (const auto layer : {Tof::kInnerTOFLayer, Tof::kOuterTOFLayer}) {
    BOOST_TEST(tofFrames[layer].lengthBC == 1);
    BOOST_TEST(tofFrames[layer].firstStartBC == 0);
    BOOST_TEST(tofDepths[layer] == 0.);
  }
  BOOST_CHECK_THROW((void)alice::trktofFrameTiming(std::span{frames}.first(3)), std::invalid_argument);
}

/// The native grid and tracker accept pT floors through the longest TOF
/// chain. The legacy width remains usable with the shorter TRK layout.
BOOST_AUTO_TEST_CASE(TheNativeGridAndFloorsCoverLongTofChains)
{
  static_assert(maxChainLayers<Base> == 15U && maxChainLayers<Tof> == 17U);
  RecoSettings settings;
  settings.iterationMinPtGeV = {0.5, 0.3};
  settings.iterationClusterMinPtGeV.assign(2U * RecoSettings::kClusterSlots, 0.);
  settings.iterationClusterMinPtGeV[17] = 0.25;
  settings.iterationClusterMinPtGeV[RecoSettings::kClusterSlots + 16] = 0.5;
  o2::base::Propagator::Instance(true)->setNominalBz(20.F);
  const auto first = workflow::iterationSettings(settings, 0);
  const auto second = workflow::iterationSettings(settings, 1);
  BOOST_TEST(first.minPtGeVByClusters[17] == 0.25);
  BOOST_TEST(second.minPtGeVByClusters[16] == 0.5);
  BOOST_TEST(second.minPtGeVByClusters[17] == 0.);

  const char* file = std::getenv("GTRK_ALICE3_GEOMETRY");
  if (file == nullptr) {
    BOOST_TEST_MESSAGE("GTRK_ALICE3_GEOMETRY unset: the geometry-dependent checks are skipped");
    return;
  }
  o2::base::GeometryManager::loadGeometry(file);
  const auto tofGeometry = std::make_unique<const Geometry<Tof>>(alice::placeChecked<Tof>());
  const std::vector<FrameTiming> trkFrames(static_cast<std::size_t>(kBaseLayers), FrameTiming{.lengthBC = 81, .firstStartBC = 0});
  const auto tofFrames = alice::trktofFrameTiming(trkFrames);
  const auto beam = workflow::nominalBeamLine();
  const auto grid = std::make_unique<const GridGeometry<Tof>>(workflow::makeGridGeometry<Tof>(settings, beam, tofFrames, 32));
  BOOST_TEST(!grid->overflows());
  BOOST_TEST(grid->axes(Tof::kInnerTOFLayer).nPhi > 0);
  const auto scattering = scatteringTable<Tof>(uniformLayerMaterial<Tof>(0.));
  BOOST_TEST(iterationParameters<Tof>(first, VertexConstraint{}, scattering).minPtByClusters[17] == 0.25F);
  BOOST_TEST(iterationParameters<Tof>(second, VertexConstraint{}, scattering).minPtByClusters[16] == 0.5F);

  settings.iterationClusterMinPtGeV.assign(2U * RecoSettings::kLegacyClusterSlots, 0.);
  settings.iterationClusterMinPtGeV[RecoSettings::kLegacyClusterSlots + 7] = 0.5;
  BOOST_TEST(workflow::iterationSettings(settings, 1).minPtGeVByClusters[7] == 0.5);
  BOOST_CHECK_THROW(workflow::makeGridGeometry<Tof>(settings, beam, tofFrames, 32), std::invalid_argument);
  const auto trkGeometry = std::make_unique<const Geometry<Base>>(alice::placeChecked<Base>());
  const auto trkGrid = std::make_unique<const GridGeometry<Base>>(workflow::makeGridGeometry<Base>(settings, beam, trkFrames, 32));
  BOOST_TEST(!trkGrid->overflows());
}

/// An IOTOF cluster as ours: its fired pixels' mean centre in its chip's
/// frame, pitch^2/12 errors, the sensor after TRKFT3's, the box of the bunch
/// crossing ending at its time; one on no chip, outside its crossing or
/// without a usable pattern left out and counted; a box past the
/// timeframe's start or end moved in and counted.
BOOST_AUTO_TEST_CASE(ATOFClusterIsItsPixelsAndTheCrossingEndingAtItsTime)
{
  alice::TofPixelFrame frame;
  for (int row = 0; row < 8; ++row) {
    frame.rowX.push_back(0.01F * static_cast<float>(row) - 0.035F);
  }
  for (int column = 0; column < 6; ++column) {
    frame.columnZ.push_back(0.025F * static_cast<float>(column) - 0.0625F);
  }
  frame.sigmaX = 0.01F / std::sqrt(12.F);
  frame.sigmaZ = 0.025F / std::sqrt(12.F);
  frame.chips = 4;
  constexpr std::uint32_t kFirstOrbit = 5;
  constexpr int kOrbits = 2;
  const std::uint64_t firstBC = std::uint64_t{kFirstOrbit} * o2::constants::lhc::LHCMaxBunches;
  using o2::iotof::Cluster;
  std::vector<Cluster> clusters{
    // one pixel, at row 2, column 3, 100 BC into the timeframe, 7.5 ns into its crossing
    Cluster(2, 3, 1, 1, 0b1, 0, 1, 0.F, firstBC + 100, 7.5F),
    // two pixels side by side across the rows (row-major 2 x 1: bits 0 and 1)
    Cluster(4, 1, 2, 1, 0b11, 0, 3, 0.F, firstBC + 2000, 0.F),
    Cluster(0, 0, 1, 1, 0b1, 0, 4, 0.F, firstBC + 10, 1.F),  // no such chip
    Cluster(0, 0, 1, 1, 0b1, 0, 0, 0.F, firstBC + 10, 30.F), // outside its crossing
    Cluster(0, 0, 1, 1, 0, 0, 0, 0.F, firstBC + 10, 1.F),    // no pixel fired
    Cluster(7, 5, 2, 1, 0b11, 0, 0, 0.F, firstBC + 10, 1.F), // a pixel past the last row
    Cluster(0, 0, 1, 1, 0b1, 0, 2, 0.F, firstBC, 1.F),       // a box from before the timeframe
  };
  std::vector<std::uint32_t> kept;
  alice::TofSkipped skipped;
  const auto ours = alice::tofClusters(clusters, frame, kFirstOrbit, kOrbits, kept, skipped);
  BOOST_REQUIRE(ours.size() == 3U);
  // In box-start order: the box moved to the timeframe's start first.
  BOOST_TEST(kept == (std::vector<std::uint32_t>{6, 0, 1}), boost::test_tools::per_element());
  BOOST_TEST(skipped.chip == 1U);
  BOOST_TEST(skipped.time == 1U);
  BOOST_TEST(skipped.pattern == 2U);
  BOOST_TEST(skipped.clamped == 1U);

  BOOST_TEST(ours[1].chipID == static_cast<SensorID>(Tof::kFirstTOFSensor + 1));
  BOOST_TEST(ours[1].xLocal == frame.rowX[2]);
  BOOST_TEST(ours[1].yLocal == frame.columnZ[3]);
  BOOST_TEST(ours[1].sigmaXLocal == frame.sigmaX);
  BOOST_TEST(ours[1].sigmaYLocal == frame.sigmaZ);
  const std::int64_t t0 = psOfBC(100) + 7500;
  BOOST_TEST(ours[1].timePS == static_cast<std::uint32_t>(t0 - psOfBC(1)));
  BOOST_TEST(ours[1].widthPS == static_cast<std::uint32_t>(psOfBC(1)));
  BOOST_TEST(ours[1].rowSpan == 1U);
  BOOST_TEST(ours[1].colSpan == 1U);

  BOOST_TEST(ours[2].chipID == static_cast<SensorID>(Tof::kFirstTOFSensor + 3));
  BOOST_TEST(ours[2].xLocal == 0.5F * (frame.rowX[4] + frame.rowX[5]), boost::test_tools::tolerance(1.e-6F));
  BOOST_TEST(ours[2].yLocal == frame.columnZ[1]);
  BOOST_TEST(ours[2].timePS == static_cast<std::uint32_t>(psOfBC(2000) - psOfBC(1)));
  BOOST_TEST(ours[2].rowSpan == 2U);

  // [t - 1 BC, t] at t = 1 ns reaches before the timeframe: moved to its start.
  BOOST_TEST(ours[0].timePS == 0U);
  BOOST_TEST(ours[0].widthPS == static_cast<std::uint32_t>(psOfBC(1)));
}
