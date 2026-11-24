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
///        layer mapping, the shape band's layers, native settings and grid,
///        frame times, measurement depths and IOTOF cluster conversion
///        (WorkflowTRK/TRKTOF.h).

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK TRKTOF
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/TimeEstPS.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Detectors/TRKFT3TOF.h"
#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Workflow/BeamLine.h"
#include "GenericTracking/Workflow/ClusterShapeBand.h"
#include "GenericTracking/Workflow/GeometryProvider.h"
#include "GenericTracking/Workflow/GridConfig.h"
#include "GenericTracking/Workflow/RecoParam.h"
#include "GenericTracking/WorkflowTRK/TRKTOF.h"

#include "BandMaterial.h"

#include "CommonConstants/LHCConstants.h"
#include "DataFormatsIOTOF/Cluster.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "IOTOFReconstruction/TopologyClassifier.h"

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <span>
#include <vector>

namespace
{
using namespace o2::gtrk;
using Base = detectors::TRKFT3;
using TOFLayout = detectors::TRKFT3TOF;
constexpr int kBaseLayers = Detector<Base>::nLayers;

} // namespace

/// TRK's shape band is TRKFT3's layers: TRKFT3TOF reads each TRK layer's,
/// and its TOF layers none (every slope).
BOOST_AUTO_TEST_CASE(TRKFT3TOFReadsTRKsShapeBandInTRKFT3sLayers)
{
  using Layers = workflow::ShapeBandLayers<TOFLayout>;
  for (int layer = 0; layer < kBaseLayers; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    BOOST_TEST(Layers::layer(detectors::trkft3tof::fromTRKFT3(id)) == id);
  }
  BOOST_TEST(Layers::layer(TOFLayout::kInnerTOFLayer) == kInvalidLayerID);
  BOOST_TEST(Layers::layer(TOFLayout::kOuterTOFLayer) == kInvalidLayerID);
  // A band calibrating TRKFT3's ML 7 alone: two rows at most, one column at slopes [0, 1).
  constexpr std::size_t kML7 = 7;
  workflow::ClusterShapeBand band;
  band.maxRows.assign(Layers::barrelLayers, 0);
  band.maxRows[kML7] = 2;
  band.slopeLow.assign(Layers::barrelLayers * workflow::ClusterShapeBand::kSpans, 0.F);
  band.slopeHigh.assign(Layers::barrelLayers * workflow::ClusterShapeBand::kSpans, 0.F);
  band.slopeHigh[(kML7 * workflow::ClusterShapeBand::kSpans) + 1] = 1.F;
  BOOST_REQUIRE(band.fits(Layers::barrelLayers));
  const auto onLayer = [&](LayerID layer) { return band.interval(Layers::layer(layer), 1, 1); };
  BOOST_TEST((onLayer(detectors::trkft3tof::fromTRKFT3(static_cast<LayerID>(kML7))) == workflow::ClusterShapeBand::Interval{0.F, 1.F}));
  BOOST_TEST((onLayer(TOFLayout::kInnerTOFLayer) == workflow::ClusterShapeBand::kEvery));
}

/// Every TRKFT3 link is TRKFT3TOF's, renumbered, with its source range; ML 6
/// -> ML 7 stays direct, beside ML 6 -> iTOF -> ML 7.
BOOST_AUTO_TEST_CASE(EveryTRKFT3LinkIsTRKFT3TOFs)
{
  for (const auto& link : Base::layerLinks) {
    const auto* const found = std::ranges::find_if(TOFLayout::layerLinks, [&](const LayerLink& other) { return other.from == detectors::trkft3tof::fromTRKFT3(link.from) && other.to == detectors::trkft3tof::fromTRKFT3(link.to); });
    BOOST_TEST_CONTEXT("link " << static_cast<int>(link.from) << " -> " << static_cast<int>(link.to))
    {
      BOOST_REQUIRE(found != TOFLayout::layerLinks.end());
      BOOST_TEST(found->sourceMin == link.sourceMin);
      BOOST_TEST(found->sourceMax == link.sourceMax);
    }
  }
  const auto has = [](LayerID from, LayerID to) { return std::ranges::any_of(TOFLayout::layerLinks, [&](const LayerLink& link) { return link.from == from && link.to == to; }); };
  BOOST_TEST(has(6, detectors::trkft3tof::fromTRKFT3(7)));
  BOOST_TEST(has(6, TOFLayout::kInnerTOFLayer));
  BOOST_TEST(has(TOFLayout::kInnerTOFLayer, detectors::trkft3tof::fromTRKFT3(7)));
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
  BOOST_REQUIRE(tofFrames.size() == static_cast<std::size_t>(Detector<TOFLayout>::nLayers));
  BOOST_REQUIRE(tofDepths.size() == static_cast<std::size_t>(Detector<TOFLayout>::nLayers));
  for (int layer = 0; layer < kBaseLayers; ++layer) {
    const auto to = detectors::trkft3tof::fromTRKFT3(static_cast<LayerID>(layer));
    BOOST_TEST(tofFrames[to].lengthBC == frames[static_cast<std::size_t>(layer)].lengthBC);
    BOOST_TEST(tofFrames[to].firstStartBC == frames[static_cast<std::size_t>(layer)].firstStartBC);
    BOOST_TEST(tofDepths[to] == depths[static_cast<std::size_t>(layer)]);
  }
  for (const auto layer : {TOFLayout::kInnerTOFLayer, TOFLayout::kOuterTOFLayer}) {
    BOOST_TEST(tofFrames[layer].lengthBC == 1);
    BOOST_TEST(tofFrames[layer].firstStartBC == 0);
    BOOST_TEST(tofDepths[layer] == 0.);
  }
  BOOST_CHECK_THROW((void)alice::trktofFrameTiming(std::span{frames}.first(3)), std::invalid_argument);
}

/// The native grid and tracker accept pT floors through the longest TOF
/// chain.
BOOST_AUTO_TEST_CASE(TheNativeGridAndFloorsCoverLongTOFChains)
{
  static_assert(maxChainLayers<Base> == 15U && maxChainLayers<TOFLayout> == 17U);
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
  const auto tofGeometry = workflow::makeGeometry<TOFLayout>();
  const std::vector<FrameTiming> trkFrames(static_cast<std::size_t>(kBaseLayers), FrameTiming{.lengthBC = 81, .firstStartBC = 0});
  const auto tofFrames = alice::trktofFrameTiming(trkFrames);
  const auto beam = workflow::nominalBeamLine();
  const auto grid = std::make_unique<const GridGeometry<TOFLayout>>(workflow::makeGridGeometry<TOFLayout>(settings, beam, tofFrames, 32, noMaterial<TOFLayout>()));
  BOOST_TEST(!grid->overflows());
  BOOST_TEST(grid->axes(TOFLayout::kInnerTOFLayer).nPhi > 0);
  const auto& scattering = noMaterial<TOFLayout>();
  BOOST_TEST(iterationParameters<TOFLayout>(first, VertexConstraint{}, scattering).minPtByClusters[17] == 0.25F);
  BOOST_TEST(iterationParameters<TOFLayout>(second, VertexConstraint{}, scattering).minPtByClusters[16] == 0.5F);
}

/// An IOTOF cluster as ours: its fired pixels' mean centre in its chip's
/// frame, pitch^2/12 errors, the sensor after TRKFT3's, the box of the bunch
/// crossing ending at its time; one on no chip, outside its crossing,
/// without a usable pattern (O2's `kHuge` among them, a full 4 x 4 pattern
/// not) or outside the timeframe left out and counted; a box starting before
/// the timeframe started at its start and counted.
BOOST_AUTO_TEST_CASE(ATOFClusterIsItsPixelsAndTheCrossingEndingAtItsTime)
{
  alice::TOFPixelFrame frame;
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
  const std::uint64_t endBC = firstBC + (std::uint64_t{kOrbits} * o2::constants::lhc::LHCMaxBunches);
  using o2::iotof::Cluster;
  std::vector<Cluster> clusters{
    // one pixel, at row 2, column 3, 100 BC into the timeframe, 7.5 ns into its crossing
    Cluster(2, 3, 1, 1, 0b1, 0, 1, 0.F, firstBC + 100, 7.5F),
    // two pixels side by side across the rows (row-major 2 x 1: bits 0 and 1)
    Cluster(4, 1, 2, 1, 0b11, 0, 3, 0.F, firstBC + 2000, 0.F),
    Cluster(0, 0, 1, 1, 0b1, 0, 4, 0.F, firstBC + 10, 1.F),                                                  // no such chip
    Cluster(0, 0, 1, 1, 0b1, 0, 0, 0.F, firstBC + 10, 30.F),                                                 // outside its crossing
    Cluster(0, 0, 1, 1, 0, 0, 0, 0.F, firstBC + 10, 1.F),                                                    // no pixel fired
    Cluster(7, 5, 2, 1, 0b11, 0, 0, 0.F, firstBC + 10, 1.F),                                                 // a pixel past the last row
    Cluster(0, 0, 1, 1, 0b1, 0, 2, 0.F, firstBC, 1.F),                                                       // a box from before the timeframe
    Cluster(0, 0, 4, 4, 0xFFFF, 0, 2, 0.F, firstBC + 50, 2.F),                                               // a full 4 x 4 pattern
    Cluster(0, 0, 4, 4, Cluster::InvalidPatternID, o2::iotof::Topologies::kHuge, 2, 0.F, firstBC + 60, 2.F), // too large for a pattern
    Cluster(0, 0, 1, 1, 0b1, 0, 2, 0.F, firstBC - 3, 1.F),                                                   // before the timeframe
    Cluster(0, 0, 1, 1, 0b1, 0, 2, 0.F, endBC + 3, 1.F),                                                     // after it
  };
  std::vector<std::uint32_t> kept;
  alice::TOFSkipped skipped;
  const auto ours = alice::tofClusters(clusters, frame, kFirstOrbit, kOrbits, kept, skipped);
  BOOST_REQUIRE(ours.size() == 4U);
  // In box-start order: the box started at the timeframe's start first.
  BOOST_TEST(kept == (std::vector<std::uint32_t>{6, 7, 0, 1}), boost::test_tools::per_element());
  BOOST_TEST(skipped.chip == 1U);
  BOOST_TEST(skipped.time == 1U);
  BOOST_TEST(skipped.pattern == 3U);
  BOOST_TEST(skipped.outside == 2U);
  BOOST_TEST(skipped.clamped == 1U);
  BOOST_TEST(ours[1].xLocal == 0.25F * (frame.rowX[0] + frame.rowX[1] + frame.rowX[2] + frame.rowX[3]), boost::test_tools::tolerance(1.e-6F));
  BOOST_TEST(ours[1].yLocal == 0.25F * (frame.columnZ[0] + frame.columnZ[1] + frame.columnZ[2] + frame.columnZ[3]), boost::test_tools::tolerance(1.e-6F));

  BOOST_TEST(ours[2].chipID == static_cast<SensorID>(TOFLayout::kFirstTOFSensor + 1));
  BOOST_TEST(ours[2].xLocal == frame.rowX[2]);
  BOOST_TEST(ours[2].yLocal == frame.columnZ[3]);
  BOOST_TEST(ours[2].sigmaXLocal == frame.sigmaX);
  BOOST_TEST(ours[2].sigmaYLocal == frame.sigmaZ);
  const std::int64_t t0 = psOfBC(100) + 7500;
  BOOST_TEST(ours[2].timePS == static_cast<std::uint32_t>(t0 - psOfBC(1)));
  BOOST_TEST(ours[2].widthPS == static_cast<std::uint32_t>(psOfBC(1)));
  BOOST_TEST(ours[2].rowSpan == 1U);
  BOOST_TEST(ours[2].colSpan == 1U);

  BOOST_TEST(ours[3].chipID == static_cast<SensorID>(TOFLayout::kFirstTOFSensor + 3));
  BOOST_TEST(ours[3].xLocal == 0.5F * (frame.rowX[4] + frame.rowX[5]), boost::test_tools::tolerance(1.e-6F));
  BOOST_TEST(ours[3].yLocal == frame.columnZ[1]);
  BOOST_TEST(ours[3].timePS == static_cast<std::uint32_t>(psOfBC(2000) - psOfBC(1)));
  BOOST_TEST(ours[3].rowSpan == 2U);

  // [t - 1 BC, t] at t = 1 ns reaches before the timeframe: moved to its start.
  BOOST_TEST(ours[0].timePS == 0U);
  BOOST_TEST(ours[0].widthPS == static_cast<std::uint32_t>(psOfBC(1)));
}
