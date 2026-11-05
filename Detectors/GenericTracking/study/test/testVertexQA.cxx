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

/// \file testVertexQA.cxx
/// \brief The vertexer's QA against the truth, on the vertexing pass's two
///        collisions.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK VertexQA
#include <boost/test/unit_test.hpp>
#include <numbers>

#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/VertexingPass.h"
#include "GenericTracking/Study/VertexQA.h"

#include "ToyBarrel8Fixture.h"
#include "UniformMaterial.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace o2::gtrk::study
{

namespace
{

using namespace toybarrel8;
using Pass = VertexingPass<Det, UniformMaterial>;

/// Two collisions of eight tracks each, 500 ns apart, every track 0.5 to 2.6
/// GeV from the origin.
constexpr std::size_t kCollisions = 2U;
constexpr std::size_t kTracksPerCollision = 8U;
constexpr std::size_t kTracks = kCollisions * kTracksPerCollision;
constexpr std::array<double, kCollisions> kStartNS{37., 537.};
constexpr double kDurationNS = 600.;

std::vector<Hit> hits()
{
  std::vector<Hit> result;
  for (std::size_t track = 0; track < kTracks; ++track) {
    const auto k = static_cast<double>(track);
    const double pt = 0.5 + (0.3 * static_cast<double>(track % kTracksPerCollision));
    const auto one =
      sensorHelixHits(std::fmod(0.3 + (2.1 * k), 2 * std::numbers::pi_v<double>), 0.03 * std::cos(1.3 * k), pt, kStartNS[track / kTracksPerCollision]);
    BOOST_REQUIRE(one.size() == kLayers);
    result.insert(result.end(), one.begin(), one.end());
  }
  return result;
}

/// The pass: the beam line, the three innermost layers, no holes, 0.3 GeV.
IterationSettings passSettings()
{
  auto settings = iterationSettings();
  settings.windows.minPtGeV = 0.3;
  settings.minClusters = 3;
  settings.maxKinks = 0;
  settings.skipLinks = false;
  settings.layerMask = 0b111U;
  settings.follow = false;
  return settings;
}

ScatteringTable<Det> noScattering()
{
  const std::array<double, kLayers> material{};
  return scatteringTable<Det>(material);
}

/// The pass's slices: three blocks, as the tracking's in these tests.
constexpr std::int32_t kSliceBlocks = 3;

} // namespace

/// The QA against the truth: each track a
/// primary of its collision, both collisions found with every findable
/// primary, nothing fake, split or merged.
BOOST_AUTO_TEST_CASE(TheQAFindsBothCollisions)
{
  const Fixture fixture{hits(), kDurationNS};
  Pass pass{passSettings(), VertexerSettings{}, Fixture::kBeamConstraint, noScattering(), std::make_shared<BoundedMemoryResource>(), 2U, 1};
  const auto result = pass.find(*fixture.frame, kSliceBlocks, 0, true);
  BOOST_REQUIRE(result.status == StageStatus::kComplete);
  BOOST_REQUIRE(!result.claims.empty());
  const auto particleOf = [](std::uint64_t track) -> std::uint64_t { return ((track / kTracksPerCollision) << 31U) | (track % kTracksPerCollision); };
  std::vector<study::TruthRow> rows(fixture.frame->size());
  study::ParticleTable particles;
  for (std::uint32_t row = 0; row < fixture.frame->size(); ++row) {
    const auto hit = fixture.frame->measurement(row).clusterIndex;
    const auto track = hit / kLayers;
    rows[row].particles[0] = particleOf(track);
    rows[row].count = 1U;
    rows[row].layer = static_cast<LayerID>(hit % kLayers);
    particles[particleOf(track)] = {.pt = 1.F, .primary = true};
  }
  const auto score = study::scoreVertices<Det>(rows, particles, result.claims, result.found);
  BOOST_TEST(score.vertices == kCollisions);
  BOOST_TEST(score.collisions == kCollisions);
  BOOST_TEST(score.found == kCollisions);
  BOOST_TEST(score.findable == kTracks);
  BOOST_TEST(score.covered == kTracks);
  BOOST_TEST(score.fakes == 0U);
  BOOST_TEST(score.splits == 0U);
  BOOST_TEST(score.merged == 0U);
  BOOST_REQUIRE(score.residualsUM.size() == kCollisions);
  for (const double residual : score.residualsUM) {
    BOOST_TEST(residual < 5.);
  }
}

} // namespace o2::gtrk::study
