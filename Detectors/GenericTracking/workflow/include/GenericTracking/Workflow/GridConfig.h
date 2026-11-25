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

/// \file GridConfig.h
/// \brief One mapping from workflow configuration to the measurement grid
///        and the tracking iterations.

#ifndef GTRK_WORKFLOWS_GRIDCONFIG_H_
#define GTRK_WORKFLOWS_GRIDCONFIG_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/TimeEstPS.h"
#include "GenericTracking/Tracking/GridSizing.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/LinkMaterialFill.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"
#include "GenericTracking/Tracking/Parallel.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/SlicePlan.h"
#include "GenericTracking/Tracking/VertexStage.h"
#include "GenericTracking/Workflow/BeamLine.h"
#include "GenericTracking/Workflow/RecoParam.h"

#include "CommonConstants/LHCConstants.h"
#include "DetectorsBase/Propagator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <iterator>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace o2::gtrk::workflow
{

/// The grid's parameters: each layer's readout \a frames (one per layer), and
/// the whole number of the longest frame nearest the configured block width
/// (one at least), the blocks' edges a
/// picosecond before that frame's starts (`GridParams::blockOrigin`), in parts
/// of a picosecond that make a bunch crossing whole. With several layers of
/// the longest length the first one's starts. The other layers' frames may
/// straddle a block's edge: correct, only not minimal. A width of 0 stays one
/// block over the timeframe.
[[nodiscard]] inline GridParams gridParams(const RecoSettings& param, std::span<const FrameTiming> frames)
{
  if (frames.size() > GridParams::kMaxLayers) {
    throw std::invalid_argument(std::format("{} layers of readout frames, the grid holds {}", frames.size(), GridParams::kMaxLayers));
  }
  GridParams params{.blockWidth = std::max<std::int64_t>(param.blockWidthPS, 0)};
  std::ranges::copy(frames, std::begin(params.frames));
  const auto longest = std::ranges::max_element(frames, {}, &FrameTiming::lengthBC);
  if (longest == frames.end() || param.blockWidthPS <= 0) {
    return params;
  }
  // A BC is kPSPerBCNumerator / kPSPerBCDenominator ps (TimeEstPS.h): the
  // denominator's parts a ps, the numerator's a BC.
  constexpr std::int64_t kPartsPerPS = kPSPerBCDenominator;
  constexpr std::int64_t kPartsPerBC = kPSPerBCNumerator;
  const double framePS = static_cast<double>(longest->lengthBC) * static_cast<double>(kPartsPerBC) / static_cast<double>(kPartsPerPS);
  const auto perBlock = std::max<std::int64_t>(1, std::llround(static_cast<double>(param.blockWidthPS) / framePS));
  params.partsPerPS = kPartsPerPS;
  params.blockWidth = perBlock * longest->lengthBC * kPartsPerBC;
  // A box's start is its frame's, rounded down by less than a picosecond.
  params.blockOrigin = (longest->firstStartBC * kPartsPerBC) - kPartsPerPS;
  return params;
}

/// The field every single-value use takes -- the windows, the pair cuts, the
/// grid, the fit's propagator, the MC truth -- in kG: O2's nominal Bz, the run's
/// map at the origin as o2::base::Propagator gives it, which is what O2's own ITS
/// tracker uses. The field must have been loaded first:
/// by GRPGeomHelper in a workflow (`fieldOnlyRequest`), by the tool itself
/// outside DPL. o2::base::Propagator is fatal otherwise.
[[nodiscard]] inline double nominalBzKiloGauss() { return static_cast<double>(o2::base::Propagator::Instance()->getNominalBz()); }

/// The window physics every iteration shares; each iteration sets its own
/// minimum pT (`iterationSettings`).
[[nodiscard]] inline SearchWindowParams searchWindowPhysics(const RecoSettings& param)
{
  return {
    .bzKiloGauss = nominalBzKiloGauss(),
    .nSigma = param.windowNSigma,
  };
}

/// Refuse a configuration \a Det's tracking cannot run with: margins of
/// another length than one per layer or none, a negative margin. Called
/// where a configuration is first read: the producer's init and the bench.
template <DetectorLayout Det>
inline void checkRecoParam(const RecoSettings& param)
{
  const auto margins = param.frameMarginBC.size();
  if (margins != 0U && margins != static_cast<std::size_t>(Detector<Det>::nLayers)) {
    throw std::invalid_argument(std::format("{}.frameMarginBC has {} entries for {}'s {} layers", RecoParamOf<Det>::type::Instance().getName(), margins, Det::name, Detector<Det>::nLayers));
  }
  for (const auto margin : param.frameMarginBC) {
    if (margin < 0) {
      throw std::invalid_argument(std::format("GTRKReco<DET>.frameMarginBC holds {}: a margin widens a frame", margin));
    }
  }
  if (!std::isfinite(param.beamZReachNSigma) || !(param.beamZReachNSigma > 0.) || !(param.windowNSigma > 0.)) {
    throw std::invalid_argument(std::format("GTRKReco<DET>.beamZReachNSigma is {} and windowNSigma {}: both must be positive", param.beamZReachNSigma, param.windowNSigma));
  }
}

/// The trackers' beam spot from \a beam, reaching `beamZReachNSigma` of its z
/// width (BeamLine.h).
[[nodiscard]] inline VertexConstraint beamSpotConstraint(const BeamLine& beam, const RecoSettings& param) { return beamSpotConstraint(beam, param.beamZReachNSigma, param.windowNSigma); }

/// The tracker's memory pool: bounded at `memoryLimitGB` if that is positive,
/// unbounded otherwise.
[[nodiscard]] inline std::shared_ptr<BoundedMemoryResource> trackerMemoryPool(const RecoSettings& param)
{
  if (!(param.memoryLimitGB > 0.)) {
    return std::make_shared<BoundedMemoryResource>();
  }
  constexpr double kBytesPerGB = 1024. * 1024. * 1024.;
  const double bytes = std::floor(param.memoryLimitGB * kBytesPerGB);
  constexpr auto kLargest = std::numeric_limits<std::size_t>::max();
  return std::make_shared<BoundedMemoryResource>(bytes >= static_cast<double>(kLargest) ? kLargest : static_cast<std::size_t>(bytes));
}

/// The iterations \a param asks for: one per entry of `iterationMinPtGeV`, at least one.
/// Each other per-iteration vector must have one entry per iteration or none
/// (`RecoSettings::kClusterSlots` per iteration for `iterationClusterMinPtGeV`); anything else is a
/// configuration error.
[[nodiscard]] inline int iterationCount(const RecoSettings& param)
{
  const auto count = param.iterationMinPtGeV.size();
  if (count == 0U) {
    throw std::invalid_argument("GTRKReco<DET>.iterationMinPtGeV has no iteration");
  }
  if (count == 0U) {
    throw std::invalid_argument("GTRKReco<DET>.iterationMinPtGeV is empty: no iteration to run");
  }
  const auto check = [&](const char* name, std::size_t size, std::size_t perIteration) {
    if (size != 0U && size != count * perIteration) {
      throw std::invalid_argument(std::format("GTRKReco<DET>.{} has {} entries; {} iterations want {} or none", name, size, count, count * perIteration));
    }
  };
  check("iterationVertexSource", param.iterationVertexSource.size(), 1U);
  check("iterationMinClusters", param.iterationMinClusters.size(), 1U);
  check("iterationMaxKinks", param.iterationMaxKinks.size(), 1U);
  check("iterationKinkSigma", param.iterationKinkSigma.size(), 1U);
  check("iterationGraphKinkSigma", param.iterationGraphKinkSigma.size(), 1U);
  check("iterationSkipLinks", param.iterationSkipLinks.size(), 1U);
  check("iterationLayerMask", param.iterationLayerMask.size(), 1U);
  check("iterationRoadEndLayers", param.iterationRoadEndLayers.size(), 1U);
  check("iterationRoadStartLayers", param.iterationRoadStartLayers.size(), 1U);
  check("iterationSharedInnermost", param.iterationSharedInnermost.size(), 1U);
  check("iterationFollow", param.iterationFollow.size(), 1U);
  check("iterationFollowHoles", param.iterationFollowHoles.size(), 1U);
  check("iterationOverlapPartners", param.iterationOverlapPartners.size(), 1U);
  check("iterationInnerTrim", param.iterationInnerTrim.size(), 1U);
  check("iterationShapeTest", param.iterationShapeTest.size(), 1U);
  check("iterationMinPGeV", param.iterationMinPGeV.size(), 1U);
  check("iterationFollowEndMaxSqrtDet", param.iterationFollowEndMaxSqrtDet.size(), 1U);
  check("iterationFollowKeepDensity", param.iterationFollowKeepDensity.size(), 1U);
  check("iterationMinClustersAfterFollow", param.iterationMinClustersAfterFollow.size(), 1U);
  check("iterationBeamConstraint", param.iterationBeamConstraint.size(), 1U);
  check("iterationFollowEarlier", param.iterationFollowEarlier.size(), 1U);
  const auto floors = param.iterationClusterMinPtGeV.size();
  if (floors != 0U && floors != count * RecoSettings::kClusterSlots) {
    throw std::invalid_argument(std::format("GTRKReco<DET>.iterationClusterMinPtGeV has {} entries; {} iterations want {} per iteration, or none", floors, count, RecoSettings::kClusterSlots));
  }
  return static_cast<int>(count);
}

/// Iteration \a iteration's settings: the grid's window physics with the
/// iteration's minimum pT, and each other knob its vector's entry or default.
[[nodiscard]] inline IterationSettings iterationSettings(const RecoSettings& param, int iteration)
{
  if (iteration < 0 || iteration >= iterationCount(param)) {
    throw std::out_of_range(std::format("iteration {} of {}", iteration, param.iterationMinPtGeV.size()));
  }
  const auto index = static_cast<std::size_t>(iteration);
  // A key without entries takes IterationSettings' default.
  const IterationSettings defaults{};
  const auto entry = [&](const std::vector<int>& values, int fallback) { return values.empty() ? fallback : values[index]; };
  const auto flag = [&](const std::vector<int>& values, bool fallback) { return values.empty() ? fallback : values[index] != 0; };
  const auto entryDouble = [&](const std::vector<double>& values, double fallback) { return values.empty() ? fallback : values[index]; };
  // The layer masks: 64 bits, one per layer of the largest layout.
  const auto entryMask = [&](const std::vector<long>& values, std::uint64_t fallback) { return values.empty() ? fallback : static_cast<std::uint64_t>(values[index]); };
  auto windows = searchWindowPhysics(param);
  windows.minPtGeV = param.iterationMinPtGeV[index];
  static_assert(RecoSettings::kClusterSlots == IterationSettings::kClusterSlots);
  std::array<double, IterationSettings::kClusterSlots> floors{};
  if (!param.iterationClusterMinPtGeV.empty()) {
    for (std::size_t count = 0; count < floors.size(); ++count) {
      floors[count] = param.iterationClusterMinPtGeV[(index * floors.size()) + count];
    }
  }
  const bool follow = flag(param.iterationFollow, defaults.follow);
  const bool followEarlier = flag(param.iterationFollowEarlier, defaults.followEarlier);
  if (followEarlier && !follow) {
    throw std::invalid_argument(std::format("GTRKReco<DET>.iterationFollowEarlier of iteration {} needs its iterationFollow", iteration));
  }
  const int followHoles = entry(param.iterationFollowHoles, defaults.followHoles);
  if (followHoles < 0) {
    throw std::invalid_argument(std::format("GTRKReco<DET>.iterationFollowHoles of iteration {} is {}", iteration, followHoles));
  }
  const bool beamConstraint = flag(param.iterationBeamConstraint, defaults.beamConstraint);
  const int innerTrim = entry(param.iterationInnerTrim, defaults.innerTrim);
  if (innerTrim < 0) {
    throw std::invalid_argument(std::format("GTRKReco<DET>.iterationInnerTrim of iteration {} is {}", iteration, innerTrim));
  }
  const double kinkSigma = entryDouble(param.iterationKinkSigma, defaults.kinkSigma);
  return {
    .windows = windows,
    .vertexSource = entry(param.iterationVertexSource, static_cast<int>(defaults.vertexSource)) == 0 ? VertexSource::kBeamLine : VertexSource::kVertexList,
    .minClusters = entry(param.iterationMinClusters, defaults.minClusters),
    .maxKinks = entry(param.iterationMaxKinks, defaults.maxKinks),
    .kinkSigma = kinkSigma,
    .graphKinkSigma = entryDouble(param.iterationGraphKinkSigma, defaults.graphKinkSigma),
    .skipLinks = flag(param.iterationSkipLinks, defaults.skipLinks),
    .layerMask = entryMask(param.iterationLayerMask, defaults.layerMask),
    .roadEndLayers = entryMask(param.iterationRoadEndLayers, defaults.roadEndLayers),
    .roadStartLayers = entryMask(param.iterationRoadStartLayers, defaults.roadStartLayers),
    .minPtGeVByClusters = floors,
    .sharedInnermost = static_cast<std::uint32_t>(std::max(entry(param.iterationSharedInnermost, static_cast<int>(defaults.sharedInnermost)), 0)),
    .follow = follow,
    .followHoles = followHoles,
    .overlapPartners = flag(param.iterationOverlapPartners, defaults.overlapPartners),
    .innerTrim = innerTrim,
    .followEndMaxSqrtDet = entryDouble(param.iterationFollowEndMaxSqrtDet, defaults.followEndMaxSqrtDet),
    .followKeepDensity = entryDouble(param.iterationFollowKeepDensity, defaults.followKeepDensity),
    .minClustersAfterFollow = std::max(entry(param.iterationMinClustersAfterFollow, defaults.minClustersAfterFollow), 0),
    .beamConstraint = beamConstraint,
    .followEarlier = followEarlier,
    .shapeTest = flag(param.iterationShapeTest, defaults.shapeTest),
    .minPGeV = entryDouble(param.iterationMinPGeV, defaults.minPGeV),
  };
}

/// The own vertexer's settings, every member
/// from the run's settings; VertexerSettings' defaults are the tests'.
[[nodiscard]] inline VertexerSettings vertexerSettings(const RecoSettings& param)
{
  return {
    .maxChi2RZ = static_cast<Scalar>(param.vertexerMaxChi2RZ),
    .maxDcaCM = static_cast<Scalar>(param.vertexerMaxDcaCM),
    .windowUM = param.vertexerWindowUM,
    .claimUM = param.vertexerClaimUM,
    .minContributors = static_cast<std::uint32_t>(std::max(param.vertexerMinContributors, 0)),
    .aloneMinContributors = static_cast<std::uint32_t>(std::max(param.vertexerAloneMinContributors, 0)),
    .significance = static_cast<Scalar>(param.vertexerSignificance),
    .clipNSigma = static_cast<Scalar>(param.windowNSigma),
  };
}

/// The vertexing pass: an iteration at the beam
/// line whose links join the three innermost layers and skip none, at
/// `vertexerMinPtGeV`. It runs its tracklets and triplets only, so the
/// road, kink and follow knobs are never read; its chains are triplets,
/// 3 clusters, which keeps its links (`takenLinks`).
[[nodiscard]] inline IterationSettings vertexingPassSettings(const RecoSettings& param)
{
  auto windows = searchWindowPhysics(param);
  windows.minPtGeV = param.vertexerMinPtGeV;
  return {
    .windows = windows,
    .vertexSource = VertexSource::kBeamLine,
    .minClusters = 3,
    .maxKinks = 0,
    .skipLinks = false,
    .layerMask = 0b111U,
    .follow = false,
    .shapeTest = param.vertexerShapeTest != 0,
  };
}

/// The iterations whose queries size the grid: the tracking's, and the
/// vertexing pass with the vertexer on.
[[nodiscard]] inline std::vector<IterationSettings> gridIterations(const RecoSettings& param)
{
  std::vector<IterationSettings> iterations;
  iterations.reserve(static_cast<std::size_t>(iterationCount(param)) + 1U);
  for (int iteration = 0; iteration < iterationCount(param); ++iteration) {
    iterations.push_back(iterationSettings(param, iteration));
  }
  if (param.vertexer) {
    iterations.push_back(vertexingPassSettings(param));
  }
  return iterations;
}

/// Whether a pass tests the measurements' cluster shapes (`iterationShapeTest`,
/// `vertexerShapeTest` with the vertexer on): the measurement producer then
/// takes the detector's cluster-shape band (ClusterShapeBand.h).
[[nodiscard]] inline bool testsShapes(const RecoSettings& param)
{
  return (param.vertexer && param.vertexerShapeTest != 0) || std::ranges::any_of(param.iterationShapeTest, [](int test) { return test != 0; });
}

/// What a timeframe's measurements are relative to, which the producer ships
/// beside them (DataSpecs.h's kMeasurementGeometry): the beam line subtracted
/// from their x and y, and the grid they are binned in, built from the
/// configuration and the timeframe's conditions; and their bounds over the
/// whole frame (`frameBounds`), which the producer takes as it places them.
/// Whoever reads the measurements takes all three from here.
template <DetectorLayout Det>
struct MeasurementGeometry {
  BeamLine beamLine;
  GridGeometry<Det> grid;
  HaloBounds<Det> bounds;
};

/// A tracker's sources as the log lines name them: "beam line" or "vertex
/// list" when every iteration takes the one, both when the iterations differ.
[[nodiscard]] inline std::string sourceLabel(std::span<const IterationSettings> settings)
{
  const auto from = [&](VertexSource source) { return std::ranges::any_of(settings, [&](const auto& s) { return s.vertexSource == source; }); };
  const bool beam = from(VertexSource::kBeamLine);
  const bool list = from(VertexSource::kVertexList);
  if (beam && list) {
    return "beam line and vertex list";
  }
  return list ? "vertex list" : "beam line";
}

/// A timeframe of \a orbits orbits' span [ps], its bunch crossings' rounded down.
[[nodiscard]] inline std::int64_t timeframeDurationPS(int orbits) { return psOfBC(static_cast<std::int64_t>(orbits) * o2::constants::lhc::LHCMaxBunches); }

/// \a table filled from the material \a lut (`fillLinkMaterial`) over the
/// vertices' z the run's passes reach about \a beamSpot (`gridIterations`:
/// every iteration's and the vertexing pass's): the trackers' scattering, and
/// what sizes the grid's cells (`makeGridGeometry`). Returns that z range.
template <DetectorLayout Det, typename Lut>
Interval fillRunLinkMaterial(LinkMaterial<Det>& table, const RecoSettings& param, const VertexConstraint& beamSpot, const Lut& lut)
{
  Interval vertexZ{.low = beamSpot.z, .high = beamSpot.z};
  for (const auto& iteration : gridIterations(param)) {
    const auto reach = vertexReach(iteration.windows, beamSpot);
    vertexZ.low = std::min(vertexZ.low, reach.low);
    vertexZ.high = std::max(vertexZ.high, reach.high);
  }
  detail::TaskArena arena{std::max(param.threads, 0)};
  fillLinkMaterial<Det>(table, lut, vertexZ, arena);
  return vertexZ;
}

/// The grid of the configured iterations: each layer's cells sized by the queries the iterations make
/// into it (GridSizing.h), the vertexing pass's among them (`gridIterations`), with the run's
/// \a material (`fillRunLinkMaterial`). A process derives it once, when it has its
/// configuration and conditions, and keeps it. \a frames: each layer's
/// readout frames, which its time blocks follow and the slices' reach reads;
/// \a orbitsPerTimeframe: the run's (`Conditions::orbitsPerTimeframe`).
template <DetectorLayout Det>
[[nodiscard]] GridGeometry<Det> makeGridGeometry(const RecoSettings& param, const BeamLine& beam, std::span<const FrameTiming> frames, int orbitsPerTimeframe, const LinkMaterial<Det>& material)
{
  const auto iterations = gridIterations(param);
  if (frames.size() != static_cast<std::size_t>(Detector<Det>::nLayers)) {
    throw std::invalid_argument(std::format("{}'s grid follows its {} layers' readout frames, and {} were given", Det::name, Detector<Det>::nLayers, frames.size()));
  }
  const auto windows = gridLayerWindows<Det>(iterations, beamSpotConstraint(beam, param), material);
  GridGeometry<Det> geometry{timeframeDurationPS(orbitsPerTimeframe), gridParams(param, frames), windows};
  if (geometry.overflows()) {
    throw std::invalid_argument("GTRKReco<DET>'s grid has more cells or buckets than its 32-bit numbering holds: wider time blocks");
  }
  return geometry;
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOWS_GRIDCONFIG_H_
