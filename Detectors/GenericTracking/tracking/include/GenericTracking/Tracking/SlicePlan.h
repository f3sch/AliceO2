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

/// \file SlicePlan.h
/// \brief A timeframe's slices, planned once on the host from the frame and
///        the iterations, and run by the host and the device trackers alike.
///
/// A slice's gate window and the halo width depend on the measurements'
/// time variances, not only on the grid offsets; the tracklet stage keeps
/// them as values. The plan computes them once, with the source maps, the
/// vertex windows and the extension's per-layer bound, and records what it
/// was planned from, so that a tracker can refuse a plan of another frame or
/// configuration.

#ifndef GTRK_TRACKING_SLICEPLAN_H_
#define GTRK_TRACKING_SLICEPLAN_H_

#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/TimeEstPS.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/LayerRowMap.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TimeFrameScratch.h"
#include "GenericTracking/Tracking/TrackletStage.h"
#include "GenericTracking/Tracking/TripletCut.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace o2::gtrk
{

/// Every layer's largest values over a range of rows: the time box's width
/// the tracklet stage's and the extension's windows and the slices' reach are
/// bounded with, and the z, r and r-phi variances the triplet stage's
/// dip-angle window is, and the larger of its two local variances, which the
/// overlap partners' search spreads by (OverlapPartners.h); and whether
/// every row's box is exactly one of the layer's readout frames
/// (`FrameTiming::isFrame`), which the slices' reach reads (`chainReach`).
///
/// The bounds over a whole frame (`frameBounds`) are the measurement
/// producer's: it adds every row as it places it and ships them with the grid,
/// so that no consumer scans the frame again.
template <DetectorLayout Det>
struct HaloBounds {
  static constexpr auto kLayerCount = static_cast<std::uint32_t>(Detector<Det>::nLayers);
  std::array<std::uint32_t, kLayerCount> widthPS{};
  std::array<TripletTargetBounds, kLayerCount> target{};
  std::array<Scalar, kLayerCount> varLoc{};
  std::array<bool, kLayerCount> rowsAreFrames{};

  /// The bounds over no rows: zero, and every layer's rows (none) frames.
  [[nodiscard]] static HaloBounds none()
  {
    HaloBounds bounds;
    bounds.rowsAreFrames.fill(true);
    return bounds;
  }

  /// One more row of \a layer, whose readout frames are \a frames: its time
  /// box \a box, its z, r and r-phi variances and its two local ones.
  void add(std::uint32_t layer, const FrameTiming& frames, const TimeEstPS& box, Scalar varZ, Scalar varR, Scalar varRPhi, Scalar varLoc0, Scalar varLoc1)
  {
    auto& bounds = target[layer];
    widthPS[layer] = std::max(widthPS[layer], box.width);
    rowsAreFrames[layer] = rowsAreFrames[layer] && frames.isFrame(box);
    bounds.varZ = std::max(bounds.varZ, varZ);
    bounds.varR = std::max(bounds.varR, varR);
    bounds.varRPhi = std::max(bounds.varRPhi, varRPhi);
    varLoc[layer] = std::max({varLoc[layer], varLoc0, varLoc1});
  }
};

template <DetectorLayout Det>
[[nodiscard]] HaloBounds<Det> boundsOver(const TimeFrame<Det>& frame, std::span<const MeasurementRange, HaloBounds<Det>::kLayerCount> rows)
{
  auto result = HaloBounds<Det>::none();
  const auto columns = frame.measurements();
  const auto& geometry = frame.grid().geometry();
  for (std::uint32_t layer = 0; layer < HaloBounds<Det>::kLayerCount; ++layer) {
    const auto& frames = geometry.frames(static_cast<LayerID>(layer));
    for (auto row = rows[layer].begin; row < rows[layer].end; ++row) {
      result.add(layer, frames, {.start = columns.timePS[row], .width = columns.widthPS[row]}, columns.varZ[row], columns.varR[row], columns.varRPhi[row], columns.varLoc0[row], columns.varLoc1[row]);
    }
  }
  return result;
}

/// The bounds over the whole frame: what every slice's stages take, so that
/// no window depends on the slicing. The dip-angle window is a cut at its
/// margin, not only a pre-filter, so a bound taken over a slice's halo would
/// make the candidates depend on which rows the halo held.
template <DetectorLayout Det>
[[nodiscard]] HaloBounds<Det> frameBounds(const TimeFrame<Det>& frame)
{
  std::array<MeasurementRange, HaloBounds<Det>::kLayerCount> layers{};
  for (std::uint32_t layer = 0; layer < layers.size(); ++layer) {
    layers[layer] = frame.grid().layerRows(static_cast<LayerID>(layer));
  }
  return boundsOver(frame, std::span<const MeasurementRange, HaloBounds<Det>::kLayerCount>{layers});
}

/// How far in time a slice's stages reach past its core, in ps, from each
/// link's pair reaches: how much later and
/// how much earlier the box of a pair's target can start than its source's.
///
/// Two overlapping boxes `a` (the source) and `b` (the target) have
/// `b.start - a.start` in `(-width_b, width_a)`: a link reaches its source
/// layer's widest box forward and its target layer's backward. When the two
/// layers have readout frames, the same ones, and the rows of both are exactly
/// those frames, two boxes are equal or disjoint, so a pair lies in one frame
/// and the link reaches 0 both ways. Pairwise overlap is not transitive once
/// boxes differ, so a chain drifts by its links' reaches: `R+` and `R-` are
/// the largest sums of the forward and of the backward reaches along a path
/// through the links.
struct SliceReach {
  std::int64_t pairForwardPS = 0;  ///< the largest link's forward reach
  std::int64_t pairBackwardPS = 0; ///< the largest link's backward reach
  std::int64_t forwardPS = 0;      ///< `R+`
  std::int64_t backwardPS = 0;     ///< `R-`

  /// How far from the core every stage works: `R+ + R-` on both sides. A road
  /// the core owns starts in it and its measurements start within `R-` before
  /// and `R+` after its first; the levels on it depend on the chains into its
  /// measurements, which start within `R+` before and `R-` after them. Within
  /// the gate everything is built, so the output does not depend on the
  /// slicing. With the same widths on every link, `2R`.
  [[nodiscard]] std::int64_t gatePS() const { return forwardPS + backwardPS; }

  /// The time halo a slice needs: the gate, plus one link's reach for the
  /// targets of the tracklets starting at its edges (the backward one before
  /// it, the forward one after; the halo is whole blocks, so the larger).
  [[nodiscard]] std::int64_t haloPS() const { return gatePS() + std::max(pairForwardPS, pairBackwardPS); }
};

/// The reach of the chains one iteration builds on the links \a taken marks,
/// with the triplets' time cut: in a chain `m0 m1 ... mk` of measurements on
/// consecutive links, every pair `m(i) m(i+1)` overlaps (the tracklet) and so
/// does every `m(i) m(i+2)` (`tripletTimeCompatible`, the triplet's first and
/// third). So `m(i+2)` reaches `m(i)` as a pair on the layers two links
/// apart, 0 both ways when those two share frames and their rows are them,
/// whatever the layer between -- a narrow box, a TOF hit's, between two
/// layers of the same frames moves neither.
///
/// Each bound is a sum of link and two-link reaches along the chain. For a
/// chain ending `x -> a -> b`, `b` drifts from `m0` by at most the bound at `a`
/// plus the link's reach and at most the bound at `x` plus the two-link
/// reach: the lesser of the two. Per link, the largest such bound over the
/// chains ending on it (`forwardTo`), and per layer over the chains ending
/// there (`forwardAt`); upper bounds of the same terms, so the lesser of
/// them bounds every chain (by induction along `layerOrder`). `R+` and `R-`
/// are the largest over the links, the forward and the backward bound alike;
/// the pair reaches are the links' own. A chain of two
/// measurements alone -- a tracklet no triplet takes -- is on no road and in
/// no level (a road holds a triplet), so a chain's first link counts in `R+`
/// and `R-` only where a taken link goes on from its target; it always counts
/// in the pair reaches, which size the halo's edge.
template <DetectorLayout Det>
[[nodiscard]] SliceReach chainReach(const HaloBounds<Det>& bounds, const GridGeometry<Det>& geometry, const std::array<bool, nLinks<Det>>& taken)
{
  constexpr auto kLayers = HaloBounds<Det>::kLayerCount;
  const auto zero = [&](LayerID from, LayerID to) {
    return geometry.frames(from) == geometry.frames(to) && bounds.rowsAreFrames[from] && bounds.rowsAreFrames[to];
  };
  const auto forwardOf = [&](LayerID from, LayerID to) -> std::int64_t { return zero(from, to) ? 0 : bounds.widthPS[from]; };
  const auto backwardOf = [&](LayerID from, LayerID to) -> std::int64_t { return zero(from, to) ? 0 : bounds.widthPS[to]; };
  SliceReach result;
  std::array<std::int64_t, nLinks<Det>> forwardTo{};
  std::array<std::int64_t, nLinks<Det>> backwardTo{};
  std::array<std::int64_t, kLayers> forwardAt{};
  std::array<std::int64_t, kLayers> backwardAt{};
  for (const auto layer : layerOrder<Det>) {
    // Every link into this layer is done: its source comes earlier in the order.
    for (auto i = inLinks<Det>.begin(layer); i < inLinks<Det>.end(layer); ++i) {
      if (taken[inLinks<Det>.links[i]]) {
        forwardAt[layer] = std::max(forwardAt[layer], forwardTo[inLinks<Det>.links[i]]);
        backwardAt[layer] = std::max(backwardAt[layer], backwardTo[inLinks<Det>.links[i]]);
      }
    }
    for (auto o = outLinks<Det>.begin(layer); o < outLinks<Det>.end(layer); ++o) {
      const auto out = outLinks<Det>.links[o];
      if (!taken[out]) {
        continue;
      }
      const LayerID target = Det::layerLinks[out].to;
      const std::int64_t forward = forwardOf(layer, target);
      const std::int64_t backward = backwardOf(layer, target);
      result.pairForwardPS = std::max(result.pairForwardPS, forward);
      result.pairBackwardPS = std::max(result.pairBackwardPS, backward);
      // The chains with a measurement before this layer: three or more.
      std::int64_t longerForward = 0;
      std::int64_t longerBackward = 0;
      for (auto i = inLinks<Det>.begin(layer); i < inLinks<Det>.end(layer); ++i) {
        const auto in = inLinks<Det>.links[i];
        if (!taken[in]) {
          continue;
        }
        const LayerID before = Det::layerLinks[in].from;
        longerForward = std::max(longerForward, std::min(forwardTo[in] + forward, forwardAt[before] + forwardOf(before, target)));
        longerBackward = std::max(longerBackward, std::min(backwardTo[in] + backward, backwardAt[before] + backwardOf(before, target)));
      }
      // A chain starting on this layer is this link's pair, the bound of the
      // longer chains it starts.
      forwardTo[out] = std::max(forward, longerForward);
      backwardTo[out] = std::max(backward, longerBackward);
      // Only a chain of three or more measurements -- one with a triplet --
      // is on a road or a level: the pair alone counts where a taken link
      // goes on from its target.
      const bool continues = [&] {
        for (auto n = outLinks<Det>.begin(target); n < outLinks<Det>.end(target); ++n) {
          if (taken[outLinks<Det>.links[n]]) {
            return true;
          }
        }
        return false;
      }();
      result.forwardPS = std::max({result.forwardPS, longerForward, continues ? forward : std::int64_t{0}});
      result.backwardPS = std::max({result.backwardPS, longerBackward, continues ? backward : std::int64_t{0}});
    }
  }
  return result;
}

/// The time a slice's seeding needs past its core: the gate and the halo
/// (`SliceReach::gatePS`, `haloPS`), each the largest of the iterations'.
struct SeedingReach {
  std::int64_t gatePS = 0;
  std::int64_t haloPS = 0;
};

/// The reach the seeding of \a iterations needs: each iteration's along the
/// links it takes (its out-link lists, `takenLinks`), with the triplets' time
/// cut (`chainReach`); of their gates and of their halos the largest. A
/// chain -- tracklet, triplet, neighbour, road, and the chains into a road its
/// levels count -- belongs to one iteration and lies on its
/// taken links, every three consecutive measurements a triplet; so a link no
/// iteration takes, or a path through two iterations' links, widens no gate.
/// The gate is each iteration's `R+ + R-`, and the largest of those -- not
/// the largest `R+` plus the largest `R-`, which two iterations can make
/// without either needing it. The follow does not work in slices.
template <DetectorLayout Det>
[[nodiscard]] SeedingReach seedingReach(const HaloBounds<Det>& bounds, const GridGeometry<Det>& geometry, std::span<const IterationParameters<Det>> iterations)
{
  SeedingReach result;
  for (const auto& parameters : iterations) {
    std::array<bool, nLinks<Det>> taken{};
    for (auto i = parameters.outLinkOffsets[0]; i < parameters.outLinkOffsets[HaloBounds<Det>::kLayerCount]; ++i) {
      taken[parameters.outLinkIndices[i]] = true;
    }
    const auto reach = chainReach(bounds, geometry, taken);
    result.gatePS = std::max(result.gatePS, reach.gatePS());
    result.haloPS = std::max(result.haloPS, reach.haloPS());
  }
  return result;
}

/// One slice of a timeframe: its blocks and rows, its tracklet sources, and
/// the gate's window of box starts around its core.
template <DetectorLayout Det>
struct SlicePlan {
  static constexpr auto kLayerCount = static_cast<std::uint32_t>(Detector<Det>::nLayers);
  using Slice = typename TimeFrame<Det>::TimeSlice;

  Slice slice;
  LayerRowMap<Det> sources;
  std::int64_t gateEarliestPS = 0;
  std::int64_t gateLatestPS = 0;

  /// The tracklet stage's inputs, with the frame's widest candidate boxes per
  /// layer \a maxCandidateWidthPS (`TimeframePlan::layerMaxWidthPS`) and the
  /// gate's core ranges pointing into this plan, which must outlive them.
  [[nodiscard]] TrackletInputs<Det> inputs(const IterationParameters<Det>& parameters, const std::array<std::uint32_t, kLayerCount>& maxCandidateWidthPS) const
  {
    return {
      .parameters = &parameters,
      .sources = &sources,
      .targetBlocks = slice.haloBlocks,
      .maxCandidateWidthPS = maxCandidateWidthPS,
      .gate = {.core = slice.core.data(), .earliestPS = gateEarliestPS, .latestPS = gateLatestPS},
      .halo = slice.halo,
    };
  }
};

/// The plan of \a slice, its gate reaching \a gatePS past the core
/// (`SliceReach::gatePS` of the frame's bounds).
template <DetectorLayout Det>
[[nodiscard]] SlicePlan<Det> planSlice(const TimeFrame<Det>& frame, const typename TimeFrame<Det>::TimeSlice& slice, std::int64_t gatePS)
{
  SlicePlan<Det> plan{.slice = slice, .sources = TimeFrameScratch<Det>::makeSources(slice.halo)};
  // The core's box starts, widened by the gate: a halo source beyond it lies
  // on no chain the core's output depends on. The reach is the frame's rows',
  // as the halo's; the core's edges are rounded outward to picoseconds.
  const auto& geometry = frame.grid().geometry();
  const std::int64_t parts = geometry.partsPerPS();
  const auto floorDiv = [](std::int64_t a, std::int64_t b) { return (a / b) - ((a % b) < 0 ? 1 : 0); };
  plan.gateEarliestPS = floorDiv(geometry.blockStart(slice.coreBlocks.begin), parts) - gatePS;
  plan.gateLatestPS = -floorDiv(-geometry.blockStart(slice.coreBlocks.end), parts) + gatePS;
  return plan;
}

/// A timeframe's plan for every iteration of a run -- its slices, the same
/// for every iteration, and each iteration's vertex windows -- with what it
/// was planned from.
template <DetectorLayout Det>
struct TimeframePlan {
  static constexpr auto kLayerCount = static_cast<std::uint32_t>(Detector<Det>::nLayers);

  // What it was planned from.
  std::uint32_t rows = 0U;
  std::int32_t blockCount = 0;
  std::int64_t blockWidth = 0; ///< [parts of a ps]
  std::vector<IterationParameters<Det>> parameters;

  VertexSeed beamLine;
  std::vector<VertexSeed> vertices; ///< the timeframe's vertex list, sorted by time
  /// The widest candidate box per layer, over the whole frame (`frameBounds`).
  std::array<std::uint32_t, kLayerCount> layerMaxWidthPS{};
  /// The triplet stage's target bounds, over the whole frame (`frameBounds`).
  std::array<TripletTargetBounds, kLayerCount> layerTargetBounds{};
  /// The largest local measurement variance per layer, over the whole frame
  /// (`frameBounds`).
  std::array<Scalar, kLayerCount> layerMaxLocVariance{};
  std::int64_t haloPS = 0;
  std::vector<SlicePlan<Det>> slices;                     ///< in order
  std::vector<VertexWindowParameters<Det>> vertexWindows; ///< per iteration

  /// Throw std::invalid_argument unless this plan was planned from \a frame
  /// -- its rows and grid -- and from \a configured. The parameters compare
  /// as bytes: every copy of them comes from the same objects by trivial
  /// copies, padding included.
  void check(const TimeFrame<Det>& frame, std::span<const IterationParameters<Det>> configured) const
  {
    const auto& geometry = frame.grid().geometry();
    if (rows != frame.size() || blockCount != geometry.blockCount() || blockWidth != geometry.blockWidth()) {
      throw std::invalid_argument("a timeframe plan of another frame");
    }
    if (configured.size() != parameters.size() || vertexWindows.size() != parameters.size()) {
      throw std::invalid_argument("a timeframe plan of another number of iterations");
    }
    for (std::uint32_t i = 0; i < parameters.size(); ++i) {
      // Padding included, deliberately: both sides are copies of the same
      // objects, made by std::copy or assign, which move their bytes as they are.
      // NOLINTNEXTLINE(bugprone-suspicious-memory-comparison)
      if (std::memcmp(&configured[i], &parameters[i], sizeof(IterationParameters<Det>)) != 0) {
        throw std::invalid_argument("a timeframe plan of other iteration parameters");
      }
    }
  }
};

/// \a plan with the vertex list \a vertices, sorted by time, and each
/// iteration's vertex windows from it and \a windowSettings, the settings its
/// parameters were derived from. The slices and the frame's bounds stay as
/// they are: they depend on the frame and iterations, so a plan made before the
/// vertexer's list was known is replanned with it without scanning the frame
/// again.
template <DetectorLayout Det>
[[nodiscard]] TimeframePlan<Det> planWithVertices(TimeframePlan<Det> plan, std::span<const SearchWindowParams> windowSettings, std::span<const VertexSeed> vertices)
{
  if (windowSettings.size() != plan.parameters.size()) {
    throw std::invalid_argument("planWithVertices needs one window setting per iteration");
  }
  plan.vertices.assign(vertices.begin(), vertices.end());
  plan.vertexWindows.clear();
  plan.vertexWindows.reserve(plan.parameters.size());
  for (std::uint32_t iteration = 0; iteration < plan.parameters.size(); ++iteration) {
    plan.vertexWindows.push_back(vertexWindowParameters<Det>(plan.parameters[iteration], windowSettings[iteration], vertices));
  }
  return plan;
}

/// Plan every iteration of \a iterations over \a frame in slices of
/// \a sliceWidthBlocks grid blocks (`GridGeometry::blocksPerSlice` of the
/// configured width) with the halo the iterations' chains need
/// (`seedingReach`). The vertexing pass uses the same planner, each
/// slice finding its core's frames' vertices (`OwnedFrames`). \a bounds are
/// the frame's (`frameBounds`, as the measurement producer ships them); \a
/// windowSettings are the settings each iteration's parameters were derived
/// from; \a beamLine the seed a beam-line iteration reads; \a vertices the
/// timeframe's vertex list, sorted by time.
template <DetectorLayout Det>
[[nodiscard]] TimeframePlan<Det> planTimeframe(const TimeFrame<Det>& frame, const HaloBounds<Det>& bounds, std::span<const IterationParameters<Det>> iterations, std::span<const SearchWindowParams> windowSettings, const VertexSeed& beamLine, std::span<const VertexSeed> vertices, std::int32_t sliceWidthBlocks)
{
  if (windowSettings.size() != iterations.size()) {
    throw std::invalid_argument("planTimeframe needs one window setting per iteration");
  }
  TimeframePlan<Det> plan;
  plan.rows = frame.size();
  plan.blockCount = frame.grid().geometry().blockCount();
  plan.blockWidth = frame.grid().geometry().blockWidth();
  plan.parameters.assign(iterations.begin(), iterations.end());
  plan.beamLine = beamLine;
  plan.vertices.assign(vertices.begin(), vertices.end());
  plan.layerMaxWidthPS = bounds.widthPS;
  plan.layerTargetBounds = bounds.target;
  plan.layerMaxLocVariance = bounds.varLoc;
  const auto reach = seedingReach(bounds, frame.grid().geometry(), iterations);
  plan.haloPS = reach.haloPS;
  const auto width = std::max<std::int32_t>(1, sliceWidthBlocks);
  for (std::int32_t begin = 0; begin < plan.blockCount; begin += width) {
    const typename TimeFrame<Det>::TimeSlice slice{frame, {.begin = begin, .end = std::min(plan.blockCount, begin + width)}, plan.haloPS};
    plan.slices.push_back(planSlice(frame, slice, reach.gatePS));
  }
  return planWithVertices<Det>(std::move(plan), windowSettings, vertices);
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_SLICEPLAN_H_
