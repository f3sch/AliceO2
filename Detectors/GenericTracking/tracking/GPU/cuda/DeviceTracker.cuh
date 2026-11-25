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

/// \file DeviceTracker.cuh
/// \brief DeviceTracker<Det> on one backend: a timeframe's slices,
///        iterations and tracks.
///
/// What lives where in the arena:
///
///   run        the iterations' parameters, the geometry, the material, until
///              configureRun again or detach; the arena itself, over the
///              attached block
///   timeframe  base end: the measurement columns, the grid offsets, the used
///              mask, the TimeFrame over them and the plan's device part,
///              until the next beginTimeframe
///   tracks     scratch end, at its bottom: every iteration's tracks, one
///              chunk per iteration, until the next beginTimeframe
///   chunks     scratch end: the open iteration's records, one chunk per
///              batch of roads fitted, until endIteration has gathered them
///   batch      scratch end, inside the chunks: the captured roads not yet
///              fitted, one chunk per slice, until they fill the device or
///              the iteration ends
///   slice      base end: the slice's outputs, until the next runSlice or
///              endIteration
///   iteration  base end, inside endIteration: the gathered records, the
///              selection's arrays, the working mask, the iteration's tracks
///              and proposals
///
/// Staging, scan and sort storage are on the scratch end, released when the
/// pass that took them returns. Included by the .cu that instantiates it for
/// a detector, in the backend's namespace, like GpuStage.cuh.

#ifndef GTRK_TRACKING_GPU_DEVICETRACKER_CUH_
#define GTRK_TRACKING_GPU_DEVICETRACKER_CUH_

#include "DeviceAlgorithms.cuh"
#include "FlatStage.cuh"
#include "GpuStage.cuh"
#include "PathGlue.cuh"
#include "Platform.cuh"
#include "RoadGlue.cuh"
#include "SelectionGlue.cuh"
#include "SliceGlue.cuh"
#include "VertexGlue.cuh"

#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/DeviceStageOutput.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Tracking/GPU/DeviceTracker.h"
#include "GenericTracking/Tracking/InactiveRows.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/NeighbourStage.h"
#include "GenericTracking/Tracking/OverlapPartners.h"
#include "GenericTracking/Tracking/RoadCount.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/SensorCells.h"
#include "GenericTracking/Tracking/SlicePlan.h"
#include "GenericTracking/Tracking/TrackJoin.h"
#include "GenericTracking/Tracking/TrackRefit.h"
#include "GenericTracking/Tracking/TrackletStage.h"
#include "GenericTracking/Tracking/TripletStage.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE
{

namespace tracker_detail
{

/// The lanes that follow one track together, on every iteration: a step's
/// window holds many
/// candidates, each an update, and the lookahead confirms every row passing
/// by the next step's search -- long windows a group shares. 16 on both
/// backends, for warps of 32 and waves of 64 alike.
inline constexpr std::uint32_t kFollowLanes = 16U;

/// \a stage's follow of \a count tracks: `kFollowLanes` lanes in step a
/// track (`FollowStage::together`).
template <typename Stage>
void follow(DeviceRuntime& runtime, DeviceStream& stream, std::uint32_t count, const Stage& stage)
{
  algorithms::forEachGroup<kFollowLanes>(runtime, stream, count, stage);
}

/// The device's launches of the follow's passes (`follow::runSequence`) on
/// \a stream: a pass one thread a track, the proposals by `follow`, the bids
/// zeroed by a fill.
struct FollowLaunch {
  DeviceRuntime& runtime;
  DeviceStream& stream;

  template <typename Body>
  void pass(profile::Phase /*phase*/, std::uint32_t count, const Body& body) const
  {
    algorithms::forEach(runtime, stream, count, body);
  }

  template <typename Stage>
  void propose(profile::Phase /*phase*/, std::uint32_t count, const Stage& stage) const
  {
    tracker_detail::follow(runtime, stream, count, stage);
  }

  void clearBids(std::uint32_t* bids, std::size_t rows) const { runtime.fill(bids, 0U, rows * sizeof(std::uint32_t), stream); }
};

/// Whether a source has a query (`TrackletStage::hasQuery`), for the
/// tracklet launch's list.
template <typename Stage>
struct HasQuery {
  Stage stage;
  GPUdi() bool operator()(std::uint32_t source) const { return stage.hasQuery(source); }
};

/// The tracks whose refit succeeded and that the iteration keeps after the
/// follow (`keptAfterFollow`).
template <DetectorLayout Det>
struct KeptAfterFollow {
  const IterationParameters<Det>* parameters = nullptr;
  const std::uint8_t* count = nullptr;
  const BoundTrackState* state = nullptr;
  const PropagationStatus* status = nullptr;
  GPUdi() bool operator()(std::uint32_t track) const { return status[track] == PropagationStatus::kOk && keptAfterFollow(*parameters, count[track], state[track]); }
};

/// The kept tracks' rows in \a mask as \a after has them: the mask before the
/// selection plus what the kept tracks hold, as the host's `dropTracks`
/// rebuilds it. Two kept tracks sharing a row store the same byte.
template <std::size_t K>
struct RestoreRows {
  const ClusterList<K>* clusters = nullptr;
  const std::uint8_t* count = nullptr;
  const std::uint8_t* after = nullptr;
  std::uint8_t* mask = nullptr;
  GPUdi() void operator()(std::uint32_t track) const
  {
    for (std::uint32_t position = 0; position < count[track]; ++position) {
      const auto row = clusters[track].index[position];
      mask[row] = after[row];
    }
  }
};

/// One value copied.
struct CopyAt {
  const std::uint32_t* from = nullptr;
  std::uint32_t* to = nullptr;
  GPUdi() void operator()(std::uint32_t /*at*/) const { *to = *from; }
};

/// Every count 0.
struct ZeroCounts {
  std::uint32_t* counts = nullptr;
  GPUdi() void operator()(std::uint32_t at) const { counts[at] = 0U; }
};

/// A proposal's sort key, its row then its chi2 (not negative, so its bits
/// order it), and its index.
struct PartnerKeys {
  const std::uint32_t* row = nullptr;
  const Scalar* chi2 = nullptr;
  std::uint64_t* keys = nullptr;
  std::uint32_t* order = nullptr;
  GPUdi() void operator()(std::uint32_t at) const
  {
    keys[at] = (std::uint64_t{row[at]} << 32U) | std::uint64_t{std::bit_cast<std::uint32_t>(chi2[at])};
    order[at] = at;
  }
};

/// Whether the sorted entry at \a at is the first of its key's high half
/// (a proposal's row, a choice's row or partner): a flag per entry, by its
/// index.
struct FirstOfRow {
  const std::uint64_t* keys = nullptr;
  const std::uint32_t* sorted = nullptr;
  std::uint8_t* won = nullptr;
  GPUdi() void operator()(std::uint32_t at) const { won[sorted[at]] = (at == 0U || (keys[at] >> 32U) != (keys[at - 1U] >> 32U)) ? 1U : 0U; }
};

struct Flagged {
  const std::uint8_t* flags = nullptr;
  GPUdi() bool operator()(std::uint32_t at) const { return flags[at] != 0U; }
};

/// Per track \a at, where its winners begin among the winners ordered by
/// track: the first whose track is not below it.
struct LowerBounds {
  const std::uint32_t* tracks = nullptr;
  std::uint32_t count = 0U;
  std::uint32_t* offsets = nullptr;
  GPUdi() void operator()(std::uint32_t at) const
  {
    std::uint32_t low = 0U;
    std::uint32_t high = count;
    while (low < high) {
      const auto middle = low + ((high - low) / 2U);
      if (tracks[middle] < at) {
        low = middle + 1U;
      } else {
        high = middle;
      }
    }
    offsets[at] = low;
  }
};

/// A track index as columns: the candidates of the next iteration's overlap
/// partners the follow of earlier tracks changed.
struct TrackIndexColumns {
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
  template <template <class> class C>
  struct RecordT {
    MEMLAYOUT_APPLY(RecordT, track)
    C<std::uint32_t> track;
  };
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
};

/// The iteration's own tracks' global indices: \a earlier on.
struct OwnTracks {
  std::uint32_t earlier = 0U;
  std::uint32_t* values = nullptr;
  GPUdi() void operator()(std::uint32_t at) const { values[at] = earlier + at; }
};

/// Every index its own value.
struct Iota {
  std::uint32_t* values = nullptr;
  GPUdi() void operator()(std::uint32_t at) const { values[at] = at; }
};

/// Whether the sorted value at \a at is the first of its run: a flag per
/// entry.
struct FirstOfValue {
  const std::uint32_t* values = nullptr;
  std::uint8_t* first = nullptr;
  GPUdi() void operator()(std::uint32_t at) const { first[at] = (at == 0U || values[at] != values[at - 1U]) ? 1U : 0U; }
};

/// Whether the follow changed a track: a row proposed either way.
template <typename Proposal>
struct FollowChanged {
  const Proposal* proposals = nullptr;
  GPUdi() bool operator()(std::uint32_t at) const { return proposals[at].outwardCount + proposals[at].inwardCount != 0U; }
};

/// Per candidate \a at, its range of the claims ordered by track: the first
/// claim of its track and the first of the next.
struct HeldRange {
  const std::uint32_t* candidates = nullptr;
  const std::uint32_t* claimTracks = nullptr; ///< ascending
  std::uint32_t claimCount = 0U;
  std::uint32_t* begin = nullptr;
  std::uint32_t* end = nullptr;
  GPUdi() void operator()(std::uint32_t at) const
  {
    const auto bound = [&](std::uint32_t track) {
      std::uint32_t low = 0U;
      std::uint32_t high = claimCount;
      while (low < high) {
        const auto middle = low + ((high - low) / 2U);
        if (claimTracks[middle] < track) {
          low = middle + 1U;
        } else {
          high = middle;
        }
      }
      return low;
    };
    begin[at] = bound(candidates[at]);
    end[at] = bound(candidates[at] + 1U);
  }
};

/// A candidate's partners to fit: those it holds and its winners, where it
/// has winners; none else.
struct EntryCounts {
  const std::uint32_t* winnerOffsets = nullptr;
  const std::uint32_t* heldBegin = nullptr;
  const std::uint32_t* heldEnd = nullptr;
  std::uint32_t* counts = nullptr;
  GPUdi() void operator()(std::uint32_t at) const
  {
    const auto winners = winnerOffsets[at + 1U] - winnerOffsets[at];
    counts[at] = winners != 0U ? (heldEnd[at] - heldBegin[at]) + winners : 0U;
  }
};

/// A candidate's partners to fit into its place: those it holds, then its
/// winners, flagged fresh, in their order of proposal.
struct FillEntries {
  const std::uint32_t* winnerOffsets = nullptr;
  const std::uint32_t* winners = nullptr;
  const std::uint32_t* proposedRow = nullptr;
  const Scalar* proposedPath = nullptr;
  const std::uint32_t* heldBegin = nullptr;
  const std::uint32_t* heldEnd = nullptr;
  const std::uint32_t* heldRows = nullptr;
  const Scalar* heldPaths = nullptr;
  const std::uint32_t* offsets = nullptr;
  std::uint32_t* rows = nullptr;
  Scalar* paths = nullptr;
  std::uint8_t* fresh = nullptr;
  GPUdi() void operator()(std::uint32_t at) const
  {
    if (winnerOffsets[at + 1U] == winnerOffsets[at]) {
      return;
    }
    auto to = offsets[at];
    for (auto from = heldBegin[at]; from < heldEnd[at]; ++from, ++to) {
      rows[to] = heldRows[from];
      paths[to] = heldPaths[from];
      fresh[to] = 0U;
    }
    for (auto from = winnerOffsets[at]; from < winnerOffsets[at + 1U]; ++from, ++to) {
      rows[to] = proposedRow[winners[from]];
      paths[to] = proposedPath[winners[from]];
      fresh[to] = 1U;
    }
  }
};

/// A candidate's new claims: its winners where its fit with them was kept.
struct ClaimCounts {
  const std::uint32_t* winnerOffsets = nullptr;
  const std::uint8_t* kept = nullptr;
  std::uint32_t* counts = nullptr;
  GPUdi() void operator()(std::uint32_t at) const { counts[at] = kept[at] != 0U ? winnerOffsets[at + 1U] - winnerOffsets[at] : 0U; }
};

/// A kept candidate's fresh partners claimed, from \a claimOffsets[at] on,
/// each taken in the used mask.
struct PlaceClaims {
  const std::uint32_t* candidates = nullptr;
  const std::uint8_t* kept = nullptr;
  const std::uint32_t* offsets = nullptr;
  const std::uint32_t* rows = nullptr;
  const Scalar* paths = nullptr;
  const std::uint8_t* fresh = nullptr;
  const std::uint32_t* claimOffsets = nullptr;
  DeviceColumns<PartnerClaimColumns::template RecordT>::View claims;
  std::uint8_t* used = nullptr;
  GPUdi() void operator()(std::uint32_t at) const
  {
    if (kept[at] == 0U) {
      return;
    }
    auto to = claimOffsets[at];
    for (auto from = offsets[at]; from < offsets[at + 1U]; ++from) {
      if (fresh[from] == 0U) {
        continue;
      }
      claims.track[to] = candidates[at];
      claims.row[to] = rows[from];
      claims.path[to] = paths[from];
      ++to;
      used[rows[from]] = kExclusive;
    }
  }
};

/// After the last iteration, per track with partners \a at: the first of its
/// claims (ordered by track) at \a starts[at], the track and its range.
struct ListedRanges {
  const std::uint32_t* starts = nullptr;
  const std::uint32_t* claimTracks = nullptr;
  std::uint32_t listedCount = 0U;
  std::uint32_t claimCount = 0U;
  std::uint32_t* listed = nullptr;
  std::uint32_t* offsets = nullptr;
  GPUdi() void operator()(std::uint32_t at) const
  {
    listed[at] = claimTracks[starts[at]];
    offsets[at] = starts[at];
    if (at + 1U == listedCount) {
      offsets[listedCount] = claimCount;
    }
  }
};

/// The rows of the partners of a track whose fit with them failed, freed.
struct ReleaseRows {
  const std::uint8_t* kept = nullptr;
  const std::uint32_t* offsets = nullptr;
  const std::uint32_t* rows = nullptr;
  std::uint8_t* used = nullptr;
  GPUdi() void operator()(std::uint32_t at) const
  {
    if (kept[at] != 0U) {
      return;
    }
    for (auto from = offsets[at]; from < offsets[at + 1U]; ++from) {
      used[rows[from]] = kFree;
    }
  }
};

/// The tracklet stage as a flat run:
/// the walk per source is `TrackletStage::forEachCandidate`, the contexts
/// its source's (`SourceContext`, the vertex cache, `BetweenSource`) and each
/// query's (`QueryContext`, the vertex intervals), the test `testCandidate`,
/// the record `record`. A pair's key is the source's slot in its chunk (24
/// bits) and its query's position among the source's out-links (8).
template <typename Stage>
struct FlatTracklets {
  using SourceContext = typename Stage::SourceContext;
  using QueryContext = typename Stage::QueryContext;
  using CompatibleVertex = typename Stage::CompatibleVertex;
  using Between = typename Stage::Between;
  static constexpr std::uint32_t kVertexCapacity = Stage::kVertexCapacity;
  static constexpr std::uint32_t kQueryBits = 8U;
  /// The lanes that walk one source together in the fill: one thread a
  /// source, or a group of lanes (`WarpGroup`).
  static constexpr std::uint32_t kFillLanes = gpu::kFillLanes;
  struct Source {
    SourceContext context;
    CompatibleVertex vertices[kVertexCapacity];
    Between between;
  };
  struct Query {
    QueryContext context;
    SecondRange intervals[kVertexCapacity];
    bool hasIntervals = false;
  };
  /// A fill's bookkeeping over one source's walk (`note`).
  struct FillState {
    bool stored = false;
    std::uint32_t query = 0U;
    std::size_t lastLink = 0U;
  };
  Stage stage;
  std::uint32_t queriesPerSource = 1U; ///< the most out-links of a layer
  Source* sources = nullptr;           ///< per chunk element
  Query* queries = nullptr;            ///< per chunk element, `queriesPerSource` each
  std::uint32_t sourceBegin = 0U;

  [[nodiscard]] GPUdi() std::uint32_t count(std::uint32_t element) const { return stage.candidateCount(element); }
  template <typename Sink>
  GPUdi() void fill(std::uint32_t element, std::uint32_t position, Sink& sink) const
  {
    const std::uint32_t slot = position - sourceBegin;
    FillState state;
    const auto one = [&](const SourceContext& source, const CompatibleVertex* cachedVertices, const Between& between, const QueryContext& context, const SecondRange* intervals, const MeasurementView&, std::uint32_t targetRow) { sink.put(targetRow, note(state, slot, true, source, cachedVertices, between, context, intervals)); };
    stage.forEachCandidate(element, one);
  }
  /// `fill` by the lanes of \a group in step (`forEachCandidate`): a step's
  /// candidates in lane order at the sink's next positions, so the list is
  /// the serial walk's. The bookkeeping on the step's ballot, alike on every
  /// lane; lane 0 stores.
  template <typename Sink, typename Group>
  GPUdi() void fill(std::uint32_t element, std::uint32_t position, Sink& sink, const Group& group) const
  {
    const std::uint32_t slot = position - sourceBegin;
    const std::uint32_t lane = group.lane();
    FillState state;
    const auto step = [&](const SourceContext& source, const CompatibleVertex* cachedVertices, const Between& between, const QueryContext& context, const SecondRange* intervals, const MeasurementView&, std::uint32_t targetRow, bool has) {
      const std::uint64_t mask = group.ballot(has);
      if (mask == 0U) {
        return;
      }
      sink.putStep(mask, lane, has, targetRow, note(state, slot, lane == 0U, source, cachedVertices, between, context, intervals));
    };
    stage.forEachCandidate(element, step, group);
  }
  /// A candidate's key, its source and query noted in \a state: the source's
  /// context stored at its first candidate, a query's at its first; a
  /// query's number counts the queries with candidates before it. \a store
  /// whether this thread writes the contexts.
  [[nodiscard]] GPUdi() std::uint32_t note(FillState& state, std::uint32_t slot, bool store, const SourceContext& source, const CompatibleVertex* cachedVertices, const Between& between, const QueryContext& context, const SecondRange* intervals) const
  {
    if (!state.stored) {
      if (store) {
        sources[slot].context = source;
        for (std::uint32_t i = 0; i < source.cachedCount; ++i) {
          sources[slot].vertices[i] = cachedVertices[i];
        }
        sources[slot].between = between;
        keep(slot, 0U, context, intervals);
      }
      state.stored = true;
      state.lastLink = context.linkIndex;
      state.query = 0U;
    } else if (context.linkIndex != state.lastLink) {
      state.lastLink = context.linkIndex;
      ++state.query;
      if (store) {
        keep(slot, state.query, context, intervals);
      }
    }
    return (slot << kQueryBits) | state.query;
  }
  GPUdi() void keep(std::uint32_t slot, std::uint32_t query, const QueryContext& context, const SecondRange* intervals) const
  {
    Query& out = queries[(std::size_t{slot} * queriesPerSource) + query];
    out.context = context;
    out.hasIntervals = intervals != nullptr;
    if (intervals != nullptr) {
      for (std::uint32_t i = 0; i < kVertexCapacity; ++i) {
        out.intervals[i] = intervals[i];
      }
    }
  }
  [[nodiscard]] GPUdi() std::uint32_t sourceOf(const flat::CandidatePair& pair) const { return (pair.key >> kQueryBits) + sourceBegin; }
  [[nodiscard]] GPUdi() bool test(std::uint32_t, const flat::CandidatePair& pair) const
  {
    const Source& source = sources[pair.key >> kQueryBits];
    const Query& query = queries[(std::size_t{pair.key >> kQueryBits} * queriesPerSource) + (pair.key & ((1U << kQueryBits) - 1U))];
    const MeasurementView columns = stage.columns();
    const PairPoint target = pairPointAt(columns, pair.row);
    return stage.testCandidate(source.context, source.vertices, source.between, query.context, query.hasIntervals ? query.intervals : nullptr, columns, pair.row, target);
  }
  [[nodiscard]] GPUdi() Tracklet emit(std::uint32_t, const flat::CandidatePair& pair, std::uint32_t) const
  {
    const Source& source = sources[pair.key >> kQueryBits];
    const Query& query = queries[(std::size_t{pair.key >> kQueryBits} * queriesPerSource) + (pair.key & ((1U << kQueryBits) - 1U))];
    const PairPoint target = pairPointAt(stage.columns(), pair.row);
    return stage.record(source.context.at.row, pair.row, query.context.linkIndex, source.context.point, target);
  }
  [[nodiscard]] bool beginChunk(DeviceArena& arena, std::uint32_t begin, std::uint32_t count, std::uint32_t)
  {
    if (count > (std::uint32_t{1U} << (32U - kQueryBits))) {
      return false; // more sources in the chunk than the key numbers: the run fails for memory, which halves the chunk
    }
    sourceBegin = begin;
    DeviceSpan<Source> sourceStorage;
    DeviceSpan<Query> queryStorage;
    if (!sourceStorage.allocate(arena, std::max<std::size_t>(count, 1U), ArenaEnd::kScratch) || !queryStorage.allocate(arena, std::max<std::size_t>(std::size_t{count} * queriesPerSource, 1U), ArenaEnd::kScratch)) {
      return false;
    }
    sources = sourceStorage.data();
    queries = queryStorage.data();
    return true;
  }
  [[nodiscard]] std::size_t contextBytes(std::uint32_t count) const { return std::size_t{count} * (sizeof(Source) + (std::size_t{queriesPerSource} * sizeof(Query))); }
  [[nodiscard]] std::size_t pairBytes() const { return 0U; }
};

/// The triplet stage as a flat run:
/// the walk per first tracklet is `TripletStage::prepareFirst` and
/// `forEachWindowCandidate`, the context its `FirstContext`, the test
/// `testCandidate` with its fit kept per pair for `record`; with \a Shape,
/// the shape-testing stage.
template <DetectorLayout Det, bool Shape = false>
struct FlatTriplets {
  using Stage = TripletStage<Det, Shape>;
  using Context = typename Stage::FirstContext;
  Stage stage;
  Context* contexts = nullptr;      ///< per chunk element
  TripletFitResult* fits = nullptr; ///< per chunk pair, where the test passed
  std::uint8_t* kinked = nullptr;
  std::uint32_t sourceBegin = 0U;

  [[nodiscard]] GPUdi() std::uint32_t count(std::uint32_t element) const { return stage.candidateCount(element); }
  template <typename Sink>
  GPUdi() void fill(std::uint32_t element, std::uint32_t position, Sink& sink) const
  {
    Context context;
    if (!stage.prepareFirst(element, context)) {
      return; // no candidates: the context is never read
    }
    contexts[position - sourceBegin] = context;
    const auto one = [&](std::uint32_t candidateIndex) { sink.put(candidateIndex, position); };
    stage.forEachWindowCandidate(context, one);
  }
  [[nodiscard]] GPUdi() std::uint32_t sourceOf(const flat::CandidatePair& pair) const { return pair.key; }
  [[nodiscard]] GPUdi() bool test(std::uint32_t k, const flat::CandidatePair& pair) const
  {
    TripletFitResult fit;
    bool isKinked = false;
    if (!stage.testCandidate(contexts[pair.key - sourceBegin], pair.row, fit, isKinked)) {
      return false;
    }
    fits[k] = fit;
    kinked[k] = isKinked ? std::uint8_t{1} : std::uint8_t{0};
    return true;
  }
  [[nodiscard]] GPUdi() Triplet emit(std::uint32_t k, const flat::CandidatePair& pair, std::uint32_t) const { return stage.record(contexts[pair.key - sourceBegin], pair.row, fits[k], kinked[k] != 0U); }
  [[nodiscard]] bool beginChunk(DeviceArena& arena, std::uint32_t begin, std::uint32_t sources, std::uint32_t pairs)
  {
    sourceBegin = begin;
    DeviceSpan<Context> contextStorage;
    DeviceSpan<TripletFitResult> fitStorage;
    DeviceSpan<std::uint8_t> kinkedStorage;
    if (!contextStorage.allocate(arena, std::max<std::size_t>(sources, 1U), ArenaEnd::kScratch) || !fitStorage.allocate(arena, std::max<std::size_t>(pairs, 1U), ArenaEnd::kScratch) || !kinkedStorage.allocate(arena, std::max<std::size_t>(pairs, 1U), ArenaEnd::kScratch)) {
      return false;
    }
    contexts = contextStorage.data();
    fits = fitStorage.data();
    kinked = kinkedStorage.data();
    return true;
  }
  [[nodiscard]] std::size_t contextBytes(std::uint32_t sources) const { return std::size_t{sources} * sizeof(Context); }
  [[nodiscard]] std::size_t pairBytes() const { return sizeof(TripletFitResult) + 1U; }
};

/// A slice's roads as the iteration's fit takes them (`RoadRows`): each
/// road's rows and their count.
template <std::size_t K>
struct CapturedRoadOf {
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
  template <template <class> class C>
  struct RecordT {
    MEMLAYOUT_APPLY(RecordT, rows, count)
    C<ClusterList<K>> rows;
    C<std::uint8_t> count;
  };
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
};

/// One road of a slice captured into \a roads, by its index.
template <DetectorLayout Det>
struct CaptureRoad {
  RoadRows<Det> rows;
  typename DeviceColumns<CapturedRoadOf<maxChainLayers<Det>>::template RecordT>::View roads;

  GPUdi() void operator()(std::uint32_t road) const
  {
    roads.count[road] = static_cast<std::uint8_t>(rows.capture(road, roads.rows[road].index));
  }
};

/// The fit of the captured roads, one element per road: \a fit
/// (`RoadRowsFit`) on its rows where they are.
template <typename Fit, std::size_t K>
struct CapturedRoadFit {
  Fit fit;
  typename DeviceColumns<CapturedRoadOf<K>::template RecordT>::ConstView roads;

  template <typename Emitter>
  GPUdi() void operator()(std::uint32_t road, Emitter& emitter) const
  {
    fit(roads.rows[road].index, roads.count[road], emitter);
  }
};

/// The fit's order of the captured roads, a key and the index per road:
/// the road's layer pattern in the high word, its first row's azimuth in the
/// low word, so that the roads a warp fits in turn have the same layers and
/// neighbouring material cells.
template <DetectorLayout Det, std::size_t K>
struct RoadOrderKey {
  const TrackingGeometry<Det>* geometry;
  const TimeFrame<Det>* frame;
  typename DeviceColumns<CapturedRoadOf<K>::template RecordT>::ConstView roads;
  std::uint64_t* keys;
  std::uint32_t* indices;

  GPUdi() void operator()(std::uint32_t road) const
  {
    const std::uint32_t* rows = roads.rows[road].index;
    const auto columns = frame->measurements();
    std::uint32_t pattern = 0U;
    for (std::uint32_t i = 0; i < roads.count[road]; ++i) {
      pattern |= 1U << geometry->surface(columns.surface[rows[i]]).layer;
    }
    // The azimuth's float bits ordered as unsigned integers: a negative
    // value's bits inverted, a positive value's sign bit set.
    const std::uint32_t bits = o2::gpu::CAMath::Float2UIntReint(atan2Of(columns.y[rows[0]], columns.x[rows[0]]));
    const std::uint32_t phi = (bits & 0x80000000U) != 0U ? ~bits : bits | 0x80000000U;
    keys[road] = (static_cast<std::uint64_t>(pattern) << 32U) | phi;
    indices[road] = road;
  }
};

/// \a x's bits as an unsigned key that sorts as the float does (-0 as +0;
/// a NaN is never a key here).
[[nodiscard]] GPUdi() std::uint32_t orderedBits(Scalar x)
{
  const auto bits = std::bit_cast<std::uint32_t>(x == 0.F ? 0.F : x);
  return (bits & 0x80000000U) != 0U ? ~bits : (bits | 0x80000000U);
}

/// Whether a track has rows and its fit stands: the join's earlier tracks.
struct HasRows {
  const PropagationStatus* status = nullptr;
  const std::uint8_t* count = nullptr;
  GPUdi() bool operator()(std::uint32_t track) const { return status[track] == PropagationStatus::kOk && count[track] != 0U; }
};

/// The listed tracks' outermost and innermost rows.
template <std::size_t K>
struct EndRows {
  const ClusterList<K>* clusters = nullptr;
  const std::uint8_t* count = nullptr;
  const std::uint32_t* listed = nullptr;
  std::uint32_t* ends = nullptr;
  std::uint32_t* starts = nullptr;
  GPUdi() void operator()(std::uint32_t at) const
  {
    const auto track = listed[at];
    ends[at] = clusters[track].index[count[track] - 1U];
    starts[at] = clusters[track].index[0];
  }
};

/// (row, track) pairs from the sorted rows and their tracks.
struct ZipEnds {
  const std::uint32_t* rows = nullptr;
  const std::uint32_t* tracks = nullptr;
  join::End* out = nullptr;
  GPUdi() void operator()(std::uint32_t at) const { out[at] = {.row = rows[at], .track = tracks[at]}; }
};

/// Whether a piece chose a partner.
struct ChosePartner {
  const join::Choice* choices = nullptr;
  GPUdi() bool operator()(std::uint32_t at) const { return choices[at].partner != join::kNone; }
};

/// A listed piece's sort key, its partner then its score, and the piece.
struct JoinKeys {
  const join::Choice* choices = nullptr;
  const std::uint32_t* listed = nullptr;
  std::uint64_t* keys = nullptr;
  std::uint32_t* pieces = nullptr;
  GPUdi() void operator()(std::uint32_t at) const
  {
    const auto piece = listed[at];
    keys[at] = (std::uint64_t{choices[piece].partner} << 32U) | std::uint64_t{orderedBits(choices[piece].score)};
    pieces[at] = piece;
  }
};

/// Whether an element is not flagged.
struct Unflagged {
  const std::uint8_t* flags = nullptr;
  GPUdi() bool operator()(std::uint32_t at) const { return flags[at] == 0U; }
};

/// The joined pieces' partners.
struct PartnersOf {
  const join::Choice* choices = nullptr;
  const std::uint32_t* joined = nullptr;
  std::uint32_t* partners = nullptr;
  GPUdi() void operator()(std::uint32_t at) const { partners[at] = choices[joined[at]].partner; }
};

/// out[at] = first + at.
struct Sequence {
  std::uint32_t first = 0U;
  std::uint32_t* out = nullptr;
  GPUdi() void operator()(std::uint32_t at) const { out[at] = first + at; }
};

/// Whether a candidate chose an inactive layer's row.
struct ChoseRow {
  const inactive::Choice* choices = nullptr;
  GPUdi() bool operator()(std::uint32_t at) const { return choices[at].row != inactive::kNone; }
};

/// A listed candidate's sort key, its row then its rise, and the candidate.
struct InactiveKeys {
  const inactive::Choice* choices = nullptr;
  const std::uint32_t* listed = nullptr;
  std::uint64_t* keys = nullptr;
  std::uint32_t* positions = nullptr;
  GPUdi() void operator()(std::uint32_t at) const
  {
    const auto position = listed[at];
    keys[at] = (std::uint64_t{choices[position].row} << 32U) | std::uint64_t{orderedBits(choices[position].rise)};
    positions[at] = position;
  }
};

} // namespace tracker_detail

template <DetectorLayout Det>
class DeviceTrackerImpl final : public DeviceTracker<Det>
{
 public:
  using Report = typename DeviceTracker<Det>::SliceReport;
  using Frame = TimeFrame<Det>;
  using Index = std::uint32_t;
  static constexpr std::size_t kLayerCount = LayerRowMap<Det>::kLayerCount;
  static constexpr std::size_t kInitialStaging = std::size_t{1} << 20U;
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using Records = DeviceStageOutput<SelectionRecordOf<kMaxClusters>::template RecordT>;
  using CapturedRoads = DeviceColumns<tracker_detail::CapturedRoadOf<kMaxClusters>::template RecordT>;
  using IterationReport = typename DeviceTracker<Det>::IterationReport;
  using Tracks = typename DeviceTracker<Det>::Tracks;
  using TrackRecord = TrackOf<kMaxClusters>;
  using Proposal = FollowProposal<kMaxClusters>;
  static_assert(std::is_trivially_copyable_v<TrackingGeometry<Det>>, "the geometry uploads as bytes");

  /// The arena on \a block.
  DeviceTrackerImpl(DeviceRuntime& runtime, DeviceBlock block)
    : mRuntime(&runtime),
      mArena(runtime, block),
      mStream(runtime.createStream()),
      mTrackletStage(runtime, mArena, mStream),
      mNeighbourStage(runtime, mArena, mStream),
      mRecordStage(runtime, mArena, mStream),
      mPartnerStage(runtime, mArena, mStream),
      mRoadBatch(static_cast<std::size_t>(runtime.properties().residentThreads))
  {
    setStagingCapacity(kInitialStaging);
  }

  ~DeviceTrackerImpl() final
  {
    // Work still queued reads the arena: finish it before the arena goes.
    mStream.synchronize();
    closeFrom(Scope::kRun);
  }

  void attach(DeviceBlock block) final
  {
    detach();
    mArena.rebind(block);
  }

  void detach() noexcept final
  {
    // Work still queued reads the block: finish it before the block goes.
    mStream.synchronize();
    closeFrom(Scope::kRun);
    mArena.rebind({});
  }
  DeviceTrackerImpl(const DeviceTrackerImpl&) = delete;
  DeviceTrackerImpl& operator=(const DeviceTrackerImpl&) = delete;
  DeviceTrackerImpl(DeviceTrackerImpl&&) = delete;
  DeviceTrackerImpl& operator=(DeviceTrackerImpl&&) = delete;

  [[nodiscard]] StageStatus configureRun(std::span<const IterationParameters<Det>> iterations, const TrackingGeometry<Det>& geometry, const Propagator& propagator) final
  {
    mStream.synchronize();
    closeFrom(Scope::kRun);
    if (!mParameters.allocate(mArena, iterations.size()) || !mGeometry.allocate(mArena)) {
      closeFrom(Scope::kRun);
      return StageStatus::kOutOfMemory;
    }
    uploadObject(geometry, mGeometry, mStream);
    // The iterations' material is the tracking geometry's (`material()`): on
    // the device, the same place in the device's copy.
    const auto offset = reinterpret_cast<const char*>(&geometry.material()) - reinterpret_cast<const char*>(&geometry);
    const auto* deviceMaterial = reinterpret_cast<const LinkMaterial<Det>*>(reinterpret_cast<const char*>(mGeometry.data()) + offset);
    std::vector<IterationParameters<Det>> deviceParameters(iterations.begin(), iterations.end());
    for (auto& parameters : deviceParameters) {
      if (parameters.material != &geometry.material()) {
        throw std::logic_error("configureRun: an iteration's material is not the tracking geometry's");
      }
      parameters.material = deviceMaterial;
    }
    upload(deviceParameters.data(), mParameters, deviceParameters.size(), mStream);
    // The sensors by place, once, where an iteration follows or the layout
    // has overlap partners.
    mSensorCells = {};
    if (overlapLayers<Det> != 0U || std::ranges::any_of(iterations, [](const auto& parameters) { return parameters.follow; })) {
      // Built on the host once per geometry: configureRun runs every timeframe.
      if (!mHostCells || std::memcmp(mHostCellsGeometry.data(), static_cast<const void*>(&geometry), sizeof(geometry)) != 0) {
        mHostCells.emplace(geometry);
        std::memcpy(mHostCellsGeometry.data(), static_cast<const void*>(&geometry), sizeof(geometry));
      }
      const SensorCells<Det>& cells = *mHostCells;
      if (!mCellLayers.allocate(mArena, cells.layers().size()) || !mCellOffsets.allocate(mArena, cells.offsets().size()) || !mCellSensors.allocate(mArena, std::max<std::size_t>(cells.sensors().size(), 1U)) || !mCellFootprints.allocate(mArena, std::max<std::size_t>(cells.footprints().size(), 1U))) {
        closeFrom(Scope::kRun);
        return StageStatus::kOutOfMemory;
      }
      upload(cells.layers().data(), mCellLayers, cells.layers().size(), mStream);
      upload(cells.offsets().data(), mCellOffsets, cells.offsets().size(), mStream);
      upload(cells.sensors().data(), mCellSensors, cells.sensors().size(), mStream);
      upload(cells.footprints().data(), mCellFootprints, cells.footprints().size(), mStream);
      mSensorCells = {.layers = mCellLayers.data(), .offsets = mCellOffsets.data(), .sensors = mCellSensors.data()};
    }
    mStream.synchronize(); // the host objects may go once this returns
    mHostParameters.assign(iterations.begin(), iterations.end());
    mPropagator = propagator;
    mIterationCount = iterations.size();
    mConfigured = true;
    return StageStatus::kComplete;
  }

  [[nodiscard]] StageStatus beginTimeframe(const Frame& frame, const TimeframePlan<Det>& plan) final
  {
    if (!mConfigured) {
      throw std::logic_error("beginTimeframe before configureRun");
    }
    plan.check(frame, mHostParameters); // before anything changes
    mStream.synchronize();
    closeFrom(Scope::kTimeframe);
    mTimeframeScope.emplace(mArena.scope());
    const std::size_t rows = frame.size();
    const auto offsets = frame.grid().offsets();
    if (!mColumns.allocate(mArena, rows) || !mOffsets.allocate(mArena, offsets.size()) || !mUsed.allocate(mArena, rows) || !mFrame.allocate(mArena) || !uploadPlan(plan)) {
      closeFrom(Scope::kTimeframe);
      return StageStatus::kOutOfMemory;
    }
    upload<MeasurementT>(frame.measurements(), rows, mColumns, mStream);
    upload(offsets.data(), mOffsets, offsets.size(), mStream);
    mRuntime->fill(mUsed.data(), std::uint8_t{kFree}, rows, mStream);
    // The frame is a view: built on the host over the device columns, and
    // copied as the value it is.
    const Frame device{std::as_const(mColumns).view(), frame.grid().geometry(), std::span<const Index>{mOffsets.data(), offsets.size()}, static_cast<Index>(rows)};
    uploadObject(device, mFrame, mStream);
    mStream.synchronize();
    // The tracks' chunks, at the bottom of the scratch end for the timeframe.
    mTracksScope.emplace(mArena.scope(ArenaEnd::kScratch));
    mRows = rows;
    mHostFrame = &frame;
    mPlan = &plan;
    return StageStatus::kComplete;
  }

  void setFlatChunkPairs(std::size_t pairs) final { mFlatOptions.maxChunkPairs = pairs; }
  void setRoadBatch(std::size_t roads) final { mRoadBatch = roads; }

  void setStagingCapacity(std::size_t records) final
  {
    mTrackletStage.prepareStaging(records);
    mNeighbourStage.prepareStaging(records);
    mRecordStage.prepareStaging(records);
  }

  void beginIteration(std::size_t iteration) final
  {
    requireTimeframe();
    if (mIterationOpen) {
      throw std::logic_error("beginIteration with an iteration open");
    }
    if (iteration >= mIterationCount) {
      throw std::invalid_argument("beginIteration of an iteration the run does not have");
    }
    // Every call before returns with the stream drained.
    // The iteration's own arena peak, in a window endIteration closes.
    mOuterArenaPeak = mArena.peak();
    mArena.resetPeak();
    mChunksScope.emplace(mArena.scope(ArenaEnd::kScratch));
    mIterationReport = {};
    mIteration = iteration;
    mIterationStatus = StageStatus::kComplete;
    mIterationOpen = true;
  }

  [[nodiscard]] StageStatus runSlice(std::size_t slice, DeviceSliceDownload<Det>* download) final
  {
    if (!mIterationOpen) {
      throw std::logic_error("runSlice outside an open iteration");
    }
    if (mIterationStatus != StageStatus::kComplete) {
      throw std::logic_error("runSlice in an iteration a slice has failed");
    }
    if (slice >= mPlan->slices.size()) {
      throw std::logic_error("runSlice of a slice the plan does not have");
    }
    openSlice();
    StageStatus status = runStages(mIteration, slice, download);
    if (status == StageStatus::kComplete && mRoadChunks.rows() != 0U && mRoadChunks.rows() >= mRoadBatch) {
      status = fitBatch();
    }
    mStream.synchronize(); // the host objects the slice reads are read until here
    if (status != StageStatus::kComplete) {
      mIterationStatus = status;
      closeFrom(Scope::kSlice);
    }
    return status;
  }

  [[nodiscard]] StageStatus endIteration(DeviceIterationDownload<Det>* download) final
  {
    if (!mIterationOpen) {
      throw std::logic_error("endIteration without an open iteration");
    }
    closeFrom(Scope::kSlice); // runSlice returns with the stream drained
    StageStatus status = mIterationStatus;
    if (status == StageStatus::kComplete && mRoadChunks.rows() != 0U) {
      // The last batch's fit, its outputs in a base scope of their own until
      // its records are appended.
      auto batch = mArena.scope();
      status = fitBatch();
    }
    if (status == StageStatus::kComplete) {
      // The records of every batch, in slice order, into one array on the
      // base end: the selection's input.
      mIterationScope.emplace(mArena.scope());
      mGathered = {};
      const std::size_t recordRows = mChunks.rows();
      mIterationReport.records = recordRows;
      if (status == StageStatus::kComplete && !mGathered.allocate(mArena, recordRows)) {
        status = StageStatus::kOutOfMemory;
      }
      if (status == StageStatus::kComplete) {
        mChunks.gather(mGathered, 0U, mStream);
        if (download != nullptr) {
          gtrk::download(mGathered, download->records, recordRows, mStream);
        }
      }
      // The gather is queued before anything that reuses the chunks' space.
      closeFrom(Scope::kChunks);
      if (status == StageStatus::kComplete) {
        status = endStages(download);
      }
      mStream.synchronize(); // the downloads land, and nothing queued reads the iteration scope
    }
    closeFrom(Scope::kChunks);
    closeFrom(Scope::kIteration);
    mIterationReport.arenaPeakBytes = mArena.peak();
    mArena.raisePeak(mOuterArenaPeak);
    return status;
  }

  void downloadTracks(Tracks& tracks) final
  {
    requireNoIteration("downloadTracks");
    // The iterations' chunks gathered into one block first, while the host
    // sizes the tracks: one copy a column instead of one a column and chunk,
    // each of them a wait for the device.
    if (mTrackChunks.chunks().size() > 1U) {
      auto scratch = mArena.scope(ArenaEnd::kScratch);
      DeviceColumns<TrackRecord::template RecordT> gathered;
      if (!gathered.allocate(mArena, mTrackChunks.rows(), ArenaEnd::kScratch)) {
        throw std::bad_alloc();
      }
      mTrackChunks.gather(gathered, 0U, mStream);
      tracks.resize(mTrackChunks.rows());
      gtrk::download(gathered, tracks, 0U, gathered.size(), mStream);
      mStream.synchronize();
      return;
    }
    tracks.resize(mTrackChunks.rows());
    mTrackChunks.download(tracks, 0U, mStream);
    mStream.synchronize();
  }

  void downloadPartners(std::vector<std::uint32_t>& offsets, std::vector<OverlapPartner>& partners) final
  {
    requireNoIteration("downloadPartners");
    offsets.clear();
    partners.clear();
    if (mFinishedListed.empty()) {
      return;
    }
    // Per track, the listed ones' partners from `finishPartners`.
    offsets.assign(mTrackChunks.rows() + 1U, 0U);
    for (std::size_t e = 0; e < mFinishedListed.size(); ++e) {
      offsets[mFinishedListed[e] + 1U] = mFinishedOffsets[e + 1U] - mFinishedOffsets[e];
    }
    for (std::size_t track = 0; track < mTrackChunks.rows(); ++track) {
      offsets[track + 1U] += offsets[track];
    }
    partners = mFinishedPartners;
  }

  [[nodiscard]] StageStatus finishPartners() final
  {
    requireNoIteration("finishPartners");
    mFinishedListed.clear();
    mFinishedOffsets.clear();
    mFinishedPartners.clear();
    if constexpr (overlapLayers<Det> == 0U) {
      return StageStatus::kComplete;
    } else {
      const auto total = static_cast<Index>(mClaimChunks.rows());
      if (total == 0U) {
        return StageStatus::kComplete;
      }
      auto scratch = mArena.scope(ArenaEnd::kScratch);
      const auto trackCount = static_cast<Index>(mTrackChunks.rows());
      // The claims ordered by track, stably: each track's in the order claimed.
      DeviceColumns<PartnerClaimColumns::template RecordT> claims;
      DeviceSpan<Index> order;
      DeviceSpan<Index> sorted;
      DeviceSpan<Index> sortedTracks;
      DeviceSpan<Index> rows;
      DeviceSpan<Scalar> paths;
      DeviceSpan<std::uint8_t> first;
      DeviceSpan<std::uint8_t> layerRowOf;
      DeviceColumns<TrackRecord::template RecordT> all;
      if (!claims.allocate(mArena, total, ArenaEnd::kScratch) || !order.allocate(mArena, total, ArenaEnd::kScratch) || !sorted.allocate(mArena, total, ArenaEnd::kScratch) || !sortedTracks.allocate(mArena, total, ArenaEnd::kScratch) || !rows.allocate(mArena, total, ArenaEnd::kScratch) || !paths.allocate(mArena, total, ArenaEnd::kScratch) || !first.allocate(mArena, total, ArenaEnd::kScratch) || !layerRowOf.allocate(mArena, total, ArenaEnd::kScratch) || !all.allocate(mArena, trackCount, ArenaEnd::kScratch)) {
        return StageStatus::kOutOfMemory;
      }
      mClaimChunks.gather(claims, 0U, mStream);
      mTrackChunks.gather(all, 0U, mStream);
      const auto claimed = std::as_const(claims).view();
      algorithms::forEach(*mRuntime, mStream, total, tracker_detail::Iota{order.data()});
      auto status = algorithms::sortPairs(*mRuntime, mArena, mStream, claimed.track, sortedTracks.data(), std::as_const(order).data(), sorted.data(), total, algorithms::keyBits(trackCount));
      if (status != StageStatus::kComplete) {
        return status;
      }
      algorithms::forEach(*mRuntime, mStream, total, algorithms::Gather<Index>{claimed.row, std::as_const(sorted).data(), rows.data()});
      algorithms::forEach(*mRuntime, mStream, total, algorithms::Gather<Scalar>{claimed.path, std::as_const(sorted).data(), paths.data()});
      // The tracks with partners and their ranges.
      algorithms::forEach(*mRuntime, mStream, total, tracker_detail::FirstOfValue{std::as_const(sortedTracks).data(), first.data()});
      DeviceSpan<Index> starts;
      Index listedCount = 0U;
      status = algorithms::selectIndices(*mRuntime, mArena, mStream, total, tracker_detail::Flagged{std::as_const(first).data()}, starts, listedCount, ArenaEnd::kScratch);
      if (status != StageStatus::kComplete) {
        return status;
      }
      DeviceSpan<Index> listed;
      DeviceSpan<Index> offsets;
      DeviceSpan<std::uint32_t> merged;
      DeviceSpan<std::uint8_t> mergedRow;
      DeviceSpan<std::uint8_t> kept;
      const std::size_t mergedSize = (std::size_t{listedCount} * kMaxClusters) + total;
      if (!listed.allocate(mArena, listedCount, ArenaEnd::kScratch) || !offsets.allocate(mArena, std::size_t{listedCount} + 1U, ArenaEnd::kScratch) || !merged.allocate(mArena, mergedSize, ArenaEnd::kScratch) || !mergedRow.allocate(mArena, mergedSize, ArenaEnd::kScratch) || !kept.allocate(mArena, listedCount, ArenaEnd::kScratch)) {
        return StageStatus::kOutOfMemory;
      }
      algorithms::forEach(*mRuntime, mStream, listedCount, tracker_detail::ListedRanges{std::as_const(starts).data(), std::as_const(sortedTracks).data(), listedCount, total, listed.data(), offsets.data()});
      // Each fitted with its partners, taking that fit, and back into the chunks.
      const FitWithPartners<Det, true> fit{.parameters = std::as_const(mParameters).data(), .frame = mFrame.data(), .geometry = mGeometry.data(), .propagator = mPropagator, .tracks = all.view(), .listed = std::as_const(listed).data(), .offsets = std::as_const(offsets).data(), .rows = rows.data(), .paths = paths.data(), .fresh = nullptr, .layerRowOf = layerRowOf.data(), .merged = merged.data(), .mergedRow = mergedRow.data(), .kept = kept.data()};
      algorithms::forEach(*mRuntime, mStream, listedCount, fit);
      mTrackChunks.scatter(all, 0U, mTrackChunks.chunks().size(), mStream);
      algorithms::forEach(*mRuntime, mStream, listedCount, tracker_detail::ReleaseRows{std::as_const(kept).data(), std::as_const(offsets).data(), std::as_const(rows).data(), mUsed.data()});
      // The kept ones' partners, for the store (`downloadPartners`).
      std::vector<Index> hostListed(listedCount);
      std::vector<Index> hostOffsets(std::size_t{listedCount} + 1U);
      std::vector<std::uint8_t> hostKept(listedCount);
      std::vector<Index> hostRows(total);
      std::vector<Scalar> hostPaths(total);
      std::vector<std::uint8_t> hostLayerRowOf(total);
      gtrk::download(listed, hostListed.data(), listedCount, mStream);
      gtrk::download(offsets, hostOffsets.data(), hostOffsets.size(), mStream);
      gtrk::download(kept, hostKept.data(), listedCount, mStream);
      gtrk::download(rows, hostRows.data(), total, mStream);
      gtrk::download(paths, hostPaths.data(), total, mStream);
      gtrk::download(layerRowOf, hostLayerRowOf.data(), total, mStream);
      mStream.synchronize();
      mFinishedOffsets.push_back(0U);
      for (Index e = 0; e < listedCount; ++e) {
        if (hostKept[e] == 0U) {
          continue;
        }
        mFinishedListed.push_back(hostListed[e]);
        for (auto at = hostOffsets[e]; at < hostOffsets[e + 1U]; ++at) {
          mFinishedPartners.push_back({.row = hostRows[at], .layerRow = hostLayerRowOf[at], .ahead = static_cast<std::uint8_t>(hostPaths[at] > 0.F ? 1U : 0U)});
        }
        mFinishedOffsets.push_back(static_cast<Index>(mFinishedPartners.size()));
      }
      return StageStatus::kComplete;
    }
  }

  [[nodiscard]] StageStatus findVertices(std::size_t pass, const VertexerSettings& settings, const VertexConstraint& beamSpot, VertexerResult& result) final
  {
    requireNoIteration("findVertices");
    if (pass >= mIterationCount) {
      throw std::invalid_argument("findVertices of an iteration the run does not have");
    }
    // Slice by slice, each its own frames' vertices; the list ordered as the
    // host orders it once the last slice is in.
    result = {};
    const auto& slices = mPlan->slices;
    StageStatus status = StageStatus::kComplete;
    for (std::size_t slice = 0; slice < slices.size() && status == StageStatus::kComplete; ++slice) {
      openSlice();
      DeviceStageOutput<TrackletT> tracklets;
      DeviceSpan<Index> permutation;
      DeviceStageOutput<TripletT> triplets;
      status = tripletStages(pass, slice, tracklets, permutation, triplets);
      VertexerResult found;
      if (status == StageStatus::kComplete) {
        status = vertex_glue::findVertices<Det>(*mRuntime, mArena, mStream, std::as_const(mColumns).view(), static_cast<Index>(mRows), constView(tracklets), constView(triplets), triplets.size(), mParameters.data() + pass, settings, beamSpot, ownedFrames(mHostFrame->grid().geometry(), slices[slice].slice.coreBlocks), found);
      }
      // A complete vertexer has synchronised after its last readback.
      if (status != StageStatus::kComplete) {
        mStream.synchronize();
      }
      closeFrom(Scope::kSlice);
      appendVertices(result, found);
    }
    if (status == StageStatus::kComplete) {
      orderVertices(result);
    }
    return status;
  }

  [[nodiscard]] StageStatus replacePlan(const TimeframePlan<Det>& plan) final
  {
    requireNoIteration("replacePlan");
    plan.check(*mHostFrame, mHostParameters); // before anything changes
    mStream.synchronize();
    closeFrom(Scope::kSlice);
    // The new device part goes on top of the timeframe's base end, the old
    // one's few kB staying until the timeframe ends.
    if (!uploadPlan(plan)) {
      return StageStatus::kOutOfMemory;
    }
    mPlan = &plan;
    return StageStatus::kComplete;
  }

  void endTimeframe() noexcept final
  {
    // Work still queued reads the timeframe's storage.
    mStream.synchronize();
    closeFrom(Scope::kTimeframe);
  }

  [[nodiscard]] const IterationReport& lastIteration() const final { return mIterationReport; }

  [[nodiscard]] const Report& lastSlice() const final { return mReport; }
  [[nodiscard]] std::size_t iterationCount() const final { return mIterationCount; }
  [[nodiscard]] std::size_t arenaBytes() const final { return mArena.capacity(); }
  [[nodiscard]] std::size_t arenaUsedBytes(ArenaEnd end) const final { return mArena.used(end); }
  [[nodiscard]] std::size_t arenaPeakBytes() const final { return mArena.peak(); }

 private:
  struct LayeredNeighbours {
    NeighbourStage<Det> body;
    const Index* layerBoundaries;
    std::array<LayerID, kLayerCount> layerOrder = ::o2::gtrk::layerDepthGroups<Det>.order;
  };

  /// The plan's device part, in the timeframe scope: every slice's source map and core ranges, once for every
  /// iteration -- the core ranges with an iteration's start-layer mask per
  /// iteration, where one has a mask --, each iteration's vertex windows, the
  /// beam line and then the vertices, the follow's bound.
  /// False, with the scope to be closed by the caller, if the arena refuses.
  [[nodiscard]] bool uploadPlan(const TimeframePlan<Det>& plan)
  {
    std::vector<LayerRowMap<Det>> sources;
    std::vector<MeasurementRange> cores;
    std::vector<MeasurementRange> maskedCores; // per iteration, with its start-layer mask (RoadCount.h, `roadCores`)
    const std::vector<VertexWindowParameters<Det>>& windows = plan.vertexWindows;
    for (const auto& slice : plan.slices) {
      sources.push_back(slice.sources);
      cores.insert(cores.end(), slice.slice.core.begin(), slice.slice.core.end());
    }
    const bool startMasks = std::ranges::any_of(mHostParameters, [](const IterationParameters<Det>& one) { return one.roadStartLayers != 0U; });
    if (startMasks) {
      for (std::size_t iteration = 0U; iteration < plan.parameters.size(); ++iteration) {
        for (const auto& slice : plan.slices) {
          const auto masked = roadCores(std::span<const MeasurementRange, kLayerCount>{slice.slice.core}, mHostParameters[iteration].roadStartLayers);
          maskedCores.insert(maskedCores.end(), masked.begin(), masked.end());
        }
      }
    }
    std::vector<VertexSeed> seeds{plan.beamLine};
    seeds.insert(seeds.end(), plan.vertices.begin(), plan.vertices.end());
    if (!mSliceSources.allocate(mArena, sources.size()) || !mSliceCores.allocate(mArena, cores.size()) || (startMasks && !mSliceRoadCores.allocate(mArena, maskedCores.size())) || !mVertexWindows.allocate(mArena, windows.size()) || !mSeeds.allocate(mArena, seeds.size()) || !mLayerMaxWidthPS.allocate(mArena, kLayerCount) || !mLayerMaxLocVariance.allocate(mArena, kLayerCount) || !mTargetBounds.allocate(mArena, kLayerCount)) {
      return false;
    }
    upload(sources.data(), mSliceSources, sources.size(), mStream);
    upload(cores.data(), mSliceCores, cores.size(), mStream);
    if (startMasks) {
      upload(maskedCores.data(), mSliceRoadCores, maskedCores.size(), mStream);
    }
    mRoadCoresMasked = startMasks;
    upload(windows.data(), mVertexWindows, windows.size(), mStream);
    upload(seeds.data(), mSeeds, seeds.size(), mStream);
    upload(plan.layerMaxWidthPS.data(), mLayerMaxWidthPS, kLayerCount, mStream);
    upload(plan.layerMaxLocVariance.data(), mLayerMaxLocVariance, kLayerCount, mStream);
    upload(plan.layerTargetBounds.data(), mTargetBounds, kLayerCount, mStream);
    mStream.synchronize(); // the host vectors go when this returns
    return true;
  }

  template <template <template <class> class> class S>
  [[nodiscard]] static typename DeviceColumns<S>::ConstView constView(const DeviceStageOutput<S>& output)
  {
    return output.view();
  }

  /// A stage's run over a flat candidate list
  /// into \a output, its report into \a run, its readbacks counted.
  template <template <template <class> class> class S, typename Adapter>
  StageStatus runFlat(const Index* listed, Index listedCount, const Adapter& adapter, DeviceStageOutput<S>& output, FlatRunReport& run)
  {
    const StageStatus status = flat::runFlat(*mRuntime, mArena, mStream, listed, listedCount, adapter, output, run, mFlatOptions);
    mReport.readbacks += run.readbacks;
    return status;
  }

  /// A stage's run, counted: its readback, and a second one after a missed
  /// prediction.
  template <template <template <class> class> class S, typename Generate>
  StageStatus runStage(GpuStage<S>& gpuStage, DeviceStageOutput<S>& output, const Generate& generate)
  {
    return counted(gpuStage, gpuStage.run(output, generate));
  }

  /// The same over the \a listedCount elements \a listed names only, when set
  /// (`GpuStage::run`).
  template <template <template <class> class> class S, typename Generate>
  StageStatus runStage(GpuStage<S>& gpuStage, DeviceStageOutput<S>& output, const Generate& generate, const Index* listed, Index listedCount)
  {
    return counted(gpuStage, gpuStage.run(output, generate, listed, listedCount));
  }

  /// \a status of \a gpuStage's last run, its readbacks counted.
  template <template <template <class> class> class S>
  StageStatus counted(const GpuStage<S>& gpuStage, StageStatus status)
  {
    const bool missed = gpuStage.stagingPredictionMissed();
    mReport.readbacks += missed ? 2U : 1U;
    mReport.stagingPredictionMissed = mReport.stagingPredictionMissed || missed;
    return status;
  }

  /// A slice's tracklets and triplets, the first stages every slice runs
  /// and all the vertexing pass runs, into \a tracklets and \a triplets in
  /// the slice scope, with the triplet stage's first-tracklet order in
  /// \a permutation, with iteration \a iteration's parameters.
  StageStatus tripletStages(std::size_t iteration, std::size_t sliceIndex, DeviceStageOutput<TrackletT>& tracklets, DeviceSpan<Index>& permutation, DeviceStageOutput<TripletT>& triplets)
  {
    const IterationParameters<Det>* parameters = mParameters.data() + iteration;
    const IterationParameters<Det>& hostParameters = mHostParameters[iteration];
    const auto& slicePlan = mPlan->slices[sliceIndex];
    // The slices' sources and cores, uploaded once for every iteration.
    const LayerRowMap<Det>* sources = mSliceSources.data() + sliceIndex;
    const MeasurementRange* core = mSliceCores.data() + (sliceIndex * kLayerCount);
    const Index sourceCount = slicePlan.sources.size();

    // The tracklet stage as the host tracker builds it from the plan, bound to
    // the device copies.
    const auto inputs = slicePlan.inputs(hostParameters, mPlan->layerMaxWidthPS);
    const bool beamLine = hostParameters.vertexSource == VertexSource::kBeamLine;
    if (!tracklets.prepare(mArena, sourceCount)) {
      return StageStatus::kOutOfMemory;
    }
    // The beam-line stage for a beam-line iteration, testing the shapes where
    // the iteration does.
    StageStatus status = withTrackletVertices(hostParameters.vertexSource, [&](auto beamLineOnly) {
      return withShapeTest(hostParameters.shapeTest, [&](auto shape) {
        return withDiscWindow<decltype(beamLineOnly)::value>(hostParameters.discWindowMinP > 0.F, [&](auto discWindow) {
          const TrackletStageFor<Det, decltype(beamLineOnly)::value, decltype(shape)::value, decltype(discWindow)::value> stage{hostParameters, slicePlan.sources, *mHostFrame, inputs.targetBlocks, mPlan->beamLine, mPlan->vertices, inputs.maxCandidateWidthPS, mPlan->vertexWindows[iteration], nullptr, inputs.gate};
          const auto bound = stage.bind({
            .parameters = parameters,
            .sources = sources,
            .frame = mFrame.data(),
            .used = mUsed.data(),
            .vertexWindows = mVertexWindows.data() + iteration,
            .vertices = beamLine ? mSeeds.data() : mSeeds.data() + 1,
            .gateCore = core,
          });
          // A thread for the sources with a query only: a dead source's thread holds its warp's lanes idle while the
          // others walk; it emits nothing, so its count is the 0 the stage sets.
          auto scratch = mArena.scope(ArenaEnd::kScratch);
          DeviceSpan<Index> listed;
          Index listedCount = 0U;
          const StageStatus selected = algorithms::selectIndices(*mRuntime, mArena, mStream, sourceCount, tracker_detail::HasQuery<decltype(bound)>{bound}, listed, listedCount, ArenaEnd::kScratch);
          ++mReport.readbacks; // the list's length
          if (selected != StageStatus::kComplete) {
            return selected;
          }
          // The flat run: the sources' candidates listed, tested one a thread.
          std::uint32_t queriesPerSource = 1U;
          for (std::size_t layer = 0; layer < kLayerCount; ++layer) {
            queriesPerSource = std::max(queriesPerSource, hostParameters.outLinkOffsets[layer + 1U] - hostParameters.outLinkOffsets[layer]);
          }
          if (queriesPerSource > (1U << tracker_detail::FlatTracklets<decltype(bound)>::kQueryBits)) {
            throw std::logic_error("a layer has more out-links than a flat pair's key numbers");
          }
          mReport.trackletSources = listedCount;
          return runFlat(listedCount > 0U ? listed.data() : nullptr, listedCount, tracker_detail::FlatTracklets<decltype(bound)>{.stage = bound, .queriesPerSource = queriesPerSource}, tracklets, mReport.trackletRun);
        });
      });
    });
    if (status != StageStatus::kComplete) {
      return status;
    }
    // Each source's tracklets in (link, lambda, target) order, which the
    // triplet stage searches.
    status = glue::sortTracklets(*mRuntime, mArena, mStream, tracklets);
    if (status != StageStatus::kComplete) {
      return status;
    }
    mReport.tracklets = tracklets.size();

    // The longest tracklet paths through each measurement, which the
    // permutation and the triplet stage read.
    glue::DevicePathBounds trackletPaths;
    status = glue::trackletPaths<Det>(*mRuntime, mArena, mStream, parameters, slicePlan.sources, constView(tracklets), tracklets.offsets().data(), trackletPaths);
    if (status != StageStatus::kComplete) {
      return status;
    }
    const auto paths = trackletPaths.view(hostParameters.minClusters);

    // The triplet stage's inputs: the first-tracklet order; the targets'
    // variance bounds are the plan's, over the whole frame, uploaded with it.
    Index firstTracklets = 0U;
    status = glue::permutation<Det>(*mRuntime, mArena, mStream, constView(tracklets), tracklets.size(), parameters, sources, permutation, firstTracklets, tracklets.offsets().data(), paths);
    ++mReport.readbacks;
    if (status != StageStatus::kComplete) {
      return status;
    }
    mReport.permutation = firstTracklets;

    if (!triplets.prepare(mArena, firstTracklets)) {
      return StageStatus::kOutOfMemory;
    }
    // The flat run: the windows'
    // candidates listed per first tracklet, tested one a thread; the
    // shape-testing stage where the iteration tests the shapes.
    status = withShapeTest(hostParameters.shapeTest, [&](auto shape) {
      const TripletStage<Det, decltype(shape)::value> tripletStage{parameters, sources, mFrame.data(), constView(tracklets), tracklets.offsets().data(), mTargetBounds.data(), permutation.data(), paths};
      return runFlat(nullptr, firstTracklets, tracker_detail::FlatTriplets<Det, decltype(shape)::value>{.stage = tripletStage}, triplets, mReport.tripletRun);
    });
    if (status != StageStatus::kComplete) {
      return status;
    }
    mReport.triplets = triplets.size();
    return StageStatus::kComplete;
  }

  StageStatus runStages(std::size_t iteration, std::size_t sliceIndex, DeviceSliceDownload<Det>* download)
  {
    const IterationParameters<Det>* parameters = mParameters.data() + iteration;
    const IterationParameters<Det>& hostParameters = mHostParameters[iteration];
    // The unmasked cores by slice; the masked ones by iteration and slice (`uploadPlan`).
    const MeasurementRange* core = mRoadCoresMasked ? mSliceRoadCores.data() + (((iteration * mPlan->slices.size()) + sliceIndex) * kLayerCount) : mSliceCores.data() + (sliceIndex * kLayerCount);
    DeviceStageOutput<TrackletT> tracklets;
    DeviceSpan<Index> permutation;
    DeviceStageOutput<TripletT> triplets;
    StageStatus status = tripletStages(iteration, sliceIndex, tracklets, permutation, triplets);
    if (status != StageStatus::kComplete) {
      return status;
    }

    // The longest triplet paths through each tracklet, from the triplets
    // alone, and the triplets kept only where they are on one:
    // what the neighbour stage reads.
    glue::DevicePathBounds tripletPaths;
    status = glue::tripletPaths<Det>(*mRuntime, mArena, mStream, parameters, mPlan->slices[sliceIndex].sources, constView(tracklets), tracklets.size(), tracklets.offsets().data(), constView(triplets), triplets.offsets().data(), permutation.data(), mReport.permutation, tripletPaths);
    if (status != StageStatus::kComplete) {
      return status;
    }
    const auto paths = tripletPaths.view(hostParameters.minClusters);
    {
      DeviceStageOutput<TripletT> kept;
      status = glue::onPathTriplets(*mRuntime, mArena, mStream, paths, triplets, kept, tripletPaths);
      ++mReport.readbacks;
      if (status != StageStatus::kComplete) {
        return status;
      }
      triplets = kept;
    }
    mReport.triplets = triplets.size();

    DeviceSpan<Index> endOffsets;
    DeviceSpan<Index> endTriplets;
    status = glue::endIndex(*mRuntime, mArena, mStream, constView(triplets), triplets.size(), tracklets.size(), endOffsets, endTriplets);
    if (status != StageStatus::kComplete) {
      return status;
    }

    DeviceStageOutput<NeighbourT> neighbours;
    if (!neighbours.prepare(mArena, triplets.size())) {
      return StageStatus::kOutOfMemory;
    }
    road_glue::OrderedCounts counts;
    status = counts.prepare(*mRuntime, mArena, mStream, triplets.size());
    if (status != StageStatus::kComplete) {
      return status;
    }
    const NeighbourStage<Det> neighbourStage{parameters, mFrame.data(), constView(tracklets), constView(triplets), endOffsets.data(), endTriplets.data(), paths, triplets.view().level, counts.view(road_glue::linkSourcesOf(parameters), core)};
    status = runStage(mNeighbourStage, neighbours, LayeredNeighbours{neighbourStage, tripletPaths.layerTriplets.data()});
    if (status != StageStatus::kComplete) {
      return status;
    }
    mReport.neighbours = neighbours.size();

    // The chained count and the overflow flag in one readback; with a
    // download, also the chained triplets' list, which the tests compare
    // with the host's.
    Index chainedCount = 0U;
    status = counts.summarize(*mRuntime, mStream, constView(triplets), triplets.size(), chainedCount);
    mReport.readbacks += triplets.size() != 0U ? 1U : 0U;
    DeviceSpan<Index> chained;
    Index chainedListed = 0U;
    if (status == StageStatus::kComplete && download != nullptr) {
      status = glue::chained(*mRuntime, mArena, mStream, neighbours.offsets().data(), triplets.size(), chained, chainedListed);
      ++mReport.readbacks;
    }
    mReport.chained = chainedCount;
    if (status != StageStatus::kComplete) {
      return status;
    }

    // Counts and successor flags were finalized with the neighbours; overflow
    // was checked together with the chained summary above.
    const auto& successorAbove = counts.successorAbove;
    const auto& roads = counts.roads;
    DeviceSpan<Index> starts;
    Index startCount = 0U;
    status = road_glue::roadStarts(*mRuntime, mArena, mStream, constView(tracklets), constView(triplets), triplets.size(), successorAbove.data(), roads.data(), road_glue::roadStartRuleOf(parameters, hostParameters), starts, startCount);
    ++mReport.readbacks;
    if (status != StageStatus::kComplete) {
      return status;
    }
    DeviceSpan<Index> roadOffsets;
    Index roadCount = 0U;
    status = road_glue::roadOffsets(*mRuntime, mArena, mStream, starts.data(), startCount, roads.data(), roadOffsets, roadCount);
    ++mReport.readbacks;
    if (status != StageStatus::kComplete) {
      return status;
    }
    mReport.roads = roadCount;

    // The roads captured for their batch's fit (fitBatch).
    if (roadCount != 0U) {
      if (!mBatchScope) {
        mBatchScope.emplace(mArena.scope(ArenaEnd::kScratch));
      }
      CapturedRoads captured;
      if (!captured.allocate(mArena, roadCount, ArenaEnd::kScratch)) {
        return StageStatus::kOutOfMemory;
      }
      const RoadRows<Det> rows{constView(tracklets), constView(triplets), RoadWalk{starts.data(), roadOffsets.data(), startCount, neighbours.offsets().data(), neighbours.view().from, triplets.view().level, roads.data()}};
      algorithms::forEach(*mRuntime, mStream, roadCount, tracker_detail::CaptureRoad<Det>{.rows = rows, .roads = captured.view()});
      mRoadChunks.adopt(captured);
    }

    if (download != nullptr) {
      gtrk::download(successorAbove, download->successorAbove, triplets.size(), mStream);
      gtrk::download(roads, download->roads, triplets.size(), mStream);
      gtrk::download(starts, download->starts, startCount, mStream);
      gtrk::download(roadOffsets, download->roadOffsets, std::size_t{startCount} + 1U, mStream);
      gtrk::download(tracklets, download->tracklets, mStream);
      gtrk::download(triplets, download->triplets, mStream);
      gtrk::download(neighbours, download->neighbours, mStream);
      gtrk::download(chained, download->chained, chainedListed, mStream);
    }
    return StageStatus::kComplete; // runSlice synchronises
  }

  /// The batch's captured roads fitted, in slice order, by one run, the
  /// records appended to the iteration's chunks, and the batch's scope
  /// closed. The fit's outputs in the open base scope; its readbacks the
  /// iteration's.
  [[nodiscard]] StageStatus fitBatch()
  {
    StageStatus status = StageStatus::kComplete;
    Records records;
    {
      // The batch's roads in one array: the one slice's where it is alone.
      CapturedRoads roads;
      const auto chunks = mRoadChunks.chunks();
      if (chunks.size() == 1U) {
        roads = chunks.front();
      } else if (!roads.allocate(mArena, mRoadChunks.rows(), ArenaEnd::kScratch)) {
        status = StageStatus::kOutOfMemory;
      } else {
        mRoadChunks.gather(roads, 0U, mStream);
      }
      // The order the roads are fitted in (`RoadOrderKey`); the records are
      // placed by road whatever the order.
      const std::size_t count = roads.size();
      DeviceSpan<std::uint64_t> keys;
      DeviceSpan<std::uint64_t> sortedKeys;
      DeviceSpan<std::uint32_t> indices;
      DeviceSpan<std::uint32_t> order;
      if (status == StageStatus::kComplete && count > std::size_t{std::numeric_limits<Index>::max()}) {
        status = StageStatus::kIndexOverflow;
      }
      if (status == StageStatus::kComplete && (!keys.allocate(mArena, count, ArenaEnd::kScratch) || !sortedKeys.allocate(mArena, count, ArenaEnd::kScratch) || !indices.allocate(mArena, count, ArenaEnd::kScratch) || !order.allocate(mArena, count, ArenaEnd::kScratch))) {
        status = StageStatus::kOutOfMemory;
      }
      if (status == StageStatus::kComplete) {
        algorithms::forEach(*mRuntime, mStream, static_cast<std::uint32_t>(count), tracker_detail::RoadOrderKey<Det, kMaxClusters>{mGeometry.data(), mFrame.data(), std::as_const(roads).view(), keys.data(), indices.data()});
        status = algorithms::sortPairs(*mRuntime, mArena, mStream, std::as_const(keys).data(), sortedKeys.data(), std::as_const(indices).data(), order.data(), static_cast<std::uint32_t>(count), 64);
      }
      if (status == StageStatus::kComplete) {
        status = fitRoads(roads, order.data(), records);
      }
    }
    // The fits are queued before anything that reuses the roads' space.
    mRoadChunks.clear();
    mBatchScope.reset();
    if (status != StageStatus::kComplete) {
      return status;
    }
    if (!mChunks.append(mArena, ArenaEnd::kScratch, records.output(), 0U, records.size(), mStream)) {
      return StageStatus::kOutOfMemory;
    }
    return StageStatus::kComplete;
  }

  /// The fit of the captured \a roads with the iteration's settings into
  /// \a records, in the open base scope: one stage run, its threads taking
  /// the roads in \a order.
  [[nodiscard]] StageStatus fitRoads(const CapturedRoads& roads, const Index* order, Records& records)
  {
    if (!records.prepare(mArena, static_cast<Index>(roads.size()))) {
      return StageStatus::kOutOfMemory;
    }
    const IterationParameters<Det>* parameters = mParameters.data() + mIteration;
    const IterationParameters<Det>& hostParameters = mHostParameters[mIteration];
    // An iteration with a kink allowance takes it on the fit's last step
    // (against the road-end fakes): its own instantiation. An iteration with
    // the beam constraint fits on the stage that also tests the beam line, as
    // the host.
    const auto run = [&]<bool BeamLine, bool LastStepAllowance>() {
      using Fit = RoadRowsFit<Det, BeamLine, LastStepAllowance>;
      const tracker_detail::CapturedRoadFit<Fit, kMaxClusters> fit{.fit = Fit{parameters, mFrame.data(), mGeometry.data(), mPropagator}, .roads = std::as_const(roads).view()};
      return mRecordStage.run(records, fit, order, static_cast<Index>(roads.size()));
    };
    const StageStatus status = withBools([&]<bool BeamLine, bool LastStepAllowance>(std::bool_constant<BeamLine>, std::bool_constant<LastStepAllowance>) { return run.template operator()<BeamLine, LastStepAllowance>(); }, hostParameters.beamConstraint, takesLastStepAllowance(hostParameters));
    mIterationReport.readbacks += mRecordStage.stagingPredictionMissed() ? 2U : 1U;
    return status;
  }

  /// The end of an iteration after its gather, in the iteration scope: the
  /// selection, the refit and the
  /// follow on a working copy of the used mask, the tracks appended to
  /// the timeframe's, and only then the working mask committed. A failure
  /// before the commit leaves the used mask and the tracks as they were.
  [[nodiscard]] StageStatus endStages(DeviceIterationDownload<Det>* download)
  {
    const Index rows = static_cast<Index>(mRows);
    const auto records = std::as_const(mGathered).view();
    const auto recordCount = static_cast<Index>(mGathered.size());
    const IterationParameters<Det>& hostParameters = mHostParameters[mIteration];
    const IterationParameters<Det>* parameters = mParameters.data() + mIteration;
    auto& report = mIterationReport;

    DeviceSpan<std::uint8_t> working;
    if (!working.allocate(mArena, rows)) {
      return StageStatus::kOutOfMemory;
    }
    mRuntime->copyOnDevice(working.data(), std::as_const(mUsed).data(), rows, mStream);

    selection_glue::Selected selected;
    StageStatus status = selection_glue::select<kMaxClusters>(*mRuntime, mArena, mStream, records, recordCount, working.data(), rows, hostParameters.sharedInnermost, selected);
    report.readbacks += selected.readbacks;
    if (status != StageStatus::kComplete) {
      return status;
    }
    report.duplicates = selected.duplicates;
    report.mismatches = selected.mismatches;
    report.rounds = selected.rounds;
    report.tracks = selected.tracks;

    // The refit, one element per selected record, into the iteration's tracks.
    const Index trackCount = selected.tracks;
    DeviceColumns<TrackRecord::template RecordT> tracks;
    DeviceSpan<Proposal> proposals;
    DeviceSpan<ClusterList<kMaxClusters>> followRows;
    DeviceSpan<std::uint8_t> followCounts;
    const Index followCount = hostParameters.follow ? trackCount : 0U;
    if (!tracks.allocate(mArena, trackCount) || !proposals.allocate(mArena, followCount) || !followRows.allocate(mArena, followCount) || !followCounts.allocate(mArena, followCount)) {
      return StageStatus::kOutOfMemory;
    }
    const TrackRefitElement<Det> refit{mFrame.data(), mGeometry.data(), mPropagator, records.clusters, records.count, records.outwardState, records.chi2, std::as_const(selected.selected).data(), tracks.view(), 0U, static_cast<std::uint8_t>(mIteration), kinkBudget(hostParameters)};
    algorithms::forEach(*mRuntime, mStream, trackCount, refit);

    // The follow, the host's sequence
    // (`follow::runSequence`), its bids in the rounds' slots. The host's
    // priority, trackCount - (trackOffset + track), is the iteration's count
    // less the track.
    if (hostParameters.follow && trackCount != 0U) {
      runFollow(mIteration, tracks, trackCount, false, proposals.data(), followRows, followCounts, selected.slots, rows, working);
    }

    // The drop after the follow (the host: Tracker::dropTracks): the
    // iteration's tracks whose refit failed or below a length or a pT floor
    // after the follow dropped, the kept ones and their proposals compacted, and the mask
    // rebuilt as the one before the selection plus the kept tracks' rows as
    // the follow left them. Any refit may fail, so it runs in every iteration.
    Index keptCount = trackCount;
    if (trackCount != 0U) {
      auto scratch = mArena.scope(ArenaEnd::kScratch);
      DeviceSpan<Index> kept;
      status = algorithms::selectIndices(*mRuntime, mArena, mStream, trackCount, tracker_detail::KeptAfterFollow<Det>{parameters, std::as_const(tracks).view().count, std::as_const(tracks).view().state, std::as_const(tracks).view().status}, kept, keptCount, ArenaEnd::kScratch);
      ++report.readbacks;
      if (status != StageStatus::kComplete) {
        return status;
      }
      if (keptCount != trackCount) {
        DeviceColumns<TrackRecord::template RecordT> keptTracks;
        DeviceSpan<Proposal> keptProposals;
        DeviceSpan<std::uint8_t> cut;
        if (!keptTracks.allocate(mArena, keptCount) || !keptProposals.allocate(mArena, followCount != 0U ? keptCount : 0U) || !cut.allocate(mArena, rows)) {
          return StageStatus::kOutOfMemory;
        }
        device_storage_detail::forEachColumnPair(std::as_const(tracks).view(), keptTracks.view(), [&](auto* to, const auto* from) { algorithms::forEach(*mRuntime, mStream, keptCount, algorithms::Gather<std::remove_pointer_t<decltype(to)>>{from, kept.data(), to}); });
        if (followCount != 0U) {
          algorithms::forEach(*mRuntime, mStream, keptCount, algorithms::Gather<Proposal>{std::as_const(proposals).data(), kept.data(), keptProposals.data()});
        }
        mRuntime->copyOnDevice(cut.data(), std::as_const(mUsed).data(), rows, mStream);
        algorithms::forEach(*mRuntime, mStream, keptCount, tracker_detail::RestoreRows<kMaxClusters>{std::as_const(keptTracks).view().clusters, std::as_const(keptTracks).view().count, std::as_const(working).data(), cut.data()});
        tracks = keptTracks;
        proposals = keptProposals;
        working = cut;
      }
    }

    // The join to the earlier iterations' tracks (the host: Tracker::joinTracks)
    // where they are, on every track so far gathered into one block; the
    // earlier tracks it changes go back into their chunks as the last step,
    // after every fallible one.
    DeviceColumns<TrackRecord::template RecordT> gathered;
    const auto earlier = static_cast<Index>(mTrackChunks.rows());
    const std::size_t earlierChunks = mTrackChunks.chunks().size();
    const bool joins = earlier != 0U && keptCount != 0U;
    DeviceSpan<Index> changed;
    Index joinedCount = 0U;
    if (joins) {
      if (!gathered.allocate(mArena, std::size_t{earlier} + keptCount)) {
        return StageStatus::kOutOfMemory;
      }
      mTrackChunks.gather(gathered, 0U, mStream);
      copyColumns(tracks, 0U, gathered, earlier, keptCount, mStream);
      status = joinTracks(gathered, earlier, tracks, proposals, keptCount, followCount != 0U, changed, joinedCount);
      if (status != StageStatus::kComplete) {
        return status;
      }
      report.tracksJoined = joinedCount;
    }

    // The inactive layers' rows (the host: InactiveRows::run) of the earlier
    // tracks the join changed, then the iteration's own.
    if constexpr (Detector<Det>::inactiveLayers != 0U) {
      if (keptCount + joinedCount != 0U) {
        Index taken = 0U;
        status = joins ? inactiveRows(parameters, gathered.view(), earlier, keptCount, changed, joinedCount, working, taken) : inactiveRows(parameters, tracks.view(), 0U, keptCount, changed, 0U, working, taken);
        if (status != StageStatus::kComplete) {
          return status;
        }
        report.inactiveRows = taken;
        if (joins && keptCount != 0U) {
          copyColumns(gathered, earlier, tracks, 0U, keptCount, mStream);
        }
      }
    }

    // The overlap partners (the host: Tracker::claimPartners), claimed into
    // the used mask, for the iteration's tracks and the earlier ones it or the
    // follow of earlier tracks before it changed.
    DeviceColumns<PartnerClaimColumns::template RecordT> claims;
    Index claimCount = 0U;
    if constexpr (overlapLayers<Det> != 0U) {
      if (keptCount + joinedCount + mCarriedChunks.rows() != 0U) {
        status = claimPartners(tracks, keptCount, joins ? &gathered : nullptr, earlier, changed, joinedCount, working, claims, claimCount);
        if (status != StageStatus::kComplete) {
          return status;
        }
        report.overlapPartners = claimCount;
      }
    }

    // The successor's recovery pass (the host: Tracker::runIteration's
    // followEarlier): every track so far followed
    // again with the successor's settings before it makes its seeds, the rows
    // it takes its own. Here, with this iteration's tracks final and before
    // the commit, nothing runs in between: the tracks so far gathered with
    // this iteration's, followed, and put back as one chunk; the old chunks'
    // storage stays until the timeframe ends.
    bool committed = false;
    if (mIteration + 1U < mIterationCount && mHostParameters[mIteration + 1U].followEarlier && mHostParameters[mIteration + 1U].follow && mTrackChunks.rows() + keptCount != 0U) {
      status = followEarlier(mIteration + 1U, tracks, keptCount, working, selected.slots, rows, joinedCount != 0U ? &gathered : nullptr);
      if (status != StageStatus::kComplete) {
        return status;
      }
      committed = true;
    }

    if (download != nullptr) {
      gtrk::download(selected.ranked, download->ranked, selected.unique, mStream);
      gtrk::download(selected.selected, download->selected, trackCount, mStream);
      gtrk::download(proposals, download->proposals, proposals.size(), mStream);
      gtrk::download(tracks, download->tracks, keptCount, mStream);
      gtrk::download(working, download->used, rows, mStream);
    }

    // The last fallible steps, then the commit.
    if (claimCount != 0U && !mClaimChunks.append(mArena, ArenaEnd::kScratch, claims, 0U, claimCount, mStream)) {
      return StageStatus::kOutOfMemory;
    }
    if (!committed && !mTrackChunks.append(mArena, ArenaEnd::kScratch, tracks, 0U, keptCount, mStream)) {
      return StageStatus::kOutOfMemory;
    }
    if (!committed && joinedCount != 0U) {
      mTrackChunks.scatter(gathered, 0U, earlierChunks, mStream);
    }
    mRuntime->copyOnDevice(mUsed.data(), std::as_const(working).data(), rows, mStream);
    return StageStatus::kComplete; // endIteration synchronises
  }

  /// The iteration's overlap partners on the device (see OverlapPartners.h),
  /// the host's `Tracker::claimPartners` and `OverlapPartners::run`. The
  /// candidates: the \a count tracks of \a tracks (after the \a earlier
  /// ones), the \a joined earlier tracks the join changed (\a changed) and
  /// those the follow of earlier tracks changed (`mCarriedChunks`, used up),
  /// ascending, each once; every track so far in \a gathered where the join
  /// gathered them. Their proposals (`ProposePartners`), a sensor holding one
  /// of a candidate's claims passed over; a row wanted by more than one to
  /// the lowest chi2, then the lowest candidate (a stable sort of the
  /// proposals, in candidate order, by row and chi2); each candidate with
  /// winners fitted on its rows, its claims and the winners
  /// (`FitWithPartners`, a test); the kept ones' winners claimed into
  /// \a claims, \a claimed of them, and taken in \a working.
  [[nodiscard]] StageStatus claimPartners(DeviceColumns<TrackRecord::template RecordT>& tracks, Index count, DeviceColumns<TrackRecord::template RecordT>* gathered, Index earlier, const DeviceSpan<Index>& changed, Index joined, DeviceSpan<std::uint8_t>& working, DeviceColumns<PartnerClaimColumns::template RecordT>& claims, Index& claimed)
  {
    auto& report = mIterationReport;
    claimed = 0U;
    // Every track so far, at its global index.
    DeviceColumns<TrackRecord::template RecordT> all = tracks;
    if (gathered != nullptr) {
      all = *gathered;
    } else if (earlier != 0U) {
      if (!all.allocate(mArena, std::size_t{earlier} + count)) {
        return StageStatus::kOutOfMemory;
      }
      mTrackChunks.gather(all, 0U, mStream);
      copyColumns(tracks, 0U, all, earlier, count, mStream);
    }
    // The candidates, ascending, each once.
    const auto carried = static_cast<Index>(mCarriedChunks.rows());
    const Index listedCount = carried + joined + count;
    DeviceColumns<tracker_detail::TrackIndexColumns::template RecordT> listed;
    DeviceSpan<Index> sortedListed;
    DeviceSpan<std::uint8_t> firstListed;
    if (!listed.allocate(mArena, listedCount) || !sortedListed.allocate(mArena, listedCount) || !firstListed.allocate(mArena, listedCount)) {
      return StageStatus::kOutOfMemory;
    }
    mCarriedChunks.gather(listed, 0U, mStream);
    mCarriedChunks.clear();
    if (joined != 0U) {
      mRuntime->copyOnDevice(listed.view().track + carried, changed.data(), std::size_t{joined} * sizeof(Index), mStream);
    }
    algorithms::forEach(*mRuntime, mStream, count, tracker_detail::OwnTracks{earlier, listed.view().track + carried + joined});
    const int trackBits = algorithms::keyBits(earlier + count);
    auto status = algorithms::sortKeys(*mRuntime, mArena, mStream, std::as_const(listed).view().track, sortedListed.data(), listedCount, trackBits);
    if (status != StageStatus::kComplete) {
      return status;
    }
    algorithms::forEach(*mRuntime, mStream, listedCount, tracker_detail::FirstOfValue{std::as_const(sortedListed).data(), firstListed.data()});
    DeviceSpan<Index> firsts;
    Index n = 0U;
    status = algorithms::selectIndices(*mRuntime, mArena, mStream, listedCount, tracker_detail::Flagged{std::as_const(firstListed).data()}, firsts, n);
    ++report.readbacks;
    if (status != StageStatus::kComplete || n == 0U) {
      return status;
    }
    DeviceSpan<Index> candidates;
    if (!candidates.allocate(mArena, n)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(*mRuntime, mStream, n, algorithms::Gather<Index>{std::as_const(sortedListed).data(), std::as_const(firsts).data(), candidates.data()});
    // The claims so far, ordered by track, and each candidate's range of them.
    const auto heldCount = static_cast<Index>(mClaimChunks.rows());
    DeviceColumns<PartnerClaimColumns::template RecordT> held;
    DeviceSpan<Index> heldOrder;
    DeviceSpan<Index> heldSorted;
    DeviceSpan<Index> heldTracks;
    DeviceSpan<Index> heldRows;
    DeviceSpan<Scalar> heldPaths;
    DeviceSpan<Index> heldBegin;
    DeviceSpan<Index> heldEnd;
    const std::size_t heldSize = std::max<std::size_t>(heldCount, 1U);
    if (!held.allocate(mArena, heldSize) || !heldOrder.allocate(mArena, heldSize) || !heldSorted.allocate(mArena, heldSize) || !heldTracks.allocate(mArena, heldSize) || !heldRows.allocate(mArena, heldSize) || !heldPaths.allocate(mArena, heldSize) || !heldBegin.allocate(mArena, n) || !heldEnd.allocate(mArena, n)) {
      return StageStatus::kOutOfMemory;
    }
    if (heldCount != 0U) {
      mClaimChunks.gather(held, 0U, mStream);
      const auto view = std::as_const(held).view();
      algorithms::forEach(*mRuntime, mStream, heldCount, tracker_detail::Iota{heldOrder.data()});
      status = algorithms::sortPairs(*mRuntime, mArena, mStream, view.track, heldTracks.data(), std::as_const(heldOrder).data(), heldSorted.data(), heldCount, trackBits);
      if (status != StageStatus::kComplete) {
        return status;
      }
      algorithms::forEach(*mRuntime, mStream, heldCount, algorithms::Gather<Index>{view.row, std::as_const(heldSorted).data(), heldRows.data()});
      algorithms::forEach(*mRuntime, mStream, heldCount, algorithms::Gather<Scalar>{view.path, std::as_const(heldSorted).data(), heldPaths.data()});
    }
    algorithms::forEach(*mRuntime, mStream, n, tracker_detail::HeldRange{std::as_const(candidates).data(), std::as_const(heldTracks).data(), heldCount, heldBegin.data(), heldEnd.data()});
    // The proposals.
    DeviceStageOutput<PartnerProposal::template RecordT> proposals;
    if (!proposals.prepare(mArena, n)) {
      return StageStatus::kOutOfMemory;
    }
    mPartnerStage.prepareStaging(std::max<std::size_t>(mPartnerStage.stagingCapacity(), n));
    const ProposePartners<Det> propose{.parameters = mParameters.data() + mIteration, .frame = mFrame.data(), .geometry = mGeometry.data(), .propagator = mPropagator, .cells = mSensorCells, .tracks = std::as_const(all).view(), .candidates = std::as_const(candidates).data(), .heldBegin = std::as_const(heldBegin).data(), .heldEnd = std::as_const(heldEnd).data(), .heldRows = std::as_const(heldRows).data(), .used = std::as_const(working).data(), .maxWidthPS = std::as_const(mLayerMaxWidthPS).data(), .measurementVariance = std::as_const(mLayerMaxLocVariance).data(), .footprints = std::as_const(mCellFootprints).data()};
    status = runStage(mPartnerStage, proposals, propose);
    if (status != StageStatus::kComplete) {
      return status;
    }
    const auto total = proposals.size();
    if (total == 0U) {
      return StageStatus::kComplete;
    }
    auto scratch = mArena.scope(ArenaEnd::kScratch);
    const auto proposed = std::as_const(proposals).view();
    // By row, then chi2: a stable sort keeps the candidate order among equals.
    DeviceSpan<std::uint64_t> keys;
    DeviceSpan<std::uint64_t> sortedKeys;
    DeviceSpan<Index> order;
    DeviceSpan<Index> sorted;
    DeviceSpan<std::uint8_t> won;
    DeviceSpan<Index> winnerOffsets;
    DeviceSpan<Index> winnerTracks;
    DeviceSpan<Index> winners;
    if (!keys.allocate(mArena, total, ArenaEnd::kScratch) || !sortedKeys.allocate(mArena, total, ArenaEnd::kScratch) || !order.allocate(mArena, total, ArenaEnd::kScratch) || !sorted.allocate(mArena, total, ArenaEnd::kScratch) || !won.allocate(mArena, total, ArenaEnd::kScratch) || !winnerOffsets.allocate(mArena, std::size_t{n} + 1U, ArenaEnd::kScratch)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(*mRuntime, mStream, total, tracker_detail::PartnerKeys{proposed.row, proposed.chi2, keys.data(), order.data()});
    status = algorithms::sortPairs(*mRuntime, mArena, mStream, std::as_const(keys).data(), sortedKeys.data(), std::as_const(order).data(), sorted.data(), total, 64);
    if (status != StageStatus::kComplete) {
      return status;
    }
    // The first of each row wins; the winners per candidate, in proposal order.
    algorithms::forEach(*mRuntime, mStream, total, tracker_detail::FirstOfRow{std::as_const(sortedKeys).data(), std::as_const(sorted).data(), won.data()});
    Index winnerCount = 0U;
    status = algorithms::selectIndices(*mRuntime, mArena, mStream, total, tracker_detail::Flagged{std::as_const(won).data()}, winners, winnerCount, ArenaEnd::kScratch);
    ++report.readbacks;
    if (status != StageStatus::kComplete) {
      return status;
    }
    if (!winnerTracks.allocate(mArena, std::max<Index>(winnerCount, 1U), ArenaEnd::kScratch)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(*mRuntime, mStream, winnerCount, algorithms::Gather<Index>{proposed.track, std::as_const(winners).data(), winnerTracks.data()});
    algorithms::forEach(*mRuntime, mStream, n + 1U, tracker_detail::LowerBounds{std::as_const(winnerTracks).data(), winnerCount, winnerOffsets.data()});
    // Each candidate with winners: its claims and winners, fitted together.
    DeviceSpan<Index> entryCounts;
    DeviceSpan<Index> entryOffsets;
    DeviceSpan<Index> entryTotal;
    if (!entryCounts.allocate(mArena, std::size_t{n} + 1U, ArenaEnd::kScratch) || !entryOffsets.allocate(mArena, std::size_t{n} + 1U, ArenaEnd::kScratch) || !entryTotal.allocate(mArena, 1U, ArenaEnd::kScratch)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(*mRuntime, mStream, n + 1U, tracker_detail::ZeroCounts{entryCounts.data()});
    algorithms::forEach(*mRuntime, mStream, n, tracker_detail::EntryCounts{std::as_const(winnerOffsets).data(), std::as_const(heldBegin).data(), std::as_const(heldEnd).data(), entryCounts.data()});
    status = algorithms::exclusiveScan(*mRuntime, mArena, mStream, std::as_const(entryCounts).data(), entryOffsets.data(), std::size_t{n} + 1U);
    if (status != StageStatus::kComplete) {
      return status;
    }
    Index entries = 0U;
    algorithms::forEach(*mRuntime, mStream, 1U, tracker_detail::CopyAt{std::as_const(entryOffsets).data() + n, entryTotal.data()});
    gtrk::download(entryTotal, &entries, 1U, mStream);
    mStream.synchronize();
    ++report.readbacks;
    if (entries == 0U) {
      return StageStatus::kComplete;
    }
    DeviceSpan<Index> entryRows;
    DeviceSpan<Scalar> entryPaths;
    DeviceSpan<std::uint8_t> fresh;
    DeviceSpan<std::uint8_t> layerRowOf;
    DeviceSpan<std::uint32_t> merged;
    DeviceSpan<std::uint8_t> mergedRow;
    DeviceSpan<std::uint8_t> kept;
    const std::size_t mergedSize = (std::size_t{n} * kMaxClusters) + entries;
    if (!entryRows.allocate(mArena, entries, ArenaEnd::kScratch) || !entryPaths.allocate(mArena, entries, ArenaEnd::kScratch) || !fresh.allocate(mArena, entries, ArenaEnd::kScratch) || !layerRowOf.allocate(mArena, entries, ArenaEnd::kScratch) || !merged.allocate(mArena, mergedSize, ArenaEnd::kScratch) || !mergedRow.allocate(mArena, mergedSize, ArenaEnd::kScratch) || !kept.allocate(mArena, n, ArenaEnd::kScratch)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(*mRuntime, mStream, n, tracker_detail::FillEntries{std::as_const(winnerOffsets).data(), std::as_const(winners).data(), proposed.row, proposed.path, std::as_const(heldBegin).data(), std::as_const(heldEnd).data(), std::as_const(heldRows).data(), std::as_const(heldPaths).data(), std::as_const(entryOffsets).data(), entryRows.data(), entryPaths.data(), fresh.data()});
    const FitWithPartners<Det, false> fit{.parameters = std::as_const(mParameters).data(), .frame = mFrame.data(), .geometry = mGeometry.data(), .propagator = mPropagator, .tracks = all.view(), .listed = std::as_const(candidates).data(), .offsets = std::as_const(entryOffsets).data(), .rows = entryRows.data(), .paths = entryPaths.data(), .fresh = fresh.data(), .layerRowOf = layerRowOf.data(), .merged = merged.data(), .mergedRow = mergedRow.data(), .kept = kept.data()};
    algorithms::forEach(*mRuntime, mStream, n, fit);
    // The kept candidates' winners claimed, and the used mask.
    DeviceSpan<Index> claimCounts;
    DeviceSpan<Index> claimOffsets;
    if (!claimCounts.allocate(mArena, std::size_t{n} + 1U, ArenaEnd::kScratch) || !claimOffsets.allocate(mArena, std::size_t{n} + 1U, ArenaEnd::kScratch)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(*mRuntime, mStream, n + 1U, tracker_detail::ZeroCounts{claimCounts.data()});
    algorithms::forEach(*mRuntime, mStream, n, tracker_detail::ClaimCounts{std::as_const(winnerOffsets).data(), std::as_const(kept).data(), claimCounts.data()});
    status = algorithms::exclusiveScan(*mRuntime, mArena, mStream, std::as_const(claimCounts).data(), claimOffsets.data(), std::size_t{n} + 1U);
    if (status != StageStatus::kComplete) {
      return status;
    }
    algorithms::forEach(*mRuntime, mStream, 1U, tracker_detail::CopyAt{std::as_const(claimOffsets).data() + n, entryTotal.data()});
    gtrk::download(entryTotal, &claimed, 1U, mStream);
    mStream.synchronize();
    ++report.readbacks;
    if (claimed == 0U) {
      return StageStatus::kComplete;
    }
    if (!claims.allocate(mArena, claimed)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(*mRuntime, mStream, n, tracker_detail::PlaceClaims{std::as_const(candidates).data(), std::as_const(kept).data(), std::as_const(entryOffsets).data(), std::as_const(entryRows).data(), std::as_const(entryPaths).data(), std::as_const(fresh).data(), std::as_const(claimOffsets).data(), claims.view(), working.data()});
    return StageStatus::kComplete;
  }

  /// The join (the host's `TrackJoin::run`) of the \a count tracks at
  /// \a earlier in \a all -- the iteration's, after the earlier iterations'
  /// -- to those before them: the earlier ones rewritten in \a all, the
  /// pieces joined removed from \a tracks, their proposals (\a followed) and
  /// from \a all (\a count then the pieces kept), the earlier tracks changed
  /// into \a changed (ascending), their number into \a joined.
  [[nodiscard]] StageStatus joinTracks(DeviceColumns<TrackRecord::template RecordT>& all, Index earlier, DeviceColumns<TrackRecord::template RecordT>& tracks, DeviceSpan<Proposal>& proposals, Index& count, bool followed, DeviceSpan<Index>& changed, Index& joined)
  {
    using View = typename DeviceColumns<TrackRecord::template RecordT>::View;
    auto& report = mIterationReport;
    joined = 0U;
    const auto constView = std::as_const(all).view();
    // The earlier tracks' ends and starts, by row then track: their rows
    // sorted stably from ascending tracks.
    DeviceSpan<Index> listed;
    Index ends = 0U;
    auto status = algorithms::selectIndices(*mRuntime, mArena, mStream, earlier, tracker_detail::HasRows{constView.status, constView.count}, listed, ends);
    ++report.readbacks;
    if (status != StageStatus::kComplete || ends == 0U) {
      return status;
    }
    DeviceSpan<std::uint32_t> endRows;
    DeviceSpan<std::uint32_t> startRows;
    DeviceSpan<std::uint32_t> sortedRows;
    DeviceSpan<std::uint32_t> sortedTracks;
    DeviceSpan<join::End> endList;
    DeviceSpan<join::End> startList;
    if (!endRows.allocate(mArena, ends) || !startRows.allocate(mArena, ends) || !sortedRows.allocate(mArena, ends) || !sortedTracks.allocate(mArena, ends) || !endList.allocate(mArena, ends) || !startList.allocate(mArena, ends)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(*mRuntime, mStream, ends, tracker_detail::EndRows<kMaxClusters>{constView.clusters, constView.count, std::as_const(listed).data(), endRows.data(), startRows.data()});
    const int rowBits = algorithms::keyBits(static_cast<std::uint32_t>(mRows));
    for (const bool outermost : {true, false}) {
      status = algorithms::sortPairs(*mRuntime, mArena, mStream, std::as_const(outermost ? endRows : startRows).data(), sortedRows.data(), std::as_const(listed).data(), sortedTracks.data(), ends, rowBits);
      if (status != StageStatus::kComplete) {
        return status;
      }
      algorithms::forEach(*mRuntime, mStream, ends, tracker_detail::ZipEnds{std::as_const(sortedRows).data(), std::as_const(sortedTracks).data(), (outermost ? endList : startList).data()});
    }
    // Each piece's pair.
    DeviceSpan<join::Choice> choices;
    if (!choices.allocate(mArena, count)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(*mRuntime, mStream, count, join::Choose<Det, View>{.parameters = std::as_const(mParameters).data(), .frame = std::as_const(mFrame).data(), .geometry = std::as_const(mGeometry).data(), .propagator = mPropagator, .tracks = all.view(), .trackOffset = earlier, .ends = std::as_const(endList).data(), .starts = std::as_const(startList).data(), .endCount = ends, .maxCandidateWidthPS = std::as_const(mLayerMaxWidthPS).data(), .choices = choices.data()});
    // Each earlier track the piece of the least score that chose it: the
    // pieces with a partner by partner, then score, stably from ascending
    // pieces, the first of each partner's run.
    DeviceSpan<Index> chose;
    Index chosen = 0U;
    status = algorithms::selectIndices(*mRuntime, mArena, mStream, count, tracker_detail::ChosePartner{std::as_const(choices).data()}, chose, chosen);
    ++report.readbacks;
    if (status != StageStatus::kComplete || chosen == 0U) {
      return status;
    }
    DeviceSpan<std::uint64_t> keys;
    DeviceSpan<std::uint64_t> sortedKeys;
    DeviceSpan<std::uint32_t> pieces;
    DeviceSpan<std::uint32_t> sortedPieces;
    DeviceSpan<std::uint8_t> joinedFlags;
    if (!keys.allocate(mArena, chosen) || !sortedKeys.allocate(mArena, chosen) || !pieces.allocate(mArena, chosen) || !sortedPieces.allocate(mArena, chosen) || !joinedFlags.allocate(mArena, count)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(*mRuntime, mStream, chosen, tracker_detail::JoinKeys{std::as_const(choices).data(), std::as_const(chose).data(), keys.data(), pieces.data()});
    status = algorithms::sortPairs(*mRuntime, mArena, mStream, std::as_const(keys).data(), sortedKeys.data(), std::as_const(pieces).data(), sortedPieces.data(), chosen, 64);
    if (status != StageStatus::kComplete) {
      return status;
    }
    mRuntime->fill(joinedFlags.data(), 0U, count, mStream);
    algorithms::forEach(*mRuntime, mStream, chosen, tracker_detail::FirstOfRow{std::as_const(sortedKeys).data(), std::as_const(sortedPieces).data(), joinedFlags.data()});
    DeviceSpan<Index> joinedPieces;
    status = algorithms::selectIndices(*mRuntime, mArena, mStream, count, tracker_detail::Flagged{std::as_const(joinedFlags).data()}, joinedPieces, joined);
    ++report.readbacks;
    if (status != StageStatus::kComplete || joined == 0U) {
      return status;
    }
    // The joined tracks over the earlier ones, those listed ascending.
    algorithms::forEach(*mRuntime, mStream, joined, join::Apply<Det, View>{.parameters = std::as_const(mParameters).data(), .frame = std::as_const(mFrame).data(), .geometry = std::as_const(mGeometry).data(), .propagator = mPropagator, .tracks = all.view(), .trackOffset = earlier, .choices = std::as_const(choices).data(), .joined = std::as_const(joinedPieces).data()});
    DeviceSpan<Index> partners;
    if (!partners.allocate(mArena, joined) || !changed.allocate(mArena, joined)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(*mRuntime, mStream, joined, tracker_detail::PartnersOf{std::as_const(choices).data(), std::as_const(joinedPieces).data(), partners.data()});
    status = algorithms::sortKeys(*mRuntime, mArena, mStream, std::as_const(partners).data(), changed.data(), joined, algorithms::keyBits(earlier));
    if (status != StageStatus::kComplete) {
      return status;
    }
    // The pieces joined removed; their rows stay used, the earlier tracks'.
    DeviceSpan<Index> kept;
    Index keptCount = 0U;
    status = algorithms::selectIndices(*mRuntime, mArena, mStream, count, tracker_detail::Unflagged{std::as_const(joinedFlags).data()}, kept, keptCount);
    ++report.readbacks;
    if (status != StageStatus::kComplete) {
      return status;
    }
    DeviceColumns<TrackRecord::template RecordT> keptTracks;
    DeviceSpan<Proposal> keptProposals;
    if (!keptTracks.allocate(mArena, keptCount) || !keptProposals.allocate(mArena, followed ? keptCount : 0U)) {
      return StageStatus::kOutOfMemory;
    }
    if (keptCount != 0U) {
      device_storage_detail::forEachColumnPair(std::as_const(tracks).view(), keptTracks.view(), [&](auto* to, const auto* from) { algorithms::forEach(*mRuntime, mStream, keptCount, algorithms::Gather<std::remove_pointer_t<decltype(to)>>{from, std::as_const(kept).data(), to}); });
      if (followed) {
        algorithms::forEach(*mRuntime, mStream, keptCount, algorithms::Gather<Proposal>{std::as_const(proposals).data(), std::as_const(kept).data(), keptProposals.data()});
      }
      copyColumns(keptTracks, 0U, all, earlier, keptCount, mStream);
    }
    tracks = keptTracks;
    proposals = keptProposals;
    count = keptCount;
    return StageStatus::kComplete;
  }

  /// The inactive layers' rows (the host's `InactiveRows::run`) with the
  /// iteration's \a parameters for the \a changedCount tracks of \a view in
  /// \a changed (ascending), then its \a ownCount from \a own: the chosen
  /// rows marked in \a working, the tracks that took one rewritten, their
  /// number into \a taken.
  template <typename View>
  [[nodiscard]] StageStatus inactiveRows(const IterationParameters<Det>* parameters, View view, Index own, Index ownCount, const DeviceSpan<Index>& changed, Index changedCount, DeviceSpan<std::uint8_t>& working, Index& taken)
  {
    auto& report = mIterationReport;
    taken = 0U;
    const Index count = changedCount + ownCount;
    DeviceSpan<Index> candidates;
    DeviceSpan<inactive::Choice> choices;
    if (!candidates.allocate(mArena, count) || !choices.allocate(mArena, count)) {
      return StageStatus::kOutOfMemory;
    }
    if (changedCount != 0U) {
      mRuntime->copyOnDevice(candidates.data(), changed.data(), std::size_t{changedCount} * sizeof(Index), mStream);
    }
    if (ownCount != 0U) {
      algorithms::forEach(*mRuntime, mStream, ownCount, tracker_detail::Sequence{own, candidates.data() + changedCount});
    }
    algorithms::forEach(*mRuntime, mStream, count, inactive::Choose<Det, View>{.iteration = parameters, .timeFrame = std::as_const(mFrame).data(), .trackingGeometry = std::as_const(mGeometry).data(), .propagator = mPropagator, .tracks = view, .candidates = std::as_const(candidates).data(), .used = std::as_const(working).data(), .maxCandidateWidthPS = std::as_const(mLayerMaxWidthPS).data(), .cells = mSensorCells, .choices = choices.data()});
    // A row two tracks want to the smaller rise, then the lower track: the
    // choices by row, then rise, stably from ascending tracks, the first of
    // each row's run.
    DeviceSpan<Index> chose;
    Index chosen = 0U;
    auto status = algorithms::selectIndices(*mRuntime, mArena, mStream, count, tracker_detail::ChoseRow{std::as_const(choices).data()}, chose, chosen);
    ++report.readbacks;
    if (status != StageStatus::kComplete || chosen == 0U) {
      return status;
    }
    DeviceSpan<std::uint64_t> keys;
    DeviceSpan<std::uint64_t> sortedKeys;
    DeviceSpan<std::uint32_t> positions;
    DeviceSpan<std::uint32_t> sortedPositions;
    DeviceSpan<std::uint8_t> winners;
    if (!keys.allocate(mArena, chosen) || !sortedKeys.allocate(mArena, chosen) || !positions.allocate(mArena, chosen) || !sortedPositions.allocate(mArena, chosen) || !winners.allocate(mArena, count)) {
      return StageStatus::kOutOfMemory;
    }
    algorithms::forEach(*mRuntime, mStream, chosen, tracker_detail::InactiveKeys{std::as_const(choices).data(), std::as_const(chose).data(), keys.data(), positions.data()});
    status = algorithms::sortPairs(*mRuntime, mArena, mStream, std::as_const(keys).data(), sortedKeys.data(), std::as_const(positions).data(), sortedPositions.data(), chosen, 64);
    if (status != StageStatus::kComplete) {
      return status;
    }
    mRuntime->fill(winners.data(), 0U, count, mStream);
    algorithms::forEach(*mRuntime, mStream, chosen, tracker_detail::FirstOfRow{std::as_const(sortedKeys).data(), std::as_const(sortedPositions).data(), winners.data()});
    DeviceSpan<Index> order;
    status = algorithms::selectIndices(*mRuntime, mArena, mStream, count, tracker_detail::Flagged{std::as_const(winners).data()}, order, taken);
    ++report.readbacks;
    if (status != StageStatus::kComplete || taken == 0U) {
      return status;
    }
    algorithms::forEach(*mRuntime, mStream, taken, inactive::Apply<Det, View>{.iteration = parameters, .timeFrame = std::as_const(mFrame).data(), .trackingGeometry = std::as_const(mGeometry).data(), .propagator = mPropagator, .tracks = view, .candidates = std::as_const(candidates).data(), .choices = std::as_const(choices).data(), .order = std::as_const(order).data(), .used = working.data()});
    return StageStatus::kComplete;
  }

  /// The recovery pass of iteration \a successor over every track so far --
  /// the chunks' and this iteration's \a tracks (\a trackCount of them) --
  /// with \a working as the used mask, its bids in \a slots (one per row,
  /// \a rows): the host's `Follow::run` with `earlier`. The tracks, updated,
  /// replace the chunks as one.
  [[nodiscard]] StageStatus followEarlier(std::size_t successor, const DeviceColumns<TrackRecord::template RecordT>& tracks, Index trackCount, DeviceSpan<std::uint8_t>& working, DeviceSpan<Index>& slots, Index rows, const DeviceColumns<TrackRecord::template RecordT>* gathered)
  {
    const IterationParameters<Det>& host = mHostParameters[successor];
    const auto earlier = static_cast<Index>(mTrackChunks.rows());
    const Index total = earlier + trackCount;
    DeviceColumns<TrackRecord::template RecordT> all;
    DeviceSpan<Proposal> proposals;
    DeviceSpan<ClusterList<kMaxClusters>> followRows;
    DeviceSpan<std::uint8_t> followCounts;
    if (!all.allocate(mArena, total) || !proposals.allocate(mArena, total) || !followRows.allocate(mArena, total) || !followCounts.allocate(mArena, total)) {
      return StageStatus::kOutOfMemory;
    }
    if (gathered != nullptr) {
      copyColumns(*gathered, 0U, all, 0U, earlier, mStream);
    } else {
      mTrackChunks.gather(all, 0U, mStream);
    }
    copyColumns(tracks, 0U, all, earlier, trackCount, mStream);
    runFollow(successor, all, total, true, proposals.data(), followRows, followCounts, slots, rows, working);
    if constexpr (overlapLayers<Det> != 0U) {
      // The tracks it changed: the successor's overlap partners' candidates.
      DeviceSpan<Index> changed;
      Index changedCount = 0U;
      const auto selected = algorithms::selectIndices(*mRuntime, mArena, mStream, total, tracker_detail::FollowChanged<Proposal>{std::as_const(proposals).data()}, changed, changedCount);
      ++mIterationReport.readbacks;
      if (selected != StageStatus::kComplete) {
        return selected;
      }
      DeviceColumns<tracker_detail::TrackIndexColumns::template RecordT> carried;
      if (changedCount != 0U) {
        if (!carried.allocate(mArena, changedCount)) {
          return StageStatus::kOutOfMemory;
        }
        mRuntime->copyOnDevice(carried.view().track, std::as_const(changed).data(), std::size_t{changedCount} * sizeof(Index), mStream);
        mCarriedChunks.clear();
        if (!mCarriedChunks.append(mArena, ArenaEnd::kScratch, carried, 0U, changedCount, mStream)) {
          return StageStatus::kOutOfMemory;
        }
      }
    }
    mTrackChunks.clear();
    if (!mTrackChunks.append(mArena, ArenaEnd::kScratch, all, 0U, total, mStream)) {
      return StageStatus::kOutOfMemory;
    }
    return StageStatus::kComplete;
  }

  /// The follow of iteration \a iteration over \a count tracks of \a tracks
  /// (`follow::runSequence`, as the host's `Follow::run`): \a earlier the
  /// earlier iterations' tracks (`followEarlier`), the proposals into
  /// \a settled, the combined rows into \a followRows and \a followCounts, the
  /// bids in \a slots (one per row, \a rows of them), the marks in
  /// \a working.
  void runFollow(std::size_t iteration, DeviceColumns<TrackRecord::template RecordT>& tracks, Index count, bool earlier, Proposal* settled, DeviceSpan<ClusterList<kMaxClusters>>& followRows, DeviceSpan<std::uint8_t>& followCounts, DeviceSpan<Index>& slots, Index rows, DeviceSpan<std::uint8_t>& working)
  {
    const IterationParameters<Det>& host = mHostParameters[iteration];
    follow::runSequence(tracker_detail::FollowLaunch{.runtime = *mRuntime, .stream = mStream}, follow::Inputs<Det>{
                                                                                                 .parameters = mParameters.data() + iteration,
                                                                                                 .choices = &host,
                                                                                                 .frame = mFrame.data(),
                                                                                                 .geometry = mGeometry.data(),
                                                                                                 .propagator = mPropagator,
                                                                                                 .tracks = tracks.view(),
                                                                                                 .constTracks = std::as_const(tracks).view(),
                                                                                                 .trackOffset = 0U,
                                                                                                 .trackCount = count,
                                                                                                 .inward = host.sharedInnermost == 0U,
                                                                                                 .earlier = earlier,
                                                                                                 .used = working.data(),
                                                                                                 .rowCount = rows,
                                                                                                 .maxCandidateWidthPS = std::as_const(mLayerMaxWidthPS).data(),
                                                                                                 .cells = mSensorCells,
                                                                                                 .proposals = settled,
                                                                                                 .bids = slots.data(),
                                                                                                 .rows = followRows.data(),
                                                                                                 .counts = followCounts.data(),
                                                                                               });
  }

  /// `requireTimeframe`, and no iteration open: \a call names the caller.
  void requireNoIteration(const char* call) const
  {
    requireTimeframe();
    if (mIterationOpen) {
      throw std::logic_error(std::string{call} + " with an iteration open");
    }
  }

  /// A slice's scope on the arena, the previous one closed and the report
  /// cleared.
  void openSlice()
  {
    // Every call before returns with the stream drained.
    closeFrom(Scope::kSlice);
    mSliceScope.emplace(mArena.scope());
    mReport = {};
  }

  void requireTimeframe() const
  {
    if (!mTimeframeScope) {
      throw std::logic_error("no timeframe has been uploaded to the device");
    }
  }

  /// The tracker's scopes: the run is the
  /// arena itself; the timeframe holds, on the base end, the slice or the
  /// iteration's -- never both at once -- and, on the scratch end, the
  /// tracks and above them the open iteration's record chunks.
  enum class Scope : std::uint8_t { kRun,
                                    kTimeframe,
                                    kChunks,
                                    kIteration,
                                    kSlice };

  /// Close \a scope and every scope inside it, innermost first: the only
  /// place scopes are closed. The slice goes with each; the chunks and the
  /// iteration's scope go with their own call or with the timeframe, the
  /// tracks with the timeframe.
  void closeFrom(Scope scope) noexcept
  {
    mSliceScope.reset();
    switch (scope) {
      case Scope::kSlice:
        return;
      case Scope::kChunks:
        closeChunks();
        return;
      case Scope::kIteration:
        mIterationScope.reset();
        mGathered = {};
        return;
      case Scope::kTimeframe:
      case Scope::kRun:
        break;
    }
    mIterationScope.reset();
    mGathered = {};
    closeChunks();
    mTrackChunks.clear();
    mClaimChunks.clear();
    mCarriedChunks.clear();
    mFinishedListed.clear();
    mFinishedOffsets.clear();
    mFinishedPartners.clear();
    mTracksScope.reset();
    mTimeframeScope.reset();
    mRows = 0U;
    mHostFrame = nullptr;
    mPlan = nullptr;
    if (scope == Scope::kTimeframe) {
      return;
    }
    mArena.reset();
    mConfigured = false;
    mIterationCount = 0U;
    mParameters = {};
    mHostParameters.clear();
    mGeometry = {};
    mCellLayers = {};
    mCellOffsets = {};
    mCellSensors = {};
    mCellFootprints = {};
    mSensorCells = {};
  }

  void closeChunks() noexcept
  {
    mRoadChunks.clear();
    mBatchScope.reset();
    mChunks.clear();
    mChunksScope.reset();
    mIterationOpen = false;
  }

  DeviceRuntime* mRuntime = nullptr;
  DeviceArena mArena;
  DeviceStream mStream;
  GpuStage<TrackletT> mTrackletStage;
  GpuStage<NeighbourT> mNeighbourStage;
  GpuStage<SelectionRecordOf<kMaxClusters>::template RecordT> mRecordStage;
  GpuStage<PartnerProposal::template RecordT> mPartnerStage;
  DeviceColumnChunks<tracker_detail::CapturedRoadOf<kMaxClusters>::template RecordT> mRoadChunks; ///< the open batch's roads, a chunk per slice
  std::size_t mRoadBatch = 0U;                                                                    ///< the roads that fill the device: a batch's fit is due at them
  DeviceColumnChunks<SelectionRecordOf<kMaxClusters>::template RecordT> mChunks;
  DeviceColumns<SelectionRecordOf<kMaxClusters>::template RecordT> mGathered;
  DeviceColumnChunks<TrackRecord::template RecordT> mTrackChunks;
  /// Beside them, each track's overlap partners: their count per track, and
  /// the partners, track after track.
  DeviceColumnChunks<PartnerClaimColumns::template RecordT> mClaimChunks; ///< the overlap partners claimed, each iteration's a chunk
  /// The earlier tracks the follow of earlier tracks changed at an
  /// iteration's end: the next iteration's overlap partners' candidates.
  DeviceColumnChunks<tracker_detail::TrackIndexColumns::template RecordT> mCarriedChunks;
  /// `finishPartners`' tracks with partners, their partners' ranges and the
  /// partners, for `downloadPartners`.
  std::vector<Index> mFinishedListed;
  std::vector<Index> mFinishedOffsets;
  std::vector<OverlapPartner> mFinishedPartners;
  IterationReport mIterationReport;
  std::size_t mOuterArenaPeak = 0U; ///< the arena's peak before the open iteration's window
  flat::FlatOptions mFlatOptions;   ///< the flat runs' test hook (`setFlatChunkPairs`)
  Report mReport;

  DeviceSpan<IterationParameters<Det>> mParameters;
  std::vector<IterationParameters<Det>> mHostParameters; ///< the same, for the values the host reads
  DeviceObject<TrackingGeometry<Det>> mGeometry;
  DeviceSpan<SensorCellLayer> mCellLayers; ///< the sensors by place, for the outward follow
  DeviceSpan<std::uint32_t> mCellOffsets;
  DeviceSpan<std::uint32_t> mCellSensors;
  DeviceSpan<SensorFootprint> mCellFootprints; ///< per sensor, for the overlap partners' search
  SensorCellsView mSensorCells;
  /// The sensor cells on the host and the geometry they were built from.
  std::optional<SensorCells<Det>> mHostCells;
  std::array<std::byte, sizeof(TrackingGeometry<Det>)> mHostCellsGeometry{};
  Propagator mPropagator;
  std::size_t mIterationCount = 0U;
  bool mConfigured = false;

  std::optional<ArenaScope> mTimeframeScope;
  DeviceColumns<MeasurementT> mColumns;
  DeviceSpan<Index> mOffsets;
  DeviceSpan<std::uint8_t> mUsed;
  DeviceObject<Frame> mFrame;
  std::size_t mRows = 0U;
  const Frame* mHostFrame = nullptr;
  const TimeframePlan<Det>* mPlan = nullptr;
  DeviceSpan<LayerRowMap<Det>> mSliceSources;
  DeviceSpan<MeasurementRange> mSliceCores;
  DeviceSpan<MeasurementRange> mSliceRoadCores; ///< mSliceCores with each iteration's start-layer mask, where one has a mask
  bool mRoadCoresMasked = false;
  DeviceSpan<VertexWindowParameters<Det>> mVertexWindows;
  DeviceSpan<VertexSeed> mSeeds; ///< the beam line, then the vertices
  DeviceSpan<std::uint32_t> mLayerMaxWidthPS;
  DeviceSpan<Scalar> mLayerMaxLocVariance;       ///< the plan's, over the whole frame
  DeviceSpan<TripletTargetBounds> mTargetBounds; ///< the plan's, over the whole frame
  std::optional<ArenaScope> mTracksScope;

  std::optional<ArenaScope> mChunksScope;
  std::optional<ArenaScope> mBatchScope; ///< the open batch's roads, inside the chunks
  std::optional<ArenaScope> mIterationScope;
  std::optional<ArenaScope> mSliceScope;
  std::size_t mIteration = 0U;
  StageStatus mIterationStatus = StageStatus::kComplete;
  bool mIterationOpen = false;
};

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE

#endif // GTRK_TRACKING_GPU_DEVICETRACKER_CUH_
