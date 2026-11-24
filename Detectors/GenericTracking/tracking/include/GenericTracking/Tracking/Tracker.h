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

/// \file Tracker.h
/// \brief The timeframe's tracking: every iteration over every slice, the
///        global selection between them, and the refit of the tracks.

#ifndef GTRK_TRACKING_TRACKER_H_
#define GTRK_TRACKING_TRACKER_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/CapacityEstimator.h"
#include "GenericTracking/Base/HostAlgorithms.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/InactiveRows.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/NeighbourStage.h"
#include "GenericTracking/Tracking/OverlapPartners.h"
#include "GenericTracking/Tracking/Parallel.h"
#include "GenericTracking/Tracking/Profile.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/RoadCount.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/Selection.h"
#include "GenericTracking/Tracking/SensorCells.h"
#include "GenericTracking/Tracking/SlicePlan.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TimeFrameScratch.h"
#include "GenericTracking/Tracking/TrackJoin.h"
#include "GenericTracking/Tracking/TrackRefit.h"
#include "GenericTracking/Tracking/TrackStore.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/TrackletSort.h"
#include "GenericTracking/Tracking/TrackletStage.h"
#include "GenericTracking/Tracking/TripletStage.h"

#include "GPUCommonAlgorithm.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace o2::gtrk
{

/// The stages whose per-thread staging the `CapacityEstimator` predicts, its
/// site indices.
enum class StagingSite : std::uint8_t {
  kTracklets = 0,
  kTriplets,
  kNeighbours,
  kRecords,
};
static_assert(static_cast<std::size_t>(StagingSite::kRecords) + 1U == CapacityEstimator::kSites);

/// The tracking of a timeframe:
/// every iteration, one after another, over every slice, the global
/// selection between them into the used mask, and the refit of the
/// selected tracks.
///
/// One slice of one iteration runs tracklets, triplets, neighbours and
/// levels (`runSlice`), then the roads and their fit. Every stage works the
/// core and whatever in the halo lies within the gate (`SliceReach::gatePS`)
/// of it, so the output the core owns does not depend on the slicing.
/// Sources and targets, and the
/// targets' variance bounds, come from the halo. The scratch and
/// backend retain their allocations across calls; a failed run returns its
/// StageStatus and no usable output.
template <DetectorLayout Det>
class Tracker
{
 public:
  using Frame = TimeFrame<Det>;
  using Scratch = TimeFrameScratch<Det>;
  using Tracklets = typename Scratch::Tracklets;
  using Triplets = typename Scratch::Triplets;
  using Neighbours = typename Scratch::Neighbours;
  using VertexWindows = VertexWindowParameters<Det>;
  using Index = typename Tracklets::Index;
  static constexpr std::size_t kLayerCount = Scratch::kLayerCount;
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using Records = StageOutput<SelectionRecordOf<kMaxClusters>::template RecordT>;

  /// A slice's tracklets, the triplets built from them with their levels and
  /// each triplet's predecessors, all valid until the next run.
  struct Result {
    StageStatus status = StageStatus::kComplete;
    const Tracklets* tracklets = nullptr;
    const Triplets* triplets = nullptr;
    const Neighbours* neighbours = nullptr;
  };

  using Phase = profile::Phase;

  /// What one iteration of a timeframe did.
  struct IterationSummary {
    StageStatus status = StageStatus::kComplete;
    std::size_t slices = 0U;
    std::size_t roads = 0U;      ///< fitted, over the slices
    std::size_t records = 0U;    ///< passing the fit, over the slices
    std::size_t duplicates = 0U; ///< records found by more than one slice
    std::size_t mismatches = 0U; ///< duplicates that differ: must be none
    std::size_t rounds = 0U;
    std::size_t tracks = 0U;
    /// Selected tracks whose two-pass refit failed, dropped after the follow
    /// (which skips them); `tracks` counts them.
    std::size_t refitFailures = 0U;
    /// The iteration's tracks its follow grew, and the rows it added outward
    /// and inward.
    std::size_t tracksFollowed = 0U;
    std::size_t clustersFollowedOutward = 0U;
    std::size_t clustersFollowedInward = 0U;
    /// With followEarlier, the earlier iterations' tracks its follow grew,
    /// and their rows.
    std::size_t tracksFollowedEarlier = 0U;
    std::size_t clustersFollowedEarlierOutward = 0U;
    std::size_t clustersFollowedEarlierInward = 0U;
    /// Tracks below the iteration's length after the follow, dropped;
    /// `tracks` counts them.
    std::size_t tracksDropped = 0U;
    /// The iteration's tracks joined to an earlier iteration's (TrackJoin.h),
    /// which takes their rows; `tracks` counts them.
    std::size_t tracksJoined = 0U;
    /// Tracks that took an inactive layer's row on their whole fit
    /// (InactiveRows.h).
    std::size_t inactiveRows = 0U;
    /// With overlapPartners, the partners attached and the tracks they went
    /// to (OverlapPartners.h).
    std::size_t overlapPartners = 0U;
    std::size_t tracksWithPartners = 0U;
    std::size_t recordBytes = 0U; ///< the selection records held for the iteration
    std::int64_t haloPS = 0;
    /// Every phase of the iteration, from `iteration` down; empty with
    /// profiling off.
    profile::Profile profile;

    /// The wall time of \a phase in the iteration, over all its slices.
    [[nodiscard]] double seconds(Phase phase) const { return profile[phase].seconds; }
  };

  /// The timeframe's tracks, all iterations' in iteration order, each
  /// iteration's best first, with each iteration's summary.
  struct TimeframeResult {
    StageStatus status = StageStatus::kComplete;
    /// The tracks committed (TrackStore.h): their fields, their references,
    /// and the tables other components add.
    const TrackStore<Det>* store = nullptr;
    std::vector<IterationSummary> iterations;
    /// The run's phases: `run`, `begin timeframe`, and the iterations' summed.
    profile::Profile profile;
  };

  /// What a run left after its first iterations, all that a later iteration
  /// reads: the tracks, the used
  /// mask and the iterations' summaries, their profiles empty (`resume`
  /// does not run them).
  struct Checkpoint {
    std::vector<typename TrackRefit<Det>::Tracks::Record> tracks;
    std::vector<std::uint8_t> used;
    std::vector<IterationSummary> summaries;
  };

  /// One tracker for every iteration of \a iterations, their parameters built
  /// once here.
  Tracker(std::span<const IterationSettings> iterations, const VertexConstraint& beamSpot, const LinkMaterial<Det>& material, std::shared_ptr<BoundedMemoryResource> memoryPool, Index chunkSize, int maxConcurrency = 0)
    : mBeamLine(vertexSeed(beamSpot)),
      mMemoryPool(std::move(memoryPool)),
      mArena(maxConcurrency),
      mTrackletStage(mArena),
      mTripletStage(mArena),
      mNeighbourStage(mArena),
      mRecordStage(mArena),
      mPartners(mArena),
      mChunkSize(chunkSize)
  {
    if (iterations.empty()) {
      throw std::invalid_argument("a tracker without iterations");
    }
    for (std::size_t iteration = 0; iteration < iterations.size(); ++iteration) {
      const auto& settings = iterations[iteration];
      // An earlier iteration's partners are references the follow of its
      // tracks would not carry along.
      if (iterations[iteration].followEarlier && std::ranges::any_of(iterations.first(iteration), [](const auto& earlier) { return earlier.overlapPartners; })) {
        throw std::invalid_argument("iteration " + std::to_string(iteration) + " follows earlier iterations' tracks after an iteration with overlap partners");
      }
      // An earlier iteration with the same windows has the same links'
      // windows and bounds: derived once (`iterationParameters`).
      const IterationParameters<Det>* alike = nullptr;
      for (std::size_t earlier = 0; earlier < iteration && alike == nullptr; ++earlier) {
        if (mWindowSettings[earlier] == settings.windows) {
          alike = &mIterations[earlier];
        }
      }
      mWindowSettings.push_back(settings.windows);
      mIterations.push_back(iterationParameters<Det>(settings, beamSpot, material, alike));
    }
    mScratch.setMemoryPool(mMemoryPool);
    mTrackletStage.setMemoryPool(mMemoryPool);
    mTripletStage.setMemoryPool(mMemoryPool);
    mNeighbourStage.setMemoryPool(mMemoryPool);
    mRecordStage.setMemoryPool(mMemoryPool);
    mRoads.setMemoryPool(mMemoryPool);
    mRecords.setMemoryPool(mMemoryPool);
    mTrackletStage.setPhases({.reset = Phase::kTrackletReset, .generate = Phase::kTrackletGenerate, .scan = Phase::kTrackletScan, .place = Phase::kTrackletPlace});
    mTripletStage.setPhases({.reset = Phase::kTripletReset, .generate = Phase::kTripletGenerate, .scan = Phase::kTripletScan, .place = Phase::kTripletPlace});
    mNeighbourStage.setPhases({.reset = Phase::kNeighbourReset, .generate = Phase::kNeighbourGenerate, .scan = Phase::kNeighbourScan, .place = Phase::kNeighbourPlace});
    mRecordStage.setPhases({.reset = Phase::kRoadFitReset, .generate = Phase::kRoadFitGenerate, .scan = Phase::kRoadFitScan, .place = Phase::kRoadFitPlace});
    mSelection.setMemoryPool(mMemoryPool);
    mRefit.setMemoryPool(mMemoryPool);
    mStore.setMemoryPool(mMemoryPool);
    mPartners.setMemoryPool(mMemoryPool);
    mFollow.setMemoryPool(mMemoryPool);
    mFollowEarlier.setMemoryPool(mMemoryPool);
    mJoin.setMemoryPool(mMemoryPool);
    mInactiveRows.setMemoryPool(mMemoryPool);
    deepVectorClear(mMemoryPool.get(), mUsed, mUsedBeforeSelection, mKeptTracks, mChangedTracks);
    mDropScratch.setMemoryResource(mMemoryPool.get());
  }

  /// A tracker of one iteration.
  Tracker(const IterationSettings& settings, const VertexConstraint& beamSpot, const LinkMaterial<Det>& material, std::shared_ptr<BoundedMemoryResource> memoryPool, Index chunkSize, int maxConcurrency = 0)
    : Tracker(std::span<const IterationSettings>{&settings, 1U}, beamSpot, material, std::move(memoryPool), chunkSize, maxConcurrency)
  {
  }

  /// What the road fit and the refit need: the sensors with their link
  /// material and the field. Required by `run(frame, ...)`, not by `runSlice`.
  void setTrackModel(const TrackingGeometry<Det>& geometry, const Propagator& propagator)
  {
    mGeometry = &geometry;
    mPropagator = propagator;
    // The sensors by place, once, where an iteration follows or takes
    // overlap partners.
    if (std::ranges::any_of(mIterations, [](const auto& parameters) { return parameters.follow || parameters.overlapPartners; })) {
      mSensorCells = SensorCells<Det>(geometry);
    }
  }

  [[nodiscard]] std::size_t iterationCount() const { return mIterations.size(); }
  [[nodiscard]] const IterationParameters<Det>& iteration(std::size_t index) const { return mIterations[index]; }

  /// Start a timeframe: every row free, and each stage's staging at what the
  /// estimator learned (`stageStaging`).
  [[nodiscard]] StageStatus beginTimeframe(const Frame& frame)
  {
    return guardOutOfMemory([&] {
      clearResizeBoundedVector(mUsed, frame.size(), mMemoryPool.get(), std::uint8_t{kFree});
      mRows = frame.size();
      stageStaging(StagingSite::kTracklets, mTrackletStage);
      stageStaging(StagingSite::kTriplets, mTripletStage);
      stageStaging(StagingSite::kNeighbours, mNeighbourStage);
      stageStaging(StagingSite::kRecords, mRecordStage);
      return StageStatus::kComplete;
    });
  }

  /// The end of a timeframe, after its caller has read `run`'s result:
  /// the estimator learns each
  /// stage's largest staging use, and everything that grew with the timeframe
  /// goes back to the pool -- what each iteration releases
  /// (`releaseIterationStorage`), the shared records, the tracks and the used
  /// masks. The result, its tracks included, is invalid after.
  void endTimeframe() noexcept
  {
    learnStaging(StagingSite::kTracklets, mTrackletStage);
    learnStaging(StagingSite::kTriplets, mTripletStage);
    learnStaging(StagingSite::kNeighbours, mNeighbourStage);
    learnStaging(StagingSite::kRecords, mRecordStage);
    releaseIterationStorage();
    mSelection.clearStorage();
    mRefit.clearStorage();
    mStore.clearStorage();
    deepVectorClear(mMemoryPool.get(), mUsed, mUsedBeforeSelection);
    mFollowStartInner = {};
    mFollowStart = {};
    mFollowStartCount = {};
    mFollowUsed = {};
    mFollowEndClusters = {};
    mFollowEndCount = {};
    mFollowEndState = {};
  }

  /// The committed tracks (TrackStore.h): a component that adds a table of
  /// its own (`addTable`) does so here, before the run.
  [[nodiscard]] TrackStore<Det>& trackStore() { return mStore; }
  [[nodiscard]] const TrackStore<Det>& trackStore() const { return mStore; }

  /// The tracks as the tracker builds and extends them, rows inline: what the
  /// store commits, and what a checkpoint holds. Valid as the store is.
  [[nodiscard]] const typename TrackRefit<Det>::Tracks& workingTracks() const { return mRefit.tracks(); }

  /// The used mask, per row: what the selection has taken so far.
  [[nodiscard]] std::span<const std::uint8_t> used() const { return mUsed; }

  /// Iteration \a iteration's graph on the used mask as it stands, without
  /// the selection, the refit or the follow (for the truth-removal study):
  /// every slice of \a plan, tracklets to the
  /// numbered roads -- with \a fit, also the road fit's records -- and \a
  /// visit(result, roads, records, fitSeconds) after each slice while its
  /// graph lives, records null and fitSeconds 0 without \a fit. After `beginTimeframe`. An iteration that
  /// follows the earlier tracks has no graph of its own: std::invalid_argument. Returns the first status short of complete.
  template <typename Visit>
  [[nodiscard]] StageStatus runGraph(std::size_t iteration, const Frame& frame, const TimeframePlan<Det>& plan, bool fit, const Visit& visit)
  {
    plan.check(frame, mIterations);
    if (mIterations[iteration].followEarlier) {
      throw std::invalid_argument("iteration " + std::to_string(iteration) + " follows the earlier tracks: no graph of its own");
    }
    for (const auto& slicePlan : plan.slices) {
      const auto sliced = runPlannedSlice(iteration, frame, plan, slicePlan, slicePlan.slice, plan.beamLine, plan.vertices, plan.vertexWindows[iteration], true);
      if (sliced.status != StageStatus::kComplete) {
        return sliced.status;
      }
      const auto counted = countRoads();
      if (counted != StageStatus::kComplete) {
        return counted;
      }
      double fitSeconds = 0.;
      if (fit) {
        const auto start = std::chrono::steady_clock::now();
        const auto fitted = guardOutOfMemory([&] { return fitRoads(frame, iteration, mRecords); });
        if (fitted != StageStatus::kComplete) {
          return fitted;
        }
        fitSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
      }
      visit(sliced, std::as_const(mRoads), fit ? &std::as_const(mRecords) : nullptr, fitSeconds);
    }
    return StageStatus::kComplete;
  }

  /// Hold \a rows exclusively, as a selection takes them: the next
  /// iteration's tracklet stage skips them.
  void holdRows(std::span<const std::uint32_t> rows)
  {
    for (const auto row : rows) {
      mUsed[row] = kExclusive;
    }
  }

  /// The used mask \a used, a checkpoint's: the rows as a run's first
  /// iterations left them, after `beginTimeframe`.
  void restoreUsed(std::span<const std::uint8_t> used)
  {
    if (used.size() != mUsed.size()) {
      throw std::invalid_argument("a used mask of " + std::to_string(used.size()) + " rows for a frame of " + std::to_string(mUsed.size()));
    }
    std::ranges::copy(used, mUsed.begin());
  }

  /// The plan of a timeframe for this tracker's iterations: every slice of
  /// \a sliceWidthBlocks blocks with the halo the frame's reach needs
  /// (`SliceReach::haloPS`), from the frame's \a bounds. \a vertices
  /// are the timeframe's vertex list, sorted by time, for the iterations that
  /// use it.
  [[nodiscard]] TimeframePlan<Det> planTimeframe(const Frame& frame, const HaloBounds<Det>& bounds, std::span<const VertexSeed> vertices, std::int32_t sliceWidthBlocks) const { return o2::gtrk::planTimeframe<Det>(frame, bounds, mIterations, mWindowSettings, mBeamLine, vertices, sliceWidthBlocks); }

  /// The whole timeframe as \a plan has it: for every iteration, every slice,
  /// then the selection, the refit and the follow. \a plan must be this
  /// tracker's plan of \a frame -- std::invalid_argument otherwise -- and
  /// outlive the run.
  [[nodiscard]] TimeframeResult run(const Frame& frame, const TimeframePlan<Det>& plan) { return runIterations(frame, plan, nullptr, mIterations.size()); }

  /// The first \a iterations iterations of `run(frame, plan)`, the same as
  /// far as they go: what `checkpoint` keeps.
  [[nodiscard]] TimeframeResult run(const Frame& frame, const TimeframePlan<Det>& plan, std::size_t iterations)
  {
    if (iterations > mIterations.size()) {
      throw std::invalid_argument(std::to_string(iterations) + " iterations of a tracker of " + std::to_string(mIterations.size()));
    }
    return runIterations(frame, plan, nullptr, iterations);
  }

  /// What \a result, this tracker's last `run` and complete, left: its
  /// tracks, the used mask and its iterations' summaries. Before
  /// `endTimeframe`, which releases them.
  [[nodiscard]] Checkpoint checkpoint(const TimeframeResult& result) const
  {
    if (result.status != StageStatus::kComplete || result.store == nullptr) {
      throw std::invalid_argument("a checkpoint of an incomplete run");
    }
    Checkpoint state;
    const auto& tracks = mRefit.tracks();
    state.tracks.reserve(tracks.size());
    for (std::size_t track = 0; track < tracks.size(); ++track) {
      state.tracks.push_back(tracks.read(track));
    }
    state.used.assign(mUsed.begin(), mUsed.end());
    state.summaries = result.iterations;
    for (auto& summary : state.summaries) {
      summary.profile = {};
    }
    return state;
  }

  /// `run(frame, plan)` from \a from: its tracks, used mask and summaries as
  /// its iterations left them, then this tracker's iterations after those.
  /// The whole run's result where those iterations were this tracker's first
  /// ones on the same frame and plan: the caller checks the settings.
  [[nodiscard]] TimeframeResult resume(const Frame& frame, const TimeframePlan<Det>& plan, const Checkpoint& from)
  {
    const std::size_t first = from.summaries.size();
    if (first > mIterations.size() || from.used.size() != frame.size()) {
      throw std::invalid_argument("a checkpoint of " + std::to_string(first) + " iterations on " + std::to_string(from.used.size()) + " rows for a tracker of " + std::to_string(mIterations.size()) + " on " + std::to_string(frame.size()));
    }
    return runIterations(frame, plan, &from, mIterations.size());
  }

  /// Profile every `run`, or leave the profiles empty (the default).
  void setProfiling(bool enabled) { mProfiling = enabled; }

  /// The tracker's arena, its stages' and loops' threads: for a pass that
  /// shares them (VertexingPass).
  [[nodiscard]] detail::TaskArena& arena() noexcept { return mArena; }

  /// Record every zone and parallel chunk of the next runs into \a trace, or
  /// nothing with nullptr, the default. The trace must outlive the runs.
  void setTrace(profile::Trace* trace) { mTrace = trace; }

  /// One slice of iteration \a iteration, tracklets to levels, reading the used
  /// mask as `beginTimeframe` and the selections since left it: \a slicePlan,
  /// one of \a plan's slices (or planned alike), with \a plan's frame-wide
  /// bounds. \a vertexWindows are \a plan's for the iteration (`planTimeframe`);
  /// a beam-line iteration does not read them.
  [[nodiscard]] Result runSlice(std::size_t iteration, const Frame& frame, const TimeframePlan<Det>& plan, const SlicePlan<Det>& slicePlan, std::span<const VertexSeed> vertices, const VertexWindows& vertexWindows)
  {
    assert(iteration < mIterations.size());
    return runPlannedSlice(iteration, frame, plan, slicePlan, slicePlan.slice, mBeamLine, vertices, vertexWindows, false);
  }

  /// A completed slice of `run`, as an observer sees it: after its roads and
  /// before the selection, with the used mask it read still in place. Every
  /// pointer is valid while the observer runs.
  struct SliceObservation {
    std::size_t iteration = 0U;
    const Frame* frame = nullptr;
    const typename Frame::TimeSlice* slice = nullptr;
    const Result* result = nullptr;
    const RoadCount* roads = nullptr; ///< the slice's road counts, starts and offsets
    const Records* records = nullptr; ///< the slice's selection records
  };
  using SliceObserver = std::function<void(const Tracker&, const SliceObservation&)>;

  /// A completed iteration of `run`: its selection, refit and follow
  /// done. The selection's records() are every slice's, appended in order.
  struct IterationObservation {
    std::size_t iteration = 0U;
    const GlobalSelection<kMaxClusters>* selection = nullptr;
    const IterationSummary* summary = nullptr;
    std::span<const std::uint8_t> used;                        ///< the mask after the iteration
    const typename TrackRefit<Det>::Tracks* tracks = nullptr;  ///< every iteration's so far
    std::uint32_t trackOffset = 0U;                            ///< this iteration's first in `tracks`
    std::span<const typename Follow<Det>::Proposal> proposals; ///< one per track of it; empty without the follow
    /// With followEarlier, one per earlier track: their follow before this
    /// iteration's seeds.
    std::span<const typename Follow<Det>::Proposal> earlierProposals;
    /// With the follow, each of the iteration's tracks as the follow started
    /// from it: its inner and outward states and its rows then (the first
    /// `followStartCount` of its clusters), for a replay;
    /// empty without the follow.
    std::span<const BoundTrackState> followStartInner;
    std::span<const BoundTrackState> followStart;
    std::span<const std::uint8_t> followStartCount;
    /// With the follow, the used mask as the follow read it, before its
    /// rows and the length cut; empty without.
    std::span<const std::uint8_t> followUsed;
    /// With the follow, each selected record's track as the follow and its
    /// refit left it, before the length cut and the pT floor after the follow
    /// drop any (one per selected record, in its order): its rows, their
    /// count and its state; empty without the follow.
    std::span<const ClusterList<kMaxClusters>> followEndClusters;
    std::span<const std::uint8_t> followEndCount;
    std::span<const BoundTrackState> followEndState;
  };
  using IterationObserver = std::function<void(const Tracker&, const IterationObservation&)>;

  /// Called by `run` after every completed slice -- QA that repeats a slice's
  /// work elsewhere, such as on a device. Nothing it does changes the tracking.
  void setSliceObserver(SliceObserver observer) { mSliceObserver = std::move(observer); }
  /// Called by `run` after every completed iteration, as the slice observer.
  void setIterationObserver(IterationObserver observer) { mIterationObserver = std::move(observer); }

  /// The source map of the last slice: its halo rows.
  [[nodiscard]] const typename Scratch::Sources& sources() const { return mScratch.sources(); }
  /// The beam line as a vertex seed, which a beam-line iteration reads.
  [[nodiscard]] const VertexSeed& beamLine() const { return mBeamLine; }
  /// The most tracklets, triplets and neighbours one slice placed.
  [[nodiscard]] std::size_t trackletHighWater() const { return mTrackletHighWater; }
  [[nodiscard]] std::size_t tripletHighWater() const { return mTripletHighWater; }
  [[nodiscard]] std::size_t neighbourHighWater() const { return mNeighbourHighWater; }
  [[nodiscard]] std::size_t stagingHighWater() const { return mStagingHighWater; }
  [[nodiscard]] std::size_t stagingPredictionMisses() const { return mStagingPredictionMisses; }
  [[nodiscard]] std::size_t memoryHighWater() const { return mMemoryPool->getPeakMemory(); }

 private:
  [[nodiscard]] const IterationParameters<Det>& parameters() const { return mIterations[mIteration]; }

  /// \a site's staging for the timeframe, in \a backend: what the estimator
  /// learned from the timeframes before.
  template <typename StageBackend>
  void stageStaging(StagingSite site, StageBackend& backend)
  {
    auto& granted = mGranted[static_cast<std::size_t>(site)];
    granted = mEstimator.perThread(static_cast<std::size_t>(site), mRows);
    backend.prepareStagingPerThread(granted);
    backend.resetLargestUse();
  }

  /// Teach the estimator what \a backend's buffers held over the timeframe.
  template <typename StageBackend>
  void learnStaging(StagingSite site, const StageBackend& backend)
  {
    mEstimator.learn(static_cast<std::size_t>(site), backend.largestUse(), mGranted[static_cast<std::size_t>(site)], mRows, backend.grewSinceReset());
    mStagingPredictionMisses += backend.grewSinceReset() ? 1U : 0U;
  }

  /// `run` from \a from (none: the start) to iteration \a last.
  [[nodiscard]] TimeframeResult runIterations(const Frame& frame, const TimeframePlan<Det>& plan, const Checkpoint* from, std::size_t last)
  {
    plan.check(frame, mIterations);
    TimeframeResult result;
    profiled(result.profile, Phase::kRun, [&] { runTimeframe(frame, plan, result, from, last); });
    return result;
  }

  /// `runIterations` without its profile's activation: the iterations after
  /// \a from's (from the first without one) to \a last into \a result.
  void runTimeframe(const Frame& frame, const TimeframePlan<Det>& plan, TimeframeResult& result, const Checkpoint* from, std::size_t last)
  {
    {
      const profile::Zone zone{Phase::kBegin};
      result.status = beginTimeframe(frame);
    }
    if (result.status != StageStatus::kComplete) {
      return;
    }
    mRefit.clear();
    mStore.clear();
    std::size_t first = 0U;
    if (from != nullptr) {
      first = from->summaries.size();
      std::ranges::copy(from->used, mUsed.begin());
      // The tracks into the pool: a limit ends the run with a status, as anywhere.
      result.status = guardOutOfMemory([&] {
        auto& tracks = mRefit.tracks();
        tracks.resize(from->tracks.size());
        for (std::size_t track = 0; track < from->tracks.size(); ++track) {
          tracks.write(track, from->tracks[track]);
        }
        return StageStatus::kComplete;
      });
      if (result.status == StageStatus::kComplete) {
        result.status = commitTracks();
      }
      if (result.status != StageStatus::kComplete) {
        return;
      }
      result.iterations = from->summaries;
    }
    for (std::size_t iteration = first; iteration < last; ++iteration) {
      IterationSummary summary;
      profiled(summary.profile, Phase::kIteration, [&] { runIteration(iteration, frame, plan, summary); });
      result.profile += summary.profile;
      result.iterations.push_back(summary);
      if (summary.status != StageStatus::kComplete) {
        result.status = summary.status;
        return;
      }
    }
    result.store = &mStore;
  }

  /// The tracks not committed yet into the store (TrackStore::commit), with
  /// \a partnered the overlap partners' last run gave them.
  [[nodiscard]] StageStatus commitTracks(bool partnered = false)
  {
    const auto count = static_cast<Index>(mRefit.size() - mStore.size());
    if (partnered) {
      return mStore.commit(spreadForEach(), spreadChunk(count), std::as_const(mRefit).tracks(), mPartners.offsets(), mPartners.partners());
    }
    return mStore.commit(spreadForEach(), spreadChunk(count), std::as_const(mRefit).tracks());
  }

  /// What one iteration grew, back to the pool at its end, after its
  /// observers: the slices' storage -- the
  /// scratch, the stages' staging and chunk lists, the road counts, the
  /// records -- and the selection's records' and the follows'. Kept from one iteration
  /// to the next, each buffer would hold the largest size any iteration gave
  /// it, whichever shape the next one has. The staging sizes the estimator
  /// granted stay; the next iteration allocates its own.
  void releaseIterationStorage() noexcept
  {
    mScratch.clearStorage();
    mTrackletStage.releaseStorage();
    mTripletStage.releaseStorage();
    mNeighbourStage.releaseStorage();
    mRecordStage.releaseStorage();
    mRoads.clearStorage();
    mRecords.clearStorage();
    mSelection.clearRecordStorage();
    mFollow.clearStorage();
    mFollowEarlier.clearStorage();
    mJoin.clearStorage();
    mInactiveRows.clearStorage();
    mPartners.clearStorage();
    deepVectorClear(mMemoryPool.get(), mKeptTracks, mChangedTracks);
    mDropScratch.setMemoryResource(mMemoryPool.get());
  }

  /// Iteration \a iteration of `run`: its slices, the selection, the refit
  /// and the follow, into \a summary, whose status says how far it got.
  void runIteration(std::size_t iteration, const Frame& frame, const TimeframePlan<Det>& plan, IterationSummary& summary)
  {
    const auto forEach = hostForEach();
    const auto fits = fitForEach();
    // The selection's sorts on the tracker's threads, as its loops.
    const auto sort = [&](std::span<std::uint32_t> indices, const auto& less) { mNeighbourStage.sortIndices(indices, less); };
    summary.haloPS = plan.haloPS;
    mSelection.clear();
    // With followEarlier the earlier iterations' tracks are followed first,
    // with this iteration's settings: the rows they take are theirs before
    // this iteration's seeds are made.
    const bool earlier = mIterations[iteration].followEarlier && mRefit.size() != 0U;
    if (earlier) {
      const profile::Zone zone{Phase::kFollowEarlier};
      summary.status = runFollow(mFollowEarlier, mIterations[iteration], frame, plan, 0U, true, summary.tracksFollowedEarlier, summary.clustersFollowedEarlierOutward, summary.clustersFollowedEarlierInward);
      if (summary.status == StageStatus::kComplete) {
        // The tracks it extended or trimmed in place, committed again.
        summary.status = mStore.recommitChanged(spreadForEach(), spreadChunk(static_cast<Index>(mStore.size())), std::as_const(mRefit).tracks(), mChangedTracks);
      }
      if (summary.status != StageStatus::kComplete) {
        return;
      }
    }
    for (const auto& slicePlan : plan.slices) {
      if (summary.status != StageStatus::kComplete) {
        break;
      }
      const auto sliced = runPlannedSlice(iteration, frame, plan, slicePlan, slicePlan.slice, plan.beamLine, plan.vertices, plan.vertexWindows[iteration], true);
      summary.status = sliced.status;
      if (summary.status == StageStatus::kComplete) {
        {
          const profile::Zone zone{Phase::kRoads};
          // Guarded whole: the fit's setup allocates before its stage's own
          // guard.
          summary.status = guardOutOfMemory([&] { return runRoads(frame, summary, iteration); });
        }
        if (mSliceObserver && summary.status == StageStatus::kComplete) {
          const profile::Zone zone{Phase::kObservers};
          mSliceObserver(*this, {
                                  .iteration = iteration,
                                  .frame = &frame,
                                  .slice = &slicePlan.slice,
                                  .result = &sliced,
                                  .roads = &mRoads,
                                  .records = &mRecords,
                                });
        }
      }
      ++summary.slices;
    }
    if (summary.status != StageStatus::kComplete) {
      return;
    }

    {
      const profile::Zone zone{Phase::kSelection};
      // The mask as the selection finds it, for the drop after the follow:
      // any selected track's refit may fail.
      summary.status = guardOutOfMemory([&] {
        // Overwritten whole: kept from iteration to iteration, copied in parallel.
        resizeBoundedVector(mUsedBeforeSelection, mUsed.size(), mMemoryPool.get());
        copyInParallel(std::span<const std::uint8_t>{mUsed}, std::span<std::uint8_t>{mUsedBeforeSelection});
        return StageStatus::kComplete;
      });
      if (summary.status != StageStatus::kComplete) {
        return;
      }
      summary.status = mSelection.select(forEach, sort, mIterations[iteration].sharedInnermost, mUsed);
    }
    summary.records = mSelection.recordCount();
    summary.recordBytes = summary.records * sizeof(typename BoundedColumns<SelectionRecordOf<kMaxClusters>::template RecordT>::Record);
    summary.duplicates = mSelection.duplicates();
    summary.mismatches = mSelection.mismatches();
    summary.rounds = mSelection.rounds();
    summary.tracks = mSelection.selected().size();
    const auto trackOffset = static_cast<std::uint32_t>(mRefit.size());
    bool followed = false;
    bool started = false;
    if (summary.status == StageStatus::kComplete) {
      {
        const profile::Zone zone{Phase::kRefit};
        summary.status = mRefit.run(fits, mSelection, mIterations[iteration], frame, *mGeometry, mPropagator, static_cast<std::uint8_t>(iteration));
      }
      const auto status = mRefit.tracks().view().status;
      for (auto track = static_cast<std::size_t>(trackOffset); track < mRefit.size(); ++track) {
        summary.refitFailures += static_cast<std::size_t>(status[track] != PropagationStatus::kOk);
      }
      if (summary.status == StageStatus::kComplete && mIterations[iteration].follow) {
        // What the follow starts from, kept only for an observer's replay.
        started = static_cast<bool>(mIterationObserver);
        if (started) {
          const auto view = std::as_const(mRefit.tracks()).view();
          mFollowStartInner.assign(view.state + trackOffset, view.state + mRefit.size());
          mFollowStart.assign(view.outwardState + trackOffset, view.outwardState + mRefit.size());
          mFollowStartCount.assign(view.count + trackOffset, view.count + mRefit.size());
          mFollowUsed.assign(mUsed.begin(), mUsed.end());
        }
        {
          const profile::Zone zone{Phase::kFollow};
          summary.status = runFollow(mFollow, mIterations[iteration], frame, plan, trackOffset, false, summary.tracksFollowed, summary.clustersFollowedOutward, summary.clustersFollowedInward);
        }
        followed = summary.status == StageStatus::kComplete;
        if (started && followed) {
          const auto view = std::as_const(mRefit.tracks()).view();
          mFollowEndClusters.assign(view.clusters + trackOffset, view.clusters + mRefit.size());
          mFollowEndCount.assign(view.count + trackOffset, view.count + mRefit.size());
          mFollowEndState.assign(view.state + trackOffset, view.state + mRefit.size());
        }
      }
    }
    if ((summary.refitFailures != 0U || dropsAfterFollow(mIterations[iteration])) && summary.status == StageStatus::kComplete) {
      const profile::Zone zone{Phase::kLengthCut};
      summary.status = guardOutOfMemory([&] {
        summary.tracksDropped = dropTracks(mIterations[iteration], trackOffset, followed, started) - summary.refitFailures;
        return StageStatus::kComplete;
      });
    }
    // The join of the iteration's tracks to the earlier iterations' ones;
    // none where an iteration attaches overlap partners, which a joined
    // earlier track's references would lose.
    if (summary.status == StageStatus::kComplete && trackOffset != 0U && std::ranges::none_of(mIterations, [](const auto& parameters) { return parameters.overlapPartners; })) {
      const profile::Zone zone{Phase::kJoin};
      summary.status = joinTracks(frame, plan, trackOffset, followed, started, summary.tracksJoined);
    }
    if constexpr (Detector<Det>::inactiveLayers != 0U) {
      const auto joined = summary.tracksJoined != 0U ? mJoin.changed() : std::span<const std::uint32_t>{};
      if (summary.status == StageStatus::kComplete) {
        summary.status = gatherInactiveCandidates(trackOffset, joined);
      }
      if (summary.status == StageStatus::kComplete) {
        const profile::Zone zone{Phase::kInactiveRows};
        const auto forEach = [&](std::uint32_t count, const auto& body) { forEachIn(spreadChunk(static_cast<Index>(count)))(static_cast<Index>(count), body); };
        summary.status = mInactiveRows.run(forEach, mIterations[iteration], frame, *mGeometry, mPropagator, mRefit.tracks(), std::span<const std::uint32_t>{mInactiveCandidates}, std::span<std::uint8_t>{mUsed}, plan.layerMaxWidthPS, mSensorCells.view());
        summary.inactiveRows = mInactiveRows.changed().size();
      }
    }
    const auto changedEarlier = [&](std::span<const std::uint32_t> changed) { return !changed.empty() && changed.front() < trackOffset; };
    if (summary.status == StageStatus::kComplete && (summary.tracksJoined != 0U || changedEarlier(mInactiveRows.changed()))) {
      // The earlier tracks the join or the inactive rows changed, committed again.
      summary.status = mStore.recommitChanged(spreadForEach(), spreadChunk(static_cast<Index>(mStore.size())), std::as_const(mRefit).tracks(), mChangedTracks);
    }
    const bool partnered = summary.status == StageStatus::kComplete && mIterations[iteration].overlapPartners;
    if (partnered) {
      const profile::Zone zone{Phase::kOverlapPartners};
      summary.status = mPartners.run(fitForEach(), fitChunk(static_cast<Index>(mRefit.size() - trackOffset)), mIterations[iteration], frame, *mGeometry, mPropagator, mSensorCells.view(), std::span<const SensorFootprint>{mSensorCells.footprints()}, mRefit.tracks(), trackOffset, std::span<std::uint8_t>{mUsed}, std::span<const std::uint32_t>{plan.layerMaxWidthPS}, std::span<const Scalar>{plan.layerMaxLocVariance});
      summary.overlapPartners = mPartners.counts().attached;
      summary.tracksWithPartners = mPartners.counts().tracks;
    }
    if (summary.status == StageStatus::kComplete) {
      // The iteration's tracks as the follow, the length cut and the
      // partners left them.
      const profile::Zone zone{Phase::kCommit};
      summary.status = commitTracks(partnered);
    }
    if (mIterationObserver && summary.status == StageStatus::kComplete) {
      const profile::Zone zone{Phase::kObservers};
      mIterationObserver(*this, {
                                  .iteration = iteration,
                                  .selection = &mSelection,
                                  .summary = &summary,
                                  .used = mUsed,
                                  .tracks = &mRefit.tracks(),
                                  .trackOffset = trackOffset,
                                  .proposals = followed ? mFollow.proposals() : std::span<const typename Follow<Det>::Proposal>{},
                                  .earlierProposals = earlier ? mFollowEarlier.proposals() : std::span<const typename Follow<Det>::Proposal>{},
                                  .followStartInner = started ? std::span<const BoundTrackState>{mFollowStartInner} : std::span<const BoundTrackState>{},
                                  .followStart = started ? std::span<const BoundTrackState>{mFollowStart} : std::span<const BoundTrackState>{},
                                  .followStartCount = started ? std::span<const std::uint8_t>{mFollowStartCount} : std::span<const std::uint8_t>{},
                                  .followUsed = started ? std::span<const std::uint8_t>{mFollowUsed} : std::span<const std::uint8_t>{},
                                  .followEndClusters = started && followed ? std::span<const ClusterList<kMaxClusters>>{mFollowEndClusters} : std::span<const ClusterList<kMaxClusters>>{},
                                  .followEndCount = started && followed ? std::span<const std::uint8_t>{mFollowEndCount} : std::span<const std::uint8_t>{},
                                  .followEndState = started && followed ? std::span<const BoundTrackState>{mFollowEndState} : std::span<const BoundTrackState>{},
                                });
    }
    releaseIterationStorage();
  }

  /// \a follow run on the tracks from \a trackOffset with \a parameters
  /// (`Follow::run`), \a earlier for the earlier iterations' tracks; its
  /// proposals tallied into \a tracks, \a outward and \a inward.
  [[nodiscard]] StageStatus runFollow(Follow<Det>& follow, const IterationParameters<Det>& parameters, const Frame& frame, const TimeframePlan<Det>& plan, std::uint32_t trackOffset, bool earlier, std::size_t& tracks, std::size_t& outward, std::size_t& inward)
  {
    const auto status = follow.run(fitForEach(), parameters, frame, *mGeometry, mPropagator, mRefit.tracks(), trackOffset, parameters.sharedInnermost == 0U, earlier, mUsed, plan.layerMaxWidthPS, mSensorCells.view());
    if (status != StageStatus::kComplete) {
      return status;
    }
    for (const auto& proposal : follow.proposals()) {
      tracks += static_cast<std::size_t>(proposal.outwardCount + proposal.inwardCount != 0U);
      outward += proposal.outwardCount;
      inward += proposal.inwardCount;
    }
    return status;
  }

  /// \a body() with its profile \a into open and \a phase's zone, where
  /// profiling is on; plain otherwise.
  template <typename Body>
  void profiled(profile::Profile& into, Phase phase, const Body& body)
  {
    if (!mProfiling) {
      body();
      return;
    }
    const profile::Activation activation{into, mMemoryPool.get(), mTrace};
    const profile::Zone zone{phase};
    body();
  }

  /// The tracks from \a trackOffset on whose refit failed, and those that
  /// \a parameters do not keep after the follow (`keptAfterFollow`: fewer
  /// clusters than its length after it, or below its pT floor), dropped: the others keep their order, and their proposals with them when
  /// \a followed, and the follow's starting states when \a started kept
  /// them for an observer. The used mask becomes the one before the
  /// selection plus the rows of every kept track, marked as they are now: a dropped track's rows
  /// no kept one holds are free again, an earlier iteration's are untouched.
  /// Returns how many were dropped.
  std::size_t dropTracks(const IterationParameters<Det>& parameters, std::uint32_t trackOffset, bool followed, bool started)
  {
    auto& tracks = mRefit.tracks();
    const auto view = std::as_const(tracks).view();
    const auto count = static_cast<Index>(tracks.size() - trackOffset);
    const auto chunk = spreadChunk(count);
    const auto forEach = forEachIn(chunk);
    // The kept tracks, as offsets from trackOffset, in order (the serial
    // loop's result, in parallel).
    host_algorithms::compactIndices(forEach, count, chunk, [&](Index track) { return view.status[trackOffset + track] == PropagationStatus::kOk && keptAfterFollow(parameters, view.count[trackOffset + track], view.state[trackOffset + track]); }, mKeptTracks);
    const auto kept = static_cast<Index>(mKeptTracks.size());
    // Their rows as they are now; two kept tracks sharing a row write the same value.
    forEachIn(spreadChunk(kept))(kept, [&](Index at) {
      const auto track = trackOffset + mKeptTracks[at];
      for (std::uint8_t i = 0U; i < view.count[track]; ++i) {
        const auto row = view.clusters[track].index[i];
        std::atomic_ref<std::uint8_t>{mUsedBeforeSelection[row]}.store(mUsed[row], std::memory_order_relaxed);
      }
    });
    const std::size_t dropped = count - kept;
    keepTracks(trackOffset, followed, started);
    copyInParallel(std::span<const std::uint8_t>{mUsedBeforeSelection}, std::span<std::uint8_t>{mUsed});
    return dropped;
  }

  /// The join of the tracks from \a trackOffset on to the earlier ones
  /// (TrackJoin.h): the joined ones removed (`keepTracks`), the earlier ones
  /// that took them committed again; how many into \a joined.
  [[nodiscard]] StageStatus joinTracks(const Frame& frame, const TimeframePlan<Det>& plan, std::uint32_t trackOffset, bool followed, bool started, std::size_t& joined)
  {
    auto& tracks = mRefit.tracks();
    const auto forEach = [&](std::uint32_t count, const auto& body) { forEachIn(spreadChunk(static_cast<Index>(count)))(static_cast<Index>(count), body); };
    auto status = mJoin.run(forEach, mIterations.data(), frame, *mGeometry, mPropagator, tracks, trackOffset, plan.layerMaxWidthPS);
    joined = mJoin.joined().size();
    if (status != StageStatus::kComplete || joined == 0U) {
      return status;
    }
    status = guardOutOfMemory([&] {
      const auto count = static_cast<Index>(tracks.size() - trackOffset);
      resizeBoundedVector(mKeptTracks, count - joined, mMemoryPool.get());
      const auto gone = mJoin.joined();
      for (Index at = 0, next = 0, kept = 0; at < count; ++at) {
        if (next < gone.size() && gone[next] == at) {
          ++next;
          continue;
        }
        mKeptTracks[kept++] = at;
      }
      keepTracks(trackOffset, followed, started);
      return StageStatus::kComplete;
    });
    return status;
  }

  /// The candidates of the inactive layers' rows: the earlier tracks the
  /// join changed (\a joined, ascending), then the tracks from
  /// \a trackOffset on.
  [[nodiscard]] StageStatus gatherInactiveCandidates(std::uint32_t trackOffset, std::span<const std::uint32_t> joined)
  {
    return guardOutOfMemory([&] {
      const auto own = static_cast<std::uint32_t>(mRefit.size()) - trackOffset;
      resizeBoundedVector(mInactiveCandidates, joined.size() + own, mMemoryPool.get());
      std::ranges::copy(joined, mInactiveCandidates.begin());
      for (std::uint32_t at = 0; at < own; ++at) {
        mInactiveCandidates[joined.size() + at] = trackOffset + at;
      }
      return StageStatus::kComplete;
    });
  }

  /// The tracks from \a trackOffset on, of them only those `mKeptTracks`
  /// lists (offsets from \a trackOffset, ascending): the others removed, the
  /// kept ones in their order, their proposals with them when \a followed,
  /// and the follow's starting states when \a started kept them for an
  /// observer.
  void keepTracks(std::uint32_t trackOffset, bool followed, bool started)
  {
    auto& tracks = mRefit.tracks();
    const auto kept = static_cast<Index>(mKeptTracks.size());
    // Up to the first dropped track every kept one stays where it is; from it
    // on each moves down, through a scratch: in place, a parallel move would
    // overwrite tracks still to be read.
    const auto first = o2::gpu::CAAlgo::partitionPoint(Index{0U}, kept, [&](Index at) { return mKeptTracks[at] == at; });
    const auto moved = kept - first;
    if (moved != 0U) {
      const auto forEachMoved = forEachIn(spreadChunk(moved));
      mDropScratch.resizeRewritten(moved);
      forEachMoved(moved, [&](Index at) { mDropScratch.write(at, tracks.read(trackOffset + mKeptTracks[first + at])); });
      forEachMoved(moved, [&](Index at) { tracks.write(trackOffset + first + at, mDropScratch.read(at)); });
      if (followed) {
        mFollow.keepProposals(std::span<const Index>{mKeptTracks}.subspan(first), first, forEachMoved);
      }
      if (started) {
        for (Index at = first; at < kept; ++at) {
          mFollowStartInner[at] = mFollowStartInner[mKeptTracks[at]];
          mFollowStart[at] = mFollowStart[mKeptTracks[at]];
          mFollowStartCount[at] = mFollowStartCount[mKeptTracks[at]];
        }
      }
    }
    tracks.resize(trackOffset + kept);
    if (followed) {
      mFollow.truncateProposals(kept);
    }
    if (started) {
      mFollowStartInner.resize(kept);
      mFollowStart.resize(kept);
      mFollowStartCount.resize(kept);
    }
  }

  /// \a from copied into \a to (as long), on the host threads.
  template <typename T>
  void copyInParallel(std::span<const T> from, std::span<T> to)
  {
    const auto count = static_cast<Index>(from.size());
    forEachIn(spreadChunk(count))(count, [&](Index at) { to[at] = from[at]; });
  }

  /// The tracklet stage over host memory of \a inputs, as `runPlannedSlice`
  /// runs it, reading \a beamLine for a beam-line iteration (the general
  /// one, or the beam-line one for a beam-line iteration, which makes the
  /// same tracklets).
  template <bool BeamLine, bool Shape, bool DiscWindow>
  [[nodiscard]] TrackletStageFor<Det, BeamLine, Shape, DiscWindow> trackletStage(const TrackletInputs<Det>& inputs, const Frame& frame, const VertexSeed& beamLine, std::span<const VertexSeed> vertices, const VertexWindows& vertexWindows) const
  {
    return {*inputs.parameters, *inputs.sources, frame, inputs.targetBlocks, beamLine, vertices, inputs.maxCandidateWidthPS, vertexWindows, mUsed.data(), inputs.gate};
  }

  struct LayeredNeighbours {
    NeighbourStage<Det> body;
    std::array<Index, kLayerCount + 1U> layerBoundaries;
    LayerDepthGroups<Det> layerGroups = ::o2::gtrk::layerDepthGroups<Det>;
    template <typename Emitter>
    void operator()(Index element, Emitter& emitter) const
    {
      body(element, emitter);
    }
  };

  /// One planned slice of iteration \a iteration, tracklets to levels, with
  /// \a plan's frame-wide bounds; the stage reads \a beamLine or \a
  /// vertices, as the iteration chose, and its gate the core ranges of \a
  /// owner, the slice the caller keeps.
  [[nodiscard]] Result runPlannedSlice(std::size_t iteration, const Frame& frame, const TimeframePlan<Det>& plan, const SlicePlan<Det>& slicePlan, const typename Frame::TimeSlice& owner, const VertexSeed& beamLine, std::span<const VertexSeed> vertices, const VertexWindows& vertexWindows, bool countRoads)
  {
    assert(iteration < mIterations.size() && mUsed.size() == frame.size());
    const profile::Zone zone{Phase::kSlices};
    mIteration = iteration;
    const StageStatus status = guardOutOfMemory([&] {
      auto result = StageStatus::kComplete;
      {
        const profile::Zone trackletZone{Phase::kTracklets};
        std::optional<profile::Zone> setup{std::in_place, Phase::kTrackletSetup};
        beginTrackletSlice(plan, slicePlan, owner);
        const auto& parameters = *mTrackletInputs.parameters;
        result = withTrackletVertices(parameters.vertexSource, [&](auto beamLineOnly) {
          return withShapeTest(parameters.shapeTest, [&](auto shape) {
            return withDiscWindow<decltype(beamLineOnly)::value>(parameters.discWindowMinP > 0.F, [&](auto discWindow) {
              const auto stage = trackletStage<decltype(beamLineOnly)::value, decltype(shape)::value, decltype(discWindow)::value>(mTrackletInputs, frame, beamLine, vertices, vertexWindows);
              setup.reset();
              return mTrackletStage.run(mScratch.tracklets(), mChunkSize, stage);
            });
          });
        });
        if (result == StageStatus::kComplete) {
          // Each source's tracklets in (link, lambda, target) order, which the
          // triplet stage searches.
          {
            const profile::Zone sortZone{Phase::kTrackletSort};
            auto& tracklets = mScratch.tracklets();
            const auto view = tracklets.view();
            const auto offsets = tracklets.offsets();
            mTrackletStage.forEachElement(tracklets.elementCount(), mChunkSize, [&](Index element) { sortTrackletSegmentOnHost(view, offsets[element], offsets[element + 1U]); });
          }
        }
      }
      if (result != StageStatus::kComplete) {
        return result;
      }
      const auto triplets = runTriplets(frame, plan);
      if (triplets != StageStatus::kComplete) {
        return triplets;
      }
      return runNeighbours(frame, owner, countRoads);
    });
    if (status == StageStatus::kComplete) {
      mTrackletHighWater = std::max<std::size_t>(mTrackletHighWater, mScratch.tracklets().size());
      mStagingHighWater = std::max(mStagingHighWater, mTrackletStage.stagingCapacity());
    }
    const bool complete = status == StageStatus::kComplete;
    return {
      .status = status,
      .tracklets = complete ? &mScratch.tracklets() : nullptr,
      .triplets = complete ? &mScratch.triplets() : nullptr,
      .neighbours = complete ? &mScratch.neighbours() : nullptr,
    };
  }

  /// The triplet stage, its staged triplets kept only where they are on a
  /// triplet path of the iteration's minClusters: the bounds built from the staged triplets, then each chunk's
  /// compacted before the place. The others make no neighbour, as element or
  /// as predecessor (NeighbourStage.h), and in a low-pT pass they are most of
  /// the triplets.
  template <typename Stage>
  struct OnPathTriplets {
    Tracker* tracker = nullptr;
    const Stage* stage = nullptr;

    template <typename Emitter>
    void operator()(Index element, Emitter& emitter) const
    {
      (*stage)(element, emitter);
    }

    template <typename Staged>
    void filterStaged(const Staged& staged) const
    {
      const profile::Zone zone{Phase::kTripletPaths};
      auto& scratch = tracker->mScratch;
      auto& paths = scratch.tripletPaths();
      paths.build(tracker->parameters(), scratch.sources(), scratch.tracklets(), scratch.permutation(), [&](Index begin, Index end, const auto& body) { staged.forEach(begin, end, [&](Index element, const auto& triplets, std::size_t at) { body(element, triplets.first[at], triplets.second[at]); }); }, tracker->hostForEach());
      staged.compact([&](Index /*element*/, const auto& triplets, std::size_t at) { return paths.onPath(triplets.first[at], triplets.second[at]); });
    }
  };

  /// The fit of the counted roads with iteration \a owner's settings into
  /// \a records, whose staging is that iteration's. An iteration with the
  /// beam constraint fits on the stage that also tests the beam line; the
  /// host alone makes it.
  [[nodiscard]] StageStatus fitRoads(const Frame& frame, std::size_t owner, Records& records)
  {
    const auto roads = mRoads.roadCount();
    mIteration = owner;
    const auto& fitted = mIterations[owner];
    std::optional<profile::Zone> setup{std::in_place, Phase::kRoadFitSetup};
    records.prepare(roads);
    const auto run = [&]<bool BeamLine, bool LastStepAllowance>() {
      const RoadFitStage<Det, BeamLine, LastStepAllowance> stage{fitted, frame, *mGeometry, mPropagator, mScratch.tracklets(), mScratch.triplets(), mRoads.walk(mScratch.triplets(), mScratch.neighbours())};
      setup.reset();
      return mRecordStage.run(records, fitChunk(static_cast<Index>(roads)), stage);
    };
    // An iteration with a kink allowance takes it on the fit's last step
    // (against the road-end fakes).
    return withBools([&]<bool BeamLine, bool LastStepAllowance>(std::bool_constant<BeamLine>, std::bool_constant<LastStepAllowance>) { return run.template operator()<BeamLine, LastStepAllowance>(); }, fitted.beamConstraint, takesLastStepAllowance(fitted));
  }

  /// The slice's roads, fitted, into the selection.
  [[nodiscard]] StageStatus runRoads(const Frame& frame, IterationSummary& summary, std::size_t iteration)
  {
    const auto counted = countRoads();
    if (counted != StageStatus::kComplete) {
      return counted;
    }
    summary.roads += mRoads.roadCount();
    const auto fitted = fitRoads(frame, iteration, mRecords);
    if (fitted != StageStatus::kComplete) {
      return fitted;
    }
    const profile::Zone zone{Phase::kSelectionAppend};
    return mSelection.append(mRecords, spreadForEach());
  }

  /// The slice's roads numbered from the counts that came with the
  /// neighbours (`runNeighbours`, countRoads): `RoadCount::runPrepared`.
  [[nodiscard]] StageStatus countRoads()
  {
    return guardOutOfMemory([&] {
      const profile::Zone zone{Phase::kRoadCount};
      return mRoads.runPrepared(hostForEach(), roadStartRule(parameters()), mScratch.tracklets(), mScratch.triplets(), mChunkSize);
    });
  }

  /// The slice's source map in the scratch, and the other inputs of its
  /// tracklet stage, from its plan and \a plan's widths, in mTrackletInputs.
  void beginTrackletSlice(const TimeframePlan<Det>& plan, const SlicePlan<Det>& slicePlan, const typename Frame::TimeSlice& owner)
  {
    mScratch.beginSlice(slicePlan.sources);
    mTrackletInputs = slicePlan.inputs(parameters(), plan.layerMaxWidthPS);
    mTrackletInputs.sources = &mScratch.sources();
    mTrackletInputs.gate.core = owner.core.data();
  }

  /// The triplet stage over the slice's sorted tracklets: the first-tracklet
  /// order, then generate/scan/place with \a plan's target bounds.
  [[nodiscard]] StageStatus runTriplets(const Frame& frame, const TimeframePlan<Det>& plan)
  {
    const profile::Zone zone{Phase::kTriplets};
    auto& permutation = mScratch.permutation();
    const auto& tracklets = mScratch.tracklets();
    auto& trackletPaths = mScratch.trackletPaths();
    {
      const profile::Zone pathsZone{Phase::kTrackletPaths};
      trackletPaths.build(parameters(), mScratch.sources(), tracklets, hostForEach());
    }
    const auto paths = trackletPaths.bounds();
    auto ordered = StageStatus::kComplete;
    {
      const profile::Zone permutationZone{Phase::kPermutation};
      ordered = permutation.build(parameters(), mScratch.sources(), tracklets, paths, hostForEach(), mChunkSize);
    }
    if (ordered != StageStatus::kComplete) {
      return ordered;
    }
    auto& triplets = mScratch.triplets();
    std::optional<profile::Zone> setup{std::in_place, Phase::kTripletSetup};
    triplets.prepare(permutation.size());
    const auto result = withShapeTest(parameters().shapeTest, [&](auto shape) {
      const TripletStage<Det, decltype(shape)::value> stage{parameters(), mScratch.sources(), frame, tracklets, plan.layerTargetBounds.data(), permutation.data(), paths};
      setup.reset();
      return mTripletStage.run(triplets, mChunkSize, OnPathTriplets<TripletStage<Det, decltype(shape)::value>>{this, &stage});
    });
    if (result == StageStatus::kComplete) {
      mTripletHighWater = std::max<std::size_t>(mTripletHighWater, triplets.size());
    }
    return result;
  }

  /// The neighbour stage over every triplet, pulling its predecessors from
  /// the triplets ending with its first tracklet.
  [[nodiscard]] StageStatus runNeighbours(const Frame& frame, const typename Frame::TimeSlice& owner, bool countRoads)
  {
    const profile::Zone zone{Phase::kNeighbours};
    auto& triplets = mScratch.triplets();
    auto& ending = mScratch.ending();
    auto indexed = StageStatus::kComplete;
    {
      const profile::Zone indexZone{Phase::kEndIndex};
      // One task per chunk or bucket of the index's two-level sort.
      indexed = ending.build(mScratch.tracklets().size(), triplets, forEachIn(1U), fitChunkCount());
    }
    if (indexed != StageStatus::kComplete) {
      return indexed;
    }
    const auto& tripletPaths = mScratch.tripletPaths();
    const auto paths = tripletPaths.bounds();
    auto& neighbours = mScratch.neighbours();
    std::optional<profile::Zone> setup{std::in_place, Phase::kNeighbourSetup};
    neighbours.prepare(triplets.size());
    // Layer by layer, each triplet's level final when it is generated, and
    // with \a countRoads the road counts with it (`OrderedRoadCounts`).
    OrderedRoadCounts counts;
    if (countRoads) {
      mRoads.prepareOrdered(hostForEach(), triplets.size());
      // A start-layer mask empties the core of the layers roads may not start on.
      const auto startLayers = parameters().roadStartLayers;
      if (startLayers != 0U) {
        mRoadCores = roadCores(std::span<const MeasurementRange, kLayerCount>{owner.core}, startLayers);
      }
      counts = mRoads.ordered(parameters().linkSources, startLayers != 0U ? mRoadCores.data() : owner.core.data());
    }
    // Each layer's first triplet: its first permuted tracklet's (`TripletPaths::layerElements`).
    std::array<Index, kLayerCount + 1U> boundaries{};
    for (std::size_t layer = 0U; layer <= kLayerCount; ++layer) {
      boundaries[layer] = triplets.offsets()[tripletPaths.layerElements()[layer]];
    }
    const NeighbourStage<Det> stage{parameters(), frame, mScratch.tracklets(), triplets, ending, paths, triplets.view().level, counts};
    setup.reset();
    const auto result = mNeighbourStage.run(neighbours, mChunkSize, LayeredNeighbours{stage, boundaries});
    if (result == StageStatus::kComplete) {
      mNeighbourHighWater = std::max<std::size_t>(mNeighbourHighWater, neighbours.size());
    }
    return result;
  }

  /// body(i) for every i in [0, count) on the host threads, in chunks of the
  /// stages' size: the host algorithms' `forEach` (HostAlgorithms.h).
  [[nodiscard]] auto hostForEach() { return forEachIn(mChunkSize); }

  /// `hostForEach` in chunks of \a chunk elements.
  [[nodiscard]] auto forEachIn(Index chunk)
  {
    return [this, chunk](Index count, const auto& body) { mNeighbourStage.forEachElement(count, chunk, body); };
  }

  /// `kFitChunksPerThread` chunks per thread of the tracker's arena.
  [[nodiscard]] Index fitChunkCount() const { return kFitChunksPerThread * static_cast<Index>(std::max(mNeighbourStage.maxConcurrency(), 1)); }

  /// The chunk of a loop over \a elements fits -- the road
  /// fit, the refit, the follow's passes: about
  /// `kFitChunksPerThread` chunks per thread, so that the threads end
  /// together, and never above the stages' chunk. The output does not depend
  /// on the chunking (StageOutput.h).
  [[nodiscard]] Index fitChunk(Index elements) const
  {
    const auto chunks = fitChunkCount();
    return std::clamp<Index>((elements + chunks - 1U) / chunks, 1U, mChunkSize);
  }

  /// `hostForEach` in fit chunks (`fitChunk`).
  [[nodiscard]] auto fitForEach()
  {
    return [this](Index count, const auto& body) { mNeighbourStage.forEachElement(count, fitChunk(count), body); };
  }

  /// The chunk of a loop over \a elements that are cheap, as copies: the
  /// elements spread over `kFitChunksPerThread` chunks per thread, however
  /// many.
  [[nodiscard]] Index spreadChunk(Index elements) const
  {
    const auto chunks = fitChunkCount();
    return std::max<Index>((elements + chunks - 1U) / chunks, 1U);
  }

  /// `hostForEach` in spread chunks (`spreadChunk`).
  [[nodiscard]] auto spreadForEach()
  {
    return [this](Index count, const auto& body) { mNeighbourStage.forEachElement(count, spreadChunk(count), body); };
  }

  std::vector<SearchWindowParams> mWindowSettings;
  std::vector<IterationParameters<Det>> mIterations;
  std::size_t mIteration = 0U;
  VertexSeed mBeamLine;
  TrackletInputs<Det> mTrackletInputs;
  SliceObserver mSliceObserver;
  IterationObserver mIterationObserver;
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  Scratch mScratch;
  /// The one arena of the tracker's stages and loops, and of a vertexing
  /// pass's: TBB's workers
  /// stay with it from one loop to the next.
  detail::TaskArena mArena;
  CpuStage<TrackletT> mTrackletStage;
  CpuStage<TripletT> mTripletStage;
  CpuStage<NeighbourT> mNeighbourStage;
  CpuStage<SelectionRecordOf<kMaxClusters>::template RecordT> mRecordStage;
  RoadCount mRoads;
  std::array<MeasurementRange, kLayerCount> mRoadCores{}; ///< the slice's core as the road counting reads it, with a start-layer mask
  Records mRecords;
  GlobalSelection<kMaxClusters> mSelection;
  TrackRefit<Det> mRefit;
  TrackStore<Det> mStore;
  OverlapPartners<Det> mPartners;
  TrackJoin<Det> mJoin;                             ///< the iteration's tracks joined to the earlier ones
  InactiveRows<Det> mInactiveRows;                  ///< inactive layers' rows on the whole track's fit
  BoundedVector<std::uint32_t> mInactiveCandidates; ///< its tracks: the iteration's and the joined earlier ones
  Follow<Det> mFollow;                              ///< the iteration's own tracks'
  Follow<Det> mFollowEarlier;                       ///< the earlier iterations' tracks', with followEarlier
  /// What the follow started from, for an iteration observer only.
  std::vector<BoundTrackState> mFollowStartInner;
  std::vector<BoundTrackState> mFollowStart;
  std::vector<std::uint8_t> mFollowStartCount;
  std::vector<std::uint8_t> mFollowUsed;
  /// What the follow left, before the length cut, for an iteration observer only.
  std::vector<ClusterList<kMaxClusters>> mFollowEndClusters;
  std::vector<std::uint8_t> mFollowEndCount;
  std::vector<BoundTrackState> mFollowEndState;
  SensorCells<Det> mSensorCells;
  BoundedVector<std::uint8_t> mUsed;
  BoundedVector<std::uint8_t> mUsedBeforeSelection; ///< the mask as the selection found it (`dropTracks`)
  BoundedVector<Index> mKeptTracks;                 ///< `dropTracks`' kept tracks
  typename TrackRefit<Det>::Tracks mDropScratch;    ///< `dropTracks`' moved tracks
  BoundedVector<Index> mChangedTracks;              ///< `TrackStore::recommitChanged`'s scratch
  const TrackingGeometry<Det>* mGeometry = nullptr;
  Propagator mPropagator;
  profile::Trace* mTrace = nullptr;
  bool mProfiling = false;
  CapacityEstimator mEstimator;
  std::size_t mRows = 0U;                                        ///< the timeframe's measurement rows
  std::array<std::size_t, CapacityEstimator::kSites> mGranted{}; ///< per stage, the staging granted per thread
  Index mChunkSize = 1U;
  static constexpr Index kFitChunksPerThread = 8U; ///< `fitChunk`
  std::size_t mTrackletHighWater = 0U;
  std::size_t mStagingHighWater = 0U;
  std::size_t mStagingPredictionMisses = 0U;
  std::size_t mTripletHighWater = 0U;
  std::size_t mNeighbourHighWater = 0U;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRACKER_H_
