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

/// \file TrackJoin.h
/// \brief The join of an iteration's tracks to the earlier iterations' ones:
///        two tracks of one particle the follow left apart, where a hard
///        scatter ended the inner one's walk, made one.
///
/// After an iteration's follow and length cut, each of its tracks (a piece)
/// looks for an earlier iteration's track it continues, both ways:
///  - inward, from its inner state, for an earlier track's outermost row;
///  - outward, from its outward state, for an earlier track's innermost row;
/// on every layer one link of the layout reaches from the piece's end layer,
/// or two (a layer between crossed without a row): the window of the
/// iteration's gate there (`follow::predictedWindow`), the grid's rows in
/// it whose time boxes overlap the piece's, cut to the window. The piece's
/// own prediction carries the scattering between the two; a scatter at or
/// near the other track's end leaves that end inside it, where the other
/// track's prediction misses the piece by the deflection times the gap.
///
/// A pair -- the inner and the outer track -- is taken where the two hold no
/// layer in common, their rows fit a chain, their time boxes overlap, their
/// charges agree, and their rows, the inner's then the outer's, refitted
/// (`fitMeasurements`) with a kink budget of the two tracks' own kinks and
/// one more (the inner's iteration's kink gate) pass the road's chi2 cut both ways at their ndf. Its score is
/// the fit's deterioration, chi2(joined) - chi2(inner) - chi2(outer), per
/// degree of freedom the join adds. Each piece keeps its pair of the least
/// score, each earlier track the piece of the least score that chose it
/// (ties by the piece's index): one join per earlier track and per piece.
/// The earlier track takes the joined rows and fit; the piece goes.
///
/// The per-piece work is three steps, each an element functor the host's
/// tracker runs here and the device's in its own stage: `join::Find` lists
/// each piece's pairs past the tests without a fit (a stage's records, in the
/// piece's order of search), `join::Score` fits each pair, `join::Pick` takes
/// each piece's least score in that order; `join::Apply` writes the joins the
/// arbitration keeps. The fits apart from the search run one pair per thread
/// on a device, not the few pieces of a group that found one.

#ifndef GTRK_TRACKING_TRACKJOIN_H_
#define GTRK_TRACKING_TRACKJOIN_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/KFTrackFitter.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackRefit.h"

#include "GPUCommonMath.h"

#include <memlayout/wrapper.h>

#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <span>
#include <type_traits>
#include <utility>

namespace o2::gtrk
{

namespace join
{

inline constexpr std::uint32_t kNone = std::numeric_limits<std::uint32_t>::max();

/// A piece's pair: its partner (kNone without), which of the two is the
/// inner track, and the score.
struct Choice {
  std::uint32_t partner = kNone;
  bool pieceInner = false;
  Scalar score = std::numeric_limits<Scalar>::infinity();
};

/// A piece's pair past the tests without a fit (`Find`): the piece (from the
/// iteration's first), its partner, and whether the piece is the inner track.
struct PairCandidate {
  // As for TrackletT: upstream's `apply` macro also uses the member names as
  // enum constants.
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
  template <template <class> class C>
  struct RecordT {
    MEMLAYOUT_APPLY(RecordT, piece, partner, pieceInner)
    C<std::uint32_t> piece;
    C<std::uint32_t> partner;
    C<std::uint8_t> pieceInner;
  };
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
};

/// The layers one link, or two, of the layout reach from \a layer in
/// \a direction, a bit each (`IterationParameters::layerLinksOut`/`In`).
template <DetectorLayout Det>
GPUhdi() std::uint64_t reached(const IterationParameters<Det>& parameters, LayerID layer, PropagationDirection direction)
{
  const auto* links = direction == PropagationDirection::kForward ? parameters.layerLinksOut : parameters.layerLinksIn;
  const std::uint64_t one = links[layer];
  std::uint64_t two = 0U;
  for (std::uint64_t left = one; left != 0U; left &= left - 1U) {
    two |= links[o2::gpu::CAMath::Ctz64(left)];
  }
  return one | two;
}

/// The joined rows of \a inner then \a outer of \a tracks into \a rows
/// (\a count of them), refitted with a kink budget of the two tracks' kinks
/// and one more (\a parameters' kink gate).
template <DetectorLayout Det, typename Tracks>
GPUhdi() TrackFitResult joinedFit(const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, const Tracks& tracks, std::uint32_t inner, std::uint32_t outer, ClusterList<maxChainLayers<Det>>& rows, std::uint8_t& count)
{
  count = 0U;
  for (std::uint8_t i = 0; i < tracks.count[inner]; ++i) {
    rows.index[count++] = tracks.clusters[inner].index[i];
  }
  for (std::uint8_t i = 0; i < tracks.count[outer]; ++i) {
    rows.index[count++] = tracks.clusters[outer].index[i];
  }
  const KinkBudget budget{.max = static_cast<int>(tracks.kinks[inner]) + static_cast<int>(tracks.kinks[outer]) + 1, .gate = parameters.kinkChi2Min};
  return fitMeasurements(propagator, geometry, detail::LinkSteps{}, RowChain<Det>{.frame = &frame, .rows = rows.index}, static_cast<int>(count), budget);
}

/// Whether the pair \a inner, \a outer of \a tracks passes the tests
/// without a fit: no layer in common, rows to fit one chain, time boxes
/// overlapping, charges agreeing.
template <DetectorLayout Det, typename Tracks>
GPUhdi() bool pairAdmitted(const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Tracks& tracks, std::uint32_t inner, std::uint32_t outer)
{
  const auto& measurements = frame.measurements();
  std::uint64_t innerLayers = 0U;
  std::uint64_t outerLayers = 0U;
  for (std::uint8_t i = 0; i < tracks.count[inner]; ++i) {
    innerLayers |= std::uint64_t{1} << geometry.surface(measurements.surface[tracks.clusters[inner].index[i]]).layer;
  }
  for (std::uint8_t i = 0; i < tracks.count[outer]; ++i) {
    outerLayers |= std::uint64_t{1} << geometry.surface(measurements.surface[tracks.clusters[outer].index[i]]).layer;
  }
  if ((innerLayers & outerLayers) != 0U || static_cast<std::size_t>(tracks.count[inner]) + tracks.count[outer] > maxChainLayers<Det>) {
    return false;
  }
  const TimeEstPS innerBox{.start = tracks.timePS[inner], .width = tracks.widthPS[inner]};
  const TimeEstPS outerBox{.start = tracks.timePS[outer], .width = tracks.widthPS[outer]};
  return innerBox.isCompatible(outerBox) && (tracks.outwardState[inner].qOverP() > 0.F) == (tracks.state[outer].qOverP() > 0.F);
}

/// The admitted pair \a inner, \a outer of \a tracks (`pairAdmitted`): its
/// score, infinity where its fit is refused. \a parameters holds every
/// iteration's, the inner track's taken.
template <DetectorLayout Det, typename Tracks>
GPUhdi() Scalar pairScore(const IterationParameters<Det>* parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, const Tracks& tracks, std::uint32_t inner, std::uint32_t outer)
{
  constexpr Scalar kRefused = std::numeric_limits<Scalar>::infinity();
  const auto& own = parameters[tracks.iteration[inner]];
  ClusterList<maxChainLayers<Det>> rows;
  std::uint8_t count = 0U;
  const auto fit = joinedFit(own, frame, geometry, propagator, tracks, inner, outer, rows, count);
  if (!follow::passesRoadChi2(own, fit)) {
    return kRefused;
  }
  const int outwardNdf = (2 * static_cast<int>(count)) - fittedParameters(propagator) - (2 * static_cast<int>(fit.outwardKinks));
  if (outwardNdf <= 0 || outwardNdf >= static_cast<int>(std::size(own.roadChi2MaxByNdf)) || !(fit.outwardChi2 <= own.roadChi2MaxByNdf[outwardNdf])) {
    return kRefused;
  }
  const int added = fit.ndf - tracks.ndf[inner] - tracks.ndf[outer];
  if (added <= 0) {
    return kRefused;
  }
  return (fit.chi2 - tracks.chi2[inner] - tracks.chi2[outer]) / static_cast<Scalar>(added);
}

/// Per piece \a at (the track \a trackOffset + \a at of \a tracks), its
/// admitted pairs (`pairAdmitted`), both ways, emitted in the order of the
/// search: inward, from its inner state, for an earlier track's outermost row
/// (\a endRows, its track \a endTracks); outward, from its outward state,
/// for an earlier track's innermost row (\a startRows, \a startTracks) --
/// each \a endCount long, ordered by row then track. A window's rows are
/// those lists' (`ListedRows`): an earlier track ends on few of a layer's
/// rows. A row where several end gives the lowest of them.
template <DetectorLayout Det, typename Tracks>
struct Find {
  using Record = typename BoundedColumns<PairCandidate::template RecordT>::Record;

  const IterationParameters<Det>* parameters = nullptr; ///< every iteration's, by the tracks' iteration
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  Propagator propagator;
  Tracks tracks;
  std::uint32_t trackOffset = 0U;
  const std::uint32_t* endRows = nullptr;
  const std::uint32_t* endTracks = nullptr;
  const std::uint32_t* startRows = nullptr;
  const std::uint32_t* startTracks = nullptr;
  std::uint32_t endCount = 0U;
  const std::uint32_t* maxCandidateWidthPS = nullptr; ///< per layer

  template <typename Emitter>
  GPUhdi() void operator()(std::uint32_t at, Emitter& emitter) const
  {
    const auto piece = trackOffset + at;
    if (tracks.status[piece] != PropagationStatus::kOk || tracks.count[piece] == 0U) {
      return;
    }
    const auto& measurements = frame->measurements();
    const auto& own = parameters[tracks.iteration[piece]];
    const TimeEstPS box{.start = tracks.timePS[piece], .width = tracks.widthPS[piece]};
    std::uint64_t pieceLayers = 0U;
    for (std::uint8_t i = 0; i < tracks.count[piece]; ++i) {
      pieceLayers |= std::uint64_t{1} << geometry->surface(measurements.surface[tracks.clusters[piece].index[i]]).layer;
    }
    for (int side = 0; side < 2; ++side) {
      const bool inward = side == 0;
      const auto direction = inward ? PropagationDirection::kBackward : PropagationDirection::kForward;
      const BoundTrackState& state = inward ? tracks.state[piece] : tracks.outwardState[piece];
      const auto factors = factorCovariance(state.cov);
      const auto endRow = inward ? tracks.clusters[piece].index[0] : tracks.clusters[piece].index[tracks.count[piece] - 1U];
      const ListedRows listed{.rows = inward ? endRows : startRows, .count = endCount};
      const std::uint32_t* partners = inward ? endTracks : startTracks;
      for (std::uint64_t left = reached(own, geometry->surface(measurements.surface[endRow]).layer, direction) & ~pieceLayers; left != 0U; left &= left - 1U) {
        const auto target = static_cast<LayerID>(o2::gpu::CAMath::Ctz64(left));
        const auto prediction = predictLayerCrossing(propagator, own, geometry->planeOf(state), state, target, direction);
        if (prediction.status != PropagationStatus::kOk) {
          continue;
        }
        Scalar phiCentre = 0.F;
        Scalar secondCentre = 0.F;
        follow::windowCentre(own, target, prediction, phiCentre, secondCentre);
        CandidateWindow window{
          .earliestPS = static_cast<std::int64_t>(box.lower()) - static_cast<std::int64_t>(maxCandidateWidthPS[target]) + 1,
          .latestPS = static_cast<std::int64_t>(box.upper()) - 1,
          .afterPS = box.lower(),
          .phiCentre = phiCentre,
          .phiHalfWidth = kPi<Scalar>,
          .secondCentre = secondCentre,
          .secondHalfWidth = follow::kOpenBound,
        };
        if (!follow::predictedWindow(propagator, own, *geometry, detail::LinkSteps{}, state, factors, target, direction, prediction, window)) {
          continue;
        }
        frame->grid().visitCandidates(
          target, window, frame->grid().geometry().allBlocks(), [&](std::uint32_t row, std::uint32_t entry) {
            const auto measurement = frame->measurement(row);
            if (!spatiallyCompatible(window, measurement.phi, secondCoordinate<Det>(target, measurement.r, measurement.z))) {
              return;
            }
            const auto partner = partners[entry];
            if (pairAdmitted(*frame, *geometry, tracks, inward ? partner : piece, inward ? piece : partner)) {
              emitter.emit(Record{at, partner, static_cast<std::uint8_t>(inward ? 0U : 1U)});
            }
          },
          listed);
      }
    }
  }
};

/// Per pair \a k (`Find`'s records: \a pieces from \a trackOffset,
/// \a partners, \a pieceInner): its score (`pairScore`) into \a scores[k].
template <DetectorLayout Det, typename Tracks>
struct Score {
  const IterationParameters<Det>* parameters = nullptr; ///< every iteration's, by the tracks' iteration
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  Propagator propagator;
  Tracks tracks;
  std::uint32_t trackOffset = 0U;
  const std::uint32_t* pieces = nullptr;
  const std::uint32_t* partners = nullptr;
  const std::uint8_t* pieceInner = nullptr;
  Scalar* scores = nullptr;

  GPUhdi() void operator()(std::uint32_t k) const
  {
    const auto piece = trackOffset + pieces[k];
    const auto inner = pieceInner[k] != 0U ? piece : partners[k];
    const auto outer = pieceInner[k] != 0U ? partners[k] : piece;
    scores[k] = pairScore(parameters, *frame, *geometry, propagator, tracks, inner, outer);
  }
};

/// Per piece \a at, its pair of the least score among its pairs
/// (\a offsets[at] up to \a offsets[at + 1] of \a partners, \a pieceInner,
/// \a scores), ties by the lower partner, then the first found, into
/// \a choices[at]; none where every score is infinite.
struct Pick {
  const std::uint32_t* offsets = nullptr;
  const std::uint32_t* partners = nullptr;
  const std::uint8_t* pieceInner = nullptr;
  const Scalar* scores = nullptr;
  Choice* choices = nullptr;

  GPUhdi() void operator()(std::uint32_t at) const
  {
    Choice best;
    for (auto k = offsets[at]; k < offsets[at + 1U]; ++k) {
      const auto partner = partners[k];
      const auto score = scores[k];
      if (score < best.score || (score == best.score && best.partner != kNone && partner < best.partner)) {
        best = {.partner = partner, .pieceInner = pieceInner[k] != 0U, .score = score};
      }
    }
    choices[at] = best;
  }
};

/// Per joined piece \a i (\a joined[i], an offset from \a trackOffset, its
/// pair in \a choices at the same offset): the joined rows and fit written
/// over its partner in \a tracks.
template <DetectorLayout Det, typename Tracks>
struct Apply {
  const IterationParameters<Det>* parameters = nullptr; ///< every iteration's, by the tracks' iteration
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  Propagator propagator;
  Tracks tracks;
  std::uint32_t trackOffset = 0U;
  const Choice* choices = nullptr;
  const std::uint32_t* joined = nullptr;

  GPUhdi() void operator()(std::uint32_t i) const
  {
    const auto at = joined[i];
    const auto& choice = choices[at];
    const auto piece = trackOffset + at;
    const auto inner = choice.pieceInner ? piece : choice.partner;
    const auto outer = choice.pieceInner ? choice.partner : piece;
    ClusterList<maxChainLayers<Det>> rows;
    std::uint8_t count = 0U;
    const auto fit = joinedFit(parameters[tracks.iteration[inner]], *frame, *geometry, propagator, tracks, inner, outer, rows, count);
    writeFit(tracks, choice.partner, fit, rows, count, trackTimeBox(frame->measurements(), rows, count));
  }
};

} // namespace join

template <DetectorLayout Det>
class TrackJoin
{
 public:
  static constexpr std::size_t kLayerCount = static_cast<std::size_t>(Detector<Det>::nLayers);
  static constexpr std::uint32_t kNone = join::kNone;
  using Tracks = typename TrackRefit<Det>::Tracks;
  using Pairs = StageOutput<join::PairCandidate::template RecordT>;

  /// On \a arena, the tracker's.
  explicit TrackJoin(detail::TaskArena& arena) : mStage(arena) {}

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool = std::move(pool);
    mStage.setMemoryPool(mMemoryPool);
    mPairs.setMemoryPool(mMemoryPool);
    clearStorage();
  }
  void clearStorage() noexcept
  {
    mStage.releaseStorage();
    mPairs.clearStorage();
    deepVectorClear(mMemoryPool.get(), mListed, mKeys, mEndRows, mEndTracks, mStartRows, mStartTracks, mScores, mChoices, mBest, mJoined, mChanged);
  }

  /// The join of \a tracks' pieces from \a trackOffset on to the tracks
  /// before it, each with its own iteration's \a parameters (every
  /// iteration's), each layer's candidates bound by \a maxCandidateWidthPS:
  /// the earlier tracks rewritten, the pieces joined listed (`joined`,
  /// offsets from \a trackOffset, ascending) for the caller to drop, the
  /// earlier tracks changed listed (`changed`, ascending). \a forEach(n, body)
  /// runs body(i) for every i below n; \a chunk is the pairs' stage's chunk.
  template <typename ForEach>
  [[nodiscard]] StageStatus run(ForEach forEach, std::uint32_t chunk, const IterationParameters<Det>* parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, Tracks& tracks, std::uint32_t trackOffset, std::span<const std::uint32_t, kLayerCount> maxCandidateWidthPS)
  {
    return guardOutOfMemory([&] {
      const auto view = tracks.view();
      using View = std::remove_const_t<decltype(view)>;
      const auto size = static_cast<std::uint32_t>(tracks.size());
      const auto pieces = size - trackOffset;
      resizeBoundedVector(mJoined, 0U, mMemoryPool.get());
      resizeBoundedVector(mChanged, 0U, mMemoryPool.get());
      if (trackOffset == 0U || pieces == 0U) {
        return StageStatus::kComplete;
      }
      // The earlier tracks' ends and starts, by row, then track: each a key
      // of the row above the track, sorted, then split.
      std::uint32_t ends = 0U;
      for (std::uint32_t track = 0; track < trackOffset; ++track) {
        ends += view.status[track] == PropagationStatus::kOk && view.count[track] != 0U ? 1U : 0U;
      }
      resizeBoundedVector(mListed, ends, mMemoryPool.get());
      for (std::uint32_t track = 0, at = 0; track < trackOffset; ++track) {
        if (view.status[track] == PropagationStatus::kOk && view.count[track] != 0U) {
          mListed[at++] = track;
        }
      }
      resizeBoundedVector(mKeys, ends, mMemoryPool.get());
      resizeBoundedVector(mEndRows, ends, mMemoryPool.get());
      resizeBoundedVector(mEndTracks, ends, mMemoryPool.get());
      resizeBoundedVector(mStartRows, ends, mMemoryPool.get());
      resizeBoundedVector(mStartTracks, ends, mMemoryPool.get());
      for (const bool outermost : {true, false}) {
        forEach(ends, [&](std::uint32_t at) {
          const auto track = mListed[at];
          const std::uint64_t row = view.clusters[track].index[outermost ? view.count[track] - 1U : 0U];
          mKeys[at] = (row << 32U) | track;
        });
        mStage.sortKeys(std::span<std::uint64_t>{mKeys});
        auto& rows = outermost ? mEndRows : mStartRows;
        auto& listedTracks = outermost ? mEndTracks : mStartTracks;
        forEach(ends, [&](std::uint32_t at) {
          rows[at] = static_cast<std::uint32_t>(mKeys[at] >> 32U);
          listedTracks[at] = static_cast<std::uint32_t>(mKeys[at]);
        });
      }
      // Each piece's pairs, their scores, its pair of the least score.
      mPairs.prepare(pieces);
      const auto status = mStage.run(mPairs, chunk, join::Find<Det, View>{.parameters = parameters, .frame = &frame, .geometry = &geometry, .propagator = propagator, .tracks = view, .trackOffset = trackOffset, .endRows = mEndRows.data(), .endTracks = mEndTracks.data(), .startRows = mStartRows.data(), .startTracks = mStartTracks.data(), .endCount = ends, .maxCandidateWidthPS = maxCandidateWidthPS.data()});
      if (status != StageStatus::kComplete) {
        return status;
      }
      const auto pairs = std::as_const(mPairs).view();
      const auto total = mPairs.size();
      resizeBoundedVector(mScores, total, mMemoryPool.get());
      forEach(total, join::Score<Det, View>{.parameters = parameters, .frame = &frame, .geometry = &geometry, .propagator = propagator, .tracks = view, .trackOffset = trackOffset, .pieces = pairs.piece, .partners = pairs.partner, .pieceInner = pairs.pieceInner, .scores = mScores.data()});
      resizeBoundedVector(mChoices, pieces, mMemoryPool.get());
      forEach(pieces, join::Pick{.offsets = mPairs.offsets().data(), .partners = pairs.partner, .pieceInner = pairs.pieceInner, .scores = mScores.data(), .choices = mChoices.data()});
      // Each earlier track the piece of the least score that chose it.
      resizeBoundedVector(mBest, trackOffset, mMemoryPool.get());
      std::ranges::fill(mBest, kNone);
      for (std::uint32_t at = 0; at < pieces; ++at) {
        const auto& choice = mChoices[at];
        if (choice.partner == kNone) {
          continue;
        }
        auto& best = mBest[choice.partner];
        if (best == kNone || choice.score < mChoices[best].score) {
          best = at;
        }
      }
      std::uint32_t joined = 0U;
      for (std::uint32_t at = 0; at < pieces; ++at) {
        joined += mChoices[at].partner != kNone && mBest[mChoices[at].partner] == at ? 1U : 0U;
      }
      resizeBoundedVector(mJoined, joined, mMemoryPool.get());
      for (std::uint32_t at = 0, i = 0; at < pieces; ++at) {
        if (mChoices[at].partner != kNone && mBest[mChoices[at].partner] == at) {
          mJoined[i++] = at;
        }
      }
      // The joined tracks, written over the earlier ones.
      forEach(joined, join::Apply<Det, View>{.parameters = parameters, .frame = &frame, .geometry = &geometry, .propagator = propagator, .tracks = view, .trackOffset = trackOffset, .choices = mChoices.data(), .joined = mJoined.data()});
      resizeBoundedVector(mChanged, joined, mMemoryPool.get());
      for (std::uint32_t i = 0; i < joined; ++i) {
        mChanged[i] = mChoices[mJoined[i]].partner;
      }
      std::ranges::sort(mChanged);
      return StageStatus::kComplete;
    });
  }

  /// The last run's joined pieces, offsets from its trackOffset, ascending.
  [[nodiscard]] std::span<const std::uint32_t> joined() const { return mJoined; }
  /// The last run's earlier tracks that took a piece, ascending.
  [[nodiscard]] std::span<const std::uint32_t> changed() const { return mChanged; }

 private:
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  CpuStage<join::PairCandidate::template RecordT> mStage;
  Pairs mPairs;                            ///< each piece's pairs past the tests without a fit
  BoundedVector<std::uint32_t> mListed;    ///< the earlier tracks with rows, ascending
  BoundedVector<std::uint64_t> mKeys;      ///< their ends (or starts) above their tracks
  BoundedVector<std::uint32_t> mEndRows;   ///< each earlier track's outermost row, by row then track
  BoundedVector<std::uint32_t> mEndTracks; ///< and its track
  BoundedVector<std::uint32_t> mStartRows; ///< each earlier track's innermost row, by row then track
  BoundedVector<std::uint32_t> mStartTracks;
  BoundedVector<Scalar> mScores;        ///< per pair
  BoundedVector<join::Choice> mChoices; ///< per piece
  BoundedVector<std::uint32_t> mBest;   ///< per earlier track, the winning piece
  BoundedVector<std::uint32_t> mJoined;
  BoundedVector<std::uint32_t> mChanged;
};

} // namespace o2::gtrk

#endif
