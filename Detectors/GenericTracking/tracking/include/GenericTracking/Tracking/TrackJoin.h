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
/// The per-piece work is the element functors `join::Choose` and
/// `join::Apply`, which the host's tracker runs here and the device's in its
/// own stage; the arbitration between them reads their outputs.

#ifndef GTRK_TRACKING_TRACKJOIN_H_
#define GTRK_TRACKING_TRACKJOIN_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/KFTrackFitter.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackRefit.h"

#include "GPUCommonMath.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <utility>

namespace o2::gtrk
{

namespace join
{

inline constexpr std::uint32_t kNone = std::numeric_limits<std::uint32_t>::max();

/// An earlier track's end: its outermost (or innermost) row and the track;
/// a list of them ordered by row, then track.
struct End {
  std::uint32_t row = 0U;
  std::uint32_t track = 0U;
  GPUhdi() bool operator<(const End& other) const { return row != other.row ? row < other.row : track < other.track; }
};

/// A piece's pair: its partner (kNone without), which of the two is the
/// inner track, and the score.
struct Choice {
  std::uint32_t partner = kNone;
  bool pieceInner = false;
  Scalar score = std::numeric_limits<Scalar>::infinity();
};

/// The first of \a count \a ends at \a row (the lowest track ending there),
/// \a count where none is.
GPUhdi() std::uint32_t firstAt(const End* ends, std::uint32_t count, std::uint32_t row)
{
  std::uint32_t low = 0U;
  std::uint32_t high = count;
  while (low < high) {
    const std::uint32_t middle = low + ((high - low) / 2U);
    if (ends[middle].row < row) {
      low = middle + 1U;
    } else {
      high = middle;
    }
  }
  return low < count && ends[low].row == row ? low : count;
}

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

/// The pair \a inner, \a outer of \a tracks: its score, infinity where it is
/// refused. \a parameters holds every iteration's, the inner track's taken.
template <DetectorLayout Det, typename Tracks>
GPUhdi() Scalar pairScore(const IterationParameters<Det>* parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, const Tracks& tracks, std::uint32_t inner, std::uint32_t outer)
{
  constexpr Scalar kRefused = std::numeric_limits<Scalar>::infinity();
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
    return kRefused;
  }
  const TimeEstPS innerBox{.start = tracks.timePS[inner], .width = tracks.widthPS[inner]};
  const TimeEstPS outerBox{.start = tracks.timePS[outer], .width = tracks.widthPS[outer]};
  if (!innerBox.isCompatible(outerBox) || (tracks.outwardState[inner].qOverP() > 0.F) != (tracks.state[outer].qOverP() > 0.F)) {
    return kRefused;
  }
  const auto& own = parameters[tracks.iteration[inner]];
  ClusterList<maxChainLayers<Det>> rows;
  std::uint8_t count = 0U;
  const auto fit = joinedFit(own, frame, geometry, propagator, tracks, inner, outer, rows, count);
  if (!follow::passesRoadChi2(own, fit)) {
    return kRefused;
  }
  const int outwardNdf = (2 * static_cast<int>(count)) - 5 - (2 * static_cast<int>(fit.outwardKinks));
  if (outwardNdf <= 0 || outwardNdf >= static_cast<int>(std::size(own.roadChi2MaxByNdf)) || !(fit.outwardChi2 <= own.roadChi2MaxByNdf[outwardNdf])) {
    return kRefused;
  }
  const int added = fit.ndf - tracks.ndf[inner] - tracks.ndf[outer];
  if (added <= 0) {
    return kRefused;
  }
  return (fit.chi2 - tracks.chi2[inner] - tracks.chi2[outer]) / static_cast<Scalar>(added);
}

/// Per piece \a at (the track \a trackOffset + \a at of \a tracks), its pair
/// of the least score, both ways, into \a choices[at]: inward, from its inner
/// state, for an earlier track's outermost row (\a ends); outward, from its
/// outward state, for an earlier track's innermost row (\a starts) -- each
/// \a endCount long, ordered by row then track.
template <DetectorLayout Det, typename Tracks>
struct Choose {
  const IterationParameters<Det>* parameters = nullptr; ///< every iteration's, by the tracks' iteration
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  Propagator propagator;
  Tracks tracks;
  std::uint32_t trackOffset = 0U;
  const End* ends = nullptr;
  const End* starts = nullptr;
  std::uint32_t endCount = 0U;
  const std::uint32_t* maxCandidateWidthPS = nullptr; ///< per layer
  Choice* choices = nullptr;

  GPUhdi() void operator()(std::uint32_t at) const { choices[at] = choose(trackOffset + at); }

  [[nodiscard]] GPUhdi() Choice choose(std::uint32_t piece) const
  {
    Choice best;
    if (tracks.status[piece] != PropagationStatus::kOk || tracks.count[piece] == 0U) {
      return best;
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
      const End* partners = inward ? ends : starts;
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
        frame->grid().visitCandidates(target, window, frame->grid().geometry().allBlocks(), [&](std::uint32_t row) {
          const auto measurement = frame->measurement(row);
          if (!spatiallyCompatible(window, measurement.phi, secondCoordinate<Det>(target, measurement.r, measurement.z))) {
            return;
          }
          const auto found = firstAt(partners, endCount, row);
          if (found == endCount) {
            return;
          }
          const auto partner = partners[found].track;
          const auto inner = inward ? partner : piece;
          const auto outer = inward ? piece : partner;
          const auto score = pairScore(parameters, *frame, *geometry, propagator, tracks, inner, outer);
          if (score < best.score || (score == best.score && best.partner != kNone && partner < best.partner)) {
            best = {.partner = partner, .pieceInner = !inward, .score = score};
          }
        });
      }
    }
    return best;
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

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool = std::move(pool);
    clearStorage();
  }
  void clearStorage() noexcept { deepVectorClear(mMemoryPool.get(), mEnds, mStarts, mChoices, mBest, mJoined, mChanged); }

  /// The join of \a tracks' pieces from \a trackOffset on to the tracks
  /// before it, each with its own iteration's \a parameters (every
  /// iteration's), each layer's candidates bound by \a maxCandidateWidthPS:
  /// the earlier tracks rewritten, the pieces joined listed (`joined`,
  /// offsets from \a trackOffset, ascending) for the caller to drop, the
  /// earlier tracks changed listed (`changed`, ascending).
  template <typename ForEach>
  [[nodiscard]] StageStatus run(ForEach forEach, const IterationParameters<Det>* parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, Tracks& tracks, std::uint32_t trackOffset, std::span<const std::uint32_t, kLayerCount> maxCandidateWidthPS)
  {
    return guardOutOfMemory([&] {
      const auto view = tracks.view();
      const auto size = static_cast<std::uint32_t>(tracks.size());
      const auto pieces = size - trackOffset;
      resizeBoundedVector(mJoined, 0U, mMemoryPool.get());
      resizeBoundedVector(mChanged, 0U, mMemoryPool.get());
      if (trackOffset == 0U || pieces == 0U) {
        return StageStatus::kComplete;
      }
      // The earlier tracks' ends and starts, by row, then track.
      std::uint32_t ends = 0U;
      for (std::uint32_t track = 0; track < trackOffset; ++track) {
        ends += view.status[track] == PropagationStatus::kOk && view.count[track] != 0U ? 1U : 0U;
      }
      resizeBoundedVector(mEnds, ends, mMemoryPool.get());
      resizeBoundedVector(mStarts, ends, mMemoryPool.get());
      for (std::uint32_t track = 0, at = 0; track < trackOffset; ++track) {
        if (view.status[track] == PropagationStatus::kOk && view.count[track] != 0U) {
          mEnds[at] = {.row = view.clusters[track].index[view.count[track] - 1U], .track = track};
          mStarts[at++] = {.row = view.clusters[track].index[0], .track = track};
        }
      }
      std::sort(mEnds.begin(), mEnds.end());
      std::sort(mStarts.begin(), mStarts.end());
      resizeBoundedVector(mChoices, pieces, mMemoryPool.get());
      forEach(pieces, join::Choose<Det, decltype(view)>{.parameters = parameters, .frame = &frame, .geometry = &geometry, .propagator = propagator, .tracks = view, .trackOffset = trackOffset, .ends = mEnds.data(), .starts = mStarts.data(), .endCount = ends, .maxCandidateWidthPS = maxCandidateWidthPS.data(), .choices = mChoices.data()});
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
      forEach(joined, join::Apply<Det, decltype(view)>{.parameters = parameters, .frame = &frame, .geometry = &geometry, .propagator = propagator, .tracks = view, .trackOffset = trackOffset, .choices = mChoices.data(), .joined = mJoined.data()});
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
  BoundedVector<join::End> mEnds;       ///< each earlier track's outermost row, by row then track
  BoundedVector<join::End> mStarts;     ///< and its innermost
  BoundedVector<join::Choice> mChoices; ///< per piece
  BoundedVector<std::uint32_t> mBest;   ///< per earlier track, the winning piece
  BoundedVector<std::uint32_t> mJoined;
  BoundedVector<std::uint32_t> mChanged;
};

} // namespace o2::gtrk

#endif
