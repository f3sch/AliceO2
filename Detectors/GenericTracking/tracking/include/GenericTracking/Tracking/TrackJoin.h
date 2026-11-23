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
/// (`fitMeasurements`, the inner track's iteration's crossing rule) with a kink
/// budget of the two tracks' own kinks and one more (the inner's iteration's
/// kink gate) pass the road's chi2 cut both ways at their ndf. Its score is
/// the fit's deterioration, chi2(joined) - chi2(inner) - chi2(outer), per
/// degree of freedom the join adds. Each piece keeps its pair of the least
/// score, each earlier track the piece of the least score that chose it
/// (ties by the piece's index): one join per earlier track and per piece.
/// The earlier track takes the joined rows and fit; the piece goes.

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

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <utility>

namespace o2::gtrk
{

template <DetectorLayout Det>
class TrackJoin
{
 public:
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  static constexpr std::size_t kLayerCount = static_cast<std::size_t>(Detector<Det>::nLayers);
  static constexpr std::uint32_t kNone = std::numeric_limits<std::uint32_t>::max();
  using Tracks = typename TrackRefit<Det>::Tracks;

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool = std::move(pool);
    clearStorage();
  }
  void clearStorage() noexcept { deepVectorClear(mMemoryPool.get(), mEnds, mStarts, mChoices, mBest, mJoined, mChanged); }

  /// The join of \a tracks' pieces from \a trackOffset on to the tracks
  /// before it, each iteration's parameters \a parametersOf(i), each layer's
  /// candidates bound by \a maxCandidateWidthPS: the earlier tracks
  /// rewritten, the pieces joined listed (`joined`, offsets from
  /// \a trackOffset, ascending) for the caller to drop, the earlier tracks
  /// changed listed (`changed`, ascending).
  template <typename ForEach, typename ParametersOf>
  [[nodiscard]] StageStatus run(ForEach forEach, const ParametersOf& parametersOf, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, Tracks& tracks, std::uint32_t trackOffset, std::span<const std::uint32_t, kLayerCount> maxCandidateWidthPS)
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
      // The earlier tracks' ends and starts: (row, track), by row.
      std::uint32_t ends = 0U;
      for (std::uint32_t track = 0; track < trackOffset; ++track) {
        ends += view.status[track] == PropagationStatus::kOk && view.count[track] != 0U ? 1U : 0U;
      }
      resizeBoundedVector(mEnds, ends, mMemoryPool.get());
      resizeBoundedVector(mStarts, ends, mMemoryPool.get());
      for (std::uint32_t track = 0, at = 0; track < trackOffset; ++track) {
        if (view.status[track] == PropagationStatus::kOk && view.count[track] != 0U) {
          mEnds[at] = {view.clusters[track].index[view.count[track] - 1U], track};
          mStarts[at++] = {view.clusters[track].index[0], track};
        }
      }
      std::ranges::sort(mEnds);
      std::ranges::sort(mStarts);
      resizeBoundedVector(mChoices, pieces, mMemoryPool.get());
      forEach(pieces, [&](std::uint32_t at) { mChoices[at] = choose(parametersOf, frame, geometry, propagator, view, trackOffset + at, maxCandidateWidthPS); });
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
      forEach(static_cast<std::uint32_t>(mJoined.size()), [&](std::uint32_t i) {
        const auto at = mJoined[i];
        const auto& choice = mChoices[at];
        const auto piece = trackOffset + at;
        const auto inner = choice.pieceInner ? piece : choice.partner;
        const auto outer = choice.pieceInner ? choice.partner : piece;
        ClusterList<kMaxClusters> rows;
        std::uint8_t count = 0U;
        const auto fit = joinedFit(parametersOf(view.iteration[inner]), frame, geometry, propagator, view, inner, outer, rows, count);
        writeFit(view, choice.partner, fit, rows, count, trackTimeBox(frame.measurements(), rows, count));
      });
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
  /// A piece's pair: its partner (kNone without), which of the two is the
  /// inner track, and the score.
  struct Choice {
    std::uint32_t partner = kNone;
    bool pieceInner = false;
    Scalar score = std::numeric_limits<Scalar>::infinity();
  };

  /// The layers one link, or two, of the layout reach from \a layer in
  /// \a direction, a bit each.
  [[nodiscard]] static std::uint64_t reachedLayers(LayerID layer, PropagationDirection direction)
  {
    const bool forward = direction == PropagationDirection::kForward;
    const auto step = [&](std::uint64_t from) {
      std::uint64_t to = 0U;
      for (const auto& link : Det::layerLinks) {
        const auto source = forward ? link.from : link.to;
        const auto target = forward ? link.to : link.from;
        if (((from >> source) & 1U) != 0U) {
          to |= std::uint64_t{1} << target;
        }
      }
      return to;
    };
    const auto one = step(std::uint64_t{1} << layer);
    return one | step(one);
  }

  /// The joined rows of \a inner then \a outer into \a rows (\a count of
  /// them), refitted.
  template <typename View>
  [[nodiscard]] static TrackFitResult joinedFit(const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, const View& view, std::uint32_t inner, std::uint32_t outer, ClusterList<kMaxClusters>& rows, std::uint8_t& count)
  {
    count = 0U;
    for (std::uint8_t i = 0; i < view.count[inner]; ++i) {
      rows.index[count++] = view.clusters[inner].index[i];
    }
    for (std::uint8_t i = 0; i < view.count[outer]; ++i) {
      rows.index[count++] = view.clusters[outer].index[i];
    }
    const KinkBudget budget{.max = static_cast<int>(view.kinks[inner]) + static_cast<int>(view.kinks[outer]) + 1, .gate = parameters.kinkChi2Min};
    return withCrossingRule(parameters.crossingRule, [&]<CrossingRule Rule>(std::integral_constant<CrossingRule, Rule>) {
      return fitMeasurements<Rule>(propagator, geometry, detail::LinkSteps{}, RowChain<Det>{.frame = &frame, .rows = rows.index}, static_cast<int>(count), budget);
    });
  }

  /// \a piece's pair of the least score, both ways.
  template <typename ParametersOf, typename View>
  [[nodiscard]] Choice choose(const ParametersOf& parametersOf, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, const View& view, std::uint32_t piece, std::span<const std::uint32_t, kLayerCount> maxCandidateWidthPS) const
  {
    Choice best;
    if (view.status[piece] != PropagationStatus::kOk || view.count[piece] == 0U) {
      return best;
    }
    const auto& measurements = frame.measurements();
    const auto& parameters = parametersOf(view.iteration[piece]);
    const TimeEstPS box{.start = view.timePS[piece], .width = view.widthPS[piece]};
    std::uint64_t pieceLayers = 0U;
    for (std::uint8_t i = 0; i < view.count[piece]; ++i) {
      pieceLayers |= std::uint64_t{1} << geometry.surface(measurements.surface[view.clusters[piece].index[i]]).layer;
    }
    for (const auto direction : {PropagationDirection::kBackward, PropagationDirection::kForward}) {
      const bool inward = direction == PropagationDirection::kBackward;
      const BoundTrackState& state = inward ? view.state[piece] : view.outwardState[piece];
      const auto factors = factorCovariance(state.cov);
      const auto endRow = inward ? view.clusters[piece].index[0] : view.clusters[piece].index[view.count[piece] - 1U];
      const auto& partners = inward ? mEnds : mStarts;
      for (std::uint64_t left = reachedLayers(geometry.surface(measurements.surface[endRow]).layer, direction) & ~pieceLayers; left != 0U; left &= left - 1U) {
        const auto target = static_cast<LayerID>(std::countr_zero(left));
        const auto prediction = predictLayerCrossing(propagator, parameters, geometry.planeOf(state), state, target, direction);
        if (prediction.status != PropagationStatus::kOk) {
          continue;
        }
        Scalar phiCentre = 0.F;
        Scalar secondCentre = 0.F;
        follow::windowCentre(parameters, target, prediction, phiCentre, secondCentre);
        CandidateWindow window{
          .earliestPS = static_cast<std::int64_t>(box.lower()) - static_cast<std::int64_t>(maxCandidateWidthPS[target]) + 1,
          .latestPS = static_cast<std::int64_t>(box.upper()) - 1,
          .afterPS = box.lower(),
          .phiCentre = phiCentre,
          .phiHalfWidth = kPi<Scalar>,
          .secondCentre = secondCentre,
          .secondHalfWidth = follow::kOpenBound,
        };
        if (!follow::predictedWindow(propagator, parameters, geometry, detail::LinkSteps{}, state, factors, target, direction, prediction, window)) {
          continue;
        }
        frame.grid().visitCandidates(target, window, frame.grid().geometry().allBlocks(), [&](std::uint32_t row) {
          const auto measurement = frame.measurement(row);
          if (!spatiallyCompatible(window, measurement.phi, secondCoordinate<Det>(target, measurement.r, measurement.z))) {
            return;
          }
          const auto found = std::ranges::lower_bound(partners, std::pair{row, std::uint32_t{0}});
          if (found == partners.end() || found->first != row) {
            return;
          }
          const auto partner = found->second;
          const auto inner = inward ? partner : piece;
          const auto outer = inward ? piece : partner;
          const auto score = test(parametersOf, frame, geometry, propagator, view, inner, outer);
          if (score < best.score || (score == best.score && best.partner != kNone && partner < best.partner)) {
            best = {.partner = partner, .pieceInner = !inward, .score = score};
          }
        });
      }
    }
    return best;
  }

  /// The pair's score, infinity where it is refused.
  template <typename ParametersOf, typename View>
  [[nodiscard]] static Scalar test(const ParametersOf& parametersOf, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, const View& view, std::uint32_t inner, std::uint32_t outer)
  {
    constexpr Scalar kRefused = std::numeric_limits<Scalar>::infinity();
    const auto& measurements = frame.measurements();
    std::uint64_t innerLayers = 0U;
    std::uint64_t outerLayers = 0U;
    for (std::uint8_t i = 0; i < view.count[inner]; ++i) {
      innerLayers |= std::uint64_t{1} << geometry.surface(measurements.surface[view.clusters[inner].index[i]]).layer;
    }
    for (std::uint8_t i = 0; i < view.count[outer]; ++i) {
      outerLayers |= std::uint64_t{1} << geometry.surface(measurements.surface[view.clusters[outer].index[i]]).layer;
    }
    if ((innerLayers & outerLayers) != 0U || static_cast<std::size_t>(view.count[inner]) + view.count[outer] > kMaxClusters) {
      return kRefused;
    }
    const TimeEstPS innerBox{.start = view.timePS[inner], .width = view.widthPS[inner]};
    const TimeEstPS outerBox{.start = view.timePS[outer], .width = view.widthPS[outer]};
    if (!innerBox.isCompatible(outerBox) || (view.outwardState[inner].qOverP() > 0.F) != (view.state[outer].qOverP() > 0.F)) {
      return kRefused;
    }
    const auto& parameters = parametersOf(view.iteration[inner]);
    ClusterList<kMaxClusters> rows;
    std::uint8_t count = 0U;
    const auto fit = joinedFit(parameters, frame, geometry, propagator, view, inner, outer, rows, count);
    if (!follow::passesRoadChi2(parameters, fit)) {
      return kRefused;
    }
    const int outwardNdf = (2 * static_cast<int>(count)) - 5 - (2 * static_cast<int>(fit.outwardKinks));
    if (outwardNdf <= 0 || outwardNdf >= static_cast<int>(std::size(parameters.roadChi2MaxByNdf)) || !(fit.outwardChi2 <= parameters.roadChi2MaxByNdf[outwardNdf])) {
      return kRefused;
    }
    const int added = fit.ndf - view.ndf[inner] - view.ndf[outer];
    if (added <= 0) {
      return kRefused;
    }
    return (fit.chi2 - view.chi2[inner] - view.chi2[outer]) / static_cast<Scalar>(added);
  }

  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  BoundedVector<std::pair<std::uint32_t, std::uint32_t>> mEnds;   ///< (row, track): each earlier track's outermost row, by row
  BoundedVector<std::pair<std::uint32_t, std::uint32_t>> mStarts; ///< and its innermost
  BoundedVector<Choice> mChoices;                                 ///< per piece
  BoundedVector<std::uint32_t> mBest;                             ///< per earlier track, the winning piece
  BoundedVector<std::uint32_t> mJoined;
  BoundedVector<std::uint32_t> mChanged;
};

} // namespace o2::gtrk

#endif
