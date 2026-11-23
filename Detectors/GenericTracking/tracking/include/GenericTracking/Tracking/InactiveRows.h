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

/// \file InactiveRows.h
/// \brief An inactive layer's row (ALICE 3's TOF) for a track that crosses
///        the layer between two of its rows without one, matched on the
///        whole track's fit.
///
/// The follow takes an inactive layer's row from the prediction of the
/// rows before it alone; where a scatter at the layer before bends the
/// track, the particle's row lies outside that prediction, and the walk goes
/// on past the layer (an inactive layer is no hole). After the iteration,
/// each of its tracks -- and each earlier one a join changed -- that holds
/// rows on both sides of an inactive layer and none on it looks for one: in
/// the windows of the iteration's gate from both sides, the forward pass of
/// its rows before the layer and the backward pass of those after it (each
/// side with three rows or more), the free rows passing the gate from
/// either; of them the one whose insertion raises the whole track's refit
/// chi2 least, kept where that rise -- the row's chi2 against the track
/// fitted on both sides -- is within the gate and the refit passes the
/// road's cut. A row two tracks want goes to the smaller rise (ties by the
/// track's index).

#ifndef GTRK_TRACKING_INACTIVEROWS_H_
#define GTRK_TRACKING_INACTIVEROWS_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/KFTrackFitter.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackRefit.h"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <utility>

namespace o2::gtrk
{

template <DetectorLayout Det>
class InactiveRows
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
  void clearStorage() noexcept { deepVectorClear(mMemoryPool.get(), mChoices, mOrder, mChanged); }

  /// The search for \a tracks' \a candidates (track indices, ascending), with
  /// \a parameters (the iteration's), the \a used mask, each layer's
  /// candidates bound by \a maxCandidateWidthPS: the tracks that took a row
  /// rewritten and listed (`changed`, ascending), their rows marked used.
  template <typename ForEach>
  [[nodiscard]] StageStatus run(ForEach forEach, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, Tracks& tracks, std::span<const std::uint32_t> candidates, std::span<std::uint8_t> used, std::span<const std::uint32_t, kLayerCount> maxCandidateWidthPS)
  {
    if constexpr (Detector<Det>::inactiveLayers == 0U) {
      resizeBoundedVector(mChanged, 0U, mMemoryPool.get());
      return StageStatus::kComplete;
    } else {
      return guardOutOfMemory([&] {
        const auto view = tracks.view();
        const auto count = static_cast<std::uint32_t>(candidates.size());
        resizeBoundedVector(mChoices, count, mMemoryPool.get());
        const follow::StepContext<Det> context{.parameters = &parameters, .frame = &frame, .geometry = &geometry, .propagator = &propagator, .used = used.data()};
        forEach(count, [&](std::uint32_t at) { mChoices[at] = choose(context, view, candidates[at], maxCandidateWidthPS); });
        // A row two tracks want to the smaller rise: by row, then rise, then track.
        std::uint32_t chosen = 0U;
        for (const auto& choice : mChoices) {
          chosen += choice.row != kNone ? 1U : 0U;
        }
        resizeBoundedVector(mOrder, chosen, mMemoryPool.get());
        for (std::uint32_t at = 0, i = 0; at < count; ++at) {
          if (mChoices[at].row != kNone) {
            mOrder[i++] = at;
          }
        }
        std::ranges::sort(mOrder, [&](std::uint32_t a, std::uint32_t b) {
          const auto& x = mChoices[a];
          const auto& y = mChoices[b];
          return x.row != y.row ? x.row < y.row : (x.rise != y.rise ? x.rise < y.rise : candidates[a] < candidates[b]);
        });
        std::uint32_t kept = 0U;
        for (std::uint32_t i = 0; i < chosen; ++i) {
          if (i == 0U || mChoices[mOrder[i]].row != mChoices[mOrder[i - 1U]].row) {
            mOrder[kept++] = mOrder[i];
          }
        }
        resizeBoundedVector(mChanged, kept, mMemoryPool.get());
        forEach(kept, [&](std::uint32_t i) {
          const auto at = mOrder[i];
          const auto track = candidates[at];
          const auto& choice = mChoices[at];
          ClusterList<kMaxClusters> rows;
          const auto n = view.count[track];
          for (std::uint32_t k = 0, j = 0; k <= n; ++k) {
            rows.index[k] = k == choice.at ? choice.row : view.clusters[track].index[j++];
          }
          const auto fit = refit(parameters, frame, geometry, propagator, rows, static_cast<std::uint8_t>(n + 1U), view.kinks[track]);
          writeFit(view, track, fit, rows, static_cast<std::uint8_t>(n + 1U), trackTimeBox(frame.measurements(), rows, static_cast<std::uint8_t>(n + 1U)));
          used[choice.row] = kExclusive;
          mChanged[i] = track;
        });
        std::ranges::sort(mChanged);
        return StageStatus::kComplete;
      });
    }
  }

  /// The last run's tracks that took a row, ascending.
  [[nodiscard]] std::span<const std::uint32_t> changed() const { return mChanged; }

 private:
  /// A track's row: the row, its place among the track's rows, the refit's
  /// chi2 rise with it.
  struct Choice {
    std::uint32_t row = kNone;
    std::uint32_t at = 0U;
    Scalar rise = 0.F;
  };

  /// \a rows (\a count of them) refitted with the iteration's crossing rule
  /// and a budget of the track's own \a kinks.
  [[nodiscard]] static TrackFitResult refit(const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, const ClusterList<kMaxClusters>& rows, std::uint8_t count, std::uint8_t kinks)
  {
    const KinkBudget budget{.max = static_cast<int>(kinks), .gate = parameters.kinkChi2Min};
    return withCrossingRule(parameters.crossingRule, [&]<CrossingRule Rule>(std::integral_constant<CrossingRule, Rule>) {
      return fitMeasurements<Rule>(propagator, geometry, detail::LinkSteps{}, RowChain<Det>{.frame = &frame, .rows = rows.index}, static_cast<int>(count), budget);
    });
  }

  [[nodiscard]] static bool linked(LayerID from, LayerID to)
  {
    return std::ranges::any_of(Det::layerLinks, [&](const auto& link) { return link.from == from && link.to == to; });
  }

  /// \a track's row on the first inactive layer it straddles without one.
  template <typename View>
  [[nodiscard]] Choice choose(const follow::StepContext<Det>& context, const View& view, std::uint32_t track, std::span<const std::uint32_t, kLayerCount> maxCandidateWidthPS) const
  {
    Choice best;
    const auto& parameters = *context.parameters;
    const auto& frame = *context.frame;
    const auto& geometry = *context.geometry;
    const auto& propagator = *context.propagator;
    const auto n = static_cast<std::uint32_t>(view.count[track]);
    if (view.status[track] != PropagationStatus::kOk || n < 2U || n + 1U > kMaxClusters) {
      return best;
    }
    const auto& measurements = frame.measurements();
    const auto& rows = view.clusters[track];
    const auto layerOf = [&](std::uint32_t row) { return geometry.surface(measurements.surface[row]).layer; };
    std::uint64_t held = 0U;
    for (std::uint32_t i = 0; i < n; ++i) {
      held |= std::uint64_t{1} << layerOf(rows.index[i]);
    }
    for (std::uint32_t i = 0; i + 1U < n; ++i) {
      const LayerID before = layerOf(rows.index[i]);
      const LayerID after = layerOf(rows.index[i + 1U]);
      for (std::uint64_t left = Detector<Det>::inactiveLayers & parameters.followLayers & ~held; left != 0U; left &= left - 1U) {
        const auto layer = static_cast<LayerID>(std::countr_zero(left));
        if (!linked(before, layer) || !linked(layer, after)) {
          continue;
        }
        const RowChain<Det> chain{.frame = &frame, .rows = rows.index};
        const TimeEstPS box{.start = view.timePS[track], .width = view.widthPS[track]};
        const auto plain = refit(parameters, frame, geometry, propagator, rows, static_cast<std::uint8_t>(n), view.kinks[track]);
        if (plain.status != PropagationStatus::kOk) {
          return best;
        }
        // The candidates: free rows passing the gate from either side.
        std::uint32_t found[2 * kMaxClusters];
        std::uint32_t nFound = 0U;
        for (const bool forward : {true, false}) {
          const auto side = forward ? i + 1U : n - i - 1U;
          if (side < 3U) {
            continue;
          }
          const auto direction = forward ? PropagationDirection::kForward : PropagationDirection::kBackward;
          const RowChain<Det> part{.frame = &frame, .rows = forward ? rows.index : rows.index + i + 1U};
          const BoundTrackState start = forward ? makeChainSeed(propagator, geometry, chain, static_cast<int>(n)) : view.outwardState[track];
          BoundTrackState state;
          Scalar chi2 = 0.F;
          std::int16_t failed = -1;
          std::int16_t kinkAt = -1;
          int kinks = 0;
          const auto span = detail::spanMaterial(geometry, detail::LinkSteps{}, part, static_cast<int>(side));
          if (detail::runPass(propagator, geometry, span, part, static_cast<int>(side), direction, start, state, chi2, failed, KinkBudget{}, kinks, kinkAt, NoFitObserver{}) != PropagationStatus::kOk) {
            continue;
          }
          const auto prediction = predictLayerCrossing(propagator, parameters, geometry.planeOf(state), state, layer, direction);
          if (prediction.status != PropagationStatus::kOk) {
            continue;
          }
          Scalar phiCentre = 0.F;
          Scalar secondCentre = 0.F;
          follow::windowCentre(parameters, layer, prediction, phiCentre, secondCentre);
          CandidateWindow window{
            .earliestPS = static_cast<std::int64_t>(box.lower()) - static_cast<std::int64_t>(maxCandidateWidthPS[layer]) + 1,
            .latestPS = static_cast<std::int64_t>(box.upper()) - 1,
            .afterPS = box.lower(),
            .phiCentre = phiCentre,
            .phiHalfWidth = kPi<Scalar>,
            .secondCentre = secondCentre,
            .secondHalfWidth = follow::kOpenBound,
          };
          const auto factors = factorCovariance(state.cov);
          detail::StepMaterial traversed;
          if (!follow::predictedWindow(propagator, parameters, geometry, detail::LinkSteps{}, state, factors, layer, direction, prediction, window, nullptr, &traversed)) {
            continue;
          }
          follow::visitInline<false, false, false, false>(
            context, state, factors, 0.F, 0.F, layer, direction, window, traversed, 0U, WorkCounter::kFollowUpdates, SoloGroup{}, [](std::uint32_t /*row*/, StepOutcome /*outcome*/, Scalar /*chi2*/) {},
            [&](const auto& candidate) {
              if (nFound < 2 * kMaxClusters && std::find(found, found + nFound, candidate.row) == found + nFound) {
                found[nFound++] = candidate.row;
              }
            });
        }
        // Of them the smallest rise of the whole track's chi2, within the gate.
        ClusterList<kMaxClusters> with;
        for (std::uint32_t k = 0; k < nFound; ++k) {
          for (std::uint32_t m = 0, j = 0; m <= n; ++m) {
            with.index[m] = m == i + 1U ? found[k] : rows.index[j++];
          }
          const auto fit = refit(parameters, frame, geometry, propagator, with, static_cast<std::uint8_t>(n + 1U), view.kinks[track]);
          if (!follow::passesRoadChi2(parameters, fit)) {
            continue;
          }
          const Scalar rise = fit.chi2 - plain.chi2;
          if (rise <= parameters.clusterChi2Max && (best.row == kNone || rise < best.rise)) {
            best = {.row = found[k], .at = i + 1U, .rise = rise};
          }
        }
        return best;
      }
    }
    return best;
  }

  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  BoundedVector<Choice> mChoices;
  BoundedVector<std::uint32_t> mOrder;
  BoundedVector<std::uint32_t> mChanged;
};

} // namespace o2::gtrk

#endif
