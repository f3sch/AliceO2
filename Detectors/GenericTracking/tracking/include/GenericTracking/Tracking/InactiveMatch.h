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

/// \file InactiveMatch.h
/// \brief After the last iteration, an inactive layer's row (ALICE 3's TOF)
///        matched to a track that ends before the layer, where the row is
///        more likely the track's than a chance one.
///
/// A track's walk ends at the first layer whose own row it refuses; its
/// inactive layers' rows further out stay free. After the last iteration,
/// each track predicts from its outward state every inactive layer beyond it
/// that it holds no row on; the free rows of that layer in the window of its
/// own iteration's gate, within its time box, are its candidates. A
/// candidate's likelihood ratio compares the density of the track's
/// prediction at the row -- exp(-chi2 / 2) / (2 pi sqrt(det S)), S the
/// update's innovation covariance -- with the density of rows there by
/// chance: the free rows of the layer on the candidate's sensor within the
/// track's time box over the sensor's area. A candidate with a ratio above
/// one is more likely the track's own row than a chance row; the track
/// takes the one of largest ratio per layer, and a row two tracks want goes
/// to the larger ratio (ties by the track's index).
///
/// The match does not enter the track's fit: the store holds it as a
/// reference of its own role (`ReferenceRole::kMatch`) after the track's
/// rows, and the row is marked used.
///
/// The per-track work is the element functor `match::Choose`; the
/// arbitration reads its outputs.

#ifndef GTRK_TRACKING_INACTIVEMATCH_H_
#define GTRK_TRACKING_INACTIVEMATCH_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/TimeEstPS.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/SensorCells.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackRefit.h"

#include "GPUCommonMath.h"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <type_traits>
#include <utility>

namespace o2::gtrk
{

namespace match
{

inline constexpr std::uint32_t kNone = std::numeric_limits<std::uint32_t>::max();

/// How many inactive layers \a Det has: a track's choices per layer.
template <DetectorLayout Det>
inline constexpr auto kInactiveCount = static_cast<std::uint32_t>(std::popcount(Detector<Det>::inactiveLayers));

/// A track's choice on one inactive layer: the row, its likelihood ratio, and
/// the squared distance from the track's outward state to the layer's
/// crossing (the order of the matches along the track).
struct Choice {
  std::uint32_t row = kNone;
  Scalar ratio = 0.F;
  Scalar distance2 = 0.F;
};

/// A track's matched row, as the store takes them: by track, then along the
/// track.
struct Matched {
  std::uint32_t track = 0U;
  std::uint32_t row = 0U;
};

/// The free rows of the inactive layers per sensor, in box-start order:
/// sensor s's from `begin[s]` to `begin[s + 1]` of `rows`, their box starts
/// in `startPS`.
struct SensorRows {
  const std::uint32_t* begin = nullptr;
  const std::uint32_t* rows = nullptr;
  const std::uint32_t* startPS = nullptr;
  std::uint32_t sensors = 0U;
};

/// How many of \a sensor's free rows in \a table have a box compatible with
/// \a box, theirs at most \a maxWidthPS wide.
template <typename Columns>
[[nodiscard]] GPUhdi() std::uint32_t inTime(const SensorRows& table, const Columns& measurements, std::uint32_t sensor, const TimeEstPS& box, std::uint32_t maxWidthPS)
{
  if (sensor >= table.sensors) {
    return 0U;
  }
  // The first start that can reach the box: after box.lower() - maxWidthPS.
  const auto earliest = static_cast<std::int64_t>(box.lower()) - static_cast<std::int64_t>(maxWidthPS);
  std::uint32_t lo = table.begin[sensor];
  std::uint32_t hi = table.begin[sensor + 1U];
  while (lo < hi) {
    const std::uint32_t mid = lo + ((hi - lo) / 2U);
    const std::int64_t start = table.startPS[mid];
    if (start <= earliest) {
      lo = mid + 1U;
    } else {
      hi = mid;
    }
  }
  std::uint32_t count = 0U;
  for (std::uint32_t at = lo; at < table.begin[sensor + 1U] && table.startPS[at] < box.upper(); ++at) {
    count += box.isCompatible(timeBoxAt(measurements, table.rows[at])) ? 1U : 0U;
  }
  return count;
}

/// Per track: on every inactive layer beyond it that it holds no row on, the
/// candidate of largest likelihood ratio above one, into \a choices
/// (`kInactiveCount` per track, in the layers' order).
template <DetectorLayout Det, typename Tracks>
struct Choose {
  static constexpr std::uint32_t kMaxClusters = maxChainLayers<Det>;
  const IterationParameters<Det>* iterations = nullptr; ///< every iteration's, the track's own read
  const TimeFrame<Det>* timeFrame = nullptr;
  const TrackingGeometry<Det>* trackingGeometry = nullptr;
  Propagator propagator;
  Tracks tracks;
  const std::uint8_t* used = nullptr;
  const std::uint32_t* maxCandidateWidthPS = nullptr;
  SensorRows freeRows;
  Choice* choices = nullptr;

  GPUhdi() void operator()(std::uint32_t track) const
  {
    Choice* mine = choices + (track * kInactiveCount<Det>);
    for (std::uint32_t k = 0; k < kInactiveCount<Det>; ++k) {
      mine[k] = Choice{};
    }
    const auto n = static_cast<std::uint32_t>(tracks.count[track]);
    if (tracks.status[track] != PropagationStatus::kOk || n == 0U) {
      return;
    }
    const auto& parameters = iterations[tracks.iteration[track]];
    const auto& frame = *timeFrame;
    const auto& geometry = *trackingGeometry;
    const auto& measurements = frame.measurements();
    std::uint64_t held = 0U;
    for (std::uint32_t i = 0; i < n; ++i) {
      held |= std::uint64_t{1} << geometry.surface(measurements.surface[tracks.clusters[track].index[i]]).layer;
    }
    const follow::StepContext<Det> context{.parameters = &parameters, .frame = timeFrame, .geometry = trackingGeometry, .propagator = &propagator, .used = used};
    const BoundTrackState& state = tracks.outwardState[track];
    const auto& from = geometry.planeOf(state);
    Scalar here[3]{};
    surfaceToGlobal(from, state.loc0(), state.loc1(), here);
    const TimeEstPS box{.start = tracks.timePS[track], .width = tracks.widthPS[track]};
    const auto factors = factorCovariance(state.cov);
    std::uint32_t slot = 0U;
    for (std::uint64_t left = Detector<Det>::inactiveLayers; left != 0U; left &= left - 1U, ++slot) {
      const auto layer = static_cast<LayerID>(o2::gpu::CAMath::Ctz64(left));
      if (((held >> layer) & 1U) != 0U) {
        continue;
      }
      const auto prediction = predictLayerCrossing(propagator, parameters, from, state, layer, PropagationDirection::kForward);
      if (prediction.status != PropagationStatus::kOk) {
        continue;
      }
      Scalar phiCentre = 0.F;
      Scalar secondCentre = 0.F;
      follow::windowCentre(parameters, layer, prediction, phiCentre, secondCentre);
      CandidateWindow window{.earliestPS = static_cast<std::int64_t>(box.lower()) - static_cast<std::int64_t>(maxCandidateWidthPS[layer]) + 1, .latestPS = static_cast<std::int64_t>(box.upper()) - 1, .afterPS = box.lower(), .phiCentre = phiCentre, .phiHalfWidth = kPi<Scalar>, .secondCentre = secondCentre, .secondHalfWidth = follow::kOpenBound};
      detail::StepMaterial traversed;
      if (!follow::predictedWindow(propagator, parameters, geometry, detail::LinkSteps{}, state, factors, layer, PropagationDirection::kForward, prediction, window, nullptr, &traversed)) {
        continue;
      }
      Choice best;
      follow::visitInline<true, false, false>(
        context, state, factors, 0.F, 0.F, layer, PropagationDirection::kForward, window, traversed, 0U, WorkCounter::kFollowUpdates, SoloGroup{}, [](std::uint32_t /*row*/, StepOutcome /*outcome*/, Scalar /*chi2*/) {},
        [&](const follow::DensityChild& candidate) {
          if (candidate.row == ClusterList<kMaxClusters>::kNoCluster || !(candidate.determinant > 0.F)) {
            return;
          }
          const auto sensor = measurements.surface[candidate.row];
          const auto& bounds = geometry.surface(sensor).bounds;
          // Itself among them: at least one.
          const auto chance = o2::gpu::CAMath::Max(inTime(freeRows, measurements, sensor, box, maxCandidateWidthPS[layer]), 1U);
          const Scalar density = static_cast<Scalar>(chance) / (4.F * bounds.halfX * bounds.halfY);
          const Scalar ratio = o2::gpu::CAMath::Exp(-0.5F * candidate.chi2) / (2.F * kPi<Scalar> * o2::gpu::CAMath::Sqrt(candidate.determinant) * density);
          if (ratio > 1.F && (best.row == kNone || ratio > best.ratio)) {
            best = {.row = candidate.row, .ratio = ratio};
          }
        });
      if (best.row == kNone) {
        continue;
      }
      Scalar point[3]{};
      follow::crossingPoint(parameters, layer, prediction, point);
      best.distance2 = ((point[0] - here[0]) * (point[0] - here[0])) + ((point[1] - here[1]) * (point[1] - here[1])) + ((point[2] - here[2]) * (point[2] - here[2]));
      mine[slot] = best;
    }
  }
};

} // namespace match

/// The match of the inactive layers' rows to the tracks after the last
/// iteration (see the file's comment).
template <DetectorLayout Det>
class InactiveMatch
{
 public:
  static constexpr auto kLayerCount = static_cast<std::uint32_t>(Detector<Det>::nLayers);
  static constexpr std::uint32_t kNone = match::kNone;
  static constexpr std::uint32_t kSlots = match::kInactiveCount<Det>;
  using Tracks = typename TrackRefit<Det>::Tracks;

  using Matched = match::Matched;

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool = std::move(pool);
    clearStorage();
  }
  void clearStorage() noexcept { deepVectorClear(mMemoryPool.get(), mBegin, mRows, mStarts, mChoices, mOrder, mMatched); }

  /// The match for every track of \a tracks, each with its own iteration's
  /// parameters (\a iterations), the \a used mask, each layer's candidates
  /// bound by \a maxCandidateWidthPS: the matches (`matched`), their rows
  /// marked used.
  template <typename ForEach>
  [[nodiscard]] StageStatus run(ForEach forEach, const IterationParameters<Det>* iterations, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, const Tracks& tracks, std::span<std::uint8_t> used, std::span<const std::uint32_t, kLayerCount> maxCandidateWidthPS)
  {
    if constexpr (Detector<Det>::inactiveLayers == 0U) {
      resizeBoundedVector(mMatched, 0U, mMemoryPool.get());
      return StageStatus::kComplete;
    } else {
      return guardOutOfMemory([&] {
        const auto& measurements = frame.measurements();
        const auto rowCount = static_cast<std::uint32_t>(used.size());
        // The free rows of the inactive layers per sensor, by box start.
        constexpr auto sensors = static_cast<std::uint32_t>(Det::nSensors);
        std::uint32_t freeRows = 0U;
        for (std::uint32_t row = 0; row < rowCount; ++row) {
          freeRows += used[row] == kFree && Detector<Det>::isInactiveLayer(geometry.surface(measurements.surface[row]).layer) ? 1U : 0U;
        }
        resizeBoundedVector(mBegin, std::size_t{sensors} + 1U, mMemoryPool.get());
        std::ranges::fill(mBegin, 0U);
        resizeBoundedVector(mRows, freeRows, mMemoryPool.get());
        resizeBoundedVector(mStarts, freeRows, mMemoryPool.get());
        for (std::uint32_t row = 0, at = 0; row < rowCount; ++row) {
          if (used[row] == kFree && Detector<Det>::isInactiveLayer(geometry.surface(measurements.surface[row]).layer)) {
            ++mBegin[measurements.surface[row] + 1U];
            mRows[at++] = row;
          }
        }
        for (std::uint32_t s = 0; s < sensors; ++s) {
          mBegin[s + 1U] += mBegin[s];
        }
        std::stable_sort(mRows.begin(), mRows.end(), [&](std::uint32_t a, std::uint32_t b) {
          const auto sa = measurements.surface[a];
          const auto sb = measurements.surface[b];
          return sa != sb ? sa < sb : measurements.timePS[a] < measurements.timePS[b];
        });
        for (std::uint32_t at = 0; at < freeRows; ++at) {
          mStarts[at] = measurements.timePS[mRows[at]];
        }
        const match::SensorRows table{.begin = mBegin.data(), .rows = mRows.data(), .startPS = mStarts.data(), .sensors = sensors};
        // Each track's choices.
        const auto count = tracks.size();
        resizeBoundedVector(mChoices, std::size_t{count} * kSlots, mMemoryPool.get());
        const auto view = tracks.view();
        forEach(count, match::Choose<Det, std::remove_const_t<decltype(view)>>{.iterations = iterations, .timeFrame = &frame, .trackingGeometry = &geometry, .propagator = propagator, .tracks = view, .used = used.data(), .maxCandidateWidthPS = maxCandidateWidthPS.data(), .freeRows = table, .choices = mChoices.data()});
        // A row two tracks want to the larger ratio: by row, then ratio
        // (larger first), then track.
        std::uint32_t chosen = 0U;
        for (const auto& choice : mChoices) {
          chosen += choice.row != kNone ? 1U : 0U;
        }
        resizeBoundedVector(mOrder, chosen, mMemoryPool.get());
        for (std::uint32_t at = 0, i = 0; at < mChoices.size(); ++at) {
          if (mChoices[at].row != kNone) {
            mOrder[i++] = at;
          }
        }
        std::ranges::sort(mOrder, [&](std::uint32_t a, std::uint32_t b) {
          const auto& ca = mChoices[a];
          const auto& cb = mChoices[b];
          if (ca.row != cb.row) {
            return ca.row < cb.row;
          }
          if (!exactlyEqual(ca.ratio, cb.ratio)) {
            return ca.ratio > cb.ratio;
          }
          return a < b;
        });
        std::uint32_t kept = 0U;
        for (std::uint32_t i = 0; i < mOrder.size(); ++i) {
          if (i == 0U || mChoices[mOrder[i]].row != mChoices[mOrder[i - 1U]].row) {
            mOrder[kept++] = mOrder[i];
          }
        }
        resizeBoundedVector(mOrder, kept, mMemoryPool.get());
        // By track, then along it.
        std::ranges::sort(mOrder, [&](std::uint32_t a, std::uint32_t b) {
          const auto ta = a / kSlots;
          const auto tb = b / kSlots;
          return ta != tb ? ta < tb : mChoices[a].distance2 < mChoices[b].distance2;
        });
        resizeBoundedVector(mMatched, kept, mMemoryPool.get());
        for (std::uint32_t i = 0; i < kept; ++i) {
          const auto& choice = mChoices[mOrder[i]];
          mMatched[i] = {.track = mOrder[i] / kSlots, .row = choice.row};
          used[choice.row] = kExclusive;
        }
        return StageStatus::kComplete;
      });
    }
  }

  /// The last run's matches, by track, then along the track.
  [[nodiscard]] std::span<const Matched> matched() const { return mMatched; }

 private:
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  BoundedVector<std::uint32_t> mBegin;  ///< per sensor, its first free row in `mRows`
  BoundedVector<std::uint32_t> mRows;   ///< the inactive layers' free rows, by sensor and box start
  BoundedVector<std::uint32_t> mStarts; ///< their box starts
  BoundedVector<match::Choice> mChoices;
  BoundedVector<std::uint32_t> mOrder;
  BoundedVector<Matched> mMatched;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_INACTIVEMATCH_H_
