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
/// on past the layer (an inactive layer is no hole); a join puts two pieces
/// together across it. After the iteration, each of its tracks -- and each
/// earlier one a join changed -- that holds rows on both sides of an
/// inactive layer (a barrel layer: two consecutive rows' radii about its
/// own; a disc: a link from the one row's layer onto it and from it onto the
/// other's) and none on it looks for one: in
/// the windows of the iteration's gate from both sides, the forward pass of
/// its rows before the layer and the backward pass of those after it (each
/// side with three rows or more), the free rows passing the gate from
/// either; of them the one whose insertion raises the whole track's refit
/// chi2 least, kept where that rise -- the row's chi2 against the track
/// fitted on both sides -- is within the gate and the refit passes the
/// road's cut. A track whose walk ended short of a barrel inactive layer
/// looks the same way beyond its outermost row, in the window of its outward
/// state alone, where its walk can have ended only for lack of a sensor: the
/// inactive layer is the next one it crosses, or one active layer lies
/// between whose gate is not wholly on silicon (a gap the particle may have
/// crossed). A row two tracks want goes to the smaller rise (ties by the
/// track's index).
///
/// The per-track work is three steps, each an element functor the host's
/// tracker runs here and the device's in its own stage: `inactive::Find`
/// lists each track's candidate rows (a stage's records, in the order of the
/// search, after one for the track's fit without a row), `inactive::Score`
/// fits each, `inactive::Pick` takes each track's least rise in that order;
/// `inactive::Apply` writes the rows the arbitration keeps. The fits apart
/// from the search run one per thread on a device, not the few tracks of a
/// group that found a row.

#ifndef GTRK_TRACKING_INACTIVEROWS_H_
#define GTRK_TRACKING_INACTIVEROWS_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/KFTrackFitter.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/SensorCells.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackRefit.h"

#include "GPUCommonMath.h"

#include <memlayout/wrapper.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <type_traits>
#include <utility>

namespace o2::gtrk
{

namespace inactive
{

inline constexpr std::uint32_t kNone = std::numeric_limits<std::uint32_t>::max();

/// A track's row: the row, its place among the track's rows, the refit's
/// chi2 rise with it.
struct Choice {
  std::uint32_t row = kNone;
  std::uint32_t at = 0U;
  Scalar rise = 0.F;
};

/// A candidate's row (`Find`): the candidate (its place among the
/// candidates), the row and its place among the track's rows; the row
/// `kNone` for the track's fit without one, its candidate's first record.
struct RowCandidate {
  // As for TrackletT: upstream's `apply` macro also uses the member names as
  // enum constants.
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
  template <template <class> class C>
  struct RecordT {
    MEMLAYOUT_APPLY(RecordT, element, row, at)
    C<std::uint32_t> element;
    C<std::uint32_t> row;
    C<std::uint8_t> at;
  };
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
};

/// \a rows (\a count of them) refitted with a budget of the track's own
/// \a kinks.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() TrackFitResult refit(const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, const ClusterList<maxChainLayers<Det>>& rows, std::uint8_t count, std::uint8_t kinks)
{
  const KinkBudget budget{.max = static_cast<int>(kinks), .gate = parameters.kinkChi2Min};
  return fitMeasurements(propagator, geometry, detail::LinkSteps{}, RowChain<Det>{.frame = &frame, .rows = rows.index}, static_cast<int>(count), budget);
}

/// Whether a link of the layout leads from \a from onto \a to.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() bool linked(const IterationParameters<Det>& parameters, LayerID from, LayerID to)
{
  return ((parameters.layerLinksOut[from] >> to) & 1U) != 0U;
}

/// Whether \a state crosses a sensor of \a layer forward with every
/// crossing its gate admits on silicon (the four corners of the gate's
/// spread about the crossing, along the track's direction there).
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() bool certainCrossing(const IterationParameters<Det>& parameters, const TrackingGeometry<Det>& geometry, const Propagator& propagator, SensorCellsView cells, const BoundTrackState& state, LayerID layer)
{
  using Math = o2::gpu::CAMath;
  const auto direction = PropagationDirection::kForward;
  const auto& from = geometry.planeOf(state);
  const auto prediction = predictLayerCrossing(propagator, parameters, from, state, layer, direction);
  if (prediction.status != PropagationStatus::kOk) {
    return false;
  }
  CandidateWindow window{.phiHalfWidth = kPi<Scalar>, .secondHalfWidth = follow::kOpenBound};
  follow::windowCentre(parameters, layer, prediction, window.phiCentre, window.secondCentre);
  Scalar reach[2]{};
  if (!follow::predictedWindow(propagator, parameters, geometry, detail::LinkSteps{}, state, factorCovariance(state.cov), layer, direction, prediction, window, reach)) {
    return false;
  }
  Scalar point[3]{};
  follow::crossingPoint(parameters, layer, prediction, point);
  const auto& carrier = pointPropagator(propagator, from, state);
  const Scalar curvature = state.qOverP() * carrier.bzKiloGauss * static_cast<Scalar>(constants::math::B2C);
  Scalar landing[3]{};
  Scalar heading[3]{};
  const bool crossed = cells.anySensor(layer, prediction.phi, reach[0], prediction.second, reach[1], [&](std::uint32_t id) {
    const TrackingSurface plane = geometry.planeAt(static_cast<SensorID>(id), point);
    BoundParameters at = static_cast<const BoundParameters&>(state);
    const auto result = propagate<Scalar>(carrier, at, from, plane, direction);
    if (result.status != PropagationStatus::kOk || follow::beyondHalfTurn(curvature * result.pathLength)) {
      return false;
    }
    surfaceToGlobal(plane, at.loc0(), at.loc1(), landing);
    Scalar sinPhi = 0.F;
    Scalar cosPhi = 0.F;
    Scalar sinTheta = 0.F;
    Scalar cosTheta = 0.F;
    sincosOf(at.phi(), sinPhi, cosPhi);
    sincosOf(at.theta(), sinTheta, cosTheta);
    heading[0] = sinTheta * cosPhi;
    heading[1] = sinTheta * sinPhi;
    heading[2] = cosTheta;
    return true;
  });
  if (!crossed) {
    return false;
  }
  const Scalar phiSpread = Math::Max(window.phiHalfWidth - reach[0], 0.F);
  const Scalar secondSpread = Math::Max(window.secondHalfWidth - reach[1], 0.F);
  const Scalar radius = Math::Sqrt((landing[0] * landing[0]) + (landing[1] * landing[1]));
  const Scalar phi = Math::ATan2(landing[1], landing[0]);
  for (int corner = 0; corner < 4; ++corner) {
    const Scalar phiAt = phi + ((corner & 1) != 0 ? phiSpread : -phiSpread);
    const Scalar zAt = landing[2] + ((corner & 2) != 0 ? secondSpread : -secondSpread);
    Scalar sinPhi = 0.F;
    Scalar cosPhi = 0.F;
    sincosOf(phiAt, sinPhi, cosPhi);
    const Scalar at[3] = {radius * cosPhi, radius * sinPhi, zAt};
    const bool onSensor = cells.anySensor(layer, phiAt, reach[0], zAt, reach[1], [&](std::uint32_t id) {
      const TrackingSurface plane = geometry.planeAt(static_cast<SensorID>(id), at);
      Scalar step = 0.F;
      Scalar loc[2]{};
      return lineCrossing(plane, at, heading, step, loc) && plane.bounds.isInside(loc[0], loc[1]);
    });
    if (!onSensor) {
      return false;
    }
  }
  return true;
}

/// Per candidate \a element (the track \a candidates[element] of \a tracks),
/// its candidate rows on the first inactive layer it straddles without one,
/// or beyond its outermost row, emitted in the order of the search after the
/// record of its fit without one: the free rows in \a used, each layer's
/// candidates bound by \a maxCandidateWidthPS, the sensors' cells \a cells.
template <DetectorLayout Det, typename Tracks>
struct Find {
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using Record = typename BoundedColumns<RowCandidate::template RecordT>::Record;
  const IterationParameters<Det>* iteration = nullptr; ///< the iteration's
  const TimeFrame<Det>* timeFrame = nullptr;
  const TrackingGeometry<Det>* trackingGeometry = nullptr;
  Propagator propagator;
  Tracks tracks;
  const std::uint32_t* candidates = nullptr;
  const std::uint8_t* used = nullptr;
  const std::uint32_t* maxCandidateWidthPS = nullptr; ///< per layer
  SensorCellsView cells;

  template <typename Emitter>
  GPUhdi() void operator()(std::uint32_t element, Emitter& emitter) const
  {
    const auto track = candidates[element];
    const auto& parameters = *iteration;
    const auto& frame = *timeFrame;
    const auto& geometry = *trackingGeometry;
    const follow::StepContext<Det> context{.parameters = iteration, .frame = timeFrame, .geometry = trackingGeometry, .propagator = &propagator, .used = used};
    const auto n = static_cast<std::uint32_t>(tracks.count[track]);
    if (tracks.status[track] != PropagationStatus::kOk || n < 2U || n + 1U > kMaxClusters) {
      return;
    }
    const auto& measurements = frame.measurements();
    const auto& rows = tracks.clusters[track];
    const auto layerOf = [&](std::uint32_t row) { return geometry.surface(measurements.surface[row]).layer; };
    std::uint64_t held = 0U;
    for (std::uint32_t i = 0; i < n; ++i) {
      held |= std::uint64_t{1} << layerOf(rows.index[i]);
    }
    for (std::uint32_t i = 0; i + 1U < n; ++i) {
      const LayerID before = layerOf(rows.index[i]);
      const LayerID after = layerOf(rows.index[i + 1U]);
      for (std::uint64_t left = Detector<Det>::inactiveLayers & parameters.followLayers & ~held; left != 0U; left &= left - 1U) {
        const auto layer = static_cast<LayerID>(o2::gpu::CAMath::Ctz64(left));
        // Crossed between the two: a barrel layer's radius between theirs,
        // else a link from the one onto it and from it onto the other.
        if constexpr (Detector<Det>::hasBarrel) {
          if (Detector<Det>::isBarrelLayer(layer)) {
            const Scalar radius = parameters.layerSpans[layer].position;
            if (!(frame.measurement(rows.index[i]).r < radius && radius < frame.measurement(rows.index[i + 1U]).r)) {
              continue;
            }
          } else if (!linked(parameters, before, layer) || !linked(parameters, layer, after)) {
            continue;
          }
        } else if (!linked(parameters, before, layer) || !linked(parameters, layer, after)) {
          continue;
        }
        const RowChain<Det> chain{.frame = &frame, .rows = rows.index};
        const TimeEstPS box{.start = tracks.timePS[track], .width = tracks.widthPS[track]};
        // The candidates: free rows passing the gate from either side.
        std::uint32_t found[2 * kMaxClusters];
        std::uint32_t nFound = 0U;
        for (int pass = 0; pass < 2; ++pass) {
          const bool forward = pass == 0;
          const auto side = forward ? i + 1U : n - i - 1U;
          if (side < 3U) {
            continue;
          }
          const auto direction = forward ? PropagationDirection::kForward : PropagationDirection::kBackward;
          const RowChain<Det> part{.frame = &frame, .rows = forward ? rows.index : rows.index + i + 1U};
          const BoundTrackState start = forward ? makeChainSeed(propagator, geometry, chain, static_cast<int>(n)) : tracks.outwardState[track];
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
          follow::visitInline<false, false, false>(
            context, state, factors, 0.F, 0.F, layer, direction, window, traversed, 0U, WorkCounter::kFollowUpdates, SoloGroup{}, [](std::uint32_t /*row*/, StepOutcome /*outcome*/, Scalar /*chi2*/) {},
            [&](const auto& candidate) {
              bool known = false;
              for (std::uint32_t k = 0; k < nFound; ++k) {
                known = known || found[k] == candidate.row;
              }
              if (nFound < 2 * kMaxClusters && !known) {
                found[nFound++] = candidate.row;
              }
            });
        }
        if (nFound != 0U) {
          emitter.emit(Record{element, kNone, 0U});
          for (std::uint32_t k = 0; k < nFound; ++k) {
            emitter.emit(Record{element, found[k], static_cast<std::uint8_t>(i + 1U)});
          }
        }
        return;
      }
    }
    // Beyond the outermost row: the walk ended for lack of a sensor only.
    if constexpr (Detector<Det>::hasBarrel) {
      const LayerID lastLayer = geometry.surface(frame.measurements().surface[rows.index[n - 1U]]).layer;
      if (Detector<Det>::isBarrelLayer(lastLayer)) {
        const Scalar lastPosition = parameters.layerSpans[lastLayer].position;
        for (std::uint64_t left = Detector<Det>::inactiveLayers & parameters.followLayers & ~held; left != 0U; left &= left - 1U) {
          const auto layer = static_cast<LayerID>(o2::gpu::CAMath::Ctz64(left));
          if (!Detector<Det>::isBarrelLayer(layer) || !(parameters.layerSpans[layer].position > lastPosition)) {
            continue;
          }
          int between = 0;
          bool certain = false;
          for (int other = 0; other < static_cast<int>(Detector<Det>::nBarrelLayers); ++other) {
            const auto position = parameters.layerSpans[other].position;
            if (Detector<Det>::isInactiveLayer(static_cast<LayerID>(other)) || ((parameters.followLayers >> other) & 1U) == 0U || !(position > lastPosition) || !(position < parameters.layerSpans[layer].position)) {
              continue;
            }
            ++between;
            certain = certain || certainCrossing(parameters, geometry, propagator, cells, tracks.outwardState[track], static_cast<LayerID>(other));
          }
          if (between > 1 || certain) {
            continue;
          }
          const BoundTrackState& state = tracks.outwardState[track];
          const auto prediction = predictLayerCrossing(propagator, parameters, geometry.planeOf(state), state, layer, PropagationDirection::kForward);
          if (prediction.status != PropagationStatus::kOk) {
            continue;
          }
          const TimeEstPS box{.start = tracks.timePS[track], .width = tracks.widthPS[track]};
          Scalar phiCentre = 0.F;
          Scalar secondCentre = 0.F;
          follow::windowCentre(parameters, layer, prediction, phiCentre, secondCentre);
          CandidateWindow window{.earliestPS = static_cast<std::int64_t>(box.lower()) - static_cast<std::int64_t>(maxCandidateWidthPS[layer]) + 1, .latestPS = static_cast<std::int64_t>(box.upper()) - 1, .afterPS = box.lower(), .phiCentre = phiCentre, .phiHalfWidth = kPi<Scalar>, .secondCentre = secondCentre, .secondHalfWidth = follow::kOpenBound};
          const auto factors = factorCovariance(state.cov);
          detail::StepMaterial traversed;
          if (!follow::predictedWindow(propagator, parameters, geometry, detail::LinkSteps{}, state, factors, layer, PropagationDirection::kForward, prediction, window, nullptr, &traversed)) {
            continue;
          }
          bool listed = false;
          follow::visitInline<false, false, false>(
            context, state, factors, 0.F, 0.F, layer, PropagationDirection::kForward, window, traversed, 0U, WorkCounter::kFollowUpdates, SoloGroup{}, [](std::uint32_t /*row*/, StepOutcome /*outcome*/, Scalar /*chi2*/) {},
            [&](const auto& candidate) {
              if (!listed) {
                emitter.emit(Record{element, kNone, 0U});
                listed = true;
              }
              emitter.emit(Record{element, candidate.row, static_cast<std::uint8_t>(n)});
            });
          return;
        }
      }
    }
  }
};

/// Per record \a k of `Find` (\a elements, \a rows, \a at): the candidate's
/// track refitted with its own kinks, with the row inserted at its place, its
/// chi2 into \a chi2[k] where the fit passes the road's cut; for the record
/// without a row, the track's rows alone, its chi2 where the fit stands.
/// Infinity otherwise.
template <DetectorLayout Det, typename Tracks>
struct Score {
  const IterationParameters<Det>* iteration = nullptr; ///< the iteration's
  const TimeFrame<Det>* timeFrame = nullptr;
  const TrackingGeometry<Det>* trackingGeometry = nullptr;
  Propagator propagator;
  Tracks tracks;
  const std::uint32_t* candidates = nullptr;
  const std::uint32_t* elements = nullptr;
  const std::uint32_t* rows = nullptr;
  const std::uint8_t* at = nullptr;
  Scalar* chi2 = nullptr;

  GPUhdi() void operator()(std::uint32_t k) const
  {
    constexpr Scalar kRefused = std::numeric_limits<Scalar>::infinity();
    const auto track = candidates[elements[k]];
    const auto n = tracks.count[track];
    const auto& own = tracks.clusters[track];
    if (rows[k] == kNone) {
      const auto plain = refit(*iteration, *timeFrame, *trackingGeometry, propagator, own, n, tracks.kinks[track]);
      chi2[k] = plain.status == PropagationStatus::kOk ? plain.chi2 : kRefused;
      return;
    }
    ClusterList<maxChainLayers<Det>> with;
    for (std::uint32_t m = 0, j = 0; m <= n; ++m) {
      with.index[m] = m == at[k] ? rows[k] : own.index[j++];
    }
    const auto fit = refit(*iteration, *timeFrame, *trackingGeometry, propagator, with, static_cast<std::uint8_t>(n + 1U), tracks.kinks[track]);
    chi2[k] = follow::passesRoadChi2(*iteration, fit) ? fit.chi2 : kRefused;
  }
};

/// Per candidate \a element, of its records (\a offsets[element] up to
/// \a offsets[element + 1], `Find`'s, their \a chi2 from `Score`) the row
/// whose fit stands and whose rise over the fit without one is the least,
/// the first found among equals, kept within the gate, into
/// \a choices[element]; none where the fit without one fails.
template <DetectorLayout Det>
struct Pick {
  const IterationParameters<Det>* iteration = nullptr; ///< the iteration's
  const std::uint32_t* offsets = nullptr;
  const std::uint32_t* rows = nullptr;
  const std::uint8_t* at = nullptr;
  const Scalar* chi2 = nullptr;
  Choice* choices = nullptr;

  GPUhdi() void operator()(std::uint32_t element) const
  {
    Choice best;
    const auto first = offsets[element];
    const auto end = offsets[element + 1U];
    constexpr Scalar kRefused = std::numeric_limits<Scalar>::infinity();
    if (first != end && chi2[first] < kRefused) {
      const Scalar plain = chi2[first];
      for (auto k = first + 1U; k < end; ++k) {
        const Scalar rise = chi2[k] - plain;
        if (chi2[k] < kRefused && rise <= iteration->clusterChi2Max && (best.row == kNone || rise < best.rise)) {
          best = {.row = rows[k], .at = at[k], .rise = rise};
        }
      }
    }
    choices[element] = best;
  }
};

/// Per kept choice \a i (\a order[i], a candidate's place): its row
/// inserted into the track's, the track refitted with its own kinks and
/// written, the row marked exclusive in \a used.
template <DetectorLayout Det, typename Tracks>
struct Apply {
  const IterationParameters<Det>* iteration = nullptr; ///< the iteration's
  const TimeFrame<Det>* timeFrame = nullptr;
  const TrackingGeometry<Det>* trackingGeometry = nullptr;
  Propagator propagator;
  Tracks tracks;
  const std::uint32_t* candidates = nullptr;
  const Choice* choices = nullptr;
  const std::uint32_t* order = nullptr;
  std::uint8_t* used = nullptr;

  GPUhdi() void operator()(std::uint32_t i) const
  {
    const auto at = order[i];
    const auto track = candidates[at];
    const auto& choice = choices[at];
    ClusterList<maxChainLayers<Det>> rows;
    const auto n = tracks.count[track];
    for (std::uint32_t k = 0, j = 0; k <= n; ++k) {
      rows.index[k] = k == choice.at ? choice.row : tracks.clusters[track].index[j++];
    }
    const auto fit = refit(*iteration, *timeFrame, *trackingGeometry, propagator, rows, static_cast<std::uint8_t>(n + 1U), tracks.kinks[track]);
    writeFit(tracks, track, fit, rows, static_cast<std::uint8_t>(n + 1U), trackTimeBox(timeFrame->measurements(), rows, static_cast<std::uint8_t>(n + 1U)));
    used[choice.row] = kExclusive;
  }
};

} // namespace inactive

template <DetectorLayout Det>
class InactiveRows
{
 public:
  static constexpr std::size_t kLayerCount = static_cast<std::size_t>(Detector<Det>::nLayers);
  static constexpr std::uint32_t kNone = inactive::kNone;
  using Tracks = typename TrackRefit<Det>::Tracks;
  using Rows = StageOutput<inactive::RowCandidate::template RecordT>;

  /// On \a arena, the tracker's.
  explicit InactiveRows(detail::TaskArena& arena) : mStage(arena) {}

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool = std::move(pool);
    mStage.setMemoryPool(mMemoryPool);
    mRows.setMemoryPool(mMemoryPool);
    clearStorage();
  }
  void clearStorage() noexcept
  {
    mStage.releaseStorage();
    mRows.clearStorage();
    deepVectorClear(mMemoryPool.get(), mChi2, mChoices, mOrder, mChanged);
  }

  /// The search for \a tracks' \a candidates (track indices, ascending), with
  /// \a parameters (the iteration's), the \a used mask, each layer's
  /// candidates bound by \a maxCandidateWidthPS: the tracks that took a row
  /// rewritten and listed (`changed`, ascending), their rows marked used.
  /// \a forEach(n, body) runs body(i) for every i below n; \a chunk is the
  /// rows' stage's chunk.
  template <typename ForEach>
  [[nodiscard]] StageStatus run(ForEach forEach, std::uint32_t chunk, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, Tracks& tracks, std::span<const std::uint32_t> candidates, std::span<std::uint8_t> used, std::span<const std::uint32_t, kLayerCount> maxCandidateWidthPS, SensorCellsView cells)
  {
    if constexpr (Detector<Det>::inactiveLayers == 0U) {
      resizeBoundedVector(mChanged, 0U, mMemoryPool.get());
      return StageStatus::kComplete;
    } else {
      return guardOutOfMemory([&] {
        const auto view = tracks.view();
        using View = std::remove_const_t<decltype(view)>;
        const auto count = static_cast<std::uint32_t>(candidates.size());
        resizeBoundedVector(mChanged, 0U, mMemoryPool.get());
        if (count == 0U) {
          return StageStatus::kComplete;
        }
        // Each candidate's rows, their fits, its row of the least rise.
        mRows.prepare(count);
        const auto status = mStage.run(mRows, chunk, inactive::Find<Det, View>{.iteration = &parameters, .timeFrame = &frame, .trackingGeometry = &geometry, .propagator = propagator, .tracks = view, .candidates = candidates.data(), .used = used.data(), .maxCandidateWidthPS = maxCandidateWidthPS.data(), .cells = cells});
        if (status != StageStatus::kComplete) {
          return status;
        }
        const auto found = std::as_const(mRows).view();
        const auto total = mRows.size();
        resizeBoundedVector(mChi2, total, mMemoryPool.get());
        forEach(total, inactive::Score<Det, View>{.iteration = &parameters, .timeFrame = &frame, .trackingGeometry = &geometry, .propagator = propagator, .tracks = view, .candidates = candidates.data(), .elements = found.element, .rows = found.row, .at = found.at, .chi2 = mChi2.data()});
        resizeBoundedVector(mChoices, count, mMemoryPool.get());
        forEach(count, inactive::Pick<Det>{.iteration = &parameters, .offsets = mRows.offsets().data(), .rows = found.row, .at = found.at, .chi2 = mChi2.data(), .choices = mChoices.data()});
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
        forEach(kept, inactive::Apply<Det, View>{.iteration = &parameters, .timeFrame = &frame, .trackingGeometry = &geometry, .propagator = propagator, .tracks = view, .candidates = candidates.data(), .choices = mChoices.data(), .order = mOrder.data(), .used = used.data()});
        resizeBoundedVector(mChanged, kept, mMemoryPool.get());
        for (std::uint32_t i = 0; i < kept; ++i) {
          mChanged[i] = candidates[mOrder[i]];
        }
        std::ranges::sort(mChanged);
        return StageStatus::kComplete;
      });
    }
  }

  /// The last run's tracks that took a row, ascending.
  [[nodiscard]] std::span<const std::uint32_t> changed() const { return mChanged; }

 private:
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  CpuStage<inactive::RowCandidate::template RecordT> mStage;
  Rows mRows;                  ///< each candidate's rows
  BoundedVector<Scalar> mChi2; ///< per row, `inactive::Score`
  BoundedVector<inactive::Choice> mChoices;
  BoundedVector<std::uint32_t> mOrder;
  BoundedVector<std::uint32_t> mChanged;
};

} // namespace o2::gtrk

#endif
