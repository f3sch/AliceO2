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

/// \file OverlapPartners.h
/// \brief A track's second row on a layer: the measurement of its crossing
///        of an overlapping sensor.
///
/// Neighbouring sensors of a layer overlap, so a particle can leave a row on
/// each. The tracker takes one row per layer; an iteration with
/// `overlapPartners` attaches the others after its follow and length cut, as
/// *overlap partners*:
///
///  - each kept track is fitted on its rows, with the link material (its
///    rows are one per layer, which the link material's steps hold), and each
///    row's filtered state kept, the inward pass's where it has one
///    (`ObservesFitStates`);
///  - from each row's state, the other sensors of its layer the track
///    crosses -- those whose footprint meets the range the track sweeps
///    within the layer's envelope (`partnerReach`), and crossed within their
///    bounds widened by the gate's spread -- and on each
///    its free rows, compatible in time, within that spread of the crossing,
///    each predicted on the sensor and tested with the follow's
///    two-coordinate cut (`clusterChi2Max`); per sensor the best is the
///    track's proposal;
///  - a row proposed by more than one track goes to the lowest chi2, then the
///    lowest track (a deterministic order);
///  - each track with partners is fitted again on its rows and partners in
///    path order, with the link material between layers and no material of
///    its own on the step between a layer's two sensors (`LayerMaterial`: the
///    bands hold the layer's whole material); where that fit fails or the rows' time
///    boxes leave none in common, the track stays as it was, without them;
///  - the attached partners are taken in the used mask, exclusively.
///
/// Only the layers the layout says can have them are searched
/// (`overlapLayers`, LayerLinks.h), and a track with no row on one is passed
/// over.
///
/// The tracks' own rows stay one per layer (`TrackOf`); the partners go to
/// the store as references (`kOverlapPartner`, TrackStore.h), in path order
/// beside the rows. A track's states, chi2 and time are then those of the
/// fit with its partners.

#ifndef GTRK_TRACKING_OVERLAPPARTNERS_H_
#define GTRK_TRACKING_OVERLAPPARTNERS_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Base/LayerLinks.h"
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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <memory_resource>
#include <numbers>
#include <span>
#include <vector>

namespace o2::gtrk
{

/// One overlap partner of a track: its row, the position among the track's
/// own rows of its row on the same layer, and whether it lies after that
/// row along the track (in its rows' order) or before.
struct OverlapPartner {
  std::uint32_t row = 0U;
  std::uint8_t layerRow = 0U;
  std::uint8_t ahead = 0U;
};

/// The device's overlap partners kept beside its tracks (DeviceTracker): a
/// count per track, and the partners as `OverlapPartner`s.
struct PartnerCount {
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
  template <template <class> class C>
  struct RecordT {
    MEMLAYOUT_APPLY(RecordT, count)
    C<std::uint32_t> count;
  };
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
};

struct PartnerRow {
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
  template <template <class> class C>
  struct RecordT {
    MEMLAYOUT_APPLY(RecordT, row, layerRow, ahead)
    C<std::uint32_t> row;
    C<std::uint8_t> layerRow;
    C<std::uint8_t> ahead;
  };
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
};

/// A track's proposal of an overlap partner (`ProposePartners`): the track
/// (from the iteration's first), the partner's row and sensor, its update's
/// chi2, its signed distance from the track's row along the track [cm], which
/// orders the partners on one side of a row, its row's position among the
/// track's rows and whether it lies after that row.
struct PartnerProposal {
  // As for TrackletT: upstream's `apply` macro also uses the member names as
  // enum constants.
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
  template <template <class> class C>
  struct RecordT {
    MEMLAYOUT_APPLY(RecordT, track, row, surface, chi2, path, layerRow, ahead)
    C<std::uint32_t> track;
    C<std::uint32_t> row;
    C<std::uint32_t> surface;
    C<Scalar> chi2;
    C<Scalar> path;
    C<std::uint8_t> layerRow;
    C<std::uint8_t> ahead;
  };
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
};

namespace overlap
{

/// Each row's filtered state of one fit and its covariance's factors, the
/// inward pass's where it has one, the outward pass's otherwise.
template <std::size_t K>
struct RowStates {
  BoundParameters state[K];
  CovarianceFactors factors[K];
  std::uint8_t pass[K]{}; ///< 0 none, 1 the outward pass's, 2 the inward pass's
};

/// The fit observer that keeps them (`ObservesFitStates`).
template <std::size_t K>
struct RowStateObserver {
  RowStates<K>* states = nullptr;
  GPUhdi() void operator()(PropagationDirection /*direction*/, int /*index*/, Scalar /*chi2*/) const {}
  GPUhdi() void state(PropagationDirection direction, int index, const BoundParameters& state, const CovarianceFactors& factors) const
  {
    const auto at = static_cast<std::size_t>(index);
    const std::uint8_t pass = direction == PropagationDirection::kBackward ? 2U : 1U;
    if (states->pass[at] < pass) {
      states->state[at] = state;
      states->factors[at] = factors;
      states->pass[at] = pass;
    }
  }
};

/// Where a track's state \a state on its row's sensor, on \a layer, may
/// cross another sensor of the layer that the later tests accept: the range
/// in azimuth about the detector's axis and in the second coordinate (z on a
/// barrel layer, the radius on a disc) that the track's line sweeps within
/// the layer's envelope (`LayerSpan`), widened by how far an accepted
/// crossing may lie off that line. That margin is a sensor's bounds widened
/// by \a spread in both local coordinates (up to sqrt(2) spread off the
/// sensor) and the helix's departure from its line, taken as the line test
/// takes it (half \a curvature times the step squared) over the line's
/// longer step to the envelope's edges. The line's sweep is exact: its
/// azimuth changes monotonically along it, so the range is that of its ends.
struct PartnerReach {
  Scalar point[3]{};     ///< the state's point [cm]
  Scalar phi = 0.F;      ///< the range's centre in azimuth
  Scalar second = 0.F;   ///< and in the second coordinate
  Scalar phiReach = 0.F; ///< its half-widths
  Scalar secondReach = 0.F;
};

template <DetectorLayout Det>
[[nodiscard]] GPUhdi() PartnerReach partnerReach(const IterationParameters<Det>& parameters, const TrackingGeometry<Det>& geometry, LayerID layer, const BoundParameters& state, Scalar spread, Scalar curvature)
{
  using Math = o2::gpu::CAMath;
  PartnerReach reach;
  surfaceToGlobal(geometry.planeOf(state), state.loc0(), state.loc1(), reach.point);
  const Scalar r = Math::Hypot(reach.point[0], reach.point[1]);
  const bool barrel = ((parameters.barrelLayers >> layer) & 1U) != 0U;
  const auto& span = parameters.layerSpans[layer];
  // The line in the transverse plane, in the frame where the point is at
  // (r, 0): direction (cosAlpha, sinAlpha), z changing by cotTheta per step.
  const Scalar alpha = state.phi() - azimuth(reach.point[0], reach.point[1]);
  const Scalar cosAlpha = cosOf(alpha);
  const Scalar sinAlpha = sinOf(alpha);
  const Scalar sinTheta = Math::Abs(sinOf(state.theta()));
  const Scalar cotTheta = cosOf(state.theta()) / sinTheta;
  const Scalar dilation = std::numbers::sqrt2_v<Scalar> * spread;
  const auto departure = [&](Scalar step) { return 0.5F * Math::Abs(curvature) * step * step; };
  // The steps back and ahead along the line to the envelope's edges.
  Scalar back = 0.F;
  Scalar ahead = 0.F;
  Scalar margin = 0.F;
  if (barrel) {
    // Radius squared along the line: (s + c)^2 + b^2.
    const Scalar c = r * cosAlpha;
    const Scalar b2 = (r * sinAlpha) * (r * sinAlpha);
    const auto edges = [&](Scalar widening) {
      const Scalar inner = Math::Max(span.positionMin - widening, 0.F);
      const Scalar outer = span.positionMax + widening;
      const auto root = [&](Scalar radius) { return sqrtClampedOf((radius * radius) - b2); };
      // Inwards, the line leaves the envelope at its inner edge if it gets
      // that close to the axis, else at its outer edge beyond.
      const bool dips = b2 < inner * inner;
      ahead = (c < 0.F && dips) ? -c - root(inner) : -c + root(outer);
      back = (c > 0.F && dips) ? c - root(inner) : c + root(outer);
    };
    // The departure over the line's longer step to the edges, as the line
    // test takes it over its step to a sensor's plane.
    edges(dilation);
    margin = dilation + departure(Math::Max(ahead, back) / sinTheta);
    edges(margin);
  } else {
    const Scalar z = reach.point[2];
    const Scalar above = span.positionMax + dilation - z;
    const Scalar below = z - (span.positionMin - dilation);
    margin = dilation + departure(Math::Max(above, below) * Math::Sqrt(1.F + (1.F / (cotTheta * cotTheta))));
    const Scalar up = (above + margin) / Math::Abs(cotTheta);
    const Scalar down = (below + margin) / Math::Abs(cotTheta);
    ahead = cotTheta >= 0.F ? up : down;
    back = cotTheta >= 0.F ? down : up;
  }
  ahead = Math::Max(ahead, 0.F);
  back = Math::Max(back, 0.F);
  const Scalar phiAhead = atan2Of(ahead * sinAlpha, r + (ahead * cosAlpha));
  const Scalar phiBack = atan2Of(-back * sinAlpha, r - (back * cosAlpha));
  const auto radiusAt = [&](Scalar step) { return sqrtClampedOf((r * r) + (2.F * step * r * cosAlpha) + (step * step)); };
  const Scalar closest = -r * cosAlpha;
  const Scalar rMin = (closest > -back && closest < ahead) ? Math::Abs(r * sinAlpha) : Math::Min(radiusAt(ahead), radiusAt(-back));
  const Scalar rMax = Math::Max(radiusAt(ahead), radiusAt(-back));
  // The margin's azimuth about the nearest the line comes to the axis.
  const Scalar widen = rMin > margin ? asinOf(margin / rMin) : kPi<Scalar>;
  reach.phi = azimuth(reach.point[0], reach.point[1]) + (0.5F * (phiAhead + phiBack));
  reach.phiReach = Math::Min((0.5F * Math::Abs(phiAhead - phiBack)) + widen, kPi<Scalar>);
  if (barrel) {
    const Scalar zAhead = reach.point[2] + (ahead * cotTheta);
    const Scalar zBack = reach.point[2] - (back * cotTheta);
    reach.second = 0.5F * (zAhead + zBack);
    reach.secondReach = Math::Min((0.5F * Math::Abs(zAhead - zBack)) + margin, follow::kOpenBound);
  } else {
    reach.second = 0.5F * (rMin + rMax);
    reach.secondReach = Math::Min((0.5F * (rMax - rMin)) + margin, follow::kOpenBound);
  }
  // A NaN range (a slope at a singular angle) looks at the whole layer.
  if (!(reach.phiReach >= 0.F) || !(Math::Abs(reach.phi) <= 4.F * kPi<Scalar>)) {
    reach.phi = 0.F;
    reach.phiReach = kPi<Scalar>;
  }
  if (!(reach.secondReach >= 0.F) || !(Math::Abs(reach.second) <= follow::kOpenBound)) {
    reach.second = 0.F;
    reach.secondReach = follow::kOpenBound;
  }
  return reach;
}

} // namespace overlap

/// One element per track of the iteration (from \a offset): its proposals
/// (see the file's comment), one per other sensor of each layer it has a row
/// on that the layout allows (`overlapLayers`). From each row's state, the
/// layer's sensors whose footprint meets the range the track can cross the
/// layer in (`partnerReach`; listed by `SensorCellsView`, then each tested
/// on its `SensorFootprint`) that the track crosses within their bounds,
/// wider by the gate's spread; on each,
/// the rows within that spread of the crossing, its own rows only. Its
/// pointers may be a device's.
template <DetectorLayout Det>
struct ProposePartners {
  static constexpr std::size_t kChainLayers = maxChainLayers<Det>;
  using TrackView = typename BoundedColumns<TrackOf<kChainLayers>::template RecordT>::ConstView;
  using Record = typename BoundedColumns<PartnerProposal::template RecordT>::Record;

  const IterationParameters<Det>* parameters = nullptr;
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  Propagator propagator;
  SensorCellsView cells;
  TrackView tracks;
  std::uint32_t offset = 0U;
  const std::uint8_t* used = nullptr;
  const std::uint32_t* maxWidthPS = nullptr;   ///< per layer, `TimeframePlan::layerMaxWidthPS`
  const Scalar* measurementVariance = nullptr; ///< per layer, `TimeframePlan::layerMaxLocVariance`
  const SensorFootprint* footprints = nullptr; ///< per sensor, `SensorCells::footprints`

  template <typename Emitter>
  GPUhdi() void operator()(std::uint32_t element, Emitter& emitter) const
  {
    using Math = o2::gpu::CAMath;
    const auto track = offset + element;
    const auto& rows = tracks.clusters[track];
    const auto count = tracks.count[track];
    const auto columns = frame->measurements();
    // A track with no row on a layer whose sensors overlap has no partner.
    bool anyLayer = false;
    for (std::uint8_t k = 0; k < count; ++k) {
      anyLayer = anyLayer || ((overlapLayers<Det> >> geometry->surface(columns.surface[rows.index[k]]).layer) & 1U) != 0U;
    }
    if (!anyLayer) {
      return;
    }
    overlap::RowStates<kChainLayers> states;
    const auto fit = fitMeasurements(propagator, *geometry, detail::LinkSteps{}, RowChain<Det>{.frame = frame, .rows = rows.index}, static_cast<int>(count), kinkBudget(*parameters), overlap::RowStateObserver<kChainLayers>{.states = &states});
    if (fit.failedPass == FitPass::kOutward) {
      return;
    }
    const TimeEstPS box{.start = tracks.timePS[track], .width = tracks.widthPS[track]};
    const Scalar nSigma = Math::Sqrt(parameters->clusterChi2Max);
    for (std::uint8_t k = 0; k < count; ++k) {
      if (states.pass[k] == 0U) {
        continue;
      }
      const auto own = rows.index[k];
      const auto ownSurface = columns.surface[own];
      const LayerID layer = geometry->surface(ownSurface).layer;
      if (((overlapLayers<Det> >> layer) & 1U) == 0U) {
        continue;
      }
      const auto& state = states.state[k];
      const auto& factors = states.factors[k];
      Scalar cov[kBoundCovSize]{};
      covarianceFromFactors(factors, cov);
      // The gate's reach about a predicted crossing: the prediction's
      // variance (`predictedWindow`) and the layer's largest measurement
      // variance, which the gate adds.
      const Scalar spread = nSigma * Math::Sqrt(Math::Max(cov[varIndex(kLoc0)], cov[varIndex(kLoc1)]) + measurementVariance[layer]);
      const Scalar sinTheta = sinOf(state.theta());
      const Scalar direction[3] = {sinTheta * cosOf(state.phi()), sinTheta * sinOf(state.phi()), cosOf(state.theta())};
      const auto& from = geometry->planeOf(state);
      const auto& stepPropagator = pointPropagator(propagator, from, state);
      const ExactScalar curvature = stepPropagator.curvature(static_cast<ExactScalar>(state.qOverP()), 1.);
      const auto reach = overlap::partnerReach<Det>(*parameters, *geometry, layer, state, spread, static_cast<Scalar>(curvature));
      cells.forEachSensorOnce(layer, reach.phi, reach.phiReach, reach.second, reach.secondReach, [&](std::uint32_t id) {
        // The cells list sensors near the range; its own footprint tells,
        // before its surface is read.
        if (id == ownSurface || !footprints[id].meets(reach.phi, reach.phiReach, reach.second, reach.secondReach)) {
          return;
        }
        const TrackingSurface& surface = geometry->surface(static_cast<SensorID>(id));
        if (surface.layer != layer) {
          return;
        }
        const TrackingSurface plane = geometry->planeAt(static_cast<SensorID>(id), reach.point);
        // The straight line along the track first: off the sensor by more
        // than the spread and the helix's departure from the line over the
        // step, it is not crossed.
        Scalar step = 0.F;
        Scalar line[2]{};
        if (!lineCrossing(plane, reach.point, direction, step, line)) {
          return;
        }
        const Scalar departure = 0.5F * static_cast<Scalar>(Math::Abs(curvature)) * step * step;
        if (!surface.bounds.isInside(line[0], line[1], spread + departure)) {
          return;
        }
        // Where the track crosses it, ahead or behind its row.
        BoundParameters at = state;
        auto crossed = propagate(stepPropagator, at, from, plane, PropagationDirection::kForward);
        const auto usable = [&](const PropagationResult& result) { return (result.status == PropagationStatus::kOk || result.status == PropagationStatus::kOutOfBounds) && !follow::beyondHalfTurn(curvature * static_cast<ExactScalar>(result.pathLength)); };
        if (!usable(crossed)) {
          at = state;
          crossed = propagate(stepPropagator, at, from, plane, PropagationDirection::kBackward);
          if (!usable(crossed)) {
            return;
          }
        }
        if (!surface.bounds.isInside(at.loc0(), at.loc1(), spread)) {
          return;
        }
        // Its rows about the crossing.
        Scalar point[3]{};
        surfaceToGlobal(plane, at.loc0(), at.loc1(), point);
        const Scalar x = point[0] - parameters->beamX;
        const Scalar y = point[1] - parameters->beamY;
        const Scalar r = Math::Hypot(x, y);
        const CandidateWindow window{
          .earliestPS = static_cast<std::int64_t>(box.lower()) - static_cast<std::int64_t>(maxWidthPS[layer]) + 1,
          .latestPS = static_cast<std::int64_t>(box.upper()) - 1,
          .afterPS = box.lower(),
          .phiCentre = azimuth(x, y),
          .phiHalfWidth = Math::Min(spread / r, kPi<Scalar>),
          .secondCentre = secondCoordinate<Det>(layer, r, point[2]),
          .secondHalfWidth = spread,
        };
        const auto towards = crossed.pathLength > 0.F ? PropagationDirection::kForward : PropagationDirection::kBackward;
        Record best{};
        bool found = false;
        frame->grid().visitCandidates(layer, window, frame->grid().geometry().allBlocks(), [&](std::uint32_t row) {
          if (columns.surface[row] != id || used[row] != kFree) {
            return;
          }
          if (!spatiallyCompatible(window, columns.phi[row], secondCoordinate<Det>(layer, columns.r[row], columns.z[row])) || !box.isCompatible(timeBoxAt(columns, row))) {
            return;
          }
          BoundParameters predicted;
          CovarianceFactors predictedFactors;
          Scalar chi2 = 0.F;
          // The transport inlined: the follow's shared copy stays its own.
          // A step between two sensors of one layer, which the link material
          // holds no step for: none of its own (`LayerMaterial`).
          if (!follow::updateCandidate<false, true>(propagator, *geometry, detail::LayerMaterial{}, measurementAt(*frame, row), towards, state, factors, predicted, predictedFactors, chi2) || !(chi2 <= parameters->clusterChi2Max)) {
            return;
          }
          if (!found || chi2 < best.chi2 || (!(best.chi2 < chi2) && row < best.row)) {
            best = Record{element, row, id, chi2, crossed.pathLength, k, static_cast<std::uint8_t>(crossed.pathLength > 0.F ? 1U : 0U)};
            found = true;
          }
        });
        if (found) {
          emitter.emit(best);
        }
      });
    }
  }
};

/// One element per track of the iteration (from \a offset) with won
/// partners: fitted again on its rows and partners in path order
/// (`detail::LayerMaterial`), written to \a tracks where the fit succeeds and
/// the rows' time boxes have one in common, `kept` 1 then. Its partners are
/// the proposals `winners` names, from `winnerOffsets[element]` to the next,
/// which this sorts into path order -- by their row's position, then their
/// distance along the track, those before the row first; `merged` holds the
/// rows it fits, from `element * kChainLayers + winnerOffsets[element]` on. Its pointers
/// may be a device's.
template <DetectorLayout Det>
struct RefitWithPartners {
  static constexpr std::size_t kChainLayers = maxChainLayers<Det>;
  using TrackView = typename BoundedColumns<TrackOf<kChainLayers>::template RecordT>::View;
  using ProposalView = typename BoundedColumns<PartnerProposal::template RecordT>::ConstView;

  const IterationParameters<Det>* parameters = nullptr;
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  Propagator propagator;
  TrackView tracks;
  std::uint32_t offset = 0U;
  ProposalView proposals;
  const std::uint32_t* winnerOffsets = nullptr;
  std::uint32_t* winners = nullptr;
  std::uint32_t* merged = nullptr;
  std::uint8_t* layerRowOf = nullptr; ///< beside `merged`, each row's own row's position
  std::uint8_t* kept = nullptr;

  GPUhdi() void operator()(std::uint32_t element) const
  {
    kept[element] = 0U;
    const auto first = winnerOffsets[element];
    const auto end = winnerOffsets[element + 1U];
    if (first == end) {
      return;
    }
    // Into path order: a few per track, sorted in place.
    const auto before = [&](std::uint32_t a, std::uint32_t b) { return proposals.layerRow[a] != proposals.layerRow[b] ? proposals.layerRow[a] < proposals.layerRow[b] : proposals.path[a] < proposals.path[b]; };
    for (auto i = first + 1U; i < end; ++i) {
      const auto moving = winners[i];
      auto at = i;
      for (; at > first && before(moving, winners[at - 1U]); --at) {
        winners[at] = winners[at - 1U];
      }
      winners[at] = moving;
    }
    const auto track = offset + element;
    const auto& rows = tracks.clusters[track];
    const auto count = tracks.count[track];
    const auto base = (static_cast<std::size_t>(element) * kChainLayers) + first;
    std::uint32_t n = 0U;
    auto next = first;
    for (std::uint8_t k = 0; k < count; ++k) {
      for (; next < end && proposals.layerRow[winners[next]] == k && proposals.ahead[winners[next]] == 0U; ++next) {
        merged[base + n] = proposals.row[winners[next]];
        layerRowOf[base + n++] = k;
      }
      merged[base + n] = rows.index[k];
      layerRowOf[base + n++] = k;
      for (; next < end && proposals.layerRow[winners[next]] == k; ++next) {
        merged[base + n] = proposals.row[winners[next]];
        layerRowOf[base + n++] = k;
      }
    }
    const auto columns = frame->measurements();
    auto box = timeBoxAt(columns, merged[base]);
    for (std::uint32_t i = 1; i < n; ++i) {
      const auto other = timeBoxAt(columns, merged[base + i]);
      if (!box.isCompatible(other)) {
        return;
      }
      box += other;
    }
    const auto fit = fitMeasurements(propagator, *geometry, detail::LayerMaterial{}, RowChain<Det>{.frame = frame, .rows = merged + base}, static_cast<int>(n), kinkBudget(*parameters));
    if (fit.status != PropagationStatus::kOk) {
      return;
    }
    tracks.state[track] = fit.state;
    tracks.outwardState[track] = fit.outwardState;
    tracks.timePS[track] = box.start;
    tracks.widthPS[track] = box.width;
    tracks.chi2[track] = fit.chi2;
    tracks.outwardChi2[track] = fit.outwardChi2;
    tracks.ndf[track] = fit.ndf;
    tracks.status[track] = fit.status;
    tracks.kinks[track] = fit.kinks;
    tracks.kinkCluster[track] = fit.kinkMeasurement < 0 ? TrackOf<kChainLayers>::kNoKink : layerRowOf[base + static_cast<std::size_t>(fit.kinkMeasurement)];
    kept[element] = 1U;
  }
};

/// An iteration's overlap partners on the host (see the file's comment): the
/// stage, and its tracks' partners until the store commits them.
template <DetectorLayout Det>
class OverlapPartners
{
 public:
  static constexpr std::size_t kChainLayers = maxChainLayers<Det>;
  using Tracks = typename TrackRefit<Det>::Tracks;
  using Proposals = StageOutput<PartnerProposal::template RecordT>;

  /// What the last run did.
  struct Counts {
    std::size_t proposed = 0U;  ///< partners proposed, one per track and sensor
    std::size_t contested = 0U; ///< of those, lost to another track
    std::size_t refused = 0U;   ///< partners of tracks whose fit with them failed
    std::size_t attached = 0U;
    std::size_t tracks = 0U; ///< tracks with partners
  };

  /// On \a arena, the tracker's.
  explicit OverlapPartners(detail::TaskArena& arena) : mStage(arena) {}

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool = std::move(pool);
    mStage.setMemoryPool(mMemoryPool);
    mProposals.setMemoryPool(mMemoryPool);
    clearStorage();
  }

  /// The run's storage back to the pool.
  void clearStorage() noexcept
  {
    mStage.releaseStorage();
    mProposals.clearStorage();
    deepVectorClear(mMemoryPool.get(), mOffsets, mPartners);
  }

  /// Partners for \a tracks from \a trackOffset on, iteration \a parameters'
  /// (see the file's comment): \a used gets the attached ones, the tracks
  /// with partners their fit with them. \a maxWidthPS per layer the widest
  /// candidate box (`TimeframePlan::layerMaxWidthPS`), \a measurementVariance
  /// the largest local variance (`TimeframePlan::layerMaxLocVariance`).
  /// \a forEach(n, body) runs body(i) for every i below n, serially or in
  /// parallel; \a chunk is the proposals' stage's chunk.
  template <typename ForEach>
  [[nodiscard]] StageStatus run(const ForEach& forEach, std::uint32_t chunk, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, const SensorCellsView& cells, std::span<const SensorFootprint> footprints, Tracks& tracks, std::uint32_t trackOffset, std::span<std::uint8_t> used, std::span<const std::uint32_t> maxWidthPS, std::span<const Scalar> measurementVariance)
  {
    return guardOutOfMemory([&] {
      mCounts = {};
      const auto n = static_cast<std::uint32_t>(tracks.size() - trackOffset);
      // Each track's proposals: count, scan, place.
      mProposals.prepare(n);
      const ProposePartners<Det> propose{.parameters = &parameters, .frame = &frame, .geometry = &geometry, .propagator = propagator, .cells = cells, .tracks = std::as_const(tracks).view(), .offset = trackOffset, .used = used.data(), .maxWidthPS = maxWidthPS.data(), .measurementVariance = measurementVariance.data(), .footprints = footprints.data()};
      const auto staged = mStage.run(mProposals, chunk, propose);
      if (staged != StageStatus::kComplete) {
        return staged;
      }
      const auto proposals = std::as_const(mProposals).view();
      const auto total = static_cast<std::uint32_t>(mProposals.size());
      mCounts.proposed = total;
      // A row proposed by more than one track: the lowest chi2, then the
      // lowest track, takes it.
      std::vector<std::uint32_t> order(total);
      for (std::uint32_t i = 0; i < total; ++i) {
        order[i] = i;
      }
      std::ranges::sort(order, [&](std::uint32_t a, std::uint32_t b) {
        if (proposals.row[a] != proposals.row[b]) {
          return proposals.row[a] < proposals.row[b];
        }
        return proposals.chi2[a] != proposals.chi2[b] ? proposals.chi2[a] < proposals.chi2[b] : proposals.track[a] < proposals.track[b];
      });
      // The winners per track, in their order of proposal.
      std::vector<std::uint32_t> winnerOffsets(std::size_t{n} + 1U, 0U);
      std::vector<std::uint8_t> won(total, 0U);
      for (std::uint32_t at = 0; at < total; ++at) {
        if (at > 0U && proposals.row[order[at]] == proposals.row[order[at - 1U]]) {
          ++mCounts.contested;
          continue;
        }
        won[order[at]] = 1U;
        ++winnerOffsets[proposals.track[order[at]] + 1U];
      }
      for (std::uint32_t i = 0; i < n; ++i) {
        winnerOffsets[i + 1U] += winnerOffsets[i];
      }
      std::vector<std::uint32_t> winners(winnerOffsets[n]);
      {
        std::vector<std::uint32_t> fill(winnerOffsets.begin(), winnerOffsets.end() - 1);
        for (std::uint32_t at = 0; at < total; ++at) {
          if (won[at] != 0U) {
            winners[fill[proposals.track[at]]++] = at;
          }
        }
      }
      // Each track with partners fitted again with them, in parallel.
      std::vector<std::uint32_t> merged((std::size_t{n} * kChainLayers) + winners.size());
      std::vector<std::uint8_t> layerRowOf(merged.size());
      std::vector<std::uint8_t> kept(n, 0U);
      const RefitWithPartners<Det> refit{.parameters = &parameters, .frame = &frame, .geometry = &geometry, .propagator = propagator, .tracks = tracks.view(), .offset = trackOffset, .proposals = proposals, .winnerOffsets = winnerOffsets.data(), .winners = winners.data(), .merged = merged.data(), .layerRowOf = layerRowOf.data(), .kept = kept.data()};
      forEach(n, refit);
      // The partners in path order per track, and the used mask.
      resizeBoundedVector(mOffsets, std::size_t{n} + 1U, mMemoryPool.get());
      mPartners.clear();
      mOffsets[0] = 0U;
      for (std::uint32_t i = 0; i < n; ++i) {
        const auto count = winnerOffsets[i + 1U] - winnerOffsets[i];
        if (kept[i] != 0U) {
          for (auto at = winnerOffsets[i]; at < winnerOffsets[i + 1U]; ++at) {
            const auto proposal = winners[at];
            mPartners.push_back({.row = proposals.row[proposal], .layerRow = proposals.layerRow[proposal], .ahead = proposals.ahead[proposal]});
            used[proposals.row[proposal]] = kExclusive;
          }
          mCounts.attached += count;
          mCounts.tracks += 1U;
        } else {
          mCounts.refused += count;
        }
        mOffsets[i + 1U] = static_cast<std::uint32_t>(mPartners.size());
      }
      return StageStatus::kComplete;
    });
  }

  /// Per track of the last run, from its first: where its partners begin in
  /// `partners`, one more than the tracks.
  [[nodiscard]] std::span<const std::uint32_t> offsets() const { return mOffsets; }
  /// Each track's partners in path order: by their row's position, those
  /// before it first, each side by its distance along the track.
  [[nodiscard]] std::span<const OverlapPartner> partners() const { return mPartners; }
  [[nodiscard]] const Counts& counts() const { return mCounts; }

 private:
  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  CpuStage<PartnerProposal::template RecordT> mStage;
  Proposals mProposals;
  BoundedVector<std::uint32_t> mOffsets;
  BoundedVector<OverlapPartner> mPartners;
  Counts mCounts;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_OVERLAPPARTNERS_H_
