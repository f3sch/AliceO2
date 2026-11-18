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
///  - each kept track is fitted on its rows, with the tracking's material
///    (its rows are one per layer, which the surface map holds), and each
///    row's filtered state kept, the inward pass's where it has one
///    (`ObservesFitStates`);
///  - from each row's state, the other sensors of its layer the track
///    crosses -- listed near enough by `SensorCells` (`partnerReach`), and
///    crossed within their bounds widened by the gate's spread -- and on each
///    its free rows, compatible in time, within that spread of the crossing,
///    each predicted on the sensor and tested with the follow's
///    two-coordinate cut (`clusterChi2Max`); per sensor the best is the
///    track's proposal;
///  - a row proposed by more than one track goes to the lowest chi2, then the
///    lowest track (a deterministic order);
///  - each track with partners is fitted again on its rows and partners in
///    path order, with the surface map between layers and no material of its
///    own on the step between a layer's two sensors (`LayerMaterial`: the
///    map's halves hold the layer's whole material); where that fit fails or the rows' time
///    boxes leave none in common, the track stays as it was, without them;
///  - the attached partners are taken in the used mask, exclusively.
///
/// The tracks' own rows stay one per layer (`TrackOf`); the partners go to
/// the store as references (`kOverlapPartner`, TrackStore.h), in path order
/// beside the rows. A track's states, chi2 and time are then those of the
/// fit with its partners.

#ifndef GTRK_TRACKING_OVERLAPPARTNERS_H_
#define GTRK_TRACKING_OVERLAPPARTNERS_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/LayerLinks.h"
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
  void operator()(PropagationDirection /*direction*/, int /*index*/, Scalar /*chi2*/) const {}
  void state(PropagationDirection direction, int index, const BoundParameters& state, const CovarianceFactors& factors) const
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
/// cross another sensor of the layer: about its point there (the detector's
/// frame, as `SensorCells` lists the sensors), as far as the layer's whole
/// envelope lets the crossing move along the track's slope, both sensors
/// anywhere in it, plus \a spread.
struct PartnerReach {
  Scalar point[3]{}; ///< the state's point [cm]
  Scalar phi = 0.F;  ///< its azimuth about the detector's axis
  Scalar second = 0.F;
  Scalar phiReach = 0.F;
  Scalar secondReach = 0.F;
};

template <DetectorLayout Det>
[[nodiscard]] PartnerReach partnerReach(const IterationParameters<Det>& parameters, const TrackingGeometry<Det>& geometry, LayerID layer, const BoundParameters& state, Scalar spread)
{
  using Math = o2::gpu::CAMath;
  PartnerReach reach;
  surfaceToGlobal(geometry.planeOf(state), state.loc0(), state.loc1(), reach.point);
  const Scalar r = Math::Hypot(reach.point[0], reach.point[1]);
  const bool barrel = ((parameters.barrelLayers >> layer) & 1U) != 0U;
  reach.phi = azimuth(reach.point[0], reach.point[1]);
  reach.second = barrel ? reach.point[2] : r;
  const Scalar alpha = state.phi() - reach.phi;
  const Scalar sinAlpha = Math::Abs(Math::Sin(alpha));
  const Scalar cosAlpha = Math::Abs(Math::Cos(alpha));
  const Scalar sinTheta = Math::Abs(Math::Sin(state.theta()));
  const Scalar cosTheta = Math::Abs(Math::Cos(state.theta()));
  // How far the crossing moves per unit of the layer's position (radius on
  // a barrel layer, z on a disc), along the track there.
  const Scalar transverse = barrel ? 1.F / cosAlpha : sinTheta / cosTheta;
  const auto& span = parameters.layerSpans[layer];
  const Scalar envelope = span.positionMax - span.positionMin;
  reach.phiReach = Math::Min(((transverse * sinAlpha * envelope) + spread) / r, std::numbers::pi_v<Scalar>);
  reach.secondReach = Math::Min((transverse * (barrel ? cosTheta / sinTheta : cosAlpha) * envelope) + spread, follow::kOpenBound);
  // A NaN reach (a slope at a singular angle) looks at the whole layer.
  if (!(reach.phiReach >= 0.F)) {
    reach.phiReach = std::numbers::pi_v<Scalar>;
  }
  if (!(reach.secondReach >= 0.F)) {
    reach.secondReach = follow::kOpenBound;
  }
  return reach;
}

} // namespace overlap

/// An iteration's overlap partners (see the file's comment): the stage, and
/// its tracks' partners until the store commits them.
template <DetectorLayout Det>
class OverlapPartners
{
 public:
  static constexpr std::size_t K = maxChainLayers<Det>;
  using Tracks = typename TrackRefit<Det>::Tracks;

  /// A track's proposal: a partner and its update's chi2; `path` its
  /// signed distance from the track's row along the track [cm], which
  /// orders the partners on one side of a row.
  struct Proposal {
    std::uint32_t track = 0U; ///< from the iteration's first track
    OverlapPartner partner;
    std::uint32_t surface = 0U;
    Scalar chi2 = 0.F;
    Scalar path = 0.F;
  };

  /// What the last run did.
  struct Counts {
    std::size_t proposed = 0U;  ///< partners proposed, one per track and sensor
    std::size_t contested = 0U; ///< of those, lost to another track
    std::size_t refused = 0U;   ///< partners of tracks whose fit with them failed
    std::size_t attached = 0U;
    std::size_t tracks = 0U; ///< tracks with partners
  };

  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool = std::move(pool);
    clearStorage();
  }

  /// The run's storage back to the pool.
  void clearStorage() noexcept
  {
    deepVectorClear(mMemoryPool.get(), mOffsets, mPartners);
    mProposals = {};
  }

  /// Partners for \a tracks from \a trackOffset on, iteration \a parameters'
  /// (see the file's comment): \a used gets the attached ones, the tracks
  /// with partners their fit with them. \a maxWidthPS per layer the widest
  /// candidate box (`TimeframePlan::layerMaxWidthPS`), \a measurementVariance
  /// the largest local variance (`TimeframePlan::layerMaxLocVariance`). \a forEach(n, body)
  /// runs body(i) for every i below n, serially or in parallel.
  template <typename Material, typename ForEach>
  [[nodiscard]] StageStatus run(const ForEach& forEach, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& lut, const Propagator& propagator, const SensorCellsView& cells, Tracks& tracks, std::uint32_t trackOffset, std::span<std::uint8_t> used, std::span<const std::uint32_t> maxWidthPS, std::span<const Scalar> measurementVariance)
  {
    return guardOutOfMemory([&] {
      mCounts = {};
      const auto n = static_cast<std::uint32_t>(tracks.size() - trackOffset);
      const auto view = tracks.view();
      const detail::LayerMaterial<Material> layered{.lut = &lut};
      // Each track's proposals, in parallel, each track its own list.
      std::vector<std::vector<Proposal>> proposals(n);
      forEach(n, [&](std::uint32_t i) { propose(parameters, frame, geometry, lut, propagator, cells, view, trackOffset + i, i, used, maxWidthPS, measurementVariance, proposals[i]); });
      // A row proposed by more than one track: the lowest chi2, then the
      // lowest track, takes it.
      for (const auto& own : proposals) {
        mProposals.insert(mProposals.end(), own.begin(), own.end());
      }
      mCounts.proposed = mProposals.size();
      std::ranges::sort(mProposals, [](const Proposal& a, const Proposal& b) {
        if (a.partner.row != b.partner.row) {
          return a.partner.row < b.partner.row;
        }
        return a.chi2 != b.chi2 ? a.chi2 < b.chi2 : a.track < b.track;
      });
      std::vector<std::vector<Proposal>> won(n);
      for (std::size_t at = 0; at < mProposals.size(); ++at) {
        if (at > 0U && mProposals[at].partner.row == mProposals[at - 1U].partner.row) {
          ++mCounts.contested;
          continue;
        }
        won[mProposals[at].track].push_back(mProposals[at]);
      }
      // Each track with partners fitted again with them, in parallel.
      std::vector<std::uint8_t> kept(n, 0U);
      forEach(n, [&](std::uint32_t i) {
        if (!won[i].empty()) {
          kept[i] = refit(parameters, frame, geometry, layered, propagator, view, trackOffset + i, won[i]) ? 1U : 0U;
        }
      });
      // The partners in path order per track, and the used mask.
      resizeBoundedVector(mOffsets, std::size_t{n} + 1U, mMemoryPool.get());
      mPartners.clear();
      mOffsets[0] = 0U;
      for (std::uint32_t i = 0; i < n; ++i) {
        if (kept[i] != 0U) {
          for (const auto& proposal : won[i]) {
            mPartners.push_back(proposal.partner);
            used[proposal.partner.row] = kExclusive;
          }
          mCounts.attached += won[i].size();
          ++mCounts.tracks;
        } else {
          mCounts.refused += won[i].size();
        }
        mOffsets[i + 1U] = static_cast<std::uint32_t>(mPartners.size());
      }
      mProposals = {};
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
  /// Track \a track's proposals into \a out, one per other sensor of each
  /// layer it has a row on (see the file's comment). From each row's state,
  /// the layer's sensors near enough (`partnerReach`, `SensorCells`) that
  /// the track crosses within their bounds, wider by the gate's spread; on
  /// each, the rows within that spread of the crossing, its own rows only.
  template <typename Material, typename View>
  static void propose(const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator, const SensorCellsView& cells, const View& view, std::uint32_t track, std::uint32_t ordinal, std::span<const std::uint8_t> used, std::span<const std::uint32_t> maxWidthPS, std::span<const Scalar> measurementVariance, std::vector<Proposal>& out)
  {
    using Math = o2::gpu::CAMath;
    const auto& rows = view.clusters[track];
    const auto count = view.count[track];
    overlap::RowStates<K> states;
    const auto fit = fitMeasurements(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = rows.index}, static_cast<int>(count), kinkBudget(parameters), overlap::RowStateObserver<K>{.states = &states});
    if (fit.failedPass == FitPass::kOutward) {
      return;
    }
    const auto columns = frame.measurements();
    const TimeEstPS box{.start = view.timePS[track], .width = view.widthPS[track]};
    const Scalar nSigma = Math::Sqrt(parameters.clusterChi2Max);
    std::vector<std::uint32_t> seen;
    for (std::uint8_t k = 0; k < count; ++k) {
      if (states.pass[k] == 0U) {
        continue;
      }
      const auto own = rows.index[k];
      const auto ownSurface = columns.surface[own];
      const LayerID layer = geometry.surface(ownSurface).layer;
      const auto& state = states.state[k];
      const auto& factors = states.factors[k];
      Scalar cov[kBoundCovSize]{};
      covarianceFromFactors(factors, cov);
      // The gate's reach about a predicted crossing: the prediction's
      // variance and the follow's position term (`predictedWindow`), and the
      // layer's largest measurement variance, which the gate adds.
      const Scalar spread = nSigma * Math::Sqrt(Math::Max(cov[varIndex(kLoc0)], cov[varIndex(kLoc1)]) + parameters.followPositionVariance + measurementVariance[layer]);
      const Scalar sinTheta = Math::Sin(state.theta());
      const Scalar direction[3] = {sinTheta * Math::Cos(state.phi()), sinTheta * Math::Sin(state.phi()), Math::Cos(state.theta())};
      const auto reach = overlap::partnerReach<Det>(parameters, geometry, layer, state, spread);
      const auto& from = geometry.planeOf(state);
      const auto& stepPropagator = pointPropagator(propagator, from, state);
      const ExactScalar curvature = stepPropagator.curvature(static_cast<ExactScalar>(state.qOverP()), 1.);
      seen.clear();
      static_cast<void>(cells.anySensor(layer, reach.phi, reach.phiReach, reach.second, reach.secondReach, [&](std::uint32_t id) {
        if (id == ownSurface || std::ranges::find(seen, id) != seen.end()) {
          return false;
        }
        seen.push_back(id);
        const TrackingSurface& surface = geometry.surface(static_cast<SensorID>(id));
        if (surface.layer != layer) {
          return false;
        }
        TrackingSurface plane = surface;
        if (surface.radius > 0.F) {
          Scalar along0 = 0.F;
          Scalar along1 = 0.F;
          for (int i = 0; i < 3; ++i) {
            along0 += (reach.point[i] - surface.origin[i]) * surface.axis0[i];
            along1 += (reach.point[i] - surface.origin[i]) * surface.axis1[i];
          }
          plane = tangentPlane(surface, Math::ATan2(along1, along0));
        }
        // The straight line along the track first: off the sensor by more than
        // the spread and the helix's departure from the line over the step,
        // it is not crossed.
        Scalar normalDirection = 0.F;
        Scalar normalOffset = 0.F;
        for (int i = 0; i < 3; ++i) {
          normalDirection += plane.normal[i] * direction[i];
          normalOffset += plane.normal[i] * (plane.origin[i] - reach.point[i]);
        }
        if (!(Math::Abs(normalDirection) > 0.F)) {
          return false;
        }
        const Scalar step = normalOffset / normalDirection;
        Scalar line0 = 0.F;
        Scalar line1 = 0.F;
        for (int i = 0; i < 3; ++i) {
          const Scalar offset = reach.point[i] + (step * direction[i]) - plane.origin[i];
          line0 += offset * plane.axis0[i];
          line1 += offset * plane.axis1[i];
        }
        const Scalar departure = 0.5F * static_cast<Scalar>(Math::Abs(curvature)) * step * step;
        if (Math::Abs(line0) > surface.bounds.halfX + spread + departure || Math::Abs(line1) > surface.bounds.halfY + spread + departure) {
          return false;
        }
        // Where the track crosses it, ahead or behind its row.
        BoundParameters at = state;
        auto crossed = propagate(stepPropagator, at, from, plane, PropagationDirection::kForward);
        const auto usable = [&](const PropagationResult& result) { return (result.status == PropagationStatus::kOk || result.status == PropagationStatus::kOutOfBounds) && !follow::beyondHalfTurn(curvature * static_cast<ExactScalar>(result.pathLength)); };
        if (!usable(crossed)) {
          at = state;
          crossed = propagate(stepPropagator, at, from, plane, PropagationDirection::kBackward);
          if (!usable(crossed)) {
            return false;
          }
        }
        if (Math::Abs(at.loc0()) > surface.bounds.halfX + spread || Math::Abs(at.loc1()) > surface.bounds.halfY + spread) {
          return false;
        }
        // Its rows about the crossing.
        Scalar point[3]{};
        surfaceToGlobal(plane, at.loc0(), at.loc1(), point);
        const Scalar x = point[0] - parameters.beamX;
        const Scalar y = point[1] - parameters.beamY;
        const Scalar r = Math::Hypot(x, y);
        const CandidateWindow window{
          .earliestPS = static_cast<std::int64_t>(box.lower()) - static_cast<std::int64_t>(maxWidthPS[layer]) + 1,
          .latestPS = static_cast<std::int64_t>(box.upper()) - 1,
          .afterPS = box.lower(),
          .phiCentre = azimuth(x, y),
          .phiHalfWidth = Math::Min(spread / r, std::numbers::pi_v<Scalar>),
          .secondCentre = secondCoordinate<Det>(layer, r, point[2]),
          .secondHalfWidth = spread,
        };
        const auto towards = crossed.pathLength > 0.F ? PropagationDirection::kForward : PropagationDirection::kBackward;
        Proposal best;
        bool found = false;
        frame.grid().visitCandidates(layer, window, frame.grid().geometry().allBlocks(), [&](std::uint32_t row) {
          if (columns.surface[row] != id || used[row] != kFree) {
            return;
          }
          if (!spatiallyCompatible(window, columns.phi[row], secondCoordinate<Det>(layer, columns.r[row], columns.z[row])) || !box.isCompatible(timeBoxAt(columns, row))) {
            return;
          }
          BoundParameters predicted;
          CovarianceFactors predictedFactors;
          Scalar chi2 = 0.F;
          if (!follow::updateCandidate<false>(propagator, geometry, material, measurementAt(frame, row), towards, state, factors, predicted, predictedFactors, chi2, 1.F) || !(chi2 <= parameters.clusterChi2Max)) {
            return;
          }
          if (!found || chi2 < best.chi2 || (!(best.chi2 < chi2) && row < best.partner.row)) {
            best = {.track = ordinal, .partner = {.row = row, .layerRow = k, .ahead = static_cast<std::uint8_t>(crossed.pathLength > 0.F ? 1U : 0U)}, .surface = id, .chi2 = chi2, .path = crossed.pathLength};
            found = true;
          }
        });
        if (found) {
          out.push_back(best);
        }
        return false;
      }));
    }
  }

  /// Track \a track fitted on its rows and \a partners, which this sorts
  /// into path order (by their row's position, then their distance along the
  /// track, those before the row first); written to \a view where the fit succeeds and the rows'
  /// time boxes have one in common. Whether it was.
  template <typename Material, typename View>
  static bool refit(const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator, const View& view, std::uint32_t track, std::vector<Proposal>& partners)
  {
    std::ranges::sort(partners, [](const Proposal& a, const Proposal& b) {
      if (a.partner.layerRow != b.partner.layerRow) {
        return a.partner.layerRow < b.partner.layerRow;
      }
      return a.path < b.path;
    });
    const auto& rows = view.clusters[track];
    const auto count = view.count[track];
    std::vector<std::uint32_t> merged;
    std::vector<std::uint8_t> layerRowOf;
    merged.reserve(count + partners.size());
    std::size_t next = 0U;
    for (std::uint8_t k = 0; k < count; ++k) {
      for (; next < partners.size() && partners[next].partner.layerRow == k && partners[next].partner.ahead == 0U; ++next) {
        merged.push_back(partners[next].partner.row);
        layerRowOf.push_back(k);
      }
      merged.push_back(rows.index[k]);
      layerRowOf.push_back(k);
      for (; next < partners.size() && partners[next].partner.layerRow == k; ++next) {
        merged.push_back(partners[next].partner.row);
        layerRowOf.push_back(k);
      }
    }
    const auto columns = frame.measurements();
    auto box = timeBoxAt(columns, merged.front());
    for (std::size_t i = 1; i < merged.size(); ++i) {
      const auto other = timeBoxAt(columns, merged[i]);
      if (!box.isCompatible(other)) {
        return false;
      }
      box += other;
    }
    const auto fit = fitMeasurements(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = merged.data()}, static_cast<int>(merged.size()), kinkBudget(parameters));
    if (fit.status != PropagationStatus::kOk) {
      return false;
    }
    view.state[track] = fit.state;
    view.outwardState[track] = fit.outwardState;
    view.timePS[track] = box.start;
    view.widthPS[track] = box.width;
    view.chi2[track] = fit.chi2;
    view.outwardChi2[track] = fit.outwardChi2;
    view.ndf[track] = fit.ndf;
    view.status[track] = fit.status;
    view.kinks[track] = fit.kinks;
    view.kinkCluster[track] = fit.kinkMeasurement < 0 ? TrackOf<K>::kNoKink : layerRowOf[static_cast<std::size_t>(fit.kinkMeasurement)];
    return true;
  }

  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  BoundedVector<std::uint32_t> mOffsets;
  BoundedVector<OverlapPartner> mPartners;
  std::vector<Proposal> mProposals;
  Counts mCounts;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_OVERLAPPARTNERS_H_
