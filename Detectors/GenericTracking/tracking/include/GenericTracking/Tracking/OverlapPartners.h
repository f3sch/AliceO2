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
///  - each kept track is fitted on its rows walking the LUT
///    (`ExactMaterial`), and each row's filtered state kept, the inward
///    pass's where it has one (`ObservesFitStates`);
///  - about each row, the rows of the same layer on other sensors within the
///    window the state's spread and the layer's envelope allow
///    (`partnerWindow`), free, compatible in time, are each predicted on
///    their own sensor and tested with the follow's two-coordinate cut
///    (`clusterChi2Max`); per sensor the best is the track's proposal;
///  - a row proposed by more than one track goes to the lowest chi2, then the
///    lowest track (a deterministic order);
///  - each track with partners is fitted again on its rows and partners in
///    path order, walking the LUT; where that fit fails or the rows' time
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

/// The window about a track's row on \a layer, its \a state there with
/// \a factors, in which its crossing of any other sensor of the layer lies:
/// the layer's whole envelope times the track's slope there (the two
/// sensors are anywhere in it), plus the gate's spread of the state, its
/// covariance and the follow's position term. In the grid's coordinates
/// about the beam line; the time range the track's \a box less the layer's
/// widest candidate box \a maxWidthPS, as the follow's.
template <DetectorLayout Det>
[[nodiscard]] CandidateWindow partnerWindow(const IterationParameters<Det>& parameters, const TrackingGeometry<Det>& geometry, LayerID layer, const BoundParameters& state, const CovarianceFactors& factors, const TimeEstPS& box, std::uint32_t maxWidthPS)
{
  using Math = o2::gpu::CAMath;
  Scalar point[3]{};
  surfaceToGlobal(geometry.planeOf(state), state.loc0(), state.loc1(), point);
  const Scalar x = point[0] - parameters.beamX;
  const Scalar y = point[1] - parameters.beamY;
  const Scalar r = Math::Hypot(x, y);
  const bool barrel = ((parameters.barrelLayers >> layer) & 1U) != 0U;
  Scalar cov[kBoundCovSize]{};
  covarianceFromFactors(factors, cov);
  const Scalar nSigma = Math::Sqrt(parameters.clusterChi2Max);
  const Scalar spread = nSigma * Math::Sqrt(Math::Max(cov[varIndex(kLoc0)], cov[varIndex(kLoc1)]) + parameters.followPositionVariance);
  const Scalar alpha = state.phi() - azimuth(x, y);
  const Scalar sinAlpha = Math::Abs(Math::Sin(alpha));
  const Scalar cosAlpha = Math::Abs(Math::Cos(alpha));
  const Scalar sinTheta = Math::Abs(Math::Sin(state.theta()));
  const Scalar cosTheta = Math::Abs(Math::Cos(state.theta()));
  // How far the crossing moves per unit of the layer's position (radius on
  // a barrel layer, z on a disc), along the track there.
  const Scalar transverse = barrel ? 1.F / cosAlpha : sinTheta / cosTheta;
  const Scalar phiSlope = transverse * sinAlpha / r;
  const Scalar secondSlope = transverse * (barrel ? cosTheta / sinTheta : cosAlpha);
  const auto& span = parameters.layerSpans[layer];
  const Scalar envelope = span.positionMax - span.positionMin;
  CandidateWindow window{
    .earliestPS = static_cast<std::int64_t>(box.lower()) - static_cast<std::int64_t>(maxWidthPS) + 1,
    .latestPS = static_cast<std::int64_t>(box.upper()) - 1,
    .afterPS = box.lower(),
    .phiCentre = azimuth(x, y),
    .phiHalfWidth = Math::Min((spread / r) + (envelope * phiSlope), std::numbers::pi_v<Scalar>),
    .secondCentre = secondCoordinate<Det>(layer, r, point[2]),
    .secondHalfWidth = Math::Min(spread + (envelope * secondSlope), follow::kOpenBound),
  };
  // A NaN half-width (a slope at a singular angle) searches the whole layer.
  if (!(window.phiHalfWidth >= 0.F)) {
    window.phiHalfWidth = std::numbers::pi_v<Scalar>;
  }
  if (!(window.secondHalfWidth >= 0.F)) {
    window.secondHalfWidth = follow::kOpenBound;
  }
  return window;
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
  /// candidate box (`TimeframePlan::layerMaxWidthPS`). \a forEach(n, body)
  /// runs body(i) for every i below n, serially or in parallel.
  template <typename Material, typename ForEach>
  [[nodiscard]] StageStatus run(const ForEach& forEach, const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& lut, const Propagator& propagator, Tracks& tracks, std::uint32_t trackOffset, std::span<std::uint8_t> used, std::span<const std::uint32_t> maxWidthPS)
  {
    return guardOutOfMemory([&] {
      mCounts = {};
      const auto n = static_cast<std::uint32_t>(tracks.size() - trackOffset);
      const auto view = tracks.view();
      const detail::ExactMaterial<Material> exact{.lut = &lut};
      // Each track's proposals, in parallel, each track its own list.
      std::vector<std::vector<Proposal>> proposals(n);
      forEach(n, [&](std::uint32_t i) { propose(parameters, frame, geometry, exact, propagator, view, trackOffset + i, i, used, maxWidthPS, proposals[i]); });
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
          kept[i] = refit(parameters, frame, geometry, exact, propagator, view, trackOffset + i, won[i]) ? 1U : 0U;
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
  /// layer it has a row on (see the file's comment).
  template <typename Material, typename View>
  static void propose(const IterationParameters<Det>& parameters, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Material& material, const Propagator& propagator, const View& view, std::uint32_t track, std::uint32_t ordinal, std::span<const std::uint8_t> used, std::span<const std::uint32_t> maxWidthPS, std::vector<Proposal>& out)
  {
    const auto& rows = view.clusters[track];
    const auto count = view.count[track];
    overlap::RowStates<K> states;
    const auto fit = fitMeasurements(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = rows.index}, static_cast<int>(count), kinkBudget(parameters), overlap::RowStateObserver<K>{.states = &states});
    if (fit.failedPass == FitPass::kOutward) {
      return;
    }
    const auto columns = frame.measurements();
    const TimeEstPS box{.start = view.timePS[track], .width = view.widthPS[track]};
    for (std::uint8_t k = 0; k < count; ++k) {
      if (states.pass[k] == 0U) {
        continue;
      }
      const auto own = rows.index[k];
      const auto ownSurface = columns.surface[own];
      const LayerID layer = geometry.surface(ownSurface).layer;
      const auto& state = states.state[k];
      const auto window = overlap::partnerWindow<Det>(parameters, geometry, layer, state, states.factors[k], box, maxWidthPS[layer]);
      // The track's direction at its row: the side of the row a candidate
      // lies on, and the way to propagate there.
      const Scalar sinTheta = o2::gpu::CAMath::Sin(state.theta());
      const Scalar direction[3] = {sinTheta * o2::gpu::CAMath::Cos(state.phi()), sinTheta * o2::gpu::CAMath::Sin(state.phi()), o2::gpu::CAMath::Cos(state.theta())};
      const std::size_t first = out.size();
      frame.grid().visitCandidates(layer, window, frame.grid().geometry().allBlocks(), [&](std::uint32_t row) {
        const auto surface = columns.surface[row];
        if (row == own || surface == ownSurface || used[row] != kFree || geometry.surface(surface).layer != layer) {
          return;
        }
        if (!spatiallyCompatible(window, columns.phi[row], secondCoordinate<Det>(layer, columns.r[row], columns.z[row])) || !box.isCompatible(timeBoxAt(columns, row))) {
          return;
        }
        const Scalar path = (direction[0] * (columns.x[row] - columns.x[own])) + (direction[1] * (columns.y[row] - columns.y[own])) + (direction[2] * (columns.z[row] - columns.z[own]));
        const auto towards = path > 0.F ? PropagationDirection::kForward : PropagationDirection::kBackward;
        BoundParameters predicted;
        CovarianceFactors predictedFactors;
        Scalar chi2 = 0.F;
        if (!follow::updateCandidate<false>(propagator, geometry, material, measurementAt(frame, row), towards, state, states.factors[k], predicted, predictedFactors, chi2, 1.F) || !(chi2 <= parameters.clusterChi2Max)) {
          return;
        }
        // The best per sensor.
        for (std::size_t at = first; at < out.size(); ++at) {
          if (out[at].surface == surface) {
            if (chi2 < out[at].chi2 || (!(out[at].chi2 < chi2) && row < out[at].partner.row)) {
              out[at] = {.track = ordinal, .partner = {.row = row, .layerRow = k, .ahead = static_cast<std::uint8_t>(path > 0.F ? 1U : 0U)}, .surface = surface, .chi2 = chi2, .path = path};
            }
            return;
          }
        }
        out.push_back({.track = ordinal, .partner = {.row = row, .layerRow = k, .ahead = static_cast<std::uint8_t>(path > 0.F ? 1U : 0U)}, .surface = surface, .chi2 = chi2, .path = path});
      });
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
