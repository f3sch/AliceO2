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

/// \file JoinStudy.h
/// \brief The track join's dry run (plans/126, S0): after each iteration, its
///        new tracks (pieces) searched inward for an earlier iteration's
///        track ending there, each pair tested and the pairs arbitrated as the
///        join would, nothing applied; every candidate pair against the truth.
///
/// Per iteration, at the point the join would run (its tracks after the
/// follow and the length cut, before the next iteration):
///  - the ends: each earlier track's outermost row;
///  - the search: from each piece's inner state, one backward step onto every
///    layer a link of the layout reaches from the piece's innermost layer, the
///    nominal gate's window there, the grid's rows in it (their time boxes
///    overlapping the piece's) cut to the window, those that are an earlier
///    track's end;
///  - the test: the two tracks' layers disjoint, the combined rows within a
///    chain's capacity, their time boxes overlapping, the same charge; the
///    combined rows (the earlier track's, then the piece's) refitted with a
///    kink budget of the two tracks' own kinks plus one, the earlier track's
///    iteration's kink gate; both passes within the road's chi2 cut at their
///    ndf. Recorded beside: the two fits' q/pT pull, the refit's kinks and
///    where the inward pass's last one updated against the junction, and the
///    deterioration, chi2(joined) - chi2(earlier) - chi2(piece) per degree of
///    freedom the join adds;
///  - the arbitration: each piece its accepted candidate of the least
///    deterioration, each earlier track the piece of the least deterioration
///    among those that chose it (ties by index): the winners.
/// Against the truth: the two tracks' particles (a good track's, else the
/// majority's), the pair the same particle, a parent and its daughter (either
/// way), two daughters of one mother, unrelated, or with a fake; and every
/// pair of one particle's good tracks, an earlier one and a piece with
/// disjoint layers, that the search did not find (kind `missed`).

#ifndef GTRK_STUDY_JOINSTUDY_H_
#define GTRK_STUDY_JOINSTUDY_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/KFTrackFitter.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackRefit.h"
#include "GenericTracking/Study/TrackComparison.h"
#include "GenericTracking/Study/TrackingInput.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <format>
#include <fstream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace o2::gtrk::study
{

template <DetectorLayout Det>
class JoinStudy
{
 public:
  JoinStudy(const std::string& path, const TrackingTruth& truth) : mOut(path), mTruth(&truth)
  {
    if (!mOut) {
      throw std::runtime_error("cannot write " + path);
    }
    mOut << "kind,iteration,piece,earlier,earlierIteration,pieceLayer,targetLayer,earlierLayer,countEarlier,countPiece,disjoint,time,charge,qptPull,"
            "status,chi2,ndf,outwardChi2,kinks,outwardKinks,kinkAt,passes,passesOutward,deterioration,addedNdf,accepted,best,winner,"
            "truth,particleEarlier,particlePiece,pdgEarlier,pdgPiece,originPiece,ptPiece,etaPiece,linked,endDistance\n";
  }

  /// The dry run after one iteration (\a iteration, \a tracker's).
  template <typename Observation, typename Tracker>
  void observe(const Observation& iteration, const Tracker& tracker, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator)
  {
    const auto& tracks = *iteration.tracks;
    const auto view = tracks.view();
    const auto size = static_cast<std::uint32_t>(tracks.size());
    const auto offset = iteration.trackOffset;
    if (mMaxWidthPS.empty()) {
      mMaxWidthPS.assign(Detector<Det>::nLayers, 0U);
      const auto columns = frame.measurements();
      for (std::uint32_t row = 0; row < mTruth->rows.size(); ++row) {
        auto& width = mMaxWidthPS[mTruth->rows[row].layer];
        width = std::max(width, timeBoxAt(columns, row).width);
      }
    }
    const auto rowsOf = [&](std::uint32_t track) { return std::span<const std::uint32_t>{view.clusters[track].index, view.count[track]}; };
    const auto layersOf = [&](std::uint32_t track) {
      std::uint64_t layers = 0U;
      for (const auto row : rowsOf(track)) {
        layers |= std::uint64_t{1} << mTruth->rows[row].layer;
      }
      return layers;
    };
    const auto particleOf = [&](std::uint32_t track) {
      const auto good = goodParticle(mTruth->rows, rowsOf(track));
      return std::pair{good != kNoTrackParticle ? good : majorityParticle(mTruth->rows, rowsOf(track)).particle, good != kNoTrackParticle};
    };
    // The ends: each earlier track's outermost row.
    std::unordered_map<std::uint32_t, std::uint32_t> endOf;
    for (std::uint32_t track = 0; track < offset; ++track) {
      if (view.status[track] == PropagationStatus::kOk && view.count[track] != 0U) {
        endOf.emplace(view.clusters[track].index[view.count[track] - 1U], track);
      }
    }
    std::vector<Pair> pairs;
    std::unordered_set<std::uint64_t> seen;
    for (std::uint32_t piece = offset; piece < size; ++piece) {
      if (view.status[piece] != PropagationStatus::kOk || view.count[piece] == 0U) {
        continue;
      }
      const auto& pieceParameters = tracker.iteration(view.iteration[piece]);
      const BoundTrackState& state = view.state[piece];
      const auto factors = factorCovariance(state.cov);
      const LayerID innermost = mTruth->rows[view.clusters[piece].index[0]].layer;
      const TimeEstPS box{.start = view.timePS[piece], .width = view.widthPS[piece]};
      for (const auto& link : Det::layerLinks) {
        if (link.to != innermost) {
          continue;
        }
        const LayerID target = link.from;
        const auto prediction = predictLayerCrossing(propagator, pieceParameters, geometry.planeOf(state), state, target, PropagationDirection::kBackward);
        if (prediction.status != PropagationStatus::kOk) {
          continue;
        }
        Scalar phiCentre = 0.F;
        Scalar secondCentre = 0.F;
        follow::windowCentre(pieceParameters, target, prediction, phiCentre, secondCentre);
        CandidateWindow window{
          .earliestPS = static_cast<std::int64_t>(box.lower()) - static_cast<std::int64_t>(mMaxWidthPS[target]) + 1,
          .latestPS = static_cast<std::int64_t>(box.upper()) - 1,
          .afterPS = box.lower(),
          .phiCentre = phiCentre,
          .phiHalfWidth = kPi<Scalar>,
          .secondCentre = secondCentre,
          .secondHalfWidth = follow::kOpenBound,
        };
        if (!follow::predictedWindow(propagator, pieceParameters, geometry, detail::LinkSteps{}, state, factors, target, PropagationDirection::kBackward, prediction, window)) {
          continue;
        }
        frame.grid().visitCandidates(target, window, frame.grid().geometry().allBlocks(), [&](std::uint32_t row) {
          const auto measurement = frame.measurement(row);
          if (!spatiallyCompatible(window, measurement.phi, secondCoordinate<Det>(target, measurement.r, measurement.z))) {
            return;
          }
          const auto end = endOf.find(row);
          if (end == endOf.end() || !seen.insert((std::uint64_t{piece} << 32U) | end->second).second) {
            return;
          }
          pairs.push_back(test(tracker, frame, geometry, propagator, view, end->second, piece, target));
        });
      }
    }
    // The arbitration: each piece its accepted pair of the least
    // deterioration, each earlier track the least among those choosing it.
    const auto better = [&](const Pair& a, const Pair& b) { return a.deterioration < b.deterioration || (a.deterioration == b.deterioration && (a.earlier < b.earlier || (a.earlier == b.earlier && a.piece < b.piece))); };
    std::unordered_map<std::uint32_t, std::size_t> bestOfPiece;
    for (std::size_t i = 0; i < pairs.size(); ++i) {
      if (!pairs[i].accepted) {
        continue;
      }
      const auto [at, inserted] = bestOfPiece.emplace(pairs[i].piece, i);
      if (!inserted && better(pairs[i], pairs[at->second])) {
        at->second = i;
      }
    }
    std::unordered_map<std::uint32_t, std::size_t> bestOfEarlier;
    for (const auto& [piece, i] : bestOfPiece) {
      pairs[i].best = true;
      const auto [at, inserted] = bestOfEarlier.emplace(pairs[i].earlier, i);
      if (!inserted && better(pairs[i], pairs[at->second])) {
        at->second = i;
      }
    }
    for (const auto& [earlier, i] : bestOfEarlier) {
      pairs[i].winner = true;
    }
    // The truth's pairs: one particle's good tracks, an earlier one and a
    // piece with disjoint layers, the search's candidates or not.
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> earlierOf;
    for (std::uint32_t track = 0; track < offset; ++track) {
      if (view.status[track] == PropagationStatus::kOk && view.count[track] != 0U) {
        if (const auto [particle, good] = particleOf(track); good) {
          earlierOf[particle].push_back(track);
        }
      }
    }
    for (std::uint32_t piece = offset; piece < size; ++piece) {
      if (view.status[piece] != PropagationStatus::kOk || view.count[piece] == 0U) {
        continue;
      }
      const auto [particle, good] = particleOf(piece);
      const auto found = earlierOf.find(particle);
      if (!good || found == earlierOf.end()) {
        continue;
      }
      for (const auto earlier : found->second) {
        if ((layersOf(earlier) & layersOf(piece)) == 0U && !seen.contains((std::uint64_t{piece} << 32U) | earlier)) {
          Pair missed{.kind = "missed", .piece = piece, .earlier = earlier};
          missed.target = kInvalidLayerID;
          reach(tracker, frame, geometry, propagator, view, missed);
          pairs.push_back(missed);
        }
      }
    }
    for (const auto& pair : pairs) {
      write(iteration.iteration, view, pair, particleOf, rowsOf);
      auto& count = mCounts[pair.kind == std::string_view{"missed"} ? 1U : 0U];
      ++count[0];
      count[1] += pair.accepted ? 1U : 0U;
      count[2] += pair.winner ? 1U : 0U;
    }
  }

  /// The candidates, accepted and winning pairs, and the truth's pairs the
  /// search missed, over the run.
  [[nodiscard]] std::vector<std::string> summary() const
  {
    return {std::format("join study: {} candidate pairs, {} accepted, {} winners; {} pairs of one particle's good tracks the search missed", mCounts[0][0], mCounts[0][1], mCounts[0][2], mCounts[1][0])};
  }

 private:
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  struct Pair {
    const char* kind = "candidate";
    std::uint32_t piece = 0U;
    std::uint32_t earlier = 0U;
    LayerID target = 0;
    bool disjoint = false;
    bool time = false;
    bool charge = false;
    double qptPull = std::numeric_limits<double>::quiet_NaN();
    PropagationStatus status = PropagationStatus::kOk;
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    int ndf = 0;
    double outwardChi2 = std::numeric_limits<double>::quiet_NaN();
    int kinks = 0;
    int outwardKinks = 0;
    int kinkAt = 99;
    bool passes = false;
    bool passesOutward = false;
    double deterioration = std::numeric_limits<double>::infinity();
    int addedNdf = 0;
    bool accepted = false;
    bool best = false;
    bool winner = false;
    /// A missed pair: whether a link of the layout joins the earlier track's
    /// outermost layer to the piece's innermost (-1 not asked), and the
    /// earlier track's end against the piece's backward window there, the
    /// larger of its distances in the window's half widths.
    int linked = -1;
    double endDistance = std::numeric_limits<double>::quiet_NaN();
  };

  /// For a missed \a pair: the link and the earlier track's end against the
  /// piece's backward window on its layer.
  template <typename Tracker, typename View>
  void reach(const Tracker& tracker, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, const View& view, Pair& pair) const
  {
    const auto end = view.clusters[pair.earlier].index[view.count[pair.earlier] - 1U];
    const LayerID endLayer = mTruth->rows[end].layer;
    const LayerID innermost = mTruth->rows[view.clusters[pair.piece].index[0]].layer;
    pair.linked = std::ranges::any_of(Det::layerLinks, [&](const auto& link) { return link.from == endLayer && link.to == innermost; }) ? 1 : 0;
    const auto& parameters = tracker.iteration(view.iteration[pair.piece]);
    const BoundTrackState& state = view.state[pair.piece];
    const auto prediction = predictLayerCrossing(propagator, parameters, geometry.planeOf(state), state, endLayer, PropagationDirection::kBackward);
    if (prediction.status != PropagationStatus::kOk) {
      return;
    }
    Scalar phiCentre = 0.F;
    Scalar secondCentre = 0.F;
    follow::windowCentre(parameters, endLayer, prediction, phiCentre, secondCentre);
    CandidateWindow window{.phiCentre = phiCentre, .phiHalfWidth = kPi<Scalar>, .secondCentre = secondCentre, .secondHalfWidth = follow::kOpenBound};
    if (!follow::predictedWindow(propagator, parameters, geometry, detail::LinkSteps{}, state, factorCovariance(state.cov), endLayer, PropagationDirection::kBackward, prediction, window)) {
      return;
    }
    const auto measurement = frame.measurement(end);
    const double dPhi = std::abs(std::remainder(static_cast<double>(measurement.phi) - static_cast<double>(window.phiCentre), 2. * kPi<double>)) / static_cast<double>(window.phiHalfWidth);
    const double dSecond = std::abs(static_cast<double>(secondCoordinate<Det>(endLayer, measurement.r, measurement.z)) - static_cast<double>(window.secondCentre)) / static_cast<double>(window.secondHalfWidth);
    pair.endDistance = std::max(dPhi, dSecond);
  }

  template <typename Tracker, typename View>
  [[nodiscard]] Pair test(const Tracker& tracker, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, const View& view, std::uint32_t earlier, std::uint32_t piece, LayerID target) const
  {
    Pair pair{.piece = piece, .earlier = earlier, .target = target};
    std::uint64_t earlierLayers = 0U;
    std::uint64_t pieceLayers = 0U;
    for (std::uint8_t i = 0; i < view.count[earlier]; ++i) {
      earlierLayers |= std::uint64_t{1} << mTruth->rows[view.clusters[earlier].index[i]].layer;
    }
    for (std::uint8_t i = 0; i < view.count[piece]; ++i) {
      pieceLayers |= std::uint64_t{1} << mTruth->rows[view.clusters[piece].index[i]].layer;
    }
    pair.disjoint = (earlierLayers & pieceLayers) == 0U;
    const TimeEstPS earlierBox{.start = view.timePS[earlier], .width = view.widthPS[earlier]};
    const TimeEstPS pieceBox{.start = view.timePS[piece], .width = view.widthPS[piece]};
    pair.time = earlierBox.isCompatible(pieceBox);
    const BoundTrackState& outer = view.outwardState[earlier];
    const BoundTrackState& inner = view.state[piece];
    pair.charge = (outer.qOverP() > 0.F) == (inner.qOverP() > 0.F);
    const auto qpt = [](const BoundTrackState& s) { return static_cast<double>(s.qOverP()) / std::sin(static_cast<double>(s.theta())); };
    const auto qptVariance = [](const BoundTrackState& s) {
      const double sinTheta = std::sin(static_cast<double>(s.theta()));
      return static_cast<double>(s.variance(kQOverP)) / (sinTheta * sinTheta);
    };
    pair.qptPull = (qpt(inner) - qpt(outer)) / std::sqrt(qptVariance(inner) + qptVariance(outer));
    const std::uint32_t count = static_cast<std::uint32_t>(view.count[earlier]) + view.count[piece];
    if (!pair.disjoint || count > kMaxClusters) {
      return pair;
    }
    std::array<std::uint32_t, kMaxClusters> rows{};
    std::copy_n(view.clusters[earlier].index, view.count[earlier], rows.begin());
    std::copy_n(view.clusters[piece].index, view.count[piece], rows.begin() + view.count[earlier]);
    const auto& parameters = tracker.iteration(view.iteration[earlier]);
    const KinkBudget budget{.max = static_cast<int>(view.kinks[earlier]) + static_cast<int>(view.kinks[piece]) + 1, .gate = parameters.kinkChi2Min};
    const auto fit = fitMeasurements(propagator, geometry, detail::LinkSteps{}, RowChain<Det>{.frame = &frame, .rows = rows.data()}, static_cast<int>(count), budget);
    pair.status = fit.status;
    if (fit.status != PropagationStatus::kOk) {
      return pair;
    }
    pair.chi2 = static_cast<double>(fit.chi2);
    pair.ndf = fit.ndf;
    pair.outwardChi2 = static_cast<double>(fit.outwardChi2);
    pair.kinks = fit.kinks;
    pair.outwardKinks = fit.outwardKinks;
    pair.kinkAt = fit.kinkMeasurement < 0 ? 99 : fit.kinkMeasurement - static_cast<int>(view.count[earlier]);
    pair.passes = follow::passesRoadChi2(parameters, fit);
    const int outwardNdf = (2 * static_cast<int>(count)) - fittedParameters(propagator) - (2 * fit.outwardKinks);
    const auto cuts = static_cast<int>(std::size(parameters.roadChi2MaxByNdf));
    pair.passesOutward = outwardNdf > 0 && outwardNdf < cuts && fit.outwardChi2 <= parameters.roadChi2MaxByNdf[outwardNdf];
    pair.addedNdf = fit.ndf - view.ndf[earlier] - view.ndf[piece];
    pair.deterioration = pair.addedNdf > 0 ? (static_cast<double>(fit.chi2) - static_cast<double>(view.chi2[earlier]) - static_cast<double>(view.chi2[piece])) / pair.addedNdf : std::numeric_limits<double>::infinity();
    pair.accepted = pair.time && pair.charge && pair.passes && pair.passesOutward;
    return pair;
  }

  template <typename View, typename ParticleOf, typename RowsOf>
  void write(std::size_t iteration, const View& view, const Pair& pair, const ParticleOf& particleOf, const RowsOf& rowsOf)
  {
    const auto [earlierParticle, earlierGood] = particleOf(pair.earlier);
    const auto [pieceParticle, pieceGood] = particleOf(pair.piece);
    const char* truth = "unrelated";
    const auto* earlierTree = mTruth->trees.find(earlierParticle);
    const auto* pieceTree = mTruth->trees.find(pieceParticle);
    if (!earlierGood || !pieceGood) {
      truth = "fake";
    } else if (earlierParticle == pieceParticle) {
      truth = "same";
    } else if (pieceTree != nullptr && pieceTree->mother == earlierParticle) {
      truth = "daughter";
    } else if (earlierTree != nullptr && earlierTree->mother == pieceParticle) {
      truth = "mother";
    } else if (pieceTree != nullptr && earlierTree != nullptr && pieceTree->mother != kNoTreeParticle && pieceTree->mother == earlierTree->mother) {
      truth = "sibling";
    }
    const auto found = mTruth->particles.find(pieceParticle);
    const bool known = found != mTruth->particles.end();
    const auto earlierRows = rowsOf(pair.earlier);
    const auto pieceRows = rowsOf(pair.piece);
    const auto pdg = [&](std::uint64_t particle) {
      const auto at = mTruth->particles.find(particle);
      return at != mTruth->particles.end() ? at->second.pdg : 0;
    };
    mOut << std::format("{},{},{},{},{},{},{},{},{},{},{:d},{:d},{:d},{:.3f},{},{:.4g},{},{:.4g},{},{},{},{:d},{:d},{:.4g},{},{:d},{:d},{:d},{},{},{},{},{},{},{:.4f},{:.4f},{},{:.3f}\n", pair.kind, iteration,
                        pair.piece, pair.earlier, static_cast<int>(view.iteration[pair.earlier]), static_cast<int>(mTruth->rows[pieceRows.front()].layer),
                        pair.target == kInvalidLayerID ? -1 : static_cast<int>(pair.target), static_cast<int>(mTruth->rows[earlierRows.back()].layer), earlierRows.size(), pieceRows.size(),
                        pair.disjoint, pair.time, pair.charge, pair.qptPull, static_cast<int>(pair.status), pair.chi2, pair.ndf, pair.outwardChi2, pair.kinks, pair.outwardKinks, pair.kinkAt,
                        pair.passes, pair.passesOutward, pair.deterioration, pair.addedNdf, pair.accepted, pair.best, pair.winner, truth, earlierParticle, pieceParticle, pdg(earlierParticle),
                        pdg(pieceParticle), known ? static_cast<int>(found->second.origin) : -1, known ? static_cast<double>(found->second.pt) : -1., known ? static_cast<double>(found->second.eta) : 0., pair.linked, pair.endDistance);
  }

  std::ofstream mOut;
  const TrackingTruth* mTruth = nullptr;
  std::vector<std::uint32_t> mMaxWidthPS; ///< per layer, the widest time box
  std::array<std::array<std::size_t, 3>, 2> mCounts{};
};

} // namespace o2::gtrk::study

#endif
