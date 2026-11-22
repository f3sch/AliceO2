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

/// \file DecayCandidateFit.h
/// \brief Decay hypotheses and the candidates a finder makes of them,
///        fitted by index into the daughter table.
///
/// A *hypothesis* is one decay written once per run: its daughters' species
/// and charge signs (its slots), the fit variant, the mother's mass and
/// whether the fit holds it. K0S -> pi+ pi-, Lambda -> p pi-, anti-Lambda ->
/// anti-p pi+ are three hypotheses. A *candidate* is one species assignment
/// of one combination of tracks: a hypothesis, one daughter-table entry per
/// slot (the
/// track under the slot's species), the vertex's start (the finder's: the
/// helices' closest approach, `decayFitStart`, or its own) and, optionally, a
/// reference point (a primary vertex) for the observables. A candidate is
/// one fit: a pair under three hypotheses is three fits, as each species carries
/// its own scattering and loss to the vertex.
///
/// The fit (DecayVertexFit.h) of each candidate gives its result and the
/// observables a finder selects on (DecayObservables.h): the daughters'
/// invariant mass and its sigma, with a point the mother's flight from it in
/// sigma (the point's covariance added) and the cosine of its pointing.
/// Candidates are checked before anything is fitted: a slot whose entry is
/// of another species or charge, a track taken twice, or a start of another
/// count is refused with its status.

#ifndef GTRK_TRACKING_DECAYCANDIDATEFIT_H_
#define GTRK_TRACKING_DECAYCANDIDATEFIT_H_

#include "GenericTracking/Base/FixedVector.h"
#include "GenericTracking/Tracking/DaughterTable.h"
#include "GenericTracking/Tracking/DecayObservables.h"
#include "GenericTracking/Tracking/DecayVertexFit.h"
#include "GenericTracking/Tracking/TrackStore.h"

#include "ReconstructionDataFormats/PID.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace o2::gtrk
{

/// One decay hypothesis (see the file's comment).
struct DecayHypothesis {
  FixedVector<o2::track::PID, decay_fit_detail::kMaxConstrainedDaughters> species; ///< per slot
  FixedVector<std::int8_t, decay_fit_detail::kMaxConstrainedDaughters> charges;    ///< per slot, +1 or -1
  double massGeV = 0.;                                                             ///< the mother's
  bool massConstraint = false;                                                     ///< the fit holds it
  DecayFitVariant variant = DecayFitVariant::kFull;
};

/// A reference point: a position and its covariance (packed, `covIndex`).
struct ReferencePoint {
  Vec3D position{};
  std::array<double, 6> covariance{};
};

/// One candidate (see the file's comment).
struct DecayCandidate {
  std::uint16_t hypothesis = 0U;
  std::int32_t point = -1;                                                        ///< its reference point, -1 for none
  FixedVector<std::uint32_t, decay_fit_detail::kMaxConstrainedDaughters> entries; ///< per slot, into the daughter table
  FixedVector<Vec3D, decay_fit_detail::kMaxConstrainedVertices> start;
};

/// A candidate's fit and observables.
struct DecayCandidateResult {
  DecayVertexFitResult fit;
  DecayKinematics kinematics; ///< all slots'
  DecayFlight flight;         ///< from the reference point; zero without one
  bool hasPoint = false;
};

namespace decay_fit_detail
{

/// Whether \a candidate fits \a hypotheses and \a table: its hypothesis
/// exists, one entry per slot of the slot's species and charge sign, no track
/// twice, at most one start per vertex (kInvalidInput otherwise), every
/// entry's refit succeeded (kRefitFailed otherwise), a start (kNoStart
/// otherwise).
template <DetectorLayout Det>
[[nodiscard]] VertexFitStatus checkCandidate(const DecayCandidate& candidate, std::span<const DecayHypothesis> hypotheses, const DaughterTable<Det>& table, std::size_t points)
{
  if (candidate.hypothesis >= hypotheses.size() || (candidate.point >= 0 && static_cast<std::size_t>(candidate.point) >= points)) {
    return VertexFitStatus::kInvalidInput;
  }
  const auto& hypothesis = hypotheses[candidate.hypothesis];
  if (candidate.entries.size() != hypothesis.species.size() || hypothesis.charges.size() != hypothesis.species.size() || candidate.start.size() > 1U) {
    return VertexFitStatus::kInvalidInput;
  }
  const auto view = table.view();
  for (std::size_t i = 0; i < candidate.entries.size(); ++i) {
    const auto e = candidate.entries[i];
    if (e >= table.size() || table.species()[view.species[e]] != hypothesis.species[i]) {
      return VertexFitStatus::kInvalidInput;
    }
    if (view.status[e] != PropagationStatus::kOk) {
      return VertexFitStatus::kRefitFailed;
    }
    if ((view.state[e].qOverP() > 0.F ? 1 : -1) != hypothesis.charges[i]) {
      return VertexFitStatus::kInvalidInput;
    }
    for (std::size_t j = 0; j < i; ++j) {
      if (view.track[candidate.entries[j]] == view.track[e]) {
        return VertexFitStatus::kInvalidInput;
      }
    }
  }
  return candidate.start.empty() ? VertexFitStatus::kNoStart : VertexFitStatus::kOk;
}

} // namespace decay_fit_detail

/// The daughters of \a candidate (checked: `checkCandidate`) as the fit
/// takes them, from \a table's entries of \a tracks in \a frame.
template <DetectorLayout Det>
[[nodiscard]] FixedVector<DecayDaughter, decay_fit_detail::kMaxConstrainedDaughters> candidateDaughters(const DecayCandidate& candidate, const DaughterTable<Det>& table, const TrackStore<Det>& tracks, const TrackingGeometry<Det>& geometry, const TimeFrame<Det>& frame, const Propagator& base)
{
  FixedVector<DecayDaughter, decay_fit_detail::kMaxConstrainedDaughters> daughters;
  const auto view = table.view();
  for (const auto e : candidate.entries) {
    const auto pid = table.species()[view.species[e]];
    auto propagator = base;
    propagator.massGeV = pid.getMass();
    DecayDaughter daughter{.innermost = view.state[e], .surface = geometry.planeOf(measurementAt(frame, tracks.rows(view.track[e]).front())), .propagator = propagator, .mass = propagator.massGeV, .chargeMagnitude = std::abs(static_cast<double>(pid.getCharge()))};
    if (view.mixture[e] >= 0) {
      daughter.mixture = table.mixture(static_cast<std::size_t>(view.mixture[e]));
    }
    daughters.push_back(daughter);
  }
  return daughters;
}

/// \a candidates fitted into \a results (one each, in their order) on
/// \a forEach(n, body), with \a hypotheses, \a table's entries of \a tracks in
/// \a frame, \a points their reference points, \a material between each
/// daughter's innermost measurement and the vertex.
template <DetectorLayout Det, typename ForEach, typename Material>
void fitDecayCandidates(const ForEach& forEach, std::span<const DecayCandidate> candidates, std::span<const DecayHypothesis> hypotheses, const DaughterTable<Det>& table, const TrackStore<Det>& tracks, const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, const Propagator& propagator, const Material& material, std::span<const ReferencePoint> points, std::span<DecayCandidateResult> results)
{
  forEach(static_cast<std::uint32_t>(candidates.size()), [&](std::uint32_t c) {
    const auto& candidate = candidates[c];
    auto& out = results[c];
    out = {};
    if (const auto status = decay_fit_detail::checkCandidate(candidate, hypotheses, table, points.size()); status != VertexFitStatus::kOk) {
      out.fit.status = status;
      return;
    }
    const auto& hypothesis = hypotheses[candidate.hypothesis];
    const auto daughters = candidateDaughters(candidate, table, tracks, geometry, frame, propagator);
    const std::span<const DecayDaughter> span{daughters.data(), daughters.size()};
    DecayConstraints constraints;
    const DecayTopology topology; // one vertex; the fitter holds it by reference
    if (hypothesis.massConstraint) {
      constraints.mass = true;
      constraints.massGeV = hypothesis.massGeV;
    }
    withDecayFitVariant(hypothesis.variant, [&](auto variant) {
      DecayVertexFitter<Det, decltype(variant)::value, Material> fitter(span, material, constraints, topology, {});
      std::vector<std::byte> workspace(fitter.workspace());
      fitter.fit({candidate.start.data(), candidate.start.size()}, workspace, out.fit);
    });
    if (out.fit.status != VertexFitStatus::kOk) {
      return;
    }
    std::array<double, decay_fit_detail::kMaxConstrainedDaughters> masses{};
    for (std::size_t i = 0; i < daughters.size(); ++i) {
      masses[i] = daughters[i].mass;
    }
    const std::span<const Vec3D> momenta{out.fit.momenta.data(), out.fit.momenta.size()};
    const std::span<const double> charges{out.fit.charges.data(), out.fit.charges.size()};
    const int vertices = static_cast<int>(out.fit.vertices.size());
    out.kinematics = decayKinematics(momenta, std::span<const double>{masses.data(), daughters.size()}, charges, [&](int a, int b) { return out.fit.covarianceAt((3 * vertices) + a, (3 * vertices) + b); });
    if (candidate.point >= 0) {
      const auto& point = points[static_cast<std::size_t>(candidate.point)];
      out.flight = decayFlight(out.fit.vertex, momenta, charges, point.position, point.covariance, hypothesis.massGeV, vertices, [&](int a, int b) { return out.fit.covarianceAt(a, b); });
      out.hasPoint = true;
    }
  });
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_DECAYCANDIDATEFIT_H_
