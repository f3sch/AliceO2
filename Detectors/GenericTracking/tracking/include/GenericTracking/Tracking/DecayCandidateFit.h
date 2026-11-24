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
/// and charge signs (its track slots), its composites (other hypotheses, of
/// tracks only: a cascade's Lambda), the fit variant, the mother's mass and
/// whether the fit holds it and the composites' masses. K0S -> pi+ pi-,
/// Lambda -> p pi-, anti-Lambda -> anti-p pi+ are three hypotheses. A
/// *candidate* is one species assignment of one combination of tracks: a
/// hypothesis, one daughter-table entry per slot (the track under the slot's
/// species), the vertices' start (the finder's: the helices' closest
/// approaches, `startDecayCandidate`, or its own) and, optionally, a
/// reference point (a primary vertex) for the observables. A candidate is
/// one fit: a pair under three hypotheses is three fits, as each species
/// carries its own scattering and loss to the vertex.
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
#include "GenericTracking/Tracking/Species.h"
#include "GenericTracking/Tracking/TrackStore.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>

namespace o2::gtrk
{

/// One decay hypothesis (see the file's comment).
struct DecayHypothesis {
  FixedVector<Species, decay_fit_detail::kMaxConstrainedDaughters> species;     ///< per own track slot
  FixedVector<std::int8_t, decay_fit_detail::kMaxConstrainedDaughters> charges; ///< per own track slot, +1 or -1
  /// Its composites' hypotheses, each of tracks only (a cascade's Lambda).
  FixedVector<std::uint16_t, decay_fit_detail::kMaxConstrainedVertices - 1> composites;
  double massGeV = 0.;         ///< the mother's
  bool massConstraint = false; ///< the fit holds it
  bool holdComposites = false; ///< the fit holds the composites' masses
  DecayFitVariant variant = DecayFitVariant::kFull;
};

/// A hypothesis as the fit's tree (`expandHypothesis`): its leaves -- each
/// composite's tracks in turn, then its own --, their species and charges,
/// the topology (vertex 0 the hypothesis's, vertex 1 + c composite c's) and
/// each composite's leaves and mass.
struct ExpandedHypothesis {
  VertexFitStatus status = VertexFitStatus::kOk; ///< kInvalidInput or kTooLarge for a hypothesis the fit cannot take
  FixedVector<Species, decay_fit_detail::kMaxConstrainedDaughters> species;
  FixedVector<std::int8_t, decay_fit_detail::kMaxConstrainedDaughters> charges;
  DecayTopology topology;
  std::array<FixedVector<int, decay_fit_detail::kMaxConstrainedDaughters>, decay_fit_detail::kMaxConstrainedVertices - 1> compositeLeaves{};
  FixedVector<double, decay_fit_detail::kMaxConstrainedVertices - 1> compositeMasses;
};

/// \a hypotheses' \a index expanded (`ExpandedHypothesis`).
[[nodiscard]] inline ExpandedHypothesis expandHypothesis(std::span<const DecayHypothesis> hypotheses, std::size_t index)
{
  ExpandedHypothesis out;
  if (index >= hypotheses.size()) {
    out.status = VertexFitStatus::kInvalidInput;
    return out;
  }
  const auto& hypothesis = hypotheses[index];
  std::size_t leaves = hypothesis.species.size();
  for (const auto c : hypothesis.composites) {
    if (c >= hypotheses.size() || !hypotheses[c].composites.empty()) {
      out.status = VertexFitStatus::kInvalidInput;
      return out;
    }
    leaves += hypotheses[c].species.size();
  }
  if (leaves > static_cast<std::size_t>(decay_fit_detail::kMaxConstrainedDaughters) || hypothesis.charges.size() != hypothesis.species.size()) {
    out.status = leaves > static_cast<std::size_t>(decay_fit_detail::kMaxConstrainedDaughters) ? VertexFitStatus::kTooLarge : VertexFitStatus::kInvalidInput;
    return out;
  }
  out.topology.massGeV.push_back(hypothesis.massGeV);
  if (!hypothesis.composites.empty()) {
    out.topology.parent.push_back(-1);
  }
  for (std::size_t c = 0; c < hypothesis.composites.size(); ++c) {
    const auto& composite = hypotheses[hypothesis.composites[c]];
    if (composite.charges.size() != composite.species.size()) {
      out.status = VertexFitStatus::kInvalidInput;
      return out;
    }
    out.topology.parent.push_back(0);
    out.topology.massGeV.push_back(composite.massGeV);
    out.compositeMasses.push_back(composite.massGeV);
    for (std::size_t s = 0; s < composite.species.size(); ++s) {
      out.compositeLeaves[c].push_back(static_cast<int>(out.species.size()));
      out.species.push_back(composite.species[s]);
      out.charges.push_back(composite.charges[s]);
      out.topology.vertexOf.push_back(static_cast<int>(c) + 1);
    }
  }
  for (std::size_t s = 0; s < hypothesis.species.size(); ++s) {
    out.species.push_back(hypothesis.species[s]);
    out.charges.push_back(hypothesis.charges[s]);
    if (!hypothesis.composites.empty()) {
      out.topology.vertexOf.push_back(0);
    }
  }
  return out;
}

/// One candidate (see the file's comment).
struct DecayCandidate {
  std::uint16_t hypothesis = 0U;
  std::int32_t point = -1; ///< its reference point, -1 for none
  /// Per leaf of its expanded hypothesis (`expandHypothesis`), into the
  /// daughter table.
  FixedVector<std::uint32_t, decay_fit_detail::kMaxConstrainedDaughters> entries;
  FixedVector<Vec3D, decay_fit_detail::kMaxConstrainedVertices> start; ///< one per vertex of its tree
};

/// A fit's daughters' kinematics (`decayKinematics`) with their \a masses.
[[nodiscard]] inline DecayKinematics fitKinematics(const DecayVertexFitResult& fit, std::span<const double> masses)
{
  const int vertices = static_cast<int>(fit.vertices.size());
  return decayKinematics(std::span<const Vec3D>{fit.momenta.data(), fit.momenta.size()}, masses, std::span<const double>{fit.charges.data(), fit.charges.size()}, [&](int a, int b) { return fit.covarianceAt((3 * vertices) + a, (3 * vertices) + b); });
}

/// The flight of a fit's root, of mass \a mass, from \a point (`decayFlight`).
[[nodiscard]] inline DecayFlight fitFlight(const DecayVertexFitResult& fit, const MeasuredPoint& point, double mass)
{
  return decayFlight(fit.vertex, std::span<const Vec3D>{fit.momenta.data(), fit.momenta.size()}, std::span<const double>{fit.charges.data(), fit.charges.size()}, point.position, point.covariance, mass, static_cast<int>(fit.vertices.size()), [&](int a, int b) { return fit.covarianceAt(a, b); });
}

/// A candidate's fit and observables.
struct DecayCandidateResult {
  DecayVertexFitResult fit;
  DecayKinematics kinematics; ///< all leaves'
  DecayFlight flight;         ///< from the reference point; zero without one
  bool hasPoint = false;
  /// Per composite: its leaves' kinematics and its flight from its parent's vertex.
  FixedVector<DecayKinematics, decay_fit_detail::kMaxConstrainedVertices - 1> composites;
  FixedVector<DecayFlight, decay_fit_detail::kMaxConstrainedVertices - 1> compositeFlights;
};

namespace decay_fit_detail
{

/// Whether \a candidate fits its expanded hypothesis \a tree and \a table:
/// one entry per leaf of the leaf's species and charge sign, no track twice,
/// no more starts than vertices (kInvalidInput otherwise), every entry's
/// refit succeeded (kRefitFailed otherwise), a start per vertex (kNoStart
/// otherwise).
template <DetectorLayout Det>
[[nodiscard]] VertexFitStatus checkCandidate(const DecayCandidate& candidate, const ExpandedHypothesis& tree, const DaughterTable<Det>& table, std::size_t points)
{
  if (tree.status != VertexFitStatus::kOk) {
    return tree.status;
  }
  const auto vertices = static_cast<std::size_t>(tree.topology.vertices());
  if ((candidate.point >= 0 && std::cmp_greater_equal(candidate.point, points)) || candidate.entries.size() != tree.species.size() || candidate.start.size() > vertices) {
    return VertexFitStatus::kInvalidInput;
  }
  const auto view = table.view();
  for (std::size_t i = 0; i < candidate.entries.size(); ++i) {
    const auto e = candidate.entries[i];
    if (e >= table.size() || table.species()[view.species[e]] != tree.species[i]) {
      return VertexFitStatus::kInvalidInput;
    }
    if (view.status[e] != PropagationStatus::kOk) {
      return VertexFitStatus::kRefitFailed;
    }
    if ((view.state[e].qOverP() > 0.F ? 1 : -1) != tree.charges[i]) {
      return VertexFitStatus::kInvalidInput;
    }
    for (std::size_t j = 0; j < i; ++j) {
      if (view.track[candidate.entries[j]] == view.track[e]) {
        return VertexFitStatus::kInvalidInput;
      }
    }
  }
  return candidate.start.size() < vertices ? VertexFitStatus::kNoStart : VertexFitStatus::kOk;
}

} // namespace decay_fit_detail

/// Where a candidate's daughters come from: the daughter table's entries of
/// the timeframe's committed tracks, and the run's field.
template <DetectorLayout Det>
struct DaughterSource {
  const DaughterTable<Det>* table = nullptr;
  const TrackStore<Det>* tracks = nullptr;
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  Propagator propagator;

  /// Committed track \a track's innermost measurement's plane: where its
  /// entries' states are.
  [[nodiscard]] TrackingSurface innermostPlane(std::uint32_t track) const { return geometry->planeOf(measurementAt(*frame, tracks->rows(track).front())); }
};

/// Daughter-table entry \a e of \a source as the fit takes it: its state,
/// its plane, the field with its species' mass and |q|, an electron's mixture.
template <DetectorLayout Det>
[[nodiscard]] DecayDaughter daughterOf(const DaughterSource<Det>& source, std::uint32_t e)
{
  const auto view = source.table->view();
  auto propagator = source.propagator;
  propagator.setSpecies(source.table->species()[view.species[e]]);
  DecayDaughter daughter{.innermost = view.state[e], .surface = source.innermostPlane(view.track[e]), .propagator = propagator, .mass = propagator.massGeV};
  if (view.mixture[e] >= 0) {
    daughter.mixture = source.table->mixture(static_cast<std::size_t>(view.mixture[e]));
  }
  return daughter;
}

/// The daughters of \a candidate (its leaves, checked: `checkCandidate`) as
/// the fit takes them (`daughterOf`).
template <DetectorLayout Det>
[[nodiscard]] FixedVector<DecayDaughter, decay_fit_detail::kMaxConstrainedDaughters> candidateDaughters(const DecayCandidate& candidate, const DaughterSource<Det>& source)
{
  FixedVector<DecayDaughter, decay_fit_detail::kMaxConstrainedDaughters> daughters;
  for (const auto e : candidate.entries) {
    daughters.push_back(daughterOf(source, e));
  }
  return daughters;
}

/// The constraints \a hypothesis holds, \a tree its expansion: the
/// mother's mass and the composites', as it asks.
[[nodiscard]] inline DecayConstraints candidateConstraints(const DecayHypothesis& hypothesis, const ExpandedHypothesis& tree)
{
  DecayConstraints constraints;
  if (hypothesis.massConstraint) {
    constraints.mass = true;
    constraints.massGeV = hypothesis.massGeV;
  }
  if (hypothesis.holdComposites) {
    for (const double m : tree.compositeMasses) {
      constraints.compositeMasses.push_back(m);
    }
  }
  return constraints;
}

/// \a candidate's start where it has none: its daughters'
/// (`candidateDaughters`) start under its hypothesis (`decayFitStart`), one
/// per vertex, if it passes `checkCandidate` with \a points reference
/// points. Its check's status: kOk with a start, kNoStart without one.
template <DetectorLayout Det>
[[nodiscard]] VertexFitStatus startDecayCandidate(DecayCandidate& candidate, std::span<const DecayHypothesis> hypotheses, const DaughterSource<Det>& source, std::size_t points)
{
  const auto tree = expandHypothesis(hypotheses, candidate.hypothesis);
  const auto status = decay_fit_detail::checkCandidate(candidate, tree, *source.table, points);
  if (status != VertexFitStatus::kNoStart || !candidate.start.empty()) {
    return status;
  }
  const auto daughters = candidateDaughters(candidate, source);
  std::array<Vec3D, decay_fit_detail::kMaxConstrainedVertices> start{};
  if (!decayFitStart<Det>(std::span<const DecayDaughter>{daughters.data(), daughters.size()}, candidateConstraints(hypotheses[candidate.hypothesis], tree), tree.topology, {}, start)) {
    return VertexFitStatus::kNoStart;
  }
  for (int w = 0; w < tree.topology.vertices(); ++w) {
    candidate.start.push_back(start[static_cast<std::size_t>(w)]);
  }
  return VertexFitStatus::kOk;
}

/// \a candidate fitted into \a out with its observables, with
/// \a hypotheses, \a source's entries, \a points its reference points,
/// \a material between each daughter's innermost measurement and its vertex;
/// the fit's workspace from \a workspace; \a mothers its charged particles'
/// own rows, if they left any (`DecayMother`: a tracked Xi); with
/// \a production, the mother made there (`DecayConstraints::production`: a
/// measurement of its helix, with its scattering on the way and the point's
/// own covariance).
template <DetectorLayout Det, typename Material>
void fitDecayCandidate(const DecayCandidate& candidate, std::span<const DecayHypothesis> hypotheses, const DaughterSource<Det>& source, const Material& material, std::span<const MeasuredPoint> points, WorkspaceBuffer& workspace, DecayCandidateResult& out, std::span<const DecayMother> mothers = {}, const MeasuredPoint* production = nullptr)
{
  out = {};
  const auto tree = expandHypothesis(hypotheses, candidate.hypothesis);
  if (const auto status = decay_fit_detail::checkCandidate(candidate, tree, *source.table, points.size()); status != VertexFitStatus::kOk) {
    out.fit.status = status;
    return;
  }
  const auto& hypothesis = hypotheses[candidate.hypothesis];
  const auto daughters = candidateDaughters(candidate, source);
  const std::span<const DecayDaughter> span{daughters.data(), daughters.size()};
  auto constraints = candidateConstraints(hypothesis, tree);
  if (production != nullptr) {
    constraints.production = true;
    constraints.point = production->position;
    constraints.pointCovariance = production->covariance;
  }
  withDecayFitVariant(hypothesis.variant, [&](auto variant) {
    DecayVertexFitter<Det, decltype(variant)::value, Material> fitter(span, material, constraints, tree.topology, mothers);
    fitter.fit({candidate.start.data(), candidate.start.size()}, workspace.take(fitter.workspace()), out.fit);
  });
  if (out.fit.status != VertexFitStatus::kOk) {
    return;
  }
  std::array<double, decay_fit_detail::kMaxConstrainedDaughters> masses{};
  for (std::size_t i = 0; i < daughters.size(); ++i) {
    masses[i] = daughters[i].mass;
  }
  const std::span<const double> massSpan{masses.data(), daughters.size()};
  const std::span<const Vec3D> momenta{out.fit.momenta.data(), out.fit.momenta.size()};
  const std::span<const double> charges{out.fit.charges.data(), out.fit.charges.size()};
  const int vertices = static_cast<int>(out.fit.vertices.size());
  const auto at = [&](int a, int b) { return out.fit.covarianceAt(a, b); };
  out.kinematics = fitKinematics(out.fit, massSpan);
  if (candidate.point >= 0) {
    out.flight = fitFlight(out.fit, points[static_cast<std::size_t>(candidate.point)], hypothesis.massGeV);
    out.hasPoint = true;
  }
  const std::span<const Vec3D> fitted{out.fit.vertices.data(), out.fit.vertices.size()};
  for (std::size_t k = 0; k < tree.compositeMasses.size(); ++k) {
    const std::span<const int> leaves{tree.compositeLeaves[k].data(), tree.compositeLeaves[k].size()};
    out.composites.push_back(subsetKinematics(momenta, massSpan, charges, leaves, vertices, at));
    out.compositeFlights.push_back(compositeFlight(fitted, static_cast<int>(k) + 1, tree.topology.parent[k + 1], momenta, charges, leaves, tree.compositeMasses[k], at));
  }
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_DECAYCANDIDATEFIT_H_
