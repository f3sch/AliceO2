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

/// \file PrimaryVertexStudy.h
/// \brief Each collision's primary vertex fitted from the chains of the
///        particles made at it (PrimaryVertexFit.h) and held to the MC.
///
/// The association is the truth's: a collision's tracks are the chains of
/// its physical primaries made at its vertex, one per particle (its
/// longest that runs with the particle: a track of its charge, a truth
/// chain's first leg -- a looper's way back, reconstructed outgoing, has the
/// opposite charge and comes from no vertex), each refitted with its
/// particle's true species (`refitRowsWithHypothesis`); nothing else -- no secondary, no fake, no
/// other collision's track, no physical primary made away from the vertex
/// (a charm or beauty decay's daughter) -- enters. Made at the vertex: a
/// production point within 1 um of it, far above the float production
/// points' rounding (~1e-3 um), far below a heavy-flavour flight. The
/// prior is the run's beam line (its MeanVertex's position and widths). The
/// start: the beam line's transverse point, z the tracks' z at their perigee
/// about it, weighted by its variance. Everything else is the fit's own: its
/// covariance, chi2 and tails.
///
/// A record per collision with a chain: the fit against the MC vertex.

#ifndef GTRK_STUDY_PRIMARYVERTEXSTUDY_H_
#define GTRK_STUDY_PRIMARYVERTEXSTUDY_H_

#include "GenericTracking/Base/Species.h"
#include "GenericTracking/Tracking/HypothesisRefit.h"
#include "GenericTracking/Tracking/PrimaryVertexFit.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Study/MCTruth.h"
#include "GenericTracking/Study/RecordWriter.h"
#include "GenericTracking/Study/StudyChains.h"
#include "GenericTracking/Study/TruthEvent.h"
#include "GenericTracking/Study/TruthFit.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace o2::gtrk::study
{

/// One collision's primary vertex fit (`PrimaryVertexStudy`).
struct PrimaryVertexRecord {
  static constexpr std::string_view kName = "primaryVertex";
  std::uint32_t timeframe = 0U;
  std::uint8_t source = 0U;       ///< ChainSource
  std::uint64_t collision = 0U;   ///< `collisionOf` its particles
  std::uint16_t tracks = 0U;      ///< fitted
  std::uint16_t refused = 0U;     ///< chains of another species or charge
  std::uint16_t refitFailed = 0U; ///< chains whose rows' refit or carry to the beam line failed
  std::uint16_t madeAway = 0U;    ///< physical primaries' chains left out, made away from the vertex
  std::uint8_t status = 0U;       ///< VertexFitStatus
  std::uint8_t iterations = 0U;
  std::uint8_t outer = 0U;
  std::uint8_t stopReason = 0U;
  float chi2 = 0.F;
  float priorChi2 = 0.F;
  std::int16_t ndf = 0;
  float microseconds = 0.F; ///< the fit's, the refits apart
  float mcX = 0.F;          ///< the MC vertex [cm]
  float mcY = 0.F;
  float mcZ = 0.F;
  float startDZ = 0.F; ///< the start's z less the MC's [um]
  float dX = 0.F;      ///< fitted less MC [um]
  float dY = 0.F;
  float dZ = 0.F;
  float sX = 0.F; ///< the fit's sigmas [um]
  float sY = 0.F;
  float sZ = 0.F;
  float rhoXY = 0.F; ///< the fit's correlations
  float rhoXZ = 0.F;
  float rhoYZ = 0.F;
  float maxTrackChi2 = 0.F;  ///< the largest track's chi2 at the solution
  std::int32_t worstPdg = 0; ///< and its particle's species
  float worstPt = 0.F;       ///< [GeV]
  /// A failed fit: the tracks whose removal alone lets it converge, and of
  /// those the one leaving the smallest chi2 (`MockPrimaryVertex::culprit`).
  std::uint16_t rescues = 0U;
  std::int32_t culpritPdg = 0;
  float culpritPt = 0.F; ///< [GeV]
  float culpritEta = 0.F;
  std::uint16_t culpritRows = 0U;
  std::uint8_t culpritChains = 0U; ///< its particle's chains (tracks) in the frame
  float culpritQOverPPull = 0.F;   ///< its refit's innermost q/p less the MC production's, in sigma

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("source", source);
    visit("collision", collision);
    visit("tracks", tracks);
    visit("refused", refused);
    visit("refitFailed", refitFailed);
    visit("madeAway", madeAway);
    visit("status", status);
    visit("iterations", iterations);
    visit("outer", outer);
    visit("stopReason", stopReason);
    visit("chi2", chi2);
    visit("priorChi2", priorChi2);
    visit("ndf", ndf);
    visit("microseconds", microseconds);
    visit("mcX", mcX);
    visit("mcY", mcY);
    visit("mcZ", mcZ);
    visit("startDZ", startDZ);
    visit("dX", dX);
    visit("dY", dY);
    visit("dZ", dZ);
    visit("sX", sX);
    visit("sY", sY);
    visit("sZ", sZ);
    visit("rhoXY", rhoXY);
    visit("rhoXZ", rhoXZ);
    visit("rhoYZ", rhoYZ);
    visit("maxTrackChi2", maxTrackChi2);
    visit("worstPdg", worstPdg);
    visit("worstPt", worstPt);
    visit("rescues", rescues);
    visit("culpritPdg", culpritPdg);
    visit("culpritPt", culpritPt);
    visit("culpritEta", culpritEta);
    visit("culpritRows", culpritRows);
    visit("culpritChains", culpritChains);
    visit("culpritQOverPPull", culpritQOverPPull);
  }
};

/// A collision's mocked primary vertex (`mockPrimaryVertices`).
struct MockPrimaryVertex {
  std::uint64_t collision = 0U;
  std::uint16_t madeAway = 0U;    ///< its physical primaries' chains made away from the vertex, left out
  std::uint16_t refused = 0U;     ///< chains of another species or charge
  std::uint16_t refitFailed = 0U; ///< chains whose rows' refit or carry to the beam line failed
  std::uint16_t tracks = 0U;      ///< fitted
  Vec3D start{};
  PrimaryVertexFitResult fit; ///< kNoStart without a fitted track
  float microseconds = 0.F;   ///< the fit's, the refits apart
  float maxTrackChi2 = 0.F;   ///< the largest track's chi2 at the solution
  std::uint64_t worstParticle = 0U;
  /// A failed fit refitted without each track in turn: how many removals let
  /// it converge, and the particle and rows of the one leaving the smallest
  /// chi2 (without a rescue: none).
  std::uint16_t rescues = 0U;
  std::uint64_t culprit = 0U;
  std::uint16_t culpritRows = 0U;
  std::uint8_t culpritChains = 0U;
  float culpritQOverPPull = 0.F;
  /// Its fitted chains' committed tracks (the chains from the tracker): its
  /// contributors.
  std::vector<std::uint32_t> contributors;
};

namespace primary_vertex_detail
{
/// Made at the vertex: within this of it [cm] (the file's comment).
inline constexpr double kAtVertexCM = 1e-4;
} // namespace primary_vertex_detail

/// Each collision's mocked primary vertex from \a chains of \a event (the
/// file's comment), the collisions in their first chain's order.
template <DetectorLayout Det>
[[nodiscard]] std::vector<MockPrimaryVertex> mockPrimaryVertices(const TruthEvent<Det>& event, std::span<const Chain> chains)
{
  std::unordered_map<std::uint64_t, const Chain*> longest;
  std::unordered_map<std::uint64_t, std::uint32_t> chainsOf;
  for (const auto& chain : chains) {
    ++chainsOf[chain.particle];
    const double charge = particleCharge(event.particles->at(chain.particle).pdg);
    const bool outgoing = chain.track ? (chain.track->state.qOverP() > 0.F) == (charge > 0.) : chain.leg == 0U;
    if (!outgoing) {
      continue;
    }
    auto& best = longest[chain.particle];
    if (best == nullptr || chain.rows.size() > best->rows.size()) {
      best = &chain;
    }
  }
  std::unordered_map<std::uint64_t, std::size_t> index;
  std::vector<MockPrimaryVertex> vertices;
  std::vector<std::vector<const Chain*>> members;
  for (const auto& chain : chains) {
    if (const auto found = longest.find(chain.particle); found == longest.end() || found->second != &chain) {
      continue;
    }
    const auto& truth = event.particles->at(chain.particle);
    if (truth.origin != ParticleOrigin::kPhysicalPrimary) {
      continue;
    }
    const auto collision = collisionOf(chain.particle);
    const auto [at, inserted] = index.try_emplace(collision, vertices.size());
    if (inserted) {
      vertices.push_back({.collision = collision, .fit = {}, .contributors = {}});
      members.emplace_back();
    }
    const auto offset = productionOffset(*event.collisions, chain.particle, Vec3D{truth.vx, truth.vy, truth.vz});
    if (!offset || *offset > primary_vertex_detail::kAtVertexCM) {
      ++vertices[at->second].madeAway;
      continue;
    }
    members[at->second].push_back(&chain);
  }
  const auto& frame = *event.frame;
  const auto& geometry = *event.geometry;
  const auto& beam = event.beamLine;
  MeasuredPoint prior{.position = {beam.x, beam.y, beam.z}};
  prior.covariance[covIndex(0, 0)] = beam.sigmaX * beam.sigmaX;
  prior.covariance[covIndex(1, 1)] = beam.sigmaY * beam.sigmaY;
  prior.covariance[covIndex(2, 2)] = beam.sigmaZ * beam.sigmaZ;
  o2::gtrk::detail::parallelFor(*event.arena, 0U, vertices.size(), 8, [&](std::size_t begin, std::size_t end) {
    for (std::size_t v = begin; v < end; ++v) {
      auto& out = vertices[v];
      constexpr o2::gtrk::detail::LayerMaterial material{};
      std::vector<VertexTrack> tracks;
      std::vector<const Chain*> fitted;
      double weightSum = 0.;
      double weightedZ = 0.;
      for (const auto* chain : members[v]) {
        const auto& truth = event.particles->at(chain->particle);
        const auto species = speciesOf(truth.pdg);
        if (!species || std::lround(std::abs(particleCharge(truth.pdg))) != 1) {
          ++out.refused;
          continue;
        }
        const RowChain<Det> rows{.frame = &frame, .rows = chain->rows.data()};
        auto refit = refitRowsWithHypothesis(event.propagator, geometry, material, rows, static_cast<std::int32_t>(chain->rows.size()), *species);
        const auto& surface = geometry.planeOf(rows[0]);
        if (refit.status != PropagationStatus::kOk || carryToPerigee(event.propagator, geometry, surface, refit, static_cast<Scalar>(beam.x), static_cast<Scalar>(beam.y)) != PropagationStatus::kOk) {
          ++out.refitFailed;
          continue;
        }
        auto propagator = event.propagator;
        propagator.setSpecies(*species);
        tracks.push_back({.innermost = refit.innermost, .surface = surface, .propagator = propagator});
        fitted.push_back(chain);
        const auto variance = static_cast<double>(refit.perigee.covAt(kLoc1, kLoc1));
        if (variance > 0.) {
          weightSum += 1. / variance;
          weightedZ += static_cast<double>(refit.perigee.loc1()) / variance;
        }
      }
      out.tracks = static_cast<std::uint16_t>(std::min<std::size_t>(tracks.size(), 65535U));
      if (tracks.empty() || !(weightSum > 0.)) {
        continue; // kNoStart
      }
      out.start = {beam.x, beam.y, weightedZ / weightSum};
      const auto started = std::chrono::steady_clock::now();
      std::vector<double> trackChi2(tracks.size());
      out.fit = fitPrimaryVertex(std::span<const VertexTrack>{tracks}, geometry.material(), prior, out.start, {}, std::span<double>{trackChi2});
      out.microseconds = static_cast<float>(std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - started).count());
      if (out.fit.status != VertexFitStatus::kOk) {
        double best = std::numeric_limits<double>::infinity();
        std::vector<VertexTrack> others;
        for (std::size_t i = 0; i < tracks.size() && tracks.size() > 1U; ++i) {
          others.clear();
          for (std::size_t j = 0; j < tracks.size(); ++j) {
            if (j != i) {
              others.push_back(tracks[j]);
            }
          }
          const auto without = fitPrimaryVertex(std::span<const VertexTrack>{others}, geometry.material(), prior, out.start);
          if (without.status == VertexFitStatus::kOk) {
            ++out.rescues;
            if (without.chi2 < best) {
              best = without.chi2;
              out.culprit = fitted[i]->particle;
              out.culpritRows = static_cast<std::uint16_t>(fitted[i]->rows.size());
              out.culpritChains = static_cast<std::uint8_t>(std::min(chainsOf.contains(out.culprit) ? chainsOf.at(out.culprit) : 0U, 255U));
              const auto& truth = event.particles->at(out.culprit);
              const double trueQOverP = particleCharge(truth.pdg) / (static_cast<double>(truth.pt) * std::cosh(static_cast<double>(truth.eta)));
              const auto& state = tracks[i].innermost;
              out.culpritQOverPPull = static_cast<float>((static_cast<double>(state.qOverP()) - trueQOverP) / std::sqrt(static_cast<double>(state.covAt(kQOverP, kQOverP))));
            }
          }
        }
        continue;
      }
      double worstChi2 = -1.;
      for (std::size_t i = 0; i < fitted.size(); ++i) {
        if (trackChi2[i] > worstChi2) {
          worstChi2 = trackChi2[i];
          out.worstParticle = fitted[i]->particle;
        }
        if (fitted[i]->track) {
          out.contributors.push_back(fitted[i]->track->track);
        }
      }
      out.maxTrackChi2 = static_cast<float>(worstChi2);
    }
  });
  return vertices;
}

/// The collisions' fitted mocked primary vertices as reference points
/// (`referencePointsOf`).
struct ReferencePoints {
  std::vector<MeasuredPoint> points;
  std::unordered_map<std::uint64_t, std::int32_t> pointOf; ///< a collision's point, none without a fitted vertex
  /// \a collision's point, -1 for none.
  [[nodiscard]] std::int32_t of(std::uint64_t collision) const
  {
    const auto found = pointOf.find(collision);
    return found == pointOf.end() ? -1 : found->second;
  }
};

/// \a vertices' fitted ones as reference points.
[[nodiscard]] inline ReferencePoints referencePointsOf(std::span<const MockPrimaryVertex> vertices)
{
  ReferencePoints out;
  for (const auto& vertex : vertices) {
    if (vertex.fit.status == VertexFitStatus::kOk) {
      out.pointOf.emplace(vertex.collision, static_cast<std::int32_t>(out.points.size()));
      out.points.push_back({.position = vertex.fit.vertex, .covariance = vertex.fit.covariance});
    }
  }
  return out;
}

/// A timeframe's committed tracks associated to the mocked primary vertices
/// (`mockAssociation`), as a vertexer's association would give them: each
/// track's collision is its good particle's (a fake's majority particle's),
/// its flags say whether it entered its vertex's fit.
struct MockAssociation {
  std::vector<std::uint64_t> particle; ///< per track: its good particle, a fake's majority one, kNoTrackParticle on noise only
  std::vector<std::uint8_t> fake;      ///< per track: 1 without a good particle
  std::vector<std::uint8_t> flags;     ///< per track: TrackFlag bits
  ReferencePoints vertices;            ///< the collisions' fitted vertices
  /// Whether track \a t entered its vertex's fit.
  [[nodiscard]] bool contributor(std::uint32_t t) const { return (flags[t] & static_cast<std::uint8_t>(TrackFlag::kPrimaryVertexContributor)) != 0U; }
};

/// \a event's tracks associated to their collisions' mocked primary vertices
/// (`mockPrimaryVertices` of their physical primaries' good tracks): every
/// track the fit of its vertex took flagged `kPrimaryVertexContributor`.
template <DetectorLayout Det>
[[nodiscard]] MockAssociation mockAssociation(const TruthEvent<Det>& event)
{
  const auto& store = *event.tracks;
  MockAssociation out;
  const auto n = store.size();
  out.particle.assign(n, kNoTrackParticle);
  out.fake.assign(n, 0U);
  out.flags.assign(n, 0U);
  o2::gtrk::detail::parallelFor(*event.arena, 0U, n, 256, [&](std::uint32_t begin, std::uint32_t end) {
    for (std::uint32_t t = begin; t < end; ++t) {
      const auto rows = store.rows(t);
      const auto good = goodParticle(event.rows, rows);
      out.particle[t] = good != kNoTrackParticle ? good : majorityParticle(event.rows, rows).particle;
      out.fake[t] = good == kNoTrackParticle ? 1U : 0U;
    }
  });
  ChainSelection primaries;
  primaries.origins = 1U << static_cast<std::uint32_t>(ParticleOrigin::kPhysicalPrimary);
  ChainCounts unusedCounts;
  const auto chains = trackChains(event, store, primaries, unusedCounts);
  const auto vertices = mockPrimaryVertices(event, std::span<const Chain>{chains});
  out.vertices = referencePointsOf(vertices);
  for (const auto& vertex : vertices) {
    if (vertex.fit.status != VertexFitStatus::kOk) {
      continue;
    }
    for (const auto t : vertex.contributors) {
      out.flags[t] |= static_cast<std::uint8_t>(TrackFlag::kPrimaryVertexContributor);
    }
  }
  return out;
}

template <DetectorLayout Det>
class PrimaryVertexStudy
{
 public:
  explicit PrimaryVertexStudy(RecordWriter& writer) : mWriter(&writer) { writer.begin<PrimaryVertexRecord>(); }

  /// \a chains of \a event, from \a source.
  void process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source)
  {
    for (const auto& vertex : mockPrimaryVertices(event, chains)) {
      mWriter->write(recordOf(event, vertex, source));
      mFitted += vertex.fit.status == VertexFitStatus::kOk ? 1U : 0U;
      ++mCollisions;
    }
  }

  [[nodiscard]] std::vector<std::string> summary() const
  {
    return {std::format("primary vertex study: {} collisions with a primary's chain, {} fitted, records into {}", mCollisions, mFitted, mWriter->path())};
  }

 private:
  [[nodiscard]] static PrimaryVertexRecord recordOf(const TruthEvent<Det>& event, const MockPrimaryVertex& vertex, ChainSource source)
  {
    const auto& fit = vertex.fit;
    PrimaryVertexRecord record{.timeframe = event.timeframe, .source = static_cast<std::uint8_t>(source), .collision = vertex.collision, .tracks = vertex.tracks, .refused = vertex.refused, .refitFailed = vertex.refitFailed, .madeAway = vertex.madeAway};
    const auto mc = event.collisions->find(vertex.collision);
    if (mc != event.collisions->end()) {
      record.mcX = static_cast<float>(mc->second[0]);
      record.mcY = static_cast<float>(mc->second[1]);
      record.mcZ = static_cast<float>(mc->second[2]);
    }
    record.status = static_cast<std::uint8_t>(fit.status);
    if (vertex.tracks == 0U || fit.status == VertexFitStatus::kNoStart) {
      return record;
    }
    record.startDZ = static_cast<float>(1e4 * (vertex.start[2] - static_cast<double>(record.mcZ)));
    record.microseconds = vertex.microseconds;
    record.iterations = fit.iterations;
    record.outer = fit.outer;
    record.stopReason = fit.stopReason;
    if (fit.status != VertexFitStatus::kOk) {
      record.rescues = vertex.rescues;
      if (vertex.rescues > 0U) {
        const auto& culprit = event.particles->at(vertex.culprit);
        record.culpritPdg = culprit.pdg;
        record.culpritPt = culprit.pt;
        record.culpritEta = culprit.eta;
        record.culpritRows = vertex.culpritRows;
        record.culpritChains = vertex.culpritChains;
        record.culpritQOverPPull = vertex.culpritQOverPPull;
      }
      return record;
    }
    record.chi2 = static_cast<float>(fit.chi2);
    record.priorChi2 = static_cast<float>(fit.priorChi2);
    record.ndf = static_cast<std::int16_t>(fit.ndf);
    const std::array<float*, 3> d{&record.dX, &record.dY, &record.dZ};
    const std::array<float*, 3> sigma{&record.sX, &record.sY, &record.sZ};
    const std::array<double, 3> truthPoint{record.mcX, record.mcY, record.mcZ};
    for (std::uint32_t k = 0; k < 3; ++k) {
      *d[k] = static_cast<float>(1e4 * (fit.vertex[k] - truthPoint[k]));
      *sigma[k] = static_cast<float>(1e4 * std::sqrt(fit.covariance[decay_fit_detail::packedIndex(k, k)]));
    }
    const auto rho = [&](std::uint32_t a, std::uint32_t b) { return static_cast<float>(fit.covariance[decay_fit_detail::packedIndex(a, b)] / std::sqrt(fit.covariance[decay_fit_detail::packedIndex(a, a)] * fit.covariance[decay_fit_detail::packedIndex(b, b)])); };
    record.rhoXY = rho(0, 1);
    record.rhoXZ = rho(0, 2);
    record.rhoYZ = rho(1, 2);
    const auto& worstTruth = event.particles->at(vertex.worstParticle);
    record.maxTrackChi2 = vertex.maxTrackChi2;
    record.worstPdg = worstTruth.pdg;
    record.worstPt = worstTruth.pt;
    return record;
  }

  RecordWriter* mWriter = nullptr;
  std::uint64_t mCollisions = 0U;
  std::uint64_t mFitted = 0U;
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_PRIMARYVERTEXSTUDY_H_
