// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DecayFinderStudy.h
/// \brief A V0 finder on the tracker's tracks with the primary vertex and the
///        association mocked from the MC, held to the MC decays.
///
/// The tracks: every committed track. Its association: the collision of its
/// good particle (the one on all its rows), a fake's of its majority
/// particle; a track of noise only is not combined. Each collision's point:
/// its mocked primary vertex (`mockPrimaryVertices`: a fit of the primaries
/// made at the MC vertex), a track of that fit tested against the vertex
/// without it.
///
/// No PID: every associated track is refitted as a pion and as a proton
/// (DaughterTable.h). Each entry's displacement is its impact parameters'
/// significance against its collision's point at its perigee,
/// S^2 = r^T (V + J C J^T)^-1 r: r the transverse and longitudinal
/// impact parameters, V their covariance, C the point's, J their derivative
/// by the point (-u across, cot(theta) t - z along). Every opposite-charge
/// pair of a collision under every V0 hypothesis (K0S -> pi+ pi-, Lambda ->
/// p pi-, anti-Lambda -> anti-p pi+) is a candidate: nothing is cut, every
/// one is fitted and written, its legs' S^2, the fit and its observables
/// beside it for the analysis to select on. Each candidate starts at its
/// helices' closest approach and is fitted (DecayCandidateFit.h); a
/// candidate whose leg's refit failed is written with that status, unfitted.
/// The collisions go in blocks, each block's records written before the
/// next is fitted, so the memory is a block's.
///
/// Records: one per candidate (`DecayFinderCandidateRecord`): its MC class,
/// whether its species are its tracks' true ones (`pidMatch`: perfect PID
/// keeps those alone) and how many of its tracks are primaries (a perfect
/// rejection of primaries keeps none), its legs' S^2, the start's and the
/// legs' innermost radii, the fit and its observables. One per MC V0 decaying into its two charged daughters
/// (`DecayFinderDecayRecord`, a second file beside the first, `-decays`):
/// the first stage that lost it, its legs' S^2 under their true species,
/// its own candidate's fit.

#ifndef GTRK_STUDY_DECAYFINDERSTUDY_H_
#define GTRK_STUDY_DECAYFINDERSTUDY_H_

#include "GenericTracking/Tracking/DaughterTable.h"
#include "GenericTracking/Tracking/DecayCandidateFit.h"
#include "GenericTracking/Tracking/HypothesisRefit.h"
#include "GenericTracking/Study/DecayCandidates.h"
#include "GenericTracking/Study/MCTruth.h"
#include "GenericTracking/Study/PrimaryVertexStudy.h"
#include "GenericTracking/Study/RecordWriter.h"
#include "GenericTracking/Study/StudyChains.h"
#include "GenericTracking/Study/TruthEvent.h"
#include "GenericTracking/Study/TruthFit.h"

#include "ReconstructionDataFormats/PID.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <format>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace o2::gtrk::study
{

/// A candidate's MC class (`DecayFinderCandidateRecord::truthClass`), one of:
enum class FinderClass : std::uint8_t {
  kTrue = 1,        ///< an MC V0's own pair under its own hypothesis, charges right
  kOtherHypothesis, ///< an MC V0's own pair under another hypothesis (a Lambda fitted as K0S)
  kWrongCharge,     ///< an MC V0's own pair, a track's reconstructed charge wrong
  kOtherDecay,      ///< both tracks from one mother that is no V0 of the hypotheses (a rho0, a conversion)
  kOneLeg,          ///< one track from an MC V0, the other not
  kCombinatorial,   ///< no common mother, every track with a good particle
  kFake             ///< a track without a good particle
};

/// Where an MC V0 was lost (`DecayFinderDecayRecord::stage`), the first of:
enum class FinderStage : std::uint8_t {
  kFound = 0,   ///< its own candidate was fitted
  kNotFound,    ///< a daughter has no good track
  kRefitFailed, ///< a daughter's refit under its true species failed
  kNoPoint,     ///< its collision has no mocked primary vertex
  kCarryFailed, ///< a daughter's carry to the point failed
  kWrongCharge, ///< both tracks of one reconstructed charge
  kNoStart,     ///< the helices give no start
  kFitFailed    ///< the fit failed
};

/// One candidate (`DecayFinderStudy`).
struct DecayFinderCandidateRecord {
  static constexpr std::string_view kName = "decayFinderCandidate";
  std::uint32_t timeframe = 0U;
  std::uint64_t collision = 0U;
  std::uint8_t hypothesis = 0U; ///< 0 K0S, 1 Lambda, 2 anti-Lambda
  std::uint8_t truthClass = 0U; ///< FinderClass
  std::uint8_t pidMatch = 0U;   ///< every slot's species its track's true one
  std::uint8_t primaries = 0U;  ///< its tracks whose particle is a physical primary
  std::uint64_t mother = 0U;    ///< the common MC mother, 0 without one
  std::int32_t motherPdg = 0;
  float s2Positive = 0.F; ///< the legs' displacement S^2
  float s2Negative = 0.F;
  float ptPositive = 0.F; ///< the legs' entries' pT [GeV]
  float ptNegative = 0.F;
  float firstPositive = 0.F; ///< the legs' innermost measurements' transverse radii [cm]
  float firstNegative = 0.F;
  float startRadius = 0.F;  ///< the start's transverse radius [cm], 0 without one
  std::uint8_t status = 0U; ///< VertexFitStatus
  std::uint8_t iterations = 0U;
  float chi2 = 0.F;
  std::int16_t ndf = 0;
  float radius = 0.F; ///< the vertex's transverse radius [cm]
  float mass = 0.F;
  float sMass = 0.F;
  float length = 0.F; ///< the flight from the point [cm]
  float sLength = 0.F;
  float cosPointing = 0.F;
  float microseconds = 0.F;

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("collision", collision);
    visit("hypothesis", hypothesis);
    visit("truthClass", truthClass);
    visit("pidMatch", pidMatch);
    visit("primaries", primaries);
    visit("mother", mother);
    visit("motherPdg", motherPdg);
    visit("s2Positive", s2Positive);
    visit("s2Negative", s2Negative);
    visit("ptPositive", ptPositive);
    visit("ptNegative", ptNegative);
    visit("firstPositive", firstPositive);
    visit("firstNegative", firstNegative);
    visit("startRadius", startRadius);
    visit("status", status);
    visit("iterations", iterations);
    visit("chi2", chi2);
    visit("ndf", ndf);
    visit("radius", radius);
    visit("mass", mass);
    visit("sMass", sMass);
    visit("length", length);
    visit("sLength", sLength);
    visit("cosPointing", cosPointing);
    visit("microseconds", microseconds);
  }
};

/// One MC V0 (`DecayFinderStudy`).
struct DecayFinderDecayRecord {
  static constexpr std::string_view kName = "decayFinderDecay";
  std::uint32_t timeframe = 0U;
  std::uint64_t particle = 0U;
  std::int32_t pdg = 0;
  float decayR = 0.F;     ///< the MC decay point's transverse radius [cm]
  float pt = 0.F;         ///< the mother's at production [GeV]
  float ptPositive = 0.F; ///< the daughters' at production [GeV]
  float ptNegative = 0.F;
  std::uint8_t stage = 0U; ///< FinderStage
  float s2Positive = std::numeric_limits<float>::quiet_NaN();
  float s2Negative = std::numeric_limits<float>::quiet_NaN();
  std::uint8_t status = 0U; ///< its own candidate's VertexFitStatus, 0 without one
  float chi2 = 0.F;
  std::int16_t ndf = 0;
  float mass = 0.F;
  float sMass = 0.F;
  float length = 0.F;
  float sLength = 0.F;
  float cosPointing = 0.F;
  float dXY = 0.F; ///< the fitted vertex's transverse distance from the MC decay point [cm]
  float dZ = 0.F;

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("particle", particle);
    visit("pdg", pdg);
    visit("decayR", decayR);
    visit("pt", pt);
    visit("ptPositive", ptPositive);
    visit("ptNegative", ptNegative);
    visit("stage", stage);
    visit("s2Positive", s2Positive);
    visit("s2Negative", s2Negative);
    visit("status", status);
    visit("chi2", chi2);
    visit("ndf", ndf);
    visit("mass", mass);
    visit("sMass", sMass);
    visit("length", length);
    visit("sLength", sLength);
    visit("cosPointing", cosPointing);
    visit("dXY", dXY);
    visit("dZ", dZ);
  }
};

namespace decay_finder_detail
{

/// The collisions fitted and written as one block (the file's comment):
/// what bounds the records held at once.
inline constexpr std::size_t kCollisionsPerBlock = 256U;

/// \a path with \a tag before its extension.
[[nodiscard]] inline std::string taggedPath(const std::string& path, std::string_view tag)
{
  const auto dot = path.rfind('.');
  return path.substr(0, dot) + std::string{tag} + path.substr(dot);
}

/// \a state's displacement S^2 against \a point (`DecayFinderStudy`'s
/// comment), its perigee about the point's transverse position; NaN if the
/// carry fails or the matrix is singular.
template <DetectorLayout Det>
[[nodiscard]] double displacement(const Propagator& base, const TrackingGeometry<Det>& geometry, const TrackingSurface& surface, const BoundTrackState& state, o2::track::PID pid, const ReferencePoint& point)
{
  HypothesisFit<> fit{.pid = pid};
  fit.innermost = state;
  if (carryToPerigee(base, geometry, surface, fit, static_cast<Scalar>(point.position[0]), static_cast<Scalar>(point.position[1])) != PropagationStatus::kOk) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const auto& perigee = fit.perigee;
  const std::array<double, 2> r{static_cast<double>(perigee.loc0()), static_cast<double>(perigee.loc1()) - point.position[2]};
  // The loc0 axis in global coordinates, from the plane itself.
  std::array<Scalar, 3> origin{};
  std::array<Scalar, 3> across{};
  surfaceToGlobal(fit.perigeePlane, 0.F, 0.F, origin.data());
  surfaceToGlobal(fit.perigeePlane, 1.F, 0.F, across.data());
  const auto phi = static_cast<double>(perigee.phi());
  const double cotTheta = 1. / std::tan(static_cast<double>(perigee.theta()));
  const std::array<std::array<double, 3>, 2> j{{{-(across[0] - origin[0]), -(across[1] - origin[1]), -(across[2] - origin[2])}, {cotTheta * std::cos(phi), cotTheta * std::sin(phi), -1.}}};
  std::array<double, 3> m{static_cast<double>(perigee.covAt(kLoc0, kLoc0)), static_cast<double>(perigee.covAt(kLoc0, kLoc1)), static_cast<double>(perigee.covAt(kLoc1, kLoc1))};
  const auto projected = [&](int a, int b) {
    double sum = 0.;
    for (int k = 0; k < 3; ++k) {
      for (int l = 0; l < 3; ++l) {
        sum += j[a][k] * point.covariance[covIndex(k, l)] * j[b][l];
      }
    }
    return sum;
  };
  m[0] += projected(0, 0);
  m[1] += projected(0, 1);
  m[2] += projected(1, 1);
  const double det = (m[0] * m[2]) - (m[1] * m[1]);
  if (!(det > 0.) || !(m[0] > 0.)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return ((r[0] * r[0] * m[2]) - (2. * r[0] * r[1] * m[1]) + (r[1] * r[1] * m[0])) / det;
}

} // namespace decay_finder_detail

template <DetectorLayout Det>
class DecayFinderStudy
{
 public:
  explicit DecayFinderStudy(RecordWriter& writer) : mWriter(&writer), mDecayWriter(std::make_unique<RecordWriter>(decay_finder_detail::taggedPath(writer.path(), "-decays")))
  {
    writer.begin<DecayFinderCandidateRecord>();
    mDecayWriter->begin<DecayFinderDecayRecord>();
  }

  /// \a event's tracks, every one, whatever the source's chains.
  void process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source);

  [[nodiscard]] std::vector<std::string> summary() const
  {
    return {std::format("decay finder study: {} tracks combined ({} without a particle), {} candidates, {} fitted (records into {}), {} MC V0s ({} found) into {}; the daughter table {} entries in {:.3f} s, the S^2 in {:.3f} s, the vertices in {:.3f} s, the fits {:.1f} thread-s",
                        mTracks, mUnassociated, mPairs, mFitted, mWriter->path(), mDecays, mFound, mDecayWriter->path(), mEntries, mFlushSeconds, mDisplacementSeconds, mVertexSeconds, mFitMicroseconds * 1e-6)};
  }

 private:
  RecordWriter* mWriter = nullptr;
  std::unique_ptr<RecordWriter> mDecayWriter;
  std::size_t mTracks = 0U;
  std::size_t mUnassociated = 0U;
  std::size_t mPairs = 0U;
  std::size_t mFitted = 0U;
  std::size_t mDecays = 0U;
  std::size_t mFound = 0U;
  std::size_t mEntries = 0U;
  double mFlushSeconds = 0.;
  double mDisplacementSeconds = 0.;
  double mVertexSeconds = 0.;
  double mFitMicroseconds = 0.;
};

template <DetectorLayout Det>
void DecayFinderStudy<Det>::process(const TruthEvent<Det>& event, std::span<const Chain> /*chains*/, ChainSource source)
{
  using namespace decay_finder_detail;
  using o2::track::PID;
  if (source != ChainSource::kTracks || event.tracks == nullptr || event.trees == nullptr) {
    throw std::runtime_error("decay finder study: it takes the tracker's tracks and the decay trees");
  }
  const auto& store = *event.tracks;
  const auto& geometry = *event.geometry;
  const auto& frame = *event.frame;
  const auto& trees = *event.trees;
  const auto committed = store.tracks().view();

  const std::array<PID, 2> species{PID{PID::Pion}, PID{PID::Proton}};
  constexpr std::array<std::int32_t, 3> kMothers{310, 3122, -3122};
  std::array<DecayHypothesis, 3> hypotheses{};
  const auto define = [&](DecayHypothesis& h, PID positive, PID negative, double mass) {
    h.species.push_back(positive);
    h.charges.push_back(1);
    h.species.push_back(negative);
    h.charges.push_back(-1);
    h.massGeV = mass;
  };
  define(hypotheses[0], PID{PID::Pion}, PID{PID::Pion}, 0.497611);
  define(hypotheses[1], PID{PID::Proton}, PID{PID::Pion}, 1.115683);
  define(hypotheses[2], PID{PID::Pion}, PID{PID::Proton}, 1.115683);
  const auto forEach = [&](std::uint32_t n, const auto& body) {
    o2::gtrk::detail::parallelFor(*event.arena, 0U, n, 16, [&](std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; ++i) {
        body(static_cast<std::uint32_t>(i));
      }
    });
  };

  // The mocked primary vertices, and the vertex without each of their tracks.
  const auto vertexStart = std::chrono::steady_clock::now();
  ChainSelection primaries;
  primaries.origins = 1U << static_cast<unsigned>(ParticleOrigin::kPhysicalPrimary);
  ChainCounts unusedCounts;
  const auto primaryChains = trackChains(event, store, primaries, unusedCounts);
  std::vector<ReferencePoint> points;
  std::unordered_map<std::uint64_t, std::int32_t> pointOf;
  std::unordered_map<std::uint32_t, ReferencePoint> withoutTrack;
  for (auto& vertex : mockPrimaryVertices(event, std::span<const Chain>{primaryChains})) {
    if (vertex.fit.status != VertexFitStatus::kOk) {
      continue;
    }
    pointOf.emplace(vertex.collision, static_cast<std::int32_t>(points.size()));
    points.push_back({.position = vertex.fit.vertex, .covariance = vertex.fit.covariance});
    for (std::size_t i = 0; i < vertex.trackIds.size(); ++i) {
      withoutTrack.emplace(vertex.trackIds[i], vertex.without[i]);
    }
  }
  mVertexSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - vertexStart).count();

  // Every committed track's particle (its good one, a fake's majority one)
  // and its collision.
  const auto n = static_cast<std::uint32_t>(store.size());
  std::vector<std::uint64_t> particleOf(n, kNoTrackParticle);
  std::vector<std::uint8_t> fake(n, 0U);
  forEach(n, [&](std::uint32_t t) {
    const auto rows = store.rows(t);
    const auto good = goodParticle(event.rows, rows);
    particleOf[t] = good != kNoTrackParticle ? good : majorityParticle(event.rows, rows).particle;
    fake[t] = good == kNoTrackParticle ? 1U : 0U;
  });
  std::unordered_map<std::uint64_t, std::array<std::vector<std::uint32_t>, 2>> byCollision; // positive, negative
  std::vector<std::uint64_t> collisions;
  for (std::uint32_t t = 0; t < n; ++t) {
    if (particleOf[t] == kNoTrackParticle) {
      ++mUnassociated;
      continue;
    }
    const auto collision = collisionOf(particleOf[t]);
    auto [at, inserted] = byCollision.try_emplace(collision);
    if (inserted) {
      collisions.push_back(collision);
    }
    at->second[committed.state[t].qOverP() > 0.F ? 0 : 1].push_back(t);
    ++mTracks;
  }

  // The table: every combined track as a pion and as a proton.
  DaughterTable<Det> table(std::span<const PID>{species});
  if (table.begin(n) != StageStatus::kComplete) {
    throw std::runtime_error("decay finder study: no daughter table");
  }
  for (const auto& [collision, lists] : byCollision) {
    for (const auto& list : lists) {
      for (const auto t : list) {
        for (std::size_t s = 0; s < species.size(); ++s) {
          table.request(t, static_cast<std::uint8_t>(s));
        }
      }
    }
  }
  const auto flushStart = std::chrono::steady_clock::now();
  if (table.flush(forEach, 4096U, event.propagator, geometry, frame, store) != StageStatus::kComplete) {
    throw std::runtime_error("decay finder study: the daughter table's flush failed");
  }
  mFlushSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - flushStart).count();
  mEntries += table.size();

  // Each entry's displacement against its collision's point (a primary
  // vertex's own track against the vertex without it).
  const auto entries = table.view();
  const auto displacementStart = std::chrono::steady_clock::now();
  std::vector<float> s2(table.size(), std::numeric_limits<float>::quiet_NaN());
  std::vector<std::uint8_t> noPoint(table.size(), 0U);
  forEach(static_cast<std::uint32_t>(table.size()), [&](std::uint32_t e) {
    const auto t = entries.track[e];
    if (entries.status[e] != PropagationStatus::kOk) {
      return;
    }
    const auto own = withoutTrack.find(t);
    const auto shared = pointOf.find(collisionOf(particleOf[t]));
    if (own == withoutTrack.end() && shared == pointOf.end()) {
      noPoint[e] = 1U;
      return;
    }
    const auto& point = own != withoutTrack.end() ? own->second : points[static_cast<std::size_t>(shared->second)];
    const auto& surface = geometry.planeOf(measurementAt(frame, store.rows(t).front()));
    s2[e] = static_cast<float>(displacement(event.propagator, geometry, surface, entries.state[e], table.species()[entries.species[e]], point));
  });
  mDisplacementSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - displacementStart).count();

  // The MC V0s decaying into their two charged daughters, by daughter.
  struct Decay {
    const TreeParticle* mother = nullptr;
    std::array<const TreeParticle*, 2> legs{}; ///< the positive, the negative daughter
    int hypothesis = 0;
    DecayFinderDecayRecord record;
    bool fitted = false;
  };
  std::vector<Decay> decays;
  std::unordered_map<std::uint64_t, std::size_t> decayOfLeg;
  for (const auto& member : trees.particles) {
    const auto h = static_cast<int>(std::ranges::find(kMothers, member.pdg) - kMothers.begin());
    if (h >= static_cast<int>(kMothers.size())) {
      continue;
    }
    const auto products = decayProducts(trees, member);
    if (products.size() != 2U) {
      continue;
    }
    const bool firstPositive = particleCharge(products[0]->pdg) > 0.;
    const std::array<const TreeParticle*, 2> legs{firstPositive ? products[0] : products[1], firstPositive ? products[1] : products[0]};
    bool matches = particleCharge(legs[0]->pdg) > 0. && particleCharge(legs[1]->pdg) < 0.;
    for (std::size_t slot = 0; slot < 2U && matches; ++slot) {
      const auto truth = truthSpecies(legs[slot]->pdg);
      matches = truth && truth->pid == hypotheses[static_cast<std::size_t>(h)].species[slot];
    }
    if (!matches) {
      continue;
    }
    const auto& decayPoint = products.front()->vertex;
    const auto& p = member.momentum;
    decayOfLeg.emplace(legs[0]->id, decays.size());
    decayOfLeg.emplace(legs[1]->id, decays.size());
    decays.push_back({.mother = &member, .legs = legs, .hypothesis = h, .record = {.timeframe = event.timeframe, .particle = member.id, .pdg = member.pdg, .decayR = static_cast<float>(std::hypot(decayPoint[0], decayPoint[1])), .pt = static_cast<float>(std::hypot(p[0], p[1])), .ptPositive = static_cast<float>(std::hypot(legs[0]->momentum[0], legs[0]->momentum[1])), .ptNegative = static_cast<float>(std::hypot(legs[1]->momentum[0], legs[1]->momentum[1]))}});
  }
  // Each MC daughter's longest good track.
  std::unordered_map<std::uint64_t, std::uint32_t> trackOfLeg;
  for (std::uint32_t t = 0; t < n; ++t) {
    if (fake[t] != 0U || !decayOfLeg.contains(particleOf[t])) {
      continue;
    }
    const auto [at, inserted] = trackOfLeg.try_emplace(particleOf[t], t);
    if (!inserted && store.rows(t).size() > store.rows(at->second).size()) {
      at->second = t;
    }
  }

  // The pairs, collision by collision; every candidate fitted on its own.
  const auto motherOf = [&](std::uint32_t t) {
    const auto* particle = fake[t] != 0U ? nullptr : trees.find(particleOf[t]);
    return particle == nullptr ? kNoTreeParticle : particle->mother;
  };
  const auto isPrimary = [&](std::uint32_t t) {
    const auto found = event.particles->find(particleOf[t]);
    return fake[t] == 0U && found != event.particles->end() && found->second.origin == ParticleOrigin::kPhysicalPrimary;
  };
  const auto one = [](std::uint32_t count, const auto& body) {
    for (std::uint32_t i = 0; i < count; ++i) {
      body(i);
    }
  };
  const auto radiusOf = [&](std::int32_t e) {
    const auto& state = entries.state[static_cast<std::size_t>(e)];
    std::array<Scalar, 3> at{};
    surfaceToGlobal(geometry.planeOf(state), state.loc0(), state.loc1(), at.data());
    return static_cast<float>(std::hypot(at[0], at[1]));
  };
  std::vector<std::vector<DecayFinderCandidateRecord>> perCollision(std::min(collisions.size(), kCollisionsPerBlock));
  for (std::size_t first = 0; first < collisions.size(); first += kCollisionsPerBlock) {
    const std::size_t last = std::min(collisions.size(), first + kCollisionsPerBlock);
    o2::gtrk::detail::parallelFor(*event.arena, first, last, 1, [&](std::size_t begin, std::size_t end) {
      for (std::size_t c = begin; c < end; ++c) {
        const auto collision = collisions[c];
        const auto& lists = byCollision.at(collision);
        const auto point = pointOf.find(collision);
        auto& out = perCollision[c - first];
        out.clear();
        for (const auto positive : lists[0]) {
          for (const auto negative : lists[1]) {
            const auto mother = motherOf(positive);
            const bool sameMother = mother != kNoTreeParticle && mother == motherOf(negative);
            const auto ownDecay = sameMother ? decayOfLeg.find(particleOf[positive]) : decayOfLeg.end();
            const bool ownPair = ownDecay != decayOfLeg.end() && trackOfLeg.at(particleOf[positive]) == positive && trackOfLeg.contains(particleOf[negative]) && trackOfLeg.at(particleOf[negative]) == negative;
            for (std::size_t h = 0; h < hypotheses.size(); ++h) {
              const auto& hypothesis = hypotheses[h];
              const auto ePositive = table.entryOf(positive, static_cast<std::uint8_t>(table.speciesIndex(hypothesis.species[0])));
              const auto eNegative = table.entryOf(negative, static_cast<std::uint8_t>(table.speciesIndex(hypothesis.species[1])));
              const float s2Positive = s2[static_cast<std::size_t>(ePositive)];
              const float s2Negative = s2[static_cast<std::size_t>(eNegative)];
              const bool own = ownPair && decays[ownDecay->second].hypothesis == static_cast<int>(h);
              DecayFinderCandidateRecord record{.timeframe = event.timeframe, .collision = collision, .hypothesis = static_cast<std::uint8_t>(h), .s2Positive = s2Positive, .s2Negative = s2Negative};
              record.primaries = static_cast<std::uint8_t>((isPrimary(positive) ? 1 : 0) + (isPrimary(negative) ? 1 : 0));
              const auto truthPositive = truthSpecies(event.particles->at(particleOf[positive]).pdg);
              const auto truthNegative = truthSpecies(event.particles->at(particleOf[negative]).pdg);
              record.pidMatch = static_cast<std::uint8_t>(truthPositive && truthNegative && truthPositive->pid == hypothesis.species[0] && truthNegative->pid == hypothesis.species[1] ? 1U : 0U);
              if (fake[positive] != 0U || fake[negative] != 0U) {
                record.truthClass = static_cast<std::uint8_t>(FinderClass::kFake);
              } else if (sameMother) {
                const auto* m = trees.find(mother);
                record.mother = mother;
                record.motherPdg = m != nullptr ? m->pdg : 0;
                if (ownDecay != decayOfLeg.end()) {
                  const bool chargesRight = particleCharge(event.particles->at(particleOf[positive]).pdg) > 0. && particleCharge(event.particles->at(particleOf[negative]).pdg) < 0.;
                  record.truthClass = static_cast<std::uint8_t>(!chargesRight ? FinderClass::kWrongCharge : (decays[ownDecay->second].hypothesis == static_cast<int>(h) ? FinderClass::kTrue : FinderClass::kOtherHypothesis));
                } else {
                  record.truthClass = static_cast<std::uint8_t>(FinderClass::kOtherDecay);
                }
              } else {
                const bool positiveLeg = decayOfLeg.contains(particleOf[positive]);
                const bool negativeLeg = decayOfLeg.contains(particleOf[negative]);
                record.truthClass = static_cast<std::uint8_t>(positiveLeg != negativeLeg ? FinderClass::kOneLeg : FinderClass::kCombinatorial);
              }
              const auto ptOf = [&](std::int32_t e) { const auto& state = entries.state[static_cast<std::size_t>(e)]; return static_cast<float>(std::abs(std::sin(state.theta()) / state.qOverP())); };
              const auto refitted = [&](std::int32_t e) { return entries.status[static_cast<std::size_t>(e)] == PropagationStatus::kOk; };
              const auto nan = std::numeric_limits<float>::quiet_NaN();
              record.ptPositive = refitted(ePositive) ? ptOf(ePositive) : nan;
              record.ptNegative = refitted(eNegative) ? ptOf(eNegative) : nan;
              record.firstPositive = refitted(ePositive) ? radiusOf(ePositive) : nan;
              record.firstNegative = refitted(eNegative) ? radiusOf(eNegative) : nan;
              DecayCandidate candidate{.hypothesis = static_cast<std::uint16_t>(h), .point = point == pointOf.end() ? -1 : point->second};
              candidate.entries.push_back(static_cast<std::uint32_t>(ePositive));
              candidate.entries.push_back(static_cast<std::uint32_t>(eNegative));
              DecayCandidateResult result;
              const auto started = std::chrono::steady_clock::now();
              if (entries.status[static_cast<std::size_t>(ePositive)] == PropagationStatus::kOk && entries.status[static_cast<std::size_t>(eNegative)] == PropagationStatus::kOk) {
                const auto daughters = candidateDaughters(candidate, table, store, geometry, frame, event.propagator);
                std::array<Vec3D, decay_fit_detail::kMaxConstrainedVertices> start{};
                if (decayFitStart<Det>(std::span<const DecayDaughter>{daughters.data(), daughters.size()}, {}, {}, {}, start)) {
                  candidate.start.push_back(start[0]);
                  record.startRadius = static_cast<float>(std::hypot(start[0][0], start[0][1]));
                }
              }
              fitDecayCandidates<Det>(one, std::span<const DecayCandidate>{&candidate, 1}, std::span<const DecayHypothesis>{hypotheses}, table, store, frame, geometry, event.propagator, geometry.material(), std::span<const ReferencePoint>{points}, std::span<DecayCandidateResult>{&result, 1});
              record.microseconds = static_cast<float>(std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - started).count());
              record.status = static_cast<std::uint8_t>(result.fit.status);
              record.iterations = result.fit.iterations;
              if (result.fit.status == VertexFitStatus::kOk) {
                record.chi2 = static_cast<float>(result.fit.chi2);
                record.ndf = static_cast<std::int16_t>(result.fit.ndf);
                record.radius = static_cast<float>(std::hypot(result.fit.vertex[0], result.fit.vertex[1]));
                record.mass = static_cast<float>(result.kinematics.mass);
                record.sMass = static_cast<float>(result.kinematics.sMass);
                if (result.hasPoint) {
                  record.length = static_cast<float>(result.flight.length);
                  record.sLength = static_cast<float>(result.flight.sLength);
                  record.cosPointing = static_cast<float>(result.flight.cosPointing);
                }
              }
              if (own) {
                // The decay's own candidate: its record's fit (one writer per decay).
                auto& decay = decays[ownDecay->second];
                auto& d = decay.record;
                decay.fitted = true;
                d.status = record.status;
                d.s2Positive = s2Positive;
                d.s2Negative = s2Negative;
                if (result.fit.status == VertexFitStatus::kOk) {
                  const auto& decayPoint = decayProducts(trees, *decay.mother).front()->vertex;
                  d.chi2 = record.chi2;
                  d.ndf = record.ndf;
                  d.mass = record.mass;
                  d.sMass = record.sMass;
                  d.length = record.length;
                  d.sLength = record.sLength;
                  d.cosPointing = record.cosPointing;
                  d.dXY = static_cast<float>(std::hypot(result.fit.vertex[0] - decayPoint[0], result.fit.vertex[1] - decayPoint[1]));
                  d.dZ = static_cast<float>(result.fit.vertex[2] - decayPoint[2]);
                }
              }
              out.push_back(record);
            }
          }
        }
      }
    });
    for (std::size_t c = first; c < last; ++c) {
      for (const auto& record : perCollision[c - first]) {
        mWriter->write(record);
        mFitted += record.status == static_cast<std::uint8_t>(VertexFitStatus::kOk) ? 1U : 0U;
        mFitMicroseconds += record.microseconds;
      }
      mPairs += perCollision[c - first].size();
    }
  }

  // Each MC V0's stage.
  for (auto& decay : decays) {
    auto& d = decay.record;
    const auto positive = trackOfLeg.find(decay.legs[0]->id);
    const auto negative = trackOfLeg.find(decay.legs[1]->id);
    const auto stageOf = [&]() {
      if (positive == trackOfLeg.end() || negative == trackOfLeg.end()) {
        return FinderStage::kNotFound;
      }
      const auto& hypothesis = hypotheses[static_cast<std::size_t>(decay.hypothesis)];
      const auto ePositive = table.entryOf(positive->second, static_cast<std::uint8_t>(table.speciesIndex(hypothesis.species[0])));
      const auto eNegative = table.entryOf(negative->second, static_cast<std::uint8_t>(table.speciesIndex(hypothesis.species[1])));
      if (ePositive < 0 || eNegative < 0 || entries.status[static_cast<std::size_t>(ePositive)] != PropagationStatus::kOk || entries.status[static_cast<std::size_t>(eNegative)] != PropagationStatus::kOk) {
        return FinderStage::kRefitFailed;
      }
      if (noPoint[static_cast<std::size_t>(ePositive)] != 0U || noPoint[static_cast<std::size_t>(eNegative)] != 0U) {
        return FinderStage::kNoPoint;
      }
      d.s2Positive = s2[static_cast<std::size_t>(ePositive)];
      d.s2Negative = s2[static_cast<std::size_t>(eNegative)];
      if (std::isnan(d.s2Positive) || std::isnan(d.s2Negative)) {
        return FinderStage::kCarryFailed;
      }
      if (!decay.fitted) {
        return FinderStage::kWrongCharge; // both its tracks of one reconstructed charge: never paired
      }
      if (d.status == static_cast<std::uint8_t>(VertexFitStatus::kNoStart)) {
        return FinderStage::kNoStart;
      }
      return d.status == static_cast<std::uint8_t>(VertexFitStatus::kOk) ? FinderStage::kFound : FinderStage::kFitFailed;
    };
    d.stage = static_cast<std::uint8_t>(stageOf());
    mFound += d.stage == static_cast<std::uint8_t>(FinderStage::kFound) ? 1U : 0U;
    mDecayWriter->write(d);
  }
  mDecays += decays.size();
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYFINDERSTUDY_H_
