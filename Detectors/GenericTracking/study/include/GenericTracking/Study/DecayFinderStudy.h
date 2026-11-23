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
/// The tracks: every committed track associated to its collision's mocked
/// primary vertex (`mockAssociation`: the collision of its good particle, a
/// fake's of its majority particle; the vertex a fit of the primaries made
/// at the MC vertex). A track flagged as its vertex's contributor
/// (`TrackFlag::kPrimaryVertexContributor`) is not combined, nor one of noise
/// only.
///
/// No PID: every associated track is refitted as a pion, a proton and a kaon
/// (DaughterTable.h). Each entry's displacement is its impact parameters'
/// significance S^2 against its collision's point at its perigee
/// (`displacement`). Every opposite-charge
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
/// legs' innermost radii, the helices' gaps at the start, the fit and its
/// observables. One per MC V0 decaying into its two charged daughters
/// (`DecayFinderDecayRecord`, a second file beside the first, `-decays`):
/// the first stage that lost it, its legs' S^2 under their true species,
/// its own candidate's fit.
///
/// Cascades (CascadeFinder.h): the Lambda candidates of the pairs that
/// `DecayParam` (GTRKDecay) selects, each with every bachelor of its
/// collision, the tracks refitted as kaons too; records per candidate
/// (`-cascades`) and per MC cascade (`-cascadeDecays`).

#ifndef GTRK_STUDY_DECAYFINDERSTUDY_H_
#define GTRK_STUDY_DECAYFINDERSTUDY_H_

#include "GenericTracking/Tracking/DaughterTable.h"
#include "GenericTracking/Tracking/DecayCandidateFit.h"
#include "GenericTracking/Study/CascadeFinder.h"
#include "GenericTracking/Study/DecayCandidates.h"
#include "GenericTracking/Study/DecayParam.h"
#include "GenericTracking/Study/DecayStudyTools.h"
#include "GenericTracking/Study/MCTruth.h"
#include "GenericTracking/Study/MotherSearch.h"
#include "GenericTracking/Study/PrimaryVertexStudy.h"
#include "GenericTracking/Study/RecordWriter.h"
#include "GenericTracking/Study/StateTruth.h"
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
#include <unordered_set>
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
  kContributor, ///< a daughter's track entered its primary vertex's fit: never combined
  kRefitFailed, ///< a daughter's refit under its true species failed
  kNoPoint,     ///< its collision has no mocked primary vertex
  kCarryFailed, ///< a daughter's carry to the point failed
  kWrongCharge, ///< a daughter's good tracks all of the wrong reconstructed charge (a looper's way back, reconstructed outgoing); taken with kNotFound
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
  float startRadius = 0.F; ///< the start's transverse radius [cm], 0 without one
  /// The legs' helices at the start (`approachOf`): the gap between their
  /// transverse circles (0 where they cross) and between their z [cm]; NaN
  /// without a start.
  float gapXY = std::numeric_limits<float>::quiet_NaN();
  float gapZ = std::numeric_limits<float>::quiet_NaN();
  std::uint8_t status = 0U; ///< VertexFitStatus
  std::uint8_t iterations = 0U;
  float chi2 = 0.F;
  std::int16_t ndf = 0;
  float radius = 0.F; ///< the vertex's transverse radius [cm]
  float mass = 0.F;
  float sMass = 0.F;
  float length = std::numeric_limits<float>::quiet_NaN(); ///< the flight from the point [cm]; NaN without a point
  float sLength = std::numeric_limits<float>::quiet_NaN();
  float cosPointing = std::numeric_limits<float>::quiet_NaN();
  std::uint8_t ahead = 0U; ///< its legs whose first measurement lies behind the fitted vertex (`leavesAhead`)
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
    visit("gapXY", gapXY);
    visit("gapZ", gapZ);
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
    visit("ahead", ahead);
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
  float length = std::numeric_limits<float>::quiet_NaN();
  float sLength = std::numeric_limits<float>::quiet_NaN();
  float cosPointing = std::numeric_limits<float>::quiet_NaN();
  float dXY = 0.F; ///< the fitted vertex's transverse distance from the MC decay point [cm]
  float dZ = 0.F;
  float opening = 0.F; ///< the daughters' MC opening angle at production [rad]
  /// Its own legs' tracks (both found): their rows; with their refits under
  /// the true species, those refits' innermost states against the MC
  /// crossing there (`stateTruth`), pulls of loc0, loc1, phi, theta, q/p;
  /// NaN without the crossing.
  std::uint8_t rowsPositive = 0U;  ///< its own leg's track: rows
  std::uint8_t wrongPositive = 0U; ///< its rows not of the daughter (`RowParticles::carries`)
  std::uint8_t layerPositive = 0U; ///< its innermost row's layer
  std::uint8_t rowsNegative = 0U;  ///< its own leg's track: rows
  std::uint8_t wrongNegative = 0U; ///< its rows not of the daughter (`RowParticles::carries`)
  std::uint8_t layerNegative = 0U; ///< its innermost row's layer
  static constexpr float nan = std::numeric_limits<float>::quiet_NaN();
  float pullLoc0Positive = nan;
  float pullLoc1Positive = nan;
  float pullPhiPositive = nan;
  float pullThetaPositive = nan;
  float pullQOverPPositive = nan; ///< charge-signed: positive where |q/p| is too large (the momentum too small)
  float dpRelPositive = nan;      ///< fitted momentum over the MC's, less 1
  float pullLoc0Negative = nan;
  float pullLoc1Negative = nan;
  float pullPhiNegative = nan;
  float pullThetaNegative = nan;
  float pullQOverPNegative = nan; ///< charge-signed: positive where |q/p| is too large (the momentum too small)
  float dpRelNegative = nan;      ///< fitted momentum over the MC's, less 1

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
    visit("opening", opening);
    visit("rowsPositive", rowsPositive);
    visit("wrongPositive", wrongPositive);
    visit("layerPositive", layerPositive);
    visit("pullLoc0Positive", pullLoc0Positive);
    visit("pullLoc1Positive", pullLoc1Positive);
    visit("pullPhiPositive", pullPhiPositive);
    visit("pullThetaPositive", pullThetaPositive);
    visit("pullQOverPPositive", pullQOverPPositive);
    visit("dpRelPositive", dpRelPositive);
    visit("rowsNegative", rowsNegative);
    visit("wrongNegative", wrongNegative);
    visit("layerNegative", layerNegative);
    visit("pullLoc0Negative", pullLoc0Negative);
    visit("pullLoc1Negative", pullLoc1Negative);
    visit("pullPhiNegative", pullPhiNegative);
    visit("pullThetaNegative", pullThetaNegative);
    visit("pullQOverPNegative", pullQOverPNegative);
    visit("dpRelNegative", dpRelNegative);
  }
};

namespace decay_finder_detail
{

/// The collisions fitted and written as one block (the file's comment):
/// what bounds the records held at once.
inline constexpr std::size_t kCollisionsPerBlock = 256U;

/// Where a leg's columns of a `DecayFinderDecayRecord` go.
struct LegTruth {
  std::uint8_t* rows = nullptr;
  std::uint8_t* wrong = nullptr;
  std::uint8_t* layer = nullptr;
  std::array<float*, 5> pulls{};
  float* dpRel = nullptr;
};

/// \a path with \a tag before its extension.
[[nodiscard]] inline std::string taggedPath(const std::string& path, std::string_view tag)
{
  const auto dot = path.rfind('.');
  return path.substr(0, dot) + std::string{tag} + path.substr(dot);
}

} // namespace decay_finder_detail

template <DetectorLayout Det>
class DecayFinderStudy
{
 public:
  explicit DecayFinderStudy(RecordWriter& writer)
    : mWriter(&writer),
      mDecayWriter(std::make_unique<RecordWriter>(decay_finder_detail::taggedPath(writer.path(), "-decays"))),
      mCascadeWriter(std::make_unique<RecordWriter>(decay_finder_detail::taggedPath(writer.path(), "-cascades"))),
      mCascadeDecayWriter(std::make_unique<RecordWriter>(decay_finder_detail::taggedPath(writer.path(), "-cascadeDecays")))
  {
    writer.begin<DecayFinderCandidateRecord>();
    mDecayWriter->begin<DecayFinderDecayRecord>();
    mCascadeWriter->begin<CascadeFinderCandidateRecord>();
    mCascadeDecayWriter->begin<CascadeFinderDecayRecord>();
  }

  /// \a event's tracks, every one, whatever the source's chains.
  void process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source);

  [[nodiscard]] std::vector<std::string> summary() const
  {
    return {std::format("decay finder study: {} tracks combined ({} primary vertices' contributors and {} without a particle not), {} candidates, {} fitted (records into {}), {} MC V0s ({} found) into {}; the daughter table {} entries in {:.3f} s, the S^2 in {:.3f} s, the vertices in {:.3f} s, the fits {:.1f} thread-s",
                        mTracks, mContributors, mUnassociated, mPairs, mFitted, mWriter->path(), mDecays, mFound, mDecayWriter->path(), mEntries, mFlushSeconds, mDisplacementSeconds, mVertexSeconds, mFitMicroseconds * 1e-6),
            std::format("cascades: {} Lambdas taken, {} candidates, {} fitted (records into {}, their fits {:.1f} thread-s), {} MC cascades ({} found) into {}", mLambdasTaken, mCascadeCandidates, mCascadesFitted, mCascadeWriter->path(), mCascadeMicroseconds * 1e-6, mCascades, mCascadesFound, mCascadeDecayWriter->path())};
  }

 private:
  RecordWriter* mWriter = nullptr;
  std::unique_ptr<RecordWriter> mDecayWriter;
  std::unique_ptr<RecordWriter> mCascadeWriter;
  std::unique_ptr<RecordWriter> mCascadeDecayWriter;
  std::size_t mLambdasTaken = 0U;
  std::size_t mCascadeCandidates = 0U;
  std::size_t mCascadesFitted = 0U;
  std::size_t mCascades = 0U;
  std::size_t mCascadesFound = 0U;
  double mCascadeMicroseconds = 0.;
  std::size_t mTracks = 0U;
  std::size_t mUnassociated = 0U;
  std::size_t mContributors = 0U;
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
  const auto& trees = *event.trees;
  const auto committed = store.tracks().view();

  const std::array<PID, 3> species{PID{PID::Pion}, PID{PID::Proton}, PID{PID::Kaon}};
  const auto hypotheses = finderHypotheses();
  const std::span<const DecayHypothesis> allHypotheses{hypotheses};
  const auto& param = DecayParam::Instance();
  const double lambdaMassGeV = hypotheses[1].massGeV;

  // The tracks associated to the mocked primary vertices; a vertex's
  // contributors are not combined.
  const auto vertexStart = std::chrono::steady_clock::now();
  const auto association = mockAssociation(event);
  mVertexSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - vertexStart).count();
  const auto& points = association.vertices.points;
  const auto& particleOf = association.particle;
  const auto& fake = association.fake;
  const auto n = static_cast<std::uint32_t>(store.size());
  std::unordered_map<std::uint64_t, std::array<std::vector<std::uint32_t>, 2>> byCollision; // positive, negative
  std::vector<std::uint64_t> collisions;
  for (std::uint32_t t = 0; t < n; ++t) {
    if (particleOf[t] == kNoTrackParticle) {
      ++mUnassociated;
      continue;
    }
    if (association.contributor(t)) {
      ++mContributors;
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

  // The table: every combined track as a pion, a proton and a kaon.
  DaughterTable<Det> table(std::span<const PID>{species});
  const auto flushStart = std::chrono::steady_clock::now();
  fillDaughterTable(table, event, "decay finder study", [&](DaughterTable<Det>& asked) {
    for (const auto& [collision, lists] : byCollision) {
      for (const auto& list : lists) {
        for (const auto t : list) {
          for (std::size_t s = 0; s < species.size(); ++s) {
            asked.request(t, static_cast<std::uint8_t>(s));
          }
        }
      }
    }
  });
  mFlushSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - flushStart).count();
  mEntries += table.size();
  const auto daughters = daughterSource(table, event);

  // Each entry's displacement against its collision's point.
  const auto entries = table.view();
  const auto displacementStart = std::chrono::steady_clock::now();
  std::vector<float> s2(table.size(), std::numeric_limits<float>::quiet_NaN());
  std::vector<std::uint8_t> noPoint(table.size(), 0U);
  eventForEach(event)(static_cast<std::uint32_t>(table.size()), [&](std::uint32_t e) {
    if (entries.status[e] != PropagationStatus::kOk) {
      return;
    }
    const auto point = association.vertices.of(collisionOf(particleOf[entries.track[e]]));
    if (point < 0) {
      noPoint[e] = 1U;
      return;
    }
    s2[e] = static_cast<float>(displacement(daughters, e, points[static_cast<std::size_t>(point)]));
  });
  mDisplacementSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - displacementStart).count();

  // The MC V0s decaying into their two charged daughters, by daughter.
  struct Decay {
    const TreeParticle* mother = nullptr;
    std::array<const TreeParticle*, 2> legs{}; ///< the positive, the negative daughter
    int hypothesis = 0;
    DecayFinderDecayRecord record;
    bool fitted = false;
    LambdaRejection rejection = LambdaRejection::kNotFitted; ///< a Lambda's own candidate's, for its cascade
  };
  std::vector<Decay> decays;
  std::unordered_map<std::uint64_t, std::size_t> decayOfLeg;
  for (const auto& member : trees.particles) {
    const auto found = std::ranges::find(kV0Mothers, member.pdg);
    if (found == kV0Mothers.end()) {
      continue;
    }
    const auto h = static_cast<int>(found - kV0Mothers.begin());
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
  // The MC cascades into a Lambda (an MC V0 above) and a charged pion
  // (Xi) or kaon (Omega), by Lambda and by bachelor.
  std::unordered_map<std::uint64_t, std::size_t> decayOfMother;
  for (std::size_t i = 0; i < decays.size(); ++i) {
    decayOfMother.emplace(decays[i].mother->id, i);
  }
  struct Cascade {
    const TreeParticle* mother = nullptr;
    std::size_t lambda = 0U; ///< into `decays`
    const TreeParticle* bachelor = nullptr;
    int sign = 0; ///< 0 a Lambda's, 1 an anti-Lambda's
    int omega = 0;
    CascadeFinderDecayRecord record;
    bool fitted = false;
  };
  std::vector<Cascade> cascades;
  std::unordered_map<std::uint64_t, std::size_t> cascadeOfLambda;
  for (const auto& member : trees.particles) {
    const int sign = member.pdg == 3312 || member.pdg == 3334 ? 0 : (member.pdg == -3312 || member.pdg == -3334 ? 1 : -1);
    if (sign < 0) {
      continue;
    }
    const int omega = std::abs(member.pdg) == 3334 ? 1 : 0;
    const auto products = decayProducts(trees, member);
    if (products.size() != 2U) {
      continue;
    }
    const bool lambdaFirst = products[0]->pdg == (sign == 0 ? 3122 : -3122);
    const auto* lambda = lambdaFirst ? products[0] : products[1];
    const auto* bachelor = lambdaFirst ? products[1] : products[0];
    const auto found = decayOfMother.find(lambda->id);
    if (lambda->pdg != (sign == 0 ? 3122 : -3122) || bachelor->pdg != (omega != 0 ? -321 : -211) * (sign == 0 ? 1 : -1) || found == decayOfMother.end()) {
      continue;
    }
    const auto& p = member.momentum;
    const auto& l = decays[found->second].legs;
    cascadeOfLambda.emplace(lambda->id, cascades.size());
    cascades.push_back({.mother = &member, .lambda = found->second, .bachelor = bachelor, .sign = sign, .omega = omega, .record = {.timeframe = event.timeframe, .particle = member.id, .pdg = member.pdg, .decayR = static_cast<float>(std::hypot(lambda->vertex[0], lambda->vertex[1])), .lambdaDecayR = decays[found->second].record.decayR, .pt = static_cast<float>(std::hypot(p[0], p[1])), .ptPositive = static_cast<float>(std::hypot(l[0]->momentum[0], l[0]->momentum[1])), .ptNegative = static_cast<float>(std::hypot(l[1]->momentum[0], l[1]->momentum[1])), .ptBachelor = static_cast<float>(std::hypot(bachelor->momentum[0], bachelor->momentum[1]))}});
  }
  std::unordered_map<std::uint64_t, std::size_t> cascadeOfBachelor;
  for (std::size_t i = 0; i < cascades.size(); ++i) {
    cascadeOfBachelor.emplace(cascades[i].bachelor->id, i);
  }
  // Each MC daughter's longest good track of its own charge (one of the
  // other charge runs against the particle: a looper's way back); the
  // daughters with good tracks of the other charge alone.
  std::unordered_map<std::uint64_t, std::uint32_t> trackOfLeg;
  std::unordered_set<std::uint64_t> reversedLeg;
  for (std::uint32_t t = 0; t < n; ++t) {
    if (fake[t] != 0U || (!decayOfLeg.contains(particleOf[t]) && !cascadeOfBachelor.contains(particleOf[t]))) {
      continue;
    }
    if ((committed.state[t].qOverP() > 0.F) != (particleCharge(event.particles->at(particleOf[t]).pdg) > 0.)) {
      reversedLeg.insert(particleOf[t]);
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
  const auto radiusOf = [&](std::int32_t e) {
    const auto& state = entries.state[static_cast<std::size_t>(e)];
    std::array<Scalar, 3> at{};
    surfaceToGlobal(geometry.planeOf(state), state.loc0(), state.loc1(), at.data());
    return std::hypot(at[0], at[1]);
  };
  std::vector<std::vector<DecayFinderCandidateRecord>> perCollision(std::min(collisions.size(), kCollisionsPerBlock));
  std::vector<std::vector<CascadeFinderCandidateRecord>> perCollisionCascades(perCollision.size());
  // A Lambda candidate a cascade takes (`DecayParam`), from the V0 pass.
  struct TakenLambda {
    std::uint32_t positive = 0U;
    std::uint32_t negative = 0U;
    int sign = 0; ///< 0 Lambda, 1 anti-Lambda
    std::int32_t ePositive = -1;
    std::int32_t eNegative = -1;
    const DecayFinderCandidateRecord* record = nullptr; ///< its V0 record (class, mother, PID, S^2, mass)
    bool own = false;                                   ///< its MC Lambda's own pair
  };
  const auto rejectionOf = [&](const DecayFinderCandidateRecord& record, int ahead) {
    if (record.status != static_cast<std::uint8_t>(VertexFitStatus::kOk)) {
      return LambdaRejection::kNotFitted;
    }
    if (!(std::abs(static_cast<double>(record.mass) - lambdaMassGeV) < param.lambdaNSigma * static_cast<double>(record.sMass))) {
      return LambdaRejection::kMass;
    }
    if (ahead != 0) {
      return LambdaRejection::kAhead;
    }
    const auto low = [&](float s2Leg) { return param.lambdaMinS2 > 0. && !(static_cast<double>(s2Leg) >= param.lambdaMinS2); };
    return low(record.s2Positive) || low(record.s2Negative) ? LambdaRejection::kS2 : LambdaRejection::kTaken;
  };
  std::vector<std::size_t> perCollisionTaken(perCollision.size());
  // The cascades' own rows (MotherSearch.h): the rows any committed track holds, flagged.
  const MotherSearch<Det> motherSearch(*event.frame, geometry, param.motherNSigma);
  std::vector<std::uint8_t> owned(event.frame->size(), 0U);
  for (std::uint32_t t = 0; t < n; ++t) {
    for (const auto row : store.rows(t)) {
      owned[row] = 1U;
    }
  }
  const double fieldScale = vertex_fit_detail::fieldScaleOf(event.propagator);
  // A collision's cascade candidates: every taken Lambda with every bachelor
  // of its charge, under Xi and Omega, the Lambda's mass free and held.
  const auto findCascades = [&](std::uint64_t collision, std::int32_t point, std::span<const TakenLambda> taken, WorkspaceBuffer& workspace, std::vector<CascadeFinderCandidateRecord>& out) {
    out.clear();
    const auto& lists = byCollision.at(collision);
    for (const auto& lambda : taken) {
      const auto& v0 = *lambda.record;
      // its MC cascade, if the pair is an MC Lambda (its own hypothesis) from one
      const auto lambdaCascade = v0.truthClass == static_cast<std::uint8_t>(FinderClass::kTrue) ? cascadeOfLambda.find(v0.mother) : cascadeOfLambda.end();
      for (const auto bachelor : lists[lambda.sign == 0 ? 1 : 0]) {
        if (bachelor == lambda.positive || bachelor == lambda.negative) {
          continue;
        }
        const bool clone = fake[bachelor] == 0U && ((fake[lambda.positive] == 0U && particleOf[bachelor] == particleOf[lambda.positive]) || (fake[lambda.negative] == 0U && particleOf[bachelor] == particleOf[lambda.negative]));
        const bool trueBachelor = lambdaCascade != cascadeOfLambda.end() && fake[bachelor] == 0U && particleOf[bachelor] == cascades[lambdaCascade->second].bachelor->id;
        const auto bachelorTruth = truthSpecies(event.particles->at(particleOf[bachelor]).pdg);
        for (int omega = 0; omega < 2; ++omega) {
          const PID bachelorPid = omega != 0 ? PID{PID::Kaon} : PID{PID::Pion};
          const auto eBachelor = table.entryOf(bachelor, static_cast<std::uint8_t>(table.speciesIndex(bachelorPid)));
          const bool bachelorRefitted = entries.status[static_cast<std::size_t>(eBachelor)] == PropagationStatus::kOk;
          float s2Bachelor = s2[static_cast<std::size_t>(eBachelor)];
          for (int held = 0; held < 2; ++held) {
            const auto h = cascadeIndex(lambda.sign, omega, held);
            CascadeFinderCandidateRecord record{.timeframe = event.timeframe, .collision = collision, .sign = static_cast<std::uint8_t>(lambda.sign), .omega = static_cast<std::uint8_t>(omega), .held = static_cast<std::uint8_t>(held), .lambdaClass = v0.truthClass, .s2Positive = v0.s2Positive, .s2Negative = v0.s2Negative, .s2Bachelor = s2Bachelor, .lambdaV0Mass = v0.mass, .lambdaV0SMass = v0.sMass};
            record.pidMatch = static_cast<std::uint8_t>(v0.pidMatch != 0U && bachelorTruth && bachelorTruth->pid == bachelorPid ? 1U : 0U);
            record.primaries = static_cast<std::uint8_t>(v0.primaries + (isPrimary(bachelor) ? 1 : 0));
            const bool lambdaFake = v0.truthClass == static_cast<std::uint8_t>(FinderClass::kFake);
            if (lambdaFake || fake[bachelor] != 0U) {
              record.truthClass = static_cast<std::uint8_t>(CascadeClass::kFake);
            } else if (clone) {
              record.truthClass = static_cast<std::uint8_t>(CascadeClass::kCloneBachelor);
            } else if (v0.truthClass == static_cast<std::uint8_t>(FinderClass::kTrue)) {
              if (trueBachelor) {
                const auto& cascade = cascades[lambdaCascade->second];
                record.truthClass = static_cast<std::uint8_t>(cascade.omega == omega ? CascadeClass::kTrue : CascadeClass::kReflection);
                record.cascade = cascade.mother->id;
              } else {
                record.truthClass = static_cast<std::uint8_t>(lambdaCascade != cascadeOfLambda.end() ? CascadeClass::kOtherBachelor : CascadeClass::kPrimaryLambda);
              }
            } else if (v0.truthClass == static_cast<std::uint8_t>(FinderClass::kOtherHypothesis) || v0.truthClass == static_cast<std::uint8_t>(FinderClass::kWrongCharge)) {
              record.truthClass = static_cast<std::uint8_t>(CascadeClass::kWrongLambda);
            } else if (v0.truthClass == static_cast<std::uint8_t>(FinderClass::kOtherDecay)) {
              record.truthClass = static_cast<std::uint8_t>(CascadeClass::kOtherMother);
            } else {
              record.truthClass = static_cast<std::uint8_t>(CascadeClass::kCombinatorial);
            }
            DecayCandidate candidate{.hypothesis = static_cast<std::uint16_t>(h), .point = point};
            candidate.entries.push_back(static_cast<std::uint32_t>(lambda.ePositive));
            candidate.entries.push_back(static_cast<std::uint32_t>(lambda.eNegative));
            candidate.entries.push_back(static_cast<std::uint32_t>(eBachelor));
            DecayCandidateResult result;
            const auto started = std::chrono::steady_clock::now();
            if (bachelorRefitted) {
              static_cast<void>(startDecayCandidate(candidate, allHypotheses, daughters, points.size()));
            }
            fitDecayCandidate(candidate, allHypotheses, daughters, geometry.material(), std::span<const MeasuredPoint>{points}, workspace, result);
            record.microseconds = static_cast<float>(std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - started).count());
            writeCandidateFit(record, result);
            const bool fitted = result.fit.status == VertexFitStatus::kOk;
            if (fitted) {
              const auto& fit = result.fit;
              record.ahead = static_cast<std::uint8_t>(leavesAhead(candidate, allHypotheses, daughters, result));
              record.lambdaMass = static_cast<float>(result.composites[0].mass);
              record.sLambdaMass = static_cast<float>(result.composites[0].sMass);
              record.radius = static_cast<float>(std::hypot(fit.vertices[0][0], fit.vertices[0][1]));
              record.lambdaRadius = static_cast<float>(std::hypot(fit.vertices[1][0], fit.vertices[1][1]));
              record.lambdaLength = static_cast<float>(result.compositeFlights[0].length);
              record.sLambdaLength = static_cast<float>(result.compositeFlights[0].sLength);
              record.lambdaCosPointing = static_cast<float>(result.compositeFlights[0].cosPointing);
            }
            if (fitted && held != 0) {
              // The mother's own rows, and the fit with them.
              const auto tree = expandHypothesis(allHypotheses, h);
              const auto start = motherStateAt(result.fit, tree.topology, 0, fieldScale);
              const std::array<std::uint32_t, 3> tracks{lambda.positive, lambda.negative, bachelor};
              TimeEstPS box{.start = committed.timePS[tracks[0]], .width = committed.widthPS[tracks[0]]};
              bool together = true;
              for (std::size_t k = 1; k < tracks.size() && together; ++k) {
                const TimeEstPS other{.start = committed.timePS[tracks[k]], .width = committed.widthPS[tracks[k]]};
                together = box.isCompatible(other);
                if (together) {
                  box += other;
                }
              }
              if (start.found && together) {
                auto propagator = event.propagator;
                propagator.massGeV = static_cast<Scalar>(hypotheses[h].massGeV);
                const auto own = [&](std::uint32_t row) { return std::ranges::any_of(tracks, [&](std::uint32_t t) { return std::ranges::find(store.rows(t), row) != store.rows(t).end(); }); };
                const auto found = motherSearch.search(start, propagator, box, own);
                record.motherLayers = found.layers;
                record.motherRows = static_cast<std::uint8_t>(found.rows.size());
                if (!found.rows.empty()) {
                  float chi2 = 0.F;
                  std::vector<TrackingSurface> planes;
                  std::vector<Vec2D> positions;
                  std::vector<Vec3D> covariance;
                  for (std::size_t k = 0; k < found.rows.size(); ++k) {
                    const auto row = found.rows[k];
                    chi2 += found.chi2[k];
                    record.motherFree += owned[row] == 0U ? 1U : 0U;
                    record.motherTrue += record.cascade != 0U && event.rows[row].carries(record.cascade) ? 1U : 0U;
                    const auto m = event.frame->measurement(row);
                    planes.push_back(geometry.planeOf(m));
                    positions.push_back({m.loc0, m.loc1});
                    covariance.push_back({m.varLoc0, m.covLoc01, m.varLoc1});
                  }
                  record.motherChi2 = chi2;
                  const std::array<DecayMother, 1> mothers{DecayMother{.vertex = 0, .planes = planes, .positions = positions, .covariance = covariance, .propagator = propagator}};
                  DecayCandidateResult tracked;
                  fitDecayCandidate(candidate, allHypotheses, daughters, geometry.material(), std::span<const MeasuredPoint>{points}, workspace, tracked, std::span<const DecayMother>{mothers});
                  record.trackedStatus = static_cast<std::uint8_t>(tracked.fit.status);
                  if (tracked.fit.status == VertexFitStatus::kOk) {
                    record.trackedChi2 = static_cast<float>(tracked.fit.chi2);
                    record.trackedNdf = static_cast<std::int16_t>(tracked.fit.ndf);
                    record.trackedMass = static_cast<float>(tracked.kinematics.mass);
                    record.trackedSMass = static_cast<float>(tracked.kinematics.sMass);
                    record.trackedRadius = static_cast<float>(std::hypot(tracked.fit.vertices[0][0], tracked.fit.vertices[0][1]));
                    if (tracked.hasPoint) {
                      record.trackedLength = static_cast<float>(tracked.flight.length);
                      record.trackedSLength = static_cast<float>(tracked.flight.sLength);
                      record.trackedCosPointing = static_cast<float>(tracked.flight.cosPointing);
                    }
                    if (record.cascade != 0U) {
                      const auto& decayPoint = decayPointOf(trees, *cascades[lambdaCascade->second].mother);
                      record.trackedDXY = static_cast<float>(std::hypot(tracked.fit.vertex[0] - decayPoint[0], tracked.fit.vertex[1] - decayPoint[1]));
                      record.trackedDZ = static_cast<float>(tracked.fit.vertex[2] - decayPoint[2]);
                    }
                  }
                }
              }
            }
            if (record.cascade != 0U) {
              auto& cascade = cascades[lambdaCascade->second];
              const auto& decayPoint = decayPointOf(trees, *cascade.mother);
              if (fitted) {
                record.dXY = static_cast<float>(std::hypot(result.fit.vertex[0] - decayPoint[0], result.fit.vertex[1] - decayPoint[1]));
                record.dZ = static_cast<float>(result.fit.vertex[2] - decayPoint[2]);
              }
              auto& d = cascade.record;
              if (record.truthClass == static_cast<std::uint8_t>(CascadeClass::kTrue) && held != 0 && fitted) {
                d.recovered = 1U;
              }
              // its own candidate: its Lambda's own pair, its bachelor's own track, its own hypothesis (one writer per cascade)
              const auto own = trackOfLeg.find(cascade.bachelor->id);
              if (record.truthClass == static_cast<std::uint8_t>(CascadeClass::kTrue) && lambda.own && own != trackOfLeg.end() && own->second == bachelor) {
                cascade.fitted = true;
                d.s2Bachelor = s2Bachelor;
                if (held != 0) {
                  writeCandidateFit(d, result);
                  d.dXY = record.dXY;
                  d.dZ = record.dZ;
                  d.motherRows = record.motherRows;
                  d.motherTrue = record.motherTrue;
                  d.trackedStatus = record.trackedStatus;
                  d.trackedMass = record.trackedMass;
                  d.trackedSMass = record.trackedSMass;
                  d.trackedDXY = record.trackedDXY;
                  d.trackedDZ = record.trackedDZ;
                } else {
                  d.statusFree = record.status;
                  d.massFree = record.mass;
                  d.sMassFree = record.sMass;
                  d.lambdaMass = record.lambdaMass;
                  d.sLambdaMass = record.sLambdaMass;
                }
              }
            }
            out.push_back(record);
          }
        }
      }
    }
  };
  for (std::size_t first = 0; first < collisions.size(); first += kCollisionsPerBlock) {
    const std::size_t last = std::min(collisions.size(), first + kCollisionsPerBlock);
    o2::gtrk::detail::parallelFor(*event.arena, first, last, 1, [&](std::size_t begin, std::size_t end) {
      WorkspaceBuffer workspace;
      for (std::size_t c = begin; c < end; ++c) {
        const auto collision = collisions[c];
        const auto& lists = byCollision.at(collision);
        const auto point = association.vertices.of(collision);
        auto& out = perCollision[c - first];
        out.clear();
        // reserved: the taken Lambdas point into `out`
        out.reserve(lists[0].size() * lists[1].size() * kV0Mothers.size());
        std::vector<TakenLambda> taken;
        for (const auto positive : lists[0]) {
          for (const auto negative : lists[1]) {
            const auto mother = motherOf(positive);
            const bool sameMother = mother != kNoTreeParticle && mother == motherOf(negative);
            const auto ownDecay = sameMother ? decayOfLeg.find(particleOf[positive]) : decayOfLeg.end();
            const bool ownPair = ownDecay != decayOfLeg.end() && trackOfLeg.contains(particleOf[positive]) && trackOfLeg.at(particleOf[positive]) == positive && trackOfLeg.contains(particleOf[negative]) && trackOfLeg.at(particleOf[negative]) == negative;
            for (std::size_t h = 0; h < kV0Mothers.size(); ++h) {
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
              DecayCandidate candidate{.hypothesis = static_cast<std::uint16_t>(h), .point = point};
              candidate.entries.push_back(static_cast<std::uint32_t>(ePositive));
              candidate.entries.push_back(static_cast<std::uint32_t>(eNegative));
              DecayCandidateResult result;
              const auto started = std::chrono::steady_clock::now();
              if (refitted(ePositive) && refitted(eNegative)) {
                static_cast<void>(startDecayCandidate(candidate, std::span<const DecayHypothesis>{hypotheses}, daughters, points.size()));
                // The legs' helices where they come closest: the start of a pair.
                const auto legs = candidateDaughters(candidate, daughters);
                const auto entryOf = [](const DecayDaughter& daughter) { return decay_fit_detail::StartEntry{.state = &daughter.innermost, .surface = &daughter.surface, .propagator = &daughter.propagator}; };
                if (const auto approach = decay_fit_detail::approachOf(entryOf(legs[0]), entryOf(legs[1])); approach.found) {
                  record.startRadius = static_cast<float>(std::hypot(approach.start[0], approach.start[1]));
                  record.gapXY = static_cast<float>(approach.gapXY);
                  record.gapZ = static_cast<float>(approach.gapZ);
                }
              }
              fitDecayCandidate(candidate, std::span<const DecayHypothesis>{hypotheses}, daughters, geometry.material(), std::span<const MeasuredPoint>{points}, workspace, result);
              record.microseconds = static_cast<float>(std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - started).count());
              writeCandidateFit(record, result);
              if (result.fit.status == VertexFitStatus::kOk) {
                record.radius = static_cast<float>(std::hypot(result.fit.vertex[0], result.fit.vertex[1]));
                record.ahead = static_cast<std::uint8_t>(leavesAhead(candidate, allHypotheses, daughters, result));
              }
              const auto rejection = h == 0 ? LambdaRejection::kNotFitted : rejectionOf(record, record.ahead);
              if (own && h != 0) {
                decays[ownDecay->second].rejection = rejection;
              }
              if (own) {
                // The decay's own candidate: its record's fit (one writer per decay).
                auto& decay = decays[ownDecay->second];
                auto& d = decay.record;
                decay.fitted = true;
                writeCandidateFit(d, result);
                d.s2Positive = s2Positive;
                d.s2Negative = s2Negative;
                if (result.fit.status == VertexFitStatus::kOk) {
                  const auto& decayPoint = decayPointOf(trees, *decay.mother);
                  d.dXY = static_cast<float>(std::hypot(result.fit.vertex[0] - decayPoint[0], result.fit.vertex[1] - decayPoint[1]));
                  d.dZ = static_cast<float>(result.fit.vertex[2] - decayPoint[2]);
                }
              }
              out.push_back(record);
              if (rejection == LambdaRejection::kTaken) {
                taken.push_back({.positive = positive, .negative = negative, .sign = static_cast<int>(h) - 1, .ePositive = ePositive, .eNegative = eNegative, .record = &out.back(), .own = own});
              }
            }
          }
        }
        perCollisionTaken[c - first] = taken.size();
        findCascades(collision, point, taken, workspace, perCollisionCascades[c - first]);
      }
    });
    for (std::size_t c = first; c < last; ++c) {
      for (const auto& record : perCollision[c - first]) {
        mWriter->write(record);
        mFitted += record.status == static_cast<std::uint8_t>(VertexFitStatus::kOk) ? 1U : 0U;
        mFitMicroseconds += record.microseconds;
      }
      mPairs += perCollision[c - first].size();
      for (const auto& record : perCollisionCascades[c - first]) {
        mCascadeWriter->write(record);
        mCascadesFitted += record.status == static_cast<std::uint8_t>(VertexFitStatus::kOk) ? 1U : 0U;
        mCascadeMicroseconds += record.microseconds;
      }
      mCascadeCandidates += perCollisionCascades[c - first].size();
      mLambdasTaken += perCollisionTaken[c - first];
    }
  }

  // Each MC V0's stage.
  for (auto& decay : decays) {
    auto& d = decay.record;
    const auto positive = trackOfLeg.find(decay.legs[0]->id);
    const auto negative = trackOfLeg.find(decay.legs[1]->id);
    const auto stageOf = [&]() {
      if (positive == trackOfLeg.end() || negative == trackOfLeg.end()) {
        const bool reversed = (positive == trackOfLeg.end() && reversedLeg.contains(decay.legs[0]->id)) || (negative == trackOfLeg.end() && reversedLeg.contains(decay.legs[1]->id));
        return reversed ? FinderStage::kWrongCharge : FinderStage::kNotFound;
      }
      if (association.contributor(positive->second) || association.contributor(negative->second)) {
        return FinderStage::kContributor;
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
        throw std::logic_error("decay finder study: an MC V0's own pair of its own charges was never fitted");
      }
      if (d.status == static_cast<std::uint8_t>(VertexFitStatus::kNoStart)) {
        return FinderStage::kNoStart;
      }
      return d.status == static_cast<std::uint8_t>(VertexFitStatus::kOk) ? FinderStage::kFound : FinderStage::kFitFailed;
    };
    d.stage = static_cast<std::uint8_t>(stageOf());
    {
      const auto& a = decay.legs[0]->momentum;
      const auto& b = decay.legs[1]->momentum;
      d.opening = static_cast<float>(std::acos(std::clamp(dotOf(a, b) / std::sqrt(dotOf(a, a) * dotOf(b, b)), -1., 1.)));
    }
    if (positive != trackOfLeg.end() && negative != trackOfLeg.end()) {
      const auto& hypothesis = hypotheses[static_cast<std::size_t>(decay.hypothesis)];
      for (std::size_t slot = 0; slot < 2U; ++slot) {
        const auto t = (slot == 0 ? positive : negative)->second;
        const auto particle = decay.legs[slot]->id;
        const auto rows = store.rows(t);
        int wrong = 0;
        for (const auto row : rows) {
          wrong += event.rows[row].carries(particle) ? 0 : 1;
        }
        const LegTruth leg = {.rows = &(slot == 0 ? d.rowsPositive : d.rowsNegative), .wrong = &(slot == 0 ? d.wrongPositive : d.wrongNegative), .layer = &(slot == 0 ? d.layerPositive : d.layerNegative), .pulls = {slot == 0 ? &d.pullLoc0Positive : &d.pullLoc0Negative, slot == 0 ? &d.pullLoc1Positive : &d.pullLoc1Negative, slot == 0 ? &d.pullPhiPositive : &d.pullPhiNegative, slot == 0 ? &d.pullThetaPositive : &d.pullThetaNegative, slot == 0 ? &d.pullQOverPPositive : &d.pullQOverPNegative}, .dpRel = slot == 0 ? &d.dpRelPositive : &d.dpRelNegative};
        *leg.rows = static_cast<std::uint8_t>(std::min<std::size_t>(rows.size(), 255U));
        *leg.wrong = static_cast<std::uint8_t>(std::min(wrong, 255));
        *leg.layer = static_cast<std::uint8_t>(event.rows[rows.front()].layer);
        const auto e = table.entryOf(t, static_cast<std::uint8_t>(table.speciesIndex(hypothesis.species[slot])));
        if (e < 0 || entries.status[static_cast<std::size_t>(e)] != PropagationStatus::kOk || event.hits == nullptr) {
          continue;
        }
        const double charge = particleCharge(decay.legs[slot]->pdg);
        const auto& state = entries.state[static_cast<std::size_t>(e)];
        const auto truth = stateTruth(geometry, *event.hits, *event.frame, state, particle, rows.front(), charge);
        if (!truth) {
          continue;
        }
        for (int k = 0; k < 5; ++k) {
          const double sign = k == 4 ? (charge > 0. ? 1. : -1.) : 1.;
          *leg.pulls[static_cast<std::size_t>(k)] = static_cast<float>(sign * truth->residual[static_cast<std::size_t>(k)] / truth->sigma[static_cast<std::size_t>(k)]);
        }
        *leg.dpRel = static_cast<float>((1. / std::abs(static_cast<double>(state.qOverP()))) / truth->momentum - 1.);
      }
    }
    mFound += d.stage == static_cast<std::uint8_t>(FinderStage::kFound) ? 1U : 0U;
    mDecayWriter->write(d);
  }
  mDecays += decays.size();

  // Each MC cascade's stage.
  for (auto& cascade : cascades) {
    auto& d = cascade.record;
    const auto& lambda = decays[cascade.lambda];
    const std::array<const TreeParticle*, 3> leaves{lambda.legs[0], lambda.legs[1], cascade.bachelor};
    std::array<std::uint32_t, 3> tracks{};
    for (std::size_t i = 0; i < leaves.size(); ++i) {
      const auto found = trackOfLeg.find(leaves[i]->id);
      if (found == trackOfLeg.end()) {
        d.missing |= static_cast<std::uint8_t>(1U << i);
        d.reversed |= static_cast<std::uint8_t>(reversedLeg.contains(leaves[i]->id) ? 1U << i : 0U);
      } else {
        tracks[i] = found->second;
      }
    }
    d.hasPoint = static_cast<std::uint8_t>(association.vertices.of(collisionOf(cascade.mother->id)) >= 0 ? 1U : 0U);
    if (const auto* own = event.chains != nullptr ? event.chains->find(cascade.mother->id) : nullptr; own != nullptr) {
      d.mcRows = static_cast<std::uint8_t>(std::min<std::size_t>(own->rows.size(), 255U));
    }
    d.lambdaRejection = static_cast<std::uint8_t>(lambda.rejection);
    const auto stageOf = [&]() {
      if (d.missing != 0U) {
        return CascadeStage::kNotFound;
      }
      if (std::ranges::any_of(tracks, [&](std::uint32_t t) { return association.contributor(t); })) {
        return CascadeStage::kContributor;
      }
      const auto& lambdaHypothesis = hypotheses[static_cast<std::size_t>(lambda.hypothesis)];
      for (std::size_t i = 0; i < 2U; ++i) {
        const auto e = table.entryOf(tracks[i], static_cast<std::uint8_t>(table.speciesIndex(lambdaHypothesis.species[i])));
        if (e < 0 || entries.status[static_cast<std::size_t>(e)] != PropagationStatus::kOk) {
          return CascadeStage::kLambdaRefitFailed;
        }
      }
      if (lambda.record.status != static_cast<std::uint8_t>(VertexFitStatus::kOk)) {
        return CascadeStage::kLambdaFitFailed;
      }
      if (lambda.rejection != LambdaRejection::kTaken) {
        return CascadeStage::kLambdaRejected;
      }
      const auto eBachelor = table.entryOf(tracks[2], static_cast<std::uint8_t>(table.speciesIndex(cascade.omega != 0 ? PID{PID::Kaon} : PID{PID::Pion})));
      if (eBachelor < 0 || entries.status[static_cast<std::size_t>(eBachelor)] != PropagationStatus::kOk) {
        return CascadeStage::kBachelorRefitFailed;
      }
      if (!cascade.fitted) {
        throw std::logic_error("decay finder study: an MC cascade's own candidate of a taken Lambda was never fitted");
      }
      if (d.status == static_cast<std::uint8_t>(VertexFitStatus::kNoStart)) {
        return CascadeStage::kNoStart;
      }
      return d.status == static_cast<std::uint8_t>(VertexFitStatus::kOk) ? CascadeStage::kFound : CascadeStage::kFitFailed;
    };
    d.stage = static_cast<std::uint8_t>(stageOf());
    mCascadesFound += d.stage == static_cast<std::uint8_t>(CascadeStage::kFound) ? 1U : 0U;
    mCascadeDecayWriter->write(d);
  }
  mCascades += cascades.size();
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYFINDERSTUDY_H_
