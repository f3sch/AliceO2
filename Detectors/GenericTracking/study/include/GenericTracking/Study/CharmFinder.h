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

/// \file CharmFinder.h
/// \brief The decay finder's charm hadrons (D0 -> K pi, D+ -> K pi pi,
///        Ds+ -> K K pi, Lambda_c+ -> p K pi, and their conjugates): their
///        hypotheses, MC classes, loss stages and records (`decayfinder`,
///        DecayFinderStudy.h, finds them).
///
/// A charm candidate is a combination of a collision's unflagged tracks of
/// the hypothesis's charges: every opposite-charge pair under D0 and its
/// conjugate, every two tracks of one charge with one of the other under
/// D+, Ds+ and Lambda_c+ (or their conjugates), each slot assignment the
/// hypothesis distinguishes (the two pions of a D+ once). One vertex each,
/// fitted on its own. A three-prong is fitted only within the pre-fit mass
/// window (`DecayParam::charmPreMassWindow`: the mass of its tracks'
/// momenta at their first measurements); the two-prongs are all fitted;
/// nothing else is cut. Flights are a few hundred um: the
/// pointing at the collision's point (`pointChi2`, the production there
/// with the point's covariance) is the candidates' handle, written beside
/// the fit.

#ifndef GTRK_STUDY_CHARMFINDER_H_
#define GTRK_STUDY_CHARMFINDER_H_

#include "GenericTracking/Tracking/DecayCandidateFit.h"
#include "GenericTracking/Tracking/Species.h"
#include "GenericTracking/Study/DecayStudyTools.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>
#include <vector>

namespace o2::gtrk::study
{

/// The charm hypotheses, in this order after the finder's others
/// (`charmIndex`): D0 -> K- pi+ and D0bar -> K+ pi- (slots pi+ K-, K+ pi-);
/// D+ -> pi+ pi+ K-, Ds+ -> K+ pi+ K-, Lambda_c+ -> p pi+ K-, each then its
/// conjugate: two slots of the first charge, then one of the other.
inline constexpr std::array<std::int32_t, 8> kCharmMothers{421, -421, 411, -411, 431, -431, 4122, -4122};

/// The weakly decaying charm hadrons a track's charm ancestor is one of.
inline constexpr std::array<std::int32_t, 7> kCharmHadrons{411, 421, 431, 4122, 4132, 4232, 4332};

/// The charm hypotheses (`kCharmMothers`' order).
[[nodiscard]] inline std::vector<DecayHypothesis> charmHypotheses()
{
  const Species pion = Species::kPion;
  const Species kaon = Species::kKaon;
  const Species proton = Species::kProton;
  const auto make = [](std::initializer_list<Species> species, std::initializer_list<std::int8_t> charges, std::int32_t pdg) {
    DecayHypothesis h;
    for (const auto s : species) {
      h.species.push_back(s);
    }
    for (const auto q : charges) {
      h.charges.push_back(q);
    }
    h.massGeV = particleMass(pdg);
    return h;
  };
  return {make({pion, kaon}, {1, -1}, 421), make({kaon, pion}, {1, -1}, -421),
          make({pion, pion, kaon}, {1, 1, -1}, 411), make({pion, pion, kaon}, {-1, -1, 1}, -411),
          make({kaon, pion, kaon}, {1, 1, -1}, 431), make({kaon, pion, kaon}, {-1, -1, 1}, -431),
          make({proton, pion, kaon}, {1, 1, -1}, 4122), make({proton, pion, kaon}, {-1, -1, 1}, -4122)};
}

/// A charm candidate's MC class (`CharmFinderCandidateRecord::truthClass`), the first of:
enum class CharmClass : std::uint8_t {
  kTrue = 1,      ///< an MC charm decay's whole set of tracks, under its own hypothesis and slots
  kReflection,    ///< the same set under another hypothesis or slot assignment (K <-> pi, D+ as Ds+)
  kFake,          ///< a track without a good particle
  kPartial,       ///< tracks of one charm decay that has more (a missing track, a pi0)
  kSomeFromCharm, ///< some tracks of a charm decay, the others not
  kCombinatorial  ///< no charm ancestor among them
};

/// Where an MC charm decay was lost (`CharmFinderDecayRecord::stage`), the first of:
enum class CharmStage : std::uint8_t {
  kFound = 0,       ///< its own candidate was fitted
  kNotFound,        ///< a leaf has no good track of its own charge (`missing`)
  kContributor,     ///< a leaf's track entered its primary vertex's fit
  kRefitFailed,     ///< a leaf's refit under its true species failed
  kOutsideWindow,   ///< a three-prong outside the pre-fit mass window
  kNoStart,         ///< the helices give no start
  kFitFailed        ///< the fit failed
};

/// One charm candidate (`decayfinder`, `-charm`).
struct CharmFinderCandidateRecord {
  static constexpr std::string_view kName = "charmFinderCandidate";
  static constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();
  std::uint32_t timeframe = 0U;
  std::uint64_t collision = 0U;
  std::uint8_t hypothesis = 0U; ///< into `kCharmMothers`
  std::uint8_t truthClass = 0U; ///< CharmClass
  std::uint8_t pidMatch = 0U;   ///< every slot's species its track's true one
  std::uint64_t mother = 0U;    ///< the tracks' common charm ancestor, 0 without one
  std::int32_t motherPdg = 0;
  float s2First = kNaN;  ///< the slots' displacement S^2 against the point
  float s2Second = kNaN;
  float s2Third = kNaN; ///< NaN for a two-prong
  float ptMin = kNaN;   ///< the slots' smallest pT [GeV]
  float preMass = kNaN; ///< from the tracks' momenta at their first measurements [GeV]
  std::uint8_t status = 0U;
  std::uint8_t iterations = 0U;
  std::uint8_t ahead = 0U; ///< its tracks whose first measurement lies behind the fitted vertex
  float chi2 = kNaN;
  std::int16_t ndf = 0;
  float mass = kNaN;
  float sMass = kNaN;
  float radius = kNaN;
  float length = kNaN; ///< the flight from the point [cm]
  float sLength = kNaN;
  float cosPointing = kNaN;
  float pointChi2 = kNaN; ///< what the production at the point (its covariance) adds, 2 dof
  float dXY = kNaN;       ///< kTrue, kReflection: the vertex less the MC decay point [cm]
  float dZ = kNaN;
  float microseconds = 0.F;

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("collision", collision);
    visit("hypothesis", hypothesis);
    visit("truthClass", truthClass);
    visit("pidMatch", pidMatch);
    visit("mother", mother);
    visit("motherPdg", motherPdg);
    visit("s2First", s2First);
    visit("s2Second", s2Second);
    visit("s2Third", s2Third);
    visit("ptMin", ptMin);
    visit("preMass", preMass);
    visit("status", status);
    visit("iterations", iterations);
    visit("ahead", ahead);
    visit("chi2", chi2);
    visit("ndf", ndf);
    visit("mass", mass);
    visit("sMass", sMass);
    visit("radius", radius);
    visit("length", length);
    visit("sLength", sLength);
    visit("cosPointing", cosPointing);
    visit("pointChi2", pointChi2);
    visit("dXY", dXY);
    visit("dZ", dZ);
    visit("microseconds", microseconds);
  }
};

/// One MC charm decay into a hypothesis's final state, through any
/// resonances (K*0, phi, Delta++, Lambda(1520)), nothing unseen
/// (`decayfinder`, `-charmDecays`).
struct CharmFinderDecayRecord {
  static constexpr std::string_view kName = "charmFinderDecay";
  static constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();
  std::uint32_t timeframe = 0U;
  std::uint64_t particle = 0U;
  std::int32_t pdg = 0;
  std::uint8_t hypothesis = 0U;
  float flight = 0.F; ///< the MC decay point's distance from the production point [cm]
  float pt = 0.F;
  float eta = 0.F; ///< the mother's pseudorapidity and rapidity at production (`rapiditiesOf`)
  float y = 0.F;
  float ptMin = 0.F;           ///< its leaves' smallest pT [GeV]
  std::uint8_t depth = 0U;     ///< resonances below it
  std::uint8_t stage = 0U;     ///< CharmStage
  std::uint8_t missing = 0U;   ///< kNotFound: the slots without a good track of their own charge, bit per slot
  std::uint8_t recovered = 0U; ///< some kTrue candidate was fitted, its own or a clone's
  float preMass = kNaN;        ///< its own tracks' pre-fit mass (all tracked and refitted), whatever the window [GeV]
  std::uint8_t status = 0U;    ///< its own candidate's
  float chi2 = kNaN;
  std::int16_t ndf = 0;
  float mass = kNaN;
  float sMass = kNaN;
  float length = kNaN;
  float sLength = kNaN;
  float cosPointing = kNaN;
  float pointChi2 = kNaN;
  float dXY = kNaN;
  float dZ = kNaN;
  /// The fitted vertex less the MC decay point along and across the true
  /// flight, and the fit's sigmas there (`writeVertexResiduals`) [cm].
  float dL = kNaN;
  float sL = kNaN;
  float dT = kNaN;
  float sT = kNaN;
  float dN = kNaN;
  float sN = kNaN;
  float productionOffset = 0.F; ///< its production point's distance from its collision's vertex [cm] (a non-prompt D)
  float mcPointChi2 = kNaN;     ///< its own candidate's production chi2 at the MC collision vertex (exact)
  float qOverPPullMax = kNaN;   ///< its own tracks' refits at their first rows against the MC: the largest |q/p pull|
  float anglePullMax = kNaN;    ///< and the largest |phi| or |theta pull|

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("particle", particle);
    visit("pdg", pdg);
    visit("hypothesis", hypothesis);
    visit("flight", flight);
    visit("pt", pt);
    visit("eta", eta);
    visit("y", y);
    visit("ptMin", ptMin);
    visit("depth", depth);
    visit("stage", stage);
    visit("missing", missing);
    visit("recovered", recovered);
    visit("preMass", preMass);
    visit("status", status);
    visit("chi2", chi2);
    visit("ndf", ndf);
    visit("mass", mass);
    visit("sMass", sMass);
    visit("length", length);
    visit("sLength", sLength);
    visit("cosPointing", cosPointing);
    visit("pointChi2", pointChi2);
    visit("dXY", dXY);
    visit("dZ", dZ);
    visit("dL", dL);
    visit("sL", sL);
    visit("dT", dT);
    visit("sT", sT);
    visit("dN", dN);
    visit("sN", sN);
    visit("productionOffset", productionOffset);
    visit("mcPointChi2", mcPointChi2);
    visit("qOverPPullMax", qOverPPullMax);
    visit("anglePullMax", anglePullMax);
  }
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_CHARMFINDER_H_
