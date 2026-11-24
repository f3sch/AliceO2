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

/// \file CascadeFinder.h
/// \brief The decay finder's cascades (Xi, Omega -> Lambda + meson): their
///        hypotheses, MC classes, loss stages and records (`decayfinder`,
///        DecayFinderStudy.h, finds them).
///
/// A cascade candidate is a selected Lambda candidate of the V0 pass
/// (`DecayParam`: its mass window, its legs' S^2, its vertex before both
/// legs' first measurements) with a bachelor: every other unflagged track of
/// the collision of the cascade's charge (a Lambda with a negative, an
/// anti-Lambda with a positive), never one of the Lambda's tracks. Each is
/// fitted under Xi (a pion bachelor) and Omega (a kaon), with the Lambda's
/// mass free and held: one tree fit of its three tracks (the V0 fit only
/// selects). Nothing is cut after the Lambda's selection.

#ifndef GTRK_STUDY_CASCADEFINDER_H_
#define GTRK_STUDY_CASCADEFINDER_H_

#include "GenericTracking/Tracking/DecayCandidateFit.h"
#include "GenericTracking/Tracking/Species.h"
#include "GenericTracking/Study/DecayStudyTools.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>
#include <vector>

namespace o2::gtrk::study
{

/// The finder's hypotheses: the V0s' (`v0Hypotheses`, `kV0Mothers`' order),
/// then the cascades' at `cascadeIndex`.
[[nodiscard]] inline std::size_t cascadeIndex(int sign, int omega, int held)
{
  return kV0Mothers.size() + (std::size_t{4} * static_cast<std::size_t>(held)) + (std::size_t{2} * static_cast<std::size_t>(omega)) + static_cast<std::size_t>(sign);
}

/// The cascade of a Lambda's \a sign (0 the Lambda's, 1 the anti-Lambda's)
/// and \a omega (0 Xi, 1 Omega).
[[nodiscard]] inline std::int32_t cascadePdg(int sign, int omega) { return (omega != 0 ? 3334 : 3312) * (sign == 0 ? 1 : -1); }

/// The V0 hypotheses, then per Lambda sign Xi and Omega, the Lambda's mass
/// free and held (`cascadeIndex`): Xi- -> Lambda pi-, Omega- -> Lambda K-
/// and their conjugates.
[[nodiscard]] inline std::vector<DecayHypothesis> finderHypotheses()
{
  const auto v0 = v0Hypotheses();
  std::vector<DecayHypothesis> hypotheses(v0.begin(), v0.end());
  hypotheses.resize(cascadeIndex(1, 1, 1) + 1U);
  for (int held = 0; held < 2; ++held) {
    for (int omega = 0; omega < 2; ++omega) {
      for (int sign = 0; sign < 2; ++sign) {
        auto& h = hypotheses[cascadeIndex(sign, omega, held)];
        h.composites.push_back(static_cast<std::uint16_t>(1 + sign)); // the Lambda's or the anti-Lambda's V0 hypothesis
        h.species.push_back(omega != 0 ? Species::kKaon : Species::kPion);
        h.charges.push_back(static_cast<std::int8_t>(sign == 0 ? -1 : 1));
        h.massGeV = particleMass(cascadePdg(sign, omega));
        h.holdComposites = held != 0;
      }
    }
  }
  return hypotheses;
}

/// Why a Lambda candidate is no cascade's (`DecayParam`); 0 taken.
enum class LambdaRejection : std::uint8_t {
  kTaken = 0,
  kNotFitted, ///< its V0 fit failed
  kMass,      ///< outside the mass window
  kAhead,     ///< its vertex lies past a leg's first measurement
  kS2         ///< a leg's S^2 below the floor
};

/// A cascade candidate's MC class (`CascadeFinderCandidateRecord::truthClass`), the first of:
enum class CascadeClass : std::uint8_t {
  kTrue = 1,      ///< an MC cascade's Lambda pair and bachelor, under its own hypothesis
  kReflection,    ///< the same under the other (a Xi fitted as Omega, or the reverse)
  kFake,          ///< a track without a good particle
  kCloneBachelor, ///< the bachelor a track of a Lambda leg's particle
  kOtherBachelor, ///< a cascade's Lambda with another bachelor
  kPrimaryLambda, ///< a Lambda from no cascade
  kWrongLambda,   ///< an MC V0's pair under another hypothesis, or a leg's charge wrong
  kOtherMother,   ///< a pair of one mother that is no V0 of the hypotheses
  kCombinatorial  ///< a pair of no common mother
};

/// Where an MC cascade was lost (`CascadeFinderDecayRecord::stage`), the first of:
enum class CascadeStage : std::uint8_t {
  kFound = 0,           ///< its own candidate (Lambda's mass held) was fitted
  kNotFound,            ///< a leaf has no good track of its own charge (`missing`, `reversed`)
  kContributor,         ///< a leaf's track entered its primary vertex's fit
  kLambdaRefitFailed,   ///< a Lambda leg's refit under its true species failed
  kLambdaFitFailed,     ///< the Lambda's own V0 fit failed
  kLambdaRejected,      ///< its own Lambda candidate not taken (`lambdaRejection`)
  kBachelorRefitFailed, ///< the bachelor's refit under its true species failed
  kNoStart,             ///< the tree gives no start
  kFitFailed            ///< the fit failed
};

/// One cascade candidate (`decayfinder`, a third file, `-cascades`).
struct CascadeFinderCandidateRecord {
  static constexpr std::string_view kName = "cascadeFinderCandidate";
  static constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();
  std::uint32_t timeframe = 0U;
  std::uint64_t collision = 0U;
  std::uint8_t sign = 0U;  ///< 0 a Lambda's (negative bachelor), 1 an anti-Lambda's
  std::uint8_t omega = 0U; ///< 0 Xi, 1 Omega
  std::uint8_t held = 0U;  ///< the Lambda's mass held
  std::uint8_t truthClass = 0U; ///< CascadeClass
  std::uint8_t lambdaClass = 0U; ///< its Lambda pair's FinderClass (DecayFinderStudy.h)
  std::uint8_t pidMatch = 0U;   ///< all three tracks' species their true ones
  std::uint8_t primaries = 0U;  ///< its tracks whose particle is a physical primary
  std::uint64_t cascade = 0U;   ///< the MC cascade of kTrue and kReflection, 0 otherwise
  float s2Positive = kNaN;      ///< the Lambda legs' S^2 (the V0 pass's)
  float s2Negative = kNaN;
  float s2Bachelor = kNaN;
  float lambdaV0Mass = kNaN; ///< the V0 fit's Lambda mass and sigma
  float lambdaV0SMass = kNaN;
  std::uint8_t status = 0U; ///< VertexFitStatus
  std::uint8_t iterations = 0U;
  std::uint8_t ahead = 0U; ///< its tracks whose first measurement lies behind their fitted vertex
  float chi2 = kNaN;
  std::int16_t ndf = 0;
  float mass = kNaN; ///< the three tracks' [GeV]
  float sMass = kNaN;
  float lambdaMass = kNaN; ///< the Lambda's in the tree fit
  float sLambdaMass = kNaN;
  float radius = kNaN;       ///< the cascade's vertex's transverse radius [cm]
  float lambdaRadius = kNaN; ///< the Lambda's
  float length = kNaN;       ///< the cascade's flight from the point [cm]
  float sLength = kNaN;
  float cosPointing = kNaN;
  float lambdaLength = kNaN; ///< the Lambda's flight from the cascade's vertex [cm]
  float sLambdaLength = kNaN;
  float lambdaCosPointing = kNaN;
  float dXY = kNaN; ///< kTrue and kReflection: the cascade's vertex less its MC decay point [cm]
  float dZ = kNaN;
  /// Lambda's mass held: the mother's own rows found (MotherSearch.h) and,
  /// with any, the fit with them (`tracked*`).
  std::uint8_t motherLayers = 0U; ///< the layers searched
  std::uint8_t motherRows = 0U;
  std::uint8_t motherFree = 0U;   ///< of those, rows no committed track holds
  std::uint8_t motherTrue = 0U;   ///< of those, rows of the MC cascade (kTrue, kReflection)
  float motherChi2 = kNaN;        ///< the rows' update chi2 summed
  std::uint8_t trackedStatus = 0U;
  float trackedChi2 = kNaN;
  std::int16_t trackedNdf = 0;
  float trackedMass = kNaN;
  float trackedSMass = kNaN;
  float trackedRadius = kNaN;
  float trackedLength = kNaN;
  float trackedSLength = kNaN;
  float trackedCosPointing = kNaN;
  float trackedDXY = kNaN;
  float trackedDZ = kNaN;
  /// Lambda's mass held, with a point: the chi2 the mother's production
  /// there adds (2 dof; `fitDecayCandidate`'s production), without and with
  /// its own rows -- how well it points back.
  float pointChi2 = kNaN;
  float trackedPointChi2 = kNaN;
  float trackedMCPointChi2 = kNaN; ///< kTrue, kReflection: the same at the MC collision's vertex
  float microseconds = 0.F;

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("collision", collision);
    visit("sign", sign);
    visit("omega", omega);
    visit("held", held);
    visit("truthClass", truthClass);
    visit("lambdaClass", lambdaClass);
    visit("pidMatch", pidMatch);
    visit("primaries", primaries);
    visit("cascade", cascade);
    visit("s2Positive", s2Positive);
    visit("s2Negative", s2Negative);
    visit("s2Bachelor", s2Bachelor);
    visit("lambdaV0Mass", lambdaV0Mass);
    visit("lambdaV0SMass", lambdaV0SMass);
    visit("status", status);
    visit("iterations", iterations);
    visit("ahead", ahead);
    visit("chi2", chi2);
    visit("ndf", ndf);
    visit("mass", mass);
    visit("sMass", sMass);
    visit("lambdaMass", lambdaMass);
    visit("sLambdaMass", sLambdaMass);
    visit("radius", radius);
    visit("lambdaRadius", lambdaRadius);
    visit("length", length);
    visit("sLength", sLength);
    visit("cosPointing", cosPointing);
    visit("lambdaLength", lambdaLength);
    visit("sLambdaLength", sLambdaLength);
    visit("lambdaCosPointing", lambdaCosPointing);
    visit("dXY", dXY);
    visit("dZ", dZ);
    visit("motherLayers", motherLayers);
    visit("motherRows", motherRows);
    visit("motherFree", motherFree);
    visit("motherTrue", motherTrue);
    visit("motherChi2", motherChi2);
    visit("trackedStatus", trackedStatus);
    visit("trackedChi2", trackedChi2);
    visit("trackedNdf", trackedNdf);
    visit("trackedMass", trackedMass);
    visit("trackedSMass", trackedSMass);
    visit("trackedRadius", trackedRadius);
    visit("trackedLength", trackedLength);
    visit("trackedSLength", trackedSLength);
    visit("trackedCosPointing", trackedCosPointing);
    visit("trackedDXY", trackedDXY);
    visit("trackedDZ", trackedDZ);
    visit("pointChi2", pointChi2);
    visit("trackedPointChi2", trackedPointChi2);
    visit("trackedMCPointChi2", trackedMCPointChi2);
    visit("microseconds", microseconds);
  }
};

/// One MC cascade decaying into a Lambda (into its proton and pion) and a
/// charged pion (Xi) or kaon (Omega) (`decayfinder`, `-cascadeDecays`).
struct CascadeFinderDecayRecord {
  static constexpr std::string_view kName = "cascadeFinderDecay";
  static constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();
  std::uint32_t timeframe = 0U;
  std::uint64_t particle = 0U;
  std::int32_t pdg = 0;
  float decayR = 0.F;       ///< the MC decay point's transverse radius [cm]
  float lambdaDecayR = 0.F; ///< the Lambda's
  float pt = 0.F;           ///< at production [GeV]
  float eta = 0.F;          ///< the mother's pseudorapidity and rapidity at production (`rapiditiesOf`)
  float y = 0.F;
  float ptPositive = 0.F;   ///< the Lambda's daughters' and the bachelor's at production [GeV]
  float ptNegative = 0.F;
  float ptBachelor = 0.F;
  std::uint8_t stage = 0U;           ///< CascadeStage
  std::uint8_t missing = 0U;         ///< kNotFound: the leaves without a good track of their own charge, bits 1 the Lambda's positive, 2 its negative, 4 the bachelor
  std::uint8_t reversed = 0U;        ///< of those, the leaves with good tracks of the other charge alone, the same bits
  std::uint8_t lambdaRejection = 0U; ///< LambdaRejection of its own Lambda candidate
  std::uint8_t hasPoint = 0U;        ///< its collision has a mocked primary vertex
  std::uint8_t recovered = 0U;       ///< some kTrue candidate (Lambda's mass held) was fitted, its own or a clone's
  float s2Bachelor = kNaN;
  std::uint8_t status = 0U;     ///< its own candidate's (Lambda's mass held) VertexFitStatus, 0 without one
  std::uint8_t statusFree = 0U; ///< and with the Lambda's mass free
  float chi2 = kNaN;
  std::int16_t ndf = 0;
  float mass = kNaN;
  float sMass = kNaN;
  float massFree = kNaN; ///< the Lambda's mass free
  float sMassFree = kNaN;
  float lambdaMass = kNaN; ///< free: the Lambda's in the tree fit
  float sLambdaMass = kNaN;
  float length = kNaN;
  float sLength = kNaN;
  float cosPointing = kNaN;
  float dXY = kNaN; ///< the fitted cascade vertex less the MC decay point [cm]
  float dZ = kNaN;
  /// The fitted vertex less the MC decay point along and across the true
  /// flight, and the fit's sigmas there (`writeVertexResiduals`) [cm].
  float dL = kNaN;
  float sL = kNaN;
  float dT = kNaN;
  float sT = kNaN;
  float dN = kNaN;
  float sN = kNaN;
  std::uint8_t mcRows = 0U;     ///< the cascade's own rows in the MC (its chain's)
  std::uint8_t motherRows = 0U; ///< its own candidate's (held): rows found, and of those its own
  std::uint8_t motherTrue = 0U;
  std::uint8_t trackedStatus = 0U;
  float trackedMass = kNaN;
  float trackedSMass = kNaN;
  float trackedDXY = kNaN;
  float trackedDZ = kNaN;

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("particle", particle);
    visit("pdg", pdg);
    visit("decayR", decayR);
    visit("lambdaDecayR", lambdaDecayR);
    visit("pt", pt);
    visit("eta", eta);
    visit("y", y);
    visit("ptPositive", ptPositive);
    visit("ptNegative", ptNegative);
    visit("ptBachelor", ptBachelor);
    visit("stage", stage);
    visit("missing", missing);
    visit("reversed", reversed);
    visit("lambdaRejection", lambdaRejection);
    visit("hasPoint", hasPoint);
    visit("recovered", recovered);
    visit("s2Bachelor", s2Bachelor);
    visit("status", status);
    visit("statusFree", statusFree);
    visit("chi2", chi2);
    visit("ndf", ndf);
    visit("mass", mass);
    visit("sMass", sMass);
    visit("massFree", massFree);
    visit("sMassFree", sMassFree);
    visit("lambdaMass", lambdaMass);
    visit("sLambdaMass", sLambdaMass);
    visit("length", length);
    visit("sLength", sLength);
    visit("cosPointing", cosPointing);
    visit("dXY", dXY);
    visit("dZ", dZ);
    visit("dL", dL);
    visit("sL", sL);
    visit("dT", dT);
    visit("sT", sT);
    visit("dN", dN);
    visit("sN", sN);
    visit("mcRows", mcRows);
    visit("motherRows", motherRows);
    visit("motherTrue", motherTrue);
    visit("trackedStatus", trackedStatus);
    visit("trackedMass", trackedMass);
    visit("trackedSMass", trackedSMass);
    visit("trackedDXY", trackedDXY);
    visit("trackedDZ", trackedDZ);
  }
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_CASCADEFINDER_H_
