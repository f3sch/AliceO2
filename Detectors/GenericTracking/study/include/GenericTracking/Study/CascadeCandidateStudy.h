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

/// \file CascadeCandidateStudy.h
/// \brief Cascades (Xi, Omega) fitted as candidates of composite hypotheses
///        (DecayCandidateFit.h) on the tracker's tracks: the Lambdas grouped
///        by the MC, the bachelors every track of their collision.
///
/// The Lambdas: every MC Lambda and anti-Lambda decaying into its proton
/// and pion, both with a good track (each its longest) of the right
/// reconstructed charge. The bachelors: every other track associated to the
/// Lambda's collision (`mockAssociation`), no contributor of its primary
/// vertex, of the cascade's charge, under the pion and the kaon: Xi- -> Lambda pi-,
/// Omega- -> Lambda K- and their conjugates, each without and with the
/// Lambda's mass held. Each candidate is one tree fit of its three tracks
/// (the Lambda's vertex and the cascade's), started at the tree's closest
/// approaches (`decayFitStart`), its reference point its collision's mocked
/// primary vertex. Nothing is cut: every candidate is written, the
/// bachelor's displacement S^2 beside its fit and observables.

#ifndef GTRK_STUDY_CASCADECANDIDATESTUDY_H_
#define GTRK_STUDY_CASCADECANDIDATESTUDY_H_

#include "GenericTracking/Tracking/DaughterTable.h"
#include "GenericTracking/Tracking/DecayCandidateFit.h"
#include "GenericTracking/Study/DecayCandidates.h"
#include "GenericTracking/Study/DecayStudyTools.h"
#include "GenericTracking/Study/MCTruth.h"
#include "GenericTracking/Study/PrimaryVertexStudy.h"
#include "GenericTracking/Study/RecordWriter.h"
#include "GenericTracking/Study/StudyChains.h"
#include "GenericTracking/Study/TruthEvent.h"
#include "GenericTracking/Study/TruthFit.h"

#include "ReconstructionDataFormats/PID.h"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <format>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace o2::gtrk::study
{

/// A cascade candidate's MC class (`CascadeCandidateRecord::truthClass`).
enum class CascadeClass : std::uint8_t {
  kTrue = 1,        ///< an MC cascade's Lambda and bachelor under its own hypothesis
  kOtherHypothesis, ///< an MC cascade's Lambda and bachelor under the other (a Xi fitted as Omega)
  kOtherBachelor,   ///< a cascade's Lambda with another track
  kPrimaryLambda    ///< a Lambda from no cascade
};

/// One cascade candidate (`CascadeCandidateStudy`).
struct CascadeCandidateRecord {
  static constexpr std::string_view kName = "cascadeCandidate";
  std::uint32_t timeframe = 0U;
  std::uint64_t lambda = 0U;                                  ///< the MC Lambda
  std::uint64_t cascade = 0U;                                 ///< its MC mother if it is a cascade, 0 otherwise
  std::int32_t cascadePdg = 0;                                ///< that mother's species, 0 for none
  std::uint8_t hypothesis = 0U;                               ///< 0 Xi, 1 Omega (the charge the Lambda's)
  std::uint8_t held = 0U;                                     ///< the Lambda's mass held
  std::uint8_t truthClass = 0U;                               ///< CascadeClass
  std::uint8_t pidMatch = 0U;                                 ///< the bachelor's species its true one
  std::uint8_t bachelorFake = 0U;                             ///< the bachelor's track has no good particle
  float s2Bachelor = std::numeric_limits<float>::quiet_NaN(); ///< its displacement S^2 against the point
  std::uint8_t status = 0U;                                   ///< VertexFitStatus
  std::uint8_t iterations = 0U;
  float chi2 = 0.F;
  std::int16_t ndf = 0;
  float mass = 0.F; ///< the three tracks' invariant mass [GeV]
  float sMass = 0.F;
  float lambdaMass = 0.F; ///< the Lambda's in the tree fit
  float sLambdaMass = 0.F;
  float radius = 0.F;       ///< the cascade's vertex's transverse radius [cm]
  float lambdaRadius = 0.F; ///< the Lambda's
  float length = 0.F;       ///< the cascade's flight from the point [cm]
  float sLength = 0.F;
  float cosPointing = 0.F;
  float lambdaLength = 0.F; ///< the Lambda's flight from the cascade's vertex [cm]
  float sLambdaLength = 0.F;
  float lambdaCosPointing = 0.F;
  float dXY = std::numeric_limits<float>::quiet_NaN(); ///< a true cascade's vertex less its MC decay point [cm]
  float dZ = std::numeric_limits<float>::quiet_NaN();
  float microseconds = 0.F;

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("lambda", lambda);
    visit("cascade", cascade);
    visit("cascadePdg", cascadePdg);
    visit("hypothesis", hypothesis);
    visit("held", held);
    visit("truthClass", truthClass);
    visit("pidMatch", pidMatch);
    visit("bachelorFake", bachelorFake);
    visit("s2Bachelor", s2Bachelor);
    visit("status", status);
    visit("iterations", iterations);
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
    visit("microseconds", microseconds);
  }
};

template <DetectorLayout Det>
class CascadeCandidateStudy
{
 public:
  explicit CascadeCandidateStudy(RecordWriter& writer) : mWriter(&writer) { writer.begin<CascadeCandidateRecord>(); }

  /// \a event's tracks (the source's chains: the Lambdas' daughters and the
  /// primaries of the mocked vertices).
  void process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source);

  [[nodiscard]] std::vector<std::string> summary() const
  {
    return {std::format("cascade candidate study: {} Lambdas, {} candidates, {} fitted, records into {}", mLambdas, mCandidates, mFitted, mWriter->path())};
  }

 private:
  RecordWriter* mWriter = nullptr;
  std::size_t mLambdas = 0U;
  std::size_t mCandidates = 0U;
  std::size_t mFitted = 0U;
};

template <DetectorLayout Det>
void CascadeCandidateStudy<Det>::process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source)
{
  using o2::track::PID;
  if (source != ChainSource::kTracks || event.tracks == nullptr || event.trees == nullptr) {
    throw std::runtime_error("cascade candidate study: it takes the tracker's tracks and the decay trees");
  }
  const auto& store = *event.tracks;
  const auto& geometry = *event.geometry;
  const auto& trees = *event.trees;
  const auto committed = store.tracks().view();

  // The hypotheses: the two Lambdas, then per Lambda sign Xi and Omega,
  // without and with the Lambda's mass held.
  const std::array<PID, 3> species{PID{PID::Pion}, PID{PID::Kaon}, PID{PID::Proton}};
  std::vector<DecayHypothesis> hypotheses(10);
  hypotheses[0] = twoProngHypothesis(PID{PID::Proton}, PID{PID::Pion}, particleMass(3122));
  hypotheses[1] = twoProngHypothesis(PID{PID::Pion}, PID{PID::Proton}, particleMass(-3122));
  // index 2 + 4 held + 2 omega + sign: sign 0 the Lambda's (negative bachelor), 1 the anti-Lambda's
  const auto cascadeIndex = [](int sign, int omega, int held) { return std::size_t{2} + (std::size_t{4} * static_cast<std::size_t>(held)) + (std::size_t{2} * static_cast<std::size_t>(omega)) + static_cast<std::size_t>(sign); };
  const auto cascadePdg = [](int sign, int omega) { return (omega != 0 ? 3334 : 3312) * (sign == 0 ? 1 : -1); };
  for (int held = 0; held < 2; ++held) {
    for (int omega = 0; omega < 2; ++omega) {
      for (int sign = 0; sign < 2; ++sign) {
        auto& h = hypotheses[cascadeIndex(sign, omega, held)];
        h.composites.push_back(static_cast<std::uint16_t>(sign));
        h.species.push_back(omega != 0 ? PID{PID::Kaon} : PID{PID::Pion});
        h.charges.push_back(static_cast<std::int8_t>(sign == 0 ? -1 : 1));
        h.massGeV = particleMass(cascadePdg(sign, omega));
        h.holdComposites = held != 0;
      }
    }
  }

  // The tracks associated to the mocked primary vertices; a vertex's
  // contributors are no bachelors.
  const auto association = mockAssociation(event);
  const auto& points = association.vertices.points;
  const auto& particleOf = association.particle;
  const auto& fake = association.fake;
  const auto n = static_cast<std::uint32_t>(store.size());
  std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> byCollision;
  for (std::uint32_t t = 0; t < n; ++t) {
    if (particleOf[t] != kNoTrackParticle && !association.contributor(t)) {
      byCollision[collisionOf(particleOf[t])].push_back(t);
    }
  }

  // The MC Lambdas with both daughters found, of the right charges.
  const auto chainOf = longestChains(chains);
  struct Lambda {
    const TreeParticle* particle = nullptr;
    int sign = 0;                          ///< 0 Lambda, 1 anti-Lambda
    std::array<std::uint32_t, 2> tracks{}; ///< the positive's, the negative's
  };
  std::vector<Lambda> lambdas;
  for (const auto& member : trees.particles) {
    if (std::abs(member.pdg) != 3122) {
      continue;
    }
    const auto products = decayProducts(trees, member);
    if (products.size() != 2U) {
      continue;
    }
    const int sign = member.pdg > 0 ? 0 : 1;
    std::array<std::uint32_t, 2> tracks{};
    bool usable = true;
    for (const auto* product : products) {
      const auto chain = chainOf.find(product->id);
      if (chain == chainOf.end()) {
        usable = false;
        break;
      }
      const auto& context = chain->second->track;
      if (!context) {
        usable = false;
        break;
      }
      const auto t = context->track;
      const int slot = committed.state[t].qOverP() > 0.F ? 0 : 1;
      const auto truth = truthSpecies(product->pdg);
      if (!truth || truth->pid != hypotheses[static_cast<std::size_t>(sign)].species[static_cast<std::size_t>(slot)] || (particleCharge(product->pdg) > 0.) != (slot == 0)) {
        usable = false;
        break;
      }
      tracks[static_cast<std::size_t>(slot)] = t;
    }
    if (usable && tracks[0] != tracks[1]) {
      lambdas.push_back({.particle = &member, .sign = sign, .tracks = tracks});
    }
  }

  // The table: the Lambdas' tracks under their species, every bachelor as a pion and a kaon.
  DaughterTable<Det> table(std::span<const PID>{species});
  fillDaughterTable(table, event, "cascade candidate study", [&](DaughterTable<Det>& asked) {
    for (const auto& lambda : lambdas) {
      const auto& h = hypotheses[static_cast<std::size_t>(lambda.sign)];
      asked.request(lambda.tracks[0], static_cast<std::uint8_t>(asked.speciesIndex(h.species[0])));
      asked.request(lambda.tracks[1], static_cast<std::uint8_t>(asked.speciesIndex(h.species[1])));
      for (const auto t : byCollision[collisionOf(lambda.particle->id)]) {
        asked.request(t, static_cast<std::uint8_t>(asked.speciesIndex(PID{PID::Pion})));
        asked.request(t, static_cast<std::uint8_t>(asked.speciesIndex(PID{PID::Kaon})));
      }
    }
  });
  const auto entries = table.view();
  const auto daughters = daughterSource(table, event);

  // The candidates, Lambda by Lambda.
  std::vector<std::vector<CascadeCandidateRecord>> perLambda(lambdas.size());
  o2::gtrk::detail::parallelFor(*event.arena, 0U, lambdas.size(), 1, [&](std::size_t begin, std::size_t end) {
    WorkspaceBuffer workspace;
    for (std::size_t l = begin; l < end; ++l) {
      const auto& lambda = lambdas[l];
      const auto collision = collisionOf(lambda.particle->id);
      const auto point = association.vertices.of(collision);
      const auto* mother = trees.find(lambda.particle->mother);
      const bool fromCascade = mother != nullptr && (std::abs(mother->pdg) == 3312 || std::abs(mother->pdg) == 3334);
      const auto& lambdaHypothesis = hypotheses[static_cast<std::size_t>(lambda.sign)];
      const auto eProton = table.entryOf(lambda.tracks[0], static_cast<std::uint8_t>(table.speciesIndex(lambdaHypothesis.species[0])));
      const auto eNegative = table.entryOf(lambda.tracks[1], static_cast<std::uint8_t>(table.speciesIndex(lambdaHypothesis.species[1])));
      for (const auto bachelor : byCollision.at(collision)) {
        if (bachelor == lambda.tracks[0] || bachelor == lambda.tracks[1] || (committed.state[bachelor].qOverP() > 0.F ? 1 : 0) != lambda.sign) {
          continue;
        }
        const bool trueBachelor = fromCascade && fake[bachelor] == 0U && trees.find(particleOf[bachelor]) != nullptr && trees.find(particleOf[bachelor])->mother == mother->id;
        for (int omega = 0; omega < 2; ++omega) {
          const PID bachelorPid = omega != 0 ? PID{PID::Kaon} : PID{PID::Pion};
          const auto eBachelor = table.entryOf(bachelor, static_cast<std::uint8_t>(table.speciesIndex(bachelorPid)));
          const auto bachelorTruth = truthSpecies(event.particles->at(particleOf[bachelor]).pdg);
          float s2 = std::numeric_limits<float>::quiet_NaN();
          if (point >= 0 && entries.status[static_cast<std::size_t>(eBachelor)] == PropagationStatus::kOk) {
            s2 = static_cast<float>(displacement(daughters, static_cast<std::uint32_t>(eBachelor), points[static_cast<std::size_t>(point)]));
          }
          for (int held = 0; held < 2; ++held) {
            const auto h = cascadeIndex(lambda.sign, omega, held);
            CascadeCandidateRecord record{.timeframe = event.timeframe, .lambda = lambda.particle->id, .cascade = fromCascade ? mother->id : 0U, .cascadePdg = fromCascade ? mother->pdg : 0, .hypothesis = static_cast<std::uint8_t>(omega), .held = static_cast<std::uint8_t>(held), .bachelorFake = fake[bachelor], .s2Bachelor = s2};
            record.pidMatch = static_cast<std::uint8_t>(bachelorTruth && bachelorTruth->pid == bachelorPid ? 1U : 0U);
            if (!fromCascade) {
              record.truthClass = static_cast<std::uint8_t>(CascadeClass::kPrimaryLambda);
            } else if (!trueBachelor) {
              record.truthClass = static_cast<std::uint8_t>(CascadeClass::kOtherBachelor);
            } else {
              record.truthClass = static_cast<std::uint8_t>(mother->pdg == cascadePdg(lambda.sign, omega) ? CascadeClass::kTrue : CascadeClass::kOtherHypothesis);
            }
            DecayCandidate candidate{.hypothesis = static_cast<std::uint16_t>(h), .point = point};
            candidate.entries.push_back(static_cast<std::uint32_t>(eProton));
            candidate.entries.push_back(static_cast<std::uint32_t>(eNegative));
            candidate.entries.push_back(static_cast<std::uint32_t>(eBachelor));
            const auto started = std::chrono::steady_clock::now();
            static_cast<void>(startDecayCandidate(candidate, std::span<const DecayHypothesis>{hypotheses}, daughters, points.size()));
            DecayCandidateResult result;
            fitDecayCandidate(candidate, std::span<const DecayHypothesis>{hypotheses}, daughters, geometry.material(), std::span<const MeasuredPoint>{points}, workspace, result);
            record.microseconds = static_cast<float>(std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - started).count());
            writeCandidateFit(record, result);
            if (result.fit.status == VertexFitStatus::kOk) {
              const auto& fit = result.fit;
              record.lambdaMass = static_cast<float>(result.composites[0].mass);
              record.sLambdaMass = static_cast<float>(result.composites[0].sMass);
              record.radius = static_cast<float>(std::hypot(fit.vertices[0][0], fit.vertices[0][1]));
              record.lambdaRadius = static_cast<float>(std::hypot(fit.vertices[1][0], fit.vertices[1][1]));
              record.lambdaLength = static_cast<float>(result.compositeFlights[0].length);
              record.sLambdaLength = static_cast<float>(result.compositeFlights[0].sLength);
              record.lambdaCosPointing = static_cast<float>(result.compositeFlights[0].cosPointing);
              if (record.truthClass == static_cast<std::uint8_t>(CascadeClass::kTrue)) {
                const auto& decayPoint = decayPointOf(trees, *mother);
                record.dXY = static_cast<float>(std::hypot(fit.vertex[0] - decayPoint[0], fit.vertex[1] - decayPoint[1]));
                record.dZ = static_cast<float>(fit.vertex[2] - decayPoint[2]);
              }
            }
            perLambda[l].push_back(record);
          }
        }
      }
    }
  });
  for (const auto& records : perLambda) {
    for (const auto& record : records) {
      mWriter->write(record);
      mFitted += record.status == static_cast<std::uint8_t>(VertexFitStatus::kOk) ? 1U : 0U;
    }
    mCandidates += records.size();
  }
  mLambdas += lambdas.size();
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_CASCADECANDIDATESTUDY_H_
