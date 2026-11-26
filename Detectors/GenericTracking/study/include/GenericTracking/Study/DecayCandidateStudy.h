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

/// \file DecayCandidateStudy.h
/// \brief The V0 hypotheses fitted as candidates (DecayCandidateFit.h) on the
///        tracker's tracks, the pairs grouped by the MC: a stand-in for a
///        finder.
///
/// The pairs: every MC decay into two final charged particles of opposite
/// charge (DecayCandidates.h's leaves) whose both particles have a good
/// track, each its longest. Each pair is a candidate under every V0
/// hypothesis whose slots' charges its tracks' reconstructed charges fit --
/// K0S -> pi+ pi-, Lambda -> p pi-, anti-Lambda -> anti-p pi+ --, each slot
/// the track's daughter-table entry under the slot's species (no PID: every
/// species the hypotheses take), the start the helices' closest approach
/// (`decayFitStart`), the reference point its collision's mocked primary
/// vertex (PrimaryVertexStudy.h). A record per candidate, with whether its
/// species are the tracks' true ones (`pidMatch`: perfect PID keeps those
/// alone) and whether it is the decay's own hypothesis.

#ifndef GTRK_STUDY_DECAYCANDIDATESTUDY_H_
#define GTRK_STUDY_DECAYCANDIDATESTUDY_H_

#include "GenericTracking/Base/Species.h"
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

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <format>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace o2::gtrk::study
{

/// One candidate's fit (`DecayCandidateStudy`).
struct DecayCandidateRecord {
  static constexpr std::string_view kName = "decayCandidate";
  std::uint32_t timeframe = 0U;
  std::uint64_t particle = 0U;      ///< the MC mother
  std::int32_t pdg = 0;             ///< its species
  std::uint8_t hypothesis = 0U;     ///< 0 K0S, 1 Lambda, 2 anti-Lambda
  std::uint8_t pidMatch = 0U;       ///< every slot's species its track's true one
  std::uint8_t trueHypothesis = 0U; ///< the mother is the hypothesis's
  std::uint8_t status = 0U;         ///< VertexFitStatus
  std::uint8_t iterations = 0U;
  std::uint8_t hasPoint = 0U; ///< its collision's mocked primary vertex was fitted
  float chi2 = 0.F;
  std::int16_t ndf = 0;
  float dX = 0.F; ///< fitted less the MC decay point [cm]
  float dY = 0.F;
  float dZ = 0.F;
  float mass = 0.F; ///< the slots' invariant mass [GeV]
  float sMass = 0.F;
  float length = 0.F; ///< the flight from the point along the momentum [cm]
  float sLength = 0.F;
  float cosPointing = 0.F;
  float microseconds = 0.F; ///< the candidate's fit, its observables with it

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("particle", particle);
    visit("pdg", pdg);
    visit("hypothesis", hypothesis);
    visit("pidMatch", pidMatch);
    visit("trueHypothesis", trueHypothesis);
    visit("status", status);
    visit("iterations", iterations);
    visit("hasPoint", hasPoint);
    visit("chi2", chi2);
    visit("ndf", ndf);
    visit("dX", dX);
    visit("dY", dY);
    visit("dZ", dZ);
    visit("mass", mass);
    visit("sMass", sMass);
    visit("length", length);
    visit("sLength", sLength);
    visit("cosPointing", cosPointing);
    visit("microseconds", microseconds);
  }
};

template <DetectorLayout Det>
class DecayCandidateStudy
{
 public:
  explicit DecayCandidateStudy(RecordWriter& writer) : mWriter(&writer) { writer.begin<DecayCandidateRecord>(); }

  /// \a chains of \a event, from \a source (the tracks only).
  void process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source)
  {
    if (source != ChainSource::kTracks || event.tracks == nullptr || event.trees == nullptr) {
      throw std::runtime_error("decay candidate study: it takes the tracker's tracks and the decay trees");
    }
    const std::array<Species, 2> species{Species::kPion, Species::kProton};
    const auto hypotheses = v0Hypotheses();

    // The collisions' mocked primary vertices.
    const auto reference = referencePointsOf(mockPrimaryVertices(event, chains));

    // The pairs, positive track first.
    const auto chainOf = longestChains(chains);
    struct Pair {
      const TreeParticle* mother = nullptr;
      std::array<std::uint32_t, 2> tracks{}; ///< the positive's, the negative's
      std::array<std::int32_t, 2> pdgs{};
    };
    std::vector<Pair> pairs;
    const auto tracks = event.tracks->tracks().view();
    for (const auto& member : event.trees->particles) {
      const auto products = decayProducts(*event.trees, member);
      if (products.size() != 2U || !isFinalCharged(products[0]->pdg) || !isFinalCharged(products[1]->pdg)) {
        continue;
      }
      const auto first = chainOf.find(products[0]->id);
      const auto second = chainOf.find(products[1]->id);
      if (first == chainOf.end() || second == chainOf.end()) {
        continue;
      }
      const auto& firstContext = first->second->track;
      const auto& secondContext = second->second->track;
      if (!firstContext || !secondContext) {
        continue;
      }
      const auto chargeOf = [&](std::uint32_t t) { return tracks.state[t].qOverP() > 0.F ? 1 : -1; };
      if (chargeOf(firstContext->track) == chargeOf(secondContext->track)) {
        continue;
      }
      const bool firstPositive = chargeOf(firstContext->track) > 0;
      pairs.push_back({.mother = &member, .tracks = {firstPositive ? firstContext->track : secondContext->track, firstPositive ? secondContext->track : firstContext->track}, .pdgs = {firstPositive ? products[0]->pdg : products[1]->pdg, firstPositive ? products[1]->pdg : products[0]->pdg}});
    }

    // Every pair's tracks under every species, refitted once.
    DaughterTable<Det> table(std::span<const Species>{species});
    const auto flushed = std::chrono::steady_clock::now();
    fillDaughterTable(table, event, "decay candidate study", [&](DaughterTable<Det>& asked) {
      for (const auto& pair : pairs) {
        for (const auto t : pair.tracks) {
          for (std::size_t s = 0; s < species.size(); ++s) {
            asked.request(t, static_cast<std::uint8_t>(s));
          }
        }
      }
    });
    mFlushSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - flushed).count();
    mEntries += table.size();
    const auto daughters = daughterSource(table, event);

    // The candidates: a pair under each hypothesis, started at its helices' closest approach.
    std::vector<DecayCandidate> candidates;
    std::vector<DecayCandidateRecord> records;
    for (const auto& pair : pairs) {
      for (std::size_t h = 0; h < hypotheses.size(); ++h) {
        const auto& hypothesis = hypotheses[h];
        DecayCandidate candidate{.hypothesis = static_cast<std::uint16_t>(h), .point = reference.of(collisionOf(pair.mother->id)), .entries = {}, .start = {}};
        bool pidMatch = true;
        for (std::uint32_t slot = 0; slot < 2U; ++slot) {
          const auto entry = table.entryOf(pair.tracks[slot], static_cast<std::uint8_t>(table.speciesIndex(hypothesis.species[slot])));
          candidate.entries.push_back(static_cast<std::uint32_t>(entry));
          const auto truthSpecies = speciesOf(pair.pdgs[slot]);
          pidMatch = pidMatch && truthSpecies && *truthSpecies == hypothesis.species[slot];
        }
        static_cast<void>(startDecayCandidate(candidate, std::span<const DecayHypothesis>{hypotheses}, daughters, static_cast<std::uint32_t>(reference.points.size())));
        candidates.push_back(candidate);
        records.push_back({.timeframe = event.timeframe, .particle = pair.mother->id, .pdg = pair.mother->pdg, .hypothesis = static_cast<std::uint8_t>(h), .pidMatch = static_cast<std::uint8_t>(pidMatch ? 1U : 0U), .trueHypothesis = static_cast<std::uint8_t>(pair.mother->pdg == kV0Mothers[h] ? 1U : 0U)});
      }
    }

    // The fits, each candidate timed with its observables.
    std::vector<DecayCandidateResult> results(candidates.size());
    std::vector<float> microseconds(candidates.size());
    o2::gtrk::detail::parallelFor(*event.arena, 0U, candidates.size(), 16, [&](std::size_t begin, std::size_t end) {
      WorkspaceBuffer workspace;
      for (std::size_t c = begin; c < end; ++c) {
        const auto started = std::chrono::steady_clock::now();
        fitDecayCandidate(candidates[c], std::span<const DecayHypothesis>{hypotheses}, daughters, event.geometry->material(), std::span<const MeasuredPoint>{reference.points}, workspace, results[c]);
        microseconds[c] = static_cast<float>(std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - started).count());
      }
    });
    for (std::size_t c = 0; c < candidates.size(); ++c) {
      auto& record = records[c];
      const auto& result = results[c];
      record.microseconds = microseconds[c];
      writeCandidateFit(record, result);
      if (result.fit.status == VertexFitStatus::kOk) {
        const auto& decay = decayPointOf(*event.trees, *event.trees->find(record.particle));
        record.dX = static_cast<float>(result.fit.vertex[0] - decay[0]);
        record.dY = static_cast<float>(result.fit.vertex[1] - decay[1]);
        record.dZ = static_cast<float>(result.fit.vertex[2] - decay[2]);
        ++mFitted;
      }
      mWriter->write(record);
    }
    mPairs += pairs.size();
    mCandidates += candidates.size();
  }

  [[nodiscard]] std::vector<std::string> summary() const
  {
    return {std::format("decay candidate study: {} pairs, {} candidates, {} fitted, records into {}; daughter table {} entries in {:.3f} s", mPairs, mCandidates, mFitted, mWriter->path(), mEntries, mFlushSeconds)};
  }

 private:
  RecordWriter* mWriter = nullptr;
  std::uint64_t mPairs = 0U;
  std::uint64_t mCandidates = 0U;
  std::uint64_t mFitted = 0U;
  std::uint64_t mEntries = 0U;
  double mFlushSeconds = 0.;
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYCANDIDATESTUDY_H_
