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

/// \file DecayCensusStudy.h
/// \brief What decays a sample holds for a decay fit, and how many a source
///        of chains (StudyChains.h) gives whole.
///
/// A record per decaying member of the timeframe's decay trees
/// (DecayTrees.h): every member that is no final charged particle and has
/// decay products (DecayCandidates.h) -- a K0S, a Lambda, a Xi, a D0, a
/// converted photon, a resonance. It says where the decay is, how deep its
/// tree goes, what it ends in, and whether the source has a chain for every
/// leaf: a decay is complete when nothing on the way is unseen and every leaf
/// has one. A charged member's own chain, if it has one, is counted too: a
/// tracked Xi or Omega.
///
/// Of a particle's several track chains (clones) the longest counts.

#ifndef GTRK_STUDY_DECAYCENSUSSTUDY_H_
#define GTRK_STUDY_DECAYCENSUSSTUDY_H_

#include "GenericTracking/Study/DecayCandidates.h"
#include "GenericTracking/Study/DecayTrees.h"
#include "GenericTracking/Study/MCTruth.h"
#include "GenericTracking/Study/RecordWriter.h"
#include "GenericTracking/Study/StateTruth.h"
#include "GenericTracking/Study/StudyChains.h"
#include "GenericTracking/Study/TruthEvent.h"

#include <algorithm>
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

/// One decaying member of a decay tree (`DecayCensusStudy`).
struct DecayCensusRecord {
  static constexpr std::string_view kName = "decayCensus";

  std::uint32_t timeframe = 0U;
  std::uint8_t source = 0U; ///< ChainSource
  std::uint64_t particle = 0U;
  std::int32_t pdg = 0;
  std::int32_t motherPdg = 0;  ///< its mother's in the trees, 0 for a root
  std::int32_t rootPdg = 0;    ///< the topmost member above it, itself for a root
  std::uint8_t generation = 0; ///< its ancestors in the trees, 0 for a root
  std::uint8_t process = 0U;   ///< how it was made: a TMCProcess
  std::int8_t charge = 0;      ///< [e]
  float p = 0.F;               ///< at production [GeV]
  float pt = 0.F;
  float eta = 0.F;
  /// Its production point's distance from its collision's vertex [cm], NaN
  /// without the collision.
  float productionOffset = std::numeric_limits<float>::quiet_NaN();
  float flight = 0.F; ///< from production to the decay point [cm]
  float decayR = 0.F; ///< the decay point's distance from the z axis [cm]
  float decayZ = 0.F; ///< [cm]
  std::uint8_t products = 0U;
  std::uint8_t productProcess = 0U; ///< the first product's TMCProcess
  std::uint8_t depth = 0U;          ///< generations of composites below it
  std::uint8_t composites = 0U;     ///< decaying members below it
  std::uint8_t leaves = 0U;         ///< final charged particles
  std::uint8_t neutralUnseen = 0U;
  std::uint8_t chargedUnseen = 0U;
  std::uint8_t leavesWithChain = 0U;
  std::uint16_t leafMinRows = 0U; ///< the fewest rows of a leaf's chain, 0 if a leaf has none
  float leafMinPt = 0.F;          ///< the lowest true pT of a leaf [GeV]
  std::uint16_t ownRows = 0U;     ///< its own chain's rows, 0 without one
  bool closed = false;            ///< nothing unseen on the way to its leaves
  bool complete = false;          ///< closed, and every leaf with a chain

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("source", source);
    visit("particle", particle);
    visit("pdg", pdg);
    visit("motherPdg", motherPdg);
    visit("rootPdg", rootPdg);
    visit("generation", generation);
    visit("process", process);
    visit("charge", charge);
    visit("p", p);
    visit("pt", pt);
    visit("eta", eta);
    visit("productionOffset", productionOffset);
    visit("flight", flight);
    visit("decayR", decayR);
    visit("decayZ", decayZ);
    visit("products", products);
    visit("productProcess", productProcess);
    visit("depth", depth);
    visit("composites", composites);
    visit("leaves", leaves);
    visit("neutralUnseen", neutralUnseen);
    visit("chargedUnseen", chargedUnseen);
    visit("leavesWithChain", leavesWithChain);
    visit("leafMinRows", leafMinRows);
    visit("leafMinPt", leafMinPt);
    visit("ownRows", ownRows);
    visit("closed", closed);
    visit("complete", complete);
  }
};

template <DetectorLayout Det>
class DecayCensusStudy
{
 public:
  explicit DecayCensusStudy(RecordWriter& writer) : mWriter(&writer) { writer.begin<DecayCensusRecord>(); }

  /// \a event's decays, against \a chains from \a source.
  void process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source);

  [[nodiscard]] std::vector<std::string> summary() const
  {
    return {std::format("decay census: {} decaying members into {}, {} closed, {} complete", mRecords, mWriter->path(), mClosed, mComplete)};
  }

 private:
  RecordWriter* mWriter = nullptr;
  std::size_t mRecords = 0U;
  std::size_t mClosed = 0U;
  std::size_t mComplete = 0U;
};

template <DetectorLayout Det>
void DecayCensusStudy<Det>::process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source)
{
  if (event.trees == nullptr) {
    throw std::runtime_error("decay census: the host has no decay trees (a truth dump of another version?)");
  }
  const auto& trees = *event.trees;
  const auto longest = longestChains(chains);
  const auto rowsOf = [&](std::uint64_t particle) -> std::size_t {
    const auto found = longest.find(particle);
    return found == longest.end() ? 0U : found->second->rows.size();
  };
  const auto clampRows = [](std::size_t rows) { return static_cast<std::uint16_t>(std::min<std::size_t>(rows, std::numeric_limits<std::uint16_t>::max())); };
  const auto ptOf = [](const TreeParticle& particle) { return std::hypot(particle.momentum[0], particle.momentum[1]); };

  for (const auto& member : trees.particles) {
    if (isFinalCharged(member.pdg)) {
      continue;
    }
    const auto products = decayProducts(trees, member);
    if (products.empty()) {
      continue;
    }
    const auto shape = decayShape(trees, member, particleCharge);
    DecayCensusRecord record{
      .timeframe = event.timeframe,
      .source = static_cast<std::uint8_t>(source),
      .particle = member.id,
      .pdg = member.pdg,
      .rootPdg = member.pdg,
      .process = member.process,
      .charge = static_cast<std::int8_t>(std::lround(particleCharge(member.pdg))),
      .products = shape.products,
      .productProcess = products.front()->process,
      .depth = shape.depth,
      .composites = static_cast<std::uint8_t>(std::min<std::size_t>(shape.composites.size(), 255U)),
      .leaves = static_cast<std::uint8_t>(std::min<std::size_t>(shape.leaves.size(), 255U)),
      .neutralUnseen = shape.neutralUnseen,
      .chargedUnseen = shape.chargedUnseen,
      .ownRows = clampRows(rowsOf(member.id)),
      .closed = shape.closed(),
    };
    if (const auto* mother = trees.find(member.mother); mother != nullptr) {
      record.motherPdg = mother->pdg;
      for (const auto* above = mother; above != nullptr; above = trees.find(above->mother)) {
        record.rootPdg = above->pdg;
        ++record.generation;
      }
    }
    const auto& momentum = member.momentum;
    const double p = std::hypot(momentum[0], momentum[1], momentum[2]);
    const double pt = ptOf(member);
    record.p = static_cast<float>(p);
    record.pt = static_cast<float>(pt);
    record.eta = static_cast<float>(pt > 0. ? std::asinh(momentum[2] / pt) : std::copysign(std::numeric_limits<double>::infinity(), momentum[2]));
    if (const auto offset = productionOffset(*event.collisions, member.id, member.vertex)) {
      record.productionOffset = static_cast<float>(*offset);
    }
    const auto& decay = products.front()->vertex;
    record.flight = static_cast<float>(std::hypot(decay[0] - member.vertex[0], decay[1] - member.vertex[1], decay[2] - member.vertex[2]));
    record.decayR = static_cast<float>(std::hypot(decay[0], decay[1]));
    record.decayZ = static_cast<float>(decay[2]);

    std::size_t withChain = 0U;
    std::size_t minRows = std::numeric_limits<std::size_t>::max();
    double minPt = std::numeric_limits<double>::infinity();
    for (const auto* leaf : shape.leaves) {
      const auto rows = rowsOf(leaf->id);
      withChain += rows != 0U ? 1U : 0U;
      minRows = std::min(minRows, rows);
      minPt = std::min(minPt, ptOf(*leaf));
    }
    record.leavesWithChain = static_cast<std::uint8_t>(std::min<std::size_t>(withChain, 255U));
    record.leafMinRows = shape.leaves.empty() ? 0U : clampRows(minRows);
    record.leafMinPt = shape.leaves.empty() ? 0.F : static_cast<float>(minPt);
    record.complete = record.closed && withChain == shape.leaves.size();

    mWriter->write(record);
    ++mRecords;
    mClosed += record.closed ? 1U : 0U;
    mComplete += record.complete ? 1U : 0U;
  }
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYCENSUSSTUDY_H_
