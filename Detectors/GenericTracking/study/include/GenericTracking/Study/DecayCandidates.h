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

/// \file DecayCandidates.h
/// \brief A decay of the MC trees (DecayTrees.h) as a fit sees it: the
///        charged particles it ends in, which the tracker measures, and what
///        keeps it from being fitted whole.
///
/// A member's decay products are the daughters its decay made: in transport
/// (kPDecay) or by the generator (kPPrimary); a photon's are its conversion's
/// pair (kPPair). A hadronic interaction's products are none.
///
/// A decay is followed down to its final charged particles -- electrons,
/// muons, pions, kaons, protons and the light nuclei d, t, 3He, 4He, which
/// cross the tracker before they decay, if they do -- through every other
/// member's decay products: a Lambda's in a cascade, a Xi's in a charm
/// baryon's. Those are the decay's leaves. A member on the way that is no
/// leaf and has no decay products in the trees -- a neutrino, a neutron, a
/// photon that did not convert, a Lambda that left the tracker, or one whose
/// products carry no rows, as the trees list no daughters of a particle that
/// is no row's ancestor -- is one the fit cannot see: the decay is
/// incomplete.
///
/// No combinatorics: the truth groups the particles.

#ifndef GTRK_STUDY_DECAYCANDIDATES_H_
#define GTRK_STUDY_DECAYCANDIDATES_H_

#include "GenericTracking/Study/DecayTrees.h"

#include <cstdint>
#include <cstdlib>
#include <vector>

namespace o2::gtrk::study
{

/// Whether species \a pdg is a final charged particle of a decay: one the
/// tracker measures, a decay's leaf.
[[nodiscard]] constexpr bool isFinalCharged(std::int32_t pdg)
{
  switch (std::abs(pdg)) {
    case 11:         // e
    case 13:         // mu
    case 211:        // pi
    case 321:        // K
    case 2212:       // p
    case 1000010020: // d
    case 1000010030: // t
    case 1000020030: // 3He
    case 1000020040: // 4He
      return true;
    default:
      return false;
  }
}

/// TMCProcess codes a decay's products are made by (TMCProcess.h).
inline constexpr std::uint8_t kProcessPrimary = 0U; ///< kPPrimary: the generator's
inline constexpr std::uint8_t kProcessDecay = 4U;   ///< kPDecay
inline constexpr std::uint8_t kProcessPair = 5U;    ///< kPPair

/// Whether \a daughter is one of \a mother's decay products.
[[nodiscard]] constexpr bool isDecayProduct(const TreeParticle& mother, const TreeParticle& daughter)
{
  if (mother.pdg == 22) {
    return daughter.process == kProcessPair;
  }
  return daughter.process == kProcessDecay || daughter.process == kProcessPrimary;
}

/// \a particle's decay products, members of \a trees, in its daughters' order.
[[nodiscard]] inline std::vector<const TreeParticle*> decayProducts(const DecayTrees& trees, const TreeParticle& particle)
{
  std::vector<const TreeParticle*> products;
  for (const auto id : trees.daughtersOf(particle)) {
    const auto* daughter = trees.find(id);
    if (daughter != nullptr && isDecayProduct(particle, *daughter)) {
      products.push_back(daughter);
    }
  }
  return products;
}

/// A decay down to its leaves.
struct DecayShape {
  std::vector<const TreeParticle*> leaves;     ///< the final charged particles, depth first
  std::vector<const TreeParticle*> composites; ///< the decaying members below the root, depth first
  std::uint8_t products = 0U;                  ///< the root's own decay products
  std::uint8_t depth = 0U;                     ///< the most generations of composites below the root
  std::uint8_t neutralUnseen = 0U;             ///< neutral members, no leaf, without decay products
  std::uint8_t chargedUnseen = 0U;             ///< charged ones (charge from \a chargeOf)
  /// Whether nothing on the way keeps the fit from seeing it whole.
  [[nodiscard]] bool closed() const { return neutralUnseen == 0U && chargedUnseen == 0U && !leaves.empty(); }
};

/// \a root's decay down to its leaves; \a chargeOf(pdg) a species' charge
/// [e].
template <typename ChargeOf>
[[nodiscard]] DecayShape decayShape(const DecayTrees& trees, const TreeParticle& root, const ChargeOf& chargeOf)
{
  DecayShape shape;
  struct Visit {
    const TreeParticle* particle;
    std::uint8_t generation;
  };
  std::vector<Visit> stack;
  const auto rootProducts = decayProducts(trees, root);
  shape.products = static_cast<std::uint8_t>(rootProducts.size());
  for (auto product = rootProducts.rbegin(); product != rootProducts.rend(); ++product) {
    stack.push_back({*product, 1U});
  }
  while (!stack.empty()) {
    const auto [particle, generation] = stack.back();
    stack.pop_back();
    if (isFinalCharged(particle->pdg)) {
      shape.leaves.push_back(particle);
      continue;
    }
    const auto products = decayProducts(trees, *particle);
    if (products.empty()) {
      if (chargeOf(particle->pdg) == 0.) {
        ++shape.neutralUnseen;
      } else {
        ++shape.chargedUnseen;
      }
      continue;
    }
    shape.composites.push_back(particle);
    if (generation > shape.depth) {
      shape.depth = generation;
    }
    for (auto product = products.rbegin(); product != products.rend(); ++product) {
      stack.push_back({*product, static_cast<std::uint8_t>(generation + 1U)});
    }
  }
  return shape;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYCANDIDATES_H_
