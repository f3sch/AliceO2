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

/// \file DecayTrees.h
/// \brief The MC decay trees the rows' particles belong to: every member's
///        mother, daughters, production point and momentum.
///
/// The particle table (`ParticleTable`) holds the particles that carry rows,
/// each with its mother's species only. A decay fit needs the tree a decay
/// belongs to: the mothers that leave no rows (a Lambda, a D0), the decay
/// points, the mothers' momenta and every daughter, measured or not. The
/// trees hold:
/// - every particle a row carries;
/// - its ancestors, up to the first that is no hadron, lepton, photon or
///   nucleus (a parton, a diquark, a string: the generator's history);
/// - every daughter of those ancestors.
///
/// An ancestor's daughters are all the event's particles that name it as
/// their mother, in the event's order; how each was made (`process`) tells a
/// decay from a hadronic interaction or a conversion. A particle that is no
/// ancestor of a row's particle lists no daughters.

#ifndef GTRK_STUDY_DECAYTREES_H_
#define GTRK_STUDY_DECAYTREES_H_

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <span>
#include <vector>

namespace o2::gtrk::study
{

/// No particle: a root's mother.
inline constexpr std::uint64_t kNoTreeParticle = std::numeric_limits<std::uint64_t>::max();

/// One particle of a decay tree.
struct TreeParticle {
  std::uint64_t id = 0U; ///< its packed identity, as a row's particles (source, event, track)
  std::int32_t pdg = 0;
  std::uint8_t process = 0U; ///< how it was made: a TMCProcess (kPPrimary for the generator's)
  std::uint64_t mother = kNoTreeParticle;
  std::array<double, 3> vertex{};   ///< its production point [cm]
  std::array<double, 3> momentum{}; ///< at production [GeV]
  std::uint32_t firstDaughter = 0U; ///< into `DecayTrees::daughters`
  std::uint32_t daughterCount = 0U;
};

/// The decay trees of one timeframe.
struct DecayTrees {
  std::vector<TreeParticle> particles; ///< sorted by identity
  std::vector<std::uint64_t> daughters;

  /// The member \a id, none if it is not one.
  [[nodiscard]] const TreeParticle* find(std::uint64_t id) const
  {
    const auto found = std::ranges::lower_bound(particles, id, {}, &TreeParticle::id);
    return found == particles.end() || found->id != id ? nullptr : &*found;
  }

  /// \a particle's daughters' identities.
  [[nodiscard]] std::span<const std::uint64_t> daughtersOf(const TreeParticle& particle) const
  {
    return std::span<const std::uint64_t>{daughters}.subspan(particle.firstDaughter, particle.daughterCount);
  }
};

/// Whether a particle of species \a pdg is a member a tree climbs through: a
/// hadron, a lepton, a photon or a nucleus. Quarks, gluons, diquarks, the
/// gauge and Higgs bosons, the generators' own codes (81-100) and the
/// diffractive states (the reggeon 110, the pomeron 990, 99xxxxx) are not.
[[nodiscard]] constexpr bool isTreeSpecies(std::int32_t pdg)
{
  const auto code = pdg < 0 ? -static_cast<std::int64_t>(pdg) : static_cast<std::int64_t>(pdg);
  if (code >= 11 && code <= 16) {
    return true; // leptons
  }
  if (code == 22) {
    return true;
  }
  if (code < 100 || code == 110 || code == 990 || (code >= 9900000 && code < 10000000)) {
    return false;
  }
  // A diquark is 1000 a + 100 b + 10 * 0 + spin: no third quark.
  return code < 1000 || code >= 1000000000 || ((code / 10) % 10) != 0;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYTREES_H_
