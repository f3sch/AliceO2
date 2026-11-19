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

/// \file testDecayCandidates.cxx
/// \brief A decay followed to its final charged particles (DecayCandidates.h):
///        through its composites, its products told from other daughters,
///        and what keeps it from being seen whole.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK decay candidates
#include "GenericTracking/Study/DecayCandidates.h"

#include <boost/test/unit_test.hpp>

#include <cstdint>
#include <cstdlib>
#include <vector>

namespace o2::gtrk::study
{

namespace
{

double chargeOf(std::int32_t pdg)
{
  switch (pdg) {
    case 211:
    case 2212:
    case -11:
      return 1.;
    case -211:
    case 3312:
    case 11:
      return -1.;
    default:
      return 0.;
  }
}

/// A Xi- (1) into a Lambda (2) and a pi- (3), the Lambda into p (4) and pi-
/// (5); a Sigma+ (6) into a neutron (7) and a pi+ (8); a photon (9) that
/// converted into e- (10) and e+ (11), with a Compton electron (12) before.
DecayTrees someTrees()
{
  return {
    .particles = {
      {.id = 1U, .pdg = 3312, .process = 0, .firstDaughter = 0, .daughterCount = 2},
      {.id = 2U, .pdg = 3122, .process = 4, .mother = 1U, .firstDaughter = 2, .daughterCount = 2},
      {.id = 3U, .pdg = -211, .process = 4, .mother = 1U},
      {.id = 4U, .pdg = 2212, .process = 4, .mother = 2U},
      {.id = 5U, .pdg = -211, .process = 4, .mother = 2U},
      {.id = 6U, .pdg = 3222, .process = 0, .firstDaughter = 4, .daughterCount = 2},
      {.id = 7U, .pdg = 2112, .process = 4, .mother = 6U},
      {.id = 8U, .pdg = 211, .process = 4, .mother = 6U},
      {.id = 9U, .pdg = 22, .process = 0, .firstDaughter = 6, .daughterCount = 3},
      {.id = 10U, .pdg = 11, .process = 5, .mother = 9U},
      {.id = 11U, .pdg = -11, .process = 5, .mother = 9U},
      {.id = 12U, .pdg = 11, .process = 6, .mother = 9U},
    },
    .daughters = {2U, 3U, 4U, 5U, 7U, 8U, 12U, 10U, 11U},
  };
}

std::vector<std::uint64_t> idsOf(const std::vector<const TreeParticle*>& particles)
{
  std::vector<std::uint64_t> ids;
  ids.reserve(particles.size());
  for (const auto* particle : particles) {
    ids.push_back(particle->id);
  }
  return ids;
}

} // namespace

BOOST_AUTO_TEST_CASE(ACascadeEndsInItsThreeChargedParticles)
{
  const auto trees = someTrees();
  const auto shape = decayShape(trees, *trees.find(1U), chargeOf);
  BOOST_TEST(shape.products == 2U);
  BOOST_TEST(idsOf(shape.leaves) == (std::vector<std::uint64_t>{4U, 5U, 3U}));
  BOOST_TEST(idsOf(shape.composites) == (std::vector<std::uint64_t>{2U}));
  BOOST_TEST(shape.depth == 1U);
  BOOST_TEST(shape.closed());
}

BOOST_AUTO_TEST_CASE(ANeutronLeavesTheDecayUnseen)
{
  const auto trees = someTrees();
  const auto shape = decayShape(trees, *trees.find(6U), chargeOf);
  BOOST_TEST(idsOf(shape.leaves) == (std::vector<std::uint64_t>{8U}));
  BOOST_TEST(shape.neutralUnseen == 1U);
  BOOST_TEST(!shape.closed());
}

BOOST_AUTO_TEST_CASE(AConversionsProductsAreItsPairOnly)
{
  const auto trees = someTrees();
  const auto shape = decayShape(trees, *trees.find(9U), chargeOf);
  BOOST_TEST(shape.products == 2U);
  BOOST_TEST(idsOf(shape.leaves) == (std::vector<std::uint64_t>{10U, 11U}));
  BOOST_TEST(shape.closed());
}

BOOST_AUTO_TEST_CASE(PartonsAndDiquarksAreNoTreeMembers)
{
  BOOST_TEST(isTreeSpecies(421));
  BOOST_TEST(isTreeSpecies(-4122));
  BOOST_TEST(isTreeSpecies(1010010030));
  BOOST_TEST(isTreeSpecies(22));
  BOOST_TEST(isTreeSpecies(-13));
  BOOST_TEST(!isTreeSpecies(4));
  BOOST_TEST(!isTreeSpecies(21));
  BOOST_TEST(!isTreeSpecies(2101));
  BOOST_TEST(!isTreeSpecies(92));
  BOOST_TEST(!isTreeSpecies(990));
  BOOST_TEST(!isTreeSpecies(9900110));
}

} // namespace o2::gtrk::study
