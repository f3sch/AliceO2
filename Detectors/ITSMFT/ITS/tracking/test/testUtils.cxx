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

#include <boost/test/tools/old/interface.hpp>
#define BOOST_TEST_MODULE ITS Utils
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "ITStracking/Utils.h"
using namespace o2::its::utils;

BOOST_AUTO_TEST_CASE(test_1)
{
  const std::array<Bracket, 2> brackets{{{594, 1188},
                                         {297, 891}}};
  const auto result = computeSmallestBracket(brackets);
  BOOST_CHECK_EQUAL(result.first, 594);
  BOOST_CHECK_EQUAL(result.second, 891);
}

BOOST_AUTO_TEST_CASE(test_2)
{
  const std::array<Bracket, 7> brackets{{{594, 1188},
                                         {594, 1188},
                                         {594, 1188},
                                         {297, 891},
                                         {297, 891},
                                         {297, 891},
                                         {891, 1485}}};
  const auto result = computeSmallestBracket(brackets);
  BOOST_CHECK_EQUAL(result.first, 594);
  BOOST_CHECK_EQUAL(result.second, 1188);
}

BOOST_AUTO_TEST_CASE(test_3)
{
  std::array<Bracket, 7> brackets;
  brackets.fill(InvalidBracket);
  brackets[0] = {500, 600};
  brackets[3] = {500, 600};
  brackets[6] = {500, 600};
  const auto result = computeSmallestBracket(brackets);
  BOOST_CHECK_EQUAL(result.first, 500);
  BOOST_CHECK_EQUAL(result.second, 600);
}

BOOST_AUTO_TEST_CASE(test_4)
{
  std::array<Bracket, 3> brackets;
  brackets.fill(InvalidBracket);
  brackets[0] = {100, 200};
  brackets[2] = {500, 600};
  const auto result = computeSmallestBracket(brackets);
  BOOST_CHECK(result == InvalidBracket);
}

BOOST_AUTO_TEST_CASE(test_5)
{
  std::array<Bracket, 7> brackets{{{5346, 5940},
                                   {5346, 5940},
                                   {5346, 5940},
                                   {5049, 5643},
                                   {4455, 5049},
                                   {5049, 5643},
                                   {5049, 5643}}};
  const auto result = computeSmallestBracket(brackets);
  BOOST_CHECK(result != InvalidBracket);
}

BOOST_AUTO_TEST_CASE(test_6)
{
  std::array<Bracket, 7> brackets{{{46926, 47520},
                                   {46926, 47520},
                                   {46926, 47520},
                                   {47817, 48411},
                                   {47817, 48411},
                                   {47817, 48411},
                                   {47817, 48411}}};
  const auto result = computeSmallestBracket(brackets);
  BOOST_CHECK(result == InvalidBracket);
}
