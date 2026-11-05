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

/// \file testDefinitions.cxx
/// \brief Exercises the shared index types, and stands as the build's canary.
///
/// The assertions here are deliberately modest -- the point of this test at
/// this stage is as much that it *links* as what it checks: it proves that the
/// GPUCommon annotations and Boost.Test resolve together against the O2
/// environment, which is the thing that actually breaks when the toolchain or
/// the aliBuild install moves.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Definitions
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Base/Definitions.h"

#include <cstdint>
#include <type_traits>

namespace o2::gtrk
{

BOOST_AUTO_TEST_CASE(IndexTypesAreTheDeclaredWidths)
{
  // These widths are part of the on-device record layout, so a silent change
  // here would change the size of every cluster and cell array.
  static_assert(std::is_same_v<SensorID, std::uint32_t>);
  static_assert(std::is_same_v<LayerID, std::uint16_t>);
  BOOST_TEST(sizeof(SensorID) == 4U);
  BOOST_TEST(sizeof(LayerID) == 2U);
}

BOOST_AUTO_TEST_CASE(SentinelsFailTheBoundsTestWithoutASpecialCase)
{
  // The whole reason the sentinels are ~0 rather than -1: one unsigned
  // comparison `id < n` has to reject both an out-of-range index and "none".
  static_assert(!(kInvalidSensorID < SensorID{14732}));
  static_assert(!(kInvalidLayerID < LayerID{46})); // 10 barrel + 36 discs
}

} // namespace o2::gtrk
