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

/// \file testDeviceMaterial.cxx
/// \brief O2's material LUT ported to the device with DeviceFlatObject, its
///        lookups against the host's.
///
/// Reads the LUT GTRK_TEST_MATLUT names and passes without
/// checking anything if there is none. Built once per backend. Needs the
/// device.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK device material
#include <boost/test/unit_test.hpp>

#include "DeviceTestFixture.h"
#include "DeviceTestKernels.h"

#include "DetectorsBase/MatLayerCylSet.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace o2::gtrk
{

using namespace device_test;

BOOST_AUTO_TEST_CASE(TheLutGivesTheHostsBudgetsOnTheDevice)
{
  const char* path = std::getenv("GTRK_TEST_MATLUT");
  if (path == nullptr || !std::filesystem::exists(path)) {
    BOOST_TEST_MESSAGE("no material LUT: nothing to compare");
    return;
  }
  const std::unique_ptr<o2::base::MatLayerCylSet> lut{o2::base::MatLayerCylSet::loadFromFile(path)};
  BOOST_REQUIRE(lut != nullptr);
  std::mt19937 random{13U};
  std::uniform_real_distribution<float> transverse{-60.F, 60.F};
  std::uniform_real_distribution<float> longitudinal{-150.F, 150.F};
  std::uniform_real_distribution<float> step{-3.F, 3.F};
  std::vector<std::array<float, 6>> segments(20000U);
  for (auto& segment : segments) {
    const float x = transverse(random);
    const float y = transverse(random);
    const float z = longitudinal(random);
    segment = {x, y, z, x + step(random), y + step(random), z + step(random)};
  }
  // The arena holds the LUT's flat buffer and the segments.
  Device device{lut->getFlatBufferSize() + kArenaBytes};
  std::vector<std::array<float, 5>> onDevice;
  BOOST_REQUIRE((kernels().matBudgets(*device.runtime, device.stream, *device.arena, *lut, segments, onDevice) == StageStatus::kComplete));
  // Per quantity (mean density, x/X0, length, the moment's x/X0 and its
  // second moment): the host's values, and the largest difference.
  std::array<std::vector<float>, 5> hosts;
  std::size_t differing = 0U;
  std::size_t nonZero = 0U;
  for (std::size_t i = 0; i < segments.size(); ++i) {
    const auto& s = segments[i];
    const auto budget = lut->getMatBudget(s[0], s[1], s[2], s[3], s[4], s[5]);
    const auto moment = lut->getMatMoment(s[0], s[1], s[2], s[3], s[4], s[5]);
    const std::array<float, 5> host{budget.meanRho, budget.meanX2X0, budget.length, moment.x2X0, moment.x2X0T2};
    nonZero += budget.meanX2X0 > 0.F ? 1U : 0U;
    bool same = true;
    for (std::size_t q = 0; q < host.size(); ++q) {
      hosts[q].push_back(host[q]);
      same = same && std::bit_cast<std::uint32_t>(host[q]) == std::bit_cast<std::uint32_t>(onDevice[i][q]);
    }
    differing += same ? 0U : 1U;
  }
  // A difference against the quantity's mean over the segments: fast math
  // moves where a segment crosses a cell boundary, which on a sliver of a
  // cell is all of a budget near zero.
  double worst = 0.;
  for (std::size_t q = 0; q < hosts.size(); ++q) {
    double mean = 0.;
    for (const float v : hosts[q]) {
      mean += std::abs(static_cast<double>(v));
    }
    mean /= static_cast<double>(hosts[q].size());
    for (std::size_t i = 0; i < segments.size(); ++i) {
      const double difference = std::abs(static_cast<double>(onDevice[i][q]) - static_cast<double>(hosts[q][i]));
      worst = std::max(worst, mean > 0. ? difference / mean : difference);
    }
  }
  BOOST_TEST_MESSAGE(runtime()->backend() << ": " << segments.size() << " segments, " << nonZero << " through material, " << differing << " differing, worst difference over the mean " << worst);
  BOOST_TEST(nonZero > 0U);
#ifdef GPUCA_DETERMINISTIC_MODE
  BOOST_TEST(differing == 0U);
#else
  // A port that is wrong is off by far more.
  BOOST_TEST(worst < 1.e-2);
#endif
}

} // namespace o2::gtrk
