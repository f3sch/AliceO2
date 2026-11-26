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

/// \file testTrackletSort.cxx
/// \brief Each source's tracklets in (link, lambda, target) order.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK TrackletSort
#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/TrackletSort.h"

#include "SerialReferences.h"
#include "StageOutputComparison.h"

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <vector>

namespace o2::gtrk
{

namespace
{

using Output = StageOutput<TrackletT>;
using Cpu = CpuStage<TrackletT>;
using Index = Output::Index;

constexpr Index kElements = 97;

// Repeated values, and -0 next to +0, so that lambda ties are common and the
// target has to decide them.
constexpr std::array<Scalar, 7> kLambdas{0.25F, -0.F, 0.F, 0.25F, -1.5F, 1.25F, 0.F};

/// Element e emits (e * 7) % 11 tracklets, zero and one among them, with
/// unique targets in scrambled order and keys that tie on link and lambda.
struct Generator {
  template <typename Emitter>
  void operator()(Index element, Emitter& emitter) const
  {
    const Index count = (element * 7U) % 11U;
    for (Index k = 0; k < count; ++k) {
      // k -> 7k + e is one-to-one modulo the prime 11, and k < count <= 10.
      const Index scrambled = ((k * 7U) + element) % 11U;
      emitter.emit(Tracklet{
        .source = element,
        .target = (element * 1000U) + scrambled,
        .link = static_cast<std::uint16_t>((k * 3U + element) % 3U),
        .lambda = kLambdas[(k + element) % kLambdas.size()],
        .phi = static_cast<Scalar>(k) - static_cast<Scalar>(element),
      });
    }
  }
};

std::vector<Tracklet> segment(const Output& output, Index element)
{
  const auto view = output.view();
  std::vector<Tracklet> result;
  for (auto i = output.offsets()[element]; i < output.offsets()[element + 1U]; ++i) {
    result.push_back({.source = view.source[i], .target = view.target[i], .link = view.link[i], .lambda = view.lambda[i], .phi = view.phi[i]});
  }
  return result;
}

void generateSerial(Output& output)
{
  output.setMemoryPool(freshPool());
  output.prepare(kElements);
  BOOST_REQUIRE((runSerial(output, Generator{}) == StageStatus::kComplete));
}

} // namespace

BOOST_AUTO_TEST_CASE(EverySegmentIsInKeyOrderAndKeepsItsTracklets)
{
  Output unsorted;
  generateSerial(unsorted);
  Output sorted;
  generateSerial(sorted);
  sortTrackletSegments(sorted);
  BOOST_REQUIRE(sorted.size() == unsorted.size());

  bool anyTie = false;
  for (Index element = 0; element < kElements; ++element) {
    const auto before = segment(unsorted, element);
    const auto after = segment(sorted, element);
    BOOST_REQUIRE(after.size() == before.size());
    for (std::size_t i = 1; i < after.size(); ++i) {
      BOOST_TEST(trackletBefore(after[i - 1U], after[i]), "element " << element << " position " << i);
      anyTie = anyTie || (after[i - 1U].link == after[i].link && !(after[i - 1U].lambda < after[i].lambda));
    }
    // The same tracklets, in the order std::sort gives with the same key.
    auto expected = before;
    std::ranges::sort(expected, [](const Tracklet& a, const Tracklet& b) { return trackletBefore(a, b); });
    BOOST_TEST(std::ranges::equal(after, expected, sameTracklet), "element " << element);
  }
  BOOST_TEST(anyTie); // the target decided at least one equal link and lambda
}

namespace
{

constexpr Index kVariedElements = 40;

/// Element e emits a source of a length from 0 to 2000 with one of three
/// lambda patterns -- one lambda for the whole source, a few repeated values
/// with -0 and +0, or all distinct -- over 12 links, its targets unique and
/// scrambled, and the same targets in every source, so that neighbouring
/// sources share keys.
struct VariedGenerator {
  static constexpr std::array<Index, 10> kLengths{0U, 1U, 2U, 3U, 31U, 32U, 33U, 129U, 512U, 2000U};

  template <typename Emitter>
  void operator()(Index element, Emitter& emitter) const
  {
    const Index count = kLengths[element % kLengths.size()];
    const auto pattern = (element / kLengths.size()) % 3U;
    for (Index k = 0; k < count; ++k) {
      Scalar lambda = 0.5F;
      if (pattern == 1U) {
        lambda = kLambdas[(std::size_t{k} * 5U) % kLambdas.size()];
      } else if (pattern == 2U) {
        lambda = static_cast<Scalar>(static_cast<std::int32_t>((k * 7919U) % 4001U) - 2000) * 1.e-3F;
      }
      emitter.emit(Tracklet{
        .source = element,
        .target = (k * 7919U) % 100003U,
        .link = static_cast<std::uint16_t>((k * 5U + element) % 12U),
        .lambda = lambda,
        .phi = static_cast<Scalar>(k),
      });
    }
  }
};

/// One source as long as a source can be, over 12 links, every lambda the same.
struct LongGenerator {
  static constexpr auto kCount = static_cast<Index>(Output::kMaxPerElement);

  template <typename Emitter>
  void operator()(Index /*element*/, Emitter& emitter) const
  {
    for (Index k = 0; k < kCount; ++k) {
      emitter.emit(Tracklet{.source = 0U, .target = (k * 7919U) % 100003U, .link = static_cast<std::uint16_t>(k % 12U), .lambda = 0.F, .phi = 0.F});
    }
  }
};

} // namespace

BOOST_AUTO_TEST_CASE(TheHostSortGivesTheReferenceOrderForAnyLengthAndTies)
{
  Output serial;
  serial.setMemoryPool(freshPool());
  serial.prepare(kVariedElements);
  BOOST_REQUIRE((runSerial(serial, VariedGenerator{}) == StageStatus::kComplete));
  sortTrackletSegments(serial);

  for (const std::int32_t threads : {1, 2, 4}) {
    for (const Index chunk : {1U, 3U, 64U}) {
      Output parallel;
      const auto parallelPool = freshPool();
      parallel.setMemoryPool(parallelPool);
      parallel.prepare(kVariedElements);
      detail::TaskArena cpuArena{threads};
      Cpu cpu{cpuArena};
      cpu.setMemoryPool(parallelPool);
      BOOST_REQUIRE((cpu.run(parallel, chunk, VariedGenerator{}) == StageStatus::kComplete));
      const auto view = parallel.view();
      const auto offsets = parallel.offsets();
      const auto sort = [&] { cpu.forEachElement(parallel.elementCount(), chunk, [&](Index element) { sortTrackletSegmentOnHost(view, offsets[element], offsets[element + 1U]); }); };
      sort();
      BOOST_TEST(sameStageOutput(parallel, serial), "threads " << threads << ", chunk " << chunk);
      // Sorting a sorted output changes nothing.
      sort();
      BOOST_TEST(sameStageOutput(parallel, serial), "threads " << threads << ", chunk " << chunk << ", twice");
    }
  }
}

BOOST_AUTO_TEST_CASE(TheHostSortTakesOneLongSourceOfOneLambda)
{
  // 65536 tracklets of one source, every lambda equal: the target alone
  // decides within a link. The serial reference would take ~10^9 moves here,
  // so the expectation is std::sort's on the same key.
  Output output;
  output.setMemoryPool(freshPool());
  output.prepare(1U);
  BOOST_REQUIRE((runSerial(output, LongGenerator{}) == StageStatus::kComplete));
  auto expected = segment(output, 0U);
  std::ranges::sort(expected, [](const Tracklet& a, const Tracklet& b) { return trackletBefore(a, b); });

  sortTrackletSegmentOnHost(output.view(), output.offsets()[0], output.offsets()[1]);
  BOOST_TEST(std::ranges::equal(segment(output, 0U), expected, sameTracklet));
}

BOOST_AUTO_TEST_CASE(TheKeysOrderLambdasWithDenormalsAsZero)
{
  // Every pair of lambdas -- both zeros, the denormals, both signs, the
  // range's ends, ties -- and links: the keys order the lambdas as their
  // values do once the denormals are zero, whether or not this process
  // flushes them, and the link-lambda keys as `trackletBefore` with the
  // target aside.
  const auto flushed = [](Scalar x) { return std::abs(x) < FLT_MIN ? 0.F : x; };
  const std::array<Scalar, 13> lambdas{-1.5707964F, -1.F, -0.25F, -1.e-30F, -1.e-45F, -0.F, 0.F, 1.e-45F, 1.e-30F, 0.25F, 0.5F, 1.F, 1.5707964F};
  const std::array<std::uint16_t, 3> links{0U, 1U, 11U};
  for (const Scalar a : lambdas) {
    for (const Scalar b : lambdas) {
      BOOST_TEST((trackletLambdaKey(a) < trackletLambdaKey(b)) == (flushed(a) < flushed(b)), a << " " << b);
      BOOST_TEST((trackletLambdaKey(a) == trackletLambdaKey(b)) == exactlyEqual(flushed(a), flushed(b)), a << " " << b);
      for (const std::uint16_t la : links) {
        for (const std::uint16_t lb : links) {
          const bool keyBefore = trackletLinkLambdaKey(la, a) < trackletLinkLambdaKey(lb, b);
          // The same target on both sides: only link and lambda decide.
          BOOST_TEST(keyBefore == trackletBefore(la, a, 0U, lb, b, 0U), la << " " << a << " " << lb << " " << b);
        }
      }
    }
  }
}

} // namespace o2::gtrk
