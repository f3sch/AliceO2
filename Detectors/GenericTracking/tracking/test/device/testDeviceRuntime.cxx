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

// The device runtime -- properties, streams, the arena, the shared functions,
// the lane groups -- against the host; built once for CUDA and once for HIP.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Device Runtime
#include <boost/test/unit_test.hpp>

#include "DeviceTestFixture.h"
#include "DeviceTestKernels.h"
#include "PropagationPrecisionFixture.h"

#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Tracking/MeasurementGrid.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <memory_resource>
#include <numbers>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace o2::gtrk
{
namespace
{

using namespace device_test;

std::shared_ptr<BoundedMemoryResource> hostPool() { return std::make_shared<BoundedMemoryResource>(); }

template <template <template <class> class> class S>
bool sameColumns(const BoundedColumns<S>& left, const BoundedColumns<S>& right)
{
  if (left.size() != right.size()) {
    return false;
  }
  const auto lhs = left.view();
  const auto rhs = right.view();
  bool same = true;
  lhs.apply([&](const auto&... a) { rhs.apply([&](const auto&... b) { same = ((std::memcmp(a, b, left.size() * sizeof(*a)) == 0) && ...); }); });
  return same;
}

std::vector<float> smokeInputs()
{
  constexpr int kSteps = 2000;
  std::vector<float> inputs;
  inputs.reserve(kSteps + 3U);
  for (int i = 0; i < kSteps; ++i) {
    inputs.push_back(-4.F + (8.F * static_cast<float>(i) / static_cast<float>(kSteps - 1)));
  }
  inputs.push_back(std::numeric_limits<float>::quiet_NaN());
  inputs.push_back(1.e10F);
  inputs.push_back(-1.e10F);
  return inputs;
}

template <typename T>
double ulps(T a, T b)
{
  return static_cast<double>(std::abs(a - b)) / static_cast<double>(std::numeric_limits<T>::epsilon() * std::max(std::abs(a), T{1}));
}

} // namespace

BOOST_AUTO_TEST_CASE(PropertiesAndCompiledWaveWidthAgree)
{
  auto device = runtime();
  const auto& properties = device->properties();
  BOOST_TEST_MESSAGE(device->backend() << ": " << properties.name << ", " << properties.architecture << ", warp " << properties.warpSize << ", "
                                       << properties.computeUnits << " compute units");
  if (std::string{device->backend()} == "CUDA") {
    BOOST_TEST(properties.warpSize == 32);
  } else if (properties.architecture.starts_with("gfx9")) {
    BOOST_TEST(properties.warpSize == 64);
  }
}

BOOST_AUTO_TEST_CASE(StreamsAndCopies)
{
  auto device = runtime();
  auto stream = device->createStream();
  DeviceStream moved = std::move(stream);

  constexpr std::size_t kBytes = 4096U;
  std::vector<std::uint8_t> input(kBytes, std::uint8_t{0xA5});
  std::vector<std::uint8_t> output(kBytes, std::uint8_t{0});
  {
    DeviceArena arena{*device, block(kBytes)};
    DeviceSpan<std::uint8_t> storage;
    BOOST_REQUIRE(storage.allocate(arena, kBytes));
    upload(input.data(), storage, kBytes, moved);
    download(storage, output.data(), kBytes, moved);
    moved.synchronize();
    BOOST_TEST(std::memcmp(input.data(), output.data(), kBytes) == 0);
  }
}

BOOST_AUTO_TEST_CASE(ArenaAlignmentMetricsResetAndRefusal)
{
  auto device = runtime();
  DeviceArena arena{*device, block(2048U)};
  const auto first = arena.allocate(7U, 16U);
  const auto second = arena.allocate(19U, 256U);
  BOOST_REQUIRE(first);
  BOOST_REQUIRE(second);
  BOOST_TEST(static_cast<std::byte*>(second.data) - static_cast<std::byte*>(first.data) == 256);
  BOOST_TEST(arena.used() == 275U);
  BOOST_TEST(arena.peak() == 275U);

  const auto used = arena.used();
  const auto peak = arena.peak();
  const auto refused = arena.allocate(arena.capacity());
  BOOST_TEST(!static_cast<bool>(refused));
  BOOST_TEST(arena.used() == used);
  BOOST_TEST(arena.peak() == peak);
  arena.reset();
  const auto again = arena.allocate(7U, 16U);
  BOOST_TEST(again.data == first.data);
  arena.resetPeak();
  BOOST_TEST(arena.peak() == arena.used());
}

BOOST_AUTO_TEST_CASE(ScopesInvalidateOnlyWhatTheyRewind)
{
  auto device = runtime();
  DeviceArena arena{*device, block(4096U)};
  DeviceSpan<std::uint32_t> timeframe;
  BOOST_REQUIRE(timeframe.allocate(arena, 8U));
  DeviceSpan<std::uint32_t> slice;
  {
    auto scope = arena.scope();
    BOOST_REQUIRE(slice.allocate(arena, 8U));
    BOOST_CHECK_NO_THROW(static_cast<void>(slice.data()));
  }
  BOOST_CHECK_THROW(static_cast<void>(slice.data()), std::logic_error);
  BOOST_CHECK_NO_THROW(static_cast<void>(timeframe.data()));

  arena.reset();
  BOOST_CHECK_THROW(static_cast<void>(timeframe.data()), std::logic_error);
}

BOOST_AUTO_TEST_CASE(ColumnAllocationIsTransactional)
{
  auto device = runtime();
  DeviceArena arena{*device, block(800U)};
  const auto used = arena.used();
  DeviceColumns<TrackletT> columns;
  BOOST_TEST(!columns.allocate(arena, 1U));
  BOOST_TEST(arena.used() == used);
}

/// A transactional allocation opens its own scope and, on success, hands its
/// columns to the scope around it: they outlive an inner scope opened later and
/// die with the one they were allocated in.
BOOST_AUTO_TEST_CASE(CommittedColumnsLiveAndDieWithTheirScope)
{
  auto device = runtime();
  DeviceArena arena{*device, block(std::size_t{16} * 1024U)};
  DeviceColumns<TrackletT> columns;
  {
    auto slice = arena.scope();
    const auto before = arena.used();
    BOOST_REQUIRE(columns.allocate(arena, 4U));
    BOOST_TEST(arena.used() > before);
    BOOST_CHECK_NO_THROW(static_cast<void>(columns.view()));
    {
      auto inner = arena.scope();
      DeviceSpan<std::uint32_t> scratch;
      BOOST_REQUIRE(scratch.allocate(arena, 16U));
    }
    BOOST_CHECK_NO_THROW(static_cast<void>(columns.view()));
    BOOST_TEST(columns.size() == 4U);
  }
  BOOST_CHECK_THROW(static_cast<void>(columns.view()), std::logic_error);
}

/// The scratch end grows down from the top: both ends allocate without
/// touching each other, meet but never cross, and each end's scopes rewind
/// only that end.
BOOST_AUTO_TEST_CASE(TheScratchEndIsIndependentOfTheBase)
{
  auto device = runtime();
  constexpr std::size_t kCapacity = 4096U;
  DeviceArena arena{*device, block(kCapacity)};
  DeviceSpan<std::uint8_t> kept;
  BOOST_REQUIRE(kept.allocate(arena, 1000U));
  // The first allocation sits at the arena's base (the device allocation is 256-aligned).
  const auto* origin = reinterpret_cast<const std::byte*>(kept.data());
  BOOST_TEST(arena.used(ArenaEnd::kBase) == 1000U);
  BOOST_TEST(arena.used(ArenaEnd::kScratch) == 0U);

  DeviceSpan<std::byte> staging;
  {
    auto scratch = arena.scope(ArenaEnd::kScratch);
    BOOST_REQUIRE(staging.allocate(arena, 1000U, ArenaEnd::kScratch));
    const std::byte* at = staging.data();
    BOOST_TEST(at >= origin + arena.used(ArenaEnd::kBase));
    BOOST_TEST(at + 1000 <= origin + kCapacity);
    BOOST_TEST(std::bit_cast<std::uintptr_t>(at) % 256U == 0U);
    BOOST_TEST(arena.used() == arena.used(ArenaEnd::kBase) + arena.used(ArenaEnd::kScratch));

    // The ends meet: what does not fit between them is refused and changes nothing.
    const auto base = arena.used(ArenaEnd::kBase);
    const auto top = arena.used(ArenaEnd::kScratch);
    BOOST_TEST(!static_cast<bool>(arena.allocate(kCapacity - base - top + 1U, 1U, ArenaEnd::kScratch)));
    BOOST_TEST(!static_cast<bool>(arena.allocate(kCapacity - base - top + 1U, 1U, ArenaEnd::kBase)));
    BOOST_TEST(arena.used(ArenaEnd::kBase) == base);
    BOOST_TEST(arena.used(ArenaEnd::kScratch) == top);

    // A base-end scope opened and closed meanwhile leaves the scratch alone.
    {
      auto slice = arena.scope();
      BOOST_REQUIRE(arena.allocate(16U));
    }
    BOOST_TEST(arena.used(ArenaEnd::kScratch) == top);
    BOOST_CHECK_NO_THROW(static_cast<void>(staging.data()));
  }
  BOOST_TEST(arena.used(ArenaEnd::kScratch) == 0U);
  BOOST_CHECK_THROW(static_cast<void>(staging.data()), std::logic_error);
  BOOST_CHECK_NO_THROW(static_cast<void>(kept.data()));
  BOOST_TEST(arena.peak() >= 2000U);

  // Transactional columns on the scratch end, and a reset that clears both ends.
  DeviceColumns<TrackletT> columns;
  BOOST_REQUIRE(columns.allocate(arena, 8U, ArenaEnd::kScratch));
  BOOST_TEST(arena.used(ArenaEnd::kScratch) > 0U);
  BOOST_CHECK_NO_THROW(static_cast<void>(columns.view()));
  arena.reset();
  BOOST_TEST(arena.used() == 0U);
  BOOST_CHECK_THROW(static_cast<void>(columns.view()), std::logic_error);
}

BOOST_AUTO_TEST_CASE(SharedFunctionsStillAgreeWithTheHost)
{
  auto device = runtime();
  auto stream = device->createStream();
  auto inputs = smokeInputs();
  std::vector<SmokeRow> rows(inputs.size());
  DeviceArena arena{*device, block((inputs.size() * (sizeof(float) + sizeof(SmokeRow))) + 1024U)};
  DeviceSpan<float> deviceInputs;
  DeviceSpan<SmokeRow> deviceRows;
  DeviceSpan<std::uint32_t> maximum;
  BOOST_REQUIRE(deviceInputs.allocate(arena, inputs.size()));
  BOOST_REQUIRE(deviceRows.allocate(arena, inputs.size()));
  BOOST_REQUIRE(maximum.allocate(arena, 1U));
  upload(inputs.data(), deviceInputs, inputs.size(), stream);
  device->fill(maximum.data(), 0U, sizeof(std::uint32_t), stream);
  kernels().smoke(*device, stream, deviceInputs.data(), deviceRows.data(), maximum.data(), static_cast<std::uint32_t>(inputs.size()));
  download(deviceRows, rows.data(), rows.size(), stream);
  std::uint32_t maximumValue = 0U;
  download(maximum, &maximumValue, 1U, stream);
  stream.synchronize();

  std::uint32_t expectedMaximum = 0U;
  for (std::uint32_t i = 0; i < inputs.size(); ++i) {
    expectedMaximum = std::max(expectedMaximum, (i * 37U) % static_cast<std::uint32_t>(inputs.size()));
  }
  BOOST_TEST(maximumValue == expectedMaximum);

  double worstSinF = 0.;
  double worstAtan2F = 0.;
  double worstSinD = 0.;
  double worstAtan2D = 0.;
  double worstSqrtD = 0.;
  for (std::size_t i = 0; i < inputs.size(); ++i) {
    const float x = inputs[i];
    BOOST_TEST(rows[i].bin == clampedBin(x * 4.F, 16));
    if (!std::isfinite(x) || std::abs(x) > 10.F) {
      continue;
    }
    worstSinF = std::max(worstSinF, static_cast<double>(std::abs(rows[i].sinF - sinOf(x))));
    worstAtan2F = std::max(worstAtan2F, static_cast<double>(std::abs(rows[i].atan2F - atan2Of(x, 1.5F))));
    worstSinD = std::max(worstSinD, ulps(rows[i].sinD, sinOf(static_cast<double>(x))));
    worstAtan2D = std::max(worstAtan2D, ulps(rows[i].atan2D, atan2Of(static_cast<double>(x), 1.5)));
    worstSqrtD = std::max(worstSqrtD, ulps(rows[i].sqrtD, o2::math_utils::detail::sqrt(static_cast<double>(std::abs(x)))));
  }
#ifdef GPUCA_DETERMINISTIC_MODE
  BOOST_TEST(worstSinF == 0.);
  BOOST_TEST(worstAtan2F == 0.);
  BOOST_TEST(worstSinD == 0.);
  BOOST_TEST(worstAtan2D == 0.);
#else
  BOOST_TEST(worstSinF <= 1.e-5);
  BOOST_TEST(worstAtan2F <= 1.e-5);
  BOOST_TEST(worstSinD <= 3.);
  BOOST_TEST(worstAtan2D <= 3.);
#endif
  BOOST_TEST(worstSqrtD == 0.);
}

/// Every function of Trig.h, in both precisions, over the propagator's turns,
/// wide arguments, the reduction's hard cases and the special values: in the
/// deterministic build the host's bits, all of them (PortableMath.h).
BOOST_AUTO_TEST_CASE(TrigGivesTheHostsBitsInTheDeterministicBuild)
{
  std::vector<double> ys;
  std::vector<double> xs;
  std::mt19937_64 generator{42U};
  std::uniform_real_distribution<double> unit{-1., 1.};
  std::uniform_real_distribution<double> turnCosine{-0.2, 1.};
  std::uniform_real_distribution<double> wide{-60., 60.};
  for (int i = 0; i < 400'000; ++i) {
    ys.push_back(unit(generator));
    xs.push_back(turnCosine(generator));
  }
  for (int i = 0; i < 400'000; ++i) {
    ys.push_back(wide(generator));
    xs.push_back(wide(generator));
  }
  // Multiples of pi/2 and the next two doubles above each, where sin and cos
  // reduce with most cancellation; not above zero, where they are subnormal.
  for (int n = -200; n <= 200; ++n) {
    double y = n * (std::numbers::pi / 2.);
    for (int step = 0; step < (n == 0 ? 1 : 3); ++step) {
      ys.push_back(y);
      xs.push_back(-y);
      y = std::nextafter(y, std::numeric_limits<double>::infinity());
    }
  }
  // No subnormals: a process that loads O2's libO2GPUDataTypes.so runs with
  // the host's FTZ and DAZ set (it is linked with crtfastmath.o), and the
  // devices keep them, so a subnormal argument differs before any function
  // sees it. 1.5e6 is below 2^20 pi/2, where sin and cos stop being portable.
  constexpr double kInf = std::numeric_limits<double>::infinity();
  const std::array special{0., -0., 1., -1., 0.5, -0.5, 0.975, 1.e-30, 1.e30, 1.5e6, kInf, -kInf, std::numeric_limits<double>::quiet_NaN()};
  for (const double y : special) {
    for (const double x : special) {
      ys.push_back(y);
      xs.push_back(x);
    }
  }
  const auto count = static_cast<std::uint32_t>(ys.size());

  auto device = runtime();
  auto stream = device->createStream();
  DeviceArena arena{*device, block((std::size_t{count} * ((2U * sizeof(double)) + sizeof(TrigRow))) + 4096U)};
  DeviceSpan<double> deviceYs;
  DeviceSpan<double> deviceXs;
  DeviceSpan<TrigRow> deviceRows;
  BOOST_REQUIRE(deviceYs.allocate(arena, count));
  BOOST_REQUIRE(deviceXs.allocate(arena, count));
  BOOST_REQUIRE(deviceRows.allocate(arena, count));
  upload(ys.data(), deviceYs, count, stream);
  upload(xs.data(), deviceXs, count, stream);
  kernels().trig(*device, stream, deviceYs.data(), deviceXs.data(), deviceRows.data(), count);
  std::vector<TrigRow> rows(count);
  download(deviceRows, rows.data(), count, stream);
  stream.synchronize();

  constexpr std::size_t kFunctions = 15U;
  constexpr std::array<const char*, kFunctions> kNames{"sin", "cos", "sincos.sin", "sincos.cos", "asin", "atan", "atan2", "sinf",
                                                       "cosf", "sincosf.sin", "sincosf.cos", "asinf", "atanf", "atan2f", "logf"};
  std::array<std::size_t, kFunctions> differ{};
  const auto compare = [&](std::size_t function, auto host, auto onDevice) {
    using Bits = std::conditional_t<sizeof(host) == sizeof(std::uint64_t), std::uint64_t, std::uint32_t>;
    differ[function] += std::bit_cast<Bits>(host) == std::bit_cast<Bits>(onDevice) ? 0U : 1U;
  };
  for (std::uint32_t i = 0; i < count; ++i) {
    const TrigRow host = trigRow(ys[i], xs[i]);
    const TrigRow& onDevice = rows[i];
    compare(0U, host.sinD, onDevice.sinD);
    compare(1U, host.cosD, onDevice.cosD);
    compare(2U, host.sincosSinD, onDevice.sincosSinD);
    compare(3U, host.sincosCosD, onDevice.sincosCosD);
    compare(4U, host.asinD, onDevice.asinD);
    compare(5U, host.atanD, onDevice.atanD);
    compare(6U, host.atan2D, onDevice.atan2D);
    compare(7U, host.sinF, onDevice.sinF);
    compare(8U, host.cosF, onDevice.cosF);
    compare(9U, host.sincosSinF, onDevice.sincosSinF);
    compare(10U, host.sincosCosF, onDevice.sincosCosF);
    compare(11U, host.asinF, onDevice.asinF);
    compare(12U, host.atanF, onDevice.atanF);
    compare(13U, host.atan2F, onDevice.atan2F);
    compare(14U, host.logF, onDevice.logF);
  }
  for (std::size_t function = 0; function < kFunctions; ++function) {
    BOOST_TEST_MESSAGE(kNames[function] << ": " << differ[function] << " of " << count << " differ from the host");
#ifdef GPUCA_DETERMINISTIC_MODE
    BOOST_TEST(differ[function] == 0U, kNames[function] << " differs from the host in " << differ[function] << " of " << count);
#endif
  }
}

BOOST_AUTO_TEST_CASE(TrackletAndMeasurementViewsRoundTrip)
{
  auto device = runtime();
  auto stream = device->createStream();
  DeviceArena arena{*device, block(std::size_t{64} * 1024U)};

  const auto trackletMemory = hostPool();
  const auto trackletBackMemory = hostPool();
  BoundedColumns<TrackletT> tracklets;
  BoundedColumns<TrackletT> trackletsBack;
  tracklets.setMemoryResource(trackletMemory.get());
  trackletsBack.setMemoryResource(trackletBackMemory.get());
  tracklets.resize(3U);
  for (std::size_t i = 0; i < tracklets.size(); ++i) {
    tracklets.write(i, Tracklet{.source = 0xA1000000U + static_cast<std::uint32_t>(i),
                                .target = 0xB2000000U + static_cast<std::uint32_t>(i),
                                .link = static_cast<std::uint16_t>(0xC300U + i),
                                .lambda = std::bit_cast<float>(0x3ED40000U + static_cast<std::uint32_t>(i)),
                                .phi = std::bit_cast<float>(0xBEE50000U + static_cast<std::uint32_t>(i))});
  }
  DeviceColumns<TrackletT> deviceTracklets;
  DeviceColumns<TrackletT> deviceTrackletsBack;
  DeviceObject<TrackletConstView> trackletView;
  BOOST_REQUIRE(deviceTracklets.allocate(arena, tracklets.size()));
  BOOST_REQUIRE(deviceTrackletsBack.allocate(arena, tracklets.size()));
  BOOST_REQUIRE(trackletView.allocate(arena));
  upload<TrackletT>(std::as_const(tracklets).view(), tracklets.size(), deviceTracklets, stream);
  const TrackletConstView trackletPointers = std::as_const(deviceTracklets).view();
  uploadObject(trackletPointers, trackletView, stream);
  kernels().copyTracklets(*device, stream, trackletView.data(), deviceTrackletsBack.view(), static_cast<std::uint32_t>(tracklets.size()));
  download(deviceTrackletsBack, trackletsBack, tracklets.size(), stream);

  const auto measurementMemory = hostPool();
  const auto measurementBackMemory = hostPool();
  BoundedColumns<MeasurementT> measurements;
  BoundedColumns<MeasurementT> measurementsBack;
  measurements.setMemoryResource(measurementMemory.get());
  measurementsBack.setMemoryResource(measurementBackMemory.get());
  measurements.resize(3U);
  for (std::size_t i = 0; i < measurements.size(); ++i) {
    const auto f = static_cast<float>(i);
    measurements.write(i, Measurement{.loc0 = 1.F + f,
                                      .loc1 = 2.F + f,
                                      .timePS = 3'250U + static_cast<std::uint32_t>(i),
                                      .widthPS = 7U + static_cast<std::uint32_t>(i),
                                      .varLoc0 = 4.F + f,
                                      .covLoc01 = -5.F - f,
                                      .varLoc1 = 6.F + f,
                                      .clusterIndex = 0x88000000U + static_cast<std::uint32_t>(i),
                                      .surface = static_cast<SensorID>(0x9900U + i),
                                      .x = 10.F + f,
                                      .y = 11.F + f,
                                      .z = 12.F + f,
                                      .r = 13.F + f,
                                      .phi = 14.F + f,
                                      .varZ = 15.F + f,
                                      .varR = 16.F + f,
                                      .varRPhi = 17.F + f});
  }
  DeviceColumns<MeasurementT> deviceMeasurements;
  DeviceColumns<MeasurementT> deviceMeasurementsBack;
  DeviceObject<MeasurementView> measurementView;
  BOOST_REQUIRE(deviceMeasurements.allocate(arena, measurements.size()));
  BOOST_REQUIRE(deviceMeasurementsBack.allocate(arena, measurements.size()));
  BOOST_REQUIRE(measurementView.allocate(arena));
  upload<MeasurementT>(std::as_const(measurements).view(), measurements.size(), deviceMeasurements, stream);
  const MeasurementView measurementPointers = std::as_const(deviceMeasurements).view();
  uploadObject(measurementPointers, measurementView, stream);
  kernels().copyMeasurements(*device, stream, measurementView.data(), deviceMeasurementsBack.view(), static_cast<std::uint32_t>(measurements.size()));
  download(deviceMeasurementsBack, measurementsBack, measurements.size(), stream);
  stream.synchronize();

  BOOST_TEST(sameColumns(tracklets, trackletsBack));
  BOOST_TEST(sameColumns(measurements, measurementsBack));
}

// The lane groups: groups of 2 lanes up to a whole warp, 37 items --
// a last block not full, groups with no trip, groups of one warp on different
// trip counts and nested collectives -- each lane's ballots and broadcasts
// against what the script gives it on the host, every lane of the group
// evaluated.
BOOST_AUTO_TEST_CASE(LaneGroupsAgreeWhateverTheOtherGroupsOfTheWarpDo)
{
  using namespace lane_group_script;
  Device device;
  constexpr std::uint32_t kItems = 37U;
  for (std::uint32_t lanes = 2U; std::cmp_less_equal(lanes, device.runtime->properties().warpSize); lanes *= 2U) {
    std::vector<std::uint64_t> folded;
    BOOST_REQUIRE(kernels().laneGroups(*device.runtime, device.stream, *device.arena, lanes, kItems, folded) == StageStatus::kComplete);
    BOOST_REQUIRE(folded.size() == std::size_t{2U} * kItems * lanes);
    std::size_t mismatches = 0U;
    for (std::uint32_t item = 0U; item < kItems; ++item) {
      for (std::uint32_t lane = 0U; lane < lanes; ++lane) {
        std::array<std::uint64_t, 2> expected{};
        for (std::uint32_t trip = 0U; trip < trips(item); ++trip) {
          std::uint64_t ballot = 0U;
          for (std::uint32_t other = 0U; other < lanes; ++other) {
            ballot |= vote(other, trip, item) ? std::uint64_t{1} << other : 0U;
          }
          expected[0] ^= ballot * (trip + 1U);
          const auto got = payload(source(trip, item, lanes), trip, item);
          expected[1] += got.a + (7U * got.b) + (13U * got.c);
          for (std::uint32_t inner = 0U; inner < innerTrips(trip, item); ++inner) {
            std::uint64_t innerBallot = 0U;
            for (std::uint32_t other = 0U; other < lanes; ++other) {
              innerBallot |= innerVote(other, inner) ? std::uint64_t{1} << other : 0U;
            }
            expected[0] ^= innerBallot << (inner + 1U);
            expected[1] += innerValue(innerSource(inner, trip, lanes), inner);
          }
        }
        const std::size_t at = 2U * ((std::size_t{item} * lanes) + lane);
        mismatches += folded[at] != expected[0] || folded[at + 1U] != expected[1] ? 1U : 0U;
      }
    }
    BOOST_TEST_CONTEXT("groups of " << lanes << " lanes") { BOOST_TEST(mismatches == 0U); }
  }
}

// Independent double helix evaluation at each device-reported flight length,
// over the same geometry distribution as the CPU precision sweep. No plane
// inversion is used by the oracle; both charges/fields, straight trajectories,
// grazing planes, off-centre targets and forward/backward steps are covered.
BOOST_AUTO_TEST_CASE(DevicePropagationClosesWithinOneMicron)
{
  Device device;
  using namespace propagation_precision;
  std::mt19937_64 random{20260928U};
  std::uniform_real_distribution<double> unit{0., 1.};
  const auto sample = [&](double low, double high) { return low + (high - low) * unit(random); };
  std::vector<PropagationProbe> probes;
  constexpr int kTrials = 1000000;
  probes.reserve(kTrials);
  for (int trial = 0; trial < kTrials; ++trial) {
    const double phi = sample(-std::numbers::pi, std::numbers::pi);
    const double theta = sample(0.01, std::numbers::pi - 0.01);
    const double pt = std::exp(sample(std::log(0.03), std::log(10000.)));
    const double charge = trial % 2 == 0 ? 1. : -1.;
    const double radius = sample(0., 90.);
    const double azimuth = sample(-std::numbers::pi, std::numbers::pi);
    PropagationProbe probe;
    probe.from =
      plane(radius * std::cos(azimuth), radius * std::sin(azimuth), sample(-300., 300.), sample(-std::numbers::pi, std::numbers::pi), trial % 3 == 0, 1);
    auto& state = probe.before;
    state.surface = probe.from.id;
    state.par[kLoc0] = static_cast<Scalar>(sample(-1.5, 1.5));
    state.par[kLoc1] = static_cast<Scalar>(sample(-3.5, 3.5));
    state.par[kPhi] = static_cast<Scalar>(phi);
    state.par[kTheta] = static_cast<Scalar>(theta);
    state.par[kQOverP] = static_cast<Scalar>(charge * std::sin(theta) / pt);
    probe.prop.bzKiloGauss = trial % 5 == 0 ? -20.F : 20.F;
    if (trial % 17 == 0) {
      probe.prop.bzKiloGauss = 0.F;
    }
    const double path = sample(0.001, 300.) * (trial % 7 == 0 ? -1. : 1.);
    probe.direction = path > 0. ? PropagationDirection::kForward : PropagationDirection::kBackward;
    const double turn = static_cast<double>(state.qOverP()) * static_cast<double>(probe.prop.bzKiloGauss) * static_cast<double>(constants::math::B2C) * path;
    if (std::abs(turn) > 3.1) {
      continue;
    }
    ExactScalar origin[3]{};
    detail::boundToGlobal(probe.from, state.loc0(), state.loc1(), origin);
    const double transverse = path * std::sin(static_cast<double>(state.theta()));
    const double along = transverse * sinc(turn);
    const double across = transverse * versine(turn);
    const double endX = origin[0] + along * std::cos(static_cast<double>(state.phi())) - across * std::sin(static_cast<double>(state.phi()));
    const double endY = origin[1] + along * std::sin(static_cast<double>(state.phi())) + across * std::cos(static_cast<double>(state.phi()));
    const double endZ = origin[2] + path * std::cos(static_cast<double>(state.theta()));
    probe.to = plane(endX, endY, endZ, sample(-std::numbers::pi, std::numbers::pi), trial % 3 == 1, 2);
    const double offsets[2]{sample(-5., 5.), sample(-5., 5.)};
    for (int axis = 0; axis < 3; ++axis) {
      probe.to.origin[axis] +=
        static_cast<Scalar>(offsets[0] * static_cast<double>(probe.to.axis0[axis]) + offsets[1] * static_cast<double>(probe.to.axis1[axis]));
    }
    probes.push_back(probe);
  }
  // Match the CPU's million-trial coverage without enlarging the device
  // arena: each callback closes its temporary scope before the next batch.
  constexpr std::size_t kBatch = 100000U;
  std::vector<PropagationAnswer> answers;
  Worst worst;
  for (std::size_t base = 0U; base < probes.size(); base += kBatch) {
    const auto count = std::min(kBatch, probes.size() - base);
    const auto batch = std::span<const PropagationProbe>{probes}.subspan(base, count);
    BOOST_REQUIRE(kernels().propagation(*device.runtime, device.stream, *device.arena, batch, answers) == StageStatus::kComplete);
    BOOST_REQUIRE(answers.size() == count);
    for (std::size_t row = 0U; row < count; ++row) {
      const auto& probe = batch[row];
      const auto& answer = answers[row];
      if (answer.result.status == PropagationStatus::kOk || answer.result.status == PropagationStatus::kOutOfBounds) {
        worst.add(closure(probe.prop, probe.before, answer.after, probe.from, probe.to, answer.result.pathLength));
      }
    }
  }
  worst.report(device.runtime->backend());
  BOOST_TEST(worst.checked > kTrials / 2);
  BOOST_TEST(worst.value.position < kMicron);
  BOOST_TEST(worst.value.plane < kMicron);
  BOOST_TEST(worst.value.azimuth < 1.e-6);
}

} // namespace o2::gtrk
