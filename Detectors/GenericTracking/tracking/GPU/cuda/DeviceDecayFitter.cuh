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

/// \file DeviceDecayFitter.cuh
/// \brief DeviceDecayFitter<Det, Material> on one backend: a batch's columns
///        uploaded, its decays' workspaces sized and the decays fitted, one
///        thread each, in launches that fit the arena.

#ifndef GTRK_TRACKING_GPU_DEVICEDECAYFITTER_CUH_
#define GTRK_TRACKING_GPU_DEVICEDECAYFITTER_CUH_

#include "DeviceAlgorithms.cuh"
#include "DeviceMaterial.cuh"
#include "Platform.cuh"

#include "GenericTracking/Tracking/DecayFitBatch.h"
#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceDecayFitter.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Tracking/GPU/DeviceStream.h"
#include "GenericTracking/Tracking/Profile.h"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE
{

namespace decay_fit_device
{

/// Decay \a index's workspace [bytes] (`batchDecayWorkspace`).
template <DetectorLayout Det, DecayFitVariant Variant, typename Material>
struct SizeWorkspace {
  DecayFitBatchView batch;
  const Material* material = nullptr;
  std::uint64_t* bytes = nullptr;
  GPUd() void operator()(std::uint32_t index) const { bytes[index] = batchDecayWorkspace<Det, Variant>(batch, index, *material); }
};

/// Decay first + i of \a batch in the workspace at offsets[i] of \a scratch,
/// up to offsets[i + 1], its result into results[i] (the launch's own).
template <DetectorLayout Det, DecayFitVariant Variant, typename Material>
struct FitDecay {
  DecayFitBatchView batch;
  const Material* material = nullptr;
  std::uint32_t first = 0U;
  const std::uint64_t* offsets = nullptr;
  std::byte* scratch = nullptr;
  DecayVertexFitResult* results = nullptr;
  GPUd() void operator()(std::uint32_t i) const
  {
    fitBatchDecay<Det, Variant>(batch, first + i, *material, std::span<std::byte>{scratch + offsets[i], static_cast<std::size_t>(offsets[i + 1U] - offsets[i])}, results[i]);
  }
};

/// The workspaces' alignment in a launch's scratch: every decay's starts at
/// a multiple of it (the fitter's own alignment divides it).
inline constexpr std::uint64_t kSliceAlignment = kArenaAlignment;

} // namespace decay_fit_device

template <DetectorLayout Det, typename Material>
class DeviceDecayFitterImpl final : public DeviceDecayFitter<Det, Material>
{
 public:
  using Timing = typename DeviceDecayFitter<Det, Material>::Timing;

  DeviceDecayFitterImpl(DeviceRuntime& runtime, DeviceBlock block) : mRuntime(&runtime), mArena(runtime, block), mStream(runtime.createStream()) {}
  ~DeviceDecayFitterImpl() final { mStream.synchronize(); }
  DeviceDecayFitterImpl(const DeviceDecayFitterImpl&) = delete;
  DeviceDecayFitterImpl& operator=(const DeviceDecayFitterImpl&) = delete;
  DeviceDecayFitterImpl(DeviceDecayFitterImpl&&) = delete;
  DeviceDecayFitterImpl& operator=(DeviceDecayFitterImpl&&) = delete;

  [[nodiscard]] bool configure(const Material& material) final
  {
    mStream.synchronize();
    mArena.reset();
    mMaterial.release();
    const bool uploaded = mMaterial.upload(mArena, material, mStream);
    mStream.synchronize(); // the host material may go once this returns
    return uploaded;
  }

  [[nodiscard]] bool fit(DecayFitVariant variant, const DecayFitBatch& batch, std::span<DecayVertexFitResult> results, Timing& timing) final
  {
    return withDecayFitVariant(variant, [&](auto one) { return fitAs<decltype(one)::value>(batch, results, timing); });
  }

 private:
  using Clock = std::chrono::steady_clock;

  template <DecayFitVariant Variant>
  [[nodiscard]] bool fitAs(const DecayFitBatch& batch, std::span<DecayVertexFitResult> results, Timing& timing)
  {
    using namespace decay_fit_device;
    timing = {};
    const auto count = batch.size();
    if (count == 0U) {
      return true;
    }
    mArena.resetPeak();
    // The inputs at the base end for this call; each launch's workspaces and
    // results at the scratch end, the results copied out when it is done.
    // The arena's peak: the inputs (decays, daughters, mothers, their rows,
    // the workspace sizes) plus the largest launch's offsets, its slices
    // (each workspace rounded up to kSliceAlignment) and its results.
    auto scope = mArena.scope(ArenaEnd::kBase);
    DeviceSpan<DecayFitInput> decays;
    DeviceSpan<DecayDaughter> daughters;
    DeviceSpan<DecayMotherRows> mothers;
    DeviceSpan<TrackingSurface> planes;
    DeviceSpan<Vec2D> positions;
    DeviceSpan<Vec3D> covariance;
    DeviceSpan<std::uint64_t> bytes;
    const auto allocated = decays.allocate(mArena, count) && daughters.allocate(mArena, batch.daughters().size()) && mothers.allocate(mArena, batch.mothers().size()) && planes.allocate(mArena, batch.planes().size()) && positions.allocate(mArena, batch.positions().size()) && covariance.allocate(mArena, batch.covariance().size()) && bytes.allocate(mArena, count);
    if (!allocated) {
      return false;
    }
    auto start = Clock::now();
    upload(batch.decays().data(), decays, count, mStream);
    upload(batch.daughters().data(), daughters, batch.daughters().size(), mStream);
    upload(batch.mothers().data(), mothers, batch.mothers().size(), mStream);
    upload(batch.planes().data(), planes, batch.planes().size(), mStream);
    upload(batch.positions().data(), positions, batch.positions().size(), mStream);
    upload(batch.covariance().data(), covariance, batch.covariance().size(), mStream);
    mStream.synchronize();
    timing.uploadSeconds = profile::secondsSince(start);

    const DecayFitBatchView view{.decays = decays.data(), .daughters = daughters.data(), .mothers = mothers.data(), .planes = planes.data(), .positions = positions.data(), .covariance = covariance.data(), .count = static_cast<std::uint32_t>(count)};
    start = Clock::now();
    algorithms::forEach(*mRuntime, mStream, static_cast<std::uint32_t>(count), SizeWorkspace<Det, Variant, Material>{.batch = view, .material = mMaterial.data(), .bytes = bytes.data()});
    std::vector<std::uint64_t> sizes;
    download(bytes, sizes, count, mStream);
    mStream.synchronize();
    timing.sizeSeconds = profile::secondsSince(start);

    // Launches of consecutive decays while their slices and results fit the
    // scratch the arena has left; a decay that alone does not is refused.
    std::vector<std::uint32_t> refusedDecays;
    std::vector<std::uint64_t> offsets;
    start = Clock::now();
    for (std::size_t first = 0; first < count;) {
      const std::size_t available = mArena.capacity() - mArena.used();
      offsets.assign(1U, 0U);
      std::size_t last = first;
      for (; last < count; ++last) {
        const std::uint64_t slice = ((sizes[last] + kSliceAlignment - 1U) / kSliceAlignment) * kSliceAlignment;
        const std::uint64_t total = offsets.back() + slice;
        const std::uint64_t offsetBytes = (offsets.size() + 1U) * sizeof(std::uint64_t) + kSliceAlignment;
        const std::uint64_t resultBytes = (offsets.size() * sizeof(DecayVertexFitResult)) + kSliceAlignment;
        if (total + offsetBytes + resultBytes + kSliceAlignment > available) {
          break;
        }
        offsets.push_back(total);
      }
      if (last == first) {
        refusedDecays.push_back(static_cast<std::uint32_t>(first));
        ++first;
        continue;
      }
      auto launch = mArena.scope(ArenaEnd::kScratch);
      const auto launched = static_cast<std::uint32_t>(last - first);
      DeviceSpan<std::uint64_t> deviceOffsets;
      DeviceSpan<std::byte> scratch;
      DeviceSpan<DecayVertexFitResult> fitted;
      if (!deviceOffsets.allocate(mArena, offsets.size(), ArenaEnd::kScratch) || !scratch.allocate(mArena, offsets.back() == 0U ? 1U : offsets.back(), ArenaEnd::kScratch) || !fitted.allocate(mArena, launched, ArenaEnd::kScratch)) {
        return false;
      }
      upload(offsets.data(), deviceOffsets, offsets.size(), mStream);
      algorithms::forEach(*mRuntime, mStream, launched, FitDecay<Det, Variant, Material>{.batch = view, .material = mMaterial.data(), .first = static_cast<std::uint32_t>(first), .offsets = deviceOffsets.data(), .scratch = scratch.data(), .results = fitted.data()});
      const auto copied = Clock::now();
      download(fitted, results.data() + first, launched, mStream);
      // The next launch takes these slices: this one must be done with them.
      mStream.synchronize();
      timing.downloadSeconds += profile::secondsSince(copied);
      ++timing.launches;
      first = last;
    }
    timing.fitSeconds = profile::secondsSince(start) - timing.downloadSeconds;
    for (const auto index : refusedDecays) {
      results[index].reset();
      results[index].status = DecayFitStatus::kTooLarge;
    }
    timing.refused = refusedDecays.size();
    timing.peakBytes = mArena.peak();
    return true;
  }

  DeviceRuntime* mRuntime = nullptr;
  DeviceArena mArena;
  DeviceStream mStream;
  DeviceMaterial<Material> mMaterial;
};

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE

#endif // GTRK_TRACKING_GPU_DEVICEDECAYFITTER_CUH_
