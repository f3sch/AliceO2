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

/// \file DeviceRuntime.h
/// \brief A GPU backend's runtime: streams, copies, fills and the device
///        trackers, vendor-neutral.
///
/// Made for the device O2's GPU reconstruction initialised, by the
/// reconstruction's backend (GPUReconstruction::GetGenericTrackingRuntime);
/// the device memory is the reconstruction's pool, of which each tracker gets
/// a block.

#ifndef GTRK_TRACKING_GPU_DEVICERUNTIME_H_
#define GTRK_TRACKING_GPU_DEVICERUNTIME_H_

#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/ITS3.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Detectors/TRKFT3TOF.h"
#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceStream.h"
#include "GenericTracking/Tracking/LinkMaterial.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace o2::gtrk
{

/// Declared only: the runtime names the trackers it makes, their users
/// include DeviceTracker.h.
template <DetectorLayout Det>
class DeviceTracker;
template <DetectorLayout Det, typename Material>
class DeviceDecayFitter;

/// One device-to-device copy of copyEachOnDevice.
struct DeviceCopy {
  void* to = nullptr;
  const void* from = nullptr;
  std::size_t bytes = 0U;
};

struct DeviceProperties {
  std::string name;
  std::string architecture;
  int warpSize = 0;
  int residentThreads = 0; ///< the multiprocessors' threads at once, all of them
  int multiprocessors = 0;
  std::size_t freeMemoryBytes = 0; ///< the device's free memory when the runtime was made
};

/// Every call runs on the device current on the calling thread: the one
/// the reconstruction initialised.
class DeviceRuntime
{
 public:
  DeviceRuntime() = default;
  virtual ~DeviceRuntime() = default;
  DeviceRuntime(const DeviceRuntime&) = delete;
  DeviceRuntime& operator=(const DeviceRuntime&) = delete;
  DeviceRuntime(DeviceRuntime&&) = delete;
  DeviceRuntime& operator=(DeviceRuntime&&) = delete;

  [[nodiscard]] virtual const char* backend() const noexcept = 0;
  [[nodiscard]] virtual const DeviceProperties& properties() const noexcept = 0;

  [[nodiscard]] DeviceStream createStream();

  virtual void copyToDevice(void* to, const void* from, std::size_t bytes, DeviceStream& stream) = 0;
  virtual void copyToHost(void* to, const void* from, std::size_t bytes, DeviceStream& stream) = 0;
  virtual void copyOnDevice(void* to, const void* from, std::size_t bytes, DeviceStream& stream) = 0;
  /// Every copy of \a copies, none overlapping another's destination or
  /// source, by one kernel (a few for many copies) instead of a copy each.
  virtual void copyEachOnDevice(std::span<const DeviceCopy> copies, DeviceStream& stream) = 0;
  virtual void fill(void* to, std::uint8_t value, std::size_t bytes, DeviceStream& stream) = 0;

  /// Backend launchers use this to recover cudaStream_t/hipStream_t without
  /// exposing either type in a public header.
  [[nodiscard]] static void* nativeStream(DeviceStream& stream);

  /// The device tracker for \a D, its arena on \a block; std::runtime_error
  /// for a layout the backends are not built for.
  template <DetectorLayout D>
  [[nodiscard]] std::unique_ptr<DeviceTracker<D>> makeTracker(DeviceBlock block)
  {
    if constexpr (std::is_same_v<D, detectors::ITS>) {
      return makeITSTracker(block);
    } else if constexpr (std::is_same_v<D, detectors::ITS3>) {
      return makeITS3Tracker(block);
    } else if constexpr (std::is_same_v<D, detectors::MFT>) {
      return makeMFTTracker(block);
    } else if constexpr (std::is_same_v<D, detectors::TRKFT3>) {
      return makeTRKFT3Tracker(block);
    } else if constexpr (std::is_same_v<D, detectors::TRKFT3TOF>) {
      return makeTRKFT3TOFTracker(block);
    } else {
      throw std::runtime_error(std::string{"GPU backend "} + backend() + ": no device tracker is built for " + D::name);
    }
  }

  /// The device decay vertex fitter for \a D (DeviceDecayFitter.h), its
  /// arena on \a block; std::runtime_error for a layout the backends do not
  /// fit decays of (ALICE 3's TRKFT3 only).
  template <DetectorLayout D>
  [[nodiscard]] std::unique_ptr<DeviceDecayFitter<D, LinkMaterial<D>>> makeDecayFitter(DeviceBlock block)
  {
    if constexpr (std::is_same_v<D, detectors::TRKFT3>) {
      return makeTRKFT3DecayFitter(block);
    } else {
      throw std::runtime_error(std::string{"GPU backend "} + backend() + ": no device decay fitter is built for " + D::name);
    }
  }

 protected:
  [[nodiscard]] virtual void* createStreamHandle() = 0;

 private:
  friend class DeviceStream;

  virtual void destroyStreamHandle(void* handle) noexcept = 0;
  virtual void synchronizeStreamHandle(void* handle) = 0;

  [[nodiscard]] virtual std::unique_ptr<DeviceTracker<detectors::ITS>> makeITSTracker(DeviceBlock block) = 0;
  [[nodiscard]] virtual std::unique_ptr<DeviceTracker<detectors::ITS3>> makeITS3Tracker(DeviceBlock block) = 0;
  [[nodiscard]] virtual std::unique_ptr<DeviceTracker<detectors::MFT>> makeMFTTracker(DeviceBlock block) = 0;
  [[nodiscard]] virtual std::unique_ptr<DeviceTracker<detectors::TRKFT3>> makeTRKFT3Tracker(DeviceBlock block) = 0;
  [[nodiscard]] virtual std::unique_ptr<DeviceTracker<detectors::TRKFT3TOF>> makeTRKFT3TOFTracker(DeviceBlock block) = 0;
  [[nodiscard]] virtual std::unique_ptr<DeviceDecayFitter<detectors::TRKFT3, LinkMaterial<detectors::TRKFT3>>> makeTRKFT3DecayFitter(DeviceBlock block) = 0;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GPU_DEVICERUNTIME_H_
