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

/// \file DeviceMaterial.cuh
/// \brief A material's copy on the device, for the device tracker and the
///        device decay fitter.

#ifndef GTRK_TRACKING_GPU_DEVICEMATERIAL_CUH_
#define GTRK_TRACKING_GPU_DEVICEMATERIAL_CUH_

#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceStorage.h"
#include "GenericTracking/Tracking/GPU/DeviceStream.h"

#include "FlatObject.h"

#include <type_traits>

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE
{

/// Where a material lives on the device: an O2 FlatObject -- the
/// MatLayerCylSet -- with its buffer, anything else as a value. Decided by
/// the base class, not by trivial copyability: O2 compiles a FlatObject's
/// copy constructor and destructor only outside GPU translation units, so in
/// this one a MatLayerCylSet looks trivially copyable, and copied as bytes it
/// keeps its host buffer pointers.
template <typename Material, bool Plain = !std::is_base_of_v<o2::gpu::FlatObject, Material>>
struct DeviceMaterial {
  static_assert(std::is_trivially_copyable_v<Material>, "a material that is not a FlatObject uploads as bytes");
  DeviceObject<Material> storage;
  [[nodiscard]] bool upload(DeviceArena& arena, const Material& material, DeviceStream& stream)
  {
    if (!storage.allocate(arena)) {
      return false;
    }
    uploadObject(material, storage, stream);
    return true;
  }
  void release() noexcept { storage = {}; }
  [[nodiscard]] const Material* data() const { return storage.data(); }
};

template <typename Material>
struct DeviceMaterial<Material, false> {
  DeviceFlatObject<Material> storage;
  [[nodiscard]] bool upload(DeviceArena& arena, const Material& material, DeviceStream& stream) { return storage.upload(arena, material, stream); }
  /// The device copy let go, the host image kept for the next run's upload.
  void release() noexcept { storage.release(); }
  [[nodiscard]] const Material* data() const { return storage.data(); }
};

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE

#endif // GTRK_TRACKING_GPU_DEVICEMATERIAL_CUH_
