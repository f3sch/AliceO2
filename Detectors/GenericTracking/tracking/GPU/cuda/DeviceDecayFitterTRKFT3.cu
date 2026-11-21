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

/// \file DeviceDecayFitterTRKFT3.cu
/// \brief The backend's device decay vertex fitter for ALICE 3's TRKFT3.

#include "DeviceDecayFitter.cuh"

#include "GenericTracking/Detectors/TRKFT3.h"

#include "DetectorsBase/MatLayerCylSet.h"

#include <memory>

namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE
{

std::unique_ptr<DeviceDecayFitter<detectors::TRKFT3, LinkMaterial<detectors::TRKFT3>>> makeTRKFT3DecayFitter(DeviceRuntime& runtime, DeviceBlock block) { return std::make_unique<DeviceDecayFitterImpl<detectors::TRKFT3, LinkMaterial<detectors::TRKFT3>>>(runtime, block); }

} // namespace o2::gtrk::gpu::GTRK_BACKEND_NAMESPACE
