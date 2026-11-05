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

/// \file DeviceCompileCheckToyDisc3.cu
/// \brief DeviceCompileCheck.cuh for ToyDisc3.

// First, as in O2's own .cu files: hipify-perl turns it into
// hip/hip_runtime.h, which O2's GPU headers expect to be included already.
#include <cuda_runtime.h>

#include "DeviceCompileCheck.cuh"

#include "ToyDisc3.h"

namespace o2::gtrk::device_check
{

void instantiateToyDisc3()
{
  instantiate<detectors::ToyDisc3>();
}

} // namespace o2::gtrk::device_check
