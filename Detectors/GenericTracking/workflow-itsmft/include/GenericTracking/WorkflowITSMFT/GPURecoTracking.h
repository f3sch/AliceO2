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

/// \file GPURecoTracking.h
/// \brief The ITS's and the MFT's tracking inside the GPU reconstruction
///        (Workflow/GPURecoTracking.h).

#ifndef GTRK_WORKFLOWITSMFT_GPURECOTRACKING_H_
#define GTRK_WORKFLOWITSMFT_GPURECOTRACKING_H_

#include "GenericTracking/Workflow/GPURecoTracking.h"

#include <memory>
#include <string_view>

namespace o2::gtrk::workflow
{

/// The tracking of \a detector, "its" or "mft"; std::invalid_argument for
/// another.
[[nodiscard]] std::unique_ptr<GPURecoTracking> makeGPURecoTracking(std::string_view detector);

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOWITSMFT_GPURECOTRACKING_H_
