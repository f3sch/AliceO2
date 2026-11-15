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

/// \file DeviceStorage.cxx
/// \brief The flat-object functions of DeviceStorage.h for the types the
///        device trackers port: their GPU translation units only declare them.

#include "GenericTracking/Tracking/GPU/DeviceStorage.h"

#include "DetectorsBase/MatLayerCylSet.h"

#include <cstddef>

namespace o2::gtrk
{

template std::size_t flatObjectBufferBytes<o2::base::MatLayerCylSet>(const o2::base::MatLayerCylSet& source);
template void flatObjectImage<o2::base::MatLayerCylSet>(const o2::base::MatLayerCylSet& source, char* futureBuffer, FlatObjectImage& image);

} // namespace o2::gtrk
