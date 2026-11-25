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

/// \file TrackerITS3.cxx
/// \brief ITS3's tracking device (TrackerDevice.h), in a unit of its own,
///        instantiating the tracker.

#include "GenericTracking/Detectors/ITS3.h"
#include "GenericTracking/Workflow/TrackerDevice.h"
#include "GenericTracking/WorkflowITSMFT/ConditionsITSMFT.h"

namespace o2::gtrk
{

template o2::framework::DataProcessorSpec workflow::getTrackerSpec<detectors::ITS3>(bool withVertices);

} // namespace o2::gtrk
