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

#include "GenericTracking/WorkflowITSMFT/Placement.h"

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Detectors/Alpide.h"

#include "ITSMFTBase/SegmentationAlpide.h"

namespace o2::gtrk::alice
{

namespace
{

using Segmentation = o2::itsmft::SegmentationAlpide;

// Alpide.h restates O2's numbers so that the layouts need no O2 detector
// library; this is where they are held to them.
static_assert(detectors::alpide::kRows == Segmentation::NRows);
static_assert(detectors::alpide::kColumns == Segmentation::NCols);
static_assert(exactlyEqual(static_cast<float>(detectors::alpide::kPitchRow), Segmentation::PitchRow));
static_assert(exactlyEqual(static_cast<float>(detectors::alpide::kPitchColumn), Segmentation::PitchCol));
static_assert(exactlyEqual(static_cast<float>(detectors::alpide::kEffectiveThickness), Segmentation::SensorLayerThicknessEff));

} // namespace

const double kActiveCentreX = 0.5 * (static_cast<double>(Segmentation::PassiveEdgeReadOut) - static_cast<double>(Segmentation::PassiveEdgeTop));

} // namespace o2::gtrk::alice
