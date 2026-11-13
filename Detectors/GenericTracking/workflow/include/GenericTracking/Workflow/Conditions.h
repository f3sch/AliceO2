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

/// \file Conditions.h
/// \brief The run's conditions for a GTRK device, per detector:
///        `Conditions<Det, WithMaterial>`, declared here; the ITS's and the
///        MFT's in WorkflowITSMFT/ConditionsITSMFT.h, ALICE 3's in
///        WorkflowTRK/ConditionsTRK.h.
///
/// Every detector's give a device, from the first timeframe on, what it
/// places its layout from (`alice::placeChips`, which reads the detector's
/// GeometryTGeo), the field as the propagator's nominal Bz (GridConfig.h's
/// `nominalBzKiloGauss`), its readout frames per layer and its layers'
/// measurement depths; they do not change within a run (a device sees one).
/// The MeanVertex does, and is the measurement producer's own input, per
/// timeframe, where the detector has one (`kMeanVertex`):
///
/// - `Conditions(inputs)` declares its condition inputs at construction;
/// - `init(ic)` in the device's init, `checkUpdates(pc)` at the top of its
///   run, `finaliseCCDB(matcher, obj)` forwarded from the device's;
/// - held by the GPU reconstruction (`ConditionsHolder::kGPUReco`), only the
///   detector's own: the GRPECS, the field and the ITS's LUT are that
///   device's (a second GRPGeomHelper request is fatal);
/// - `frameTiming()`; `measurementDepthsCM()`, empty for 0 on every layer;
///   `origin()`, the detector's data origin; `kMeanVertex`, whether the
///   detector's beam line is the run's MeanVertex (`GLO/Calib/MeanVertex`);
/// - with `WithMaterial`, the material LUT (MaterialLUT.h) as one more
///   condition, `materialLut()`: the devices that track take it, the
///   measurement producer and the MC vertex reader (the default) fetch none.

#ifndef GTRK_WORKFLOWS_CONDITIONS_H_
#define GTRK_WORKFLOWS_CONDITIONS_H_

#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Workflow/DataSpecs.h"
#include "GenericTracking/Workflow/MaterialLUT.h"

#include "DetectorsBase/GRPGeomHelper.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/InitContext.h"
#include "Framework/InputRecord.h"
#include "Framework/InputSpec.h"
#include "Framework/ProcessingContext.h"
#include "Headers/DataHeader.h"

#include <fairlogger/Logger.h>

#include <concepts>
#include <cstdint>
#include <format>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace o2::gtrk::workflow
{

/// Who holds a `Conditions`: a device of its own, which requests the run's
/// GRPECS and field itself (GRPGeomHelper), or the GPU reconstruction, whose
/// device requests them -- and GLO/Param/MatLUT, the ITS's material LUT --
/// for everything it runs, and fetches and finalises them.
enum class ConditionsHolder : std::uint8_t { kOwnDevice,
                                             kGPUReco };

/// The conditions a device of \a Det reads (see above); with \a WithMaterial
/// also the material LUT; held by \a Holder.
template <DetectorLayout Det, bool WithMaterial = false, ConditionsHolder Holder = ConditionsHolder::kOwnDevice>
class Conditions;

/// The material LUT's condition a `Conditions` holds: the LUT's, GRPGeomHelper's
/// for the ITS inside the GPU reconstruction, or none.
template <DetectorLayout Det, bool WithMaterial, ConditionsHolder Holder = ConditionsHolder::kOwnDevice>
using MaterialConditionOf = std::conditional_t<!WithMaterial, NoMaterialCondition, std::conditional_t<Holder == ConditionsHolder::kGPUReco && std::same_as<Det, detectors::ITS>, GRPGeomMaterialCondition, MaterialCondition<Det>>>;

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOWS_CONDITIONS_H_
