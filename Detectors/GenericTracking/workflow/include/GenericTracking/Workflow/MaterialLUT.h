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

/// \file MaterialLUT.h
/// \brief The tracking's material LUT, a CCDB condition per detector.
///
/// | detector | CCDB path |
/// |---|---|
/// | ITS | `ITS/Param/MatLUT` (`o2-gtrk-build-matlut -d ITS`) |
/// | MFT | `MFT/Param/MatLUT` (`o2-gtrk-build-matlut -d MFT`) |
/// | TRKFT3 | `TRK/Param/MatLUT` (`o2-gtrk-build-matlut -d TRKFT3`, one LUT for TRK and FT3) |
/// | TRKFT3TOF | `TRK/Param/MatLUTTOF` (`o2-gtrk-build-matlut -d TRKFT3TOF`, to the TOF, its layers bands of their own) |
///
/// The builder takes a LUT's extent from the layers it serves and makes each
/// barrel layer's material band a run of shells of its own; the tracking
/// refuses a LUT without such runs or one that does not reach a layer's
/// sensors (`detail::lutBands`). O2's global LUT, `GLO/Param/MatLUT`, has no
/// such runs and stays the rest of the reconstruction's.
///
/// `MaterialCondition<Det>` is the condition input `Conditions<Det, true>`
/// holds (Conditions.h): the devices that track fetch it, the measurement
/// producer and the MC vertex reader (`Conditions<Det>`) hold
/// `NoMaterialCondition` in its place and fetch nothing. The LUT arrives
/// through `finaliseCCDB`, rectified as GRPGeomHelper rectifies
/// `GLO/Param/MatLUT`, and the CCDB cache owns it; inside the GPU
/// reconstruction too, as one of the detector's own inputs. Offline, from a
/// CCDB snapshot through DPL's `--condition-remap`; a program without DPL
/// (the bench) reads the same path through `o2::ccdb::BasicCCDBManager`.

#ifndef GTRK_WORKFLOWS_MATERIALLUT_H_
#define GTRK_WORKFLOWS_MATERIALLUT_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Detectors/TRKFT3TOF.h"
#include "GenericTracking/Workflow/DataSpecs.h"

#include "DetectorsBase/MatLayerCylSet.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/InputRecord.h"
#include "Framework/InputSpec.h"
#include "Framework/ProcessingContext.h"
#include "Headers/DataHeader.h"

#include <fairlogger/Logger.h>

#include <concepts>
#include <string>
#include <vector>

namespace o2::gtrk::workflow
{

/// The CCDB path of \a Det's material LUT: TRKFT3TOF's its own, which
/// reaches the TOF, so that TRKFT3's stays TRKFT3's.
template <DetectorLayout Det>
[[nodiscard]] std::string materialLutPath()
{
  return detectorOrigin<Det>().template as<std::string>() + (std::same_as<Det, detectors::TRKFT3TOF> ? "/Param/MatLUTTOF" : "/Param/MatLUT");
}

/// \a Det's material LUT as a condition input (see above).
template <DetectorLayout Det>
class MaterialCondition
{
 public:
  /// Declares the input: called last by the holding `Conditions`'
  /// constructor, so that the LUT's input follows the others.
  static void declare(std::vector<o2::framework::InputSpec>& inputs) { inputs.emplace_back(inputBinding<Det>(kBinding), detectorOrigin<Det>(), "MATLUT", 0, o2::framework::Lifetime::Condition, o2::framework::ccdbParamSpec(materialLutPath<Det>())); }

  /// At the top of the device's run: delivers a new one through finaliseCCDB.
  static void checkUpdates(o2::framework::ProcessingContext& pc) { pc.inputs().get<o2::base::MatLayerCylSet*>(inputBinding<Det>(kBinding).c_str()); }

  /// Whether \a matcher is the LUT, which is then taken.
  bool finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj)
  {
    if (!(matcher == o2::framework::ConcreteDataMatcher(detectorOrigin<Det>(), "MATLUT", 0))) {
      return false;
    }
    mLut = o2::base::MatLayerCylSet::rectifyPtrFromFile(static_cast<o2::base::MatLayerCylSet*>(obj));
    LOGP(info, "{}'s material LUT from the CCDB ({})", Det::name, materialLutPath<Det>());
    return true;
  }

  /// The CCDB cache's; null before the first timeframe.
  [[nodiscard]] const o2::base::MatLayerCylSet* lut() const { return mLut; }

 private:
  static constexpr const char* kBinding = "gtrkMatLUT";

  const o2::base::MatLayerCylSet* mLut = nullptr;
};

/// What a `Conditions` without the material LUT holds in its place: no
/// input, nothing fetched.
struct NoMaterialCondition {
  static void declare(std::vector<o2::framework::InputSpec>& /*inputs*/) {}
  static void checkUpdates(o2::framework::ProcessingContext& /*pc*/) {}
  static bool finaliseCCDB(o2::framework::ConcreteDataMatcher& /*matcher*/, void* /*obj*/) { return false; }
};

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOWS_MATERIALLUT_H_
