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

/// \file ConditionsITSMFT.h
/// \brief The ITS's, ITS3's and the MFT's run conditions for a GTRK device
///        (`Conditions`, Workflow/Conditions.h).
///
/// The ITS's and MFT's are what O2's own trackers take by default: not the
/// full TGeo but the detector's lightweight `GeometryTGeo` object
/// (`ITS/Config/Geometry`, `MFT/Config/Geometry`), with its chip matrices
/// filled, which the placement reads like a loaded geometry; the field
/// (`GLO/Config/GRPMagField`) through GRPGeomHelper; and `<DET>/Config/AlpideParam`
/// for the frames. ITS3's are the ITS's (O2's ITS3 comes with the ITS's
/// origin and GeometryTGeo), with its tiles' measurement depths. Offline,
/// from CCDB snapshots through DPL's `--condition-remap`.

#ifndef GTRK_WORKFLOWITSMFT_CONDITIONSITSMFT_H_
#define GTRK_WORKFLOWITSMFT_CONDITIONSITSMFT_H_

#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/ITS3.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Workflow/Conditions.h"
#include "GenericTracking/Workflow/DataSpecs.h"
#include "GenericTracking/Workflow/MaterialLUT.h"

#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/InitContext.h"
#include "Framework/InputRecord.h"
#include "Framework/InputSpec.h"
#include "Framework/ProcessingContext.h"
#include "Headers/DataHeader.h"
#include "ITSBase/GeometryTGeo.h"
#include "MFTBase/GeometryTGeo.h"

#include <fairlogger/Logger.h>

#include <concepts>
#include <cstddef>
#include <format>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace o2::gtrk::workflow
{

/// The ITS's, ITS3's and the MFT's.
template <DetectorLayout Det, bool WithMaterial, ConditionsHolder Holder>
  requires std::same_as<Det, detectors::ITS> || std::same_as<Det, detectors::ITS3> || std::same_as<Det, detectors::MFT>
class Conditions<Det, WithMaterial, Holder>
{
  static constexpr bool kITS3 = std::same_as<Det, detectors::ITS3>;
  static constexpr bool kITS = std::same_as<Det, detectors::ITS> || kITS3; ///< the ITS's GeometryTGeo and frames
  static constexpr bool kOwnDevice = Holder == ConditionsHolder::kOwnDevice;

 public:
  explicit Conditions(std::vector<o2::framework::InputSpec>& inputs)
  {
    if constexpr (kOwnDevice) {
      mRequest = std::make_shared<o2::base::GRPGeomRequest>(false, true, false, true, false, o2::base::GRPGeomRequest::None, inputs);
    }
    const auto detector = origin().template as<std::string>();
    inputs.emplace_back(inputBinding<Det>("geomTGeo"), origin(), "GEOMTGEO", 0, o2::framework::Lifetime::Condition, o2::framework::ccdbParamSpec(detector + "/Config/Geometry"));
    inputs.emplace_back(inputBinding<Det>("alpideParam"), origin(), "ALPIDEPARAM", 0, o2::framework::Lifetime::Condition, o2::framework::ccdbParamSpec(detector + "/Config/AlpideParam"));
    MaterialConditionOf<Det, WithMaterial>::declare(inputs);
  }

  static constexpr bool kMeanVertex = true;

  /// The origin of O2's data and conditions: ITS3's are the ITS's.
  [[nodiscard]] static o2::header::DataOrigin origin()
  {
    if constexpr (kITS3) {
      return o2::header::gDataOriginITS;
    } else {
      return detectorOrigin<Det>();
    }
  }

  void init(o2::framework::InitContext& /*ic*/)
  {
    if constexpr (kOwnDevice) {
      o2::base::GRPGeomHelper::instance().setRequest(mRequest);
    }
  }

  void checkUpdates(o2::framework::ProcessingContext& pc) const
  {
    if constexpr (kOwnDevice) {
      o2::base::GRPGeomHelper::instance().checkUpdates(pc);
    }
    pc.inputs().get<GeometryTGeoType*>(inputBinding<Det>("geomTGeo").c_str());
    pc.inputs().get<AlpideParam*>(inputBinding<Det>("alpideParam").c_str()); // the CCDB fetcher updates the singleton
    MaterialConditionOf<Det, WithMaterial>::checkUpdates(pc);
  }

  /// Whether \a matcher was one of these conditions, then taken.
  bool finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj)
  {
    if (mMaterial.finaliseCCDB(matcher, obj)) {
      return true;
    }
    if constexpr (kOwnDevice) {
      if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
        return true;
      }
    }
    if (matcher == o2::framework::ConcreteDataMatcher(origin(), "GEOMTGEO", 0)) {
      GeometryTGeoType::adopt(static_cast<GeometryTGeoType*>(obj));
      LOGP(info, "{} GeometryTGeo from the CCDB adopted", Det::name);
      return true;
    }
    return false;
  }

  /// Frame k starts at k * length + delay + bias BC, independently per layer.
  [[nodiscard]] std::vector<FrameTiming> frameTiming() const
  {
    const auto& param = AlpideParam::Instance();
    std::vector<FrameTiming> result;
    result.reserve(Detector<Det>::nLayers);
    for (int layer = 0; layer < Detector<Det>::nLayers; ++layer) {
      const auto length = param.getROFLengthInBC(layer);
      if (length <= 0) {
        throw std::runtime_error(std::format("{} layer {}: AlpideParam gives readout frames of {} BC", Det::name, layer, length));
      }
      result.push_back({.lengthBC = length, .firstStartBC = param.getROFDelayInBC(layer) + param.getROFBiasInBC(layer)});
    }
    return result;
  }

  /// The run's orbits per timeframe, its GRPECS's, as O2's ITS tracking
  /// takes them.
  [[nodiscard]] static int orbitsPerTimeframe() { return o2::base::GRPGeomHelper::getNHBFPerTF(); }

  /// None: the ALPIDE's clusters are at the effective sensitive layer's
  /// centre, 0 on every layer; on ITS3's tiles O2's response depth.
  [[nodiscard]] std::vector<double> measurementDepthsCM() const
  {
    if constexpr (kITS3) {
      std::vector<double> depths(Detector<Det>::nLayers, 0.);
      for (std::size_t layer = 0; layer < detectors::mosaix::kNominalRadius.size(); ++layer) {
        depths[layer] = detectors::mosaix::measurementDepth(static_cast<int>(layer));
      }
      return depths;
    } else {
      return {};
    }
  }

  /// The run's material LUT, the CCDB cache's; null before the first timeframe.
  [[nodiscard]] const o2::base::MatLayerCylSet* materialLut() const
    requires WithMaterial
  {
    return mMaterial.lut();
  }

 private:
  using GeometryTGeoType = std::conditional_t<kITS, o2::its::GeometryTGeo, o2::mft::GeometryTGeo>;
  using AlpideParam = o2::itsmft::DPLAlpideParam<kITS ? o2::detectors::DetID::ITS : o2::detectors::DetID::MFT>;

  std::shared_ptr<o2::base::GRPGeomRequest> mRequest; ///< the run-wide conditions' request, in a device of its own
  [[no_unique_address]] MaterialConditionOf<Det, WithMaterial> mMaterial;
};

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOWITSMFT_CONDITIONSITSMFT_H_
