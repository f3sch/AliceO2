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

/// \file ClusterConverterITS3.h
/// \brief ITS3's cluster converter (ITS3.cxx): the library's own, no
///        installed header.

#ifndef GTRK_WORKFLOWITS3_CLUSTERCONVERTERITS3_H_
#define GTRK_WORKFLOWITS3_CLUSTERCONVERTERITS3_H_

#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Detectors/ITS3.h"
#include "GenericTracking/Workflow/DataSpecs.h"
#include "GenericTracking/WorkflowITS3/Clusters.h"

#include "ITSMFTInput.h"

#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Logger.h"
#include "Framework/Task.h"
#include "Framework/TimingInfo.h"
#include "ITS3Reconstruction/TopologyDictionary.h"
#include "ITSBase/GeometryTGeo.h"
#include "SimulationDataFormat/ConstMCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

#include <cstddef>
#include <cstdint>
#include <format>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace o2::gtrk::workflow
{

/// O2's ITS3 clusters, as its clusterer and cluster reader publish them --
/// the ITS's origin, ITS/COMPCLUSTERS and the rest, one stream or one per
/// layer -- turned into GTRK/IT3CLUSTERS and GTRK/IT3CLUSMCTR (Clusters.h),
/// in O2's order, so that a measurement's cluster index addresses the
/// streams concatenated in layer order. With ITS3's topology dictionary
/// (IT3/Calib/ClusterDictionary), the ITS's frames (ITS/Config/AlpideParam)
/// and geometry (ITS/Config/Geometry, which builds ITS3's), whose tiles'
/// orientation the conversion takes. No frame selection.
class ITS3ClusterConverterTask : public o2::framework::Task
{
 public:
  using AlpideParam = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>;

  ITS3ClusterConverterTask(bool useMC, std::int32_t streams) : mUseMC(useMC), mStreams(streams) {}

  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj) final
  {
    if (matcher == o2::framework::ConcreteDataMatcher(o2::header::gDataOriginITS, "ALPIDEPARAM", 0)) {
      return; // the CCDB fetch updates the singleton
    }
    if (matcher == o2::framework::ConcreteDataMatcher(o2::header::gDataOriginIT3, "CLUSDICT", 0)) {
      mDictionary = static_cast<const o2::its3::TopologyDictionary*>(obj);
      LOGP(info, "ITS3 topology dictionary from the CCDB: {} inner-layer, {} outer-barrel patterns", mDictionary->getSize(true), mDictionary->getSize(false));
      return;
    }
    if (matcher == o2::framework::ConcreteDataMatcher(o2::header::gDataOriginITS, "GEOMTGEO", 0)) {
      o2::its::GeometryTGeo::adopt(static_cast<o2::its::GeometryTGeo*>(obj));
      mTiles.reset();
      LOGP(info, "ITS3 GeometryTGeo from the CCDB adopted");
    }
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    // Each delivers a new object through finaliseCCDB.
    pc.inputs().get<o2::its3::TopologyDictionary*>("dictionary");
    pc.inputs().get<o2::its::GeometryTGeo*>("geomTGeo");
    pc.inputs().get<AlpideParam*>("alpideParam");
    if (!mTiles) {
      mTiles.emplace(*o2::its::GeometryTGeo::Instance());
    }
    const auto input = readITSMFTInput(pc, mUseMC, mStreams);
    const auto& param = AlpideParam::Instance();
    const auto orbit = pc.services().get<o2::framework::TimingInfo>().firstTForbit;
    std::vector<Cluster> clusters;
    clusters.reserve(input.size());
    std::size_t rofCount = 0U;
    for (std::uint32_t stream = 0; stream < input.streams().size(); ++stream) {
      const auto& source = input.streams()[stream];
      // One stream holds every layer, whose frames must then be alike.
      const std::int32_t lengthBC = param.getROFLengthInBC(static_cast<std::int32_t>(stream));
      for (std::int32_t layer = 0; mStreams == 1 && layer < Detector<detectors::ITS3>::nLayers; ++layer) {
        if (param.getROFLengthInBC(layer) != lengthBC) {
          throw std::runtime_error("staggered ITS3 timing needs per-layer cluster input (--enable-its-staggering)");
        }
      }
      if (source.rofs.size() > 1U && alice::readoutFrames(source.rofs).lengthBC != lengthBC) {
        throw std::runtime_error(std::format("ITS3 stream {}: frames other than the run's length of {} BC", stream, lengthBC));
      }
      const auto local = alice::its3LocalClusters(source.compact, source.patterns, *mDictionary, *mTiles);
      const auto converted = alice::its3Clusters(local, source.rofs, {.lengthBC = lengthBC}, orbit);
      clusters.insert(clusters.end(), converted.begin(), converted.end());
      rofCount += source.rofs.size();
    }
    LOGP(info, "ITS3: {} clusters in {} readout frames from {} streams for the timeframe at orbit {}", clusters.size(), rofCount, mStreams, orbit);
    pc.outputs().snapshot(o2::framework::Output{gDataOrigin, description<detectors::ITS3>(kClusters), 0}, clusters);
    if (mUseMC) {
      auto& flat = pc.outputs().make<o2::dataformats::ConstMCTruthContainer<o2::MCCompLabel>>(o2::framework::Output{gDataOrigin, description<detectors::ITS3>(kClustersMCTruth), 0});
      concatenateLabels(input).flatten_to(flat);
    }
  }

  /// The converter's device; \a streams 1, or one per layer (staggered).
  [[nodiscard]] static o2::framework::DataProcessorSpec spec(bool useMC, std::int32_t streams)
  {
    using namespace o2::framework;
    std::vector<InputSpec> inputs;
    addITSMFTClusterInputs(inputs, o2::header::gDataOriginITS, useMC, streams);
    inputs.emplace_back("dictionary", o2::header::gDataOriginIT3, "CLUSDICT", 0, Lifetime::Condition, ccdbParamSpec("IT3/Calib/ClusterDictionary"));
    inputs.emplace_back("alpideParam", o2::header::gDataOriginITS, "ALPIDEPARAM", 0, Lifetime::Condition, ccdbParamSpec("ITS/Config/AlpideParam"));
    inputs.emplace_back("geomTGeo", o2::header::gDataOriginITS, "GEOMTGEO", 0, Lifetime::Condition, ccdbParamSpec("ITS/Config/Geometry"));
    std::vector<OutputSpec> outputs{OutputSpec{gDataOrigin, description<detectors::ITS3>(kClusters), 0, Lifetime::Timeframe}};
    if (useMC) {
      outputs.emplace_back(gDataOrigin, description<detectors::ITS3>(kClustersMCTruth), 0, Lifetime::Timeframe);
    }
    return DataProcessorSpec{
      .name = deviceName<detectors::ITS3>("cluster-converter"),
      .inputs = inputs,
      .outputs = outputs,
      .algorithm = AlgorithmSpec{adaptFromTask<ITS3ClusterConverterTask>(useMC, streams)},
      .options = Options{},
    };
  }

 private:
  bool mUseMC;
  std::int32_t mStreams;
  const o2::its3::TopologyDictionary* mDictionary = nullptr; ///< the CCDB cache's
  std::optional<alice::ITS3Tiles> mTiles;                    ///< from the adopted geometry, at its first timeframe
};

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOWITS3_CLUSTERCONVERTERITS3_H_
