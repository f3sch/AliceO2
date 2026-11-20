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

/// \file TrkClusterConverter.h
/// \brief ALICE 3's cluster converter, TRKFT3's and TRKFT3TOF's (TRK.cxx,
///        TRKTOF.cxx): the library's own, no installed header.

#ifndef GTRK_WORKFLOWTRK_TRKCLUSTERCONVERTER_H_
#define GTRK_WORKFLOWTRK_TRKCLUSTERCONVERTER_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Detectors/TRKFT3TOF.h"
#include "GenericTracking/Clustering/Cluster.h"
#include "GenericTracking/Workflow/DataSpecs.h"
#include "GenericTracking/WorkflowTRK/ConditionsTRK.h"
#include "GenericTracking/WorkflowTRK/TRK.h"
#include "GenericTracking/WorkflowTRK/TRKTOF.h"

#include "DataFormatsIOTOF/Cluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsTRKFT3/Cluster.h"
#include "DataFormatsTRKFT3/ROFRecord.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DetectorsRaw/HBFUtils.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Logger.h"
#include "Framework/Task.h"
#include "Framework/TimingInfo.h"
#include "Headers/DataHeader.h"
#include "IOTOFBase/GeometryTGeo.h"
#include "SimulationDataFormat/ConstMCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "TRKBase/GeometryTGeo.h"

#include <gsl/span>

#include <concepts>
#include <cstdint>
#include <format>
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace o2::gtrk::workflow
{

/// ALICE 3's TRK and FT3 clusters, one stream per barrel layer and one per
/// disc as O2's cluster readers publish them (TRK/COMPCLUSTERS/<layer>,
/// FT3/COMPCLUSTERS/<disc> and the rest), turned into GTRK/TRKCLUSTERS and
/// GTRK/TRKCLUSMCTR: each layer's
/// stream converted with its own readout frames, each cluster with its
/// position and errors from the topology dictionaries, CCDB conditions,
/// then the layers
/// concatenated in order, barrel first, so a measurement's cluster index is
/// its position in that concatenation and the labels follow it. For
/// TRKFT3TOF the TOF's clusters (TF3/CLUSTERS/0, TRKTOF.h's `tofClusters`)
/// follow the TRK's, with their labels.
template <DetectorLayout Det>
class TrkClusterConverterTask : public o2::framework::Task
{
 public:
  static constexpr bool kTOF = std::same_as<Det, detectors::TRKFT3TOF>;
  /// TRK's and FT3's streams, TRKFT3's layers.
  static constexpr int kLayers = Detector<detectors::TRKFT3>::nLayers;
  static constexpr int kBarrelLayers = Detector<detectors::TRKFT3>::nBarrelLayers;

  /// TRKFT3's layer \a layer's stream: TRK's barrel layer, or FT3's disc after them.
  [[nodiscard]] static o2::header::DataOrigin originOf(int layer) { return layer < kBarrelLayers ? o2::header::DataOrigin{"TRK"} : o2::header::DataOrigin{"FT3"}; }
  [[nodiscard]] static o2::header::DataHeader::SubSpecificationType subspecOf(int layer) { return static_cast<o2::header::DataHeader::SubSpecificationType>(layer < kBarrelLayers ? layer : layer - kBarrelLayers); }

  TrkClusterConverterTask(bool useMC, TrkDictionaryConditions dictionaries, std::shared_ptr<o2::base::GRPGeomRequest> geometryRequest)
    : mUseMC(useMC), mDictionaries(dictionaries), mGeometryRequest(std::move(geometryRequest)) {}

  void init(o2::framework::InitContext& /*ic*/) final { o2::base::GRPGeomHelper::instance().setRequest(mGeometryRequest); }

  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj) final
  {
    if (!o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
      mDictionaries.finaliseCCDB(matcher, obj);
    }
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    o2::base::GRPGeomHelper::instance().checkUpdates(pc);
    // Dictionary coordinates use the VD segmentation's petal count, which
    // GeometryTGeo derives from the loaded geometry in this device's process.
    o2::trk::GeometryTGeo::Instance();
    TrkDictionaryConditions::checkUpdates(pc);
    const auto& timing = pc.services().get<o2::framework::TimingInfo>();
    std::vector<Cluster> clusters;
    o2::dataformats::MCTruthContainer<o2::MCCompLabel> labels;
    for (int layer = 0; layer < kLayers; ++layer) {
      const auto patterns = pc.inputs().get<gsl::span<unsigned char>>(std::format("patterns_{}", layer));
      const auto rofs = pc.inputs().get<gsl::span<o2::trkft3::ROFRecord>>(std::format("rofs_{}", layer));
      if (rofs.empty()) {
        continue;
      }
      // The timeframe's first orbit as DPL announces it, held to the one the
      // frames imply, as for the ITS.
      const auto implied = o2::raw::HBFUtils::Instance().getFirstIRofTF(rofs.front().getBCData()).orbit;
      if (implied != timing.firstTForbit) {
        LOGP(fatal, "TRKFT3 layer {}: the readout frames belong to the timeframe starting at orbit {}, DPL says {}", layer, implied, timing.firstTForbit);
      }
      const std::span<const unsigned char> patternBytes{patterns.data(), patterns.size()};
      const auto local = [&]() {
        if (layer < kBarrelLayers) {
          const auto compact = pc.inputs().get<gsl::span<o2::trkft3::TRKCluster>>(std::format("compclusters_{}", layer));
          return alice::trkLocalClusters(std::span{compact.data(), compact.size()}, patternBytes, *mDictionaries.trk());
        }
        const auto compact = pc.inputs().get<gsl::span<o2::trkft3::FT3Cluster>>(std::format("compclusters_{}", layer));
        return alice::ft3LocalClusters(std::span{compact.data(), compact.size()}, patternBytes, *mDictionaries.ft3());
      }();
      const std::span<const o2::trkft3::ROFRecord> frames{rofs.data(), rofs.size()};
      const auto converted = alice::toClusters(local, frames, alice::readoutFrames(frames), timing.firstTForbit);
      clusters.insert(clusters.end(), converted.begin(), converted.end());
      if (mUseMC) {
        const auto layerLabels = pc.inputs().get<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>(std::format("labels_{}", layer));
        if (layerLabels->getIndexedSize() != converted.size()) {
          LOGP(fatal, "TRKFT3 layer {}: {} labelled clusters for {} clusters", layer, layerLabels->getIndexedSize(), converted.size());
        }
        labels.mergeAtBack(*layerLabels);
      }
    }
    if constexpr (kTOF) {
      appendTOF(pc, timing.firstTForbit, clusters, labels);
    }
    LOGP(info, "{}: {} clusters on {} layers for the timeframe at orbit {}", Det::name, clusters.size(), Detector<Det>::nLayers, timing.firstTForbit);
    pc.outputs().snapshot(o2::framework::Output{gDataOrigin, description<Det>(kClusters), 0}, clusters);
    if (mUseMC) {
      auto& flat = pc.outputs().make<o2::dataformats::ConstMCTruthContainer<o2::MCCompLabel>>(o2::framework::Output{gDataOrigin, description<Det>(kClustersMCTruth), 0});
      labels.flatten_to(flat);
    }
  }

  /// The converter's inputs, outputs and device (`getClusterConverterSpec`).
  [[nodiscard]] static o2::framework::DataProcessorSpec spec(bool useMC)
  {
    using namespace o2::framework;
    std::vector<InputSpec> inputs;
    std::vector<OutputSpec> outputs{OutputSpec{gDataOrigin, description<Det>(kClusters), 0, Lifetime::Timeframe}};
    for (int layer = 0; layer < kLayers; ++layer) {
      const auto origin = originOf(layer);
      const auto subspec = subspecOf(layer);
      inputs.emplace_back(std::format("compclusters_{}", layer), origin, "COMPCLUSTERS", subspec, Lifetime::Timeframe);
      inputs.emplace_back(std::format("patterns_{}", layer), origin, "PATTERNS", subspec, Lifetime::Timeframe);
      inputs.emplace_back(std::format("rofs_{}", layer), origin, "CLUSTERSROF", subspec, Lifetime::Timeframe);
      if (useMC) {
        inputs.emplace_back(std::format("labels_{}", layer), origin, "CLUSTERSMCTR", subspec, Lifetime::Timeframe);
      }
    }
    if constexpr (kTOF) {
      inputs.emplace_back("tofclusters", o2::header::gDataOriginTF3, "CLUSTERS", 0, Lifetime::Timeframe);
      inputs.emplace_back("tofrofs", o2::header::gDataOriginTF3, "CLUSTERSROF", 0, Lifetime::Timeframe);
      if (useMC) {
        inputs.emplace_back("toflabels", o2::header::gDataOriginTF3, "CLUSTERSMCTR", 0, Lifetime::Timeframe);
      }
    }
    if (useMC) {
      outputs.emplace_back(gDataOrigin, description<Det>(kClustersMCTruth), 0, Lifetime::Timeframe);
    }
    const TrkDictionaryConditions dictionaries(inputs);
    auto geometryRequest = std::make_shared<o2::base::GRPGeomRequest>(false, false, false, false, false, o2::base::GRPGeomRequest::Aligned, inputs);
    return DataProcessorSpec{
      .name = deviceName<Det>("cluster-converter"),
      .inputs = inputs,
      .outputs = outputs,
      .algorithm = AlgorithmSpec{adaptFromTask<TrkClusterConverterTask>(useMC, dictionaries, geometryRequest)},
      .options = Options{},
    };
  }

 private:
  /// The TOF's clusters after the TRK's, and their labels: those
  /// `tofClusters` keeps, the others counted in one warning per timeframe.
  void appendTOF(o2::framework::ProcessingContext& pc, std::uint32_t tfFirstOrbit, std::vector<Cluster>& clusters, o2::dataformats::MCTruthContainer<o2::MCCompLabel>& labels)
  {
    if (!mTOFFrame) {
      auto& geometry = *o2::iotof::GeometryTGeo::Instance();
      mTOFFrame = alice::TOFPixelFrame::build(geometry);
    }
    // The timeframe's first orbit, held to the one the TOF's readout frame
    // implies, as for the TRK's layers.
    const auto tofFrames = pc.inputs().get<gsl::span<o2::itsmft::ROFRecord>>("tofrofs");
    if (!tofFrames.empty()) {
      const auto implied = o2::raw::HBFUtils::Instance().getFirstIRofTF(tofFrames.front().getBCData()).orbit;
      if (implied != tfFirstOrbit) {
        LOGP(fatal, "TRKFT3TOF: the TOF's readout frames belong to the timeframe starting at orbit {}, DPL says {}", implied, tfFirstOrbit);
      }
    }
    const auto tof = pc.inputs().get<gsl::span<o2::iotof::Cluster>>("tofclusters");
    std::vector<std::uint32_t> kept;
    alice::TOFSkipped skipped;
    const auto converted = alice::tofClusters(std::span{tof.data(), tof.size()}, *mTOFFrame, tfFirstOrbit, o2::raw::HBFUtils::Instance().getNOrbitsPerTF(), kept, skipped);
    if (skipped.total() != 0 || skipped.clamped != 0) {
      LOGP(warn, "TRKFT3TOF: of {} TOF clusters {} left out (no chip {}, time in its crossing {}, pattern {}, outside the timeframe {}), {} boxes started at the timeframe's start", tof.size(), skipped.total(), skipped.chip, skipped.time, skipped.pattern, skipped.outside, skipped.clamped);
    }
    if (mUseMC) {
      const auto tofLabels = pc.inputs().get<const o2::dataformats::MCTruthContainer<o2::MCCompLabel>*>("toflabels");
      if (tofLabels->getIndexedSize() != tof.size()) {
        LOGP(fatal, "TRKFT3TOF: {} labelled TOF clusters for {} clusters", tofLabels->getIndexedSize(), tof.size());
      }
      const auto first = clusters.size();
      for (std::size_t i = 0; i < kept.size(); ++i) {
        const auto own = tofLabels->getLabels(kept[i]);
        labels.addElements(first + i, gsl::span<const o2::MCCompLabel>{own.data(), own.size()});
      }
    }
    clusters.insert(clusters.end(), converted.begin(), converted.end());
  }

  bool mUseMC;
  TrkDictionaryConditions mDictionaries;
  std::shared_ptr<o2::base::GRPGeomRequest> mGeometryRequest;
  std::optional<alice::TOFPixelFrame> mTOFFrame; ///< TRKFT3TOF's, from the first timeframe's geometry
};

/// `getClusterConverterSpec` of TRKFT3 and TRKFT3TOF.
template <DetectorLayout Det>
[[nodiscard]] o2::framework::DataProcessorSpec trkClusterConverterSpec(bool useMC)
{
  return TrkClusterConverterTask<Det>::spec(useMC);
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOWTRK_TRKCLUSTERCONVERTER_H_
