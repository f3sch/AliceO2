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

/// \file TRK.cxx
/// \brief ALICE 3's own, TRK's barrel and FT3's discs: its chips placed, its
///        cluster format, its dictionaries' depths, frames and field (TRK.h),
///        its cluster converter (TrkClusterConverter.h) and its measurement
///        producer.

#include "GenericTracking/WorkflowTRK/TRK.h"

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Alice/Clusters.h"
#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Clustering/Cluster.h"
#include "GenericTracking/Workflow/ClusterConverterSpec.h"
#include "GenericTracking/Workflow/DataSpecs.h"
#include "GenericTracking/Workflow/MeasurementProducer.h"
#include "GenericTracking/Workflow/MeasurementProducerSpec.h"
#include "GenericTracking/WorkflowTRK/ConditionsTRK.h"

#include "TrkClusterConverter.h"

#include "DataFormatsTRKFT3/Cluster.h"
#include "DataFormatsTRKFT3/ROFRecord.h"
#include "DetectorsRaw/HBFUtils.h"
#include "FT3Base/GeometryTGeo.h"
#include "Field/ALICE3MagneticField.h"
#include "Framework/ConcreteDataMatcher.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Logger.h"
#include "Framework/Task.h"
#include "Framework/TimingInfo.h"
#include "MathUtils/Cartesian.h"
#include "MathUtils/Utils.h"
#include "SimulationDataFormat/ConstMCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "TRKBase/AlmiraParam.h"
#include "TRKBase/GeometryTGeo.h"
#include "TRKBase/SegmentationChip.h"
#include "TRKBase/Specs.h"
#include "TRKReconstruction/TopologyDictionary.h"

#include <gsl/span>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <initializer_list>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace o2::gtrk::alice
{

namespace
{

// TRKFT3.h restates O2's TRK sensors; here they are held to O2's, to the bit.
namespace specs = o2::trk::constants;
// Our VD sensors are four petals (`trkft3::vdWidth`); O2 takes the petal
// count from the loaded geometry, so `placeChips` holds the widths to O2's.
static_assert(detectors::trkft3::kVDRadius == specs::VD::petal::layer::radii);
static_assert(detectors::trkft3::kVDPetalGap == specs::VD::petal::layer::gaps);
static_assert(exactlyEqual(detectors::trkft3::kVDLength, specs::VD::petal::layer::length));
static_assert(exactlyEqual(detectors::trkft3::kVDThickness, specs::VD::silicon::thickness));
static_assert(exactlyEqual(detectors::trkft3::kVDPitch, specs::VD::petal::layer::pitchX) && exactlyEqual(detectors::trkft3::kVDPitch, specs::VD::petal::layer::pitchZ));
static_assert(exactlyEqual(detectors::trkft3::kChipWidth, specs::moduleMLOT::chip::width));
static_assert(exactlyEqual(detectors::trkft3::kChipPassiveEdgeReadOut, specs::moduleMLOT::chip::passiveEdgeReadOut));
static_assert(exactlyEqual(detectors::trkft3::kChipLength, specs::moduleMLOT::chip::length));
static_assert(exactlyEqual(detectors::trkft3::kChipThickness, specs::moduleMLOT::silicon::thickness));
static_assert(exactlyEqual(detectors::trkft3::kChipPitch, specs::moduleMLOT::chip::pitchX) && exactlyEqual(detectors::trkft3::kChipPitch, specs::moduleMLOT::chip::pitchZ));

} // namespace

template <>
std::vector<Sensor> placeChips<detectors::TRKFT3>()
{
  // The VD's curved sensors (sub-detector 0) and the ML/OT chips, whose
  // sensor volume is the active matrix, so their origin needs no shift; then
  // FT3's chips, the same chip with the same volume, their layers after the
  // barrel's.
  auto& trk = *o2::trk::GeometryTGeo::Instance();
  auto& ft3 = *o2::ft3::GeometryTGeo::Instance();
  trk.fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  ft3.fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  // O2's VD widths follow the loaded geometry's petal count; ours are four petals'.
  for (std::size_t vdLayer = 0; vdLayer < detectors::trkft3::kVDRadius.size(); ++vdLayer) {
    const int layer = static_cast<int>(vdLayer);
    const float width = o2::trk::SegmentationChip::getWidthVD(layer);
    if (!exactlyEqual(static_cast<float>(detectors::trkft3::vdWidth(layer)), width)) {
      throw std::runtime_error(std::format("VD layer {}: the geometry's sensors are {} cm wide, ours (four petals) {} cm", layer, width, detectors::trkft3::vdWidth(layer)));
    }
  }
  std::vector<Sensor> sensors;
  sensors.reserve(static_cast<std::size_t>(trk.getNumberOfChips()) + static_cast<std::size_t>(ft3.getNumberOfChips()));
  for (int chip = 0; chip < trk.getNumberOfChips(); ++chip) {
    const auto& l2g = trk.getMatrixL2G(chip);
    const int layer = trk.getLayerTRK(chip);
    sensors.push_back(trk.getSubDetID(chip) == 0 ? curvedSensor(l2g, layer) : flatSensor(l2g, layer, 0.));
  }
  const int firstDiscLayer = Detector<detectors::TRKFT3>::nBarrelLayers;
  for (int chip = 0; chip < ft3.getNumberOfChips(); ++chip) {
    sensors.push_back(flatSensor(ft3.getMatrixL2G(chip), firstDiscLayer + ft3.getLayer(chip), 0.));
  }
  return sensors;
}

template <>
std::array<SensorSpec, Detector<detectors::TRKFT3>::nLayers> placeSpecs<detectors::TRKFT3>()
{
  return Geometry<detectors::TRKFT3>::layoutSpecs();
}

namespace
{

/// The local clusters of TRK's or FT3's \a clusters, each chip numbered
/// \a firstSensor plus O2's chip ID; \a subDetectors, the sub-detectors their
/// clusters may name (a cluster naming another is refused).
template <typename ClusterType, typename Dictionary>
std::vector<LocalCluster> alice3LocalClusters(std::span<const ClusterType> clusters, std::span<const unsigned char> patterns, const Dictionary& dictionary, std::uint32_t firstSensor, std::initializer_list<int> subDetectors)
{
  std::vector<LocalCluster> result;
  result.reserve(clusters.size());
  std::size_t offset = 0;
  for (const auto& cluster : clusters) {
    if (std::ranges::find(subDetectors, static_cast<int>(cluster.subDetID)) == subDetectors.end()) {
      throw std::runtime_error(std::format("a cluster on chip {} names sub-detector {}", cluster.chipID, cluster.subDetID));
    }
    // The pattern: its row span, its column span, then its pixels' bits row
    // by row, the first pixel in the highest bit (O2's Clusterer<TRK>).
    if (offset + 2U > patterns.size()) {
      throw std::runtime_error("TRK cluster patterns ran out");
    }
    const auto* pattern = patterns.data() + offset;
    const auto rowSpan = pattern[0];
    const auto colSpan = pattern[1];
    const auto bytes = 2U + ((static_cast<std::size_t>(rowSpan) * colSpan + 7U) / 8U);
    if (offset + bytes > patterns.size()) {
      throw std::runtime_error("TRK cluster patterns ran out");
    }
    offset += bytes;
    int fired = 0;
    for (std::size_t byte = 2; byte < bytes; ++byte) {
      fired += std::popcount(pattern[byte]);
    }
    LocalCluster local{.chip = firstSensor + cluster.chipID, .nPixels = static_cast<std::uint32_t>(fired), .rowSpan = rowSpan, .colSpan = colSpan};
    const int layer = Dictionary::layerOf(cluster);
    if (layer < 0 || layer >= dictionary.getNumberOfLayers()) {
      throw std::runtime_error(std::format("a cluster on chip {} is on no layer of the topology dictionary", cluster.chipID));
    }
    const int id = dictionary.findID(layer, pattern);
    const auto position = dictionary.getFlatCoordinates(cluster, pattern, id);
    const auto& entry = dictionary.getEntry(layer, id);
    local.x = position.X();
    local.z = position.Z();
    local.sigma2X = entry.xWidth * entry.xWidth;
    local.sigma2Z = entry.zWidth * entry.zWidth;
    result.push_back(local);
  }
  if (offset != patterns.size()) {
    throw std::runtime_error(std::format("{} bytes of TRK cluster patterns left over", patterns.size() - offset));
  }
  return result;
}

} // namespace

std::vector<LocalCluster> trkLocalClusters(std::span<const o2::trkft3::TRKCluster> clusters, std::span<const unsigned char> patterns, const TrkDictionary& dictionary) { return alice3LocalClusters(clusters, patterns, dictionary, 0U, {0, 1}); }

std::vector<LocalCluster> ft3LocalClusters(std::span<const o2::trkft3::FT3Cluster> clusters, std::span<const unsigned char> patterns, const Ft3Dictionary& dictionary) { return alice3LocalClusters(clusters, patterns, dictionary, static_cast<std::uint32_t>(detectors::TRKFT3::kFirstDiscSensor), {2}); }

std::vector<Cluster> toClusters(std::span<const LocalCluster> local, std::span<const o2::trkft3::ROFRecord> rofs, ReadoutFrames frames, std::uint32_t tfFirstOrbit)
{
  // The TRK's sensor volume is its active matrix: our origin is O2's (TRK.cxx's placeChips).
  return detail::toClusters(local, rofs, frames, tfFirstOrbit, 0.);
}

ReadoutFrames readoutFrames(std::span<const o2::trkft3::ROFRecord> rofs) { return detail::readoutFrames(rofs); }

std::vector<double> measurementDepthsCM(const TrkDictionary& trk, const Ft3Dictionary& ft3)
{
  std::vector<double> depths;
  depths.reserve(static_cast<std::size_t>(trk.getNumberOfLayers()) + static_cast<std::size_t>(ft3.getNumberOfLayers()));
  for (int layer = 0; layer < trk.getNumberOfLayers(); ++layer) {
    depths.push_back(-static_cast<double>(trk.getLayer(layer).depth));
  }
  for (int layer = 0; layer < ft3.getNumberOfLayers(); ++layer) {
    depths.push_back(-static_cast<double>(ft3.getLayer(layer).depth));
  }
  return depths;
}

std::vector<FrameTiming> trkFrameTiming()
{
  const auto& param = o2::trk::AlmiraParam::Instance();
  constexpr int kLayers = Detector<detectors::TRKFT3>::nLayers;
  // Our layers are the parameter's entries: TRK's barrel layers, then FT3's
  // discs, which take entries 11-22.
  static_assert(kLayers == static_cast<int>(o2::trk::AlmiraParam::kNLayers));
  std::vector<FrameTiming> frames;
  frames.reserve(kLayers);
  for (int layer = 0; layer < kLayers; ++layer) {
    frames.push_back({.lengthBC = param.getROFLengthInBC(layer), .firstStartBC = param.getROFBiasInBC(layer) + param.getROFDelayInBC(layer)});
  }
  return frames;
}

double alice3BzKiloGauss()
{
  const char* macro = std::getenv("ALICE3_MAGFIELD_MACRO");
  if (macro == nullptr || !std::filesystem::exists(macro)) {
    throw std::runtime_error(std::string{"ALICE3_MAGFIELD_MACRO names no file ("} + (macro == nullptr ? "unset" : macro) + "): the ALICE 3 field is unknown, and O2 would take a toy one");
  }
  o2::field::ALICE3MagneticField field;
  const double origin[3] = {0., 0., 0.};
  double b[3] = {0., 0., 0.};
  field.Field(origin, b);
  return b[2];
}

} // namespace o2::gtrk::alice

namespace o2::gtrk::workflow
{

template <>
o2::framework::DataProcessorSpec getClusterConverterSpec<detectors::TRKFT3>(const o2::framework::ConfigContext& /*context*/, bool useMC)
{
  return trkClusterConverterSpec<detectors::TRKFT3>(useMC);
}

template o2::framework::DataProcessorSpec getMeasurementProducerSpec<detectors::TRKFT3>();

} // namespace o2::gtrk::workflow
