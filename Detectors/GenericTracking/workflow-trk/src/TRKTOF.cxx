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

/// \file TRKTOF.cxx
/// \brief ALICE 3's tracker with its TOF layers (TRKTOF.h): TRKFT3TOF's chips,
///        the IOTOF clusters as ours, frames and depths, its cluster
///        converter and its measurement producer.

#include "GenericTracking/WorkflowTRK/TRKTOF.h"

#include "GenericTracking/Base/TimeEstPS.h"
#include "GenericTracking/Detectors/TRKFT3TOF.h"
#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Workflow/ClusterConverterSpec.h"
#include "GenericTracking/Workflow/GridConfig.h"
#include "GenericTracking/Workflow/MeasurementProducer.h"
#include "GenericTracking/Workflow/MeasurementProducerSpec.h"
#include "GenericTracking/WorkflowTRK/ConditionsTRK.h"

#include "TrkClusterConverter.h"

#include "CommonConstants/LHCConstants.h"
#include "DataFormatsIOTOF/Cluster.h"
#include "IOTOFBase/GeometryTGeo.h"
#include "IOTOFBase/IOTOFBaseParam.h"
#include "MathUtils/Cartesian.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace o2::gtrk::alice
{

template <>
std::vector<Sensor> placeChips<detectors::TRKFT3TOF>()
{
  using Layout = detectors::TRKFT3TOF;
  auto sensors = placeChips<detectors::TRKFT3>();
  for (auto& sensor : sensors) {
    sensor.layer = detectors::trkft3tof::fromTRKFT3(sensor.layer);
  }
  // The IOTOF chips: the plane of O2's sensor volume through the active
  // matrix's centre, loc0 along its x (rows), loc1 along its z (columns), the
  // normal x cross z, its -y: `tofClusters`' frame (TofPixelFrame).
  auto& tof = *o2::iotof::GeometryTGeo::Instance();
  const auto& specs = o2::iotof::ChipSpecificsParam::Instance();
  const auto first = tof.detectorToLocal(0, 0, 0);
  const auto last = tof.detectorToLocal(specs.NRows - 1, specs.NCols - 1, 0);
  const double cx = (first.X() + last.X()) / 2.F;
  const double cz = (first.Z() + last.Z()) / 2.F;
  const int chips = tof.getITOFNumberOfChips() + tof.getOTOFNumberOfChips();
  for (int chip = 0; chip < chips; ++chip) {
    const auto* placement = tof.extractMatrixSensor(chip);
    if (placement == nullptr) {
      throw std::runtime_error("TRKFT3TOF: IOTOF chip " + std::to_string(chip) + " has no sensor matrix");
    }
    const o2::math_utils::Transform3D matrix(*placement);
    using Point = o2::math_utils::Point3D<double>;
    const auto origin = matrix * Point{cx, 0., cz};
    const auto axis0 = matrix * Point{cx + 1., 0., cz};
    const auto axis1 = matrix * Point{cx, 0., cz + 1.};
    const double u[3] = {axis0.X() - origin.X(), axis0.Y() - origin.Y(), axis0.Z() - origin.Z()};
    const double v[3] = {axis1.X() - origin.X(), axis1.Y() - origin.Y(), axis1.Z() - origin.Z()};
    const double n[3] = {(u[1] * v[2]) - (u[2] * v[1]), (u[2] * v[0]) - (u[0] * v[2]), (u[0] * v[1]) - (u[1] * v[0])};
    Sensor sensor{};
    for (int i = 0; i < 3; ++i) {
      sensor.rot[(3 * i) + 0] = u[i];
      sensor.rot[(3 * i) + 1] = v[i];
      sensor.rot[(3 * i) + 2] = n[i];
    }
    sensor.trans[0] = origin.X();
    sensor.trans[1] = origin.Y();
    sensor.trans[2] = origin.Z();
    sensor.layer = tof.getIOTOFLayer(chip) == 0 ? Layout::kInnerTOFLayer : Layout::kOuterTOFLayer;
    sensors.push_back(sensor);
  }
  return sensors;
}

TofPixelFrame TofPixelFrame::build(const o2::iotof::GeometryTGeo& geometry)
{
  const auto& specs = o2::iotof::ChipSpecificsParam::Instance();
  TofPixelFrame frame;
  while (geometry.isValidIOTOFChipIndex(static_cast<int>(frame.chips))) {
    ++frame.chips;
  }
  if (frame.chips == 0 || specs.NRows <= 0 || specs.NCols <= 0) {
    throw std::runtime_error(std::format("IOTOF: {} chips of {} x {} pixels in the geometry", frame.chips, specs.NRows, specs.NCols));
  }
  // The pixel centres in the chip frame: x depends on the row alone, z on
  // the column alone (GeometryTGeo::detectorToLocal), about the active
  // matrix's centre (`placeChips<TRKFT3TOF>`).
  const auto first = geometry.detectorToLocal(0, 0, 0);
  const auto last = geometry.detectorToLocal(specs.NRows - 1, specs.NCols - 1, 0);
  const double cx = (first.X() + last.X()) / 2.;
  const double cz = (first.Z() + last.Z()) / 2.;
  frame.rowX.resize(static_cast<std::size_t>(specs.NRows));
  frame.columnZ.resize(static_cast<std::size_t>(specs.NCols));
  for (int row = 0; row < specs.NRows; ++row) {
    frame.rowX[static_cast<std::size_t>(row)] = static_cast<float>(geometry.detectorToLocal(row, 0, 0).X() - cx);
  }
  for (int column = 0; column < specs.NCols; ++column) {
    frame.columnZ[static_cast<std::size_t>(column)] = static_cast<float>(geometry.detectorToLocal(0, column, 0).Z() - cz);
  }
  frame.sigmaX = static_cast<float>(specs.PitchRow / std::sqrt(12.));
  frame.sigmaZ = static_cast<float>(specs.PitchCol / std::sqrt(12.));
  return frame;
}

std::vector<Cluster> tofClusters(std::span<const o2::iotof::Cluster> clusters, const TofPixelFrame& frame, std::uint32_t tfFirstOrbit, int orbitsPerTimeframe, std::vector<std::uint32_t>& kept, TofSkipped& skipped)
{
  constexpr double kBunchNS = o2::constants::lhc::LHCBunchSpacingNS;
  const auto firstBC = static_cast<std::int64_t>(tfFirstOrbit) * o2::constants::lhc::LHCMaxBunches;
  const std::int64_t durationPS = workflow::timeframeDurationPS(orbitsPerTimeframe);
  const std::int64_t widthPS = psOfBC(1);
  std::vector<Cluster> result;
  result.reserve(clusters.size());
  kept.clear();
  kept.reserve(clusters.size());
  for (std::size_t index = 0; index < clusters.size(); ++index) {
    const auto& cluster = clusters[index];
    const auto chip = cluster.getChipID();
    if (chip >= frame.chips) {
      ++skipped.chip;
      continue;
    }
    const float fineNS = cluster.getTimeWrtBc();
    if (!(fineNS >= 0.F && static_cast<double>(fineNS) < kBunchNS)) {
      ++skipped.time;
      continue;
    }
    // The fired pixels' mean centre: the pattern is the cluster's
    // rows x columns, row-major, from its first row and column.
    const unsigned rows = cluster.getRowSpan();
    const unsigned columns = cluster.getColSpan();
    const unsigned mask = cluster.getPattern();
    if (mask == o2::iotof::Cluster::InvalidPatternID || rows == 0 || columns == 0 || rows * columns > 16 || cluster.getRow() + rows > frame.rowX.size() || cluster.getCol() + columns > frame.columnZ.size()) {
      ++skipped.pattern;
      continue;
    }
    double x = 0.;
    double z = 0.;
    int pixels = 0;
    for (unsigned bit = 0; bit < rows * columns; ++bit) {
      if ((mask & (1U << bit)) != 0U) {
        x += frame.rowX[cluster.getRow() + (bit / columns)];
        z += frame.columnZ[cluster.getCol() + (bit % columns)];
        ++pixels;
      }
    }
    if (pixels == 0) {
      ++skipped.pattern;
      continue;
    }
    // Its time from the timeframe's start, the bunch crossings' difference
    // taken first; the box the crossing that ends there, in the timeframe.
    const std::int64_t timePS = psOfBC(static_cast<std::int64_t>(cluster.getBc()) - firstBC) + std::llround(static_cast<double>(fineNS) * 1000.);
    std::int64_t start = timePS - widthPS;
    if (start < 0 || timePS > durationPS) {
      ++skipped.clamped;
      start = std::clamp<std::int64_t>(start, 0, durationPS - widthPS);
    }
    result.push_back({
      .chipID = static_cast<SensorID>(detectors::TRKFT3TOF::kFirstTOFSensor + static_cast<int>(chip)),
      .xLocal = static_cast<float>(x / pixels),
      .yLocal = static_cast<float>(z / pixels),
      .sigmaXLocal = frame.sigmaX,
      .sigmaYLocal = frame.sigmaZ,
      .timePS = static_cast<std::uint32_t>(start),
      .widthPS = static_cast<std::uint32_t>(widthPS),
      .rowSpan = static_cast<std::uint8_t>(rows),
      .colSpan = static_cast<std::uint8_t>(columns),
    });
    kept.push_back(static_cast<std::uint32_t>(index));
  }
  // In box-start order, as the measurement producer takes a layer timed per
  // hit (MeasurementProducer.h): its stable counting sort keeps every
  // bucket's rows in it. The clusters come in the clusterer's order.
  std::vector<std::uint32_t> order(result.size());
  std::iota(order.begin(), order.end(), 0U);
  std::ranges::stable_sort(order, {}, [&](std::uint32_t i) { return result[i].timePS; });
  std::vector<Cluster> sorted(result.size());
  std::vector<std::uint32_t> sortedKept(kept.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    sorted[i] = result[order[i]];
    sortedKept[i] = kept[order[i]];
  }
  kept = std::move(sortedKept);
  return sorted;
}

std::vector<FrameTiming> trktofFrameTiming(std::span<const FrameTiming> trk)
{
  using Layout = detectors::TRKFT3TOF;
  if (trk.size() != static_cast<std::size_t>(Detector<detectors::TRKFT3>::nLayers)) {
    throw std::invalid_argument(std::format("TRKFT3TOF: {} TRKFT3 layers' frames, the layout has {}", trk.size(), Detector<detectors::TRKFT3>::nLayers));
  }
  std::vector<FrameTiming> frames(static_cast<std::size_t>(Detector<Layout>::nLayers), FrameTiming{.lengthBC = 1, .firstStartBC = 0});
  for (std::size_t layer = 0; layer < trk.size(); ++layer) {
    frames[detectors::trkft3tof::fromTRKFT3(static_cast<LayerID>(layer))] = trk[layer];
  }
  return frames;
}

std::vector<double> trktofMeasurementDepthsCM(std::span<const double> trk)
{
  using Layout = detectors::TRKFT3TOF;
  std::vector<double> depths(static_cast<std::size_t>(Detector<Layout>::nLayers), 0.);
  for (std::size_t layer = 0; layer < trk.size() && std::cmp_less(layer, Detector<detectors::TRKFT3>::nLayers); ++layer) {
    depths[detectors::trkft3tof::fromTRKFT3(static_cast<LayerID>(layer))] = trk[layer];
  }
  return depths;
}

} // namespace o2::gtrk::alice

namespace o2::gtrk::workflow
{

template <>
o2::framework::DataProcessorSpec getClusterConverterSpec<detectors::TRKFT3TOF>(const o2::framework::ConfigContext& /*context*/, bool useMC)
{
  return trkClusterConverterSpec<detectors::TRKFT3TOF>(useMC);
}

template o2::framework::DataProcessorSpec getMeasurementProducerSpec<detectors::TRKFT3TOF>();

} // namespace o2::gtrk::workflow
