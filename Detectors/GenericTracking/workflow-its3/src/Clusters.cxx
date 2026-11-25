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

#include "GenericTracking/WorkflowITS3/Clusters.h"

#include "GenericTracking/WorkflowITSMFT/Placement.h"

#include "DataFormatsITSMFT/ClusterPattern.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITS3Base/SpecsV2.h"
#include "ITS3Reconstruction/IOUtils.h"
#include "ITS3Reconstruction/TopologyDictionary.h"
#include "ITSBase/GeometryTGeo.h"
#include "MathUtils/Utils.h"

#include <gsl/span>

#include <cmath>
#include <cstddef>
#include <format>
#include <stdexcept>

namespace o2::gtrk::alice
{

namespace
{
namespace specs = o2::its3::constants;
namespace mosaix = detectors::mosaix;
} // namespace

ITS3Tiles::ITS3Tiles(o2::its::GeometryTGeo& geometry) : mRightHanded(specs::detID::nChips, 0U)
{
  geometry.fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  for (std::uint32_t chip = 0; chip < specs::detID::nChips; ++chip) {
    double r[9];
    double t[3];
    geometry.getMatrixL2G(static_cast<int>(chip)).GetComponents(r[0], r[1], r[2], t[0], r[3], r[4], r[5], t[1], r[6], r[7], r[8], t[2]);
    const double determinant = (r[0] * ((r[4] * r[8]) - (r[5] * r[7]))) - (r[1] * ((r[3] * r[8]) - (r[5] * r[6]))) + (r[2] * ((r[3] * r[7]) - (r[4] * r[6])));
    if (std::abs(std::abs(determinant) - 1.) > 1e-9) {
      throw std::runtime_error(std::format("ITS3 tile {}: its L2G matrix has determinant {}", chip, determinant));
    }
    mRightHanded[chip] = determinant > 0. ? 1U : 0U;
  }
}

std::array<double, 3> tileUnrolled(const ITS3Tiles& tiles, std::uint32_t chip, double x, double y, double z)
{
  const int layer = specs::detID::getDetID2Layer(static_cast<int>(chip));
  const double inner = mosaix::innerRadius(layer);
  // Our angle from the arc's middle runs with O2's on a right-handed tile.
  const double angle = std::atan2(y, x) - mosaix::middleAngle(layer);
  return {(tiles.rightHanded(chip) ? angle : -angle) * inner, z, std::hypot(x, y) - inner};
}

std::vector<LocalCluster> its3LocalClusters(std::span<const o2::itsmft::CompClusterExt> clusters, std::span<const unsigned char> patterns, const o2::its3::TopologyDictionary& dictionary, const ITS3Tiles& tiles)
{
  std::vector<LocalCluster> result;
  result.reserve(clusters.size());
  const gsl::span<const unsigned char> stream{patterns.data(), patterns.size()};
  auto pattern = stream.begin();
  for (const auto& cluster : clusters) {
    const auto chip = static_cast<std::uint32_t>(cluster.getSensorID());
    const bool tile = specs::detID::isDetITS3(static_cast<int>(chip));
    const auto id = cluster.getPatternID();
    if (id != o2::itsmft::CompCluster::InvalidPatternID && id >= dictionary.getSize(tile)) {
      throw std::runtime_error(std::format("ITS3 chip {}: pattern ID {} beyond the dictionary's {}", chip, id, dictionary.getSize(tile)));
    }
    const auto read = pattern;
    float sigma2X = 0.F;
    float sigma2Z = 0.F;
    const auto point = o2::its3::ioutils::extractClusterData(cluster, pattern, &dictionary, sigma2X, sigma2Z);
    LocalCluster local{.chip = chip, .sigma2X = sigma2X, .sigma2Z = sigma2Z};
    // The pattern: a group's and an invalid ID's from the stream, read
    // again where the decode read it; a common topology's from the dictionary.
    const auto take = [&local](const o2::itsmft::ClusterPattern& patt) {
      local.nPixels = static_cast<std::uint32_t>(patt.getNPixels());
      local.rowSpan = static_cast<std::uint8_t>(patt.getRowSpan());
      local.colSpan = static_cast<std::uint8_t>(patt.getColumnSpan());
    };
    if (pattern != read) {
      auto again = read;
      take(o2::itsmft::ClusterPattern{again});
    } else {
      take(dictionary.getPattern(id, tile));
    }
    if (tile) {
      const auto unrolled = tileUnrolled(tiles, chip, point.X(), point.Y(), point.Z());
      const int layer = specs::detID::getDetID2Layer(static_cast<int>(chip));
      const double scale = mosaix::innerRadius(layer) / mosaix::middleRadius(layer);
      local.x = static_cast<float>(unrolled[0]);
      local.z = static_cast<float>(unrolled[1]);
      local.sigma2X = static_cast<float>(static_cast<double>(sigma2X) * scale * scale);
    } else {
      local.x = static_cast<float>(static_cast<double>(point.X()) - kActiveCentreX);
      local.z = point.Z();
    }
    result.push_back(local);
  }
  return result;
}

std::vector<Cluster> its3Clusters(std::span<const LocalCluster> local, std::span<const o2::itsmft::ROFRecord> rofs, ReadoutFrames frames, std::uint32_t tfFirstOrbit) { return detail::toClusters(local, rofs, frames, tfFirstOrbit, 0.); }

std::array<double, Detector<detectors::ITS3>::nLayers> its3MeasurementDepthsCM()
{
  std::array<double, Detector<detectors::ITS3>::nLayers> depths{};
  for (int layer = 0; layer < static_cast<int>(specs::nLayers); ++layer) {
    depths[static_cast<std::size_t>(layer)] = mosaix::middleRadius(layer) - mosaix::innerRadius(layer) + specs::pixelarray::pixels::apts::responseYShift;
  }
  return depths;
}

} // namespace o2::gtrk::alice
