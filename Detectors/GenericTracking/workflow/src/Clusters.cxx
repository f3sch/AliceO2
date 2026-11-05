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

#include "GenericTracking/Alice/Clusters.h"
#include "GenericTracking/Alice/Placement.h"

#include "DataFormatsITSMFT/ClusterPattern.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "ITSMFTBase/SegmentationAlpide.h"

#include <cstdint>
#include <stdexcept>

namespace o2::gtrk::alice
{

namespace
{

using Segmentation = o2::itsmft::SegmentationAlpide;

// O2's trackers' errors for a cluster without a valid pattern ID: half a pitch.
constexpr float kDefaultError2Row = (Segmentation::PitchRow * 0.5F) * (Segmentation::PitchRow * 0.5F);
constexpr float kDefaultError2Column = (Segmentation::PitchCol * 0.5F) * (Segmentation::PitchCol * 0.5F);

} // namespace

std::vector<LocalCluster> localClusters(std::span<const o2::itsmft::CompClusterExt> clusters, std::span<const unsigned char> patterns,
                                        const o2::itsmft::TopologyDictionary& dictionary)
{
  std::vector<LocalCluster> result;
  result.reserve(clusters.size());
  auto pattern = patterns.begin();
  const auto readPattern = [&pattern, &patterns]() {
    if (pattern == patterns.end()) {
      throw std::runtime_error("cluster patterns ran out");
    }
    return o2::itsmft::ClusterPattern{pattern};
  };
  for (const auto& cluster : clusters) {
    LocalCluster local{.chip = cluster.getSensorID(), .x = 0.F, .z = 0.F, .sigma2X = kDefaultError2Row, .sigma2Z = kDefaultError2Column, .nPixels = 0U};
    o2::math_utils::Point3D<float> position;
    // The pattern's spans: a common topology's from the dictionary, a group's
    // and an invalid ID's from the stream.
    const auto spans = [&local](const o2::itsmft::ClusterPattern& patt) {
      local.rowSpan = static_cast<std::uint8_t>(patt.getRowSpan());
      local.colSpan = static_cast<std::uint8_t>(patt.getColumnSpan());
    };
    const auto id = cluster.getPatternID();
    if (id != o2::itsmft::CompCluster::InvalidPatternID) {
      local.sigma2X = dictionary.getErr2X(id);
      local.sigma2Z = dictionary.getErr2Z(id);
      if (!dictionary.isGroup(id)) {
        position = dictionary.getClusterCoordinates(cluster);
        local.nPixels = static_cast<std::uint32_t>(dictionary.getNpixels(id));
        spans(dictionary.getPattern(id));
      } else {
        const auto patt = readPattern();
        position = o2::itsmft::TopologyDictionary::getClusterCoordinates(cluster, patt);
        local.nPixels = static_cast<std::uint32_t>(patt.getNPixels());
        spans(patt);
      }
    } else {
      const auto patt = readPattern();
      position = o2::itsmft::TopologyDictionary::getClusterCoordinates(cluster, patt, false);
      local.nPixels = static_cast<std::uint32_t>(patt.getNPixels());
      spans(patt);
    }
    local.x = position.X();
    local.z = position.Z();
    result.push_back(local);
  }
  return result;
}

std::vector<Cluster> toClusters(std::span<const LocalCluster> local, std::span<const o2::itsmft::ROFRecord> rofs, ReadoutFrames frames,
                                std::uint32_t tfFirstOrbit)
{
  return detail::toClusters(local, rofs, frames, tfFirstOrbit, kActiveCentreX);
}

ReadoutFrames readoutFrames(std::span<const o2::itsmft::ROFRecord> rofs) { return detail::readoutFrames(rofs); }

} // namespace o2::gtrk::alice
