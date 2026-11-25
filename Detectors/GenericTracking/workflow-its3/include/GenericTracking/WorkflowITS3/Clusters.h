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

/// \file Clusters.h
/// \brief ITS3's compact clusters, `o2::itsmft::CompClusterExt` with ITS3's
///        topology dictionary, as our clusters: the two steps of
///        Alice/Clusters.h, in our sensors' frames.
///
/// A tile's cluster comes from O2's dictionary as a point of the tile's
/// cylinder frame (`its3::ioutils::extractClusterData`); ours is its unrolled
/// position on our tile (ITS3.cxx `tileSensor`): the arc at the inner radius
/// from the tile's middle, against O2's flat coordinate, which O2 measures at
/// the middle radius from the far end, and with O2's angle on a mirrored tile
/// against ours. An outer-barrel cluster is the ITS's: ALICE's local x less
/// the active matrix's offset (WorkflowITSMFT/Placement.h).

#ifndef GTRK_WORKFLOWITS3_CLUSTERS_H_
#define GTRK_WORKFLOWITS3_CLUSTERS_H_

#include "GenericTracking/Detectors/ITS3.h"
#include "GenericTracking/Alice/Clusters.h"
#include "GenericTracking/Clustering/Cluster.h"

#include <array>
#include <cstdint>
#include <span>
#include <vector>

namespace o2::itsmft
{
class CompClusterExt;
class ROFRecord;
} // namespace o2::itsmft

namespace o2::its3
{
class TopologyDictionary;
} // namespace o2::its3

namespace o2::its
{
class GeometryTGeo;
} // namespace o2::its

namespace o2::gtrk::alice
{

/// Each tile's orientation in O2's geometry, the loaded one's: whether its
/// L2G is right-handed (O2 mirrors every other tile).
class ITS3Tiles
{
 public:
  /// From \a geometry's L2G matrices; throws std::runtime_error if a tile's
  /// is not a rotation or a reflection.
  explicit ITS3Tiles(o2::its::GeometryTGeo& geometry);

  /// Whether tile \a chip's frame is right-handed.
  [[nodiscard]] bool rightHanded(std::uint32_t chip) const { return mRightHanded[chip] != 0U; }

 private:
  std::vector<std::uint8_t> mRightHanded; ///< by O2's chip ID
};

/// Our unrolled position on tile \a chip of a point of O2's tile frame,
/// (\a x, \a y) across its cylinder and \a z along it: across, along, and
/// off the inner surface [cm].
[[nodiscard]] std::array<double, 3> tileUnrolled(const ITS3Tiles& tiles, std::uint32_t chip, double x, double y, double z);

/// Every compact cluster's position and errors in our sensor's frame, in
/// order, reading the patterns of the group and invalid ones as it goes. An
/// entry without errors in the dictionary (O2's ITS3 dictionary has rare
/// entries at zero) takes O2's default, half a pitch.
/// Throws std::runtime_error for a pattern ID beyond \a dictionary.
[[nodiscard]] std::vector<LocalCluster> its3LocalClusters(std::span<const o2::itsmft::CompClusterExt> clusters, std::span<const unsigned char> patterns, const o2::its3::TopologyDictionary& dictionary, const ITS3Tiles& tiles);

/// Our clusters from those of `its3LocalClusters`, in the same order, as
/// Alice/Clusters.h's `toClusters`.
[[nodiscard]] std::vector<Cluster> its3Clusters(std::span<const LocalCluster> local, std::span<const o2::itsmft::ROFRecord> rofs, ReadoutFrames frames, std::uint32_t tfFirstOrbit);

} // namespace o2::gtrk::alice

#endif // GTRK_WORKFLOWITS3_CLUSTERS_H_
