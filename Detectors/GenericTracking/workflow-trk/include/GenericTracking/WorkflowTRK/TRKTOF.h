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

/// \file TRKTOF.h
/// \brief ALICE 3's tracker with its TOF layers (Detectors/TRKFT3TOF.h): the
///        IOTOF clusters as ours, the layout's frames and depths.
///
/// **Clusters.** `o2::iotof::Cluster`, one stream (TF3/CLUSTERS/0), each with
/// its absolute bunch crossing and its time in it (`getBc`, `getTimeWrtBc`):
/// `tofClusters` gives each its pixels' mean centre in its chip's frame -- the
/// active matrix's centre, x across the rows, z along the columns, as
/// `placeChips<TRKFT3TOF>` places the chip -- pitch^2/12 variances, and the
/// time box of the bunch crossing that ends at its time, [t - 1 BC, t].
///
/// **Frames and depths.** TRKFT3's on its layers (TRK.h), one bunch crossing
/// at measurement depth 0 on the TOF's: a TOF row's box is its own, never a
/// readout frame.

#ifndef GTRK_ALICE_TRKTOF_H_
#define GTRK_ALICE_TRKTOF_H_

#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Clustering/Cluster.h"

#include <cstdint>
#include <span>
#include <vector>

namespace o2::iotof
{
class Cluster;
class GeometryTGeo;
} // namespace o2::iotof

namespace o2::gtrk::alice
{

/// Where an IOTOF chip's pixels lie in its frame (`tofClusters`): every
/// chip's the same, from O2's IOTOF GeometryTGeo and ChipSpecificsParam.
struct TOFPixelFrame {
  std::vector<float> rowX;    ///< a row's pixel centre, along the frame's x [cm]
  std::vector<float> columnZ; ///< a column's, along its z [cm]
  float sigmaX = 0.F;         ///< the row pitch over sqrt(12) [cm]
  float sigmaZ = 0.F;         ///< the column pitch's
  std::uint32_t chips = 0;    ///< the IOTOF chips, iTOF and oTOF

  /// From \a geometry: the frame's origin the active matrix's centre. Throws
  /// std::runtime_error if the geometry has no chips or pixels.
  [[nodiscard]] static TOFPixelFrame build(const o2::iotof::GeometryTGeo& geometry);
};

/// The clusters `tofClusters` left out, by why.
struct TOFSkipped {
  std::uint32_t chip = 0;    ///< on no IOTOF chip
  std::uint32_t time = 0;    ///< a time in its bunch crossing outside [0, 1 BC)
  std::uint32_t pattern = 0; ///< too large for a pattern (O2's `kHuge`), no pixel fired, or one outside the chip
  std::uint32_t outside = 0; ///< its time outside the timeframe
  std::uint32_t clamped = 0; ///< kept: its box's start moved to the timeframe's start
  [[nodiscard]] std::uint32_t total() const { return chip + time + pattern + outside; }
};

/// Every usable cluster of \a clusters as ours, in its box's start order
/// (stable), and its index in \a clusters into \a kept: sensor `TRKFT3TOF::kFirstTOFSensor` plus its
/// chip, its position and errors in \a frame, its time box [t - 1 BC, t] in
/// the timeframe starting at orbit \a tfFirstOrbit of \a orbitsPerTimeframe
/// orbits, moved into the timeframe where it reaches past one of its ends
/// (counted). A cluster on no chip, with a time outside its bunch crossing or
/// without a usable pattern is left out and counted in \a skipped.
[[nodiscard]] std::vector<Cluster> tofClusters(std::span<const o2::iotof::Cluster> clusters, const TOFPixelFrame& frame, std::uint32_t tfFirstOrbit, std::int32_t orbitsPerTimeframe, std::vector<std::uint32_t>& kept, TOFSkipped& skipped);

/// TRKFT3TOF's readout frames: TRKFT3's (\a trk, `trkFrameTiming`) on its
/// layers, one bunch crossing from BC 0 on the TOF's.
[[nodiscard]] std::vector<FrameTiming> trktofFrameTiming(std::span<const FrameTiming> trk);

/// TRKFT3TOF's measurement depths [cm]: TRKFT3's (\a trk,
/// `measurementDepthsCM`) on its layers, 0 on the TOF's, whose clusters lie
/// on the sensor's mid-plane.
[[nodiscard]] std::vector<double> trktofMeasurementDepthsCM(std::span<const double> trk);

} // namespace o2::gtrk::alice

#endif // GTRK_ALICE_TRKTOF_H_
