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
/// \brief The ITS's and MFT's compact clusters, `o2::itsmft::CompClusterExt`,
///        as our clusters: the two steps of Alice/Clusters.h.

#ifndef GTRK_WORKFLOWITSMFT_CLUSTERS_H_
#define GTRK_WORKFLOWITSMFT_CLUSTERS_H_

#include "GenericTracking/Alice/Clusters.h"
#include "GenericTracking/Clustering/Cluster.h"

#include <cstdint>
#include <span>
#include <vector>

namespace o2::itsmft
{
class CompClusterExt;
class ROFRecord;
class TopologyDictionary;
} // namespace o2::itsmft

namespace o2::gtrk::alice
{

/// Every compact cluster's local position and errors, in order, reading the
/// patterns of the group and invalid ones as it goes. Throws
/// std::runtime_error for a pattern ID beyond \a dictionary; patterns that
/// run out before a cluster's fail gsl::span's bounds check.
[[nodiscard]] std::vector<LocalCluster> localClusters(std::span<const o2::itsmft::CompClusterExt> clusters, std::span<const unsigned char> patterns, const o2::itsmft::TopologyDictionary& dictionary);

/// Our clusters from the local ones, in the same order. \a tfFirstOrbit is
/// the timeframe's first orbit, which a cluster's time is measured from.
/// Throws std::runtime_error if a ROF does not cover its clusters exactly and
/// in order, or starts before the timeframe.
[[nodiscard]] std::vector<Cluster> toClusters(std::span<const LocalCluster> local, std::span<const o2::itsmft::ROFRecord> rofs, ReadoutFrames frames, std::uint32_t tfFirstOrbit);

/// The readout frames as \a rofs lay them out: O2 writes one record per frame
/// of the timeframe, empty or not, so their spacing is the frame length. That
/// is the length the data were taken with, whatever a configuration dump says:
/// in the MC the ITS's comes from CCDB (ITS/Config/AlpideParam, 198 BC) while
/// the workflows' .ini files print the parameter's default, 891. Throws
/// std::runtime_error unless there are two frames or more, evenly spaced.
[[nodiscard]] ReadoutFrames readoutFrames(std::span<const o2::itsmft::ROFRecord> rofs);

} // namespace o2::gtrk::alice

#endif // GTRK_WORKFLOWITSMFT_CLUSTERS_H_
