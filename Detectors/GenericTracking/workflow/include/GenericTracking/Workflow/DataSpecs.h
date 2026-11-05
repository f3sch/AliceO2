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

/// \file DataSpecs.h
/// \brief The origin and descriptions the GTRK devices exchange.
///
/// One place, so that a producer and a consumer cannot disagree about a
/// spelling that the framework only checks at run time.
///
/// The timeframe's start is o2::framework::TimingInfo's, its length
/// o2::raw::HBFUtils'; the rest is the detector's RecoParam.

#ifndef GTRK_WORKFLOW_DATASPECS_H_
#define GTRK_WORKFLOW_DATASPECS_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Detectors/TRKFT3.h"

#include "Headers/DataHeader.h"

#include <cctype>
#include <concepts>
#include <cstdint>
#include <string>
#include <string_view>

namespace o2::gtrk::workflow
{

inline constexpr o2::header::DataOrigin gDataOrigin{"GTRK"};

/// \a Det's data origin: its readers' and its conditions' (`Conditions<Det>::origin()`).
template <DetectorLayout Det>
[[nodiscard]] o2::header::DataOrigin detectorOrigin()
{
  if constexpr (std::same_as<Det, detectors::ITS>) {
    return o2::header::DataOrigin{"ITS"};
  } else if constexpr (std::same_as<Det, detectors::MFT>) {
    return o2::header::DataOrigin{"MFT"};
  } else {
    static_assert(std::same_as<Det, detectors::TRKFT3>, "the ITS, the MFT and TRKFT3 only");
    return o2::header::DataOrigin{"TRK"};
  }
}

// A detector's messages are GTRK/<DET><kind>, <DET> its data origin
// (`detectorOrigin<Det>()`): GTRK/ITSCLUSTERS, GTRK/MFTMEAS, ..., so that
// the chains of several detectors run in one workflow.

/// A flat array of gtrk::Cluster, one timeframe per message.
inline constexpr std::string_view kClusters{"CLUSTERS"};

/// The MC labels of that array, as an MCTruthContainer; with MC only.
inline constexpr std::string_view kClustersMCTruth{"CLUSMCTR"};

/// One SoA column of Measurement per message part, the subspecification being
/// gtrk::MeasurementColumn, followed by the grid prefix-offset part. DPL
/// frames the parts, so no additional message header is needed.
inline constexpr std::string_view kMeasurements{"MEAS"};

/// What the measurements are relative to (GridConfig.h's MeasurementGeometry):
/// the beam line subtracted from their x and y and the grid they are binned
/// in, the producer's for this timeframe.
inline constexpr std::string_view kMeasurementGeometry{"MEASGEO"};

/// The native vertex seeds, an optional input of the tracking: a flat array of
/// gtrk::VertexSeed sorted by `timePS`, `timePS` and `widthPS` in ps from the
/// timeframe's start, z in cm, `varZ` and `varTransverse` in cm^2 in the frame
/// of the shipped beam line. The study's MC vertex reader fills it.
inline constexpr std::string_view kVertices{"VERTICES"};

// The tracking device's products (TrackerSpec.h), one message each, every
// timeframe, empty when there is nothing.

/// The tracks: a flat array of the tracker's records (`TrackRecord<Det>`),
/// in the tracker's order, every record with its refit's status. A record's
/// clusters are rows of the same timeframe's measurement message.
inline constexpr std::string_view kTracks{"TRACKS"};

/// The own vertexer's vertices: a flat array of gtrk::VertexSeed, sorted by
/// time, then z; empty without the vertexer.
inline constexpr std::string_view kVertexSeeds{"VTX"};

/// Each of those vertices' contributors, a std::uint32_t per vertex.
inline constexpr std::string_view kVertexContributors{"VTXCONTR"};

/// The input binding of the measurement message's part \a part: a
/// gtrk::MeasurementColumn, or the grid's prefix offsets (`kColGridOffsets`).
[[nodiscard]] inline std::string columnBinding(std::uint32_t part) { return "col" + std::to_string(part); }

/// Detector \a Det's description of \a kind.
template <DetectorLayout Det>
[[nodiscard]] o2::header::DataDescription description(std::string_view kind)
{
  o2::header::DataDescription result;
  result.runtimeInit((detectorOrigin<Det>().template as<std::string>() + std::string{kind}).c_str());
  return result;
}

/// Detector \a Det's device \a role: gtrk-<det>-<role>.
template <DetectorLayout Det>
[[nodiscard]] std::string deviceName(std::string_view role)
{
  auto det = detectorOrigin<Det>().template as<std::string>();
  for (auto& c : det) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  return "gtrk-" + det + "-" + std::string{role};
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_DATASPECS_H_
