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

/// \file ShapeShares.h
/// \brief How often the cluster spans occur, a timeframe's: for a particle's
///        own crossing, from the detector's cluster-shape calibration, and for
///        any row of a layer, from the timeframe's own rows.

#ifndef GTRK_TRACKING_SHAPESHARES_H_
#define GTRK_TRACKING_SHAPESHARES_H_

#include "GenericTracking/Base/LayerSpec.h"

#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace o2::gtrk
{

/// The shares of the cluster spans the follow's last-row trim reads
/// (`logSizeRatio`), per layer of the layout, as the logarithms of the
/// shares, spans from `kSpans` - 1 on together:
///  - a particle's own crossing (the detector's cluster-shape calibration,
///    ClusterShapeBand.h): per bin of the slope (`edges`, `bins` of them),
///    its column span at its slope along the columns and its row span at its
///    slope along the rows; NaN on a layer or in a bin the calibration does
///    not know;
///  - any row of the layer in the timeframe, its column span and its row
///    span.
/// `bins` 0: no calibration (a detector without a band), no shares read.
/// The measurement producer makes them per timeframe and ships them with its
/// grid (`MeasurementGeometry`); the frame carries them (TimeFrame.h).
template <DetectorLayout Det>
struct ShapeShares {
  static constexpr std::uint32_t kSpans = 16U;
  static constexpr std::uint32_t kMaxBins = 64U;
  static constexpr std::size_t kLayers = static_cast<std::size_t>(Detector<Det>::nLayers);

  std::uint32_t bins = 0U;
  float edges[kMaxBins + 1U]{};
  float ownColumn[kLayers][kMaxBins][kSpans]{};
  float ownRow[kLayers][kMaxBins][kSpans]{};
  float chanceColumn[kLayers][kSpans]{};
  float chanceRow[kLayers][kSpans]{};
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_SHAPESHARES_H_
