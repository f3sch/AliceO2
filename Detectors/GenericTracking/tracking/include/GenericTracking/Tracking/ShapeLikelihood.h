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

/// \file ShapeLikelihood.h
/// \brief How much more likely a row's cluster spans are for the track's own
///        crossing than for any row of its layer (`ShapeShares`).

#ifndef GTRK_TRACKING_SHAPELIKELIHOOD_H_
#define GTRK_TRACKING_SHAPELIKELIHOOD_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/ShapeShares.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include "GPUCommonMath.h"

#include <cstdint>

namespace o2::gtrk
{

/// log P(spans | the track's own crossing) - log P(spans | a row of the
/// layer) for a row of \a spans (`MeasurementT::spans`) on \a layer crossed
/// at the slopes \a columnSlope and \a rowSlope, the column's span and the
/// row's taken as independent; a side the calibration does not know, a row
/// without a pattern or a frame without shares adds nothing.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() Scalar logSizeRatio(const ShapeShares<Det>& shares, LayerID layer, std::uint32_t spans, Scalar columnSlope, Scalar rowSlope)
{
  constexpr std::uint32_t kLast = ShapeShares<Det>::kSpans - 1U;
  const std::uint32_t row = o2::gpu::CAMath::Min(measurementRowSpan(spans), kLast);
  const std::uint32_t column = o2::gpu::CAMath::Min(measurementColumnSpan(spans), kLast);
  if (shares.bins == 0U || row == 0U || column == 0U) {
    return 0.F;
  }
  const auto binOf = [&](Scalar slope) -> std::int32_t {
    if (!(slope >= shares.edges[0]) || !(slope < shares.edges[shares.bins])) {
      return -1;
    }
    std::uint32_t bin = 0U;
    while (slope >= shares.edges[bin + 1U]) {
      ++bin;
    }
    return static_cast<std::int32_t>(bin);
  };
  Scalar ratio = 0.F;
  const std::int32_t columnBin = binOf(columnSlope);
  if (columnBin >= 0) {
    const float own = shares.ownColumn[layer][columnBin][column];
    ratio += own == own ? own - shares.chanceColumn[layer][column] : 0.F;
  }
  const std::int32_t rowBin = binOf(rowSlope);
  if (rowBin >= 0) {
    const float own = shares.ownRow[layer][rowBin][row];
    ratio += own == own ? own - shares.chanceRow[layer][row] : 0.F;
  }
  return ratio;
}

/// `logSizeRatio` for row \a row of \a columns on \a layer, crossed along
/// \a direction (any length): its slopes along the columns and the rows over
/// its length along the normal -- on the barrel the band's (z and r-phi, the
/// radial normal at the row), on a disc the chip's own axes
/// (ClusterShapeBand.h).
template <DetectorLayout Det, typename Columns>
[[nodiscard]] GPUhdi() Scalar logSizeRatioAlong(const ShapeShares<Det>& shares, const TrackingGeometry<Det>& geometry, const Columns& columns, std::uint32_t row, LayerID layer, const Scalar (&direction)[3])
{
  if (shares.bins == 0U) {
    return 0.F;
  }
  Scalar alongNormal = 0.F;
  Scalar alongColumns = 0.F;
  Scalar alongRows = 0.F;
  if (Detector<Det>::isBarrelLayer(layer)) {
    const Scalar x = columns.x[row];
    const Scalar y = columns.y[row];
    const Scalar r = columns.r[row];
    if (!(r > 0.F)) {
      return 0.F;
    }
    alongNormal = o2::gpu::CAMath::Abs((direction[0] * x) + (direction[1] * y)) / r;
    alongColumns = o2::gpu::CAMath::Abs(direction[2]);
    alongRows = o2::gpu::CAMath::Abs((direction[1] * x) - (direction[0] * y)) / r;
  } else {
    const auto& sensor = geometry.surface(columns.surface[row]);
    const auto dot = [&](const Scalar(&axis)[3]) { return o2::gpu::CAMath::Abs((direction[0] * axis[0]) + (direction[1] * axis[1]) + (direction[2] * axis[2])); };
    alongNormal = dot(sensor.normal);
    alongColumns = dot(sensor.axis1);
    alongRows = dot(sensor.axis0);
  }
  return alongNormal > 0.F ? logSizeRatio(shares, layer, columns.spans[row], alongColumns / alongNormal, alongRows / alongNormal) : 0.F;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_SHAPELIKELIHOOD_H_
