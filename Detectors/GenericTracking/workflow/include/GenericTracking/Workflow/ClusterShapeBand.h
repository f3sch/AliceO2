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

/// \file ClusterShapeBand.h
/// \brief A detector's calibrated cluster-shape band, a CCDB object
///        (`<DET>/Calib/ClusterShapeBand`): per cluster shape the slope
///        interval a measurement carries.
///
/// Made from the simulation: per
/// calibrated layer, the most rows a particle above the calibration's floor
/// leaves, and per column span the hull of the local slopes (ShapeCut.h) at
/// which such a particle leaves that many columns -- the calibration's band
/// inverted, one interval per shape. A cluster
/// without a pattern, or on a layer the band does not calibrate, gets every
/// slope; one of more rows than its layer's most, or of a span no slope
/// gives, none. The slope is the barrel's (the normal radial, the columns
/// along z), so a band calibrates barrel layers only.
///
/// Beside the band, how often such a particle leaves each span at each
/// slope: per layer of the detector's whole layout -- the barrel's, then the
/// discs' (`ShapeBandLayers::layers`) -- and bin of the slope (`slopeEdges`),
/// the share of the calibration's crossings leaving each column span at its
/// slope along the columns and each row span at its slope along the rows,
/// spans from `kFrequencySpans` - 1 on together; none (every share 0) on a
/// layer or in a bin with too few crossings to say. On the barrel the slopes
/// are the band's (along z and along r-phi, over the length along the radial
/// normal); on a disc the chip's own (along its columns, loc1, and its rows,
/// loc0, over the length along its normal). The follow's last-row trim reads
/// them (`follow::ShapeFrequencies`) against the timeframe's own rows' spans.
///
/// A band is its detector's sensors', so its layers are numbered as the
/// detector's own layout numbers them, whatever layout tracks the detector
/// (`ShapeBandLayers`): TRK's band serves TRKFT3 and TRKFT3TOF alike.

#ifndef GTRK_WORKFLOW_CLUSTERSHAPEBAND_H_
#define GTRK_WORKFLOW_CLUSTERSHAPEBAND_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/ShapeCut.h"
#include "GenericTracking/Tracking/ShapeShares.h"
#ifdef ENABLE_UPGRADES
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Detectors/TRKFT3TOF.h"
#endif

#include <Rtypes.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace o2::gtrk::workflow
{

struct ClusterShapeBand {
  /// A measurement's slope interval [first, second).
  using Interval = std::pair<Scalar, Scalar>;
  static constexpr Interval kEvery{0.F, kSlopeUnbounded};
  static constexpr Interval kNone{std::numeric_limits<Scalar>::infinity(), 0.F};
  /// The column spans a cluster can have, 0 to 255.
  static constexpr std::uint32_t kSpans = 256;

  std::vector<std::int32_t> maxRows; ///< per layer from the first, its most rows; 0 on a layer the band does not calibrate, as every one after
  std::vector<float> slopeLow;  ///< per layer and column span, at layer * kSpans + span: the interval's start
  std::vector<float> slopeHigh; ///< and its end

  /// The spans the frequencies tell apart: 1 to kFrequencySpans - 1, the
  /// last one and every longer one together (0: no pattern).
  static constexpr std::uint32_t kFrequencySpans = 16;
  std::vector<float> slopeEdges;      ///< the slope bins' edges, ascending, from 0
  std::vector<float> columnFrequency; ///< per layer of the layout, slope bin and column span, at (layer * bins + bin) * kFrequencySpans + span
  std::vector<float> rowFrequency;    ///< per layer, slope bin and row span, the same

  [[nodiscard]] std::uint32_t slopeBins() const { return slopeEdges.empty() ? 0U : static_cast<std::uint32_t>(slopeEdges.size()) - 1U; }
  /// The layers the span shares cover.
  [[nodiscard]] std::uint32_t shareLayers() const { return slopeBins() == 0U ? 0U : static_cast<std::uint32_t>(columnFrequency.size()) / (slopeBins() * kFrequencySpans); }

  /// Whether it is a band of a layout of \a barrelLayers barrel layers: of
  /// its first layers at most, every layer's spans there --
  /// and its span shares every one of the layout's \a layers.
  [[nodiscard]] bool fits(std::uint32_t barrelLayers, std::uint32_t layers) const
  {
    const auto frequencies = layers * slopeBins() * kFrequencySpans;
    return maxRows.size() <= barrelLayers && slopeLow.size() == maxRows.size() * kSpans && slopeHigh.size() == maxRows.size() * kSpans && slopeBins() != 0U && columnFrequency.size() == frequencies && rowFrequency.size() == frequencies;
  }

  /// The slope interval of a cluster of \a rowSpan rows and \a colSpan
  /// columns on layer \a layer.
  [[nodiscard]] Interval interval(std::uint32_t layer, std::uint8_t rowSpan, std::uint8_t colSpan) const
  {
    if (colSpan == 0U || layer >= maxRows.size() || maxRows[layer] == 0) {
      return kEvery;
    }
    if (std::cmp_greater(rowSpan, maxRows[layer])) {
      return kNone;
    }
    const auto at = (layer * kSpans) + colSpan;
    return {slopeLow[at], slopeHigh[at]};
  }

  /// One line: the layers it calibrates, with their most rows.
  [[nodiscard]] std::string describe() const
  {
    std::string layers;
    for (std::size_t layer = 0; layer < maxRows.size(); ++layer) {
      if (maxRows[layer] != 0) {
        layers += std::format("{}layer {} ({} rows)", layers.empty() ? "" : ", ", layer, maxRows[layer]);
      }
    }
    return "cluster-shape band: " + (layers.empty() ? std::string{"no layer"} : layers);
  }

  ClassDefNV(ClusterShapeBand, 2)
};

/// How \a Det's layers read its detector's band: layer `layer(l)` of a band
/// of `barrelLayers` layers at most; a layout's own numbering by default.
template <DetectorLayout Det>
struct ShapeBandLayers {
  static constexpr auto barrelLayers = static_cast<std::uint32_t>(Detector<Det>::nBarrelLayers);
  static constexpr auto layers = static_cast<std::uint32_t>(Detector<Det>::nLayers);
  [[nodiscard]] static constexpr LayerID layer(LayerID layer) { return layer; }
};

#ifdef ENABLE_UPGRADES
/// TRKFT3TOF reads TRK's band in TRKFT3's numbering; its TOF layers are none
/// of the band's (kInvalidLayerID, past every band layer: every slope).
template <>
struct ShapeBandLayers<detectors::TRKFT3TOF> {
  static constexpr auto barrelLayers = static_cast<std::uint32_t>(Detector<detectors::TRKFT3>::nBarrelLayers);
  static constexpr auto layers = static_cast<std::uint32_t>(Detector<detectors::TRKFT3>::nLayers);
  [[nodiscard]] static constexpr LayerID layer(LayerID layer) { return detectors::trkft3tof::toTRKFT3(layer); }
};
#endif

/// \a shares' calibration side from \a band, per layer of \a Det through
/// `ShapeBandLayers`: the logarithms of its span shares, NaN where it has
/// none; no calibration (0 bins) without a band. Throws std::runtime_error
/// for a band of more slope bins than the shares hold.
template <DetectorLayout Det>
void fillOwnShares(const ClusterShapeBand* band, ShapeShares<Det>& shares)
{
  using Shares = ShapeShares<Det>;
  static_assert(Shares::kSpans == ClusterShapeBand::kFrequencySpans);
  shares.bins = 0U;
  if (band == nullptr) {
    return;
  }
  const auto bins = band->slopeBins();
  if (bins > Shares::kMaxBins) {
    throw std::runtime_error(std::format("the cluster-shape band's {} slope bins exceed the {} the shares hold", bins, Shares::kMaxBins));
  }
  shares.bins = bins;
  for (std::uint32_t edge = 0; edge <= bins; ++edge) {
    shares.edges[edge] = band->slopeEdges[edge];
  }
  constexpr float kNone = std::numeric_limits<float>::quiet_NaN();
  for (std::uint32_t layer = 0; layer < Shares::kLayers; ++layer) {
    const std::uint32_t own = ShapeBandLayers<Det>::layer(static_cast<LayerID>(layer));
    for (std::uint32_t bin = 0; bin < bins; ++bin) {
      for (std::uint32_t span = 0; span < Shares::kSpans; ++span) {
        const auto at = (((own * bins) + bin) * Shares::kSpans) + span;
        const bool known = own < band->shareLayers();
        shares.ownColumn[layer][bin][span] = known && band->columnFrequency[at] > 0.F ? std::log(band->columnFrequency[at]) : kNone;
        shares.ownRow[layer][bin][span] = known && band->rowFrequency[at] > 0.F ? std::log(band->rowFrequency[at]) : kNone;
      }
    }
  }
}

/// \a shares' timeframe side from each layer's counts of its rows' column
/// and row spans (\a columns, \a rows: layer * kSpans + span, spans from
/// kSpans - 1 on together): the logarithms of their shares, each span counted
/// once more (no share 0).
template <DetectorLayout Det>
void fillChanceShares(std::span<const std::uint32_t> columns, std::span<const std::uint32_t> rows, ShapeShares<Det>& shares)
{
  using Shares = ShapeShares<Det>;
  for (std::uint32_t layer = 0; layer < Shares::kLayers; ++layer) {
    std::uint64_t columnTotal = Shares::kSpans;
    std::uint64_t rowTotal = Shares::kSpans;
    for (std::uint32_t span = 0; span < Shares::kSpans; ++span) {
      columnTotal += columns[(layer * Shares::kSpans) + span];
      rowTotal += rows[(layer * Shares::kSpans) + span];
    }
    for (std::uint32_t span = 0; span < Shares::kSpans; ++span) {
      shares.chanceColumn[layer][span] = static_cast<float>(std::log(static_cast<double>(columns[(layer * Shares::kSpans) + span] + 1U) / static_cast<double>(columnTotal)));
      shares.chanceRow[layer][span] = static_cast<float>(std::log(static_cast<double>(rows[(layer * Shares::kSpans) + span] + 1U) / static_cast<double>(rowTotal)));
    }
  }
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_CLUSTERSHAPEBAND_H_
