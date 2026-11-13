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

/// \file FrameSelection.h
/// \brief A detector's selection of readout frames in the devices that track:
///        `FrameSelectionOf<Det>::type`, none by default (`NoFrameSelection`);
///        the ITS's in WorkflowITSMFT/ITSROFSelection.h.
///
/// What the measurement producer and the tracker ask of it:
///
/// - `configured()`, whether the run selects frames at all; nothing below is
///   used without;
/// - `addMaskInputs(inputs)`, the frames' mask among a device's inputs, and
///   `receiveMask(pc)`, the timeframe's; `Mask::enabled(layer, frame)`, the
///   producer's test of a cluster's frame;
/// - the tracker's selection, from the mask and the layers' frames:
///   `cut(vertices, contributors)`, the cut by the vertices' contributors,
///   true if configured; `filter(measurements, offsets, layerOf)`, the
///   measurements of the frames left, none if every one is; `clockMask()`,
///   the frames selected on the output clock, the adaptor's IRFRAMES.

#ifndef GTRK_WORKFLOW_FRAMESELECTION_H_
#define GTRK_WORKFLOW_FRAMESELECTION_H_

#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/PairCut.h"

#include "Framework/InputSpec.h"

#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace o2::framework
{
class ProcessingContext;
}

namespace o2::gtrk::workflow
{

/// The measurements of the selected frames, copied, in their bucket order.
struct SelectedMeasurements {
  MeasurementStore columns;
  std::vector<std::uint32_t> offsets;
  std::vector<std::uint32_t> originalRows; ///< selected row -> the shipped measurement's row
};

/// Every frame selected.
class NoFrameSelection
{
 public:
  struct Mask {
    [[nodiscard]] static bool enabled(int /*layer*/, std::uint32_t /*frame*/) { return true; }
  };

  [[nodiscard]] static bool configured() { return false; }
  static void addMaskInputs(std::vector<o2::framework::InputSpec>& /*inputs*/) {}
  [[nodiscard]] static Mask receiveMask(o2::framework::ProcessingContext& /*pc*/) { return {}; }

  NoFrameSelection(const Mask& /*mask*/, std::span<const FrameTiming> /*timing*/) {}

  [[nodiscard]] static bool cut(std::span<const VertexSeed> /*seeds*/, std::span<const std::uint32_t> /*contributors*/) { return false; }
  template <class LayerOf>
  [[nodiscard]] static std::optional<SelectedMeasurements> filter(MeasurementView /*source*/, std::span<const std::uint32_t> /*offsets*/, const LayerOf& /*layerOf*/)
  {
    return std::nullopt;
  }
  [[nodiscard]] static std::vector<std::uint8_t> clockMask() { return {}; }
};

/// Detector \a Det's frame selection.
template <DetectorLayout Det>
struct FrameSelectionOf {
  using type = NoFrameSelection; // NOLINT(readability-identifier-naming): a type trait's member
};

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_FRAMESELECTION_H_
