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

/// \file ITSROFSelection.h
/// \brief The ITS's readout frame selection, O2's own (`fastMultConfig`,
///        `FastMultEst`), over the three devices that track:
///
/// - the cluster converter, which reads O2's clusters, selects the frames
///   by the clusters' multiplicity, the random rejection and the triggers
///   (`selectFrames`) and ships the frames' mask (`kFrameMask`);
/// - the measurement producer, which writes no measurement of a frame the
///   mask disables (`FrameMask::enabled`);
/// - the tracker, which applies the vertex-multiplicity cut to the mask,
///   drops the measurements of the frames it disables
///   (`ITSFrameSelection`) and publishes the clock layer's frames
///   (`kSelectedROFs`), the adaptor's IRFRAMES.
///
/// The producer and the tracker are the detectors' own (Workflow/
/// MeasurementProducer.h, Workflow/TrackerDevice.h): they take it as the
/// ITS's `FrameSelectionOf` (Workflow/FrameSelection.h).
///
/// Without a selection configured (`frameSelection()`) none of them takes
/// part: the converter ships no mask and the tracker an empty list, every
/// frame selected.

#ifndef GTRK_WORKFLOW_ITSROFSELECTION_H_
#define GTRK_WORKFLOW_ITSROFSELECTION_H_

#include "ITSMFTInput.h"

#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Workflow/FrameSelection.h"
#include "Framework/InputSpec.h"
#include "Framework/OutputSpec.h"

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <vector>

namespace o2::framework
{
class ProcessingContext;
}

namespace o2::gtrk::workflow
{

/// Whether the ITS selects frames: a cut of `fastMultConfig` or a trigger
/// (`GTRKRecoITS.selectWithTriggers`). Read where the devices are defined and
/// in each device, from the same configuration.
[[nodiscard]] bool frameSelection();

/// The frames' mask as the converter ships it: each layer's frame count, then
/// one byte per frame, layer by layer (O2's `ROFMaskTable` layout).
struct FrameMask {
  std::vector<std::uint32_t> counts; ///< per layer
  std::vector<std::uint8_t> mask;    ///< 1 for an enabled frame

  /// Whether frame \a frame of layer \a layer is enabled; throws
  /// std::out_of_range for a frame the mask does not have.
  [[nodiscard]] bool enabled(int layer, std::uint32_t frame) const;
};

/// The converter's inputs and outputs for the selection: the trigger stream,
/// if one is configured, and the mask.
void addFrameSelectionInputs(std::vector<o2::framework::InputSpec>& inputs);
void addFrameMaskOutputs(std::vector<o2::framework::OutputSpec>& outputs);
/// The mask's inputs, the producer's and the tracker's.
void addFrameMaskInputs(std::vector<o2::framework::InputSpec>& inputs);

/// The converter's selection of \a input's frames, \a timing each layer's.
[[nodiscard]] FrameMask selectFrames(o2::framework::ProcessingContext& pc, const ITSMFTInput& input, std::span<const FrameTiming> timing);
/// The converter's mask into its outputs.
void shipFrameMask(o2::framework::ProcessingContext& pc, const FrameMask& mask);
/// The mask from the inputs of `addFrameMaskInputs`.
[[nodiscard]] FrameMask receiveFrameMask(o2::framework::ProcessingContext& pc);

/// The ITS's frame selection as the producer and the tracker take it
/// (Workflow/FrameSelection.h); the tracker's part, the vertex-multiplicity
/// cut on the converter's mask and the frames it leaves, is the object.
class ITSFrameSelection
{
 public:
  [[nodiscard]] static bool configured() { return frameSelection(); }
  static void addMaskInputs(std::vector<o2::framework::InputSpec>& inputs) { addFrameMaskInputs(inputs); }
  [[nodiscard]] static FrameMask receiveMask(o2::framework::ProcessingContext& pc) { return receiveFrameMask(pc); }

  ITSFrameSelection(const FrameMask& mask, std::span<const FrameTiming> timing);
  ~ITSFrameSelection();
  ITSFrameSelection(const ITSFrameSelection&) = delete;
  ITSFrameSelection& operator=(const ITSFrameSelection&) = delete;
  ITSFrameSelection(ITSFrameSelection&&) = delete;
  ITSFrameSelection& operator=(ITSFrameSelection&&) = delete;

  /// Disables the frames overlapping a vertex that fails the cut; the
  /// vertices stay in the output, as in O2's ITS. Returns whether the cut is
  /// configured.
  bool cut(std::span<const VertexSeed> seeds, std::span<const std::uint32_t> contributors);
  /// The clock layer's frames, 1 for a selected one.
  [[nodiscard]] std::vector<std::uint8_t> clockMask() const;

  /// The measurements of \a source whose frames are still enabled, in their
  /// bucket order, \a layerOf a measurement's layer from its surface; none
  /// when every one is.
  template <class LayerOf>
  [[nodiscard]] std::optional<SelectedMeasurements> filter(MeasurementView source, std::span<const std::uint32_t> offsets, const LayerOf& layerOf) const
  {
    std::vector<std::uint8_t> keep(offsets.back());
    bool all = true;
    for (std::uint32_t row = 0; row < offsets.back(); ++row) {
      keep[row] = enabled(layerOf(source.surface[row]), source.timePS[row], source.widthPS[row]) ? 1U : 0U;
      all = all && keep[row] != 0U;
    }
    if (all) {
      return std::nullopt;
    }
    return select(source, offsets, keep);
  }

 private:
  /// Whether the frame of a measurement on \a layer, its time box
  /// [\a timePS, + \a widthPS), is enabled.
  [[nodiscard]] bool enabled(int layer, std::uint32_t timePS, std::uint32_t widthPS) const;
  [[nodiscard]] static SelectedMeasurements select(MeasurementView source, std::span<const std::uint32_t> offsets, std::span<const std::uint8_t> keep);

  struct State;
  std::unique_ptr<State> mState;
};

template <>
struct FrameSelectionOf<detectors::ITS> {
  using type = ITSFrameSelection; // NOLINT(readability-identifier-naming): a type trait's member
};

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_ITSROFSELECTION_H_
