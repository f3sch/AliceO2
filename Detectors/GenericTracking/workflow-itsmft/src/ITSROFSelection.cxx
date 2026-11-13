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

#include "ITSROFSelection.h"

#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Workflow/DataSpecs.h"
#include "GenericTracking/Workflow/RecoParam.h"
#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "DataFormatsITSMFT/PhysTrigger.h"
#include "DataFormatsTRD/TriggerRecord.h"
#include "Framework/DataAllocator.h"
#include "Framework/InputRecord.h"
#include "Framework/Logger.h"
#include "Framework/ProcessingContext.h"
#include "Framework/TimingInfo.h"
#include "ITStracking/FastMultEst.h"

#include <algorithm>
#include <array>
#include <format>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace o2::gtrk::workflow
{

namespace
{
using Estimator = o2::its::FastMultEst;
using Overlap = Estimator::ROFOverlapTableN;
using Mask = Estimator::ROFMaskTableN;
using Timing = o2::itsmft::tracking::ROFTimingLayer;
using Alpide = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>;
constexpr int kLayers = Estimator::NLayers;

int triggerType()
{
  const auto type = recoParam<detectors::ITS>().selectWithTriggers;
  if (type < 0 || type > 2) {
    throw std::invalid_argument("GTRKRecoITS.selectWithTriggers must be 0 (none), 1 (ITS physics) or 2 (TRD)");
  }
  return type;
}

/// O2's frame tables of \a counts frames per layer: the layer's length, its
/// delay and bias from AlpideParam, its margin.
Overlap overlapTable(std::span<const std::uint32_t> counts, std::span<const FrameTiming> timing)
{
  if (counts.size() != kLayers || timing.size() != kLayers) {
    throw std::runtime_error("ITS frame selection: a frame count and a frame timing per layer");
  }
  const auto& alpide = Alpide::Instance();
  const auto& margins = recoParam<detectors::ITS>().frameMarginBC;
  Overlap overlap;
  for (int layer = 0; layer < kLayers; ++layer) {
    overlap.defineLayer(layer, {static_cast<Timing::BCType>(counts[layer]), static_cast<Timing::BCType>(timing[layer].lengthBC), static_cast<Timing::BCType>(alpide.getROFDelayInBC(layer)), static_cast<Timing::BCType>(alpide.getROFBiasInBC(layer)), static_cast<Timing::BCType>(margins.empty() ? 0 : margins[layer])});
  }
  overlap.init();
  return overlap;
}

std::string binding(std::string_view name) { return inputBinding<detectors::ITS>(name); }
} // namespace

bool frameSelection() { return o2::its::FastMultEstConfig::Instance().isRequested() || triggerType() != 0; }

bool FrameMask::enabled(int layer, std::uint32_t frame) const
{
  const auto at = static_cast<std::size_t>(layer);
  if (at >= counts.size() || frame >= counts[at]) {
    throw std::out_of_range(std::format("ITS frame selection: frame {} of layer {} is not in the mask", frame, layer));
  }
  return mask[std::accumulate(counts.begin(), counts.begin() + layer, std::size_t{0}) + frame] != 0U;
}

void addFrameSelectionInputs(std::vector<o2::framework::InputSpec>& inputs)
{
  using o2::framework::Lifetime;
  if (triggerType() == 1) {
    inputs.emplace_back(binding("triggers"), "ITS", "PHYSTRIG", 0, Lifetime::Timeframe);
  } else if (triggerType() == 2) {
    inputs.emplace_back(binding("triggers"), "TRD", "TRKTRGRD", 0, Lifetime::Timeframe);
  }
}

void addFrameMaskOutputs(std::vector<o2::framework::OutputSpec>& outputs)
{
  using o2::framework::Lifetime;
  outputs.emplace_back(gDataOrigin, description<detectors::ITS>(kFrameMask), 0, Lifetime::Timeframe);
  outputs.emplace_back(gDataOrigin, description<detectors::ITS>(kFrameMask), 1, Lifetime::Timeframe);
}

void addFrameMaskInputs(std::vector<o2::framework::InputSpec>& inputs)
{
  using o2::framework::Lifetime;
  if (o2::its::FastMultEstConfig::Instance().isVtxMultCutRequested() && !recoParam<detectors::ITS>().vertexer) {
    throw std::invalid_argument("fastMultConfig vertex cuts require GTRKRecoITS.vertexer=true (upstream VertexSeed inputs have no contributor counts)");
  }
  inputs.emplace_back(binding("frameMaskCounts"), gDataOrigin, description<detectors::ITS>(kFrameMask), 0, Lifetime::Timeframe);
  inputs.emplace_back(binding("frameMask"), gDataOrigin, description<detectors::ITS>(kFrameMask), 1, Lifetime::Timeframe);
}

FrameMask selectFrames(o2::framework::ProcessingContext& pc, const ITSMFTInput& input, std::span<const FrameTiming> timing)
{
  const auto streams = input.streams();
  if (streams.size() != 1U && streams.size() != static_cast<std::size_t>(kLayers)) {
    throw std::runtime_error("ITS frame selection: one stream or one per layer");
  }
  const bool staggered = streams.size() != 1U;
  std::array<gsl::span<const o2::itsmft::CompClusterExt>, kLayers> clusters{};
  std::array<gsl::span<const o2::itsmft::ROFRecord>, kLayers> rofs{};
  FrameMask result;
  result.counts.resize(kLayers);
  for (std::size_t stream = 0; stream < streams.size(); ++stream) {
    clusters[stream] = {streams[stream].compact.data(), streams[stream].compact.size()};
    rofs[stream] = {streams[stream].rofs.data(), streams[stream].rofs.size()};
  }
  // Every layer's frames as its records lay them out: O2 writes one record
  // per frame of the timeframe.
  for (int layer = 0; layer < kLayers; ++layer) {
    result.counts[layer] = static_cast<std::uint32_t>(rofs[staggered ? layer : 0].size());
  }
  auto overlap = overlapTable(result.counts, timing);
  Mask mask{overlap};

  const auto firstOrbit = pc.services().get<o2::framework::TimingInfo>().firstTForbit;
  const o2::InteractionRecord firstIR{0, firstOrbit};
  std::vector<o2::itsmft::PhysTrigger> triggers;
  if (triggerType() == 1) {
    const auto received = pc.inputs().get<gsl::span<o2::itsmft::PhysTrigger>>(binding("triggers").c_str());
    triggers.assign(received.begin(), received.end());
  } else if (triggerType() == 2) {
    auto lastIR = firstIR;
    const auto received = pc.inputs().get<gsl::span<o2::trd::TriggerRecord>>(binding("triggers").c_str());
    for (const auto& trigger : received) {
      if (trigger.getBCData() >= lastIR && trigger.getNumberOfTracklets() != 0) {
        lastIR = trigger.getBCData();
        triggers.push_back({.ir = lastIR, .data = 0});
      }
    }
  }
  // The estimator indexes the selection layer's frames with a trigger's
  // frame: the triggers outside this input's frames are not handed to it.
  const auto& mult = o2::its::FastMultEstConfig::Instance();
  const int layer = mult.isMultCutRequested() ? std::clamp(mult.cutMultClusLayer, 0, kLayers - 1) : overlap.getView().getClock();
  const auto count = rofs[staggered ? layer : 0].size();
  std::erase_if(triggers, [&](const auto& trigger) {
    const auto bc = trigger.ir.differenceInBC(firstIR);
    return bc < 0 || bc > std::numeric_limits<Timing::BCType>::max() || static_cast<std::size_t>(overlap.getLayer(layer).getROF(static_cast<Timing::BCType>(bc))) >= count;
  });
  const auto selected = Estimator{}.selectROFs(rofs, clusters, triggers, firstOrbit, staggered, overlap.getView(), mask);
  LOGP(info, "ITS frame selection: {} of {} frames on selection layer {}", selected, count, layer);

  const auto view = mask.getView();
  result.mask.assign(view.mFlatMask, view.mFlatMask + mask.getFlatMaskSize());
  return result;
}

void shipFrameMask(o2::framework::ProcessingContext& pc, const FrameMask& mask)
{
  using o2::framework::Output;
  pc.outputs().snapshot(Output{gDataOrigin, description<detectors::ITS>(kFrameMask), 0}, mask.counts);
  pc.outputs().snapshot(Output{gDataOrigin, description<detectors::ITS>(kFrameMask), 1}, mask.mask);
}

FrameMask receiveFrameMask(o2::framework::ProcessingContext& pc)
{
  const auto counts = pc.inputs().get<gsl::span<std::uint32_t>>(binding("frameMaskCounts").c_str());
  const auto mask = pc.inputs().get<gsl::span<std::uint8_t>>(binding("frameMask").c_str());
  FrameMask result{.counts = {counts.begin(), counts.end()}, .mask = {mask.begin(), mask.end()}};
  if (result.counts.size() != static_cast<std::size_t>(kLayers) || result.mask.size() != std::accumulate(result.counts.begin(), result.counts.end(), std::size_t{0})) {
    throw std::runtime_error("ITS frame selection: the mask does not hold its frame counts");
  }
  return result;
}

struct ITSFrameSelection::State {
  Overlap overlap;
  Mask mask;
  std::vector<FrameTiming> timing;
  std::vector<int> margins;
};

ITSFrameSelection::ITSFrameSelection(const FrameMask& mask, std::span<const FrameTiming> timing)
  : mState(std::make_unique<State>())
{
  auto& state = *mState;
  state.overlap = overlapTable(mask.counts, timing);
  state.mask = Mask{state.overlap};
  state.timing.assign(timing.begin(), timing.end());
  const auto& margins = recoParam<detectors::ITS>().frameMarginBC;
  state.margins = margins.empty() ? std::vector<int>(kLayers, 0) : margins;
  std::size_t at = 0;
  for (int layer = 0; layer < kLayers; ++layer) {
    for (std::uint32_t frame = 0; frame < mask.counts[layer]; ++frame, ++at) {
      state.mask.setROFEnabled(layer, static_cast<int>(frame), mask.mask[at]);
    }
  }
}

ITSFrameSelection::~ITSFrameSelection() = default;

bool ITSFrameSelection::cut(std::span<const VertexSeed> seeds, std::span<const std::uint32_t> contributors)
{
  if (!o2::its::FastMultEstConfig::Instance().isVtxMultCutRequested()) {
    return false;
  }
  if (seeds.size() != contributors.size()) {
    throw std::runtime_error("ITS vertex selection: contributor counts do not match the vertices");
  }
  std::vector<o2::its::Vertex> vertices;
  vertices.reserve(seeds.size());
  for (std::size_t v = 0; v < seeds.size(); ++v) {
    auto& vertex = vertices.emplace_back();
    vertex.setNContributors(static_cast<ushort>(std::min<std::uint32_t>(contributors[v], std::numeric_limits<ushort>::max())));
    const auto start = bcOfPS(seeds[v].timePS);
    const auto end = bcOfPS(static_cast<std::int64_t>(seeds[v].timePS) + seeds[v].widthPS);
    vertex.setTimeStamp({static_cast<std::uint32_t>(start), static_cast<std::uint32_t>(end - start)});
  }
  Estimator{}.selectROFsWithVertices(vertices, mState->overlap.getView(), mState->mask);
  return true;
}

std::vector<std::uint8_t> ITSFrameSelection::clockMask() const
{
  const auto clock = mState->overlap.getView().getClock();
  std::vector<std::uint8_t> result(mState->overlap.getLayer(clock).mNROFsTF);
  const auto mask = mState->mask.getView();
  for (std::size_t r = 0; r < result.size(); ++r) {
    result[r] = mask.isROFEnabled(clock, static_cast<int>(r)) ? 1U : 0U;
  }
  return result;
}

bool ITSFrameSelection::enabled(int layer, std::uint32_t timePS, std::uint32_t widthPS) const
{
  // The producer's box (`measurementTimeBox`): the frame [start, start +
  // length) widened by the layer's margin, clamped at the timeframe's start
  // below, so the frame is read from the upper edge; both edges exact in BC.
  const auto& timing = mState->timing[static_cast<std::size_t>(layer)];
  const auto upperBC = bcOfPS(timePS) + bcOfPS(widthPS);
  const auto fromFirst = upperBC - mState->margins[static_cast<std::size_t>(layer)] - timing.lengthBC - timing.firstStartBC;
  if (fromFirst < 0 || fromFirst % timing.lengthBC != 0) {
    throw std::runtime_error(std::format("ITS vertex selection: a measurement's box on layer {} is not one of the layer's frames", layer));
  }
  return mState->mask.getView().isROFEnabled(layer, static_cast<int>(fromFirst / timing.lengthBC));
}

SelectedMeasurements ITSFrameSelection::select(MeasurementView source, std::span<const std::uint32_t> offsets, std::span<const std::uint8_t> keep)
{
  std::vector<std::uint32_t> rows;
  for (std::uint32_t row = 0; row < keep.size(); ++row) {
    if (keep[row] != 0U) {
      rows.push_back(row);
    }
  }
  SelectedMeasurements selected;
  resizeMeasurementColumns(selected.columns, rows.size());
  selected.originalRows = rows;
  for (std::size_t i = 0; i < rows.size(); ++i) {
    Measurement measurement;
    measurement.apply(source[rows[i]], [](auto& to, const auto& from) -> auto& { to = from; return to; });
    writeMeasurement(selected.columns, i, measurement);
  }
  selected.offsets.reserve(offsets.size());
  auto at = rows.begin();
  for (const auto offset : offsets) {
    at = std::lower_bound(at, rows.end(), offset);
    selected.offsets.push_back(static_cast<std::uint32_t>(at - rows.begin()));
  }
  return selected;
}

} // namespace o2::gtrk::workflow
