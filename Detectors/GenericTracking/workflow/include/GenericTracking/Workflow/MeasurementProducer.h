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

/// \file MeasurementProducer.h
/// \brief The measurement producer, every detector's: instantiated in the
///        detector's own unit (ITS.cxx, MFT.cxx, TRK.cxx), where its
///        `Conditions<Det>` is.

#ifndef GTRK_WORKFLOW_MEASUREMENTPRODUCER_H_
#define GTRK_WORKFLOW_MEASUREMENTPRODUCER_H_

#include "GenericTracking/Base/FrameTiming.h"
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Base/TimeEstPS.h"
#include "GenericTracking/Clustering/Cluster.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/MeasurementGridBuilder.h"
#include "GenericTracking/Workflow/BeamLine.h"
#include "GenericTracking/Workflow/ClusterShapeBand.h"
#include "GenericTracking/Workflow/Conditions.h"
#include "GenericTracking/Workflow/DataSpecs.h"
#include "GenericTracking/Workflow/FrameSelection.h"
#include "GenericTracking/Workflow/GeometryProvider.h"
#include "GenericTracking/Workflow/GridConfig.h"
#include "GenericTracking/Workflow/MaterialLUT.h"
#include "GenericTracking/Workflow/MeasurementProducerSpec.h"
#include "GenericTracking/Workflow/RecoParam.h"

#include "DataFormatsCalibration/MeanVertexObject.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Logger.h"
#include "Framework/Task.h"

#include <gsl/span>

#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace o2::gtrk::workflow
{

/// A cluster's time box: the readout frame, as the converter wrote it from the
/// frame's record, moved to where the frame starts by the delay and bias of
/// \a frames, the layer's, and widened by the layer's margin, every edge
/// rounded once from BC.
[[nodiscard]] inline TimeEstPS measurementTimeBox(const Cluster& cluster, LayerID layer, const RecoSettings& param, const FrameTiming& frames)
{
  const auto margin = param.frameMarginBC.empty() ? 0 : param.frameMarginBC[layer];
  const auto startBC = bcOfPS(cluster.timePS) + frames.firstStartBC;
  return frameBox(startBC, startBC + bcOfPS(cluster.widthPS), margin);
}

/// Turns the cluster array into one message part per Measurement column,
/// ordered by (layer, time block, phi cell, second cell, box start), plus the
/// prefix offsets for those buckets.
///
/// Clusters arrive ordered by their readout frame -- on a layer timed per hit
/// (LayerSpec.h, `hitTimedLayers`) by their own box's start, which they keep
/// as their converter wrote it. The counting sort is
/// stable, so the box starts inside every bucket remain ordered without
/// another sort. The
/// geometry transform runs once per cluster, from its sensor ID and local
/// position: the float global position, with its transverse radius, azimuth
/// and projected variances, is kept beside the bucket number until the scatter
/// pass writes both out. It also ships the grid's layout and the beam line it
/// was derived from (`MeasurementGeometry`, GTRK/<DET>MEASGEO), and subtracts
/// the beam line from every measurement's x and y.
template <DetectorLayout Det>
class MeasurementProducerTask : public o2::framework::Task
{
 public:
  using Selection = typename FrameSelectionOf<Det>::type;

  MeasurementProducerTask(std::shared_ptr<Conditions<Det, true>> conditions, bool testsShapes) : mConditions(std::move(conditions)), mTestsShapes(testsShapes) {}

  void init(o2::framework::InitContext& ic) final
  {
    checkRecoParam<Det>(recoParam<Det>());
    mConditions->init(ic);
  }

  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj) final
  {
    if (matcher == o2::framework::ConcreteDataMatcher(Conditions<Det>::origin(), "SHAPEBAND", 0)) {
      mShapeBand = static_cast<const ClusterShapeBand*>(obj);
      if (!mShapeBand->fits(ShapeBandLayers<Det>::barrelLayers, ShapeBandLayers<Det>::layers)) {
        throw std::runtime_error(std::format("{}'s cluster-shape band from the CCDB is not one of its detector's {} barrel layers with span shares for its {} layers", Det::name, ShapeBandLayers<Det>::barrelLayers, ShapeBandLayers<Det>::layers));
      }
      LOGP(info, "{}'s {} from the CCDB", Det::name, mShapeBand->describe());
      return;
    }
    mConditions->finaliseCCDB(matcher, obj);
  }

  void run(o2::framework::ProcessingContext& pc) final
  {
    // The run's conditions: the geometry a placed layout is placed from, and
    // the field the grid's cells come from through the windows.
    mConditions->checkUpdates(pc);
    if (mTestsShapes) {
      pc.inputs().get<ClusterShapeBand*>("shapeBand"); // delivers a new one through finaliseCCDB
    }
    if (!mGeometry) {
      mGeometry = makeGeometry<Det>(mConditions->measurementDepthsCM());
      mFrames = mConditions->frameTiming();
    }
    // The timeframe's beam line, its MeanVertex's where the detector has one,
    // and the grid's layout, derived from the configuration, it and the
    // material LUT (the queries' scattering sizes the cells): rebuilt when it
    // changes, whatever order the timeframes come in.
    auto beamLine = nominalBeamLine();
    if constexpr (Conditions<Det>::kMeanVertex) {
      beamLine = beamLineOf(*pc.inputs().get<o2::dataformats::MeanVertexObject*>("meanVertex"));
    }
    if (!mGrid || beamLine != mBeamLine) {
      mBeamLine = beamLine;
      LOGP(info, "the beam line, subtracted from every measurement's x and y: {}", describe(mBeamLine));
      const auto* lut = mConditions->materialLut();
      if (lut == nullptr) {
        throw std::runtime_error(std::format("{}: no material LUT from the CCDB ({})", Det::name, materialLutPath<Det>()));
      }
      // The trackers' link material, held only while the grid is derived.
      const auto material = std::make_unique<LinkMaterial<Det>>();
      fillRunLinkMaterial<Det>(*material, recoParam<Det>(), beamSpotConstraint(mBeamLine, recoParam<Det>()), *lut);
      mGrid.emplace(makeGridGeometry<Det>(recoParam<Det>(), mBeamLine, mFrames, mConditions->orbitsPerTimeframe(), *material));
      std::string cells;
      std::string frames;
      for (std::size_t layer = 0; layer < GridGeometry<Det>::kLayerCount; ++layer) {
        const auto& axes = mGrid->axes(static_cast<LayerID>(layer));
        cells += std::format("{}{}x{}", layer == 0U ? "" : " ", axes.nPhi, axes.nSecond);
        frames += std::format("{}{}@{}", layer == 0U ? "" : " ", mFrames[layer].lengthBC, mFrames[layer].firstStartBC);
      }
      LOGP(info,
           "the grid, derived from the iterations: {} blocks of {} from {} in {} parts of a ps (frames of length@first start in BC by layer: {}); "
           "phi x second cells by layer: {}; {} buckets",
           mGrid->blockCount(), mGrid->blockWidth(), mGrid->blockOrigin(), mGrid->partsPerPS(), frames, cells, mGrid->bucketCount());
    }
    const auto clusters = pc.inputs().get<gsl::span<Cluster>>("clusters");
    if (clusters.size() > std::numeric_limits<std::uint32_t>::max()) {
      LOGP(fatal, "{} measurements exceed the uint32 row-index range", clusters.size());
      return;
    }
    for (const auto& cluster : clusters) {
      if (cluster.chipID >= static_cast<SensorID>(Det::nSensors)) {
        LOGP(fatal, "chip {} is outside {}", cluster.chipID, Det::name);
        return;
      }
    }

    // The detector's frame selection (FrameSelection.h): no measurement of a
    // frame its mask disables; a measurement keeps its cluster's index. A
    // cluster's record is a whole number of frames from BC 0.
    std::vector<std::uint32_t> kept;
    const auto mask = Selection::receiveMask(pc);
    const bool selected = !mask.allEnabled();
    if (selected) {
      kept.reserve(clusters.size());
      for (std::uint32_t i = 0; i < clusters.size(); ++i) {
        const auto layer = mGeometry->sensor(clusters[i].chipID).layer;
        if (mask.enabled(layer, static_cast<std::uint32_t>(bcOfPS(clusters[i].timePS) / mFrames[layer].lengthBC))) {
          kept.push_back(i);
        }
      }
    }
    const std::size_t n = selected ? kept.size() : clusters.size();
    const auto clusterOf = [&](std::size_t row) { return selected ? kept[row] : static_cast<std::uint32_t>(row); };
    const auto& param = recoParam<Det>();
    const auto& grid = *mGrid;
    const auto systematic2 = static_cast<Scalar>(param.clusterSystematicSigmaCM * param.clusterSystematicSigmaCM);
    const auto localVariance = [systematic2](float sigma) { return (sigma * sigma) + systematic2; };
    const auto start = std::chrono::steady_clock::now();
    const std::int64_t durationPS = timeframeDurationPS(mConditions->orbitsPerTimeframe());

    auto offsets = pc.outputs().make<std::uint32_t>(o2::framework::Output{gDataOrigin, description<Det>(kMeasurements), kColGridOffsets}, static_cast<std::size_t>(grid.bucketCount()) + 1U);
    MeasurementGridBuilder builder{std::span<std::uint32_t>{offsets.data(), offsets.size()}};

    std::vector<std::uint32_t> buckets(n);
    std::vector<TimeEstPS> boxes(n);
    std::vector<std::array<Scalar, 8>> positions(n); // x, y, z, r, phi, varZ, varR, varRPhi
    // The frame's bounds, row by row as the columns get them (`frameBounds`).
    auto bounds = HaloBounds<Det>::none();
    // The slope intervals only with a band: without one every measurement takes kEvery.
    const auto* band = mShapeBand;
    std::vector<ClusterShapeBand::Interval> slopes(band != nullptr ? n : 0U);
    for (std::size_t i = 0; i < n; ++i) {
      const auto& cluster = clusters[clusterOf(i)];
      // The tracking takes every measurement as a point with a positive
      // covariance: none of its steps tests for one that is not.
      const Scalar varLoc0 = localVariance(cluster.sigmaXLocal);
      const Scalar varLoc1 = localVariance(cluster.sigmaYLocal);
      if (!std::isfinite(cluster.xLocal) || !std::isfinite(cluster.yLocal) || !std::isfinite(varLoc0) || !std::isfinite(varLoc1) || !(varLoc0 > 0.F) || !(varLoc1 > 0.F)) {
        LOGP(fatal, "a {} cluster on chip {} is not a measurement: local position ({}, {}), sigmas ({}, {})", Det::name, cluster.chipID, cluster.xLocal, cluster.yLocal, cluster.sigmaXLocal, cluster.sigmaYLocal);
        return;
      }
      const auto point = gridPoint<Det>(*mGeometry, cluster.chipID, static_cast<double>(cluster.xLocal), static_cast<double>(cluster.yLocal), mBeamLine.x, mBeamLine.y);
      const auto variances = pointVariances<Det>(*mGeometry, cluster.chipID, static_cast<double>(cluster.xLocal), point, varLoc0, 0.F, varLoc1);
      // Every frame's record where its layer's frames have theirs, every
      // lengthBC from BC 0: the boxes, the blocks and the slices' reach follow them.
      // A layer timed per hit keeps the cluster's own box.
      const auto& frames = mFrames[point.layer];
      if (Detector<Det>::isHitTimedLayer(point.layer)) {
        if (cluster.widthPS == 0U || static_cast<std::int64_t>(cluster.timePS) + cluster.widthPS > durationPS) {
          LOGP(fatal, "a {} cluster's time box [{}, {}) ps on layer {} is empty or outside the timeframe", Det::name, cluster.timePS, static_cast<std::int64_t>(cluster.timePS) + cluster.widthPS, point.layer);
          return;
        }
        boxes[i] = {.start = cluster.timePS, .width = cluster.widthPS};
      } else {
        if (bcOfPS(cluster.timePS) % frames.lengthBC != 0 || bcOfPS(cluster.widthPS) != frames.lengthBC) {
          LOGP(fatal, "a {} cluster's readout frame record [{}, {}) ps on layer {} is not one of the layer's frames of {} BC", Det::name, cluster.timePS, static_cast<std::int64_t>(cluster.timePS) + cluster.widthPS, point.layer, frames.lengthBC);
          return;
        }
        boxes[i] = measurementTimeBox(cluster, point.layer, param, frames);
      }
      buckets[i] = grid.bucketAt(point.layer, boxes[i].start, point.phi, point.second);
      positions[i] = {point.x, point.y, point.z, point.r, point.phi, variances.z, variances.r, variances.rPhi};
      bounds.add(static_cast<std::size_t>(point.layer), grid.frames(point.layer), boxes[i], positions[i][5], positions[i][6], positions[i][7], varLoc0, varLoc1);
      if (band != nullptr) {
        slopes[i] = band->interval(static_cast<std::size_t>(ShapeBandLayers<Det>::layer(point.layer)), cluster.rowSpan, cluster.colSpan);
      }
      builder.count(buckets[i]);
    }
    builder.prefixSum();

    auto columns = makeMeasurementColumns<gsl::span>([&]<class T>(std::type_identity<T>, MeasurementColumn column) { return pc.outputs().make<T>(o2::framework::Output{gDataOrigin, description<Det>(kMeasurements), static_cast<std::uint32_t>(column)}, n); });

    for (std::size_t i = 0; i < n; ++i) {
      const auto& cluster = clusters[clusterOf(i)];
      const auto at = builder.place(buckets[i]);
      const auto slope = band != nullptr ? slopes[i] : ClusterShapeBand::kEvery;
      writeMeasurement(columns, at,
                       Measurement{
                         .loc0 = cluster.xLocal,
                         .loc1 = cluster.yLocal,
                         .timePS = boxes[i].start,
                         .widthPS = boxes[i].width,
                         .varLoc0 = localVariance(cluster.sigmaXLocal),
                         .covLoc01 = 0.F,
                         .varLoc1 = localVariance(cluster.sigmaYLocal),
                         .clusterIndex = clusterOf(i),
                         .surface = cluster.chipID,
                         .x = positions[i][0],
                         .y = positions[i][1],
                         .z = positions[i][2],
                         .r = positions[i][3],
                         .phi = positions[i][4],
                         .varZ = positions[i][5],
                         .varR = positions[i][6],
                         .varRPhi = positions[i][7],
                         .slopeLow = slope.first,
                         .slopeHigh = slope.second,
                         .spans = measurementSpans(cluster.rowSpan, cluster.colSpan),
                       });
    }
    builder.finish();
    pc.outputs().snapshot(o2::framework::Output{gDataOrigin, description<Det>(kMeasurementGeometry), 0}, MeasurementGeometry<Det>{mBeamLine, grid, bounds});

    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
    const auto payloadBytes = n * kMeasurementRowBytes + offsets.size() * sizeof(std::uint32_t);
    LOGP(info,
         "{} measurements in {} blocks and {} buckets over {} layers, {} columns, "
         "{} MB, {:.3f} s, about the beam line {}",
         n, grid.blockCount(), grid.bucketCount(), Detector<Det>::nLayers, kNMeasurementColumns, payloadBytes / 1000000U, elapsed.count(), describe(mBeamLine));
  }

 private:
  std::unique_ptr<const Geometry<Det>> mGeometry;     ///< placed on the first timeframe, once the conditions are in
  BeamLine mBeamLine;                                 ///< the last timeframe's
  std::optional<GridGeometry<Det>> mGrid;             ///< the grid's layout, derived from the configuration and mBeamLine
  std::vector<FrameTiming> mFrames;                   ///< each layer's readout frames, from the conditions
  std::shared_ptr<Conditions<Det, true>> mConditions; ///< the run's field, geometry, frames and material LUT (Conditions.h)
  bool mTestsShapes = false;                          ///< whether a pass tests the cluster shapes, which takes the detector's band
  const ClusterShapeBand* mShapeBand = nullptr;       ///< the CCDB cache's, with a pass that tests the shapes
};

template <DetectorLayout Det>
o2::framework::DataProcessorSpec getMeasurementProducerSpec()
{
  using namespace o2::framework;
  std::vector<OutputSpec> outputs;
  outputs.reserve(kNMeasurementParts);
  for (std::uint32_t part = 0; part < kNMeasurementParts; ++part) {
    outputs.emplace_back(gDataOrigin, description<Det>(kMeasurements), part, Lifetime::Timeframe);
  }
  outputs.emplace_back(gDataOrigin, description<Det>(kMeasurementGeometry), 0, Lifetime::Timeframe);
  static_assert(std::is_trivially_copyable_v<MeasurementGeometry<Det>>, "shipped as it is");

  std::vector<InputSpec> inputs{InputSpec{"clusters", gDataOrigin, description<Det>(kClusters), 0, Lifetime::Timeframe}};
  auto conditions = std::make_shared<Conditions<Det, true>>(inputs);
  if constexpr (Conditions<Det>::kMeanVertex) {
    inputs.emplace_back("meanVertex", "GLO", "MEANVERTEX", 0, Lifetime::Condition, ccdbParamSpec("GLO/Calib/MeanVertex", {}, 1));
  }
  FrameSelectionOf<Det>::type::addMaskInputs(inputs);
  // A pass testing the shapes takes the detector's band.
  const bool shapes = testsShapes(recoParam<Det>());
  if (shapes) {
    const auto origin = Conditions<Det>::origin();
    inputs.emplace_back("shapeBand", origin, "SHAPEBAND", 0, Lifetime::Condition, ccdbParamSpec(origin.template as<std::string>() + "/Calib/ClusterShapeBand"));
  }

  return DataProcessorSpec{
    .name = deviceName<Det>("measurement-producer"),
    .inputs = inputs,
    .outputs = outputs,
    .algorithm = AlgorithmSpec{adaptFromTask<MeasurementProducerTask<Det>>(conditions, shapes)},
    .options = Options{},
  };
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_MEASUREMENTPRODUCER_H_
