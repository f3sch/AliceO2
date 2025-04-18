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
///

#ifndef ITSTRACKINGGPU_TRACKERTRAITSGPU_H_
#define ITSTRACKINGGPU_TRACKERTRAITSGPU_H_

#include "ITStracking/Configuration.h"
#include "ITStracking/Definitions.h"
#include "ITStracking/TrackerTraits.h"
#include "ITStrackingGPU/TimeFrameGPU.h"
#include "Framework/Logger.h"

namespace o2
{
namespace its
{

template <int nLayers = 7>
class TrackerTraitsGPU final : public TrackerTraits
{
 public:
  TrackerTraitsGPU() = default;
  ~TrackerTraitsGPU() override = default;

  void adoptTimeFrame(TimeFrame* tf) final;
  void initialiseTimeFrame(const int iteration) final;
  void setBz(float) final;

  void computeLayerTracklets(const int iteration, int, int) final { LOGP(fatal, "computeLayerTracklers must never be called from Hybrid traits!"); };
  void computeLayerCells(const int iteration) final { LOGP(fatal, "computeLayerCells must never be called from Hybrid traits!"); };
  void findCellsNeighbours(const int iteration) final { LOGP(fatal, "findCellsNeighbours must never be called from Hybrid traits!"); };
  void findRoads(const int iteration) final { LOGP(fatal, "findRoads must never be called from Hybrid traits!"); };
  void extendTracks(const int iteration) final { LOGP(fatal, "extendTracks must never be called from Hybrid traits!"); };
  void findShortPrimaries() final { LOGP(fatal, "findShortPrimaries must never be called from Hybrid traits!"); };

  void initialiseTimeFrameHybrid(const int iteration) final { initialiseTimeFrame(iteration); };
  void computeTrackletsHybrid(const int iteration, int, int) final;
  void computeCellsHybrid(const int iteration) final;
  void findCellsNeighboursHybrid(const int iteration) final;
  void findRoadsHybrid(const int iteration) final;

  // TimeFrameGPU information forwarding
  int getTFNumberOfClusters() const override;
  int getTFNumberOfTracklets() const override;
  int getTFNumberOfCells() const override;

 private:
  IndexTableUtils* mDeviceIndexTableUtils;
  gpu::TimeFrameGPU<7>* mTimeFrameGPU;
};

template <int nLayers>
inline void TrackerTraitsGPU<nLayers>::adoptTimeFrame(TimeFrame* tf)
{
  mTimeFrameGPU = static_cast<gpu::TimeFrameGPU<nLayers>*>(tf);
  mTimeFrame = static_cast<TimeFrame*>(tf);
}
} // namespace its
} // namespace o2

#endif