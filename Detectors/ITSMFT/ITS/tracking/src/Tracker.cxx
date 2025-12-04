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
/// \file Tracker.cxx
/// \brief
///

#include "ITStracking/Tracker.h"

#include "ITStracking/BoundedAllocator.h"
#include "ITStracking/Configuration.h"
#include "ITStracking/Constants.h"
#include "ITStracking/ROFLookupTables.h"
#include "ITStracking/TrackerTraits.h"
#include "ITStracking/TrackingConfigParam.h"

#include <cassert>
#include <format>
#include <cstdlib>
#include <string>

namespace o2::its
{
using o2::its::constants::GB;

template <int NLayers>
Tracker<NLayers>::Tracker(TrackerTraits<NLayers>* traits) : mTraits(traits)
{
  if (traits->isGPU()) {
    ITSGpuTrackingParamConfig::Instance().maybeOverride();
    ITSGpuTrackingParamConfig::Instance().printKeyValues(true, true);
  }
}

template <int NLayers>
void Tracker<NLayers>::clustersToTracks(const LogFunc& logger, const LogFunc& error)
{
  LogFunc evalLog = [](const std::string&) {};

  double total{0};
  mTraits->updateTrackingParameters(mRecoParams);

  int iteration{0}, iSlice{0}, iVertex{0}, maxNvertices{-1};
  auto handleException = [&](const auto& err) {
    LOGP(error, "Too much memory used during {} in slice {} in iteration {} iVtx={}: {:.2f} GB. Current limit is {:.2f} GB, check the detector status and/or the selections.",
         StateNames[mCurState], iSlice, iteration, iVertex, (double)mTimeFrame->getArtefactsMemory() / GB, (double)mRecoParams[iteration].params.MaxMemory / GB);
    if (typeid(err) != typeid(std::bad_alloc)) { // only print if the exceptions is different from what is expected
      LOGP(error, "Exception: {}", err.what());
    }
    if (mRecoParams[iteration].params.DropTFUponFailure) {
      mMemoryPool->print();
      mTimeFrame->wipe();
      ++mNumberOfDroppedTFs;
      error("...Dropping Timeframe...");
    } else {
      throw err;
    }
  };

  try {
    for (iteration = 0; iteration < (int)mRecoParams.size(); ++iteration) {
      const auto& reco = mRecoParams[iteration];
      mMemoryPool->setMaxMemory(reco.params.MaxMemory);
      // rebuilt rof -> vertex lookup table
      if (reco.steps[RecoIterationSteps::kUpdateVertexTable] && mTimeFrame->getROFVertexLookupTable().needsUpdate()) {
        mTimeFrame->updateROFVertexLookupTable();
      }
      if (reco.params.PerPrimaryVertexProcessing) { // find the largest span of seeding vertices found
        maxNvertices = mTimeFrame->getROFVertexLookupTableView().getMaxVerticesPerROF();
      }

      // FIXME:
      // if (iteration == 3 && mTrkParams[0].DoUPCIteration) {
      // UPC: mark all already used ROFs and
      // then only use those where no tracjs where found.
      // has to be the last iteration
      // }

      double timeTracklets{0.}, timeCells{0.}, timeCellSeeds{0.}, timeNeighbours{0.}, timeRoads{0.};
      int nTracklets{0}, nCells{0}, nCellSeeds{-static_cast<int>(mTimeFrame->getPrimaryVerticesNum())}, nNeighbours{0}, nTracks{-static_cast<int>(mTimeFrame->getNumberOfTracks())};
      iVertex = std::min(maxNvertices, 0);
      logger(std::format("==== ITS {} Tracking iteration {} summary ==== ({})", mTraits->getName(), iteration, reco.name));

      if (reco.steps[RecoIterationSteps::kRunTruthSeeding]) {
        evaluateTask(&Tracker::computeTruthSeeding, StateNames[mCurState = kTruthSeeding], iteration, evalLog);
        continue;
      }

      total += evaluateTask(&Tracker::initialiseTimeFrame, StateNames[mCurState = kTFInit], iteration, logger);
      do {
        for (iSlice = 0; iSlice < reco.params.NTimeSlices; ++iSlice) {
          if (reco.steps[RecoIterationSteps::kRunTrackleting]) {
            timeTracklets += evaluateTask(&Tracker::computeTracklets, StateNames[mCurState = kTrackleting], iteration, evalLog, iSlice, iVertex);
            nTracklets += mTraits->getTFNumberOfTracklets();
          }
          if (reco.steps[RecoIterationSteps::kRunCellFinding]) {
            timeCells += evaluateTask(&Tracker::computeCells, StateNames[mCurState = kCelling], iteration, evalLog);
            nCells += mTraits->getTFNumberOfCells();
          }
          if (reco.steps[RecoIterationSteps::kRunCellSeeding]) {
            timeCellSeeds += evaluateTask(&Tracker::findCellSeeds, StateNames[mCurState = kSeeding], iteration, evalLog);
            nCellSeeds += mTimeFrame->getPrimaryVerticesNum();
          }
          if (reco.steps[RecoIterationSteps::kRunCellNeighborFinding]) {
            timeNeighbours += evaluateTask(&Tracker::findCellsNeighbours, StateNames[mCurState = kNeighbouring], iteration, evalLog);
            nNeighbours += mTimeFrame->getNumberOfNeighbours();
          }
          if (reco.steps[RecoIterationSteps::kRunRoadFinding]) {
            timeRoads += evaluateTask(&Tracker::findRoads, StateNames[mCurState = kRoading], iteration, evalLog);
          }
        }
      } while (++iVertex < maxNvertices);
      if (reco.steps[RecoIterationSteps::kRunTrackleting]) {
        logger(std::format(" - Tracklet finding: {} tracklets found in {:.2f} ms", nTracklets, timeTracklets));
      }
      if (reco.steps[RecoIterationSteps::kRunCellFinding]) {
        logger(std::format(" - Cell finding: {} cells found in {:.2f} ms", nCells, timeCells));
      }
      if (reco.steps[RecoIterationSteps::kRunCellSeeding]) {
        logger(std::format(" - Cell seeding: {} seeds found in {:.2f} ms", nCellSeeds, timeCellSeeds));
      }
      if (reco.steps[RecoIterationSteps::kRunCellNeighborFinding]) {
        logger(std::format(" - Neighbours finding: {} neighbours found in {:.2f} ms", nNeighbours, timeNeighbours));
      }
      if (reco.steps[RecoIterationSteps::kRunRoadFinding]) {
        logger(std::format(" - Track finding: {} tracks found in {:.2f} ms", nTracks + mTimeFrame->getNumberOfTracks(), timeRoads));
      }
      total += timeTracklets + timeCells + timeNeighbours + timeRoads + timeCellSeeds;
      // FIXME:
      // if (mTraits->supportsExtendTracks() && mTrkParams[iteration].UseTrackFollower) {
      //   int nExtendedTracks{-mTimeFrame->mNExtendedTracks}, nExtendedClusters{-mTimeFrame->mNExtendedUsedClusters};
      //   auto timeExtending = evaluateTask(&Tracker::extendTracks, "Extending tracks", iteration, evalLog, iteration);
      //   total += timeExtending;
      //   logger(std::format(" - Extending Tracks: {} extended tracks using {} clusters found in {:.2f} ms", nExtendedTracks + mTimeFrame->mNExtendedTracks, nExtendedClusters + mTimeFrame->mNExtendedUsedClusters, timeExtending));
      // }
    }
    if (mTraits->supportsFindShortPrimaries() && mRecoParams[0].params.FindShortTracks) {
      auto nTracksB = mTimeFrame->getNumberOfTracks();
      total += evaluateTask(&Tracker::findShortPrimaries, "Short primaries finding", 0, logger);
      auto nTracksA = mTimeFrame->getNumberOfTracks();
      logger(std::format("  `-> found {} additional tracks", nTracksA - nTracksB));
    }
    if constexpr (constants::DoTimeBenchmarks) {
      logger(std::format("=== TimeFrame {} processing completed in: {:.2f} ms using {} thread(s) ===", mTimeFrameCounter, total, mTraits->getNThreads()));
    }
  } catch (const BoundedMemoryResource::MemoryLimitExceeded& err) {
    handleException(err);
    return;
  } catch (const std::bad_alloc& err) {
    handleException(err);
    return;
  } catch (...) {
    error("Uncaught exception, all bets are off...");
  }

  if (mTimeFrame->hasMCinformation()) {
    computeTracksMClabels();
  }
  rectifyClusterIndices();
  sortTracks();

  ++mTimeFrameCounter;
  mTotalTime += total;
}

template <int NLayers>
void Tracker<NLayers>::computeTracksMClabels()
{
  for (auto& track : mTimeFrame->getTracks()) {
    std::vector<std::pair<MCCompLabel, size_t>> occurrences;
    occurrences.clear();

    for (int iCluster = 0; iCluster < TrackITSExt::MaxClusters; ++iCluster) {
      const int index = track.getClusterIndex(iCluster);
      if (index == constants::UnusedIndex) {
        continue;
      }
      auto labels = mTimeFrame->getClusterLabels(iCluster, index);
      bool found{false};
      for (size_t iOcc{0}; iOcc < occurrences.size(); ++iOcc) {
        std::pair<o2::MCCompLabel, size_t>& occurrence = occurrences[iOcc];
        for (const auto& label : labels) {
          if (label == occurrence.first) {
            ++occurrence.second;
            found = true;
          }
        }
      }
      if (!found) {
        for (const auto& label : labels) {
          occurrences.emplace_back(label, 1);
        }
      }
    }
    std::sort(std::begin(occurrences), std::end(occurrences), [](auto e1, auto e2) {
      return e1.second > e2.second;
    });

    auto maxOccurrencesValue = occurrences[0].first;
    uint32_t pattern = track.getPattern();
    // set fake clusters pattern
    for (int ic{TrackITSExt::MaxClusters}; ic--;) {
      auto clid = track.getClusterIndex(ic);
      if (clid != constants::UnusedIndex) {
        auto labelsSpan = mTimeFrame->getClusterLabels(ic, clid);
        for (const auto& currentLabel : labelsSpan) {
          if (currentLabel == maxOccurrencesValue) {
            pattern |= 0x1 << (16 + ic); // set bit if correct
            break;
          }
        }
      }
    }
    track.setPattern(pattern);
    if (occurrences[0].second < track.getNumberOfClusters()) {
      maxOccurrencesValue.setFakeFlag();
    }
    mTimeFrame->getTracksLabel().emplace_back(maxOccurrencesValue);
  }
}

template <int NLayers>
void Tracker<NLayers>::rectifyClusterIndices()
{
  for (auto& track : mTimeFrame->getTracks()) {
    for (int iCluster = 0; iCluster < TrackITSExt::MaxClusters; ++iCluster) {
      const int index = track.getClusterIndex(iCluster);
      if (index != constants::UnusedIndex) {
        track.setExternalClusterIndex(iCluster, mTimeFrame->getClusterExternalIndex(iCluster, index));
      }
    }
  }
}

template <int NLayers>
void Tracker<NLayers>::sortTracks()
{
  auto& trks = mTimeFrame->getTracks();
  std::vector<size_t> indices(trks.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&trks](size_t i, size_t j) {
    const auto& a = trks[i];
    const auto& b = trks[j];
    const auto at = a.getTimeStamp();
    const auto bt = b.getTimeStamp();
    if (at.getTimeStamp() != bt.getTimeStamp()) { // sort first in time
      return at.getTimeStamp() < bt.getTimeStamp();
    }
    return a.isBetter(b, 1e9); // then sort tracks in quality
  });
  bounded_vector<TrackITSExt> sortedTrks(trks.get_allocator());
  sortedTrks.reserve(trks.size());
  for (size_t idx : indices) {
    sortedTrks.push_back(trks[idx]);
  }
  trks.swap(sortedTrks);
  if (mTimeFrame->hasMCinformation()) {
    auto& trksLabels = mTimeFrame->getTracksLabel();
    bounded_vector<MCCompLabel> sortedLabels(trksLabels.get_allocator());
    sortedLabels.reserve(trksLabels.size());
    for (size_t idx : indices) {
      sortedLabels.push_back(trksLabels[idx]);
    }
    trksLabels.swap(sortedLabels);
  }
}

template <int NLayers>
void Tracker<NLayers>::adoptTimeFrame(TimeFrame<NLayers>& tf)
{
  mTimeFrame = &tf;
  mTraits->adoptTimeFrame(&tf);
}

template <int NLayers>
void Tracker<NLayers>::printSummary() const
{
  auto avgTF = mTotalTime * 1.e-3 / ((mTimeFrameCounter > 0) ? (double)mTimeFrameCounter : -1.0);
  auto avgTFwithDropped = mTotalTime * 1.e-3 / (((mTimeFrameCounter + mNumberOfDroppedTFs) > 0) ? (double)(mTimeFrameCounter + mNumberOfDroppedTFs) : -1.0);
  LOGP(info, "Tracker summary: Processed {} TFs (dropped {}) in TOT={:.2f} s, AVG/TF={:.2f} ({:.2f}) s", mTimeFrameCounter, mNumberOfDroppedTFs, mTotalTime * 1.e-3, avgTF, avgTFwithDropped);
}

template class Tracker<7>;

} // namespace o2::its
