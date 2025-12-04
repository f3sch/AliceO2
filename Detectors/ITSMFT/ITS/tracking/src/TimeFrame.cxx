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
/// \file TimeFrame.cxx
/// \brief
///

#include <numeric>
#include <utility>

#include "Framework/Logger.h"
#include "DetectorsRaw/HBFUtils.h"
#include "ITStracking/TimeFrame.h"
#include "ITStracking/Configuration.h"
#include "ITStracking/MathUtils.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "ITStracking/BoundedAllocator.h"

namespace
{
struct ClusterHelper {
  float phi;
  float r;
  int bin;
  int ind;
};
} // namespace

namespace o2::its
{

constexpr float DefClusErrorRow = o2::itsmft::SegmentationAlpide::PitchRow * 0.5;
constexpr float DefClusErrorCol = o2::itsmft::SegmentationAlpide::PitchCol * 0.5;
constexpr float DefClusError2Row = DefClusErrorRow * DefClusErrorRow;
constexpr float DefClusError2Col = DefClusErrorCol * DefClusErrorCol;

template <int NLayers>
void TimeFrame<NLayers>::loadROFrameData(gsl::span<const o2::itsmft::ROFRecord> rofs,
                                         gsl::span<const itsmft::CompClusterExt> clusters,
                                         gsl::span<const unsigned char>::iterator& pattIt,
                                         const itsmft::TopologyDictionary* dict,
                                         int layer,
                                         const dataformats::MCTruthContainer<MCCompLabel>* mcLabels)
{
  GeometryTGeo* geom = GeometryTGeo::Instance();
  geom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G));
  resetROFrameData(layer);
  prepareROFrameData(clusters, layer);

  // check for missing/empty/unset rofs
  // the code requires consistent monotonically increasing input without gaps
  const auto& timing = mROFOverlapTableView.getLayer(layer);
  if (timing.mNROFsTF != rofs.size()) {
    LOGP(fatal, "Received inconsistent number of rofs on layer:{} expected:{} received:{}", layer, timing.mNROFsTF, rofs.size());
  }

  for (int32_t iRof{0}; iRof < rofs.size(); ++iRof) {
    const auto& rof = rofs[iRof];
    for (int clusterId{rof.getFirstEntry()}; clusterId < rof.getFirstEntry() + rof.getNEntries(); ++clusterId) {
      const auto& c = clusters[clusterId];
      int lay = geom->getLayer(c.getSensorID());
      auto pattID = c.getPatternID();
      o2::math_utils::Point3D<float> locXYZ;
      float sigmaY2 = DefClusError2Row, sigmaZ2 = DefClusError2Col, sigmaYZ = 0; // Dummy COG errors (about half pixel size)
      unsigned int clusterSize{0};
      if (pattID != itsmft::CompCluster::InvalidPatternID) {
        sigmaY2 = dict->getErr2X(pattID);
        sigmaZ2 = dict->getErr2Z(pattID);
        if (!dict->isGroup(pattID)) {
          locXYZ = dict->getClusterCoordinates(c);
          clusterSize = dict->getNpixels(pattID);
        } else {
          o2::itsmft::ClusterPattern patt(pattIt);
          locXYZ = dict->getClusterCoordinates(c, patt);
          clusterSize = patt.getNPixels();
        }
      } else {
        o2::itsmft::ClusterPattern patt(pattIt);
        locXYZ = dict->getClusterCoordinates(c, patt, false);
        clusterSize = patt.getNPixels();
      }
      mClusterSize[layer][clusterId] = std::clamp(clusterSize, 0u, 255u);
      auto sensorID = c.getSensorID();
      // Inverse transformation to the local --> tracking
      auto trkXYZ = geom->getMatrixT2L(sensorID) ^ locXYZ;
      // Transformation to the local --> global
      auto gloXYZ = geom->getMatrixL2G(sensorID) * locXYZ;
      addTrackingFrameInfoToLayer(layer, gloXYZ.x(), gloXYZ.y(), gloXYZ.z(), trkXYZ.x(), geom->getSensorRefAlpha(sensorID),
                                  std::array<float, 2>{trkXYZ.y(), trkXYZ.z()},
                                  std::array<float, 3>{sigmaY2, sigmaYZ, sigmaZ2});
      /// Rotate to the global frame
      addClusterToLayer(layer, gloXYZ.x(), gloXYZ.y(), gloXYZ.z(), mUnsortedClusters[layer].size());
      addClusterExternalIndexToLayer(layer, clusterId);
    }
    mROFramesClusters[layer][iRof + 1] = mUnsortedClusters[layer].size(); // effectively calculating an exclusive sum
  }

  if (mcLabels != nullptr) {
    mClusterLabels[layer] = mcLabels;
  }
}

template <int NLayers>
void TimeFrame<NLayers>::resetROFrameData(int layer)
{
  deepVectorClear(mUnsortedClusters[layer], getMaybeFrameworkHostResource());
  deepVectorClear(mTrackingFrameInfo[layer], getMaybeFrameworkHostResource());
  deepVectorClear(mClusterExternalIndices[layer], mMemoryPool.get());
  clearResizeBoundedVector(mROFramesClusters[layer], mROFOverlapTableView.getLayer(layer).mNROFsTF + 1, getMaybeFrameworkHostResource());
}

template <int NLayers>
void TimeFrame<NLayers>::prepareROFrameData(gsl::span<const itsmft::CompClusterExt> clusters, int layer)
{
  GeometryTGeo* geom = GeometryTGeo::Instance();
  std::array<int, NLayers> clusterCountPerLayer{};
  for (const auto& clus : clusters) {
    auto lay = geom->getLayer(clus.getSensorID());
    if (lay != layer) {
      LOGP(fatal, "received layer from cluster {} while preparing data from {}!", lay, layer);
    }
    ++clusterCountPerLayer[lay];
  }
  mUnsortedClusters[layer].reserve(clusterCountPerLayer[layer]);
  mTrackingFrameInfo[layer].reserve(clusterCountPerLayer[layer]);
  mClusterExternalIndices[layer].reserve(clusterCountPerLayer[layer]);
  clearResizeBoundedVector(mClusterSize[layer], clusterCountPerLayer[layer], mMemoryPool.get());
}

template <int NLayers>
void TimeFrame<NLayers>::prepareClusters(const TrackingParameters& trkParam, const int maxLayers)
{

  const int numBins{trkParam.PhiBins * trkParam.ZBins};
  const int stride{numBins + 1};
  bounded_vector<ClusterHelper> cHelper(mMemoryPool.get());
  bounded_vector<int> clsPerBin(numBins, 0, mMemoryPool.get());
  bounded_vector<int> lutPerBin(numBins, 0, mMemoryPool.get());
  for (int iLayer{0}, stopLayer = std::min(trkParam.NLayers, maxLayers); iLayer < stopLayer; ++iLayer) {
    for (int rof{0}; rof < getNrof(iLayer); ++rof) {
      // TODO how to deal with mult mask?
      // if ((int)mMultiplicityCutMask.size() == mNrof && !mMultiplicityCutMask[rof]) {
      //   continue;
      // }
      const auto& unsortedClusters{getUnsortedClustersOnLayer(rof, iLayer)};
      const int clustersNum{static_cast<int>(unsortedClusters.size())};
      auto* tableBase = mIndexTables[iLayer].data() + rof * stride;

      cHelper.resize(clustersNum);

      for (int iCluster{0}; iCluster < clustersNum; ++iCluster) {
        const Cluster& c = unsortedClusters[iCluster];
        ClusterHelper& h = cHelper[iCluster];

        float x = c.xCoordinate;
        float y = c.yCoordinate;
        if (hasMeanVertex()) {
          x -= mMeanVertex->getX();
          y -= mMeanVertex->getY();
        }
        const float z = c.zCoordinate;

        float phi = math_utils::computePhi(x, y);
        int zBin{mIndexTableUtils.getZBinIndex(iLayer, z)};
        if (zBin < 0 || zBin >= trkParam.ZBins) {
          zBin = std::clamp(zBin, 0, trkParam.ZBins - 1);
          mBogusClusters[iLayer]++;
        }
        int bin = mIndexTableUtils.getBinIndex(zBin, mIndexTableUtils.getPhiBinIndex(phi));
        h.phi = phi;
        h.r = math_utils::hypot(x, y);
        mMinR[iLayer] = o2::gpu::GPUCommonMath::Min(h.r, mMinR[iLayer]);
        mMaxR[iLayer] = o2::gpu::GPUCommonMath::Max(h.r, mMaxR[iLayer]);
        h.bin = bin;
        h.ind = clsPerBin[bin]++;
      }
      std::exclusive_scan(clsPerBin.begin(), clsPerBin.end(), lutPerBin.begin(), 0);

      auto clusters2beSorted{getClustersOnLayer(rof, iLayer)};
      for (int iCluster{0}; iCluster < clustersNum; ++iCluster) {
        const ClusterHelper& h = cHelper[iCluster];
        Cluster& c = clusters2beSorted[lutPerBin[h.bin] + h.ind];

        c = unsortedClusters[iCluster];
        c.phi = h.phi;
        c.radius = h.r;
        c.indexTableBinIndex = h.bin;
      }
      std::copy_n(lutPerBin.data(), clsPerBin.size(), tableBase);
      std::fill_n(tableBase + clsPerBin.size(), stride - clsPerBin.size(), clustersNum);

      std::fill(clsPerBin.begin(), clsPerBin.end(), 0);
      cHelper.clear();
    }
  }
}

template <int NLayers>
void TimeFrame<NLayers>::initialise(const RecoIteration& reco)
{
  const auto& trkParam = reco.params;
  if (reco.steps[RecoIterationSteps::kInitMemory]) {
    mBogusClusters.fill(0); // reset bogus counter for this TF
    clearResizeBoundedVector(mCells, trkParam.CellsPerRoad(), mMemoryPool.get());
    clearResizeBoundedVector(mCellsLookupTable, trkParam.CellsPerRoad() - 1, mMemoryPool.get());
    clearResizeBoundedVector(mCellsNeighbours, trkParam.CellsPerRoad() - 1, mMemoryPool.get());
    clearResizeBoundedVector(mCellsNeighboursLUT, trkParam.CellsPerRoad() - 1, mMemoryPool.get());
    clearResizeBoundedVector(mCellLabels, trkParam.CellsPerRoad(), mMemoryPool.get());
    clearResizeBoundedVector(mTracklets, trkParam.TrackletsPerRoad(), mMemoryPool.get());
    clearResizeBoundedVector(mTrackletLabels, trkParam.TrackletsPerRoad(), mMemoryPool.get());
    clearResizeBoundedVector(mTrackletsLookupTable, trkParam.TrackletsPerRoad(), mMemoryPool.get());
    mIndexTableUtils.setTrackingParameters(trkParam);
    for (int iLayer{0}; iLayer < trkParam.NLayers; ++iLayer) {
      clearResizeBoundedVector(mClusters[iLayer], mUnsortedClusters[iLayer].size(), getMaybeFrameworkHostResource());
      clearResizeBoundedVector(mUsedClusters[iLayer], mUnsortedClusters[iLayer].size(), getMaybeFrameworkHostResource());
      if (iLayer < (int)mCells.size()) {
        mTrackletsLookupTable[iLayer].resize(mClusters[iLayer + 1].size() + 1, 0);
      }
    }
    for (int iLayer{0}; iLayer < NLayers; ++iLayer) {
      clearResizeBoundedVector(mIndexTables[iLayer], getNrof(iLayer) * ((trkParam.ZBins * trkParam.PhiBins) + 1), getMaybeFrameworkHostResource());
    }
    for (int iLayer{0}; iLayer < trkParam.NLayers; ++iLayer) {
      if (trkParam.SystErrorY2[iLayer] > 0.f || trkParam.SystErrorZ2[iLayer] > 0.f) {
        for (auto& tfInfo : mTrackingFrameInfo[iLayer]) {
          /// Account for alignment systematics in the cluster covariance matrix
          tfInfo.covarianceTrackingFrame[0] += trkParam.SystErrorY2[iLayer];
          tfInfo.covarianceTrackingFrame[2] += trkParam.SystErrorZ2[iLayer];
        }
      }
    }
    mMinR.fill(std::numeric_limits<float>::max());
    mMaxR.fill(std::numeric_limits<float>::min());
  }

  if (reco.steps[RecoIterationSteps::kUpdateClusters]) {
    prepareClusters(trkParam, reco.params.NLayers);
  }

  // these change dynamically with the given tracking parameters
  mMSangles.fill(0.f);
  mPhiCuts.fill(0.f);
  mPositionResolution.fill(0.f);
  float minCurvature{o2::gpu::CAMath::Abs(mBz * o2::constants::math::B2C) / trkParam.TrackletMinPt};
  for (int iLayer{0}; iLayer < trkParam.NLayers; ++iLayer) {
    mMSangles[iLayer] = math_utils::MSangle(0.14f, trkParam.TrackletMinPt, trkParam.LayerxX0[iLayer]);
    mPositionResolution[iLayer] = o2::gpu::CAMath::Sqrt((0.5f * (trkParam.SystErrorZ2[iLayer] + trkParam.SystErrorY2[iLayer])) + math_utils::Sq(trkParam.LayerResolution[iLayer]));
    if (iLayer < trkParam.NLayers - 1) {
      const float r1 = trkParam.LayerRadii[iLayer];
      const float r2 = trkParam.LayerRadii[iLayer + 1];
      const float res1 = mPositionResolution[iLayer];
      const float res2 = mPositionResolution[iLayer + 1];
      const float cosTheta1half = o2::gpu::CAMath::Sqrt(1.f - math_utils::Sq(0.5f * r1 * minCurvature));
      const float cosTheta2half = o2::gpu::CAMath::Sqrt(1.f - math_utils::Sq(0.5f * r2 * minCurvature));
      float x = (r2 * cosTheta1half) - (r1 * cosTheta2half);
      float delta = o2::gpu::CAMath::Sqrt(1.f / (1.f - 0.25f * math_utils::Sq(x * minCurvature)) * (math_utils::Sq((0.25f * r1 * r2 * math_utils::Sq(minCurvature) / cosTheta2half) + cosTheta1half) * math_utils::Sq(res1) + math_utils::Sq((0.25f * r1 * r2 * math_utils::Sq(minCurvature) / cosTheta1half) + cosTheta2half) * math_utils::Sq(res2)));
      mPhiCuts[iLayer] = std::min(o2::gpu::CAMath::ASin(0.5f * x * minCurvature) + 2.f * mMSangles[iLayer] + delta, o2::constants::math::PI * 0.5f);
    }
  }

  if (static bool initOnce{false}; !initOnce) {
    initOnce = true;
    // initialise the rolling vertex once with large weights
    mRollingVertex.setX(0.f);
    mRollingVertex.setY(0.f);
    mRollingVertex.setZ(0.f);
    mRollingVertex.setCov(1e6, dataformats::VertexBase::kCovXX);
    mRollingVertex.setCov(1e6, dataformats::VertexBase::kCovYY);
    mRollingVertex.setCov(1e6, dataformats::VertexBase::kCovZZ);
  }
}

template <int NLayers>
void TimeFrame<NLayers>::setMeanVertex(const dataformats::MeanVertexObject* mv, float extraErr2)
{
  mMeanVertex = mv;
  float ex2 = mMeanVertex->getSigmaX2() + extraErr2, ey2 = mMeanVertex->getSigmaY2() + extraErr2, exy = mMeanVertex->getSigmaXY();
  float det = (ex2 * ey2) - (exy * exy);
  if (det < constants::Tolerance || ex2 < constants::Tolerance || ey2 < constants::Tolerance) {
    LOGP(fatal, "Singular matrix for mean vertex sxx={:+.4e} syy={:+.4e} sxy={:+.4e}", ex2, ey2, exy);
  }
  mMeanVertexXYInvErr[0] = ey2 / det;
  mMeanVertexXYInvErr[1] = -exy / det;
  mMeanVertexXYInvErr[2] = -ex2 / det;
}

template <int NLayers>
unsigned long TimeFrame<NLayers>::getArtefactsMemory() const
{
  unsigned long size{0};
  for (const auto& trkl : mTracklets) {
    size += sizeof(Tracklet) * trkl.size();
  }
  for (const auto& cells : mCells) {
    size += sizeof(CellSeedN) * cells.size();
  }
  for (const auto& cellsN : mCellsNeighbours) {
    size += sizeof(int) * cellsN.size();
  }
  return size + sizeof(Road<NLayers - 2>) * mRoads.size();
}

template <int NLayers>
void TimeFrame<NLayers>::printArtefactsMemory() const
{
  LOGP(info, "TimeFrame: Artefacts occupy {:.2f} MB", getArtefactsMemory() / constants::MB);
}

template <int NLayers>
void TimeFrame<NLayers>::fillPrimaryVerticesXandAlpha()
{
  deepVectorClear(mPValphaX);
  mPValphaX.reserve(mPrimaryVertices.size());
  for (auto& pv : mPrimaryVertices) {
    mPValphaX.emplace_back(std::array<float, 2>{o2::gpu::CAMath::Hypot(pv.getX(), pv.getY()), math_utils::computePhi(pv.getX(), pv.getY())});
  }
}

template <int NLayers>
void TimeFrame<NLayers>::setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
{
  mMemoryPool = std::move(pool);

  auto initVector = [&]<typename T>(bounded_vector<T>& vec, bool useExternal = false) {
    std::pmr::memory_resource* mr = (useExternal) ? mExtMemoryPool.get() : mMemoryPool.get();
    deepVectorClear(vec, mr);
  };

  auto initContainers = [&]<typename Container>(Container& container, bool useExternal = false) {
    for (auto& v : container) {
      initVector(v, useExternal);
    }
  };

  // these will only reside on the host for the cpu part
  initContainers(mClusterExternalIndices);
  initVector(mPrimaryVertices);
  initVector(mRoads);
  initContainers(mClusterSize);
  initVector(mPValphaX);
  initVector(mTracks);
  initContainers(mTracklets);
  initContainers(mCells);
  initContainers(mCellsNeighbours);
  initContainers(mCellsLookupTable);
  // MC info (we don't know if we have MC (yet))
  initContainers(mTrackletLabels);
  initContainers(mCellLabels);
  initVector(mTracksLabel);
  initVector(mPrimaryVerticesLabels);
  // these will use possibly an externally provided allocator
  initContainers(mClusters, hasFrameworkAllocator());
  initContainers(mUsedClusters, hasFrameworkAllocator());
  initContainers(mUnsortedClusters, hasFrameworkAllocator());
  initContainers(mIndexTables, hasFrameworkAllocator());
  initContainers(mTrackingFrameInfo, hasFrameworkAllocator());
  initContainers(mROFramesClusters, hasFrameworkAllocator());
}

template <int nLayers>
void TimeFrame<nLayers>::setFrameworkAllocator(ExternalAllocator* ext)
{
  mExternalAllocator = ext;
  mExtMemoryPool = std::make_shared<BoundedMemoryResource>(mExternalAllocator);
}

template <int NLayers>
void TimeFrame<NLayers>::wipe()
{
  deepVectorClear(mTracks);
  deepVectorClear(mTracklets);
  deepVectorClear(mCells);
  deepVectorClear(mRoads);
  deepVectorClear(mCellsNeighbours);
  deepVectorClear(mCellsLookupTable);
  deepVectorClear(mPrimaryVertices);
  deepVectorClear(mTrackletsLookupTable);
  deepVectorClear(mClusterExternalIndices);
  deepVectorClear(mClusterSize);
  deepVectorClear(mPValphaX);
  // if we use the external host allocator then the assumption is that we
  // don't clear the memory ourself
  if (!hasFrameworkAllocator()) {
    deepVectorClear(mClusters);
    deepVectorClear(mUsedClusters);
    deepVectorClear(mUnsortedClusters);
    deepVectorClear(mIndexTables);
    deepVectorClear(mTrackingFrameInfo);
    deepVectorClear(mROFramesClusters);
  }
  // only needed to clear if we have MC info
  if (hasMCinformation()) {
    deepVectorClear(mTrackletLabels);
    deepVectorClear(mCellLabels);
    deepVectorClear(mTracksLabel);
    deepVectorClear(mPrimaryVerticesLabels);
  }
}

template class TimeFrame<7>;

} // namespace o2::its
