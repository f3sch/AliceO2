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
/// \file TrackerTraits.cxx
/// \brief
///

#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_sort.h>

#include <algorithm>
#include <iterator>
#include <limits>
#include <ranges>
#include <utility>
#include <cmath>

#include "DetectorsRaw/HBFUtils.h"
#include "Steer/MCKinematicsReader.h"
#include "SimulationDataFormat/O2DatabasePDG.h"
#include "CommonConstants/MathConstants.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "DetectorsBase/Propagator.h"
#include "ITSMFTBase/DPLAlpideParam.h"
#include "GPUCommonMath.h"
#include "ITStracking/Cell.h"
#include "ITStracking/MathUtils.h"
#include "ITStracking/Constants.h"
#include "ITStracking/Seeding.h"
#include "ITStracking/TrackerTraits.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITStracking/IndexTableUtils.h"
#include "ITStracking/Tracklet.h"
#include "ReconstructionDataFormats/Track.h"

/// optimization output
// TODO should this go somewhere else
#define OPTIMISATION_TRACKLETS (1 << 0)
#define OPTIMISATION_CELLS (1 << 1)
#define OPTIMISATION_CELLSNEIGH (1 << 2)
// #define OPTIMISATION (OPTIMISATION_CELLS | OPTIMISATION_TRACKLETS | OPTIMISATION_CELLSNEIGH)
#define OPTIMISATION (0)
#define OPTIMISATION_ANY (OPTIMISATION != 0)
#define OPTIMISATION_SET(flag) ((OPTIMISATION & (flag)) != 0)
#define OPTIMISATION_NOT_SET(flag) ((OPTIMISATION & (flag)) == 0)
#define OPTIMISATION_DOWNSAMPLE 1 // down-sample fraction

namespace o2::its
{
namespace
{
utils::TreeStreamRedirector* sDBGOut{nullptr};
}

struct PassMode {
  using OnePass = std::integral_constant<int, 0>;
  using TwoPassCount = std::integral_constant<int, 1>;
  using TwoPassInsert = std::integral_constant<int, 2>;
};

template <int NLayers>
TrackerTraits<NLayers>::TrackerTraits()
{
#if OPTIMISATION_ANY
  sDBGOut = new utils::TreeStreamRedirector("its_debug.root");
#endif
}

template <int NLayers>
TrackerTraits<NLayers>::~TrackerTraits()
{
#if OPTIMISATION_ANY
  sDBGOut->Close();
  delete sDBGOut;
  sDBGOut = nullptr;
#endif
}

template <int NLayers>
void TrackerTraits<NLayers>::computeLayerTracklets(const int iteration, const int iSlice, const int iVertex)
{
  for (int iLayer = 0; iLayer < mRecoParams[iteration].params.TrackletsPerRoad(); ++iLayer) {
    mTimeFrame->getTracklets()[iLayer].clear();
    mTimeFrame->getTrackletsLabel(iLayer).clear();
    if (iLayer > 0) {
      std::fill(mTimeFrame->getTrackletsLookupTable()[iLayer - 1].begin(), mTimeFrame->getTrackletsLookupTable()[iLayer - 1].end(), 0);
    }
  }

  const Vertex diamondVert(mRecoParams[iteration].params.Diamond, mRecoParams[iteration].params.DiamondCov, 1, 1.f);
  gsl::span<const Vertex> diamondSpan(&diamondVert, 1);

  mTaskArena->execute([&] {
    auto forTracklets = [&](auto Tag, int iLayer, int pivotROF, int base, int& offset) -> int {
      // do we even have clusters for this rof on this layer
      const auto& layer0 = mTimeFrame->getClustersOnLayer(pivotROF, iLayer);
      if (layer0.empty()) {
        return 0;
      }

      // if (!mTimeFrame->mMultiplicityCutMask[iLayer][pivotROF]) {
      //   return 0;
      // }

      // get seen seeding vertices
      gsl::span<const Vertex> primaryVertices = mRecoParams[iteration].params.UseDiamond ? diamondSpan : mTimeFrame->getPrimaryVertices(iLayer, pivotROF);
      if (primaryVertices.empty()) {
        return 0;
      }

      const int startVtx = iVertex >= 0 ? iVertex : 0;
      const int endVtx = iVertex >= 0 ? o2::gpu::CAMath::Min(iVertex + 1, int(primaryVertices.size())) : int(primaryVertices.size());
      if (endVtx <= startVtx || (iVertex + 1) > primaryVertices.size()) {
        return 0;
      }

      // does this layer have any overlap with the next layer
      const auto& rofOverlap = mTimeFrame->getROFOverlapTableView().getOverlap(iLayer, iLayer + 1, pivotROF);
      if (!rofOverlap.getEntries()) {
        return 0;
      }

      int localCount = 0;
      const float meanDeltaR = mRecoParams[iteration].params.LayerRadii[iLayer + 1] - mRecoParams[iteration].params.LayerRadii[iLayer];
      const float meanDeltaR2 = math_utils::Sq(meanDeltaR);

      auto& tracklets = mTimeFrame->getTracklets()[iLayer];
      for (int iCluster = 0; iCluster < int(layer0.size()); ++iCluster) {
        const Cluster& currentCluster = layer0[iCluster];
        const int currentSortedIndex = mTimeFrame->getSortedIndex(pivotROF, iLayer, iCluster);
        if (iteration && mTimeFrame->isClusterUsed(iLayer, currentCluster.clusterId)) {
          continue;
        }

        const float inverseR0 = 1.f / currentCluster.radius;
        const float inverseR02 = math_utils::Sq(inverseR0);

        for (int iV = startVtx; iV < endVtx; ++iV) {
          const auto& pv = primaryVertices[iV];

          // vertex resolution
          const float pvX2 = pv.getSigmaX2() / float(pv.getNContributors());
          const float pvY2 = pv.getSigmaY2() / float(pv.getNContributors());
          const float pvR2 = pvX2 + pvY2;
          const float pvZ2 = pv.getSigmaZ2() / float(pv.getNContributors());
          // xyz-position resolution
          const float res2 = math_utils::Sq(mTimeFrame->getPositionResolution(iLayer));
          const float deltaZ = currentCluster.zCoordinate - pv.getZ();
          /// calculate the lookup window for the next layer in z&phi
          // the assumption is that the clusters are already in the vertex coordinate system
          // Note: this is only approximate can we account for the difference somehow?
          // phi-window:
          // phi=arctan((y_c - y_pv)/(x_c - x_pv)); in this case this is just determined by the cluster
          // Var[phi] = fac^2 * [ S*(s_x_c^2 + s_x_pv^2) + T*(s_y_c^2+s_y_pv^2)]
          // where fac = 1/(S+T), S=(y_c - y_pv)^2, T=(x_c - x_pv)^2
          // additionally, adding the phi cuts due to the stave inclination which also accounts MS
          // NOTE: the cluster x&y positions are already in an approximate frame of the pv
          const float ms2 = math_utils::Sq(mTimeFrame->getMSangle(iLayer));
          const float phiExtrap = currentCluster.phi;
          const float phiS = math_utils::Sq(currentCluster.yCoordinate);
          const float phiT = math_utils::Sq(currentCluster.xCoordinate);
          const float phiFac = 1.f / (phiS + phiT);
          const float phiVar = (math_utils::Sq(phiFac) * (phiS * (res2 + pvX2) + phiT * (res2 + pvY2)));
          const float phiNSigma = mRecoParams[iteration].params.NSigmaCut * (o2::gpu::CAMath::Sqrt(phiVar) + mTimeFrame->getPhiCut(iLayer));
          // z-window:
          // z(r) = z_c + (z_pv - z_c)/(r_pv - r_c) * (r - r_c)
          // define S = r_pv - r_c, r_n = r_c + meanDeltaR
          // a := d z(r_n) / d z_c = 1 - (meanDeltaR / S)
          // b := d z(r_n) / d z_pv = meanDeltaR / S
          // d z(r_n) / d r_c = - (r_n - r_pv) * (z_c - z_pv) / S^2
          // d z(r_n) / d r_pv = - (r_c - r_n) * (z_c - z_pv) / S^2
          const float tanLambda = deltaZ * inverseR0;
          const float zExtrap = currentCluster.zCoordinate + (tanLambda * meanDeltaR);
          const float a = 1.f - (meanDeltaR * inverseR0);
          const float rA2 = math_utils::Sq(a);
          const float b = meanDeltaR * inverseR0;
          const float rB2 = math_utils::Sq(b);
          const float S = 1.f / inverseR0;
          const float dzdrc = -(meanDeltaR - S) * deltaZ * inverseR02; // = - (r_n - r_pv) * deltaZ / S^2
          const float dzdrpv = meanDeltaR * deltaZ * inverseR02;       // = - (r_c - r_n) * deltaZ / S^2
          const float zAtRmin = tanLambda * (mTimeFrame->getMinR(iLayer + 1) - currentCluster.radius) + currentCluster.zCoordinate;
          const float zAtRmax = tanLambda * (mTimeFrame->getMaxR(iLayer + 1) - currentCluster.radius) + currentCluster.zCoordinate;
          const float zWidth2 = (1.f / 12.f) * math_utils::Sq(zAtRmax - zAtRmin); // accounting for stave inclination
          const float zVar = (rA2 * res2) + (rB2 * pvZ2) + (dzdrc * dzdrc * res2) + (dzdrpv * dzdrpv * pvR2) + zWidth2;
          const float zNSigma = mRecoParams[iteration].params.NSigmaCut * o2::gpu::CAMath::Sqrt(zVar);

          const auto bins = getBinsRect(iteration, iLayer + 1, phiExtrap, phiNSigma, zExtrap, zNSigma);
          if (bins.x == 0 && bins.y == 0 && bins.z == 0 && bins.w == 0) {
            continue;
          }
          int phiBinsNum = bins.w - bins.y + 1;
          if (phiBinsNum < 0) { // handle wrap around
            phiBinsNum += mRecoParams[iteration].params.PhiBins;
          }

          for (int targetROF = rofOverlap.getFirstEntry(); targetROF < rofOverlap.getEntriesBound(); ++targetROF) {
            // if (!mTimeFrame->mMultiplicityCutMask[iLayer + 1][targetROF]) {
            //   continue;
            // }
            const auto& layer1 = mTimeFrame->getClustersOnLayer(targetROF, iLayer + 1);
            if (layer1.empty()) {
              continue;
            }
            const auto& targetIndexTable = mTimeFrame->getIndexTable(targetROF, iLayer + 1);
            const int zBinRange = (bins.z - bins.x) + 1;
            for (int iPhi = 0; iPhi < phiBinsNum; ++iPhi) {
              const int iPhiBin = (bins.y + iPhi) % mRecoParams[iteration].params.PhiBins;
              const int firstBinIdx = mTimeFrame->getIndexTableUtils().getBinIndex(bins.x, iPhiBin);
              const int maxBinIdx = firstBinIdx + zBinRange;
              const int firstRow = targetIndexTable[firstBinIdx];
              const int lastRow = targetIndexTable[maxBinIdx];
              for (int iNext = firstRow; iNext < lastRow; ++iNext) {
                if (iNext >= int(layer1.size())) {
                  break;
                }
                const Cluster& nextCluster = layer1[iNext];
                if (iteration && mTimeFrame->isClusterUsed(iLayer + 1, nextCluster.clusterId)) {
                  continue;
                }

                const float deltaPhi = o2::gpu::CAMath::Abs(o2::math_utils::toPMPi(currentCluster.phi - nextCluster.phi));
                const float deltaZ = o2::gpu::CAMath::Abs((tanLambda * (nextCluster.radius - currentCluster.radius)) + currentCluster.zCoordinate - nextCluster.zCoordinate);

                const bool accept = deltaZ < zNSigma && deltaPhi < phiNSigma;
                debugComputeLayerTracklets(iteration, iLayer, currentCluster, nextCluster, pv, zNSigma, phiNSigma, accept);

                if (accept) {
                  const float phi{o2::gpu::CAMath::ATan2(currentCluster.yCoordinate - nextCluster.yCoordinate, currentCluster.xCoordinate - nextCluster.xCoordinate)};
                  const float tanL = (currentCluster.zCoordinate - nextCluster.zCoordinate) / (currentCluster.radius - nextCluster.radius);

                  if constexpr (decltype(Tag)::value == PassMode::OnePass::value) {
                    tracklets.emplace_back(currentSortedIndex, mTimeFrame->getSortedIndex(targetROF, iLayer + 1, iNext), tanL, phi, pivotROF, targetROF);
                  } else if constexpr (decltype(Tag)::value == PassMode::TwoPassCount::value) {
                    ++localCount;
                  } else if constexpr (decltype(Tag)::value == PassMode::TwoPassInsert::value) {
                    const int idx = base + offset++;
                    tracklets[idx] = Tracklet(currentSortedIndex, mTimeFrame->getSortedIndex(targetROF, iLayer + 1, iNext), tanL, phi, pivotROF, targetROF);
                  }
                }
              }
            }
          }
        }
      }
      return localCount;
    };

    int dummy{0};
    if (mTaskArena->max_concurrency() <= 1) {
      for (int iLayer{0}; iLayer < mRecoParams[iteration].params.TrackletsPerRoad(); ++iLayer) {
        const auto& rofSlices = mTimeFrame->getROFTimeSliceTableView().getROFSlice(iLayer, iSlice);
        const auto& timeMask = mTimeFrame->getROFTimeSliceTableView().getMask(iLayer, iSlice);
        for (int pivotROF{rofSlices.getFirstEntry()}; pivotROF < rofSlices.getEntriesBound(); ++pivotROF) {
          if (!timeMask[pivotROF - rofSlices.getFirstEntry()]) {
            continue;
          }
          forTracklets(PassMode::OnePass{}, iLayer, pivotROF, 0, dummy);
        }
      }
    } else {
      tbb::parallel_for(0, mRecoParams[iteration].params.TrackletsPerRoad(), [&](const int iLayer) {
        const auto& rofSlices = mTimeFrame->getROFTimeSliceTableView().getROFSlice(iLayer, iSlice);
        const auto& timeMask = mTimeFrame->getROFTimeSliceTableView().getMask(iLayer, iSlice);
        bounded_vector<int> perROFCount(rofSlices.getEntries() + 1, mMemoryPool.get());
        tbb::parallel_for(rofSlices.getFirstEntry(), rofSlices.getEntriesBound(), [&](const int pivotROF) {
          if (!timeMask[pivotROF - rofSlices.getFirstEntry()]) {
            return;
          }
          perROFCount[pivotROF - rofSlices.getFirstEntry()] = forTracklets(PassMode::TwoPassCount{}, iLayer, pivotROF, 0, dummy);
        });
        std::exclusive_scan(perROFCount.begin(), perROFCount.end(), perROFCount.begin(), 0);
        const int nTracklets = perROFCount.back();
        mTimeFrame->getTracklets()[iLayer].resize(nTracklets);
        if (nTracklets == 0) {
          return;
        }
        tbb::parallel_for(rofSlices.getFirstEntry(), rofSlices.getEntriesBound(), [&](const int pivotROF) {
          if (!timeMask[pivotROF - rofSlices.getFirstEntry()]) {
            return;
          }
          int baseIdx = perROFCount[pivotROF - rofSlices.getFirstEntry()];
          if (baseIdx == perROFCount[pivotROF + 1 - rofSlices.getFirstEntry()]) {
            return;
          }
          int localIdx = 0;
          forTracklets(PassMode::TwoPassInsert{}, iLayer, pivotROF, baseIdx, localIdx);
        });
      });
    }

    tbb::parallel_for(0, mRecoParams[iteration].params.TrackletsPerRoad(), [&](const int iLayer) {
      /// Sort tracklets
      auto& trkl{mTimeFrame->getTracklets()[iLayer]};
      tbb::parallel_sort(trkl.begin(), trkl.end(), [](const Tracklet& a, const Tracklet& b) -> bool {
        if (a.firstClusterIndex != b.firstClusterIndex) {
          return a.firstClusterIndex < b.firstClusterIndex;
        }
        return a.secondClusterIndex < b.secondClusterIndex;
      });
      /// Remove duplicates
      trkl.erase(std::unique(trkl.begin(), trkl.end(), [](const Tracklet& a, const Tracklet& b) -> bool {
                   return a.firstClusterIndex == b.firstClusterIndex && a.secondClusterIndex == b.secondClusterIndex;
                 }),
                 trkl.end());
      if (iLayer > 0) { /// recalculate lut
        auto& lut{mTimeFrame->getTrackletsLookupTable()[iLayer - 1]};
        clearResizeBoundedVector(lut, mTimeFrame->getNumberOfClusters(iLayer) + 1, mMemoryPool.get(), 0);
        for (const auto& tkl : trkl) {
          lut[tkl.firstClusterIndex + 1]++;
        }
        std::inclusive_scan(lut.begin(), lut.end(), lut.begin());
      }
    });

    /// Create tracklets labels
    if ((mTimeFrame->hasMCinformation() && mRecoParams[iteration].params.createArtefactLabels) || OPTIMISATION_ANY) {
      tbb::parallel_for(0, mRecoParams[iteration].params.TrackletsPerRoad(), [&](const int iLayer) {
        for (auto& trk : mTimeFrame->getTracklets()[iLayer]) {
          MCCompLabel label;
          int currentId{mTimeFrame->getClusters()[iLayer][trk.firstClusterIndex].clusterId};
          int nextId{mTimeFrame->getClusters()[iLayer + 1][trk.secondClusterIndex].clusterId};
          for (const auto& lab1 : mTimeFrame->getClusterLabels(iLayer, currentId)) {
            for (const auto& lab2 : mTimeFrame->getClusterLabels(iLayer + 1, nextId)) {
              if (lab1 == lab2 && lab1.isValid()) {
                label = lab1;
                break;
              }
            }
            if (label.isValid()) {
              break;
            }
          }
          mTimeFrame->getTrackletsLabel(iLayer).emplace_back(label);
        }
      });
    }
  });
}

template <int NLayers>
void TrackerTraits<NLayers>::computeLayerCells(const int iteration)
{
  for (int iLayer = 0; iLayer < mRecoParams[iteration].params.CellsPerRoad(); ++iLayer) {
    deepVectorClear(mTimeFrame->getCells()[iLayer]);
    if (iLayer > 0) {
      deepVectorClear(mTimeFrame->getCellsLookupTable()[iLayer - 1]);
    }
    if ((mTimeFrame->hasMCinformation() && mRecoParams[iteration].params.createArtefactLabels) || OPTIMISATION_ANY) {
      deepVectorClear(mTimeFrame->getCellsLabel(iLayer));
    }
  }

  mTaskArena->execute([&] {
    auto forTrackletCells = [&](auto Tag, int iLayer, bounded_vector<CellSeedN>& layerCells, int iTracklet, int offset = 0) -> int {
      const Tracklet& currentTracklet{mTimeFrame->getTracklets()[iLayer][iTracklet]};
      const int nextLayerClusterIndex{currentTracklet.secondClusterIndex};
      const int nextLayerFirstTrackletIndex{mTimeFrame->getTrackletsLookupTable()[iLayer][nextLayerClusterIndex]};
      const int nextLayerLastTrackletIndex{mTimeFrame->getTrackletsLookupTable()[iLayer][nextLayerClusterIndex + 1]};

      // properties for the middle layer
      const float resMid2 = math_utils::Sq(mTimeFrame->getPositionResolution(iLayer + 1));
      const float msMid2 = math_utils::Sq(mTimeFrame->getMSangle(iLayer + 1));
      // the allowed tgl variance is entirely determined by the middle layer
      const float tglNSigma = o2::gpu::CAMath::Sqrt(resMid2 + msMid2) * mRecoParams[iteration].params.NSigmaCut;

      int foundCells{0};
      for (int iNextTracklet{nextLayerFirstTrackletIndex}; iNextTracklet < nextLayerLastTrackletIndex; ++iNextTracklet) {
        // ensure that the tracklets are sharing the middle cluster
        const Tracklet& nextTracklet{mTimeFrame->getTracklets()[iLayer + 1][iNextTracklet]};
        if (mTimeFrame->getTracklets()[iLayer + 1][iNextTracklet].firstClusterIndex != nextLayerClusterIndex) {
          break;
        }

        // need to ensure that the clusters in layer (iLayer) and (iLayer+2) are compatible
        if (!mTimeFrame->getROFOverlapTableView().isCompatible(iLayer, currentTracklet.rof[0], iLayer + 2, nextTracklet.rof[1])) {
          continue;
        }

        debugComputeLayerCells(iteration, iLayer, iTracklet, iNextTracklet);

        // calculate their compatibility in TgL
        const float deltaTanLambda = o2::gpu::CAMath::Abs(currentTracklet.tanLambda - nextTracklet.tanLambda);
        if (deltaTanLambda >= tglNSigma) {
          continue;
        }

        const auto& cls1 = mTimeFrame->getClusters()[iLayer][currentTracklet.firstClusterIndex];
        const auto& cls2 = mTimeFrame->getClusters()[iLayer + 1][nextTracklet.firstClusterIndex];
        const auto& cls3 = mTimeFrame->getClusters()[iLayer + 2][nextTracklet.secondClusterIndex];
        // curvature consistency cut arxiv:2401.16046 4.1.2
        // the transverse cluster positions are in the presumed beamspot frame
        const float k123 = math_utils::computeCurvature(cls1.xCoordinate, cls1.yCoordinate,
                                                        cls2.xCoordinate, cls2.yCoordinate,
                                                        cls3.xCoordinate, cls3.yCoordinate);
        const float k013 = math_utils::computeCurvature(0.f, 0.f,
                                                        cls1.xCoordinate, cls1.yCoordinate,
                                                        cls3.xCoordinate, cls3.yCoordinate);
        const float dk = o2::gpu::CAMath::Abs(k123 - k013);
        const float dR2 = 0.25f * (math_utils::Sq(cls1.xCoordinate - cls3.xCoordinate) + math_utils::Sq(cls1.yCoordinate - cls3.yCoordinate));
        const float tgl2 = math_utils::Sq(0.5 * (currentTracklet.tanLambda + nextTracklet.tanLambda));
        const float snl2 = tgl2 / (1.f + tgl2);
        const float kSigma2Pos = 6.f * resMid2 / math_utils::Sq(dR2);
        const float kSigma2MS = msMid2 / (dR2 * snl2);
        const float kNSigma = mRecoParams[iteration].params.NSigmaCut * o2::gpu::CAMath::Sqrt(kSigma2Pos + kSigma2MS);
        if (dk >= kNSigma) {
          continue;
        }

        /// Track seed preparation. Clusters are numbered progressively from the innermost going outward.
        const int clusId[3] = {cls1.clusterId, cls2.clusterId, cls3.clusterId};
        const auto& cluster1_glo = mTimeFrame->getUnsortedClusters()[iLayer][clusId[0]];
        const auto& cluster2_glo = mTimeFrame->getUnsortedClusters()[iLayer + 1][clusId[1]];
        const auto& cluster3_tf = mTimeFrame->getTrackingFrameInfoOnLayer(iLayer + 2)[clusId[2]];
        auto track{buildTrackSeed(cluster1_glo, cluster2_glo, cluster3_tf)};

        float chi2{0.f};
        bool good{false};
        for (int iC{2}; iC--;) {
          const TrackingFrameInfo& trackingHit = mTimeFrame->getTrackingFrameInfoOnLayer(iLayer + iC)[clusId[iC]];

          if (!track.rotate(trackingHit.alphaTrackingFrame)) {
            break;
          }

          if (!track.propagateTo(trackingHit.xTrackingFrame, getBz())) {
            break;
          }

          if (!track.correctForMaterial(mRecoParams[iteration].params.LayerxX0[iLayer + iC], mRecoParams[iteration].params.LayerxX0[iLayer] * constants::Radl * constants::Rho, true)) {
            break;
          }

          const auto predChi2{track.getPredictedChi2Quiet(trackingHit.positionTrackingFrame, trackingHit.covarianceTrackingFrame)};
          if (!iC && predChi2 > mRecoParams[iteration].params.MaxChi2ClusterAttachment) {
            break;
          }

          if (!track.o2::track::TrackParCov::update(trackingHit.positionTrackingFrame, trackingHit.covarianceTrackingFrame)) {
            break;
          }

          good = !iC;
          chi2 += predChi2;
        }
        if (good) {
          if constexpr (decltype(Tag)::value == PassMode::OnePass::value) {
            layerCells.emplace_back(iLayer, clusId[0], clusId[1], clusId[2], iTracklet, iNextTracklet, track, chi2);
            ++foundCells;
          } else if constexpr (decltype(Tag)::value == PassMode::TwoPassCount::value) {
            ++foundCells;
          } else if constexpr (decltype(Tag)::value == PassMode::TwoPassInsert::value) {
            layerCells[offset++] = CellSeedN(iLayer, clusId[0], clusId[1], clusId[2], iTracklet, iNextTracklet, track, chi2);
          } else {
            static_assert(false, "Unknown mode!");
          }
        }
      }
      return foundCells;
    };

    tbb::parallel_for(0, mRecoParams[iteration].params.CellsPerRoad(), [&](const int iLayer) {
      if (mTimeFrame->getTracklets()[iLayer + 1].empty() ||
          mTimeFrame->getTracklets()[iLayer].empty()) {
        return;
      }

      auto& layerCells = mTimeFrame->getCells()[iLayer];
      const int currentLayerTrackletsNum{static_cast<int>(mTimeFrame->getTracklets()[iLayer].size())};
      bounded_vector<int> perTrackletCount(currentLayerTrackletsNum + 1, 0, mMemoryPool.get());
      if (mTaskArena->max_concurrency() <= 1) {
        for (int iTracklet{0}; iTracklet < currentLayerTrackletsNum; ++iTracklet) {
          perTrackletCount[iTracklet] = forTrackletCells(PassMode::OnePass{}, iLayer, layerCells, iTracklet);
        }
        std::exclusive_scan(perTrackletCount.begin(), perTrackletCount.end(), perTrackletCount.begin(), 0);
      } else {
        tbb::parallel_for(0, currentLayerTrackletsNum, [&](const int iTracklet) {
          perTrackletCount[iTracklet] = forTrackletCells(PassMode::TwoPassCount{}, iLayer, layerCells, iTracklet);
        });

        std::exclusive_scan(perTrackletCount.begin(), perTrackletCount.end(), perTrackletCount.begin(), 0);
        auto totalCells{perTrackletCount.back()};
        if (totalCells == 0) {
          return;
        }
        layerCells.resize(totalCells);

        tbb::parallel_for(0, currentLayerTrackletsNum, [&](const int iTracklet) {
          int offset = perTrackletCount[iTracklet];
          if (offset == perTrackletCount[iTracklet + 1]) {
            return;
          }
          forTrackletCells(PassMode::TwoPassInsert{}, iLayer, layerCells, iTracklet, offset);
        });
      }

      if (iLayer > 0) {
        auto& lut = mTimeFrame->getCellsLookupTable()[iLayer - 1];
        lut.resize(currentLayerTrackletsNum + 1);
        std::copy_n(perTrackletCount.begin(), currentLayerTrackletsNum + 1, lut.begin());
      }
    });

    /// Create cells labels
    if ((mTimeFrame->hasMCinformation() && mRecoParams[iteration].params.createArtefactLabels) || OPTIMISATION_ANY) {
      tbb::parallel_for(0, mRecoParams[iteration].params.CellsPerRoad(), [&](const int iLayer) {
        mTimeFrame->getCellsLabel(iLayer).reserve(mTimeFrame->getCells()[iLayer].size());
        for (const auto& cell : mTimeFrame->getCells()[iLayer]) {
          MCCompLabel currentLab{mTimeFrame->getTrackletsLabel(iLayer)[cell.getFirstTrackletIndex()]};
          MCCompLabel nextLab{mTimeFrame->getTrackletsLabel(iLayer + 1)[cell.getSecondTrackletIndex()]};
          mTimeFrame->getCellsLabel(iLayer).emplace_back(currentLab == nextLab ? currentLab : MCCompLabel());
        }
      });
    }
  });
}

template <int NLayers>
void TrackerTraits<NLayers>::findCellSeeds(const int iteration)
{
  const auto& propagator = o2::base::Propagator::Instance();
  auto& cells = mTimeFrame->getCells()[0];
  LOGP(info, "received {} | {} tracklets", mTimeFrame->getTracklets()[0].size(), mTimeFrame->getTracklets()[1].size());
  LOGP(info, "received {} cells", cells.size());

  mTaskArena->execute([&] {
    // we get cells from the first three layers
    bounded_vector<LinearizedTrack> ltracks;
    ltracks.resize(cells.size());
    tbb::parallel_for(size_t(0), cells.size(), [&](size_t iCell) {
      // relate each cell to the imposed mean vertex or assumed beamline if the former is not provided
      auto& cell = cells[iCell];
      // impose general quality cuts
      if (cell.getPt() < mRecoParams[iteration].params.SeedingMinPtTrk) {
        return;
      }
      dataformats::VertexBase vtx;
      dataformats::DCA dca;
      float z = cell.getZAt(0., getBz());
      if (z < -999.f) { // if outside of acceptance impose mean vertex z
        z = mTimeFrame->getMeanVertex().getZ();
      }
      if (mTimeFrame->hasMeanVertex()) {
        mTimeFrame->getMeanVertexConstraint()->setMeanXYVertexAtZ(vtx, z);
      }
      if (!propagator->propagateToDCA(vtx, cell, getBz(), 2.0f, o2::base::Propagator::MatCorrType::USEMatCorrLUT, &dca, nullptr, 0, mRecoParams[iteration].params.SeedingDCATolerance)) {
        ltracks[iCell].markDead();
        return;
      }
      if (dca.getY() * dca.getY() / (dca.getSigmaY2()) >= mRecoParams[iteration].params.SeedingDCAMaxPull) {
        ltracks[iCell].markDead();
        return;
      }
      ltracks[iCell] = LinearizedTrack(cell, iCell);
      if (ltracks[iCell].isDead()) {
        return;
      }
      const auto& cls = cell.getClusters();
      const int sta = cell.getUserField();
      int startBC = std::numeric_limits<int>::max(), endBC = std::numeric_limits<int>::min();
      for (int i{sta}; i < 3; ++i) {
        int rof = mTimeFrame->getClusterROF(i, cls[i]);
        int rofStartBC = mTimeFrame->getROFOverlapTableView().getLayer(i).getROFStartInBC(rof);
        int rofEndBC = mTimeFrame->getROFOverlapTableView().getLayer(i).getROFEndInBC(rof);
        startBC = o2::gpu::CAMath::Min(startBC, rofStartBC);
        endBC = o2::gpu::CAMath::Max(endBC, rofEndBC);
      }
      if (endBC - startBC < 0) { // this should not happen
        ltracks[iCell].markDead();
      }
      ltracks[iCell].time.setTimeStamp(startBC);
      ltracks[iCell].time.setTimeStampError(endBC - startBC);
    });
    tbb::parallel_sort(ltracks.begin(), ltracks.end(), [](const LinearizedTrack& a, const LinearizedTrack& b) {
      const bool aDead = a.isDead();
      const bool bDead = b.isDead();
      // all dead tracks are sorted to the end
      if (aDead != bDead) {
        return !aDead; // a < b only if a is alive and b is dead
      }
      // sort them in time and then increasing z
      if (!aDead) {
        const auto ta = a.time.getTimeStamp();
        const auto tb = b.time.getTimeStamp();
        if (ta != tb) {
          return ta < tb;
        }
        return a.z < b.z;
      }
      return false;
    });
    // drop all dead tracks
    auto firstDead = std::partition_point(ltracks.begin(), ltracks.end(), [](const LinearizedTrack& t) { return !t.isDead(); });
    ltracks.erase(firstDead, ltracks.end());
    // do DBscan
    // scan gives association in time and z
    const auto dbRes = mDBScan.cluster(ltracks.data(), ltracks.size());
    bounded_vector<VertexSeed> vtxSeeds(dbRes.nClusters, mMemoryPool.get());
    bounded_vector<VertexLabel> vtxSeedsLbl(dbRes.nClusters, mMemoryPool.get());
    // we only care about the source&event of the tracks, not the trackId
    auto composeVtxLabel = [](o2::MCCompLabel& lbl) -> void {
      lbl.set(o2::MCCompLabel::maxTrackID(), lbl.getEventID(), lbl.getSourceID(), lbl.isFake());
    };
    tbb::parallel_for(0, dbRes.nClusters, [&](const int32_t iCls) {
      // create vertex seeds based on meanvertex and scanned z which one can take as start
      auto& seed = vtxSeeds[iCls];
      seed.setZ(dbRes.zCentroids[iCls]);
      if (mTimeFrame->hasMeanVertex()) {
        mTimeFrame->getMeanVertexConstraint()->setMeanXYVertexAtZ(seed, seed.getZ());
      }
      seed.idx = iCls;
      seed.tukeyC = mRecoParams[iteration].params.SeedingTukeyStartIter;
      seed.scaleSig2ITuk2I = 1.f / (seed.tukeyC * seed.tukeyC);
      const auto& range = dbRes.ranges[iCls];
      // fit iteratively
      for (int iter{0}; iter < mRecoParams[iteration].params.SeedingMaxFitIter; ++iter) {
        seed.iteration = iter;
        // 1. update tukey scaling
        seed.updateTukeyScale(ltracks.data(), dbRes.labels.data(), range);
        // 2. reset for new iteration with new weights
        seed.resetForNewIteration();
        // 3. account all tracks
        for (uint32_t entry{range.getFirstEntry()}; entry < range.getEntriesBound(); ++entry) {
          const auto& lt = ltracks[entry];
          if (lt.isDead() || dbRes.labels[entry] != iCls) {
            continue;
          }
          seed.accountTrack(lt);
        }
        // 4. if needed relax scaling
        if (seed.getNContributors() < mRecoParams[iteration].params.SeedingMinTracksIter) {
          seed.status = VertexSeed::kIterateRelaxScale;
          continue;
        }
        // 5. impose mean vertex constraint if wanted
        if (mTimeFrame->hasMeanVertex()) {
          const auto& err = mTimeFrame->getMeanVertexInvErr();
          seed.C(0, 0) += err[0]; // cxx
          seed.C(0, 1) += err[1]; // cxy
          seed.C(1, 1) += err[2]; // cyy
          float x = mTimeFrame->getMeanVertexConstraint()->getXAtZ(seed.getZ());
          float y = mTimeFrame->getMeanVertexConstraint()->getYAtZ(seed.getZ());
          seed.b(0) += err[0] * x + err[1] * y;
          seed.b(1) += err[1] * x + err[2] * y;
        }
        // 6. solve LS fit
        seed.solveVertex();
        // 7. check for convergence
        float avgChi2 = seed.wghChi2 / o2::gpu::CAMath::Max(1.f, seed.wghSum);
        if (avgChi2 < mRecoParams[iteration].params.SeedingMaxChi2Iter) {
          seed.status = VertexSeed::kConverged;
          break;
        }
        seed.status = VertexSeed::kIterateFurther;
      }
      // find time bracket of vertex by looking at used tracks
      uint32_t startBC = std::numeric_limits<uint32_t>::max(), endBC = 0u;
      for (uint32_t entry{range.getFirstEntry()}; entry < range.getEntriesBound(); ++entry) {
        const auto& lt = ltracks[entry];
        if (lt.isDead() || dbRes.labels[entry] != iCls) {
          continue;
        }
        float chi2Red = lt.getChi2(seed);
        if (chi2Red > mRecoParams[iteration].params.SeedingMaxChi2Iter) {
          continue;
        }
        startBC = o2::gpu::CAMath::Min(startBC, lt.time.getTimeStamp());
        endBC = o2::gpu::CAMath::Max(endBC, lt.time.getTimeStamp() + lt.time.getTimeStampError());
      }
      seed.getTimeStamp().setTimeStamp(startBC);
      seed.getTimeStamp().setTimeStampError(endBC - startBC);
      // add additional errors
      for (int i{0}; i < dataformats::VertexBase::kNCov; ++i) {
        seed.setCov(seed.getCov(i) + mRecoParams[iteration].params.SeedingVertexExtraErr2[i], i);
      }
      // if mc is present calculate labels and purity
      if (mTimeFrame->hasMCinformation()) {
        // use boyer-moore voting
        int accepted{0}, weight{0};
        o2::MCCompLabel lbl;
        for (uint32_t entry{range.getFirstEntry()}; entry < range.getEntriesBound(); ++entry) {
          const auto& lt = ltracks[entry];
          if (lt.isDead() || dbRes.labels[entry] != iCls) {
            continue;
          }
          if (seed.acceptTrack(lt)) {
            ++accepted;
            // create cell label
            const auto& cell = mTimeFrame->getCells()[0][lt.cellIdx];
            o2::MCCompLabel cl;
            for (const auto& lab1 : mTimeFrame->getClusterLabels(0, cell.getFirstClusterIndex())) {
              for (const auto& lab2 : mTimeFrame->getClusterLabels(1, cell.getSecondClusterIndex())) {
                if (lab1 == lab2 && lab1.isValid()) {
                  cl = lab1;
                  break;
                }
              }
              if (cl.isValid()) {
                break;
              }
            }
            composeVtxLabel(cl); // normalize label
            if (weight == 0) {
              lbl = cl;
              weight = 1;
            } else {
              (cl == lbl) ? ++weight : --weight;
            }
          }
        }
        if (accepted > 0) {
          vtxSeedsLbl[iCls] = {lbl, static_cast<float>(weight) / static_cast<float>(accepted)};
        } else {
          vtxSeedsLbl[iCls] = {lbl, 0};
        }
      }
    });
    // TODO reduce debris vertices (eliminate small mult vertices within vicinity of high multiplicity ones)
    // Sort vertices in time and multiplicity
    for (int i{0}; i < dbRes.nClusters; ++i) {
      if (vtxSeeds[i].status == VertexSeed::kKilled) {
        continue;
      }
      mTimeFrame->addPrimaryVertex(vtxSeeds[i]);
      if (mTimeFrame->hasMCinformation()) {
        mTimeFrame->addPrimaryVertexLabel(vtxSeedsLbl[i]);
      }
    }
    // if no mean vertex constraint imposed, use rolling weighted average (KF-filter)
    // this should theoretical stabelize very quickly
    if (!mTimeFrame->hasMeanVertex()) {
      ROOT::Math::SMatrix<float, 3, 3> I;
      I(0, 0) = 1.f;
      I(1, 1) = 1.f;
      I(2, 2) = 1.f;
      auto& avg = mTimeFrame->getMeanVertexRolling();
      ROOT::Math::SVector<float, 3> x;
      x(0) = avg.getX();
      x(1) = avg.getY();
      x(2) = avg.getZ();
      ROOT::Math::SMatrix<float, 3, 3> P;
      P(0, 0) = avg.getCov(dataformats::VertexBase::kCovXX);
      P(0, 1) = avg.getCov(dataformats::VertexBase::kCovXY);
      P(1, 0) = avg.getCov(dataformats::VertexBase::kCovXY);
      P(1, 1) = avg.getCov(dataformats::VertexBase::kCovYY);
      P(0, 2) = avg.getCov(dataformats::VertexBase::kCovXZ);
      P(2, 0) = avg.getCov(dataformats::VertexBase::kCovXZ);
      P(1, 2) = avg.getCov(dataformats::VertexBase::kCovYZ);
      P(2, 1) = avg.getCov(dataformats::VertexBase::kCovYZ);
      P(2, 2) = avg.getCov(dataformats::VertexBase::kCovZZ);
      for (const auto& vtx : vtxSeeds) {
        if (vtx.getNContributors() < mRecoParams[iteration].params.SeedingMinContrib) {
          continue;
        }
        ROOT::Math::SVector<float, 3> z;
        z(0) = vtx.getX();
        z(1) = vtx.getY();
        z(2) = vtx.getZ();
        ROOT::Math::SMatrix<float, 3, 3> R;
        R(0, 0) = vtx.getCov(dataformats::VertexBase::kCovXX);
        R(0, 1) = vtx.getCov(dataformats::VertexBase::kCovXY);
        R(1, 0) = vtx.getCov(dataformats::VertexBase::kCovXY);
        R(1, 1) = vtx.getCov(dataformats::VertexBase::kCovYY);
        R(0, 2) = vtx.getCov(dataformats::VertexBase::kCovXZ);
        R(2, 0) = vtx.getCov(dataformats::VertexBase::kCovXZ);
        R(1, 2) = vtx.getCov(dataformats::VertexBase::kCovYZ);
        R(2, 1) = vtx.getCov(dataformats::VertexBase::kCovYZ);
        R(2, 2) = vtx.getCov(dataformats::VertexBase::kCovZZ);
        ROOT::Math::SMatrix<float, 3, 3> S = P + R; // innovation
        if (!S.Invert()) {
          continue;
        }
        ROOT::Math::SMatrix<float, 3, 3> K = P * S; // gain
        ROOT::Math::SVector<float, 3> y = z - x;    // innovation
        x += K * y;                                 // state update
        P = (I - K) * P;                            // cov update
      }
      // update avg
      // do we need to symmetrize the matrix, probably not errors are anyways small
      avg.setX(x[0]);
      avg.setY(x[1]);
      avg.setZ(x[2]);
      avg.setCov(P(0, 0), dataformats::VertexBase::kCovXX);
      avg.setCov(P(0, 1), dataformats::VertexBase::kCovXY);
      avg.setCov(P(0, 2), dataformats::VertexBase::kCovXZ);
      avg.setCov(P(1, 1), dataformats::VertexBase::kCovYY);
      avg.setCov(P(1, 2), dataformats::VertexBase::kCovYZ);
      avg.setCov(P(2, 2), dataformats::VertexBase::kCovZZ);
    }
  });
  sortSeeds();
}

template <int NLayers>
void TrackerTraits<NLayers>::findCellsNeighbours(const int iteration)
{
  mTaskArena->execute([&] {
    for (int iLayer{0}; iLayer < mRecoParams[iteration].params.NeighboursPerRoad(); ++iLayer) {
      deepVectorClear(mTimeFrame->getCellsNeighbours()[iLayer]);
      deepVectorClear(mTimeFrame->getCellsNeighboursLUT()[iLayer]);
      if (mTimeFrame->getCells()[iLayer + 1].empty() ||
          mTimeFrame->getCellsLookupTable()[iLayer].empty()) {
        continue;
      }

      int nCells{static_cast<int>(mTimeFrame->getCells()[iLayer].size())};
      bounded_vector<Neighbor> cellsNeighbours(mMemoryPool.get());

      auto forCellNeighbour = [&](auto Tag, int iCell, int offset = 0) -> int {
        const auto& currentCellSeed{mTimeFrame->getCells()[iLayer][iCell]};
        const int nextLayerTrackletIndex{currentCellSeed.getSecondTrackletIndex()};
        const int nextLayerFirstCellIndex{mTimeFrame->getCellsLookupTable()[iLayer][nextLayerTrackletIndex]};
        const int nextLayerLastCellIndex{mTimeFrame->getCellsLookupTable()[iLayer][nextLayerTrackletIndex + 1]};
        int foundNextCells{0};
        for (int iNextCell{nextLayerFirstCellIndex}; iNextCell < nextLayerLastCellIndex; ++iNextCell) {
          auto nextCellSeed{mTimeFrame->getCells()[iLayer + 1][iNextCell]}; /// copy
          if (nextCellSeed.getFirstTrackletIndex() != nextLayerTrackletIndex) {
            break;
          }

          // by construction for each cell we already know that their tracklets are compatible in time
          // and above we check that the cells use the same middle tracklet so we just have to make sure that the
          // outermost tracklets are compatible in time
          const auto& trkl00 = mTimeFrame->getTracklets()[iLayer][currentCellSeed.getFirstTrackletIndex()];
          const auto& trkl11 = mTimeFrame->getTracklets()[iLayer + 2][nextCellSeed.getSecondTrackletIndex()];
          if (!mTimeFrame->getROFOverlapTableView().isCompatible(iLayer, trkl00.rof[0], iLayer + 2, trkl11.rof[1])) {
            continue;
          }

          if (!nextCellSeed.rotate(currentCellSeed.getAlpha()) ||
              !nextCellSeed.propagateTo(currentCellSeed.getX(), getBz())) {
            continue;
          }
          const float chi2 = currentCellSeed.getPredictedChi2(nextCellSeed); /// TODO: switch to the chi2 wrt cluster to avoid correlation
          const bool accept = chi2 <= mRecoParams[iteration].params.MaxChi2ClusterAttachment;

          debugFindCellsNeighbours(iteration, iLayer, iCell, iNextCell, accept);

          if (accept) {
            if constexpr (decltype(Tag)::value == PassMode::OnePass::value) {
              cellsNeighbours.emplace_back(iCell, iNextCell, currentCellSeed.getLevel() + 1);
            } else if constexpr (decltype(Tag)::value == PassMode::TwoPassCount::value) {
              ++foundNextCells;
            } else if constexpr (decltype(Tag)::value == PassMode::TwoPassInsert::value) {
              cellsNeighbours[offset++] = {iCell, iNextCell, currentCellSeed.getLevel() + 1};
            } else {
              static_assert(false, "Unknown mode!");
            }
          }
        }
        return foundNextCells;
      };

      if (mTaskArena->max_concurrency() <= 1) {
        for (int iCell{0}; iCell < nCells; ++iCell) {
          forCellNeighbour(PassMode::OnePass{}, iCell);
        }
      } else {
        bounded_vector<int> perCellCount(nCells + 1, 0, mMemoryPool.get());
        tbb::parallel_for(0, nCells, [&](const int iCell) {
          perCellCount[iCell] = forCellNeighbour(PassMode::TwoPassCount{}, iCell);
        });

        std::exclusive_scan(perCellCount.begin(), perCellCount.end(), perCellCount.begin(), 0);
        int totalCellNeighbours = perCellCount.back();
        if (totalCellNeighbours == 0) {
          deepVectorClear(mTimeFrame->getCellsNeighbours()[iLayer]);
          continue;
        }
        cellsNeighbours.resize(totalCellNeighbours);

        tbb::parallel_for(0, nCells, [&](const int iCell) {
          int offset = perCellCount[iCell];
          if (offset == perCellCount[iCell + 1]) {
            return;
          }
          forCellNeighbour(PassMode::TwoPassInsert{}, iCell, offset);
        });
      }

      if (cellsNeighbours.empty()) {
        continue;
      }

      tbb::parallel_sort(cellsNeighbours.begin(), cellsNeighbours.end(), [](const auto& a, const auto& b) {
        return a.nextCell < b.nextCell;
      });

      auto& cellsNeighbourLUT = mTimeFrame->getCellsNeighboursLUT()[iLayer];
      cellsNeighbourLUT.assign(mTimeFrame->getCells()[iLayer + 1].size(), 0);
      for (const auto& neigh : cellsNeighbours) {
        ++cellsNeighbourLUT[neigh.nextCell];
      }
      std::inclusive_scan(cellsNeighbourLUT.begin(), cellsNeighbourLUT.end(), cellsNeighbourLUT.begin());

      mTimeFrame->getCellsNeighbours()[iLayer].reserve(cellsNeighbours.size());
      std::ranges::transform(cellsNeighbours, std::back_inserter(mTimeFrame->getCellsNeighbours()[iLayer]), [](const auto& neigh) { return neigh.cell; });

      for (auto it = cellsNeighbours.begin(); it != cellsNeighbours.end();) {
        int cellIdx = it->nextCell;
        int maxLvl = it->level;
        while (++it != cellsNeighbours.end() && it->nextCell == cellIdx) {
          maxLvl = std::max(maxLvl, it->level);
        }
        mTimeFrame->getCells()[iLayer + 1][cellIdx].setLevel(maxLvl);
      }
    }
  });
}

template <int NLayers>
void TrackerTraits<NLayers>::processNeighbours(int iteration, int iLayer, int iLevel, const bounded_vector<CellSeedN>& currentCellSeed, const bounded_vector<int>& currentCellId, bounded_vector<CellSeedN>& updatedCellSeeds, bounded_vector<int>& updatedCellsIds)
{
  const auto& propagator = o2::base::Propagator::Instance();

  mTaskArena->execute([&] {
    auto forCellNeighbours = [&](auto Tag, int iCell, int offset = 0) -> int {
      const auto& currentCell{currentCellSeed[iCell]};

      if constexpr (decltype(Tag)::value != PassMode::TwoPassInsert::value) {
        if (currentCell.getLevel() != iLevel) {
          return 0;
        }
        if (currentCellId.empty() && (mTimeFrame->isClusterUsed(iLayer, currentCell.getFirstClusterIndex()) ||
                                      mTimeFrame->isClusterUsed(iLayer + 1, currentCell.getSecondClusterIndex()) ||
                                      mTimeFrame->isClusterUsed(iLayer + 2, currentCell.getThirdClusterIndex()))) {
          return 0; /// this we do only on the first iteration, hence the check on currentCellId
        }
      }

      const int cellId = currentCellId.empty() ? iCell : currentCellId[iCell];
      const int startNeighbourId{cellId ? mTimeFrame->getCellsNeighboursLUT()[iLayer - 1][cellId - 1] : 0};
      const int endNeighbourId{mTimeFrame->getCellsNeighboursLUT()[iLayer - 1][cellId]};
      int foundSeeds{0};
      for (int iNeighbourCell{startNeighbourId}; iNeighbourCell < endNeighbourId; ++iNeighbourCell) {
        const int neighbourCellId = mTimeFrame->getCellsNeighbours()[iLayer - 1][iNeighbourCell];
        const auto& neighbourCell = mTimeFrame->getCells()[iLayer - 1][neighbourCellId];
        if (neighbourCell.getSecondTrackletIndex() != currentCell.getFirstTrackletIndex()) {
          continue;
        }
        if (mTimeFrame->isClusterUsed(iLayer - 1, neighbourCell.getFirstClusterIndex())) {
          continue;
        }
        if (currentCell.getLevel() - 1 != neighbourCell.getLevel()) {
          continue;
        }

        /// Let's start the fitting procedure
        CellSeedN seed{currentCell};
        const auto& trHit = mTimeFrame->getTrackingFrameInfoOnLayer(iLayer - 1)[neighbourCell.getFirstClusterIndex()];

        if (!seed.rotate(trHit.alphaTrackingFrame)) {
          continue;
        }

        if (!propagator->propagateToX(seed, trHit.xTrackingFrame, getBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, mRecoParams[iteration].params.CorrType)) {
          continue;
        }

        if (mRecoParams[iteration].params.CorrType == o2::base::PropagatorF::MatCorrType::USEMatCorrNONE) {
          if (!seed.correctForMaterial(mRecoParams[iteration].params.LayerxX0[iLayer - 1], mRecoParams[iteration].params.LayerxX0[iLayer - 1] * constants::Radl * constants::Rho, true)) {
            continue;
          }
        }

        auto predChi2{seed.getPredictedChi2Quiet(trHit.positionTrackingFrame, trHit.covarianceTrackingFrame)};
        if ((predChi2 > mRecoParams[iteration].params.MaxChi2ClusterAttachment) || predChi2 < 0.f) {
          continue;
        }
        seed.setChi2(seed.getChi2() + predChi2);
        if (!seed.o2::track::TrackParCov::update(trHit.positionTrackingFrame, trHit.covarianceTrackingFrame)) {
          continue;
        }

        if constexpr (decltype(Tag)::value != PassMode::TwoPassCount::value) {
          seed.getClusters()[iLayer - 1] = neighbourCell.getFirstClusterIndex();
          seed.setLevel(neighbourCell.getLevel());
          seed.setFirstTrackletIndex(neighbourCell.getFirstTrackletIndex());
          seed.setSecondTrackletIndex(neighbourCell.getSecondTrackletIndex());
        }

        if constexpr (decltype(Tag)::value == PassMode::OnePass::value) {
          updatedCellSeeds.push_back(seed);
          updatedCellsIds.push_back(neighbourCellId);
        } else if constexpr (decltype(Tag)::value == PassMode::TwoPassCount::value) {
          ++foundSeeds;
        } else if constexpr (decltype(Tag)::value == PassMode::TwoPassInsert::value) {
          updatedCellSeeds[offset] = seed;
          updatedCellsIds[offset++] = neighbourCellId;
        } else {
          static_assert(false, "Unknown mode!");
        }
      }
      return foundSeeds;
    };

    const int nCells = static_cast<int>(currentCellSeed.size());
    if (mTaskArena->max_concurrency() <= 1) {
      for (int iCell{0}; iCell < nCells; ++iCell) {
        forCellNeighbours(PassMode::OnePass{}, iCell);
      }
    } else {
      bounded_vector<int> perCellCount(nCells + 1, 0, mMemoryPool.get());
      tbb::parallel_for(0, nCells, [&](const int iCell) {
        perCellCount[iCell] = forCellNeighbours(PassMode::TwoPassCount{}, iCell);
      });

      std::exclusive_scan(perCellCount.begin(), perCellCount.end(), perCellCount.begin(), 0);
      auto totalNeighbours{perCellCount.back()};
      if (totalNeighbours == 0) {
        return;
      }
      updatedCellSeeds.resize(totalNeighbours);
      updatedCellsIds.resize(totalNeighbours);

      tbb::parallel_for(0, nCells, [&](const int iCell) {
        int offset = perCellCount[iCell];
        if (offset == perCellCount[iCell + 1]) {
          return;
        }
        forCellNeighbours(PassMode::TwoPassInsert{}, iCell, offset);
      });
    }
  });
}

template <int NLayers>
void TrackerTraits<NLayers>::findRoads(const int iteration)
{
  // we don't need tracklets anymore
  for (int iLayer{0}; iLayer < mRecoParams[iteration].params.TrackletsPerRoad(); ++iLayer) {
    deepVectorClear(mTimeFrame->getTracklets()[iLayer]);
    deepVectorClear(mTimeFrame->getTrackletsLookupTable()[iLayer]);
    if ((mTimeFrame->hasMCinformation() && mRecoParams[iteration].params.createArtefactLabels) || OPTIMISATION_ANY) {
      deepVectorClear(mTimeFrame->getTrackletsLabel(iLayer));
    }
  }

  bounded_vector<bounded_vector<int>> firstClusters(mRecoParams[iteration].params.NLayers, bounded_vector<int>(mMemoryPool.get()), mMemoryPool.get());
  bounded_vector<bounded_vector<int>> sharedFirstClusters(mRecoParams[iteration].params.NLayers, bounded_vector<int>(mMemoryPool.get()), mMemoryPool.get());
  firstClusters.resize(mRecoParams[iteration].params.NLayers);
  sharedFirstClusters.resize(mRecoParams[iteration].params.NLayers);
  for (int startLevel{mRecoParams[iteration].params.CellsPerRoad()}; startLevel >= mRecoParams[iteration].params.CellMinimumLevel(); --startLevel) {

    auto seedFilter = [&](const auto& seed) {
      return seed.getQ2Pt() <= 1.e3 && seed.getChi2() <= mRecoParams[iteration].params.MaxChi2NDF * (((startLevel + 2) * 2) - 5);
    };

    bounded_vector<CellSeedN> trackSeeds(mMemoryPool.get());
    for (int startLayer{mRecoParams[iteration].params.NeighboursPerRoad()}; startLayer >= startLevel - 1; --startLayer) {
      if ((mRecoParams[iteration].params.StartLayerMask & (1 << (startLayer + 2))) == 0) {
        continue;
      }

      bounded_vector<int> lastCellId(mMemoryPool.get()), updatedCellId(mMemoryPool.get());
      bounded_vector<CellSeedN> lastCellSeed(mMemoryPool.get()), updatedCellSeed(mMemoryPool.get());

      processNeighbours(iteration, startLayer, startLevel, mTimeFrame->getCells()[startLayer], lastCellId, updatedCellSeed, updatedCellId);

      int level = startLevel;
      for (int iLayer{startLayer - 1}; iLayer > 0 && level > 2; --iLayer) {
        lastCellSeed.swap(updatedCellSeed);
        lastCellId.swap(updatedCellId);
        deepVectorClear(updatedCellSeed); /// tame the memory peaks
        deepVectorClear(updatedCellId);   /// tame the memory peaks
        processNeighbours(iteration, iLayer, --level, lastCellSeed, lastCellId, updatedCellSeed, updatedCellId);
      }
      deepVectorClear(lastCellId);   /// tame the memory peaks
      deepVectorClear(lastCellSeed); /// tame the memory peaks

      if (!updatedCellSeed.empty()) {
        trackSeeds.reserve(trackSeeds.size() + std::count_if(updatedCellSeed.begin(), updatedCellSeed.end(), seedFilter));
        std::copy_if(updatedCellSeed.begin(), updatedCellSeed.end(), std::back_inserter(trackSeeds), seedFilter);
      }
    }

    if (trackSeeds.empty()) {
      continue;
    }

    bounded_vector<TrackITSExt> tracks(mMemoryPool.get());
    mTaskArena->execute([&] {
      auto forSeed = [&](auto Tag, int iSeed, int offset = 0) {
        TrackITSExt temporaryTrack = seedTrackForRefit(trackSeeds[iSeed]);
        o2::track::TrackPar linRef{temporaryTrack};
        bool fitSuccess = fitTrack(temporaryTrack, 0, mRecoParams[iteration].params.NLayers, 1, mRecoParams[iteration].params.MaxChi2ClusterAttachment, mRecoParams[iteration].params.MaxChi2NDF, o2::constants::math::VeryBig, 0, &linRef);
        if (!fitSuccess) {
          return 0;
        }
        temporaryTrack.getParamOut() = temporaryTrack.getParamIn();
        linRef = temporaryTrack.getParamOut(); // use refitted track as lin.reference
        temporaryTrack.resetCovariance();
        temporaryTrack.setCov(temporaryTrack.getQ2Pt() * temporaryTrack.getQ2Pt() * temporaryTrack.getCov()[o2::track::CovLabels::kSigQ2Pt2], o2::track::CovLabels::kSigQ2Pt2);
        temporaryTrack.setChi2(0);
        fitSuccess = fitTrack(temporaryTrack, mRecoParams[iteration].params.NLayers - 1, -1, -1, mRecoParams[iteration].params.MaxChi2ClusterAttachment, mRecoParams[iteration].params.MaxChi2NDF, 50.f, 0, &linRef);
        if (!fitSuccess || temporaryTrack.getPt() < mRecoParams[iteration].params.MinPt[mRecoParams[iteration].params.NLayers - temporaryTrack.getNClusters()]) {
          return 0;
        }
        if constexpr (decltype(Tag)::value == PassMode::OnePass::value) {
          tracks.push_back(temporaryTrack);
        } else if constexpr (decltype(Tag)::value == PassMode::TwoPassCount::value) {
          // nothing to do
        } else if constexpr (decltype(Tag)::value == PassMode::TwoPassInsert::value) {
          tracks[offset] = temporaryTrack;
        } else {
          static_assert(false, "Unknown mode!");
        }
        return 1;
      };

      const int nSeeds = static_cast<int>(trackSeeds.size());
      if (mTaskArena->max_concurrency() <= 1) {
        for (int iSeed{0}; iSeed < nSeeds; ++iSeed) {
          forSeed(PassMode::OnePass{}, iSeed);
        }
      } else {
        bounded_vector<int> perSeedCount(nSeeds + 1, 0, mMemoryPool.get());
        tbb::parallel_for(0, nSeeds, [&](const int iSeed) {
          perSeedCount[iSeed] = forSeed(PassMode::TwoPassCount{}, iSeed);
        });

        std::exclusive_scan(perSeedCount.begin(), perSeedCount.end(), perSeedCount.begin(), 0);
        auto totalTracks{perSeedCount.back()};
        if (totalTracks == 0) {
          return;
        }
        tracks.resize(totalTracks);

        tbb::parallel_for(0, nSeeds, [&](const int iSeed) {
          if (perSeedCount[iSeed] == perSeedCount[iSeed + 1]) {
            return;
          }
          forSeed(PassMode::TwoPassInsert{}, iSeed, perSeedCount[iSeed]);
        });
      }

      deepVectorClear(trackSeeds);
      tbb::parallel_sort(tracks.begin(), tracks.end(), [](const auto& a, const auto& b) {
        return a.getChi2() < b.getChi2();
      });
    });

    for (auto& track : tracks) {
      int nShared = 0;
      bool isFirstShared{false};
      int firstLayer{-1}, firstCluster{-1};
      for (int iLayer{0}; iLayer < mRecoParams[iteration].params.NLayers; ++iLayer) {
        if (track.getClusterIndex(iLayer) == constants::UnusedIndex) {
          continue;
        }
        bool isShared = mTimeFrame->isClusterUsed(iLayer, track.getClusterIndex(iLayer));
        nShared += int(isShared);
        if (firstLayer < 0) {
          firstCluster = track.getClusterIndex(iLayer);
          isFirstShared = isShared && mRecoParams[iteration].params.AllowSharingFirstCluster && std::find(firstClusters[iLayer].begin(), firstClusters[iLayer].end(), firstCluster) != firstClusters[iLayer].end();
          firstLayer = iLayer;
        }
      }

      /// do not account for the first cluster in the shared clusters number if it is allowed
      if (nShared - int(isFirstShared && mRecoParams[iteration].params.AllowSharingFirstCluster) > mRecoParams[iteration].params.ClusterSharing) {
        continue;
      }

      // here we can do the calculation of the time bracket simply
      // by checkig in which rofs the clusters are
      int bcStart{0}, bcEnd{std::numeric_limits<int>::max()};
      for (int iLayer{0}; iLayer < mRecoParams[iteration].params.NLayers; ++iLayer) {
        if (track.getClusterIndex(iLayer) == constants::UnusedIndex) {
          continue;
        }
        mTimeFrame->markUsedCluster(iLayer, track.getClusterIndex(iLayer));
        int currentROF = mTimeFrame->getClusterROF(iLayer, track.getClusterIndex(iLayer));
        int bcClsSta = mTimeFrame->getROFOverlapTableView().getLayer(iLayer).getROFStartInBC(currentROF);
        int bcClsEnd = mTimeFrame->getROFOverlapTableView().getLayer(iLayer).getROFEndInBC(currentROF);
        bcStart = std::max(bcStart, bcClsSta);
        bcEnd = std::min(bcEnd, bcClsEnd);
      }
      track.getTimeStamp().setTimeStamp(bcStart);
      track.getTimeStamp().setTimeStampError(bcEnd - bcStart + 1);
      track.setUserField(0);
      track.getParamOut().setUserField(0);
      mTimeFrame->getTracks().emplace_back(track);

      firstClusters[firstLayer].push_back(firstCluster);
      if (isFirstShared) {
        sharedFirstClusters[firstLayer].push_back(firstCluster);
      }
    }
  }

  if (mRecoParams[iteration].params.AllowSharingFirstCluster) {
    /// Now we have to set the shared cluster flag
    for (int iLayer{0}; iLayer < mRecoParams[iteration].params.NLayers; ++iLayer) {
      std::sort(sharedFirstClusters[iLayer].begin(), sharedFirstClusters[iLayer].end());
    }
    for (auto& track : mTimeFrame->getTracks()) {
      int firstLayer{mRecoParams[iteration].params.NLayers}, firstCluster{constants::UnusedIndex};
      for (int iLayer{0}; iLayer < mRecoParams[iteration].params.NLayers; ++iLayer) {
        if (track.getClusterIndex(iLayer) == constants::UnusedIndex) {
          continue;
        }
        firstLayer = iLayer;
        firstCluster = track.getClusterIndex(iLayer);
        break;
      }
      if (std::binary_search(sharedFirstClusters[firstLayer].begin(), sharedFirstClusters[firstLayer].end(), firstCluster)) {
        track.setSharedClusters();
      }
    }
  }

  // remove roads and cells
  int nCells = mTimeFrame->getCells().size();
  for (int iLayer{0}; iLayer < nCells; ++iLayer) {
    deepVectorClear(mTimeFrame->getTrackletsLookupTable()[iLayer]);
    if (iLayer < nCells - 1) {
      deepVectorClear(mTimeFrame->getCellsLookupTable()[iLayer]);
      deepVectorClear(mTimeFrame->getCellsNeighbours()[iLayer]);
      deepVectorClear(mTimeFrame->getCellsNeighboursLUT()[iLayer]);
    }
    if (iLayer == 0 && mRecoParams[iteration].params.FindShortTracks) {
      continue;
    }
    deepVectorClear(mTimeFrame->getCells()[iLayer]);
    if ((mTimeFrame->hasMCinformation() && mRecoParams[iteration].params.createArtefactLabels) || OPTIMISATION_ANY) {
      deepVectorClear(mTimeFrame->getCellsLabel(iLayer));
    }
  }
  deepVectorClear(mTimeFrame->getRoads());
}

template <int NLayers>
void TrackerTraits<NLayers>::extendTracks(const int iteration)
{
  // TODO fix
  // for (auto& track : mTimeFrame->getTracks()) {
  //   auto backup{track};
  //   bool success{false};
  //   // the order here biases towards top extension, tracks should probably be fitted separately in the directions and then compared.
  //   if ((mRecoParams[iteration].params.UseTrackFollowerMix || mRecoParams[iteration].params.UseTrackFollowerTop) && track.getLastClusterLayer() != mRecoParams[iteration].params.NLayers - 1) {
  //     success = success || trackFollowing(&track, rof, true, iteration);
  //   }
  //   if ((mRecoParams[iteration].params.UseTrackFollowerMix || (mRecoParams[iteration].params.UseTrackFollowerBot && !success)) && track.getFirstClusterLayer() != 0) {
  //     success = success || trackFollowing(&track, rof, false, iteration);
  //   }
  //   if (success) {
  //     /// We have to refit the track
  //     track.resetCovariance();
  //     track.setChi2(0);
  //     bool fitSuccess = fitTrack(track, 0, mRecoParams[iteration].params.NLayers, 1, mRecoParams[iteration].params.MaxChi2ClusterAttachment, mRecoParams[0].params.MaxChi2NDF);
  //     if (!fitSuccess) {
  //       track = backup;
  //       continue;
  //     }
  //     track.getParamOut() = track;
  //     track.resetCovariance();
  //     track.setChi2(0);
  //     fitSuccess = fitTrack(track, mRecoParams[iteration].params.NLayers - 1, -1, -1, mRecoParams[iteration].params.MaxChi2ClusterAttachment, mRecoParams[0].params.MaxChi2NDF, 50.);
  //     if (!fitSuccess) {
  //       track = backup;
  //       continue;
  //     }
  //     mTimeFrame->mNExtendedTracks++;
  //     mTimeFrame->mNExtendedUsedClusters += track.getNClusters() - backup.getNClusters();
  //     auto pattern = track.getPattern();
  //     auto diff = (pattern & ~backup.getPattern()) & 0xff;
  //     pattern |= (diff << 24);
  //     track.setPattern(pattern);
  //     /// Make sure that the newly attached clusters get marked as used
  //     for (int iLayer{0}; iLayer < mRecoParams[iteration].params.NLayers; ++iLayer) {
  //       if (track.getClusterIndex(iLayer) == constants::UnusedIndex) {
  //         continue;
  //       }
  //       mTimeFrame->markUsedCluster(iLayer, track.getClusterIndex(iLayer));
  //     }
  //   }
  // }
}

template <int NLayers>
void TrackerTraits<NLayers>::findShortPrimaries(const int iteration)
{
  // TODO fix
  // const auto propagator = o2::base::Propagator::Instance();
  // mTimeFrame->fillPrimaryVerticesXandAlpha();
  //
  // for (auto& cell : mTimeFrame->getCells()[0]) {
  //   auto& cluster3_glo = mTimeFrame->getClusters()[2][cell.getThirdClusterIndex()];
  //   auto& cluster2_glo = mTimeFrame->getClusters()[1][cell.getSecondClusterIndex()];
  //   auto& cluster1_glo = mTimeFrame->getClusters()[0][cell.getFirstClusterIndex()];
  //   if (mTimeFrame->isClusterUsed(2, cluster1_glo.clusterId) ||
  //       mTimeFrame->isClusterUsed(1, cluster2_glo.clusterId) ||
  //       mTimeFrame->isClusterUsed(0, cluster3_glo.clusterId)) {
  //     continue;
  //   }
  //
  //   std::array<int, 3> rofs{
  //     mTimeFrame->getClusterROF(2, cluster3_glo.clusterId),
  //     mTimeFrame->getClusterROF(1, cluster2_glo.clusterId),
  //     mTimeFrame->getClusterROF(0, cluster1_glo.clusterId)};
  //   if (rofs[0] != rofs[1] && rofs[1] != rofs[2] && rofs[0] != rofs[2]) {
  //     continue;
  //   }
  //
  //   int rof{rofs[0]};
  //   if (rofs[1] == rofs[2]) {
  //     rof = rofs[2];
  //   }
  //
  //   auto pvs{mTimeFrame->getPrimaryVertices(0, rof)};
  //   auto pvsXAlpha{mTimeFrame->getPrimaryVerticesXAlpha(0, rof)};
  //
  //   const auto& cluster3_tf = mTimeFrame->getTrackingFrameInfoOnLayer(2)[cluster3_glo.clusterId];
  //   TrackITSExt temporaryTrack{buildTrackSeed(cluster1_glo, cluster2_glo, cluster3_tf)};
  //   temporaryTrack.setExternalClusterIndex(0, cluster1_glo.clusterId, true);
  //   temporaryTrack.setExternalClusterIndex(1, cluster2_glo.clusterId, true);
  //   temporaryTrack.setExternalClusterIndex(2, cluster3_glo.clusterId, true);
  //
  //   /// add propagation to the primary vertices compatible with the ROF(s) of the cell
  //   bool fitSuccess = fitTrack(temporaryTrack, 1, -1, -1);
  //   if (!fitSuccess) {
  //     continue;
  //   }
  //   fitSuccess = false;
  //
  //   TrackITSExt bestTrack{temporaryTrack}, backup{temporaryTrack};
  //   float bestChi2{std::numeric_limits<float>::max()};
  //   for (int iV{0}; iV < (int)pvs.size(); ++iV) {
  //     temporaryTrack = backup;
  //     if (!temporaryTrack.rotate(pvsXAlpha[iV][1])) {
  //       continue;
  //     }
  //     if (!propagator->propagateTo(temporaryTrack, pvsXAlpha[iV][0], true)) {
  //       continue;
  //     }
  //
  //     float pvRes{mRecoParams[0].params.PVres / o2::gpu::CAMath::Sqrt(float(pvs[iV].getNContributors()))};
  //     const float posVtx[2]{0.f, pvs[iV].getZ()};
  //     const float covVtx[3]{pvRes, 0.f, pvRes};
  //     float chi2 = temporaryTrack.getPredictedChi2Quiet(posVtx, covVtx);
  //     if (chi2 < bestChi2) {
  //       if (!temporaryTrack.track::TrackParCov::update(posVtx, covVtx)) {
  //         continue;
  //       }
  //       bestTrack = temporaryTrack;
  //       bestChi2 = chi2;
  //     }
  //   }
  //
  //   bestTrack.resetCovariance();
  //   bestTrack.setChi2(0.f);
  //   fitSuccess = fitTrack(bestTrack, 0, mRecoParams[iteration].params.NLayers, 1, mRecoParams[iteration].params.MaxChi2ClusterAttachment, mRecoParams[iteration].params.MaxChi2NDF);
  //   if (!fitSuccess) {
  //     continue;
  //   }
  //   bestTrack.getParamOut() = bestTrack;
  //   bestTrack.resetCovariance();
  //   bestTrack.setChi2(0.f);
  //   fitSuccess = fitTrack(bestTrack, mRecoParams[iteration].params.NLayers - 1, -1, -1, mRecoParams[iteration].params.MaxChi2ClusterAttachment, mRecoParams[iteration].params.MaxChi2NDF, 50.);
  //   if (!fitSuccess) {
  //     continue;
  //   }
  //   mTimeFrame->markUsedCluster(0, bestTrack.getClusterIndex(0));
  //   mTimeFrame->markUsedCluster(1, bestTrack.getClusterIndex(1));
  //   mTimeFrame->markUsedCluster(2, bestTrack.getClusterIndex(2));
  //   mTimeFrame->getTracks().emplace_back(bestTrack);
  // }
}

template <int NLayers>
bool TrackerTraits<NLayers>::fitTrack(TrackITSExt& track, int start, int end, int step, float chi2clcut, float chi2ndfcut, float maxQoverPt, int nCl, o2::track::TrackPar* linRef)
{
  auto propInstance = o2::base::Propagator::Instance();

  for (int iLayer{start}; iLayer != end; iLayer += step) {
    if (track.getClusterIndex(iLayer) == constants::UnusedIndex) {
      continue;
    }
    const TrackingFrameInfo& trackingHit = mTimeFrame->getTrackingFrameInfoOnLayer(iLayer)[track.getClusterIndex(iLayer)];
    if (linRef) {
      if (!track.rotate(trackingHit.alphaTrackingFrame, *linRef, getBz())) {
        return false;
      }
      if (!propInstance->propagateToX(track, *linRef, trackingHit.xTrackingFrame, getBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, mRecoParams[0].params.CorrType)) {
        return false;
      }
      if (mRecoParams[0].params.CorrType == o2::base::PropagatorF::MatCorrType::USEMatCorrNONE) {
        if (!track.correctForMaterial(*linRef, mRecoParams[0].params.LayerxX0[iLayer], mRecoParams[0].params.LayerxX0[iLayer] * constants::Radl * constants::Rho, true)) {
          continue;
        }
      }
    } else {
      if (!track.rotate(trackingHit.alphaTrackingFrame)) {
        return false;
      }
      if (!propInstance->propagateToX(track, trackingHit.xTrackingFrame, getBz(), o2::base::PropagatorImpl<float>::MAX_SIN_PHI, o2::base::PropagatorImpl<float>::MAX_STEP, mRecoParams[0].params.CorrType)) {
        return false;
      }
      if (mRecoParams[0].params.CorrType == o2::base::PropagatorF::MatCorrType::USEMatCorrNONE) {
        if (!track.correctForMaterial(mRecoParams[0].params.LayerxX0[iLayer], mRecoParams[0].params.LayerxX0[iLayer] * constants::Radl * constants::Rho, true)) {
          continue;
        }
      }
    }
    auto predChi2{track.getPredictedChi2Quiet(trackingHit.positionTrackingFrame, trackingHit.covarianceTrackingFrame)};
    if ((nCl >= 3 && predChi2 > chi2clcut) || predChi2 < 0.f) {
      return false;
    }
    track.setChi2(track.getChi2() + predChi2);
    if (!track.o2::track::TrackParCov::update(trackingHit.positionTrackingFrame, trackingHit.covarianceTrackingFrame)) {
      return false;
    }
    if (linRef && mRecoParams[0].params.ShiftRefToCluster) { // displace the reference to the last updated cluster
      linRef->setY(trackingHit.positionTrackingFrame[0]);
      linRef->setZ(trackingHit.positionTrackingFrame[1]);
    }
    nCl++;
  }
  return std::abs(track.getQ2Pt()) < maxQoverPt && track.getChi2() < chi2ndfcut * (nCl * 2 - 5);
}

template <int NLayers>
bool TrackerTraits<NLayers>::trackFollowing(TrackITSExt* track, int rof, bool outward, const int iteration)
{
  auto propInstance = o2::base::Propagator::Instance();
  const int step = -1 + (static_cast<int>(outward) * 2);
  const int end = outward ? mRecoParams[iteration].params.NLayers - 1 : 0;
  bounded_vector<TrackITSExt> hypotheses(1, *track, mMemoryPool.get()); // possibly avoid reallocation
  for (size_t iHypo{0}; iHypo < hypotheses.size(); ++iHypo) {
    auto hypo{hypotheses[iHypo]};
    int iLayer = static_cast<int>(outward ? hypo.getLastClusterLayer() : hypo.getFirstClusterLayer());
    // per layer we add new hypotheses
    while (iLayer != end) {
      iLayer += step; // step through all layers until we reach the end, this allows for skipping on empty layers
      const float r = mRecoParams[iteration].params.LayerRadii[iLayer];
      // get an estimate of the trackinf-frame x for the next step
      float x{-999};
      if (!hypo.getXatLabR(r, x, mTimeFrame->getBz(), o2::track::DirAuto) || x <= 0.f) {
        continue;
      }
      // estimate hypo's trk parameters at that x
      auto& hypoParam{outward ? hypo.getParamOut() : hypo.getParamIn()};
      if (!propInstance->propagateToX(hypoParam, x, mTimeFrame->getBz(), base::PropagatorF::MAX_SIN_PHI,
                                      base::PropagatorF::MAX_STEP, mRecoParams[iteration].params.CorrType)) {
        continue;
      }

      if (mRecoParams[iteration].params.CorrType == base::PropagatorF::MatCorrType::USEMatCorrNONE) { // account for material affects if propagator does not
        if (!hypoParam.correctForMaterial(mRecoParams[iteration].params.LayerxX0[iLayer], mRecoParams[iteration].params.LayerxX0[iLayer] * constants::Radl * constants::Rho, true)) {
          continue;
        }
      }

      // calculate the search window on this layer
      const float phi{hypoParam.getPhi()};
      const float ePhi{o2::gpu::CAMath::Sqrt(hypoParam.getSigmaSnp2() / hypoParam.getCsp2())};
      const float z{hypoParam.getZ()};
      const float eZ{o2::gpu::CAMath::Sqrt(hypoParam.getSigmaZ2())};
      const int4 selectedBinsRect{getBinsRect(iteration, iLayer, phi, mRecoParams[iteration].params.NSigmaCut * ePhi, z, mRecoParams[iteration].params.NSigmaCut * eZ)};
      if (selectedBinsRect.x == 0 && selectedBinsRect.y == 0 && selectedBinsRect.z == 0 && selectedBinsRect.w == 0) {
        continue;
      }

      int phiBinsNum{selectedBinsRect.w - selectedBinsRect.y + 1};

      if (phiBinsNum < 0) {
        phiBinsNum += mRecoParams[iteration].params.PhiBins;
      }

      gsl::span<const Cluster> layer1 = mTimeFrame->getClustersOnLayer(rof, iLayer);
      if (layer1.empty()) {
        continue;
      }

      // check all clusters in search windows for possible new hypotheses
      for (int iPhiCount = 0; iPhiCount < phiBinsNum; iPhiCount++) {
        int iPhiBin = (selectedBinsRect.y + iPhiCount) % mRecoParams[iteration].params.PhiBins;
        const int firstBinIndex{mTimeFrame->getIndexTableUtils().getBinIndex(selectedBinsRect.x, iPhiBin)};
        const int maxBinIndex{firstBinIndex + selectedBinsRect.z - selectedBinsRect.x + 1};
        const int firstRowClusterIndex = mTimeFrame->getIndexTable(rof, iLayer)[firstBinIndex];
        const int maxRowClusterIndex = mTimeFrame->getIndexTable(rof, iLayer)[maxBinIndex];

        for (int iNextCluster{firstRowClusterIndex}; iNextCluster < maxRowClusterIndex; ++iNextCluster) {
          if (iNextCluster >= (int)layer1.size()) {
            break;
          }
          const Cluster& nextCluster{layer1[iNextCluster]};

          if (mTimeFrame->isClusterUsed(iLayer, nextCluster.clusterId)) {
            continue;
          }

          const TrackingFrameInfo& trackingHit = mTimeFrame->getTrackingFrameInfoOnLayer(iLayer)[nextCluster.clusterId];

          auto tbupdated{hypo};
          auto& tbuParams = outward ? tbupdated.getParamOut() : tbupdated.getParamIn();
          if (!tbuParams.rotate(trackingHit.alphaTrackingFrame)) {
            continue;
          }

          if (!propInstance->propagateToX(tbuParams, trackingHit.xTrackingFrame, mTimeFrame->getBz(),
                                          base::PropagatorF::MAX_SIN_PHI, base::PropagatorF::MAX_STEP, base::PropagatorF::MatCorrType::USEMatCorrNONE)) {
            continue;
          }

          auto predChi2{tbuParams.getPredictedChi2Quiet(trackingHit.positionTrackingFrame, trackingHit.covarianceTrackingFrame)};
          if (predChi2 >= track->getChi2() * mRecoParams[iteration].params.NSigmaCut) {
            continue;
          }

          if (!tbuParams.o2::track::TrackParCov::update(trackingHit.positionTrackingFrame, trackingHit.covarianceTrackingFrame)) {
            continue;
          }
          tbupdated.setChi2(tbupdated.getChi2() + predChi2); /// This is wrong for outward propagation as the chi2 refers to inward parameters
          tbupdated.setExternalClusterIndex(iLayer, nextCluster.clusterId, true);
          hypotheses.emplace_back(tbupdated);
        }
      }
    }
  }

  TrackITSExt* bestHypo{track};
  bool swapped{false};
  for (auto& hypo : hypotheses) {
    if (hypo.isBetter(*bestHypo, track->getChi2() * mRecoParams[iteration].params.NSigmaCut)) {
      bestHypo = &hypo;
      swapped = true;
    }
  }
  *track = *bestHypo;
  return swapped;
}

// create a new seed either from the existing track inner param or reseed from the edgepointd and cluster in the middle
template <int nLayers>
TrackITSExt TrackerTraits<nLayers>::seedTrackForRefit(const CellSeedN& seed)
{
  TrackITSExt temporaryTrack(seed);
  int lrMin = nLayers, lrMax = 0, lrMid = 0;
  for (int iL = 0; iL < nLayers; ++iL) {
    const int idx = seed.getCluster(iL);
    temporaryTrack.setExternalClusterIndex(iL, idx, idx != constants::UnusedIndex);
    if (idx != constants::UnusedIndex) {
      lrMin = o2::gpu::CAMath::Min(lrMin, iL);
      lrMax = o2::gpu::CAMath::Max(lrMax, iL);
    }
  }
  int ncl = temporaryTrack.getNClusters();
  if (ncl < mRecoParams[0].params.ReseedIfShorter) { // reseed with circle passing via edges and the midpoint
    if (ncl == mRecoParams[0].params.NLayers) {
      lrMin = 0;
      lrMax = mRecoParams[0].params.NLayers - 1;
      lrMid = (lrMin + lrMax) / 2;
    } else {
      lrMid = lrMin + 1;
      float midR = 0.5 * (mRecoParams[0].params.LayerRadii[lrMax] + mRecoParams[0].params.LayerRadii[lrMin]), dstMidR = o2::gpu::GPUCommonMath::Abs(midR - mRecoParams[0].params.LayerRadii[lrMid]);
      for (int iL = lrMid + 1; iL < lrMax; ++iL) { // find the midpoint as closest to the midR
        auto dst = o2::gpu::GPUCommonMath::Abs(midR - mRecoParams[0].params.LayerRadii[iL]);
        if (dst < dstMidR) {
          lrMid = iL;
          dstMidR = dst;
        }
      }
    }
    const auto& cluster0_tf = mTimeFrame->getTrackingFrameInfoOnLayer(lrMin)[seed.getCluster(lrMin)]; // if the sensor frame!
    const auto& cluster1_gl = mTimeFrame->getUnsortedClusters()[lrMid][seed.getCluster(lrMid)];       // global frame
    const auto& cluster2_gl = mTimeFrame->getUnsortedClusters()[lrMax][seed.getCluster(lrMax)];       // global frame
    temporaryTrack.getParamIn() = buildTrackSeed(cluster2_gl, cluster1_gl, cluster0_tf, true);
  }
  temporaryTrack.resetCovariance();
  temporaryTrack.setCov(temporaryTrack.getQ2Pt() * temporaryTrack.getQ2Pt() * temporaryTrack.getCov()[o2::track::CovLabels::kSigQ2Pt2], o2::track::CovLabels::kSigQ2Pt2);
  return temporaryTrack;
}

/// Clusters are given from inside outward (cluster3 is the outermost). The outermost cluster is given in the tracking
/// frame coordinates whereas the others are referred to the global frame.
template <int nLayers>
track::TrackParCov TrackerTraits<nLayers>::buildTrackSeed(const Cluster& cluster1, const Cluster& cluster2, const TrackingFrameInfo& tf3, bool reverse)
{
  const float sign = reverse ? -1.f : 1.f;

  float ca = NAN, sa = NAN;
  o2::gpu::CAMath::SinCos(tf3.alphaTrackingFrame, sa, ca);

  const float x1 = (cluster1.xCoordinate * ca) + (cluster1.yCoordinate * sa);
  const float y1 = (-cluster1.xCoordinate * sa) + (cluster1.yCoordinate * ca);
  const float x2 = (cluster2.xCoordinate * ca) + (cluster2.yCoordinate * sa);
  const float y2 = (-cluster2.xCoordinate * sa) + (cluster2.yCoordinate * ca);
  const float x3 = tf3.xTrackingFrame;
  const float y3 = tf3.positionTrackingFrame[0];

  float snp = NAN, q2pt = NAN, q2pt2 = NAN;
  if (mIsZeroField) {
    const float tgp = o2::gpu::CAMath::ATan2(y3 - y1, x3 - x1);
    snp = sign * tgp / o2::gpu::CAMath::Sqrt(1.f + (tgp * tgp));
    q2pt = sign / track::kMostProbablePt;
    q2pt2 = 1.f;
  } else {
    const float crv = math_utils::computeCurvature(x3, y3, x2, y2, x1, y1);
    snp = sign * crv * (x3 - math_utils::computeCurvatureCentreX(x3, y3, x2, y2, x1, y1));
    q2pt = sign * crv / (mBz * o2::constants::math::B2C);
    q2pt2 = crv * crv;
  }

  const float tgl = 0.5f * (math_utils::computeTanDipAngle(x1, y1, x2, y2, cluster1.zCoordinate, cluster2.zCoordinate) +
                            math_utils::computeTanDipAngle(x2, y2, x3, y3, cluster2.zCoordinate, tf3.positionTrackingFrame[1]));
  const float sg2q2pt = track::kC1Pt2max * (q2pt2 > 0.0005f ? (q2pt2 < 1.f ? q2pt2 : 1.f) : 0.0005f);

  return {x3, tf3.alphaTrackingFrame, {y3, tf3.positionTrackingFrame[1], snp, tgl, q2pt}, {tf3.covarianceTrackingFrame[0], tf3.covarianceTrackingFrame[1], tf3.covarianceTrackingFrame[2], 0.f, 0.f, track::kCSnp2max, 0.f, 0.f, 0.f, track::kCTgl2max, 0.f, 0.f, 0.f, 0.f, sg2q2pt}};
}

template <int NLayers>
void TrackerTraits<NLayers>::computeTruthSeeding()
{
  LOGP(info, "Using truth seeds as vertices; will skip computations");
  const auto dc = o2::steer::DigitizationContext::loadFromFile("collisioncontext.root");
  const auto irs = dc->getEventRecords();
  // TODO in principle need to account for the bias which is anyways not well defined now
  int64_t roFrameBiasInBC{0};
  for (int iLayer{0}; iLayer < NLayers; ++iLayer) {
    roFrameBiasInBC += o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>::Instance().getROFBiasInBC(iLayer);
  }
  roFrameBiasInBC /= 7;
  o2::steer::MCKinematicsReader mcReader(dc);
  const int iSrc = 0; // take only events from collision generator
  auto eveId2colId = dc->getCollisionIndicesForSource(iSrc);
  int nVerts{0};
  for (int iEve{0}; iEve < mcReader.getNEvents(iSrc); ++iEve) {
    const auto& ir = irs[eveId2colId[iEve]];
    if (!ir.isDummy()) { // do we need this, is this for diffractive events?
      const auto& eve = mcReader.getMCEventHeader(iSrc, iEve);
      int64_t bc = ((ir - o2::raw::HBFUtils::Instance().getFirstSampledTFIR()).toLong() - roFrameBiasInBC);
      // find closest lower border in bcs
      Vertex vert;
      vert.getTimeStamp().setTimeStamp(bc);
      vert.getTimeStamp().setTimeStampError(1); // place it exactly at the BC
      // count potential contributors
      int nCont{0};
      // for (const auto& trk : mcReader.getTracks(iSrc, iEve)) {
      //   if (!trk.leftTrace(o2::detectors::DetID::ITS, eve.getDetId2HitBitLUT())) {
      //     continue;
      //   }
      //   if (!trk.isPrimary() || trk.GetPt() < 0.2 || std::abs(trk.GetEta()) > 1.1) {
      //     continue;
      //   }
      //   if (auto pdg = o2::O2DatabasePDG::Instance()->GetParticle(trk.GetPdgCode()); pdg && pdg->Charge() == 0) {
      //     continue;
      //   }
      //   ++nCont;
      // }
      // set minimum to 1 sometimes for diffractive events there is nothing acceptance
      vert.setNContributors(std::max(1, nCont));
      vert.setXYZ((float)eve.GetX(), (float)eve.GetY(), (float)eve.GetZ());
      vert.setChi2(1); // not used as constraint
      vert.setCov(25e-4, 25e-4, 25e-4, 25e-4, 125e-4, 100e-4);
      vert.print();
      o2::MCCompLabel lbl(o2::MCCompLabel::maxTrackID(), iEve, iSrc, false);
      mTimeFrame->addPrimaryVertex(vert);
      VertexLabel poll{lbl, 1.f};
      mTimeFrame->addPrimaryVertexLabel(poll);
      ++nVerts;
    }
    mcReader.releaseTracksForSourceAndEvent(iSrc, iEve);
  }
  sortSeeds();
  LOGP(info, "Imposed {} MC vertices", nVerts);
}

template <int NLayers>
void TrackerTraits<NLayers>::sortSeeds()
{
  // sort seeding vertices: 1. time, 2. error, 3. multiplicity
  // this by definition ensures that LUT is correct
  // and the timeslicing works as indented
  auto& seeds = mTimeFrame->getPrimaryVertices();
  auto& lbls = mTimeFrame->getPrimaryVerticesLabels();
  bounded_vector<size_t> indices(seeds.size(), mMemoryPool.get());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&seeds](const size_t a, const size_t b) {
    const auto& seedA = seeds[a];
    const auto& seedB = seeds[b];
    const auto& tA = seedA.getTimeStamp();
    const auto& tB = seedB.getTimeStamp();
    if (tA.getTimeStamp() != tB.getTimeStamp()) {
      return tA.getTimeStamp() < tB.getTimeStamp();
    }
    if (tA.getTimeStampError() != tB.getTimeStampError()) {
      return tA.getTimeStampError() < tB.getTimeStampError();
    }
    return seedA.getNContributors() < seedB.getNContributors();
  });
  bounded_vector<Vertex> seedsSorted(indices.size(), mMemoryPool.get());
  bounded_vector<VertexLabel> lblsSorted(mMemoryPool.get());
  if (mTimeFrame->hasMCinformation()) {
    lblsSorted.resize(indices.size());
  }
  for (size_t i{0}; i < indices.size(); ++i) {
    seedsSorted[i] = seeds[indices[i]];
    if (mTimeFrame->hasMCinformation()) {
      lblsSorted[i] = lbls[indices[i]];
    }
  }
  std::copy(seedsSorted.begin(), seedsSorted.end(), seeds.begin());
  if (mTimeFrame->hasMCinformation()) {
    std::copy(lblsSorted.begin(), lblsSorted.end(), lbls.begin());
  }
}

template <int NLayers>
void TrackerTraits<NLayers>::setBz(float bz)
{
  mBz = bz;
  mIsZeroField = std::abs(mBz) < 0.01;
  mTimeFrame->setBz(bz);
}

template <int NLayers>
bool TrackerTraits<NLayers>::isMatLUT() const
{
  return o2::base::Propagator::Instance()->getMatLUT() && (mRecoParams[0].params.CorrType == o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrLUT);
}

template <int NLayers>
void TrackerTraits<NLayers>::updateTrackingParameters(const std::vector<RecoIteration>& recoPars)
{
  mRecoParams = recoPars;
  if (static bool onceDone{false}; !onceDone) {
    onceDone = true;
    dbscan::DBSCANParams p{};
    p.minPts = mRecoParams[0].params.SeedingDBScanMinPt;
    p.eps[0] = mRecoParams[0].params.SeedingDBScanEpsZ;
    p.eps[1] = mRecoParams[0].params.SeedingDBScanEpsT;
    mDBScan = dbscan::DBSCAN(p);
  }
}

template <int NLayers>
void TrackerTraits<NLayers>::setNThreads(int n, std::shared_ptr<tbb::task_arena>& arena)
{
#if OPTIMISATION_ANY
  LOGP(info, "Debug output enabled; enforcing single-threaded");
  mTaskArena = std::make_shared<tbb::task_arena>(1);
#else
  if (arena == nullptr) {
    mTaskArena = std::make_shared<tbb::task_arena>(std::abs(n));
    LOGP(info, "Setting tracker with {} threads.", n);
  } else {
    mTaskArena = arena;
    LOGP(info, "Attaching tracker to calling thread's arena");
  }
#endif
}

/// Debug dumps
class MCStat
{
  struct Record {
    using MCMap = std::unordered_map<o2::MCCompLabel, uint32_t>;
    int32_t mIter{-1}, mLayer{-1};
    MCMap mMap;
    void print() const
    {
      uint32_t good{0};
      for (const auto& [lbl, s] : mMap) {
        if (lbl.isValid() && s > 0) {
          ++good;
        }
      }
      double frac = ((double)good / (double)mMap.size()) * 100.;
      LOGP(info, "\titer:{} layer:{} -> {}/{} (DS:{}) ({:.2f}%)", mIter, mLayer, good, mMap.size(), OPTIMISATION_DOWNSAMPLE, frac);
    }
  };
  Record* mCurRec{nullptr};
  std::vector<Record> mRecords;
  const char* mName;

 public:
  MCStat(const MCStat&) = default;
  MCStat(MCStat&&) = delete;
  MCStat& operator=(const MCStat&) = default;
  MCStat& operator=(MCStat&&) = delete;
  MCStat(const char* name) : mName(name) {}
  ~MCStat()
  {
    LOGP(info, "Stats for {}:", mName);
    for (const auto& rec : mRecords) {
      rec.print();
    }
  }
  void account(int iteration, int layer, const o2::MCCompLabel& lbl, bool acc)
  {
    if (!lbl.isValid()) {
      return;
    }
    if ((mCurRec == nullptr) || mCurRec->mIter != iteration || mCurRec->mLayer != layer) {
      bool found{false};
      for (auto& rec : mRecords) {
        if (rec.mIter == iteration && rec.mLayer == layer) {
          mCurRec = &rec;
          found = true;
        }
      }
      if (!found) {
        mRecords.emplace_back();
        mRecords.back().mIter = iteration;
        mRecords.back().mLayer = layer;
        mCurRec = &mRecords.back();
      }
    }
    auto& s = mCurRec->mMap[lbl];
    if (acc) {
      ++s;
    }
  }
};

template <int NLayers>
inline void TrackerTraits<NLayers>::debugComputeLayerTracklets(int iteration, int layer, const Cluster& currentCls, const Cluster& nextCls, const Vertex& pv, float sigmaZ, float sigmaPhi, bool accepted)
{
#if OPTIMISATION_NOT_SET(OPTIMISATION_TRACKLETS)
  return; // no-op
#endif

  static LogLogThrottler logger;
  LOG_IF(info, logger.needToLog(iteration, layer)) << "debugTree: LayerTracklets:" << iteration << ":" << layer << " (1:" << OPTIMISATION_DOWNSAMPLE << ") dumped entries " << logger.evCount;
  if (OPTIMISATION_DOWNSAMPLE > 1 && ((logger.evCount) % OPTIMISATION_DOWNSAMPLE) != 0) {
    return;
  }

  MCCompLabel label;
  if (mTimeFrame->hasMCinformation()) {
    int currentId{currentCls.clusterId};
    int nextId{nextCls.clusterId};
    for (auto& lab1 : mTimeFrame->getClusterLabels(layer, currentId)) {
      for (auto& lab2 : mTimeFrame->getClusterLabels(layer + 1, nextId)) {
        if (lab1 == lab2 && lab1.isValid()) {
          label = lab1;
          break;
        }
      }
      if (label.isValid()) {
        break;
      }
    }
  }

  const float deltaPvZ = currentCls.zCoordinate - pv.getZ();
  const float tanLambda = deltaPvZ / currentCls.radius;
  const float deltaZ = o2::gpu::GPUCommonMath::Abs((tanLambda * (nextCls.radius - currentCls.radius)) + currentCls.zCoordinate - nextCls.zCoordinate);
  const float tglNext = (nextCls.zCoordinate - pv.getZ()) / nextCls.radius;
  const float deltaZ2Pv = nextCls.zCoordinate + ((nextCls.zCoordinate - currentCls.zCoordinate) / (nextCls.radius - currentCls.radius) * (nextCls.radius - pv.getR()));
  const float deltaPhi = o2::gpu::CAMath::Abs(o2::math_utils::toPMPi(currentCls.phi - nextCls.phi));

  static MCStat stats("trackleting");
  stats.account(iteration, layer, label, accepted);

  (*sDBGOut) << "tracklets"
             << "iter=" << iteration
             << "lay=" << layer
             << "lbl=" << label
             << "pv=" << pv
             << "curCls=" << currentCls
             << "nextCls=" << nextCls
             << "sigZ=" << sigmaZ
             << "delZ=" << deltaZ
             << "sigPhi=" << sigmaPhi
             << "delPhi=" << deltaPhi
             << "tgl=" << tanLambda
             << "tglNext=" << tglNext
             << "delZPv=" << deltaZ2Pv
             << "acc=" << accepted
             << "\n";
}

template <int NLayers>
inline void TrackerTraits<NLayers>::debugComputeLayerCells(int iteration, int layer, int currentTrkl, int nextTrkl)
{
#if !(OPTIMISATION_SET(OPTIMISATION_CELLS))
  return; // no-op
#endif

  static LogLogThrottler logger;
  LOG_IF(info, logger.needToLog(iteration, layer)) << "debugTree: LayerCells:" << iteration << ":" << layer << " (1:" << OPTIMISATION_DOWNSAMPLE << ") dumped entries " << logger.evCount;
  if (OPTIMISATION_DOWNSAMPLE > 1 && ((logger.evCount) % OPTIMISATION_DOWNSAMPLE) != 0) {
    return;
  }

  const Tracklet& currentTracklet{mTimeFrame->getTracklets()[layer][currentTrkl]};
  const auto currentLbl = mTimeFrame->getTrackletsLabel(layer)[currentTrkl];
  const Tracklet& nextTracklet{mTimeFrame->getTracklets()[layer + 1][nextTrkl]};
  const auto nextLbl = mTimeFrame->getTrackletsLabel(layer + 1)[nextTrkl];
  const auto lbl = (currentLbl == nextLbl) ? currentLbl : o2::MCCompLabel();

  const float delTgl = std::abs(currentTracklet.tanLambda - nextTracklet.tanLambda);
  const float delPhi = o2::gpu::CAMath::Abs(o2::math_utils::toPMPi(nextTracklet.phi - currentTracklet.phi));

  const auto& cls1 = mTimeFrame->getClusters()[layer][currentTracklet.firstClusterIndex];
  const auto& cls2 = mTimeFrame->getClusters()[layer + 1][nextTracklet.firstClusterIndex];
  const auto& cls3 = mTimeFrame->getClusters()[layer + 2][nextTracklet.secondClusterIndex];
  const float k123 = math_utils::computeCurvature(cls1.xCoordinate, cls1.yCoordinate,
                                                  cls2.xCoordinate, cls2.yCoordinate,
                                                  cls3.xCoordinate, cls3.yCoordinate);
  const float k013 = math_utils::computeCurvature(0.f, 0.f,
                                                  cls1.xCoordinate, cls1.yCoordinate,
                                                  cls3.xCoordinate, cls3.yCoordinate);
  const float dk = o2::gpu::CAMath::Abs(k123 - k013);
  const float dR2 = 0.25f * (math_utils::Sq(cls1.xCoordinate - cls3.xCoordinate) + math_utils::Sq(cls1.yCoordinate - cls3.yCoordinate));
  const float tgl2 = math_utils::Sq(0.5 * (currentTracklet.tanLambda + nextTracklet.tanLambda));
  const float snl2 = tgl2 / (1.f + tgl2);
  const float resMid2 = math_utils::Sq(mTimeFrame->getPositionResolution(layer + 1));
  const float msMid2 = math_utils::Sq(mTimeFrame->getMSangle(layer + 1));
  const float tglNSigma = o2::gpu::CAMath::Sqrt(resMid2 + msMid2) * mRecoParams[iteration].params.NSigmaCut;
  const float kSigma2Pos = 6.f * resMid2 / math_utils::Sq(dR2);
  const float kSigma2MS = msMid2 / (dR2 * snl2);
  const float kNSigma = mRecoParams[iteration].params.NSigmaCut * o2::gpu::CAMath::Sqrt(kSigma2Pos + kSigma2MS);

  bool accepted = dk < kNSigma && delTgl < tglNSigma;

  static MCStat stats("celling");
  stats.account(iteration, layer, lbl, accepted);

  (*sDBGOut) << "cells"
             << "iter=" << iteration
             << "lay=" << layer
             << "lbl=" << lbl
             << "curLbl=" << currentLbl
             << "nextLbl=" << nextLbl
             << "delTgl=" << delTgl
             << "delPhi=" << delPhi
             << "curTrkl=" << currentTracklet
             << "nextTrkl=" << nextTracklet
             << "k123=" << k123
             << "k013=" << k013
             << "dk=" << dk
             << "kNSigma=" << kNSigma
             << "tglNSigma=" << tglNSigma
             << "acc=" << accepted
             << "\n";
}

template <int NLayers>
inline void TrackerTraits<NLayers>::debugFindCellsNeighbours(int iteration, int layer, int currentCell, int nextCell, bool accepted)
{
#if OPTIMISATION_NOT_SET(OPTIMISATION_CELLSNEIGH)
  return; // no-op
#endif

  static LogLogThrottler logger;
  LOG_IF(info, logger.needToLog(iteration, layer)) << "debugTree: CellsNeighbours:" << iteration << ":" << layer << " (1:" << OPTIMISATION_DOWNSAMPLE << ") dumped entries " << logger.evCount;
  if (OPTIMISATION_DOWNSAMPLE > 1 && ((logger.evCount) % OPTIMISATION_DOWNSAMPLE) != 0) {
    return;
  }

  const auto& currentCellSeed{mTimeFrame->getCells()[layer][currentCell]};
  bool good{true};
  auto nextCellSeed{mTimeFrame->getCells()[layer + 1][nextCell]};
  if (!nextCellSeed.rotate(currentCellSeed.getAlpha()) ||
      !nextCellSeed.propagateTo(currentCellSeed.getX(), getBz())) {
    good = false;
  }
  const auto& currentLbl = mTimeFrame->getCellsLabel(layer)[currentCell];
  const auto& nextLbl = mTimeFrame->getCellsLabel(layer + 1)[nextCell];
  const auto lbl = (currentLbl == nextLbl) ? currentLbl : o2::MCCompLabel();
  const auto chi2 = currentCellSeed.getPredictedChi2(nextCellSeed);

  (*sDBGOut) << "cneigh"
             << "iter=" << iteration
             << "lay=" << layer
             << "lbl=" << lbl
             << "curLbl=" << currentLbl
             << "nextLbl=" << nextLbl
             << "chi2=" << chi2
             << "curCell=" << currentCellSeed
             << "nextCell=" << nextCellSeed
             << "prop=" << good
             << "acc=" << accepted
             << "\n";
}

// explicitly instaniate the ITS2/ITS3 tracker functions
template class TrackerTraits<7>;

} // namespace o2::its
