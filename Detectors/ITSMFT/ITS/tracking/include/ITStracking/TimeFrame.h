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

#ifndef TRACKINGITSU_INCLUDE_TIMEFRAME_H_
#define TRACKINGITSU_INCLUDE_TIMEFRAME_H_

#include <array>
#include <cstdint>
#include <vector>
#include <utility>
#include <algorithm>
#include <numeric>
#include <gsl/gsl>

#include "ITStracking/Cell.h"
#include "ITStracking/Cluster.h"
#include "ITStracking/Configuration.h"
#include "ITStracking/Definitions.h"
#include "ITStracking/Road.h"
#include "ITStracking/Tracklet.h"
#include "ITStracking/IndexTableUtils.h"
#include "ITStracking/ExternalAllocator.h"
#include "ITStracking/BoundedAllocator.h"
#include "ITStracking/ROFLookupTables.h"
#include "DataFormatsITS/TrackITS.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "DetectorsBase/Propagator.h"
#include "DataFormatsCalibration/MeanVertexObject.h"

namespace o2
{
namespace gpu
{
class GPUChainITS;
}

namespace itsmft
{
class CompClusterExt;
class TopologyDictionary;
class ROFRecord;
} // namespace itsmft

namespace its
{
namespace gpu
{
template <int>
class TimeFrameGPU;
}

template <int NLayers = 7>
struct TimeFrame {
  using IndexTableUtilsN = IndexTableUtils<NLayers>;
  using ROFOverlapTableN = ROFOverlapTable<NLayers>;
  using ROFVertexLookupTableN = ROFVertexLookupTable<NLayers>;
  using ROFTimeSliceTableN = ROFTimeSliceTable<NLayers>;
  using CellSeedN = CellSeed<NLayers>;
  friend class gpu::TimeFrameGPU<NLayers>;

  TimeFrame() = default;
  TimeFrame(const TimeFrame&) = delete;
  TimeFrame(TimeFrame&&) = delete;
  TimeFrame& operator=(const TimeFrame&) = delete;
  TimeFrame& operator=(TimeFrame&&) = delete;
  virtual ~TimeFrame() = default;

  void initialise(const RecoIteration& reco);

  const Vertex& getPrimaryVertex(const int ivtx) const { return mPrimaryVertices[ivtx]; }
  auto& getPrimaryVertices() { return mPrimaryVertices; };
  auto getPrimaryVerticesNum() { return mPrimaryVertices.size(); };
  const auto& getPrimaryVertices() const { return mPrimaryVertices; };
  auto& getPrimaryVerticesLabels() { return mPrimaryVerticesLabels; };
  gsl::span<const Vertex> getPrimaryVertices(int layer, int rofId) const;
  gsl::span<const std::array<float, 2>> getPrimaryVerticesXAlpha(int layer, int rofId) const;
  void fillPrimaryVerticesXandAlpha();
  void addPrimaryVertex(const Vertex& vertex) { mPrimaryVertices.emplace_back(vertex); }
  void addPrimaryVertexLabel(const VertexLabel& label) { mPrimaryVerticesLabels.push_back(label); }

  int loadROFrameData(const o2::itsmft::ROFRecord& rof, gsl::span<const itsmft::Cluster> clusters,
                      const dataformats::MCTruthContainer<MCCompLabel>* mcLabels = nullptr);

  void loadROFrameData(gsl::span<const o2::itsmft::ROFRecord> rofs,
                       gsl::span<const itsmft::CompClusterExt> clusters,
                       gsl::span<const unsigned char>::iterator& pattIt,
                       const itsmft::TopologyDictionary* dict,
                       int layer,
                       const dataformats::MCTruthContainer<MCCompLabel>* mcLabels = nullptr);
  void resetROFrameData(int iLayer);
  void prepareROFrameData(gsl::span<const itsmft::CompClusterExt> clusters, int layer);

  int getTotalClusters() const;
  int getTotalClustersPerROFrange(int rofMin, int range, int layerId) const;
  int getSortedIndex(int rofId, int layer, int idx) const { return mROFramesClusters[layer][rofId] + idx; }
  int getSortedStartIndex(const int rofId, const int layer) const { return mROFramesClusters[layer][rofId]; }
  int getNrof(int layer) const { return mROFramesClusters[layer].size() - 1; }

  auto& getMinRs() { return mMinR; }
  auto& getMaxRs() { return mMaxR; }
  float getMinR(int layer) const { return mMinR[layer]; }
  float getMaxR(int layer) const { return mMaxR[layer]; }
  float getMSangle(int layer) const { return mMSangles[layer]; }
  auto& getMSangles() { return mMSangles; }
  float getPhiCut(int layer) const { return mPhiCuts[layer]; }
  auto& getPhiCuts() { return mPhiCuts; }
  float getPositionResolution(int layer) const { return mPositionResolution[layer]; }
  auto& getPositionResolutions() { return mPositionResolution; }

  // seeding vertex constraint or current best estimate
  void setMeanVertex(const dataformats::MeanVertexObject* mv, float extraErr2 = 0.f);
  const dataformats::MeanVertexObject* getMeanVertexConstraint() const { return mMeanVertex; }
  auto& getMeanVertexRolling() { return mRollingVertex; }
  const dataformats::VertexBase& getMeanVertex() const { return (hasMeanVertex()) ? mMeanVertex->getMeanVertex() : mRollingVertex; }
  const auto& getMeanVertexInvErr() const { return mMeanVertexXYInvErr; }
  bool hasMeanVertex() const noexcept { return mMeanVertex != nullptr; }

  auto& getClusterLabelsContainer() { return mClusterLabels; }
  gsl::span<Cluster> getClustersOnLayer(int rofId, int layerId);
  gsl::span<const Cluster> getClustersOnLayer(int rofId, int layerId) const;
  gsl::span<const Cluster> getClustersPerROFrange(int rofMin, int range, int layerId) const;
  gsl::span<const Cluster> getUnsortedClustersOnLayer(int rofId, int layerId) const;
  gsl::span<uint8_t> getUsedClustersROF(int rofId, int layerId);
  gsl::span<const uint8_t> getUsedClustersROF(int rofId, int layerId) const;
  auto& getROFrameClusters(int layerId, int rofId) { return mROFramesClusters[layerId][rofId]; }
  gsl::span<const int> getROFrameClusters(int layerId) const;
  gsl::span<const int> getROFramesClustersPerROFrange(int rofMin, int range, int layerId) const;
  gsl::span<const int> getNClustersROFrange(int rofMin, int range, int layerId) const;
  gsl::span<const int> getIndexTablePerROFrange(int rofMin, int range, int layerId) const;
  gsl::span<int> getIndexTable(int rofId, int layerId);
  const auto& getTrackingFrameInfoOnLayer(int layerId) const { return mTrackingFrameInfo[layerId]; }

  // navigation tables
  const auto& getIndexTableUtils() const { return mIndexTableUtils; }
  const auto& getROFOverlapTable() const { return mROFOverlapTable; }
  const auto& getROFOverlapTableView() const { return mROFOverlapTableView; }
  void setROFOverlapTable(ROFOverlapTableN& table)
  {
    mROFOverlapTable = std::move(table);
    mROFOverlapTableView = mROFOverlapTable.getView();
  }
  const auto& getROFVertexLookupTable() const { return mROFVertexLookupTable; }
  const auto& getROFVertexLookupTableView() const { return mROFVertexLookupTableView; }
  void setROFVertexLookupTable(ROFVertexLookupTableN& table)
  {
    mROFVertexLookupTable = std::move(table);
    mROFVertexLookupTableView = mROFVertexLookupTable.getView();
  }
  void updateROFVertexLookupTable() { mROFVertexLookupTable.update(mPrimaryVertices.data(), mPrimaryVertices.size()); }
  const auto& getROFTimeSliceTable() const { return mROFTimeSliceTable; }
  const auto& getROFTimeSliceTableView() const { return mROFTimeSliceTableView; }
  void setROFTimeSliceTable(ROFTimeSliceTableN& table)
  {
    mROFTimeSliceTable = std::move(table);
    mROFTimeSliceTableView = mROFTimeSliceTable.getView();
  }

  // cluster information
  const TrackingFrameInfo& getClusterTrackingFrameInfo(int layerId, const Cluster& cl) const;
  gsl::span<const MCCompLabel> getClusterLabels(int layerId, const Cluster& cl) const { return getClusterLabels(layerId, cl.clusterId); }
  gsl::span<const MCCompLabel> getClusterLabels(int layerId, const int clId) const { return mClusterLabels[layerId]->getLabels(mClusterExternalIndices[layerId][clId]); }
  int getClusterExternalIndex(int layerId, const int clId) const { return mClusterExternalIndices[layerId][clId]; }
  int getClusterSize(int layer, int clusterId) const { return mClusterSize[layer][clusterId]; }
  void setClusterSize(int layer, bounded_vector<uint8_t>& v) { mClusterSize[layer] = std::move(v); }
  bool isClusterUsed(int layer, int clusterId) const { return mUsedClusters[layer][clusterId]; }
  void markUsedCluster(int layer, int clusterId) { mUsedClusters[layer][clusterId] = true; }
  gsl::span<unsigned char> getUsedClusters(const int layer);
  auto& getClusters() { return mClusters; }
  auto& getUnsortedClusters() { return mUnsortedClusters; }
  int getClusterROF(int iLayer, int iCluster);
  int getNumberOfClusters(int layer = -1) const;
  size_t getNumberOfUsedClusters() const;
  template <typename... T>
  void addClusterToLayer(int layer, T&&... args);
  template <typename... T>
  void addTrackingFrameInfoToLayer(int layer, T&&... args);
  void addClusterExternalIndexToLayer(int layer, const int idx) { mClusterExternalIndices[layer].push_back(idx); }

  // mc information
  bool hasMCinformation() const { return mClusterLabels[0] != nullptr; }

  // tracklet information
  auto& getTracklets() { return mTracklets; }
  auto& getTrackletsLabel(int layer) { return mTrackletLabels[layer]; }
  auto& getTrackletsLookupTable() { return mTrackletsLookupTable; }
  virtual int getNumberOfTracklets() const;

  // cells information
  auto& getCells() { return mCells; }
  auto& getCellsLabel(int layer) { return mCellLabels[layer]; }
  auto& getCellsLookupTable() { return mCellsLookupTable; }
  auto& getCellsNeighbours() { return mCellsNeighbours; }
  auto& getCellsNeighboursLUT() { return mCellsNeighboursLUT; }
  virtual int getNumberOfCells() const;

  // roads information
  auto& getRoads() { return mRoads; }
  virtual int getNumberOfNeighbours() const;

  // tracks information
  auto& getTracks() { return mTracks; }
  auto& getTracksLabel() { return mTracksLabel; }
  size_t getNumberOfTracks() const noexcept { return mTracks.size(); };
  auto getNumberOfExtendedTracks() const { return mNExtendedTracks; }
  auto getNumberOfUsedExtendedClusters() const { return mNExtendedUsedClusters; }

  /// memory management
  virtual void wipe();
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool);
  auto& getMemoryPool() const noexcept { return mMemoryPool; }
  bool checkMemory(unsigned long max) { return getArtefactsMemory() < max; }
  unsigned long getArtefactsMemory() const;
  void printArtefactsMemory() const;

  int hasBogusClusters() const
  {
    return std::accumulate(mBogusClusters.begin(), mBogusClusters.end(), 0);
  }

  /// State if memory will be externally managed by the GPU framework
  ExternalAllocator* mExternalAllocator{nullptr};
  std::shared_ptr<BoundedMemoryResource> mExtMemoryPool; // host memory pool managed by the framework
  auto getFrameworkAllocator() { return mExternalAllocator; };
  void setFrameworkAllocator(ExternalAllocator* ext);
  bool hasFrameworkAllocator() const noexcept { return mExternalAllocator != nullptr; }
  std::pmr::memory_resource* getMaybeFrameworkHostResource() { return hasFrameworkAllocator() ? mExtMemoryPool.get() : mMemoryPool.get(); }

  // magnetic field
  void setBz(float bz) { mBz = bz; }
  float getBz() const { return mBz; }

  // Propagator
  const o2::base::PropagatorImpl<float>* getDevicePropagator() const { return mPropagatorDevice; }
  virtual void setDevicePropagator(const o2::base::PropagatorImpl<float>* /*unused*/) {};

  // interface
  virtual bool isGPU() const noexcept { return false; }
  virtual const char* getName() const noexcept { return "CPU"; }

 private:
  void prepareClusters(const TrackingParameters& trkParam, const int maxLayers = NLayers);

  // tf data
  std::array<bounded_vector<Cluster>, NLayers> mClusters;
  std::array<bounded_vector<TrackingFrameInfo>, NLayers> mTrackingFrameInfo;
  std::array<bounded_vector<int>, NLayers> mClusterExternalIndices;
  std::array<bounded_vector<int>, NLayers> mROFramesClusters;
  std::array<const dataformats::MCTruthContainer<MCCompLabel>*, NLayers> mClusterLabels{nullptr};
  std::array<bounded_vector<int>, NLayers> mNClustersPerROF;
  std::array<bounded_vector<int>, NLayers> mIndexTables;
  std::vector<bounded_vector<int>> mTrackletsLookupTable;
  std::array<bounded_vector<uint8_t>, NLayers> mUsedClusters;
  int mNExtendedTracks{0};
  int mNExtendedUsedClusters{0};
  bounded_vector<Vertex> mPrimaryVertices;
  bounded_vector<VertexLabel> mPrimaryVerticesLabels;

  std::array<bounded_vector<Cluster>, NLayers> mUnsortedClusters;
  std::vector<bounded_vector<Tracklet>> mTracklets;
  std::vector<bounded_vector<CellSeedN>> mCells;
  bounded_vector<Road<NLayers - 2>> mRoads;
  bounded_vector<TrackITSExt> mTracks;
  bounded_vector<MCCompLabel> mTracksLabel;
  std::vector<bounded_vector<int>> mCellsNeighbours;
  std::vector<bounded_vector<int>> mCellsLookupTable;

  const o2::base::PropagatorImpl<float>* mPropagatorDevice = nullptr; // Needed only for GPU
  float mBz = 999.;
  const dataformats::MeanVertexObject* mMeanVertex{nullptr};
  dataformats::VertexBase mRollingVertex;
  std::array<float, 3> mMeanVertexXYInvErr{};
  std::array<float, NLayers> mMinR;
  std::array<float, NLayers> mMaxR;
  std::array<float, NLayers> mMSangles;
  std::array<float, NLayers - 1> mPhiCuts;
  std::array<float, NLayers> mPositionResolution;
  std::array<bounded_vector<uint8_t>, NLayers> mClusterSize;

  bounded_vector<std::array<float, 2>> mPValphaX; /// PV x and alpha for track propagation
  std::vector<bounded_vector<MCCompLabel>> mTrackletLabels;
  std::vector<bounded_vector<MCCompLabel>> mCellLabels;
  std::vector<bounded_vector<int>> mCellsNeighboursLUT;
  std::array<uint32_t, NLayers> mBogusClusters; /// keep track of clusters with wild coordinates

  // lookup tables
  IndexTableUtilsN mIndexTableUtils;
  ROFOverlapTableN mROFOverlapTable;
  ROFOverlapTableN::View mROFOverlapTableView;
  ROFVertexLookupTableN mROFVertexLookupTable;
  ROFVertexLookupTableN::View mROFVertexLookupTableView;
  ROFTimeSliceTableN mROFTimeSliceTable;
  ROFTimeSliceTableN::View mROFTimeSliceTableView;

  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
};

template <int NLayers>
inline gsl::span<const Vertex> TimeFrame<NLayers>::getPrimaryVertices(int layer, int rofId) const
{
  const auto& pvs = mROFVertexLookupTableView.getVertices(layer, rofId);
  return {&mPrimaryVertices[pvs.getFirstEntry()], static_cast<gsl::span<const Vertex>::size_type>(pvs.getEntries())};
}

template <int NLayers>
inline gsl::span<const std::array<float, 2>> TimeFrame<NLayers>::getPrimaryVerticesXAlpha(int layer, int rofId) const
{
  const auto& pvs = mROFVertexLookupTableView.getVertices(layer, rofId);
  return {&(mPValphaX[pvs.getFirstEntry()]), static_cast<gsl::span<const std::array<float, 2>>::size_type>(pvs.getEntries())};
}

template <int NLayers>
inline gsl::span<const int> TimeFrame<NLayers>::getROFrameClusters(int layerId) const
{
  return {&mROFramesClusters[layerId][0], static_cast<gsl::span<const int>::size_type>(mROFramesClusters[layerId].size())};
}

template <int NLayers>
inline gsl::span<Cluster> TimeFrame<NLayers>::getClustersOnLayer(int rofId, int layerId)
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mClusters[layerId][startIdx], static_cast<gsl::span<Cluster>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<const Cluster> TimeFrame<NLayers>::getClustersOnLayer(int rofId, int layerId) const
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mClusters[layerId][startIdx], static_cast<gsl::span<const Cluster>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<uint8_t> TimeFrame<NLayers>::getUsedClustersROF(int rofId, int layerId)
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mUsedClusters[layerId][startIdx], static_cast<gsl::span<uint8_t>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<const uint8_t> TimeFrame<NLayers>::getUsedClustersROF(int rofId, int layerId) const
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mUsedClusters[layerId][startIdx], static_cast<gsl::span<const uint8_t>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<const Cluster> TimeFrame<NLayers>::getClustersPerROFrange(int rofMin, int range, int layerId) const
{
  if (rofMin < 0 || rofMin >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofMin]}; // First cluster of rofMin
  int endIdx{mROFramesClusters[layerId][o2::gpu::CAMath::Min(rofMin + range, getNrof(layerId))]};
  return {&mClusters[layerId][startIdx], static_cast<gsl::span<Cluster>::size_type>(endIdx - startIdx)};
}

template <int NLayers>
inline gsl::span<const int> TimeFrame<NLayers>::getROFramesClustersPerROFrange(int rofMin, int range, int layerId) const
{
  int chkdRange{o2::gpu::CAMath::Min(range, getNrof(layerId) - rofMin)};
  return {&mROFramesClusters[layerId][rofMin], static_cast<gsl::span<int>::size_type>(chkdRange)};
}

template <int NLayers>
inline gsl::span<const int> TimeFrame<NLayers>::getNClustersROFrange(int rofMin, int range, int layerId) const
{
  int chkdRange{o2::gpu::CAMath::Min(range, getNrof(layerId) - rofMin)};
  return {&mNClustersPerROF[layerId][rofMin], static_cast<gsl::span<int>::size_type>(chkdRange)};
}

template <int NLayers>
inline int TimeFrame<NLayers>::getTotalClustersPerROFrange(int rofMin, int range, int layerId) const
{
  int startIdx{rofMin}; // First cluster of rofMin
  int endIdx{o2::gpu::CAMath::Min(rofMin + range, getNrof(layerId))};
  return mROFramesClusters[layerId][endIdx] - mROFramesClusters[layerId][startIdx];
}

template <int NLayers>
inline gsl::span<const int> TimeFrame<NLayers>::getIndexTablePerROFrange(int rofMin, int range, int layerId) const
{
  const int iTableSize{mIndexTableUtils.getNphiBins() * mIndexTableUtils.getNzBins() + 1};
  int chkdRange{o2::gpu::CAMath::Min(range, getNrof(layerId) - rofMin)};
  return {&mIndexTables[layerId][rofMin * iTableSize], static_cast<gsl::span<int>::size_type>(chkdRange * iTableSize)};
}

template <int NLayers>
inline int TimeFrame<NLayers>::getClusterROF(int iLayer, int iCluster)
{
  return std::lower_bound(mROFramesClusters[iLayer].begin(), mROFramesClusters[iLayer].end(), iCluster + 1) - mROFramesClusters[iLayer].begin() - 1;
}

template <int NLayers>
inline gsl::span<const Cluster> TimeFrame<NLayers>::getUnsortedClustersOnLayer(int rofId, int layerId) const
{
  if (rofId < 0 || rofId >= getNrof(layerId)) {
    return {};
  }
  int startIdx{mROFramesClusters[layerId][rofId]};
  return {&mUnsortedClusters[layerId][startIdx], static_cast<gsl::span<Cluster>::size_type>(mROFramesClusters[layerId][rofId + 1] - startIdx)};
}

template <int NLayers>
inline gsl::span<int> TimeFrame<NLayers>::getIndexTable(int rofId, int layer)
{
  if (rofId < 0 || rofId >= getNrof(layer)) {
    return {};
  }
  const int tableSize = mIndexTableUtils.getNphiBins() * mIndexTableUtils.getNzBins() + 1;
  return {&mIndexTables[layer][rofId * tableSize], static_cast<gsl::span<int>::size_type>(tableSize)};
}

template <int NLayers>
template <typename... T>
void TimeFrame<NLayers>::addClusterToLayer(int layer, T&&... values)
{
  mUnsortedClusters[layer].emplace_back(std::forward<T>(values)...);
}

template <int NLayers>
template <typename... T>
void TimeFrame<NLayers>::addTrackingFrameInfoToLayer(int layer, T&&... values)
{
  mTrackingFrameInfo[layer].emplace_back(std::forward<T>(values)...);
}

template <int NLayers>
inline gsl::span<uint8_t> TimeFrame<NLayers>::getUsedClusters(const int layer)
{
  return {&mUsedClusters[layer][0], static_cast<gsl::span<uint8_t>::size_type>(mUsedClusters[layer].size())};
}

template <int NLayers>
inline int TimeFrame<NLayers>::getTotalClusters() const
{
  size_t totalClusters{0};
  for (const auto& clusters : mUnsortedClusters) {
    totalClusters += clusters.size();
  }
  return int(totalClusters);
}

template <int NLayers>
inline int TimeFrame<NLayers>::getNumberOfClusters(int layer) const
{
  if (layer >= 0) {
    return mClusters[layer].size();
  }
  int nClusters = 0;
  for (const auto& layer : mClusters) {
    nClusters += layer.size();
  }
  return nClusters;
}

template <int NLayers>
inline int TimeFrame<NLayers>::getNumberOfCells() const
{
  int nCells = 0;
  for (const auto& layer : mCells) {
    nCells += layer.size();
  }
  return nCells;
}

template <int NLayers>
inline int TimeFrame<NLayers>::getNumberOfTracklets() const
{
  int nTracklets = 0;
  for (const auto& layer : mTracklets) {
    nTracklets += layer.size();
  }
  return nTracklets;
}

template <int NLayers>
inline int TimeFrame<NLayers>::getNumberOfNeighbours() const
{
  int n{0};
  for (const auto& l : mCellsNeighbours) {
    n += l.size();
  }
  return n;
}

template <int NLayers>
inline size_t TimeFrame<NLayers>::getNumberOfUsedClusters() const
{
  size_t nClusters = 0;
  for (const auto& layer : mUsedClusters) {
    nClusters += std::count(layer.begin(), layer.end(), true);
  }
  return nClusters;
}

} // namespace its
} // namespace o2

#endif
