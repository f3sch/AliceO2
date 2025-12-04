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

#ifndef TRACKINGITSU_INCLUDE_ROFOVERLAPTABLE_H_
#define TRACKINGITSU_INCLUDE_ROFOVERLAPTABLE_H_

#include <cstddef>
#include <cstdint>
#include <ranges>
#include <limits>
#include <numeric>
#include <vector>

#include "CommonConstants/LHCConstants.h"
#include "CommonDataFormat/RangeReference.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "ITStracking/Definitions.h"
#include "GPUCommonLogger.h"
#include "GPUCommonMath.h"
#include "GPUCommonDef.h"

namespace o2::its
{

// Layer timing definition
struct LayerTiming {
  using BCType = int32_t;
  using ROFRange = dataformats::RangeReference<BCType, BCType>;
  BCType mNROFsTF{0};   // number of ROFs per timeframe
  BCType mROFLength{0}; // ROF length in BC
  BCType mROFDelay{0};  // delay of ROFs wrt LHC orbit
  BCType mROFDelta{0};  // added delta in BC for compatibility

  // return start of ROF in BC
  GPUhdi() BCType getROFStartInBC(BCType rofId, bool actual = false) const noexcept
  {
    assert(rofId < mNROFsTF);
    return (mROFLength * rofId) + ((actual) ? mROFDelay : 0);
  }

  // return end of ROF in BCs
  GPUhdi() BCType getROFEndInBC(BCType rofId, bool actual = false) const noexcept
  {
    assert(rofId < mNROFsTF);
    return getROFStartInBC(rofId, actual) + mROFLength;
  }

  // return time-interval of rof [start, end)
  GPUhdi() ROFRange getROFTimeBounds(BCType rofId, bool actual = false) const noexcept
  {
    auto start = getROFStartInBC(rofId, actual);
    auto end = getROFEndInBC(rofId, actual);
    if (actual) {
      start -= mROFDelta;
      end += mROFDelta;
    }
    return {start, end - start + 1};
  }

  GPUhdi() BCType getROFIdForBC(BCType bc) const noexcept
  {
    const BCType rofId = (bc - mROFDelay) / mROFLength;
    assert(rofId >= 0);
    return rofId;
  }

  GPUhi() o2::InteractionRecord getIRForROFId(BCType rofId) const noexcept
  {
    return o2::InteractionRecord::long2IR(getROFStartInBC(rofId));
  }

  GPUh() std::string asString() const
  {
    return std::format("NROFsPerTF {:4} ROFLength {:4} ({:4} per Orbit) ROFDelay {:4} ROFDelta {:4}", mNROFsTF, mROFLength, (o2::constants::lhc::LHCMaxBunches / mROFLength), mROFDelay, mROFDelta);
  }

  GPUh() void print() const
  {
    LOG(info) << asString();
  }
};

// Base class for lookup to define layers
template <int32_t NLayers>
class LayerTimingBase
{
 protected:
  LayerTiming mLayers[NLayers];

 public:
  using T = LayerTiming::BCType;
  GPUdDefault() LayerTimingBase() = default;

  // Define the time structure for one layer
  // layer: which layer
  // nROFsTF: the total number rofs in a TF (the pattern is repeating per orbit so this is just a multiple, handling
  //          both edges correctly)
  // rofLength: ROF length in BC
  // rofDelay: delay of ROFs wrt LHC orbit
  // deltaROF: added delta in BC considering compatibility among rofs
  GPUh() void defineLayer(int32_t layer, int32_t nROFsTF, int32_t rofLength, int32_t rofDelay, int32_t deltaROF)
  {
    assert(std::numeric_limits<T>::max() < nROFsTF);
    mLayers[layer] = {nROFsTF, rofLength, rofDelay, deltaROF};
  }

  GPUh() void defineLayer(int32_t layer, const LayerTiming& timing)
  {
    mLayers[layer] = timing;
  }

  GPUhdi() const LayerTiming& getLayer(int32_t layer)
  {
    assert(layer >= 0 && layer < NLayers);
    return mLayers[layer];
  }

  GPUhdi() constexpr int32_t getEntries() noexcept { return NLayers; }
};

// GPU friendly view of the table below
template <int32_t NLayers, typename TableEntry, typename TableIndex>
struct ROFOverlapTableView {
  using ROFRange = LayerTiming::ROFRange;
  const TableEntry* mFlatTable{nullptr};
  const TableIndex* mIndices{nullptr};
  const LayerTiming* mLayers{nullptr};

  GPUhdi() const TableEntry& getOverlap(int32_t from, int32_t to, size_t rofIdx) const noexcept
  {
    assert(from < NLayers && to < NLayers);
    const size_t linearIdx = (from * NLayers) + to;
    const auto& idx = mIndices[linearIdx];
    assert(rofIdx < idx.getEntries());
    return mFlatTable[idx.getFirstEntry() + rofIdx];
  }

  GPUhdi() bool isCompatible(int32_t layer0, size_t rof0, int32_t layer1, size_t rof1) const noexcept
  {
    if (layer0 == layer1) { // layer is compatible with itself
      return rof0 == rof1;
    }

    assert(layer0 < NLayers && layer1 < NLayers);
    const size_t linearIdx = (layer0 * NLayers) + layer1;
    const auto& idx = mIndices[linearIdx];

    if (rof0 >= idx.getEntries()) {
      return false;
    }

    const auto& overlap = mFlatTable[idx.getFirstEntry() + rof0];

    if (overlap.getEntries() == 0) {
      return false;
    }

    const size_t firstCompatible = overlap.getFirstEntry();
    const size_t lastCompatible = firstCompatible + overlap.getEntries() - 1;
    return rof1 >= firstCompatible && rof1 <= lastCompatible;
  }

  GPUhdi() const LayerTiming& getLayer(int32_t layer) const noexcept
  {
    assert(layer >= 0 && layer < NLayers);
    return mLayers[layer];
  }

  GPUh() void printAll() const
  {
    for (int32_t i = 0; i < NLayers; ++i) {
      for (int32_t j = 0; j < NLayers; ++j) {
        if (i != j) {
          printMapping(i, j);
        }
      }
    }
    printSummary();
  }

  GPUh() void printMapping(int32_t from, int32_t to) const
  {
    if (from == to) {
      LOGP(error, "No self-lookup supported");
      return;
    }

    constexpr int w_index = 10;
    constexpr int w_first = 12;
    constexpr int w_last = 12;
    constexpr int w_count = 10;

    LOGF(info, "Overlap mapping: Layer %d -> Layer %d", from, to);
    LOGP(info, "From: {}", mLayers[from].asString());
    LOGP(info, "To  : {}", mLayers[to].asString());
    LOGF(info, "%*s | %*s | %*s | %*s", w_index, "ROF.index", w_first, "First.ROF", w_last, "Last.ROF", w_count, "Count");
    LOGF(info, "%.*s-+-%.*s-+-%.*s-+-%.*s", w_index, "----------", w_first, "------------", w_last, "------------", w_count, "----------");

    const size_t linearIdx = (from * NLayers) + to;
    const auto& idx = mIndices[linearIdx];
    for (int32_t i = 0; i < idx.getEntries(); ++i) {
      const auto& overlap = getOverlap(from, to, i);
      LOGF(info, "%*d | %*d | %*d | %*d", w_index, i, w_first, overlap.getFirstEntry(), w_last, overlap.getEntriesBound() - 1, w_count, overlap.getEntries());
    }
  }

  GPUh() void printSummary() const
  {
    uint32_t totalEntries{0};
    size_t flatTableSize{0};

    for (int32_t i = 0; i < NLayers; ++i) {
      for (int32_t j = 0; j < NLayers; ++j) {
        if (i != j) {
          const size_t linearIdx = i * NLayers + j;
          const auto& idx = mIndices[linearIdx];
          totalEntries += idx.getEntries();
          flatTableSize += idx.getEntries();
        }
      }
    }

    const uint32_t totalBytes = (flatTableSize * sizeof(TableEntry)) + (NLayers * NLayers * sizeof(TableIndex));
    LOGF(info, "------------------------------------------------------------");
    LOGF(info, "Total overlap table size: %u entries", totalEntries);
    LOGF(info, "Flat table size: %zu entries", flatTableSize);
    LOGF(info, "Total view size: %u bytes", totalBytes);
    LOGF(info, "------------------------------------------------------------");
  }
};

// Precalculated lookup table to find overlapping ROFs in another layer given a ROF index in the current layer
template <int32_t NLayers>
class ROFOverlapTable : public LayerTimingBase<NLayers>
{
 public:
  using T = LayerTimingBase<NLayers>::T;
  using TableEntry = dataformats::RangeReference<T, T>;
  using TableIndex = dataformats::RangeReference<T, T>;

  using View = ROFOverlapTableView<NLayers, TableEntry, TableIndex>;
  GPUdDefault() ROFOverlapTable() = default;

  GPUh() void init()
  {
    std::vector<TableEntry> table[NLayers][NLayers];
    for (int32_t i{0}; i < NLayers; ++i) {
      for (int32_t j{0}; j < NLayers; ++j) {
        if (i != j) { // we do not need self-lookup
          buildMapping(i, j, table[i][j]);
        }
      }
    }
    flatten(table);
  }

  GPUh() View getView() const
  {
    View view;
    view.mFlatTable = mFlatTable.data();
    view.mIndices = mIndices;
    view.mLayers = this->mLayers;
    return view;
  }

  GPUh() View getDeviceView(const TableEntry* deviceFlatTablePtr, const TableIndex* deviceIndicesPtr, const LayerTiming* deviceLayerTimingPtr) const
  {
    View view;
    view.mFlatTable = deviceFlatTablePtr;
    view.mIndices = deviceIndicesPtr;
    view.mLayers = deviceLayerTimingPtr;
    return view;
  }

  GPUh() size_t getFlatTableSize() const noexcept { return mFlatTable.size(); }
  static GPUh() constexpr size_t getIndicesSize() { return NLayers * NLayers; }

 private:
  GPUh() void buildMapping(int32_t from, int32_t to, std::vector<TableEntry>& table)
  {
    const auto& layerFrom = this->mLayers[from];
    const auto& layerTo = this->mLayers[to];
    table.resize(layerFrom.mNROFsTF);

    for (int32_t iROF{0}; iROF < layerFrom.mNROFsTF; ++iROF) {
      int32_t startFrom = layerFrom.mROFDelay + (iROF * layerFrom.mROFLength);
      int32_t endFrom = startFrom + layerFrom.mROFLength;
      startFrom -= layerFrom.mROFDelta;
      endFrom += layerFrom.mROFDelta;
      int32_t firstROFTo = o2::gpu::CAMath::Max(0, (startFrom - layerTo.mROFDelay) / layerTo.mROFLength);
      int32_t lastROFTo = (endFrom - layerTo.mROFDelay - 1) / layerTo.mROFLength;
      firstROFTo = o2::gpu::CAMath::Max(0, firstROFTo);
      lastROFTo = o2::gpu::CAMath::Min(layerTo.mNROFsTF - 1, lastROFTo);

      // verify overlap
      while (firstROFTo <= lastROFTo) {
        int32_t startTo = layerTo.mROFDelay + (firstROFTo * layerTo.mROFLength);
        int32_t endTo = startTo + layerTo.mROFLength;
        if (endTo > startFrom && startTo < endFrom) {
          break;
        }
        ++firstROFTo;
      }
      while (lastROFTo >= firstROFTo) {
        int32_t startTo = layerTo.mROFDelay + (lastROFTo * layerTo.mROFLength);
        int32_t endTo = startTo + layerTo.mROFLength;
        if (endTo > startFrom && startTo < endFrom) {
          break;
        }
        --lastROFTo;
      }
      int32_t count = (firstROFTo <= lastROFTo) ? (lastROFTo - firstROFTo + 1) : 0;
      table[iROF] = {static_cast<T>(firstROFTo), static_cast<T>(count)};
    }
  }

  GPUh() void flatten(const std::vector<TableEntry> table[NLayers][NLayers])
  {
    size_t total{0};
    for (int32_t i{0}; i < NLayers; ++i) {
      for (int32_t j{0}; j < NLayers; ++j) {
        if (i != j) { // we don not need self-lookup
          total += table[i][j].size();
        }
      }
    }

    mFlatTable.reserve(total);

    for (int32_t i{0}; i < NLayers; ++i) {
      for (int32_t j{0}; j < NLayers; ++j) {
        size_t idx = (i * NLayers) + j;
        if (i != j) {
          mIndices[idx].setFirstEntry(static_cast<T>(mFlatTable.size()));
          mIndices[idx].setEntries(static_cast<T>(table[i][j].size()));
          mFlatTable.insert(mFlatTable.end(), table[i][j].begin(), table[i][j].end());
        } else {
          mIndices[idx] = {0, 0};
        }
      }
    }
  }

  TableIndex mIndices[NLayers * NLayers];
  std::vector<TableEntry> mFlatTable;
};

// GPU friendly view of the table below
template <int32_t NLayers, typename TableEntry, typename TableIndex>
struct ROFVertexLookupTableView {
  const TableEntry* mFlatTable{nullptr};
  const TableIndex* mIndices{nullptr};
  const LayerTiming* mLayers{nullptr};

  GPUhdi() const LayerTiming& getLayer(int32_t layer) const noexcept
  {
    assert(layer >= 0 && layer < NLayers);
    return mLayers[layer];
  }

  GPUhdi() const TableEntry& getVertices(int32_t layer, size_t rofIdx) const noexcept
  {
    assert(layer < NLayers);
    const auto& idx = mIndices[layer];
    assert(rofIdx < idx.getEntries());
    return mFlatTable[idx.getFirstEntry() + rofIdx];
  }

  GPUh() int32_t getMaxVerticesPerROF() const noexcept
  {
    int32_t maxCount = 0;
    for (int32_t layer = 0; layer < NLayers; ++layer) {
      const auto& idx = mIndices[layer];
      for (int32_t i = 0; i < idx.getEntries(); ++i) {
        const auto& entry = mFlatTable[idx.getFirstEntry() + i];
        maxCount = o2::gpu::CAMath::Max(maxCount, static_cast<int32_t>(entry.getEntries()));
      }
    }
    return maxCount;
  }

  // Check if a specific vertex is compatible with a given ROF
  GPUhdi() bool isVertexCompatible(int32_t layer, size_t rofIdx, size_t vertexIdx) const noexcept
  {
    assert(layer < NLayers);
    const auto& idx = mIndices[layer];

    if (rofIdx >= idx.getEntries()) {
      return false;
    }

    const auto& entry = mFlatTable[idx.getFirstEntry() + rofIdx];

    if (entry.getEntries() == 0) {
      return false;
    }

    const size_t firstVertex = entry.getFirstEntry();
    const size_t lastVertex = firstVertex + entry.getEntries() - 1;
    return vertexIdx >= firstVertex && vertexIdx <= lastVertex;
  }

  GPUh() void printAll() const
  {
    for (int32_t i = 0; i < NLayers; ++i) {
      printLayer(i);
    }
    printSummary();
  }

  GPUh() void printLayer(int32_t layer) const
  {
    constexpr int w_rof = 10;
    constexpr int w_first = 12;
    constexpr int w_last = 12;
    constexpr int w_count = 10;

    LOGF(info, "Vertex lookup: Layer %d", layer);
    LOGF(info, "%*s | %*s | %*s | %*s", w_rof, "ROF.index", w_first, "First.Vtx", w_last, "Last.Vtx", w_count, "Count");
    LOGF(info, "%.*s-+-%.*s-+-%.*s-+-%.*s", w_rof, "----------", w_first, "------------", w_last, "------------", w_count, "----------");

    const auto& idx = mIndices[layer];
    for (int32_t i = 0; i < idx.getEntries(); ++i) {
      const auto& entry = mFlatTable[idx.getFirstEntry() + i];
      int first = entry.getFirstEntry();
      int count = entry.getEntries();
      int last = first + count - 1;
      LOGF(info, "%*d | %*d | %*d | %*d", w_rof, i, w_first, first, w_last, last, w_count, count);
    }
  }

  GPUh() void printSummary() const
  {
    uint32_t totalROFs{0};
    uint32_t totalVertexRefs{0};

    for (int32_t i = 0; i < NLayers; ++i) {
      const auto& idx = mIndices[i];
      totalROFs += idx.getEntries();

      for (int32_t j = 0; j < idx.getEntries(); ++j) {
        const auto& entry = mFlatTable[idx.getFirstEntry() + j];
        totalVertexRefs += entry.getEntries();
      }
    }

    const uint32_t totalBytes = (totalROFs * sizeof(TableEntry)) + (NLayers * sizeof(TableIndex));
    LOGF(info, "------------------------------------------------------------");
    LOGF(info, "Total ROFs in table: %u", totalROFs);
    LOGF(info, "Total vertex references: %u", totalVertexRefs);
    LOGF(info, "Total view size: %u bytes", totalBytes);
    LOGF(info, "------------------------------------------------------------");
  }
};

// Precalculated lookup table to find vertices compatible with ROFs
// Given a layer and ROF index, returns the range of vertices that overlap in time.
// The vertex time is defined as asymmetrical, it provides the beginning and range
// from the lowest common time bracket (BCs).
// e.g., [beginning, beginning+range)
template <int32_t NLayers>
class ROFVertexLookupTable : public LayerTimingBase<NLayers>
{
 public:
  using T = LayerTimingBase<NLayers>::T;
  using BCType = LayerTiming::BCType;
  using TableEntry = dataformats::RangeReference<T, T>;
  using TableIndex = dataformats::RangeReference<T, T>;

  using View = ROFVertexLookupTableView<NLayers, TableEntry, TableIndex>;

  GPUdDefault() ROFVertexLookupTable() = default;

  GPUh() size_t getFlatTableSize() const noexcept { return mFlatTable.size(); }
  static GPUh() constexpr size_t getIndicesSize() { return NLayers; }

  // Build the lookup table given a sorted array of vertices
  // vertices must be sorted by timestamp, then by error (secondary)
  GPUh() void init(const Vertex* vertices, size_t nVertices)
  {
    if (nVertices > std::numeric_limits<T>::max()) {
      LOGF(fatal, "too many vertices %zu, max supported is %u", nVertices, std::numeric_limits<T>::max());
    }

    std::vector<TableEntry> table[NLayers];
    for (int32_t layer{0}; layer < NLayers; ++layer) {
      buildMapping(layer, vertices, nVertices, table[layer]);
    }
    flatten(table);
  }

  // Pre-allocated needed memory, then use update(...)
  GPUh() void init()
  {
    size_t total{0};
    for (int32_t layer{0}; layer < NLayers; ++layer) {
      total += this->mLayers[layer].mNROFsTF;
    }
    mFlatTable.resize(total, {0, 0});
    size_t offset = 0;
    for (int32_t layer{0}; layer < NLayers; ++layer) {
      size_t nROFs = this->mLayers[layer].mNROFsTF;
      mIndices[layer].setFirstEntry(static_cast<T>(offset));
      mIndices[layer].setEntries(static_cast<T>(nROFs));
      offset += nROFs;
    }
    mNeedsUpdate = true;
  }

  GPUh() bool needsUpdate() const noexcept { return mNeedsUpdate; }

  // Recalculate lookup table with new vertices
  GPUh() void update(const Vertex* vertices, size_t nVertices)
  {
    size_t offset = 0;
    for (int32_t layer{0}; layer < NLayers; ++layer) {
      const auto& idx = mIndices[layer];
      size_t nROFs = idx.getEntries();
      for (size_t iROF = 0; iROF < nROFs; ++iROF) {
        updateROFMapping(layer, iROF, vertices, nVertices, offset + iROF);
      }
      offset += nROFs;
    }
    mNeedsUpdate = true;
  }

  GPUh() View getView() const
  {
    View view;
    view.mFlatTable = mFlatTable.data();
    view.mIndices = mIndices;
    view.mLayers = this->mLayers;
    return view;
  }

  GPUh() View getDeviceView(const TableEntry* deviceFlatTablePtr, const TableIndex* deviceIndicesPtr, const LayerTiming* deviceLayerTimingPtr) const
  {
    View view;
    view.mFlatTable = deviceFlatTablePtr;
    view.mIndices = deviceIndicesPtr;
    view.mLayers = deviceLayerTimingPtr;
    return view;
  }

 private:
  // Build the mapping for one layer
  GPUh() void buildMapping(int32_t layer, const Vertex* vertices, size_t nVertices, std::vector<TableEntry>& table)
  {
    const auto& layerDef = this->mLayers[layer];
    table.resize(layerDef.mNROFsTF);

    size_t vertexSearchStart = 0;

    for (int32_t iROF{0}; iROF < layerDef.mNROFsTF; ++iROF) {
      BCType rofStart = layerDef.mROFDelay + (iROF * layerDef.mROFLength);
      BCType rofEnd = rofStart + layerDef.mROFLength;
      rofStart -= layerDef.mROFDelta;
      rofEnd += layerDef.mROFDelta;
      size_t firstVertex = binarySearchFirst(vertices, nVertices, vertexSearchStart, rofStart);
      size_t lastVertex = firstVertex;
      while (lastVertex < nVertices && vertices[lastVertex].getTimeStamp().getTimeStamp() < rofEnd) {
        ++lastVertex;
      }
      size_t count = (lastVertex > firstVertex) ? (lastVertex - firstVertex) : 0;
      table[iROF] = {static_cast<T>(firstVertex), static_cast<T>(count)};
      vertexSearchStart = firstVertex;
    }
  }

  // Update a single ROF's vertex mapping
  GPUh() void updateROFMapping(int32_t layer, size_t iROF, const Vertex* vertices, size_t nVertices, size_t flatTableIdx)
  {
    const auto& layerDef = this->mLayers[layer];
    BCType rofStart = layerDef.mROFDelay + (iROF * layerDef.mROFLength);
    BCType rofEnd = rofStart + layerDef.mROFLength;
    rofStart -= layerDef.mROFDelta;
    rofEnd += layerDef.mROFDelta;
    size_t firstVertex = binarySearchFirst(vertices, nVertices, 0, rofStart);
    size_t lastVertex = firstVertex;
    while (lastVertex < nVertices && static_cast<int32_t>(vertices[lastVertex].getTimeStamp().getTimeStamp()) < rofEnd) {
      ++lastVertex;
    }
    size_t count = (lastVertex > firstVertex) ? (lastVertex - firstVertex) : 0;
    mFlatTable[flatTableIdx].setFirstEntry(static_cast<T>(firstVertex));
    mFlatTable[flatTableIdx].setEntries(static_cast<T>(count));
  }

  // Binary search for first vertex where maxBC >= targetBC
  GPUh() size_t binarySearchFirst(const Vertex* vertices, size_t nVertices, size_t searchStart, BCType targetBC) const
  {
    size_t left = searchStart;
    size_t right = nVertices;

    while (left < right) {
      size_t mid = left + ((right - left) / 2);
      if (static_cast<int32_t>(vertices[mid].getTimeStamp().getTimeStamp() + vertices[mid].getTimeStamp().getTimeStampError()) <= targetBC) {
        left = mid + 1;
      } else {
        right = mid;
      }
    }

    return left;
  }

  // Compress the temporary table into a single flat table
  GPUh() void flatten(const std::vector<TableEntry> table[NLayers])
  {
    // Count total entries
    size_t total{0};
    for (int32_t i{0}; i < NLayers; ++i) {
      total += table[i].size();
    }

    mFlatTable.reserve(total);

    // Build flat table and indices
    for (int32_t i{0}; i < NLayers; ++i) {
      mIndices[i].setFirstEntry(static_cast<T>(mFlatTable.size()));
      mIndices[i].setEntries(static_cast<T>(table[i].size()));
      mFlatTable.insert(mFlatTable.end(), table[i].begin(), table[i].end());
    }
  }

  bool mNeedsUpdate{false};
  TableIndex mIndices[NLayers];
  std::vector<TableEntry> mFlatTable;
};

// GPU friendly view of the table below
template <int32_t NLayers, typename BCRange, typename ROFRange>
struct ROFTimeSliceTableView {
  const BCRange* mFlatTable{nullptr};
  const uint8_t* mFlatMask{nullptr};
  const int32_t* mLayerROFOffsets{nullptr};
  const LayerTiming* mLayers{nullptr};
  int32_t mSlices{0};

  GPUhdi() const uint8_t* getMask(int32_t layer, int32_t slice) const
  {
    assert(layer < NLayers);
    assert(slice < mSlices);
    const auto& range = getSlice(layer, slice);
    return &mFlatMask[mLayerROFOffsets[layer] + range.getFirstEntry()];
  }

  // Get the ROF range for a given layer and slice
  GPUhdi() ROFRange getROFSlice(int32_t layer, int32_t slice) const
  {
    const auto& bcSlices = getSlice(layer, slice);
    const auto rofIdStart = mLayers[layer].getROFIdForBC(bcSlices.getFirstEntry());
    const auto rofIdEnd = mLayers[layer].getROFIdForBC(bcSlices.getEntriesBound());
    assert(rofIdEnd <= mLayers[layer].mNROFsTF);
    return {rofIdStart, rofIdEnd - rofIdStart};
  }

  // Get the BC range for a given layer and slice
  GPUhdi() const BCRange& getSlice(int32_t layer, int32_t slice) const
  {
    assert(layer < NLayers);
    assert(slice < mSlices);
    return mFlatTable[(layer * mSlices) + slice];
  }

  // Check if a BC falls within a specific slice
  GPUhdi() bool isInSlice(int32_t layer, int32_t slice, int32_t bc) const
  {
    const auto& range = getSlice(layer, slice);
    int32_t start = range.getFirstEntry();
    int32_t end = range.getEntriesBound();
    return bc >= start && bc < end;
  }

  // Find which slice a BC belongs to (returns -1 if not found)
  GPUhdi() int32_t findSlice(int32_t layer, int32_t bc) const
  {
    assert(layer < NLayers);
    int32_t left = 0;
    int32_t right = mSlices - 1;

    while (left <= right) {
      int32_t mid = left + ((right - left) / 2);
      const auto& range = getSlice(layer, mid);
      int32_t start = range.getFirstEntry();
      int32_t end = range.getEntriesBound();

      if (bc < start) {
        right = mid - 1;
      } else if (bc >= end) {
        left = mid + 1;
      } else {
        return mid;
      }
    }
    return -1;
  }

  GPUh() void printAll() const
  {
    for (int32_t i = 0; i < NLayers; ++i) {
      printLayer(i);
    }
  }

  GPUh() void printLayer(int32_t layer) const
  {
    constexpr int w_slice = 10;
    constexpr int w_start = 12;
    constexpr int w_end = 12;
    constexpr int w_range = 10;
    constexpr int w_rofs = 12;
    constexpr int w_active = 12;
    LOGF(info, "Slice table: Layer %d", layer);
    LOGF(info, "%*s | %*s | %*s | %*s | %*s | %*s", w_slice, "Slice", w_start, "BC.Start", w_end, "BC.End", w_range, "Range", w_rofs, "ROFs", w_active, "ActiveROFs");
    LOGF(info, "%.*s-+-%.*s-+-%.*s-+-%.*s-+-%.*s-+-%.*s", w_slice, "----------", w_start, "------------", w_end, "------------", w_range, "----------", w_rofs, "----------", w_active, "------------");
    for (int32_t i = 0; i < mSlices; ++i) {
      const auto& range = getSlice(layer, i);
      int32_t start = range.getFirstEntry();
      int32_t rangeLen = range.getEntries();
      int32_t end = range.getEntriesBound();
      const auto& rofRange = getROFSlice(layer, i);

      // Count active ROFs in this slice
      const uint8_t* mask = getMask(layer, i);
      int32_t activeCount = 0;
      for (int32_t j = 0; j < rangeLen; ++j) {
        activeCount += mask[j];
      }

      LOGF(info, "%*d | %*d | %*d | %*d | %*d-%d | %*d/%d", w_slice, i, w_start, start, w_end, end, w_range, rangeLen, w_rofs, rofRange.getFirstEntry(), rofRange.getEntriesBound(), w_active, activeCount, rangeLen);
    }
  }
};

template <int32_t NLayers>
class ROFTimeSliceTable : public LayerTimingBase<NLayers>
{
 public:
  using ROFRange = LayerTiming::ROFRange;
  using BCRange = dataformats::RangeReference<int32_t, int32_t>;
  using View = ROFTimeSliceTableView<NLayers, ROFRange, ROFRange>;

  GPUdDefault() ROFTimeSliceTable() = default;

  GPUh() void init(int32_t nSlices)
  {
    if (nSlices < 1) {
      LOGP(fatal, "Using {} time slices makes no sense", nSlices);
    }
    mSlices = nSlices;
    std::vector<ROFRange> table[NLayers];
    int32_t totalROFs = 0;
    for (int32_t layer{0}; layer < NLayers; ++layer) {
      buildSlices(layer, table[layer]);
      mLayerROFOffsets[layer] = totalROFs;
      totalROFs += this->getLayer(layer).mNROFsTF;
    }
    mFlatMask.resize(totalROFs, 1);
    flatten(table);
    mNeedsUpdate = true;
  }

  GPUh() size_t getFlatTableSize() const noexcept { return mFlatTable.size(); }
  GPUh() size_t getFlatMaskSize() const noexcept { return mFlatMask.size(); }

  GPUh() bool needsUpdate() const noexcept { return mNeedsUpdate; }

  GPUh() void selectROFs(const std::vector<BCRange>& ts)
  {
    resetMask();
    for (const auto& t : ts) {
      selectROFs(t.getFirstEntry(), t.getEntriesBound());
    }
  }

  GPUh() void selectROFs(int32_t bcStart, int32_t bcEnd)
  {
    for (int32_t layer{0}; layer < NLayers; ++layer) {
      const auto& lay = this->getLayer(layer);
      int32_t offset = mLayerROFOffsets[layer];
      for (int32_t rofId{0}; rofId < lay.mNROFsTF; ++rofId) {
        auto bounds = lay.getROFTimeBounds(rofId);
        int32_t rofStart = bounds.getFirstEntry();
        int32_t rofEnd = bounds.getEntriesBound();
        bool isCompatible = rofStart < bcEnd && rofEnd > bcStart;
        mFlatMask[offset + rofId] = isCompatible ? 1 : 0;
      }
    }
    mNeedsUpdate = true;
  }

  GPUh() void resetMask(int32_t s = 1)
  {
    std::memset(mFlatMask.data(), s, mFlatMask.size());
    mNeedsUpdate = true;
  }

  GPUh() void invertMask()
  {
    std::ranges::transform(mFlatMask, mFlatMask.begin(), [](uint8_t x) { return 1 - x; });
    mNeedsUpdate = true;
  }

  GPUh() View getView() const
  {
    View view;
    view.mFlatTable = mFlatTable.data();
    view.mFlatMask = mFlatMask.data();
    view.mLayerROFOffsets = mLayerROFOffsets;
    view.mLayers = this->mLayers;
    view.mSlices = mSlices;
    return view;
  }

  GPUh() View getDeviceView(const ROFRange* deviceFlatTablePtr, const uint8_t* deviceFlatMaskPtr, const int32_t* deviceOffsetPtr, const LayerTiming* deviceLayerTimingPtr) const
  {
    View view;
    view.mFlatTable = deviceFlatTablePtr;
    view.mFlatMask = deviceFlatMaskPtr;
    view.mLayerROFOffsets = deviceOffsetPtr;
    view.mLayers = deviceLayerTimingPtr;
    view.mSlices = mSlices;
    return view;
  }

 private:
  GPUh() void buildSlices(int32_t layer, std::vector<ROFRange>& table)
  {
    table.reserve(mSlices);
    const auto& lay = this->getLayer(layer);
    int32_t nBCsPerTF = (lay.mNROFsTF * lay.mROFLength);
    for (int32_t iSlice{0}; iSlice < mSlices; ++iSlice) {
      int32_t bcStart = (iSlice * nBCsPerTF) / mSlices;
      int32_t bcEnd = (((iSlice + 1) * nBCsPerTF) / mSlices) + 1;
      table.emplace_back(bcStart, bcEnd - bcStart);
    }
  }

  GPUh() void flatten(const std::vector<ROFRange> table[NLayers])
  {
    auto total = static_cast<size_t>(mSlices * NLayers);
    mFlatTable.reserve(total);
    for (int32_t layer = 0; layer < NLayers; ++layer) {
      mFlatTable.insert(mFlatTable.end(), table[layer].begin(), table[layer].end());
    }
  }

  bool mNeedsUpdate{false};
  int32_t mSlices{0};
  int32_t mLayerROFOffsets[NLayers]{};
  std::vector<ROFRange> mFlatTable;
  std::vector<uint8_t> mFlatMask;
};

} // namespace o2::its

#endif
