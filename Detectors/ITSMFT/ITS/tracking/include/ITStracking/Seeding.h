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

#ifndef TRACKINGITSU_INCLUDE_SEEDING_H_
#define TRACKINGITSU_INCLUDE_SEEDING_H_

#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <oneapi/tbb/parallel_for.h>

#include <cstdint>
#include <vector>
#include <array>

#include "Framework/Logger.h"
#include "ReconstructionDataFormats/Track.h"
#include "CommonDataFormat/RangeReference.h"
#include "ITStracking/Definitions.h"
#include "ITStracking/Constants.h"
#include "GPUCommonDef.h"
#include "GPUCommonMath.h"

namespace o2::its
{

using LTRef = dataformats::RangeReference<uint32_t, uint32_t>;

// Linearized track parametrisation near a vertex
// follows the definition of the pvertexer
struct LinearizedTrack {
  enum Status : uint8_t {
    kDead = (1 << 0),
  };
  GPUhdDefault() LinearizedTrack() = default;
  GPUhd() LinearizedTrack(const o2::track::TrackParCov& trk, int32_t cell) : cellIdx(cell), x(trk.getX()), y(trk.getY()), z(trk.getZ()), tgL(trk.getTgl()), tgP(trk.getSnp() / o2::gpu::CAMath::Sqrt(1.f - trk.getSnp()) * (1.f + trk.getSnp()))
  {
    o2::gpu::CAMath::SinCos(trk.getAlpha(), sinAlp, cosAlp);
    const float syy = trk.getSigmaY2(), szz = trk.getSigmaZ2(), syz = trk.getSigmaZY();
    const float det = (syy * szz) - (syz * syz);
    if (det <= constants::Tolerance) {
      markDead();
      return;
    }
    const float detI = 1.f / det;
    sig2YI = szz * detI;
    sig2ZI = syy * detI;
    sigYZI = -syz * detI;
  }

  float getChi2(const dataformats::VertexBase&) const;

  int32_t cellIdx{-1}; ///< index of attached cell
  float x{0.f};        ///< reference X
  float y{0.f};        ///< Y at X
  float z{0.f};        ///< Z at X
  float sig2YI{0.f};   ///< YY component of inverse cov.matrix
  float sig2ZI{0.f};   ///< ZZ component of inverse cov.matrix
  float sigYZI{0.f};   ///< YZ component of inverse cov.matrix
  float tgP{0.f};      ///< tangent(phi) in tracking frame
  float tgL{0.f};      ///< tangent(lambda)
  float cosAlp{0.f};   ///< cos of alpha frame
  float sinAlp{0.f};   ///< sin of alpha frame
  TimeEstBC time;      ///< timestamp; exclusive range of BCs that are spanned [time, time+error)
  uint8_t status{0};   ///< status bits

  GPUhdi() void markDead() noexcept { status |= kDead; }
  GPUhdi() bool isDead() const noexcept { return ((status & kDead) == kDead) || cellIdx < 0; }

  GPUh() void print() const
  {
    LOGP(info, "LT: cell={} x={} y={} z={} t={}+{} status={:b}", cellIdx, x, y, z, time.getTimeStamp(), time.getTimeStampError(), status);
  }

  ClassDefNV(LinearizedTrack, 1);
};

struct VertexSeed : public Vertex {
  enum FitStatus : uint8_t {
    kKilled,
    kConverged,
    kIterateFirst,
    kIterateFurther,
    kIterateRelaxScale,
  };
  uint8_t status = kIterateFirst;                                      // current fit status
  uint8_t iteration = 0;                                               // current iteration
  int32_t idx = -1;                                                    // db cluster index
  float wghSum = 0.f;                                                  // sum of tracks weights
  float wghChi2 = 0.f;                                                 // sum of tracks weighted chi2's
  float scaleSig2ITuk2I = 0.f;                                         // inverted scaled tukey weight
  float tukeyC = 0.f;                                                  // tukey constant
  float wghMin = 0.f;                                                  // minimum weight to account track
  ROOT::Math::SMatrix<float, 3, 3, ROOT::Math::MatRepSym<float, 3>> C; // C matrix
  ROOT::Math::SVector<float, 3> b;                                     // b vector

  void updateTukeyScale(const LinearizedTrack*, const int*, const LTRef&);
  void resetForNewIteration();
  bool acceptTrack(const LinearizedTrack&);
  void accountTrack(const LinearizedTrack&);
  void solveVertex();
};

namespace dbscan
{

// Do a scan in z and time
constexpr int32_t NDim{2};

// Configuration parameters
struct DBSCANParams {
  std::array<float, NDim> eps; // Maximum distance per dimension (z,t)
  int32_t minPts;              // Minimum points to form a dense region
};

// Clustering result
struct DBSCANResult {
  std::vector<int32_t> labels;
  std::vector<float> zCentroids;
  std::vector<LTRef> ranges;
  int32_t nClusters = 0;
  int32_t nNoise = 0;
};

// Neighbor list
using NeighborList = std::vector<std::vector<size_t>>;

// Point classification
constexpr int32_t DB_NOISE = -1;
constexpr int32_t DB_UNVISITED = -2;

// #define ITS_DB_MEASURE_TIMING
#ifdef ITS_DB_MEASURE_TIMING
class ScopedTimer
{
  std::string_view name;
  std::chrono::high_resolution_clock::time_point start;

 public:
  ScopedTimer(const ScopedTimer&) = default;
  ScopedTimer(ScopedTimer&&) = delete;
  ScopedTimer& operator=(const ScopedTimer&) = default;
  ScopedTimer& operator=(ScopedTimer&&) = delete;
  explicit ScopedTimer(std::string_view name)
    : name(name), start(std::chrono::high_resolution_clock::now()) {}

  ~ScopedTimer()
  {
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
    LOGP(info, "{} : {:.2f} ms", name, elapsed_ms);
  }
};
#define SCOPED_TIMER(name) ScopedTimer _timer##__LINE__(name)
#else
#define SCOPED_TIMER(name) ((void)0)
#endif

class DBSCANDistance
{
  using EPS = decltype(DBSCANParams::eps);

 public:
  DBSCANDistance() = default;
  DBSCANDistance(const EPS& eps) : mEps(eps) {}

  // Check if two linearized tracks are neighbors using L-infinity distance
  // Returns true if ALL dimensions are within their respective thresholds
  // Important to remember is that time is given as an asymmetric bracket
  bool areNeighbours(const LinearizedTrack& l0, const LinearizedTrack& l1) const
  {
    if (l0.isDead() || l1.isDead()) {
      return false; // dead tracks can never be neighbors of anything
    }
    const float diffZ = o2::gpu::CAMath::Abs(l0.z - l1.z);
    if (diffZ > mEps[0]) {
      return false;
    }
    return isTimeCompatible(l0, l1);
  }

  float getDistance(const LinearizedTrack& l0, const LinearizedTrack& l1) const
  {
    if (l0.isDead() || l1.isDead()) {
      return std::numeric_limits<float>::max();
    }
    if (isTimeCompatible(l0, l1)) {
      return std::numeric_limits<float>::max();
    }
    return o2::gpu::CAMath::Abs(l0.z - l1.z);
  }

 private:
  bool isTimeCompatible(const LinearizedTrack& l0, const LinearizedTrack& l1) const
  {
    const auto end0 = l0.time.getTimeStamp() + l0.time.getTimeStampError();
    const auto end1 = l1.time.getTimeStamp() + l1.time.getTimeStampError();
    const int diffT = o2::gpu::CAMath::Max(0, o2::gpu::CAMath::Max((int)l0.time.getTimeStamp(), (int)l1.time.getTimeStamp())) - (int)o2::gpu::CAMath::Min(end0, end1);
    return diffT <= (int)mEps[1];
  }

  EPS mEps;
};

class DBSCAN
{
 public:
  DBSCAN() = default;
  DBSCAN(const DBSCANParams& p);

  DBSCANResult cluster(const LinearizedTrack* tracks, size_t n) const;

  const auto& getParams() const { return mParams; }
  const auto& getDistance() const { return mDistance; }

 private:
  void findNeighbors(const LinearizedTrack* tracks, size_t n, NeighborList& neighbors) const;
  void classify(size_t n, const NeighborList& neighbors, std::vector<int32_t>& labels) const;
  void finalize(const LinearizedTrack* tracks, size_t n, DBSCANResult& result) const;

  DBSCANParams mParams;
  DBSCANDistance mDistance;
};

// Grid cell for spatial partitioning
using GridCell = std::vector<size_t>;
//  Grid coordinates
using GridCoord = std::array<int32_t, NDim>;
/// General index grid
class Grid
{
 public:
  Grid(const LinearizedTrack* tracks, size_t n, const std::array<float, NDim>& cellSizes)
    : mTracks(tracks), mNPoints(n), mCellSizes(cellSizes)
  {
    SCOPED_TIMER("Grid construction");
    computeBounds();
    computeGridDimensions();
    allocateCells();
    assignCells();
  }

  // Get grid coordinates for a point
  [[nodiscard]] GridCoord getGridCoords(size_t idx) const
  {
    GridCoord coords{};
    const auto& trk = mTracks[idx];
    // z
    coords[0] = static_cast<int32_t>((trk.z - mMinBounds[0]) / mCellSizes[0]);
    coords[0] = std::clamp(coords[0], 0, static_cast<int32_t>(mGridDims[0]) - 1);
    // time
    coords[1] = static_cast<int32_t>(((trk.time.getTimeStamp() + (trk.time.getTimeStampError() / 2)) - mMinBounds[1]) / mCellSizes[1]);
    coords[1] = std::clamp(coords[1], 0, static_cast<int32_t>(mGridDims[1]) - 1);
    return coords;
  }

  // Get flat index from grid coordinates
  [[nodiscard]] size_t getCellIndex(const GridCoord& coords) const
  {
    int32_t index = 0, stride = 1;
#pragma unroll(NDim)
    for (size_t d = 0; d < NDim; ++d) {
      index += coords[d] * stride;
      stride *= mGridDims[d];
    }
    return static_cast<size_t>(index);
  }

  // Get cell at grid coordinates
  [[nodiscard]] const GridCell* getCell(const GridCoord& coords) const
  {
#pragma unroll(NDim)
    for (size_t d = 0; d < NDim; ++d) {
      if (coords[d] < 0 || coords[d] >= static_cast<int32_t>(mGridDims[d])) {
        return nullptr;
      }
    }
    return &mCells[getCellIndex(coords)];
  }

  // Get neighboring cells (including the cell itself)
  void getNeighborCells(const GridCoord& coords, std::vector<const GridCell*>& neighbors) const
  {
    neighbors.clear();
    neighbors.reserve(static_cast<size_t>(std::pow(3, NDim)));
    GridCoord offset{};
    enumerateNeighborOffsets<0>(coords, offset, neighbors);
  }

 private:
  template <int32_t Dim>
  void enumerateNeighborOffsets(const GridCoord& base, GridCoord& offset, std::vector<const GridCell*>& output) const
  {
    if constexpr (Dim == NDim) {
      GridCoord nbr;
#pragma unroll(NDim)
      for (size_t d = 0; d < NDim; ++d) {
        nbr[d] = base[d] + offset[d];
      }
      const GridCell* cell = getCell(nbr);
      if (cell) {
        output.push_back(cell);
      }
      return;
    } else {
      for (int32_t v = -1; v <= 1; ++v) {
        offset[Dim] = v;
        enumerateNeighborOffsets<Dim + 1>(base, offset, output);
      }
    }
  }

  void computeBounds()
  {
    mMinBounds.fill(std::numeric_limits<float>::max());
    mMaxBounds.fill(std::numeric_limits<float>::lowest());
    for (size_t i{0}; i < mNPoints; ++i) {
      const auto& trk = mTracks[i];
      // z
      mMinBounds[0] = std::min(mMinBounds[0], trk.z);
      mMaxBounds[0] = std::max(mMaxBounds[0], trk.z);
      // time
      mMinBounds[1] = std::min(mMinBounds[1], (float)trk.time.getTimeStamp());
      mMaxBounds[1] = std::max(mMaxBounds[1], (float)(trk.time.getTimeStamp() + (trk.time.getTimeStampError() / 2)));
    }
  }

  void computeGridDimensions()
  {
    mGridDims.fill(1);
#pragma unroll(NDim)
    for (size_t d = 0; d < NDim; ++d) {
      float range = mMaxBounds[d] - mMinBounds[d];
      mGridDims[d] = std::max(size_t(1), static_cast<size_t>(std::ceil(range / mCellSizes[d])));
    }
  }

  void allocateCells()
  {
    auto total_cells = static_cast<size_t>(mGridDims[0] * mGridDims[1]);
    mCells.resize(total_cells);
  }

  void assignCells()
  {
    for (size_t i = 0; i < mNPoints; ++i) {
      mCells[getCellIndex(getGridCoords(i))].push_back(i);
    }
    tbb::parallel_for(size_t(0), mCells.size(), [&](size_t c) {
      std::sort(mCells[c].begin(), mCells[c].end(), [&](size_t a, size_t b) {
        return mTracks[a].z < mTracks[b].z;
      });
    });
  }

  const LinearizedTrack* mTracks;
  size_t mNPoints;
  std::array<float, NDim> mCellSizes;
  std::array<float, NDim> mMinBounds{};
  std::array<float, NDim> mMaxBounds{};
  std::array<size_t, NDim> mGridDims{};
  std::vector<GridCell> mCells;
};

} // namespace dbscan

} // namespace o2::its

#endif
