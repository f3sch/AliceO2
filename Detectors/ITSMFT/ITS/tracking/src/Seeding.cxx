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

#include <oneapi/tbb/parallel_for.h>

#include <cstdint>
#include <limits>
#include <vector>

#include "GPUCommonMath.h"
#include "ITStracking/Seeding.h"

namespace o2::its
{

float LinearizedTrack::getChi2(const dataformats::VertexBase& v) const
{
  // get residuals (Y and Z DCA in track frame) and calculate chi2
  float dx = (v.getX() * cosAlp) + (v.getY() * sinAlp) - x; // VX rotated to track frame - trackX
  float dy = y + (tgP * dx) - (-v.getX() * sinAlp + v.getY() * cosAlp);
  float dz = z + (tgL * dx) - v.getZ();
  float chi2 = (dy * dy * sig2YI + dz * dz * sig2ZI) + (2.f * dy * dz * sigYZI);
  return chi2 / 2.f; // not using time
}

void VertexSeed::updateTukeyScale(const LinearizedTrack* ltracks, const int* labels, const LTRef& ref)
{
  constexpr float sigma2Gaus{1.4826f};
  if (status == kIterateFirst) {
    scaleSig2ITuk2I = 1.f / (tukeyC * tukeyC);
    return;
  }
  if (status == kIterateRelaxScale) {
    // relax by 1/3
    scaleSig2ITuk2I *= 2.f / 3.f;
    return;
  }

  std::vector<float> chi2s;
  chi2s.reserve(getNContributors());

  for (uint32_t entry{ref.getFirstEntry()}; entry < ref.getEntriesBound(); ++entry) {
    const auto& lt = ltracks[entry];
    if (lt.isDead() || labels[entry] != idx) {
      continue;
    }
    chi2s.push_back(lt.getChi2((const dataformats::VertexBase&)*this));
  }
  if (chi2s.empty()) {
    return;
  }
  // get MAD scale
  std::sort(chi2s.begin(), chi2s.end());
  float median = chi2s[chi2s.size() / 2];
  std::vector<float> dev;
  dev.reserve(chi2s.size());
  for (float c : chi2s) {
    dev.push_back(o2::gpu::CAMath::Abs(c - median));
  }
  std::sort(dev.begin(), dev.end());
  float mad = dev[dev.size() / 2];
  float scale = sigma2Gaus * mad;
  scaleSig2ITuk2I = 1.f / (scale * tukeyC * tukeyC);
}

void VertexSeed::resetForNewIteration()
{
  wghSum = wghChi2 = 0.f;
  std::fill(C.begin(), C.end(), 0.f);
  std::fill(b.begin(), b.end(), 0.f);
}

bool VertexSeed::acceptTrack(const LinearizedTrack& lt)
{
  float chi2Red = lt.getChi2((const dataformats::VertexBase&)*this);
  float wghT = (1.f - (chi2Red * scaleSig2ITuk2I)); // weighted distance to vertex
  return wghT >= wghMin;
}

void VertexSeed::accountTrack(const LinearizedTrack& lt)
{
  float chi2Red = lt.getChi2((const dataformats::VertexBase&)*this);
  float wghT = (1.f - (chi2Red * scaleSig2ITuk2I)); // weighted distance to vertex
  if (wghT < wghMin) {
    return;
  }
  wghT *= wghT;
  wghSum += wghT;
  wghChi2 += wghT * chi2Red;
  float syyI(lt.sig2YI), szzI(lt.sig2ZI), syzI(lt.sigYZI); // reweighted inverse cov
  syyI *= wghT;
  syzI *= wghT;
  szzI *= wghT;
  // solve line equation
  float tmpSP = lt.sinAlp * lt.tgP, tmpCP = lt.cosAlp * lt.tgP,
        tmpSC = lt.sinAlp + tmpCP, tmpCS = -lt.cosAlp + tmpSP,
        tmpCL = lt.cosAlp * lt.tgL, tmpSL = lt.sinAlp * lt.tgL,
        tmpYXP = lt.y - (lt.tgP * lt.x), tmpZXL = lt.z - (lt.tgL * lt.x),
        tmpCLzz = tmpCL * szzI, tmpSLzz = tmpSL * szzI, tmpSCyz = tmpSC * syzI,
        tmpCSyz = tmpCS * syzI, tmpCSyy = tmpCS * syyI, tmpSCyy = tmpSC * syyI,
        tmpSLyz = tmpSL * syzI, tmpCLyz = tmpCL * syzI;
  // symmetric matrix equation
  C(0, 0) += tmpCL * (tmpCLzz + tmpSCyz + tmpSCyz) + tmpSC * tmpSCyy;         // dchi^2/dx/dx
  C(0, 1) += tmpCL * (tmpSLzz + tmpCSyz) + tmpSL * tmpSCyz + tmpSC * tmpCSyy; // dchi^2/dx/dy
  C(0, 2) += -lt.sinAlp * syzI - tmpCLzz - tmpCP * syzI;                      // dchi^2/dx/dz
  C(1, 1) += tmpSL * (tmpSLzz + tmpCSyz + tmpCSyz) + tmpCS * tmpCSyy;         // dchi^2/dy/dy
  C(1, 2) += -(tmpCSyz + tmpSLzz);                                            // dchi^2/dy/dz
  C(2, 2) += szzI;                                                            // dchi^2/dz/dz
  // RHS
  b(0) += -(tmpCLyz + tmpSCyy) * tmpYXP - (tmpCLzz + tmpSCyz) * tmpZXL;
  b(1) += -tmpYXP * (tmpCSyy + tmpSLyz) - tmpZXL * (tmpCSyz + tmpSLzz);
  b(2) += tmpZXL * szzI + tmpYXP * syzI;
  // account as contributor
  addContributor();
}

void VertexSeed::solveVertex()
{
  // solve C*a=b by inversion a=C^-1*b
  if (!C.InvertFast()) {
    status = kIterateFurther;
    return;
  }
  auto sol = C * b;
  setXYZ(sol(0), sol(1), sol(2));
  setCov(C(0, 0), C(1, 0), C(1, 1), C(2, 0), C(2, 1), C(2, 2));
  setChi2(wghChi2 / o2::gpu::CAMath::Max(1.f, wghSum - 3.f));
}

namespace dbscan
{

DBSCAN::DBSCAN(const DBSCANParams& p) : mParams(p), mDistance(mParams.eps) {}

DBSCANResult DBSCAN::cluster(const LinearizedTrack* tracks, size_t n) const
{
  DBSCANResult result;
  result.labels.resize(n, DB_UNVISITED);
  if (n == 0) {
    return result;
  }

  // Step 1: Find neighbors for all points using grid
  NeighborList neighbors;
  findNeighbors(tracks, n, neighbors);
  // Step 2: Classify points and form clusters
  classify(n, neighbors, result.labels);
  // Step 3: finalize results
  finalize(tracks, n, result);

  return result;
}

void DBSCAN::findNeighbors(const LinearizedTrack* tracks, size_t n, NeighborList& neighbors) const
{
  Grid grid(tracks, n, mParams.eps);

  // Parallel neighbor finding
  SCOPED_TIMER("find neighbors");
  neighbors.resize(n);
  tbb::parallel_for(
    tbb::blocked_range<size_t>(0, n), [&](const tbb::blocked_range<size_t>& range) {
      std::vector<const GridCell*> neighborCells;
      neighborCells.reserve(9);

      for (size_t i = range.begin(); i < range.end(); ++i) {
        const auto& query = tracks[i];
        if (query.isDead()) {
          continue;
        }
        auto coords = grid.getGridCoords(i);
        grid.getNeighborCells(coords, neighborCells);
        for (const GridCell* cell : neighborCells) {
          for (auto idx : *cell) {
            if (idx == i) {
              continue;
            }
            const float diffZ = tracks[idx].z - query.z;
            if (diffZ > mParams.eps[0]) { // sorted in z the rest is not reachable
              break;
            }
            if (diffZ < -mParams.eps[0]) {
              continue;
            }
            if (mDistance.areNeighbours(query, tracks[idx])) {
              neighbors[i].push_back(idx);
            }
          }
        }
      }
    });
}

void DBSCAN::classify(size_t n, const NeighborList& neighbors, std::vector<int32_t>& labels) const
{
  SCOPED_TIMER("classify");

  std::vector<std::atomic<size_t>> parent(n);
  std::vector<bool> isCore(n, false);

  tbb::parallel_for(size_t(0), n, [&](size_t i) {
    parent[i].store(i, std::memory_order_relaxed);
    isCore[i] = static_cast<int32_t>(neighbors[i].size()) >= mParams.minPts;
  });

  auto find = [&](size_t i) -> size_t {
    size_t root = i;
    while (parent[root].load(std::memory_order_relaxed) != root) {
      root = parent[root].load(std::memory_order_relaxed);
    }
    while (parent[i].load(std::memory_order_relaxed) != root) {
      size_t next = parent[i].load(std::memory_order_relaxed);
      parent[i].compare_exchange_weak(next, root, std::memory_order_relaxed);
      i = next;
    }
    return root;
  };

  auto unite = [&](size_t a, size_t b) {
    while (true) {
      size_t ra = find(a), rb = find(b);
      if (ra == rb) {
        return;
      }
      // Always make smaller index the root (deterministic)
      if (ra > rb) {
        std::swap(ra, rb);
      }
      size_t expected = rb;
      if (parent[rb].compare_exchange_weak(expected, ra, std::memory_order_relaxed)) {
        return;
      }
    }
  };

  // Union core-core neighbors
  tbb::parallel_for(size_t(0), n, [&](size_t i) {
    if (!isCore[i]) {
      return;
    }
    for (size_t j : neighbors[i]) {
      if (isCore[j]) {
        unite(i, j);
      }
    }
  });

  // Collect unique roots and sort for deterministic cluster IDs
  std::vector<size_t> roots;
  roots.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    if (isCore[i] && find(i) == i) {
      roots.push_back(i);
    }
  }
  std::sort(roots.begin(), roots.end());

  // Map root -> cluster ID (smallest root = cluster 0, etc.)
  std::vector<int32_t> rootToCluster(n, -1);
  for (int32_t c = 0; c < static_cast<int32_t>(roots.size()); ++c) {
    rootToCluster[roots[c]] = c;
  }

  // Label points
  tbb::parallel_for(size_t(0), n, [&](size_t i) {
    if (isCore[i]) {
      labels[i] = rootToCluster[find(i)];
    } else {
      // Border: assign to cluster of smallest-index core neighbor
      size_t minCoreNeighbor = std::numeric_limits<size_t>::max();
      for (size_t j : neighbors[i]) {
        if (isCore[j] && j < minCoreNeighbor) {
          minCoreNeighbor = j;
        }
      }
      if (minCoreNeighbor != std::numeric_limits<size_t>::max()) {
        labels[i] = rootToCluster[find(minCoreNeighbor)];
      } else {
        labels[i] = DB_NOISE;
      }
    }
  });
}

void DBSCAN::finalize(const LinearizedTrack* tracks, size_t n, DBSCANResult& result) const
{
  SCOPED_TIMER("finalizing");
  // Phase 1: Count clusters and noise points
  int32_t max_label = *std::ranges::max_element(result.labels);
  result.nClusters = max_label + 1;
  result.nNoise = static_cast<int32_t>(std::count(result.labels.begin(), result.labels.end(), DB_NOISE));
  // Phase 2: compute z-centroids and also set the ranges
  if (result.nClusters > 0) {
    result.zCentroids.resize(result.nClusters);
    result.ranges.resize(result.nClusters);
    tbb::parallel_for(0, result.nClusters, [&](const int32_t iCls) {
      float count{0.f}, zSum = {0.f};
      uint32_t first{std::numeric_limits<uint32_t>::max()}, last{0};
      for (size_t iTrk{0}; iTrk < n; ++iTrk) { // FIXME we rescan everything here, this can be done smarter
        if (tracks[iTrk].isDead()) {
          continue;
        }
        if (result.labels[iTrk] == iCls) {
          first = o2::gpu::CAMath::Min(first, (uint32_t)iTrk);
          last = o2::gpu::CAMath::Max(last, (uint32_t)iTrk);
          zSum += tracks[iTrk].z;
          ++count;
        }
      }
      result.ranges[iCls].set(first, last - first + 1);
      result.zCentroids[iCls] = zSum / count;
    });
  }
}

} // namespace dbscan
} // namespace o2::its
