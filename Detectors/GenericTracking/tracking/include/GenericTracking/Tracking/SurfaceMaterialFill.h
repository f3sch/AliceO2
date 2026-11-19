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

/// \file SurfaceMaterialFill.h
/// \brief The surface map's tables (SurfaceMaterial.h) filled from the
/// material LUT and the sensor planes, on the host once per run.

#ifndef GTRK_TRACKING_SURFACEMATERIALFILL_H_
#define GTRK_TRACKING_SURFACEMATERIALFILL_H_

#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Parallel.h"
#include "GenericTracking/Tracking/SearchWindow.h"
#include "GenericTracking/Tracking/SensorCells.h"
#include "GenericTracking/Tracking/SurfaceMaterial.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include "DetectorsBase/MatLayerCylSet.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <span>
#include <vector>

namespace o2::gtrk
{

/// Where the ray from the nominal vertex along \a u crosses one of layer
/// \a id's sensor planes, from a grid of the sensors over the layer's
/// azimuth and second coordinate: the farthest crossing with \a last, else
/// the nearest; the ray's parameter t (the point is t u), or <= 0 if it
/// crosses no sensor there.
template <DetectorLayout Det>
class SensorPlaneGrid
{
 public:
  static constexpr int kPhi = SurfaceMaterialMap<Det>::kPhiBins;
  static constexpr int kSecond = SurfaceMaterialMap<Det>::kSecondBins;

  SensorPlaneGrid(std::span<const TrackingSurface> surfaces, LayerID id) : mBarrel(Detector<Det>::isBarrelLayer(id)), mSecondMin(secondMin<Det>(id)), mPerCm(kSecond / (secondMax<Det>(id) - mSecondMin))
  {
    mCells.resize(static_cast<std::size_t>(kPhi) * kSecond);
    for (const auto& surface : surfaces) {
      if (surface.layer != id) {
        continue;
      }
      const auto footprint = SensorCells<Det>::footprint(surface, mBarrel);
      const double phiLo = footprint.phiCentre + footprint.phiLow;
      const double phiHi = footprint.phiCentre + footprint.phiHigh;
      const double sLo = footprint.secondLow;
      const double sHi = footprint.secondHigh;
      const int s0 = std::max(0, static_cast<int>(std::floor((sLo - mSecondMin) * mPerCm)));
      const int s1 = std::min(kSecond - 1, static_cast<int>(std::floor((sHi - mSecondMin) * mPerCm)));
      const int p0 = static_cast<int>(std::floor(phiLo * kPhi / kTwoPi<double>));
      const int p1 = static_cast<int>(std::floor(phiHi * kPhi / kTwoPi<double>));
      for (int ip = p0; ip <= p1; ++ip) {
        const int wrapped = ((ip % kPhi) + kPhi) % kPhi;
        for (int is = s0; is <= s1; ++is) {
          mCells[(static_cast<std::size_t>(wrapped) * kSecond) + is].push_back(&surface);
        }
      }
    }
  }

  /// The ray's first and last sensor crossings, t <= 0 for none: the ray
  /// along \a u from the nominal vertex, which crosses the layer's nominal
  /// surface in the grid's azimuth bin \a phiBin at \a second.
  void crossings(const double u[3], int phiBin, double second, double& first, double& last) const
  {
    first = -1.;
    last = -1.;
    const int ip = std::clamp(phiBin, 0, kPhi - 1);
    const int is = std::clamp(static_cast<int>((second - mSecondMin) * mPerCm), 0, kSecond - 1);
    const auto take = [&](double t) {
      first = first <= 0. || t < first ? t : first;
      last = last <= 0. || t > last ? t : last;
    };
    for (const TrackingSurface* surface : mCells[(static_cast<std::size_t>(ip) * kSecond) + is]) {
      if (surface->radius > 0.F) {
        // |w - (w.n) n| = radius for w = t u - origin, n the axis.
        const double o[3] = {surface->origin[0], surface->origin[1], surface->origin[2]};
        const double un = dotOf(u, surface->normal);
        const double on = dotOf(o, surface->normal);
        double up[3];
        double op[3];
        for (int k = 0; k < 3; ++k) {
          up[k] = u[k] - (un * surface->normal[k]);
          op[k] = o[k] - (on * surface->normal[k]);
        }
        const double qa = dotOf(up, up);
        const double qb = -2. * dotOf(up, op);
        const double reach = double(surface->radius) + surface->depth;
        const double qc = dotOf(op, op) - (reach * reach);
        const double disc = (qb * qb) - (4. * qa * qc);
        if (qa <= 0. || disc < 0.) {
          continue;
        }
        const double t = (-qb + std::sqrt(disc)) / (2. * qa);
        if (!(t > 0.)) {
          continue;
        }
        const double point[3] = {t * u[0], t * u[1], t * u[2]};
        double local[3];
        surfaceToLocal(*surface, point, local);
        if (!surface->bounds.isInside(static_cast<Scalar>(surface->radius * std::atan2(local[1], local[0])), static_cast<Scalar>(local[2]))) {
          continue;
        }
        take(t);
        continue;
      }
      constexpr double kVertex[3] = {0., 0., 0.};
      double t = 0.;
      double loc[2];
      if (std::abs(cosIncidence(*surface, u)) < 1e-9 || !lineCrossing(*surface, kVertex, u, t, loc) || !(t > 0.)) {
        continue;
      }
      if (!surface->bounds.isInside(static_cast<Scalar>(loc[0]), static_cast<Scalar>(loc[1]))) {
        continue;
      }
      take(t);
    }
  }

 private:
  bool mBarrel;
  double mSecondMin = 0.;
  double mPerCm = 0.;
  std::vector<std::vector<const TrackingSurface*>> mCells;
};

/// The linked pairs' material (SurfaceMaterial.h) from the LUT: for each cell
/// of a layer, rays from the nominal vertex through points spread over the
/// cell (thin walls fall between cell centres), each walked from where it
/// crosses a sensor plane of the cell's layer to the midpoint between that
/// and its crossing of a linked layer -- the sensor planes, not the layers'
/// nominal radii or z, since a stave's own material splits between the steps
/// in and out at its sensor. Where a ray crosses two sensors of a layer
/// (staggered or overlapping staves: a step may leave or reach either), the
/// mean of the walks from both and the midpoint between the mean crossings;
/// where it crosses none, the nominal surface. A step between a pair starts
/// and ends on a sensor and had no row on a layer in between: a cell takes
/// the rays that cross a sensor at both ends and none between, else those
/// that cross none between, else all.
///
/// A ray serves every pair of its layer: its crossings of all layers are
/// found once, and it is walked once outward and once inward, cut at the
/// pairs' midpoints in order. The layers' azimuth bins in parallel on
/// \a arena's threads, the tracker's.
template <DetectorLayout Det>
inline void fillSurfaceMaterial(SurfaceMaterialMap<Det>& map, const o2::base::MatLayerCylSet& lut, std::span<const TrackingSurface> surfaces, detail::TaskArena& arena)
{
  using Map = SurfaceMaterialMap<Det>;
  constexpr int kL = Map::kLayers;
  map.pairOf.fill(0);
  const auto surfaceOf = [](LayerID id) {
    const double sMin = secondMin<Det>(id);
    const double sMax = secondMax<Det>(id);
    return MaterialSurface{.secondMin = static_cast<Scalar>(sMin), .secondBinsPerCm = static_cast<Scalar>(Map::kSecondBins / (sMax - sMin)), .barrel = static_cast<std::uint8_t>(Detector<Det>::isBarrelLayer(id) ? 1U : 0U)};
  };
  // Each layer's halves: the pair, which half, the pair's other layer.
  struct Half {
    int pair;
    int half;
    LayerID other;
  };
  std::array<std::vector<Half>, kL> halves;
  for (int pair = 0; pair < Map::kPairs; ++pair) {
    const auto& link = Det::layerLinks[pair];
    map.surfaces[pair] = {surfaceOf(link.from), surfaceOf(link.to)};
    map.pairOf[(static_cast<std::size_t>(link.from) * kL) + link.to] = static_cast<std::int16_t>(pair + 1);
    map.pairOf[(static_cast<std::size_t>(link.to) * kL) + link.from] = static_cast<std::int16_t>(-(pair + 1));
    halves[link.from].push_back({.pair = pair, .half = 0, .other = link.to});
    halves[link.to].push_back({.pair = pair, .half = 1, .other = link.from});
  }
  std::vector<std::unique_ptr<SensorPlaneGrid<Det>>> grids;
  grids.reserve(kL);
  for (int layer = 0; layer < kL; ++layer) {
    grids.push_back(std::make_unique<SensorPlaneGrid<Det>>(surfaces, static_cast<LayerID>(layer)));
  }
  // Where the vertex ray along u crosses layer id's nominal surface: t, or <= 0.
  const auto nominal = [](LayerID id, const double u[3]) {
    if (Detector<Det>::isBarrelLayer(id)) {
      const double r = std::hypot(u[0], u[1]);
      return r > 0. ? layerDistance<Det>(id) / r : -1.;
    }
    return u[2] != 0. ? discZ<Det>(id) / u[2] : -1.;
  };
  constexpr int kSubPhi = 4;
  constexpr int kSubSecond = 2;
  detail::parallelFor(arena, 0U, static_cast<std::size_t>(kL) * Map::kPhiBins, 1U, [&](std::size_t begin, std::size_t end) {
    // A layer's crossing of the ray: the first and last sensor (or the
    // nominal surface twice), whether on a sensor.
    struct Crossing {
      double first = -1.;
      double last = -1.;
      double nominalT = -1.;
      bool onSensor = false;
    };
    struct Sum {
      double x2X0 = 0.;
      double xRho = 0.;
      int n = 0;
    };
    // A half's walk end and its result for the ray.
    struct Cut {
      double t;
      std::size_t index;
    };
    std::array<Crossing, kL> crossing;
    std::vector<std::array<Sum, 3>> sums; // sensed, clear, all
    std::vector<double> mid;
    std::vector<double> walkedX2X0;
    std::vector<double> walkedXRho;
    std::vector<Cut> outward;
    std::vector<Cut> inward;
    for (std::size_t job = begin; job < end; ++job) {
      const auto own = static_cast<LayerID>(job / Map::kPhiBins);
      const int iPhi = static_cast<int>(job % Map::kPhiBins);
      const auto& mine = halves[own];
      if (mine.empty()) {
        continue;
      }
      const bool barrel = Detector<Det>::isBarrelLayer(own);
      const MaterialSurface surface = surfaceOf(own);
      mid.resize(mine.size());
      walkedX2X0.resize(mine.size());
      walkedXRho.resize(mine.size());
      for (int iSecond = 0; iSecond < Map::kSecondBins; ++iSecond) {
        sums.assign(mine.size(), {});
        for (int a = 0; a < kSubPhi; ++a) {
          const double phi = kTwoPi<double> * (iPhi + ((a + 0.5) / kSubPhi)) / Map::kPhiBins;
          const double c = std::cos(phi);
          const double s = std::sin(phi);
          for (int b = 0; b < kSubSecond; ++b) {
            const double second = surface.secondMin + ((iSecond + ((b + 0.5) / kSubSecond)) / surface.secondBinsPerCm);
            double p[3];
            if (barrel) {
              const double r = layerDistance<Det>(own);
              p[0] = r * c;
              p[1] = r * s;
              p[2] = second;
            } else {
              p[0] = second * c;
              p[1] = second * s;
              p[2] = discZ<Det>(own);
            }
            const double norm = std::sqrt(dotOf(p, p));
            const double u[3] = {p[0] / norm, p[1] / norm, p[2] / norm};
            // The ray keeps its azimuth; the sensors are looked up only on the
            // layers between the own layer's pairs' ends, where a pair reads
            // them.
            const int phiBin = static_cast<int>(phi * Map::kPhiBins / kTwoPi<double>);
            const double axial = std::hypot(u[0], u[1]);
            for (int layer = 0; layer < kL; ++layer) {
              auto& at = crossing[layer];
              at.nominalT = nominal(static_cast<LayerID>(layer), u);
              at.onSensor = false;
              at.first = at.nominalT;
              at.last = at.nominalT;
            }
            double lowest = crossing[own].nominalT;
            double highest = crossing[own].nominalT;
            for (const Half& h : mine) {
              if (crossing[h.other].nominalT > 0.) {
                lowest = std::min(lowest, crossing[h.other].nominalT);
                highest = std::max(highest, crossing[h.other].nominalT);
              }
            }
            for (int layer = 0; layer < kL; ++layer) {
              auto& at = crossing[layer];
              if (!(at.nominalT > 0.) || at.nominalT < lowest || at.nominalT > highest) {
                continue;
              }
              double first = -1.;
              double last = -1.;
              grids[layer]->crossings(u, phiBin, Detector<Det>::isBarrelLayer(static_cast<LayerID>(layer)) ? at.nominalT * u[2] : at.nominalT * axial, first, last);
              if (first > 0.) {
                at.first = first;
                at.last = last;
                at.onSensor = true;
              }
            }
            const Crossing& here = crossing[own];
            const double tOwn = 0.5 * (here.first + here.last);
            // The walk ends, outward and inward of the own layer.
            outward.clear();
            inward.clear();
            for (std::size_t h = 0; h < mine.size(); ++h) {
              const Crossing& there = crossing[mine[h].other];
              mid[h] = -1.;
              if (!(there.nominalT > 0.)) {
                continue;
              }
              mid[h] = 0.5 * (tOwn + (0.5 * (there.first + there.last)));
              (mid[h] > tOwn ? outward : inward).push_back({.t = mid[h], .index = h});
            }
            // Walked from the own layer's crossing nearer each end, and the
            // stretch between its two crossings, half of which each half takes.
            const auto walk = [&](double from, double to) {
              return lut.getMatBudget(static_cast<float>(from * u[0]), static_cast<float>(from * u[1]), static_cast<float>(from * u[2]), static_cast<float>(to * u[0]), static_cast<float>(to * u[1]), static_cast<float>(to * u[2]));
            };
            double betweenX2X0 = 0.;
            double betweenXRho = 0.;
            if (here.last > here.first) {
              const auto between = walk(here.first, here.last);
              betweenX2X0 = between.meanX2X0;
              betweenXRho = between.meanRho * between.length;
            }
            std::ranges::sort(outward, {}, &Cut::t);
            std::ranges::sort(inward, std::ranges::greater{}, &Cut::t);
            for (auto* cuts : {&outward, &inward}) {
              double at = cuts == &outward ? here.last : here.first;
              double x2X0 = 0.;
              double xRho = 0.;
              for (const Cut& cut : *cuts) {
                if (std::abs(cut.t - at) > 0.) {
                  const auto step = walk(at, cut.t);
                  x2X0 += step.meanX2X0;
                  xRho += step.meanRho * step.length;
                  at = cut.t;
                }
                walkedX2X0[cut.index] = x2X0 + (0.5 * betweenX2X0);
                walkedXRho[cut.index] = xRho + (0.5 * betweenXRho);
              }
            }
            for (std::size_t h = 0; h < mine.size(); ++h) {
              if (!(mid[h] > 0.)) {
                continue;
              }
              const int from = mine[h].half == 0 ? own : mine[h].other;
              const int to = mine[h].half == 0 ? mine[h].other : own;
              const double tFrom = 0.5 * (crossing[from].first + crossing[from].last);
              const double tTo = 0.5 * (crossing[to].first + crossing[to].last);
              bool crosses = false;
              for (int layer = 0; layer < kL && !crosses; ++layer) {
                const auto& at = crossing[layer];
                crosses = layer != from && layer != to && at.onSensor && at.nominalT > std::min(tFrom, tTo) && at.nominalT < std::max(tFrom, tTo);
              }
              const bool sensed = !crosses && crossing[from].onSensor && crossing[to].onSensor;
              for (int set = 0; set < 3; ++set) {
                if ((set == 0 && !sensed) || (set == 1 && crosses)) {
                  continue;
                }
                auto& sum = sums[h][set];
                sum.x2X0 += walkedX2X0[h];
                sum.xRho += walkedXRho[h];
                ++sum.n;
              }
            }
          }
        }
        for (std::size_t h = 0; h < mine.size(); ++h) {
          const auto& set = sums[h];
          const Sum& taken = set[0].n > 0 ? set[0] : (set[1].n > 0 ? set[1] : set[2]);
          auto& cell = map.cells[(((static_cast<std::size_t>(mine[h].pair) * 2U) + static_cast<std::size_t>(mine[h].half)) * Map::kCells) + (static_cast<std::size_t>(iPhi) * Map::kSecondBins) + static_cast<std::size_t>(iSecond)];
          cell = taken.n > 0 ? MaterialCell{.x2X0 = static_cast<float>(taken.x2X0 / taken.n), .xRho = static_cast<float>(taken.xRho / taken.n)} : MaterialCell{};
        }
      }
    }
  });
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_SURFACEMATERIALFILL_H_
