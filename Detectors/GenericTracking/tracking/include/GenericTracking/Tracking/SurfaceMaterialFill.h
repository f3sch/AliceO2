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

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Parallel.h"
#include "GenericTracking/Tracking/SearchWindow.h"
#include "GenericTracking/Tracking/SurfaceMaterial.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include "DetectorsBase/MatLayerCylSet.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numbers>
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

  SensorPlaneGrid(std::span<const TrackingSurface> surfaces, LayerID id) : mBarrel(Detector<Det>::isBarrelLayer(id))
  {
    mSecondMin = secondMin<Det>(id);
    mPerCm = kSecond / (secondMax<Det>(id) - mSecondMin);
    mCells.resize(static_cast<std::size_t>(kPhi) * kSecond);
    for (const auto& surface : surfaces) {
      if (surface.layer != id) {
        continue;
      }
      double phiLo = 1e9, phiHi = -1e9, sLo = 1e9, sHi = -1e9;
      const bool curved = surface.radius > 0.F;
      double centre[3];
      for (int k = 0; k < 3; ++k) {
        centre[k] = curved ? surface.origin[k] + (surface.radius * surface.axis0[k]) : surface.origin[k];
      }
      const double reference = std::atan2(centre[1], centre[0]);
      for (int i = -1; i <= 1; i += 2) {
        for (int j = -1; j <= 1; j += 2) {
          double g[3];
          for (int k = 0; k < 3; ++k) {
            if (curved) {
              // A curved sensor: its unrolled x is the arc at its radius about
              // its axis (the normal), x along axis0 at angle 0.
              const double angle = i * surface.bounds.halfX / surface.radius;
              g[k] = surface.origin[k] + (surface.radius * ((std::cos(angle) * surface.axis0[k]) + (std::sin(angle) * surface.axis1[k]))) + (j * surface.bounds.halfY * surface.normal[k]);
            } else {
              g[k] = surface.origin[k] + (i * surface.bounds.halfX * surface.axis0[k]) + (j * surface.bounds.halfY * surface.axis1[k]);
            }
          }
          double phi = std::atan2(g[1], g[0]) - reference;
          phi -= 2. * std::numbers::pi * std::round(phi / (2. * std::numbers::pi));
          phiLo = std::min(phiLo, reference + phi);
          phiHi = std::max(phiHi, reference + phi);
          const double second = mBarrel ? g[2] : std::hypot(g[0], g[1]);
          sLo = std::min(sLo, second);
          sHi = std::max(sHi, second);
        }
      }
      const int s0 = std::max(0, static_cast<int>(std::floor((sLo - mSecondMin) * mPerCm)));
      const int s1 = std::min(kSecond - 1, static_cast<int>(std::floor((sHi - mSecondMin) * mPerCm)));
      const int p0 = static_cast<int>(std::floor(phiLo * kPhi / (2. * std::numbers::pi)));
      const int p1 = static_cast<int>(std::floor(phiHi * kPhi / (2. * std::numbers::pi)));
      for (int ip = p0; ip <= p1; ++ip) {
        const int wrapped = ((ip % kPhi) + kPhi) % kPhi;
        for (int is = s0; is <= s1; ++is) {
          mCells[(static_cast<std::size_t>(wrapped) * kSecond) + is].push_back(&surface);
        }
      }
    }
  }

  [[nodiscard]] double cross(const double u[3], double nominalT, bool last) const
  {
    const double p[3] = {nominalT * u[0], nominalT * u[1], nominalT * u[2]};
    double phi = std::atan2(p[1], p[0]);
    phi = phi < 0. ? phi + (2. * std::numbers::pi) : phi;
    const int ip = std::clamp(static_cast<int>(phi * kPhi / (2. * std::numbers::pi)), 0, kPhi - 1);
    const double second = mBarrel ? p[2] : std::hypot(p[0], p[1]);
    const int is = std::clamp(static_cast<int>((second - mSecondMin) * mPerCm), 0, kSecond - 1);
    double best = -1.;
    for (const TrackingSurface* surface : mCells[(static_cast<std::size_t>(ip) * kSecond) + is]) {
      if (surface->radius > 0.F) {
        // |w - (w.n) n| = radius for w = t u - origin, n the axis.
        const double o[3] = {surface->origin[0], surface->origin[1], surface->origin[2]};
        const double un = (u[0] * surface->normal[0]) + (u[1] * surface->normal[1]) + (u[2] * surface->normal[2]);
        const double on = (o[0] * surface->normal[0]) + (o[1] * surface->normal[1]) + (o[2] * surface->normal[2]);
        double up[3], op[3];
        for (int k = 0; k < 3; ++k) {
          up[k] = u[k] - (un * surface->normal[k]);
          op[k] = o[k] - (on * surface->normal[k]);
        }
        const double qa = (up[0] * up[0]) + (up[1] * up[1]) + (up[2] * up[2]);
        const double qb = -2. * ((up[0] * op[0]) + (up[1] * op[1]) + (up[2] * op[2]));
        const double reach = double(surface->radius) + surface->depth;
        const double qc = (op[0] * op[0]) + (op[1] * op[1]) + (op[2] * op[2]) - (reach * reach);
        const double disc = (qb * qb) - (4. * qa * qc);
        if (qa <= 0. || disc < 0.) {
          continue;
        }
        const double t = (-qb + std::sqrt(disc)) / (2. * qa);
        if (!(t > 0.)) {
          continue;
        }
        double w[3];
        for (int k = 0; k < 3; ++k) {
          w[k] = (t * u[k]) - o[k];
        }
        const double x = (w[0] * surface->axis0[0]) + (w[1] * surface->axis0[1]) + (w[2] * surface->axis0[2]);
        const double y = (w[0] * surface->axis1[0]) + (w[1] * surface->axis1[1]) + (w[2] * surface->axis1[2]);
        const double along = (w[0] * surface->normal[0]) + (w[1] * surface->normal[1]) + (w[2] * surface->normal[2]);
        if (!surface->bounds.isInside(static_cast<Scalar>(surface->radius * std::atan2(y, x)), static_cast<Scalar>(along))) {
          continue;
        }
        if (best <= 0. || (last ? t > best : t < best)) {
          best = t;
        }
        continue;
      }
      const double nu = (surface->normal[0] * u[0]) + (surface->normal[1] * u[1]) + (surface->normal[2] * u[2]);
      if (std::abs(nu) < 1e-9) {
        continue;
      }
      const double t = ((surface->normal[0] * surface->origin[0]) + (surface->normal[1] * surface->origin[1]) + (surface->normal[2] * surface->origin[2])) / nu;
      if (!(t > 0.)) {
        continue;
      }
      double d[3];
      for (int k = 0; k < 3; ++k) {
        d[k] = (t * u[k]) - surface->origin[k];
      }
      const double loc0 = (d[0] * surface->axis0[0]) + (d[1] * surface->axis0[1]) + (d[2] * surface->axis0[2]);
      const double loc1 = (d[0] * surface->axis1[0]) + (d[1] * surface->axis1[1]) + (d[2] * surface->axis1[2]);
      if (!surface->bounds.isInside(static_cast<Scalar>(loc0), static_cast<Scalar>(loc1))) {
        continue;
      }
      if (best <= 0. || (last ? t > best : t < best)) {
        best = t;
      }
    }
    return best;
  }

 private:
  bool mBarrel;
  double mSecondMin = 0.;
  double mPerCm = 0.;
  std::vector<std::vector<const TrackingSurface*>> mCells;
};

/// The linked pairs' material (SurfaceMaterial.h) from the LUT: for each cell
/// of a half's layer, rays from the nominal vertex through points spread over
/// the cell (thin walls fall between cell centres), each walked from where it
/// crosses a sensor plane of the cell's layer to the midpoint between that
/// and its crossing of the pair's other layer -- the sensor planes, not the
/// layers' nominal radii or z, since a stave's own material splits between
/// the steps in and out at its sensor. Where a ray crosses two sensors of a
/// layer, both; where it crosses none, the nominal surface. A step between
/// the pair starts and ends on a sensor and had no row on a layer in between:
/// a cell takes the rays that cross a sensor at both ends and none between,
/// else those that cross none between, else all. The halves in parallel on \a arena's
/// threads, the tracker's.
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
  for (int pair = 0; pair < Map::kPairs; ++pair) {
    const auto& link = Det::layerLinks[pair];
    map.surfaces[pair] = {surfaceOf(link.from), surfaceOf(link.to)};
    map.pairOf[(static_cast<std::size_t>(link.from) * kL) + link.to] = static_cast<std::int16_t>(pair + 1);
    map.pairOf[(static_cast<std::size_t>(link.to) * kL) + link.from] = static_cast<std::int16_t>(-(pair + 1));
  }
  std::vector<std::unique_ptr<SensorPlaneGrid<Det>>> grids;
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
  // On a sensor plane where the ray crosses one (\a onSensor), else the
  // nominal surface.
  const auto onLayer = [&](LayerID id, const double u[3], bool last, bool& onSensor) {
    onSensor = false;
    const double t = nominal(id, u);
    if (!(t > 0.)) {
      return t;
    }
    const double plane = grids[id]->cross(u, t, last);
    onSensor = plane > 0.;
    return onSensor ? plane : t;
  };
  detail::parallelFor(arena, 0U, 2U * Map::kPairs, 1U, [&](std::size_t first, std::size_t end) {
    for (std::size_t job = first; job < end; ++job) {
      const int pair = static_cast<int>(job / 2U);
      const int half = static_cast<int>(job % 2U);
      const auto& link = Det::layerLinks[pair];
      const LayerID own = half == 0 ? link.from : link.to;
      const bool barrel = Detector<Det>::isBarrelLayer(own);
      const auto& surface = map.surfaces[pair][half];
      constexpr int kSubPhi = 4;
      constexpr int kSubSecond = 2;
      for (int iPhi = 0; iPhi < Map::kPhiBins; ++iPhi) {
        for (int iSecond = 0; iSecond < Map::kSecondBins; ++iSecond) {
          // The cell's rays in three sets, the first non-empty taken: those
          // that cross a sensor at both ends and none between, those that
          // cross none between, all of them.
          struct Sum {
            double x2X0 = 0.;
            double xRho = 0.;
            int n = 0;
          };
          Sum sensed;
          Sum clear;
          Sum all;
          for (int a = 0; a < kSubPhi; ++a) {
            const double phi = 2. * std::numbers::pi * (iPhi + ((a + 0.5) / kSubPhi)) / Map::kPhiBins;
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
              const double norm = std::sqrt((p[0] * p[0]) + (p[1] * p[1]) + (p[2] * p[2]));
              const double u[3] = {p[0] / norm, p[1] / norm, p[2] / norm};
              // A layer's sensors may overlap (staggered or overlapping
              // staves): a step may leave or reach either, so the own layer is
              // taken at both its first and its last crossing, the other at
              // their mean.
              bool fromSensor = false;
              bool toSensor = false;
              const double fromFirst = onLayer(link.from, u, false, fromSensor);
              const double fromLast = onLayer(link.from, u, true, fromSensor);
              const double toFirst = onLayer(link.to, u, false, toSensor);
              const double toLast = onLayer(link.to, u, true, toSensor);
              if (!(fromFirst > 0.) || !(toFirst > 0.)) {
                continue;
              }
              const double tFrom = 0.5 * (fromFirst + fromLast);
              const double tTo = 0.5 * (toFirst + toLast);
              double x2X0 = 0.;
              double xRho = 0.;
              for (const double tOwn : {half == 0 ? fromFirst : toFirst, half == 0 ? fromLast : toLast}) {
                const double tMid = 0.5 * (tOwn + (half == 0 ? tTo : tFrom));
                const auto walked = lut.getMatBudget(static_cast<float>(tOwn * u[0]), static_cast<float>(tOwn * u[1]), static_cast<float>(tOwn * u[2]), static_cast<float>(tMid * u[0]), static_cast<float>(tMid * u[1]), static_cast<float>(tMid * u[2]));
                x2X0 += 0.5 * walked.meanX2X0;
                xRho += 0.5 * walked.meanRho * walked.length;
              }
              bool crosses = false;
              for (int other = 0; other < kL && !crosses; ++other) {
                if (other == link.from || other == link.to) {
                  continue;
                }
                const double t = nominal(static_cast<LayerID>(other), u);
                if (t > std::min(tFrom, tTo) && t < std::max(tFrom, tTo)) {
                  crosses = grids[other]->cross(u, t, false) > 0.;
                }
              }
              for (Sum* sum : {&all, crosses ? nullptr : &clear, (crosses || !fromSensor || !toSensor) ? nullptr : &sensed}) {
                if (sum != nullptr) {
                  sum->x2X0 += x2X0;
                  sum->xRho += xRho;
                  ++sum->n;
                }
              }
            }
          }
          auto& cell = map.cells[(((static_cast<std::size_t>(pair) * 2U) + half) * Map::kCells) + (static_cast<std::size_t>(iPhi) * Map::kSecondBins) + iSecond];
          const Sum& taken = sensed.n > 0 ? sensed : (clear.n > 0 ? clear : all);
          cell = taken.n > 0 ? MaterialCell{.x2X0 = static_cast<float>(taken.x2X0 / taken.n), .xRho = static_cast<float>(taken.xRho / taken.n)} : MaterialCell{};
        }
      }
    }
  });
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_SURFACEMATERIALFILL_H_
