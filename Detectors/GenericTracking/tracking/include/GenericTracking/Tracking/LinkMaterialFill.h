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

/// \file LinkMaterialFill.h
/// \brief The seeding's material tables (LinkMaterial.h) filled from the
/// material LUT and the placed layers, on the host once per run and beam
/// line.

#ifndef GTRK_TRACKING_LINKMATERIALFILL_H_
#define GTRK_TRACKING_LINKMATERIALFILL_H_

#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/LinkMaterial.h"
#include "GenericTracking/Tracking/Parallel.h"
#include "GenericTracking/Tracking/SearchWindow.h"
#include "GenericTracking/Tracking/SensorCells.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include "DetectorsBase/MatLayerCylSet.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
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
  /// The grid's cells: a stave's width in azimuth, a sensor's length along the layer.
  static constexpr int kPhi = 256;
  static constexpr int kSecond = 128;

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

namespace detail
{

/// A layer's envelope faces: a barrel layer's least and largest radius, a
/// disc's faces toward and away from the vertex (signed z).
template <DetectorLayout Det>
[[nodiscard]] inline double innerFace(LayerID layer)
{
  const auto& e = PlacedEnvelopes<Det>::get(layer);
  return Detector<Det>::isBarrelLayer(layer) ? e.rMin : (e.position > 0. ? e.zMin : e.zMax);
}
template <DetectorLayout Det>
[[nodiscard]] inline double outerFace(LayerID layer)
{
  const auto& e = PlacedEnvelopes<Det>::get(layer);
  return Detector<Det>::isBarrelLayer(layer) ? e.rMax : (e.position > 0. ? e.zMax : e.zMin);
}

/// Where origin + l u crosses a face of \a layer at \a face (a radius or a
/// signed z): the farther crossing of a cylinder; NaN for none.
template <DetectorLayout Det>
[[nodiscard]] inline double faceCrossing(LayerID layer, double face, const Vec3D& origin, const Vec3D& u)
{
  if (Detector<Det>::isBarrelLayer(layer)) {
    const double a = (u[0] * u[0]) + (u[1] * u[1]);
    const double b = 2. * ((origin[0] * u[0]) + (origin[1] * u[1]));
    const double c = (origin[0] * origin[0]) + (origin[1] * origin[1]) - (face * face);
    const double discriminant = (b * b) - (4. * a * c);
    if (!(a > 0.) || discriminant < 0.) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return (-b + std::sqrt(discriminant)) / (2. * a);
  }
  return u[2] != 0. ? (face - origin[2]) / u[2] : std::numeric_limits<double>::quiet_NaN();
}

/// A layer's point at its nominal position, azimuth \a phi, second coordinate \a second.
template <DetectorLayout Det>
[[nodiscard]] inline Vec3D nominalPoint(LayerID layer, double second, double phi)
{
  const double position = PlacedEnvelopes<Det>::get(layer).position;
  if (Detector<Det>::isBarrelLayer(layer)) {
    return {position * std::cos(phi), position * std::sin(phi), second};
  }
  return {second * std::cos(phi), second * std::sin(phi), position};
}

/// The LUT's x/X0 along origin + l u for l in [l0, l1] in \a pieces, with
/// its x-weighted sums of the fraction l / \a length and its square.
struct RaySums {
  double x = 0.;
  double xf = 0.;
  double xf2 = 0.;
  double xRho = 0.;
};

[[nodiscard]] inline RaySums raySums(const o2::base::MatLayerCylSet& lut, const Vec3D& origin, const Vec3D& u, double l0, double l1, double length, int pieces)
{
  RaySums sums;
  if (!std::isfinite(l0) || !std::isfinite(l1) || !(l1 > l0)) {
    return sums;
  }
  for (int i = 0; i < pieces; ++i) {
    const double a = l0 + ((l1 - l0) * i / pieces);
    const double b = l0 + ((l1 - l0) * (i + 1) / pieces);
    const auto budget = lut.getMatBudget(static_cast<float>(origin[0] + (a * u[0])), static_cast<float>(origin[1] + (a * u[1])), static_cast<float>(origin[2] + (a * u[2])), static_cast<float>(origin[0] + (b * u[0])), static_cast<float>(origin[1] + (b * u[1])), static_cast<float>(origin[2] + (b * u[2])));
    const auto x = static_cast<double>(budget.meanX2X0);
    const double f = 0.5 * (a + b) / length;
    sums.xRho += static_cast<double>(budget.meanRho) * static_cast<double>(budget.length);
    sums.x += x;
    sums.xf += x * f;
    sums.xf2 += x * f * f;
  }
  return sums;
}

[[nodiscard]] inline Vec3D unitFrom(const Vec3D& from, const Vec3D& to, double& length)
{
  const Vec3D d{to[0] - from[0], to[1] - from[1], to[2] - from[2]};
  length = std::sqrt((d[0] * d[0]) + (d[1] * d[1]) + (d[2] * d[2]));
  return {d[0] / length, d[1] / length, d[2] / length};
}

/// The gap of the rays' sums \a total over \a rays rays (those through no
/// material included).
[[nodiscard]] inline MaterialGap gapOf(const RaySums& total, int rays)
{
  if (!(total.x > 0.) || rays == 0) {
    return {};
  }
  const double x = total.x / rays;
  return {.x = static_cast<float>(x), .m1 = static_cast<float>(total.xf / total.x), .m2 = static_cast<float>(total.xf2 / total.x), .logX = static_cast<float>(std::log(x)), .xRho = static_cast<float>(total.xRho / rays)};
}

} // namespace detail

/// \a table from \a lut and the placed layers (PlacedEnvelopes, their sensor
/// planes \a surfaces), the vertices' z over \a vertexZ: every bin's material
/// averaged over rays spread across its azimuth, in parallel on \a arena.
///
/// A node's rays come from the nominal vertex through the layer's nominal
/// surface; its band is cut where the ray crosses a sensor plane of the
/// layer (the nominal surface where it crosses none) and scaled to normal
/// incidence. A gap's rays run between the two layers' nominal surfaces at
/// one azimuth, from where they leave the source's band to where they enter
/// the target's; a source and target whose bands overlap on the ray leave an
/// empty gap. A vertex gap's run from the vertex to the layer's band.
template <DetectorLayout Det>
inline void fillLinkMaterial(LinkMaterial<Det>& table, const o2::base::MatLayerCylSet& lut, std::span<const TrackingSurface> surfaces, Interval vertexZ, detail::TaskArena& arena)
{
  using Table = LinkMaterial<Det>;
  constexpr int kL = Table::kLayers;
  constexpr int kRays = 4;     // per azimuth bin
  constexpr int kPieces = 32;  // per gap ray
  constexpr int kNodePieces = 4; // per side of a node's sensor
  table.linkOf.fill(-1);
  for (int link = 0; link < Table::kLinks; ++link) {
    const auto& l = Det::layerLinks[link];
    table.linkOf[(static_cast<std::size_t>(l.from) * kL) + l.to] = static_cast<std::int16_t>(link);
    table.linkFrom[link] = l.from;
    table.linkTo[link] = l.to;
    table.linkSource[link] = {.low = static_cast<float>(l.sourceMin), .perUnit = static_cast<float>(Table::kSecond / (l.sourceMax - l.sourceMin))};
  }
  for (int layer = 0; layer < kL; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    table.layerSecond[layer] = {.low = static_cast<float>(secondMin<Det>(id)), .perUnit = static_cast<float>(Table::kSecond / (secondMax<Det>(id) - secondMin<Det>(id)))};
  }
  const double zLow = vertexZ.low;
  const double zHigh = vertexZ.high > vertexZ.low ? vertexZ.high : vertexZ.low + 1.;
  table.vertexZ = {.low = static_cast<float>(zLow), .perUnit = static_cast<float>(Table::kVertexZ / (zHigh - zLow))};

  std::vector<std::unique_ptr<SensorPlaneGrid<Det>>> grids;
  grids.reserve(kL);
  for (int layer = 0; layer < kL; ++layer) {
    grids.push_back(std::make_unique<SensorPlaneGrid<Det>>(surfaces, static_cast<LayerID>(layer)));
  }
  const auto centre = [](const MaterialAxis& axis, int bin) { return axis.low + ((bin + 0.5) / axis.perUnit); };
  const auto phiOfRay = [](int bin, int bins, int ray) { return kTwoPi<double> * (bin + ((ray + 0.5) / kRays)) / bins; };

  // Jobs: a layer's nodes and vertex gaps, a link's gaps.
  const std::size_t jobs = static_cast<std::size_t>(kL + Table::kLinks);
  detail::parallelFor(arena, 0U, jobs, 1U, [&](std::size_t first, std::size_t last) {
    for (std::size_t job = first; job < last; ++job) {
      if (job < static_cast<std::size_t>(kL)) {
        const auto id = static_cast<LayerID>(job);
        const bool barrel = Detector<Det>::isBarrelLayer(id);
        for (int phi = 0; phi < Table::kNodePhi; ++phi) {
          for (int s = 0; s < Table::kSecond; ++s) {
            const double second = centre(table.layerSecond[job], s);
            double front = 0.;
            double behind = 0.;
            double xRho = 0.;
            for (int ray = 0; ray < kRays; ++ray) {
              const double angle = phiOfRay(phi, Table::kNodePhi, ray);
              const Vec3D origin{0., 0., 0.};
              double length = 0.;
              const Vec3D u = detail::unitFrom(origin, detail::nominalPoint<Det>(id, second, angle), length);
              const double l0 = detail::faceCrossing<Det>(id, detail::innerFace<Det>(id), origin, u);
              const double l1 = detail::faceCrossing<Det>(id, detail::outerFace<Det>(id), origin, u);
              const double uu[3] = {u[0], u[1], u[2]};
              double firstCrossing = -1.;
              double lastCrossing = -1.;
              grids[job]->crossings(uu, static_cast<int>(angle * SensorPlaneGrid<Det>::kPhi / kTwoPi<double>), second, firstCrossing, lastCrossing);
              const double split = std::clamp(firstCrossing > 0. ? firstCrossing : length, l0, l1);
              // To normal incidence: the ray's path factor through the band.
              const double incidence = barrel ? std::hypot(u[0], u[1]) : std::abs(u[2]);
              const auto inFront = detail::raySums(lut, origin, u, l0, split, 1., kNodePieces);
              const auto inBehind = detail::raySums(lut, origin, u, split, l1, 1., kNodePieces);
              front += inFront.x * incidence;
              behind += inBehind.x * incidence;
              xRho += (inFront.xRho + inBehind.xRho) * incidence;
            }
            const double whole = (front + behind) / kRays;
            table.nodes[(((job * Table::kNodePhi) + static_cast<std::size_t>(phi)) * Table::kSecond) + static_cast<std::size_t>(s)] = {.front = static_cast<float>(front / kRays), .behind = static_cast<float>(behind / kRays), .logX = whole > 0. ? static_cast<float>(std::log(whole)) : 0.F, .xRho = static_cast<float>(xRho / kRays)};
          }
        }
        for (int phi = 0; phi < Table::kGapPhi; ++phi) {
          for (int v = 0; v < Table::kVertexZ; ++v) {
            const double z = centre(table.vertexZ, v);
            for (int s = 0; s < Table::kSecond; ++s) {
              const double second = centre(table.layerSecond[job], s);
              detail::RaySums total;
              for (int ray = 0; ray < kRays; ++ray) {
                const Vec3D origin{0., 0., z};
                double length = 0.;
                const Vec3D u = detail::unitFrom(origin, detail::nominalPoint<Det>(id, second, phiOfRay(phi, Table::kGapPhi, ray)), length);
                const auto sums = detail::raySums(lut, origin, u, 0., detail::faceCrossing<Det>(id, detail::innerFace<Det>(id), origin, u), length, kPieces);
                total.x += sums.x;
                total.xf += sums.xf;
                total.xf2 += sums.xf2;
                total.xRho += sums.xRho;
              }
              table.vertexGaps[(((((job * Table::kGapPhi) + static_cast<std::size_t>(phi)) * Table::kVertexZ) + static_cast<std::size_t>(v)) * Table::kSecond) + static_cast<std::size_t>(s)] = detail::gapOf(total, kRays);
            }
          }
        }
        continue;
      }
      const std::size_t link = job - static_cast<std::size_t>(kL);
      const auto& l = Det::layerLinks[link];
      for (int phi = 0; phi < Table::kGapPhi; ++phi) {
        for (int s = 0; s < Table::kSecond; ++s) {
          const double sourceSecond = centre(table.linkSource[link], s);
          float rowMax = 0.F;
          for (int t = 0; t < Table::kSecond; ++t) {
            const double targetSecond = centre(table.layerSecond[l.to], t);
            detail::RaySums total;
            for (int ray = 0; ray < kRays; ++ray) {
              const double angle = phiOfRay(phi, Table::kGapPhi, ray);
              const Vec3D source = detail::nominalPoint<Det>(l.from, sourceSecond, angle);
              double length = 0.;
              const Vec3D u = detail::unitFrom(source, detail::nominalPoint<Det>(l.to, targetSecond, angle), length);
              const double l0 = detail::faceCrossing<Det>(l.from, detail::outerFace<Det>(l.from), source, u);
              const double l1 = detail::faceCrossing<Det>(l.to, detail::innerFace<Det>(l.to), source, u);
              const auto sums = detail::raySums(lut, source, u, l0, l1, length, kPieces);
              total.x += sums.x;
              total.xf += sums.xf;
              total.xf2 += sums.xf2;
              total.xRho += sums.xRho;
            }
            const auto gap = detail::gapOf(total, kRays);
            table.gaps[(((((link * Table::kGapPhi) + static_cast<std::size_t>(phi)) * Table::kSecond) + static_cast<std::size_t>(s)) * Table::kSecond) + static_cast<std::size_t>(t)] = gap;
            rowMax = std::max(rowMax, gap.x * (1.F - (2.F * gap.m1) + gap.m2));
          }
          table.gapRowMax[(((link * Table::kGapPhi) + static_cast<std::size_t>(phi)) * Table::kSecond) + static_cast<std::size_t>(s)] = rowMax;
        }
      }
    }
  });
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_LINKMATERIALFILL_H_
