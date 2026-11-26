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

/// \file BandMaterial.h
/// \brief The seeding's link material (LinkMaterial.h) of material in the
/// layers' bands alone: the tests' model of it.
///
/// The cuts' oracles and the statistical tests of the cuts need material
/// they can place and count by hand, which a binary LUT cannot give. Host
/// only: the table is filled on the host and read where the tracking reads
/// it.

#ifndef GTRK_TEST_BANDMATERIAL_H_
#define GTRK_TEST_BANDMATERIAL_H_

#include "GenericTracking/Base/LayerEnvelopes.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/LinkMaterial.h"
#include "GenericTracking/Tracking/LinkMaterialFill.h"
#include "GenericTracking/Tracking/SearchWindow.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <utility>

namespace o2::gtrk
{

/// A thin cylinder about the beam that scatters and measures nothing, as a
/// beam pipe: \a x2X0 at normal incidence.
struct ThinCylinder {
  double r = 0.;
  double x2X0 = 0.;
};

/// A link material (LinkMaterial.h) of layers whose envelope bands each hold
/// `x2X0[layer]` at normal incidence, evenly over their depth, and of
/// \a cylinders, and nothing else: what `fillLinkMaterial` makes of a LUT with that material, every
/// crossing of a band along a bin's central ray integrated exactly, the
/// layout being symmetric in azimuth. A skip link's gap holds the bands it
/// passes over, a vertex gap the bands inside the layer. With \a xRho, each
/// band's rho x [g/cm^2] at normal incidence, for the energy loss (none
/// without). The layout's geometry must be loaded (PlacedEnvelopes). On the
/// heap: a layout's table is megabytes.
template <DetectorLayout Det>
[[nodiscard]] std::unique_ptr<LinkMaterial<Det>> bandMaterial(const std::array<double, static_cast<std::size_t>(Detector<Det>::nLayers)>& x2X0, std::span<const ThinCylinder> cylinders = {}, Interval vertexZ = {.low = -20., .high = 20.}, std::span<const double> xRho = {})
{
  using Table = LinkMaterial<Det>;
  constexpr auto kL = static_cast<std::uint32_t>(Table::kLayers);
  auto table = std::make_unique<Table>();
  table->linkOf.fill(-1);
  for (std::uint32_t link = 0; link < Table::kLinks; ++link) {
    const auto& l = Det::layerLinks[link];
    table->linkOf[(std::uint32_t{l.from} * kL) + l.to] = static_cast<std::int16_t>(link);
    table->linkFrom[link] = l.from;
    table->linkTo[link] = l.to;
    table->linkSource[link] = {.low = static_cast<float>(secondMin<Det>(l.from)), .perUnit = static_cast<float>(Table::kSourceSecond / (secondMax<Det>(l.from) - secondMin<Det>(l.from)))};
  }
  for (std::uint32_t layer = 0; layer < kL; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    table->layerSecond[layer] = {.low = static_cast<float>(secondMin<Det>(id)), .perUnit = static_cast<float>(Table::kSecond / (secondMax<Det>(id) - secondMin<Det>(id)))};
    const double inner = std::abs(detail::innerFace<Det>(detail::envelopeBands<Det>(), id));
    const double outer = std::abs(detail::outerFace<Det>(detail::envelopeBands<Det>(), id));
    table->layerDepth[layer] = {.low = static_cast<float>(inner), .perUnit = static_cast<float>(1. / (outer - inner))};
    table->layerSide[layer] = Detector<Det>::isBarrelLayer(id) ? std::int8_t{0} : (PlacedEnvelopes<Det>::get(id).position > 0. ? std::int8_t{1} : std::int8_t{-1});
  }
  table->vertexZ = {.low = static_cast<float>(vertexZ.low), .perUnit = static_cast<float>(Table::kVertexZ / (vertexZ.high - vertexZ.low))};
  const auto centre = [](const MaterialAxis& axis, std::uint32_t bin) { return static_cast<double>(axis.low) + ((bin + 0.5) / static_cast<double>(axis.perUnit)); };

  // The bands' x/X0 along origin + l u for l in [l0, l1], with its sums of
  // l / length and its square weighted by x: per band, the ray's intervals
  // inside the band's hull in r and z.
  const auto sums = [&](const Vec3D& origin, const Vec3D& u, double l0, double l1, double length) {
    detail::RaySums total;
    if (!std::isfinite(l0) || !std::isfinite(l1) || !(l1 > l0)) {
      return total;
    }
    const double a = (u[0] * u[0]) + (u[1] * u[1]);
    const double b = 2. * ((origin[0] * u[0]) + (origin[1] * u[1]));
    const double c = (origin[0] * origin[0]) + (origin[1] * origin[1]);
    // Where r^2(l) <= radius^2: [low, high], empty if low > high.
    const auto within = [&](double radius) {
      if (!(a > 0.)) {
        return c <= radius * radius ? std::pair{-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()} : std::pair{1., 0.};
      }
      const double discriminant = (b * b) - (4. * a * (c - (radius * radius)));
      if (discriminant < 0.) {
        return std::pair{1., 0.};
      }
      const double root = std::sqrt(discriminant);
      return std::pair{(-b - root) / (2. * a), (-b + root) / (2. * a)};
    };
    const auto add = [&](double density, double p, double q) {
      if (q > p) {
        total.x += density * (q - p);
        total.xf += density * ((q * q) - (p * p)) / (2. * length);
        total.xf2 += density * ((q * q * q) - (p * p * p)) / (3. * length * length);
        for (std::uint32_t k = 0; k < kGapProfile; ++k) {
          const double low = std::max(p, length * k / kGapProfile);
          const double high = std::min(q, length * (k + 1) / kGapProfile);
          if (high > low) {
            total.byFraction[k] += density * (high - low);
          }
        }
      }
    };
    for (std::uint32_t layer = 0; layer < kL; ++layer) {
      const auto id = static_cast<LayerID>(layer);
      const auto& e = PlacedEnvelopes<Det>::get(id);
      const double thickness = Detector<Det>::isBarrelLayer(id) ? e.rMax - e.rMin : e.zMax - e.zMin;
      if (!(x2X0[layer] > 0.) || !(thickness > 0.)) {
        continue;
      }
      const double density = x2X0[layer] / thickness;
      double low = l0;
      double high = l1;
      if (!exactlyEqual(u[2], 0.)) {
        const double z0 = (e.zMin - origin[2]) / u[2];
        const double z1 = (e.zMax - origin[2]) / u[2];
        low = std::max(low, std::min(z0, z1));
        high = std::min(high, std::max(z0, z1));
      } else if (origin[2] < e.zMin || origin[2] > e.zMax) {
        continue;
      }
      const auto [outerLow, outerHigh] = within(e.rMax);
      low = std::max(low, outerLow);
      high = std::min(high, outerHigh);
      if (!(high > low)) {
        continue;
      }
      // Less the part inside the band's inner radius.
      const auto [innerLow, innerHigh] = within(e.rMin);
      if (innerHigh > innerLow) {
        add(density, low, std::min(high, innerLow));
        add(density, std::max(low, innerHigh), high);
      } else {
        add(density, low, high);
      }
    }
    for (const auto& cylinder : cylinders) {
      const auto [enter, leave] = within(cylinder.r);
      for (const double l : {enter, leave}) {
        if (l <= l0 || l >= l1) {
          continue;
        }
        // The path factor: one over the cosine of the ray to the cylinder's normal.
        const double radial = ((origin[0] + (l * u[0])) * u[0]) + ((origin[1] + (l * u[1])) * u[1]);
        const double x = cylinder.x2X0 * cylinder.r / std::abs(radial);
        const double f = l / length;
        total.x += x;
        total.xf += x * f;
        total.xf2 += x * f * f;
        total.byFraction[static_cast<std::uint32_t>(std::clamp(static_cast<std::int32_t>(f * kGapProfile), 0, kGapProfile - 1))] += x;
      }
    }
    return total;
  };

  for (std::uint32_t layer = 0; layer < kL; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    MaterialNode node;
    node.x = static_cast<float>(x2X0[layer]);
    node.xRho = xRho.empty() ? 0.F : static_cast<float>(xRho[layer]);
    for (std::int32_t k = 0; k <= kNodeDepths; ++k) {
      node.depth[k] = static_cast<std::uint8_t>(std::lround(255. * k / kNodeDepths));
    }
    for (std::uint32_t bin = 0; bin < std::uint32_t{Table::kNodePhi} * Table::kSecond; ++bin) {
      table->nodes[(layer * Table::kNodePhi * Table::kSecond) + bin] = node;
    }
    for (std::uint32_t v = 0; v < Table::kVertexZ; ++v) {
      for (std::uint32_t s = 0; s < Table::kSecond; ++s) {
        const Vec3D origin{0., 0., centre(table->vertexZ, v)};
        double length = 0.;
        const Vec3D u = detail::unitFrom(origin, detail::nominalPoint<Det>(id, centre(table->layerSecond[layer], s), 0.), length);
        const double cosine = Detector<Det>::isBarrelLayer(id) ? std::hypot(u[0], u[1]) : std::abs(u[2]);
        const auto total = sums(origin, u, 0., detail::faceCrossing<Det>(id, detail::innerFace<Det>(detail::envelopeBands<Det>(), id), origin, u), length).scaled(cosine);
        const auto gap = detail::gapOf(total, 1);
        const auto profile = detail::profileOf(total);
        for (std::uint32_t phi = 0; phi < Table::kGapPhi; ++phi) {
          const auto index = (((((layer * Table::kGapPhi) + phi) * Table::kVertexZ) + v) * Table::kSecond) + s;
          table->vertexGaps[index] = gap;
          table->vertexProfiles[index] = profile;
        }
      }
    }
  }
  for (std::uint32_t link = 0; link < Table::kLinks; ++link) {
    const auto& l = Det::layerLinks[link];
    for (std::uint32_t s = 0; s < Table::kSourceSecond; ++s) {
      const Vec3D source = detail::nominalPoint<Det>(l.from, centre(table->linkSource[link], s), 0.);
      float rowMax = 0.F;
      for (std::uint32_t t = 0; t < Table::kSecond; ++t) {
        double length = 0.;
        const Vec3D u = detail::unitFrom(source, detail::nominalPoint<Det>(l.to, centre(table->layerSecond[l.to], t), 0.), length);
        const double l0 = detail::faceCrossing<Det>(l.from, detail::outerFace<Det>(detail::envelopeBands<Det>(), l.from), source, u);
        const double l1 = detail::faceCrossing<Det>(l.to, detail::innerFace<Det>(detail::envelopeBands<Det>(), l.to), source, u);
        const auto gap = detail::gapOf(sums(source, u, l0, l1, length), 1);
        table->gapMax[link] = std::max(table->gapMax[link], gap.x);
        rowMax = std::max(rowMax, gap.x * (1.F - (2.F * gap.m1) + gap.m2));
        for (std::uint32_t phi = 0; phi < Table::kGapPhi; ++phi) {
          table->gaps[(((((link * Table::kGapPhi) + phi) * Table::kSourceSecond) + s) * Table::kSecond) + t] = gap;
        }
      }
      for (std::uint32_t phi = 0; phi < Table::kGapPhi; ++phi) {
        table->gapRowMax[(((link * Table::kGapPhi) + phi) * Table::kSourceSecond) + s] = rowMax;
      }
    }
  }
  return table;
}

/// The same with \a x2X0 in every layer's band.
template <DetectorLayout Det>
[[nodiscard]] std::unique_ptr<LinkMaterial<Det>> bandMaterial(double x2X0, std::span<const ThinCylinder> cylinders = {}, Interval vertexZ = {.low = -20., .high = 20.})
{
  std::array<double, static_cast<std::size_t>(Detector<Det>::nLayers)> layers;
  layers.fill(x2X0);
  return bandMaterial<Det>(layers, cylinders, vertexZ);
}

/// The share of layer \a layer's band in front of \a point (its depth in the
/// band, the band even across it: `bandMaterial`'s) and the path factor of a
/// step along \a d there, as LinkMaterial::step takes them.
template <DetectorLayout Det>
[[nodiscard]] std::pair<Scalar, Scalar> bandShare(LayerID layer, const Scalar point[3], const Scalar d[3])
{
  const auto& e = PlacedEnvelopes<Det>::get(layer);
  const double length = std::hypot(d[0], d[1], d[2]);
  if (Detector<Det>::isBarrelLayer(layer)) {
    const double r = std::hypot(point[0], point[1]);
    const double cosine = static_cast<double>(std::abs((d[0] * point[0]) + (d[1] * point[1]))) / (r * length);
    return {static_cast<Scalar>(std::clamp((r - e.rMin) / (e.rMax - e.rMin), 0., 1.)), static_cast<Scalar>(1. / cosine)};
  }
  const double depth = (static_cast<double>(std::abs(point[2])) - std::min(std::abs(e.zMin), std::abs(e.zMax))) / (e.zMax - e.zMin);
  return {static_cast<Scalar>(std::clamp(depth, 0., 1.)), static_cast<Scalar>(length / static_cast<double>(std::abs(d[2])))};
}

/// \a table filled from \a medium -- a test's UniformMaterial, anything with
/// a LUT's `getMatBudget` -- as the tracking fills it from its LUT
/// (`fillLinkMaterial`), the vertices over \a vertexZ.
template <DetectorLayout Det, typename Medium>
void fillFromMedium(LinkMaterial<Det>& table, const Medium& medium, Interval vertexZ = {.low = -20., .high = 20.})
{
  detail::TaskArena arena{0};
  fillLinkMaterial<Det>(table, medium, vertexZ, arena);
}

/// The link material without any, made once per layout: the iterations'
/// parameters point at it.
template <DetectorLayout Det>
[[nodiscard]] const LinkMaterial<Det>& noMaterial()
{
  static const auto table = bandMaterial<Det>(0.);
  return *table;
}

} // namespace o2::gtrk

#endif // GTRK_TEST_BANDMATERIAL_H_
