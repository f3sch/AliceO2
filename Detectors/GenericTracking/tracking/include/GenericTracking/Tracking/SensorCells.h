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

/// \file SensorCells.h
/// \brief The sensors of each layer by where they are, for the follow's hole
///        test: which of a layer's sensors
///        a track may cross near a predicted crossing.
///
/// Each layer is cut into cells of azimuth and its second coordinate (z on
/// a barrel layer, r on a disc), and each sensor is listed in every cell its
/// active area's footprint reaches: the bounding box of its outline, which a
/// flat sensor's corners and the points of its edges nearest the beam axis
/// bound, and a curved one's corners. A query over a range of both
/// coordinates visits the sensors of the cells it covers, a superset of
/// those whose footprint it meets (each sensor is listed one cell beyond its
/// footprint too, against the query's rounding). The cells are as wide as the layer's
/// narrowest footprint in either coordinate, so a cell holds a few sensors.
///
/// Built from the tracking geometry when the track model is set; the view is
/// what a device reads.

#ifndef GTRK_TRACKING_SENSORCELLS_H_
#define GTRK_TRACKING_SENSORCELLS_H_

#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include "GPUCommonMath.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace o2::gtrk
{

/// One layer's cells: \a nPhi azimuth cells of `kTwoPi<Scalar> / nPhi` from 0, and
/// \a nSecond cells of \a secondBin from \a secondLow; its cells are the
/// CSR's from \a cellOffset on, the second coordinate's index major.
struct SensorCellLayer {
  Scalar phiBin = 0.F;
  Scalar secondLow = 0.F;
  Scalar secondBin = 0.F;
  std::uint32_t nPhi = 0U;
  std::uint32_t nSecond = 0U;
  std::uint32_t cellOffset = 0U;
};

/// A sensor's footprint about the detector's axis: its azimuth from
/// \a phiLow to \a phiHigh, continuous (\a phiHigh may pass pi), and its
/// second coordinate's range (`SensorCells::footprint`, rounded outwards).
struct SensorFootprint {
  Scalar phiLow = 0.F;
  Scalar phiHigh = 0.F;
  Scalar secondLow = 0.F;
  Scalar secondHigh = 0.F;

  /// Whether it meets the range \a phi +- \a phiReach, \a second +-
  /// \a secondReach.
  [[nodiscard]] GPUhdi() bool meets(Scalar phi, Scalar phiReach, Scalar second, Scalar secondReach) const
  {
    using Math = o2::gpu::CAMath;
    if (secondLow > second + secondReach || secondHigh < second - secondReach) {
      return false;
    }
    if (phiReach >= kPi<Scalar>) {
      return true;
    }
    // Its start about phi in [-pi, pi); past pi its end continues from -pi.
    const Scalar low = (phiLow - phi) - (kTwoPi<Scalar> * Math::Floor(((phiLow - phi) / kTwoPi<Scalar>)+0.5F));
    const Scalar high = low + (phiHigh - phiLow);
    return (low <= phiReach && high >= -phiReach) || (high - kTwoPi<Scalar> >= -phiReach);
  }
};

/// What a stage reads: per layer its cells, then per cell its sensors.
struct SensorCellsView {
  const SensorCellLayer* layers = nullptr;
  const std::uint32_t* offsets = nullptr;
  const std::uint32_t* sensors = nullptr;

  /// Every sensor of \a layer listed in a cell that the range \a phi +-
  /// \a phiReach, \a second +- \a secondReach covers, until \a visit returns
  /// true; whether one did. A sensor listed in several of those cells is
  /// visited once per cell.
  template <typename Visit>
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool anySensor(LayerID layer, Scalar phi, Scalar phiReach, Scalar second, Scalar secondReach, const Visit& visit) const
  {
    using Math = o2::gpu::CAMath;
    const auto& cells = layers[layer];
    if (cells.nPhi == 0U || cells.nSecond == 0U) {
      return false;
    }
    const auto nPhi = static_cast<int>(cells.nPhi);
    const auto nSecond = static_cast<int>(cells.nSecond);
    const int phiFirst = static_cast<int>(Math::Floor((phi - phiReach) / cells.phiBin));
    // At most every azimuth cell once.
    const int phiCount = Math::Min(static_cast<int>(Math::Floor((phi + phiReach) / cells.phiBin)) - phiFirst + 1, nPhi);
    // A range beyond every sensor's footprint meets none; one touching the
    // outermost footprints' edges, at either end, reads the outermost cells,
    // which list them.
    const Scalar secondHigh = cells.secondLow + (static_cast<Scalar>(nSecond) * cells.secondBin);
    if (second + secondReach < cells.secondLow || second - secondReach > secondHigh) {
      return false;
    }
    const int secondFirst = Math::Min(Math::Max(static_cast<int>(Math::Floor((second - secondReach - cells.secondLow) / cells.secondBin)), 0), nSecond - 1);
    const int secondLast = Math::Min(Math::Max(static_cast<int>(Math::Floor((second + secondReach - cells.secondLow) / cells.secondBin)), 0), nSecond - 1);
    for (int s = secondFirst; s <= secondLast; ++s) {
      for (int i = 0; i < phiCount; ++i) {
        const int p = (((phiFirst + i) % nPhi) + nPhi) % nPhi;
        const auto cell = cells.cellOffset + static_cast<std::uint32_t>((s * nPhi) + p);
        for (auto at = offsets[cell]; at < offsets[cell + 1U]; ++at) {
          if (visit(sensors[at])) {
            return true;
          }
        }
      }
    }
    return false;
  }

  /// \a visit(sensor) for every sensor of \a layer listed in a cell that
  /// the range \a phi +- \a phiReach, \a second +- \a secondReach covers,
  /// each once: a sensor is listed in a rectangle of cells, so it is visited
  /// in the first cell of the range's part of that rectangle, the one with
  /// no listing cell before it in either coordinate. A range all round in
  /// azimuth may visit a sensor twice, at its two ends.
  template <typename Visit>
  GPUhdi() void forEachSensorOnce(LayerID layer, Scalar phi, Scalar phiReach, Scalar second, Scalar secondReach, const Visit& visit) const
  {
    using Math = o2::gpu::CAMath;
    const auto& cells = layers[layer];
    if (cells.nPhi == 0U || cells.nSecond == 0U) {
      return;
    }
    const auto nPhi = static_cast<int>(cells.nPhi);
    const auto nSecond = static_cast<int>(cells.nSecond);
    const int phiFirst = static_cast<int>(Math::Floor((phi - phiReach) / cells.phiBin));
    const int phiCount = Math::Min(static_cast<int>(Math::Floor((phi + phiReach) / cells.phiBin)) - phiFirst + 1, nPhi);
    const Scalar secondHigh = cells.secondLow + (static_cast<Scalar>(nSecond) * cells.secondBin);
    if (second + secondReach < cells.secondLow || second - secondReach > secondHigh) {
      return;
    }
    const int secondFirst = Math::Min(Math::Max(static_cast<int>(Math::Floor((second - secondReach - cells.secondLow) / cells.secondBin)), 0), nSecond - 1);
    const int secondLast = Math::Min(Math::Max(static_cast<int>(Math::Floor((second + secondReach - cells.secondLow) / cells.secondBin)), 0), nSecond - 1);
    const auto cellOf = [&](int s, int i) { return cells.cellOffset + static_cast<std::uint32_t>((s * nPhi) + ((((phiFirst + i) % nPhi) + nPhi) % nPhi)); };
    const auto lists = [&](std::uint32_t cell, std::uint32_t sensor) {
      for (auto at = offsets[cell]; at < offsets[cell + 1U]; ++at) {
        if (sensors[at] == sensor) {
          return true;
        }
      }
      return false;
    };
    for (int s = secondFirst; s <= secondLast; ++s) {
      for (int i = 0; i < phiCount; ++i) {
        const auto cell = cellOf(s, i);
        for (auto at = offsets[cell]; at < offsets[cell + 1U]; ++at) {
          const auto sensor = sensors[at];
          if ((s > secondFirst && lists(cellOf(s - 1, i), sensor)) || (i > 0 && lists(cellOf(s, i - 1), sensor))) {
            continue;
          }
          visit(sensor);
        }
      }
    }
  }
};

/// The host's cells of every layer of \a Det, from \a geometry.
template <DetectorLayout Det>
class SensorCells
{
 public:
  static constexpr auto kLayers = static_cast<std::size_t>(Detector<Det>::nLayers);

  SensorCells() = default;
  explicit SensorCells(const TrackingGeometry<Det>& geometry) { build(geometry); }

  [[nodiscard]] SensorCellsView view() const { return {.layers = mLayers.data(), .offsets = mOffsets.data(), .sensors = mSensors.data()}; }
  [[nodiscard]] const std::vector<SensorCellLayer>& layers() const { return mLayers; }
  [[nodiscard]] const std::vector<std::uint32_t>& offsets() const { return mOffsets; }
  [[nodiscard]] const std::vector<std::uint32_t>& sensors() const { return mSensors; }
  [[nodiscard]] bool empty() const { return mLayers.empty(); }
  /// Per sensor, its footprint (`SensorFootprint`).
  [[nodiscard]] const std::vector<SensorFootprint>& footprints() const { return mFootprints; }

  /// A surface's footprint: its azimuth range about \a phiCentre, and its
  /// second coordinate's range.
  struct Footprint {
    double phiCentre = 0.;
    double phiLow = 0.;
    double phiHigh = 0.;
    double secondLow = 0.;
    double secondHigh = 0.;
  };

  /// The footprint of \a surface, a sensor of a barrel layer when \a barrel.
  [[nodiscard]] static Footprint footprint(const TrackingSurface& surface, bool barrel)
  {
    std::vector<Vec3D> points;
    const auto add = [&](const Vec3D& point) { points.push_back(point); };
    if (surface.radius > 0.F) {
      // A curved sensor: its unrolled corners on the cylinder at its
      // measurement depth, and the middle of each long edge.
      const auto halfX = static_cast<double>(surface.bounds.halfX);
      for (const double u : {-halfX, 0., halfX}) {
        const double angle = u / static_cast<double>(surface.radius);
        for (const double v : {-static_cast<double>(surface.bounds.halfY), static_cast<double>(surface.bounds.halfY)}) {
          Vec3D point{};
          for (std::size_t i = 0; i < 3; ++i) {
            const double radial = (std::cos(angle) * static_cast<double>(surface.axis0[i])) + (std::sin(angle) * static_cast<double>(surface.axis1[i]));
            point[i] = static_cast<double>(surface.origin[i]) + (static_cast<double>(surface.radius + surface.depth) * radial) + (v * static_cast<double>(surface.normal[i]));
          }
          add(point);
        }
      }
    } else {
      const auto at = [&](double u, double v) {
        Vec3D point{};
        surfaceToGlobal(surface, u, v, point.data());
        return point;
      };
      const auto halfY = static_cast<double>(surface.bounds.halfY);
      const std::array<Vec3D, 4> corners = {at(-static_cast<double>(surface.bounds.halfX), -halfY), at(static_cast<double>(surface.bounds.halfX), -halfY), at(static_cast<double>(surface.bounds.halfX), halfY), at(-static_cast<double>(surface.bounds.halfX), halfY)};
      for (std::size_t i = 0; i < corners.size(); ++i) {
        const auto& a = corners[i];
        const auto& b = corners[(i + 1U) % corners.size()];
        add(a);
        // The edge's point nearest the beam axis: the least radius on it.
        const double dx = b[0] - a[0];
        const double dy = b[1] - a[1];
        const double length2 = (dx * dx) + (dy * dy);
        if (length2 > 0.) {
          const double t = std::clamp(-((a[0] * dx) + (a[1] * dy)) / length2, 0., 1.);
          add({a[0] + (t * dx), a[1] + (t * dy), a[2] + (t * (b[2] - a[2]))});
        }
      }
    }
    Vec3D centre{};
    for (const auto& point : points) {
      for (std::size_t i = 0; i < 3; ++i) {
        centre[i] += point[i] / static_cast<double>(points.size());
      }
    }
    Footprint result;
    result.phiCentre = std::atan2(centre[1], centre[0]);
    result.phiLow = std::numeric_limits<double>::infinity();
    result.phiHigh = -std::numeric_limits<double>::infinity();
    result.secondLow = std::numeric_limits<double>::infinity();
    result.secondHigh = -std::numeric_limits<double>::infinity();
    for (const auto& point : points) {
      const double phi = toPMPiOf(std::atan2(point[1], point[0]) - result.phiCentre);
      const double second = barrel ? point[2] : std::hypot(point[0], point[1]);
      result.phiLow = std::min(result.phiLow, phi);
      result.phiHigh = std::max(result.phiHigh, phi);
      result.secondLow = std::min(result.secondLow, second);
      result.secondHigh = std::max(result.secondHigh, second);
    }
    return result;
  }

 private:
  void build(const TrackingGeometry<Det>& geometry)
  {
    const auto& surfaces = geometry.surfaces();
    std::vector<Footprint> footprints(surfaces.size());
    for (std::size_t id = 0; id < surfaces.size(); ++id) {
      footprints[id] = footprint(surfaces[id], Detector<Det>::isBarrelLayer(surfaces[id].layer));
    }
    const auto down = [](double value) {
      const auto rounded = static_cast<Scalar>(value);
      return static_cast<double>(rounded) > value ? std::nextafter(rounded, -std::numeric_limits<Scalar>::infinity()) : rounded;
    };
    const auto up = [](double value) {
      const auto rounded = static_cast<Scalar>(value);
      return static_cast<double>(rounded) < value ? std::nextafter(rounded, std::numeric_limits<Scalar>::infinity()) : rounded;
    };
    mFootprints.resize(surfaces.size());
    for (std::size_t id = 0; id < surfaces.size(); ++id) {
      const auto& at = footprints[id];
      mFootprints[id] = {.phiLow = down(at.phiCentre + at.phiLow), .phiHigh = up(at.phiCentre + at.phiHigh), .secondLow = down(at.secondLow), .secondHigh = up(at.secondHigh)};
    }
    mLayers.assign(kLayers, SensorCellLayer{});
    std::vector<std::vector<std::uint32_t>> lists;
    for (std::size_t layer = 0; layer < kLayers; ++layer) {
      double phiWidth = kTwoPi<double>;
      double secondWidth = std::numeric_limits<double>::infinity();
      double low = std::numeric_limits<double>::infinity();
      double high = -std::numeric_limits<double>::infinity();
      for (std::size_t id = 0; id < surfaces.size(); ++id) {
        if (surfaces[id].layer != layer) {
          continue;
        }
        const auto& at = footprints[id];
        phiWidth = std::min(phiWidth, at.phiHigh - at.phiLow);
        secondWidth = std::min(secondWidth, at.secondHigh - at.secondLow);
        low = std::min(low, at.secondLow);
        high = std::max(high, at.secondHigh);
      }
      auto& cells = mLayers[layer];
      cells.cellOffset = static_cast<std::uint32_t>(lists.size());
      if (!(low <= high) || !(phiWidth > 0.) || !(secondWidth > 0.)) {
        continue;
      }
      cells.nPhi = static_cast<std::uint32_t>(std::max(1., std::floor(kTwoPi<double> / phiWidth)));
      cells.phiBin = static_cast<Scalar>(kTwoPi<double> / cells.nPhi);
      cells.nSecond = static_cast<std::uint32_t>(std::max(1., std::ceil((high - low) / secondWidth)));
      cells.secondLow = static_cast<Scalar>(low);
      cells.secondBin = static_cast<Scalar>((high - low) / cells.nSecond);
      lists.resize(lists.size() + (static_cast<std::size_t>(cells.nPhi) * cells.nSecond));
      const auto nPhi = static_cast<long>(cells.nPhi);
      for (std::size_t id = 0; id < surfaces.size(); ++id) {
        if (surfaces[id].layer != layer) {
          continue;
        }
        const auto& at = footprints[id];
        // The cells the footprint reaches and one more on each side, so that
        // a query's float rounding at a cell's edge cannot miss it.
        const auto phiFirst = static_cast<long>(std::floor((at.phiCentre + at.phiLow) / static_cast<double>(cells.phiBin))) - 1L;
        const auto phiLast = static_cast<long>(std::floor((at.phiCentre + at.phiHigh) / static_cast<double>(cells.phiBin))) + 1L;
        const auto bin = [&](double second) { return std::clamp(static_cast<long>(std::floor((second - static_cast<double>(cells.secondLow)) / static_cast<double>(cells.secondBin))), 0L, static_cast<long>(cells.nSecond) - 1L); };
        for (auto s = std::max(bin(at.secondLow) - 1L, 0L); s <= std::min(bin(at.secondHigh) + 1L, static_cast<long>(cells.nSecond) - 1L); ++s) {
          for (auto p = phiFirst; p <= std::min(phiLast, phiFirst + nPhi - 1L); ++p) {
            const auto cell = static_cast<std::size_t>(cells.cellOffset) + static_cast<std::size_t>((s * nPhi) + (((p % nPhi) + nPhi) % nPhi));
            lists[cell].push_back(static_cast<std::uint32_t>(id));
          }
        }
      }
    }
    mOffsets.assign(lists.size() + 1U, 0U);
    for (std::size_t cell = 0; cell < lists.size(); ++cell) {
      mOffsets[cell + 1U] = mOffsets[cell] + static_cast<std::uint32_t>(lists[cell].size());
    }
    mSensors.clear();
    mSensors.reserve(mOffsets.back());
    for (const auto& list : lists) {
      mSensors.insert(mSensors.end(), list.begin(), list.end());
    }
  }

  std::vector<SensorCellLayer> mLayers;
  std::vector<std::uint32_t> mOffsets;
  std::vector<std::uint32_t> mSensors;
  std::vector<SensorFootprint> mFootprints;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_SENSORCELLS_H_
