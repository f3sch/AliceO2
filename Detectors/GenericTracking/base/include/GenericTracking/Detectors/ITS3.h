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

/// \file ITS3.h
/// \brief ITS3: three truly cylindrical inner layers of MOSAIX tiles and the
///        ITS's four outer-barrel layers of ALPIDE chips, placed by O2.
///
/// A placed layout (PlacedLayerSpec) like ITS.h, placed from O2's
/// o2::its::GeometryTGeo, which builds ITS3 when the loaded geometry holds
/// it. A sensor is O2's chip: on the inner layers a tile's pixel array (a
/// piece of a cylinder about the beam, 2 half-layers x 3, 4 or 5 segments x
/// 12 RSUs x 12 tiles), then the outer barrel's ALPIDEs, numbered after the
/// tiles. Between the tiles of a layer lie their periphery, the gaps between
/// segments and the equatorial gap between the half-layers: no silicon, so a
/// crossing there leaves no cluster.
///
/// The tile is restated from O2's ITS3Base/SpecsV2.h and SegmentationMosaix.h,
/// which exist only with ENABLE_UPGRADES while the base library builds this
/// layout in every build; workflow-its3/src/ITS3.cxx holds them to O2's. The
/// tables are generated: `o2-gtrk-import-geometry -d ITS3` prints the chip
/// counts, `o2-gtrk-layer-links -d ITS3` the links.

#ifndef GTRK_DETECTORS_ITS3_H_
#define GTRK_DETECTORS_ITS3_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/Alpide.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace o2::gtrk::detectors
{

namespace mosaix
{

/// A tile's pixel array: 442 rows across it, 156 columns along the beam.
inline constexpr std::int32_t kRows = 442;
inline constexpr std::int32_t kColumns = 156;
inline constexpr double kWidth = 9.197 * kMm;  ///< across, the arc at kMiddleRadius
inline constexpr double kLength = 3.571 * kMm; ///< along the beam
inline constexpr double kPitchRow = kWidth / static_cast<double>(kRows);
inline constexpr double kPitchColumn = kLength / static_cast<double>(kColumns);

/// The layers' nominal radii, and the silicon about them: 45 um, its inner
/// surface 25 um below the nominal radius (half the silicon and the 5 um
/// metal stack), its outer 20 um above.
inline constexpr std::array<double, 3> kNominalRadius{19.0006 * kMm, 25.228 * kMm, 31.4554 * kMm};
inline constexpr double kSiliconThickness = 45 * kUm;
inline constexpr double kMetalStackThickness = 5 * kUm;
inline constexpr double kInnerOffset = (kSiliconThickness + kMetalStackThickness) / 2.;
inline constexpr double kOuterOffset = (kSiliconThickness - kMetalStackThickness) / 2.;

[[nodiscard]] constexpr double innerRadius(std::int32_t layer) { return kNominalRadius[static_cast<std::uint32_t>(layer)] - kInnerOffset; }
[[nodiscard]] constexpr double outerRadius(std::int32_t layer) { return kNominalRadius[static_cast<std::uint32_t>(layer)] + kOuterOffset; }
/// The radius O2 measures a tile's flat coordinate at, the silicon's middle.
[[nodiscard]] constexpr double middleRadius(std::int32_t layer) { return (innerRadius(layer) + outerRadius(layer)) / 2.; }

/// Where O2's dictionary puts a tile's cluster: this far above the silicon's
/// middle (ITS3Base's APTS response depth), which the MC hits confirm.
inline constexpr double kResponseShift = 15.5 * kUm;

/// A tile layer's measurement depth (Geometry::measurementDepth), from our
/// inner surface [cm].
[[nodiscard]] constexpr double measurementDepth(std::int32_t layer) { return middleRadius(layer) - innerRadius(layer) + kResponseShift; }

/// The angle of the middle of a tile's arc in O2's tile frame, from its
/// local +x axis, where O2's flat coordinate starts at half the width.
[[nodiscard]] constexpr double middleAngle(std::int32_t layer) { return 0.5 * kWidth / middleRadius(layer); }

/// The tile of layer \a layer as our curved sensor: its arc and pitch across
/// at the inner surface, where SensorSpec::radius measures them, scaled
/// from O2's at the middle radius.
[[nodiscard]] constexpr SensorSpec tile(std::int32_t layer)
{
  const double scale = innerRadius(layer) / middleRadius(layer);
  return {
    .shape = ShapeType::kCylinder,
    .halfX = 0.5 * kWidth * scale,
    .halfY = 0.5 * kLength,
    .halfThickness = 0.5 * kSiliconThickness,
    .pitchX = kPitchRow * scale,
    .pitchY = kPitchColumn,
    .radius = innerRadius(layer),
  };
}

} // namespace mosaix

struct ITS3 {
  static constexpr const char* name = "ITS3";
  /// The layers on which a particle can leave two rows, on two overlapping
  /// sensors (LayerLinks.h `overlapLayers`): the outer barrel's. The tiles of
  /// an inner layer abut without overlapping.
  static constexpr std::uint64_t overlapLayers = 0b1111000U;
  // clang-format off
  // Chip counts measured by o2-gtrk-import-geometry from an ITS3 (ALICE2.1) o2sim_geometry-aligned.root: 27144 chips.
  static constexpr std::array<PlacedLayerSpec, 7> barrelLayers{{
    {.sensorCount = 864, .sensor = mosaix::tile(0)},
    {.sensorCount = 1152, .sensor = mosaix::tile(1)},
    {.sensorCount = 1440, .sensor = mosaix::tile(2)},
    {.sensorCount = 2688, .sensor = alpide::kSensor},
    {.sensorCount = 3360, .sensor = alpide::kSensor},
    {.sensorCount = 8232, .sensor = alpide::kSensor},
    {.sensorCount = 9408, .sensor = alpide::kSensor},
  }};
  // Traced by o2-gtrk-layer-links for ITS3: pT >= 0.1 GeV at Bz = 5.00675 kG,
  // beam spot z in [-18, 18] cm, |eta| <= 6, source intervals padded by 0.5 cm.
  // The comment on each link is the fraction of its source layer's crossings that take it.
  static constexpr std::array<LayerLink, 11> layerLinks{{
    {.from = 0, .to = 1, .sourceMin = -13.49, .sourceMax = 13.50}, // 0.9884
    {.from = 0, .to = 3, .sourceMin = -13.50, .sourceMax = 13.50}, // 0.0076
    {.from = 0, .to = 5, .sourceMin = -13.09, .sourceMax = 13.09}, // 0.0040
    {.from = 1, .to = 2, .sourceMin = -13.50, .sourceMax = 13.50}, // 0.9880
    {.from = 1, .to = 3, .sourceMin = -13.50, .sourceMax = 13.50}, // 0.0120
    {.from = 2, .to = 3, .sourceMin = -13.50, .sourceMax = 13.50}, // 0.9943
    {.from = 2, .to = 5, .sourceMin = -12.71, .sourceMax = 12.71}, // 0.0057
    {.from = 3, .to = 4, .sourceMin = -37.84, .sourceMax = 37.84}, // 0.8916
    {.from = 3, .to = 5, .sourceMin = -42.71, .sourceMax = 42.71}, // 0.1084
    {.from = 4, .to = 5, .sourceMin = -42.71, .sourceMax = 42.71}, // 1.0000
    {.from = 5, .to = 6, .sourceMin = -67.33, .sourceMax = 67.33}, // 1.0000
  }};
  // clang-format on
  static constexpr std::array<PlacedLayerSpec, 0> discLayers{};
  static constexpr std::int32_t nSensors = countSensors(barrelLayers, discLayers);
  /// The first outer-barrel chip's ID, after every tile.
  static constexpr std::int32_t kFirstOuterSensor = countSensors(std::array{barrelLayers[0], barrelLayers[1], barrelLayers[2]}, discLayers);
};

} // namespace o2::gtrk::detectors

#endif // GTRK_DETECTORS_ITS3_H_
