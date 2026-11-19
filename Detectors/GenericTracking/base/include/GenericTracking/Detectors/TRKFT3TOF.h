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

/// \file TRKFT3TOF.h
/// \brief TRKFT3 with ALICE 3's inner and outer TOF as two inactive barrel
///        layers: the TOF clusters found, owned and fitted by the tracking.
///
/// TRKFT3's layers in their original order, the iTOF (r 21.4 cm)
/// inserted as barrel layer 7 between ML layers 6 and 7, and the oTOF (r 92 cm)
/// as barrel layer 12 after the last OT layer: TRKFT3's layers 7-10 are 8-11
/// here, its discs 11-22 are 13-24. Both TOF layers are inactive (LayerSpec.h).
/// Missing hits consume zero holes, and the TRK material slab boundaries
/// stay in place. The tracking schedule selects the iTOF seeding passes.
/// The default seeding mask leaves oTOF to the follower.
/// Sensor IDs are TRKFT3's, then the IOTOF chips at `kFirstTOFSensor` plus
/// O2's IOTOF chip ID (iTOF chips first).
///
/// The TRK links are TRKFT3's renumbered: ML 6 -> ML 7 stays a direct link,
/// beside ML 6 -> iTOF -> ML 7. The links into and out of the TOF layers are
/// traced by `o2-gtrk-layer-links -d TRKFT3TOF`. The passives are TRKFT3's:
/// the iTOF's material stays passive 7's.
///
/// A TOF cluster's time box is the bunch crossing that ends at its time,
/// [t - 1 BC, t] (WorkflowTRK/TRK.h, `tofClusters`): two TOF rows of a track
/// lie within 1 BC of each other, and the box never narrows a TRK readout
/// frame by more than the flight time to the TOF.

#ifndef GTRK_DETECTORS_TRKFT3TOF_H_
#define GTRK_DETECTORS_TRKFT3TOF_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/TRKFT3.h"

#include <array>
#include <cstdint>

namespace o2::gtrk::detectors
{

namespace trkft3tof
{

/// An IOTOF chip: O2's `ChipSpecificsParam` active matrix, 240 rows of
/// 100 um by 120 columns of 250 um, its frame centred on it; the plane at the
/// sensor's mid-plane (measurement depth 0).
inline constexpr SensorSpec kTOFChip{
  .shape = ShapeType::kRectangle,
  .halfX = 1.2,
  .halfY = 1.5,
  .halfThickness = 0.0025,
  .pitchX = 0.01,
  .pitchY = 0.025,
};
inline constexpr SensorSpec kInnerTOFChip = kTOFChip;
inline constexpr SensorSpec kOuterTOFChip = kTOFChip;

/// TRKFT3's layer \a layer in TRKFT3TOF: 0-6 stay, 7-10 move past the iTOF,
/// the discs past the oTOF too.
[[nodiscard]] constexpr LayerID fromTRKFT3(LayerID layer) { return layer <= 6 ? layer : static_cast<LayerID>(layer <= 10 ? layer + 1 : layer + 2); }

/// TRKFT3's overlap layers (LayerLinks.h `overlapLayers`) renumbered; none on
/// the TOF.
inline constexpr std::uint64_t kOverlapLayers = [] {
  std::uint64_t layers = 0U;
  for (int layer = 0; layer < 64; ++layer) {
    if (((TRKFT3::overlapLayers >> layer) & 1U) != 0U) {
      layers |= std::uint64_t{1} << fromTRKFT3(static_cast<LayerID>(layer));
    }
  }
  return layers;
}();

} // namespace trkft3tof

struct TRKFT3TOF {
  static constexpr const char* name = "TRKFT3TOF";
  static constexpr LayerID kInnerTOFLayer = 7;
  static constexpr LayerID kOuterTOFLayer = 12;
  static constexpr std::uint64_t inactiveLayers = (std::uint64_t{1} << kInnerTOFLayer) | (std::uint64_t{1} << kOuterTOFLayer);
  /// The TOF's clusters carry their own time box, [t - 1 BC, t]
  /// (WorkflowTRK/TRKTOF.h), not a readout frame's record.
  static constexpr std::uint64_t hitTimedLayers = inactiveLayers;
  static constexpr std::uint64_t defaultSeedingMask = ~(std::uint64_t{1} << kOuterTOFLayer);
  /// The layers with overlapping sensors: TRKFT3's (`trkft3tof::kOverlapLayers`).
  static constexpr std::uint64_t overlapLayers = trkft3tof::kOverlapLayers;
  // clang-format off
  static constexpr std::array<PlacedLayerSpec, 13> barrelLayers{{
    {.sensorCount = 4, .sensor = trkft3::vdSensor(0)},
    {.sensorCount = 4, .sensor = trkft3::vdSensor(1)},
    {.sensorCount = 4, .sensor = trkft3::vdSensor(2)},
    {.sensorCount = 880, .sensor = trkft3::kChip},
    {.sensorCount = 1232, .sensor = trkft3::kChip},
    {.sensorCount = 1584, .sensor = trkft3::kChip},
    {.sensorCount = 2288, .sensor = trkft3::kChip},
    {.sensorCount = 2352, .sensor = trkft3tof::kInnerTOFChip},
    {.sensorCount = 3344, .sensor = trkft3::kChip},
    {.sensorCount = 10560, .sensor = trkft3::kChip},
    {.sensorCount = 13376, .sensor = trkft3::kChip},
    {.sensorCount = 17600, .sensor = trkft3::kChip},
    {.sensorCount = 56544, .sensor = trkft3tof::kOuterTOFChip},
  }};
  static constexpr std::array<PlacedLayerSpec, 12> discLayers{{
    {.sensorCount = 440, .sensor = trkft3::kChip},
    {.sensorCount = 440, .sensor = trkft3::kChip},
    {.sensorCount = 440, .sensor = trkft3::kChip},
    {.sensorCount = 1864, .sensor = trkft3::kChip},
    {.sensorCount = 1864, .sensor = trkft3::kChip},
    {.sensorCount = 1864, .sensor = trkft3::kChip},
    {.sensorCount = 440, .sensor = trkft3::kChip},
    {.sensorCount = 440, .sensor = trkft3::kChip},
    {.sensorCount = 440, .sensor = trkft3::kChip},
    {.sensorCount = 1864, .sensor = trkft3::kChip},
    {.sensorCount = 1864, .sensor = trkft3::kChip},
    {.sensorCount = 1864, .sensor = trkft3::kChip},
  }};
  static constexpr std::array<LayerLink, 170> layerLinks{{
    {.from = 0, .to = 1, .sourceMin = -21.42, .sourceMax = 21.42, .skipsLayer = false}, // 0.9664
    {.from = 0, .to = 3, .sourceMin = -21.75, .sourceMax = 21.75, .skipsLayer = false}, // 0.0041
    {.from = 0, .to = 14, .sourceMin = -22.64, .sourceMax = -18.67, .skipsLayer = false}, // 0.0023
    {.from = 0, .to = 15, .sourceMin = -23.75, .sourceMax = -16.92, .skipsLayer = false}, // 0.0067
    {.from = 0, .to = 16, .sourceMin = -22.73, .sourceMax = -18.67, .skipsLayer = true}, // 0.0001
    {.from = 0, .to = 17, .sourceMin = -23.69, .sourceMax = -18.67, .skipsLayer = true}, // 0.0009
    {.from = 0, .to = 18, .sourceMin = -24.96, .sourceMax = -15.17, .skipsLayer = false}, // 0.0046
    {.from = 0, .to = 20, .sourceMin = 18.67, .sourceMax = 22.64, .skipsLayer = false}, // 0.0023
    {.from = 0, .to = 21, .sourceMin = 16.92, .sourceMax = 23.75, .skipsLayer = false}, // 0.0067
    {.from = 0, .to = 22, .sourceMin = 18.67, .sourceMax = 22.73, .skipsLayer = true}, // 0.0001
    {.from = 0, .to = 23, .sourceMin = 18.67, .sourceMax = 23.69, .skipsLayer = true}, // 0.0009
    {.from = 0, .to = 24, .sourceMin = 15.17, .sourceMax = 24.96, .skipsLayer = false}, // 0.0046
    {.from = 1, .to = 2, .sourceMin = -21.86, .sourceMax = 21.86, .skipsLayer = false}, // 0.8459
    {.from = 1, .to = 3, .sourceMin = -25.50, .sourceMax = 25.50, .skipsLayer = false}, // 0.0796
    {.from = 1, .to = 14, .sourceMin = -25.50, .sourceMax = -11.50, .skipsLayer = false}, // 0.0069
    {.from = 1, .to = 15, .sourceMin = -25.50, .sourceMax = -8.38, .skipsLayer = false}, // 0.0183
    {.from = 1, .to = 16, .sourceMin = -25.50, .sourceMax = -13.05, .skipsLayer = true}, // 0.0004
    {.from = 1, .to = 17, .sourceMin = -25.50, .sourceMax = -9.93, .skipsLayer = true}, // 0.0021
    {.from = 1, .to = 18, .sourceMin = -25.50, .sourceMax = -6.82, .skipsLayer = false}, // 0.0095
    {.from = 1, .to = 20, .sourceMin = 11.50, .sourceMax = 25.50, .skipsLayer = false}, // 0.0069
    {.from = 1, .to = 21, .sourceMin = 8.38, .sourceMax = 25.50, .skipsLayer = false}, // 0.0183
    {.from = 1, .to = 22, .sourceMin = 13.05, .sourceMax = 25.50, .skipsLayer = true}, // 0.0004
    {.from = 1, .to = 23, .sourceMin = 9.93, .sourceMax = 25.50, .skipsLayer = true}, // 0.0021
    {.from = 1, .to = 24, .sourceMin = 6.82, .sourceMax = 25.50, .skipsLayer = false}, // 0.0095
    {.from = 2, .to = 3, .sourceMin = -25.50, .sourceMax = 25.50, .skipsLayer = false}, // 0.9432
    {.from = 2, .to = 14, .sourceMin = -25.50, .sourceMax = -10.19, .skipsLayer = false}, // 0.0031
    {.from = 2, .to = 15, .sourceMin = -25.50, .sourceMax = -10.19, .skipsLayer = false}, // 0.0165
    {.from = 2, .to = 16, .sourceMin = -25.50, .sourceMax = -18.11, .skipsLayer = true}, // 0.0001
    {.from = 2, .to = 17, .sourceMin = -25.50, .sourceMax = -10.19, .skipsLayer = true}, // 0.0016
    {.from = 2, .to = 18, .sourceMin = -25.50, .sourceMax = -10.19, .skipsLayer = false}, // 0.0071
    {.from = 2, .to = 20, .sourceMin = 10.19, .sourceMax = 25.50, .skipsLayer = false}, // 0.0031
    {.from = 2, .to = 21, .sourceMin = 10.19, .sourceMax = 25.50, .skipsLayer = false}, // 0.0165
    {.from = 2, .to = 22, .sourceMin = 18.11, .sourceMax = 25.50, .skipsLayer = true}, // 0.0001
    {.from = 2, .to = 23, .sourceMin = 10.19, .sourceMax = 25.50, .skipsLayer = true}, // 0.0016
    {.from = 2, .to = 24, .sourceMin = 10.19, .sourceMax = 25.50, .skipsLayer = false}, // 0.0071
    {.from = 3, .to = 4, .sourceMin = -54.45, .sourceMax = 54.45, .skipsLayer = false}, // 0.8386
    {.from = 3, .to = 13, .sourceMin = -60.88, .sourceMax = -45.34, .skipsLayer = false}, // 0.0100
    {.from = 3, .to = 14, .sourceMin = -64.49, .sourceMax = -45.34, .skipsLayer = false}, // 0.0484
    {.from = 3, .to = 15, .sourceMin = -64.49, .sourceMax = -45.34, .skipsLayer = true}, // 0.0189
    {.from = 3, .to = 16, .sourceMin = -64.49, .sourceMax = -45.46, .skipsLayer = true}, // 0.0005
    {.from = 3, .to = 17, .sourceMin = -64.49, .sourceMax = -45.37, .skipsLayer = true}, // 0.0009
    {.from = 3, .to = 18, .sourceMin = -64.49, .sourceMax = -45.54, .skipsLayer = true}, // 0.0020
    {.from = 3, .to = 19, .sourceMin = 45.34, .sourceMax = 60.88, .skipsLayer = false}, // 0.0100
    {.from = 3, .to = 20, .sourceMin = 45.34, .sourceMax = 64.49, .skipsLayer = false}, // 0.0484
    {.from = 3, .to = 21, .sourceMin = 45.34, .sourceMax = 64.49, .skipsLayer = true}, // 0.0189
    {.from = 3, .to = 22, .sourceMin = 45.46, .sourceMax = 64.49, .skipsLayer = true}, // 0.0005
    {.from = 3, .to = 23, .sourceMin = 45.37, .sourceMax = 64.49, .skipsLayer = true}, // 0.0009
    {.from = 3, .to = 24, .sourceMin = 45.54, .sourceMax = 64.49, .skipsLayer = true}, // 0.0020
    {.from = 4, .to = 5, .sourceMin = -53.12, .sourceMax = 53.12, .skipsLayer = false}, // 0.7777
    {.from = 4, .to = 13, .sourceMin = -64.49, .sourceMax = -42.63, .skipsLayer = false}, // 0.0820
    {.from = 4, .to = 14, .sourceMin = -64.49, .sourceMax = -42.63, .skipsLayer = true}, // 0.0279
    {.from = 4, .to = 15, .sourceMin = -64.49, .sourceMax = -42.67, .skipsLayer = true}, // 0.0012
    {.from = 4, .to = 16, .sourceMin = -64.49, .sourceMax = -43.27, .skipsLayer = true}, // 0.0000
    {.from = 4, .to = 17, .sourceMin = -64.49, .sourceMax = -43.40, .skipsLayer = true}, // 0.0000
    {.from = 4, .to = 18, .sourceMin = -64.48, .sourceMax = -48.06, .skipsLayer = true}, // 0.0000
    {.from = 4, .to = 19, .sourceMin = 42.63, .sourceMax = 64.49, .skipsLayer = false}, // 0.0820
    {.from = 4, .to = 20, .sourceMin = 42.63, .sourceMax = 64.49, .skipsLayer = true}, // 0.0279
    {.from = 4, .to = 21, .sourceMin = 42.67, .sourceMax = 64.49, .skipsLayer = true}, // 0.0012
    {.from = 4, .to = 22, .sourceMin = 43.27, .sourceMax = 64.49, .skipsLayer = true}, // 0.0000
    {.from = 4, .to = 23, .sourceMin = 43.40, .sourceMax = 64.49, .skipsLayer = true}, // 0.0000
    {.from = 4, .to = 24, .sourceMin = 48.06, .sourceMax = 64.48, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 6, .sourceMin = -46.99, .sourceMax = 46.99, .skipsLayer = false}, // 0.5203
    {.from = 5, .to = 9, .sourceMin = -48.04, .sourceMax = 48.04, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 12, .sourceMin = -59.09, .sourceMax = 59.09, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 13, .sourceMin = -64.49, .sourceMax = -30.29, .skipsLayer = false}, // 0.2372
    {.from = 5, .to = 14, .sourceMin = -64.49, .sourceMax = -30.32, .skipsLayer = true}, // 0.0026
    {.from = 5, .to = 15, .sourceMin = -64.49, .sourceMax = -31.27, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 16, .sourceMin = -64.48, .sourceMax = -32.10, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 17, .sourceMin = -63.91, .sourceMax = -34.79, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 18, .sourceMin = -58.92, .sourceMax = -32.99, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 19, .sourceMin = 30.29, .sourceMax = 64.49, .skipsLayer = false}, // 0.2372
    {.from = 5, .to = 20, .sourceMin = 30.32, .sourceMax = 64.49, .skipsLayer = true}, // 0.0026
    {.from = 5, .to = 21, .sourceMin = 31.27, .sourceMax = 64.49, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 22, .sourceMin = 32.10, .sourceMax = 64.48, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 23, .sourceMin = 34.79, .sourceMax = 63.91, .skipsLayer = true}, // 0.0000
    {.from = 5, .to = 24, .sourceMin = 32.99, .sourceMax = 58.92, .skipsLayer = true}, // 0.0000
    {.from = 6, .to = 7, .sourceMin = -59.76, .sourceMax = 59.76, .skipsLayer = false}, // 0.7981
    {.from = 6, .to = 8, .sourceMin = -47.81, .sourceMax = 47.81, .skipsLayer = false}, // 0.1654
    {.from = 6, .to = 9, .sourceMin = -64.49, .sourceMax = 64.49, .skipsLayer = false}, // 0.1386
    {.from = 6, .to = 12, .sourceMin = -64.46, .sourceMax = 64.46, .skipsLayer = true}, // 0.0000
    {.from = 6, .to = 13, .sourceMin = -64.49, .sourceMax = -25.80, .skipsLayer = false}, // 0.3430
    {.from = 6, .to = 14, .sourceMin = -64.49, .sourceMax = -26.95, .skipsLayer = true}, // 0.0011
    {.from = 6, .to = 15, .sourceMin = -64.45, .sourceMax = -37.92, .skipsLayer = true}, // 0.0000
    {.from = 6, .to = 16, .sourceMin = -64.49, .sourceMax = -30.59, .skipsLayer = false}, // 0.0038
    {.from = 6, .to = 17, .sourceMin = -64.48, .sourceMax = -33.05, .skipsLayer = true}, // 0.0001
    {.from = 6, .to = 18, .sourceMin = -64.46, .sourceMax = -36.83, .skipsLayer = true}, // 0.0000
    {.from = 6, .to = 19, .sourceMin = 25.80, .sourceMax = 64.49, .skipsLayer = false}, // 0.3430
    {.from = 6, .to = 20, .sourceMin = 26.95, .sourceMax = 64.49, .skipsLayer = true}, // 0.0011
    {.from = 6, .to = 21, .sourceMin = 37.92, .sourceMax = 64.45, .skipsLayer = true}, // 0.0000
    {.from = 6, .to = 22, .sourceMin = 30.59, .sourceMax = 64.49, .skipsLayer = false}, // 0.0038
    {.from = 6, .to = 23, .sourceMin = 33.05, .sourceMax = 64.48, .skipsLayer = true}, // 0.0001
    {.from = 6, .to = 24, .sourceMin = 36.83, .sourceMax = 64.46, .skipsLayer = true}, // 0.0000
    {.from = 7, .to = 8, .sourceMin = -50.76, .sourceMax = 50.76, .skipsLayer = false}, // 0.2092
    {.from = 7, .to = 9, .sourceMin = -63.91, .sourceMax = 63.91, .skipsLayer = false}, // 0.1751
    {.from = 7, .to = 12, .sourceMin = -63.91, .sourceMax = 63.91, .skipsLayer = true}, // 0.0000
    {.from = 7, .to = 13, .sourceMin = -63.91, .sourceMax = -31.08, .skipsLayer = false}, // 0.3028
    {.from = 7, .to = 14, .sourceMin = -63.91, .sourceMax = -32.29, .skipsLayer = true}, // 0.0005
    {.from = 7, .to = 15, .sourceMin = -63.91, .sourceMax = -44.64, .skipsLayer = true}, // 0.0000
    {.from = 7, .to = 16, .sourceMin = -63.91, .sourceMax = -35.97, .skipsLayer = false}, // 0.0045
    {.from = 7, .to = 17, .sourceMin = -63.91, .sourceMax = -38.60, .skipsLayer = true}, // 0.0001
    {.from = 7, .to = 18, .sourceMin = -63.91, .sourceMax = -42.70, .skipsLayer = true}, // 0.0000
    {.from = 7, .to = 19, .sourceMin = 31.08, .sourceMax = 63.92, .skipsLayer = false}, // 0.3028
    {.from = 7, .to = 20, .sourceMin = 32.29, .sourceMax = 63.92, .skipsLayer = true}, // 0.0005
    {.from = 7, .to = 21, .sourceMin = 44.64, .sourceMax = 63.91, .skipsLayer = true}, // 0.0000
    {.from = 7, .to = 22, .sourceMin = 35.97, .sourceMax = 63.92, .skipsLayer = false}, // 0.0045
    {.from = 7, .to = 23, .sourceMin = 38.60, .sourceMax = 63.91, .skipsLayer = true}, // 0.0001
    {.from = 7, .to = 24, .sourceMin = 42.70, .sourceMax = 63.91, .skipsLayer = true}, // 0.0000
    {.from = 8, .to = 9, .sourceMin = -64.49, .sourceMax = 64.49, .skipsLayer = false}, // 0.9844
    {.from = 8, .to = 13, .sourceMin = -64.49, .sourceMax = -47.69, .skipsLayer = false}, // 0.0068
    {.from = 8, .to = 16, .sourceMin = -64.48, .sourceMax = -59.59, .skipsLayer = false}, // 0.0010
    {.from = 8, .to = 19, .sourceMin = 47.69, .sourceMax = 64.49, .skipsLayer = false}, // 0.0068
    {.from = 8, .to = 22, .sourceMin = 59.59, .sourceMax = 64.48, .skipsLayer = false}, // 0.0010
    {.from = 9, .to = 10, .sourceMin = -101.31, .sourceMax = 101.31, .skipsLayer = false}, // 0.1297
    {.from = 9, .to = 12, .sourceMin = -128.88, .sourceMax = 128.88, .skipsLayer = false}, // 0.0490
    {.from = 9, .to = 16, .sourceMin = -128.88, .sourceMax = -70.15, .skipsLayer = false}, // 0.4250
    {.from = 9, .to = 17, .sourceMin = -128.88, .sourceMax = -77.47, .skipsLayer = true}, // 0.0101
    {.from = 9, .to = 18, .sourceMin = -128.88, .sourceMax = -87.08, .skipsLayer = true}, // 0.0001
    {.from = 9, .to = 22, .sourceMin = 70.15, .sourceMax = 128.88, .skipsLayer = false}, // 0.4250
    {.from = 9, .to = 23, .sourceMin = 77.47, .sourceMax = 128.88, .skipsLayer = true}, // 0.0101
    {.from = 9, .to = 24, .sourceMin = 87.08, .sourceMax = 128.88, .skipsLayer = true}, // 0.0001
    {.from = 10, .to = 11, .sourceMin = -101.25, .sourceMax = 101.25, .skipsLayer = false}, // 0.1171
    {.from = 10, .to = 12, .sourceMin = -128.88, .sourceMax = 128.88, .skipsLayer = false}, // 0.3932
    {.from = 10, .to = 16, .sourceMin = -128.88, .sourceMax = -91.29, .skipsLayer = false}, // 0.4406
    {.from = 10, .to = 17, .sourceMin = -128.88, .sourceMax = -110.80, .skipsLayer = true}, // 0.0009
    {.from = 10, .to = 22, .sourceMin = 91.29, .sourceMax = 128.88, .skipsLayer = false}, // 0.4406
    {.from = 10, .to = 23, .sourceMin = 110.80, .sourceMax = 128.88, .skipsLayer = true}, // 0.0009
    {.from = 11, .to = 12, .sourceMin = -128.88, .sourceMax = 128.88, .skipsLayer = false}, // 1.0000
    {.from = 13, .to = 9, .sourceMin = 23.87, .sourceMax = 35.17, .skipsLayer = false}, // 0.0952
    {.from = 13, .to = 12, .sourceMin = 19.67, .sourceMax = 30.44, .skipsLayer = true}, // 0.0001
    {.from = 13, .to = 14, .sourceMin = 9.50, .sourceMax = 33.31, .skipsLayer = false}, // 0.7973
    {.from = 13, .to = 15, .sourceMin = 9.50, .sourceMax = 30.66, .skipsLayer = true}, // 0.0066
    {.from = 13, .to = 16, .sourceMin = 9.56, .sourceMax = 35.14, .skipsLayer = false}, // 0.0969
    {.from = 13, .to = 17, .sourceMin = 9.53, .sourceMax = 34.13, .skipsLayer = true}, // 0.0039
    {.from = 13, .to = 18, .sourceMin = 9.96, .sourceMax = 32.62, .skipsLayer = true}, // 0.0001
    {.from = 14, .to = 9, .sourceMin = 33.24, .sourceMax = 35.16, .skipsLayer = false}, // 0.0000
    {.from = 14, .to = 12, .sourceMin = 23.52, .sourceMax = 35.04, .skipsLayer = true}, // 0.0000
    {.from = 14, .to = 15, .sourceMin = 9.50, .sourceMax = 33.83, .skipsLayer = false}, // 0.8259
    {.from = 14, .to = 16, .sourceMin = 10.57, .sourceMax = 35.17, .skipsLayer = false}, // 0.1645
    {.from = 14, .to = 17, .sourceMin = 9.57, .sourceMax = 35.16, .skipsLayer = true}, // 0.0088
    {.from = 14, .to = 18, .sourceMin = 9.55, .sourceMax = 34.99, .skipsLayer = true}, // 0.0008
    {.from = 15, .to = 12, .sourceMin = 29.04, .sourceMax = 34.16, .skipsLayer = true}, // 0.0000
    {.from = 15, .to = 16, .sourceMin = 12.06, .sourceMax = 35.17, .skipsLayer = false}, // 0.6962
    {.from = 15, .to = 17, .sourceMin = 9.90, .sourceMax = 35.16, .skipsLayer = false}, // 0.1845
    {.from = 15, .to = 18, .sourceMin = 9.51, .sourceMax = 35.05, .skipsLayer = false}, // 0.1193
    {.from = 16, .to = 12, .sourceMin = 36.91, .sourceMax = 71.13, .skipsLayer = false}, // 0.0563
    {.from = 16, .to = 17, .sourceMin = 15.18, .sourceMax = 68.75, .skipsLayer = false}, // 0.9588
    {.from = 16, .to = 18, .sourceMin = 15.18, .sourceMax = 63.29, .skipsLayer = true}, // 0.0412
    {.from = 17, .to = 12, .sourceMin = 45.30, .sourceMax = 71.14, .skipsLayer = false}, // 0.0630
    {.from = 17, .to = 18, .sourceMin = 15.18, .sourceMax = 68.29, .skipsLayer = false}, // 1.0000
    {.from = 18, .to = 12, .sourceMin = 56.58, .sourceMax = 71.14, .skipsLayer = false}, // 1.0000
    {.from = 19, .to = 9, .sourceMin = 23.87, .sourceMax = 35.17, .skipsLayer = false}, // 0.0952
    {.from = 19, .to = 12, .sourceMin = 19.67, .sourceMax = 30.44, .skipsLayer = true}, // 0.0001
    {.from = 19, .to = 20, .sourceMin = 9.50, .sourceMax = 33.31, .skipsLayer = false}, // 0.7973
    {.from = 19, .to = 21, .sourceMin = 9.50, .sourceMax = 30.66, .skipsLayer = true}, // 0.0066
    {.from = 19, .to = 22, .sourceMin = 9.56, .sourceMax = 35.14, .skipsLayer = false}, // 0.0969
    {.from = 19, .to = 23, .sourceMin = 9.53, .sourceMax = 34.13, .skipsLayer = true}, // 0.0039
    {.from = 19, .to = 24, .sourceMin = 9.96, .sourceMax = 32.62, .skipsLayer = true}, // 0.0001
    {.from = 20, .to = 9, .sourceMin = 33.24, .sourceMax = 35.16, .skipsLayer = false}, // 0.0000
    {.from = 20, .to = 12, .sourceMin = 23.52, .sourceMax = 35.04, .skipsLayer = true}, // 0.0000
    {.from = 20, .to = 21, .sourceMin = 9.50, .sourceMax = 33.83, .skipsLayer = false}, // 0.8259
    {.from = 20, .to = 22, .sourceMin = 10.57, .sourceMax = 35.17, .skipsLayer = false}, // 0.1645
    {.from = 20, .to = 23, .sourceMin = 9.57, .sourceMax = 35.16, .skipsLayer = true}, // 0.0088
    {.from = 20, .to = 24, .sourceMin = 9.55, .sourceMax = 34.99, .skipsLayer = true}, // 0.0008
    {.from = 21, .to = 12, .sourceMin = 29.04, .sourceMax = 34.16, .skipsLayer = true}, // 0.0000
    {.from = 21, .to = 22, .sourceMin = 12.06, .sourceMax = 35.17, .skipsLayer = false}, // 0.6962
    {.from = 21, .to = 23, .sourceMin = 9.90, .sourceMax = 35.16, .skipsLayer = false}, // 0.1845
    {.from = 21, .to = 24, .sourceMin = 9.51, .sourceMax = 35.05, .skipsLayer = false}, // 0.1193
    {.from = 22, .to = 23, .sourceMin = 15.18, .sourceMax = 68.75, .skipsLayer = false}, // 0.9588
    {.from = 22, .to = 24, .sourceMin = 15.18, .sourceMax = 63.29, .skipsLayer = true}, // 0.0412
    {.from = 23, .to = 24, .sourceMin = 15.18, .sourceMax = 68.29, .skipsLayer = false}, // 1.0000
  }};
  static constexpr std::array<PassiveSpec, 11> passives{{
    {.kind = PassiveSpec::Kind::kBarrel, .position = 0.4899, .bandLow = 0.0000, .bandHigh = 0.4999, .low = -25.00, .high = 25.00, .x2X0 = 2.224e-03},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 0.8002, .bandLow = 0.5096, .bandHigh = 1.1999, .low = -25.00, .high = 25.00, .x2X0 = 1.228e-03},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 1.7708, .bandLow = 1.2096, .bandHigh = 2.5000, .low = -25.00, .high = 25.00, .x2X0 = 9.985e-04},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 3.7955, .bandLow = 2.5096, .bandHigh = 6.8657, .low = -63.98, .high = 63.98, .x2X0 = 1.132e-02},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 7.9710, .bandLow = 7.8977, .bandHigh = 8.8055, .low = -63.98, .high = 63.98, .x2X0 = 3.692e-04},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 10.0133, .bandLow = 9.8542, .bandHigh = 11.7623, .low = -63.98, .high = 63.98, .x2X0 = 4.094e-04},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 13.8492, .bandLow = 12.7722, .bandHigh = 18.7300, .low = -63.98, .high = 63.98, .x2X0 = 4.863e-04},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 21.2092, .bandLow = 20.1379, .bandHigh = 29.9999, .low = -63.98, .high = 63.98, .x2X0 = 2.027e-02},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 38.2067, .bandLow = 31.0792, .bandHigh = 44.8000, .low = -128.38, .high = 128.38, .x2X0 = 2.557e-02},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 52.1226, .bandLow = 45.3527, .bandHigh = 59.8000, .low = -128.38, .high = 128.38, .x2X0 = 4.899e-04},
    {.kind = PassiveSpec::Kind::kBarrel, .position = 69.6851, .bandLow = 60.2898, .bandHigh = 79.8000, .low = -128.38, .high = 128.38, .x2X0 = 6.376e-04},
  }};
  // clang-format on
  static constexpr int nSensors = TRKFT3::nSensors + 2352 + 56544;
  static constexpr int kFirstDiscSensor = TRKFT3::kFirstDiscSensor;
  /// The first TOF sensor's ID: O2's IOTOF chip 0, after every TRKFT3 chip.
  static constexpr int kFirstTOFSensor = TRKFT3::nSensors;
};

} // namespace o2::gtrk::detectors

#endif // GTRK_DETECTORS_TRKFT3TOF_H_
