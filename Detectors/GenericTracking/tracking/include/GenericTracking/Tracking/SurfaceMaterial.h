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

/// \file SurfaceMaterial.h
/// \brief A step's material from two table reads instead of a walk through
/// the LUT.
///
/// For every linked pair of layers (Det::layerLinks), the LUT's material in
/// two halves, tabulated once per run (SurfaceMaterialFill.h): from the pair's
/// first layer (the link's `from`) to the midpoint between its crossings of
/// the two layers, over the first layer's surface; from the midpoint to the
/// second layer, over the second's (azimuth and z on a cylinder, azimuth and r
/// on a disc). A step between the two layers reads each half where it crosses
/// that half's layer, scaled from the tabulated ray's angle to the layer to
/// its own, so that a track that bends or does not come from the vertex reads
/// the material about each of its ends. Steps between layers that are not
/// linked are not tabulated (`SurfaceStep::found`); they walk the LUT.

#ifndef GTRK_TRACKING_SURFACEMATERIAL_H_
#define GTRK_TRACKING_SURFACEMATERIAL_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"

#include "GPUCommonMath.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace o2::gtrk
{

/// One layer's surface as a pair's half is tabulated over it.
struct MaterialSurface {
  Scalar secondMin = 0.F; ///< z on a cylinder, r on a disc [cm]
  Scalar secondBinsPerCm = 0.F;
  std::uint8_t barrel = 1U;
};

struct MaterialCell {
  float x2X0 = 0.F; ///< x/X0 along the tabulated ray
  float xRho = 0.F; ///< rho * length along the tabulated ray [g/cm^2]
};

/// A step's material by halves about its midpoint; `found` false where the
/// pair is not tabulated.
struct SurfaceStep {
  Scalar nearerX2X0 = 0.F;
  Scalar fartherX2X0 = 0.F;
  Scalar xRho = 0.F;
  bool found = false;
};

template <DetectorLayout Det>
struct SurfaceMaterialMap {
  static constexpr int kLayers = Detector<Det>::nLayers;
  static constexpr int kPairs = static_cast<int>(std::size(Det::layerLinks));
  static constexpr int kPhiBins = 256;
  static constexpr int kSecondBins = 128;
  static constexpr std::size_t kCells = static_cast<std::size_t>(kPhiBins) * kSecondBins;

  /// Per ordered layer pair: 1 + the pair's index, negative where the step
  /// goes from the pair's second layer to its first; 0 if not linked.
  std::array<std::int16_t, static_cast<std::size_t>(kLayers) * kLayers> pairOf{};
  /// Per pair: the first layer's surface [0] and the second's [1].
  std::array<std::array<MaterialSurface, 2>, kPairs> surfaces{};
  /// Per pair and half (0: about the first layer, 1: about the second), the
  /// cells.
  std::array<MaterialCell, static_cast<std::size_t>(kPairs) * 2U * kCells> cells{};

  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::size_t index(int pair, int half, Scalar phi, Scalar second) const
  {
    constexpr Scalar kTwoPi = 6.283185307F;
    const MaterialSurface& s = surfaces[pair][half];
    int iPhi = static_cast<int>((phi < 0.F ? phi + kTwoPi : phi) * (kPhiBins / kTwoPi));
    iPhi = iPhi < 0 ? 0 : (iPhi >= kPhiBins ? kPhiBins - 1 : iPhi);
    int iSecond = static_cast<int>((second - s.secondMin) * s.secondBinsPerCm);
    iSecond = iSecond < 0 ? 0 : (iSecond >= kSecondBins ? kSecondBins - 1 : iSecond);
    return (((static_cast<std::size_t>(pair) * 2U) + static_cast<std::size_t>(half)) * kCells) + (static_cast<std::size_t>(iPhi) * kSecondBins) + static_cast<std::size_t>(iSecond);
  }

  /// The step from \a start on layer \a fromLayer to \a end on \a toLayer.
  [[nodiscard]] GPUhdi() SurfaceStep step(LayerID fromLayer, LayerID toLayer, const Scalar start[3], const Scalar end[3]) const
  {
    using Math = o2::gpu::GPUCommonMath;
    SurfaceStep out;
    if (fromLayer >= kLayers || toLayer >= kLayers) {
      return out;
    }
    const int code = pairOf[(static_cast<std::size_t>(fromLayer) * kLayers) + toLayer];
    if (code == 0) {
      return out;
    }
    out.found = true;
    const bool reversed = code < 0;
    const int pair = (reversed ? -code : code) - 1;
    const Scalar dx = end[0] - start[0];
    const Scalar dy = end[1] - start[1];
    const Scalar dz = end[2] - start[2];
    const Scalar length = Math::Sqrt((dx * dx) + (dy * dy) + (dz * dz));
    if (!(length > 1e-5F)) {
      return out;
    }
    const Scalar d[3] = {dx / length, dy / length, dz / length};
    // A half read at \a p on its surface, scaled from the tabulated ray's angle
    // to the surface to the step's.
    const auto read = [&](int half, const Scalar p[3]) {
      const MaterialSurface& s = surfaces[pair][half];
      const Scalar r = Math::Sqrt((p[0] * p[0]) + (p[1] * p[1]));
      const Scalar distance = Math::Sqrt((r * r) + (p[2] * p[2]));
      Scalar rayCosine = 1.F;
      Scalar stepCosine = 1.F;
      if (s.barrel != 0U) {
        rayCosine = distance > 0.F ? r / distance : 1.F;
        stepCosine = r > 0.F ? Math::Abs(((d[0] * p[0]) + (d[1] * p[1])) / r) : 1.F;
      } else {
        rayCosine = distance > 0.F ? Math::Abs(p[2]) / distance : 1.F;
        stepCosine = Math::Abs(d[2]);
      }
      Scalar scale = stepCosine > 0.F ? rayCosine / stepCosine : 5.F;
      scale = scale < 0.2F ? 0.2F : (scale > 5.F ? 5.F : scale);
      const MaterialCell& c = cells[index(pair, half, Math::ATan2(p[1], p[0]), s.barrel != 0U ? p[2] : r)];
      return MaterialCell{.x2X0 = scale * c.x2X0, .xRho = scale * c.xRho};
    };
    const MaterialCell first = read(0, reversed ? end : start);
    const MaterialCell second = read(1, reversed ? start : end);
    out.nearerX2X0 = reversed ? second.x2X0 : first.x2X0;
    out.fartherX2X0 = reversed ? first.x2X0 : second.x2X0;
    out.xRho = first.xRho + second.xRho;
    return out;
  }
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_SURFACEMATERIAL_H_
