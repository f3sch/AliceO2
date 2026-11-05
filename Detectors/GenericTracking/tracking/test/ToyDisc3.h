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

/// \file ToyDisc3.h
/// \brief Three disc layers and no barrel at all.
///
/// Its purpose is the empty barrel array.  Keeping an all-disc detector
/// instantiated from the start means `if constexpr (Detector<Det>::hasBarrel)` has a
/// tested false branch long before the forward tracking of stage 6 arrives,
/// rather than that branch being compiled for the first time months later.

#ifndef GTRK_DETECTORS_TOYDISC3_H_
#define GTRK_DETECTORS_TOYDISC3_H_

#include "AnalyticLayout.h"

#include <array>

namespace o2::gtrk::detectors
{

struct ToyDisc3 {
  static constexpr const char* name = "ToyDisc3";

  static constexpr SensorSpec sensor{.shape = ShapeType::kRectangle, .halfX = 1.0, .halfY = 2.0, .halfThickness = 0.005, .pitchX = 0.0025, .pitchY = 0.0025};

  static constexpr RingSpec ring{.r = 5., .nModules = 6, .phi0 = 0., .zOffset = 0., .gap = 0.06};

  static constexpr std::array<BarrelLayerSpec, 0> barrelSpecs{};

  static constexpr std::array<DiscLayerSpec, 3> discSpecs{{
    {.z = 10., .nRings = 1, .rings = {ring}, .sensor = sensor},
    {.z = 20., .nRings = 1, .rings = {ring}, .sensor = sensor},
    {.z = 30., .nRings = 1, .rings = {ring}, .sensor = sensor},
  }};

  // Traced by o2-gtrk-layer-links' former toy mode for ToyDisc3, now kept by hand:
  // pT >= 0.1 GeV at Bz = 20 kG, beam spot z in [-1, 1] cm, |eta| <= 6, source intervals padded by 0.5 cm.
  // The comment on each link is the fraction of its source layer's crossings that take it.
  static constexpr std::array<LayerLink, 2> layerLinks{{
    {.from = 0, .to = 1, .sourceMin = 2.50, .sourceMax = 4.23, .skipsLayer = false}, // 1.0000
    {.from = 1, .to = 2, .sourceMin = 2.50, .sourceMax = 5.31, .skipsLayer = false}, // 1.0000
  }};

  static constexpr auto barrelLayers = placedLayers(barrelSpecs);
  static constexpr auto discLayers = placedLayers(discSpecs);
  static constexpr int nSensors = countSensors(barrelLayers, discLayers);
};

inline const bool kToyDisc3Published = publishToy<ToyDisc3>();

} // namespace o2::gtrk::detectors

#endif // GTRK_DETECTORS_TOYDISC3_H_
