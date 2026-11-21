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

/// \file BetweenBoundCheck.h
/// \brief The kink between a link's layers against the link's bound.
///
/// The tracklet stage tests a pair at its link's bounds -- on the kink of the
/// link's gap (`IterationParameters::betweenBound`) and on the source term's
/// rescaling to the pair's span (`betweenSourceScale`) -- before computing
/// the pair's own (`betweenScattering`), and refuses it if that test fails.
/// That refuses only pairs the own test would refuse too as long as no
/// pair's kink or rescaling exceeds its bound, which this checks.

#ifndef GTRK_TEST_BETWEENBOUNDCHECK_H_
#define GTRK_TEST_BETWEENBOUNDCHECK_H_

#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Base/LayerEnvelopes.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/LinkMaterial.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/SearchWindow.h"
#include "GenericTracking/Tracking/TripletFit.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

namespace o2::gtrk::test
{

/// One link with material between its layers: the largest of its pairs'
/// kinks over its bound, and of their source terms' rescalings.
struct BetweenBoundRatio {
  LayerID from = 0;
  LayerID to = 0;
  double worst = 0.;
  double worstSourceScale = 0.;
};

/// The value a fraction \a u of the way from \a low to \a high.
[[nodiscard]] inline double interpolate(double low, double high, double u) { return low + (u * (high - low)); }

/// Per link with material between its layers, `betweenScattering` of pairs
/// a prompt iteration makes -- a source over the source layer's extent, the
/// target where the line from a vertex in the beam spot's reach through the
/// source meets the target's nominal surface, both at an azimuth drawn
/// uniformly -- over the link's `betweenBound` and `betweenSourceScale`: at
/// most one where the bounds hold.
template <DetectorLayout Det>
std::vector<BetweenBoundRatio> betweenOverBound(const LinkMaterial<Det>& material, const SearchWindowParams& windows, const VertexConstraint& beamSpot)
{
  constexpr int kDraws = 8192;
  IterationSettings settings;
  settings.windows = windows;
  const auto parameters = iterationParameters<Det>(settings, beamSpot, material);
  const double reach = windows.nSigma * beamSpot.sigmaZ;
  std::mt19937 rng(20261010);
  std::uniform_real_distribution<double> uniform(0., 1.);
  std::vector<BetweenBoundRatio> ratios;
  for (std::size_t link = 0; link < Det::layerLinks.size(); ++link) {
    if (!(parameters.betweenBound[link] > 0.F)) {
      continue;
    }
    const auto& pair = Det::layerLinks[link];
    BetweenBoundRatio ratio{.from = pair.from, .to = pair.to};
    const bool barrelFrom = Detector<Det>::isBarrelLayer(pair.from);
    const bool barrelTo = Detector<Det>::isBarrelLayer(pair.to);
    const double positionTo = PlacedEnvelopes<Det>::get(pair.to).position;
    for (int draw = 0; draw < kDraws; ++draw) {
      auto vertex = vertexSeed(beamSpot);
      vertex.z = static_cast<Scalar>(interpolate(beamSpot.z - reach, beamSpot.z + reach, uniform(rng)));
      const double phi = kTwoPi<double> * uniform(rng);
      const double r = barrelFrom ? PlacedEnvelopes<Det>::get(pair.from).position : interpolate(radiusMin<Det>(pair.from), radiusMax<Det>(pair.from), uniform(rng));
      const double z = barrelFrom ? interpolate(std::max(pair.sourceMin, zMin<Det>(pair.from)), std::min(pair.sourceMax, zMax<Det>(pair.from)), uniform(rng)) : PlacedEnvelopes<Det>::get(pair.from).position;
      // The line from the vertex through the source to the target's surface.
      const double scale = barrelTo ? positionTo / r : (positionTo - vertex.z) / (z - vertex.z);
      const double rTarget = barrelTo ? positionTo : r * scale;
      const double zTarget = barrelTo ? vertex.z + ((z - vertex.z) * scale) : positionTo;
      if (!(scale > 1.) || rTarget < radiusMin<Det>(pair.to) || rTarget > radiusMax<Det>(pair.to) || zTarget < zMin<Det>(pair.to) || zTarget > zMax<Det>(pair.to)) {
        continue;
      }
      const auto point = [&](double pr, double pz) { return PairPoint{.x = static_cast<Scalar>(pr * std::cos(phi)), .y = static_cast<Scalar>(pr * std::sin(phi)), .z = static_cast<Scalar>(pz), .r = static_cast<Scalar>(pr)}; };
      const auto source = point(r, z);
      const auto target = point(rTarget, zTarget);
      Scalar spanX = 0.F;
      static_cast<void>(sourceScattering<Det>(material, parameters.pairCut, vertex, pair.from, source, atan2Of(source.y, source.x), spanX));
      const auto own = betweenScattering<Det>(material, parameters.pairCut, betweenSource<Det>(source, pair.from), link, pair.to, vertex, source, target, spanX);
      ratio.worst = std::max(ratio.worst, static_cast<double>(own.between) / static_cast<double>(parameters.betweenBound[link]));
      ratio.worstSourceScale = std::max(ratio.worstSourceScale, static_cast<double>(own.sourceScale) / static_cast<double>(parameters.betweenSourceScale[link]));
    }
    ratios.push_back(ratio);
  }
  return ratios;
}

} // namespace o2::gtrk::test

#endif // GTRK_TEST_BETWEENBOUNDCHECK_H_
