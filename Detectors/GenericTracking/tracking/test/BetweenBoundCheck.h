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
/// The tracklet stage tests a pair at its link's bound on the kink of a
/// passive scatterer between the layers (`IterationParameters::betweenBound`)
/// before computing the pair's own (`betweenScattering`), and refuses it if
/// that test fails. That refuses only pairs the own kink would refuse too as
/// long as no pair's kink exceeds the bound, which this checks.

#ifndef GTRK_TEST_BETWEENBOUNDCHECK_H_
#define GTRK_TEST_BETWEENBOUNDCHECK_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/SearchWindow.h"
#include "GenericTracking/Tracking/TripletFit.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

namespace o2::gtrk::test
{

/// One link with a passive between its layers: the largest of its pairs'
/// kinks over its bound.
struct BetweenBoundRatio {
  LayerID from = 0;
  LayerID to = 0;
  double worst = 0.;
  /// Pairs whose kink over the link's `betweenPassives` is not, bit for bit,
  /// the one over every passive: none.
  int maskMismatches = 0;
  /// Pairs whose kink from the source's part (`betweenSource`, the beam-line
  /// tracklet stage) is not, bit for bit, the
  /// masked one: none.
  int sourcePartMismatches = 0;
};

/// The value a fraction \a u of the way from \a low to \a high.
[[nodiscard]] inline double interpolate(double low, double high, double u) { return low + (u * (high - low)); }

/// The point of layer \a layer at azimuth 0, fractions \a u of the way
/// through its radii and \a v through its z.
template <DetectorLayout Det>
[[nodiscard]] PairPoint extentPoint(LayerID layer, double u, double v)
{
  const auto r = static_cast<Scalar>(interpolate(radiusMin<Det>(layer), radiusMax<Det>(layer), u));
  return PairPoint{.x = r, .z = static_cast<Scalar>(interpolate(zMin<Det>(layer), zMax<Det>(layer), v)), .r = r};
}

/// Per link with a passive between, `betweenScattering` at points over both
/// layers' extents -- every corner, then uniform draws -- and vertices over
/// the beam spot's reach, over the link's `betweenBound`: at most one where
/// the bound holds.
template <DetectorLayout Det>
std::vector<BetweenBoundRatio> betweenOverBound(const ScatteringTable<Det>& scattering, const SearchWindowParams& windows, const VertexConstraint& beamSpot)
{
  constexpr int kCorners = 32;
  constexpr int kDraws = 4096;
  IterationSettings settings;
  settings.windows = windows;
  const auto parameters = iterationParameters<Det>(settings, beamSpot, scattering);
  const double reach = windows.nSigma * beamSpot.sigmaZ;
  std::mt19937 rng(20261010);
  std::uniform_real_distribution<double> uniform(0., 1.);
  std::vector<BetweenBoundRatio> ratios;
  for (std::size_t link = 0; link < Det::layerLinks.size(); ++link) {
    if (!parameters.passiveBetween[link]) {
      continue;
    }
    const auto& pair = Det::layerLinks[link];
    BetweenBoundRatio ratio{.from = pair.from, .to = pair.to};
    for (int draw = 0; draw < kCorners + kDraws; ++draw) {
      std::array<double, 5> u{};
      for (std::size_t i = 0; i < u.size(); ++i) {
        u[i] = draw < kCorners ? static_cast<double>((draw >> i) & 1) : uniform(rng);
      }
      auto vertex = vertexSeed(beamSpot);
      vertex.z = static_cast<Scalar>(interpolate(beamSpot.z - reach, beamSpot.z + reach, u[4]));
      const auto source = extentPoint<Det>(pair.from, u[0], u[1]);
      const auto target = extentPoint<Det>(pair.to, u[2], u[3]);
      const Scalar kink = betweenScattering<Det>(parameters.scattering, parameters.pairCut, vertex, source, target);
      const Scalar masked = betweenScattering<Det>(parameters.scattering, parameters.pairCut, vertex, source, target, parameters.betweenPassives[link]);
      ratio.maskMismatches += sameBits(masked, kink) ? 0 : 1;
      const Scalar fromSource = betweenScattering<Det>(betweenSource<Det>(parameters.scattering, vertex, source, parameters.pairCut.logBeta2), parameters.pairCut, vertex, target, parameters.betweenPassives[link]);
      ratio.sourcePartMismatches += sameBits(fromSource, masked) ? 0 : 1;
      ratio.worst = std::max(ratio.worst, static_cast<double>(kink) / static_cast<double>(parameters.betweenBound[link]));
    }
    ratios.push_back(ratio);
  }
  return ratios;
}

/// Per layout link, segments between points over both layers' extents --
/// every corner, then uniform draws -- whose passive visits (`visitPassivesBetween`)
/// differ between every passive and the table's `outOf` the source layer or
/// `into` the target: none.
template <DetectorLayout Det>
int segmentMaskMismatches(const ScatteringTable<Det>& table)
{
  constexpr int kCorners = 16;
  constexpr int kDraws = 4096;
  std::mt19937 rng(20261015);
  std::uniform_real_distribution<double> uniform(0., 1.);
  const PairCut cut{.invBetaMomentum2 = 1.F};
  int mismatches = 0;
  for (const auto& link : Det::layerLinks) {
    for (int draw = 0; draw < kCorners + kDraws; ++draw) {
      std::array<double, 4> u{};
      for (std::size_t i = 0; i < u.size(); ++i) {
        u[i] = draw < kCorners ? static_cast<double>((draw >> i) & 1) : uniform(rng);
      }
      const auto from = extentPoint<Det>(link.from, u[0], u[1]);
      const auto to = extentPoint<Det>(link.to, u[2], u[3]);
      const auto visits = [&](std::uint64_t passives) {
        std::vector<std::uint32_t> seen;
        visitPassivesBetween(table, cut, from, to, 1.F, [&](Scalar variance, Scalar t) { seen.push_back(std::bit_cast<std::uint32_t>(variance * t)); }, passives);
        return seen;
      };
      const auto all = visits(ScatteringTable<Det>::kAllPassives);
      mismatches += visits(table.outOf(link.from)) != all ? 1 : 0;
      mismatches += visits(table.into(link.to)) != all ? 1 : 0;
    }
  }
  return mismatches;
}

} // namespace o2::gtrk::test

#endif // GTRK_TEST_BETWEENBOUNDCHECK_H_
