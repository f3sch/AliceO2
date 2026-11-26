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

/// \file testPassiveScattering.cxx
/// \brief A passive scatterer among the layers -- the beam pipe with layers
///        on both sides of it, as ALICE 3's -- in the seeding's cuts and
///        windows, as the link material holds it: in the gap of the link
///        across it and in the vertex gaps beyond it.
///
/// The layout is a test fixture: three barrel layers inside a pipe and three
/// outside it. Tracks are followed as helices from the beam line and kinked
/// at every layer and at the pipe with the scattering the cut models --
/// Highland's theta0 at normal incidence, the path factor from the track's
/// dip -- so that the cut's true pairs can be counted against its n sigma.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Passive Scattering
#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/Statistics.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/LinkMaterial.h"
#include "GenericTracking/Tracking/NeighbourCut.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/SearchWindow.h"
#include "GenericTracking/Tracking/TripletCut.h"

#include "AnalyticLayout.h"
#include "BandMaterial.h"
#include "BetweenBoundCheck.h"
#include "Oracles.h"

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

namespace o2::gtrk
{
namespace
{

struct ToyPipe {
  static constexpr const char* name = "ToyPipe";

  static constexpr SensorSpec sensor{.shape = ShapeType::kRectangle, .halfX = 0.5, .halfY = 10., .halfThickness = 0.005, .pitchX = 0.0025, .pitchY = 0.0025};

  // ALICE 3's shape: three layers inside the pipe, three outside it.
  static constexpr std::array<BarrelLayerSpec, 6> barrelSpecs{{
    {.r = 0.5, .nStaves = 4, .nModulesPerStave = 1, .phi0 = 0., .tilt = 0., .moduleGap = 0., .sensor = sensor},
    {.r = 1.2, .nStaves = 8, .nModulesPerStave = 1, .phi0 = 0., .tilt = 0., .moduleGap = 0., .sensor = sensor},
    {.r = 2.5, .nStaves = 16, .nModulesPerStave = 1, .phi0 = 0., .tilt = 0., .moduleGap = 0., .sensor = sensor},
    {.r = 7., .nStaves = 45, .nModulesPerStave = 1, .phi0 = 0., .tilt = 0., .moduleGap = 0., .sensor = sensor},
    {.r = 10., .nStaves = 64, .nModulesPerStave = 1, .phi0 = 0., .tilt = 0., .moduleGap = 0., .sensor = sensor},
    {.r = 13., .nStaves = 83, .nModulesPerStave = 1, .phi0 = 0., .tilt = 0., .moduleGap = 0., .sensor = sensor},
  }};

  static constexpr std::array<DiscLayerSpec, 0> discSpecs{};
  static constexpr auto barrelLayers = placedLayers(barrelSpecs);
  static constexpr auto discLayers = placedLayers(discSpecs);

  static constexpr std::array<LayerLink, 5> layerLinks{{
    {.from = 0, .to = 1, .sourceMin = -10., .sourceMax = 10.},
    {.from = 1, .to = 2, .sourceMin = -10., .sourceMax = 10.},
    {.from = 2, .to = 3, .sourceMin = -10., .sourceMax = 10.},
    {.from = 3, .to = 4, .sourceMin = -10., .sourceMax = 10.},
    {.from = 4, .to = 5, .sourceMin = -10., .sourceMax = 10.},
  }};

  static constexpr std::int32_t nSensors = countSensors(barrelLayers, discLayers);
};

[[maybe_unused]] const bool kToyPipePublished = publishToy<ToyPipe>();

using Det = ToyPipe;

constexpr double kLayerX2X0 = 1.e-3;
/// The pipe between layers 2 and 3, 4 cm short of layer 3.
constexpr double kPipeRadius = 3.;
constexpr double kBzKiloGauss = 20.;
constexpr double kMinPtGeV = 0.3;
constexpr std::uint32_t kCrossingLink = 2; // 2 -> 3
constexpr std::uint32_t kBeyondLink = 3;   // 3 -> 4, the pipe inside its source

SearchWindowParams windowParams()
{
  SearchWindowParams params;
  params.minPtGeV = kMinPtGeV;
  params.bzKiloGauss = kBzKiloGauss;
  return params;
}

const VertexConstraint kBeamSpot{.z = 0., .sigmaZ = 1., .sigmaTransverse = 0.002};

std::unique_ptr<LinkMaterial<Det>> table(double pipeX2X0)
{
  const std::array<ThinCylinder, 1> pipe{{{.r = kPipeRadius, .x2X0 = pipeX2X0}}};
  return bandMaterial<Det>(kLayerX2X0, pipe);
}

/// The pair's source term and its kink between, as the tracklet stage takes
/// them: the source's at the pair's span (`PairScattering::sourceScale`);
/// \a withBetween false leaves the kink between out.
struct PairTerms {
  Scalar source = 0.F;
  Scalar between = 0.F;
};

PairTerms pairTerms(const IterationParameters<Det>& parameters, std::uint32_t link, const VertexSeed& vertex, const PairPoint& source, const PairPoint& target, bool withBetween)
{
  const auto& pair = Det::layerLinks[link];
  Scalar spanX = 0.F;
  PairTerms terms{.source = sourceScattering<Det>(*parameters.material, parameters.pairCut, vertex, pair.from, source, atan2Of(source.y, source.x), spanX)};
  if (withBetween && parameters.betweenBound[link] > 0.F) {
    const auto between = betweenScattering<Det>(*parameters.material, parameters.pairCut, betweenSource<Det>(source, pair.from), link, pair.to, vertex, source, target, spanX);
    terms.source *= between.sourceScale;
    terms.between = between.between;
  }
  return terms;
}

/// One track's crossings of the six layers, kinked on the way.
struct Crossings {
  std::array<PairPoint, 6> points{};
  bool complete = false;
};

/// Follow a helix of transverse momentum \a pt from the beam line, kinking
/// it at every radius in \a scatterers with its Highland theta0 at the
/// crossing, and record where it crosses each layer.
Crossings follow(double pt, double phi, double cotTheta, double charge, double vertexZ, std::span<const double> scatterers, std::span<const double> scattererX2X0, std::mt19937& rng)
{
  std::normal_distribution<double> gauss(0., 1.);
  const double radius = curvatureRadius(pt, kBzKiloGauss);
  const auto& layers = Det::barrelSpecs;
  // Every radius the track has to reach, in order, and what happens there.
  struct Stop {
    double r;
    std::int32_t layer; // -1: a scatterer only
    double x2X0; // scattering at the crossing
  };
  std::vector<Stop> stops;
  stops.reserve(layers.size() + scatterers.size());
  std::int32_t layer = 0;
  for (const auto& spec : layers) {
    stops.push_back({.r = spec.r, .layer = layer++, .x2X0 = kLayerX2X0});
  }
  for (std::size_t i = 0; i < scatterers.size(); ++i) {
    stops.push_back({.r = scatterers[i], .layer = -1, .x2X0 = scattererX2X0[i]});
  }
  std::ranges::sort(stops, {}, &Stop::r);

  Crossings result;
  double x = 0.;
  double y = 0.;
  double z = vertexZ;
  for (const auto& stop : stops) {
    // The helix from (x, y) at azimuth phi, by transverse path s:
    // x + (R/q)(sin(phi + q s/R) - sin phi), y - (R/q)(cos(phi + q s/R) - cos phi).
    const auto at = [&](double path, double& px, double& py) {
      const double turned = phi + (charge * path / radius);
      px = x + (radius / charge * (std::sin(turned) - std::sin(phi)));
      py = y - (radius / charge * (std::cos(turned) - std::cos(phi)));
    };
    double low = 0.;
    double high = 0.;
    double px = x;
    double py = y;
    // Grow the bracket in quarter-radius steps; turning back before the stop ends the track.
    double previousR = std::hypot(x, y);
    while (true) {
      high += 0.25 * std::min(radius, stop.r);
      at(high, px, py);
      const double r = std::hypot(px, py);
      if (r >= stop.r) {
        break;
      }
      if (r <= previousR) {
        return result;
      }
      previousR = r;
      low = high;
    }
    for (std::int32_t i = 0; i < 60; ++i) {
      const double mid = 0.5 * (low + high);
      at(mid, px, py);
      (std::hypot(px, py) < stop.r ? low : high) = mid;
    }
    const double path = 0.5 * (low + high);
    at(path, px, py);
    x = px;
    y = py;
    z += cotTheta * path;
    phi += charge * path / radius;
    if (stop.layer >= 0) {
      const auto r = static_cast<Scalar>(std::hypot(x, y));
      result.points[static_cast<std::uint32_t>(stop.layer)] = {
        .x = static_cast<Scalar>(x),
        .y = static_cast<Scalar>(y),
        .z = static_cast<Scalar>(z),
        .r = r,
        .varZ = 2.5e-7F,
        .varR = 2.5e-7F,
        .varRPhi = 2.5e-7F,
      };
    }
    // The kink: Highland at normal incidence times the path factor, at the
    // track's momentum, in two projections.
    const double sinTheta = 1. / std::sqrt(1. + (cotTheta * cotTheta));
    const double pathFactor = 1. / sinTheta;
    const double p = pt / sinTheta;
    const double theta0 = std::sqrt(detail::highlandFactor(stop.x2X0, 1.) * pathFactor) / p;
    phi += gauss(rng) * theta0 / sinTheta;
    const double theta = std::atan2(1., cotTheta) + (gauss(rng) * theta0);
    cotTheta = std::cos(theta) / std::sin(theta);
  }
  result.complete = true;
  return result;
}

struct LinkCount {
  std::int32_t pairs = 0;
  std::int32_t accepted = 0;
  std::int32_t inPhiWindow = 0;
  std::int32_t inSecondWindow = 0;
};

/// Pairs of consecutive layers of tracks at pT_min -- where the curvature test
/// is the n-sigma test it states -- tested with the cut
/// as the tracklet stage runs it -- \a withBetween false leaves the kink
/// between source and target out, \a cutTable is the table the cut is built
/// with.
std::array<LinkCount, 5> countPairs(const LinkMaterial<Det>& cutTable, bool withBetween, double pipeX2X0, std::int32_t tracks)
{
  IterationSettings settings;
  settings.windows = windowParams();
  const auto parameters = iterationParameters<Det>(settings, kBeamSpot, cutTable);
  const auto beamLine = vertexSeed(kBeamSpot);
  std::mt19937 rng(20260930);
  std::uniform_real_distribution<double> uniform(0., 1.);
  std::normal_distribution<double> gauss(0., 1.);
  const std::array<double, 1> scatterers{kPipeRadius};
  const std::array<double, 1> scattererX2X0{pipeX2X0};
  std::array<LinkCount, 5> counts{};
  for (std::int32_t track = 0; track < tracks; ++track) {
    const double phi = kTwoPi<double> * uniform(rng);
    const double cotTheta = -0.5 + uniform(rng);
    const double charge = track % 2 == 0 ? 1. : -1.;
    const double vertexZ = gauss(rng) * kBeamSpot.sigmaZ;
    const auto crossings = follow(kMinPtGeV, phi, cotTheta, charge, vertexZ, scatterers, scattererX2X0, rng);
    if (!crossings.complete) {
      continue;
    }
    for (std::uint32_t link = 0; link < Det::layerLinks.size(); ++link) {
      const auto& source = crossings.points[Det::layerLinks[link].from];
      const auto& target = crossings.points[Det::layerLinks[link].to];
      const auto terms = pairTerms(parameters, link, beamLine, source, target, withBetween);
      auto& count = counts[link];
      ++count.pairs;
      if (!pairCompatible(parameters.pairCut, beamLine, source, terms.source, target, terms.between)) {
        continue;
      }
      ++count.accepted;
      const auto& window = parameters.windows[link];
      const auto dPhi = toPMPiOf<double>(std::atan2(target.y, target.x) - std::atan2(source.y, source.x));
      count.inPhiWindow += std::abs(dPhi) <= static_cast<double>(window.phiHalfWidth) ? 1 : 0;
      const double centre = (static_cast<double>(window.secondScale) * static_cast<double>(source.z)) + static_cast<double>(window.secondOffset);
      count.inSecondWindow += std::abs(static_cast<double>(target.z) - centre) <= static_cast<double>(window.secondHalfWidth) ? 1 : 0;
    }
  }
  return counts;
}

/// Consecutive triplets and neighbours of tracks at pT_min, through the
/// stages' own cuts with \a cutTable: a triplet counted when both its pairs
/// pass, a neighbour when both its triplets do -- what the stages see.
struct ChainCounts {
  std::array<LinkCount, 4> triplets{}; ///< first link 0..3
  std::array<LinkCount, 3> neighbours{};
};

ChainCounts countChains(const LinkMaterial<Det>& cutTable, double pipeX2X0, std::int32_t tracks)
{
  IterationSettings settings;
  settings.windows = windowParams();
  const auto parameters = iterationParameters<Det>(settings, kBeamSpot, cutTable);
  const auto beamLine = vertexSeed(kBeamSpot);
  std::mt19937 rng(20261001);
  std::uniform_real_distribution<double> uniform(0., 1.);
  std::normal_distribution<double> gauss(0., 1.);
  const std::array<double, 1> scatterers{kPipeRadius};
  const std::array<double, 1> scattererX2X0{pipeX2X0};
  const TripletTargetBounds bounds{.varZ = 2.5e-7F, .varR = 2.5e-7F, .varRPhi = 2.5e-7F};
  ChainCounts counts;
  for (std::int32_t track = 0; track < tracks; ++track) {
    const double phi = kTwoPi<double> * uniform(rng);
    const double cotTheta = -0.5 + uniform(rng);
    const double charge = track % 2 == 0 ? 1. : -1.;
    const double vertexZ = gauss(rng) * kBeamSpot.sigmaZ;
    const auto crossings = follow(kMinPtGeV, phi, cotTheta, charge, vertexZ, scatterers, scattererX2X0, rng);
    if (!crossings.complete) {
      continue;
    }
    const auto timed = [&](std::uint32_t layer) { return TimedPoint{.point = crossings.points[layer], .time = {.start = 0U, .width = 1U}}; };
    std::array<bool, 5> pairPasses{};
    for (std::uint32_t link = 0; link < Det::layerLinks.size(); ++link) {
      const auto& source = crossings.points[Det::layerLinks[link].from];
      const auto& target = crossings.points[Det::layerLinks[link].to];
      const auto terms = pairTerms(parameters, link, beamLine, source, target, true);
      pairPasses[link] = pairCompatible(parameters.pairCut, beamLine, source, terms.source, target, terms.between);
    }
    std::array<bool, 4> tripletPasses{};
    std::array<TripletFitResult, 4> fits{};
    std::array<KinkTerms, 4> terms{};
    for (std::uint32_t first = 0; first + 1U < Det::layerLinks.size(); ++first) {
      if (!pairPasses[first] || !pairPasses[first + 1U]) {
        continue;
      }
      const auto a = timed(first);
      const auto b = timed(first + 1U);
      const auto c = timed(first + 2U);
      auto& count = counts.triplets[first];
      ++count.pairs;
      const auto dip = [](const PairPoint& from, const PairPoint& to) { return std::atan2(static_cast<double>(to.z - from.z), std::hypot(static_cast<double>(to.x - from.x), static_cast<double>(to.y - from.y))); };
      const auto layers = std::array{static_cast<LayerID>(first), static_cast<LayerID>(first + 1U), static_cast<LayerID>(first + 2U)};
      const auto halfWidth = static_cast<double>(tripletLambdaHalfWidth(parameters, first, first + 1U, a.point, b.point, bounds));
      if (std::abs(dip(b.point, c.point) - dip(a.point, b.point)) > halfWidth) {
        continue;
      }
      bool kinked = false;
      if (!tripletCompatible(parameters, layers, a, b, c, fits[first], terms[first], kinked)) {
        continue;
      }
      tripletPasses[first] = true;
      ++count.accepted;
    }
    // GCC 14 addresses `neighbours[first]` as `triplets[4 + first]` once
    // the cuts are inlined and finds that subscript out of bounds;
    // `first` stays below `neighbours.size()`.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
    for (std::uint32_t first = 0; first + 2U < Det::layerLinks.size(); ++first) {
      if (!tripletPasses[first] || !tripletPasses[first + 1U]) {
        continue;
      }
      auto& count = counts.neighbours[first];
      ++count.pairs;
      if (!neighbourCurvatureCompatible(parameters, fits[first].curvature, fits[first].varCurvature, fits[first + 1U].curvature, fits[first + 1U].varCurvature)) {
        continue;
      }
      TripletFitResult fit;
      KinkTerms bcd;
      bool made = false;
      bool valid = false;
      if (neighbourCompatible(parameters, first, first + 1U, first + 2U, timed(first), timed(first + 1U), timed(first + 2U), timed(first + 3U), bcd, made, valid, fit, 1.F, 1.F)) {
        ++count.accepted;
      }
    }
#pragma GCC diagnostic pop
  }
  return counts;
}

} // namespace

/// Only the link across the pipe has it between its layers; a source beyond
/// the pipe has it inside, a source inside the pipe does not.
BOOST_AUTO_TEST_CASE(ThePipeIsBetweenOrInsideByRadius)
{
  IterationSettings settings;
  settings.windows = windowParams();
  const auto withPipe = table(0.01);
  const auto withoutPipe = table(0.);
  const auto parameters = iterationParameters<Det>(settings, kBeamSpot, *withPipe);
  for (std::uint32_t link = 0; link < Det::layerLinks.size(); ++link) {
    BOOST_TEST_CONTEXT("link " << link) { BOOST_TEST((parameters.betweenBound[link] > 0.F) == (link == kCrossingLink)); }
  }

  const auto cut = pairCut(settings.windows);
  const auto beamLine = vertexSeed(kBeamSpot);
  const PairPoint inside{.x = 2.5F, .z = 0.5F, .r = 2.5F};
  const PairPoint beyond{.x = 7.F, .z = 1.4F, .r = 7.F};
  const PairPoint far{.x = 10.F, .z = 2.F, .r = 10.F};
  BOOST_TEST(sourceScattering<Det>(*withPipe, cut, beamLine, 2, inside) == sourceScattering<Det>(*withoutPipe, cut, beamLine, 2, inside));
  BOOST_TEST(sourceScattering<Det>(*withPipe, cut, beamLine, 3, beyond) > sourceScattering<Det>(*withoutPipe, cut, beamLine, 3, beyond));
  const auto between = [&](const LinkMaterial<Det>& material, std::uint32_t link, LayerID from, const PairPoint& source, const PairPoint& target) {
    Scalar spanX = 0.F;
    static_cast<void>(sourceScattering<Det>(material, cut, beamLine, from, source, atan2Of(source.y, source.x), spanX));
    return betweenScattering<Det>(material, cut, betweenSource<Det>(source, from), link, static_cast<LayerID>(from + 1), beamLine, source, target, spanX).between;
  };
  BOOST_TEST(between(*withPipe, kCrossingLink, 2, inside, beyond) > 0.F);
  BOOST_TEST(between(*withoutPipe, kCrossingLink, 2, inside, beyond) == 0.F);
  BOOST_TEST(between(*withPipe, kBeyondLink, 3, beyond, far) == 0.F);
}

/// Tracks scattered as the cut models them pass it at its n-sigma share on
/// every link, and the windows hold the pairs the cut accepts. Without the
/// kink between source and target, the pairs across the pipe fail; without
/// the pipe inside the source, the pairs beyond it.
/// No pair's kink between its layers exceeds the link's bound: the tracklet
/// stage's test at the bound refuses only pairs the pair's own kink would
/// refuse too (ALICE 3's barrel: testTrkChainPulls).
BOOST_AUTO_TEST_CASE(NoPairsKinkExceedsItsLinksBound)
{
  const auto ratios = test::betweenOverBound<Det>(*table(0.01), windowParams(), kBeamSpot);
  BOOST_TEST(ratios.size() == 1U);
  for (const auto& ratio : ratios) {
    BOOST_TEST(ratio.worst <= 1.);
    BOOST_TEST(ratio.worstSourceScale <= 1.);
  }
}

BOOST_AUTO_TEST_CASE(TruePairsPassAcrossAndBeyondThePipe)
{
  constexpr std::int32_t kTracks = 20000;
  constexpr double kPipeX2X0 = 0.01;
  const auto full = countPairs(*table(kPipeX2X0), true, kPipeX2X0, kTracks);
  for (std::uint32_t link = 0; link < full.size(); ++link) {
    const auto& count = full[link];
    BOOST_TEST_MESSAGE("link " << link << ": " << count.accepted << "/" << count.pairs << " accepted, " << count.inPhiWindow << " and " << count.inSecondWindow << " of them in the phi and second windows");
    BOOST_TEST_CONTEXT("link " << link)
    {
      BOOST_TEST(count.pairs > kTracks / 2);
      // 3 sigma: 0.13 % one-sided on the curvature, 0.27 % on the pitch.
      BOOST_TEST(fraction(count.accepted, count.pairs) > 0.993);
      BOOST_TEST(count.inPhiWindow == count.accepted);
      BOOST_TEST(count.inSecondWindow == count.accepted);
    }
  }

  const auto noBetween = countPairs(*table(kPipeX2X0), false, kPipeX2X0, kTracks);
  BOOST_TEST_MESSAGE("without the kink between: link 2 " << noBetween[kCrossingLink].accepted << "/" << noBetween[kCrossingLink].pairs);
  BOOST_TEST(fraction(noBetween[kCrossingLink].accepted, noBetween[kCrossingLink].pairs) < 0.98);

  const auto noPipe = countPairs(*table(0.), true, kPipeX2X0, kTracks);
  BOOST_TEST_MESSAGE("without the pipe in the table: link 3 " << noPipe[kBeyondLink].accepted << "/" << noPipe[kBeyondLink].pairs);
  BOOST_TEST(fraction(noPipe[kBeyondLink].accepted, noPipe[kBeyondLink].pairs) < 0.99);
}

/// Triplets and neighbours with the pipe between their points, of every kind
/// -- between the first two, the middle two, the last two -- pass their cuts
/// at the n-sigma share when the cuts take its kink where it happens, and
/// fail without it.
BOOST_AUTO_TEST_CASE(TripletsAndNeighboursAcrossThePipe)
{
  constexpr std::int32_t kTracks = 20000;
  constexpr double kPipeX2X0 = 0.01;
  const auto full = countChains(*table(kPipeX2X0), kPipeX2X0, kTracks);
  const auto without = countChains(*table(0.), kPipeX2X0, kTracks);
  for (std::uint32_t first = 0; first < full.triplets.size(); ++first) {
    BOOST_TEST_MESSAGE("triplet " << first << "-" << first + 2U << ": " << full.triplets[first].accepted << "/" << full.triplets[first].pairs << ", without the pipe " << without.triplets[first].accepted << "/" << without.triplets[first].pairs);
    BOOST_TEST_CONTEXT("triplet from layer " << first) { BOOST_TEST(fraction(full.triplets[first].accepted, full.triplets[first].pairs) > 0.99); }
  }
  for (std::uint32_t first = 0; first < full.neighbours.size(); ++first) {
    BOOST_TEST_MESSAGE("neighbour " << first << "-" << first + 3U << ": " << full.neighbours[first].accepted << "/" << full.neighbours[first].pairs << ", without the pipe " << without.neighbours[first].accepted << "/" << without.neighbours[first].pairs);
    BOOST_TEST_CONTEXT("neighbour from layer " << first) { BOOST_TEST(fraction(full.neighbours[first].accepted, full.neighbours[first].pairs) > 0.99); }
  }
  // The pipe lies 0.5 cm beyond layer 2 of the 4.5 between layers 2 and 3.
  // Between the middle points it matters: the triplet 1-3 and the neighbour
  // 1-4 fail without it. Just past the first point, t = 0.11, it turns the
  // angle at `b` by 0.11 theta_p and hardly shows (triplet 2-4).
  BOOST_TEST(fraction(without.triplets[1].accepted, without.triplets[1].pairs) < 0.9);
  BOOST_TEST(fraction(without.neighbours[1].accepted, without.neighbours[1].pairs) < 0.95);
}

} // namespace o2::gtrk
