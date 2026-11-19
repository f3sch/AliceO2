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

/// \file testHypothesisRefit.cxx
/// \brief A refit under a species' hypothesis: the PID's fitter and mass, the state carried to the vertex.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Hypothesis Refit
#include "GenericTracking/Tracking/HypothesisRefit.h"

#include "ElectronChain.h"
#include "UniformMaterial.h"

#include "MathUtils/detail/trigonometric.h"

#include <boost/test/unit_test.hpp>
#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <random>
#include <vector>

namespace o2::gtrk
{
namespace
{

using test::ElectronChain;
using test::medianWidth;
constexpr int kLayers = static_cast<int>(ElectronChain::kLayers);

bool sameBits(const Scalar* a, const Scalar* b, int n) { return std::memcmp(a, b, sizeof(Scalar) * static_cast<std::size_t>(n)) == 0; }

double standardDeviation(const std::vector<double>& values)
{
  double mean = 0.;
  for (const double v : values) {
    mean += v / static_cast<double>(values.size());
  }
  double variance = 0.;
  for (const double v : values) {
    variance += (v - mean) * (v - mean) / static_cast<double>(values.size() - 1);
  }
  return std::sqrt(variance);
}

} // namespace

/// Any species but an electron: the Kalman fit at its mass, the same states to
/// the bit; a species of charge two is refused.
BOOST_AUTO_TEST_CASE(AHadronIsTheKalmanFitAtItsMass)
{
  const UniformMaterial material{.x2X0PerCm = 2.5e-3F, .rho = 0.055F};
  ElectronChain chain;
  std::mt19937 rng(1051);
  for (const auto id : {o2::track::PID::Muon, o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton}) {
    const o2::track::PID pid{id};
    while (!chain.generate(2.F, material, true, rng, true)) {
    }
    Propagator prop = chain.propagator;
    prop.massGeV = pid.getMass();
    const auto kalman = fitMeasurements(prop, chain.tracking, material, chain.measurements.data(), kLayers);
    const auto refit = refitWithHypothesis(chain.propagator, chain.tracking, material, chain.measurements.data(), kLayers, pid, 0.F, 0.F);
    BOOST_TEST_CONTEXT(pid.getName())
    {
      BOOST_REQUIRE(refit.status == PropagationStatus::kOk);
      BOOST_TEST(sameBits(refit.innermost.par, kalman.state.par, kNBoundParams));
      BOOST_TEST(sameBits(refit.innermost.cov, kalman.state.cov, kBoundCovSize));
      BOOST_TEST(refit.chi2 == kalman.chi2);
      BOOST_TEST(refit.perigeeMixture.size == 0);
    }
  }
  const auto refused = refitWithHypothesis(chain.propagator, chain.tracking, material, chain.measurements.data(), kLayers, o2::track::PID{o2::track::PID::Alpha}, 0.F, 0.F);
  BOOST_TEST((refused.status == PropagationStatus::kNotFinite));
  BOOST_TEST((refused.failedPass == FitPass::kNone));
}

/// From a vertex through uniform material: pions (no radiation), the Kalman
/// fit carried to the perigee plane about the vertex, unit pulls of d0, z0,
/// phi, theta and q/p there, the vertex at the origin and displaced; electrons (Bethe and
/// Heitler's losses sampled), the mixture's q/p marginal covers the true q/p
/// at the vertex (the PIT flat) and the mode beats the electron-mass Kalman
/// fit carried alike.
BOOST_AUTO_TEST_CASE(TheVertexStateIsCalibrated)
{
  const UniformMaterial material{.x2X0PerCm = 1e-2F, .rho = 0.22F}; // 2 % X0 between layers, from the vertex too
  constexpr int kTracks = 1000;
  ElectronChain chain;
  std::mt19937 rng(1052);
  // Pions, from the origin and from a vertex displaced in the transverse plane.
  chain.propagator.massGeV = o2::track::PID::getMass(o2::track::PID::Pion);
  constexpr std::array<const char*, kNBoundParams> kNames{"d0", "z0", "phi", "theta", "q/p"};
  for (const auto& [x, y] : {std::pair{0.F, 0.F}, std::pair{0.3F, -0.2F}}) {
    chain.vertexX = x;
    chain.vertexY = y;
    std::array<std::vector<double>, kNBoundParams> pulls;
    for (int track = 0; track < kTracks; ++track) {
      const Scalar qOverP = (track % 2 == 0 ? 2.F : -2.F);
      while (!chain.generate(qOverP, material, false, rng, true)) {
      }
      const auto refit = refitWithHypothesis(chain.propagator, chain.tracking, material, chain.measurements.data(), kLayers, o2::track::PID{o2::track::PID::Pion}, x, y);
      BOOST_REQUIRE(refit.status == PropagationStatus::kOk);
      // The truth on the refit's perigee plane: through the vertex, as the
      // track; its direction and q/p there the vertex's.
      const double truth[kNBoundParams] = {0., 0., chain.vertex.phi(), chain.vertex.theta(), chain.vertex.qOverP()};
      for (int i = 0; i < kNBoundParams; ++i) {
        double d = static_cast<double>(refit.perigee.par[i]) - truth[i];
        if (i == kPhi) {
          math_utils::detail::bringToPMPi<double>(d);
        }
        pulls[static_cast<std::size_t>(i)].push_back(d / std::sqrt(static_cast<double>(refit.perigee.variance(i))));
      }
    }
    std::string line = fmt::format("pions from the vertex at ({}, {}) cm, at the vertex, pull mean and standard deviation:", x, y);
    for (int i = 0; i < kNBoundParams; ++i) {
      const auto& p = pulls[static_cast<std::size_t>(i)];
      double mean = 0.;
      for (const double v : p) {
        mean += v / static_cast<double>(p.size());
      }
      const double width = standardDeviation(p);
      line += fmt::format("  {} {:+.3f} {:.3f}", kNames[static_cast<std::size_t>(i)], mean, width);
      BOOST_TEST_CONTEXT("vertex (" << x << ", " << y << "), " << kNames[static_cast<std::size_t>(i)] << " mean " << mean << " width " << width)
      {
        BOOST_TEST(std::abs(mean) < 0.1);
        BOOST_TEST(std::abs(width - 1.) < 0.1);
      }
    }
    BOOST_TEST_MESSAGE(line);
  }
  chain.vertexX = 0.F;
  chain.vertexY = 0.F;
  // Electrons.
  chain.propagator.massGeV = o2::track::PID::getMass(o2::track::PID::Electron);
  std::vector<double> pit;
  std::vector<double> mode;
  std::vector<double> kalman;
  std::vector<double> d0Pull;
  for (int track = 0; track < kTracks; ++track) {
    const Scalar qOverP = (track % 2 == 0 ? 2.F : -2.F);
    while (!chain.generate(qOverP, material, true, rng, true)) {
    }
    const auto* measurements = chain.measurements.data();
    const auto electron = refitWithHypothesis(chain.propagator, chain.tracking, material, measurements, kLayers, o2::track::PID{o2::track::PID::Electron}, 0.F, 0.F);
    // The electron-mass Kalman fit carried alike: the hadron path at the electron's mass.
    const auto kalmanFit = fitMeasurements(chain.propagator, chain.tracking, material, measurements, kLayers);
    if (electron.status != PropagationStatus::kOk || kalmanFit.status != PropagationStatus::kOk) {
      continue;
    }
    BoundTrackState carried = kalmanFit.state;
    auto factors = factorCovariance(carried.cov);
    const auto& innermostSurface = chain.tracking.planeOf(measurements[0]);
    const auto plane = detail::hypothesis::perigeePlaneOf(chain.propagator, carried, innermostSurface, 0.F, 0.F);
    if (detail::hypothesis::carryInward(chain.propagator, material, innermostSurface, plane, carried, factors) != PropagationStatus::kOk) {
      continue;
    }
    const double u = chain.vertex.qOverP();
    pit.push_back(marginalCdf(electron.perigeeMixture, u));
    mode.push_back(electron.perigeeMomentum * std::abs(u));
    kalman.push_back(std::abs(u / static_cast<double>(carried.qOverP())));
    d0Pull.push_back(static_cast<double>(electron.perigee.loc0()) / std::sqrt(static_cast<double>(electron.perigee.variance(kLoc0))));
  }
  const auto tails = [](const std::vector<double>& v) { return static_cast<double>(std::ranges::count_if(v, [](double p) { return p < 0.05 || p > 0.95; })) / static_cast<double>(v.size()); };
  const auto low = [](const std::vector<double>& v) { return static_cast<double>(std::ranges::count_if(v, [](double r) { return r < 0.9; })) / static_cast<double>(v.size()); };
  double pitMean = 0.;
  for (const double p : pit) {
    pitMean += p / static_cast<double>(pit.size());
  }
  const auto [modeMedian, modeWidth] = medianWidth(mode);
  const auto [kalmanMedian, kalmanWidth] = medianWidth(kalman);
  const auto [d0Median, d0Width] = medianWidth(d0Pull);
  BOOST_TEST_MESSAGE(fmt::format(
    "electrons at the vertex: {} fits, q/p PIT mean {:.3f}, outside 5-95 % {:.1f} %; p / p_true median, width, < 0.9: mode {:.3f} {:.3f} {:.1f} %, "
    "Kalman {:.3f} {:.3f} {:.1f} %; d0 pull median {:+.3f}, robust width {:.3f}",
    pit.size(), pitMean, 100. * tails(pit), modeMedian, modeWidth, 100. * low(mode), kalmanMedian, kalmanWidth, 100. * low(kalman), d0Median, d0Width));
  BOOST_TEST(pit.size() > 950U);
  BOOST_TEST(std::abs(pitMean - 0.5) < 0.03);
  BOOST_TEST(std::abs(tails(pit) - 0.1) < 0.035);
  BOOST_TEST(low(mode) < 0.8 * low(kalman));
  BOOST_TEST(modeWidth < 0.9 * kalmanWidth);
}

} // namespace o2::gtrk
