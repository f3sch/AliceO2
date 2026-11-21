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

/// \file testBetheHeitler.cxx
/// \brief The bremsstrahlung loss mixture against Bethe and Heitler's model, and one component's loss applied to a state.
///
/// Bethe and Heitler's model is exact here: -ln z is Gamma-distributed with
/// shape c = t / ln 2, so E[z^k] = (1 + k)^-c and P(z < z0) = Q(c, -ln z0),
/// the regularised upper incomplete gamma function.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Bethe Heitler
#include "GenericTracking/Tracking/BetheHeitler.h"
#include "GenericTracking/Tracking/MeasurementUpdate.h"

#include <boost/math/special_functions/gamma.hpp>
#include <boost/test/unit_test.hpp>
#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>
#include <vector>

namespace o2::gtrk
{
namespace
{

/// The step thicknesses the mixture is checked at: below a layer's to
/// beyond the thickest services, the table's edges included.
std::vector<double> thicknesses()
{
  std::vector<double> t;
  double x = 2e-5;
  while (x < 0.2) {
    t.push_back(x);
    x *= 1.25;
  }
  for (const double edge : {1e-4, 0.2}) {
    t.push_back(edge);
  }
  std::ranges::sort(t);
  return t;
}

double betheHeitlerShape(double t) { return t / std::numbers::ln2; }

/// P(z < z0) under Bethe and Heitler's model.
double betheHeitlerCdf(double t, double z0) { return boost::math::gamma_q(betheHeitlerShape(t), -std::log(z0)); }

double normalCdf(double x) { return 0.5 * std::erfc(-x / std::numbers::sqrt2); }

struct Mixture {
  int size = 0;
  BetheHeitlerComponent components[kBetheHeitlerComponents];

  explicit Mixture(double t) { size = betheHeitlerMixture(static_cast<Scalar>(t), components); }

  [[nodiscard]] double moment(int k) const
  {
    double sum = 0.;
    for (int i = 0; i < size; ++i) {
      const double m = components[i].mean;
      const double v = components[i].variance;
      sum += components[i].weight * (k == 1 ? m : (v + (m * m)));
    }
    return sum;
  }
  /// P(z < z0).
  [[nodiscard]] double cdf(double z0) const
  {
    double sum = 0.;
    for (int i = 0; i < size; ++i) {
      sum += components[i].weight * normalCdf((z0 - components[i].mean) / std::sqrt(components[i].variance));
    }
    return sum;
  }
  /// P(1 / z > y0) of the components carried to 1 / z to first order, as a
  /// loss along the momentum takes them (applyBetheHeitler).
  [[nodiscard]] double linearisedInverseTail(double y0) const
  {
    double sum = 0.;
    for (int i = 0; i < size; ++i) {
      const double m = components[i].mean;
      sum += components[i].weight * normalCdf(-(y0 - (1. / m)) / (std::sqrt(components[i].variance) / (m * m)));
    }
    return sum;
  }
};

std::array<double, kBoundCovSize> randomCovariance(std::mt19937& rng)
{
  std::normal_distribution<double> gauss;
  const double scale[kNBoundParams] = {1e-2, 3e-2, 1e-3, 2e-3, 5e-2};
  double b[kNBoundParams][kNBoundParams];
  for (auto& row : b) {
    for (double& x : row) {
      x = gauss(rng);
    }
  }
  std::array<double, kBoundCovSize> cov{};
  for (int i = 0; i < kNBoundParams; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = i == j ? 0.5 : 0.;
      for (int k = 0; k < kNBoundParams; ++k) {
        sum += b[i][k] * b[j][k];
      }
      cov[covIndex(i, j)] = sum * scale[i] * scale[j];
    }
  }
  return cov;
}

} // namespace

BOOST_AUTO_TEST_CASE(TheRankOneUpdateIsTheDenseSum)
{
  std::mt19937 rng(105);
  std::normal_distribution<double> gauss;
  double worst = 0.;
  for (int trial = 0; trial < 200; ++trial) {
    const auto start = randomCovariance(rng);
    double cov[kBoundCovSize];
    std::ranges::copy(start, cov);
    double upper[kBoundUpperSize];
    double diagonal[kNBoundParams];
    detail::factorUD(cov, upper, diagonal);
    double a[kNBoundParams];
    for (double& x : a) {
      x = gauss(rng) * 0.05;
    }
    const double c = std::exp(gauss(rng));
    double expected[kBoundCovSize];
    for (int i = 0; i < kNBoundParams; ++i) {
      for (int j = 0; j <= i; ++j) {
        expected[covIndex(i, j)] = start[covIndex(i, j)] + (c * a[i] * a[j]);
      }
    }
    detail::rankOneUpdateUD(upper, diagonal, c, a);
    double result[kBoundCovSize];
    detail::covarianceFromUD(upper, diagonal, result);
    for (int i = 0; i < kNBoundParams; ++i) {
      for (int j = 0; j <= i; ++j) {
        const double norm = std::sqrt(expected[covIndex(i, i)] * expected[covIndex(j, j)]);
        worst = std::max(worst, std::abs(result[covIndex(i, j)] - expected[covIndex(i, j)]) / norm);
      }
    }
  }
  BOOST_TEST_MESSAGE(fmt::format("rank-one update: worst |difference| / sqrt(P_ii P_jj) {:.2e}", worst));
  BOOST_TEST(worst < 1e-13);
}

BOOST_AUTO_TEST_CASE(TheMixtureIsAProbability)
{
  for (const double t : thicknesses()) {
    const Mixture mixture(t);
    BOOST_TEST_CONTEXT("t " << t)
    {
      BOOST_TEST(mixture.size == (t < 1e-4 ? 1 : kBetheHeitlerComponents));
      double sum = 0.;
      for (int i = 0; i < mixture.size; ++i) {
        const auto& c = mixture.components[i];
        BOOST_TEST(c.weight > 0.F);
        BOOST_TEST(c.mean > 0.F);
        BOOST_TEST(c.mean <= 1.F);
        BOOST_TEST(c.variance >= 0.F);
        BOOST_TEST((c.variance > 0.F) == (t >= 1e-4));
        sum += c.weight;
      }
      BOOST_TEST(sum == 1., boost::test_tools::tolerance(1e-6));
    }
  }
  BetheHeitlerComponent unused[kBetheHeitlerComponents];
  for (const Scalar t : {0.2001F, -1e-6F, std::numeric_limits<Scalar>::quiet_NaN(), std::numeric_limits<Scalar>::infinity()}) {
    BOOST_TEST(betheHeitlerMixture(t, unused) == 0);
  }
}

/// The mixture has Bethe and Heitler's mean and second moment at every
/// thickness, to float precision; its CDF follows the model's within what
/// the table's fit and its interpolation leave, printed per thickness.
BOOST_AUTO_TEST_CASE(TheMixtureFollowsBetheHeitler)
{
  const double tails[] = {0.5, 0.8, 0.9, 0.95, 0.99};
  BOOST_TEST_MESSAGE("t        n  E[1-z] rel  E[z^2] diff  P(z<z0) mixture / Bethe-Heitler, z0 = 0.5 0.8 0.9 0.95 0.99; P(1/z > 2, 1.25) linearised / Bethe-Heitler");
  constexpr double kFloatEpsilon = std::numeric_limits<Scalar>::epsilon();
  double worstTail = 0.;
  for (const double t : thicknesses()) {
    const Mixture mixture(t);
    const double c = betheHeitlerShape(t);
    const double meanLoss = (1. - mixture.moment(1)) / -std::expm1(-t);
    const double second = mixture.moment(2) - std::pow(3., -c);
    std::string line = fmt::format("{:<8.4g} {}  {:+.2e}  {:+.2e}  ", t, mixture.size, meanLoss - 1., second);
    BOOST_TEST_CONTEXT("t " << t)
    {
      if (t < 1e-4) {
        BOOST_TEST(mixture.moment(1) == 1.);
      } else {
        BOOST_TEST(std::abs(mixture.moment(1) - std::exp(-t)) < 4. * kFloatEpsilon);
        BOOST_TEST(std::abs(second) < 8. * kFloatEpsilon);
        for (const double z0 : tails) {
          const double ratio = mixture.cdf(z0) / betheHeitlerCdf(t, z0);
          line += fmt::format(" {:.3f}", ratio);
          worstTail = std::max(worstTail, std::abs(std::log(ratio)));
        }
        for (const double y0 : {2., 1.25}) {
          line += fmt::format("  {:.2f}", mixture.linearisedInverseTail(y0) / betheHeitlerCdf(t, 1. / y0));
        }
      }
    }
    BOOST_TEST_MESSAGE(line);
  }
  BOOST_TEST_MESSAGE(fmt::format("six components: tail ratio within a factor {:.3f}", std::exp(worstTail)));
  BOOST_TEST(worstTail < std::log(1.25));
}

/// A loss component applied to a state against the dense formula: q/p
/// scaled by s = 1 / z along the momentum, z against it, the covariance S P
/// S plus the loss's variance on q/p.
BOOST_AUTO_TEST_CASE(ALossComponentScalesQOverPAndItsCovariance)
{
  std::mt19937 rng(1050);
  std::uniform_real_distribution<double> uniform;
  double worst = 0.;
  for (int trial = 0; trial < 200; ++trial) {
    const auto cov = randomCovariance(rng);
    BoundTrackState state;
    state.par[kQOverP] = static_cast<Scalar>((uniform(rng) < 0.5 ? -1. : 1.) / (0.05 + uniform(rng)));
    for (int i = 0; i < kBoundCovSize; ++i) {
      state.cov[i] = static_cast<Scalar>(cov[i]);
    }
    auto factors = factorCovariance(state.cov);
    const BetheHeitlerComponent component{.weight = 1.F, .mean = static_cast<Scalar>(0.2 + (0.8 * uniform(rng))), .variance = static_cast<Scalar>(0.05 * uniform(rng))};
    const bool along = trial % 2 == 0;
    const double u = state.qOverP();
    const double z = component.mean;
    const double s = along ? 1. / z : z;
    const double slope = along ? u / (z * z) : u;
    double expected[kBoundCovSize];
    for (int i = 0; i < kNBoundParams; ++i) {
      for (int j = 0; j <= i; ++j) {
        expected[covIndex(i, j)] = static_cast<double>(state.cov[covIndex(i, j)]) * (i == kQOverP ? s : 1.) * (j == kQOverP ? s : 1.);
      }
    }
    expected[covIndex(kQOverP, kQOverP)] += slope * slope * component.variance;
    BOOST_TEST(applyBetheHeitler(component, along ? PropagationDirection::kForward : PropagationDirection::kBackward, state, factors));
    BOOST_TEST(state.qOverP() == u * s, boost::test_tools::tolerance(1e-6));
    Scalar result[kBoundCovSize];
    covarianceFromFactors(factors, result);
    for (int i = 0; i < kNBoundParams; ++i) {
      for (int j = 0; j <= i; ++j) {
        const double norm = std::sqrt(expected[covIndex(i, i)] * expected[covIndex(j, j)]);
        worst = std::max(worst, std::abs(result[covIndex(i, j)] - expected[covIndex(i, j)]) / norm);
      }
    }
  }
  BOOST_TEST_MESSAGE(fmt::format("loss on a state: worst |difference| / sqrt(P_ii P_jj) {:.2e}", worst));
  BOOST_TEST(worst < 1e-5);
}

/// A loss along the momentum, then the same loss against it, without
/// spread: the state and its covariance come back.
BOOST_AUTO_TEST_CASE(AnUnspreadLossIsUndoneAgainstTheMomentum)
{
  std::mt19937 rng(10500);
  const auto cov = randomCovariance(rng);
  BoundTrackState state;
  state.par[kQOverP] = -2.5F;
  for (int i = 0; i < kBoundCovSize; ++i) {
    state.cov[i] = static_cast<Scalar>(cov[i]);
  }
  auto factors = factorCovariance(state.cov);
  const auto before = factors;
  const BetheHeitlerComponent component{.weight = 1.F, .mean = 0.37F, .variance = 0.F};
  BOOST_TEST(applyBetheHeitler(component, PropagationDirection::kForward, state, factors));
  BOOST_TEST(applyBetheHeitler(component, PropagationDirection::kBackward, state, factors));
  BOOST_TEST(state.qOverP() == -2.5F, boost::test_tools::tolerance(1e-6F));
  for (int i = 0; i < kBoundUpperSize; ++i) {
    BOOST_TEST(factors.upper[i] == before.upper[i], boost::test_tools::tolerance(1e-6F));
  }
  for (int i = 0; i < kNBoundParams; ++i) {
    BOOST_TEST(factors.diagonal[i] == before.diagonal[i], boost::test_tools::tolerance(1e-6F));
  }
}

} // namespace o2::gtrk
