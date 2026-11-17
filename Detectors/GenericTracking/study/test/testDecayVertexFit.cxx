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

/// \file testDecayVertexFit.cxx
/// \brief A decay vertex fitted from its daughters' states (DecayVertexFit.h)
///        on exact helices without material: the truth back from exact
///        states, unit pulls and the chi2's mean from smeared ones.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK decay vertex fit
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Study/DecayVertexFit.h"

#include "DetectorsBase/MatCell.h"

#include <array>
#include <cmath>
#include <random>
#include <vector>

namespace o2::gtrk::study
{

namespace
{

/// No material anywhere.
struct NoMaterial {
  [[nodiscard]] o2::base::MatBudget getMatBudget(float, float, float, float, float, float) const { return {}; }
};

constexpr std::array<double, 3> kVertex{0.3, -0.2, 1.};
constexpr std::array<double, 5> kSigmas{1e-3, 1e-3, 1e-3, 1e-3, 0.01}; ///< loc0, loc1, phi, theta [cm, rad]; q/p relative

struct Truth {
  double phi, theta, qOverP;
};
const std::array<Truth, 3> kDaughters{Truth{0.3, 1.2, 2.}, Truth{-0.4, 1.6, -1.25}, Truth{0.1, 1.4, -0.8}};

/// A daughter of \a truth from kVertex, exact on a plane 3 cm downstream, its
/// covariance kSigmas's.
DecayDaughter daughterOf(const Truth& truth, const Propagator& prop)
{
  const auto start = beamLinePlane(static_cast<Scalar>(kVertex[0]), static_cast<Scalar>(kVertex[1]), static_cast<Scalar>(truth.phi));
  BoundTrackState state;
  state.par[kLoc0] = 0.F;
  state.par[kLoc1] = static_cast<Scalar>(kVertex[2]);
  state.par[kPhi] = static_cast<Scalar>(truth.phi);
  state.par[kTheta] = static_cast<Scalar>(truth.theta);
  state.par[kQOverP] = static_cast<Scalar>(truth.qOverP);
  TrackingSurface target;
  const double c = std::cos(truth.phi);
  const double s = std::sin(truth.phi);
  target.origin[0] = static_cast<Scalar>(kVertex[0] + (3. * c));
  target.origin[1] = static_cast<Scalar>(kVertex[1] + (3. * s));
  target.origin[2] = static_cast<Scalar>(kVertex[2]);
  target.axis0[0] = static_cast<Scalar>(-s);
  target.axis0[1] = static_cast<Scalar>(c);
  target.axis1[2] = 1.F;
  target.normal[0] = static_cast<Scalar>(c);
  target.normal[1] = static_cast<Scalar>(s);
  target.id = 1;
  SolenoidJacobian jacobian;
  BoundTrackState out = state;
  const auto status = jacobianAnalytic(prop, state, start, target, PropagationDirection::kForward, jacobian, &out);
  BOOST_REQUIRE(status == PropagationStatus::kOk || status == PropagationStatus::kOutOfBounds);
  for (auto& entry : out.cov) {
    entry = 0.F;
  }
  for (int i = 0; i < 4; ++i) {
    out.covAt(i, i) = static_cast<Scalar>(kSigmas[i] * kSigmas[i]);
  }
  out.covAt(kQOverP, kQOverP) = static_cast<Scalar>(std::pow(kSigmas[4] * truth.qOverP, 2));
  out.surface = target.id;
  return {.innermost = out, .surface = target, .propagator = prop, .mass = prop.massGeV};
}

} // namespace

BOOST_AUTO_TEST_CASE(ExactStatesGiveTheTruthBack)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  std::vector<DecayDaughter> daughters;
  for (const auto& truth : kDaughters) {
    daughters.push_back(daughterOf(truth, prop));
  }
  for (const auto variant : {DecayFitVariant::kFull, DecayFitVariant::kNoMaterial}) {
    const auto fit = fitDecayVertex(std::span<const DecayDaughter>{daughters}, NoMaterial{}, variant);
    BOOST_REQUIRE(fit.status == DecayFitStatus::kOk);
    for (int k = 0; k < 3; ++k) {
      BOOST_TEST(std::abs(fit.vertex[k] - kVertex[k]) < 2e-5); // the states are floats
    }
    for (std::size_t i = 0; i < kDaughters.size(); ++i) {
      BOOST_TEST(std::abs(fit.momenta[i][0] - kDaughters[i].phi) < 1e-5);
      BOOST_TEST(std::abs(fit.momenta[i][1] - kDaughters[i].theta) < 1e-5);
      BOOST_TEST(std::abs((fit.momenta[i][2] / kDaughters[i].qOverP) - 1.) < 1e-5);
    }
    BOOST_TEST(fit.ndf == 3);
  }
}

BOOST_AUTO_TEST_CASE(SmearedStatesGiveUnitPulls)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  std::vector<DecayDaughter> exact;
  for (const auto& truth : kDaughters) {
    exact.push_back(daughterOf(truth, prop));
  }
  std::mt19937 random{12345U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 2000;
  std::array<double, 3> sum{};
  std::array<double, 3> sum2{};
  double chi2 = 0.;
  int fitted = 0;
  for (int trial = 0; trial < kTrials; ++trial) {
    auto daughters = exact;
    for (auto& daughter : daughters) {
      for (int k = 0; k < 5; ++k) {
        daughter.innermost.par[k] += static_cast<Scalar>(gauss(random) * std::sqrt(static_cast<double>(daughter.innermost.covAt(k, k))));
      }
    }
    const auto fit = fitDecayVertex(std::span<const DecayDaughter>{daughters}, NoMaterial{}, DecayFitVariant::kFull);
    if (fit.status != DecayFitStatus::kOk) {
      continue;
    }
    ++fitted;
    const int n = 3 + (3 * static_cast<int>(daughters.size()));
    for (int k = 0; k < 3; ++k) {
      const double pull = (fit.vertex[k] - kVertex[k]) / std::sqrt(fit.covariance[k * n + k]);
      sum[k] += pull;
      sum2[k] += pull * pull;
    }
    chi2 += fit.chi2;
  }
  BOOST_TEST(fitted == kTrials);
  for (int k = 0; k < 3; ++k) {
    const double mean = sum[k] / fitted;
    const double width = std::sqrt((sum2[k] / fitted) - (mean * mean));
    BOOST_TEST_MESSAGE("vertex " << k << ": pull mean " << mean << ", width " << width);
    BOOST_TEST(std::abs(mean) < 0.1);
    BOOST_TEST(std::abs(width - 1.) < 0.05);
  }
  BOOST_TEST_MESSAGE("chi2 mean " << chi2 / fitted << " at ndf 3");
  BOOST_TEST(std::abs((chi2 / fitted) - 3.) < 0.2);
}

} // namespace o2::gtrk::study
