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
///        on exact helices without material: the truth back from exact inputs, unit pulls and the
///        chi2's mean from smeared ones.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK decay vertex fit
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Study/DecayConstraints.h"
#include "GenericTracking/Study/DecayProfile.h"
#include "GenericTracking/Study/DecayVertexFit.h"

#include "DetectorsBase/MatCell.h"

#include <array>
#include <cmath>
#include <random>
#include <tuple>
#include <vector>

namespace o2::gtrk::study
{

namespace
{

/// No material anywhere.
struct NoMaterial {
  [[nodiscard]] o2::base::MatBudget getMatBudget(float, float, float, float, float, float) const { return {}; }
  [[nodiscard]] o2::base::MatMoment getMatMoment(float, float, float, float, float, float) const { return {}; }
};

constexpr std::array<double, 3> kVertex{0.3, -0.2, 1.};
constexpr std::array<double, 5> kSigmas{1e-3, 1e-3, 1e-3, 1e-3, 0.01}; ///< loc0, loc1, phi, theta [cm, rad]; q/p relative

struct Truth {
  double phi, theta, qOverP;
};
const std::array<Truth, 3> kDaughters{Truth{0.3, 1.2, 2.}, Truth{-0.4, 1.6, -1.25}, Truth{0.1, 1.4, -0.8}};

/// A daughter of \a truth from kVertex, exact on a plane 3 cm downstream, its
/// covariance kSigmas's.
DecayDaughter daughterOf(const Truth& truth, const Propagator& prop, const std::array<double, 3>& vertex = kVertex)
{
  const auto start = beamLinePlane(static_cast<Scalar>(vertex[0]), static_cast<Scalar>(vertex[1]), static_cast<Scalar>(truth.phi));
  BoundTrackState state;
  state.par[kLoc0] = 0.F;
  state.par[kLoc1] = static_cast<Scalar>(vertex[2]);
  state.par[kPhi] = static_cast<Scalar>(truth.phi);
  state.par[kTheta] = static_cast<Scalar>(truth.theta);
  state.par[kQOverP] = static_cast<Scalar>(truth.qOverP);
  TrackingSurface target;
  const double c = std::cos(truth.phi);
  const double s = std::sin(truth.phi);
  target.origin[0] = static_cast<Scalar>(vertex[0] + (3. * c));
  target.origin[1] = static_cast<Scalar>(vertex[1] + (3. * s));
  target.origin[2] = static_cast<Scalar>(vertex[2]);
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

/// A K0S-like decay: a neutral mother of mass kMass made at kPoint flies to
/// kDecay; its two pions' momenta sum to its momentum along the flight.
constexpr double kMass = 0.497611;
constexpr double kPion = 0.13957;
constexpr std::array<double, 3> kPoint{0.01, -0.02, 0.5};
constexpr std::array<double, 3> kDecay{1.21, 0.98, 1.7};

std::array<Truth, 2> pionPair()
{
  std::array<double, 3> dir{kDecay[0] - kPoint[0], kDecay[1] - kPoint[1], kDecay[2] - kPoint[2]};
  const double length = std::hypot(dir[0], dir[1], dir[2]);
  for (auto& d : dir) {
    d /= length;
  }
  // A unit vector across the flight, in the transverse plane.
  const double across = std::hypot(dir[0], dir[1]);
  const std::array<double, 3> normal{-dir[1] / across, dir[0] / across, 0.};
  // p1 = a dir + b normal, p2 = c dir - b normal; c from the mass.
  const double a = 0.7;
  const double b = 0.15;
  const auto massOf = [&](double c) {
    const double e1 = std::sqrt((a * a) + (b * b) + (kPion * kPion));
    const double e2 = std::sqrt((c * c) + (b * b) + (kPion * kPion));
    return std::sqrt(((e1 + e2) * (e1 + e2)) - ((a + c) * (a + c)));
  };
  double low = 0.05;
  double high = 5.;
  for (int i = 0; i < 200; ++i) {
    const double mid = 0.5 * (low + high);
    // the mass falls as c grows past a: bisect on the side of the root
    (massOf(mid) > kMass ? low : high) = mid;
  }
  const double c = 0.5 * (low + high);
  const auto truthOf = [](const std::array<double, 3>& p, double charge) {
    const double total = std::hypot(p[0], p[1], p[2]);
    return Truth{std::atan2(p[1], p[0]), std::atan2(std::hypot(p[0], p[1]), p[2]), charge / total};
  };
  const std::array<double, 3> p1{(a * dir[0]) + (b * normal[0]), (a * dir[1]) + (b * normal[1]), (a * dir[2]) + (b * normal[2])};
  const std::array<double, 3> p2{(c * dir[0]) - (b * normal[0]), (c * dir[1]) - (b * normal[1]), (c * dir[2]) - (b * normal[2])};
  return {truthOf(p1, 1.), truthOf(p2, -1.)};
}

} // namespace

/// Every constraint set: the constraints hold, the vertex pulls
/// stay at one and the chi2's mean is its new ndf.
BOOST_AUTO_TEST_CASE(ConstrainedFitsHoldTheirConstraints)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  prop.massGeV = static_cast<Scalar>(kPion);
  const auto pions = pionPair();
  std::vector<DecayDaughter> exactStates;
  for (const auto& truth : pions) {
    exactStates.push_back(daughterOf(truth, prop, kDecay));
  }
  std::mt19937 random{777U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 1000;
  const double fieldScale = 20. * static_cast<double>(constants::math::B2C);
  const std::array<double, 2> masses{kPion, kPion};
  for (const int set : {1, 2, 3}) {
    const DecayConstraints constraints{.mass = (set & 1) != 0, .massGeV = kMass, .production = (set & 2) != 0, .point = kPoint};
    {
      std::array<double, 3> sum{};
      std::array<double, 3> sum2{};
      double chi2 = 0.;
      double worst = 0.;
      int fitted = 0;
      int ndf = 0;
      for (int trial = 0; trial < kTrials; ++trial) {
        auto states = exactStates;
        for (auto& daughter : states) {
          for (int k = 0; k < 5; ++k) {
            daughter.innermost.par[k] += static_cast<Scalar>(gauss(random) * std::sqrt(static_cast<double>(daughter.innermost.covAt(k, k))));
          }
        }
        const auto fit = fitDecayVertex(std::span<const DecayDaughter>{states}, NoMaterial{}, DecayFitVariant::kNoMaterial, constraints);
        if (fit.status != DecayFitStatus::kOk) {
          continue;
        }
        ++fitted;
        ndf = fit.ndf;
        const int n = 3 + 6;
        for (int k = 0; k < 3; ++k) {
          const double pull = (fit.vertex[k] - kDecay[k]) / std::sqrt(fit.covariance[k * n + k]);
          sum[k] += pull;
          sum2[k] += pull * pull;
        }
        chi2 += fit.chi2;
        std::vector<double> values;
        std::vector<double> unused;
        BOOST_REQUIRE(decay_fit_detail::evaluateConstraints(constraints, DecayTopology{}, std::vector{fit.vertex}, fit.momenta, masses, std::array<double, 2>{1., -1.}, fieldScale, values, unused));
        for (const auto value : values) {
          worst = std::max(worst, std::abs(value));
        }
      }
      BOOST_TEST_MESSAGE("constraints " << set << " states" << ": fitted " << fitted << ", chi2 mean " << chi2 / fitted << " at ndf " << ndf << ", largest constraint residual " << worst);
      BOOST_TEST(fitted >= kTrials - 5);
      BOOST_TEST(worst < 1e-5); // a tenth of a micron, 1e-5 GeV^2: float rounding, far below any resolution
      BOOST_TEST(std::abs((chi2 / fitted / ndf) - 1.) < 0.1);
      for (int k = 0; k < 3; ++k) {
        const double mean = sum[k] / fitted;
        const double width = std::sqrt((sum2[k] / fitted) - (mean * mean));
        BOOST_TEST_MESSAGE("  vertex " << k << ": pull mean " << mean << ", width " << width);
        BOOST_TEST(std::abs(mean) < 0.12);
        BOOST_TEST(std::abs(width - 1.) < 0.08);
      }
    }
  }
}

/// A conversion-like pair -- one direction, opposite charges -- under the
/// collinear constraint, profiled along its line (DecayProfile.h) over cells
/// of equal prior about the truth, so that the likelihood alone decides: the
/// constraints hold, the chi2's mean is its ndf, the vertex
/// pulls -- along the flight the profile's spread -- are one.
BOOST_AUTO_TEST_CASE(CollinearDaughtersHoldTheirConstraint)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  prop.massGeV = static_cast<Scalar>(kPion);
  const std::array<Truth, 2> pair{Truth{0.4, 1.3, 1. / 0.3}, Truth{0.4, 1.3, -1. / 0.5}};
  std::vector<DecayDaughter> exactStates;
  for (const auto& truth : pair) {
    exactStates.push_back(daughterOf(truth, prop, kDecay));
  }
  std::mt19937 random{4242U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 1000;
  const double fieldScale = 20. * static_cast<double>(constants::math::B2C);
  const std::array<double, 2> masses{kPion, kPion};
  std::array<double, 3> direction{kDecay[0] - kPoint[0], kDecay[1] - kPoint[1], kDecay[2] - kPoint[2]};
  const double distance = std::hypot(direction[0], direction[1], direction[2]);
  for (auto& d : direction) {
    d /= distance;
  }
  // Cells of equal prior, 1 cm either side of the truth (the likelihood's
  // spread along the flight is about 0.5 mm here), 200 um long.
  std::vector<LineCell> cells;
  for (int k = -50; k < 50; ++k) {
    cells.push_back({.s = distance + ((k + 0.5) * 0.02), .halfLength = 0.01, .probability = 1.});
  }
  {
    std::array<double, 3> sum{};
    std::array<double, 3> sum2{};
    double chi2 = 0.;
    double worst = 0.;
    int fitted = 0;
    int ndf = 0;
    for (int trial = 0; trial < kTrials; ++trial) {
      auto states = exactStates;
      for (auto& daughter : states) {
        for (int k = 0; k < 5; ++k) {
          daughter.innermost.par[k] += static_cast<Scalar>(gauss(random) * std::sqrt(static_cast<double>(daughter.innermost.covAt(k, k))));
        }
      }
      const DecayConstraints constraints{.collinear = true};
      const auto states_ = std::span<const DecayDaughter>{states};
      const auto statesFit = [&](const DecayConstraints& set) { return fitDecayVertex(states_, NoMaterial{}, DecayFitVariant::kNoMaterial, set); };
      const auto fit = profileAlongLine(statesFit, constraints, kPoint, direction, cells).fit;
      if (fit.status != DecayFitStatus::kOk) {
        continue;
      }
      ++fitted;
      ndf = fit.ndf;
      const int n = 3 + 6;
      for (int k = 0; k < 3; ++k) {
        const double pull = (fit.vertex[k] - kDecay[k]) / std::sqrt(fit.covariance[k * n + k]);
        sum[k] += pull;
        sum2[k] += pull * pull;
      }
      chi2 += fit.chi2;
      std::vector<double> values;
      std::vector<double> unused;
      BOOST_REQUIRE(decay_fit_detail::evaluateConstraints(constraints, DecayTopology{}, std::vector{fit.vertex}, fit.momenta, masses, std::array<double, 2>{1., -1.}, fieldScale, values, unused));
      for (const auto value : values) {
        worst = std::max(worst, std::abs(value));
      }
    }
    BOOST_TEST_MESSAGE("collinear" << " states" << ": fitted " << fitted << ", chi2 mean " << chi2 / fitted << " at ndf " << ndf << ", largest constraint residual " << worst);
    BOOST_TEST(fitted >= kTrials - 5);
    BOOST_TEST(worst < 1e-5); // a tenth of a micron, 1e-5 GeV^2: float rounding, far below any resolution
    BOOST_TEST(std::abs((chi2 / fitted / ndf) - 1.) < 0.1);
    for (int k = 0; k < 3; ++k) {
      const double mean = sum[k] / fitted;
      const double width = std::sqrt((sum2[k] / fitted) - (mean * mean));
      BOOST_TEST_MESSAGE("  vertex " << k << ": pull mean " << mean << ", width " << width);
      BOOST_TEST(std::abs(mean) < 0.12);
      BOOST_TEST(std::abs(width - 1.) < 0.1);
    }
  }
}

/// A cascade: a Xi- at kVertex into a Lambda and a pi-, the Lambda 2 cm on
/// along its momentum into p pi-. Fitted as one tree (the Lambda's link
/// constrained): exact inputs give both vertices back; smeared ones unit
/// pulls on both and the chi2's mean at its ndf.
BOOST_AUTO_TEST_CASE(ACascadeFitsAsOneTree)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  const Truth proton{0.25, 1.2, 1. / 1.0};
  const Truth pion{0.45, 1.25, -1. / 0.25};
  const Truth bachelor{-0.3, 1.4, -1. / 0.4};
  const auto vectorOf = [](const Truth& t) {
    const double p = 1. / std::abs(t.qOverP);
    return std::array<double, 3>{p * std::cos(t.phi) * std::sin(t.theta), p * std::sin(t.phi) * std::sin(t.theta), p * std::cos(t.theta)};
  };
  std::array<double, 3> lambda{};
  for (const auto& t : {proton, pion}) {
    const auto v = vectorOf(t);
    for (int c = 0; c < 3; ++c) {
      lambda[c] += v[c];
    }
  }
  const double norm = std::hypot(lambda[0], lambda[1], lambda[2]);
  const std::array<double, 3> second{kVertex[0] + (2. * lambda[0] / norm), kVertex[1] + (2. * lambda[1] / norm), kVertex[2] + (2. * lambda[2] / norm)};
  const DecayTopology topology{.vertexOf = {1, 1, 0}, .parent = {-1, 0}};
  const std::vector<DecayDaughter> exactStates{daughterOf(proton, prop, second), daughterOf(pion, prop, second), daughterOf(bachelor, prop, kVertex)};
  const std::array<std::array<double, 3>, 2> truth{kVertex, second};

  const auto exact = fitDecayVertex(std::span<const DecayDaughter>{exactStates}, NoMaterial{}, DecayFitVariant::kNoMaterial, {}, topology);
  BOOST_REQUIRE(exact.status == DecayFitStatus::kOk);
  for (int w = 0; w < 2; ++w) {
    for (int c = 0; c < 3; ++c) {
      BOOST_TEST(std::abs(exact.vertices[w][c] - truth[w][c]) < 1e-4); // float states: well below a micron
    }
  }

  std::mt19937 random{9090U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 1000;
  {
    std::array<double, 6> sum{};
    std::array<double, 6> sum2{};
    double chi2 = 0.;
    int fitted = 0;
    int ndf = 0;
    for (int trial = 0; trial < kTrials; ++trial) {
      auto states = exactStates;
      for (auto& daughter : states) {
        for (int c = 0; c < 5; ++c) {
          daughter.innermost.par[c] += static_cast<Scalar>(gauss(random) * std::sqrt(static_cast<double>(daughter.innermost.covAt(c, c))));
        }
      }
      const auto fit = fitDecayVertex(std::span<const DecayDaughter>{states}, NoMaterial{}, DecayFitVariant::kNoMaterial, {}, topology);
      if (fit.status != DecayFitStatus::kOk) {
        continue;
      }
      ++fitted;
      ndf = fit.ndf;
      const int n = 6 + 9;
      for (int k = 0; k < 6; ++k) {
        const double pull = (fit.vertices[k / 3][k % 3] - truth[k / 3][k % 3]) / std::sqrt(fit.covariance[k * n + k]);
        sum[k] += pull;
        sum2[k] += pull * pull;
      }
      chi2 += fit.chi2;
    }
    BOOST_TEST_MESSAGE("cascade" << " states" << ": fitted " << fitted << ", chi2 mean " << chi2 / fitted << " at ndf " << ndf);
    BOOST_TEST(fitted >= kTrials - 5);
    BOOST_TEST(std::abs((chi2 / fitted / ndf) - 1.) < 0.1);
    for (int k = 0; k < 6; ++k) {
      const double mean = sum[k] / fitted;
      const double width = std::sqrt((sum2[k] / fitted) - (mean * mean));
      BOOST_TEST_MESSAGE("  vertex " << k / 3 << " coordinate " << k % 3 << ": pull mean " << mean << ", width " << width);
      BOOST_TEST(std::abs(mean) < 0.12);
      BOOST_TEST(std::abs(width - 1.) < 0.08);
    }
  }
  // Sequential: the Lambda fitted first, then the Xi from the bachelor and the
  // Lambda as a neutral pseudo-track.
  std::array<double, 3> sum{};
  std::array<double, 3> sum2{};
  double chi2 = 0.;
  int fitted = 0;
  int ndf = 0;
  for (int trial = 0; trial < kTrials; ++trial) {
    auto states = exactStates;
    for (auto& daughter : states) {
      for (int c = 0; c < 5; ++c) {
        daughter.innermost.par[c] += static_cast<Scalar>(gauss(random) * std::sqrt(static_cast<double>(daughter.innermost.covAt(c, c))));
      }
    }
    const auto lambdaFit = fitDecayVertex(std::span<const DecayDaughter>{states}.first(2), NoMaterial{}, DecayFitVariant::kNoMaterial);
    if (lambdaFit.status != DecayFitStatus::kOk) {
      continue;
    }
    const std::vector<DecayDaughter> xi{states[2], pseudoTrackOf(lambdaFit, 1.115683, 0., prop)};
    const auto fit = fitDecayVertex(std::span<const DecayDaughter>{xi}, NoMaterial{}, DecayFitVariant::kNoMaterial);
    if (fit.status != DecayFitStatus::kOk) {
      continue;
    }
    ++fitted;
    ndf = fit.ndf;
    const int n = 3 + 6;
    for (int c = 0; c < 3; ++c) {
      const double pull = (fit.vertex[c] - kVertex[c]) / std::sqrt(fit.covariance[c * n + c]);
      sum[c] += pull;
      sum2[c] += pull * pull;
    }
    chi2 += fit.chi2;
  }
  BOOST_TEST_MESSAGE("cascade sequential: fitted " << fitted << ", chi2 mean " << chi2 / fitted << " at ndf " << ndf);
  BOOST_TEST(fitted >= kTrials - 5);
  BOOST_TEST(std::abs((chi2 / fitted / ndf) - 1.) < 0.15);
  for (int c = 0; c < 3; ++c) {
    const double mean = sum[c] / fitted;
    const double width = std::sqrt((sum2[c] / fitted) - (mean * mean));
    BOOST_TEST_MESSAGE("  Xi vertex coordinate " << c << ": pull mean " << mean << ", width " << width);
    BOOST_TEST(std::abs(mean) < 0.12);
    BOOST_TEST(std::abs(width - 1.) < 0.08);
  }
  // Tracked: the Xi's crossings of planes 1, 2, ... cm before its decay
  // vertex (the outermost first), one, two or four of them, with and without
  // its production point 6 cm before.
  std::array<double, 3> xi{};
  for (const auto& t : {proton, pion, bachelor}) {
    const auto v = vectorOf(t);
    for (int c = 0; c < 3; ++c) {
      xi[c] += v[c];
    }
  }
  const double xiP = std::hypot(xi[0], xi[1], xi[2]);
  const Truth xiTruth{std::atan2(xi[1], xi[0]), std::atan2(std::hypot(xi[0], xi[1]), xi[2]), -1. / xiP};
  Propagator xiProp = prop;
  xiProp.massGeV = 1.32171F;
  const auto at = beamLinePlane(static_cast<Scalar>(kVertex[0]), static_cast<Scalar>(kVertex[1]), static_cast<Scalar>(xiTruth.phi));
  BoundTrackState xiState;
  xiState.par[kLoc0] = 0.F;
  xiState.par[kLoc1] = static_cast<Scalar>(kVertex[2]);
  xiState.par[kPhi] = static_cast<Scalar>(xiTruth.phi);
  xiState.par[kTheta] = static_cast<Scalar>(xiTruth.theta);
  xiState.par[kQOverP] = static_cast<Scalar>(xiTruth.qOverP);
  const auto planeBefore = [&](double distance, int id) {
    const double c = std::cos(xiTruth.phi);
    const double sn = std::sin(xiTruth.phi);
    TrackingSurface plane;
    plane.origin[0] = static_cast<Scalar>(kVertex[0] - (distance * c));
    plane.origin[1] = static_cast<Scalar>(kVertex[1] - (distance * sn));
    plane.axis0[0] = static_cast<Scalar>(-sn);
    plane.axis0[1] = static_cast<Scalar>(c);
    plane.axis1[2] = 1.F;
    plane.normal[0] = static_cast<Scalar>(c);
    plane.normal[1] = static_cast<Scalar>(sn);
    plane.id = static_cast<SensorID>(id);
    return plane;
  };
  const auto crossing = [&](const TrackingSurface& plane) {
    SolenoidJacobian jacobian;
    BoundTrackState out = xiState;
    const auto status = jacobianAnalytic(xiProp, xiState, at, plane, PropagationDirection::kBackward, jacobian, &out);
    BOOST_REQUIRE(status == PropagationStatus::kOk || status == PropagationStatus::kOutOfBounds);
    return out;
  };
  Scalar madeAt[3]{};
  {
    const auto origin = planeBefore(6., 30);
    const auto out = crossing(origin);
    surfaceToGlobal(origin, out.loc0(), out.loc1(), madeAt);
  }
  for (const int rowCount : {1, 2, 4}) {
    DecayMother tracked{.propagator = xiProp};
    for (int r = 0; r < rowCount; ++r) {
      const auto plane = planeBefore(1. + r, 20 + r);
      const auto out = crossing(plane);
      tracked.planes.push_back(plane);
      tracked.positions.push_back({out.loc0(), out.loc1()});
      tracked.covariance.push_back({25e-8, 0., 25e-8});
    }
    for (const bool withPoint : {false, true}) {
      const DecayConstraints set{.production = withPoint, .point = {madeAt[0], madeAt[1], madeAt[2]}};
      std::array<double, 6> tsum{};
      std::array<double, 6> tsum2{};
      double tchi2 = 0.;
      int tfitted = 0;
      int tndf = 0;
      for (int trial = 0; trial < kTrials; ++trial) {
        auto states = exactStates;
        for (auto& daughter : states) {
          for (int c = 0; c < 5; ++c) {
            daughter.innermost.par[c] += static_cast<Scalar>(gauss(random) * std::sqrt(static_cast<double>(daughter.innermost.covAt(c, c))));
          }
        }
        auto mother = tracked;
        for (auto& position : mother.positions) {
          position[0] += gauss(random) * 5e-4;
          position[1] += gauss(random) * 5e-4;
        }
        const auto fit = fitDecayVertex(std::span<const DecayDaughter>{states}, NoMaterial{}, DecayFitVariant::kNoMaterial, set, topology, std::span<const DecayMother>{&mother, 1});
        if (fit.status != DecayFitStatus::kOk) {
          continue;
        }
        ++tfitted;
        tndf = fit.ndf;
        const int n = 6 + 9;
        for (int k = 0; k < 6; ++k) {
          const double pull = (fit.vertices[k / 3][k % 3] - truth[k / 3][k % 3]) / std::sqrt(fit.covariance[k * n + k]);
          tsum[k] += pull;
          tsum2[k] += pull * pull;
        }
        tchi2 += fit.chi2;
      }
      BOOST_TEST_MESSAGE("cascade tracked, " << rowCount << " rows" << (withPoint ? ", from its point" : "") << ": fitted " << tfitted << ", chi2 mean " << tchi2 / tfitted << " at ndf " << tndf);
      BOOST_TEST(tfitted >= kTrials - 5);
      BOOST_TEST(std::abs((tchi2 / tfitted / tndf) - 1.) < 0.1);
      for (int k = 0; k < 6; ++k) {
        const double mean = tsum[k] / tfitted;
        const double width = std::sqrt((tsum2[k] / tfitted) - (mean * mean));
        BOOST_TEST_MESSAGE("  vertex " << k / 3 << " coordinate " << k % 3 << ": pull mean " << mean << ", width " << width);
        BOOST_TEST(std::abs(mean) < 0.12);
        BOOST_TEST(std::abs(width - 1.) < 0.08);
      }
    }
  }
}

/// A tree with a charged composite: a Xi_c+-like root at kVertex into
/// pi+ pi+ and a Xi-, which curves 2 cm to its own decay into pi- and a
/// Lambda, which flies 3 cm in a line into p pi-. The Xi's link is its helix
/// from its vertex back to the root's (the fit's chain): exact inputs give
/// the three vertices back, smeared ones unit pulls and the chi2's mean at
/// its ndf.
BOOST_AUTO_TEST_CASE(AChargedCompositeLinksOnItsHelix)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  const auto vectorOf = [](const Truth& t) {
    const double p = 1. / std::abs(t.qOverP);
    return std::array<double, 3>{p * std::cos(t.phi) * std::sin(t.theta), p * std::sin(t.phi) * std::sin(t.theta), p * std::cos(t.theta)};
  };
  const auto truthOf = [](const std::array<double, 3>& v, double charge) {
    const double p = std::hypot(v[0], v[1], v[2]);
    return Truth{std::atan2(v[1], v[0]), std::atan2(std::hypot(v[0], v[1]), v[2]), (charge != 0. ? charge : 1.) / p};
  };
  // the Xi- at the root, carried 2 cm along its helix: its decay vertex and
  // its momentum there
  const Truth xiAtRoot{0.2, 1.3, -1. / 1.2};
  Propagator xiProp = prop;
  xiProp.massGeV = 1.32171F;
  const auto start = beamLinePlane(static_cast<Scalar>(kVertex[0]), static_cast<Scalar>(kVertex[1]), static_cast<Scalar>(xiAtRoot.phi));
  BoundTrackState xiState;
  xiState.par[kLoc0] = 0.F;
  xiState.par[kLoc1] = static_cast<Scalar>(kVertex[2]);
  xiState.par[kPhi] = static_cast<Scalar>(xiAtRoot.phi);
  xiState.par[kTheta] = static_cast<Scalar>(xiAtRoot.theta);
  xiState.par[kQOverP] = static_cast<Scalar>(xiAtRoot.qOverP);
  TrackingSurface on;
  on.origin[0] = static_cast<Scalar>(kVertex[0] + (2. * std::cos(xiAtRoot.phi)));
  on.origin[1] = static_cast<Scalar>(kVertex[1] + (2. * std::sin(xiAtRoot.phi)));
  on.axis0[0] = static_cast<Scalar>(-std::sin(xiAtRoot.phi));
  on.axis0[1] = static_cast<Scalar>(std::cos(xiAtRoot.phi));
  on.axis1[2] = 1.F;
  on.normal[0] = static_cast<Scalar>(std::cos(xiAtRoot.phi));
  on.normal[1] = static_cast<Scalar>(std::sin(xiAtRoot.phi));
  SolenoidJacobian jacobian;
  BoundTrackState there = xiState;
  const auto status = jacobianAnalytic(xiProp, xiState, start, on, PropagationDirection::kForward, jacobian, &there);
  BOOST_REQUIRE(status == PropagationStatus::kOk || status == PropagationStatus::kOutOfBounds);
  Scalar xiDecay[3]{};
  surfaceToGlobal(on, there.loc0(), there.loc1(), xiDecay);
  const std::array<double, 3> second{xiDecay[0], xiDecay[1], xiDecay[2]};
  const auto xi = vectorOf(Truth{static_cast<double>(there.phi()), static_cast<double>(there.theta()), static_cast<double>(there.qOverP())});
  // its daughters at its decay: a pi-, the Lambda the rest
  const Truth bachelor{0.5, 1.25, -1. / 0.3};
  const auto b = vectorOf(bachelor);
  const std::array<double, 3> lambda{xi[0] - b[0], xi[1] - b[1], xi[2] - b[2]};
  const double lambdaP = std::hypot(lambda[0], lambda[1], lambda[2]);
  const std::array<double, 3> third{second[0] + (3. * lambda[0] / lambdaP), second[1] + (3. * lambda[1] / lambdaP), second[2] + (3. * lambda[2] / lambdaP)};
  const Truth proton{0.1, 1.2, 1. / 0.8};
  const auto pr = vectorOf(proton);
  const Truth lambdaPion = truthOf({lambda[0] - pr[0], lambda[1] - pr[1], lambda[2] - pr[2]}, -1.);
  const Truth piPlusA{-0.3, 1.5, 1. / 0.6};
  const Truth piPlusB{0.7, 1.1, 1. / 0.4};
  const DecayTopology topology{.vertexOf = {0, 0, 1, 2, 2}, .parent = {-1, 0, 1}};
  const std::vector<DecayDaughter> exactStates{daughterOf(piPlusA, prop, kVertex), daughterOf(piPlusB, prop, kVertex), daughterOf(bachelor, prop, second), daughterOf(proton, prop, third), daughterOf(lambdaPion, prop, third)};
  const std::array<std::array<double, 3>, 3> truth{kVertex, second, third};

  const auto exact = fitDecayVertex(std::span<const DecayDaughter>{exactStates}, NoMaterial{}, DecayFitVariant::kNoMaterial, {}, topology);
  BOOST_REQUIRE(exact.status == DecayFitStatus::kOk);
  for (int w = 0; w < 3; ++w) {
    for (int c = 0; c < 3; ++c) {
      BOOST_TEST(std::abs(exact.vertices[w][c] - truth[w][c]) < 1e-4); // float states
    }
  }
  std::mt19937 random{1234U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 1000;
  std::array<double, 9> sum{};
  std::array<double, 9> sum2{};
  double chi2 = 0.;
  int fitted = 0;
  int ndf = 0;
  for (int trial = 0; trial < kTrials; ++trial) {
    auto states = exactStates;
    for (auto& daughter : states) {
      for (int c = 0; c < 5; ++c) {
        daughter.innermost.par[c] += static_cast<Scalar>(gauss(random) * std::sqrt(static_cast<double>(daughter.innermost.covAt(c, c))));
      }
    }
    const auto fit = fitDecayVertex(std::span<const DecayDaughter>{states}, NoMaterial{}, DecayFitVariant::kNoMaterial, {}, topology);
    if (fit.status != DecayFitStatus::kOk) {
      continue;
    }
    ++fitted;
    ndf = fit.ndf;
    const int n = 9 + 15;
    for (int k = 0; k < 9; ++k) {
      const double pull = (fit.vertices[k / 3][k % 3] - truth[k / 3][k % 3]) / std::sqrt(fit.covariance[k * n + k]);
      sum[k] += pull;
      sum2[k] += pull * pull;
    }
    chi2 += fit.chi2;
  }
  BOOST_TEST_MESSAGE("charged composite: fitted " << fitted << ", chi2 mean " << chi2 / fitted << " at ndf " << ndf);
  BOOST_TEST(fitted >= kTrials - 5);
  BOOST_TEST(std::abs((chi2 / fitted / ndf) - 1.) < 0.1);
  for (int k = 0; k < 9; ++k) {
    const double mean = sum[k] / fitted;
    const double width = std::sqrt((sum2[k] / fitted) - (mean * mean));
    BOOST_TEST_MESSAGE("  vertex " << k / 3 << " coordinate " << k % 3 << ": pull mean " << mean << ", width " << width);
    BOOST_TEST(std::abs(mean) < 0.12);
    BOOST_TEST(std::abs(width - 1.) < 0.08);
  }
}

/// Electrons whose states are two-component mixtures in q/p, as a GSF gives
/// them: a core (75 %, 1 %) and a radiative tail (25 %, |q/p| 15 % up, 5 %
/// wide). Each trial draws a component for each electron and smears its
/// state with it; the mixture handed to the fit is the two corrections of
/// that draw. Under the pair's mass, which tells the combinations apart,
/// the Gaussian sum (`kGaussianSum`) covers the momenta (mean pull^2 near 1)
/// and resolves them better than the collapsed fit (kFull).
BOOST_AUTO_TEST_CASE(ElectronMixturesSumOverTheirComponents)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  prop.massGeV = 0.000511F;
  const std::array<Truth, 2> pair{Truth{0.3, 1.2, 1. / 0.6}, Truth{0.9, 1.35, -1. / 0.4}};
  std::vector<DecayDaughter> exactStates;
  for (const auto& truth : pair) {
    exactStates.push_back(daughterOf(truth, prop, kVertex));
  }
  constexpr std::array<double, 2> kWeight{0.75, 0.25};
  constexpr std::array<double, 2> kShift{0., 0.15}; ///< relative, |q/p| up
  constexpr std::array<double, 2> kWidth{0.01, 0.05};
  double energy = 0.;
  std::array<double, 3> total{};
  for (const auto& truth : pair) {
    const double p = 1. / std::abs(truth.qOverP);
    energy += std::hypot(p, 0.000511);
    total[0] += p * std::cos(truth.phi) * std::sin(truth.theta);
    total[1] += p * std::sin(truth.phi) * std::sin(truth.theta);
    total[2] += p * std::cos(truth.theta);
  }
  const DecayConstraints massOfPair{.mass = true, .massGeV = std::sqrt((energy * energy) - (total[0] * total[0]) - (total[1] * total[1]) - (total[2] * total[2]))};
  for (auto& daughter : exactStates) {
    daughter.mass = 0.000511;
  }
  std::mt19937 random{2468U};
  std::normal_distribution<double> gauss;
  std::uniform_real_distribution<double> uniform;
  constexpr int kTrials = 1000;
  double collapsedRms = 0.;
  for (const auto variant : {DecayFitVariant::kFull, DecayFitVariant::kGaussianSum}) {
    double pull2 = 0.;
    double residual2 = 0.;
    int fitted = 0;
    std::mt19937 local = random;
    for (int trial = 0; trial < kTrials; ++trial) {
      auto states = exactStates;
      for (std::size_t i = 0; i < states.size(); ++i) {
        auto& daughter = states[i];
        const double truthQ = pair[i].qOverP;
        const int drawn = uniform(local) < kWeight[0] ? 0 : 1;
        auto measured = daughter.innermost;
        for (int c = 0; c < 4; ++c) {
          measured.par[c] += static_cast<Scalar>(gauss(local) * std::sqrt(static_cast<double>(measured.covAt(c, c))));
        }
        measured.par[kQOverP] = static_cast<Scalar>(truthQ * (1. + kShift[drawn] + (kWidth[drawn] * gauss(local))));
        GSFMixture<BetheHeitlerLoss::kComponents> mixture;
        for (int j = 0; j < 2; ++j) {
          auto& component = mixture.components[j];
          component.logWeight = static_cast<Scalar>(std::log(kWeight[j]));
          component.state = measured;
          component.state.par[kQOverP] = static_cast<Scalar>(static_cast<double>(measured.par[kQOverP]) / (1. + kShift[j]));
          component.state.covAt(kQOverP, kQOverP) = static_cast<Scalar>(std::pow(kWidth[j] * truthQ, 2));
          component.factors = factorCovariance(component.state.cov);
        }
        mixture.size = 2;
        BoundTrackState collapsed;
        Scalar mode = 0.F;
        detail::gsf::finishMixture(mixture, collapsed, mode);
        daughter.innermost = collapsed;
        daughter.mixture = mixture;
      }
      const auto fit = fitDecayVertex(std::span<const DecayDaughter>{states}, NoMaterial{}, variant, massOfPair);
      if (fit.status != DecayFitStatus::kOk) {
        continue;
      }
      ++fitted;
      const int n = 3 + 6;
      for (int i = 0; i < 2; ++i) {
        const int q = 3 + (3 * i) + 2;
        const double residual = fit.momenta[i][2] - pair[i].qOverP;
        pull2 += residual * residual / fit.covariance[q * n + q];
        residual2 += std::pow(residual / pair[i].qOverP, 2);
      }
    }
    const double meanPull2 = pull2 / (2. * fitted);
    const double rms = std::sqrt(residual2 / (2. * fitted));
    BOOST_TEST_MESSAGE((variant == DecayFitVariant::kGaussianSum ? "gaussian sum" : "collapsed") << ": fitted " << fitted << ", q/p mean pull^2 " << meanPull2 << ", relative q/p rms " << rms);
    BOOST_TEST(fitted >= kTrials - 5);
    if (variant == DecayFitVariant::kGaussianSum) {
      BOOST_TEST(std::abs(meanPull2 - 1.) < 0.15);
      BOOST_TEST(rms < collapsedRms);
    } else {
      collapsedRms = rms;
    }
  }
}

/// A D0-like decay into K- e+ and an unseen neutrino, made 1 mm before its
/// vertex: the neutrino an unmeasured daughter, its momentum fitted under
/// the mother's mass and production point (`missingMomentumStarts`, two
/// branches). Exact inputs: one branch gives the neutrino back. Smeared: on
/// the branch nearer the truth the vertex pulls are one and the chi2's mean
/// is its ndf (one).
BOOST_AUTO_TEST_CASE(AnUnseenNeutrinoIsFittedFromItsMotherMassAndPoint)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  const auto vectorOf = [](const Truth& t) {
    const double p = 1. / std::abs(t.qOverP);
    return std::array<double, 3>{p * std::cos(t.phi) * std::sin(t.theta), p * std::sin(t.phi) * std::sin(t.theta), p * std::cos(t.theta)};
  };
  const Truth kaon{0.35, 1.25, -1. / 1.1};
  const Truth positron{0.05, 1.4, 1. / 0.7};
  const std::array<double, 3> neutrino{0.3, 0.25, 0.1};
  constexpr double kKaon = 0.493677;
  constexpr double kElectron = 0.000511;
  std::array<double, 3> total = neutrino;
  double energy = std::hypot(neutrino[0], neutrino[1], neutrino[2]);
  for (const auto& [truth, mass] : {std::pair{kaon, kKaon}, std::pair{positron, kElectron}}) {
    const auto v = vectorOf(truth);
    energy += std::sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]) + (mass * mass));
    for (int c = 0; c < 3; ++c) {
      total[c] += v[c];
    }
  }
  const double pMother = std::hypot(total[0], total[1], total[2]);
  const double massMother = std::sqrt((energy * energy) - (pMother * pMother));
  const std::array<double, 3> made{kVertex[0] - (0.1 * total[0] / pMother), kVertex[1] - (0.1 * total[1] / pMother), kVertex[2] - (0.1 * total[2] / pMother)};
  Propagator kaonProp = prop;
  kaonProp.massGeV = static_cast<Scalar>(kKaon);
  Propagator positronProp = prop;
  positronProp.massGeV = static_cast<Scalar>(kElectron);
  auto exactK = daughterOf(kaon, kaonProp, kVertex);
  exactK.mass = kKaon;
  auto exactE = daughterOf(positron, positronProp, kVertex);
  exactE.mass = kElectron;
  const std::vector<DecayDaughter> visible{exactK, exactE};
  const std::array<double, 2> visibleMasses{kKaon, kElectron};
  const DecayConstraints constraints{.mass = true, .massGeV = massMother, .production = true, .point = made};
  const auto fitFrom = [&](const std::vector<DecayDaughter>& states, const std::array<double, 3>& start) {
    auto all = states;
    DecayDaughter unseen{.mass = 0., .neutral = true, .unmeasured = true};
    const double p = std::hypot(start[0], start[1], start[2]);
    unseen.innermost.par[kPhi] = static_cast<Scalar>(std::atan2(start[1], start[0]));
    unseen.innermost.par[kTheta] = static_cast<Scalar>(std::atan2(std::hypot(start[0], start[1]), start[2]));
    unseen.innermost.par[kQOverP] = static_cast<Scalar>(1. / p);
    all.push_back(unseen);
    return fitDecayVertex(std::span<const DecayDaughter>{all}, NoMaterial{}, DecayFitVariant::kNoMaterial, constraints);
  };
  // exact: one branch gives the neutrino back
  {
    const auto free = fitDecayVertex(std::span<const DecayDaughter>{visible}, NoMaterial{}, DecayFitVariant::kNoMaterial);
    BOOST_REQUIRE(free.status == DecayFitStatus::kOk);
    const auto starts = missingMomentumStarts(free, visibleMasses, made, massMother, 0.);
    BOOST_REQUIRE(!starts.empty());
    double best = std::numeric_limits<double>::infinity();
    for (const auto& start : starts) {
      const auto fit = fitFrom(visible, start);
      if (fit.status == DecayFitStatus::kOk) {
        const auto& q = fit.momenta[2];
        const double p = 1. / q[2];
        const std::array<double, 3> v{p * std::cos(q[0]) * std::sin(q[1]), p * std::sin(q[0]) * std::sin(q[1]), p * std::cos(q[1])};
        best = std::min(best, std::hypot(v[0] - neutrino[0], v[1] - neutrino[1], v[2] - neutrino[2]));
      }
    }
    BOOST_TEST_MESSAGE("unseen neutrino, exact: nearest branch off by " << best << " GeV");
    // the fit stops within what the float inputs resolve, and the neutrino's
    // momentum is fixed loosely (its sigma of order 0.1 GeV)
    BOOST_TEST(best < 1e-3 * std::hypot(neutrino[0], neutrino[1], neutrino[2]));
  }
  std::mt19937 random{1357U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 1000;
  std::array<double, 3> sum{};
  std::array<double, 3> sum2{};
  double chi2 = 0.;
  int fitted = 0;
  int ndf = 0;
  for (int trial = 0; trial < kTrials; ++trial) {
    auto states = visible;
    for (auto& daughter : states) {
      for (int c = 0; c < 5; ++c) {
        daughter.innermost.par[c] += static_cast<Scalar>(gauss(random) * std::sqrt(static_cast<double>(daughter.innermost.covAt(c, c))));
      }
    }
    const auto free = fitDecayVertex(std::span<const DecayDaughter>{states}, NoMaterial{}, DecayFitVariant::kNoMaterial);
    if (free.status != DecayFitStatus::kOk) {
      continue;
    }
    DecayVertexFitResult nearest;
    double best = std::numeric_limits<double>::infinity();
    for (const auto& start : missingMomentumStarts(free, visibleMasses, made, massMother, 0.)) {
      const auto fit = fitFrom(states, start);
      if (fit.status != DecayFitStatus::kOk) {
        continue;
      }
      const auto& q = fit.momenta[2];
      const double p = 1. / q[2];
      const double off = std::hypot((p * std::cos(q[0]) * std::sin(q[1])) - neutrino[0], (p * std::sin(q[0]) * std::sin(q[1])) - neutrino[1], (p * std::cos(q[1])) - neutrino[2]);
      if (off < best) {
        best = off;
        nearest = fit;
      }
    }
    if (!std::isfinite(best)) {
      continue;
    }
    ++fitted;
    ndf = nearest.ndf;
    const int n = 3 + 9;
    for (int c = 0; c < 3; ++c) {
      const double pull = (nearest.vertex[c] - kVertex[c]) / std::sqrt(nearest.covariance[c * n + c]);
      sum[c] += pull;
      sum2[c] += pull * pull;
    }
    chi2 += nearest.chi2;
  }
  BOOST_TEST_MESSAGE("unseen neutrino: fitted " << fitted << ", chi2 mean " << chi2 / fitted << " at ndf " << ndf);
  BOOST_TEST(fitted >= kTrials - 50);
  BOOST_TEST(std::abs((chi2 / fitted / ndf) - 1.) < 0.15);
  for (int c = 0; c < 3; ++c) {
    const double mean = sum[c] / fitted;
    const double width = std::sqrt((sum2[c] / fitted) - (mean * mean));
    BOOST_TEST_MESSAGE("  vertex " << c << ": pull mean " << mean << ", width " << width);
    BOOST_TEST(std::abs(mean) < 0.12);
    BOOST_TEST(std::abs(width - 1.) < 0.1);
  }
}

/// A D+-like K- pi+ pi+ whose K- pi+ pair is a K*0 at its mass: its mass
/// taken softly by its Breit-Wigner. Exact inputs give the vertex back;
/// smeared ones unit vertex pulls.
BOOST_AUTO_TEST_CASE(AResonanceMassIsSoft)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  constexpr double kKaon = 0.493677;
  constexpr double kKStar = 0.89555;
  constexpr double kWidth = 0.0473;
  // the K- and the pi+ of the K*0: their momenta at the K*'s mass, the
  // pion's magnitude solved for it
  const Truth kaon{0.2, 1.3, -1. / 0.9};
  const double phiPi = 0.55;
  const double thetaPi = 1.2;
  const auto unit = [](double phi, double theta) { return std::array<double, 3>{std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta)}; };
  const auto k = unit(kaon.phi, kaon.theta);
  const double pk = 1. / std::abs(kaon.qOverP);
  const double ek = std::hypot(pk, kKaon);
  const auto u = unit(phiPi, thetaPi);
  const double cosine = (k[0] * u[0]) + (k[1] * u[1]) + (k[2] * u[2]);
  // M^2 = mK^2 + mpi^2 + 2 (ek epi - pk p cos): bisect p
  double low = 1e-3;
  double high = 10.;
  for (int i = 0; i < 200; ++i) {
    const double p = 0.5 * (low + high);
    const double m2 = (kKaon * kKaon) + (kPion * kPion) + (2. * ((ek * std::hypot(p, kPion)) - (pk * p * cosine)));
    (m2 < kKStar * kKStar ? low : high) = p;
  }
  const Truth pion{phiPi, thetaPi, 1. / (0.5 * (low + high))};
  const Truth other{-0.4, 1.5, 1. / 0.6};
  Propagator kaonProp = prop;
  kaonProp.massGeV = static_cast<Scalar>(kKaon);
  Propagator pionProp = prop;
  pionProp.massGeV = static_cast<Scalar>(kPion);
  auto a = daughterOf(kaon, kaonProp, kVertex);
  a.mass = kKaon;
  auto b = daughterOf(pion, pionProp, kVertex);
  b.mass = kPion;
  auto c = daughterOf(other, pionProp, kVertex);
  c.mass = kPion;
  const std::vector<DecayDaughter> exactStates{a, b, c};
  std::mt19937 random{97531U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 1000;
  {
    const DecayConstraints set{.resonances = {Resonance{.daughters = {0, 1}, .massGeV = kKStar, .widthGeV = kWidth}}};
    const auto exact = fitDecayVertex(std::span<const DecayDaughter>{exactStates}, NoMaterial{}, DecayFitVariant::kNoMaterial, set);
    BOOST_REQUIRE(exact.status == DecayFitStatus::kOk);
    for (int i = 0; i < 3; ++i) {
      BOOST_TEST(std::abs(exact.vertex[i] - kVertex[i]) < 1e-4); // float states
    }
    std::array<double, 3> sum{};
    std::array<double, 3> sum2{};
    int fitted = 0;
    for (int trial = 0; trial < kTrials; ++trial) {
      auto states = exactStates;
      for (auto& daughter : states) {
        for (int i = 0; i < 5; ++i) {
          daughter.innermost.par[i] += static_cast<Scalar>(gauss(random) * std::sqrt(static_cast<double>(daughter.innermost.covAt(i, i))));
        }
      }
      const auto fit = fitDecayVertex(std::span<const DecayDaughter>{states}, NoMaterial{}, DecayFitVariant::kNoMaterial, set);
      if (fit.status != DecayFitStatus::kOk) {
        continue;
      }
      ++fitted;
      const int n = 3 + 9;
      for (int i = 0; i < 3; ++i) {
        const double pull = (fit.vertex[i] - kVertex[i]) / std::sqrt(fit.covariance[i * n + i]);
        sum[i] += pull;
        sum2[i] += pull * pull;
      }
    }
    BOOST_TEST_MESSAGE("resonance: fitted " << fitted);
    BOOST_TEST(fitted >= kTrials - 5);
    for (int i = 0; i < 3; ++i) {
      const double mean = sum[i] / fitted;
      const double width = std::sqrt((sum2[i] / fitted) - (mean * mean));
      BOOST_TEST_MESSAGE("  vertex " << i << ": pull mean " << mean << ", width " << width);
      BOOST_TEST(std::abs(mean) < 0.12);
      BOOST_TEST(std::abs(width - 1.) < 0.1);
    }
  }
}

/// A hypertriton-like decay into 3He pi-: the 3He's q/p is 2 / p
/// (`chargeMagnitude`). Exact inputs give the mother's mass back; under it,
/// smeared ones unit vertex pulls and the chi2's mean at its ndf.
BOOST_AUTO_TEST_CASE(ADoublyChargedDaughterHasItsMomentum)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  constexpr double kHelium3 = 2.80839;
  const Truth helium{0.3, 1.25, 2. / 2.4}; // q / p, |q| = 2
  const Truth pion{0.7, 1.4, -1. / 0.25};
  double energy = 0.;
  std::array<double, 3> total{};
  for (const auto& [truth, mass, charge] : {std::tuple{helium, kHelium3, 2.}, std::tuple{pion, kPion, 1.}}) {
    const double p = charge / std::abs(truth.qOverP);
    energy += std::hypot(p, mass);
    total[0] += p * std::cos(truth.phi) * std::sin(truth.theta);
    total[1] += p * std::sin(truth.phi) * std::sin(truth.theta);
    total[2] += p * std::cos(truth.theta);
  }
  const double massMother = std::sqrt((energy * energy) - (total[0] * total[0]) - (total[1] * total[1]) - (total[2] * total[2]));
  Propagator heliumProp = prop;
  heliumProp.massGeV = static_cast<Scalar>(kHelium3);
  Propagator pionProp = prop;
  pionProp.massGeV = static_cast<Scalar>(kPion);
  auto he = daughterOf(helium, heliumProp, kVertex);
  he.mass = kHelium3;
  he.chargeMagnitude = 2.;
  auto pi = daughterOf(pion, pionProp, kVertex);
  pi.mass = kPion;
  const std::vector<DecayDaughter> exactStates{he, pi};
  {
    const auto fit = fitDecayVertex(std::span<const DecayDaughter>{exactStates}, NoMaterial{}, DecayFitVariant::kNoMaterial);
    BOOST_REQUIRE(fit.status == DecayFitStatus::kOk);
    double e = 0.;
    std::array<double, 3> sum{};
    for (int i = 0; i < 2; ++i) {
      const auto& q = fit.momenta[i];
      const double p = std::abs(fit.charges[i]) / std::abs(q[2]);
      e += std::hypot(p, exactStates[i].mass);
      sum[0] += p * std::cos(q[0]) * std::sin(q[1]);
      sum[1] += p * std::sin(q[0]) * std::sin(q[1]);
      sum[2] += p * std::cos(q[1]);
    }
    const double mass = std::sqrt((e * e) - (sum[0] * sum[0]) - (sum[1] * sum[1]) - (sum[2] * sum[2]));
    BOOST_TEST_MESSAGE("doubly charged: mass " << mass << " against " << massMother);
    BOOST_TEST(std::abs(mass - massMother) < 1e-4); // float states
  }
  const DecayConstraints set{.mass = true, .massGeV = massMother};
  std::mt19937 random{8642U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 1000;
  std::array<double, 3> sum{};
  std::array<double, 3> sum2{};
  double chi2 = 0.;
  int fitted = 0;
  int ndf = 0;
  for (int trial = 0; trial < kTrials; ++trial) {
    auto states = exactStates;
    for (auto& daughter : states) {
      for (int c = 0; c < 5; ++c) {
        daughter.innermost.par[c] += static_cast<Scalar>(gauss(random) * std::sqrt(static_cast<double>(daughter.innermost.covAt(c, c))));
      }
    }
    const auto fit = fitDecayVertex(std::span<const DecayDaughter>{states}, NoMaterial{}, DecayFitVariant::kNoMaterial, set);
    if (fit.status != DecayFitStatus::kOk) {
      continue;
    }
    ++fitted;
    ndf = fit.ndf;
    const int n = 3 + 6;
    for (int c = 0; c < 3; ++c) {
      const double pull = (fit.vertex[c] - kVertex[c]) / std::sqrt(fit.covariance[c * n + c]);
      sum[c] += pull;
      sum2[c] += pull * pull;
    }
    chi2 += fit.chi2;
  }
  BOOST_TEST_MESSAGE("doubly charged under its mass: fitted " << fitted << ", chi2 mean " << chi2 / fitted << " at ndf " << ndf);
  BOOST_TEST(fitted >= kTrials - 5);
  BOOST_TEST(std::abs((chi2 / fitted / ndf) - 1.) < 0.1);
  for (int c = 0; c < 3; ++c) {
    const double mean = sum[c] / fitted;
    const double width = std::sqrt((sum2[c] / fitted) - (mean * mean));
    BOOST_TEST_MESSAGE("  vertex " << c << ": pull mean " << mean << ", width " << width);
    BOOST_TEST(std::abs(mean) < 0.12);
    BOOST_TEST(std::abs(width - 1.) < 0.08);
  }
}

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
