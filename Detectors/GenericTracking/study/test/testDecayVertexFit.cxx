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
///        and from their rows (DecayRowsFit.h) on exact helices without
///        material: the truth back from exact inputs, unit pulls and the
///        chi2's mean from smeared ones.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK decay vertex fit
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Study/DecayConstraints.h"
#include "GenericTracking/Study/DecayProfile.h"
#include "GenericTracking/Study/DecayRowsFit.h"
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

/// A daughter of \a truth from kVertex as rows: its exact crossings of six
/// planes 2-7 cm downstream (the library's propagation), 5 um errors.
RowsDaughter rowsOf(const Truth& truth, const Propagator& prop, const std::array<double, 3>& vertex = kVertex)
{
  const auto start = beamLinePlane(static_cast<Scalar>(vertex[0]), static_cast<Scalar>(vertex[1]), static_cast<Scalar>(truth.phi));
  BoundTrackState state;
  state.par[kLoc0] = 0.F;
  state.par[kLoc1] = static_cast<Scalar>(vertex[2]);
  state.par[kPhi] = static_cast<Scalar>(truth.phi);
  state.par[kTheta] = static_cast<Scalar>(truth.theta);
  state.par[kQOverP] = static_cast<Scalar>(truth.qOverP);
  RowsDaughter rows{.propagator = prop};
  for (int k = 0; k < 6; ++k) {
    TrackingSurface plane;
    const double c = std::cos(truth.phi);
    const double s = std::sin(truth.phi);
    const double distance = 2. + k;
    plane.origin[0] = static_cast<Scalar>(vertex[0] + (distance * c));
    plane.origin[1] = static_cast<Scalar>(vertex[1] + (distance * s));
    plane.axis0[0] = static_cast<Scalar>(-s);
    plane.axis0[1] = static_cast<Scalar>(c);
    plane.axis1[2] = 1.F;
    plane.normal[0] = static_cast<Scalar>(c);
    plane.normal[1] = static_cast<Scalar>(s);
    plane.id = static_cast<SensorID>(k + 1);
    SolenoidJacobian jacobian;
    BoundTrackState out = state;
    const auto status = jacobianAnalytic(prop, state, start, plane, PropagationDirection::kForward, jacobian, &out);
    BOOST_REQUIRE(status == PropagationStatus::kOk || status == PropagationStatus::kOutOfBounds);
    rows.planes.push_back(plane);
    rows.positions.push_back({out.loc0(), out.loc1()});
    rows.covariance.push_back({25e-8, 0., 25e-8});
  }
  return rows;
}

/// The rows' fit's start: the truth, as the states' fit would give it.
DecayVertexFitResult truthStart()
{
  DecayVertexFitResult start;
  start.status = DecayFitStatus::kOk;
  start.vertex = kVertex;
  start.vertices = {kVertex};
  for (const auto& truth : kDaughters) {
    start.momenta.push_back({truth.phi, truth.theta, truth.qOverP});
  }
  return start;
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

/// Every constraint set, on both fits: the constraints hold, the vertex pulls
/// stay at one and the chi2's mean is its new ndf.
BOOST_AUTO_TEST_CASE(ConstrainedFitsHoldTheirConstraints)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  prop.massGeV = static_cast<Scalar>(kPion);
  const auto pions = pionPair();
  std::vector<DecayDaughter> exactStates;
  std::vector<RowsDaughter> exactRows;
  for (const auto& truth : pions) {
    exactStates.push_back(daughterOf(truth, prop, kDecay));
    exactRows.push_back(rowsOf(truth, prop, kDecay));
  }
  std::mt19937 random{777U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 1000;
  const double fieldScale = 20. * static_cast<double>(constants::math::B2C);
  const std::array<double, 2> masses{kPion, kPion};
  for (const int set : {1, 2, 3}) {
    const DecayConstraints constraints{.mass = (set & 1) != 0, .massGeV = kMass, .production = (set & 2) != 0, .point = kPoint};
    for (const bool rows : {false, true}) {
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
        auto own = exactRows;
        for (auto& daughter : own) {
          for (auto& position : daughter.positions) {
            position[0] += gauss(random) * 5e-4;
            position[1] += gauss(random) * 5e-4;
          }
        }
        const auto free = fitDecayVertex(std::span<const DecayDaughter>{states}, NoMaterial{}, DecayFitVariant::kNoMaterial);
        const auto fit = rows ? fitDecayRows(std::span<const RowsDaughter>{own}, NoMaterial{}, free, constraints) : fitDecayVertex(std::span<const DecayDaughter>{states}, NoMaterial{}, DecayFitVariant::kNoMaterial, constraints);
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
      BOOST_TEST_MESSAGE("constraints " << set << (rows ? " rows" : " states") << ": fitted " << fitted << ", chi2 mean " << chi2 / fitted << " at ndf " << ndf << ", largest constraint residual " << worst);
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
/// of equal prior about the truth, so that the likelihood alone decides: on
/// both fits the constraints hold, the chi2's mean is its ndf, the vertex
/// pulls -- along the flight the profile's spread -- are one.
BOOST_AUTO_TEST_CASE(CollinearDaughtersHoldTheirConstraint)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  prop.massGeV = static_cast<Scalar>(kPion);
  const std::array<Truth, 2> pair{Truth{0.4, 1.3, 1. / 0.3}, Truth{0.4, 1.3, -1. / 0.5}};
  std::vector<DecayDaughter> exactStates;
  std::vector<RowsDaughter> exactRows;
  for (const auto& truth : pair) {
    exactStates.push_back(daughterOf(truth, prop, kDecay));
    exactRows.push_back(rowsOf(truth, prop, kDecay));
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
  for (const bool rows : {false, true}) {
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
      auto own = exactRows;
      for (auto& daughter : own) {
        for (auto& position : daughter.positions) {
          position[0] += gauss(random) * 5e-4;
          position[1] += gauss(random) * 5e-4;
        }
      }
      const DecayConstraints constraints{.collinear = true};
      const auto states_ = std::span<const DecayDaughter>{states};
      const auto rows_ = std::span<const RowsDaughter>{own};
      const auto statesFit = [&](const DecayConstraints& set) { return fitDecayVertex(states_, NoMaterial{}, DecayFitVariant::kNoMaterial, set); };
      const auto start = profileAlongLine(statesFit, DecayConstraints{}, kPoint, direction, cells);
      const auto rowsFit = [&](const DecayConstraints& set) { return fitDecayRows(rows_, NoMaterial{}, start.fit, set); };
      const auto fit = rows ? profileAlongLine(rowsFit, constraints, kPoint, direction, cells).fit : profileAlongLine(statesFit, constraints, kPoint, direction, cells).fit;
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
    BOOST_TEST_MESSAGE("collinear" << (rows ? " rows" : " states") << ": fitted " << fitted << ", chi2 mean " << chi2 / fitted << " at ndf " << ndf << ", largest constraint residual " << worst);
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
/// along its momentum into p pi-. Both fits as one tree (the Lambda's link
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
  const std::vector<RowsDaughter> exactRows{rowsOf(proton, prop, second), rowsOf(pion, prop, second), rowsOf(bachelor, prop, kVertex)};
  const std::array<std::array<double, 3>, 2> truth{kVertex, second};

  const auto exact = fitDecayVertex(std::span<const DecayDaughter>{exactStates}, NoMaterial{}, DecayFitVariant::kNoMaterial, {}, topology);
  BOOST_REQUIRE(exact.status == DecayFitStatus::kOk);
  const auto exactFromRows = fitDecayRows(std::span<const RowsDaughter>{exactRows}, NoMaterial{}, exact, {}, topology);
  BOOST_REQUIRE(exactFromRows.status == DecayFitStatus::kOk);
  for (int w = 0; w < 2; ++w) {
    for (int c = 0; c < 3; ++c) {
      BOOST_TEST(std::abs(exact.vertices[w][c] - truth[w][c]) < 1e-4);         // float states: well below a micron
      BOOST_TEST(std::abs(exactFromRows.vertices[w][c] - truth[w][c]) < 1e-4); // the rows' float crossings
    }
  }

  std::mt19937 random{9090U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 1000;
  for (const bool rows : {false, true}) {
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
      auto own = exactRows;
      for (auto& daughter : own) {
        for (auto& position : daughter.positions) {
          position[0] += gauss(random) * 5e-4;
          position[1] += gauss(random) * 5e-4;
        }
      }
      const auto free = fitDecayVertex(std::span<const DecayDaughter>{states}, NoMaterial{}, DecayFitVariant::kNoMaterial, {}, topology);
      const auto fit = rows ? fitDecayRows(std::span<const RowsDaughter>{own}, NoMaterial{}, free, {}, topology) : free;
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
    BOOST_TEST_MESSAGE("cascade" << (rows ? " rows" : " states") << ": fitted " << fitted << ", chi2 mean " << chi2 / fitted << " at ndf " << ndf);
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
    const std::vector<DecayDaughter> xi{states[2], pseudoTrackOf(lambdaFit, 1.115683, prop)};
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
  // Tracked: the Xi's own state 1 cm before its decay vertex (its helix from
  // the vertex carried back), one more measurement of the tree.
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
  DecayMother tracked{.propagator = xiProp, .vertex = 0};
  {
    const auto at = beamLinePlane(static_cast<Scalar>(kVertex[0]), static_cast<Scalar>(kVertex[1]), static_cast<Scalar>(xiTruth.phi));
    BoundTrackState state;
    state.par[kLoc0] = 0.F;
    state.par[kLoc1] = static_cast<Scalar>(kVertex[2]);
    state.par[kPhi] = static_cast<Scalar>(xiTruth.phi);
    state.par[kTheta] = static_cast<Scalar>(xiTruth.theta);
    state.par[kQOverP] = static_cast<Scalar>(xiTruth.qOverP);
    TrackingSurface before;
    const double c = std::cos(xiTruth.phi);
    const double sn = std::sin(xiTruth.phi);
    before.origin[0] = static_cast<Scalar>(kVertex[0] - c);
    before.origin[1] = static_cast<Scalar>(kVertex[1] - sn);
    before.origin[2] = static_cast<Scalar>(kVertex[2]);
    before.axis0[0] = static_cast<Scalar>(-sn);
    before.axis0[1] = static_cast<Scalar>(c);
    before.axis1[2] = 1.F;
    before.normal[0] = static_cast<Scalar>(c);
    before.normal[1] = static_cast<Scalar>(sn);
    before.id = 7;
    SolenoidJacobian jacobian;
    BoundTrackState out = state;
    const auto status = jacobianAnalytic(xiProp, state, at, before, PropagationDirection::kBackward, jacobian, &out);
    BOOST_REQUIRE(status == PropagationStatus::kOk || status == PropagationStatus::kOutOfBounds);
    for (auto& entry : out.cov) {
      entry = 0.F;
    }
    for (int i = 0; i < 4; ++i) {
      out.covAt(i, i) = static_cast<Scalar>(kSigmas[i] * kSigmas[i]);
    }
    out.covAt(kQOverP, kQOverP) = static_cast<Scalar>(std::pow(kSigmas[4] * xiTruth.qOverP, 2));
    out.surface = before.id;
    tracked.outermost = out;
    tracked.surface = before;
  }
  {
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
      for (int c = 0; c < 5; ++c) {
        mother.outermost.par[c] += static_cast<Scalar>(gauss(random) * std::sqrt(static_cast<double>(mother.outermost.covAt(c, c))));
      }
      const std::array<DecayMother, 1> mothers{mother};
      const auto fit = fitDecayVertex(std::span<const DecayDaughter>{states}, NoMaterial{}, DecayFitVariant::kNoMaterial, {}, topology, mothers);
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
    BOOST_TEST_MESSAGE("cascade tracked: fitted " << tfitted << ", chi2 mean " << tchi2 / tfitted << " at ndf " << tndf);
    BOOST_TEST(tfitted >= kTrials - 5);
    BOOST_TEST(std::abs((tchi2 / tfitted / tndf) - 1.) < 0.1);
    for (int k = 0; k < 6; ++k) {
      const double mean = tsum[k] / tfitted;
      const double width = std::sqrt((tsum2[k] / tfitted) - (mean * mean));
      BOOST_TEST_MESSAGE("  tracked vertex " << k / 3 << " coordinate " << k % 3 << ": pull mean " << mean << ", width " << width);
      BOOST_TEST(std::abs(mean) < 0.12);
      BOOST_TEST(std::abs(width - 1.) < 0.08);
    }
  }
  // Tracked, from rows: the Xi's crossings of four planes 1-4 cm before its
  // decay vertex, nearest first, in the rows' fit.
  RowsMother trackedRows{.rows = {.propagator = xiProp}, .vertex = 0};
  {
    const auto at = beamLinePlane(static_cast<Scalar>(kVertex[0]), static_cast<Scalar>(kVertex[1]), static_cast<Scalar>(xiTruth.phi));
    BoundTrackState state;
    state.par[kLoc0] = 0.F;
    state.par[kLoc1] = static_cast<Scalar>(kVertex[2]);
    state.par[kPhi] = static_cast<Scalar>(xiTruth.phi);
    state.par[kTheta] = static_cast<Scalar>(xiTruth.theta);
    state.par[kQOverP] = static_cast<Scalar>(xiTruth.qOverP);
    for (int k = 0; k < 4; ++k) {
      TrackingSurface plane;
      const double c = std::cos(xiTruth.phi);
      const double sn = std::sin(xiTruth.phi);
      const double distance = -(1. + k);
      plane.origin[0] = static_cast<Scalar>(kVertex[0] + (distance * c));
      plane.origin[1] = static_cast<Scalar>(kVertex[1] + (distance * sn));
      plane.axis0[0] = static_cast<Scalar>(-sn);
      plane.axis0[1] = static_cast<Scalar>(c);
      plane.axis1[2] = 1.F;
      plane.normal[0] = static_cast<Scalar>(c);
      plane.normal[1] = static_cast<Scalar>(sn);
      plane.id = static_cast<SensorID>(20 + k);
      SolenoidJacobian jacobian;
      BoundTrackState out = state;
      const auto status = jacobianAnalytic(xiProp, state, at, plane, PropagationDirection::kBackward, jacobian, &out);
      BOOST_REQUIRE(status == PropagationStatus::kOk || status == PropagationStatus::kOutOfBounds);
      trackedRows.rows.planes.push_back(plane);
      trackedRows.rows.positions.push_back({out.loc0(), out.loc1()});
      trackedRows.rows.covariance.push_back({25e-8, 0., 25e-8});
    }
  }
  {
    const std::array<RowsMother, 1> exactMothers{trackedRows};
    const auto exactTracked = fitDecayRows(std::span<const RowsDaughter>{exactRows}, NoMaterial{}, exact, {}, topology, exactMothers);
    BOOST_REQUIRE(exactTracked.status == DecayFitStatus::kOk);
    for (int w = 0; w < 2; ++w) {
      for (int c = 0; c < 3; ++c) {
        BOOST_TEST(std::abs(exactTracked.vertices[w][c] - truth[w][c]) < 1e-4);
      }
    }
    std::array<double, 6> tsum{};
    std::array<double, 6> tsum2{};
    double tchi2 = 0.;
    int tfitted = 0;
    int tndf = 0;
    for (int trial = 0; trial < kTrials; ++trial) {
      auto own = exactRows;
      for (auto& daughter : own) {
        for (auto& position : daughter.positions) {
          position[0] += gauss(random) * 5e-4;
          position[1] += gauss(random) * 5e-4;
        }
      }
      auto mother = trackedRows;
      for (auto& position : mother.rows.positions) {
        position[0] += gauss(random) * 5e-4;
        position[1] += gauss(random) * 5e-4;
      }
      const std::array<RowsMother, 1> mothers{mother};
      const auto fit = fitDecayRows(std::span<const RowsDaughter>{own}, NoMaterial{}, exact, {}, topology, mothers);
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
    BOOST_TEST_MESSAGE("cascade tracked rows: fitted " << tfitted << ", chi2 mean " << tchi2 / tfitted << " at ndf " << tndf);
    BOOST_TEST(tfitted >= kTrials - 5);
    BOOST_TEST(std::abs((tchi2 / tfitted / tndf) - 1.) < 0.1);
    for (int k = 0; k < 6; ++k) {
      const double mean = tsum[k] / tfitted;
      const double width = std::sqrt((tsum2[k] / tfitted) - (mean * mean));
      BOOST_TEST_MESSAGE("  tracked rows vertex " << k / 3 << " coordinate " << k % 3 << ": pull mean " << mean << ", width " << width);
      BOOST_TEST(std::abs(mean) < 0.12);
      BOOST_TEST(std::abs(width - 1.) < 0.08);
    }
  }
  // Tracked from rows, made at a point 6 cm before its decay vertex: the
  // walk ends there, the point measured with 5 um.
  {
    const auto at = beamLinePlane(static_cast<Scalar>(kVertex[0]), static_cast<Scalar>(kVertex[1]), static_cast<Scalar>(xiTruth.phi));
    BoundTrackState state;
    state.par[kLoc0] = 0.F;
    state.par[kLoc1] = static_cast<Scalar>(kVertex[2]);
    state.par[kPhi] = static_cast<Scalar>(xiTruth.phi);
    state.par[kTheta] = static_cast<Scalar>(xiTruth.theta);
    state.par[kQOverP] = static_cast<Scalar>(xiTruth.qOverP);
    TrackingSurface origin = trackedRows.rows.planes[0];
    const double c = std::cos(xiTruth.phi);
    const double sn = std::sin(xiTruth.phi);
    origin.origin[0] = static_cast<Scalar>(kVertex[0] - (6. * c));
    origin.origin[1] = static_cast<Scalar>(kVertex[1] - (6. * sn));
    SolenoidJacobian jacobian;
    BoundTrackState out = state;
    const auto status = jacobianAnalytic(xiProp, state, at, origin, PropagationDirection::kBackward, jacobian, &out);
    BOOST_REQUIRE(status == PropagationStatus::kOk || status == PropagationStatus::kOutOfBounds);
    Scalar made[3]{};
    surfaceToGlobal(origin, out.loc0(), out.loc1(), made);
    auto mother = trackedRows;
    mother.pointCovariance = {25e-8, 0., 25e-8};
    std::array<double, 6> psum{};
    std::array<double, 6> psum2{};
    double pchi2 = 0.;
    int pfitted = 0;
    int pndf = 0;
    for (int trial = 0; trial < kTrials; ++trial) {
      auto own = exactRows;
      for (auto& daughter : own) {
        for (auto& position : daughter.positions) {
          position[0] += gauss(random) * 5e-4;
          position[1] += gauss(random) * 5e-4;
        }
      }
      auto smeared = mother;
      for (auto& position : smeared.rows.positions) {
        position[0] += gauss(random) * 5e-4;
        position[1] += gauss(random) * 5e-4;
      }
      const double across = gauss(random) * 5e-4;
      DecayConstraints production{.production = true, .point = {made[0] - (across * sn), made[1] + (across * c), made[2] + (gauss(random) * 5e-4)}};
      const std::array<RowsMother, 1> mothers{smeared};
      const auto fit = fitDecayRows(std::span<const RowsDaughter>{own}, NoMaterial{}, exact, production, topology, mothers);
      if (fit.status != DecayFitStatus::kOk) {
        continue;
      }
      ++pfitted;
      pndf = fit.ndf;
      const int n = 6 + 9;
      for (int k = 0; k < 6; ++k) {
        const double pull = (fit.vertices[k / 3][k % 3] - truth[k / 3][k % 3]) / std::sqrt(fit.covariance[k * n + k]);
        psum[k] += pull;
        psum2[k] += pull * pull;
      }
      pchi2 += fit.chi2;
    }
    BOOST_TEST_MESSAGE("cascade tracked rows from a point: fitted " << pfitted << ", chi2 mean " << pchi2 / pfitted << " at ndf " << pndf);
    BOOST_TEST(pfitted >= kTrials - 5);
    BOOST_TEST(std::abs((pchi2 / pfitted / pndf) - 1.) < 0.1);
    for (int k = 0; k < 6; ++k) {
      const double mean = psum[k] / pfitted;
      const double width = std::sqrt((psum2[k] / pfitted) - (mean * mean));
      BOOST_TEST_MESSAGE("  from a point vertex " << k / 3 << " coordinate " << k % 3 << ": pull mean " << mean << ", width " << width);
      BOOST_TEST(std::abs(mean) < 0.12);
      BOOST_TEST(std::abs(width - 1.) < 0.08);
    }
  }
}

BOOST_AUTO_TEST_CASE(ExactRowsGiveTheTruthBack)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  std::vector<RowsDaughter> daughters;
  for (const auto& truth : kDaughters) {
    daughters.push_back(rowsOf(truth, prop));
  }
  auto start = truthStart();
  start.vertices[0][0] += 0.01; // away from the truth, to find it again
  const auto fit = fitDecayRows(std::span<const RowsDaughter>{daughters}, NoMaterial{}, start);
  BOOST_REQUIRE(fit.status == DecayFitStatus::kOk);
  for (int k = 0; k < 3; ++k) {
    BOOST_TEST(std::abs(fit.vertex[k] - kVertex[k]) < 2e-5);
  }
  for (std::size_t i = 0; i < kDaughters.size(); ++i) {
    BOOST_TEST(std::abs((fit.momenta[i][2] / kDaughters[i].qOverP) - 1.) < 1e-3);
  }
  BOOST_TEST(fit.ndf == (2 * 18) - 12);
}

BOOST_AUTO_TEST_CASE(SmearedRowsGiveUnitPulls)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  std::vector<RowsDaughter> exact;
  for (const auto& truth : kDaughters) {
    exact.push_back(rowsOf(truth, prop));
  }
  std::mt19937 random{54321U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 2000;
  std::array<double, 3> sum{};
  std::array<double, 3> sum2{};
  double chi2 = 0.;
  int fitted = 0;
  for (int trial = 0; trial < kTrials; ++trial) {
    auto daughters = exact;
    for (auto& daughter : daughters) {
      for (auto& position : daughter.positions) {
        position[0] += gauss(random) * 5e-4;
        position[1] += gauss(random) * 5e-4;
      }
    }
    const auto fit = fitDecayRows(std::span<const RowsDaughter>{daughters}, NoMaterial{}, truthStart());
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
    BOOST_TEST_MESSAGE("rows: vertex " << k << ": pull mean " << mean << ", width " << width);
    BOOST_TEST(std::abs(mean) < 0.1);
    BOOST_TEST(std::abs(width - 1.) < 0.05);
  }
  BOOST_TEST_MESSAGE("rows: chi2 mean " << chi2 / fitted << " at ndf 24");
  BOOST_TEST(std::abs((chi2 / fitted) - 24.) < 0.6);
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
