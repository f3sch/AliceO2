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
#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Tracking/DecayFitBatch.h"
#include "GenericTracking/Tracking/DecayFitConstraints.h"
#include "GenericTracking/Tracking/DecayVertexFit.h"
#include "GenericTracking/Tracking/PrimaryVertexFit.h"

#include "UniformMaterial.h"

#include "DetectorsBase/MatCell.h"

#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <vector>

namespace o2::gtrk::study
{

using detectors::TRKFT3;

namespace
{

constexpr Vec3D kVertex{0.3, -0.2, 1.};
constexpr std::array<double, 5> kSigmas{1e-3, 1e-3, 1e-3, 1e-3, 0.01}; ///< loc0, loc1, phi, theta [cm, rad]; q/p relative

struct Truth {
  double phi, theta, qOverP;
};

/// \a truth's momentum vector [GeV], computed apart from the fit's.
Vec3D vectorOf(const Truth& truth)
{
  const double p = 1. / std::abs(truth.qOverP);
  return {p * std::cos(truth.phi) * std::sin(truth.theta), p * std::sin(truth.phi) * std::sin(truth.theta), p * std::cos(truth.theta)};
}

/// The truth of momentum vector \a v and \a charge (q/p 1/p for a neutral one).
Truth truthOf(const Vec3D& v, double charge)
{
  const double p = std::hypot(v[0], v[1], v[2]);
  return Truth{.phi = phiOf(v), .theta = thetaOf(v), .qOverP = (charge != 0. ? charge : 1.) / p};
}
const std::array<Truth, 3> kDaughters{Truth{.phi = 0.3, .theta = 1.2, .qOverP = 2.}, Truth{.phi = -0.4, .theta = 1.6, .qOverP = -1.25}, Truth{.phi = 0.1, .theta = 1.4, .qOverP = -0.8}};

/// A daughter of \a truth from kVertex, exact on a plane 3 cm downstream, its
/// covariance kSigmas's.
DecayDaughter daughterOf(const Truth& truth, const Propagator& prop, const Vec3D& vertex = kVertex)
{
  const auto start = beamLinePlane(static_cast<Scalar>(vertex[0]), static_cast<Scalar>(vertex[1]), static_cast<Scalar>(truth.phi));
  BoundTrackState state;
  state.par[kLoc0] = 0.F;
  state.par[kLoc1] = static_cast<Scalar>(vertex[2]);
  state.par[kPhi] = static_cast<Scalar>(truth.phi);
  state.par[kTheta] = static_cast<Scalar>(truth.theta);
  state.par[kQOverP] = static_cast<Scalar>(truth.qOverP);
  auto target = beamLinePlane(static_cast<Scalar>(vertex[0] + (3. * std::cos(truth.phi))), static_cast<Scalar>(vertex[1] + (3. * std::sin(truth.phi))), static_cast<Scalar>(truth.phi));
  target.origin[2] = static_cast<Scalar>(vertex[2]);
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
constexpr Vec3D kPoint{0.01, -0.02, 0.5};
constexpr Vec3D kDecay{1.21, 0.98, 1.7};

std::array<Truth, 2> pionPair()
{
  Vec3D dir{kDecay[0] - kPoint[0], kDecay[1] - kPoint[1], kDecay[2] - kPoint[2]};
  const double length = std::hypot(dir[0], dir[1], dir[2]);
  for (auto& d : dir) {
    d /= length;
  }
  // A unit vector across the flight, in the transverse plane.
  const double across = std::hypot(dir[0], dir[1]);
  const Vec3D normal{-dir[1] / across, dir[0] / across, 0.};
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
  const Vec3D p1{(a * dir[0]) + (b * normal[0]), (a * dir[1]) + (b * normal[1]), (a * dir[2]) + (b * normal[2])};
  const Vec3D p2{(c * dir[0]) - (b * normal[0]), (c * dir[1]) - (b * normal[1]), (c * dir[2]) - (b * normal[2])};
  return {truthOf(p1, 1.), truthOf(p2, -1.)};
}

/// \a name's fit in double and in the reduced-precision configurations
/// (`DecayVertexFitter`'s scalars): the double fit as `fitDecayVertex` gives
/// it (asserted), and each one's Delta (the double fit's chi2 under its last
/// model at that fit's parameters, less at its own) over the double fit's R,
/// as a message.
template <DecayFitVariant Variant = DecayFitVariant::kFull>
void reportPrecision(const char* name, std::span<const DecayDaughter> daughters, const DecayConstraints& constraints = {}, const DecayTopology& topology = {})
{
  const UniformMaterial material;
  std::array<Vec3D, decay_fit_detail::kMaxConstrainedVertices> starts{};
  const std::span<const Vec3D> start{starts.data(), decayFitStart<TRKFT3>(daughters, constraints, topology, {}, starts) ? static_cast<std::size_t>(topology.vertices()) : 0U};
  DecayVertexFitter<TRKFT3, Variant, UniformMaterial> reference(daughters, material, constraints, topology, {});
  std::vector<std::byte> workspace(reference.workspace());
  DecayVertexFitResult fit;
  reference.fit(start, workspace, fit);
  const auto plain = fitDecayVertex<TRKFT3, Variant>(daughters, material, constraints, topology);
  BOOST_TEST(static_cast<int>(fit.status) == static_cast<int>(plain.status));
  BOOST_TEST(fit.chi2 == plain.chi2);
  BOOST_TEST((fit.vertex == plain.vertex));
  if (fit.status != VertexFitStatus::kOk) {
    BOOST_TEST_MESSAGE(name << ": the double fit itself ends with status " << static_cast<int>(fit.status));
    return;
  }
  const std::vector<double> solution(reference.solution().begin(), reference.solution().end());
  double chi2 = 0.;
  double violation = 0.;
  BOOST_REQUIRE(reference.frozenMerit(solution, chi2, violation));
  const auto one = [&]<typename Model, typename Solve>(const char* label) {
    DecayVertexFitter<TRKFT3, Variant, UniformMaterial, Model, Solve> fitter(daughters, material, constraints, topology, {});
    std::vector<std::byte> fitterWorkspace(fitter.workspace());
    DecayVertexFitResult result;
    fitter.fit(start, fitterWorkspace, result);
    double at = std::numeric_limits<double>::quiet_NaN();
    double atViolation = std::numeric_limits<double>::quiet_NaN();
    if (result.status == VertexFitStatus::kOk && !reference.frozenMerit(fitter.solution(), at, atViolation)) {
      at = std::numeric_limits<double>::quiet_NaN();
    }
    BOOST_TEST_MESSAGE(name << ", " << label << ": status " << static_cast<int>(result.status) << ", Delta / R " << (at - chi2) / static_cast<double>(fit.stopBound) << " (R " << fit.stopBound << "), |c| " << violation << " -> " << atViolation << ", vertex shift [um] "
                            << 1e4 * std::hypot(result.vertex[0] - fit.vertex[0], result.vertex[1] - fit.vertex[1], result.vertex[2] - fit.vertex[2]));
  };
  one.template operator()<float, double>("model float");
  one.template operator()<double, float>("solve float");
  one.template operator()<float, float>("both float");
}

/// Pulls' mean and width (the population's: the second moment less the
/// mean's square) over \a n of them.
struct PullMoments {
  double sum = 0.;
  double sum2 = 0.;
  void add(double pull)
  {
    sum += pull;
    sum2 += pull * pull;
  }
  [[nodiscard]] double mean(int n) const { return sum / n; }
  [[nodiscard]] double width(int n) const { return std::sqrt((sum2 / n) - (mean(n) * mean(n))); }
};

/// \a states smeared by their covariances' diagonals, drawn from \a random.
std::vector<DecayDaughter> smeared(std::vector<DecayDaughter> states, std::mt19937& random, std::normal_distribution<double>& gauss)
{
  for (auto& daughter : states) {
    if (daughter.unmeasured) {
      continue;
    }
    for (int k = 0; k < 5; ++k) {
      daughter.innermost.par[k] += static_cast<Scalar>(gauss(random) * std::sqrt(static_cast<double>(daughter.innermost.covAt(k, k))));
    }
  }
  return states;
}

/// \a states smeared by their covariances' diagonals, the same draw for a
/// seed.
std::vector<DecayDaughter> smeared(std::vector<DecayDaughter> states, unsigned seed)
{
  std::mt19937 random{seed};
  std::normal_distribution<double> gauss;
  return smeared(std::move(states), random, gauss);
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
  exactStates.reserve(pions.size());
  for (const auto& truth : pions) {
    exactStates.push_back(daughterOf(truth, prop, kDecay));
  }
  std::mt19937 random{777U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 1000;
  const Vec2D masses{kPion, kPion};
  for (const int set : {1, 2, 3}) {
    const DecayConstraints constraints{.mass = (set & 1) != 0, .massGeV = kMass, .production = (set & 2) != 0, .point = kPoint};
    {
      std::array<PullMoments, 3> pulls{};
      double chi2 = 0.;
      double worst = 0.;
      int fitted = 0;
      int ndf = 0;
      for (int trial = 0; trial < kTrials; ++trial) {
        auto states = smeared(exactStates, random, gauss);
        const auto fit = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{states}, UniformMaterial{}, constraints);
        if (fit.status != VertexFitStatus::kOk) {
          continue;
        }
        ++fitted;
        ndf = fit.ndf;
        for (int k = 0; k < 3; ++k) {
          const double pull = (fit.vertex[k] - kDecay[k]) / std::sqrt(fit.covarianceAt(k, k));
          pulls[k].add(pull);
        }
        chi2 += fit.chi2;
        std::array<double, decay_fit_detail::kMaxHardRows> values{};
        std::array<double, static_cast<std::size_t>(decay_fit_detail::kMaxHardRows) * decay_fit_detail::kConstraintInputs> unused{};
        int written = 0;
        BOOST_REQUIRE(decay_fit_detail::evaluateConstraints(constraints, DecayTopology{}, std::vector{fit.vertex}, fit.momenta, masses, Vec2D{1., -1.}, values.data(), unused.data(), decay_fit_detail::kConstraintInputs, decay_fit_detail::kMaxHardRows, written));
        for (int j = 0; j < written; ++j) {
          worst = std::max(worst, std::abs(values[j]));
        }
      }
      BOOST_TEST_MESSAGE("constraints " << set << " states" << ": fitted " << fitted << ", chi2 mean " << chi2 / fitted << " at ndf " << ndf << ", largest constraint residual " << worst);
      BOOST_TEST(fitted >= kTrials - 5);
      BOOST_TEST(worst < 1e-5); // a tenth of a micron, 1e-5 GeV^2: float rounding, far below any resolution
      BOOST_TEST(std::abs((chi2 / fitted / ndf) - 1.) < 0.1);
      for (int k = 0; k < 3; ++k) {
        const double mean = pulls[k].mean(fitted);
        const double width = pulls[k].width(fitted);
        BOOST_TEST_MESSAGE("  vertex " << k << ": pull mean " << mean << ", width " << width);
        BOOST_TEST(std::abs(mean) < 0.12);
        BOOST_TEST(std::abs(width - 1.) < 0.08);
      }
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
  const Truth proton{.phi = 0.25, .theta = 1.2, .qOverP = 1.};
  const Truth pion{.phi = 0.45, .theta = 1.25, .qOverP = -1. / 0.25};
  const Truth bachelor{.phi = -0.3, .theta = 1.4, .qOverP = -1. / 0.4};
  Vec3D lambda{};
  for (const auto& t : {proton, pion}) {
    const auto v = vectorOf(t);
    for (int c = 0; c < 3; ++c) {
      lambda[c] += v[c];
    }
  }
  const double norm = std::hypot(lambda[0], lambda[1], lambda[2]);
  const Vec3D second{kVertex[0] + (2. * lambda[0] / norm), kVertex[1] + (2. * lambda[1] / norm), kVertex[2] + (2. * lambda[2] / norm)};
  const DecayTopology topology{.vertexOf = {1, 1, 0}, .parent = {-1, 0}};
  const std::vector<DecayDaughter> exactStates{daughterOf(proton, prop, second), daughterOf(pion, prop, second), daughterOf(bachelor, prop, kVertex)};
  const std::array<Vec3D, 2> truth{kVertex, second};

  const auto exact = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{exactStates}, UniformMaterial{}, {}, topology);
  BOOST_REQUIRE(exact.status == VertexFitStatus::kOk);
  for (int w = 0; w < 2; ++w) {
    for (int c = 0; c < 3; ++c) {
      BOOST_TEST(std::abs(exact.vertices[w][c] - truth[w][c]) < 1e-4); // float states: well below a micron
    }
  }

  std::mt19937 random{9090U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 1000;
  {
    std::array<PullMoments, 6> pulls{};
    double chi2 = 0.;
    int fitted = 0;
    int ndf = 0;
    for (int trial = 0; trial < kTrials; ++trial) {
      auto states = smeared(exactStates, random, gauss);
      const auto fit = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{states}, UniformMaterial{}, {}, topology);
      if (fit.status != VertexFitStatus::kOk) {
        continue;
      }
      ++fitted;
      ndf = fit.ndf;
      for (int k = 0; k < 6; ++k) {
        const double pull = (fit.vertices[k / 3][k % 3] - truth[k / 3][k % 3]) / std::sqrt(fit.covarianceAt(k, k));
        pulls[k].add(pull);
      }
      chi2 += fit.chi2;
    }
    BOOST_TEST_MESSAGE("cascade" << " states" << ": fitted " << fitted << ", chi2 mean " << chi2 / fitted << " at ndf " << ndf);
    BOOST_TEST(fitted >= kTrials - 5);
    BOOST_TEST(std::abs((chi2 / fitted / ndf) - 1.) < 0.1);
    for (int k = 0; k < 6; ++k) {
      const double mean = pulls[k].mean(fitted);
      const double width = pulls[k].width(fitted);
      BOOST_TEST_MESSAGE("  vertex " << k / 3 << " coordinate " << k % 3 << ": pull mean " << mean << ", width " << width);
      BOOST_TEST(std::abs(mean) < 0.12);
      BOOST_TEST(std::abs(width - 1.) < 0.08);
    }
  }
  // Tracked: the Xi's crossings of planes 1, 2, ... cm before its decay
  // vertex (the outermost first), one, two or four of them, with and without
  // its production point 6 cm before.
  Vec3D xi{};
  for (const auto& t : {proton, pion, bachelor}) {
    const auto v = vectorOf(t);
    for (int c = 0; c < 3; ++c) {
      xi[c] += v[c];
    }
  }
  const double xiP = std::hypot(xi[0], xi[1], xi[2]);
  const Truth xiTruth{.phi = phiOf(xi), .theta = thetaOf(xi), .qOverP = -1. / xiP};
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
    auto plane = beamLinePlane(static_cast<Scalar>(kVertex[0] - (distance * std::cos(xiTruth.phi))), static_cast<Scalar>(kVertex[1] - (distance * std::sin(xiTruth.phi))), static_cast<Scalar>(xiTruth.phi));
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
    std::vector<TrackingSurface> planes;
    std::vector<Vec2D> exactPositions;
    std::vector<Vec3D> covariance;
    for (int r = 0; r < rowCount; ++r) {
      const auto plane = planeBefore(1. + r, 20 + r);
      const auto out = crossing(plane);
      planes.push_back(plane);
      exactPositions.push_back({out.loc0(), out.loc1()});
      covariance.push_back({25e-8, 0., 25e-8});
    }
    for (const bool withPoint : {false, true}) {
      const DecayConstraints set{.production = withPoint, .point = {madeAt[0], madeAt[1], madeAt[2]}};
      std::array<PullMoments, 6> tpulls{};
      double tchi2 = 0.;
      int tfitted = 0;
      int tndf = 0;
      for (int trial = 0; trial < kTrials; ++trial) {
        auto states = smeared(exactStates, random, gauss);
        auto positions = exactPositions;
        for (auto& position : positions) {
          position[0] += gauss(random) * 5e-4;
          position[1] += gauss(random) * 5e-4;
        }
        const DecayMother mother{.planes = planes, .positions = positions, .covariance = covariance, .propagator = xiProp};
        const auto fit = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{states}, UniformMaterial{}, set, topology, std::span<const DecayMother>{&mother, 1});
        if (fit.status != VertexFitStatus::kOk) {
          continue;
        }
        ++tfitted;
        tndf = fit.ndf;
        for (int k = 0; k < 6; ++k) {
          const double pull = (fit.vertices[k / 3][k % 3] - truth[k / 3][k % 3]) / std::sqrt(fit.covarianceAt(k, k));
          tpulls[k].add(pull);
        }
        tchi2 += fit.chi2;
      }
      BOOST_TEST_MESSAGE("cascade tracked, " << rowCount << " rows" << (withPoint ? ", from its point" : "") << ": fitted " << tfitted << ", chi2 mean " << tchi2 / tfitted << " at ndf " << tndf);
      BOOST_TEST(tfitted >= kTrials - 5);
      BOOST_TEST(std::abs((tchi2 / tfitted / tndf) - 1.) < 0.1);
      for (int k = 0; k < 6; ++k) {
        const double mean = tpulls[k].mean(tfitted);
        const double width = tpulls[k].width(tfitted);
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
  // the Xi- at the root, carried 2 cm along its helix: its decay vertex and
  // its momentum there
  const Truth xiAtRoot{.phi = 0.2, .theta = 1.3, .qOverP = -1. / 1.2};
  Propagator xiProp = prop;
  xiProp.massGeV = 1.32171F;
  const auto start = beamLinePlane(static_cast<Scalar>(kVertex[0]), static_cast<Scalar>(kVertex[1]), static_cast<Scalar>(xiAtRoot.phi));
  BoundTrackState xiState;
  xiState.par[kLoc0] = 0.F;
  xiState.par[kLoc1] = static_cast<Scalar>(kVertex[2]);
  xiState.par[kPhi] = static_cast<Scalar>(xiAtRoot.phi);
  xiState.par[kTheta] = static_cast<Scalar>(xiAtRoot.theta);
  xiState.par[kQOverP] = static_cast<Scalar>(xiAtRoot.qOverP);
  const Scalar onOrigin[3]{static_cast<Scalar>(kVertex[0] + (2. * std::cos(xiAtRoot.phi))), static_cast<Scalar>(kVertex[1] + (2. * std::sin(xiAtRoot.phi))), 0.F};
  const TrackingSurface on = radialPlane(onOrigin, static_cast<Scalar>(std::sin(xiAtRoot.phi)), static_cast<Scalar>(std::cos(xiAtRoot.phi)));
  SolenoidJacobian jacobian;
  BoundTrackState there = xiState;
  const auto status = jacobianAnalytic(xiProp, xiState, start, on, PropagationDirection::kForward, jacobian, &there);
  BOOST_REQUIRE(status == PropagationStatus::kOk || status == PropagationStatus::kOutOfBounds);
  Scalar xiDecay[3]{};
  surfaceToGlobal(on, there.loc0(), there.loc1(), xiDecay);
  const Vec3D second{xiDecay[0], xiDecay[1], xiDecay[2]};
  const auto xi = vectorOf(Truth{.phi = static_cast<double>(there.phi()), .theta = static_cast<double>(there.theta()), .qOverP = static_cast<double>(there.qOverP())});
  // its daughters at its decay: a pi-, the Lambda the rest
  const Truth bachelor{.phi = 0.5, .theta = 1.25, .qOverP = -1. / 0.3};
  const auto b = vectorOf(bachelor);
  const Vec3D lambda{xi[0] - b[0], xi[1] - b[1], xi[2] - b[2]};
  const double lambdaP = std::hypot(lambda[0], lambda[1], lambda[2]);
  const Vec3D third{second[0] + (3. * lambda[0] / lambdaP), second[1] + (3. * lambda[1] / lambdaP), second[2] + (3. * lambda[2] / lambdaP)};
  const Truth proton{.phi = 0.1, .theta = 1.2, .qOverP = 1. / 0.8};
  const auto pr = vectorOf(proton);
  const Truth lambdaPion = truthOf({lambda[0] - pr[0], lambda[1] - pr[1], lambda[2] - pr[2]}, -1.);
  const Truth piPlusA{.phi = -0.3, .theta = 1.5, .qOverP = 1. / 0.6};
  const Truth piPlusB{.phi = 0.7, .theta = 1.1, .qOverP = 1. / 0.4};
  const DecayTopology topology{.vertexOf = {0, 0, 1, 2, 2}, .parent = {-1, 0, 1}};
  const std::vector<DecayDaughter> exactStates{daughterOf(piPlusA, prop, kVertex), daughterOf(piPlusB, prop, kVertex), daughterOf(bachelor, prop, second), daughterOf(proton, prop, third), daughterOf(lambdaPion, prop, third)};
  const std::array<Vec3D, 3> truth{kVertex, second, third};

  const auto exact = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{exactStates}, UniformMaterial{}, {}, topology);
  BOOST_REQUIRE(exact.status == VertexFitStatus::kOk);
  for (int w = 0; w < 3; ++w) {
    for (int c = 0; c < 3; ++c) {
      BOOST_TEST(std::abs(exact.vertices[w][c] - truth[w][c]) < 1e-4); // float states
    }
  }
  std::mt19937 random{1234U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 1000;
  std::array<PullMoments, 9> pulls{};
  double chi2 = 0.;
  int fitted = 0;
  int ndf = 0;
  for (int trial = 0; trial < kTrials; ++trial) {
    auto states = smeared(exactStates, random, gauss);
    const auto fit = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{states}, UniformMaterial{}, {}, topology);
    if (fit.status != VertexFitStatus::kOk) {
      continue;
    }
    ++fitted;
    ndf = fit.ndf;
    for (int k = 0; k < 9; ++k) {
      const double pull = (fit.vertices[k / 3][k % 3] - truth[k / 3][k % 3]) / std::sqrt(fit.covarianceAt(k, k));
      pulls[k].add(pull);
    }
    chi2 += fit.chi2;
  }
  BOOST_TEST_MESSAGE("charged composite: fitted " << fitted << ", chi2 mean " << chi2 / fitted << " at ndf " << ndf);
  BOOST_TEST(fitted >= kTrials - 5);
  BOOST_TEST(std::abs((chi2 / fitted / ndf) - 1.) < 0.1);
  for (int k = 0; k < 9; ++k) {
    const double mean = pulls[k].mean(fitted);
    const double width = pulls[k].width(fitted);
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
  const std::array<Truth, 2> pair{Truth{.phi = 0.3, .theta = 1.2, .qOverP = 1. / 0.6}, Truth{.phi = 0.9, .theta = 1.35, .qOverP = -1. / 0.4}};
  std::vector<DecayDaughter> exactStates;
  exactStates.reserve(pair.size());
  for (const auto& truth : pair) {
    exactStates.push_back(daughterOf(truth, prop, kVertex));
  }
  constexpr Vec2D kWeight{0.75, 0.25};
  constexpr Vec2D kShift{0., 0.15}; ///< relative, |q/p| up
  constexpr Vec2D kWidth{0.01, 0.05};
  double energy = 0.;
  Vec3D total{};
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
  std::mt19937 const random{2468U};
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
      const std::span<const DecayDaughter> span{states};
      const auto fit = withDecayFitVariant(variant, [&](auto one) { return fitDecayVertex<TRKFT3, decltype(one)::value>(span, UniformMaterial{}, massOfPair); });
      if (fit.status != VertexFitStatus::kOk) {
        continue;
      }
      ++fitted;
      for (int i = 0; i < 2; ++i) {
        const int q = 3 + (3 * i) + 2;
        const double residual = fit.momenta[i][2] - pair[i].qOverP;
        pull2 += residual * residual / fit.covarianceAt(q, q);
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
  const Truth kaon{.phi = 0.35, .theta = 1.25, .qOverP = -1. / 1.1};
  const Truth positron{.phi = 0.05, .theta = 1.4, .qOverP = 1. / 0.7};
  const Vec3D neutrino{0.3, 0.25, 0.1};
  constexpr double kKaon = 0.493677;
  constexpr double kElectron = 0.000511;
  Vec3D total = neutrino;
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
  const Vec3D made{kVertex[0] - (0.1 * total[0] / pMother), kVertex[1] - (0.1 * total[1] / pMother), kVertex[2] - (0.1 * total[2] / pMother)};
  Propagator kaonProp = prop;
  kaonProp.massGeV = static_cast<Scalar>(kKaon);
  Propagator positronProp = prop;
  positronProp.massGeV = static_cast<Scalar>(kElectron);
  auto exactK = daughterOf(kaon, kaonProp, kVertex);
  exactK.mass = kKaon;
  auto exactE = daughterOf(positron, positronProp, kVertex);
  exactE.mass = kElectron;
  const std::vector<DecayDaughter> visible{exactK, exactE};
  const Vec2D visibleMasses{kKaon, kElectron};
  const DecayConstraints constraints{.mass = true, .massGeV = massMother, .production = true, .point = made};
  const auto fitFrom = [&](const std::vector<DecayDaughter>& states, const Vec3D& start) {
    auto all = states;
    DecayDaughter unseen{.mass = 0., .neutral = true, .unmeasured = true};
    startUnmeasured(unseen, start);
    all.push_back(unseen);
    return fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{all}, UniformMaterial{}, constraints);
  };
  // exact: one branch gives the neutrino back
  {
    const auto free = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{visible}, UniformMaterial{});
    BOOST_REQUIRE(free.status == VertexFitStatus::kOk);
    const auto starts = missingMomentumStarts(free, visibleMasses, made, massMother, 0.);
    BOOST_REQUIRE(!starts.empty());
    double best = std::numeric_limits<double>::infinity();
    for (const auto& start : starts) {
      const auto fit = fitFrom(visible, start);
      if (fit.status == VertexFitStatus::kOk) {
        const auto& q = fit.momenta[2];
        const double p = 1. / q[2];
        const Vec3D v{p * std::cos(q[0]) * std::sin(q[1]), p * std::sin(q[0]) * std::sin(q[1]), p * std::cos(q[1])};
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
  std::array<PullMoments, 3> pulls{};
  double chi2 = 0.;
  int fitted = 0;
  int ndf = 0;
  for (int trial = 0; trial < kTrials; ++trial) {
    auto states = smeared(visible, random, gauss);
    const auto free = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{states}, UniformMaterial{});
    if (free.status != VertexFitStatus::kOk) {
      continue;
    }
    DecayVertexFitResult nearest;
    double best = std::numeric_limits<double>::infinity();
    for (const auto& start : missingMomentumStarts(free, visibleMasses, made, massMother, 0.)) {
      const auto fit = fitFrom(states, start);
      if (fit.status != VertexFitStatus::kOk) {
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
    for (int c = 0; c < 3; ++c) {
      const double pull = (nearest.vertex[c] - kVertex[c]) / std::sqrt(nearest.covarianceAt(c, c));
      pulls[c].add(pull);
    }
    chi2 += nearest.chi2;
  }
  BOOST_TEST_MESSAGE("unseen neutrino: fitted " << fitted << ", chi2 mean " << chi2 / fitted << " at ndf " << ndf);
  BOOST_TEST(fitted >= kTrials - 50);
  BOOST_TEST(std::abs((chi2 / fitted / ndf) - 1.) < 0.15);
  for (int c = 0; c < 3; ++c) {
    const double mean = pulls[c].mean(fitted);
    const double width = pulls[c].width(fitted);
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
  const Truth kaon{.phi = 0.2, .theta = 1.3, .qOverP = -1. / 0.9};
  const double phiPi = 0.55;
  const double thetaPi = 1.2;
  const auto unit = [](double phi, double theta) { return Vec3D{std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta)}; };
  const auto k = unit(kaon.phi, kaon.theta);
  const double pk = 1. / std::abs(kaon.qOverP);
  const double ek = std::hypot(pk, kKaon);
  const auto u = unit(phiPi, thetaPi);
  const double cosine = dotOf(k, u);
  // M^2 = mK^2 + mpi^2 + 2 (ek epi - pk p cos): bisect p
  double low = 1e-3;
  double high = 10.;
  for (int i = 0; i < 200; ++i) {
    const double p = 0.5 * (low + high);
    const double m2 = (kKaon * kKaon) + (kPion * kPion) + (2. * ((ek * std::hypot(p, kPion)) - (pk * p * cosine)));
    (m2 < kKStar * kKStar ? low : high) = p;
  }
  const Truth pion{.phi = phiPi, .theta = thetaPi, .qOverP = 1. / (0.5 * (low + high))};
  const Truth other{.phi = -0.4, .theta = 1.5, .qOverP = 1. / 0.6};
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
    const auto exact = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{exactStates}, UniformMaterial{}, set);
    BOOST_REQUIRE(exact.status == VertexFitStatus::kOk);
    for (int i = 0; i < 3; ++i) {
      BOOST_TEST(std::abs(exact.vertex[i] - kVertex[i]) < 1e-4); // float states
    }
    std::array<PullMoments, 3> pulls{};
    int fitted = 0;
    for (int trial = 0; trial < kTrials; ++trial) {
      auto states = smeared(exactStates, random, gauss);
      const auto fit = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{states}, UniformMaterial{}, set);
      if (fit.status != VertexFitStatus::kOk) {
        continue;
      }
      ++fitted;
      for (int i = 0; i < 3; ++i) {
        const double pull = (fit.vertex[i] - kVertex[i]) / std::sqrt(fit.covarianceAt(i, i));
        pulls[i].add(pull);
      }
    }
    BOOST_TEST_MESSAGE("resonance: fitted " << fitted);
    BOOST_TEST(fitted >= kTrials - 5);
    for (int i = 0; i < 3; ++i) {
      const double mean = pulls[i].mean(fitted);
      const double width = pulls[i].width(fitted);
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
  const Truth helium{.phi = 0.3, .theta = 1.25, .qOverP = 2. / 2.4}; // q / p, |q| = 2
  const Truth pion{.phi = 0.7, .theta = 1.4, .qOverP = -1. / 0.25};
  double energy = 0.;
  Vec3D total{};
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
    const auto fit = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{exactStates}, UniformMaterial{});
    BOOST_REQUIRE(fit.status == VertexFitStatus::kOk);
    double e = 0.;
    Vec3D sum{};
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
  std::array<PullMoments, 3> pulls{};
  double chi2 = 0.;
  int fitted = 0;
  int ndf = 0;
  for (int trial = 0; trial < kTrials; ++trial) {
    auto states = smeared(exactStates, random, gauss);
    const auto fit = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{states}, UniformMaterial{}, set);
    if (fit.status != VertexFitStatus::kOk) {
      continue;
    }
    ++fitted;
    ndf = fit.ndf;
    for (int c = 0; c < 3; ++c) {
      const double pull = (fit.vertex[c] - kVertex[c]) / std::sqrt(fit.covarianceAt(c, c));
      pulls[c].add(pull);
    }
    chi2 += fit.chi2;
  }
  BOOST_TEST_MESSAGE("doubly charged under its mass: fitted " << fitted << ", chi2 mean " << chi2 / fitted << " at ndf " << ndf);
  BOOST_TEST(fitted >= kTrials - 5);
  BOOST_TEST(std::abs((chi2 / fitted / ndf) - 1.) < 0.1);
  for (int c = 0; c < 3; ++c) {
    const double mean = pulls[c].mean(fitted);
    const double width = pulls[c].width(fitted);
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
  daughters.reserve(kDaughters.size());
  for (const auto& truth : kDaughters) {
    daughters.push_back(daughterOf(truth, prop));
  }
  {
    const auto fit = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{daughters}, UniformMaterial{});
    BOOST_REQUIRE(fit.status == VertexFitStatus::kOk);
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
  exact.reserve(kDaughters.size());
  for (const auto& truth : kDaughters) {
    exact.push_back(daughterOf(truth, prop));
  }
  std::mt19937 random{12345U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 2000;
  std::array<PullMoments, 3> pulls{};
  double chi2 = 0.;
  int fitted = 0;
  for (int trial = 0; trial < kTrials; ++trial) {
    auto daughters = smeared(exact, random, gauss);
    const auto fit = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{daughters}, UniformMaterial{});
    if (fit.status != VertexFitStatus::kOk) {
      continue;
    }
    ++fitted;
    for (int k = 0; k < 3; ++k) {
      const double pull = (fit.vertex[k] - kVertex[k]) / std::sqrt(fit.covarianceAt(k, k));
      pulls[k].add(pull);
    }
    chi2 += fit.chi2;
  }
  BOOST_TEST(fitted == kTrials);
  for (int k = 0; k < 3; ++k) {
    const double mean = pulls[k].mean(fitted);
    const double width = pulls[k].width(fitted);
    BOOST_TEST_MESSAGE("vertex " << k << ": pull mean " << mean << ", width " << width);
    BOOST_TEST(std::abs(mean) < 0.1);
    BOOST_TEST(std::abs(width - 1.) < 0.05);
  }
  BOOST_TEST_MESSAGE("chi2 mean " << chi2 / fitted << " at ndf 3");
  BOOST_TEST(std::abs((chi2 / fitted) - 3.) < 0.2);
}

/// The primary vertex fit's block solve is the decay fit's dense one: on the
/// same two and three tracks, one vertex, no constraint, and a prior so wide
/// that it weighs nothing, from the same start, the same vertex, momenta and
/// chi2 to rounding.
BOOST_AUTO_TEST_CASE(APrimaryVertexIsTheDecayFitsVertex)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  for (const unsigned seed : {1U, 2U, 3U}) {
    for (const std::size_t count : {std::size_t{2}, std::size_t{3}}) {
      std::vector<DecayDaughter> exact;
      for (std::size_t i = 0; i < count; ++i) {
        exact.push_back(daughterOf(kDaughters[i], prop));
      }
      const auto daughters = smeared(exact, seed);
      const auto decay = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{daughters}, UniformMaterial{});
      BOOST_REQUIRE(decay.status == VertexFitStatus::kOk);
      std::vector<VertexTrack> tracks;
      for (const auto& daughter : daughters) {
        tracks.push_back({.innermost = daughter.innermost, .surface = daughter.surface, .propagator = daughter.propagator});
      }
      BeamPrior wide{.position = {0., 0., 0.}};
      for (int k = 0; k < 3; ++k) {
        wide.covariance[covIndex(k, k)] = 1e12;
      }
      std::vector<Vec3D> momenta(count);
      std::array<Vec3D, decay_fit_detail::kMaxConstrainedVertices> start{};
      BOOST_REQUIRE(decayFitStart<TRKFT3>(std::span<const DecayDaughter>{daughters}, {}, {}, {}, start));
      const auto primary = fitPrimaryVertex(std::span<const VertexTrack>{tracks}, UniformMaterial{}, wide, start[0], std::span<Vec3D>{momenta});
      BOOST_TEST_MESSAGE("primary against decay fit: iterations " << int(primary.iterations) << " / " << int(decay.iterations) << ", outer " << int(primary.outer) << " / " << int(decay.outer) << ", dz " << primary.vertex[2] - decay.vertex[2] << ", dchi2 " << (primary.chi2 - primary.priorChi2) - decay.chi2);
      BOOST_REQUIRE(primary.status == VertexFitStatus::kOk);
      BOOST_TEST(primary.ndf == decay.ndf + 3);
      for (int k = 0; k < 3; ++k) {
        BOOST_TEST(std::abs(primary.vertex[k] - decay.vertex[k]) < 1e-9);
        BOOST_TEST(std::abs(primary.covariance[covIndex(k, k)] - decay.covarianceAt(k, k)) < 1e-6 * decay.covarianceAt(k, k));
      }
      for (std::size_t i = 0; i < count; ++i) {
        for (int c = 0; c < 3; ++c) {
          BOOST_TEST(std::abs(momenta[i][c] - decay.momenta[i][c]) < 1e-9 * std::max(1., std::abs(decay.momenta[i][c])));
        }
      }
      BOOST_TEST(std::abs((primary.chi2 - primary.priorChi2) - decay.chi2) < 1e-6 * std::max(1., decay.chi2));
    }
  }
}

/// A primary vertex of three tracks under a beam prior about the true
/// vertex: unit pulls, the chi2's mean its ndf; the vertex without a track
/// is the fit without it, to the linearisation.
BOOST_AUTO_TEST_CASE(APrimaryVertexHasUnitPullsAndDropsATrack)
{
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  std::vector<VertexTrack> exact;
  for (const auto& truth : kDaughters) {
    const auto daughter = daughterOf(truth, prop);
    exact.push_back({.innermost = daughter.innermost, .surface = daughter.surface, .propagator = daughter.propagator});
  }
  constexpr std::array<double, 3> kWidths{0.002, 0.002, 0.05}; // [cm]
  std::mt19937 random{4321U};
  std::normal_distribution<double> gauss;
  constexpr int kTrials = 2000;
  std::array<PullMoments, 3> pulls{};
  double chi2 = 0.;
  int fitted = 0;
  double worstDrop = 0.;
  for (int trial = 0; trial < kTrials; ++trial) {
    auto tracks = exact;
    for (auto& track : tracks) {
      for (int k = 0; k < 5; ++k) {
        track.innermost.par[k] += static_cast<Scalar>(gauss(random) * std::sqrt(static_cast<double>(track.innermost.covAt(k, k))));
      }
    }
    BeamPrior prior;
    for (int k = 0; k < 3; ++k) {
      prior.position[k] = kVertex[k] + (kWidths[k] * gauss(random));
      prior.covariance[covIndex(k, k)] = kWidths[k] * kWidths[k];
    }
    PrimaryVertexFitter<UniformMaterial> fitter(std::span<const VertexTrack>{tracks}, UniformMaterial{}, prior);
    std::vector<std::max_align_t> storage((fitter.workspace() + sizeof(std::max_align_t) - 1U) / sizeof(std::max_align_t));
    PrimaryVertexFitResult fit;
    fitter.fit(prior.position, std::as_writable_bytes(std::span{storage}).first(fitter.workspace()), fit);
    if (fit.status != VertexFitStatus::kOk) {
      continue;
    }
    ++fitted;
    for (int k = 0; k < 3; ++k) {
      pulls[k].add((fit.vertex[k] - kVertex[k]) / std::sqrt(fit.covariance[covIndex(k, k)]));
    }
    chi2 += fit.chi2;
    if (trial < 50) {
      Vec3D dropped{};
      std::array<double, 6> droppedCovariance{};
      BOOST_REQUIRE(fitter.withoutTrack(2, dropped, droppedCovariance));
      const auto refit = fitPrimaryVertex(std::span<const VertexTrack>{tracks.data(), 2}, UniformMaterial{}, prior, prior.position);
      BOOST_REQUIRE(refit.status == VertexFitStatus::kOk);
      for (int k = 0; k < 3; ++k) {
        worstDrop = std::max(worstDrop, std::abs(dropped[k] - refit.vertex[k]) / std::sqrt(refit.covariance[covIndex(k, k)]));
      }
    }
  }
  BOOST_TEST(fitted == kTrials);
  for (int k = 0; k < 3; ++k) {
    BOOST_TEST_MESSAGE("primary vertex " << k << ": pull mean " << pulls[k].mean(fitted) << ", width " << pulls[k].width(fitted));
    BOOST_TEST(std::abs(pulls[k].mean(fitted)) < 0.1);
    BOOST_TEST(std::abs(pulls[k].width(fitted) - 1.) < 0.05);
  }
  BOOST_TEST_MESSAGE("primary chi2 mean " << chi2 / fitted << " at ndf 6; the dropped track's vertex within " << worstDrop << " sigma of the refit");
  BOOST_TEST(std::abs((chi2 / fitted) - 6.) < 0.3);
  BOOST_TEST(worstDrop < 0.05);
}

BOOST_AUTO_TEST_CASE(TheMixtureCollapseIsTheTwoPassOne)
{
  // Components of 3 parameters with log-weights spread over 700 (rising and
  // falling maxima; exp of the spread underflows), streamed, merged from
  // halves, against the normalised two-pass collapse.
  using decay_fit_detail::MixtureCollapse;
  constexpr int n = 3;
  constexpr std::size_t kEntries = static_cast<std::size_t>(n) * n; // a full n x n matrix
  std::mt19937 random(7U);
  std::normal_distribution<double> gauss;
  struct Component {
    double logWeight;
    std::vector<double> p;
    std::vector<double> sigma;
    double chi2;
  };
  std::vector<Component> mixture;
  for (const double logWeight : {-3., 1.5, -700., 2., 0.5, -1., 2.5, -20.}) {
    Component c{.logWeight = logWeight, .p = {gauss(random), gauss(random), gauss(random)}, .sigma = std::vector<double>(kEntries, 0.), .chi2 = 3. + gauss(random)};
    std::array<double, kEntries> b{};
    for (auto& v : b) {
      v = gauss(random);
    }
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
          c.sigma[(i * n) + j] += b[(i * n) + k] * b[(j * n) + k];
        }
      }
    }
    mixture.push_back(c);
  }
  double largest = -std::numeric_limits<double>::infinity();
  for (const auto& c : mixture) {
    largest = std::max(largest, c.logWeight);
  }
  double total = 0.;
  for (const auto& c : mixture) {
    total += std::exp(c.logWeight - largest);
  }
  std::vector<double> mean(n, 0.);
  double chi2 = 0.;
  for (const auto& c : mixture) {
    const double w = std::exp(c.logWeight - largest) / total;
    for (int i = 0; i < n; ++i) {
      mean[i] += w * c.p[i];
    }
    chi2 += w * c.chi2;
  }
  std::vector<double> covariance(kEntries, 0.);
  for (const auto& c : mixture) {
    const double w = std::exp(c.logWeight - largest) / total;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        covariance[(i * n) + j] += w * (c.sigma[(i * n) + j] + ((c.p[i] - mean[i]) * (c.p[j] - mean[j])));
      }
    }
  }

  // each collapse's mean and scatter (packed) in storage of its own; the
  // components' covariances packed as the collapse takes them
  std::vector<std::vector<double>> storage(12, std::vector<double>(covSize(n), 0.));
  for (auto& c : mixture) {
    std::vector<double> packedSigma(covSize(n));
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j <= i; ++j) {
        packedSigma[covIndex(i, j)] = c.sigma[(i * n) + j];
      }
    }
    c.sigma = packedSigma;
  }
  const auto collapseIn = [&](int at) { return MixtureCollapse(std::span<double>{storage[2U * static_cast<std::size_t>(at)].data(), n}, std::span<double>{storage[(2U * static_cast<std::size_t>(at)) + 1]}); };
  auto stream = collapseIn(0);
  auto first = collapseIn(1);
  auto second = collapseIn(2);
  for (std::size_t c = 0; c < mixture.size(); ++c) {
    stream.add(mixture[c].logWeight, mixture[c].p, mixture[c].sigma, mixture[c].chi2);
    (c % 2 == 0 ? first : second).add(mixture[c].logWeight, mixture[c].p, mixture[c].sigma, mixture[c].chi2);
  }
  auto merged = collapseIn(3);
  merged.merge(collapseIn(4)); // an empty one is the identity
  merged.merge(second);
  merged.merge(first);
  for (const auto* one : {&stream, &merged}) {
    BOOST_TEST(one->count == static_cast<int>(mixture.size()));
    BOOST_TEST(std::abs(one->meanChi2() - chi2) < 1e-12);
    std::vector<double> c(covSize(n));
    one->covarianceInto(c.data());
    for (int i = 0; i < n; ++i) {
      BOOST_TEST(std::abs(one->mean[i] - mean[i]) < 1e-12);
      for (int j = 0; j < n; ++j) {
        BOOST_TEST(std::abs(c[covIndex(i, j)] - covariance[(i * n) + j]) < 1e-12);
      }
    }
  }
  auto empty = collapseIn(5);
  empty.merge(collapseIn(4));
  BOOST_TEST(empty.count == 0);
}

BOOST_AUTO_TEST_CASE(InputsOutsideTheCapacitiesAreRefused)
{
  // Two daughters of a two-vertex tree whose second vertex hangs from
  // itself, a resonance of a third daughter, a tracked mother at a vertex
  // that is not there, one with more rows than a track holds.
  const std::vector<DecayDaughter> two(2);
  const DecayTopology loop{.vertexOf = {0, 1}, .parent = {-1, 1}};
  BOOST_TEST(static_cast<int>(fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{two}, UniformMaterial{}, {}, loop).status) == static_cast<int>(VertexFitStatus::kInvalidInput));
  DecayConstraints resonance;
  resonance.resonances.push_back({.daughters = {0, 2}, .massGeV = 1., .widthGeV = 0.1});
  BOOST_TEST(static_cast<int>(fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{two}, UniformMaterial{}, resonance).status) == static_cast<int>(VertexFitStatus::kInvalidInput));
  const DecayMother elsewhere{.vertex = 1};
  BOOST_TEST(static_cast<int>(fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{two}, UniformMaterial{}, {}, {}, std::span<const DecayMother>{&elsewhere, 1}).status) == static_cast<int>(VertexFitStatus::kInvalidInput));
  const std::vector<TrackingSurface> planes(maxChainLayers<TRKFT3> + 1U);
  const std::vector<Vec2D> positions(planes.size());
  const std::vector<Vec3D> covariance(planes.size());
  const DecayMother longer{.planes = planes, .positions = positions, .covariance = covariance};
  BOOST_TEST(static_cast<int>(fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{two}, UniformMaterial{}, {}, {}, std::span<const DecayMother>{&longer, 1}).status) == static_cast<int>(VertexFitStatus::kTooLarge));
  const std::vector<DecayDaughter> six(6);
  BOOST_TEST(static_cast<int>(fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{six}, UniformMaterial{}).status) == static_cast<int>(VertexFitStatus::kTooLarge));
}

BOOST_AUTO_TEST_CASE(InputsBuiltBeyondTheirCapacityOrInconsistentAreRefused)
{
  const std::vector<DecayDaughter> two(2);
  const auto statusOf = [](const DecayVertexFitResult& fit) { return static_cast<int>(fit.status); };
  // A fixed-capacity input built with more than it holds keeps what fits and
  // says so; the fit refuses it as too large.
  const FixedVector<double, 3> four{1., 2., 3., 4.};
  BOOST_TEST(four.size() == 3U);
  BOOST_TEST(four.overflowed());
  BOOST_TEST(four[2] == 3.);
  const DecayTopology crowded{.vertexOf = {0, 0, 0, 0, 0, 0}};
  BOOST_TEST(statusOf(fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{two}, UniformMaterial{}, {}, crowded)) == static_cast<int>(VertexFitStatus::kTooLarge));
  // A mixture that claims more components than it holds.
  auto mixed = two;
  mixed[0].mixture.size = std::remove_cvref_t<decltype(mixed[0].mixture)>::kCapacity + 1;
  BOOST_TEST(statusOf(fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{mixed}, UniformMaterial{})) == static_cast<int>(VertexFitStatus::kInvalidInput));
  // A composite's mass where the tree has no composite: no equation, no
  // degree of freedom either.
  DecayConstraints orphan;
  orphan.compositeMasses.push_back(1.);
  BOOST_TEST(statusOf(fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{two}, UniformMaterial{}, orphan)) == static_cast<int>(VertexFitStatus::kInvalidInput));
  // A tracked mother with fewer positions than planes, through a batch: the
  // batch keeps each array's length, so the fit's check sees the mismatch.
  const std::vector<TrackingSurface> planes(2);
  const std::vector<Vec2D> positions(1);
  const std::vector<Vec3D> covariance(2);
  const DecayMother ragged{.planes = planes, .positions = positions, .covariance = covariance};
  DecayFitBatch batch;
  batch.add(std::span<const DecayDaughter>{two}, {}, {}, std::span<const DecayMother>{&ragged, 1}, {});
  std::vector<DecayVertexFitResult> results(1);
  fitDecayBatch<TRKFT3, DecayFitVariant::kFull>(batch.view(), UniformMaterial{}, std::span<DecayVertexFitResult>{results});
  BOOST_TEST(statusOf(results[0]) == static_cast<int>(VertexFitStatus::kInvalidInput));
}

BOOST_AUTO_TEST_CASE(AnAzimuthAnyTurnsAwayIsReduced)
{
  // Newton's steps are not bounded in phi: an azimuth whole turns away
  // reduces as std::remainder does.
  for (const int turns : {-7, -5, -1, 0, 1, 5, 7}) {
    for (const double offset : {-3., -0.3, 0., 0.3, 3.}) {
      const double x = offset + (kTwoPi<double> * turns);
      BOOST_TEST(std::abs(toPMPiOf(x) - std::remainder(x, kTwoPi<double>)) < 1e-12);
    }
  }
}

BOOST_AUTO_TEST_CASE(AComponentWhoseWeightUnderflowsLeavesNoTrace)
{
  // (log w = -1000, mean 1e16) then (log w = 0, mean 1): the first one's
  // relative weight is zero in double, so the mean is 1.
  using decay_fit_detail::MixtureCollapse;
  std::array<double, 1> mean{};
  std::array<double, 1> scatter{};
  MixtureCollapse<double> collapse(std::span<double>{mean}, std::span<double>{scatter});
  const std::array<double, 1> far{1e16};
  const std::array<double, 1> near{1.};
  const std::array<double, 1> sigma{1.};
  collapse.add(-1000., far, sigma, 0.);
  collapse.add(0., near, sigma, 0.);
  BOOST_TEST(collapse.mean[0] == 1.);
  std::array<double, 1> c{};
  collapse.covarianceInto(c.data());
  BOOST_TEST(c[0] == 1.);
}

BOOST_AUTO_TEST_CASE(TheFitsCrossingIsThePropagatorsOne)
{
  // The fit's crossing (`crossToPlane`, double, on dual numbers) against the
  // propagator's (`jacobianAnalytic`) on the planes both take, a plane along
  // the beam and one across it; its derivatives against central differences
  // of itself.
  using namespace decay_fit_detail;
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  prop.massGeV = static_cast<Scalar>(kPion);
  const double fieldScale = 20. * static_cast<double>(constants::math::B2C);
  const auto set = [](Scalar* v, auto x, auto y, auto z) {
    v[0] = static_cast<Scalar>(x);
    v[1] = static_cast<Scalar>(y);
    v[2] = static_cast<Scalar>(z);
  };
  const auto crossing = [&](const std::array<double, 6>& start, const TrackingSurface& plane, std::array<Dual<6>, 3>& out) {
    FreeState<6, double> from{};
    for (int c = 0; c < 6; ++c) {
      from[c] = dual::variable<6>(start[c], c);
    }
    FreeState<6, double> to{};
    std::array<Dual<6>, 2> local{};
    const bool reached = crossToPlane(from, plane, fieldScale, to, local);
    out = {local[0], local[1], to[3]};
    return reached;
  };
  for (const Truth& truth : {Truth{.phi = 0.3, .theta = 1.2, .qOverP = 1. / 0.3}, Truth{.phi = -0.4, .theta = 1.9, .qOverP = -1.}, Truth{.phi = 2.5, .theta = 0.6, .qOverP = 0.2}}) {
    const Vec3D vertex{0.2, -0.1, 0.4};
    const auto start = beamLinePlane(static_cast<Scalar>(vertex[0]), static_cast<Scalar>(vertex[1]), static_cast<Scalar>(truth.phi));
    BoundTrackState state;
    state.par[kLoc1] = static_cast<Scalar>(vertex[2]);
    state.par[kPhi] = static_cast<Scalar>(truth.phi);
    state.par[kTheta] = static_cast<Scalar>(truth.theta);
    state.par[kQOverP] = static_cast<Scalar>(truth.qOverP);
    const std::array<double, 6> free{vertex[0], vertex[1], vertex[2], static_cast<double>(state.phi()), static_cast<double>(state.theta()), static_cast<double>(state.qOverP())};
    // along the beam, 2 cm out along the direction; across it, 3 cm downstream
    auto along = beamLinePlane(static_cast<Scalar>(vertex[0] + (2. * std::cos(truth.phi))), static_cast<Scalar>(vertex[1] + (2. * std::sin(truth.phi))), static_cast<Scalar>(truth.phi));
    along.origin[2] = static_cast<Scalar>(vertex[2]);
    TrackingSurface across;
    const double side = std::cos(truth.theta) > 0. ? 1. : -1.;
    set(across.origin, 0.F, 0.F, static_cast<Scalar>(vertex[2] + (3. * side)));
    set(across.axis0, 1.F, 0.F, 0.F);
    set(across.axis1, 0.F, static_cast<Scalar>(side), 0.F);
    set(across.normal, 0.F, 0.F, static_cast<Scalar>(side));
    for (const auto* plane : {&along, &across}) {
      SolenoidJacobian jacobian;
      BoundTrackState reference = state;
      const auto status = jacobianAnalytic(prop, state, start, *plane, PropagationDirection::kForward, jacobian, &reference);
      BOOST_REQUIRE(status == PropagationStatus::kOk || status == PropagationStatus::kOutOfBounds);
      std::array<Dual<6>, 3> fit{};
      BOOST_REQUIRE(crossing(free, *plane, fit));
      // the propagator's state is float: its rounding, relative to the
      // lever arm of a few cm
      const double tolerance = 64. * std::numeric_limits<Scalar>::epsilon() * 10.;
      BOOST_TEST(std::abs(fit[0].v - static_cast<double>(reference.loc0())) < tolerance);
      BOOST_TEST(std::abs(fit[1].v - static_cast<double>(reference.loc1())) < tolerance);
      BOOST_TEST(std::abs(toPMPiOf(fit[2].v - static_cast<double>(reference.phi()))) < tolerance);
      for (int k = 0; k < 6; ++k) {
        const double h = 1e-6 * std::max(std::abs(free[k]), 1.);
        auto up = free;
        auto down = free;
        up[k] += h;
        down[k] -= h;
        std::array<Dual<6>, 3> above{};
        std::array<Dual<6>, 3> below{};
        BOOST_REQUIRE(crossing(up, *plane, above));
        BOOST_REQUIRE(crossing(down, *plane, below));
        for (int r = 0; r < 3; ++r) {
          const double difference = (above[r].v - below[r].v) / (2. * h);
          BOOST_TEST(std::abs(fit[r].d[k] - difference) < 1e-6 * std::max(std::abs(difference), 1.));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(ReducedPrecisionOnHardCases)
{
  // A boosted, nearly collinear V0: two pions of 25 GeV, 1 mrad apart, with
  // and without its mass.
  Propagator prop;
  prop.bzKiloGauss = 20.F;
  prop.massGeV = static_cast<Scalar>(kPion);
  for (const double opening : {1e-3, 1e-2}) {
    constexpr double kP = 25.;
    const double half = 0.5 * opening;
    const Truth plus{.phi = 0.4 + half, .theta = 1.3, .qOverP = 1. / kP};
    const Truth minus{.phi = 0.4 - half, .theta = 1.3, .qOverP = -1. / kP};
    const auto pair = smeared({daughterOf(plus, prop, kDecay), daughterOf(minus, prop, kDecay)}, 11U);
    const double e = std::hypot(kP, kPion);
    const double mass = std::sqrt((2. * kPion * kPion) + (2. * ((e * e) - (kP * kP * std::cos(opening)))));
    const auto name = std::string("boosted V0, opening ") + std::to_string(opening);
    reportPrecision(name.c_str(), pair);
    reportPrecision((name + ", under its mass").c_str(), pair, DecayConstraints{.mass = true, .massGeV = mass});
  }
  // A resonance: the pion pair's mass as a Breit-Wigner of 50 MeV.
  {
    const auto pions = pionPair();
    const auto pair = smeared({daughterOf(pions[0], prop, kDecay), daughterOf(pions[1], prop, kDecay)}, 12U);
    DecayConstraints resonance;
    resonance.resonances.push_back({.daughters = {0, 1}, .massGeV = kMass, .widthGeV = 0.05});
    reportPrecision("resonance", pair, resonance);
  }
  // An unseen neutral: a K0S's two pions with a massless third, its
  // momentum fitted from the mother's mass and production point.
  {
    const auto pions = pionPair();
    const auto visible = smeared({daughterOf(pions[0], prop, kDecay), daughterOf(pions[1], prop, kDecay)}, 13U);
    const Vec2D masses{kPion, kPion};
    const auto free = fitDecayVertex<TRKFT3, DecayFitVariant::kFull>(std::span<const DecayDaughter>{visible}, UniformMaterial{});
    BOOST_REQUIRE(free.status == VertexFitStatus::kOk);
    // a mother heavier than the pair: room for the unseen one
    const auto starts = missingMomentumStarts(free, masses, kPoint, kMass + 0.1, 0.);
    BOOST_REQUIRE(!starts.empty());
    auto all = visible;
    DecayDaughter unseen{.mass = 0., .neutral = true, .unmeasured = true};
    startUnmeasured(unseen, starts.front());
    all.push_back(unseen);
    reportPrecision("missing momentum", all, DecayConstraints{.mass = true, .massGeV = kMass + 0.1, .production = true, .point = kPoint});
  }
  // A crossing near the quarter turn: a helix of radius R from the origin
  // along x to the plane x = d along the beam, d / R towards 1, in double and
  // in float: whether both accept it, and the float crossing's offset.
  {
    using namespace decay_fit_detail;
    TrackingSurface plane;
    plane.normal[0] = 1.F;
    plane.axis0[1] = 1.F;
    plane.axis1[2] = 1.F;
    constexpr double kD = 3.;
    plane.origin[0] = static_cast<Scalar>(kD);
    const double fieldScale = 20. * static_cast<double>(constants::math::B2C);
    for (const double ratio : {0.9, 0.99, 0.999, 0.9999, 0.99999}) {
      const double radius = kD / ratio;
      const double qOverPt = 1. / (fieldScale * radius); // rho = fieldScale q/pT
      const auto crossing = [&]<typename T>(Vec2D& at) {
        FreeState<6, T> from{};
        for (int c = 0; c < 6; ++c) {
          from[c].d[c] = T(1.);
        }
        from[3].v = T(0.);
        from[4].v = static_cast<T>(kPi<double> / 2.);
        from[5].v = static_cast<T>(qOverPt);
        FreeState<6, T> to{};
        std::array<Dual<6, T>, 2> local{};
        const bool reached = crossToPlane(from, plane, fieldScale, to, local);
        at = {static_cast<double>(local[0].v), static_cast<double>(local[1].v)};
        return reached;
      };
      Vec2D inDouble{};
      Vec2D inFloat{};
      const bool reachedDouble = crossing.template operator()<double>(inDouble);
      const bool reachedFloat = crossing.template operator()<float>(inFloat);
      BOOST_TEST_MESSAGE("quarter turn, d / R " << ratio << ": reached in double " << reachedDouble << ", in float " << reachedFloat << ", loc0 " << inDouble[0] << ", float less double [um] " << 1e4 * (inFloat[0] - inDouble[0]));
    }
  }
}

} // namespace o2::gtrk::study
