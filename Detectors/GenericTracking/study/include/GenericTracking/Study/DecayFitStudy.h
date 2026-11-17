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

/// \file DecayFitStudy.h
/// \brief Single-vertex decays fitted from their daughters' chains
///        (DecayVertexFit.h) and held to the MC.
///
/// The decays: members of the decay trees with a vertex of their own (the
/// products made away from the member's production point), whose products are
/// all final charged particles (no composite below: one vertex), two to
/// `kMaxDaughters` of them -- a photon's conversion pair among them --, every
/// one with a chain from the source that is findable: its rows -- the chain
/// the fit takes, a truth chain's first outward leg -- on a chain of at least
/// four linked layers (`chainLayers`, `isFindable`, no pT floor): a track of
/// three is barely defined. An electron daughter
/// is refitted with the Gaussian-sum filter and carried as its mixture
/// (DecayVertexFit.h); the rows' fit takes its radiation as a Gaussian. Each daughter's rows are refitted
/// with its true species (`refitWithHypothesis`); the vertex is fitted under
/// every `DecayFitVariant`, a record each: from the daughters' states, and
/// from their rows (`kRows`, DecayRowsFit.h, the information limit); and
/// both fits (kFull, kRows) under the mother's mass, its production at the
/// MC collision vertex, and both (DecayConstraints.h); a photon's under
/// collinear daughters in place of the mass. Two nearly parallel daughters
/// measure a conversion's point along the flight only through their
/// curvature, weakly and far from linearly, so every photon fit is profiled
/// along the line from the collision towards the pair's first rows over the
/// material table's cells on it (`conversionCells`, `profileAlongLine`).
///
/// Truth: the vertex against the MC decay point (the products' production
/// point), in the global frame and along the true flight; the daughters'
/// momenta against their MC momenta at production; the invariant mass of the
/// refitted momenta and of the carried ones (without the vertex's
/// information) against the true daughters' (with the species' masses), the
/// summed momentum against theirs. The decay length L is the vertex's
/// distance from the collision vertex along the mother's direction, ct = L M /
/// p with the mother's PDG mass, the truth's from the true decay point and
/// the true daughters' momenta; the mother's direction is the summed
/// momentum's.

#ifndef GTRK_STUDY_DECAYFITSTUDY_H_
#define GTRK_STUDY_DECAYFITSTUDY_H_

#include "GenericTracking/Study/DecayCandidates.h"
#include "GenericTracking/Study/DecayConstraints.h"
#include "GenericTracking/Study/DecayProfile.h"
#include "GenericTracking/Study/DecayRowsFit.h"
#include "GenericTracking/Study/TrackComparison.h"
#include "GenericTracking/Study/DecayTrees.h"
#include "GenericTracking/Study/DecayVertexFit.h"
#include "GenericTracking/Study/RecordWriter.h"
#include "GenericTracking/Study/StateTruth.h"
#include "GenericTracking/Study/StudyChains.h"
#include "GenericTracking/Study/TruthEvent.h"
#include "GenericTracking/Study/TruthFit.h"
#include "GenericTracking/Tracking/HypothesisRefit.h"
#include "GenericTracking/Tracking/Parallel.h"
#include "GenericTracking/Tracking/RoadFit.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <format>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace o2::gtrk::study
{

/// One decay fitted under one variant (`DecayFitStudy`).
struct DecayFitRecord {
  static constexpr std::string_view kName = "decayFit";
  static constexpr std::size_t kDaughterColumns = 3; ///< daughters with their own columns
  static constexpr std::uint8_t kRefitFailed = 254U; ///< `status`: a daughter's refit failed

  std::uint32_t timeframe = 0U;
  std::uint8_t source = 0U; ///< ChainSource
  std::uint64_t particle = 0U;
  std::int32_t pdg = 0;
  std::uint8_t variant = 0U;     ///< DecayFitVariant
  std::uint8_t constraints = 0U; ///< 1 the mass, 2 the production point, 4 collinear daughters, 8 a conversion's profile along its line, 16 a tree's composites' masses, 32 the tracked mother's own rows
  std::uint8_t daughters = 0U;
  std::uint8_t status = 0U; ///< DecayFitStatus, or kRefitFailed
  std::uint8_t iterations = 0U;
  std::uint8_t outer = 0U;
  std::uint8_t stopReason = 0U;
  float stopMetric = 0.F;
  float modelShift = 0.F;
  float microseconds = 0.F; ///< the fit's time (with a profile's or a sequential fit's every fit)
  float chi2 = 0.F;
  std::int16_t ndf = 0;
  float flight = 0.F; ///< the true one [cm]
  /// The mother's production point's distance from its collision's vertex
  /// [cm]: 0 for a mother made at the collision, NaN without the collision.
  float productionOffset = std::numeric_limits<float>::quiet_NaN();
  float decayR = 0.F;        ///< the true decay point's distance from the z axis [cm]
  float pTrue = 0.F;         ///< the true daughters' summed |p| [GeV]
  float minDaughterPt = 0.F; ///< the lowest true daughter pT [GeV]
  float maxCarryX2X0 = 0.F;  ///< the most material a daughter's last carry crossed
  std::uint8_t minDaughterRows = 0U;
  std::uint8_t motherRows = 0U;  ///< a charged root's own rows (a tracked Xi or Omega), 0 without
  float pointSigma = 0.F;        ///< the production point's scattering across the flight [cm] the fit took (0: held exactly)
  std::uint8_t refitPass = 0U;   ///< a failed daughter refit's FitPass
  std::uint8_t refitStatus = 0U; ///< and its PropagationStatus
  float dX = 0.F;                ///< fitted less true vertex [cm]
  float dY = 0.F;
  float dZ = 0.F;
  float sX = 0.F; ///< the fit's sigmas [cm]
  float sY = 0.F;
  float sZ = 0.F;
  float dL = 0.F; ///< along the true flight [cm]
  float sL = 0.F;
  float dT = 0.F; ///< across it, in the transverse plane [cm]
  float sT = 0.F;
  float dN = 0.F; ///< across it, the third axis [cm]
  float sN = 0.F;
  float mass = 0.F; ///< of the refitted momenta [GeV]
  float sMass = 0.F;
  float massCarried = 0.F; ///< of the last carried states' momenta [GeV]
  float sMassCarried = 0.F;
  float massTrue = 0.F;
  float p = 0.F; ///< the refitted momenta's summed |p| [GeV]
  float sP = 0.F;
  float decayLength = 0.F; ///< L [cm]
  float sDecayLength = 0.F;
  float decayLengthTrue = 0.F;
  float ct = 0.F; ///< [cm]
  float sCt = 0.F;
  float ctTrue = 0.F;
  /// A conversion's profile along its line: the true point's distance along
  /// it less the profile's mean, over its spread, and the spread [cm]; NaN
  /// for another decay.
  float profilePull = std::numeric_limits<float>::quiet_NaN();
  float profileSigma = std::numeric_limits<float>::quiet_NaN();
  /// A tree's first composite below the root (a cascade's Lambda): its
  /// species, its vertex against its true decay point along and across its
  /// true flight, its decay length from its parent's vertex along its
  /// momentum, its mass from its daughters; NaN for a one-vertex decay.
  std::int32_t childPdg = 0;
  float dChildL = std::numeric_limits<float>::quiet_NaN();
  float sChildL = std::numeric_limits<float>::quiet_NaN();
  float dChildT = std::numeric_limits<float>::quiet_NaN();
  float sChildT = std::numeric_limits<float>::quiet_NaN();
  float childLength = std::numeric_limits<float>::quiet_NaN();
  float sChildLength = std::numeric_limits<float>::quiet_NaN();
  float childLengthTrue = std::numeric_limits<float>::quiet_NaN();
  float childMass = std::numeric_limits<float>::quiet_NaN();
  float sChildMass = std::numeric_limits<float>::quiet_NaN();
  float childMassTrue = std::numeric_limits<float>::quiet_NaN();
  float dPhiMother = 0.F; ///< the mother's direction, fitted less true
  float dThetaMother = 0.F;
  float sPhiMother = 0.F;
  float sThetaMother = 0.F;
  /// The first daughters' refitted phi, theta, q/p less the true ones, and
  /// their sigmas.
  std::array<float, kDaughterColumns> dPhi{};
  std::array<float, kDaughterColumns> dTheta{};
  std::array<float, kDaughterColumns> dQOverP{};
  std::array<float, kDaughterColumns> sPhi{};
  std::array<float, kDaughterColumns> sTheta{};
  std::array<float, kDaughterColumns> sQOverP{};
  std::array<float, kDaughterColumns> dPRel{}; ///< fitted |p| over the true one, less 1
  /// The energy loss from the decay point to the first row, as fractions of
  /// |p|: the MC's (|p| entering the first row's sensor over |p| at
  /// production) and the model's over the straight segment between (the LUT's
  /// material, `energyLoss`); and the daughter's Kalman fit's |p| on its
  /// first row over the MC's there, less 1. NaN without the MC hits.
  std::array<float, kDaughterColumns> lossTrue{};
  std::array<float, kDaughterColumns> lossModel{};
  std::array<float, kDaughterColumns> dPRelInner{};
  /// The daughter's MC hits (sensitive-volume crossings) before its first
  /// row's: where it crossed a sensor and left no row of its chain; -1
  /// without the MC hits.
  std::array<std::int16_t, kDaughterColumns> hitsBefore{-1, -1, -1};

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("source", source);
    visit("particle", particle);
    visit("pdg", pdg);
    visit("variant", variant);
    visit("constraints", constraints);
    visit("daughters", daughters);
    visit("status", status);
    visit("iterations", iterations);
    visit("outer", outer);
    visit("modelShift", modelShift);
    visit("microseconds", microseconds);
    visit("stopReason", stopReason);
    visit("stopMetric", stopMetric);
    visit("chi2", chi2);
    visit("ndf", ndf);
    visit("flight", flight);
    visit("productionOffset", productionOffset);
    visit("decayR", decayR);
    visit("pTrue", pTrue);
    visit("minDaughterPt", minDaughterPt);
    visit("maxCarryX2X0", maxCarryX2X0);
    visit("minDaughterRows", minDaughterRows);
    visit("motherRows", motherRows);
    visit("pointSigma", pointSigma);
    visit("refitPass", refitPass);
    visit("refitStatus", refitStatus);
    visit("dX", dX);
    visit("dY", dY);
    visit("dZ", dZ);
    visit("sX", sX);
    visit("sY", sY);
    visit("sZ", sZ);
    visit("dL", dL);
    visit("sL", sL);
    visit("dT", dT);
    visit("sT", sT);
    visit("dN", dN);
    visit("sN", sN);
    visit("mass", mass);
    visit("sMass", sMass);
    visit("massCarried", massCarried);
    visit("sMassCarried", sMassCarried);
    visit("massTrue", massTrue);
    visit("p", p);
    visit("sP", sP);
    visit("decayLength", decayLength);
    visit("sDecayLength", sDecayLength);
    visit("decayLengthTrue", decayLengthTrue);
    visit("ct", ct);
    visit("sCt", sCt);
    visit("ctTrue", ctTrue);
    visit("profilePull", profilePull);
    visit("profileSigma", profileSigma);
    visit("childPdg", childPdg);
    visit("dChildL", dChildL);
    visit("sChildL", sChildL);
    visit("dChildT", dChildT);
    visit("sChildT", sChildT);
    visit("childLength", childLength);
    visit("sChildLength", sChildLength);
    visit("childLengthTrue", childLengthTrue);
    visit("childMass", childMass);
    visit("sChildMass", sChildMass);
    visit("childMassTrue", childMassTrue);
    visit("dPhiMother", dPhiMother);
    visit("dThetaMother", dThetaMother);
    visit("sPhiMother", sPhiMother);
    visit("sThetaMother", sThetaMother);
    static constexpr std::array<std::string_view, kDaughterColumns> kSuffix{"0", "1", "2"};
    for (std::size_t d = 0; d < kDaughterColumns; ++d) {
      // The names must outlive the writer's first call: kept in static storage.
      static const std::array<std::array<std::string, 11>, kDaughterColumns> kNames = [] {
        std::array<std::array<std::string, 11>, kDaughterColumns> names;
        for (std::size_t k = 0; k < kDaughterColumns; ++k) {
          const std::string s{kSuffix[k]};
          names[k] = {"dPhi" + s, "dTheta" + s, "dQOverP" + s, "sPhi" + s, "sTheta" + s, "sQOverP" + s, "dPRel" + s, "lossTrue" + s, "lossModel" + s, "dPRelInner" + s, "hitsBefore" + s};
        }
        return names;
      }();
      visit(kNames[d][0], dPhi[d]);
      visit(kNames[d][1], dTheta[d]);
      visit(kNames[d][2], dQOverP[d]);
      visit(kNames[d][3], sPhi[d]);
      visit(kNames[d][4], sTheta[d]);
      visit(kNames[d][5], sQOverP[d]);
      visit(kNames[d][6], dPRel[d]);
      visit(kNames[d][7], lossTrue[d]);
      visit(kNames[d][8], lossModel[d]);
      visit(kNames[d][9], dPRelInner[d]);
      visit(kNames[d][10], hitsBefore[d]);
    }
  }
};

namespace decay_fit_study_detail
{

/// A daughter's momentum vector from phi, theta, q/p.
[[nodiscard]] inline std::array<double, 3> momentumOf(const std::array<double, 3>& q)
{
  const double p = 1. / std::abs(q[2]);
  return {p * std::cos(q[0]) * std::sin(q[1]), p * std::sin(q[0]) * std::sin(q[1]), p * std::cos(q[1])};
}

/// d(px, py, pz) / d(phi, theta, q/p), row-major 3 x 3.
[[nodiscard]] inline std::array<double, 9> momentumJacobian(const std::array<double, 3>& q)
{
  const double p = 1. / std::abs(q[2]);
  const double dpdk = -p / q[2]; // d|p| / d(q/p)
  const double cp = std::cos(q[0]);
  const double sp = std::sin(q[0]);
  const double ct = std::cos(q[1]);
  const double st = std::sin(q[1]);
  return {-p * sp * st, p * cp * ct, dpdk * cp * st,
          p * cp * st, p * sp * ct, dpdk * sp * st,
          0., -p * st, dpdk * ct};
}

/// The invariant mass of \a momenta with \a masses and its sigma from
/// \a covariance (\a stride wide, the momenta's phi, theta, q/p at
/// \a offset + 3 i), and the summed |p| and its sigma.
struct Kinematics {
  double mass = 0.;
  double sMass = 0.;
  double p = 0.;
  double sP = 0.;
};

[[nodiscard]] inline Kinematics kinematicsOf(std::span<const std::array<double, 3>> momenta, std::span<const double> masses, std::span<const double> covariance, int stride, int offset)
{
  const auto count = momenta.size();
  std::array<double, 3> total{};
  double energy = 0.;
  std::vector<std::array<double, 3>> vectors(count);
  std::vector<double> energies(count);
  for (std::size_t i = 0; i < count; ++i) {
    vectors[i] = momentumOf(momenta[i]);
    const double p2 = (vectors[i][0] * vectors[i][0]) + (vectors[i][1] * vectors[i][1]) + (vectors[i][2] * vectors[i][2]);
    energies[i] = std::sqrt(p2 + (masses[i] * masses[i]));
    energy += energies[i];
    for (int k = 0; k < 3; ++k) {
      total[k] += vectors[i][k];
    }
  }
  const double p = std::hypot(total[0], total[1], total[2]);
  const double m2 = (energy * energy) - (p * p);
  Kinematics out{.mass = std::sqrt(std::max(m2, 0.)), .p = p};
  // Gradients by every daughter's phi, theta, q/p.
  std::vector<double> gMass(3 * count);
  std::vector<double> gP(3 * count);
  for (std::size_t i = 0; i < count; ++i) {
    const auto jacobian = momentumJacobian(momenta[i]);
    for (int a = 0; a < 3; ++a) {
      double dm2 = 0.;
      double dp = 0.;
      for (int k = 0; k < 3; ++k) {
        const double dvk = jacobian[k * 3 + a];
        // d(E_i)/dp_ik = p_ik / E_i; dM^2 = 2 E dE - 2 P . dP
        dm2 += 2. * ((energy * vectors[i][k] / energies[i]) - total[k]) * dvk;
        dp += (total[k] / p) * dvk;
      }
      gMass[3 * i + a] = out.mass > 0. ? dm2 / (2. * out.mass) : 0.;
      gP[3 * i + a] = dp;
    }
  }
  double varMass = 0.;
  double varP = 0.;
  for (std::size_t a = 0; a < 3 * count; ++a) {
    for (std::size_t b = 0; b < 3 * count; ++b) {
      const double c = covariance[(offset + a) * stride + (offset + b)];
      varMass += gMass[a] * c * gMass[b];
      varP += gP[a] * c * gP[b];
    }
  }
  out.sMass = std::sqrt(std::max(varMass, 0.));
  out.sP = std::sqrt(std::max(varP, 0.));
  return out;
}

/// The mother's decay length from \a point along its direction, ct with
/// \a mass, and its direction, each with its sigma from \a covariance (the
/// fit's, V and the momenta).
struct MotherQuantities {
  double length = 0.;
  double sLength = 0.;
  double ct = 0.;
  double sCt = 0.;
  double phi = 0.;
  double sPhi = 0.;
  double theta = 0.;
  double sTheta = 0.;
};

[[nodiscard]] inline MotherQuantities motherQuantities(const std::array<double, 3>& vertex, std::span<const std::array<double, 3>> momenta, const std::array<double, 3>& point, double mass, std::span<const double> covariance, int vertices)
{
  using namespace decay_fit_detail;
  using D = Dual<kConstraintInputs>;
  const int count = static_cast<int>(momenta.size());
  const int n = (3 * vertices) + (3 * count);
  std::array<D, 3> v{};
  for (int k = 0; k < 3; ++k) {
    v[k].v = vertex[k];
    v[k].d[k] = 1.;
  }
  std::array<D, 3> total{};
  for (int i = 0; i < count; ++i) {
    std::array<D, 3> q{};
    for (int k = 0; k < 3; ++k) {
      q[k].v = momenta[i][k];
      q[k].d[(3 * vertices) + (3 * i) + k] = 1.;
    }
    const D p = D{q[2].v > 0. ? 1. : -1.} / q[2];
    const D pt = p * sin(q[1]);
    total[0] = total[0] + (pt * cos(q[0]));
    total[1] = total[1] + (pt * sin(q[0]));
    total[2] = total[2] + (p * cos(q[1]));
  }
  const D pt = sqrt((total[0] * total[0]) + (total[1] * total[1]));
  const D p = sqrt((pt * pt) + (total[2] * total[2]));
  const D length = (((v[0] + (-point[0])) * total[0]) + ((v[1] + (-point[1])) * total[1]) + ((v[2] + (-point[2])) * total[2])) / p;
  const D ct = (mass * length) / p;
  const D phi = atan2(total[1], total[0]);
  const D theta = atan2(pt, total[2]);
  const auto sigma = [&](const D& value) {
    double variance = 0.;
    for (int a = 0; a < n; ++a) {
      for (int b = 0; b < n; ++b) {
        variance += value.d[a] * covariance[a * n + b] * value.d[b];
      }
    }
    return std::sqrt(std::max(variance, 0.));
  };
  return {.length = length.v, .sLength = sigma(length), .ct = ct.v, .sCt = sigma(ct), .phi = phi.v, .sPhi = sigma(phi), .theta = theta.v, .sTheta = sigma(theta)};
}

/// The cells of \a lut the line from \a origin along \a direction crosses up
/// to \a end, where a photon may have converted: each cylindrical layer's
/// stretch of it, split where the line crosses the layer's z bins; each with
/// the probability of a pair there, 1 - exp(-7/9 x/X0). Cells without
/// material are left out: no pair is made there.
/// The covariance (loc0, loc0 loc1, loc1) [cm^2] on the plane about the
/// production point \a point of a charged mother of \a momentum [GeV] and \a mass
/// flying straight to \a vertex: a kink theta_k at the fraction t_k of the
/// way, at distance t_k L from the point, moves its crossing by theta_k t_k L
/// across the flight, by that over sin(theta) along the beam (the plane
/// contains the beam's direction); Highland's formula for the segment's whole
/// x/X0, shared in proportion to x/X0 (`beamScattering`'s way):
/// theta0^2 L^2 (int t^2 dx) / X.
[[nodiscard]] inline std::array<double, 3> productionScattering(const o2::base::MatLayerCylSet& lut, const std::array<double, 3>& point, const std::array<double, 3>& vertex, const std::array<double, 3>& momentum, double mass)
{
  const auto moment = lut.getMatMoment(static_cast<float>(point[0]), static_cast<float>(point[1]), static_cast<float>(point[2]), static_cast<float>(vertex[0]), static_cast<float>(vertex[1]), static_cast<float>(vertex[2]));
  const double p = std::hypot(momentum[0], momentum[1], momentum[2]);
  const double sinTheta = std::hypot(momentum[0], momentum[1]) / p;
  if (!(moment.x2X0 > 0.F) || !(sinTheta > 0.)) {
    return {};
  }
  const double theta0 = highlandTheta0(moment.x2X0, moment.x2X0, static_cast<Scalar>(p), static_cast<Scalar>(mass));
  const double length2 = std::pow(vertex[0] - point[0], 2) + std::pow(vertex[1] - point[1], 2) + std::pow(vertex[2] - point[2], 2);
  const double across = theta0 * theta0 * length2 * static_cast<double>(moment.x2X0T2) / static_cast<double>(moment.x2X0);
  return {across, 0., across / (sinTheta * sinTheta)};
}

[[nodiscard]] inline std::vector<LineCell> conversionCells(const o2::base::MatLayerCylSet& lut, const std::array<double, 3>& origin, const std::array<double, 3>& direction, double end)
{
  std::vector<LineCell> cells;
  // r^2(s) = a s^2 + 2 b s + c over the transverse plane.
  const double a = (direction[0] * direction[0]) + (direction[1] * direction[1]);
  const double b = (origin[0] * direction[0]) + (origin[1] * direction[1]);
  const double c = (origin[0] * origin[0]) + (origin[1] * origin[1]);
  // The stretch of s with r(s) below \a r: [low, high], empty if low > high.
  const auto inside = [&](double r) {
    const double discriminant = (b * b) - (a * (c - (r * r)));
    if (!(a > 0.) || discriminant < 0.) {
      return std::array<double, 2>{1., 0.};
    }
    const double root = std::sqrt(discriminant);
    return std::array<double, 2>{(-b - root) / a, (-b + root) / a};
  };
  const auto at = [&](double s) { return std::array<float, 3>{static_cast<float>(origin[0] + (s * direction[0])), static_cast<float>(origin[1] + (s * direction[1])), static_cast<float>(origin[2] + (s * direction[2]))}; };
  for (int i = 0; i < lut.getNLayers(); ++i) {
    const auto& layer = lut.getLayer(i);
    const auto outer = inside(layer.getRMax());
    const auto inner = inside(layer.getRMin());
    // The layer's stretches: inside its outer radius, outside its inner one.
    std::vector<std::array<double, 2>> stretches;
    if (inner[0] > inner[1]) {
      stretches.push_back(outer);
    } else {
      stretches.push_back({outer[0], std::min(outer[1], inner[0])});
      stretches.push_back({std::max(outer[0], inner[1]), outer[1]});
    }
    for (auto stretch : stretches) {
      stretch[0] = std::max(stretch[0], 0.);
      stretch[1] = std::min(stretch[1], end);
      if (!(stretch[1] > stretch[0])) {
        continue;
      }
      // Split at the layer's z bins.
      std::vector<double> edges{stretch[0], stretch[1]};
      if (direction[2] != 0.) {
        for (int k = 0; k <= layer.getNZBins(); ++k) {
          const double s = (static_cast<double>(layer.getZBinMin(k)) - origin[2]) / direction[2];
          if (s > stretch[0] && s < stretch[1]) {
            edges.push_back(s);
          }
        }
      }
      std::ranges::sort(edges);
      for (std::size_t k = 0; k + 1 < edges.size(); ++k) {
        const auto from = at(edges[k]);
        const auto to = at(edges[k + 1]);
        const double x2X0 = lut.getMatBudget(from[0], from[1], from[2], to[0], to[1], to[2]).meanX2X0;
        if (x2X0 > 0.) {
          cells.push_back({.s = 0.5 * (edges[k] + edges[k + 1]), .halfLength = 0.5 * (edges[k + 1] - edges[k]), .probability = 1. - std::exp(-(7. / 9.) * x2X0)});
        }
      }
    }
  }
  return cells;
}

} // namespace decay_fit_study_detail

template <DetectorLayout Det>
class DecayFitStudy
{
 public:
  static constexpr std::size_t kMaxDaughters = 5;

  explicit DecayFitStudy(RecordWriter& writer) : mWriter(&writer) { writer.begin<DecayFitRecord>(); }

  /// \a event's single-vertex decays, from \a chains of \a source.
  void process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source);

  [[nodiscard]] std::vector<std::string> summary() const
  {
    return {std::format("decay fit study: {} decays, {} records into {}; {} with a failed daughter refit; {} left out for a daughter not findable", mDecays, mRecords, mWriter->path(), mRefitFailed, mUnfindable)};
  }

 private:
  RecordWriter* mWriter = nullptr;
  std::size_t mDecays = 0U;
  std::size_t mRecords = 0U;
  std::size_t mRefitFailed = 0U;
  std::size_t mUnfindable = 0U;
};

template <DetectorLayout Det>
void DecayFitStudy<Det>::process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source)
{
  using namespace decay_fit_study_detail;
  if (event.trees == nullptr) {
    throw std::runtime_error("decay fit study: the host has no decay trees (a truth dump of another version?)");
  }
  const auto& trees = *event.trees;
  std::unordered_map<std::uint64_t, const Chain*> chainOf;
  for (const auto& chain : chains) {
    auto& best = chainOf[chain.particle];
    if (best == nullptr || chain.rows.size() > best->rows.size()) {
      best = &chain;
    }
  }

  // A daughter's chain must be findable: its rows on a chain of four linked
  // layers, at any pT.
  const auto rule = comparisonSettings<Det>(4, 0.);
  const auto findable = [&](const TreeParticle& leaf) {
    const auto chain = chainOf.find(leaf.id);
    return chain != chainOf.end() && isFindable(std::numeric_limits<float>::max(), chainLayers(event.rows, chain->second->rows, rule), rule);
  };

  struct Work {
    const TreeParticle* mother = nullptr;
    std::vector<const TreeParticle*> leaves;
    std::vector<const TreeParticle*> composites; ///< a tree's, vertex 1 on
    DecayTopology topology;
    std::array<double, 3> decay{}; ///< the root's true decay point
    std::vector<DecayFitRecord> records;
    bool refitFailed = false;
    std::uint8_t refitPass = 0U;
    std::uint8_t refitStatus = 0U;
  };
  std::vector<Work> work;
  for (const auto& member : trees.particles) {
    if (isFinalCharged(member.pdg)) {
      continue;
    }
    const auto products = decayProducts(trees, member);
    if (products.size() < 2) {
      continue; // no decay vertex: a K0 becoming a K0S, which is a root of its own
    }
    const auto& decay = products.front()->vertex;
    if (decay == member.vertex) {
      continue; // no vertex of its own: a resonance
    }
    const auto shape = decayShape(trees, member, particleCharge);
    if (!shape.closed() || shape.leaves.size() < 2 || shape.leaves.size() > kMaxDaughters) {
      continue;
    }
    // A tree: vertex 0 the root, then the composites depth first. A
    // composite with a single product (a K0 becoming a K0S) has no vertex:
    // its product takes its place. Every other composite has its own, two
    // products or more, is neutral (its link is a line) and flies.
    std::vector<const TreeParticle*> composites;
    std::unordered_map<std::uint64_t, int> vertexOf{{member.id, 0}};
    bool linkable = true;
    for (const auto* composite : shape.composites) {
      const auto own = decayProducts(trees, *composite);
      if (own.size() == 1U) {
        vertexOf[composite->id] = vertexOf.at(composite->mother);
        continue;
      }
      if (particleCharge(composite->pdg) != 0. || own.front()->vertex == composite->vertex) {
        linkable = false;
        break;
      }
      vertexOf[composite->id] = static_cast<int>(composites.size()) + 1;
      composites.push_back(composite);
    }
    if (!linkable || composites.size() + 1U > static_cast<std::size_t>(decay_fit_detail::kMaxConstrainedVertices)) {
      continue;
    }
    DecayTopology topology;
    if (!composites.empty()) {
      topology.parent.push_back(-1);
      for (const auto* composite : composites) {
        topology.parent.push_back(vertexOf.at(composite->mother));
      }
      for (const auto* leaf : shape.leaves) {
        topology.vertexOf.push_back(vertexOf.at(leaf->mother));
      }
    }
    const bool usable = std::ranges::all_of(shape.leaves, [&](const TreeParticle* leaf) {
      return truthSpecies(leaf->pdg).has_value() && chainOf.contains(leaf->id);
    });
    if (!usable) {
      continue;
    }
    if (!std::ranges::all_of(shape.leaves, [&](const TreeParticle* leaf) { return findable(*leaf); })) {
      ++mUnfindable;
      continue;
    }
    work.push_back({.mother = &member, .leaves = shape.leaves, .composites = std::move(composites), .topology = std::move(topology), .decay = decay});
  }

  const auto& frame = *event.frame;
  const auto& geometry = *event.geometry;
  const auto& material = *event.material;
  o2::gtrk::detail::parallelFor(*event.arena, 0U, work.size(), 4, [&](std::size_t begin, std::size_t end) {
    for (std::size_t index = begin; index < end; ++index) {
      auto& one = work[index];
      const auto count = one.leaves.size();
      std::vector<DecayDaughter> daughters(count);
      std::vector<RowsDaughter> rowsDaughters(count);
      std::vector<double> masses(count);
      std::uint8_t minRows = std::numeric_limits<std::uint8_t>::max();
      for (const auto* leaf : one.leaves) {
        minRows = std::min<std::uint8_t>(minRows, static_cast<std::uint8_t>(std::min<std::size_t>(chainOf.at(leaf->id)->rows.size(), 255U)));
      }
      for (std::size_t i = 0; i < count; ++i) {
        const auto* leaf = one.leaves[i];
        const auto& chain = *chainOf.at(leaf->id);
        const auto species = *truthSpecies(leaf->pdg);
        const RowChain<Det> rows{.frame = &frame, .rows = chain.rows.data()};
        const auto refit = refitWithHypothesis(event.propagator, geometry, material, rows, static_cast<int>(chain.rows.size()), species.pid, 0.F, 0.F);
        if (refit.failedPass != FitPass::kNone || (refit.status != PropagationStatus::kOk && refit.innermost.surface == kInvalidSensorID)) {
          one.refitFailed = true;
          one.refitPass = static_cast<std::uint8_t>(refit.failedPass);
          one.refitStatus = static_cast<std::uint8_t>(refit.status);
          break;
        }
        auto propagator = event.propagator;
        propagator.massGeV = species.mass;
        daughters[i] = {.innermost = refit.innermost, .surface = geometry.planeOf(rows[0]), .propagator = propagator, .mass = species.mass};
        if (species.pid == o2::track::PID::Electron) {
          daughters[i].mixture = refit.innermostMixture;
        }
        auto& own = rowsDaughters[i];
        own.propagator = propagator;
        own.electron = species.pid == o2::track::PID::Electron;
        for (const auto at : chain.rows) {
          const auto m = frame.measurement(at);
          own.planes.push_back(geometry.planeOf(m));
          own.positions.push_back({m.loc0, m.loc1});
          own.covariance.push_back({m.varLoc0, m.covLoc01, m.varLoc1});
        }
        masses[i] = species.mass;
      }

      // The truth.
      const auto& mother = *one.mother;
      const auto& decay = one.decay;
      std::array<double, 3> pTrue{};
      double eTrue = 0.;
      double minPt = std::numeric_limits<double>::infinity();
      std::vector<std::array<double, 3>> qTrue(count);
      for (std::size_t i = 0; i < count; ++i) {
        const auto& m = one.leaves[i]->momentum;
        const double p = std::hypot(m[0], m[1], m[2]);
        eTrue += std::sqrt((p * p) + (masses[i] * masses[i]));
        for (int k = 0; k < 3; ++k) {
          pTrue[k] += m[k];
        }
        minPt = std::min(minPt, std::hypot(m[0], m[1]));
        qTrue[i] = {std::atan2(m[1], m[0]), std::atan2(std::hypot(m[0], m[1]), m[2]), particleCharge(one.leaves[i]->pdg) / p};
      }
      const double pTotal = std::hypot(pTrue[0], pTrue[1], pTrue[2]);
      // The flight's frame: along it, across it in the transverse plane, and
      // the third axis.
      std::array<double, 3> along{decay[0] - mother.vertex[0], decay[1] - mother.vertex[1], decay[2] - mother.vertex[2]};
      const double flight = std::hypot(along[0], along[1], along[2]);
      for (auto& component : along) {
        component /= flight;
      }
      std::array<double, 3> across{-along[1], along[0], 0.};
      const double acrossNorm = std::hypot(across[0], across[1]);
      if (acrossNorm > 0.) {
        across = {across[0] / acrossNorm, across[1] / acrossNorm, 0.};
      } else {
        across = {1., 0., 0.};
      }
      const std::array<double, 3> third{(along[1] * across[2]) - (along[2] * across[1]), (along[2] * across[0]) - (along[0] * across[2]), (along[0] * across[1]) - (along[1] * across[0])};

      DecayFitRecord base{
        .timeframe = event.timeframe,
        .source = static_cast<std::uint8_t>(0),
        .particle = mother.id,
        .pdg = mother.pdg,
        .daughters = static_cast<std::uint8_t>(count),
        .flight = static_cast<float>(flight),
        .decayR = static_cast<float>(std::hypot(decay[0], decay[1])),
        .pTrue = static_cast<float>(pTotal),
        .minDaughterPt = static_cast<float>(minPt),
        .minDaughterRows = minRows,
        .refitPass = one.refitPass,
        .refitStatus = one.refitStatus,
        .massTrue = static_cast<float>(std::sqrt(std::max((eTrue * eTrue) - (pTotal * pTotal), 0.))),
      };
      // The daughters' energy loss before their first row, true and modelled,
      // and their Kalman fits' momenta there.
      for (std::size_t i = 0; i < std::min(count, DecayFitRecord::kDaughterColumns); ++i) {
        base.lossTrue[i] = base.lossModel[i] = base.dPRelInner[i] = std::numeric_limits<float>::quiet_NaN();
        if (one.refitFailed || event.hits == nullptr) {
          continue;
        }
        const auto* leaf = one.leaves[i];
        const auto first = frame.measurement(chainOf.at(leaf->id)->rows.front());
        const auto* segment = event.hits->crossing(leaf->id, first.surface, {static_cast<double>(first.x), static_cast<double>(first.y), static_cast<double>(first.z)});
        if (segment == nullptr) {
          continue;
        }
        const double pProduction = std::hypot(leaf->momentum[0], leaf->momentum[1], leaf->momentum[2]);
        const double pFirst = std::hypot(segment->momentum[0], segment->momentum[1], segment->momentum[2]);
        base.lossTrue[i] = static_cast<float>(1. - (pFirst / pProduction));
        const auto budget = material.getMatBudget(static_cast<float>(leaf->vertex[0]), static_cast<float>(leaf->vertex[1]), static_cast<float>(leaf->vertex[2]), static_cast<float>(segment->start[0]), static_cast<float>(segment->start[1]), static_cast<float>(segment->start[2]));
        const auto loss = energyLoss(static_cast<Scalar>(1. / pProduction), static_cast<Scalar>(masses[i]), -budget.meanRho * budget.length);
        if (loss.ok) {
          base.lossModel[i] = static_cast<float>(1. - (1. / static_cast<double>(loss.momentumScale)));
        }
        base.dPRelInner[i] = static_cast<float>((1. / std::abs(static_cast<double>(daughters[i].innermost.qOverP())) / pFirst) - 1.);
        base.hitsBefore[i] = static_cast<std::int16_t>(std::ranges::count_if(event.hits->segments(leaf->id), [&](const HitSegment* other) { return other->time < segment->time; }));
      }

      // The mother: its PDG mass, its collision's vertex, the true decay length.
      const auto* species = o2::O2DatabasePDG::Instance()->GetParticle(mother.pdg);
      const double motherMass = species == nullptr ? 0. : species->Mass();
      const auto collision = event.collisions->find(collisionOf(mother.id));
      const bool hasCollision = collision != event.collisions->end();
      const std::array<double, 3> point = hasCollision ? collision->second : std::array<double, 3>{};
      // Made at the collision: its production point is the collision's to
      // the precision the kinematics store it (float); otherwise (a Lambda
      // from a Xi) a production point at the collision is false.
      const bool madeAtCollision = hasCollision && std::ranges::equal(mother.vertex, point, [](double a, double b) { return static_cast<float>(a) == static_cast<float>(b); });
      if (hasCollision) {
        base.productionOffset = static_cast<float>(std::hypot(mother.vertex[0] - point[0], mother.vertex[1] - point[1], mother.vertex[2] - point[2]));
        base.decayLengthTrue = static_cast<float>((((decay[0] - point[0]) * pTrue[0]) + ((decay[1] - point[1]) * pTrue[1]) + ((decay[2] - point[2]) * pTrue[2])) / pTotal);
        base.ctTrue = static_cast<float>(base.decayLengthTrue * motherMass / pTotal);
      }
      // A conversion: its line from the collision towards its daughters'
      // first rows, and the material table's cells on it up to the nearer one.
      struct ConversionLine {
        std::array<double, 3> direction{};
        std::vector<LineCell> cells;
        double along = 0.; ///< the true point's distance along it
      };
      std::optional<ConversionLine> conversion;
      if (mother.pdg == 22 && hasCollision && !one.refitFailed) {
        std::array<double, 3> firstRows{};
        std::vector<std::array<double, 3>> firsts;
        for (std::size_t i = 0; i < count; ++i) {
          const auto m = frame.measurement(chainOf.at(one.leaves[i]->id)->rows.front());
          firsts.push_back({m.x, m.y, m.z});
          for (int k = 0; k < 3; ++k) {
            firstRows[k] += firsts.back()[k] / static_cast<double>(count);
          }
        }
        ConversionLine line;
        const double length = std::hypot(firstRows[0] - point[0], firstRows[1] - point[1], firstRows[2] - point[2]);
        for (int k = 0; k < 3; ++k) {
          line.direction[k] = (firstRows[k] - point[k]) / length;
        }
        double nearest = std::numeric_limits<double>::infinity();
        for (const auto& first : firsts) {
          nearest = std::min(nearest, ((first[0] - point[0]) * line.direction[0]) + ((first[1] - point[1]) * line.direction[1]) + ((first[2] - point[2]) * line.direction[2]));
        }
        line.cells = conversionCells(material, point, line.direction, nearest);
        for (int k = 0; k < 3; ++k) {
          line.along += (decay[k] - point[k]) * line.direction[k];
        }
        conversion = std::move(line);
      }
      // The fits: every variant alone, then kFull and kRows under the mass,
      // the production point, and both.
      struct Job {
        DecayFitVariant variant;
        std::uint8_t constraints;
      };
      const bool tree = !one.composites.empty();
      // A charged root's own rows (a tracked Xi, Omega): its Kalman fit with its
      // mass, the outward pass's end on its last row.
      // For the rows' fit its rows directly, any number, its last row first.
      std::optional<DecayMother> trackedMother;
      std::optional<RowsMother> trackedRows;
      if (tree && !one.refitFailed && particleCharge(mother.pdg) != 0. && motherMass > 0.) {
        if (const auto found = chainOf.find(mother.id); found != chainOf.end()) {
          const auto& own = *found->second;
          base.motherRows = static_cast<std::uint8_t>(std::min<std::size_t>(own.rows.size(), 255U));
          auto propagator = event.propagator;
          propagator.massGeV = static_cast<Scalar>(motherMass);
          if (!own.rows.empty()) {
            RowsMother byRows{.rows = {.propagator = propagator}, .vertex = 0};
            for (auto at = own.rows.rbegin(); at != own.rows.rend(); ++at) {
              const auto m = frame.measurement(*at);
              byRows.rows.planes.push_back(geometry.planeOf(m));
              byRows.rows.positions.push_back({m.loc0, m.loc1});
              byRows.rows.covariance.push_back({m.varLoc0, m.covLoc01, m.varLoc1});
            }
            trackedRows = std::move(byRows);
          }
          const RowChain<Det> rows{.frame = &frame, .rows = own.rows.data()};
          const auto ownFit = fitMeasurements(propagator, geometry, material, rows, static_cast<int>(own.rows.size()));
          if (ownFit.status == PropagationStatus::kOk) {
            trackedMother = DecayMother{.outermost = ownFit.outwardState, .surface = geometry.planeOf(ownFit.outwardState), .propagator = propagator, .vertex = 0};
          }
        }
      }
      std::vector<Job> jobs;
      for (const auto variant : {DecayFitVariant::kFull, DecayFitVariant::kCarryOnce, DecayFitVariant::kNoMaterial, DecayFitVariant::kRows}) {
        jobs.push_back({variant, 0U});
      }
      if (tree) {
        // A tree: its composites' masses, with the root's mass, its
        // production point, both; the sequential fit (a composite first) for
        // a tree of two vertices.
        for (const auto variant : {DecayFitVariant::kFull, DecayFitVariant::kRows}) {
          for (const std::uint8_t set : {std::uint8_t{16}, std::uint8_t{17}, std::uint8_t{18}, std::uint8_t{19}}) {
            jobs.push_back({variant, set});
          }
        }
        if (one.topology.vertices() == 2) {
          for (const std::uint8_t set : {std::uint8_t{0}, std::uint8_t{16}, std::uint8_t{18}, std::uint8_t{19}}) {
            jobs.push_back({DecayFitVariant::kSequential, set});
          }
        }
        // A tracked mother: its own rows in the tree, and in the sequential
        // fit's root (strangeness tracking as a fit).
        for (const std::uint8_t set : {std::uint8_t{48}, std::uint8_t{50}, std::uint8_t{51}}) {
          if (trackedMother) {
            jobs.push_back({DecayFitVariant::kFull, set});
            if (one.topology.vertices() == 2) {
              jobs.push_back({DecayFitVariant::kSequential, set});
            }
          }
          if (trackedRows) {
            jobs.push_back({DecayFitVariant::kRows, set});
          }
        }
      } else {
        // A photon: collinear daughters in place of the mass.
        const std::uint8_t first = mother.pdg == 22 ? std::uint8_t{4} : std::uint8_t{1};
        for (const auto variant : {DecayFitVariant::kFull, DecayFitVariant::kRows}) {
          for (const std::uint8_t set : {first, std::uint8_t{2}, static_cast<std::uint8_t>(first | 2U)}) {
            jobs.push_back({variant, set});
          }
        }
      }
      // A tree's composites: their PDG masses, the first one's truth.
      std::vector<double> compositeMasses;
      for (const auto* composite : one.composites) {
        const auto* entry = o2::O2DatabasePDG::Instance()->GetParticle(composite->pdg);
        compositeMasses.push_back(entry == nullptr ? 0. : entry->Mass());
      }
      std::vector<std::size_t> childLeaves; // the first composite's leaves
      std::array<double, 3> childDecay{};
      std::array<double, 3> childAxis{}; // its true flight direction
      std::array<double, 3> childAcross{};
      if (tree) {
        base.childPdg = one.composites.front()->pdg;
        childDecay = decayProducts(trees, *one.composites.front()).front()->vertex;
        std::array<double, 3> childP{};
        double childE = 0.;
        for (std::size_t i = 0; i < count; ++i) {
          if (one.topology.under(static_cast<int>(i), 1)) {
            childLeaves.push_back(i);
            const auto& m = one.leaves[i]->momentum;
            const double p = std::hypot(m[0], m[1], m[2]);
            childE += std::sqrt((p * p) + (masses[i] * masses[i]));
            for (int c = 0; c < 3; ++c) {
              childP[c] += m[c];
            }
          }
        }
        const double childPTotal = std::hypot(childP[0], childP[1], childP[2]);
        for (int c = 0; c < 3; ++c) {
          childAxis[c] = childP[c] / childPTotal;
        }
        const double t = std::hypot(childAxis[0], childAxis[1]);
        childAcross = t > 0. ? std::array<double, 3>{-childAxis[1] / t, childAxis[0] / t, 0.} : std::array<double, 3>{1., 0., 0.};
        base.childLengthTrue = static_cast<float>(((childDecay[0] - decay[0]) * childAxis[0]) + ((childDecay[1] - decay[1]) * childAxis[1]) + ((childDecay[2] - decay[2]) * childAxis[2]));
        base.childMassTrue = static_cast<float>(std::sqrt(std::max((childE * childE) - (childPTotal * childPTotal), 0.)));
      }
      std::array<DecayVertexFitResult, kDecayFitVariants> fits{};
      for (const auto& job : jobs) {
        DecayFitRecord record = base;
        record.variant = static_cast<std::uint8_t>(job.variant);
        record.constraints = job.constraints;
        if (one.refitFailed) {
          record.status = DecayFitRecord::kRefitFailed;
          one.records.push_back(record);
          continue;
        }
        const bool wantMass = (job.constraints & 1U) != 0U;
        const bool wantProduction = (job.constraints & 2U) != 0U;
        const bool wantCollinear = (job.constraints & 4U) != 0U;
        const bool wantComposites = (job.constraints & 16U) != 0U;
        if ((wantMass && !(motherMass > 0.)) || (wantProduction && !madeAtCollision) || (wantCollinear && count != 2)) {
          record.status = static_cast<std::uint8_t>(DecayFitStatus::kNoStart); // nothing to constrain to
          one.records.push_back(record);
          continue;
        }
        DecayConstraints constraints{.mass = wantMass, .massGeV = motherMass, .production = wantProduction, .point = point, .collinear = wantCollinear};
        std::optional<RowsMother> jobMother = (job.constraints & 32U) != 0U ? trackedRows : std::nullopt;
        // A charged mother's production point measured, with its scattering
        // on the way from it to the unconstrained fit's vertex at its
        // momentum; a neutral one does not scatter: held exactly.
        if (const auto& free = fits[static_cast<std::size_t>(DecayFitVariant::kFull)]; wantProduction && particleCharge(mother.pdg) != 0. && free.status == DecayFitStatus::kOk && !free.vertices.empty()) {
          std::array<double, 3> total{};
          for (const auto& m : free.momenta) {
            const double pm = 1. / std::abs(m[2]);
            total[0] += pm * std::sin(m[1]) * std::cos(m[0]);
            total[1] += pm * std::sin(m[1]) * std::sin(m[0]);
            total[2] += pm * std::cos(m[1]);
          }
          constraints.pointCovariance = productionScattering(material, point, free.vertices[0], total, motherMass > 0. ? motherMass : static_cast<double>(base.massTrue));
          record.pointSigma = static_cast<float>(std::sqrt(constraints.pointCovariance[0]));
          // A tracked mother's rows' walk ends at the point: the scattering
          // between its innermost row and the point.
          if (jobMother) {
            Scalar innermost[3]{};
            surfaceToGlobal(jobMother->rows.planes.back(), static_cast<Scalar>(jobMother->rows.positions.back()[0]), static_cast<Scalar>(jobMother->rows.positions.back()[1]), innermost);
            jobMother->pointCovariance = productionScattering(material, point, {innermost[0], innermost[1], innermost[2]}, total, motherMass);
            if (job.variant == DecayFitVariant::kRows) {
              record.pointSigma = static_cast<float>(std::sqrt(jobMother->pointCovariance[0])); // the walk's end's, which this fit takes
            }
          }
        }
        if (wantComposites && job.variant != DecayFitVariant::kSequential) {
          constraints.compositeMasses = compositeMasses;
        }
        const auto& full = fits[static_cast<std::size_t>(DecayFitVariant::kFull)];
        const bool wantMother = (job.constraints & 32U) != 0U;
        const auto mothers = wantMother && trackedMother ? std::span<const DecayMother>{&*trackedMother, 1} : std::span<const DecayMother>{};
        const auto rowsMothers = jobMother ? std::span<const RowsMother>{&*jobMother, 1} : std::span<const RowsMother>{};
        const auto fitUnder = [&](const DecayConstraints& set) {
          return job.variant == DecayFitVariant::kRows ? fitDecayRows(std::span<const RowsDaughter>{rowsDaughters}, material, full, set, one.topology, rowsMothers) : fitDecayVertex(std::span<const DecayDaughter>{daughters}, material, job.variant, set, one.topology, mothers);
        };
        DecayVertexFitResult outcome;
        const auto started = std::chrono::steady_clock::now();
        // The sequential fit: the composite (vertex 1) from its daughters,
        // under its mass if asked; then the root from its own daughters and
        // the composite as a pseudo-track.
        DecayVertexFitResult childFit;
        std::vector<double> pseudoJacobian; // the pseudo-track by the child fit's parameters
        std::vector<double> rootMasses;
        const bool sequential = job.variant == DecayFitVariant::kSequential;
        if (sequential) {
          std::vector<DecayDaughter> childDaughters;
          std::vector<DecayDaughter> rootDaughters;
          for (std::size_t i = 0; i < count; ++i) {
            if (one.topology.vertexOfDaughter(static_cast<int>(i)) == 1) {
              childDaughters.push_back(daughters[i]);
            } else {
              rootDaughters.push_back(daughters[i]);
              rootMasses.push_back(masses[i]);
            }
          }
          const DecayConstraints childSet{.mass = wantComposites, .massGeV = compositeMasses.front()};
          childFit = fitDecayVertex(std::span<const DecayDaughter>{childDaughters}, material, DecayFitVariant::kFull, childSet);
          if (childFit.status != DecayFitStatus::kOk) {
            outcome.status = childFit.status;
          } else {
            rootDaughters.push_back(pseudoTrackOf(childFit, compositeMasses.front(), event.propagator, &pseudoJacobian));
            rootMasses.push_back(compositeMasses.front());
            outcome = fitDecayVertex(std::span<const DecayDaughter>{rootDaughters}, material, DecayFitVariant::kFull, constraints, {}, mothers);
          }
        } else if (conversion) {
          record.constraints |= 8U;
          const auto profile = profileAlongLine(fitUnder, constraints, point, conversion->direction, conversion->cells);
          outcome = profile.fit;
          if (outcome.status != DecayFitStatus::kNoProfile) {
            record.profilePull = static_cast<float>((conversion->along - profile.mean) / profile.sigma);
            record.profileSigma = static_cast<float>(profile.sigma);
          }
        } else {
          outcome = fitUnder(constraints);
        }
        record.microseconds = static_cast<float>(std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - started).count());
        if (job.constraints == 0U) {
          fits[static_cast<std::size_t>(job.variant)] = outcome;
        }
        const auto& fit = outcome;
        record.status = static_cast<std::uint8_t>(fit.status);
        record.iterations = fit.iterations;
        record.stopReason = fit.stopReason;
        record.outer = fit.outer;
        record.modelShift = fit.modelShift;
        record.stopMetric = fit.stopMetric;
        if (fit.status != DecayFitStatus::kOk) {
          one.records.push_back(record);
          continue;
        }
        const int vertexCount = static_cast<int>(fit.vertices.size());
        const std::size_t fitCount = fit.momenta.size();
        const int n = (3 * vertexCount) + (3 * static_cast<int>(fitCount));
        const auto& c = fit.covariance;
        const auto& fitMasses = sequential ? rootMasses : masses;
        record.chi2 = static_cast<float>(fit.chi2);
        record.ndf = static_cast<std::int16_t>(fit.ndf);
        record.maxCarryX2X0 = static_cast<float>(*std::ranges::max_element(fit.carryX2X0));
        const std::array<double, 3> d{fit.vertex[0] - decay[0], fit.vertex[1] - decay[1], fit.vertex[2] - decay[2]};
        record.dX = static_cast<float>(d[0]);
        record.dY = static_cast<float>(d[1]);
        record.dZ = static_cast<float>(d[2]);
        record.sX = static_cast<float>(std::sqrt(c[0]));
        record.sY = static_cast<float>(std::sqrt(c[n + 1]));
        record.sZ = static_cast<float>(std::sqrt(c[2 * n + 2]));
        const auto project = [&](const std::array<double, 3>& axis, float& residual, float& sigma) {
          residual = static_cast<float>((d[0] * axis[0]) + (d[1] * axis[1]) + (d[2] * axis[2]));
          double variance = 0.;
          for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
              variance += axis[a] * c[a * n + b] * axis[b];
            }
          }
          sigma = static_cast<float>(std::sqrt(std::max(variance, 0.)));
        };
        project(along, record.dL, record.sL);
        project(across, record.dT, record.sT);
        project(third, record.dN, record.sN);
        const auto fitted = kinematicsOf(fit.momenta, fitMasses, c, n, 3 * vertexCount);
        record.mass = static_cast<float>(fitted.mass);
        record.sMass = static_cast<float>(fitted.sMass);
        record.p = static_cast<float>(fitted.p);
        record.sP = static_cast<float>(fitted.sP);
        if (hasCollision) {
          const auto motherFit = motherQuantities(fit.vertex, fit.momenta, point, motherMass, c, vertexCount);
          record.decayLength = static_cast<float>(motherFit.length);
          record.sDecayLength = static_cast<float>(motherFit.sLength);
          record.ct = static_cast<float>(motherFit.ct);
          record.sCt = static_cast<float>(motherFit.sCt);
          record.dPhiMother = static_cast<float>(std::remainder(motherFit.phi - std::atan2(pTrue[1], pTrue[0]), 2. * std::numbers::pi));
          record.dThetaMother = static_cast<float>(motherFit.theta - std::atan2(std::hypot(pTrue[0], pTrue[1]), pTrue[2]));
          record.sPhiMother = static_cast<float>(motherFit.sPhi);
          record.sThetaMother = static_cast<float>(motherFit.sTheta);
        }
        // The child: its vertex, decay length and mass.
        if (tree) {
          const auto& child = sequential ? childFit.vertex : fit.vertices[1];
          const int childN = sequential ? 3 + (3 * static_cast<int>(childFit.momenta.size())) : n;
          const int childAt = sequential ? 0 : 3;
          const auto& childC = sequential ? childFit.covariance : c;
          const std::array<double, 3> e{child[0] - childDecay[0], child[1] - childDecay[1], child[2] - childDecay[2]};
          const auto childVariance = [&](const std::array<double, 3>& axis) {
            double variance = 0.;
            for (int a = 0; a < 3; ++a) {
              for (int b = 0; b < 3; ++b) {
                variance += axis[a] * childC[((childAt + a) * childN) + childAt + b] * axis[b];
              }
            }
            return variance;
          };
          record.dChildL = static_cast<float>((e[0] * childAxis[0]) + (e[1] * childAxis[1]) + (e[2] * childAxis[2]));
          record.sChildL = static_cast<float>(std::sqrt(std::max(childVariance(childAxis), 0.)));
          record.dChildT = static_cast<float>((e[0] * childAcross[0]) + (e[1] * childAcross[1]) + (e[2] * childAcross[2]));
          record.sChildT = static_cast<float>(std::sqrt(std::max(childVariance(childAcross), 0.)));
          // The decay length from the parent's vertex along the true flight;
          // its variance the two vertices' with their correlation (one fit:
          // its covariance; sequential: the parent's vertex depends on the
          // child's through the pseudo-track, the parent fit's last daughter:
          // cov(V_p, V_c) = G J C_c, G the parent's gain on it, J the
          // pseudo-track by the child's parameters).
          record.childLength = static_cast<float>(((child[0] - fit.vertex[0]) * childAxis[0]) + ((child[1] - fit.vertex[1]) * childAxis[1]) + ((child[2] - fit.vertex[2]) * childAxis[2]));
          double lengthVariance = childVariance(childAxis);
          std::array<double, 9> crossCovariance{};
          if (sequential) {
            const auto& gain = fit.vertexGain.back();
            for (int a = 0; a < 3; ++a) {
              for (int b = 0; b < 3; ++b) {
                double value = 0.;
                for (int r = 0; r < 5; ++r) {
                  for (int q = 0; q < childN; ++q) {
                    value += gain[(a * 5) + r] * pseudoJacobian[(r * childN) + q] * childC[(q * childN) + b];
                  }
                }
                crossCovariance[(a * 3) + b] = value;
              }
            }
          } else {
            for (int a = 0; a < 3; ++a) {
              for (int b = 0; b < 3; ++b) {
                crossCovariance[(a * 3) + b] = c[(a * n) + 3 + b];
              }
            }
          }
          for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
              lengthVariance += childAxis[a] * c[(a * n) + b] * childAxis[b];
              lengthVariance -= 2. * childAxis[a] * crossCovariance[(a * 3) + b] * childAxis[b];
            }
          }
          record.sChildLength = static_cast<float>(std::sqrt(std::max(lengthVariance, 0.)));
          if (sequential) {
            std::vector<double> childMasses;
            for (const auto i : childLeaves) {
              childMasses.push_back(masses[i]);
            }
            const auto childKinematics = kinematicsOf(childFit.momenta, childMasses, childFit.covariance, childN, 3);
            record.childMass = static_cast<float>(childKinematics.mass);
            record.sChildMass = static_cast<float>(childKinematics.sMass);
          } else {
            // the child's leaves' momenta and their covariance block
            std::vector<std::array<double, 3>> childMomenta;
            std::vector<double> childMasses;
            for (const auto i : childLeaves) {
              childMomenta.push_back(fit.momenta[i]);
              childMasses.push_back(masses[i]);
            }
            const int m = 3 * static_cast<int>(childLeaves.size());
            std::vector<double> block(static_cast<std::size_t>(m * m));
            for (int a = 0; a < m; ++a) {
              for (int b = 0; b < m; ++b) {
                const int ia = (3 * vertexCount) + (3 * static_cast<int>(childLeaves[a / 3])) + (a % 3);
                const int ib = (3 * vertexCount) + (3 * static_cast<int>(childLeaves[b / 3])) + (b % 3);
                block[(a * m) + b] = c[(ia * n) + ib];
              }
            }
            const auto childKinematics = kinematicsOf(childMomenta, childMasses, block, m, 0);
            record.childMass = static_cast<float>(childKinematics.mass);
            record.sChildMass = static_cast<float>(childKinematics.sMass);
          }
        }
        if (sequential) {
          one.records.push_back(record); // no per-daughter columns: the composite's daughters are not in the root's fit
          continue;
        }
        // The carried momenta: independent of one another, each its own 3 x 3.
        std::vector<double> carried(static_cast<std::size_t>(9 * count * count), 0.);
        const int stride = 3 * static_cast<int>(count);
        for (std::size_t i = 0; i < count; ++i) {
          for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
              carried[(3 * i + a) * stride + (3 * i + b)] = fit.carriedCovariance[i][a * 3 + b];
            }
          }
        }
        const auto measured = kinematicsOf(fit.carried, masses, carried, stride, 0);
        record.massCarried = static_cast<float>(measured.mass);
        record.sMassCarried = static_cast<float>(measured.sMass);
        for (std::size_t i = 0; i < std::min(count, DecayFitRecord::kDaughterColumns); ++i) {
          const int q = (3 * vertexCount) + (3 * static_cast<int>(i));
          record.dPhi[i] = static_cast<float>(std::remainder(fit.momenta[i][0] - qTrue[i][0], 2. * std::numbers::pi));
          record.dTheta[i] = static_cast<float>(fit.momenta[i][1] - qTrue[i][1]);
          record.dQOverP[i] = static_cast<float>(fit.momenta[i][2] - qTrue[i][2]);
          record.sPhi[i] = static_cast<float>(std::sqrt(c[q * n + q]));
          record.sTheta[i] = static_cast<float>(std::sqrt(c[(q + 1) * n + (q + 1)]));
          record.sQOverP[i] = static_cast<float>(std::sqrt(c[(q + 2) * n + (q + 2)]));
          record.dPRel[i] = static_cast<float>((std::abs(qTrue[i][2]) / std::abs(fit.momenta[i][2])) - 1.);
        }
        one.records.push_back(record);
      }
    }
  });
  for (auto& one : work) {
    for (auto& record : one.records) {
      record.source = static_cast<std::uint8_t>(source);
      mWriter->write(record);
    }
    mRecords += one.records.size();
    mRefitFailed += one.refitFailed ? 1U : 0U;
  }
  mDecays += work.size();
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYFITSTUDY_H_
