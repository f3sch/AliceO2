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
/// `kMaxDaughters` of them (a photon's conversion is parked), every
/// one with a chain from the source that is findable: its rows -- the chain
/// the fit takes, a truth chain's first outward leg -- on a chain of at least
/// four linked layers (`chainLayers`, `isFindable`, no pT floor): a track of
/// three is barely defined. An electron daughter
/// is refitted with the Gaussian-sum filter and carried as its mixture
/// (DecayVertexFit.h). Each daughter's rows are refitted with its true
/// species (`refitWithHypothesis`); the vertex is fitted from the daughters'
/// states alone, under the mother's mass, its production at the MC collision
/// vertex, and both (DecayConstraints.h), a record each.
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
  std::uint8_t constraints = 0U; ///< 1 the mass, 2 the production point, 16 a tree's composites' masses, 32 the tracked mother's own rows, 64 an unseen neutral's momentum fitted
  std::uint8_t daughters = 0U;
  std::uint8_t status = 0U; ///< DecayFitStatus, or kRefitFailed
  std::uint8_t iterations = 0U;
  std::uint8_t outer = 0U;
  std::uint8_t stopReason = 0U;
  float stopMetric = 0.F;
  float modelShift = 0.F;
  float microseconds = 0.F;        ///< the fit's time
  std::uint16_t combinations = 0U; ///< the Gaussian sum's component combinations
  std::int8_t branch = -1;         ///< an unseen neutral's start (0, 1), -1 without one
  std::uint8_t resonance = 0U;     ///< 1 the resonances' masses taken (Breit-Wigner)
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
  float pointSigma = 0.F;        ///< the production point's equations' sigma across the flight [cm] (0: held exactly)
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
    visit("resonance", resonance);
    visit("branch", branch);
    visit("combinations", combinations);
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

/// q/p's numerator for a daughter of \a charge: |q|, 1 if neutral (its
/// q/p is 1/p).
[[nodiscard]] inline double numeratorOf(double charge) { return charge != 0. ? std::abs(charge) : 1.; }

/// A daughter's momentum vector from phi, theta, q/p and its \a charge.
[[nodiscard]] inline std::array<double, 3> momentumOf(const std::array<double, 3>& q, double charge)
{
  const double p = numeratorOf(charge) / std::abs(q[2]);
  return {p * std::cos(q[0]) * std::sin(q[1]), p * std::sin(q[0]) * std::sin(q[1]), p * std::cos(q[1])};
}

/// d(px, py, pz) / d(phi, theta, q/p), row-major 3 x 3.
[[nodiscard]] inline std::array<double, 9> momentumJacobian(const std::array<double, 3>& q, double charge)
{
  const double p = numeratorOf(charge) / std::abs(q[2]);
  const double dpdk = -p / q[2]; // d|p| / d(q/p)
  const double cp = std::cos(q[0]);
  const double sp = std::sin(q[0]);
  const double ct = std::cos(q[1]);
  const double st = std::sin(q[1]);
  return {-p * sp * st, p * cp * ct, dpdk * cp * st,
          p * cp * st, p * sp * ct, dpdk * sp * st,
          0., -p * st, dpdk * ct};
}

/// The invariant mass of \a momenta with \a masses (and \a charges) and its sigma from
/// \a covariance (\a stride wide, the momenta's phi, theta, q/p at
/// \a offset + 3 i), and the summed |p| and its sigma.
struct Kinematics {
  double mass = 0.;
  double sMass = 0.;
  double p = 0.;
  double sP = 0.;
};

[[nodiscard]] inline Kinematics kinematicsOf(std::span<const std::array<double, 3>> momenta, std::span<const double> masses, std::span<const double> charges, std::span<const double> covariance, int stride, int offset)
{
  const auto count = momenta.size();
  std::array<double, 3> total{};
  double energy = 0.;
  std::vector<std::array<double, 3>> vectors(count);
  std::vector<double> energies(count);
  for (std::size_t i = 0; i < count; ++i) {
    vectors[i] = momentumOf(momenta[i], charges[i]);
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
    const auto jacobian = momentumJacobian(momenta[i], charges[i]);
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

[[nodiscard]] inline MotherQuantities motherQuantities(const std::array<double, 3>& vertex, std::span<const std::array<double, 3>> momenta, std::span<const double> charges, const std::array<double, 3>& point, double mass, std::span<const double> covariance, int vertices)
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
    const D p = D{q[2].v > 0. ? numeratorOf(charges[i]) : -numeratorOf(charges[i])} / q[2];
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
    std::array<double, 3> decay{};               ///< the root's true decay point
    const TreeParticle* unseen = nullptr;        ///< the root's one unseen neutral product, if any
    std::vector<const TreeParticle*> resonances; ///< composites that decay where they are made: transparent
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
    if (member.pdg == 22) {
      continue; // a conversion: parked (its point along the flight needs a profile)
    }
    const auto& decay = products.front()->vertex;
    if (decay == member.vertex) {
      continue; // no vertex of its own: a resonance
    }
    const auto shape = decayShape(trees, member, particleCharge);
    // Closed, or with one unseen neutral straight from the root (a neutrino,
    // a photon, a pi0 as a whole): its momentum the fit's, from the mother's
    // mass and production point.
    const bool missing = shape.chargedUnseen == 0U && shape.unseen.size() == 1U && shape.unseen.front()->mother == member.id;
    if (!(shape.closed() || missing) || shape.leaves.size() < 2 || shape.leaves.size() + (missing ? 1U : 0U) > kMaxDaughters) {
      continue;
    }
    // A tree: vertex 0 the root, then the composites depth first. A
    // composite with a single product (a K0 becoming a K0S) has no vertex:
    // its product takes its place. Every other composite has its own, two
    // products or more, and flies: a neutral one's link a line, a charged
    // one's its helix (the fit's chains).
    std::vector<const TreeParticle*> composites;
    std::unordered_map<std::uint64_t, int> vertexOf{{member.id, 0}};
    bool linkable = true;
    std::vector<const TreeParticle*> resonances;
    for (const auto* composite : shape.composites) {
      const auto own = decayProducts(trees, *composite);
      if (own.size() == 1U) {
        vertexOf[composite->id] = vertexOf.at(composite->mother);
        continue;
      }
      if (own.front()->vertex == composite->vertex) {
        // a resonance: its products at its parent's vertex, its mass soft
        vertexOf[composite->id] = vertexOf.at(composite->mother);
        resonances.push_back(composite);
        continue;
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
    if (missing && !topology.vertexOf.empty()) {
      topology.vertexOf.push_back(0); // the unseen one, last, at the root
    }
    work.push_back({.mother = &member, .leaves = shape.leaves, .composites = std::move(composites), .topology = std::move(topology), .decay = decay, .unseen = missing ? shape.unseen.front() : nullptr, .resonances = std::move(resonances)});
  }

  const auto& frame = *event.frame;
  const auto& geometry = *event.geometry;
  const auto& material = *event.material;
  o2::gtrk::detail::parallelFor(*event.arena, 0U, work.size(), 4, [&](std::size_t begin, std::size_t end) {
    for (std::size_t index = begin; index < end; ++index) {
      auto& one = work[index];
      const auto count = one.leaves.size();
      std::vector<DecayDaughter> daughters(count);
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
        daughters[i] = {.innermost = refit.innermost, .surface = geometry.planeOf(rows[0]), .propagator = propagator, .mass = species.mass, .chargeMagnitude = std::abs(particleCharge(leaf->pdg))};
        if (species.pid == o2::track::PID::Electron) {
          daughters[i].mixture = refit.innermostMixture;
        }
        masses[i] = species.mass;
      }
      // The unseen one: an unmeasured daughter, last, its start per job.
      double unseenMass = 0.;
      if (one.unseen != nullptr) {
        const auto* entry = o2::O2DatabasePDG::Instance()->GetParticle(one.unseen->pdg);
        unseenMass = entry == nullptr ? 0. : entry->Mass();
        daughters.push_back({.mass = unseenMass, .neutral = true, .unmeasured = true});
        masses.push_back(unseenMass);
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
      if (one.unseen != nullptr) {
        const auto& m = one.unseen->momentum;
        eTrue += std::sqrt((m[0] * m[0]) + (m[1] * m[1]) + (m[2] * m[2]) + (unseenMass * unseenMass));
        for (int k = 0; k < 3; ++k) {
          pTrue[k] += m[k];
        }
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
        base.dPRelInner[i] = static_cast<float>((daughters[i].momentumNumerator() / std::abs(static_cast<double>(daughters[i].innermost.qOverP())) / pFirst) - 1.);
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
      // The fits: every variant alone, then kFull under the mass, the
      // production point, and both.
      struct Job {
        DecayFitVariant variant;
        std::uint8_t constraints;
        int branch = -1;        ///< an unseen neutral's start: the quadratic's smaller (0) or larger (1) root
        bool resonance = false; ///< the resonances' masses taken
      };
      const bool tree = !one.composites.empty();
      // A tree's charged particles' own rows -- the root's (a tracked Xi,
      // Omega), a charged composite's (a tracked Xi in a Xi_c) --, any number,
      // the outermost first, with its mass (DecayVertexFit.h's chains).
      std::vector<DecayMother> trackedMothers;
      if (tree && !one.refitFailed) {
        std::vector<std::pair<const TreeParticle*, int>> charged{{&mother, 0}};
        for (std::size_t w = 0; w < one.composites.size(); ++w) {
          charged.emplace_back(one.composites[w], static_cast<int>(w) + 1);
        }
        for (const auto& [particle, vertex] : charged) {
          const auto* entry = o2::O2DatabasePDG::Instance()->GetParticle(particle->pdg);
          const auto found = chainOf.find(particle->id);
          if (particleCharge(particle->pdg) == 0. || entry == nullptr || found == chainOf.end() || found->second->rows.empty()) {
            continue;
          }
          const auto& own = *found->second;
          if (vertex == 0) {
            base.motherRows = static_cast<std::uint8_t>(std::min<std::size_t>(own.rows.size(), 255U));
          }
          auto propagator = event.propagator;
          propagator.massGeV = static_cast<Scalar>(entry->Mass());
          DecayMother tracked{.vertex = vertex, .propagator = propagator};
          for (auto at = own.rows.rbegin(); at != own.rows.rend(); ++at) {
            const auto m = frame.measurement(*at);
            tracked.planes.push_back(geometry.planeOf(m));
            tracked.positions.push_back({m.loc0, m.loc1});
            tracked.covariance.push_back({m.varLoc0, m.covLoc01, m.varLoc1});
          }
          trackedMothers.push_back(std::move(tracked));
        }
      }
      std::vector<Job> jobs;
      jobs.push_back({DecayFitVariant::kFull, 0U});
      const bool missing = one.unseen != nullptr;
      if (missing) {
        // The visible daughters alone (a tree's composites' masses hold),
        // then the unseen one's momentum under the mother's mass and
        // production point, from either start.
        if (tree) {
          jobs.push_back({DecayFitVariant::kFull, 16U});
        }
        for (const int branch : {0, 1}) {
          jobs.push_back({DecayFitVariant::kFull, static_cast<std::uint8_t>(64U | 3U | (tree ? 16U : 0U)), branch});
        }
      } else if (tree) {
        // A tree: its composites' masses, with the root's mass, its
        // production point, both.
        for (const std::uint8_t set : {std::uint8_t{16}, std::uint8_t{17}, std::uint8_t{18}, std::uint8_t{19}}) {
          jobs.push_back({DecayFitVariant::kFull, set});
        }
        // A tracked mother: its own rows in the tree.
        if (!trackedMothers.empty()) {
          for (const std::uint8_t set : {std::uint8_t{48}, std::uint8_t{50}, std::uint8_t{51}}) {
            jobs.push_back({DecayFitVariant::kFull, set});
          }
        }
      } else {
        for (const std::uint8_t set : {std::uint8_t{1}, std::uint8_t{2}, std::uint8_t{3}}) {
          jobs.push_back({DecayFitVariant::kFull, set});
        }
      }
      // Resonances' masses: alone, and with the mother's mass and
      // production point (and a tree's composites' masses).
      if (!one.resonances.empty() && !missing) {
        for (const std::uint8_t set : {std::uint8_t{0}, static_cast<std::uint8_t>(tree ? 19U : 3U)}) {
          jobs.push_back({.variant = DecayFitVariant::kFull, .constraints = set, .resonance = true});
        }
      }
      // Electron daughters' mixtures kept: alone and under the masses.
      const bool electrons = std::ranges::any_of(daughters, [](const DecayDaughter& d) { return d.mixture.has_value(); });
      if (electrons) {
        for (const std::uint8_t set : {std::uint8_t{0}, tree ? std::uint8_t{16} : std::uint8_t{1}}) {
          if (!missing || set != 1U) {
            jobs.push_back({DecayFitVariant::kGaussianSum, set});
          }
        }
      }
      // A tree's composites: their PDG masses, the first one's truth.
      std::vector<double> compositeMasses;
      for (const auto* composite : one.composites) {
        const auto* entry = o2::O2DatabasePDG::Instance()->GetParticle(composite->pdg);
        compositeMasses.push_back(entry == nullptr ? 0. : entry->Mass());
      }
      // each vertex's particle's mass, for its flight's material in the fit
      one.topology.massGeV.assign(1, motherMass);
      one.topology.massGeV.insert(one.topology.massGeV.end(), compositeMasses.begin(), compositeMasses.end());
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
      DecayVertexFitResult free; // kFull without constraints: the missing momentum's start
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
        const bool wantComposites = (job.constraints & 16U) != 0U;
        if ((wantMass && !(motherMass > 0.)) || (wantProduction && !madeAtCollision)) {
          record.status = static_cast<std::uint8_t>(DecayFitStatus::kNoStart); // nothing to constrain to
          one.records.push_back(record);
          continue;
        }
        DecayConstraints constraints{.mass = wantMass, .massGeV = motherMass, .production = wantProduction, .point = point};
        if (wantComposites) {
          constraints.compositeMasses = compositeMasses;
        }
        record.resonance = job.resonance ? 1U : 0U;
        if (job.resonance) {
          for (const auto* resonance : one.resonances) {
            const auto* entry = o2::O2DatabasePDG::Instance()->GetParticle(resonance->pdg);
            Resonance taken{.massGeV = entry == nullptr ? 0. : entry->Mass(), .widthGeV = entry == nullptr ? 0. : entry->Width()};
            for (std::size_t i = 0; i < count; ++i) {
              for (const auto* up = one.leaves[i]; up != nullptr; up = trees.find(up->mother)) {
                if (up->id == resonance->id) {
                  taken.daughters.push_back(static_cast<int>(i));
                  break;
                }
              }
            }
            if (taken.massGeV > 0. && taken.widthGeV > 0. && taken.daughters.size() >= 2U) {
              constraints.resonances.push_back(std::move(taken));
            }
          }
        }
        const bool wantMother = (job.constraints & 32U) != 0U;
        const auto tracked = wantMother ? std::span<const DecayMother>{trackedMothers} : std::span<const DecayMother>{};
        // The unseen one with the visible ones only under the missing fit,
        // from the visible fit's starts (`missingMomentumStarts`).
        record.branch = static_cast<std::int8_t>(job.branch);
        auto jobDaughters = daughters;
        const bool wantMissing = (job.constraints & 64U) != 0U;
        if (wantMissing) {
          std::optional<std::array<double, 3>> start;
          if (free.status == DecayFitStatus::kOk) {
            const auto starts = missingMomentumStarts(free, masses, point, motherMass, unseenMass);
            if (static_cast<std::size_t>(job.branch) < starts.size()) {
              start = starts[static_cast<std::size_t>(job.branch)];
            }
          }
          const double unseenP = start ? std::hypot((*start)[0], (*start)[1], (*start)[2]) : 0.;
          if (!(unseenP > 0.)) {
            record.status = static_cast<std::uint8_t>(DecayFitStatus::kNoStart);
            one.records.push_back(record);
            continue;
          }
          auto& unseen = jobDaughters.back();
          unseen.innermost.par[kPhi] = static_cast<Scalar>(std::atan2((*start)[1], (*start)[0]));
          unseen.innermost.par[kTheta] = static_cast<Scalar>(std::atan2(std::hypot((*start)[0], (*start)[1]), (*start)[2]));
          unseen.innermost.par[kQOverP] = static_cast<Scalar>(1. / unseenP);
        } else if (missing) {
          jobDaughters.pop_back();
        }
        const auto fitUnder = [&](const DecayConstraints& set) { return fitDecayVertex(std::span<const DecayDaughter>{jobDaughters}, material, job.variant, set, one.topology, tracked); };
        const auto started = std::chrono::steady_clock::now();
        const auto outcome = fitUnder(constraints);
        record.microseconds = static_cast<float>(std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - started).count());
        if (job.constraints == 0U && job.variant == DecayFitVariant::kFull) {
          free = outcome;
        }
        const auto& fit = outcome;
        record.status = static_cast<std::uint8_t>(fit.status);
        record.iterations = fit.iterations;
        record.stopReason = fit.stopReason;
        record.outer = fit.outer;
        record.modelShift = fit.modelShift;
        record.pointSigma = static_cast<float>(fit.pointSigma);
        record.combinations = fit.combinations;
        record.stopMetric = fit.stopMetric;
        if (fit.status != DecayFitStatus::kOk) {
          one.records.push_back(record);
          continue;
        }
        const int vertexCount = static_cast<int>(fit.vertices.size());
        const std::size_t fitCount = fit.momenta.size();
        const int n = (3 * vertexCount) + (3 * static_cast<int>(fitCount));
        const auto& c = fit.covariance;
        const auto& fitMasses = masses;
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
        const auto fitted = kinematicsOf(fit.momenta, fitMasses, fit.charges, c, n, 3 * vertexCount);
        record.mass = static_cast<float>(fitted.mass);
        record.sMass = static_cast<float>(fitted.sMass);
        record.p = static_cast<float>(fitted.p);
        record.sP = static_cast<float>(fitted.sP);
        if (hasCollision) {
          const auto motherFit = motherQuantities(fit.vertex, fit.momenta, fit.charges, point, motherMass, c, vertexCount);
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
          const auto& child = fit.vertices[1];
          const int childN = n;
          const int childAt = 3;
          const auto& childC = c;
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
          // its variance the two vertices' with their correlation.
          record.childLength = static_cast<float>(((child[0] - fit.vertex[0]) * childAxis[0]) + ((child[1] - fit.vertex[1]) * childAxis[1]) + ((child[2] - fit.vertex[2]) * childAxis[2]));
          double lengthVariance = childVariance(childAxis);
          for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
              lengthVariance += childAxis[a] * c[(a * n) + b] * childAxis[b];
              lengthVariance -= 2. * childAxis[a] * c[(a * n) + 3 + b] * childAxis[b];
            }
          }
          record.sChildLength = static_cast<float>(std::sqrt(std::max(lengthVariance, 0.)));
          {
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
            std::vector<double> childCharges;
            for (const auto i : childLeaves) {
              childCharges.push_back(fit.charges[i]);
            }
            const auto childKinematics = kinematicsOf(childMomenta, childMasses, childCharges, block, m, 0);
            record.childMass = static_cast<float>(childKinematics.mass);
            record.sChildMass = static_cast<float>(childKinematics.sMass);
          }
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
        const auto measured = kinematicsOf(fit.carried, masses, fit.charges, carried, stride, 0);
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
