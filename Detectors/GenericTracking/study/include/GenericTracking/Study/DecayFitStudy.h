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
/// `kMaxDaughters` of them, every one findable -- its rows on a chain of at
/// least four linked layers (`isFindable`, no pT floor): a track of three is
/// barely defined -- and with a chain from the source; no
/// photon's conversion and no electron daughter (their fit is the
/// Gaussian-sum filter's, a later step). Each daughter's rows are refitted
/// with its true species (`refitWithHypothesis`); the vertex is fitted under
/// every `DecayFitVariant`, a record each: from the daughters' states, and
/// from their rows (`kRows`, DecayRowsFit.h, the information limit).
///
/// Truth: the vertex against the MC decay point (the products' production
/// point), in the global frame and along the true flight; the daughters'
/// momenta against their MC momenta at production; the invariant mass of the
/// refitted momenta and of the carried ones (without the vertex's
/// information) against the true daughters' (with the species' masses), the
/// summed momentum against theirs.

#ifndef GTRK_STUDY_DECAYFITSTUDY_H_
#define GTRK_STUDY_DECAYFITSTUDY_H_

#include "GenericTracking/Study/DecayCandidates.h"
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
  std::uint8_t variant = 0U; ///< DecayFitVariant
  std::uint8_t daughters = 0U;
  std::uint8_t status = 0U; ///< DecayFitStatus, or kRefitFailed
  std::uint8_t iterations = 0U;
  float chi2 = 0.F;
  std::int16_t ndf = 0;
  float flight = 0.F;        ///< the true one [cm]
  float decayR = 0.F;        ///< the true decay point's distance from the z axis [cm]
  float pTrue = 0.F;         ///< the true daughters' summed |p| [GeV]
  float minDaughterPt = 0.F; ///< the lowest true daughter pT [GeV]
  float maxCarryX2X0 = 0.F;  ///< the most material a daughter's last carry crossed
  std::uint8_t minDaughterRows = 0U;
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
  /// The first daughters' refitted phi, theta, q/p less the true ones, and
  /// their sigmas.
  std::array<float, kDaughterColumns> dPhi{};
  std::array<float, kDaughterColumns> dTheta{};
  std::array<float, kDaughterColumns> dQOverP{};
  std::array<float, kDaughterColumns> sPhi{};
  std::array<float, kDaughterColumns> sTheta{};
  std::array<float, kDaughterColumns> sQOverP{};

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("source", source);
    visit("particle", particle);
    visit("pdg", pdg);
    visit("variant", variant);
    visit("daughters", daughters);
    visit("status", status);
    visit("iterations", iterations);
    visit("chi2", chi2);
    visit("ndf", ndf);
    visit("flight", flight);
    visit("decayR", decayR);
    visit("pTrue", pTrue);
    visit("minDaughterPt", minDaughterPt);
    visit("maxCarryX2X0", maxCarryX2X0);
    visit("minDaughterRows", minDaughterRows);
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
    static constexpr std::array<std::string_view, kDaughterColumns> kSuffix{"0", "1", "2"};
    for (std::size_t d = 0; d < kDaughterColumns; ++d) {
      // The names must outlive the writer's first call: kept in static storage.
      static const std::array<std::array<std::string, 6>, kDaughterColumns> kNames = [] {
        std::array<std::array<std::string, 6>, kDaughterColumns> names;
        for (std::size_t k = 0; k < kDaughterColumns; ++k) {
          const std::string s{kSuffix[k]};
          names[k] = {"dPhi" + s, "dTheta" + s, "dQOverP" + s, "sPhi" + s, "sTheta" + s, "sQOverP" + s};
        }
        return names;
      }();
      visit(kNames[d][0], dPhi[d]);
      visit(kNames[d][1], dTheta[d]);
      visit(kNames[d][2], dQOverP[d]);
      visit(kNames[d][3], sPhi[d]);
      visit(kNames[d][4], sTheta[d]);
      visit(kNames[d][5], sQOverP[d]);
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

  // A daughter must be findable: on a chain of four linked layers, at any pT.
  const auto rule = comparisonSettings<Det>(4, 0.);
  const auto layers = particleLayers(event.rows, rule);
  const auto findable = [&](const TreeParticle& leaf) {
    const auto found = layers.find(leaf.id);
    return found != layers.end() && isFindable(std::numeric_limits<float>::max(), found->second, rule);
  };

  struct Work {
    const TreeParticle* mother = nullptr;
    std::vector<const TreeParticle*> leaves;
    std::vector<DecayFitRecord> records;
    bool refitFailed = false;
    std::uint8_t refitPass = 0U;
    std::uint8_t refitStatus = 0U;
  };
  std::vector<Work> work;
  for (const auto& member : trees.particles) {
    if (member.pdg == 22 || isFinalCharged(member.pdg)) {
      continue;
    }
    const auto products = decayProducts(trees, member);
    if (products.empty()) {
      continue;
    }
    const auto& decay = products.front()->vertex;
    if (decay == member.vertex) {
      continue; // no vertex of its own: a resonance
    }
    const auto shape = decayShape(trees, member, particleCharge);
    if (!shape.closed() || !shape.composites.empty() || shape.leaves.size() < 2 || shape.leaves.size() > kMaxDaughters) {
      continue;
    }
    const bool usable = std::ranges::all_of(shape.leaves, [&](const TreeParticle* leaf) {
      return std::abs(leaf->pdg) != 11 && truthSpecies(leaf->pdg).has_value() && chainOf.contains(leaf->id);
    });
    if (!usable) {
      continue;
    }
    if (!std::ranges::all_of(shape.leaves, [&](const TreeParticle* leaf) { return findable(*leaf); })) {
      ++mUnfindable;
      continue;
    }
    work.push_back({.mother = &member, .leaves = shape.leaves});
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
        auto& own = rowsDaughters[i];
        own.propagator = propagator;
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
      const auto& decay = one.leaves.front()->vertex; // a leaf of a one-vertex decay is a product
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
      std::array<DecayVertexFitResult, kDecayFitVariants> fits{};
      for (std::size_t v = 0; v < kDecayFitVariants; ++v) {
        DecayFitRecord record = base;
        record.variant = static_cast<std::uint8_t>(v);
        if (one.refitFailed) {
          record.status = DecayFitRecord::kRefitFailed;
          one.records.push_back(record);
          continue;
        }
        const auto variant = static_cast<DecayFitVariant>(v);
        fits[v] = variant == DecayFitVariant::kRows ? fitDecayRows(std::span<const RowsDaughter>{rowsDaughters}, material, fits[static_cast<std::size_t>(DecayFitVariant::kFull)]) : fitDecayVertex(std::span<const DecayDaughter>{daughters}, material, variant);
        const auto& fit = fits[v];
        record.status = static_cast<std::uint8_t>(fit.status);
        record.iterations = fit.iterations;
        if (fit.status != DecayFitStatus::kOk) {
          one.records.push_back(record);
          continue;
        }
        const int n = 3 + (3 * static_cast<int>(count));
        const auto& c = fit.covariance;
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
        const auto fitted = kinematicsOf(fit.momenta, masses, c, n, 3);
        record.mass = static_cast<float>(fitted.mass);
        record.sMass = static_cast<float>(fitted.sMass);
        record.p = static_cast<float>(fitted.p);
        record.sP = static_cast<float>(fitted.sP);
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
          const int q = 3 + (3 * static_cast<int>(i));
          record.dPhi[i] = static_cast<float>(std::remainder(fit.momenta[i][0] - qTrue[i][0], 2. * std::numbers::pi));
          record.dTheta[i] = static_cast<float>(fit.momenta[i][1] - qTrue[i][1]);
          record.dQOverP[i] = static_cast<float>(fit.momenta[i][2] - qTrue[i][2]);
          record.sPhi[i] = static_cast<float>(std::sqrt(c[q * n + q]));
          record.sTheta[i] = static_cast<float>(std::sqrt(c[(q + 1) * n + (q + 1)]));
          record.sQOverP[i] = static_cast<float>(std::sqrt(c[(q + 2) * n + (q + 2)]));
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
