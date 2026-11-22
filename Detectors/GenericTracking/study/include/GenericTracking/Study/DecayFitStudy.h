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
/// vertex, and both (DecayFitConstraints.h), a record each.
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

#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Tracking/DaughterTable.h"
#include "GenericTracking/Tracking/DecayCandidateFit.h"
#include "GenericTracking/Tracking/DecayFitBatch.h"
#include "GenericTracking/Tracking/DecayFitConstraints.h"
#include "GenericTracking/Tracking/DecayObservables.h"
#include "GenericTracking/Tracking/DecayVertexFit.h"
#include "GenericTracking/Tracking/GPU/DeviceContext.h"
#include "GenericTracking/Tracking/GPU/DeviceDecayFitter.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"
#include "GenericTracking/Tracking/HypothesisRefit.h"
#include "GenericTracking/Tracking/Parallel.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Study/DecayCandidates.h"
#include "GenericTracking/Study/DecayStudyTools.h"
#include "GenericTracking/Study/DecayTrees.h"
#include "GenericTracking/Study/RecordWriter.h"
#include "GenericTracking/Study/StateTruth.h"
#include "GenericTracking/Study/StudyChains.h"
#include "GenericTracking/Study/TrackComparison.h"
#include "GenericTracking/Study/TruthEvent.h"
#include "GenericTracking/Study/TruthFit.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <memory>
#include <optional>
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
  std::uint8_t status = 0U; ///< VertexFitStatus, or kRefitFailed
  std::uint8_t iterations = 0U;
  std::uint8_t outer = 0U;
  std::uint8_t stopReason = 0U;
  float stopMetric = 0.F;
  float stopBound = 0.F; ///< what the last step's gain was held to (the inputs' resolution R)
  float modelShift = 0.F;
  float microseconds = 0.F;             ///< the fit's time
  std::uint16_t combinations = 0U;      ///< the Gaussian sum's component combinations fitted and taken
  std::uint16_t combinationsTried = 0U; ///< and all it tried
  std::int8_t branch = -1;              ///< an unseen neutral's start (0, 1), -1 without one
  std::uint8_t resonance = 0U;          ///< 1 the resonances' masses taken (Breit-Wigner)
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
  std::uint8_t motherRows = 0U;        ///< a charged root's own rows fitted (a tracked Xi or Omega), 0 without
  std::uint8_t motherTrimmed = 0U;     ///< its rows left out, one row per layer kept
  std::uint8_t compositesTrimmed = 0U; ///< the charged composites' rows left out, one row per layer kept
  float pointSigma = 0.F;              ///< the production point's equations' sigma across the flight [cm] (0: held exactly)
  std::uint8_t refitPass = 0U;         ///< a failed daughter refit's FitPass
  std::uint8_t refitStatus = 0U;       ///< and its PropagationStatus
  float dX = 0.F;                      ///< fitted less true vertex [cm]
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
    visit("combinationsTried", combinationsTried);
    visit("microseconds", microseconds);
    visit("stopReason", stopReason);
    visit("stopMetric", stopMetric);
    visit("stopBound", stopBound);
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
    visit("motherTrimmed", motherTrimmed);
    visit("compositesTrimmed", compositesTrimmed);
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

/// One decay fit in reduced precision against the same fit in double
/// (`DecayFitStudy` as `decayprecision`): the model, the solve or both in
/// float (`DecayVertexFitter`'s scalars), judged in the double fit's own
/// terms -- its chi2 under its last model at the float fit's parameters less
/// at its own, against its stopping bound R.
struct DecayPrecisionRecord {
  static constexpr std::string_view kName = "decayPrecision";
  std::uint32_t timeframe = 0U;
  std::uint8_t source = 0U; ///< ChainSource
  std::uint64_t particle = 0U;
  std::int32_t pdg = 0;
  std::uint8_t variant = 0U;     ///< DecayFitVariant
  std::uint8_t constraints = 0U; ///< as DecayFitRecord's
  std::int8_t branch = -1;
  std::uint8_t resonance = 0U;
  std::uint8_t daughters = 0U;
  std::uint8_t configuration = 0U; ///< 0 the double fit itself, 1 the model in float, 2 the solve in float, 3 both
  std::uint8_t statusDouble = 0U;  ///< VertexFitStatus of the double fit
  std::uint8_t status = 0U;        ///< and of this one
  std::uint8_t stopReasonDouble = 0U;
  std::uint8_t iterations = 0U;
  double stopMetricDouble = 0.;
  double stopBoundDouble = 0.; ///< R
  /// The double fit's chi2 under its last model at this fit's parameters,
  /// less at its own; NaN unless both converged.
  double delta = std::numeric_limits<double>::quiet_NaN();
  double violationDouble = std::numeric_limits<double>::quiet_NaN(); ///< |c|_1 in double at the double fit's parameters
  double violation = std::numeric_limits<double>::quiet_NaN();       ///< and at this one's
  float dVertex = std::numeric_limits<float>::quiet_NaN();           ///< the root vertex's shift from the double fit's [um]
  float dMass = std::numeric_limits<float>::quiet_NaN();             ///< the fitted mass's [GeV]
  float sTRel = std::numeric_limits<float>::quiet_NaN();             ///< sT's (across the true flight, transverse) relative difference
  float chi2Rel = std::numeric_limits<float>::quiet_NaN();           ///< the fit's own chi2's relative difference
  float microseconds = 0.F;

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("source", source);
    visit("particle", particle);
    visit("pdg", pdg);
    visit("variant", variant);
    visit("constraints", constraints);
    visit("branch", branch);
    visit("resonance", resonance);
    visit("daughters", daughters);
    visit("configuration", configuration);
    visit("statusDouble", statusDouble);
    visit("status", status);
    visit("stopReasonDouble", stopReasonDouble);
    visit("iterations", iterations);
    visit("stopMetricDouble", stopMetricDouble);
    visit("stopBoundDouble", stopBoundDouble);
    visit("delta", delta);
    visit("violationDouble", violationDouble);
    visit("violation", violation);
    visit("dVertex", dVertex);
    visit("dMass", dMass);
    visit("sTRel", sTRel);
    visit("chi2Rel", chi2Rel);
    visit("microseconds", microseconds);
  }
};

/// `decaydevice`'s record: one fit, made directly (`fitDecayVertex`),
/// through the packed batch on the host (`fitDecayBatch`) and, with a GPU
/// backend, on the device (DeviceDecayFitter.h).
struct DecayDeviceRecord {
  static constexpr std::string_view kName = "decayDevice";
  static constexpr std::uint8_t kNone = 255U; ///< `statusDevice` without a device
  std::uint32_t timeframe = 0U;
  std::uint8_t source = 0U; ///< ChainSource
  std::uint64_t particle = 0U;
  std::int32_t pdg = 0;
  std::uint8_t variant = 0U;     ///< DecayFitVariant
  std::uint8_t constraints = 0U; ///< as DecayFitRecord's
  std::int8_t branch = -1;
  std::uint8_t resonance = 0U;
  std::uint8_t daughters = 0U;
  std::uint8_t statusDirect = 0U; ///< VertexFitStatus
  std::uint8_t statusHost = 0U;
  std::uint8_t statusDevice = kNone;
  std::uint8_t hostSame = 0U; ///< 1: the batch's host fit equals the direct one to the bit (statuses, vertices, momenta, chi2, covariance)
  /// The device fit against the direct one, both converged: the root
  /// vertex's shift [um], the mass's difference [GeV], chi2's relative
  /// difference, the largest relative difference of a covariance diagonal.
  float dVertex = std::numeric_limits<float>::quiet_NaN();
  float dMass = std::numeric_limits<float>::quiet_NaN();
  float chi2Rel = std::numeric_limits<float>::quiet_NaN();
  float diagonalRel = std::numeric_limits<float>::quiet_NaN();

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("source", source);
    visit("particle", particle);
    visit("pdg", pdg);
    visit("variant", variant);
    visit("constraints", constraints);
    visit("branch", branch);
    visit("resonance", resonance);
    visit("daughters", daughters);
    visit("statusDirect", statusDirect);
    visit("statusHost", statusHost);
    visit("statusDevice", statusDevice);
    visit("hostSame", hostSame);
    visit("dVertex", dVertex);
    visit("dMass", dMass);
    visit("chi2Rel", chi2Rel);
    visit("diagonalRel", diagonalRel);
  }
};

namespace decay_fit_study_detail
{

using decay_fit_detail::momentumOf;
using decay_fit_detail::numeratorOf;

using decay_fit_detail::quadraticForm;

/// \a chain's rows, one per layer: on a layer holding more (overlapping
/// sensors, split clusters), the row of the smallest global radius, the
/// smaller row ID on a tie; the kept rows in the chain's order.
[[nodiscard]] inline std::vector<std::uint32_t> onePerLayer(std::span<const TruthRow> rows, std::span<const std::uint32_t> chain)
{
  const auto better = [&](std::uint32_t a, std::uint32_t b) {
    return rows[a].r < rows[b].r || (rows[a].r == rows[b].r && a < b);
  };
  std::vector<std::uint32_t> kept;
  for (const auto row : chain) {
    const bool beaten = std::ranges::any_of(chain, [&](std::uint32_t other) {
      return other != row && rows[other].layer == rows[row].layer && better(other, row);
    });
    if (!beaten) {
      kept.push_back(row);
    }
  }
  return kept;
}

/// The fit of \a daughters in double and in the three reduced-precision
/// configurations (`DecayPrecisionRecord`), a record each from \a base into
/// \a out; \a masses the daughters' (the fitted mass), \a across the true
/// flight's transverse normal (sT). Returns the double fit's result.
template <DetectorLayout Det, DecayFitVariant Variant, typename Material>
DecayVertexFitResult measurePrecision(std::span<const DecayDaughter> daughters, const Material& material, const DecayConstraints& constraints, const DecayTopology& topology, std::span<const DecayMother> mothers, std::span<const double> masses, const Vec3D& across,
                                      const DecayPrecisionRecord& base, std::vector<DecayPrecisionRecord>& out)
{
  const auto massOf = [&](const DecayVertexFitResult& fit) { return fitKinematics(fit, masses).mass; };
  const auto sigmaT = [&](const DecayVertexFitResult& fit) { return sqrtClampedOf(quadraticForm(across, across, 3, [&](int a, int b) { return fit.covarianceAt(a, b); })); };
  std::array<Vec3D, decay_fit_detail::kMaxConstrainedVertices> starts{};
  const std::span<const Vec3D> start{starts.data(), decayFitStart<Det>(daughters, constraints, topology, mothers, starts) ? static_cast<std::size_t>(topology.vertices()) : 0U};
  using Reference = DecayVertexFitter<Det, Variant, Material, double, double>;
  Reference reference(daughters, material, constraints, topology, mothers);
  WorkspaceBuffer workspace; // the reference's, read again by its frozenMerit
  auto started = std::chrono::steady_clock::now();
  DecayVertexFitResult fit;
  reference.fit(start, workspace.take(reference.workspace()), fit);
  auto record = base;
  record.microseconds = static_cast<float>(std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - started).count());
  record.statusDouble = record.status = static_cast<std::uint8_t>(fit.status);
  record.stopReasonDouble = fit.stopReason;
  record.iterations = fit.iterations;
  record.stopMetricDouble = fit.stopMetric;
  record.stopBoundDouble = fit.stopBound;
  const bool ok = fit.status == VertexFitStatus::kOk;
  const std::vector<double> solution(reference.solution().begin(), reference.solution().end());
  double chi2 = 0.;
  double violation = 0.;
  const bool merit = ok && reference.frozenMerit(solution, chi2, violation);
  if (merit) {
    record.delta = 0.;
    record.violationDouble = record.violation = violation;
    record.dVertex = record.dMass = record.sTRel = record.chi2Rel = 0.F;
  }
  out.push_back(record);
  const double mass = ok ? massOf(fit) : 0.;
  const double sT = ok ? sigmaT(fit) : 0.;
  const auto reduced = [&]<typename Model, typename Solve>(std::uint8_t configuration) {
    DecayVertexFitter<Det, Variant, Material, Model, Solve> fitter(daughters, material, constraints, topology, mothers);
    WorkspaceBuffer fitterWorkspace;
    const auto begin = std::chrono::steady_clock::now();
    DecayVertexFitResult result;
    fitter.fit(start, fitterWorkspace.take(fitter.workspace()), result);
    auto one = record;
    one.configuration = configuration;
    one.microseconds = static_cast<float>(std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - begin).count());
    one.status = static_cast<std::uint8_t>(result.status);
    one.iterations = result.iterations;
    one.delta = one.violation = std::numeric_limits<double>::quiet_NaN();
    one.dVertex = one.dMass = one.sTRel = one.chi2Rel = std::numeric_limits<float>::quiet_NaN();
    if (merit && result.status == VertexFitStatus::kOk) {
      const auto at = fitter.solution();
      double atChi2 = 0.;
      double atViolation = 0.;
      if (reference.frozenMerit(at, atChi2, atViolation)) {
        one.delta = atChi2 - chi2;
        one.violation = atViolation;
      }
      one.dVertex = static_cast<float>(1e4 * std::hypot(result.vertex[0] - fit.vertex[0], result.vertex[1] - fit.vertex[1], result.vertex[2] - fit.vertex[2]));
      one.dMass = static_cast<float>(massOf(result) - mass);
      one.sTRel = static_cast<float>((sigmaT(result) / sT) - 1.);
      one.chi2Rel = static_cast<float>((result.chi2 / fit.chi2) - 1.);
    }
    out.push_back(one);
  };
  reduced.template operator()<float, double>(1U);
  reduced.template operator()<double, float>(2U);
  reduced.template operator()<float, float>(3U);
  return fit;
}

/// What `decaydevice` compares of a fit: its status, vertices, momenta,
/// chi2, covariance diagonal, and the mass of its momenta with \a
/// daughters' masses.
struct FitSummary {
  std::uint8_t status = 0U;
  FixedVector<Vec3D, DecayVertexFitResult::kVertices> vertices;
  FixedVector<Vec3D, DecayVertexFitResult::kDaughters> momenta;
  FixedVector<double, (3 * DecayVertexFitResult::kVertices) + (3 * DecayVertexFitResult::kDaughters)> diagonal;
  double chi2 = 0.;
  double mass = 0.;
};

[[nodiscard]] inline FitSummary summarise(const DecayVertexFitResult& fit, std::span<const DecayDaughter> daughters)
{
  FitSummary summary{.status = static_cast<std::uint8_t>(fit.status), .vertices = fit.vertices, .momenta = fit.momenta, .chi2 = fit.chi2};
  if (fit.status != VertexFitStatus::kOk) {
    return summary;
  }
  const int n = (3 * static_cast<int>(fit.vertices.size())) + (3 * static_cast<int>(fit.momenta.size()));
  for (int a = 0; a < n; ++a) {
    summary.diagonal.push_back(fit.covarianceAt(a, a));
  }
  std::vector<double> masses;
  masses.reserve(fit.momenta.size());
  for (std::size_t i = 0; i < fit.momenta.size(); ++i) {
    masses.push_back(daughters[i].mass);
  }
  summary.mass = fitKinematics(fit, masses).mass;
  return summary;
}

[[nodiscard]] inline bool sameBits(const FitSummary& a, const FitSummary& b)
{
  if (a.status != b.status || a.vertices.size() != b.vertices.size() || a.momenta.size() != b.momenta.size() || a.diagonal.size() != b.diagonal.size() || !o2::gtrk::sameBits(a.chi2, b.chi2) || !o2::gtrk::sameBits(a.mass, b.mass)) {
    return false;
  }
  for (std::size_t w = 0; w < a.vertices.size(); ++w) {
    for (int c = 0; c < 3; ++c) {
      if (!o2::gtrk::sameBits(a.vertices[w][c], b.vertices[w][c])) {
        return false;
      }
    }
  }
  for (std::size_t i = 0; i < a.momenta.size(); ++i) {
    for (int c = 0; c < 3; ++c) {
      if (!o2::gtrk::sameBits(a.momenta[i][c], b.momenta[i][c])) {
        return false;
      }
    }
  }
  for (std::size_t j = 0; j < a.diagonal.size(); ++j) {
    if (!o2::gtrk::sameBits(a.diagonal[j], b.diagonal[j])) {
      return false;
    }
  }
  return true;
}

} // namespace decay_fit_study_detail

template <DetectorLayout Det>
class DecayFitStudy
{
 public:
  static constexpr std::size_t kMaxDaughters = 5;

  /// What the study writes: the fits (`decayfit`, DecayFitRecord); every fit
  /// also in reduced precision (`decayprecision`, DecayPrecisionRecord); every
  /// fit also through a packed batch on the host and, with a GPU backend, on
  /// the device (`decaydevice`, DecayDeviceRecord).
  enum class Mode : std::uint8_t {
    kFit,
    kPrecision,
    kDevice
  };

  explicit DecayFitStudy(RecordWriter& writer, Mode mode = Mode::kFit, const TrackingDevice& device = {}) : mWriter(&writer), mPrecision(mode == Mode::kPrecision), mMode(mode), mDeviceSettings(device)
  {
    if (mode == Mode::kPrecision) {
      writer.begin<DecayPrecisionRecord>();
    } else if (mode == Mode::kDevice) {
      writer.begin<DecayDeviceRecord>();
    } else {
      writer.begin<DecayFitRecord>();
    }
    if (mode == Mode::kDevice && !device.backend.empty()) {
      const std::array<std::size_t, 1> blocks{device.arenaMB * 1024U * 1024U};
      mContext = std::make_unique<DeviceContext>(device.backend, device.device, blocks);
      mDevice = mContext->runtime().template makeDecayFitter<Det>(mContext->block(0));
    }
  }

  /// \a event's single-vertex decays, from \a chains of \a source.
  void process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source);

  [[nodiscard]] std::vector<std::string> summary() const
  {
    std::vector<std::string> lines{std::format("decay fit study: {} decays, {} records into {}; {} with a failed daughter refit; {} left out for a daughter not findable", mDecays, mRecords, mWriter->path(), mRefitFailed, mUnfindable)};
    if (mMode == Mode::kDevice) {
      const auto perEvent = [&](double seconds) { return mEvents > 0U ? 1e3 * seconds / static_cast<double>(mEvents) : 0.; };
      lines.push_back(std::format("decay device: {} fits of {} events in {} timeframes; the batch on the host (one thread) {:.4f} ms per event ({:.4f} s); equal to the direct fit to the bit: {} of {}",
                                  mBatchFits, mEvents, mTimeframes, perEvent(mHostSeconds), mHostSeconds, mHostSame, mBatchFits));
      if (mDevice) {
        lines.push_back(std::format("decay device on {}: {:.4f} ms per event ({:.4f} s: upload {:.4f}, sizes {:.4f}, fits {:.4f} in {} launches, download {:.4f}); statuses equal to the direct fit's: {} of {}; {} refused (workspace); arena peak {} B",
                                    mDeviceSettings.backend, perEvent(mDeviceTotal.uploadSeconds + mDeviceTotal.sizeSeconds + mDeviceTotal.fitSeconds + mDeviceTotal.downloadSeconds), mDeviceTotal.uploadSeconds + mDeviceTotal.sizeSeconds + mDeviceTotal.fitSeconds + mDeviceTotal.downloadSeconds,
                                    mDeviceTotal.uploadSeconds, mDeviceTotal.sizeSeconds, mDeviceTotal.fitSeconds, mDeviceTotal.launches, mDeviceTotal.downloadSeconds, mDeviceStatusSame, mBatchFits, mDeviceTotal.refused, mDeviceTotal.peakBytes));
      }
    }
    return lines;
  }

 private:
  using Device = DeviceDecayFitter<Det, LinkMaterial<Det>>;
  /// `decaydevice`: the batches of one source's timeframe fitted on the host
  /// and on the device, compared with the direct fits, written.
  void fitBatches(const TruthEvent<Det>& event, std::array<DecayFitBatch, 2>& batches, std::array<std::vector<DecayDeviceRecord>, 2>& records, std::array<std::vector<decay_fit_study_detail::FitSummary>, 2>& direct);

  RecordWriter* mWriter = nullptr;
  bool mPrecision = false;
  Mode mMode = Mode::kFit;
  TrackingDevice mDeviceSettings; ///< where `decaydevice` fits beside the host batch
  std::unique_ptr<DeviceContext> mContext;
  std::unique_ptr<Device> mDevice;
  const LinkMaterial<Det>* mConfigured = nullptr; ///< the material the device has
  std::size_t mTimeframes = 0U;
  std::size_t mEvents = 0U;
  std::size_t mBatchFits = 0U;
  std::size_t mHostSame = 0U;
  std::size_t mDeviceStatusSame = 0U;
  double mHostSeconds = 0.;
  typename Device::Timing mDeviceTotal;
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
  const auto chainOf = longestChains(chains);

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
    Vec3D decay{};                               ///< the root's true decay point
    const TreeParticle* unseen = nullptr;        ///< the root's one unseen neutral product, if any
    std::vector<const TreeParticle*> resonances; ///< composites that decay where they are made: transparent
    std::vector<DecayFitRecord> records;
    std::vector<DecayPrecisionRecord> precision; ///< `decayprecision`'s
    /// `decaydevice`'s: per variant the fits' inputs, their records' identity
    /// and the direct fits.
    std::array<DecayFitBatch, 2> batches;
    std::array<std::vector<DecayDeviceRecord>, 2> device;
    std::array<std::vector<FitSummary>, 2> direct;
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
    if (composites.size() + 1U > static_cast<std::size_t>(decay_fit_detail::kMaxConstrainedVertices)) {
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
  constexpr o2::gtrk::detail::LayerMaterial material{};
  // The tracker's tracks: each daughter's (track, species) refit from the
  // timeframe's daughter table, once whatever the decays that take it.
  std::optional<DaughterTable<Det>> table;
  std::optional<DaughterSource<Det>> tableSource;
  if (source == ChainSource::kTracks) {
    const std::array<o2::track::PID, 5> species{o2::track::PID{o2::track::PID::Electron}, o2::track::PID{o2::track::PID::Muon}, o2::track::PID{o2::track::PID::Pion}, o2::track::PID{o2::track::PID::Kaon}, o2::track::PID{o2::track::PID::Proton}};
    table.emplace(std::span<const o2::track::PID>{species});
    fillDaughterTable(*table, event, "decay fit study", [&](DaughterTable<Det>& asked) {
      for (const auto& one : work) {
        for (const auto* leaf : one.leaves) {
          const auto& context = chainOf.at(leaf->id)->track;
          const auto species = truthSpecies(leaf->pdg);
          if (!context || !species) {
            throw std::logic_error("decay fit study: a track chain without its track, or a leaf without its species");
          }
          asked.request(context->track, static_cast<std::uint8_t>(asked.speciesIndex(species->pid)));
        }
      }
    });
    tableSource = daughterSource(*table, event);
  }
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
        // The tracks' entry in the table, or the truth chain refitted.
        HypothesisFit<> refit;
        std::uint32_t entry = 0U;
        if (table) {
          entry = static_cast<std::uint32_t>(table->entryOf(chain.track->track, static_cast<std::uint8_t>(table->speciesIndex(species.pid))));
          refit.status = table->view().status[entry];
          refit.failedPass = table->view().failedPass[entry];
        } else {
          refit = refitRowsWithHypothesis(event.propagator, geometry, material, rows, static_cast<int>(chain.rows.size()), species.pid);
        }
        if (refit.status != PropagationStatus::kOk) {
          one.refitFailed = true;
          one.refitPass = static_cast<std::uint8_t>(refit.failedPass);
          one.refitStatus = static_cast<std::uint8_t>(refit.status);
          break;
        }
        if (table) {
          daughters[i] = daughterOf(*tableSource, entry);
        } else {
          auto propagator = event.propagator;
          propagator.massGeV = species.mass;
          daughters[i] = {.innermost = refit.innermost, .surface = geometry.planeOf(rows[0]), .propagator = propagator, .mass = species.mass, .chargeMagnitude = std::abs(particleCharge(leaf->pdg))};
          if (species.pid == o2::track::PID::Electron) {
            daughters[i].mixture = refit.innermostMixture;
          }
        }
        masses[i] = species.mass;
      }
      // The unseen one: an unmeasured daughter, last, its start per job.
      double unseenMass = 0.;
      if (one.unseen != nullptr) {
        unseenMass = particleMass(one.unseen->pdg);
        daughters.push_back({.mass = unseenMass, .neutral = true, .unmeasured = true});
        masses.push_back(unseenMass);
      }

      // The truth.
      const auto& mother = *one.mother;
      const auto& decay = one.decay;
      Vec3D pTrue{};
      double eTrue = 0.;
      double minPt = std::numeric_limits<double>::infinity();
      std::vector<Vec3D> qTrue(count);
      for (std::size_t i = 0; i < count; ++i) {
        const auto& m = one.leaves[i]->momentum;
        const double p = std::hypot(m[0], m[1], m[2]);
        eTrue += energyOf(decay_fit_detail::Momentum<double>{m, p}, masses[i]);
        for (int k = 0; k < 3; ++k) {
          pTrue[k] += m[k];
        }
        minPt = std::min(minPt, std::hypot(m[0], m[1]));
        qTrue[i] = {phiOf(m), thetaOf(m), particleCharge(one.leaves[i]->pdg) / p};
      }
      if (one.unseen != nullptr) {
        const auto& m = one.unseen->momentum;
        eTrue += energyOf(decay_fit_detail::Momentum<double>{m, std::hypot(m[0], m[1], m[2])}, unseenMass);
        for (int k = 0; k < 3; ++k) {
          pTrue[k] += m[k];
        }
      }
      const double pTotal = std::hypot(pTrue[0], pTrue[1], pTrue[2]);
      // The flight's frame: along it, across it in the transverse plane, and
      // the third axis.
      Vec3D along{decay[0] - mother.vertex[0], decay[1] - mother.vertex[1], decay[2] - mother.vertex[2]};
      const double flight = std::hypot(along[0], along[1], along[2]);
      for (auto& component : along) {
        component /= flight;
      }
      Vec3D across{-along[1], along[0], 0.};
      const double acrossNorm = std::hypot(across[0], across[1]);
      if (acrossNorm > 0.) {
        across = {across[0] / acrossNorm, across[1] / acrossNorm, 0.};
      } else {
        across = {1., 0., 0.};
      }
      const Vec3D third = crossOf(along, across);

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
        .massTrue = static_cast<float>(sqrtClampedOf((eTrue * eTrue) - (pTotal * pTotal))),
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
        const auto* segment = event.hits->crossing(leaf->id, first);
        if (segment == nullptr) {
          continue;
        }
        const double pProduction = std::hypot(leaf->momentum[0], leaf->momentum[1], leaf->momentum[2]);
        const double pFirst = std::hypot(segment->momentum[0], segment->momentum[1], segment->momentum[2]);
        base.lossTrue[i] = static_cast<float>(1. - (pFirst / pProduction));
        const auto budget = event.material->getMatBudget(static_cast<float>(leaf->vertex[0]), static_cast<float>(leaf->vertex[1]), static_cast<float>(leaf->vertex[2]), static_cast<float>(segment->start[0]), static_cast<float>(segment->start[1]), static_cast<float>(segment->start[2]));
        const auto loss = energyLoss(static_cast<Scalar>(1. / pProduction), static_cast<Scalar>(masses[i]), -budget.meanRho * budget.length);
        if (loss.ok) {
          base.lossModel[i] = static_cast<float>(1. - (1. / static_cast<double>(loss.momentumScale)));
        }
        base.dPRelInner[i] = static_cast<float>((numeratorOf(daughters[i].charge()) / std::abs(static_cast<double>(daughters[i].innermost.qOverP())) / pFirst) - 1.);
        base.hitsBefore[i] = static_cast<std::int16_t>(std::ranges::count_if(event.hits->segments(leaf->id), [&](const HitSegment* other) { return other->time < segment->time; }));
      }

      // The mother: its PDG mass, its collision's vertex, the true decay length.
      const double motherMass = particleMass(mother.pdg);
      const auto collision = event.collisions->find(collisionOf(mother.id));
      const bool hasCollision = collision != event.collisions->end();
      const Vec3D point = hasCollision ? collision->second : Vec3D{};
      // Made at the collision: its production point is the collision's to
      // the precision the kinematics store it (float); otherwise (a Lambda
      // from a Xi) a production point at the collision is false.
      const bool madeAtCollision = hasCollision && std::ranges::equal(mother.vertex, point, [](double a, double b) { return static_cast<float>(a) == static_cast<float>(b); });
      if (hasCollision) {
        base.productionOffset = static_cast<float>(*productionOffset(*event.collisions, mother.id, mother.vertex));
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
      // Omega), a charged composite's (a tracked Xi in a Xi_c) --, one per
      // layer as a track holds them, the outermost first, with its mass
      // (DecayVertexFit.h's chains).
      // The fit's DecayMother holds spans of its rows: their storage here,
      // reserved so that the spans stay valid.
      struct OwnRows {
        std::vector<TrackingSurface> planes;
        std::vector<Vec2D> positions;
        std::vector<Vec3D> covariance;
      };
      std::vector<OwnRows> motherRows;
      std::vector<DecayMother> trackedMothers;
      if (tree && !one.refitFailed) {
        std::vector<std::pair<const TreeParticle*, int>> charged{{&mother, 0}};
        for (std::size_t w = 0; w < one.composites.size(); ++w) {
          charged.emplace_back(one.composites[w], static_cast<int>(w) + 1);
        }
        motherRows.reserve(charged.size());
        for (const auto& [particle, vertex] : charged) {
          const auto found = chainOf.find(particle->id);
          // a species the PDG table lacks has no charge either
          if (particleCharge(particle->pdg) == 0. || found == chainOf.end() || found->second->rows.empty()) {
            continue;
          }
          const auto own = onePerLayer(event.rows, found->second->rows);
          const auto trimmed = static_cast<std::uint8_t>(std::min<std::size_t>(found->second->rows.size() - own.size(), 255U));
          if (vertex == 0) {
            base.motherRows = static_cast<std::uint8_t>(std::min<std::size_t>(own.size(), 255U));
            base.motherTrimmed = trimmed;
          } else {
            base.compositesTrimmed = static_cast<std::uint8_t>(std::min(base.compositesTrimmed + trimmed, 255));
          }
          auto propagator = event.propagator;
          propagator.massGeV = static_cast<Scalar>(particleMass(particle->pdg));
          auto& rows = motherRows.emplace_back();
          for (auto at = own.rbegin(); at != own.rend(); ++at) {
            const auto m = frame.measurement(*at);
            rows.planes.push_back(geometry.planeOf(m));
            rows.positions.push_back({m.loc0, m.loc1});
            rows.covariance.push_back({m.varLoc0, m.covLoc01, m.varLoc1});
          }
          trackedMothers.push_back({.vertex = vertex, .planes = rows.planes, .positions = rows.positions, .covariance = rows.covariance, .propagator = propagator});
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
      const bool electrons = std::ranges::any_of(daughters, [](const DecayDaughter& d) { return d.mixture.size > 0; });
      if (electrons) {
        for (const std::uint8_t set : {std::uint8_t{0}, tree ? std::uint8_t{16} : std::uint8_t{1}}) {
          if (!missing || set != 1U) {
            jobs.push_back({DecayFitVariant::kGaussianSum, set});
          }
        }
      }
      // A tree's composites: their PDG masses, the first one's truth.
      std::vector<double> compositeMasses;
      compositeMasses.reserve(one.composites.size());
      for (const auto* composite : one.composites) {
        compositeMasses.push_back(particleMass(composite->pdg));
      }
      // each vertex's particle's mass, for its flight's material in the fit
      one.topology.massGeV.clear();
      one.topology.massGeV.push_back(motherMass);
      for (const double m : compositeMasses) {
        one.topology.massGeV.push_back(m);
      }
      std::vector<int> childLeaves; // the first composite's leaves
      Vec3D childDecay{};
      Vec3D childAxis{}; // its true flight direction
      Vec3D childAcross{};
      if (tree) {
        base.childPdg = one.composites.front()->pdg;
        childDecay = decayPointOf(trees, *one.composites.front());
        Vec3D childP{};
        double childE = 0.;
        for (std::size_t i = 0; i < count; ++i) {
          if (one.topology.under(static_cast<int>(i), 1)) {
            childLeaves.push_back(static_cast<int>(i));
            const auto& m = one.leaves[i]->momentum;
            const double p = std::hypot(m[0], m[1], m[2]);
            childE += energyOf(decay_fit_detail::Momentum<double>{m, p}, masses[i]);
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
        childAcross = t > 0. ? Vec3D{-childAxis[1] / t, childAxis[0] / t, 0.} : Vec3D{1., 0., 0.};
        base.childLengthTrue = static_cast<float>(((childDecay[0] - decay[0]) * childAxis[0]) + ((childDecay[1] - decay[1]) * childAxis[1]) + ((childDecay[2] - decay[2]) * childAxis[2]));
        base.childMassTrue = static_cast<float>(sqrtClampedOf((childE * childE) - (childPTotal * childPTotal)));
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
          record.status = static_cast<std::uint8_t>(VertexFitStatus::kNoStart); // nothing to constrain to
          one.records.push_back(record);
          continue;
        }
        DecayConstraints constraints{.mass = wantMass, .massGeV = motherMass, .production = wantProduction, .point = point};
        if (wantComposites) {
          for (const double m : compositeMasses) {
            constraints.compositeMasses.push_back(m);
          }
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
              constraints.resonances.push_back(taken);
            }
          }
        }
        const bool wantMother = (job.constraints & 32U) != 0U;
        const auto tracked = wantMother ? std::span<const DecayMother>{trackedMothers} : std::span<const DecayMother>{};
        // The unseen one with the visible ones only under the missing fit,
        // from the visible fit's starts (`missingMomentumStarts`).
        record.branch = static_cast<std::int8_t>(job.branch);
        auto jobDaughters = daughters;
        auto jobTopology = one.topology;
        const bool wantMissing = (job.constraints & 64U) != 0U;
        if (wantMissing) {
          std::optional<Vec3D> start;
          if (free.status == VertexFitStatus::kOk) {
            const auto starts = missingMomentumStarts(free, masses, point, motherMass, unseenMass);
            if (static_cast<std::size_t>(job.branch) < starts.size()) {
              start = starts[static_cast<std::size_t>(job.branch)];
            }
          }
          const double unseenP = start ? std::hypot((*start)[0], (*start)[1], (*start)[2]) : 0.;
          if (!(unseenP > 0.)) {
            record.status = static_cast<std::uint8_t>(VertexFitStatus::kNoStart);
            one.records.push_back(record);
            continue;
          }
          startUnmeasured(jobDaughters.back(), *start);
        } else if (missing) {
          // the visible daughters alone: the unseen one out of the tree too
          jobDaughters.pop_back();
          if (!jobTopology.vertexOf.empty()) {
            jobTopology.vertexOf.pop_back();
          }
        }
        const auto fitUnder = [&](const DecayConstraints& set) {
          const std::span<const DecayDaughter> span{jobDaughters};
          return withDecayFitVariant(job.variant, [&](auto variant) { return fitDecayVertex<Det, decltype(variant)::value>(span, event.geometry->material(), set, jobTopology, tracked); });
        };
        const auto started = std::chrono::steady_clock::now();
        const auto inPrecisions = [&]() {
          const std::span<const DecayDaughter> span{jobDaughters};
          const DecayPrecisionRecord identity{.timeframe = record.timeframe, .particle = record.particle, .pdg = record.pdg, .variant = record.variant, .constraints = record.constraints, .branch = record.branch, .resonance = record.resonance, .daughters = static_cast<std::uint8_t>(jobDaughters.size())};
          return withDecayFitVariant(job.variant, [&](auto variant) { return measurePrecision<Det, decltype(variant)::value>(span, event.geometry->material(), constraints, jobTopology, tracked, masses, across, identity, one.precision); });
        };
        const auto outcome = mPrecision ? inPrecisions() : fitUnder(constraints);
        record.microseconds = static_cast<float>(std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - started).count());
        if (job.constraints == 0U && job.variant == DecayFitVariant::kFull) {
          free = outcome;
        }
        if (mMode == Mode::kDevice) {
          const auto v = static_cast<std::size_t>(job.variant);
          std::array<Vec3D, decay_fit_detail::kMaxConstrainedVertices> start{};
          const bool startFound = decayFitStart<Det>(jobDaughters, constraints, jobTopology, tracked, start);
          one.batches[v].add(jobDaughters, constraints, jobTopology, tracked, {start.data(), startFound ? static_cast<std::size_t>(jobTopology.vertices()) : 0U});
          one.device[v].push_back({.timeframe = record.timeframe, .particle = record.particle, .pdg = record.pdg, .variant = record.variant, .constraints = record.constraints, .branch = record.branch, .resonance = record.resonance, .daughters = static_cast<std::uint8_t>(jobDaughters.size()), .statusDirect = static_cast<std::uint8_t>(outcome.status)});
          one.direct[v].push_back(summarise(outcome, jobDaughters));
          continue;
        }
        if (mPrecision) {
          continue;
        }
        const auto& fit = outcome;
        record.status = static_cast<std::uint8_t>(fit.status);
        record.iterations = fit.iterations;
        record.stopReason = fit.stopReason;
        record.outer = fit.outer;
        record.modelShift = fit.modelShift;
        record.pointSigma = static_cast<float>(fit.pointSigma);
        record.combinations = fit.combinations;
        record.combinationsTried = fit.combinationsTried;
        record.stopMetric = fit.stopMetric;
        record.stopBound = fit.stopBound;
        if (fit.status != VertexFitStatus::kOk) {
          one.records.push_back(record);
          continue;
        }
        const int vertexCount = static_cast<int>(fit.vertices.size());
        const auto at = [&](int a, int b) { return fit.covarianceAt(a, b); };
        const auto& fitMasses = masses;
        record.chi2 = static_cast<float>(fit.chi2);
        record.ndf = static_cast<std::int16_t>(fit.ndf);
        record.maxCarryX2X0 = static_cast<float>(*std::ranges::max_element(fit.carryX2X0));
        const Vec3D d{fit.vertex[0] - decay[0], fit.vertex[1] - decay[1], fit.vertex[2] - decay[2]};
        record.dX = static_cast<float>(d[0]);
        record.dY = static_cast<float>(d[1]);
        record.dZ = static_cast<float>(d[2]);
        record.sX = static_cast<float>(std::sqrt(at(0, 0)));
        record.sY = static_cast<float>(std::sqrt(at(1, 1)));
        record.sZ = static_cast<float>(std::sqrt(at(2, 2)));
        const auto project = [&](const Vec3D& axis, float& residual, float& sigma) {
          residual = static_cast<float>(dotOf(d, axis));
          sigma = static_cast<float>(sqrtClampedOf(quadraticForm(axis, axis, 3, at)));
        };
        project(along, record.dL, record.sL);
        project(across, record.dT, record.sT);
        project(third, record.dN, record.sN);
        const auto fitted = fitKinematics(fit, fitMasses);
        record.mass = static_cast<float>(fitted.mass);
        record.sMass = static_cast<float>(fitted.sMass);
        record.p = static_cast<float>(fitted.p);
        record.sP = static_cast<float>(fitted.sP);
        if (hasCollision) {
          const auto motherFit = fitFlight(fit, MeasuredPoint{.position = point}, motherMass);
          record.decayLength = static_cast<float>(motherFit.length);
          record.sDecayLength = static_cast<float>(motherFit.sLength);
          record.ct = static_cast<float>(motherFit.ct);
          record.sCt = static_cast<float>(motherFit.sCt);
          record.dPhiMother = static_cast<float>(toPMPiOf(motherFit.phi - phiOf(pTrue)));
          record.dThetaMother = static_cast<float>(motherFit.theta - thetaOf(pTrue));
          record.sPhiMother = static_cast<float>(motherFit.sPhi);
          record.sThetaMother = static_cast<float>(motherFit.sTheta);
        }
        // The child: its vertex, decay length and mass.
        if (tree) {
          const auto& child = fit.vertices[1];
          const Vec3D e{child[0] - childDecay[0], child[1] - childDecay[1], child[2] - childDecay[2]};
          // the covariance of projections of the vertices from \a first and \a second (0: the parent's, 3: the child's)
          const auto vertexCovariance = [&](const Vec3D& axis, int first, int second) { return quadraticForm(axis, axis, 3, [&](int a, int b) { return at(first + a, second + b); }); };
          const auto childVariance = [&](const Vec3D& axis) { return vertexCovariance(axis, 3, 3); };
          record.dChildL = static_cast<float>(dotOf(e, childAxis));
          record.sChildL = static_cast<float>(sqrtClampedOf(childVariance(childAxis)));
          record.dChildT = static_cast<float>(dotOf(e, childAcross));
          record.sChildT = static_cast<float>(sqrtClampedOf(childVariance(childAcross)));
          // The decay length from the parent's vertex along the true flight;
          // its variance the two vertices' with their correlation.
          record.childLength = static_cast<float>(((child[0] - fit.vertex[0]) * childAxis[0]) + ((child[1] - fit.vertex[1]) * childAxis[1]) + ((child[2] - fit.vertex[2]) * childAxis[2]));
          const double lengthVariance = childVariance(childAxis) + vertexCovariance(childAxis, 0, 0) - (2. * vertexCovariance(childAxis, 0, 3));
          record.sChildLength = static_cast<float>(sqrtClampedOf(lengthVariance));
          const auto childKinematics = subsetKinematics(std::span<const Vec3D>{fit.momenta.data(), fit.momenta.size()}, masses, std::span<const double>{fit.charges.data(), fit.charges.size()}, childLeaves, vertexCount, at);
          record.childMass = static_cast<float>(childKinematics.mass);
          record.sChildMass = static_cast<float>(childKinematics.sMass);
        }
        // The carried momenta, of the seen daughters (an unseen one has no
        // carried state): independent of one another, each its own 3 x 3.
        const auto carried = [&](int a, int b) { return a / 3 == b / 3 ? fit.carriedCovariance[a / 3][covIndex(a % 3, b % 3)] : 0.; };
        const auto measured = decayKinematics(std::span(fit.carried).first(count), std::span(masses).first(count), std::span(fit.charges).first(count), carried);
        record.massCarried = static_cast<float>(measured.mass);
        record.sMassCarried = static_cast<float>(measured.sMass);
        for (std::size_t i = 0; i < std::min(count, DecayFitRecord::kDaughterColumns); ++i) {
          const int q = (3 * vertexCount) + (3 * static_cast<int>(i));
          record.dPhi[i] = static_cast<float>(toPMPiOf(fit.momenta[i][0] - qTrue[i][0]));
          record.dTheta[i] = static_cast<float>(fit.momenta[i][1] - qTrue[i][1]);
          record.dQOverP[i] = static_cast<float>(fit.momenta[i][2] - qTrue[i][2]);
          record.sPhi[i] = static_cast<float>(std::sqrt(at(q, q)));
          record.sTheta[i] = static_cast<float>(std::sqrt(at(q + 1, q + 1)));
          record.sQOverP[i] = static_cast<float>(std::sqrt(at(q + 2, q + 2)));
          record.dPRel[i] = static_cast<float>((std::abs(qTrue[i][2]) / std::abs(fit.momenta[i][2])) - 1.);
        }
        one.records.push_back(record);
      }
    }
  });
  if (mMode == Mode::kDevice) {
    std::array<DecayFitBatch, 2> batches;
    std::array<std::vector<DecayDeviceRecord>, 2> records;
    std::array<std::vector<FitSummary>, 2> direct;
    for (auto& one : work) {
      for (std::size_t v = 0; v < 2U; ++v) {
        batches[v].append(one.batches[v]);
        records[v].insert(records[v].end(), one.device[v].begin(), one.device[v].end());
        direct[v].insert(direct[v].end(), one.direct[v].begin(), one.direct[v].end());
      }
    }
    for (auto& variant : records) {
      for (auto& record : variant) {
        record.source = static_cast<std::uint8_t>(source);
      }
    }
    fitBatches(event, batches, records, direct);
  }
  for (auto& one : work) {
    if (mPrecision) {
      for (auto& record : one.precision) {
        record.source = static_cast<std::uint8_t>(source);
        mWriter->write(record);
      }
      mRecords += one.precision.size();
    } else if (mMode == Mode::kFit) {
      for (auto& record : one.records) {
        record.source = static_cast<std::uint8_t>(source);
        mWriter->write(record);
      }
      mRecords += one.records.size();
    }
    mRefitFailed += one.refitFailed ? 1U : 0U;
  }
  mDecays += work.size();
}

template <DetectorLayout Det>
void DecayFitStudy<Det>::fitBatches(const TruthEvent<Det>& event, std::array<DecayFitBatch, 2>& batches, std::array<std::vector<DecayDeviceRecord>, 2>& records, std::array<std::vector<decay_fit_study_detail::FitSummary>, 2>& direct)
{
  using namespace decay_fit_study_detail;
  const auto& material = event.geometry->material();
  ++mTimeframes;
  mEvents += event.collisions != nullptr ? event.collisions->size() : 0U;
  if (mDevice && mConfigured != &material) {
    if (!mDevice->configure(material)) {
      throw std::runtime_error("decay device: the material does not fit the device arena");
    }
    mConfigured = &material;
  }
  for (std::size_t v = 0; v < 2U; ++v) {
    auto& batch = batches[v];
    if (batch.empty()) {
      continue;
    }
    const auto variant = static_cast<DecayFitVariant>(v);
    std::vector<DecayVertexFitResult> host(batch.size());
    const auto started = std::chrono::steady_clock::now();
    withDecayFitVariant(variant, [&](auto one) { fitDecayBatch<Det, decltype(one)::value>(batch.view(), material, host); });
    mHostSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    std::vector<DecayVertexFitResult> device;
    if (mDevice) {
      device.resize(batch.size());
      typename Device::Timing timing;
      if (!mDevice->fit(variant, batch, device, timing)) {
        throw std::runtime_error("decay device: a timeframe's decays do not fit the device arena");
      }
      mDeviceTotal.uploadSeconds += timing.uploadSeconds;
      mDeviceTotal.sizeSeconds += timing.sizeSeconds;
      mDeviceTotal.fitSeconds += timing.fitSeconds;
      mDeviceTotal.downloadSeconds += timing.downloadSeconds;
      mDeviceTotal.launches += timing.launches;
      mDeviceTotal.refused += timing.refused;
      mDeviceTotal.peakBytes = std::max(mDeviceTotal.peakBytes, timing.peakBytes);
    }
    const auto view = batch.view();
    for (std::size_t i = 0; i < batch.size(); ++i) {
      auto& record = records[v][i];
      const auto& input = view.decays[i];
      const std::span<const DecayDaughter> daughters{view.daughters + input.firstDaughter, input.daughters};
      const auto hostSummary = summarise(host[i], daughters);
      record.statusHost = hostSummary.status;
      record.hostSame = sameBits(hostSummary, direct[v][i]) ? 1U : 0U;
      mHostSame += record.hostSame;
      if (mDevice) {
        const auto deviceSummary = summarise(device[i], daughters);
        const auto& reference = direct[v][i];
        record.statusDevice = deviceSummary.status;
        mDeviceStatusSame += deviceSummary.status == reference.status ? 1U : 0U;
        if (deviceSummary.status == 0U && reference.status == 0U) {
          const auto& a = deviceSummary.vertices[0];
          const auto& b = reference.vertices[0];
          record.dVertex = static_cast<float>(1e4 * std::sqrt(((a[0] - b[0]) * (a[0] - b[0])) + ((a[1] - b[1]) * (a[1] - b[1])) + ((a[2] - b[2]) * (a[2] - b[2]))));
          record.dMass = static_cast<float>(deviceSummary.mass - reference.mass);
          record.chi2Rel = static_cast<float>((deviceSummary.chi2 - reference.chi2) / reference.chi2);
          double largest = 0.;
          for (std::size_t j = 0; j < reference.diagonal.size() && j < deviceSummary.diagonal.size(); ++j) {
            largest = std::max(largest, std::abs(deviceSummary.diagonal[j] - reference.diagonal[j]) / std::abs(reference.diagonal[j]));
          }
          record.diagonalRel = static_cast<float>(largest);
        }
      }
      mWriter->write(record);
    }
    mBatchFits += batch.size();
    mRecords += batch.size();
  }
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYFITSTUDY_H_
