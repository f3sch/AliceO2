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

/// \file SpeciesFitStudy.h
/// \brief A chain's rows (StudyChains.h) fitted with its particle's true
///        species (TruthFit.h) and held to the MC.
///
/// For every chain of an electron, muon, pion, kaon or proton, every
/// `TruthFitVariant`: a record each. Each fit's state on the chain's
/// comparison row -- a track's row on its fitted state's surface, a truth
/// chain's innermost row -- is held to the particle's MC crossing there
/// (`stateTruth`); the true PID's refit (`kHypothesis`) also at the perigee
/// about the MC production point (`perigeeTruth`). The inward-only variants
/// start from the chain's outward state: a track's own, a truth chain's pion
/// fit's (when its outward pass succeeded).
///
/// A record's `status` is the fit's PropagationStatus at the rows, 254 for
/// an inward variant without an outward state to start from, 255 for a fit
/// without the MC crossing to hold it to; its residuals only with kOk.

#ifndef GTRK_STUDY_SPECIESFITSTUDY_H_
#define GTRK_STUDY_SPECIESFITSTUDY_H_

#include "GenericTracking/Tracking/HypothesisRefit.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Study/MCTruth.h"
#include "GenericTracking/Study/RecordWriter.h"
#include "GenericTracking/Study/StateTruth.h"
#include "GenericTracking/Study/StudyChains.h"
#include "GenericTracking/Study/TruthEvent.h"
#include "GenericTracking/Study/TruthFit.h"

#include <bit>
#include <cmath>
#include <cstdint>
#include <format>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace o2::gtrk::study
{

/// One fit of a chain (`SpeciesFitStudy`).
struct SpeciesFitRecord {
  static constexpr std::string_view kName = "speciesFit";
  static constexpr std::uint8_t kNoSeed = 254U;     ///< `status`: an inward variant with no outward state
  static constexpr std::uint8_t kNoCrossing = 255U; ///< `status`: a fit without the MC crossing
  std::uint32_t timeframe = 0U;
  std::uint8_t source = 0U; ///< ChainSource
  std::int32_t track = -1;  ///< the tracker's track, -1 for a truth chain
  std::uint64_t particle = 0U;
  std::uint8_t leg = 0U;
  std::int32_t pdg = 0;
  std::int8_t charge = 0; ///< [e]
  std::uint8_t origin = 0U;
  float eta = 0.F;            ///< the truth's at production
  std::uint8_t layer = 0U;    ///< the comparison row's
  std::uint8_t kinks = 0U;    ///< a track's own refit's; 0 for a truth chain
  std::uint8_t clusters = 0U; ///< the chain's rows
  std::uint8_t layers = 0U;   ///< the distinct layers they are on
  /// The production point's distance from its collision's vertex [cm]: 0
  /// for a particle made at the collision.
  float productionOffset = std::numeric_limits<float>::quiet_NaN();
  float p = 0.F; ///< |p| at the comparison row's crossing [GeV], 0 without one
  /// |p| at the crossing of the chain's outermost row [GeV]: the momentum lost
  /// inside the chain, 1 - pOuter / p.
  float pOuter = std::numeric_limits<float>::quiet_NaN();
  std::uint8_t variant = 0U; ///< TruthFitVariant
  std::uint8_t status = 0U;
  float chi2 = 0.F;
  std::int16_t ndf = 0;
  float dLoc0 = 0.F; ///< fitted less true at the comparison row
  float dLoc1 = 0.F;
  float dPhi = 0.F;
  float dTheta = 0.F;
  float dQOverP = 0.F;
  float sLoc0 = 0.F; ///< the fit's sigmas there
  float sLoc1 = 0.F;
  float sPhi = 0.F;
  float sTheta = 0.F;
  float sQOverP = 0.F;
  float dPtRel = 0.F; ///< fitted pT over the crossing's, less 1
  /// The inward variants' outward state: 0 the track's own, 1 the chain's
  /// pion fit's, 2 none.
  std::uint8_t seed = 0U;
  /// `kHypothesis` only: the carry to the perigee's PropagationStatus (255
  /// for the other variants), and with kOk the state there less the MC at
  /// production, and its sigmas.
  std::uint8_t perigeeStatus = 255U;
  float perigeeD0 = std::numeric_limits<float>::quiet_NaN();
  float perigeeZ = std::numeric_limits<float>::quiet_NaN();
  float perigeePhi = std::numeric_limits<float>::quiet_NaN();
  float perigeeTheta = std::numeric_limits<float>::quiet_NaN();
  float perigeeQOverP = std::numeric_limits<float>::quiet_NaN();
  float perigeeSD0 = std::numeric_limits<float>::quiet_NaN();
  float perigeeSZ = std::numeric_limits<float>::quiet_NaN();
  float perigeeSPhi = std::numeric_limits<float>::quiet_NaN();
  float perigeeSTheta = std::numeric_limits<float>::quiet_NaN();
  float perigeeSQOverP = std::numeric_limits<float>::quiet_NaN();

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("source", source);
    visit("track", track);
    visit("particle", particle);
    visit("leg", leg);
    visit("pdg", pdg);
    visit("charge", charge);
    visit("origin", origin);
    visit("eta", eta);
    visit("layer", layer);
    visit("kinks", kinks);
    visit("clusters", clusters);
    visit("layers", layers);
    visit("productionOffset", productionOffset);
    visit("p", p);
    visit("pOuter", pOuter);
    visit("variant", variant);
    visit("status", status);
    visit("chi2", chi2);
    visit("ndf", ndf);
    visit("dLoc0", dLoc0);
    visit("dLoc1", dLoc1);
    visit("dPhi", dPhi);
    visit("dTheta", dTheta);
    visit("dQOverP", dQOverP);
    visit("sLoc0", sLoc0);
    visit("sLoc1", sLoc1);
    visit("sPhi", sPhi);
    visit("sTheta", sTheta);
    visit("sQOverP", sQOverP);
    visit("dPtRel", dPtRel);
    visit("seed", seed);
    visit("perigeeStatus", perigeeStatus);
    visit("perigeeD0", perigeeD0);
    visit("perigeeZ", perigeeZ);
    visit("perigeePhi", perigeePhi);
    visit("perigeeTheta", perigeeTheta);
    visit("perigeeQOverP", perigeeQOverP);
    visit("perigeeSD0", perigeeSD0);
    visit("perigeeSZ", perigeeSZ);
    visit("perigeeSPhi", perigeeSPhi);
    visit("perigeeSTheta", perigeeSTheta);
    visit("perigeeSQOverP", perigeeSQOverP);
  }
};

template <DetectorLayout Det>
class SpeciesFitStudy
{
 public:
  explicit SpeciesFitStudy(RecordWriter& writer) : mWriter(&writer) { writer.begin<SpeciesFitRecord>(); }

  /// \a chains of \a event, from \a source.
  void process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source);

  [[nodiscard]] std::vector<std::string> summary() const
  {
    return {std::format("species fit study: {} chains, {} records into {}; refused {} of another species or charge; {} without an outward state for the inward variants", mChains, mRecords, mWriter->path(), mUnsupported, mNoSeed)};
  }

 private:
  /// One chain's records, in the variants' order.
  void fitChain(const TruthEvent<Det>& event, const Chain& chain, const Species& species, double charge, const ChainIdentity& identity, std::vector<SpeciesFitRecord>& records, bool& noSeed) const;

  RecordWriter* mWriter = nullptr;
  std::size_t mChains = 0U;
  std::size_t mUnsupported = 0U;
  std::size_t mNoSeed = 0U;
  std::size_t mRecords = 0U;
};

template <DetectorLayout Det>
void SpeciesFitStudy<Det>::process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source)
{
  struct Work {
    const Chain* chain = nullptr;
    Species species;
    double charge = 0.;
    std::vector<SpeciesFitRecord> records;
    bool noSeed = false;
  };
  std::vector<Work> work;
  for (const auto& chain : chains) {
    const auto& truth = event.particles->at(chain.particle);
    const auto species = truthSpecies(truth.pdg);
    const double charge = particleCharge(truth.pdg);
    if (!species || std::lround(std::abs(charge)) != 1) {
      ++mUnsupported;
      continue;
    }
    work.push_back({.chain = &chain, .species = *species, .charge = charge});
  }
  o2::gtrk::detail::parallelFor(*event.arena, 0U, work.size(), 16, [&](std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      const auto& chain = *work[i].chain;
      fitChain(event, chain, work[i].species, work[i].charge, identityOf(chain, event.timeframe, source), work[i].records, work[i].noSeed);
    }
  });
  for (const auto& one : work) {
    for (const auto& record : one.records) {
      mWriter->write(record);
    }
    mRecords += one.records.size();
    mNoSeed += one.noSeed ? 1U : 0U;
  }
  mChains += work.size();
}

template <DetectorLayout Det>
void SpeciesFitStudy<Det>::fitChain(const TruthEvent<Det>& event, const Chain& chain, const Species& species, double charge, const ChainIdentity& identity, std::vector<SpeciesFitRecord>& records, bool& noSeed) const
{
  const auto& truth = event.particles->at(chain.particle);
  const auto& frame = *event.frame;
  const auto& geometry = *event.geometry;
  // A truth chain may hold two rows on one layer (overlapping sensors): no
  // material on the step between them, which the layer's band holds already.
  constexpr o2::gtrk::detail::LayerMaterial material{};
  const auto count = static_cast<int>(chain.rows.size());
  const auto row = chain.track ? chain.track->comparisonRow : chain.rows.front();
  const RowChain<Det> rows{.frame = &frame, .rows = chain.rows.data()};
  const auto crossing = [&](std::uint32_t at) {
    const auto m = frame.measurement(at);
    return event.hits == nullptr ? nullptr : event.hits->crossing(chain.particle, m);
  };
  SpeciesFitRecord base{
    .timeframe = identity.timeframe,
    .source = identity.source,
    .track = identity.track,
    .particle = identity.particle,
    .leg = identity.leg,
    .pdg = truth.pdg,
    .charge = static_cast<std::int8_t>(std::lround(charge)),
    .origin = static_cast<std::uint8_t>(truth.origin),
    .eta = truth.eta,
    .layer = static_cast<std::uint8_t>(event.rows[row].layer),
    .kinks = chain.track ? chain.track->kinks : std::uint8_t{0},
    .clusters = static_cast<std::uint8_t>(count),
  };
  std::uint64_t layerBits = 0U;
  for (const auto at : chain.rows) {
    layerBits |= std::uint64_t{1} << event.rows[at].layer;
  }
  base.layers = static_cast<std::uint8_t>(std::popcount(layerBits));
  if (event.collisions != nullptr) {
    if (const auto offset = productionOffset(*event.collisions, chain.particle, {static_cast<double>(truth.vx), static_cast<double>(truth.vy), static_cast<double>(truth.vz)})) {
      base.productionOffset = static_cast<float>(*offset);
    }
  }
  // The comparison row's crossing |p| for every record, failed fits too.
  if (const auto* at = crossing(row); at != nullptr) {
    base.p = static_cast<float>(std::hypot(at->momentum[0], at->momentum[1], at->momentum[2]));
  }
  if (const auto* outer = crossing(chain.rows.back()); outer != nullptr) {
    base.pOuter = static_cast<float>(std::hypot(outer->momentum[0], outer->momentum[1], outer->momentum[2]));
  }
  // Each fit's state at the comparison row against the MC crossing there.
  const auto held = [&](SpeciesFitRecord& record, PropagationStatus status, const BoundTrackState& state) {
    record.status = static_cast<std::uint8_t>(status);
    if (status != PropagationStatus::kOk) {
      return;
    }
    const auto at = event.hits == nullptr ? std::nullopt : stateTruth(geometry, *event.hits, frame, state, chain.particle, row, charge);
    if (!at) {
      record.status = SpeciesFitRecord::kNoCrossing;
      return;
    }
    record.p = static_cast<float>(at->momentum);
    record.dPtRel = static_cast<float>((at->fittedPt / at->pt) - 1.);
    record.dLoc0 = static_cast<float>(at->residual[0]);
    record.dLoc1 = static_cast<float>(at->residual[1]);
    record.dPhi = static_cast<float>(at->residual[2]);
    record.dTheta = static_cast<float>(at->residual[3]);
    record.dQOverP = static_cast<float>(at->residual[4]);
    record.sLoc0 = static_cast<float>(at->sigma[0]);
    record.sLoc1 = static_cast<float>(at->sigma[1]);
    record.sPhi = static_cast<float>(at->sigma[2]);
    record.sTheta = static_cast<float>(at->sigma[3]);
    record.sQOverP = static_cast<float>(at->sigma[4]);
  };
  // The inward variants' outward state: a track's own, else the pion fit's.
  std::optional<BoundTrackState> outward;
  std::uint8_t seed = 2U;
  if (chain.track) {
    outward = chain.track->outward;
    seed = 0U;
  }
  std::array<TrackFitResult, static_cast<std::size_t>(TruthFitVariant::kHypothesis)> fits{};
  for (std::size_t v = 0; v < fits.size(); ++v) {
    const auto variant = static_cast<TruthFitVariant>(v);
    const bool inwardOnly = variant == TruthFitVariant::kTrueInward || variant == TruthFitVariant::kElectronInward;
    if (inwardOnly && !chain.track) {
      const auto& pion = fits[static_cast<std::size_t>(TruthFitVariant::kPion)];
      if (pion.failedPass != FitPass::kOutward) {
        outward = pion.outwardState;
        seed = 1U;
      }
    }
    SpeciesFitRecord record = base;
    record.variant = static_cast<std::uint8_t>(v);
    if (inwardOnly && !outward) {
      record.status = SpeciesFitRecord::kNoSeed;
      noSeed = true;
      records.push_back(record);
      continue;
    }
    fits[v] = truthFit(variant, species.mass, event.propagator, geometry, material, rows, count, outward ? *outward : BoundTrackState{});
    record.chi2 = fits[v].chi2;
    record.ndf = fits[v].ndf;
    held(record, fits[v].status, fits[v].state);
    records.push_back(record);
  }
  // The true PID's refit, at the rows and at the perigee about the production
  // point: a failed carry to the perigee leaves the fit at the rows standing.
  const auto refit = refitWithHypothesis(event.propagator, geometry, material, rows, count, species.pid, truth.vx, truth.vy);
  const bool fitted = refit.status == PropagationStatus::kOk || refit.failedPass == FitPass::kNone;
  SpeciesFitRecord record = base;
  record.variant = static_cast<std::uint8_t>(TruthFitVariant::kHypothesis);
  record.chi2 = refit.chi2;
  record.ndf = refit.ndf;
  held(record, fitted ? PropagationStatus::kOk : refit.status, refit.innermost);
  if (fitted) {
    record.perigeeStatus = static_cast<std::uint8_t>(refit.status);
  }
  if (refit.status == PropagationStatus::kOk) {
    const auto perigee = perigeeTruth(refit.perigee, truth, charge);
    record.perigeeD0 = static_cast<float>(perigee.residual[0]);
    record.perigeeZ = static_cast<float>(perigee.residual[1]);
    record.perigeePhi = static_cast<float>(perigee.residual[2]);
    record.perigeeTheta = static_cast<float>(perigee.residual[3]);
    record.perigeeQOverP = static_cast<float>(perigee.residual[4]);
    record.perigeeSD0 = static_cast<float>(perigee.sigma[0]);
    record.perigeeSZ = static_cast<float>(perigee.sigma[1]);
    record.perigeeSPhi = static_cast<float>(perigee.sigma[2]);
    record.perigeeSTheta = static_cast<float>(perigee.sigma[3]);
    record.perigeeSQOverP = static_cast<float>(perigee.sigma[4]);
  }
  records.push_back(record);
  for (auto& one : records) {
    one.seed = seed;
  }
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_SPECIESFITSTUDY_H_
