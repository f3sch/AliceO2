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

/// \file OverlapFitStudy.h
/// \brief What a chain's second row on a layer (overlapping sensors) is
///        worth in the fit: the chain fitted without and with it.
///
/// Of each chain of an electron, muon, pion, kaon or proton that has two rows
/// on one layer on two sensors, two fits under the true PID
/// (`refitWithHypothesis`): the first row on each layer (variant 0), and the
/// first row on each sensor (variant 1), the overlap partners added. A row on
/// a sensor that already has one of the chain's rows is left out of both:
/// two clusters of one crossing, not a second measurement. Both fits start
/// on the same innermost row and take the link material as the final refit
/// does, none on the step between a layer's two sensors (`LayerMaterial`):
/// the bands hold a layer's whole material about either sensor of an
/// overlap. Chains on fewer layers than a pass starts with are refused. Each is held to the MC at that row
/// (`stateTruth`) and, carried to the perigee about the production point,
/// to the MC at production (`perigeeTruth`). A record per fit; `status` as
/// SpeciesFitStudy's.

#ifndef GTRK_STUDY_OVERLAPFITSTUDY_H_
#define GTRK_STUDY_OVERLAPFITSTUDY_H_

#include "GenericTracking/Tracking/HypothesisRefit.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Study/MCTruth.h"
#include "GenericTracking/Study/RecordWriter.h"
#include "GenericTracking/Study/StateTruth.h"
#include "GenericTracking/Study/StudyChains.h"
#include "GenericTracking/Study/TruthEvent.h"
#include "GenericTracking/Study/TruthFit.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <format>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace o2::gtrk::study
{

/// One fit of a chain with an overlap pair (`OverlapFitStudy`).
struct OverlapFitRecord {
  static constexpr std::string_view kName = "overlapFit";
  static constexpr std::uint8_t kNoCrossing = 255U; ///< `status`: a fit without the MC crossing
  std::uint32_t timeframe = 0U;
  std::uint8_t source = 0U; ///< ChainSource
  std::uint64_t particle = 0U;
  std::int32_t pdg = 0;
  std::uint8_t origin = 0U;
  float pt = 0.F;           ///< the truth's at production [GeV]
  float eta = 0.F;          ///< the truth's at production
  std::uint8_t layers = 0U; ///< the distinct layers of the chain's rows
  std::uint8_t pairs = 0U;  ///< the overlap partners added in variant 1
  /// The layers with an overlap pair, a bit each.
  std::uint64_t pairLayers = 0U;
  std::uint8_t variant = 0U; ///< 0 a row per layer, 1 with the overlap partners
  std::uint8_t rows = 0U;    ///< the fit's rows
  std::uint8_t status = 0U;
  float chi2 = 0.F;
  std::int16_t ndf = 0;
  float p = 0.F;     ///< |p| at the innermost row's crossing [GeV]
  float dLoc0 = 0.F; ///< fitted less true at the innermost row
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
  /// The carry to the perigee's PropagationStatus (255 without a fit), and
  /// with kOk the state there less the MC at production, and its sigmas.
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
    visit("particle", particle);
    visit("pdg", pdg);
    visit("origin", origin);
    visit("pt", pt);
    visit("eta", eta);
    visit("layers", layers);
    visit("pairs", pairs);
    visit("pairLayers", pairLayers);
    visit("variant", variant);
    visit("rows", rows);
    visit("status", status);
    visit("chi2", chi2);
    visit("ndf", ndf);
    visit("p", p);
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
class OverlapFitStudy
{
 public:
  explicit OverlapFitStudy(RecordWriter& writer) : mWriter(&writer) { writer.begin<OverlapFitRecord>(); }

  /// \a chains of \a event, from \a source.
  void process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source);

  [[nodiscard]] std::vector<std::string> summary() const
  {
    return {std::format("overlap fit study: {} chains with an overlap pair fitted, {} records into {}; {} chains without one, {} on fewer than {} layers, {} of another species or charge", mChains, mRecords, mWriter->path(), mNoPair, mFewLayers, o2::gtrk::detail::kPassStart, mUnsupported)};
  }

 private:
  struct Work {
    const Chain* chain = nullptr;
    Species species;
    double charge = 0.;
    std::vector<std::uint32_t> perLayer;  ///< variant 0's rows
    std::vector<std::uint32_t> perSensor; ///< variant 1's rows
    std::uint8_t layers = 0U;
    std::uint64_t pairLayers = 0U;
    std::vector<OverlapFitRecord> records;
  };

  void fitChain(const TruthEvent<Det>& event, ChainSource source, Work& work) const;

  RecordWriter* mWriter = nullptr;
  std::size_t mChains = 0U;
  std::size_t mNoPair = 0U;
  std::size_t mFewLayers = 0U;
  std::size_t mUnsupported = 0U;
  std::size_t mRecords = 0U;
};

template <DetectorLayout Det>
void OverlapFitStudy<Det>::process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source)
{
  std::vector<Work> work;
  for (const auto& chain : chains) {
    const auto& truth = event.particles->at(chain.particle);
    const auto species = truthSpecies(truth.pdg);
    const double charge = particleCharge(truth.pdg);
    if (!species || std::lround(std::abs(charge)) != 1) {
      ++mUnsupported;
      continue;
    }
    Work one{.chain = &chain, .species = *species, .charge = charge};
    std::uint64_t layers = 0U;
    std::vector<std::uint32_t> surfaces;
    for (const auto at : chain.rows) {
      const auto bit = std::uint64_t{1} << event.rows[at].layer;
      const auto surface = event.frame->measurement(at).surface;
      if (std::ranges::find(surfaces, surface) != surfaces.end()) {
        continue;
      }
      surfaces.push_back(surface);
      one.perSensor.push_back(at);
      if ((layers & bit) == 0U) {
        one.perLayer.push_back(at);
        layers |= bit;
      } else {
        one.pairLayers |= bit;
      }
    }
    if (one.pairLayers == 0U) {
      ++mNoPair;
      continue;
    }
    one.layers = static_cast<std::uint8_t>(std::popcount(layers));
    if (one.layers < o2::gtrk::detail::kPassStart) {
      ++mFewLayers;
      continue;
    }
    work.push_back(std::move(one));
  }
  o2::gtrk::detail::parallelFor(*event.arena, 0U, work.size(), 16, [&](std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      fitChain(event, source, work[i]);
    }
  });
  for (const auto& one : work) {
    for (const auto& record : one.records) {
      mWriter->write(record);
    }
    mRecords += one.records.size();
  }
  mChains += work.size();
}

template <DetectorLayout Det>
void OverlapFitStudy<Det>::fitChain(const TruthEvent<Det>& event, ChainSource source, Work& work) const
{
  const auto& chain = *work.chain;
  const auto& truth = event.particles->at(chain.particle);
  const auto& frame = *event.frame;
  const auto& geometry = *event.geometry;
  const o2::gtrk::detail::LayerMaterial<o2::base::MatLayerCylSet> material{.inner = event.material};
  const OverlapFitRecord base{
    .timeframe = event.timeframe,
    .source = static_cast<std::uint8_t>(source),
    .particle = chain.particle,
    .pdg = truth.pdg,
    .origin = static_cast<std::uint8_t>(truth.origin),
    .pt = truth.pt,
    .eta = truth.eta,
    .layers = work.layers,
    .pairs = static_cast<std::uint8_t>(work.perSensor.size() - work.perLayer.size()),
    .pairLayers = work.pairLayers,
  };
  for (std::uint8_t variant = 0U; variant < 2U; ++variant) {
    const auto& rows = variant == 0U ? work.perLayer : work.perSensor;
    const auto count = static_cast<int>(rows.size());
    OverlapFitRecord record = base;
    record.variant = variant;
    record.rows = static_cast<std::uint8_t>(count);
    const RowChain<Det> measurements{.frame = &frame, .rows = rows.data()};
    const auto refit = refitWithHypothesis(event.propagator, geometry, material, measurements, count, work.species.pid, truth.vx, truth.vy);
    const bool fitted = refit.status == PropagationStatus::kOk || refit.failedPass == FitPass::kNone;
    record.chi2 = refit.chi2;
    record.ndf = refit.ndf;
    record.status = static_cast<std::uint8_t>(fitted ? PropagationStatus::kOk : refit.status);
    if (fitted) {
      record.perigeeStatus = static_cast<std::uint8_t>(refit.status);
      const auto at = event.hits == nullptr ? std::nullopt : stateTruth(geometry, *event.hits, frame, refit.innermost, chain.particle, rows.front(), work.charge);
      if (!at) {
        record.status = OverlapFitRecord::kNoCrossing;
      } else {
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
      }
    }
    if (refit.status == PropagationStatus::kOk) {
      const auto perigee = perigeeTruth(refit.perigee, truth, work.charge);
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
    work.records.push_back(record);
  }
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_OVERLAPFITSTUDY_H_
