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

/// \file SpeciesLikelihoodStudy.h
/// \brief A chain's rows refitted under every species of the tracking's list
///        (Species.h), each fit's likelihood of the rows.
///
/// For every chain (StudyChains.h; the tracker's good tracks or the
/// particles' own rows) and every species: a record of the Gaussian-sum fit
/// of its rows at the species' mass and |q| (`fitGaussianSum`) -- with Bethe
/// and Heitler's radiation for the electron, with none (`NoRadiationLoss`,
/// the Kalman fit) for every other species. Each record holds the fit's
/// status, chi2 and ndf and both passes' -2 ln L, the likelihood of the rows
/// under the fit's predictions (the full Gaussian: chi2 plus ln det of each
/// prediction's covariance), the start's rows by the diffuse likelihood. On
/// the same rows the species compare by their -2 ln L alone: its difference
/// between two species is -2 ln of their likelihood ratio.

#ifndef GTRK_STUDY_SPECIESLIKELIHOODSTUDY_H_
#define GTRK_STUDY_SPECIESLIKELIHOODSTUDY_H_

#include "GenericTracking/Base/Species.h"
#include "GenericTracking/Tracking/GSFTrackFitter.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Study/MCTruth.h"
#include "GenericTracking/Study/RecordWriter.h"
#include "GenericTracking/Study/StateTruth.h"
#include "GenericTracking/Study/StudyChains.h"
#include "GenericTracking/Study/TruthEvent.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <format>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace o2::gtrk::study
{

/// One species' fit of a chain (`SpeciesLikelihoodStudy`).
struct SpeciesLikelihoodRecord {
  static constexpr std::string_view kName = "speciesLikelihood";
  std::uint32_t timeframe = 0U;
  std::uint8_t source = 0U; ///< ChainSource
  std::int32_t track = -1;  ///< the tracker's track, -1 for a truth chain
  std::uint64_t particle = 0U;
  std::uint8_t leg = 0U;
  std::int32_t pdg = 0;
  std::int8_t charge = 0; ///< [e]
  std::uint8_t origin = 0U;
  float eta = 0.F;            ///< the truth's at production
  std::uint8_t clusters = 0U; ///< the chain's rows
  std::uint8_t layers = 0U;   ///< the distinct layers they are on
  /// The production point's distance from its collision's vertex [cm]: 0
  /// for a particle made at the collision.
  float productionOffset = std::numeric_limits<float>::quiet_NaN();
  /// |p| at the MC crossing of the chain's innermost row [GeV], NaN without one.
  float p = std::numeric_limits<float>::quiet_NaN();
  std::uint8_t species = 0U; ///< the hypothesis (Species)
  std::uint8_t status = 0U;  ///< PropagationStatus of the fit
  std::uint8_t failedPass = 0U;
  float chi2 = 0.F; ///< the inward pass's
  std::int16_t ndf = 0;
  float m2LogL = std::numeric_limits<float>::quiet_NaN();        ///< the inward pass's -2 ln L
  float outwardM2LogL = std::numeric_limits<float>::quiet_NaN(); ///< the outward pass's
  /// The fitted |p| on the innermost row [GeV]: |q| / |q/p|, the electron's
  /// the mode of its mixture's momentum density.
  float fittedP = std::numeric_limits<float>::quiet_NaN();

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
    visit("clusters", clusters);
    visit("layers", layers);
    visit("productionOffset", productionOffset);
    visit("p", p);
    visit("species", species);
    visit("status", status);
    visit("failedPass", failedPass);
    visit("chi2", chi2);
    visit("ndf", ndf);
    visit("m2LogL", m2LogL);
    visit("outwardM2LogL", outwardM2LogL);
    visit("fittedP", fittedP);
  }
};

template <DetectorLayout Det>
class SpeciesLikelihoodStudy
{
 public:
  explicit SpeciesLikelihoodStudy(RecordWriter& writer) : mWriter(&writer) { writer.begin<SpeciesLikelihoodRecord>(); }

  /// \a chains of \a event, from \a source.
  void process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source);

  [[nodiscard]] std::vector<std::string> summary() const
  {
    return {std::format("species likelihood study: {} chains x {} species, {} records into {}; {} fits failed", mChains, kNSpecies, mRecords, mWriter->path(), mFailed)};
  }

 private:
  RecordWriter* mWriter = nullptr;
  std::uint64_t mChains = 0U;
  std::uint64_t mRecords = 0U;
  std::uint64_t mFailed = 0U;
};

template <DetectorLayout Det>
void SpeciesLikelihoodStudy<Det>::process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source)
{
  const auto& frame = *event.frame;
  const auto& geometry = *event.geometry;
  constexpr o2::gtrk::detail::LayerMaterial material{};
  std::vector<std::vector<SpeciesLikelihoodRecord>> records(chains.size());
  o2::gtrk::detail::parallelFor(*event.arena, 0U, chains.size(), 16, [&](std::size_t begin, std::size_t end) {
    for (std::size_t c = begin; c < end; ++c) {
      const auto& chain = chains[c];
      const auto& truth = event.particles->at(chain.particle);
      const auto identity = identityOf(chain, event.timeframe, source);
      const auto count = static_cast<std::int32_t>(chain.rows.size());
      SpeciesLikelihoodRecord base{
        .timeframe = identity.timeframe,
        .source = identity.source,
        .track = identity.track,
        .particle = identity.particle,
        .leg = identity.leg,
        .pdg = truth.pdg,
        .charge = static_cast<std::int8_t>(std::lround(particleCharge(truth.pdg))),
        .origin = static_cast<std::uint8_t>(truth.origin),
        .eta = truth.eta,
        .clusters = static_cast<std::uint8_t>(std::min(count, 255)),
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
      if (event.hits != nullptr) {
        const auto m = frame.measurement(chain.rows.front());
        if (const auto* at = event.hits->crossing(chain.particle, m); at != nullptr) {
          base.p = static_cast<float>(std::hypot(at->momentum[0], at->momentum[1], at->momentum[2]));
        }
      }
      const RowChain<Det> rows{.frame = &frame, .rows = chain.rows.data()};
      auto& out = records[c];
      for (std::int32_t s = 0; s < kNSpecies; ++s) {
        const auto species = static_cast<Species>(s);
        auto prop = event.propagator;
        prop.setSpecies(species);
        SpeciesLikelihoodRecord record = base;
        record.species = static_cast<std::uint8_t>(s);
        const auto take = [&](const auto& gsf, Scalar momentum) {
          record.status = static_cast<std::uint8_t>(gsf.fit.status);
          record.failedPass = static_cast<std::uint8_t>(gsf.fit.failedPass);
          record.chi2 = gsf.fit.chi2;
          record.ndf = gsf.fit.ndf;
          record.outwardM2LogL = gsf.outwardMinusTwoLogLikelihood;
          if (gsf.fit.status == PropagationStatus::kOk) {
            record.m2LogL = gsf.minusTwoLogLikelihood;
            record.fittedP = momentum;
          }
        };
        if (species == Species::kElectron) {
          const auto gsf = fitGaussianSum(prop, geometry, material, rows, count);
          take(gsf, gsf.momentumMode);
        } else {
          const auto gsf = fitGaussianSum<NoRadiationLoss>(prop, geometry, material, rows, count);
          take(gsf, prop.momentum(gsf.fit.state.qOverP()));
        }
        out.push_back(record);
      }
    }
  });
  for (const auto& one : records) {
    for (const auto& record : one) {
      mWriter->write(record);
      mFailed += record.status != static_cast<std::uint8_t>(PropagationStatus::kOk) ? 1U : 0U;
    }
    mRecords += one.size();
  }
  mChains += chains.size();
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_SPECIESLIKELIHOODSTUDY_H_
