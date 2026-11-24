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

/// \file GSFStudy.h
/// \brief The Gaussian-sum filter against the Kalman fit and the MC, both
///        under the electron hypothesis, on chains (StudyChains.h).
///
/// Of each chain it takes the first row on each layer, at most
/// `maxChainLayers` rows, and at least `kPassStart` of them, of a charged
/// particle. Every prefix of at least `kPassStart` rows is fitted, the
/// outward pass's end held to the MC momentum where the particle enters the
/// prefix's last row, and every such suffix, the inward pass's end held to
/// it at the suffix's first row: a record per row and pass, with each fit's
/// q/p and standard deviation, the mixture's mode and median, and the
/// probability integral transform of the true q/p under the mixture's
/// marginal. Then the production point (pass 2, k and layer -1): the refit
/// under the electron hypothesis (`refitWithHypothesis`) carried to the
/// perigee plane about the MC production point, and the Kalman fit carried
/// alike, against the MC momentum at production, with the d0 pulls.

#ifndef GTRK_STUDY_GSFSTUDY_H_
#define GTRK_STUDY_GSFSTUDY_H_

#include "GenericTracking/Tracking/GSFTrackFitter.h"
#include "GenericTracking/Tracking/HypothesisRefit.h"
#include "GenericTracking/Tracking/KFTrackFitter.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/Species.h"
#include "GenericTracking/Study/RecordWriter.h"
#include "GenericTracking/Study/StateTruth.h"
#include "GenericTracking/Study/StudyChains.h"
#include "GenericTracking/Study/TruthEvent.h"

#include "CommonConstants/PhysicsConstants.h"

#include <algorithm>
#include <chrono>
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

/// One row and pass of a chain's fits (`GSFStudy`).
struct GSFRecord {
  static constexpr std::string_view kName = "gsf";
  std::uint32_t timeframe = 0U;
  std::uint8_t source = 0U; ///< ChainSource
  std::int32_t track = -1;  ///< the tracker's track, -1 for a truth chain
  std::uint64_t particle = 0U;
  std::uint8_t leg = 0U;
  std::int32_t pdg = 0;
  float pt = 0.F; ///< the truth's [GeV]
  float eta = 0.F;
  std::int32_t n = 0;      ///< the fitted rows
  std::int32_t k = -1;     ///< the row, -1 at the production point
  std::int32_t layer = -1; ///< its layer
  std::int32_t pass = 0;   ///< 0 outward, 1 inward, 2 at the production point
  double x2X0 = 0.;        ///< the row's plane's radiation length fraction; the segment to the perigee's at production
  double pTrue = 0.;       ///< the MC |p| [GeV]
  double qOverPTrue = 0.;
  double legLoss = 0.; ///< 1 - pTrue at the last row / pTrue at the first
  std::uint8_t kfOk = 0U;
  double kfQOverP = 0.;
  double kfSigma = 0.;
  double kfChi2 = 0.;
  std::uint8_t gsfOk = 0U;
  double gsfQOverP = 0.;
  double gsfSigma = 0.;
  double gsfModeP = 0.;
  double gsfMedianP = 0.;
  double gsfPit = 0.;
  std::int32_t gsfComponents = 0;
  double gsfM2LogL = 0.;
  double kfD0Pull = 0.;
  double gsfD0Pull = 0.;

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("source", source);
    visit("track", track);
    visit("particle", particle);
    visit("leg", leg);
    visit("pdg", pdg);
    visit("pt", pt);
    visit("eta", eta);
    visit("n", n);
    visit("k", k);
    visit("layer", layer);
    visit("pass", pass);
    visit("x2X0", x2X0);
    visit("pTrue", pTrue);
    visit("qOverPTrue", qOverPTrue);
    visit("legLoss", legLoss);
    visit("kfOk", kfOk);
    visit("kfQOverP", kfQOverP);
    visit("kfSigma", kfSigma);
    visit("kfChi2", kfChi2);
    visit("gsfOk", gsfOk);
    visit("gsfQOverP", gsfQOverP);
    visit("gsfSigma", gsfSigma);
    visit("gsfModeP", gsfModeP);
    visit("gsfMedianP", gsfMedianP);
    visit("gsfPit", gsfPit);
    visit("gsfComponents", gsfComponents);
    visit("gsfM2LogL", gsfM2LogL);
    visit("kfD0Pull", kfD0Pull);
    visit("gsfD0Pull", gsfD0Pull);
  }
};

template <DetectorLayout Det>
class GSFStudy
{
 public:
  explicit GSFStudy(RecordWriter& writer) : mWriter(&writer) { writer.begin<GSFRecord>(); }

  /// \a chains of \a event, from \a source.
  void process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source);

  /// The fits' count and their time.
  [[nodiscard]] std::vector<std::string> summary() const
  {
    const auto perFit = [&](double seconds) { return 1e3 * seconds / static_cast<double>(std::max<std::size_t>(mFits, 1U)); };
    return {std::format("GSF study: {} chains, {} records into {}; refused {} neutral, {} with fewer than {} layers; per fit (both passes, one thread) the Gaussian-sum fit {:.3f} ms, the Kalman fit {:.4f} ms", mChains, mRecords, mWriter->path(), mNeutral, mFewLayers, o2::gtrk::detail::kPassStart, perFit(mGsfSeconds), perFit(mKfSeconds))};
  }

 private:
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  /// The sign relating the state's q/p to the charge: q/p = sign x charge / p
  /// (Propagator's curvature convention).
  static constexpr double kChargeSign = 1.;

  /// One chain's rows as fitted, its truth, and its records.
  struct Leg {
    const Chain* chain = nullptr;
    std::vector<std::uint32_t> rows;
    int charge = 0;
    std::vector<double> pTrue;
    double legLoss = 0.;
    std::vector<GSFRecord> records;
    double gsfSeconds = 0.;
    double kfSeconds = 0.;
    std::size_t fits = 0U;
  };

  void fitLeg(Leg& leg, const TruthEvent<Det>& event, const Propagator& prop, const ChainIdentity& identity) const;

  RecordWriter* mWriter = nullptr;
  std::size_t mChains = 0U;
  std::size_t mNeutral = 0U;   ///< chains refused: no charge
  std::size_t mFewLayers = 0U; ///< chains refused: fewer layers than a pass needs
  std::size_t mRecords = 0U;
  std::size_t mFits = 0U;
  double mGsfSeconds = 0.;
  double mKfSeconds = 0.;
};

template <DetectorLayout Det>
void GSFStudy<Det>::process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source)
{
  Propagator prop = event.propagator;
  prop.massGeV = static_cast<Scalar>(o2::constants::physics::MassElectron);
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  // Each chain's rows and truth serially; then the fits in parallel, each
  // chain's records apart, written in the chains' order.
  std::vector<Leg> legs;
  for (const auto& chain : chains) {
    Leg leg{.chain = &chain};
    std::uint64_t layers = 0U;
    for (std::size_t k = 0; k < chain.rows.size() && leg.rows.size() < kMaxClusters; ++k) {
      const auto bit = std::uint64_t{1} << event.rows[chain.rows[k]].layer;
      if ((layers & bit) == 0U) {
        leg.rows.push_back(chain.rows[k]);
        layers |= bit;
      }
    }
    leg.charge = static_cast<int>(std::lround(particleCharge(event.particles->at(chain.particle).pdg)));
    if (leg.charge == 0) {
      ++mNeutral;
      continue;
    }
    if (static_cast<int>(leg.rows.size()) < o2::gtrk::detail::kPassStart) {
      ++mFewLayers;
      continue;
    }
    leg.pTrue.assign(leg.rows.size(), nan);
    for (std::size_t k = 0; k < leg.rows.size(); ++k) {
      const auto point = event.frame->measurement(leg.rows[k]);
      const auto* hit = event.hits == nullptr ? nullptr : event.hits->crossing(chain.particle, point);
      if (hit != nullptr) {
        leg.pTrue[k] = std::hypot(hit->momentum[0], hit->momentum[1], hit->momentum[2]);
      }
    }
    leg.legLoss = 1. - (leg.pTrue.back() / leg.pTrue.front());
    legs.push_back(std::move(leg));
  }
  o2::gtrk::detail::parallelFor(*event.arena, 0U, legs.size(), 4, [&](std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      fitLeg(legs[i], event, prop, identityOf(*legs[i].chain, event.timeframe, source));
    }
  });
  for (const auto& leg : legs) {
    for (const auto& record : leg.records) {
      mWriter->write(record);
    }
    mRecords += leg.records.size();
    mGsfSeconds += leg.gsfSeconds;
    mKfSeconds += leg.kfSeconds;
    mFits += leg.fits;
  }
  mChains += legs.size();
}

template <DetectorLayout Det>
void GSFStudy<Det>::fitLeg(Leg& leg, const TruthEvent<Det>& event, const Propagator& prop, const ChainIdentity& identity) const
{
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  const auto& truth = event.particles->at(identity.particle);
  const auto& geometry = *event.geometry;
  constexpr o2::gtrk::detail::LayerMaterial material{};
  const auto& frame = *event.frame;
  const auto& chain = leg.rows;
  const int n = static_cast<int>(chain.size());
  const RowChain<Det> whole{.frame = &frame, .rows = chain.data()};
  std::vector<double> x2X0(chain.size());
  for (int k = 0; k < n; ++k) {
    x2X0[static_cast<std::size_t>(k)] = o2::gtrk::detail::gsf::planeRadiation(geometry, whole, n, k);
  }
  GSFRecord base{.timeframe = identity.timeframe, .source = identity.source, .track = identity.track, .particle = identity.particle, .leg = identity.leg};
  base.pdg = truth.pdg;
  base.pt = truth.pt;
  base.eta = truth.eta;
  base.n = n;
  base.legLoss = leg.legLoss;
  for (int k = 0; k < n; ++k) {
    for (const auto direction : {PropagationDirection::kForward, PropagationDirection::kBackward}) {
      const bool outward = direction == PropagationDirection::kForward;
      const int first = outward ? 0 : k;
      const int count = outward ? k + 1 : n - k;
      const double pTrue = leg.pTrue[static_cast<std::size_t>(k)];
      if (count < o2::gtrk::detail::kPassStart || std::isnan(pTrue)) {
        continue;
      }
      const RowChain<Det> rows{.frame = &frame, .rows = chain.data() + first};
      const auto kfStart = std::chrono::steady_clock::now();
      const auto kalman = fitMeasurements(prop, geometry, material, rows, count);
      const auto gsfStart = std::chrono::steady_clock::now();
      const auto gsf = fitGaussianSum(prop, geometry, material, rows, count);
      const auto gsfEnd = std::chrono::steady_clock::now();
      leg.kfSeconds += std::chrono::duration<double>(gsfStart - kfStart).count();
      leg.gsfSeconds += std::chrono::duration<double>(gsfEnd - gsfStart).count();
      ++leg.fits;
      // A pass's end: the outward one stands if the inward one failed.
      const auto passOk = [&](const TrackFitResult& fit) { return outward ? fit.failedPass != FitPass::kOutward : fit.status == PropagationStatus::kOk; };
      const bool kfOk = passOk(kalman);
      const auto& kf = outward ? kalman.outwardState : kalman.state;
      const bool gsfOk = passOk(gsf.fit);
      const auto& gsfState = outward ? gsf.fit.outwardState : gsf.fit.state;
      const auto& mixture = outward ? gsf.outwardMixture : gsf.mixture;
      const double truthQOverP = static_cast<double>(leg.charge) * kChargeSign / pTrue;
      // The q/p marginal's distribution function; the PIT at the truth, the
      // median by bisection between the components' 8-sigma ends.
      const auto cdf = [&](double u) { return marginalCdf(mixture, u); };
      double pit = nan;
      double medianP = nan;
      if (gsfOk) {
        pit = cdf(truthQOverP);
        double low = std::numeric_limits<double>::max();
        double high = -std::numeric_limits<double>::max();
        for (int c = 0; c < mixture.size; ++c) {
          const auto& state = mixture.components[c].state;
          const double m = state.qOverP();
          const double sigma = std::sqrt(static_cast<double>(state.variance(kQOverP)));
          low = std::min(low, m - (8. * sigma));
          high = std::max(high, m + (8. * sigma));
        }
        for (int step = 0; step < 100 && high - low > 1e-12 * std::abs(high); ++step) {
          const double middle = 0.5 * (low + high);
          (cdf(middle) < 0.5 ? low : high) = middle;
        }
        medianP = 1. / std::abs(0.5 * (low + high));
      }
      GSFRecord record = base;
      record.k = k;
      record.layer = static_cast<std::int32_t>(event.rows[chain[static_cast<std::size_t>(k)]].layer);
      record.pass = outward ? 0 : 1;
      record.x2X0 = x2X0[static_cast<std::size_t>(k)];
      record.pTrue = pTrue;
      record.qOverPTrue = truthQOverP;
      record.kfOk = kfOk ? 1U : 0U;
      record.kfQOverP = kfOk ? static_cast<double>(kf.qOverP()) : nan;
      record.kfSigma = kfOk ? std::sqrt(static_cast<double>(kf.variance(kQOverP))) : nan;
      record.kfChi2 = static_cast<double>(outward ? kalman.outwardChi2 : kalman.chi2);
      record.gsfOk = gsfOk ? 1U : 0U;
      record.gsfQOverP = gsfOk ? static_cast<double>(gsfState.qOverP()) : nan;
      record.gsfSigma = gsfOk ? std::sqrt(static_cast<double>(gsfState.variance(kQOverP))) : nan;
      record.gsfModeP = gsfOk ? static_cast<double>(outward ? gsf.outwardMomentumMode : gsf.momentumMode) : nan;
      record.gsfMedianP = medianP;
      record.gsfPit = pit;
      record.gsfComponents = mixture.size;
      record.gsfM2LogL = static_cast<double>(outward ? gsf.outwardMinusTwoLogLikelihood : gsf.minusTwoLogLikelihood);
      record.kfD0Pull = nan;
      record.gsfD0Pull = nan;
      leg.records.push_back(record);
    }
  }
  // The production point: the refit under the electron hypothesis carried to
  // the perigee plane about the MC production point, and the Kalman fit
  // carried alike.
  const RowChain<Det> rows{.frame = &frame, .rows = chain.data()};
  const auto refit = refitWithHypothesis(prop, geometry, material, rows, n, Species::kElectron, truth.vx, truth.vy);
  const auto kalman = fitMeasurements(prop, geometry, material, rows, n);
  const auto& innermostSurface = geometry.planeOf(rows[0]);
  bool kfOk = kalman.status == PropagationStatus::kOk;
  BoundTrackState kf = kalman.state;
  double segmentX2X0 = nan;
  if (kfOk) {
    const auto plane = o2::gtrk::detail::hypothesis::perigeePlaneOf(prop, kf, innermostSurface, truth.vx, truth.vy);
    segmentX2X0 = o2::gtrk::detail::hypothesis::segmentRadiation(prop, geometry.material(), kf, innermostSurface, plane);
    kfOk = o2::gtrk::detail::hypothesis::carryStateInward(prop, geometry.material(), innermostSurface, plane, kf) == PropagationStatus::kOk;
  }
  const bool gsfOk = refit.status == PropagationStatus::kOk;
  const double pTrue = static_cast<double>(truth.pt) * std::cosh(static_cast<double>(truth.eta));
  const double truthQOverP = static_cast<double>(leg.charge) * kChargeSign / pTrue;
  const auto d0Pull = [](const BoundTrackState& state) { return static_cast<double>(state.loc0()) / std::sqrt(static_cast<double>(state.variance(kLoc0))); };
  GSFRecord record = base;
  record.k = -1;
  record.layer = -1;
  record.pass = 2;
  record.x2X0 = segmentX2X0;
  record.pTrue = pTrue;
  record.qOverPTrue = truthQOverP;
  record.kfOk = kfOk ? 1U : 0U;
  record.kfQOverP = kfOk ? static_cast<double>(kf.qOverP()) : nan;
  record.kfSigma = kfOk ? std::sqrt(static_cast<double>(kf.variance(kQOverP))) : nan;
  record.kfChi2 = static_cast<double>(kalman.chi2);
  record.gsfOk = gsfOk ? 1U : 0U;
  record.gsfQOverP = gsfOk ? static_cast<double>(refit.perigee.qOverP()) : nan;
  record.gsfSigma = gsfOk ? std::sqrt(static_cast<double>(refit.perigee.variance(kQOverP))) : nan;
  record.gsfModeP = gsfOk ? static_cast<double>(refit.perigeeMomentum) : nan;
  record.gsfMedianP = nan;
  record.gsfPit = gsfOk ? marginalCdf(refit.perigeeMixture, truthQOverP) : nan;
  record.gsfComponents = refit.perigeeMixture.size;
  record.gsfM2LogL = nan;
  record.kfD0Pull = kfOk ? d0Pull(kf) : nan;
  record.gsfD0Pull = gsfOk ? d0Pull(refit.perigee) : nan;
  leg.records.push_back(record);
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_GSFSTUDY_H_
