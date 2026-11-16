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

/// \file testGSFTrackFitter.cxx
/// \brief The Gaussian-sum filter: its Kalman limit, its closure against an enumeration of every loss path.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK GSF Track Fitter
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Tracking/GSFTrackFitter.h"
#include "GenericTracking/Tracking/KFTrackFitter.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "ElectronChain.h"
#include "ToyBarrel4Disc8.h"
#include "UniformMaterial.h"

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstring>
#include <random>
#include <vector>

namespace o2::gtrk
{
namespace
{

using test::ElectronChain;
using test::medianWidth;

/// No radiation: one component, z = 1, no variance.
struct NoLoss {
  static constexpr int kComponents = 1;
  static constexpr Scalar kMaxX2X0 = FLT_MAX;
  static int mixture(Scalar /*x2X0*/, BetheHeitlerComponent (&components)[kComponents])
  {
    components[0] = {.weight = 1.F, .mean = 1.F, .variance = 0.F};
    return 1;
  }
};

/// A loss of two components whatever the thickness: small and frequent, or
/// large and rare.
struct TwoLosses {
  static constexpr int kComponents = 2;
  static constexpr Scalar kMaxX2X0 = FLT_MAX;
  static constexpr BetheHeitlerComponent kLosses[kComponents] = {{.weight = 0.7F, .mean = 0.97F, .variance = 1e-4F},
                                                                 {.weight = 0.3F, .mean = 0.7F, .variance = 4e-3F}};
  static int mixture(Scalar /*x2X0*/, BetheHeitlerComponent (&components)[kComponents])
  {
    components[0] = kLosses[0];
    components[1] = kLosses[1];
    return kComponents;
  }
};

/// A chain over the first \a N barrel layers of ToyBarrel4Disc8, smeared, on
/// uniform material: generated without radiation (the closure and the
/// Kalman limit do not need the chain to follow the model).
template <std::size_t N>
struct Chain {
  using Mixed = detectors::ToyBarrel4Disc8;
  Geometry<Mixed> geometry = toyGeometry<Mixed>();
  TrackingGeometry<Mixed> tracking{geometry};
  Propagator propagator;
  std::array<Measurement, N> measurements;
  bool built = false;

  Chain(Scalar qOverP, std::mt19937& rng)
  {
    propagator.massGeV = 0.000511F;
    BoundTrackState state{};
    state.par[kPhi] = 0.12F;
    state.par[kTheta] = 1.4F;
    state.par[kQOverP] = qOverP;
    state.surface = firstSensorOfLayer<Mixed>(0);
    std::normal_distribution<double> gauss;
    constexpr Scalar kSigma = 5e-4F;
    for (std::size_t layer = 0; layer < N; ++layer) {
      if (layer > 0) {
        const auto& from = tracking.surface(state.surface);
        bool crossed = false;
        for (SensorID id = firstSensorOfLayer<Mixed>(static_cast<LayerID>(layer)); id < firstSensorOfLayer<Mixed>(static_cast<LayerID>(layer + 1)) && !crossed; ++id) {
          BoundTrackState probe = state;
          if (propagate(propagator, probe, from, tracking.surface(id), PropagationDirection::kForward).status == PropagationStatus::kOk) {
            state = probe;
            crossed = true;
          }
        }
        if (!crossed) {
          return;
        }
      }
      auto& m = measurements[layer];
      m = Measurement{};
      m.surface = state.surface;
      m.loc0 = state.loc0() + static_cast<Scalar>(gauss(rng)) * kSigma;
      m.loc1 = state.loc1() + static_cast<Scalar>(gauss(rng)) * kSigma;
      m.varLoc0 = kSigma * kSigma;
      m.varLoc1 = kSigma * kSigma;
      m.clusterIndex = static_cast<std::uint32_t>(layer);
    }
    built = true;
  }
};

/// What one configuration of fits gives against the truth on the vertex
/// end (the inward pass on the first plane).
struct Coverage {
  int fits = 0;
  int failed = 0;       ///< the filter failed, the Kalman fit did not
  int kalmanFailed = 0; ///< both may: skipped
  std::vector<double> pit;
  std::vector<double> outwardPit;  ///< the outward pass's on the last plane
  std::vector<double> kalmanRatio; ///< p / p_true of the Kalman fit
  std::vector<double> modeRatio;   ///< of the momentum's mode
  std::vector<double> meanRatio;   ///< of the collapsed mean q/p
  std::vector<double> collapsedPull;
  std::vector<double> outwardKalmanRatio; ///< on the last plane
  std::vector<double> outwardModeRatio;
  std::vector<double> innerLoss; ///< radiated between rows 0 and 2: the outward pass's start
  std::vector<double> outerLoss; ///< between rows 5 and 7: the inward pass's start
  double seconds = 0.;
  double kalmanSeconds = 0.;

  static double meanOf(const std::vector<double>& values)
  {
    double sum = 0.;
    for (const double p : values) {
      sum += p;
    }
    return sum / static_cast<double>(values.size());
  }
  /// The share outside [0.05, 0.95]: 10 % for a calibrated marginal.
  static double tailsOf(const std::vector<double>& values)
  {
    double n = 0.;
    for (const double p : values) {
      n += (p < 0.05 || p > 0.95) ? 1. : 0.;
    }
    return n / static_cast<double>(values.size());
  }
  /// The share of momenta more than 10 % off, either way.
  static double offOf(const std::vector<double>& ratios)
  {
    double n = 0.;
    for (const double r : ratios) {
      n += std::abs(r - 1.) > 0.1 ? 1. : 0.;
    }
    return n / static_cast<double>(ratios.size());
  }
  /// The share of momenta more than 10 % low.
  static double lowOf(const std::vector<double>& ratios)
  {
    double n = 0.;
    for (const double r : ratios) {
      n += r < 0.9 ? 1. : 0.;
    }
    return n / static_cast<double>(ratios.size());
  }
};

template <typename Loss = BetheHeitlerLoss, int K = Loss::kComponents>
Coverage coverage(double pGeV, const UniformMaterial& material, bool radiate, int tracks, std::uint32_t seed)
{
  ElectronChain chain;
  std::mt19937 rng(seed);
  Coverage out;
  for (int track = 0; track < tracks; ++track) {
    const auto qOverP = static_cast<Scalar>((track % 2 == 0 ? 1. : -1.) / pGeV);
    while (!chain.generate(qOverP, material, radiate, rng)) {
    }
    const auto begin = std::chrono::steady_clock::now();
    const auto gsf = fitGaussianSum<Loss, K>(chain.propagator, chain.tracking, material, chain.measurements.data(), static_cast<int>(ElectronChain::kLayers));
    out.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    const auto kalmanBegin = std::chrono::steady_clock::now();
    const auto kalman = fitMeasurements(chain.propagator, chain.tracking, material, chain.measurements.data(), static_cast<int>(ElectronChain::kLayers));
    out.kalmanSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - kalmanBegin).count();
    if (kalman.status != PropagationStatus::kOk) {
      ++out.kalmanFailed; // a track the fit's propagation cannot follow (a curler's turning point)
      continue;
    }
    if (gsf.fit.status != PropagationStatus::kOk) {
      ++out.failed;
      continue;
    }
    ++out.fits;
    const double u = chain.truth[0].qOverP();
    const double p = 1. / std::abs(u);
    out.pit.push_back(marginalCdf(gsf.mixture, u));
    out.outwardPit.push_back(marginalCdf(gsf.outwardMixture, chain.truth[ElectronChain::kLayers - 1].qOverP()));
    out.kalmanRatio.push_back(1. / std::abs(static_cast<double>(kalman.state.qOverP())) / p);
    out.modeRatio.push_back(static_cast<double>(gsf.momentumMode) / p);
    out.meanRatio.push_back(1. / std::abs(static_cast<double>(gsf.fit.state.qOverP())) / p);
    out.collapsedPull.push_back((static_cast<double>(gsf.fit.state.qOverP()) - u) / std::sqrt(static_cast<double>(gsf.fit.state.variance(kQOverP))));
    const double pLast = 1. / std::abs(static_cast<double>(chain.truth[ElectronChain::kLayers - 1].qOverP()));
    out.outwardKalmanRatio.push_back(1. / std::abs(static_cast<double>(kalman.outwardState.qOverP())) / pLast);
    out.outwardModeRatio.push_back(static_cast<double>(gsf.outwardMomentumMode) / pLast);
    out.innerLoss.push_back(chain.radiatedBetween(0, 2));
    out.outerLoss.push_back(chain.radiatedBetween(ElectronChain::kLayers - 3, ElectronChain::kLayers - 1));
  }
  return out;
}

std::string describe(const std::string& label, const Coverage& c)
{
  const auto [kalman, kalmanWidth] = medianWidth(c.kalmanRatio);
  const auto [mode, modeWidth] = medianWidth(c.modeRatio);
  const auto [mean, meanWidth] = medianWidth(c.meanRatio);
  return fmt::format("{:<34} {:5d} {:2d}/{:<2d} {:.3f} {:5.1f} %  {:.3f} {:5.1f} %   {:.3f} {:.3f} {:5.1f} %   {:.3f} {:.3f} {:5.1f} %   {:.3f} {:.3f}   {:6.3f} {:6.3f}",
                     label, c.fits, c.failed, c.kalmanFailed, Coverage::meanOf(c.pit), 100. * Coverage::tailsOf(c.pit), Coverage::meanOf(c.outwardPit), 100. * Coverage::tailsOf(c.outwardPit),
                     kalman, kalmanWidth, 100. * Coverage::lowOf(c.kalmanRatio), mode, modeWidth, 100. * Coverage::lowOf(c.modeRatio), mean, meanWidth,
                     1e3 * c.seconds / std::max(1, c.fits), 1e3 * c.kalmanSeconds / std::max(1, c.fits));
}

bool sameBits(const Scalar* a, const Scalar* b, int n) { return std::memcmp(a, b, sizeof(Scalar) * static_cast<std::size_t>(n)) == 0; }

} // namespace

/// One component that never loses: the filter is the Kalman fit, both
/// passes, to the bit -- the same start, the same steps, the same arithmetic.
BOOST_AUTO_TEST_CASE(WithoutLossItIsTheKalmanFit)
{
  std::mt19937 rng(105);
  const UniformMaterial material{.x2X0PerCm = 2e-3F, .rho = 0.3F};
  int fits = 0;
  for (int trial = 0; trial < 50; ++trial) {
    const Chain<4> chain(trial % 2 == 0 ? 3.F : -5.F, rng);
    BOOST_REQUIRE(chain.built);
    const auto kalman = fitMeasurements(chain.propagator, chain.tracking, material, chain.measurements.data(), 4);
    const auto gsf = fitGaussianSum<NoLoss>(chain.propagator, chain.tracking, material, chain.measurements.data(), 4);
    BOOST_REQUIRE(kalman.status == PropagationStatus::kOk);
    BOOST_REQUIRE(gsf.fit.status == PropagationStatus::kOk);
    BOOST_TEST(gsf.outwardMixture.size == 1);
    BOOST_TEST(gsf.mixture.size == 1);
    BOOST_TEST(sameBits(gsf.fit.outwardState.par, kalman.outwardState.par, kNBoundParams));
    BOOST_TEST(sameBits(gsf.fit.outwardState.cov, kalman.outwardState.cov, kBoundCovSize));
    BOOST_TEST(sameBits(gsf.fit.state.par, kalman.state.par, kNBoundParams));
    BOOST_TEST(sameBits(gsf.fit.state.cov, kalman.state.cov, kBoundCovSize));
    BOOST_TEST(gsf.fit.outwardChi2 == kalman.outwardChi2);
    BOOST_TEST(gsf.fit.chi2 == kalman.chi2);
    BOOST_TEST(gsf.fit.ndf == kalman.ndf);
    ++fits;
  }
  BOOST_TEST(fits == 50);
}

/// Every loss path enumerated: each plane the filter branches at -- the two
/// the start leaves, the one after it and the last -- takes one chosen loss
/// component, each path is filtered alone with the same steps (the diffuse
/// start's, then the filter's; its reference q/p the seed's scaled with each
/// loss), and its weight is the product of its loss
/// weights and its measurements' likelihoods (the start's diffuse one). The
/// unreduced filter (capacity for every path) holds exactly these paths, in
/// the order of the branching, with the normalised weights.
BOOST_AUTO_TEST_CASE(UnreducedItIsEveryLossPath)
{
  constexpr int kLayers = 4;     // ToyBarrel4Disc8's barrel
  constexpr int kBranchings = 4; // every plane
  constexpr int kPaths = 16;
  constexpr auto kRule = CrossingRule::kFirstAhead;
  std::mt19937 rng(1050);
  const UniformMaterial material{.x2X0PerCm = 2e-3F, .rho = 0.3F};
  double worstState = 0.;
  double worstWeight = 0.;
  for (int trial = 0; trial < 20; ++trial) {
    const Chain<kLayers> chain(trial % 2 == 0 ? 2.F : -4.F, rng);
    BOOST_REQUIRE(chain.built);
    const auto* measurements = chain.measurements.data();
    const auto gsf = fitGaussianSum<TwoLosses, kPaths>(chain.propagator, chain.tracking, material, measurements, kLayers);
    BOOST_REQUIRE(gsf.fit.status == PropagationStatus::kOk);
    const auto seed = makeChainSeed(chain.propagator, chain.tracking, measurements, kLayers);
    const auto span = detail::spanMaterial(chain.tracking, material, measurements, kLayers);
    for (const auto direction : {PropagationDirection::kForward, PropagationDirection::kBackward}) {
      const bool outward = direction == PropagationDirection::kForward;
      const auto& mixture = outward ? gsf.outwardMixture : gsf.mixture;
      BOOST_REQUIRE(mixture.size == kPaths);
      const int first = outward ? 0 : kLayers - 1;
      const int increment = outward ? 1 : -1;
      double logWeight[kPaths];
      BoundTrackState state[kPaths];
      CovarianceFactors factors[kPaths];
      for (int path = 0; path < kPaths; ++path) {
        const auto lossOf = [&](int b) { return TwoLosses::kLosses[(path >> (kBranchings - 1 - b)) & 1]; };
        state[path] = outward ? seed : gsf.fit.outwardState;
        logWeight[path] = 0.;
        Scalar chi2 = 0.F;
        // The pass's reference (PassReference), scaled with each loss.
        detail::PassReference reference{.qOverP = state[path].qOverP()};
        const auto scaleReference = [&](const BetheHeitlerComponent& loss) {
          reference.qOverP *= betheHeitlerScale(loss, direction);
          reference.lossScale *= betheHeitlerScale(loss, direction);
        };
        const auto referenceTo = [&](int index) { return detail::referenceFor(reference, kLayers, direction, [&] { return detail::leaveOneOutQOverP(chain.propagator, chain.tracking, measurements, kLayers, index); }); };
        detail::DiffuseFactors diffuse;
        BOOST_REQUIRE(detail::startDiffuse(measurements[first], state[path], diffuse));
        int plane = first;
        for (int b = 0; b < detail::kPassStart - 1; ++b, plane += increment) {
          const auto loss = lossOf(b);
          logWeight[path] += std::log(static_cast<double>(loss.weight));
          BOOST_REQUIRE(detail::gsf::applyBetheHeitler(loss, direction, state[path], diffuse));
          scaleReference(loss);
          Scalar density = 0.F;
          const Scalar stepReference = referenceTo(plane + increment);
          Scalar carried = stepReference;
          BOOST_REQUIRE((detail::diffuseStepBody<kRule, true>(chain.propagator, chain.tracking, span, measurements[plane + increment], direction, state[path], diffuse, chi2, carried, &density) == PropagationStatus::kOk));
          detail::advanceReference(reference, kLayers, direction, stepReference, carried);
          logWeight[path] -= 0.5 * static_cast<double>(density);
        }
        BOOST_REQUIRE(diffuse.rank == 0);
        factors[path] = diffuse.star;
        for (int b = detail::kPassStart - 1; b < kBranchings; ++b, plane += increment) {
          const auto loss = lossOf(b);
          logWeight[path] += std::log(static_cast<double>(loss.weight));
          BOOST_REQUIRE(applyBetheHeitler(loss, direction, state[path], factors[path]));
          scaleReference(loss);
          if (b + 1 < kBranchings) {
            Scalar logLikelihood = 0.F;
            const Scalar stepReference = referenceTo(plane + increment);
            Scalar carried = stepReference;
            BOOST_REQUIRE((detail::filterStepBody<false, kRule, true>(chain.propagator, chain.tracking, span, measurements[plane + increment], direction, state[path], factors[path], chi2, carried, 1.F, &logLikelihood) == PropagationStatus::kOk));
            detail::advanceReference(reference, kLayers, direction, stepReference, carried);
            logWeight[path] += logLikelihood;
          }
        }
        covarianceFromFactors(factors[path], state[path].cov);
      }
      double maxLog = -1e300;
      for (const double l : logWeight) {
        maxLog = std::max(maxLog, l);
      }
      double sum = 0.;
      for (const double l : logWeight) {
        sum += std::exp(l - maxLog);
      }
      for (int path = 0; path < kPaths; ++path) {
        const double expected = std::exp(logWeight[path] - maxLog) / sum;
        const double got = std::exp(static_cast<double>(mixture.components[path].logWeight));
        worstWeight = std::max(worstWeight, std::abs(got - expected));
        for (int i = 0; i < kNBoundParams; ++i) {
          const double sigma = std::sqrt(static_cast<double>(state[path].variance(i)));
          worstState = std::max(worstState, std::abs(static_cast<double>(mixture.components[path].state.par[i] - state[path].par[i])) / sigma);
        }
        for (int i = 0; i < kBoundCovSize; ++i) {
          const double scale = std::abs(static_cast<double>(state[path].cov[i])) + 1e-30;
          worstState = std::max(worstState, std::abs(static_cast<double>(mixture.components[path].state.cov[i] - state[path].cov[i])) / scale);
        }
      }
    }
  }
  BOOST_TEST_MESSAGE(fmt::format("every loss path: worst weight difference {:.2e}, worst state difference {:.2e} (sigma or relative)", worstWeight, worstState));
  BOOST_TEST(worstWeight < 1e-5);
  BOOST_TEST(worstState < 1e-5);
}

/// Sampled Bethe-Heitler losses, independent of the mixture: each pass's q/p
/// marginal covers the true q/p on its end plane (the probability integral
/// transform flat: its mean a half, a tenth of it outside 5-95 %). The
/// momentum's posterior is skewed, so neither its mode nor its mean is
/// unbiased; both are reported beside the Kalman fit's momentum, which takes
/// no radiation and which the mode beats on thick material (asserted: fewer
/// momenta more than 10 % low, a narrower core). Without radiation the
/// filter is as calibrated. One component (K 1) is not: reported only, with
/// the convergence in K and the cost per fit.
BOOST_AUTO_TEST_CASE(TheMixtureCoversSampledLosses)
{
  const UniformMaterial thin{.x2X0PerCm = 2.5e-3F, .rho = 0.055F}; // silicon's ratio, 0.5 % X0 per gap
  const UniformMaterial thick{.x2X0PerCm = 1e-2F, .rho = 0.22F};   // 2 % X0 per gap
  constexpr int kTracks = 1000;
  struct Row {
    std::string label;
    Coverage coverage;
    bool calibrated; ///< asserted: flat PIT on both ends
  };
  std::vector<Row> rows;
  for (const double p : {0.2, 0.5, 2.0}) {
    rows.push_back({fmt::format("thin, p {:.1f} GeV", p), coverage<>(p, thin, true, kTracks, 105), true});
  }
  for (const double p : {0.2, 0.5}) {
    rows.push_back({fmt::format("thick, p {:.1f} GeV", p), coverage<>(p, thick, true, kTracks, 106), true});
  }
  rows.push_back({"thick, p 0.5 GeV, K 1", coverage<BetheHeitlerLoss, 1>(0.5, thick, true, kTracks, 106), false});
  rows.push_back({"thick, p 0.5 GeV, K 4", coverage<BetheHeitlerLoss, 4>(0.5, thick, true, kTracks, 106), false});
  rows.push_back({"thick, p 0.5 GeV, K 24", coverage<BetheHeitlerLoss, 24>(0.5, thick, true, kTracks, 106), false});
  // Who owns which loss: the ionisation alone (the no-loss model on tracks without radiation, the Kalman fit's
  // own q/p coverage), the radiation alone (no ionisation), and a radiation the tracks do not have.
  rows.push_back({"thick, p 0.5 GeV, ionisation only", coverage<NoLoss>(0.5, thick, false, kTracks, 107), false});
  const UniformMaterial thickNoIonisation{.x2X0PerCm = thick.x2X0PerCm, .rho = 0.F};
  rows.push_back({"thick, p 0.5 GeV, radiation only", coverage<>(0.5, thickNoIonisation, true, kTracks, 108), true});
  rows.push_back({"thick, p 0.5 GeV, tracks without radiation", coverage<>(0.5, thick, false, kTracks, 107), false});
  BOOST_TEST_MESSAGE("                                    fits fail/KF inward PIT       outward PIT      Kalman p/p_true          mode p/p_true            mean p/p_true   ms/fit");
  BOOST_TEST_MESSAGE("                                                 mean  5-95 %     mean  5-95 %     median width  <0.9       median width  <0.9       median width    GSF  Kalman");
  for (const auto& row : rows) {
    BOOST_TEST_MESSAGE(describe(row.label, row.coverage));
    if (!row.calibrated) {
      continue;
    }
    const auto& c = row.coverage;
    BOOST_TEST_CONTEXT(row.label)
    {
      BOOST_TEST(c.failed == 0);
      // 1000 tracks: the PIT mean's error 0.009, the tail share's 0.95 %.
      BOOST_TEST(std::abs(Coverage::meanOf(c.pit) - 0.5) < 0.03);
      BOOST_TEST(std::abs(Coverage::meanOf(c.outwardPit) - 0.5) < 0.03);
      BOOST_TEST(std::abs(Coverage::tailsOf(c.pit) - 0.1) < 0.035);
      BOOST_TEST(std::abs(Coverage::tailsOf(c.outwardPit) - 0.1) < 0.035);
      if (row.label.starts_with("thick")) {
        BOOST_TEST(Coverage::lowOf(c.modeRatio) < 0.8 * Coverage::lowOf(c.kalmanRatio));
        BOOST_TEST(medianWidth(c.modeRatio).second < 0.9 * medianWidth(c.kalmanRatio).second);
      }
    }
  }
}

/// (iv) The start: tracks that radiated a fifth of their momentum or more
/// between the rows a pass's start takes (rows 0-2 the outward pass's, near
/// the vertex; rows 5-7 the inward pass's), reported apart. Selected by
/// their true loss, so the probability integral transforms are not flat for
/// them (they sit in the posterior's low tail) and are reported only; asserted:
/// the momentum at each end beats the Kalman fit's on these tracks (fewer more
/// than 10 % off: the Kalman fit errs low at the vertex end and high at the
/// last plane). A loss inside a pass's own start leaves that pass calibrated
/// at its far end (the start branches); a loss on the last rows a pass takes
/// is beyond what its remaining measurements show (the inward pass at the
/// vertex end, after a loss between the innermost rows): both fits keep most
/// of it there.
BOOST_AUTO_TEST_CASE(LargeLossesInsideAStart)
{
  const UniformMaterial thick{.x2X0PerCm = 1e-2F, .rho = 0.22F};
  const auto all = coverage<>(0.5, thick, true, 4000, 109);
  const auto subset = [&](const std::vector<double>& loss, const std::vector<double>& values) {
    std::vector<double> out;
    for (std::size_t i = 0; i < values.size(); ++i) {
      if (loss[i] > 0.2) {
        out.push_back(values[i]);
      }
    }
    return out;
  };
  BOOST_TEST_MESSAGE("                       tracks  inward PIT      outward PIT     vertex end p/p_true: Kalman / GSF mode       last plane p/p_true: Kalman / GSF mode");
  BOOST_TEST_MESSAGE("                               mean 5-95 %     mean 5-95 %     median width >10 % off                       median width >10 % off");
  for (const auto& [label, loss] : {std::pair<const char*, const std::vector<double>*>{"all", nullptr}, {"lost > 20 % in rows 0-2", &all.innerLoss}, {"lost > 20 % in rows 5-7", &all.outerLoss}}) {
    const auto pick = [&, selected = loss](const std::vector<double>& values) { return selected == nullptr ? values : subset(*selected, values); };
    const auto pit = pick(all.pit);
    const auto outwardPit = pick(all.outwardPit);
    const auto kalman = pick(all.kalmanRatio);
    const auto mode = pick(all.modeRatio);
    const auto outwardKalman = pick(all.outwardKalmanRatio);
    const auto outwardMode = pick(all.outwardModeRatio);
    const auto row = [&](const std::vector<double>& r) {
      const auto [median, width] = medianWidth(r);
      return fmt::format("{:.3f} {:.3f} {:4.1f} %", median, width, 100. * Coverage::offOf(r));
    };
    BOOST_TEST_MESSAGE(fmt::format("{:<24} {:5d}  {:.3f} {:4.1f} %    {:.3f} {:4.1f} %    {} / {}    {} / {}", label, pit.size(), Coverage::meanOf(pit), 100. * Coverage::tailsOf(pit),
                                   Coverage::meanOf(outwardPit), 100. * Coverage::tailsOf(outwardPit), row(kalman), row(mode), row(outwardKalman), row(outwardMode)));
    BOOST_TEST_CONTEXT(label)
    {
      BOOST_TEST(pit.size() > 100U);
      BOOST_TEST(Coverage::offOf(mode) < Coverage::offOf(kalman));
      BOOST_TEST(Coverage::offOf(outwardMode) < Coverage::offOf(outwardKalman));
    }
  }
}

/// The reduction: K = 1, 4, the default (6, one component per loss
/// component) and 12 against K = 24 on the same tracks, per track, in the
/// reference's standard deviations (the collapsed vertex-end state) and in
/// the momentum's mode, relative. The merge cost reads q/p alone: what it
/// gives up about the other parameters shows here. Asserted on the default.
BOOST_AUTO_TEST_CASE(TheReductionConverges)
{
  const UniformMaterial thick{.x2X0PerCm = 1e-2F, .rho = 0.22F};
  ElectronChain chain;
  std::mt19937 rng(110);
  constexpr int kTracks = 300;
  constexpr int kLayers = static_cast<int>(ElectronChain::kLayers);
  std::array<std::array<std::vector<double>, kNBoundParams + 1>, 4> difference; // per K: the parameters, the mode
  for (int track = 0; track < kTracks; ++track) {
    const Scalar qOverP = (track % 2 == 0 ? 2.F : -2.F);
    while (!chain.generate(qOverP, thick, true, rng)) {
    }
    const auto* measurements = chain.measurements.data();
    const auto reference = fitGaussianSum<BetheHeitlerLoss, 24>(chain.propagator, chain.tracking, thick, measurements, kLayers);
    if (reference.fit.status != PropagationStatus::kOk) {
      continue;
    }
    const auto compare = [&](int slot, const auto& gsf) {
      BOOST_REQUIRE(gsf.fit.status == PropagationStatus::kOk);
      for (int i = 0; i < kNBoundParams; ++i) {
        double d = static_cast<double>(gsf.fit.state.par[i]) - static_cast<double>(reference.fit.state.par[i]);
        if (i == kPhi) {
          math_utils::detail::bringToPMPi<double>(d);
        }
        difference[slot][i].push_back(std::abs(d) / std::sqrt(static_cast<double>(reference.fit.state.variance(i))));
      }
      difference[slot][kNBoundParams].push_back(std::abs(static_cast<double>(gsf.momentumMode / reference.momentumMode) - 1.));
    };
    compare(0, fitGaussianSum<BetheHeitlerLoss, 1>(chain.propagator, chain.tracking, thick, measurements, kLayers));
    compare(1, fitGaussianSum<BetheHeitlerLoss, 4>(chain.propagator, chain.tracking, thick, measurements, kLayers));
    compare(2, fitGaussianSum(chain.propagator, chain.tracking, thick, measurements, kLayers));
    compare(3, fitGaussianSum<BetheHeitlerLoss, 12>(chain.propagator, chain.tracking, thick, measurements, kLayers));
  }
  const auto quantile = [](std::vector<double> v, double q) {
    std::ranges::sort(v);
    return v[static_cast<std::size_t>(q * static_cast<double>(v.size() - 1))];
  };
  BOOST_TEST_MESSAGE("against K 24, |difference| / sigma (median, 95 %) of loc0 loc1 phi theta q/p; the mode, relative");
  const int ks[] = {1, 4, BetheHeitlerLoss::kComponents, 12};
  for (int slot = 0; slot < 4; ++slot) {
    std::string line = fmt::format("K {:2d}:", ks[slot]);
    for (int i = 0; i <= kNBoundParams; ++i) {
      line += fmt::format("  {:.3f} {:.3f}", quantile(difference[slot][i], 0.5), quantile(difference[slot][i], 0.95));
    }
    BOOST_TEST_MESSAGE(line);
  }
  for (int i = 0; i < kNBoundParams; ++i) {
    BOOST_TEST_CONTEXT("the default K, parameter " << i) { BOOST_TEST(quantile(difference[2][i], 0.95) < 0.1); }
  }
  BOOST_TEST(quantile(difference[2][kNBoundParams], 0.95) < 0.02); // a fifth of the resolution here
}

} // namespace o2::gtrk
