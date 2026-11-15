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
#include "ToyBarrel4Disc8.h"
#include "ToyBarrel8.h"
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

/// An electron through ToyBarrel8's eight layers (2-16 cm) on uniform
/// material, as the filter models it: on each step the scattering in two
/// kicks of half the chord's Highland theta0 (behind the sensor left, in
/// front of the sensor reached, the logarithm at the span's material), the
/// ionisation over the chord on reaching the next sensor (its mean and the
/// straggling the filter takes, Gaussian in q/p), and at each
/// plane the radiation of the half-gap before it, then of the half-gap after
/// it, each fraction z drawn from Bethe and Heitler's model itself (-ln z
/// Gamma-distributed with shape t / ln 2), not from the mixture. truth[i] is
/// the state on plane i after the half-gap before it: what the outward pass
/// ends with on the last plane and the inward one on the first. The chords
/// are those of the track without material.
struct ElectronChain {
  static constexpr std::size_t kLayers = 8;
  using Barrel = detectors::ToyBarrel8;
  Geometry<Barrel> geometry = toyGeometry<Barrel>();
  TrackingGeometry<Barrel> tracking{geometry};
  Propagator propagator;
  std::array<BoundTrackState, kLayers> truth;
  std::array<Measurement, kLayers> measurements;
  Scalar sigma = 5e-4F;

  ElectronChain() { propagator.massGeV = 0.000511F; }

  [[nodiscard]] bool cross(BoundTrackState& state, LayerID layer) const
  {
    const auto& from = tracking.surface(state.surface);
    for (SensorID id = firstSensorOfLayer<Barrel>(layer); id < firstSensorOfLayer<Barrel>(static_cast<LayerID>(layer + 1)); ++id) {
      BoundTrackState probe = state;
      if (propagate(propagator, probe, from, tracking.surface(id), PropagationDirection::kForward).status == PropagationStatus::kOk) {
        state = probe;
        return true;
      }
    }
    return false;
  }

  static Scalar fraction(double x2X0, std::mt19937& rng)
  {
    if (!(x2X0 > 0.)) {
      return 1.F;
    }
    std::gamma_distribution<double> gamma(x2X0 / std::numbers::ln2, 1.);
    return static_cast<Scalar>(std::exp(-gamma(rng)));
  }

  /// A track of \a qOverP; false if it does not reach every layer.
  bool generate(Scalar qOverP, const UniformMaterial& material, bool radiate, std::mt19937& rng)
  {
    BoundTrackState start{};
    start.par[kPhi] = 0.05F;
    start.par[kTheta] = 1.55F;
    start.par[kQOverP] = qOverP;
    start.surface = firstSensorOfLayer<Barrel>(0);
    // The chords without material.
    std::array<Scalar, kLayers> chord{};
    BoundTrackState state = start;
    Scalar previous[3];
    surfaceToGlobal(tracking.surface(state.surface), state.loc0(), state.loc1(), previous);
    for (std::size_t layer = 1; layer < kLayers; ++layer) {
      if (!cross(state, static_cast<LayerID>(layer))) {
        return false;
      }
      Scalar point[3];
      surfaceToGlobal(tracking.surface(state.surface), state.loc0(), state.loc1(), point);
      chord[layer] = std::hypot(point[0] - previous[0], point[1] - previous[1], point[2] - previous[2]);
      std::copy(point, point + 3, previous);
    }
    Scalar span = 0.F;
    for (const Scalar c : chord) {
      span += material.x2X0PerCm * c;
    }
    std::normal_distribution<double> gauss;
    const auto kick = [&](BoundTrackState& s, Scalar theta0) {
      const auto sinTheta = std::sin(s.theta());
      s.par[kTheta] += static_cast<Scalar>(gauss(rng)) * theta0;
      s.par[kPhi] += static_cast<Scalar>(gauss(rng)) * theta0 / sinTheta;
    };
    const auto halfGap = [&](std::size_t layer) { return 0.5 * static_cast<double>(material.x2X0PerCm * chord[layer]); };
    state = start;
    for (std::size_t layer = 0; layer < kLayers; ++layer) {
      if (layer > 0) {
        const Scalar theta0 = highlandTheta0(0.5F * material.x2X0PerCm * chord[layer], span, 1.F / std::abs(state.qOverP()), propagator.massGeV);
        kick(state, theta0);
        if (!cross(state, static_cast<LayerID>(layer))) {
          return false;
        }
        Scalar variance = 0.F;
        if (!detail::applyEnergyLoss(propagator, material.rho * chord[layer], PropagationDirection::kForward, state, variance)) {
          return false;
        }
        state.par[kQOverP] += static_cast<Scalar>(gauss(rng)) * std::sqrt(variance); // the straggling the filter takes
        kick(state, theta0);
        if (radiate) {
          state.par[kQOverP] /= fraction(halfGap(layer), rng);
        }
      }
      truth[layer] = state;
      auto& m = measurements[layer];
      m = Measurement{};
      m.surface = state.surface;
      m.loc0 = state.loc0() + static_cast<Scalar>(gauss(rng)) * sigma;
      m.loc1 = state.loc1() + static_cast<Scalar>(gauss(rng)) * sigma;
      m.varLoc0 = sigma * sigma;
      m.varLoc1 = sigma * sigma;
      m.clusterIndex = static_cast<std::uint32_t>(layer);
      if (radiate && layer + 1 < kLayers) {
        state.par[kQOverP] /= fraction(halfGap(layer + 1), rng);
      }
      if (!isFiniteScalar(state.qOverP()) || std::abs(state.qOverP()) > 50.F) {
        return false;
      }
    }
    return true;
  }
};

/// The q/p marginal's distribution function at \a u: the probability
/// integral transform of the truth.
template <int K>
double marginalCdf(const GSFMixture<K>& mixture, double u)
{
  double sum = 0.;
  for (int k = 0; k < mixture.size; ++k) {
    const double m = mixture.state[k].qOverP();
    const double s = std::sqrt(static_cast<double>(mixture.state[k].cov[covIndex(kQOverP, kQOverP)]));
    sum += std::exp(static_cast<double>(mixture.logWeight[k])) * 0.5 * std::erfc(-(u - m) / (s * std::numbers::sqrt2));
  }
  return sum;
}

/// The median and the robust width (half the 16-84 % range) of \a values.
std::pair<double, double> medianWidth(std::vector<double> values)
{
  std::ranges::sort(values);
  const auto at = [&](double q) { return values[static_cast<std::size_t>(q * static_cast<double>(values.size() - 1))]; };
  return {at(0.5), 0.5 * (at(0.84) - at(0.16))};
}

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

template <int K, typename Loss = BetheHeitlerLoss>
Coverage coverage(double pGeV, const UniformMaterial& material, bool radiate, int tracks, std::uint32_t seed)
{
  ElectronChain chain;
  std::mt19937 rng(seed);
  Coverage out;
  for (int track = 0; track < tracks; ++track) {
    const Scalar qOverP = static_cast<Scalar>((track % 2 == 0 ? 1. : -1.) / pGeV);
    while (!chain.generate(qOverP, material, radiate, rng)) {
    }
    const auto begin = std::chrono::steady_clock::now();
    const auto gsf = fitGaussianSum<K, Loss>(chain.propagator, chain.tracking, material, chain.measurements.data(), static_cast<int>(ElectronChain::kLayers));
    out.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    const auto kalmanBegin = std::chrono::steady_clock::now();
    const auto kalman = fitMeasurements(chain.propagator, chain.tracking, material, chain.measurements.data(), static_cast<int>(ElectronChain::kLayers));
    out.kalmanSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - kalmanBegin).count();
    if (kalman.status != PropagationStatus::kOk) {
      ++out.kalmanFailed; // a track the fit's propagation cannot follow (a curler's turning point)
      continue;
    }
    if (gsf.inward.status != PropagationStatus::kOk) {
      ++out.failed;
      continue;
    }
    ++out.fits;
    const double u = chain.truth[0].qOverP();
    const double p = 1. / std::abs(u);
    out.pit.push_back(marginalCdf(gsf.inward.mixture, u));
    out.outwardPit.push_back(marginalCdf(gsf.outward.mixture, chain.truth[ElectronChain::kLayers - 1].qOverP()));
    out.kalmanRatio.push_back(1. / std::abs(static_cast<double>(kalman.state.qOverP())) / p);
    out.modeRatio.push_back(static_cast<double>(gsf.inward.momentumMode) / p);
    out.meanRatio.push_back(1. / std::abs(static_cast<double>(gsf.inward.collapsed.qOverP())) / p);
    out.collapsedPull.push_back((static_cast<double>(gsf.inward.collapsed.qOverP()) - u) / std::sqrt(static_cast<double>(gsf.inward.collapsed.variance(kQOverP))));
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
    const auto gsf = fitGaussianSum<4, NoLoss>(chain.propagator, chain.tracking, material, chain.measurements.data(), 4);
    BOOST_REQUIRE(kalman.status == PropagationStatus::kOk);
    BOOST_REQUIRE(gsf.outward.status == PropagationStatus::kOk);
    BOOST_REQUIRE(gsf.inward.status == PropagationStatus::kOk);
    BOOST_TEST(gsf.outward.mixture.size == 1);
    BOOST_TEST(gsf.inward.mixture.size == 1);
    BOOST_TEST(sameBits(gsf.outward.collapsed.par, kalman.outwardState.par, kNBoundParams));
    BOOST_TEST(sameBits(gsf.outward.collapsed.cov, kalman.outwardState.cov, kBoundCovSize));
    BOOST_TEST(sameBits(gsf.inward.collapsed.par, kalman.state.par, kNBoundParams));
    BOOST_TEST(sameBits(gsf.inward.collapsed.cov, kalman.state.cov, kBoundCovSize));
    BOOST_TEST(gsf.outward.mixture.chi2[0] == kalman.outwardChi2);
    BOOST_TEST(gsf.inward.mixture.chi2[0] == kalman.chi2);
    BOOST_TEST(gsf.inward.qOverPMode == kalman.state.qOverP());
    ++fits;
  }
  BOOST_TEST(fits == 50);
}

/// Every loss path enumerated: each plane the filter branches at -- the two
/// the start leaves, the one after it and the last -- takes one chosen loss
/// component, each path is filtered alone with the same steps (the diffuse
/// start's, then the filter's), and its weight is the product of its loss
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
    const auto gsf = fitGaussianSum<kPaths, TwoLosses>(chain.propagator, chain.tracking, material, measurements, kLayers);
    BOOST_REQUIRE(gsf.outward.status == PropagationStatus::kOk);
    BOOST_REQUIRE(gsf.inward.status == PropagationStatus::kOk);
    const auto seed = makeChainSeed(chain.propagator, chain.tracking, measurements, kLayers);
    const auto span = detail::spanMaterial(chain.tracking, material, measurements, kLayers);
    for (const auto direction : {PropagationDirection::kForward, PropagationDirection::kBackward}) {
      const bool outward = direction == PropagationDirection::kForward;
      const auto& pass = outward ? gsf.outward : gsf.inward;
      BOOST_REQUIRE(pass.mixture.size == kPaths);
      const int first = outward ? 0 : kLayers - 1;
      const int increment = outward ? 1 : -1;
      double logWeight[kPaths];
      BoundTrackState state[kPaths];
      CovarianceFactors factors[kPaths];
      for (int path = 0; path < kPaths; ++path) {
        const auto lossOf = [&](int b) { return TwoLosses::kLosses[(path >> (kBranchings - 1 - b)) & 1]; };
        state[path] = outward ? seed : gsf.outward.collapsed;
        logWeight[path] = 0.;
        Scalar chi2 = 0.F;
        detail::DiffuseFactors diffuse;
        BOOST_REQUIRE(detail::startDiffuse(measurements[first], state[path], diffuse));
        int plane = first;
        for (int b = 0; b < detail::kPassStart - 1; ++b, plane += increment) {
          const auto loss = lossOf(b);
          logWeight[path] += std::log(static_cast<double>(loss.weight));
          BOOST_REQUIRE(detail::gsf::applyLossDiffuse(loss, direction, state[path], diffuse));
          Scalar density = 0.F;
          BOOST_REQUIRE((detail::gsf::diffuseBranchStep<kRule>(chain.propagator, chain.tracking, span, measurements[plane + increment], direction, state[path], diffuse, chi2, density) == PropagationStatus::kOk));
          logWeight[path] -= 0.5 * static_cast<double>(density);
        }
        BOOST_REQUIRE(diffuse.rank == 0);
        factors[path] = diffuse.star;
        for (int b = detail::kPassStart - 1; b < kBranchings; ++b, plane += increment) {
          const auto loss = lossOf(b);
          logWeight[path] += std::log(static_cast<double>(loss.weight));
          BOOST_REQUIRE(applyBetheHeitler(loss, direction, state[path], factors[path]));
          if (b + 1 < kBranchings) {
            Scalar logLikelihood = 0.F;
            BOOST_REQUIRE((detail::gsf::componentStep<kRule>(chain.propagator, chain.tracking, span, measurements[plane + increment], direction, state[path], factors[path], chi2, logLikelihood) == PropagationStatus::kOk));
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
        const double got = std::exp(static_cast<double>(pass.mixture.logWeight[path]));
        worstWeight = std::max(worstWeight, std::abs(got - expected));
        for (int i = 0; i < kNBoundParams; ++i) {
          const double sigma = std::sqrt(static_cast<double>(state[path].variance(i)));
          worstState = std::max(worstState, std::abs(static_cast<double>(pass.mixture.state[path].par[i] - state[path].par[i])) / sigma);
        }
        for (int i = 0; i < kBoundCovSize; ++i) {
          const double scale = std::abs(static_cast<double>(state[path].cov[i])) + 1e-30;
          worstState = std::max(worstState, std::abs(static_cast<double>(pass.mixture.state[path].cov[i] - state[path].cov[i])) / scale);
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
    rows.push_back({fmt::format("thin, p {:.1f} GeV, K 12", p), coverage<12>(p, thin, true, kTracks, 105), true});
  }
  for (const double p : {0.2, 0.5}) {
    rows.push_back({fmt::format("thick, p {:.1f} GeV, K 12", p), coverage<12>(p, thick, true, kTracks, 106), true});
  }
  rows.push_back({"thick, p 0.5 GeV, K 1", coverage<1>(0.5, thick, true, kTracks, 106), false});
  rows.push_back({"thick, p 0.5 GeV, K 4", coverage<4>(0.5, thick, true, kTracks, 106), false});
  rows.push_back({"thick, p 0.5 GeV, K 24", coverage<24>(0.5, thick, true, kTracks, 106), false});
  // Who owns which loss: the ionisation alone (the no-loss model on tracks without radiation, the Kalman fit's
  // own q/p coverage), the radiation alone (no ionisation), and a radiation the tracks do not have.
  rows.push_back({"thick, p 0.5 GeV, ionisation only", coverage<12, NoLoss>(0.5, thick, false, kTracks, 107), false});
  const UniformMaterial thickNoIonisation{.x2X0PerCm = thick.x2X0PerCm, .rho = 0.F};
  rows.push_back({"thick, p 0.5 GeV, radiation only", coverage<12>(0.5, thickNoIonisation, true, kTracks, 108), true});
  rows.push_back({"thick, p 0.5 GeV, tracks without radiation", coverage<12>(0.5, thick, false, kTracks, 107), false});
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

} // namespace o2::gtrk
