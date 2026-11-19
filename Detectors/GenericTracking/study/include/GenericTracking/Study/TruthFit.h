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

/// \file TruthFit.h
/// \brief A chain's rows (StudyChains.h) fitted with its particle's true
///        species: the mass for the scattering and for the energy loss
///        together and apart, and the PID's own refit (HypothesisRefit.h).
///
/// The fit takes one mass hypothesis for both (`Propagator::massGeV`). An
/// `OraclePropagator` carries the energy loss's own (`lossMassGeV`);
/// `stepPropagator` hands it on whole, so the fit's steps see it. The
/// tracking never instantiates it.

#ifndef GTRK_STUDY_TRUTHFIT_H_
#define GTRK_STUDY_TRUTHFIT_H_

#include "GenericTracking/Tracking/KFTrackFitter.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/RoadFit.h"

#include "CommonConstants/PhysicsConstants.h"
#include "ReconstructionDataFormats/PID.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <optional>
#include <stdexcept>

namespace o2::gtrk::study
{

/// A propagator whose energy loss takes \a lossMass, its scattering
/// `massGeV`.
struct OraclePropagator : Propagator {
  Scalar lossMass = static_cast<Scalar>(o2::constants::physics::MassPionCharged);
};

[[nodiscard]] inline Scalar lossMassGeV(const OraclePropagator& prop) { return prop.lossMass; }

/// The oracle's propagator for every step: itself, as `stepPropagator`
/// gives a `Propagator`.
template <typename Geometry, typename Chain>
[[nodiscard]] const OraclePropagator& stepPropagator(const OraclePropagator& prop, const Geometry& /*geometry*/, const Chain& /*measurements*/, int /*from*/, int /*to*/)
{
  return prop; // NOLINT(bugprone-return-const-ref-from-parameter)
}

/// A propagator whose fit steps take Bethe-Heitler's bremsstrahlung as one
/// Gaussian, for every species alike: over t = x/X0 the retained energy
/// fraction z has the mean mu = e^-t and the variance
/// v = e^(-t ln3 / ln2) - mu^2; outward q/p becomes q/p / z, inward q/p z.
/// With `mean` the state takes mu (outward q/p / mu, inward q/p mu); with
/// `variance` the step's q/p variance grows by the linearised v (outward
/// (q/p)^2 v / mu^4, inward (q/p)^2 v).
struct RadiationPropagator : Propagator {
  bool mean = true;
  bool variance = true;
};

inline void applyRadiation(const RadiationPropagator& prop, Scalar x2X0, PropagationDirection direction, BoundParameters& state, Scalar& varQOverP)
{
  constexpr Scalar kLn3OverLn2 = 1.5849625F;
  const Scalar mu = std::exp(-x2X0);
  const Scalar v = std::exp(-x2X0 * kLn3OverLn2) - (mu * mu);
  const Scalar u = state.par[kQOverP];
  const bool outward = direction == PropagationDirection::kForward;
  if (prop.variance) {
    varQOverP += outward ? u * u * v / (mu * mu * mu * mu) : u * u * v;
  }
  if (prop.mean) {
    state.par[kQOverP] = outward ? u / mu : u * mu;
  }
}

/// The radiation propagator for every step: itself, as `stepPropagator`
/// gives a `Propagator`.
template <typename Geometry, typename Chain>
[[nodiscard]] const RadiationPropagator& stepPropagator(const RadiationPropagator& prop, const Geometry& /*geometry*/, const Chain& /*measurements*/, int /*from*/, int /*to*/)
{
  return prop; // NOLINT(bugprone-return-const-ref-from-parameter)
}

/// The fits of a chain, each without kinks: the pion's two passes (the
/// control); the true mass's, for both, for the scattering alone, for the
/// energy loss alone; the true mass's inward pass alone from the chain's
/// outward state (a track's own, as the tracking's refit; a truth chain's
/// pion fit's); the electron hypothesis for every species, the inward pass
/// from the same outward state. Then the true PID's refit
/// (`refitWithHypothesis`: the Gaussian-sum filter for an electron, the
/// Kalman fit otherwise), carried to the perigee about the production point;
/// `truthFit` makes the first six.
enum class TruthFitVariant : std::uint8_t { kPion,
                                            kTrue,
                                            kTrueScattering,
                                            kTrueLoss,
                                            kTrueInward,
                                            kElectronInward,
                                            kHypothesis,
                                            kCount };
inline constexpr std::size_t kTruthFitVariants = static_cast<std::size_t>(TruthFitVariant::kCount);

/// A species' mass, as the fits take it, and its PID.
struct Species {
  Scalar mass = 0.F; ///< [GeV]
  o2::track::PID pid;
};

/// The species of \a pdg: e, mu, pi, K, p; none for another.
[[nodiscard]] inline std::optional<Species> truthSpecies(std::int32_t pdg)
{
  namespace physics = o2::constants::physics;
  using o2::track::PID;
  switch (std::abs(pdg)) {
    case 11:
      return Species{.mass = static_cast<Scalar>(physics::MassElectron), .pid = PID{PID::Electron}};
    case 13:
      return Species{.mass = static_cast<Scalar>(physics::MassMuon), .pid = PID{PID::Muon}};
    case 211:
      return Species{.mass = static_cast<Scalar>(physics::MassPionCharged), .pid = PID{PID::Pion}};
    case 321:
      return Species{.mass = static_cast<Scalar>(physics::MassKaonCharged), .pid = PID{PID::Kaon}};
    case 2212:
      return Species{.mass = static_cast<Scalar>(physics::MassProton), .pid = PID{PID::Proton}};
    default:
      return std::nullopt;
  }
}

/// \a variant's fit (one of the first six) of the \a n rows of \a chain, with
/// \a base's field and tolerances and \a mass, the particle's; \a outward the
/// outward state the inward-only variants start from.
template <DetectorLayout Det, typename Material, typename Chain>
[[nodiscard]] TrackFitResult truthFit(TruthFitVariant variant, Scalar mass, const Propagator& base, const TrackingGeometry<Det>& geometry, const Material& material, const Chain& chain, int n, const BoundTrackState& outward)
{
  const auto pion = static_cast<Scalar>(o2::constants::physics::MassPionCharged);
  const auto electron = static_cast<Scalar>(o2::constants::physics::MassElectron);
  const auto with = [&](Scalar scattering, Scalar loss) {
    OraclePropagator prop;
    static_cast<Propagator&>(prop) = base;
    prop.massGeV = scattering;
    prop.lossMass = loss;
    return prop;
  };
  switch (variant) {
    case TruthFitVariant::kPion:
      return fitMeasurements(with(pion, pion), geometry, material, chain, n);
    case TruthFitVariant::kTrue:
      return fitMeasurements(with(mass, mass), geometry, material, chain, n);
    case TruthFitVariant::kTrueScattering:
      return fitMeasurements(with(mass, pion), geometry, material, chain, n);
    case TruthFitVariant::kTrueLoss:
      return fitMeasurements(with(pion, mass), geometry, material, chain, n);
    case TruthFitVariant::kTrueInward:
      return fitInward(with(mass, mass), geometry, material, chain, n, outward);
    case TruthFitVariant::kElectronInward:
      return fitInward(with(electron, electron), geometry, material, chain, n, outward);
    case TruthFitVariant::kHypothesis:
    case TruthFitVariant::kCount:
      break;
  }
  throw std::invalid_argument("truthFit: kHypothesis is refitWithHypothesis's, kCount no variant");
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_TRUTHFIT_H_
