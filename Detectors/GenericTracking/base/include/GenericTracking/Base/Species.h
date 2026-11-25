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

/// \file Species.h
/// \brief The species a track can be fitted as: a mass and an |q|.
///
/// The charged particles that leave rows of their own in the tracker: the
/// stable ones and the weakly decaying ones that live long enough (Xi,
/// Omega, Sigma, the hypernuclei). One entry covers a particle and its
/// antiparticle: the curvature measures the sign, a hypothesis is the mass
/// and |q|. `o2::track::PID` lacks the lithium ions and the Sigmas, so the
/// tracking keeps its own list; the list is fixed and its order is policy.

#ifndef GTRK_BASE_SPECIES_H_
#define GTRK_BASE_SPECIES_H_

#include "GenericTracking/Base/Definitions.h"

#include "CommonConstants/PhysicsConstants.h"

#include <cstdint>

#ifndef GPUCA_GPUCODE_DEVICE
#include <optional>
#include <string_view>
#endif

namespace o2::gtrk
{

enum class Species : std::uint8_t { kElectron,
                                    kMuon,
                                    kPion,
                                    kKaon,
                                    kProton,
                                    kDeuteron,
                                    kTriton,
                                    kHelium3,
                                    kHelium4,
                                    kLithium6,
                                    kLithium7,
                                    kXi,
                                    kOmega,
                                    kSigmaPlus,
                                    kSigmaMinus,
                                    kHyperTriton,
                                    kHyperHydrogen4,
                                    kHyperHelium4,
                                    kHyperHelium5,
                                    kCount };
inline constexpr int kNSpecies = static_cast<int>(Species::kCount);

/// The mass [GeV] of \a species: O2's (`PhysicsConstants.h`); the lithium
/// ions' from O2's ion table (`extra_ions_pdg_table.dat`).
GPUhdi() constexpr float massOf(Species species)
{
  namespace physics = o2::constants::physics;
  switch (species) {
    case Species::kElectron:
      return static_cast<float>(physics::MassElectron);
    case Species::kMuon:
      return static_cast<float>(physics::MassMuon);
    case Species::kPion:
      return static_cast<float>(physics::MassPionCharged);
    case Species::kKaon:
      return static_cast<float>(physics::MassKaonCharged);
    case Species::kProton:
      return static_cast<float>(physics::MassProton);
    case Species::kDeuteron:
      return static_cast<float>(physics::MassDeuteron);
    case Species::kTriton:
      return static_cast<float>(physics::MassTriton);
    case Species::kHelium3:
      return static_cast<float>(physics::MassHelium3);
    case Species::kHelium4:
      return static_cast<float>(physics::MassAlpha);
    case Species::kLithium6:
      return 5.60152F;
    case Species::kLithium7:
      return 6.53383F;
    case Species::kXi:
      return static_cast<float>(physics::MassXiMinus);
    case Species::kOmega:
      return static_cast<float>(physics::MassOmegaMinus);
    case Species::kSigmaPlus:
      return static_cast<float>(physics::MassSigmaPlus);
    case Species::kSigmaMinus:
      return static_cast<float>(physics::MassSigmaMinus);
    case Species::kHyperTriton:
      return static_cast<float>(physics::MassHyperTriton);
    case Species::kHyperHydrogen4:
      return static_cast<float>(physics::MassHyperHydrogen4);
    case Species::kHyperHelium4:
      return static_cast<float>(physics::MassHyperHelium4);
    case Species::kHyperHelium5:
      return static_cast<float>(physics::MassHyperHelium5);
    case Species::kCount:
      break;
  }
  return 0.F;
}

/// The |q| [e] of \a species.
GPUhdi() constexpr float absChargeOf(Species species)
{
  switch (species) {
    case Species::kHelium3:
    case Species::kHelium4:
    case Species::kHyperHelium4:
    case Species::kHyperHelium5:
      return 2.F;
    case Species::kLithium6:
    case Species::kLithium7:
      return 3.F;
    default:
      return 1.F;
  }
}

#ifndef GPUCA_GPUCODE_DEVICE
/// The PDG code of \a species' particle (Sigma+ and Sigma- each their own).
[[nodiscard]] constexpr std::int32_t pdgOf(Species species)
{
  switch (species) {
    case Species::kElectron:
      return 11;
    case Species::kMuon:
      return 13;
    case Species::kPion:
      return 211;
    case Species::kKaon:
      return 321;
    case Species::kProton:
      return 2212;
    case Species::kDeuteron:
      return 1000010020;
    case Species::kTriton:
      return 1000010030;
    case Species::kHelium3:
      return 1000020030;
    case Species::kHelium4:
      return 1000020040;
    case Species::kLithium6:
      return 1000030060;
    case Species::kLithium7:
      return 1000030070;
    case Species::kXi:
      return 3312;
    case Species::kOmega:
      return 3334;
    case Species::kSigmaPlus:
      return 3222;
    case Species::kSigmaMinus:
      return 3112;
    case Species::kHyperTriton:
      return 1010010030;
    case Species::kHyperHydrogen4:
      return 1010010040;
    case Species::kHyperHelium4:
      return 1010020040;
    case Species::kHyperHelium5:
      return 1010020050;
    case Species::kCount:
      break;
  }
  return 0;
}

/// The species of PDG code \a pdg, either sign; none for a particle not in
/// the list.
[[nodiscard]] constexpr std::optional<Species> speciesOf(std::int32_t pdg)
{
  const std::int32_t code = pdg < 0 ? -pdg : pdg;
  for (int s = 0; s < kNSpecies; ++s) {
    if (pdgOf(static_cast<Species>(s)) == code) {
      return static_cast<Species>(s);
    }
  }
  return std::nullopt;
}

/// A short name of \a species, for reports and record columns.
[[nodiscard]] constexpr std::string_view nameOf(Species species)
{
  constexpr std::string_view kNames[kNSpecies] = {"e", "mu", "pi", "K", "p", "d", "t", "3He", "4He", "6Li", "7Li", "Xi", "Omega", "Sigma+", "Sigma-", "3LH", "4LH", "4LHe", "5LHe"};
  return species == Species::kCount ? std::string_view{} : kNames[static_cast<int>(species)];
}
#endif

} // namespace o2::gtrk

#endif // GTRK_BASE_SPECIES_H_
