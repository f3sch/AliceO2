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

/// \file StateTruth.h
/// \brief A fitted state against its particle's MC crossing of the state's
///        sensor, and a particle's charge from its PDG code.

#ifndef GTRK_STUDY_STATETRUTH_H_
#define GTRK_STUDY_STATETRUTH_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Study/Alice/Hits.h"
#include "GenericTracking/Study/TrackComparison.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include "SimulationDataFormat/O2DatabasePDG.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <optional>

namespace o2::gtrk::study
{

/// The charge [e] of the species \a pdg, 0 for one O2's PDG table does not
/// have.
[[nodiscard]] inline double particleCharge(std::int32_t pdg)
{
  const auto* particle = o2::O2DatabasePDG::Instance()->GetParticle(pdg);
  return particle == nullptr ? 0. : particle->Charge() / 3.;
}

/// The mass [GeV] of the species \a pdg, 0 for one O2's PDG table does not
/// have.
[[nodiscard]] inline double particleMass(std::int32_t pdg)
{
  const auto* particle = o2::O2DatabasePDG::Instance()->GetParticle(pdg);
  return particle == nullptr ? 0. : particle->Mass();
}

/// A fitted state against its particle's MC crossing.
struct StateTruth {
  std::array<double, 5> residual{}; ///< fitted less true: loc0, loc1, phi, theta, q/p
  std::array<double, 5> sigma{};    ///< the fit's
  double pt = 0.;                   ///< the truth's [GeV]
  double momentum = 0.;             ///< the truth's [GeV]
  double fittedPt = 0.;             ///< [GeV]
  double ptSigma = 0.;              ///< of fittedPt, J C J^T over theta and q/p [GeV]
};

/// Where the straight MC segment from \a start to \a end crosses \a plane:
/// its loc0, loc1 there; none for a segment along the plane.
[[nodiscard]] inline std::optional<std::array<double, 2>> segmentCrossing(const TrackingSurface& plane, const std::array<double, 3>& start, const std::array<double, 3>& end)
{
  const auto local = [&](const std::array<double, 3>& global) {
    std::array<double, 3> l{};
    for (std::size_t i = 0; i < 3U; ++i) {
      const double d = global[i] - static_cast<double>(plane.origin[i]);
      l[0] += d * static_cast<double>(plane.axis0[i]);
      l[1] += d * static_cast<double>(plane.axis1[i]);
      l[2] += d * static_cast<double>(plane.normal[i]);
    }
    return l;
  };
  const auto a = local(start);
  const auto b = local(end);
  if (std::abs(a[2] - b[2]) < 1e-9) {
    return std::nullopt;
  }
  const double t = a[2] / (a[2] - b[2]);
  return std::array<double, 2>{a[0] + (t * (b[0] - a[0])), a[1] + (t * (b[1] - a[1]))};
}

/// \a state of particle \a particle (of charge \a charge) against the
/// particle's MC crossing of the state's sensor, the one nearest the
/// measurement at \a row, on the plane the state is bound to
/// (TrackingGeometry::planeOf: a flat sensor's own at its measurement depth, a
/// curved sensor's tangent plane at the state's angle): where the MC
/// segment's normal coordinate is zero; its direction and momentum there.
/// None without the segment, for a neutral particle or a segment along the
/// plane.
template <DetectorLayout Det>
[[nodiscard]] std::optional<StateTruth> stateTruth(const TrackingGeometry<Det>& geometry, const HitIndex& hits, const TimeFrame<Det>& frame, const BoundTrackState& state, std::uint64_t particle, std::uint32_t row, double charge)
{
  const auto m = frame.measurement(row);
  const auto* segment = hits.crossing(particle, state.surface, {static_cast<double>(m.x), static_cast<double>(m.y), static_cast<double>(m.z)});
  if (segment == nullptr || exactlyEqual(charge, 0.)) {
    return std::nullopt;
  }
  const TrackingSurface plane = geometry.planeOf(state); // a copy: a curved sensor's is made on the spot
  const auto crossing = segmentCrossing(plane, segment->start, segment->end);
  if (!crossing) {
    return std::nullopt;
  }
  const auto& p = segment->momentum;
  StateTruth truth;
  truth.pt = std::hypot(p[0], p[1]);
  const double momentum = std::hypot(truth.pt, p[2]);
  truth.momentum = momentum;
  truth.residual = {static_cast<double>(state.par[kLoc0]) - (*crossing)[0], static_cast<double>(state.par[kLoc1]) - (*crossing)[1], std::remainder(static_cast<double>(state.par[kPhi]) - std::atan2(p[1], p[0]), 2. * std::numbers::pi), static_cast<double>(state.par[kTheta]) - std::atan2(truth.pt, p[2]), static_cast<double>(state.par[kQOverP]) - (charge / momentum)};
  for (int i = 0; i < 5; ++i) {
    truth.sigma[static_cast<std::size_t>(i)] = std::sqrt(std::max(static_cast<double>(state.covAt(i, i)), 0.));
  }
  // pT = sin(theta) / |q/p|: d/dtheta = cos(theta) / |q/p|, d/d(q/p) = -pT / (q/p).
  const auto theta = static_cast<double>(state.par[kTheta]);
  const auto qOverP = static_cast<double>(state.par[kQOverP]);
  truth.fittedPt = std::sin(theta) / std::abs(qOverP);
  const double dTheta = std::cos(theta) / std::abs(qOverP);
  const double dQOverP = -truth.fittedPt / qOverP;
  const double variance = (dTheta * dTheta * static_cast<double>(state.covAt(kTheta, kTheta))) + (2. * dTheta * dQOverP * static_cast<double>(state.covAt(kTheta, kQOverP))) + (dQOverP * dQOverP * static_cast<double>(state.covAt(kQOverP, kQOverP)));
  truth.ptSigma = std::sqrt(std::max(variance, 0.));
  return truth;
}

/// A state on the perigee plane about a particle's MC production point
/// against its MC kinematics there.
struct PerigeeTruth {
  std::array<double, 5> residual{}; ///< fitted less true: d0, z, phi, theta, q/p
  std::array<double, 5> sigma{};    ///< the fit's
};

/// \a state on the perigee plane about the production point of the particle
/// of \a truth and \a charge (`perigeePlaneOf`, BeamConstraint.h: loc0 the
/// signed transverse distance, the d0, loc1 the z) against its production
/// point and momentum.
[[nodiscard]] inline PerigeeTruth perigeeTruth(const BoundTrackState& state, const ParticleTruth& truth, double charge)
{
  const auto eta = static_cast<double>(truth.eta);
  const double momentum = static_cast<double>(truth.pt) * std::cosh(eta);
  PerigeeTruth result;
  result.residual = {static_cast<double>(state.par[kLoc0]), static_cast<double>(state.par[kLoc1]) - static_cast<double>(truth.vz), std::remainder(static_cast<double>(state.par[kPhi]) - static_cast<double>(truth.phi), 2. * std::numbers::pi), static_cast<double>(state.par[kTheta]) - (2. * std::atan(std::exp(-eta))), static_cast<double>(state.par[kQOverP]) - (charge / momentum)};
  for (int i = 0; i < 5; ++i) {
    result.sigma[static_cast<std::size_t>(i)] = std::sqrt(std::max(static_cast<double>(state.covAt(i, i)), 0.));
  }
  return result;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_STATETRUTH_H_
