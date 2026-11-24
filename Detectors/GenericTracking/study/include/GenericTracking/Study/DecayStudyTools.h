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

/// \file DecayStudyTools.h
/// \brief What the decay studies on the tracker's tracks share: the V0
///        hypotheses, the daughter table filled on the event's threads, an
///        entry's displacement against a reference point, a candidate's fit
///        into a record, the fitted vertex against the true decay point, the
///        Armenteros-Podolanski variables, a mother's rapidities.

#ifndef GTRK_STUDY_DECAYSTUDYTOOLS_H_
#define GTRK_STUDY_DECAYSTUDYTOOLS_H_

#include "GenericTracking/Tracking/DaughterTable.h"
#include "GenericTracking/Tracking/DecayCandidateFit.h"
#include "GenericTracking/Tracking/HypothesisRefit.h"
#include "GenericTracking/Study/StateTruth.h"
#include "GenericTracking/Study/TruthEvent.h"

#include "ReconstructionDataFormats/PID.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <span>
#include <stdexcept>
#include <string_view>

namespace o2::gtrk::study
{

/// The hypothesis of a decay into \a positive and \a negative, the mother's
/// mass \a massGeV.
[[nodiscard]] inline DecayHypothesis twoProngHypothesis(o2::track::PID positive, o2::track::PID negative, double massGeV)
{
  DecayHypothesis h;
  h.species.push_back(positive);
  h.charges.push_back(1);
  h.species.push_back(negative);
  h.charges.push_back(-1);
  h.massGeV = massGeV;
  return h;
}

/// The V0s' species: K0S, Lambda, anti-Lambda (`v0Hypotheses`).
inline constexpr std::array<std::int32_t, 3> kV0Mothers{310, 3122, -3122};

/// The V0 hypotheses in `kV0Mothers`' order: K0S -> pi+ pi-, Lambda -> p pi-,
/// anti-Lambda -> anti-p pi+.
[[nodiscard]] inline std::array<DecayHypothesis, 3> v0Hypotheses()
{
  using o2::track::PID;
  return {twoProngHypothesis(PID{PID::Pion}, PID{PID::Pion}, particleMass(kV0Mothers[0])), twoProngHypothesis(PID{PID::Proton}, PID{PID::Pion}, particleMass(kV0Mothers[1])), twoProngHypothesis(PID{PID::Pion}, PID{PID::Proton}, particleMass(kV0Mothers[2]))};
}

/// The daughter table's forEach on \a event's threads: body(i) for every i
/// below n.
template <DetectorLayout Det>
[[nodiscard]] auto eventForEach(const TruthEvent<Det>& event)
{
  return [&event](std::uint32_t n, const auto& body) {
    o2::gtrk::detail::parallelFor(*event.arena, 0U, n, 16, [&](std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; ++i) {
        body(static_cast<std::uint32_t>(i));
      }
    });
  };
}

/// \a table begun for \a event's tracks, its pairs asked for by
/// \a request(table) and refitted on the event's threads;
/// std::runtime_error naming \a study where it cannot be.
template <DetectorLayout Det, typename Request>
void fillDaughterTable(DaughterTable<Det>& table, const TruthEvent<Det>& event, std::string_view study, const Request& request)
{
  if (event.tracks == nullptr || table.begin(event.tracks->size()) != StageStatus::kComplete) {
    throw std::runtime_error(std::format("{}: no daughter table", study));
  }
  request(table);
  if (table.flush(eventForEach(event), 4096U, event.propagator, *event.geometry, *event.frame, *event.tracks) != StageStatus::kComplete) {
    throw std::runtime_error(std::format("{}: the daughter table's flush failed", study));
  }
}

/// \a table's entries of \a event's tracks as the candidates' daughters.
template <DetectorLayout Det>
[[nodiscard]] DaughterSource<Det> daughterSource(const DaughterTable<Det>& table, const TruthEvent<Det>& event)
{
  return {.table = &table, .tracks = event.tracks, .frame = event.frame, .geometry = event.geometry, .propagator = event.propagator};
}

/// Entry \a e's displacement against \a point: its impact parameters'
/// significance at its perigee about the point's transverse position,
/// S^2 = r^T (V + J C J^T)^-1 r -- r the transverse and longitudinal impact
/// parameters, V their covariance, C the point's, J their derivative by the
/// point (-u across, cot(theta) t - z along); NaN if the carry fails or the
/// matrix is singular.
template <DetectorLayout Det>
[[nodiscard]] double displacement(const DaughterSource<Det>& source, std::uint32_t e, const MeasuredPoint& point)
{
  const auto entries = source.table->view();
  HypothesisFit<> fit{.pid = source.table->species()[entries.species[e]]};
  fit.innermost = entries.state[e];
  if (carryToPerigee(source.propagator, *source.geometry, source.innermostPlane(entries.track[e]), fit, static_cast<Scalar>(point.position[0]), static_cast<Scalar>(point.position[1])) != PropagationStatus::kOk) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const auto& perigee = fit.perigee;
  const std::array<double, 2> r{static_cast<double>(perigee.loc0()), static_cast<double>(perigee.loc1()) - point.position[2]};
  // The loc0 axis in global coordinates, from the plane itself.
  std::array<Scalar, 3> origin{};
  std::array<Scalar, 3> across{};
  surfaceToGlobal(fit.perigeePlane, 0.F, 0.F, origin.data());
  surfaceToGlobal(fit.perigeePlane, 1.F, 0.F, across.data());
  const auto phi = static_cast<double>(perigee.phi());
  const double cotTheta = 1. / std::tan(static_cast<double>(perigee.theta()));
  const std::array<std::array<double, 3>, 2> j{{{-(across[0] - origin[0]), -(across[1] - origin[1]), -(across[2] - origin[2])}, {cotTheta * std::cos(phi), cotTheta * std::sin(phi), -1.}}};
  std::array<double, 3> m{static_cast<double>(perigee.covAt(kLoc0, kLoc0)), static_cast<double>(perigee.covAt(kLoc0, kLoc1)), static_cast<double>(perigee.covAt(kLoc1, kLoc1))};
  const auto projected = [&](int a, int b) {
    double sum = 0.;
    for (int k = 0; k < 3; ++k) {
      for (int l = 0; l < 3; ++l) {
        sum += j[a][k] * point.covariance[covIndex(k, l)] * j[b][l];
      }
    }
    return sum;
  };
  m[0] += projected(0, 0);
  m[1] += projected(0, 1);
  m[2] += projected(1, 1);
  const double det = (m[0] * m[2]) - (m[1] * m[1]);
  if (!(det > 0.) || !(m[0] > 0.)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return ((r[0] * r[0] * m[2]) - (2. * r[0] * r[1] * m[1]) + (r[1] * r[1] * m[0])) / det;
}

/// The leaves of \a candidate, fitted into \a result (converged), whose
/// first measurement lies behind their vertex in the fit (`planeAhead`: the
/// fit carried them forward to it); a decay's own lie before it.
template <DetectorLayout Det>
[[nodiscard]] int leavesAhead(const DecayCandidate& candidate, std::span<const DecayHypothesis> hypotheses, const DaughterSource<Det>& source, const DecayCandidateResult& result)
{
  const auto tree = expandHypothesis(hypotheses, candidate.hypothesis);
  const auto daughters = candidateDaughters(candidate, source);
  int ahead = 0;
  for (std::size_t i = 0; i < daughters.size(); ++i) {
    const auto& vertex = result.fit.vertices[static_cast<std::size_t>(tree.topology.vertexOfDaughter(static_cast<int>(i)))];
    const auto& daughter = daughters[i];
    const auto plane = detail::hypothesis::perigeePlaneOf(daughter.propagator, daughter.innermost, daughter.surface, static_cast<Scalar>(vertex[0]), static_cast<Scalar>(vertex[1]));
    ahead += detail::hypothesis::planeAhead(daughter.innermost, daughter.surface, plane) ? 1 : 0;
  }
  return ahead;
}

/// The frame of a true flight from \a production to \a decay: along it,
/// across it in the transverse plane, and the third axis (along x across).
struct FlightFrame {
  Vec3D along{};
  Vec3D across{};
  Vec3D third{};
};

[[nodiscard]] inline FlightFrame flightFrame(const Vec3D& production, const Vec3D& decay)
{
  FlightFrame frame{.along = {decay[0] - production[0], decay[1] - production[1], decay[2] - production[2]}};
  const double flight = std::hypot(frame.along[0], frame.along[1], frame.along[2]);
  for (auto& component : frame.along) {
    component /= flight;
  }
  const double transverse = std::hypot(frame.along[0], frame.along[1]);
  frame.across = transverse > 0. ? Vec3D{-frame.along[1] / transverse, frame.along[0] / transverse, 0.} : Vec3D{1., 0., 0.};
  frame.third = crossOf(frame.along, frame.across);
  return frame;
}

/// \a fit's root vertex less the true \a decay point of a mother made at
/// \a production into \a record: the z difference (dZ), the transverse
/// distance (dXY) if the record has it, and, if the record has them, the residuals and the
/// fit's sigmas along and across the true flight (`FlightFrame`: dL, sL, dT,
/// sT, dN, sN) for the pulls.
template <typename Record>
void writeVertexResiduals(Record& record, const DecayVertexFitResult& fit, const Vec3D& production, const Vec3D& decay)
{
  const Vec3D d{fit.vertex[0] - decay[0], fit.vertex[1] - decay[1], fit.vertex[2] - decay[2]};
  if constexpr (requires { record.dXY; }) {
    record.dXY = static_cast<float>(std::hypot(d[0], d[1]));
  }
  record.dZ = static_cast<float>(d[2]);
  if constexpr (requires { record.dL; }) {
    const auto frame = flightFrame(production, decay);
    const auto project = [&](const Vec3D& axis, float& residual, float& sigma) {
      residual = static_cast<float>(dotOf(d, axis));
      sigma = static_cast<float>(sqrtClampedOf(decay_fit_detail::quadraticForm(axis, axis, 3, [&](int a, int b) { return fit.covarianceAt(a, b); })));
    };
    project(frame.along, record.dL, record.sL);
    project(frame.across, record.dT, record.sT);
    project(frame.third, record.dN, record.sN);
  }
}

/// The Armenteros-Podolanski variables of a two-prong fit's daughters at
/// their vertex, the first daughter the positive: alpha = (pL+ - pL-) /
/// (pL+ + pL-), pL each daughter's momentum along their sum, and qT the
/// first daughter's momentum across the sum [GeV].
struct Armenteros {
  float alpha = 0.F;
  float qT = 0.F;
};

[[nodiscard]] inline Armenteros armenteros(const DecayVertexFitResult& fit)
{
  std::array<Vec3D, 2> p{};
  for (std::size_t i = 0; i < p.size(); ++i) {
    const auto& m = fit.momenta[i];
    p[i] = decay_fit_detail::momentumOf<double>(m[0], m[1], m[2], decay_fit_detail::numeratorOf(fit.charges[i])).vector;
  }
  const Vec3D sum{p[0][0] + p[1][0], p[0][1] + p[1][1], p[0][2] + p[1][2]};
  const double norm = std::sqrt(dotOf(sum, sum));
  const double first = dotOf(p[0], sum) / norm;
  const double second = dotOf(p[1], sum) / norm;
  const auto across = crossOf(p[0], sum);
  return {.alpha = static_cast<float>((first - second) / (first + second)), .qT = static_cast<float>(std::sqrt(dotOf(across, across)) / norm)};
}

/// A particle's pseudorapidity and rapidity, its \a momentum and \a mass.
struct Rapidities {
  float eta = 0.F;
  float y = 0.F;
};

[[nodiscard]] inline Rapidities rapiditiesOf(const Vec3D& momentum, double mass)
{
  const double pt = std::hypot(momentum[0], momentum[1]);
  const double energy = std::sqrt(dotOf(momentum, momentum) + (mass * mass));
  return {.eta = static_cast<float>(pt > 0. ? etaOf(pt, momentum[2]) : std::copysign(std::numeric_limits<double>::infinity(), momentum[2])), .y = static_cast<float>(0.5 * std::log((energy + momentum[2]) / (energy - momentum[2])))};
}

/// \a result's fit into \a record's columns of the same names: its status
/// (and iterations, if the record has them); with a converged fit its chi2,
/// ndf, mass (and whether it has a point) and, with a point, its flight.
template <typename Record>
void writeCandidateFit(Record& record, const DecayCandidateResult& result)
{
  record.status = static_cast<std::uint8_t>(result.fit.status);
  if constexpr (requires { record.iterations; }) {
    record.iterations = result.fit.iterations;
  }
  if (result.fit.status != VertexFitStatus::kOk) {
    return;
  }
  record.chi2 = static_cast<float>(result.fit.chi2);
  record.ndf = static_cast<std::int16_t>(result.fit.ndf);
  record.mass = static_cast<float>(result.kinematics.mass);
  record.sMass = static_cast<float>(result.kinematics.sMass);
  if constexpr (requires { record.hasPoint; }) {
    record.hasPoint = result.hasPoint ? 1U : 0U;
  }
  if (result.hasPoint) {
    record.length = static_cast<float>(result.flight.length);
    record.sLength = static_cast<float>(result.flight.sLength);
    record.cosPointing = static_cast<float>(result.flight.cosPointing);
  }
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYSTUDYTOOLS_H_
