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
///        into a record.

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
