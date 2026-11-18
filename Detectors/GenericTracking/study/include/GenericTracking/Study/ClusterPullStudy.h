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

/// \file ClusterPullStudy.h
/// \brief Every row of a chain held to its particle's MC crossing of the
///        row's sensor: the cluster's position less the crossing's on the
///        sensor's plane, in the cluster's own errors.
///
/// For each chain (StudyChains.h) and each of its rows: the particle's MC
/// segment through the row's sensor nearest the cluster (`HitIndex`), crossed
/// with the plane the fits bind the row to (TrackingGeometry::planeOf); the
/// residuals in loc0 and loc1 and the cluster's sigmas (the measurement's
/// varLoc0, varLoc1), the incidence (the angle between the MC direction and
/// the plane's normal), the true momentum there, the cluster's radius and z.
/// A row whose particle has no segment there, or a segment along the plane,
/// gives no record.
///
/// With the next row's segment too, the scattering between: theta's change
/// from this segment's start to the next one's (a uniform solenoid field
/// leaves theta alone; the material does not), the material table's x/X0 on
/// the straight chord between the two starts, and Highland's theta0 for it
/// at the particle's momentum and mass (the projected angle's width the fits
/// take): their ratio is a unit Gaussian if the model holds.

#ifndef GTRK_STUDY_CLUSTERPULLSTUDY_H_
#define GTRK_STUDY_CLUSTERPULLSTUDY_H_

#include "GenericTracking/Study/RecordWriter.h"
#include "GenericTracking/Study/StateTruth.h"
#include "GenericTracking/Study/StudyChains.h"
#include "GenericTracking/Study/TruthEvent.h"
#include "GenericTracking/Tracking/Material.h"

#include "SimulationDataFormat/O2DatabasePDG.h"

#include <cmath>
#include <cstdint>
#include <format>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace o2::gtrk::study
{

/// One row against the MC.
struct ClusterPullRecord {
  static constexpr std::string_view kName = "clusterPull";
  std::uint32_t timeframe = 0U;
  std::uint64_t particle = 0U;
  std::int32_t pdg = 0;
  std::uint8_t row = 0U;      ///< its place in the chain, 0 the innermost
  std::uint32_t surface = 0U; ///< the sensor
  float r = 0.F;              ///< the cluster's transverse radius [cm]
  float z = 0.F;              ///< [cm]
  float p = 0.F;              ///< the true momentum at the crossing [GeV]
  float cosIncidence = 0.F;   ///< |cos| of the angle between the MC direction and the plane's normal
  float dLoc0 = 0.F;          ///< the cluster less the crossing [cm]
  float dLoc1 = 0.F;
  float sLoc0 = 0.F; ///< the cluster's [cm]
  float sLoc1 = 0.F;
  float dTheta = 0.F;    ///< theta at the next row's segment's start less here's [rad]; 0 at the last row
  float chordX2X0 = 0.F; ///< the material table's x/X0 between the two starts
  float theta0 = 0.F;    ///< Highland's for it [rad]; 0 without a next row

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("particle", particle);
    visit("pdg", pdg);
    visit("row", row);
    visit("surface", surface);
    visit("r", r);
    visit("z", z);
    visit("p", p);
    visit("cosIncidence", cosIncidence);
    visit("dLoc0", dLoc0);
    visit("dLoc1", dLoc1);
    visit("sLoc0", sLoc0);
    visit("sLoc1", sLoc1);
    visit("dTheta", dTheta);
    visit("chordX2X0", chordX2X0);
    visit("theta0", theta0);
  }
};

template <DetectorLayout Det>
class ClusterPullStudy
{
 public:
  explicit ClusterPullStudy(RecordWriter& writer) : mWriter(&writer) { writer.begin<ClusterPullRecord>(); }

  /// \a chains of \a event.
  void process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource /*source*/)
  {
    if (event.hits == nullptr) {
      return;
    }
    const auto& frame = *event.frame;
    const auto& geometry = *event.geometry;
    const auto segmentOf = [&](std::uint64_t particle, std::uint32_t row) {
      const auto m = frame.measurement(row);
      return event.hits->crossing(particle, m.surface, {static_cast<double>(m.x), static_cast<double>(m.y), static_cast<double>(m.z)});
    };
    for (const auto& chain : chains) {
      const auto& truth = event.particles->at(chain.particle);
      const auto* species = o2::O2DatabasePDG::Instance()->GetParticle(truth.pdg);
      const double mass = species == nullptr ? 0. : species->Mass();
      for (std::size_t i = 0; i < chain.rows.size(); ++i) {
        const auto m = frame.measurement(chain.rows[i]);
        const auto* segment = segmentOf(chain.particle, chain.rows[i]);
        if (segment == nullptr) {
          ++mMissing;
          continue;
        }
        const auto plane = geometry.planeOf(m);
        const auto crossing = segmentCrossing(plane, segment->start, segment->end);
        if (!crossing) {
          ++mMissing;
          continue;
        }
        const auto& mom = segment->momentum;
        const double momentum = std::hypot(mom[0], mom[1], mom[2]);
        const double along = std::abs((mom[0] * static_cast<double>(plane.normal[0])) + (mom[1] * static_cast<double>(plane.normal[1])) + (mom[2] * static_cast<double>(plane.normal[2])));
        float dTheta = 0.F;
        float chordX2X0 = 0.F;
        float theta0 = 0.F;
        if (i + 1 < chain.rows.size()) {
          if (const auto* next = segmentOf(chain.particle, chain.rows[i + 1]); next != nullptr) {
            const auto thetaOf = [](const std::array<double, 3>& q) { return std::atan2(std::hypot(q[0], q[1]), q[2]); };
            dTheta = static_cast<float>(thetaOf(next->momentum) - thetaOf(mom));
            const auto& a = segment->start;
            const auto& b = next->start;
            chordX2X0 = event.material->getMatBudget(static_cast<float>(a[0]), static_cast<float>(a[1]), static_cast<float>(a[2]), static_cast<float>(b[0]), static_cast<float>(b[1]), static_cast<float>(b[2])).meanX2X0;
            theta0 = highlandTheta0(chordX2X0, chordX2X0, static_cast<Scalar>(momentum), static_cast<Scalar>(mass));
          }
        }
        mWriter->write(ClusterPullRecord{
          .timeframe = event.timeframe,
          .particle = chain.particle,
          .pdg = truth.pdg,
          .row = static_cast<std::uint8_t>(std::min<std::size_t>(i, 255U)),
          .surface = static_cast<std::uint32_t>(m.surface),
          .r = static_cast<float>(std::hypot(static_cast<double>(m.x), static_cast<double>(m.y))),
          .z = m.z,
          .p = static_cast<float>(momentum),
          .cosIncidence = static_cast<float>(along / momentum),
          .dLoc0 = static_cast<float>(static_cast<double>(m.loc0) - (*crossing)[0]),
          .dLoc1 = static_cast<float>(static_cast<double>(m.loc1) - (*crossing)[1]),
          .sLoc0 = std::sqrt(static_cast<float>(m.varLoc0)),
          .sLoc1 = std::sqrt(static_cast<float>(m.varLoc1)),
          .dTheta = dTheta,
          .chordX2X0 = chordX2X0,
          .theta0 = theta0,
        });
        ++mRecords;
      }
    }
  }

  [[nodiscard]] std::vector<std::string> summary() const
  {
    return {std::format("cluster pull study: {} rows into {}; {} without the particle's MC crossing", mRecords, mWriter->path(), mMissing)};
  }

 private:
  RecordWriter* mWriter = nullptr;
  std::size_t mRecords = 0U;
  std::size_t mMissing = 0U;
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_CLUSTERPULLSTUDY_H_
