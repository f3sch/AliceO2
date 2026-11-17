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

/// \file ParticleChains.h
/// \brief Every particle's own measurements in the order it crossed them, from
///        the MC truth: its labels, and with the MC hits their time of flight.
///
/// A particle's rows are the measurements carrying it among their labels (the
/// one truth rule, `rowsOfParticles`). They are first ordered by their r-z
/// distance from the particle's production vertex, the row as the
/// tie-break. That order falls on a looper's way back, which comes nearer the
/// vertex in r while it goes on in z and in time. With the hits, a particle
/// every one of whose rows has its crossing (`HitIndex::crossing` on the row's
/// sensor) is ordered by the crossings' time of flight instead, the row as the
/// tie-break, and is *timed*; one with a row the hits do not have keeps the
/// distance order as a whole. A timed particle's rows are split into legs
/// where the radial direction of its crossings turns: leg 0 is the first run
/// of one direction (outward for a particle that starts outward), 1 the next.

#ifndef GTRK_STUDY_PARTICLECHAINS_H_
#define GTRK_STUDY_PARTICLECHAINS_H_

#include "GenericTracking/Study/Alice/Hits.h"
#include "GenericTracking/Study/TrackingInput.h"
#include "GenericTracking/Tracking/Parallel.h"
#include "GenericTracking/Tracking/TimeFrame.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace o2::gtrk::study
{

/// One particle's rows, in path order.
struct ParticleRows {
  std::vector<std::uint32_t> rows;
  /// Timed only: each row's leg, and its crossing.
  std::vector<std::uint8_t> legs;
  std::vector<const HitSegment*> crossings;
  bool timed = false;
  /// With the hits: 0 timed, the distance order kept; 1 timed and reordered;
  /// 2 not timed; plus 4 where the first row on some layer changed.
  std::uint8_t outcome = 0U;
};

/// Every particle's `ParticleRows` of one timeframe, sorted by the particles'
/// packed identity.
class ParticleChains
{
 public:
  /// The chains of \a rows, \a frame's, whose particles' truth is
  /// \a particles; timed with \a hits if given. \a arena runs the ordering.
  template <DetectorLayout Det>
  [[nodiscard]] static ParticleChains build(std::span<const TruthRow> rows, const ParticleTable& particles, const TimeFrame<Det>& frame, const HitIndex* hits, o2::gtrk::detail::TaskArena& arena);

  /// The particles, in increasing packed identity.
  [[nodiscard]] std::span<const std::uint64_t> particles() const { return mParticles; }
  /// Whether the order came from the hits (`build` had them).
  [[nodiscard]] bool fromHits() const { return mFromHits; }

  /// \a particle's rows; none (nullptr) if no row carries it.
  [[nodiscard]] const ParticleRows* find(std::uint64_t particle) const
  {
    const auto at = std::ranges::lower_bound(mParticles, particle);
    return at != mParticles.end() && *at == particle ? &mChains[static_cast<std::size_t>(at - mParticles.begin())] : nullptr;
  }
  [[nodiscard]] bool contains(std::uint64_t particle) const { return find(particle) != nullptr; }
  /// \a particle's rows; throws std::out_of_range if no row carries it.
  [[nodiscard]] const ParticleRows& at(std::uint64_t particle) const
  {
    const auto* found = find(particle);
    if (found == nullptr) {
      throw std::out_of_range("ParticleChains: no rows carry the particle");
    }
    return *found;
  }

  /// The production vertex in r-z [cm] of a particle with \a truth, the
  /// origin for one without (nullptr): the reference of `pathDistance`.
  [[nodiscard]] static std::array<double, 2> vertexRZ(const ParticleTruth* truth)
  {
    if (truth == nullptr) {
      return {0., 0.};
    }
    return {std::hypot(static_cast<double>(truth->vx), static_cast<double>(truth->vy)), static_cast<double>(truth->vz)};
  }
  /// How far \a row is from a production \a vertex (`vertexRZ`) in r-z [cm].
  [[nodiscard]] static double pathDistance(std::span<const TruthRow> rows, std::uint32_t row, const std::array<double, 2>& vertex)
  {
    return std::hypot(static_cast<double>(rows[row].r) - vertex[0], static_cast<double>(rows[row].z) - vertex[1]);
  }

  /// The radial direction of \a segment: +1 outward, -1 inward; from its ends,
  /// or where they are at one radius (a crossing at the turning point) from
  /// its momentum. A leg's crossings all have its sign.
  [[nodiscard]] static int radialSign(const HitSegment& segment)
  {
    const double dr = std::hypot(segment.end[0], segment.end[1]) - std::hypot(segment.start[0], segment.start[1]);
    if (std::abs(dr) > 1e-4) {
      return dr > 0. ? 1 : -1;
    }
    const double along = (0.5 * (segment.start[0] + segment.end[0]) * segment.momentum[0]) + (0.5 * (segment.start[1] + segment.end[1]) * segment.momentum[1]);
    return along >= 0. ? 1 : -1;
  }

 private:
  template <DetectorLayout Det>
  static void timeOrder(std::uint64_t particle, ParticleRows& chain, std::span<const TruthRow> rows, const TimeFrame<Det>& frame, const HitIndex& hits);

  std::vector<std::uint64_t> mParticles;
  std::vector<ParticleRows> mChains;
  bool mFromHits = false;
};

template <DetectorLayout Det>
ParticleChains ParticleChains::build(std::span<const TruthRow> rows, const ParticleTable& particles, const TimeFrame<Det>& frame, const HitIndex* hits, o2::gtrk::detail::TaskArena& arena)
{
  ParticleChains result;
  result.mFromHits = hits != nullptr;
  auto rowsOf = rowsOfParticles(rows);
  result.mParticles.reserve(rowsOf.size());
  for (const auto& entry : rowsOf) {
    result.mParticles.push_back(entry.first);
  }
  std::ranges::sort(result.mParticles);
  result.mChains.resize(result.mParticles.size());
  for (std::size_t i = 0; i < result.mParticles.size(); ++i) {
    result.mChains[i].rows = std::move(rowsOf.at(result.mParticles[i]));
  }
  o2::gtrk::detail::parallelFor(arena, 0U, result.mParticles.size(), 1024, [&](std::size_t a, std::size_t b) {
    for (std::size_t i = a; i < b; ++i) {
      const auto particle = result.mParticles[i];
      auto& chain = result.mChains[i];
      const auto found = particles.find(particle);
      const auto vertex = vertexRZ(found != particles.end() ? &found->second : nullptr);
      const auto distance = [&](std::uint32_t row) { return pathDistance(rows, row, vertex); };
      std::ranges::sort(chain.rows, [&](std::uint32_t x, std::uint32_t y) { return distance(x) < distance(y) || (!(distance(y) < distance(x)) && x < y); });
      if (hits != nullptr) {
        timeOrder(particle, chain, rows, frame, *hits);
      }
    }
  });
  return result;
}

template <DetectorLayout Det>
void ParticleChains::timeOrder(std::uint64_t particle, ParticleRows& chain, std::span<const TruthRow> rows, const TimeFrame<Det>& frame, const HitIndex& hits)
{
  constexpr auto kNone = std::numeric_limits<std::uint32_t>::max();
  constexpr auto kLayers = static_cast<std::size_t>(Detector<Det>::nLayers);
  auto& own = chain.rows;
  std::vector<std::pair<double, std::uint32_t>> keyed;
  std::vector<const HitSegment*> segments;
  keyed.reserve(own.size());
  segments.reserve(own.size());
  for (const auto row : own) {
    const auto point = frame.measurement(row);
    const auto* segment = hits.crossing(particle, point.surface, {static_cast<double>(point.x), static_cast<double>(point.y), static_cast<double>(point.z)});
    if (segment == nullptr) {
      chain.outcome = 2U;
      return;
    }
    keyed.emplace_back(segment->time, row);
    segments.push_back(segment);
  }
  const auto firstOnLayers = [&] {
    std::array<std::uint32_t, kLayers> first{};
    first.fill(kNone);
    for (const auto row : own) {
      auto& at = first[rows[row].layer];
      at = at == kNone ? row : at;
    }
    return first;
  };
  const auto before = firstOnLayers();
  std::vector<std::size_t> order(own.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::ranges::sort(order, [&](std::size_t x, std::size_t y) { return keyed[x] < keyed[y]; });
  bool same = true;
  chain.crossings.resize(own.size());
  for (std::size_t k = 0; k < own.size(); ++k) {
    same = same && own[k] == keyed[order[k]].second;
    own[k] = keyed[order[k]].second;
    chain.crossings[k] = segments[order[k]];
  }
  chain.outcome = static_cast<std::uint8_t>((same ? 0U : 1U) | (before != firstOnLayers() ? 4U : 0U));
  chain.legs.resize(own.size());
  std::uint8_t leg = 0U;
  for (std::size_t k = 0; k < own.size(); ++k) {
    if (k != 0U && radialSign(*chain.crossings[k]) != radialSign(*chain.crossings[k - 1U]) && leg < 255U) {
      ++leg;
    }
    chain.legs[k] = leg;
  }
  chain.timed = true;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_PARTICLECHAINS_H_
