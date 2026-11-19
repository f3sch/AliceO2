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

/// \file StudyChains.h
/// \brief What a truth study fits: chains, each a particle and its rows, from
///        the MC truth or from the tracker's output.
///
/// A study does not know where its chains come from. A truth chain is a
/// particle's own rows (`ParticleChains`): its first leg, when the hits timed
/// it and the leg runs outward, so that its rows in time order are its
/// motion, innermost first, which the fits' seed and energy loss assume. A
/// track chain is a good track's rows in the track's order, with what the
/// tracker made of them (`TrackContext`). One selection (`ChainSelection`)
/// picks the particles for both.

#ifndef GTRK_STUDY_STUDYCHAINS_H_
#define GTRK_STUDY_STUDYCHAINS_H_

#include "GenericTracking/Tracking/TrackRefit.h"
#include "GenericTracking/Tracking/TrackStore.h"
#include "GenericTracking/Study/ParticleChains.h"
#include "GenericTracking/Study/StateTruth.h"
#include "GenericTracking/Study/TrackComparison.h"
#include "GenericTracking/Study/TruthEvent.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace o2::gtrk::study
{

/// Where a study's chains come from.
enum class ChainSource : std::uint8_t { kTruth,
                                        kTracks };

/// What a chain from the tracker's output carries besides its rows.
struct TrackContext {
  std::uint32_t track = 0U;         ///< into the tracker's output
  BoundTrackState state;            ///< the track's fitted state, on its innermost surface
  BoundTrackState outward;          ///< its outward pass's state
  std::uint32_t comparisonRow = 0U; ///< the row on `state`'s surface
  std::uint8_t kinks = 0U;
};

/// A particle and its rows, in the order of its motion, innermost first.
struct Chain {
  std::uint64_t particle = 0U; ///< the truth chain's particle, or the track's good particle
  std::span<const std::uint32_t> rows;
  std::uint8_t leg = 0U;             ///< a truth chain's leg; 0 for a track
  std::optional<TrackContext> track; ///< a chain from the tracker's output only
};

/// Each particle's longest chain of \a chains, the first of equal length.
[[nodiscard]] inline std::unordered_map<std::uint64_t, const Chain*> longestChains(std::span<const Chain> chains)
{
  std::unordered_map<std::uint64_t, const Chain*> longest;
  for (const auto& chain : chains) {
    auto& best = longest[chain.particle];
    if (best == nullptr || chain.rows.size() > best->rows.size()) {
      best = &chain;
    }
  }
  return longest;
}

/// What names a chain in a study's records.
struct ChainIdentity {
  std::uint32_t timeframe = 0U;
  std::uint8_t source = 0U; ///< ChainSource
  std::int32_t track = -1;  ///< the tracker's track, -1 for a truth chain
  std::uint64_t particle = 0U;
  std::uint8_t leg = 0U;
};

/// \a chain's identity in timeframe \a timeframe, from \a source.
[[nodiscard]] inline ChainIdentity identityOf(const Chain& chain, std::uint32_t timeframe, ChainSource source)
{
  return {.timeframe = timeframe, .source = static_cast<std::uint8_t>(source), .track = chain.track ? static_cast<std::int32_t>(chain.track->track) : -1, .particle = chain.particle, .leg = chain.leg};
}

/// The particles a study takes, from either source.
struct ChainSelection {
  /// A bit per `ParticleOrigin`: physical primaries and decay daughters.
  std::uint8_t origins = (1U << static_cast<unsigned>(ParticleOrigin::kPhysicalPrimary)) | (1U << static_cast<unsigned>(ParticleOrigin::kDecay));
  std::vector<std::int32_t> pdgs; ///< |PDG codes|; empty, every species
  double minPtGeV = 0.;
  std::uint32_t minRows = 0U; ///< on the chain
  std::uint64_t layers = 0U;  ///< a bit per layer the particle must have a row on
  bool findable = false;      ///< the QA's findable rule too (`isFindable`, a chain of four layers above `minPtGeV`)
};

/// The selection \a text names: `key=value` items separated by `;`, every
/// key optional: `origins=primary+decay+material`, `pdg=11+211`,
/// `minPt=0.1`, `minRows=3`, `layers=<mask>`, and `findable` alone. Throws
/// std::invalid_argument for another key or a value it cannot read.
[[nodiscard]] ChainSelection parseChainSelection(std::string_view text);

/// What a producer took and left.
struct ChainCounts {
  std::size_t selected = 0U; ///< particles the selection takes
  std::size_t chains = 0U;
  std::size_t untimed = 0U;  ///< truth: a row without its crossing
  std::size_t inward = 0U;   ///< truth: the first leg runs inward
  std::size_t tooShort = 0U; ///< fewer rows than `minRows`
  std::size_t noTruth = 0U;  ///< tracks: no good particle, no row on the state's surface, or no crossing there
};

/// \a counts in one line, \a label first.
[[nodiscard]] std::string describe(const ChainCounts& counts, std::string_view label);

namespace study_chains_detail
{
/// Whether \a selection takes the particle of \a truth with \a rows, its rows
/// in \a all; \a findable its chains by the QA's rule, if asked.
[[nodiscard]] inline bool selected(const ChainSelection& selection, const ParticleTruth& truth, std::span<const std::uint32_t> rows, std::span<const TruthRow> all, const std::unordered_map<std::uint64_t, ParticleLayers>* findable, const ComparisonSettings* rule, std::uint64_t particle)
{
  if ((selection.origins & (1U << static_cast<unsigned>(truth.origin))) == 0U || static_cast<double>(truth.pt) < selection.minPtGeV) {
    return false;
  }
  if (!selection.pdgs.empty() && std::ranges::find(selection.pdgs, std::abs(truth.pdg)) == selection.pdgs.end()) {
    return false;
  }
  if (selection.layers != 0U) {
    std::uint64_t layers = 0U;
    for (const auto row : rows) {
      layers |= std::uint64_t{1} << all[row].layer;
    }
    if ((layers & selection.layers) != selection.layers) {
      return false;
    }
  }
  if (findable != nullptr) {
    const auto chain = findable->find(particle);
    return chain != findable->end() && isFindable(truth.pt, chain->second, *rule);
  }
  return true;
}
} // namespace study_chains_detail

/// The truth chains of \a event's particles \a selection takes, in particle
/// order; what was taken and left into \a counts.
template <DetectorLayout Det>
[[nodiscard]] std::vector<Chain> truthChains(const TruthEvent<Det>& event, const ChainSelection& selection, ChainCounts& counts)
{
  const auto rule = comparisonSettings<Det>(4, selection.minPtGeV);
  const auto findable = selection.findable ? particleLayers(event.rows, rule) : std::unordered_map<std::uint64_t, ParticleLayers>{};
  std::vector<Chain> chains;
  for (const auto particle : event.chains->particles()) {
    const auto found = event.particles->find(particle);
    if (found == event.particles->end()) {
      continue;
    }
    const auto& own = event.chains->at(particle);
    if (!study_chains_detail::selected(selection, found->second, own.rows, event.rows, selection.findable ? &findable : nullptr, &rule, particle)) {
      continue;
    }
    ++counts.selected;
    if (!own.timed) {
      ++counts.untimed;
      continue;
    }
    if (ParticleChains::radialSign(*own.crossings.front()) < 0) {
      ++counts.inward;
      continue;
    }
    std::size_t n = 0U;
    while (n < own.rows.size() && own.legs[n] == 0U) {
      ++n;
    }
    if (n < selection.minRows) {
      ++counts.tooShort;
      continue;
    }
    chains.push_back({.particle = particle, .rows = std::span<const std::uint32_t>{own.rows.data(), n}});
  }
  counts.chains += chains.size();
  return chains;
}

/// The chains of \a tracks' good tracks whose particle \a selection takes,
/// in the tracks' order: those with a row on the fitted state's surface and
/// the particle's MC crossing there (`stateTruth`); what was taken and left
/// into \a counts. The chains point into \a tracks.
template <DetectorLayout Det>
[[nodiscard]] std::vector<Chain> trackChains(const TruthEvent<Det>& event, const TrackStore<Det>& tracks, const ChainSelection& selection, ChainCounts& counts)
{
  const auto rule = comparisonSettings<Det>(4, selection.minPtGeV);
  const auto findable = selection.findable ? particleLayers(event.rows, rule) : std::unordered_map<std::uint64_t, ParticleLayers>{};
  const auto view = tracks.tracks().view();
  std::vector<Chain> chains;
  for (std::size_t track = 0; track < tracks.size(); ++track) {
    const auto rows = tracks.rows(track);
    const auto particle = goodParticle(event.rows, rows);
    const auto found = particle == kNoTrackParticle ? event.particles->end() : event.particles->find(particle);
    if (found == event.particles->end()) {
      ++counts.noTruth;
      continue;
    }
    const auto* own = event.chains->find(particle);
    if (own == nullptr || !study_chains_detail::selected(selection, found->second, own->rows, event.rows, selection.findable ? &findable : nullptr, &rule, particle)) {
      continue;
    }
    ++counts.selected;
    const auto& state = view.state[track];
    const auto innermost = std::ranges::find_if(rows, [&](std::uint32_t row) { return event.frame->measurement(row).surface == state.surface; });
    if (innermost == rows.end() || event.hits == nullptr || !stateTruth(*event.geometry, *event.hits, *event.frame, state, particle, *innermost, particleCharge(found->second.pdg))) {
      ++counts.noTruth;
      continue;
    }
    if (rows.size() < selection.minRows) {
      ++counts.tooShort;
      continue;
    }
    chains.push_back({
      .particle = particle,
      .rows = rows,
      .track = TrackContext{
        .track = static_cast<std::uint32_t>(track),
        .state = state,
        .outward = view.outwardState[track],
        .comparisonRow = *innermost,
        .kinks = view.kinks[track],
      },
    });
  }
  counts.chains += chains.size();
  return chains;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_STUDYCHAINS_H_
