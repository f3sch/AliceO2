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

/// \file OverlapStudy.h
/// \brief A particle's second measurement on one layer -- an overlap of two
///        sensors -- and what the tracker's output makes of it.
///
/// An *overlap pair* is two consecutive rows of one particle's leg on the
/// same layer (its rows in path order, ParticleChains.h); on two sensors
/// where neighbouring sensors overlap, on one sensor for a split cluster or
/// a delta ray's. One record per pair of each chain (StudyChains.h, the
/// first outward leg): where it is, how far apart the two rows are, and the
/// angle between the two sensors' planes. Before a chain's pairs, one record
/// of the chain itself (`kind` 0: its rows in `group`, its distinct layers in
/// `indexA`), the denominator of every share.
///
/// With the tracker's output (`processTracks`), the same records say who
/// holds each row: nobody, the particle's best track (its good track with
/// the most rows), another good track of it (a clone), a good track of
/// another particle (a shared cluster) or a fake. A row is *claimable* where
/// a particle's best track holds a row of the same leg and layer on another
/// sensor: a track that attached its overlap partners would hold it. The
/// summary counts the fakes and clones holding claimable rows, by whether
/// the best track came from an earlier iteration than the holder (a claim at
/// the end of each iteration would have kept the row from it), the same or
/// a later one.

#ifndef GTRK_STUDY_OVERLAPSTUDY_H_
#define GTRK_STUDY_OVERLAPSTUDY_H_

#include "GenericTracking/Study/RecordWriter.h"
#include "GenericTracking/Study/StudyChains.h"
#include "GenericTracking/Study/TruthEvent.h"
#include "GenericTracking/Tracking/TrackRefit.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <format>
#include <limits>
#include <numbers>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace o2::gtrk::study
{

/// Who holds a row of an overlap pair.
enum class OverlapHolder : std::uint8_t {
  kFree,          ///< no track
  kBest,          ///< the particle's best track
  kClone,         ///< another good track of the particle
  kOtherParticle, ///< a good track of another particle: the cluster carries both
  kFake,          ///< a track with no good particle
  kNoTracks,      ///< the study ran on the truth alone
};

/// One overlap pair.
struct OverlapRecord {
  static constexpr std::string_view kName = "overlaps";
  std::uint32_t timeframe = 0U;
  std::uint8_t kind = 1U; ///< 0 a chain, 1 a pair
  std::uint64_t particle = 0U;
  std::int32_t pdg = 0;
  std::uint8_t origin = 0U; ///< ParticleOrigin
  float pt = 0.F;           ///< [GeV]
  float eta = 0.F;
  std::uint8_t layer = 0U;
  std::uint8_t group = 0U;  ///< the leg's rows on the layer
  std::uint8_t indexA = 0U; ///< the pair's first row's place in the chain, 0 the innermost
  std::uint8_t indexB = 0U;
  std::uint32_t surfaceA = 0U;
  std::uint32_t surfaceB = 0U;
  float rA = 0.F; ///< the first row's transverse radius [cm]
  float zA = 0.F; ///< [cm]
  float phiA = 0.F;
  float dR = 0.F;          ///< the second row less the first [cm]
  float dRPhi = 0.F;       ///< the azimuth's difference times the first radius [cm]
  float dZ = 0.F;          ///< [cm]
  float dist = 0.F;        ///< in space [cm]
  float cosNormals = 0.F;  ///< between the two sensors' planes' normals
  std::uint8_t found = 0U; ///< whether the particle has a good track; 0 without tracks
  std::int32_t bestTrack = -1;
  std::uint8_t bestIteration = 0xFFU;
  std::int32_t holderA = -1; ///< a track holding the row, -1 for none
  std::uint8_t classA = 0U;  ///< OverlapHolder
  std::uint8_t iterationA = 0xFFU;
  std::int32_t holderB = -1;
  std::uint8_t classB = 0U;
  std::uint8_t iterationB = 0xFFU;

  template <typename Visit>
  void columns(const Visit& visit) const
  {
    visit("timeframe", timeframe);
    visit("kind", kind);
    visit("particle", particle);
    visit("pdg", pdg);
    visit("origin", origin);
    visit("pt", pt);
    visit("eta", eta);
    visit("layer", layer);
    visit("group", group);
    visit("indexA", indexA);
    visit("indexB", indexB);
    visit("surfaceA", surfaceA);
    visit("surfaceB", surfaceB);
    visit("rA", rA);
    visit("zA", zA);
    visit("phiA", phiA);
    visit("dR", dR);
    visit("dRPhi", dRPhi);
    visit("dZ", dZ);
    visit("dist", dist);
    visit("cosNormals", cosNormals);
    visit("found", found);
    visit("bestTrack", bestTrack);
    visit("bestIteration", bestIteration);
    visit("holderA", holderA);
    visit("classA", classA);
    visit("iterationA", iterationA);
    visit("holderB", holderB);
    visit("classB", classB);
    visit("iterationB", iterationB);
  }
};

template <DetectorLayout Det>
class OverlapStudy
{
  static constexpr std::uint32_t kNone = std::numeric_limits<std::uint32_t>::max();

 public:
  OverlapStudy(RecordWriter& writer, ChainSelection selection) : mWriter(&writer), mSelection(std::move(selection)) { writer.begin<OverlapRecord>(); }

  /// The truth chains' pairs; the tracks' chains are not read (`processTracks` takes the output whole).
  void process(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source)
  {
    if (source == ChainSource::kTruth) {
      write(event, chains, nullptr);
    }
  }

  /// The truth chains' pairs with their rows' holders among \a tracks, and the
  /// fakes' and clones' claimable rows.
  void processTracks(const TruthEvent<Det>& event, const TrackStore<Det>& tracks)
  {
    const auto view = tracks.tracks().view();
    const auto rowsOf = [&](std::uint32_t track) { return tracks.rows(track); };
    Holders holders;
    holders.holder.assign(event.rows.size(), kNone);
    holders.good.assign(tracks.size(), kNoTrackParticle);
    holders.iteration.resize(tracks.size());
    for (std::uint32_t track = 0; track < tracks.size(); ++track) {
      const auto rows = rowsOf(track);
      holders.good[track] = goodParticle(event.rows, rows);
      holders.iteration[track] = view.iteration[track];
      for (const auto row : rows) {
        holders.holder[row] = holders.holder[row] == kNone ? track : holders.holder[row];
      }
      if (holders.good[track] != kNoTrackParticle) {
        auto [at, inserted] = holders.best.try_emplace(holders.good[track], track);
        if (!inserted && rowsOf(at->second).size() < rowsOf(track).size()) {
          at->second = track;
        }
      }
    }

    // Claimable rows: of any particle (the tracker claims regardless of the
    // truth), per leg and layer, the rows on sensors other than the best
    // track's rows there.
    std::vector<std::uint32_t> claimer(event.rows.size(), kNone);
    for (const auto particle : event.chains->particles()) {
      const auto best = holders.best.find(particle);
      const auto& own = event.chains->at(particle);
      if (best == holders.best.end() || !own.timed) {
        continue;
      }
      const auto bestRows = rowsOf(best->second);
      for (std::size_t i = 0; i < own.rows.size(); ++i) {
        if (std::ranges::find(bestRows, own.rows[i]) != bestRows.end()) {
          continue;
        }
        const auto surface = event.frame->measurement(own.rows[i]).surface;
        for (std::size_t j = 0; j < own.rows.size(); ++j) {
          if (j != i && own.legs[j] == own.legs[i] && event.rows[own.rows[j]].layer == event.rows[own.rows[i]].layer && event.frame->measurement(own.rows[j]).surface != surface && std::ranges::find(bestRows, own.rows[j]) != bestRows.end()) {
            claimer[own.rows[i]] = claimer[own.rows[i]] == kNone ? best->second : claimer[own.rows[i]];
            break;
          }
        }
      }
    }
    std::size_t fakes = 0U;
    for (std::uint32_t track = 0; track < tracks.size(); ++track) {
      const bool fake = holders.good[track] == kNoTrackParticle;
      const bool clone = !fake && holders.best.at(holders.good[track]) != track;
      fakes += fake ? 1U : 0U;
      mClones += clone ? 1U : 0U;
      if (!fake && !clone) {
        continue;
      }
      std::size_t claimable = 0U;
      std::array<bool, 3> when{}; // the claimer earlier, the same iteration, later
      for (const auto row : rowsOf(track)) {
        const auto by = claimer[row];
        if (by == kNone || by == track) {
          continue;
        }
        ++claimable;
        const auto a = holders.iteration[by];
        const auto b = holders.iteration[track];
        when[a < b ? 0 : (a == b ? 1 : 2)] = true;
      }
      if (claimable == 0U) {
        continue;
      }
      auto& counts = fake ? mFakeClaims : mCloneClaims;
      ++counts.tracks;
      counts.rows += claimable;
      counts.earlier += when[0] ? 1U : 0U;
      counts.same += !when[0] && when[1] ? 1U : 0U;
      counts.later += !when[0] && !when[1] ? 1U : 0U;
      counts.byRows[std::min<std::size_t>(claimable, counts.byRows.size() - 1U)] += 1U;
      counts.lengthLeft += rowsOf(track).size() - claimable;
    }
    mFakes += fakes;
    mTracks += tracks.size();

    ChainCounts counts;
    const auto chains = truthChains(event, mSelection, counts);
    write(event, chains, &holders);
  }

  [[nodiscard]] std::vector<std::string> summary() const
  {
    std::vector<std::string> lines;
    lines.push_back(std::format("overlap study: {} pairs into {} ({} on two sensors) on {} chains; {} untimed chains skipped", mPairs, mWriter->path(), mTwoSensors, mChains, mUntimed));
    if (mTracks > 0U) {
      const auto line = [&](std::string_view what, std::size_t all, const Claims& c) {
        return std::format("overlap study: {} {} of {} tracks; {} hold claimable rows ({} rows; by rows 1/2/3+: {}/{}/{}; their other rows {:.2f} on average), the best track from an earlier iteration for {}, the same for {}, a later one for {}", all, what, mTracks, c.tracks, c.rows, c.byRows[1], c.byRows[2], c.byRows[3], c.tracks > 0U ? static_cast<double>(c.lengthLeft) / static_cast<double>(c.tracks) : 0., c.earlier, c.same, c.later);
      };
      lines.push_back(line("fakes", mFakes, mFakeClaims));
      lines.push_back(line("clones", mClones, mCloneClaims));
    }
    return lines;
  }

 private:
  struct Holders {
    std::vector<std::uint32_t> holder; ///< per row, the first track holding it
    std::vector<std::uint64_t> good;   ///< per track, its good particle
    std::vector<std::uint8_t> iteration;
    std::unordered_map<std::uint64_t, std::uint32_t> best; ///< per particle, its good track with the most rows
  };

  struct Claims {
    std::size_t tracks = 0U;
    std::size_t rows = 0U;
    std::size_t earlier = 0U;
    std::size_t same = 0U;
    std::size_t later = 0U;
    std::size_t lengthLeft = 0U;
    std::array<std::size_t, 4> byRows{};
  };

  void write(const TruthEvent<Det>& event, std::span<const Chain> chains, const Holders* holders)
  {
    const auto& frame = *event.frame;
    const auto& geometry = *event.geometry;
    for (const auto& chain : chains) {
      const auto* own = event.chains->find(chain.particle);
      if (own == nullptr || !own->timed) {
        ++mUntimed;
        continue;
      }
      ++mChains;
      const auto& truth = event.particles->at(chain.particle);
      std::uint32_t bestTrack = kNone;
      if (holders != nullptr) {
        const auto best = holders->best.find(chain.particle);
        bestTrack = best != holders->best.end() ? best->second : kNone;
      }
      const bool found = bestTrack != kNone;
      const auto classify = [&](std::uint32_t row, OverlapRecord& record, bool first) {
        auto& holder = first ? record.holderA : record.holderB;
        auto& kind = first ? record.classA : record.classB;
        auto& iteration = first ? record.iterationA : record.iterationB;
        if (holders == nullptr) {
          kind = static_cast<std::uint8_t>(OverlapHolder::kNoTracks);
          return;
        }
        const auto track = holders->holder[row];
        if (track == kNone) {
          kind = static_cast<std::uint8_t>(OverlapHolder::kFree);
          return;
        }
        holder = static_cast<std::int32_t>(track);
        iteration = holders->iteration[track];
        const auto good = holders->good[track];
        const auto value = good == kNoTrackParticle ? OverlapHolder::kFake : (good != chain.particle ? OverlapHolder::kOtherParticle : (bestTrack == track ? OverlapHolder::kBest : OverlapHolder::kClone));
        kind = static_cast<std::uint8_t>(value);
      };
      std::uint64_t layers = 0U;
      for (const auto row : chain.rows) {
        layers |= std::uint64_t{1} << event.rows[row].layer;
      }
      mWriter->write(OverlapRecord{
        .timeframe = event.timeframe,
        .kind = 0U,
        .particle = chain.particle,
        .pdg = truth.pdg,
        .origin = static_cast<std::uint8_t>(truth.origin),
        .pt = truth.pt,
        .eta = truth.eta,
        .group = static_cast<std::uint8_t>(std::min<std::size_t>(chain.rows.size(), 255U)),
        .indexA = static_cast<std::uint8_t>(std::popcount(layers)),
        .found = static_cast<std::uint8_t>(found ? 1U : 0U),
        .bestTrack = found ? static_cast<std::int32_t>(bestTrack) : -1,
        .bestIteration = found ? holders->iteration[bestTrack] : static_cast<std::uint8_t>(0xFFU),
        .classA = static_cast<std::uint8_t>(holders == nullptr ? OverlapHolder::kNoTracks : OverlapHolder::kFree),
        .classB = static_cast<std::uint8_t>(holders == nullptr ? OverlapHolder::kNoTracks : OverlapHolder::kFree),
      });
      for (std::size_t i = 0; i < chain.rows.size(); ++i) {
        const auto layer = event.rows[chain.rows[i]].layer;
        bool firstOnLayer = true;
        for (std::size_t k = 0; k < i; ++k) {
          firstOnLayer = firstOnLayer && event.rows[chain.rows[k]].layer != layer;
        }
        if (!firstOnLayer) {
          continue;
        }
        std::vector<std::size_t> group;
        for (std::size_t k = i; k < chain.rows.size(); ++k) {
          if (event.rows[chain.rows[k]].layer == layer) {
            group.push_back(k);
          }
        }
        for (std::size_t g = 1; g < group.size(); ++g) {
          const auto a = frame.measurement(chain.rows[group[g - 1]]);
          const auto b = frame.measurement(chain.rows[group[g]]);
          const double ra = std::hypot(static_cast<double>(a.x), static_cast<double>(a.y));
          const double rb = std::hypot(static_cast<double>(b.x), static_cast<double>(b.y));
          const double phia = std::atan2(static_cast<double>(a.y), static_cast<double>(a.x));
          double dphi = std::atan2(static_cast<double>(b.y), static_cast<double>(b.x)) - phia;
          dphi -= dphi > std::numbers::pi ? 2. * std::numbers::pi : (dphi < -std::numbers::pi ? -2. * std::numbers::pi : 0.);
          const auto& na = geometry.planeOf(a).normal;
          const auto& nb = geometry.planeOf(b).normal;
          OverlapRecord record{
            .timeframe = event.timeframe,
            .kind = 1U,
            .particle = chain.particle,
            .pdg = truth.pdg,
            .origin = static_cast<std::uint8_t>(truth.origin),
            .pt = truth.pt,
            .eta = truth.eta,
            .layer = static_cast<std::uint8_t>(layer),
            .group = static_cast<std::uint8_t>(std::min<std::size_t>(group.size(), 255U)),
            .indexA = static_cast<std::uint8_t>(std::min<std::size_t>(group[g - 1], 255U)),
            .indexB = static_cast<std::uint8_t>(std::min<std::size_t>(group[g], 255U)),
            .surfaceA = static_cast<std::uint32_t>(a.surface),
            .surfaceB = static_cast<std::uint32_t>(b.surface),
            .rA = static_cast<float>(ra),
            .zA = a.z,
            .phiA = static_cast<float>(phia),
            .dR = static_cast<float>(rb - ra),
            .dRPhi = static_cast<float>(ra * dphi),
            .dZ = b.z - a.z,
            .dist = static_cast<float>(std::hypot(static_cast<double>(b.x - a.x), static_cast<double>(b.y - a.y), static_cast<double>(b.z - a.z))),
            .cosNormals = static_cast<float>(std::abs((static_cast<double>(na[0]) * nb[0]) + (static_cast<double>(na[1]) * nb[1]) + (static_cast<double>(na[2]) * nb[2]))),
            .found = static_cast<std::uint8_t>(found ? 1U : 0U),
            .bestTrack = found ? static_cast<std::int32_t>(bestTrack) : -1,
            .bestIteration = found ? holders->iteration[bestTrack] : static_cast<std::uint8_t>(0xFFU),
          };
          classify(chain.rows[group[g - 1]], record, true);
          classify(chain.rows[group[g]], record, false);
          mWriter->write(record);
          ++mPairs;
          mTwoSensors += a.surface != b.surface ? 1U : 0U;
        }
      }
    }
  }

  RecordWriter* mWriter = nullptr;
  ChainSelection mSelection;
  std::size_t mPairs = 0U;
  std::size_t mTwoSensors = 0U;
  std::size_t mChains = 0U;
  std::size_t mUntimed = 0U;
  std::size_t mTracks = 0U;
  std::size_t mFakes = 0U;
  std::size_t mClones = 0U;
  Claims mFakeClaims;
  Claims mCloneClaims;
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_OVERLAPSTUDY_H_
