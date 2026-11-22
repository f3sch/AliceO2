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

/// \file DaughterTable.h
/// \brief The timeframe's committed tracks refitted under the species a
///        decay finder asks for: one entry per (track, species), shared by
///        every candidate that takes it.
///
/// The tracking fits every track as a pion; a decay's fit needs each
/// daughter under its own species (its scattering and loss), and a finder
/// without PID tries several species per track. An *entry* is one committed
/// track's rows (its references in the store: its layers' rows and its
/// overlap partners, in path order) refitted under one species
/// (`refitRowsWithHypothesis`, an electron by the Gaussian-sum filter): its
/// state on the innermost measurement, the fit's chi2 and outcome, an
/// electron's mixture. Every species is made the same way, the pion too, so
/// two hypotheses of one track differ by the mass only. The innermost
/// measurement's plane, the propagator, the mass and |q| are not stored:
/// they follow from the track's first row, the run's field and the species.
///
/// The table is filled on request, in phases: a finder marks the
/// (track, species) pairs it needs (`request`, safe from any number of
/// threads), `flush` numbers the new marks after the entries held -- in
/// (track, species) order, so the numbering does not depend on the order of
/// the requests -- and refits them, and the entries are read by index until
/// the next flush grows the columns. An entry is immutable once filled; a
/// failed refit stays an entry with its status, pass and measurement, never
/// retried. A species is the table's by its index in the run's list, which
/// the table is made with. The memory: tracks x species map entries, one
/// entry per pair asked for, one mixture per electron entry, all from the
/// caller's resource; `begin` starts a timeframe.

#ifndef GTRK_TRACKING_DAUGHTERTABLE_H_
#define GTRK_TRACKING_DAUGHTERTABLE_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/Columns.h"
#include "GenericTracking/Base/HostAlgorithms.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/GSFTrackFitter.h"
#include "GenericTracking/Tracking/HypothesisRefit.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/TrackFitterUtils.h"
#include "GenericTracking/Tracking/TrackStore.h"

#include "ReconstructionDataFormats/PID.h"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory_resource>
#include <span>
#include <stdexcept>
#include <vector>

namespace o2::gtrk
{

/// One entry of a `DaughterTable`.
struct DaughterEntry {
  // As for TrackletT: upstream's `apply` macro also uses the member names as
  // enum constants.
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
  template <template <class> class C>
  struct RecordT {
    MEMLAYOUT_APPLY(RecordT, track, species, state, chi2, ndf, status, failedPass, failedMeasurement, mixture)
    C<std::uint32_t> track;   ///< the committed track
    C<std::uint8_t> species;  ///< its index in the table's species
    C<BoundTrackState> state; ///< on the innermost measurement; an electron's its mixture's collapse
    C<Scalar> chi2;           ///< the inward pass's
    C<std::int16_t> ndf;
    C<PropagationStatus> status; ///< the refit's: kOk, or why it failed
    C<FitPass> failedPass;       ///< where it failed
    C<std::int16_t> failedMeasurement;
    C<std::int32_t> mixture; ///< its index in the mixtures, -1 for another species than the electron
  };
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
};

/// The timeframe's (track, species) refits (see the file's comment).
template <DetectorLayout Det>
class DaughterTable
{
 public:
  static constexpr int kComponents = BetheHeitlerLoss::kComponents;
  using Mixture = GSFMixture<kComponents>;
  using Entries = BoundedColumns<DaughterEntry::template RecordT>;
  using Entry = typename Entries::Record;
  /// `entryOf`: no entry and none asked for.
  static constexpr std::int32_t kNone = -1;

  /// A table of \a species, the run's list (at most 255, each once);
  /// std::invalid_argument otherwise.
  explicit DaughterTable(std::span<const o2::track::PID> species) : mSpecies(species.begin(), species.end())
  {
    if (mSpecies.empty() || mSpecies.size() > std::numeric_limits<std::uint8_t>::max()) {
      throw std::invalid_argument("daughter table: one to 255 species");
    }
    for (std::size_t a = 0; a < mSpecies.size(); ++a) {
      for (std::size_t b = 0; b < a; ++b) {
        if (mSpecies[a] == mSpecies[b]) {
          throw std::invalid_argument("daughter table: a species listed twice");
        }
      }
    }
  }

  /// Every allocation on \a resource, what was held released.
  void setMemoryResource(std::pmr::memory_resource* resource)
  {
    mResource = resource;
    deepVectorClear(mEntryOf, resource);
    deepVectorClear(mMixtures, resource);
    mEntries.setMemoryResource(resource);
  }

  /// A timeframe of \a tracks committed tracks: no entries, nothing asked
  /// for. kIndexOverflow where tracks x species exceeds the entries' index.
  [[nodiscard]] StageStatus begin(std::size_t tracks)
  {
    mEntries.resize(0U);
    mMixtures.clear();
    mTracks = 0U;
    const std::uint64_t pairs = static_cast<std::uint64_t>(tracks) * mSpecies.size();
    if (pairs > static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max())) {
      mEntryOf.clear();
      return StageStatus::kIndexOverflow;
    }
    return guardOutOfMemory([&] {
      mEntryOf.assign(static_cast<std::size_t>(pairs), kNone);
      mTracks = tracks;
      return StageStatus::kComplete;
    });
  }

  [[nodiscard]] std::span<const o2::track::PID> species() const { return mSpecies; }
  /// \a pid's index in the species, -1 if the table has none of it.
  [[nodiscard]] int speciesIndex(o2::track::PID pid) const
  {
    const auto found = std::ranges::find(mSpecies, pid);
    return found == mSpecies.end() ? -1 : static_cast<int>(found - mSpecies.begin());
  }
  [[nodiscard]] std::size_t tracks() const { return mTracks; }

  /// Ask for \a track under species \a s; from any thread. Nothing if its
  /// entry exists or is asked for.
  void request(std::uint32_t track, std::uint8_t s)
  {
    const std::atomic_ref<std::int32_t> slot{mEntryOf[key(track, s)]};
    std::int32_t none = kNone;
    slot.compare_exchange_strong(none, kRequested, std::memory_order_relaxed);
  }

  /// The pairs asked for since the last flush, numbered after the entries in
  /// (track, species) order and refitted from \a tracks' rows in \a frame
  /// with \a propagator's field on \a forEach(n, body) (body(i) for every i
  /// below n, serially or in parallel; \a chunk the compaction's chunk).
  /// kOutOfMemory where the pool refuses the columns, with no entry added.
  template <typename ForEach>
  [[nodiscard]] StageStatus flush(const ForEach& forEach, std::uint32_t chunk, const Propagator& propagator, const TrackingGeometry<Det>& geometry, const TimeFrame<Det>& frame, const TrackStore<Det>& tracks)
  {
    return guardOutOfMemory([&] {
      BoundedVector<std::uint32_t> asked{std::pmr::polymorphic_allocator<std::uint32_t>{mResource}};
      host_algorithms::compactIndices(forEach, static_cast<std::uint32_t>(mEntryOf.size()), chunk, [&](std::uint32_t k) { return mEntryOf[k] == kRequested; }, asked);
      const auto n = static_cast<std::uint32_t>(asked.size());
      // The electrons' mixtures, one each, in the entries' order.
      BoundedVector<std::uint32_t> mixtureOffsets(std::size_t{n} + 1U, 0U, std::pmr::polymorphic_allocator<std::uint32_t>{mResource});
      const auto isElectron = [&](std::uint32_t i) { return mSpecies[asked[i] % mSpecies.size()] == o2::track::PID::Electron ? 1U : 0U; };
      if (!host_algorithms::exclusiveScan(forEach, n, chunk, isElectron, std::span<std::uint32_t>{mixtureOffsets}, mResource)) {
        return StageStatus::kIndexOverflow;
      }
      const std::size_t base = mEntries.size();
      const std::size_t mixtureBase = mMixtures.size();
      mEntries.resize(base + n);
      try {
        mMixtures.resize(mixtureBase + mixtureOffsets[n]);
      } catch (...) {
        mEntries.resize(base);
        throw;
      }
      const auto view = mEntries.view();
      // A chain of a track's rows may hold two on one layer (its overlap
      // partners): no material on the step between them, which the layer's
      // band holds already.
      constexpr detail::LayerMaterial material{};
      forEach(n, [&](std::uint32_t i) {
        const std::uint32_t k = asked[i];
        const auto track = static_cast<std::uint32_t>(k / mSpecies.size());
        const auto s = static_cast<std::uint8_t>(k % mSpecies.size());
        const auto at = base + i;
        const auto rows = tracks.rows(track);
        const auto fit = refitRowsWithHypothesis<kComponents>(propagator, geometry, material, RowChain<Det>{.frame = &frame, .rows = rows.data()}, static_cast<int>(rows.size()), mSpecies[s]);
        view.track[at] = track;
        view.species[at] = s;
        view.state[at] = fit.innermost;
        view.chi2[at] = fit.chi2;
        view.ndf[at] = fit.ndf;
        view.status[at] = fit.status;
        view.failedPass[at] = fit.failedPass;
        view.failedMeasurement[at] = fit.failedMeasurement;
        view.mixture[at] = kNone;
        if (isElectron(i) != 0U) {
          const auto mixture = mixtureBase + mixtureOffsets[i];
          mMixtures[mixture] = fit.innermostMixture;
          view.mixture[at] = static_cast<std::int32_t>(mixture);
        }
        mEntryOf[k] = static_cast<std::int32_t>(at);
      });
      return StageStatus::kComplete;
    });
  }

  /// \a track's entry under species \a s; kNone without one (none asked for,
  /// or asked for and not flushed yet: kNone too).
  [[nodiscard]] std::int32_t entryOf(std::uint32_t track, std::uint8_t s) const
  {
    const auto entry = mEntryOf[key(track, s)];
    return entry >= 0 ? entry : kNone;
  }
  [[nodiscard]] std::size_t size() const { return mEntries.size(); }
  [[nodiscard]] Entry entry(std::size_t i) const { return mEntries.read(i); }
  [[nodiscard]] typename Entries::ConstView view() const { return mEntries.view(); }
  /// Electron entry \a i's mixture on its innermost measurement.
  [[nodiscard]] const Mixture& mixture(std::size_t i) const { return mMixtures[i]; }

 private:
  /// A pair asked for, not flushed yet.
  static constexpr std::int32_t kRequested = -2;

  [[nodiscard]] std::size_t key(std::uint32_t track, std::uint8_t s) const { return (static_cast<std::size_t>(track) * mSpecies.size()) + s; }

  std::vector<o2::track::PID> mSpecies;
  std::pmr::memory_resource* mResource = std::pmr::get_default_resource();
  std::size_t mTracks = 0U;
  BoundedVector<std::int32_t> mEntryOf;
  Entries mEntries;
  BoundedVector<Mixture> mMixtures;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_DAUGHTERTABLE_H_
