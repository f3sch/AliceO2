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

/// \file TrackStore.h
/// \brief The timeframe's committed tracks: what each refers to, and the
///        tables other sources add to them.
///
/// A *track table* holds, for every committed track, a range of records of
/// one memlayout struct: one track to any number of records, the records in
/// struct-of-arrays columns, each track's range `{first, count}` into them.
/// It is filled in batches by count-then-fill (count per track, exclusive
/// scan, place), so its size is exact for any input. A track committed again
/// gets a new range after the records held; its old one is no longer read and
/// goes with the rest at the timeframe's end.
///
/// The store holds one table of its own, the tracks' *references*: per
/// reference the row it names, the source whose measurements that row indexes
/// (`kOwnRows` for the tracker's frame, or a source the store names with
/// `addSource`) and its role (`ReferenceRole`). Any other component may add a
/// table of its own record type (`addTable`): a matching step's matches to its
/// detector's measurements, say, keyed by the committed track. The tracker
/// does not know those types; it grows every table with the tracks it commits.
///
/// Each table's ranges and records are plain arrays: a device holds the same
/// columns and reads them through the same views.

#ifndef GTRK_TRACKING_TRACKSTORE_H_
#define GTRK_TRACKING_TRACKSTORE_H_

#include "GenericTracking/Base/BoundedAllocator.h"
#include "GenericTracking/Base/Columns.h"
#include "GenericTracking/Base/HostAlgorithms.h"
#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/StageOutput.h"
#include "GenericTracking/Tracking/TrackRefit.h"

#include "GPUCommonDef.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <memory_resource>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace o2::gtrk
{

/// Which measurements a reference's row indexes.
using SourceID = std::uint8_t;
/// The tracker's own frame: the rows its tracks are made of.
inline constexpr SourceID kOwnRows = 0U;

/// One reference of a committed track.
struct TrackReference {
  // As for TrackletT: upstream's `apply` macro also uses the member names as
  // enum constants.
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
  template <template <class> class C>
  struct RecordT {
    MEMLAYOUT_APPLY(RecordT, row, source, role)
    C<std::uint32_t> row; ///< the row in its source's measurements
    C<SourceID> source;
    C<std::uint8_t> role; ///< ReferenceRole bits
  };
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members,cppcoreguidelines-missing-std-forward,readability-identifier-naming)
};

/// A committed track's records in a track table: [first, first + count).
struct TrackRange {
  std::uint32_t first = 0U;
  std::uint32_t count = 0U;
};

namespace detail
{

/// What the store does to each table it holds, whatever its records.
class TrackTableBase
{
 public:
  TrackTableBase() = default;
  TrackTableBase(const TrackTableBase&) = delete;
  TrackTableBase& operator=(const TrackTableBase&) = delete;
  TrackTableBase(TrackTableBase&&) = delete;
  TrackTableBase& operator=(TrackTableBase&&) = delete;
  virtual ~TrackTableBase() = default;

  /// Every allocation reseated on \a resource, what was held released.
  virtual void setMemoryResource(std::pmr::memory_resource* resource) = 0;
  /// No tracks and no records; the storage kept.
  virtual void clear() = 0;
  /// Ranges for \a tracks tracks, the new ones empty.
  virtual void resizeTracks(std::size_t tracks) = 0;
};

/// A table type's key in a store: the address of its own variable.
template <template <template <class> class> class S>
inline constexpr char kTableKey = 0;

} // namespace detail

/// One track to any number of records of \a S (see the file's comment).
template <template <template <class> class> class S>
class TrackTable final : public detail::TrackTableBase
{
 public:
  using Records = BoundedColumns<S>;
  using Record = typename Records::Record;

  /// What a kernel reads: each track's range, and the records' columns.
  struct ConstView {
    const TrackRange* ranges = nullptr;
    typename Records::ConstView records;
  };

  void setMemoryResource(std::pmr::memory_resource* resource) override
  {
    mResource = resource;
    deepVectorClear(mRanges, resource);
    mRecords.setMemoryResource(resource);
  }

  void clear() override
  {
    mRanges.clear();
    mRecords.resize(0U);
  }

  void resizeTracks(std::size_t tracks) override { mRanges.resize(tracks, TrackRange{}); }

  [[nodiscard]] std::size_t trackCount() const { return mRanges.size(); }
  /// Every record held, those of superseded ranges too.
  [[nodiscard]] std::size_t recordCount() const { return mRecords.size(); }
  [[nodiscard]] TrackRange range(std::size_t track) const { return mRanges[track]; }
  /// Record \a k of \a track's range.
  [[nodiscard]] Record read(std::size_t track, std::uint32_t k) const
  {
    assert(k < mRanges[track].count);
    return mRecords.read(mRanges[track].first + k);
  }
  [[nodiscard]] ConstView view() const { return {.ranges = mRanges.data(), .records = mRecords.view()}; }
  [[nodiscard]] const Records& records() const { return mRecords; }

  /// New ranges for \a n tracks, the i-th \a trackOf(i) (all distinct and
  /// below `trackCount`), after the records held: \a count(i) records each,
  /// written by \a fill(i, view, at) to the columns' rows from `at` on.
  /// Count, scan and place on \a forEach(n, body), which runs body(i) for
  /// every i below n, serially or in parallel; \a chunk is the scan's chunk.
  /// kIndexOverflow where the records would not fit a 32-bit index.
  template <typename ForEach, typename TrackOf, typename Count, typename Fill>
  [[nodiscard]] StageStatus append(const ForEach& forEach, std::uint32_t chunk, std::uint32_t n, const TrackOf& trackOf, const Count& count, const Fill& fill)
  {
    return guardOutOfMemory([&] {
      BoundedVector<std::uint32_t> offsets(std::size_t{n} + 1U, std::uint32_t{0}, std::pmr::polymorphic_allocator<std::uint32_t>{mResource});
      if (!host_algorithms::exclusiveScan(forEach, n, chunk, count, std::span<std::uint32_t>{offsets}, mResource)) {
        return StageStatus::kIndexOverflow;
      }
      const std::uint64_t base = mRecords.size();
      if (base + offsets[n] > std::numeric_limits<std::uint32_t>::max()) {
        return StageStatus::kIndexOverflow;
      }
      mRecords.resize(static_cast<std::size_t>(base + offsets[n]));
      const auto view = mRecords.view();
      forEach(n, [&](std::uint32_t i) {
        const auto at = static_cast<std::uint32_t>(base + offsets[i]);
        const std::uint32_t track = trackOf(i);
        assert(track < mRanges.size());
        mRanges[track] = TrackRange{.first = at, .count = offsets[i + 1U] - offsets[i]};
        fill(i, view, at);
      });
      return StageStatus::kComplete;
    });
  }

 private:
  std::pmr::memory_resource* mResource = std::pmr::get_default_resource();
  BoundedVector<TrackRange> mRanges;
  Records mRecords;
};

/// The references of a track's inline rows (`TrackOf`'s `clusters`, the first
/// `count`): each its own detector's row, held as its layer's. A device
/// stage fills the same way.
template <std::size_t K>
struct LayerRowReferences {
  const ClusterList<K>* clusters = nullptr;
  const std::uint8_t* count = nullptr;

  [[nodiscard]] GPUhdi() std::uint32_t references(std::uint32_t track) const { return count[track]; }

  template <typename View>
  GPUhdi() void write(std::uint32_t track, const View& out, std::uint32_t at) const
  {
    for (std::uint32_t k = 0; k < count[track]; ++k) {
      out.row[at + k] = clusters[track].index[k];
      out.source[at + k] = kOwnRows;
      out.role[at + k] = static_cast<std::uint8_t>(ReferenceRole::kLayerRow);
    }
  }
};

/// The timeframe's committed tracks, \a Det's (see the file's comment): the
/// tracker commits each iteration's tracks once its follow and length cut are
/// done, and again any that a later step changed in place.
template <DetectorLayout Det>
class TrackStore
{
 public:
  static constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  using References = TrackTable<TrackReference::template RecordT>;
  using Tracks = typename TrackRefit<Det>::Tracks;

  TrackStore() { mSources.emplace_back("rows"); }

  /// Every table's storage on \a pool, the tables that will be added too.
  void setMemoryPool(std::shared_ptr<BoundedMemoryResource> pool)
  {
    mMemoryPool = std::move(pool);
    clearStorage();
  }

  /// Start a timeframe: no committed tracks; the tables and sources stay.
  void clear()
  {
    mTracks = 0U;
    mReferences.clear();
    for (auto& table : mTables) {
      table.second->clear();
    }
  }

  /// Everything the timeframe grew back to the pool.
  void clearStorage() noexcept
  {
    mTracks = 0U;
    mReferences.setMemoryResource(resource());
    for (auto& table : mTables) {
      table.second->setMemoryResource(resource());
    }
  }

  /// The committed tracks: those of the store's TrackRefit, from 0.
  [[nodiscard]] std::size_t size() const { return mTracks; }
  [[nodiscard]] const References& references() const { return mReferences; }
  [[nodiscard]] References& references() { return mReferences; }

  /// A table of \a S's records per committed track, empty for those committed
  /// so far; std::logic_error if the store has one of \a S already.
  template <template <template <class> class> class S>
  TrackTable<S>& addTable()
  {
    if (findTable<S>() != nullptr) {
      throw std::logic_error("a track store holds one table per record type");
    }
    auto table = std::make_unique<TrackTable<S>>();
    table->setMemoryResource(resource());
    table->resizeTracks(mTracks);
    auto& added = *table;
    mTables.emplace_back(&detail::kTableKey<S>, std::move(table));
    return added;
  }

  /// The table of \a S's records, nullptr without one.
  template <template <template <class> class> class S>
  [[nodiscard]] TrackTable<S>* findTable()
  {
    for (auto& table : mTables) {
      if (table.first == &detail::kTableKey<S>) {
        return static_cast<TrackTable<S>*>(table.second.get());
      }
    }
    return nullptr;
  }

  template <template <template <class> class> class S>
  [[nodiscard]] const TrackTable<S>* findTable() const
  {
    return const_cast<TrackStore*>(this)->template findTable<S>();
  }

  /// A source of measurements that references name, \a name for reports:
  /// its id. std::length_error past the byte's 255 sources.
  SourceID addSource(std::string name)
  {
    if (mSources.size() > std::numeric_limits<SourceID>::max()) {
      throw std::length_error("a track store names at most 256 sources");
    }
    mSources.push_back(std::move(name));
    return static_cast<SourceID>(mSources.size() - 1U);
  }

  [[nodiscard]] std::size_t sourceCount() const { return mSources.size(); }
  [[nodiscard]] std::string_view sourceName(SourceID source) const { return mSources.at(source); }

  /// \a tracks from the first not committed yet to their end, committed with
  /// their rows (`LayerRowReferences`); every table grows to hold them, with
  /// empty ranges in the others. \a forEach as `TrackTable::append`'s.
  template <typename ForEach>
  [[nodiscard]] StageStatus commit(const ForEach& forEach, std::uint32_t chunk, const Tracks& tracks)
  {
    assert(tracks.size() >= mTracks);
    const auto begin = static_cast<std::uint32_t>(mTracks);
    const auto n = static_cast<std::uint32_t>(tracks.size() - mTracks);
    const auto status = guardOutOfMemory([&] {
      grow(tracks.size());
      return StageStatus::kComplete;
    });
    if (status != StageStatus::kComplete) {
      return status;
    }
    const auto view = tracks.view();
    const LayerRowReferences<kMaxClusters> rows{.clusters = view.clusters, .count = view.count};
    return mReferences.append(
      forEach, chunk, n, [begin](std::uint32_t i) { return begin + i; }, [&](std::uint32_t i) { return rows.references(begin + i); }, [&](std::uint32_t i, const auto& out, std::uint32_t at) { rows.write(begin + i, out, at); });
  }

  /// The committed tracks whose rows in \a tracks are no longer their
  /// references -- changed in place since their commit -- committed again.
  /// \a changed is the scratch for their indices.
  template <typename ForEach>
  [[nodiscard]] StageStatus recommitChanged(const ForEach& forEach, std::uint32_t chunk, const Tracks& tracks, BoundedVector<std::uint32_t>& changed)
  {
    assert(tracks.size() >= mTracks);
    const auto view = tracks.view();
    const auto status = guardOutOfMemory([&] {
      host_algorithms::compactIndices(forEach, static_cast<std::uint32_t>(mTracks), chunk, [&](std::uint32_t track) { return !holdsRows(view, track); }, changed);
      return StageStatus::kComplete;
    });
    if (status != StageStatus::kComplete || changed.empty()) {
      return status;
    }
    const LayerRowReferences<kMaxClusters> rows{.clusters = view.clusters, .count = view.count};
    return mReferences.append(
      forEach, chunk, static_cast<std::uint32_t>(changed.size()), [&](std::uint32_t i) { return changed[i]; }, [&](std::uint32_t i) { return rows.references(changed[i]); }, [&](std::uint32_t i, const auto& out, std::uint32_t at) { rows.write(changed[i], out, at); });
  }

  /// Whether every committed track's references are its rows in \a tracks,
  /// in order, each its own detector's and held as its layer's.
  [[nodiscard]] bool holdsRowsOf(const Tracks& tracks) const
  {
    if (tracks.size() != mTracks || mReferences.trackCount() != mTracks) {
      return false;
    }
    const auto view = tracks.view();
    for (std::size_t track = 0; track < mTracks; ++track) {
      if (!holdsRows(view, track)) {
        return false;
      }
    }
    return true;
  }

 private:
  /// The pool, or the default resource before one is set.
  [[nodiscard]] std::pmr::memory_resource* resource() const noexcept { return mMemoryPool ? static_cast<std::pmr::memory_resource*>(mMemoryPool.get()) : std::pmr::get_default_resource(); }

  /// Whether \a track's references are its rows in \a view.
  template <typename View>
  [[nodiscard]] bool holdsRows(const View& view, std::size_t track) const
  {
    const auto range = mReferences.range(track);
    if (range.count != view.count[track]) {
      return false;
    }
    const auto references = mReferences.view().records;
    for (std::uint32_t k = 0; k < range.count; ++k) {
      const auto at = range.first + k;
      if (references.row[at] != view.clusters[track].index[k] || references.source[at] != kOwnRows || references.role[at] != static_cast<std::uint8_t>(ReferenceRole::kLayerRow)) {
        return false;
      }
    }
    return true;
  }

  void grow(std::size_t tracks)
  {
    mReferences.resizeTracks(tracks);
    for (auto& table : mTables) {
      table.second->resizeTracks(tracks);
    }
    mTracks = tracks;
  }

  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  std::size_t mTracks = 0U;
  References mReferences;
  std::vector<std::pair<const void*, std::unique_ptr<detail::TrackTableBase>>> mTables;
  std::vector<std::string> mSources;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRACKSTORE_H_
