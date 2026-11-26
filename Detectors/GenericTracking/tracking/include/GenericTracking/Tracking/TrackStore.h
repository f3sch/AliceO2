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
/// The store holds each committed track's own fields (`CommittedTrack`: its
/// states, time, chi2, iteration) and one table of its own, the tracks'
/// *references*: per reference the row it names, the source whose
/// measurements that row indexes (`kOwnRows` for the tracker's frame, or a
/// source the store names with `addSource`) and its role (`ReferenceRole`).
/// Any other component may add a table of its own record type (`addTable`): a
/// matching step's matches to its detector's measurements, say, keyed by the
/// committed track. The tracker does not know those types; it grows every
/// table with the tracks it commits.
///
/// The tracker builds and extends its tracks with their rows inline (`TrackOf`,
/// one per layer at most) and commits them here; what reads tracks after the
/// tracking reads the store. `flatten` and `assign` take it to and from the
/// tracking device's products; `inlineTracks` gives back the tracker's layout
/// for what compares or digests in it.
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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <memory_resource>
#include <new>
#include <numeric>
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
  // As for MeasurementT: an SoA element proxy holds references.
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
  template <template <class> class C>
  struct RecordT {
    MEMLAYOUT_APPLY(RecordT, row, source, role)
    C<std::uint32_t> row; ///< the row in its source's measurements
    C<SourceID> source;
    C<std::uint8_t> role; ///< ReferenceRole bits
  };
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
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
  virtual ~TrackTableBase() = default;

  /// Every allocation reseated on \a resource, what was held released.
  virtual void setMemoryResource(std::pmr::memory_resource* resource) = 0;
  /// No tracks and no records; the storage kept.
  virtual void clear() = 0;
  /// Ranges for \a tracks tracks, the new ones empty.
  virtual void resizeTracks(std::uint32_t tracks) = 0;

 protected:
  // A table moves as its own type only, never through the base.
  TrackTableBase() = default;
  TrackTableBase(const TrackTableBase&) = default;
  TrackTableBase& operator=(const TrackTableBase&) = default;
  TrackTableBase(TrackTableBase&&) = default;
  TrackTableBase& operator=(TrackTableBase&&) = default;
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

  void resizeTracks(std::uint32_t tracks) override { mRanges.resize(tracks, TrackRange{}); }

  [[nodiscard]] std::uint32_t trackCount() const { return static_cast<std::uint32_t>(mRanges.size()); }
  /// Every record held, those of superseded ranges too.
  [[nodiscard]] std::uint32_t recordCount() const { return mRecords.size(); }
  [[nodiscard]] TrackRange range(std::uint32_t track) const { return mRanges[track]; }
  /// Record \a k of \a track's range.
  [[nodiscard]] Record read(std::uint32_t track, std::uint32_t k) const
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
      mRecords.resize(static_cast<std::uint32_t>(base + offsets[n]));
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

/// A committed track's own fields: `TrackOf`'s without its rows and their
/// count, which are its references and their range.
struct CommittedTrack {
  // As for MeasurementT: an SoA element proxy holds references.
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
  template <template <class> class C>
  struct RecordT {
    MEMLAYOUT_APPLY(RecordT, state, outwardState, timePS, widthPS, chi2, outwardChi2, ndf, status, iteration, addedOutward, addedInward, kinks, kinkCluster)
    C<BoundTrackState> state;
    C<BoundTrackState> outwardState;
    C<std::uint32_t> timePS;
    C<std::uint32_t> widthPS;
    C<Scalar> chi2;        ///< the inward pass's
    C<Scalar> outwardChi2; ///< the outward pass's, which `outwardState` ends
    C<std::int16_t> ndf;
    C<PropagationStatus> status; ///< the refit's: kOk, or why it failed
    C<std::uint8_t> iteration;   ///< the iteration that found it
    C<std::uint8_t> addedOutward;
    C<std::uint8_t> addedInward;
    C<std::uint8_t> kinks;       ///< the inward pass's
    C<std::uint8_t> kinkCluster; ///< its last kink's position among its layer rows; kNoKink without
  };
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
  static constexpr std::uint8_t kNoKink = 0xFFU;
};

/// Track \a from's fields of \a tracks (a `TrackOf` view) into track \a to
/// of \a committed (a `CommittedTrack` view), or back with \a Back.
template <bool Back = false, typename TrackView, typename CommittedView>
GPUhdi() void copyTrackFields(const TrackView& tracks, std::uint32_t from, const CommittedView& committed, std::uint32_t to)
{
  const auto copy = [](auto&& field, auto&& trackField) {
    if constexpr (Back) {
      trackField = field;
    } else {
      field = trackField;
    }
  };
  copy(committed.state[to], tracks.state[from]);
  copy(committed.outwardState[to], tracks.outwardState[from]);
  copy(committed.timePS[to], tracks.timePS[from]);
  copy(committed.widthPS[to], tracks.widthPS[from]);
  copy(committed.chi2[to], tracks.chi2[from]);
  copy(committed.outwardChi2[to], tracks.outwardChi2[from]);
  copy(committed.ndf[to], tracks.ndf[from]);
  copy(committed.status[to], tracks.status[from]);
  copy(committed.iteration[to], tracks.iteration[from]);
  copy(committed.addedOutward[to], tracks.addedOutward[from]);
  copy(committed.addedInward[to], tracks.addedInward[from]);
  copy(committed.kinks[to], tracks.kinks[from]);
  copy(committed.kinkCluster[to], tracks.kinkCluster[from]);
}

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

/// A track's rows \a rows (\a count of them) and its overlap partners
/// \a partners to \a end (OverlapPartners.h), each partner's `layerRow` its
/// row's position, written to \a out from \a at in path order: for each row,
/// its partners before it, the row, its partners after it.
template <std::size_t K, typename Partner, typename View>
GPUhdi() void writePartneredRows(const ClusterList<K>& rows, std::uint32_t count, const Partner* partners, const Partner* end, const View& out, std::uint32_t at)
{
  const auto emit = [&](std::uint32_t row, ReferenceRole role) {
    out.row[at] = row;
    out.source[at] = kOwnRows;
    out.role[at] = static_cast<std::uint8_t>(role);
    ++at;
  };
  for (std::uint32_t k = 0; k < count; ++k) {
    for (; partners != end && partners->layerRow == k && partners->ahead == 0U; ++partners) {
      emit(partners->row, ReferenceRole::kOverlapPartner);
    }
    emit(rows.index[k], ReferenceRole::kLayerRow);
    for (; partners != end && partners->layerRow == k; ++partners) {
      emit(partners->row, ReferenceRole::kOverlapPartner);
    }
  }
}

/// The references of a track's inline rows and of its overlap partners
/// (`writePartneredRows`). \a offsets and \a partners are per track from
/// \a first, the first of the tracks committed with them.
template <std::size_t K, typename Partner>
struct PartneredReferences {
  LayerRowReferences<K> rows;
  std::uint32_t first = 0U;
  const std::uint32_t* offsets = nullptr;
  const Partner* partners = nullptr;

  [[nodiscard]] GPUhdi() std::uint32_t references(std::uint32_t track) const { return rows.references(track) + offsets[track - first + 1U] - offsets[track - first]; }

  template <typename View>
  GPUhdi() void write(std::uint32_t track, const View& out, std::uint32_t at) const
  {
    writePartneredRows(rows.clusters[track], rows.count[track], partners + offsets[track - first], partners + offsets[track - first + 1U], out, at);
  }
};

/// The timeframe's committed tracks, \a Det's (see the file's comment): the
/// tracker commits each iteration's tracks once its follow and length cut are
/// done, and again any that a later step changed in place.
template <DetectorLayout Det>
class TrackStore
{
 public:
  static constexpr std::uint32_t kMaxClusters = maxChainLayers<Det>;
  using References = TrackTable<TrackReference::template RecordT>;
  /// The tracker's tracks as it builds them, rows inline (`TrackOf`).
  using Tracks = typename TrackRefit<Det>::Tracks;
  /// The committed tracks' own fields.
  using TrackColumns = BoundedColumns<CommittedTrack::template RecordT>;
  using TrackRecord = typename TrackColumns::Record;
  using ReferenceRecord = typename References::Record;

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
    mTrackColumns.resize(0U);
    mReferences.clear();
    for (auto& table : mTables) {
      table.second->clear();
    }
  }

  /// Everything the timeframe grew back to the pool.
  void clearStorage() noexcept
  {
    mTracks = 0U;
    mTrackColumns.setMemoryResource(resource());
    mReferences.setMemoryResource(resource());
    for (auto& table : mTables) {
      table.second->setMemoryResource(resource());
    }
  }

  /// The committed tracks, in commit order: every iteration's after the
  /// ones before, each iteration's best first.
  [[nodiscard]] std::uint32_t size() const { return mTracks; }
  [[nodiscard]] const TrackColumns& tracks() const { return mTrackColumns; }
  /// \a track's rows in its own detector's frame, in path order: its
  /// references' rows -- its layers' and its overlap partners -- every
  /// reference being one of the tracker's rows.
  [[nodiscard]] std::span<const std::uint32_t> rows(std::uint32_t track) const
  {
    const auto range = mReferences.range(track);
    return {mReferences.view().records.row + range.first, range.count};
  }
  [[nodiscard]] const References& references() const { return mReferences; }
  [[nodiscard]] References& references() { return mReferences; }
  /// \a track's layer rows alone, one per layer, in path order, into \a out
  /// (a chain's length): how many.
  [[nodiscard]] std::uint32_t layerRows(std::uint32_t track, std::span<std::uint32_t, kMaxClusters> out) const
  {
    const auto range = mReferences.range(track);
    const auto records = mReferences.view().records;
    std::uint32_t count = 0U;
    for (std::uint32_t at = range.first; at < range.first + range.count; ++at) {
      if ((records.role[at] & static_cast<std::uint8_t>(ReferenceRole::kLayerRow)) != 0U) {
        out[count++] = records.row[at];
      }
    }
    return count;
  }

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
    return const_cast<TrackStore*>(this)->template findTable<S>(); // NOLINT(cppcoreguidelines-pro-type-const-cast): the same lookup, read only
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

  [[nodiscard]] std::uint32_t sourceCount() const { return static_cast<std::uint32_t>(mSources.size()); }
  [[nodiscard]] std::string_view sourceName(SourceID source) const { return mSources.at(source); }

  /// \a tracks from the first not committed yet to their end, committed with
  /// their rows (`LayerRowReferences`); every table grows to hold them, with
  /// empty ranges in the others. \a forEach as `TrackTable::append`'s.
  template <typename ForEach>
  [[nodiscard]] StageStatus commit(const ForEach& forEach, std::uint32_t chunk, const Tracks& tracks)
  {
    const auto view = tracks.view();
    return commitWith(forEach, chunk, tracks, LayerRowReferences<kMaxClusters>{.clusters = view.clusters, .count = view.count});
  }

  /// The same, with their overlap partners (OverlapPartners.h): per track
  /// from the first not committed yet, \a offsets where its \a partners
  /// begin, one more than the tracks; each track's in path order.
  template <typename ForEach, typename Partner>
  [[nodiscard]] StageStatus commit(const ForEach& forEach, std::uint32_t chunk, const Tracks& tracks, std::span<const std::uint32_t> offsets, std::span<const Partner> partners)
  {
    assert(offsets.size() == tracks.size() - mTracks + 1U);
    const auto view = tracks.view();
    return commitWith(forEach, chunk, tracks, PartneredReferences<kMaxClusters, Partner>{.rows = {.clusters = view.clusters, .count = view.count}, .first = mTracks, .offsets = offsets.data(), .partners = partners.data()});
  }

  /// The committed tracks as \a tracks holds them after a step changed them
  /// in place (the follow of earlier iterations' tracks): every one's fields
  /// again, and new references for those whose rows are no longer theirs.
  /// \a changed is the scratch for their indices.
  template <typename ForEach>
  [[nodiscard]] StageStatus recommitChanged(const ForEach& forEach, std::uint32_t chunk, const Tracks& tracks, BoundedVector<std::uint32_t>& changed)
  {
    assert(tracks.size() >= mTracks);
    const auto view = tracks.view();
    const auto committed = mTrackColumns.view();
    forEach(mTracks, [&](std::uint32_t track) { copyTrackFields(view, track, committed, track); });
    const auto status = guardOutOfMemory([&] {
      host_algorithms::compactIndices(forEach, mTracks, chunk, [&](std::uint32_t track) { return !holdsRows(view, track); }, changed);
      return StageStatus::kComplete;
    });
    if (status != StageStatus::kComplete || changed.empty()) {
      return status;
    }
    const LayerRowReferences<kMaxClusters> rows{.clusters = view.clusters, .count = view.count};
    return mReferences.append(
      forEach, chunk, static_cast<std::uint32_t>(changed.size()), [&](std::uint32_t i) { return changed[i]; }, [&](std::uint32_t i) { return rows.references(changed[i]); }, [&](std::uint32_t i, const auto& out, std::uint32_t at) { rows.write(changed[i], out, at); });
  }

  /// The committed tracks \a listed (ascending) as \a tracks holds them,
  /// with their overlap partners (OverlapPartners.h): each one's fields
  /// again, and its rows and partners as new references -- \a listed[i]'s
  /// partners from \a offsets[i] to the next, in path order.
  template <typename ForEach, typename Partner>
  [[nodiscard]] StageStatus recommitPartnered(const ForEach& forEach, std::uint32_t chunk, const Tracks& tracks, std::span<const std::uint32_t> listed, std::span<const std::uint32_t> offsets, std::span<const Partner> partners)
  {
    assert(offsets.size() == listed.size() + 1U);
    const auto view = tracks.view();
    const auto committed = mTrackColumns.view();
    const auto n = static_cast<std::uint32_t>(listed.size());
    forEach(n, [&](std::uint32_t i) { copyTrackFields(view, listed[i], committed, listed[i]); });
    return mReferences.append(
      forEach, chunk, n, [&](std::uint32_t i) { return listed[i]; }, [&](std::uint32_t i) { return view.count[listed[i]] + offsets[i + 1U] - offsets[i]; },
      [&](std::uint32_t i, const auto& out, std::uint32_t at) { writePartneredRows(view.clusters[listed[i]], view.count[listed[i]], partners.data() + offsets[i], partners.data() + offsets[i + 1U], out, at); });
  }

  /// The committed tracks of \a matched (by track, then along it; a
  /// `track` and a `row` each) with their matches after their references
  /// (`ReferenceRole::kMatch`, InactiveMatch.h): each track's references
  /// again, as the store holds them, then its matched rows.
  template <typename ForEach, typename Matched>
  [[nodiscard]] StageStatus recommitMatched(const ForEach& forEach, std::uint32_t chunk, std::span<const Matched> matched)
  {
    // Each matched track's first entry in \a matched.
    std::vector<std::uint32_t> firsts;
    for (std::uint32_t i = 0; i < matched.size(); ++i) {
      if (i == 0U || matched[i].track != matched[i - 1U].track) {
        firsts.push_back(i);
      }
    }
    firsts.push_back(static_cast<std::uint32_t>(matched.size()));
    const auto n = static_cast<std::uint32_t>(firsts.size() - 1U);
    std::vector<TrackRange> before(n);
    for (std::uint32_t i = 0; i < n; ++i) {
      before[i] = mReferences.range(matched[firsts[i]].track);
    }
    // The fill reads the old references through \a out: the append resizes
    // the columns first, and their earlier records stay where they were.
    return mReferences.append(
      forEach, chunk, n, [&](std::uint32_t i) { return matched[firsts[i]].track; }, [&](std::uint32_t i) { return before[i].count + (firsts[i + 1U] - firsts[i]); },
      [&](std::uint32_t i, const auto& out, std::uint32_t at) {
        for (std::uint32_t k = 0; k < before[i].count; ++k) {
          out.row[at + k] = out.row[before[i].first + k];
          out.source[at + k] = out.source[before[i].first + k];
          out.role[at + k] = out.role[before[i].first + k];
        }
        at += before[i].count;
        for (std::uint32_t m = firsts[i]; m < firsts[i + 1U]; ++m) {
          out.row[at] = matched[m].row;
          out.source[at] = kOwnRows;
          out.role[at] = static_cast<std::uint8_t>(ReferenceRole::kMatch);
          ++at;
        }
      });
  }

  /// Whether the committed tracks are \a tracks: the same fields, bit for
  /// bit, and as the references held as their layers' its rows in order,
  /// each its own detector's.
  [[nodiscard]] bool holds(const Tracks& tracks) const
  {
    if (tracks.size() != mTracks || mReferences.trackCount() != mTracks) {
      return false;
    }
    const auto view = tracks.view();
    for (std::uint32_t track = 0; track < mTracks; ++track) {
      if (!holdsRows(view, track) || !sameFields(view, track)) {
        return false;
      }
    }
    return true;
  }

  /// The committed tracks as flat products: \a records their fields,
  /// \a offsets (one more than the tracks) where each one's references begin
  /// in \a references. Each vector is cleared first. The products' order is
  /// by time (the box's start), then the more clusters first (references to
  /// the tracker's own rows that are not outliers), then the smaller inward
  /// chi2 per degree of freedom; equal ones in the order committed. The
  /// store keeps the order committed, which the tracker's own tables follow.
  template <class Records, class Offsets, class ReferenceVector>
  void flatten(Records& records, Offsets& offsets, ReferenceVector& references) const
  {
    records.clear();
    offsets.clear();
    references.clear();
    const auto columns = mTrackColumns.view();
    std::vector<std::uint32_t> clusters(mTracks, 0U);
    std::vector<Scalar> chi2PerDof(mTracks);
    for (std::uint32_t track = 0; track < mTracks; ++track) {
      const auto range = mReferences.range(track);
      for (std::uint32_t k = 0; k < range.count; ++k) {
        const auto reference = mReferences.read(track, k);
        const bool outlier = (reference.role & static_cast<std::uint8_t>(ReferenceRole::kOutlier)) != 0U;
        clusters[track] += reference.source == kOwnRows && !outlier ? 1U : 0U;
      }
      // No degree of freedom, or a refit that left no finite chi2: last.
      const auto ratio = columns.ndf[track] > 0 ? columns.chi2[track] / static_cast<Scalar>(columns.ndf[track]) : std::numeric_limits<Scalar>::infinity();
      chi2PerDof[track] = std::isfinite(ratio) ? ratio : std::numeric_limits<Scalar>::infinity();
    }
    std::vector<std::uint32_t> order(mTracks);
    std::iota(order.begin(), order.end(), 0U);
    std::ranges::stable_sort(order, [&](std::uint32_t a, std::uint32_t b) {
      if (columns.timePS[a] != columns.timePS[b]) {
        return columns.timePS[a] < columns.timePS[b];
      }
      if (clusters[a] != clusters[b]) {
        return clusters[a] > clusters[b];
      }
      return chi2PerDof[a] < chi2PerDof[b];
    });
    records.reserve(mTracks);
    offsets.reserve(mTracks + 1U);
    offsets.push_back(0U);
    for (const auto track : order) {
      records.push_back(mTrackColumns.read(track));
      const auto range = mReferences.range(track);
      for (std::uint32_t k = 0; k < range.count; ++k) {
        references.push_back(mReferences.read(track, k));
      }
      offsets.push_back(static_cast<std::uint32_t>(references.size()));
    }
  }

  /// The store holding flat products (`flatten`'s; all three empty for no
  /// tracks) and nothing else; the added tables stay, empty for every track. std::invalid_argument where
  /// the offsets do not fit the records and references.
  void assign(std::span<const TrackRecord> records, std::span<const std::uint32_t> offsets, std::span<const ReferenceRecord> references)
  {
    if (offsets.empty() && records.empty() && references.empty()) {
      clear();
      return;
    }
    if (offsets.size() != records.size() + 1U || offsets.front() != 0U || offsets.back() != references.size() || !std::ranges::is_sorted(offsets)) {
      throw std::invalid_argument("track products whose offsets do not fit their tracks and references");
    }
    clear();
    grow(static_cast<std::uint32_t>(records.size()));
    for (std::uint32_t track = 0; track < records.size(); ++track) {
      mTrackColumns.write(track, records[track]);
    }
    const auto status = mReferences.append(
      [](std::uint32_t n, const auto& body) {
        for (std::uint32_t i = 0; i < n; ++i) {
          body(i);
        }
      },
      1024U, static_cast<std::uint32_t>(records.size()), [](std::uint32_t i) { return i; }, [&](std::uint32_t i) { return offsets[i + 1U] - offsets[i]; },
      [&](std::uint32_t i, const auto& out, std::uint32_t at) {
        for (std::uint32_t k = offsets[i]; k < offsets[i + 1U]; ++k) {
          out.row[at + k - offsets[i]] = references[k].row;
          out.source[at + k - offsets[i]] = references[k].source;
          out.role[at + k - offsets[i]] = references[k].role;
        }
      });
    if (status != StageStatus::kComplete) {
      throw std::bad_alloc();
    }
  }

 private:
  /// The pool, or the default resource before one is set.
  [[nodiscard]] std::pmr::memory_resource* resource() const noexcept { return mMemoryPool ? static_cast<std::pmr::memory_resource*>(mMemoryPool.get()) : std::pmr::get_default_resource(); }

  /// Whether \a track's references held as its layers' are its rows in
  /// \a view, in order, its own detector's.
  template <typename View>
  [[nodiscard]] bool holdsRows(const View& view, std::uint32_t track) const
  {
    const auto range = mReferences.range(track);
    const auto references = mReferences.view().records;
    std::uint32_t k = 0;
    for (std::uint32_t at = range.first; at < range.first + range.count; ++at) {
      if ((references.role[at] & static_cast<std::uint8_t>(ReferenceRole::kLayerRow)) == 0U) {
        continue;
      }
      if (k == view.count[track] || references.row[at] != view.clusters[track].index[k] || references.source[at] != kOwnRows) {
        return false;
      }
      ++k;
    }
    return k == view.count[track];
  }

  /// Whether \a track's fields are those of \a view's track, bit for bit.
  template <typename View>
  [[nodiscard]] bool sameFields(const View& view, std::uint32_t track) const
  {
    const auto mine = mTrackColumns.read(track);
    TrackRecord theirs{};
    const typename TrackColumns::View one{CommittedTrack::RecordT<memlayout::pointer>{.state = &theirs.state, .outwardState = &theirs.outwardState, .timePS = &theirs.timePS, .widthPS = &theirs.widthPS, .chi2 = &theirs.chi2, .outwardChi2 = &theirs.outwardChi2, .ndf = &theirs.ndf, .status = &theirs.status, .iteration = &theirs.iteration, .addedOutward = &theirs.addedOutward, .addedInward = &theirs.addedInward, .kinks = &theirs.kinks, .kinkCluster = &theirs.kinkCluster}};
    copyTrackFields(view, track, one, 0U);
    bool same = true;
    // Bit for bit is the comparison asked for: a column's bytes as stored.
    // NOLINTNEXTLINE(bugprone-suspicious-memory-comparison)
    mine.apply([&](const auto&... a) { theirs.apply([&](const auto&... b) { ((same = same && std::memcmp(&a, &b, sizeof(a)) == 0), ...); }); });
    return same;
  }

  /// `commit` with \a references, which counts and writes each track's.
  template <typename ForEach, typename References>
  [[nodiscard]] StageStatus commitWith(const ForEach& forEach, std::uint32_t chunk, const Tracks& tracks, const References& references)
  {
    assert(tracks.size() >= mTracks);
    const auto begin = mTracks;
    const auto n = tracks.size() - mTracks;
    const auto status = guardOutOfMemory([&] {
      grow(tracks.size());
      return StageStatus::kComplete;
    });
    if (status != StageStatus::kComplete) {
      return status;
    }
    const auto view = tracks.view();
    const auto committed = mTrackColumns.view();
    forEach(n, [&](std::uint32_t i) { copyTrackFields(view, begin + i, committed, begin + i); });
    return mReferences.append(
      forEach, chunk, n, [begin](std::uint32_t i) { return begin + i; }, [&](std::uint32_t i) { return references.references(begin + i); }, [&](std::uint32_t i, const auto& out, std::uint32_t at) { references.write(begin + i, out, at); });
  }

  void grow(std::uint32_t tracks)
  {
    mTrackColumns.resize(tracks);
    mReferences.resizeTracks(tracks);
    for (auto& table : mTables) {
      table.second->resizeTracks(tracks);
    }
    mTracks = tracks;
  }

  std::shared_ptr<BoundedMemoryResource> mMemoryPool;
  std::uint32_t mTracks = 0U;
  TrackColumns mTrackColumns;
  References mReferences;
  std::vector<std::pair<const void*, std::unique_ptr<detail::TrackTableBase>>> mTables;
  std::vector<std::string> mSources;
};

/// The time box of a track's \a rows (`trackTimeBox`'s, from its committed
/// rows).
template <typename Columns>
[[nodiscard]] TimeEstPS trackTimeBox(const Columns& columns, std::span<const std::uint32_t> rows)
{
  auto box = timeBoxAt(columns, rows[0]);
  for (std::uint32_t i = 1; i < rows.size(); ++i) {
    box += timeBoxAt(columns, rows[i]);
  }
  return box;
}

/// \a store's tracks as the tracker builds them (`TrackOf`), on \a resource:
/// the fields, each track's layer rows inline and the slots after them
/// `kNoCluster`; its other references are not in the layout. For what reads that layout -- the digests, a comparison with
/// a device's tracks; std::length_error for a track with more rows than the
/// layout holds.
template <DetectorLayout Det>
[[nodiscard]] typename TrackStore<Det>::Tracks inlineTracks(const TrackStore<Det>& store, std::pmr::memory_resource* resource)
{
  constexpr std::uint32_t kClusters = TrackStore<Det>::kMaxClusters;
  typename TrackStore<Det>::Tracks tracks;
  tracks.setMemoryResource(resource);
  tracks.resize(store.size());
  const auto view = tracks.view();
  const auto committed = store.tracks().view();
  for (std::uint32_t track = 0; track < store.size(); ++track) {
    copyTrackFields<true>(view, track, committed, track);
    std::fill(std::begin(view.clusters[track].index), std::end(view.clusters[track].index), ClusterList<kClusters>::kNoCluster);
    const auto range = store.references().range(track);
    const auto references = store.references().view().records;
    std::uint8_t count = 0U;
    for (std::uint32_t at = range.first; at < range.first + range.count; ++at) {
      if ((references.role[at] & static_cast<std::uint8_t>(ReferenceRole::kLayerRow)) == 0U) {
        continue;
      }
      if (count == kClusters) {
        throw std::length_error("a committed track with more layer rows than a chain holds");
      }
      view.clusters[track].index[count++] = references.row[at];
    }
    view.count[track] = count;
  }
  return tracks;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRACKSTORE_H_
