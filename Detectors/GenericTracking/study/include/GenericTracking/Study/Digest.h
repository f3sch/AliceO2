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

/// \file Digest.h
/// \brief A 64-bit fingerprint of output columns, as their values' bits.
///
/// A change to how the host reaches the tracks may not change a single bit
/// of them, and comparing two runs' digests is how that is checked without
/// keeping either run's output.
///
/// Every value is hashed as the bits of its fields, never as the bytes of its
/// object, whose padding nobody writes. A column type is accepted only when its
/// object representation is its value (integers, `ClusterList`), when it is
/// floating point, or when it is decomposed here field by field; anything else
/// fails to compile rather than hash its padding.

#ifndef GTRK_STUDY_DIGEST_H_
#define GTRK_STUDY_DIGEST_H_

#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackRefit.h"
#include "GenericTracking/Tracking/TrackStore.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory_resource>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace o2::gtrk::study
{

/// A running 64-bit hash over words: stable across runs, builds and
/// machines of one byte order. Not cryptographic; a fingerprint.
class Digest
{
 public:
  /// One value, as the bits of its fields.
  template <typename T>
  void add(const T& value)
  {
    if constexpr (std::is_same_v<T, BoundTrackState>) {
      for (const auto parameter : value.par) {
        add(parameter);
      }
      for (const auto covariance : value.cov) {
        add(covariance);
      }
      add(value.surface);
      // 0 on every flat sensor: hashed only when set, so that a layout without
      // curved sensors digests without it.
      if (!exactlyEqual(value.surfacePhi, Scalar{0})) {
        add(value.surfacePhi);
      }
    } else if constexpr (std::is_same_v<T, VertexSeed>) {
      add(value.timePS);
      add(value.widthPS);
      add(value.z);
      add(value.varZ);
      add(value.varTransverse);
    } else {
      static_assert(std::is_floating_point_v<T> || std::has_unique_object_representations_v<T>, "a digest reads a value's bits; decompose types with padding field by field");
      bytes(&value, sizeof(T));
    }
  }

  /// \a count values at \a values, preceded by their count.
  template <typename T>
  void column(const T* values, std::size_t count)
  {
    mix(count);
    if constexpr (std::is_floating_point_v<T> || std::has_unique_object_representations_v<T>) {
      bytes(values, count * sizeof(T));
    } else {
      for (std::size_t index = 0; index < count; ++index) {
        add(values[index]);
      }
    }
  }

  template <typename T>
  void column(std::span<T> values)
  {
    column(values.data(), values.size());
  }

  /// Every column of a memlayout view, \a rows rows each, in schema order.
  template <typename View>
  void columns(const View& view, std::size_t rows)
  {
    view.apply([&](const auto*... columns) { (column(columns, rows), ...); });
  }

  /// The fingerprint so far, avalanched; the digest itself is unchanged.
  [[nodiscard]] std::uint64_t result() const
  {
    std::uint64_t hash = mState ^ mWords;
    hash ^= hash >> 33U;
    hash *= 0xff51afd7ed558ccdULL;
    hash ^= hash >> 33U;
    hash *= 0xc4ceb9fe1a85ec53ULL;
    hash ^= hash >> 33U;
    return hash;
  }

 private:
  /// The bytes \a data, \a size of them, in order.
  void bytes(const void* data, std::size_t size)
  {
    const auto* at = static_cast<const unsigned char*>(data);
    while (size >= sizeof(std::uint64_t)) {
      std::uint64_t word = 0U;
      std::memcpy(&word, at, sizeof(word));
      mix(word);
      at += sizeof(word);
      size -= sizeof(word);
    }
    if (size != 0U) {
      std::uint64_t word = 0U;
      std::memcpy(&word, at, size);
      mix(word ^ (std::uint64_t{size} << 56U));
    }
  }

  void mix(std::uint64_t word)
  {
    mState = std::rotl(mState ^ (word * 0x9e3779b97f4a7c15ULL), 31) * 0xbf58476d1ce4e5b9ULL;
    ++mWords;
  }

  std::uint64_t mState = 0x243f6a8885a308d3ULL;
  std::uint64_t mWords = 0U;
};

/// The counts of every iteration's summary of \a result (a tracker's
/// `TimeframeResult`): part of both digests. The counts of the earlier
/// iterations' follow and the length cut (`tracksFollowedEarlier`,
/// `tracksDropped`...) are left out; the tracks they change are in both.
template <typename Result>
void addSummaries(Digest& digest, const Result& result)
{
  for (const auto& summary : result.iterations) {
    for (const std::size_t count : {static_cast<std::size_t>(summary.status), summary.slices, summary.roads, summary.records, summary.duplicates, summary.mismatches, summary.rounds, summary.tracks, summary.refitFailures, summary.tracksFollowed, summary.clustersFollowedOutward, summary.clustersFollowedInward}) {
      digest.add(static_cast<std::uint64_t>(count));
    }
    digest.add(summary.haloPS);
  }
}

/// \a store's references other than its tracks' layer rows (overlap
/// partners, other sources) into \a digest, where it has any -- nothing
/// otherwise, so that a run without them digests as before: per such
/// reference its track, row, source and role, in track order, or with
/// \a inClusters as cluster indices of \a frame and sorted, with the track
/// named by its first row's cluster, which no track order changes.
template <DetectorLayout Det>
void addOtherReferences(Digest& digest, const TrackStore<Det>& store, const TimeFrame<Det>* inClusters = nullptr)
{
  const auto& references = store.references();
  const auto view = references.view().records;
  std::vector<std::array<std::uint32_t, 4>> others;
  for (std::size_t track = 0; track < store.size(); ++track) {
    const auto range = references.range(track);
    const std::size_t mine = others.size();
    std::uint32_t firstRow = 0U;
    bool haveFirst = false;
    for (std::uint32_t at = range.first; at < range.first + range.count; ++at) {
      if ((view.role[at] & static_cast<std::uint8_t>(ReferenceRole::kLayerRow)) != 0U && view.source[at] == kOwnRows) {
        if (!haveFirst) {
          firstRow = view.row[at];
          haveFirst = true;
        }
        continue;
      }
      others.push_back({static_cast<std::uint32_t>(track), view.row[at], view.source[at], view.role[at]});
    }
    if (inClusters != nullptr) {
      const auto columns = inClusters->measurements();
      for (std::size_t at = mine; at < others.size(); ++at) {
        others[at][0] = haveFirst ? columns.clusterIndex[firstRow] : std::numeric_limits<std::uint32_t>::max();
        others[at][1] = others[at][2] == kOwnRows ? columns.clusterIndex[others[at][1]] : others[at][1];
      }
    }
  }
  if (others.empty()) {
    return;
  }
  if (inClusters != nullptr) {
    std::ranges::sort(others);
  }
  digest.add(static_cast<std::uint64_t>(others.size()));
  for (const auto& other : others) {
    digest.add(other);
  }
}

/// The fingerprint of a complete run's output: every column of every track,
/// the used mask \a used the run left (the tracker's `used()`), and the counts
/// of every iteration's summary -- nothing that is a time, nor a size in
/// bytes, which a record's layout changes without changing the output. Two
/// runs that agree on it produced the same tracks bit for bit.
template <typename Result>
[[nodiscard]] std::uint64_t trackingDigest(const Result& result, std::span<const std::uint8_t> used)
{
  Digest digest;
  digest.add(static_cast<std::uint32_t>(result.status));
  if (result.store != nullptr) {
    // The tracker's layout, rows inline: the same bytes as before the store.
    const auto tracks = inlineTracks(*result.store, std::pmr::get_default_resource());
    digest.columns(tracks.view(), tracks.size());
    addOtherReferences(digest, *result.store);
  }
  digest.column(used);
  addSummaries(digest, result);
  return digest.result();
}

/// \a tracks in cluster numbers into \a digest: each track's rows as cluster
/// indices, the tracks in the order of their cluster lists, and every other
/// column as is. The part of `clusterDigest` that reads only the tracks, so
/// that a device's tracks digest alike.
template <DetectorLayout Det>
void addTracksInClusters(Digest& digest, const typename TrackRefit<Det>::Tracks& tracks, const TimeFrame<Det>& frame)
{
  constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  const auto columns = frame.measurements();
  const auto clusterOf = [&](std::uint32_t row) { return row == ClusterList<kMaxClusters>::kNoCluster ? row : columns.clusterIndex[row]; };
  const auto view = tracks.view();
  const std::size_t n = tracks.size();
  std::vector<ClusterList<kMaxClusters>> clusters(n);
  for (std::size_t track = 0; track < n; ++track) {
    for (std::size_t at = 0; at < kMaxClusters; ++at) {
      clusters[track].index[at] = at < view.count[track] ? clusterOf(view.clusters[track].index[at]) : ClusterList<kMaxClusters>::kNoCluster;
    }
  }
  std::vector<std::size_t> order(n);
  for (std::size_t track = 0; track < n; ++track) {
    order[track] = track;
  }
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    const auto& first = clusters[a].index;
    const auto& second = clusters[b].index;
    const auto differ = std::mismatch(std::begin(first), std::end(first), std::begin(second));
    if (differ.first != std::end(first)) {
      return *differ.first < *differ.second;
    }
    return view.iteration[a] < view.iteration[b];
  });
  digest.add(static_cast<std::uint64_t>(n));
  for (const auto track : order) {
    view.apply([&](const auto*... column) {
      (
        [&](const auto* values) {
          if constexpr (std::is_same_v<std::remove_cv_t<std::remove_pointer_t<decltype(values)>>, ClusterList<kMaxClusters>>) {
            digest.add(clusters[track]);
          } else {
            digest.add(values[track]);
          }
        }(column),
        ...);
    });
  }
}

/// The same output in the producer's cluster numbers instead of the frame's
/// rows. Rows are numbered in grid order, so another grid layout renumbers
/// them without changing a track; this does not see the renumbering. Every
/// column of every track, its clusters as cluster numbers, the tracks ordered
/// by those numbers (then by iteration); which cluster the used mask \a used
/// marks how; and every count of the iteration summaries. \a frame is the one
/// \a result was tracked on.
template <DetectorLayout Det, typename Result>
[[nodiscard]] std::uint64_t clusterDigest(const Result& result, const TimeFrame<Det>& frame, std::span<const std::uint8_t> used)
{
  const auto columns = frame.measurements();
  Digest digest;
  digest.add(static_cast<std::uint32_t>(result.status));
  if (result.store != nullptr) {
    addTracksInClusters<Det>(digest, inlineTracks(*result.store, std::pmr::get_default_resource()), frame);
    addOtherReferences(digest, *result.store, &frame);
  }
  std::uint32_t clusterEnd = 0U;
  for (std::uint32_t row = 0; row < frame.size(); ++row) {
    clusterEnd = std::max(clusterEnd, columns.clusterIndex[row] + 1U);
  }
  std::vector<std::uint8_t> usedByCluster(clusterEnd, std::uint8_t{0xFFU});
  for (std::uint32_t row = 0; row < frame.size(); ++row) {
    usedByCluster[columns.clusterIndex[row]] = used[row];
  }
  digest.column(std::span<const std::uint8_t>{usedByCluster});
  addSummaries(digest, result);
  return digest.result();
}

/// \a tracks' columns, as "digest, tracks" of a device's tracks.
template <class Tracks>
[[nodiscard]] std::uint64_t tracksDigest(const Tracks& tracks)
{
  Digest digest;
  digest.columns(std::as_const(tracks).view(), tracks.size());
  return digest.result();
}

/// \a store's tracks as "digest, tracks": the tracker's layout's columns
/// (`inlineTracks`).
template <DetectorLayout Det>
[[nodiscard]] std::uint64_t tracksDigest(const TrackStore<Det>& store)
{
  const auto tracks = inlineTracks(store, std::pmr::get_default_resource());
  Digest digest;
  digest.columns(tracks.view(), tracks.size());
  addOtherReferences(digest, store);
  return digest.result();
}

/// A vertex list and its contributor counts, as the vertexer gives them.
[[nodiscard]] inline std::uint64_t vertexDigest(std::span<const VertexSeed> vertices, std::span<const std::uint32_t> contributors)
{
  Digest digest;
  digest.add(vertices.size());
  for (const auto& vertex : vertices) {
    digest.add(vertex);
  }
  digest.column(contributors.data(), contributors.size());
  return digest.result();
}

/// \a records (the tracking device's product) as the tracker's columns, on
/// \a resource: what the digests of a device's tracks take.
template <DetectorLayout Det>
[[nodiscard]] typename TrackRefit<Det>::Tracks recordsTracks(std::span<const typename TrackRefit<Det>::Tracks::Record> records, std::pmr::memory_resource* resource)
{
  typename TrackRefit<Det>::Tracks tracks;
  tracks.setMemoryResource(resource);
  tracks.resize(records.size());
  for (std::size_t track = 0; track < records.size(); ++track) {
    tracks.write(track, records[track]);
  }
  return tracks;
}

/// \a tracks in cluster numbers alone: unchanged by a grid layout, and the
/// same for a device's tracks as for the host's where they are the same.
template <DetectorLayout Det>
[[nodiscard]] std::uint64_t tracksClusterDigest(const typename TrackRefit<Det>::Tracks& tracks, const TimeFrame<Det>& frame)
{
  Digest digest;
  addTracksInClusters<Det>(digest, tracks, frame);
  return digest.result();
}

template <DetectorLayout Det>
[[nodiscard]] std::uint64_t tracksClusterDigest(const TrackStore<Det>& store, const TimeFrame<Det>& frame)
{
  Digest digest;
  addTracksInClusters<Det>(digest, inlineTracks(store, std::pmr::get_default_resource()), frame);
  addOtherReferences(digest, store, &frame);
  return digest.result();
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DIGEST_H_
