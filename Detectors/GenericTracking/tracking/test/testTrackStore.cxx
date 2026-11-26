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

/// \file testTrackStore.cxx
/// \brief The committed tracks' store: each track's references are its rows,
///        a track changed in place is committed again, a table another
///        component adds grows with the tracks, and the timeframe's storage
///        goes back to the pool.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK TrackStore
#include "GenericTracking/Tracking/OverlapPartners.h"
#include "GenericTracking/Tracking/TrackStore.h"

#include "ToyBarrel8Fixture.h"

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <thread>
#include <vector>

namespace o2::gtrk
{

namespace
{

using Det = toybarrel8::Det;
using Store = TrackStore<Det>;
using Tracks = Store::Tracks;
constexpr std::uint32_t kClusters = Store::kMaxClusters;

/// A matching step's record, as another component would add: per committed
/// track, any number of hits of its detector with their chi2.
struct ToyMatch {
  // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
  template <template <class> class C>
  struct RecordT {
    MEMLAYOUT_APPLY(RecordT, hit, chi2)
    C<std::uint32_t> hit;
    C<float> chi2;
  };
  // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
};

/// body(i) for i below n, in order.
const auto serial = [](std::uint32_t n, const auto& body) {
  for (std::uint32_t i = 0; i < n; ++i) {
    body(i);
  }
};

/// body(i) for i below n, on four threads, each a contiguous quarter.
const auto threaded = [](std::uint32_t n, const auto& body) {
  constexpr std::uint32_t kThreads = 4U;
  std::vector<std::jthread> threads;
  threads.reserve(kThreads);
  for (std::uint32_t t = 0; t < kThreads; ++t) {
    threads.emplace_back([&, t] {
      for (std::uint32_t i = (n * t) / kThreads; i < (n * (t + 1U)) / kThreads; ++i) {
        body(i);
      }
    });
  }
};

/// \a counts.size() more tracks in \a tracks, track t's rows 100 t + k.
void appendTracks(Tracks& tracks, const std::vector<std::uint8_t>& counts)
{
  const auto begin = tracks.size();
  const auto added = static_cast<std::uint32_t>(counts.size());
  tracks.resize(begin + added);
  const auto view = tracks.view();
  for (std::uint32_t i = 0; i < added; ++i) {
    const auto track = begin + i;
    view.count[track] = counts[i];
    std::ranges::fill(view.clusters[track].index, ClusterList<kClusters>::kNoCluster);
    for (std::uint8_t k = 0; k < counts[i]; ++k) {
      view.clusters[track].index[k] = (100U * track) + k;
    }
  }
}

std::vector<std::uint32_t> referencedRows(const Store& store, std::uint32_t track)
{
  std::vector<std::uint32_t> rows;
  for (std::uint32_t k = 0; k < store.references().range(track).count; ++k) {
    rows.push_back(store.references().read(track, k).row);
  }
  return rows;
}

} // namespace

/// Two iterations' tracks committed one after the other: each track's range
/// holds its rows in order, its own detector's and as its layer's, the
/// ranges one after another and exactly as long as the rows.
BOOST_AUTO_TEST_CASE(EachCommittedTrackRefersToItsRows)
{
  auto pool = std::make_shared<BoundedMemoryResource>();
  Store store;
  store.setMemoryPool(pool);
  Tracks tracks;
  tracks.setMemoryResource(pool.get());
  appendTracks(tracks, {3U, 0U, static_cast<std::uint8_t>(kClusters)});
  BOOST_REQUIRE(store.commit(serial, 2U, tracks) == StageStatus::kComplete);
  appendTracks(tracks, {4U, 2U});
  BOOST_REQUIRE(store.commit(serial, 2U, tracks) == StageStatus::kComplete);
  BOOST_TEST(store.size() == 5U);
  BOOST_TEST(store.holds(tracks));
  BOOST_TEST(store.references().recordCount() == 3U + 0U + kClusters + 4U + 2U);
  std::uint32_t next = 0U;
  for (std::uint32_t track = 0; track < store.size(); ++track) {
    const auto range = store.references().range(track);
    BOOST_TEST(range.first == next);
    BOOST_TEST(range.count == tracks.view().count[track]);
    next += range.count;
    for (std::uint32_t k = 0; k < range.count; ++k) {
      const auto reference = store.references().read(track, k);
      BOOST_TEST(reference.row == (100U * track) + k);
      BOOST_TEST(reference.source == kOwnRows);
      BOOST_TEST(reference.role == static_cast<std::uint8_t>(ReferenceRole::kLayerRow));
    }
  }
}

/// The products in time order, then the more clusters, then the smaller
/// chi2 per degree of freedom, each track with its own references; the store
/// keeps the order committed.
BOOST_AUTO_TEST_CASE(TheProductsAreOrderedByTimeThenQuality)
{
  Store store;
  Tracks tracks;
  appendTracks(tracks, {3U, 5U, 5U, 4U});
  const auto view = tracks.view();
  const std::uint32_t timePS[] = {200U, 100U, 100U, 100U};
  const float chi2[] = {1.F, 9.F, 5.F, 0.F};
  const std::int16_t ndf[] = {1, 5, 5, 0};
  for (std::uint32_t track = 0; track < tracks.size(); ++track) {
    view.timePS[track] = timePS[track];
    view.chi2[track] = chi2[track];
    view.ndf[track] = ndf[track];
  }
  BOOST_REQUIRE(store.commit(serial, 2U, tracks) == StageStatus::kComplete);
  std::vector<Store::TrackRecord> records;
  std::vector<std::uint32_t> offsets;
  std::vector<Store::ReferenceRecord> references;
  store.flatten(records, offsets, references);
  const std::uint32_t expected[] = {2U, 1U, 3U, 0U};
  BOOST_REQUIRE(records.size() == 4U);
  for (std::size_t place = 0; place < records.size(); ++place) {
    const auto track = expected[place];
    BOOST_TEST_CONTEXT("place " << place)
    {
      BOOST_TEST(records[place].timePS == timePS[track]);
      BOOST_TEST(records[place].chi2 == chi2[track]);
      BOOST_TEST(offsets[place + 1U] - offsets[place] == tracks.view().count[track]);
      BOOST_TEST(references[offsets[place]].row == 100U * track);
    }
  }
  BOOST_TEST(store.holds(tracks));
}

/// The same tracks committed on four threads: the same ranges and records as
/// the serial commit.
BOOST_AUTO_TEST_CASE(AParallelCommitIsTheSerialOne)
{
  Tracks tracks;
  std::vector<std::uint8_t> counts;
  counts.reserve(101U);
  for (std::uint32_t track = 0; track < 101U; ++track) {
    counts.push_back(static_cast<std::uint8_t>((track * 7U) % (kClusters + 1U)));
  }
  appendTracks(tracks, counts);
  Store one;
  Store four;
  BOOST_REQUIRE(one.commit(serial, 8U, tracks) == StageStatus::kComplete);
  BOOST_REQUIRE(four.commit(threaded, 8U, tracks) == StageStatus::kComplete);
  BOOST_TEST(four.holds(tracks));
  for (std::uint32_t track = 0; track < tracks.size(); ++track) {
    BOOST_TEST(one.references().range(track).first == four.references().range(track).first);
    BOOST_TEST(one.references().range(track).count == four.references().range(track).count);
  }
}

/// A committed track whose rows change in place (as the follow of earlier
/// iterations' tracks extends or trims them) no longer matches; committed
/// again, it gets a new range after the records held and matches, and the
/// unchanged tracks keep theirs.
BOOST_AUTO_TEST_CASE(ATrackChangedInPlaceIsCommittedAgain)
{
  Store store;
  Tracks tracks;
  appendTracks(tracks, {3U, 5U, 2U});
  BOOST_REQUIRE(store.commit(serial, 2U, tracks) == StageStatus::kComplete);
  const auto before = store.references().recordCount();
  const auto view = tracks.view();
  view.clusters[0].index[3] = 7777U; // extended by a row
  view.count[0] = 4U;
  view.clusters[2].index[1] = 8888U; // its outer end replaced
  BOOST_TEST(!store.holds(tracks));
  BoundedVector<std::uint32_t> changed;
  BOOST_REQUIRE(store.recommitChanged(serial, 2U, tracks, changed) == StageStatus::kComplete);
  BOOST_TEST((std::vector<std::uint32_t>(changed.begin(), changed.end()) == std::vector<std::uint32_t>{0U, 2U}));
  BOOST_TEST(store.holds(tracks));
  BOOST_TEST(store.references().recordCount() == before + 4U + 2U);
  BOOST_TEST(store.references().range(0).first == before);
  BOOST_TEST(store.references().range(1).first == 3U);
  BOOST_TEST((referencedRows(store, 0) == std::vector<std::uint32_t>{0U, 1U, 2U, 7777U}));
  BOOST_TEST((referencedRows(store, 2) == std::vector<std::uint32_t>{200U, 8888U}));
  // Nothing changed since: nothing committed again.
  BOOST_REQUIRE(store.recommitChanged(serial, 2U, tracks, changed) == StageStatus::kComplete);
  BOOST_TEST(changed.empty());
  BOOST_TEST(store.references().recordCount() == before + 6U);
}

/// A table another component adds -- here a matching step's hits per track --
/// starts with an empty range for every track committed so far, grows with
/// the tracks committed after, and holds any number of records per track.
/// One table per record type; the store's sources are named, its own first.
BOOST_AUTO_TEST_CASE(AnAddedTableGrowsWithTheCommittedTracks)
{
  Store store;
  Tracks tracks;
  appendTracks(tracks, {3U, 4U});
  BOOST_REQUIRE(store.commit(serial, 2U, tracks) == StageStatus::kComplete);
  auto& matches = store.addTable<ToyMatch::RecordT>();
  BOOST_TEST(matches.trackCount() == 2U);
  BOOST_TEST(matches.range(1).count == 0U);
  BOOST_CHECK_THROW(store.addTable<ToyMatch::RecordT>(), std::logic_error);
  BOOST_TEST(store.findTable<ToyMatch::RecordT>() == &matches);
  appendTracks(tracks, {5U});
  BOOST_REQUIRE(store.commit(serial, 2U, tracks) == StageStatus::kComplete);
  BOOST_TEST(matches.trackCount() == 3U);
  // Tracks 2 and 0 matched, with three and one hits.
  const std::vector<std::uint32_t> matched{2U, 0U};
  const std::vector<std::uint32_t> hits{3U, 1U};
  const auto status = matches.append(
    serial, 2U, 2U, [&](std::uint32_t i) { return matched[i]; }, [&](std::uint32_t i) { return hits[i]; },
    [&](std::uint32_t i, const auto& out, std::uint32_t at) {
      for (std::uint32_t k = 0; k < hits[i]; ++k) {
        out.hit[at + k] = (1000U * matched[i]) + k;
        out.chi2[at + k] = static_cast<float>(k);
      }
    });
  BOOST_REQUIRE(status == StageStatus::kComplete);
  BOOST_TEST(matches.range(2).count == 3U);
  BOOST_TEST(matches.range(1).count == 0U);
  BOOST_TEST(matches.read(2, 2).hit == 2002U);
  BOOST_TEST(matches.read(0, 0).hit == 0U);
  BOOST_TEST(matches.view().records.chi2[matches.range(2).first + 1U] == 1.F);
  // Sources: the tracker's rows, then those added.
  BOOST_TEST(store.sourceCount() == 1U);
  BOOST_TEST(store.sourceName(kOwnRows) == "rows");
  const auto tof = store.addSource("TOF");
  BOOST_TEST(tof == 1U);
  BOOST_TEST(store.sourceName(tof) == "TOF");
}

/// Tracks committed with overlap partners: each partner beside its row on
/// the same layer, before or after it as it lies along the track, held as a
/// partner; the layer rows still the tracks' own rows (`holds`), the tracker's
/// layout without the partners (`inlineTracks`), `rows` with them.
BOOST_AUTO_TEST_CASE(PartnersGoBesideTheirRows)
{
  Store store;
  Tracks tracks;
  appendTracks(tracks, {2U, 3U});
  BOOST_REQUIRE(store.commit(serial, 2U, tracks) == StageStatus::kComplete);
  appendTracks(tracks, {3U, 2U});
  // Track 2 (rows 200-202): a partner after row 0 and one before row 2;
  // track 3 none.
  const std::vector<std::uint32_t> offsets{0U, 2U, 2U};
  const std::vector<OverlapPartner> partners{{.row = 9000U, .layerRow = 0U, .ahead = 1U}, {.row = 9001U, .layerRow = 2U, .ahead = 0U}};
  BOOST_REQUIRE(store.commit(serial, 2U, tracks, std::span<const std::uint32_t>{offsets}, std::span<const OverlapPartner>{partners}) == StageStatus::kComplete);
  const auto rows = store.rows(2);
  BOOST_TEST((std::vector<std::uint32_t>(rows.begin(), rows.end()) == std::vector<std::uint32_t>{200U, 9000U, 201U, 9001U, 202U}));
  const std::vector<std::uint8_t> roles{1U, 2U, 1U, 2U, 1U};
  for (std::uint32_t k = 0; k < 5U; ++k) {
    BOOST_TEST(store.references().read(2, k).role == static_cast<std::uint8_t>(roles[k] == 1U ? ReferenceRole::kLayerRow : ReferenceRole::kOverlapPartner));
  }
  BOOST_TEST(store.rows(3).size() == 2U);
  BOOST_TEST(store.holds(tracks));
  const auto inlined = inlineTracks(store, std::pmr::get_default_resource());
  BOOST_TEST(inlined.view().count[2] == 3U);
  BOOST_TEST(inlined.view().clusters[2].index[1] == 201U);
  BOOST_TEST(inlined.view().clusters[2].index[3] == ClusterList<kClusters>::kNoCluster);
}

/// Committed tracks given their overlap partners after the run
/// (`recommitPartnered`): the listed ones' references are their rows and
/// partners in path order and their fields the tracks' again; the others keep
/// theirs.
BOOST_AUTO_TEST_CASE(PartnersAddedToCommittedTracks)
{
  Store store;
  Tracks tracks;
  appendTracks(tracks, {3U, 2U, 3U});
  BOOST_REQUIRE(store.commit(serial, 2U, tracks) == StageStatus::kComplete);
  tracks.view().chi2[2] = 7.F;
  // Track 2 (rows 200-202): two partners before row 1, one after row 2.
  const std::vector<std::uint32_t> listed{2U};
  const std::vector<std::uint32_t> offsets{0U, 3U};
  const std::vector<OverlapPartner> partners{{.row = 9000U, .layerRow = 1U, .ahead = 0U}, {.row = 9001U, .layerRow = 1U, .ahead = 0U}, {.row = 9002U, .layerRow = 2U, .ahead = 1U}};
  BOOST_REQUIRE(store.recommitPartnered(serial, 2U, tracks, std::span<const std::uint32_t>{listed}, std::span<const std::uint32_t>{offsets}, std::span<const OverlapPartner>{partners}) == StageStatus::kComplete);
  BOOST_TEST((referencedRows(store, 2) == std::vector<std::uint32_t>{200U, 9000U, 9001U, 201U, 202U, 9002U}));
  BOOST_TEST((referencedRows(store, 0) == std::vector<std::uint32_t>{0U, 1U, 2U}));
  BOOST_TEST(store.tracks().view().chi2[2] == 7.F);
  BOOST_TEST(store.holds(tracks));
}

/// A new timeframe starts with no tracks in any table, the tables still
/// there; the timeframe's end gives every table's storage back to the pool.
BOOST_AUTO_TEST_CASE(TheTimeframesStorageGoesBackToThePool)
{
  auto pool = std::make_shared<BoundedMemoryResource>();
  Store store;
  store.setMemoryPool(pool);
  auto& matches = store.addTable<ToyMatch::RecordT>();
  Tracks tracks;
  appendTracks(tracks, {3U, 4U, 5U});
  BOOST_REQUIRE(store.commit(serial, 2U, tracks) == StageStatus::kComplete);
  BOOST_TEST(pool->getUsedMemory() > 0U);
  store.clear();
  BOOST_TEST(store.size() == 0U);
  BOOST_TEST(store.references().recordCount() == 0U);
  BOOST_TEST(matches.trackCount() == 0U);
  BOOST_TEST(store.findTable<ToyMatch::RecordT>() == &matches);
  BOOST_REQUIRE(store.commit(serial, 2U, tracks) == StageStatus::kComplete);
  BOOST_TEST(store.holds(tracks));
  store.clearStorage();
  BOOST_TEST(pool->getUsedMemory() == 0U);
}

} // namespace o2::gtrk
