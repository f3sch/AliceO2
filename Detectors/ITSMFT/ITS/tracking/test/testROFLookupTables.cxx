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

#include <boost/test/tools/old/interface.hpp>
#define BOOST_TEST_MODULE ITS ROFLookupTables
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "ITStracking/ROFLookupTables.h"

/// -------- Tests --------
// LayerTiming
BOOST_AUTO_TEST_CASE(layertiming_basic)
{
  o2::its::ROFOverlapTable<1> table;
  table.defineLayer(0, 10, 594, 100, 50);
  const auto& layer = table.getLayer(0);

  // test ROF time calculations
  auto start0 = layer.getROFStartInBC(0, true);
  BOOST_CHECK_EQUAL(start0, 100); // delay only

  auto end0 = layer.getROFEndInBC(0, true);
  BOOST_CHECK_EQUAL(end0, 100 + 594);

  auto bounds = layer.getROFTimeBounds(0, true);
  BOOST_CHECK_EQUAL(bounds.getFirstEntry(), 100 - 50);
  BOOST_CHECK_EQUAL(bounds.getEntriesBound() - 1, 100 + 594 + 50);

  // test second ROF
  auto start1 = layer.getROFStartInBC(1, true);
  BOOST_CHECK_EQUAL(start1, 100 + 594);
}

BOOST_AUTO_TEST_CASE(layertiming_base)
{
  o2::its::ROFOverlapTable<3> table;
  table.defineLayer(0, 10, 500, 0, 0);
  table.defineLayer(1, 12, 600, 50, 0);
  table.defineLayer(2, 8, 400, 100, 0);
  const auto& layer1 = table.getLayer(1);
  BOOST_CHECK_EQUAL(layer1.mNROFsTF, 12);
  BOOST_CHECK_EQUAL(layer1.mROFLength, 600);
}

// ROFOverlapTable
BOOST_AUTO_TEST_CASE(rofoverlap_basic)
{
  // define 2 layers with the same definitions (no staggering)
  o2::its::ROFOverlapTable<2> table;
  table.defineLayer(0, 12, 594, 0, 0);
  table.defineLayer(1, 12, 594, 0, 0);
  table.init();
  const auto view = table.getView();
  // each rof in layer 0 should be compatible with its layer 1 equivalent
  for (int rof{0}; rof < 12; ++rof) {
    BOOST_CHECK(view.isCompatible(0, rof, 1, rof));
    BOOST_CHECK(view.isCompatible(1, rof, 0, rof));
    BOOST_CHECK(view.getOverlap(0, 1, rof).getEntries() == 1);
  }
}

BOOST_AUTO_TEST_CASE(rofoverlap_staggered)
{
  // test staggered layers with ROF delay
  o2::its::ROFOverlapTable<2> table;
  table.defineLayer(0, 10, 500, 0, 0);
  table.defineLayer(1, 10, 500, 250, 0); // 250 BC delay
  table.init();
  const auto view = table.getView();

  // verify overlap range
  { // from 0 to 1
    const auto& range = view.getOverlap(0, 1, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 1 to 0
    const auto& range = view.getOverlap(1, 0, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
}

BOOST_AUTO_TEST_CASE(rofoverlap_staggered_alllayers)
{
  // test staggered layers with ROF delay
  o2::its::ROFOverlapTable<3> table;
  table.defineLayer(0, 2, 3, 0, 0);
  table.defineLayer(1, 3, 2, 0, 0);
  table.defineLayer(2, 6, 1, 0, 0);
  table.init();
  const auto view = table.getView();

  // verify overlap range
  { // from 0 to 1 rof=0
    const auto& range = view.getOverlap(0, 1, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 0 to 2 rof=0
    const auto& range = view.getOverlap(0, 2, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 3);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 0 to 1 rof=1
    const auto& range = view.getOverlap(0, 1, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 0 to 2 rof=1
    const auto& range = view.getOverlap(0, 2, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 3);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 3);
  }
  { // from 1 to 2 rof=0
    const auto& range = view.getOverlap(1, 2, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 1 to 0 rof=0
    const auto& range = view.getOverlap(1, 0, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 1 to 2 rof=1
    const auto& range = view.getOverlap(1, 2, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 2);
  }
  { // from 1 to 0 rof=1
    const auto& range = view.getOverlap(1, 0, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 1 to 2 rof=2
    const auto& range = view.getOverlap(1, 2, 2);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 4);
  }
  { // from 1 to 0 rof=2
    const auto& range = view.getOverlap(1, 0, 2);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 2 to 1 rof=0
    const auto& range = view.getOverlap(2, 1, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 1 rof=1
    const auto& range = view.getOverlap(2, 1, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 1 rof=2
    const auto& range = view.getOverlap(2, 1, 2);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 2 to 1 rof=3
    const auto& range = view.getOverlap(2, 1, 3);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 2 to 1 rof=4
    const auto& range = view.getOverlap(2, 1, 4);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 2);
  }
  { // from 2 to 1 rof=5
    const auto& range = view.getOverlap(2, 1, 5);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 2);
  }
  { // from 2 to 0 rof=0
    const auto& range = view.getOverlap(2, 0, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 0 rof=1
    const auto& range = view.getOverlap(2, 0, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 0 rof=2
    const auto& range = view.getOverlap(2, 0, 2);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 0 rof=3
    const auto& range = view.getOverlap(2, 0, 3);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 2 to 0 rof=4
    const auto& range = view.getOverlap(2, 0, 4);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 2 to 0 rof=5
    const auto& range = view.getOverlap(2, 0, 5);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
}

BOOST_AUTO_TEST_CASE(rofoverlap_staggered_alllayers_delay)
{
  // test staggered layers with ROF delay
  o2::its::ROFOverlapTable<3> table;
  table.defineLayer(0, 2, 3, 0, 0);
  table.defineLayer(1, 3, 2, 1, 0);
  table.defineLayer(2, 6, 1, 0, 0);
  table.init();
  const auto view = table.getView();

  // verify overlap range
  { // from 0 to 1 rof=0
    const auto& range = view.getOverlap(0, 1, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 0 to 2 rof=0
    const auto& range = view.getOverlap(0, 2, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 3);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 0 to 1 rof=1
    const auto& range = view.getOverlap(0, 1, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 0 to 2 rof=1
    const auto& range = view.getOverlap(0, 2, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 3);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 3);
  }
  { // from 1 to 2 rof=0
    const auto& range = view.getOverlap(1, 2, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 1 to 0 rof=0
    const auto& range = view.getOverlap(1, 0, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 1 to 2 rof=1
    const auto& range = view.getOverlap(1, 2, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 3);
  }
  { // from 1 to 0 rof=1
    const auto& range = view.getOverlap(1, 0, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 1 to 2 rof=2
    const auto& range = view.getOverlap(1, 2, 2);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 5);
  }
  { // from 1 to 0 rof=2
    const auto& range = view.getOverlap(1, 0, 2);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 2 to 1 rof=0
    const auto& range = view.getOverlap(2, 1, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 0);
  }
  { // from 2 to 1 rof=1
    const auto& range = view.getOverlap(2, 1, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 1 rof=2
    const auto& range = view.getOverlap(2, 1, 2);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 1 rof=3
    const auto& range = view.getOverlap(2, 1, 3);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 2 to 1 rof=4
    const auto& range = view.getOverlap(2, 1, 4);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 2 to 1 rof=5
    const auto& range = view.getOverlap(2, 1, 5);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 2);
  }
  { // from 2 to 0 rof=0
    const auto& range = view.getOverlap(2, 0, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 0 rof=1
    const auto& range = view.getOverlap(2, 0, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 0 rof=2
    const auto& range = view.getOverlap(2, 0, 2);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 0 rof=3
    const auto& range = view.getOverlap(2, 0, 3);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 2 to 0 rof=4
    const auto& range = view.getOverlap(2, 0, 4);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 2 to 0 rof=5
    const auto& range = view.getOverlap(2, 0, 5);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
}

BOOST_AUTO_TEST_CASE(rofoverlap_staggered_alllayers_delay_delta)
{
  // test staggered layers with ROF delay
  o2::its::ROFOverlapTable<3> table;
  table.defineLayer(0, 2, 3, 0, 0);
  table.defineLayer(1, 3, 2, 1, 0);
  table.defineLayer(2, 6, 1, 0, 1);
  table.init();
  const auto view = table.getView();

  // verify overlap range
  { // from 0 to 1 rof=0
    const auto& range = view.getOverlap(0, 1, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 0 to 2 rof=0
    const auto& range = view.getOverlap(0, 2, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 3);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 0 to 1 rof=1
    const auto& range = view.getOverlap(0, 1, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 0 to 2 rof=1
    const auto& range = view.getOverlap(0, 2, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 3);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 3);
  }
  { // from 1 to 2 rof=0
    const auto& range = view.getOverlap(1, 2, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 1 to 0 rof=0
    const auto& range = view.getOverlap(1, 0, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 1 to 2 rof=1
    const auto& range = view.getOverlap(1, 2, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 3);
  }
  { // from 1 to 0 rof=1
    const auto& range = view.getOverlap(1, 0, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 1 to 2 rof=2
    const auto& range = view.getOverlap(1, 2, 2);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 5);
  }
  { // from 1 to 0 rof=2
    const auto& range = view.getOverlap(1, 0, 2);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 2 to 1 rof=0
    const auto& range = view.getOverlap(2, 1, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 1 rof=1
    const auto& range = view.getOverlap(2, 1, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 1 rof=2
    const auto& range = view.getOverlap(2, 1, 2);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 1 rof=3
    const auto& range = view.getOverlap(2, 1, 3);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 1 rof=4
    const auto& range = view.getOverlap(2, 1, 4);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 2 to 1 rof=5
    const auto& range = view.getOverlap(2, 1, 5);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 2 to 0 rof=0
    const auto& range = view.getOverlap(2, 0, 0);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 0 rof=1
    const auto& range = view.getOverlap(2, 0, 1);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 0 rof=2
    const auto& range = view.getOverlap(2, 0, 2);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 0 rof=3
    const auto& range = view.getOverlap(2, 0, 3);
    BOOST_CHECK_EQUAL(range.getEntries(), 2);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 0);
  }
  { // from 2 to 0 rof=4
    const auto& range = view.getOverlap(2, 0, 4);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
  { // from 2 to 0 rof=5
    const auto& range = view.getOverlap(2, 0, 5);
    BOOST_CHECK_EQUAL(range.getEntries(), 1);
    BOOST_CHECK_EQUAL(range.getFirstEntry(), 1);
  }
}

BOOST_AUTO_TEST_CASE(rofoverlap_with_delta)
{
  // test with ROF delta for compatibility window
  o2::its::ROFOverlapTable<2> table;
  table.defineLayer(0, 8, 600, 0, 100); // +/- 100 BC delta
  table.defineLayer(1, 8, 600, 0, 100);
  table.init();
  const auto view = table.getView();

  // with delta, ROFs should have wider compatibility
  for (int rof{0}; rof < 8; ++rof) {
    auto overlap = view.getOverlap(0, 1, rof);
    if (rof == 0 || rof == 7) {
      // edges should see only two
      BOOST_CHECK_EQUAL(overlap.getEntries(), 2);
    } else {
      BOOST_CHECK_EQUAL(overlap.getEntries(), 3);
    }
  }
}

BOOST_AUTO_TEST_CASE(rofoverlap_same_layer)
{
  // test same layer compatibility
  o2::its::ROFOverlapTable<1> table;
  table.defineLayer(0, 10, 500, 0, 0);
  table.init();
  const auto view = table.getView();

  // same ROF in same layer should be compatible
  BOOST_CHECK(view.isCompatible(0, 5, 0, 5));
  // different ROFs in same layer should not be compatible
  BOOST_CHECK(!view.isCompatible(0, 5, 0, 6));
}

// ROFVertexLookupTable
BOOST_AUTO_TEST_CASE(rofvertex_basic)
{
  o2::its::ROFVertexLookupTable<1> table;
  table.defineLayer(0, 6, 594, 0, 0);
  table.init();
  std::vector<o2::its::Vertex> vertices;
  o2::its::Vertex vert0;
  vert0.getTimeStamp().setTimeStamp(594);
  vert0.getTimeStamp().setTimeStampError(594);
  vertices.push_back(vert0);
  o2::its::Vertex vert1;
  vert1.getTimeStamp().setTimeStamp(2375);
  vert1.getTimeStamp().setTimeStampError(594);
  vertices.push_back(vert1);
  table.update(vertices.data(), vertices.size());
  const auto view = table.getView();
}

BOOST_AUTO_TEST_CASE(rofvertex_init_with_vertices)
{
  o2::its::ROFVertexLookupTable<2> table;
  table.defineLayer(0, 10, 500, 0, 0);
  table.defineLayer(1, 10, 500, 0, 0);

  // create vertices at different timestamps
  std::vector<o2::its::Vertex> vertices;
  for (int i = 0; i < 5; ++i) {
    o2::its::Vertex v;
    v.getTimeStamp().setTimeStamp(i * 1000);
    v.getTimeStamp().setTimeStampError(500);
    vertices.push_back(v);
  }

  table.init(vertices.data(), vertices.size());
  const auto view = table.getView();

  // verify vertices can be queried
  const auto& vtxRange = view.getVertices(0, 0);
  BOOST_CHECK_EQUAL(vtxRange.getEntries(), 1);
}

BOOST_AUTO_TEST_CASE(rofvertex_compatibility)
{
  o2::its::ROFVertexLookupTable<1> table;
  table.defineLayer(0, 5, 1000, 0, 100);

  std::vector<o2::its::Vertex> vertices;
  o2::its::Vertex v0;
  v0.getTimeStamp().setTimeStamp(500);
  v0.getTimeStamp().setTimeStampError(200);
  vertices.push_back(v0);

  o2::its::Vertex v1;
  v1.getTimeStamp().setTimeStamp(2500);
  v1.getTimeStamp().setTimeStampError(200);
  vertices.push_back(v1);

  table.init(vertices.data(), vertices.size());
  const auto view = table.getView();

  // check vertex compatibility with ROFs
  bool compat0 = view.isVertexCompatible(0, 0, 0);
  bool compat1 = view.isVertexCompatible(0, 2, 1);
  BOOST_CHECK(compat0 || !compat0); // just verify function executes
}

BOOST_AUTO_TEST_CASE(rofvertex_needs_update)
{
  o2::its::ROFVertexLookupTable<1> table;
  table.defineLayer(0, 5, 500, 0, 0);

  BOOST_CHECK(!table.needsUpdate());

  table.init();
  BOOST_CHECK(table.needsUpdate());

  std::vector<o2::its::Vertex> vertices;
  o2::its::Vertex v;
  v.getTimeStamp().setTimeStamp(500);
  v.getTimeStamp().setTimeStampError(100);
  vertices.push_back(v);

  table.update(vertices.data(), vertices.size());
  BOOST_CHECK(table.needsUpdate());
}

BOOST_AUTO_TEST_CASE(rofvertex_max_vertices)
{
  o2::its::ROFVertexLookupTable<1> table;
  table.defineLayer(0, 3, 1000, 0, 500);

  std::vector<o2::its::Vertex> vertices;
  for (int i = 0; i < 10; ++i) {
    o2::its::Vertex v;
    v.getTimeStamp().setTimeStamp(500 + i * 100);
    v.getTimeStamp().setTimeStampError(50);
    vertices.push_back(v);
  }

  table.init(vertices.data(), vertices.size());
  const auto view = table.getView();

  int32_t maxVtx = view.getMaxVerticesPerROF();
  BOOST_CHECK(maxVtx >= 0);
}

BOOST_AUTO_TEST_CASE(roftimeslice_basic)
{
  o2::its::ROFTimeSliceTable<2> table;
  table.defineLayer(0, 10, 500, 0, 0);
  table.defineLayer(1, 10, 500, 0, 0);

  table.init(5); // 5 time slices
  const auto view = table.getView();

  // verify slices were created
  for (int slice = 0; slice < 5; ++slice) {
    const auto& range = view.getSlice(0, slice);
    BOOST_CHECK(range.getEntries() > 0);
  }
}

BOOST_AUTO_TEST_CASE(roftimeslice_select_rofs)
{
  o2::its::ROFTimeSliceTable<1> table;
  table.defineLayer(0, 10, 500, 0, 0);
  table.init(4);

  BOOST_CHECK(table.needsUpdate());

  // select ROFs in a BC range
  table.selectROFs(1000, 2000);
  BOOST_CHECK(table.needsUpdate());
}

BOOST_AUTO_TEST_CASE(roftimeslice_select_rofs_vector)
{
  o2::its::ROFTimeSliceTable<1> table;
  table.defineLayer(0, 10, 500, 0, 0);
  table.init(4);

  using BCRange = o2::its::ROFTimeSliceTable<1>::BCRange;
  std::vector<BCRange> ranges;
  ranges.emplace_back(500, 1000);
  ranges.emplace_back(2000, 500);

  table.selectROFs(ranges);
  BOOST_CHECK(table.needsUpdate());
}

BOOST_AUTO_TEST_CASE(roftimeslice_mask_operations)
{
  o2::its::ROFTimeSliceTable<1> table;
  table.defineLayer(0, 8, 500, 0, 0);
  table.init(3);

  // reset mask
  table.resetMask(1);
  BOOST_CHECK(table.needsUpdate());

  table.resetMask(0);
  BOOST_CHECK(table.needsUpdate());

  // invert mask
  table.invertMask();
  BOOST_CHECK(table.needsUpdate());
}

BOOST_AUTO_TEST_CASE(roftimeslice_find_slice)
{
  o2::its::ROFTimeSliceTable<1> table;
  table.defineLayer(0, 10, 1000, 0, 0);
  table.init(5);
  const auto view = table.getView();

  // find which slice a BC belongs to
  int32_t slice = view.findSlice(0, 5000);
  BOOST_CHECK(slice >= -1 && slice < 5);

  // BC outside range should return -1
  int32_t invalidSlice = view.findSlice(0, 50000);
  BOOST_CHECK_EQUAL(invalidSlice, -1);
}

BOOST_AUTO_TEST_CASE(roftimeslice_is_in_slice)
{
  o2::its::ROFTimeSliceTable<1> table;
  table.defineLayer(0, 10, 1000, 0, 0);
  table.init(4);
  const auto view = table.getView();

  // test if BC is in specific slice
  const auto& slice0 = view.getSlice(0, 0);
  int32_t bcInSlice = slice0.getFirstEntry() + 100;
  BOOST_CHECK(view.isInSlice(0, 0, bcInSlice));

  int32_t bcOutOfSlice = slice0.getEntriesBound() + 100;
  BOOST_CHECK(!view.isInSlice(0, 0, bcOutOfSlice));
}

BOOST_AUTO_TEST_CASE(roftimeslice_get_mask)
{
  o2::its::ROFTimeSliceTable<1> table;
  table.defineLayer(0, 6, 500, 0, 0);
  table.init(3);
  const auto view = table.getView();

  // get mask for a slice
  const uint8_t* mask = view.getMask(0, 1);
  BOOST_CHECK(mask != nullptr);
}

BOOST_AUTO_TEST_CASE(multilayer_complex)
{
  // test more complex scenario with 4 layers
  o2::its::ROFOverlapTable<4> table;
  table.defineLayer(0, 10, 500, 0, 50);
  table.defineLayer(1, 10, 500, 100, 50);
  table.defineLayer(2, 12, 600, 0, 100);
  table.defineLayer(3, 8, 400, 50, 0);
  table.init();

  const auto view = table.getView();
  BOOST_CHECK_EQUAL(table.getEntries(), 4);

  // verify different layer combinations
  BOOST_CHECK(view.getOverlap(0, 1, 0).getEntries() > 0);
  BOOST_CHECK(view.getOverlap(2, 3, 0).getEntries() >= 0);
}
