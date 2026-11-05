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

/// \file testQualifiers.cxx
/// \brief The tracking headers' unforced `GPUhdi()` functions are the listed stage bodies.
///
/// Every helper of the tracking headers is `GPUhdi() GPUAlwaysInline()`, so that
/// a caller growing past GCC's inlining limits cannot push one out of line --
/// it happened twice, each time a call per candidate. The compiler holds a forced function to
/// that; what it cannot see is a new helper declared `GPUhdi()` alone. This test
/// reads the headers and fails on any such declaration not listed below,
/// and on a listed one that is gone.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Qualifiers
#include <boost/test/unit_test.hpp>

#include <cctype>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <utility>

namespace
{

using Declarations = std::map<std::pair<std::string, std::string>, int>;

/// The functions left to the compiler, by header and name, with their count:
/// what a compiler keeps out of line today (GCC 14 -O2, nvcc, hipcc), apart
/// from the small helpers among it.
const Declarations kStageBodies = {
  // a stage's per-element body
  {{"Follow.h", "operator()"}, 9},
  // the follow's per-track body for a group of lanes
  {{"Follow.h", "together"}, 1},
  {{"RoadCount.h", "operator()"}, 1},
  {{"RoadFit.h", "operator()"}, 1},
  {{"Selection.h", "operator()"}, 3},
  {{"TrackRefit.h", "operator()"}, 1},
  {{"TrackletStage.h", "operator()"}, 1},
  {{"TripletStage.h", "operator()"}, 1},
  // the stages' walks and tests apart from their bodies
  {{"TrackletStage.h", "forEachCandidate"}, 1},
  {{"TrackletStage.h", "testCandidate"}, 1},
  {{"TripletStage.h", "prepareFirst"}, 1},
  {{"TripletStage.h", "forEachWindowCandidate"}, 1},
  {{"TripletStage.h", "testCandidate"}, 1},
  // the walks and drivers around the candidates
  {{"Follow.h", "confirm"}, 1},
  {{"Follow.h", "follow"}, 1},
  {{"Follow.h", "inwardStart"}, 1},
  {{"Follow.h", "kinkStep"}, 1},
  {{"Follow.h", "nominal"}, 1},
  {{"Follow.h", "predictedWindow"}, 1},
  {{"Follow.h", "targetsOf"}, 1},
  // the nominal pass as a confirmation sees it
  {{"Follow.h", "witnessOf"}, 1},
  {{"MeasurementGrid.h", "cellSpan"}, 1},
  {{"MeasurementGrid.h", "visitBuckets"}, 2},
  {{"MeasurementGrid.h", "visitCandidates"}, 2},
  {{"MeasurementGrid.h", "visitIfCompatible"}, 1},
  {{"TrackletStage.h", "accepted"}, 1},
  // forced, its three visitor sites take HIP's vertex-list tracklet kernels
  // from 136 to 248 VGPRs, +14 %
  {{"TrackletStage.h", "expensive"}, 1},
  // large computations; the longitudinal pair test is faster as a call
  //
  {{"CovarianceTransport.h", "transportFactors"}, 1},
  {{"DisplacedPair.h", "displacedPair"}, 1},
  {{"Material.h", "energyLoss"}, 1},
  {{"PairCut.h", "addLongitudinalPairTerms"}, 1},
  {{"PairCut.h", "betweenScattering"}, 1},
  {{"PairCut.h", "betweenSource"}, 1},
  {{"TripletFit.h", "fitAdjacentTriplets"}, 1},
  {{"TripletFit.h", "segmentBackward"}, 1},
  {{"TripletFit.h", "segmentForward"}, 1},
};

/// The declared function's name: the identifier before the parameter list,
/// or `operator...` up to it.
std::string declaredName(const std::string& declaration)
{
  const auto op = declaration.find("operator");
  if (op != std::string::npos) {
    return declaration.substr(op, declaration.find('(', op + 9) - op);
  }
  const auto open = declaration.find('(');
  auto end = open;
  while (end > 0 && declaration[end - 1] == ' ') {
    --end;
  }
  auto begin = end;
  while (begin > 0 && (std::isalnum(static_cast<unsigned char>(declaration[begin - 1])) != 0 || declaration[begin - 1] == '_')) {
    --begin;
  }
  return declaration.substr(begin, end - begin);
}

/// The leaves never inlined on the device, `GPUhdnii()`: the
/// propagation of a state to a plane with the material, which the device otherwise copies
/// into every search and window of every follow variant.
const Declarations kDeviceLeaves = {
  // the propagation of a state to a plane with the material
  {{"Follow.h", "transportCandidate"}, 1},
  // the search of one window, kept out of line so that the kernels do not absorb it once the
  // transport is gone (the TRKFT3 unit doubled its compile otherwise)
  {{"Follow.h", "visit"}, 1},
};

/// The device-only leaves never inlined on the device, O2's `GPUdnii()`.
const Declarations kDeviceOnlyLeaves = {
  // a fit's step: the state carried to the next plane and updated
  {{"TrackFitter.h", "filterStep"}, 1},
};

/// The tracking headers' declarations with \a qualifier (a `GPU...()` word), by header and
/// name; with \a skipForced, those followed by `GPUAlwaysInline()` left out.
Declarations declarationsWith(const std::string& qualifier, bool skipForced)
{
  Declarations found;
  for (const auto& entry : std::filesystem::recursive_directory_iterator(GTRK_TRACKING_HEADERS)) {
    if (entry.path().extension() != ".h") {
      continue;
    }
    std::ifstream in(entry.path());
    std::string line;
    while (std::getline(in, line)) {
      const auto at = line.find(qualifier);
      if (at == std::string::npos || line.find("//") < at) {
        continue;
      }
      // The return type may end the line, the name opening the next one.
      std::string declaration = line.substr(at + qualifier.size());
      std::string next;
      if (declaration.find('(') == std::string::npos && std::getline(in, next)) {
        declaration += " " + next;
      }
      // A helper, forced.
      if (skipForced && declaration.starts_with(" GPUAlwaysInline()")) {
        continue;
      }
      ++found[{std::filesystem::relative(entry.path(), GTRK_TRACKING_HEADERS).string(), declaredName(declaration)}];
    }
  }
  return found;
}

/// \a qualifier is declared on the leaves of \a leaves only, as often as listed.
void checkLeaves(const std::string& qualifier, const Declarations& leaves)
{
  const auto found = declarationsWith(qualifier, false);
  for (const auto& [key, count] : found) {
    const auto listed = leaves.find(key);
    BOOST_TEST_INFO(key.first << " " << key.second << ": " << qualifier << " is for a leaf the device must not copy, listed in this test");
    BOOST_TEST(count == (listed == leaves.end() ? 0 : listed->second));
  }
  for (const auto& [key, count] : leaves) {
    BOOST_TEST_INFO(key.first << " " << key.second << ": listed, but not declared " << qualifier << " that often");
    BOOST_TEST(found.contains(key));
  }
}

} // namespace

BOOST_AUTO_TEST_CASE(GPUhdniiIsOnlyOnTheListedLeaves) { checkLeaves("GPUhdnii()", kDeviceLeaves); }

BOOST_AUTO_TEST_CASE(GPUdniiIsOnlyOnTheListedLeaves) { checkLeaves("GPUdnii()", kDeviceOnlyLeaves); }

BOOST_AUTO_TEST_CASE(GPUhdiIsOnlyOnTheListedStageBodies)
{
  const auto found = declarationsWith("GPUhdi()", true);
  for (const auto& [key, count] : found) {
    const auto listed = kStageBodies.find(key);
    BOOST_TEST_INFO(key.first << " " << key.second << ": a helper is GPUhdi() GPUAlwaysInline(); GPUhdi() alone is for a stage body, listed in this test");
    BOOST_TEST(count == (listed == kStageBodies.end() ? 0 : listed->second));
  }
  for (const auto& [key, count] : kStageBodies) {
    BOOST_TEST_INFO(key.first << " " << key.second << ": listed, but not declared GPUhdi() that often");
    BOOST_TEST(found.contains(key));
  }
}
