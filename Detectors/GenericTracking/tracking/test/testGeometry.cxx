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

/// \file testGeometry.cxx
/// \brief Exercises the compile-time layout model and the placements it builds.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK Geometry
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Base/LayerEnvelopes.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "ToyBarrel3.h"
#include "ToyBarrel4Disc8.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <set>
#include <span>
#include <stdexcept>

namespace o2::gtrk
{

using detectors::ToyBarrel3;
using detectors::ToyBarrel4Disc8;

// The concept is the contract; if a layout stops satisfying it the failure
// should be here and not in some template instantiation deep in the tracker.
/// Positions recovered by differencing global coordinates carry an absolute
/// error of order (detector size) x (epsilon). In double that is below 1e-13 cm
/// at |z| = 300 cm, so the tolerance below is generous by nine orders of
/// magnitude and exists only to avoid comparing floating point for exact
/// equality.
///
/// It stays absolute rather than relative. The quantities being checked are
/// small differences of large coordinates -- a sub-mm ring stagger read off a
/// 300 cm position -- and a relative tolerance on those measures the wrong
/// thing, which cost two test failures back when this was float.
inline constexpr double kPositionTolerance = 1e-9;

static_assert(DetectorLayout<ToyBarrel3>);
static_assert(DetectorLayout<ToyBarrel4Disc8>);

// The topology each toy exists to represent.
static_assert(Detector<ToyBarrel3>::hasBarrel && !Detector<ToyBarrel3>::hasDiscs);
static_assert(Detector<ToyBarrel4Disc8>::hasBarrel && Detector<ToyBarrel4Disc8>::hasDiscs);
static_assert(Detector<ToyBarrel3>::nLayers == 3);
static_assert(Detector<ToyBarrel4Disc8>::nBarrelLayers == 4);
static_assert(Detector<ToyBarrel4Disc8>::nDiscLayers == 8);
static_assert(Detector<ToyBarrel4Disc8>::isBarrelLayer(0));
static_assert(Detector<ToyBarrel4Disc8>::isBarrelLayer(3));
static_assert(!Detector<ToyBarrel4Disc8>::isBarrelLayer(4));
static_assert(Detector<ToyBarrel4Disc8>::discIndex(4) == 0U);

// Sensor counts are derived, never typed in twice, so this cannot drift.
static_assert(ToyBarrel3::nSensors == 12);

BOOST_AUTO_TEST_CASE(BarrelStaveLengthAccountsForGaps)
{
  constexpr BarrelLayerSpec threeUp{.r = 5., .nStaves = 1, .nModulesPerStave = 3, .phi0 = 0., .tilt = 0., .moduleGap = 0.1, .sensor = ToyBarrel3::sensor};
  // three 2 cm sensors and the two gaps between them
  static_assert(threeUp.nSensors() == 3);
  BOOST_TEST(threeUp.halfLength() == 0.5 * (3. * 2. + 2. * 0.1), boost::test_tools::tolerance(1e-12));
}

/// Chip IDs must be dense, unique, and cover exactly [0, nSensors).
template <DetectorLayout Det>
static void checkIdsAreDense()
{
  const auto geo = toyGeometry<Det>();
  std::set<SensorID> seen;
  for (SensorID id = 0; id < static_cast<SensorID>(Det::nSensors); ++id) {
    seen.insert(id);
    BOOST_TEST(geo.sensor(id).layer < Detector<Det>::nLayers);
  }
  BOOST_TEST(seen.size() == static_cast<std::size_t>(Det::nSensors));
}

BOOST_AUTO_TEST_CASE(ChipIdsAreDenseAndCarryALayer)
{
  checkIdsAreDense<ToyBarrel3>();
}

/// A rotation matrix, not merely an invertible one: orthonormal rows and
/// determinant +1. A left-handed frame would round-trip perfectly and still
/// mirror every cluster position.
template <DetectorLayout Det>
static void checkRotationsAreProper()
{
  const auto geo = toyGeometry<Det>();
  for (SensorID id = 0; id < static_cast<SensorID>(Det::nSensors); ++id) {
    const auto& r = geo.sensor(id).rot;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        const double dot = r[i * 3 + 0] * r[j * 3 + 0] + r[i * 3 + 1] * r[j * 3 + 1] + r[i * 3 + 2] * r[j * 3 + 2];
        BOOST_TEST(dot == (i == j ? 1. : 0.), boost::test_tools::tolerance(1e-12));
      }
    }
    const double det = r[0] * (r[4] * r[8] - r[5] * r[7]) - r[1] * (r[3] * r[8] - r[5] * r[6]) + r[2] * (r[3] * r[7] - r[4] * r[6]);
    BOOST_TEST(det == 1., boost::test_tools::tolerance(1e-12));
  }
}

BOOST_AUTO_TEST_CASE(RotationsAreProperAndRightHanded)
{
  checkRotationsAreProper<ToyBarrel3>();
}

BOOST_AUTO_TEST_CASE(BarrelSensorsSitOnTheirLayerRadius)
{
  const auto geo = toyGeometry<ToyBarrel3>();
  for (SensorID id = 0; id < static_cast<SensorID>(ToyBarrel3::nSensors); ++id) {
    const auto& s = geo.sensor(id);
    BOOST_TEST(std::hypot(s.trans[0], s.trans[1]) == ToyBarrel3::barrelSpecs[s.layer].r, boost::test_tools::tolerance(1e-12));
    BOOST_TEST(s.trans[2] == 0., boost::test_tools::tolerance(1e-12)); // one module per stave
  }
}

BOOST_AUTO_TEST_CASE(BarrelNormalIsRadialAndLongAxisIsBeamwise)
{
  const auto geo = toyGeometry<ToyBarrel3>();
  for (SensorID id = 0; id < static_cast<SensorID>(ToyBarrel3::nSensors); ++id) {
    const auto& s = geo.sensor(id);
    // untilted: the normal (third column) points straight out along the centre
    const double nx = s.rot[2];
    const double ny = s.rot[5];
    BOOST_TEST(nx * s.trans[0] + ny * s.trans[1] == std::hypot(s.trans[0], s.trans[1]), boost::test_tools::tolerance(1e-12));
    // the long axis (second column) is the beam direction
    BOOST_TEST(s.rot[1] == 0., boost::test_tools::tolerance(1e-12));
    BOOST_TEST(s.rot[4] == 0., boost::test_tools::tolerance(1e-12));
    BOOST_TEST(s.rot[7] == 1., boost::test_tools::tolerance(1e-12));
  }
}

BOOST_AUTO_TEST_CASE(SensorsAreOrderedLayerByLayer)
{
  // The ID ordering is a promise the TGeo copy numbers will depend on, so it is
  // pinned here rather than left to whatever the build loops happen to do.
  const auto geo = toyGeometry<ToyBarrel3>();
  LayerID previous = 0;
  for (SensorID id = 0; id < static_cast<SensorID>(ToyBarrel3::nSensors); ++id) {
    BOOST_TEST(geo.sensor(id).layer >= previous);
    previous = geo.sensor(id).layer;
  }
  BOOST_TEST(firstSensorOfLayer<ToyBarrel3>(0) == 0U);
  BOOST_TEST(firstSensorOfLayer<ToyBarrel3>(1) == 4U);
  BOOST_TEST(firstSensorOfLayer<ToyBarrel3>(2) == 8U);
}

// ---------------------------------------------------------------------------
// Barrel and discs together
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(BarrelAndDiscsAreHermeticInAzimuth)
{
  // A tilted stave design covers azimuth with some overlap; so do the
  // trapezoids of a ring, at both radial edges.
  for (const auto& layer : ToyBarrel4Disc8::barrelSpecs) {
    const double covered = static_cast<double>(layer.nStaves) * 2. * layer.sensor.halfX * std::cos(layer.tilt);
    const double circumference = 2. * std::numbers::pi * layer.r;
    const double coverage = covered / circumference;
    BOOST_TEST(coverage > 1.);  // no gaps
    BOOST_TEST(coverage < 1.3); // and not a nonsensical pile-up of silicon
  }
  for (const auto& layer : ToyBarrel4Disc8::discSpecs) {
    for (int i = 0; i < layer.nRings; ++i) {
      const auto& ring = layer.rings[i];
      const double rInner = ring.r - layer.sensor.halfY;
      const double rOuter = ring.r + layer.sensor.halfY;
      const double inner = static_cast<double>(ring.nModules) * 2. * layer.sensor.halfX / (2. * std::numbers::pi * rInner);
      const double outer = static_cast<double>(ring.nModules) * 2. * layer.sensor.halfX / (2. * std::numbers::pi * rOuter);
      BOOST_TEST(inner > 1.);
      BOOST_TEST(outer > 1.);
    }
  }
}

BOOST_AUTO_TEST_CASE(RingModulesAlternateFrontAndBack)
{
  // Azimuthal coverage above 1 means the trapezoids of a ring overlap; if they
  // were also coplanar they would physically intersect, which is what the TGeo
  // overlap check reported (61.7 um per adjacent pair) before the ring gap was
  // taken from the XML. Neighbours must sit a full gap apart in z.
  const auto geo = toyGeometry<ToyBarrel4Disc8>();
  for (LayerID layer = 0; layer < Detector<ToyBarrel4Disc8>::nDiscLayers; ++layer) {
    const auto& disc = ToyBarrel4Disc8::discSpecs[layer];
    SensorID first = firstSensorOfLayer<ToyBarrel4Disc8>(static_cast<LayerID>(Detector<ToyBarrel4Disc8>::nBarrelLayers + layer));
    for (int iRing = 0; iRing < disc.nRings; ++iRing) {
      const auto& ring = disc.rings[iRing];
      // an even module count is what lets the alternation close on itself
      BOOST_TEST(ring.nModules % 2 == 0);
      for (int i = 0; i < ring.nModules; ++i) {
        const double zHere = geo.sensor(first + static_cast<SensorID>(i)).trans[2];
        const double zNext = geo.sensor(first + static_cast<SensorID>((i + 1) % ring.nModules)).trans[2];
        BOOST_TEST(std::abs(std::abs(zHere - zNext) - ring.gap) < kPositionTolerance);
        // and the separation must clear the silicon, or they still touch
        BOOST_TEST(ring.gap > 2. * disc.sensor.halfThickness);
      }
      first += static_cast<SensorID>(ring.nModules);
    }
  }
}

BOOST_AUTO_TEST_CASE(DiscsRunMonotonicallyAlongTheBeam)
{
  // Chip numbering follows the layer order, so the discs must be sorted by z
  // rather than grouped by endcap.
  double previous = -1.e6;
  for (const auto& layer : ToyBarrel4Disc8::discSpecs) {
    BOOST_TEST(layer.z > previous);
    previous = layer.z;
  }
}

BOOST_AUTO_TEST_CASE(BarrelAndDiscGeometryBuildsAndIsSelfConsistent)
{
  checkIdsAreDense<ToyBarrel4Disc8>();
  checkRotationsAreProper<ToyBarrel4Disc8>();
}

BOOST_AUTO_TEST_CASE(SensorsSitWhereTheirLayerSaysTheyDo)
{
  const auto geo = toyGeometry<ToyBarrel4Disc8>();
  for (SensorID id = 0; id < static_cast<SensorID>(ToyBarrel4Disc8::nSensors); ++id) {
    const auto& s = geo.sensor(id);
    if (s.layer < Detector<ToyBarrel4Disc8>::nBarrelLayers) {
      const auto& barrel = ToyBarrel4Disc8::barrelSpecs[s.layer];
      BOOST_TEST(std::hypot(s.trans[0], s.trans[1]) == barrel.r, boost::test_tools::tolerance(1e-12));
      BOOST_TEST(std::abs(s.trans[2]) <= barrel.halfLength());
    } else {
      const auto& disc = ToyBarrel4Disc8::discSpecs[s.layer - Detector<ToyBarrel4Disc8>::nBarrelLayers];
      // The sensor's ring is the one at its radius.
      const double radius = std::hypot(s.trans[0], s.trans[1]);
      const auto rings = std::span{disc.rings}.first(static_cast<std::size_t>(disc.nRings));
      const auto ring = std::ranges::find_if(rings, [&](const RingSpec& candidate) { return std::abs(candidate.r - radius) < 1e-9; });
      BOOST_REQUIRE(ring != rings.end());
      BOOST_TEST(std::abs(std::abs(s.trans[2] - (disc.z + ring->zOffset)) - 0.5 * ring->gap) < kPositionTolerance);
    }
  }
}

namespace
{

/// A curved sensor like the ALICE 3 VD's: a 71 degree arc of a 0.5 cm cylinder,
/// 20 um thick and 50 cm long.
constexpr SensorSpec
  kCurved{.shape = ShapeType::kCylinder, .halfX = 0.3112, .halfY = 25., .halfThickness = 0.001, .pitchX = 0.001, .pitchY = 0.001, .radius = 0.5};

/// A sensor rotated by \a angle about the beam axis and moved by \a t.
Sensor rotatedAboutBeam(double angle, std::array<double, 3> t)
{
  Sensor sensor{};
  const double c = std::cos(angle);
  const double s = std::sin(angle);
  const double rot[9] = {c, -s, 0., s, c, 0., 0., 0., 1.};
  std::ranges::copy(rot, std::begin(sensor.rot));
  std::ranges::copy(t, std::begin(sensor.trans));
  return sensor;
}

} // namespace

BOOST_AUTO_TEST_CASE(CurvedSensorUnrollsOntoItsCylinder)
{
  // The middle of the arc is on the local +x axis, the unrolled z outwards.
  double frame[3];
  kCurved.unrolledToFrame(0., 3., 0.001, frame);
  BOOST_TEST(std::abs(frame[0] - 0.501) < kPositionTolerance);
  BOOST_TEST(std::abs(frame[1]) < kPositionTolerance);
  BOOST_TEST(std::abs(frame[2] - 3.) < kPositionTolerance);
  // The unrolled x is the arc length on the inner surface, counter-clockwise.
  kCurved.unrolledToFrame(0.5 * std::numbers::pi / 4., 0., 0., frame);
  BOOST_TEST(std::abs(frame[0] - frame[1]) < kPositionTolerance);
  BOOST_TEST(frame[1] > 0.);
  static_assert(exactlyEqual(kCurved.surfaceMin(), 0.) && exactlyEqual(kCurved.surfaceMax(), 2. * kCurved.halfThickness));
}

BOOST_AUTO_TEST_CASE(CurvedEnvelopeIsTheHullOfItsSilicon)
{
  // On the beam axis, and moved off it so far that the nearest and farthest
  // points lie inside the arc rather than at its ends.
  for (const auto& t : {std::array<double, 3>{0., 0., 0.}, std::array<double, 3>{0.03, -0.02, 1.}, std::array<double, 3>{0.2, 0.25, -2.}}) {
    const auto sensor = rotatedAboutBeam(std::numbers::pi / 4., t);
    const std::array<SensorSpec, 1> specs{kCurved};
    const auto envelope = measureEnvelopes(std::span(&sensor, 1), specs, 1).front();
    // The silicon sampled densely lies within the hull and reaches it.
    double rMin = 1e9;
    double rMax = 0.;
    for (int i = 0; i <= 2000; ++i) {
      const double x = -kCurved.halfX + (2. * kCurved.halfX * i / 2000.);
      for (const double z : {kCurved.surfaceMin(), kCurved.surfaceMax()}) {
        double frame[3];
        kCurved.unrolledToFrame(x, 0., z, frame);
        double global[3];
        sensor.toGlobal(frame, global);
        rMin = std::min(rMin, std::hypot(global[0], global[1]));
        rMax = std::max(rMax, std::hypot(global[0], global[1]));
      }
    }
    BOOST_TEST(envelope.rMin <= rMin + kPositionTolerance);
    BOOST_TEST(envelope.rMax >= rMax - kPositionTolerance);
    // 2000 steps over the arc miss an extremum by at most a micron.
    BOOST_TEST(rMin - envelope.rMin < 1e-4);
    BOOST_TEST(envelope.rMax - rMax < 1e-4);
    BOOST_TEST(std::abs(envelope.zMin - (t[2] - 25.)) < kPositionTolerance);
    BOOST_TEST(std::abs(envelope.zMax - (t[2] + 25.)) < kPositionTolerance);
  }
  // On the axis the position is the middle of the silicon.
  const auto onAxis = rotatedAboutBeam(1., {0., 0., 0.});
  const std::array<SensorSpec, 1> specs{kCurved};
  BOOST_TEST(std::abs(measureEnvelopes(std::span(&onAxis, 1), specs, 1).front().position - 0.501) < kPositionTolerance);
}

BOOST_AUTO_TEST_CASE(CurvedSensorMustBeAboutTheBeam)
{
  auto tilted = rotatedAboutBeam(0., {0., 0., 0.});
  const double c = std::cos(0.01);
  const double s = std::sin(0.01);
  const double rot[9] = {1., 0., 0., 0., c, -s, 0., s, c};
  std::ranges::copy(rot, std::begin(tilted.rot));
  const std::array<SensorSpec, 1> specs{kCurved};
  BOOST_CHECK_THROW(static_cast<void>(measureEnvelopes(std::span(&tilted, 1), specs, 1)), std::runtime_error);
}

} // namespace o2::gtrk
