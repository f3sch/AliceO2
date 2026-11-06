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

/// \file testAliceLayouts.cxx
/// \brief The ITS and MFT as placed layouts: their tables, the extents the
///        tracking measures from the placements, and the tracker instantiated
///        for both (the envelopes as run-time data).
///
/// With `GTRK_ALICE_GEOMETRY` naming an o2sim geometry file, the placements
/// are also taken from O2's GeometryTGeo and held to its L2G matrices, and the
/// envelopes measured from the pp sample's ideal geometry are held to the
/// table that was compiled in until 2026-10-04; without it those cases pass
/// without checking.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK ALICE layouts
#include <boost/test/unit_test.hpp>

#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Base/LayerEnvelopes.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Tracking/CpuBackend.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/MeasurementGridBuilder.h"
#include "GenericTracking/Tracking/Tracker.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include "DetectorsBase/GeometryManager.h"
#include "ITSBase/GeometryTGeo.h"
#include "MFTBase/GeometryTGeo.h"
#include "MathUtils/Cartesian.h"
#include "MathUtils/Utils.h"

#include "GridTestSupport.h"
#include "UniformMaterial.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace o2::gtrk
{

namespace
{

using detectors::ITS;
using detectors::MFT;

static_assert(Detector<ITS>::nBarrelLayers == 7 && Detector<ITS>::nDiscLayers == 0);
static_assert(Detector<MFT>::nBarrelLayers == 0 && Detector<MFT>::nDiscLayers == 10);
static_assert(ITS::nSensors == 24120 && MFT::nSensors == 936);

/// Every link joins two layers of its detector, and none joins a layer to itself.
template <DetectorLayout Det>
void checkLinks()
{
  for (const auto& link : Det::layerLinks) {
    BOOST_TEST(link.from < Detector<Det>::nLayers);
    BOOST_TEST(link.to < Detector<Det>::nLayers);
    BOOST_TEST(link.from != link.to);
    BOOST_TEST(link.sourceMin < link.sourceMax);
  }
}

/// Where a stand-in layer sits: the ITS's radii and the MFT's disc z, roughly.
template <DetectorLayout Det>
double standInPosition(LayerID layer)
{
  return Detector<Det>::isBarrelLayer(layer) ? 2.3 + (6.2 * static_cast<double>(layer)) : -45.3 - (3.6 * static_cast<double>(layer));
}

/// A stand-in placement: every chip on its layer at the same place, facing
/// along x on a barrel layer and along z on a disc, with the disc chips 8 cm
/// from the beam. Enough to build the tracking geometry; not a detector.
template <DetectorLayout Det>
std::vector<Sensor> standInPlacements(double shiftZ = 0.)
{
  std::vector<Sensor> sensors;
  for (int layer = 0; layer < Detector<Det>::nLayers; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    for (int chip = 0; chip < Detector<Det>::sensorCount(id); ++chip) {
      Sensor sensor{};
      if (Detector<Det>::isBarrelLayer(id)) {
        sensor.rot[3] = sensor.rot[7] = sensor.rot[2] = 1.; // local x, y, z along global y, z, x
        sensor.trans[0] = standInPosition<Det>(id);
      } else {
        sensor.rot[0] = sensor.rot[4] = sensor.rot[8] = 1.; // local axes along the global ones
        sensor.trans[0] = 8.;
        sensor.trans[2] = standInPosition<Det>(id);
      }
      sensor.trans[2] += shiftZ;
      sensor.layer = id;
      sensors.push_back(sensor);
    }
  }
  return sensors;
}

/// The whole tracker, run on an empty timeframe: that it compiles and runs
/// for a barrel-only and a disc-only placed layout is the point.
template <DetectorLayout Det>
void runEmpty()
{
  const auto sensors = standInPlacements<Det>();
  const auto geometry = std::make_unique<const Geometry<Det>>(sensors);
  const auto tracking = std::make_unique<const TrackingGeometry<Det>>(*geometry);
  BOOST_TEST(tracking->surface(0).bounds.halfX == static_cast<Scalar>(detectors::alpide::kSensor.halfX));

  const SearchWindowParams windows{.minPtGeV = 0.3, .bzKiloGauss = 5.};
  const VertexConstraint beamSpot{.z = 0., .sigmaZ = 6., .sigmaTransverse = 0.002};
  const GridGeometry<Det> grid{20'000'000, test::withToyFrames(GridParams{}), test::oneIterationWindows<Det>(windows, beamSpot)};
  std::vector<std::uint32_t> offsets(static_cast<std::size_t>(grid.bucketCount()) + 1U, 0U);
  const MeasurementStore columns;
  const TimeFrame<Det> frame{measurementView(columns), grid, std::span<const std::uint32_t>{offsets}};

  const IterationSettings settings{.windows = windows, .vertexSource = VertexSource::kBeamLine};
  Tracker<Det, UniformMaterial> tracker{settings, beamSpot, scatteringTable<Det>(uniformLayerMaterial<Det>(0.01)),
                                        std::make_shared<BoundedMemoryResource>(), 1U};
  tracker.setTrackModel(*tracking, UniformMaterial{}, Propagator{.bzKiloGauss = 5.F});
  const auto result = tracker.run(frame, tracker.planTimeframe(frame, {}, 1, 0.));
  BOOST_TEST((result.status == StageStatus::kComplete));
  BOOST_REQUIRE(result.tracks != nullptr);
  BOOST_TEST(result.tracks->size() == 0U);
}

/// Our sensor frame against O2's L2G: ALICE's local (x, 0, z) is our
/// (x - kActiveCentreX, z, 0) on the same chip.
template <typename GeometryTGeo>
void checkAgainstL2G(GeometryTGeo& o2Geometry, std::span<const Sensor> sensors)
{
  o2Geometry.fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  double worst = 0.;
  for (std::size_t chip = 0; chip < sensors.size(); chip += 7U) {
    for (const double xa : {-0.6, 0.0581, 0.7}) {
      for (const double za : {-1.4, 0., 1.4}) {
        const auto expected = o2Geometry.getMatrixL2G(static_cast<int>(chip))(o2::math_utils::Point3D<double>{xa, 0., za});
        const double local[3] = {xa - alice::kActiveCentreX, za, 0.};
        double global[3];
        sensors[chip].toGlobal(local, global);
        worst = std::max({worst, std::abs(global[0] - expected.X()), std::abs(global[1] - expected.Y()), std::abs(global[2] - expected.Z())});
      }
    }
    BOOST_TEST(static_cast<int>(sensors[chip].layer) == o2Geometry.getLayer(static_cast<int>(chip)));
  }
  BOOST_TEST(worst < 1e-9);
}

/// The envelope tables compiled into ITS.h and MFT.h until 2026-10-04,
/// measured by o2-gtrk-import-geometry from the pp sample's ideal
/// o2sim_geometry-aligned.root. The measurement at load must give them back
/// to the bit on that geometry, which is what keeps every digest in place.
// Measured values that happen to lie near mathematical constants are not them.
// NOLINTBEGIN(modernize-use-std-numbers)
constexpr std::array<LayerEnvelope, 7> kPastedIts{{
  {.position = 2.3260, .rMin = 2.2331, .rMax = 2.6036, .zMin = -13.5571, .zMax = 13.5571},
  {.position = 3.1354, .rMin = 3.0066, .rMax = 3.4012, .zMin = -13.5571, .zMax = 13.5571},
  {.position = 3.9162, .rMin = 3.7820, .rMax = 4.1614, .zMin = -13.5571, .zMax = 13.5571},
  {.position = 19.5882, .rMin = 19.3518, .rMax = 19.8957, .zMin = -42.2071, .zMax = 42.2071},
  {.position = 24.5272, .rMin = 24.3018, .rMax = 24.8096, .zMin = -42.2071, .zMax = 42.2071},
  {.position = 34.3546, .rMin = 34.1418, .rMax = 34.6084, .zMin = -73.8721, .zMax = 73.8721},
  {.position = 39.3106, .rMin = 39.1018, .rMax = 39.5554, .zMin = -73.8721, .zMax = 73.8721},
}};
constexpr std::array<LayerEnvelope, 10> kPastedMft{{
  {.position = -45.2889, .rMin = 2.3029, .rMax = 12.2259, .zMin = -45.2903, .zMax = -45.2874},
  {.position = -46.7111, .rMin = 2.3029, .rMax = 12.2259, .zMin = -46.7126, .zMax = -46.7097},
  {.position = -48.5889, .rMin = 2.3029, .rMax = 12.2259, .zMin = -48.5903, .zMax = -48.5874},
  {.position = -50.0111, .rMin = 2.3029, .rMax = 12.2259, .zMin = -50.0126, .zMax = -50.0097},
  {.position = -52.3889, .rMin = 2.3029, .rMax = 13.6984, .zMin = -52.3903, .zMax = -52.3874},
  {.position = -53.8111, .rMin = 2.3029, .rMax = 13.6984, .zMin = -53.8126, .zMax = -53.8097},
  {.position = -67.6889, .rMin = 3.4029, .rMax = 16.7089, .zMin = -67.6903, .zMax = -67.6874},
  {.position = -69.1111, .rMin = 3.4029, .rMax = 16.7089, .zMin = -69.1126, .zMax = -69.1097},
  {.position = -76.0889, .rMin = 3.8029, .rMax = 17.3105, .zMin = -76.0903, .zMax = -76.0874},
  {.position = -77.5111, .rMin = 3.8029, .rMax = 17.3105, .zMin = -77.5126, .zMax = -77.5097},
}};
// NOLINTEND(modernize-use-std-numbers)

template <DetectorLayout Det, std::size_t N>
void checkPasted(const std::array<LayerEnvelope, N>& pasted)
{
  static_assert(N == Detector<Det>::nLayers);
  for (std::size_t layer = 0; layer < N; ++layer) {
    const auto& e = PlacedEnvelopes<Det>::all()[layer];
    BOOST_TEST(e.position == pasted[layer].position);
    BOOST_TEST(e.rMin == pasted[layer].rMin);
    BOOST_TEST(e.rMax == pasted[layer].rMax);
    BOOST_TEST(e.zMin == pasted[layer].zMin);
    BOOST_TEST(e.zMax == pasted[layer].zMax);
  }
}

const char* geometryFile() { return std::getenv("GTRK_ALICE_GEOMETRY"); }

} // namespace

// First, before any case constructs a geometry: the extents of a placed layout
// are the loaded geometry's, and there is none yet.
BOOST_AUTO_TEST_CASE(ExtentsNeedAGeometry)
{
  BOOST_CHECK_THROW(static_cast<void>(layerDistance<ITS>(0)), std::logic_error);
  BOOST_CHECK_THROW(static_cast<void>(discZ<MFT>(0)), std::logic_error);
}

BOOST_AUTO_TEST_CASE(LinkTablesJoinLayersOfTheirDetector)
{
  checkLinks<ITS>();
  checkLinks<MFT>();
}

BOOST_AUTO_TEST_CASE(ExtentsAreMeasuredFromThePlacements)
{
  const auto& half = detectors::alpide::kSensor;
  {
    const Geometry<ITS> geometry{standInPlacements<ITS>()};
    for (LayerID layer = 0; layer < Detector<ITS>::nLayers; ++layer) {
      const auto& e = PlacedEnvelopes<ITS>::all()[layer];
      const double radius = standInPosition<ITS>(layer);
      BOOST_TEST(layerDistance<ITS>(layer) == e.position);
      BOOST_TEST(e.position == radius, boost::test_tools::tolerance(1e-9));
      // The chip faces along x: its thickness spans r, its columns z.
      BOOST_TEST(radiusMin<ITS>(layer) == e.rMin);
      BOOST_TEST(e.rMin <= radius - half.halfThickness);
      BOOST_TEST(e.rMin > radius - half.halfThickness - 2. * kEnvelopeResolution);
      BOOST_TEST(radiusMax<ITS>(layer) == e.rMax);
      BOOST_TEST(e.rMax >= std::hypot(radius + half.halfThickness, half.halfX));
      BOOST_TEST(secondMin<ITS>(layer) == e.zMin);
      BOOST_TEST(e.zMin <= -half.halfY);
      BOOST_TEST(secondMax<ITS>(layer) == e.zMax);
      BOOST_TEST(e.zMax >= half.halfY);
      BOOST_TEST(discZ<ITS>(layer) == 0.);
    }
  }
  {
    const Geometry<MFT> geometry{standInPlacements<MFT>()};
    for (LayerID layer = 0; layer < Detector<MFT>::nLayers; ++layer) {
      const auto& e = PlacedEnvelopes<MFT>::all()[layer];
      const double z = standInPosition<MFT>(layer);
      BOOST_TEST(discZ<MFT>(layer) == e.position);
      BOOST_TEST(e.position == z, boost::test_tools::tolerance(1e-9));
      BOOST_TEST(layerDistance<MFT>(layer) == -e.position); // the MFT is on the C side
      BOOST_TEST(secondMin<MFT>(layer) == e.rMin);
      BOOST_TEST(e.rMin <= 8. - half.halfX);
      BOOST_TEST(secondMax<MFT>(layer) == e.rMax);
      BOOST_TEST(e.rMax >= std::hypot(8. + half.halfX, half.halfY));
      BOOST_TEST(zMin<MFT>(layer) == e.zMin);
      BOOST_TEST(e.zMin <= z - half.halfThickness);
      BOOST_TEST(zMax<MFT>(layer) == e.zMax);
      BOOST_TEST(e.zMax >= z + half.halfThickness);
    }
  }
  // A later geometry of the same layout replaces the table: shifted by 1 cm
  // along z, every z extent follows to the micron.
  {
    const Geometry<MFT> shifted{standInPlacements<MFT>(1.)};
    for (LayerID layer = 0; layer < Detector<MFT>::nLayers; ++layer) {
      BOOST_TEST(discZ<MFT>(layer) == standInPosition<MFT>(layer) + 1., boost::test_tools::tolerance(1e-9));
      BOOST_TEST(zMin<MFT>(layer) == PlacedEnvelopes<MFT>::all()[layer].zMin);
    }
  }
}

// The beam line off the axis: the producer subtracts it, so a geometry moved
// with the beam gives every row the same position, radius, azimuth and
// variances about the beam line as the unmoved one about the axis.
BOOST_AUTO_TEST_CASE(GridPointsAreAboutTheBeamLine)
{
  constexpr double kBeamX = -0.057;
  constexpr double kBeamY = -0.0565;
  auto moved = standInPlacements<ITS>();
  for (auto& sensor : moved) {
    sensor.trans[0] += kBeamX;
    sensor.trans[1] += kBeamY;
  }
  const Geometry<ITS> reference{standInPlacements<ITS>()};
  const Geometry<ITS> shifted{moved};
  const auto tolerance = boost::test_tools::tolerance(1e-5F);
  for (const SensorID chip : {SensorID{0}, SensorID{107}, SensorID{24119}}) {
    for (const auto& [loc0, loc1] : {std::pair{0.3, -1.1}, std::pair{-0.6, 0.9}}) {
      const auto a = gridPoint<ITS>(reference, chip, loc0, loc1);
      const auto b = gridPoint<ITS>(shifted, chip, loc0, loc1, kBeamX, kBeamY);
      BOOST_TEST(b.layer == a.layer);
      BOOST_TEST(b.x == a.x, tolerance);
      BOOST_TEST(b.y == a.y, tolerance);
      BOOST_TEST(b.z == a.z, tolerance);
      BOOST_TEST(b.r == a.r, tolerance);
      BOOST_TEST(b.phi == a.phi, tolerance);
      BOOST_TEST(b.second == a.second, tolerance);
      const auto va = pointVariances<ITS>(reference, chip, loc0, a, 1e-6F, 0.F, 4e-6F);
      const auto vb = pointVariances<ITS>(shifted, chip, loc0, b, 1e-6F, 0.F, 4e-6F);
      BOOST_TEST(vb.z == va.z, tolerance);
      BOOST_TEST(vb.r == va.r, tolerance);
      BOOST_TEST(vb.rPhi == va.rPhi, tolerance);
      // And without the subtraction the moved geometry is off by the beam.
      const auto c = gridPoint<ITS>(shifted, chip, loc0, loc1);
      BOOST_TEST(static_cast<double>(c.x - a.x) == kBeamX, boost::test_tools::tolerance(1e-3));
      BOOST_TEST(static_cast<double>(c.y - a.y) == kBeamY, boost::test_tools::tolerance(1e-3));
    }
  }
}

BOOST_AUTO_TEST_CASE(PlacedGeometryRejectsAWrongCount)
{
  auto sensors = standInPlacements<MFT>();
  sensors.pop_back();
  BOOST_CHECK_THROW(Geometry<MFT>{sensors}, std::invalid_argument);
  sensors = standInPlacements<MFT>();
  sensors.back().layer = static_cast<LayerID>(Detector<MFT>::nLayers);
  BOOST_CHECK_THROW(Geometry<MFT>{sensors}, std::invalid_argument);
  // Every chip moved onto layer 0: the count agrees, but a layer is empty.
  sensors = standInPlacements<MFT>();
  for (auto& sensor : sensors) {
    sensor.layer = 0;
  }
  BOOST_CHECK_THROW(Geometry<MFT>{sensors}, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(TrackerRunsOnTheIts) { runEmpty<ITS>(); }

BOOST_AUTO_TEST_CASE(TrackerRunsOnTheMft) { runEmpty<MFT>(); }

BOOST_AUTO_TEST_CASE(PlacementsFollowO2)
{
  if (geometryFile() == nullptr) {
    BOOST_TEST_MESSAGE("GTRK_ALICE_GEOMETRY not set: O2's placements are not checked");
    return;
  }
  o2::base::GeometryManager::loadGeometry(geometryFile());
  const auto its = alice::placeChecked<ITS>();
  const auto mft = alice::placeChecked<MFT>();
  checkAgainstL2G(*o2::its::GeometryTGeo::Instance(), its);
  checkAgainstL2G(*o2::mft::GeometryTGeo::Instance(), mft);

  // The normal convention the frame mapping gives: away from the beam on the
  // ITS, alternating along z on the MFT.
  for (const auto& sensor : its) {
    const double outward = (sensor.rot[2] * sensor.trans[0]) + (sensor.rot[5] * sensor.trans[1]);
    BOOST_TEST_REQUIRE(outward > 0.);
  }
  for (const auto& sensor : mft) {
    BOOST_TEST_REQUIRE(std::abs(sensor.rot[8]) > 0.999999);
  }

  // The envelopes measured at load, against the tables that were compiled in.
  // Only the pp sample's ideal geometry reproduces them; an aligned one does
  // not, and is not meant to.
  if (std::getenv("GTRK_ALICE_GEOMETRY_IDEAL") != nullptr) {
    [[maybe_unused]] const Geometry<ITS> itsGeometry{its};
    [[maybe_unused]] const Geometry<MFT> mftGeometry{mft};
    checkPasted<ITS>(kPastedIts);
    checkPasted<MFT>(kPastedMft);
  } else {
    BOOST_TEST_MESSAGE("GTRK_ALICE_GEOMETRY_IDEAL not set: the measured envelopes are not held to the pasted tables");
  }
}

} // namespace o2::gtrk
