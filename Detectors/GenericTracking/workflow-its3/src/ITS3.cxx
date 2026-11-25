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

/// \file ITS3.cxx
/// \brief ITS3's own: its chips placed, its cluster converter and its
///        measurement producer.

#include "GenericTracking/Detectors/ITS3.h"

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Workflow/ClusterConverterSpec.h"
#include "GenericTracking/Workflow/MeasurementProducer.h"
#include "GenericTracking/Workflow/MeasurementProducerSpec.h"
#include "GenericTracking/WorkflowITSMFT/ConditionsITSMFT.h"
#include "GenericTracking/WorkflowITSMFT/ITSMFTOptions.h"
#include "GenericTracking/WorkflowITSMFT/Placement.h"

#include "ClusterConverterITS3.h"

#include "ITS3Base/SegmentationMosaix.h"
#include "ITS3Base/SpecsV2.h"
#include "ITSBase/GeometryTGeo.h"
#include "MathUtils/Cartesian.h"
#include "MathUtils/Utils.h"

#include <cmath>
#include <cstddef>
#include <format>
#include <stdexcept>
#include <utility>
#include <vector>

namespace o2::gtrk::alice
{

namespace
{

// ITS3.h restates O2's tile; here it is held to O2's, to the bit.
namespace specs = o2::its3::constants;
using o2::its3::SegmentationMosaix;
namespace mosaix = detectors::mosaix;
static_assert(mosaix::kRows == SegmentationMosaix::NRows && mosaix::kColumns == SegmentationMosaix::NCols);
static_assert(exactlyEqual(mosaix::kWidth, specs::pixelarray::width) && exactlyEqual(mosaix::kLength, specs::pixelarray::length));
static_assert(exactlyEqual(mosaix::kPitchRow, specs::pixelarray::pixels::mosaix::pitchX) && exactlyEqual(mosaix::kPitchColumn, specs::pixelarray::pixels::mosaix::pitchZ));
static_assert(exactlyEqual(mosaix::kSiliconThickness, specs::silicon::thickness) && exactlyEqual(mosaix::kMetalStackThickness, specs::metalstack::thickness));
static_assert(mosaix::kNominalRadius == specs::radii);
static_assert(exactlyEqual(mosaix::kResponseShift, specs::pixelarray::pixels::apts::responseYShift));
static_assert([] {
  for (int layer = 0; std::cmp_less(layer, specs::nLayers); ++layer) {
    const auto i = static_cast<std::size_t>(layer);
    if (!exactlyEqual(mosaix::innerRadius(layer), specs::radiiInner[i]) || !exactlyEqual(mosaix::outerRadius(layer), specs::radiiOuter[i]) || !exactlyEqual(mosaix::middleRadius(layer), specs::radiiMiddle[i])) {
      return false;
    }
  }
  return true;
}());
// The tiles are O2's chips 0 to nChips - 1, layer by layer, the outer barrel's after them.
static_assert(detectors::ITS3::barrelLayers[0].sensorCount == static_cast<int>(specs::detID::l0IDTot));
static_assert(detectors::ITS3::barrelLayers[1].sensorCount == static_cast<int>(specs::detID::l1IDTot));
static_assert(detectors::ITS3::barrelLayers[2].sensorCount == static_cast<int>(specs::detID::l2IDTot));
static_assert(detectors::ITS3::kFirstOuterSensor == static_cast<int>(specs::detID::nChips));

/// Our curved sensor on tile layer \a layer from the tile's L2G matrix.
///
/// O2's tile frame is its cylinder's: the origin on the axis, z along it, the
/// pixel array an arc from angle 0 at the local +x axis, its rows counted
/// from there; every other tile's frame is mirrored (a left-handed L2G). Ours
/// is right-handed whatever O2's: the origin and z O2's, x through the arc's
/// middle (SensorSpec::unrolledToFrame), y = z x x. On a mirrored tile our
/// angle runs against O2's.
Sensor tileSensor(const o2::math_utils::Transform3D& l2g, int layer)
{
  double m[3][3];
  double t[3];
  l2g.GetComponents(m[0][0], m[0][1], m[0][2], t[0], m[1][0], m[1][1], m[1][2], t[1], m[2][0], m[2][1], m[2][2], t[2]);
  const double angle = mosaix::middleAngle(layer);
  Vec3D x{};
  Vec3D z{};
  for (std::size_t i = 0; i < 3; ++i) {
    x[i] = (std::cos(angle) * m[i][0]) + (std::sin(angle) * m[i][1]);
    z[i] = m[i][2];
  }
  const Vec3D y{(z[1] * x[2]) - (z[2] * x[1]), (z[2] * x[0]) - (z[0] * x[2]), (z[0] * x[1]) - (z[1] * x[0])};
  Sensor sensor{};
  for (std::size_t i = 0; i < 3; ++i) {
    sensor.rot[(3 * i) + 0] = x[i];
    sensor.rot[(3 * i) + 1] = y[i];
    sensor.rot[(3 * i) + 2] = z[i];
    sensor.trans[i] = t[i];
  }
  sensor.layer = static_cast<LayerID>(layer);
  return sensor;
}

} // namespace

template <>
std::vector<Sensor> placeChips<detectors::ITS3>()
{
  auto& geometry = *o2::its::GeometryTGeo::Instance();
  geometry.fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  std::vector<Sensor> sensors;
  sensors.reserve(static_cast<std::size_t>(geometry.getNumberOfChips()));
  for (int chip = 0; chip < geometry.getNumberOfChips(); ++chip) {
    const int layer = geometry.getLayer(chip);
    const bool tile = specs::detID::isDetITS3(chip);
    if (tile != std::cmp_less(layer, specs::nLayers)) {
      throw std::runtime_error(std::format("ITS3: chip {} on layer {} is {}a tile: the loaded geometry is not ITS3's", chip, layer, tile ? "" : "not "));
    }
    sensors.push_back(tile ? tileSensor(geometry.getMatrixL2G(chip), layer) : flatSensor(geometry.getMatrixL2G(chip), layer, kActiveCentreX));
  }
  return sensors;
}

template <>
std::array<SensorSpec, Detector<detectors::ITS3>::nLayers> placeSpecs<detectors::ITS3>()
{
  return Geometry<detectors::ITS3>::layoutSpecs();
}

} // namespace o2::gtrk::alice

namespace o2::gtrk::workflow
{

template <>
o2::framework::DataProcessorSpec getClusterConverterSpec<detectors::ITS3>(const o2::framework::ConfigContext& context, bool useMC)
{
  return ITS3ClusterConverterTask::spec(useMC, staggeredInput<detectors::ITS3>(context) ? Detector<detectors::ITS3>::nLayers : 1);
}

template o2::framework::DataProcessorSpec getMeasurementProducerSpec<detectors::ITS3>();

} // namespace o2::gtrk::workflow
