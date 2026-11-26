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

#include "GenericTracking/Tools/Material.h"

#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Base/LayerEnvelopes.h"

#include "DetectorsBase/MatLayerCylSet.h"
#include "IOTOFBase/GeometryTGeo.h"
#include "TRKBase/GeometryTGeo.h"

#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMatrix.h>
#include <TGeoNavigator.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TGeoVoxelFinder.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

namespace o2::gtrk::tools
{

std::int32_t fillThreads(std::int32_t threads) { return std::max(threads < 0 ? o2::base::MatLayerCylSet::getNThreadsFromEnv() : threads, 1); }

void forEachOnNavigators(std::size_t count, std::int32_t threads, const std::function<void(std::size_t, TGeoNavigator*)>& work)
{
  threads = std::max(threads, 1);
  if (threads == 1) {
    TGeoNavigator* navigator = gGeoManager->GetCurrentNavigator();
    for (std::size_t i = 0; i < count; ++i) {
      work(i, navigator);
    }
    return;
  }
  gGeoManager->SetMaxThreads(threads);
  // TGeo voxelizes a volume whose voxels need a rebuild when it first looks
  // into it (TGeoVoxelFinder::GetCheckList); several threads would rebuild the
  // same voxels at once.
  TIter volumes(gGeoManager->GetListOfVolumes());
  while (auto* volume = dynamic_cast<TGeoVolume*>(volumes())) {
    auto* voxels = volume->GetVoxels();
    if (voxels != nullptr && voxels->NeedRebuild()) {
      voxels->Voxelize();
      volume->FindOverlaps();
    }
  }
  std::atomic<std::size_t> next{0};
  std::vector<std::thread> workers;
  workers.reserve(static_cast<std::size_t>(threads));
  for (std::int32_t thread = 0; thread < threads; ++thread) {
    workers.emplace_back([&] {
      TGeoNavigator* navigator = gGeoManager->AddNavigator();
      for (std::size_t i = next++; i < count; i = next++) {
        work(i, navigator);
      }
      gGeoManager->RemoveNavigator(navigator);
    });
  }
  for (auto& worker : workers) {
    worker.join();
  }
}

SideWalls trkSideWalls()
{
  if (gGeoManager == nullptr || gGeoManager->GetTopVolume() == nullptr) {
    throw std::runtime_error("trkSideWalls: no geometry loaded");
  }
  SideWalls walls;
  TGeoIterator next(gGeoManager->GetTopVolume());
  while (const TGeoNode* node = next()) {
    if (std::string_view{node->GetVolume()->GetName()} != "VD_SideWall") {
      continue;
    }
    const auto* box = dynamic_cast<const TGeoBBox*>(node->GetVolume()->GetShape());
    if (box == nullptr) {
      throw std::runtime_error("trkSideWalls: a VD_SideWall is not a box");
    }
    const TGeoMatrix* matrix = next.GetCurrentMatrix();
    const auto radius = [&](double x) {
      const Vec3D local{x, 0., 0.};
      Vec3D global{};
      matrix->LocalToMaster(local.data(), global.data());
      return std::hypot(global[0], global[1]);
    };
    const double inner = std::min(radius(-box->GetDX()), radius(box->GetDX()));
    const double outer = std::max(radius(-box->GetDX()), radius(box->GetDX()));
    const double thickness = 2. * box->GetDY();
    walls.thickness = walls.thickness > 0. ? std::min(walls.thickness, thickness) : thickness;
    walls.rMin = walls.rMax > 0. ? std::min(walls.rMin, inner) : inner;
    walls.rMax = std::max(walls.rMax, outer);
  }
  return walls;
}

namespace
{
/// A gas, below this density, is no material: air is 0.0012 g/cm^3.
constexpr double kGasDensity = 0.01; // g/cm^3
} // namespace

std::vector<double> beamLineCylinderFaces(double rMax, double zMax, double dz, std::int32_t threads)
{
  if (gGeoManager == nullptr || gGeoManager->GetTopVolume() == nullptr) {
    throw std::runtime_error("beamLineCylinderFaces: no geometry loaded");
  }
  constexpr std::int32_t kAzimuths = 256;
  constexpr double kResolution = 1e-4; // cm
  const double golden = kPi<double> * (3. - std::sqrt(5.));
  const auto nz = static_cast<std::int32_t>(std::floor(2. * zMax / dz)) + 1;
  const auto facesAt = [&](double z, TGeoNavigator* navigator) {
    std::map<std::int64_t, std::int32_t> seen; // per radius, the azimuths meeting it
    for (std::int32_t ia = 0; ia < kAzimuths; ++ia) {
      const double phi = ia * golden;
      const Vec3D direction{std::cos(phi), std::sin(phi), 0.};
      Vec3D start{0., 0., z};
      navigator->InitTrack(start.data(), direction.data());
      std::set<std::int64_t> mine;
      double r = 0.;
      bool solid = false;
      while (r < rMax && !navigator->IsOutside()) {
        const bool inSolid = navigator->GetCurrentVolume()->GetMaterial()->GetDensity() >= kGasDensity;
        if (inSolid != solid) {
          mine.insert(std::lround(r / kResolution));
          solid = inSolid;
        }
        navigator->FindNextBoundaryAndStep(rMax - r);
        const double step = navigator->GetStep();
        if (step <= 0.) {
          break;
        }
        r += step;
      }
      for (const std::int64_t face : mine) {
        ++seen[face];
      }
    }
    std::vector<std::int64_t> found;
    for (const auto& [face, count] : seen) {
      if (2 * count > kAzimuths && face > 0) {
        found.push_back(face);
      }
    }
    return found;
  };
  std::vector<std::vector<std::int64_t>> byZ(static_cast<std::size_t>(nz));
  forEachOnNavigators(byZ.size(), threads, [&](std::size_t iz, TGeoNavigator* navigator) { byZ[iz] = facesAt(-zMax + (static_cast<double>(iz) * dz), navigator); });
  std::set<std::int64_t> faces;
  for (const auto& found : byZ) {
    faces.insert(found.begin(), found.end());
  }
  std::vector<double> result;
  result.reserve(faces.size());
  for (const std::int64_t face : faces) {
    result.push_back(static_cast<double>(face) * kResolution);
  }
  return result;
}

namespace
{
/// Along each barrel sensor's normal from the centre of its silicon, both
/// ways: the centre must lie in a volume \a isSensor accepts (else \a who
/// throws), then the length of the volumes beyond the sensor that \a beyond
/// accepts, per barrel layer the thickest; both read gGeoManager's current
/// volume.
template <typename IsSensor, typename Beyond>
std::vector<ChipStack> alongNormals(const char* who, std::span<const Sensor> sensors, std::span<const SensorSpec> specs, std::int32_t nBarrelLayers, const IsSensor& isSensor, const Beyond& beyond)
{
  if (gGeoManager == nullptr) {
    throw std::runtime_error(std::format("{}: no geometry loaded", who));
  }
  std::vector<ChipStack> thickness(static_cast<std::size_t>(nBarrelLayers));
  for (const auto& sensor : sensors) {
    if (std::cmp_greater_equal(sensor.layer, nBarrelLayers)) {
      continue;
    }
    // The centre of the silicon, and the normal away from the beam: a flat
    // sensor's local z, a curved one's radius from its axis.
    const auto& spec = specs[sensor.layer];
    const Vec3D centre = siliconCentreOf(sensor, spec);
    double normal[3];
    if (spec.curved()) {
      const double dx = centre[0] - sensor.trans[0];
      const double dy = centre[1] - sensor.trans[1];
      const double norm = std::hypot(dx, dy);
      normal[0] = dx / norm;
      normal[1] = dy / norm;
      normal[2] = 0.;
    } else {
      normal[0] = sensor.rot[2];
      normal[1] = sensor.rot[5];
      normal[2] = sensor.rot[8];
      if ((normal[0] * centre[0]) + (normal[1] * centre[1]) < 0.) {
        for (double& n : normal) {
          n = -n;
        }
      }
    }
    // Out of the sensor along the normal, then through what lies beyond it.
    const auto lengthAlong = [&](double sign) {
      const double direction[3] = {sign * normal[0], sign * normal[1], sign * normal[2]};
      gGeoManager->InitTrack(centre.data(), direction);
      if (!isSensor()) {
        throw std::runtime_error(std::format("{}: the centre of a sensor of layer {} is in {}, not a sensor", who, sensor.layer, gGeoManager->GetCurrentVolume()->GetName()));
      }
      gGeoManager->FindNextBoundaryAndStep();
      double length = 0.;
      while (!gGeoManager->IsOutside() && beyond()) {
        gGeoManager->FindNextBoundaryAndStep();
        length += gGeoManager->GetStep();
      }
      return length;
    };
    auto& layer = thickness[sensor.layer];
    layer.inner = std::max(layer.inner, lengthAlong(-1.));
    layer.outer = std::max(layer.outer, lengthAlong(1.));
  }
  return thickness;
}

/// Whether gGeoManager's current volume's name holds \a pattern.
bool currentNamed(std::string_view pattern)
{
  return std::string_view{gGeoManager->GetCurrentVolume()->GetName()}.find(pattern) != std::string_view::npos;
}
} // namespace

std::vector<ChipStack> trkStackThickness(std::span<const Sensor> sensors, std::span<const SensorSpec> specs, std::int32_t nBarrelLayers)
{
  const std::string sensorPattern = o2::trk::GeometryTGeo::getTRKSensorPattern();
  const std::string stackPattern = o2::trk::GeometryTGeo::getTRKMetalStackPattern();
  return alongNormals("trkStackThickness", sensors, specs, nBarrelLayers, [&] { return currentNamed(sensorPattern); }, [&] { return currentNamed(stackPattern); });
}

std::vector<ChipStack> tofSolidThickness(std::span<const Sensor> sensors, std::span<const SensorSpec> specs, std::int32_t nBarrelLayers)
{
  const std::string inner = o2::iotof::GeometryTGeo::getITOFSensorPattern();
  const std::string outer = o2::iotof::GeometryTGeo::getOTOFSensorPattern();
  return alongNormals(
    "tofSolidThickness", sensors, specs, nBarrelLayers, [&] { return currentNamed(inner) || currentNamed(outer); },
    [] { return gGeoManager->GetCurrentVolume()->GetMaterial()->GetDensity() >= kGasDensity; });
}

std::vector<LayerEnvelope> withStacks(std::span<const LayerEnvelope> envelopes, std::span<const Sensor> sensors, std::span<const SensorSpec> specs, std::span<const ChipStack> stack)
{
  // A flat sensor's silicon is centred on its frame: a face moves by its
  // stack with the half-thickness grown by all of it, and each face is
  // measured on its own. A curved one's starts at its radius, which the inner
  // stack moves in.
  const auto thicker = [&](bool outer) {
    std::vector<SensorSpec> result(specs.begin(), specs.end());
    for (std::uint32_t layer = 0; layer < stack.size(); ++layer) {
      auto& spec = result[layer];
      if (spec.curved()) {
        spec.radius -= stack[layer].inner;
        spec.halfThickness += 0.5 * (stack[layer].inner + stack[layer].outer);
      } else {
        spec.halfThickness += outer ? stack[layer].outer : stack[layer].inner;
      }
    }
    return measureEnvelopes(sensors, result, static_cast<std::int32_t>(stack.size()), kEnvelopeResolution);
  };
  const auto inward = thicker(false);
  const auto outward = thicker(true);
  std::vector<LayerEnvelope> result(envelopes.begin(), envelopes.end());
  for (std::uint32_t layer = 0; layer < stack.size(); ++layer) {
    result[layer].rMin = inward[layer].rMin;
    result[layer].rMax = outward[layer].rMax;
  }
  return result;
}

} // namespace o2::gtrk::tools
