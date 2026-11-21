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

#include "GenericTracking/Alice/Placement.h"

#include <TClass.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoNavigator.h>
#include <TGeoNode.h>
#include <TGeoPcon.h>
#include <TGeoTube.h>
#include <TGeoVolume.h>

#include <algorithm>
#include <bitset>
#include <cmath>
#include <limits>
#include <map>
#include <set>

namespace o2::gtrk::alice
{

namespace
{

/// The rotation and translation of one L2G matrix.
struct Matrix {
  double m[3][3];
  double t[3];
};

Matrix componentsOf(const o2::math_utils::Transform3D& l2g)
{
  Matrix matrix{};
  auto& m = matrix.m;
  auto& t = matrix.t;
  l2g.GetComponents(m[0][0], m[0][1], m[0][2], t[0], m[1][0], m[1][1], m[1][2], t[1], m[2][0], m[2][1], m[2][2], t[2]);
  return matrix;
}

} // namespace

Sensor flatSensor(const o2::math_utils::Transform3D& l2g, int layer, double activeCentreX)
{
  const auto [m, t] = componentsOf(l2g);
  Sensor sensor{};
  for (int i = 0; i < 3; ++i) {
    // Row-major local -> global, so column j is the global direction of our
    // local axis j.
    sensor.rot[(3 * i) + 0] = m[i][0];
    sensor.rot[(3 * i) + 1] = m[i][2];
    sensor.rot[(3 * i) + 2] = -m[i][1];
    sensor.trans[i] = t[i] + (activeCentreX * m[i][0]);
  }
  sensor.layer = static_cast<LayerID>(layer);
  return sensor;
}

Sensor curvedSensor(const o2::math_utils::Transform3D& l2g, int layer)
{
  const auto [m, t] = componentsOf(l2g);
  Sensor sensor{};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      sensor.rot[(3 * i) + j] = m[i][j];
    }
    sensor.trans[i] = t[i];
  }
  sensor.layer = static_cast<LayerID>(layer);
  return sensor;
}

namespace
{

/// The extent of a volume's shape placed by \a matrix: its radius (from the
/// beam axis) and z.
struct ShapeExtent {
  double rMin = std::numeric_limits<double>::infinity();
  double rMax = 0.;
  double zMin = std::numeric_limits<double>::infinity();
  double zMax = -std::numeric_limits<double>::infinity();
};

/// The least distance from the origin to the convex hull of \a points in the
/// plane: 0 if the hull holds the origin.
double leastRadius(const std::vector<std::array<double, 2>>& points)
{
  // The hull's edges by angle around its centroid: the corners of a placed box.
  std::array<double, 2> centre{0., 0.};
  for (const auto& p : points) {
    centre[0] += p[0] / static_cast<double>(points.size());
    centre[1] += p[1] / static_cast<double>(points.size());
  }
  auto hull = points;
  std::ranges::sort(hull, [&](const auto& a, const auto& b) { return std::atan2(a[1] - centre[1], a[0] - centre[0]) < std::atan2(b[1] - centre[1], b[0] - centre[0]); });
  bool inside = true;
  double least = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < hull.size(); ++i) {
    const auto& a = hull[i];
    const auto& b = hull[(i + 1) % hull.size()];
    const double ex = b[0] - a[0];
    const double ey = b[1] - a[1];
    // The origin left of every edge (counter-clockwise) is inside.
    if ((ex * (-a[1])) - (ey * (-a[0])) < 0.) {
      inside = false;
    }
    const double length2 = (ex * ex) + (ey * ey);
    const double t = length2 > 0. ? std::clamp(-((a[0] * ex) + (a[1] * ey)) / length2, 0., 1.) : 0.;
    least = std::min(least, std::hypot(a[0] + (t * ex), a[1] + (t * ey)));
  }
  return inside ? 0. : least;
}

ShapeExtent extentOf(const TGeoShape& shape, const TGeoHMatrix& matrix)
{
  ShapeExtent out;
  const double* rotation = matrix.GetRotationMatrix();
  const double* translation = matrix.GetTranslation();
  const bool alongBeam = std::abs(std::abs(rotation[8]) - 1.) < 1e-9 && std::hypot(translation[0], translation[1]) < kEnvelopeResolution;
  const auto& box = static_cast<const TGeoBBox&>(shape);
  const double* origin = box.GetOrigin();
  std::vector<std::array<double, 2>> corners;
  for (const double sx : {-1., 1.}) {
    for (const double sy : {-1., 1.}) {
      for (const double sz : {-1., 1.}) {
        const double local[3] = {origin[0] + (sx * box.GetDX()), origin[1] + (sy * box.GetDY()), origin[2] + (sz * box.GetDZ())};
        double global[3];
        matrix.LocalToMaster(local, global);
        corners.push_back({global[0], global[1]});
        out.rMax = std::max(out.rMax, std::hypot(global[0], global[1]));
        out.zMin = std::min(out.zMin, global[2]);
        out.zMax = std::max(out.zMax, global[2]);
      }
    }
  }
  out.rMin = leastRadius(corners);
  // A shape about the beam axis: its own radii (a box about the axis holds it).
  if (alongBeam && shape.IsA()->InheritsFrom(TGeoTube::Class())) {
    const auto& tube = static_cast<const TGeoTube&>(shape);
    out.rMin = tube.GetRmin();
    out.rMax = tube.GetRmax();
  } else if (alongBeam && shape.IsA()->InheritsFrom(TGeoPcon::Class())) {
    const auto& pcon = static_cast<const TGeoPcon&>(shape);
    double least = std::numeric_limits<double>::infinity();
    for (int i = 0; i < pcon.GetNz(); ++i) {
      least = std::min(least, pcon.GetRmin(i));
    }
    out.rMin = least;
  } else if (out.rMin == 0.) {
    // A shape whose box holds the axis but not a shape about it: the distance
    // from the axis, at the box's ends and middle, to the shape.
    double least = std::numeric_limits<double>::infinity();
    for (const double z : {out.zMin, 0.5 * (out.zMin + out.zMax), out.zMax}) {
      const double global[3] = {0., 0., z};
      double local[3];
      matrix.MasterToLocal(global, local);
      least = std::min(least, shape.Contains(local) ? 0. : shape.Safety(local, false));
    }
    out.rMin = least;
  }
  return out;
}

} // namespace

std::vector<LayerBand> materialBands(std::span<const Sensor> sensors, std::span<const SensorSpec> specs, std::span<const LayerEnvelope> envelopes, int nBarrelLayers)
{
  if (gGeoManager == nullptr) {
    throw std::runtime_error("material bands: no geometry loaded");
  }
  using Layers = std::bitset<64>;
  const auto layerCount = envelopes.size();
  if (layerCount > Layers{}.size()) {
    throw std::runtime_error("material bands: more layers than the band search holds");
  }
  // Each sensor's node path, and for every volume on it the layers it holds.
  auto* navigator = gGeoManager->GetCurrentNavigator();
  std::vector<std::vector<std::string>> paths(sensors.size());
  std::map<std::string, Layers> held;
  for (std::size_t i = 0; i < sensors.size(); ++i) {
    const auto centre = siliconCentreOf(sensors[i], specs[sensors[i].layer]);
    if (navigator->FindNode(centre[0], centre[1], centre[2]) == nullptr) {
      throw std::runtime_error("material bands: sensor " + std::to_string(i) + "'s centre is in no volume");
    }
    const std::string path = navigator->GetPath();
    std::size_t at = 0;
    while ((at = path.find('/', at + 1)) != std::string::npos) {
      paths[i].push_back(path.substr(0, at));
    }
    paths[i].push_back(path);
    for (const auto& prefix : paths[i]) {
      held[prefix].set(sensors[i].layer);
    }
  }
  // Each layer's outermost volumes holding its sensors alone.
  std::vector<std::set<std::string>> roots(layerCount);
  for (std::size_t i = 0; i < sensors.size(); ++i) {
    for (const auto& prefix : paths[i]) {
      if (held[prefix].count() == 1U) {
        roots[sensors[i].layer].insert(prefix);
        break;
      }
    }
  }
  // The density of the medium a volume sits in: its nearest enclosing volume's
  // that is not an assembly.
  const auto mediumAbove = [](const std::string& path) {
    std::string at = path;
    while (true) {
      const auto slash = at.rfind('/');
      if (slash == 0U || slash == std::string::npos) {
        return 0.;
      }
      at = at.substr(0, slash);
      gGeoManager->cd(at.c_str());
      const auto* volume = gGeoManager->GetCurrentVolume();
      if (!volume->IsAssembly() && volume->GetMaterial() != nullptr) {
        return volume->GetMaterial()->GetDensity();
      }
    }
  };
  std::vector<LayerBand> bands(layerCount);
  for (std::size_t layer = 0; layer < layerCount; ++layer) {
    const bool barrel = static_cast<int>(layer) < nBarrelLayers;
    const auto& e = envelopes[layer];
    double low = barrel ? e.rMin : e.zMin;
    double high = barrel ? e.rMax : e.zMax;
    if (!barrel) {
      bands[layer] = {.low = low, .high = high};
      continue;
    }
    for (const auto& root : roots[layer]) {
      const double medium = mediumAbove(root);
      if (!gGeoManager->cd(root.c_str())) {
        throw std::runtime_error("material bands: no volume at " + root);
      }
      const TGeoHMatrix top = *gGeoManager->GetCurrentMatrix();
      // Volumes of the medium the layer sits in (its air) hold nothing.
      const auto take = [&](const TGeoVolume& volume, const TGeoHMatrix& matrix) {
        if (volume.IsAssembly() || volume.GetMaterial() == nullptr || !(volume.GetMaterial()->GetDensity() > medium)) {
          return;
        }
        const auto extent = extentOf(*volume.GetShape(), matrix);
        low = std::min(low, barrel ? extent.rMin : extent.zMin);
        high = std::max(high, barrel ? extent.rMax : extent.zMax);
      };
      auto* volume = gGeoManager->GetCurrentVolume();
      take(*volume, top);
      TGeoIterator next(volume);
      while (const TGeoNode* node = next()) {
        TGeoHMatrix matrix = top;
        matrix.Multiply(next.GetCurrentMatrix());
        take(*node->GetVolume(), matrix);
      }
    }
    bands[layer] = {.low = low, .high = high};
  }
  // Two barrel layers whose bands meet over a common z end half-way between
  // their envelopes.
  const auto barrels = static_cast<std::size_t>(nBarrelLayers);
  for (std::size_t a = 0; a < barrels; ++a) {
    for (std::size_t b = 0; b < barrels; ++b) {
      const auto& ea = envelopes[a];
      const auto& eb = envelopes[b];
      const double aHigh = ea.rMax;
      const double bLow = eb.rMin;
      const bool across = ea.zMin < eb.zMax && eb.zMin < ea.zMax;
      if (a == b || !across || !(aHigh <= bLow) || !(bands[a].high > bands[b].low)) {
        continue;
      }
      const double middle = 0.5 * (aHigh + bLow);
      bands[a].high = std::max(aHigh, std::min(bands[a].high, middle));
      bands[b].low = std::min(bLow, std::max(bands[b].low, middle));
    }
  }
  return bands;
}

} // namespace o2::gtrk::alice
