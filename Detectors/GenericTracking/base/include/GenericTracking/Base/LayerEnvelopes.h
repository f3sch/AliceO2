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

/// \file LayerEnvelopes.h
/// \brief Where a layout's layers are: measured from the sensors it was
///        loaded with, and held for the process.
///
/// A layout (the ITS, the MFT, the TRK) has no compile-time positions: where its
/// chips sit is the loaded geometry's business, aligned as that geometry is.
/// Constructing its Geometry measures each layer's envelope from the sensors
/// and publishes the table here; the extent functions of
/// Tracking/SearchWindow.h, and with them the grid, the search windows, the
/// scattering table and the material sampling, read it. Nothing on the device
/// does. Each layer's sensor shape is published beside it: a layout's own
/// where O2 states it at compile time, O2's runtime parameters where O2 holds
/// it there (Alice/Placement.h, `placeSpecs`). The compiled layout keeps
/// topology only: the layer count, each layer's kind and sensor count, the
/// links.

#ifndef GTRK_BASE_LAYERENVELOPES_H_
#define GTRK_BASE_LAYERENVELOPES_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Base/Sensor.h"

#include <array>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace o2::gtrk
{

/// Where a sensor's silicon is centred, in the global frame: a flat sensor's
/// origin, a curved sensor's middle of the silicon at the middle of its arc.
[[nodiscard]] inline Vec3D siliconCentreOf(const Sensor& sensor, const SensorSpec& spec)
{
  double local[3];
  spec.unrolledToFrame(0., 0., 0.5 * (spec.surfaceMin() + spec.surfaceMax()), local);
  Vec3D global{};
  sensor.toGlobal(local, global.data());
  return global;
}

/// The rounding a geometry measures its envelopes with [cm]: 1 um, the
/// minima down, the maxima up, the position to the nearest.
inline constexpr double kEnvelopeResolution = 1e-4;

/// The hull in r and z of one sensor's active silicon, with its effective
/// thickness. A curved sensor's axis must be the beam's (std::runtime_error).
struct SensorExtent {
  double rMin = std::numeric_limits<double>::infinity();
  double rMax = 0.;
  double zMin = std::numeric_limits<double>::infinity();
  double zMax = -std::numeric_limits<double>::infinity();
};

[[nodiscard]] SensorExtent extentOf(const Sensor& sensor, const SensorSpec& spec);

/// The envelope of each layer's active silicon: the hull in r and z of every
/// sensor's active matrix with its effective thickness, and as position the
/// mean radius of its sensors' centres on a barrel or their mean z on a disc.
/// \a specs holds each layer's sensor shape, layers numbered barrel first:
/// the first \a nBarrelLayers are barrel layers. A curved sensor's centre is
/// the middle of its silicon at the middle of its arc; its axis must be the
/// beam's. With a \a resolution, the minima are rounded down to it, the maxima
/// up and the position to the nearest; 0 leaves the values exact. Throws
/// std::runtime_error for a sensor beyond the layers, a layer without one or a
/// curved sensor whose axis is not the beam's.
[[nodiscard]] std::vector<LayerEnvelope> measureEnvelopes(std::span<const Sensor> sensors, std::span<const SensorSpec> specs, int nBarrelLayers, double resolution = 0.);

/// One line per layer, for a log; layers numbered barrel first.
[[nodiscard]] std::string describeEnvelopes(std::span<const LayerEnvelope> envelopes, int nBarrelLayers);

/// One entry per layer of the layout \a Det loaded in this process: its
/// envelope (`PlacedEnvelopes`) or its sensor shape (`PlacedSensorSpecs`).
///
/// One table per layout per process, as O2's `GeometryTGeo::Instance()` holds
/// the one loaded detector for everything that asks: `Geometry<Det>` publishes
/// at construction, a later geometry of the same layout replaces the table,
/// and reading before any exists throws std::logic_error. Published once at
/// set-up, before any tracking thread runs, and read-only after.
template <DetectorLayout Det, typename Entry>
class PlacedLayers
{
 public:
  static void publish(std::span<const Entry> entries)
  {
    if (entries.size() != sTable.size()) {
      throw std::invalid_argument(std::string{Det::name} + ": " + std::to_string(entries.size()) + " placed layers for " + std::to_string(sTable.size()));
    }
    for (std::size_t layer = 0; layer < sTable.size(); ++layer) {
      sTable[layer] = entries[layer];
    }
    sPublished = true;
  }

  [[nodiscard]] static const Entry& get(LayerID layer)
  {
    if (!sPublished) {
      throw std::logic_error(std::string{Det::name} + ": no geometry loaded, so no placed layers; construct Geometry<" + Det::name + "> first");
    }
    if (layer >= sTable.size()) {
      throw std::out_of_range(std::string{Det::name} + ": layer " + std::to_string(layer) + " of " + std::to_string(sTable.size()));
    }
    return sTable[layer];
  }

  [[nodiscard]] static std::span<const Entry> all()
  {
    if (!sPublished) {
      throw std::logic_error(std::string{Det::name} + ": no geometry loaded, so no placed layers; construct Geometry<" + Det::name + "> first");
    }
    return sTable;
  }

 private:
  inline static std::array<Entry, Detector<Det>::nLayers> sTable{};
  inline static bool sPublished = false;
};

template <DetectorLayout Det>
using PlacedEnvelopes = PlacedLayers<Det, LayerEnvelope>;

template <DetectorLayout Det>
using PlacedSensorSpecs = PlacedLayers<Det, SensorSpec>;

/// A layer's material band along its depth -- the radius on a barrel layer,
/// the signed z on a disc [cm]: the extent of the volumes that hold the
/// layer's sensors and none of another layer's (its staves with their
/// supports, cooling and cables), the sensors' envelope within it. Measured
/// from the geometry by the material LUT's builder (alice::materialBands),
/// which makes each barrel layer's band a run of the LUT's shells of its own;
/// the link material reads it back from the LUT (LinkMaterialFill.h).
struct LayerBand {
  double low = 0.;
  double high = 0.;
};

/// One line per layer, for a log; layers numbered barrel first.
[[nodiscard]] std::string describeBands(std::span<const LayerBand> bands, int nBarrelLayers);

} // namespace o2::gtrk

#endif // GTRK_BASE_LAYERENVELOPES_H_
