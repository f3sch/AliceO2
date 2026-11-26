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

/// \file layerLinks.cxx
/// \brief Trace a detector's layer links from its layout.
///
/// **Trace.** A barrel layer is phi-symmetric, the segment `r = R,
/// z in [-L, L]` of its stated envelope in (r, z). A helix is a straight line
/// in (s, z) with `r(s) = 2 R_c sin(s / 2 R_c)`, so its barrel crossings are
/// closed form. Trajectories run over z0 in the given beam-spot
/// range, eta over the given acceptance and curvature from zero to the floor
/// pT's, each followed over its outgoing half-turn only. Consecutive crossings
/// form links; a link's source interval is the hull of the source coordinate
/// over its crossings, padded by the margin. Whether a pair on a link leaves a
/// hole is the pair's own and the iteration's (Tracking/LinkHoles.h).
///
/// The printed table is pasted into the detector's header. Nothing checks the
/// two against each other afterwards; regenerating the table is part of
/// changing a layout.
///
/// **Discs.** Every layout is placed from the loaded geometry (ITS, MFT,
/// TRKFT3), and a disc layer need not be a ring: an MFT face covers only part
/// of its disc, so a track crosses one face of a disc or the other. The discs'
/// crossings are therefore found in three dimensions, over azimuth and both
/// charges, against a raster of the sensors' active silicon placed from the
/// geometry.

#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/ITS3.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Detectors/TRKFT3TOF.h"
#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Tracking/SearchWindow.h"
#include "GenericTracking/Workflow/GeometryProvider.h"
#include "GenericTracking/WorkflowTRK/TRK.h"
#include "GenericTracking/Tools/SensorFrame.h"

#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

#include <boost/program_options.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace bpo = boost::program_options;

namespace o2::gtrk::tools
{
namespace
{

/// The trace's settings, all but the field's Bz the options' (their defaults in `run`).
struct Options {
  std::string detector;
  double beamSpotZ;
  double beamSpotHalfRange;
  double minPtGeV;
  std::string field;       ///< the run's GRPMagField file
  double bzKiloGauss = 0.; ///< its nominal Bz, set once the field is loaded
  double etaMax;
  std::int32_t etaSteps;
  std::int32_t z0Steps;
  std::int32_t curvatureSteps;
  double margin;
  /// Links carrying less than this fraction of their source layer's crossings
  /// are left out of the printed table (the shares are those of the full trace).
  double minShare;
  std::string geometry;
  std::int32_t phiSteps;
};

/// One barrel layer's surface in (r, z): its stated envelope; the discs are
/// traced against their placed silicon (`discCoverage`).
struct Surface {
  LayerID layer = 0;
  double position = 0.; ///< R [cm]
  double low = 0.;      ///< extent along z [cm]
  double high = 0.;
};

template <DetectorLayout Det>
std::vector<Surface> surfaces()
{
  std::vector<Surface> result;
  for (std::int32_t layer = 0; layer < Detector<Det>::nLayers; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    if (Detector<Det>::isBarrelLayer(id)) {
      result.push_back({.layer = id, .position = layerDistance<Det>(id), .low = secondMin<Det>(id), .high = secondMax<Det>(id)});
    }
  }
  return result;
}

/// A crossing of one layer by one trajectory.
struct Crossing {
  double path = 0.; ///< transverse path length from the origin
  LayerID layer = 0;
  double second = 0.; ///< z on a barrel layer, r on a disc [cm]
};

/// What is known about one link from either source.
struct LinkStats {
  std::uint32_t count = 0;
  double low = std::numeric_limits<double>::infinity();
  double high = -std::numeric_limits<double>::infinity();

  void add(double second)
  {
    ++count;
    low = std::min(low, second);
    high = std::max(high, second);
  }
};

using LinkMap = std::map<std::pair<LayerID, LayerID>, LinkStats>;

/// Visit the consecutive transitions of one list of crossings, taking each
/// layer's first crossing only: the link, and its two crossings.
template <typename Visitor>
void forEachTransition(std::vector<Crossing>& crossings, std::vector<bool>& seen, Visitor visit)
{
  std::ranges::sort(crossings, {}, &Crossing::path);
  seen.assign(seen.size(), false);
  const Crossing* previous = nullptr;
  for (const auto& crossing : crossings) {
    if (seen[crossing.layer]) {
      continue; // a second ring of the same disc, or a sensor overlap
    }
    seen[crossing.layer] = true;
    if (previous != nullptr) {
      visit(std::pair{previous->layer, crossing.layer}, *previous, crossing);
    }
    previous = &crossing;
  }
}

/// Which points of a placed disc layer's plane are on active silicon, as a
/// raster of kCell squares centred on the beam line.
struct DiscCoverage {
  static constexpr double kCell = 0.01; // [cm]

  LayerID layer = 0;
  double z = 0.;    ///< [cm]
  double rMin = 0.; ///< [cm]
  double rMax = 0.;
  std::int32_t cells = 0; ///< per side
  std::vector<bool> covered;

  [[nodiscard]] double half() const { return 0.5 * kCell * cells; }

  [[nodiscard]] bool contains(double x, double y) const
  {
    const double ix = std::floor((x + half()) / kCell);
    const double iy = std::floor((y + half()) / kCell);
    if (ix < 0. || iy < 0. || ix >= cells || iy >= cells) {
      return false;
    }
    return covered[(static_cast<std::uint32_t>(iy) * static_cast<std::uint32_t>(cells)) + static_cast<std::uint32_t>(ix)];
  }
};

/// The coverage of every disc layer of a placed geometry. A cell is covered
/// when its centre, on the sensor's plane, is inside the sensor's bounds.
template <DetectorLayout Det>
std::vector<DiscCoverage> discCoverage(const Geometry<Det>& geometry)
{
  std::vector<DiscCoverage> result;
  for (std::int32_t layer = Detector<Det>::nBarrelLayers; layer < Detector<Det>::nLayers; ++layer) {
    const auto id = static_cast<LayerID>(layer);
    DiscCoverage coverage{.layer = id, .z = discZ<Det>(id), .rMin = secondMin<Det>(id), .rMax = secondMax<Det>(id), .cells = 0, .covered = {}};
    coverage.cells = 2 * static_cast<std::int32_t>(std::ceil(coverage.rMax / DiscCoverage::kCell)) + 2;
    coverage.covered.assign(static_cast<std::size_t>(coverage.cells) * static_cast<std::size_t>(coverage.cells), false);
    const double half = coverage.half();
    for (const auto& sensor : geometry.sensors()) {
      if (sensor.layer != id) {
        continue;
      }
      const auto& spec = geometry.layerSpec(id);
      // The sensor's bounding box in the raster, from its corners.
      double low[2] = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
      double high[2] = {-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
      for (const double u : {-spec.halfX, spec.halfX}) {
        for (const double v : {-spec.halfY, spec.halfY}) {
          const double local[3] = {u, v, 0.};
          double global[3];
          sensor.toGlobal(local, global);
          for (std::int32_t i = 0; i < 2; ++i) {
            low[i] = std::min(low[i], global[i]);
            high[i] = std::max(high[i], global[i]);
          }
        }
      }
      const auto cell = [half](double value) { return static_cast<std::int32_t>(std::floor((value + half) / DiscCoverage::kCell)); };
      for (std::int32_t iy = std::max(cell(low[1]), 0); iy <= std::min(cell(high[1]), coverage.cells - 1); ++iy) {
        for (std::int32_t ix = std::max(cell(low[0]), 0); ix <= std::min(cell(high[0]), coverage.cells - 1); ++ix) {
          const double global[3] = {((ix + 0.5) * DiscCoverage::kCell) - half, ((iy + 0.5) * DiscCoverage::kCell) - half, sensor.trans[2]};
          double local[3];
          toLocal(sensor, global, local);
          if (contains(spec, local[0], local[1])) {
            coverage.covered[static_cast<std::uint32_t>((iy * coverage.cells) + ix)] = true;
          }
        }
      }
    }
    result.push_back(std::move(coverage));
  }
  return result;
}

template <DetectorLayout Det>
LinkMap trace(const Options& options, std::span<const DiscCoverage> placedDiscs)
{
  const auto layout = surfaces<Det>();
  const double maxCurvature = 1. / curvatureRadius(options.minPtGeV, options.bzKiloGauss);

  LinkMap links;
  std::vector<Crossing> crossings;
  std::vector<Crossing> all;
  std::vector<std::pair<const DiscCoverage*, double>> reached; // a placed disc and the transverse path to it
  std::vector<bool> seen(static_cast<std::size_t>(Detector<Det>::nLayers));
  const auto step = [](double low, double high, std::int32_t steps, std::int32_t i) { return steps > 1 ? low + ((high - low) * i / (steps - 1)) : 0.5 * (low + high); };
  for (std::int32_t ie = 0; ie < options.etaSteps; ++ie) {
    const double cotTheta = std::sinh(step(-options.etaMax, options.etaMax, options.etaSteps, ie));
    for (std::int32_t iz = 0; iz < options.z0Steps; ++iz) {
      const double z0 = step(options.beamSpotZ - options.beamSpotHalfRange, options.beamSpotZ + options.beamSpotHalfRange, options.z0Steps, iz);
      for (std::int32_t ik = 0; ik < options.curvatureSteps; ++ik) {
        const double k = options.curvatureSteps > 1 ? maxCurvature * ik / (options.curvatureSteps - 1) : 0.;
        const double halfTurn = k > 0. ? kPi<double> / k : std::numeric_limits<double>::infinity();
        crossings.clear();
        for (const auto& surface : layout) {
          if (k * surface.position > 2.) {
            continue; // turns back before it gets there
          }
          const double s = k > 0. ? 2. / k * std::asin(0.5 * k * surface.position) : surface.position;
          const double z = z0 + (s * cotTheta);
          if (z >= surface.low && z <= surface.high) {
            crossings.push_back({.path = s, .layer = surface.layer, .second = z});
          }
        }
        const auto count = [&links](auto link, const Crossing& from, const Crossing& /*to*/) { links[link].add(from.second); };

        // Where a placed disc's plane is reached does not depend on azimuth
        // or charge; whether it is on silicon there does.
        reached.clear();
        if (!exactlyEqual(cotTheta, 0.)) {
          for (const auto& disc : placedDiscs) {
            const double s = (disc.z - z0) / cotTheta;
            const double r = k > 0. ? 2. / k * std::sin(0.5 * k * s) : s;
            if (s > 0. && s <= halfTurn && r >= disc.rMin && r <= disc.rMax) {
              reached.emplace_back(&disc, s);
            }
          }
        }
        if (reached.empty()) {
          forEachTransition(crossings, seen, count);
          continue;
        }
        for (const double charge : {-1., 1.}) {
          for (std::int32_t ip = 0; ip < options.phiSteps; ++ip) {
            const double phi0 = kTwoPi<double> * ip / options.phiSteps;
            all = crossings;
            for (const auto& [disc, s] : reached) {
              const double turn = charge * k;
              const double x = k > 0. ? (std::sin(phi0 + (turn * s)) - std::sin(phi0)) / turn : s * std::cos(phi0);
              const double y = k > 0. ? (std::cos(phi0) - std::cos(phi0 + (turn * s))) / turn : s * std::sin(phi0);
              if (disc->contains(x, y)) {
                all.push_back({.path = s, .layer = disc->layer, .second = std::hypot(x, y)});
              }
            }
            forEachTransition(all, seen, count);
          }
          if (!(k > 0.)) {
            break; // a straight line has no charge
          }
        }
      }
    }
  }
  return links;
}

template <DetectorLayout Det>
std::string layerName(LayerID layer)
{
  std::ostringstream name;
  if (Detector<Det>::isBarrelLayer(layer)) {
    name << layer << " (barrel r " << std::defaultfloat << std::setprecision(4) << layerDistance<Det>(layer) << ")";
  } else {
    name << layer << " (disc z " << std::showpos << std::defaultfloat << std::setprecision(4) << discZ<Det>(layer) << std::noshowpos << ")";
  }
  return name.str();
}

/// The fraction of each link's source layer's transitions that take it.
std::map<std::pair<LayerID, LayerID>, double> weights(const LinkMap& links)
{
  std::map<LayerID, std::uint32_t> leaving;
  for (const auto& [link, stats] : links) {
    leaving[link.first] += stats.count;
  }
  std::map<std::pair<LayerID, LayerID>, double> result;
  for (const auto& [link, stats] : links) {
    result[link] = static_cast<double>(stats.count) / static_cast<double>(leaving[link.first]);
  }
  return result;
}

/// The links at or above \a minShare of their source layer's crossings, the
/// shares those of the full map.
LinkMap pruned(const LinkMap& links, double minShare)
{
  const auto weight = weights(links);
  LinkMap result;
  for (const auto& [link, stats] : links) {
    if (weight.at(link) >= minShare) {
      result.emplace(link, stats);
    }
  }
  return result;
}

template <DetectorLayout Det>
void printTable(const LinkMap& links, const Options& options)
{
  const auto weight = weights(links);
  const auto kept = pruned(links, options.minShare);
  std::cout << std::defaultfloat << "  // Traced by o2-gtrk-layer-links for " << Det::name << ": pT >= " << options.minPtGeV << " GeV at Bz = " << options.bzKiloGauss << " kG,\n"
            << "  // beam spot z in [" << options.beamSpotZ - options.beamSpotHalfRange << ", " << options.beamSpotZ + options.beamSpotHalfRange << "] cm, |eta| <= " << options.etaMax << ", source intervals padded by " << options.margin << " cm.\n"
            << "  // The comment on each link is the fraction of its source layer's crossings that take it.\n";
  if (options.minShare > 0.) {
    std::cout << "  // Links below " << options.minShare << " of their source layer's crossings are left out (" << links.size() - kept.size() << " of " << links.size() << ").\n";
  }
  std::cout << "  static constexpr std::array<LayerLink, " << kept.size() << "> layerLinks{{\n";
  for (const auto& [link, stats] : kept) {
    const double low = std::floor((stats.low - options.margin) * 100.) / 100.;
    const double high = std::ceil((stats.high + options.margin) * 100.) / 100.;
    std::cout << "    {.from = " << link.first << ", .to = " << link.second << std::fixed << std::setprecision(2) << ", .sourceMin = " << low << ", .sourceMax = " << high << "}, // " << std::setprecision(4) << weight.at(link) << '\n';
  }
  std::cout << "  }};\n";
}

template <DetectorLayout Det>
void printSummary(const LinkMap& links)
{
  std::uint32_t twoLinkPaths = 0;
  std::map<LayerID, std::uint32_t> outDegree;
  for (const auto& [link, stats] : links) {
    ++outDegree[link.first];
  }
  std::uint32_t maxOut = 0;
  for (const auto& [link, stats] : links) {
    twoLinkPaths += outDegree.contains(link.second) ? outDegree[link.second] : 0U;
    maxOut = std::max(maxOut, outDegree[link.first]);
  }
  std::cerr << Det::name << ": " << Detector<Det>::nLayers << " layers, " << links.size() << " links, " << twoLinkPaths << " two-link paths, at most " << maxOut << " links leaving a layer\n";
}

template <DetectorLayout Det>
std::int32_t runFor(const Options& options)
{
  // A layout's extents are its loaded geometry's, on the ITS as well as on the MFT.
  o2::base::GeometryManager::loadGeometry(options.geometry);
  const auto geometry = workflow::makeGeometry<Det>();
  std::vector<DiscCoverage> coverage;
  if constexpr (Detector<Det>::hasDiscs) {
    coverage = discCoverage<Det>(*geometry);
  }
  const auto traced = trace<Det>(options, coverage);
  printSummary<Det>(traced);
  printTable<Det>(traced, options);
  return 0;
}

std::int32_t run(std::int32_t argc, char** argv)
{
  Options options;
  bpo::options_description description("Trace a detector's layer links");
  description.add_options()                                                                                      //
    ("help,h", "this message")                                                                                   //
    ("detector,d", bpo::value(&options.detector)->required(), "ITS, ITS3, MFT, TRKFT3 or TRKFT3TOF")             //
    ("beam-spot-z", bpo::value(&options.beamSpotZ)->default_value(0.), "centre of the traced vertex range [cm]") //
    ("beam-spot-range", bpo::value(&options.beamSpotHalfRange)->default_value(18.),
     "half-width of the traced vertex range [cm]; 3 sigma of a 6 cm spot by default") //
    ("min-pt", bpo::value(&options.minPtGeV)->default_value(0.1),
     "floor pT [GeV]: the lowest any iteration will ask for") //
    ("grp-mag-field", bpo::value(&options.field)->default_value(""),
     "the run's GRPMagField, e.g. <sample>/ccdb/GLO/Config/GRPMagField/snapshot.root; the trace takes its nominal Bz. Not for TRKFT3, whose field is "
     "ALICE3_MAGFIELD_MACRO's")                                                       //
    ("eta-max", bpo::value(&options.etaMax)->default_value(6.), "traced |eta| range") //
    ("eta-steps", bpo::value(&options.etaSteps)->default_value(6001), "eta samples")  //
    ("z0-steps", bpo::value(&options.z0Steps)->default_value(13), "vertex samples")   //
    ("curvature-steps", bpo::value(&options.curvatureSteps)->default_value(201),
     "curvature samples, from a straight line to the floor pT") //
    ("min-share", bpo::value(&options.minShare)->default_value(0.),
     "leave links carrying less than this fraction of their source layer's crossings out of the printed table") //
    ("margin", bpo::value(&options.margin)->default_value(0.5),
     "padding on each end of a source interval [cm]") //
    ("geometry,g", bpo::value(&options.geometry)->default_value("o2sim_geometry-aligned.root"),
     "o2sim geometry a placed layout's discs are traced against") //
    ("phi-steps", bpo::value(&options.phiSteps)->default_value(180), "azimuth samples, per charge, where a placed disc is reached");

  bpo::variables_map vm;
  try {
    bpo::store(bpo::parse_command_line(argc, argv, description), vm);
    if (vm.contains("help")) {
      std::cout << description << '\n';
      return 0;
    }
    bpo::notify(vm);
  } catch (const bpo::error& error) {
    std::cerr << error.what() << "\n\n"
              << description << '\n';
    return 1;
  }
  // The field as the tracking takes it: O2's nominal Bz of the run's map;
  // for ALICE 3, whose GRP records none, its field object's Bz at the origin.
  if (options.detector == "TRKFT3" || options.detector == "TRKFT3TOF") {
    options.bzKiloGauss = alice::alice3BzKiloGauss();
  } else {
    if (options.field.empty()) {
      std::cerr << "no field: --grp-mag-field names the run's GRPMagField\n";
      return 1;
    }
    const std::unique_ptr<o2::parameters::GRPMagField> grp{o2::parameters::GRPMagField::loadFrom(options.field)};
    o2::base::Propagator::initFieldFromGRP(grp.get());
    options.bzKiloGauss = static_cast<double>(o2::base::Propagator::Instance()->getNominalBz());
  }

  using namespace detectors;
  if (options.detector == "ITS") {
    return runFor<ITS>(options);
  }
  if (options.detector == "ITS3") {
    return runFor<ITS3>(options);
  }
  if (options.detector == "MFT") {
    return runFor<MFT>(options);
  }
  if (options.detector == "TRKFT3") {
    return runFor<TRKFT3>(options);
  }
  if (options.detector == "TRKFT3TOF") {
    return runFor<TRKFT3TOF>(options);
  }
  std::cerr << "unknown detector '" << options.detector << "'\n";
  return 1;
}

} // namespace
} // namespace o2::gtrk::tools

int main(int argc, char** argv)
{
  try {
    return o2::gtrk::tools::run(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
