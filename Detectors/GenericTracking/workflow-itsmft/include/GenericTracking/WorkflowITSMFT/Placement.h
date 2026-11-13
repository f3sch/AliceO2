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

/// \file Placement.h
/// \brief The ITS's and MFT's ALPIDE chips as our sensors (Alice/Placement.h).
///
/// **The origin.** In the plane, the active matrix is not centred on the chip:
/// the readout edge is wider than the opposite one, which moves the matrix
/// centre along local x by `(PassiveEdgeReadOut - PassiveEdgeTop) / 2`,
/// 0.58 mm. A sensor's origin is put there, so that the bounds of
/// `detectors::alpide::kSensor` are symmetric about it, and a cluster's loc0 is
/// ALICE's local x less `kActiveCentreX`..

#ifndef GTRK_WORKFLOWITSMFT_PLACEMENT_H_
#define GTRK_WORKFLOWITSMFT_PLACEMENT_H_

#include "GenericTracking/Alice/Placement.h"
#include "GenericTracking/Base/Sensor.h"

#include "MathUtils/Utils.h"

#include <cstddef>
#include <vector>

namespace o2::gtrk::alice
{

/// Where the active matrix's centre sits in ALICE's local x [cm].
extern const double kActiveCentreX;

/// Every chip of \a geometry as a flat sensor with the ALPIDE's active
/// matrix (the ITS's and MFT's), indexed by O2's chip ID, `layer`
/// `GeometryTGeo::getLayer`, `ring` 0.
template <typename GeometryTGeo>
[[nodiscard]] std::vector<Sensor> placeAlpides(GeometryTGeo& geometry)
{
  geometry.fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));
  std::vector<Sensor> sensors;
  sensors.reserve(static_cast<std::size_t>(geometry.getNumberOfChips()));
  for (int chip = 0; chip < geometry.getNumberOfChips(); ++chip) {
    sensors.push_back(flatSensor(geometry.getMatrixL2G(chip), geometry.getLayer(chip), kActiveCentreX));
  }
  return sensors;
}

} // namespace o2::gtrk::alice

#endif // GTRK_WORKFLOWITSMFT_PLACEMENT_H_
