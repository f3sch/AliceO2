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
#ifndef GTRK_TEST_PROPAGATIONPRECISIONFIXTURE_H_
#define GTRK_TEST_PROPAGATIONPRECISIONFIXTURE_H_
#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Tracking/Propagator.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
namespace o2::gtrk::propagation_precision
{

inline TrackingSurface plane(double x, double y, double z, double angle, bool disc, SensorID id)
{
  const Scalar origin[3]{static_cast<Scalar>(x), static_cast<Scalar>(y), static_cast<Scalar>(z)};
  const auto sinAngle = static_cast<Scalar>(std::sin(angle));
  const auto cosAngle = static_cast<Scalar>(std::cos(angle));
  TrackingSurface result = disc ? discPlane(origin, sinAngle, cosAngle) : radialPlane(origin, sinAngle, cosAngle);
  result.id = id;
  result.bounds = {.halfX = 1.e5F, .halfY = 1.e5F};
  return result;
}

/// Slide \a surface's origin by \a offset0 and \a offset1 along its own
/// axes: the same plane, its centre elsewhere.
inline void slideOrigin(TrackingSurface& surface, double offset0, double offset1)
{
  for (std::int32_t i = 0; i < 3; ++i) {
    surface.origin[i] += static_cast<Scalar>(offset0 * static_cast<double>(surface.axis0[i]) + offset1 * static_cast<double>(surface.axis1[i]));
  }
}

/// How far the stored answer is from the helix it claims to be on.
struct Closure {
  double position = 0.; ///< stored local position against the helix point [cm]
  double plane = 0.;    ///< that point's distance from the target plane [cm]
  double azimuth = 0.;  ///< stored azimuth against the accumulated turn [rad]
};

/// Re-evaluate the helix in double at the path length the propagation
/// reported; its float rounding is 1e-7 relative, 30 nm over 3 m.
inline Closure closure(const Propagator& prop, const BoundTrackState& before, const BoundTrackState& after, const TrackingSurface& from, const TrackingSurface& to, Scalar reportedPath)
{
  const auto path = static_cast<double>(reportedPath);
  const auto phi = static_cast<double>(before.phi());
  const auto theta = static_cast<double>(before.theta());
  // The turn per unit path carries no theta: q/pT's sin(theta) cancels against
  // the transverse part of the step.
  const double turn = static_cast<double>(before.qOverP()) * static_cast<double>(prop.bzKiloGauss) * static_cast<double>(constants::math::B2C) * path;
  const double transverse = path * std::sin(theta);
  const double along = transverse * sinc(turn);
  const double across = transverse * versine(turn);
  const double step[3] = {along * std::cos(phi) - across * std::sin(phi), along * std::sin(phi) + across * std::cos(phi), path * std::cos(theta)};

  double point[3];
  surfaceToGlobal(from, before.loc0(), before.loc1(), point);
  for (std::int32_t i = 0; i < 3; ++i) {
    point[i] += step[i];
  }
  double local[3];
  surfaceToLocal(to, point, local);

  Closure out;
  out.position = std::hypot(local[0] - static_cast<double>(after.loc0()), local[1] - static_cast<double>(after.loc1()));
  out.plane = std::abs(local[2]);
  out.azimuth = std::abs(toPMPiOf(phi + turn - static_cast<double>(after.phi())));
  return out;
}

struct Worst {
  Closure value;
  std::int32_t checked = 0;

  void add(const Closure& next)
  {
    value.position = std::max(value.position, next.position);
    value.plane = std::max(value.plane, next.plane);
    value.azimuth = std::max(value.azimuth, next.azimuth);
    ++checked;
  }

  void report(const char* label) const { std::cout << label << ": checked=" << checked << " position um=" << value.position / kUm << " plane um=" << value.plane / kUm << " azimuth urad=" << value.azimuth * 1.e6 << '\n'; }
};

} // namespace o2::gtrk::propagation_precision
#endif
