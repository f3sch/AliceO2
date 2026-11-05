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

/// \file Perigee.h
/// \brief The direction of a uniform-solenoid helix at its closest approach
///        to the beam line.

#ifndef GTRK_TRACKING_PERIGEE_H_
#define GTRK_TRACKING_PERIGEE_H_

#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/Trig.h"

#include "MathUtils/detail/trigonometric.h"

namespace o2::gtrk
{

/// The direction [rad] of a bound state's helix at its transverse closest
/// approach to the z axis, without a virtual beam-line surface.
///
/// With A = r.dot(direction) and B = r.dot(leftNormal), closest approach obeys
///
///   kappa A cos(psi) + (1 + kappa B) sin(psi) = 0,
///
/// solved for the turn psi; at zero curvature it is 0.
template <typename Prop>
[[nodiscard]] GPUdi() Scalar perigeePhi(const Prop& prop, const BoundTrackState& state, const TrackingSurface& surface)
{
  ExactScalar position[3]{};
  detail::boundToGlobal(surface, static_cast<ExactScalar>(state.loc0()), static_cast<ExactScalar>(state.loc1()), position);

  const auto phi = static_cast<ExactScalar>(state.phi());
  const auto theta = static_cast<ExactScalar>(state.theta());
  ExactScalar sinPhi = 0.;
  ExactScalar cosPhi = 0.;
  ExactScalar sinTheta = 0.;
  ExactScalar cosTheta = 0.;
  sincosOf(phi, sinPhi, cosPhi);
  sincosOf(theta, sinTheta, cosTheta);

  const ExactScalar curvature = pointPropagator(prop, surface, state).curvature(static_cast<ExactScalar>(state.qOverP()), sinTheta);
  const ExactScalar along = position[0] * cosPhi + position[1] * sinPhi;
  const ExactScalar across = -position[0] * sinPhi + position[1] * cosPhi;
  const ExactScalar numerator = -curvature * along;
  const ExactScalar denominator = 1. + curvature * across;
  const ExactScalar turn = atan2Of(numerator, denominator);
  ExactScalar phi0 = phi + turn;
  math_utils::detail::bringToPMPi<ExactScalar>(phi0);
  return static_cast<Scalar>(phi0);
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_PERIGEE_H_
