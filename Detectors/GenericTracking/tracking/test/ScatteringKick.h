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

/// \file ScatteringKick.h
/// \brief The multiple-scattering kick the simulated chains of the fit tests
///        give a true state on a layer.
///
/// Theta and phi each move by a Gaussian draw times the plane angle
/// \a theta0, phi's divided by sin(theta) so that both are the same angle in
/// space; theta's draw comes first. The two forms differ only in which theta
/// that sine is of.

#ifndef GTRK_TEST_SCATTERINGKICK_H_
#define GTRK_TEST_SCATTERINGKICK_H_

#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/Definitions.h"

#include <cmath>
#include <random>

namespace o2::gtrk::test
{

/// The kick, phi's scaled by the sine of the theta before it.
inline void scatteringKick(BoundTrackState& state, Scalar theta0, std::normal_distribution<double>& gauss, std::mt19937& rng)
{
  const auto sinTheta = std::sin(state.theta());
  state.par[kTheta] += static_cast<Scalar>(gauss(rng)) * theta0;
  state.par[kPhi] += static_cast<Scalar>(gauss(rng)) * theta0 / sinTheta;
}

/// The kick, phi's scaled by the sine of the theta already kicked.
inline void scatteringKickKickedTheta(BoundTrackState& state, Scalar theta0, std::normal_distribution<double>& gauss, std::mt19937& rng)
{
  state.par[kTheta] += static_cast<Scalar>(gauss(rng)) * theta0;
  state.par[kPhi] += static_cast<Scalar>(gauss(rng)) * theta0 / std::sin(state.theta());
}

} // namespace o2::gtrk::test

#endif // GTRK_TEST_SCATTERINGKICK_H_
