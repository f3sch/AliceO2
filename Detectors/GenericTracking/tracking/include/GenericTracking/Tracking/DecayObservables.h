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

/// \file DecayObservables.h
/// \brief What a decay's fit gives a finder to select on: the invariant mass
///        of its daughters and its sigma, and with a reference point (a
///        primary vertex) the mother's flight from it and its direction.
///
/// Each quantity is a function of the fit's parameters (the vertices, the
/// daughters' phi, theta, q/p at theirs), differentiated forward (dual
/// numbers); its sigma is the gradient's quadratic form with the fit's
/// covariance, given as an accessor at(a, b) over the parameters the
/// quantity takes (so a packed result, a full matrix or a block of one
/// serve alike). A reference point's covariance adds to the flight's
/// variance, taken as independent of the decay (the point must not contain
/// the daughters).

#ifndef GTRK_TRACKING_DECAYOBSERVABLES_H_
#define GTRK_TRACKING_DECAYOBSERVABLES_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/DenseSolve.h"
#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Tracking/DecayFitConstraints.h"
#include "GenericTracking/Tracking/DecayFitMath.h"

#include "GPUCommonDef.h"

#include <array>
#include <span>

namespace o2::gtrk
{

namespace decay_fit_detail
{

/// u^T C v over \a n entries of u and v, C's entry (a, b) \a at(a, b): a
/// projection's variance, or two projections' covariance.
template <typename U, typename V, typename At>
[[nodiscard]] GPUhdi() double quadraticForm(const U& u, const V& v, int n, const At& at)
{
  double sum = 0.;
  for (int a = 0; a < n; ++a) {
    for (int b = 0; b < n; ++b) {
      sum += u[a] * at(a, b) * v[b];
    }
  }
  return sum;
}

} // namespace decay_fit_detail

/// The invariant mass of a set of daughters and its sigma, and their summed
/// |p| and its sigma.
struct DecayKinematics {
  double mass = 0.;
  double sMass = 0.;
  double p = 0.;
  double sP = 0.;
};

/// \a momenta's (phi, theta, q/p each) kinematics with \a masses and
/// \a charges; \a at(a, b) their covariance, momentum i's components at
/// 3 i + c.
template <typename At>
[[nodiscard]] GPUhdi() DecayKinematics decayKinematics(std::span<const Vec3D> momenta, std::span<const double> masses, std::span<const double> charges, const At& at)
{
  using namespace decay_fit_detail;
  using D = dual::Dual<3 * kMaxConstrainedDaughters>;
  const int count = static_cast<int>(momenta.size());
  MomentumSum<D> sum;
  for (int i = 0; i < count; ++i) {
    const auto momentum = momentumOfInputs<D::kInputs>(momenta[i][0], momenta[i][1], momenta[i][2], 3 * i, charges[i]);
    sum.add(momentum.vector, energyOf(momentum, masses[i]));
  }
  const D p = sqrtOf(sum.squared());
  const D m2 = sum.massSquared();
  const auto sigma = [&](const D& value) { return sqrtClampedOf(quadraticForm(value.d, value.d, 3 * count, at)); };
  const double mass = sqrtClampedOf(m2.v);
  return {.mass = mass, .sMass = mass > 0. ? sigma(m2) / (2. * mass) : 0., .p = p.v, .sP = sigma(p)};
}

/// The mother's flight from a reference point: its decay length along its
/// momentum, ct with its mass, its direction, the cosine of the angle
/// between its momentum and the line from the point (pointing), each with
/// its sigma where it has one.
struct DecayFlight {
  double length = 0.;
  double sLength = 0.;
  double ct = 0.;
  double sCt = 0.;
  double phi = 0.;
  double sPhi = 0.;
  double theta = 0.;
  double sTheta = 0.;
  double cosPointing = 0.;
};

/// The flight of the mother decaying at \a vertex into \a momenta (with
/// \a charges) from \a point, its mass \a mass; \a at(a, b) the fit's
/// covariance over its parameters, the vertices' first (\a vertices of
/// them, the mother's first), then the momenta; \a pointCovariance the
/// point's, packed (`covIndex`), zero for an exact point.
template <typename At>
[[nodiscard]] GPUhdi() DecayFlight decayFlight(const Vec3D& vertex, std::span<const Vec3D> momenta, std::span<const double> charges, const Vec3D& point, const std::array<double, 6>& pointCovariance, double mass, int vertices, const At& at)
{
  using namespace decay_fit_detail;
  using D = Dual<kConstraintInputs>;
  const int count = static_cast<int>(momenta.size());
  const int n = (3 * vertices) + (3 * count);
  std::array<D, 3> v{};
  for (int k = 0; k < 3; ++k) {
    v[k].v = vertex[k];
    v[k].d[k] = 1.;
  }
  MomentumSum<D> sum;
  for (int i = 0; i < count; ++i) {
    sum.add(momentumOfInputs<kConstraintInputs>(momenta[i][0], momenta[i][1], momenta[i][2], (3 * vertices) + (3 * i), charges[i]).vector);
  }
  const auto& total = sum.vector;
  const auto direction = directionOfMomentum(total);
  const D& p = direction.p;
  const D length = (((v[0] + (-point[0])) * total[0]) + ((v[1] + (-point[1])) * total[1]) + ((v[2] + (-point[2])) * total[2])) / p;
  const D ct = (mass * length) / p;
  const D& phi = direction.phi;
  const D& theta = direction.theta;
  // The point's share: d length / d point = -P / |P|, d ct = (mass / |P|) d length.
  const std::array<double, 3> dLength{-total[0].v / p.v, -total[1].v / p.v, -total[2].v / p.v};
  const double pointLength = quadraticForm(dLength, dLength, 3, [&](int a, int b) { return pointCovariance[covIndex(a, b)]; });
  const double pointCt = (mass / p.v) * (mass / p.v) * pointLength;
  const auto variance = [&](const D& value) { return quadraticForm(value.d, value.d, n, at); };
  const Vec3D flight{vertex[0] - point[0], vertex[1] - point[1], vertex[2] - point[2]};
  const double flightNorm = sqrtOf(dotOf(flight, flight));
  const double cosPointing = flightNorm > 0. ? ((flight[0] * total[0].v) + (flight[1] * total[1].v) + (flight[2] * total[2].v)) / (flightNorm * p.v) : 0.;
  return {.length = length.v, .sLength = sqrtClampedOf(variance(length) + pointLength), .ct = ct.v, .sCt = sqrtClampedOf(variance(ct) + pointCt), .phi = phi.v, .sPhi = sqrtClampedOf(variance(phi)), .theta = theta.v, .sTheta = sqrtClampedOf(variance(theta)), .cosPointing = cosPointing};
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_DECAYOBSERVABLES_H_
