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

/// \file DecayConstraints.h
/// \brief Hard constraints on a decay vertex fit's parameters -- the vertex V
///        and each daughter's phi, theta, q/p at it -- and the constrained
///        Gauss-Newton step both fits take (DecayVertexFit.h, DecayRowsFit.h).
///
/// - The mass: the daughters' summed four-momentum at V has the mother's mass,
///   m^2 - M^2 = 0, each daughter's energy from its species' mass.
/// - The production point: the mother -- the daughters' summed momentum at V,
///   their summed charge -- came from a known point C. Its trajectory from V,
///   a helix in the uniform field or a line if neutral, crosses the perigee
///   plane about C (through C's transverse point along the beam, normal to the
///   mother's direction at its closest approach to it) at C: loc0 = 0 and
///   loc1 = C_z, two equations (`predictOnPlane`). The plane is the current
///   iterate's, as a measurement's.
///
/// - Collinear daughters (a photon's conversion): the two daughters' momenta
///   parallel at V, phi1 - phi2 = 0 and theta1 - theta2 = 0. A pair has a
///   mass of at least 2 m_e, so a mass constraint cannot hold one; parallel
///   momenta is the approximation (DecayTreeFitter's, CMS's colinearity).
///
/// - V's position along a line, s = (V - origin) . direction: one equation;
///   a conversion's profile fixes its point along the flight with it
///   (DecayProfile.h).
///
/// The values and their derivatives by (V, the momenta) come from forward
/// differentiation. A constrained step solves the normal equations with the
/// constraints linearised, by Lagrange multipliers; its covariance is the
/// unconstrained one less what the constraints fix.

#ifndef GTRK_STUDY_DECAYCONSTRAINTS_H_
#define GTRK_STUDY_DECAYCONSTRAINTS_H_

#include "GenericTracking/Study/DecayFitMath.h"

#include <array>
#include <cmath>
#include <numbers>
#include <optional>
#include <span>
#include <vector>

namespace o2::gtrk::study
{

/// V at \a s along a line: (V - origin) . direction = s.
struct LinePosition {
  std::array<double, 3> origin{};
  std::array<double, 3> direction{}; ///< unit
  double s = 0.;                     ///< [cm]
};

/// The constraints a decay vertex fit takes.
struct DecayConstraints {
  bool mass = false;
  double massGeV = 0.; ///< the mother's
  bool production = false;
  std::array<double, 3> point{}; ///< where the mother was made [cm]
  bool collinear = false;        ///< two daughters only
  std::optional<LinePosition> position;
  [[nodiscard]] int count() const { return (mass ? 1 : 0) + (production ? 2 : 0) + (collinear ? 2 : 0) + (position ? 1 : 0); }
};

namespace decay_fit_detail
{

/// The most daughters the constraints take, and so their derivatives' width.
inline constexpr int kMaxConstrainedDaughters = 5;
inline constexpr int kConstraintInputs = 3 + (3 * kMaxConstrainedDaughters);

template <int N>
Dual<N> atan2(const Dual<N>& y, const Dual<N>& x)
{
  Dual<N> r{std::atan2(y.v, x.v)};
  const double norm = (x.v * x.v) + (y.v * y.v);
  for (int i = 0; i < N; ++i) {
    r.d[i] = ((x.v * y.d[i]) - (y.v * x.d[i])) / norm;
  }
  return r;
}

/// \a constraints' values at \a vertex and \a momenta (phi, theta, q/p each;
/// unit charges) with the daughters' \a masses, and their rows of derivatives
/// by (V, the momenta), 3 + 3N wide; false if the mother misses its plane or
/// there are too many daughters.
[[nodiscard]] inline bool evaluateConstraints(const DecayConstraints& constraints, const std::array<double, 3>& vertex, std::span<const std::array<double, 3>> momenta, std::span<const double> masses, double fieldScale, std::vector<double>& values, std::vector<double>& rows)
{
  using D = Dual<kConstraintInputs>;
  const int count = static_cast<int>(momenta.size());
  const int width = 3 + (3 * count);
  if (count > kMaxConstrainedDaughters || (constraints.collinear && count != 2)) {
    return false;
  }
  values.clear();
  rows.clear();
  std::array<D, 3> v{};
  for (int k = 0; k < 3; ++k) {
    v[k].v = vertex[k];
    v[k].d[k] = 1.;
  }
  std::array<D, 3> total{};
  D energy{};
  double charge = 0.;
  for (int i = 0; i < count; ++i) {
    std::array<D, 3> q{};
    for (int k = 0; k < 3; ++k) {
      q[k].v = momenta[i][k];
      q[k].d[3 + (3 * i) + k] = 1.;
    }
    const double sign = q[2].v > 0. ? 1. : -1.;
    charge += sign;
    const D p = D{sign} / q[2]; // |p| = 1 / |q/p|
    const D pt = p * sin(q[1]);
    const std::array<D, 3> vector{pt * cos(q[0]), pt * sin(q[0]), p * cos(q[1])};
    for (int k = 0; k < 3; ++k) {
      total[k] = total[k] + vector[k];
    }
    energy = energy + sqrt((p * p) + (masses[i] * masses[i]));
  }
  const auto push = [&](const D& value) {
    values.push_back(value.v);
    for (int c = 0; c < width; ++c) {
      rows.push_back(value.d[c]);
    }
  };
  if (constraints.mass) {
    const D p2 = (total[0] * total[0]) + (total[1] * total[1]) + (total[2] * total[2]);
    push((energy * energy) - p2 + (-(constraints.massGeV * constraints.massGeV)));
  }
  if (constraints.production) {
    const D pt = sqrt((total[0] * total[0]) + (total[1] * total[1]));
    const D p = sqrt((pt * pt) + (total[2] * total[2]));
    const D phi = atan2(total[1], total[0]);
    const D theta = atan2(pt, total[2]);
    const D qOverP = D{charge} / p;
    // The plane about C: normal to the mother's direction at its closest
    // approach to C, from the current values (`perigeePhi`).
    const double curvature = fieldScale * charge / pt.v;
    const double dx = vertex[0] - constraints.point[0];
    const double dy = vertex[1] - constraints.point[1];
    const double along = (dx * std::cos(phi.v)) + (dy * std::sin(phi.v));
    const double across = (-dx * std::sin(phi.v)) + (dy * std::cos(phi.v));
    const PerigeePlane plane{.x = constraints.point[0], .y = constraints.point[1], .phi = phi.v + std::atan2(-curvature * along, 1. + (curvature * across))};
    std::array<D, 5> crossing{};
    if (!predictOnPlane(v[0], v[1], v[2], phi, theta, qOverP, fieldScale, plane, crossing)) {
      return false;
    }
    push(crossing[0]);
    push(crossing[1] + (-constraints.point[2]));
  }
  if (constraints.collinear) {
    D dPhi{};
    dPhi.v = std::remainder(momenta[0][0] - momenta[1][0], 2. * std::numbers::pi);
    dPhi.d[3] = 1.;
    dPhi.d[6] = -1.;
    D dTheta{};
    dTheta.v = momenta[0][1] - momenta[1][1];
    dTheta.d[4] = 1.;
    dTheta.d[7] = -1.;
    push(dPhi);
    push(dTheta);
  }
  if (constraints.position) {
    const auto& line = *constraints.position;
    D along{};
    for (int k = 0; k < 3; ++k) {
      along.v += (vertex[k] - line.origin[k]) * line.direction[k];
      along.d[k] = line.direction[k];
    }
    along.v -= line.s;
    push(along);
  }
  return true;
}

/// The Gauss-Newton step under the hard constraints C step + c = 0 (\a rows
/// C, m x n, \a values c) with the normal matrix \a normal (n x n) and its
/// right-hand side \a gradient: the Lagrange system [[N, C^T], [C, 0]]
/// [step; lambda] = [g; -c], solved by LU with partial pivoting, as N alone
/// may be singular where the constraints fix its null direction (two parallel
/// daughters along the flight). With \a covariance, the constrained
/// covariance: the system's inverse's N block; with \a sensitivity (n x m),
/// the parameters' derivatives by each constraint's value, the inverse's
/// C^T block. False if the system is singular (an exactly zero pivot).
[[nodiscard]] inline bool constrainedStep(const std::vector<double>& normal, int n, const std::vector<double>& gradient, const std::vector<double>& rows, const std::vector<double>& values, std::vector<double>& step, std::vector<double>* covariance, std::vector<double>* sensitivity = nullptr)
{
  const int m = static_cast<int>(values.size());
  const int size = n + m;
  std::vector<double> system(static_cast<std::size_t>(size * size), 0.);
  for (int a = 0; a < n; ++a) {
    for (int b = 0; b < n; ++b) {
      system[a * size + b] = normal[a * n + b];
    }
  }
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      system[(n + j) * size + i] = rows[j * n + i];
      system[i * size + (n + j)] = rows[j * n + i];
    }
  }
  // LU with partial pivoting, in place.
  std::vector<int> pivot(static_cast<std::size_t>(size));
  for (int k = 0; k < size; ++k) {
    int best = k;
    for (int i = k + 1; i < size; ++i) {
      if (std::abs(system[i * size + k]) > std::abs(system[best * size + k])) {
        best = i;
      }
    }
    if (system[best * size + k] == 0.) {
      return false;
    }
    pivot[k] = best;
    if (best != k) {
      for (int j = 0; j < size; ++j) {
        std::swap(system[k * size + j], system[best * size + j]);
      }
    }
    for (int i = k + 1; i < size; ++i) {
      const double factor = system[i * size + k] / system[k * size + k];
      system[i * size + k] = factor;
      for (int j = k + 1; j < size; ++j) {
        system[i * size + j] -= factor * system[k * size + j];
      }
    }
  }
  const auto solve = [&](std::vector<double>& b) {
    for (int k = 0; k < size; ++k) {
      std::swap(b[k], b[pivot[k]]);
    }
    for (int i = 0; i < size; ++i) {
      for (int k = 0; k < i; ++k) {
        b[i] -= system[i * size + k] * b[k];
      }
    }
    for (int i = size - 1; i >= 0; --i) {
      for (int k = i + 1; k < size; ++k) {
        b[i] -= system[i * size + k] * b[k];
      }
      b[i] /= system[i * size + i];
    }
  };
  std::vector<double> rhs(static_cast<std::size_t>(size), 0.);
  std::copy(gradient.begin(), gradient.end(), rhs.begin());
  for (int j = 0; j < m; ++j) {
    rhs[n + j] = -values[j];
  }
  solve(rhs);
  step.assign(rhs.begin(), rhs.begin() + n);
  if (covariance != nullptr) {
    covariance->assign(static_cast<std::size_t>(n * n), 0.);
    for (int c = 0; c < n; ++c) {
      std::vector<double> unit(static_cast<std::size_t>(size), 0.);
      unit[c] = 1.;
      solve(unit);
      for (int r = 0; r < n; ++r) {
        (*covariance)[r * n + c] = unit[r];
      }
    }
  }
  if (sensitivity != nullptr) {
    sensitivity->assign(static_cast<std::size_t>(n * m), 0.);
    for (int j = 0; j < m; ++j) {
      std::vector<double> unit(static_cast<std::size_t>(size), 0.);
      unit[n + j] = 1.;
      solve(unit);
      for (int r = 0; r < n; ++r) {
        (*sensitivity)[r * m + j] = unit[r];
      }
    }
  }
  return true;
}

} // namespace decay_fit_detail

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYCONSTRAINTS_H_
