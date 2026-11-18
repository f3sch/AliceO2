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
///        Gauss-Newton step the fit takes (DecayVertexFit.h).
///
/// - The mass: the daughters' summed four-momentum at V has the mother's mass,
///   m^2 - M^2 = 0, each daughter's energy from its species' mass.
/// - The production point: the mother -- the daughters' summed momentum at V,
///   their summed charge -- came from a known point C. Its trajectory from V,
///   a helix in the uniform field or a line if neutral, crosses the perigee
///   plane about C (through C's transverse point along the beam, normal to the
///   mother's direction at its closest approach to it) at C: loc0 = 0 and
///   loc1 = C_z, two equations. The plane is the current iterate's, as a
///   measurement's. The fit takes them with a tracked mother's own rows as
///   one chain walked from V (DecayVertexFit.h, `DecayMother`), the
///   scattering on the way their nuisances' (`ScatteredEquations`); they are
///   not among the values here.
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
#include <limits>
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

/// A decay tree's shape: which vertex each daughter leaves, and how the
/// vertices hang together. Vertex 0 is the root's decay point; every other
/// vertex v is the decay point of a composite made at its parent's vertex
/// (`parent[v] < v`), whose momentum is the sum of its subtree's daughters'.
/// Empty, one vertex for every daughter.
struct DecayTopology {
  std::vector<int> vertexOf; ///< per daughter
  std::vector<int> parent;   ///< per vertex, -1 for the root
  /// Per vertex, the decaying particle's mass [GeV]: its flight's material
  /// (DecayVertexFit.h's chains); empty or 0 where unknown.
  std::vector<double> massGeV;
  [[nodiscard]] int vertices() const { return parent.empty() ? 1 : static_cast<int>(parent.size()); }
  [[nodiscard]] int vertexOfDaughter(int daughter) const { return vertexOf.empty() ? 0 : vertexOf[daughter]; }
  /// Whether \a daughter's vertex is \a vertex or below it.
  [[nodiscard]] bool under(int daughter, int vertex) const
  {
    for (int v = vertexOfDaughter(daughter); v >= 0; v = parent.empty() ? -1 : parent[v]) {
      if (v == vertex) {
        return true;
      }
    }
    return false;
  }
};

/// A resonance among the daughters: it decays where it is made, so its
/// products sit at its parent's vertex; their invariant mass m is its own
/// by the relativistic Breit-Wigner: -2 ln BW = 2 ln(1 + u^2),
/// u = (m^2 - M^2) / (M Gamma) a standard Cauchy variable, solved by Fisher
/// scoring (its information 1/2 per unit u^2). (A Gaussian of the
/// Breit-Wigner's FWHM in m was measured worse: it pulls the tails' decays.)
struct Resonance {
  std::vector<int> daughters; ///< its products among the fit's daughters
  double massGeV = 0.;
  double widthGeV = 0.;
};

/// The constraints a decay vertex fit takes. Besides these, a tree's links
/// always hold: every composite below the root flies from its parent's
/// vertex to its own along its momentum, two equations -- a neutral one in a
/// line (here), a charged one on its helix with its scattering
/// (DecayVertexFit.h's chains).
struct DecayConstraints {
  bool mass = false;
  double massGeV = 0.; ///< the mother's
  bool production = false;
  std::array<double, 3> point{}; ///< where the mother was made [cm]
  bool collinear = false;        ///< two daughters only
  std::optional<LinePosition> position;
  /// Per vertex below the root, its composite's mass to hold [GeV]; 0 for none.
  std::vector<double> compositeMasses;
  /// Resonances' masses, measured softly (a term of the chi2 each).
  std::vector<Resonance> resonances;
  /// What the constraints add to the degrees of freedom.
  [[nodiscard]] int count(const DecayTopology& topology = {}) const
  {
    int composites = 0;
    for (const auto m : compositeMasses) {
      composites += m > 0. ? 1 : 0;
    }
    return (mass ? 1 : 0) + (production ? 2 : 0) + (collinear ? 2 : 0) + (position ? 1 : 0) + composites + (2 * (topology.vertices() - 1)) + static_cast<int>(resonances.size());
  }
};

namespace decay_fit_detail
{

/// The most vertices and daughters the constraints take, and so their
/// derivatives' width.
inline constexpr int kMaxConstrainedVertices = 3;
inline constexpr int kMaxConstrainedDaughters = 5;
inline constexpr int kConstraintInputs = (3 * kMaxConstrainedVertices) + (3 * kMaxConstrainedDaughters);

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

/// \a constraints' values at \a vertices (\a topology's) and \a momenta (phi,
/// theta, q/p each at its daughter's vertex; a neutral daughter's q/p its
/// 1/p) with the daughters' \a masses and \a charges, and their rows of
/// derivatives by (the vertices, the
/// momenta), 3K + 3N wide; false if the tree is too large.
[[nodiscard]] inline bool evaluateConstraints(const DecayConstraints& constraints, const DecayTopology& topology, std::span<const std::array<double, 3>> vertices, std::span<const std::array<double, 3>> momenta, std::span<const double> masses, std::span<const double> charges, double fieldScale, std::vector<double>& values, std::vector<double>& rows)
{
  using D = Dual<kConstraintInputs>;
  const int count = static_cast<int>(momenta.size());
  const int k = topology.vertices();
  const int width = (3 * k) + (3 * count);
  if (count > kMaxConstrainedDaughters || k > kMaxConstrainedVertices || (constraints.collinear && count != 2)) {
    return false;
  }
  values.clear();
  rows.clear();
  std::vector<std::array<D, 3>> v(static_cast<std::size_t>(k));
  for (int w = 0; w < k; ++w) {
    for (int c = 0; c < 3; ++c) {
      v[w][c].v = vertices[w][c];
      v[w][c].d[(3 * w) + c] = 1.;
    }
  }
  // Every daughter's momentum vector, energy and charge.
  std::vector<std::array<D, 3>> vectors(static_cast<std::size_t>(count));
  std::vector<D> energies(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    std::array<D, 3> q{};
    for (int c = 0; c < 3; ++c) {
      q[c].v = momenta[i][c];
      q[c].d[(3 * k) + (3 * i) + c] = 1.;
    }
    const D p = D{q[2].v > 0. ? 1. : -1.} / q[2]; // |p| = 1 / |q/p|
    const D pt = p * sin(q[1]);
    vectors[i] = {pt * cos(q[0]), pt * sin(q[0]), p * cos(q[1])};
    energies[i] = sqrt((p * p) + (masses[i] * masses[i]));
  }
  // A vertex's composite: its subtree's momentum, energy, charge.
  const auto composite = [&](int vertex, std::array<D, 3>& total, D& energy, double& charge) {
    total = {};
    energy = D{};
    charge = 0.;
    for (int i = 0; i < count; ++i) {
      if (topology.under(i, vertex)) {
        for (int c = 0; c < 3; ++c) {
          total[c] = total[c] + vectors[i][c];
        }
        energy = energy + energies[i];
        charge += charges[i];
      }
    }
  };
  const auto push = [&](const D& value) {
    values.push_back(value.v);
    for (int c = 0; c < width; ++c) {
      rows.push_back(value.d[c]);
    }
  };
  const auto massOf = [&](const std::array<D, 3>& total, const D& energy, double m) {
    const D p2 = (total[0] * total[0]) + (total[1] * total[1]) + (total[2] * total[2]);
    return (energy * energy) - p2 + (-(m * m));
  };
  std::array<D, 3> total{};
  D energy{};
  double charge = 0.;
  composite(0, total, energy, charge);
  if (constraints.mass) {
    push(massOf(total, energy, constraints.massGeV));
  }
  // The links: each composite below the root flies from its parent's vertex
  // to its own along its momentum P: the offset's components across P, on
  // u1 = z x P / |z x P| (x where P is along z) and u2 = P x u1 / |P|.
  for (int w = 1; w < k; ++w) {
    composite(w, total, energy, charge);
    if (charge != 0.) {
      continue; // a charged composite's link: its chain's (DecayVertexFit.h)
    }
    const int up = topology.parent[w];
    std::array<D, 3> offset{};
    for (int c = 0; c < 3; ++c) {
      offset[c] = v[w][c] - v[up][c];
    }
    const D p = sqrt((total[0] * total[0]) + (total[1] * total[1]) + (total[2] * total[2]));
    std::array<D, 3> u1{-total[1], total[0], D{}};
    D across = sqrt((u1[0] * u1[0]) + (u1[1] * u1[1]));
    if (across.v == 0.) {
      u1 = {D{1.}, D{}, D{}};
      across = D{1.};
    }
    for (auto& c : u1) {
      c = c / across;
    }
    const std::array<D, 3> u2{((total[1] * u1[2]) - (total[2] * u1[1])) / p, ((total[2] * u1[0]) - (total[0] * u1[2])) / p, ((total[0] * u1[1]) - (total[1] * u1[0])) / p};
    push((offset[0] * u1[0]) + (offset[1] * u1[1]) + (offset[2] * u1[2]));
    push((offset[0] * u2[0]) + (offset[1] * u2[1]) + (offset[2] * u2[2]));
  }
  for (int w = 1; w < k && w - 1 < static_cast<int>(constraints.compositeMasses.size()); ++w) {
    if (constraints.compositeMasses[w - 1] > 0.) {
      composite(w, total, energy, charge);
      push(massOf(total, energy, constraints.compositeMasses[w - 1]));
    }
  }
  if (constraints.collinear) {
    D dPhi{};
    dPhi.v = std::remainder(momenta[0][0] - momenta[1][0], 2. * std::numbers::pi);
    dPhi.d[(3 * k)] = 1.;
    dPhi.d[(3 * k) + 3] = -1.;
    D dTheta{};
    dTheta.v = momenta[0][1] - momenta[1][1];
    dTheta.d[(3 * k) + 1] = 1.;
    dTheta.d[(3 * k) + 4] = -1.;
    push(dPhi);
    push(dTheta);
  }
  if (constraints.position) {
    const auto& line = *constraints.position;
    D along{};
    for (int c = 0; c < 3; ++c) {
      along.v += (vertices[0][c] - line.origin[c]) * line.direction[c];
      along.d[c] = line.direction[c];
    }
    along.v -= line.s;
    push(along);
  }
  return true;
}

/// Equations predicted - measured - L eta = 0 whose measurements carry a
/// covariance S = L L^T (`semidefiniteFactor`): eta, as many of the fit's
/// parameters from \a at, carry a unit Gaussian prior (`addPrior`); `apply`
/// subtracts L eta from the equations' values from \a row and puts -L in
/// their columns. Minimising eta^T eta under them is the chi2 of the
/// measurements with S, never inverting S: a measurement S leaves exact (a
/// zero row) holds exactly, one with a small S is no ill-conditioned weight.
struct ScatteredEquations {
  int at = -1;  ///< eta's first parameter; -1: none
  int row = -1; ///< the first equation's row
  int count = 0;
  std::vector<double> factor; ///< L, count x count, lower triangular
  void apply(const std::vector<double>& p, int n, std::vector<double>& values, std::vector<double>& rows) const
  {
    for (int a = 0; a < count; ++a) {
      for (int b = 0; b <= a; ++b) {
        const double l = factor[(a * count) + b];
        values[row + a] -= l * p[at + b];
        rows[((row + a) * n) + at + b] = -l;
      }
    }
  }
  void addPrior(const std::vector<double>& p, int n, double& chi2, std::vector<double>* normal, std::vector<double>* gradient) const
  {
    for (int b = 0; b < count; ++b) {
      chi2 += p[at + b] * p[at + b];
      if (normal != nullptr) {
        (*normal)[((at + b) * n) + at + b] += 1.;
        (*gradient)[at + b] -= p[at + b];
      }
    }
  }
};

/// The Gauss-Newton step under the hard constraints C step + c = 0 (\a rows
/// C, m x n, \a values c) with the normal matrix \a normal (n x n) and its
/// right-hand side \a gradient: the Lagrange system [[N, C^T], [C, 0]]
/// [step; lambda] = [g; -c], solved by LU with partial pivoting, as N alone
/// may be singular where the constraints fix its null direction (two parallel
/// daughters along the flight). With \a covariance, the constrained
/// covariance: the system's inverse's N block; with \a sensitivity (n x m),
/// the parameters' derivatives by each constraint's value, the inverse's
/// C^T block; with \a logDeterminant, ln |det| of the system. False if the
/// system is singular (an exactly zero pivot).
[[nodiscard]] inline bool constrainedStep(const std::vector<double>& normal, int n, const std::vector<double>& gradient, const std::vector<double>& rows, const std::vector<double>& values, std::vector<double>& step, std::vector<double>* covariance, std::vector<double>* sensitivity = nullptr, std::vector<double>* multipliers = nullptr, double* logDeterminant = nullptr)
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
  if (logDeterminant != nullptr) {
    *logDeterminant = 0.;
    for (int k = 0; k < size; ++k) {
      *logDeterminant += std::log(std::abs(system[k * size + k]));
    }
  }
  solve(rhs);
  step.assign(rhs.begin(), rhs.begin() + n);
  if (multipliers != nullptr) {
    multipliers->assign(rhs.begin() + n, rhs.end());
  }
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
