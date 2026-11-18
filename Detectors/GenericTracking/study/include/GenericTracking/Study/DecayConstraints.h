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
/// The values and their derivatives by (V, the momenta) come from forward
/// differentiation. A constrained step solves the normal equations with the
/// constraints linearised, by Lagrange multipliers; its covariance is the
/// unconstrained one less what the constraints fix.

#ifndef GTRK_STUDY_DECAYCONSTRAINTS_H_
#define GTRK_STUDY_DECAYCONSTRAINTS_H_

#include "GenericTracking/Study/DecayFitMath.h"

#include <array>
#include <cmath>
#include <span>
#include <type_traits>
#include <vector>

namespace o2::gtrk::study
{

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
    return (mass ? 1 : 0) + (production ? 2 : 0) + composites + (2 * (topology.vertices() - 1)) + static_cast<int>(resonances.size());
  }
};

namespace decay_fit_detail
{

/// The most vertices and daughters the constraints take, and so their
/// derivatives' width.
inline constexpr int kMaxConstrainedVertices = 3;
inline constexpr int kMaxConstrainedDaughters = 5;
inline constexpr int kConstraintInputs = (3 * kMaxConstrainedVertices) + (3 * kMaxConstrainedDaughters);
/// The most equations `evaluateConstraints` gives: the mother's mass, and per
/// composite below the root its link (two) and its mass.
inline constexpr int kMaxHardRows = 1 + (3 * (kMaxConstrainedVertices - 1));

/// \a constraints' values at \a vertices (\a topology's) and \a momenta (phi,
/// theta, q/p each at its daughter's vertex; a neutral daughter's q/p its
/// 1/p) with the daughters' \a masses and \a charges, and their rows of
/// derivatives by (the vertices, the momenta), 3K + 3N wide: \a written
/// values into \a values, row j at \a rows + j \a stride, at most
/// \a capacity of them; false if the tree is too large or they are more.
/// Evaluated in the arithmetic M, stored in S.
template <typename M = double, typename S>
[[nodiscard]] bool evaluateConstraints(const DecayConstraints& constraints, const DecayTopology& topology, std::span<const std::array<double, 3>> vertices, std::span<const std::array<double, 3>> momenta, std::span<const double> masses, std::span<const double> charges, double fieldScale, S* values, S* rows, int stride, int capacity, int& written)
{
  using D = Dual<kConstraintInputs, M>;
  const int count = static_cast<int>(momenta.size());
  const int k = topology.vertices();
  const int width = (3 * k) + (3 * count);
  written = 0;
  if (count > kMaxConstrainedDaughters || k > kMaxConstrainedVertices) {
    return false;
  }
  std::array<std::array<D, 3>, kMaxConstrainedVertices> v{};
  for (int w = 0; w < k; ++w) {
    for (int c = 0; c < 3; ++c) {
      v[w][c].v = static_cast<M>(vertices[w][c]);
      v[w][c].d[(3 * w) + c] = M(1.);
    }
  }
  // Every daughter's momentum vector, energy and charge.
  std::array<std::array<D, 3>, kMaxConstrainedDaughters> vectors{};
  std::array<D, kMaxConstrainedDaughters> energies{};
  for (int i = 0; i < count; ++i) {
    const int q = (3 * k) + (3 * i);
    const auto momentum = momentumOf(dual::variable<kConstraintInputs, M>(momenta[i][0], q), dual::variable<kConstraintInputs, M>(momenta[i][1], q + 1), dual::variable<kConstraintInputs, M>(momenta[i][2], q + 2), numeratorOf(charges[i]));
    vectors[i] = momentum.vector;
    energies[i] = sqrt((momentum.magnitude * momentum.magnitude) + (masses[i] * masses[i]));
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
  bool room = true;
  const auto push = [&](const D& value) {
    if (written >= capacity) {
      room = false;
      return;
    }
    values[written] = static_cast<S>(value.v);
    for (int c = 0; c < width; ++c) {
      rows[(written * stride) + c] = static_cast<S>(value.d[c]);
    }
    ++written;
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
  return room;
}

/// Equations predicted - measured - L eta = 0 whose measurements carry a
/// covariance S = L L^T (`semidefiniteFactor`): eta, as many of the fit's
/// parameters from \a at, carry a unit Gaussian prior (`addPrior`); `apply`
/// subtracts L eta from the equations' values from \a row and puts -L in
/// their columns. Minimising eta^T eta under them is the chi2 of the
/// measurements with S, never inverting S: a measurement S leaves exact (a
/// zero row) holds exactly, one with a small S is no ill-conditioned weight.
template <typename S = double>
struct ScatteredEquations {
  int at = -1;  ///< eta's first parameter; -1: none
  int row = -1; ///< the first equation's row
  int count = 0;
  const S* factor = nullptr; ///< L, count x count (`semidefiniteFactor`: not triangular)
  void apply(std::span<const double> p, int n, S* values, S* rows) const
  {
    for (int a = 0; a < count; ++a) {
      for (int b = 0; b < count; ++b) {
        const S l = factor[(a * count) + b];
        values[row + a] -= l * static_cast<S>(p[at + b]);
        rows[((row + a) * n) + at + b] = -l;
      }
    }
  }
  void addPrior(std::span<const double> p, int n, double& chi2, S* normal, S* gradient) const
  {
    for (int b = 0; b < count; ++b) {
      chi2 += p[at + b] * p[at + b];
      if (normal != nullptr) {
        normal[((at + b) * n) + at + b] += S(1.);
        gradient[at + b] -= static_cast<S>(p[at + b]);
      }
    }
  }
};

/// The Lagrange system's storage, the caller's: the matrix ((n + m)^2), the
/// solution (n + m), the pivots (n + m).
template <typename S = double>
struct LagrangeSpace {
  S* system = nullptr;
  S* solution = nullptr;
  int* pivots = nullptr;
};

/// The Gauss-Newton step under the hard constraints C step + c = 0 (\a rows
/// C, m x n, \a values c) with the normal matrix \a normal (n x n) and its
/// right-hand side \a gradient: the Lagrange system [[N, C^T], [C, 0]]
/// [step; lambda] = [g; -c], solved by LU with partial pivoting, as N alone
/// may be singular where the constraints fix its null direction (two parallel
/// daughters along the flight). With \a step, the step; with \a covariance,
/// the constrained covariance: the system's inverse's N block (n x n); with
/// \a multipliers, lambda; with \a logDeterminant, ln |det| of the system.
/// False if the system is singular (an exactly zero pivot).
template <typename S>
[[nodiscard]] bool constrainedStep(const LagrangeSpace<S>& space, const S* normal, int n, const S* gradient, const S* rows, std::span<const S> values, std::type_identity_t<S>* step, std::type_identity_t<S>* covariance, std::type_identity_t<S>* multipliers = nullptr, std::type_identity_t<S>* logDeterminant = nullptr)
{
  const int m = static_cast<int>(values.size());
  const int size = n + m;
  S* system = space.system;
  std::fill_n(system, size * size, S(0.));
  for (int a = 0; a < n; ++a) {
    for (int b = 0; b < n; ++b) {
      system[(a * size) + b] = normal[(a * n) + b];
    }
  }
  for (int j = 0; j < m; ++j) {
    for (int b = 0; b < n; ++b) {
      system[((n + j) * size) + b] = rows[(j * n) + b];
      system[(b * size) + n + j] = rows[(j * n) + b];
    }
  }
  int* pivots = space.pivots;
  const auto lu = viewOf(system, size);
  if (!dense::luFactor(lu, size, pivots)) {
    return false;
  }
  S* solution = space.solution;
  for (int a = 0; a < n; ++a) {
    solution[a] = gradient[a];
  }
  for (int j = 0; j < m; ++j) {
    solution[n + j] = -values[j];
  }
  dense::luSolve(lu, size, pivots, solution);
  if (step != nullptr) {
    std::copy_n(solution, n, step);
  }
  if (multipliers != nullptr) {
    std::copy_n(solution + n, m, multipliers);
  }
  if (logDeterminant != nullptr) {
    *logDeterminant = dense::luLogAbsDeterminant(lu, size);
  }
  if (covariance != nullptr) {
    // The system's inverse's N block, a column at a time.
    for (int b = 0; b < n; ++b) {
      std::fill_n(solution, size, S(0.));
      solution[b] = S(1.);
      dense::luSolve(lu, size, pivots, solution);
      for (int a = 0; a < n; ++a) {
        covariance[(a * n) + b] = solution[a];
      }
    }
  }
  return true;
}

} // namespace decay_fit_detail

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYCONSTRAINTS_H_
