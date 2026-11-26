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

/// \file DecayFitConstraints.h
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

#ifndef GTRK_TRACKING_DECAYFITCONSTRAINTS_H_
#define GTRK_TRACKING_DECAYFITCONSTRAINTS_H_

#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Base/FixedVector.h"
#include "GenericTracking/Tracking/DecayFitMath.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>

namespace o2::gtrk
{

namespace decay_fit_detail
{
/// The most vertices and daughters the constraints take, and so their
/// derivatives' width.
inline constexpr std::int32_t kMaxConstrainedVertices = 3;
inline constexpr std::int32_t kMaxConstrainedDaughters = 5;
/// The most resonances among the daughters: each has two products or more,
/// so a tree of them over the daughters has at most one fewer.
inline constexpr std::int32_t kMaxResonances = kMaxConstrainedDaughters - 1;
} // namespace decay_fit_detail

/// A decay tree's shape: which vertex each daughter leaves, and how the
/// vertices hang together. Vertex 0 is the root's decay point; every other
/// vertex v is the decay point of a composite made at its parent's vertex
/// (`parent[v] < v`), whose momentum is the sum of its subtree's daughters'.
/// Empty, one vertex for every daughter.
struct DecayTopology {
  FixedVector<std::int32_t, decay_fit_detail::kMaxConstrainedDaughters> vertexOf; ///< per daughter
  FixedVector<std::int32_t, decay_fit_detail::kMaxConstrainedVertices> parent;    ///< per vertex, -1 for the root
  /// Per vertex, the decaying particle's mass [GeV]: its flight's material
  /// (DecayVertexFit.h's chains); empty or 0 where unknown.
  FixedVector<double, decay_fit_detail::kMaxConstrainedVertices> massGeV;
  [[nodiscard]] GPUhdi() std::uint32_t vertices() const { return parent.empty() ? 1U : parent.size(); }
  [[nodiscard]] GPUhdi() std::int32_t vertexOfDaughter(std::uint32_t daughter) const { return vertexOf.empty() ? 0 : vertexOf[daughter]; }
  /// Whether \a daughter's vertex is \a vertex or below it.
  [[nodiscard]] GPUhdi() bool under(std::uint32_t daughter, std::uint32_t vertex) const
  {
    for (std::int32_t v = vertexOfDaughter(daughter); v >= 0; v = parent.empty() ? -1 : parent[static_cast<std::uint32_t>(v)]) {
      if (std::cmp_equal(v, vertex)) {
        return true;
      }
    }
    return false;
  }
  /// The summed charge of the daughters (their \a charges) at \a vertex or
  /// below it.
  template <typename Charges>
  [[nodiscard]] GPUhdi() double chargeUnder(const Charges& charges, std::uint32_t vertex) const
  {
    double sum = 0.;
    for (std::uint32_t i = 0; i < charges.size(); ++i) {
      sum += under(i, vertex) ? charges[i] : 0.;
    }
    return sum;
  }
};

/// A resonance among the daughters: it decays where it is made, so its
/// products sit at its parent's vertex; their invariant mass m is its own
/// by the relativistic Breit-Wigner: -2 ln BW = 2 ln(1 + u^2),
/// u = (m^2 - M^2) / (M Gamma) a standard Cauchy variable, solved by Fisher
/// scoring (its information 1/2 per unit u^2). The Breit-Wigner, not a
/// Gaussian in m: a Gaussian's quadratic penalty pulls the decays in the
/// tails, which the Cauchy's logarithmic one leaves.
struct Resonance {
  FixedVector<std::int32_t, decay_fit_detail::kMaxConstrainedDaughters> daughters; ///< its products among the fit's daughters
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
  Vec3D point{}; ///< where the mother was made [cm]
  /// The point's own covariance (x, y, z, packed: `covIndex`) [cm^2]: a
  /// fitted primary vertex's; zero for a point known exactly.
  std::array<double, 6> pointCovariance{};
  /// Per vertex below the root, its composite's mass to hold [GeV]; 0 for none.
  FixedVector<double, decay_fit_detail::kMaxConstrainedVertices> compositeMasses;
  /// Resonances' masses, measured softly (a term of the chi2 each).
  FixedVector<Resonance, decay_fit_detail::kMaxResonances> resonances;
  /// The composites whose mass is held.
  [[nodiscard]] GPUhdi() std::uint32_t heldMasses() const
  {
    std::uint32_t held = 0;
    for (const auto m : compositeMasses) {
      held += m > 0. ? 1U : 0U;
    }
    return held;
  }
  /// What the constraints add to the degrees of freedom: besides the
  /// equations, the links' two each (a charged composite's through its
  /// chain) and the resonances' terms.
  [[nodiscard]] GPUhdi() std::uint32_t count(const DecayTopology& topology = {}) const
  {
    return (mass ? 1U : 0U) + (production ? 2U : 0U) + heldMasses() + (2 * (topology.vertices() - 1)) + resonances.size();
  }
};

namespace decay_fit_detail
{

inline constexpr std::int32_t kConstraintInputs = (3 * kMaxConstrainedVertices) + (3 * kMaxConstrainedDaughters);
/// The most equations `evaluateConstraints` gives: the mother's mass, and per
/// composite below the root its link (two) and its mass.
inline constexpr std::int32_t kMaxHardRows = 1 + (3 * (kMaxConstrainedVertices - 1));

/// Equations written as values and rows of derivatives: equation j's value
/// into \a values[j], its derivatives by the first \a width parameters
/// (the vertices', the momenta's: 3K + 3N) into row j at \a rows + j
/// \a stride; at most \a capacity of them, beyond which `room` turns false.
template <typename S>
struct EquationWriter {
  S* values = nullptr;
  S* rows = nullptr;
  std::uint32_t stride = 0;
  std::uint32_t width = 0;
  std::uint32_t capacity = 0;
  std::uint32_t written = 0;
  bool room = true;
  template <std::int32_t N, typename M>
  GPUhdi() void push(const Dual<N, M>& value)
  {
    if (written >= capacity) {
      room = false;
      return;
    }
    values[written] = static_cast<S>(value.v);
    for (std::uint32_t c = 0; c < width; ++c) {
      rows[(written * stride) + c] = static_cast<S>(value.d[c]);
    }
    ++written;
  }
};

/// \a constraints' values at \a vertices (\a topology's) and \a momenta (phi,
/// theta, q/p each at its daughter's vertex; a neutral daughter's q/p its
/// 1/p) with the daughters' \a masses and \a charges, and their rows of
/// derivatives by (the vertices, the momenta), after \a out's equations;
/// false if the tree is too large or \a out has no room for them.
/// Evaluated in the arithmetic M, stored in S.
template <typename M = double, typename S>
[[nodiscard]] GPUhdi() bool evaluateConstraints(const DecayConstraints& constraints, const DecayTopology& topology, std::span<const Vec3D> vertices, std::span<const Vec3D> momenta, std::span<const double> masses, std::span<const double> charges, EquationWriter<S>& out)
{
  using D = Dual<kConstraintInputs, M>;
  const auto count = static_cast<std::uint32_t>(momenta.size());
  const std::uint32_t k = topology.vertices();
  if (count > kMaxConstrainedDaughters || k > kMaxConstrainedVertices) {
    return false;
  }
  std::array<std::array<D, 3>, kMaxConstrainedVertices> v{};
  for (std::uint32_t w = 0; w < k; ++w) {
    v[w] = vertexOfInputs<kConstraintInputs, M>(vertices[w], w);
  }
  // Every daughter's momentum vector and energy.
  std::array<std::array<D, 3>, kMaxConstrainedDaughters> vectors{};
  std::array<D, kMaxConstrainedDaughters> energies{};
  for (std::uint32_t i = 0; i < count; ++i) {
    const auto momentum = momentumOfInputs<kConstraintInputs, M>(momenta[i][0], momenta[i][1], momenta[i][2], (3 * k) + (3 * i), charges[i]);
    vectors[i] = momentum.vector;
    energies[i] = energyOf(momentum, masses[i]);
  }
  // A vertex's composite: its subtree's momentum and energy.
  const auto composite = [&](std::uint32_t vertex) {
    MomentumSum<D> sum;
    for (std::uint32_t i = 0; i < count; ++i) {
      if (topology.under(i, vertex)) {
        sum.add(vectors[i], energies[i]);
      }
    }
    return sum;
  };
  const auto massOf = [&](const MomentumSum<D>& sum, double m) { return sum.massSquared() + static_cast<M>(-(m * m)); };
  if (constraints.mass) {
    out.push(massOf(composite(0), constraints.massGeV));
  }
  // The links: each composite below the root flies from its parent's vertex
  // to its own along its momentum P: the offset's components across P, on
  // u1 = z x P / |z x P| (x where P is along z) and u2 = P x u1 / |P|.
  for (std::uint32_t w = 1; w < k; ++w) {
    if (!exactlyEqual(topology.chargeUnder(charges, w), 0.)) {
      continue; // a charged composite's link: its chain's (DecayVertexFit.h)
    }
    const auto sum = composite(w);
    const auto& total = sum.vector;
    const auto up = static_cast<std::uint32_t>(topology.parent[w]);
    std::array<D, 3> offset{};
    for (std::uint32_t c = 0; c < 3; ++c) {
      offset[c] = v[w][c] - v[up][c];
    }
    const D p = sqrtOf(sum.squared());
    std::array<D, 3> u1{-total[1], total[0], D{}};
    D across = sqrtOf((u1[0] * u1[0]) + (u1[1] * u1[1]));
    if (exactlyEqual(static_cast<double>(across.v), 0.)) {
      u1 = {D{1.}, D{}, D{}};
      across = D{1.};
    }
    for (auto& c : u1) {
      c = c / across;
    }
    std::array<D, 3> u2 = crossOf(total, u1);
    for (auto& c : u2) {
      c = c / p;
    }
    out.push(dotOf(offset, u1));
    out.push(dotOf(offset, u2));
  }
  for (std::uint32_t w = 1; w < k && w - 1 < constraints.compositeMasses.size(); ++w) {
    if (constraints.compositeMasses[w - 1] > 0.) {
      out.push(massOf(composite(w), constraints.compositeMasses[w - 1]));
    }
  }
  return out.room;
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
  std::uint32_t at = 0;  ///< eta's first parameter
  std::uint32_t row = 0; ///< the first equation's row
  std::uint32_t count = 0;
  const S* factor = nullptr; ///< L, count x count (`semidefiniteFactor`: not triangular)
  GPUhdi() void apply(std::span<const double> p, std::uint32_t n, S* values, S* rows) const
  {
    for (std::uint32_t a = 0; a < count; ++a) {
      for (std::uint32_t b = 0; b < count; ++b) {
        const S l = factor[(a * count) + b];
        values[row + a] -= l * static_cast<S>(p[at + b]);
        rows[((row + a) * n) + at + b] = -l;
      }
    }
  }
  /// \a normal packed (`dense::PackedView`).
  GPUhdi() void addPrior(std::span<const double> p, double& chi2, S* normal, S* gradient) const
  {
    for (std::uint32_t b = 0; b < count; ++b) {
      chi2 += p[at + b] * p[at + b];
      if (normal != nullptr) {
        normal[packedIndex(at + b, at + b)] += S(1.);
        gradient[at + b] -= static_cast<S>(p[at + b]);
      }
    }
  }
};

/// The Lagrange system's storage, the caller's: the matrix packed
/// (`covSize(n + m)`), the solution (n + m), the factorisation's
/// exchanges and block sizes (n + m each).
template <typename S = double>
struct LagrangeSpace {
  S* system = nullptr;
  S* solution = nullptr;
  std::int32_t* swaps = nullptr;
  std::int32_t* blocks = nullptr;
};

/// The Gauss-Newton step under the hard constraints C step + c = 0 (\a rows
/// C, m x n, \a values c) with the normal matrix \a normal (n x n, packed)
/// and its right-hand side \a gradient: the Lagrange system [[N, C^T], [C,
/// 0]] [step; lambda] = [g; -c], symmetric and indefinite, solved by L D L^T
/// with Bunch-Kaufman pivoting, as N alone may be singular where the
/// constraints fix its null direction (two parallel daughters along the
/// flight). With \a step, the step; with \a covariance, the constrained
/// covariance: the system's inverse's N block (packed); with \a multipliers,
/// lambda; with \a logDeterminant, ln |det| of the system. False if the
/// system is singular (an exactly singular pivot block).
template <typename S>
[[nodiscard]] GPUhdi() bool constrainedStep(const LagrangeSpace<S>& space, const S* normal, std::uint32_t n, const S* gradient, const S* rows, std::span<const S> values, std::type_identity_t<S>* step, std::type_identity_t<S>* covariance, std::type_identity_t<S>* multipliers = nullptr, std::type_identity_t<S>* logDeterminant = nullptr)
{
  const auto m = static_cast<std::uint32_t>(values.size());
  const std::uint32_t size = n + m;
  const dense::PackedView<S> system{space.system};
  for (std::uint32_t a = 0; a < n; ++a) {
    for (std::uint32_t b = 0; b <= a; ++b) {
      space.system[packedIndex(a, b)] = normal[packedIndex(a, b)];
    }
  }
  for (std::uint32_t j = 0; j < m; ++j) {
    for (std::uint32_t b = 0; b < n; ++b) {
      space.system[packedIndex(n + j, b)] = rows[(j * n) + b];
    }
    for (std::uint32_t k = 0; k <= j; ++k) {
      space.system[packedIndex(n + j, n + k)] = S(0.);
    }
  }
  if (!dense::ldltFactor(system, size, space.swaps, space.blocks)) {
    return false;
  }
  S* solution = space.solution;
  for (std::uint32_t a = 0; a < n; ++a) {
    solution[a] = gradient[a];
  }
  for (std::uint32_t j = 0; j < m; ++j) {
    solution[n + j] = -values[j];
  }
  dense::ldltSolve(system, size, space.swaps, space.blocks, solution);
  if (step != nullptr) {
    copyN(solution, n, step);
  }
  if (multipliers != nullptr) {
    copyN(solution + n, m, multipliers);
  }
  if (logDeterminant != nullptr) {
    *logDeterminant = dense::ldltLogAbsDeterminant(system, size, space.blocks);
  }
  if (covariance != nullptr) {
    // The system's inverse's N block, a column at a time, its lower triangle.
    for (std::uint32_t b = 0; b < n; ++b) {
      fillN(solution, size, S(0.));
      solution[b] = S(1.);
      dense::ldltSolve(system, size, space.swaps, space.blocks, solution);
      for (std::uint32_t a = b; a < n; ++a) {
        covariance[packedIndex(a, b)] = solution[a];
      }
    }
  }
  return true;
}

} // namespace decay_fit_detail

} // namespace o2::gtrk

#endif // GTRK_TRACKING_DECAYFITCONSTRAINTS_H_
