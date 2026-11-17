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

/// \file DecayRowsFit.h
/// \brief One decay vertex fitted from its daughters' rows directly: the
///        information limit the fit from the daughters' states
///        (DecayVertexFit.h) is held to.
///
/// One least-squares fit of every daughter's rows with the vertex. The
/// parameters: the vertex V; each daughter's direction and q/p as it was made
/// at V; and a kick -- in phi, in theta and in q/p -- at V and at each of its
/// rows but the last, the material's there, each with its prior: zero mean,
/// the scattering's and the straggling's variance (a General Broken Lines
/// fit). The material model is the Kalman fit's: a step's material split in
/// halves, each half's scattering kicked at its end of the step -- the first
/// step's nearer half at V, so that the momentum parameters are the
/// daughters' as they were made, which the truth and the constraints mean --,
/// Highland's formula with its logarithm at the span of the rows
/// (`SpanMaterial`), the step from the vertex to the first row at its own; the
/// mean ionisation over every step (`energyLoss`), its straggling the q/p
/// kick's variance. Between planes the helix is exact in the uniform field, in
/// double, its derivatives by forward differentiation of the six-component
/// free state (position, phi, theta, q/p); the material and the loss are
/// re-evaluated every iteration, the loss's derivative by q/p taken (by a
/// central difference: left out, the fit's fixed point is not the minimum and
/// slow daughters come out too slow), the material's left out.
///
/// Gauss-Newton from (a)'s solution. The kicks start on a trajectory that
/// follows the rows (`followRows`): at every row the kick that turns the
/// direction onto the helix, at the current curvature, through the next
/// row's measured point; with kicks at zero a slow daughter's helix from the
/// vertex may miss its later planes altogether. A trial step is
/// evaluated with the material along its own trajectory, and halved while
/// its trajectory misses a plane or its chi2 rises against the current
/// point's under the same material: a slow daughter's reach depends on its
/// loss. Converged when the vertex and the momenta stop moving
/// (in their sigmas); the kicks are nuisance parameters.

#ifndef GTRK_STUDY_DECAYROWSFIT_H_
#define GTRK_STUDY_DECAYROWSFIT_H_

#include "GenericTracking/Study/DecayConstraints.h"
#include "GenericTracking/Study/DecayVertexFit.h"
#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/TrackFitterUtils.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <span>
#include <vector>

namespace o2::gtrk::study
{

/// A daughter as the rows' fit takes it: its rows' planes and measurements,
/// innermost first.
struct RowsDaughter {
  std::vector<TrackingSurface> planes;
  std::vector<std::array<double, 2>> positions;  ///< loc0, loc1 on its plane [cm]
  std::vector<std::array<double, 3>> covariance; ///< var loc0, cov, var loc1 [cm^2]
  Propagator propagator;                         ///< the field, with the daughter's species' mass
  /// An electron: its radiation on every step, the Bethe-Heitler mean loss
  /// and, as the q/p kick's prior, its variance (a Gaussian in place of the
  /// loss's density, which the Gaussian-sum fit keeps).
  bool electron = false;
};

namespace decay_rows_detail
{

using Free = std::array<decay_fit_detail::Dual<6>, 6>; ///< x, y, z, phi, theta, q/p, by the step's start

/// \a from carried along its helix to \a plane (a plane along the beam, or
/// across it): the free state at the crossing nearest along the helix,
/// ahead or behind, and its loc0, loc1 there; false if the helix never
/// reaches the plane or the plane is tilted.
///
/// On a plane along the beam with normal angle phiL the crossing has
/// sin(phi_end - phiL) = sin(phi - phiL) - rho u, u the start's distance
/// along the normal; of the two solutions phi_end and every turn by 2 pi, the
/// smallest |rho s| = |phi_end - phi| is the nearest crossing. In double the
/// difference loses nothing that matters at small rho.
[[nodiscard]] inline bool crossToPlane(const Free& from, const TrackingSurface& plane, double fieldScale, Free& to, std::array<decay_fit_detail::Dual<6>, 2>& local)
{
  using namespace decay_fit_detail;
  using D = Dual<6>;
  const double nx = plane.normal[0];
  const double ny = plane.normal[1];
  const double nz = plane.normal[2];
  const D& phi = from[3];
  const D& theta = from[4];
  const D sinTheta = sin(theta);
  const D cotTheta = cos(theta) / sinTheta;
  const D rho = fieldScale * (from[5] / sinTheta);
  D s;
  if (std::abs(nz) < 1e-6) {
    const double phiL = std::atan2(ny, nx);
    const D u = (nx * (from[0] + (-static_cast<double>(plane.origin[0])))) + (ny * (from[1] + (-static_cast<double>(plane.origin[1]))));
    const D delta = phi + (-phiL);
    if (rho.v == 0.) {
      const D cosDelta = cos(delta);
      if (!(std::abs(cosDelta.v) > 0.)) {
        return false;
      }
      s = (-u) / cosDelta;
    } else {
      const D sinEnd = sin(delta) - (rho * u);
      if (!(std::abs(sinEnd.v) <= 1.)) {
        return false;
      }
      const D first = chain(sinEnd, std::asin(sinEnd.v), 1. / std::sqrt(std::max(1e-300, 1. - (sinEnd.v * sinEnd.v))));
      // The turn to each solution, brought within (-pi, pi]; the smaller wins.
      const auto turnTo = [&](const D& end) {
        D turn = end - delta;
        turn.v = std::remainder(turn.v, 2. * std::numbers::pi);
        return turn;
      };
      const D one = turnTo(first);
      const D two = turnTo((-first) + std::numbers::pi);
      s = (std::abs(one.v) <= std::abs(two.v) ? one : two) / rho;
    }
  } else if (std::abs(nx) < 1e-6 && std::abs(ny) < 1e-6) {
    if (!(std::abs(cotTheta.v) > 0.)) {
      return false;
    }
    s = (((-1.) * from[2]) + static_cast<double>(plane.origin[2])) / cotTheta;
  } else {
    return false;
  }
  const D half = 0.5 * (rho * s);
  const D chord = s * sinc(half);
  to[0] = from[0] + (chord * cos(phi + half));
  to[1] = from[1] + (chord * sin(phi + half));
  to[2] = from[2] + (s * cotTheta);
  to[3] = phi + (rho * s);
  to[4] = theta;
  to[5] = from[5];
  const std::array<D, 3> d{to[0] + (-static_cast<double>(plane.origin[0])), to[1] + (-static_cast<double>(plane.origin[1])), to[2] + (-static_cast<double>(plane.origin[2]))};
  local[0] = (static_cast<double>(plane.axis0[0]) * d[0]) + (static_cast<double>(plane.axis0[1]) * d[1]) + (static_cast<double>(plane.axis0[2]) * d[2]);
  local[1] = (static_cast<double>(plane.axis1[0]) * d[0]) + (static_cast<double>(plane.axis1[1]) * d[1]) + (static_cast<double>(plane.axis1[2]) * d[2]);
  return true;
}

} // namespace decay_rows_detail

/// The vertex of \a daughters fitted from their rows, starting from \a start
/// (the states' fit's solution), with \a material's material.
template <typename Material>
[[nodiscard]] DecayVertexFitResult fitDecayRows(std::span<const RowsDaughter> daughters, const Material& material, const DecayVertexFitResult& start, const DecayConstraints& constraints = {}, int maxIterations = 20)
{
  using namespace decay_fit_detail;
  using namespace decay_rows_detail;
  using D = Dual<6>;
  DecayVertexFitResult result;
  result.status = DecayFitStatus::kNotConverged;
  const int count = static_cast<int>(daughters.size());
  if (count < 2 || start.status != DecayFitStatus::kOk) {
    result.status = DecayFitStatus::kNoStart;
    return result;
  }

  // The parameter layout: V, then per daughter its phi, theta, q/p as it was
  // made, its kick at V and its kicks at its rows but the last (three each).
  std::vector<int> firstOf(static_cast<std::size_t>(count));
  int n = 3;
  for (int i = 0; i < count; ++i) {
    firstOf[i] = n;
    n += 6 + (3 * (static_cast<int>(daughters[i].planes.size()) - 1));
  }
  std::vector<double> parameters(static_cast<std::size_t>(n), 0.);
  std::copy(start.vertex.begin(), start.vertex.end(), parameters.begin());
  for (int i = 0; i < count; ++i) {
    for (int k = 0; k < 3; ++k) {
      parameters[firstOf[i] + k] = start.momenta[i][k];
    }
  }

  // Per daughter and row the material model at the current trajectory:
  // the kick's variances after the row, the loss's factor on q/p before it.
  struct RowMaterial {
    std::array<double, 3> kickVariance{}; ///< phi, theta, q/p
    double lossScale = 1.;
    double lossSlope = 0.; ///< d lossScale / d(q/p): a slow daughter's loss depends on its momentum
  };
  std::vector<std::vector<RowMaterial>> model(static_cast<std::size_t>(count));
  std::vector<std::array<double, 3>> vertexKick(static_cast<std::size_t>(count)); ///< the kick at V's variances

  // The chi2 at \a p, and with \a normal its normal equations; with \a refresh
  // the material model is evaluated along the trajectory first.
  const auto evaluate = [&](const std::vector<double>& p, bool refresh, double& chi2, std::vector<double>* normal, std::vector<double>* gradient) {
    chi2 = 0.;
    if (normal != nullptr) {
      std::ranges::fill(*normal, 0.);
      std::ranges::fill(*gradient, 0.);
    }
    for (int i = 0; i < count; ++i) {
      const auto& daughter = daughters[i];
      const int rows = static_cast<int>(daughter.planes.size());
      const double fieldScale = static_cast<double>(daughter.propagator.bzKiloGauss) * static_cast<double>(constants::math::B2C);
      const double mass = daughter.propagator.massGeV;
      const int base = firstOf[i];
      const int local = 9 + (3 * (rows - 1)); // V, the momentum, the kicks: this daughter's parameters
      // index of local parameter -> global
      std::vector<int> index(static_cast<std::size_t>(local));
      index[0] = 0;
      index[1] = 1;
      index[2] = 2;
      for (int k = 0; k < local - 3; ++k) {
        index[3 + k] = base + k;
      }
      auto& rowModel = model[i];
      if (refresh) {
        rowModel.assign(static_cast<std::size_t>(rows), {});
      }
      // The free state leaving V, the momentum kicked there, and its
      // derivatives by the local parameters (6 x local).
      std::array<double, 6> state{p[0], p[1], p[2], p[base] + p[base + 3], p[base + 1] + p[base + 4], p[base + 2] + p[base + 5]};
      std::vector<double> derivative(static_cast<std::size_t>(6 * local), 0.);
      for (int k = 0; k < 6; ++k) {
        derivative[k * local + k] = 1.;
      }
      for (int k = 3; k < 6; ++k) {
        derivative[k * local + k + 3] = 1.;
      }
      double spanX2X0 = 0.;
      if (refresh) {
        // The span's material along the rows, for Highland's logarithm.
        for (int row = 1; row < rows; ++row) {
          Scalar a[3]{};
          Scalar b[3]{};
          surfaceToGlobal(daughter.planes[row - 1], static_cast<Scalar>(daughter.positions[row - 1][0]), static_cast<Scalar>(daughter.positions[row - 1][1]), a);
          surfaceToGlobal(daughter.planes[row], static_cast<Scalar>(daughter.positions[row][0]), static_cast<Scalar>(daughter.positions[row][1]), b);
          spanX2X0 += material.getMatBudget(a[0], a[1], a[2], b[0], b[1], b[2]).meanX2X0;
        }
      }
      for (int row = 0; row < rows; ++row) {
        Free from{};
        for (int k = 0; k < 6; ++k) {
          from[k].v = state[k];
          from[k].d[k] = 1.;
        }
        Free to{};
        std::array<D, 2> loc{};
        if (!crossToPlane(from, daughter.planes[row], fieldScale, to, loc)) {
          return false;
        }
        if (refresh) {
          // This step's material, along the chord from where it starts.
          const Scalar a[3] = {static_cast<Scalar>(state[0]), static_cast<Scalar>(state[1]), static_cast<Scalar>(state[2])};
          const Scalar b[3] = {static_cast<Scalar>(to[0].v), static_cast<Scalar>(to[1].v), static_cast<Scalar>(to[2].v)};
          auto step = o2::gtrk::detail::stepMaterial(material, a, b);
          if (row != 0) {
            step.nearerLogX2X0 = static_cast<Scalar>(spanX2X0);
            step.fartherLogX2X0 = static_cast<Scalar>(spanX2X0);
          }
          const double momentum = 1. / std::abs(state[5]);
          const double farther = highlandTheta0(step.fartherX2X0, step.fartherLogX2X0, static_cast<Scalar>(momentum), static_cast<Scalar>(mass));
          const auto loss = energyLoss(static_cast<Scalar>(state[5]), static_cast<Scalar>(mass), -step.xRho);
          // The step's nearer half kicks at its start -- the previous row, or
          // the vertex --, its farther half at this row, with the loss's
          // straggling.
          const double nearer = highlandTheta0(step.nearerX2X0, step.nearerLogX2X0, static_cast<Scalar>(momentum), static_cast<Scalar>(mass));
          if (row != 0) {
            rowModel[row - 1].kickVariance[0] += nearer * nearer;
            rowModel[row - 1].kickVariance[1] += nearer * nearer;
          } else {
            vertexKick[i] = {nearer * nearer, nearer * nearer, 0.};
          }
          auto& here = rowModel[row];
          here.lossScale = loss.ok ? loss.momentumScale : 1.;
          double radiationVariance = 0.;
          if (daughter.electron) {
            // z = E after / E before: E[z] = 2^-c, E[z^2] = 3^-c, c = t / ln 2;
            // q/p scales by 1 / E[z], its variance (q/p)^2 var(z) / E[z]^4.
            const double c = static_cast<double>(step.nearerX2X0 + step.fartherX2X0) / std::numbers::ln2;
            const double meanZ = std::pow(2., -c);
            const double varZ = std::pow(3., -c) - (meanZ * meanZ);
            const double kAfter = state[5] * here.lossScale;
            here.lossScale /= meanZ;
            radiationVariance = kAfter * kAfter * varZ / (meanZ * meanZ * meanZ * meanZ);
          }
          // Its slope by q/p, by a central difference of 1 %.
          const double h = 0.01 * state[5];
          const auto up = energyLoss(static_cast<Scalar>(state[5] + h), static_cast<Scalar>(mass), -step.xRho);
          const auto down = energyLoss(static_cast<Scalar>(state[5] - h), static_cast<Scalar>(mass), -step.xRho);
          here.lossSlope = loss.ok && up.ok && down.ok ? (static_cast<double>(up.momentumScale) - static_cast<double>(down.momentumScale)) / (2. * h) : 0.;
          here.kickVariance = {farther * farther, farther * farther, (loss.ok ? loss.varQOverP : 0.) + radiationVariance};
        }
        // Chain the derivatives: d(to)/d(local) = d(to)/d(from) d(from)/d(local).
        std::vector<double> next(static_cast<std::size_t>(6 * local), 0.);
        for (int r = 0; r < 6; ++r) {
          for (int c = 0; c < local; ++c) {
            double value = 0.;
            for (int k = 0; k < 6; ++k) {
              value += to[r].d[k] * derivative[k * local + c];
            }
            next[r * local + c] = value;
          }
        }
        // The residual on this row and its derivatives.
        std::array<std::vector<double>, 2> h{std::vector<double>(static_cast<std::size_t>(local), 0.), std::vector<double>(static_cast<std::size_t>(local), 0.)};
        for (int a = 0; a < 2; ++a) {
          for (int c = 0; c < local; ++c) {
            double value = 0.;
            for (int k = 0; k < 6; ++k) {
              value += loc[a].d[k] * derivative[k * local + c];
            }
            h[a][c] = value;
          }
        }
        const std::array<double, 2> residual{daughter.positions[row][0] - loc[0].v, daughter.positions[row][1] - loc[1].v};
        const auto& cov = daughter.covariance[row];
        const double det = (cov[0] * cov[2]) - (cov[1] * cov[1]);
        if (!(det > 0.)) {
          return false;
        }
        const std::array<double, 4> weight{cov[2] / det, -cov[1] / det, -cov[1] / det, cov[0] / det};
        const std::array<double, 2> weighted{(weight[0] * residual[0]) + (weight[1] * residual[1]), (weight[2] * residual[0]) + (weight[3] * residual[1])};
        chi2 += (residual[0] * weighted[0]) + (residual[1] * weighted[1]);
        if (normal != nullptr) {
          for (int a = 0; a < local; ++a) {
            const double g = (h[0][a] * weighted[0]) + (h[1][a] * weighted[1]);
            (*gradient)[index[a]] += g;
            for (int b = 0; b < local; ++b) {
              const double value = (h[0][a] * ((weight[0] * h[0][b]) + (weight[1] * h[1][b]))) + (h[1][a] * ((weight[2] * h[0][b]) + (weight[3] * h[1][b])));
              (*normal)[index[a] * n + index[b]] += value;
            }
          }
        }
        // The state after the row: its loss, then its kick (none after the last).
        for (int k = 0; k < 6; ++k) {
          state[k] = to[k].v;
        }
        derivative.swap(next);
        // q/p after the loss: k scale(k), its derivative scale + k scale'.
        const double lossDerivative = rowModel[row].lossScale + (state[5] * rowModel[row].lossSlope);
        state[5] *= rowModel[row].lossScale;
        for (int c = 0; c < local; ++c) {
          derivative[5 * local + c] *= lossDerivative;
        }
        if (row + 1 < rows) {
          const int kick = 9 + (3 * row); // local index of this row's kick
          for (int k = 0; k < 3; ++k) {
            state[3 + k] += p[index[kick + k]];
            derivative[(3 + k) * local + kick + k] += 1.;
          }
        }
      }
    }
    // The kicks' priors (with the model of this evaluation): at V, then at
    // the rows.
    for (int i = 0; i < count; ++i) {
      const int rows = static_cast<int>(daughters[i].planes.size());
      for (int row = -1; row + 1 < rows; ++row) {
        for (int k = 0; k < 3; ++k) {
          const int at = firstOf[i] + 6 + (3 * row) + k;
          double variance = row < 0 ? vertexKick[i][k] : model[i][row].kickVariance[k];
          if (k == 0) {
            const double s = std::sin(p[firstOf[i] + 1]);
            variance /= s * s; // phi's kick: theta0 / sin(theta)
          }
          // A kick without material is held at zero.
          const double w = variance > 0. ? 1. / variance : 1e30;
          chi2 += p[at] * p[at] * w;
          if (normal != nullptr) {
            (*normal)[at * n + at] += w;
            (*gradient)[at] -= p[at] * w;
          }
        }
      }
    }
    if (constraints.prior) {
      chi2 += addLinePrior(*constraints.prior, std::span<const double>{p.data(), 3}, normal, gradient, n);
    }
    return true;
  };

  // The kicks' start: each daughter walked along its rows, its direction at
  // every row turned towards the next row's measured point.
  const auto followRows = [&](std::vector<double>& p) {
    for (int i = 0; i < count; ++i) {
      const auto& daughter = daughters[i];
      const int rows = static_cast<int>(daughter.planes.size());
      const double fieldScale = static_cast<double>(daughter.propagator.bzKiloGauss) * static_cast<double>(constants::math::B2C);
      const double mass = daughter.propagator.massGeV;
      const int base = firstOf[i];
      std::array<double, 6> state{p[0], p[1], p[2], p[base] + p[base + 3], p[base + 1] + p[base + 4], p[base + 2] + p[base + 5]};
      for (int row = 0; row + 1 < rows; ++row) {
        Free from{};
        for (int k = 0; k < 6; ++k) {
          from[k].v = state[k];
        }
        Free to{};
        std::array<D, 2> loc{};
        if (!crossToPlane(from, daughter.planes[row], fieldScale, to, loc)) {
          return; // the start stays where it got to
        }
        const Scalar a[3] = {static_cast<Scalar>(state[0]), static_cast<Scalar>(state[1]), static_cast<Scalar>(state[2])};
        const Scalar b[3] = {static_cast<Scalar>(to[0].v), static_cast<Scalar>(to[1].v), static_cast<Scalar>(to[2].v)};
        const auto step = o2::gtrk::detail::stepMaterial(material, a, b);
        const auto loss = energyLoss(static_cast<Scalar>(state[5]), static_cast<Scalar>(mass), -step.xRho);
        for (int k = 0; k < 6; ++k) {
          state[k] = to[k].v;
        }
        state[5] *= loss.ok ? loss.momentumScale : 1.;
        // The next row's measured point, and the helix towards it.
        Scalar next[3]{};
        surfaceToGlobal(daughter.planes[row + 1], static_cast<Scalar>(daughter.positions[row + 1][0]), static_cast<Scalar>(daughter.positions[row + 1][1]), next);
        const double cx = next[0] - state[0];
        const double cy = next[1] - state[1];
        const double chord = std::hypot(cx, cy);
        const double rho = fieldScale * state[5] / std::sin(state[4]);
        const double halfTurn = rho * chord / 2.;
        if (!(chord > 0.) || !(std::abs(halfTurn) < 1.)) {
          continue; // no helix of this curvature reaches it: no kick
        }
        const double turn = 2. * std::asin(halfTurn); // rho s
        const double arc = std::abs(rho) > 0. ? turn / rho : chord;
        const double phi = std::atan2(cy, cx) - (turn / 2.);
        const double theta = std::atan2(arc, next[2] - state[2]);
        const int kick = base + 6 + (3 * row);
        p[kick] = std::remainder(phi - state[3], 2. * std::numbers::pi);
        p[kick + 1] = theta - state[4];
        state[3] += p[kick];
        state[4] += p[kick + 1];
      }
    }
  };
  followRows(parameters);

  // The constraints' values and rows at \a p, in this fit's layout (their
  // columns are V and the momenta; the kicks' are zero).
  const bool constrained = constraints.count() > 0;
  std::vector<double> masses(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    masses[i] = daughters[i].propagator.massGeV;
  }
  const double fieldScale = static_cast<double>(daughters[0].propagator.bzKiloGauss) * static_cast<double>(constants::math::B2C);
  std::vector<double> values;
  std::vector<double> rows;
  const auto constraintsAt = [&](const std::vector<double>& p) {
    std::vector<std::array<double, 3>> momenta(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
      momenta[i] = {p[firstOf[i]], p[firstOf[i] + 1], p[firstOf[i] + 2]};
    }
    std::vector<double> core;
    if (!evaluateConstraints(constraints, {p[0], p[1], p[2]}, momenta, masses, fieldScale, values, core)) {
      return false;
    }
    const int width = 3 + (3 * count);
    const int m = static_cast<int>(values.size());
    rows.assign(static_cast<std::size_t>(m * n), 0.);
    for (int a = 0; a < m; ++a) {
      for (int k = 0; k < 3; ++k) {
        rows[a * n + k] = core[a * width + k];
      }
      for (int i = 0; i < count; ++i) {
        for (int k = 0; k < 3; ++k) {
          rows[a * n + firstOf[i] + k] = core[a * width + 3 + (3 * i) + k];
        }
      }
    }
    return true;
  };

  std::vector<double> normal(static_cast<std::size_t>(n * n));
  std::vector<double> gradient(static_cast<std::size_t>(n));
  std::vector<double> factor;
  double chi2 = 0.;
  bool converged = false;
  for (int iteration = 0; iteration < maxIterations; ++iteration) {
    result.iterations = static_cast<std::uint8_t>(iteration + 1);
    if (!evaluate(parameters, true, chi2, &normal, &gradient)) {
      result.status = DecayFitStatus::kNoCrossing;
      return result;
    }
    factor = normal;
    if (!cholesky(factor, n)) {
      result.status = DecayFitStatus::kSingular;
      return result;
    }
    std::vector<double> step = gradient;
    if (constrained) {
      if (!constraintsAt(parameters)) {
        result.status = DecayFitStatus::kNoCrossing;
        return result;
      }
      if (!constrainedStep(factor, n, gradient, rows, values, step, nullptr)) {
        result.status = DecayFitStatus::kSingular;
        return result;
      }
    } else {
      choleskySolve(factor, n, step);
    }
    std::vector<double> trial(parameters.size());
    double trialChi2 = 0.;
    double scale = 1.;
    bool taken = false;
    const auto kept = model;
    for (int halving = 0; halving < 10; ++halving, scale *= 0.5) {
      for (int k = 0; k < n; ++k) {
        trial[k] = parameters[k] + (scale * step[k]);
      }
      // The trial and the current point, both under the trial trajectory's
      // material: the priors move with the model.
      double currentChi2 = 0.;
      if (evaluate(trial, true, trialChi2, nullptr, nullptr) && (constrained || (evaluate(parameters, false, currentChi2, nullptr, nullptr) && trialChi2 <= currentChi2 * (1. + 1e-12) + 1e-12))) {
        taken = true;
        break;
      }
    }
    if (!taken) {
      model = kept;
      converged = !constrained;
      break;
    }
    const auto inverse = choleskyInverse(factor, n);
    double largest = 0.;
    const auto watch = [&](int k) { largest = std::max(largest, std::abs(scale * step[k]) / std::sqrt(inverse[k * n + k])); };
    for (int k = 0; k < 3; ++k) {
      watch(k);
    }
    for (int i = 0; i < count; ++i) {
      for (int k = 0; k < 3; ++k) {
        watch(firstOf[i] + k);
      }
    }
    parameters = trial;
    chi2 = trialChi2;
    if (largest < 1e-3) {
      converged = true;
      break;
    }
  }
  result.status = converged ? DecayFitStatus::kOk : DecayFitStatus::kNotConverged;
  if (!evaluate(parameters, true, chi2, &normal, &gradient)) {
    result.status = DecayFitStatus::kNoCrossing;
    return result;
  }
  factor = normal;
  if (!cholesky(factor, n)) {
    result.status = DecayFitStatus::kSingular;
    return result;
  }
  // The covariance of V and the momenta: the full inverse's rows of them.
  std::vector<double> inverse;
  if (constrained) {
    std::vector<double> unused;
    if (!constraintsAt(parameters) || !constrainedStep(factor, n, gradient, rows, values, unused, &inverse)) {
      result.status = DecayFitStatus::kSingular;
      return result;
    }
  } else {
    inverse = choleskyInverse(factor, n);
  }
  const int m = 3 + (3 * count);
  std::vector<int> kept{0, 1, 2};
  for (int i = 0; i < count; ++i) {
    for (int k = 0; k < 3; ++k) {
      kept.push_back(firstOf[i] + k);
    }
  }
  result.covariance.assign(static_cast<std::size_t>(m * m), 0.);
  for (int a = 0; a < m; ++a) {
    for (int b = 0; b < m; ++b) {
      result.covariance[a * m + b] = inverse[kept[a] * n + kept[b]];
    }
  }
  std::copy_n(parameters.begin(), 3, result.vertex.begin());
  result.momenta.resize(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    for (int k = 0; k < 3; ++k) {
      result.momenta[i][k] = parameters[firstOf[i] + k];
    }
  }
  result.chi2 = chi2;
  int measurements = 0;
  for (const auto& daughter : daughters) {
    measurements += 2 * static_cast<int>(daughter.planes.size());
  }
  result.ndf = measurements - (3 + (3 * count)) + constraints.count() + (constraints.prior ? 1 : 0);
  result.carried = start.carried;
  result.carriedCovariance = start.carriedCovariance;
  result.carryX2X0 = start.carryX2X0;
  return result;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYROWSFIT_H_
