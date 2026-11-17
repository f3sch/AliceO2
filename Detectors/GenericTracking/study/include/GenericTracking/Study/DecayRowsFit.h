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
/// parameters: the vertex V; each daughter's direction and q/p at V; and at
/// each of its rows but the last a kick -- in phi, in theta and in q/p --
/// that its material gives it there, each with its prior: zero mean, the
/// scattering's and the straggling's variance (a General Broken Lines fit).
/// The material model is the Kalman fit's: a step's material split in halves,
/// each half's scattering kicked at its end of the step, Highland's formula
/// with its logarithm at the span of the rows (`SpanMaterial`), the step
/// from the vertex to the first row at its own; the mean ionisation applied
/// over every step (`energyLoss`); the vertex's half of the first step is
/// absorbed by the free direction at V. Between planes the helix is exact in
/// the uniform field, in double, its derivatives by forward differentiation
/// of the six-component free state (position, phi, theta, q/p); the
/// derivatives of the material and the loss by the state are left out (they
/// are re-evaluated every iteration).
///
/// Gauss-Newton, the start (a)'s solution, kicks at zero.

#ifndef GTRK_STUDY_DECAYROWSFIT_H_
#define GTRK_STUDY_DECAYROWSFIT_H_

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
};

namespace decay_rows_detail
{

using Free = std::array<decay_fit_detail::Dual<6>, 6>; ///< x, y, z, phi, theta, q/p, by the step's start

/// \a from carried along its helix to \a plane (a plane along the beam, or
/// across it): the free state there and its loc0, loc1; false if the helix
/// does not reach it ahead or the plane is tilted.
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
    double phiL = std::atan2(ny, nx);
    D u = (nx * (from[0] + (-static_cast<double>(plane.origin[0])))) + (ny * (from[1] + (-static_cast<double>(plane.origin[1]))));
    if (std::cos(phi.v - phiL) < 0.) {
      phiL += std::numbers::pi;
      u = -u;
    }
    const D delta = phi + (-phiL);
    const D sinDelta = sin(delta);
    const D cosDelta = cos(delta);
    const D sinEnd = sinDelta - (rho * u);
    if (!(std::abs(sinEnd.v) < 1.)) {
      return false;
    }
    const D rootEnd = sqrt((-(sinEnd * sinEnd)) + 1.);
    const D wOverRho = (-(u * cosDelta)) - ((sinDelta * u) * (sinEnd + sinDelta) / (cosDelta + rootEnd));
    const D w = rho * wOverRho;
    if (!(std::abs(w.v) < 1.)) {
      return false;
    }
    s = wOverRho * asinc(w);
  } else if (std::abs(nx) < 1e-6 && std::abs(ny) < 1e-6) {
    if (!(std::abs(cotTheta.v) > 1e-9)) {
      return false;
    }
    s = (((-1.) * from[2]) + static_cast<double>(plane.origin[2])) / cotTheta;
  } else {
    return false;
  }
  if (!(s.v > 0.)) {
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
[[nodiscard]] DecayVertexFitResult fitDecayRows(std::span<const RowsDaughter> daughters, const Material& material, const DecayVertexFitResult& start, int maxIterations = 20)
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

  // The parameter layout: V, then per daughter its phi, theta, q/p and its
  // kicks (three per row but the last).
  std::vector<int> firstOf(static_cast<std::size_t>(count));
  int n = 3;
  for (int i = 0; i < count; ++i) {
    firstOf[i] = n;
    n += 3 + (3 * (static_cast<int>(daughters[i].planes.size()) - 1));
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
  };
  std::vector<std::vector<RowMaterial>> model(static_cast<std::size_t>(count));

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
      const int local = 6 + (3 * (rows - 1)); // V, the momentum, the kicks: this daughter's parameters
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
      // The free state and its derivatives by the local parameters (6 x local).
      std::array<double, 6> state{p[0], p[1], p[2], p[base], p[base + 1], p[base + 2]};
      std::vector<double> derivative(static_cast<std::size_t>(6 * local), 0.);
      for (int k = 0; k < 6; ++k) {
        derivative[k * local + k] = 1.;
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
          // The step's nearer half kicks at its start -- the previous row; at
          // the vertex the free direction absorbs it --, its farther half at
          // this row, with the loss's straggling.
          if (row != 0) {
            const double nearer = highlandTheta0(step.nearerX2X0, step.nearerLogX2X0, static_cast<Scalar>(momentum), static_cast<Scalar>(mass));
            rowModel[row - 1].kickVariance[0] += nearer * nearer;
            rowModel[row - 1].kickVariance[1] += nearer * nearer;
          }
          auto& here = rowModel[row];
          here.lossScale = loss.ok ? loss.momentumScale : 1.;
          here.kickVariance = {farther * farther, farther * farther, loss.ok ? loss.varQOverP : 0.};
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
        state[5] *= rowModel[row].lossScale;
        for (int c = 0; c < local; ++c) {
          derivative[5 * local + c] *= rowModel[row].lossScale;
        }
        if (row + 1 < rows) {
          const int kick = 6 + (3 * row); // local index of this row's kick
          for (int k = 0; k < 3; ++k) {
            state[3 + k] += p[index[kick + k]];
            derivative[(3 + k) * local + kick + k] += 1.;
          }
        }
      }
    }
    // The kicks' priors (with the model of this evaluation).
    for (int i = 0; i < count; ++i) {
      const int rows = static_cast<int>(daughters[i].planes.size());
      for (int row = 0; row + 1 < rows; ++row) {
        for (int k = 0; k < 3; ++k) {
          const int at = firstOf[i] + 3 + (3 * row) + k;
          double variance = model[i][row].kickVariance[k];
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
    choleskySolve(factor, n, step);
    std::vector<double> trial(parameters.size());
    double trialChi2 = 0.;
    double scale = 1.;
    bool taken = false;
    for (int halving = 0; halving < 10; ++halving, scale *= 0.5) {
      for (int k = 0; k < n; ++k) {
        trial[k] = parameters[k] + (scale * step[k]);
      }
      if (evaluate(trial, false, trialChi2, nullptr, nullptr) && trialChi2 <= chi2 * (1. + 1e-12) + 1e-12) {
        taken = true;
        break;
      }
    }
    if (!taken) {
      converged = true;
      break;
    }
    const auto inverse = choleskyInverse(factor, n);
    double largest = 0.;
    for (int k = 0; k < n; ++k) {
      largest = std::max(largest, std::abs(scale * step[k]) / std::sqrt(inverse[k * n + k]));
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
  const auto inverse = choleskyInverse(factor, n);
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
  result.ndf = measurements - (3 + (3 * count));
  result.carried = start.carried;
  result.carriedCovariance = start.carriedCovariance;
  result.carryX2X0 = start.carryX2X0;
  return result;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYROWSFIT_H_
