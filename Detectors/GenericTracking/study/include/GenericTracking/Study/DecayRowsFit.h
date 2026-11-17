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
/// re-evaluated every iteration, their derivatives by the state left out. A
/// kick whose material is nil is no parameter: it is held at zero.
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
#include <limits>
#include <functional>
#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <ranges>
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

/// A charged composite's own rows (a tracked Xi or Omega), from its decay
/// vertex back -- its last row first --: walked from the vertex against its
/// momentum, its state there its daughters' summed momentum and charge, a
/// kick at the vertex and after every row but the last, the loss regained
/// on the way back. The root's walk, under a production point, ends there:
/// a kick after its innermost row too, then the plane about the point
/// (through it along the beam, normal to the walk's direction at its closest
/// approach), measured at the point with `pointCovariance`, the scattering
/// between the innermost row and the point; the point is then not predicted
/// from the vertex.
struct RowsMother {
  RowsDaughter rows;                       ///< its rows, its last row first; `propagator` with its mass
  int vertex = 0;                          ///< the vertex it decayed at
  std::array<double, 3> pointCovariance{}; ///< (loc0, loc0 loc1, loc1) [cm^2] on the plane about the production point
};

namespace decay_rows_detail
{

using Free = std::array<decay_fit_detail::Dual<6>, 6>; ///< x, y, z, phi, theta, q/p, by the step's start

/// \a from carried along its helix to \a plane (a plane along the beam, or
/// across it): the free state at the crossing nearest along the helix,
/// ahead or behind, and its loc0, loc1 there; false if the helix never
/// reaches the plane.
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
  if (nz == 0.) {
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
      if (!(std::abs(sinEnd.v) < 1.)) {
        return false; // missed, or grazing: no crossing to bind to
      }
      const D first = chain(sinEnd, std::asin(sinEnd.v), 1. / std::sqrt(1. - (sinEnd.v * sinEnd.v)));
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
  } else if (nx == 0. && ny == 0.) {
    if (!(std::abs(cotTheta.v) > 0.)) {
      return false;
    }
    s = (((-1.) * from[2]) + static_cast<double>(plane.origin[2])) / cotTheta;
  } else {
    // A tilted plane: f(s) = n . (x(s) - origin) = 0 by Newton's method from
    // the straight line's crossing, while each step shrinks; s's derivatives
    // by the implicit function, -df / f'(s).
    const auto f = [&](const D& at) {
      const D half = 0.5 * (rho * at);
      const D chord = at * sinc(half);
      return (nx * ((from[0] + (chord * cos(phi + half))) + (-static_cast<double>(plane.origin[0])))) + (ny * ((from[1] + (chord * sin(phi + half))) + (-static_cast<double>(plane.origin[1])))) + (nz * ((from[2] + (at * cotTheta)) + (-static_cast<double>(plane.origin[2]))));
    };
    const auto slope = [&](double at) { return (nx * std::cos(phi.v + (rho.v * at))) + (ny * std::sin(phi.v + (rho.v * at))) + (nz * cotTheta.v); };
    const double ahead = slope(0.);
    if (!(std::abs(ahead) > 0.)) {
      return false;
    }
    double at = -f(D{0.}).v / ahead;
    double last = std::numeric_limits<double>::infinity();
    for (int step = 0; step < kMaxSteps; ++step) {
      const double fs = slope(at);
      if (!(std::abs(fs) > 0.)) {
        return false;
      }
      const double change = f(D{at}).v / fs;
      if (!(std::abs(change) < last)) {
        break; // no longer shrinking: at the root to rounding
      }
      last = std::abs(change);
      at -= change;
    }
    const D value = f(D{at});
    // Closed to the plane's own precision (its geometry is floats), or no
    // crossing: f's terms' size times float epsilon.
    {
      const double half = 0.5 * rho.v * at;
      const double chord = at * (half == 0. ? 1. : std::sin(half) / half);
      const std::array<double, 3> x{from[0].v + (chord * std::cos(phi.v + half)), from[1].v + (chord * std::sin(phi.v + half)), from[2].v + (at * cotTheta.v)};
      double size = 0.;
      for (int c = 0; c < 3; ++c) {
        size += std::abs(static_cast<double>(plane.normal[c])) * (std::abs(x[c]) + std::abs(static_cast<double>(plane.origin[c])));
      }
      if (!(std::abs(value.v) <= static_cast<double>(std::numeric_limits<Scalar>::epsilon()) * size)) {
        return false;
      }
    }
    s = D{at};
    const double fs = slope(at);
    for (int k = 0; k < 6; ++k) {
      s.d[k] = -value.d[k] / fs;
    }
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
[[nodiscard]] DecayVertexFitResult fitDecayRows(std::span<const RowsDaughter> daughters, const Material& material, const DecayVertexFitResult& start, const DecayConstraints& requested = {}, const DecayTopology& topology = {}, std::span<const RowsMother> mothers = {})
{
  using namespace decay_fit_detail;
  using namespace decay_rows_detail;
  using D = Dual<6>;
  // A tracked root's walk ends at the production point: the point's
  // equations leave the constraints.
  const auto endsAtPoint = [&](const RowsMother& mother) { return requested.production && mother.vertex == 0; };
  DecayConstraints constraints = requested;
  constraints.production = requested.production && std::ranges::none_of(mothers, endsAtPoint);
  DecayVertexFitResult result;
  result.status = DecayFitStatus::kNotConverged;
  const int count = static_cast<int>(daughters.size());
  if (count < 2 || start.status != DecayFitStatus::kOk) {
    result.status = DecayFitStatus::kNoStart;
    return result;
  }

  // The parameter layout: the vertices (a tree's, DecayTopology), then per
  // daughter its phi, theta, q/p as it was made at its vertex, its kick there
  // and its kicks at its rows but the last (three each).
  const int vertexCount = topology.vertices();
  if (static_cast<int>(start.vertices.size()) != vertexCount) {
    result.status = DecayFitStatus::kNoStart;
    return result;
  }
  std::vector<int> firstOf(static_cast<std::size_t>(count));
  int n = 3 * vertexCount;
  for (int i = 0; i < count; ++i) {
    firstOf[i] = n;
    n += 6 + (3 * (static_cast<int>(daughters[i].planes.size()) - 1));
  }
  // Each mother's kicks: at its vertex, and after every row but the last.
  std::vector<int> firstOfMother(mothers.size());
  for (std::size_t j = 0; j < mothers.size(); ++j) {
    firstOfMother[j] = n;
    n += 3 * (static_cast<int>(mothers[j].rows.planes.size()) + (endsAtPoint(mothers[j]) ? 1 : 0));
  }
  // A production point's two nuisances (`ScatteredPoint`): of the point's
  // equations, or of the tracked root's walk's end at it.
  const int etaAt = requested.production ? n : -1;
  n += requested.production ? 2 : 0;
  const auto ending = std::ranges::find_if(mothers, endsAtPoint);
  const bool walkEnds = ending != mothers.end();
  ScatteredPoint scattered{.at = etaAt, .factor = ScatteredPoint::factorOf(walkEnds ? ending->pointCovariance : requested.pointCovariance)};
  std::vector<double> parameters(static_cast<std::size_t>(n), 0.);
  for (int w = 0; w < vertexCount; ++w) {
    std::copy(start.vertices[w].begin(), start.vertices[w].end(), parameters.begin() + (3 * w));
  }
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
  std::vector<std::vector<RowMaterial>> motherModel(mothers.size());
  std::vector<std::array<double, 3>> motherVertexKick(mothers.size());
  std::vector<std::array<double, 3>> vertexKick(static_cast<std::size_t>(count)); ///< the kick at V's variances

  // The kicks without material: held at zero, out of the solve.
  std::vector<char> held(static_cast<std::size_t>(n), 0);
  const auto holdFixed = [&](std::vector<double>& p, std::vector<double>& matrix, std::vector<double>& rhs) {
    for (int k = 0; k < n; ++k) {
      if (held[k] != 0) {
        p[k] = 0.;
        for (int j = 0; j < n; ++j) {
          matrix[k * n + j] = 0.;
          matrix[j * n + k] = 0.;
        }
        matrix[k * n + k] = 1.;
        rhs[k] = 0.;
      }
    }
  };

  // The chi2 at \a p, and with \a normal its normal equations; with \a refresh
  // the material model is evaluated along the trajectory first.
  // The inputs' precision (`minimiseFrozen`), a walk's exact end's
  // equations, and whether a trajectory stopped in the material, of the last
  // evaluation.
  constexpr double kEpsilon = std::numeric_limits<Scalar>::epsilon(); // the measurements' precision
  double rounding = 0.;
  std::vector<double> terminalValues;
  std::vector<double> terminalRows;
  bool stopped = false;
  const auto evaluate = [&](const std::vector<double>& p, bool refresh, double& chi2, std::vector<double>* normal, std::vector<double>* gradient) {
    chi2 = 0.;
    rounding = 0.;
    terminalValues.clear();
    terminalRows.clear();
    if (normal != nullptr) {
      std::ranges::fill(*normal, 0.);
      std::ranges::fill(*gradient, 0.);
    }
    // One chain's walk: from \a state (its derivatives by the \a local
    // parameters, \a index into the global ones) through its rows, the kick
    // after row r at local index \a kickBase + 3 r; \a backward walks against
    // the momentum (a mother's), regaining the loss; with \a end (a mother's)
    // it ends at the production point, kicked after its last row too.
    const auto walk = [&](const RowsDaughter& chain, int local, const std::vector<int>& index, std::array<double, 6> state, std::vector<double> derivative, int kickBase, bool backward, std::vector<RowMaterial>& rowModel, std::array<double, 3>& firstKick, const RowsMother* end) {
      const int rows = static_cast<int>(chain.planes.size());
      const double fieldScale = static_cast<double>(chain.propagator.bzKiloGauss) * static_cast<double>(constants::math::B2C);
      const double mass = chain.propagator.massGeV;
      // A measured position on a plane: its chi2 and normal equations, the
      // crossing \a loc's derivatives chained through \a from's.
      const auto measure = [&](const std::array<D, 2>& loc, const std::vector<double>& from, const std::array<double, 2>& position, const std::array<double, 3>& cov) {
        std::array<std::vector<double>, 2> h{std::vector<double>(static_cast<std::size_t>(local), 0.), std::vector<double>(static_cast<std::size_t>(local), 0.)};
        for (int a = 0; a < 2; ++a) {
          for (int c = 0; c < local; ++c) {
            double value = 0.;
            for (int k = 0; k < 6; ++k) {
              value += loc[a].d[k] * from[k * local + c];
            }
            h[a][c] = value;
          }
        }
        const std::array<double, 2> residual{position[0] - loc[0].v, position[1] - loc[1].v};
        const double det = (cov[0] * cov[2]) - (cov[1] * cov[1]);
        if (!(cov[0] > 0.) || !(det > 0.)) {
          return false;
        }
        const std::array<double, 4> weight{cov[2] / det, -cov[1] / det, -cov[1] / det, cov[0] / det};
        const std::array<double, 2> weighted{(weight[0] * residual[0]) + (weight[1] * residual[1]), (weight[2] * residual[0]) + (weight[3] * residual[1])};
        chi2 += (residual[0] * weighted[0]) + (residual[1] * weighted[1]);
        rounding += (std::pow(kEpsilon * position[0], 2) * weight[0]) + (std::pow(kEpsilon * position[1], 2) * weight[3]);
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
        return true;
      };
      if (refresh) {
        rowModel.assign(static_cast<std::size_t>(rows), {});
      }
      double spanX2X0 = 0.;
      if (refresh) {
        // The span's material along the rows, for Highland's logarithm.
        for (int row = 1; row < rows; ++row) {
          Scalar a[3]{};
          Scalar b[3]{};
          surfaceToGlobal(chain.planes[row - 1], static_cast<Scalar>(chain.positions[row - 1][0]), static_cast<Scalar>(chain.positions[row - 1][1]), a);
          surfaceToGlobal(chain.planes[row], static_cast<Scalar>(chain.positions[row][0]), static_cast<Scalar>(chain.positions[row][1]), b);
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
        if (!crossToPlane(from, chain.planes[row], fieldScale, to, loc)) {
          return false;
        }
        if (refresh) {
          // This step's material, along the chord from where it starts.
          const Scalar a[3] = {static_cast<Scalar>(state[0]), static_cast<Scalar>(state[1]), static_cast<Scalar>(state[2])};
          const Scalar b[3] = {static_cast<Scalar>(to[0].v), static_cast<Scalar>(to[1].v), static_cast<Scalar>(to[2].v)};
          auto step = o2::gtrk::detail::stepMaterial(material, a, b);
          // Highland's logarithm at the span's x/X0, the first step's at its
          // own whole (its halves are one segment's).
          const Scalar logX2X0 = row != 0 ? static_cast<Scalar>(spanX2X0) : step.nearerX2X0 + step.fartherX2X0;
          step.nearerLogX2X0 = logX2X0;
          step.fartherLogX2X0 = logX2X0;
          const double momentum = 1. / std::abs(state[5]);
          const double farther = highlandTheta0(step.fartherX2X0, step.fartherLogX2X0, static_cast<Scalar>(momentum), static_cast<Scalar>(mass));
          const auto loss = energyLoss(static_cast<Scalar>(state[5]), static_cast<Scalar>(mass), backward ? step.xRho : -step.xRho);
          // The step's nearer half kicks at its start -- the previous row, or
          // the vertex --, its farther half at this row, with the loss's
          // straggling.
          const double nearer = highlandTheta0(step.nearerX2X0, step.nearerLogX2X0, static_cast<Scalar>(momentum), static_cast<Scalar>(mass));
          if (row != 0) {
            rowModel[row - 1].kickVariance[0] += nearer * nearer;
            rowModel[row - 1].kickVariance[1] += nearer * nearer;
          } else {
            firstKick = {nearer * nearer, nearer * nearer, 0.};
          }
          if (!loss.ok) {
            stopped = true; // the particle stops in this step's material
            return false;
          }
          auto& here = rowModel[row];
          here.lossScale = loss.momentumScale;
          double radiationVariance = 0.;
          if (chain.electron) {
            // z = E after / E before: E[z] = 2^-c, E[z^2] = 3^-c, c = t / ln 2;
            // q/p scales by 1 / E[z], its variance (q/p)^2 var(z) / E[z]^4.
            const double c = static_cast<double>(step.nearerX2X0 + step.fartherX2X0) / std::numbers::ln2;
            const double meanZ = std::pow(2., -c);
            const double varZ = std::pow(3., -c) - (meanZ * meanZ);
            const double kAfter = state[5] * here.lossScale;
            here.lossScale /= meanZ;
            radiationVariance = kAfter * kAfter * varZ / (meanZ * meanZ * meanZ * meanZ);
          }
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
        if (!measure(loc, derivative, chain.positions[row], chain.covariance[row])) {
          return false;
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
        if (row + 1 < rows || end != nullptr) {
          const int kick = kickBase + (3 * row); // local index of this row's kick
          for (int k = 0; k < 3; ++k) {
            state[3 + k] += p[index[kick + k]];
            derivative[(3 + k) * local + kick + k] += 1.;
          }
        }
      }
      if (end != nullptr) {
        // The plane about the production point, normal to the walk's
        // direction at its closest approach (`perigeePhi`), from the current
        // values; measured at the point.
        const auto& point = requested.point;
        const double curvature = fieldScale * state[5] / std::sin(state[4]);
        const double dx = state[0] - point[0];
        const double dy = state[1] - point[1];
        const double along = (dx * std::cos(state[3])) + (dy * std::sin(state[3]));
        const double across = (-dx * std::sin(state[3])) + (dy * std::cos(state[3]));
        const double phi = state[3] + std::atan2(-curvature * along, 1. + (curvature * across));
        const auto plane = beamLinePlane(static_cast<Scalar>(point[0]), static_cast<Scalar>(point[1]), static_cast<Scalar>(phi));
        Free from{};
        for (int k = 0; k < 6; ++k) {
          from[k].v = state[k];
          from[k].d[k] = 1.;
        }
        Free to{};
        std::array<D, 2> loc{};
        if (!crossToPlane(from, plane, fieldScale, to, loc)) {
          return false;
        }
        // The walk ends at the point: two equations, its scattering between
        // the innermost row and the point their nuisances' (`constraintsAt`).
        const std::array<double, 2> target{0., point[2]};
        for (int a = 0; a < 2; ++a) {
          terminalValues.push_back(loc[a].v - target[a]);
          std::vector<double> row(static_cast<std::size_t>(n), 0.);
          for (int c = 0; c < local; ++c) {
            double value = 0.;
            for (int k = 0; k < 6; ++k) {
              value += loc[a].d[k] * derivative[k * local + c];
            }
            row[index[c]] += value;
          }
          terminalRows.insert(terminalRows.end(), row.begin(), row.end());
        }
      }
      return true;
    };

    for (int i = 0; i < count; ++i) {
      const auto& daughter = daughters[i];
      const int rows = static_cast<int>(daughter.planes.size());
      const int base = firstOf[i];
      const int local = 9 + (3 * (rows - 1)); // V, the momentum, the kicks: this daughter's parameters
      // index of local parameter -> global
      std::vector<int> index(static_cast<std::size_t>(local));
      const int at = 3 * topology.vertexOfDaughter(i);
      index[0] = at;
      index[1] = at + 1;
      index[2] = at + 2;
      for (int k = 0; k < local - 3; ++k) {
        index[3 + k] = base + k;
      }
      // The free state leaving V, the momentum kicked there, and its
      // derivatives by the local parameters (6 x local).
      const std::array<double, 6> state{p[at], p[at + 1], p[at + 2], p[base] + p[base + 3], p[base + 1] + p[base + 4], p[base + 2] + p[base + 5]};
      std::vector<double> derivative(static_cast<std::size_t>(6 * local), 0.);
      for (int k = 0; k < 6; ++k) {
        derivative[k * local + k] = 1.;
      }
      for (int k = 3; k < 6; ++k) {
        derivative[k * local + k + 3] = 1.;
      }
      if (!walk(daughter, local, index, state, std::move(derivative), 9, false, model[i], vertexKick[i], nullptr)) {
        return false;
      }
    }

    // The mothers: from their vertex back through their rows, their state
    // there their daughters' summed momentum and charge, kicked there.
    for (std::size_t j = 0; j < mothers.size(); ++j) {
      const auto& mother = mothers[j];
      const int w = mother.vertex;
      std::vector<int> under;
      for (int i = 0; i < count; ++i) {
        if (topology.under(i, w)) {
          under.push_back(i);
        }
      }
      const int m = static_cast<int>(under.size());
      const int rows = static_cast<int>(mother.rows.planes.size());
      const int kicks = rows + (endsAtPoint(mother) ? 1 : 0);
      const int local = 3 + (3 * m) + (3 * kicks); // V, the daughters' momenta, the kicks
      std::vector<int> index(static_cast<std::size_t>(local));
      for (int k = 0; k < 3; ++k) {
        index[k] = (3 * w) + k;
      }
      for (int u = 0; u < m; ++u) {
        for (int k = 0; k < 3; ++k) {
          index[3 + (3 * u) + k] = firstOf[under[u]] + k;
        }
      }
      for (int k = 0; k < 3 * kicks; ++k) {
        index[3 + (3 * m) + k] = firstOfMother[j] + k;
      }
      using DM = Dual<kConstraintInputs>;
      std::array<DM, 3> total{};
      double charge = 0.;
      for (int u = 0; u < m; ++u) {
        std::array<DM, 3> q{};
        for (int k = 0; k < 3; ++k) {
          q[k].v = p[firstOf[under[u]] + k];
          q[k].d[(3 * u) + k] = 1.;
        }
        const double sign = q[2].v > 0. ? 1. : -1.;
        charge += sign;
        const DM pu = DM{sign} / q[2];
        const DM ptu = pu * sin(q[1]);
        total[0] = total[0] + (ptu * cos(q[0]));
        total[1] = total[1] + (ptu * sin(q[0]));
        total[2] = total[2] + (pu * cos(q[1]));
      }
      const DM pt = sqrt((total[0] * total[0]) + (total[1] * total[1]));
      const DM pTotal = sqrt((pt * pt) + (total[2] * total[2]));
      const std::array<DM, 3> direction{atan2(total[1], total[0]), atan2(pt, total[2]), DM{charge} / pTotal};
      const int kickV = 3 + (3 * m);
      std::array<double, 6> state{p[3 * w], p[(3 * w) + 1], p[(3 * w) + 2], 0., 0., 0.};
      std::vector<double> derivative(static_cast<std::size_t>(6 * local), 0.);
      for (int k = 0; k < 3; ++k) {
        derivative[k * local + k] = 1.;
        state[3 + k] = direction[k].v + p[index[kickV + k]];
        for (int c = 0; c < 3 * m; ++c) {
          derivative[(3 + k) * local + 3 + c] = direction[k].d[c];
        }
        derivative[(3 + k) * local + kickV + k] = 1.;
      }
      if (!walk(mother.rows, local, index, state, std::move(derivative), kickV + 3, true, motherModel[j], motherVertexKick[j], endsAtPoint(mother) ? &mother : nullptr)) {
        return false;
      }
      // Its kicks' priors: at the vertex, then after its rows (but the last,
      // unless it ends at the production point).
      for (int row = -1; row + 1 < kicks; ++row) {
        for (int k = 0; k < 3; ++k) {
          const int at = firstOfMother[j] + 3 + (3 * row) + k;
          double variance = row < 0 ? motherVertexKick[j][k] : motherModel[j][row].kickVariance[k];
          if (k == 0) {
            const double s = std::sin(direction[1].v);
            variance /= s * s;
          }
          if (!(variance > 0.)) {
            held[at] = 1;
            continue;
          }
          held[at] = 0;
          const double weightKick = 1. / variance;
          chi2 += p[at] * p[at] * weightKick;
          if (normal != nullptr) {
            (*normal)[at * n + at] += weightKick;
            (*gradient)[at] -= p[at] * weightKick;
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
          if (!(variance > 0.)) {
            held[at] = 1; // no material: the kick is held at zero
            continue;
          }
          held[at] = 0;
          const double w = 1. / variance;
          chi2 += p[at] * p[at] * w;
          if (normal != nullptr) {
            (*normal)[at * n + at] += w;
            (*gradient)[at] -= p[at] * w;
          }
        }
      }
    }
    scattered.addPrior(p, n, chi2, normal, gradient);
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
      const int at = 3 * topology.vertexOfDaughter(i);
      std::array<double, 6> state{p[at], p[at + 1], p[at + 2], p[base] + p[base + 3], p[base + 1] + p[base + 4], p[base + 2] + p[base + 5]};
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
  const bool constrained = constraints.count(topology) > 0 || walkEnds;
  std::vector<double> masses(static_cast<std::size_t>(count));
  std::vector<double> charges(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    masses[i] = daughters[i].propagator.massGeV;
    charges[i] = start.momenta[i][2] > 0. ? 1. : -1.;
  }
  const double fieldScale = static_cast<double>(daughters[0].propagator.bzKiloGauss) * static_cast<double>(constants::math::B2C);
  std::vector<double> values;
  std::vector<double> rows;
  const auto constraintsAt = [&](const std::vector<double>& p) {
    std::vector<std::array<double, 3>> momenta(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
      momenta[i] = {p[firstOf[i]], p[firstOf[i] + 1], p[firstOf[i] + 2]};
    }
    std::vector<std::array<double, 3>> vertices(static_cast<std::size_t>(vertexCount));
    for (int w = 0; w < vertexCount; ++w) {
      vertices[w] = {p[3 * w], p[(3 * w) + 1], p[(3 * w) + 2]};
    }
    std::vector<double> core;
    if (!evaluateConstraints(constraints, topology, vertices, momenta, masses, charges, fieldScale, values, core)) {
      return false;
    }
    const int width = (3 * vertexCount) + (3 * count);
    const int m = static_cast<int>(values.size());
    rows.assign(static_cast<std::size_t>(m * n), 0.);
    for (int a = 0; a < m; ++a) {
      for (int c = 0; c < 3 * vertexCount; ++c) {
        rows[a * n + c] = core[a * width + c];
      }
      for (int i = 0; i < count; ++i) {
        for (int c = 0; c < 3; ++c) {
          rows[a * n + firstOf[i] + c] = core[a * width + (3 * vertexCount) + (3 * i) + c];
        }
      }
    }
    // A walk's end at the point, from the evaluation at \a p (after the
    // position, which no mother's fit takes); the point's equations, its or
    // the production constraint's, smeared by their nuisances.
    scattered.row = walkEnds ? static_cast<int>(values.size()) : (constraints.mass ? 1 : 0);
    values.insert(values.end(), terminalValues.begin(), terminalValues.end());
    rows.insert(rows.end(), terminalRows.begin(), terminalRows.end());
    scattered.apply(p, n, values, rows);
    // A held kick is no parameter: none of its columns (`holdFixed`).
    for (std::size_t a = 0; a < values.size(); ++a) {
      for (int k = 0; k < n; ++k) {
        rows[(a * n) + k] = held[k] != 0 ? 0. : rows[(a * n) + k];
      }
    }
    return true;
  };

  std::vector<double> normal(static_cast<std::size_t>(n * n));
  std::vector<double> gradient(static_cast<std::size_t>(n));
  std::vector<double> factor;
  double chi2 = 0.;
  int measurementCount = 0;
  for (const auto& daughter : daughters) {
    measurementCount += 2 * static_cast<int>(daughter.planes.size());
  }
  for (const auto& mother : mothers) {
    measurementCount += 2 * (static_cast<int>(mother.rows.planes.size()) + (endsAtPoint(mother) ? 1 : 0));
  }
  const int ndf = measurementCount - ((3 * vertexCount) + (3 * count)) + constraints.count(topology);
  // Evaluated at \a p with the material along it, the held kicks at zero.
  const auto evaluateHeld = [&](std::vector<double>& p) {
    if (!evaluate(p, true, chi2, &normal, &gradient)) {
      return false;
    }
    const bool moved = std::ranges::any_of(std::views::iota(0, n), [&](int k) { return held[k] != 0 && p[k] != 0.; });
    if (moved) {
      for (int k = 0; k < n; ++k) {
        p[k] = held[k] != 0 ? 0. : p[k];
      }
      if (!evaluate(p, true, chi2, &normal, &gradient)) {
        return false;
      }
    }
    holdFixed(p, normal, gradient);
    return true;
  };
  // The material model along the trajectory is the model: minimised under
  // it, evaluated anew along the solution, until that no longer moves the
  // solution by what the inputs resolve (the minimiser takes no step).

  const auto frozen = [&](const std::vector<double>& p, double& value, double& resolution, std::vector<double>* matrix, std::vector<double>* rhs) {
    const bool ok = evaluate(p, false, value, matrix, rhs);
    resolution = rounding;
    return ok;
  };
  double mu = 0.;
  double metric = 0.;
  int iterations = 0;
  int steps = 0;
  double shift = 0.;
  double lastShift = std::numeric_limits<double>::infinity();
  result.status = DecayFitStatus::kNotConverged;
  for (int outer = 0; outer < kMaxSteps; ++outer) {
    stopped = false;
    if (!evaluateHeld(parameters)) {
      result.status = stopped ? DecayFitStatus::kStopped : DecayFitStatus::kNoCrossing;
      return result;
    }
    result.outer = static_cast<std::uint8_t>(outer + 1);
    const auto inner = minimiseFrozen(n, ndf, constrained, frozen, constraintsAt, holdFixed, values, rows, parameters, chi2, mu, metric, iterations, steps, shift);
    if (inner != DecayFitStatus::kOk) {
      result.status = inner;
      break;
    }
    result.modelShift = static_cast<float>(shift);
    if (steps == 0) {
      // The model anew along the solution moves it less than the inputs
      // resolve: converged.
      result.status = DecayFitStatus::kOk;
      result.stopReason = 1U;
      break;
    }
    if (outer > 1 && shift >= lastShift) {
      // The model's moves no longer shrink: the material table's
      // discretisation (piecewise constant along the trajectory) leaves no
      // exact fixed point. Each solution here is a converged minimum under a
      // model along a trajectory close to it; taken if the model moves it by
      // less than its error (`kModelTolerance`).
      result.status = shift <= kModelTolerance ? DecayFitStatus::kOk : DecayFitStatus::kNotConverged;
      result.stopReason = shift <= kModelTolerance ? 2U : 0U;
      break;
    }
    if (outer != 0) {
      lastShift = shift;
    }
  }
  result.iterations = static_cast<std::uint8_t>(std::min(iterations, 255));
  result.stopMetric = static_cast<float>(metric);
  if (!evaluateHeld(parameters)) {
    result.status = DecayFitStatus::kNoCrossing;
    return result;
  }
  factor = normal;
  if (!constrained && !cholesky(factor, n)) {
    result.status = DecayFitStatus::kSingular;
    return result;
  }
  // The covariance of V and the momenta: the full inverse's rows of them.
  std::vector<double> inverse;
  std::vector<double> sensitivity;
  if (constrained) {
    std::vector<double> unused;
    if (!constraintsAt(parameters) || !constrainedStep(normal, n, gradient, rows, values, unused, &inverse, &sensitivity)) {
      result.status = DecayFitStatus::kSingular;
      return result;
    }
  } else {
    inverse = choleskyInverse(factor, n);
  }
  const int m = (3 * vertexCount) + (3 * count);
  std::vector<int> kept;
  for (int c = 0; c < 3 * vertexCount; ++c) {
    kept.push_back(c);
  }
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
  if (constraints.position) {
    // the position is the last constraint (`evaluateConstraints`)
    const int constraintCount = static_cast<int>(values.size());
    result.positionSensitivity.resize(static_cast<std::size_t>(m));
    for (int a = 0; a < m; ++a) {
      result.positionSensitivity[a] = sensitivity[kept[a] * constraintCount + (constraintCount - 1)];
    }
  }
  std::copy_n(parameters.begin(), 3, result.vertex.begin());
  result.vertices.resize(static_cast<std::size_t>(vertexCount));
  for (int w = 0; w < vertexCount; ++w) {
    std::copy_n(parameters.begin() + (3 * w), 3, result.vertices[w].begin());
  }
  result.momenta.resize(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    for (int k = 0; k < 3; ++k) {
      result.momenta[i][k] = parameters[firstOf[i] + k];
    }
  }
  result.chi2 = chi2;
  result.ndf = ndf;
  result.carried = start.carried;
  result.carriedCovariance = start.carriedCovariance;
  result.carryX2X0 = start.carryX2X0;
  return result;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYROWSFIT_H_
