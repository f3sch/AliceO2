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

/// \file DecayVertexFit.h
/// \brief One decay vertex fitted from its charged daughters' fitted tracks:
///        the vertex and every daughter's momentum at it.
///
/// The parameters are the vertex V and, per daughter, its direction and q/p
/// at V (phi, theta, q/p; Billoir's choice). A daughter enters as its fitted
/// state on its innermost measurement, which holds all its rows'
/// information, carried inward with its species' mass and the material
/// between (`carryInward`) to the perigee plane about the current vertex: the
/// plane through the vertex's transverse point along the beam, normal to the
/// daughter's direction at its closest approach (`perigeePlaneOf`). Its five
/// parameters there and their covariance are the measurement. The
/// prediction is the helix from V with the daughter's parameters to that
/// plane, exact in a uniform field: V and the plane are close, so no
/// material lies between. Its derivatives come from forward differentiation
/// of the same formula, exact to rounding.
///
/// Gauss-Newton on all parameters at once (one dense normal matrix, 3 + 3N
/// wide): every iteration carries the daughters anew to the plane about the
/// latest vertex, so the material and the transport follow the vertex. A
/// step that raises the chi2 is halved. The start is the closest approach of
/// the first two daughters' helices in the transverse plane.
///
/// The fit is in double: the states come in as floats.

#ifndef GTRK_STUDY_DECAYVERTEXFIT_H_
#define GTRK_STUDY_DECAYVERTEXFIT_H_

#include "GenericTracking/Study/DecayConstraints.h"
#include "GenericTracking/Study/DecayFitMath.h"
#include "GenericTracking/Tracking/BeamConstraint.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/GSFTrackFitter.h"
#include "GenericTracking/Tracking/HypothesisRefit.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/TrackFitterUtils.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include <algorithm>
#include <functional>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <vector>

namespace o2::gtrk::study
{

/// How a decay vertex fit takes its daughters.
enum class DecayFitVariant : std::uint8_t {
  kFull,       ///< carried with their material to the plane about every iteration's vertex
  kCarryOnce,  ///< carried once, to the plane about the start; later iterations re-linearise only
  kNoMaterial, ///< carried every iteration, without material
  kRows,       ///< not this fit's: the rows' fit (DecayRowsFit.h), from kFull's solution
  kSequential, ///< not one fit: a tree's composite fitted first, a pseudo-track in its parent's fit (`pseudoTrackOf`)
  kCount
};
inline constexpr std::size_t kDecayFitVariants = static_cast<std::size_t>(DecayFitVariant::kCount);

/// Why a decay vertex fit stopped.
enum class DecayFitStatus : std::uint8_t {
  kOk,
  kNoStart,     ///< the first two daughters' helices give no start
  kCarryFailed, ///< a daughter could not be carried to its plane
  kNoCrossing,  ///< the prediction's helix does not reach a plane
  kSingular,    ///< the normal matrix or a measurement's covariance is not positive definite
  kNotConverged,
  kNoProfile, ///< no point along a conversion's line gave a fit (DecayFitStudy.h)
  kStopped,   ///< a daughter's trajectory stops in the material (the loss has no solution)
  kCount
};

/// A daughter as the fit takes it.
struct DecayDaughter {
  BoundTrackState innermost; ///< fitted, on `surface`; an electron's its mixture's collapse
  TrackingSurface surface;   ///< its innermost measurement's plane
  Propagator propagator;     ///< the field, with the daughter's species' mass
  double mass = 0.;          ///< [GeV]
  /// An electron's Gaussian-sum fit there: carried as a mixture, its
  /// radiation on the way taken (`takeRadiation`), and collapsed on the plane.
  std::optional<GSFMixture<BetheHeitlerLoss::kComponents>> mixture;
  /// A neutral daughter -- a composite fitted before, as a pseudo-track
  /// (`pseudoTrackOf`): a line, its q/p its 1/p, carried without material.
  bool neutral = false;
  /// The daughter's charge [e], 0 for a neutral one.
  [[nodiscard]] double charge() const { return neutral ? 0. : (innermost.qOverP() > 0.F ? 1. : -1.); }
};

/// A charged composite that left its own rows before it decayed (a tracked Xi
/// or Omega): its rows' fit, the outward pass's end on its last row, carried
/// along its momentum with the material to the plane about its decay vertex,
/// is one more measurement, predicted from the vertex and its daughters'
/// summed momentum there with their summed charge.
struct DecayMother {
  BoundTrackState outermost; ///< its own fit's outward state, on `surface`
  TrackingSurface surface;   ///< its last row's plane
  Propagator propagator;     ///< the field, with its mass
  int vertex = 0;            ///< the vertex it decayed at
};

/// A decay vertex fit's result.
struct DecayVertexFitResult {
  DecayFitStatus status = DecayFitStatus::kNoStart;
  std::uint8_t iterations = 0U;
  std::uint8_t outer = 0U;                     ///< the model's evaluations (carries, material along the trajectory)
  std::uint8_t stopReason = 0U;                ///< 1 the model anew moved nothing, 2 its moves stopped shrinking within `kModelTolerance`, 0 not converged
  float stopMetric = 0.F;                      ///< the last step's predicted chi2 gain
  float modelShift = 0.F;                      ///< what the last model evaluation moved the solution by (its first step's predicted gain)
  std::array<double, 3> vertex{};              ///< the root's, `vertices[0]`
  std::vector<std::array<double, 3>> vertices; ///< every vertex of the tree (DecayTopology)
  /// Per daughter: phi, theta, q/p at its vertex.
  std::vector<std::array<double, 3>> momenta;
  /// The covariance of (the vertices, the daughters' momenta in order),
  /// (3K + 3N)^2.
  std::vector<double> covariance;
  double chi2 = 0.;
  int ndf = 0;
  /// The daughters' last carried states' phi, theta, q/p: their momenta
  /// without the vertex's information.
  std::vector<std::array<double, 3>> carried;
  /// The covariance of each carried state's phi, theta, q/p, 3 x 3.
  std::vector<std::array<double, 9>> carriedCovariance;
  /// The material [x/X0] each daughter's last carry crossed.
  std::vector<double> carryX2X0;
  /// Per daughter: the first vertex's derivatives by its measurement (the
  /// five parameters on its plane), 3 x 5: how a pseudo-track's errors reach
  /// its parent's vertex.
  std::vector<std::array<double, 15>> vertexGain;
  /// With a `LinePosition` constraint: the parameters' (as `covariance`)
  /// derivatives by its s.
  std::vector<double> positionSensitivity;
};

namespace decay_fit_detail
{

/// The Gauss-Newton minimum under a fixed model -- the carried measurements,
/// the material along a trajectory --: smooth, in double. Each step solves
/// the normal equations, under hard constraints the Lagrange system; it is
/// halved until it lowers the merit, the chi2 plus 2 mu |c|_1 with mu the
/// largest |lambda| seen (the exact penalty: the constrained step is a descent
/// direction of it for mu >= |lambda|_inf). Converged when the step's
/// predicted gain, step^T N step + mu |c|_1 - lambda.c (the merit's halved
/// first-order decrease), is no larger than what the inputs resolve: they are
/// floats, each measured value m known to epsilon_float |m|, which moves the
/// solution by R = sum (epsilon_float |m|)^2 w in the same metric (\a
/// evaluate's `rounding`); plus double's rounding of the chi2, epsilon
/// (chi2 + ndf). A step that no halving makes lower the merit above that is a
/// failure, not a minimum.
///
/// \a evaluate(p, chi2, rounding, normal*, gradient*); \a constraintsAt(p)
/// fills \a values and \a rows; \a prepare(p, normal, gradient) after each
/// full evaluation (held parameters). \a steps: the steps taken; \a first:
/// the first step's predicted gain, what the model moved the solution by.
template <typename Evaluate, typename ConstraintsAt, typename Prepare>
[[nodiscard]] DecayFitStatus minimiseFrozen(int n, int ndf, bool constrained, const Evaluate& evaluate, const ConstraintsAt& constraintsAt, const Prepare& prepare, const std::vector<double>& values, const std::vector<double>& rows, std::vector<double>& parameters, double& chi2, double& mu, double& metric, int& iterations, int& steps, double& first)
{
  using namespace decay_fit_detail;
  std::vector<double> normal(static_cast<std::size_t>(n * n));
  std::vector<double> gradient(static_cast<std::size_t>(n));
  std::vector<double> factor;
  std::vector<double> step;
  std::vector<double> multipliers;
  std::vector<double> trial(static_cast<std::size_t>(n));
  const auto violationOf = [&]() {
    double sum = 0.;
    for (const double c : values) {
      sum += std::abs(c);
    }
    return sum;
  };
  steps = 0;
  for (int iteration = 0; iteration < kMaxSteps; ++iteration) {
    ++iterations;
    double rounding = 0.;
    if (!evaluate(parameters, chi2, rounding, &normal, &gradient)) {
      return DecayFitStatus::kNoCrossing;
    }
    prepare(parameters, normal, gradient);
    step = gradient;
    double violation = 0.;
    double lambdaC = 0.;
    if (constrained) {
      if (!constraintsAt(parameters)) {
        return DecayFitStatus::kNoCrossing;
      }
      if (!constrainedStep(normal, n, gradient, rows, values, step, nullptr, nullptr, &multipliers)) {
        return DecayFitStatus::kSingular;
      }
      for (std::size_t j = 0; j < multipliers.size(); ++j) {
        mu = std::max(mu, std::abs(multipliers[j]));
        lambdaC += multipliers[j] * values[j];
      }
      violation = violationOf();
    } else {
      factor = normal;
      if (!cholesky(factor, n)) {
        return DecayFitStatus::kSingular;
      }
      choleskySolve(factor, n, step);
    }
    metric = (mu * violation) - lambdaC;
    for (int a = 0; a < n; ++a) {
      for (int b = 0; b < n; ++b) {
        metric += step[a] * normal[a * n + b] * step[b];
      }
    }
    if (iteration == 0) {
      first = metric;
    }
    const double resolvable = (std::numeric_limits<double>::epsilon() * (chi2 + std::max(ndf, 0))) + rounding;
    if (metric <= resolvable) {
      return DecayFitStatus::kOk;
    }
    const double merit = chi2 + (2. * mu * violation);
    bool reached = false;
    bool taken = false;
    double scale = 1.;
    for (int halving = 0; halving < kMaxSteps && !taken; ++halving, scale *= 0.5) {
      for (int j = 0; j < n; ++j) {
        trial[j] = parameters[j] + (scale * step[j]);
      }
      double trialChi2 = 0.;
      double trialRounding = 0.;
      if (!evaluate(trial, trialChi2, trialRounding, nullptr, nullptr) || (constrained && !constraintsAt(trial))) {
        continue;
      }
      reached = true;
      taken = trialChi2 + (constrained ? 2. * mu * violationOf() : 0.) < merit;
    }
    if (!taken) {
      return reached ? DecayFitStatus::kNotConverged : DecayFitStatus::kNoCrossing;
    }
    parameters = trial;
    ++steps;
  }
  return DecayFitStatus::kNotConverged;
}

/// The global position of \a state on \a surface [cm].
[[nodiscard]] inline std::array<double, 3> positionOf(const BoundTrackState& state, const TrackingSurface& surface)
{
  Scalar position[3]{};
  surfaceToGlobal(surface, state.loc0(), state.loc1(), position);
  return {position[0], position[1], position[2]};
}

/// The start: where the transverse circles of the first two daughters come
/// closest -- their crossing with the smaller z gap, or the point between
/// them on the line of their centres --, at the mean of their z there.
[[nodiscard]] inline std::optional<std::array<double, 3>> closestApproach(const DecayDaughter& a, const DecayDaughter& b)
{
  struct Circle {
    double cx, cy, radius, rho, px, py, pz, cotTheta;
  };
  const auto circleOf = [](const DecayDaughter& daughter) {
    const auto position = positionOf(daughter.innermost, daughter.surface);
    const double phi = daughter.innermost.phi();
    const double theta = daughter.innermost.theta();
    const double rho = daughter.propagator.curvature(daughter.innermost.qOverP(), std::sin(theta));
    return Circle{position[0] - (std::sin(phi) / rho), position[1] + (std::cos(phi) / rho), 1. / std::abs(rho), rho, position[0], position[1], position[2], std::cos(theta) / std::sin(theta)};
  };
  const auto one = circleOf(a);
  const auto two = circleOf(b);
  if (!std::isfinite(one.radius) || !std::isfinite(two.radius)) {
    return std::nullopt;
  }
  const double ux = two.cx - one.cx;
  const double uy = two.cy - one.cy;
  const double d = std::hypot(ux, uy);
  if (!(d > 0.)) {
    return std::nullopt;
  }
  const double ex = ux / d;
  const double ey = uy / d;
  std::vector<std::array<double, 2>> points;
  if (d > one.radius + two.radius) {
    const double along = one.radius + (0.5 * (d - one.radius - two.radius));
    points.push_back({one.cx + (along * ex), one.cy + (along * ey)});
  } else if (d < std::abs(one.radius - two.radius)) {
    const double sign = one.radius > two.radius ? 1. : -1.;
    const double along = 0.5 * ((sign * one.radius) + d + (sign * two.radius));
    points.push_back({one.cx + (along * ex), one.cy + (along * ey)});
  } else {
    const double along = ((one.radius * one.radius) - (two.radius * two.radius) + (d * d)) / (2. * d);
    const double across = std::sqrt(std::max(0., (one.radius * one.radius) - (along * along)));
    points.push_back({one.cx + (along * ex) - (across * ey), one.cy + (along * ey) + (across * ex)});
    points.push_back({one.cx + (along * ex) + (across * ey), one.cy + (along * ey) - (across * ex)});
  }
  // z on a circle's helix at a transverse point: the nearer way round from
  // the daughter's innermost position.
  const auto zAt = [](const Circle& circle, const std::array<double, 2>& point) {
    const double from = std::atan2(circle.py - circle.cy, circle.px - circle.cx);
    const double to = std::atan2(point[1] - circle.cy, point[0] - circle.cx);
    const double turn = std::remainder(to - from, 2. * std::numbers::pi);
    return circle.pz + ((turn / circle.rho) * circle.cotTheta);
  };
  std::optional<std::array<double, 3>> best;
  double bestGap = std::numeric_limits<double>::infinity();
  for (const auto& point : points) {
    const double z1 = zAt(one, point);
    const double z2 = zAt(two, point);
    if (std::abs(z1 - z2) < bestGap) {
      bestGap = std::abs(z1 - z2);
      best = std::array<double, 3>{point[0], point[1], 0.5 * (z1 + z2)};
    }
  }
  return best;
}

/// The closest points of two lines (a point and a unit direction each):
/// their middle; none if they are parallel.
[[nodiscard]] inline std::optional<std::array<double, 3>> closestOfLines(const std::array<double, 3>& a, const std::array<double, 3>& u, const std::array<double, 3>& b, const std::array<double, 3>& w)
{
  const std::array<double, 3> d{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  const double uw = (u[0] * w[0]) + (u[1] * w[1]) + (u[2] * w[2]);
  const double ud = (u[0] * d[0]) + (u[1] * d[1]) + (u[2] * d[2]);
  const double wd = (w[0] * d[0]) + (w[1] * d[1]) + (w[2] * d[2]);
  const double det = 1. - (uw * uw);
  if (!(det > 0.)) {
    return std::nullopt;
  }
  const double s = ((uw * wd) - ud) / det;
  const double t = (wd - (uw * ud)) / det;
  return std::array<double, 3>{0.5 * (a[0] + (s * u[0]) + b[0] + (t * w[0])), 0.5 * (a[1] + (s * u[1]) + b[1] + (t * w[1])), 0.5 * (a[2] + (s * u[2]) + b[2] + (t * w[2]))};
}

/// A daughter's direction (unit) and momentum vector from its innermost state.
[[nodiscard]] inline std::array<double, 3> directionOf(const BoundTrackState& state)
{
  const double st = std::sin(state.theta());
  return {std::cos(state.phi()) * st, std::sin(state.phi()) * st, std::cos(state.theta())};
}

/// Each vertex's start, the deepest first: two charged daughters' helices'
/// closest approach (`closestApproach`); otherwise the closest points of two
/// lines -- a daughter's from its innermost point along its direction (a
/// neutral one's is its line), a composite's from its own start along its
/// daughters' summed momentum.
[[nodiscard]] inline std::optional<std::vector<std::array<double, 3>>> treeStart(std::span<const DecayDaughter> daughters, const DecayTopology& topology)
{
  const int k = topology.vertices();
  const int count = static_cast<int>(daughters.size());
  std::vector<std::array<double, 3>> starts(static_cast<std::size_t>(k));
  for (int w = k - 1; w >= 0; --w) {
    std::vector<int> tracks;
    for (int i = 0; i < count; ++i) {
      if (topology.vertexOfDaughter(i) == w) {
        tracks.push_back(i);
      }
    }
    if (tracks.size() >= 2 && !daughters[tracks[0]].neutral && !daughters[tracks[1]].neutral) {
      const auto start = closestApproach(daughters[tracks[0]], daughters[tracks[1]]);
      if (!start) {
        return std::nullopt;
      }
      starts[w] = *start;
      continue;
    }
    struct Line {
      std::array<double, 3> point;
      std::array<double, 3> direction;
    };
    std::vector<Line> lines;
    for (const auto i : tracks) {
      lines.push_back({positionOf(daughters[i].innermost, daughters[i].surface), directionOf(daughters[i].innermost)});
    }
    for (int child = w + 1; child < k; ++child) {
      if (topology.parent[child] != w) {
        continue;
      }
      std::array<double, 3> total{};
      for (int i = 0; i < count; ++i) {
        if (topology.under(i, child)) {
          const auto u = directionOf(daughters[i].innermost);
          const double p = 1. / std::abs(static_cast<double>(daughters[i].innermost.qOverP()));
          for (int c = 0; c < 3; ++c) {
            total[c] += p * u[c];
          }
        }
      }
      const double norm = std::hypot(total[0], total[1], total[2]);
      lines.push_back({starts[child], {total[0] / norm, total[1] / norm, total[2] / norm}});
    }
    if (lines.size() < 2) {
      return std::nullopt;
    }
    const auto start = closestOfLines(lines[0].point, lines[0].direction, lines[1].point, lines[1].direction);
    if (!start) {
      return std::nullopt;
    }
    starts[w] = *start;
  }
  return starts;
}

} // namespace decay_fit_detail

/// \a state, its covariance \a factors, carried from \a from to \a to along
/// the momentum, as a Kalman step carries a state between two measurements,
/// with \a lut's material on the straight segment between (`carryInward`'s
/// way, the other direction: the loss taken, not undone).
template <typename Material>
[[nodiscard]] PropagationStatus carryOutward(const Propagator& prop, const Material& lut, const TrackingSurface& from, const TrackingSurface& to, BoundTrackState& state, CovarianceFactors& factors)
{
  Scalar start[3]{};
  surfaceToGlobal(from, state.loc0(), state.loc1(), start);
  const Scalar leaving = state.theta();
  SolenoidJacobian jacobian;
  const auto status = jacobianAnalytic(prop, state, from, to, PropagationDirection::kForward, jacobian, &state);
  if (status != PropagationStatus::kOk && status != PropagationStatus::kOutOfBounds) {
    return status;
  }
  state.surface = to.id;
  Scalar end[3]{};
  surfaceToGlobal(to, state.loc0(), state.loc1(), end);
  const auto step = o2::gtrk::detail::stepMaterial(lut, start, end);
  const Scalar momentum = 1.F / o2::gpu::GPUCommonMath::Abs(state.qOverP());
  const auto before = scatteringNoise(highlandTheta0(step.nearerX2X0, step.nearerLogX2X0, momentum, prop.massGeV), leaving);
  return o2::gtrk::detail::transportStep<false>(prop, step, momentum, before, PropagationDirection::kForward, jacobian, state, factors, 1.F);
}

/// A composite fitted before (\a fit, a one-vertex fit) as a neutral
/// pseudo-track of mass \a mass, for the next level's fit: on the plane
/// through its vertex along the beam, normal to its momentum's transverse
/// direction, loc0 0, loc1 the vertex's z, the summed momentum's phi, theta
/// and 1/p; their covariance the fit's carried by forward differentiation
/// (its \a jacobian, 5 x the fit's parameters, if asked).
/// \a propagator's field is set to none: the composite flies in a line.
[[nodiscard]] inline DecayDaughter pseudoTrackOf(const DecayVertexFitResult& fit, double mass, Propagator propagator, std::vector<double>* jacobian = nullptr)
{
  using namespace decay_fit_detail;
  using D = Dual<kConstraintInputs>;
  const int count = static_cast<int>(fit.momenta.size());
  const int n = 3 + (3 * count);
  std::array<D, 3> total{};
  for (int i = 0; i < count; ++i) {
    std::array<D, 3> q{};
    for (int c = 0; c < 3; ++c) {
      q[c].v = fit.momenta[i][c];
      q[c].d[3 + (3 * i) + c] = 1.;
    }
    const D p = D{q[2].v > 0. ? 1. : -1.} / q[2];
    const D pt = p * sin(q[1]);
    total[0] = total[0] + (pt * cos(q[0]));
    total[1] = total[1] + (pt * sin(q[0]));
    total[2] = total[2] + (p * cos(q[1]));
  }
  const D pt = sqrt((total[0] * total[0]) + (total[1] * total[1]));
  const D p = sqrt((pt * pt) + (total[2] * total[2]));
  const D phi = atan2(total[1], total[0]);
  const D theta = atan2(pt, total[2]);
  const D inverse = D{1.} / p;
  // loc0, loc1 where the line through the vertex crosses the plane through
  // its transverse point (normal n along its transverse direction): a shift
  // of the vertex along axis0 moves loc0; along n it slides the crossing
  // along the line, z by -cot(theta) per unit, which leaves loc0.
  const double cotTheta = std::cos(theta.v) / std::sin(theta.v);
  D loc0{};
  loc0.d[0] = -std::sin(phi.v);
  loc0.d[1] = std::cos(phi.v);
  D loc1{fit.vertex[2]};
  loc1.d[0] = -cotTheta * std::cos(phi.v);
  loc1.d[1] = -cotTheta * std::sin(phi.v);
  loc1.d[2] = 1.;
  const std::array<D, 5> state{loc0, loc1, phi, theta, inverse};
  DecayDaughter daughter{.propagator = propagator, .mass = mass, .neutral = true};
  daughter.propagator.bzKiloGauss = 0.F;
  daughter.propagator.massGeV = static_cast<Scalar>(mass);
  daughter.surface = beamLinePlane(static_cast<Scalar>(fit.vertex[0]), static_cast<Scalar>(fit.vertex[1]), static_cast<Scalar>(phi.v));
  for (int r = 0; r < 5; ++r) {
    daughter.innermost.par[r] = static_cast<Scalar>(state[r].v);
    for (int c = 0; c <= r; ++c) {
      double value = 0.;
      for (int a = 0; a < n; ++a) {
        for (int b = 0; b < n; ++b) {
          value += state[r].d[a] * fit.covariance[(a * n) + b] * state[c].d[b];
        }
      }
      daughter.innermost.covAt(r, c) = static_cast<Scalar>(value);
    }
  }
  daughter.innermost.surface = daughter.surface.id;
  if (jacobian != nullptr) {
    jacobian->assign(static_cast<std::size_t>(5 * n), 0.);
    for (int r = 0; r < 5; ++r) {
      for (int a = 0; a < n; ++a) {
        (*jacobian)[(r * n) + a] = state[r].d[a];
      }
    }
  }
  return daughter;
}

/// The vertex of \a daughters (at least two) under \a variant, with
/// \a material's material between each daughter's innermost measurement and
/// the vertex, and under \a constraints: then every step is the constrained
/// one (DecayConstraints.h), taken whole unless a helix misses its plane;
/// a tree (\a topology) has a vertex per composite, each daughter carried to
/// the plane about its own, and the links always constrained.
template <typename Material>
[[nodiscard]] DecayVertexFitResult fitDecayVertex(std::span<const DecayDaughter> daughters, const Material& material, DecayFitVariant variant, const DecayConstraints& constraints = {}, const DecayTopology& topology = {}, std::span<const DecayMother> mothers = {})
{
  using namespace decay_fit_detail;
  constexpr int kDerivatives = 6; // vx, vy, vz, phi, theta, q/p
  DecayVertexFitResult result;
  const int count = static_cast<int>(daughters.size());
  const int k = topology.vertices();
  const int core = (3 * k) + (3 * count);
  // A production point: two nuisance parameters after the momenta, its
  // scattering's (`ScatteredPoint`).
  const int n = core + (constraints.production ? 2 : 0);
  const ScatteredPoint scattered{.at = constraints.production ? core : -1, .row = constraints.mass ? 1 : 0, .factor = ScatteredPoint::factorOf(constraints.pointCovariance)};
  if (count < 2) {
    return result;
  }
  const auto start = treeStart(daughters, topology);
  if (!start) {
    return result;
  }

  struct Measurement {
    PerigeePlane plane;
    std::array<double, 5> value{};
    std::array<double, 25> weight{}; ///< the inverse covariance
  };
  std::vector<Measurement> measurements(static_cast<std::size_t>(count));
  std::vector<Measurement> motherMeasurements(mothers.size());
  result.carried.resize(static_cast<std::size_t>(count));
  result.carriedCovariance.resize(static_cast<std::size_t>(count));
  result.carryX2X0.resize(static_cast<std::size_t>(count));

  // Each daughter carried to the perigee plane about its vertex in \a at.
  const auto carry = [&](const std::vector<double>& at) {
    for (int i = 0; i < count; ++i) {
      const auto& daughter = daughters[i];
      const int w = topology.vertexOfDaughter(i);
      // A neutral pseudo-track is measured on its own plane: a line without
      // material, which the prediction reaches exactly; nothing to carry.
      const auto plane = daughter.neutral ? daughter.surface : detail::hypothesis::perigeePlaneOf(daughter.propagator, daughter.innermost, daughter.surface, static_cast<Scalar>(at[3 * w]), static_cast<Scalar>(at[(3 * w) + 1]));
      BoundTrackState state = daughter.innermost;
      auto factors = factorCovariance(state.cov);
      PropagationStatus status = PropagationStatus::kOk;
      const bool bare = variant == DecayFitVariant::kNoMaterial || daughter.neutral;
      result.carryX2X0[i] = bare ? 0. : detail::hypothesis::segmentRadiation(daughter.propagator, material, daughter.innermost, daughter.surface, plane);
      if (daughter.mixture) {
        // The electron: its mixture branched by the segment's radiation, every
        // component carried, the collapse on the plane.
        auto mixture = *daughter.mixture;
        if (variant != DecayFitVariant::kNoMaterial && (!(result.carryX2X0[i] >= 0.) || !detail::gsf::takeRadiation<BetheHeitlerLoss>(static_cast<Scalar>(result.carryX2X0[i]), PropagationDirection::kBackward, mixture))) {
          return DecayFitStatus::kCarryFailed;
        }
        status = detail::gsf::stepEach(mixture, [&](GSFComponent<CovarianceFactors>& component) {
          return variant == DecayFitVariant::kNoMaterial ? detail::hypothesis::carryInward(daughter.propagator, o2::gtrk::detail::GivenMaterial{}, daughter.surface, plane, component.state, component.factors) : detail::hypothesis::carryInward(daughter.propagator, material, daughter.surface, plane, component.state, component.factors);
        });
        if (status != PropagationStatus::kOk || mixture.size == 0) {
          return DecayFitStatus::kCarryFailed;
        }
        Scalar mode = 0.F;
        detail::gsf::finishMixture(mixture, state, mode);
      } else if (!daughter.neutral) {
        if (bare) {
          status = detail::hypothesis::carryInward(daughter.propagator, o2::gtrk::detail::GivenMaterial{}, daughter.surface, plane, state, factors);
        } else {
          status = detail::hypothesis::carryInward(daughter.propagator, material, daughter.surface, plane, state, factors);
        }
        if (status != PropagationStatus::kOk) {
          return DecayFitStatus::kCarryFailed;
        }
        covarianceFromFactors(factors, state.cov);
      }
      std::vector<double> covariance(25);
      for (int r = 0; r < 5; ++r) {
        for (int c = 0; c < 5; ++c) {
          covariance[r * 5 + c] = state.covAt(r, c);
        }
      }
      if (!cholesky(covariance, 5)) {
        return DecayFitStatus::kSingular;
      }
      const auto inverse = choleskyInverse(covariance, 5);
      auto& m = measurements[i];
      m.plane = perigeePlaneFrame(plane);
      for (int r = 0; r < 5; ++r) {
        m.value[r] = state.par[r];
      }
      std::ranges::copy(inverse, m.weight.begin());
      result.carried[i] = {state.phi(), state.theta(), state.qOverP()};
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          result.carriedCovariance[i][r * 3 + c] = state.covAt(r + 2, c + 2);
        }
      }
    }
    // Each tracked mother carried along its momentum to the plane about its
    // decay vertex.
    for (std::size_t j = 0; j < mothers.size(); ++j) {
      const auto& mother = mothers[j];
      const int w = mother.vertex;
      const auto plane = detail::hypothesis::perigeePlaneOf(mother.propagator, mother.outermost, mother.surface, static_cast<Scalar>(at[3 * w]), static_cast<Scalar>(at[(3 * w) + 1]));
      BoundTrackState state = mother.outermost;
      auto factors = factorCovariance(state.cov);
      const auto status = variant == DecayFitVariant::kNoMaterial ? carryOutward(mother.propagator, o2::gtrk::detail::GivenMaterial{}, mother.surface, plane, state, factors) : carryOutward(mother.propagator, material, mother.surface, plane, state, factors);
      if (status != PropagationStatus::kOk) {
        return DecayFitStatus::kCarryFailed;
      }
      covarianceFromFactors(factors, state.cov);
      std::vector<double> covariance(25);
      for (int r = 0; r < 5; ++r) {
        for (int c = 0; c < 5; ++c) {
          covariance[r * 5 + c] = state.covAt(r, c);
        }
      }
      if (!cholesky(covariance, 5)) {
        return DecayFitStatus::kSingular;
      }
      const auto inverse = choleskyInverse(covariance, 5);
      auto& m = motherMeasurements[j];
      m.plane = perigeePlaneFrame(plane);
      for (int r = 0; r < 5; ++r) {
        m.value[r] = state.par[r];
      }
      std::ranges::copy(inverse, m.weight.begin());
    }
    return DecayFitStatus::kOk;
  };

  // The parameters: the vertices, then each daughter's phi, theta, q/p.
  std::vector<double> parameters(static_cast<std::size_t>(n));
  for (int w = 0; w < k; ++w) {
    std::copy((*start)[w].begin(), (*start)[w].end(), parameters.begin() + (3 * w));
  }
  if (const auto status = carry(parameters); status != DecayFitStatus::kOk) {
    result.status = status;
    return result;
  }
  for (int i = 0; i < count; ++i) {
    for (int c = 0; c < 3; ++c) {
      parameters[(3 * k) + (3 * i) + c] = result.carried[i][c];
    }
  }

  // The chi2 at \a p and, with \a normal and \a gradient, the normal
  // equations' matrix and right-hand side; false if a helix misses its plane.
  // The field a charged mother flies in: a charged daughter's (a neutral
  // pseudo-track's propagator has none).
  double physicalField = 0.;
  for (const auto& daughter : daughters) {
    if (!daughter.neutral) {
      physicalField = static_cast<double>(daughter.propagator.bzKiloGauss) * static_cast<double>(constants::math::B2C);
      break;
    }
  }
  // Each daughter's prediction's derivatives (5 x its 6 parameters) and
  // their indices, kept at the last evaluation for the vertex's gains.
  std::vector<std::array<double, 30>> jacobians(static_cast<std::size_t>(count));
  std::vector<std::array<int, 6>> jacobianIndex(static_cast<std::size_t>(count));
  constexpr double kEpsilon = std::numeric_limits<Scalar>::epsilon(); // the measurements' precision
  const auto evaluate = [&](const std::vector<double>& p, double& chi2, double& rounding, std::vector<double>* normal, std::vector<double>* gradient) {
    chi2 = 0.;
    rounding = 0.;
    if (normal != nullptr) {
      std::ranges::fill(*normal, 0.);
      std::ranges::fill(*gradient, 0.);
    }
    for (int i = 0; i < count; ++i) {
      const auto& daughter = daughters[i];
      const double fieldScale = daughter.neutral ? 0. : static_cast<double>(daughter.propagator.bzKiloGauss) * static_cast<double>(constants::math::B2C);
      using D = Dual<kDerivatives>;
      std::array<D, kDerivatives> in{};
      const int q = (3 * k) + (3 * i);
      const int w = 3 * topology.vertexOfDaughter(i);
      const std::array<int, kDerivatives> index{w, w + 1, w + 2, q, q + 1, q + 2};
      for (int c = 0; c < kDerivatives; ++c) {
        in[c].v = p[index[c]];
        in[c].d[c] = 1.;
      }
      std::array<D, 5> predicted{};
      if (!predictOnPlane(in[0], in[1], in[2], in[3], in[4], in[5], fieldScale, measurements[i].plane, predicted)) {
        return false;
      }
      std::array<double, 5> residual{};
      for (int r = 0; r < 5; ++r) {
        residual[r] = measurements[i].value[r] - predicted[r].v;
      }
      residual[2] = std::remainder(residual[2], 2. * std::numbers::pi);
      const auto& weight = measurements[i].weight;
      for (int r = 0; r < 5; ++r) {
        rounding += std::pow(kEpsilon * measurements[i].value[r], 2) * weight[r * 5 + r];
      }
      for (int r = 0; r < 5; ++r) {
        for (int a = 0; a < kDerivatives; ++a) {
          jacobians[i][r * 6 + a] = predicted[r].d[a];
        }
      }
      jacobianIndex[i] = index;
      std::array<double, 5> weighted{};
      for (int r = 0; r < 5; ++r) {
        for (int c = 0; c < 5; ++c) {
          weighted[r] += weight[r * 5 + c] * residual[c];
        }
        chi2 += residual[r] * weighted[r];
      }
      if (normal == nullptr) {
        continue;
      }
      // H^T W H and H^T W r over this daughter's six parameters.
      for (int a = 0; a < kDerivatives; ++a) {
        double g = 0.;
        for (int r = 0; r < 5; ++r) {
          g += predicted[r].d[a] * weighted[r];
        }
        (*gradient)[index[a]] += g;
        for (int b = 0; b < kDerivatives; ++b) {
          double h = 0.;
          for (int r = 0; r < 5; ++r) {
            for (int c = 0; c < 5; ++c) {
              h += predicted[r].d[a] * weight[r * 5 + c] * predicted[c].d[b];
            }
          }
          (*normal)[index[a] * n + index[b]] += h;
        }
      }
    }
    // The tracked mothers: predicted from their vertex and their daughters'
    // summed momentum there, with their summed charge.
    for (std::size_t j = 0; j < mothers.size(); ++j) {
      using DM = Dual<kConstraintInputs>;
      const int w = mothers[j].vertex;
      std::array<DM, 3> v{};
      for (int c = 0; c < 3; ++c) {
        v[c].v = p[(3 * w) + c];
        v[c].d[(3 * w) + c] = 1.;
      }
      std::array<DM, 3> total{};
      double charge = 0.;
      for (int i = 0; i < count; ++i) {
        if (!topology.under(i, w)) {
          continue;
        }
        const int q = (3 * k) + (3 * i);
        std::array<DM, 3> m{};
        for (int c = 0; c < 3; ++c) {
          m[c].v = p[q + c];
          m[c].d[q + c] = 1.;
        }
        charge += daughters[i].charge();
        const DM pi = DM{m[2].v > 0. ? 1. : -1.} / m[2];
        const DM pti = pi * sin(m[1]);
        total[0] = total[0] + (pti * cos(m[0]));
        total[1] = total[1] + (pti * sin(m[0]));
        total[2] = total[2] + (pi * cos(m[1]));
      }
      const DM pt = sqrt((total[0] * total[0]) + (total[1] * total[1]));
      const DM pTotal = sqrt((pt * pt) + (total[2] * total[2]));
      const DM phi = atan2(total[1], total[0]);
      const DM theta = atan2(pt, total[2]);
      const DM qOverP = DM{charge} / pTotal;
      const double motherField = static_cast<double>(mothers[j].propagator.bzKiloGauss) * static_cast<double>(constants::math::B2C);
      std::array<DM, 5> predicted{};
      if (!predictOnPlane(v[0], v[1], v[2], phi, theta, qOverP, motherField, motherMeasurements[j].plane, predicted)) {
        return false;
      }
      std::array<double, 5> residual{};
      for (int r = 0; r < 5; ++r) {
        residual[r] = motherMeasurements[j].value[r] - predicted[r].v;
      }
      residual[2] = std::remainder(residual[2], 2. * std::numbers::pi);
      const auto& weight = motherMeasurements[j].weight;
      for (int r = 0; r < 5; ++r) {
        rounding += std::pow(kEpsilon * motherMeasurements[j].value[r], 2) * weight[r * 5 + r];
      }
      std::array<double, 5> weighted{};
      for (int r = 0; r < 5; ++r) {
        for (int c = 0; c < 5; ++c) {
          weighted[r] += weight[r * 5 + c] * residual[c];
        }
        chi2 += residual[r] * weighted[r];
      }
      if (normal == nullptr) {
        continue;
      }
      for (int a = 0; a < n; ++a) {
        double g = 0.;
        for (int r = 0; r < 5; ++r) {
          g += predicted[r].d[a] * weighted[r];
        }
        (*gradient)[a] += g;
        for (int b = 0; b < n; ++b) {
          double h = 0.;
          for (int r = 0; r < 5; ++r) {
            for (int c = 0; c < 5; ++c) {
              h += predicted[r].d[a] * weight[r * 5 + c] * predicted[c].d[b];
            }
          }
          (*normal)[a * n + b] += h;
        }
      }
    }
    scattered.addPrior(p, n, chi2, normal, gradient);
    return true;
  };

  // The constraints' values and rows at \a p.
  const bool constrained = constraints.count(topology) > 0;
  std::vector<double> masses(static_cast<std::size_t>(count));
  std::vector<double> charges(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    masses[i] = daughters[i].mass;
    charges[i] = daughters[i].charge();
  }
  std::vector<double> values;
  std::vector<double> rows;
  const auto constraintsAt = [&](const std::vector<double>& p) {
    std::vector<std::array<double, 3>> vertices(static_cast<std::size_t>(k));
    for (int w = 0; w < k; ++w) {
      vertices[w] = {p[3 * w], p[(3 * w) + 1], p[(3 * w) + 2]};
    }
    std::vector<std::array<double, 3>> momenta(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
      const int q = (3 * k) + (3 * i);
      momenta[i] = {p[q], p[q + 1], p[q + 2]};
    }
    std::vector<double> coreRows;
    if (!evaluateConstraints(constraints, topology, vertices, momenta, masses, charges, physicalField, values, coreRows)) {
      return false;
    }
    const int m = static_cast<int>(values.size());
    rows.assign(static_cast<std::size_t>(m * n), 0.);
    for (int a = 0; a < m; ++a) {
      std::copy_n(coreRows.begin() + (a * core), core, rows.begin() + (a * n));
    }
    scattered.apply(p, n, values, rows);
    return true;
  };

  std::vector<double> normal(static_cast<std::size_t>(n * n));
  std::vector<double> gradient(static_cast<std::size_t>(n));
  std::vector<double> factor;
  double chi2 = 0.;
  const int ndf = (2 * count) - (3 * k) + constraints.count(topology) + (5 * static_cast<int>(mothers.size()));
  // The carried measurements are the model: minimised under them, carried
  // anew to the planes about the solution, until a carry no longer moves the
  // solution by what the inputs resolve (the minimiser takes no step); once
  // under kCarryOnce.

  const auto noPreparation = [](std::vector<double>&, std::vector<double>&, std::vector<double>&) {};
  double mu = 0.;
  double metric = 0.;
  int iterations = 0;
  int steps = 0;
  double shift = 0.;
  double lastShift = std::numeric_limits<double>::infinity();
  result.status = DecayFitStatus::kNotConverged;
  for (int outer = 0; outer < kMaxSteps; ++outer) {
    if (outer != 0) {
      if (variant == DecayFitVariant::kCarryOnce) {
        result.status = DecayFitStatus::kOk;
        result.stopReason = 1U;
        break;
      }
      if (const auto status = carry(parameters); status != DecayFitStatus::kOk) {
        result.status = status;
        return result;
      }
    }
    result.outer = static_cast<std::uint8_t>(outer + 1);
    const auto inner = minimiseFrozen(n, ndf, constrained, evaluate, constraintsAt, noPreparation, values, rows, parameters, chi2, mu, metric, iterations, steps, shift);
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
  // The covariance at the solution, with the measurements of the last carry.
  double rounding = 0.;
  if (!evaluate(parameters, chi2, rounding, &normal, &gradient)) {
    result.status = DecayFitStatus::kNoCrossing;
    return result;
  }
  factor = normal;
  if (!constrained && !cholesky(factor, n)) {
    result.status = DecayFitStatus::kSingular;
    return result;
  }
  if (constrained) {
    std::vector<double> unused;
    std::vector<double> sensitivity;
    if (!constraintsAt(parameters) || !constrainedStep(normal, n, gradient, rows, values, unused, &result.covariance, &sensitivity)) {
      result.status = DecayFitStatus::kSingular;
      return result;
    }
    if (constraints.position) {
      // the position is the last constraint (`evaluateConstraints`)
      const int m = static_cast<int>(values.size());
      result.positionSensitivity.resize(static_cast<std::size_t>(n));
      for (int r = 0; r < n; ++r) {
        result.positionSensitivity[r] = sensitivity[r * m + (m - 1)];
      }
    }
  } else {
    result.covariance = choleskyInverse(factor, n);
  }
  std::copy_n(parameters.begin(), 3, result.vertex.begin());
  result.vertices.resize(static_cast<std::size_t>(k));
  for (int w = 0; w < k; ++w) {
    std::copy_n(parameters.begin() + (3 * w), 3, result.vertices[w].begin());
  }
  result.momenta.resize(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    for (int c = 0; c < 3; ++c) {
      result.momenta[i][c] = parameters[(3 * k) + (3 * i) + c];
    }
  }
  // Each daughter's gain on the first vertex: d V_0 / d m_i = (C J_i^T W_i)
  // over V_0's rows, with the constrained covariance under constraints.
  result.vertexGain.assign(static_cast<std::size_t>(count), {});
  for (int i = 0; i < count; ++i) {
    const auto& weight = measurements[i].weight;
    for (int v = 0; v < 3; ++v) {
      std::array<double, 5> cj{}; // (C J^T)_{v, s}
      for (int sIndex = 0; sIndex < 5; ++sIndex) {
        for (int a = 0; a < 6; ++a) {
          cj[sIndex] += result.covariance[(v * n) + jacobianIndex[i][a]] * jacobians[i][(sIndex * 6) + a];
        }
      }
      for (int r = 0; r < 5; ++r) {
        double g = 0.;
        for (int sIndex = 0; sIndex < 5; ++sIndex) {
          g += cj[sIndex] * weight[(sIndex * 5) + r];
        }
        result.vertexGain[i][(v * 5) + r] = g;
      }
    }
  }
  // The vertices' and momenta's block: the nuisances left out.
  if (n != core) {
    std::vector<double> block(static_cast<std::size_t>(core * core));
    for (int a = 0; a < core; ++a) {
      std::copy_n(result.covariance.begin() + (a * n), core, block.begin() + (a * core));
    }
    result.covariance = std::move(block);
    if (!result.positionSensitivity.empty()) {
      result.positionSensitivity.resize(static_cast<std::size_t>(core));
    }
  }
  result.chi2 = chi2;
  result.ndf = ndf;
  return result;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYVERTEXFIT_H_
