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

#ifndef GTRK_TRACKING_DECAYVERTEXFIT_H_
#define GTRK_TRACKING_DECAYVERTEXFIT_H_

#include "GenericTracking/Base/LayerLinks.h"
#include "GenericTracking/Tracking/BeamConstraint.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/DecayFitConstraints.h"
#include "GenericTracking/Tracking/DecayFitMath.h"
#include "GenericTracking/Tracking/GSFTrackFitter.h"
#include "GenericTracking/Tracking/HypothesisRefit.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/TrackFitterUtils.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include "DetectorsBase/MatCell.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <new>
#include <numbers>
#include <numeric>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace o2::gtrk
{

/// A daughter as the fit takes it.
struct DecayDaughter {
  BoundTrackState innermost; ///< fitted, on `surface`; an electron's its mixture's collapse
  TrackingSurface surface;   ///< its innermost measurement's plane
  Propagator propagator;     ///< the field, with the daughter's species' mass
  double mass = 0.;          ///< [GeV]
  /// An electron's Gaussian-sum fit there: carried as a mixture, its
  /// radiation on the way taken (`takeRadiation`), and collapsed on the plane;
  /// no components for any other daughter.
  GSFMixture<BetheHeitlerLoss::kComponents> mixture{};
  /// A neutral daughter: a line, its q/p its 1/p (an unmeasured one, below).
  bool neutral = false;
  /// A particle not seen at all (a neutrino, a photon, a pi0 as a whole),
  /// neutral: no measurement, its momentum fitted through the constraints
  /// (its mother's mass and production point fix three components); its
  /// start the caller's, `innermost`'s phi, theta and 1/p.
  bool unmeasured = false;
  /// A charged daughter's |q| [e]: 2 for 3He, 4He.
  double chargeMagnitude = 1.;
  /// The daughter's charge [e], 0 for a neutral one.
  [[nodiscard]] GPUhdi() double charge() const { return neutral ? 0. : (innermost.qOverP() > 0.F ? chargeMagnitude : -chargeMagnitude); }
};

/// A charged particle of the decay -- the root, or a composite below it --
/// that left rows of its own before it decayed (a tracked Xi or Omega): its
/// rows, the outermost first, as the fit walks them from its decay vertex
/// against the momentum (`DecayVertexFitter`'s chains), any number of them.
struct DecayMother {
  int vertex = 0;                                    ///< the vertex it decayed at: the root's, or a charged composite's
  std::span<const TrackingSurface> planes;           ///< the caller's storage, as the others
  std::span<const std::array<double, 2>> positions;  ///< loc0, loc1 [cm]
  std::span<const std::array<double, 3>> covariance; ///< loc0, loc0 loc1, loc1 [cm^2]
  Propagator propagator;                             ///< the field, with its mass
};

/// A decay vertex fit's result.
struct DecayVertexFitResult {
  DecayFitStatus status = DecayFitStatus::kNoStart;
  std::uint8_t iterations = 0U;
  std::uint8_t outer = 0U;              ///< the model's evaluations (carries, material along the trajectory)
  std::uint8_t stopReason = 0U;         ///< 1 the model anew moved nothing, 2 its moves stopped shrinking within `kModelTolerance`, 0 not converged
  float stopMetric = 0.F;               ///< the last step's predicted chi2 gain
  float stopBound = 0.F;                ///< and what it was held to: the inputs' resolution R (`minimiseFrozen`)
  float modelShift = 0.F;               ///< what the last model evaluation moved the solution by (its first step's predicted gain)
  double pointSigma = 0.;               ///< the production point's equations' sigma across the flight [cm], with the scattering onto it
  std::uint16_t combinations = 0U;      ///< `kGaussianSum`: the component combinations fitted and taken
  std::uint16_t combinationsTried = 0U; ///< and all combinations tried
  static constexpr std::size_t kVertices = decay_fit_detail::kMaxConstrainedVertices;
  static constexpr std::size_t kDaughters = decay_fit_detail::kMaxConstrainedDaughters;
  std::array<double, 3> vertex{};                         ///< the root's, `vertices[0]`
  FixedVector<std::array<double, 3>, kVertices> vertices; ///< every vertex of the tree (DecayTopology)
  /// Per daughter: its charge [e] (|p| = |q| / |q/p|; a neutral one's q/p is 1/p).
  FixedVector<double, kDaughters> charges;
  /// Per daughter: phi, theta, q/p at its vertex.
  FixedVector<std::array<double, 3>, kDaughters> momenta;
  /// The covariance of (the vertices, the daughters' momenta in order),
  /// 3K + 3N square, packed (`covIndex`; `covarianceAt`).
  FixedVector<double, covSize((3 * kVertices) + (3 * kDaughters))> covariance;
  double chi2 = 0.;
  int ndf = 0;
  /// The daughters' last carried states' phi, theta, q/p: their momenta
  /// without the vertex's information.
  FixedVector<std::array<double, 3>, kDaughters> carried;
  /// The covariance of each carried state's phi, theta, q/p, packed.
  FixedVector<std::array<double, 6>, kDaughters> carriedCovariance;
  /// The material [x/X0] each daughter's last carry crossed.
  FixedVector<double, kDaughters> carryX2X0;

  /// Element (\a row, \a column) of `covariance`.
  [[nodiscard]] GPUhdi() double covarianceAt(int row, int column) const { return covariance[static_cast<std::size_t>(covIndex(row, column))]; }

  /// Back to a fresh result, in place (no temporary of this size).
  GPUhdi() void reset()
  {
    status = DecayFitStatus::kNoStart;
    iterations = outer = stopReason = 0U;
    stopMetric = stopBound = modelShift = 0.F;
    pointSigma = 0.;
    combinations = combinationsTried = 0U;
    vertex = {};
    vertices.clear();
    charges.clear();
    momenta.clear();
    covariance.clear();
    chi2 = 0.;
    ndf = 0;
    carried.clear();
    carriedCovariance.clear();
    carryX2X0.clear();
  }
};

namespace decay_fit_detail
{

/// The minimiser's buffers in a fitter's workspace, for n parameters and m
/// constraint equations: the solve's in S, the trial parameters in double.
template <typename S = double>
struct SolveSpace {
  S* normal = nullptr;       ///< n x n, packed
  S* gradient = nullptr;     ///< n
  S* factor = nullptr;       ///< max(n, m), packed: a Cholesky factor
  S* step = nullptr;         ///< n
  double* trial = nullptr;   ///< n
  S* multipliers = nullptr;  ///< m
  LagrangeSpace<S> lagrange; ///< the Lagrange system's
};

/// The Gauss-Newton minimum under a fixed model -- the carried measurements,
/// the material along a trajectory --: smooth; the model evaluated in M, the
/// solves in S, the parameters, chi2, merit and the stopping test in double.
/// Each step solves
/// the normal equations, under hard constraints the Lagrange system; it is
/// halved until it lowers the merit, the chi2 plus 2 mu |c|_1 with mu the
/// largest |lambda| seen (the exact penalty: the constrained step is a descent
/// direction of it for mu >= |lambda|_inf). Converged when the step's
/// predicted gain, step^T N step + mu |c|_1 - lambda.c (the merit's halved
/// first-order decrease), is no larger than what the inputs resolve: they are
/// floats, each measured value m known to epsilon_float |m|, which moves the
/// solution by R = sum (epsilon_float |m|)^2 w in the same metric (\a
/// evaluate's `rounding`); plus the model arithmetic's rounding of the chi2,
/// epsilon_M (chi2 + ndf). A step that no halving makes lower the merit above
/// that is a failure, not a minimum.
///
/// \a evaluate(p, chi2, rounding, normal*, gradient*); \a constraintsAt(p)
/// fills \a values and \a rows (m x n). \a steps: the steps taken; \a
/// first: the first step's predicted gain, what the model moved the
/// solution by. Its buffers are \a space's.
template <typename M, typename S, typename Evaluate, typename ConstraintsAt>
[[nodiscard]] GPUd() DecayFitStatus minimiseFrozen(int n, int ndf, bool constrained, const Evaluate& evaluate, const ConstraintsAt& constraintsAt, std::span<const S> values, const S* rows, const SolveSpace<S>& space, std::span<double> parameters, double& chi2, double& mu, double& metric, double& bound, int& iterations, int& steps, double& first)
{
  using namespace decay_fit_detail;
  S* normal = space.normal;
  S* gradient = space.gradient;
  Cholesky<S> llt{.l = space.factor};
  S* step = space.step;
  S* multipliers = space.multipliers;
  const std::span<double> trial{space.trial, static_cast<std::size_t>(n)};
  const auto violationOf = [&]() {
    double sum = 0.;
    for (const S c : values) {
      sum += o2::gpu::CAMath::Abs(static_cast<double>(c));
    }
    return sum;
  };
  steps = 0;
  for (int iteration = 0; iteration < kMaxSteps; ++iteration) {
    ++iterations;
    double rounding = 0.;
    if (!evaluate(parameters, chi2, rounding, normal, gradient)) {
      return DecayFitStatus::kNoCrossing;
    }
    decay_fit_detail::copyN(gradient, n, step);
    double violation = 0.;
    double lambdaC = 0.;
    if (constrained) {
      if (!constraintsAt(parameters)) {
        return DecayFitStatus::kNoCrossing;
      }
      if (!constrainedStep(space.lagrange, normal, n, gradient, rows, values, step, nullptr, multipliers)) {
        return DecayFitStatus::kSingular;
      }
      for (std::size_t j = 0; j < values.size(); ++j) {
        mu = o2::gpu::CAMath::Max(mu, o2::gpu::CAMath::Abs(static_cast<double>(multipliers[j])));
        lambdaC += static_cast<double>(multipliers[j]) * static_cast<double>(values[j]);
      }
      violation = violationOf();
    } else {
      if (!llt.factor(normal, n)) {
        return DecayFitStatus::kSingular;
      }
      llt.solve(step);
    }
    metric = (mu * violation) - lambdaC;
    for (int a = 0; a < n; ++a) {
      for (int b = 0; b < n; ++b) {
        metric += static_cast<double>(step[a]) * static_cast<double>(normal[covIndex(a, b)]) * static_cast<double>(step[b]);
      }
    }
    if (iteration == 0) {
      first = metric;
    }
    const double resolvable = (static_cast<double>(std::numeric_limits<M>::epsilon()) * (chi2 + o2::gpu::CAMath::Max(ndf, 0))) + rounding;
    bound = resolvable;
    if (metric <= resolvable) {
      return DecayFitStatus::kOk;
    }
    const double merit = chi2 + (2. * mu * violation);
    bool reached = false;
    bool taken = false;
    double scale = 1.;
    for (int halving = 0; halving < kMaxSteps && !taken; ++halving, scale *= 0.5) {
      for (int j = 0; j < n; ++j) {
        trial[j] = parameters[j] + (scale * static_cast<double>(step[j]));
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
    decay_fit_detail::copyN(trial.data(), parameters.size(), parameters.data());
    ++steps;
  }
  return DecayFitStatus::kNotConverged;
}

/// The global position of \a state on \a surface [cm].
[[nodiscard]] GPUhdi() std::array<double, 3> positionOf(const BoundTrackState& state, const TrackingSurface& surface)
{
  std::array<double, 3> position{};
  o2::gtrk::detail::boundToGlobal(surface, static_cast<ExactScalar>(state.loc0()), static_cast<ExactScalar>(state.loc1()), position.data());
  return position;
}

/// An entry of a vertex's start (`treeStart`): a daughter's state, or one
/// made for a composite; a charged one's propagator.
struct StartEntry {
  const BoundTrackState* state = nullptr;
  const TrackingSurface* surface = nullptr;
  const Propagator* propagator = nullptr;
  bool neutral = false;
};

/// The start: where the transverse circles of the first two entries come
/// closest -- their crossing with the smaller z gap, or the point between
/// them on the line of their centres --, at the mean of their z there.
[[nodiscard]] GPUhdi() bool closestApproach(const StartEntry& a, const StartEntry& b, std::array<double, 3>& start)
{
  struct Circle {
    double cx, cy, radius, rho, px, py, pz, cotTheta;
  };
  const auto circleOf = [](const StartEntry& entry) {
    const auto position = positionOf(*entry.state, *entry.surface);
    const double phi = entry.state->phi();
    const double theta = entry.state->theta();
    const double rho = entry.propagator->curvature(entry.state->qOverP(), sinOf(theta));
    return Circle{.cx = position[0] - (sinOf(phi) / rho), .cy = position[1] + (cosOf(phi) / rho), .radius = 1. / o2::gpu::CAMath::Abs(rho), .rho = rho, .px = position[0], .py = position[1], .pz = position[2], .cotTheta = cosOf(theta) / sinOf(theta)};
  };
  const auto one = circleOf(a);
  const auto two = circleOf(b);
  if (!isFiniteOf(one.radius) || !isFiniteOf(two.radius)) {
    return false;
  }
  const double ux = two.cx - one.cx;
  const double uy = two.cy - one.cy;
  const double d = hypotOf(ux, uy);
  if (!(d > 0.)) {
    return false;
  }
  const double ex = ux / d;
  const double ey = uy / d;
  FixedVector<std::array<double, 2>, 2> points;
  if (d > one.radius + two.radius) {
    const double along = one.radius + (0.5 * (d - one.radius - two.radius));
    points.push_back({one.cx + (along * ex), one.cy + (along * ey)});
  } else if (d < o2::gpu::CAMath::Abs(one.radius - two.radius)) {
    const double sign = one.radius > two.radius ? 1. : -1.;
    const double along = 0.5 * ((sign * one.radius) + d + (sign * two.radius));
    points.push_back({one.cx + (along * ex), one.cy + (along * ey)});
  } else {
    const double along = ((one.radius * one.radius) - (two.radius * two.radius) + (d * d)) / (2. * d);
    const double across = sqrtOf(o2::gpu::CAMath::Max(0., (one.radius * one.radius) - (along * along)));
    points.push_back({one.cx + (along * ex) - (across * ey), one.cy + (along * ey) + (across * ex)});
    points.push_back({one.cx + (along * ex) + (across * ey), one.cy + (along * ey) - (across * ex)});
  }
  // z on a circle's helix at a transverse point: the nearer way round from
  // the daughter's innermost position.
  const auto zAt = [](const Circle& circle, const std::array<double, 2>& point) {
    const double from = atan2Of(circle.py - circle.cy, circle.px - circle.cx);
    const double to = atan2Of(point[1] - circle.cy, point[0] - circle.cx);
    const double turn = toPMPiOf(to - from);
    return circle.pz + ((turn / circle.rho) * circle.cotTheta);
  };
  bool found = false;
  double bestGap = std::numeric_limits<double>::infinity();
  for (const auto& point : points) {
    const double z1 = zAt(one, point);
    const double z2 = zAt(two, point);
    if (o2::gpu::CAMath::Abs(z1 - z2) < bestGap) {
      bestGap = o2::gpu::CAMath::Abs(z1 - z2);
      start = {point[0], point[1], 0.5 * (z1 + z2)};
      found = true;
    }
  }
  return found;
}

/// The closest points of two lines (a point and a unit direction each):
/// their middle; none if they are parallel.
[[nodiscard]] GPUhdi() bool closestOfLines(const std::array<double, 3>& a, const std::array<double, 3>& u, const std::array<double, 3>& b, const std::array<double, 3>& w, std::array<double, 3>& start)
{
  const std::array<double, 3> d{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  const double uw = (u[0] * w[0]) + (u[1] * w[1]) + (u[2] * w[2]);
  const double ud = (u[0] * d[0]) + (u[1] * d[1]) + (u[2] * d[2]);
  const double wd = (w[0] * d[0]) + (w[1] * d[1]) + (w[2] * d[2]);
  const double det = 1. - (uw * uw);
  if (!(det > 0.)) {
    return false;
  }
  const double s = ((uw * wd) - ud) / det;
  const double t = (wd - (uw * ud)) / det;
  start = {0.5 * (a[0] + (s * u[0]) + b[0] + (t * w[0])), 0.5 * (a[1] + (s * u[1]) + b[1] + (t * w[1])), 0.5 * (a[2] + (s * u[2]) + b[2] + (t * w[2]))};
  return true;
}

/// A state's direction (unit).
[[nodiscard]] GPUhdi() std::array<double, 3> directionOf(const BoundTrackState& state)
{
  return momentumOf<double>(state.phi(), state.theta(), 1., 1.).vector;
}

/// Each vertex's start, the deepest first, from its entries -- its daughters,
/// and its composites, each from its own start along its daughters' summed
/// momentum with their summed charge --: two charged entries' helices'
/// closest approach (`closestApproach`); otherwise the closest points of two
/// lines, an entry's from its point along its direction. Into \a starts
/// (one per vertex); false if a vertex has none.
[[nodiscard]] GPUhdi() bool treeStart(std::span<const DecayDaughter> daughters, const DecayTopology& topology, std::span<std::array<double, 3>> starts)
{
  const int k = topology.vertices();
  const int count = static_cast<int>(daughters.size());
  for (int w = k - 1; w >= 0; --w) {
    // The vertex's entries: its daughters, then its composites, each as a
    // daughter from its own start along its daughters' summed momentum with
    // their summed charge (a helix if charged, a line if not).
    FixedVector<StartEntry, kMaxConstrainedDaughters + kMaxConstrainedVertices> entries;
    for (int i = 0; i < count; ++i) {
      if (topology.vertexOfDaughter(i) == w && !daughters[i].unmeasured) {
        entries.push_back({.state = &daughters[i].innermost, .surface = &daughters[i].surface, .propagator = &daughters[i].propagator, .neutral = daughters[i].neutral});
      }
    }
    std::array<BoundTrackState, kMaxConstrainedVertices> compositeStates{};
    std::array<TrackingSurface, kMaxConstrainedVertices> compositeSurfaces{};
    int made = 0;
    for (int child = w + 1; child < k; ++child) {
      if (topology.parent[child] != w) {
        continue;
      }
      MomentumSum<double> sum;
      double charge = 0.;
      const DecayDaughter* field = nullptr;
      for (int i = 0; i < count; ++i) {
        if (topology.under(i, child)) {
          const auto& state = daughters[i].innermost;
          sum.add(momentumOf<double>(state.phi(), state.theta(), state.qOverP(), numeratorOf(daughters[i].charge())).vector);
          charge += daughters[i].charge();
          field = daughters[i].neutral ? field : &daughters[i];
        }
      }
      const bool neutral = charge == 0. || field == nullptr;
      const auto direction = directionOfMomentum(sum.vector);
      auto& state = compositeStates[made];
      auto& surface = compositeSurfaces[made];
      ++made;
      surface = beamLinePlane(static_cast<Scalar>(starts[child][0]), static_cast<Scalar>(starts[child][1]), static_cast<Scalar>(direction.phi));
      state.par[kLoc0] = 0.F;
      state.par[kLoc1] = static_cast<Scalar>(starts[child][2]);
      state.par[kPhi] = static_cast<Scalar>(direction.phi);
      state.par[kTheta] = static_cast<Scalar>(direction.theta);
      state.par[kQOverP] = static_cast<Scalar>((neutral ? 1. : charge) / direction.p);
      entries.push_back({.state = &state, .surface = &surface, .propagator = neutral ? nullptr : &field->propagator, .neutral = neutral});
    }
    if (entries.size() < 2) {
      return false;
    }
    // Two charged entries: their helices' closest approach; else two lines.
    FixedVector<int, kMaxConstrainedDaughters + kMaxConstrainedVertices> charged;
    for (std::size_t e = 0; e < entries.size(); ++e) {
      if (!entries[e].neutral) {
        charged.push_back(static_cast<int>(e));
      }
    }
    if (charged.size() >= 2) {
      if (!closestApproach(entries[charged[0]], entries[charged[1]], starts[w])) {
        return false;
      }
      continue;
    }
    if (!closestOfLines(positionOf(*entries[0].state, *entries[0].surface), directionOf(*entries[0].state), positionOf(*entries[1].state, *entries[1].surface), directionOf(*entries[1].state), starts[w])) {
      return false;
    }
  }
  return true;
}

} // namespace decay_fit_detail

/// An unseen neutral's starts for the missing-momentum fit, from the fit of
/// the visible daughters (\a visible, its momenta with \a masses): the
/// mother flies along the line from its production point \a point to that
/// vertex, |P_M| = x on it (u its direction), the unseen one x u - P, P and
/// E the visible ones' momentum and energy: M^2 = (E + E_n)^2 - x^2 is
/// 2 E E_n = M^2 - m^2 - E^2 - P^2 + 2 x (u.P), squared a quadratic in x;
/// its roots with x > 0 that solve it unsquared (the smaller first). Smeared past the touching point (no real root), its
/// nearest. Each start the unseen one's momentum vector [GeV].
[[nodiscard]] GPUhdi() FixedVector<std::array<double, 3>, 2> missingMomentumStarts(const DecayVertexFitResult& visible, std::span<const double> masses, const std::array<double, 3>& point, double motherMass, double unseenMass)
{
  decay_fit_detail::MomentumSum<double> sum;
  for (std::size_t i = 0; i < visible.momenta.size(); ++i) {
    const auto& q = visible.momenta[i];
    const auto momentum = decay_fit_detail::momentumOf(q[0], q[1], q[2], decay_fit_detail::numeratorOf(visible.charges[i]));
    sum.add(momentum.vector, decay_fit_detail::energyOf(momentum, masses[i]));
  }
  const auto& total = sum.vector;
  const double energy = sum.energy;
  std::array<double, 3> u{visible.vertex[0] - point[0], visible.vertex[1] - point[1], visible.vertex[2] - point[2]};
  const double length = hypotOf(u[0], u[1], u[2]);
  if (!(length > 0.)) {
    return {};
  }
  for (int c = 0; c < 3; ++c) {
    u[c] /= length;
  }
  const double p2 = sum.squared();
  const double b = (u[0] * total[0]) + (u[1] * total[1]) + (u[2] * total[2]);
  const double kappa = (motherMass * motherMass) - (unseenMass * unseenMass) - (energy * energy) - p2; // 2 E E_n = kappa + 2 x b
  const double qa = 4. * ((b * b) - (energy * energy));
  const double qb = 4. * ((kappa * b) + (2. * energy * energy * b));
  const double qc = (kappa * kappa) - (4. * energy * energy * (p2 + (unseenMass * unseenMass)));
  const double discriminant = (qb * qb) - (4. * qa * qc);
  FixedVector<double, 2> roots;
  if (qa == 0.) {
    return {};
  }
  if (discriminant >= 0.) {
    for (const double sign : {-1., 1.}) {
      const double x = (-qb + (sign * sqrtOf(discriminant))) / (2. * qa);
      if (x > 0. && kappa + (2. * x * b) >= 0.) {
        roots.push_back(x);
      }
    }
    if (roots.size() == 2 && roots[1] < roots[0]) {
      const double smaller = roots[1];
      roots[1] = roots[0];
      roots[0] = smaller;
    }
  } else if (-qb / (2. * qa) > 0.) {
    roots.push_back(-qb / (2. * qa));
  }
  FixedVector<std::array<double, 3>, 2> starts;
  for (const double x : roots) {
    const std::array<double, 3> q{(x * u[0]) - total[0], (x * u[1]) - total[1], (x * u[2]) - total[2]};
    if (hypotOf(q[0], q[1], q[2]) > 0.) {
      starts.push_back(q);
    }
  }
  return starts;
}

/// An unseen neutral's start from its momentum vector \a start [GeV] (one of
/// `missingMomentumStarts`): \a daughter's phi, theta and 1/p.
GPUhdi() void startUnmeasured(DecayDaughter& daughter, const std::array<double, 3>& start)
{
  const auto direction = decay_fit_detail::directionOfMomentum(start);
  daughter.innermost.par[kPhi] = static_cast<Scalar>(direction.phi);
  daughter.innermost.par[kTheta] = static_cast<Scalar>(direction.theta);
  daughter.innermost.par[kQOverP] = static_cast<Scalar>(1. / direction.p);
}

namespace decay_fit_detail
{

/// Whether the inputs fit the fit's capacities (`kTooLarge`: more vertices
/// or daughters than the constraints take, a tracked particle with more rows
/// than a track of Det holds) and refer only inside the decay
/// (`kInvalidInput`); checked before anything walks them.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() DecayFitStatus validateInputs(std::span<const DecayDaughter> daughters, const DecayConstraints& constraints, const DecayTopology& topology, std::span<const DecayMother> mothers)
{
  const int count = static_cast<int>(daughters.size());
  const int k = topology.vertices();
  if (count > kMaxConstrainedDaughters || k > kMaxConstrainedVertices) {
    return DecayFitStatus::kTooLarge;
  }
  // an input built beyond its fixed capacity kept only what fits
  bool overflowed = topology.vertexOf.overflowed() || topology.parent.overflowed() || topology.massGeV.overflowed() || constraints.compositeMasses.overflowed() || constraints.resonances.overflowed();
  for (const auto& resonance : constraints.resonances) {
    overflowed = overflowed || resonance.daughters.overflowed();
  }
  if (overflowed) {
    return DecayFitStatus::kTooLarge;
  }
  for (const auto& daughter : daughters) {
    if (daughter.mixture.size < 0 || daughter.mixture.size > std::remove_cvref_t<decltype(daughter.mixture)>::kCapacity) {
      return DecayFitStatus::kInvalidInput;
    }
  }
  const int vertexOfCount = static_cast<int>(topology.vertexOf.size());
  if (vertexOfCount != 0 && vertexOfCount != count) {
    return DecayFitStatus::kInvalidInput;
  }
  for (const int v : topology.vertexOf) {
    if (v < 0 || v >= k) {
      return DecayFitStatus::kInvalidInput;
    }
  }
  const int parentCount = static_cast<int>(topology.parent.size());
  for (int v = 0; v < parentCount; ++v) {
    if (v == 0 ? topology.parent[v] != -1 : (topology.parent[v] < 0 || topology.parent[v] >= v)) {
      return DecayFitStatus::kInvalidInput;
    }
  }
  // composite masses are per vertex below the root: k - 1 of them at most
  const int vertexMasses = static_cast<int>(topology.massGeV.size());
  const int compositeMasses = static_cast<int>(constraints.compositeMasses.size());
  if (vertexMasses > k || compositeMasses > k - 1) {
    return DecayFitStatus::kInvalidInput;
  }
  for (const auto& resonance : constraints.resonances) {
    for (const int i : resonance.daughters) {
      if (i < 0 || i >= count) {
        return DecayFitStatus::kInvalidInput;
      }
    }
  }
  for (std::size_t a = 0; a < mothers.size(); ++a) {
    const auto& mother = mothers[a];
    if (mother.vertex < 0 || mother.vertex >= k || mother.positions.size() != mother.planes.size() || mother.covariance.size() != mother.planes.size()) {
      return DecayFitStatus::kInvalidInput;
    }
    for (std::size_t b = 0; b < a; ++b) {
      if (mothers[b].vertex == mother.vertex) {
        return DecayFitStatus::kInvalidInput;
      }
    }
    if (mother.planes.size() > maxChainLayers<Det>) {
      return DecayFitStatus::kTooLarge;
    }
  }
  return DecayFitStatus::kOk;
}

} // namespace decay_fit_detail

/// The vertex of \a daughters (at least two) under \a Variant, with
/// \a material's material between each daughter's innermost measurement and
/// the vertex, and under \a constraints: then every step is the constrained
/// one (DecayFitConstraints.h), taken whole unless a helix misses its plane;
/// a tree (\a topology) has a vertex per composite, each daughter carried to
/// the plane about its own, and the links always constrained. A fitter fits
/// one decay, once (`fit`).
///
/// It allocates nothing, and holds only what is small and fixed: the counts,
/// per vertex and per daughter a few scalars, the chains' descriptions.
/// Everything that grows with the decay -- the parameters (n), the
/// constraint equations (m), the chains' equations (c) and rows, the matrices
/// over them, the carried measurements, the Gaussian sum's components, the
/// chains' model's scratch -- lives in one workspace the caller hands in, of
/// `workspace()`'s bytes, carved by `layout`; the result goes to the caller's
/// `DecayVertexFitResult`. On a device both can be slices of an arena, sized
/// per decay.
///
/// Its arithmetic: ModelScalar evaluates the model from the parameters --
/// the carry's crossings, the predictions, residuals, chain equations,
/// constraint values and all their derivatives (dual numbers), the chains'
/// Jacobians and the scattering they carry --; SolveScalar holds and solves
/// what the model gives -- the measurements' weights, the normal matrix, the
/// gradient, the constraint rows and values, the factors, the step, the
/// multipliers, the covariances, the log-determinants, the Gaussian sum's
/// collapse's scatter --. Always double, O(n) or O(1) on any device: the
/// parameters and their update, the chi2 and merit sums, mu, the stopping
/// test and its bound -- epsilon_ModelScalar (chi2 + ndf) plus the inputs'
/// rounding --, the outer loop's shift. Both double hold the fit's precision
/// for every input; float in either loses it for some (constrained masses of
/// nearly collinear daughters, the covariance), which `DecayPrecisionStudy`
/// shows decay by decay.
template <DetectorLayout Det, DecayFitVariant Variant, typename Material, typename ModelScalar = double, typename SolveScalar = double>
class DecayVertexFitter
{
  using M = ModelScalar;
  using S = SolveScalar;
  static constexpr int kMaxVertices = decay_fit_detail::kMaxConstrainedVertices;
  static constexpr int kMaxDaughters = decay_fit_detail::kMaxConstrainedDaughters;
  static constexpr int kMaxRows = static_cast<int>(maxChainLayers<Det>); ///< a tracked particle's
  using Mixture = decltype(DecayDaughter::mixture);

  // The chains (below): one per charged particle that decayed here -- the
  // root, every charged composite -- and the root if neutral under a
  // production point; their equations two per row and two at the end, and as
  // many nuisances after the momenta (`ScatteredEquations`).
  struct Chain {
    int vertex = 0;                    ///< where it decayed
    const DecayMother* rows = nullptr; ///< its own rows, if it left any
    int end = -2;                      ///< where it was made: its parent's vertex; -1 the production point; -2 nowhere
    double charge = 0.;                ///< its daughters' sum [e]
    double mass = 0.;                  ///< [GeV], for its flight's material
    int first = 0;                     ///< its first equation among the chains'
    int count = 0;                     ///< two per row and two at its end
    int rowCount = 0;                  ///< its rows
    M* loss = nullptr;                 ///< the model: q/p's factor after each row (in the workspace)
    TrackingSurface endPlane;          ///< the model: the plane about its end point
  };
  struct Measurement {
    TrackingSurface plane;
    std::array<double, 5> value{};
    std::array<S, kBoundCovSize> weight{}; ///< the inverse covariance, packed
  };
  // An electron's carried mixture, its components kept apart for the
  // Gaussian sum (`kGaussianSum`): ln of its normalised weight, its
  // measurement, ln det of its covariance.
  struct Component {
    double logWeight = 0.;
    Measurement measurement;
    double logDetCovariance = 0.;
  };
  /// A daughter's kept components, in the workspace (at most its mixture's
  /// capacity).
  struct Components {
    Component* items = nullptr;
    std::size_t count = 0;
    [[nodiscard]] GPUhdi() std::size_t size() const { return count; }
    [[nodiscard]] GPUhdi() bool empty() const { return count == 0U; }
    GPUhdi() void clear() { count = 0U; }
    GPUhdi() void push_back(const Component& component) { items[count++] = component; } // NOLINT(readability-identifier-naming): the std containers' name
    [[nodiscard]] GPUhdi() const Component& operator[](std::size_t i) const { return items[i]; }
  };
  /// The chains' model's scratch, for the longest chain (`refreshChains`).
  struct ChainScratch {
    std::array<M, 36>* step = nullptr;                 ///< each crossing's step Jacobian (rows + 1)
    std::array<M, 12>* local = nullptr;                ///< and its local one (rows + 1)
    o2::gtrk::detail::StepMaterial* segment = nullptr; ///< the segments' material (rows)
    std::array<M, 6>* at = nullptr;                    ///< the states at the scatterers (rows + 1)
    M* straggling = nullptr;                           ///< the loss's straggling (rows)
    M* jac = nullptr;                                  ///< a kink's sensitivity of the equations (3 (2 rows + 2))
  };
  /// Where everything that grows with the decay lives in the workspace.
  struct Space {
    double* parameters = nullptr;          ///< n
    decay_fit_detail::SolveSpace<S> solve; ///< the minimiser's and the solves' buffers
    S* rows = nullptr;                     ///< the constraints' rows of derivatives, m x n
    S* values = nullptr;                   ///< the constraints' values, m
    S* covariance = nullptr;               ///< n x n
    S* chainCovariance = nullptr;          ///< c x c
    S* chainFactor = nullptr;              ///< c x c, the scattered equations' L
    int* factorOrder = nullptr;            ///< c
    M* loss = nullptr;                     ///< the chains' rows' losses
    ChainScratch scratch;
    Measurement* measurements = nullptr; ///< one per daughter
    Component* components = nullptr;     ///< the Gaussian sum's: the electrons' kept components
    Measurement* saved = nullptr;        ///< and the collapsed measurements, while combinations replace them
    double* from = nullptr;              ///< its start (n)
    double* fitParameters = nullptr;     ///< a combination's parameters (n)
    S* fitCovariance = nullptr;          ///< and covariance (n x n)
    S* square = nullptr;                 ///< C C^T (m x m)
    double* collapseMean = nullptr;      ///< the collapse's mean (n)
    S* collapseScatter = nullptr;        ///< and scatter (n x n)
  };

 public:
  /// What a decay's workspace is sized by.
  struct Counts {
    int parameters = 0;     ///< n
    int values = 0;         ///< m, the constraint equations, the chains' among them
    int chainEquations = 0; ///< c
    int daughters = 0;
    int electrons = 0;    ///< the daughters with a mixture
    int chainRows = 0;    ///< all chains' rows
    int longestChain = 0; ///< the most rows of one chain
  };

  GPUd() DecayVertexFitter(std::span<const DecayDaughter> daughters, const Material& material, const DecayConstraints& constraints, const DecayTopology& topology, std::span<const DecayMother> mothers)
    : mDaughters(daughters), mMaterial(material), mConstraints(constraints), mTopology(topology), mMothers(mothers)
  {
  }

  /// The bytes of workspace a decay of \a counts takes (`layout`).
  [[nodiscard]] GPUhdi() static std::size_t workspaceSize(const Counts& counts)
  {
    Space unused;
    return layout(counts, nullptr, unused);
  }

  /// This decay's workspace [bytes]; none for inputs `fit` refuses.
  [[nodiscard]] GPUd() std::size_t workspace()
  {
    if (!prepare()) {
      return 0U;
    }
    return workspaceSize(counts());
  }

  /// The fit into \a result, in \a workspace (at least `workspace()` bytes,
  /// aligned to kAlignment; else refused as `kTooLarge`).
  GPUd() void fit(std::span<std::byte> workspace, DecayVertexFitResult& result)
  {
    using namespace decay_fit_detail;
    mResult = &result;
    result.reset();
    if (!prepare()) {
      result.status = mStatus;
      return;
    }
    result.charges = mCharges;
    if (workspace.size() < workspaceSize(counts()) || reinterpret_cast<std::uintptr_t>(workspace.data()) % kAlignment != 0U) { // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast): the address's alignment
      result.status = DecayFitStatus::kTooLarge;
      return;
    }
    bind(workspace.data());
    if (mCount < 2) {
      return;
    }
    std::array<std::array<double, 3>, kMaxVertices> start{};
    if (!treeStart(mDaughters, mTopology, start)) {
      return;
    }
    mResult->carried.resize(static_cast<std::size_t>(mCount));
    mResult->carriedCovariance.resize(static_cast<std::size_t>(mCount));
    mResult->carryX2X0.resize(static_cast<std::size_t>(mCount));

    // The field a charged mother flies in: a charged daughter's (a neutral
    // pseudo-track's propagator has none).
    for (const auto& daughter : mDaughters) {
      if (!daughter.neutral) {
        mPhysicalField = static_cast<double>(daughter.propagator.bzKiloGauss) * static_cast<double>(constants::math::B2C);
        break;
      }
    }

    // The parameters: the vertices, then each daughter's phi, theta, q/p.
    decay_fit_detail::fillN(mSpace.parameters, mParameterCount, 0.);
    for (int w = 0; w < mVertexCount; ++w) {
      decay_fit_detail::copyN(start[w].data(), 3, mSpace.parameters + (3 * w));
    }
    if (const auto status = carry(parameters()); status != DecayFitStatus::kOk) {
      mResult->status = status;
      return;
    }
    for (int i = 0; i < mCount; ++i) {
      for (int c = 0; c < 3; ++c) {
        mSpace.parameters[(3 * mVertexCount) + (3 * i) + c] = mResult->carried[i][c];
      }
    }
    // The chain's model, along the start (it needs the momenta).
    if (mChainCount > 0) {
      if (const auto status = refreshChains(parameters()); status != DecayFitStatus::kOk) {
        mResult->status = status;
        return;
      }
    }

    mConstrained = mConstraints.count(mTopology) > 0 || mChainCount > 0;
    mMasses.resize(static_cast<std::size_t>(mCount));
    for (int i = 0; i < mCount; ++i) {
      mMasses[i] = mDaughters[i].mass;
    }
    int measured = 0;
    for (const auto& daughter : mDaughters) {
      measured += daughter.unmeasured ? -3 : 2; // five measured against three; none against three
    }
    mNdf = measured - (3 * mVertexCount) + mConstraints.count(mTopology) + (2 * mChainRows);
    if (!modelLoop() || !covariance()) {
      return;
    }
    if constexpr (Variant == DecayFitVariant::kGaussianSum) {
      gaussianSum();
    }
    decay_fit_detail::copyN(mSpace.parameters, 3, mResult->vertex.begin());
    mResult->vertices.resize(static_cast<std::size_t>(mVertexCount));
    for (int w = 0; w < mVertexCount; ++w) {
      decay_fit_detail::copyN(mSpace.parameters + (3 * w), 3, mResult->vertices[w].begin());
    }
    mResult->momenta.resize(static_cast<std::size_t>(mCount));
    for (int i = 0; i < mCount; ++i) {
      for (int c = 0; c < 3; ++c) {
        mResult->momenta[i][c] = mSpace.parameters[(3 * mVertexCount) + (3 * i) + c];
      }
    }
    // The vertices' and momenta's block, the nuisances left out: packed by
    // rows, it is the covariance's first entries.
    mResult->covariance.resize(static_cast<std::size_t>(covSize(mCore)));
    decay_fit_detail::copyN(mSpace.covariance, covSize(mCore), mResult->covariance.data());
    mResult->chi2 = mChi2;
    mResult->ndf = mNdf;
  }

  /// After `fit`: its parameters (n: the vertices, the momenta, the chains'
  /// nuisances), as the solution or the Gaussian sum's collapse left them.
  [[nodiscard]] GPUd() std::span<const double> solution() const { return {mSpace.parameters, static_cast<std::size_t>(mParameterCount)}; }

  /// After `fit`: the chi2 at \a p (n parameters, as `solution`) under the
  /// fit's last model -- its last carry's measurements, its chains' last
  /// model --, and the constraints' violation |c|_1 there; false if a helix
  /// misses its plane.
  [[nodiscard]] GPUd() bool frozenMerit(std::span<const double> p, double& chi2, double& violation)
  {
    double rounding = 0.;
    violation = 0.;
    if (!evaluate(p, chi2, rounding, nullptr, nullptr)) {
      return false;
    }
    if (mConstrained) {
      if (!constraintsAt(p)) {
        return false;
      }
      for (const S c : values()) {
        violation += o2::gpu::CAMath::Abs(static_cast<double>(c));
      }
    }
    return true;
  }

 private:
  static constexpr int kDerivatives = 6;                                     // vx, vy, vz, phi, theta, q/p
  static constexpr double kEpsilon = std::numeric_limits<Scalar>::epsilon(); // the measurements' precision

  /// The inputs checked, the counts and the chains set up; false (and
  /// `mResult`'s status) for inputs the fit refuses. Once.
  [[nodiscard]] GPUd() bool prepare()
  {
    using namespace decay_fit_detail;
    if (mPrepared) {
      return mStatus == DecayFitStatus::kOk;
    }
    mPrepared = true;
    mStatus = validateInputs<Det>(mDaughters, mConstraints, mTopology, mMothers);
    if (mStatus != DecayFitStatus::kOk) {
      return false;
    }
    mCount = static_cast<int>(mDaughters.size());
    mVertexCount = mTopology.vertices();
    mCore = (3 * mVertexCount) + (3 * mCount);
    mCharges.resize(static_cast<std::size_t>(mCount));
    for (int i = 0; i < mCount; ++i) {
      mCharges[i] = mDaughters[i].charge();
    }
    for (int i = 0; i < mCount; ++i) {
      mElectrons += mDaughters[i].mixture.size > 0 && !mDaughters[i].unmeasured ? 1 : 0;
    }
    setUpChains();
    mParameterCount = mCore + mChainCount;
    mScattered = ScatteredEquations<S>{.at = mCore, .count = mChainCount};
    mValueCount = hardRows() + mChainCount;
    return true;
  }

  /// The equations `evaluateConstraints` gives for these constraints and
  /// this tree: the mother's mass, a neutral composite's link (two), a
  /// composite's mass.
  [[nodiscard]] GPUd() int hardRows() const
  {
    int rows = mConstraints.mass ? 1 : 0;
    for (int w = 1; w < mVertexCount; ++w) {
      rows += mTopology.chargeUnder(mCharges, w) == 0. ? 2 : 0;
    }
    return rows + mConstraints.heldMasses();
  }

  /// The alignment the workspace starts at: every array's in it.
  static constexpr std::size_t kAlignment = alignof(std::max_align_t);

  [[nodiscard]] GPUd() Counts counts() const
  {
    return {.parameters = mParameterCount, .values = mValueCount, .chainEquations = mChainCount, .daughters = mCount, .electrons = mElectrons, .chainRows = mChainRows, .longestChain = mLongestChain};
  }

  /// An array of \a count Ts at \a at (aligned for T) of \a base, \a at
  /// moved past it; none without \a base.
  template <typename T>
  [[nodiscard]] GPUhdi() static T* place(std::byte* base, std::size_t& at, std::size_t count)
  {
    at = (at + alignof(T) - 1U) / alignof(T) * alignof(T);
    T* first = base != nullptr ? reinterpret_cast<T*>(base + at) : nullptr; // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast): the workspace's bytes, constructed as Ts below
    if (first != nullptr) {
      for (std::size_t i = 0; i < count; ++i) {
        new (first + i) T{};
      }
    }
    at += count * sizeof(T);
    return first;
  }

  /// The workspace's layout for \a k: its arrays one after the other in
  /// \a base (none: only its size), into \a space; its size [bytes].
  [[nodiscard]] GPUhdi() static std::size_t layout(const Counts& k, std::byte* base, Space& space)
  {
    const auto n = static_cast<std::size_t>(k.parameters);
    const auto m = static_cast<std::size_t>(k.values);
    const auto c = static_cast<std::size_t>(k.chainEquations);
    const std::size_t system = n + m;
    const std::size_t larger = o2::gpu::CAMath::Max(n, m);
    const auto packed = [](std::size_t size) { return (size * (size + 1U)) / 2U; };
    const auto rows = static_cast<std::size_t>(k.longestChain);
    std::size_t at = 0;
    space.parameters = place<double>(base, at, n);
    space.solve.trial = place<double>(base, at, n);
    space.solve.normal = place<S>(base, at, packed(n));
    space.solve.gradient = place<S>(base, at, n);
    space.solve.factor = place<S>(base, at, packed(larger));
    space.solve.step = place<S>(base, at, n);
    space.solve.multipliers = place<S>(base, at, m);
    space.rows = place<S>(base, at, m * n);
    space.values = place<S>(base, at, m);
    space.solve.lagrange.system = place<S>(base, at, packed(system));
    space.solve.lagrange.solution = place<S>(base, at, system);
    space.covariance = place<S>(base, at, packed(n));
    space.chainCovariance = place<S>(base, at, packed(c));
    space.chainFactor = place<S>(base, at, c * c);
    space.solve.lagrange.swaps = place<int>(base, at, system);
    space.solve.lagrange.blocks = place<int>(base, at, system);
    space.factorOrder = place<int>(base, at, c);
    space.loss = place<M>(base, at, static_cast<std::size_t>(k.chainRows));
    if (c > 0U) {
      space.scratch.step = place<std::array<M, 36>>(base, at, rows + 1U);
      space.scratch.local = place<std::array<M, 12>>(base, at, rows + 1U);
      space.scratch.segment = place<o2::gtrk::detail::StepMaterial>(base, at, rows);
      space.scratch.at = place<std::array<M, 6>>(base, at, rows + 1U);
      space.scratch.straggling = place<M>(base, at, rows);
      space.scratch.jac = place<M>(base, at, 3U * ((2U * rows) + 2U));
    }
    space.measurements = place<Measurement>(base, at, static_cast<std::size_t>(k.daughters));
    if constexpr (Variant == DecayFitVariant::kGaussianSum) {
      space.components = place<Component>(base, at, static_cast<std::size_t>(k.electrons) * static_cast<std::size_t>(Mixture::kCapacity));
      space.saved = place<Measurement>(base, at, static_cast<std::size_t>(k.daughters));
      space.from = place<double>(base, at, n);
      space.fitParameters = place<double>(base, at, n);
      space.collapseMean = place<double>(base, at, n);
      space.fitCovariance = place<S>(base, at, packed(n));
      space.square = place<S>(base, at, packed(m));
      space.collapseScatter = place<S>(base, at, packed(n));
    }
    return at;
  }

  /// Carve the workspace at \a base (`layout`), the chains' losses at 1 and
  /// each electron's share of the components.
  GPUd() void bind(std::byte* base)
  {
    static_cast<void>(layout(counts(), base, mSpace));
    int first = 0;
    for (auto& chain : mChains) {
      chain.loss = mSpace.loss + first;
      decay_fit_detail::fillN(chain.loss, chain.rowCount, M(1.));
      first += chain.rowCount;
    }
    int electron = 0; // NOLINT(misc-const-correctness): incremented under if constexpr, which the check does not see
    for (int i = 0; i < mCount; ++i) {
      mComponents[i] = {};
      if constexpr (Variant == DecayFitVariant::kGaussianSum) {
        if (mDaughters[i].mixture.size > 0 && !mDaughters[i].unmeasured) {
          mComponents[i].items = mSpace.components + (electron * Mixture::kCapacity);
          ++electron;
        }
      }
    }
  }

  [[nodiscard]] GPUd() std::span<double> parameters() const { return {mSpace.parameters, static_cast<std::size_t>(mParameterCount)}; }
  [[nodiscard]] GPUd() std::span<const S> values() const { return {mSpace.values, static_cast<std::size_t>(mValueCount)}; }

  [[nodiscard]] GPUd() const DecayMother* rowsAt(int w) const
  {
    using namespace decay_fit_detail;
    for (const auto& mother : mMothers) {
      if (mother.vertex == w) {
        return &mother;
      }
    }
    return nullptr;
  }

  [[nodiscard]] GPUd() double massAt(int w) const
  {
    using namespace decay_fit_detail;
    if (const auto* own = rowsAt(w); own != nullptr) {
      return static_cast<double>(own->propagator.massGeV);
    }
    return static_cast<std::size_t>(w) < mTopology.massGeV.size() ? mTopology.massGeV[w] : (w == 0 ? mConstraints.massGeV : 0.);
  }

  GPUd() void setUpChains()
  {
    if (mConstraints.production || rowsAt(0) != nullptr) {
      mChains.push_back({.vertex = 0, .rows = rowsAt(0), .end = mConstraints.production ? -1 : -2, .charge = mTopology.chargeUnder(mCharges, 0), .mass = massAt(0)});
    }
    for (int w = 1; w < mVertexCount; ++w) {
      if (mTopology.chargeUnder(mCharges, w) != 0.) {
        mChains.push_back({.vertex = w, .rows = rowsAt(w), .end = mTopology.parent[w], .charge = mTopology.chargeUnder(mCharges, w), .mass = massAt(w)});
      }
    }
    for (auto& chain : mChains) {
      const int rows = chain.rows != nullptr ? static_cast<int>(chain.rows->planes.size()) : 0;
      chain.first = mChainCount;
      chain.count = (2 * rows) + (chain.end != -2 ? 2 : 0);
      chain.rowCount = rows;
      mLongestChain = o2::gpu::CAMath::Max(mLongestChain, rows);
      mChainCount += chain.count;
      mChainRows += rows;
    }
  }

  // Each chain: from its vertex V with its daughters' summed momentum and
  // charge, against the momentum through its rows (outermost first), then to
  // the plane about its end point -- the production point C for the root, its
  // parent's vertex for a composite: its link -- (through the point along the
  // beam, normal to the direction at the closest approach to it). Its model,
  // evaluated along the current solution with each carry (`refreshChains`):
  // the mean loss regained after each row, the end's plane, and the
  // equations' covariance S = R + sum_k J_k Q_k J_k^T + S_end: R the rows'
  // own; a scatterer at V and at each row, the Kalman's halves of the
  // segments' material on either side (Highland with the logarithm at the
  // whole path's x/X0; the loss's straggling on q/p), J_k how a kink there
  // moves every later equation; S_end the last segment's scattering onto the
  // end point, exact from the material table's moment (theta0^2 L^2
  // int t^2 dx / X across, over sin^2(theta) along the beam). A neutral
  // particle crosses nothing: S = R. The chains are separate particles: S is
  // block diagonal. The material model takes unit charge (Highland's p, the
  // loss): a doubly charged chain's flight is not modelled.
  [[nodiscard]] GPUd() decay_fit_detail::FreeState<decay_fit_detail::kConstraintInputs, M> chainStart(std::span<const double> p, const Chain& chain) const
  {
    using namespace decay_fit_detail;
    using DC = Dual<kConstraintInputs, M>;
    FreeState<kConstraintInputs, M> state{};
    for (int c = 0; c < 3; ++c) {
      state[c].v = static_cast<M>(p[(3 * chain.vertex) + c]);
      state[c].d[(3 * chain.vertex) + c] = M(1.);
    }
    MomentumSum<DC> sum;
    for (int i = 0; i < mCount; ++i) {
      if (mTopology.under(i, chain.vertex)) {
        const int q = (3 * mVertexCount) + (3 * i);
        sum.add(momentumOfInputs<kConstraintInputs, M>(p[q], p[q + 1], p[q + 2], q, mCharges[i]).vector);
      }
    }
    const auto direction = directionOfMomentum(sum.vector);
    state[3] = direction.phi;
    state[4] = direction.theta;
    state[5] = DC{static_cast<M>(chain.charge)} / direction.p;
    return state;
  }

  // A chain's end point at \a p: C, or its parent's vertex (a parameter).
  [[nodiscard]] GPUd() std::array<decay_fit_detail::Dual<decay_fit_detail::kConstraintInputs, M>, 3> endPoint(std::span<const double> p, const Chain& chain) const
  {
    using namespace decay_fit_detail;
    std::array<Dual<kConstraintInputs, M>, 3> point{};
    for (int c = 0; c < 3; ++c) {
      if (chain.end == -1) {
        point[c].v = static_cast<M>(mConstraints.point[c]);
      } else {
        point[c].v = static_cast<M>(p[(3 * chain.end) + c]);
        point[c].d[(3 * chain.end) + c] = M(1.);
      }
    }
    return point;
  }

  // The plane about \a point for a state at (x, y) with direction phi.
  [[nodiscard]] GPUd() TrackingSurface endPlaneOf(const std::array<double, 3>& point, double x, double y, double phi, double theta, double qOverP) const
  {
    using namespace decay_fit_detail;
    const double curvature = mPhysicalField * qOverP / sinOf(theta);
    return beamLinePlane(static_cast<Scalar>(point[0]), static_cast<Scalar>(point[1]), static_cast<Scalar>(perigeeAzimuth(x - point[0], y - point[1], phi, curvature)));
  }

  // The chains' equations at \a p, predicted less measured, their rows of
  // derivatives by (the vertices, the momenta), under the current model:
  // after the \a written equations already in `mSpace.values` and `mSpace.rows`.
  [[nodiscard]] GPUd() bool chainEquations(std::span<const double> p, int& written)
  {
    using namespace decay_fit_detail;
    bool room = true;
    const auto push = [&](const Dual<kConstraintInputs, M>& value) {
      if (written >= mValueCount) {
        room = false;
        return;
      }
      mSpace.values[written] = static_cast<S>(value.v);
      for (int c = 0; c < mCore; ++c) {
        mSpace.rows[(written * mParameterCount) + c] = static_cast<S>(value.d[c]);
      }
      ++written;
    };
    for (const auto& chain : mChains) {
      auto state = chainStart(p, chain);
      const int rows = chain.rowCount;
      for (int j = 0; j < rows; ++j) {
        FreeState<kConstraintInputs, M> to{};
        std::array<Dual<kConstraintInputs, M>, 2> local{};
        if (!crossToPlane(state, chain.rows->planes[j], mPhysicalField, to, local)) {
          return false;
        }
        push(local[0] + (-chain.rows->positions[j][0]));
        push(local[1] + (-chain.rows->positions[j][1]));
        state = to;
        state[5] = chain.loss[j] * state[5];
      }
      if (chain.end != -2) {
        FreeState<kConstraintInputs, M> to{};
        std::array<Dual<kConstraintInputs, M>, 2> local{};
        if (!crossToPlane(state, chain.endPlane, mPhysicalField, to, local)) {
          return false;
        }
        // the end point's own place on the plane
        const auto point = endPoint(p, chain);
        const auto& plane = chain.endPlane;
        std::array<Dual<kConstraintInputs, M>, 2> target{};
        for (int c = 0; c < 3; ++c) {
          const auto d = point[c] + (-static_cast<double>(plane.origin[c]));
          target[0] = target[0] + (static_cast<double>(plane.axis0[c]) * d);
          target[1] = target[1] + (static_cast<double>(plane.axis1[c]) * d);
        }
        push(local[0] - target[0]);
        push(local[1] - target[1]);
      }
    }
    return room;
  }

  // The chains' model along \a p's trajectories.
  [[nodiscard]] GPUd() DecayFitStatus refreshChains(std::span<const double> p)
  {
    using namespace decay_fit_detail;
    S* covariance = mSpace.chainCovariance;
    decay_fit_detail::fillN(covariance, covSize(mChainCount), S(0.));
    for (auto& chain : mChains) {
      const auto begin = chainStart(p, chain);
      std::array<M, 6> state{};
      for (int c = 0; c < 6; ++c) {
        state[c] = begin[c].v;
      }
      const int rows = chain.rowCount;
      const bool scatters = chain.charge != 0.;
      const int planes = rows + (chain.end != -2 ? 1 : 0);
      std::array<double, 3> end{};
      if (chain.end != -2) {
        const auto point = endPoint(p, chain);
        end = {point[0].v, point[1].v, point[2].v};
      }
      // Each crossing's step and local Jacobians by the state before it, the
      // segments' material and the states at the scatterers.
      auto* step = mSpace.scratch.step;
      auto* local = mSpace.scratch.local;
      auto* segment = mSpace.scratch.segment;
      auto* at = mSpace.scratch.at;
      M* straggling = mSpace.scratch.straggling;
      for (int j = 0; j <= rows; ++j) {
        step[j] = {};
        local[j] = {};
        at[j] = {};
      }
      for (int j = 0; j < rows; ++j) {
        segment[j] = {};
        straggling[j] = M(0.);
      }
      at[0] = state;
      double totalX2X0 = 0.;
      o2::base::MatMoment last{};
      double lastLength = 0.;
      for (int j = 0; j < planes; ++j) {
        const bool isRow = j < rows;
        if (!isRow) {
          chain.endPlane = endPlaneOf(end, state[0], state[1], state[3], state[4], state[5]);
        }
        const auto& plane = isRow ? chain.rows->planes[j] : chain.endPlane;
        FreeState<6, M> from{};
        for (int c = 0; c < 6; ++c) {
          from[c].v = state[c];
          from[c].d[c] = M(1.);
        }
        FreeState<6, M> to{};
        std::array<Dual<6, M>, 2> loc{};
        if (!crossToPlane(from, plane, mPhysicalField, to, loc)) {
          return DecayFitStatus::kNoCrossing;
        }
        for (int r = 0; r < 6; ++r) {
          for (int c = 0; c < 6; ++c) {
            step[j][(r * 6) + c] = to[r].d[c];
          }
        }
        for (int r = 0; r < 2; ++r) {
          for (int c = 0; c < 6; ++c) {
            local[j][(r * 6) + c] = loc[r].d[c];
          }
        }
        const Scalar a[3] = {static_cast<Scalar>(state[0]), static_cast<Scalar>(state[1]), static_cast<Scalar>(state[2])};
        const Scalar b[3] = {static_cast<Scalar>(to[0].v), static_cast<Scalar>(to[1].v), static_cast<Scalar>(to[2].v)};
        if (isRow) {
          M scale = 1.;
          if (scatters) {
            segment[j] = o2::gtrk::detail::stepMaterial(mMaterial, a, b);
            totalX2X0 += static_cast<double>(segment[j].nearerX2X0 + segment[j].fartherX2X0);
            const auto loss = energyLoss(static_cast<Scalar>(to[5].v), static_cast<Scalar>(chain.mass), segment[j].xRho);
            if (!loss.ok) {
              return DecayFitStatus::kStopped;
            }
            scale = loss.momentumScale;
            straggling[j] = loss.varQOverP;
          }
          chain.loss[j] = scale;
          for (int c = 0; c < 6; ++c) {
            state[c] = to[c].v;
          }
          state[5] *= scale;
          at[j + 1] = state;
        } else if (scatters) {
          // the last segment, from the innermost row (or V) to the end point, t from it
          last = mMaterial.getMatMoment(static_cast<Scalar>(end[0]), static_cast<Scalar>(end[1]), static_cast<Scalar>(end[2]), a[0], a[1], a[2]);
          lastLength = hypotOf(state[0] - end[0], state[1] - end[1], state[2] - end[2]);
          totalX2X0 += static_cast<double>(last.x2X0);
        }
      }
      // This chain's block (packed, its lower triangle): its rows' own, then
      // each scatterer's kink carried to every later equation, then the last
      // segment's onto the end.
      const auto entry = [&](int a, int b) -> S& { return covariance[covIndex(chain.first + a, chain.first + b)]; };
      for (int j = 0; j < rows; ++j) {
        const int e = 2 * j;
        entry(e, e) = chain.rows->covariance[j][0];
        entry(e + 1, e) = chain.rows->covariance[j][1];
        entry(e + 1, e + 1) = chain.rows->covariance[j][2];
      }
      if (!scatters) {
        continue;
      }
      const auto theta0 = [&](Scalar x2X0, const std::array<M, 6>& where) {
        return static_cast<double>(highlandTheta0(x2X0, static_cast<Scalar>(totalX2X0), static_cast<Scalar>(o2::gpu::CAMath::Abs(chain.charge / where[5])), static_cast<Scalar>(chain.mass)));
      };
      // scatterer k sits at at[k]: V (k = 0) or after row k - 1
      for (int kink = 0; kink <= rows; ++kink) {
        const auto& where = at[kink];
        M angle = 0.;
        if (kink < rows) {
          const auto nearer = static_cast<M>(theta0(segment[kink].nearerX2X0, where));
          angle += nearer * nearer;
        }
        if (kink > 0) {
          const auto farther = static_cast<M>(theta0(segment[kink - 1].fartherX2X0, where));
          angle += farther * farther;
        }
        const M sinTheta = sinOf(where[4]);
        const std::array<M, 3> variance{angle / (sinTheta * sinTheta), angle, kink > 0 ? straggling[kink - 1] : M(0.)};
        // the kink's sensitivity, 6 x 3, carried through the later crossings
        std::array<M, 18> g{};
        for (int c = 0; c < 3; ++c) {
          g[((3 + c) * 3) + c] = M(1.);
        }
        M* jac = mSpace.scratch.jac;
        decay_fit_detail::fillN(jac, 3 * chain.count, M(0.));
        for (int j = kink; j < planes; ++j) {
          for (int r = 0; r < 2; ++r) {
            for (int c = 0; c < 3; ++c) {
              M value = 0.;
              for (int q = 0; q < 6; ++q) {
                value += local[j][(r * 6) + q] * g[(q * 3) + c];
              }
              jac[(((2 * j) + r) * 3) + c] = value;
            }
          }
          std::array<M, 18> next{};
          for (int r = 0; r < 6; ++r) {
            for (int c = 0; c < 3; ++c) {
              for (int q = 0; q < 6; ++q) {
                next[(r * 3) + c] += step[j][(r * 6) + q] * g[(q * 3) + c];
              }
            }
          }
          if (j < rows) {
            for (int c = 0; c < 3; ++c) {
              next[(5 * 3) + c] *= chain.loss[j];
            }
          }
          g = next;
        }
        for (int a = 0; a < chain.count; ++a) {
          for (int b = 0; b <= a; ++b) {
            for (int c = 0; c < 3; ++c) {
              entry(a, b) += static_cast<S>(jac[(a * 3) + c] * variance[c] * jac[(b * 3) + c]);
            }
          }
        }
      }
      if (chain.end != -2 && last.x2X0 > 0.F) {
        const auto t0 = static_cast<double>(highlandTheta0(last.x2X0, static_cast<Scalar>(totalX2X0), static_cast<Scalar>(o2::gpu::CAMath::Abs(chain.charge / state[5])), static_cast<Scalar>(chain.mass)));
        const double across = t0 * t0 * lastLength * lastLength * static_cast<double>(last.x2X0T2) / static_cast<double>(last.x2X0);
        const double sinTheta = sinOf(static_cast<double>(state[4]));
        const int e = 2 * rows;
        entry(e, e) += static_cast<S>(across);
        entry(e + 1, e + 1) += static_cast<S>(across / (sinTheta * sinTheta));
      }
    }
    for (const auto& chain : mChains) {
      if (chain.end == -1) {
        const int e = chain.first + chain.count - 2;
        mResult->pointSigma = sqrtOf(static_cast<double>(covariance[covIndex(e, e)]));
      }
    }
    // the covariance factored in place (not needed after)
    dense::semidefiniteFactor(dense::PackedView<S>{covariance}, mChainCount, mSpace.factorOrder, dense::MatrixView<S>{mSpace.chainFactor, mChainCount});
    mScattered.factor = mSpace.chainFactor;
    return DecayFitStatus::kOk;
  }

  // Each daughter carried to the perigee plane about its vertex in \a at.
  [[nodiscard]] GPUd() DecayFitStatus carry(std::span<const double> at)
  {
    using namespace decay_fit_detail;
    for (int i = 0; i < mCount; ++i) {
      const auto& daughter = mDaughters[i];
      if (daughter.unmeasured) {
        mResult->carried[i] = {daughter.innermost.phi(), daughter.innermost.theta(), daughter.innermost.qOverP()};
        continue;
      }
      const auto first = 3U * static_cast<std::size_t>(mTopology.vertexOfDaughter(i)); // its vertex's first coordinate
      const auto plane = detail::hypothesis::perigeePlaneOf(daughter.propagator, daughter.innermost, daughter.surface, static_cast<Scalar>(at[first]), static_cast<Scalar>(at[first + 1]));
      BoundTrackState state = daughter.innermost;
      mResult->carryX2X0[i] = detail::hypothesis::segmentRadiation(daughter.propagator, mMaterial, daughter.innermost, daughter.surface, plane);
      std::array<S, kBoundCovSize> factorStorage{};
      if (daughter.mixture.size > 0) {
        // The electron: its mixture branched by the segment's radiation, every
        // component carried, the collapse on the plane.
        auto mixture = daughter.mixture;
        const auto status = detail::hypothesis::carryMixtureInward(daughter.propagator, mMaterial, daughter.surface, plane, static_cast<Scalar>(mResult->carryX2X0[i]), mixture);
        if (status != PropagationStatus::kOk || mixture.size == 0) {
          return DecayFitStatus::kCarryFailed;
        }
        detail::gsf::collapseMixture(mixture, state);
        // its components kept apart for the Gaussian sum
        if constexpr (Variant == DecayFitVariant::kGaussianSum) {
          auto& kept = mComponents[i];
          kept.clear();
          for (int j = 0; j < mixture.size; ++j) {
            const auto& one = mixture.components[j].state;
            Cholesky<S> llt{.l = factorStorage.data()};
            if (!llt.factor(covarianceOf<S>(one).data(), 5)) {
              continue; // a component without a covariance is no measurement
            }
            Component component{.logWeight = static_cast<double>(mixture.components[j].logWeight)}; // normalised by the collapse
            component.measurement.plane = plane;
            for (int r = 0; r < 5; ++r) {
              component.measurement.value[r] = one.par[r];
            }
            component.logDetCovariance = static_cast<double>(llt.logDeterminant());
            llt.inverse(component.measurement.weight.data());
            kept.push_back(component);
          }
        }
      } else if (detail::hypothesis::carryStateInward(daughter.propagator, mMaterial, daughter.surface, plane, state) != PropagationStatus::kOk) {
        return DecayFitStatus::kCarryFailed;
      }
      Cholesky<S> llt{.l = factorStorage.data()};
      if (!llt.factor(covarianceOf<S>(state).data(), 5)) {
        return DecayFitStatus::kSingular;
      }
      auto& m = mSpace.measurements[i];
      m.plane = plane;
      for (int r = 0; r < 5; ++r) {
        m.value[r] = state.par[r];
      }
      llt.inverse(m.weight.data());
      mResult->carried[i] = {state.phi(), state.theta(), state.qOverP()};
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c <= r; ++c) {
          mResult->carriedCovariance[i][covIndex(r, c)] = state.covAt(r + 2, c + 2);
        }
      }
    }
    return DecayFitStatus::kOk;
  }

  // The chi2 at \a p and, with \a normal (n x n) and \a gradient (n), the
  // normal equations' matrix and right-hand side; false if a helix misses
  // its plane.
  [[nodiscard]] GPUd() bool evaluate(std::span<const double> p, double& chi2, double& rounding, S* normal, S* gradient) const
  {
    using namespace decay_fit_detail;
    chi2 = 0.;
    rounding = 0.;
    if (normal != nullptr) {
      decay_fit_detail::fillN(normal, covSize(mParameterCount), S(0.));
      decay_fit_detail::fillN(gradient, mParameterCount, S(0.));
    }
    for (int i = 0; i < mCount; ++i) {
      const auto& daughter = mDaughters[i];
      if (daughter.unmeasured) {
        continue; // no measurement
      }
      const double fieldScale = daughter.neutral ? 0. : static_cast<double>(daughter.propagator.bzKiloGauss) * static_cast<double>(constants::math::B2C);
      using D = Dual<kDerivatives, M>;
      const int q = (3 * mVertexCount) + (3 * i);
      const int w = 3 * mTopology.vertexOfDaughter(i);
      const std::array<int, kDerivatives> index{w, w + 1, w + 2, q, q + 1, q + 2};
      FreeState<kDerivatives, M> from{};
      for (int c = 0; c < kDerivatives; ++c) {
        from[c] = dual::variable<kDerivatives, M>(p[index[c]], c);
      }
      FreeState<kDerivatives, M> to{};
      std::array<D, 2> local{};
      if (!crossToPlane(from, mSpace.measurements[i].plane, fieldScale, to, local)) {
        return false;
      }
      const std::array<D, 5> predicted{local[0], local[1], to[3], to[4], to[5]};
      std::array<M, 5> residual{};
      for (int r = 0; r < 5; ++r) {
        residual[r] = static_cast<M>(mSpace.measurements[i].value[r]) - predicted[r].v;
      }
      residual[2] = toPMPiOf(residual[2]);
      const auto& weight = mSpace.measurements[i].weight;
      for (int r = 0; r < 5; ++r) {
        rounding += ((kEpsilon * mSpace.measurements[i].value[r]) * (kEpsilon * mSpace.measurements[i].value[r])) * static_cast<double>(weight[covIndex(r, r)]);
      }
      // W r in double: the chi2 is a sum in double
      std::array<double, 5> weighted{};
      for (int r = 0; r < 5; ++r) {
        for (int c = 0; c < 5; ++c) {
          weighted[r] += static_cast<double>(weight[covIndex(r, c)]) * static_cast<double>(residual[c]);
        }
        chi2 += static_cast<double>(residual[r]) * weighted[r];
      }
      if (normal == nullptr) {
        continue;
      }
      // H^T W H (its lower triangle: the indices rise with a) and H^T W r
      // over this daughter's six parameters.
      for (int a = 0; a < kDerivatives; ++a) {
        S g = 0.;
        for (int r = 0; r < 5; ++r) {
          g += static_cast<S>(predicted[r].d[a]) * static_cast<S>(weighted[r]);
        }
        gradient[index[a]] += g;
        for (int b = 0; b <= a; ++b) {
          S h = 0.;
          for (int r = 0; r < 5; ++r) {
            for (int c = 0; c < 5; ++c) {
              h += static_cast<S>(predicted[r].d[a]) * weight[covIndex(r, c)] * static_cast<S>(predicted[c].d[b]);
            }
          }
          normal[covIndex(index[a], index[b])] += h;
        }
      }
    }
    // The chains' rows: the inputs' precision (`minimiseFrozen`).
    for (const auto& mother : mMothers) {
      for (std::size_t j = 0; j < mother.planes.size(); ++j) {
        rounding += ((kEpsilon * mother.positions[j][0]) * (kEpsilon * mother.positions[j][0]) / mother.covariance[j][0]) + ((kEpsilon * mother.positions[j][1]) * (kEpsilon * mother.positions[j][1]) / mother.covariance[j][2]);
      }
    }
    mScattered.addPrior(p, chi2, normal, gradient);
    // The resonances' masses (`Resonance`): Breit-Wigner terms.
    for (const auto& resonance : mConstraints.resonances) {
      using DR = Dual<kConstraintInputs, M>;
      MomentumSum<DR> sum;
      for (const int i : resonance.daughters) {
        const int q = (3 * mVertexCount) + (3 * i);
        const auto momentum = momentumOfInputs<kConstraintInputs, M>(p[q], p[q + 1], p[q + 2], q, mCharges[i]);
        sum.add(momentum.vector, energyOf(momentum, mDaughters[i].mass));
      }
      const DR mass2 = sum.massSquared();
      if (!(mass2.v > 0.)) {
        return false;
      }
      const double big = resonance.massGeV;
      const DR u = (mass2 + (-(big * big))) / DR{static_cast<M>(big * resonance.widthGeV)};
      const M u2 = u.v * u.v;
      chi2 += 2. * log1pOf(static_cast<double>(u2));
      const M weight = 0.5;                       // the Fisher information of chi2 / 2 per unit u^2
      const M slope = M(2.) * u.v / (M(1.) + u2); // d(chi2 / 2) / du
      if (normal == nullptr) {
        continue;
      }
      for (int a = 0; a < mCore; ++a) {
        gradient[a] -= static_cast<S>(slope * u.d[a]);
        for (int b = 0; b <= a; ++b) {
          normal[covIndex(a, b)] += static_cast<S>(weight * u.d[a] * u.d[b]);
        }
      }
    }
    return true;
  }

  // The constraints' values and rows at \a p, into `mSpace.values` and `mSpace.rows`
  // (m x n).
  [[nodiscard]] GPUd() bool constraintsAt(std::span<const double> p)
  {
    using namespace decay_fit_detail;
    std::array<std::array<double, 3>, kMaxVertices> vertices{};
    for (int w = 0; w < mVertexCount; ++w) {
      const auto first = 3U * static_cast<std::size_t>(w);
      vertices[w] = {p[first], p[first + 1], p[first + 2]};
    }
    std::array<std::array<double, 3>, kMaxDaughters> momenta{};
    for (int i = 0; i < mCount; ++i) {
      const int q = (3 * mVertexCount) + (3 * i);
      momenta[i] = {p[q], p[q + 1], p[q + 2]};
    }
    int written = 0;
    if (!evaluateConstraints<M>(mConstraints, mTopology, std::span<const std::array<double, 3>>{vertices.data(), static_cast<std::size_t>(mVertexCount)}, std::span<const std::array<double, 3>>{momenta.data(), static_cast<std::size_t>(mCount)}, std::span<const double>{mMasses}, std::span<const double>{mCharges}, mSpace.values, mSpace.rows, mParameterCount, mValueCount, written)) {
      return false;
    }
    mScattered.row = written;
    if (mChainCount > 0 && !chainEquations(p, written)) {
      return false;
    }
    if (written != mValueCount) {
      return false;
    }
    for (int a = 0; a < mValueCount; ++a) {
      decay_fit_detail::fillN(mSpace.rows + (a * mParameterCount) + mCore, static_cast<std::size_t>(mParameterCount - mCore), S(0.));
    }
    if (mChainCount > 0) {
      mScattered.apply(p, mParameterCount, mSpace.values, mSpace.rows);
    }
    return true;
  }

  /// `minimiseFrozen` under the current model, from \a parameters.
  [[nodiscard]] GPUd() DecayFitStatus minimise(std::span<double> parameters, double& chi2, double& mu, double& metric, double& bound, int& iterations, int& steps, double& first)
  {
    const auto evaluateAt = [this](std::span<const double> p, double& value, double& rounding, S* normal, S* gradient) { return evaluate(p, value, rounding, normal, gradient); };
    const auto constraintsOf = [this](std::span<const double> p) { return constraintsAt(p); };
    return decay_fit_detail::minimiseFrozen<M>(mParameterCount, mNdf, mConstrained, evaluateAt, constraintsOf, values(), mSpace.rows, mSpace.solve, parameters, chi2, mu, metric, bound, iterations, steps, first);
  }

  // The carried measurements are the model: minimised under them, carried
  // anew to the planes about the solution, until a carry no longer moves the
  // solution by what the inputs resolve (the minimiser takes no step).
  /// False if the fit ends here (`mResult` is final).
  [[nodiscard]] GPUd() bool modelLoop()
  {
    using namespace decay_fit_detail;
    double lastShift = std::numeric_limits<double>::infinity();
    mResult->status = DecayFitStatus::kNotConverged;
    for (int outer = 0; outer < kMaxSteps; ++outer) {
      if (outer != 0) {
        if (const auto status = carry(parameters()); status != DecayFitStatus::kOk) {
          mResult->status = status;
          return false;
        }
        if (mChainCount > 0) {
          if (const auto status = refreshChains(parameters()); status != DecayFitStatus::kOk) {
            mResult->status = status;
            return false;
          }
        }
      }
      mResult->outer = static_cast<std::uint8_t>(outer + 1);
      const auto inner = minimise(parameters(), mChi2, mMu, mMetric, mBound, mIterations, mSteps, mShift);
      if (inner != DecayFitStatus::kOk) {
        mResult->status = inner;
        break;
      }
      mResult->modelShift = static_cast<float>(mShift);
      if (mSteps == 0) {
        // The model anew along the solution moves it less than the inputs
        // resolve: converged.
        mResult->status = DecayFitStatus::kOk;
        mResult->stopReason = 1U;
        break;
      }
      if (outer > 1 && mShift >= lastShift) {
        // The model's moves no longer shrink: the material table's
        // discretisation (piecewise constant along the trajectory) leaves no
        // exact fixed point. Each solution here is a converged minimum under a
        // model along a trajectory close to it; taken if the model moves it by
        // less than its error (`kModelTolerance`).
        mResult->status = mShift <= kModelTolerance ? DecayFitStatus::kOk : DecayFitStatus::kNotConverged;
        mResult->stopReason = mShift <= kModelTolerance ? 2U : 0U;
        break;
      }
      if (outer != 0) {
        lastShift = mShift;
      }
    }
    mResult->iterations = static_cast<std::uint8_t>(o2::gpu::CAMath::Min(mIterations, 255));
    mResult->stopMetric = static_cast<float>(mMetric);
    mResult->stopBound = static_cast<float>(mBound);
    return true;
  }

  // The covariance at the solution, with the measurements of the last carry,
  // into `mSpace.covariance` (packed).
  /// False if the fit ends here (`mResult` is final).
  [[nodiscard]] GPUd() bool covariance()
  {
    using namespace decay_fit_detail;
    double rounding = 0.;
    if (!evaluate(parameters(), mChi2, rounding, mSpace.solve.normal, mSpace.solve.gradient)) {
      mResult->status = DecayFitStatus::kNoCrossing;
      return false;
    }
    if (mConstrained) {
      if (!constraintsAt(parameters()) || !constrainedStep(mSpace.solve.lagrange, mSpace.solve.normal, mParameterCount, mSpace.solve.gradient, mSpace.rows, values(), nullptr, mSpace.covariance)) {
        mResult->status = DecayFitStatus::kSingular;
        return false;
      }
    } else {
      Cholesky<S> llt{.l = mSpace.solve.factor};
      if (!llt.factor(mSpace.solve.normal, mParameterCount)) {
        mResult->status = DecayFitStatus::kSingular;
        return false;
      }
      llt.inverse(mSpace.covariance);
    }
    return true;
  }

  // The Gaussian sum: the electrons' carried mixtures kept, a fit per
  // combination of their components (one each) from this solution under this
  // model, weighted by the components' weights times the combination's
  // evidence -- exp(-chi2 / 2) times the posterior's volume over the
  // measurements', sqrt(det C) / prod sqrt(det V), C over the constraints'
  // surface (det(Z^T N Z) = |det KKT| / det(C C^T)) --; the result the
  // posterior mixture's mean and covariance, its chi2 the weighted mean.
  GPUd() void gaussianSum()
  {
    using namespace decay_fit_detail;
    const auto n = static_cast<std::size_t>(mParameterCount);
    FixedVector<int, kMaxDaughters> electrons;
    for (int i = 0; i < mCount; ++i) {
      if (!mComponents[i].empty()) {
        electrons.push_back(i);
      }
    }
    if (!electrons.empty()) {
      decay_fit_detail::copyN(mSpace.measurements, mCount, mSpace.saved);
      decay_fit_detail::copyN(mSpace.parameters, n, mSpace.from);
      // Combination k picks, per electron, its component k's digit in the
      // mixed radix of the mixtures' sizes, the first electron's the fastest.
      std::size_t tried = 1U;
      for (const int e : electrons) {
        tried *= mComponents[e].size();
      }
      MixtureCollapse<S> collapse(std::span<double>{mSpace.collapseMean, n}, std::span<S>{mSpace.collapseScatter, static_cast<std::size_t>(covSize(mParameterCount))});
      const std::span<double> fitParameters{mSpace.fitParameters, n};
      for (std::size_t combination = 0; combination < tried; ++combination) {
        double logWeight = 0.;
        std::size_t digits = combination;
        for (const int e : electrons) {
          const auto& component = mComponents[e][digits % mComponents[e].size()];
          digits /= mComponents[e].size();
          mSpace.measurements[e] = component.measurement;
          logWeight += component.logWeight - (0.5 * component.logDetCovariance);
        }
        decay_fit_detail::copyN(mSpace.from, n, mSpace.fitParameters);
        double fitChi2 = 0.;
        double combinationMu = mMu;
        double combinationMetric = 0.;
        double combinationBound = 0.;
        int combinationIterations = 0;
        int combinationSteps = 0;
        double combinationShift = 0.;
        if (minimise(fitParameters, fitChi2, combinationMu, combinationMetric, combinationBound, combinationIterations, combinationSteps, combinationShift) != DecayFitStatus::kOk) {
          continue;
        }
        double unusedRounding = 0.;
        S* matrix = mSpace.solve.normal;
        S* rhs = mSpace.solve.gradient;
        bool ok = evaluate(fitParameters, fitChi2, unusedRounding, matrix, rhs);
        double logVolume = 0.;
        if (ok && mConstrained) {
          S logKKT = 0.;
          ok = constraintsAt(fitParameters) && constrainedStep(mSpace.solve.lagrange, matrix, mParameterCount, rhs, mSpace.rows, values(), nullptr, mSpace.fitCovariance, nullptr, &logKKT);
          // C C^T, m x m, packed.
          const int m = mValueCount;
          S* square = mSpace.square;
          for (int j = 0; j < m; ++j) {
            for (int c = 0; c <= j; ++c) {
              S sum = 0.;
              for (int b = 0; b < mParameterCount; ++b) {
                sum += mSpace.rows[(j * mParameterCount) + b] * mSpace.rows[(c * mParameterCount) + b];
              }
              square[covIndex(j, c)] = sum;
            }
          }
          Cholesky<S> llt{.l = mSpace.solve.factor};
          ok = ok && llt.factor(square, m);
          logVolume = ok ? static_cast<double>(logKKT) - static_cast<double>(llt.logDeterminant()) : 0.;
        } else if (ok) {
          Cholesky<S> llt{.l = mSpace.solve.factor};
          ok = llt.factor(matrix, mParameterCount);
          if (ok) {
            logVolume = static_cast<double>(llt.logDeterminant());
            llt.inverse(mSpace.fitCovariance);
          }
        }
        if (ok) {
          collapse.add(logWeight - (0.5 * fitChi2) - (0.5 * logVolume), fitParameters, std::span<const S>{mSpace.fitCovariance, static_cast<std::size_t>(covSize(mParameterCount))}, fitChi2);
        }
      }
      decay_fit_detail::copyN(mSpace.saved, mCount, mSpace.measurements);
      mResult->combinationsTried = static_cast<std::uint16_t>(o2::gpu::CAMath::Min<std::size_t>(tried, 65535U));
      if (collapse.count == 0) {
        mResult->status = DecayFitStatus::kNotConverged;
      } else {
        decay_fit_detail::copyN(collapse.mean.data(), collapse.mean.size(), mSpace.parameters);
        mChi2 = collapse.meanChi2();
        collapse.covarianceInto(mSpace.covariance);
        mResult->combinations = static_cast<std::uint16_t>(o2::gpu::CAMath::Min(collapse.count, 65535));
      }
    }
  }

  std::span<const DecayDaughter> mDaughters;
  // The fitter lives for one fit over its caller's inputs.
  const Material& mMaterial;            // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
  const DecayConstraints& mConstraints; // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
  const DecayTopology& mTopology;       // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
  std::span<const DecayMother> mMothers;
  DecayVertexFitResult* mResult = nullptr;      ///< the caller's, during `fit`
  DecayFitStatus mStatus = DecayFitStatus::kOk; ///< `prepare`'s verdict
  bool mPrepared = false;
  int mCount = 0;          ///< the daughters
  int mVertexCount = 0;    ///< the tree's vertices
  int mCore = 0;           ///< the vertices' and momenta's parameters, 3K + 3N
  int mParameterCount = 0; ///< with the chains' nuisances after them: n
  int mChainCount = 0;     ///< the chains' equations: c
  int mChainRows = 0;      ///< the chains' rows
  int mLongestChain = 0;   ///< the most rows of one chain
  int mElectrons = 0;      ///< the daughters with a mixture
  int mValueCount = 0;     ///< the constraint equations, the chains' among them: m
  FixedVector<double, kMaxDaughters> mCharges;
  FixedVector<Chain, kMaxVertices> mChains;
  decay_fit_detail::ScatteredEquations<S> mScattered;
  std::array<Components, kMaxDaughters> mComponents{};
  double mPhysicalField = 0.;
  FixedVector<double, kMaxDaughters> mMasses;
  bool mConstrained = false;
  double mChi2 = 0.;
  int mNdf = 0;
  // the minimiser's state, across the model's evaluations
  double mMu = 0.;
  double mMetric = 0.;
  double mBound = 0.;
  int mIterations = 0;
  int mSteps = 0;
  double mShift = 0.;
  Space mSpace; ///< the workspace (`bind`)
};

/// A `DecayVertexFitter`'s fit of one decay, in a workspace of its own.
template <DetectorLayout Det, DecayFitVariant Variant, typename ModelScalar = double, typename SolveScalar = double, typename Material>
[[nodiscard]] DecayVertexFitResult fitDecayVertex(std::span<const DecayDaughter> daughters, const Material& material, const DecayConstraints& constraints = {}, const DecayTopology& topology = {}, std::span<const DecayMother> mothers = {})
{
  DecayVertexFitter<Det, Variant, Material, ModelScalar, SolveScalar> fitter(daughters, material, constraints, topology, mothers);
  std::vector<std::byte> workspace(fitter.workspace());
  DecayVertexFitResult result;
  fitter.fit(workspace, result);
  return result;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_DECAYVERTEXFIT_H_
