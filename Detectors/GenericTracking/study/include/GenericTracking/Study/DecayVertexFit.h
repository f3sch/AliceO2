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
#include "GenericTracking/Base/LayerLinks.h"
#include "DetectorsBase/MatCell.h"
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
#include <utility>
#include <numeric>
#include <numbers>
#include <optional>
#include <span>
#include <vector>

namespace o2::gtrk::study
{

inline constexpr std::size_t kDecayFitVariants = static_cast<std::size_t>(DecayFitVariant::kCount);

/// A daughter as the fit takes it.
struct DecayDaughter {
  BoundTrackState innermost; ///< fitted, on `surface`; an electron's its mixture's collapse
  TrackingSurface surface;   ///< its innermost measurement's plane
  Propagator propagator;     ///< the field, with the daughter's species' mass
  double mass = 0.;          ///< [GeV]
  /// An electron's Gaussian-sum fit there: carried as a mixture, its
  /// radiation on the way taken (`takeRadiation`), and collapsed on the plane.
  std::optional<GSFMixture<BetheHeitlerLoss::kComponents>> mixture;
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
  [[nodiscard]] double charge() const { return neutral ? 0. : (innermost.qOverP() > 0.F ? chargeMagnitude : -chargeMagnitude); }
  /// q/p's numerator: |p| = it / |q/p| (a neutral one's q/p is 1/p).
  [[nodiscard]] double momentumNumerator() const { return neutral ? 1. : chargeMagnitude; }
};

/// A charged particle of the decay -- the root, or a composite below it --
/// that left rows of its own before it decayed (a tracked Xi or Omega): its
/// rows, the outermost first, as the fit walks them from its decay vertex
/// against the momentum (`DecayVertexFitter`'s chains), any number of them.
struct DecayMother {
  int vertex = 0; ///< the vertex it decayed at: the root's, or a charged composite's
  std::vector<TrackingSurface> planes;
  std::vector<std::array<double, 2>> positions;  ///< loc0, loc1 [cm]
  std::vector<std::array<double, 3>> covariance; ///< loc0, loc0 loc1, loc1 [cm^2]
  Propagator propagator;                         ///< the field, with its mass
};

/// A decay vertex fit's result.
struct DecayVertexFitResult {
  DecayFitStatus status = DecayFitStatus::kNoStart;
  std::uint8_t iterations = 0U;
  std::uint8_t outer = 0U;                     ///< the model's evaluations (carries, material along the trajectory)
  std::uint8_t stopReason = 0U;                ///< 1 the model anew moved nothing, 2 its moves stopped shrinking within `kModelTolerance`, 0 not converged
  float stopMetric = 0.F;                      ///< the last step's predicted chi2 gain
  float stopBound = 0.F;                       ///< and what it was held to: the inputs' resolution R (`minimiseFrozen`)
  float modelShift = 0.F;                      ///< what the last model evaluation moved the solution by (its first step's predicted gain)
  double pointSigma = 0.;                      ///< the production point's equations' sigma across the flight [cm], with the scattering onto it
  std::uint16_t combinations = 0U;             ///< `kGaussianSum`: the component combinations fitted and taken
  std::uint16_t combinationsTried = 0U;        ///< and all combinations tried
  std::array<double, 3> vertex{};              ///< the root's, `vertices[0]`
  std::vector<std::array<double, 3>> vertices; ///< every vertex of the tree (DecayTopology)
  /// Per daughter: its charge [e] (|p| = |q| / |q/p|; a neutral one's q/p is 1/p).
  std::vector<double> charges;
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
/// fills \a values and \a rows. \a steps: the steps taken; \a first:
/// the first step's predicted gain, what the model moved the solution by.
template <typename Evaluate, typename ConstraintsAt>
[[nodiscard]] DecayFitStatus minimiseFrozen(int n, int ndf, bool constrained, const Evaluate& evaluate, const ConstraintsAt& constraintsAt, const std::vector<double>& values, const std::vector<double>& rows, std::vector<double>& parameters, double& chi2, double& mu, double& metric, double& bound, int& iterations, int& steps, double& first)
{
  using namespace decay_fit_detail;
  std::vector<double> normal(static_cast<std::size_t>(n * n));
  std::vector<double> gradient(static_cast<std::size_t>(n));
  Cholesky llt;
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
    step = gradient;
    double violation = 0.;
    double lambdaC = 0.;
    if (constrained) {
      if (!constraintsAt(parameters)) {
        return DecayFitStatus::kNoCrossing;
      }
      if (!constrainedStep(normal, n, gradient, rows, values, step, nullptr, &multipliers)) {
        return DecayFitStatus::kSingular;
      }
      for (std::size_t j = 0; j < multipliers.size(); ++j) {
        mu = std::max(mu, std::abs(multipliers[j]));
        lambdaC += multipliers[j] * values[j];
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
        metric += step[a] * normal[a * n + b] * step[b];
      }
    }
    if (iteration == 0) {
      first = metric;
    }
    const double resolvable = (std::numeric_limits<double>::epsilon() * (chi2 + std::max(ndf, 0))) + rounding;
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

/// A state's direction (unit).
[[nodiscard]] inline std::array<double, 3> directionOf(const BoundTrackState& state)
{
  return momentumOf<double>(state.phi(), state.theta(), 1., 1.).vector;
}

/// Each vertex's start, the deepest first, from its entries -- its daughters,
/// and its composites, each from its own start along its daughters' summed
/// momentum with their summed charge --: two charged entries' helices'
/// closest approach (`closestApproach`); otherwise the closest points of two
/// lines, an entry's from its point along its direction.
[[nodiscard]] inline std::optional<std::vector<std::array<double, 3>>> treeStart(std::span<const DecayDaughter> daughters, const DecayTopology& topology)
{
  const int k = topology.vertices();
  const int count = static_cast<int>(daughters.size());
  std::vector<std::array<double, 3>> starts(static_cast<std::size_t>(k));
  for (int w = k - 1; w >= 0; --w) {
    std::vector<int> tracks;
    for (int i = 0; i < count; ++i) {
      if (topology.vertexOfDaughter(i) == w && !daughters[i].unmeasured) {
        tracks.push_back(i);
      }
    }
    // The vertex's entries: its daughters, then its composites, each as a
    // daughter from its own start along its daughters' summed momentum with
    // their summed charge (a helix if charged, a line if not).
    std::vector<DecayDaughter> entries;
    for (const auto i : tracks) {
      entries.push_back(daughters[i]);
    }
    for (int child = w + 1; child < k; ++child) {
      if (topology.parent[child] != w) {
        continue;
      }
      std::array<double, 3> total{};
      double charge = 0.;
      const DecayDaughter* field = nullptr;
      for (int i = 0; i < count; ++i) {
        if (topology.under(i, child)) {
          const auto& state = daughters[i].innermost;
          const auto momentum = momentumOf<double>(state.phi(), state.theta(), state.qOverP(), daughters[i].momentumNumerator());
          for (int c = 0; c < 3; ++c) {
            total[c] += momentum.vector[c];
          }
          charge += daughters[i].charge();
          field = daughters[i].neutral ? field : &daughters[i];
        }
      }
      const double p = std::hypot(total[0], total[1], total[2]);
      const double phi = std::atan2(total[1], total[0]);
      DecayDaughter composite{.neutral = charge == 0. || field == nullptr};
      if (!composite.neutral) {
        composite.propagator = field->propagator;
      }
      composite.surface = beamLinePlane(static_cast<Scalar>(starts[child][0]), static_cast<Scalar>(starts[child][1]), static_cast<Scalar>(phi));
      composite.innermost.par[kLoc0] = 0.F;
      composite.innermost.par[kLoc1] = static_cast<Scalar>(starts[child][2]);
      composite.innermost.par[kPhi] = static_cast<Scalar>(phi);
      composite.innermost.par[kTheta] = static_cast<Scalar>(std::atan2(std::hypot(total[0], total[1]), total[2]));
      composite.innermost.par[kQOverP] = static_cast<Scalar>((composite.neutral ? 1. : charge) / p);
      entries.push_back(composite);
    }
    if (entries.size() < 2) {
      return std::nullopt;
    }
    // Two charged entries: their helices' closest approach; else two lines.
    std::vector<const DecayDaughter*> charged;
    for (const auto& entry : entries) {
      if (!entry.neutral) {
        charged.push_back(&entry);
      }
    }
    if (charged.size() >= 2) {
      const auto start = closestApproach(*charged[0], *charged[1]);
      if (!start) {
        return std::nullopt;
      }
      starts[w] = *start;
      continue;
    }
    const auto start = closestOfLines(positionOf(entries[0].innermost, entries[0].surface), directionOf(entries[0].innermost), positionOf(entries[1].innermost, entries[1].surface), directionOf(entries[1].innermost));
    if (!start) {
      return std::nullopt;
    }
    starts[w] = *start;
  }
  return starts;
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
[[nodiscard]] inline std::vector<std::array<double, 3>> missingMomentumStarts(const DecayVertexFitResult& visible, std::span<const double> masses, const std::array<double, 3>& point, double motherMass, double unseenMass)
{
  std::array<double, 3> total{};
  double energy = 0.;
  for (std::size_t i = 0; i < visible.momenta.size(); ++i) {
    const auto& q = visible.momenta[i];
    const auto momentum = decay_fit_detail::momentumOf(q[0], q[1], q[2], decay_fit_detail::numeratorOf(visible.charges[i]));
    energy += std::sqrt((momentum.magnitude * momentum.magnitude) + (masses[i] * masses[i]));
    for (int c = 0; c < 3; ++c) {
      total[c] += momentum.vector[c];
    }
  }
  std::array<double, 3> u{visible.vertex[0] - point[0], visible.vertex[1] - point[1], visible.vertex[2] - point[2]};
  const double length = std::hypot(u[0], u[1], u[2]);
  if (!(length > 0.)) {
    return {};
  }
  for (auto& c : u) {
    c /= length;
  }
  const double p2 = (total[0] * total[0]) + (total[1] * total[1]) + (total[2] * total[2]);
  const double b = (u[0] * total[0]) + (u[1] * total[1]) + (u[2] * total[2]);
  const double kappa = (motherMass * motherMass) - (unseenMass * unseenMass) - (energy * energy) - p2; // 2 E E_n = kappa + 2 x b
  const double qa = 4. * ((b * b) - (energy * energy));
  const double qb = 4. * ((kappa * b) + (2. * energy * energy * b));
  const double qc = (kappa * kappa) - (4. * energy * energy * (p2 + (unseenMass * unseenMass)));
  const double discriminant = (qb * qb) - (4. * qa * qc);
  std::vector<double> roots;
  if (qa == 0.) {
    return {};
  }
  if (discriminant >= 0.) {
    for (const double sign : {-1., 1.}) {
      const double x = (-qb + (sign * std::sqrt(discriminant))) / (2. * qa);
      if (x > 0. && kappa + (2. * x * b) >= 0.) {
        roots.push_back(x);
      }
    }
    std::ranges::sort(roots);
  } else if (-qb / (2. * qa) > 0.) {
    roots.push_back(-qb / (2. * qa));
  }
  std::vector<std::array<double, 3>> starts;
  for (const double x : roots) {
    const std::array<double, 3> q{(x * u[0]) - total[0], (x * u[1]) - total[1], (x * u[2]) - total[2]};
    if (std::hypot(q[0], q[1], q[2]) > 0.) {
      starts.push_back(q);
    }
  }
  return starts;
}

namespace decay_fit_detail
{

/// Whether the inputs fit the fit's capacities (`kTooLarge`: more vertices
/// or daughters than the constraints take, a tracked particle with more rows
/// than a track of Det holds) and refer only inside the decay
/// (`kInvalidInput`); checked before anything walks them.
template <DetectorLayout Det>
[[nodiscard]] DecayFitStatus validateInputs(std::span<const DecayDaughter> daughters, const DecayConstraints& constraints, const DecayTopology& topology, std::span<const DecayMother> mothers)
{
  const int count = static_cast<int>(daughters.size());
  const int k = topology.vertices();
  if (count > kMaxConstrainedDaughters || k > kMaxConstrainedVertices) {
    return DecayFitStatus::kTooLarge;
  }
  if (!topology.vertexOf.empty() && static_cast<int>(topology.vertexOf.size()) != count) {
    return DecayFitStatus::kInvalidInput;
  }
  for (const int v : topology.vertexOf) {
    if (v < 0 || v >= k) {
      return DecayFitStatus::kInvalidInput;
    }
  }
  for (int v = 0; v < static_cast<int>(topology.parent.size()); ++v) {
    if (v == 0 ? topology.parent[v] != -1 : (topology.parent[v] < 0 || topology.parent[v] >= v)) {
      return DecayFitStatus::kInvalidInput;
    }
  }
  if (static_cast<int>(topology.massGeV.size()) > k || static_cast<int>(constraints.compositeMasses.size()) > k) {
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
/// one (DecayConstraints.h), taken whole unless a helix misses its plane;
/// a tree (\a topology) has a vertex per composite, each daughter carried to
/// the plane about its own, and the links always constrained. A fitter fits
/// one decay, once (`fit`).
template <DetectorLayout Det, DecayFitVariant Variant, typename Material>
class DecayVertexFitter
{
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
    std::vector<double> loss;          ///< the model: q/p's factor after each row
    TrackingSurface endPlane;          ///< the model: the plane about its end point
  };
  struct Measurement {
    TrackingSurface plane;
    std::array<double, 5> value{};
    std::array<double, 25> weight{}; ///< the inverse covariance
  };
  // An electron's carried mixture, its components kept apart for the
  // Gaussian sum (`kGaussianSum`): ln of its normalised weight, its
  // measurement, ln det of its covariance.
  struct Component {
    double logWeight = 0.;
    Measurement measurement;
    double logDetCovariance = 0.;
  };

 public:
  DecayVertexFitter(std::span<const DecayDaughter> daughters, const Material& material, const DecayConstraints& constraints, const DecayTopology& topology, std::span<const DecayMother> mothers)
    : mDaughters(daughters), mMaterial(material), mConstraints(constraints), mTopology(topology), mMothers(mothers)
  {
  }

  [[nodiscard]] DecayVertexFitResult fit()
  {
    using namespace decay_fit_detail;
    if (const auto status = validateInputs<Det>(mDaughters, mConstraints, mTopology, mMothers); status != DecayFitStatus::kOk) {
      mResult.status = status;
      return mResult;
    }
    mCount = static_cast<int>(mDaughters.size());
    mVertexCount = mTopology.vertices();
    mCore = (3 * mVertexCount) + (3 * mCount);
    mCharges.resize(static_cast<std::size_t>(mCount));
    for (int i = 0; i < mCount; ++i) {
      mCharges[i] = mDaughters[i].charge();
    }
    mResult.charges = mCharges;
    setUpChains();
    mParameterCount = mCore + mChainCount;
    mScattered = ScatteredEquations{.at = mCore, .count = mChainCount};
    if (mCount < 2) {
      return mResult;
    }
    const auto start = treeStart(mDaughters, mTopology);
    if (!start) {
      return mResult;
    }
    mMeasurements.resize(static_cast<std::size_t>(mCount));
    mComponents.resize(static_cast<std::size_t>(mCount));
    mResult.carried.resize(static_cast<std::size_t>(mCount));
    mResult.carriedCovariance.resize(static_cast<std::size_t>(mCount));
    mResult.carryX2X0.resize(static_cast<std::size_t>(mCount));

    // The field a charged mother flies in: a charged daughter's (a neutral
    // pseudo-track's propagator has none).
    for (const auto& daughter : mDaughters) {
      if (!daughter.neutral) {
        mPhysicalField = static_cast<double>(daughter.propagator.bzKiloGauss) * static_cast<double>(constants::math::B2C);
        break;
      }
    }

    // The parameters: the vertices, then each daughter's phi, theta, q/p.
    mParameters.resize(static_cast<std::size_t>(mParameterCount));
    for (int w = 0; w < mVertexCount; ++w) {
      std::copy((*start)[w].begin(), (*start)[w].end(), mParameters.begin() + (3 * w));
    }
    if (const auto status = carry(mParameters); status != DecayFitStatus::kOk) {
      mResult.status = status;
      return mResult;
    }
    for (int i = 0; i < mCount; ++i) {
      for (int c = 0; c < 3; ++c) {
        mParameters[(3 * mVertexCount) + (3 * i) + c] = mResult.carried[i][c];
      }
    }
    // The chain's model, along the start (it needs the momenta).
    if (mChainCount > 0) {
      if (const auto status = refreshChains(mParameters); status != DecayFitStatus::kOk) {
        mResult.status = status;
        return mResult;
      }
    }

    mConstrained = mConstraints.count(mTopology) > 0 || mChainCount > 0;
    mMasses.resize(static_cast<std::size_t>(mCount));
    for (int i = 0; i < mCount; ++i) {
      mMasses[i] = mDaughters[i].mass;
    }
    mNormal.resize(static_cast<std::size_t>(mParameterCount * mParameterCount));
    mGradient.resize(static_cast<std::size_t>(mParameterCount));
    int measured = 0;
    for (const auto& daughter : mDaughters) {
      measured += daughter.unmeasured ? -3 : 2; // five measured against three; none against three
    }
    mNdf = measured - (3 * mVertexCount) + mConstraints.count(mTopology) + (2 * mChainRows);
    if (!modelLoop() || !covariance()) {
      return mResult;
    }
    if constexpr (Variant == DecayFitVariant::kGaussianSum) {
      gaussianSum();
    }
    std::copy_n(mParameters.begin(), 3, mResult.vertex.begin());
    mResult.vertices.resize(static_cast<std::size_t>(mVertexCount));
    for (int w = 0; w < mVertexCount; ++w) {
      std::copy_n(mParameters.begin() + (3 * w), 3, mResult.vertices[w].begin());
    }
    mResult.momenta.resize(static_cast<std::size_t>(mCount));
    for (int i = 0; i < mCount; ++i) {
      for (int c = 0; c < 3; ++c) {
        mResult.momenta[i][c] = mParameters[(3 * mVertexCount) + (3 * i) + c];
      }
    }
    // The vertices' and momenta's block: the nuisances left out.
    if (mParameterCount != mCore) {
      std::vector<double> block(static_cast<std::size_t>(mCore * mCore));
      for (int a = 0; a < mCore; ++a) {
        std::copy_n(mResult.covariance.begin() + (a * mParameterCount), mCore, block.begin() + (a * mCore));
      }
      mResult.covariance = std::move(block);
    }
    mResult.chi2 = mChi2;
    mResult.ndf = mNdf;
    return mResult;
  }

 private:
  static constexpr int kDerivatives = 6;                                     // vx, vy, vz, phi, theta, q/p
  static constexpr double kEpsilon = std::numeric_limits<Scalar>::epsilon(); // the measurements' precision

  [[nodiscard]] double chargeUnder(int w) const
  {
    using namespace decay_fit_detail;
    double sum = 0.;
    for (int i = 0; i < mCount; ++i) {
      sum += mTopology.under(i, w) ? mCharges[i] : 0.;
    }
    return sum;
  }

  [[nodiscard]] const DecayMother* rowsAt(int w) const
  {
    using namespace decay_fit_detail;
    const auto found = std::ranges::find_if(mMothers, [&](const DecayMother& m) { return m.vertex == w; });
    return found == mMothers.end() ? nullptr : &*found;
  }

  [[nodiscard]] double massAt(int w) const
  {
    using namespace decay_fit_detail;
    if (const auto* own = rowsAt(w); own != nullptr) {
      return static_cast<double>(own->propagator.massGeV);
    }
    return static_cast<std::size_t>(w) < mTopology.massGeV.size() ? mTopology.massGeV[w] : (w == 0 ? mConstraints.massGeV : 0.);
  }

  void setUpChains()
  {
    if (mConstraints.production || rowsAt(0) != nullptr) {
      mChains.push_back({.vertex = 0, .rows = rowsAt(0), .end = mConstraints.production ? -1 : -2, .charge = chargeUnder(0), .mass = massAt(0)});
    }
    for (int w = 1; w < mVertexCount; ++w) {
      if (chargeUnder(w) != 0.) {
        mChains.push_back({.vertex = w, .rows = rowsAt(w), .end = mTopology.parent[w], .charge = chargeUnder(w), .mass = massAt(w)});
      }
    }
    for (auto& chain : mChains) {
      const int rows = chain.rows != nullptr ? static_cast<int>(chain.rows->planes.size()) : 0;
      chain.first = mChainCount;
      chain.count = (2 * rows) + (chain.end != -2 ? 2 : 0);
      chain.loss.assign(static_cast<std::size_t>(rows), 1.);
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
  [[nodiscard]] decay_fit_detail::FreeState<decay_fit_detail::kConstraintInputs> chainStart(const std::vector<double>& p, const Chain& chain) const
  {
    using namespace decay_fit_detail;
    using DC = Dual<kConstraintInputs>;
    FreeState<kConstraintInputs> state{};
    for (int c = 0; c < 3; ++c) {
      state[c].v = p[(3 * chain.vertex) + c];
      state[c].d[(3 * chain.vertex) + c] = 1.;
    }
    std::array<DC, 3> total{};
    for (int i = 0; i < mCount; ++i) {
      if (!mTopology.under(i, chain.vertex)) {
        continue;
      }
      const int q = (3 * mVertexCount) + (3 * i);
      const auto momentum = momentumOf(dual::variable<kConstraintInputs>(p[q], q), dual::variable<kConstraintInputs>(p[q + 1], q + 1), dual::variable<kConstraintInputs>(p[q + 2], q + 2), numeratorOf(mCharges[i]));
      for (int c = 0; c < 3; ++c) {
        total[c] = total[c] + momentum.vector[c];
      }
    }
    const DC pt = sqrt((total[0] * total[0]) + (total[1] * total[1]));
    const DC pTotal = sqrt((pt * pt) + (total[2] * total[2]));
    state[3] = atan2(total[1], total[0]);
    state[4] = atan2(pt, total[2]);
    state[5] = DC{chain.charge} / pTotal;
    return state;
  }

  // A chain's end point at \a p: C, or its parent's vertex (a parameter).
  [[nodiscard]] std::array<decay_fit_detail::Dual<decay_fit_detail::kConstraintInputs>, 3> endPoint(const std::vector<double>& p, const Chain& chain) const
  {
    using namespace decay_fit_detail;
    std::array<Dual<kConstraintInputs>, 3> point{};
    for (int c = 0; c < 3; ++c) {
      if (chain.end == -1) {
        point[c].v = mConstraints.point[c];
      } else {
        point[c].v = p[(3 * chain.end) + c];
        point[c].d[(3 * chain.end) + c] = 1.;
      }
    }
    return point;
  }

  // The plane about \a point for a state at (x, y) with direction phi.
  [[nodiscard]] TrackingSurface endPlaneOf(const std::array<double, 3>& point, double x, double y, double phi, double theta, double qOverP) const
  {
    using namespace decay_fit_detail;
    const double curvature = mPhysicalField * qOverP / std::sin(theta);
    const double dx = x - point[0];
    const double dy = y - point[1];
    const double along = (dx * std::cos(phi)) + (dy * std::sin(phi));
    const double across = (-dx * std::sin(phi)) + (dy * std::cos(phi));
    return beamLinePlane(static_cast<Scalar>(point[0]), static_cast<Scalar>(point[1]), static_cast<Scalar>(phi + std::atan2(-curvature * along, 1. + (curvature * across))));
  }

  // The chains' equations at \a p, predicted less measured, their rows of
  // derivatives by (the vertices, the momenta), under the current model.
  [[nodiscard]] bool chainEquations(const std::vector<double>& p, std::vector<double>& equations, std::vector<double>& coreRows) const
  {
    using namespace decay_fit_detail;
    const auto push = [&](const Dual<kConstraintInputs>& value) {
      equations.push_back(value.v);
      for (int c = 0; c < mCore; ++c) {
        coreRows.push_back(value.d[c]);
      }
    };
    for (const auto& chain : mChains) {
      auto state = chainStart(p, chain);
      const int rows = static_cast<int>(chain.loss.size());
      for (int j = 0; j < rows; ++j) {
        FreeState<kConstraintInputs> to{};
        std::array<Dual<kConstraintInputs>, 2> local{};
        if (!crossToPlane(state, chain.rows->planes[j], mPhysicalField, to, local)) {
          return false;
        }
        push(local[0] + (-chain.rows->positions[j][0]));
        push(local[1] + (-chain.rows->positions[j][1]));
        state = to;
        state[5] = chain.loss[j] * state[5];
      }
      if (chain.end != -2) {
        FreeState<kConstraintInputs> to{};
        std::array<Dual<kConstraintInputs>, 2> local{};
        if (!crossToPlane(state, chain.endPlane, mPhysicalField, to, local)) {
          return false;
        }
        // the end point's own place on the plane
        const auto point = endPoint(p, chain);
        const auto& plane = chain.endPlane;
        std::array<Dual<kConstraintInputs>, 2> target{};
        for (int c = 0; c < 3; ++c) {
          const auto d = point[c] + (-static_cast<double>(plane.origin[c]));
          target[0] = target[0] + (static_cast<double>(plane.axis0[c]) * d);
          target[1] = target[1] + (static_cast<double>(plane.axis1[c]) * d);
        }
        push(local[0] - target[0]);
        push(local[1] - target[1]);
      }
    }
    return true;
  }

  // The chains' model along \a p's trajectories.
  [[nodiscard]] DecayFitStatus refreshChains(const std::vector<double>& p)
  {
    using namespace decay_fit_detail;
    std::vector<double> covariance(static_cast<std::size_t>(mChainCount * mChainCount), 0.);
    for (auto& chain : mChains) {
      const auto begin = chainStart(p, chain);
      std::array<double, 6> state{};
      for (int c = 0; c < 6; ++c) {
        state[c] = begin[c].v;
      }
      const int rows = static_cast<int>(chain.loss.size());
      const bool scatters = chain.charge != 0.;
      const int planes = rows + (chain.end != -2 ? 1 : 0);
      std::array<double, 3> end{};
      if (chain.end != -2) {
        const auto point = endPoint(p, chain);
        end = {point[0].v, point[1].v, point[2].v};
      }
      // Each crossing's step and local Jacobians by the state before it, the
      // segments' material and the states at the scatterers.
      std::vector<std::array<double, 36>> step(static_cast<std::size_t>(planes));
      std::vector<std::array<double, 12>> local(static_cast<std::size_t>(planes));
      std::vector<o2::gtrk::detail::StepMaterial> segment(static_cast<std::size_t>(rows));
      std::vector<std::array<double, 6>> at(static_cast<std::size_t>(rows + 1));
      std::vector<double> straggling(static_cast<std::size_t>(rows), 0.);
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
        FreeState<6> from{};
        for (int c = 0; c < 6; ++c) {
          from[c].v = state[c];
          from[c].d[c] = 1.;
        }
        FreeState<6> to{};
        std::array<Dual<6>, 2> loc{};
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
          double scale = 1.;
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
          lastLength = std::hypot(state[0] - end[0], state[1] - end[1], state[2] - end[2]);
          totalX2X0 += static_cast<double>(last.x2X0);
        }
      }
      // This chain's block: its rows' own, then each scatterer's kink carried
      // to every later equation, then the last segment's onto the end.
      const auto entry = [&](int a, int b) -> double& { return covariance[((chain.first + a) * mChainCount) + chain.first + b]; };
      for (int j = 0; j < rows; ++j) {
        const int e = 2 * j;
        entry(e, e) = chain.rows->covariance[j][0];
        entry(e, e + 1) = chain.rows->covariance[j][1];
        entry(e + 1, e) = chain.rows->covariance[j][1];
        entry(e + 1, e + 1) = chain.rows->covariance[j][2];
      }
      if (!scatters) {
        continue;
      }
      const auto theta0 = [&](Scalar x2X0, const std::array<double, 6>& where) {
        return static_cast<double>(highlandTheta0(x2X0, static_cast<Scalar>(totalX2X0), static_cast<Scalar>(std::abs(chain.charge / where[5])), static_cast<Scalar>(chain.mass)));
      };
      // scatterer k sits at at[k]: V (k = 0) or after row k - 1
      for (int kink = 0; kink <= rows; ++kink) {
        const auto& where = at[kink];
        double angle = 0.;
        if (kink < rows) {
          const double nearer = theta0(segment[kink].nearerX2X0, where);
          angle += nearer * nearer;
        }
        if (kink > 0) {
          const double farther = theta0(segment[kink - 1].fartherX2X0, where);
          angle += farther * farther;
        }
        const double sinTheta = std::sin(where[4]);
        const std::array<double, 3> variance{angle / (sinTheta * sinTheta), angle, kink > 0 ? straggling[kink - 1] : 0.};
        // the kink's sensitivity, 6 x 3, carried through the later crossings
        std::array<double, 18> g{};
        for (int c = 0; c < 3; ++c) {
          g[((3 + c) * 3) + c] = 1.;
        }
        std::vector<double> jac(static_cast<std::size_t>(chain.count * 3), 0.);
        for (int j = kink; j < planes; ++j) {
          for (int r = 0; r < 2; ++r) {
            for (int c = 0; c < 3; ++c) {
              double value = 0.;
              for (int q = 0; q < 6; ++q) {
                value += local[j][(r * 6) + q] * g[(q * 3) + c];
              }
              jac[(((2 * j) + r) * 3) + c] = value;
            }
          }
          std::array<double, 18> next{};
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
          for (int b = 0; b < chain.count; ++b) {
            for (int c = 0; c < 3; ++c) {
              entry(a, b) += jac[(a * 3) + c] * variance[c] * jac[(b * 3) + c];
            }
          }
        }
      }
      if (chain.end != -2 && last.x2X0 > 0.F) {
        const double t0 = static_cast<double>(highlandTheta0(last.x2X0, static_cast<Scalar>(totalX2X0), static_cast<Scalar>(std::abs(chain.charge / state[5])), static_cast<Scalar>(chain.mass)));
        const double across = t0 * t0 * lastLength * lastLength * static_cast<double>(last.x2X0T2) / static_cast<double>(last.x2X0);
        const double sinTheta = std::sin(state[4]);
        const int e = 2 * rows;
        entry(e, e) += across;
        entry(e + 1, e + 1) += across / (sinTheta * sinTheta);
      }
    }
    for (const auto& chain : mChains) {
      if (chain.end == -1) {
        const int e = chain.first + chain.count - 2;
        mResult.pointSigma = std::sqrt(covariance[(e * mChainCount) + e]);
      }
    }
    mScattered.factor = semidefiniteFactor(covariance, mChainCount);
    return DecayFitStatus::kOk;
  }

  // Each daughter carried to the perigee plane about its vertex in \a at.
  [[nodiscard]] DecayFitStatus carry(const std::vector<double>& at)
  {
    using namespace decay_fit_detail;
    for (int i = 0; i < mCount; ++i) {
      const auto& daughter = mDaughters[i];
      if (daughter.unmeasured) {
        mResult.carried[i] = {daughter.innermost.phi(), daughter.innermost.theta(), daughter.innermost.qOverP()};
        continue;
      }
      const int w = mTopology.vertexOfDaughter(i);
      const auto plane = detail::hypothesis::perigeePlaneOf(daughter.propagator, daughter.innermost, daughter.surface, static_cast<Scalar>(at[3 * w]), static_cast<Scalar>(at[(3 * w) + 1]));
      BoundTrackState state = daughter.innermost;
      auto factors = factorCovariance(state.cov);
      PropagationStatus status = PropagationStatus::kOk;
      mResult.carryX2X0[i] = detail::hypothesis::segmentRadiation(daughter.propagator, mMaterial, daughter.innermost, daughter.surface, plane);
      if (daughter.mixture) {
        // The electron: its mixture branched by the segment's radiation, every
        // component carried, the collapse on the plane.
        auto mixture = *daughter.mixture;
        if (!(mResult.carryX2X0[i] >= 0.) || !detail::gsf::takeRadiation<BetheHeitlerLoss>(static_cast<Scalar>(mResult.carryX2X0[i]), PropagationDirection::kBackward, mixture)) {
          return DecayFitStatus::kCarryFailed;
        }
        status = detail::gsf::stepEach(mixture, [&](GSFComponent<CovarianceFactors>& component) {
          return detail::hypothesis::carryInward(daughter.propagator, mMaterial, daughter.surface, plane, component.state, component.factors);
        });
        if (status != PropagationStatus::kOk || mixture.size == 0) {
          return DecayFitStatus::kCarryFailed;
        }
        Scalar mode = 0.F;
        detail::gsf::finishMixture(mixture, state, mode);
        auto& kept = mComponents[i];
        kept.clear();
        double largest = -std::numeric_limits<double>::infinity();
        for (int j = 0; j < mixture.size; ++j) {
          largest = std::max(largest, static_cast<double>(mixture.components[j].logWeight));
        }
        double total = 0.;
        for (int j = 0; j < mixture.size; ++j) {
          total += std::exp(static_cast<double>(mixture.components[j].logWeight) - largest);
        }
        for (int j = 0; j < mixture.size; ++j) {
          const auto& one = mixture.components[j].state;
          Cholesky llt;
          if (!llt.factor(covarianceOf(one), 5)) {
            continue; // a component without a covariance is no measurement
          }
          Component component{.logWeight = static_cast<double>(mixture.components[j].logWeight) - largest - std::log(total)};
          component.measurement.plane = plane;
          for (int r = 0; r < 5; ++r) {
            component.measurement.value[r] = one.par[r];
          }
          component.logDetCovariance = llt.logDeterminant();
          weightOf(llt, component.measurement.weight);
          kept.push_back(component);
        }
      } else {
        status = detail::hypothesis::carryInward(daughter.propagator, mMaterial, daughter.surface, plane, state, factors);
        if (status != PropagationStatus::kOk) {
          return DecayFitStatus::kCarryFailed;
        }
        covarianceFromFactors(factors, state.cov);
      }
      Cholesky llt;
      if (!llt.factor(covarianceOf(state), 5)) {
        return DecayFitStatus::kSingular;
      }
      auto& m = mMeasurements[i];
      m.plane = plane;
      for (int r = 0; r < 5; ++r) {
        m.value[r] = state.par[r];
      }
      weightOf(llt, m.weight);
      mResult.carried[i] = {state.phi(), state.theta(), state.qOverP()};
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          mResult.carriedCovariance[i][r * 3 + c] = state.covAt(r + 2, c + 2);
        }
      }
    }
    return DecayFitStatus::kOk;
  }

  // The chi2 at \a p and, with \a normal and \a gradient, the normal
  // equations' matrix and right-hand side; false if a helix misses its plane.
  [[nodiscard]] bool evaluate(const std::vector<double>& p, double& chi2, double& rounding, std::vector<double>* normal, std::vector<double>* gradient) const
  {
    using namespace decay_fit_detail;
    chi2 = 0.;
    rounding = 0.;
    if (normal != nullptr) {
      std::ranges::fill(*normal, 0.);
      std::ranges::fill(*gradient, 0.);
    }
    for (int i = 0; i < mCount; ++i) {
      const auto& daughter = mDaughters[i];
      if (daughter.unmeasured) {
        continue; // no measurement
      }
      const double fieldScale = daughter.neutral ? 0. : static_cast<double>(daughter.propagator.bzKiloGauss) * static_cast<double>(constants::math::B2C);
      using D = Dual<kDerivatives>;
      const int q = (3 * mVertexCount) + (3 * i);
      const int w = 3 * mTopology.vertexOfDaughter(i);
      const std::array<int, kDerivatives> index{w, w + 1, w + 2, q, q + 1, q + 2};
      FreeState<kDerivatives> from{};
      for (int c = 0; c < kDerivatives; ++c) {
        from[c] = dual::variable<kDerivatives>(p[index[c]], c);
      }
      FreeState<kDerivatives> to{};
      std::array<D, 2> local{};
      if (!crossToPlane(from, mMeasurements[i].plane, fieldScale, to, local)) {
        return false;
      }
      const std::array<D, 5> predicted{local[0], local[1], to[3], to[4], to[5]};
      std::array<double, 5> residual{};
      for (int r = 0; r < 5; ++r) {
        residual[r] = mMeasurements[i].value[r] - predicted[r].v;
      }
      residual[2] = std::remainder(residual[2], 2. * std::numbers::pi);
      const auto& weight = mMeasurements[i].weight;
      for (int r = 0; r < 5; ++r) {
        rounding += std::pow(kEpsilon * mMeasurements[i].value[r], 2) * weight[r * 5 + r];
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
          (*normal)[index[a] * mParameterCount + index[b]] += h;
        }
      }
    }
    // The chains' rows: the inputs' precision (`minimiseFrozen`).
    for (const auto& mother : mMothers) {
      for (std::size_t j = 0; j < mother.planes.size(); ++j) {
        rounding += (std::pow(kEpsilon * mother.positions[j][0], 2) / mother.covariance[j][0]) + (std::pow(kEpsilon * mother.positions[j][1], 2) / mother.covariance[j][2]);
      }
    }
    mScattered.addPrior(p, mParameterCount, chi2, normal, gradient);
    // The resonances' masses (`Resonance`): Breit-Wigner terms.
    for (const auto& resonance : mConstraints.resonances) {
      using DR = Dual<kConstraintInputs>;
      DR energy{};
      std::array<DR, 3> total{};
      for (const int i : resonance.daughters) {
        const int q = (3 * mVertexCount) + (3 * i);
        const auto momentum = momentumOf(dual::variable<kConstraintInputs>(p[q], q), dual::variable<kConstraintInputs>(p[q + 1], q + 1), dual::variable<kConstraintInputs>(p[q + 2], q + 2), numeratorOf(mCharges[i]));
        for (int c = 0; c < 3; ++c) {
          total[c] = total[c] + momentum.vector[c];
        }
        energy = energy + sqrt((momentum.magnitude * momentum.magnitude) + (mDaughters[i].mass * mDaughters[i].mass));
      }
      const DR mass2 = (energy * energy) - ((total[0] * total[0]) + (total[1] * total[1]) + (total[2] * total[2]));
      if (!(mass2.v > 0.)) {
        return false;
      }
      const double big = resonance.massGeV;
      const DR u = (mass2 + (-(big * big))) / DR{big * resonance.widthGeV};
      const double u2 = u.v * u.v;
      chi2 += 2. * std::log1p(u2);
      const double weight = 0.5;                 // the Fisher information of chi2 / 2 per unit u^2
      const double slope = 2. * u.v / (1. + u2); // d(chi2 / 2) / du
      if (normal == nullptr) {
        continue;
      }
      for (int a = 0; a < mCore; ++a) {
        (*gradient)[a] -= slope * u.d[a];
        for (int b = 0; b < mCore; ++b) {
          (*normal)[a * mParameterCount + b] += weight * u.d[a] * u.d[b];
        }
      }
    }
    return true;
  }

  // The constraints' values and rows at \a p.
  [[nodiscard]] bool constraintsAt(const std::vector<double>& p)
  {
    using namespace decay_fit_detail;
    std::vector<std::array<double, 3>> vertices(static_cast<std::size_t>(mVertexCount));
    for (int w = 0; w < mVertexCount; ++w) {
      vertices[w] = {p[3 * w], p[(3 * w) + 1], p[(3 * w) + 2]};
    }
    std::vector<std::array<double, 3>> momenta(static_cast<std::size_t>(mCount));
    for (int i = 0; i < mCount; ++i) {
      const int q = (3 * mVertexCount) + (3 * i);
      momenta[i] = {p[q], p[q + 1], p[q + 2]};
    }
    std::vector<double> coreRows;
    if (!evaluateConstraints(mConstraints, mTopology, vertices, momenta, mMasses, mCharges, mPhysicalField, mValues, coreRows)) {
      return false;
    }
    mScattered.row = static_cast<int>(mValues.size());
    if (mChainCount > 0 && !chainEquations(p, mValues, coreRows)) {
      return false;
    }
    const int m = static_cast<int>(mValues.size());
    mRows.assign(static_cast<std::size_t>(m * mParameterCount), 0.);
    for (int a = 0; a < m; ++a) {
      std::copy_n(coreRows.begin() + (a * mCore), mCore, mRows.begin() + (a * mParameterCount));
    }
    if (mChainCount > 0) {
      mScattered.apply(p, mParameterCount, mValues, mRows);
    }
    return true;
  }

  /// `minimiseFrozen` under the current model, from \a parameters.
  [[nodiscard]] DecayFitStatus minimise(std::vector<double>& parameters, double& chi2, double& mu, double& metric, double& bound, int& iterations, int& steps, double& first)
  {
    const auto evaluateAt = [this](const std::vector<double>& p, double& value, double& rounding, std::vector<double>* normal, std::vector<double>* gradient) { return evaluate(p, value, rounding, normal, gradient); };
    const auto constraintsOf = [this](const std::vector<double>& p) { return constraintsAt(p); };
    return decay_fit_detail::minimiseFrozen(mParameterCount, mNdf, mConstrained, evaluateAt, constraintsOf, mValues, mRows, parameters, chi2, mu, metric, bound, iterations, steps, first);
  }

  // The carried measurements are the model: minimised under them, carried
  // anew to the planes about the solution, until a carry no longer moves the
  // solution by what the inputs resolve (the minimiser takes no step).
  /// False if the fit ends here (`mResult` is final).
  [[nodiscard]] bool modelLoop()
  {
    using namespace decay_fit_detail;
    double lastShift = std::numeric_limits<double>::infinity();
    mResult.status = DecayFitStatus::kNotConverged;
    for (int outer = 0; outer < kMaxSteps; ++outer) {
      if (outer != 0) {
        if (const auto status = carry(mParameters); status != DecayFitStatus::kOk) {
          mResult.status = status;
          return false;
        }
        if (mChainCount > 0) {
          if (const auto status = refreshChains(mParameters); status != DecayFitStatus::kOk) {
            mResult.status = status;
            return false;
          }
        }
      }
      mResult.outer = static_cast<std::uint8_t>(outer + 1);
      const auto inner = minimise(mParameters, mChi2, mMu, mMetric, mBound, mIterations, mSteps, mShift);
      if (inner != DecayFitStatus::kOk) {
        mResult.status = inner;
        break;
      }
      mResult.modelShift = static_cast<float>(mShift);
      if (mSteps == 0) {
        // The model anew along the solution moves it less than the inputs
        // resolve: converged.
        mResult.status = DecayFitStatus::kOk;
        mResult.stopReason = 1U;
        break;
      }
      if (outer > 1 && mShift >= lastShift) {
        // The model's moves no longer shrink: the material table's
        // discretisation (piecewise constant along the trajectory) leaves no
        // exact fixed point. Each solution here is a converged minimum under a
        // model along a trajectory close to it; taken if the model moves it by
        // less than its error (`kModelTolerance`).
        mResult.status = mShift <= kModelTolerance ? DecayFitStatus::kOk : DecayFitStatus::kNotConverged;
        mResult.stopReason = mShift <= kModelTolerance ? 2U : 0U;
        break;
      }
      if (outer != 0) {
        lastShift = mShift;
      }
    }
    mResult.iterations = static_cast<std::uint8_t>(std::min(mIterations, 255));
    mResult.stopMetric = static_cast<float>(mMetric);
    mResult.stopBound = static_cast<float>(mBound);
    return true;
  }

  // The covariance at the solution, with the measurements of the last carry.
  /// False if the fit ends here (`mResult` is final).
  [[nodiscard]] bool covariance()
  {
    using namespace decay_fit_detail;
    double rounding = 0.;
    if (!evaluate(mParameters, mChi2, rounding, &mNormal, &mGradient)) {
      mResult.status = DecayFitStatus::kNoCrossing;
      return false;
    }
    if (mConstrained) {
      std::vector<double> unused;
      if (!constraintsAt(mParameters) || !constrainedStep(mNormal, mParameterCount, mGradient, mRows, mValues, unused, &mResult.covariance)) {
        mResult.status = DecayFitStatus::kSingular;
        return false;
      }
    } else {
      Cholesky llt;
      if (!llt.factor(mNormal, mParameterCount)) {
        mResult.status = DecayFitStatus::kSingular;
        return false;
      }
      mResult.covariance = llt.inverse();
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
  void gaussianSum()
  {
    using namespace decay_fit_detail;
    std::vector<int> electrons;
    for (int i = 0; i < mCount; ++i) {
      if (!mComponents[i].empty()) {
        electrons.push_back(i);
      }
    }
    if (!electrons.empty()) {
      const auto collapsed = mMeasurements;
      const auto from = mParameters;
      // Combination k picks, per electron, its component k's digit in the
      // mixed radix of the mixtures' sizes, the first electron's the fastest.
      std::size_t tried = 1U;
      for (const int e : electrons) {
        tried *= mComponents[e].size();
      }
      MixtureCollapse collapse(mParameterCount);
      std::vector<double> fitParameters(static_cast<std::size_t>(mParameterCount));
      std::vector<double> fitCovariance;
      for (std::size_t combination = 0; combination < tried; ++combination) {
        double logWeight = 0.;
        std::size_t digits = combination;
        for (const int e : electrons) {
          const auto& component = mComponents[e][digits % mComponents[e].size()];
          digits /= mComponents[e].size();
          mMeasurements[e] = component.measurement;
          logWeight += component.logWeight - (0.5 * component.logDetCovariance);
        }
        fitParameters = from;
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
        std::vector<double> matrix(static_cast<std::size_t>(mParameterCount * mParameterCount));
        std::vector<double> rhs(static_cast<std::size_t>(mParameterCount));
        bool ok = evaluate(fitParameters, fitChi2, unusedRounding, &matrix, &rhs);
        double logVolume = 0.;
        if (ok && mConstrained) {
          std::vector<double> unused;
          double logKKT = 0.;
          ok = constraintsAt(fitParameters) && constrainedStep(matrix, mParameterCount, rhs, mRows, mValues, unused, &fitCovariance, nullptr, &logKKT);
          // C C^T, m x m.
          const int m = static_cast<int>(mValues.size());
          std::vector<double> square(static_cast<std::size_t>(m * m), 0.);
          for (int j = 0; j < m; ++j) {
            for (int c = 0; c <= j; ++c) {
              double sum = 0.;
              for (int b = 0; b < mParameterCount; ++b) {
                sum += mRows[(j * mParameterCount) + b] * mRows[(c * mParameterCount) + b];
              }
              square[(j * m) + c] = sum;
              square[(c * m) + j] = sum;
            }
          }
          Cholesky llt;
          ok = ok && llt.factor(square, m);
          logVolume = ok ? logKKT - llt.logDeterminant() : 0.;
        } else if (ok) {
          Cholesky llt;
          ok = llt.factor(matrix, mParameterCount);
          if (ok) {
            logVolume = llt.logDeterminant();
            fitCovariance = llt.inverse();
          }
        }
        if (ok) {
          collapse.add(logWeight - (0.5 * fitChi2) - (0.5 * logVolume), fitParameters, fitCovariance, fitChi2);
        }
      }
      mMeasurements = collapsed;
      mResult.combinationsTried = static_cast<std::uint16_t>(std::min<std::size_t>(tried, 65535U));
      if (collapse.count == 0) {
        mResult.status = DecayFitStatus::kNotConverged;
      } else {
        mParameters = collapse.mean;
        mChi2 = collapse.meanChi2();
        mResult.covariance = collapse.covariance();
        mResult.combinations = static_cast<std::uint16_t>(std::min(collapse.count, 65535));
      }
    }
  }

  std::span<const DecayDaughter> mDaughters;
  const Material& mMaterial;
  const DecayConstraints& mConstraints;
  const DecayTopology& mTopology;
  std::span<const DecayMother> mMothers;
  DecayVertexFitResult mResult;
  int mCount = 0;          ///< the daughters
  int mVertexCount = 0;    ///< the tree's vertices
  int mCore = 0;           ///< the vertices' and momenta's parameters, 3K + 3N
  int mParameterCount = 0; ///< with the chains' nuisances after them
  int mChainCount = 0;     ///< the chains' equations
  int mChainRows = 0;      ///< the chains' rows
  std::vector<double> mCharges;
  std::vector<Chain> mChains;
  decay_fit_detail::ScatteredEquations mScattered;
  std::vector<Measurement> mMeasurements;
  std::vector<std::vector<Component>> mComponents;
  double mPhysicalField = 0.;
  std::vector<double> mParameters;
  std::vector<double> mMasses;
  std::vector<double> mValues; ///< the constraints' values
  std::vector<double> mRows;   ///< and their rows of derivatives
  bool mConstrained = false;
  std::vector<double> mNormal;
  std::vector<double> mGradient;
  double mChi2 = 0.;
  int mNdf = 0;
  // the minimiser's state, across the model's evaluations
  double mMu = 0.;
  double mMetric = 0.;
  double mBound = 0.;
  int mIterations = 0;
  int mSteps = 0;
  double mShift = 0.;
};

/// A `DecayVertexFitter`'s fit of one decay.
template <DetectorLayout Det, DecayFitVariant Variant, typename Material>
[[nodiscard]] DecayVertexFitResult fitDecayVertex(std::span<const DecayDaughter> daughters, const Material& material, const DecayConstraints& constraints = {}, const DecayTopology& topology = {}, std::span<const DecayMother> mothers = {})
{
  return DecayVertexFitter<Det, Variant, Material>(daughters, material, constraints, topology, mothers).fit();
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYVERTEXFIT_H_
