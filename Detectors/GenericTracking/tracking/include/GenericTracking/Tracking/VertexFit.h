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

/// \file VertexFit.h
/// \brief What every vertex fit of tracks shares: a track's measurement at
///        the vertex, its term in the normal equations, the minimiser under a
///        fixed model and the loop that renews the model.
///
/// A vertex fit's parameters are its vertices and, per track, its direction
/// and q/p at its vertex (phi, theta, q/p; Billoir's choice). A track enters
/// as its fitted state carried with its species' mass and the material
/// between to the perigee plane about the current vertex; its five
/// parameters there and their weight are the measurement
/// (`TrackMeasurement`). The prediction is the helix from the vertex with the
/// track's parameters to that plane, exact in a uniform field, its
/// derivatives by forward differentiation (`trackTerm`).
///
/// The minimiser (`minimiseFrozen`) takes the linear algebra of the normal
/// equations as a type: the decay fit's is dense (any tree, constraints by
/// their Lagrange system, DecayVertexFit.h), the primary vertex fit's is
/// block-arrow -- each track's momentum couples only to the vertex --, solved
/// by eliminating each track's block (PrimaryVertexFit.h). Both run the same
/// steps, damping and stop, report the same diagnostics
/// (`VertexFitDiagnostics`) and carve their workspaces alike (`place`).

#ifndef GTRK_TRACKING_VERTEXFIT_H_
#define GTRK_TRACKING_VERTEXFIT_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/DecayFitMath.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/HypothesisRefit.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include "GPUCommonDef.h"
#include "GPUCommonMath.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <span>
#include <vector>

namespace o2::gtrk
{

/// A point known to its covariance: the beam line as the primary vertex
/// fit's prior, a primary vertex as a decay's reference point.
struct MeasuredPoint {
  Vec3D position{};                   ///< [cm]
  std::array<double, 6> covariance{}; ///< x, y, z; packed (`covIndex`) [cm^2]
};

/// How a vertex fit ended: the first part of its result (`modelLoop`).
struct VertexFitDiagnostics {
  VertexFitStatus status = VertexFitStatus::kNoStart;
  std::uint8_t iterations = 0U;
  std::uint8_t outer = 0U;      ///< the model's evaluations (carries, material along the trajectory)
  std::uint8_t stopReason = 0U; ///< 1 the model anew moved nothing, 2 its moves stopped shrinking within `kModelTolerance`, 0 not converged
  float stopMetric = 0.F;       ///< the last step's predicted chi2 gain
  float stopBound = 0.F;        ///< and what it was held to: the inputs' resolution R (`minimiseFrozen`)
  float modelShift = 0.F;       ///< what the last model evaluation moved the solution by (its first step's predicted gain)
};

/// A fit's workspace on the host: aligned to
/// `vertex_fit_detail::kWorkspaceAlignment`, grown
/// to the largest asked for and reused.
class WorkspaceBuffer
{
 public:
  /// \a bytes of it.
  [[nodiscard]] std::span<std::byte> take(std::size_t bytes)
  {
    const std::size_t words = (bytes + sizeof(std::max_align_t) - 1U) / sizeof(std::max_align_t);
    if (mWords.size() < words) {
      mWords.resize(words);
    }
    return std::as_writable_bytes(std::span{mWords}).first(bytes);
  }

 private:
  std::vector<std::max_align_t> mWords;
};

} // namespace o2::gtrk

namespace o2::gtrk::vertex_fit_detail
{

/// The alignment a fit's workspace starts at: every array's in it.
inline constexpr std::size_t kWorkspaceAlignment = alignof(std::max_align_t);

/// An array of \a count Ts at \a at (aligned for T) of a workspace's
/// \a base, \a at moved past it; none without \a base (the size alone).
template <typename T>
[[nodiscard]] GPUhdi() T* place(std::byte* base, std::size_t& at, std::size_t count)
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

/// Whether \a workspace holds \a bytes, aligned to kWorkspaceAlignment.
[[nodiscard]] GPUhdi() bool holdsWorkspace(std::span<std::byte> workspace, std::size_t bytes)
{
  return workspace.size() >= bytes && reinterpret_cast<std::uintptr_t>(workspace.data()) % kWorkspaceAlignment == 0U; // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast): the address's alignment
}

/// The measurements' precision: their states are floats.
inline constexpr double kMeasurementPrecision = std::numeric_limits<Scalar>::epsilon();

/// An input \a value's rounding squared, (epsilon_float value)^2: times
/// its weight, its share of the inputs' resolution R (`minimiseFrozen`).
[[nodiscard]] GPUhdi() double roundingOf(double value) { return (kMeasurementPrecision * value) * (kMeasurementPrecision * value); }

/// The helix's curvature per q/pT in \a propagator's field, Bz B2C.
[[nodiscard]] GPUhdi() double fieldScaleOf(const Propagator& propagator) { return static_cast<double>(propagator.bzKiloGauss) * static_cast<double>(constants::math::B2C); }

/// A track's parameters in a vertex fit: its vertex's x, y, z, then its phi,
/// theta, q/p there.
inline constexpr int kTrackParameters = 6;

/// A track's measurement: its state carried to \a plane, the perigee plane
/// about its vertex; the inverse covariance packed.
template <typename S>
struct TrackMeasurement {
  TrackingSurface plane;
  std::array<double, 5> value{};
  std::array<S, kBoundCovSize> weight{};
};

/// \a state on \a plane as a measurement into \a m, and ln det of its
/// covariance into \a logDet if given; false if the covariance is not
/// positive definite.
template <typename S>
[[nodiscard]] GPUd() bool measurementOf(const BoundTrackState& state, const TrackingSurface& plane, TrackMeasurement<S>& m, double* logDet = nullptr)
{
  using namespace decay_fit_detail;
  std::array<S, kBoundCovSize> factor{};
  Cholesky<S> llt{.l = factor.data()};
  if (!llt.factor(covarianceOf<S>(state).data(), 5)) {
    return false;
  }
  m.plane = plane;
  for (int r = 0; r < 5; ++r) {
    m.value[r] = state.par[r];
  }
  if (logDet != nullptr) {
    *logDet = static_cast<double>(llt.logDeterminant());
  }
  llt.inverse(m.weight.data());
  return true;
}

/// \a state carried from \a surface to \a plane with \a propagator's mass
/// and \a material's material between (`carryStateInward`), its measurement
/// there into \a m; kCarryFailed or kSingular where either fails.
template <typename Material, typename S>
[[nodiscard]] GPUd() VertexFitStatus carriedMeasurement(const Propagator& propagator, const Material& material, const TrackingSurface& surface, const TrackingSurface& plane, BoundTrackState& state, TrackMeasurement<S>& m)
{
  if (detail::hypothesis::carryStateInward(propagator, material, surface, plane, state) != PropagationStatus::kOk) {
    return VertexFitStatus::kCarryFailed;
  }
  return measurementOf(state, plane, m) ? VertexFitStatus::kOk : VertexFitStatus::kSingular;
}

/// A track's term in the normal equations over its six parameters
/// (`kTrackParameters`): the lower triangle of H^T W H, packed, and H^T W r.
template <typename S>
struct TrackTerm {
  std::array<S, covSize(kTrackParameters)> normal{};
  std::array<S, kTrackParameters> gradient{};
};

/// The track of measurement \a m at \a at (its six parameters): its chi2
/// added to \a chi2, the inputs' rounding in its metric to \a rounding
/// (`minimiseFrozen`), and with \a term its term; \a fieldScale Bz B2C, 0 for
/// a neutral track's line. False if the helix misses the plane.
template <typename M, typename S>
[[nodiscard]] GPUd() bool trackTerm(const TrackMeasurement<S>& m, double fieldScale, const std::array<double, kTrackParameters>& at, double& chi2, double& rounding, TrackTerm<S>* term)
{
  using namespace decay_fit_detail;
  using D = Dual<kTrackParameters, M>;
  FreeState<kTrackParameters, M> from{};
  for (int c = 0; c < kTrackParameters; ++c) {
    from[c] = dual::variable<kTrackParameters, M>(at[c], c);
  }
  FreeState<kTrackParameters, M> to{};
  std::array<D, 2> local{};
  if (!crossToPlane(from, m.plane, fieldScale, to, local)) {
    return false;
  }
  const std::array<D, 5> predicted{local[0], local[1], to[3], to[4], to[5]};
  std::array<M, 5> residual{};
  for (int r = 0; r < 5; ++r) {
    residual[r] = static_cast<M>(m.value[r]) - predicted[r].v;
  }
  residual[2] = toPMPiOf(residual[2]);
  const auto& weight = m.weight;
  for (int r = 0; r < 5; ++r) {
    rounding += roundingOf(m.value[r]) * static_cast<double>(weight[covIndex(r, r)]);
  }
  // W r in double: the chi2 is a sum in double
  std::array<double, 5> weighted{};
  for (int r = 0; r < 5; ++r) {
    for (int c = 0; c < 5; ++c) {
      weighted[r] += static_cast<double>(weight[covIndex(r, c)]) * static_cast<double>(residual[c]);
    }
    chi2 += static_cast<double>(residual[r]) * weighted[r];
  }
  if (term == nullptr) {
    return true;
  }
  for (int a = 0; a < kTrackParameters; ++a) {
    S g = 0.;
    for (int r = 0; r < 5; ++r) {
      g += static_cast<S>(predicted[r].d[a]) * static_cast<S>(weighted[r]);
    }
    term->gradient[a] = g;
    for (int b = 0; b <= a; ++b) {
      S h = 0.;
      for (int r = 0; r < 5; ++r) {
        for (int c = 0; c < 5; ++c) {
          h += static_cast<S>(predicted[r].d[a]) * weight[covIndex(r, c)] * static_cast<S>(predicted[c].d[b]);
        }
      }
      term->normal[covIndex(a, b)] = h;
    }
  }
  return true;
}

/// The minimiser's state across a fit's model evaluations
/// (`minimiseFrozen`, `modelLoop`).
struct MinimiserState {
  double mu = 0.;     ///< the largest |lambda| of the Gauss-Newton steps seen: the merit's penalty
  double metric = 0.; ///< the last step's predicted gain
  double bound = 0.;  ///< and what it was held to
  int iterations = 0; ///< all minimisations'
  int steps = 0;      ///< the last minimisation's steps taken
  double shift = 0.;  ///< the last minimisation's first predicted gain: what the model moved the solution by
};

/// The Gauss-Newton minimum under a fixed model -- the carried measurements,
/// the material along a trajectory --: smooth; the model evaluated in M, the
/// solves in S, the parameters, chi2, merit and the stopping test in double.
///
/// Converged when the Gauss-Newton step's predicted gain, step^T N step +
/// mu |c|_1 - lambda.c (the merit's halved first-order decrease), is no larger
/// than what the inputs resolve: they are floats, each measured value m known
/// to epsilon_float |m|, which moves the solution by R = sum (epsilon_float
/// |m|)^2 w in the same metric (the evaluation's `rounding`); plus the model
/// arithmetic's rounding of the chi2, epsilon_M (chi2 + ndf). The merit is the
/// chi2 plus 2 mu |c|_1, mu the largest |lambda| of the Gauss-Newton steps seen
/// (the exact penalty: the constrained step is a descent direction of it for
/// mu >= |lambda|_inf); without constraints the chi2.
///
/// The step taken is Levenberg-Marquardt's: (N + tau D) d = g, D = diag N,
/// under hard constraints the Lagrange system with N + tau D. Where N is
/// nearly singular (two nearly parallel tracks: their vertex along their
/// line) the Gauss-Newton step overshoots, as the curvature there is the
/// residuals' second-order terms it leaves out; the damping tau stands in for
/// them. tau starts at 0 (Gauss-Newton) and follows Nielsen's update from the
/// full step's gain ratio rho = (merit - merit(full)) / predicted: accepted,
/// tau x max(1/3, 1 - (2 rho - 1)^3) and nu = 2; rejected, tau x nu (at least
/// the curvature the model missed along the step, (merit(full) - model(full))
/// / d^T D d, which also gives tau's first value), nu doubled. A damped step
/// that is no descent of the merit by its predicted gain (the constraints'
/// restoring part costs more chi2 than the penalty gains) is replaced by the
/// Gauss-Newton step, tau back to 0. Each step is then halved until it lowers
/// the merit; a step that no halving makes lower it above the resolution is a
/// failure, not a minimum.
///
/// \a system, the normal equations' linear algebra, over n parameters:
/// - `bool evaluate(p, chi2, rounding, bool normal)`: the chi2 at p, with
///   `normal` the normal equations there; false if a helix misses its plane;
/// - `bool constrained()`;
/// - `VertexFitStatus newton(p, S* step, double& mu, double& violation,
///   double& lambdaC)`: the Gauss-Newton step from the last evaluation, at p
///   (its constraints' values there, mu raised to its multipliers' largest,
///   |c|_1 and lambda.c);
/// - `bool damped(double tau, S* d)`: the step of N + tau D;
/// - `double quadratic(const S* v, double start)`: start + v^T N v;
/// - `void forms(const S* v, double& gd, double& dNd, double& dDd)`;
/// - `bool trial(p, chi2, violation)`: the chi2 and |c|_1 at a trial point.
/// \a state: mu, the iterations, and this minimisation's steps, its last
/// step's gain and bound and its first step's gain (`MinimiserState`).
/// \a trial, \a step and \a dampedStep are n each.
template <typename M, typename S, typename System>
[[nodiscard]] GPUd() VertexFitStatus minimiseFrozen(System& system, int n, int ndf, std::span<double> parameters, double* trial, S* step, S* dampedStep, double& chi2, MinimiserState& state)
{
  double& mu = state.mu;
  double& metric = state.metric;
  int& steps = state.steps;
  const bool constrained = system.constrained();
  const std::span<double> trialPoint{trial, static_cast<std::size_t>(n)};
  double damping = 0.;
  double growth = 2.;
  steps = 0;
  for (int iteration = 0; iteration < decay_fit_detail::kMaxSteps; ++iteration) {
    ++state.iterations;
    double rounding = 0.;
    if (!system.evaluate(parameters, chi2, rounding, true)) {
      return VertexFitStatus::kNoCrossing;
    }
    double violation = 0.;
    double lambdaC = 0.;
    if (const auto status = system.newton(parameters, step, mu, violation, lambdaC); status != VertexFitStatus::kOk) {
      return status;
    }
    metric = system.quadratic(step, (mu * violation) - lambdaC);
    if (iteration == 0) {
      state.shift = metric;
    }
    const double resolvable = (static_cast<double>(std::numeric_limits<M>::epsilon()) * (chi2 + o2::gpu::CAMath::Max(ndf, 0))) + rounding;
    state.bound = resolvable;
    if (metric <= resolvable) {
      return VertexFitStatus::kOk;
    }
    const double merit = chi2 + (2. * mu * violation);
    // The step d, damped while tau > 0; its quadratic forms and predicted gain.
    S* d = step;
    if (damping > 0.) {
      d = dampedStep;
      if (!system.damped(damping, d)) {
        return VertexFitStatus::kSingular;
      }
    }
    double gd = 0.;
    double dNd = 0.;
    double dDd = 0.;
    const auto formsOf = [&](const S* v) {
      system.forms(v, gd, dNd, dDd);
      return (2. * gd) - dNd + (constrained ? 2. * mu * violation : 0.);
    };
    double predicted = formsOf(d);
    if (damping > 0. && !(predicted > 0.)) {
      d = step;
      damping = 0.;
      predicted = formsOf(d);
    }
    bool reached = false;
    bool taken = false;
    double scale = 1.;
    for (int halving = 0; halving < decay_fit_detail::kMaxSteps && !taken; ++halving, scale *= 0.5) {
      for (int j = 0; j < n; ++j) {
        trialPoint[j] = parameters[j] + (scale * static_cast<double>(d[j]));
      }
      double trialChi2 = 0.;
      double trialViolation = 0.;
      const bool evaluated = system.trial(trialPoint, trialChi2, trialViolation);
      const double trialMerit = evaluated ? trialChi2 + (constrained ? 2. * mu * trialViolation : 0.) : 0.;
      if (halving == 0) {
        // Nielsen's update from the full step.
        if (evaluated && trialMerit < merit) {
          const double t = (2. * (merit - trialMerit) / predicted) - 1.;
          damping *= o2::gpu::CAMath::Max(1. / 3., 1. - (t * t * t));
          growth = 2.;
        } else {
          const double missed = !(dDd > 0.) ? 0. : (evaluated ? (trialMerit - (merit - predicted)) / dDd : dNd / dDd);
          damping = o2::gpu::CAMath::Max(damping * growth, missed);
          growth *= 2.;
        }
      }
      if (!evaluated) {
        continue;
      }
      reached = true;
      taken = trialMerit < merit;
    }
    if (!taken) {
      return reached ? VertexFitStatus::kNotConverged : VertexFitStatus::kNoCrossing;
    }
    decay_fit_detail::copyN(trialPoint.data(), parameters.size(), parameters.data());
    ++steps;
  }
  return VertexFitStatus::kNotConverged;
}

/// The carried measurements are the model: minimised under them
/// (\a minimise(), which updates \a state), the model renewed along the
/// solution (\a refresh()), until a renewal no longer moves the solution by
/// what the inputs resolve (the minimiser takes no step). How it ended into
/// \a out. False if a renewal failed (its status in \a out): the fit ends
/// there; true otherwise, \a out's status the loop's, with the minimiser's
/// iterations and last bound.
template <typename Refresh, typename Minimise>
[[nodiscard]] GPUd() bool modelLoop(const Refresh& refresh, const Minimise& minimise, const MinimiserState& state, VertexFitDiagnostics& out)
{
  double lastShift = std::numeric_limits<double>::infinity();
  out.status = VertexFitStatus::kNotConverged;
  for (int outer = 0; outer < decay_fit_detail::kMaxSteps; ++outer) {
    if (outer != 0) {
      if (const auto status = refresh(); status != VertexFitStatus::kOk) {
        out.status = status;
        return false;
      }
    }
    out.outer = static_cast<std::uint8_t>(outer + 1);
    const auto inner = minimise();
    if (inner != VertexFitStatus::kOk) {
      out.status = inner;
      break;
    }
    out.modelShift = static_cast<float>(state.shift);
    if (state.steps == 0) {
      // The model anew along the solution moves it less than the inputs
      // resolve: converged.
      out.status = VertexFitStatus::kOk;
      out.stopReason = 1U;
      break;
    }
    if (outer > 1 && state.shift >= lastShift) {
      // The model's moves no longer shrink: the material table's
      // discretisation (piecewise constant along the trajectory) leaves no
      // exact fixed point. Each solution here is a converged minimum under a
      // model along a trajectory close to it; taken if the model moves it by
      // less than its error (`kModelTolerance`).
      out.status = state.shift <= decay_fit_detail::kModelTolerance ? VertexFitStatus::kOk : VertexFitStatus::kNotConverged;
      out.stopReason = state.shift <= decay_fit_detail::kModelTolerance ? 2U : 0U;
      break;
    }
    if (outer != 0) {
      lastShift = state.shift;
    }
  }
  out.iterations = static_cast<std::uint8_t>(o2::gpu::CAMath::Min(state.iterations, 255));
  out.stopMetric = static_cast<float>(state.metric);
  out.stopBound = static_cast<float>(state.bound);
  return true;
}

} // namespace o2::gtrk::vertex_fit_detail

#endif // GTRK_TRACKING_VERTEXFIT_H_
