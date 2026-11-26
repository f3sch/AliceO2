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

/// \file PrimaryVertexFit.h
/// \brief A primary vertex fitted from its tracks and the beam line: the
///        vertex and every track's momentum at it.
///
/// The decay fit's mathematics (VertexFit.h) for one vertex of any number of
/// tracks: the parameters are the vertex V and, per track, its phi, theta,
/// q/p at V; each track's fitted state is carried with its species' mass and
/// the material to the perigee plane about the current vertex, the exact
/// helix from V its prediction; the beam line is a Gaussian prior on V
/// (its position and widths, chi2 (V - b)^T B^-1 (V - b)). Gauss-Newton with
/// Levenberg-Marquardt damping under a fixed model, the model renewed along
/// the solution, as the decay fit.
///
/// No constraint couples two tracks, so the normal matrix is block-arrow:
/// each track's momentum block C_i couples only to the vertex (E_i). Each
/// step eliminates the tracks' blocks (Billoir's fit): with A and a the
/// vertex's block and right-hand side, c_i the track's,
///
///     S = A - sum E_i^T C_i^-1 E_i,  h = a - sum E_i^T C_i^-1 c_i,
///     S dV = h,  dq_i = C_i^-1 (c_i - E_i dV),
///
/// the same step as one dense solve, at O(N) per iteration; the damping is
/// added to the blocks' diagonals before the elimination. The vertex's
/// covariance is S^-1 at the solution. A track's own share of S and h at
/// the solution (S_i, h_i) gives the vertex without that track at the same
/// linearisation: V + (S - S_i)^-1 (h - h_i), covariance (S - S_i)^-1.
///
/// The fit is in double: the states come in as floats. It allocates nothing:
/// everything that grows with the tracks lives in the caller's workspace.

#ifndef GTRK_TRACKING_PRIMARYVERTEXFIT_H_
#define GTRK_TRACKING_PRIMARYVERTEXFIT_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/HypothesisRefit.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/TrackingSurface.h"
#include "GenericTracking/Tracking/VertexFit.h"

#include "GPUCommonDef.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>

namespace o2::gtrk
{

/// A track as the primary vertex fit takes it.
struct VertexTrack {
  BoundTrackState innermost; ///< fitted, on `surface` (its species' refit)
  TrackingSurface surface;   ///< its innermost measurement's plane
  Propagator propagator;     ///< the field, with the track's species' mass
};

/// A primary vertex fit's result: how it ended, then the vertex.
struct PrimaryVertexFitResult : VertexFitDiagnostics {
  Vec3D vertex{};
  std::array<double, 6> covariance{}; ///< the vertex's, packed
  double chi2 = 0.;                   ///< the tracks' and the prior's
  double priorChi2 = 0.;              ///< the prior's share
  std::int32_t ndf = 0;               ///< 2 per track: 5 measured, 3 fitted; the prior's 3 against the vertex's 3
};

/// The vertex of \a tracks (at least one) and \a prior (the beam line),
/// with \a material's material between each track's innermost measurement
/// and the vertex. A fitter fits one vertex, once (`fit`); after it, each track's momentum at
/// the vertex, its chi2 and the vertex without it.
template <typename Material>
class PrimaryVertexFitter
{
  using M = double;
  using S = double;
  using Measurement = vertex_fit_detail::TrackMeasurement<S>;
  using Term = vertex_fit_detail::TrackTerm<S>;
  /// A track's block eliminated: the factor of its (damped) C_i, C_i^-1 E_i
  /// (row-major, momentum by vertex) and C_i^-1 c_i.
  struct Eliminated {
    std::array<S, 6> factor{};
    std::array<S, 9> inverseE{};
    std::array<S, 3> inverseC{};
  };
  /// Where everything that grows with the tracks lives in the workspace.
  struct Space {
    double* parameters = nullptr; ///< n = 3 + 3N
    double* trial = nullptr;      ///< n
    S* step = nullptr;            ///< n
    S* dampedStep = nullptr;      ///< n
    Measurement* measurements = nullptr;
    Term* terms = nullptr;
    Eliminated* eliminated = nullptr;
    double* trackChi2 = nullptr;
  };

 public:
  GPUd() PrimaryVertexFitter(std::span<const VertexTrack> tracks, const Material& material, const MeasuredPoint& prior)
    : mTracks(tracks), mMaterial(material), mPrior(prior), mCount(static_cast<std::uint32_t>(tracks.size())), mParameterCount(3 + (3 * mCount))
  {
  }

  /// The bytes of workspace a fit of \a tracks tracks takes.
  [[nodiscard]] GPUhdi() static std::size_t workspaceSize(std::uint32_t tracks)
  {
    Space unused;
    return layout(tracks, nullptr, unused);
  }
  [[nodiscard]] GPUd() std::size_t workspace() const { return workspaceSize(mCount); }

  /// The fit into \a result from the vertex \a start (the momenta from each
  /// track's carry to the plane about it), in \a workspace (at least
  /// `workspace()` bytes, aligned to `kWorkspaceAlignment`; else
  /// kTooLarge). No track, or a prior whose covariance is not positive
  /// definite: kInvalidInput.
  GPUd() void fit(const Vec3D& start, std::span<std::byte> workspace, PrimaryVertexFitResult& result)
  {
    using namespace decay_fit_detail;
    result = {};
    if (mCount < 1 || !priorWeight()) {
      result.status = VertexFitStatus::kInvalidInput;
      return;
    }
    if (!vertex_fit_detail::holdsWorkspace(workspace, workspaceSize(mCount))) {
      result.status = VertexFitStatus::kTooLarge;
      return;
    }
    static_cast<void>(layout(mCount, workspace.data(), mSpace));
    fillN(mSpace.parameters, mParameterCount, 0.);
    copyN(start.data(), 3, mSpace.parameters);
    if (const auto status = carry(parameters()); status != VertexFitStatus::kOk) {
      result.status = status;
      return;
    }
    mNdf = 2 * static_cast<std::int32_t>(mCount);
    const auto refresh = [this]() { return carry(parameters()); };
    const auto minimiseOnce = [this]() { return minimise(); };
    if (!vertex_fit_detail::modelLoop(refresh, minimiseOnce, mMinimiser, result)) {
      return;
    }
    // The covariance at the solution, the measurements of the last carry.
    double rounding = 0.;
    if (!evaluate(parameters(), mChi2, rounding, true)) {
      result.status = VertexFitStatus::kNoCrossing;
      return;
    }
    if (!reduce(0.)) {
      result.status = VertexFitStatus::kSingular;
      return;
    }
    Cholesky<S> llt{.l = mFactor.data()};
    if (!llt.factor(mReduced.data(), 3)) {
      result.status = VertexFitStatus::kSingular;
      return;
    }
    llt.inverse(result.covariance.data());
    copyN(mSpace.parameters, 3, result.vertex.data());
    result.chi2 = mChi2;
    result.priorChi2 = mPriorChi2;
    result.ndf = mNdf;
    mFitted = true;
  }

  /// After `fit`: track \a i's phi, theta, q/p at the vertex.
  [[nodiscard]] GPUd() Vec3D momentum(std::uint32_t i) const { return {mSpace.parameters[3 + (3 * i)], mSpace.parameters[4 + (3 * i)], mSpace.parameters[5 + (3 * i)]}; }

  /// After `fit`: track \a i's chi2 at the solution.
  [[nodiscard]] GPUd() double trackChi2(std::uint32_t i) const { return mSpace.trackChi2[i]; }

  /// After a fit that reached its covariance: the vertex without track \a i
  /// and its covariance (packed), at the solution's linearisation; false
  /// if that leaves no positive-definite S (never with the prior).
  [[nodiscard]] GPUd() bool withoutTrack(std::uint32_t i, Vec3D& vertex, std::array<double, 6>& covariance) const
  {
    using namespace decay_fit_detail;
    if (!mFitted) {
      return false;
    }
    std::array<S, 6> reduced = mReduced;
    std::array<S, 3> rhs = mRhs;
    std::array<S, 6> own{};
    std::array<S, 3> ownRhs{};
    trackShare(i, own, ownRhs);
    for (std::uint32_t k = 0; k < 6; ++k) {
      reduced[k] -= own[k];
    }
    for (std::uint32_t k = 0; k < 3; ++k) {
      rhs[k] -= ownRhs[k];
    }
    std::array<S, 6> factor{};
    Cholesky<S> llt{.l = factor.data()};
    if (!llt.factor(reduced.data(), 3)) {
      return false;
    }
    llt.solve(rhs.data());
    for (std::uint32_t k = 0; k < 3; ++k) {
      vertex[k] = mSpace.parameters[k] + static_cast<double>(rhs[k]);
    }
    llt.inverse(covariance.data());
    return true;
  }

 private:
  // The linear algebra of the block-arrow normal equations for
  // `minimiseFrozen`.
  struct BlockSystem {
    PrimaryVertexFitter* fitter = nullptr;
    [[nodiscard]] GPUd() bool constrained() const { return false; }
    [[nodiscard]] GPUd() bool evaluate(std::span<const double> p, double& chi2, double& rounding, bool normal) const { return fitter->evaluate(p, chi2, rounding, normal); }
    [[nodiscard]] GPUd() VertexFitStatus newton(std::span<const double> /*p*/, S* step, double& /*mu*/, double& /*violation*/, double& /*lambdaC*/) const { return fitter->solve(0., step) ? VertexFitStatus::kOk : VertexFitStatus::kSingular; }
    [[nodiscard]] GPUd() bool damped(double tau, S* d) const { return fitter->solve(tau, d); }
    [[nodiscard]] GPUd() double quadratic(const S* v, double start) const { return fitter->quadratic(v, start); }
    GPUd() void forms(const S* v, double& gd, double& dNd, double& dDd) const { fitter->forms(v, gd, dNd, dDd); }
    [[nodiscard]] GPUd() bool trial(std::span<const double> p, double& chi2, double& violation) const
    {
      double rounding = 0.;
      violation = 0.;
      return fitter->evaluate(p, chi2, rounding, false);
    }
  };

  [[nodiscard]] GPUhdi() static std::size_t layout(std::uint32_t tracks, std::byte* base, Space& space)
  {
    const std::size_t n = std::size_t{3} + (std::size_t{3} * tracks);
    using vertex_fit_detail::place;
    std::size_t at = 0U;
    space.parameters = place<double>(base, at, n);
    space.trial = place<double>(base, at, n);
    space.step = place<S>(base, at, n);
    space.dampedStep = place<S>(base, at, n);
    space.measurements = place<Measurement>(base, at, tracks);
    space.terms = place<Term>(base, at, tracks);
    space.eliminated = place<Eliminated>(base, at, tracks);
    space.trackChi2 = place<double>(base, at, tracks);
    return at;
  }

  [[nodiscard]] GPUd() std::span<double> parameters() const { return {mSpace.parameters, mParameterCount}; }

  // The prior's weight B^-1 into mPriorWeight; false if B is not positive
  // definite.
  [[nodiscard]] GPUd() bool priorWeight()
  {
    using namespace decay_fit_detail;
    std::array<S, 6> factor{};
    Cholesky<S> llt{.l = factor.data()};
    std::array<S, 6> covariance{};
    copyN(mPrior.covariance.data(), 6, covariance.data());
    if (!llt.factor(covariance.data(), 3)) {
      return false;
    }
    llt.inverse(mPriorWeight.data());
    return true;
  }

  // Every track carried to the perigee plane about the vertex of \a at, its
  // measurement there; the momenta's start from the carried states.
  [[nodiscard]] GPUd() VertexFitStatus carry(std::span<const double> at)
  {
    for (std::uint32_t i = 0; i < mCount; ++i) {
      const auto& track = mTracks[i];
      const auto plane = detail::hypothesis::perigeePlaneOf(track.propagator, track.innermost, track.surface, static_cast<Scalar>(at[0]), static_cast<Scalar>(at[1]));
      BoundTrackState state = track.innermost;
      if (const auto status = vertex_fit_detail::carriedMeasurement(track.propagator, mMaterial, track.surface, plane, state, mSpace.measurements[i]); status != VertexFitStatus::kOk) {
        return status;
      }
      if (!mStarted) {
        mSpace.parameters[3 + (3 * i)] = state.phi();
        mSpace.parameters[4 + (3 * i)] = state.theta();
        mSpace.parameters[5 + (3 * i)] = state.qOverP();
      }
    }
    mStarted = true;
    return VertexFitStatus::kOk;
  }

  // The chi2 at \a p (the tracks' and the prior's, each track's into
  // `trackChi2`) and, with \a normal, each track's term and the vertex's
  // block A and right-hand side a; false if a helix misses its plane.
  [[nodiscard]] GPUd() bool evaluate(std::span<const double> p, double& chi2, double& rounding, bool normal)
  {
    using namespace decay_fit_detail;
    chi2 = 0.;
    rounding = 0.;
    for (std::uint32_t i = 0; i < mCount; ++i) {
      const double fieldScale = vertex_fit_detail::fieldScaleOf(mTracks[i].propagator);
      const std::array<double, vertex_fit_detail::kTrackParameters> at{p[0], p[1], p[2], p[3 + (3 * i)], p[4 + (3 * i)], p[5 + (3 * i)]};
      const double before = chi2;
      if (!vertex_fit_detail::trackTerm<M>(mSpace.measurements[i], fieldScale, at, chi2, rounding, normal ? &mSpace.terms[i] : nullptr)) {
        return false;
      }
      mSpace.trackChi2[i] = chi2 - before;
    }
    // The prior: r = b - V, weight B^-1.
    std::array<double, 3> r{};
    for (std::uint32_t k = 0; k < 3; ++k) {
      r[k] = mPrior.position[k] - p[k];
    }
    double prior = 0.;
    for (std::uint32_t k = 0; k < 3; ++k) {
      rounding += vertex_fit_detail::roundingOf(mPrior.position[k]) * static_cast<double>(mPriorWeight[packedIndex(k, k)]);
      for (std::uint32_t l = 0; l < 3; ++l) {
        prior += r[k] * static_cast<double>(mPriorWeight[packedIndex(k, l)]) * r[l];
      }
    }
    mPriorChi2 = prior;
    chi2 += prior;
    if (!normal) {
      return true;
    }
    for (std::uint32_t k = 0; k < 3; ++k) {
      S g = 0.;
      for (std::uint32_t l = 0; l < 3; ++l) {
        g += mPriorWeight[packedIndex(k, l)] * static_cast<S>(r[l]);
      }
      mGradient[k] = g;
      for (std::uint32_t l = 0; l <= k; ++l) {
        mVertexBlock[packedIndex(k, l)] = mPriorWeight[packedIndex(k, l)];
      }
    }
    for (std::uint32_t i = 0; i < mCount; ++i) {
      const auto& term = mSpace.terms[i];
      for (std::uint32_t k = 0; k < 3; ++k) {
        mGradient[k] += term.gradient[k];
        for (std::uint32_t l = 0; l <= k; ++l) {
          mVertexBlock[packedIndex(k, l)] += term.normal[packedIndex(k, l)];
        }
      }
    }
    return true;
  }

  // Track i's E_i (momentum row r, vertex column c) and C_i, a_i, c_i.
  [[nodiscard]] GPUd() S coupling(std::uint32_t i, std::uint32_t r, std::uint32_t c) const { return mSpace.terms[i].normal[decay_fit_detail::packedIndex(3 + r, c)]; }
  [[nodiscard]] GPUd() S momentumBlock(std::uint32_t i, std::uint32_t r, std::uint32_t c) const { return mSpace.terms[i].normal[decay_fit_detail::packedIndex(3 + r, 3 + c)]; }

  // The tracks' blocks eliminated under damping \a tau (N + tau diag N):
  // each track's `Eliminated`, the reduced vertex block S (`mReduced`) and
  // right-hand side h (`mRhs`); false if a block is not positive definite.
  [[nodiscard]] GPUd() bool reduce(double tau)
  {
    using namespace decay_fit_detail;
    const auto t = static_cast<S>(tau);
    for (std::uint32_t k = 0; k < 3; ++k) {
      mRhs[k] = mGradient[k];
      for (std::uint32_t l = 0; l <= k; ++l) {
        mReduced[packedIndex(k, l)] = mVertexBlock[packedIndex(k, l)] + (k == l ? t * mVertexBlock[packedIndex(k, k)] : S(0.));
      }
    }
    for (std::uint32_t i = 0; i < mCount; ++i) {
      auto& e = mSpace.eliminated[i];
      std::array<S, 6> block{};
      for (std::uint32_t r = 0; r < 3; ++r) {
        for (std::uint32_t c = 0; c <= r; ++c) {
          block[packedIndex(r, c)] = momentumBlock(i, r, c) + (r == c ? t * momentumBlock(i, r, r) : S(0.));
        }
      }
      Cholesky<S> llt{.l = e.factor.data()};
      if (!llt.factor(block.data(), 3)) {
        return false;
      }
      for (std::uint32_t c = 0; c < 3; ++c) {
        std::array<S, 3> column{coupling(i, 0, c), coupling(i, 1, c), coupling(i, 2, c)};
        llt.solve(column.data());
        for (std::uint32_t r = 0; r < 3; ++r) {
          e.inverseE[(3 * r) + c] = column[r];
        }
      }
      for (std::uint32_t r = 0; r < 3; ++r) {
        e.inverseC[r] = mSpace.terms[i].gradient[3 + r];
      }
      llt.solve(e.inverseC.data());
      eliminate(i, mReduced, mRhs);
    }
    return true;
  }

  // Track i's eliminated block taken from \a reduced and \a rhs: E_i^T
  // C_i^-1 E_i and E_i^T C_i^-1 c_i.
  GPUd() void eliminate(std::uint32_t i, std::array<S, 6>& reduced, std::array<S, 3>& rhs) const
  {
    using namespace decay_fit_detail;
    const auto& e = mSpace.eliminated[i];
    for (std::uint32_t k = 0; k < 3; ++k) {
      for (std::uint32_t r = 0; r < 3; ++r) {
        rhs[k] -= coupling(i, r, k) * e.inverseC[r];
      }
      for (std::uint32_t l = 0; l <= k; ++l) {
        S sum = 0.;
        for (std::uint32_t r = 0; r < 3; ++r) {
          sum += coupling(i, r, k) * e.inverseE[(3 * r) + l];
        }
        reduced[packedIndex(k, l)] -= sum;
      }
    }
  }

  // Track i's share of the undamped S and h (`reduce(0)`): A_i - E_i^T
  // C_i^-1 E_i and a_i - E_i^T C_i^-1 c_i.
  GPUd() void trackShare(std::uint32_t i, std::array<S, 6>& share, std::array<S, 3>& rhs) const
  {
    using namespace decay_fit_detail;
    const auto& term = mSpace.terms[i];
    for (std::uint32_t k = 0; k < 3; ++k) {
      rhs[k] = term.gradient[k];
      for (std::uint32_t l = 0; l <= k; ++l) {
        share[packedIndex(k, l)] = term.normal[packedIndex(k, l)];
      }
    }
    eliminate(i, share, rhs);
  }

  // The step of N + tau diag N into \a d (n): the vertex's from S, each
  // track's back-substituted; false if a block or S is not positive definite.
  [[nodiscard]] GPUd() bool solve(double tau, S* d)
  {
    using namespace decay_fit_detail;
    if (!reduce(tau)) {
      return false;
    }
    Cholesky<S> llt{.l = mFactor.data()};
    if (!llt.factor(mReduced.data(), 3)) {
      return false;
    }
    std::array<S, 3> dv = mRhs;
    llt.solve(dv.data());
    copyN(dv.data(), 3, d);
    for (std::uint32_t i = 0; i < mCount; ++i) {
      const auto& e = mSpace.eliminated[i];
      for (std::uint32_t r = 0; r < 3; ++r) {
        S value = e.inverseC[r];
        for (std::uint32_t c = 0; c < 3; ++c) {
          value -= e.inverseE[(3 * r) + c] * dv[c];
        }
        d[3 + (3 * i) + r] = value;
      }
    }
    return true;
  }

  // start + v^T N v over the blocks (the undamped N of the last evaluation).
  [[nodiscard]] GPUd() double quadratic(const S* v, double start) const
  {
    double sum = start;
    for (std::uint32_t k = 0; k < 3; ++k) {
      for (std::uint32_t l = 0; l < 3; ++l) {
        sum += static_cast<double>(v[k]) * static_cast<double>(mVertexBlock[decay_fit_detail::packedIndex(k, l)]) * static_cast<double>(v[l]);
      }
    }
    for (std::uint32_t i = 0; i < mCount; ++i) {
      const S* q = v + 3 + (std::size_t{3} * i);
      for (std::uint32_t r = 0; r < 3; ++r) {
        for (std::uint32_t c = 0; c < 3; ++c) {
          sum += 2. * static_cast<double>(q[r]) * static_cast<double>(coupling(i, r, c)) * static_cast<double>(v[c]);
          sum += static_cast<double>(q[r]) * static_cast<double>(momentumBlock(i, r, c)) * static_cast<double>(q[c]);
        }
      }
    }
    return sum;
  }

  GPUd() void forms(const S* v, double& gd, double& dNd, double& dDd) const
  {
    gd = dDd = 0.;
    for (std::uint32_t k = 0; k < 3; ++k) {
      gd += static_cast<double>(mGradient[k]) * static_cast<double>(v[k]);
      dDd += static_cast<double>(mVertexBlock[decay_fit_detail::packedIndex(k, k)]) * static_cast<double>(v[k]) * static_cast<double>(v[k]);
    }
    for (std::uint32_t i = 0; i < mCount; ++i) {
      const S* q = v + 3 + (std::size_t{3} * i);
      for (std::uint32_t r = 0; r < 3; ++r) {
        gd += static_cast<double>(mSpace.terms[i].gradient[3 + r]) * static_cast<double>(q[r]);
        dDd += static_cast<double>(momentumBlock(i, r, r)) * static_cast<double>(q[r]) * static_cast<double>(q[r]);
      }
    }
    dNd = quadratic(v, 0.);
  }

  [[nodiscard]] GPUd() VertexFitStatus minimise()
  {
    BlockSystem system{.fitter = this};
    return vertex_fit_detail::minimiseFrozen<M>(system, mParameterCount, mNdf, parameters(), mSpace.trial, mSpace.step, mSpace.dampedStep, mChi2, mMinimiser);
  }

  std::span<const VertexTrack> mTracks;
  const Material& mMaterial; // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
  MeasuredPoint mPrior;
  std::uint32_t mCount = 0;
  std::uint32_t mParameterCount = 0;
  std::int32_t mNdf = 0;
  bool mStarted = false;
  bool mFitted = false;
  std::array<S, 6> mPriorWeight{};
  std::array<S, 6> mVertexBlock{}; ///< A, the prior's weight in it
  std::array<S, 3> mGradient{};    ///< a
  std::array<S, 6> mReduced{};     ///< S
  std::array<S, 3> mRhs{};         ///< h
  std::array<S, 6> mFactor{};
  double mChi2 = 0.;
  double mPriorChi2 = 0.;
  vertex_fit_detail::MinimiserState mMinimiser;
  Space mSpace;
};

/// A `PrimaryVertexFitter`'s fit of one vertex from \a start, in a workspace
/// of its own; \a momenta and \a trackChi2, if not empty (one per track),
/// get each track's momentum at the vertex and its chi2.
template <typename Material>
[[nodiscard]] PrimaryVertexFitResult fitPrimaryVertex(std::span<const VertexTrack> tracks, const Material& material, const MeasuredPoint& prior, const Vec3D& start, std::span<Vec3D> momenta = {}, std::span<double> trackChi2 = {})
{
  PrimaryVertexFitter<Material> fitter(tracks, material, prior);
  WorkspaceBuffer storage;
  PrimaryVertexFitResult result;
  fitter.fit(start, storage.take(fitter.workspace()), result);
  if (result.status == VertexFitStatus::kOk) {
    for (std::uint32_t i = 0; i < momenta.size(); ++i) {
      momenta[i] = fitter.momentum(i);
    }
    for (std::uint32_t i = 0; i < trackChi2.size(); ++i) {
      trackChi2[i] = fitter.trackChi2(i);
    }
  }
  return result;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_PRIMARYVERTEXFIT_H_
