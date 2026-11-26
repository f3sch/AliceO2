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

/// \file MotherSearch.h
/// \brief A decay's charged mother's own rows found from its fitted decay
///        ("strangeness tracking": a cascade's Xi or Omega on the layers it
///        crossed before it decayed).
///
/// The mother's state at its fitted vertex (`motherStateAt`): its helix from
/// the vertex with its daughters' summed momentum and charge, bound to the
/// plane through the vertex normal to its transverse direction, its
/// covariance the fit's propagated (J C J^T over the vertex and the
/// momenta). From there inward, against the momentum, onto every layer --
/// barrel and disc -- the straight line back along the mother's direction to
/// the beam axis crosses within the layer's span, the nearest to the vertex
/// first, the follow's own step: the nominal crossing
/// (`predictLayerCrossing`), the window of the carried covariance with the
/// material on the way (`follow::predictedWindow`), the rows in the time
/// box, each updated on its own sensor (`follow::updateCandidate`); the
/// smallest chi2 under the gate taken, the state going on from it. The
/// first step starts from the vertex's plane, of no layer: its material the
/// table's gap share (`stepMaterial`). One row per layer; rows of the
/// decay's own tracks never.

#ifndef GTRK_STUDY_MOTHERSEARCH_H_
#define GTRK_STUDY_MOTHERSEARCH_H_

#include "GenericTracking/Tracking/BeamConstraint.h"
#include "GenericTracking/Tracking/DecayFitMath.h"
#include "GenericTracking/Tracking/DecayVertexFit.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/LinkHoles.h"
#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Tracking/TimeFrame.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace o2::gtrk::study
{

/// A decay's charged mother at its fitted vertex (`motherStateAt`).
struct MotherState {
  BoundTrackState state; ///< on `plane`, no sensor
  TrackingSurface plane; ///< through the vertex, normal to the mother's transverse direction; no layer
  bool found = false;
};

/// The mother decaying at vertex \a w of \a fit (\a topology's): its
/// daughters under \a w summed, their momentum and charge, from the vertex,
/// bound to the plane through it normal to that transverse direction
/// (`crossToPlane`'s dual numbers over the fit's parameters), its covariance
/// J C J^T from the fit's; none for a neutral one or where the fit holds no
/// covariance. \a fieldScale: Bz B2C.
[[nodiscard]] inline MotherState motherStateAt(const DecayVertexFitResult& fit, const DecayTopology& topology, std::uint32_t w, double fieldScale)
{
  using namespace decay_fit_detail;
  constexpr std::int32_t kInputs = (3 * kMaxConstrainedVertices) + (3 * kMaxConstrainedDaughters);
  using D = Dual<kInputs>;
  MotherState out;
  const auto vertices = fit.vertices.size();
  const auto count = fit.momenta.size();
  const auto inputs = static_cast<std::int32_t>((3 * vertices) + (3 * count));
  if (std::cmp_not_equal(fit.covariance.size(), covSize(inputs))) {
    return out;
  }
  FreeState<kInputs> from{};
  const auto vertex = vertexOfInputs<kInputs>(fit.vertices[w], w);
  for (std::uint32_t c = 0; c < 3; ++c) {
    from[c] = vertex[c];
  }
  MomentumSum<D> sum;
  double charge = 0.;
  for (std::uint32_t i = 0; i < count; ++i) {
    if (topology.under(i, w)) {
      const auto& m = fit.momenta[i];
      sum.add(momentumOfInputs<kInputs>(m[0], m[1], m[2], (3 * vertices) + (3 * i), fit.charges[i]).vector);
      charge += fit.charges[i];
    }
  }
  if (exactlyEqual(charge, 0.)) {
    return out;
  }
  const auto direction = directionOfMomentum(sum.vector);
  from[3] = direction.phi;
  from[4] = direction.theta;
  from[5] = D{charge} / direction.p;
  const auto& at = fit.vertices[w];
  out.plane = beamLinePlane(static_cast<Scalar>(at[0]), static_cast<Scalar>(at[1]), static_cast<Scalar>(direction.phi.v));
  FreeState<kInputs> to{};
  std::array<D, 2> local{};
  if (!crossToPlane<kInputs>(from, out.plane, fieldScale, to, local)) {
    return out;
  }
  const std::array<const D*, kNBoundParams> bound{local.data(), &local[1], &to[3], &to[4], &to[5]};
  for (std::int32_t a = 0; a < kNBoundParams; ++a) {
    out.state.par[a] = static_cast<Scalar>(bound[static_cast<std::uint32_t>(a)]->v);
  }
  for (std::int32_t a = 0; a < kNBoundParams; ++a) {
    for (std::int32_t b = 0; b <= a; ++b) {
      double value = 0.;
      for (std::int32_t k = 0; k < inputs; ++k) {
        double row = 0.;
        for (std::int32_t l = 0; l < inputs; ++l) {
          row += fit.covarianceAt(k, l) * bound[static_cast<std::uint32_t>(b)]->d[static_cast<std::uint32_t>(l)];
        }
        value += bound[static_cast<std::uint32_t>(a)]->d[static_cast<std::uint32_t>(k)] * row;
      }
      out.state.covAt(a, b) = static_cast<Scalar>(value);
    }
  }
  out.found = true;
  return out;
}

/// A mother's rows found (`MotherSearch::search`), the outermost first.
struct MotherRows {
  std::vector<std::uint32_t> rows;
  std::vector<float> chi2;  ///< each row's update chi2, before the update
  std::uint8_t layers = 0U; ///< the layers searched: crossed inside the vertex
};

/// The search's fixed parts for one timeframe: the layers' spans and the
/// gate, the widest time box per layer.
template <DetectorLayout Det>
class MotherSearch
{
 public:
  static constexpr std::size_t kLayerCount = static_cast<std::size_t>(Detector<Det>::nLayers);

  /// \a nSigma: the gate on each row, `chiSquareQuantile(2, nSigma)`.
  MotherSearch(const TimeFrame<Det>& frame, const TrackingGeometry<Det>& geometry, double nSigma) : mFrame(&frame), mGeometry(&geometry)
  {
    for (std::uint32_t layer = 0; layer < kLayerCount; ++layer) {
      mParameters.layerSpans[layer] = layerSpan<Det>(static_cast<LayerID>(layer));
      const auto range = frame.grid().layerRows(static_cast<LayerID>(layer));
      for (auto row = range.begin; row < range.end; ++row) {
        mMaxWidthPS[layer] = std::max(mMaxWidthPS[layer], timeBoxAt(frame.measurements(), row).width);
      }
    }
    mParameters.barrelLayers = barrelLayerBits<Det>();
    mParameters.clusterChi2Max = static_cast<Scalar>(chiSquareQuantile(2, nSigma));
  }

  /// From \a start (`motherStateAt`) with \a propagator (the mother's mass),
  /// in the time box \a box (the decay's tracks'), every row but those
  /// \a skip(row) names.
  template <typename Skip>
  [[nodiscard]] MotherRows search(const MotherState& start, const Propagator& propagator, TimeEstPS box, const Skip& skip) const
  {
    MotherRows out;
    BoundParameters state = static_cast<const BoundParameters&>(start.state);
    auto factors = factorCovariance(start.state.cov);
    TrackingSurface from = start.plane;
    const auto& frame = *mFrame;
    const auto& geometry = *mGeometry;
    for (const auto layer : layersBehind(start)) {
      const auto id = static_cast<LayerID>(layer);
      const auto nominal = predictLayerCrossing(propagator, mParameters, from, state, id, PropagationDirection::kBackward);
      if (nominal.status != PropagationStatus::kOk) {
        continue;
      }
      Scalar phiCentre = 0.F;
      Scalar secondCentre = 0.F;
      follow::windowCentre(mParameters, id, nominal, phiCentre, secondCentre);
      CandidateWindow window{
        .earliestPS = static_cast<std::int64_t>(box.lower()) - static_cast<std::int64_t>(mMaxWidthPS[static_cast<std::uint32_t>(layer)]) + 1,
        .latestPS = static_cast<std::int64_t>(box.upper()) - 1,
        .afterPS = box.lower(),
        .phiCentre = phiCentre,
        .phiHalfWidth = kPi<Scalar>,
        .secondCentre = secondCentre,
        .secondHalfWidth = follow::kOpenBound,
      };
      if (!follow::predictedWindow(propagator, mParameters, geometry, detail::LinkSteps{}, from, state, factors, id, PropagationDirection::kBackward, nominal, window)) {
        continue;
      }
      ++out.layers;
      std::uint32_t best = std::numeric_limits<std::uint32_t>::max();
      Scalar bestChi2 = mParameters.clusterChi2Max;
      BoundParameters bestState;
      CovarianceFactors bestFactors;
      frame.grid().visitCandidates(id, window, frame.grid().geometry().allBlocks(), [&](std::uint32_t row) {
        if (skip(row)) {
          return;
        }
        const auto m = frame.measurement(row);
        if (!spatiallyCompatible(window, m.phi, secondCoordinate<Det>(id, m.r, m.z))) {
          return;
        }
        BoundParameters updated;
        CovarianceFactors updatedFactors;
        Scalar chi2 = 0.F;
        if (!follow::updateCandidate(propagator, geometry, detail::LinkSteps{}, from, measurementAt(frame, row), PropagationDirection::kBackward, state, factors, updated, updatedFactors, chi2)) {
          return;
        }
        if (chi2 < bestChi2) {
          bestChi2 = chi2;
          best = row;
          bestState = updated;
          bestFactors = updatedFactors;
        }
      });
      if (best == std::numeric_limits<std::uint32_t>::max()) {
        continue;
      }
      out.rows.push_back(best);
      out.chi2.push_back(bestChi2);
      state = bestState;
      factors = bestFactors;
      from = geometry.planeOf(state);
      box += timeBoxAt(frame.measurements(), best);
    }
    return out;
  }

 private:
  /// The layers between \a start's vertex and the beam axis, on the straight
  /// line back along the mother's direction, crossed within their span (a
  /// barrel layer whose sensors begin inside the vertex's radius, a disc
  /// between the vertex's z and the line's at the axis), the nearest to the
  /// vertex first.
  [[nodiscard]] std::vector<std::int32_t> layersBehind(const MotherState& start) const
  {
    using Math = o2::gpu::CAMath;
    const Scalar vertexR = Math::Sqrt((start.plane.origin[0] * start.plane.origin[0]) + (start.plane.origin[1] * start.plane.origin[1]));
    const Scalar vertexZ = start.state.loc1();
    const Scalar cotTheta = Math::Cos(start.state.theta()) / Math::Sin(start.state.theta());
    const Scalar axisZ = vertexZ - (vertexR * cotTheta);
    std::vector<std::pair<Scalar, std::int32_t>> crossed;
    for (std::int32_t layer = 0; layer < Detector<Det>::nLayers; ++layer) {
      const auto& span = mParameters.layerSpans[layer];
      Scalar fraction = 0.F; // from the vertex towards the axis
      if (((mParameters.barrelLayers >> layer) & 1U) != 0U) {
        if (!(span.positionMin < vertexR)) {
          continue;
        }
        fraction = Math::Max(vertexR - span.position, Scalar{0}) / vertexR;
      } else {
        const Scalar dz = vertexZ - axisZ;
        if (!(Math::Abs(dz) > 0.F) || !((span.position - axisZ) / dz > 0.F && (span.position - axisZ) / dz < 1.F)) {
          continue;
        }
        fraction = (vertexZ - span.position) / dz;
      }
      const Scalar second = ((mParameters.barrelLayers >> layer) & 1U) != 0U ? vertexZ - (fraction * vertexR * cotTheta) : vertexR * (1.F - fraction);
      if (second < span.low || second > span.high) {
        continue;
      }
      crossed.emplace_back(fraction, layer);
    }
    std::ranges::sort(crossed);
    std::vector<std::int32_t> layers;
    layers.reserve(crossed.size());
    for (const auto& [fraction, layer] : crossed) {
      layers.push_back(layer);
    }
    return layers;
  }

  const TimeFrame<Det>* mFrame = nullptr;
  const TrackingGeometry<Det>* mGeometry = nullptr;
  IterationParameters<Det> mParameters{};
  std::array<std::uint32_t, kLayerCount> mMaxWidthPS{};
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_MOTHERSEARCH_H_
