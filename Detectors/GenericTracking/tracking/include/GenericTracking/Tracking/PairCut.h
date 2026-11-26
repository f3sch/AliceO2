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

/// \file PairCut.h
/// \brief The exact spatial cuts a tracklet passes.
///
/// Two tests on the measured positions, in the frame of the corrected beam
/// line: a common helix pitch from the vertex to both measurements, and the
/// curvature of the circle through the beam point and both measurements. One
/// expression serves every kind of layer pair, and each measurement's own
/// position and errors enter rather than a layer's nominal ones.
///
/// The cheap transverse test runs first. Only its survivors pay for the
/// division and arcsines of the helix-pitch test. The loop over layers that
/// multiple scattering needs happens once per source and vertex, in
/// `sourceScattering`.
///
/// The iteration's cut and the vertex are separate arguments: one source is
/// tested against every vertex seed compatible with it in time, so the vertex
/// changes within one source's work while the cut stays the iteration's.

#ifndef GTRK_TRACKING_PAIRCUT_H_
#define GTRK_TRACKING_PAIRCUT_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/LinkMaterial.h"
#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Tracking/SearchWindow.h"
#include "GenericTracking/Tracking/WorkCounters.h"

#include "GPUCommonMath.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>

namespace o2::gtrk
{

/// `(0.0136 GeV)^2 x (1 + 0.038 ln(x / beta^2))^2`: Highland's theta0^2
/// (beta p)^2 at normal incidence (`highlandTheta0`), for \a x2X0 = x/X0 and
/// \a logRatio = ln(x / beta^2), in T; 0 without material or where the
/// logarithm turns it negative.
template <typename T>
  requires std::is_floating_point_v<T>
[[nodiscard]] GPUhdi() GPUAlwaysInline() T highlandK(T x2X0, std::type_identity_t<T> logRatio)
{
  const T highland = T(1) + (T(0.038) * logRatio);
  return x2X0 > T(0) && highland > T(0) ? T(0.0136) * T(0.0136) * x2X0 * highland * highland : T(0);
}

namespace detail
{

/// `highlandK` for \a x = x/X0 and \a beta2 = beta^2, in double.
[[nodiscard]] inline double highlandFactor(double x, double beta2) { return x > 0. && beta2 > 0. ? highlandK(x, std::log(x / beta2)) : 0.; }

} // namespace detail

/// The cut parameters of one iteration, the same for every layer pair and
/// every vertex.
struct PairCut {
  Scalar twoRMin = 0.F;          ///< twice the softest track's transverse radius [cm], the vertexer's; 0 with no field
  Scalar nSigma = 3.F;           ///< how many sigma of the Gaussian terms are admitted
  Scalar invBetaMomentum2 = 0.F; ///< 1 / (beta pT)^2 at pT_min [GeV^-2], for the scattering sum (`invBetaMomentumMin2`)
  /// The largest transverse curvature above pT_min [1/cm]: zero with no field,
  /// unbounded with no pT cut.
  Scalar maxCurvature = std::numeric_limits<Scalar>::infinity();
  /// A track's momentum times its 3D curvature, |B2C Bz| [GeV/cm]: what turns
  /// a fitted curvature into the momentum its scattering is taken at; 0 with
  /// no field.
  Scalar momentumCurvature = 0.F;
  Scalar mass2 = 0.F; ///< the windows' mass squared [GeV^2] (`kWindowMassGeV`), for beta
  /// ln(beta^2) of the momentum `invBetaMomentum2` is at, for the Highland
  /// factors' logarithm (`highlandK`); 0, beta = 1, with no pT cut.
  Scalar logBeta2 = 0.F;
};

/// ln(beta^2) at pT_min for the windows' mass: the least of an iteration's
/// tracks, at which its cuts take the scattering (`invBetaMomentumMin2`).
[[nodiscard]] inline double logBeta2Min(const SearchWindowParams& params)
{
  const double pt2 = params.minPtGeV * params.minPtGeV;
  return params.minPtGeV > 0. ? std::log(pt2 / (pt2 + (kWindowMassGeV * kWindowMassGeV))) : 0.;
}

[[nodiscard]] inline PairCut pairCut(const SearchWindowParams& params)
{
  const double radius = curvatureRadius(params.minPtGeV, params.bzKiloGauss);
  return {
    .twoRMin = static_cast<Scalar>(2. * radius),
    .nSigma = static_cast<Scalar>(params.nSigma),
    .invBetaMomentum2 = params.minPtGeV > 0. ? static_cast<Scalar>(invBetaMomentumMin2(params)) : 0.F,
    .maxCurvature = static_cast<Scalar>(maxCurvature(params)),
    .momentumCurvature = static_cast<Scalar>(std::abs(static_cast<double>(o2::constants::math::B2C) * params.bzKiloGauss)),
    .mass2 = static_cast<Scalar>(kWindowMassGeV * kWindowMassGeV),
    .logBeta2 = static_cast<Scalar>(logBeta2Min(params)),
  };
}

/// A candidate origin of the tracks in one collision: the beam line, one
/// entry for the whole timeframe, or a vertex from a list sorted by time.
///
/// Transverse positions are in the frame of the corrected beam line, so a
/// vertex contributes only its transverse spread there.
struct VertexSeed {
  /// The time box `[timePS, timePS + widthPS)` [ps from the timeframe's start]:
  /// a point for an MC vertex, the span of
  /// its contributors' boxes for a found one, every time for the beam line.
  std::uint32_t timePS = 0U;
  std::uint32_t widthPS = 0U;
  Scalar z = 0.F;             ///< [cm]
  Scalar varZ = 0.F;          ///< [cm^2]
  Scalar varTransverse = 0.F; ///< transverse spread squared [cm^2]
};

/// A seed from a spatial constraint alone, compatible with every time: the
/// beam line, from the beam spot's z and widths.
[[nodiscard]] inline VertexSeed vertexSeed(const VertexConstraint& constraint)
{
  return {
    .timePS = 0U,
    .widthPS = std::numeric_limits<std::uint32_t>::max(),
    .z = static_cast<Scalar>(constraint.z),
    .varZ = static_cast<Scalar>(constraint.sigmaZ * constraint.sigmaZ),
    .varTransverse = static_cast<Scalar>(constraint.sigmaTransverse * constraint.sigmaTransverse),
  };
}

/// What the tests read of one measurement: its shipped position and its
/// variances along z, r and r-phi.
struct PairPoint {
  Scalar x = 0.F;
  Scalar y = 0.F;
  Scalar z = 0.F;
  Scalar r = 0.F;
  Scalar varZ = 0.F;
  Scalar varR = 0.F;
  Scalar varRPhi = 0.F;
};

template <typename Row>
[[nodiscard]] GPUhdi() GPUAlwaysInline() PairPoint pairPoint(const Row& row)
{
  return {.x = row.x, .y = row.y, .z = row.z, .r = row.r, .varZ = row.varZ, .varR = row.varR, .varRPhi = row.varRPhi};
}

/// The same for row \a row of \a columns, reading only these seven columns: a
/// whole row through a view's `operator[]` is a call that reads all of them.
template <typename Columns>
[[nodiscard]] GPUhdi() GPUAlwaysInline() PairPoint pairPointAt(const Columns& columns, std::uint32_t row)
{
  return {
    .x = columns.x[row],
    .y = columns.y[row],
    .z = columns.z[row],
    .r = columns.r[row],
    .varZ = columns.varZ[row],
    .varR = columns.varR[row],
    .varRPhi = columns.varRPhi[row],
  };
}

/// `S = sum_k theta_k^2 s_k^2` over the material the source's line from the
/// vertex crosses up to and including the source's own layer, for the
/// iteration's minimum pT (LinkMaterial.h): the source layer's vertex gap,
/// a share x of it at the fraction f of the way, and the layer's node at the
/// source.
///
/// A kink at distance `s_k` from the vertex moves the collinearity by
/// `theta_k s_k (s_j - s_i)` and the curvature cross product by
/// `theta_k s_k d`, so this one sum serves both tests. With `p = pT g` and
/// `s = f L`, `L = r g`: `S = r^2 / (beta pT)^2 [K(x_gap) <f^2> + K(x_node)]`,
/// each part's Highland logarithm at its own x along the line, at the cut's
/// beta (`PairCut::logBeta2`). The node's path: the line crosses a barrel
/// layer radially, a disc along z.
/// \a sourcePhi: the source's azimuth, its measurement's. Highland's
/// logarithm is the span's from the vertex through the source, its x/X0 in
/// \a spanX: a pair adds its gap to it (`betweenScattering`).
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar sourceScattering(const LinkMaterial<Det>& material, const PairCut& cut, const VertexSeed& vertex, LayerID sourceLayer, const PairPoint& source, Scalar sourcePhi, Scalar& spanX)
{
  spanX = 0.F;
  if (!(source.r > 0.F)) {
    return 0.F;
  }
  const Scalar dz = source.z - vertex.z;
  const Scalar second = LinkMaterial<Det>::second(sourceLayer, source.r, source.z);
  const auto& gap = material.vertexGapAt(sourceLayer, sourcePhi, vertex.z, second);
  const auto& node = material.nodeAt(sourceLayer, sourcePhi, second);
  const Scalar length = o2::gpu::GPUCommonMath::Sqrt((source.r * source.r) + (dz * dz));
  const Scalar absDz = o2::gpu::GPUCommonMath::Abs(dz);
  const Scalar incidence = Detector<Det>::isBarrelLayer(sourceLayer) ? length / source.r : (absDz > 0.F ? length / absDz : 0.F);
  // The gap and the node at normal incidence, both at the line's.
  const Scalar nodeX = node.x * incidence;
  const Scalar gapX = gap.x * incidence;
  spanX = gapX + nodeX;
  if (!(spanX > 0.F)) {
    return 0.F;
  }
  const Scalar logSpan = logOf(spanX);
  return source.r * source.r * cut.invBetaMomentum2 * ((materialKick(gapX, logSpan, cut.logBeta2) * gap.m2) + materialKick(nodeX, logSpan, cut.logBeta2));
}

/// The same without the span, the azimuth from the point.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar sourceScattering(const LinkMaterial<Det>& material, const PairCut& cut, const VertexSeed& vertex, LayerID sourceLayer, const PairPoint& source)
{
  Scalar spanX = 0.F;
  return sourceScattering(material, cut, vertex, sourceLayer, source, atan2Of(source.y, source.x), spanX);
}

/// Highland's factor `(1 + 0.038 ln(x / beta^2))^2` of a span \a x over the
/// same of the span \a from it grew from: what a scattering taken at \a from's
/// logarithm becomes at \a x's; 1 where \a x adds nothing.
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar spanRatio(Scalar from, Scalar x, Scalar logBeta2)
{
  if (!(from > 0.F) || !(x > from)) {
    return 1.F;
  }
  const Scalar low = 1.F + (0.038F * (logOf(from) - logBeta2));
  const Scalar high = 1.F + (0.038F * (logOf(x) - logBeta2));
  return low > 0.F ? (high * high) / (low * low) : 1.F;
}

/// What `betweenScattering` takes of a source once for all its links and
/// candidates: its gap azimuth bin and its second coordinate.
template <DetectorLayout Det>
struct BetweenSource {
  std::uint32_t phiBin = 0U;
  Scalar second = 0.F;
};

template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() BetweenSource<Det> betweenSource(const PairPoint& source, LayerID sourceLayer, Scalar sourcePhi)
{
  return {.phiBin = static_cast<std::uint32_t>(LinkMaterial<Det>::phiBinOf(sourcePhi, LinkMaterial<Det>::kGapPhi)), .second = LinkMaterial<Det>::second(sourceLayer, source.r, source.z)};
}

template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() BetweenSource<Det> betweenSource(const PairPoint& source, LayerID sourceLayer)
{
  return betweenSource<Det>(source, sourceLayer, atan2Of(source.y, source.x));
}

/// A pair's between term and the factor its source's term takes at the
/// pair's span (`betweenScattering`).
struct PairScattering {
  Scalar between = 0.F;
  Scalar sourceScale = 1.F;
};

/// `T = sum_p theta_p^2 (s_j - s_p)^2` over the link's gap between the source
/// and the target (LinkMaterial.h), for the iteration's minimum pT: a kink
/// there moves the target alone. With the gap's share of x at fraction f of
/// the chord L from the source, `T = K(x) (1 - 2 <f> + <f^2>) L^2 / (beta p)^2`,
/// `p = pT g` from the source's line from the vertex. The source's own layer
/// is `sourceScattering`'s, the target's own does not move the target.
/// Highland's logarithm is the pair's span's, the source's \a sourceSpanX
/// (`sourceScattering`) and the gap: `sourceScale` carries the source's term
/// there.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() PairScattering betweenScattering(const LinkMaterial<Det>& material, const PairCut& cut, const BetweenSource<Det>& from, std::uint32_t link, LayerID targetLayer, const VertexSeed& vertex, const PairPoint& source, const PairPoint& target, Scalar sourceSpanX, Scalar sourceHighland)
{
  if (!(source.r > 0.F)) {
    return {};
  }
  const Scalar slope = (source.z - vertex.z) / source.r;
  const Scalar g2 = 1.F + (slope * slope);
  const auto& gap = material.gap(material.gapRowAt(link, from.phiBin, from.second), targetLayer, LinkMaterial<Det>::second(targetLayer, target.r, target.z));
  const Scalar dx = target.x - source.x;
  const Scalar dy = target.y - source.y;
  const Scalar dz = target.z - source.z;
  const Scalar lever = 1.F - (2.F * gap.m1) + gap.m2;
  const Scalar spanX = sourceSpanX + gap.x;
  if (!(spanX > 0.F)) {
    return {};
  }
  const Scalar logSpan = logOf(spanX);
  const Scalar highland = 1.F + (0.038F * (logSpan - cut.logBeta2));
  const Scalar ratio = sourceHighland > 0.F && sourceSpanX > 0.F ? highland / sourceHighland : 1.F;
  return {
    .between = materialKick(gap.x, logSpan, cut.logBeta2) * lever * ((dx * dx) + (dy * dy) + (dz * dz)) * cut.invBetaMomentum2 / g2,
    .sourceScale = ratio * ratio,
  };
}

/// The same, the source's Highland factor from its span.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() PairScattering betweenScattering(const LinkMaterial<Det>& material, const PairCut& cut, const BetweenSource<Det>& from, std::uint32_t link, LayerID targetLayer, const VertexSeed& vertex, const PairPoint& source, const PairPoint& target, Scalar sourceSpanX)
{
  return betweenScattering(material, cut, from, link, targetLayer, vertex, source, target, sourceSpanX, sourceSpanX > 0.F ? 1.F + (0.038F * (logOf(sourceSpanX) - cut.logBeta2)) : 0.F);
}

/// The least x/X0 a source's span can have on layer \a from (mean over
/// azimuth, as `largestSourceSpan` takes the largest): where a gap's share of
/// a pair's logarithm is largest.
template <DetectorLayout Det>
[[nodiscard]] inline double leastSourceSpan(const LinkMaterial<Det>& material, LayerID from)
{
  using Table = LinkMaterial<Det>;
  double least = std::numeric_limits<double>::infinity();
  for (std::int32_t s = 0; s < Table::kSecond; ++s) {
    double node = std::numeric_limits<double>::infinity();
    for (std::int32_t phi = 0; phi < Table::kNodePhi; ++phi) {
      node = std::min(node, static_cast<double>(material.nodes[(((std::uint32_t{from} * Table::kNodePhi) + static_cast<std::uint32_t>(phi)) * Table::kSecond) + static_cast<std::uint32_t>(s)].x));
    }
    for (std::int32_t v = 0; v < Table::kVertexZ; ++v) {
      for (std::int32_t phi = 0; phi < Table::kGapPhi; ++phi) {
        least = std::min(least, node + static_cast<double>(material.vertexGaps[(((((std::uint32_t{from} * Table::kGapPhi) + static_cast<std::uint32_t>(phi)) * Table::kVertexZ) + static_cast<std::uint32_t>(v)) * Table::kSecond) + static_cast<std::uint32_t>(s)].x));
      }
    }
  }
  return least;
}

/// The largest path factor through layer \a from's band (path over the band's
/// depth along the line from the vertex) of the sources in its second
/// coordinate bin \a s from vertices in vertex bin \a v within \a reach: over
/// the corners of the bin, the vertices' interval and the band's depth, where
/// the factor, monotonic in each, is largest. Zero if no vertex of the reach
/// is in the bin.
template <DetectorLayout Det>
[[nodiscard]] inline double largestIncidence(const LinkMaterial<Det>& material, LayerID from, std::int32_t s, std::int32_t v, const Interval& reach)
{
  const double v0 = std::max(reach.low, static_cast<double>(material.vertexZ.low) + static_cast<double>(static_cast<float>(v) / material.vertexZ.perUnit));
  const double v1 = std::min(reach.high, static_cast<double>(material.vertexZ.low) + static_cast<double>(static_cast<float>(v + 1) / material.vertexZ.perUnit));
  if (v0 > v1) {
    return 0.;
  }
  const bool barrel = Detector<Det>::isBarrelLayer(from);
  const auto& axis = material.layerSecond[from];
  const double s0 = axis.low + (static_cast<float>(s) / axis.perUnit);
  const double s1 = axis.low + (static_cast<float>(s + 1) / axis.perUnit);
  // The band's faces (LinkMaterial::layerDepth holds them unsigned).
  const auto& depth = material.layerDepth[from];
  const double inner = depth.low;
  const double outer = static_cast<double>(depth.low) + (1. / static_cast<double>(depth.perUnit));
  const bool forward = barrel || PlacedEnvelopes<Det>::get(from).position > 0.;
  const double d0 = forward ? inner : -outer;
  const double d1 = forward ? outer : -inner;
  double largest = 0.;
  for (const double second : {s0, s1}) {
    for (const double face : {d0, d1}) {
      for (const double zv : {v0, v1}) {
        const double r = barrel ? face : second;
        const double dz = (barrel ? second : face) - zv;
        const double length = std::hypot(r, dz);
        const double across = barrel ? r : std::abs(dz);
        if (across > 0.) {
          largest = std::max(largest, length / across);
        }
      }
    }
  }
  return largest;
}

/// The largest x/X0 a source's span can have on layer \a from
/// (`sourceScattering`'s \a spanX) for vertices in the reach: what bounds a
/// pair's logarithm before it is read.
template <DetectorLayout Det>
[[nodiscard]] inline double largestSourceSpan(const LinkMaterial<Det>& material, const SearchWindowParams& params, LayerID from, const VertexConstraint& vertex)
{
  using Table = LinkMaterial<Det>;
  const Interval reach = vertexReach(params, vertex);
  double largest = 0.;
  for (std::int32_t s = 0; s < Table::kSecond; ++s) {
    double node = 0.;
    for (std::int32_t phi = 0; phi < Table::kNodePhi; ++phi) {
      node = std::max(node, static_cast<double>(material.nodes[(((std::uint32_t{from} * Table::kNodePhi) + static_cast<std::uint32_t>(phi)) * Table::kSecond) + static_cast<std::uint32_t>(s)].x));
    }
    for (std::int32_t v = 0; v < Table::kVertexZ; ++v) {
      const double incidence = largestIncidence<Det>(material, from, s, v, reach);
      if (!(incidence > 0.)) {
        continue;
      }
      for (std::int32_t phi = 0; phi < Table::kGapPhi; ++phi) {
        const auto& gap = material.vertexGaps[(((((std::uint32_t{from} * Table::kGapPhi) + static_cast<std::uint32_t>(phi)) * Table::kVertexZ) + static_cast<std::uint32_t>(v)) * Table::kSecond) + static_cast<std::uint32_t>(s)];
        largest = std::max(largest, (static_cast<double>(gap.x) + node) * incidence);
      }
    }
  }
  return largest;
}

/// One standard deviation of a source's transverse direction on \a from
/// [rad], at its largest over the layer and the vertex range: `sqrt(S)/r_i`,
/// with `S` the `sourceScattering` sum at pT_min, and the beam width over the
/// radius. The transverse test takes both times the pair's chord.
///
/// `S / r^2 = [K(x_gap) <f^2> + K(x_node)] / (beta pT)^2` holds the layer's
/// radius out: the largest bracket over the layer's coordinates and the
/// vertices within the reach, its gap and node each the mean over azimuth --
/// a window holds the material the tracks see, as the per-layer averages did,
/// not every stave edge's.
template <DetectorLayout Det>
[[nodiscard]] inline double sourceKinkBound(const LinkMaterial<Det>& material, const SearchWindowParams& params, LayerID from, const VertexConstraint& vertex)
{
  using Table = LinkMaterial<Det>;
  const double rMin = radiusMin<Det>(from);
  if (!(rMin > 0.) || !(params.minPtGeV > 0.)) {
    return std::numeric_limits<double>::infinity();
  }
  const auto logBeta2 = static_cast<Scalar>(logBeta2Min(params));
  const Interval reach = vertexReach(params, vertex);
  double largest = 0.;
  for (std::int32_t s = 0; s < Table::kSecond; ++s) {
    double node = 0.;
    for (std::int32_t phi = 0; phi < Table::kNodePhi; ++phi) {
      const auto& n = material.nodes[(((std::uint32_t{from} * Table::kNodePhi) + static_cast<std::uint32_t>(phi)) * Table::kSecond) + static_cast<std::uint32_t>(s)];
      node += static_cast<double>(n.x) / Table::kNodePhi;
    }
    for (std::int32_t v = 0; v < Table::kVertexZ; ++v) {
      const double incidence = largestIncidence<Det>(material, from, s, v, reach);
      if (!(incidence > 0.)) {
        continue;
      }
      double gapX = 0.;
      double gapM2 = 0.;
      for (std::int32_t phi = 0; phi < Table::kGapPhi; ++phi) {
        const auto& gap = material.vertexGaps[(((((std::uint32_t{from} * Table::kGapPhi) + static_cast<std::uint32_t>(phi)) * Table::kVertexZ) + static_cast<std::uint32_t>(v)) * Table::kSecond) + static_cast<std::uint32_t>(s)];
        gapX += static_cast<double>(gap.x) * incidence / Table::kGapPhi;
        gapM2 += static_cast<double>(gap.x) * incidence * static_cast<double>(gap.m2) / Table::kGapPhi;
      }
      const double nodeX = node * incidence;
      const double spanX = gapX + nodeX;
      if (spanX > 0.) {
        const auto logSpan = static_cast<Scalar>(std::log(spanX));
        // K is linear in x at one logarithm: the gap's x <f^2> and the node's x.
        largest = std::max(largest, static_cast<double>(materialKick(static_cast<Scalar>(gapM2 + nodeX), logSpan, logBeta2)));
      }
    }
  }
  const double beam = vertex.sigmaTransverse / rMin;
  return std::sqrt((largest / (params.minPtGeV * params.minPtGeV) * betaMomentumFactor(params)) + (beam * beam));
}

/// How a link's bound takes its gaps' azimuth bins: the largest, a bound
/// every pair keeps (`IterationParameters::betweenBound`, the test before the
/// pair's own term), or the mean, what a window holds (as the per-layer
/// averages did).
enum class AzimuthBound : std::uint8_t { kLargest,
                                         kMean };

/// A bound on `sqrt(T)` [cm] (`betweenScattering`) for the pairs on the link
/// \a from -> \a to: the largest `K(x) (1 - 2 <f> + <f^2>) L^2` over the
/// link's source and target bins a pair can join, each with the longest chord
/// between them (their corners on the layers' nominal surfaces), at pT_min
/// (`g >= 1`), its azimuth bins taken by \a azimuth. The pairs leave from the
/// link's source interval toward the line from the vertex reach through the
/// source, a bin either side. Zero if its gaps hold nothing.
template <DetectorLayout Det>
[[nodiscard]] inline double betweenKinkBound(const LinkMaterial<Det>& material, const SearchWindowParams& params, LayerID from, LayerID to, const VertexConstraint& vertex, AzimuthBound azimuth)
{
  using Table = LinkMaterial<Det>;
  const std::int32_t link = material.link(from, to);
  if (link < 0 || !(params.minPtGeV > 0.)) {
    return 0.;
  }
  const auto logBeta2 = static_cast<Scalar>(logBeta2Min(params));
  // The pair's logarithm at its span, the source's at its largest (bounded) or mean.
  const double sourceSpan = largestSourceSpan(material, params, from, vertex);
  const bool barrelFrom = Detector<Det>::isBarrelLayer(from);
  const bool barrelTo = Detector<Det>::isBarrelLayer(to);
  const double positionFrom = PlacedEnvelopes<Det>::get(from).position;
  const double positionTo = PlacedEnvelopes<Det>::get(to).position;
  const auto& sourceAxis = material.linkSource[static_cast<std::uint32_t>(link)];
  const auto& targetAxis = material.layerSecond[to];
  // A layer point at second coordinate c, as (r, z).
  const auto rz = [](bool barrel, double position, double c) { return barrel ? std::pair{position, c} : std::pair{c, position}; };
  const auto& pair = Det::layerLinks[static_cast<std::uint32_t>(link)];
  const Interval reach = vertexReach(params, vertex);
  // The target coordinate the line from a vertex at zv through a source point reaches.
  const auto projected = [&](double cs, double zv) {
    const auto [ra, za] = rz(barrelFrom, positionFrom, cs);
    if (barrelTo) {
      return ra > 0. ? zv + ((za - zv) * positionTo / ra) : std::numeric_limits<double>::quiet_NaN();
    }
    return std::abs(za - zv) > 0. ? ra * (positionTo - zv) / (za - zv) : std::numeric_limits<double>::quiet_NaN();
  };
  double largest = 0.;
  for (std::int32_t s = 0; s < Table::kSourceSecond; ++s) {
    const double s0 = sourceAxis.low + (static_cast<float>(s) / sourceAxis.perUnit);
    const double s1 = sourceAxis.low + (static_cast<float>(s + 1) / sourceAxis.perUnit);
    if (s1 < pair.sourceMin || s0 > pair.sourceMax) {
      continue;
    }
    double low = std::numeric_limits<double>::infinity();
    double high = -std::numeric_limits<double>::infinity();
    for (const double cs : {s0, s1}) {
      for (const double zv : {reach.low, reach.high}) {
        const double c = projected(cs, zv);
        if (std::isfinite(c)) {
          low = std::min(low, c);
          high = std::max(high, c);
        }
      }
    }
    if (!(low <= high)) {
      continue;
    }
    const std::int32_t tLow = std::max(0, targetAxis.bin(static_cast<Scalar>(low), Table::kSecond) - 1);
    const std::int32_t tHigh = std::min(Table::kSecond - 1, targetAxis.bin(static_cast<Scalar>(high), Table::kSecond) + 1);
    for (std::int32_t t = tLow; t <= tHigh; ++t) {
      const double t0 = targetAxis.low + (static_cast<float>(t) / targetAxis.perUnit);
      const double t1 = targetAxis.low + (static_cast<float>(t + 1) / targetAxis.perUnit);
      double chord2 = 0.;
      for (const double cs : {s0, s1}) {
        for (const double ct : {t0, t1}) {
          const auto [ra, za] = rz(barrelFrom, positionFrom, cs);
          const auto [rb, zb] = rz(barrelTo, positionTo, ct);
          chord2 = std::max(chord2, ((rb - ra) * (rb - ra)) + ((zb - za) * (zb - za)));
        }
      }
      double cell = 0.;
      for (std::int32_t phi = 0; phi < Table::kGapPhi; ++phi) {
        const auto& gap = material.gaps[(((((static_cast<std::uint32_t>(link) * Table::kGapPhi) + static_cast<std::uint32_t>(phi)) * Table::kSourceSecond) + static_cast<std::uint32_t>(s)) * Table::kSecond) + static_cast<std::uint32_t>(t)];
        const auto logSpan = static_cast<Scalar>(std::log(sourceSpan + static_cast<double>(gap.x) > 0. ? sourceSpan + static_cast<double>(gap.x) : 1.));
        const double term = static_cast<double>(materialKick(gap.x, logSpan, logBeta2)) * (1. - (2. * static_cast<double>(gap.m1)) + static_cast<double>(gap.m2)) * chord2;
        cell = azimuth == AzimuthBound::kLargest ? std::max(cell, term) : cell + (term / Table::kGapPhi);
      }
      largest = std::max(largest, cell);
    }
  }
  return std::sqrt(largest * betaMomentumFactor(params)) / params.minPtGeV;
}

/// The window a pair query uses, and the one that sizes the grid's cells
/// (GridSizing.h): its phi and second-coordinate half-widths hold what the
/// pair tests accept with \a material's scattering.
template <DetectorLayout Det>
[[nodiscard]] inline LayerPairWindow layerPairWindow(LayerID from, LayerID to, const SearchWindowParams& params, const VertexConstraint& vertex, const LinkMaterial<Det>& material)
{
  return layerPairWindow<Det>(from, to, secondMax<Det>(from), params, vertex, sourceKinkBound<Det>(material, params, from, vertex), betweenKinkBound<Det>(material, params, from, to, vertex, AzimuthBound::kMean));
}

/// The window an iteration seeding from \a source queries with, and sizes the
/// grid's cells with: without an origin any straight leg's
/// (`originlessPairWindow`), else the beam line's (`layerPairWindow`).
template <DetectorLayout Det>
[[nodiscard]] inline LayerPairWindow iterationPairWindow(LayerID from, LayerID to, const SearchWindowParams& params, VertexSource source, const VertexConstraint& vertex, const LinkMaterial<Det>& material)
{
  return source == VertexSource::kNone ? originlessPairWindow<Det>(from, to) : layerPairWindow<Det>(from, to, params, vertex, material);
}

/// The two tests' ingredients for one source, vertex and candidate.
struct PairTerms {
  bool beforeTurn = false;      ///< the circle's turning point is not between the two
  Scalar curvature = 0.F;       ///< 2 |P_i x P_j|
  Scalar bound = 0.F;           ///< kappa_max r_i r_j d: the curvature bound, in the same units
  Scalar varCurvature = 0.F;    ///< (2 sigma_X)^2
  Scalar collinearity2 = 0.F;   ///< D^2
  Scalar varCollinearity = 0.F; ///< sigma_D^2
};

namespace detail
{

struct ArcLengthDerivatives {
  Scalar value;
  Scalar chord;
  Scalar curvature;
};

/// Arc length and its derivatives with respect to chord and signed curvature.
[[nodiscard]] GPUhdi() GPUAlwaysInline() ArcLengthDerivatives arcLengthDerivatives(Scalar curvature, Scalar chord)
{
  const Scalar halfTurnSine = o2::gpu::GPUCommonMath::Clamp(0.5F * curvature * chord, -1.F, 1.F);
  const auto quotient = asinOverArgumentWithDerivative(halfTurnSine);
  const Scalar ratio = quotient.value;
  const Scalar derivative = quotient.derivative;
  return {.value = chord * ratio, .chord = ratio + (halfTurnSine * derivative), .curvature = 0.5F * chord * chord * derivative};
}

struct ArcPairGeometry {
  Scalar curvature;
  Scalar curvatureScale;
  Scalar chord;
  bool bounded; ///< the curvature is the cut's bound, which the positions do not move
  ArcLengthDerivatives source;
  ArcLengthDerivatives target;
};

/// The circle through the beam point and the pair, with its curvature held to
/// what a track above pT_min can have. The positions can put it further out
/// -- up to 2/r when the chord is unresolved -- and the arc length there, and
/// its derivatives, would describe a track the cut does not admit. The bound
/// never exceeds a curvature the measured one already reached, so the arcsine
/// argument stays within [-1, 1]. Arc length is even in the curvature, so the
/// bound's sign does not matter.
[[nodiscard]] GPUhdi() GPUAlwaysInline() ArcPairGeometry arcPairGeometry(const PairPoint& source, const PairPoint& target, Scalar cross, Scalar chord, Scalar maxCurvature)
{
  const Scalar curvatureScale = 2.F / (source.r * target.r * chord);
  const Scalar measured = curvatureScale * cross;
  // Also false for the NaN of a zero chord.
  const bool bounded = !(o2::gpu::GPUCommonMath::Abs(measured) <= maxCurvature);
  const Scalar curvature = bounded ? maxCurvature : measured;
  return {
    .curvature = curvature,
    .curvatureScale = curvatureScale,
    .chord = chord,
    .bounded = bounded,
    .source = arcLengthDerivatives(curvature, source.r),
    .target = arcLengthDerivatives(curvature, target.r),
  };
}

} // namespace detail

/// What the pair tests compute of the two positions alone: the transverse
/// chord and the cross product, once a pair for all of them.
struct PairChord {
  Scalar dx;
  Scalar dy;
  Scalar chord2;
  Scalar chord;
  Scalar cross;
};

[[nodiscard]] GPUhdi() GPUAlwaysInline() PairChord pairChord(const PairPoint& source, const PairPoint& target)
{
  const Scalar dx = target.x - source.x;
  const Scalar dy = target.y - source.y;
  const Scalar chord2 = (dx * dx) + (dy * dy);
  return {.dx = dx, .dy = dy, .chord2 = chord2, .chord = o2::gpu::GPUCommonMath::Sqrt(chord2), .cross = (source.x * target.y) - (target.x * source.y)};
}

/// The circle through the beam point turns where it is farthest from it,
/// opposite the beam point across its centre. The pair lies on one side of
/// that diameter -- no turning point between them -- exactly when
/// `P_i . P_j` lies between `r_i^2` and `r_j^2`: the step from one to the
/// other leads outwards at both ends. That also keeps them in one hemisphere.
/// A straight line from anywhere (a pass without an origin) the same: it
/// leads outward from the one to the other.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool noTurnBetween(const PairPoint& source, const PairPoint& target)
{
  const Scalar dot = (source.x * target.x) + (source.y * target.y);
  const Scalar sourceR2 = source.r * source.r;
  const Scalar targetR2 = target.r * target.r;
  const Scalar lowR2 = sourceR2 < targetR2 ? sourceR2 : targetR2;
  const Scalar highR2 = sourceR2 < targetR2 ? targetR2 : sourceR2;
  return dot >= lowR2 && dot <= highR2;
}

/// \a betweenScatter is `betweenScattering` for the pair, zero where its
/// link's gap holds no material.
///
/// Always inlined, as `pairCompatible`: the tracklet stage's grid walk runs
/// the transverse test on every candidate, and GCC keeps the two out of
/// line once the walk grows past its inlining limits -- a call each per
/// candidate.
[[nodiscard]] GPUhdi() GPUAlwaysInline() PairTerms transversePairTerms(const PairCut& cut, const VertexSeed& vertex, const PairPoint& source, Scalar sourceScatter, const PairPoint& target, const PairChord& pair, Scalar betweenScatter = 0.F)
{
  const Scalar chord2 = pair.chord2;
  const Scalar chord = pair.chord;
  const Scalar cross = pair.cross;
  // A kink between the two moves P_j against the side of length r_i.
  const Scalar varCross = (source.r * source.r * (target.varRPhi + betweenScatter)) + (target.r * target.r * source.varRPhi) + (chord2 * (vertex.varTransverse + sourceScatter));

  return {
    .beforeTurn = noTurnBetween(source, target),
    .curvature = 2.F * o2::gpu::GPUCommonMath::Abs(cross),
    .bound = cut.maxCurvature * source.r * target.r * chord,
    .varCurvature = 4.F * varCross,
  };
}

GPUhdi() void addLongitudinalPairTerms(PairTerms& terms, const PairCut& cut, const VertexSeed& vertex, const PairPoint& source, Scalar sourceScatter, const PairPoint& target, const PairChord& pair, Scalar betweenScatter = 0.F)
{
  const auto geometry = detail::arcPairGeometry(source, target, pair.cross, pair.chord, cut.maxCurvature);
  const Scalar dzSource = source.z - vertex.z;
  const Scalar dzTarget = target.z - vertex.z;
  const Scalar collinearity = (geometry.source.value * dzTarget) - (geometry.target.value * dzSource);

  const Scalar invSourceR = 1.F / source.r;
  const Scalar invTargetR = 1.F / target.r;
  // The collinearity's gradient in each point's radial and r-phi directions. At a fixed curvature only
  // the arcs' radii move it, radially.
  Scalar sourceRadialD = dzTarget * geometry.source.chord;
  Scalar targetRadialD = -dzSource * geometry.target.chord;
  Scalar sourceTangentD = 0.F;
  Scalar targetTangentD = 0.F;
  // A measured curvature follows the positions: kappa = 2 sin(h) / c, h the pair's azimuth gap and c its
  // chord. A point moved radially changes the chord alone, d kappa / d r_i = kappa (P_i . P_j - r_i^2) /
  // (r_i c^2); moved along r-phi by a length t, h and the chord, d kappa / d t_i = 2 ((P_i x P_j)^2 / c^2 -
  // P_i . P_j) / (r_i^2 r_j c), and the opposite sign for the target. A bounded curvature does not follow.
  if (!geometry.bounded) {
    const Scalar curvatureD = (dzTarget * geometry.source.curvature) - (dzSource * geometry.target.curvature);
    const Scalar dot = (source.x * target.x) + (source.y * target.y);
    const Scalar invChord2 = 1.F / pair.chord2;
    const Scalar radialScale = curvatureD * geometry.curvature * invChord2;
    sourceRadialD += radialScale * (dot - (source.r * source.r)) * invSourceR;
    targetRadialD += radialScale * (dot - (target.r * target.r)) * invTargetR;
    const Scalar turnD = curvatureD * geometry.curvatureScale * ((pair.cross * pair.cross * invChord2) - dot);
    sourceTangentD = turnD * invSourceR;
    targetTangentD = -turnD * invTargetR;
  }
  // The vertex's gradient is minus the sum of the two (D depends on the positions relative to it).
  const Scalar sourceRadialX = sourceRadialD * invSourceR;
  const Scalar sourceTangentX = sourceTangentD * invSourceR;
  const Scalar targetRadialX = targetRadialD * invTargetR;
  const Scalar targetTangentX = targetTangentD * invTargetR;
  const Scalar vertexXD = -((sourceRadialX * source.x) - (sourceTangentX * source.y) + (targetRadialX * target.x) - (targetTangentX * target.y));
  const Scalar vertexYD = -((sourceRadialX * source.y) + (sourceTangentX * source.x) + (targetRadialX * target.y) + (targetTangentX * target.x));
  const Scalar arcGap = geometry.target.value - geometry.source.value;
  const Scalar dz = target.z - source.z;

  terms.collinearity2 = collinearity * collinearity;
  terms.varCollinearity = (source.varR * o2::gpu::CAMath::Square(sourceRadialD)) + (o2::gpu::CAMath::Square(sourceTangentD) * source.varRPhi) + (o2::gpu::CAMath::Square(targetRadialD) * target.varR) + (o2::gpu::CAMath::Square(targetTangentD) * target.varRPhi) + ((o2::gpu::CAMath::Square(vertexXD) + o2::gpu::CAMath::Square(vertexYD)) * vertex.varTransverse) + (o2::gpu::CAMath::Square(geometry.target.value) * source.varZ) + (o2::gpu::CAMath::Square(geometry.source.value) * target.varZ) + (o2::gpu::CAMath::Square(arcGap) * vertex.varZ) + ((o2::gpu::CAMath::Square(arcGap) + o2::gpu::CAMath::Square(dz)) * sourceScatter) +
                          // A kink between the two moves the target's (s, z) corner against (s_i, z_i - z_v).
                          ((o2::gpu::CAMath::Square(geometry.source.value) + o2::gpu::CAMath::Square(dzSource)) * betweenScatter);
}

/// Curvature at most the iteration's bound plus n sigma, with no turning point
/// between the two: a looper's legs are found separately. The bound is
/// `maxCurvature`: without a pT cut infinite, and every such pair passes; at
/// zero field 0, every track straight, and the pair points back to the beam
/// line within its errors and scattering.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool transverseCompatible(const PairCut& cut, const PairTerms& terms)
{
  const Scalar excess = terms.curvature - terms.bound;
  return terms.beforeTurn && (!(excess > 0.F) || excess * excess <= cut.nSigma * cut.nSigma * terms.varCurvature);
}

/// Collinear with the vertex within n sigma.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool longitudinalCompatible(const PairCut& cut, const PairTerms& terms) { return terms.collinearity2 <= cut.nSigma * cut.nSigma * terms.varCollinearity; }

/// Always inlined (see `transversePairTerms`); the longitudinal half stays a
/// call.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool pairCompatible(const PairCut& cut, const VertexSeed& vertex, const PairPoint& source, Scalar sourceScatter, const PairPoint& target, const PairChord& pair, Scalar betweenScatter = 0.F)
{
  countWork(WorkCounter::kPairTests);
  auto terms = transversePairTerms(cut, vertex, source, sourceScatter, target, pair, betweenScatter);
  if (!transverseCompatible(cut, terms)) {
    return false;
  }
  countWork(WorkCounter::kPairTransverse);
  addLongitudinalPairTerms(terms, cut, vertex, source, sourceScatter, target, pair, betweenScatter);
  const bool compatible = longitudinalCompatible(cut, terms);
  if (compatible) {
    countWork(WorkCounter::kPairLongitudinal);
  }
  return compatible;
}

/// The same three without a computed chord, one pair each: the harness's and
/// the tests' (StageFunnel, the sink).
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool pairCompatible(const PairCut& cut, const VertexSeed& vertex, const PairPoint& source, Scalar sourceScatter, const PairPoint& target, Scalar betweenScatter = 0.F) { return pairCompatible(cut, vertex, source, sourceScatter, target, pairChord(source, target), betweenScatter); }

[[nodiscard]] GPUhdi() GPUAlwaysInline() PairTerms transversePairTerms(const PairCut& cut, const VertexSeed& vertex, const PairPoint& source, Scalar sourceScatter, const PairPoint& target, Scalar betweenScatter = 0.F) { return transversePairTerms(cut, vertex, source, sourceScatter, target, pairChord(source, target), betweenScatter); }

GPUhdi() GPUAlwaysInline() void addLongitudinalPairTerms(PairTerms& terms, const PairCut& cut, const VertexSeed& vertex, const PairPoint& source, Scalar sourceScatter, const PairPoint& target, Scalar betweenScatter = 0.F) { addLongitudinalPairTerms(terms, cut, vertex, source, sourceScatter, target, pairChord(source, target), betweenScatter); }

} // namespace o2::gtrk

#endif // GTRK_TRACKING_PAIRCUT_H_
