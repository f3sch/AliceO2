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

/// \file LinkMaterial.h
/// \brief The material the seeding's cuts take, as layer nodes and the gaps
/// between them.
///
/// Space is cut once into parts, each holding its own material: every
/// layer's envelope band (a node, at the layer's crossing), and between the
/// band a link leaves and the band it reaches, the link's gap. A link that
/// passes over a layer has that layer's band inside its gap, at its own
/// place, so a pair on the direct link and the same track through the layer
/// between cross the same material. Nothing is compiled in: the parts follow
/// from the placed layers' envelopes, their material from the LUT
/// (LinkMaterialFill.h).
///
/// - A node holds its band's x/X0 at normal incidence, in front of and behind
///   the layer's sensors, by azimuth and the layer's second coordinate.
/// - A link's gap holds the x/X0 along the chord between the two layers'
///   nominal surfaces from the source band's outer face to the target band's
///   inner face, with the mean and the mean square of the fraction of the
///   chord the material lies at, by the source's azimuth and both ends'
///   second coordinates: what the lever arms of a pair window or a triplet
///   kink need.
/// - A layer's vertex gap holds the same from a vertex on the beam line to
///   the layer's inner face, by azimuth, the vertex's z and the layer's
///   second coordinate.

#ifndef GTRK_TRACKING_LINKMATERIAL_H_
#define GTRK_TRACKING_LINKMATERIAL_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/Definitions.h"

#include "GPUCommonMath.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace o2::gtrk
{

/// A gap's material: x/X0 along its chord and the x-weighted mean and mean
/// square of the fraction of the chord, from the source, it lies at.
struct MaterialGap {
  float x = 0.F;
  float m1 = 0.F;
  float m2 = 0.F;
  float half = 0.F; ///< the share of x in the chord's first half, from the source
  float xRho = 0.F; ///< rho x along the chord [g/cm^2], for the energy loss
};

/// The parts a vertex gap's profile resolves its ray in.
inline constexpr int kGapProfile = 16;

/// How a vertex gap's material lies along its ray from the vertex to the
/// layer's nominal point: the share before each of `kGapProfile + 1` equal
/// fractions of the ray, in 1/255. A segment from a point inside the layer's
/// reach (a decay vertex) takes the share beyond it.
struct GapProfile {
  std::uint8_t share[kGapProfile + 1]{};
  /// The share of the gap's x between the fractions \a f0 < \a f1 of the
  /// ray, and the x-weighted mean and mean square of the fraction over it
  /// (\a mean, \a square), each part's material even across it.
  GPUhdi() Scalar between(Scalar f0, Scalar f1, Scalar& mean, Scalar& square) const
  {
    Scalar total = 0.F;
    Scalar first = 0.F;
    Scalar second = 0.F;
    for (int k = 0; k < kGapProfile; ++k) {
      const Scalar low = static_cast<Scalar>(k) / static_cast<Scalar>(kGapProfile);
      const Scalar high = static_cast<Scalar>(k + 1) / static_cast<Scalar>(kGapProfile);
      const Scalar a = f0 > low ? f0 : low;
      const Scalar b = f1 < high ? f1 : high;
      if (!(b > a)) {
        continue;
      }
      const Scalar w = (static_cast<Scalar>(share[k + 1]) - static_cast<Scalar>(share[k])) * (1.F / 255.F) * (b - a) * static_cast<Scalar>(kGapProfile);
      total += w;
      first += w * 0.5F * (a + b);
      second += w * ((a * a) + (a * b) + (b * b)) * (1.F / 3.F);
    }
    mean = total > 0.F ? first / total : 0.5F * (f0 + f1);
    square = total > 0.F ? second / total : mean * mean;
    return total;
  }
};

/// The depths a node's profile resolves its band in.
inline constexpr int kNodeDepths = 32;

/// A layer band's x/X0 at normal incidence and how it is spread across the
/// band's depth: the share in front of each of `kNodeDepths + 1` equal depths
/// from the band's inner face to its outer one, in 1/255. A measurement's
/// depth in its band splits the band in what lies in front of it and behind
/// it, whichever of a layer's staggered or tilted staves it is on.
struct MaterialNode {
  float x = 0.F;    ///< across the band
  float xRho = 0.F; ///< rho x across the band at normal incidence [g/cm^2]
  std::uint8_t depth[kNodeDepths + 1]{};
  /// The share of the band in front of the fraction \a t of its depth.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar frontShare(Scalar t) const
  {
    const Scalar at = (t < 0.F ? 0.F : (t > 1.F ? 1.F : t)) * static_cast<Scalar>(kNodeDepths);
    const int k = at >= static_cast<Scalar>(kNodeDepths) ? kNodeDepths - 1 : static_cast<int>(at);
    const auto low = static_cast<Scalar>(depth[k]);
    return (low + ((static_cast<Scalar>(depth[k + 1]) - low) * (at - static_cast<Scalar>(k)))) * (1.F / 255.F);
  }
};

/// One step's material between two linked layers' crossings, in the step's
/// direction (`LinkMaterial::step`): the band it leaves (the part on its far
/// side from where the track came) and the gap as one kick at the mean place
/// with a spread, the band it reaches (the part in front of its sensor) at
/// the end.
struct LinkStep {
  Scalar before = 0.F;       ///< x/X0 along the step: the band left and the gap
  Scalar after = 0.F;        ///< x/X0 along the step: the band reached, before its sensor
  Scalar lever = 1.F;        ///< 1 - <f> of `before`, f the fraction of the step
  Scalar spread = 0.F;       ///< <f^2> - <f>^2 of `before`
  Scalar xRho = 0.F;         ///< the whole step's [g/cm^2]
  Scalar left = 0.F;         ///< of `before`, the band left
  Scalar gapFirstHalf = 0.F; ///< of `before` less `left`, what lies in the step's first half
  bool found = false;        ///< the layers are a link (either way)
};

/// One binned coordinate: bin = (value - low) * perUnit, clamped.
struct MaterialAxis {
  float low = 0.F;
  float perUnit = 0.F;
  [[nodiscard]] GPUhdi() GPUAlwaysInline() int bin(Scalar value, int bins) const
  {
    const int b = static_cast<int>((value - low) * perUnit);
    return b < 0 ? 0 : (b >= bins ? bins - 1 : b);
  }
};

template <DetectorLayout Det>
struct LinkMaterial {
  static constexpr int kLayers = Detector<Det>::nLayers;
  static constexpr int kLinks = static_cast<int>(std::size(Det::layerLinks));
  /// Azimuth bins of the nodes: as fine as a stave, whose sensors split its
  /// band in front and behind differently across its width.
  static constexpr int kNodePhi = 1024;
  /// Azimuth bins of the gaps: supports and services are not uniform in it.
  static constexpr int kGapPhi = 32;
  /// Bins of a layer's second coordinate, of a link's source layer's (finer:
  /// a link's material changes fast with where it leaves, past a stave's
  /// end), and the vertex's z.
  static constexpr int kSecond = 16;
  static constexpr int kSourceSecond = 32;
  static constexpr int kVertexZ = 8;

  [[nodiscard]] static constexpr std::array<std::int16_t, static_cast<std::size_t>(kLayers) * kLayers> noLinks()
  {
    std::array<std::int16_t, static_cast<std::size_t>(kLayers) * kLayers> none{};
    for (auto& link : none) {
      link = -1;
    }
    return none;
  }

  /// The link from layer `from` to layer `to`, -1 for none (an unfilled
  /// table has none); each link's layers.
  std::array<std::int16_t, static_cast<std::size_t>(kLayers) * kLayers> linkOf = noLinks();
  std::array<LayerID, kLinks> linkFrom{};
  std::array<LayerID, kLinks> linkTo{};
  /// Each layer's second coordinate over its envelope, the source interval of
  /// each link, the vertex's z.
  std::array<MaterialAxis, kLayers> layerSecond{};
  std::array<MaterialAxis, kLinks> linkSource{};
  MaterialAxis vertexZ{};
  /// Each layer's band in depth: the fraction of the way from its inner face
  /// to its outer one, in radius on a barrel layer, in |z| on a disc.
  std::array<MaterialAxis, kLayers> layerDepth{};
  /// Each disc's side of the vertex, +1 or -1 (0 on a barrel layer).
  std::array<std::int8_t, kLayers> layerSide{};
  std::array<MaterialNode, static_cast<std::size_t>(kLayers) * kNodePhi * kSecond> nodes{};
  std::array<MaterialGap, static_cast<std::size_t>(kLinks) * kGapPhi * kSourceSecond * kSecond> gaps{};
  /// Per row of a link's gaps (one source azimuth and coordinate), the
  /// largest `x (1 - 2 m1 + m2)` over the targets: what a triplet's dip-angle
  /// window takes for its second link before the third point is known.
  std::array<float, static_cast<std::size_t>(kLinks) * kGapPhi * kSourceSecond> gapRowMax{};
  /// Per link, its gaps' largest x/X0: what bounds a pair's span before its gap is read.
  std::array<float, kLinks> gapMax{};
  /// Each layer's vertex gaps, from the vertex to its band, at normal
  /// incidence (radial on a barrel layer, along z on a disc): a ray's x/X0
  /// times its depth over its length. A use takes its own line's path.
  std::array<MaterialGap, static_cast<std::size_t>(kLayers) * kGapPhi * kVertexZ * kSecond> vertexGaps{};
  /// Each vertex gap's profile along its ray (`segment`).
  std::array<GapProfile, static_cast<std::size_t>(kLayers) * kGapPhi * kVertexZ * kSecond> vertexProfiles{};

  /// The azimuth bin of \a phi [rad, either -pi..pi or 0..2pi].
  [[nodiscard]] GPUhdi() GPUAlwaysInline() static int phiBinOf(Scalar phi, int bins)
  {
    phi = phi < 0.F ? phi + kTwoPi<Scalar> : phi;
    const int b = static_cast<int>(phi * (static_cast<Scalar>(bins) / kTwoPi<Scalar>));
    return b < 0 ? 0 : (b >= bins ? bins - 1 : b);
  }
  [[nodiscard]] GPUhdi() GPUAlwaysInline() static int phiBin(Scalar x, Scalar y, int bins) { return phiBinOf(atan2Of(y, x), bins); }

  [[nodiscard]] GPUhdi() GPUAlwaysInline() int link(LayerID from, LayerID to) const { return linkOf[(static_cast<std::size_t>(from) * kLayers) + to]; }

  /// Layer \a layer's node at azimuth \a phi with second coordinate \a second.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() const MaterialNode& nodeAt(LayerID layer, Scalar phi, Scalar second) const
  {
    const auto bin = static_cast<std::size_t>(phiBinOf(phi, kNodePhi));
    const auto s = static_cast<std::size_t>(layerSecond[layer].bin(second, kSecond));
    return nodes[(((static_cast<std::size_t>(layer) * kNodePhi) + bin) * kSecond) + s];
  }

  /// Layer \a layer's node at the point (\a x, \a y) with second coordinate \a second.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() const MaterialNode& node(LayerID layer, Scalar x, Scalar y, Scalar second) const
  {
    const auto phi = static_cast<std::size_t>(phiBin(x, y, kNodePhi));
    const auto s = static_cast<std::size_t>(layerSecond[layer].bin(second, kSecond));
    return nodes[(((static_cast<std::size_t>(layer) * kNodePhi) + phi) * kSecond) + s];
  }

  /// The first cell of link \a link's row for a source at (\a x, \a y) with
  /// second coordinate \a coordinate.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint32_t gapRow(std::size_t link, Scalar x, Scalar y, Scalar coordinate) const
  {
    const auto phi = static_cast<std::uint32_t>(phiBin(x, y, kGapPhi));
    const auto s = static_cast<std::uint32_t>(linkSource[link].bin(coordinate, kSourceSecond));
    return ((((static_cast<std::uint32_t>(link) * kGapPhi) + phi) * kSourceSecond) + s) * kSecond;
  }

  /// The same for a source in gap azimuth bin \a phiBin (`phiBin` with `kGapPhi`).
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint32_t gapRowAt(std::size_t link, std::uint32_t phiBin, Scalar coordinate) const
  {
    const auto s = static_cast<std::uint32_t>(linkSource[link].bin(coordinate, kSourceSecond));
    return ((((static_cast<std::uint32_t>(link) * kGapPhi) + phiBin) * kSourceSecond) + s) * kSecond;
  }

  /// The fraction of layer \a layer's band's depth a point at \a r, \a z lies at.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar depthOf(LayerID layer, Scalar r, Scalar z) const
  {
    const auto& axis = layerDepth[layer];
    return ((Detector<Det>::isBarrelLayer(layer) ? r : o2::gpu::GPUCommonMath::Abs(z)) - axis.low) * axis.perUnit;
  }

  /// A point's second coordinate on layer \a layer: z on a barrel layer, r on a disc.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() static Scalar second(LayerID layer, Scalar r, Scalar z) { return Detector<Det>::isBarrelLayer(layer) ? z : r; }

  /// The gap of row \a row (`gapRow`) toward a target on layer \a target
  /// with second coordinate \a second.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() const MaterialGap& gap(std::uint32_t row, LayerID target, Scalar second) const { return gaps[row + static_cast<std::uint32_t>(layerSecond[target].bin(second, kSecond))]; }

  [[nodiscard]] GPUhdi() GPUAlwaysInline() float rowMax(std::uint32_t row) const { return gapRowMax[row / kSecond]; }

  /// The step from \a start on layer \a fromLayer to \a end on \a toLayer:
  /// along a link outward, against one inward; between two barrel layers that
  /// are no link, the run of links over the barrel layers the chord crosses
  /// between them, which the track crossed without a measurement
  /// (`composedStep`); `found` false where neither is. A node takes the
  /// step's incidence on its layer, radial on a barrel layer, along z on a
  /// disc; a barrel-to-barrel gap the chord's obliquity in the transverse
  /// plane over the table's radial chord.
  [[nodiscard]] GPUhdi() LinkStep step(LayerID fromLayer, LayerID toLayer, const Scalar start[3], const Scalar end[3]) const
  {
    if (fromLayer >= static_cast<LayerID>(kLayers) || toLayer >= static_cast<LayerID>(kLayers)) {
      return {};
    }
    if (link(fromLayer, toLayer) >= 0 || link(toLayer, fromLayer) >= 0) {
      return linkStep(fromLayer, toLayer, start, end);
    }
    return composedStep(fromLayer, toLayer, start, end);
  }

  /// `step` between two barrel layers that are no link: the chord from
  /// \a start to \a end in parts at the barrel layers between them whose
  /// band's middle radius it crosses within their z (layers numbered by
  /// radius), each part a link's step (`linkStep`) and each crossed layer's
  /// band whole at its crossing; one step of them, its `before` everything up
  /// to the band reached, at its parts' places. `found` false where a part is
  /// no link.
  [[nodiscard]] GPUhdi() LinkStep composedStep(LayerID fromLayer, LayerID toLayer, const Scalar start[3], const Scalar end[3]) const
  {
    using Math = o2::gpu::GPUCommonMath;
    LinkStep out;
    if (!Detector<Det>::isBarrelLayer(fromLayer) || !Detector<Det>::isBarrelLayer(toLayer) || fromLayer == toLayer) {
      return out;
    }
    const Scalar d[3] = {end[0] - start[0], end[1] - start[1], end[2] - start[2]};
    const Scalar a = (d[0] * d[0]) + (d[1] * d[1]);
    const Scalar b = (start[0] * d[0]) + (start[1] * d[1]);
    const Scalar c0 = (start[0] * start[0]) + (start[1] * start[1]);
    if (!(a > 0.F)) {
      return out;
    }
    const int sign = toLayer > fromLayer ? 1 : -1;
    LayerID previous = fromLayer;
    Scalar previousPoint[3] = {start[0], start[1], start[2]};
    Scalar previousT = 0.F;
    // The sums of the parts before the band reached: x, x <f>, x <f^2>.
    Scalar sum = 0.F;
    Scalar first = 0.F;
    Scalar second = 0.F;
    // What lies in the step's first half (the parts' own split at their
    // middles where the step's middle falls inside one, linear elsewhere).
    Scalar firstHalf = 0.F;
    const auto add = [&](const LinkStep& part, Scalar t0, Scalar t1, bool last) {
      const Scalar span = t1 - t0;
      const Scalar mean = 1.F - part.lever;
      const Scalar square = part.spread + (mean * mean);
      sum += part.before;
      first += part.before * (t0 + (span * mean));
      second += part.before * ((t0 * t0) + (2.F * t0 * span * mean) + (span * span * square));
      const Scalar middle = part.before - part.left;
      if (t0 > 0.F && t0 < 0.5F) {
        firstHalf += part.left;
      }
      if (t1 <= 0.5F) {
        firstHalf += middle;
      } else if (t0 < 0.5F) {
        firstHalf += middle * (0.5F - t0) / span;
      }
      if (!last) {
        sum += part.after;
        first += part.after * t1;
        second += part.after * t1 * t1;
        if (t1 < 0.5F) {
          firstHalf += part.after;
        }
      }
      out.xRho += part.xRho;
    };
    for (int layer = static_cast<int>(fromLayer) + sign; std::cmp_not_equal(layer, toLayer); layer += sign) {
      const auto crossed = static_cast<LayerID>(layer);
      if (!Detector<Det>::isBarrelLayer(crossed)) {
        return out;
      }
      // The chord's crossing of the band's middle radius beyond the last part.
      const Scalar radius = layerDepth[crossed].low + (0.5F / layerDepth[crossed].perUnit);
      const Scalar discriminant = (b * b) - (a * (c0 - (radius * radius)));
      if (!(discriminant >= 0.F)) {
        continue;
      }
      const Scalar root = Math::Sqrt(discriminant);
      Scalar t = (-b - root) / a;
      if (!(t > previousT)) {
        t = (-b + root) / a;
      }
      if (!(t > previousT) || !(t < 1.F)) {
        continue;
      }
      const Scalar point[3] = {start[0] + (t * d[0]), start[1] + (t * d[1]), start[2] + (t * d[2])};
      const auto& z = layerSecond[crossed];
      if (point[2] < z.low || point[2] > z.low + (static_cast<Scalar>(kSecond) / z.perUnit)) {
        continue;
      }
      if (link(previous, crossed) < 0 && link(crossed, previous) < 0) {
        return out;
      }
      const auto part = linkStep(previous, crossed, previousPoint, point);
      if (previous == fromLayer) {
        out.left = part.left;
      }
      add(part, previousT, t, false);
      previous = crossed;
      previousPoint[0] = point[0];
      previousPoint[1] = point[1];
      previousPoint[2] = point[2];
      previousT = t;
    }
    if (previous == fromLayer || (link(previous, toLayer) < 0 && link(toLayer, previous) < 0)) {
      return out;
    }
    const auto lastPart = linkStep(previous, toLayer, previousPoint, end);
    add(lastPart, previousT, 1.F, true);
    out.found = true;
    out.before = sum;
    out.after = lastPart.after;
    out.gapFirstHalf = firstHalf;
    if (sum > 0.F) {
      const Scalar mean = first / sum;
      out.lever = 1.F - mean;
      const Scalar spread = (second / sum) - (mean * mean);
      out.spread = spread > 0.F ? spread : 0.F;
    }
    return out;
  }

  /// `step` along the link between \a fromLayer and \a toLayer, either way.
  [[nodiscard]] GPUhdi() LinkStep linkStep(LayerID fromLayer, LayerID toLayer, const Scalar start[3], const Scalar end[3]) const
  {
    using Math = o2::gpu::GPUCommonMath;
    LinkStep out;
    int l = link(fromLayer, toLayer);
    const bool inward = l < 0;
    if (inward) {
      l = link(toLayer, fromLayer);
      if (l < 0) {
        return out;
      }
    }
    out.found = true;
    const Scalar dx = end[0] - start[0];
    const Scalar dy = end[1] - start[1];
    const Scalar dz = end[2] - start[2];
    const Scalar length = Math::Sqrt((dx * dx) + (dy * dy) + (dz * dz));
    if (!(length > 0.F)) {
      return out;
    }
    const Scalar* source = inward ? end : start;
    const Scalar* target = inward ? start : end;
    const LayerID sourceLayer = linkFrom[l];
    const LayerID targetLayer = linkTo[l];
    const Scalar rSource = Math::Sqrt((source[0] * source[0]) + (source[1] * source[1]));
    const Scalar rTarget = Math::Sqrt((target[0] * target[0]) + (target[1] * target[1]));
    const auto& g = gaps[gapRow(static_cast<std::size_t>(l), source[0], source[1], second(sourceLayer, rSource, source[2])) + static_cast<std::uint32_t>(layerSecond[targetLayer].bin(second(targetLayer, rTarget, target[2]), kSecond))];
    Scalar oblique = 1.F;
    if (Detector<Det>::isBarrelLayer(sourceLayer) && Detector<Det>::isBarrelLayer(targetLayer)) {
      oblique = Math::Sqrt((dx * dx) + (dy * dy)) / Math::Abs(rTarget - rSource);
    }
    // The cosine of the step's incidence on a layer at the point p.
    const auto cosine = [&](LayerID layer, const Scalar* p) {
      if (Detector<Det>::isBarrelLayer(layer)) {
        const Scalar r = Math::Sqrt((p[0] * p[0]) + (p[1] * p[1]));
        return r > 0.F ? Math::Abs((dx * p[0]) + (dy * p[1])) / (r * length) : 1.F;
      }
      return Math::Abs(dz) / length;
    };
    const Scalar rStart = Math::Sqrt((start[0] * start[0]) + (start[1] * start[1]));
    const Scalar rEnd = Math::Sqrt((end[0] * end[0]) + (end[1] * end[1]));
    const auto& left = node(fromLayer, start[0], start[1], second(fromLayer, rStart, start[2]));
    const auto& reached = node(toLayer, end[0], end[1], second(toLayer, rEnd, end[2]));
    // The share of each band in front of the point the step leaves or reaches.
    const Scalar leftFront = left.frontShare(depthOf(fromLayer, rStart, start[2]));
    const Scalar reachedFront = reached.frontShare(depthOf(toLayer, rEnd, end[2]));
    const Scalar cosLeft = cosine(fromLayer, start);
    const Scalar cosReached = cosine(toLayer, end);
    const Scalar pathLeft = cosLeft > 0.F ? 1.F / cosLeft : 1.F;
    const Scalar pathReached = cosReached > 0.F ? 1.F / cosReached : 1.F;
    // Outward a band's part behind the point is left and the part in front reached; inward the reverse.
    const Scalar leftShare = inward ? leftFront : 1.F - leftFront;
    const Scalar reachedShare = inward ? 1.F - reachedFront : reachedFront;
    const Scalar leftX = left.x * leftShare * pathLeft;
    const Scalar reachedX = reached.x * reachedShare * pathReached;
    const Scalar gapX = g.x * oblique;
    // The gap's moments in the step's direction.
    const Scalar m1 = inward ? 1.F - g.m1 : g.m1;
    const Scalar m2 = inward ? 1.F - (2.F * g.m1) + g.m2 : g.m2;
    out.before = leftX + gapX;
    out.after = reachedX;
    out.left = leftX;
    out.gapFirstHalf = gapX * (inward ? 1.F - g.half : g.half);
    if (out.before > 0.F) {
      const Scalar mean = gapX * m1 / out.before;
      const Scalar square = gapX * m2 / out.before;
      out.lever = 1.F - mean;
      const Scalar spread = square - (mean * mean);
      out.spread = spread > 0.F ? spread : 0.F;
    }
    out.xRho = (left.xRho * leftShare * pathLeft) + (g.xRho * oblique) + (reached.xRho * reachedShare * pathReached);
    return out;
  }

  /// Layer \a layer's vertex gap for a vertex at \a z and a point of the
  /// layer at azimuth \a phi with second coordinate \a second.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() const MaterialGap& vertexGapAt(LayerID layer, Scalar phi, Scalar z, Scalar second) const
  {
    const auto bin = static_cast<std::size_t>(phiBinOf(phi, kGapPhi));
    const auto v = static_cast<std::size_t>(vertexZ.bin(z, kVertexZ));
    const auto s = static_cast<std::size_t>(layerSecond[layer].bin(second, kSecond));
    return vertexGaps[(((((static_cast<std::size_t>(layer) * kGapPhi) + bin) * kVertexZ) + v) * kSecond) + s];
  }

  /// Layer \a layer's vertex gap for a vertex at \a z and a point (\a x,
  /// \a y) of the layer with second coordinate \a second.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() const MaterialGap& vertexGap(LayerID layer, Scalar x, Scalar y, Scalar z, Scalar second) const { return vertexGaps[vertexIndex(layer, x, y, z, second)]; }

  /// The index of that vertex gap (and its profile).
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::size_t vertexIndex(LayerID layer, Scalar x, Scalar y, Scalar z, Scalar second) const
  {
    const auto phi = static_cast<std::size_t>(phiBin(x, y, kGapPhi));
    const auto v = static_cast<std::size_t>(vertexZ.bin(z, kVertexZ));
    const auto s = static_cast<std::size_t>(layerSecond[layer].bin(second, kSecond));
    return (((((static_cast<std::size_t>(layer) * kGapPhi) + phi) * kVertexZ) + v) * kSecond) + s;
  }

  /// The middle of layer \a layer's band: a radius on a barrel layer, |z| on a disc.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar bandMiddle(LayerID layer) const { return layerDepth[layer].low + (0.5F / layerDepth[layer].perUnit); }

  /// The material of the straight segment from \a a to \a b, each on layer
  /// \a aLayer or \a bLayer or, kInvalidLayerID, a point inside the
  /// detector off its layers (a vertex, a decay's): between two layers
  /// `step`; otherwise a share of a vertex gap -- of the layer whose reach
  /// (the volume inside its band) holds the free points, the layer of a point
  /// on one if the other lies inside it -- between the points' fractions of
  /// its ray (`GapProfile`), the fractions radial on a barrel layer and along
  /// z on a disc, its vertex where the segment's line passes the beam axis
  /// closest; and the band in front of a point on the layer of that gap. In
  /// the segment's direction, `after` that band at \a b. `found` false where
  /// no layer's reach holds the free points.
  [[nodiscard]] GPUhdi() LinkStep segment(LayerID aLayer, const Scalar a[3], LayerID bLayer, const Scalar b[3]) const
  {
    using Math = o2::gpu::GPUCommonMath;
    if (aLayer != kInvalidLayerID && bLayer != kInvalidLayerID) {
      return step(aLayer, bLayer, a, b);
    }
    LinkStep out;
    const Scalar d[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const Scalar length = Math::Sqrt((d[0] * d[0]) + (d[1] * d[1]) + (d[2] * d[2]));
    if (!(length > 0.F)) {
      out.found = true;
      return out;
    }
    const Scalar ra = Math::Sqrt((a[0] * a[0]) + (a[1] * a[1]));
    const Scalar rb = Math::Sqrt((b[0] * b[0]) + (b[1] * b[1]));
    // Whether a point lies inside layer l's reach (inward of its band).
    const auto inside = [&](LayerID l, Scalar r, Scalar z) {
      const Scalar inner = layerDepth[l].low;
      if (Detector<Det>::isBarrelLayer(l)) {
        const auto& axis = layerSecond[l];
        return r < inner && z >= axis.low && z <= axis.low + (static_cast<Scalar>(kSecond) / axis.perUnit);
      }
      return Math::Abs(z) < inner && z * static_cast<Scalar>(layerSide[l]) >= 0.F;
    };
    // The layer of the gap: a point's on a layer if the other is inside it, else the first holding both.
    LayerID layer = kInvalidLayerID;
    bool aOn = false;
    bool bOn = false;
    if (aLayer != kInvalidLayerID && inside(aLayer, rb, b[2])) {
      layer = aLayer;
      aOn = true;
    } else if (bLayer != kInvalidLayerID && inside(bLayer, ra, a[2])) {
      layer = bLayer;
      bOn = true;
    } else {
      // the innermost barrel layer holding both, else the nearest disc
      for (int l = 0; l < kLayers; ++l) {
        const auto candidate = static_cast<LayerID>(l);
        if (!inside(candidate, ra, a[2]) || !inside(candidate, rb, b[2])) {
          continue;
        }
        const bool better = layer == kInvalidLayerID || (Detector<Det>::isBarrelLayer(candidate) != Detector<Det>::isBarrelLayer(layer) ? Detector<Det>::isBarrelLayer(candidate) : bandMiddle(candidate) < bandMiddle(layer));
        if (better) {
          layer = candidate;
        }
      }
    }
    if (layer == kInvalidLayerID) {
      return out;
    }
    const bool barrel = Detector<Det>::isBarrelLayer(layer);
    // The vertex: the line's closest approach to the beam axis.
    const Scalar dt2 = (d[0] * d[0]) + (d[1] * d[1]);
    const Scalar tAxis = dt2 > 0.F ? -((a[0] * d[0]) + (a[1] * d[1])) / dt2 : 0.F;
    const Scalar z0 = a[2] + (tAxis * d[2]);
    // The ray's end: the layer's band middle along the line, else the far point's direction.
    const Scalar middle = bandMiddle(layer);
    const Scalar* far = (barrel ? rb > ra : Math::Abs(b[2]) > Math::Abs(a[2])) ? b : a;
    const Scalar rFar = far == b ? rb : ra;
    Scalar n[3];
    if (barrel) {
      const Scalar scale = rFar > 0.F ? middle / rFar : 1.F;
      n[0] = far[0] * scale;
      n[1] = far[1] * scale;
      n[2] = z0 + ((far[2] - z0) * scale);
    } else {
      const Scalar zEnd = far[2] >= 0.F ? middle : -middle;
      const Scalar dz = far[2] - z0;
      const Scalar scale = Math::Abs(dz) > 0.F ? (zEnd - z0) / dz : 1.F;
      n[0] = far[0] * scale;
      n[1] = far[1] * scale;
      n[2] = zEnd;
    }
    const Scalar rn = Math::Sqrt((n[0] * n[0]) + (n[1] * n[1]));
    const auto fraction = [&](Scalar r, Scalar z) {
      if (barrel) {
        return rn > 0.F ? r / rn : 0.F;
      }
      const Scalar span = n[2] - z0;
      return Math::Abs(span) > 0.F ? (z - z0) / span : 0.F;
    };
    const Scalar fa = aOn ? 1.F : fraction(ra, a[2]);
    const Scalar fb = bOn ? 1.F : fraction(rb, b[2]);
    const auto index = vertexIndex(layer, n[0], n[1], z0, second(layer, rn, n[2]));
    const auto& g = vertexGaps[index];
    Scalar mean = 0.F;
    Scalar square = 0.F;
    const Scalar lo = fa < fb ? fa : fb;
    const Scalar hi = fa < fb ? fb : fa;
    const Scalar share = vertexProfiles[index].between(lo < 0.F ? 0.F : lo, hi > 1.F ? 1.F : hi, mean, square);
    // The segment's path per unit of depth (radius or z) it covers: the
    // gap is at normal incidence.
    const Scalar df = Math::Abs(fb - fa);
    const Scalar depth = df * (barrel ? rn : Math::Abs(n[2] - z0));
    const Scalar path = depth > 0.F ? length / depth : 1.F;
    const Scalar gapX = g.x * share * path;
    // The fraction of the segment, from a, a ray fraction f lies at.
    const Scalar span = fb - fa;
    const Scalar meanU = Math::Abs(span) > 0.F ? (mean - fa) / span : 0.5F;
    const Scalar squareU = Math::Abs(span) > 0.F ? (square - (2.F * fa * mean) + (fa * fa)) / (span * span) : 0.25F;
    // The band in front of a point on the gap's layer.
    Scalar nodeX = 0.F;
    Scalar nodeRho = 0.F;
    if (aOn || bOn) {
      const Scalar* p = aOn ? a : b;
      const Scalar r = aOn ? ra : rb;
      const auto& band = node(layer, p[0], p[1], second(layer, r, p[2]));
      const Scalar front = band.frontShare(depthOf(layer, r, p[2]));
      const Scalar cosine = barrel ? (r > 0.F ? Math::Abs((d[0] * p[0]) + (d[1] * p[1])) / (r * length) : 1.F) : Math::Abs(d[2]) / length;
      const Scalar factor = cosine > 0.F ? 1.F / cosine : 1.F;
      nodeX = band.x * front * factor;
      nodeRho = band.xRho * front * factor;
    }
    out.found = true;
    out.xRho = (g.xRho * share * path) + nodeRho;
    out.gapFirstHalf = gapX * (meanU < 0.5F ? 1.F : 0.F);
    if (bOn) {
      out.before = gapX;
      out.after = nodeX;
      out.lever = 1.F - meanU;
      const Scalar spread = squareU - (meanU * meanU);
      out.spread = spread > 0.F ? spread : 0.F;
    } else {
      // the band at a first, then the gap
      out.left = nodeX;
      out.before = nodeX + gapX;
      if (out.before > 0.F) {
        const Scalar gapMean = gapX * meanU / out.before;
        const Scalar gapSquare = gapX * squareU / out.before;
        out.lever = 1.F - gapMean;
        const Scalar spread = gapSquare - (gapMean * gapMean);
        out.spread = spread > 0.F ? spread : 0.F;
      }
    }
    return out;
  }
};

/// A step's x/X0 and the integral of t^2 dx/X0 over it, t from 0 at its start
/// to 1 at its end: `before` at its mean place and spread, `after` at the end.
struct StepMoment {
  Scalar x2X0 = 0.F;
  Scalar x2X0T2 = 0.F;
};

[[nodiscard]] GPUhdi() GPUAlwaysInline() StepMoment momentOf(const LinkStep& step)
{
  const Scalar mean = 1.F - step.lever;
  return {.x2X0 = step.before + step.after, .x2X0T2 = (step.before * (step.spread + (mean * mean))) + step.after};
}

/// Highland's `theta0^2 (beta p)^2` of \a x = x/X0 along the path, its
/// logarithm at \a logX = ln x' -- the table's own material, a node's at
/// normal incidence, a gap's along its chord, so that no lookup takes a
/// logarithm -- and \a logBeta2 = ln(beta^2); 0 without material.
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar materialKick(Scalar x, Scalar logX, Scalar logBeta2)
{
  if (!(x > 0.F)) {
    return 0.F;
  }
  const Scalar highland = 1.F + (0.038F * (logX - logBeta2));
  return highland > 0.F ? 0.0136F * 0.0136F * x * highland * highland : 0.F;
}

/// The same with the logarithm at \a x itself.
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar materialKick(Scalar x, Scalar logBeta2) { return x > 0.F ? materialKick(x, logOf(x), logBeta2) : 0.F; }

} // namespace o2::gtrk

#endif // GTRK_TRACKING_LINKMATERIAL_H_
