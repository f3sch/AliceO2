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

namespace o2::gtrk
{

/// A gap's material: x/X0 along its chord and the x-weighted mean and mean
/// square of the fraction of the chord, from the source, it lies at.
struct MaterialGap {
  float x = 0.F;
  float m1 = 0.F;
  float m2 = 0.F;
  float logX = 0.F; ///< ln x, Highland's logarithm (0 without material)
  float xRho = 0.F; ///< rho x along the chord [g/cm^2], for the energy loss
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
  float logX = 0.F; ///< ln x, Highland's logarithm at normal incidence (0 without material)
  float xRho = 0.F; ///< rho x across the band at normal incidence [g/cm^2]
  std::uint8_t depth[kNodeDepths + 1]{};
  /// The share of the band in front of the fraction \a t of its depth.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar frontShare(Scalar t) const
  {
    const Scalar at = (t < 0.F ? 0.F : (t > 1.F ? 1.F : t)) * static_cast<Scalar>(kNodeDepths);
    const int k = at >= static_cast<Scalar>(kNodeDepths) ? kNodeDepths - 1 : static_cast<int>(at);
    const Scalar low = static_cast<Scalar>(depth[k]);
    return (low + ((static_cast<Scalar>(depth[k + 1]) - low) * (at - static_cast<Scalar>(k)))) * (1.F / 255.F);
  }
};

/// One step's material between two linked layers' crossings, in the step's
/// direction (`LinkMaterial::step`): the band it leaves (the part on its far
/// side from where the track came) and the gap as one kick at the mean place
/// with a spread, the band it reaches (the part in front of its sensor) at
/// the end.
struct LinkStep {
  Scalar before = 0.F; ///< x/X0 along the step: the band left and the gap
  Scalar after = 0.F;  ///< x/X0 along the step: the band reached, before its sensor
  Scalar lever = 1.F;  ///< 1 - <f> of `before`, f the fraction of the step
  Scalar spread = 0.F; ///< <f^2> - <f>^2 of `before`
  Scalar xRho = 0.F;   ///< the whole step's [g/cm^2]
  bool found = false;  ///< the layers are a link (either way)
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

  /// The link from layer `from` to layer `to`, -1 for none; each link's layers.
  std::array<std::int16_t, static_cast<std::size_t>(kLayers) * kLayers> linkOf{};
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
  std::array<MaterialNode, static_cast<std::size_t>(kLayers) * kNodePhi * kSecond> nodes{};
  std::array<MaterialGap, static_cast<std::size_t>(kLinks) * kGapPhi * kSourceSecond * kSecond> gaps{};
  /// Per row of a link's gaps (one source azimuth and coordinate), the
  /// largest `x (1 - 2 m1 + m2)` over the targets: what a triplet's dip-angle
  /// window takes for its second link before the third point is known.
  std::array<float, static_cast<std::size_t>(kLinks) * kGapPhi * kSourceSecond> gapRowMax{};
  /// Per link, its gaps' largest x/X0: what bounds a pair's span before its gap is read.
  std::array<float, kLinks> gapMax{};
  std::array<MaterialGap, static_cast<std::size_t>(kLayers) * kGapPhi * kVertexZ * kSecond> vertexGaps{};

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
  /// second coordinate \a second.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint32_t gapRow(std::size_t link, Scalar x, Scalar y, Scalar second) const
  {
    const auto phi = static_cast<std::uint32_t>(phiBin(x, y, kGapPhi));
    const auto s = static_cast<std::uint32_t>(linkSource[link].bin(second, kSourceSecond));
    return ((((static_cast<std::uint32_t>(link) * kGapPhi) + phi) * kSourceSecond) + s) * kSecond;
  }

  /// The same for a source in gap azimuth bin \a phiBin (`phiBin` with `kGapPhi`).
  [[nodiscard]] GPUhdi() GPUAlwaysInline() std::uint32_t gapRowAt(std::size_t link, std::uint32_t phiBin, Scalar second) const
  {
    const auto s = static_cast<std::uint32_t>(linkSource[link].bin(second, kSourceSecond));
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
  /// along a link outward, against one inward; `found` false where the layers
  /// are no link. A node takes the step's incidence on its layer, radial on a
  /// barrel layer, along z on a disc; a barrel-to-barrel gap the chord's
  /// obliquity in the transverse plane over the table's radial chord.
  [[nodiscard]] GPUhdi() LinkStep step(LayerID fromLayer, LayerID toLayer, const Scalar start[3], const Scalar end[3]) const
  {
    using Math = o2::gpu::GPUCommonMath;
    LinkStep out;
    if (fromLayer >= static_cast<LayerID>(kLayers) || toLayer >= static_cast<LayerID>(kLayers)) {
      return out;
    }
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
  [[nodiscard]] GPUhdi() GPUAlwaysInline() const MaterialGap& vertexGap(LayerID layer, Scalar x, Scalar y, Scalar z, Scalar second) const
  {
    const auto phi = static_cast<std::size_t>(phiBin(x, y, kGapPhi));
    const auto v = static_cast<std::size_t>(vertexZ.bin(z, kVertexZ));
    const auto s = static_cast<std::size_t>(layerSecond[layer].bin(second, kSecond));
    return vertexGaps[(((((static_cast<std::size_t>(layer) * kGapPhi) + phi) * kVertexZ) + v) * kSecond) + s];
  }
};

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
