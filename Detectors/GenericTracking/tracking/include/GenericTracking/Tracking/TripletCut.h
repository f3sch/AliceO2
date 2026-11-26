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

/// \file TripletCut.h
/// \brief The cuts a triplet passes.
///
/// A triplet is a tracklet a -> b on one link and b -> c on the next. It
/// passes, cheapest first:
///
/// 1. the dip-angle cut, `|lambda_bc - lambda_ab|` within a half-width that
///    depends on `a -> b` alone, so the triplet stage finds its candidates
///    by a binary search in `b`'s lambda-sorted tracklets;
/// 2. the boxes of `a` and `c` overlap;
/// 3. the circle's curvature against pT_min;
/// 4. the triplet fit's chi2;
/// 5. the fit again with the scattering at the triplet's own momentum
///    (`tripletOwnMomentumCompatible`).
///
/// A triplet failing 3-5 is tried once more with the iteration's kink
/// allowance on the middle layer (`tripletCompatible`).
///
/// Every tolerance follows from the scattering, the measurements' errors and
/// pT_min, with the iteration's `nSigma`.

#ifndef GTRK_TRACKING_TRIPLETCUT_H_
#define GTRK_TRACKING_TRIPLETCUT_H_

#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Tracking/Iteration.h"
#include "GenericTracking/Tracking/MathUtils.h"
#include "GenericTracking/Tracking/Measurement.h"
#include "GenericTracking/Tracking/PairCut.h"
#include "GenericTracking/Tracking/TripletFit.h"
#include "GenericTracking/Tracking/WorkCounters.h"

#include "GPUCommonMath.h"

#include <cstddef>

namespace o2::gtrk
{

/// A measurement as the triplet cuts read it.
struct TimedPoint {
  PairPoint point;
  TimeEstPS time; ///< the time box
};

/// Row \a row of \a columns as a `TimedPoint`, reading only the columns it holds.
template <typename Columns>
[[nodiscard]] GPUhdi() GPUAlwaysInline() TimedPoint timedPointAt(const Columns& columns, std::uint32_t row)
{
  return {.point = pairPointAt(columns, row), .time = timeBoxAt(columns, row)};
}

/// Bounds on the third measurement's variances over the rows a triplet can
/// read on its layer: the dip-angle half-width is fixed before they are seen.
struct TripletTargetBounds {
  Scalar varZ = 0.F;
  Scalar varR = 0.F;
  Scalar varRPhi = 0.F;
};

namespace detail
{

/// How far a chord's dip angle can lie from the helix's, at most, for the
/// chord's half-turn sine \a halfTurnSine and a chord slope \a slope.
///
/// `tan(lambda_chord) = tan(lambda) f(x)`, `f(x) = asin(x)/x`, so the chord
/// is steeper by `atan(t f) - atan(t)`. That grows with the true slope `t` up
/// to `t = 1/sqrt(f)`, where it is `atan((f - 1)/(2 sqrt(f)))`, and falls
/// beyond. The chord's slope bounds the true one from above.
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar chordDipExcess(Scalar halfTurnSine, Scalar slope, Scalar atanSlope)
{
  const Scalar f = asinOverArgument(o2::gpu::CAMath::Min(halfTurnSine, 1.F));
  const Scalar peak = 1.F / o2::gpu::CAMath::Sqrt(f);
  if (slope < peak) {
    return atanOf(slope * f) - atanSlope;
  }
  return atanOf((f - 1.F) / (2.F * o2::gpu::CAMath::Sqrt(f)));
}

/// A measurement's variance across the chord with direction (\a ux, \a uy):
/// its r and r-phi errors projected on the chord's normal.
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar varianceAcross(const PairPoint& point, Scalar ux, Scalar uy)
{
  const Scalar invR = point.r > 0.F ? 1.F / point.r : 0.F;
  const Scalar radialAcross = ((point.x * uy) - (point.y * ux)) * invR;
  const Scalar tangentialAcross = ((point.x * ux) + (point.y * uy)) * invR;
  return (point.varR * radialAcross * radialAcross) + (point.varRPhi * tangentialAcross * tangentialAcross);
}

} // namespace detail

/// What `tripletLambdaHalfWidth` takes from the first tracklet `a -> b` alone,
/// made once per first tracklet (`tripletFirstTerms`) and shared by the
/// out-links of `b`: the chord's dip, the scattering of `b`'s node and of the
/// first gap, where `b` reads the second links' gaps, and the first two
/// measurements' error terms.
template <DetectorLayout Det>
struct TripletFirstTerms {
  bool valid = false;           ///< the tracklet has a length
  Scalar chord = 0.F;           ///< its transverse chord
  Scalar sine = 0.F;            ///< |sin(lambda)|
  Scalar cosine = 0.F;          ///< cos(lambda)
  Scalar slope = 0.F;           ///< tan(lambda), 0 where cos(lambda) is 0
  Scalar atanSlope = 0.F;       ///< atan(slope)
  Scalar scattering = 0.F;      ///< `b`'s node with the kink allowance, then the gap between `a` and `b` by its fraction
  std::uint32_t middlePhi = 0U; ///< `b`'s gap azimuth bin
  Scalar firstSpanX = 0.F;      ///< x/X0 of `b`'s node and the first gap
  Scalar middleSecond = 0.F;    ///< `b`'s second coordinate
  Scalar inFirst = 0.F;         ///< 1 / L_ab
  Scalar zFirst = 0.F;          ///< varZ_a / L_ab^2
  Scalar chordFirst = 0.F;      ///< (varR_a + varRPhi_a) / L_ab^2
  Scalar varZMiddle = 0.F;      ///< varZ_b
  Scalar alongMiddle = 0.F;     ///< varR_b + varRPhi_b
};

/// \a firstLink: the first tracklet's link.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() TripletFirstTerms<Det> tripletFirstTerms(const IterationParameters<Det>& parameters, std::uint32_t firstLink, const PairPoint& first, const PairPoint& middle)
{
  TripletFirstTerms<Det> terms;
  terms.chord = pairChord(first, middle).chord;
  const Scalar dz = middle.z - first.z;
  const Scalar length = o2::gpu::CAMath::Hypot(terms.chord, dz);
  if (!(length > 0.F)) {
    return terms;
  }
  terms.valid = true;
  terms.sine = o2::gpu::CAMath::Abs(dz) / length; // |sin(lambda)|
  terms.cosine = terms.chord / length;            // cos(lambda)
  terms.slope = terms.cosine > 0.F ? terms.sine / terms.cosine : 0.F;
  terms.atanSlope = atanOf(terms.slope);
  // The third measurement is not known yet: the first gap enters with its
  // fraction, the second link's with its largest from `b` (the half-width).
  // The middle node takes the iteration's kink allowance, so the window holds
  // the kinked triplets too (1 without one).
  using Table = LinkMaterial<Det>;
  const auto& material = *parameters.material;
  const LayerID firstLayer = material.linkFrom[firstLink];
  const LayerID middleLayer = material.linkTo[firstLink];
  const Scalar middleSecond = Table::second(middleLayer, middle.r, middle.z);
  const auto& node = material.node(middleLayer, middle.x, middle.y, middleSecond);
  TripletMaterial partial;
  partial.node = node.x;

  partial.barrel = Detector<Det>::isBarrelLayer(middleLayer);
  partial.first = material.gap(material.gapRow(firstLink, first.x, first.y, Table::second(firstLayer, first.r, first.z)), middleLayer, middleSecond);
  partial = partial.withSpan();
  terms.scattering = (middleScattering(partial, parameters.pairCut, terms.cosine) * parameters.graphKinkScale2) + gapScattering(partial.first, partial.first.m2, partial.spanLog, parameters.pairCut, terms.cosine);
  terms.firstSpanX = partial.node + partial.first.x;
  terms.middlePhi = static_cast<std::uint32_t>(LinkMaterial<Det>::phiBin(middle.x, middle.y, LinkMaterial<Det>::kGapPhi));
  terms.middleSecond = middleSecond;
  terms.inFirst = 1.F / length;
  terms.zFirst = first.varZ * terms.inFirst * terms.inFirst;
  const Scalar alongFirst = first.varR + first.varRPhi;
  terms.chordFirst = alongFirst * terms.inFirst * terms.inFirst;
  terms.varZMiddle = middle.varZ;
  terms.alongMiddle = middle.varR + middle.varRPhi;
  return terms;
}

/// The half-width of the dip-angle cut for the first tracklet `a -> b` and a
/// second tracklet on \a secondLink, whose target's variances are bounded by
/// \a target.
///
/// The sum of
/// - the chords' excess over the helix's dip angle at the curvature bound,
///   for the longer of `a -> b` and the second link's longest chord;
/// - `n` sigma of the scattering at `b`, at pT_min, and of the three
///   measurements' z and along-chord errors, with the second link's
///   shortest length and the target's largest variances:
///
/// ```
/// var = theta_0^2 + cos^2(lambda) [varZ_a/L_ab^2 + varZ_b (1/L_ab + 1/L_bc)^2 + varZ_c/L_bc^2]
///                 + sin^2(lambda) [the same along the chord]
/// ```
///
/// with `theta_0^2` the middle node's and the gaps' (`tripletScattering`): the
/// first gap by its fraction, the second link's at its largest from `b`.
///
/// The along-chord error of a measurement is bounded by its r plus r-phi
/// variance. Everything the target does not fix is taken where it is largest,
/// so the half-width holds for every candidate on the link.
///
/// \a first is the first tracklet's part (`tripletFirstTerms`), the same for
/// every out-link of `b`; this adds the link's.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar tripletLambdaHalfWidth(const IterationParameters<Det>& parameters, std::uint32_t secondLink, const TripletFirstTerms<Det>& first, const TripletTargetBounds& target)
{
  const auto& bounds = parameters.tripletBounds[secondLink];
  if (!first.valid || !(bounds.lengthMin > 0.F)) {
    return kPi<Scalar>;
  }
  const Scalar longest = o2::gpu::CAMath::Max(first.chord, bounds.chordMax);
  const Scalar excess = detail::chordDipExcess(0.5F * parameters.pairCut.maxCurvature * longest, first.slope, first.atanSlope);
  const auto& material = *parameters.material;
  const Scalar largest = material.rowMax(material.gapRowAt(secondLink, first.middlePhi, first.middleSecond));
  const Scalar spanX = first.firstSpanX + largest;
  const Scalar scattering = first.scattering + (materialKick(largest, spanX > 0.F ? logOf(spanX) : 0.F, parameters.pairCut.logBeta2) * first.cosine * first.cosine * parameters.pairCut.invBetaMomentum2);
  const Scalar inSecond = 1.F / bounds.lengthMin;
  const Scalar inMiddle = first.inFirst + inSecond;
  const Scalar alongTarget = target.varR + target.varRPhi;
  const Scalar zTerm = first.zFirst + (first.varZMiddle * inMiddle * inMiddle) + (target.varZ * inSecond * inSecond);
  const Scalar chordTerm = first.chordFirst + (first.alongMiddle * inMiddle * inMiddle) + (alongTarget * inSecond * inSecond);
  const Scalar variance = scattering + (first.cosine * first.cosine * zTerm) + (first.sine * first.sine * chordTerm);
  return excess + (parameters.pairCut.nSigma * o2::gpu::CAMath::Sqrt(variance));
}

/// The same for one link from the points, the first tracklet's part made here.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() Scalar tripletLambdaHalfWidth(const IterationParameters<Det>& parameters, std::uint32_t firstLink, std::uint32_t secondLink, const PairPoint& first, const PairPoint& middle, const TripletTargetBounds& target)
{
  return tripletLambdaHalfWidth(parameters, secondLink, tripletFirstTerms(parameters, firstLink, first, middle), target);
}

/// The first and the third measurement's boxes overlap: with the links' two
/// pairs, the three boxes share a time.
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool tripletTimeCompatible(const TimedPoint& first, const TimedPoint& third) { return first.time.isCompatible(third.time); }

/// The circle through the three measurements curves at most as the softest
/// track of the iteration, within `n` sigma.
///
/// Its curvature `2 (P_ab x P_bc) / (d_ab d_bc d_ac)` changes by `2/d_ac` times
/// any change of the angle between the two chords, which moves `a` by its
/// displacement over `d_ab`, `c` over `d_bc`, and `b` over both. The
/// displacements are the measurements' errors across the chord `a -> c`, and
/// the scattering at `b` turns the second chord by `theta_0 / cos(lambda)`.
/// \a middleScale multiplies
/// the middle node's scattering: the kink allowance (1 for the nominal cut).
/// \a chords are the points' (`tripletChords`), \a material the triplet's
/// (`tripletMaterial`).
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool tripletCurvatureCompatible(const IterationParameters<Det>& parameters, const TripletMaterial& material, const TripletChords& chords, const PairPoint& first, const PairPoint& middle, const PairPoint& third, Scalar middleScale = 1.F)
{
  const Scalar dxFirst = middle.x - first.x;
  const Scalar dyFirst = middle.y - first.y;
  const Scalar dxSecond = third.x - middle.x;
  const Scalar dySecond = third.y - middle.y;
  const Scalar dxWhole = third.x - first.x;
  const Scalar dyWhole = third.y - first.y;
  const Scalar chordFirst = chords.length01;
  const Scalar chordSecond = chords.length12;
  const Scalar chordWhole = chords.length02;
  if (!(chordFirst > 0.F) || !(chordSecond > 0.F) || !(chordWhole > 0.F)) {
    return false;
  }
  const Scalar cross = (dxFirst * dySecond) - (dyFirst * dxSecond);
  const Scalar curvature = 2.F * o2::gpu::CAMath::Abs(cross) / (chordFirst * chordSecond * chordWhole);
  const Scalar excess = curvature - parameters.pairCut.maxCurvature;
  if (!(excess > 0.F)) {
    return true;
  }

  const Scalar dz = third.z - first.z;
  const Scalar cosine2 = (chordWhole * chordWhole) / ((chordWhole * chordWhole) + (dz * dz)); // cos^2(lambda)
  const Scalar scattering = tripletScattering(material, parameters.pairCut, o2::gpu::CAMath::Sqrt(cosine2), middleScale) / cosine2;
  const Scalar ux = dxWhole / chordWhole;
  const Scalar uy = dyWhole / chordWhole;
  const Scalar inFirst = 1.F / chordFirst;
  const Scalar inSecond = 1.F / chordSecond;
  const Scalar inMiddle = inFirst + inSecond;
  const Scalar angleVariance = scattering + (detail::varianceAcross(first, ux, uy) * inFirst * inFirst) + (detail::varianceAcross(middle, ux, uy) * inMiddle * inMiddle) + (detail::varianceAcross(third, ux, uy) * inSecond * inSecond);
  const Scalar scale = 2.F / chordWhole;
  const Scalar nSigma = parameters.pairCut.nSigma;
  return excess * excess <= nSigma * nSigma * scale * scale * angleVariance;
}

/// The same, the chords made here.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool tripletCurvatureCompatible(const IterationParameters<Det>& parameters, const TripletMaterial& material, const PairPoint& first, const PairPoint& middle, const PairPoint& third, Scalar middleScale = 1.F)
{
  return tripletCurvatureCompatible(parameters, material, tripletChords(first, middle, third), first, middle, third, middleScale);
}

namespace detail
{

/// The fit's chi2 cut again, every scatterer's variance taken at the
/// triplet's own momentum instead of pT_min's:
/// the momentum of the fitted curvature raised by `nSigma` of its sigmas, the
/// variance never above pT_min's. A triplet as soft as pT_min keeps the fit
/// just made; a stiffer one is refitted from its kink terms. The window and
/// the first fit take pT_min's variance, as a track at pT_min needs; this
/// takes from a stiff-looking triplet the allowance only a soft track needs.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool tripletOwnMomentumCompatible(const IterationParameters<Det>& parameters, const TripletMaterial& material, const PairPoint& first, const PairPoint& middle, const PairPoint& third, Scalar middleScale, const TripletFitResult& fit, const KinkTerms& terms)
{
  const PairCut& nominal = parameters.pairCut;
  const Scalar raised = o2::gpu::CAMath::Abs(fit.curvature) + (nominal.nSigma * o2::gpu::CAMath::Sqrt(fit.varCurvature));
  if (!(nominal.momentumCurvature > 0.F) || !(raised > 0.F) || !(terms.sineTheta > 0.F)) {
    return true;
  }
  const Scalar momentum = nominal.momentumCurvature / raised;
  const Scalar pt = momentum * terms.sineTheta;
  const Scalar own = (1.F + (nominal.mass2 / (momentum * momentum))) / (pt * pt);
  if (!(own < nominal.invBetaMomentum2)) {
    return true;
  }
  PairCut cut = nominal;
  cut.invBetaMomentum2 = own;
  cut.logBeta2 = logOf(momentum * momentum / ((momentum * momentum) + nominal.mass2));
  TripletFitFactor factor;
  TripletFitResult refit;
  return makeTripletFitFactor(terms, first, middle, third, tripletScattering(material, cut, terms.sineTheta, middleScale), factor) && fitTriplet(factor, refit) && refit.chi2 <= nominal.nSigma * nominal.nSigma;
}

/// The curvature test and the fit with the middle layer's scattering times
/// \a middleScale: the nominal cuts at 1, the kink allowance at the
/// iteration's `graphKinkScale2`. The kink terms depend on the three points alone:
/// the first try to fit makes them into \a terms (\a made, and \a valid with
/// `makeKinkTerms`' result) and the other reuses them. Both read the points'
/// \a chords.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool tripletSpatiallyCompatible(const IterationParameters<Det>& parameters, const TripletMaterial& material, const TimedPoint& first, const TimedPoint& middle, const TimedPoint& third, const TripletChords& chords, Scalar middleScale, TripletFitResult& fit, KinkTerms& terms, bool& made, bool& valid)
{
  if (!tripletCurvatureCompatible(parameters, material, chords, first.point, middle.point, third.point, middleScale)) {
    countWork(WorkCounter::kTripletCurvature);
    return false;
  }
  countWork(WorkCounter::kTripletFits);
  TripletFitResult result;
  bool fitted = false;
  {
    const WorkTimer timer{WorkCounter::kTripletFitNS};
    if (!made) {
      valid = makeKinkTerms(chords, first.point, middle.point, third.point, terms);
      made = true;
    }
    fitted = valid && fitTriplet(terms, first.point, middle.point, third.point, material, parameters.pairCut, result, middleScale);
  }
  if (!fitted) {
    countWork(WorkCounter::kTripletFitFailed);
    return false;
  }
  const Scalar nSigma = parameters.pairCut.nSigma;
  if (!(result.chi2 <= nSigma * nSigma)) {
    countWork(WorkCounter::kTripletChi2);
    return false;
  }
  if (!tripletOwnMomentumCompatible(parameters, material, first.point, middle.point, third.point, middleScale, result, terms)) {
    countWork(WorkCounter::kTripletOwnMomentum);
    return false;
  }
  fit = result;
  return true;
}

} // namespace detail

/// Every cut after the dip angle; \a fit receives the triplet fit and \a
/// terms the triplet's kink terms when all pass. A triplet that fails the
/// curvature test, the fit or its refit at its own momentum -- not the time
/// cut -- is tried once more with the iteration's kink allowance on the middle
/// layer's scattering (`graphKinkScale2`), its kink
/// terms made once for both; if it passes then, \a kinked is set and \a fit is
/// the kinked fit, its variance enlarged. Without an allowance there is no
/// second try.
template <DetectorLayout Det>
[[nodiscard]] GPUhdi() GPUAlwaysInline() bool tripletCompatible(const IterationParameters<Det>& parameters, const TripletMaterial& material, const TimedPoint& first, const TimedPoint& middle, const TimedPoint& third, TripletFitResult& fit, KinkTerms& terms, bool& kinked)
{
  kinked = false;
  if (!tripletTimeCompatible(first, third)) {
    countWork(WorkCounter::kTripletTime);
    return false;
  }
  const auto chords = tripletChords(first.point, middle.point, third.point);
  bool made = false;
  bool valid = false;
  if (detail::tripletSpatiallyCompatible(parameters, material, first, middle, third, chords, 1.F, fit, terms, made, valid)) {
    return true;
  }
  if (!(parameters.graphKinkScale2 > 1.F)) {
    return false;
  }
  if (!detail::tripletSpatiallyCompatible(parameters, material, first, middle, third, chords, parameters.graphKinkScale2, fit, terms, made, valid)) {
    return false;
  }
  countWork(WorkCounter::kTripletKinked);
  kinked = true;
  return true;
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRIPLETCUT_H_
