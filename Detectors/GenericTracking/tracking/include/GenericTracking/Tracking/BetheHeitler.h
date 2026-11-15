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

/// \file BetheHeitler.h
/// \brief An electron's bremsstrahlung loss over a step as a Gaussian mixture, and one component's loss applied to a state.
///
/// Bremsstrahlung leaves an electron the fraction z = E_after / E_before of
/// its energy. Bethe and Heitler's model over t = x / X0 radiation lengths
/// makes -ln z Gamma-distributed with shape t / ln 2 and unit scale: E[z^k]
/// = (1 + k)^(-t / ln 2), so E[z] = e^-t. The density is far from Gaussian
/// (most tracks lose little, a few lose most of their energy), so a
/// Gaussian-sum filter takes it as a mixture of Gaussians in z, each a
/// component's weight, mean and variance; each fit component is carried on
/// with each loss component.
///
/// The momentum is taken as the energy (an electron's mass is far below the
/// tracking's momenta), so z scales the momentum.

#ifndef GTRK_TRACKING_BETHEHEITLER_H_
#define GTRK_TRACKING_BETHEHEITLER_H_

#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/MathUtils.h"

#include "GPUCommonMath.h"

namespace o2::gtrk
{

/// One Gaussian of the loss mixture over the retained energy fraction z.
struct BetheHeitlerComponent {
  Scalar weight = 0.F;
  Scalar mean = 0.F;
  Scalar variance = 0.F;
};

/// The most components a step's mixture has.
inline constexpr int kBetheHeitlerComponents = 6;

/// The thickest step the mixture describes [x / X0]; a thicker step's loss is
/// the product of the losses of its parts (the shapes of the Gamma add), so
/// the caller applies it in parts no thicker than this.
inline constexpr Scalar kBetheHeitlerMaxX2X0 = 0.2F;

namespace detail
{

/// Below this thickness the step takes no loss.
inline constexpr Scalar kBetheHeitlerNoLossX2X0 = 1e-4F;
/// Below this thickness the loss is one Gaussian with Bethe-Heitler's mean
/// and variance.
inline constexpr Scalar kBetheHeitlerOneGaussianX2X0 = 0.002F;

/// The six-component parameterisation ACTS's Gaussian-sum filter takes
/// (Examples/Configs/betheHeitler_geantSim_cdf_nC6_O5.json, fitted to
/// Geant4's losses), with the two thresholds above as ACTS sets them: in
/// each of two ranges of t, per component, the polynomials in t of order
/// five, highest power first, of the logit of the weight, the logit of the
/// mean and the logarithm of the variance; the weights are normalised after.
/// Up to a few percent of a radiation length it is far closer to Bethe and
/// Heitler's tails and mean than ACTS's built-in default (ATLAS's
/// BetheHeitler_cdf_nC6_O5), which overstates them several times there.
/// On a GPU the table is in constant memory, which host code in a device
/// compilation unit cannot read: only device functions read it.
inline constexpr Scalar kBetheHeitlerRangeSplitX2X0 = 0.1F;
GPUconstexpr() Scalar kBetheHeitlerCoefficients[2][kBetheHeitlerComponents][3][6] = {
  // 0.0 <= t <= 0.1
  {{{3.21010941e+06F, -9.47056307e+05F, 1.08298759e+05F, -6.18635540e+03F, 2.01353611e+02F, -8.07615150e+00F},
    {4.50794410e+05F, -1.06471811e+05F, 9.01391148e+03F, -3.60813661e+02F, 7.03767376e+00F, -2.36102850e+00F},
    {7.45862413e+05F, -1.66341886e+05F, 1.28384366e+04F, -4.51208283e+02F, 8.18425268e+00F, -5.91975367e+00F}},
   {{3.58843797e+06F, -1.05317566e+06F, 1.19346411e+05F, -6.70611908e+03F, 2.10950205e+02F, -6.55684191e+00F},
    {0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, -1.50620196e+00F, -3.15693080e-01F},
    {0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, -1.30276933e+00F, -3.54004817e+00F}},
   {{2.31554526e+06F, -7.18915013e+05F, 8.65851581e+04F, -5.17997620e+03F, 1.74413154e+02F, -5.33349363e+00F},
    {0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, -3.61142838e+00F, 1.59829983e+00F},
    {0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, -9.86977145e-01F, -4.02968219e+00F}},
   {{2.54927417e+06F, -7.85993455e+05F, 9.37375917e+04F, -5.53264239e+03F, 1.81186194e+02F, -5.23901537e+00F},
    {-3.10981530e+05F, 1.10072484e+05F, -1.42801441e+04F, 8.41090255e+02F, -2.59709043e+01F, 3.18923109e+00F},
    {-5.55773386e+05F, 1.31543581e+05F, -1.16947471e+04F, 5.61054360e+02F, -1.55306777e+01F, -6.16186123e+00F}},
   {{2.69321543e+06F, -8.22500526e+05F, 9.68999597e+04F, -5.62972948e+03F, 1.79511682e+02F, -4.69839512e+00F},
    {0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, -5.48123172e+00F, 4.75133157e+00F},
    {-5.28365651e+05F, 1.24095313e+05F, -1.09310906e+04F, 5.25092571e+02F, -1.33013679e+01F, -9.16216201e+00F}},
   {{-2.70113597e+06F, 8.21801931e+05F, -9.66584667e+04F, 5.62268084e+03F, -1.87722128e+02F, 3.86177487e+00F},
    {-5.20768893e+05F, 1.59496294e+05F, -1.89990449e+04F, 1.10221067e+03F, -3.44364661e+01F, 7.43762383e+00F},
    {6.47667253e+04F, -2.69949965e+04F, 2.08766018e+03F, 7.66090494e+01F, 2.69959344e+00F, -1.63386708e+01F}}},
  // 0.1 <= t <= 0.2
  {{{0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, 3.41770310e+00F, -4.50233705e+00F},
    {0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, -3.94570530e+00F, -2.07634548e+00F},
    {0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, -7.14571219e+00F, -5.36289790e+00F}},
   {{0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, 1.14516875e+00F, -2.84954253e+00F},
    {0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, -5.29481145e+00F, -9.81800690e-02F},
    {0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, -6.26827496e+00F, -3.18603845e+00F}},
   {{-1.61214395e+04F, 1.53848432e+04F, -5.66893389e+03F, 9.88616066e+02F, -7.47160038e+01F, -3.86654892e-04F},
    {0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, -5.38920094e+00F, 1.66461070e+00F},
    {8.54188653e+04F, -8.47298452e+04F, 3.25098743e+04F, -6.03602323e+03F, 5.47032072e+02F, -2.32697644e+01F}},
   {{-2.17533635e+03F, 2.54805371e+03F, -1.08001532e+03F, 1.98017276e+02F, -1.10046967e+01F, -1.91224118e+00F},
    {-7.23982224e+04F, 7.16920265e+04F, -2.75401692e+04F, 5.12730379e+03F, -4.67418036e+02F, 1.91089862e+01F},
    {1.19909189e+05F, -1.18220722e+05F, 4.51916880e+04F, -8.36935920e+03F, 7.57724284e+02F, -3.27470300e+01F}},
   {{3.93011343e+03F, -3.75414183e+03F, 1.44257955e+03F, -2.93772911e+02F, 3.35529243e+01F, -3.11305862e+00F},
    {0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, 0.00000000e+00F, -6.11381645e+00F, 4.60839585e+00F},
    {1.57534242e+05F, -1.56440302e+05F, 6.02868488e+04F, -1.12640816e+04F, 1.03039079e+03F, -4.54372232e+01F}},
   {{3.23971641e+04F, -3.21141345e+04F, 1.22941040e+04F, -2.24702096e+03F, 1.86863452e+02F, -5.67098495e+00F},
    {-5.36173406e+04F, 5.43028466e+04F, -2.14355586e+04F, 4.12704029e+03F, -3.97505376e+02F, 2.17893222e+01F},
    {1.72155367e+05F, -1.72774707e+05F, 6.71741467e+04F, -1.26027591e+04F, 1.15227854e+03F, -5.38272103e+01F}}}};

GPUdi() GPUAlwaysInline() Scalar betheHeitlerPolynomial(int range, int component, int quantity, Scalar t)
{
  Scalar sum = 0.F;
  for (int power = 0; power < 6; ++power) {
    sum = (sum * t) + kBetheHeitlerCoefficients[range][component][quantity][power];
  }
  return sum;
}

GPUhdi() GPUAlwaysInline() Scalar logistic(Scalar x) { return 1.F / (1.F + o2::gpu::GPUCommonMath::Exp(-x)); }

} // namespace detail

/// The loss mixture of a step through \a x2X0 radiation lengths, in
/// \a mixture; returns its number of components: 1 or
/// kBetheHeitlerComponents, 0 for a step thicker than kBetheHeitlerMaxX2X0
/// (or a negative or NaN thickness).
GPUdi() int betheHeitlerMixture(Scalar x2X0, BetheHeitlerComponent (&mixture)[kBetheHeitlerComponents])
{
  if (!(x2X0 >= 0.F && x2X0 <= kBetheHeitlerMaxX2X0)) {
    return 0;
  }
  if (x2X0 < detail::kBetheHeitlerNoLossX2X0) {
    mixture[0] = {.weight = 1.F, .mean = 1.F, .variance = 0.F};
    return 1;
  }
  if (x2X0 < detail::kBetheHeitlerOneGaussianX2X0) {
    // E[z] = 2^-c = e^-t, E[z^2] = 3^-c, c = t / ln 2
    constexpr Scalar kLn3OverLn2 = 1.5849625F;
    const Scalar mean = o2::gpu::GPUCommonMath::Exp(-x2X0);
    mixture[0] = {.weight = 1.F, .mean = mean, .variance = o2::gpu::GPUCommonMath::Exp(-kLn3OverLn2 * x2X0) - (mean * mean)};
    return 1;
  }
  const int range = x2X0 < detail::kBetheHeitlerRangeSplitX2X0 ? 0 : 1;
  Scalar sum = 0.F;
  for (int k = 0; k < kBetheHeitlerComponents; ++k) {
    mixture[k] = {.weight = detail::logistic(detail::betheHeitlerPolynomial(range, k, 0, x2X0)),
                  .mean = detail::logistic(detail::betheHeitlerPolynomial(range, k, 1, x2X0)),
                  .variance = o2::gpu::GPUCommonMath::Exp(detail::betheHeitlerPolynomial(range, k, 2, x2X0))};
    sum += mixture[k].weight;
  }
  for (int k = 0; k < kBetheHeitlerComponents; ++k) {
    mixture[k].weight /= sum;
  }
  return kBetheHeitlerComponents;
}

/// One loss component applied to \a state and its covariance \a factors:
/// along the momentum q/p -> q/p / z, against it (an inward pass undoing the
/// loss) q/p -> q/p z, with z the component's mean and its variance carried
/// to q/p to first order about it. The covariance's q/p row and column are
/// scaled by the same factor, and the loss's variance added on q/p: along
/// the momentum (q/p)^2 var(z) / z^4, against it (q/p)^2 var(z). Against the
/// momentum q/p is linear in z, so the mean is exact there. False if the
/// result is not finite or the factors not positive.
GPUhdi() bool applyBetheHeitler(const BetheHeitlerComponent& component, PropagationDirection direction, BoundTrackState& state, CovarianceFactors& factors)
{
  const bool along = direction == PropagationDirection::kForward;
  const Scalar scale = along ? 1.F / component.mean : component.mean;
  const Scalar qOverP = state.par[kQOverP] * scale;
  state.par[kQOverP] = qOverP;
  for (int row = 0; row < kQOverP; ++row) {
    factors.upper[detail::strictUpperIndex(row, kQOverP)] /= scale;
  }
  factors.diagonal[kQOverP] *= scale * scale;
  const Scalar slope = along ? qOverP * scale : state.par[kQOverP] / scale;
  Scalar unit[kNBoundParams]{};
  unit[kQOverP] = 1.F;
  detail::rankOneUpdateUD(factors.upper, factors.diagonal, slope * slope * component.variance, unit);
  return isFiniteScalar(qOverP) && positiveDiagonal(factors);
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_BETHEHEITLER_H_
