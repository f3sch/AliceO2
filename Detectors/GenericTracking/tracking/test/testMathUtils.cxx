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

/// \file testMathUtils.cxx
/// \brief Exercises scalar helpers shared by the tracking code.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test GTRK MathUtils
#include "GenericTracking/Tracking/MathUtils.h"

#include "CommonConstants/MathConstants.h"

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

namespace o2::gtrk
{

BOOST_AUTO_TEST_CASE(TransverseArcHasTheStraightAndCurvedLimits)
{
  static_assert(std::is_same_v<decltype(transverseArcAtRadius(1.F, 2.F)), float>);
  static_assert(std::is_same_v<decltype(transverseArcAtRadius(1., 2.)), double>);

  BOOST_TEST(transverseArcAtRadius(0., 4.) == 4.);
  BOOST_TEST(transverseArcAtRadius(-1., 4.) == 4.);

  constexpr double curvature = 0.2;
  constexpr double radius = 4.;
  const double expected = 2. * std::asin(0.5 * curvature * radius) / curvature;
  BOOST_TEST(transverseArcAtRadius(curvature, radius) == expected, boost::test_tools::tolerance(1e-14));
}

/// The quotient from one asin with its derivative is, bit for bit, the
/// quotient alone, over the domain and at its ends.
BOOST_AUTO_TEST_CASE(AsinOverArgumentWithDerivativeIsTheQuotient)
{
  std::int32_t mismatches = 0;
  for (std::int32_t i = -20000; i <= 20000; ++i) {
    const float x = static_cast<float>(i) / 20000.F;
    const auto both = detail::asinOverArgumentWithDerivative(x);
    mismatches += sameBits(both.value, detail::asinOverArgument(x)) ? 0 : 1;
  }
  for (const float x : {0.05F, -0.05F, std::nextafter(0.05F, 0.F), 1.F, -1.F}) {
    const auto both = detail::asinOverArgumentWithDerivative(x);
    mismatches += sameBits(both.value, detail::asinOverArgument(x)) ? 0 : 1;
  }
  BOOST_TEST(mismatches == 0);
}

BOOST_AUTO_TEST_CASE(ChiSquareQuantilesMatchKnownValues)
{
  // Three degrees of freedom at the 68.27 percent and 99.73 percent
  // two-sided Gaussian confidence levels.
  BOOST_TEST(chiSquareQuantile(3, 1.) == 3.5267, boost::test_tools::tolerance(1.e-4));
  BOOST_TEST(chiSquareQuantile(3, 3.) == 14.1564, boost::test_tools::tolerance(1.e-4));
  // One degree of freedom is n^2 exactly; twenty was checked by numerical
  // integration of the density.
  BOOST_TEST(chiSquareQuantile(1, 3.) == 9., boost::test_tools::tolerance(1.e-6));
  BOOST_TEST(chiSquareQuantile(20, 3.) == 42.0802, boost::test_tools::tolerance(1.e-4));
}

namespace
{

// References in long double (64-bit mantissa, about 2000 times double's
// precision): closed forms where they do not cancel, long series where they
// do.
using Long = long double;

Long sincReference(Long x) { return exactlyEqual(x, 0.L) ? 1.L : std::sin(x) / x; }

Long dSincReference(Long x)
{
  if (std::abs(x) < 0.25L) {
    // -x/3 sum_n (-x^2)^n / prod 2k (2k + 3)
    Long term = 1.L;
    Long sum = 1.L;
    for (std::int32_t k = 1; k < 40; ++k) {
      term *= -(x * x) / (2.L * k * ((2.L * k) + 3.L));
      sum += term;
    }
    return -(x / 3.L) * sum;
  }
  return (std::cos(x) - (std::sin(x) / x)) / x;
}

Long versineReference(Long x)
{
  const Long h = sincReference(x / 2.L);
  return 0.5L * x * h * h;
}

Long dVersineReference(Long x)
{
  const Long h = sincReference(x / 2.L);
  return h * (std::cos(x / 2.L) - (0.5L * h));
}

Long asinOverArgumentReference(Long x) { return exactlyEqual(x, 0.L) ? 1.L : std::asin(x) / x; }

Long dAsinOverArgumentReference(Long x)
{
  if (std::abs(x) < 0.25L) {
    // sum_n 2n a_n x^{2n-1}, a_n = C(2n, n) / (4^n (2n + 1))
    Long a = 1.L;
    Long sum = 0.L;
    Long power = x;
    for (std::int32_t n = 1; n < 60; ++n) {
      a *= ((2.L * n) - 1.L) * ((2.L * n) - 1.L) / ((2.L * n) * ((2.L * n) + 1.L));
      sum += 2.L * n * a * power;
      power *= x * x;
    }
    return sum;
  }
  return ((x / std::sqrt((1.L - x) * (1.L + x))) - std::asin(x)) / (x * x);
}

Long atanOverArgumentReference(Long x) { return exactlyEqual(x, 0.L) ? 1.L : std::atan(x) / x; }

// The arguments: log-spaced from 1e-30 to the end, both signs, and a
// thousand around dSinc's switch.
template <typename T>
std::vector<T> argumentsUpTo(double end)
{
  std::vector<T> xs;
  const auto steps = static_cast<std::int32_t>(std::floor((std::log10(end) + 30.) / 0.005));
  for (std::int32_t step = 0; step <= steps; ++step) {
    const T x = static_cast<T>(std::pow(10., -30. + (step * 0.005)));
    if (x > T(0) && static_cast<double>(x) <= end) {
      xs.push_back(x);
      xs.push_back(-x);
    }
  }
  const auto at = static_cast<double>(detail::DSincSeries<T>::kSwitch);
  for (std::int32_t i = -500; i <= 500; ++i) {
    xs.push_back(static_cast<T>(at * (1. + (i * 1e-4))));
  }
  return xs;
}

// The worst relative error, in units of T's epsilon, of f against reference.
template <typename T, typename F, typename R>
double worstError(const std::vector<T>& xs, F f, R reference)
{
  double worst = 0.;
  for (const T x : xs) {
    const Long truth = reference(static_cast<Long>(x));
    if (exactlyEqual(truth, 0.L)) {
      continue;
    }
    const double error = static_cast<double>(std::abs((static_cast<Long>(f(x)) - truth) / truth)) / static_cast<double>(std::numeric_limits<T>::epsilon());
    worst = std::max(worst, error);
  }
  return worst;
}

template <typename T>
void checkAccuracy(const char* type)
{
  // sinc and its relatives up to x = 2 (dSinc's first zero is beyond 4); the
  // asin family up to 1 - 1/16 here, its derivative's 1 / r growing
  // towards 1 (checked apart below).
  const auto xs = argumentsUpTo<T>(2.);
  const auto inside = argumentsUpTo<T>(0.9375);
  const double sincError = worstError(xs, [](T x) { return sinc(x); }, sincReference);
  const double versineError = worstError(xs, [](T x) { return versine(x); }, versineReference);
  const double sincVersineError = std::max(worstError(xs, [](T x) { T s; T v; sincVersine(x, s, v); return s; }, sincReference), worstError(xs, [](T x) { T s; T v; sincVersine(x, s, v); return v; }, versineReference));
  const double dSincError = worstError(xs, [](T x) { return dSinc(x); }, dSincReference);
  const double dVersineError = worstError(xs, [](T x) { return dVersine(x); }, dVersineReference);
  const double asinError = worstError(inside, [](T x) { return detail::asinOverArgument(x); }, asinOverArgumentReference);
  const double dAsinError = worstError(inside, [](T x) { return detail::asinOverArgumentWithDerivative(x).derivative; }, dAsinOverArgumentReference);
  const double atanError = worstError(xs, [](T x) { return detail::atanOverArgument(x); }, atanOverArgumentReference);
  BOOST_TEST_MESSAGE(type << ": dSinc " << detail::DSincSeries<T>::kTerms << " terms below " << static_cast<double>(detail::DSincSeries<T>::kSwitch) << "; worst errors [eps]: sinc " << sincError << ", versine " << versineError << ", sincVersine " << sincVersineError << ", dSinc " << dSincError << ", dVersine " << dVersineError << ", asin(x)/x " << asinError << ", its derivative " << dAsinError << ", atan(x)/x " << atanError);
  // The values: a transcendental and a quotient, a few eps.
  BOOST_TEST(sincError <= 4.);
  BOOST_TEST(versineError <= 4.);
  BOOST_TEST(sincVersineError <= 4.);
  BOOST_TEST(asinError <= 4.);
  BOOST_TEST(atanError <= 4.);
  // dSinc within its target; dVersine and asin's derivative add a value's
  // few eps to it (asin's: A^2 and the root).
  BOOST_TEST(dSincError <= detail::kDSincTarget);
  BOOST_TEST(dVersineError <= 4.);
  BOOST_TEST(dAsinError <= detail::kDSincTarget + 8.);
  // Towards 1, A' = -sinc'(theta) A^2 / r: r's factored form keeps its
  // digits, so the error stays bounded where 1 - x^2 would lose them.
  std::vector<T> nearOne;
  for (std::int32_t k = 1; k < std::numeric_limits<T>::digits; ++k) {
    nearOne.push_back(T(1) - std::ldexp(T(1), -k));
  }
  const double nearOneError = worstError(nearOne, [](T x) { return detail::asinOverArgumentWithDerivative(x).derivative; }, dAsinOverArgumentReference);
  BOOST_TEST_MESSAGE(type << ": asin(x)/x's derivative towards 1: worst " << nearOneError << " eps");
  BOOST_TEST(nearOneError <= detail::kDSincTarget + 8.);
}

} // namespace

BOOST_AUTO_TEST_CASE(TheSincFamilyIsAccurate)
{
  checkAccuracy<float>("float");
  checkAccuracy<double>("double");
}

BOOST_AUTO_TEST_CASE(DualSincAndAsinChainTheScalarOnes)
{
  using D = dual::Dual<1>;
  for (const double x : {0., 1e-9, 1e-3, 0.3, 0.8, 1.7}) {
    const auto s = sinc(dual::variable<1>(x, 0));
    BOOST_TEST(s.v == sinc(x));
    BOOST_TEST(s.d[0] == dSinc(x));
  }
  for (const double x : {0., 1e-9, 1e-3, 0.3, 0.8, 0.99}) {
    const D a = detail::asinOverArgument(dual::variable<1>(x, 0));
    const auto both = detail::asinOverArgumentWithDerivative(x);
    BOOST_TEST(a.v == both.value);
    BOOST_TEST(a.d[0] == both.derivative);
  }
}

} // namespace o2::gtrk
