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

/// \file Dual.h
/// \brief Forward differentiation: a value carried with its derivatives by N
///        inputs (a dual number), the arithmetic and the functions on it.

#ifndef GTRK_STUDY_DUAL_H_
#define GTRK_STUDY_DUAL_H_

#include <array>
#include <cmath>
#include <limits>

namespace o2::gtrk::study::dual
{

/// A value and its derivatives by N inputs: forward differentiation.
template <int N>
struct Dual {
  double v = 0.;
  std::array<double, N> d{};
};

template <int N>
Dual<N> operator+(const Dual<N>& a, const Dual<N>& b)
{
  Dual<N> r{a.v + b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a.d[i] + b.d[i];
  }
  return r;
}
template <int N>
Dual<N> operator-(const Dual<N>& a, const Dual<N>& b)
{
  Dual<N> r{a.v - b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a.d[i] - b.d[i];
  }
  return r;
}
template <int N>
Dual<N> operator-(const Dual<N>& a)
{
  Dual<N> r{-a.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = -a.d[i];
  }
  return r;
}
template <int N>
Dual<N> operator*(const Dual<N>& a, const Dual<N>& b)
{
  Dual<N> r{a.v * b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a.d[i] * b.v + a.v * b.d[i];
  }
  return r;
}
template <int N>
Dual<N> operator/(const Dual<N>& a, const Dual<N>& b)
{
  Dual<N> r{a.v / b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = (a.d[i] - r.v * b.d[i]) / b.v;
  }
  return r;
}
template <int N>
Dual<N> operator*(double a, const Dual<N>& b)
{
  Dual<N> r{a * b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a * b.d[i];
  }
  return r;
}
template <int N>
Dual<N> operator+(const Dual<N>& a, double b)
{
  Dual<N> r = a;
  r.v += b;
  return r;
}
/// f(a) with f'(a) given.
template <int N>
Dual<N> chain(const Dual<N>& a, double value, double derivative)
{
  Dual<N> r{value};
  for (int i = 0; i < N; ++i) {
    r.d[i] = derivative * a.d[i];
  }
  return r;
}
template <int N>
Dual<N> sin(const Dual<N>& a)
{
  return chain(a, std::sin(a.v), std::cos(a.v));
}
template <int N>
Dual<N> cos(const Dual<N>& a)
{
  return chain(a, std::cos(a.v), -std::sin(a.v));
}
template <int N>
Dual<N> sqrt(const Dual<N>& a)
{
  const double root = std::sqrt(a.v);
  return chain(a, root, 0.5 / root);
}
/// Below these arguments sinc and asinc take their series, through the x^6
/// term: there its first omitted term is below double's rounding, where the
/// closed forms' derivatives lose digits to cancellation. sinc's omitted term
/// is x^8 / 9!, asinc's 35 x^8 / 1152.
inline const double kSincSeries = std::pow(362880. * std::numeric_limits<double>::epsilon(), 1. / 8.);
inline const double kAsincSeries = std::pow((1152. / 35.) * std::numeric_limits<double>::epsilon(), 1. / 8.);

/// sin(x) / x, 1 at 0.
template <int N>
Dual<N> sinc(const Dual<N>& a)
{
  const double x = a.v;
  if (std::abs(x) < kSincSeries) {
    const double x2 = x * x;
    return chain(a, 1. - (x2 / 6.) + (x2 * x2 / 120.) - (x2 * x2 * x2 / 5040.), x * (-(1. / 3.) + (x2 / 30.) - (x2 * x2 / 840.)));
  }
  return chain(a, std::sin(x) / x, (std::cos(x) / x) - (std::sin(x) / (x * x)));
}
/// asin(x) / x, 1 at 0.
template <int N>
Dual<N> asinc(const Dual<N>& a)
{
  const double x = a.v;
  if (std::abs(x) < kAsincSeries) {
    const double x2 = x * x;
    return chain(a, 1. + (x2 / 6.) + (3. * x2 * x2 / 40.) + (5. * x2 * x2 * x2 / 112.), x * ((1. / 3.) + (3. * x2 / 10.) + (15. * x2 * x2 / 56.)));
  }
  const double root = std::sqrt(1. - (x * x));
  return chain(a, std::asin(x) / x, (1. / (x * root)) - (std::asin(x) / (x * x)));
}
template <int N>
Dual<N> atan2(const Dual<N>& y, const Dual<N>& x)
{
  Dual<N> r{std::atan2(y.v, x.v)};
  const double norm = (x.v * x.v) + (y.v * y.v);
  for (int i = 0; i < N; ++i) {
    r.d[i] = ((x.v * y.d[i]) - (y.v * x.d[i])) / norm;
  }
  return r;
}

} // namespace o2::gtrk::study::dual

#endif // GTRK_STUDY_DUAL_H_
