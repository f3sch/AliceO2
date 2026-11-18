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

#ifndef GTRK_BASE_DUAL_H_
#define GTRK_BASE_DUAL_H_

#include <array>
#include <cmath>
#include <limits>
#include <type_traits>

namespace o2::gtrk::dual
{

/// A value and its derivatives by N inputs: forward differentiation, in
/// the arithmetic T (double or float).
template <int N, typename T = double>
struct Dual {
  static constexpr int kInputs = N;
  using Value = T;
  T v = 0.;
  std::array<T, N> d{};
};

/// A scalar or a dual number's arithmetic.
template <typename T>
struct ScalarOf {
  using type = T;
};
template <int N, typename T>
struct ScalarOf<Dual<N, T>> {
  using type = T;
};

template <int N, typename T>
Dual<N, T> operator+(const Dual<N, T>& a, const Dual<N, T>& b)
{
  Dual<N, T> r{a.v + b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a.d[i] + b.d[i];
  }
  return r;
}
template <int N, typename T>
Dual<N, T> operator-(const Dual<N, T>& a, const Dual<N, T>& b)
{
  Dual<N, T> r{a.v - b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a.d[i] - b.d[i];
  }
  return r;
}
template <int N, typename T>
Dual<N, T> operator-(const Dual<N, T>& a)
{
  Dual<N, T> r{-a.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = -a.d[i];
  }
  return r;
}
template <int N, typename T>
Dual<N, T> operator*(const Dual<N, T>& a, const Dual<N, T>& b)
{
  Dual<N, T> r{a.v * b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a.d[i] * b.v + a.v * b.d[i];
  }
  return r;
}
template <int N, typename T>
Dual<N, T> operator/(const Dual<N, T>& a, const Dual<N, T>& b)
{
  Dual<N, T> r{a.v / b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = (a.d[i] - r.v * b.d[i]) / b.v;
  }
  return r;
}
/// A constant times a dual number: the constant taken in the dual's
/// arithmetic.
template <int N, typename T>
Dual<N, T> operator*(std::type_identity_t<T> a, const Dual<N, T>& b)
{
  Dual<N, T> r{a * b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a * b.d[i];
  }
  return r;
}
template <int N, typename T>
Dual<N, T> operator+(const Dual<N, T>& a, std::type_identity_t<T> b)
{
  Dual<N, T> r = a;
  r.v += b;
  return r;
}
/// f(a) with f'(a) given.
template <int N, typename T>
Dual<N, T> chain(const Dual<N, T>& a, T value, T derivative)
{
  Dual<N, T> r{value};
  for (int i = 0; i < N; ++i) {
    r.d[i] = derivative * a.d[i];
  }
  return r;
}
template <int N, typename T>
Dual<N, T> sin(const Dual<N, T>& a)
{
  return chain(a, std::sin(a.v), std::cos(a.v));
}
template <int N, typename T>
Dual<N, T> cos(const Dual<N, T>& a)
{
  return chain(a, std::cos(a.v), -std::sin(a.v));
}
template <int N, typename T>
Dual<N, T> sqrt(const Dual<N, T>& a)
{
  const T root = std::sqrt(a.v);
  return chain(a, root, T(0.5) / root);
}
/// Below these arguments sinc and asinc take their series, through the x^6
/// term: there its first omitted term is below T's rounding, where the
/// closed forms' derivatives lose digits to cancellation. sinc's omitted term
/// is x^8 / 9!, asinc's 35 x^8 / 1152.
template <typename T>
inline const T kSincSeries = std::pow(T(362880.) * std::numeric_limits<T>::epsilon(), T(1.) / T(8.));
template <typename T>
inline const T kAsincSeries = std::pow((T(1152.) / T(35.)) * std::numeric_limits<T>::epsilon(), T(1.) / T(8.));

/// sin(x) / x, 1 at 0.
template <int N, typename T>
Dual<N, T> sinc(const Dual<N, T>& a)
{
  const T x = a.v;
  if (std::abs(x) < kSincSeries<T>) {
    const T x2 = x * x;
    return chain(a, T(1.) - (x2 / T(6.)) + (x2 * x2 / T(120.)) - (x2 * x2 * x2 / T(5040.)), x * (-(T(1.) / T(3.)) + (x2 / T(30.)) - (x2 * x2 / T(840.))));
  }
  return chain(a, std::sin(x) / x, (std::cos(x) / x) - (std::sin(x) / (x * x)));
}
/// asin(x) / x, 1 at 0.
template <int N, typename T>
Dual<N, T> asinc(const Dual<N, T>& a)
{
  const T x = a.v;
  if (std::abs(x) < kAsincSeries<T>) {
    const T x2 = x * x;
    return chain(a, T(1.) + (x2 / T(6.)) + (T(3.) * x2 * x2 / T(40.)) + (T(5.) * x2 * x2 * x2 / T(112.)), x * ((T(1.) / T(3.)) + (T(3.) * x2 / T(10.)) + (T(15.) * x2 * x2 / T(56.))));
  }
  const T root = std::sqrt(T(1.) - (x * x));
  return chain(a, std::asin(x) / x, (T(1.) / (x * root)) - (std::asin(x) / (x * x)));
}
template <int N, typename T>
Dual<N, T> atan2(const Dual<N, T>& y, const Dual<N, T>& x)
{
  Dual<N, T> r{std::atan2(y.v, x.v)};
  const T norm = (x.v * x.v) + (y.v * y.v);
  for (int i = 0; i < N; ++i) {
    r.d[i] = ((x.v * y.d[i]) - (y.v * x.d[i])) / norm;
  }
  return r;
}

/// An input: \a value (taken in T), its derivative by itself (input
/// \a index) one.
template <int N, typename T = double>
[[nodiscard]] Dual<N, T> variable(std::type_identity_t<T> value, int index)
{
  Dual<N, T> r{value};
  r.d[index] = T(1.);
  return r;
}

/// A value, of a dual number or of a scalar.
template <int N, typename T>
[[nodiscard]] T valueOf(const Dual<N, T>& a)
{
  return a.v;
}
[[nodiscard]] inline double valueOf(double a) { return a; }
[[nodiscard]] inline float valueOf(float a) { return a; }

} // namespace o2::gtrk::dual

#endif // GTRK_BASE_DUAL_H_
