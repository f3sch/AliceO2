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

#include "GenericTracking/Base/ElementaryMath.h"

#include "GPUCommonDef.h"

#include <array>
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
GPUhdi() Dual<N, T> operator+(const Dual<N, T>& a, const Dual<N, T>& b)
{
  Dual<N, T> r{a.v + b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a.d[i] + b.d[i];
  }
  return r;
}
template <int N, typename T>
GPUhdi() Dual<N, T> operator-(const Dual<N, T>& a, const Dual<N, T>& b)
{
  Dual<N, T> r{a.v - b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a.d[i] - b.d[i];
  }
  return r;
}
template <int N, typename T>
GPUhdi() Dual<N, T> operator-(const Dual<N, T>& a)
{
  Dual<N, T> r{-a.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = -a.d[i];
  }
  return r;
}
template <int N, typename T>
GPUhdi() Dual<N, T> operator*(const Dual<N, T>& a, const Dual<N, T>& b)
{
  Dual<N, T> r{a.v * b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a.d[i] * b.v + a.v * b.d[i];
  }
  return r;
}
template <int N, typename T>
GPUhdi() Dual<N, T> operator/(const Dual<N, T>& a, const Dual<N, T>& b)
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
GPUhdi() Dual<N, T> operator*(std::type_identity_t<T> a, const Dual<N, T>& b)
{
  Dual<N, T> r{a * b.v};
  for (int i = 0; i < N; ++i) {
    r.d[i] = a * b.d[i];
  }
  return r;
}
template <int N, typename T>
GPUhdi() Dual<N, T> operator+(const Dual<N, T>& a, std::type_identity_t<T> b)
{
  Dual<N, T> r = a;
  r.v += b;
  return r;
}
/// f(a) with f'(a) given.
template <int N, typename T>
GPUhdi() Dual<N, T> chain(const Dual<N, T>& a, T value, T derivative)
{
  Dual<N, T> r{value};
  for (int i = 0; i < N; ++i) {
    r.d[i] = derivative * a.d[i];
  }
  return r;
}
template <int N, typename T>
GPUhdi() Dual<N, T> sin(const Dual<N, T>& a)
{
  return chain(a, sinOf(a.v), cosOf(a.v));
}
template <int N, typename T>
GPUhdi() Dual<N, T> cos(const Dual<N, T>& a)
{
  return chain(a, cosOf(a.v), -sinOf(a.v));
}
template <int N, typename T>
GPUhdi() Dual<N, T> sqrt(const Dual<N, T>& a)
{
  const T root = sqrtOf(a.v);
  return chain(a, root, T(0.5) / root);
}
template <int N, typename T>
GPUhdi() Dual<N, T> asin(const Dual<N, T>& a)
{
  return chain(a, asinOf(a.v), T(1.) / sqrtOf(T(1.) - (a.v * a.v)));
}
template <int N, typename T>
GPUhdi() Dual<N, T> atan2(const Dual<N, T>& y, const Dual<N, T>& x)
{
  Dual<N, T> r{atan2Of(y.v, x.v)};
  const T norm = (x.v * x.v) + (y.v * y.v);
  for (int i = 0; i < N; ++i) {
    r.d[i] = ((x.v * y.d[i]) - (y.v * x.d[i])) / norm;
  }
  return r;
}

/// An input: \a value (taken in T), its derivative by itself (input
/// \a index) one.
template <int N, typename T = double>
[[nodiscard]] GPUhdi() Dual<N, T> variable(std::type_identity_t<T> value, int index)
{
  Dual<N, T> r{value};
  r.d[index] = T(1.);
  return r;
}

/// A value, of a dual number or of a scalar.
template <int N, typename T>
[[nodiscard]] GPUhdi() T valueOf(const Dual<N, T>& a)
{
  return a.v;
}
[[nodiscard]] GPUhdi() double valueOf(double a) { return a; }
[[nodiscard]] GPUhdi() float valueOf(float a) { return a; }

} // namespace o2::gtrk::dual

#endif // GTRK_BASE_DUAL_H_
