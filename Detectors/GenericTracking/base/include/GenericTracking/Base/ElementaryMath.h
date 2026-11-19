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

/// \file ElementaryMath.h
/// \brief The elementary functions of shared code, overloaded for float and
///        double, for host and device alike, and the finiteness test.
///
/// Outside O2's deterministic mode they are `CAMath`'s: `Sin`, `ATan2`,
/// `SinCos`, ... in float and `Sind`, `ATan2d`, `SinCosd`, ... in double.
///
/// O2's deterministic mode is different. It computes the float functions in
/// double and rounds, which almost always hides the host's and the device's
/// libraries disagreeing in the last place of the double; the double
/// functions it leaves native, and there nothing hides it (PortableMath.h). So
/// in that mode every function here is PortableMath's, in both precisions: the
/// same arithmetic on the host and on every device, the float ones in double
/// and rounded, as O2 does. The square root is correctly rounded everywhere
/// (IEEE) and needs no variant.

#ifndef GTRK_BASE_ELEMENTARYMATH_H_
#define GTRK_BASE_ELEMENTARYMATH_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/PortableMath.h"

#include "GPUCommonMath.h"

#include <cfloat>
#include <cstdint>
#include <limits>
#include <numbers>

// The documented libdevice float implementation, called directly so CUDA's
// --use-fast-math does not substitute the lower-accuracy intrinsic.
#if defined(__CUDACC__) && !defined(__HIPCC__)
extern "C" __device__ void __nv_sincosf(float, float*, float*);
#endif

namespace o2::gtrk
{

/// pi and 2 pi in T's own precision, the one definition for GenericTracking:
/// `kPi<Scalar>` for the tracking's float, `kPi<double>` for the exact
/// solves, the fits and the studies. O2's `constants::math::PI` is a float:
/// in double it is 9e-8 off pi.
template <typename T>
inline constexpr T kPi = std::numbers::pi_v<T>;
template <typename T>
inline constexpr T kTwoPi = T(2) * std::numbers::pi_v<T>;

GPUhdi() GPUAlwaysInline() float sinOf(float x)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  return static_cast<float>(portable::sin(static_cast<double>(x)));
#else
  return o2::gpu::CAMath::Sin(x);
#endif
}

GPUhdi() GPUAlwaysInline() float cosOf(float x)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  return static_cast<float>(portable::cos(static_cast<double>(x)));
#else
  return o2::gpu::CAMath::Cos(x);
#endif
}

GPUhdi() GPUAlwaysInline() float asinOf(float x)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  return static_cast<float>(portable::asin(static_cast<double>(x)));
#else
  return o2::gpu::CAMath::ASin(x);
#endif
}

GPUhdi() GPUAlwaysInline() float atanOf(float x)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  return static_cast<float>(portable::atan(static_cast<double>(x)));
#else
  return o2::gpu::CAMath::ATan(x);
#endif
}

GPUhdi() GPUAlwaysInline() float logOf(float x)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  return static_cast<float>(portable::log(static_cast<double>(x)));
#else
  return o2::gpu::CAMath::Log(x);
#endif
}

GPUhdi() GPUAlwaysInline() float atan2Of(float y, float x)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  return static_cast<float>(portable::atan2(static_cast<double>(y), static_cast<double>(x)));
#else
  return o2::gpu::CAMath::ATan2(y, x);
#endif
}

#ifdef GPUCA_DETERMINISTIC_MODE
GPUhdi() GPUAlwaysInline() double sinOf(double x) { return portable::sin(x); }
GPUhdi() GPUAlwaysInline() double cosOf(double x) { return portable::cos(x); }
GPUhdi() GPUAlwaysInline() double asinOf(double x) { return portable::asin(x); }
GPUhdi() GPUAlwaysInline() double atanOf(double x) { return portable::atan(x); }
GPUhdi() GPUAlwaysInline() double atan2Of(double y, double x) { return portable::atan2(y, x); }
GPUhdi() GPUAlwaysInline() double logOf(double x) { return portable::log(x); }
#else
GPUhdi() GPUAlwaysInline() double sinOf(double x) { return o2::gpu::CAMath::Sind(x); }
GPUhdi() GPUAlwaysInline() double cosOf(double x) { return o2::gpu::CAMath::Cosd(x); }
GPUhdi() GPUAlwaysInline() double asinOf(double x) { return o2::gpu::CAMath::ASind(x); }
GPUhdi() GPUAlwaysInline() double atanOf(double x) { return o2::gpu::CAMath::ATand(x); }
GPUhdi() GPUAlwaysInline() double atan2Of(double y, double x) { return o2::gpu::CAMath::ATan2d(y, x); }
GPUhdi() GPUAlwaysInline() double logOf(double x) { return o2::gpu::CAMath::Logd(x); }
#endif

GPUhdi() GPUAlwaysInline() float sqrtOf(float x) { return o2::gpu::CAMath::Sqrt(x); }
GPUhdi() GPUAlwaysInline() double sqrtOf(double x) { return o2::gpu::CAMath::Sqrtd(x); }

/// The length of a vector of two or three components. CAMath's `Hypot` is
/// float's only; double's is its sum of squares under `Sqrtd`, the same form.
GPUhdi() GPUAlwaysInline() float hypotOf(float a, float b) { return o2::gpu::CAMath::Hypot(a, b); }
GPUhdi() GPUAlwaysInline() float hypotOf(float a, float b, float c) { return o2::gpu::CAMath::Hypot(a, b, c); }
GPUhdi() GPUAlwaysInline() double hypotOf(double a, double b) { return sqrtOf((a * a) + (b * b)); }
GPUhdi() GPUAlwaysInline() double hypotOf(double a, double b, double c) { return sqrtOf((a * a) + (b * b) + (c * c)); }

/// ln(1 + x) without losing x's digits to the sum: ln(u) x / (u - 1) with
/// u = 1 + x rounded, exact where u is 1 (Goldberg's form). CAMath has no
/// log1p; `logOf` keeps the deterministic mode's.
GPUhdi() GPUAlwaysInline() double log1pOf(double x)
{
  const double u = 1. + x;
  return u == 1. ? x : logOf(u) * x / (u - 1.);
}

GPUhdi() GPUAlwaysInline() void sincosOf(float x, float& s, float& c)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  s = sinOf(x);
  c = cosOf(x);
#else
  o2::gpu::CAMath::SinCos(x, s, c);
#endif
}

GPUhdi() GPUAlwaysInline() void sincosOf(double x, double& s, double& c)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  s = portable::sin(x);
  c = portable::cos(x);
#else
  o2::gpu::CAMath::SinCosd(x, s, c);
#endif
}

/// sin and cos of a direction for the double plane solve. Under CUDA's
/// -use_fast_math (outside the deterministic mode) `sincosf` becomes the
/// intrinsic, too coarse here; libdevice's precise one is called instead.
GPUhdi() GPUAlwaysInline() void directionSincosOf(float angle, float& sine, float& cosine)
{
#if defined(__CUDA_ARCH__) && !defined(GPUCA_DETERMINISTIC_MODE)
  __nv_sincosf(angle, &sine, &cosine);
#else
  sincosOf(angle, sine, cosine);
#endif
}

GPUhdi() GPUAlwaysInline() void directionSincosOf(float angle, double& sine, double& cosine)
{
  float s = 0.F;
  float c = 0.F;
#if defined(__CUDA_ARCH__) && !defined(GPUCA_DETERMINISTIC_MODE)
  __nv_sincosf(angle, &s, &c);
#else
  sincosOf(angle, s, c);
#endif
  sine = s;
  cosine = c;
}

/// An angle brought to [-pi, pi] by whole turns, any number of them, with pi
/// in T's own precision (as std::remainder by 2 pi): the nearest whole number
/// of turns through an integer, exact while it is below 2^(digits - 1); an
/// angle beyond that (or not finite) has no meaningful turn and is returned
/// as it is. O2's `toPMPi` takes one turn and float's pi; this one any, for
/// an angle a fit's unbounded steps may carry turns away.
template <typename T>
[[nodiscard]] GPUhdi() T toPMPiOf(T x)
{
  constexpr T kExact = static_cast<T>(std::uint64_t{1} << (std::numeric_limits<T>::digits - 1));
  const T turns = x / kTwoPi<T>;
  if (!(turns < kExact && turns > -kExact)) {
    return x;
  }
  const T whole = static_cast<T>(static_cast<std::int64_t>(turns + (turns >= T(0) ? T(0.5) : T(-0.5))));
  T reduced = x - (kTwoPi<T> * whole);
  // the rounded quotient can leave the result a hair outside
  if (reduced > kPi<T>) {
    reduced -= kTwoPi<T>;
  } else if (reduced < -kPi<T>) {
    reduced += kTwoPi<T>;
  }
  return reduced;
}

/// Whether \a x is finite, in every build: O2's `GPUCommonMath::Finite`
/// and `IsNaN` are constants outside the deterministic mode, on the host as
/// on the device. NaN fails both comparisons, an infinity one of them. The
/// bounds are the C macros, which a device compiler accepts where a
/// `std::numeric_limits` call may not.
GPUhdi() GPUAlwaysInline() constexpr bool isFiniteOf(float x) { return x >= -FLT_MAX && x <= FLT_MAX; }
GPUhdi() GPUAlwaysInline() constexpr bool isFiniteOf(double x) { return x >= -DBL_MAX && x <= DBL_MAX; }

} // namespace o2::gtrk

#endif // GTRK_BASE_ELEMENTARYMATH_H_
