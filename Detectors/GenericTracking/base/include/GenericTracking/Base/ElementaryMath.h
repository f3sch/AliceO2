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
///        double, for host and device alike.
///
/// Outside O2's deterministic mode they are `CAMath`'s: `Sin`, `ATan2`, ... in
/// float and `Sind`, `ATan2d`, ... in double.
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

namespace o2::gtrk
{

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

} // namespace o2::gtrk

#endif // GTRK_BASE_ELEMENTARYMATH_H_
