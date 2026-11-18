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

/// \file Trig.h
/// \brief The elementary functions of shared tracking code (ElementaryMath.h)
///        and the sincos that needs O2's MathUtils: `sincosOf` is
///        `o2::math_utils::detail::sincos`, one call on hosts and devices, and
///        the spelling O2's own track code uses (in the deterministic mode
///        PortableMath's sin and cos).

#ifndef GTRK_TRACKING_TRIG_H_
#define GTRK_TRACKING_TRIG_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Base/PortableMath.h"

#include "GPUCommonMath.h"
#include "MathUtils/detail/trigonometric.h"

// The documented libdevice float implementation, called directly so CUDA's
// --use-fast-math does not substitute the lower-accuracy intrinsic.
#if defined(__CUDACC__) && !defined(__HIPCC__)
extern "C" __device__ void __nv_sincosf(float, float*, float*);
#endif

namespace o2::gtrk
{

GPUhdi() GPUAlwaysInline() void sincosOf(float x, float& s, float& c)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  s = sinOf(x);
  c = cosOf(x);
#else
  math_utils::detail::sincos<float>(x, s, c);
#endif
}

GPUhdi() GPUAlwaysInline() void sincosOf(double x, double& s, double& c)
{
#ifdef GPUCA_DETERMINISTIC_MODE
  s = portable::sin(x);
  c = portable::cos(x);
#else
  math_utils::detail::sincos<double>(x, s, c);
#endif
}

/// sin and cos of a direction for the double plane solve. Under CUDA's
/// -use_fast_math (outside the deterministic mode) `sincosf` becomes the
/// intrinsic, too coarse here; libdevice's precise one is called instead.
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

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRIG_H_
