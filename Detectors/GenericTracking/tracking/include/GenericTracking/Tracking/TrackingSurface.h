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

/// \file TrackingSurface.h
/// \brief The sensor planes as the tracking sees them, and the conversion
///        between a bound position and the global frame.
///
/// This is the device side of the host-to-tracking boundary: the double
/// Geometry<Det> narrows here, once, into fixed-size float records;
/// Geometry<Det> is the host side.
///
/// A surface's three axes are the **columns** of Sensor::rot, which is stored
/// row-major and maps local to global.  Writing the rows instead is a
/// transpose that a same-surface round trip cannot see, because it would use
/// the same wrong matrix in both directions.  testTrackingSurface.cxx
/// therefore checks against Geometry<Det>::toGlobal, an independent
/// implementation, and not against this file's own inverse.

#ifndef GTRK_TRACKING_TRACKINGSURFACE_H_
#define GTRK_TRACKING_TRACKINGSURFACE_H_

#include "GenericTracking/Base/Definitions.h"
#include "GenericTracking/Base/Dual.h"
#include "GenericTracking/Base/ElementaryMath.h"
#include "GenericTracking/Base/Geometry.h"
#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Tracking/BoundTrackState.h"
#include "GenericTracking/Tracking/Definitions.h"
#include "GenericTracking/Tracking/LinkMaterial.h"

#include "GPUCommonMath.h"

#include <array>
#include <cstdint>
#include <type_traits>

namespace o2::gtrk
{

/// The active area of a sensor, a rectangle in its own local frame.
struct SurfaceBounds {
  Scalar halfX = 0.F; ///< half width along loc0
  Scalar halfY = 0.F; ///< half length along loc1

  /// Written so that a NaN coordinate is outside.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool isInside(Scalar loc0, Scalar loc1) const { return loc1 >= -halfY && loc1 <= halfY && loc0 >= -halfX && loc0 <= halfX; }
  /// Inside the rectangle grown by \a margin on every side; NaN is outside.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() bool isInside(Scalar loc0, Scalar loc1, Scalar margin) const { return loc1 >= -(halfY + margin) && loc1 <= halfY + margin && loc0 >= -(halfX + margin) && loc0 <= halfX + margin; }
};

/// One sensor plane.
///
/// Held as a record rather than as columns of a struct-of-arrays because a
/// propagation reads a whole surface at once -- origin, both in-plane axes and
/// the normal -- so the fifteen floats travel together.
struct TrackingSurface {
  Scalar origin[3]{}; ///< sensor centre in the global frame [cm]
  Scalar axis0[3]{};  ///< global direction of local loc0
  Scalar axis1[3]{};  ///< global direction of local loc1
  Scalar normal[3]{}; ///< global plane normal, axis0 x axis1

  SurfaceBounds bounds{};

  /// The surface's own chip id.  A propagation is handed the surface and not a
  /// geometry to look it up in, so without this the resulting state could not
  /// say which plane it is bound to.
  SensorID id = kInvalidSensorID;

  LayerID layer = kInvalidLayerID;

  /// A curved sensor's own entry is its cylinder, not a plane: the origin a
  /// point of its axis, axis0 and axis1 the sensor's local x and y, the normal
  /// its axis, the bounds unrolled, and this the radius its unrolled x is
  /// measured at (SensorSpec::radius). 0 on a flat sensor, and on the planes
  /// tangentPlane() makes. [cm]
  Scalar radius = 0.F;
  /// On a plane tangentPlane() made: its angle in the sensor, which a state
  /// bound to it carries (BoundTrackState::surfacePhi). 0 otherwise. [rad]
  Scalar surfacePhi = 0.F;
  /// A curved sensor's measurement depth (Geometry::measurementDepth): its
  /// tangent planes touch at radius + depth, while its unrolled x stays
  /// measured at the radius. 0 on a flat sensor, whose origin already lies at
  /// its depth, and on the planes tangentPlane() makes. [cm]
  Scalar depth = 0.F;
};

/// The plane touching the curved sensor \a cylinder at its local angle
/// \a angle, at radius + depth: normal
/// radial, axis1 along the cylinder's axis, axis0 = axis1 x normal along
/// increasing angle. Its origin is placed so that the touching point is at
/// loc0 = radius * angle, loc1 = 0 at the axis point: a point of the sensor near there has its unrolled
/// coordinates as its plane coordinates, and a measurement keeps its own. The
/// bounds are the sensor's unrolled ones.
GPUhdi() GPUAlwaysInline() TrackingSurface tangentPlane(const TrackingSurface& cylinder, Scalar angle)
{
  Scalar sinAngle = 0.F;
  Scalar cosAngle = 0.F;
  sincosOf(angle, sinAngle, cosAngle);
  TrackingSurface plane = cylinder;
  const Scalar arc = cylinder.radius * angle;
  for (std::int32_t i = 0; i < 3; ++i) {
    const Scalar radial = cosAngle * cylinder.axis0[i] + sinAngle * cylinder.axis1[i];
    const Scalar along = -sinAngle * cylinder.axis0[i] + cosAngle * cylinder.axis1[i];
    plane.origin[i] = cylinder.origin[i] + (cylinder.radius + cylinder.depth) * radial - arc * along;
    plane.axis0[i] = along;
    plane.axis1[i] = cylinder.normal[i];
    plane.normal[i] = radial;
  }
  plane.radius = 0.F;
  plane.depth = 0.F;
  plane.surfacePhi = angle;
  return plane;
}

/// global = origin + loc0 * axis0 + loc1 * axis1, in the precision T of \a g
/// (Scalar, ExactScalar, or a Dual of the fits).
template <typename T>
GPUhdi() GPUAlwaysInline() void surfaceToGlobal(const TrackingSurface& s, std::type_identity_t<T> loc0, std::type_identity_t<T> loc1, T g[3])
{
  // The surface's floats as plain numbers: a Dual multiplies them as
  // scalars; in float or double the same sum as origin + loc0 a0 + loc1 a1.
  using Plain = typename dual::ScalarOf<T>::type;
  for (std::int32_t i = 0; i < 3; ++i) {
    g[i] = ((static_cast<Plain>(s.axis0[i]) * loc0) + static_cast<Plain>(s.origin[i])) + (static_cast<Plain>(s.axis1[i]) * loc1);
  }
}

/// The inverse, including the out-of-plane component (the height above the
/// plane along its normal).
///
/// l[2] is returned rather than discarded so that a caller can check that a
/// point it believes is on the surface actually is, instead of assuming it.
template <typename T>
GPUhdi() GPUAlwaysInline() void surfaceToLocal(const TrackingSurface& s, const T g[3], T l[3])
{
  // As in surfaceToGlobal, the surface's floats enter as plain numbers.
  using Plain = typename dual::ScalarOf<T>::type;
  const T d0 = g[0] + static_cast<Plain>(-s.origin[0]);
  const T d1 = g[1] + static_cast<Plain>(-s.origin[1]);
  const T d2 = g[2] + static_cast<Plain>(-s.origin[2]);
  l[0] = static_cast<Plain>(s.axis0[0]) * d0 + static_cast<Plain>(s.axis0[1]) * d1 + static_cast<Plain>(s.axis0[2]) * d2;
  l[1] = static_cast<Plain>(s.axis1[0]) * d0 + static_cast<Plain>(s.axis1[1]) * d1 + static_cast<Plain>(s.axis1[2]) * d2;
  l[2] = static_cast<Plain>(s.normal[0]) * d0 + static_cast<Plain>(s.normal[1]) * d1 + static_cast<Plain>(s.normal[2]) * d2;
}

/// \a n dot \a d, the cosine of incidence on a surface.  Its magnitude is what
/// the grazing check thresholds on.
template <typename T>
GPUhdi() GPUAlwaysInline() T cosIncidence(const TrackingSurface& s, const T d[3])
{
  using Plain = typename dual::ScalarOf<T>::type;
  return static_cast<Plain>(s.normal[0]) * d[0] + static_cast<Plain>(s.normal[1]) * d[1] + static_cast<Plain>(s.normal[2]) * d[2];
}

/// Where the straight line through \a point along \a direction crosses plane
/// \a s: the line's parameter \a step (in units of \a direction) and the local
/// position \a loc there. False, and nothing written, for a line parallel to
/// the plane; a caller that needs a margin from grazing checks
/// cosIncidence first.
template <typename T>
GPUhdi() GPUAlwaysInline() bool lineCrossing(const TrackingSurface& s, const T point[3], const T direction[3], T& step, T loc[2])
{
  T normalDirection = T(0);
  T normalOffset = T(0);
  for (std::int32_t i = 0; i < 3; ++i) {
    normalDirection += static_cast<T>(s.normal[i]) * direction[i];
    normalOffset += static_cast<T>(s.normal[i]) * (static_cast<T>(s.origin[i]) - point[i]);
  }
  if (!(normalDirection > T(0) || normalDirection < T(0))) {
    return false;
  }
  step = normalOffset / normalDirection;
  loc[0] = T(0);
  loc[1] = T(0);
  for (std::int32_t i = 0; i < 3; ++i) {
    const T offset = point[i] + (step * direction[i]) - static_cast<T>(s.origin[i]);
    loc[0] += offset * static_cast<T>(s.axis0[i]);
    loc[1] += offset * static_cast<T>(s.axis1[i]);
  }
  return true;
}

/// The plane through \a origin facing the transverse direction (cos phi,
/// sin phi), given as \a sinPhi, \a cosPhi: loc0 along (-sin phi, cos phi),
/// loc1 along z, the frame of a barrel layer at azimuth phi. No sensor, no
/// layer, empty bounds: the caller sets what it needs.
GPUhdi() GPUAlwaysInline() TrackingSurface radialPlane(const Scalar origin[3], Scalar sinPhi, Scalar cosPhi)
{
  TrackingSurface plane;
  for (std::int32_t i = 0; i < 3; ++i) {
    plane.origin[i] = origin[i];
  }
  plane.axis0[0] = -sinPhi;
  plane.axis0[1] = cosPhi;
  plane.axis1[2] = 1.F;
  plane.normal[0] = cosPhi;
  plane.normal[1] = sinPhi;
  return plane;
}

/// The plane through \a origin facing +z: loc0 along the transverse direction
/// (cos phi, sin phi), loc1 along (-sin phi, cos phi), the frame of a disc at
/// azimuth phi. As radialPlane, no sensor, no layer, empty bounds.
GPUhdi() GPUAlwaysInline() TrackingSurface discPlane(const Scalar origin[3], Scalar sinPhi, Scalar cosPhi)
{
  TrackingSurface plane;
  for (std::int32_t i = 0; i < 3; ++i) {
    plane.origin[i] = origin[i];
  }
  plane.axis0[0] = cosPhi;
  plane.axis0[1] = sinPhi;
  plane.axis1[0] = -sinPhi;
  plane.axis1[1] = cosPhi;
  plane.normal[2] = 1.F;
  return plane;
}

/// The angle in curved sensor \a cylinder of global point \a g, as
/// tangentPlane takes it: of g's offset from the axis point, in the plane of
/// axis0 and axis1, from axis0.
GPUhdi() GPUAlwaysInline() Scalar angleOn(const TrackingSurface& cylinder, const Scalar g[3])
{
  Scalar along0 = 0.F;
  Scalar along1 = 0.F;
  for (std::int32_t i = 0; i < 3; ++i) {
    along0 += (g[i] - cylinder.origin[i]) * cylinder.axis0[i];
    along1 += (g[i] - cylinder.origin[i]) * cylinder.axis1[i];
  }
  return atan2Of(along1, along0);
}

/// Every sensor plane of one detector, as the tracking sees them.
///
/// Fixed extent, so the whole thing is one contiguous block sized at compile
/// time and uploads as a single buffer.  Built from Geometry<Det>, which is
/// double and stays on the host; this is where that narrows.
template <DetectorLayout Det>
class TrackingGeometry
{
 public:
  explicit TrackingGeometry(const Geometry<Det>& geo) { build(geo); }

  /// A sensor's own entry: its plane, or a curved sensor's cylinder
  /// (TrackingSurface::radius). Its layer is the sensor's either way.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() const TrackingSurface& surface(SensorID id) const { return mSurfaces[id]; }

  /// What planeOf returns: the table entry itself where every sensor is flat,
  /// a plane built on the spot where some are curved.
  using Plane = std::conditional_t<Detector<Det>::hasCurvedSensors, TrackingSurface, const TrackingSurface&>;

  /// The plane a state is bound to: its sensor's, or on a curved sensor the
  /// tangent plane at the state's angle.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() Plane planeOf(const BoundParameters& state) const { return plane(state.surface, state.surfacePhi); }

  /// The plane a measurement is on: its sensor's, or on a curved sensor the
  /// tangent plane where it sits, at the angle its unrolled loc0 gives.
  template <typename M>
    requires(!std::is_base_of_v<BoundParameters, M>)
  [[nodiscard]] GPUhdi() GPUAlwaysInline() Plane planeOf(const M& measurement) const
  {
    if constexpr (Detector<Det>::hasCurvedSensors) {
      const TrackingSurface& s = mSurfaces[measurement.surface];
      return s.radius > 0.F ? tangentPlane(s, measurement.loc0 / s.radius) : s;
    } else {
      return mSurfaces[measurement.surface];
    }
  }

  /// Sensor \a id's plane, or a curved sensor's tangent plane where global
  /// point \a g lies (angleOn).
  [[nodiscard]] GPUhdi() GPUAlwaysInline() Plane planeAt(SensorID id, [[maybe_unused]] const Scalar g[3]) const
  {
    if constexpr (Detector<Det>::hasCurvedSensors) {
      const TrackingSurface& s = mSurfaces[id];
      return s.radius > 0.F ? tangentPlane(s, angleOn(s, g)) : s;
    } else {
      return mSurfaces[id];
    }
  }

  /// Sensor \a id's plane, or a curved sensor's tangent plane at \a angle.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() Plane plane(SensorID id, [[maybe_unused]] Scalar angle) const
  {
    if constexpr (Detector<Det>::hasCurvedSensors) {
      const TrackingSurface& s = mSurfaces[id];
      return s.radius > 0.F ? tangentPlane(s, angle) : s;
    } else {
      return mSurfaces[id];
    }
  }

  [[nodiscard]] const std::array<TrackingSurface, Det::nSensors>& surfaces() const { return mSurfaces; }

  /// The material of the layers' bands and the links' gaps (LinkMaterial.h),
  /// filled from the LUT once per run and beam line (`fillLinkMaterial`):
  /// the seeding's cuts and the follow's and the fits' steps read it.
  [[nodiscard]] GPUhdi() GPUAlwaysInline() const LinkMaterial<Det>& material() const { return mMaterial; }
  [[nodiscard]] LinkMaterial<Det>& material() { return mMaterial; }

 private:
  void build(const Geometry<Det>& geo);

  std::array<TrackingSurface, Det::nSensors> mSurfaces{};
  LinkMaterial<Det> mMaterial{};
};

template <DetectorLayout Det>
void TrackingGeometry<Det>::build(const Geometry<Det>& geo)
{
  for (SensorID id = 0; id < static_cast<SensorID>(Det::nSensors); ++id) {
    const Sensor& src = geo.sensor(id);
    const SensorSpec& spec = geo.specOf(id);
    TrackingSurface& dst = mSurfaces[id];
    // A flat sensor's plane at its layer's measurement depth, along its normal.
    const double depth = geo.measurementDepth(src.layer);
    const double shift = spec.curved() ? 0. : depth;

    for (std::int32_t i = 0; i < 3; ++i) {
      dst.origin[i] = static_cast<Scalar>(src.trans[i] + (shift * src.rot[3 * i + 2]));
      // Columns, not rows: rot is row-major local -> global, so column j is
      // the global direction of local axis j.
      dst.axis0[i] = static_cast<Scalar>(src.rot[3 * i + 0]);
      dst.axis1[i] = static_cast<Scalar>(src.rot[3 * i + 1]);
      dst.normal[i] = static_cast<Scalar>(src.rot[3 * i + 2]);
    }

    dst.bounds.halfX = static_cast<Scalar>(spec.halfX);
    dst.bounds.halfY = static_cast<Scalar>(spec.halfY);
    dst.radius = spec.curved() ? static_cast<Scalar>(spec.radius) : 0.F;
    dst.depth = spec.curved() ? static_cast<Scalar>(depth) : 0.F;

    dst.id = id;
    dst.layer = src.layer;
  }
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_TRACKINGSURFACE_H_
