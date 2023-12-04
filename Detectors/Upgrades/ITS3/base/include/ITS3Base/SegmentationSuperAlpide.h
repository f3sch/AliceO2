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

/// \file SegmentationSuperAlpide.h
/// \brief Definition of the SegmentationSuperAlpide class
/// \author Fabrizio Grosa <fgrosa@cern.ch>
/// \author felix.schlepper@cern.ch

#ifndef ALICEO2_ITS3_SEGMENTATIONSUPERALPIDE_H_
#define ALICEO2_ITS3_SEGMENTATIONSUPERALPIDE_H_

#include "MathUtils/Cartesian.h"
#include "CommonConstants/MathConstants.h"
#include "ITS3Base/Specs.h"

#include <type_traits>

namespace o2
{
namespace its3
{

/// Segmentation and response for pixels in ITS3 upgrade
class SegmentationSuperAlpide
{
  // This class defines the segmenation of the pixelArray in the tile. We define
  // two coordinate systems, one width x,z detector local coordianates (cm) and
  // the more natural row,col layout: Also all the transformation between these
  // two. The class provides the transformation from the tile to TGeo
  // coordinates.

  // row,col=0
  // |
  // v
  // x----------------------x
  // |           |          |
  // |           |          |
  // |           |          |
  // |           |          |                        ^ x
  // |           |          |                        |
  // |           |          |                        |
  // |           |          |                        |
  // |-----------X----------|  X marks (x,z)=(0,0)   X----> z
  // |           |          |
  // |           |          |
  // |           |          |
  // |           |          |
  // |           |          |
  // |           |          |
  // x----------------------x
 public:
  static constexpr unsigned int mNCols{156};
  static constexpr unsigned int mNRows{440};
  static constexpr unsigned int nPixels{mNCols * mNRows};
  static constexpr double mPitchCol{constants::pixelarray::width / static_cast<double>(mNCols)};
  static constexpr double mPitchRow{constants::pixelarray::length / static_cast<double>(mNRows)};
  static constexpr double mSensorLayerThicknessEff{66 * constants::mu};

 public:
  SegmentationSuperAlpide(int layer = 0) : mLayer{layer} {}

  /// Transformation from the curved surface to a flat surface
  /// It works only if the detector is not rototraslated
  /// \param xCurved Detector local curved coordinate x in cm with respect to
  /// the center of the sensitive volume.
  /// \param yCurved Detector local curved coordinate y in cm with respect to
  /// the center of the sensitive volume.
  /// \param xFlat Detector local flat coordinate x in cm with respect to
  /// the center of the sensitive volume.
  /// \param yFlat Detector local flat coordinate y in cm with respect to
  /// the center of the sensitive volume.
  template <typename T = double>
  void curvedToFlat(T xCurved, T yCurved, T& xFlat, T& yFlat)
  {
    T dist = std::hypot(xCurved, yCurved);
    yFlat = dist - mEffRadius;
    // phi is the angle between the x axis and the center of the pixel array
    T phi = std::atan2(yCurved, xCurved);
    // phiOffset is the angle between the x axis and the center of the pixel array
    T phiOffset = std::asin((constants::pixelarray::length / 2.) / dist);
    T actualPhi = phi - phiOffset;
  }

  /// Transformation from the flat surface to a curved surface
  /// It works only if the detector is not rototraslated
  /// \param xFlat Detector local flat coordinate x in cm with respect to
  /// the center of the sensitive volume.
  /// \param yFlat Detector local flat coordinate y in cm with respect to
  /// the center of the sensitive volume.
  /// \param xCurved Detector local curved coordinate x in cm with respect to
  /// the center of the sensitive volume.
  /// \param yCurved Detector local curved coordinate y in cm with respect to
  /// the center of the sensitive volume.
  template <typename T = double>
  void flatToCurved(T xFlat, T yFlat, T& xCurved, T& yCurved)
  {
    T dist = yFlat + mEffRadius;
    T phi = xFlat / dist;
    // phiOffset is the angle between the x axis and the center of the pixel array
    T phiOffset = std::asin((constants::pixelarray::length / 2.) / dist);
    xCurved = dist * std::cos(phi + phiOffset);
    yCurved = dist * std::sin(phi + phiOffset);
  }

  /// Transformation from Geant detector centered local coordinates (cm) to
  /// Pixel cell numbers iRow and iCol.
  /// Returns true if point x,z is inside sensitive volume, false otherwise.
  /// A value of -1 for iRow or iCol indicates that this point is outside of the
  /// detector segmentation as defined.
  /// \param float x Detector local coordinate x in cm with respect to
  /// the center of the sensitive volume.
  /// \param float z Detector local coordinate z in cm with respect to
  /// the center of the sensitive volume.
  /// \param int iRow Detector x cell coordinate.
  /// \param int iCol Detector z cell coordinate.
  template <typename T = double, typename TI = int>
  bool localToDetector(T const xRow, T const zCol, TI& iRow, TI& iCol) const noexcept
  {
    localToDetectorUnchecked(xRow, zCol, iRow, iCol);
    if (!isValid<TI>(iRow, iCol)) {
      iRow = iCol = -1;
      return false;
    }
    return true;
  }

  // Same as localToDetector w.o. checks.
  template <typename T = double, typename TI = int>
  void localToDetectorUnchecked(T const xRow, T const zCol, TI& iRow, TI& iCol) const noexcept
  {
    namespace cp = constants::pixelarray;
    T x = cp::length / 2. - xRow; // transformation to upper edge of pixelarray
    T z = zCol + cp::width / 2.;  // transformation to left edge of pixelarray
    iRow = std::floor(x / mPitchRow);
    iCol = std::floor(z / mPitchCol);
  }

  /// Transformation from Detector cell coordinates to Geant detector centered
  /// local coordinates (cm)
  /// \param int iRow Detector x cell coordinate.
  /// \param int iCol Detector z cell coordinate.
  /// \param float x Detector local coordinate x in cm with respect to the
  /// center of the sensitive volume.
  /// \param float z Detector local coordinate z in cm with respect to the
  /// center of the sensitive volume.
  /// If iRow and or iCol is outside of the segmentation range a value of -0.5*Dx()
  /// or -0.5*Dz() is returned.
  template <typename TI = int, typename T = double>
  bool detectorToLocal(TI const iRow, TI const iCol, T& xRow, T& zCol) const noexcept
  {
    detectorToLocalUnchecked(iRow, iCol, xRow, zCol);
    return isValid<T>(xRow, zCol);
  }

  // Same as detectorToLocal w.o. checks.
  // We position ourself in the middle of the pixel.
  template <typename TI = int, typename T = double>
  void detectorToLocalUnchecked(TI const iRow, TI const iCol, T& xRow, T& zCol) const noexcept
  {
    namespace cp = constants::pixelarray;
    xRow = -(iRow + 0.5) * mPitchRow + cp::length / 2.;
    zCol = (iCol + 0.5) * mPitchCol - cp::width / 2.;
  }

  template <typename T = double, typename TI = int>
  bool detectorToLocal(T const row, T const col, T& xRow, T& zCol) const noexcept
  {
    return detectorToLocal(static_cast<TI>(row), static_cast<TI>(col), xRow, zCol);
  }

  template <typename T = double, typename TI = int>
  void detectorToLocalUnchecked(T const row, T const col, T& xRow, T& zCol) const noexcept
  {
    detectorToLocalUnchecked(static_cast<TI>(row), static_cast<TI>(col), xRow, zCol);
  }

  template <typename T = double, typename TI = int>
  bool detectorToLocal(T const row, T const col, math_utils::Point3D<T>& loc) const noexcept
  {
    T xRow, zCol;
    if (detectorToLocal<TI>(row, col, xRow, zCol)) {
      return false;
    }
    loc.SetCoordinates(xRow, 0., zCol);
    return true;
  }

  template <typename T = double, typename TP = T>
  void detectorToLocalUnchecked(T const row, T const col, math_utils::Point3D<TP>& loc) const noexcept
  {
    T xRow, zCol;
    detectorToLocalUnchecked<T, TP>(row, col, xRow, zCol);
    loc.SetCoordinates(xRow, 0., zCol);
  }

  double getEffRadius() const noexcept
  {
    return mEffRadius;
  }

 private:
  template <typename T>
  bool isValid(T const row, T const col) const noexcept
  {
    namespace cp = constants::pixelarray;
    if constexpr (std::is_floating_point_v<T>) {
      return !static_cast<bool>(row < 0. || row >= cp::length || col < 0. || col >= cp::width);
    } else {
      return !static_cast<bool>(row < 0 || row >= mNRows || col < 0 || col >= mNCols);
    }
  }

  const int mLayer; ///< chip layer
  const double mEffRadius{constants::radii[mLayer] + constants::thickness / 2.};
};
} // namespace its3
} // namespace o2

#endif
