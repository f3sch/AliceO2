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
///
/// \file IndexTableUtils.h
/// \brief
///

#ifndef TRACKINGITSU_INCLUDE_INDEXTABLEUTILS_H_
#define TRACKINGITSU_INCLUDE_INDEXTABLEUTILS_H_

#include <array>

#include "ITStracking/Constants.h"
#include "ITStracking/Configuration.h"
#include "ITStracking/Definitions.h"
#include "CommonConstants/MathConstants.h"
#include "GPUCommonMath.h"
#include "GPUCommonDef.h"

namespace o2::its
{

template <int nLayers>
class IndexTableUtils
{
 public:
  template <class T>
  void setTrackingParameters(const T& params);
  float getInverseZCoordinate(const int layerIndex) const noexcept;
  GPUhdi() int getZBinIndex(const int, const float) const noexcept;
  GPUhdi() int getPhiBinIndex(const float) const noexcept;
  GPUhdi() int getBinIndex(const int, const int) const noexcept;
  GPUhdi() int countRowSelectedBins(const int*, const int, const int, const int) const noexcept;
  GPUhdi() void print() const;

  GPUhdi() int getNzBins() const noexcept { return mNzBins; }
  GPUhdi() int getNphiBins() const noexcept { return mNphiBins; }
  GPUhdi() float getLayerZ(int i) const noexcept { return mLayerZ[i]; }
  GPUhdi() void setNzBins(const int zBins) noexcept { mNzBins = zBins; }
  GPUhdi() void setNphiBins(const int phiBins) noexcept { mNphiBins = phiBins; }

 private:
  int mNzBins = 0;
  int mNphiBins = 0;
  float mInversePhiBinSize = 0.f;
  std::array<float, nLayers> mLayerZ{};
  std::array<float, nLayers> mInverseZBinSize{};
};

template <int nLayers>
template <class T>
inline void IndexTableUtils<nLayers>::setTrackingParameters(const T& params)
{
  mInversePhiBinSize = params.PhiBins / o2::constants::math::TwoPI;
  mNzBins = params.ZBins;
  mNphiBins = params.PhiBins;
  for (int iLayer{0}; iLayer < params.LayerZ.size(); ++iLayer) {
    mLayerZ[iLayer] = params.LayerZ[iLayer];
  }
  for (unsigned int iLayer{0}; iLayer < params.LayerZ.size(); ++iLayer) {
    mInverseZBinSize[iLayer] = 0.5f * params.ZBins / params.LayerZ[iLayer];
  }
}

template <int nLayers>
inline float IndexTableUtils<nLayers>::getInverseZCoordinate(const int layerIndex) const noexcept
{
  return 0.5f * mNzBins / mLayerZ[layerIndex];
}

template <int nLayers>
GPUhdi() int IndexTableUtils<nLayers>::getZBinIndex(const int layerIndex, const float zCoordinate) const noexcept
{
  return (zCoordinate + mLayerZ[layerIndex]) * mInverseZBinSize[layerIndex];
}

template <int nLayers>
GPUhdi() int IndexTableUtils<nLayers>::getPhiBinIndex(const float currentPhi) const noexcept
{
  return (currentPhi * mInversePhiBinSize);
}

template <int nLayers>
GPUhdi() int IndexTableUtils<nLayers>::getBinIndex(const int zIndex, const int phiIndex) const noexcept
{
  return o2::gpu::GPUCommonMath::Min((phiIndex * mNzBins) + zIndex, (mNzBins * mNphiBins) - 1);
}

template <int nLayers>
GPUhdi() int IndexTableUtils<nLayers>::countRowSelectedBins(const int* indexTable, const int phiBinIndex, const int minZBinIndex, const int maxZBinIndex) const noexcept
{
  const int firstBinIndex{getBinIndex(minZBinIndex, phiBinIndex)};
  const int maxBinIndex{firstBinIndex + maxZBinIndex - minZBinIndex + 1};
  return indexTable[maxBinIndex] - indexTable[firstBinIndex];
}

template <int nLayers>
GPUhdi() void IndexTableUtils<nLayers>::print() const
{
  printf("NzBins: %d, NphiBins: %d, InversePhiBinSize: %f\n", mNzBins, mNphiBins, mInversePhiBinSize);
  for (int iLayer{0}; iLayer < nLayers; ++iLayer) {
    printf("Layer %d: Z: %f, InverseZBinSize: %f\n", iLayer, mLayerZ[iLayer], mInverseZBinSize[iLayer]);
  }
}

} // namespace o2::its
#endif /* TRACKINGITSU_INCLUDE_INDEXTABLEUTILS_H_ */
