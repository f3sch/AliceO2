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

/// \file GeometryTGeo.h
/// \brief Definition of the GeometryTGeo class for ITS3
/// \author felix.schlepper@cern.ch

#ifndef ALICEO2_ITS3_GEOMETRYTGEO_H_
#define ALICEO2_ITS3_GEOMETRYTGEO_H_

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "DetectorsBase/GeometryManager.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsCommonDataFormats/DetMatrixCache.h"
#include "MathUtils/Utils.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITS3Base/Specs.h"

namespace o2::its3
{
class GeometryTGeo final : public o2::detectors::DetMatrixCache
{
  using ITS2TGeo = o2::its::GeometryTGeo;

 public:
  using Mat3D = o2::math_utils::Transform3D;
  using DetMatrixCache::getMatrixL2G;
  using DetMatrixCache::getMatrixT2G;
  using DetMatrixCache::getMatrixT2GRot;
  using DetMatrixCache::getMatrixT2L;
  using o2::detectors::DetMatrixCache::fillMatrixCache;

  GeometryTGeo() = default;
  GeometryTGeo(const o2::detectors::DetID& detid) : o2::detectors::DetMatrixCache(detid) {}

  GeometryTGeo(const GeometryTGeo& src) = delete;
  GeometryTGeo(GeometryTGeo&& src) = delete;
  GeometryTGeo& operator=(const GeometryTGeo& geom) = delete;
  GeometryTGeo& operator=(GeometryTGeo&& geom) = delete;
  ~GeometryTGeo() final = default;

  static GeometryTGeo* Instance()
  {
    // get (create if needed) a unique instance of the object
    if (!mInstance) {
      mInstance = std::make_unique<GeometryTGeo>();
    }
    return mInstance.get();
  }

  // Extracts ITS3 Parameters from TGeo.
  void Build(int mask);

  void fillMatrixCache(int mask) override;

  // Get Layer number from Index
  [[nodiscard]] int getLayer(int index) const noexcept;
  // Get Chip ID from Index in ITS3
  [[nodiscard]] bool getChipIdITS3(int index, int lay, int& chip, int& segment, int& rsu, int& tile) const;
  // Get Last Chip ID of layer
  [[nodiscard]] int getLastChipIndex(int lay) const { return mLastChipIndex[lay]; }
  // Get First Chip ID of layer
  [[nodiscard]] int getFirstChipIndex(int lay) const { return (lay == 0) ? 0 : mLastChipIndex[lay - 1] + 1; }

  // TODO FS better to query this from the TGeo
  static int getNumberOfLayers() { return 7; }
  static int getNumberOfChipsPerLayer(int layer) { return 2; }

  // TODO FS placeholders
  float getAlphaFromGlobalITS3(int isn, o2::math_utils::Point3D<float> gloXYZ) { return 0; }
  float getSensorRefAlpha(int isn) const { return mCacheRefAlpha[isn]; }
  const Mat3D getT2LMatrixITS3(int isn, float alpha) { return {}; }

  // create matrix for transformation from sensor local frame to global one
  TGeoHMatrix& createT2LMatrix(int isn);

 private:
  static std::unique_ptr<GeometryTGeo> mInstance; ///< singletone instance
  ITS2TGeo* mITS2Instance;

  static constexpr std::array<bool, constants::nTotLayers> mIsITS3Layer{true, true, true, false, false, false, false}; ///< mask indicating a new layer
  std::array<int, constants::nITS2Layers> mITS2NumberOfStaves{};                                                       ///< ITS2 number of staves/layer(layer)
  std::array<int, constants::nITS2Layers> mITS2NumberOfHalfStaves{};                                                   ///< ITS2 the number of substaves/stave(layer)
  std::array<int, constants::nITS2Layers> mITS2NumberOfModules{};                                                      ///< ITS2 number of modules/substave(layer)
  std::array<int, constants::nITS2Layers> mITS2NumberOfChipsPerModule{};                                               ///< ITS2 number of chips per module (group of chips on substaves)
  std::array<int, constants::nITS2Layers> mITS2NumberOfChipRowsPerModule{};                                            ///< ITS2 number of chips rows per module (relevant for OB modules)
  std::array<int, constants::nITS2Layers> mITS2NumberOfChipsPerHalfStave{};                                            ///< ITS2 number of chips per substave
  std::array<int, constants::nITS2Layers> mITS2NumberOfChipsPerStave{};                                                ///< ITS2 number of chips per stave
  std::array<int, constants::nITS2Layers> mITS2NumberOfChipsPerHalfBarrel{};                                           ///< ITS2 number of chips per halfbarrel
  std::array<int, constants::nITS2Layers> mITS2NumberOfChipsPerLayer{};                                                ///< ITS2 number of chips per layer
  std::array<int, constants::nITS3Layers> mITS3NumberOfTilesPerLayer{};                                                ///< ITS3 number of chips/tiles per layer
  std::array<int, constants::nITS3Layers> mITS3NumberOfTilesPerChip{};                                                 ///< ITS3 number of chips/tiles per chip
  std::array<int, constants::nTotLayers> mLastChipIndex{};                                                             ///< max ID of the detctor in the layer

  // Cache tracking frame for every sensor
  void fillTrackingFramesCache();
  [[nodiscard]] bool isTrackingFrameCached() const { return !mCacheRefX.empty(); }
  std::vector<double> mCacheRefX{}; ///< sensors tracking plane reference X
  void extractSensorXAlpha(int isn, double& x, double& alp);
  std::vector<double> mCacheRefAlpha{}; ///< sensors tracking plane reference alpha

  // Extract Matrix of Sensor
  [[nodiscard]] TGeoHMatrix* extractMatrixSensor(int index) const;

  // Name Patterns for the IB
  static std::string sLayerNameITS3;      ///< ITS3 Layer name
  static std::string sCarbonFormNameITS3; ///< ITS3 CarbonForm name
  static std::string sChipNameITS3;       ///< ITS3 Chip name
  static std::string sSegmentNameITS3;    ///< ITS3 Segment name
  static std::string sRSUNameITS3;        ///< ITS3 RSU name
  static std::string sTileNameITS3;       ///< ITS3 Tile name
  static std::string sPixelArrayNameITS3; ///< ITS3 PixelArray name

 public:
  static const char* getITS3LayerPatternRaw() { return sLayerNameITS3.c_str(); };
  static const char* getITS3LayerPattern(int layer = 0) { return Form("%s_%d", getITS3LayerPatternRaw(), layer); };
  static const char* getITS3CarbonFormPatternRaw() { return sCarbonFormNameITS3.c_str(); };
  static const char* getITS3CarbonFormPattern(int layer = 0) { return Form("%s_%d", getITS3CarbonFormPatternRaw(), layer); };
  static const char* getITS3ChipPatternRaw() { return sChipNameITS3.c_str(); };
  static const char* getITS3ChipPattern(int layer = 0) { return Form("%s_%d", getITS3ChipPatternRaw(), layer); };
  static const char* getITS3SegmentPatternRaw() { return sSegmentNameITS3.c_str(); };
  static const char* getITS3SegmentPattern(int layer = 0) { return Form("%s_%d", getITS3SegmentPatternRaw(), layer); };
  static const char* getITS3RSUPatternRaw() { return sRSUNameITS3.c_str(); };
  static const char* getITS3RSUPattern(int layer = 0) { return Form("%s_%d", getITS3RSUPatternRaw(), layer); };
  static const char* getITS3TilePatternRaw() { return sTileNameITS3.c_str(); };
  static const char* getITS3TilePattern(int layer = 0) { return Form("%s_%d", getITS3TilePatternRaw(), layer); };
  static const char* getITS3PixelArrayPatternRaw() { return sPixelArrayNameITS3.c_str(); };
  static const char* getITS3PixelArrayPattern(int layer = 0) { return Form("%s_%d", getITS3PixelArrayPatternRaw(), layer); };

 private:
  ClassDefOverride(GeometryTGeo, 0);
};
} // namespace o2::its3

#endif
