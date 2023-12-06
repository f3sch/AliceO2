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

/// \file GeometryTGeo.cxx
/// \brief Implementation of the GeometryTGeo class
/// \author felix.schlepper@cern.ch

#include "ITS3Base/GeometryTGeo.h"
#include "ITS3Base/SegmentationSuperAlpide.h"
#include "ITS3Base/Specs.h"
#include "CommonConstants/MathConstants.h"

namespace o2m = o2::constants::math;

ClassImp(o2::its3::GeometryTGeo);

namespace o2::its3
{
std::unique_ptr<GeometryTGeo> GeometryTGeo::mInstance;

std::string GeometryTGeo::sLayerNameITS3 = "ITS3Layer";
std::string GeometryTGeo::sCarbonFormNameITS3 = "ITS3CarbonForm";
std::string GeometryTGeo::sChipNameITS3 = "ITS3Chip";
std::string GeometryTGeo::sSegmentNameITS3 = "ITS3Segment";
std::string GeometryTGeo::sRSUNameITS3 = "ITS3RSU";
std::string GeometryTGeo::sTileNameITS3 = "ITS3Tile";
std::string GeometryTGeo::sPixelArrayNameITS3 = "ITS3PixelArray";

void GeometryTGeo::Build(int mask)
{
  if (isBuilt()) {
    LOGP(warn, "ITS3 alreay built, skipping...");
    return;
  }

  if (gGeoManager == nullptr) {
    LOGP(fatal, "gGeoManager is not loaded!");
  }

  // we require a pointer to the ITS2 TGeo instance but do not build it
  mITS2Instance = ITS2TGeo::Instance(false);
  if (mITS2Instance == nullptr) {
    LOGP(fatal, "ITS2::Instance ref is null");
  }

  // Count the number of tiles/chips
  int numberOfChips{0};
  for (int iLayer{0}; iLayer < constants::nTotLayers; ++iLayer) {
    if (mIsITS3Layer[iLayer]) {
      mITS3NumberOfTilesPerChip[iLayer] = constants::rsu::nTiles * constants::nSegments[iLayer] * constants::segment::nRSUs;
      numberOfChips += mITS3NumberOfTilesPerLayer[iLayer] = 2 * mITS3NumberOfTilesPerChip[iLayer];
    } else {
      // TODO FS
    }
    mLastChipIndex[iLayer] = numberOfChips - 1;
  }
  setSize(numberOfChips);
  fillTrackingFramesCache();
  fillMatrixCache(mask);
}

void GeometryTGeo::fillTrackingFramesCache()
{
  // fill for every sensor its tracking frame parameteres
  if (!isTrackingFrameCached()) {
    // special cache for sensors tracking frame X and alpha params
    mCacheRefX.resize(mSize);
    mCacheRefAlpha.resize(mSize);
    for (int i = 0; i < mSize; i++) {
      extractSensorXAlpha(i, mCacheRefX[i], mCacheRefAlpha[i]);
    }
  }
}

void GeometryTGeo::extractSensorXAlpha(int isn, double& x, double& alp)
{
  // calculate r and phi of the impact of the normal on the sensor
  // (i.e. phi of the tracking frame alpha and X of the sensor in this frame)
  const auto matL2G = extractMatrixSensor(isn);
  double locA[3] = {-100., 0., 0.}, locB[3] = {100., 0., 0.}, gloA[3], gloB[3];
  int iLayer = getLayer(isn);

  if (mIsITS3Layer[iLayer]) {
    // in this case we need the line tangent to the circumference
    double radius = 0.;
    SegmentationSuperAlpide seg(iLayer);
    radius = seg.getEffRadius();
    locA[1] = radius;
    locB[1] = radius;
  }

  matL2G->LocalToMaster(locA, gloA);
  matL2G->LocalToMaster(locB, gloB);
  double dx = gloB[0] - gloA[0], dy = gloB[1] - gloA[1];
  double t = (gloB[0] * dx + gloB[1] * dy) / (dx * dx + dy * dy);
  double xp = gloB[0] - dx * t, yp = gloB[1] - dy * t;
  x = std::hypot(xp, yp);
  alp = std::atan2(yp, xp);
  o2::math_utils::bringTo02Pid(alp);
}

TGeoHMatrix* GeometryTGeo::extractMatrixSensor(int index) const
{
  int lay = getLayer(index);
  if (!mIsITS3Layer[lay]) {
    if (const auto mat = mITS2Instance->extractMatrixSensor(index); mat == nullptr) {
      LOGP(fatal, "Failed to extract ITS2 Matrix from index {}", index);
    } else {
      return mat;
    }
  }

  int chip, segment, rsu, tile;
  TString path = Form("/cave_1/barrel_1/%s_2/", mITS2Instance->getITSVolPattern());
  if (getChipIdITS3(index, lay, chip, segment, rsu, tile)) {
    path += Form("%s_0/%s_%d/%s_0/%s_%d/%s_%d/%s_%d/%s_0", getITS3LayerPattern(lay), getITS3CarbonFormPattern(lay), chip, getITS3ChipPattern(lay),
                 getITS3SegmentPattern(lay), segment, getITS3RSUPattern(lay), rsu, getITS3TilePattern(lay), tile, getITS3PixelArrayPattern(lay));
  } else {
    LOGP(fatal, "Failed to extract ITS3 Matrix from Geometry with index {} -> layer={} chip={} segment={} rsu={} tile={}", index, lay, chip, segment, rsu, tile);
  }
  if (!gGeoManager->CheckPath(path.Data())) {
    LOGP(error, "Pathquery failed: {}", path.Data());
    return nullptr;
  }

  static TGeoHMatrix matTmp;
  gGeoManager->PushPath();
  if (!gGeoManager->cd(path.Data())) {
    gGeoManager->PopPath();
    LOGP(error, "Error in cd-ing to {}", path.Data());
    return nullptr;
  }

  matTmp = *gGeoManager->GetCurrentMatrix(); // matrix will change after cd;
  gGeoManager->PopPath();

  return &matTmp;
}

void GeometryTGeo::fillMatrixCache(int mask)
{
  if (mSize < 1) {
    LOG(warning) << "The method Build was not called yet";
    Build(mask);
    return;
  }

  // build matrices
  if ((mask & o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G)) && !getCacheL2G().isFilled()) {
    // Matrices for Local (Sensor!!! rather than the full chip) to Global frame transformation
    LOG(info) << "Loading ITS2 and ITS3 L2G matrices from TGeo";
    auto& cacheL2G = getCacheL2G();
    cacheL2G.setSize(mSize);

    for (int i = 0; i < mSize; i++) {
      TGeoHMatrix* hm = extractMatrixSensor(i);
      cacheL2G.setMatrix(Mat3D(*hm), i);
    }
  }

  if ((mask & o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L)) && !getCacheT2L().isFilled()) {
    // matrices for Tracking to Local (Sensor!!! rather than the full chip) frame transformation
    LOG(info) << "Loading ITS2 and ITS3 T2L matrices from TGeo";
    auto& cacheT2L = getCacheT2L();
    cacheT2L.setSize(mSize);
    for (int i = 0; i < mSize; i++) {
      TGeoHMatrix& hm = createT2LMatrix(i);
      cacheT2L.setMatrix(Mat3D(hm), i);
    }
  }

  if ((mask & o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2G)) && !getCacheT2G().isFilled()) {
    LOG(warning) << "It is faster to use 2D rotation for T2G instead of full Transform3D matrices";
    // matrices for Tracking to Global frame transformation
    LOG(info) << "Loading ITS2 and ITS3 T2G matrices from TGeo";
    auto& cacheT2G = getCacheT2G();
    cacheT2G.setSize(mSize);
    for (int i = 0; i < mSize; i++) {
      TGeoHMatrix& mat = createT2LMatrix(i);
      mat.MultiplyLeft(extractMatrixSensor(i));
      cacheT2G.setMatrix(Mat3D(mat), i);
    }
  }

  if ((mask & o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2GRot)) && !getCacheT2GRot().isFilled()) {
    // 2D rotation matrices for Tracking frame to Global rotations
    LOG(info) << "Loading ITS2 and ITS3 T2G rotation 2D matrices";
    auto& cacheT2Gr = getCacheT2GRot();
    cacheT2Gr.setSize(mSize);
    for (int i = 0; i < mSize; i++) {
      cacheT2Gr.setMatrix(Rot2D(getSensorRefAlpha(i)), i);
    }
  }
}

bool GeometryTGeo::getChipIdITS3(int index, int lay, int& chip, int& segment, int& rsu, int& tile) const
{
  // TODO FS has to reflect filling in fillMatrixCache
  if (!mIsITS3Layer[lay]) {
    return false;
  }
  index -= getFirstChipIndex(lay);
  chip = index / mITS3NumberOfTilesPerChip[lay];
  index %= mITS3NumberOfTilesPerChip[lay];
  segment = index / constants::segment::nTilesPerSegment;
  index -= segment * constants::segment::nTilesPerSegment;
  rsu = index / constants::segment::nRSUs;
  tile = index % constants::rsu::nTiles;
  return true;
}

int GeometryTGeo::getLayer(int index) const noexcept
{
  int lay = 0;
  while (index > mLastChipIndex[lay]) {
    lay++;
  }
  return lay;
}

TGeoHMatrix& GeometryTGeo::createT2LMatrix(int isn)
{
  // create for sensor isn the TGeo matrix for Tracking to Local frame transformations
  static TGeoHMatrix t2l;
  double x = 0.f, alp = 0.f;
  extractSensorXAlpha(isn, x, alp);
  t2l.Clear();
  t2l.RotateZ(alp * o2m::Rad2Deg); // rotate in direction of normal to the sensor plane
  const TGeoHMatrix* matL2G = extractMatrixSensor(isn);
  const TGeoHMatrix& matL2Gi = matL2G->Inverse();
  t2l.MultiplyLeft(&matL2Gi);
  return t2l;
}

} // namespace o2::its3
