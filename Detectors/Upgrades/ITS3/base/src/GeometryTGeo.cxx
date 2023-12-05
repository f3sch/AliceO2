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
// #include "ITSBase/GeometryTGeo.h"

// using ITS2TGeo = o2::its::GeometryTGeo;

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

  // Count the number of tiles/chips
  int numberOfChips{0};
  for (int iLayer{0}; iLayer < constants::nTotLayers; ++iLayer) {
    if (mIsITS3Layer[iLayer]) {
      mITS3NumberOfTilesPerChip[iLayer] = constants::rsu::nTiles * constants::nSegments[iLayer] * constants::segment::nRSUs;
      numberOfChips += mITS3NumberOfTilesPerLayer[iLayer] = 2 * mITS3NumberOfTilesPerChip[iLayer];
    } else {
      // TODO
    }
    mLastChipIndex[iLayer] = numberOfChips - 1;
  }
  setSize(numberOfChips);

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

TGeoHMatrix* GeometryTGeo::extractMatrixSensor(unsigned int index) const
{

  // extract matrix transforming from the PHYSICAL sensor frame to global one
  // Note, the if the effective sensitive layer thickness is smaller than the
  // total physical sensor tickness, this matrix is biased and connot be used
  // directly for transformation from sensor frame to global one. Therefore, we
  // need to add a shift.

  unsigned int lay = getLayer(index), hba, stav, sstav, mod, chipInMod, hemi, segment, rsu, tile;
  // TString path = Form("/cave_1/barrel_1/%s_2/", GeometryTGeo::getITSVolPattern());
  if (mIsITS3Layer[lay]) {
    if (!getChipIdITS3(index, lay, hemi, segment, rsu, tile)) {
      return nullptr;
    }
  } else {
    if (!getChipIdITS2(index, lay, hba, stav, sstav, mod, chipInMod)) {
      return nullptr;
    }
  }
  return nullptr;
}

void GeometryTGeo::fillMatrixCache(int mask)
{
  // populate matrix cache for requested transformations
}

bool GeometryTGeo::getChipIdITS2(unsigned int index, unsigned int lay, unsigned int& hba, unsigned int& sta, unsigned int& hsta, unsigned int& mod, unsigned int& chip) const
{
  if (mIsITS3Layer[lay]) {
    return false;
  }
  index -= getFirstChipIndex(lay);
  auto its2layer = lay - constants::nITS3Layers;
  hba = index / mITS2NumberOfChipsPerHalfBarrel[its2layer];
  index %= mITS2NumberOfChipsPerHalfBarrel[its2layer];
  sta = index / mITS2NumberOfChipsPerStave[its2layer];
  index %= mITS2NumberOfChipsPerStave[its2layer];
  hsta = mITS2NumberOfHalfStaves[its2layer] > 0 ? index / mITS2NumberOfChipsPerHalfStave[its2layer] : -1;
  index %= mITS2NumberOfChipsPerHalfStave[its2layer];
  mod = mITS2NumberOfModules[its2layer] > 0 ? index / mITS2NumberOfChipsPerModule[its2layer] : -1;
  chip = index % mITS2NumberOfChipsPerModule[its2layer];

  return true;
}

bool GeometryTGeo::getChipIdITS3(unsigned int index, unsigned int lay, unsigned int& chip, unsigned int& segment, unsigned int& rsu, unsigned int& tile) const
{
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

unsigned int GeometryTGeo::getLayer(unsigned int index) const noexcept
{
  unsigned int lay = 0;
  while (index > mLastChipIndex[lay]) {
    lay++;
  }
  return lay;
}

} // namespace o2::its3
