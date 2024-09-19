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

#include "ITS3Simulation/DescriptorInnerBarrelITS3.h"
#include "ITS3Base/SpecsV2.h"
#include "fairlogger/Logger.h"
#include "ITSBase/GeometryTGeo.h"
#include "TGeoManager.h"

using namespace o2::its3;

ClassImp(DescriptorInnerBarrelITS3);

void DescriptorInnerBarrelITS3::createLayer(int iLayer, TGeoVolume* dest)
{
  LOGP(info, "ITS3-IB: Creating Layer {}", iLayer);
  mIBLayers[iLayer] = std::make_unique<ITS3Layer>(iLayer);
  mIBLayers[iLayer]->createLayer(dest);
}

void DescriptorInnerBarrelITS3::createServices(TGeoVolume* dest)
{
  LOGP(info, "ITS3-IB: Creating Services");
  mServices = std::make_unique<ITS3Services>();
  mServices->createCYSSAssembly(dest);
}

void DescriptorInnerBarrelITS3::addAlignableVolumesLayer(int idLayer, int wrapperLayerId, TString& parentPath, int& lastUID) const
{
  const TString wrpV = Form("%s%d_1", o2::its::GeometryTGeo::getITSWrapVolPattern(), wrapperLayerId);
  TString path = Form("%s/%s/%s%d_0", parentPath.Data(), wrpV.Data(), o2::its::GeometryTGeo::getITS3LayerPattern(), idLayer);

  for (int iHalfBarrel{0}; iHalfBarrel < 2; ++iHalfBarrel) {
    addAlignableVolumesCarbonForm(idLayer, iHalfBarrel, path, lastUID);
  }
}

void DescriptorInnerBarrelITS3::addAlignableVolumesCarbonForm(int idLayer, int iHalfBarrel, TString& parentPath, int& lastUID) const
{
  TString path = parentPath;
  path = Form("%s/%s_%d", parentPath.Data(), o2::its::GeometryTGeo::getITS3CarbonFormPattern(idLayer), iHalfBarrel);
  const TString sname = o2::its::GeometryTGeo::composeSymNameHalfBarrel(idLayer, iHalfBarrel, true);
  if (!gGeoManager->SetAlignableEntry(sname.Data(), path.Data())) {
    LOG(fatal) << "Unable to set alignable entry ! " << sname << " : " << path;
  }

  addAlignableVolumesChip(idLayer, iHalfBarrel, parentPath, lastUID);
}

void DescriptorInnerBarrelITS3::addAlignableVolumesChip(int idLayer, int iHalfBarrel, TString& parentPath, int& lastUID) const
{
  using gITS = o2::its::GeometryTGeo;
  for (unsigned int iSegment{0}; iSegment < constants::nSegments[idLayer]; ++iSegment) {
    for (unsigned int iRSU{0}; iRSU < constants::segment::nRSUs; ++iRSU) {
      for (unsigned int iTile{0}; iTile < constants::rsu::nTiles; ++iTile) {
        TString path = parentPath;
        path += Form("/%s_%d/", gITS::getITS3CarbonFormPattern(idLayer), iHalfBarrel);
        path += Form("%s_0/", gITS::getITS3ChipPattern(idLayer));
        path += Form("%s_%d/", gITS::getITS3SegmentPattern(idLayer), iSegment);
        path += Form("%s_%d/", gITS::getITS3RSUPattern(idLayer), iRSU);
        path += Form("%s_%d/", gITS::getITS3TilePattern(idLayer), iTile);
        path += Form("%s_0", gITS::getITS3PixelArrayPattern(idLayer));
        if (!gGeoManager->CheckPath(path.Data())) {
          LOG(fatal) << "Path does not exist: " << path;
        }

        TString sname = o2::its::GeometryTGeo::composeSymNameChip(idLayer, iHalfBarrel, iSegment, iRSU, iTile, 0, true);
        int modUID = o2::base::GeometryManager::getSensID(o2::detectors::DetID::IT3, lastUID++);

        if (!gGeoManager->SetAlignableEntry(sname, path.Data(), modUID)) {
          LOG(fatal) << "Unable to set alignable entry ! " << sname << " : " << path;
        }
      }
    }
  }
}
