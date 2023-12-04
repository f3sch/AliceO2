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

/// \file DescriptorInnerBarrelITS3.h
/// \brief Definition of the DescriptorInnerBarrelITS3 class
/// \author felix.schlepper@cern.ch

#ifndef ALICEO2_ITS3_DESCRIPTORINNERBARRELITS3_H
#define ALICEO2_ITS3_DESCRIPTORINNERBARRELITS3_H

#include <string>
#include <vector>

#include "ITSBase/DescriptorInnerBarrel.h"
#include "ITS3Simulation/ITS3Layer.h"

namespace o2
{
namespace its3
{

class DescriptorInnerBarrelITS3 : public o2::its::DescriptorInnerBarrel
{
 public:
  ITS3Layer* createLayer(int idLayer, TGeoVolume* dest);
  void createServices(TGeoVolume* dest);
  void configure() {}

 private:
  ClassDefNV(DescriptorInnerBarrelITS3, 0); /// ITS3 inner barrel geometry descriptor
};
} // namespace its3
} // namespace o2

#endif
