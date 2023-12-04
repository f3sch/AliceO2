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

/// \file Detector.h
/// \brief Definition of the Detector class
/// \author felix.schlepper@cern.ch

#ifndef ALICEO2_ITS3_DETECTOR_H_
#define ALICEO2_ITS3_DETECTOR_H_

#include "Rtypes.h"      // for Int_t, Double_t, Float_t, Bool_t, etc
#include "TGeoManager.h" // for gGeoManager, TGeoManager (ptr only)

#include "ITSMFTSimulation/Hit.h" // for Hit
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Detector.h"

namespace o2::its3
{

class Detector : public o2::base::DetImpl<Detector>
{
 public:
  Detector();
  ~Detector() override;

  Bool_t ProcessHits(FairVolume* v = nullptr) override;
  void Register() override;

 public:
  void Reset() override;
  void ConstructGeometry() override;
  void InitGeometry() override {}
  void GeneratePrimaries() override {}
  void BeginEvent() override {}
  void BeginPrimary() override {}
  void PreTrack() override {}
  void Stepping() override {}
  void PostTrack() override {}
  void FinishPrimary() override {}
  void FinishEvent() override {}

 private:
  void createMaterials();
  void createDetectorGeometry();

 private:
  std::vector<o2::itsmft::Hit>* mHits{o2::utils::createSimVector<o2::itsmft::Hit>()}; // Container for hit data
};

} // namespace o2::its3

#endif
