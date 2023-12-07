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

#include "Rtypes.h"
#include "TGeoManager.h"
#include "TLorentzVector.h"

#include "ITSMFTSimulation/Hit.h"
#include "DetectorsBase/Detector.h"
#include "ITS3Base/GeometryTGeo.h"
#include "ITS3Base/Specs.h"

namespace o2::its3
{

class Detector : public o2::base::DetImpl<Detector>
{
 public:
  Detector(bool isActive = true);
  ~Detector() override;
  Detector(const Detector&) = default;
  Detector& operator=(const Detector&) = default;

  Bool_t ProcessHits(FairVolume* vol = nullptr) override;
  void Register() override;

 public:
  void InitializeO2Detector() override;
  void Reset() override
  {
    if (!o2::utils::ShmManager::Instance().isOperational()) {
      mHits->clear();
    }
  }
  void ConstructGeometry() override;
  void BeginPrimary() override {}
  void PreTrack() override {}
  void PostTrack() override {}
  void FinishPrimary() override {}

  /// Gets the produced collections
  std::vector<o2::itsmft::Hit>* mHits{o2::utils::createSimVector<o2::itsmft::Hit>()}; // Container for hit data
  std::vector<o2::itsmft::Hit>* getHits(Int_t iColl) const
  {
    if (iColl == 0) {
      return mHits;
    }
    return nullptr;
  }

 private:
  /// this is transient data about track passing the sensor
  struct TrackData {               // this is transient
    bool mHitStarted;              //! hit creation started
    unsigned char mTrkStatusStart; //! track status flag
    TLorentzVector mPositionStart; //! position at entrance
    TLorentzVector mMomentumStart; //! momentum
    double mEnergyLoss;            //! energy loss
  } mTrackData;                    //!

  /// Creates all the materials
  void createMaterials();

  /// Creates the Detector Geometry
  void createDetectorGeometry();

  /// Define the sensitive volumes of the Geometry
  void defineSensitiveVolumes();

  o2::itsmft::Hit* addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos,
                          const TVector3& startMom, double startE, double endTime, double eLoss,
                          unsigned char startStatus, unsigned char endStatus);

  GeometryTGeo* mTGeo{nullptr};                            //! GeometryTGeo instance
  std::array<int, constants::nTotLayers> mLayerID{};       //! layer identifiers, e.g. the pixelarray/sensor in each layer
  std::array<TString, constants::nTotLayers> mLayerName{}; //! layer names, e.g. the pixelarray/sensor in each layer

  template <typename Det>
  friend class o2::base::DetImpl;
  ClassDefOverride(Detector, 0);
};
} // namespace o2::its3

#ifdef USESHM
namespace o2
{
namespace base
{
template <>
struct UseShm<o2::its3::Detector> {
  static constexpr bool value = true;
};
} // namespace base
} // namespace o2
#endif

#endif
