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

/// \author felix.schlepper@cern.ch

#include "ITS3Simulation/Detector.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITSMFTSimulation/Hit.h"
#include "ITS3Base/GeometryTGeo.h"
#include "DetectorsBase/Stack.h"
#include "SimulationDataFormat/TrackReference.h"
#include "TVirtualMC.h"
#include "TVirtualMCStack.h" // for TVirtualMCStack
#include "Rtypes.h"

#include "FairDetector.h"      // for FairDetector
#include <fairlogger/Logger.h> // for LOG, LOG_IF
#include "FairRootManager.h"   // for FairRootManager
#include "FairRun.h"           // for FairRun
#include "FairRuntimeDb.h"     // for FairRuntimeDb
#include "FairVolume.h"        // for FairVolume

#include "fairlogger/Logger.h"

#include <cstddef>

using ITS2g = o2::its::GeometryTGeo;
using ITS3g = o2::its3::GeometryTGeo;
using o2::itsmft::Hit;

ClassImp(o2::its3::Detector);

namespace o2::its3
{
using namespace constants;

Detector::Detector(bool isActive) : o2::base::DetImpl<Detector>("IT3", isActive)
{
  for (int iLayer{0}; iLayer < nTotLayers; ++iLayer) {
    if (isITS3Layer[iLayer]) {
      mLayerName[iLayer] = ITS3g::getITS3PixelArrayPattern(iLayer);
    } else {
      mLayerName[iLayer].Form("%s%d", ITS2g::getITSSensorPattern(), iLayer);
    }
  }
}

Detector::~Detector()
{
  if (mHits != nullptr) {
    o2::utils::freeSimVector(mHits);
  }
}

void Detector::InitializeO2Detector()
{
  defineSensitiveVolumes();

  for (int iLayer{0}; iLayer < nTotLayers; ++iLayer) {
    mLayerID[iLayer] = gMC ? TVirtualMC::GetMC()->VolId(mLayerName[iLayer]) : 0;
  }

  mTGeo = ITS3g::Instance();
}

void Detector::Register()
{
  if (FairRootManager::Instance() != nullptr) {
    FairRootManager::Instance()->RegisterAny(addNameTo("Hit").data(), mHits, kTRUE);
  }
}

Bool_t Detector::ProcessHits(FairVolume* vol)
{
  if (!(fMC->TrackCharge())) {
    return kFALSE;
  }
  Int_t layer{0}, volID{vol->getMCid()};
  bool notSens{false};
  while ((layer < nTotLayers) && (notSens = (volID != mLayerID[layer]))) {
    ++layer;
  }
  if (notSens) {
    LOGP(error, "Sensor/PixelArray in Layer {} is not sensitive!", layer);
    return kFALSE; // This should never happen in principle!
  }

  auto stack{(o2::data::Stack*)fMC->GetStack()};
  if (fMC->IsTrackExiting()) {
    // Keep the track refs for the innermost and outermost layers only
    o2::TrackReference tr(*fMC, GetDetId());
    tr.setTrackID(stack->GetCurrentTrackNumber());
    tr.setUserId(layer);
    stack->addTrackReference(tr);
  }

  bool startHit{false}, stopHit{false};
  unsigned char status{0};
  if (fMC->IsTrackEntering()) {
    status |= Hit::kTrackEntering;
  }
  if (fMC->IsTrackInside()) {
    status |= Hit::kTrackInside;
  }
  if (fMC->IsTrackExiting()) {
    status |= Hit::kTrackExiting;
  }
  if (fMC->IsTrackOut()) {
    status |= Hit::kTrackOut;
  }
  if (fMC->IsTrackStop()) {
    status |= Hit::kTrackStopped;
  }
  if (fMC->IsTrackAlive()) {
    status |= Hit::kTrackAlive;
  }

  // track is entering or created in the volume
  if ((status & Hit::kTrackEntering) || (status & Hit::kTrackInside && !mTrackData.mHitStarted)) {
    startHit = true;
  } else if ((status & (Hit::kTrackExiting | Hit::kTrackOut | Hit::kTrackStopped))) {
    stopHit = true;
  }

  if (!startHit) {
    // increment energy loss at all steps except entrance
    mTrackData.mEnergyLoss += fMC->Edep();
  }

  if (!(startHit | stopHit)) {
    return kFALSE; // do nothing
  }

  if (startHit) {
    mTrackData.mEnergyLoss = 0.;
    fMC->TrackMomentum(mTrackData.mMomentumStart);
    fMC->TrackPosition(mTrackData.mPositionStart);
    mTrackData.mTrkStatusStart = status;
    mTrackData.mHitStarted = true;
  }

  if (stopHit) {
    TLorentzVector positionStop;
    fMC->TrackPosition(positionStop);
    int index{0}, l1{0}, l2{0}, l3{0}, l4{0}, l5{0};
    fMC->CurrentVolOffID(1, l1);
    fMC->CurrentVolOffID(2, l2);
    fMC->CurrentVolOffID(3, l3);
    fMC->CurrentVolOffID(4, l4);
    fMC->CurrentVolOffID(5, l5);
    index = mTGeo->getChipIndex(layer, l5, l4, l3, l2, l1);

    Hit* p = addHit(stack->GetCurrentTrackNumber(), index, mTrackData.mPositionStart.Vect(), positionStop.Vect(),
                    mTrackData.mMomentumStart.Vect(), mTrackData.mMomentumStart.E(), positionStop.T(),
                    mTrackData.mEnergyLoss, mTrackData.mTrkStatusStart, status);
    // Increment number of Detector det points in TParticle
    stack->addHit(GetDetId());
  }

  return kTRUE;
}

void Detector::createMaterials()
{
  Int_t ifield = 2;
  Float_t fieldm = 10.0;
  o2::base::Detector::initFieldTrackingParams(ifield, fieldm);
  Float_t tmaxfd = 0.1;   // 1.0; // Degree
  Float_t stemax = 1.0;   // cm
  Float_t deemax = 0.1;   // 30.0; // Fraction of particle's energy 0<deemax<=1
  Float_t epsil = 1.0E-4; // 1.0; // cm
  Float_t stmin = 0.0;    // cm "Default value used"

  Float_t tmaxfdSi = 0.1;    // .10000E+01; // Degree
  Float_t stemaxSi = 0.0075; //  .10000E+01; // cm
  Float_t deemaxSi = 0.1;    // 0.30000E-02; // Fraction of particle's energy 0<deemax<=1
  Float_t epsilSi = 1.0E-4;  // .10000E+01;
  Float_t stminSi = 0.0;     // cm "Default value used"

  Float_t tmaxfdAir = 0.1;        // .10000E+01; // Degree
  Float_t stemaxAir = .10000E+01; // cm
  Float_t deemaxAir = 0.1;        // 0.30000E-02; // Fraction of particle's energy 0<deemax<=1
  Float_t epsilAir = 1.0E-4;      // .10000E+01;
  Float_t stminAir = 0.0;         // cm "Default value used"

  o2::base::Detector::Material(1, "SI$", 0.28086E+02, 0.14000E+02, 0.23300E+01, 0.93600E+01, 0.99900E+03);
  o2::base::Detector::Medium(1, "SI$", 3, 0, ifield, fieldm, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);
}

void Detector::createDetectorGeometry()
{
  // we attach everything to the volume named ‘barrel’
  auto volBarrel = gGeoManager->GetVolume("barrel");
  if (volBarrel == nullptr) {
    LOGP(fatal, "Could not get top node: volume 'barrel'!");
  }

  auto volITS = new TGeoVolumeAssembly(ITS2g::getITSVolPattern());
  volBarrel->AddNode(volITS, 2, new TGeoTranslation(0, 30., 0)); // Copy number is 2 to cheat AliGeoManager::CheckSymNamesLUT
}

void Detector::ConstructGeometry()
{
  // First the materials are created.
  createMaterials();
  // Now the actual geometry.
  createDetectorGeometry();
}

void Detector::defineSensitiveVolumes()
{
  if (gGeoManager == nullptr) {
    LOGP(fatal, "gGeoManager nullptr");
  }

  TGeoVolume* vol{nullptr};

  for (int iLayer{0}; iLayer < nTotLayers; ++iLayer) {
    if (isITS3Layer[iLayer]) {
      vol = gGeoManager->GetVolume(ITS3g::getITS3PixelArrayPattern(iLayer));
    } else {
      TString volName = ITS2g::getITSSensorPattern() + TString::Itoa(iLayer, 10);
      vol = gGeoManager->GetVolume(volName.Data());
    }

    if (vol == nullptr) {
      LOGP(fatal, "Failed to define sensitive Volume of ITS3 detector");
    }
    AddSensitiveVolume(vol);
  }
}

Hit* Detector::addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos,
                      const TVector3& startMom, double startE, double endTime, double eLoss, unsigned char startStatus,
                      unsigned char endStatus)
{
  mHits->emplace_back(trackID, detID, startPos, endPos, startMom, startE, endTime, eLoss, startStatus, endStatus);
  return &(mHits->back());
}

} // namespace o2::its3
