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

ClassImp(o2::its3::Detector);

namespace o2::its3
{

Detector::Detector(bool isActive) : o2::base::DetImpl<Detector>("IT3", isActive) {}

Detector::~Detector()
{
  if (mHits != nullptr) {
    o2::utils::freeSimVector(mHits);
  }
}

void Detector::Register()
{
  // This will create a branch in the output tree called Hit, setting the last
  // parameter to kFALSE means that this collection will not be written to the file,
  // it will exist only during the simulation

  if (FairRootManager::Instance()) {
    FairRootManager::Instance()->RegisterAny(addNameTo("Hit").data(), mHits, kTRUE);
  }
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

} // namespace o2::its3
