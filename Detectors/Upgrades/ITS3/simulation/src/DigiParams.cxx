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

/// \file DigiParams.cxx
/// \brief Implementation of the ITS3 digitization steering params

#include <fairlogger/Logger.h> // for LOG
#include "ITS3Simulation/DigiParams.h"
#include <cassert>

ClassImp(o2::its3::DigiParams);

using namespace o2::its3;

//______________________________________________
void DigiParams::print() const
{
  // print settings
  printf("ITS3 DigiParams settings:\n");
  printf("Continuous readout             : %s\n", isContinuous() ? "ON" : "OFF");
  printf("Readout Frame Length(ns)       : %f\n", getROFrameLength());
  printf("Strobe delay (ns)              : %f\n", getStrobeDelay());
  printf("Strobe length (ns)             : %f\n", getStrobeLength());
  printf("Threshold (N electrons)        : %d\n", getChargeThreshold());
  printf("Min N electrons to account     : %d\n", getMinChargeToAccount());
  printf("Number of charge sharing steps : %d\n", getNSimSteps());
  printf("ELoss to N electrons factor    : %e\n", getEnergyToNElectrons());
  printf("Noise level per pixel          : %e\n", getNoisePerPixel());
  printf("Charge time-response:\n");
  getSignalShape().print();
}
