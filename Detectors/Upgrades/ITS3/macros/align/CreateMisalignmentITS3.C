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

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "TRandom.h"
#include "TMatrixD.h"

#include "ITS3Align/MisalignmentParameters.h"
#endif

void CreateMisalignmentITS3(bool dummy = false)
{
  gRandom->SetSeed(42);

  // Legendre coeff.
  constexpr int nOrder{2};
  auto getRandom = []() {
    constexpr double scale{100.e-4};
    return scale * gRandom->Uniform(-1.0, 1.0);
  };

  o2::its3::align::MisalignmentParameters params;

  for (int sensorID{0}; sensorID < 6; ++sensorID) {
    TMatrixD coeffX(nOrder + 1, nOrder + 1);
    TMatrixD coeffY(nOrder + 1, nOrder + 1);
    TMatrixD coeffZ(nOrder + 1, nOrder + 1);
    for (int i{0}; i <= nOrder; ++i) {
      for (int j{0}; j <= i; ++j) {
        // some random scaling as higher order parameters have higher influence
        coeffX(i, j) = (dummy) ? 0.0 : getRandom() / (1.0 + i * j * 2.0);
        coeffZ(i, j) = (dummy) ? 0.0 : getRandom() / (1.0 + i * j * 2.0);
        coeffY(i, j) = (dummy) ? 0.0 : getRandom() / (1.0 + i * j * 2.0);
      }
    }

    params.setLegendreCoeffX(sensorID, coeffX);
    params.setLegendreCoeffY(sensorID, coeffY);
    params.setLegendreCoeffZ(sensorID, coeffZ);
    params.printLegendreParams(sensorID);
  }

  params.store("misparams.root");
}
