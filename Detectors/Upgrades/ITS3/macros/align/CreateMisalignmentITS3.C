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
#include "DetectorsCommonDataFormats/AlignParam.h"
#endif

void CreateMisalignmentITS3(bool dummy = false, bool manual = false)
{
  gRandom->SetSeed(42);

  // Legendre coeff.
  constexpr int nOrder{2};
  auto getRandom = []() {
    constexpr double scale{50.e-4};
    return scale * gRandom->Uniform(-1.0, 1.0);
  };

  auto getSign = []() { return gRandom->Uniform() ? -1.0 : 1.0; };

  o2::its3::align::MisalignmentParameters params;

  for (int sensorID{0}; sensorID < 6; ++sensorID) {
    o2::detectors::AlignParam p;
    p.setTranslation(getRandom(), getRandom(), getRandom());
    params.setGlobal(sensorID, p);
  }

  if (dummy) {
    const TMatrixD coeffNull(0 + 1, 0 + 1);
    for (int sensorID{0}; sensorID < 6; ++sensorID) {
      params.setLegendreCoeff(sensorID, coeffNull);
    }
  } else if (manual) {
    for (int sensorID{0}; sensorID < 6; ++sensorID) {
      constexpr double scale{100e-4};
      int nOrder{2};
      TMatrixD coeffMinus(nOrder + 1, nOrder + 1);
      TMatrixD coeffPlus(nOrder + 1, nOrder + 1);
      /* coeffMinus(0, 0) = -scale; */
      coeffMinus(2, 2) = -scale;
      coeffPlus(2, 2) = -scale;
      if (sensorID % 2 == 0) {
        params.setLegendreCoeff(sensorID, coeffPlus);
      } else {
        params.setLegendreCoeff(sensorID, coeffMinus);
      }
    }
  } else {
    for (int sensorID{0}; sensorID < 6; ++sensorID) {
      TMatrixD coeff(nOrder + 1, nOrder + 1);
      for (int i{0}; i <= nOrder; ++i) {
        for (int j{0}; j <= i; ++j) {
          // some random scaling as higher order parameters have higher influence
          coeff(i, j) = getRandom() / (1.0 + i * j * 2.0);
        }
      }
      params.setLegendreCoeff(sensorID, coeff);
    }
  }

  for (int sensorID{0}; sensorID < 6; ++sensorID) {
    params.getLegendreCoeff(sensorID).Print();
  }

  params.store("misparams.root");
}
