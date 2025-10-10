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

/// \file MisalignmentParameters.h
/// \brief Definition of the MisalignmentParameters class

#ifndef ITS3_MISALIGNMENTPARAMETERS_H_
#define ITS3_MISALIGNMENTPARAMETERS_H_

#include "ITS3Base/SpecsV2.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

#include "TNamed.h"
#include "TFile.h"
#include "TMatrixD.h"

#include <array>
#include <string>

namespace o2::its3::align
{

class MisalignmentParameters : public TNamed
{
 public:
  MisalignmentParameters();
  MisalignmentParameters(const MisalignmentParameters&) = default;
  MisalignmentParameters(MisalignmentParameters&&) = delete;
  MisalignmentParameters& operator=(MisalignmentParameters&&) = delete;
  MisalignmentParameters& operator=(const MisalignmentParameters& o);

  // IO
  bool store(const std::string& file) const;
  static MisalignmentParameters* load(const std::string& file);

  /// Global getters
  const o2::detectors::AlignParam& getGlobal(unsigned int detID) const { return mGlobal[detID]; }
  void setGlobal(unsigned int detID, const o2::detectors::AlignParam& param) { mGlobal[detID] = param; }

  const TMatrixD& getLegendreCoeff(unsigned int sensorID) const { return mLegCoeff[sensorID]; }
  void setLegendreCoeff(unsigned int sensorID, const TMatrixD& m) { setMatrix(mLegCoeff[sensorID], m); }

  void printParams(unsigned int detID) const;

 private:
  inline void setMatrix(TMatrixD& o, const TMatrixD& n)
  {
    o.ResizeTo(n.GetNrows(), n.GetNcols());
    o = n;
  }

  static constexpr unsigned int nDetectors{constants::detID::nChips}; ///! for now just the IB
  std::array<o2::detectors::AlignParam, nDetectors> mGlobal;          ///< Array to hold the global misalignment
  std::array<TMatrixD, constants::nSensorsIB> mLegCoeff;              // Legendre Polynominals coefficients

  ClassDefOverride(MisalignmentParameters, 2);
};

} // namespace o2::its3::align

#endif
