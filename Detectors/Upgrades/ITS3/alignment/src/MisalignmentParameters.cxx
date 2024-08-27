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

/// \file MisalignmentParameter.cxx
/// \brief Implementation of the MisalignmentParameter class

#include "ITS3Align/MisalignmentParameters.h"
#include "Framework/Logger.h"

#include "TFile.h"

#include <memory>

ClassImp(o2::its3::align::MisalignmentParameters);

namespace o2::its3::align
{

MisalignmentParameters::MisalignmentParameters()
{
  SetName("MisalignmentParameters");
  SetTitle("ITS3 MisalignmentParameters");
}

MisalignmentParameters& MisalignmentParameters::operator=(const MisalignmentParameters& o)
{
  if (this != &o) {
    SetName(o.GetName());
    SetTitle(o.GetTitle());
    mGlobal = o.mGlobal;
    for (unsigned int s{0}; s < constants::nSensorsIB; ++s) {
      o.mLegCoeff[s].Print();
      setMatrix(mLegCoeff[s], o.mLegCoeff[s]);
    }
  }
  return *this;
}

bool MisalignmentParameters::store(const std::string& file) const
{
  std::unique_ptr<TFile> fOut(TFile::Open(file.c_str(), "RECREATE"));
  if (fOut == nullptr || fOut->IsZombie()) {
    LOGP(info, "Unable to save misalignment parameters");
    return false;
  }
  fOut->WriteObjectAny(this, "o2::its3::align::MisalignmentParameters", "ccdb_object");
  return true;
}

MisalignmentParameters* MisalignmentParameters::load(const std::string& file)
{
  LOGP(info, "Loading Parameters from {}", file);
  std::unique_ptr<TFile> fIn(TFile::Open(file.c_str(), "READ"));
  auto p = fIn->Get<MisalignmentParameters>("ccdb_object");
  if (p == nullptr) {
    LOGP(fatal, "Unable to load parameters from file!");
  }
  return p;
}

void MisalignmentParameters::printParams(unsigned int detID) const
{
  LOGP(info, "Parameters for ID={}:", detID);
  const auto& glo = mGlobal[detID];
  LOGP(info, " - Global Trans: X={} Y={} Z={}", glo.getX(), glo.getY(), glo.getZ());
  LOGP(info, " - Global Rots: Phi={} Psi={} Theta={}", glo.getPhi(), glo.getPsi(), glo.getTheta());
  if (constants::detID::isDetITS3(detID)) {
    auto sensorID = constants::detID::getSensorID(detID);
    mLegCoeff[sensorID].Print();
  }
}

} // namespace o2::its3::align
