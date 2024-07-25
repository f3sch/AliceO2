// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file MisAlignHits.C
/// \brief Misalign on-the-fly hits for ITS3
/// \author felix.schlepper@cern.ch

#include "ITS3Base/MisAlignHelpers.h"

void MisAlignHitsITS3(const std::string& configFilePath = "")
{
  auto& params = MisAlignGlobalParams::Instance();
  params.writeINI("default_parameters_local.ini", "MisAlignGlobalParams");
  LOGP(info, "{:*^90}", " ITS3 LOCAL MISALIGNMENT START ");
  if (configFilePath.empty()) {
    LOGP(info, "No user config provided using defaults");
  } else {
    LOGP(info, "User config at {}", configFilePath);
    params.updateFromFile(configFilePath);
  }
  params.writeINI("used_parameters_local.ini", "MisAlignGlobalParams");
  params.printKeyValues(true, true);

  std::array<Trafo3D, its3c::nSensorsIB> gRotoTranslations{};
  for (int iSensor{0}; iSensor < (int)its3c::nSensorsIB; ++iSensor) {
    Euler3D rot{
      params.getValueAs<float>(fmt::format("MisAlignGlobalParams.Sensor{}Phi", iSensor)),
      params.getValueAs<float>(fmt::format("MisAlignGlobalParams.Sensor{}Theta", iSensor)),
      params.getValueAs<float>(fmt::format("MisAlignGlobalParams.Sensor{}Psi", iSensor)),
    };
    Trans3D trans{
      params.getValueAs<float>(fmt::format("MisAlignGlobalParams.Sensor{}Dx", iSensor)),
      params.getValueAs<float>(fmt::format("MisAlignGlobalParams.Sensor{}Dy", iSensor)),
      params.getValueAs<float>(fmt::format("MisAlignGlobalParams.Sensor{}Dz", iSensor)),
    };
    gRotoTranslations[iSensor] = Trafo3D(rot, trans);
  }

  const fs::path oldHitFileSrc{fs::current_path().string() + "/o2sim_HitsIT3.root"};
  const fs::path oldHitFileDest{fs::current_path().string() + "/o2sim_HitsIT3_Orig.root"};
  backup_file(oldHitFileSrc, oldHitFileDest);

  std::unique_ptr<TFile> origFile{TFile::Open(oldHitFileDest.c_str(), "READ")};
  if (origFile == nullptr || origFile->IsZombie()) {
    LOGP(fatal, "Original file {} cannot be opened", oldHitFileDest.c_str());
  }

  std::unique_ptr<TTree> origTree{origFile->Get<TTree>("o2sim")};
  if (origTree == nullptr) {
    LOGP(fatal, "Cannot get hit-tree from orignal file");
  }
  std::vector<o2::itsmft::Hit> origHits, *origHitsPtr{&origHits};
  origTree->SetBranchAddress("IT3Hit", &origHitsPtr);

  std::unique_ptr<TFile> newFile{TFile::Open(oldHitFileSrc.c_str(), "RECREATE")};
  if (newFile == nullptr || newFile->IsZombie()) {
    LOGP(fatal, "New file {} cannot be opened", oldHitFileSrc.c_str());
  }

  auto newTree = std::make_unique<TTree>("o2sim", "o2sim");
  std::vector<o2::itsmft::Hit> newHits, *newHitsPtr{nullptr};
  newTree->Branch("IT3Hit", &newHitsPtr);

  auto misalignHit = [&](const o2::itsmft::Hit& origHit) -> std::optional<o2::itsmft::Hit> {
    if (!its3c::detID::isDetITS3(origHit.GetDetectorID())) { // only modify IB
      return origHit;
    }
    int sensorID{its3c::detID::getSensorID(origHit.GetDetectorID())};

    auto hit{origHit};
    Point3D hitXYZ{hit.GetX(), hit.GetY(), hit.GetZ()};

    // Global RotoTranslation
    hitXYZ = gRotoTranslations[sensorID] * hitXYZ;
    hit.SetXYZ(hitXYZ.X(), hitXYZ.Y(), hitXYZ.Z());
    return hit;
  };

  ULong64_t totalOrigHits{0}, totalNewHits{0};
  for (int iEntry{0}; origTree->LoadTree(iEntry) >= 0; ++iEntry) {
    if (origTree->GetEntry(iEntry) <= 0) {
      continue;
    }

    newHits.clear();
    newHits.reserve(origHits.size());
    for (const auto& origHit : origHits) {
      if (auto newHit = misalignHit(origHit)) {
        newHits.emplace_back(*newHit);
      }
    }

    newHitsPtr = &newHits;
    newTree->Fill();

    totalNewHits += newHits.size();
    totalOrigHits += origHits.size();
  }

  newFile->WriteTObject(newTree.get());

  LOGP(info, "Summary: Total orignal Hits {}", totalOrigHits);
  LOGP(info, "Summary: Total misaligned Hits {}", totalNewHits);
  LOGP(info, "Summary: Discarded Hits {}", totalOrigHits - totalNewHits);

  LOGP(info, "{:*^90}", " ITS3 LOCAL MISALIGNMENT END ");
}
