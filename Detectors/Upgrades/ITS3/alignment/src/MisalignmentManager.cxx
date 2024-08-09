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

#include "Framework/Logger.h"
#include "ITS3Align/MisalignmentManager.h"
#include "ITS3Align/MisalignmentHits.h"
#include "SimConfig/DigiParams.h"

#include "TFile.h"
#include "TStopwatch.h"

#include <string>
#include <filesystem>
#include <memory>
#include <algorithm>
#include <optional>
#include <vector>

namespace fs = std::filesystem;

namespace o2::its3::align
{

void MisalignmentManager::createBackup(const fs::path& src, const fs::path& dest)
{
  if (fs::exists(dest)) {
    LOGP(info, "Previous orignal file found, using this as src");
  } else {
    if (!fs::exists(src)) {
      LOGP(fatal, "File {} does not exist", src.c_str());
    }
    LOGP(info, "Trying to backup file to {}", dest.c_str());
    try {
      fs::rename(src, dest);
    } catch (const fs::filesystem_error& err) {
      LOGP(fatal, "Cannot create backup file for Hit-File: {}", err.what());
    }
  }
}

void MisalignmentManager::misalignHits()
{
  LOGP(info, "{:*^90}", " ITS3 LOCAL MISALIGNMENT START ");

  TStopwatch timer;
  timer.Start();

  MisAlignmentHits MisAligner;
  MisAligner.init();

  const fs::path oldHitFileSrc{fs::current_path().string() + "/" + o2::conf::DigiParams::Instance().digitizationgeometry_prefix + "_HitsIT3.root"};
  const fs::path oldHitFileDest{fs::current_path().string() + "/" + o2::conf::DigiParams::Instance().digitizationgeometry_prefix + "_HitsIT3_Orig.root"};
  createBackup(oldHitFileSrc, oldHitFileDest);

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

  LOGP(info, "Preparations done; starting hit loop");
  auto nEntries = origTree->GetEntries();
  ULong64_t totalOrigHits{0}, totalNewHits{0};
  for (Long64_t iEntry{0}; origTree->LoadTree(iEntry) >= 0; ++iEntry) {
    if (origTree->GetEntry(iEntry) <= 0) {
      continue;
    }

    const auto progress = (iEntry * 100) / nEntries;
    LOG_IF(info, progress % 10 == 0) << "Processing event " << iEntry << " / " << nEntries;

    newHits.clear();
    newHits.reserve(origHits.size());
    for (const auto& origHit : origHits) {
      if (auto newHit = MisAligner.processHit(iEntry, origHit)) {
        newHits.emplace_back(*newHit);
      }
    }

    newHitsPtr = &newHits;
    newTree->Fill();

    totalNewHits += newHits.size();
    totalOrigHits += origHits.size();
  }

  newFile->WriteTObject(newTree.get());

  timer.Stop();

  MisAligner.printStats();

  auto totalDiscardedHits = totalOrigHits - totalNewHits;
  LOGP(info, "Summary: Total orignal Hits {}", totalOrigHits);
  LOGP(info, "Summary: Total misaligned Hits {} ({:.2f}%)", totalNewHits, static_cast<float>(totalNewHits) / static_cast<float>(totalOrigHits) * 100);
  LOGP(info, "Summary: Total discarded Hits {} ({:.2f}%)", totalDiscardedHits, static_cast<float>(totalDiscardedHits) / static_cast<float>(totalOrigHits) * 100);
  LOGP(info, "Summary: Misalignment took {}s", timer.CpuTime());
  LOGP(info, "{:*^90}", " ITS3 LOCAL MISALIGNMENT END ");
}

} // namespace o2::its3::align
