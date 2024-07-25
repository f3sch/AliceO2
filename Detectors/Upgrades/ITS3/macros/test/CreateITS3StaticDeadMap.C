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

/// \file CreateITS3StaticDeadMap.C
/// \brief Macro to create a static dead channel map for ITS3 from a file

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "TFile.h"
#include "TError.h"

#include "DataFormatsITSMFT/NoiseMap.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <string>
#include <vector>
#include <ranges>
#include <memory>
#endif

void CreateITS3StaticDeadMap(const std::string& filename = "input.txt")
{
  constexpr int nChips{27144};
  o2::itsmft::NoiseMap map{nChips};

  std::set<int> intSet;
  std::ifstream infile(filename);

  if (!infile) {
    Error("", "Failed to open file %s", filename.c_str());
    return;
  }

  std::string line, item;
  while (std::getline(infile, line)) {
    std::stringstream ss(line);
    for (auto item : line | std::views::split(',') | std::views::transform([](auto&& str) { return std::stoi(std::string(str.begin(), str.end())); })) {
      intSet.insert(item);
    }
  }

  for (const auto& iChip : intSet) {
    map.maskFullChip(iChip);
  }

  std::unique_ptr<TFile> f{TFile::Open("snapshot.root", "RECREATE")};
  f->WriteObjectAny(&map, "o2::itsmft::NoiseMap", "ccdb_object");
}
