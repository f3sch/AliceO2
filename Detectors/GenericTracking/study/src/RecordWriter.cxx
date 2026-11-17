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

#include "GenericTracking/Study/RecordWriter.h"

#include <TFile.h>

#include <stdexcept>

namespace o2::gtrk::study
{

RecordWriter::RecordWriter(const std::string& path) : mPath(path)
{
  if (path.ends_with(".csv")) {
    mCsv = std::make_unique<std::ofstream>(path);
    if (!*mCsv) {
      throw std::runtime_error("cannot write " + path);
    }
  } else if (path.ends_with(".root")) {
    mTree = std::make_unique<o2::utils::TreeStreamRedirector>(path.c_str(), "recreate");
    const auto* file = mTree->GetFile();
    if (file == nullptr || file->IsZombie() || !file->IsWritable()) {
      throw std::runtime_error("cannot write " + path);
    }
  } else {
    throw std::invalid_argument("a study's output is a .csv or a .root file, not " + path);
  }
}

RecordWriter::~RecordWriter()
{
  if (mTree) {
    mTree->Close();
  }
}

} // namespace o2::gtrk::study
