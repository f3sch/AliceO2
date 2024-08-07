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

#ifndef ITS3_MISALIGNMENTMANAGER_H_
#define ITS3_MISALIGNMENTMANAGER_H_

#include <filesystem>

namespace o2::its3::align
{

class MisalignmentManager
{
 public:
  void misalignHits();

  static void createBackup(const std::filesystem::path& src, const std::filesystem::path& dest);
 private:
};

} // namespace o2::its3::align

#endif
