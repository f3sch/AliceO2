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

#include "SimConfig/DetectorVersions.h"

namespace o2::conf
{

const DetectorMap_t DetectorVersions{
#ifdef ENABLE_UPGRADES
  {"ALICE3", {
               // clang-format off
               // Active
               "TRK",
               "FT3",
               "FCT",
               "RCH",
               "MI3",
               "ECL",
               // Passive
               "HALL",
               "MAG",
               "A3IP",
               "A3ABSO",
               "A3MAG",
               // clang-format on
             }},
  // ALICE 2.1
  {"ALICE2.1", {
                 // clang-format off
                 // Active
                 "IT3",
                 "TPC",
                 "TRD",
                 "TOF",
                 "PHS",
                 "EMC",
                 "HMP",
                 "MFT",
                 "MCH",
                 "MID",
                 "ZDC",
                 "FT0",
                 "FV0",
                 "FDD",
                 "CTP",
                 "FOC",
                 // Passive
                 "HALL",
                 "MAG",
                 "DIPO",
                 "COMP",
                 "PIPE",
                 "ABSO",
                 "SHIL"
                 // clang-format on
                 // Note for IT3 the pipe is automatically replaced with another version
               }},
#endif
  // ALICE 2
  {"ALICE2", {
               // clang-format off
               // Active
               "ITS",
               "TPC",
               "TRD",
               "TOF",
               "PHS",
               "EMC",
               "HMP",
               "MFT",
               "MCH",
               "MID",
               "ZDC",
               "FT0",
               "FV0",
               "FDD",
               "CTP",
               // Passive
               "HALL",
               "MAG",
               "DIPO",
               "COMP",
               "PIPE",
               "ABSO",
               "SHIL",
               // clang-format on
             }},
};

void printDetMap(const DetectorMap_t& map)
{
  LOGP(error, "List of available versions including their detectors:");
  for (int i{0}; const auto& [version, elements] : map) {
    LOGP(error, " - {: >2d}. {}:", i++, version);
    for (int j{0}; const auto& element : elements) {
      LOGP(error, "\t\t* {: >2d}.\t{}", j++, element);
    }
  }
}

} // namespace o2::conf
