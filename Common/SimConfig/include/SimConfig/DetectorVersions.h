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

#ifndef O2_DETECTORVERSIONS_H_
#define O2_DETECTORVERSIONS_H_

#include <string>
#include <unordered_map>
#include <vector>

#include "Framework/Logger.h"

namespace o2::conf
{
using DetectorElements_t = std::vector<std::string>;
using DetectorMap_t = std::unordered_map<std::string, DetectorElements_t>;

// Container defining different general evolutions of the ALICE experiment. Each
// evolution is given a name and a list defining the names of the detectors and
// passive elements present.
extern const DetectorMap_t DetectorVersions;

void printDetMap(const DetectorMap_t& map);
} // namespace o2::conf

#endif // O2_DETECTORVERSIONS_H_
