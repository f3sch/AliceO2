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
///
/// \file Definitions.h
/// \brief

#ifndef TRACKINGITS_DEFINITIONS_H_
#define TRACKINGITS_DEFINITIONS_H_

#include <type_traits>
#include <cstdint>
#include <tuple>

#include "SimulationDataFormat/MCCompLabel.h"
#include "CommonDataFormat/TimeStamp.h"
#include "ReconstructionDataFormats/Vertex.h"

namespace o2::its
{

template <bool IsConst, typename T>
using maybe_const = std::conditional_t<IsConst, const T, T>;

// Time estimates are given in BC
// error needs to cover maximum 1 orbit
// this is an asymmetric error defining an interval [time, time+error)
using TimeEstBC = o2::dataformats::TimeStampWithError<uint32_t, uint16_t>;
using Vertex = o2::dataformats::Vertex<TimeEstBC>;
// MC vertex label with purity
using VertexLabel = std::pair<o2::MCCompLabel, float>;

// simple implemnetion of logging with exp. backoff
struct LogLogThrottler {
  uint64_t evCount{0};
  uint64_t nextLog{1};
  int32_t iteration{-1};
  int32_t layer{-1};
  bool needToLog(int32_t iter, int32_t lay)
  {
    if (iteration != iter || layer != lay) {
      iteration = iter;
      layer = lay;
      evCount = 0;
      nextLog = 1;
    }
    if (++evCount > nextLog) {
      nextLog *= 2;
      return true;
    }
    return false;
  }
};

} // namespace o2::its

#endif
