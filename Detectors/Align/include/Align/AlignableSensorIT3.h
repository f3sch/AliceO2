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

/// @file   AlignableSensorIT3.h
/// @brief  IT3 sensor

#ifndef ALIGNABLESENSORIT3_H
#define ALIGNABLESENSORIT3_H

#include "Align/AlignableSensor.h"

namespace o2::align
{

class AlignableSensorIT3 final : public AlignableSensor
{
 public:
  AlignableSensorIT3() = default;
  AlignableSensorIT3(const AlignableSensorIT3&) = default;
  AlignableSensorIT3(AlignableSensorIT3&&) = delete;
  AlignableSensorIT3& operator=(const AlignableSensorIT3&) = default;
  AlignableSensorIT3& operator=(AlignableSensorIT3&&) = delete;
  AlignableSensorIT3(const char* name, int vid, int iid, Controller* ctr);
  ~AlignableSensorIT3() final = default;
  void prepareMatrixT2L() final;

 protected:

  ClassDefOverride(AlignableSensorIT3, 1)
};

} // namespace o2
#endif
