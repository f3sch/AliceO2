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

/// \file   TRDCalibGainReaderSpec.h
/// \brief Calibration Reader for gain
/// \author Felix Schlepper

#ifndef O2_TRD_CALIBGAINREADER
#define O2_TRD_CALIBGAINREADER

#include "TFile.h"
#include "TTree.h"

#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"

namespace o2
{
namespace trd
{

class TRDCalibGainReader : public o2::framework::Task
{
 public:
  TRDCalibGainReader() = default;
  ~TRDCalibGainReader() override = default;
  void init(o2::framework::InitContext& ic) final;
  void run(o2::framework::ProcessingContext& pc) final;

 private:
  void connectTree();
  std::unique_ptr<TFile> mFile;
  std::unique_ptr<TTree> mTree;
  std::string mInFileName{"trdgain.root"};
  std::string mInTreeName{"calibdata"};
};

/// create a processor spec
/// read TRD calibration data from a root file
framework::DataProcessorSpec getTRDCalibGainReaderSpec();

} // namespace trd
} // namespace o2

#endif /* O2_TRD_CALIBGAINREADER */
