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

#ifndef O2_CTP_RAWDECODER_H
#define O2_CTP_RAWDECODER_H

#include <vector>
#include <deque>

#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "DataFormatsCTP/Digits.h"
#include "DataFormatsCTP/LumiInfo.h"

namespace o2
{
namespace ctp
{
namespace reco_workflow
{

/// \class RawDecoderSpec
/// \brief Coverter task for Raw data to CTP digits
/// \author Roman Lietava from CPV example
///
class RawDecoderSpec : public framework::Task
{
 public:
  /// \brief Constructor
  /// \param propagateMC If true the MCTruthContainer is propagated to the output
  RawDecoderSpec(bool digits, bool lumi) : mDoDigits(digits), mDoLumi(lumi) {}

  /// \brief Destructor
  ~RawDecoderSpec() override = default;

  /// \brief Initializing the RawDecoderSpec
  /// \param ctx Init context
  void init(framework::InitContext& ctx) final;

  /// \brief Run conversion of raw data to cells
  /// \param ctx Processing context
  ///
  /// The following branches are linked:
  /// Input RawData: {"ROUT", "RAWDATA", 0, Lifetime::Timeframe}
  /// Output HW errors: {"CTP", "RAWHWERRORS", 0, Lifetime::Timeframe} -later
  void run(framework::ProcessingContext& ctx) final;
  static void makeGBTWordInverse(std::vector<gbtword80_t>& diglets, gbtword80_t& GBTWord, gbtword80_t& remnant, uint32_t& size_gbt, uint32_t Npld);

 protected:
 private:
  // for digits
  bool mDoDigits = true;
  std::vector<CTPDigit> mOutputDigits;
  // for lumi
  bool mDoLumi = true;
  gbtword80_t mTVXMask = 0x4; // TVX is 3rd input
  LumiInfo mOutputLumiInfo;
  bool mVerbose = false;
  uint64_t mCounts = 0;
  uint32_t mNTFToIntegrate = 1;
  uint32_t mNHBIntegrated = 0;
  std::deque<size_t> mHistory;
};

/// \brief Creating DataProcessorSpec for the CTP
///
o2::framework::DataProcessorSpec getRawDecoderSpec(bool askSTFDist, bool digits, bool lumi);

} // namespace reco_workflow

} // namespace ctp

} // namespace o2

#endif