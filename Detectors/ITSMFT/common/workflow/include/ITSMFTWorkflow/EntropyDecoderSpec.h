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

/// @file   EntropyDecoderSpec.h
/// @brief  Convert CTF (EncodedBlocks) to clusters streams

#ifndef O2_ITSMFT_ENTROPYDECODER_SPEC
#define O2_ITSMFT_ENTROPYDECODER_SPEC

#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "Headers/DataHeader.h"
#include "ITSMFTReconstruction/CTFCoder.h"
#include "DataFormatsITSMFT/NoiseMap.h"
#include "ITSMFTBase/DPLAlpideParam.h"
#include "ITSMFTReconstruction/LookUp.h"
#include <TStopwatch.h>

namespace o2::itsmft
{

template <int N>
class EntropyDecoderSpec : public o2::framework::Task
{
 public:
  static constexpr o2::header::DataOrigin Origin{N == o2::detectors::DetID::ITS ? o2::header::gDataOriginITS : o2::header::gDataOriginMFT};
  static constexpr int NLayers{o2::itsmft::DPLAlpideParam<N>::getNLayers()};
  static constexpr const char* DeviceName{N == o2::detectors::DetID::ITS ? "its-entropy-decoder" : "mft-entropy-decoder"};

  EntropyDecoderSpec(int verbosity, bool getDigits = false);
  ~EntropyDecoderSpec() override = default;
  void init(o2::framework::InitContext& ic) final;
  void run(o2::framework::ProcessingContext& pc) final;
  void endOfStream(o2::framework::EndOfStreamContext& ec) final;
  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj) final;

 private:
  void updateTimeDependentParams(o2::framework::ProcessingContext& pc);

  o2::itsmft::CTFCoder mCTFCoder;
  const NoiseMap* mNoiseMap = nullptr;
  LookUp mPattIdConverter;
  bool mGetDigits{false};
  bool mMaskNoise{false};
  bool mUseClusterDictionary{true};
  std::string mDetPrefix;
  std::string mCTFDictPath;
  TStopwatch mTimer;
};

framework::DataProcessorSpec getITSEntropyDecoderSpec(int verbosity, bool getDigits, unsigned int sspec);
framework::DataProcessorSpec getMFTEntropyDecoderSpec(int verbosity, bool getDigits, unsigned int sspec);

} // namespace o2::itsmft

#endif
