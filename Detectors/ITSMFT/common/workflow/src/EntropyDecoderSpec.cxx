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

/// @file   EntropyDecoderSpec.cxx

#include <vector>

#include "Framework/ControlService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/CCDBParamSpec.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "ITSMFTWorkflow/EntropyDecoderSpec.h"
#include "ITSMFTReconstruction/ClustererParam.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "DataFormatsITSMFT/PhysTrigger.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"

using namespace o2::framework;

namespace o2
{
namespace itsmft
{

template <int N>
EntropyDecoderSpec<N>::EntropyDecoderSpec(int verbosity, bool getDigits)
  : mCTFCoder(o2::ctf::CTFCoderBase::OpType::Decoder, N), mGetDigits(getDigits)
{
  assert(orig == o2::header::gDataOriginITS || orig == o2::header::gDataOriginMFT);
  mDetPrefix = Origin == o2::header::gDataOriginITS ? "_ITS" : "_MFT";
  mTimer.Stop();
  mTimer.Reset();
  mCTFCoder.setVerbosity(verbosity);
  mCTFCoder.setDictBinding(std::string("ctfdict") + mDetPrefix);
}

template <int N>
void EntropyDecoderSpec<N>::init(o2::framework::InitContext& ic)
{
  mCTFCoder.init<CTF>(ic);
  mMaskNoise = ic.options().get<bool>("mask-noise");
  mUseClusterDictionary = !ic.options().get<bool>("ignore-cluster-dictionary");
}

template <int N>
void EntropyDecoderSpec<N>::run(ProcessingContext& pc)
{
  if (pc.services().get<o2::framework::TimingInfo>().globalRunNumberChanged) {
    mTimer.Reset();
  }
  auto cput = mTimer.CpuTime();
  mTimer.Start(false);
  o2::ctf::CTFIOSize iosize;
  updateTimeDependentParams(pc);
  auto buff = pc.inputs().get<gsl::span<o2::ctf::BufferType>>(std::string("ctf") + mDetPrefix);
  // since the buff is const, we cannot use EncodedBlocks::relocate directly, instead we wrap its data to another flat object
  //  const auto ctfImage = o2::itsmft::CTF::getImage(buff.data());

  // this produces weird memory problems in unrelated devices, to be understood
  // auto& trigs = pc.outputs().make<std::vector<o2::itsmft::PhysTrigger>>(OutputRef{"phystrig"}); // dummy output

  if constexpr (DPLAlpideParam<N>::supportsStaggering()) {
    // for now we need to 'mock' the staggered output and sort ordering ourselves
    std::vector<o2::itsmft::ROFRecord> rofs;
    std::vector<o2::itsmft::Digit> digits;
    std::vector<o2::itsmft::CompClusterExt> clusters;
    std::vector<unsigned char> patterns;
    // do the actual read
    if (mGetDigits) {
      if (buff.size()) {
        iosize = mCTFCoder.decode(o2::itsmft::CTF::getImage(buff.data()), rofs, digits, mNoiseMap, mPattIdConverter);
      }
      mTimer.Stop();
      LOG(info) << "Decoded " << digits.size() << " digits in " << rofs.size() << " RO frames, (" << iosize.asString() << ") in " << mTimer.CpuTime() - cput << " s";
    } else {
      if (buff.size()) {
        iosize = mCTFCoder.decode(o2::itsmft::CTF::getImage(buff.data()), rofs, clusters, patterns, mNoiseMap, mPattIdConverter);
      }
      mTimer.Stop();
      LOG(info) << "Decoded " << clusters.size() << " clusters in " << rofs.size() << " RO frames, (" << iosize.asString() << ") in " << mTimer.CpuTime() - cput << " s";
    }
    std::array<std::vector<o2::itsmft::ROFRecord>, NLayers> rofsPerLayer;
    std::array<std::vector<o2::itsmft::Digit>, NLayers> digitsPerLayer;
    std::array<std::vector<o2::itsmft::CompClusterExt>, NLayers> clustersPerLayer;
    std::array<std::vector<unsigned char>, NLayers> patternsPerLayer;
    std::array<std::vector<int>, NLayers> firstEntries;
    std::array<std::vector<int>, NLayers> nEntries;
    for (int iLayer{0}; iLayer < NLayers; ++iLayer) {
      rofsPerLayer[iLayer] = rofs;
      firstEntries[iLayer].resize(rofs.size(), 0);
      nEntries[iLayer].resize(rofs.size(), 0);
    }
    // now we need to filter the data per layer
    // TODO implement also for cluster input
    for (size_t iROF{0}; iROF < rofs.size(); ++iROF) {
      const auto& rof = rofs[iROF];
      for (int iEntry{rof.getFirstEntry()}; iEntry < (rof.getFirstEntry() + rof.getNEntries()); ++iEntry) {
        const auto& dig = digits[iEntry];
        int lay = ChipMappingITS::getLayer(dig.getChipIndex());
        digitsPerLayer[lay].push_back(dig);
        ++(nEntries[lay][iROF]);
      }
    }
    for (int iLayer{0}; iLayer < NLayers; ++iLayer) {
      std::exclusive_scan(nEntries[iLayer].begin(), nEntries[iLayer].end(), firstEntries[iLayer].begin(), 0);
      for (int iROF{0}; iROF < rofs.size(); ++iROF) {
        rofsPerLayer[iLayer][iROF].setFirstEntry(firstEntries[iLayer][iROF]);
        rofsPerLayer[iLayer][iROF].setNEntries(nEntries[iLayer][iROF]);
      }
    }
    for (uint32_t iLayer{0}; iLayer < NLayers; ++iLayer) {
      pc.outputs().snapshot(OutputRef{"ROframes", iLayer}, rofsPerLayer[iLayer]);
      if (mGetDigits) {
        pc.outputs().snapshot(OutputRef{"Digits", iLayer}, digitsPerLayer[iLayer]);
      } else {
        pc.outputs().snapshot(OutputRef{"compClusters", iLayer}, clustersPerLayer[iLayer]);
        pc.outputs().snapshot(OutputRef{"patterns", iLayer}, patternsPerLayer[iLayer]);
      }
    }
  } else {
    auto& rofs = pc.outputs().make<std::vector<o2::itsmft::ROFRecord>>(OutputRef{"ROframes", 0});
    if (mGetDigits) {
      auto& digits = pc.outputs().make<std::vector<o2::itsmft::Digit>>(OutputRef{"Digits", 0});
      if (buff.size()) {
        iosize = mCTFCoder.decode(o2::itsmft::CTF::getImage(buff.data()), rofs, digits, mNoiseMap, mPattIdConverter);
      }
      mTimer.Stop();
      LOG(info) << "Decoded " << digits.size() << " digits in " << rofs.size() << " RO frames, (" << iosize.asString() << ") in " << mTimer.CpuTime() - cput << " s";
    } else {
      auto& compcl = pc.outputs().make<std::vector<o2::itsmft::CompClusterExt>>(OutputRef{"compClusters", 0});
      auto& patterns = pc.outputs().make<std::vector<unsigned char>>(OutputRef{"patterns", 0});
      if (buff.size()) {
        iosize = mCTFCoder.decode(o2::itsmft::CTF::getImage(buff.data()), rofs, compcl, patterns, mNoiseMap, mPattIdConverter);
      }
      mTimer.Stop();
      LOG(info) << "Decoded " << compcl.size() << " clusters in " << rofs.size() << " RO frames, (" << iosize.asString() << ") in " << mTimer.CpuTime() - cput << " s";
    }
    // hack: output empty messages to avoid dropping the TF
    for (uint32_t iLayer{1}; iLayer < NLayers; ++iLayer) {
      pc.outputs().make<std::vector<o2::itsmft::ROFRecord>>(OutputRef{"ROframes", iLayer});
      if (mGetDigits) {
        pc.outputs().make<std::vector<o2::itsmft::Digit>>(OutputRef{"Digits", iLayer});
      } else {
        pc.outputs().make<std::vector<o2::itsmft::CompClusterExt>>(OutputRef{"compClusters", iLayer});
        pc.outputs().make<std::vector<unsigned char>>(OutputRef{"patterns", iLayer});
      }
    }
  }
  pc.outputs().snapshot({"ctfrep", 0}, iosize);
} // namespace itsmft

template <int N>
void EntropyDecoderSpec<N>::endOfStream(EndOfStreamContext& ec)
{
  LOGF(info, "%s Entropy Decoding total timing: Cpu: %.3e Real: %.3e s in %d slots",
       Origin.as<std::string>(), mTimer.CpuTime(), mTimer.RealTime(), mTimer.Counter() - 1);
}

template <int N>
void EntropyDecoderSpec<N>::updateTimeDependentParams(ProcessingContext& pc)
{
  if (pc.services().get<o2::framework::TimingInfo>().globalRunNumberChanged) { // this params need to be queried only once
    if (mMaskNoise) {
      pc.inputs().get<o2::itsmft::NoiseMap*>(std::string("noise") + mDetPrefix);
    }
    if (mGetDigits || mMaskNoise) {
      pc.inputs().get<o2::itsmft::TopologyDictionary*>(std::string("cldict") + mDetPrefix);
    }
  }
  mCTFCoder.updateTimeDependentParams(pc, true);
}

template <int N>
void EntropyDecoderSpec<N>::finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj)
{
  if (matcher == ConcreteDataMatcher(Origin, "NOISEMAP", 0)) {
    mNoiseMap = (o2::itsmft::NoiseMap*)obj;
    LOG(info) << Origin.as<std::string>() << " noise map updated";
    return;
  }
  if (matcher == ConcreteDataMatcher(Origin, "CLUSDICT", 0)) {
    LOG(info) << Origin.as<std::string>() << " cluster dictionary updated" << (!mUseClusterDictionary ? " but its using is disabled" : "");
    mPattIdConverter.setDictionary((const TopologyDictionary*)obj);
    return;
  }
  if (mCTFCoder.finaliseCCDB<CTF>(matcher, obj)) {
    return;
  }
}

template <int N>
DataProcessorSpec getEntropyDecoderSpec(int verbosity, bool getDigits, unsigned int sspec)
{
  using EntropyDecoder = EntropyDecoderSpec<N>;

  std::string det = EntropyDecoder::Origin.template as<std::string>();
  std::string nm = "_" + det;
  std::vector<InputSpec> inputs;
  inputs.emplace_back(std::string("ctf") + nm, EntropyDecoder::Origin, "CTFDATA", sspec, Lifetime::Timeframe);
  inputs.emplace_back(std::string("noise") + nm, EntropyDecoder::Origin, "NOISEMAP", 0, Lifetime::Condition, ccdbParamSpec(fmt::format("{}/Calib/NoiseMap", det)));
  inputs.emplace_back(std::string("cldict") + nm, EntropyDecoder::Origin, "CLUSDICT", 0, Lifetime::Condition, ccdbParamSpec(fmt::format("{}/Calib/ClusterDictionary", det)));
  inputs.emplace_back(std::string("ctfdict") + nm, EntropyDecoder::Origin, "CTFDICT", 0, Lifetime::Condition, ccdbParamSpec(fmt::format("{}/Calib/CTFDictionaryTree", det)));
  inputs.emplace_back(std::string("trigoffset"), "CTP", "Trig_Offset", 0, Lifetime::Condition, ccdbParamSpec("CTP/Config/TriggerOffsets"));

  std::vector<OutputSpec> outputs;
  // this is a special dummy input which makes sense only in sync workflows

  // this produces weird memory problems in unrelated devices, to be understood
  // outputs.emplace_back(OutputSpec{{"phystrig"}, orig, "PHYSTRIG", 0, Lifetime::Timeframe});

  for (uint32_t iLayer{0}; iLayer < EntropyDecoder::NLayers; ++iLayer) {
    if (getDigits) {
      outputs.emplace_back(OutputSpec{{"Digits"}, EntropyDecoder::Origin, "DIGITS", iLayer, Lifetime::Timeframe});
      outputs.emplace_back(OutputSpec{{"ROframes"}, EntropyDecoder::Origin, "DIGITSROF", iLayer, Lifetime::Timeframe});
    } else {
      outputs.emplace_back(OutputSpec{{"compClusters"}, EntropyDecoder::Origin, "COMPCLUSTERS", iLayer, Lifetime::Timeframe});
      outputs.emplace_back(OutputSpec{{"ROframes"}, EntropyDecoder::Origin, "CLUSTERSROF", iLayer, Lifetime::Timeframe});
      outputs.emplace_back(OutputSpec{{"patterns"}, EntropyDecoder::Origin, "PATTERNS", iLayer, Lifetime::Timeframe});
    }
  }
  outputs.emplace_back(OutputSpec{{"ctfrep"}, EntropyDecoder::Origin, "CTFDECREP", 0, Lifetime::Timeframe});

  return DataProcessorSpec{
    .name = EntropyDecoder::DeviceName,
    .inputs = inputs,
    .outputs = outputs,
    .algorithm = AlgorithmSpec{adaptFromTask<EntropyDecoder>(verbosity, getDigits)},
    .options = Options{
      {"ctf-dict", VariantType::String, "ccdb", {"CTF dictionary: empty or ccdb=CCDB, none=no external dictionary otherwise: local filename"}},
      {"mask-noise", VariantType::Bool, false, {"apply noise mask to digits or clusters (involves reclusterization)"}},
      {"ignore-cluster-dictionary", VariantType::Bool, false, {"do not use cluster dictionary, always store explicit patterns"}},
      {"and-version", VariantType::String, {"version of and entropy coder implementation to use"}}}};
}

framework::DataProcessorSpec getITSEntropyDecoderSpec(int verbosity, bool getDigits, unsigned int sspec) { return getEntropyDecoderSpec<o2::detectors::DetID::ITS>(verbosity, getDigits, sspec); }
framework::DataProcessorSpec getMFTEntropyDecoderSpec(int verbosity, bool getDigits, unsigned int sspec) { return getEntropyDecoderSpec<o2::detectors::DetID::MFT>(verbosity, getDigits, sspec); }

} // namespace itsmft
} // namespace o2
