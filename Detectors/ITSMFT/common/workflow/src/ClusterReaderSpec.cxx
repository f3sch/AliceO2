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

/// @file   ClusterReaderSpec.cxx

#include <vector>
#include <cassert>

#include <TTree.h>

#include "Framework/ControlService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/Logger.h"
#include "ITSMFTWorkflow/ClusterReaderSpec.h"
#include "ITSMFTBase/DPLAlpideParam.h"
#include "DataFormatsITSMFT/PhysTrigger.h"
#include "CommonUtils/NameConf.h"

using namespace o2::framework;
using namespace o2::itsmft;

namespace o2
{
namespace itsmft
{

template <int N>
ClusterReader<N>::ClusterReader(bool useMC, bool usePatterns, bool triggerOut) : mUseMC(useMC), mUsePatterns(usePatterns), mTriggerOut(triggerOut), mDetName(Origin.as<std::string>()), mDetNameLC(mDetName)
{
  std::transform(mDetNameLC.begin(), mDetNameLC.end(), mDetNameLC.begin(), ::tolower);
}

template <int N>
void ClusterReader<N>::init(InitContext& ic)
{
  mFileName = o2::utils::Str::concat_string(o2::utils::Str::rectifyDirectory(ic.options().get<std::string>("input-dir")),
                                            ic.options().get<std::string>((mDetNameLC + "-cluster-infile").c_str()));
  connectTree(mFileName);
}

template <int N>
void ClusterReader<N>::run(ProcessingContext& pc)
{
  auto ent = mTree->GetReadEntry() + 1;
  assert(ent < mTree->GetEntries()); // this should not happen
  mTree->GetEntry(ent);

  for (uint32_t iLayer = 0; iLayer < NLayers; ++iLayer) {
    LOG(info) << mDetName << "ClusterReader:" << iLayer << " pushes " << mClusROFRec[iLayer]->size() << " ROFRecords, " << mClusterCompArray[iLayer]->size() << " compact clusters at entry " << ent;
    pc.outputs().snapshot(Output{Origin, "CLUSTERSROF", iLayer}, *mClusROFRec[iLayer]);
    pc.outputs().snapshot(Output{Origin, "COMPCLUSTERS", iLayer}, *mClusterCompArray[iLayer]);
    if (mUsePatterns) {
      pc.outputs().snapshot(Output{Origin, "PATTERNS", iLayer}, *mPatternsArray[iLayer]);
    }
    if (mUseMC) {
      pc.outputs().snapshot(Output{Origin, "CLUSTERSMCTR", iLayer}, *mClusterMCTruth[iLayer]);
      pc.outputs().snapshot(Output{Origin, "CLUSTERSMC2ROF", iLayer}, *mClusMC2ROFs[iLayer]);
    }
  }
  if (mTriggerOut) {
    std::vector<o2::itsmft::PhysTrigger> dummyTrig;
    pc.outputs().snapshot(Output{Origin, "PHYSTRIG", 0}, dummyTrig);
  }
  if (mTree->GetReadEntry() + 1 >= mTree->GetEntries()) {
    pc.services().get<ControlService>().endOfStream();
    pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
  }
}

template <int N>
void ClusterReader<N>::connectTree(const std::string& filename)
{
  mTree.reset(nullptr); // in case it was already loaded
  mFile.reset(TFile::Open(filename.c_str()));
  assert(mFile && !mFile->IsZombie());
  mTree.reset((TTree*)mFile->Get(mClusTreeName.c_str()));
  assert(mTree);

  for (uint32_t iLayer = 0; iLayer < NLayers; ++iLayer) {
    mTree->SetBranchAddress(getBranchName(mClusROFBranchName, iLayer).c_str(), &mClusROFRec[iLayer]);
    mTree->SetBranchAddress(getBranchName(mClusterCompBranchName, iLayer).c_str(), &mClusterCompArray[iLayer]);
    if (mUsePatterns) {
      mTree->SetBranchAddress(getBranchName(mClusterPattBranchName, iLayer).c_str(), &mPatternsArray[iLayer]);
    }
    if (mUseMC) {
      if (mTree->GetBranch(getBranchName(mClustMCTruthBranchName, iLayer).c_str()) &&
          mTree->GetBranch(getBranchName(mClustMC2ROFBranchName, iLayer).c_str())) {
        mTree->SetBranchAddress(getBranchName(mClustMCTruthBranchName, iLayer).c_str(), &mClusterMCTruth[iLayer]);
        mTree->SetBranchAddress(getBranchName(mClustMC2ROFBranchName, iLayer).c_str(), &mClusMC2ROFs[iLayer]);
      } else {
        LOG(info) << "MC-truth is missing";
        mUseMC = false;
      }
    }
  }
  LOG(info) << "Loaded tree from " << filename << " with " << mTree->GetEntries() << " entries";
}

template <int N>
std::string ClusterReader<N>::getBranchName(const std::string& base, int index) const
{
  return mDetName + base + "_" + std::to_string(index);
}

namespace
{
template <int N>
std::vector<OutputSpec> makeOutChannels(o2::header::DataOrigin detOrig, bool mctruth, bool usePatterns, bool triggerOut)
{
  std::vector<OutputSpec> outputs;
  for (int iLayer = 0; iLayer < o2::itsmft::DPLAlpideParam<N>::getNLayers(); ++iLayer) {
    outputs.emplace_back(detOrig, "CLUSTERSROF", iLayer, Lifetime::Timeframe);
    outputs.emplace_back(detOrig, "COMPCLUSTERS", iLayer, Lifetime::Timeframe);
    if (usePatterns) {
      outputs.emplace_back(detOrig, "PATTERNS", iLayer, Lifetime::Timeframe);
    }
    if (mctruth) {
      outputs.emplace_back(detOrig, "CLUSTERSMCTR", iLayer, Lifetime::Timeframe);
      outputs.emplace_back(detOrig, "CLUSTERSMC2ROF", iLayer, Lifetime::Timeframe);
    }
  }
  if (triggerOut) {
    outputs.emplace_back(detOrig, "PHYSTRIG", 0, Lifetime::Timeframe);
  }
  return outputs;
}
} // namespace

DataProcessorSpec getITSClusterReaderSpec(bool useMC, bool usePatterns, bool triggerOut)
{
  return DataProcessorSpec{
    .name = "its-cluster-reader",
    .inputs = Inputs{},
    .outputs = makeOutChannels<o2::detectors::DetID::ITS>("ITS", useMC, usePatterns, triggerOut),
    .algorithm = AlgorithmSpec{adaptFromTask<ITSClusterReader>(useMC, usePatterns, triggerOut)},
    .options = Options{
      {"its-cluster-infile", VariantType::String, "o2clus_its.root", {"Name of the input cluster file"}},
      {"input-dir", VariantType::String, "none", {"Input directory"}}}};
}

DataProcessorSpec getMFTClusterReaderSpec(bool useMC, bool usePatterns, bool triggerOut)
{
  return DataProcessorSpec{
    .name = "mft-cluster-reader",
    .inputs = Inputs{},
    .outputs = makeOutChannels<o2::detectors::DetID::MFT>("MFT", useMC, usePatterns, triggerOut),
    .algorithm = AlgorithmSpec{adaptFromTask<MFTClusterReader>(useMC, usePatterns, triggerOut)},
    .options = Options{
      {"mft-cluster-infile", VariantType::String, "mftclusters.root", {"Name of the input cluster file"}},
      {"input-dir", VariantType::String, "none", {"Input directory"}}}};
}

} // namespace itsmft
} // namespace o2
