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

#include "CommonUtils/TreeStreamRedirector.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/Task.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITSStudies/Helpers.h"
#include "ITSStudies/TrackExtension.h"
#include "SimulationDataFormat/MCTrack.h"
#include "Steer/MCKinematicsReader.h"

#include <bitset>

#include "TFile.h"

namespace o2::its::study
{
using namespace o2::framework;
using namespace o2::globaltracking;

using GTrackID = o2::dataformats::GlobalTrackID;
using o2::steer::MCKinematicsReader;
class TrackExtensionStudy : public Task
{
  struct ParticleInfo {
    int event;
    int pdg;
    float pt;
    float eta;
    float phi;
    int mother;
    int first;
    float vx;
    float vy;
    float vz;
    uint8_t clusters = 0u;
    uint8_t fakeClusters = 0u;
    uint8_t isReco = 0u;
    uint8_t isFake = 0u;
    bool isPrimary = false;
    unsigned char storedStatus = 2; /// not stored = 2, fake = 1, good = 0
    int prodProcess;
    o2::its::TrackITS track;
    MCTrack mcTrack;
  };

 public:
  TrackExtensionStudy(std::shared_ptr<DataRequest> dr,
                      mask_t src,
                      std::shared_ptr<o2::steer::MCKinematicsReader> kineReader,
                      std::shared_ptr<o2::base::GRPGeomRequest> gr) : mDataRequest(dr), mTracksSrc(src), mKineReader(kineReader), mGGCCDBRequest(gr)
  {
    LOGP(info, "Read MCKine reader with {} sources", mKineReader->getNSources());
  }

  ~TrackExtensionStudy() final = default;
  void init(InitContext& /*ic*/) final;
  void run(ProcessingContext& /*pc*/) final;
  void endOfStream(EndOfStreamContext& /*ec*/) final;
  void process();

 private:
  static constexpr std::array<uint8_t, 9> mBitPatternsBefore{15, 30, 31, 60, 62, 63, 120, 124, 126};
  static constexpr std::array<uint8_t, 16> mBitPatternsAfter{31, 47, 61, 62, 63, 79, 94, 95, 111, 121, 122, 123, 124, 125, 126, 127};

  void updateTimeDependentParams(ProcessingContext& pc);
  std::string mOutFileName = "TrackExtensionStudy.root";
  std::shared_ptr<MCKinematicsReader> mKineReader;
  GeometryTGeo* mGeometry{};

  gsl::span<const o2::itsmft::ROFRecord> mTracksROFRecords;
  gsl::span<const o2::its::TrackITS> mTracks;
  gsl::span<const o2::MCCompLabel> mTracksMCLabels;
  gsl::span<const o2::itsmft::CompClusterExt> mClusters;
  gsl::span<const int> mInputITSidxs;
  const o2::dataformats::MCLabelContainer* mClustersMCLCont{};

  GTrackID::mask_t mTracksSrc{};
  std::shared_ptr<DataRequest> mDataRequest;
  std::vector<std::vector<std::vector<ParticleInfo>>> mParticleInfo; // src/event/track
  unsigned short mMask = 0x7f;

  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  std::unique_ptr<utils::TreeStreamRedirector> mStream;
  bool mWithTree{false};

  std::unique_ptr<TH1D> mHTrackCounts;
  std::unique_ptr<TH1D> mHLengthAny, mHLengthGood, mHLengthFake;
  std::unique_ptr<TH1D> mHChi2Any, mHChi2Good, mHChi2Fake;
  std::unique_ptr<TH1D> mHPtAny, mHPtGood, mHPtFake;
  std::unique_ptr<TH1D> mHExtensionAny, mHExtensionGood, mHExtensionFake;
  std::unique_ptr<TH2D> mHExtensionPatternsAny, mHExtensionPatternsGood, mHExtensionPatternsFake, mHExtensionPatternsGoodMissed, mHExtensionPatternsGoodEmpty;
  ;

  template <class T, typename... C, typename... F>
  std::unique_ptr<T> createHistogram(C... n, F... b)
  {
    auto t = std::make_unique<T>(n..., b...);
    mHistograms.push_back(static_cast<TH1*>(t.get()));
    return std::move(t);
  }

  std::vector<TH1*> mHistograms;
};

void TrackExtensionStudy::init(InitContext& ic)
{
  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);
  mWithTree = ic.options().get<bool>("with-tree");

  constexpr size_t effHistBins = 100;
  constexpr float effPtCutLow = 0.01;
  constexpr float effPtCutHigh = 10.;
  auto xbins = helpers::makeLogBinning(effHistBins, effPtCutLow, effPtCutHigh);

  // Track Counting
  mHTrackCounts = createHistogram<TH1D>("hTrackCounts", "Track Stats", 10, 0, 10);
  mHTrackCounts->GetXaxis()->SetBinLabel(1, "Total Tracks");
  mHTrackCounts->GetXaxis()->SetBinLabel(2, "Normal ANY Tracks");
  mHTrackCounts->GetXaxis()->SetBinLabel(3, "Normal GOOD Tracks");
  mHTrackCounts->GetXaxis()->SetBinLabel(4, "Normal FAKE Tracks");
  mHTrackCounts->GetXaxis()->SetBinLabel(5, "Extended ANY Tracks");
  mHTrackCounts->GetXaxis()->SetBinLabel(6, "Extended GOOD Tracks");
  mHTrackCounts->GetXaxis()->SetBinLabel(7, "Extended FAKE Tracks");
  mHTrackCounts->GetXaxis()->SetBinLabel(8, "Extended FAKE BEFORE Tracks");
  mHTrackCounts->GetXaxis()->SetBinLabel(9, "Extended FAKE AFTER Tracks");
  mHTrackCounts->GetXaxis()->SetBinLabel(10, "Extended FAKE BEFORE&AFTER Tracks");

  // Length
  mHLengthAny = createHistogram<TH1D>("hLengthAny", "Extended Tracks Length (ANY);NCluster;Entries", 5, 3, 8);
  mHLengthGood = createHistogram<TH1D>("hLengthGood", "Extended Tracks Length (GOOD);NCluster;Entries", 5, 3, 8);
  mHLengthFake = createHistogram<TH1D>("hLengthFake", "Extended Tracks Length (FAKE);NCluster;Entries", 5, 3, 8);

  // Chi2
  mHChi2Any = createHistogram<TH1D>("hChi2Any", "Extended Tracks Length (ANY);#chi^{2};Entries", 50, 0, 100);
  mHChi2Good = createHistogram<TH1D>("hChi2Good", "Extended Tracks Length (GOOD);#chi^{2};Entries", 50, 0, 100);
  mHChi2Fake = createHistogram<TH1D>("hChi2Fake", "Extended Tracks Length (FAKE);#chi^{2};Entries", 50, 0, 100);

  // Pt
  mHPtAny = createHistogram<TH1D>("hPtAny", "Extended Tracks Length (ANY);#it{p}_{T};Entries", xbins.size(), effPtCutLow, effPtCutHigh);
  mHPtGood = createHistogram<TH1D>("hPtGood", "Extended Tracks Length (GOOD);#it{p}_{T};Entries", xbins.size(), effPtCutLow, effPtCutHigh);
  mHPtFake = createHistogram<TH1D>("hPtFake", "Extended Tracks Length (FAKE);#it{p}_{T};Entries", xbins.size(), effPtCutLow, effPtCutHigh);

  // Length
  mHExtensionAny = createHistogram<TH1D>("hExtensionAny", "Extended Tracks Length (ANY);Extended Layer;Entries", 7, 0, 7);
  mHExtensionGood = createHistogram<TH1D>("hExtensionGood", "Extended Tracks Length (GOOD);Extended Layer;Entries", 7, 0, 7);
  mHExtensionFake = createHistogram<TH1D>("hExtensionFake", "Extended Tracks Length (FAKE);Extended Layer;Entries", 7, 0, 7);

  // Patterns
  auto makePatternAxisLabels = [&](TH1* h, bool xBefore = true) {
    for (int i{1}; i <= h->GetXaxis()->GetNbins(); ++i) {
      if (xBefore) {
        h->GetXaxis()->SetBinLabel(i, fmt::format("{:07b}", mBitPatternsBefore[i - 1]).c_str());
      } else {
        h->GetXaxis()->SetBinLabel(i, fmt::format("{:07b}", mBitPatternsAfter[i - 1]).c_str());
      }
    }
    for (int i{1}; i <= h->GetYaxis()->GetNbins(); ++i) {
      h->GetYaxis()->SetBinLabel(i, fmt::format("{:07b}", mBitPatternsAfter[i - 1]).c_str());
    }
  };
  mHExtensionPatternsAny = createHistogram<TH2D>("hExtensionPatternsAny", "Extended Tracks Pattern (ANY);Before;After;Entries", mBitPatternsBefore.size(), 0, mBitPatternsBefore.size(), mBitPatternsAfter.size(), 0, mBitPatternsAfter.size());
  makePatternAxisLabels(mHExtensionPatternsAny.get());
  mHExtensionPatternsGood = createHistogram<TH2D>("hExtensionPatternsGood", "Extended Tracks Pattern (GOOD);Before;After;Entries", mBitPatternsBefore.size(), 0, mBitPatternsBefore.size(), mBitPatternsAfter.size(), 0, mBitPatternsAfter.size());
  makePatternAxisLabels(mHExtensionPatternsGood.get());
  mHExtensionPatternsFake = createHistogram<TH2D>("hExtensionPatternsFake", "Extended Tracks Pattern (FAKE);Before;After;Entries", mBitPatternsBefore.size(), 0, mBitPatternsBefore.size(), mBitPatternsAfter.size(), 0, mBitPatternsAfter.size());
  makePatternAxisLabels(mHExtensionPatternsFake.get());
  mHExtensionPatternsGoodMissed = createHistogram<TH2D>("hExtensionPatternsGoodMissed", "Extended Tracks Pattern (GOOD) Missed Clusters;After;Missed;Entries", mBitPatternsAfter.size(), 0, mBitPatternsAfter.size(), mBitPatternsAfter.size(), 0, mBitPatternsAfter.size());
  makePatternAxisLabels(mHExtensionPatternsGoodMissed.get(), false);
  mHExtensionPatternsGoodEmpty = createHistogram<TH2D>("hExtensionPatternsGoodEmpty", "Extended Tracks Pattern (GOOD) Empty Clusters;Before;After;Entries", mBitPatternsAfter.size(), 0, mBitPatternsAfter.size(), mBitPatternsAfter.size(), 0, mBitPatternsAfter.size());
  makePatternAxisLabels(mHExtensionPatternsGoodEmpty.get(), false);

  mStream = std::make_unique<utils::TreeStreamRedirector>(mOutFileName.c_str(), "RECREATE");
}

void TrackExtensionStudy::run(ProcessingContext& pc)
{
  o2::globaltracking::RecoContainer recoData;
  recoData.collectData(pc, *mDataRequest);
  updateTimeDependentParams(pc);

  mTracksROFRecords = recoData.getITSTracksROFRecords();
  mTracks = recoData.getITSTracks();
  mTracksMCLabels = recoData.getITSTracksMCLabels();
  mClusters = recoData.getITSClusters();
  mClustersMCLCont = recoData.getITSClustersMCLabels();
  mInputITSidxs = recoData.getITSTracksClusterRefs();

  LOGP(info, "** Found in {} rofs:\n\t- {} clusters with {} labels\n\t- {} tracks with {} labels",
       mTracksROFRecords.size(), mClusters.size(), mClustersMCLCont->getIndexedSize(), mTracks.size(), mTracksMCLabels.size());
  LOGP(info, "** Found {} sources from kinematic files", mKineReader->getNSources());

  process();
}

void TrackExtensionStudy::process()
{
  LOGP(info, "** Filling particle table ... ");
  mParticleInfo.resize(mKineReader->getNSources()); // sources
  for (int iSource{0}; iSource < mKineReader->getNSources(); ++iSource) {
    mParticleInfo[iSource].resize(mKineReader->getNEvents(iSource)); // events
    for (int iEvent{0}; iEvent < mKineReader->getNEvents(iSource); ++iEvent) {
      mParticleInfo[iSource][iEvent].resize(mKineReader->getTracks(iSource, iEvent).size()); // tracks
      for (auto iPart{0}; iPart < mKineReader->getTracks(iEvent).size(); ++iPart) {
        auto& part = mKineReader->getTracks(iSource, iEvent)[iPart];
        mParticleInfo[iSource][iEvent][iPart].event = iEvent;
        mParticleInfo[iSource][iEvent][iPart].pdg = part.GetPdgCode();
        mParticleInfo[iSource][iEvent][iPart].pt = part.GetPt();
        mParticleInfo[iSource][iEvent][iPart].phi = part.GetPhi();
        mParticleInfo[iSource][iEvent][iPart].eta = part.GetEta();
        mParticleInfo[iSource][iEvent][iPart].vx = part.Vx();
        mParticleInfo[iSource][iEvent][iPart].vy = part.Vy();
        mParticleInfo[iSource][iEvent][iPart].vz = part.Vz();
        mParticleInfo[iSource][iEvent][iPart].isPrimary = part.isPrimary();
        mParticleInfo[iSource][iEvent][iPart].mother = part.getMotherTrackId();
        mParticleInfo[iSource][iEvent][iPart].prodProcess = part.getProcess();
      }
    }
  }
  LOGP(info, "** Creating particle/clusters correspondance ... ");
  for (auto iSource{0}; iSource < mParticleInfo.size(); ++iSource) {
    for (auto iCluster{0}; iCluster < mClusters.size(); ++iCluster) {
      auto labs = mClustersMCLCont->getLabels(iCluster); // ideally I can have more than one label per cluster
      for (auto& lab : labs) {
        if (!lab.isValid()) {
          continue; // We want to skip channels related to noise, e.g. sID = 99: QED
        }
        int trackID, evID, srcID;
        bool fake;
        lab.get(trackID, evID, srcID, fake);
        auto& cluster = mClusters[iCluster];
        auto layer = mGeometry->getLayer(cluster.getSensorID());
        mParticleInfo[srcID][evID][trackID].clusters |= (1 << layer);
        if (fake) {
          mParticleInfo[srcID][evID][trackID].fakeClusters |= (1 << layer);
        }
      }
    }
  }

  LOGP(info, "** Analysing tracks ... ");
  int unaccounted{0}, good{0}, fakes{0}, extended{0};
  for (auto iTrack{0}; iTrack < mTracks.size(); ++iTrack) {
    const auto& lab = mTracksMCLabels[iTrack];
    if (!lab.isValid()) {
      unaccounted++;
      continue;
    }
    int trackID, evID, srcID;
    bool fake;
    lab.get(trackID, evID, srcID, fake);

    if (srcID == 99) { // skip QED
      unaccounted++;
      continue;
    }

    for (int iLayer{0}; iLayer < 7; ++iLayer) {
      if (mTracks[iTrack].isExtendedOnLayer(iLayer)) {
        ++extended;
        break;
      }
    }

    mParticleInfo[srcID][evID][trackID].isReco += !fake;
    mParticleInfo[srcID][evID][trackID].isFake += fake;
    if (mTracks[iTrack].isBetter(mParticleInfo[srcID][evID][trackID].track, 1.e9)) {
      mParticleInfo[srcID][evID][trackID].storedStatus = fake;
      mParticleInfo[srcID][evID][trackID].track = mTracks[iTrack];
      mParticleInfo[srcID][evID][trackID].mcTrack = *mKineReader->getTrack(lab);
    }
    fakes += fake;
    good += !fake;
  }
  LOGP(info, "** Some statistics:");
  LOGP(info, "\t- Total number of tracks: {}", mTracks.size());
  LOGP(info, "\t- Total number of tracks not corresponding to particles: {} ({:.2f} %)", unaccounted, unaccounted * 100. / mTracks.size());
  LOGP(info, "\t- Total number of fakes: {} ({:.2f} %)", fakes, fakes * 100. / mTracks.size());
  LOGP(info, "\t- Total number of good: {} ({:.2f} %)", good, good * 100. / mTracks.size());
  LOGP(info, "\t- Total number of extensions: {} ({:.2f} %)", extended, extended * 100. / mTracks.size());

  LOGP(info, "** Filling histograms ... ");
  for (auto iTrack{0}; iTrack < mTracks.size(); ++iTrack) {
    auto& lab = mTracksMCLabels[iTrack];
    if (!lab.isValid()) {
      unaccounted++;
      continue;
    }
    int trackID, evID, srcID;
    bool fake;
    lab.get(trackID, evID, srcID, fake);
    const auto& part = mParticleInfo[srcID][evID][trackID];
    if (!part.isPrimary) {
      continue;
    }
    const auto& trk = part.track;
    bool isGood = part.isReco && !part.isFake;
    mHTrackCounts->Fill(0);

    std::bitset<7> extPattern{0};
    for (int iLayer{0}; iLayer < 7; ++iLayer) {
      if (trk.isExtendedOnLayer(iLayer)) {
        extPattern.set(iLayer);
      }
    }

    if (!extPattern.any()) {
      mHTrackCounts->Fill(1);
      if (part.isReco || !part.isFake) {
        mHTrackCounts->Fill(2);
      } else {
        mHTrackCounts->Fill(3);
      }
      continue;
    }

    mHTrackCounts->Fill(4);
    mHLengthAny->Fill(trk.getNClusters());
    mHChi2Any->Fill(trk.getChi2());
    mHPtAny->Fill(trk.getPt());
    if (isGood) {
      mHTrackCounts->Fill(5);
      mHLengthGood->Fill(trk.getNClusters());
      mHChi2Good->Fill(trk.getChi2());
      mHPtGood->Fill(trk.getPt());
    } else {
      mHTrackCounts->Fill(6);
      mHLengthFake->Fill(trk.getNClusters());
      mHChi2Fake->Fill(trk.getChi2());
      mHPtFake->Fill(trk.getPt());
    }

    std::bitset<7> clusPattern{static_cast<uint8_t>(trk.getPattern())};
    for (int iLayer{0}; iLayer < 7; ++iLayer) {
      if (extPattern.test(iLayer)) {
        extPattern.set(iLayer);
        mHExtensionAny->Fill(iLayer);
        if (isGood) {
          mHExtensionGood->Fill(iLayer);
        } else {
          mHExtensionFake->Fill(iLayer);
        }
      }
    }
    std::bitset<7> oldPattern{clusPattern & ~extPattern}, holePattern{clusPattern};
    holePattern.flip();
    auto clusN = clusPattern.to_ulong();
    auto clusIdx = std::distance(std::begin(mBitPatternsAfter), std::find(std::begin(mBitPatternsAfter), std::end(mBitPatternsAfter), clusN));
    auto oldN = oldPattern.to_ulong();
    auto oldIdx = std::distance(std::begin(mBitPatternsBefore), std::find(std::begin(mBitPatternsBefore), std::end(mBitPatternsBefore), oldN));
    mHExtensionPatternsAny->Fill(oldIdx, clusIdx);
    if (isGood) {
      mHExtensionPatternsGood->Fill(oldIdx, clusIdx);
    } else {
      mHExtensionPatternsFake->Fill(oldIdx, clusIdx);
    }

    // old pattern
    bool oldFake{false}, newFake{false};
    for (int iLayer{0}; iLayer < 7; ++iLayer) {
      if (trk.isFakeOnLayer(iLayer)) {
        if (oldPattern.test(iLayer)) {
          oldFake = true;
        } else if (extPattern.test(iLayer)) {
          newFake = true;
        }
      }
    }
    if (oldFake && newFake) {
      mHTrackCounts->Fill(9);
    } else if (oldFake) {
      mHTrackCounts->Fill(7);
    } else if (newFake) {
      mHTrackCounts->Fill(8);
    }

    // Check if we missed some clusters
    if (isGood && holePattern.any()) {
      auto missPattern{clusPattern}, emptyPattern{clusPattern};
      for (int iLayer{0}; iLayer < 7; ++iLayer) {
        if (!holePattern.test(iLayer)) {
          continue;
        }

        // Check if there was actually a cluster that we missed
        if ((part.clusters & (1 << iLayer)) != 0) {
          missPattern.set(iLayer);
        } else {
          emptyPattern.set(iLayer);
        }
      }

      if (missPattern != clusPattern) {
        auto missN = missPattern.to_ulong();
        auto missIdx = std::distance(std::begin(mBitPatternsAfter), std::find(std::begin(mBitPatternsAfter), std::end(mBitPatternsAfter), missN));
        mHExtensionPatternsGoodMissed->Fill(clusIdx, missIdx);
      }
      if (emptyPattern != clusPattern) {
        auto emptyN = emptyPattern.to_ulong();
        auto emptyIdx = std::distance(std::begin(mBitPatternsAfter), std::find(std::begin(mBitPatternsAfter), std::end(mBitPatternsAfter), emptyN));
        mHExtensionPatternsGoodEmpty->Fill(clusIdx, emptyIdx);
      }
    }

    if (mWithTree) {
      std::array<float, 3> xyz{(float)part.mcTrack.GetStartVertexCoordinatesX(), (float)part.mcTrack.GetStartVertexCoordinatesY(), (float)part.mcTrack.GetStartVertexCoordinatesZ()};
      std::array<float, 3> pxyz{(float)part.mcTrack.GetStartVertexMomentumX(), (float)part.mcTrack.GetStartVertexMomentumY(), (float)part.mcTrack.GetStartVertexMomentumZ()};
      auto pdg = O2DatabasePDG::Instance()->GetParticle(part.pdg);
      if (pdg == nullptr) {
        LOGP(fatal, "MC info not available");
      }
      auto mcTrk = o2::track::TrackPar(xyz, pxyz, TMath::Nint(pdg->Charge() / 3.), true);
      (*mStream) << "tree"
                 << "trk=" << trk
                 << "mcTrk=" << mcTrk
                 << "\n";
    }
  }
}

void TrackExtensionStudy::updateTimeDependentParams(ProcessingContext& pc)
{
  static bool initOnceDone = false;
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  if (!initOnceDone) { // this params need to be queried only once
    initOnceDone = true;
    mGeometry = GeometryTGeo::Instance();
    mGeometry->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::T2GRot, o2::math_utils::TransformType::T2G));
  }
}

void TrackExtensionStudy::endOfStream(EndOfStreamContext& ec)
{
  LOGP(info, "Writing results to {}", mOutFileName);
  mStream->GetFile()->cd();
  for (const auto h : mHistograms) {
    h->Write();
  }
  mStream->Close();
}

DataProcessorSpec getTrackExtensionStudy(mask_t srcTracksMask, mask_t srcClustersMask, std::shared_ptr<o2::steer::MCKinematicsReader> kineReader)
{
  std::vector<OutputSpec> outputs;
  auto dataRequest = std::make_shared<DataRequest>();
  dataRequest->requestTracks(srcTracksMask, true);
  dataRequest->requestClusters(srcClustersMask, true);

  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
                                                              true,                              // GRPECS=true
                                                              false,                             // GRPLHCIF
                                                              true,                              // GRPMagField
                                                              true,                              // askMatLUT
                                                              o2::base::GRPGeomRequest::Aligned, // geometry
                                                              dataRequest->inputs,
                                                              true);

  return DataProcessorSpec{
    "its-study-track-extension",
    dataRequest->inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<TrackExtensionStudy>(dataRequest, srcTracksMask, kineReader, ggRequest)},
    Options{{"with-tree", o2::framework::VariantType::Bool, false, {"Produce in addition a tree"}}}};
}

} // namespace o2::its::study
