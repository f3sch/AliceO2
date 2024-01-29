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

/// \file MatchTPCITS.h
/// \brief Class to perform TPC ITS matching
/// \author ruben.shahoyan@cern.ch

#ifndef ALICEO2_GLOBTRACKING_MATCHTPCITS_QC_
#define ALICEO2_GLOBTRACKING_MATCHTPCITS_QC_

#include <TH1F.h>
#include <TH2F.h>
#include <TEfficiency.h>
#include <TObjArray.h>
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "Framework/ProcessingContext.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "Steer/MCKinematicsReader.h"
#include "ReconstructionDataFormats/PID.h"
#include <unordered_map>
#include <vector>
#include <array>
#include <set>

namespace o2::globaltracking
{

using GID = o2::dataformats::GlobalTrackID;

struct LblInfo {
  int mIdx = -1;
  bool mIsPhysicalPrimary = false;
  bool mIsDuplicate = false;
};

class MatchITSTPCQC
{
 public:
  enum matchType : int { TPC = 0,
                         ITS,
                         SIZE };

  MatchITSTPCQC() = default;
  ~MatchITSTPCQC();

  bool init();
  void initDataRequest();
  void run(o2::framework::ProcessingContext& ctx);
  void setDataRequest(const std::shared_ptr<o2::globaltracking::DataRequest>& dr) { mDataRequest = dr; }
  void finalize();
  void reset();
  void getHistos(TObjArray& objar);

  /// \brief Publishes the histograms to the publisher e.g. the one provided by the QC task
  /// \tparam T type of the publisher
  /// \param publisher the publisher e.g. getObjectsManager()
  template <typename T>
  void publishHistograms(const std::shared_ptr<T>& publisher)
  {
    for (int i = 0; i < matchType::SIZE; ++i) {
      publisher->startPublishing(mPtNum[i]);
      publisher->startPublishing(mPtDen[i]);
      publisher->startPublishing(mFractionITSTPCmatch[i]);

      publisher->startPublishing(mPtNum_noEta0[i]);
      publisher->startPublishing(mPtDen_noEta0[i]);
      publisher->startPublishing(mFractionITSTPCmatch_noEta0[i]);

      publisher->startPublishing(mPhiNum[i]);
      publisher->startPublishing(mPhiDen[i]);
      publisher->startPublishing(mFractionITSTPCmatchPhi[i]);

      if (mUseTrkPID) { // Vs Tracking PID hypothesis
        for (int j = 0; j < o2::track::PID::NIDs; ++j) {
          // Pt
          publisher->startPublishing(mPtNumVsTrkPID[i][j]);
          publisher->startPublishing(mPtDenVsTrkPID[i][j]);
          publisher->startPublishing(mFractionITSTPCmatchPtVsTrkPID[i][j]);
          // Phi
          publisher->startPublishing(mPhiNumVsTrkPID[i][j]);
          publisher->startPublishing(mPhiDenVsTrkPID[i][j]);
          publisher->startPublishing(mFractionITSTPCmatchPhiVsTrkPID[i][j]);
          // Eta
          publisher->startPublishing(mEtaNumVsTrkPID[i][j]);
          publisher->startPublishing(mEtaDenVsTrkPID[i][j]);
          publisher->startPublishing(mFractionITSTPCmatchEtaVsTrkPID[i][j]);
        }
      }

      publisher->startPublishing(mPhiVsPtNum[i]);
      publisher->startPublishing(mPhiVsPtDen[i]);
      publisher->startPublishing(mFractionITSTPCmatchPhiVsPt[i]);

      publisher->startPublishing(mEtaNum[i]);
      publisher->startPublishing(mEtaDen[i]);
      publisher->startPublishing(mFractionITSTPCmatchEta[i]);

      publisher->startPublishing(mEtaVsPtNum[i]);
      publisher->startPublishing(mEtaVsPtDen[i]);
      publisher->startPublishing(mFractionITSTPCmatchEtaVsPt[i]);

      publisher->startPublishing(mClsVsPtNum[i]);
      publisher->startPublishing(mClsVsPtDen[i]);
      publisher->startPublishing(mFractionITSTPCmatchClsVsPt[i]);

      publisher->startPublishing(mChi2VsPtNum[i]);
      publisher->startPublishing(mChi2VsPtDen[i]);
      publisher->startPublishing(mFractionITSTPCmatchChi2VsPt[i]);

      publisher->startPublishing(m1OverPtNum[i]);
      publisher->startPublishing(m1OverPtDen[i]);
      publisher->startPublishing(mFractionITSTPCmatch1OverPt[i]);

      // MC
      if (mUseMC) {
        publisher->startPublishing(mPtPhysPrimNum[i]);
        publisher->startPublishing(mPtPhysPrimDen[i]);
        publisher->startPublishing(mFractionITSTPCmatchPhysPrim[i]);

        publisher->startPublishing(m1OverPtPhysPrimNum[i]);
        publisher->startPublishing(m1OverPtPhysPrimDen[i]);
        publisher->startPublishing(mFractionITSTPCmatchPhysPrim1OverPt[i]);

        publisher->startPublishing(mPhiPhysPrimNum[i]);
        publisher->startPublishing(mPhiPhysPrimDen[i]);
        publisher->startPublishing(mFractionITSTPCmatchPhiPhysPrim[i]);

        publisher->startPublishing(mEtaPhysPrimNum[i]);
        publisher->startPublishing(mEtaPhysPrimDen[i]);
        publisher->startPublishing(mFractionITSTPCmatchEtaPhysPrim[i]);

        publisher->startPublishing(mDupPtNum[i]);
        publisher->startPublishing(mDupPtDen[i]);
        publisher->startPublishing(mFractionITSTPCmatchDupPt[i]);

        publisher->startPublishing(mQPtNum[i]);
        publisher->startPublishing(mQPtDen[i]);
        publisher->startPublishing(mFractionITSTPCmatchQPt[i]);
      }
    }
    publisher->startPublishing(mChi2Matching);
    publisher->startPublishing(mChi2Refit);
    publisher->startPublishing(mTimeResVsPt);
    publisher->startPublishing(mResidualPt);
    publisher->startPublishing(mResidualPhi);
    publisher->startPublishing(mResidualEta);
    publisher->startPublishing(mDCAr);
    publisher->startPublishing(mDCArVsPtNum);
    publisher->startPublishing(mDCArVsPtDen);
    publisher->startPublishing(mFractionITSTPCmatchDCArVsPt);
  }

  void setSources(GID::mask_t src) { mSrc = src; }
  void setUseTrkPID(bool b) { mUseTrkPID = b; }
  bool getUseTrkPID() const { return mUseTrkPID; }
  void setUseMC(bool b) { mUseMC = b; }
  bool getUseMC() const { return mUseMC; }
  void setDebugHist(bool b) { mDebugHist = b; }
  bool getDebugHist() const { return mDebugHist; }
  void deleteHistograms();
  void setBz(float bz) { mBz = bz; }

  // ITS track
  void setMinPtITSCut(float v) { mPtITSCut = v; };
  void setEtaITSCut(float v) { mEtaITSCut = v; }; // TODO: define 2 different values for min and max (**)
  void setMinNClustersITS(float v) { mMinNClustersITS = v; }
  void setMaxChi2PerClusterITS(float v) { mMaxChi2PerClusterITS = v; }
  // TO DO: define an agreed way to implement the setter for ITS matching (min. # layers, which layers)
  // [...] --> exploit the method TrackCuts::setRequireHitsInITSLayers(...)
  // TPC track
  void setMinPtTPCCut(float v) { mPtTPCCut = v; };
  void setEtaTPCCut(float v) { mEtaTPCCut = v; }; // TODO: define 2 different values for min and max (***)
  void setMinNTPCClustersCut(float v) { mNTPCClustersCut = v; }
  void setMinDCAtoBeamPipeCut(std::array<float, 2> v)
  {
    setMinDCAtoBeamPipeDistanceCut(v[0]);
    setMinDCAtoBeamPipeYCut(v[1]);
  }
  void setMinDCAtoBeamPipeDistanceCut(float v) { mDCATPCCut = v; }
  void setMinDCAtoBeamPipeYCut(float v) { mDCATPCCutY = v; }
  // ITS-TPC kinematics
  void setPtCut(float v) { mPtCut = v; }
  void setMaxPtCut(float v) { mPtMaxCut = v; }
  void setEtaCut(float v) { mEtaCut = v; } // TODO: define 2 different values for min and max (*)

 private:
  std::shared_ptr<o2::globaltracking::DataRequest> mDataRequest;
  o2::globaltracking::RecoContainer mRecoCont;
  GID::mask_t mSrc = GID::getSourcesMask("ITS,TPC,ITS-TPC");
  GID::mask_t mAllowedSources = GID::getSourcesMask("ITS,TPC,ITS-TPC");
  // TPC
  gsl::span<const o2::tpc::TrackTPC> mTPCTracks;
  // ITS
  gsl::span<const o2::its::TrackITS> mITSTracks;
  // ITS-TPC
  gsl::span<const o2::dataformats::TrackTPCITS> mITSTPCTracks;
  bool mUseMC = false;                                                                     // Usage of the MC information
  bool mDebugHist = false;                                                                 // Debug Histograms
  bool mUseTrkPID = false;                                                                 // Usage of the PID hypothesis in tracking
  float mBz = 0;                                                                           ///< nominal Bz
  std::array<std::unordered_map<o2::MCCompLabel, LblInfo>, matchType::SIZE> mMapLabels;    // map with labels that have been found for the matched ITSTPC tracks; key is the label,
                                                                                           // value is the LbLinfo with the id of the track with the highest pT found with that label so far,
                                                                                           // and the flag to say if it is a physical primary or not
  std::array<std::unordered_map<o2::MCCompLabel, LblInfo>, matchType::SIZE> mMapRefLabels; // map with labels that have been found for the unmatched TPC tracks; key is the label,
                                                                                           // value is the LblInfo with the id of the track with the highest number of TPC clusters found
                                                                                           // with that label so far, and the flag to say if it is a physical primary or not
  o2::steer::MCKinematicsReader mcReader;                                                  // reader of MC information
  TDatabasePDG* mPDGDB;                                                                    // PDG Database

  // Pt
  TH1F* mPtNum[matchType::SIZE] = {};
  TH1F* mPtDen[matchType::SIZE] = {};
  TEfficiency* mFractionITSTPCmatch[matchType::SIZE] = {};
  TH1F* mPtNum_noEta0[matchType::SIZE] = {};
  TH1F* mPtDen_noEta0[matchType::SIZE] = {};
  TEfficiency* mFractionITSTPCmatch_noEta0[matchType::SIZE] = {};
  // Pt split per PID hypothesis in tracking
  TH1F* mPtNumVsTrkPID[matchType::SIZE][track::PID::NIDs] = {};
  TH1F* mPtDenVsTrkPID[matchType::SIZE][track::PID::NIDs] = {};
  TEfficiency* mFractionITSTPCmatchPtVsTrkPID[matchType::SIZE][track::PID::NIDs] = {};
  // Phi
  TH1F* mPhiNum[matchType::SIZE] = {};
  TH1F* mPhiDen[matchType::SIZE] = {};
  TEfficiency* mFractionITSTPCmatchPhi[matchType::SIZE] = {};
  TH2F* mPhiVsPtNum[matchType::SIZE] = {};
  TH2F* mPhiVsPtDen[matchType::SIZE] = {};
  TEfficiency* mFractionITSTPCmatchPhiVsPt[matchType::SIZE] = {};
  // Phi split per PID hypothesis in tracking
  TH1F* mPhiNumVsTrkPID[matchType::SIZE][track::PID::NIDs] = {};
  TH1F* mPhiDenVsTrkPID[matchType::SIZE][track::PID::NIDs] = {};
  TEfficiency* mFractionITSTPCmatchPhiVsTrkPID[matchType::SIZE][track::PID::NIDs] = {};
  // Eta
  TH1F* mEtaNum[matchType::SIZE] = {};
  TH1F* mEtaDen[matchType::SIZE] = {};
  TEfficiency* mFractionITSTPCmatchEta[matchType::SIZE] = {};
  TH2F* mEtaVsPtNum[matchType::SIZE] = {};
  TH2F* mEtaVsPtDen[matchType::SIZE] = {};
  TEfficiency* mFractionITSTPCmatchEtaVsPt[matchType::SIZE] = {};
  // Clusters
  TH2F* mClsVsPtNum[matchType::SIZE] = {};
  TH2F* mClsVsPtDen[matchType::SIZE] = {};
  TEfficiency* mFractionITSTPCmatchClsVsPt[matchType::SIZE] = {};
  // Chi2
  TH2F* mChi2VsPtNum[matchType::SIZE] = {};
  TH2F* mChi2VsPtDen[matchType::SIZE] = {};
  TEfficiency* mFractionITSTPCmatchChi2VsPt[matchType::SIZE] = {};
  // Eta split per PID hypothesis in tracking
  TH1F* mEtaNumVsTrkPID[matchType::SIZE][track::PID::NIDs] = {};
  TH1F* mEtaDenVsTrkPID[matchType::SIZE][track::PID::NIDs] = {};
  TEfficiency* mFractionITSTPCmatchEtaVsTrkPID[matchType::SIZE][track::PID::NIDs] = {};
  // Residuals
  TH2F* mResidualPt = nullptr;
  TH2F* mResidualPhi = nullptr;
  TH2F* mResidualEta = nullptr;
  // Others
  TH1F* mChi2Matching = nullptr;
  TH1F* mChi2Refit = nullptr;
  TH2F* mTimeResVsPt = nullptr;
  TH1F* mDCAr = nullptr;
  TH2F* mDCArVsPtNum = nullptr;
  TH2F* mDCArVsPtDen = nullptr;
  TEfficiency* mFractionITSTPCmatchDCArVsPt = nullptr;
  // 1/Pt
  TH1F* m1OverPtNum[matchType::SIZE] = {};
  TH1F* m1OverPtDen[matchType::SIZE] = {};
  TEfficiency* mFractionITSTPCmatch1OverPt[matchType::SIZE] = {};

  /// MC
  // Pt
  TH1F* mPtPhysPrimNum[matchType::SIZE] = {};
  TH1F* mPtPhysPrimDen[matchType::SIZE] = {};
  TEfficiency* mFractionITSTPCmatchPhysPrim[matchType::SIZE] = {};
  // 1/Pt
  TH1F* m1OverPtPhysPrimNum[matchType::SIZE] = {};
  TH1F* m1OverPtPhysPrimDen[matchType::SIZE] = {};
  TEfficiency* mFractionITSTPCmatchPhysPrim1OverPt[matchType::SIZE] = {};
  // Phi
  TH1F* mPhiPhysPrimNum[matchType::SIZE] = {};
  TH1F* mPhiPhysPrimDen[matchType::SIZE] = {};
  TEfficiency* mFractionITSTPCmatchPhiPhysPrim[matchType::SIZE] = {};
  // Eta
  TH1F* mEtaPhysPrimNum[matchType::SIZE] = {};
  TH1F* mEtaPhysPrimDen[matchType::SIZE] = {};
  TEfficiency* mFractionITSTPCmatchEtaPhysPrim[matchType::SIZE] = {};
  // MC Duplicates
  TH1F* mDupPtNum[matchType::SIZE] = {};
  TH1F* mDupPtDen[matchType::SIZE] = {};
  TEfficiency* mFractionITSTPCmatchDupPt[matchType::SIZE] = {};
  // MC Wrong Charge
  TH1F* mQPtNum[matchType::SIZE] = {};
  TH1F* mQPtDen[matchType::SIZE] = {};
  TEfficiency* mFractionITSTPCmatchQPt[matchType::SIZE] = {};

  void setEfficiency(TEfficiency* eff, TH1* hnum, TH1* hden, bool is2D = false);

  int mNTPCSelectedTracks = 0;
  int mNITSSelectedTracks = 0;
  int mNITSTPCSelectedTracks[matchType::SIZE] = {0, 0};

  // cut values
  // ITS track
  float mPtITSCut = 0.f;                                                // min pT for ITS track
  float mEtaITSCut = 1e10f;                                             // eta window for ITS track --> TODO: define 2 different values for min and max (**)
  int mMinNClustersITS = 0;                                             // min number of ITS clusters
  float mMaxChi2PerClusterITS{1e10f};                                   // max its fit chi2 per ITS cluster
  std::vector<std::pair<int8_t, std::set<uint8_t>>> mRequiredITSHits{}; // vector of ITS requirements (minNRequiredHits in specific requiredLayers)
  // TPC track
  float mPtTPCCut = 0.1f;        // min pT for TPC track
  float mEtaTPCCut = 1.4f;       // eta window for TPC track --> TODO: define 2 different values for min and max (***)
  int32_t mNTPCClustersCut = 60; // minimum number of TPC clusters for TPC track
  float mDCATPCCut = 100.f;      // max DCA 3D to PV for TPC track
  float mDCATPCCutY = 10.f;      // max DCA xy to PV for TPC track
  // ITS-TPC kinematics
  float mPtCut = 0.1f;
  float mPtMaxCut = 1e10f;
  float mEtaCut = 1e10f; // 1e10f as defaults of Detectors/GlobalTracking/include/GlobalTracking/TrackCuts.h
                         // TODO: define 2 different values for min and max (*)

  ClassDefNV(MatchITSTPCQC, 2);
};
} // namespace o2::globaltracking

#endif
