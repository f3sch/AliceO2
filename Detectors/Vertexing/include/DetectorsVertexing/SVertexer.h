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

/// \file SVertexer.h
/// \brief Secondary vertex finder
/// \author ruben.shahoyan@cern.ch
#ifndef O2_S_VERTEXER_H
#define O2_S_VERTEXER_H

#include "gsl/span"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "ReconstructionDataFormats/Decay3Body.h"
#include "ReconstructionDataFormats/DecayNBodyIndex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "ReconstructionDataFormats/VtxTrackRef.h"
#include "CommonDataFormat/RangeReference.h"
#include "DCAFitter/DCAFitterN.h"
#include "DetectorsVertexing/SVertexerParams.h"
#include "DetectorsVertexing/SVertexHypothesis.h"
#include "StrangenessTracking/StrangenessTracker.h"
#include "DataFormatsTPC/TrackTPC.h"
#include <numeric>
#include <algorithm>
#include <tuple>
#include "GPUO2InterfaceRefit.h"
#include "TPCFastTransform.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "SimulationDataFormat/MCUtils.h"
#include "Steer/MCKinematicsReader.h"
#include "boost/functional/hash.hpp"

namespace o2
{
namespace tpc
{
class VDriftCorrFact;
}
namespace gpu
{
class CorrectionMapsHelper;
}

namespace vertexing
{

namespace o2d = o2::dataformats;

class SVertexer
{
 public:
  using GIndex = o2::dataformats::VtxTrackIndex;
  using VRef = o2::dataformats::VtxTrackRef;
  using PVertex = const o2::dataformats::PrimaryVertex;
  using V0 = o2::dataformats::V0;
  using V0Index = o2::dataformats::V0Index;
  using Cascade = o2::dataformats::Cascade;
  using CascadeIndex = o2::dataformats::CascadeIndex;
  using Decay3Body = o2::dataformats::Decay3Body;
  using Decay3BodyIndex = o2::dataformats::Decay3BodyIndex;
  using RRef = o2::dataformats::RangeReference<int, int>;
  using VBracket = o2::math_utils::Bracket<int>;

  enum HypV0 { Photon,
               K0,
               Lambda,
               AntiLambda,
               HyperTriton,
               AntiHyperTriton,
               Hyperhydrog4,
               AntiHyperhydrog4,
               NHypV0 };

  enum HypCascade {
    XiMinus,
    OmegaMinus,
    NHypCascade
  };

  enum Hyp3body {
    H3L3body,
    AntiH3L3body,
    He4L3body,
    AntiHe4L3body,
    He5L3body,
    AntiHe5L3body,
    NHyp3body
  };

  static constexpr int POS = 0, NEG = 1;
  struct TrackCand : o2::track::TrackParCov {
    GIndex gid{};
    VBracket vBracket{};
    float minR = 0; // track lowest point r
  };

  SVertexer(bool enabCascades = true, bool enab3body = false) : mEnableCascades{enabCascades}, mEnable3BodyDecays{enab3body}
  {
  }

  void setEnableCascades(bool v) { mEnableCascades = v; }
  void setEnable3BodyDecays(bool v) { mEnable3BodyDecays = v; }
  void init();
  void process(const o2::globaltracking::RecoContainer& recoTracks, o2::framework::ProcessingContext& pc);
  int getNV0s() const { return mNV0s; }
  int getNCascades() const { return mNCascades; }
  int getN3Bodies() const { return mN3Bodies; }
  int getNStrangeTracks() const { return mNStrangeTracks; }
  auto& getMeanVertex() const { return mMeanVertex; }
  void setMeanVertex(const o2d::VertexBase& v) { mMeanVertex = v; }
  void setNThreads(int n);
  int getNThreads() const { return mNThreads; }
  void setUseMC(bool v) { mUseMC = v; }
  bool getUseMC() const { return mUseMC; }
  void setUseDebug(bool v) { mUseDebug = v; }
  bool getUseDebug() const { return mUseDebug; }

  void setTPCTBin(int nbc)
  {
    // set TPC time bin in BCs
    mMUS2TPCBin = 1.f / (nbc * o2::constants::lhc::LHCBunchSpacingMUS);
  }
  void setTPCVDrift(const o2::tpc::VDriftCorrFact& v);
  void setTPCCorrMaps(o2::gpu::CorrectionMapsHelper* maph);
  void initTPCTransform();
  void setStrangenessTracker(o2::strangeness_tracking::StrangenessTracker* tracker) { mStrTracker = tracker; }
  o2::strangeness_tracking::StrangenessTracker* getStrangenessTracker() { return mStrTracker; }

 private:
  template <class TVI, class TCI, class T3I, class TR>
  void extractPVReferences(const TVI& v0s, TR& vtx2V0Refs, const TCI& cascades, TR& vtx2CascRefs, const T3I& vtxs3, TR& vtx2body3Refs);
  bool checkV0(const TrackCand& seed0, const TrackCand& seed1, int iP, int iN, int ithread);
  int checkCascades(const V0Index& v0Idx, const V0& v0, float rv0, std::array<float, 3> pV0, float p2V0, int avoidTrackID, int posneg, VBracket v0vlist, int ithread);
  int check3bodyDecays(const V0Index& v0Idx, const V0& v0, float rv0, std::array<float, 3> pV0, float p2V0, int avoidTrackID, int posneg, VBracket v0vlist, int ithread);
  void setupThreads();
  void buildT2V(const o2::globaltracking::RecoContainer& recoTracks);
  void updateTimeDependentParams();
  bool acceptTrack(GIndex gid, const o2::track::TrackParCov& trc) const;
  bool processTPCTrack(const o2::tpc::TrackTPC& trTPC, GIndex gid, int vtxid);
  float correctTPCTrack(o2::track::TrackParCov& trc, const o2::tpc::TrackTPC tTPC, float tmus, float tmusErr) const;

  uint64_t getPairIdx(GIndex id1, GIndex id2) const
  {
    return (uint64_t(id1) << 32) | id2;
  }

  // at the moment not used
  o2::gpu::CorrectionMapsHelper* mTPCCorrMapsHelper = nullptr;
  std::unique_ptr<o2::gpu::GPUO2InterfaceRefit> mTPCRefitter; ///< TPC refitter used for TPC tracks refit during the reconstruction
  o2::strangeness_tracking::StrangenessTracker* mStrTracker = nullptr;
  gsl::span<const PVertex> mPVertices;
  std::vector<std::vector<V0>> mV0sTmp;
  std::vector<std::vector<Cascade>> mCascadesTmp;
  std::vector<std::vector<Decay3Body>> m3bodyTmp;
  std::vector<std::vector<V0Index>> mV0sIdxTmp;
  std::vector<std::vector<CascadeIndex>> mCascadesIdxTmp;
  std::vector<std::vector<Decay3BodyIndex>> m3bodyIdxTmp;
  std::array<std::vector<TrackCand>, 2> mTracksPool{}; // pools of positive and negative seeds sorted in min VtxID
  std::array<std::vector<int>, 2> mVtxFirstTrack{};    // 1st pos. and neg. track of the pools for each vertex

  o2d::VertexBase mMeanVertex{{0., 0., 0.}, {0.1 * 0.1, 0., 0.1 * 0.1, 0., 0., 6. * 6.}};
  const SVertexerParams* mSVParams = nullptr;
  std::array<SVertexHypothesis, NHypV0> mV0Hyps;
  std::array<SVertexHypothesis, NHypCascade> mCascHyps;
  std::array<SVertex3Hypothesis, NHyp3body> m3bodyHyps;

  std::vector<DCAFitterN<2>> mFitterV0;
  std::vector<DCAFitterN<2>> mFitterCasc;
  std::vector<DCAFitterN<3>> mFitter3body;
  int mNThreads = 1;
  int mNV0s = 0, mNCascades = 0, mN3Bodies = 0, mNStrangeTracks = 0;
  float mMinR2ToMeanVertex = 0;
  float mMaxDCAXY2ToMeanVertex = 0;
  float mMaxDCAXY2ToMeanVertexV0Casc = 0;
  float mMaxDCAXY2ToMeanVertex3bodyV0 = 0;
  float mMinR2DiffV0Casc = 0;
  float mMaxR2ToMeanVertexCascV0 = 0;
  float mMinPt2V0 = 1e-6;
  float mMaxTgl2V0 = 2. * 2.;
  float mMinPt2Casc = 1e-4;
  float mMaxTgl2Casc = 2. * 2.;
  float mMinPt23Body = 1e-4;
  float mMaxTgl23Body = 2.f * 2.f;
  float mMUS2TPCBin = 1.f / (8 * o2::constants::lhc::LHCBunchSpacingMUS);
  float mTPCBin2Z = 0;
  float mTPCVDrift = 0;
  float mTPCVDriftCorrFact = 1.; ///< TPC nominal correction factort (wrt ref)
  float mTPCVDriftRef = 0;
  float mTPCDriftTimeOffset = 0; ///< drift time offset in mus

  bool mEnableCascades = true;
  bool mEnable3BodyDecays = false;
  bool mUseMC = false;
  bool mUseDebug = false;

  long mCounterTPConly{0};
  long mCounterTrueGammas{0};
  long mCounterTrueGammasITS{0};
  long mCounterTrueGammasTPC{0};
  long mCounterTrueGammasITSTPC{0};
  std::vector<double> mTrueGammasITSPt;
  std::vector<double> mTrueGammasTPCPt;
  std::vector<double> mTrueGammasITSTPCPt;
  gsl::span<const o2::MCCompLabel> mITSTrkLabels;    ///< input ITS Track MC labels
  gsl::span<const o2::MCCompLabel> mTPCTrkLabels;    ///< input TPC Track MC labels
  gsl::span<const o2::MCCompLabel> mITSTPCTrkLabels; ///< input ITS-TPC Track MC labels
  std::unique_ptr<o2::utils::TreeStreamRedirector> mDebugStream;
  static constexpr std::string_view mDebugStreamName{"candTracks"};
  o2::steer::MCKinematicsReader mcReader; // reader of MC information
  void writeDebugV0Candidates(o2::tpc::TrackTPC const& trk, GIndex gid, int vtxid, o2::track::TrackParCov const& candTrk);
  void writeDebugBTrackPools(const o2::globaltracking::RecoContainer& recoData);
  void writeDebugATrackPools(const o2::globaltracking::RecoContainer& recoData);
  template <class TVI, class RECO>
  void writeDebugV0Found(TVI const& v0s, RECO const& recoData);

  using key_t = std::tuple<int, int, int>;

  struct key_hash {
    size_t operator()(const key_t& k) const
    {
      const auto& [eve, src, trkid] = k;
      size_t seed = 0;
      boost::hash_combine(seed, eve);
      boost::hash_combine(seed, src);
      boost::hash_combine(seed, trkid);
      return seed;
    }
  };
  using map_before_t = std::unordered_map<std::tuple<int, int, int>, std::tuple<MCTrack, MCTrack, MCTrack, bool>, key_hash>;
  using map_after_t = std::unordered_map<std::tuple<int, int, int>, std::tuple<TrackCand, MCTrack, TrackCand, MCTrack, MCTrack>, key_hash>;

  bool checkMother(o2::MCTrack const* mother, const std::vector<o2::MCTrack>& pcontainer)
  {
    if (mother == nullptr) {
      return false;
    }
    bool isPhysicalPrimary = o2::mcutils::MCTrackNavigator::isPhysicalPrimary(*mother, pcontainer);
    if (mother->GetPdgCode() != 22 || !mother->isPrimary() || !isPhysicalPrimary) {
      return false;
    }
    return true;
  }

  o2::MCCompLabel getLabel(GIndex const& gid, o2::globaltracking::RecoContainer const& recoData)
  {
    if (gid.getSource() == GIndex::ITS && recoData.isTrackSourceLoaded(GIndex::ITS)) {
      return mITSTrkLabels[gid.getIndex()];
    } else if (gid.getSource() == GIndex::ITSTPC && recoData.isTrackSourceLoaded(GIndex::ITSTPC)) {
      return mITSTPCTrkLabels[gid.getIndex()];
    } else if (gid.getSource() == GIndex::TPC && recoData.isTrackSourceLoaded(GIndex::TPC)) {
      return mTPCTrkLabels[gid.getIndex()];
    }
    return o2::MCCompLabel(true);
  }

  bool checkLabels(o2::MCCompLabel const& lbl0, o2::MCCompLabel const& lbl1)
  {
    if (!lbl0.isValid() || !lbl1.isValid() || lbl0.isFake() || lbl1.isFake()) {
      return false;
    }
    if (lbl0.getEventID() != lbl1.getEventID()) {
      return false;
    }
    if (lbl0.getSourceID() != lbl1.getSourceID()) {
      return false;
    }
    if (lbl0.getTrackID() == lbl1.getTrackID()) {
      return false; // possible looper?
    }
    return true;
  }

  bool checkPair(o2::MCTrack const* mcTrk0, o2::MCTrack const* mcTrk1)
  {
    if (mcTrk0 == nullptr || mcTrk1 == nullptr || mcTrk0 == mcTrk1) {
      return false;
    }
    if (auto mcMotherId0 = mcTrk0->getMotherTrackId(), mcMotherId1 = mcTrk1->getMotherTrackId();
        mcMotherId0 != mcMotherId1 || mcMotherId0 == -1 || mcMotherId1 == -1) {
      return false;
    }
    if (mcTrk0->getProcess() != kPPair || mcTrk1->getProcess() != kPPair) {
      return false;
    }
    return true;
  }

  bool checkITS(GIndex const& gid0, GIndex const& gid1)
  {
    auto gid0ITS = gid0.includesDet(o2::detectors::DetID::ITS);
    auto gid1ITS = gid1.includesDet(o2::detectors::DetID::ITS);
    if (gid0ITS && gid1ITS) {
      return true;
    }
    return false;
  }

  bool checkTPC(GIndex const& gid0, GIndex const& gid1)
  {
    if (checkITS(gid0, gid1)) {
      return false;
    }
    auto gid0TPC = gid0.includesDet(o2::detectors::DetID::TPC);
    auto gid1TPC = gid1.includesDet(o2::detectors::DetID::TPC);
    if (gid0TPC && gid1TPC) {
      return true;
    }
    return false;
  }

  bool checkITSTPC(GIndex const& gid0, GIndex const& gid1)
  {
    auto gid0ITS = gid0.includesDet(o2::detectors::DetID::ITS);
    auto gid1ITS = gid1.includesDet(o2::detectors::DetID::ITS);
    auto gid0TPC = gid0.includesDet(o2::detectors::DetID::TPC);
    auto gid1TPC = gid1.includesDet(o2::detectors::DetID::TPC);
    if (gid0ITS && gid1ITS && gid0TPC && gid1TPC) {
      return true;
    }
    return false;
  }
};
} // namespace vertexing
} // namespace o2

#endif
