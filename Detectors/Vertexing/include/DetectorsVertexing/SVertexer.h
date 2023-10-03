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
#include <utility>
#include "GPUO2InterfaceRefit.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "SimulationDataFormat/MCUtils.h"
#include "Steer/MCKinematicsReader.h"
#include "boost/functional/hash.hpp"
#include "TPDGCode.h"
#include "TArrayD.h"

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
  using Vec3D = ROOT::Math::SVector<double, 3>;

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
  void setMeanVertex(const o2d::MeanVertexObject* v)
  {
    if (v == nullptr) {
      return;
    }
    mMeanVertex = v->getMeanVertex();
  }
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

  gsl::span<const o2::MCCompLabel> mITSTrkLabels;
  gsl::span<const o2::MCCompLabel> mTPCTrkLabels;
  gsl::span<const o2::MCCompLabel> mITSTPCTrkLabels;
  gsl::span<const o2::MCCompLabel> mITSTPCTOFTrkLabels;
  gsl::span<const o2::MCCompLabel> mITSTPCTRDTrkLabels;
  gsl::span<const o2::MCCompLabel> mITSTPCTRDTOFTrkLabels;
  gsl::span<const o2::MCCompLabel> mTPCTRDTOFTrkLabels;
  gsl::span<const o2::MCCompLabel> mTPCTRDTrkLabels;
  gsl::span<const o2::MCCompLabel> mTPCTOFTrkLabels;
  o2::utils::TreeStreamRedirector mDebugStream{"svertexer-debug.root", "recreate"};
  o2::steer::MCKinematicsReader mcReader; // reader of MC information
  void writeDebugV0Candidates(o2::tpc::TrackTPC const& trk, GIndex gid, int vtxid, o2::track::TrackParCov const& candTrk);
  void writeDebugWithoutTiming(const o2::globaltracking::RecoContainer& recoData);
  void writeDebugWithTiming(const o2::globaltracking::RecoContainer& recoData);
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
  using map_timing_t = std::unordered_map<key_t, std::tuple<GIndex, GIndex, bool>, key_hash>;
  using map_before_t = std::unordered_map<key_t, std::tuple<MCTrack, MCTrack, MCTrack, bool>, key_hash>;
  using map_after_t = std::unordered_map<key_t, std::tuple<TrackCand, MCTrack, TrackCand, MCTrack, MCTrack>, key_hash>;
  using map_mc_t = std::unordered_map<key_t, std::pair<key_t, key_t>, key_hash>;
  using map_mc_particle_t = std::unordered_map<key_t, bool, key_hash>;

  map_mc_t mD0V0Map;
  map_mc_t mD1V0Map;
  map_mc_t mMotherV0Map;
  map_mc_particle_t mMCParticle;

  using key_dup_t = std::tuple<GIndex, GIndex>;

  struct key_dup_hash {
    size_t operator()(const key_dup_t& k) const
    {
      const auto& [gid0, gid1] = k;
      size_t seed = 0;
      boost::hash_combine(seed, gid0.getRaw());
      boost::hash_combine(seed, gid1.getRaw());
      return seed;
    }
  };
  using map_dup_t = std::unordered_map<key_dup_t, ULong64_t, key_dup_hash>;

  static bool checkMother(o2::MCTrack const* mother, const std::vector<o2::MCTrack>& pcontainer)
  {
    if (mother == nullptr) {
      return false;
    }
    bool isPhysicalPrimary = o2::mcutils::MCTrackNavigator::isPhysicalPrimary(*mother, pcontainer);
    return mother->GetPdgCode() == v0Type && mother->isPrimary() && isPhysicalPrimary;
  }

  o2::MCCompLabel getLabel(GIndex const& gid)
  {
    if (gid.getSource() == GIndex::ITSTPCTRDTOF) {
      return mITSTPCTRDTOFTrkLabels[gid.getIndex()];
    } else if (gid.getSource() == GIndex::ITSTPCTOF) {
      return mITSTPCTOFTrkLabels[gid.getIndex()];
    } else if (gid.getSource() == GIndex::ITSTPCTRD) {
      return mITSTPCTRDTrkLabels[gid.getIndex()];
    } else if (gid.getSource() == GIndex::TPCTRD) {
      return mTPCTRDTrkLabels[gid.getIndex()];
    } else if (gid.getSource() == GIndex::TPCTOF) {
      return mTPCTOFTrkLabels[gid.getIndex()];
    } else if (gid.getSource() == GIndex::ITSTPC) {
      return mITSTPCTrkLabels[gid.getIndex()];
    } else if (gid.getSource() == GIndex::ITS) {
      return mITSTrkLabels[gid.getIndex()];
    } else if (gid.getSource() == GIndex::TPC) {
      return mTPCTrkLabels[gid.getIndex()];
    }
    return {};
  }

  static bool checkLabels(o2::MCCompLabel const& lbl0, o2::MCCompLabel const& lbl1)
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

  static bool checkPair(o2::MCTrack const* mcTrk0, o2::MCTrack const* mcTrk1)
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

  template <typename Enum>
  struct Counter_t {
    const std::array<std::string_view, static_cast<size_t>(Enum::NSIZE)>& _names;
    const std::string _treeName;
    Counter_t(std::string_view const& treeName, std::array<std::string_view, static_cast<size_t>(Enum::NSIZE)> const& names) : _treeName{treeName}, _names{names}
    {
      mTotCounters.fill(0);
      for (auto& c : mCounters) {
        c.fill(0);
      }
    }

    std::array<ULong64_t, static_cast<size_t>(Enum::NSIZE)> mTotCounters{};
    std::array<std::array<ULong64_t, static_cast<size_t>(Enum::NSIZE)>, 55> mCounters{};
    std::array<std::array<ULong64_t, static_cast<size_t>(Enum::NSIZE)>, 55> mCountersDup{};
    std::array<std::array<std::unordered_map<GIndex, bool>, static_cast<size_t>(Enum::NSIZE)>, 55> mDupMap{};
    std::array<std::array<ULong64_t, static_cast<size_t>(Enum::NSIZE)>, 55> mCountersV0{};
    std::array<std::array<ULong64_t, static_cast<size_t>(Enum::NSIZE)>, 55> mCountersV0Dup{};
    std::array<std::array<ULong64_t, static_cast<size_t>(Enum::NSIZE)>, 55> mCountersMC{};
    std::array<std::array<map_dup_t, static_cast<size_t>(Enum::NSIZE)>, 55> mDup2Map{};

    ULong64_t mCounterMC{0};

    void inc(Enum e, GIndex const& gid0, GIndex const& gid1)
    {
      auto c = static_cast<unsigned int>(e);
      auto gid0ITS = gid0.includesDet(o2::detectors::DetID::ITS);
      auto gid0TPC = gid0.includesDet(o2::detectors::DetID::TPC);
      auto gid0TRD = gid0.includesDet(o2::detectors::DetID::TRD);
      auto gid0TOF = gid0.includesDet(o2::detectors::DetID::TOF);
      auto gid1ITS = gid1.includesDet(o2::detectors::DetID::ITS);
      auto gid1TPC = gid1.includesDet(o2::detectors::DetID::TPC);
      auto gid1TRD = gid1.includesDet(o2::detectors::DetID::TRD);
      auto gid1TOF = gid1.includesDet(o2::detectors::DetID::TOF);
      int i = getCombination(gid0ITS, gid0TPC, gid0TRD, gid0TOF, gid1ITS, gid1TPC, gid1TRD, gid1TOF);
      ++mCounters[i][c];
      ++mTotCounters[c];
    }
    void inc(Enum e, PVertex const& pvertex, std::array<float, 3> const& svertex, const o2::track::TrackParCov& seedP, const o2::track::TrackParCov& seedN, GIndex const& gid0, GIndex const& gid1, o2::MCCompLabel const& lbl0, o2::MCCompLabel const& lbl1, bool checkLabels, map_mc_t const& d0, map_mc_t const& d1, o2::steer::MCKinematicsReader& mcReader, utils::TreeStreamRedirector& mDebugStream)
    {
      // foundV0
      auto c = static_cast<unsigned int>(e);
      bool duplicate{false};
      bool trueV0{false};
      auto gid0ITS = gid0.includesDet(o2::detectors::DetID::ITS);
      auto gid0TPC = gid0.includesDet(o2::detectors::DetID::TPC);
      auto gid0TRD = gid0.includesDet(o2::detectors::DetID::TRD);
      auto gid0TOF = gid0.includesDet(o2::detectors::DetID::TOF);
      auto gid1ITS = gid1.includesDet(o2::detectors::DetID::ITS);
      auto gid1TPC = gid1.includesDet(o2::detectors::DetID::TPC);
      auto gid1TRD = gid1.includesDet(o2::detectors::DetID::TRD);
      auto gid1TOF = gid1.includesDet(o2::detectors::DetID::TOF);
      int i = getCombination(gid0ITS, gid0TPC, gid0TRD, gid0TOF, gid1ITS, gid1TPC, gid1TRD, gid1TOF);
      ++mCounters[i][c];
      ++mTotCounters[c];
      auto idx = std::make_tuple(gid0, gid1);
      auto dup = mDup2Map[i][c].find(idx) != mDup2Map[i][c].end();
      if (dup) {
        duplicate = true;
        ++mCountersDup[i][c];
      } else {
        ++mDup2Map[i][c][idx];
      }
      if (checkLabels) {
        auto idx0 = std::make_tuple(lbl0.getSourceID(), lbl0.getEventID(), lbl0.getTrackID());
        auto it0 = d0.find(idx0);
        auto it0T = it0 != d0.end();
        auto idx1 = std::make_tuple(lbl1.getSourceID(), lbl1.getEventID(), lbl1.getTrackID());
        auto it1 = d1.find(idx1);
        auto it1T = it1 != d1.end();
        if (it0T && it1T) {
          if ((*it0).second.second == (*it1).second.second) {
            trueV0 = true;
            ++mCountersV0[i][c];
            if (dup) {
              ++mCountersV0Dup[i][c];
            }
          }
        }
      }
      const MCTrack *mcTrk0 = mcReader.getTrack(lbl0), *mcTrk1 = mcReader.getTrack(lbl1);
      if (mcTrk0 == nullptr || mcTrk1 == nullptr) {
        return;
      }
      std::array<float, 3> xyz0{(float)mcTrk0->GetStartVertexCoordinatesX(), (float)mcTrk0->GetStartVertexCoordinatesY(), (float)mcTrk0->GetStartVertexCoordinatesZ()};
      std::array<float, 3> pxyz0{(float)mcTrk0->GetStartVertexMomentumX(), (float)mcTrk0->GetStartVertexMomentumY(), (float)mcTrk0->GetStartVertexMomentumZ()};
      std::array<float, 3> xyz1{(float)mcTrk1->GetStartVertexCoordinatesX(), (float)mcTrk1->GetStartVertexCoordinatesY(), (float)mcTrk1->GetStartVertexCoordinatesZ()};
      std::array<float, 3> pxyz1{(float)mcTrk1->GetStartVertexMomentumX(), (float)mcTrk1->GetStartVertexMomentumY(), (float)mcTrk1->GetStartVertexMomentumZ()};
      auto pPDG0 = TDatabasePDG::Instance()->GetParticle(mcTrk0->GetPdgCode());
      auto pPDG1 = TDatabasePDG::Instance()->GetParticle(mcTrk1->GetPdgCode());
      if (pPDG0 == nullptr || pPDG1 == nullptr) {
        return;
      }
      o2::track::TrackPar mctr0(xyz0, pxyz0, TMath::Nint(pPDG0->Charge() / 3), false);
      o2::track::TrackPar mctr1(xyz1, pxyz1, TMath::Nint(pPDG1->Charge() / 3), false);
      if (!mctr0.rotate(seedP.getAlpha()) || !o2::base::Propagator::Instance()->PropagateToXBxByBz(mctr0, seedP.getX()) ||
          !mctr1.rotate(seedN.getAlpha()) || !o2::base::Propagator::Instance()->PropagateToXBxByBz(mctr1, seedN.getX())) {
        return;
      }

      mDebugStream << _treeName.c_str()
                   << "recoSeedP=" << seedP
                   << "recoSeedN=" << seedN
                   << "mcSeedP=" << mctr0
                   << "mcSeedN=" << mctr1
                   << "svertex=" << svertex
                   << "pvertex=" << pvertex
                   << "case=" << c
                   << "isDuplicate=" << duplicate
                   << "isV0=" << trueV0
                   << "\n";
    }

    void inc(Enum e, PVertex const& pvertex, Vec3D const& svertex, const TrackCand& seedP, const TrackCand& seedN, o2::MCCompLabel const& lbl0, o2::MCCompLabel const& lbl1, bool checkLabels, map_mc_t const& d0, map_mc_t const& d1, o2::steer::MCKinematicsReader& mcReader, utils::TreeStreamRedirector& mDebugStream, bool write = true, bool useMC = false)
    {
      // checkV0
      auto c = static_cast<unsigned int>(e);
      bool duplicate{false};
      bool trueV0{false};
        TArrayD sv{3};
        sv.SetAt(svertex[0], 0);
        sv.SetAt(svertex[1], 1);
        sv.SetAt(svertex[2], 2);
      auto gid0ITS = seedP.gid.includesDet(o2::detectors::DetID::ITS);
      auto gid0TPC = seedP.gid.includesDet(o2::detectors::DetID::TPC);
      auto gid0TRD = seedP.gid.includesDet(o2::detectors::DetID::TRD);
      auto gid0TOF = seedP.gid.includesDet(o2::detectors::DetID::TOF);
      auto gid1ITS = seedN.gid.includesDet(o2::detectors::DetID::ITS);
      auto gid1TPC = seedN.gid.includesDet(o2::detectors::DetID::TPC);
      auto gid1TRD = seedN.gid.includesDet(o2::detectors::DetID::TRD);
      auto gid1TOF = seedN.gid.includesDet(o2::detectors::DetID::TOF);
      int i = getCombination(gid0ITS, gid0TPC, gid0TRD, gid0TOF, gid1ITS, gid1TPC, gid1TRD, gid1TOF);
      ++mCounters[i][c];
      ++mTotCounters[c];
      auto idx = std::make_tuple(seedP.gid, seedN.gid);
      auto dup = mDup2Map[i][c].find(idx) != mDup2Map[i][c].end();
      if (dup) {
        duplicate = true;
        ++mCountersDup[i][c];
      } else {
        ++mDup2Map[i][c][idx];
      }
      if (checkLabels) {
        auto idx0 = std::make_tuple(lbl0.getSourceID(), lbl0.getEventID(), lbl0.getTrackID());
        auto it0 = d0.find(idx0);
        auto it0T = it0 != d0.end();
        auto idx1 = std::make_tuple(lbl1.getSourceID(), lbl1.getEventID(), lbl1.getTrackID());
        auto it1 = d1.find(idx1);
        auto it1T = it1 != d1.end();
        if (it0T && it1T) {
          if ((*it0).second.second == (*it1).second.second) {
            trueV0 = true;
            ++mCountersV0[i][c];
            if (dup) {
              ++mCountersV0Dup[i][c];
            }
          }
        }
      }
      if (!write) {
        return;
      }

      mDebugStream << _treeName.c_str()
                   << "recoSeedP=" << (o2::track::TrackParCov)seedP
                   << "recoSeedN=" << (o2::track::TrackParCov)seedN
                   << "svertex=" << sv
                   << "pvertex=" << pvertex
                   << "case=" << c
                   << "isDuplicate=" << duplicate
                   << "isV0=" << trueV0
                   << "\n";
    }

    bool inc(Enum e, bool gid0ITS, bool gid0TPC, bool gid0TRD, bool gid0TOF)
    {
      auto c = static_cast<unsigned int>(e);
      int i = getType(gid0ITS, gid0TPC, gid0TRD, gid0TOF);
      ++mCountersMC[i][c];
      ++mCounterMC;
      return i == 9;
    }

    bool inc(Enum e, bool gid0ITS, bool gid0TPC, bool gid0TRD, bool gid0TOF, bool gid1ITS, bool gid1TPC, bool gid1TRD, bool gid1TOF)
    {
      auto c = static_cast<unsigned int>(e);
      int i = getCombination(gid0ITS, gid0TPC, gid0TRD, gid0TOF, gid1ITS, gid1TPC, gid1TRD, gid1TOF);
      ++mCounters[i][c];
      ++mTotCounters[c];
      return i == 54;
    }

    int getCombination(bool gid0ITS, bool gid0TPC, bool gid0TRD, bool gid0TOF, bool gid1ITS, bool gid1TPC, bool gid1TRD, bool gid1TOF)
    {
      int i = -1;
      if (gid0ITS && gid0TPC && gid0TRD && gid0TOF && gid1ITS && gid1TPC && gid1TRD && gid1TOF) {
        i = 0;
      } else if ((gid0ITS && gid0TPC && !gid0TRD && gid0TOF && gid1ITS && gid1TPC && gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && gid0TOF && gid1ITS && gid1TPC && !gid1TRD && gid1TOF)) {
        i = 1;
      } else if ((gid0ITS && gid0TPC && gid0TRD && !gid0TOF && gid1ITS && gid1TPC && gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && gid0TOF && gid1ITS && gid1TPC && gid1TRD && !gid1TOF)) {
        i = 2;
      } else if ((!gid0ITS && gid0TPC && gid0TRD && gid0TOF && gid1ITS && gid1TPC && gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && gid0TOF && !gid1ITS && gid1TPC && gid1TRD && gid1TOF)) {
        i = 3;
      } else if ((!gid0ITS && gid0TPC && gid0TRD && !gid0TOF && gid1ITS && gid1TPC && gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && gid0TOF && !gid1ITS && gid1TPC && gid1TRD && !gid1TOF)) {
        i = 4;
      } else if ((!gid0ITS && gid0TPC && !gid0TRD && gid0TOF && gid1ITS && gid1TPC && gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && gid1TOF)) {
        i = 5;
      } else if ((gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && gid0TOF && gid1ITS && gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 6;
      } else if ((gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && gid0TOF && gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 7;
      } else if ((!gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 8;
      } else if ((!gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && gid0TOF && !gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 9;
      } else if (gid0ITS && gid0TPC && !gid0TRD && gid0TOF && gid1ITS && gid1TPC && !gid1TRD && gid1TOF) {
        i = 10;
      } else if ((gid0ITS && gid0TPC && gid0TRD && !gid0TOF && gid1ITS && gid1TPC && !gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && !gid0TRD && gid0TOF && gid1ITS && gid1TPC && gid1TRD && !gid1TOF)) {
        i = 11;
      } else if ((!gid0ITS && gid0TPC && gid0TRD && gid0TOF && gid1ITS && gid1TPC && !gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && !gid0TRD && gid0TOF && !gid1ITS && gid1TPC && gid1TRD && gid1TOF)) {
        i = 12;
      } else if ((!gid0ITS && gid0TPC && gid0TRD && !gid0TOF && gid1ITS && gid1TPC && !gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && !gid0TRD && gid0TOF && !gid1ITS && gid1TPC && gid1TRD && !gid1TOF)) {
        i = 13;
      } else if ((!gid0ITS && gid0TPC && !gid0TRD && gid0TOF && gid1ITS && gid1TPC && !gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && !gid0TRD && gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && gid1TOF)) {
        i = 14;
      } else if ((gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && !gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && !gid0TRD && gid0TOF && gid1ITS && gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 15;
      } else if ((gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && !gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && !gid0TRD && gid0TOF && gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 16;
      } else if ((!gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && !gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && !gid0TRD && gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 17;
      } else if ((!gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && !gid1TRD && gid1TOF) ||
                 (gid0ITS && gid0TPC && !gid0TRD && gid0TOF && !gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 18;
      } else if (gid0ITS && gid0TPC && gid0TRD && !gid0TOF && gid1ITS && gid1TPC && gid1TRD && !gid1TOF) {
        i = 19;
      } else if ((!gid0ITS && gid0TPC && gid0TRD && gid0TOF && gid1ITS && gid1TPC && gid1TRD && !gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && gid1TRD && gid1TOF)) {
        i = 20;
      } else if ((!gid0ITS && gid0TPC && gid0TRD && !gid0TOF && gid1ITS && gid1TPC && gid1TRD && !gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && gid1TRD && !gid1TOF)) {
        i = 21;
      } else if ((!gid0ITS && gid0TPC && !gid0TRD && gid0TOF && gid1ITS && gid1TPC && gid1TRD && !gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && gid1TOF)) {
        i = 22;
      } else if ((gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && gid1TRD && !gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && !gid0TOF && gid1ITS && gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 23;
      } else if ((gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && gid1TRD && !gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && !gid0TOF && gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 24;
      } else if ((!gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && gid1TRD && !gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 25;
      } else if ((!gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && gid1TRD && !gid1TOF) ||
                 (gid0ITS && gid0TPC && gid0TRD && !gid0TOF && !gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 26;
      } else if (!gid0ITS && gid0TPC && gid0TRD && gid0TOF && !gid1ITS && gid1TPC && gid1TRD && gid1TOF) {
        i = 27;
      } else if ((!gid0ITS && gid0TPC && gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && gid1TRD && gid1TOF) ||
                 (!gid0ITS && gid0TPC && gid0TRD && gid0TOF && !gid1ITS && gid1TPC && gid1TRD && !gid1TOF)) {
        i = 28;
      } else if ((!gid0ITS && gid0TPC && !gid0TRD && gid0TOF && !gid1ITS && gid1TPC && gid1TRD && gid1TOF) ||
                 (!gid0ITS && gid0TPC && gid0TRD && gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && gid1TOF)) {
        i = 29;
      } else if ((gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && gid1TRD && gid1TOF) ||
                 (!gid0ITS && gid0TPC && gid0TRD && gid0TOF && gid1ITS && gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 30;
      } else if ((gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && gid1TRD && gid1TOF) ||
                 (!gid0ITS && gid0TPC && gid0TRD && gid0TOF && gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 31;
      } else if ((!gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && gid1TRD && gid1TOF) ||
                 (!gid0ITS && gid0TPC && gid0TRD && gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 32;
      } else if ((!gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && gid1TRD && gid1TOF) ||
                 (!gid0ITS && gid0TPC && gid0TRD && gid0TOF && !gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 33;
      } else if (!gid0ITS && gid0TPC && gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && gid1TRD && !gid1TOF) {
        i = 34;
      } else if ((!gid0ITS && gid0TPC && !gid0TRD && gid0TOF && !gid1ITS && gid1TPC && gid1TRD && !gid1TOF) ||
                 (!gid0ITS && gid0TPC && gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && gid1TOF)) {
        i = 35;
      } else if ((gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && gid1TRD && !gid1TOF) ||
                 (!gid0ITS && gid0TPC && gid0TRD && !gid0TOF && gid1ITS && gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 36;
      } else if ((gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && gid1TRD && !gid1TOF) ||
                 (!gid0ITS && gid0TPC && gid0TRD && !gid0TOF && gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 37;
      } else if ((!gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && gid1TRD && !gid1TOF) ||
                 (!gid0ITS && gid0TPC && gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 38;
      } else if ((!gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && gid1TRD && !gid1TOF) ||
                 (!gid0ITS && gid0TPC && gid0TRD && !gid0TOF && !gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 39;
      } else if (!gid0ITS && gid0TPC && !gid0TRD && gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && gid1TOF) {
        i = 40;
      } else if ((gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && gid1TOF) ||
                 (!gid0ITS && gid0TPC && !gid0TRD && gid0TOF && gid1ITS && gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 41;
      } else if ((gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && gid1TOF) ||
                 (!gid0ITS && gid0TPC && !gid0TRD && gid0TOF && gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 42;
      } else if ((!gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && gid1TOF) ||
                 (!gid0ITS && gid0TPC && !gid0TRD && gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 43;
      } else if ((!gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && gid1TOF) ||
                 (!gid0ITS && gid0TPC && !gid0TRD && gid0TOF && !gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 44;
      } else if (gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && !gid1TRD && !gid1TOF) {
        i = 45;
      } else if ((gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && !gid1TRD && !gid1TOF) ||
                 (gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 46;
      } else if ((!gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && !gid1TRD && !gid1TOF) ||
                 (gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 47;
      } else if ((!gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && gid1TPC && !gid1TRD && !gid1TOF) ||
                 (gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 48;
      } else if (gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF) {
        i = 49;
      } else if ((!gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF) ||
                 (gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 50;
      } else if ((!gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF) ||
                 (gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 51;
      } else if (!gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && !gid1TOF) {
        i = 52;
      } else if ((!gid0ITS && !gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && gid1TPC && !gid1TRD && !gid1TOF) ||
                 (!gid0ITS && gid0TPC && !gid0TRD && !gid0TOF && !gid1ITS && !gid1TPC && !gid1TRD && !gid1TOF)) {
        i = 53;
      } else {
        i = 54;
      }
      return i;
    }

    int getType(bool gidITS, bool gidTPC, bool gidTRD, bool gidTOF)
    {
      int i = -1;
      if (gidITS && gidTPC && gidTRD && gidTOF) {
        i = 0;
      } else if (gidITS && gidTPC && !gidTRD && gidTOF) {
        i = 1;
      } else if (gidITS && gidTPC && gidTRD && !gidTOF) {
        i = 2;
      } else if (!gidITS && gidTPC && gidTRD && gidTOF) {
        i = 3;
      } else if (!gidITS && gidTPC && gidTRD && !gidTOF) {
        i = 4;
      } else if (!gidITS && gidTPC && !gidTRD && gidTOF) {
        i = 5;
      } else if (gidITS && gidTPC && !gidTRD && !gidTOF) {
        i = 6;
      } else if (gidITS && !gidTPC && !gidTRD && !gidTOF) {
        i = 7;
      } else if (!gidITS && gidTPC && !gidTRD && !gidTOF) {
        i = 8;
      } else {
        i = 9;
      }
      return i;
    }

    void inc(Enum e, GIndex const& gid, o2::MCCompLabel const& lbl, map_mc_t const& d0, map_mc_t const& d1, map_mc_particle_t const& mcparticles)
    {
      auto c = static_cast<unsigned int>(e);
      auto gidITS = gid.includesDet(o2::detectors::DetID::ITS);
      auto gidTPC = gid.includesDet(o2::detectors::DetID::TPC);
      auto gidTRD = gid.includesDet(o2::detectors::DetID::TRD);
      auto gidTOF = gid.includesDet(o2::detectors::DetID::TOF);
      int i = getType(gidITS, gidTPC, gidTRD, gidTOF);
      ++mCounters[i][c];
      ++mTotCounters[c];
      auto dup = mDupMap[i][c].find(gid) != mDupMap[i][c].end();
      if (dup) {
        ++mCountersDup[i][c];
      } else {
        mDupMap[i][c][gid] = true;
      }
      if (lbl.isFake() || !lbl.isValid()) {
        return;
      }
      auto idx = std::make_tuple(lbl.getSourceID(), lbl.getEventID(), lbl.getTrackID());
      auto it0 = d0.find(idx);
      auto it1 = d1.find(idx);
      auto it0T = it0 != d0.end();
      auto it1T = it1 != d1.end();
      if (it0T || it1T) {
        ++mCountersV0[i][c];
        if (dup) {
          ++mCountersV0Dup[i][c];
        }
      }
      auto mcgen = mcparticles.find(idx) != mcparticles.end();
      if (mcgen) {
        ++mCountersMC[i][c];
      }
    }

    void print()
    {
      for (int i{0}; i < static_cast<int>(Enum::NSIZE); ++i) {
        LOGP(info, "{} - CHECK: {}: {}", i, _names[i], mTotCounters[i]);
        LOGP(info, "                 `--> Track           {:>10} / {:>10} / {:>10} / {:>10} / {:>10}", "true V0", "mcparticle", "reco", "V0 dup", "dup");
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF {:>10} / {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[0][i], mCountersMC[0][i], mCounters[0][i], mCountersV0Dup[0][i], mCountersDup[0][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF {:>10} / {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[1][i], mCountersMC[1][i], mCounters[1][i], mCountersV0Dup[1][i], mCountersDup[1][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     {:>10} / {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[2][i], mCountersMC[2][i], mCounters[2][i], mCountersV0Dup[2][i], mCountersDup[2][i]);
        LOGP(info, "                 `-->     TPC-TRD-TOF {:>10} / {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[3][i], mCountersMC[3][i], mCounters[3][i], mCountersV0Dup[3][i], mCountersDup[3][i]);
        LOGP(info, "                 `-->     TPC-TRD     {:>10} / {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[4][i], mCountersMC[4][i], mCounters[4][i], mCountersV0Dup[4][i], mCountersDup[4][i]);
        LOGP(info, "                 `-->     TPC-   -TOF {:>10} / {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[5][i], mCountersMC[5][i], mCounters[5][i], mCountersV0Dup[5][i], mCountersDup[5][i]);
        LOGP(info, "                 `--> ITS-TPC         {:>10} / {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[6][i], mCountersMC[6][i], mCounters[6][i], mCountersV0Dup[6][i], mCountersDup[6][i]);
        LOGP(info, "                 `--> ITS             {:>10} / {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[7][i], mCountersMC[7][i], mCounters[7][i], mCountersV0Dup[7][i], mCountersDup[7][i]);
        LOGP(info, "                 `-->     TPC         {:>10} / {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[8][i], mCountersMC[8][i], mCounters[8][i], mCountersV0Dup[8][i], mCountersDup[8][i]);
        LOGP(info, "                 `--> ???             {:>10} / {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[9][i], mCountersMC[9][i], mCounters[9][i], mCountersV0Dup[9][i], mCountersDup[9][i]);
      }
    }
    void print2()
    {
      for (int i{0}; i < static_cast<int>(Enum::NSIZE); ++i) {
        LOGP(info, "{} - CHECK: {}: {}", i, _names[i], mTotCounters[i]);
        LOGP(info, "                 `--> Track           & Track           {:>10} / {:>10} / {:>10} / {:>10}", "true V0s", "reco", "V0 dup", "dup");
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF & ITS-TPC-TRD-TOF {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[0][i], mCounters[0][i], mCountersV0Dup[0][i], mCountersDup[0][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF & ITS-TPC-   -TOF {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[1][i], mCounters[1][i], mCountersV0Dup[1][i], mCountersDup[1][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF & ITS-TPC-TRD     {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[2][i], mCounters[2][i], mCountersV0Dup[2][i], mCountersDup[2][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF &     TPC-TRD-TOF {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[3][i], mCounters[3][i], mCountersV0Dup[3][i], mCountersDup[3][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF &     TPC-TRD     {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[4][i], mCounters[4][i], mCountersV0Dup[4][i], mCountersDup[4][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF &     TPC-   -TOF {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[5][i], mCounters[5][i], mCountersV0Dup[5][i], mCountersDup[5][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF & ITS-TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[6][i], mCounters[6][i], mCountersV0Dup[6][i], mCountersDup[6][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF & ITS             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[7][i], mCounters[7][i], mCountersV0Dup[7][i], mCountersDup[7][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF &     TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[8][i], mCounters[8][i], mCountersV0Dup[8][i], mCountersDup[8][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF & ???             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[9][i], mCounters[9][i], mCountersV0Dup[9][i], mCountersDup[9][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF & ITS-TPC-   -TOF {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[10][i], mCounters[10][i], mCountersV0Dup[10][i], mCountersDup[10][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF & ITS-TPC-TRD     {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[11][i], mCounters[11][i], mCountersV0Dup[11][i], mCountersDup[11][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF &     TPC-TRD-TOF {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[12][i], mCounters[12][i], mCountersV0Dup[12][i], mCountersDup[12][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF &     TPC-TRD     {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[13][i], mCounters[13][i], mCountersV0Dup[13][i], mCountersDup[13][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF &     TPC-   -TOF {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[14][i], mCounters[14][i], mCountersV0Dup[14][i], mCountersDup[14][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF & ITS-TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[15][i], mCounters[15][i], mCountersV0Dup[15][i], mCountersDup[15][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF & ITS             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[16][i], mCounters[16][i], mCountersV0Dup[16][i], mCountersDup[16][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF &     TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[17][i], mCounters[17][i], mCountersV0Dup[17][i], mCountersDup[17][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF & ???             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[18][i], mCounters[18][i], mCountersV0Dup[18][i], mCountersDup[18][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     & ITS-TPC-TRD     {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[19][i], mCounters[19][i], mCountersV0Dup[19][i], mCountersDup[19][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     &     TPC-TRD-TOF {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[20][i], mCounters[20][i], mCountersV0Dup[20][i], mCountersDup[20][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     &     TPC-TRD     {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[21][i], mCounters[21][i], mCountersV0Dup[21][i], mCountersDup[21][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     &     TPC-   -TOF {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[22][i], mCounters[22][i], mCountersV0Dup[22][i], mCountersDup[22][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     & ITS-TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[23][i], mCounters[23][i], mCountersV0Dup[23][i], mCountersDup[23][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     & ITS             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[24][i], mCounters[24][i], mCountersV0Dup[24][i], mCountersDup[24][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     &     TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[25][i], mCounters[25][i], mCountersV0Dup[25][i], mCountersDup[25][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     & ???             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[26][i], mCounters[26][i], mCountersV0Dup[26][i], mCountersDup[26][i]);
        LOGP(info, "                 `-->     TPC-TRD-TOF &     TPC-TRD-TOF {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[27][i], mCounters[27][i], mCountersV0Dup[27][i], mCountersDup[27][i]);
        LOGP(info, "                 `-->     TPC-TRD-TOF &     TPC-TRD     {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[28][i], mCounters[28][i], mCountersV0Dup[28][i], mCountersDup[28][i]);
        LOGP(info, "                 `-->     TPC-TRD-TOF &     TPC-   -TOF {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[29][i], mCounters[29][i], mCountersV0Dup[29][i], mCountersDup[29][i]);
        LOGP(info, "                 `-->     TPC-TRD-TOF & ITS-TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[30][i], mCounters[30][i], mCountersV0Dup[30][i], mCountersDup[30][i]);
        LOGP(info, "                 `-->     TPC-TRD-TOF & ITS             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[31][i], mCounters[31][i], mCountersV0Dup[31][i], mCountersDup[31][i]);
        LOGP(info, "                 `-->     TPC-TRD-TOF &     TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[32][i], mCounters[32][i], mCountersV0Dup[32][i], mCountersDup[32][i]);
        LOGP(info, "                 `-->     TPC-TRD-TOF & ???             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[33][i], mCounters[33][i], mCountersV0Dup[33][i], mCountersDup[33][i]);
        LOGP(info, "                 `-->     TPC-TRD     &     TPC-TRD     {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[34][i], mCounters[34][i], mCountersV0Dup[34][i], mCountersDup[34][i]);
        LOGP(info, "                 `-->     TPC-TRD     &     TPC-   -TOF {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[35][i], mCounters[35][i], mCountersV0Dup[35][i], mCountersDup[35][i]);
        LOGP(info, "                 `-->     TPC-TRD     & ITS-TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[36][i], mCounters[36][i], mCountersV0Dup[36][i], mCountersDup[36][i]);
        LOGP(info, "                 `-->     TPC-TRD     & ITS             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[37][i], mCounters[37][i], mCountersV0Dup[37][i], mCountersDup[37][i]);
        LOGP(info, "                 `-->     TPC-TRD     &     TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[38][i], mCounters[38][i], mCountersV0Dup[38][i], mCountersDup[38][i]);
        LOGP(info, "                 `-->     TPC-TRD     & ???             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[39][i], mCounters[39][i], mCountersV0Dup[39][i], mCountersDup[39][i]);
        LOGP(info, "                 `-->     TPC-   -TOF &     TPC-   -TOF {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[40][i], mCounters[40][i], mCountersV0Dup[40][i], mCountersDup[40][i]);
        LOGP(info, "                 `-->     TPC-   -TOF & ITS-TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[41][i], mCounters[41][i], mCountersV0Dup[41][i], mCountersDup[41][i]);
        LOGP(info, "                 `-->     TPC-   -TOF & ITS             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[42][i], mCounters[42][i], mCountersV0Dup[42][i], mCountersDup[42][i]);
        LOGP(info, "                 `-->     TPC-   -TOF &     TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[43][i], mCounters[43][i], mCountersV0Dup[43][i], mCountersDup[43][i]);
        LOGP(info, "                 `-->     TPC-   -TOF & ???             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[44][i], mCounters[44][i], mCountersV0Dup[44][i], mCountersDup[44][i]);
        LOGP(info, "                 `--> ITS-TPC         & ITS-TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[45][i], mCounters[45][i], mCountersV0Dup[45][i], mCountersDup[45][i]);
        LOGP(info, "                 `--> ITS-TPC         & ITS             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[46][i], mCounters[46][i], mCountersV0Dup[46][i], mCountersDup[46][i]);
        LOGP(info, "                 `--> ITS-TPC         &     TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[47][i], mCounters[47][i], mCountersV0Dup[47][i], mCountersDup[47][i]);
        LOGP(info, "                 `--> ITS-TPC         & ???             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[48][i], mCounters[48][i], mCountersV0Dup[48][i], mCountersDup[48][i]);
        LOGP(info, "                 `--> ITS             & ITS             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[49][i], mCounters[49][i], mCountersV0Dup[49][i], mCountersDup[49][i]);
        LOGP(info, "                 `--> ITS             &     TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[50][i], mCounters[50][i], mCountersV0Dup[50][i], mCountersDup[50][i]);
        LOGP(info, "                 `--> ITS             & ???             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[51][i], mCounters[51][i], mCountersV0Dup[51][i], mCountersDup[51][i]);
        LOGP(info, "                 `-->     TPC         &     TPC         {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[52][i], mCounters[52][i], mCountersV0Dup[52][i], mCountersDup[52][i]);
        LOGP(info, "                 `-->     TPC         & ???             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[53][i], mCounters[53][i], mCountersV0Dup[53][i], mCountersDup[53][i]);
        LOGP(info, "                 `--> ???             & ???             {:>10} / {:>10} / {:>10} / {:>10}", mCountersV0[54][i], mCounters[54][i], mCountersV0Dup[54][i], mCountersDup[54][i]);
      }
    }

    void printMC()
    {
      int i{0};
      LOGP(info, "{} - CHECK: {}: {}", i, "MC Particles generated", mCounterMC);
      LOGP(info, "                 `--> Track           {:>10}", "MCParticles");
      LOGP(info, "                 `--> ITS-TPC-TRD-TOF {:>10}", mCountersMC[0][i]);
      LOGP(info, "                 `--> ITS-TPC-   -TOF {:>10}", mCountersMC[1][i]);
      LOGP(info, "                 `--> ITS-TPC-TRD     {:>10}", mCountersMC[2][i]);
      LOGP(info, "                 `-->     TPC-TRD-TOF {:>10}", mCountersMC[3][i]);
      LOGP(info, "                 `-->     TPC-TRD     {:>10}", mCountersMC[4][i]);
      LOGP(info, "                 `-->     TPC-   -TOF {:>10}", mCountersMC[5][i]);
      LOGP(info, "                 `--> ITS-TPC         {:>10}", mCountersMC[6][i]);
      LOGP(info, "                 `--> ITS             {:>10}", mCountersMC[7][i]);
      LOGP(info, "                 `-->     TPC         {:>10}", mCountersMC[8][i]);
      LOGP(info, "                 `--> ???             {:>10}", mCountersMC[9][i]);
    }

    void printMC2()
    {
      for (int i{0}; i < static_cast<int>(Enum::NSIZE); ++i) {
        LOGP(info, "{} - CHECK: {}: {}", i, _names[i], mTotCounters[i]);
        LOGP(info, "                 `--> Track           & Track           {:>15}", "true V0s pairs");
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF & ITS-TPC-TRD-TOF {:>15}", mCounters[0][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF & ITS-TPC-   -TOF {:>15}", mCounters[1][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF & ITS-TPC-TRD     {:>15}", mCounters[2][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF &     TPC-TRD-TOF {:>15}", mCounters[3][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF &     TPC-TRD     {:>15}", mCounters[4][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF &     TPC-   -TOF {:>15}", mCounters[5][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF & ITS-TPC         {:>15}", mCounters[6][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF & ITS             {:>15}", mCounters[7][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF &     TPC         {:>15}", mCounters[8][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD-TOF & ???             {:>15}", mCounters[9][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF & ITS-TPC-   -TOF {:>15}", mCounters[10][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF & ITS-TPC-TRD     {:>15}", mCounters[11][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF &     TPC-TRD-TOF {:>15}", mCounters[12][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF &     TPC-TRD     {:>15}", mCounters[13][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF &     TPC-   -TOF {:>15}", mCounters[14][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF & ITS-TPC         {:>15}", mCounters[15][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF & ITS             {:>15}", mCounters[16][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF &     TPC         {:>15}", mCounters[17][i]);
        LOGP(info, "                 `--> ITS-TPC-   -TOF & ???             {:>15}", mCounters[18][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     & ITS-TPC-TRD     {:>15}", mCounters[19][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     &     TPC-TRD-TOF {:>15}", mCounters[20][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     &     TPC-TRD     {:>15}", mCounters[21][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     &     TPC-   -TOF {:>15}", mCounters[22][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     & ITS-TPC         {:>15}", mCounters[23][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     & ITS             {:>15}", mCounters[24][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     &     TPC         {:>15}", mCounters[25][i]);
        LOGP(info, "                 `--> ITS-TPC-TRD     & ???             {:>15}", mCounters[26][i]);
        LOGP(info, "                 `-->     TPC-TRD-TOF &     TPC-TRD-TOF {:>15}", mCounters[27][i]);
        LOGP(info, "                 `-->     TPC-TRD-TOF &     TPC-TRD     {:>15}", mCounters[28][i]);
        LOGP(info, "                 `-->     TPC-TRD-TOF &     TPC-   -TOF {:>15}", mCounters[29][i]);
        LOGP(info, "                 `-->     TPC-TRD-TOF & ITS-TPC         {:>15}", mCounters[30][i]);
        LOGP(info, "                 `-->     TPC-TRD-TOF & ITS             {:>15}", mCounters[31][i]);
        LOGP(info, "                 `-->     TPC-TRD-TOF &     TPC         {:>15}", mCounters[32][i]);
        LOGP(info, "                 `-->     TPC-TRD-TOF & ???             {:>15}", mCounters[33][i]);
        LOGP(info, "                 `-->     TPC-TRD     &     TPC-TRD     {:>15}", mCounters[34][i]);
        LOGP(info, "                 `-->     TPC-TRD     &     TPC-   -TOF {:>15}", mCounters[35][i]);
        LOGP(info, "                 `-->     TPC-TRD     & ITS-TPC         {:>15}", mCounters[36][i]);
        LOGP(info, "                 `-->     TPC-TRD     & ITS             {:>15}", mCounters[37][i]);
        LOGP(info, "                 `-->     TPC-TRD     &     TPC         {:>15}", mCounters[38][i]);
        LOGP(info, "                 `-->     TPC-TRD     & ???             {:>15}", mCounters[39][i]);
        LOGP(info, "                 `-->     TPC-   -TOF &     TPC-   -TOF {:>15}", mCounters[40][i]);
        LOGP(info, "                 `-->     TPC-   -TOF & ITS-TPC         {:>15}", mCounters[41][i]);
        LOGP(info, "                 `-->     TPC-   -TOF & ITS             {:>15}", mCounters[42][i]);
        LOGP(info, "                 `-->     TPC-   -TOF &     TPC         {:>15}", mCounters[43][i]);
        LOGP(info, "                 `-->     TPC-   -TOF & ???             {:>15}", mCounters[44][i]);
        LOGP(info, "                 `--> ITS-TPC         & ITS-TPC         {:>15}", mCounters[45][i]);
        LOGP(info, "                 `--> ITS-TPC         & ITS             {:>15}", mCounters[46][i]);
        LOGP(info, "                 `--> ITS-TPC         &     TPC         {:>15}", mCounters[47][i]);
        LOGP(info, "                 `--> ITS-TPC         & ???             {:>15}", mCounters[48][i]);
        LOGP(info, "                 `--> ITS             & ITS             {:>15}", mCounters[49][i]);
        LOGP(info, "                 `--> ITS             &     TPC         {:>15}", mCounters[50][i]);
        LOGP(info, "                 `--> ITS             & ???             {:>15}", mCounters[51][i]);
        LOGP(info, "                 `-->     TPC         &     TPC         {:>15}", mCounters[52][i]);
        LOGP(info, "                 `-->     TPC         & ???             {:>15}", mCounters[53][i]);
        LOGP(info, "                 `--> ???             & ???             {:>15}", mCounters[54][i]);
      }
    }
  };

  enum class CHECKV0 : unsigned int {
    FPROCESS = 0,
    MINR2TOMEANVERTEX,
    REJCAUSALITY,
    PROPVTX,
    REJPT2,
    REJTGL,
    REJCPA,
    CALLED,
    NSIZE,
  };
  static constexpr std::array<std::string_view, static_cast<size_t>(CHECKV0::NSIZE)> checkV0Names{
    "Fitter Processing",
    "Min R2 to mean Vertex",
    "Rejection Causality",
    "Propagating to vertex",
    "Rejection Pt2",
    "Rejection TgL",
    "Rejection Cos Pointing Angle",
    "#CALLED",
  };
  static constexpr std::string_view checkV0TreeName{"checkV0"};
  Counter_t<CHECKV0> mCounterV0{checkV0TreeName, checkV0Names};

  enum class BUILDT2V : unsigned int {
    NOTLOADED = 0,
    TPCTRACK,
    TPCEXCLUDE,
    TPCSPROCESS,
    TPCFPROCESS,
    AMBIGIOUS,
    ACCOUNT,
    REJECTED,
    HEAVY,
    BACCEPT,
    NACCEPT,
    NACCEPTAMBI,
    AACCEPT,
    CALLED,
    NSIZE,
  };
  static constexpr std::array<std::string_view, static_cast<size_t>(BUILDT2V::NSIZE)> buildT2VNames{
    "Track source not loaded",
    "TPC track",
    "Excluded TPC track",
    "TPC track successfully processed",
    "TPC track constrained processed",
    "Ambigious track",
    "Ambigious: Already accounted track (latter one)",
    "Ambigious: Already rejected track",
    "Heavy ionising particles",
    "Before acceptTrack",
    "Not acceptTrack",
    "Not acceptTrack Ambigious",
    "After acceptTrack",
    "#CALLED",
  };
  static constexpr std::string_view buildT2VTreeName{"buildT2V"};
  Counter_t<BUILDT2V> mCounterBuildT2V{buildT2VTreeName, buildT2VNames};

  enum class MCGEN : unsigned int {
    GEN = 0,
    NSIZE,
  };
  static constexpr std::array<std::string_view, static_cast<size_t>(MCGEN::NSIZE)> mcNames{
    "V0s generated",
  };
  static constexpr std::string_view mcTreeName{"alksjdlk"};
  Counter_t<MCGEN> mCounterMC{mcTreeName, mcNames};

  enum class Reco : unsigned int {
    WITHOUT = 0,
    WITH,
    NSIZE,
  };
  static constexpr std::array<std::string_view, static_cast<size_t>(Reco::NSIZE)> recoNames{
    "Reconstructed true V0s findable without timing information (e.g., in RecoContainer)",
    "Reconstructed true V0s findable with timing information (e.g., in TrackPool)",
  };
  static constexpr std::string_view recoTreeName{"recoTracks"};
  Counter_t<Reco> mCounterReco{recoTreeName, recoNames};

  enum class V0Found : unsigned int {
    FOUND,
    NSIZE,
  };
  static constexpr std::array<std::string_view, static_cast<size_t>(V0Found::NSIZE)> v0Names{
    "Found V0s"};
  static constexpr std::string_view v0TreeName{"v0Found"};
  Counter_t<V0Found> mCounterV0Found{v0TreeName, v0Names};

  static constexpr auto v0Type = kGamma;
};
} // namespace vertexing
} // namespace o2

#endif
