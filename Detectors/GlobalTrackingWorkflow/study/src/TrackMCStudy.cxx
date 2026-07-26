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

#include <vector>
#include <TStopwatch.h>
#include <TMCProcess.h>
#include <TPDGCode.h>
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DataFormatsGlobalTracking/RecoContainerCreateTracksVariadic.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "TPCCalibration/VDriftHelper.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"
#include "ITStracking/IOUtils.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "ITSBase/GeometryTGeo.h"
#include "SimulationDataFormat/MCEventLabel.h"
#include "SimulationDataFormat/MCUtils.h"
#include "SimulationDataFormat/O2DatabasePDG.h"
#include "SimulationDataFormat/TrackReference.h"
#include "CommonDataFormat/BunchFilling.h"
#include "CommonUtils/NameConf.h"
#include "DataFormatsFT0/RecPoints.h"
#include "DataFormatsITSMFT/TrkClusRef.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/CCDBParamSpec.h"
#include "FT0Reconstruction/InteractionTag.h"
#include "DataFormatsITSMFT/DPLAlpideParam.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "GlobalTrackingStudy/TrackMCStudy.h"
#include "GlobalTrackingStudy/TrackMCStudyConfig.h"
#include "GlobalTrackingStudy/TrackMCStudyTypes.h"
#include "GlobalTracking/MatchTPCITSParams.h"
#include "TPCBase/ParameterElectronics.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/PrimaryVertexExt.h"
#include "DataFormatsFT0/RecPoints.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "ReconstructionDataFormats/VtxTrackRef.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Steer/MCKinematicsReader.h"
#include "DCAFitter/DCAFitterN.h"
#include "DetectorsVertexing/SVertexerParams.h"
#include "DetectorsVertexing/SVertexer.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include "GPUO2InterfaceRefit.h"
#include "GPUParam.h"
#include "GPUParam.inc"
#include "MathUtils/fit.h"
#include "MathUtils/Primitive2D.h"
#include "TPCFastTransformPOD.h"
#include <TRandom.h>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <utility>
#include <gsl/span>

// workflow to study relation of reco tracks to MCTruth
// o2-trackmc-study-workflow --device-verbosity 3 -b --run

namespace o2::trackstudy
{

using namespace o2::framework;
using DetID = o2::detectors::DetID;
using DataRequest = o2::globaltracking::DataRequest;

using PVertex = o2::dataformats::PrimaryVertex;
using V2TRef = o2::dataformats::VtxTrackRef;
using VTIndex = o2::dataformats::VtxTrackIndex;
using VTIndexV = std::pair<int, o2::dataformats::VtxTrackIndex>;
using GTrackID = o2::dataformats::GlobalTrackID;
using TBracket = o2::math_utils::Bracketf_t;

using timeEst = o2::dataformats::TimeStampWithError<float, float>;

class TrackMCStudy final : public Task
{
 public:
  TrackMCStudy(std::shared_ptr<DataRequest> dr, std::shared_ptr<o2::base::GRPGeomRequest> gr, GTrackID::mask_t src, bool checkSV, bool useCCDBParams, bool enableCasc, bool enable3body)
    : mDataRequest(dr), mGGCCDBRequest(gr), mTracksSrc(src), mCheckSV(checkSV), mUseCCDBParams(useCCDBParams), mEnableCascades(enableCasc), mEnable3BodyDecays(enable3body) {}
  ~TrackMCStudy() final = default;
  void init(InitContext& ic) final;
  void run(ProcessingContext& pc) final;
  void endOfStream(EndOfStreamContext& ec) final;
  void finaliseCCDB(ConcreteDataMatcher& matcher, void* obj) final;
  void process(const o2::globaltracking::RecoContainer& recoData);

 private:
  void processTPCTrackRefs();
  void processITSTracks(const o2::globaltracking::RecoContainer& recoData);
  void loadTPCOccMap(const o2::globaltracking::RecoContainer& recoData);
  void fillMCClusterInfo(const o2::globaltracking::RecoContainer& recoData);
  void prepareITSData(const o2::globaltracking::RecoContainer& recoData);
  bool processMCParticle(int src, int ev, int trid);
  bool addMCParticle(const MCTrack& mctr, const o2::MCCompLabel& lb, TParticlePDG* pPDG = nullptr);
  bool acceptMCCharged(const MCTrack& tr, const o2::MCCompLabel& lb, int followDec = -1);
  bool propagateToRefX(o2::track::TrackParCov& trcTPC, o2::track::TrackParCov& trcITS);
  bool refitV0(int i, o2::dataformats::V0& v0, const o2::globaltracking::RecoContainer& recoData);
  void checkSVertexer(const o2::globaltracking::RecoContainer& recoData);
  void evalTPCPhotonTune(const o2::globaltracking::RecoContainer& recoData);
  void processRecSVs(const o2::globaltracking::RecoContainer& recoData);
  void collectGammaConversions(int src, int ev);
  void processGammaConversions(const o2::globaltracking::RecoContainer& recoData);
  void updateTimeDependentParams(ProcessingContext& pc);
  float getDCAYCut(float pt) const;

  const std::vector<o2::MCTrack>* mCurrMCTracks = nullptr;
  TVector3 mCurrMCVertex;
  o2::tpc::VDriftHelper mTPCVDriftHelper{};
  const o2::gpu::TPCFastTransformPOD* mTPCCorrMaps{nullptr};
  std::shared_ptr<DataRequest> mDataRequest;
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  std::unique_ptr<o2::utils::TreeStreamRedirector> mDBGOut;
  std::vector<float> mTBinClOcc;                            ///< TPC occupancy histo: i-th entry is the integrated occupancy for ~1 orbit starting from the TB = i*mNTPCOccBinLength
  std::vector<float> mTBinClOccHist;                        //< original occupancy
  std::vector<long> mIntBC;                                 ///< interaction global BC wrt TF start
  std::vector<float> mTPCOcc;                               ///< TPC occupancy for this interaction time
  std::vector<int> mITSOcc;                                 //< N ITS clusters in the ROF containing collision
  std::vector<o2::BaseCluster<float>> mITSClustersArray;    ///< ITS clusters created in run() method from compact clusters
  const o2::itsmft::TopologyDictionary* mITSDict = nullptr; ///< cluster patterns dictionary

  bool mCheckSV = false;          //< check SV binding (apart from prongs availability)
  bool mUseCCDBParams = false;    //< fetch SVertexerParams from CCDB, as the reconstruction does
  bool mEnableCascades = true;    //< must match the reconstruction: it changes the V0 selection itself
  bool mEnable3BodyDecays = true; //< idem, notably the cosPA floor applied to V0 candidates
  bool mRecProcStage = false;     //< flag that the MC particle was added only at the stage of reco tracks processing
  int mNTPCOccBinLength = 0;      ///< TPC occ. histo bin length in TBs
  float mNTPCOccBinLengthInv = -1.f;
  int mVerbose = 0;
  float mITSTimeBiasMUS = 0.f;
  float mITSROFrameLengthMUS = 0.f; ///< ITS RO frame in mus
  float mTPCTBinMUS = 0.;           ///< TPC time bin duration in microseconds

  int mNCheckDecays = 0;

  GTrackID::mask_t mTracksSrc{};
  o2::steer::MCKinematicsReader mcReader; // reader of MC information
  std::vector<int> mITSROF;
  std::vector<TBracket> mITSROFBracket;
  std::vector<o2::MCCompLabel> mDecProdLblPool; // labels of decay products to watch, added to MC map
  std::vector<MCVertex> mMCVtVec{};

  struct DecayRef {
    o2::MCCompLabel mother{};
    o2::track::TrackPar parent{};
    int pdg = 0;
    int daughterFirst = -1;
    int daughterLast = -1;
    int foundSVID = -1;
    SVCheck svCheck{};
  };
  std::vector<std::vector<DecayRef>> mDecaysMaps; // for every parent particle to watch, store its label and entries of 1st/last decay product labels in mDecProdLblPool
  std::vector<GammaConvInfo> mGammaConversions;
  std::unordered_set<o2::MCCompLabel> mGammaDaughterLabels;
  std::unordered_map<VTIndex, std::vector<int>> mGammaTrackPVRefs;
  std::unordered_map<o2::MCCompLabel, TrackFamily> mSelMCTracks;
  std::unordered_map<o2::MCCompLabel, std::pair<int, int>> mSelTRefIdx;
  std::vector<o2::track::TrackPar> mSelTRefs;
  o2::vertexing::DCAFitterN<2> mFitterV0;
  o2::vertexing::SVertexer mSVertexer{};
  static constexpr float MaxSnp = 0.9; // max snp of ITS or TPC track at xRef to be matched
};

void TrackMCStudy::init(InitContext& ic)
{
  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);
  mcReader.initFromDigitContext("collisioncontext.root");

  mDBGOut = std::make_unique<o2::utils::TreeStreamRedirector>("trackMCStudy.root", "recreate");
  mVerbose = ic.options().get<int>("device-verbosity");

  const auto& params = o2::trackstudy::TrackMCStudyConfig::Instance();
  params.printKeyValues(true, true);
  mNCheckDecays = (int)params.decayPDG.size();
  mDecaysMaps.resize(mNCheckDecays);
  if (params.storeTPCPhotonTune) {
    // the signal sample is defined by mGammaDaughterLabels, which only collectGammaConversions fills
    if (!params.storeGammaConversions) {
      LOGP(fatal, "trmcconf.storeTPCPhotonTune requires trmcconf.storeGammaConversions, otherwise every "
                  "recorded TPC track would be classified as background");
    }
    if (!mCheckSV) {
      LOGP(fatal, "trmcconf.storeTPCPhotonTune needs the SVertexer seeds pool, do not pass --ignore-sv-check");
    }
  }
  if (mCheckSV) {
    mSVertexer.setEnableCascades(mEnableCascades);
    mSVertexer.setEnable3BodyDecays(mEnable3BodyDecays);
    mSVertexer.setNThreads(1);
    mSVertexer.setUseMC(false);
    mSVertexer.setCollectSeedRejections(true);
  }
}

void TrackMCStudy::run(ProcessingContext& pc)
{
  o2::globaltracking::RecoContainer recoData;
  for (int i = 0; i < mNCheckDecays; i++) {
    mDecaysMaps[i].clear();
  }
  mDecProdLblPool.clear();
  mGammaConversions.clear();
  mGammaDaughterLabels.clear();
  mGammaTrackPVRefs.clear();
  mMCVtVec.clear();
  mCurrMCTracks = nullptr;

  recoData.collectData(pc, *mDataRequest.get()); // select tracks of needed type, with minimal cuts, the real selected will be done in the vertexer
  updateTimeDependentParams(pc);                 // Make sure this is called after recoData.collectData, which may load some conditions
  mRecProcStage = false;
  process(recoData);
}

void TrackMCStudy::updateTimeDependentParams(ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  mTPCVDriftHelper.extractCCDBInputs(pc);
  auto const& raw = pc.inputs().get<const char*>("corrMap");
  mTPCCorrMaps = &o2::gpu::TPCFastTransformPOD::get(raw);
  static bool initOnceDone = false;
  if (!initOnceDone) { // this params need to be queried only once
    initOnceDone = true;
    if (mCheckSV && mUseCCDBParams) {
      pc.inputs().get<o2::vertexing::SVertexerParams*>("SVParam");
      pc.inputs().get<o2::dataformats::MeanVertexObject*>("meanvtx");
    }
    const auto& alpParamsITS = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>::Instance();
    mITSROFrameLengthMUS = o2::base::GRPGeomHelper::instance().getGRPECS()->isDetContinuousReadOut(o2::detectors::DetID::ITS) ? alpParamsITS.roFrameLengthInBC * o2::constants::lhc::LHCBunchSpacingMUS : alpParamsITS.roFrameLengthTrig * 1.e-3;
    LOGP(info, "VertexTrackMatcher ITSROFrameLengthMUS:{}", mITSROFrameLengthMUS);

    auto& elParam = o2::tpc::ParameterElectronics::Instance();
    mTPCTBinMUS = elParam.ZbinWidth;
    o2::its::GeometryTGeo::Instance()->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2GRot) | o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L));
    if (mCheckSV) {
      const auto& svparam = o2::vertexing::SVertexerParams::Instance();
      mFitterV0.setBz(o2::base::Propagator::Instance()->getNominalBz());
      mFitterV0.setUseAbsDCA(svparam.useAbsDCA);
      mFitterV0.setPropagateToPCA(false);
      mFitterV0.setMaxR(svparam.maxRIni);
      mFitterV0.setMinParamChange(svparam.minParamChange);
      mFitterV0.setMinRelChi2Change(svparam.minRelChi2Change);
      mFitterV0.setMaxDZIni(svparam.maxDZIni);
      mFitterV0.setMaxDXYIni(svparam.maxDXYIni);
      mFitterV0.setMaxChi2(svparam.maxChi2);
      mFitterV0.setMatCorrType(o2::base::Propagator::MatCorrType(svparam.matCorr));
      mFitterV0.setUsePropagator(svparam.usePropagator);
      mFitterV0.setRefitWithMatCorr(svparam.refitWithMatCorr);
      mFitterV0.setMaxStep(svparam.maxStep);
      mFitterV0.setMaxSnp(svparam.maxSnp);
      mFitterV0.setMinXSeed(svparam.minXSeed);
      mSVertexer.init();
    }
  }
  if (mCheckSV) {
    mSVertexer.setTPCCorrMaps(mTPCCorrMaps);
    if (mTPCVDriftHelper.isUpdated()) {
      mSVertexer.setTPCVDrift(mTPCVDriftHelper.getVDriftObject());
    }
  }
}

void TrackMCStudy::process(const o2::globaltracking::RecoContainer& recoData)
{
  constexpr float SQRT12Inv = 0.288675f;
  const auto& params = o2::trackstudy::TrackMCStudyConfig::Instance();
  auto pvvec = recoData.getPrimaryVertices();
  auto pvvecLbl = recoData.getPrimaryVertexMCLabels();
  auto trackIndex = recoData.getPrimaryVertexMatchedTracks(); // Global ID's for associated tracks
  auto vtxRefs = recoData.getPrimaryVertexMatchedTrackRefs(); // references from vertex to these track IDs
  auto prop = o2::base::Propagator::Instance();
  int nv = vtxRefs.size();
  float vdriftTB = mTPCVDriftHelper.getVDriftObject().getVDrift() * o2::tpc::ParameterElectronics::Instance().ZbinWidth;                                                     // VDrift expressed in cm/TimeBin
  float itsBias = (0.5 * mITSROFrameLengthMUS) + o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>::Instance().roFrameBiasInBC * o2::constants::lhc::LHCBunchSpacingMUS; // ITS time is supplied in \mus as beginning of ROF

  prepareITSData(recoData);
  loadTPCOccMap(recoData);
  auto getITSPatt = [&](GTrackID gid, uint8_t& ncl) {
    int8_t patt = 0;
    if (gid.getSource() == VTIndex::ITSAB) {
      const auto& itsTrf = recoData.getITSABRefs()[gid];
      ncl = itsTrf.getNClusters();
      for (int il = 0; il < 7; il++) {
        if (itsTrf.hasHitOnLayer(il)) {
          patt |= 0x1 << il;
        }
      }
      patt |= 0x1 << 7;
    } else {
      const auto& itsTr = recoData.getITSTrack(gid);
      for (int il = 0; il < 7; il++) {
        if (itsTr.hasHitOnLayer(il)) {
          patt |= 0x1 << il;
          ncl++;
        }
      }
    }
    return patt;
  };

  auto fillTPCClusterInfo = [&recoData](const o2::tpc::TrackTPC& trc, RecTrack& tref) {
    if (recoData.inputsTPCclusters) {
      uint8_t clSect = 0, clRow = 0, lowestR = -1;
      uint32_t clIdx = 0;
      const auto clRefs = recoData.getTPCTracksClusterRefs();
      const auto tpcClusAcc = recoData.getTPCClusters();
      const auto shMap = recoData.clusterShMapTPC;
      for (int ic = 0; ic < trc.getNClusterReferences(); ic++) { // outside -> inside ordering, but on the sector boundaries backward jumps are possible
        trc.getClusterReference(clRefs, ic, clSect, clRow, clIdx);
        if (clRow < lowestR) {
          tref.rowCountTPC++;
          lowestR = clRow;
        }
        unsigned int absoluteIndex = tpcClusAcc.clusterOffset[clSect][clRow] + clIdx;
        if (shMap[absoluteIndex] & o2::gpu::GPUTPCGMMergedTrackHit::flagShared) {
          tref.nClTPCShared++;
        }
      }
      tref.lowestPadRow = lowestR;
      const auto& clus = tpcClusAcc.clusters[clSect][clRow][clIdx];
      int padFromEdge = int(clus.getPad()), npads = o2::gpu::GPUTPCGeometry::NPads(clRow);
      if (padFromEdge > npads / 2) {
        padFromEdge = npads - 1 - padFromEdge;
      }
      tref.padFromEdge = uint8_t(padFromEdge);
      trc.getClusterReference(clRefs, 0, clSect, clRow, clIdx);
      tref.rowMaxTPC = clRow;
    }
  };

  auto flagTPCClusters = [&recoData](const o2::tpc::TrackTPC& trc, o2::MCCompLabel lbTrc) {
    if (recoData.inputsTPCclusters) {
      const auto clRefs = recoData.getTPCTracksClusterRefs();
      const auto* TPCClMClab = recoData.inputsTPCclusters->clusterIndex.clustersMCTruth;
      const auto& TPCClusterIdxStruct = recoData.inputsTPCclusters->clusterIndex;
      for (int ic = 0; ic < trc.getNClusterReferences(); ic++) {
        uint8_t clSect = 0, clRow = 0;
        uint32_t clIdx = 0;
        trc.getClusterReference(clRefs, ic, clSect, clRow, clIdx);
        auto labels = TPCClMClab->getLabels(clIdx + TPCClusterIdxStruct.clusterOffset[clSect][clRow]);
        for (auto& lbl : labels) {
          if (lbl == lbTrc) {
            const_cast<o2::MCCompLabel&>(lbl).setFakeFlag(true); // actually, in this way we are flagging that this cluster was correctly attached
            break;
          }
        }
      }
    }
  };

  {
    const auto* digconst = mcReader.getDigitizationContext();
    const auto& mcEvRecords = digconst->getEventRecords(false);
    int ITSTimeBias = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>::Instance().roFrameBiasInBC;
    int ITSROFLen = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>::Instance().roFrameLengthInBC;
    unsigned int rofCount = 0;
    const auto ITSClusROFRec = recoData.getITSClustersROFRecords();
    for (const auto& mcIR : mcEvRecords) {
      long tbc = mcIR.differenceInBC(recoData.startIR);
      auto& mcVtx = mMCVtVec.emplace_back();
      mcVtx.ts = (tbc * o2::constants::lhc::LHCBunchSpacingMUS) + mcIR.getTimeOffsetWrtBC() * 1e-3;
      mcVtx.ID = mIntBC.size();
      mIntBC.push_back(tbc);
      int occBin = tbc / 8 * mNTPCOccBinLengthInv;
      mTPCOcc.push_back(occBin < 0 ? mTBinClOcc[0] : (occBin >= mTBinClOcc.size() ? mTBinClOcc.back() : mTBinClOcc[occBin]));
      // fill ITS occupancy
      long gbc = mcIR.toLong();
      while (rofCount < ITSClusROFRec.size()) {
        long rofbcMin = ITSClusROFRec[rofCount].getBCData().toLong() + ITSTimeBias, rofbcMax = rofbcMin + ITSROFLen;
        if (gbc < rofbcMin) { // IRs and ROFs are sorted, so this IR is prior of all ROFs
          mITSOcc.push_back(0);
        } else if (gbc < rofbcMax) {
          mITSOcc.push_back(ITSClusROFRec[rofCount].getNEntries());
        } else {
          rofCount++; // test next ROF
          continue;
        }
        break;
      }
      if (mNTPCOccBinLengthInv > 0.f) {
        mcVtx.occTPCV.resize(params.nOccBinsDrift);
        int grp = TMath::Max(1, TMath::Nint(params.nTBPerOccBin * mNTPCOccBinLengthInv));
        for (int ib = 0; ib < params.nOccBinsDrift; ib++) {
          float smb = 0;
          int tbs = occBin + TMath::Nint(ib * params.nTBPerOccBin * mNTPCOccBinLengthInv);
          for (int ig = 0; ig < grp; ig++) {
            if (tbs >= 0 && tbs < int(mTBinClOccHist.size())) {
              smb += mTBinClOccHist[tbs];
            }
            tbs++;
          }
          mcVtx.occTPCV[ib] = smb;
        }
      }
      if (rofCount >= ITSClusROFRec.size()) {
        mITSOcc.push_back(0); // IR after the last ROF
      }
    }
  }
  // collect interesting MC particle (tracks and parents)
  int curSrcMC = 0, curEvMC = 0;
  for (curSrcMC = 0; curSrcMC < (int)mcReader.getNSources(); curSrcMC++) {
    if (mVerbose > 1) {
      LOGP(info, "Source {}", curSrcMC);
    }
    int nev = mcReader.getNEvents(curSrcMC);
    bool okAccVtx = true;
    if (nev != (int)mMCVtVec.size()) {
      LOGP(debug, "source {} has {} events while {} MC vertices were booked", curSrcMC, nev, mMCVtVec.size());
      okAccVtx = false;
      if (nev > (int)mMCVtVec.size()) { // QED
        continue;
      }
    }
    for (curEvMC = 0; curEvMC < nev; curEvMC++) {
      if (mVerbose > 1) {
        LOGP(info, "Event {}", curEvMC);
      }
      mCurrMCTracks = &mcReader.getTracks(curSrcMC, curEvMC);
      const_cast<o2::dataformats::MCEventHeader&>(mcReader.getMCEventHeader(curSrcMC, curEvMC)).GetVertex(mCurrMCVertex);
      if (okAccVtx) {
        auto& pos = mMCVtVec[curEvMC].pos;
        if (pos[2] < -999) {
          pos[0] = mCurrMCVertex.X();
          pos[1] = mCurrMCVertex.Y();
          pos[2] = mCurrMCVertex.Z();
        }
      }
      if (params.storeGammaConversions) {
        collectGammaConversions(curSrcMC, curEvMC);
      }
      for (int itr = 0; itr < mCurrMCTracks->size(); itr++) {
        processMCParticle(curSrcMC, curEvMC, itr);
      }
    }
  }
  if (mVerbose > 0) {
    for (int id = 0; id < mNCheckDecays; id++) {
      LOGP(info, "Decay PDG={} : {} entries", params.decayPDG[id], mDecaysMaps[id].size());
    }
  }

  // add reconstruction info to MC particles. If MC particle was not selected before but was reconstrected, account MC info
  mRecProcStage = true; // MC particles accepted only at this stage will be flagged
  for (int iv = 0; iv < nv; iv++) {
    if (mVerbose > 1) {
      LOGP(info, "processing PV {} of {}", iv, nv);
    }
    o2::MCEventLabel pvLbl;
    int pvID = -1;
    if (iv < (int)pvvecLbl.size()) {
      pvLbl = pvvecLbl[iv];
      pvID = iv;
      if (pvLbl.isSet() && pvLbl.getEventID() < mMCVtVec.size()) {
        mMCVtVec[pvLbl.getEventID()].recVtx.emplace_back(RecPV{.pv = pvvec[iv], .mcEvLbl = pvLbl});
      }
    }
    const auto& vtref = vtxRefs[iv];
    for (int is = GTrackID::NSources; is--;) {
      DetID::mask_t dm = GTrackID::getSourceDetectorsMask(is);
      if (!mTracksSrc[is] || !recoData.isTrackSourceLoaded(is) || !(dm[DetID::ITS] || dm[DetID::TPC])) {
        continue;
      }
      int idMin = vtref.getFirstEntryOfSource(is), idMax = idMin + vtref.getEntriesOfSource(is);
      for (int i = idMin; i < idMax; i++) {
        auto vid = trackIndex[i];
        const auto& trc = recoData.getTrackParam(vid);
        auto lbl = recoData.getTrackMCLabel(vid);
        bool isGammaDaughter = false;
        if (lbl.isValid()) {
          lbl.setFakeFlag(false);
          isGammaDaughter = mGammaDaughterLabels.contains(lbl);
        }
        if (!isGammaDaughter && (trc.getPt() < params.minPt || std::abs(trc.getTgl()) > params.maxTgl)) {
          continue;
        }
        if (lbl.isValid()) {
          if (isGammaDaughter && iv < (int)pvvec.size()) {
            auto& pvIDs = mGammaTrackPVRefs[vid];
            if (std::find(pvIDs.begin(), pvIDs.end(), iv) == pvIDs.end()) {
              pvIDs.push_back(iv);
            }
          }
          auto entry = mSelMCTracks.find(lbl);
          if (entry == mSelMCTracks.end()) { // add the track which was not added during MC scan
            if (lbl.getSourceID() != curSrcMC || lbl.getEventID() != curEvMC) {
              curSrcMC = lbl.getSourceID();
              curEvMC = lbl.getEventID();
              mCurrMCTracks = &mcReader.getTracks(curSrcMC, curEvMC);
              const_cast<o2::dataformats::MCEventHeader&>(mcReader.getMCEventHeader(curSrcMC, curEvMC)).GetVertex(mCurrMCVertex);
            }
            if (!acceptMCCharged((*mCurrMCTracks)[lbl.getTrackID()], lbl)) {
              continue;
            }
            entry = mSelMCTracks.find(lbl);
          }
          auto& trackFamily = entry->second;
          if (vid.isAmbiguous()) { // do not repeat ambiguous tracks
            if (trackFamily.contains(vid)) {
              continue;
            }
          }
          auto& trf = trackFamily.recTracks.emplace_back();
          trf.gid = vid; //  account(iv, vid);
          trf.pvID = pvID;
          trf.pvLabel = pvLbl;
          while (dm[DetID::ITS] && dm[DetID::TPC]) { // this track should have both ITS and TPC parts, if ITS was mismatched, fill it to its proper MC track slot
            auto gidSet = recoData.getSingleDetectorRefs(vid);
            if (!gidSet[GTrackID::ITS].isSourceSet()) {
              break; // AB track, nothing to check
            }
            auto lblITS = recoData.getTrackMCLabel(gidSet[GTrackID::ITS]);
            if (lblITS == trackFamily.mcTrackInfo.label) {
              break; // correct match, no need for special treatment
            }
            const auto& trcITSF = recoData.getTrackParam(gidSet[GTrackID::ITS]);
            if (trcITSF.getPt() < params.minPt || std::abs(trcITSF.getTgl()) > params.maxTgl) {
              break; // ignore this track
            }
            auto entryOfFake = mSelMCTracks.find(lblITS);
            if (entryOfFake == mSelMCTracks.end()) { // this MC track was not selected
              break;
            }
            auto& trackFamilyOfFake = entryOfFake->second;
            auto& trfOfFake = trackFamilyOfFake.recTracks.emplace_back();
            trfOfFake.gid = gidSet[GTrackID::ITS]; //  account(iv, vid);
            break;
          }
          if (mVerbose > 1) {
            LOGP(info, "Matched rec track {} to MC track {}", vid.asString(), entry->first.asString());
          }
        } else {
          continue;
        }
      }
    }
  }

  LOGP(info, "collected {} MC tracks", mSelMCTracks.size());
  if (params.minTPCRefsToExtractClRes > 0 || params.storeTPCTrackRefs) { // prepare MC trackrefs for TPC
    processTPCTrackRefs();
  }

  int mcnt = 0;
  for (auto& entry : mSelMCTracks) {
    auto& trackFam = entry.second;
    auto& tracks = trackFam.recTracks;
    mcnt++;
    if (tracks.empty()) {
      continue;
    }
    if (mVerbose > 1) {
      LOGP(info, "Processing MC track#{} {} -> {} reconstructed tracks", mcnt - 1, entry.first.asString(), tracks.size());
    }
    // sort according to the gid complexity (in principle, should be already sorted due to the backwards loop over NSources above
    std::sort(tracks.begin(), tracks.end(), [](const RecTrack& lhs, const RecTrack& rhs) {
      const auto mskL = lhs.gid.getSourceDetectorsMask();
      const auto mskR = rhs.gid.getSourceDetectorsMask();
      bool itstpcL = mskL[DetID::ITS] && mskL[DetID::TPC], itstpcR = mskR[DetID::ITS] && mskR[DetID::TPC];
      if (itstpcL && !itstpcR) { // to avoid TPC/TRD or TPC/TOF shadowing ITS/TPC
        return true;
      }
      return lhs.gid.getSource() > rhs.gid.getSource();
    });
    if (params.storeTPCTrackRefs) {
      auto rft = mSelTRefIdx.find(entry.first);
      if (rft != mSelTRefIdx.end()) {
        auto rfent = rft->second;
        for (int irf = rfent.first; irf < rfent.second; irf++) {
          trackFam.mcTrackInfo.trackRefsTPC.push_back(mSelTRefs[irf]);
        }
      }
    }
    // fill track params
    int tcnt = 0;
    for (auto& tref : tracks) {
      if (tref.gid.isSourceSet()) {
        auto gidSet = recoData.getSingleDetectorRefs(tref.gid);
        tref.track = recoData.getTrackParam(tref.gid);
        if (recoData.getTrackMCLabel(tref.gid).isFake()) {
          tref.flags |= RecTrack::FakeGLO;
        }
        auto msk = tref.gid.getSourceDetectorsMask();
        if (msk[DetID::ITS]) {
          if (gidSet[GTrackID::ITS].isSourceSet()) { // has ITS track rather than AB tracklet
            tref.pattITS = getITSPatt(gidSet[GTrackID::ITS], tref.nClITS);
            if (trackFam.entITS < 0) {
              trackFam.entITS = tcnt;
            }
            auto lblITS = recoData.getTrackMCLabel(gidSet[GTrackID::ITS]);
            if (lblITS.isFake()) {
              tref.flags |= RecTrack::FakeITS;
            }
            if (lblITS == trackFam.mcTrackInfo.label) {
              trackFam.entITSFound = tcnt;
            }
          } else { // AB ITS tracklet
            tref.pattITS = getITSPatt(gidSet[GTrackID::ITSAB], tref.nClITS);
            if (recoData.getTrackMCLabel(gidSet[GTrackID::ITSAB]).isFake()) {
              tref.flags |= RecTrack::FakeITS;
            }
          }
          if (msk[DetID::TPC]) {
            if (trackFam.entITSTPC < 0) { // has both ITS and TPC contribution
              trackFam.entITSTPC = tcnt;
            }
            if (recoData.getTrackMCLabel(gidSet[GTrackID::ITSTPC]).isFake()) {
              tref.flags |= RecTrack::FakeITSTPC;
            }

            if (msk[DetID::TRD]) {
              if (recoData.getTrackMCLabel(gidSet[GTrackID::ITSTPCTRD]).isFake()) {
                tref.flags |= RecTrack::FakeTRD;
              }
              if (msk[DetID::TOF]) {
                if (recoData.getTrackMCLabel(gidSet[GTrackID::ITSTPCTRDTOF]).isFake()) {
                  tref.flags |= RecTrack::FakeTOF;
                }
              }
            } else {
              if (msk[DetID::TOF]) {
                if (recoData.getTrackMCLabel(gidSet[GTrackID::ITSTPCTOF]).isFake()) {
                  tref.flags |= RecTrack::FakeTOF;
                }
              }
            }
          }
        }
        if (msk[DetID::TPC]) {
          const auto& trtpc = recoData.getTPCTrack(gidSet[GTrackID::TPC]);
          tref.nClTPC = trtpc.getNClusters();
          if (trtpc.hasBothSidesClusters()) {
            tref.flags |= RecTrack::HASACSides;
          }
          fillTPCClusterInfo(trtpc, tref);
          flagTPCClusters(trtpc, entry.first);
          if (trackFam.entTPC < 0) {
            trackFam.entTPC = tcnt;
            trackFam.tpcT0 = trtpc.getTime0();
          }
          if (recoData.getTrackMCLabel(gidSet[GTrackID::TPC]).isFake()) {
            tref.flags |= RecTrack::FakeTPC;
          }
          if (!msk[DetID::ITS]) {
            if (msk[DetID::TRD]) {
              if (recoData.getTrackMCLabel(gidSet[GTrackID::TPCTRD]).isFake()) {
                tref.flags |= RecTrack::FakeTRD;
              }
              if (msk[DetID::TOF]) {
                if (recoData.getTrackMCLabel(gidSet[GTrackID::TPCTRDTOF]).isFake()) {
                  tref.flags |= RecTrack::FakeTOF;
                }
              }
            } else {
              if (msk[DetID::TOF]) {
                if (recoData.getTrackMCLabel(gidSet[GTrackID::TPCTOF]).isFake()) {
                  tref.flags |= RecTrack::FakeTOF;
                }
              }
            }
          }
        }
        float ts = 0, terr = 0;
        if (tref.gid.getSource() != GTrackID::ITS) {
          recoData.getTrackTime(tref.gid, ts, terr);
          tref.ts = timeEst{ts, terr};
        } else {
          const auto& itsBra = mITSROFBracket[mITSROF[tref.gid.getIndex()]];
          tref.ts = timeEst{itsBra.mean(), itsBra.delta() * SQRT12Inv};
        }
      } else {
        LOGP(info, "Invalid entry {} of {} getTrackMCLabel {}", tcnt, tracks.size(), tref.gid.asString());
      }
      tcnt++;
    }
    if (trackFam.entITS > -1 && trackFam.entTPC > -1) { // ITS and TPC were found but matching failed
      auto vidITS = recoData.getITSContributorGID(tracks[trackFam.entITS].gid);
      auto vidTPC = recoData.getTPCContributorGID(tracks[trackFam.entTPC].gid);
      auto trcTPC = recoData.getTrackParam(vidTPC);
      auto trcITS = recoData.getTrackParamOut(vidITS);
      if (propagateToRefX(trcTPC, trcITS)) {
        trackFam.trackITSProp = trcITS;
        trackFam.trackTPCProp = trcTPC;
      } else {
        trackFam.trackITSProp.invalidate();
        trackFam.trackTPCProp.invalidate();
      }
    } else {
      trackFam.trackITSProp.invalidate();
      trackFam.trackTPCProp.invalidate();
    }
  }

  // SVertices (V0s)
  if (mCheckSV) {
    auto v0s = recoData.getV0sIdx();
    auto prpr = [](o2::trackstudy::TrackFamily& f) {
      std::string s;
      s += fmt::format(" par {} Ntpccl={} Nitscl={} ", f.mcTrackInfo.pdgParent, f.mcTrackInfo.nTPCCl, f.mcTrackInfo.nITSCl);
      for (auto& t : f.recTracks) {
        s += t.gid.asString();
        s += " ";
      }
      return s;
    };
    for (int svID = 0; svID < (int)v0s.size(); svID++) {
      const auto& v0idx = v0s[svID];
      int nOKProngs = 0, realMCSVID = -1;
      int8_t decTypeID = -1;
      for (int ipr = 0; ipr < v0idx.getNProngs(); ipr++) {
        auto mcl = recoData.getTrackMCLabel(v0idx.getProngID(ipr)); // was this MC particle selected?
        mcl.setFakeFlag(false);                                     // mSelMCTracks is keyed by labels without the fake flag
        auto itl = mSelMCTracks.find(mcl);
        if (itl == mSelMCTracks.end()) {
          nOKProngs = -1; // was not selected as interesting one, ignore
          break;
        }
        auto& trackFamily = itl->second;
        int decayParentIndex = trackFamily.mcTrackInfo.parentEntry;
        if (decayParentIndex < 0) { // does not come from decay
          break;
        }
        if (ipr == 0) {
          realMCSVID = decayParentIndex;
          decTypeID = trackFamily.mcTrackInfo.parentDecID;
          nOKProngs = 1;
          LOGP(debug, "Prong{} {} comes from {}/{}", ipr, prpr(trackFamily), decTypeID, realMCSVID);
          continue;
        }
        if (realMCSVID != decayParentIndex || decTypeID != trackFamily.mcTrackInfo.parentDecID) {
          break;
        }
        LOGP(debug, "Prong{} {} comes from {}/{}", ipr, prpr(trackFamily), decTypeID, realMCSVID);
        nOKProngs++;
      }
      if (nOKProngs == v0idx.getNProngs()) { // all prongs are from the decay of MC parent which deemed to be interesting, flag it
        LOGP(debug, "Decay {}/{} was found", decTypeID, realMCSVID);
        mDecaysMaps[decTypeID][realMCSVID].foundSVID = svID;
      }
    }
    // build the seeds pool once, both consumers below read the same prepared state. Calling
    // prepareSeeds twice would rebuild it and wipe the seed rejection map collected on the way.
    if (params.checkSVertexerCuts || params.storeTPCPhotonTune) {
      mSVertexer.prepareSeeds(recoData);
    }
    if (params.checkSVertexerCuts) {
      checkSVertexer(recoData);
    }
    if (params.storeTPCPhotonTune) {
      evalTPCPhotonTune(recoData);
    }
    if (params.storeRecSV) {
      processRecSVs(recoData);
    }
  }
  if (params.storeGammaConversions) {
    processGammaConversions(recoData);
  }

  // collect ITS/TPC cluster info for selected MC particles
  fillMCClusterInfo(recoData);

  // single tracks
  for (auto& entry : mSelMCTracks) {
    auto& trackFam = entry.second;
    (*mDBGOut) << "tracks" << "tr=" << trackFam << "\n";
  }

  // decays
  std::vector<TrackFamily> decFam;
  for (int id = 0; id < mNCheckDecays; id++) {
    std::string decTreeName = fmt::format("dec{}", params.decayPDG[id]);
    for (const auto& dec : mDecaysMaps[id]) {
      decFam.clear();
      bool skip = false;
      for (int idd = dec.daughterFirst; idd <= dec.daughterLast; idd++) {
        auto dtLbl = mDecProdLblPool[idd]; // daughter MC label
        const auto& dtFamily = mSelMCTracks[dtLbl];
        if (dtFamily.mcTrackInfo.pdgParent != dec.pdg) {
          LOGP(error, "{}-th decay (pdg={}): {} in {}:{} range refers to MC track with pdgParent = {}", id, params.decayPDG[id], idd, dec.daughterFirst, dec.daughterLast, dtFamily.mcTrackInfo.pdgParent);
          skip = true;
          break;
        }
        decFam.push_back(dtFamily);
      }
      if (!skip) {
        o2::dataformats::V0 v0;
        if (dec.foundSVID >= 0 && !refitV0(dec.foundSVID, v0, recoData)) {
          v0.invalidate();
        }
        (*mDBGOut) << decTreeName.c_str()
                   << "pdgPar=" << dec.pdg
                   << "trPar=" << dec.parent
                   << "prod=" << decFam
                   << "found=" << dec.foundSVID
                   << "sv=" << v0
                   << "svCheck=" << dec.svCheck
                   << "\n";
      }
    }
  }

  for (auto& mcVtx : mMCVtVec) { // sort rec.vertices in mult. order
    std::sort(mcVtx.recVtx.begin(), mcVtx.recVtx.end(), [](const RecPV& lhs, const RecPV& rhs) {
      return lhs.pv.getNContributors() > rhs.pv.getNContributors();
    });
    (*mDBGOut) << "mcVtxTree" << "mcVtx=" << mcVtx << "\n";
  }

  if (params.storeITSInfo) {
    processITSTracks(recoData);
  }
}

void TrackMCStudy::processTPCTrackRefs()
{
  constexpr float alpsec[18] = {0.174533, 0.523599, 0.872665, 1.221730, 1.570796, 1.919862, 2.268928, 2.617994, 2.967060, 3.316126, 3.665191, 4.014257, 4.363323, 4.712389, 5.061455, 5.410521, 5.759587, 6.108652};
  constexpr float sinAlp[18] = {0.173648, 0.500000, 0.766044, 0.939693, 1.000000, 0.939693, 0.766044, 0.500000, 0.173648, -0.173648, -0.500000, -0.766044, -0.939693, -1.000000, -0.939693, -0.766044, -0.500000, -0.173648};
  constexpr float cosAlp[18] = {0.984808, 0.866025, 0.642788, 0.342020, 0.000000, -0.342020, -0.642788, -0.866025, -0.984808, -0.984808, -0.866025, -0.642788, -0.342020, -0.000000, 0.342020, 0.642788, 0.866025, 0.984808};
  const auto& params = o2::trackstudy::TrackMCStudyConfig::Instance();
  for (auto& entry : mSelMCTracks) {
    auto lb = entry.first;
    auto trspan = mcReader.getTrackRefs(lb.getSourceID(), lb.getEventID(), lb.getTrackID());
    int q = entry.second.mcTrackInfo.track.getCharge();
    if (q * q != 1) {
      continue;
    }
    int ref0entry = mSelTRefs.size(), nrefsSel = 0;
    for (const auto& trf : trspan) {
      if (trf.getDetectorId() != 1) { // process TPC only
        continue;
      }
      float pT = std::hypot(trf.Px(), trf.Py());
      if (pT < 0.05) {
        continue;
      }
      float secX, secY, phi = std::atan2(trf.Y(), trf.X());
      int sector = o2::math_utils::angle2Sector(phi);
      o2::math_utils::rotateZInv(trf.X(), trf.Y(), secX, secY, sinAlp[sector], cosAlp[sector]); // sector coordinates
      float phiPt = std::atan2(trf.Py(), trf.Px());
      o2::math_utils::bringTo02Pi(phiPt);
      auto dphiPt = phiPt - alpsec[sector];
      if (dphiPt > o2::constants::math::PI) { // account for wraps
        dphiPt -= o2::constants::math::TwoPI;
      } else if (dphiPt < -o2::constants::math::PI) {
        dphiPt += o2::constants::math::TwoPI;
      } else if (std::abs(dphiPt) > o2::constants::math::PIHalf * 0.8) {
        continue; // ignore backward going or parallel to padrows tracks
      }
      float tgL = trf.Pz() / pT;
      std::array<float, 5> pars = {secY, trf.Z(), std::sin(dphiPt), tgL, q / pT};
      auto& refTrack = mSelTRefs.emplace_back(secX, alpsec[sector], pars);
      refTrack.setUserField(uint16_t(sector));
      nrefsSel++;
    }
    if (nrefsSel < params.minTPCRefsToExtractClRes) {
      mSelTRefs.resize(ref0entry); // discard unused tracks
      continue;
    } else {
      mSelTRefIdx[lb] = std::make_pair(ref0entry, ref0entry + nrefsSel);
    }
  }
}

void TrackMCStudy::fillMCClusterInfo(const o2::globaltracking::RecoContainer& recoData)
{
  // TPC clusters info
  const auto& TPCClusterIdxStruct = recoData.inputsTPCclusters->clusterIndex;
  const auto* TPCClMClab = recoData.inputsTPCclusters->clusterIndex.clustersMCTruth;
  const auto& params = o2::trackstudy::TrackMCStudyConfig::Instance();

  ClResTPC clRes{};
  for (uint8_t row = 0; row < 152; row++) { // we need to go in increasing row, so this should be the outer loop
    for (uint8_t sector = 0; sector < 36; sector++) {
      unsigned int offs = TPCClusterIdxStruct.clusterOffset[sector][row];
      for (unsigned int icl0 = 0; icl0 < TPCClusterIdxStruct.nClusters[sector][row]; icl0++) {
        const auto labels = TPCClMClab->getLabels(icl0 + offs);
        int ncontLb = 0; // number of real contrubutors to this label (w/o noise)
        for (const auto& lbl : labels) {
          if (!lbl.isValid()) {
            continue;
          }
          ncontLb++;
        }
        const auto& clus = TPCClusterIdxStruct.clusters[sector][row][icl0];
        int tbinH = int(clus.getTime() * mNTPCOccBinLengthInv); // time bin converted to slot of the occ. histo
        clRes.contTracks.clear();
        bool doClusRes = (params.minTPCRefsToExtractClRes > 0) && (params.rejectClustersResStat <= 0. || gRandom->Rndm() < params.rejectClustersResStat);
        for (auto lbl : labels) {
          bool corrAttach = lbl.isFake(); // was this flagged in the flagTPCClusters called from process ?
          lbl.setFakeFlag(false);
          auto entry = mSelMCTracks.find(lbl);
          if (entry == mSelMCTracks.end()) { // not selected
            continue;
          }
          auto& mctr = entry->second.mcTrackInfo;
          mctr.nTPCCl++;
          if (row > mctr.maxTPCRow) {
            mctr.maxTPCRow = row;
            mctr.maxTPCRowSect = sector;
            mctr.nUsedPadRows++;
          } else if (row == 0 && mctr.nUsedPadRows == 0) {
            mctr.nUsedPadRows++;
          }
          if (row < mctr.minTPCRow) {
            mctr.minTPCRow = row;
            mctr.minTPCRowSect = sector;
          }
          if (mctr.minTPCRowSect == sector && row > mctr.maxTPCRowInner) {
            mctr.maxTPCRowInner = row;
          }
          if (ncontLb > 1) {
            mctr.nTPCClShared++;
          }
          // try to extract ideal track position
          if (doClusRes) {
            auto entTRefIDsIt = mSelTRefIdx.find(lbl);
            if (entTRefIDsIt == mSelTRefIdx.end()) {
              continue;
            }
            float xc, yc, zc;
            mTPCCorrMaps->Transform(sector, row, clus.getPad(), clus.getTime(), xc, yc, zc, mctr.bcInTF / 8.); // nominal time of the track

            const auto& entTRefIDs = entTRefIDsIt->second;
            // find bracketing TRef params
            int entIDBelow = -1, entIDAbove = -1;
            float xBelow = -1e6, xAbove = 1e6;

            for (int entID = entTRefIDs.first; entID < entTRefIDs.second; entID++) {
              const auto& refTr = mSelTRefs[entID];
              if (refTr.getUserField() != sector % 18) {
                continue;
              }
              if ((refTr.getX() < xc) && (refTr.getX() > xBelow) && (refTr.getX() > xc - params.maxTPCRefExtrap)) {
                xBelow = refTr.getX();
                entIDBelow = entID;
              }
              if ((refTr.getX() > xc) && (refTr.getX() < xAbove) && (refTr.getX() < xc + params.maxTPCRefExtrap)) {
                xAbove = refTr.getX();
                entIDAbove = entID;
              }
            }
            if ((entIDBelow < 0 && entIDAbove < 0) || (params.requireTopBottomRefs && (entIDBelow < 0 || entIDAbove < 0))) {
              continue;
            }
            auto prop = o2::base::Propagator::Instance();
            o2::track::TrackPar tparAbove, tparBelow;
            bool okBelow = entIDBelow >= 0 && prop->PropagateToXBxByBz((tparBelow = mSelTRefs[entIDBelow]), xc, 0.99, 2.);
            bool okAbove = entIDAbove >= 0 && prop->PropagateToXBxByBz((tparAbove = mSelTRefs[entIDAbove]), xc, 0.99, 2.);
            if ((!okBelow && !okAbove) || (params.requireTopBottomRefs && (!okBelow || !okAbove))) {
              continue;
            }

            int nmeas = 0;
            auto& clCont = clRes.contTracks.emplace_back();
            clCont.corrAttach = corrAttach;
            if (okBelow) {
              clCont.below = {mSelTRefs[entIDBelow].getX(), tparBelow.getY(), tparBelow.getZ()};
              clCont.snp += tparBelow.getSnp();
              clCont.tgl += tparBelow.getTgl();
              clCont.q2pt += tparBelow.getQ2Pt();
              nmeas++;
            }
            if (okAbove) {
              clCont.above = {mSelTRefs[entIDAbove].getX(), tparAbove.getY(), tparAbove.getZ()};
              clCont.snp += tparAbove.getSnp();
              clCont.tgl += tparAbove.getTgl();
              clCont.q2pt += tparAbove.getQ2Pt();
              nmeas++;
            }
            if (nmeas) {
              if (clRes.contTracks.size() == 1) {
                int occBin = mctr.bcInTF / 8 * mNTPCOccBinLengthInv;
                clRes.occ = occBin < 0 ? mTBinClOcc[0] : (occBin >= mTBinClOcc.size() ? mTBinClOcc.back() : mTBinClOcc[occBin]);
              }
              clCont.xyz = {xc, yc, zc};
              if (nmeas > 1) {
                clCont.snp *= 0.5;
                clCont.tgl *= 0.5;
                clCont.q2pt *= 0.5;
              }
            } else {
              clRes.contTracks.pop_back();
            }
          }
        }
        if (clRes.getNCont()) {
          clRes.sect = sector;
          clRes.row = row;
          clRes.qtot = clus.getQtot();
          clRes.qmax = clus.getQmax();
          clRes.flags = clus.getFlags();
          clRes.sigmaTimePacked = clus.sigmaTimePacked;
          clRes.sigmaPadPacked = clus.sigmaPadPacked;
          clRes.ncont = ncontLb;
          clRes.sortCont();

          if (tbinH < 0) {
            tbinH = 0;
          } else if (tbinH >= int(mTBinClOccHist.size())) {
            tbinH = (int)mTBinClOccHist.size() - 1;
          }
          clRes.occBin = mTBinClOccHist[tbinH];

          (*mDBGOut) << "clres" << "clr=" << clRes << "\n";
        }
      }
    }
  }
  // fill ITS cluster info
  const auto* mcITSClusters = recoData.getITSClustersMCLabels();
  const auto& ITSClusters = recoData.getITSClusters();
  for (unsigned int icl = 0; icl < ITSClusters.size(); icl++) {
    const auto labels = mcITSClusters->getLabels(icl);
    for (const auto& lbl : labels) {
      auto entry = mSelMCTracks.find(lbl);
      if (entry == mSelMCTracks.end()) { // not selected
        continue;
      }
      auto& mctr = entry->second.mcTrackInfo;
      mctr.nITSCl++;
      mctr.pattITSCl |= 0x1 << o2::itsmft::ChipMappingITS::getLayer(ITSClusters[icl].getChipID());
    }
  }
}

bool TrackMCStudy::propagateToRefX(o2::track::TrackParCov& trcTPC, o2::track::TrackParCov& trcITS)
{
  bool refReached = false;
  constexpr float TgHalfSector = 0.17632698f;
  const auto& par = o2::globaltracking::MatchTPCITSParams::Instance();
  int trialsLeft = 2;
  while (o2::base::Propagator::Instance()->PropagateToXBxByBz(trcTPC, par.XMatchingRef, MaxSnp, 2., par.matCorr)) {
    if (refReached) {
      break;
    }
    // make sure the track is indeed within the sector defined by alpha
    if (fabs(trcTPC.getY()) < par.XMatchingRef * TgHalfSector) {
      refReached = true;
      break; // ok, within
    }
    if (!trialsLeft--) {
      break;
    }
    auto alphaNew = o2::math_utils::angle2Alpha(trcTPC.getPhiPos());
    if (!trcTPC.rotate(alphaNew) != 0) {
      break; // failed (RS: check effect on matching tracks to neighbouring sector)
    }
  }
  if (!refReached) {
    return false;
  }
  refReached = false;
  float alp = trcTPC.getAlpha();
  return !(trcITS.rotate(alp) == 0) && o2::base::Propagator::Instance()->PropagateToXBxByBz(trcITS, par.XMatchingRef, MaxSnp, 2., par.matCorr);
}

void TrackMCStudy::endOfStream(EndOfStreamContext& ec)
{
  mDBGOut.reset();
}

void TrackMCStudy::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  if (mTPCVDriftHelper.accountCCDBInputs(matcher, obj)) {
    return;
  }
  if (matcher == ConcreteDataMatcher("ITS", "ALPIDEPARAM", 0)) {
    LOG(info) << "ITS Alpide param updated";
    const auto& par = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>::Instance();
    par.printKeyValues(true, true);
    mITSTimeBiasMUS = par.roFrameBiasInBC * o2::constants::lhc::LHCBunchSpacingNS * 1e-3;
    mITSROFrameLengthMUS = par.roFrameLengthInBC * o2::constants::lhc::LHCBunchSpacingNS * 1e-3;
    return;
  }
  if (matcher == ConcreteDataMatcher("ITS", "CLUSDICT", 0)) {
    LOG(info) << "cluster dictionary updated";
    mITSDict = (const o2::itsmft::TopologyDictionary*)obj;
    return;
  }
  if (matcher == ConcreteDataMatcher("GLO", "MEANVERTEX", 0)) {
    LOG(info) << "Imposing new MeanVertex: " << ((const o2::dataformats::MeanVertexObject*)obj)->asString();
    mSVertexer.setMeanVertex((const o2::dataformats::MeanVertexObject*)obj);
    return;
  }
  if (matcher == ConcreteDataMatcher("GLO", "SVPARAM", 0)) {
    LOG(info) << "SVertexer Params updated from ccdb";
    const auto& par = o2::vertexing::SVertexerParams::Instance();
    par.printKeyValues(true, true);
    return;
  }
}

//_____________________________________________________
void TrackMCStudy::prepareITSData(const o2::globaltracking::RecoContainer& recoData)
{
  const auto ITSTracksArray = recoData.getITSTracks();
  const auto ITSTrackROFRec = recoData.getITSTracksROFRecords();
  int nROFs = ITSTrackROFRec.size();
  mITSROF.clear();
  mITSROFBracket.clear();
  mITSROF.reserve(ITSTracksArray.size());
  mITSROFBracket.reserve(ITSTracksArray.size());
  for (int irof = 0; irof < nROFs; irof++) {
    const auto& rofRec = ITSTrackROFRec[irof];
    long nBC = rofRec.getBCData().differenceInBC(recoData.startIR);
    float tMin = nBC * o2::constants::lhc::LHCBunchSpacingMUS + mITSTimeBiasMUS;
    float tMax = tMin + mITSROFrameLengthMUS;
    mITSROFBracket.emplace_back(tMin, tMax);
    for (int it = 0; it < rofRec.getNEntries(); it++) {
      mITSROF.push_back(irof);
    }
  }
}
/*
float TrackMCStudy::getDCAYCut(float pt) const
{
  static TF1 fun("dcayvspt", mDCAYFormula.c_str(), 0, 20);
  return fun.Eval(pt);
}
*/

void TrackMCStudy::collectGammaConversions(int src, int ev)
{
  const auto& params = o2::trackstudy::TrackMCStudyConfig::Instance();
  const auto& tracks = *mCurrMCTracks;
  for (int photonID = 0; photonID < (int)tracks.size(); photonID++) {
    const auto& photon = tracks[photonID];
    if (photon.GetPdgCode() != kGamma) {
      continue;
    }
    const bool physicalPrimary = o2::mcutils::MCTrackNavigator::isPhysicalPrimary(photon, tracks);
    const bool producedByGenerator = photon.getProcess() == TMCProcess::kPPrimary;
    if (params.gammaRequirePrimaryOrGenerator && !physicalPrimary && !producedByGenerator) {
      continue;
    }
    if (photon.GetPt() < params.gammaMinPt || photon.GetPt() > params.gammaMaxPt ||
        std::abs(photon.GetEta()) > params.gammaMaxEta) {
      continue;
    }

    const int firstDaughter = photon.getFirstDaughterTrackId();
    const int lastDaughter = photon.getLastDaughterTrackId();
    if (firstDaughter < 0 || lastDaughter < firstDaughter || lastDaughter >= (int)tracks.size()) {
      continue;
    }
    int positronID = -1, electronID = -1;
    std::array<float, 3> conversionXYZ{};
    for (int daughterID = firstDaughter; daughterID <= lastDaughter; daughterID++) {
      const auto& daughter = tracks[daughterID];
      if (daughter.getMotherTrackId() != photonID || daughter.getProcess() != TMCProcess::kPPair) {
        continue;
      }
      if (daughter.GetPdgCode() == kPositron && positronID < 0) {
        positronID = daughterID;
        conversionXYZ = {(float)daughter.Vx(), (float)daughter.Vy(), (float)daughter.Vz()};
      } else if (daughter.GetPdgCode() == kElectron && electronID < 0) {
        electronID = daughterID;
        conversionXYZ = {(float)daughter.Vx(), (float)daughter.Vy(), (float)daughter.Vz()};
      }
    }
    if (positronID < 0 || electronID < 0) {
      continue;
    }

    const float conversionRadius = std::hypot(conversionXYZ[0], conversionXYZ[1]);
    if (conversionRadius < params.gammaMinR || conversionRadius > params.gammaMaxR) {
      continue;
    }
    const float rzSlope = std::tan(2.f * std::atan(std::exp(-params.gammaMaxEta)));
    if (params.gammaApplyRZLineCut &&
        conversionRadius < std::abs(conversionXYZ[2]) * rzSlope - params.gammaRZMargin) {
      continue;
    }

    const std::array<o2::MCCompLabel, GammaConvInfo::NLegs> daughterLabels{
      o2::MCCompLabel(positronID, ev, src), o2::MCCompLabel(electronID, ev, src)};
    bool daughtersRegistered = true;
    for (int leg = 0; leg < GammaConvInfo::NLegs; leg++) {
      const auto label = daughterLabels[leg];
      const auto& daughter = tracks[label.getTrackID()];
      auto* pdg = O2DatabasePDG::Instance()->GetParticle(daughter.GetPdgCode());
      if (!pdg || !addMCParticle(daughter, label, pdg)) {
        LOGP(error, "Failed to register selected conversion daughter {}", label.asString());
        daughtersRegistered = false;
        break;
      }
    }
    if (!daughtersRegistered) {
      continue;
    }

    auto& conversion = mGammaConversions.emplace_back();
    conversion.photonLabel = o2::MCCompLabel(photonID, ev, src);
    conversion.daughterLabels = daughterLabels;
    conversion.trueConversionXYZ = conversionXYZ;
    conversion.photonPt = photon.GetPt();
    conversion.photonEta = photon.GetEta();
    conversion.photonPhi = photon.GetPhi();
    for (const auto& label : daughterLabels) {
      mGammaDaughterLabels.insert(label);
    }
  }
}

bool TrackMCStudy::processMCParticle(int src, int ev, int trid)
{
  const auto& mcPart = (*mCurrMCTracks)[trid];
  int pdg = mcPart.GetPdgCode();
  bool res = false;
  while (true) {
    auto lbl = o2::MCCompLabel(trid, ev, src);
    int decay = -1; // is this decay to watch?
    const auto& params = o2::trackstudy::TrackMCStudyConfig::Instance();
    if (mcPart.T() < params.decayMotherMaxT) {
      for (int id = 0; id < mNCheckDecays; id++) {
        if (params.decayPDG[id] == std::abs(pdg)) {
          decay = id;
          break;
        }
      }
      if (decay >= 0) { // check if decay and kinematics is acceptable
        auto& decayPool = mDecaysMaps[decay];
        int idd0 = mcPart.getFirstDaughterTrackId(), idd1 = mcPart.getLastDaughterTrackId(); // we want only charged and trackable daughters
        int dtStart = mDecProdLblPool.size(), dtEnd = -1;
        if (idd0 < 0) {
          break;
        }
        for (int idd = idd0; idd <= idd1; idd++) {
          const auto& product = (*mCurrMCTracks)[idd];
          auto lbld = o2::MCCompLabel(idd, ev, src);
          if (!acceptMCCharged(product, lbld, decay)) {
            decay = -1; // discard decay
            mDecProdLblPool.resize(dtStart);
            break;
          }
          mDecProdLblPool.push_back(lbld); // register prong entry and label
        }
        if (decay >= 0) {
          // account decay
          dtEnd = mDecProdLblPool.size();
          for (int dtid = dtStart; dtid < dtEnd; dtid++) { // flag selected decay parent entry in the prongs MCs
            mSelMCTracks[mDecProdLblPool[dtid]].mcTrackInfo.parentEntry = decayPool.size();
            mSelMCTracks[mDecProdLblPool[dtid]].mcTrackInfo.parentDecID = int8_t(decay);
          }
          dtEnd--;
          std::array<float, 3> xyz{(float)mcPart.GetStartVertexCoordinatesX(), (float)mcPart.GetStartVertexCoordinatesY(), (float)mcPart.GetStartVertexCoordinatesZ()};
          std::array<float, 3> pxyz{(float)mcPart.GetStartVertexMomentumX(), (float)mcPart.GetStartVertexMomentumY(), (float)mcPart.GetStartVertexMomentumZ()};
          decayPool.emplace_back(DecayRef{.mother = lbl,
                                          .parent = o2::track::TrackPar(xyz, pxyz, TMath::Nint(O2DatabasePDG::Instance()->GetParticle(mcPart.GetPdgCode())->Charge() / 3), false),
                                          .pdg = mcPart.GetPdgCode(),
                                          .daughterFirst = dtStart,
                                          .daughterLast = dtEnd});
          if (mVerbose > 1) {
            LOGP(info, "Adding MC parent pdg={} {}, with prongs in {}:{} range", pdg, lbl.asString(), dtStart, dtEnd);
          }
          res = true; // Accept!
        }
        break;
      }
    }
    // check if this is a charged which should be processed but was not accounted as a decay product
    if (mSelMCTracks.find(lbl) == mSelMCTracks.end()) {
      res = acceptMCCharged(mcPart, lbl);
    }
    break;
  }
  return res;
}

bool TrackMCStudy::acceptMCCharged(const MCTrack& tr, const o2::MCCompLabel& lb, int followDecay)
{
  const auto& params = o2::trackstudy::TrackMCStudyConfig::Instance();
  // Species filter, applied before anything expensive. Daughters of a watched decay bypass it, so
  // restricting the species can never make a watched decay disappear from the denominator.
  if (followDecay < 0 && !params.selectMCPDG.empty() &&
      std::find(params.selectMCPDG.begin(), params.selectMCPDG.end(), tr.GetPdgCode()) == params.selectMCPDG.end()) {
    return false;
  }
  if (tr.GetPt() < params.minPtMC ||
      std::abs(tr.GetTgl()) > params.maxTglMC ||
      tr.R2() > params.maxRMC * params.maxRMC) {
    if (mVerbose > 1 && followDecay > -1) {
      LOGP(info, "rejecting decay {} prong : pdg={}, pT={}, tgL={}, r={}", followDecay, tr.GetPdgCode(), tr.GetPt(), tr.GetTgl(), std::sqrt(tr.R2()));
    }
    return false;
  }
  float dx = tr.GetStartVertexCoordinatesX() - mCurrMCVertex.X(), dy = tr.GetStartVertexCoordinatesY() - mCurrMCVertex.Y(), dz = tr.GetStartVertexCoordinatesZ() - mCurrMCVertex.Z();
  float r2 = (dx * dx) + (dy * dy);
  float posTgl2 = r2 > 1 && std::abs(dz) < 20 ? dz * dz / r2 : 0;
  if (posTgl2 > params.maxPosTglMC * params.maxPosTglMC) {
    if (mVerbose > 1 && followDecay > -1) {
      LOGP(info, "rejecting decay {} prong : pdg={}, pT={}, tgL={}, dr={}, dz={} r={}, z={}, posTgl={}", followDecay, tr.GetPdgCode(), tr.GetPt(), tr.GetTgl(), std::sqrt(r2), dz, std::sqrt(tr.R2()), tr.GetStartVertexCoordinatesZ(), std::sqrt(posTgl2));
    }
    return false;
  }
  if (params.requireITSorTPCTrackRefs) {
    auto trspan = mcReader.getTrackRefs(lb.getSourceID(), lb.getEventID(), lb.getTrackID());
    bool ok = false;
    for (const auto& trf : trspan) {
      if (trf.getDetectorId() == DetID::ITS || trf.getDetectorId() == DetID::TPC) {
        ok = true;
        break;
      }
    }
    if (!ok) {
      return false;
    }
  }
  TParticlePDG* pPDG = O2DatabasePDG::Instance()->GetParticle(tr.GetPdgCode());
  if (!pPDG) {
    LOGP(debug, "Unknown particle {}", tr.GetPdgCode());
    return false;
  }
  if (pPDG->Charge() == 0.) {
    return false;
  }
  return addMCParticle(tr, lb, pPDG);
}

bool TrackMCStudy::addMCParticle(const MCTrack& mcPart, const o2::MCCompLabel& lb, TParticlePDG* pPDG)
{
  std::array<float, 3> xyz{(float)mcPart.GetStartVertexCoordinatesX(), (float)mcPart.GetStartVertexCoordinatesY(), (float)mcPart.GetStartVertexCoordinatesZ()};
  std::array<float, 3> pxyz{(float)mcPart.GetStartVertexMomentumX(), (float)mcPart.GetStartVertexMomentumY(), (float)mcPart.GetStartVertexMomentumZ()};
  if (!pPDG && !(pPDG = O2DatabasePDG::Instance()->GetParticle(mcPart.GetPdgCode()))) {
    LOGP(debug, "Unknown particle {}", mcPart.GetPdgCode());
    return false;
  }
  auto& mcEntry = mSelMCTracks[lb];
  mcEntry.mcTrackInfo.pdg = mcPart.GetPdgCode();
  mcEntry.mcTrackInfo.track = o2::track::TrackPar(xyz, pxyz, TMath::Nint(pPDG->Charge() / 3), true);
  mcEntry.mcTrackInfo.label = lb;
  mcEntry.mcTrackInfo.bcInTF = mIntBC[lb.getEventID()];
  mcEntry.mcTrackInfo.occTPC = mTPCOcc[lb.getEventID()];
  mcEntry.mcTrackInfo.occITS = mITSOcc[lb.getEventID()];
  mcEntry.mcTrackInfo.occTPCV = mMCVtVec[lb.getEventID()].occTPCV;
  if (mRecProcStage) {
    mcEntry.mcTrackInfo.setAddedAtRecStage();
  }
  if (o2::mcutils::MCTrackNavigator::isPhysicalPrimary(mcPart, *mCurrMCTracks)) {
    mcEntry.mcTrackInfo.setPrimary();
  }
  int moth = -1;
  o2::MCCompLabel mclbPar;
  if ((moth = mcPart.getMotherTrackId()) >= 0) {
    const auto& mcPartPar = (*mCurrMCTracks)[moth];
    mcEntry.mcTrackInfo.pdgParent = mcPartPar.GetPdgCode();
  }
  if (mcPart.isPrimary() && mcReader.getNEvents(lb.getSourceID()) == mMCVtVec.size()) {
    mMCVtVec[lb.getEventID()].nTrackSel++;
  }
  if (mVerbose > 1) {
    LOGP(info, "Adding charged MC pdg={} {} ", mcPart.GetPdgCode(), lb.asString());
  }
  return true;
}

bool TrackMCStudy::refitV0(int i, o2::dataformats::V0& v0, const o2::globaltracking::RecoContainer& recoData)
{
  const auto& id = recoData.getV0sIdx()[i];
  auto seedP = recoData.getTrackParam(id.getProngID(0));
  auto seedN = recoData.getTrackParam(id.getProngID(1));
  bool isTPConly = (id.getProngID(0).getSource() == GTrackID::TPC) || (id.getProngID(1).getSource() == GTrackID::TPC);
  const auto& svparam = o2::vertexing::SVertexerParams::Instance();
  if (svparam.mTPCTrackPhotonTune && isTPConly) {
    mFitterV0.setMaxDZIni(svparam.mTPCTrackMaxDZIni);
    mFitterV0.setMaxDXYIni(svparam.mTPCTrackMaxDXYIni);
    mFitterV0.setMaxChi2(svparam.mTPCTrackMaxChi2);
    mFitterV0.setCollinear(true);
  }
  int nCand = mFitterV0.process(seedP, seedN);
  if (svparam.mTPCTrackPhotonTune && isTPConly) { // restore
    // Reset immediately to the defaults
    mFitterV0.setMaxDZIni(svparam.maxDZIni);
    mFitterV0.setMaxDXYIni(svparam.maxDXYIni);
    mFitterV0.setMaxChi2(svparam.maxChi2);
    mFitterV0.setCollinear(false);
  }
  if (nCand == 0) { // discard this pair
    return false;
  }
  const int cand = 0;
  if (!mFitterV0.isPropagateTracksToVertexDone(cand) && !mFitterV0.propagateTracksToVertex(cand)) {
    return false;
  }
  const auto& trPProp = mFitterV0.getTrack(0, cand);
  const auto& trNProp = mFitterV0.getTrack(1, cand);
  std::array<float, 3> pP{}, pN{};
  trPProp.getPxPyPzGlo(pP);
  trNProp.getPxPyPzGlo(pN);
  std::array<float, 3> pV0 = {pP[0] + pN[0], pP[1] + pN[1], pP[2] + pN[2]};
  auto p2V0 = (pV0[0] * pV0[0]) + (pV0[1] * pV0[1]) + (pV0[2] * pV0[2]);
  const auto& pv = recoData.getPrimaryVertex(id.getVertexID());
  const auto v0XYZ = mFitterV0.getPCACandidatePos(cand);
  float dx = v0XYZ[0] - pv.getX(), dy = v0XYZ[1] - pv.getY(), dz = v0XYZ[2] - pv.getZ(), prodXYZv0 = (dx * pV0[0]) + (dy * pV0[1]) + (dz * pV0[2]);
  float cosPA = prodXYZv0 / std::sqrt((dx * dx + dy * dy + dz * dz) * p2V0);
  new (&v0) o2::dataformats::V0(v0XYZ, pV0, mFitterV0.calcPCACovMatrixFlat(cand), trPProp, trNProp);
  v0.setDCA(mFitterV0.getChi2AtPCACandidate(cand));
  v0.setCosPA(cosPA);
  return true;
}

//_____________________________________________________
// Replay the SVertexer V0 selection on the prongs of every watched MC decay, to learn at which
// stage a decay which was in principle reconstructable was actually lost.
void TrackMCStudy::checkSVertexer(const o2::globaltracking::RecoContainer& recoData)
{
  using SV = o2::vertexing::SVertexer;
  const auto& pools = mSVertexer.getTracksPool(); // prepareSeeds() is called once by process()
  const auto& seedRejMap = mSVertexer.getSeedRejMap();

  // prong pairs of the V0s the reconstruction actually produced. This, and not the label-based
  // foundSVID, is the reference for the replay: it answers directly whether the production
  // SVertexer built a V0 out of these two tracks, without depending on any MC matching.
  auto pairKey = [](VTIndex a, VTIndex b) { return (uint64_t(a) << 32) | uint64_t(b); };
  std::unordered_set<uint64_t> recV0Pairs;
  for (const auto& v0i : recoData.getV0sIdx()) {
    recV0Pairs.insert(pairKey(v0i.getProngID(0), v0i.getProngID(1)));
  }

  // reco track -> position in the SVertexer seeds pool
  std::unordered_map<VTIndex, std::pair<int8_t, int32_t>> poolMap;
  for (int side = 0; side < 2; side++) {
    for (size_t i = 0; i < pools[side].size(); i++) {
      poolMap[pools[side][i].gid] = {int8_t(side), int32_t(i)};
    }
  }

  // recTracks are sorted best-first, so the first one which is in the pool is the one the
  // SVertexer would have used for this MC particle
  auto fillProng = [&pools, &poolMap, &seedRejMap](const TrackFamily& fam, SVProngInfo& pr) {
    for (const auto& rt : fam.recTracks) {
      auto ent = poolMap.find(rt.gid);
      if (ent == poolMap.end()) {
        continue;
      }
      pr.gid = rt.gid;
      pr.poolSide = ent->second.first;
      pr.poolEntry = ent->second.second;
      const auto& seed = pools[pr.poolSide][pr.poolEntry];
      pr.vBrMin = seed.vBracket.getMin();
      pr.vBrMax = seed.vBracket.getMax();
      pr.minR = seed.minR;
      pr.hasTPC = seed.hasTPC;
      pr.nITSclu = seed.nITSclu;
      return;
    }
    // reconstructed, but never made it into the pool: report the best track the SVertexer
    // explicitly dropped, falling back to the best one overall if none was recorded
    if (!fam.recTracks.empty()) {
      pr.gid = fam.recTracks.front().gid;
      for (const auto& rt : fam.recTracks) {
        auto rj = seedRejMap.find(rt.gid);
        if (rj != seedRejMap.end()) {
          pr.gid = rt.gid;
          pr.seedRej = rj->second;
          break;
        }
      }
    }
  };

  std::array<int, o2::vertexing::SVertexer::NRejV0> disagrByRej{};
  int nChecked = 0, nInconsistent = 0;
  for (int id = 0; id < mNCheckDecays; id++) {
    for (auto& dec : mDecaysMaps[id]) {
      auto& chk = dec.svCheck;
      chk.foundSVID = dec.foundSVID;
      if (dec.daughterLast - dec.daughterFirst != 1) {
        continue; // not a 2-prong decay, leave as NotChecked
      }
      for (int ip = 0; ip < 2; ip++) {
        auto ent = mSelMCTracks.find(mDecProdLblPool[dec.daughterFirst + ip]);
        if (ent != mSelMCTracks.end()) {
          fillProng(ent->second, chk.prongs[ip]);
        }
      }
      if (!chk.bothProngsReconstructed()) {
        chk.stage = SVCheck::NoProngs;
        continue;
      }
      if (!chk.bothProngsSeeded()) {
        chk.stage = SVCheck::NotSeeded;
        continue;
      }
      int ipP = -1, ipN = -1;
      for (int ip = 0; ip < 2; ip++) {
        if (chk.prongs[ip].poolSide == SV::POS) {
          ipP = ip;
        } else {
          ipN = ip;
        }
      }
      if (ipP < 0 || ipN < 0) {
        chk.stage = SVCheck::SameCharge;
        continue;
      }
      const auto& seedP = pools[SV::POS][chk.prongs[ipP].poolEntry];
      const auto& seedN = pools[SV::NEG][chk.prongs[ipN].poolEntry];
      if (seedP.vBracket.getOverlap(seedN.vBracket).isInvalid()) {
        chk.stage = SVCheck::NoBracketOverlap; // the pair is never even tried
        continue;
      }
      auto rej = mSVertexer.checkV0(seedP, seedN, chk.prongs[ipP].poolEntry, chk.prongs[ipN].poolEntry, 0);
      chk.rejV0 = uint8_t(rej);
      chk.stage = rej == SV::RejNone ? SVCheck::Found : SVCheck::CutRejected;
      // does the replay agree with the V0 list the reconstruction actually produced for this pair?
      bool inReco = recV0Pairs.contains(pairKey(chk.prongs[ipP].gid, chk.prongs[ipN].gid));
      chk.pairInReco = inReco;
      chk.replayConsistent = (rej == SV::RejNone) == inReco;
      nChecked++;
      if (!chk.replayConsistent) {
        nInconsistent++;
        disagrByRej[rej]++;
        LOGP(debug, "SVertexer replay disagrees for decay pdg={}: reconstruction {} this pair, replay returned {}",
             dec.pdg, inReco ? "kept" : "dropped", o2::vertexing::SVertexer::V0RejNames[rej]);
      }
    }
  }
  // A disagreement does not invalidate one decay, it invalidates the comparison: it means the replay
  // is not seeing the conditions the reconstruction saw. Report the rate rather than aborting, the
  // per-decay flag is in the tree so the affected entries can be cut away offline.
  if (nInconsistent) {
    std::string brk;
    for (int i = 0; i < (int)disagrByRej.size(); i++) {
      if (disagrByRej[i]) {
        // RejNone means the reconstruction dropped a pair the replay kept, any other code means the
        // reconstruction kept a pair the replay rejected on that particular cut
        bool replayKept = (i == o2::vertexing::SVertexer::RejNone);
        brk += fmt::format(" [reco {} / replay {}: {}]", replayKept ? "dropped" : "kept",
                           replayKept ? "kept" : o2::vertexing::SVertexer::V0RejNames[i], disagrByRej[i]);
      }
    }
    LOGP(warn, "SVertexer replay disagrees with the reconstruction for {}/{} pairs ({:.2f}%):{}",
         nInconsistent, nChecked, nChecked ? 100. * nInconsistent / nChecked : 0., brk);
  }
}

//_____________________________________________________
// Build the conversion-specific, clone-aware efficiency record. Unlike the generic decay output,
// this path deliberately applies no daughter MC or reconstructed-track kinematic quality cuts.
void TrackMCStudy::processGammaConversions(const o2::globaltracking::RecoContainer& recoData)
{
  if (mGammaConversions.empty()) {
    return;
  }
  const auto& params = o2::trackstudy::TrackMCStudyConfig::Instance();
  const auto pvvec = recoData.getPrimaryVertices();
  const auto pvvecLbl = recoData.getPrimaryVertexMCLabels();
  const auto v0IDs = recoData.getV0sIdx();

  std::unordered_map<o2::MCCompLabel, std::vector<std::array<VTIndex, 2>>> gammaV0Pairs;
  gammaV0Pairs.reserve(v0IDs.size());
  for (const auto& v0ID : v0IDs) {
    const auto positronTrackID = v0ID.getProngID(0);
    const auto electronTrackID = v0ID.getProngID(1);
    auto positronLabel = recoData.getTrackMCLabel(positronTrackID);
    auto electronLabel = recoData.getTrackMCLabel(electronTrackID);
    if (!positronLabel.isValid() || !electronLabel.isValid()) {
      continue;
    }
    positronLabel.setFakeFlag(false);
    electronLabel.setFakeFlag(false);
    if (positronLabel == electronLabel ||
        positronLabel.getSourceID() != electronLabel.getSourceID() ||
        positronLabel.getEventID() != electronLabel.getEventID()) {
      continue;
    }
    const auto* positronMC = mcReader.getTrack(positronLabel);
    const auto* electronMC = mcReader.getTrack(electronLabel);
    if (!positronMC || !electronMC ||
        positronMC->GetPdgCode() != kPositron || electronMC->GetPdgCode() != kElectron ||
        positronMC->getProcess() != TMCProcess::kPPair || electronMC->getProcess() != TMCProcess::kPPair ||
        positronMC->getMotherTrackId() < 0 ||
        positronMC->getMotherTrackId() != electronMC->getMotherTrackId()) {
      continue;
    }
    const int motherID = positronMC->getMotherTrackId();
    const auto* motherMC = mcReader.getTrack(positronLabel.getSourceID(), positronLabel.getEventID(), motherID);
    if (!motherMC || motherMC->GetPdgCode() != kGamma) {
      continue;
    }
    const o2::MCCompLabel photonLabel(motherID, positronLabel.getEventID(), positronLabel.getSourceID());
    gammaV0Pairs[photonLabel].push_back({positronTrackID, electronTrackID});
  }

  // Mirror the central-barrel labelMask filled by AODProducerWorkflowDPL::fillMCTrackLabelsTable.
  auto getAODMCTrackMask = [&recoData](GTrackID trackID) {
    uint16_t mask = 0;
    const auto mcTruth = recoData.getTrackMCLabel(trackID);
    if (mcTruth.isValid()) {
      if (mcTruth.isFake()) {
        mask |= uint16_t(1) << 15;
      }
      if (trackID.includesDet(DetID::TPC) && trackID.getSource() != GTrackID::TPC) {
        const auto contributorIDs = recoData.getSingleDetectorRefs(trackID);
        if (contributorIDs[GTrackID::ITSTPC].isIndexSet() &&
            recoData.getTrackMCLabel(contributorIDs[GTrackID::ITSTPC]).isFake()) {
          mask |= uint16_t(1) << 13;
        }
      }
      if (trackID.includesDet(DetID::ITS)) {
        const auto itsID = recoData.getITSContributorGID(trackID);
        if (itsID.getSource() == GTrackID::ITS) {
          const auto& itsTrack = recoData.getITSTrack(itsID);
          for (int layer = 0; layer < 7; layer++) {
            if (itsTrack.isFakeOnLayer(layer)) {
              mask |= uint16_t(1) << layer;
            }
          }
        } else if (itsID.getSource() == GTrackID::ITSAB &&
                   recoData.getTrackMCLabel(itsID).isFake()) {
          mask |= uint16_t(1) << 12;
        }
      }
    } else if (mcTruth.isNoise()) {
      mask |= uint16_t(1) << 14;
    }
    return mask;
  };
  auto getPairType = [](VTIndex positronID, VTIndex electronID) {
    const auto positronMask = positronID.getSourceDetectorsMask();
    const auto electronMask = electronID.getSourceDetectorsMask();
    const bool positronITSTPC = positronMask[DetID::ITS] && positronMask[DetID::TPC];
    const bool electronITSTPC = electronMask[DetID::ITS] && electronMask[DetID::TPC];
    const bool positronTPCOnly = !positronMask[DetID::ITS] && positronMask[DetID::TPC];
    const bool electronTPCOnly = !electronMask[DetID::ITS] && electronMask[DetID::TPC];
    const bool positronITSOnly = positronMask[DetID::ITS] && !positronMask[DetID::TPC];
    const bool electronITSOnly = electronMask[DetID::ITS] && !electronMask[DetID::TPC];

    if (positronITSTPC && electronITSTPC) {
      return GammaConvInfo::PairITSTPCITSTPC;
    }
    if ((positronITSTPC && electronTPCOnly) || (positronTPCOnly && electronITSTPC)) {
      return GammaConvInfo::PairITSTPCTPCOnly;
    }
    if (positronTPCOnly && electronTPCOnly) {
      return GammaConvInfo::PairTPCOnlyTPCOnly;
    }
    if (positronITSOnly || electronITSOnly) {
      return GammaConvInfo::PairContainsITSOnly;
    }
    return GammaConvInfo::PairOther;
  };

  for (auto& conversion : mGammaConversions) {
    conversion.referencePV = -1;
    conversion.referencePVNumContrib = -1;
    for (int pvID = 0; pvID < (int)pvvec.size() && pvID < (int)pvvecLbl.size(); pvID++) {
      const auto& label = pvvecLbl[pvID];
      if (!label.isSet() || label.getSourceID() != conversion.photonLabel.getSourceID() ||
          label.getEventID() != conversion.photonLabel.getEventID()) {
        continue;
      }
      const int nContributors = pvvec[pvID].getNContributors();
      if (nContributors > conversion.referencePVNumContrib) {
        conversion.referencePV = pvID;
        conversion.referencePVNumContrib = nContributors;
      }
    }
    if (conversion.referencePV < 0) {
      continue; // match O2Physics: the efficiency sample is conditional on an accepted reference collision
    }

    std::array<std::unordered_set<VTIndex>, GammaConvInfo::NLegs> eligibleIDs;
    std::array<std::unordered_set<VTIndex>, GammaConvInfo::NLegs> referenceEligibleIDs;
    conversion.referenceEligiblePairTypes.fill(false);
    conversion.rawV0FoundUsingReferenceEligiblePairTypes.fill(false);
    for (int leg = 0; leg < GammaConvInfo::NLegs; leg++) {
      conversion.tracks[leg].clear();
      auto familyIt = mSelMCTracks.find(conversion.daughterLabels[leg]);
      if (familyIt == mSelMCTracks.end()) {
        continue;
      }
      const int expectedSign = leg == GammaConvInfo::Positron ? 1 : -1;
      for (const auto& recTrack : familyIt->second.recTracks) {
        auto& trackInfo = conversion.tracks[leg].emplace_back();
        trackInfo.track = recTrack;
        auto pvIt = mGammaTrackPVRefs.find(recTrack.gid);
        if (pvIt != mGammaTrackPVRefs.end()) {
          trackInfo.pvIDs = pvIt->second;
        }
        // AOD stores an ambiguous track once, with the first collision assignment encountered.
        trackInfo.inReferencePV = recTrack.pvID == conversion.referencePV;
        trackInfo.mcMask = getAODMCTrackMask(recTrack.gid);

        const auto detectorMask = recTrack.gid.getSourceDetectorsMask();
        if (params.gammaRequireCleanTrackLabel && trackInfo.mcMask != 0) {
          trackInfo.rejection = GammaConvTrackInfo::DirtyLabel;
        } else if (recTrack.track.getSign() * expectedSign <= 0) {
          trackInfo.rejection = GammaConvTrackInfo::WrongSign;
        } else if (!detectorMask[DetID::ITS] && !detectorMask[DetID::TPC]) {
          trackInfo.rejection = GammaConvTrackInfo::NoITSorTPC;
        } else if (params.gammaRequireTPC && !detectorMask[DetID::TPC]) {
          trackInfo.rejection = GammaConvTrackInfo::NoTPC;
        } else {
          trackInfo.rejection = GammaConvTrackInfo::Accepted;
          eligibleIDs[leg].insert(recTrack.gid);
          if (trackInfo.inReferencePV) {
            referenceEligibleIDs[leg].insert(recTrack.gid);
          }
        }
      }
    }

    conversion.bothLegsFoundAnywhere = !eligibleIDs[GammaConvInfo::Positron].empty() &&
                                       !eligibleIDs[GammaConvInfo::Electron].empty();
    conversion.bothLegsFoundInReferencePV = !referenceEligibleIDs[GammaConvInfo::Positron].empty() &&
                                            !referenceEligibleIDs[GammaConvInfo::Electron].empty();
    for (const auto positronID : referenceEligibleIDs[GammaConvInfo::Positron]) {
      for (const auto electronID : referenceEligibleIDs[GammaConvInfo::Electron]) {
        conversion.referenceEligiblePairTypes[getPairType(positronID, electronID)] = true;
      }
    }

    conversion.rawV0FoundUsingAnywhereEligiblePair = false;
    conversion.rawV0FoundUsingReferenceEligiblePair = false;
    auto gammaV0It = gammaV0Pairs.find(conversion.photonLabel);
    if (gammaV0It != gammaV0Pairs.end()) {
      for (const auto& pair : gammaV0It->second) {
        if (eligibleIDs[GammaConvInfo::Positron].contains(pair[GammaConvInfo::Positron]) &&
            eligibleIDs[GammaConvInfo::Electron].contains(pair[GammaConvInfo::Electron])) {
          conversion.rawV0FoundUsingAnywhereEligiblePair = true;
        }
        if (referenceEligibleIDs[GammaConvInfo::Positron].contains(pair[GammaConvInfo::Positron]) &&
            referenceEligibleIDs[GammaConvInfo::Electron].contains(pair[GammaConvInfo::Electron])) {
          conversion.rawV0FoundUsingReferenceEligiblePair = true;
          const auto pairType = getPairType(pair[GammaConvInfo::Positron], pair[GammaConvInfo::Electron]);
          conversion.rawV0FoundUsingReferenceEligiblePairTypes[pairType] = true;
        }
      }
    }

    if (conversion.rawV0FoundUsingAnywhereEligiblePair) {
      conversion.terminalReason = GammaConvInfo::V0Stored;
    } else if (!conversion.bothLegsFoundAnywhere) {
      const bool positronFound = !eligibleIDs[GammaConvInfo::Positron].empty();
      const bool electronFound = !eligibleIDs[GammaConvInfo::Electron].empty();
      conversion.terminalReason = !positronFound && !electronFound ? GammaConvInfo::NeitherLegEligible :
                                  !positronFound                   ? GammaConvInfo::NoEligiblePositron :
                                                                    GammaConvInfo::NoEligibleElectron;
    } else {
      conversion.terminalReason = GammaConvInfo::BothLegsFoundNoV0;
    }

    (*mDBGOut) << "gammaConv"
               << "conv=" << conversion
               << "\n";
  }
}

//_____________________________________________________
// Reproduce the TPC-only photon-tune cuts of SVertexer::processTPCTrack and record the continuous
// variables behind them, so that the thresholds can be studied instead of only their outcome.
//
// This duplicates production logic on purpose, to keep the SVertexer untouched. It is affordable
// because correctTPCTrack only shifts Z (X, alpha, snp, tgl and q2pt are left alone), so the helix
// circle is identical on the raw track and only the drift shift has to be redone here. The recomputed
// Z is compared against the SVertexer seed of every accepted track, which makes a future divergence
// of correctTPCTrack show up as a warning rather than as quietly wrong numbers.
void TrackMCStudy::evalTPCPhotonTune(const o2::globaltracking::RecoContainer& recoData)
{
  const auto& params = o2::trackstudy::TrackMCStudyConfig::Instance();
  const auto& svparam = o2::vertexing::SVertexerParams::Instance();
  if (!recoData.isTrackSourceLoaded(GTrackID::TPC) || svparam.mExcludeTPCtracks) {
    return;
  }
  const auto pvvec = recoData.getPrimaryVertices();
  const auto pvvecLbl = recoData.getPrimaryVertexMCLabels();
  const auto trackIndex = recoData.getPrimaryVertexMatchedTracks();
  const auto vtxRefs = recoData.getPrimaryVertexMatchedTrackRefs();
  const float bz = o2::base::Propagator::Instance()->getNominalBz();

  // SVertexer::setTPCTBin is never called anywhere in O2, so mMUS2TPCBin keeps its default of
  // 1/(8 BC). Mirrored here rather than assumed silently.
  const float mus2TPCBin = 1.f / (8 * o2::constants::lhc::LHCBunchSpacingMUS);
  const float tpcBin2Z = mTPCVDriftHelper.getVDriftObject().getVDrift() / mus2TPCBin;

  // Z of every TPC clone the SVertexer actually kept, keyed by track and vertex, for the cross-check
  std::map<std::pair<uint32_t, int>, float> poolZ;
  for (const auto& pool : mSVertexer.getTracksPool()) {
    for (const auto& seed : pool) {
      if (seed.gid.getSource() == GTrackID::TPC) {
        poolZ[{uint32_t(seed.gid), seed.vBracket.getMin()}] = seed.getZ();
      }
    }
  }

  std::vector<TPCTuneInfo> recs;
  int nZMismatch = 0, nZChecked = 0;
  const int nv = (int)vtxRefs.size() - 1; // the last entry holds the unassigned tracks
  for (int iv = 0; iv < nv; iv++) {
    const auto& vtref = vtxRefs[iv];
    const auto& vtx = pvvec[iv];
    int idMin = vtref.getFirstEntryOfSource(GTrackID::TPC);
    int idMax = idMin + vtref.getEntriesOfSource(GTrackID::TPC);
    for (int i = idMin; i < idMax; i++) {
      const auto gid = trackIndex[i];
      if (gid.getSource() != GTrackID::TPC) {
        continue;
      }
      const auto& tTPC = recoData.getTPCTrack(gid);

      TPCTuneInfo rec;
      rec.gid = gid;
      rec.vtxID = iv;
      rec.x = tTPC.getX();
      rec.nClusters = (int16_t)tTPC.getNClusters();

      auto lbl = recoData.getTrackMCLabel(gid);
      rec.labelFake = lbl.isFake();
      if (lbl.isValid()) {
        lbl.setFakeFlag(false);
        rec.mcLabel = lbl;
        rec.isSignal = mGammaDaughterLabels.contains(lbl);
        // is the vertex under test the collision this track actually comes from?
        if (iv < (int)pvvecLbl.size()) {
          const auto& pvLbl = pvvecLbl[iv];
          rec.isCorrectPV = pvLbl.isSet() && pvLbl.getSourceID() == lbl.getSourceID() &&
                            pvLbl.getEventID() == lbl.getEventID();
        }
        if (const auto* mcTr = mcReader.getTrack(lbl)) {
          rec.mcPdg = mcTr->GetPdgCode();
          rec.mcPt = mcTr->GetPt();
          rec.mcR = std::hypot(mcTr->Vx(), mcTr->Vy());
          const int moth = mcTr->getMotherTrackId();
          if (moth >= 0) {
            if (const auto* mcMoth = mcReader.getTrack(lbl.getSourceID(), lbl.getEventID(), moth)) {
              rec.mcMotherPdg = mcMoth->GetPdgCode();
            }
          }
        }
      }
      // signal is always kept, the background is the bulk and gets sampled
      if (!rec.isSignal && params.tpcTuneBkgSamplingFrac < 1.f &&
          gRandom->Rndm() > params.tpcTuneBkgSamplingFrac) {
        continue;
      }

      if (svparam.mTPCTrackMaxX > 0. && tTPC.getX() > svparam.mTPCTrackMaxX) {
        rec.stage = TPCTuneInfo::RejMaxX;
        recs.push_back(rec);
        continue;
      }
      if (tTPC.hasBothSidesClusters()) { // effectively constrained, handled as a normal track
        rec.stage = TPCTuneInfo::BothSides;
        recs.push_back(rec);
        continue;
      }

      const auto twe = vtx.getTimeStamp();
      const float tTB = twe.getTimeStamp() * mus2TPCBin;
      const float driftErr = twe.getTimeStampError() * mus2TPCBin * tpcBin2Z;
      if (driftErr < 0.f) {
        rec.stage = TPCTuneInfo::RejTimeCorr;
        recs.push_back(rec);
        continue;
      }
      const float dDrift = (tTB - tTPC.getTime0()) * tpcBin2Z;
      rec.zCorr = tTPC.getZ() + (tTPC.hasASideClustersOnly() ? dDrift : -dDrift);
      rec.stage = TPCTuneInfo::Evaluated;

      // dDPV: extrapolation of the track back to the beam line in Z
      rec.dz2Beam = std::abs(tTPC.getX() * tTPC.getTgl() - rec.zCorr + vtx.getZ());
      // dRD2: correctTPCTrack does not touch the transverse parameters, so the circle of the raw
      // track is the circle of the clone
      float sna = 0.f, csa = 0.f;
      o2::math_utils::CircleXYf_t circle;
      tTPC.getCircleParams(bz, circle, sna, csa);
      rec.cR = std::hypot(circle.xC, circle.yC);
      rec.rC = circle.rC;
      rec.drd2Sq = (rec.cR * rec.cR) - (rec.rC * rec.rC);

      const bool dCls = rec.nClusters < svparam.mTPCTrackMinNClusters;
      const bool dDPV = rec.dz2Beam > svparam.mTPCTrack2Beam;
      // must track SVertexer::processTPCTrack: a negative argument means the helix encloses the
      // beam line, which is rejected rather than passed through a NaN comparison
      const bool dRD2 = rec.drd2Sq < 0.f || std::sqrt(rec.drd2Sq) > svparam.mTPCTrackXY2Radius;
      rec.accepted = !(dCls || dDPV || dRD2);

      auto poolIt = poolZ.find({uint32_t(gid), iv});
      if (poolIt != poolZ.end()) {
        rec.zPool = poolIt->second;
        nZChecked++;
        if (std::abs(rec.zPool - rec.zCorr) > 1e-3) {
          nZMismatch++;
        }
      }
      recs.push_back(rec);
    }
  }

  if (nZMismatch) {
    LOGP(warn, "TPC photon tune: recomputed Z differs from the SVertexer seed for {}/{} tracks. The local "
               "copy of correctTPCTrack has drifted from SVertexer::correctTPCTrack, the tpcTune output is not usable",
         nZMismatch, nZChecked);
  }

  (*mDBGOut) << "tpcTune" << "orbit=" << recoData.startIR.orbit
             << "minNCl=" << svparam.mTPCTrackMinNClusters
             << "cutZ2Beam=" << svparam.mTPCTrack2Beam
             << "cutXY2Radius=" << svparam.mTPCTrackXY2Radius
             << "cutMaxX=" << svparam.mTPCTrackMaxX
             // needed to undo the sampling: signal is kept in full, background is not, so any
             // absolute rate or signal fraction computed from the stored records is otherwise wrong
             << "bkgSampFrac=" << params.tpcTuneBkgSamplingFrac
             << "trk=" << recs << "\n";
}

//_____________________________________________________
// Store every reconstructed V0 together with the MC origin of its prongs, to study the composition
// of the sample (true decays vs combinatorial).
void TrackMCStudy::processRecSVs(const o2::globaltracking::RecoContainer& recoData)
{
  const auto& params = o2::trackstudy::TrackMCStudyConfig::Instance();
  auto v0IDs = recoData.getV0sIdx();
  std::vector<RecSVInfo> recSVs;
  for (int iv = 0; iv < (int)v0IDs.size(); iv++) {
    RecSVInfo inf;
    inf.v0ID = v0IDs[iv];
    const o2::MCTrack* mcTr[2] = {nullptr, nullptr};
    o2::MCCompLabel lbClean[2]; // labels with the fake flag cleared, as used to key mSelMCTracks
    bool fake = false;
    for (int ip = 0; ip < 2; ip++) {
      auto lbl = recoData.getTrackMCLabel(inf.v0ID.getProngID(ip));
      inf.prongLbl[ip] = lbl; // kept with the fake flag, it is part of the information
      if (!lbl.isValid()) {
        continue;
      }
      fake = fake || lbl.isFake();
      lbl.setFakeFlag(false);
      lbClean[ip] = lbl;
      if ((mcTr[ip] = mcReader.getTrack(lbl))) {
        inf.prongPDG[ip] = mcTr[ip]->GetPdgCode();
      }
    }
    if (!mcTr[0] || !mcTr[1]) {
      inf.kind = RecSVInfo::Unknown;
    } else if (fake) {
      inf.kind = RecSVInfo::FakeProng;
    } else if (lbClean[0].getSourceID() == lbClean[1].getSourceID() &&
               lbClean[0].getEventID() == lbClean[1].getEventID() &&
               mcTr[0]->getMotherTrackId() == mcTr[1]->getMotherTrackId() &&
               mcTr[0]->getMotherTrackId() >= 0) {
      inf.kind = RecSVInfo::TrueDecay;
      const auto* moth = mcReader.getTrack(lbClean[0].getSourceID(), lbClean[0].getEventID(), mcTr[0]->getMotherTrackId());
      if (moth) {
        inf.mcMotherPDG = moth->GetPdgCode();
      }
      // was this the decay of a particle we are watching? if so, record where to find it
      auto ent = mSelMCTracks.find(lbClean[0]);
      if (ent != mSelMCTracks.end()) {
        inf.mcMotherEntry = ent->second.mcTrackInfo.parentEntry;
        inf.mcMotherDecID = ent->second.mcTrackInfo.parentDecID;
      }
    } else {
      inf.kind = RecSVInfo::DifferentMothers;
    }
    // true decays are rare and always kept, the combinatorial ones are the bulk and get sampled
    if (!inf.isTrueDecay() && params.recSVSamplingFrac < 1.f && gRandom->Rndm() > params.recSVSamplingFrac) {
      continue;
    }
    if (!refitV0(iv, inf.v0, recoData)) {
      inf.v0.invalidate();
    }
    recSVs.push_back(inf);
  }
  (*mDBGOut) << "recSV"
             << "orbit=" << recoData.startIR.orbit
             << "sv=" << recSVs
             << "\n";
}

void TrackMCStudy::loadTPCOccMap(const o2::globaltracking::RecoContainer& recoData)
{
  auto NHBPerTF = o2::base::GRPGeomHelper::instance().getGRPECS()->getNHBFPerTF();
  const auto& TPCOccMap = recoData.occupancyMapTPC;
  auto prop = o2::base::Propagator::Instance();
  auto TPCRefitter = std::make_unique<o2::gpu::GPUO2InterfaceRefit>(&recoData.inputsTPCclusters->clusterIndex, mTPCCorrMaps, prop->getNominalBz(),
                                                                    recoData.getTPCTracksClusterRefs().data(), 0, recoData.clusterShMapTPC.data(), TPCOccMap.data(), TPCOccMap.size(), nullptr, prop);
  mNTPCOccBinLength = TPCRefitter->getParam()->rec.tpc.occupancyMapTimeBins;
  mTBinClOcc.clear();
  if (mNTPCOccBinLength > 1 && TPCOccMap.size()) {
    mNTPCOccBinLengthInv = 1. / mNTPCOccBinLength;
    int nTPCBins = NHBPerTF * o2::constants::lhc::LHCMaxBunches / 8, ninteg = 0;
    int nTPCOccBins = nTPCBins * mNTPCOccBinLengthInv, sumBins = std::max(1, int(o2::constants::lhc::LHCMaxBunches / 8 * mNTPCOccBinLengthInv));
    mTBinClOcc.resize(nTPCOccBins);
    mTBinClOccHist.resize(nTPCOccBins);
    float sm = 0., tb = 0.5 * mNTPCOccBinLength;
    for (int i = 0; i < nTPCOccBins; i++) {
      mTBinClOccHist[i] = TPCRefitter->getParam()->GetUnscaledMult(tb);
      tb += mNTPCOccBinLength;
    }
    for (int i = nTPCOccBins; i--;) {
      sm += mTBinClOccHist[i];
      if (i + sumBins < nTPCOccBins) {
        sm -= mTBinClOccHist[i + sumBins];
      }
      mTBinClOcc[i] = sm;
    }
  } else {
    mTBinClOcc.resize(1);
    mTBinClOccHist.resize(1);
  }
}

void TrackMCStudy::processITSTracks(const o2::globaltracking::RecoContainer& recoData)
{
  if (!mITSDict) {
    LOGP(warn, "ITS data is not loaded");
    return;
  }
  const auto itsTracks = recoData.getITSTracks();
  const auto itsLbls = recoData.getITSTracksMCLabels();
  const auto itsClRefs = recoData.getITSTracksClusterRefs();
  const auto clusITS = recoData.getITSClusters();
  const auto patterns = recoData.getITSClustersPatterns();
  const auto& params = o2::trackstudy::TrackMCStudyConfig::Instance();
  auto pattIt = patterns.begin();
  mITSClustersArray.clear();
  mITSClustersArray.reserve(clusITS.size());

  o2::its::ioutils::convertCompactClusters(clusITS, pattIt, mITSClustersArray, mITSDict);
  auto geom = o2::its::GeometryTGeo::Instance();
  int ntr = itsLbls.size();
  LOGP(info, "We have {} ITS clusters and the number of patterns is {}, ITSdict:{} NMCLabels: {}", clusITS.size(), patterns.size(), mITSDict != nullptr, itsLbls.size());

  std::vector<int> evord(ntr);
  std::iota(evord.begin(), evord.end(), 0);
  std::sort(evord.begin(), evord.end(), [&](int i, int j) { return itsLbls[i] < itsLbls[j]; });
  std::vector<ITSHitInfo> outHitInfo;
  std::array<int, 7> cl2arr{};

  for (int itr0 = 0; itr0 < ntr; itr0++) {
    auto itr = evord[itr0];
    const auto& itsTr = itsTracks[itr];
    const auto& itsLb = itsLbls[itr];
    //    LOGP(info,"proc {} {} {}",itr0, itr, itsLb.asString());
    int nCl = itsTr.getNClusters();
    if (itsLb.isFake() || nCl < params.minITSClForITSoutput) {
      continue;
    }
    auto entrySel = mSelMCTracks.find(itsLb);
    if (entrySel == mSelMCTracks.end()) {
      continue;
    }
    outHitInfo.clear();
    cl2arr.fill(-1);
    auto clEntry = itsTr.getFirstClusterEntry();
    for (int iCl = nCl; iCl--;) { // clusters are stored from outer to inner layers
      const auto& cls = mITSClustersArray[itsClRefs[clEntry + iCl]];
      int hpos = outHitInfo.size();
      auto& hinf = outHitInfo.emplace_back();
      hinf.clus = cls;
      hinf.clus.setCount(geom->getLayer(cls.getSensorID()));
      geom->getSensorXAlphaRefPlane(cls.getSensorID(), hinf.chipX, hinf.chipAlpha);
      cl2arr[hinf.clus.getCount()] = hpos; // to facilitate finding the cluster of the layer
    }
    auto trspan = mcReader.getTrackRefs(itsLb.getSourceID(), itsLb.getEventID(), itsLb.getTrackID());
    int ilrc = -1, nrefAcc = 0;
    for (const auto& trf : trspan) {
      if (trf.getDetectorId() != 0) { // process ITS only
        continue;
      }
      int lrt = trf.getUserId(); // layer of the reference, but there might be multiple hits on the same layer
      int clEnt = cl2arr[lrt];
      if (clEnt < 0) {
        continue;
      }
      auto& hinf = outHitInfo[clEnt];
      float traX, traY;
      o2::math_utils::rotateZInv(trf.X(), trf.Y(), traX, traY, std::sin(hinf.chipAlpha), std::cos(hinf.chipAlpha)); // tracking coordinates of the reference
      if (hinf.trefXT < 1 || std::abs(traX - hinf.chipX) < std::abs(hinf.trefXT - hinf.chipX)) {
        if (hinf.trefXT < 1) {
          nrefAcc++;
        }
        hinf.tref = trf;
        hinf.trefXT = traX;
        hinf.trefYT = traY;
      }
    }
    (*mDBGOut) << "itsTree" << "hits=" << outHitInfo << "trIn=" << ((o2::track::TrackParCov&)itsTr) << "trOut=" << itsTr.getParamOut() << "mcTr=" << entrySel->second.mcTrackInfo.track << "mcPDG=" << entrySel->second.mcTrackInfo.pdg << "nTrefs=" << nrefAcc << "\n";
  }
}

DataProcessorSpec getTrackMCStudySpec(GTrackID::mask_t srcTracks, GTrackID::mask_t srcClusters, bool checkSV, bool useCCDBParams, bool enableCasc, bool enable3body)
{
  std::vector<OutputSpec> outputs;
  Options opts{
    {"device-verbosity", VariantType::Int, 0, {"Verbosity level"}},
    {"dcay-vs-pt", VariantType::String, "0.0105 + 0.0350 / pow(x, 1.1)", {"Formula for global tracks DCAy vs pT cut"}},
    {"min-tpc-clusters", VariantType::Int, 60, {"Cut on TPC clusters"}},
    {"max-tpc-dcay", VariantType::Float, 2.f, {"Cut on TPC dcaY"}},
    {"max-tpc-dcaz", VariantType::Float, 2.f, {"Cut on TPC dcaZ"}},
    {"min-x-prop", VariantType::Float, 6.f, {"track should be propagated to this X at least"}}};
  auto dataRequest = std::make_shared<DataRequest>();
  bool useMC = true;
  dataRequest->requestTracks(srcTracks, useMC);
  dataRequest->requestClusters(srcClusters, useMC);
  dataRequest->requestPrimaryVertices(useMC);
  if (checkSV) {
    dataRequest->requestSecondaryVertices(useMC);
    dataRequest->inputs.emplace_back("meanvtx", "GLO", "MEANVERTEX", 0, Lifetime::Condition, ccdbParamSpec("GLO/Calib/MeanVertex", {}, 1));
    if (useCCDBParams) { // the same object the SecondaryVertexingSpec reads, so the cuts are identical
      dataRequest->inputs.emplace_back("SVParam", "GLO", "SVPARAM", 0, Lifetime::Condition, ccdbParamSpec("GLO/Config/SVertexerParam"));
    }
  }
  o2::tpc::VDriftHelper::requestCCDBInputs(dataRequest->inputs);
  dataRequest->inputs.emplace_back("corrMap", o2::header::gDataOriginTPC, "TPCCORRMAP", 0, Lifetime::Timeframe);
  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
                                                              true,                              // GRPECS=true
                                                              true,                              // GRPLHCIF
                                                              true,                              // GRPMagField
                                                              true,                              // askMatLUT
                                                              o2::base::GRPGeomRequest::Aligned, // geometry
                                                              dataRequest->inputs,
                                                              true);

  return DataProcessorSpec{
    .name = "track-mc-study",
    .inputs = dataRequest->inputs,
    .outputs = outputs,
    .algorithm = AlgorithmSpec{adaptFromTask<TrackMCStudy>(dataRequest, ggRequest, srcTracks, checkSV, useCCDBParams, enableCasc, enable3body)},
    .options = opts};
}

} // namespace o2::trackstudy
