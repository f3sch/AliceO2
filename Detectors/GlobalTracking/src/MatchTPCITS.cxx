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

#include <TTree.h>
#include <cassert>
#include <algorithm>

#include <fairlogger/Logger.h>
#include "Field/MagneticField.h"
#include "Field/MagFieldFast.h"
#include "ITSBase/GeometryTGeo.h"
#include "CommonUtils/TreeStream.h"

#include "DataFormatsTPC/Defs.h"
#include "TPCBase/ParameterElectronics.h"
#include "TPCBase/ParameterDetector.h"
#include "MathUtils/Cartesian.h"
#include "MathUtils/Utils.h"
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "CommonConstants/GeomConstants.h"
#include "DetectorsBase/GeometryManager.h"

#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <TFile.h>
#include <TGeoGlobalMagField.h>
#include "SimulationDataFormat/MCTruthContainer.h"
#include "TPCReconstruction/TPCFastTransformHelperO2.h"
#include "CommonUtils/NameConf.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "GlobalTracking/MatchTPCITS.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DataFormatsGlobalTracking/RecoContainerCreateTracksVariadic.h"
#include "DataFormatsGlobalTracking/TrackTuneParams.h"
#include "DataFormatsTPC/WorkflowHelper.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "ITStracking/IOUtils.h"

#ifdef ENABLE_UPGRADES
#include "ITS3Reconstruction/IOUtils.h"
#endif

#include "GPUO2Interface.h" // Needed for propper settings in GPUParam.h
#include "GPUParam.h"
#include "GPUParam.inc"
#ifdef WITH_OPENMP
#include <omp.h>
#endif

using namespace o2::globaltracking;

using MatrixDSym4 = ROOT::Math::SMatrix<double, 4, 4, ROOT::Math::MatRepSym<double, 4>>;
using MatrixD4 = ROOT::Math::SMatrix<double, 4, 4, ROOT::Math::MatRepStd<double, 4>>;
using NAMES = o2::base::NameConf;
using GTrackID = o2::dataformats::GlobalTrackID;
using TrackTunePar = o2::globaltracking::TrackTuneParams;
constexpr float MatchTPCITS::Tan70, MatchTPCITS::Cos70I2, MatchTPCITS::MaxSnp, MatchTPCITS::MaxTgp;

LinksPoolMT* TPCABSeed::gLinksPool = nullptr;

//______________________________________________
MatchTPCITS::MatchTPCITS() = default;

//______________________________________________
MatchTPCITS::~MatchTPCITS() = default;

//______________________________________________
void MatchTPCITS::run(const o2::globaltracking::RecoContainer& inp,
                      pmr::vector<o2::dataformats::TrackTPCITS>& matchedTracks,
                      pmr::vector<o2::itsmft::TrkClusRef>& ABTrackletRefs,
                      pmr::vector<int>& ABTrackletClusterIDs,
                      pmr::vector<o2::MCCompLabel>& matchLabels,
                      pmr::vector<o2::MCCompLabel>& ABTrackletLabels,
                      pmr::vector<o2::dataformats::Triplet<float, float, float>>& calib)
{
  ///< perform matching for provided input
  if (!mInitDone) {
    LOG(fatal) << "init() was not done yet";
  }
  clear();
  mRecoCont = &inp;
  mStartIR = inp.startIR;
  updateTimeDependentParams();

  mTimer[SWTot].Start(false);

  while (1) {
    if (!prepareITSData() || !prepareTPCData() || !prepareFITData()) {
      break;
    }
    if (mVDriftCalibOn) { // in the beginning of the output vector we send the full and reference VDrift used for this TF
      calib.emplace_back(mTPCVDrift, mTPCDrift.refVDrift, -999.);
      calib.emplace_back(mTPCDriftTimeOffset, mTPCDrift.refTimeOffset, -999.);
    }

    mTimer[SWDoMatching].Start(false);
    for (int sec = o2::constants::math::NSectors; sec--;) {
      doMatching(sec);
    }
    mTimer[SWDoMatching].Stop();
    if constexpr (false) { // enabling this creates very verbose output
      mTimer[SWTot].Stop();
      printCandidatesTPC();
      printCandidatesITS();
      mTimer[SWTot].Start(false);
    }

    selectBestMatches();

    bool fullMatchRefitDone = false;
    if (mUseFT0 && mParams->runAfterBurner) {
      fullMatchRefitDone = runAfterBurner(matchedTracks, matchLabels, ABTrackletLabels, ABTrackletClusterIDs, ABTrackletRefs, calib);
    }
    if (!fullMatchRefitDone) {
      refitWinners(matchedTracks, matchLabels, calib); // it afterburner is active, full matches refit will be done by it
    }
    if (mParams->verbosity > 0) {
      reportSizes(matchedTracks, ABTrackletRefs, ABTrackletClusterIDs, matchLabels, ABTrackletLabels, calib);
    }
#ifdef _ALLOW_DEBUG_TREES_
    if (mDBGOut && isDebugFlag(WinnerMatchesTree)) {
      dumpWinnerMatches();
    }
#endif
    break;
  }
  mTimer[SWTot].Stop();

  if (mParams->verbosity > 0) {
    reportTiming();
  }
  mTFCount++;
}

//______________________________________________
void MatchTPCITS::reportTiming()
{
  for (int i = 0; i < NStopWatches; i++) {
    LOGF(info, "Timing for %15s: Cpu: %.3e Real: %.3e s in %d slots of TF#%d", TimerName[i], mTimer[i].CpuTime(), mTimer[i].RealTime(), mTimer[i].Counter() - 1, mTFCount);
  }
}

//______________________________________________
void MatchTPCITS::end()
{
#ifdef _ALLOW_DEBUG_TREES_
  mDBGOut.reset();
#endif
}

//______________________________________________
void MatchTPCITS::clear()
{
  ///< clear results of previous TF reconstruction
  mMatchRecordsTPC.clear();
  mMatchRecordsITS.clear();
  mWinnerChi2Refit.clear();
  mITSWork.clear();
  mTPCWork.clear();
  mInteractions.clear();
  mITSROFTimes.clear();
  mITSTrackROFContMapping.clear();
  mITSClustersArray.clear();
  mITSClusterSizes.clear();
  mTPCABSeeds.clear();
  mTPCABIndexCache.clear();
  mABWinnersIDs.clear();
  mABClusterLinkIndex.clear();
  mNMatchesControl = 0;

  for (int sec = o2::constants::math::NSectors; sec--;) {
    mITSSectIndexCache[sec].clear();
    mITSTimeStart[sec].clear();
    mTPCSectIndexCache[sec].clear();
    mTPCTimeStart[sec].clear();
  }

  if (mMCTruthON) {
    mTPCLblWork.clear();
    mITSLblWork.clear();
  }
  for (int i = 0; i < mNThreads; i++) {
    mABLinksPool.threadPool[i].clear();
  }
}

//______________________________________________
void MatchTPCITS::setTPCVDrift(const o2::tpc::VDriftCorrFact& v)
{
  mTPCDrift = v;
  mTPCVDrift = v.getVDrift();
  mTPCDriftTimeOffset = v.getTimeOffset();
}

//______________________________________________
void MatchTPCITS::setTPCCorrMaps(o2::gpu::CorrectionMapsHelper* maph)
{
  mTPCCorrMapsHelper = maph;
}

//______________________________________________
void MatchTPCITS::init()
{
  ///< perform initizalizations, precalculate what is needed
  if (mInitDone) {
    LOG(error) << "Initialization was already done";
    return;
  }
  for (int i = NStopWatches; i--;) {
    mTimer[i].Stop();
    mTimer[i].Reset();
  }
  mParams = &Params::Instance();
  YMaxAtXMatchingRef = mParams->XMatchingRef * 0.17632698; ///< max Y in the sector at reference X
  mParams->printKeyValues();
  mFT0Params = &o2::ft0::InteractionTag::Instance();
  setUseMatCorrFlag(mParams->matCorr);
  auto* prop = o2::base::Propagator::Instance();
  if (!prop->getMatLUT() && mParams->matCorr == o2::base::Propagator::MatCorrType::USEMatCorrLUT) {
    LOG(warning) << "Requested material LUT is not loaded, switching to TGeo usage";
    setUseMatCorrFlag(o2::base::Propagator::MatCorrType::USEMatCorrTGeo);
  }

  // make sure T2GRot matrices are loaded into ITS geometry helper
  o2::its::GeometryTGeo::Instance()->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2GRot) | o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L));

  mSectEdgeMargin = mParams->crudeAbsDiffCut[o2::track::kY] / std::sqrt(Cos70I2);

#ifdef _ALLOW_DEBUG_TREES_
  // debug streamer
  if (mDBGFlags) {
    mDBGOut = std::make_unique<o2::utils::TreeStreamRedirector>(mDebugTreeFileName.data(), "recreate");
  }
#endif

  if (mParams->runAfterBurner) { // only used in AfterBurner
    mRGHelper.init();    // prepare helper for TPC track / ITS clusters matching
  }

  clear();

  mInitDone = true;

  if (fair::Logger::Logging(fair::Severity::info)) {
    print();
  }
}

//______________________________________________
void MatchTPCITS::updateTimeDependentParams()
{
  ///< update parameters depending on time (once per TF)
  auto& elParam = o2::tpc::ParameterElectronics::Instance();
  auto& detParam = o2::tpc::ParameterDetector::Instance();
  mTPCTBinMUS = elParam.ZbinWidth;
  mTPCTBinNS = mTPCTBinMUS * 1e3;
  mTPCZMax = detParam.TPClength;
  mTPCTBinMUSInv = 1. / mTPCTBinMUS;
  assert(mITSROFrameLengthMUS > 0.0f);
  mTPCBin2Z = mTPCTBinMUS * mTPCVDrift;
  mZ2TPCBin = 1. / mTPCBin2Z;
  mTPCVDriftInv = 1. / mTPCVDrift;
  mNTPCBinsFullDrift = mTPCZMax * mZ2TPCBin;
  mTPCTimeEdgeTSafeMargin = z2TPCBin(mParams->safeMarginTPCTimeEdge);
  mTPCExtConstrainedNSigmaInv = 1.f / mParams->tpcExtConstrainedNSigma;
  mBz = o2::base::Propagator::Instance()->getNominalBz();
  mFieldON = std::abs(mBz) > 0.01;

  mMinTPCTrackPtInv = (mFieldON && mParams->minTPCTrackR > 0) ? 1. / std::abs(mParams->minTPCTrackR * mBz * o2::constants::math::B2C) : 999.;
  mMinITSTrackPtInv = (mFieldON && mParams->minITSTrackR > 0) ? 1. / std::abs(mParams->minITSTrackR * mBz * o2::constants::math::B2C) : 999.;

  o2::math_utils::Point3D<float> p0(90., 1., 1), p1(90., 100., 100.);
  auto matbd = o2::base::Propagator::Instance()->getMatBudget(mParams->matCorr, p0, p1);
  mTPCmeanX0Inv = matbd.meanX2X0 / matbd.length;

  const auto& trackTune = TrackTuneParams::Instance();
  float scale = mTPCCorrMapsHelper->getInstLumiCTP();
  if (scale < 0.f) {
    scale = 0.f;
  }
  mCovDiagInner = trackTune.getCovInnerTotal(scale);
  mCovDiagOuter = trackTune.getCovOuterTotal(scale);
}

//______________________________________________
void MatchTPCITS::selectBestMatches()
{
  ///< loop over match records and select the ones with best chi2
  mTimer[SWSelectBest].Start(false);
  int nValidated = 0, iter = 0;
  mNMatches = 0;
  mNCalibPrelim = 0;
  do {
    nValidated = 0;
    int ntpc = mTPCWork.size(), nremaining = 0;
    for (int it = 0; it < ntpc; it++) {
      auto& tTPC = mTPCWork[it];
      if (isDisabledTPC(tTPC) || isValidatedTPC(tTPC)) {
        continue;
      }
      nremaining++;
      if (validateTPCMatch(it)) {
        nValidated++;
        if (mVDriftCalibOn && (!mFieldON || std::abs(tTPC.getQ2Pt()) < mParams->maxVDriftTrackQ2Pt)) {
          mNCalibPrelim++;
        }
        continue;
      }
    }
    if (mParams->verbosity > 0) {
      LOGP(info, "iter {}: Validated {} of {} remaining matches", iter, nValidated, nremaining);
    }
    iter++;
    mNMatches += nValidated;
  } while (nValidated);

  mTimer[SWSelectBest].Stop();
  LOGP(info, "Validated {} matches out of {} for {} TPC and {} ITS tracks in {} iterations", mNMatches, mNMatchesControl, mTPCWork.size(), mITSWork.size(), iter);
}

//______________________________________________
bool MatchTPCITS::validateTPCMatch(int iTPC)
{
  auto& tTPC = mTPCWork[iTPC];
  auto& rcTPC = mMatchRecordsTPC[tTPC.matchID]; // best TPC->ITS match
  // check if it is consistent with corresponding ITS->TPC match
  auto& tITS = mITSWork[rcTPC.partnerID];       //  partner ITS track
  auto& rcITS = mMatchRecordsITS[tITS.matchID]; // best ITS->TPC match record
  if (rcITS.nextRecID == Validated) {
    return false;
  }
  if (rcITS.partnerID == iTPC) {        // is best matching TPC track for this ITS track actually iTPC?
    int cloneID = tITS.getCloneShift(); // check if there is a clone of tITS
    while (cloneID) {
      cloneID += rcTPC.partnerID;
      auto& tITSClone = mITSWork[cloneID];
      if (isDisabledITS(tITSClone)) { // ignore clone
        break;
      }
      int nextITSCloneMatchID = tITSClone.matchID;
      if (rcITS.isBetter(mMatchRecordsITS[nextITSCloneMatchID])) { // best ITS->TPC match record for the clone is worse than the rcITS
        LOGP(debug, "Suppressing clone cloneID={} of winner clone {} of source ITS {}", cloneID, rcTPC.partnerID, tITS.sourceID);
        while (nextITSCloneMatchID > MinusOne) {
          auto& rcITSClone = mMatchRecordsITS[nextITSCloneMatchID];
          removeITSfromTPC(cloneID, rcITSClone.partnerID);
          nextITSCloneMatchID = rcITSClone.nextRecID;
        }
        tITSClone.matchID = MinusTen; // disable
        break;
      }
      return false; // ignore match at this iteration
    }
    // unlink winner TPC track from all ITS candidates except winning one
    int nextTPC = rcTPC.nextRecID;
    while (nextTPC > MinusOne) {
      auto& rcTPCrem = mMatchRecordsTPC[nextTPC];
      removeTPCfromITS(iTPC, rcTPCrem.partnerID); // remove references on mtID from ITS match=rcTPCrem.partnerID
      nextTPC = rcTPCrem.nextRecID;
    }
    rcTPC.nextRecID = Validated;
    int itsWinID = rcTPC.partnerID;

    // unlink winner ITS match from all TPC matches using it
    int nextITS = rcITS.nextRecID;
    while (nextITS > MinusOne) {
      auto& rcITSrem = mMatchRecordsITS[nextITS];
      removeITSfromTPC(itsWinID, rcITSrem.partnerID); // remove references on itsWinID from TPC match=rcITSrem.partnerID
      nextITS = rcITSrem.nextRecID;
    }
    rcITS.nextRecID = Validated;
    tTPC.gid.setBit(0); // Flag full match
    return true;
  }
  return false;
}

//______________________________________________
int MatchTPCITS::getNMatchRecordsTPC(const TrackLocTPC& tTPC) const
{
  ///< get number of matching records for TPC track
  int count = 0, recID = tTPC.matchID;
  while (recID > MinusOne) {
    recID = mMatchRecordsTPC[recID].nextRecID;
    count++;
  }
  return count;
}

//______________________________________________
int MatchTPCITS::getNMatchRecordsITS(const TrackLocITS& tTPC) const
{
  ///< get number of matching records for ITS track
  int count = 0, recID = tTPC.matchID;
  while (recID > MinusOne) {
    auto& itsRecord = mMatchRecordsITS[recID];
    recID = itsRecord.nextRecID;
    count++;
  }
  return count;
}

//______________________________________________
void MatchTPCITS::addTPCSeed(const o2::track::TrackParCov& _tr, float t0, float terr, GTrackID srcGID, int tpcID)
{
  // account single TPC seed, can be from standalone TPC track or constrained track from match to TRD and/or TOF
  const float SQRT12DInv = 2. / sqrt(12.);
  if (_tr.getX() > o2::constants::geom::XTPCInnerRef + 0.1 || std::abs(_tr.getQ2Pt()) > mMinTPCTrackPtInv) {
    return;
  }
  const auto& tpcOrig = mTPCTracksArray[tpcID];
  // discard tracks w/o certain number of total or innermost pads (last cluster is innermost one)
  if (tpcOrig.getNClusterReferences() < mParams->minTPCClusters) {
    return;
  }
  uint8_t clSect = 0, clRow = 0;
  uint32_t clIdx = 0;
  tpcOrig.getClusterReference(mTPCTrackClusIdx, tpcOrig.getNClusterReferences() - 1, clSect, clRow, clIdx);
  if (clRow > mParams->askMinTPCRow[clSect]) {
    return;
  }
  // create working copy of track param
  bool extConstrained = srcGID.getSource() != GTrackID::TPC;
  if (extConstrained) {
    terr *= mParams->tpcExtConstrainedNSigma;
  } else {
    terr += tpcTimeBin2MUS(tpcOrig.hasBothSidesClusters() ? mParams->safeMarginTPCITSTimeBin : mTPCTimeEdgeTSafeMargin);
  }
  auto& trc = mTPCWork.emplace_back(
    TrackLocTPC{_tr, {t0 - terr, t0 + terr}, extConstrained ? t0 : tpcTimeBin2MUS(tpcOrig.getTime0()),
                // for A/C constrained tracks the terr is half-interval, for externally constrained tracks it is sigma*Nsigma
                terr * (extConstrained ? mTPCExtConstrainedNSigmaInv : SQRT12DInv),
                tpcID,
                srcGID,
                MinusOne,
                (extConstrained || tpcOrig.hasBothSidesClusters()) ? TrackLocTPC::Constrained : (tpcOrig.hasASideClustersOnly() ? TrackLocTPC::ASide : TrackLocTPC::CSide)});
  // propagate to matching Xref
  const auto& trackTune = TrackTuneParams::Instance();
  // only TPC standalone need to be corrected on the input, provided they were not corrected at the source level,
  // other inputs are corrected in respective upstream matching processes
  if (srcGID.getSource() == GTrackID::TPC && !trackTune.sourceLevelTPC) {
    if (trackTune.useTPCInnerCorr) {
      trc.updateParams(trackTune.tpcParInner);
    }
    if (trackTune.tpcCovInnerType != TrackTuneParams::AddCovType::Disable) {
      trc.updateCov(mCovDiagInner, trackTune.tpcCovInnerType == TrackTuneParams::AddCovType::WithCorrelations);
    }
  }
  if (!propagateToRefX(trc)) {
    mTPCWork.pop_back(); // discard track whose propagation to XMatchingRef failed
    return;
  }
  if (mMCTruthON) {
    mTPCLblWork.emplace_back(mTPCTrkLabels[tpcID]);
  }
  // cache work track index
  mTPCSectIndexCache[o2::math_utils::angle2Sector(trc.getAlpha())].push_back(mTPCWork.size() - 1);
}

//______________________________________________
bool MatchTPCITS::prepareTPCData()
{
  ///< load TPC data and prepare for matching
  mTimer[SWPrepTPC].Start(false);
  const auto& inp = *mRecoCont;

  mTPCTracksArray = inp.getTPCTracks();
  mTPCTrackClusIdx = inp.getTPCTracksClusterRefs();
  mTPCClusterIdxStruct = &inp.inputsTPCclusters->clusterIndex;
  mTPCRefitterShMap = inp.clusterShMapTPC;
  mTPCRefitterOccMap = inp.occupancyMapTPC;

  if (mMCTruthON) {
    mTPCTrkLabels = inp.getTPCTracksMCLabels();
  }

  int ntr = mTPCTracksArray.size(), ntrW = 0.7 * ntr;
  mMatchRecordsTPC.reserve(mParams->maxMatchCandidates * ntrW); // number of records might be actually more than N tracks!
  mTPCWork.reserve(ntrW);
  if (mMCTruthON) {
    mTPCLblWork.reserve(ntrW);
  }
  for (int sec = o2::constants::math::NSectors; sec--;) {
    mTPCSectIndexCache[sec].reserve(100 + 1.2 * ntrW / o2::constants::math::NSectors);
  }

  auto creator = [this](auto& trk, GTrackID gid, float time0, float terr) {
    if constexpr (isITSTrack<decltype(trk)>()) {
      // do nothing, ITS tracks will be processed in a direct loop over ROFs
    }
    if constexpr (!std::is_base_of_v<o2::track::TrackParCov, std::decay_t<decltype(trk)>>) {
      return true;
    } else if (std::abs(trk.getQ2Pt()) > mMinTPCTrackPtInv) {
      return true;
    }
    if constexpr (isTPCTrack<decltype(trk)>()) {
      // unconstrained TPC track, with t0 = TrackTPC.getTime0+0.5*(DeltaFwd-DeltaBwd) and terr = 0.5*(DeltaFwd+DeltaBwd) in TimeBins
      if (!this->mSkipTPCOnly && trk.getNClusters() > 0) {
        this->addTPCSeed(trk, this->tpcTimeBin2MUS(time0), this->tpcTimeBin2MUS(terr), gid, gid.getIndex());
      }
    }
    if constexpr (isTPCTOFTrack<decltype(trk)>()) {
      // TPC track constrained by TOF time, time and its error in \mus
      this->addTPCSeed(trk, time0, terr, gid, this->mRecoCont->getTPCContributorGID(gid));
    }
    if constexpr (isTRDTrack<decltype(trk)>()) {
      // TPC track constrained by TRD trigger time, time and its error in \mus
      this->addTPCSeed(trk, time0, terr, gid, this->mRecoCont->getTPCContributorGID(gid));
    }
    // note: TPCTRDTPF tracks are actually TRD track with extra TOF cluster
    return true;
  };
  mRecoCont->createTracksVariadic(creator);

  float maxTime = 0;
  int nITSROFs = mITSROFTimes.size();
  // sort tracks in each sector according to their timeMax
  for (int sec = o2::constants::math::NSectors; sec--;) {
    auto& indexCache = mTPCSectIndexCache[sec];
    if (mParams->verbosity > 0) {
      LOG(info) << "Sorting sector" << sec << " | " << indexCache.size() << " TPC tracks";
    }
    if (!indexCache.size()) {
      continue;
    }
    std::sort(indexCache.begin(), indexCache.end(), [this](int a, int b) {
      auto& trcA = mTPCWork[a];
      auto& trcB = mTPCWork[b];
      return (trcA.tBracket.getMax() - trcB.tBracket.getMax()) < 0.;
    });

    // build array of 1st entries with tmax corresponding to each ITS ROF (or trigger),
    // TPC tracks below this entry cannot match to ITS tracks of this and higher ROFs

    float tmax = mTPCWork[indexCache.back()].tBracket.getMax();
    if (maxTime < tmax) {
      maxTime = tmax;
    }
    int nbins = 1 + (mITSTriggered ? time2ITSROFrameTrig(tmax, 0) : time2ITSROFrameCont(tmax));
    auto& timeStart = mTPCTimeStart[sec];
    timeStart.resize(nbins, -1);
    int itsROF = 0;

    timeStart[0] = 0;
    for (int itr = 0; itr < (int)indexCache.size(); itr++) {
      auto& trc = mTPCWork[indexCache[itr]];
      while (itsROF < nITSROFs && !(trc.tBracket < mITSROFTimes[itsROF])) { // 1st ITS frame after max allowed time for this TPC track
        itsROF++;
      }
      int itsROFMatch = itsROF;
      if (itsROFMatch && timeStart[--itsROFMatch] == -1) { // register ITSrof preceding the one which exceeds the TPC track tmax
        timeStart[itsROFMatch] = itr;
      }
    }
    for (int i = 1; i < nbins; i++) {
      if (timeStart[i] == -1) { // fill gaps with preceding indices
        timeStart[i] = timeStart[i - 1];
      }
    }
  } // loop over tracks of single sector

  // FIXME We probably don't need this
  /*
  // create mapping from TPC time to ITS ROFs
  if (mITSROFTimes.back() < maxTime) {
    maxTime = mITSROFTimes.back().getMax();
  }
  int nb = int(maxTime) + 1;
  mITSROFofTPCBin.resize(nb, -1);
  int itsROF = 0;
  for (int ib = 0; ib < nb; ib++) {
    while (itsROF < nITSROFs && ib < mITSROFTimes[itsROF].getMin()) {
      itsROF++;
    }
    mITSROFofTPCBin[ib] = itsROF;
  }
*/
  mTPCRefitter = std::make_unique<o2::gpu::GPUO2InterfaceRefit>(mTPCClusterIdxStruct, mTPCCorrMapsHelper, mBz, mTPCTrackClusIdx.data(), 0, mTPCRefitterShMap.data(), mTPCRefitterOccMap.data(), mTPCRefitterOccMap.size(), nullptr, o2::base::Propagator::Instance());
  mNTPCOccBinLength = mTPCRefitter->getParam()->rec.tpc.occupancyMapTimeBins;
  mTBinClOcc.clear();
  if (mNTPCOccBinLength > 1 && mTPCRefitterOccMap.size()) {
    mNTPCOccBinLengthInv = 1. / mNTPCOccBinLength;
    int nTPCBins = mNHBPerTF * o2::constants::lhc::LHCMaxBunches / 8, ninteg = 0;
    int nTPCOccBins = nTPCBins * mNTPCOccBinLengthInv, sumBins = std::max(1, int(o2::constants::lhc::LHCMaxBunches / 8 * mNTPCOccBinLengthInv));
    mTBinClOcc.resize(nTPCOccBins);
    float sm = 0., tb = (nTPCOccBins - 0.5) * mNTPCOccBinLength, mltPrev = 0.;
    for (int i = nTPCOccBins; i--;) {
      float mlt = mTPCRefitter->getParam()->GetUnscaledMult(tb);
      sm += mlt;
      mTBinClOcc.push_back(sm);
      if (ninteg++ > mNTPCOccBinLength) {
        sm -= mltPrev;
      }
      //      LOGP(info, "BIN {} of {} -> {} with inst val {} (prev = {}) BL={} nInt={} tb={}", i, nTPCOccBins, sm, mlt, mltPrev, mNTPCOccBinLength, ninteg, tb);
      mltPrev = mlt;
      tb -= mNTPCOccBinLength;
    }
  } else {
    mTBinClOcc.resize(1);
  }
  mInteractionMUSLUT.clear();
  mInteractionMUSLUT.resize(maxTime + 3 * o2::constants::lhc::LHCOrbitMUS, -1);
  mTimer[SWPrepTPC].Stop();
  return mTPCWork.size() > 0;
}

//_____________________________________________________
bool MatchTPCITS::prepareITSData()
{
  static size_t errCount = 0;
  constexpr size_t MaxErrors2Report = 10;
  // Do preparatory work for matching
  mTimer[SWPrepITS].Start(false);
  const auto& inp = *mRecoCont;

  // ITS clusters
  mITSClusterROFRec = inp.getITSClustersROFRecords();
  const auto clusITS = inp.getITSClusters();
  if (mITSClusterROFRec.empty() || clusITS.empty()) {
    LOG(info) << "No ITS clusters";
    return false;
  }
  const auto patterns = inp.getITSClustersPatterns();
  auto pattIt = patterns.begin();
  mITSClustersArray.reserve(clusITS.size());
#ifdef ENABLE_UPGRADES
  if (mWithIT3) {
    o2::its3::ioutils::convertCompactClusters(clusITS, pattIt, mITSClustersArray, mIT3Dict);
  } else {
    o2::its::ioutils::convertCompactClusters(clusITS, pattIt, mITSClustersArray, mITSDict);
  }
#else
  o2::its::ioutils::convertCompactClusters(clusITS, pattIt, mITSClustersArray, mITSDict);
#endif

  // ITS clusters sizes
  mITSClusterSizes.reserve(clusITS.size());
  auto pattIt2 = patterns.begin();
  for (auto& clus : clusITS) {
    auto pattID = clus.getPatternID();
    unsigned int npix;
#ifdef ENABLE_UPGRADES
    if ((pattID == o2::itsmft::CompCluster::InvalidPatternID) || ((mWithIT3) ? mIT3Dict->isGroup(pattID) : mITSDict->isGroup(pattID))) { // braces guarantee evaluation order
#else
    if (pattID == o2::itsmft::CompCluster::InvalidPatternID || mITSDict->isGroup(pattID)) {
#endif
      o2::itsmft::ClusterPattern patt;
      patt.acquirePattern(pattIt2);
      npix = patt.getNPixels();
    } else {
#ifdef ENABLE_UPGRADES
      if (mWithIT3) {
        npix = mIT3Dict->getNpixels(pattID);
      } else {
        npix = mITSDict->getNpixels(pattID);
      }
#else
      npix = mITSDict->getNpixels(pattID);
#endif
    }
    mITSClusterSizes.push_back(std::clamp(npix, 0u, 255u));
  }

  if (mMCTruthON) {
    mITSClsLabels = inp.mcITSClusters.get();
  }

  // ITS tracks
  mITSTracksArray = inp.getITSTracks();
  mITSTrackClusIdx = inp.getITSTracksClusterRefs();
  mITSTrackROFRec = inp.getITSTracksROFRecords();
  if (mMCTruthON) {
    mITSTrkLabels = inp.getITSTracksMCLabels();
  }
  int nROFs = mITSTrackROFRec.size();
  mITSWork.reserve(mITSTracksArray.size());

  // total N ITS clusters in TF
  const auto& lastClROF = mITSClusterROFRec[nROFs - 1];
  int nITSClus = lastClROF.getFirstEntry() + lastClROF.getNEntries();
  mABClusterLinkIndex.resize(nITSClus, MinusOne);
  for (int sec = o2::constants::math::NSectors; sec--;) {
    mITSTimeStart[sec].resize(nROFs, -1); // start of ITS work tracks in every sector
  }
  long nHBF = o2::base::GRPGeomHelper::getNHBFPerTF();
  long maxBCs = nHBF * long(o2::constants::lhc::LHCMaxBunches);
  o2::track::TrackLTIntegral trackLTInt;
  trackLTInt.setTimeNotNeeded();

  for (int irof = 0; irof < nROFs; irof++) {
    const auto& rofRec = mITSTrackROFRec[irof];
    long nBC = rofRec.getBCData().differenceInBC(mStartIR);
    if (nBC > maxBCs || nBC < 0) {
      if (++errCount < MaxErrors2Report) {
        LOGP(alarm, "ITS ROF#{} start {} is not compatible with TF 1st orbit {} or TF length of {} HBFs",
             irof, rofRec.getBCData().asString(), mStartIR.asString(), nHBF);
      }
      break;
    }
    float tMin = nBC * o2::constants::lhc::LHCBunchSpacingMUS + mITSTimeBiasMUS;
    float tMax = (nBC + mITSROFrameLengthInBC) * o2::constants::lhc::LHCBunchSpacingMUS + mITSTimeBiasMUS;
    if (!mITSTriggered) {
      size_t irofCont = nBC / mITSROFrameLengthInBC;
      if (mITSTrackROFContMapping.size() <= irofCont) { // there might be gaps in the non-empty rofs, this will map continuous ROFs index to non empty ones
        mITSTrackROFContMapping.resize((1 + irofCont / 128) * 128, 0);
      }
      mITSTrackROFContMapping[irofCont] = irof;
    }

    mITSROFTimes.emplace_back(tMin, tMax); // ITS ROF min/max time

    for (int sec = o2::constants::math::NSectors; sec--;) {      // start of sector's tracks for this ROF
      mITSTimeStart[sec][irof] = mITSSectIndexCache[sec].size(); // The sorting does not affect this
    }

    int trlim = rofRec.getFirstEntry() + rofRec.getNEntries();
    for (int it = rofRec.getFirstEntry(); it < trlim; it++) {
      const auto& trcOrig = mITSTracksArray[it];
      if (mParams->runAfterBurner) {
        flagUsedITSClusters(trcOrig);
      }
      if (trcOrig.getParamOut().getX() < 1.) {
        continue; // backward refit failed
      }
      if (std::abs(trcOrig.getQ2Pt()) > mMinITSTrackPtInv) {
        continue;
      }
      int nWorkTracks = mITSWork.size();
      // working copy of outer track param
      auto& trc = mITSWork.emplace_back(TrackLocITS{trcOrig.getParamOut(), {tMin, tMax}, it, irof, MinusOne});
      if (!trc.rotate(o2::math_utils::angle2Alpha(trc.getPhiPos()))) {
        mITSWork.pop_back(); // discard failed track
        continue;
      }
      // make sure the track is at the ref. radius
      trackLTInt.clearFast();
      if (!propagateToRefX(trc, &trackLTInt)) {
        mITSWork.pop_back(); // discard failed track
        continue;            // add to cache only those ITS tracks which reached ref.X and have reasonable snp
      }
      trc.xrho = trackLTInt.getXRho(); // we collect seen x*rho and distance to the reference X for further PID correcrions
      trc.dL = trackLTInt.getL();

      if (mMCTruthON) {
        mITSLblWork.emplace_back(mITSTrkLabels[it]);
      }
      trc.setUserField(0);
      // cache work track index
      int sector = o2::math_utils::angle2Sector(trc.getAlpha());
      mITSSectIndexCache[sector].push_back(nWorkTracks);

      // If the ITS track is very close to the sector edge, it may match also to a TPC track in the neighb. sector.
      // For a track with Yr and Phir at Xr the distance^2 between the poisition of this track in the neighb. sector
      // when propagated to Xr (in this neighbouring sector) and the edge will be (neglecting the curvature)
      // [(Xr*tg(10)-Yr)/(tgPhir+tg70)]^2  / cos(70)^2  // for the next sector
      // [(Xr*tg(10)+Yr)/(tgPhir-tg70)]^2  / cos(70)^2  // for the prev sector
      // Distances to the sector edges in neighbourings sectors (at Xref in their proper frames)
      float trcY = trc.getY(), tgp = trc.getSnp();
      tgp /= std::sqrt((1.f - tgp) * (1.f + tgp)); // tan of track direction XY

      float dyUpDn[2] = {std::abs((YMaxAtXMatchingRef - trcY) / (tgp + Tan70)), std::abs((YMaxAtXMatchingRef + trcY) / (tgp - Tan70))}; // sector up, down edge distances
      // we do the cloning for closest edge only
      int sel = dyUpDn[0] < dyUpDn[1] ? 0 : 1;
      if (dyUpDn[sel] < mSectEdgeMargin) { // need to check this track for matching in sector up or down
        int sectNeib = sel == 0 ? (sector < (o2::constants::math::NSectors - 1) ? sector + 1 : 0) : (sector > 1 ? sector - 1 : o2::constants::math::NSectors - 1);
        addLastTrackCloneForNeighbourSector(sectNeib, &trackLTInt);
      }
    }
  }

  if (!mITSTriggered) { // fill the gaps;
    int nr = mITSTrackROFContMapping.size();
    for (int i = 1; i < nr; i++) {
      if (mITSTrackROFContMapping[i] < mITSTrackROFContMapping[i - 1]) {
        mITSTrackROFContMapping[i] = mITSTrackROFContMapping[i - 1];
      }
    }
  }

  // sort tracks in each sector according to their min time, then tgl
  // RSTODO: sorting in tgl will be dangerous once the tracks with different time uncertaincies will be added
  for (int sec = o2::constants::math::NSectors; sec--;) {
    auto& indexCache = mITSSectIndexCache[sec];
    if (mParams->verbosity > 0) {
      LOG(info) << "Sorting sector" << sec << " | " << indexCache.size() << " ITS tracks";
    }
    if (!indexCache.size()) {
      continue;
    }
    std::sort(indexCache.begin(), indexCache.end(), [this](int a, int b) {
      auto& trackA = mITSWork[a];
      auto& trackB = mITSWork[b];
      if (trackA.tBracket.getMin() < trackB.tBracket.getMin()) {
        return true;
      } else if (trackA.tBracket.getMin() > trackB.tBracket.getMin()) {
        return false;
      }
      return trackA.getTgl() < trackB.getTgl();
    });
  } // loop over tracks of single sector
  mMatchRecordsITS.reserve(mITSWork.size() * mParams->maxMatchCandidates);
  mTimer[SWPrepITS].Stop();

  return nITSClus > 0;
}

//_____________________________________________________
bool MatchTPCITS::prepareFITData()
{
  // If available, read FIT Info
  if (mUseFT0) {
    mFITInfo = mRecoCont->getFT0RecPoints();
    prepareInteractionTimes();
  }
  return true;
}

//_____________________________________________________
void MatchTPCITS::doMatching(int sec)
{
  ///< run matching for currently cached ITS data for given TPC sector
  auto& cacheITS = mITSSectIndexCache[sec]; // array of cached ITS track indices for this sector
  auto& cacheTPC = mTPCSectIndexCache[sec]; // array of cached ITS track indices for this sector
  auto& timeStartTPC = mTPCTimeStart[sec];  // array of 1st TPC track with timeMax in ITS ROFrame
  auto& timeStartITS = mITSTimeStart[sec];
  int nTracksTPC = cacheTPC.size(), nTracksITS = cacheITS.size();
  if (!nTracksTPC || !nTracksITS) {
    if (mParams->verbosity > 0) {
      LOG(info) << "Matchng sector " << sec << " : N tracks TPC:" << nTracksTPC << " ITS:" << nTracksITS << " in sector " << sec;
    }
    return;
  }

  /// full drift time + safety margin
  float maxTDriftSafe = tpcTimeBin2MUS(mNTPCBinsFullDrift + mParams->safeMarginTPCITSTimeBin + mTPCTimeEdgeTSafeMargin);
  float vdErrT = tpcTimeBin2MUS(mZ2TPCBin * mParams->maxVDriftUncertainty);

  // get min ROFrame of ITS tracks currently in cache
  auto minROFITS = mITSWork[cacheITS.front()].roFrame;

  if (minROFITS >= int(timeStartTPC.size())) {
    LOG(info) << "ITS min ROFrame " << minROFITS << " exceeds all cached TPC track ROF eqiuvalent " << cacheTPC.size() - 1;
    return;
  }

  int nCheckTPCControl = 0, nCheckITSControl = 0, nMatchesControl = 0; // temporary
  int idxMinTPC = timeStartTPC[minROFITS];                             // index of 1st cached TPC track within cached ITS ROFrames
  auto t2nbs = tpcTimeBin2MUS(mZ2TPCBin * mParams->tpcTimeICMatchingNSigma);
  bool checkInteractionCandidates = mUseFT0 && mParams->validateMatchByFIT != MatchTPCITSParams::Disable;

  int itsROBin = 0;
  for (int itpc = idxMinTPC; itpc < nTracksTPC; itpc++) {
    auto& trefTPC = mTPCWork[cacheTPC[itpc]];
    // estimate ITS 1st ROframe bin this track may match to: TPC track are sorted according to their
    // timeMax, hence the timeMax - MaxmNTPCBinsFullDrift are non-decreasing
    auto tmn = trefTPC.tBracket.getMax() - maxTDriftSafe;
    itsROBin = mITSTriggered ? time2ITSROFrameTrig(tmn, itsROBin) : time2ITSROFrameCont(tmn);

    if (itsROBin >= int(timeStartITS.size())) { // time of TPC track exceeds the max time of ITS in the cache
      break;
    }
    int iits0 = timeStartITS[itsROBin];
    nCheckTPCControl++;
    for (auto iits = iits0; iits < nTracksITS; iits++) {
      auto& trefITS = mITSWork[cacheITS[iits]];
      // compare if the ITS and TPC tracks may overlap in time
      LOG(debug) << "TPC bracket: " << trefTPC.tBracket.asString() << " ITS bracket: " << trefITS.tBracket.asString() << " TPCtgl: " << trefTPC.getTgl() << " ITStgl: " << trefITS.getTgl();
      if (trefTPC.tBracket < trefITS.tBracket) { // since TPC tracks are sorted in timeMax and ITS tracks are sorted in timeMin all following ITS tracks also will not match
        break;
      }
      if (trefTPC.tBracket > trefITS.tBracket) { // its bracket precedes TPC bracket
        continue;
      }

      // is corrected TPC track time compatible with ITS ROF expressed
      auto deltaT = (trefITS.getZ() - trefTPC.getZ()) * mTPCVDriftInv;                                                    // drift time difference corresponding to Z differences
      auto timeCorr = trefTPC.getCorrectedTime(deltaT);                                                                   // TPC time required to match to Z of ITS track
      auto timeCorrErr = std::sqrt(trefITS.getSigmaZ2() + trefTPC.getSigmaZ2()) * t2nbs + mParams->safeMarginTimeCorrErr; // nsigma*error
      if (mVDriftCalibOn) {
        timeCorrErr += vdErrT * (250. - abs(trefITS.getZ())); // account for the extra error from TPC VDrift uncertainty
      }
      o2::math_utils::Bracketf_t trange(timeCorr - timeCorrErr, timeCorr + timeCorrErr);
      LOG(debug) << "TPC range: " << trange.asString() << " ITS bracket: " << trefITS.tBracket.asString() << " DZ: " << (trefITS.getZ() - trefTPC.getZ()) << " DT: " << timeCorr;
      // check if the assigned time is strictly within the ITS bracket
      auto cmpITSBracket = trefITS.tBracket.isOutside(timeCorr);
      if (cmpITSBracket) { // no, check if brackets are overlapping at all
        if (trefITS.tBracket.isOutside(trange)) {
          continue;
        }
        if (mParams->ITSTimeOutliersPolicy == MatchTPCITSParams::TimeOutliersPolicy::Adjust) {
          if (cmpITSBracket == o2::math_utils::Bracketf_t::Below) {
            timeCorr = trefITS.tBracket.getMin();
            trange.setMin(timeCorr);
          } else {
            timeCorr = trefITS.tBracket.getMax();
            trange.setMax(timeCorr);
          }
        } else if (mParams->ITSTimeOutliersPolicy == MatchTPCITSParams::TimeOutliersPolicy::Reject) {
          continue;
        }
      }

      nCheckITSControl++;
      float chi2 = -1;
      int rejFlag = compareTPCITSTracks(trefITS, trefTPC, chi2);

#ifdef _ALLOW_DEBUG_TREES_
      if (mDBGOut && ((rejFlag == Accept && isDebugFlag(MatchTreeAccOnly)) || isDebugFlag(MatchTreeAll))) {
        fillTPCITSmatchTree(cacheITS[iits], cacheTPC[itpc], rejFlag, chi2, timeCorr);
      }
#endif
      /*
      // RS: this might be dangerous for ITS tracks with different time coverages.
      if (rejFlag == RejectOnTgl) {
        // ITS tracks in each ROFrame are ordered in Tgl, hence if this check failed on Tgl check
        // (i.e. tgl_its>tgl_tpc+tolerance), then all other ITS tracks in this ROFrame will also have tgl too large.
        // Jump on the 1st ITS track of the next ROFrame
        int rof = trefITS.roFrame;
        bool stop = false;
        do {
          if (++rof >= int(timeStartITS.size())) {
            stop = true;
            break; // no more ITS ROFrames in cache
          }
          iits = timeStartITS[rof] - 1;                  // next track to be checked -1
        } while (iits <= timeStartITS[trefITS.roFrame]); // skip empty bins
        if (stop) {
          break;
        }
        continue;
      }
      */
      if (rejFlag != Accept) {
        continue;
      }
      int matchedIC = MinusOne;
      if (!isCosmics()) {
        // validate by bunch filling scheme
        if (mUseBCFilling) {
          auto irBracket = tBracket2IRBracket(trange);
          if (irBracket.isInvalid()) {
            continue;
          }
        }
        if (checkInteractionCandidates && mInteractions.size()) {
          // check if corrected TPC track time is compatible with any of interaction times
          int tmus = trange.getMin();
          if (tmus < 0) {
            tmus = 0;
          }
          auto entStart = tmus < int(mInteractionMUSLUT.size()) ? mInteractionMUSLUT[tmus] : (mInteractionMUSLUT.size() ? mInteractionMUSLUT.back() : 0);
          for (int ent = entStart; ent < (int)mInteractions.size(); ent++) {
            auto cmp = mInteractions[ent].tBracket.isOutside(trange);
            if (cmp == o2::math_utils::Bracketf_t::Above) { // trange is above this interaction candidate, the following ones may match
              continue;
            }
            if (cmp == o2::math_utils::Bracketf_t::Inside) {
              matchedIC = ent;
            }
            break; // we loop till 1st matching IC or the one above the trange (since IC are ordered, all others will be above too)
          }
        }
        if (mParams->validateMatchByFIT == MatchTPCITSParams::Require && matchedIC == MinusOne) {
          continue;
        }
      }
      registerMatchRecordTPC(cacheITS[iits], cacheTPC[itpc], chi2, matchedIC); // register matching candidate
      nMatchesControl++;
    }
  }
  if (mParams->verbosity > 0) {
    LOG(info) << "Match sector " << sec << " N tracks TPC:" << nTracksTPC << " ITS:" << nTracksITS
              << " N TPC tracks checked: " << nCheckTPCControl << " (starting from " << idxMinTPC
              << "), checks: " << nCheckITSControl << ", matches:" << nMatchesControl;
  }
  mNMatchesControl += nMatchesControl;
}

//______________________________________________
void MatchTPCITS::suppressMatchRecordITS(int itsID, int tpcID)
{
  ///< suppress the reference on the tpcID in the list of matches recorded for itsID
  auto& tITS = mITSWork[itsID];
  int topID = MinusOne, recordID = tITS.matchID; // 1st entry in mMatchRecordsITS
  while (recordID > MinusOne) {                  // navigate over records for given ITS track
    if (mMatchRecordsITS[recordID].partnerID == tpcID) {
      // unlink this record, connecting its child to its parrent
      if (topID < 0) {
        tITS.matchID = mMatchRecordsITS[recordID].nextRecID;
      } else {
        mMatchRecordsITS[topID].nextRecID = mMatchRecordsITS[recordID].nextRecID;
      }
      return;
    }
    topID = recordID;
    recordID = mMatchRecordsITS[recordID].nextRecID; // check next record
  }
}

//______________________________________________
bool MatchTPCITS::registerMatchRecordTPC(int iITS, int iTPC, float chi2, int candIC)
{
  ///< record matching candidate, making sure that number of ITS candidates per TPC track, sorted
  ///< in matching chi2 does not exceed allowed number
  auto& tTPC = mTPCWork[iTPC];                                   // get MatchRecord structure of this TPC track, create if none
  if (tTPC.matchID < 0) {                                        // no matches yet, just add new record
    registerMatchRecordITS(iITS, iTPC, chi2, candIC);            // register TPC track in the ITS records
    tTPC.matchID = mMatchRecordsTPC.size();                      // new record will be added in the end
    mMatchRecordsTPC.emplace_back(iITS, chi2, MinusOne, candIC); // create new record with empty reference on next match
    return true;
  }

  int count = 0, nextID = tTPC.matchID, topID = MinusOne;
  do {
    auto& nextMatchRec = mMatchRecordsTPC[nextID];
    count++;
    if (!nextMatchRec.isBetter(chi2, candIC)) { // need to insert new record before nextMatchRec?
      if (count < mParams->maxMatchCandidates) {
        break;                                                // will insert in front of nextID
      } else {                                                // max number of candidates reached, will overwrite the last one
        suppressMatchRecordITS(nextMatchRec.partnerID, iTPC); // flag as disabled the overriden ITS match
        registerMatchRecordITS(iITS, iTPC, chi2, candIC);     // register TPC track entry in the ITS records
        // reuse the record of suppressed ITS match to store better one
        nextMatchRec.chi2 = chi2;
        nextMatchRec.matchedIC = candIC;
        nextMatchRec.partnerID = iITS;
        return true;
      }
    }
    topID = nextID; // check next match record
    nextID = nextMatchRec.nextRecID;
  } while (nextID > MinusOne);

  // if count == mParams->maxMatchCandidates, the max number of candidates was already reached, and the
  // new candidated was either discarded (if its chi2 is worst one) or has overwritten worst
  // existing candidate. Otherwise, we need to add new entry
  if (count < mParams->maxMatchCandidates) {
    if (topID < 0) {                                                       // the new match is top candidate
      topID = tTPC.matchID = mMatchRecordsTPC.size();                      // register new record as top one
    } else {                                                               // there are better candidates
      topID = mMatchRecordsTPC[topID].nextRecID = mMatchRecordsTPC.size(); // register to his parent
    }
    // nextID==-1 will mean that the while loop run over all candidates->the new one is the worst (goes to the end)
    registerMatchRecordITS(iITS, iTPC, chi2, candIC);          // register TPC track in the ITS records
    mMatchRecordsTPC.emplace_back(iITS, chi2, nextID, candIC); // create new record with empty reference on next match
    // make sure that after addition the number of candidates don't exceed allowed number
    count++;
    while (nextID > MinusOne) {
      if (count > mParams->maxMatchCandidates) {
        suppressMatchRecordITS(mMatchRecordsTPC[nextID].partnerID, iTPC);
        // exclude nextID record, w/o changing topID (which becomes the last record)
        nextID = mMatchRecordsTPC[topID].nextRecID = mMatchRecordsTPC[nextID].nextRecID;
        continue;
      }
      count++;
      topID = nextID;
      nextID = mMatchRecordsTPC[nextID].nextRecID;
    }
    return true;
  } else {
    return false; // unless nextID was assigned OverrideExisting, new candidate was discarded
  }
}

//______________________________________________
void MatchTPCITS::registerMatchRecordITS(const int iITS, int iTPC, float chi2, int candIC)
{
  ///< register TPC match in ITS tracks match records, ordering them in quality
  auto& tITS = mITSWork[iITS];
  int idnew = mMatchRecordsITS.size();
  auto& newRecord = mMatchRecordsITS.emplace_back(iTPC, chi2, MinusOne, candIC); // associate iTPC with this record
  if (tITS.matchID < 0) {
    tITS.matchID = idnew;
    return;
  }
  // there are other matches for this ITS track, insert the new record preserving quality order
  // navigate till last record or the one with worse chi2
  int topID = MinusOne, nextRecord = tITS.matchID;
  do {
    auto& nextMatchRec = mMatchRecordsITS[nextRecord];
    if (!nextMatchRec.isBetter(chi2, candIC)) { // need to insert new record before nextMatchRec?
      newRecord.nextRecID = nextRecord;         // new one will refer to old one it overtook
      if (topID < 0) {
        tITS.matchID = idnew; // the new one is the best match, track will refer to it
      } else {
        mMatchRecordsITS[topID].nextRecID = idnew; // new record will follow existing better one
      }
      return;
    }
    topID = nextRecord;
    nextRecord = mMatchRecordsITS[nextRecord].nextRecID;
  } while (nextRecord > MinusOne);

  // if we reached here, the new record should be added in the end
  mMatchRecordsITS[topID].nextRecID = idnew; // register new link
}

//______________________________________________
int MatchTPCITS::compareTPCITSTracks(const TrackLocITS& tITS, const TrackLocTPC& tTPC, float& chi2) const
{
  ///< compare pair of ITS and TPC tracks
  chi2 = -1.f;
  int rejFlag = Accept;
  float diff; // make rough check differences and their nsigmas

  // start with check on Tgl, since rjection on it will allow to profit from sorting
  diff = tITS.getParam(o2::track::kTgl) - tTPC.getParam(o2::track::kTgl);
  if ((rejFlag = roughCheckDif(diff, mParams->crudeAbsDiffCut[o2::track::kTgl], RejectOnTgl))) {
    return rejFlag;
  }
  auto err2Tgl = tITS.getDiagError2(o2::track::kTgl) + tTPC.getDiagError2(o2::track::kTgl);
  if (mVDriftCalibOn) {
    auto addErr = tITS.getParam(o2::track::kTgl) * mParams->maxVDriftUncertainty;
    err2Tgl += addErr * addErr; // account for VDrift uncertainty
  }
  diff *= diff / err2Tgl;
  if ((rejFlag = roughCheckDif(diff, mParams->crudeNSigma2Cut[o2::track::kTgl], RejectOnTgl + NSigmaShift))) {
    return rejFlag;
  }
  // do we need to account for different PID hypotheses used for ITS and TPC tracks propagation to ref. X?
  bool testOtherPID = false;
  float itsParam[5] = {tITS.getY(), tITS.getZ(), tITS.getSnp(), tITS.getTgl(), tITS.getQ2Pt()};
  if (tTPC.getPID() > tITS.getPID() && tITS.dL > 0.f && tTPC.getP2() / tTPC.getPID().getMass2() < mParams->minBetaGammaForPIDDiff) {
    o2::track::TrackPar tPID(mITSTracksArray[tITS.sourceID].getParamOut()); // clone original ITS track at highest update point
    tPID.setPID(tTPC.getPID(), true);
    if (!tPID.correctForELoss(tITS.xrho)) {
      return RejectoOnPIDCorr;
    }
    float dCurv = (tPID.getQ2Pt() - tITS.getQ2Pt()) * mBz * o2::constants::math::B2C, dLEff = tITS.dL * mParams->ITSStepEffFraction, dCurvL = dCurv * dLEff;
    itsParam[o2::track::kQ2Pt] = tPID.getQ2Pt();
    itsParam[o2::track::kSnp] += dCurvL;
    if (std::abs(itsParam[o2::track::kSnp]) >= 1.) {
      itsParam[o2::track::kSnp] = std::copysign(0.99, itsParam[o2::track::kSnp]);
    }
    itsParam[o2::track::kY] += dCurvL * dLEff * 0.5;
    testOtherPID = true;
  }

  diff = itsParam[o2::track::kY] - tTPC.getParam(o2::track::kY);
  if ((rejFlag = roughCheckDif(diff, mParams->crudeAbsDiffCut[o2::track::kY], RejectOnY))) {
    return rejFlag;
  }
  diff *= diff / (tITS.getDiagError2(o2::track::kY) + tTPC.getDiagError2(o2::track::kY));
  if ((rejFlag = roughCheckDif(diff, mParams->crudeNSigma2Cut[o2::track::kY], RejectOnY + NSigmaShift))) {
    return rejFlag;
  }

  if (tTPC.constraint == TrackLocTPC::Constrained) { // in continuous only constrained tracks can be compared in Z
    diff = itsParam[o2::track::kZ] - tTPC.getParam(o2::track::kZ);
    if ((rejFlag = roughCheckDif(diff, mParams->crudeAbsDiffCut[o2::track::kZ], RejectOnZ))) {
      return rejFlag;
    }
    diff *= diff / (tITS.getDiagError2(o2::track::kZ) + tTPC.getDiagError2(o2::track::kZ));
    if ((rejFlag = roughCheckDif(diff, mParams->crudeNSigma2Cut[o2::track::kZ], RejectOnZ + NSigmaShift))) {
      return rejFlag;
    }
  }

  diff = itsParam[o2::track::kSnp] - tTPC.getParam(o2::track::kSnp);
  if ((rejFlag = roughCheckDif(diff, mParams->crudeAbsDiffCut[o2::track::kSnp], RejectOnSnp))) {
    return rejFlag;
  }
  diff *= diff / (tITS.getDiagError2(o2::track::kSnp) + tTPC.getDiagError2(o2::track::kSnp));
  if ((rejFlag = roughCheckDif(diff, mParams->crudeNSigma2Cut[o2::track::kSnp], RejectOnSnp + NSigmaShift))) {
    return rejFlag;
  }

  diff = itsParam[o2::track::kQ2Pt] - tTPC.getParam(o2::track::kQ2Pt);
  if ((rejFlag = roughCheckDif(diff, mParams->crudeAbsDiffCut[o2::track::kQ2Pt], RejectOnQ2Pt))) {
    return rejFlag;
  }
  diff *= diff / (tITS.getDiagError2(o2::track::kQ2Pt) + tTPC.getDiagError2(o2::track::kQ2Pt));
  if ((rejFlag = roughCheckDif(diff, mParams->crudeNSigma2Cut[o2::track::kQ2Pt], RejectOnQ2Pt + NSigmaShift))) {
    return rejFlag;
  }
  // calculate mutual chi2 excluding Z in continuous mode
  if (testOtherPID) { // temporarily substitute pion params by alternative ones
    auto tITSAlt = tITS;
    tITSAlt.setPID(tTPC.getPID());
    tITSAlt.setParam(itsParam[o2::track::kY], o2::track::kY);
    tITSAlt.setParam(itsParam[o2::track::kSnp], o2::track::kSnp);
    tITSAlt.setParam(itsParam[o2::track::kQ2Pt], o2::track::kQ2Pt);
    chi2 = getPredictedChi2NoZ(tITSAlt, tTPC);
  } else {
    chi2 = getPredictedChi2NoZ(tITS, tTPC);
  }
  if (chi2 > mParams->cutMatchingChi2 || chi2 < 0.) { // sometimes due to the numerical stability the chi2 is negative, reject it.
    return RejectOnChi2;
  }

  return Accept;
}

//______________________________________________
void MatchTPCITS::printCandidatesTPC() const
{
  ///< print mathing records
  int ntpc = mTPCWork.size();
  printf("\n\nPrinting all TPC -> ITS matches for %d TPC tracks\n", ntpc);
  for (int i = 0; i < ntpc; i++) {
    const auto& tTPC = mTPCWork[i];
    int nm = getNMatchRecordsTPC(tTPC);
    printf("*** trackTPC#%6d %6d : Ncand = %d Best = %d\n", i, tTPC.sourceID, nm, tTPC.matchID);
    int count = 0, recID = tTPC.matchID;
    while (recID > MinusOne) {
      const auto& rcTPC = mMatchRecordsTPC[recID];
      const auto& tITS = mITSWork[rcTPC.partnerID];
      printf("  * cand %2d : ITS track %6d(src:%6d) Chi2: %.2f\n", count, rcTPC.partnerID, tITS.sourceID, rcTPC.chi2);
      count++;
      recID = rcTPC.nextRecID;
    }
  }
}

//______________________________________________
void MatchTPCITS::printCandidatesITS() const
{
  ///< print mathing records
  int nits = mITSWork.size();
  printf("\n\nPrinting all ITS -> TPC matches for %d ITS tracks\n", nits);

  for (int i = 0; i < nits; i++) {
    const auto& tITS = mITSWork[i];
    printf("*** trackITS#%6d %6d : Ncand = %d Best = %d\n", i, tITS.sourceID, getNMatchRecordsITS(tITS), tITS.matchID);
    int count = 0, recID = tITS.matchID;
    while (recID > MinusOne) {
      const auto& rcITS = mMatchRecordsITS[recID];
      const auto& tTPC = mTPCWork[rcITS.partnerID];
      printf("  * cand %2d : TPC track %6d(src:%6d) Chi2: %.2f\n", count, rcITS.partnerID, tTPC.sourceID, rcITS.chi2);
      count++;
      recID = rcITS.nextRecID;
    }
  }
}

//______________________________________________
float MatchTPCITS::getPredictedChi2NoZ(const o2::track::TrackParCov& trITS, const o2::track::TrackParCov& trTPC) const
{
  /// get chi2 between 2 tracks, neglecting Z parameter.
  /// 2 tracks must be defined at the same parameters X,alpha (check is currently commented)

  //  if (std::abs(trITS.getAlpha() - trTPC.getAlpha()) > FLT_EPSILON) {
  //    LOG(error) << "The reference Alpha of the tracks differ: "
  //        << trITS.getAlpha() << " : " << trTPC.getAlpha();
  //    return 2. * o2::track::HugeF;
  //  }
  //  if (std::abs(trITS.getX() - trTPC.getX()) > FLT_EPSILON) {
  //    LOG(error) << "The reference X of the tracks differ: "
  //        << trITS.getX() << " : " << trTPC.getX();
  //    return 2. * o2::track::HugeF;
  //  }
  MatrixDSym4 covMat;
  covMat(0, 0) = static_cast<double>(trITS.getSigmaY2()) + static_cast<double>(trTPC.getSigmaY2());
  covMat(1, 0) = static_cast<double>(trITS.getSigmaSnpY()) + static_cast<double>(trTPC.getSigmaSnpY());
  covMat(1, 1) = static_cast<double>(trITS.getSigmaSnp2()) + static_cast<double>(trTPC.getSigmaSnp2());
  covMat(2, 0) = static_cast<double>(trITS.getSigmaTglY()) + static_cast<double>(trTPC.getSigmaTglY());
  covMat(2, 1) = static_cast<double>(trITS.getSigmaTglSnp()) + static_cast<double>(trTPC.getSigmaTglSnp());
  covMat(2, 2) = static_cast<double>(trITS.getSigmaTgl2()) + static_cast<double>(trTPC.getSigmaTgl2());
  if (mVDriftCalibOn) {
    auto addErr = trITS.getParam(o2::track::kTgl) * mParams->maxVDriftUncertainty;
    covMat(2, 2) += addErr * addErr;
  }
  covMat(3, 0) = static_cast<double>(trITS.getSigma1PtY()) + static_cast<double>(trTPC.getSigma1PtY());
  covMat(3, 1) = static_cast<double>(trITS.getSigma1PtSnp()) + static_cast<double>(trTPC.getSigma1PtSnp());
  covMat(3, 2) = static_cast<double>(trITS.getSigma1PtTgl()) + static_cast<double>(trTPC.getSigma1PtTgl());
  covMat(3, 3) = static_cast<double>(trITS.getSigma1Pt2()) + static_cast<double>(trTPC.getSigma1Pt2());
  if (!covMat.Invert()) {
    LOG(error) << "Cov.matrix inversion failed: " << covMat;
    return 2. * o2::track::HugeF;
  }
  double chi2diag = 0., chi2ndiag = 0.,
         diff[o2::track::kNParams - 1] = {trITS.getParam(o2::track::kY) - trTPC.getParam(o2::track::kY),
                                          trITS.getParam(o2::track::kSnp) - trTPC.getParam(o2::track::kSnp),
                                          trITS.getParam(o2::track::kTgl) - trTPC.getParam(o2::track::kTgl),
                                          trITS.getParam(o2::track::kQ2Pt) - trTPC.getParam(o2::track::kQ2Pt)};
  for (int i = o2::track::kNParams - 1; i--;) {
    chi2diag += diff[i] * diff[i] * covMat(i, i);
    for (int j = i; j--;) {
      chi2ndiag += diff[i] * diff[j] * covMat(i, j);
    }
  }
  return chi2diag + 2. * chi2ndiag;
}

//______________________________________________
void MatchTPCITS::addLastTrackCloneForNeighbourSector(int sector, o2::track::TrackLTIntegral* trackLTInt)
{
  // add clone of the src ITS track cache, propagate it to ref.X in requested sector
  // and register its index in the sector cache. Used for ITS tracks which are so close
  // to their setctor edge that their matching should be checked also in the neighbouring sector
  mITSWork.push_back(mITSWork.back()); // clone the last track defined in given sector
  auto& trc = mITSWork.back();
  if (trc.rotate(o2::math_utils::sector2Angle(sector)) &&
      o2::base::Propagator::Instance()->PropagateToXBxByBz(trc, mParams->XMatchingRef, MaxSnp, 2., mUseMatCorrFlag, trackLTInt)) {
    // TODO: use faster prop here, no 3d field, materials
    mITSSectIndexCache[sector].push_back(mITSWork.size() - 1); // register track CLONE
    // flag clone
    mITSWork.back().setCloneBefore();
    if (trackLTInt) {
      mITSWork.back().xrho = trackLTInt->getXRho(); // we collect seen x*rho and distance to the reference X for further PID correcrions
      mITSWork.back().dL = trackLTInt->getL();
    }
    mITSWork[mITSWork.size() - 2].setCloneAfter();
    if (mMCTruthON) {
      mITSLblWork.emplace_back(mITSTrkLabels[trc.sourceID]);
    }
  } else {
    mITSWork.pop_back(); // rotation / propagation failed
  }
}

//______________________________________________
bool MatchTPCITS::propagateToRefX(o2::track::TrackParCov& trc, o2::track::TrackLTIntegral* lti)
{
  // propagate track to matching reference X, making sure its assigned alpha
  // is consistent with TPC sector
  constexpr float TgHalfSector = 0.17632698f;
  bool refReached = false;
  refReached = mParams->XMatchingRef < 10.; // RS: tmp, to cover XMatchingRef~0
  int trialsLeft = 2;
  while (o2::base::Propagator::Instance()->PropagateToXBxByBz(trc, mParams->XMatchingRef, MaxSnp, 2., mUseMatCorrFlag, lti)) {
    if (refReached) {
      break;
    }
    // make sure the track is indeed within the sector defined by alpha
    if (fabs(trc.getY()) < mParams->XMatchingRef * TgHalfSector) {
      refReached = true;
      break; // ok, within
    }
    if (!trialsLeft--) {
      break;
    }
    auto alphaNew = o2::math_utils::angle2Alpha(trc.getPhiPos());
    if (!trc.rotate(alphaNew) != 0) {
      break; // failed (RS: check effect on matching tracks to neighbouring sector)
    }
  }
  return refReached && std::abs(trc.getSnp()) < MaxSnp;
}

//______________________________________________
void MatchTPCITS::print() const
{
  ///< print all settings
  printf("\n******************** TPC-ITS matching component ********************\n");
  if (!mInitDone) {
    printf("init is not done yet\n");
    return;
  }

  printf("MC truth: %s\n", mMCTruthON ? "on" : "off");
  printf("Matching reference X: %.3f\n", mParams->XMatchingRef);
  printf("Account Z dimension: %s\n", mCompareTracksDZ ? "on" : "off");
  printf("Cut on matching chi2: %.3f\n", mParams->cutMatchingChi2);
  printf("Max number ITS candidates per TPC track: %d\n", mParams->maxMatchCandidates);
  printf("Crude cut on track params: ");
  for (int i = 0; i < o2::track::kNParams; i++) {
    printf(" %.3e", mParams->crudeAbsDiffCut[i]);
  }
  printf("\n");

  printf("NSigma^2 cut on track params: ");
  for (int i = 0; i < o2::track::kNParams; i++) {
    printf(" %6.2f", mParams->crudeNSigma2Cut[i]);
  }
  printf("\n");

  printf("TPC Z->time(bins) bracketing safety margin: %6.2f\n", mParams->safeMarginTPCTimeEdge);

#ifdef _ALLOW_DEBUG_TREES_

  printf("Output debug tree (%s) file: %s\n", mDBGFlags ? "on" : "off", mDebugTreeFileName.data());
  if (getDebugFlags()) {
    printf("Debug stream flags:\n");
    if (isDebugFlag(MatchTreeAll | MatchTreeAccOnly)) {
      printf("* matching canditate pairs: %s\n", isDebugFlag(MatchTreeAccOnly) ? "accepted" : "all");
    }
    if (isDebugFlag(WinnerMatchesTree)) {
      printf("* winner matches\n");
    }
  }
#endif

  printf("**********************************************************************\n");
}

//______________________________________________
void MatchTPCITS::refitWinners(pmr::vector<o2::dataformats::TrackTPCITS>& matchedTracks, pmr::vector<o2::MCCompLabel>& matchLabels, pmr::vector<o2::dataformats::Triplet<float, float, float>>& calib)
{
  ///< refit winning tracks
  mTimer[SWRefit].Start(false);
  matchedTracks.reserve(mNMatches + mABWinnersIDs.size());
  matchedTracks.resize(mNMatches);
  if (mMCTruthON) {
    matchLabels.reserve(mNMatches + mABWinnersIDs.size());
    matchLabels.resize(mNMatches);
  }
  if (mVDriftCalibOn) {
    calib.reserve(mNCalibPrelim * 1.2 + 1);
  }
  std::vector<int> tpcToFit;
  tpcToFit.reserve(mNMatches);
  for (int iTPC = 0; iTPC < (int)mTPCWork.size(); iTPC++) {
    const auto& tTPC = mTPCWork[iTPC];
    if (!isDisabledTPC(tTPC) && tTPC.gid.testBit(0)) {
      tpcToFit.push_back(iTPC);
    }
  }
  LOG(debug) << "Refitting winner matches";
  mWinnerChi2Refit.resize(mITSWork.size(), -1.f);
  int nToFit = (int)tpcToFit.size();
  unsigned int nFailedRefit{0};

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(mNThreads) \
  reduction(+                                                     \
            : nFailedRefit)
#endif
  for (int ifit = 0; ifit < nToFit; ifit++) {
    int iTPC = tpcToFit[ifit], iITS;
    const auto& tTPC = mTPCWork[iTPC];
    if (refitTrackTPCITS(ifit, iTPC, iITS, matchedTracks, matchLabels, calib)) {
      mWinnerChi2Refit[iITS] = matchedTracks.back().getChi2Refit();
    } else {
      ++nFailedRefit;
    }
  }
  LOGP(info, "Failed {} TPC-ITS refits out of {}", nFailedRefit, nToFit);

  // suppress tracks failed on refit and fill calib/debug data (if needed)
  int last = nToFit;
  mNMatches = 0;
  for (int ifit = 0; ifit < nToFit; ifit++) {
    int itpc = tpcToFit[ifit];
    if (!matchedTracks[ifit].isValid()) { // move the last good track from the back to the slot to delete
      while (--last > ifit && !matchedTracks[last].isValid()) {
      } // find the highest valid track
      if (last > ifit) {
        matchedTracks[ifit] = matchedTracks[last];
        matchedTracks[last].invalidate();
        itpc = tpcToFit[last];
        if (mMCTruthON) {
          matchLabels[ifit] = matchLabels[last];
        }
      } else {
        break;
      }
    }
    if (mDBGOut || mVDriftCalibOn) {
      fillCalibDebug(ifit, itpc, matchedTracks[ifit], calib);
    }
    mNMatches++;
  }
  // adjust sizes
  matchedTracks.resize(mNMatches);
  if (mMCTruthON) {
    matchLabels.resize(mNMatches);
  }
  mTimer[SWRefit].Stop();
}

//______________________________________________
void MatchTPCITS::fillCalibDebug(int ifit, int iTPC, const o2::dataformats::TrackTPCITS& match, pmr::vector<o2::dataformats::Triplet<float, float, float>>& calib)
{
  const auto& tTPC = mTPCWork[iTPC];
  int iITS = mMatchRecordsTPC[tTPC.matchID].partnerID;
  const auto& tITS = mITSWork[iITS];
  float minDiffFT0 = -999., timeC = 0.f;
  std::vector<float> dtimes;
  bool fillVDCalib = mVDriftCalibOn && (!mFieldON || std::abs(match.getQ2Pt()) < mParams->maxVDriftTrackQ2Pt);
  if (fillVDCalib || mDBGOut) {
    timeC = match.getTimeMUS().getTimeStamp(); // find closest FIT record
    float minDiffA = mParams->maxVDritTimeOffset;
    if (mInteractions.size()) {
      int timeC0 = timeC - minDiffA;
      if (timeC0 < 0) {
        timeC0 = 0;
      }
      auto entStart = timeC0 < int(mInteractionMUSLUT.size()) ? mInteractionMUSLUT[timeC0] : (mInteractionMUSLUT.size() ? mInteractionMUSLUT.back() : 0);
      for (int ent = entStart; ent < (int)mInteractions.size(); ent++) {
        float diff = mInteractions[ent].tBracket.mean() - timeC;
        if (diff > minDiffA) {
          break; // all following will be the same
        }
        if (diff < -minDiffA) {
          continue;
        }
        dtimes.push_back(diff);
        minDiffFT0 = diff;
        minDiffA = std::abs(minDiffFT0);
      }
    }
  }
  if (fillVDCalib) {
    calib.emplace_back(tITS.getTgl(), tTPC.getTgl(), minDiffFT0);
  }
#ifdef _ALLOW_DEBUG_TREES_
  if (mDBGOut) {
    o2::track::TrackPar itsRefPIDCorr(tITS);   // version ad hoc corrected for TPC PID being different from the pion, as it is done in the matching
    o2::track::TrackParCov itsRefAltPID(tITS); // version with full propagation to account for TPC PID being different from the pion
    itsRefPIDCorr.setX(0);
    if (tTPC.getPID() > tITS.getPID() && tITS.dL > 0.f && tTPC.getP2() / tTPC.getPID().getMass2() < mParams->minBetaGammaForPIDDiff) {
      itsRefAltPID = mITSTracksArray[tITS.sourceID].getParamOut(); // clone original ITS track at highest update point
      itsRefPIDCorr.setPID(tTPC.getPID(), true);
      itsRefPIDCorr = itsRefAltPID;
      // fast correction
      if (!itsRefPIDCorr.correctForELoss(tITS.xrho)) {
        itsRefPIDCorr.setX(-10);
      } else {
        float q2ptPID = itsRefPIDCorr.getQ2Pt();
        float dCurv = (q2ptPID - tITS.getQ2Pt()) * mBz * o2::constants::math::B2C, dLEff = tITS.dL * mParams->ITSStepEffFraction, dCurvL = dCurv * dLEff;
        itsRefPIDCorr = tITS;
        itsRefPIDCorr.setPID(tTPC.getPID(), true);
        itsRefPIDCorr.setQ2Pt(q2ptPID);
        auto snp = tITS.getSnp() + dCurvL;
        if (std::abs(snp) >= 1.) {
          snp = std::copysign(0.99, snp);
        }
        itsRefPIDCorr.setSnp(snp);
        itsRefPIDCorr.setY(tITS.getY() + dCurvL * dLEff * 0.5);
      }
      // full propagation
      if (!itsRefAltPID.rotate(tTPC.getAlpha()) || !o2::base::Propagator::Instance()->PropagateToXBxByBz(itsRefAltPID, mParams->XMatchingRef, MaxSnp, 2., mUseMatCorrFlag)) {
        itsRefAltPID.setX(-10);
      }
    }
    int tb = mTPCTracksArray[tTPC.sourceID].getTime0() * mNTPCOccBinLengthInv;
    float mltTPC = tb < 0 ? mTBinClOcc[0] : (tb >= mTBinClOcc.size() ? mTBinClOcc.back() : mTBinClOcc[tb]);
    (*mDBGOut) << "refit"
               << "tpcOrig=" << mTPCTracksArray[tTPC.sourceID] << "itsOrig=" << mITSTracksArray[tITS.sourceID] << "itsRef=" << tITS << "tpcRef=" << tTPC << "matchRefit=" << match
               << "timeCorr=" << timeC << "dTimeFT0=" << minDiffFT0 << "dTimes=" << dtimes
               << "itsRefAltPID=" << itsRefAltPID << "itsRefPIDCorr=" << itsRefPIDCorr;
    if (mMCTruthON) {
      (*mDBGOut) << "refit"
                 << "itsLbl=" << mITSLblWork[iITS] << "tpcLbl=" << mTPCLblWork[iTPC];
    }
    (*mDBGOut) << "refit"
               << "multTPC=" << mltTPC
               << "multITSTr=" << mITSTrackROFRec[tITS.roFrame]
               << "multITSCl=" << mITSClusterROFRec[tITS.roFrame]
               << "tf=" << mTFCount << "\n";
  }
#endif
}

//______________________________________________
bool MatchTPCITS::refitTrackTPCITS(int slot, int iTPC, int& iITS, pmr::vector<o2::dataformats::TrackTPCITS>& matchedTracks, pmr::vector<o2::MCCompLabel>& matchLabels, pmr::vector<o2::dataformats::Triplet<float, float, float>>& calib)
{
  ///< refit in inward direction the pair of TPC and ITS tracks

  const float maxStep = 2.f; // max propagation step (TODO: tune)
  const auto& tTPC = mTPCWork[iTPC];
  const auto& tpcMatchRec = mMatchRecordsTPC[tTPC.matchID];
  iITS = tpcMatchRec.partnerID;
  const auto& tITS = mITSWork[iITS];
  const auto& itsTrOrig = mITSTracksArray[tITS.sourceID];
  auto& trfit = matchedTracks[slot];
  ((o2::track::TrackParCov&)trfit) = (const o2::track::TrackParCov&)tTPC;
  trfit.getParamOut() = (const o2::track::TrackParCov&)tITS; // create a copy of TPC track at xRef
  trfit.getParamOut().setUserField(0);                       // reset eventual clones flag
  trfit.setPID(tTPC.getPID(), true);
  trfit.getParamOut().setPID(tTPC.getPID(), true);
  // in continuos mode the Z of TPC track is meaningless, unless it is CE crossing
  // track (currently absent, TODO)
  if (!mCompareTracksDZ) {
    trfit.setZ(tITS.getZ()); // fix the seed Z
  }
  float deltaT = (trfit.getZ() - tTPC.getZ()) * mTPCVDriftInv;                                                                                   // time correction in \mus
  float timeErr = tTPC.constraint == TrackLocTPC::Constrained ? tTPC.timeErr : std::sqrt(tITS.getSigmaZ2() + tTPC.getSigmaZ2()) * mTPCVDriftInv; // estimate the error on time
  if (timeErr > mITSTimeResMUS && tTPC.constraint != TrackLocTPC::Constrained) {
    timeErr = mITSTimeResMUS; // chose smallest error
    deltaT = tTPC.constraint == TrackLocTPC::ASide ? tITS.tBracket.mean() - tTPC.time0 : tTPC.time0 - tITS.tBracket.mean();
  }
  timeErr += mParams->globalTimeExtraErrorMUS;
  float timeC = tTPC.getCorrectedTime(deltaT) + mParams->globalTimeBiasMUS; /// precise time estimate, optionally corrected for bias
  o2::math_utils::Bracketf_t::Relation relITS;
  if (mParams->ITSTimeOutliersPolicy != MatchTPCITSParams::TimeOutliersPolicy::Tolerate && (relITS = tITS.tBracket.isOutside(timeC))) { /// track time is outside of the nominal ITS time
    if (mParams->ITSTimeOutliersPolicy == MatchTPCITSParams::TimeOutliersPolicy::Adjust) {
      timeC = relITS == o2::math_utils::Bracketf_t::Below ? tITS.tBracket.getMin() : tITS.tBracket.getMax(); /// assign ITS boundary
    } else {                                                                                                 // == MatchTPCITSParams::TimeOutliersPolicy::Reject
      trfit.invalidate();
      return false;
    }
  }
  if (timeC < 0) { // RS TODO similar check is needed for other edge of TF
    if (timeC + std::min(timeErr, mParams->tfEdgeTimeToleranceMUS * mTPCTBinMUSInv) < 0) {
      trfit.invalidate(); // destroy failed track
      return false;
    }
    timeC = 0.;
  }

  // refit TPC track inward into the ITS
  int nclRefit = 0, ncl = itsTrOrig.getNumberOfClusters();
  float chi2 = 0.f;
  auto geom = o2::its::GeometryTGeo::Instance();
  auto propagator = o2::base::Propagator::Instance();
  int clEntry = itsTrOrig.getFirstClusterEntry();

  float addErr2 = 0;
  // extra error on tgl due to the assumed vdrift uncertainty
  if (mVDriftCalibOn) {
    addErr2 = tITS.getParam(o2::track::kTgl) * mParams->maxVDriftUncertainty;
    addErr2 *= addErr2;
    trfit.updateCov(addErr2, o2::track::kSigTgl2);
  }

  for (int icl = 0; icl < ncl; icl++) {
    const auto& clus = mITSClustersArray[mITSTrackClusIdx[clEntry++]];
    float alpha = geom->getSensorRefAlpha(clus.getSensorID()), x = clus.getX();
    if (!trfit.rotate(alpha) ||
        // note: here we also calculate the L,T integral (in the inward direction, but this is irrelevant)
        // note: we should eventually use TPC pid in the refit (TODO)
        // note: since we are at small R, we can use field BZ component at origin rather than 3D field
        !propagator->propagateToX(trfit, x, propagator->getNominalBz(),
                                  MaxSnp, maxStep, mUseMatCorrFlag, &trfit.getLTIntegralOut())) {
      break;
    }
    chi2 += trfit.getPredictedChi2(clus);
    if (!trfit.update(clus)) {
      break;
    }
    nclRefit++;
  }
  if (nclRefit != ncl) {
    LOGP(debug, "Refit in ITS failed after ncl={}, match between TPC track #{} and ITS track #{}", nclRefit, tTPC.sourceID, tITS.sourceID);
    LOGP(debug, "{:s}", trfit.asString());
    trfit.invalidate(); // destroy failed track
    return false;
  }

  // We need to update the LTOF integral by the distance to the "primary vertex"
  // We want to leave the track at the the position of its last update, so we do a fast propagation on the TrackPar copy of trfit,
  // and since for the LTOF calculation the material effects are irrelevant, we skip material corrections
  const o2::dataformats::VertexBase vtxDummy; // at the moment using dummy vertex: TODO use MeanVertex constraint instead
  o2::track::TrackPar trpar(trfit);
  if (!propagator->propagateToDCA(vtxDummy.getXYZ(), trpar, propagator->getNominalBz(),
                                  maxStep, MatCorrType::USEMatCorrNONE, nullptr, &trfit.getLTIntegralOut())) {
    LOG(error) << "LTOF integral might be incorrect";
  }
  auto& tofL = trfit.getLTIntegralOut(); // this is TL integral calculated from the RefX to the DCA to the beamline, invert material integrals for outward propagation
  tofL.setX2X0(-tofL.getX2X0());
  tofL.setXRho(-tofL.getXRho());

  // outward refit
  auto& tracOut = trfit.getParamOut(); // this is a clone of ITS outward track already at the matching reference X
  if (tTPC.getPID() > tITS.getPID() && tTPC.getP2() / tTPC.getPID().getMass2() < mParams->minBetaGammaForPIDDiff) {
    // in case the TPC track hypothesis is not pion, we redo the outward propagation to ref.x with TPC PID
    tracOut = mITSTracksArray[tITS.sourceID].getParamOut();
    tracOut.setPID(tTPC.getPID(), true);
    if (!tracOut.rotate(tTPC.getAlpha()) || !o2::base::Propagator::Instance()->PropagateToXBxByBz(tracOut, mParams->XMatchingRef, MaxSnp, 2., mUseMatCorrFlag)) {
      LOGP(debug, "Failed to rotate ITSouter with imposed PID to TPC alpha {} or propagate to X={}: {:s}", tTPC.getAlpha(), mParams->XMatchingRef, tracOut.asString());
      trfit.invalidate(); // destroy failed track
      return false;
    }
  }
  {
    float xtogo = 0;
    if (!tracOut.getXatLabR(o2::constants::geom::XTPCInnerRef, xtogo, mBz, o2::track::DirOutward) ||
        !propagator->PropagateToXBxByBz(tracOut, xtogo, MaxSnp, 10., mUseMatCorrFlag, &tofL)) {
      LOG(debug) << "Propagation to inner TPC boundary X=" << xtogo << " failed, Xtr=" << tracOut.getX() << " snp=" << tracOut.getSnp();
      trfit.invalidate(); // destroy failed track
      return false;
    }
    if (mVDriftCalibOn) {
      tracOut.updateCov(addErr2, o2::track::kSigTgl2);
    }
    float chi2Out = 0;
    auto posStart = tracOut.getXYZGlo();
    auto tImposed = timeC * mTPCTBinMUSInv;
    if (std::abs(tImposed - mTPCTracksArray[tTPC.sourceID].getTime0()) > 550) {
      LOGP(alarm, "Impossible imposed timebin {} for TPC track time0:{}, dBwd:{} dFwd:{} TB | ZShift:{}, TShift:{}", tImposed, mTPCTracksArray[tTPC.sourceID].getTime0(),
           mTPCTracksArray[tTPC.sourceID].getDeltaTBwd(), mTPCTracksArray[tTPC.sourceID].getDeltaTFwd(), trfit.getZ() - tTPC.getZ(), deltaT);
      LOGP(info, "Trc: {}", mTPCTracksArray[tTPC.sourceID].asString());
      trfit.invalidate(); // destroy failed track
      return false;
    }
    int retVal = mTPCRefitter->RefitTrackAsTrackParCov(tracOut, mTPCTracksArray[tTPC.sourceID].getClusterRef(), tImposed, &chi2Out, true, false); // outward refit
    if (retVal < 0) {
      LOG(debug) << "Refit failed";
      trfit.invalidate(); // destroy failed track
      return false;
    }
    auto posEnd = tracOut.getXYZGlo();
    auto lInt = propagator->estimateLTIncrement(tracOut, posStart, posEnd);
    tofL.addStep(lInt, tracOut.getP2Inv());
    tofL.addX2X0(lInt * mTPCmeanX0Inv);
    propagator->PropagateToXBxByBz(tracOut, o2::constants::geom::XTPCOuterRef, MaxSnp, 10., mUseMatCorrFlag, &tofL);

    const auto& trackTune = TrackTuneParams::Instance();
    if (trackTune.useTPCOuterCorr) {
      tracOut.updateParams(trackTune.tpcParOuter);
    }
    if (trackTune.tpcCovOuterType != TrackTuneParams::AddCovType::Disable) {
      tracOut.updateCov(mCovDiagOuter, trackTune.tpcCovOuterType == TrackTuneParams::AddCovType::WithCorrelations);
    }
  }
  trfit.setChi2Match(tpcMatchRec.chi2);
  trfit.setChi2Refit(chi2);
  trfit.setTimeMUS(timeC, timeErr);
  trfit.setRefTPC({unsigned(tTPC.sourceID), o2::dataformats::GlobalTrackID::TPC});
  trfit.setRefITS({unsigned(tITS.sourceID), o2::dataformats::GlobalTrackID::ITS});

  if (mMCTruthON) { // store MC info: we assign TPC track label and declare the match fake if the ITS and TPC labels are different (their fake flag is ignored)
    matchLabels[slot] = mTPCLblWork[iTPC];
    matchLabels[slot].setFakeFlag(mITSLblWork[iITS] != mTPCLblWork[iTPC]);
  }

  return true;
}

//______________________________________________
bool MatchTPCITS::refitABTrack(int iITSAB, const TPCABSeed& seed, pmr::vector<o2::dataformats::TrackTPCITS>& matchedTracks, pmr::vector<int>& ABTrackletClusterIDs, pmr::vector<o2::itsmft::TrkClusRef>& ABTrackletRefs)
{
  ///< refit AfterBurner track

  const float maxStep = 2.f; // max propagation step (TODO: tune)
  const auto& tTPC = mTPCWork[seed.tpcWID];
  const auto& winLink = seed.getLink(seed.winLinkID);
  auto& newtr = matchedTracks.emplace_back(winLink, winLink); // create a copy of winner param at innermost ITS cluster
  auto& tracOut = newtr.getParamOut();
  auto& tofL = newtr.getLTIntegralOut();
  auto geom = o2::its::GeometryTGeo::Instance();
  auto propagator = o2::base::Propagator::Instance();
  tracOut.resetCovariance();
  propagator->estimateLTFast(tofL, winLink); // guess about initial value for the track integral from the origin
  // refit track outward in the ITS
  const auto& itsClRefs = ABTrackletRefs[iITSAB];
  int nclRefit = 0, ncl = itsClRefs.getNClusters();

  float chi2 = 0.f;
  // NOTE: the ITS cluster absolute indices are stored from inner to outer layers
  for (int icl = itsClRefs.getFirstEntry(); icl < itsClRefs.getEntriesBound(); icl++) {
    const auto& clus = mITSClustersArray[ABTrackletClusterIDs[icl]];
    float alpha = geom->getSensorRefAlpha(clus.getSensorID()), x = clus.getX();
    if (!tracOut.rotate(alpha) ||
        // note: here we also calculate the L,T integral
        // note: we should eventually use TPC pid in the refit (TODO)
        // note: since we are at small R, we can use field BZ component at origin rather than 3D field
        !propagator->propagateToX(tracOut, x, propagator->getNominalBz(), MaxSnp, maxStep, mUseMatCorrFlag, &tofL)) {
      break;
    }
    chi2 += tracOut.getPredictedChi2(clus);
    if (!tracOut.update(clus)) {
      break;
    }
    nclRefit++;
  }
  if (nclRefit != ncl) {
    LOGP(debug, "AfterBurner refit in ITS failed after ncl={}, match between TPC track #{} and ITS tracklet #{}", nclRefit, tTPC.sourceID, iITSAB);
    LOGP(debug, "{:s}", tracOut.asString());
    matchedTracks.pop_back(); // destroy failed track
    return false;
  }
  // perform TPC refit with interaction time constraint
  float timeC = mInteractions[seed.ICCanID].tBracket.mean();
  float timeErr = mInteractions[seed.ICCanID].tBracket.delta(); // RS FIXME shall we use gaussian error?
  {
    float xtogo = 0;
    if (!tracOut.getXatLabR(o2::constants::geom::XTPCInnerRef, xtogo, mBz, o2::track::DirOutward) ||
        !propagator->PropagateToXBxByBz(tracOut, xtogo, MaxSnp, 10., mUseMatCorrFlag, &tofL)) {
      LOG(debug) << "Propagation to inner TPC boundary X=" << xtogo << " failed, Xtr=" << tracOut.getX() << " snp=" << tracOut.getSnp();
      matchedTracks.pop_back(); // destroy failed track
      return false;
    }
    float chi2Out = 0;
    auto posStart = tracOut.getXYZGlo();
    int retVal = mTPCRefitter->RefitTrackAsTrackParCov(tracOut, mTPCTracksArray[tTPC.sourceID].getClusterRef(), timeC * mTPCTBinMUSInv, &chi2Out, true, false); // outward refit
    if (retVal < 0) {
      LOG(debug) << "Refit failed";
      matchedTracks.pop_back(); // destroy failed track
      return false;
    }
    auto posEnd = tracOut.getXYZGlo();
    auto lInt = propagator->estimateLTIncrement(tracOut, posStart, posEnd);
    tofL.addStep(lInt, tracOut.getP2Inv());
    tofL.addX2X0(lInt * mTPCmeanX0Inv);
    propagator->PropagateToXBxByBz(tracOut, o2::constants::geom::XTPCOuterRef, MaxSnp, 10., mUseMatCorrFlag, &tofL);
    const auto& trackTune = TrackTuneParams::Instance();
    if (trackTune.useTPCOuterCorr) {
      tracOut.updateParams(trackTune.tpcParOuter);
    }
    if (trackTune.tpcCovOuterType != TrackTuneParams::AddCovType::Disable) {
      tracOut.updateCov(mCovDiagOuter, trackTune.tpcCovOuterType == TrackTuneParams::AddCovType::WithCorrelations);
    }
  }

  newtr.setChi2Match(winLink.chi2Norm());
  newtr.setChi2Refit(chi2);
  newtr.setTimeMUS(timeC, timeErr);
  newtr.setRefTPC({unsigned(tTPC.sourceID), o2::dataformats::GlobalTrackID::TPC});
  newtr.setRefITS({unsigned(iITSAB), o2::dataformats::GlobalTrackID::ITSAB});

  return true;
}

//______________________________________________
bool MatchTPCITS::refitTPCInward(o2::track::TrackParCov& trcIn, float& chi2, float xTgt, int trcID, float timeTB) const
{
  // inward refit
  const auto& tpcTrOrig = mTPCTracksArray[trcID];

  trcIn = tpcTrOrig.getOuterParam();
  chi2 = 0;

  auto propagator = o2::base::Propagator::Instance();
  int retVal = mTPCRefitter->RefitTrackAsTrackParCov(trcIn, tpcTrOrig.getClusterRef(), timeTB, &chi2, false, true); // inward refit with matrix reset
  if (retVal < 0) {
    LOG(warning) << "Refit failed";
    LOG(warning) << trcIn.asString();
    return false;
  }
  //
  // propagate to the inner edge of the TPC
  // Note: it is allowed to not reach the requested radius
  if (!propagator->PropagateToXBxByBz(trcIn, xTgt, MaxSnp, 2., mUseMatCorrFlag)) {
    LOG(debug) << "Propagation to target X=" << xTgt << " failed, Xtr=" << trcIn.getX() << " snp=" << trcIn.getSnp() << " pT=" << trcIn.getPt();
    LOG(debug) << trcIn.asString();
    return false;
  }
  return true;
}

//>>============================= AfterBurner for TPC-track / ITS cluster matching ===================>>
//______________________________________________
int MatchTPCITS::prepareABSeeds()
{
  ///< select TPC tracks to be considered in afterburner, clone them as seeds for every matching interaction candidate
  const auto& outerLr = mRGHelper.layers.back();
  // to avoid difference between 3D field propagation and Bz-bazed getXatLabR we propagate RMax+margin
  const float ROuter = outerLr.rRange.getMax() + 0.5f;
  auto propagator = o2::base::Propagator::Instance();

  for (int iTPC = 0; iTPC < (int)mTPCWork.size(); iTPC++) {
    auto& tTPC = mTPCWork[iTPC];
    if (isDisabledTPC(tTPC)) {
      // Popagate to the vicinity of the out layer. Note: the Z of the track might be uncertain,
      // in this case the material corrections will be correct only in the limit of their uniformity in Z,
      // which should be good assumption....
      float xTgt;
      if (!tTPC.getXatLabR(ROuter, xTgt, propagator->getNominalBz(), o2::track::DirInward) ||
          !propagator->PropagateToXBxByBz(tTPC, xTgt, MaxSnp, 2., mUseMatCorrFlag)) {
        continue;
      }
      mTPCABIndexCache.push_back(iTPC);
    }
  }
  // sort tracks according to their timeMin
  std::sort(mTPCABIndexCache.begin(), mTPCABIndexCache.end(), [this](int a, int b) {
    auto& trcA = mTPCWork[a];
    auto& trcB = mTPCWork[b];
    return (trcA.tBracket.getMin() - trcB.tBracket.getMin()) < 0.;
  });

  float maxTDriftSafe = tpcTimeBin2MUS(mNTPCBinsFullDrift + mParams->safeMarginTPCITSTimeBin + mTPCTimeEdgeTSafeMargin);
  int nIntCand = mInteractions.size(), nTPCCand = mTPCABIndexCache.size();
  int tpcStart = 0;
  for (int ic = 0; ic < nIntCand; ic++) {
    int icFirstSeed = mTPCABSeeds.size();
    auto& intCand = mInteractions[ic];
    auto tic = intCand.tBracket.mean();
    for (int it = tpcStart; it < nTPCCand; it++) {
      auto& trc = mTPCWork[mTPCABIndexCache[it]];
      auto cmp = trc.tBracket.isOutside(intCand.tBracket);
      if (cmp < 0) {
        break; // all other TPC tracks will be also in future wrt the interaction
      }
      if (cmp > 0) {
        if (trc.tBracket.getMin() + maxTDriftSafe < intCand.tBracket.getMin()) {
          tpcStart++; // all following int.candidates would be in future wrt this track
        }
        continue;
      }
      // we beed to create seed from this TPC track and interaction candidate
      float dt = trc.getSignedDT(tic - trc.time0);
      float dz = dt * mTPCVDrift, z = trc.getZ() + dz;
      if (outerLr.zRange.isOutside(z, std::sqrt(trc.getSigmaZ2()) + 2.)) { // RS FIXME introduce margin as parameter?
        continue;
      }
      // make sure seed crosses the outer ITS layer (with some margin)
      auto& seed = mTPCABSeeds.emplace_back(mTPCABIndexCache[it], ic, trc);
      seed.track.setZ(z); // RS FIXME : in case of distortions and large dz the track must be refitted
    }
    intCand.seedsRef.set(icFirstSeed, mTPCABSeeds.size() - icFirstSeed);
  }
  return mTPCABIndexCache.size();
}

//______________________________________________
int MatchTPCITS::prepareInteractionTimes()
{
  // guess interaction times from various sources and relate with ITS rofs
  const float ft0Uncertainty = 0.5e-3;
  int nITSROFs = mITSROFTimes.size();
  if (mFITInfo.size()) {
    int rof = 0;
    for (const auto& ft : mFITInfo) {
      if (!mFT0Params->isSelected(ft)) {
        continue;
      }
      auto fitTime = ft.getInteractionRecord().differenceInBCMUS(mStartIR);
      if (fitTime < 0) { // should not happen
        continue;
      }
      if (size_t(fitTime) >= mInteractionMUSLUT.size()) {
        mInteractionMUSLUT.resize(size_t(fitTime) + 1, -1);
      }
      if (mInteractionMUSLUT[fitTime] < 0) {
        mInteractionMUSLUT[fitTime] = mInteractions.size();
      }
      for (; rof < nITSROFs; rof++) {
        if (mITSROFTimes[rof] < fitTime) {
          continue;
        }
        break;
      }
      if (rof >= nITSROFs) {
        break;
      }
      mInteractions.emplace_back(ft.getInteractionRecord(), fitTime, ft0Uncertainty, rof, o2::detectors::DetID::FT0);
    }
  }
  int ent = 0;
  for (auto& val : mInteractionMUSLUT) {
    if (val < 0) { // was not assigned == no interactions in this mus, assign previous one
      val = ent;
    } else {
      ent = val > 0 ? val - 1 : val;
    }
  }
  return mInteractions.size();
}

//______________________________________________
bool MatchTPCITS::runAfterBurner(pmr::vector<o2::dataformats::TrackTPCITS>& matchedTracks, pmr::vector<o2::MCCompLabel>& matchLabels, pmr::vector<o2::MCCompLabel>& ABTrackletLabels,
                                 pmr::vector<int>& ABTrackletClusterIDs, pmr::vector<o2::itsmft::TrkClusRef>& ABTrackletRefs, pmr::vector<o2::dataformats::Triplet<float, float, float>>& calib)
{
  mTimer[SWABSeeds].Start(false);
  mNABRefsClus = 0;
  prepareABSeeds();
  int nIntCand = mInteractions.size(), nABSeeds = mTPCABSeeds.size();
  LOGP(info, "AfterBurner will check {} seeds from {} TPC tracks and {} interaction candidates with {} threads", nABSeeds, mTPCABIndexCache.size(), nIntCand, mNThreads); // TMP
  mTimer[SWABSeeds].Stop();
  if (!nIntCand || !mTPCABSeeds.size()) {
    return false;
  }
  mTimer[SWABMatch].Start(false);

  std::vector<ITSChipClustersRefs> itsChipClRefsBuff(mNThreads);
#ifdef ENABLE_UPGRADES
  // with upgrades the datatype changed, hence we need to initialize
  // each element individually
  std::generate(itsChipClRefsBuff.begin(), itsChipClRefsBuff.end(), []() {
    return ITSChipClustersRefs(o2::its::GeometryTGeo::Instance()->getNumberOfChips());
  });
#endif

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(mNThreads)
#endif
  for (int ic = 0; ic < nIntCand; ic++) {
    const auto& intCand = mInteractions[ic];
    LOGP(debug, "cand T {} Entries: {} : {} : {} | ITS ROF: {}", intCand.tBracket.mean(), intCand.seedsRef.getEntries(), intCand.seedsRef.getFirstEntry(), intCand.seedsRef.getEntriesBound(), intCand.rofITS);
    if (!intCand.seedsRef.getEntries()) {
      continue;
    }
#ifdef WITH_OPENMP
    uint8_t tid = (uint8_t)omp_get_thread_num();
#else
    uint8_t tid = 0;
#endif
    fillClustersForAfterBurner(intCand.rofITS, 1, itsChipClRefsBuff[tid]);                           // RS FIXME account for possibility of filling 2 ROFs
    for (int is = intCand.seedsRef.getFirstEntry(); is < intCand.seedsRef.getEntriesBound(); is++) { // loop over all seeds of this interaction candidate
      processABSeed(is, itsChipClRefsBuff[tid], tid);
    }
  }
  mTimer[SWABMatch].Stop();
  mTimer[SWABWinners].Start(false);
  int nwin = 0;
  // select winners
  struct SID {
    int seedID = -1;
    float chi2 = 1e9;
  };
  std::vector<SID> candAB;
  candAB.reserve(nABSeeds);
  mABWinnersIDs.reserve(mTPCABIndexCache.size());

  for (int i = 0; i < nABSeeds; i++) {
    auto& ABSeed = mTPCABSeeds[i];
    if (ABSeed.isDisabled()) {
      continue;
    }
    candAB.emplace_back(SID{i, ABSeed.getLink(ABSeed.getBestLinkID()).chi2Norm()});
  }
  std::sort(candAB.begin(), candAB.end(), [](SID a, SID b) { return a.chi2 < b.chi2; });
  for (int i = 0; i < (int)candAB.size(); i++) {
    auto& ABSeed = mTPCABSeeds[candAB[i].seedID];
    if (ABSeed.isDisabled()) {
      // RSTMP      LOG(info) << "Iter: " << iter << " seed is disabled: " << i << "[" << candAB[i].seedID << "/" << candAB[i].chi2 << "]"  << " last lr: " << int(ABSeed.lowestLayer);
      continue;
    }
    auto& tTPC = mTPCWork[ABSeed.tpcWID];
    if (tTPC.matchID > MinusOne) { // this tracks was already validated with other IC
      ABSeed.disable();
      // RSTMP      LOG(info) << "Iter: " << iter << " disabling seed " << i << "[" << candAB[i].seedID << "/" << candAB[i].chi2 << "]" << " TPC track " << ABSeed.tpcWID << " already validated"  << " last lr: " << int(ABSeed.lowestLayer);
      continue;
    }
    auto bestID = ABSeed.getBestLinkID();
    if (ABSeed.checkLinkHasUsedClusters(bestID, mABClusterLinkIndex)) {
      ABSeed.setNeedAlternative(); // flag for later processing
      // RSTMP      LOG(info) << "Iter: " << iter << " seed has used clusters " << i << "[" << candAB[i].seedID << "/" << candAB[i].chi2 << "]"  << " last lr: " << int(ABSeed.lowestLayer) << " Ncont: " << int(link.nContLayers);;
      continue;
    }
    ABSeed.validate(bestID);
    ABSeed.flagLinkUsedClusters(bestID, mABClusterLinkIndex);
    mABWinnersIDs.push_back(tTPC.matchID = candAB[i].seedID);
    mNABRefsClus += ABSeed.getNLayers();
    nwin++;
    // RSTMP      LOG(info) << "Iter: " << iter << " validated seed " << i << "[" << candAB[i].seedID << "/" << candAB[i].chi2 << "] for TPC track " << ABSeed.tpcWID << " last lr: " << int(ABSeed.lowestLayer) << " Ncont: " << int(link.nContLayers);
  }
  mTimer[SWABWinners].Stop();
  mTimer[SWABRefit].Start(false);
  refitABWinners(matchedTracks, matchLabels, ABTrackletLabels, ABTrackletClusterIDs, ABTrackletRefs, calib);
  mTimer[SWABRefit].Stop();
  return true;
}

//______________________________________________
void MatchTPCITS::refitABWinners(pmr::vector<o2::dataformats::TrackTPCITS>& matchedTracks, pmr::vector<o2::MCCompLabel>& matchLabels, pmr::vector<o2::MCCompLabel>& ABTrackletLabels,
                                 pmr::vector<int>& ABTrackletClusterIDs, pmr::vector<o2::itsmft::TrkClusRef>& ABTrackletRefs, pmr::vector<o2::dataformats::Triplet<float, float, float>>& calib)
{
  // refit normal matches
  refitWinners(matchedTracks, matchLabels, calib);

  ABTrackletClusterIDs.reserve(mNABRefsClus);
  ABTrackletRefs.reserve(mABWinnersIDs.size());
  if (mMCTruthON) {
    ABTrackletLabels.reserve(mABWinnersIDs.size());
  }
  if (matchedTracks.capacity() < mABWinnersIDs.size() + matchedTracks.size()) {
    LOGP(warn, "need to expand matched tracks container from {} to {}", matchedTracks.capacity(), mABWinnersIDs.size() + matchedTracks.size());
    matchedTracks.reserve(mABWinnersIDs.size() + matchedTracks.size());
    if (mMCTruthON) {
      matchLabels.reserve(mABWinnersIDs.size() + matchedTracks.size());
    }
  }

  std::map<o2::MCCompLabel, int> labelOccurence;
  auto accountClusterLabel = [&labelOccurence, itsClLabs = mITSClsLabels](int clID) {
    auto labels = itsClLabs->getLabels(clID);
    for (auto lab : labels) { // check all labels of the cluster
      if (lab.isSet()) {
        labelOccurence[lab]++;
      }
    }
  };

  for (auto wid : mABWinnersIDs) {
    const auto& ABSeed = mTPCABSeeds[wid];
    int start = ABTrackletClusterIDs.size();
    int lID = ABSeed.winLinkID, ncl = 0;
    auto& clref = ABTrackletRefs.emplace_back(start, ncl);
    while (lID > MinusOne) {
      const auto& winL = ABSeed.getLink(lID);
      if (winL.clID > MinusOne) {
        ABTrackletClusterIDs.push_back(winL.clID);
        ncl++;
        clref.pattern |= 0x1 << winL.layerID;
        clref.setClusterSize(winL.layerID, mITSClusterSizes[winL.clID]);
        if (mMCTruthON) {
          accountClusterLabel(winL.clID);
        }
      }
      lID = winL.parentID;
    }
    clref.setEntries(ncl);
    if (!refitABTrack(ABTrackletRefs.size() - 1, ABSeed, matchedTracks, ABTrackletClusterIDs, ABTrackletRefs)) { // on failure, destroy added tracklet reference
      ABTrackletRefs.pop_back();
      ABTrackletClusterIDs.resize(start); // RSS
      if (mMCTruthON) {
        labelOccurence.clear();
      }
      continue;
    }
    if (mMCTruthON) {
      o2::MCCompLabel lab;
      int maxL = 0; // find most encountered label
      for (auto [label, count] : labelOccurence) {
        if (count > maxL) {
          maxL = count;
          lab = label;
        }
      }
      if (maxL < ncl) {
        lab.setFakeFlag();
      }
      labelOccurence.clear();
      ABTrackletLabels.push_back(lab); // ITSAB tracklet label
      auto& lblGlo = matchLabels.emplace_back(mTPCLblWork[ABSeed.tpcWID]);
      lblGlo.setFakeFlag(lab != lblGlo);
      LOG(debug) << "ABWinner ncl=" << ncl << " mcLBAB " << lab << " mcLBGlo " << lblGlo << " chi2=" << ABSeed.getLink(ABSeed.winLinkID).chi2Norm() << " pT = " << ABSeed.track.getPt();
    }
    // build MC label
  }
  LOG(info) << "AfterBurner validated " << ABTrackletRefs.size() << " tracks";
}

//______________________________________________
void MatchTPCITS::processABSeed(int sid, const ITSChipClustersRefs& itsChipClRefs, uint8_t tID)
{
  // prepare matching hypothesis tree for given seed
  auto& ABSeed = mTPCABSeeds[sid];
  ABSeed.threadID = tID;
  ABSeed.linksEntry = mABLinksPool.threadPool[tID].size();
  followABSeed(ABSeed.track, itsChipClRefs, MinusTen, NITSLayers - 1, ABSeed); // check matches on outermost layer
  for (int ilr = NITSLayers - 1; ilr > mParams->lowestLayerAB; ilr--) {
    int nextLinkID = ABSeed.firstInLr[ilr];
    if (nextLinkID < 0) {
      break;
    }
    while (nextLinkID > MinusOne) {
      const auto& seedLink = ABSeed.getLink(nextLinkID);
      if (seedLink.isDisabled()) {
        nextLinkID = seedLink.nextOnLr;
        continue;
      }
      int next2nextLinkID = seedLink.nextOnLr;                            // fetch now since the seedLink may change due to the relocation
      followABSeed(seedLink, itsChipClRefs, nextLinkID, ilr - 1, ABSeed); // check matches on the next layer
      nextLinkID = next2nextLinkID;
    }
  }
  // is this seed has chance to be validated?
  auto candID = ABSeed.getBestLinkID();
  if (ABSeed.isDisabled() ||
      ABSeed.lowestLayer > mParams->requireToReachLayerAB ||
      candID < 0 ||
      ABSeed.getLink(candID).nContLayers < mParams->minContributingLayersAB) { // free unused links
    ABSeed.disable();
    mABLinksPool.threadPool[tID].resize(size_t(ABSeed.linksEntry));
  }

  /* // RS FIXME remove on final clean-up
  auto bestLinkID = ABSeed.getBestLinkID();
  if (bestLinkID>MinusOne) {
    const auto& bestL = ABSeed.getLink(bestLinkID);
    LOG(info) << "seed " << sid << " last lr: " << int(ABSeed.lowestLayer) << " Ncont: " << int(bestL.nContLayers) << " chi2 " << bestL.chi2;
  }
  else {
    LOG(info) << "seed " << sid << " : NONE";
  }
  */
}

//______________________________________________
int MatchTPCITS::followABSeed(const o2::track::TrackParCov& seed, const ITSChipClustersRefs& itsChipClRefs, int seedID, int lrID, TPCABSeed& ABSeed)
{

  auto propagator = o2::base::Propagator::Instance();
  float xTgt;
  const auto& lr = mRGHelper.layers[lrID];
  auto seedC = seed;
  if (!seedC.getXatLabR(lr.rRange.getMax(), xTgt, propagator->getNominalBz(), o2::track::DirInward) ||
      !propagator->propagateToX(seedC, xTgt, propagator->getNominalBz(), MaxSnp, 2., mUseMatCorrFlag)) { // Bz-propagation only in ITS
    return -1;
  }

  float zDRStep = -seedC.getTgl() * lr.rRange.delta(); // approximate Z span when going from layer rMin to rMax
  float errZ = std::sqrt(seedC.getSigmaZ2() + mParams->err2ABExtraZ);
  if (lr.zRange.isOutside(seedC.getZ(), mParams->nABSigmaZ * errZ + std::abs(zDRStep))) {
    return 0;
  }
  std::vector<int> chipSelClusters;        // preliminary cluster candidates //RS TODO do we keep this local / consider array instead of vector
  o2::math_utils::CircleXYf_t trcCircle;   // circle parameters for B ON data
  o2::math_utils::IntervalXYf_t trcLinPar; // line parameters for B OFF data
  float sna, csa;
  // approximate errors
  float errY = std::sqrt(seedC.getSigmaY2() + mParams->err2ABExtraY), errYFrac = errY * mRGHelper.ladderWidthInv();
  if (mFieldON) {
    seedC.getCircleParams(propagator->getNominalBz(), trcCircle, sna, csa);
  } else {
    seedC.getLineParams(trcLinPar, sna, csa);
  }
  float xCurr, yCurr;
  o2::math_utils::rotateZ(seedC.getX(), seedC.getY(), xCurr, yCurr, sna, csa); // lab X,Y
  float phi = std::atan2(yCurr, xCurr);                                        // RS: TODO : can we use fast atan2 here?
  // find approximate ladder and chip_in_ladder corresponding to this track extrapolation
  int nLad2Check = 0, ladIDguess = lr.getLadderID(phi), chipIDguess = lr.getChipID(seedC.getZ() + 0.5 * zDRStep);
  std::array<int, MaxLadderCand> lad2Check;
  nLad2Check = mFieldON ? findLaddersToCheckBOn(lrID, ladIDguess, trcCircle, errYFrac, lad2Check) : findLaddersToCheckBOff(lrID, ladIDguess, trcLinPar, errYFrac, lad2Check);

  for (int ilad = nLad2Check; ilad--;) {
    int ladID = lad2Check[ilad];
    const auto& lad = lr.ladders[ladID];

    // we assume that close chips on the same ladder will have close xyEdges, so it is enough to calculate track-chip crossing
    // coordinates xCross,yCross,zCross for this central chipIDguess, although we are going to check also neighbours
    float t = 1e9, xCross, yCross;
    const auto& chipC = lad.chips[chipIDguess];
    if (mFieldON) {
      chipC.xyEdges.circleCrossParam(trcCircle, t);
    } else {
      chipC.xyEdges.lineCrossParam(trcLinPar, t);
    }
    chipC.xyEdges.eval(t, xCross, yCross);
    float dx = xCross - xCurr, dy = yCross - yCurr, dst2 = dx * dx + dy * dy, dst = sqrtf(dst2);
    // Z-step sign depends on radius decreasing or increasing during the propagation
    float zCross = seedC.getZ() + seedC.getTgl() * (dst2 < 2 * (dx * xCurr + dy * yCurr) ? dst : -dst);

    for (int ich = -1; ich < 2; ich++) {
      int chipID = chipIDguess + ich;

      if (chipID < 0 || chipID >= static_cast<int>(lad.chips.size())) {
        continue;
      }
      if (lad.chips[chipID].zRange.isOutside(zCross, mParams->nABSigmaZ * errZ)) {
        continue;
      }
      const auto& clRange = itsChipClRefs.chipRefs[lad.chips[chipID].id];
      if (!clRange.getEntries()) {
        LOG(debug) << "No clusters in chip range";
        continue;
      }
      // track Y error in chip frame
      float errYcalp = errY * (csa * chipC.csAlp + sna * chipC.snAlp); // sigY_rotate(from alpha0 to alpha1) = sigY * cos(alpha1 - alpha0);
      float tolerZ = errZ * mParams->nABSigmaZ, tolerY = errYcalp * mParams->nABSigmaY;
      float yTrack = -xCross * chipC.snAlp + yCross * chipC.csAlp;                                           // track-chip crossing Y in chip frame
      if (!preselectChipClusters(chipSelClusters, clRange, itsChipClRefs, yTrack, zCross, tolerY, tolerZ)) { // select candidate clusters for this chip
        LOG(debug) << "No compatible clusters found";
        continue;
      }
      o2::track::TrackParCov trcLC = seedC;

      if (!trcLC.rotate(chipC.alp) || !trcLC.propagateTo(chipC.xRef, propagator->getNominalBz())) {
        LOG(debug) << " failed to rotate to alpha=" << chipC.alp << " or prop to X=" << chipC.xRef;
        // trcLC.print();
        break; // the chips of the ladder are practically on the same X and alpha
      }

      for (auto clID : chipSelClusters) {
        const auto& cls = mITSClustersArray[clID];
        auto chi2 = trcLC.getPredictedChi2(cls);
        if (chi2 > mParams->cutABTrack2ClChi2) {
          continue;
        }
        int lnkID = registerABTrackLink(ABSeed, trcLC, clID, seedID, lrID, ladID, chi2); // add new link with track copy
        if (lnkID > MinusOne) {
          auto& link = ABSeed.getLink(lnkID);
          link.update(cls);
          if (seedID >= MinusOne) {
            ABSeed.getLink(seedID).nDaughters++; // RS FIXME : do we need this?
          }
          if (lrID < ABSeed.lowestLayer) {
            ABSeed.lowestLayer = lrID; // update lowest layer reached
          }
        }
      }
    }
  }
  return 0;
}

//______________________________________________
void MatchTPCITS::accountForOverlapsAB(int lrSeed)
{
  // TODO
  LOG(warning) << "TODO";
}

//______________________________________________
int MatchTPCITS::findLaddersToCheckBOn(int ilr, int lad0, const o2::math_utils::CircleXYf_t& circle, float errYFrac,
                                       std::array<int, MatchTPCITS::MaxLadderCand>& lad2Check) const
{
  // check if ladder lad0 and at most +-MaxUpDnLadders around it are compatible with circular track of
  // r^2 = r2 and centered at xC,yC
  const auto& lr = mRGHelper.layers[ilr];
  int nacc = 0, jmp = 0;
  if (lr.ladders[lad0].xyEdges.seenByCircle(circle, errYFrac)) {
    lad2Check[nacc++] = lad0;
  }
  bool doUp = true, doDn = true;
  while ((doUp || doDn) && jmp++ < MaxUpDnLadders) {
    if (doUp) {
      int ldID = (lad0 + jmp) % lr.nLadders;
      if (lr.ladders[ldID].xyEdges.seenByCircle(circle, errYFrac)) {
        lad2Check[nacc++] = ldID;
      } else {
        doUp = false;
      }
    }
    if (doDn) {
      int ldID = lad0 - jmp;
      if (ldID < 0) {
        ldID += lr.nLadders;
      }
      if (lr.ladders[ldID].xyEdges.seenByCircle(circle, errYFrac)) {
        lad2Check[nacc++] = ldID;
      } else {
        doDn = false;
      }
    }
  }
  return nacc;
}

//______________________________________________
int MatchTPCITS::findLaddersToCheckBOff(int ilr, int lad0, const o2::math_utils::IntervalXYf_t& trcLinPar, float errYFrac,
                                        std::array<int, MatchTPCITS::MaxLadderCand>& lad2Check) const
{
  // check if ladder lad0 and at most +-MaxUpDnLadders around it are compatible with linear track

  const auto& lr = mRGHelper.layers[ilr];
  int nacc = 0, jmp = 0;
  if (lr.ladders[lad0].xyEdges.seenByLine(trcLinPar, errYFrac)) {
    lad2Check[nacc++] = lad0;
  }
  bool doUp = true, doDn = true;
  while ((doUp || doDn) && jmp++ < MaxUpDnLadders) {
    if (doUp) {
      int ldID = (lad0 + jmp) % lr.nLadders;
      if (lr.ladders[ldID].xyEdges.seenByLine(trcLinPar, errYFrac)) {
        lad2Check[nacc++] = ldID;
      } else {
        doUp = false;
      }
    }
    if (doDn) {
      int ldID = lad0 - jmp;
      if (ldID < 0) {
        ldID += lr.nLadders;
      }
      if (lr.ladders[ldID].xyEdges.seenByLine(trcLinPar, errYFrac)) {
        lad2Check[nacc++] = ldID;
      } else {
        doDn = false;
      }
    }
  }
  return nacc;
}

//______________________________________________
int MatchTPCITS::registerABTrackLink(TPCABSeed& ABSeed, const o2::track::TrackParCov& trc, int clID, int parentID, int lr, int laddID, float chi2Cl)
{
  // registers new ABLink on the layer, assigning provided kinematics. The link will be registered in a
  // way preserving the quality ordering of the links on the layer
  int lnkID = ABSeed.getNLinks(), nextID = ABSeed.firstInLr[lr], nc = 1 + (parentID > MinusOne ? ABSeed.getLink(parentID).nContLayers : 0);
  float chi2 = chi2Cl + (parentID > MinusOne ? ABSeed.getLink(parentID).chi2 : 0.);
  // LOG(info) << "Reg on lr "  << lr << " nc = " << nc << " chi2cl=" << chi2Cl << " -> " << chi2; // RSTMP

  if (ABSeed.firstInLr[lr] == MinusOne) { // no links on this layer yet
    ABSeed.firstInLr[lr] = lnkID;
    ABSeed.addLink(trc, clID, parentID, MinusOne, lr, nc, laddID, chi2);
    return lnkID;
  }
  // add new link sorting links of this layer in quality

  int count = 0, topID = MinusOne;
  do {
    auto& nextLink = ABSeed.getLink(nextID);
    count++;
    bool newIsBetter = parentID <= MinusOne ? isBetter(chi2, nextLink.chi2) : isBetter(ABSeed.getLink(parentID).chi2NormPredict(chi2Cl), nextLink.chi2Norm());
    if (newIsBetter) {                          // need to insert new link before nextLink
      if (count < mParams->maxABLinksOnLayer) { // will insert in front of nextID
        ABSeed.addLink(trc, clID, parentID, nextID, lr, nc, laddID, chi2);
        if (topID == MinusOne) {        // are we comparing new link with best link on the layer?
          ABSeed.firstInLr[lr] = lnkID; // flag as best on the layer
        } else {
          ABSeed.getLink(topID).nextOnLr = lnkID; // point from previous one
        }
        return lnkID;
      } else { // max number of candidates reached, will overwrite the last one
        nextLink = ABTrackLink(trc, clID, parentID, MinusOne, lr, nc, laddID, chi2);
        return nextID;
      }
    }
    topID = nextID;
    nextID = nextLink.nextOnLr;
  } while (nextID > MinusOne);
  // new link is worse than all others, add it only if there is a room to expand
  if (count < mParams->maxABLinksOnLayer) {
    ABSeed.addLink(trc, clID, parentID, MinusOne, lr, nc, laddID, chi2);
    if (topID > MinusOne) {
      ABSeed.getLink(topID).nextOnLr = lnkID; // point from previous one
    }
    return lnkID;
  }
  return MinusOne; // link to be ignored
}

//______________________________________________
float MatchTPCITS::correctTPCTrack(o2::track::TrackParCov& trc, const TrackLocTPC& tTPC, const InteractionCandidate& cand) const
{
  // Correct the track copy trc of the working TPC track tTPC in continuous RO mode for the assumed interaction time
  // return extra uncertainty in Z due to the interaction time incertainty
  // TODO: at the moment, apply simple shift, but with Z-dependent calibration we may
  // need to do corrections on TPC cluster level and refit
  if (tTPC.constraint == TrackLocTPC::Constrained) {
    return 0.f;
  }
  auto tpcTrOrig = mTPCTracksArray[tTPC.sourceID];
  float timeIC = cand.tBracket.mean();
  float driftErr = cand.tBracket.delta() * mTPCBin2Z;

  // we use this for refit, at the moment it is not done ...
  /*
  {
    float r = std::sqrt(trc.getX()*trc.getX() + trc.getY()*trc.getY());
    float chi2 = 0;
    bool res = refitTPCInward(trc, chi2, r, tTPC.sourceID(), timeIC );
    if (!res) {
      return -1;
    }
    float xTgt;
    auto propagator = o2::base::Propagator::Instance();
    if (!trc.getXatLabR(r, xTgt, propagator->getNominalBz(), o2::track::DirInward) ||
 !propagator->PropagateToXBxByBz(trc, xTgt, MaxSnp, 2., mUseMatCorrFlag)) {
      return -1;
    }
  }
  */
  // if interaction time precedes the initial assumption on t0 (i.e. timeIC < timeTrc),
  // the track actually was drifting longer, i.e. tracks should be shifted closer to the CE
  float dDrift = (timeIC - tTPC.time0) * mTPCBin2Z;

  // float zz = tTPC.getZ() + (tpcTrOrig.hasASideClustersOnly() ? dDrift : -dDrift);                                 // tmp
  // LOG(info) << "CorrTrack Z=" << trc.getZ() << " (zold= " << zz << ") at TIC= " << timeIC << " Ttr= " << tTPC.time0; // tmp

  // we use this w/o refit
  //
  {
    trc.setZ(tTPC.getZ() + (tTPC.constraint == TrackLocTPC::ASide ? dDrift : -dDrift));
  }
  //
  trc.setCov(trc.getSigmaZ2() + driftErr * driftErr, o2::track::kSigZ2);

  return driftErr;
}

//______________________________________________
void MatchTPCITS::fillClustersForAfterBurner(int rofStart, int nROFs, ITSChipClustersRefs& itsChipClRefs)
{
  // Prepare unused clusters of given ROFs range for matching in the afterburner
  // Note: normally only 1 ROF needs to be filled (nROFs==1 ) unless we want
  // to account for interaction on the boundary of 2 rofs, which then may contribute to both ROFs.
  int first = mITSClusterROFRec[rofStart].getFirstEntry(), last = first;
  for (int ir = nROFs; ir--;) {
    last += mITSClusterROFRec[rofStart + ir].getNEntries();
  }
  itsChipClRefs.clear();
  auto& idxSort = itsChipClRefs.clusterID;
  for (int icl = first; icl < last; icl++) {
    if (mABClusterLinkIndex[icl] != MinusTen) { // clusters with MinusOne are used in main matching
      idxSort.push_back(icl);
    }
  }
  // sort in chip, Z
  const auto& clusArr = mITSClustersArray;
  std::sort(idxSort.begin(), idxSort.end(), [&clusArr](int i, int j) {
    const auto &clI = clusArr[i], &clJ = clusArr[j];
    if (clI.getSensorID() < clJ.getSensorID()) {
      return true;
    }
    if (clI.getSensorID() == clJ.getSensorID()) {
      return clI.getZ() < clJ.getZ();
    }
    return false;
  });

  int ncl = idxSort.size();
  int lastSens = -1, nClInSens = 0;
  ClusRange* chipClRefs = nullptr;
  for (int icl = 0; icl < ncl; icl++) {
    const auto& clus = mITSClustersArray[idxSort[icl]];
    int sens = clus.getSensorID();
    if (sens != lastSens) {
      if (chipClRefs) { // finalize chip reference
        chipClRefs->setEntries(nClInSens);
        nClInSens = 0;
      }
      chipClRefs = &itsChipClRefs.chipRefs[(lastSens = sens)];
      chipClRefs->setFirstEntry(icl);
    }
    nClInSens++;
  }
  if (chipClRefs) {
    chipClRefs->setEntries(nClInSens); // finalize last chip reference
  }
}

//______________________________________________
void MatchTPCITS::setITSTimeBiasInBC(int n)
{
  mITSTimeBiasInBC = n;
  mITSTimeBiasMUS = mITSTimeBiasInBC * o2::constants::lhc::LHCBunchSpacingNS * 1e-3;
}

//______________________________________________
void MatchTPCITS::setITSROFrameLengthMUS(float fums)
{
  mITSROFrameLengthMUS = fums;
  mITSTimeResMUS = mITSROFrameLengthMUS / std::sqrt(12.f);
  mITSROFrameLengthMUSInv = 1. / mITSROFrameLengthMUS;
  mITSROFrameLengthInBC = std::max(1, int(mITSROFrameLengthMUS / (o2::constants::lhc::LHCBunchSpacingNS * 1e-3)));
}

//______________________________________________
void MatchTPCITS::setITSROFrameLengthInBC(int nbc)
{
  mITSROFrameLengthInBC = nbc;
  mITSROFrameLengthMUS = nbc * o2::constants::lhc::LHCBunchSpacingNS * 1e-3;
  mITSTimeResMUS = mITSROFrameLengthMUS / std::sqrt(12.f);
  mITSROFrameLengthMUSInv = 1. / mITSROFrameLengthMUS;
}

//___________________________________________________________________
void MatchTPCITS::setBunchFilling(const o2::BunchFilling& bf)
{
  mBunchFilling = bf;
  // find closest (from above) filled bunch
  int minBC = bf.getFirstFilledBC(), maxBC = bf.getLastFilledBC();
  if (minBC < 0 && mUseBCFilling) {
    if (mUseBCFilling) {
      mUseBCFilling = false;
      LOG(warning) << "Disabling match validation by BunchFilling as no interacting bunches found";
    }
    return;
  }
  int bcAbove = minBC;
  for (int i = o2::constants::lhc::LHCMaxBunches; i--;) {
    if (bf.testBC(i)) {
      bcAbove = i;
    }
    mClosestBunchAbove[i] = bcAbove;
  }
  int bcBelow = maxBC;
  for (int i = 0; i < o2::constants::lhc::LHCMaxBunches; i++) {
    if (bf.testBC(i)) {
      bcBelow = i;
    }
    mClosestBunchBelow[i] = bcBelow;
  }
}

//___________________________________________________________________
MatchTPCITS::BracketIR MatchTPCITS::tBracket2IRBracket(const BracketF tbrange)
{
  // convert time bracket to IR bracket
  o2::InteractionRecord irMin(mStartIR), irMax(mStartIR);
  if (tbrange.getMin() > 0) {
    irMin += o2::InteractionRecord(tbrange.getMin() * 1000.f); // time in ns is needed
  }
  irMax += o2::InteractionRecord(tbrange.getMax() * 1000.f);
  irMax++; // to account for rounding
  int bc = mClosestBunchAbove[irMin.bc];
  if (bc < irMin.bc) {
    irMin.orbit++;
  }
  irMin.bc = bc;
  bc = mClosestBunchBelow[irMax.bc];
  if (bc > irMax.bc) {
    if (irMax.orbit == 0) {
      bc = 0;
    } else {
      irMax.orbit--;
    }
  }
  irMax.bc = bc;
  return {irMin, irMax};
}

//______________________________________________
void MatchTPCITS::removeTPCfromITS(int tpcID, int itsID)
{
  ///< remove reference to tpcID track from itsID track matches
  auto& tITS = mITSWork[itsID];
  if (isValidatedITS(tITS)) {
    return;
  }
  int topID = MinusOne, next = tITS.matchID; // ITS MatchRecord
  while (next > MinusOne) {
    auto& rcITS = mMatchRecordsITS[next];
    if (rcITS.partnerID == tpcID) {
      if (topID < 0) {
        tITS.matchID = rcITS.nextRecID;
      } else {
        mMatchRecordsITS[topID].nextRecID = rcITS.nextRecID;
      }
      return;
    }
    topID = next;
    next = rcITS.nextRecID;
  }
}

//______________________________________________
void MatchTPCITS::removeITSfromTPC(int itsID, int tpcID)
{
  ///< remove reference to itsID track from matches of tpcID track
  auto& tTPC = mTPCWork[tpcID];
  if (isValidatedTPC(tTPC)) {
    return;
  }
  int topID = MinusOne, next = tTPC.matchID;
  while (next > MinusOne) {
    auto& rcTPC = mMatchRecordsTPC[next];
    if (rcTPC.partnerID == itsID) {
      if (topID < 0) {
        tTPC.matchID = rcTPC.nextRecID;
      } else {
        mMatchRecordsTPC[topID].nextRecID = rcTPC.nextRecID;
      }
      return;
    }
    topID = next;
    next = rcTPC.nextRecID;
  }
}

//______________________________________________
void MatchTPCITS::flagUsedITSClusters(const o2::its::TrackITS& track)
{
  // flag clusters used by this track
  int clEntry = track.getFirstClusterEntry();
  for (int icl = track.getNumberOfClusters(); icl--;) {
    mABClusterLinkIndex[mITSTrackClusIdx[clEntry++]] = MinusTen;
  }
}

//__________________________________________________________
int MatchTPCITS::preselectChipClusters(std::vector<int>& clVecOut, const ClusRange& clRange, const ITSChipClustersRefs& itsChipClRefs,
                                       float trackY, float trackZ, float tolerY, float tolerZ) const
{
  clVecOut.clear();
  int icID = clRange.getFirstEntry();
  for (int icl = clRange.getEntries(); icl--;) { // note: clusters within a chip are sorted in Z
    int clID = itsChipClRefs.clusterID[icID++];  // so, we go in clusterID increasing direction
    const auto& cls = mITSClustersArray[clID];
    float dz = cls.getZ() - trackZ;
    LOG(debug) << "cl" << icl << '/' << clID << " "
               << " dZ: " << dz << " [" << tolerZ << "| dY: " << trackY - cls.getY() << " [" << tolerY << "]";
    if (dz > tolerZ) {
      float clsZ = cls.getZ();
      LOG(debug) << "Skip the rest since " << trackZ << " < " << clsZ << "\n";
      break;
    } else if (dz < -tolerZ) {
      LOG(debug) << "Skip cluster dz=" << dz << " Ztr=" << trackZ << " zCl=" << cls.getZ();
      continue;
    }
    if (fabs(trackY - cls.getY()) > tolerY) {
      LOG(debug) << "Skip cluster dy= " << trackY - cls.getY() << " Ytr=" << trackY << " yCl=" << cls.getY();
      continue;
    }
    clVecOut.push_back(clID);
  }
  return clVecOut.size();
}

//__________________________________________________________
void MatchTPCITS::reportSizes(pmr::vector<o2::dataformats::TrackTPCITS>& matchedTracks,
                              pmr::vector<o2::itsmft::TrkClusRef>& ABTrackletRefs,
                              pmr::vector<int>& ABTrackletClusterIDs,
                              pmr::vector<o2::MCCompLabel>& matchLabels,
                              pmr::vector<o2::MCCompLabel>& ABTrackletLabels,
                              pmr::vector<o2::dataformats::Triplet<float, float, float>>& calib)
{
  size_t sizTotShm = 0, capTotShm = 0, sizTot = 0, capTot = 0, siz = 0, cap = 0, cnt = 0, cntCap = 0;
  {
    siz = matchedTracks.size() * sizeof(o2::dataformats::TrackTPCITS);
    cap = matchedTracks.capacity() * sizeof(o2::dataformats::TrackTPCITS);
    sizTotShm += siz;
    capTotShm += cap;
    LOGP(info, "Size SHM, matchedTracks         : size {:9} cap {:9}", siz, cap);
    //
    siz = ABTrackletRefs.size() * sizeof(o2::itsmft::TrkClusRef);
    cap = ABTrackletRefs.capacity() * sizeof(o2::itsmft::TrkClusRef);
    sizTotShm += siz;
    capTotShm += cap;
    LOGP(info, "Size SHM, ABTrackletRefs        : size {:9} cap {:9}", siz, cap);
    //
    siz = ABTrackletClusterIDs.size() * sizeof(int);
    cap = ABTrackletClusterIDs.capacity() * sizeof(int);
    sizTotShm += siz;
    capTotShm += cap;
    LOGP(info, "Size SHM, ABTrackletClusterIDs  : size {:9} cap {:9}", siz, cap);
    //
    siz = matchLabels.size() * sizeof(o2::MCCompLabel);
    cap = matchLabels.capacity() * sizeof(o2::MCCompLabel);
    sizTotShm += siz;
    capTotShm += cap;
    LOGP(info, "Size SHM, matchLabels           : size {:9} cap {:9}", siz, cap);
    //
    siz = ABTrackletLabels.size() * sizeof(o2::MCCompLabel);
    cap = ABTrackletLabels.capacity() * sizeof(o2::MCCompLabel);
    sizTotShm += siz;
    capTotShm += cap;
    LOGP(info, "Size SHM, ABTrackletLabels      : size {:9} cap {:9}", siz, cap);
    //
    siz = calib.size() * sizeof(o2::dataformats::Triplet<float, float, float>);
    cap = calib.capacity() * sizeof(o2::dataformats::Triplet<float, float, float>);
    sizTotShm += siz;
    capTotShm += cap;
    LOGP(info, "Size SHM, calib                 : size {:9} cap {:9}", siz, cap);
  }
  {
    siz = mITSClustersArray.size() * sizeof(ITSCluster);
    cap = mITSClustersArray.capacity() * sizeof(ITSCluster);
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, mITSClustersArray     : size {:9} cap {:9}", siz, cap);
    //
    siz = mMatchRecordsTPC.size() * sizeof(MatchRecord);
    cap = mMatchRecordsTPC.capacity() * sizeof(MatchRecord);
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, mMatchRecordsTPC      : size {:9} cap {:9}", siz, cap);
    //
    siz = mMatchRecordsITS.size() * sizeof(MatchRecord);
    cap = mMatchRecordsITS.capacity() * sizeof(MatchRecord);
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, mMatchRecordsITS      : size {:9} cap {:9}", siz, cap);
    //
    siz = mITSROFTimes.size() * sizeof(BracketF);
    cap = mITSROFTimes.capacity() * sizeof(BracketF);
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, mITSROFTimes          : size {:9} cap {:9}", siz, cap);
    //
    siz = mTPCWork.size() * sizeof(TrackLocTPC);
    cap = mTPCWork.capacity() * sizeof(TrackLocTPC);
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, mTPCWork              : size {:9} cap {:9}", siz, cap);
    //
    siz = mITSWork.size() * sizeof(TrackLocITS);
    cap = mITSWork.capacity() * sizeof(TrackLocITS);
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, mITSWork              : size {:9} cap {:9}", siz, cap);
    //
    siz = mWinnerChi2Refit.size() * sizeof(float);
    cap = mWinnerChi2Refit.capacity() * sizeof(float);
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, mWinnerChi2Refit      : size {:9} cap {:9}", siz, cap);
    //
    siz = mTPCABSeeds.size() * sizeof(float);
    cap = mTPCABSeeds.capacity() * sizeof(float);
    cnt = 0;
    cntCap = 0;
    for (const auto& a : mTPCABSeeds) {
      siz += a.sizeInternal();
      cap += a.capInternal();
      cnt += a.getNLinks();
      cntCap += a.getNLinks();
    }
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, mTPCABSeeds           : size {:9} cap {:9} | internals size:{}/capacity:{} for {} elements", siz, cap, cnt, cntCap, mTPCABSeeds.size());
    //
    siz = mTPCABIndexCache.size() * sizeof(int);
    cap = mTPCABIndexCache.capacity() * sizeof(int);
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, mTPCABIndexCache      : size {:9} cap {:9}", siz, cap);
    //
    siz = mABWinnersIDs.size() * sizeof(int);
    cap = mABWinnersIDs.capacity() * sizeof(int);
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, mABWinnersIDs         : size {:9} cap {:9}", siz, cap);
    //
    siz = mABClusterLinkIndex.size() * sizeof(int);
    cap = mABClusterLinkIndex.capacity() * sizeof(int);
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, mABClusterLinkIndex   : size {:9} cap {:9}", siz, cap);
    //
    for (int is = 0; is < o2::constants::math::NSectors; is++) {
      siz += mTPCSectIndexCache[is].size() * sizeof(int);
      cap += mTPCSectIndexCache[is].capacity() * sizeof(int);
    }
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, mTPCSectIndexCache    : size {:9} cap {:9}", siz, cap);
    //
    for (int is = 0; is < o2::constants::math::NSectors; is++) {
      siz += mITSSectIndexCache[is].size() * sizeof(int);
      cap += mITSSectIndexCache[is].capacity() * sizeof(int);
    }
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, mITSSectIndexCache    : size {:9} cap {:9}", siz, cap);
    //
    for (int is = 0; is < o2::constants::math::NSectors; is++) {
      siz += mTPCTimeStart[is].size() * sizeof(int);
      cap += mTPCTimeStart[is].capacity() * sizeof(int);
    }
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, mTPCTimeStart         : size {:9} cap {:9}", siz, cap);
    //
    for (int is = 0; is < o2::constants::math::NSectors; is++) {
      siz += mITSTimeStart[is].size() * sizeof(int);
      cap += mITSTimeStart[is].capacity() * sizeof(int);
    }
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, mITSTimeStart         : size {:9} cap {:9}", siz, cap);
    //
    siz = mITSTrackROFContMapping.size() * sizeof(int);
    cap = mITSTrackROFContMapping.capacity() * sizeof(int);
    sizTot += siz;
    capTot += cap;
    LOGP(info, "Size RSS, ITSTrackROFContMapping: size {:9} cap {:9}", siz, cap);
  }
  LOGP(info, "TotalSizes/Capacities: SHM: {}/{}  Heap: {}/{}", sizTotShm, capTotShm, sizTot, capTot);
}

//__________________________________________________________
void MatchTPCITS::setNThreads(int n)
{
#ifdef WITH_OPENMP
  mNThreads = n > 0 ? n : 1;
#else
  LOG(warning) << "Multithreading is not supported, imposing single thread";
  mNThreads = 1;
#endif
  mABLinksPool.threadPool.resize(mNThreads);
  TPCABSeed::gLinksPool = &mABLinksPool;
}

//<<============================= AfterBurner for TPC-track / ITS cluster matching ===================<<

#ifdef _ALLOW_DEBUG_TREES_
//______________________________________________
void MatchTPCITS::setDebugFlag(UInt_t flag, bool on)
{
  ///< set debug stream flag
  if (on) {
    mDBGFlags |= flag;
  } else {
    mDBGFlags &= ~flag;
  }
}

//_________________________________________________________
void MatchTPCITS::fillTPCITSmatchTree(int itsID, int tpcID, int rejFlag, float chi2, float tCorr)
{
  ///< fill debug tree for ITS TPC tracks matching check

  mTimer[SWDBG].Start(false);

  auto& trackITS = mITSWork[itsID];
  auto& trackTPC = mTPCWork[tpcID];
  if (chi2 < 0.) { // need to recalculate
    chi2 = getPredictedChi2NoZ(trackITS, trackTPC);
  }
  (*mDBGOut) << "match"
             << "tf=" << mTFCount << "chi2Match=" << chi2 << "its=" << trackITS << "tpc=" << trackTPC << "tcorr=" << tCorr;
  if (mMCTruthON) {
    (*mDBGOut) << "match"
               << "itsLbl=" << mITSLblWork[itsID] << "tpcLbl=" << mTPCLblWork[tpcID];
  }
  int tb = mTPCTracksArray[trackTPC.sourceID].getTime0() * mNTPCOccBinLengthInv;
  float mltTPC = tb < 0 ? mTBinClOcc[0] : (tb >= mTBinClOcc.size() ? mTBinClOcc.back() : mTBinClOcc[tb]);
  (*mDBGOut) << "match"
             << "rejFlag=" << rejFlag
             << "multTPC=" << tb
             << "multITSTr=" << mITSTrackROFRec[trackITS.roFrame]
             << "multITSCl=" << mITSClusterROFRec[trackITS.roFrame]
             << "\n";

  mTimer[SWDBG].Stop();
}

//______________________________________________
void MatchTPCITS::dumpWinnerMatches()
{
  ///< write winner matches into debug tree

  mTimer[SWDBG].Start(false);

  LOG(info) << "Dumping debug tree for winner matches";
  for (int iits = 0; iits < int(mITSWork.size()); iits++) {
    auto& tITS = mITSWork[iits];
    if (isDisabledITS(tITS)) {
      continue;
    }
    auto& itsMatchRec = mMatchRecordsITS[tITS.matchID];
    int itpc = itsMatchRec.partnerID;
    auto& tTPC = mTPCWork[itpc];

    (*mDBGOut) << "matchWin"
               << "tf=" << mTFCount << "chi2Match=" << itsMatchRec.chi2 << "chi2Refit=" << mWinnerChi2Refit[iits] << "its=" << tITS << "tpc=" << tTPC;

    if (mMCTruthON) {
      (*mDBGOut) << "matchWin"
                 << "itsLbl=" << mITSLblWork[iits] << "tpcLbl=" << mTPCLblWork[itpc];
    }
    int tb = mTPCTracksArray[tTPC.sourceID].getTime0() * mNTPCOccBinLengthInv;
    float mltTPC = tb < 0 ? mTBinClOcc[0] : (tb >= mTBinClOcc.size() ? mTBinClOcc.back() : mTBinClOcc[tb]);
    (*mDBGOut) << "matchWin"
               << "multTPC=" << mTPCRefitter->getParam()->GetUnscaledMult(mTPCTracksArray[tTPC.sourceID].getTime0())
               << "multITSTr=" << mITSTrackROFRec[tITS.roFrame]
               << "multITSCl=" << mITSClusterROFRec[tITS.roFrame]
               << "\n";
  }
  mTimer[SWDBG].Stop();
}

#endif
