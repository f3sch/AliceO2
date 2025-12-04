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

#include <cassert>
#include <memory>

#include <oneapi/tbb/task_arena.h>

#include "ITSMFTBase/DPLAlpideParam.h"
#include "ITSBase/GeometryTGeo.h"

#include "ITSReconstruction/FastMultEstConfig.h"
#include "ITSReconstruction/FastMultEst.h"

#include "ITStracking/Configuration.h"
#include "ITStracking/TrackingConfigParam.h"
#include "ITStracking/TrackingInterface.h"

#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/PhysTrigger.h"
#include "DataFormatsTRD/TriggerRecord.h"
#include "CommonDataFormat/IRFrame.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "ITStracking/BoundedAllocator.h"
#include "Framework/InputRecordWalker.h"
#include "Framework/DataRefUtils.h"
#include "Framework/DeviceSpec.h"

using namespace o2::framework;
using namespace o2::its;

void ITSTrackingInterface::initialise()
{
  // get parameters
  const auto& trackConf = o2::its::TrackerParamConfig::Instance();
  if (auto parmode = (TrackingMode::Type)trackConf.trackingMode; mMode == TrackingMode::Unset || (parmode != TrackingMode::Unset && mMode != parmode)) {
    LOGP(info, "Tracking mode overwritten by configurable params from {} to {}", TrackingMode::toString(mMode), TrackingMode::toString(parmode));
    mMode = parmode;
  }
  auto iterations = TrackingMode::getRecoIterations(mMode);
  LOGP(info, "Initializing tracker in {} phase reconstruction with {} passes for tracking", TrackingMode::toString(mMode), iterations.size());
  mTracker->setParameters(iterations);

  if (mMode == TrackingMode::Cosmics) {
    mRunVertexer = false;
    mCosmicsProcessing = true;
    LOGP(info, "Cosmic mode enabled, will skip vertexing");
  }

  // threading
  bool clamped{false};
  int nThreads = trackConf.nThreads;
  if (nThreads > 0) {
    const int hw = std::thread::hardware_concurrency();
    const int maxThreads = (hw == 0 ? 1 : hw);
    nThreads = std::clamp(nThreads, 1, maxThreads);
    clamped = trackConf.nThreads > maxThreads;
  }
  LOGP(info, "Tracker has {} thread(s){}", nThreads, (clamped) ? " (clamped)" : "");
  mTaskArena = std::make_shared<tbb::task_arena>(std::abs(nThreads));
  mTracker->setNThreads(trackConf.nThreads, mTaskArena);

  // prepare data filter
  for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
    mFilter.emplace_back("compClusters", "ITS", "COMPCLUSTERS", iLayer, Lifetime::Timeframe);
    mFilter.emplace_back("patterns", "ITS", "PATTERNS", iLayer, Lifetime::Timeframe);
    mFilter.emplace_back("ROframe", "ITS", "CLUSTERSROF", iLayer, Lifetime::Timeframe);
    if (mIsMC) {
      mFilter.emplace_back("itsmclabels", "ITS", "CLUSTERSMCTR", iLayer, Lifetime::Timeframe);
    }
  }
}

void ITSTrackingInterface::run(framework::ProcessingContext& pc)
{
  if (static bool doneOnce{false}; !doneOnce) {
    doneOnce = true;

    // prepare rof lookup table(s)
    // has to be done here to ensure we get the right number of HB per TF
    const int nOrbitsPerTF = o2::base::GRPGeomHelper::getNHBFPerTF();
    TimeFrameN::ROFOverlapTableN rofTable;
    TimeFrameN::ROFVertexLookupTableN vtxTable;
    TimeFrameN::ROFTimeSliceTableN sliceTable;
    const auto& par = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>::Instance();
    const auto& trackParams = mTracker->getParameters();
    for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
      const int nROFsPerOrbit = o2::constants::lhc::LHCMaxBunches / par.getROFLengthInBC(iLayer);
      const LayerTiming timing{.mNROFsTF = (nROFsPerOrbit * nOrbitsPerTF), .mROFLength = par.getROFLengthInBC(iLayer), .mROFDelay = par.getROFDelayInBC(iLayer), .mROFDelta = trackParams[0].params.DeltaROF[iLayer]};
      rofTable.defineLayer(iLayer, timing);
      vtxTable.defineLayer(iLayer, timing);
      sliceTable.defineLayer(iLayer, timing);
    }
    rofTable.init();
    mTimeFrame->setROFOverlapTable(rofTable);
    vtxTable.init();
    mTimeFrame->setROFVertexLookupTable(vtxTable);
    sliceTable.init(trackParams[0].params.NTimeSlices);
    mTimeFrame->setROFTimeSliceTable(sliceTable);
  }

  std::array<gsl::span<const itsmft::CompClusterExt>, NLayers> compClusters;
  std::array<gsl::span<const unsigned char>, NLayers> patterns;
  std::array<gsl::span<const itsmft::ROFRecord>, NLayers> rofsinput;
  std::array<const dataformats::MCTruthContainer<MCCompLabel>*, NLayers> labels{};

  // filter input and compose
  for (const DataRef& ref : framework::InputRecordWalker{pc.inputs(), mFilter}) {
    auto const* dh = DataRefUtils::getHeader<o2::header::DataHeader*>(ref);
    if (framework::DataRefUtils::match(ref, {"compClusters", framework::ConcreteDataTypeMatcher{"ITS", "COMPCLUSTERS"}})) {
      compClusters[dh->subSpecification] = pc.inputs().get<gsl::span<o2::itsmft::CompClusterExt>>(ref);
    }
    if (framework::DataRefUtils::match(ref, {"patterns", framework::ConcreteDataTypeMatcher{"ITS", "PATTERNS"}})) {
      patterns[dh->subSpecification] = pc.inputs().get<gsl::span<unsigned char>>(ref);
    }
    if (framework::DataRefUtils::match(ref, {"ROframes", framework::ConcreteDataTypeMatcher{"ITS", "CLUSTERSROF"}})) {
      rofsinput[dh->subSpecification] = pc.inputs().get<gsl::span<o2::itsmft::ROFRecord>>(ref);
    }
    if (framework::DataRefUtils::match(ref, {"itsmclabels", framework::ConcreteDataTypeMatcher{"ITS", "CLUSTERSMCTR"}})) {
      labels[dh->subSpecification] = pc.inputs().get<const dataformats::MCTruthContainer<MCCompLabel>*>(ref).release();
    }
  }
  const auto& alpParams = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>::Instance();
  for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
    LOGP(info, "ITSTracker:{} pulled {} clusters, {} RO frames", iLayer, compClusters[iLayer].size(), rofsinput[iLayer].size());
    if (compClusters[iLayer].empty()) {
      LOGP(warn, " -> received no processable data on layer {}", iLayer);
    }
    if (mIsMC) {
      LOG(info) << " -> " << labels[iLayer]->getIndexedSize() << " MC label objects";
    }
  }

  // trigger
  gsl::span<const o2::itsmft::PhysTrigger> physTriggers;
  std::vector<o2::itsmft::PhysTrigger> fromTRD;
  if (mUseTriggers == 2) { // use TRD triggers
    o2::InteractionRecord ir{0, pc.services().get<o2::framework::TimingInfo>().firstTForbit};
    auto trdTriggers = pc.inputs().get<gsl::span<o2::trd::TriggerRecord>>("phystrig");
    for (const auto& trig : trdTriggers) {
      if (trig.getBCData() >= ir && trig.getNumberOfTracklets()) {
        ir = trig.getBCData();
        fromTRD.emplace_back(o2::itsmft::PhysTrigger{ir, 0});
      }
    }
    physTriggers = gsl::span<const o2::itsmft::PhysTrigger>(fromTRD.data(), fromTRD.size());
  } else if (mUseTriggers == 1) { // use Phys triggers from ITS stream
    physTriggers = pc.inputs().get<gsl::span<o2::itsmft::PhysTrigger>>("phystrig");
  }

  // prepare output data
  auto& irFrames = pc.outputs().make<std::vector<o2::dataformats::IRFrame>>(Output{"ITS", "IRFRAMES", 0});
  // FIXME
  // irFrames.reserve(trackROFvec.size());
  auto& allClusIdx = pc.outputs().make<std::vector<int>>(Output{"ITS", "TRACKCLSID", 0});
  auto& allTracks = pc.outputs().make<std::vector<o2::its::TrackITS>>(Output{"ITS", "TRACKS", 0});
  auto& vertices = pc.outputs().make<std::vector<Vertex>>(Output{"ITS", "VERTICES", 0});

  // MC
  static pmr::vector<o2::MCCompLabel> dummyMCLabTracks, dummyMCLabVerts;
  static pmr::vector<float> dummyMCPurVerts;
  auto& allTrackLabels = mIsMC ? pc.outputs().make<std::vector<o2::MCCompLabel>>(Output{"ITS", "TRACKSMCTR", 0}) : dummyMCLabTracks;
  auto& allVerticesLabels = mIsMC ? pc.outputs().make<std::vector<o2::MCCompLabel>>(Output{"ITS", "VERTICESMCTR", 0}) : dummyMCLabVerts;
  auto& allVerticesPurities = mIsMC ? pc.outputs().make<std::vector<float>>(Output{"ITS", "VERTICESMCPUR", 0}) : dummyMCPurVerts;

  bool continuous = o2::base::GRPGeomHelper::instance().getGRPECS()->isDetContinuousReadOut(o2::detectors::DetID::ITS);
  LOG(info) << "ITSTracker RO: continuous=" << continuous;

  if (mOverrideBeamEstimation) {
    mTimeFrame->setMeanVertex(mMeanVertex, TrackerParamConfig::Instance().seedingMeanVertexExtraErr2);
  }

  mTracker->setBz(o2::base::Propagator::Instance()->getNominalBz());

  for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
    gsl::span<const unsigned char>::iterator pattIt = patterns[iLayer].begin();
    loadROF(rofsinput[iLayer], compClusters[iLayer], pattIt, iLayer, labels[iLayer]);
  }

  auto logger = [&](const std::string& s) { LOG(info) << s; };
  auto fatalLogger = [&](const std::string& s) { LOG(fatal) << s; };
  auto errorLogger = [&](const std::string& s) { LOG(error) << s; };

  FastMultEst multEst; // mult estimator
  // std::array<std::vector<uint8_t>, NLayers> processingMask, processUPCMask;
  // std::array<int, NLayers> cutRandomMult{};
  // for (int iLayer = 0; iLayer < NLayers; ++iLayer) {
  //   cutRandomMult[iLayer] = int(rofsinput[iLayer].size()) - multEst.selectROFs(rofsinput[iLayer], compClusters[iLayer], physTriggers, processingMask[iLayer]);
  //   processUPCMask[iLayer].resize(processingMask.size(), 0u);
  // }
  // int cutVertexMult{0}, cutUPCVertex{0};
  // mTimeFrame->setMultiplicityCutMask(processingMask);
  // float vertexerElapsedTime{0.f};
  // if (mRunVertexer) {
  // Run seeding vertexer
  // vertexerElapsedTime = mVertexer->clustersToVertices(logger);
  // } else { // cosmics
  //   mTimeFrame->resetRofPV();
  // }
  const auto& multEstConf = FastMultEstConfig::Instance(); // parameters for mult estimation and cuts
  gsl::span<const std::pair<MCCompLabel, float>> vMCRecInfo;
  gsl::span<const MCCompLabel> vMCContLabels;
  // for (auto iRof{0}; iRof < rofsinput[1].size(); ++iRof) {
  //   bounded_vector<Vertex> vtxVecLoc;
  //   auto& vtxROF = vertROFvec.emplace_back(rofsinput[1][iRof]);
  //   vtxROF.setFirstEntry(vertices.size());
  //   if (mRunVertexer) {
  //     auto vtxSpan = mTimeFrame->getPrimaryVertices(iRof);
  //     if (mIsMC) {
  //       vMCRecInfo = mTimeFrame->getPrimaryVerticesMCRecInfo(iRof);
  //     }
  //     if (o2::its::TrackerParamConfig::Instance().doUPCIteration) {
  //       if (!vtxSpan.empty()) {
  //         if (vtxSpan[0].isFlagSet(Vertex::UPCMode) == 1) { // at least one vertex in this ROF and it is from second vertex iteration
  //           LOGP(debug, "ROF {} rejected as vertices are from the UPC iteration", iRof);
  //           // processUPCMask[1][iRof] = 1;
  //           // cutUPCVertex++;
  //           vtxROF.setFlag(o2::itsmft::ROFRecord::VtxUPCMode);
  //         } else { // in all cases except if as standard mode vertex was found, the ROF was processed with UPC settings
  //           vtxROF.setFlag(o2::itsmft::ROFRecord::VtxStdMode);
  //         }
  //       } else {
  //         vtxROF.setFlag(o2::itsmft::ROFRecord::VtxUPCMode);
  //       }
  //     } else {
  //       vtxROF.setFlag(o2::itsmft::ROFRecord::VtxStdMode);
  //     }
  //     vtxROF.setNEntries(vtxSpan.size());
  //     bool selROF = vtxSpan.empty();
  //     for (int iV{0}, iVC{0}; iV < vtxSpan.size(); ++iV) {
  //       const auto& v = vtxSpan[iV];
  //       if (multEstConf.isVtxMultCutRequested() && !multEstConf.isPassingVtxMultCut(v.getNContributors())) {
  //         iVC += v.getNContributors();
  //         continue; // skip vertex of unwanted multiplicity
  //       }
  //       selROF = true;
  //       vertices.push_back(v);
  //       if (mIsMC && !VertexerParamConfig::Instance().useTruthSeeding) {
  //         allVerticesLabels.push_back(vMCRecInfo[iV].first);
  //         allVerticesPurities.push_back(vMCRecInfo[iV].second);
  //       }
  //       iVC += v.getNContributors();
  //     }
  //     // FIXME
  //     // if (processingMask[iRof] && !selROF) { // passed selection in clusters and not in vertex multiplicity
  //     //   LOGP(info, "ROF {} rejected by the vertex multiplicity selection [{},{}]", iRof, multEstConf.cutMultVtxLow, multEstConf.cutMultVtxHigh);
  //     //   processingMask[iRof] = selROF;
  //     //   cutVertexMult++;
  //     // }
  //   } else { // cosmics
  //     vtxVecLoc.emplace_back();
  //     vtxVecLoc.back().setNContributors(1);
  //     vtxROF.setNEntries(vtxVecLoc.size());
  //     for (auto& v : vtxVecLoc) {
  //       vertices.push_back(v);
  //     }
  //     mTimeFrame->addPrimaryVertices(vtxVecLoc, 0);
  //   }
  // }
  // if (mRunVertexer) {
  //   LOG(info) << fmt::format(" - Vertex seeding total elapsed time: {} ms for {} ({} + {}) vertices found in {}/{} ROFs",
  //                            vertexerElapsedTime,
  //                            mTimeFrame->getPrimaryVerticesNum(),
  //                            mTimeFrame->getTotVertIteration()[0],
  //                            o2::its::VertexerParamConfig::Instance().nIterations > 1 ? mTimeFrame->getTotVertIteration()[1] : 0,
  //                            rofsinput[1].size() - mTimeFrame->getNoVertexROF(),
  //                            rofsinput[1].size());
  //   // LOG(info) << fmt::format(" - FastMultEst: rejected {}/{} ROFs: random/mult.sel:{} (seed {}), vtx.sel:{}", cutRandomMult + cutVertexMult, trackROFspan.size(), cutRandomMult, multEst.lastRandomSeed, cutVertexMult);
  // }
  // if (mOverrideBeamEstimation) {
  //   LOG(info) << fmt::format(" - Beam position set to: {}, {} from meanvertex object", mTimeFrame->getBeamX(), mTimeFrame->getBeamY());
  // } else {
  //   LOG(info) << fmt::format(" - Beam position computed for the TF: {}, {}", mTimeFrame->getBeamX(), mTimeFrame->getBeamY());
  // }
  // mTimeFrame->setMultiplicityCutMask(processingMask);
  // mTimeFrame->setROFMask(processUPCMask);
  // Run CA tracker
  if (mMode == o2::its::TrackingMode::Async && o2::its::TrackerParamConfig::Instance().fataliseUponFailure) {
    mTracker->clustersToTracks(logger, fatalLogger);
  } else {
    mTracker->clustersToTracks(logger, errorLogger);
  }

  if (mTimeFrame->hasBogusClusters()) {
    LOG(warning) << fmt::format(" - The processed timeframe had {} clusters with wild z coordinates, check the dictionaries", mTimeFrame->hasBogusClusters());
  }

  // vertices are already sorted in time
  if (size_t totVertices{mTimeFrame->getPrimaryVerticesNum()}; totVertices) {
    allVerticesLabels.reserve(totVertices);
    if (mIsMC) {
      allVerticesLabels.reserve(totVertices);
      allVerticesPurities.reserve(totVertices);
    }
    for (size_t iVtx{0}; iVtx < totVertices; ++iVtx) {
      vertices.emplace_back(mTimeFrame->getPrimaryVertex(iVtx));
      if (mIsMC) {
        const auto& lbl = mTimeFrame->getPrimaryVerticesLabels()[iVtx];
        allVerticesLabels.emplace_back(lbl.first);
        allVerticesPurities.emplace_back(lbl.second);
      }
    }
  }

  // tracks are already sorted in time
  if (size_t totTracks{mTimeFrame->getNumberOfTracks()}; totTracks) {
    allTracks.reserve(totTracks);
    allClusIdx.reserve(mTimeFrame->getNumberOfUsedClusters());
    if (mIsMC) {
      allTrackLabels.reserve(totTracks);
    }
    const auto& lbls = mTimeFrame->getTracksLabel();
    auto& tracks = mTimeFrame->getTracks();
    for (size_t iTrk{0}; iTrk < totTracks; ++iTrk) {
      auto& trc{tracks[iTrk]};
      trc.setFirstClusterEntry(allClusIdx.size()); // before adding tracks, create final cluster indices
      int ncl = trc.getNumberOfClusters(), nclf = 0;
      for (int ic = TrackITSExt::MaxClusters; ic--;) { // track internally keeps in->out cluster indices, but we want to store the references as out->in!!!
        auto clid = trc.getClusterIndex(ic);
        if (clid >= 0) {
          trc.setClusterSize(ic, mTimeFrame->getClusterSize(ic, clid));
          allClusIdx.emplace_back(clid);
          nclf++;
        }
      }
      assert(ncl == nclf);
      allTracks.emplace_back(trc);
      if (mIsMC) {
        allTrackLabels.emplace_back(mTimeFrame->getTracksLabel()[iTrk]);
      }
    }
  }

  LOGP(info, "ITSTracker pushed {} tracks and {} vertices", allTracks.size(), vertices.size());
  if (mIsMC) {
    LOGP(info, "ITSTracker pushed {} track labels", allTrackLabels.size());
    LOGP(info, "ITSTracker pushed {} vertex labels and purities", allVerticesLabels.size());
  }
  mTimeFrame->wipe();
}

void ITSTrackingInterface::updateTimeDependentParams(framework::ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  static bool initOnceDone = false;
  if (mOverrideBeamEstimation) {
    pc.inputs().get<o2::dataformats::MeanVertexObject*>("meanvtx");
  }
  if (!initOnceDone) { // this params need to be queried only once
    initOnceDone = true;
    pc.inputs().get<o2::itsmft::TopologyDictionary*>("itscldict"); // just to trigger the finaliseCCDB
    pc.inputs().get<o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>*>("itsalppar");
    if (pc.inputs().getPos("itsTGeo") >= 0) {
      pc.inputs().get<o2::its::GeometryTGeo*>("itsTGeo");
    }
    GeometryTGeo* geom = GeometryTGeo::Instance();
    geom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::T2GRot, o2::math_utils::TransformType::T2G));
    initialise();

    if (pc.services().get<const o2::framework::DeviceSpec>().inputTimesliceId == 0) { // print settings only for the 1st pipeling
      o2::its::TrackerParamConfig::Instance().printKeyValues(true, true);
      for (const auto& par : getTracker()->getParameters()) {
        LOGP(info, "{}", par.asString());
      }
    }
  }
}

void ITSTrackingInterface::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  if (matcher == ConcreteDataMatcher("ITS", "CLUSDICT", 0)) {
    LOG(info) << "cluster dictionary updated";
    setClusterDictionary((const o2::itsmft::TopologyDictionary*)obj);
    return;
  }
  // Note: strictly speaking, for Configurable params we don't need finaliseCCDB check, the singletons are updated at the CCDB fetcher level
  if (matcher == ConcreteDataMatcher("ITS", "ALPIDEPARAM", 0)) {
    LOG(info) << "Alpide param updated";
    const auto& par = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>::Instance();
    par.printKeyValues(true, true);
    return;
  }
  if (matcher == ConcreteDataMatcher("GLO", "MEANVERTEX", 0)) {
    LOGP(info, "Mean vertex acquired");
    setMeanVertex((const o2::dataformats::MeanVertexObject*)obj);
    return;
  }
  if (matcher == ConcreteDataMatcher("ITS", "GEOMTGEO", 0)) {
    LOG(info) << "ITS GeometryTGeo loaded from ccdb";
    o2::its::GeometryTGeo::adopt((o2::its::GeometryTGeo*)obj);
    return;
  }
}

void ITSTrackingInterface::printSummary() const
{
  mTracker->printSummary();
}

void ITSTrackingInterface::setTraitsFromProvider(TrackerTraitsN* trackerTraits,
                                                 TimeFrameN* frame)
{
  mTracker = std::make_unique<TrackerN>(trackerTraits);
  mTimeFrame = frame;
  mTracker->adoptTimeFrame(*mTimeFrame);

  // set common memory resource
  if (!mMemoryPool) {
    mMemoryPool = std::make_shared<BoundedMemoryResource>();
  }
  trackerTraits->setMemoryPool(mMemoryPool);
  mTimeFrame->setMemoryPool(mMemoryPool);
  mTracker->setMemoryPool(mMemoryPool);
}

void ITSTrackingInterface::loadROF(gsl::span<const itsmft::ROFRecord>& trackROFspan,
                                   gsl::span<const itsmft::CompClusterExt> clusters,
                                   gsl::span<const unsigned char>::iterator& pattIt,
                                   int layer,
                                   const dataformats::MCTruthContainer<MCCompLabel>* mcLabels)
{
  mTimeFrame->loadROFrameData(trackROFspan, clusters, pattIt, mDict, layer, mcLabels);
}
