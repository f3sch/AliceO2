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

#include <utility>
#include <vector>

#include "TDatabasePDG.h"
#include "TStopwatch.h"

#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "DetectorsBase/Propagator.h"
#include "SimulationDataFormat/MCEventLabel.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "GlobalTrackingStudy/TPCV0Study.h"
#include "ReconstructionDataFormats/V0.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "Steer/MCKinematicsReader.h"
#include "Framework/Task.h"

namespace o2::trackstudy
{
using namespace o2::framework;
using DetID = o2::detectors::DetID;
using GTrackID = o2::dataformats::GlobalTrackID;
using DataRequest = o2::globaltracking::DataRequest;

using timeEst = o2::dataformats::TimeStampWithError<float, float>;

class TPCV0StudySpec final : public Task
{
 public:
  TPCV0StudySpec(const TPCV0StudySpec&) = delete;
  TPCV0StudySpec(TPCV0StudySpec&&) = delete;
  TPCV0StudySpec& operator=(const TPCV0StudySpec&) = delete;
  TPCV0StudySpec& operator=(TPCV0StudySpec&&) = delete;
  TPCV0StudySpec(std::shared_ptr<DataRequest> dr, std::shared_ptr<o2::base::GRPGeomRequest> gr, bool useMC)
    : mDataRequest(std::move(dr)), mGGCCDBRequest(std::move(gr)), mUseMC(useMC) {}
  ~TPCV0StudySpec() final = default;
  void init(InitContext& ic) final;
  void run(ProcessingContext& pc) final;
  void endOfStream(EndOfStreamContext& ec) final;
  void finaliseCCDB(ConcreteDataMatcher& matcher, void* obj) final;
  void process(o2::globaltracking::RecoContainer& recoData);

 private:
  static void updateTimeDependentParams(ProcessingContext& pc);
  TStopwatch mTimer;
  std::shared_ptr<DataRequest> mDataRequest{};
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest{};
  bool mUseMC{false};
  std::unique_ptr<o2::utils::TreeStreamRedirector> mTree;
  o2::steer::MCKinematicsReader mcReader;

  gsl::span<const o2::tpc::TrackTPC> mTPCTracksArray;
  gsl::span<const o2::MCCompLabel> mTPCTrkLabels;
  gsl::span<const o2::dataformats::PrimaryVertex> mPrimVertices;
  gsl::span<const o2::dataformats::VtxTrackRef> mPrimVer2TRefs;
  gsl::span<const o2::MCEventLabel> mMCPrimVertices;
  gsl::span<const o2::dataformats::V0> mV0s;
  gsl::span<const o2::dataformats::V0Index> mV0sIdx;
};

void TPCV0StudySpec::init(InitContext& /*ic*/)
{
  mTimer.Stop();
  mTimer.Reset();
  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);
  mTree = std::make_unique<o2::utils::TreeStreamRedirector>("tpc-v0-study.root", "recreate");
}

void TPCV0StudySpec::run(ProcessingContext& pc)
{
  mTimer.Start(false);
  o2::globaltracking::RecoContainer recoData;
  recoData.collectData(pc, *mDataRequest); // select tracks of needed type, with minimal cuts, the real selected will be done in the vertexer
  updateTimeDependentParams(pc);           // Make sure this is called after recoData.collectData, which may load some conditions
  process(recoData);
  mTimer.Stop();
}

void TPCV0StudySpec::updateTimeDependentParams(ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  static bool initOnceDone = false;
  if (!initOnceDone) { // this params need to be queried only once
    initOnceDone = true;
    // none at the moment
  }
}

void TPCV0StudySpec::process(o2::globaltracking::RecoContainer& recoData)
{
  mTPCTracksArray = recoData.getTPCTracks();
  mPrimVertices = recoData.getPrimaryVertices();
  mPrimVer2TRefs = recoData.getPrimaryVertexMatchedTrackRefs();
  mMCPrimVertices = recoData.getPrimaryVertexMCLabels();
  mV0s = recoData.getV0s();
  mV0sIdx = recoData.getV0sIdx();
  if (!mV0sIdx.empty() && (!mV0s.empty() && mV0sIdx.size() != mV0s.size())) {
    LOGP(fatal, "Mismatch between input SVertices indices and kinematics (not requested?): V0: {}/{}", mV0sIdx.size(), mV0s.size());
  }
  LOGP(info, "Found {} reconstructed V0", mV0sIdx.size());

  std::vector<o2::InteractionTimeRecord> intRecs;
  if (mUseMC) { // extract MC tracks
    const o2::steer::DigitizationContext* digCont = nullptr;
    if (!mcReader.initFromDigitContext("collisioncontext.root")) {
      LOGP(fatal, "Initialization of MCKinematicsReader failed!");
    }
    digCont = mcReader.getDigitizationContext();
    intRecs = digCont->getEventRecords();
    mTPCTrkLabels = recoData.getTPCTracksMCLabels();
  }

  for (const auto& v0Idx : mV0sIdx) {
    // auto v0 = mV0s[iv];
    // auto recoPosTrk = v0.getProng(0);
    // auto recoEleTrk = v0.getProng(1);
    const auto& recoPosTrk = recoData.getTPCTrack(v0Idx.getProngID(0));
    const auto& recoEleTrk = recoData.getTPCTrack(v0Idx.getProngID(1));
    const auto& recoPVtx = recoData.getPrimaryVertex(v0Idx.getVertexID());
    if (mUseMC) {
      auto mcPosLab = mTPCTrkLabels[v0Idx.getProngID(0)];
      auto mcEleLab = mTPCTrkLabels[v0Idx.getProngID(1)];
      if (!mcPosLab.isValid() || !mcEleLab.isValid() || mcPosLab.getEventID() != mcEleLab.getEventID()) {
        continue;
      }
      auto mcPosTrk = mcReader.getTrack(mcPosLab);
      auto mcEleTrk = mcReader.getTrack(mcEleLab);
      if ((mcPosTrk == nullptr) || (mcEleTrk == nullptr)) {
        continue;
      }
      if (auto mcPosMotherId = mcPosTrk->getMotherTrackId(), mcEleMotherId = mcEleTrk->getMotherTrackId();
          mcPosMotherId != mcEleMotherId || mcPosMotherId == -1 || mcEleMotherId == -1) {
        continue;
      }
      auto mcMother = mcReader.getTrack(mcPosLab.getEventID(), mcPosTrk->getMotherTrackId());
      // Propagate mc tracks to same x as reco tracks
      std::array<float, 3> mcPosXYZ{static_cast<float>(mcPosTrk->GetStartVertexCoordinatesX()), static_cast<float>(mcPosTrk->GetStartVertexCoordinatesY()), static_cast<float>(mcPosTrk->GetStartVertexCoordinatesZ())};
      std::array<float, 3> mcPosMomXYZ{static_cast<float>(mcPosTrk->GetStartVertexMomentumX()), static_cast<float>(mcPosTrk->GetStartVertexMomentumY()), static_cast<float>(mcPosTrk->GetStartVertexMomentumZ())};
      TParticlePDG* mcPosPDG = TDatabasePDG::Instance()->GetParticle(mcPosTrk->GetPdgCode());
      std::array<float, 3> mcEleXYZ{static_cast<float>(mcEleTrk->GetStartVertexCoordinatesX()), static_cast<float>(mcEleTrk->GetStartVertexCoordinatesY()), static_cast<float>(mcEleTrk->GetStartVertexCoordinatesZ())};
      std::array<float, 3> mcEleMomXYZ{static_cast<float>(mcEleTrk->GetStartVertexMomentumX()), static_cast<float>(mcEleTrk->GetStartVertexMomentumY()), static_cast<float>(mcEleTrk->GetStartVertexMomentumZ())};
      TParticlePDG* mcElePDG = TDatabasePDG::Instance()->GetParticle(mcEleTrk->GetPdgCode());
      if ((mcPosPDG == nullptr) || (mcElePDG == nullptr)) {
        continue;
      }
      o2::track::TrackPar mcPosTrkProp(mcPosXYZ, mcPosMomXYZ, TMath::Nint(mcPosPDG->Charge() / 3), false);
      o2::track::TrackPar mcEleTrkProp(mcEleXYZ, mcEleMomXYZ, TMath::Nint(mcElePDG->Charge() / 3), false);
      if (auto prop = o2::base::Propagator::Instance(); !mcPosTrkProp.rotate(recoPosTrk.getAlpha()) ||
                                                        !mcEleTrkProp.rotate(recoEleTrk.getAlpha()) ||
                                                        !prop->PropagateToXBxByBz(mcPosTrkProp, recoPosTrk.getX()) ||
                                                        !prop->PropagateToXBxByBz(mcEleTrkProp, recoEleTrk.getX())) {
        continue;
      }

      (*mTree) << "mcTPCTracks"
               << "mcPosTrk" << mcPosTrk
               << "mcEleTrk" << mcEleTrk
               << "mcPosTrkProp" << mcPosTrkProp
               << "mcEleTrkProp" << mcEleTrkProp
               << "mcMother" << mcMother
               << "\n";
    }
    LOGP(info, "WRITING STUFF");
    (*mTree) << "tpcTracks"
             << "recoPosTrk" << recoPosTrk
             << "recoEleTrk" << recoEleTrk
             << "recoPVtx" << recoPVtx
             << "\n";
  }
}

void TPCV0StudySpec::endOfStream(EndOfStreamContext& /*ec*/)
{
  mTree.reset();
  LOGF(info, "TPC V0 Study total timing: Cpu: %.3e Real: %.3e s in %d slots",
       mTimer.CpuTime(), mTimer.RealTime(), mTimer.Counter() - 1);
}

void TPCV0StudySpec::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
}

DataProcessorSpec getTPCV0StudySpec(GTrackID::mask_t srcTracks, bool useMC)
{
  auto dataRequest = std::make_shared<DataRequest>();

  dataRequest->requestTracks(srcTracks, useMC);
  dataRequest->requestPrimaryVertertices(useMC);
  dataRequest->requestSecondaryVertices(useMC);
  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
                                                              true,                              // GRPECS=true
                                                              false,                             // GRPLHCIF
                                                              true,                              // GRPMagField
                                                              true,                              // askMatLUT
                                                              o2::base::GRPGeomRequest::Aligned, // geometry
                                                              dataRequest->inputs,
                                                              true);

  return DataProcessorSpec{
    "tpc-v0-study",
    dataRequest->inputs,
    {},
    AlgorithmSpec{adaptFromTask<TPCV0StudySpec>(dataRequest, ggRequest, useMC)},
    {}};
}

} // namespace o2::trackstudy
