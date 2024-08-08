// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "ITS3Align/MisalignmentHits.h"
#include "ITS3Base/SegmentationSuperAlpide.h"
#include "ITS3Base/ITS3Params.h"
#include "Framework/Logger.h"

#include "Math/Factory.h"
#include "Math/UnaryOperators.h"
#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "TString.h"

#include <memory>
#include <string>
#include <cstring>

namespace o2::its3::align
{

void MisAlignmentHits::init()
{
  if (o2::its3::ITS3Params::Instance().misalignmentHitsUseProp) {
    mMethod = PropMethod::Propagator;
  } else {
    mMethod = PropMethod::Line;
  }

  mGeo = o2::its::GeometryTGeo::Instance();

  mMinimizer.reset(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
  if (mMinimizer == nullptr) {
    LOGP(fatal, "Cannot create minimizer");
  }
  mMinimizer->SetMaxFunctionCalls(10'000'000);
  mMinimizer->SetStrategy(1);
  mMinimizer->SetPrintLevel(0);

  if (mMethod == PropMethod::Propagator) {
    mMCReader = std::make_unique<o2::steer::MCKinematicsReader>("collisioncontext.root");
    LOGP(info, "Using propagator to find intersection");
  } else {
    LOGP(info, "Using local straight-line to find intersection");
    mMinimizer->SetFunction(mLine);
  }

  resetStats();

  if (auto file = o2::its3::ITS3Params::Instance().misalignmentHitsParams; file.empty()) {
    LOGP(fatal, "No parameter file specified");
  } else {
    mDeformations.init(file);
  }
}

std::optional<o2::itsmft::Hit> MisAlignmentHits::processHit(const o2::itsmft::Hit& hit)
{
  LOGP(debug, "Procesing hit on {} at {},{},{} phi={}", hit.GetDetectorID(), hit.GetPos().X(), hit.GetPos().Y(), hit.GetPos().Z(), o2::math_utils::to02Pi(hit.GetPos().phi()));
  ++mStats[Stats::kHitTotal];

  if (!constants::detID::isDetITS3(hit.GetDetectorID())) {
    ++mStats[Stats::kHitIsOB];
    return hit;
  }
  ++mStats[Stats::kHitIsIB];

  // Set the working hits
  mCurHit = hit;
  mCurWorkingHits[WorkingHit::kEntering] = WorkingHit(WorkingHit::kEntering, hit);
  mCurWorkingHits[WorkingHit::kExiting] = WorkingHit(WorkingHit::kExiting, hit);

  // Do work
  if (!deformHit(WorkingHit::kEntering) || !deformHit(WorkingHit::kExiting)) {
    ++mStats[Stats::kHitDead];
    return std::nullopt;
  }
  ++mStats[Stats::kHitAlive];

  // Set the possibly new detectorIDs with mid point approximation
  auto midPointOrig = mCurWorkingHits[WorkingHit::kEntering].mPoint + (mCurWorkingHits[WorkingHit::kExiting].mPoint - mCurWorkingHits[WorkingHit::kEntering].mPoint) * 0.5;
  auto midPointDef = mCurWorkingHits[WorkingHit::kEntering].mPointDef + (mCurWorkingHits[WorkingHit::kExiting].mPointDef - mCurWorkingHits[WorkingHit::kEntering].mPointDef) * 0.5;
  const int idDef = getDetID(midPointDef), idOrig = getDetID(midPointOrig);
  if (idDef == -1) {
    return std::nullopt;
  }

  if (idDef != idOrig) {
    ++mStats[Stats::kHitMigrated];
  } else {
    ++mStats[Stats::kHitNotMigrated];
  }

  // Get new postion
  mCurHit.SetPosStart(mCurWorkingHits[WorkingHit::kEntering].mPointDef);
  mCurHit.SetPos(mCurWorkingHits[WorkingHit::kExiting].mPointDef);
  mCurHit.SetDetectorID(idDef);

  ++mStats[Stats::kHitSuccess];
  return mCurHit;
}

bool MisAlignmentHits::deformHit(WorkingHit::HitType t)
{
  auto& wHit = mCurWorkingHits[t];

  mMinimizer->Clear(); // clear for next iteration
  if (mMethod == PropMethod::Line) {
    prepareLineMethod(t);
  } else {
    preparePropagtorMethod(t);
  }
  constexpr double minStep{1e-5};
  constexpr double phiMargin{0.3};
  constexpr double zMargin{4.0};
  mMinimizer->SetLimitedVariable(0, "t", 0.0, minStep, -2.5, 2.5);
  mMinimizer->SetLimitedVariable(1, "phiStar", wHit.mPhi, minStep,
                                 std::max(static_cast<double>(wHit.mPhiBorder1), static_cast<double>(wHit.mPhi) - phiMargin),
                                 std::min(static_cast<double>(wHit.mPhiBorder2), static_cast<double>(wHit.mPhi) + phiMargin));
  mMinimizer->SetLimitedVariable(2, "zStar", wHit.mPoint.Z(), minStep,
                                 std::max(static_cast<double>(-constants::segment::lengthSensitive / 2.f), static_cast<double>(wHit.mPoint.Z()) - zMargin),
                                 std::min(static_cast<double>(constants::segment::lengthSensitive / 2.f), static_cast<double>(wHit.mPoint.Z()) + zMargin));

  mMinimizer->Minimize(); // perform the actual minimization

  auto ss = mMinimizer->Status();
  if (ss == 1) {
    ++mStats[Stats::kMinimizerCovPos];
  } else if (ss == 2) {
    ++mStats[Stats::kMinimizerHesse];
  } else if (ss == 3) {
    ++mStats[Stats::kMinimizerEDM];
  } else if (ss == 4) {
    ++mStats[Stats::kMinimizerLimit];
  } else if (ss == 5) {
    ++mStats[Stats::kMinimizerOther];
  } else {
    ++mStats[Stats::kMinimizerConverged];
  }

  if (ss == 0 || ss == 1) { // for Minuit2 0=ok, 1=ok with pos. forced hesse
    ++mStats[Stats::kMinimizerStatusOk];
    if (mMinimizer->MinValue() < 1e-4) { // within 1 um considering the pixel pitch this good enough
      ++mStats[Stats::kMinimizerValueOk];
    } else {
      ++mStats[Stats::kMinimizerValueBad];
      return false;
    }
  } else {
    ++mStats[Stats::kMinimizerStatusBad];
    return false;
  }

  // Valid solution found; calculate new position on ideal geo
  wHit.recalculateIdeal(static_cast<float>(mMinimizer->X()[1]), static_cast<float>(mMinimizer->X()[2]));

  return true;
}

int MisAlignmentHits::getDetID(const o2::math_utils::Point3D<float>& point)
{
  gGeoManager->PushPath();
  auto id = getDetIDFromCords(point);
  gGeoManager->PopPath();
  return id;
}

int MisAlignmentHits::getDetIDFromCords(const o2::math_utils::Point3D<float>& point)
{
  // retrive if any the node which constains the point
  const auto node = gGeoManager->FindNode(point.X(), point.Y(), point.Z());
  if (node == nullptr) {
    ++mStats[Stats::kFindNodeFailed];
    return -1;
  }
  ++mStats[Stats::kFindNodeSuccess];

  // check if this node is a sensitive volume
  const std::string path = gGeoManager->GetPath();
  if (path.find(o2::its::GeometryTGeo::getITS3SensorPattern()) == std::string::npos) {
    ++mStats[Stats::kProjNonSensitive];
    return -1;
  }
  ++mStats[Stats::kProjSensitive];

  return getDetIDFromPath(path);
}

int MisAlignmentHits::getDetIDFromPath(const std::string& path) const
{
  static const std::regex pattern{R"(/cave_1/barrel_1/ITSV_2/ITSUWrapVol0_1/ITS3Layer(\d+)_(\d+)/ITS3CarbonForm(\d+)_(\d+)/ITS3Chip(\d+)_(\d+)/ITS3Segment(\d+)_(\d+)/ITS3RSU(\d+)_(\d+)/ITS3Tile(\d+)_(\d+)/ITS3PixelArray(\d+)_(\d+))"};
  if (std::smatch matches; std::regex_search(path, matches, pattern)) {
    if (matches.size() == 15) {
      int iLayer = std::stoi(matches[1]);
      int iCarbonForm = std::stoi(matches[4]);
      int iSegment = std::stoi(matches[8]);
      int iRSU = std::stoi(matches[10]);
      int iTile = std::stoi(matches[12]);
      return mGeo->getChipIndex(iLayer, iCarbonForm, 0, iSegment, iRSU, iTile);
    } else {
      LOGP(fatal, "Path did not contain expected number of matches ({})!", matches.size());
    }
  } else {
    LOGP(fatal, "Path was not matched ({})!", path);
  }
  __builtin_unreachable();
}

void MisAlignmentHits::printStats() const
{
  LOGP(info, "Processed {} Hits (IB:{}; OB:{}):", mStats[Stats::kHitTotal], mStats[Stats::kHitIsIB], mStats[Stats::kHitIsOB]);
  LOGP(info, "  - Minimizer Status: {} ok {} bad (2x)", mStats[Stats::kMinimizerStatusOk], mStats[Stats::kMinimizerStatusBad]);
  LOGP(info, "  - Minimizer Value: {} ok {} bad (2x)", mStats[Stats::kMinimizerValueOk], mStats[Stats::kMinimizerValueBad]);
  LOGP(info, "  - Minimizer Detailed: {} Converged {} pos. forced Hesse (2x)", mStats[Stats::kMinimizerConverged], mStats[Stats::kMinimizerHesse]);
  LOGP(info, "  - Minimizer Detailed: {} EDM {} call limit {} other (2x)", mStats[Stats::kMinimizerEDM], mStats[Stats::kMinimizerLimit], mStats[Stats::kMinimizerOther]);
  LOGP(info, "  - FindNode: {} ok {} failed", mStats[Stats::kFindNodeSuccess], mStats[Stats::kFindNodeFailed]);
  LOGP(info, "  - IsSensitve: {} yes {} no", mStats[Stats::kProjSensitive], mStats[Stats::kProjNonSensitive]);
  LOGP(info, "  - IsAlive: {} yes {} no", mStats[Stats::kHitAlive], mStats[Stats::kHitDead]);
  LOGP(info, "  - HasMigrated: {} yes {} no", mStats[Stats::kHitMigrated], mStats[Stats::kHitNotMigrated]);
  LOGP(info, "  --> Good Hits {}", mStats[Stats::kHitSuccess]);
}

void MisAlignmentHits::prepareLineMethod(WorkingHit::HitType from)
{
  // Set the starint point and radius
  // always start from the entering hit that way t is always pos. defined
  mLine.mStart = mCurWorkingHits[WorkingHit::kEntering].mPoint;
  mLine.mRadius = mCurWorkingHits[from].mRadius;
  mLine.mSensorID = mCurWorkingHits[from].mSensorID;
  // Calculate the direction vector
  mLine.mD[0] = mCurWorkingHits[WorkingHit::kExiting].mPoint.X() - mCurWorkingHits[WorkingHit::kEntering].mPoint.X();
  mLine.mD[1] = mCurWorkingHits[WorkingHit::kExiting].mPoint.Y() - mCurWorkingHits[WorkingHit::kEntering].mPoint.Y();
  mLine.mD[2] = mCurWorkingHits[WorkingHit::kExiting].mPoint.Z() - mCurWorkingHits[WorkingHit::kEntering].mPoint.Z();
}

} // namespace o2::its3::align
