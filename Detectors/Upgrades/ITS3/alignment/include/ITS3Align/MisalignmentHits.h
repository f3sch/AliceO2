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

#ifndef ITS3_MISALIGNMENTHITS_H_
#define ITS3_MISALIGNMENTHITS_H_

#include "Math/IFunction.h"
#include "Math/Minimizer.h"

#include "ITS3Align/Deformations.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITSMFTSimulation/Hit.h"
#include "MathUtils/Cartesian.h"
#include "MathUtils/Utils.h"
#include "Steer/MCKinematicsReader.h"

#include <regex>
#include <memory>
#include <array>
#include <optional>
#include <tuple>

namespace o2::its3::align
{

class MisAlignmentHits
{
 public:
  enum class PropMethod {
    Propagator,
    Line,
  };

  void init();

  std::optional<o2::itsmft::Hit> processHit(const o2::itsmft::Hit& hit);

  void resetStats() { mStats.fill(0ull); }
  void printStats() const;

 private:
  Deformations mDeformations;
  std::unique_ptr<ROOT::Math::Minimizer> mMinimizer;
  PropMethod mMethod{PropMethod::Line};
  o2::its::GeometryTGeo* mGeo{nullptr};
  std::unique_ptr<o2::steer::MCKinematicsReader> mMCReader;

  int getDetID(const o2::math_utils::Point3D<float>& point);
  int getDetIDFromCords(const o2::math_utils::Point3D<float>& point);
  int getDetIDFromPath(const std::string& path) const;

  // We treat each hit as two separate hits', one for the entering and one for the exiting hit
  struct WorkingHit {
    enum HitType : uint8_t {
      kEntering = 0,
      kExiting,
      kTypes,
    };

    WorkingHit() = default;

    WorkingHit(HitType t, const o2::itsmft::Hit& hit) : mType(t),
                                                        mDetID(hit.GetDetectorID()),
                                                        mLayerID(constants::detID::getDetID2Layer(mDetID)),
                                                        mSensorID(constants::detID::getSensorID(mDetID))
    {
      if (mType == kEntering) {
        mRadius = constants::radiiInner[mLayerID];
        mPoint = hit.GetPosStart();
      } else {
        mRadius = constants::radiiOuter[mLayerID];
        mPoint = hit.GetPos();
      }

      // Pre-calculate the normalized u,v coordinates
      const bool isTop = mSensorID % 2 == 0;
      mPhi = o2::math_utils::to02Pi(std::atan2(mPoint.Y(), mPoint.X()));
      mPhiBorder1 = o2::math_utils::to02Pi(((isTop) ? 0.f : 1.f) * TMath::Pi() + std::asin(constants::equatorialGap / 2.f / mRadius));
      mPhiBorder2 = o2::math_utils::to02Pi(((isTop) ? 1.f : 2.f) * TMath::Pi() - std::asin(constants::equatorialGap / 2.f / mRadius));
      mU = ((mPhi - mPhiBorder1) * 2.f) / (mPhiBorder2 - mPhiBorder1) - 1.f;
      mV = (2.f * mPoint.Z() + constants::segment::lengthSensitive) / constants::segment::lengthSensitive - 1.f;
    }

    void recalculateIdeal(float phi, float z)
    {
      mPointDef.SetX(mRadius * std::cos(phi));
      mPointDef.SetY(mRadius * std::sin(phi));
      mPointDef.SetZ(z);
    }

    HitType mType;
    short mDetID;
    int mLayerID;
    int mSensorID;
    float mRadius;
    float mPhi;
    o2::math_utils::Point3D<float> mPoint;
    o2::math_utils::Point3D<float> mPointDef;
    float mU; // u is normalized phi
    float mV; // u is normalized z

    float mPhiBorder1;
    float mPhiBorder2;
  };
  std::array<WorkingHit, WorkingHit::kTypes> mCurWorkingHits;
  o2::itsmft::Hit mCurHit;

  bool deformHit(WorkingHit::HitType t);

  auto getDeformation(unsigned int id, double u, double v) const
  {
    return mDeformations.getDeformation(id, u, v);
  }

  // Mimize function assuming a straight line
  // given in the parametric representation by y_v = t * d_x + x_s
  // assuming no offset is needed
  class StraightLine : public ROOT::Math::IBaseFunctionMultiDim
  {
   public:
    StraightLine(const MisAlignmentHits* m) : mMis(m) {}

    std::array<double, 3> mD;
    o2::math_utils::Point3D<float> mStart;
    unsigned int mSensorID;
    double mRadius;
    const MisAlignmentHits* mMis;

    unsigned int NDim() const override { return 3; }
    ROOT::Math::IBaseFunctionMultiDim* Clone() const override { return nullptr; }

   private:
    double DoEval(const double* x) const override
    {
      const double t = x[0];
      const double phi = x[1];
      const double z = x[2];

      /// Find the point along the line given current t
      double xline = mStart.X() + t * mD[0],
             yline = mStart.Y() + t * mD[1],
             zline = mStart.Z() + t * mD[2];

      // Find the point of the deformed geometry given a certain phi' and z'
      double xideal = mRadius * std::cos(phi), yideal = mRadius * std::sin(phi),
             zideal = z;
      double u, v;
      auto [dx, dy, dz] = mMis->getDeformation(mSensorID, u, v);
      double xdef = xideal + dx, ydef = yideal + dy, zdef = zideal + dz;

      // Minimize the euclidean distance of the line point and the deformed point
      return std::hypot(xline - xdef, yline - ydef, zline - zdef);
    }
  };
  StraightLine mLine{this};

  void prepareLineMethod(WorkingHit::HitType from);
  void preparePropagtorMethod(WorkingHit::HitType from) {}

  enum Stats : uint8_t {
    kHitTotal = 0,
    kHitIsOB,
    kHitIsIB,
    kHitDead,
    kHitAlive,
    kHitSuccess,
    kHitMigrated,
    kHitNotMigrated,
    kHitEntBoundary,
    kHitExtBoundary,
    kHitNoBoundary,
    kHitSameBoundary,
    kFindNodeFailed,
    kFindNodeSuccess,
    kProjSensitive,
    kProjNonSensitive,
    kDetIDOk,
    kDetIDBad,
    kMinimizerStatusOk,
    kMinimizerStatusBad,
    kMinimizerValueOk,
    kMinimizerValueBad,
    kMinimizerConverged,
    kMinimizerCovPos,
    kMinimizerHesse,
    kMinimizerEDM,
    kMinimizerLimit,
    kMinimizerOther,
    kALL,
  };
  std::array<ULong64_t, Stats::kALL> mStats;
};

} // namespace o2::its3::align

#endif
