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

/// @file   AlignableDetectorIT3.cxx
/// @brief  IT3 detector wrapper

#include "Align/AlignableDetectorIT3.h"
#include "Align/AlignableVolume.h"
#include "Align/AlignableSensorIT3.h"
#include "Align/Controller.h"
#include "Align/AlignConfig.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITS3Reconstruction/TopologyDictionary.h"
#include "ITS3Reconstruction/IOUtils.h"
#include "DataFormatsITSMFT/TrkClusRef.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include <TMath.h>
#include <cstdio>

using namespace TMath;
using namespace o2::align::utils;

namespace o2
{
namespace align
{

AlignableDetectorIT3::AlignableDetectorIT3(Controller* ctr) : AlignableDetector(DetID::IT3, ctr) {}

void AlignableDetectorIT3::defineVolumes()
{
  auto geom = o2::its::GeometryTGeo::Instance();

  AlignableVolume *volITS = nullptr, *volLr = nullptr, *volHB = nullptr, *volSt = nullptr, *volHSt = nullptr, *volMod = nullptr;
  AlignableSensorIT3* sens = nullptr;
  //
  std::unordered_map<std::string, AlignableVolume*> sym2vol;
  addVolume(volITS = new AlignableVolume(geom->composeSymNameITS(), getDetLabel(), mController));
  sym2vol[volITS->getSymName()] = volITS;
  //
  int nonSensCnt = 0;
  for (int ilr = 0; ilr < geom->getNumberOfLayers(); ilr++) {
    for (int ihb = 0; ihb < geom->getNumberOfHalfBarrels(); ihb++) {
      addVolume(volLr = new AlignableVolume(geom->composeSymNameHalfBarrel(ilr, ihb, ilr < 3), getNonSensLabel(nonSensCnt++), mController));
      sym2vol[volLr->getSymName()] = volLr;
      volLr->setParent(volITS);

      if (ilr > 2) {
        int nstavesHB = geom->getNumberOfStaves(ilr) / 2;
        for (int ist = 0; ist < nstavesHB; ist++) {
          addVolume(volSt = new AlignableVolume(geom->composeSymNameStave(ilr, ihb, ist), getNonSensLabel(nonSensCnt++), mController));
          sym2vol[volSt->getSymName()] = volSt;
          volSt->setParent(volLr);
          for (int ihst = 0; ihst < geom->getNumberOfHalfStaves(ilr); ihst++) {
            addVolume(volHSt = new AlignableVolume(geom->composeSymNameHalfStave(ilr, ihb, ist, ihst), getNonSensLabel(nonSensCnt++), mController));
            sym2vol[volHSt->getSymName()] = volHSt;
            volHSt->setParent(volSt);
            for (int imd = 0; imd < geom->getNumberOfModules(ilr); imd++) {
              addVolume(volMod = new AlignableVolume(geom->composeSymNameModule(ilr, ihb, ist, ihst, imd), getNonSensLabel(nonSensCnt++), mController));
              sym2vol[volMod->getSymName()] = volMod;
              volMod->setParent(volHSt);
            } // module
          } // halfstave
        } // stave
      }
    } // layer halfBarrel
  } // layer

  for (int ich = 0; ich < geom->getNumberOfChips(); ich++) {
    int chID = o2::base::GeometryManager::getSensID(mDetID, ich);
    addVolume(sens = new AlignableSensorIT3(o2::base::GeometryManager::getSymbolicName(mDetID, ich), chID, getSensLabel(ich), mController));
    int lay = 0, hba, sta = 0, ssta = 0, modd = 0, chip = 0;
    geom->getChipId(ich, lay, hba, sta, ssta, modd, chip);
    // We se the parent to be the carbon form
    AlignableVolume* parVol = nullptr;
    if (lay < 3) {
      parVol = sym2vol[geom->composeSymNameHalfBarrel(lay, hba, true)];
    } else {
      auto parent = modd < 0 ? geom->composeSymNameStave(lay, hba, sta) : geom->composeSymNameModule(lay, hba, sta, ssta, modd);
      parVol = sym2vol[parent];
    }
    if (!parVol) {
      throw std::runtime_error(fmt::format("did not find parent for chip {}", chID));
    }
    sens->setParent(parVol);
  }
  //
}

//____________________________________________
int AlignableDetectorIT3::processPoints(GIndex gid, int npntCut, bool inv)
{
  // Extract the points corresponding to this detector, recalibrate/realign them to the
  // level of the "starting point" for the alignment/calibration session.
  // If inv==true, the track propagates in direction of decreasing tracking X
  // (i.e. upper leg of cosmic track)
  //
  auto algTrack = mController->getAlgTrack();
  auto recoData = mController->getRecoContainer();
  const auto& algConf = AlignConfig::Instance();
  int npoints = 0;
  auto procClus = [this, inv, &npoints, &algTrack](const ClusterD& clus) {
    auto* sensor = this->getSensor(clus.getSensorID());
    auto& pnt = algTrack->addDetectorPoint();
    const auto* sysE = sensor->getAddError(); // additional syst error
    pnt.setYZErrTracking(clus.getSigmaY2() + sysE[0] * sysE[0], clus.getSigmaYZ(), clus.getSigmaZ2() + sysE[1] * sysE[1]);
    if (this->getUseErrorParam()) { // errors will be calculated just before using the point in the fit, using track info
      pnt.setNeedUpdateFromTrack();
    }
    pnt.setXYZTracking(clus.getX(), clus.getY(), clus.getZ());
    pnt.setSensor(sensor);
    pnt.setAlphaSens(sensor->getAlpTracking());
    pnt.setXSens(sensor->getXTracking());
    pnt.setDetID(this->mDetID);
    pnt.setSID(sensor->getSID());
    pnt.setContainsMeasurement();
    pnt.setInvDir(inv);
    pnt.init();
    npoints++;
  };
  std::array<int, 7> clusIDs{};
  int nOverlaps = 0;
  if (gid.getSource() == GIndex::ITS) {
    const auto tracks = recoData->getITSTracks();
    if (tracks.empty()) {
      return -1; // source not loaded?
    }
    const auto& track = tracks[gid.getIndex()];
    if (track.getNClusters() < npntCut) {
      return -1;
    }
    const auto& clusIdx = recoData->getITSTracksClusterRefs();
    // do we want to apply some cuts?
    int clEntry = track.getFirstClusterEntry();
    int preevSensID = -1;
    bool errReported = false;
    for (int icl = track.getNumberOfClusters(); icl--;) { // clusters refs are stored from outer to inner layers, we loop in inner -> outer direction
      const auto& clus = mITSClustersArray[(clusIDs[npoints] = clusIdx[clEntry + icl])];
      if (clus.getSensorID() < preevSensID && !errReported) { // clusters are ordered from outer to inner layer, hence decreasing sensorID
        std::string errstr{};
        for (int ie = track.getNumberOfClusters(); ie--;) {
          errstr += fmt::format(" {}", mITSClustersArray[clusIdx[clEntry + ie]].getSensorID());
        }
        LOGP(error, "wrong ITS clusters order? : chips {}", errstr);
        errReported = true;
      }
      preevSensID = clus.getSensorID();
      procClus(clus);
    }
  } else { // ITSAB
    LOGP(fatal, "ITSAB tracks not supported yet");
  }
  if (npoints < npntCut) { // reset points to original start
    algTrack->suppressLastPoints(npoints);
    npoints = 0;
    return 0;
  }

  // do we need to process overlaps?
  mNPoints += npoints;
  return npoints;
}

//____________________________________________
bool AlignableDetectorIT3::prepareDetectorData()
{
  // prepare TF data for processing: convert clusters
  const auto& algConf = AlignConfig::Instance();
  auto recoData = mController->getRecoContainer();
  const auto& clusITS = recoData->getITSClusters();
  const auto& clusITSROF = recoData->getITSClustersROFRecords();
  const auto& patterns = recoData->getITSClustersPatterns();
  auto pattIt = patterns.begin();
  mITSClustersArray.clear();
  mITSClustersArray.reserve(clusITS.size());
  LOGP(info, "There are {} clusters in {} ROFs", clusITS.size(), clusITSROF.size());
  int ROFCount = 0;
  int16_t curSensID = -1;
  struct ROFChipEntry {
    int16_t rofCount = -1;
    int chipFirstEntry = -1;
  };

  for (const auto& rof : clusITSROF) {
    int maxic = rof.getFirstEntry() + rof.getNEntries();
    for (int ic = rof.getFirstEntry(); ic < maxic; ic++) {
      const auto& c = clusITS[ic];
      int16_t sensID = c.getSensorID();
      auto* sensor = getSensor(sensID);
      double sigmaY2, sigmaZ2, sigmaYZ = 0, locXYZC[3], traXYZ[3];
      auto pattItCopy = pattIt;
      auto locXYZP = o2::its3::ioutils::extractClusterData(c, pattIt, mITSDict, sigmaY2, sigmaZ2); // local ideal coordinates
      std::array<double, 3> locXYZ{locXYZP.X(), locXYZP.Y(), locXYZP.Z()};
      const auto& matAlg = sensor->getMatrixClAlg(); // local alignment matrix !!! RS FIXME
      matAlg.LocalToMaster(locXYZ.data(), locXYZC);  // aligned point in the local frame
      const auto& mat = sensor->getMatrixT2L();      // RS FIXME check if correct
      mat.MasterToLocal(locXYZC, traXYZ);
      mITSClustersArray.emplace_back(sensID, traXYZ[0], traXYZ[1], traXYZ[2], sigmaY2, sigmaZ2, sigmaYZ); // local --> tracking
    } // clusters of ROF

    ROFCount++;
  } // loop over ROFs
  return true;
}

//____________________________________________
void AlignableDetectorIT3::Print(const Option_t* opt) const
{
  AlignableDetector::Print(opt);
}

//____________________________________________
void AlignableDetectorIT3::SetAddErrorLr(int ilr, double sigY, double sigZ)
{
  // set syst. errors for specific layer
  auto geom = o2::its::GeometryTGeo::Instance();
  int chMin = geom->getFirstChipIndex(ilr), chMax = geom->getLastChipIndex(ilr);
  for (int isn = chMin; isn <= chMax; isn++) {
    getSensor(isn)->setAddError(sigY, sigZ);
  }
}

//____________________________________________
void AlignableDetectorIT3::SetSkipLr(int ilr)
{
  // exclude sensor of the layer from alignment
  auto geom = o2::its::GeometryTGeo::Instance();
  int chMin = geom->getFirstChipIndex(ilr), chMax = geom->getLastChipIndex(ilr);
  for (int isn = chMin; isn <= chMax; isn++) {
    getSensor(isn)->setSkip();
  }
}

//_________________________________________________
void AlignableDetectorIT3::setUseErrorParam(int v)
{
  // set type of points error parameterization // RS DO WE NEED THIS?
  mUseErrorParam = v;
}

} // namespace align
} // namespace o2
