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

#ifndef O2_TRACKING_STUDY_CONFIG_H
#define O2_TRACKING_STUDY_CONFIG_H
#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include <vector>

namespace o2::trackstudy
{
struct TrackMCStudyConfig : o2::conf::ConfigurableParamHelper<TrackMCStudyConfig> {
  float minPt = 0.05;
  float maxTgl = 1.5;
  float minPtMC = 0.05;
  float maxTglMC = 1.5;
  float maxRMC = 33.;
  float maxPosTglMC = 2.;
  float maxPVZOffset = 15.;
  float decayMotherMaxT = 1.0f; // max TOF in ns for mother particles to study
  bool requireITSorTPCTrackRefs = true;
  bool requireTopBottomRefs = false;
  bool storeTPCTrackRefs = false;
  bool storeITSInfo = true;
  int minTPCRefsToExtractClRes = 2;
  int nOccBinsDrift = 10; // number of bins for TPC max drift time, where we integrate the occupancies
  int nTBPerOccBin = 48;  // number of TB per occ bin
  float rejectClustersResStat = 0.1;
  float maxTPCRefExtrap = 2;                         // max dX to extrapolate the track ref when extrapolating track true posions
  int minITSClForITSoutput = 7;                      // create special ITS otput only for long enough tracks
  std::vector<int> decayPDG = {310, 3122, 411, 421}; // decays to study, matched on |PDG|
  std::vector<int> selectMCPDG = {};                 // if non-empty, only MC tracks with exactly these PDG codes are collected
  bool checkSVertexerCuts = true;                    // replay the SVertexer V0 selection on the prongs of decays
  bool storeRecSV = true;                            // store all reconstructed V0s with the MC origin of their prongs
  float recSVSamplingFrac = 1.f;                     // fraction of combinatorial reconstructed V0s to store, true decays are always kept

  /// Specific output to study gamma conversions
  bool storeGammaConversions = false;         // store the dedicated gammaConv tree
  bool gammaRequirePrimaryOrGenerator = true; // require a physical-primary or generator photon
  float gammaMinPt = 0.01f;                   // minimum truth photon pT in GeV/c
  float gammaMaxPt = 9999999.f;               // maximum truth photon pT in GeV/c
  float gammaMaxEta = 0.9f;                   // maximum absolute truth photon eta
  float gammaMinR = 2.f;                      // minimum truth conversion radius in cm
  float gammaMaxR = 90.f;                     // maximum truth conversion radius in cm
  bool gammaApplyRZLineCut = true;            // apply the conversion R-Z fiducial line
  float gammaRZMargin = 7.f;                  // R-Z line offset in cm
  bool gammaRequireTPC = true;                // reject ITS-only daughter representations
  bool gammaRequireCleanTrackLabel = true;    // require the AOD-equivalent MC track-label mask to be zero

  /// Cuts applied to TPC-only seeds by SVertexer::processTPCTrack under mTPCTrackPhotonTune.
  /// Requires storeGammaConversions, which is what defines the signal sample.
  bool storeTPCPhotonTune = false;      // store the dedicated tpcTune tree
  float tpcTuneBkgSamplingFrac = 0.02f; // fraction of non-signal TPC tracks to store, signal is always kept

  O2ParamDef(TrackMCStudyConfig, "trmcconf");
};
} // namespace o2::trackstudy

#endif
