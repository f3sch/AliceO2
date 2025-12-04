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

#ifndef ALICEO2_ITSDPLTRACKINGPARAM_H_
#define ALICEO2_ITSDPLTRACKINGPARAM_H_

#include <limits>
#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"

namespace o2::its
{

struct TrackerParamConfig : public o2::conf::ConfigurableParamHelper<TrackerParamConfig> {
  // Use TGeo for mat. budget
  static const int MaxIter = 4;
  static const int MinTrackLength = 4;
  static const int MaxTrackLength = 7;
  bool useMatCorrTGeo = false;                                              // use full geometry to correct for material budget accounting in the fits. Default is to use the material budget LUT.
  bool useFastMaterial = false;                                             // use faster material approximation for material budget accounting in the fits.
  int deltaROF[MaxTrackLength] = {0};                                       // configure the width of the window in ROFs to be considered for the tracking.
  int minTrackLgtIter[MaxIter] = {};                                        // minimum track length at each iteration, used only if >0, otherwise use code defaults
  uint8_t startLayerMask[MaxIter] = {};                                     // mask of start layer for this iteration (if >0)
  float minPtIterLgt[MaxIter * (MaxTrackLength - MinTrackLength + 1)] = {}; // min.pT for given track length at this iteration, used only if >0, otherwise use code defaults
  float sysErrY2[MaxTrackLength] = {0};                                     // systematic error^2 in Y per layer
  float sysErrZ2[MaxTrackLength] = {0};                                     // systematic error^2 in Z per layer
  int nTimeSlices{1};
  float maxChi2ClusterAttachment = -1.f;
  float maxChi2NDF = -1.f;
  float nSigmaCut = -1.f;
  float minPt = -1.f;
  int LUTbinsPhi = -1;
  int LUTbinsZ = -1;
  float diamondPos[3] = {0.f, 0.f, 0.f};                         // override the position of the vertex
  float diamondCov[6] = {25.e-6f, 0.f, 0.f, 25.e-6f, 0.f, 36.f}; // cov
  bool useDiamond = false;                                       // enable overriding the vertex position
  int useTrackFollower = -1;                                     // bit 0: allow mixing implies bits 1&2; bit 1: topwards; bit2: downwards; => 0 off
  int findShortTracks = -1;
  int nROFsPerIterations = 0;              // size of the slice of ROFs to be processed at a time, preferably integer divisors of nROFs per TF, to balance the iterations.
  int nOrbitsPerIterations = 0;            // not implemented: size of the slice of ROFs to be processed at a time, computed using the number of ROFs per orbit.
  bool perPrimaryVertexProcessing = false; // perform the full tracking considering the vertex hypotheses one at the time.
  bool saveTimeBenchmarks = false;         // dump metrics on file
  bool overrideBeamEstimation = false;     // use beam position from meanVertex CCDB object
  int trackingMode = -1;                   // -1: unset, 0=sync, 1=async, 2=cosmics used by gpuwf only
  bool doUPCIteration = false;             // Perform an additional iteration for UPC events on tagged vertices. You want to combine this config with VertexerParamConfig.nIterations=2
  int nIterations = MaxIter;               // overwrite the number of iterations
  int reseedIfShorter = 6;                 // for the final refit reseed the track with circle if they are shorter than this value
  bool shiftRefToCluster{true};            // TrackFit: after update shift the linearization reference to cluster

  /// seeding
  float seedingMeanVertexExtraErr2{0.f};   // additional error imposed on mean vertex cov.
  float seedingDCATolerance{0.8f};         // maximum allowed DCA to meanvertex for track to enter pool
  float seedingDCAMaxPull{3.f};            // maximum allowed initial pull on DCA to meanvertex
  float seedingMaxChi2Iter{3.f};           // maximum chi2 change required to end iteration
  float seedingTukeyStartIter{5.f};        // start value for tukey scaling for iteration
  float seedingMinWghTrk{0.01f};           // minimum weight a track has to have to be accounted
  int seedingMaxFitIter{3};                // maximum iterations for fit
  int seedingMinTracksIter{20};            // minimum tracks needed for one iteration
  int seedingDBScanMinPt{50};              // DBSCAN: minimum number of cluster points
  float seedingDBScanEpsZ{0.02f};          // DBSCAN: maximum epsilon z
  float seedingDBScanEpsT{10.f};           // DBSCAN: maximum epsilon t (BC)
  float seedingVertexExtraErr2[6] = {0.f}; // impose additional errors on seeding vertices
  bool seedingUseMCTruth{false};           // skip seeding and impose MC event information

  bool createArtefactLabels{false}; // create on-the-fly labels for the artefacts

  int nThreads = 1;
  bool printMemory = false;
  size_t maxMemory = std::numeric_limits<size_t>::max();
  bool dropTFUponFailure = false;
  bool fataliseUponFailure = true;       // granular management of the fatalisation in async mode
  bool allowSharingFirstCluster = false; // allow first cluster sharing among tracks

  O2ParamDef(TrackerParamConfig, "ITSCATrackerParam");
};

struct ITSGpuTrackingParamConfig : public o2::conf::ConfigurableParamHelper<ITSGpuTrackingParamConfig> {
  /// Set nBlocks/nThreads to summarily override all kernel launch parameters in each iteration.
  /// Parameters must start with nBlocks/nThreads.
  static constexpr int OverrideValue{-1};
  static constexpr char const* BlocksName = "nBlocks";
  static constexpr char const* ThreadsName = "nThreads";
  int nBlocks = OverrideValue;
  int nThreads = OverrideValue;
  void maybeOverride() const;

  /// Individual kernel launch parameter for each iteration
  int nBlocksLayerTracklets = 60;
  int nThreadsLayerTracklets = 256;

  int nBlocksLayerCells = 60;
  int nThreadsLayerCells = 256;

  int nBlocksFindNeighbours = 60;
  int nThreadsFindNeighbours = 256;

  int nBlocksProcessNeighbours = 60;
  int nThreadsProcessNeighbours = 256;

  int nBlocksTracksSeeds = 60;
  int nThreadsTracksSeeds = 256;

  O2ParamDef(ITSGpuTrackingParamConfig, "ITSGpuTrackingParam");
};

} // namespace o2::its
#endif
