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
///
/// \file Configuration.h
/// \brief
///

#ifndef TRACKINGITSU_INCLUDE_CONFIGURATION_H_
#define TRACKINGITSU_INCLUDE_CONFIGURATION_H_

#include <limits>
#include <vector>

#include "DetectorsBase/Propagator.h"
#include "CommonUtils/EnumFlags.h"
#include "ITStracking/Constants.h"

namespace o2::its
{

struct TrackingParameters {
  int CellMinimumLevel() const noexcept { return MinTrackLength - constants::ClustersPerCell + 1; }
  int NeighboursPerRoad() const noexcept { return NLayers - 3; }
  int CellsPerRoad() const noexcept { return NLayers - 2; }
  int TrackletsPerRoad() const noexcept { return NLayers - 1; }
  std::string asString() const;

  int NLayers = 7;
  std::vector<int> DeltaROF = {0, 0, 0, 0, 0, 0, 0}; // Delta in BC to define search window
  std::vector<float> LayerZ = {16.333f + 1, 16.333f + 1, 16.333f + 1, 42.140f + 1, 42.140f + 1, 73.745f + 1, 73.745f + 1};
  std::vector<float> LayerRadii = {2.33959f, 3.14076f, 3.91924f, 19.6213f, 24.5597f, 34.388f, 39.3329f};
  std::vector<float> LayerxX0 = {5.e-3f, 5.e-3f, 5.e-3f, 1.e-2f, 1.e-2f, 1.e-2f, 1.e-2f};
  std::vector<float> LayerResolution = {5.e-4f, 5.e-4f, 5.e-4f, 5.e-4f, 5.e-4f, 5.e-4f, 5.e-4f};
  std::vector<float> SystErrorY2 = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
  std::vector<float> SystErrorZ2 = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
  int ZBins{256};
  int PhiBins{128};
  bool UseDiamond = false;
  float Diamond[3] = {0.f, 0.f, 0.f};
  float DiamondCov[6] = {25.e-6f, 0.f, 0.f, 25.e-6f, 0.f, 36.f};

  /// General parameters
  bool AllowSharingFirstCluster = false;
  int ClusterSharing = 0;
  int MinTrackLength = 7;
  float NSigmaCut = 5;
  int NTimeSlices{1};
  /// Trackleting cuts
  float TrackletMinPt = 0.3f;
  /// Fitter parameters
  o2::base::PropagatorImpl<float>::MatCorrType CorrType = o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrNONE;
  float MaxChi2ClusterAttachment = 60.f;
  float MaxChi2NDF = 30.f;
  int ReseedIfShorter = 6; // reseed for the final fit track with the length shorter than this
  std::vector<float> MinPt = {0.f, 0.f, 0.f, 0.f};
  unsigned char StartLayerMask = 0x7F;
  bool ShiftRefToCluster = true; // TrackFit: after update shift the linearization reference to cluster
  bool FindShortTracks = false;
  bool PerPrimaryVertexProcessing = false;
  bool SaveTimeBenchmarks = false;
  bool DoUPCIteration = false;
  bool FataliseUponFailure = true;
  /// Cluster attachment
  bool UseTrackFollower = false;
  bool UseTrackFollowerTop = false;
  bool UseTrackFollowerBot = false;
  bool UseTrackFollowerMix = false;

  /// Seeding parameters
  float SeedingDCATolerance{0.8f};         // maximum allowed DCA to meanvertex for track to enter pool
  float SeedingDCAMaxPull{3.f};            // maximum allowed initial pull on DCA to meanvertex
  float SeedingMaxChi2Iter{3.f};           // maximum chi2 change required to end iteration
  float SeedingTukeyStartIter{5.f};        // start value for tukey scaling for iteration
  float SeedingMinWghTrk{0.01f};           // minimum weight a track has to have to be accounted
  float SeedingMinPtTrk{0.05};             // minimum pt for cells to enter pool
  int SeedingMinContrib{8};                // minimum number of contributors to account seeding vertex for rolling average
  int SeedingMaxFitIter{3};                // maximum iterations for fit
  int SeedingMinTracksIter{20};            // minimum tracks needed for one iteration
  int SeedingDBScanMinPt{50};              // DBSCAN: minimum number of cluster points
  float SeedingDBScanEpsZ{0.01f};          // DBSCAN: maximum epsilon z
  float SeedingDBScanEpsT{10.f};           // DBSCAN: maximum epsilon t (BC)
  float SeedingVertexExtraErr2[6] = {0.f}; // impose additional errors on seeding vertices

  bool createArtefactLabels{false};

  bool PrintMemory = false; // print allocator usage in epilog report
  size_t MaxMemory = std::numeric_limits<size_t>::max();
  bool DropTFUponFailure = false;
};

enum class RecoIterationSteps : uint16_t {
  /// which steps should be run
  kRunTrackleting,         // run the trackleting step
  kRunCellFinding,         // run the cell finding step (connect tracklets)
  kRunCellSeeding,         // run the cell seeding step (use cells to find seeding vertices)
  kRunCellNeighborFinding, // run the cell neighbor finding step (connect cells)
  kRunRoadFinding,         // run the road finding step (resolve ambiguities to find best roads/tracks)
  kRunTruthSeeding,        // run truth seeding (imposing MC event information as seeding vertices)
  /// extra steps
  kUpdateVertexTable, // update the vertex table for the current pool of vertices
  kUpdateClusters,    // update the cluster position wrt current beam constraint
  kInitMemory,        // initialize all vectors to use memory resource
};

struct RecoIteration {
  TrackingParameters params;
  utils::EnumFlags<RecoIterationSteps> steps;
  std::string name;
  std::string asString() const;
};

namespace TrackingMode
{
enum Type : int8_t {
  Unset = -1, // Special value to leave a default in case we want to override via Configurable Params
  Sync = 0,
  Async = 1,
  Cosmics = 2,
  Off = 3,
};

Type fromString(std::string_view str);
std::string toString(Type mode);

std::vector<RecoIteration> getRecoIterations(Type mode);
}; // namespace TrackingMode

} // namespace o2::its

#endif /* TRACKINGITSU_INCLUDE_CONFIGURATION_H_ */
