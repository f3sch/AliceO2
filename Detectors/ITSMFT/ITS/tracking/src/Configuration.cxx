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

#include <algorithm>
#include <format>
#include <limits>
#include <string_view>
#include <vector>

#include "Framework/Logger.h"
#include "ITStracking/Constants.h"
#include "ITStracking/Configuration.h"
#include "ITStracking/TrackingConfigParam.h"

using namespace o2::its;

std::string TrackingParameters::asString() const
{
  std::string str = std::format("NZb:{} NPhB:{} PerVtx:{} DropFail:{} ClSh:{} TrklMinPt:{:.2f} MinCl:{}", ZBins, PhiBins, PerPrimaryVertexProcessing, DropTFUponFailure, ClusterSharing, TrackletMinPt, MinTrackLength);
  bool first = true;
  for (int il = NLayers; il >= MinTrackLength; il--) {
    int slot = NLayers - il;
    if (slot < (int)MinPt.size() && MinPt[slot] > 0) {
      if (first) {
        first = false;
        str += " MinPt: ";
      }
      str += std::format("L{}:{:.2f} ", il, MinPt[slot]);
    }
  }
  str += " SystErrY/Z:";
  for (size_t i = 0; i < SystErrorY2.size(); i++) {
    str += std::format("{:.2e}/{:.2e} ", SystErrorY2[i], SystErrorZ2[i]);
  }
  str += std::format(" TimeSlices:{}", NTimeSlices);
  if (std::any_of(DeltaROF.begin(), DeltaROF.end(), [](auto b) { return b != 0; })) {
    str += " DeltaROF:[";
    for (int il = 0; il < NLayers; ++il) {
      str += std::to_string(DeltaROF[il]) + ",";
    }
    str[str.length() - 1] = ']';
  }
  if (std::numeric_limits<size_t>::max() != MaxMemory) {
    str += std::format(" MemLimit {:.2f} GB", double(MaxMemory) / constants::GB);
  }
  return str;
}

namespace
{
constexpr bool iequals(std::string_view a, std::string_view b)
{
  return std::equal(a.begin(), a.end(), b.begin(), b.end(),
                    [](char x, char y) { return std::tolower(x) == std::tolower(y); });
}
} // namespace

TrackingMode::Type TrackingMode::fromString(std::string_view str)
{
  constexpr std::array smodes = {
    std::pair{"sync", Sync},
    std::pair{"async", Async},
    std::pair{"cosmics", Cosmics},
    std::pair{"unset", Unset},
    std::pair{"off", Off}};

  auto it = std::find_if(smodes.begin(), smodes.end(), [&str](const auto& pair) {
    return iequals(str, pair.first);
  });
  if (it == smodes.end()) {
    LOGP(fatal, "Unrecognized tracking mode '{}'", str);
  }
  return it->second;
}

std::string TrackingMode::toString(TrackingMode::Type mode)
{
  if (mode == TrackingMode::Sync) {
    return "sync";
  } else if (mode == TrackingMode::Async) {
    return "async";
  } else if (mode == TrackingMode::Cosmics) {
    return "cosmics";
  } else if (mode == TrackingMode::Unset) {
    return "unset";
  } else if (mode == TrackingMode::Off) {
    return "off";
  }
  LOGP(fatal, "Unrecognized tracking mode '{}'", (int)mode);
  return ""; // not reachable
}

std::vector<RecoIteration> TrackingMode::getRecoIterations(TrackingMode::Type mode)
{
  const auto& tc = o2::its::TrackerParamConfig::Instance();
  std::vector<RecoIteration> recoIterations;

  // set the size and name
  if (mode == TrackingMode::Async) {
    recoIterations.resize(tc.doUPCIteration ? 5 : 4);
    recoIterations[0].name = "ASYNC_SEED";
    recoIterations[1].name = "ASYNC_LONG";
    recoIterations[2].name = "ASYNC_LONG_RETRY";
    recoIterations[3].name = "ASYNC_REST";
    if (tc.doUPCIteration) {
      recoIterations[4].name = "ASYNC_UPC";
    }
  } else if (mode == TrackingMode::Sync) {
    recoIterations.resize(2);
    recoIterations[0].name = "SYNC_TIGHT";
  } else if (mode == TrackingMode::Cosmics) {
    // in case of cosmics we do not do the seeding step
    recoIterations.resize(1);
    recoIterations[0].name = "COSMICS";
  } else {
    LOGP(fatal, "Unsupported ITS tracking mode {} ", toString(mode));
  }

  // always done in any first two iterations (unless cosmics)
  recoIterations[0].steps.set(RecoIterationSteps::kInitMemory);
  if (mode != TrackingMode::Cosmics) {
    recoIterations[1].steps.set(RecoIterationSteps::kInitMemory);

    // standards steps to configure seeding
    recoIterations[0].steps.set(RecoIterationSteps::kRunTrackleting, RecoIterationSteps::kRunCellFinding, RecoIterationSteps::kRunCellSeeding, RecoIterationSteps::kUpdateClusters);
    recoIterations[0].params.NLayers = 3;       // only do cell finding up until the third layer
    recoIterations[0].params.UseDiamond = true; // use the blown up diamond constrain to (e.g. luminous region)
    recoIterations[0].params.NSigmaCut = 5.f;
    recoIterations[0].params.CorrType = o2::base::PropagatorF::MatCorrType::USEMatCorrNONE; // do not use material
    recoIterations[0].params.SeedingDCATolerance = tc.seedingDCATolerance;
    recoIterations[0].params.SeedingDCAMaxPull = tc.seedingDCAMaxPull;
    recoIterations[0].params.SeedingMaxChi2Iter = tc.seedingMaxChi2Iter;
    recoIterations[0].params.SeedingTukeyStartIter = tc.seedingTukeyStartIter;
    recoIterations[0].params.SeedingMinWghTrk = tc.seedingMinWghTrk;
    recoIterations[0].params.SeedingMaxFitIter = tc.seedingMaxFitIter;
    recoIterations[0].params.SeedingMinTracksIter = tc.seedingMinTracksIter;
    recoIterations[0].params.SeedingDBScanMinPt = tc.seedingDBScanMinPt;
    recoIterations[0].params.SeedingDBScanEpsZ = tc.seedingDBScanEpsZ;
    recoIterations[0].params.SeedingDBScanEpsT = tc.seedingDBScanEpsT;
    recoIterations[0].params.PerPrimaryVertexProcessing = false;

    if (tc.seedingUseMCTruth) {
      recoIterations[0].name = "MC_SEEDING";
      recoIterations[0].steps.reset();
      recoIterations[0].steps.set(RecoIterationSteps::kRunTruthSeeding);
    }
  }

  if (mode == TrackingMode::Async) {
    recoIterations[2].params.TrackletMinPt = 0.2f;
    recoIterations[3].params.TrackletMinPt = 0.1f;

    recoIterations[1].params.MinPt[0] = 1.f / 12; // 7cl
    recoIterations[2].params.MinPt[0] = 1.f / 12; // 7cl

    recoIterations[3].params.MinTrackLength = 4;
    recoIterations[3].params.MinPt[0] = 1.f / 12; // 7cl
    recoIterations[3].params.MinPt[1] = 1.f / 5;  // 6cl
    recoIterations[3].params.MinPt[2] = 1.f / 1;  // 5cl
    recoIterations[3].params.MinPt[3] = 1.f / 6;  // 4cl

    recoIterations[3].params.StartLayerMask = (1 << 6) + (1 << 3);

    if (tc.doUPCIteration) {
      recoIterations[4].params.MinTrackLength = 4;
      recoIterations[4].params.TrackletMinPt = 0.1f;
    }
    for (size_t ip = 0; ip < recoIterations.size(); ip++) {
      // the seeding step is configured outside of this loop beforehand
      if (ip > 0) {
        recoIterations[ip].steps.set(RecoIterationSteps::kRunTrackleting, RecoIterationSteps::kRunCellFinding, RecoIterationSteps::kRunCellNeighborFinding, RecoIterationSteps::kRunRoadFinding);
        if (ip == 1) {
          recoIterations[ip].steps.set(RecoIterationSteps::kUpdateClusters, RecoIterationSteps::kUpdateVertexTable);
        }
      }
      recoIterations[ip].params.ZBins = 64;
      recoIterations[ip].params.PhiBins = 32;
      // check if something was overridden via configurable params
      if (ip < tc.MaxIter) {
        if (tc.startLayerMask[ip] > 0) {
          recoIterations[2].params.StartLayerMask = tc.startLayerMask[ip];
        }
        if (tc.minTrackLgtIter[ip] > 0) {
          recoIterations[ip].params.MinTrackLength = tc.minTrackLgtIter[ip];
        }
        for (int ilg = tc.MaxTrackLength; ilg >= tc.MinTrackLength; ilg--) {
          int lslot0 = (tc.MaxTrackLength - ilg), lslot = lslot0 + (static_cast<int>(ip) * (tc.MaxTrackLength - tc.MinTrackLength + 1));
          if (tc.minPtIterLgt[lslot] > 0.) {
            recoIterations[ip].params.MinPt[lslot0] = tc.minPtIterLgt[lslot];
          }
        }
      }
    }
  } else if (mode == TrackingMode::Sync) {
    recoIterations[0].params.ZBins = 64;
    recoIterations[0].params.PhiBins = 32;
    recoIterations[0].params.MinTrackLength = 4;
  } else if (mode == TrackingMode::Cosmics) {
    recoIterations[0].params.MinTrackLength = 4;
    recoIterations[0].params.PhiBins = 4;
    recoIterations[0].params.ZBins = 16;
    recoIterations[0].params.MaxChi2ClusterAttachment = 60.;
    recoIterations[0].params.MaxChi2NDF = 40.;
  }

  // global parameters set for every iteration
  float bFactor = std::abs(o2::base::Propagator::Instance()->getNominalBz()) / 5.0066791f;
  float bFactorTracklets = bFactor < 0.01f ? 1.f : bFactor; // for tracklets only
  for (auto& reco : recoIterations) {
    auto& p = reco.params;

    // adjust pT settings to actual mag. field
    p.TrackletMinPt *= bFactorTracklets;
    for (int ilg = tc.MaxTrackLength; ilg >= tc.MinTrackLength; ilg--) {
      int lslot = tc.MaxTrackLength - ilg;
      p.MinPt[lslot] *= bFactor;
    }
    p.ReseedIfShorter = tc.reseedIfShorter;
    p.ShiftRefToCluster = tc.shiftRefToCluster;
    p.createArtefactLabels = tc.createArtefactLabels;
    p.NTimeSlices = tc.nTimeSlices;
    p.PrintMemory = tc.printMemory;
    p.MaxMemory = tc.maxMemory;
    p.DropTFUponFailure = tc.dropTFUponFailure;
    p.SaveTimeBenchmarks = tc.saveTimeBenchmarks;
    p.FataliseUponFailure = tc.fataliseUponFailure;
    p.AllowSharingFirstCluster = tc.allowSharingFirstCluster;

    if (tc.useMatCorrTGeo) {
      p.CorrType = o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrTGeo;
    } else if (tc.useFastMaterial) {
      p.CorrType = o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrNONE;
    } else {
      p.CorrType = o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrLUT;
    }

    for (int i{0}; i < tc.MaxTrackLength; ++i) {
      p.SystErrorY2[i] = tc.sysErrY2[i] > 0 ? tc.sysErrY2[i] : p.SystErrorY2[i];
      p.SystErrorZ2[i] = tc.sysErrZ2[i] > 0 ? tc.sysErrZ2[i] : p.SystErrorZ2[i];
      p.DeltaROF[i] = tc.deltaROF[i];
    }

    p.DoUPCIteration = tc.doUPCIteration;
    p.MaxChi2ClusterAttachment = tc.maxChi2ClusterAttachment > 0 ? tc.maxChi2ClusterAttachment : p.MaxChi2ClusterAttachment;
    p.MaxChi2NDF = tc.maxChi2NDF > 0 ? tc.maxChi2NDF : p.MaxChi2NDF;
    p.PhiBins = tc.LUTbinsPhi > 0 ? tc.LUTbinsPhi : p.PhiBins;
    p.ZBins = tc.LUTbinsZ > 0 ? tc.LUTbinsZ : p.ZBins;
    p.NSigmaCut *= tc.nSigmaCut > 0 ? tc.nSigmaCut : 1.f;
    p.TrackletMinPt *= tc.minPt > 0 ? tc.minPt : 1.f;
    p.PerPrimaryVertexProcessing = tc.perPrimaryVertexProcessing;
    std::copy(tc.diamondPos, tc.diamondPos + 3, p.Diamond);
    std::copy(tc.diamondCov, tc.diamondCov + 6, p.DiamondCov);
    if (tc.useTrackFollower > 0) {
      p.UseTrackFollower = true;
      // Bit 0: Allow for mixing of top&bot extension --> implies Bits 1&2 set
      // Bit 1: Allow for top extension
      // Bit 2: Allow for bot extension
      p.UseTrackFollowerMix = ((tc.useTrackFollower & (1 << 0)) != 0);
      p.UseTrackFollowerTop = ((tc.useTrackFollower & (1 << 1)) != 0);
      p.UseTrackFollowerBot = ((tc.useTrackFollower & (1 << 2)) != 0);
    }
    if (tc.findShortTracks >= 0) {
      p.FindShortTracks = tc.findShortTracks;
    }
  }

  // opt. suppress layer iterations
  if (tc.nIterations >= 0 && tc.nIterations < recoIterations.size()) {
    recoIterations.resize(tc.nIterations);
  }

  return recoIterations;
}

std::string RecoIteration::asString() const
{
  return std::format("recoIter:{}[{}] {}", name, steps.string(), params.asString());
}
