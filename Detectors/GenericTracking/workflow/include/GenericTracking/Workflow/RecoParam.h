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

/// \file RecoParam.h
/// \brief Each detector's tracking settings as a ConfigurableParam.
///
/// A ConfigurableParam, where `--configKeyValues "GTRKReco<DET>.<field>=<value>"`
/// reaches them and every device of the chain sees the same value.
///
/// The layouts share the `RecoSettings` fields. Each layout reads its own
/// parameter through `recoParam<Det>()`: `GTRKRecoITS`, `GTRKRecoMFT`,
/// `GTRKRecoTRK` or `GTRKRecoTRKTOF`. Layer-indexed keys use the selected
/// layout's numbering.

#ifndef GTRK_WORKFLOW_RECOPARAM_H_
#define GTRK_WORKFLOW_RECOPARAM_H_

#include "GenericTracking/Base/LayerSpec.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include "DetectorsCommonDataFormats/DetID.h"

// ENABLE_UPGRADES comes with DetID.h (UpgradesStatus.h).
#ifdef ENABLE_UPGRADES
#include "GenericTracking/Detectors/ITS3.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Detectors/TRKFT3TOF.h"
#endif

#include <cstdint>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace o2::gtrk
{

struct RecoSettings {
  /// O2's `addTimeError` for the ITS and the MFT [BC]: every readout frame's
  /// time box widened by this much at both ends, for a hit recorded in the
  /// neighbouring frame. One entry per layer, or
  /// none for 0 on every layer, O2's default.
  std::vector<std::int32_t> frameMarginBC;

  /// Width of the outer time block in the measurement grid [ps], 0 for one
  /// block. This is the minimum slice quantum and is chosen from the worker
  /// memory budget; a detector with readout frames takes the whole number of
  /// frames nearest it.
  std::int64_t blockWidthPS = 10'000'000;

  // The grid's cells are derived from the iterations' queries with the
  // material LUT's scattering, with no key of their own (GridSizing.h).
  /// Added in quadrature to every cluster's local errors [cm], for what the
  /// cluster shape does not know: the alignment's residuals on real data. O2
  /// adds 10 um per layer in its async reconstruction (ITSCATrackerParam.sysErrY2,
  /// sysErrZ2 = 100e-8 cm^2); 0, the default, is the MC's dictionary error alone.
  double clusterSystematicSigmaCM = 0.;
  /// The windows' and the cuts' nSigma, the run's (every iteration's).
  double windowNSigma = 3.;
  /// How many sigma of the beam's z width the tracking reaches about the
  /// beam's z (`beamSpotConstraint`): the z range it is to be efficient in,
  /// for the beam-line windows and pair cut, the vertexing pass, the grid and
  /// the vertex-list iterations alike.
  double beamZReachNSigma = 3.;

  /// The own vertexer: before the tracking, a pass at
  /// the beam line on the three innermost layers, one slice over the
  /// timeframe, tracklets and triplets at `vertexerMinPtGeV`; its triplets'
  /// z at the beam line, frame by frame, give the vertex list the
  /// vertex-list iterations read (`iterationVertexSource`). Off by default.
  /// The rest are the vertexer's knobs (VertexStage.h, VertexerSettings).
  bool vertexer = false;
  double vertexerMinPtGeV = 0.5;
  double vertexerMaxChi2RZ = 5.;   ///< a triplet's r-z kink chi2 cut, with its scattering at its own momentum; 0 for none
  double vertexerMaxDcaCM = 0.02;  ///< its circle's transverse distance to the beam axis [cm], 0 for none
  std::int32_t vertexerWindowUM = 150;      ///< a peak's window [um]
  std::int32_t vertexerClaimUM = 2000;      ///< the half-width a vertex claims around its peak [um]
  std::int32_t vertexerMinContributors = 4; ///< the least candidates a peak's window holds
  /// The least in a frame where none holds vertexerMinContributors: that
  /// frame's rounds run again at this minimum; 0, none, by default.
  std::int32_t vertexerAloneMinContributors = 0;
  double vertexerSignificance = 5.; ///< a peak's excess over the sidebands, in their Poisson sigmas
  /// Whether the vertexing pass tests its pairs and triplets against their
  /// measurements' cluster shapes (`iterationShapeTest`); 0, off, by default.
  std::int32_t vertexerShapeTest = 0;

  /// Tracking iterations, run one after another: one per entry of
  /// `iterationMinPtGeV`. Each other per-iteration vector holds one entry per
  /// iteration, or none for its default in every iteration (GridConfig.h,
  /// `iterationCount`). They are set whole, as O2's parameters take vectors:
  /// `GTRKRecoITS.iterationMinPtGeV=[1,0.5];GTRKRecoITS.iterationMinClusters=[7,4]`.
  /// What an iteration loosens is its minimum pT; nSigma (`windowNSigma`) is
  /// the run's, the same in every iteration.
  ///
  /// Three iterations at 1.0, 0.5 and 0.3 GeV by default; a schedule sets its
  /// own.
  std::vector<double> iterationMinPtGeV{1., 0.5, 0.3};
  /// Where each iteration's vertex seeds come from: 0 the beam line (the
  /// default), 1 the timeframe's vertex list, 2 none: straight legs from
  /// anywhere, cosmic muons (barrel links only, no beam constraint).
  std::vector<std::int32_t> iterationVertexSource;
  /// The shortest road each iteration fits, in clusters, 4 by default; a road
  /// is at least 4 clusters, two triplets, whatever is set. An iteration takes
  /// only the links on a path of this many layers (Iteration.h, `takenLinks`).
  /// ITS asks 7 of its first two iterations and 4 of the third.
  std::vector<std::int32_t> iterationMinClusters;
  /// The kinks each iteration's fit may take per pass;
  /// 0 for none, 1 by default.
  std::vector<std::int32_t> iterationMaxKinks;
  /// The kink allowance of each iteration's triplet and neighbour cuts: a
  /// triplet failing the nominal cuts is kept, flagged kinked, if it passes with the middle
  /// layer's scattering at pT_min scaled by this factor; 1 (the default) for
  /// none.
  std::vector<double> iterationKinkSigma;
  /// The triplet and neighbour cuts' own kink allowance, each iteration's:
  /// above 0 it replaces
  /// iterationKinkSigma there, which then sets only the road fit's and the
  /// follow's; 0, the default, iterationKinkSigma everywhere.
  std::vector<double> iterationGraphKinkSigma;
  /// Whether each iteration makes tracklets that leave a hole, whose segment
  /// crosses the envelope of another layer of its mask (the default); 0 for
  /// none, so the tracks past a shorter layer's end are kept and those
  /// through a gap are not.
  std::vector<std::int32_t> iterationSkipLinks;
  /// The layers each iteration's seeds take rows on, a bit per layer; 0 for
  /// every layer (the default): its links join two of them, and of the layers
  /// a link passes over those among them are holes (`iterationSkipLinks`).
  /// The vertexing pass's 7 takes layers 0-2 alone.
  std::vector<long> iterationLayerMask;
  /// The layers each iteration's roads may end on, a bit per layer; 0 for any
  /// (the default). O2's ITS async asks layer 3 or 6 of its third,
  /// `(1 << 6) + (1 << 3)` = 72.
  std::vector<long> iterationRoadEndLayers;
  /// The layers each iteration's roads may start on, their innermost
  /// measurement's, a bit per layer; 0 for any (the default).
  std::vector<long> iterationRoadStartLayers;
  /// The lowest pT of a track by its cluster count [GeV], 0 for none -- the
  /// final track's: after the follow where the iteration follows, else the
  /// road's. Iteration
  /// `i`'s floor for `n` clusters is entry `18 i + n`, so 18 entries per
  /// iteration (a slot for every chain of the longest layout), or none for no
  /// floor at all. O2's ITS async asks 1/12 of 7 clusters in every iteration,
  /// and 1/5, 1 and 1/6 of 6, 5 and 4 in its third (entries 7, 25, 43, 42,
  /// 41, 40).
  static constexpr std::uint32_t kClusterSlots = 18U;
  std::vector<double> iterationClusterMinPtGeV;
  /// How many of a track's innermost clusters each iteration lets it share
  /// with another; 0, strict, by default.
  std::vector<std::int32_t> iterationSharedInnermost;
  /// Whether each iteration follows its selected tracks after their refit,
  /// outward and then inward, as a CKF does, before the next iteration;
  /// 1, on, by default.
  std::vector<std::int32_t> iterationFollow;
  /// The holes each iteration's follow lets a track cross in each direction,
  /// a sensor crossed without a row; 0 by default.
  std::vector<std::int32_t> iterationFollowHoles;
  /// The innermost clusters each iteration's road fit may drop from a road it
  /// refuses, one more each time, down to the iteration's minimum, the first
  /// it takes the road's record; 0, off, by
  /// default.
  std::vector<std::int32_t> iterationInnerTrim;
  /// Whether each iteration tests its pairs and triplets against their
  /// measurements' cluster shapes: a candidate's local slope at each of its measurements -- its line's
  /// |dz| over its length along the sensor's normal -- must lie in the
  /// measurement's slope interval, which the measurement producer makes from
  /// its detector's calibrated band (ClusterShapeBand.h); 0, off, by default.
  std::vector<std::int32_t> iterationShapeTest;
  /// Whether each iteration's road fit also weighs its rows' cluster spans:
  /// -2 log of how much more likely each row's spans are for the road's
  /// crossing than for any row of its layer (the measurement producer's
  /// shares, from the detector's calibrated band) added to the fit's chi2,
  /// held to the road's cut and carried into the record's chi2; nothing on a
  /// layout without a band; 0, off, by default.
  std::vector<std::int32_t> iterationRoadSize;
  /// Each iteration's lowest momentum [GeV] of a disc->disc pair's track,
  /// beside its minimum pT. Above 0, on a prompt beam line only, the
  /// iteration searches its disc->disc pairs in a window per source -- from
  /// the source's direction from the beam line, the curvatures and momenta
  /// the floors admit, the scattering at them and the layers' resolutions --
  /// instead of the link's, the window then part of the pair's acceptance
  /// (TrackletStage.h, `discPairWindow`); 0, the link's window, by default.
  std::vector<double> iterationMinPGeV;
  /// Each iteration's follow's end decisions on the last row's own update:
  /// inward, the sqrt|V| [cm^2]
  /// above which a last row at the layout's end is refused; outward, the
  /// density [cm^-2] above which the unconfirmed drop keeps a last row; 0, off, by default.
  std::vector<double> iterationFollowEndMaxSqrtDet;
  std::vector<double> iterationFollowKeepDensity;
  /// The clusters each iteration's tracks must hold after the follow; 0 for
  /// none (the default). A shorter track is dropped and its rows given back.
  std::vector<std::int32_t> iterationMinClustersAfterFollow;
  /// Whether each iteration updates its roads with the beam line on their
  /// innermost surface and refuses them above the update's gate, to start
  /// the follow inward from the constrained state; 0, off, by default.
  std::vector<std::int32_t> iterationBeamConstraint;
  /// Whether each iteration follows every earlier iteration's tracks too,
  /// with its own settings, their rows kept, as a recovery pass does; 0, off,
  /// by default. Needs the iteration's
  /// iterationFollow; iterationMinClustersAfterFollow applies to the
  /// iteration's own tracks only.
  std::vector<std::int32_t> iterationFollowEarlier;

  /// Requested width of one time slice [ps]. The actual width is rounded up
  /// to whole grid blocks, and zero means one slice over the timeframe.
  std::int64_t sliceWidthPS = 10'000'000;

  /// Maximum CPU workers; zero uses TBB's default.
  std::int32_t threads = 0;

  /// The host tracker's memory [GB]: a stage that would hold more fails with
  /// kOutOfMemory, and the timeframe's tracking with it, instead of growing.
  /// Zero, the default, for no limit.
  double memoryLimitGB = 0.;

  /// A timeframe whose tracking runs out of memory (memoryLimitGB) gives
  /// empty products and is counted, instead of failing the run. Off, the
  /// default, as ITS's DropTFUponFailure: the run fails.
  bool dropTFUponFailure = false;
};

/// Detector \a N's parameter name: `GTRKReco<DET>`, ITS3's `GTRKRecoITS3`.
[[nodiscard]] constexpr std::string_view recoParamName(std::int32_t n)
{
#ifdef ENABLE_UPGRADES
  if (n == o2::detectors::DetID::IT3) {
    return "GTRKRecoITS3";
  }
#endif
  return n == o2::detectors::DetID::ITS ? "GTRKRecoITS" : (n == o2::detectors::DetID::MFT ? "GTRKRecoMFT" : "GTRKRecoTRK");
}

/// Detector \a N's settings (`o2::detectors::DetID`): `GTRKReco<DET>`.
// O2ParamDef deletes the copy constructor and supplies the rest of the
// boilerplate, so the rule about declaring all five has nothing to act on here.
template <std::int32_t N>
// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions)
struct RecoParam : public o2::conf::ConfigurableParamPromoter<RecoParam<N>, RecoSettings> {
  static constexpr std::string_view getParamName() { return recoParamName(N); }
  O2ParamDef(RecoParam, getParamName().data())
};

template <std::int32_t N>
RecoParam<N> RecoParam<N>::sInstance;

/// A detector layout's O2 detector and data origin (DataSpecs.h).
template <DetectorLayout Det>
struct DetIDOf;
template <>
struct DetIDOf<detectors::ITS> : std::integral_constant<o2::detectors::DetID::ID, o2::detectors::DetID::ITS> {
};
template <>
struct DetIDOf<detectors::MFT> : std::integral_constant<o2::detectors::DetID::ID, o2::detectors::DetID::MFT> {
};
#ifdef ENABLE_UPGRADES
/// ITS3 is O2's IT3: its own messages (GTRK/IT3...) and keys; O2's ITS3
/// clusters, geometry and frames come with the ITS's origin
/// (WorkflowITS3/ConditionsITS3.h).
template <>
struct DetIDOf<detectors::ITS3> : std::integral_constant<o2::detectors::DetID::ID, o2::detectors::DetID::IT3> {
};
template <>
struct DetIDOf<detectors::TRKFT3> : std::integral_constant<o2::detectors::DetID::ID, o2::detectors::DetID::TRK> {
};
/// TRKFT3TOF publishes tracks with TRK's data origin.
template <>
struct DetIDOf<detectors::TRKFT3TOF> : std::integral_constant<o2::detectors::DetID::ID, o2::detectors::DetID::TRK> {
};
#endif

/// A detector layout's reconstruction parameters.
template <DetectorLayout Det>
struct RecoParamOf {
  using type = RecoParam<DetIDOf<Det>::value>; // NOLINT(readability-identifier-naming): a type trait's member
};

#ifdef ENABLE_UPGRADES
/// TRKFT3TOF's settings, `GTRKRecoTRKTOF`: TRK's detector with its own keys.
// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions): as RecoParam's
struct RecoParamTRKTOF : public o2::conf::ConfigurableParamPromoter<RecoParamTRKTOF, RecoSettings> {
  O2ParamDef(RecoParamTRKTOF, "GTRKRecoTRKTOF")
};

template <>
struct RecoParamOf<detectors::TRKFT3TOF> {
  using type = RecoParamTRKTOF; // NOLINT(readability-identifier-naming): a type trait's member
};
#endif

template <DetectorLayout Det>
[[nodiscard]] const RecoSettings& recoParam()
{
  return RecoParamOf<Det>::type::Instance();
}

} // namespace o2::gtrk

#endif // GTRK_WORKFLOW_RECOPARAM_H_
