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

#ifndef O2_TRACKING_STUDY_TYPES_H
#define O2_TRACKING_STUDY_TYPES_H
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "ReconstructionDataFormats/Track.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCEventLabel.h"
#include "CommonConstants/LHCConstants.h"
#include "CommonDataFormat/TimeStamp.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/V0.h"
#include "SimulationDataFormat/TrackReference.h"
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

namespace o2::trackstudy
{
struct MCTrackInfo {

  inline float getMCTimeMUS() const { return bcInTF * o2::constants::lhc::LHCBunchSpacingMUS; }
  inline bool hasITSHitOnLr(int i) const { return (pattITSCl & ((0x1 << i) & 0x7f)) != 0; }
  int getNITSClusCont() const;
  int getNITSClusForAB() const;
  int getLowestITSLayer() const;
  int getHighestITSLayer() const;
  std::vector<float> occTPCV{};
  std::vector<o2::track::TrackPar> trackRefsTPC{};
  o2::track::TrackPar track{};
  o2::MCCompLabel label{};
  float occTPC = -1.f;
  int occITS = -1.f;
  int bcInTF = -1;
  int pdg = 0;
  int pdgParent = 0;
  int parentEntry = -1;
  int16_t nTPCCl = 0;
  int16_t nTPCClShared = 0;
  int8_t parentDecID = -1;
  uint8_t minTPCRow = -1;
  uint8_t maxTPCRow = 0;
  uint8_t nUsedPadRows = 0;
  uint8_t maxTPCRowInner = 0; // highest row in the sector containing the lowest one
  uint8_t minTPCRowSect = -1;
  uint8_t maxTPCRowSect = -1;
  int8_t nITSCl = 0;
  int8_t pattITSCl = 0;
  uint8_t flags = 0;

  enum Flags : uint32_t { Primary = 0,
                          AddedAtRecStage = 2,
                          BitMask = 0xff };

  bool isPrimary() const { return isBitSet(Primary); }
  bool isAddedAtRecStage() const { return isBitSet(AddedAtRecStage); }
  void setPrimary() { setBit(Primary); }
  void setAddedAtRecStage() { setBit(AddedAtRecStage); }

  uint8_t getBits() const { return flags; }
  bool isBitSet(int bit) const { return flags & (0xff & (0x1 << bit)); }
  void setBits(std::uint8_t b) { flags = b; }
  void setBit(int bit) { flags |= BitMask & (0x1 << bit); }
  void resetBit(int bit) { flags &= ~(BitMask & (0x1 << bit)); }

  o2::track::TrackPar getTrackParTPC(float b, float x = 90) const;
  float getTrackParTPCPar(int i, float b, float x = 90) const;
  float getTrackParTPCPhiSec(float b, float x = 90) const;

  ClassDefNV(MCTrackInfo, 8);
};

struct RecTrack {
  enum FakeFlag {
    FakeITS = 0x1 << 0,
    FakeTPC = 0x1 << 1,
    FakeTRD = 0x1 << 2,
    FakeTOF = 0x1 << 3,
    FakeITSTPC = 0x1 << 4,
    FakeITSTPCTRD = 0x1 << 5,
    HASACSides = 0x1 << 6,
    FakeGLO = 0x1 << 7
  };
  o2::track::TrackParCov track{};
  o2::dataformats::VtxTrackIndex gid{};
  o2::dataformats::TimeStampWithError<float, float> ts{};
  o2::MCEventLabel pvLabel{};
  short pvID = -1;
  uint8_t nClTPCShared = 0;
  uint8_t flags = 0;
  uint8_t nClITS = 0;
  uint8_t nClTPC = 0;
  uint8_t pattITS = 0;
  int8_t lowestPadRow = -1;
  int8_t padFromEdge = -1;
  uint8_t rowMaxTPC = 0;
  uint8_t rowCountTPC = 0;

  bool isFakeGLO() const { return flags & FakeGLO; }
  bool isFakeITS() const { return flags & FakeITS; }
  bool isFakeTPC() const { return flags & FakeTPC; }
  bool isFakeTRD() const { return flags & FakeTRD; }
  bool isFakeTOF() const { return flags & FakeTOF; }
  bool isFakeITSTPC() const { return flags & FakeITSTPC; }
  bool hasACSides() const { return flags & HASACSides; }

  ClassDefNV(RecTrack, 3);
};

struct TrackPairInfo {
  RecTrack tr0;
  RecTrack tr1;
  uint8_t nshTPC = 0;
  uint8_t nshTPCRow = 0;

  int getComb() const { return tr0.track.getSign() != tr1.track.getSign() ? 0 : (tr0.track.getSign() > 0 ? 1 : 2); }
  float getDPhi() const
  {
    float dphi = tr0.track.getPhi() - tr1.track.getPhi();
    if (dphi < -o2::constants::math::PI) {
      dphi += o2::constants::math::TwoPI;
    } else if (dphi > o2::constants::math::PI) {
      dphi -= o2::constants::math::TwoPI;
    }
    return dphi;
  }
  float getDTgl() const { return tr0.track.getTgl() - tr1.track.getTgl(); }

  ClassDefNV(TrackPairInfo, 1)
};

struct TrackFamily { // set of tracks related to the same MC label
  MCTrackInfo mcTrackInfo{};
  std::vector<RecTrack> recTracks{};
  o2::track::TrackParCov trackITSProp{};
  o2::track::TrackParCov trackTPCProp{};
  int8_t entITS = -1;
  int8_t entTPC = -1;
  int8_t entITSTPC = -1;
  int8_t entITSFound = -1; // ITS track for this MC track, regardless if it was matched to TPC of another track
  int8_t flags = 0;
  float tpcT0 = -999.;

  bool contains(const o2::dataformats::VtxTrackIndex& ref) const
  {
    for (const auto& tr : recTracks) {
      if (ref == tr.gid) {
        return true;
      }
    }
    return false;
  }
  const RecTrack& getTrackWithITS() const { return entITS < 0 ? dummyRecTrack : recTracks[entITS]; }
  const RecTrack& getTrackWithTPC() const { return entTPC < 0 ? dummyRecTrack : recTracks[entTPC]; }
  const RecTrack& getTrackWithITSTPC() const { return entITSTPC < 0 ? dummyRecTrack : recTracks[entITSTPC]; }
  const RecTrack& getTrackWithITSFound() const { return entITSFound < 0 ? dummyRecTrack : recTracks[entITSFound]; }
  const RecTrack& getLongestTPCTrack() const
  {
    int n = getLongestTPCTrackEntry();
    return n < 0 ? dummyRecTrack : recTracks[n];
  }
  int getLongestTPCTrackEntry() const;
  int getNTPCClones() const;
  static RecTrack dummyRecTrack; //

  ClassDefNV(TrackFamily, 1);
};

struct ClResTPCCont {
  // contributor to TPC Cluster
  std::array<float, 3> xyz{};
  std::array<float, 3> below{};
  std::array<float, 3> above{};
  float snp = 0.;
  float tgl = 0.;
  float q2pt = 0.;
  bool corrAttach = false;

  int getNExt() const { return (below[0] > 1.) + (above[0] > 1.); }

  float getClX() const { return xyz[0]; }
  float getClY() const { return xyz[1]; }
  float getClZ() const { return xyz[2]; }

  float getDY() const { return xyz[1] - getYRef(); }
  float getDZ() const { return xyz[2] - getZRef(); }

  float getYRef() const
  {
    float y = 0;
    int n = 0;
    if (below[0] > 1.) {
      y += below[1];
      n++;
    }
    if (above[0] > 1.) {
      y += above[1];
      n++;
    }
    return n == 1 ? y : 0.5 * y;
  }

  float getZRef() const
  {
    float z = 0;
    int n = 0;
    if (below[0] > 1.) {
      z += below[2];
      n++;
    }
    if (above[0] > 1.) {
      z += above[2];
      n++;
    }
    return n == 1 ? z : 0.5 * z;
  }

  float getDXMin() const
  {
    float adxA = 1e9, adxB = 1e9;
    if (above[0] > 1.) {
      adxA = xyz[0] - above[0];
    }
    if (below[0] > 1.) {
      adxB = xyz[1] - below[0];
    }
    return std::abs(adxA) < std::abs(adxB) ? adxA : adxB;
  }

  float getDXMax() const
  {
    float adxA = 0, adxB = 0;
    if (above[0] > 1.) {
      adxA = xyz[0] - above[0];
    }
    if (below[0] > 1.) {
      adxB = xyz[0] - below[0];
    }
    return std::abs(adxA) > std::abs(adxB) ? adxA : adxB;
  }

  float getEY() const { return getNExt() > 1 ? below[1] - above[1] : -999; }
  float getEZ() const { return getNExt() > 1 ? below[2] - above[2] : -999; }

  ClassDefNV(ClResTPCCont, 1);
};

struct ClResTPC {
  uint8_t sect = 0;
  uint8_t row = 0;
  uint8_t ncont = 0;
  uint8_t flags = 0;
  uint8_t sigmaTimePacked;
  uint8_t sigmaPadPacked;
  float qmax = 0;
  float qtot = 0;
  float occ = 0;
  float occBin = 0;
  float getSigmaPad() const { return float(sigmaPadPacked) * (1.f / 32); }
  float getSigmaTime() const { return float(sigmaTimePacked) * (1.f / 32); }

  std::vector<ClResTPCCont> contTracks;
  int getNCont() const { return contTracks.size(); }

  float getDY(int i) const { return i < getNCont() ? contTracks[i].getDY() : -999.; }
  float getDZ(int i) const { return i < getNCont() ? contTracks[i].getDZ() : -999.; }
  float getYRef(int i) const { return i < getNCont() ? contTracks[i].getYRef() : -999.; }
  float getZRef(int i) const { return i < getNCont() ? contTracks[i].getZRef() : -999.; }
  float getDXMin(int i) const { return i < getNCont() ? contTracks[i].getDXMin() : -999.; }
  float getDXMax(int i) const { return i < getNCont() ? contTracks[i].getDXMax() : -999.; }
  float getEY(int i) const { return i < getNCont() ? contTracks[i].getEY() : -999.; }
  float getEZ(int i) const { return i < getNCont() ? contTracks[i].getEZ() : -999.; }

  void sortCont()
  {
    std::sort(contTracks.begin(), contTracks.end(), [](const ClResTPCCont& a, const ClResTPCCont& b) {
      float dya = a.getDY(), dyb = b.getDY(), dza = a.getDZ(), dzb = b.getDZ();
      return dya * dya + dza * dza < dyb * dyb + dzb * dzb;
    });
  }

  ClassDefNV(ClResTPC, 2);
};

struct ITSHitInfo {
  o2::BaseCluster<float> clus{};
  o2::TrackReference tref{};
  float trefXT = 0; // track ref tracking frame coordinates
  float trefYT = 0;
  float chipX = 0;
  float chipAlpha = 0;
  ClassDefNV(ITSHitInfo, 1);
};

struct RecPV {
  o2::dataformats::PrimaryVertex pv{};
  o2::MCEventLabel mcEvLbl{};
  ClassDefNV(RecPV, 1);
};

struct MCVertex {
  float getX() const { return pos[0]; }
  float getY() const { return pos[1]; }
  float getZ() const { return pos[2]; }

  std::array<float, 3> pos{0., 0., -1999.f};
  float ts = 0;
  int nTrackSel = 0; // number of selected MC charged tracks
  int ID = -1;
  std::vector<RecPV> recVtx{};
  std::vector<float> occTPCV{};
  ClassDefNV(MCVertex, 2);
};

/// State of one prong of a MC decay with respect to the SVertexer seeds pool.
/// A prong which was reconstructed but never made it into the pool can never form a V0,
/// whatever the pair cuts do, so this has to be checked before interpreting SVCheck::rejV0.
struct SVProngInfo {
  o2::dataformats::VtxTrackIndex gid; // reco track used as this prong, unset if not reconstructed
  int32_t poolEntry = -1;             // entry in the SVertexer seeds pool, -1 if not seeded
  int8_t poolSide = -1;               // SVertexer::POS / SVertexer::NEG
  int32_t vBrMin = -1;                // vertex bracket of the seed
  int32_t vBrMax = -1;
  float minR = -1.f;   // lowest radial point of the seed, used by the causality cut
  uint8_t seedRej = 0; // SVertexer::SeedRej, why the track never entered the pool
  int8_t nITSclu = -1;
  bool hasTPC = false;

  bool isReconstructed() const { return gid.isSourceSet(); }
  bool isSeeded() const { return poolEntry >= 0; }

  ClassDefNV(SVProngInfo, 1);
};

/// Why a MC decay was or was not reconstructed as a V0 by the SVertexer.
/// stage says how far the decay got, and only if it reached CutRejected is rejV0 meaningful.
struct SVCheck {
  enum Stage : int8_t {
    NotChecked = -1,  // SV checking disabled or decay not eligible
    NoProngs,         // at least one prong has no reconstructed track at all
    NotSeeded,        // a prong was reconstructed but rejected before the seeds pool, see prong seedRej
    SameCharge,       // prongs did not end up as one positive and one negative seed
    NoBracketOverlap, // seeds share no primary vertex, so the pair is never even tried
    CutRejected,      // the pair was tried and rejected, see rejV0
    Found             // the pair passes the SVertexer selection
  };
  std::array<SVProngInfo, 2> prongs{};
  int foundSVID = -1;           // reconstructed V0 matched to this decay by MC labels, -1 if none
  int8_t stage = NotChecked;    // Stage
  uint8_t rejV0 = 0;            // SVertexer::V0Rej, only if stage == CutRejected
  bool pairInReco = false;      // the reconstruction built a V0 out of exactly the two prongs replayed
  bool replayConsistent = true; // false if the replay verdict differs from pairInReco

  bool isReconstructed() const { return foundSVID >= 0; }

  bool isFound() const { return stage == Found; }
  bool bothProngsReconstructed() const { return prongs[0].isReconstructed() && prongs[1].isReconstructed(); }
  bool bothProngsSeeded() const { return prongs[0].isSeeded() && prongs[1].isSeeded(); }

  ClassDefNV(SVCheck, 1);
};

/// A reconstructed V0 together with the MC origin of its prongs, to study the composition of
/// the sample: which V0s are real decays and which are combinatorial.
struct RecSVInfo {
  enum Kind : int8_t {
    Unknown = -1,     // MC information unavailable for at least one prong
    TrueDecay,        // both prongs are daughters of the same MC mother
    DifferentMothers, // prongs come from unrelated MC particles, i.e. combinatorial
    FakeProng         // at least one prong track has a fake MC label
  };
  o2::dataformats::V0 v0;
  o2::dataformats::V0Index v0ID{};
  std::array<o2::MCCompLabel, 2> prongLbl{};
  std::array<int, 2> prongPDG{0, 0};
  int mcMotherPDG = 0;       // PDG of the common mother, 0 if there is none
  int mcMotherEntry = -1;    // entry in the decays pool of mcMotherDecID, -1 if not checked
  int8_t mcMotherDecID = -1; // which decay type the pool of mcMotherEntry refers to
  int8_t kind = Unknown;     // Kind

  bool isTrueDecay() const { return kind == TrueDecay; }
  bool isCheckedDecay() const { return mcMotherEntry >= 0; }

  ClassDefNV(RecSVInfo, 1);
};

/// One reconstructed representation of an exact truth conversion daughter.
/// rejection is exclusive: Accepted means that this track contributes to the efficiency
/// denominator, while every other value identifies the first failed eligibility requirement.
struct GammaConvTrackInfo {
  enum Rejection : uint8_t {
    Accepted,
    WrongSign,
    NoITSorTPC,
    NoTPC,
    DirtyLabel
  };
  RecTrack track{};
  std::vector<int> pvIDs{}; // all reconstructed collisions compatible with this exact track ID
  uint16_t mcMask = 0;      // AOD-equivalent MC track-label mask
  uint8_t rejection = Accepted;
  bool inReferencePV = false;

  bool isEligible() const { return rejection == Accepted; }

  ClassDefNV(GammaConvTrackInfo, 2);
};

/// One row per selected truth photon conversion. The row contains the exact daughter labels,
/// every reconstructed representation of each daughter, and the exact raw-V0 matching result.
struct GammaConvInfo {
  enum Leg : uint8_t {
    Positron,
    Electron,
    NLegs
  };
  enum PairType : uint8_t {
    PairITSTPCITSTPC,
    PairITSTPCTPCOnly,
    PairTPCOnlyTPCOnly,
    PairContainsITSOnly,
    PairOther,
    NPairTypes
  };
  enum TerminalReason : uint8_t {
    V0Stored,
    BothLegsFoundNoV0,
    NeitherLegEligible,
    NoEligiblePositron,
    NoEligibleElectron
  };

  o2::MCCompLabel photonLabel{};
  std::array<o2::MCCompLabel, NLegs> daughterLabels{};
  std::array<std::vector<GammaConvTrackInfo>, NLegs> tracks{};
  std::array<float, 3> trueConversionXYZ{};
  float photonPt = 0.f;
  float photonEta = 0.f;
  float photonPhi = 0.f;
  int referencePV = -1;
  int referencePVNumContrib = -1;
  uint8_t terminalReason = NeitherLegEligible;
  bool bothLegsFoundAnywhere = false;
  bool bothLegsFoundInReferencePV = false;
  bool rawV0FoundUsingAnywhereEligiblePair = false;
  bool rawV0FoundUsingReferenceEligiblePair = false;
  // These flags are filled once per pair category, matching the O2Physics
  // EligiblePairType histograms without selecting a single track clone.
  std::array<bool, NPairTypes> referenceEligiblePairTypes{};
  std::array<bool, NPairTypes> rawV0FoundUsingReferenceEligiblePairTypes{};

  float getTrueConversionRadius() const { return std::hypot(trueConversionXYZ[0], trueConversionXYZ[1]); }

  ClassDefNV(GammaConvInfo, 2);
};

/// The three cuts SVertexer::processTPCTrack applies to TPC-only seeds under mTPCTrackPhotonTune,
/// recorded as the continuous variables behind them so that the thresholds can be studied rather
/// than only their outcome. One record per (TPC track, primary vertex): the SVertexer builds a
/// separate time-constrained clone of the track for every compatible vertex, and dz2Beam depends on
/// the vertex Z, so the same gid legitimately appears several times with different values.
struct TPCTuneInfo {
  enum Stage : uint8_t {
    RejMaxX,     // dropped by mTPCTrackMaxX before anything else is computed
    BothSides,   // has clusters on both sides, treated as constrained, the tune never applies
    RejTimeCorr, // the drift correction to the vertex time failed
    Evaluated    // all three cut variables are filled
  };
  o2::dataformats::VtxTrackIndex gid{};
  // A soft conversion electron loops in the TPC and is usually reconstructed as several track
  // segments, all carrying this same label. Efficiencies per track and per MC particle therefore
  // differ a lot, and only the latter is what a V0 sees: the prong survives if any segment does.
  o2::MCCompLabel mcLabel{};
  int vtxID = -1;
  float x = 0.f;         // track X, the mTPCTrackMaxX variable
  float zCorr = -999.f;  // Z after the drift correction to the vertex time
  float zPool = -999.f;  // Z of the corresponding SVertexer seed, -999 if the track was rejected
  float dz2Beam = -1.f;  // |x*tgl - zCorr + vtxZ|, cut against mTPCTrack2Beam       (dDPV)
  float cR = -1.f;       // distance of the helix centre from the beam line
  float rC = -1.f;       // helix radius
  float drd2Sq = 0.f;    // cR^2 - rC^2, cut against mTPCTrackXY2Radius^2            (dRD2)
  int16_t nClusters = -1;// cut against mTPCTrackMinNClusters                        (dCls)
  int mcPdg = 0;
  int mcMotherPdg = 0;
  // Truth kinematics of the track, needed to compare this sample like-for-like against the dec22
  // prongs: collectGammaConversions applies no daughter cuts, while the dec22 daughters must pass
  // acceptMCCharged and have a reconstructed partner, so the two samples are not the same tracks.
  float mcPt = -1.f;
  float mcR = -1.f; // production radius of the track, i.e. the conversion radius for a prong
  uint8_t stage = RejMaxX;
  bool accepted = false;  // the track survived all three cuts and stayed in the seeds pool
  bool isSignal = false;  // daughter of a photon passing the trmcconf.gamma* fiducial selection
  bool labelFake = false;
  // vtxID is the MC collision the track really comes from. A TPC-only track has a wide time
  // bracket and is tried against many vertices, so most records of a genuine conversion prong are
  // wrong-collision hypotheses which dz2Beam is meant to reject. Efficiencies must be quoted on
  // the correct-vertex records only, otherwise the association combinatorics look like a loss.
  bool isCorrectPV = false;

  // SVertexer takes sqrt(cR^2 - rC^2), which is NaN whenever the helix encloses the beam line.
  // NaN > threshold is false, so such tracks silently pass the cut. Kept signed here so that the
  // population can be counted instead of disappearing.
  bool isDrd2Undefined() const { return drd2Sq < 0.f; }
  float getDrd2() const { return drd2Sq >= 0.f ? std::sqrt(drd2Sq) : -1.f; }
  bool isEvaluated() const { return stage == Evaluated; }

  ClassDefNV(TPCTuneInfo, 4);
};

} // namespace o2::trackstudy
#endif
