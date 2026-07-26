// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CheckTPCPhotonTune.C
/// \brief Distributions of the TPC-only photon-tune cuts of SVertexer::processTPCTrack
///
/// Signal are TPC tracks of e+/e- from a fiducial-selected conversion, background is everything
/// else. For each of the three cuts the macro shows the underlying continuous variable, where the
/// current threshold sits, and what a different threshold would buy in signal efficiency versus
/// background rejection. It also quantifies the drd2 = sqrt(cR^2 - rC^2) pathology: the argument is
/// negative whenever the helix encloses the beam line, production takes the sqrt anyway and the
/// resulting NaN compares false, so those tracks pass a cut that was presumably meant to reject them.
///
/// Usage:
///   root -l -b -q 'CheckTPCPhotonTune.C+("trackMCStudy.root")'

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVirtualPad.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "GlobalTrackingStudy/TrackMCStudyTypes.h"
#endif

namespace
{
const char* stageNames[] = {"rej. maxX", "both sides", "rej. time corr", "evaluated"};
constexpr int NStages = sizeof(stageNames) / sizeof(char*);

void drawCut(TH1* h, double cut)
{
  gPad->Update();
  auto* l = new TLine(cut, gPad->GetUymin(), cut, gPad->GetUymax());
  l->SetLineColor(kBlack);
  l->SetLineStyle(2);
  l->SetLineWidth(2);
  l->Draw();
}

void styleSB(TH1* sig, TH1* bkg)
{
  sig->SetLineColor(kRed + 1);
  sig->SetLineWidth(2);
  bkg->SetLineColor(kAzure + 2);
  bkg->SetLineWidth(2);
  for (auto* h : {sig, bkg}) {
    if (h->Integral() > 0.) {
      h->Scale(1. / h->Integral());
    }
  }
  sig->SetMaximum(1.25 * std::max(sig->GetMaximum(), bkg->GetMaximum()));
}

// Signal efficiency and background rejection as a function of a "keep if value <= cut" threshold,
// scanned over the bins of the input distributions.
void makeScan(const TH1F* sig, const TH1F* bkg, TH1F*& eff, TH1F*& rej, const char* name, const char* xt)
{
  const int nb = sig->GetNbinsX();
  eff = new TH1F(Form("%sEff", name), Form(";%s threshold;fraction", xt), nb, sig->GetXaxis()->GetXmin(), sig->GetXaxis()->GetXmax());
  rej = new TH1F(Form("%sRej", name), Form(";%s threshold;fraction", xt), nb, sig->GetXaxis()->GetXmin(), sig->GetXaxis()->GetXmax());
  // include the overflow: a track above the last bin edge is kept only by an infinite threshold
  const double sTot = sig->Integral(0, nb + 1);
  const double bTot = bkg->Integral(0, nb + 1);
  double sCum = sig->GetBinContent(0), bCum = bkg->GetBinContent(0);
  for (int i = 1; i <= nb; i++) {
    sCum += sig->GetBinContent(i);
    bCum += bkg->GetBinContent(i);
    eff->SetBinContent(i, sTot > 0. ? sCum / sTot : 0.);
    rej->SetBinContent(i, bTot > 0. ? 1. - bCum / bTot : 0.);
  }
  eff->SetLineColor(kRed + 1);
  eff->SetLineWidth(2);
  rej->SetLineColor(kAzure + 2);
  rej->SetLineWidth(2);
  eff->SetMinimum(0.);
  eff->SetMaximum(1.05);
}
} // namespace

/// \param minMCPt truth pT floor applied to the signal, to reproduce the acceptMCCharged selection
///        the dec22 prongs went through (trmcconf.minPtMC). Set to 0 for the unbiased sample.
void CheckTPCPhotonTune(const std::string& inpFile = "trackMCStudy.root",
                        const std::string& outFile = "CheckTPCPhotonTune.root",
                        float minMCPt = 0.f)
{
  using namespace o2::trackstudy;
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(0);

  std::unique_ptr<TFile> fl{TFile::Open(inpFile.c_str())};
  if (!fl || fl->IsZombie()) {
    printf("ERROR: cannot open %s\n", inpFile.c_str());
    return;
  }
  auto* tree = (TTree*)fl->Get("tpcTune");
  if (!tree || !tree->GetBranch("trk")) {
    printf("ERROR: no tpcTune/trk output in %s, run with trmcconf.storeTPCPhotonTune=true\n", inpFile.c_str());
    return;
  }

  std::vector<TPCTuneInfo>* recs = nullptr;
  int minNCl = -1;
  float cutZ2Beam = -1.f, cutXY2Radius = -1.f, cutMaxX = -1.f;
  float bkgSampFrac = 1.f;
  tree->SetBranchAddress("trk", &recs);
  tree->SetBranchAddress("minNCl", &minNCl);
  tree->SetBranchAddress("cutZ2Beam", &cutZ2Beam);
  tree->SetBranchAddress("cutXY2Radius", &cutXY2Radius);
  tree->SetBranchAddress("cutMaxX", &cutMaxX);
  tree->SetBranchAddress("bkgSampFrac", &bkgSampFrac);

  auto* hStage = new TH2F("hStage", "TPC seed attrition;;", NStages, -0.5, NStages - 0.5, 2, -0.5, 1.5);
  for (int i = 0; i < NStages; i++) {
    hStage->GetXaxis()->SetBinLabel(i + 1, stageNames[i]);
  }
  hStage->GetYaxis()->SetBinLabel(1, "background");
  hStage->GetYaxis()->SetBinLabel(2, "signal");

  auto* hNClS = new TH1F("hNClSig", "N TPC clusters;N clusters;norm.", 160, 0., 160.);
  auto* hNClB = new TH1F("hNClBkg", "N TPC clusters;N clusters;norm.", 160, 0., 160.);
  auto* hDzS = new TH1F("hDz2BeamSigCorrPV", "extrapolation to the beam line;|x#upointtg#lambda - z + z_{PV}| (cm);norm.", 120, 0., 120.);
  auto* hDzW = new TH1F("hDz2BeamSigWrongPV", "extrapolation to the beam line;|x#upointtg#lambda - z + z_{PV}| (cm);norm.", 120, 0., 120.);
  auto* hDzB = new TH1F("hDz2BeamBkg", "extrapolation to the beam line;|x#upointtg#lambda - z + z_{PV}| (cm);norm.", 120, 0., 120.);
  auto* hDzBc = new TH1F("hDz2BeamBkgCorrPV", ";|x#upointtg#lambda - z + z_{PV}| (cm);norm.", 120, 0., 120.);
  auto* hDzBw = new TH1F("hDz2BeamBkgWrongPV", ";|x#upointtg#lambda - z + z_{PV}| (cm);norm.", 120, 0., 120.);
  auto* hDrS = new TH1F("hDrd2Sig", "tangent length to the helix;#sqrt{c_{R}^{2}-r_{C}^{2}} (cm);norm.", 150, 0., 300.);
  auto* hDrB = new TH1F("hDrd2Bkg", "tangent length to the helix;#sqrt{c_{R}^{2}-r_{C}^{2}} (cm);norm.", 150, 0., 300.);
  auto* hUndef = new TH1F("hDrd2Undefined", "c_{R}^{2}-r_{C}^{2} < 0, sqrt is NaN and the cut passes;;tracks", 2, -0.5, 1.5);
  hUndef->GetXaxis()->SetBinLabel(1, "background");
  hUndef->GetXaxis()->SetBinLabel(2, "signal");
  // seeding efficiency of the correct-collision records against truth kinematics: this is what has
  // to be compared with the dec22 prongs, which live behind a minPtMC floor and a reconstructed partner
  auto* hPtAll = new TH1F("hPtAllSig", "signal, correct PV;true p_{T} (GeV/c);prongs", 50, 0., 1.);
  auto* hPtAcc = new TH1F("hPtAccSig", "signal, correct PV, kept;true p_{T} (GeV/c);prongs", 50, 0., 1.);
  auto* hRAll = new TH1F("hRAllSig", "signal, correct PV;true R_{conv} (cm);prongs", 45, 0., 90.);
  auto* hRAcc = new TH1F("hRAccSig", "signal, correct PV, kept;true R_{conv} (cm);prongs", 45, 0., 90.);
  auto* hCRvsRC = new TH2F("hCRvsRC", "signal;r_{C} (cm);c_{R} (cm)", 100, 0., 400., 100, 0., 400.);
  auto* hDzVsDr = new TH2F("hDzVsDr", "signal;#sqrt{c_{R}^{2}-r_{C}^{2}} (cm);|#Deltaz to beam| (cm)",
                           100, 0., 300., 100, 0., 120.);

  long nRec = 0, nSig = 0, nEval = 0, nEvalSig = 0;
  long nUndefS = 0, nUndefB = 0, nAccB = 0;
  long nUndefAccS = 0, nUndefAccB = 0; // undefined drd2 among tracks the tune actually accepted
  // records split by whether the vertex under test is the track's true collision
  long nCorr = 0, nAccCorr = 0, nWrong = 0, nAccWrong = 0;
  long nCorrCut = 0, nAccCorrCut = 0; // same, restricted to mcPt >= minMCPt
  long nBkgCorr = 0, nAccBkgCorr = 0; // background sitting on its own collision
  // per MC particle rather than per reconstructed track: this is the prong-level quantity that a V0
  // actually sees, and the one directly comparable with eps1 from CheckSVGamma
  long nMCPart = 0, nMCPartAnyAcc = 0, nMCPartCorrAcc = 0, nMCPartTracks = 0;
  auto* hNSeg = new TH1F("hNSegPerProng", "TPC segments per conversion prong;segments;prongs", 12, 0.5, 12.5);
  auto* hNSegAcc = new TH1F("hNSegAccPerProng", "accepted segments per prong;accepted segments;prongs", 12, -0.5, 11.5);
  // per distinct track, so that the vertex-independent variables are not weighted by how many
  // vertices the track happens to be compatible with
  long nTrkSig = 0, nTrkSigWithCorr = 0, nTrkSigCorrAcc = 0, nTrkSigAnyAcc = 0;

  struct TrackAgg {
    const TPCTuneInfo* rep = nullptr; // representative record, the correct-PV one when there is one
    bool isSignal = false;
    bool hasCorr = false;
    bool corrAcc = false;
    bool anyAcc = false;
  };

  const auto nEnt = tree->GetEntries();
  for (Long64_t ie = 0; ie < nEnt; ie++) {
    tree->GetEntry(ie);
    if (!recs) {
      continue;
    }
    std::map<uint32_t, TrackAgg> perTrack; // gids are unique within a TF
    struct PartAgg {
      int nSeg = 0;
      int nSegAcc = 0;
      bool anyAcc = false;
      bool corrAcc = false;
    };
    std::map<uint64_t, PartAgg> perPart; // keyed by MC label, signal only
    for (const auto& r : *recs) {
      nRec++;
      const bool s = r.isSignal;
      nSig += s;
      hStage->Fill(r.stage, s ? 1 : 0);
      if (!r.isEvaluated()) {
        continue;
      }
      nEval++;
      nEvalSig += s;

      // dz2Beam is the only vertex-dependent variable, so it is filled per record and split
      // three ways: a wrong-collision record of a genuine prong is something the cut should reject
      if (!s) {
        hDzB->Fill(r.dz2Beam);
        nAccB += r.accepted;
        // background splits the same way: a non-conversion track sitting on its own collision also
        // extrapolates back to it, which is where the peak at 0 in the background comes from
        (r.isCorrectPV ? hDzBc : hDzBw)->Fill(r.dz2Beam);
        if (r.isCorrectPV) {
          nBkgCorr++;
          nAccBkgCorr += r.accepted;
        }
      } else if (r.isCorrectPV) {
        hDzS->Fill(r.dz2Beam);
        nCorr++;
        nAccCorr += r.accepted;
        hPtAll->Fill(r.mcPt);
        hRAll->Fill(r.mcR);
        if (r.accepted) {
          hPtAcc->Fill(r.mcPt);
          hRAcc->Fill(r.mcR);
        }
        if (r.mcPt >= minMCPt) {
          nCorrCut++;
          nAccCorrCut += r.accepted;
        }
      } else {
        hDzW->Fill(r.dz2Beam);
        nWrong++;
        nAccWrong += r.accepted;
      }

      auto& agg = perTrack[uint32_t(r.gid)];
      agg.isSignal = s;
      if (s && r.mcLabel.isValid()) {
        auto& pa = perPart[r.mcLabel.getRawValue()];
        pa.anyAcc = pa.anyAcc || r.accepted;
        if (r.isCorrectPV) {
          pa.corrAcc = pa.corrAcc || r.accepted;
        }
      }
      agg.anyAcc = agg.anyAcc || r.accepted;
      if (r.isCorrectPV) {
        agg.hasCorr = true;
        agg.corrAcc = r.accepted;
        agg.rep = &r; // prefer the correct-collision record
      } else if (!agg.rep) {
        agg.rep = &r;
      }
    }

    for (const auto& [gid, agg] : perTrack) {
      if (!agg.rep) {
        continue;
      }
      const auto& r = *agg.rep;
      const bool s = agg.isSignal;
      (s ? hNClS : hNClB)->Fill(r.nClusters);
      if (r.isDrd2Undefined()) {
        hUndef->Fill(s ? 1 : 0);
        (s ? nUndefS : nUndefB)++;
        if (r.accepted) {
          (s ? nUndefAccS : nUndefAccB)++;
        }
      } else {
        (s ? hDrS : hDrB)->Fill(r.getDrd2());
      }
      if (s) {
        nTrkSig++;
        nTrkSigWithCorr += agg.hasCorr;
        nTrkSigCorrAcc += agg.corrAcc;
        nTrkSigAnyAcc += agg.anyAcc;
        hCRvsRC->Fill(r.rC, r.cR);
        hDzVsDr->Fill(r.getDrd2(), r.dz2Beam);
        // one distinct gid is one TPC segment of this MC particle
        if (r.mcLabel.isValid()) {
          auto& pa = perPart[r.mcLabel.getRawValue()];
          pa.nSeg++;
          pa.nSegAcc += agg.anyAcc;
        }
      }
    }

    for (const auto& [lbl, pa] : perPart) {
      if (!pa.nSeg) {
        continue;
      }
      nMCPart++;
      nMCPartTracks += pa.nSeg;
      nMCPartAnyAcc += pa.anyAcc;
      nMCPartCorrAcc += pa.corrAcc;
      hNSeg->Fill(std::min(pa.nSeg, 12));
      hNSegAcc->Fill(std::min(pa.nSegAcc, 11));
    }
  }
  const long nAccS = nAccCorr; // "signal accepted" means the correct-collision record survived

  auto pct = [](double a, double b) { return b > 0. ? 100. * a / b : 0.; };
  // background was sampled, signal was not: every count below has to be put back on a common
  // footing before any absolute rate or signal fraction means anything
  const double w = bkgSampFrac > 0.f ? 1. / bkgSampFrac : 1.;
  const long nEvalBkg = nEval - nEvalSig;
  const double nEvalBkgT = nEvalBkg * w, nAccBkgT = nAccB * w, nUndefBkgT = nUndefB * w;

  printf("\n=== TPC photon tune, %s ===\n", inpFile.c_str());
  printf("thresholds: minNClusters=%d  z2Beam=%.2f cm  xy2Radius=%.2f cm  maxX=%.2f\n",
         minNCl, cutZ2Beam, cutXY2Radius, cutMaxX);
  if (minNCl < 0) {
    printf("  NOTE: minNClusters < 0, so the cluster cut can never fire (counts are >= 0).\n");
  }
  if (cutMaxX < 0.f) {
    printf("  NOTE: maxX < 0, so the X cut is disabled. Only z2Beam and xy2Radius are live.\n");
  }
  printf("\nbackground sampling fraction %.4f, background counts below are scaled by %.1f\n", bkgSampFrac, w);
  printf("stored records %ld (signal %ld, background %ld)\n", nRec, nSig, nRec - nSig);
  printf("reaching the tune: signal %ld, background %ld (%.0f before sampling)\n",
         nEvalSig, nEvalBkg, nEvalBkgT);
  printf("  signal fraction of the true TPC seed sample: %.2f%%\n",
         pct(nEvalSig, nEvalSig + nEvalBkgT));

  printf("\nsignal records split by collision hypothesis:\n");
  printf("  correct collision %8ld  accepted %8ld  (%5.2f%%)   <- the real efficiency\n",
         nCorr, nAccCorr, pct(nAccCorr, nCorr));
  printf("  wrong collision   %8ld  accepted %8ld  (%5.2f%%)   <- should be rejected, not a loss\n",
         nWrong, nAccWrong, pct(nAccWrong, nWrong));
  printf("  %.1f wrong-collision hypotheses per genuine prong, so a per-record efficiency of %.2f%%\n"
         "  is mostly association combinatorics rather than lost signal\n",
         nCorr > 0 ? double(nWrong) / nCorr : 0., pct(nAccCorr + nAccWrong, nCorr + nWrong));

  if (minMCPt > 0.f) {
    printf("\nrestricted to true p_{T} >= %.3f GeV/c, i.e. the acceptMCCharged floor the dec22\n"
           "prongs live behind (trmcconf.minPtMC):\n", minMCPt);
    printf("  correct collision %8ld  accepted %8ld  (%5.2f%%)\n",
           nCorrCut, nAccCorrCut, pct(nAccCorrCut, nCorrCut));
    printf("  compare with eps1 from CheckSVGamma. A remaining gap is not the pT floor but the\n"
           "  other dec22 requirements: both daughters trackable and both reconstructed.\n");
  } else {
    printf("\nNOTE: no pT floor applied. The dec22 prongs pass acceptMCCharged (minPtMC) and need a\n"
           "reconstructed partner, so eps1 there is measured on an easier sample than this one.\n"
           "Rerun with the third argument set to trmcconf.minPtMC to compare like for like.\n");
  }

  printf("\nper MC particle, i.e. per conversion prong (the quantity a V0 actually sees):\n");
  printf("  distinct conversion prongs       %8ld\n", nMCPart);
  printf("  TPC segments per prong           %8.2f\n", nMCPart ? double(nMCPartTracks) / nMCPart : 0.);
  printf("  at least one segment kept        %8ld  (%5.2f%%)  <- compare with eps1 in CheckSVGamma\n",
         nMCPartAnyAcc, pct(nMCPartAnyAcc, nMCPart));
  printf("  correct-PV segment kept          %8ld  (%5.2f%%)\n", nMCPartCorrAcc, pct(nMCPartCorrAcc, nMCPart));
  printf("  A prong survives if any of its segments does, so the per-track efficiency above is not\n"
         "  the conversion efficiency: with n segments the two differ by roughly 1-(1-eps)^n.\n");

  printf("\nper distinct signal track (%ld tracks):\n", nTrkSig);
  printf("  has a correct-collision record   %8ld  (%5.2f%%)\n", nTrkSigWithCorr, pct(nTrkSigWithCorr, nTrkSig));
  printf("  correct-collision record kept    %8ld  (%5.2f%% of those)\n",
         nTrkSigCorrAcc, pct(nTrkSigCorrAcc, nTrkSigWithCorr));
  printf("  kept under at least one vertex   %8ld  (%5.2f%%)\n", nTrkSigAnyAcc, pct(nTrkSigAnyAcc, nTrkSig));

  printf("\nacceptance of the tune:\n");
  printf("  signal (correct PV) %8ld / %-8ld (%5.2f%%)\n", nAccCorr, nCorr, pct(nAccCorr, nCorr));
  printf("  background          %8.0f / %-8.0f (%5.2f%%)\n", nAccBkgT, nEvalBkgT, pct(nAccBkgT, nEvalBkgT));
  printf("  background rejection factor %.1fx, signal/background improvement %.2fx\n",
         nAccBkgT > 0. ? nEvalBkgT / nAccBkgT : 0.,
         pct(nAccBkgT, nEvalBkgT) > 0. ? pct(nAccCorr, nCorr) / pct(nAccBkgT, nEvalBkgT) : 0.);
  printf("  background on its OWN collision %8.0f accepted %8.0f (%5.2f%%)\n",
         nBkgCorr * w, nAccBkgCorr * w, pct(nAccBkgCorr, nBkgCorr));
  printf("  Most of the background rejection is wrong-collision hypotheses. Against background that\n"
         "  really belongs to the vertex under test the cut has far less power, compare the two\n"
         "  numbers above with the signal efficiency.\n");

  printf("\ncR^2-rC^2 < 0: the helix encloses the beam line, sqrt is NaN and the cut passes\n");
  printf("  signal     %8ld  (%5.2f%% of evaluated signal)\n", nUndefS, pct(nUndefS, nEvalSig));
  printf("  background %8.0f  (%5.2f%% of evaluated background)\n", nUndefBkgT, pct(nUndefBkgT, nEvalBkgT));
  printf("  of which still accepted by the other cuts, i.e. the real cost of fixing the sqrt:\n");
  printf("    signal     %8ld  (%5.2f%% of accepted signal would be lost)\n",
         nUndefAccS, pct(nUndefAccS, nAccS));
  printf("    background %8.0f  (%5.2f%% of accepted background would be removed)\n",
         nUndefAccB * w, pct(nUndefAccB * w, nAccBkgT));

  TH1F *effDz = nullptr, *rejDz = nullptr, *effDr = nullptr, *rejDr = nullptr;
  // scan before normalising, the scan needs the raw counts
  makeScan(hDzS, hDzB, effDz, rejDz, "hDz2Beam", "|#Deltaz to beam|");
  makeScan(hDrS, hDrB, effDr, rejDr, "hDrd2", "#sqrt{c_{R}^{2}-r_{C}^{2}}");
  styleSB(hNClS, hNClB);
  styleSB(hDzS, hDzB);
  styleSB(hDrS, hDrB);
  hDzW->SetLineColor(kOrange + 7);
  hDzW->SetLineWidth(2);
  hDzW->SetLineStyle(2);
  if (hDzW->Integral() > 0.) {
    hDzW->Scale(1. / hDzW->Integral());
  }

  auto* leg = new TLegend(0.55, 0.70, 0.95, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hDzS, "conversion e^{#pm}, correct PV", "l");
  leg->AddEntry(hDzW, "conversion e^{#pm}, wrong PV", "l");
  leg->AddEntry(hDzBc, "background, correct PV", "l");
  leg->AddEntry(hDzBw, "background, wrong PV", "l");

  auto* cnv = new TCanvas("cTPCTune", "TPC-only photon tune", 1500, 900);
  cnv->Divide(3, 2);
  cnv->cd(1);
  hNClS->Draw("hist");
  hNClB->Draw("hist same");
  leg->Draw();
  if (minNCl >= 0) {
    drawCut(hNClS, minNCl);
  }
  cnv->cd(2);
  gPad->SetLogy();
  for (auto* h : {hDzBc, hDzBw}) {
    if (h->Integral() > 0.) {
      h->Scale(1. / h->Integral());
    }
    h->SetLineColor(kAzure + 2);
    h->SetLineWidth(2);
  }
  hDzBw->SetLineStyle(2);
  hDzS->SetMaximum(3. * std::max({hDzS->GetMaximum(), hDzW->GetMaximum(), hDzBc->GetMaximum()}));
  hDzS->Draw("hist");
  hDzW->Draw("hist same");
  hDzBc->Draw("hist same");
  hDzBw->Draw("hist same");
  drawCut(hDzS, cutZ2Beam);
  leg->Draw();
  cnv->cd(3);
  gPad->SetLogy();
  hDrS->Draw("hist");
  hDrB->Draw("hist same");
  drawCut(hDrS, cutXY2Radius);
  cnv->cd(4);
  hUndef->Draw("hist");
  // legend for the scans: without it the two curves are easy to read backwards
  auto* legScan = new TLegend(0.30, 0.18, 0.92, 0.34);
  legScan->SetBorderSize(0);
  legScan->SetFillStyle(0);
  legScan->SetTextSize(0.033);
  legScan->AddEntry(effDz, "signal kept (efficiency), rises as the cut loosens", "l");
  legScan->AddEntry(rejDz, "background removed (rejection), falls as the cut loosens", "l");
  cnv->cd(5);
  effDz->Draw("hist");
  rejDz->Draw("hist same");
  drawCut(effDz, cutZ2Beam);
  legScan->Draw();
  cnv->cd(6);
  effDr->Draw("hist");
  rejDr->Draw("hist same");
  drawCut(effDr, cutXY2Radius);
  legScan->Draw();

  auto* hEffPt = (TH1F*)hPtAcc->Clone("hEffPt");
  hEffPt->SetTitle("signal, correct PV;true p_{T} (GeV/c);seeding efficiency");
  hEffPt->Divide(hPtAcc, hPtAll, 1., 1., "B");
  auto* hEffR = (TH1F*)hRAcc->Clone("hEffRconv");
  hEffR->SetTitle("signal, correct PV;true R_{conv} (cm);seeding efficiency");
  hEffR->Divide(hRAcc, hRAll, 1., 1., "B");
  for (auto* h : {hEffPt, hEffR}) {
    h->SetLineColor(kRed + 1);
    h->SetLineWidth(2);
    h->SetMinimum(0.);
    h->SetMaximum(1.05);
  }

  auto* cnv2 = new TCanvas("cTPCTuneCorr", "TPC-only photon tune, correlations", 1600, 900);
  cnv2->Divide(3, 2);
  cnv2->cd(1);
  hStage->Draw("colz text");
  cnv2->cd(2);
  hCRvsRC->Draw("colz");
  cnv2->cd(3);
  hDzVsDr->Draw("colz");
  // where the sample difference against dec22 lives: if the efficiency collapses at low pT, the
  // gap is the minPtMC floor rather than anything about the cut itself
  cnv2->cd(4);
  hEffPt->Draw("e");
  if (minMCPt > 0.f) {
    drawCut(hEffPt, minMCPt);
  }
  cnv2->cd(5);
  hEffR->Draw("e");
  // the segment multiplicity is what turns a ~50% per-track efficiency into a ~90% per-prong one
  cnv2->cd(6);
  gPad->SetLogy();
  hNSeg->SetLineColor(kAzure + 2);
  hNSeg->SetLineWidth(2);
  hNSeg->Draw("hist");
  hNSegAcc->SetLineColor(kRed + 1);
  hNSegAcc->SetLineWidth(2);
  hNSegAcc->Draw("hist same");

  std::unique_ptr<TFile> flOut{TFile::Open(outFile.c_str(), "recreate")};
  for (auto* h : {hNClS, hNClB, hDzS, hDzW, hDzB, hDrS, hDrB, hUndef, effDz, rejDz, effDr, rejDr}) {
    h->Write();
  }
  hStage->Write();
  hCRvsRC->Write();
  hDzVsDr->Write();
  hNSeg->Write();
  hNSegAcc->Write();
  hPtAll->Write();
  hPtAcc->Write();
  hRAll->Write();
  hRAcc->Write();
  hEffPt->Write();
  hEffR->Write();
  cnv->Write();
  cnv2->Write();
  const auto base = outFile.substr(0, outFile.find_last_of('.'));
  cnv->SaveAs((base + ".png").c_str());
  cnv2->SaveAs((base + "_corr.png").c_str());
  printf("\nWritten %s\n", outFile.c_str());
}
