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

/// \file CheckSVGamma.C
/// \brief Why are reconstructed prongs of gamma conversions rejected before the SVertexer seeds pool
///
/// Reads the dec<PDG> tree written by o2-trackMC-study-workflow and, for every watched decay whose
/// prongs were reconstructed as tracks, asks whether those tracks made it into the SVertexer seeds
/// pool and, if not, which check in SVertexer::buildT2V dropped them. A prong which never becomes a
/// seed can never form a V0, whatever the pair cuts do, so this stage has to be understood first.
///
/// Usage:
///   root -l -b -q 'CheckSVGamma.C+("trackMCStudy.root", 22)'

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVirtualPad.h>

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "GlobalTrackingStudy/TrackMCStudyTypes.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#endif

namespace
{
// Mirrors o2::vertexing::SVertexer::SeedRej. Kept local on purpose: including SVertexer.h drags in
// the GPU and strangeness-tracking headers, which is more than a plotting macro should need.
// The first bin is not a rejection: buildT2V only records a reason for tracks it actually examined,
// so a prong with no reason was never offered to the SVertexer at all, i.e. it is not attached to
// any primary vertex and therefore never appears in the PV-matched track list buildT2V loops over.
const char* seedRejNames[] = {"not in PV list",  // SeedRejNone
                              "source not loaded",
                              "TPC excluded",
                              "TPC maxX",
                              "TPC time corr",
                              "TPC photon tune",
                              "acceptTrack (DCA/PV contrib)",
                              "short ITS-only"};
constexpr int NSeedRej = sizeof(seedRejNames) / sizeof(char*);

// Mirrors o2::trackstudy::SVCheck::Stage, shifted by one so that NotChecked = -1 lands in bin 0.
const char* stageNames[] = {"not checked", "no prongs", "not seeded", "same charge",
                            "no bracket overlap", "cut rejected", "found"};
constexpr int NStages = sizeof(stageNames) / sizeof(char*);

// Mirrors o2::vertexing::SVertexer::V0Rej. The last two are stages rather than cuts: with cascades
// or 3-body decays enabled a pair failing an earlier cut is demoted instead of dropped, and if it
// dies there the original cut is no longer known.
const char* v0RejNames[] = {"accepted", "tgl diff", "D2R", "DR", "DCAFitter", "minR2 to mean vtx",
                            "causality", "propagation", "pt2", "tgl V0", "hypothesis",
                            "DCAXY casc V0", "DCAXY/cosPAXY", "cosPA", "after 3body", "not cascade"};
constexpr int NV0Rej = sizeof(v0RejNames) / sizeof(char*);

// one colour per SeedRej reason, used consistently across every plot in this macro
const int seedRejColors[] = {kGray + 1, kOrange + 7, kRed + 1, kMagenta + 1,
                             kViolet + 1, kBlue + 1, kAzure + 7, kGreen + 2};

void setLabels(TAxis* ax, const char** names, int n)
{
  for (int i = 0; i < n; i++) {
    ax->SetBinLabel(i + 1, names[i]);
  }
}

// Minimum number of lost prongs in a (pT, R) cell before naming a dominant cut for it, so that the
// map shows structure rather than the argmax of one or two entries.
constexpr double MinLossForDominant = 5.;
} // namespace

void CheckSVGamma(const std::string& inpFile = "trackMCStudy.root",
                  int decayPDG = 22,
                  const std::string& outFile = "CheckSVGamma.root")
{
  using GTrackID = o2::dataformats::GlobalTrackID;
  constexpr int NSrc = GTrackID::NSources;

  TH1::AddDirectory(kFALSE); // keep the histograms alive when the input file goes away

  std::unique_ptr<TFile> fl{TFile::Open(inpFile.c_str())};
  if (!fl || fl->IsZombie()) {
    printf("ERROR: cannot open %s\n", inpFile.c_str());
    return;
  }
  const std::string trName = "dec" + std::to_string(decayPDG);
  auto* tree = (TTree*)fl->Get(trName.c_str());
  if (!tree) {
    printf("ERROR: no tree %s in %s, was the workflow run with trmcconf.decayPDG=[%d] ?\n",
           trName.c_str(), inpFile.c_str(), decayPDG);
    return;
  }
  if (!tree->GetBranch("svCheck")) {
    printf("ERROR: tree %s has no svCheck branch, the workflow ran with trmcconf.checkSVertexerCuts=false\n",
           trName.c_str());
    return;
  }

  o2::trackstudy::SVCheck* svCheck = nullptr;
  std::vector<o2::trackstudy::TrackFamily>* prod = nullptr;
  tree->SetBranchAddress("svCheck", &svCheck);
  tree->SetBranchAddress("prod", &prod);

  // decay level, for the denominators
  auto* hStage = new TH1F("hStage", "Decay stage;;decays", NStages, -0.5, NStages - 0.5);
  setLabels(hStage->GetXaxis(), stageNames, NStages);

  // prong level, the actual question
  auto* hSeedRej = new TH1F("hSeedRej", "Reconstructed prongs missing from the seeds pool;;prongs",
                            NSeedRej, -0.5, NSeedRej - 0.5);
  setLabels(hSeedRej->GetXaxis(), seedRejNames, NSeedRej);

  auto* hSeedRejVsR = new TH2F("hSeedRejVsR", "Missing prongs vs MC conversion radius;R_{conv} (cm);",
                               100, 0., 100., NSeedRej, -0.5, NSeedRej - 0.5);
  setLabels(hSeedRejVsR->GetYaxis(), seedRejNames, NSeedRej);

  auto* hSeedRejVsPt = new TH2F("hSeedRejVsPt", "Missing prongs vs MC p_{T};p_{T} (GeV/c);",
                                100, 0., 2., NSeedRej, -0.5, NSeedRej - 0.5);
  setLabels(hSeedRejVsPt->GetYaxis(), seedRejNames, NSeedRej);

  auto* hSeedRejVsSrc = new TH2F("hSeedRejVsSrc", "Missing prongs vs track source;;",
                                 NSrc, -0.5, NSrc - 0.5, NSeedRej, -0.5, NSeedRej - 0.5);
  for (int i = 0; i < NSrc; i++) {
    hSeedRejVsSrc->GetXaxis()->SetBinLabel(i + 1, GTrackID::getSourceName(i).c_str());
  }
  setLabels(hSeedRejVsSrc->GetYaxis(), seedRejNames, NSeedRej);

  // reference distributions of the prongs which did become seeds, so the shapes above can be read
  // as an inefficiency rather than just as the shape of the parent sample
  auto* hSeededR = new TH1F("hSeededR", "Seeded prongs;R_{conv} (cm);prongs", 100, 0., 100.);
  auto* hSeededPt = new TH1F("hSeededPt", "Seeded prongs;p_{T} (GeV/c);prongs", 100, 0., 2.);
  auto* hAllR = new TH1F("hAllR", "Reconstructed prongs;R_{conv} (cm);prongs", 100, 0., 100.);
  auto* hAllPt = new TH1F("hAllPt", "Reconstructed prongs;p_{T} (GeV/c);prongs", 100, 0., 2.);

  // ---- findable gammas: both prongs have a reconstructed track, so the only thing that can still
  // ---- stop the V0 from being formed is the seeding. This is the sample the maps below are built on.
  const int nbPt = 20, nbR = 25;
  const double ptMax = 2., rMax = 100.;
  auto* hFndAllPtR = new TH2F("hFndAllPtR", "Findable-gamma prongs;p_{T} (GeV/c);R_{conv} (cm)",
                              nbPt, 0., ptMax, nbR, 0., rMax);
  auto* hFndAllPt = new TH1F("hFndAllPt", "Findable-gamma prongs;p_{T} (GeV/c);prongs", nbPt, 0., ptMax);
  auto* hFndAllR = new TH1F("hFndAllR", "Findable-gamma prongs;R_{conv} (cm);prongs", nbR, 0., rMax);

  std::vector<TH2F*> hFndRejPtR(NSeedRej, nullptr);
  std::vector<TH1F*> hFndRejPt(NSeedRej, nullptr), hFndRejR(NSeedRej, nullptr);
  for (int i = 0; i < NSeedRej; i++) {
    hFndRejPtR[i] = new TH2F(Form("hFndRejPtR_%d", i), Form("%s;p_{T} (GeV/c);R_{conv} (cm)", seedRejNames[i]),
                             nbPt, 0., ptMax, nbR, 0., rMax);
    hFndRejPt[i] = new TH1F(Form("hFndRejPt_%d", i), Form("%s;p_{T} (GeV/c);lost fraction", seedRejNames[i]),
                            nbPt, 0., ptMax);
    hFndRejR[i] = new TH1F(Form("hFndRejR_%d", i), Form("%s;R_{conv} (cm);lost fraction", seedRejNames[i]),
                           nbR, 0., rMax);
  }
  long nFindable = 0, nFndProngs = 0, nFndMissing = 0;

  // ---- joint vs single-prong seeding ---------------------------------------------------------
  // The pair-level efficiency is not the square of the per-prong one: the two prongs share a
  // conversion point and a collision, and dz2Beam depends on the same vertex Z for both, so they
  // succeed and fail together. Measure the joint probability instead of assuming independence.
  // Split by prong composition because the TPC-only photon tune only acts on TPC standalone prongs.
  enum Comp { BothTPC,
              OneTPC,
              NoTPC,
              NComp };
  const char* compNames[NComp] = {"both prongs TPC-only", "one prong TPC-only", "no TPC-only prong"};
  std::array<long, NComp> nCompFnd{}, nCompProngSeeded{}, nCompBothSeeded{};
  // which cut in checkV0 rejected the pair, for the decays that reached the pairing stage
  auto* hRejV0 = new TH1F("hRejV0", "Cut rejecting the pair in checkV0;;gammas", NV0Rej, -0.5, NV0Rej - 0.5);
  setLabels(hRejV0->GetXaxis(), v0RejNames, NV0Rej);
  auto* hRejV0VsR = new TH2F("hRejV0VsR", "Cut rejecting the pair;R_{conv} (cm);",
                             nbR, 0., rMax, NV0Rej, -0.5, NV0Rej - 0.5);
  setLabels(hRejV0VsR->GetYaxis(), v0RejNames, NV0Rej);

  auto* hFndDecR = new TH1F("hFndDecR", "findable gammas;R_{conv} (cm);gammas", nbR, 0., rMax);
  auto* hBothSeededR = new TH1F("hBothSeededR", "both prongs seeded;R_{conv} (cm);gammas", nbR, 0., rMax);
  auto* hProngSeededR = new TH1F("hProngSeededR", "seeded prongs;R_{conv} (cm);prongs", nbR, 0., rMax);

  long nDec = 0, nProngsRec = 0, nProngsSeeded = 0, nProngsMissing = 0, nInconsistent = 0;
  const long nEnt = tree->GetEntries();
  for (long ient = 0; ient < nEnt; ient++) {
    tree->GetEntry(ient);
    nDec++;
    hStage->Fill(svCheck->stage + 1); // NotChecked = -1 -> bin 0
    if (!svCheck->replayConsistent) {
      nInconsistent++;
    }
    // prod[ip] and svCheck->prongs[ip] are filled from the same daughterFirst..daughterLast range,
    // so they are index aligned; bail out rather than mismatch them if that ever changes
    if (!prod || prod->size() != svCheck->prongs.size()) {
      continue;
    }
    // a gamma is findable once both prongs exist as tracks: from here on only the seeding and the
    // pair cuts can lose it, so this is the right denominator for a seeding inefficiency
    const bool findable = svCheck->bothProngsReconstructed();
    if (findable) {
      nFindable++;
      // gid source is filled for every reconstructed prong, seeded or not, so the TPC-only
      // classification survives the prongs that never entered the pool
      int nTPCOnly = 0, nSeeded = 0;
      for (const auto& pr : svCheck->prongs) {
        nTPCOnly += (pr.gid.getSource() == GTrackID::TPC);
        nSeeded += pr.isSeeded();
      }
      const int comp = nTPCOnly == 2 ? BothTPC : (nTPCOnly == 1 ? OneTPC : NoTPC);
      nCompFnd[comp]++;
      nCompProngSeeded[comp] += nSeeded;
      nCompBothSeeded[comp] += (nSeeded == 2);
      // both prongs are born at the same point, so either one gives the conversion radius
      const auto& mcTr0 = (*prod)[0].mcTrackInfo.track;
      const float rDec = std::hypot(mcTr0.getX(), mcTr0.getY());
      hFndDecR->Fill(rDec);
      hProngSeededR->Fill(rDec, nSeeded);
      if (nSeeded == 2) {
        hBothSeededR->Fill(rDec);
      }
      using SVC = o2::trackstudy::SVCheck;
      if (svCheck->stage == SVC::CutRejected || svCheck->stage == SVC::Found) {
        hRejV0->Fill(svCheck->rejV0);
        hRejV0VsR->Fill(rDec, svCheck->rejV0);
      }
    }
    for (size_t ip = 0; ip < svCheck->prongs.size(); ip++) {
      const auto& pr = svCheck->prongs[ip];
      if (!pr.isReconstructed()) {
        continue; // nothing was reconstructed for this MC particle, a tracking loss not a seeding one
      }
      nProngsRec++;
      const auto& mcTr = (*prod)[ip].mcTrackInfo.track;
      // daughters are born at the conversion point, and R is invariant under the rotation into the
      // tracking frame, so this is the MC conversion radius
      const float rConv = std::hypot(mcTr.getX(), mcTr.getY());
      const float pt = mcTr.getPt();
      hAllR->Fill(rConv);
      hAllPt->Fill(pt);
      if (findable) {
        nFndProngs++;
        hFndAllPtR->Fill(pt, rConv);
        hFndAllPt->Fill(pt);
        hFndAllR->Fill(rConv);
      }
      if (pr.isSeeded()) {
        nProngsSeeded++;
        hSeededR->Fill(rConv);
        hSeededPt->Fill(pt);
        continue;
      }
      nProngsMissing++;
      const int rej = pr.seedRej < NSeedRej ? pr.seedRej : 0;
      hSeedRej->Fill(rej);
      hSeedRejVsR->Fill(rConv, rej);
      hSeedRejVsPt->Fill(pt, rej);
      hSeedRejVsSrc->Fill(pr.gid.getSource(), rej);
      if (findable) {
        nFndMissing++;
        hFndRejPtR[rej]->Fill(pt, rConv);
        hFndRejPt[rej]->Fill(pt);
        hFndRejR[rej]->Fill(rConv);
      }
    }
  }

  auto pct = [](long a, long b) { return b ? 100. * a / b : 0.; };
  printf("\n=== %s: %ld decays from %s ===\n", trName.c_str(), nDec, inpFile.c_str());
  printf("\nDecay stages:\n");
  for (int i = 0; i < NStages; i++) {
    const auto n = (long)hStage->GetBinContent(i + 1);
    if (n) {
      printf("  %-22s %8ld  (%5.2f%%)\n", stageNames[i], n, pct(n, nDec));
    }
  }
  if (nInconsistent) {
    printf("\n  NOTE: %ld decays (%.2f%%) have replayConsistent==false. The replay did not see the same\n"
           "        conditions as the reconstruction for those, exclude them before drawing conclusions.\n",
           nInconsistent, pct(nInconsistent, nDec));
  }

  printf("\nProngs with a reconstructed track: %ld\n", nProngsRec);
  printf("  became SVertexer seeds:          %8ld  (%5.2f%%)\n", nProngsSeeded, pct(nProngsSeeded, nProngsRec));
  printf("  missing from the seeds pool:     %8ld  (%5.2f%%)\n", nProngsMissing, pct(nProngsMissing, nProngsRec));
  if (nProngsMissing) {
    printf("\nWhy the missing prongs never became seeds:\n");
    for (int i = 0; i < NSeedRej; i++) {
      const auto n = (long)hSeedRej->GetBinContent(i + 1);
      if (n) {
        printf("  %-30s %8ld  (%5.2f%% of missing)\n", seedRejNames[i], n, pct(n, nProngsMissing));
      }
    }
  }

  // ---- findable gammas: build the (pT, R) maps ----------------------------------------------
  // total loss fraction per cell, and which cut is responsible for most of it
  auto* hFndLossPtR = (TH2F*)hFndAllPtR->Clone("hFndLossPtR");
  hFndLossPtR->SetTitle("Findable-gamma prongs lost before the seeds pool;p_{T} (GeV/c);R_{conv} (cm)");
  hFndLossPtR->Reset();
  auto* hFndDomPtR = (TH2F*)hFndAllPtR->Clone("hFndDomPtR");
  hFndDomPtR->SetTitle("Dominant cut killing the prong;p_{T} (GeV/c);R_{conv} (cm)");
  hFndDomPtR->Reset();

  for (int ix = 1; ix <= nbPt; ix++) {
    for (int iy = 1; iy <= nbR; iy++) {
      const double den = hFndAllPtR->GetBinContent(ix, iy);
      double tot = 0., best = 0.;
      int dom = -1;
      for (int i = 0; i < NSeedRej; i++) {
        const double n = hFndRejPtR[i]->GetBinContent(ix, iy);
        tot += n;
        if (n > best) {
          best = n;
          dom = i;
        }
      }
      if (den > 0.) {
        hFndLossPtR->SetBinContent(ix, iy, tot / den);
      }
      // +1 so that 0 reads as "no loss / no stats" in the palette
      hFndDomPtR->SetBinContent(ix, iy, (tot >= MinLossForDominant && dom >= 0) ? dom + 1 : 0.);
    }
  }
  hFndLossPtR->SetMinimum(0.);
  hFndLossPtR->SetMaximum(1.);
  hFndDomPtR->SetMinimum(-0.5);
  hFndDomPtR->SetMaximum(NSeedRej + 0.5);
  hFndDomPtR->SetContour(NSeedRej + 1);

  // 1D breakdowns: fraction of findable-gamma prongs lost to each cut, stacked so the total height
  // is the overall seeding loss
  auto* stPt = new THStack("stFndLossPt", "Findable-gamma prongs lost before the seeds pool;p_{T} (GeV/c);lost fraction");
  auto* stR = new THStack("stFndLossR", "Findable-gamma prongs lost before the seeds pool;R_{conv} (cm);lost fraction");
  auto* leg = new TLegend(0.55, 0.55, 0.98, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.032);
  for (int i = 0; i < NSeedRej; i++) {
    hFndRejPt[i]->Divide(hFndRejPt[i], hFndAllPt, 1., 1., "B");
    hFndRejR[i]->Divide(hFndRejR[i], hFndAllR, 1., 1., "B");
    for (auto* h : {(TH1*)hFndRejPt[i], (TH1*)hFndRejR[i]}) {
      h->SetFillColor(seedRejColors[i]);
      h->SetLineColor(seedRejColors[i]);
      h->SetMarkerColor(seedRejColors[i]);
    }
    if (hFndRejPt[i]->Integral() > 0. || hFndRejR[i]->Integral() > 0.) {
      stPt->Add(hFndRejPt[i], "hist");
      stR->Add(hFndRejR[i], "hist");
      leg->AddEntry(hFndRejPt[i], seedRejNames[i], "f");
    }
  }

  printf("\nFindable gammas (both prongs reconstructed): %ld of %ld decays (%5.2f%%)\n",
         nFindable, nDec, pct(nFindable, nDec));
  if (nFndProngs) {
    printf("  their prongs:                    %8ld\n", nFndProngs);
    printf("  lost before the seeds pool:      %8ld  (%5.2f%%)\n", nFndMissing, pct(nFndMissing, nFndProngs));
  }

  printf("\nJoint vs single-prong seeding (the prongs are correlated, do not square eps1):\n");
  printf("  %-22s %8s %8s %8s %8s %8s\n", "composition", "gammas", "eps1", "eps2", "eps1^2", "eps2/eps1^2");
  long totFnd = 0, totProng = 0, totBoth = 0;
  for (int ic = 0; ic < NComp; ic++) {
    if (!nCompFnd[ic]) {
      continue;
    }
    totFnd += nCompFnd[ic];
    totProng += nCompProngSeeded[ic];
    totBoth += nCompBothSeeded[ic];
    const double e1 = double(nCompProngSeeded[ic]) / (2. * nCompFnd[ic]);
    const double e2 = double(nCompBothSeeded[ic]) / nCompFnd[ic];
    printf("  %-22s %8ld %7.2f%% %7.2f%% %7.2f%% %8.2f\n", compNames[ic], nCompFnd[ic],
           100. * e1, 100. * e2, 100. * e1 * e1, e1 > 0. ? e2 / (e1 * e1) : 0.);
  }
  if (totFnd) {
    const double e1 = double(totProng) / (2. * totFnd);
    const double e2 = double(totBoth) / totFnd;
    printf("  %-22s %8ld %7.2f%% %7.2f%% %7.2f%% %8.2f\n", "all", totFnd,
           100. * e1, 100. * e2, 100. * e1 * e1, e1 > 0. ? e2 / (e1 * e1) : 0.);
  }
  printf("  eps1 = per-prong seeding, eps2 = both prongs seeded. eps2/eps1^2 > 1 means the prongs\n"
         "  fail together, so pair-level estimates built on eps1^2 are too pessimistic by that factor.\n");

  const double nPaired = hRejV0->Integral();
  if (nPaired > 0.) {
    printf("\nOf the %.0f gammas whose prongs were both seeded and tried as a pair, checkV0 says:\n", nPaired);
    for (int i = 0; i < NV0Rej; i++) {
      const auto n = (long)hRejV0->GetBinContent(i + 1);
      if (n) {
        printf("  %-22s %8ld  (%5.2f%%)\n", v0RejNames[i], n, pct(n, (long)nPaired));
      }
    }
  }

  gStyle->SetOptStat(0);
  auto* cnv = new TCanvas("cSVGamma", "SVertexer seeding of conversion prongs", 1400, 900);
  cnv->Divide(3, 2);
  cnv->cd(1);
  hStage->Draw("hist");
  cnv->cd(2);
  gPad->SetBottomMargin(0.28);
  hSeedRej->LabelsOption("v");
  hSeedRej->Draw("hist");
  cnv->cd(3);
  gPad->SetLeftMargin(0.28);
  hSeedRejVsSrc->Draw("colz");
  cnv->cd(4);
  gPad->SetLeftMargin(0.28);
  hSeedRejVsR->Draw("colz");
  cnv->cd(5);
  gPad->SetLeftMargin(0.28);
  hSeedRejVsPt->Draw("colz");
  cnv->cd(6);
  // seeding efficiency vs radius, the one plot that says where in the detector prongs are lost
  auto* hEffR = (TH1F*)hSeededR->Clone("hEffR");
  hEffR->SetTitle("Seeding efficiency of reconstructed prongs;R_{conv} (cm);seeded / reconstructed");
  hEffR->Divide(hSeededR, hAllR, 1., 1., "B");
  hEffR->SetMinimum(0.);
  hEffR->SetMaximum(1.05);
  hEffR->Draw("e");

  // second canvas: the findable-gamma seeding map, which cut kills the prong where
  auto* cnv2 = new TCanvas("cSVGammaFindable", "Seeding losses of findable gammas", 1400, 900);
  cnv2->Divide(2, 2);
  cnv2->cd(1);
  gPad->SetRightMargin(0.15);
  hFndLossPtR->Draw("colz");
  cnv2->cd(2);
  gPad->SetRightMargin(0.15);
  {
    // discrete palette so each cell reads as a cut rather than as a number
    std::vector<int> pal{kWhite};
    for (int i = 0; i < NSeedRej; i++) {
      pal.push_back(seedRejColors[i]);
    }
    gStyle->SetPalette((int)pal.size(), pal.data());
    hFndDomPtR->Draw("col");
    auto* legDom = new TLegend(0.62, 0.55, 0.98, 0.92);
    legDom->SetBorderSize(0);
    legDom->SetFillStyle(0);
    legDom->SetTextSize(0.030);
    for (int i = 0; i < NSeedRej; i++) {
      if (hFndRejPtR[i]->GetEntries() > 0) {
        auto* mk = new TH1F(Form("hLegDom_%d", i), "", 1, 0., 1.);
        mk->SetFillColor(seedRejColors[i]);
        legDom->AddEntry(mk, seedRejNames[i], "f");
      }
    }
    legDom->Draw();
  }
  cnv2->cd(3);
  gPad->SetBottomMargin(0.30);
  hRejV0->LabelsOption("v");
  hRejV0->Draw("hist");
  cnv2->cd(4);
  // joint seeding efficiency against the independent expectation, as a function of radius
  auto* hEps2R = (TH1F*)hBothSeededR->Clone("hEps2R");
  hEps2R->SetTitle("Seeding of findable gammas;R_{conv} (cm);efficiency");
  hEps2R->Divide(hBothSeededR, hFndDecR, 1., 1., "B");
  auto* hEps1R = (TH1F*)hProngSeededR->Clone("hEps1R");
  hEps1R->Divide(hProngSeededR, hFndDecR, 0.5, 1.); // 0.5 because each gamma has two prongs
  auto* hEps1SqR = (TH1F*)hEps1R->Clone("hEps1SqR");
  hEps1SqR->Multiply(hEps1R);
  hEps2R->SetLineColor(kRed + 1);
  hEps2R->SetLineWidth(2);
  hEps1SqR->SetLineColor(kAzure + 2);
  hEps1SqR->SetLineWidth(2);
  hEps1SqR->SetLineStyle(2);
  hEps2R->SetMinimum(0.);
  hEps2R->SetMaximum(1.05);
  hEps2R->Draw("e");
  hEps1SqR->Draw("hist same");
  auto* legJ = new TLegend(0.45, 0.75, 0.98, 0.92);
  legJ->SetBorderSize(0);
  legJ->SetFillStyle(0);
  legJ->AddEntry(hEps2R, "both prongs seeded (measured)", "l");
  legJ->AddEntry(hEps1SqR, "#varepsilon_{1}^{2} (independence assumed)", "l");
  legJ->Draw();

  std::unique_ptr<TFile> flOut{TFile::Open(outFile.c_str(), "recreate")};
  hFndAllPtR->Write();
  hFndLossPtR->Write();
  hFndDomPtR->Write();
  stPt->Write();
  stR->Write();
  for (int i = 0; i < NSeedRej; i++) {
    hFndRejPtR[i]->Write();
  }
  hRejV0->Write();
  hRejV0VsR->Write();
  stPt->Write();
  hFndDecR->Write();
  hBothSeededR->Write();
  hProngSeededR->Write();
  hEps2R->Write();
  hEps1R->Write();
  hEps1SqR->Write();
  cnv2->Write();
  cnv2->SaveAs((outFile.substr(0, outFile.find_last_of('.')) + "_findable.png").c_str());
  hStage->Write();
  hSeedRej->Write();
  hSeedRejVsR->Write();
  hSeedRejVsPt->Write();
  hSeedRejVsSrc->Write();
  hSeededR->Write();
  hSeededPt->Write();
  hAllR->Write();
  hAllPt->Write();
  hEffR->Write();
  cnv->Write();
  cnv->SaveAs((outFile.substr(0, outFile.find_last_of('.')) + ".png").c_str());
  printf("\nWritten %s\n", outFile.c_str());
}
