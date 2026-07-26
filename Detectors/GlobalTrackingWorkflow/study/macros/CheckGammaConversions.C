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

/// \file CheckGammaConversions.C
/// \brief Photon-conversion efficiencies from the dedicated TrackMCStudy gammaConv tree
///
/// Usage:
///   root -l -b -q 'CheckGammaConversions.C+("trackMCStudy.root")'

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>

#include <array>
#include <memory>
#include <string>

#include "GlobalTrackingStudy/TrackMCStudyTypes.h"
#endif

namespace
{
constexpr int NRBins = 90;
constexpr float RMin = 0.f;
constexpr float RMax = 90.f;
constexpr std::array<float, 8> LayerRadii{2.33959f, 3.14076f, 3.91924f, 19.6213f,
                                          24.5597f, 34.388f, 39.3329f, 83.5f};
constexpr std::array<float, 7> PtEdges{0.01f, 0.1f, 0.2f, 0.5f, 1.f, 2.f, 9999999.f};

const char* outcomeNames[] = {
  "V0 stored", "both eligible legs, no V0", "neither leg eligible",
  "no eligible positron", "no eligible electron"};
constexpr int NOutcomes = sizeof(outcomeNames) / sizeof(char*);

int getPtSlice(float pt)
{
  for (int i = 0; i + 1 < (int)PtEdges.size(); i++) {
    if (pt >= PtEdges[i] &&
        (pt < PtEdges[i + 1] || (i + 2 == (int)PtEdges.size() && pt <= PtEdges[i + 1]))) {
      return i;
    }
  }
  return -1;
}

void setLabels(TAxis* axis, const char* const* labels, int n)
{
  for (int i = 0; i < n; i++) {
    axis->SetBinLabel(i + 1, labels[i]);
  }
}

void addLayerLines()
{
  gPad->Update();
  const double ymin = gPad->GetUymin();
  const double ymax = gPad->GetUymax();
  for (size_t i = 0; i < LayerRadii.size(); i++) {
    auto* line = new TLine(LayerRadii[i], ymin, LayerRadii[i], ymax);
    line->SetLineColor(i + 1 == LayerRadii.size() ? kBlue + 2 : kGray + 1);
    line->SetLineStyle(i + 1 == LayerRadii.size() ? 2 : 3);
    line->Draw();
  }
}
} // namespace

void CheckGammaConversions(const std::string& inpFile = "trackMCStudy.root",
                           const std::string& outFile = "CheckGammaConversions.root")
{
  using namespace o2::trackstudy;
  constexpr int NPtSlices = PtEdges.size() - 1;
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(0);

  std::unique_ptr<TFile> input{TFile::Open(inpFile.c_str())};
  if (!input || input->IsZombie()) {
    printf("ERROR: cannot open %s\n", inpFile.c_str());
    return;
  }
  auto* tree = static_cast<TTree*>(input->Get("gammaConv"));
  if (!tree || !tree->GetBranch("conv")) {
    printf("ERROR: no gammaConv/conv output in %s\n", inpFile.c_str());
    return;
  }

  std::array<TH1F*, NPtSlices> hTruth{}, hBothAnywhere{}, hBothReference{}, hRawAnywhere{};
  std::array<TH1F*, NPtSlices> hBothAnywhereEff{}, hBothReferenceEff{}, hHeadline{};
  std::array<TH2F*, NPtSlices> hOutcome{};
  for (int ipt = 0; ipt < NPtSlices; ipt++) {
    const std::string suffix = "Pt" + std::to_string(ipt);
    const std::string ptTitle = std::to_string(PtEdges[ipt]) + " <= pT(gamma) < " +
                                std::to_string(PtEdges[ipt + 1]);
    hTruth[ipt] = new TH1F(("hTruthPhotons" + suffix).c_str(),
                           (ptTitle + ";true R_{conv} (cm);selected photons with reference collision").c_str(),
                           NRBins, RMin, RMax);
    hBothAnywhere[ipt] = new TH1F(("hBothLegsFoundAnywhere" + suffix).c_str(), "",
                                  NRBins, RMin, RMax);
    hBothReference[ipt] = new TH1F(("hBothLegsFound" + suffix).c_str(), "",
                                   NRBins, RMin, RMax);
    hRawAnywhere[ipt] = new TH1F(("hRawV0ExactAnywhere" + suffix).c_str(), "",
                                 NRBins, RMin, RMax);
    hBothAnywhereEff[ipt] = new TH1F(("hBothLegsFoundAnywhereEfficiency" + suffix).c_str(),
                                     (ptTitle + ";true R_{conv} (cm);eligible pair anywhere / truth photon").c_str(),
                                     NRBins, RMin, RMax);
    hBothReferenceEff[ipt] = new TH1F(("hBothLegsFoundEfficiency" + suffix).c_str(),
                                      (ptTitle + ";true R_{conv} (cm);eligible pair in reference collision / truth photon").c_str(),
                                      NRBins, RMin, RMax);
    hHeadline[ipt] = new TH1F(("hEfficiencyV0Anywhere" + suffix).c_str(),
                              (ptTitle + ";true R_{conv} (cm);exact raw V0 / eligible pair anywhere").c_str(),
                              NRBins, RMin, RMax);
    hOutcome[ipt] = new TH2F(("hTruthPairOutcome" + suffix).c_str(),
                             (ptTitle + ";true R_{conv} (cm);exclusive outcome").c_str(),
                             NRBins, RMin, RMax, NOutcomes, -0.5, NOutcomes - 0.5);
    setLabels(hOutcome[ipt]->GetYaxis(), outcomeNames, NOutcomes);
  }

  GammaConvInfo* conversion = nullptr;
  tree->SetBranchAddress("conv", &conversion);
  const auto nEntries = tree->GetEntries();
  for (Long64_t entry = 0; entry < nEntries; entry++) {
    tree->GetEntry(entry);
    if (!conversion) {
      continue;
    }
    const int ipt = getPtSlice(conversion->photonPt);
    if (ipt < 0) {
      continue;
    }
    const float radius = conversion->getTrueConversionRadius();
    hTruth[ipt]->Fill(radius);
    if (conversion->bothLegsFoundAnywhere) {
      hBothAnywhere[ipt]->Fill(radius);
    }
    if (conversion->bothLegsFoundInReferencePV) {
      hBothReference[ipt]->Fill(radius);
    }
    if (conversion->rawV0FoundUsingAnywhereEligiblePair) {
      hRawAnywhere[ipt]->Fill(radius);
    }
    hOutcome[ipt]->Fill(radius, conversion->terminalReason);
  }

  std::unique_ptr<TFile> output{TFile::Open(outFile.c_str(), "recreate")};
  for (int ipt = 0; ipt < NPtSlices; ipt++) {
    hBothAnywhereEff[ipt]->Divide(hBothAnywhere[ipt], hTruth[ipt], 1., 1., "B");
    hBothReferenceEff[ipt]->Divide(hBothReference[ipt], hTruth[ipt], 1., 1., "B");
    hHeadline[ipt]->Divide(hRawAnywhere[ipt], hBothAnywhere[ipt], 1., 1., "B");

    for (auto* hist : {hTruth[ipt], hBothAnywhere[ipt], hBothReference[ipt], hRawAnywhere[ipt],
                       hBothAnywhereEff[ipt], hBothReferenceEff[ipt], hHeadline[ipt]}) {
      hist->Write();
    }
    hOutcome[ipt]->Write();

    auto* canvas = new TCanvas(("cGammaConversionPt" + std::to_string(ipt)).c_str(),
                               "Gamma conversion efficiencies", 1600, 900);
    canvas->Divide(2, 2);
    canvas->cd(1);
    hBothAnywhereEff[ipt]->SetMinimum(0.);
    hBothAnywhereEff[ipt]->SetMaximum(1.05);
    hBothAnywhereEff[ipt]->Draw("hist");
    addLayerLines();
    canvas->cd(2);
    hBothReferenceEff[ipt]->SetMinimum(0.);
    hBothReferenceEff[ipt]->SetMaximum(1.05);
    hBothReferenceEff[ipt]->Draw("hist");
    addLayerLines();
    canvas->cd(3);
    hHeadline[ipt]->SetMinimum(0.);
    hHeadline[ipt]->SetMaximum(1.05);
    hHeadline[ipt]->Draw("hist");
    addLayerLines();
    canvas->cd(4);
    hOutcome[ipt]->Draw("colz");
    addLayerLines();
    canvas->Write();
    const auto pngName = outFile.substr(0, outFile.find_last_of('.')) + "_pt" +
                         std::to_string(ipt) + ".png";
    canvas->SaveAs(pngName.c_str());
  }
  printf("Written %s\n", outFile.c_str());
}
