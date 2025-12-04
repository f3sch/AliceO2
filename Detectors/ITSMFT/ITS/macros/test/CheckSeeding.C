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

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <vector>
#include <unistd.h>

#include <TF1.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>

#include "Framework/Logger.h"
#include "ITStracking/Definitions.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "Steer/MCKinematicsReader.h"
#endif

namespace fs = std::filesystem;
void fitGaus(TH1*);

void CheckSeeding(const std::string& trkFileName = "o2trac_its.root", const std::string& colFileName = "collisioncontext.root")
{
  fs::path base = fs::current_path();
  std::vector<fs::path> dirs;
  for (auto& entry : fs::recursive_directory_iterator(base)) {
    if (!entry.is_regular_file()) {
      continue;
    }
    if (entry.path().filename() == trkFileName) {
      dirs.push_back(entry.path().parent_path());
    }
  }
  std::sort(dirs.begin(), dirs.end());
  dirs.erase(std::unique(dirs.begin(), dirs.end()), dirs.end());

  LOG(info) << "Found " << dirs.size() << " directories containing " << trkFileName << "\n";
  if (dirs.empty()) {
    return;
  }

  const int nBins{100};
  const float yRange{0.5}, zRange{18};
  auto hVtxX = new TH1F("hVtxX", ";x (cm)", nBins, -yRange, yRange);
  auto hVtxY = new TH1F("hVtxY", ";y (cm)", nBins, -yRange, yRange);
  auto hVtxZ = new TH1F("hVtxZ", ";z (cm)", nBins, -zRange, zRange);
  const float resRange{0.1};
  auto hVtxResX = new TH1F("hVtxResX", ";x (cm)", nBins, -resRange, resRange);
  auto hVtxResY = new TH1F("hVtxResY", ";y (cm)", nBins, -resRange, resRange);
  auto hVtxResZ = new TH1F("hVtxResZ", ";z (cm)", nBins, -resRange, resRange);
  const float pullRange{5};
  auto hVtxPullX = new TH1F("hVtxPullX", "pull x", nBins, -pullRange, pullRange);
  auto hVtxPullY = new TH1F("hVtxPullY", "pull y", nBins, -pullRange, pullRange);
  auto hVtxPullZ = new TH1F("hVtxPullZ", "pull z", nBins, -pullRange, pullRange);

  auto hVtxEffNumZ = new TH1F("hVtxEffNumZ", ";z (cm);efficiency", nBins, -zRange, zRange);
  auto hVtxEffDenZ = new TH1F("hVtxEffDenZ", ";z (cm)", nBins, -zRange, zRange);
  auto hVtxPurity = new TH1F("hVtxPurity", ";purity", nBins, 0, 1);
  auto hVtxPurityVsNContrib = new TH2F("hVtxPurityVsNContrib", ";purity;contributors", nBins, 0, 1, 150, 0, 150);
  auto hVtxNContribVsNPrim = new TH2F("hVtxNContribVsNPrim", ";contributors;nprim", 150, 0, 150, 150, 0, 150);

  std::vector<o2::its::Vertex>* vertices{nullptr};
  std::vector<o2::MCCompLabel>* verticesLbl{nullptr};
  std::vector<float>* verticesPurity{nullptr};
  o2::steer::MCKinematicsReader* mcReader{nullptr};

  for (const auto& dirPath : dirs) {
    std::string dirStr = dirPath.string();
    LOG(info) << "Processing directory: " << dirStr << "\n";

    // switch working directory to the directory containing the files
    if (chdir(dirStr.c_str()) != 0) {
      perror("chdir failed");
      LOG(info) << "Skipping directory " << dirStr << " due to chdir failure\n";
      continue;
    }

    // open the track ROOT file from this directory
    fs::path trkPath = dirPath / trkFileName;
    if (!fs::exists(trkPath)) {
      LOG(info) << "Missing track file " << trkPath.string() << " (skipping)\n";
      continue;
    }

    TFile* f = TFile::Open(trkPath.string().c_str(), "READ");
    if (!f || f->IsZombie()) {
      LOG(info) << "Failed to open " << trkPath.string() << " (skipping)\n";
      if (f) {
        f->Close();
        delete f;
      }
      continue;
    }

    TTree* trkTree = dynamic_cast<TTree*>(f->Get("o2sim"));
    if (!trkTree) {
      LOG(info) << "Tree 'o2sim' not found in " << trkPath.string() << " (skipping)\n";
      f->Close();
      delete f;
      continue;
    }

    vertices = nullptr;
    verticesLbl = nullptr;
    verticesPurity = nullptr;
    trkTree->SetBranchAddress("ITSVertices", &vertices);
    trkTree->SetBranchAddress("ITSVertexMCTruth", &verticesLbl);
    trkTree->SetBranchAddress("ITSVertexMCPurity", &verticesPurity);
    fs::path colPath = dirPath / colFileName;
    if (!fs::exists(colPath)) {
      LOG(info) << "Collision file not found: " << colPath.string() << " (mcReader will be nullptr)\n";
      mcReader = nullptr;
    } else {
      delete mcReader;
      mcReader = new o2::steer::MCKinematicsReader(colPath.string().c_str());
    }

    Long64_t nEntries = trkTree->GetEntries();
    LOG(info) << "Entries in file: " << nEntries << "\n";

    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
      trkTree->GetEntry(iEntry);

      for (size_t iVtx{0}; iVtx < vertices->size(); ++iVtx) {
        const auto& vtx = vertices->at(iVtx);
        const auto& vtxLbl = verticesLbl->at(iVtx);
        const auto& purity = verticesPurity->at(iVtx);
        LOGP(info, "{}: {}", iVtx, vtx.asString());
        if (vtxLbl.isValid() && vtxLbl.isCorrect() && mcReader) {
          const auto& head = mcReader->getMCEventHeader(vtxLbl.getSourceID(), vtxLbl.getEventID());
          LOGP(info, "\t-{}", vtxLbl.asString());
          LOGP(info, "\t-Purity:{}", purity);
          LOGP(info, "\t-MC: x={} y={} z={} prim={}", head.GetX(), head.GetY(), head.GetZ(), head.GetNPrim());

          hVtxX->Fill(vtx.getX());
          hVtxY->Fill(vtx.getY());
          hVtxZ->Fill(vtx.getZ());

          hVtxResX->Fill(vtx.getX() - head.GetX());
          hVtxResY->Fill(vtx.getY() - head.GetY());
          hVtxResZ->Fill(vtx.getZ() - head.GetZ());

          hVtxPullX->Fill((vtx.getX() - head.GetX()) / vtx.getSigmaX());
          hVtxPullY->Fill((vtx.getY() - head.GetY()) / vtx.getSigmaY());
          hVtxPullZ->Fill((vtx.getZ() - head.GetZ()) / vtx.getSigmaZ());

          hVtxEffNumZ->Fill(head.GetZ());

          hVtxNContribVsNPrim->Fill(vtx.getNContributors(), head.GetNPrim());

        } else {
          LOGP(info, "\t-FAKE");
        }

        hVtxPurity->Fill(purity);
        hVtxPurityVsNContrib->Fill(purity, vtx.getNContributors());
      }
    }

    // fill den
    for (int iEve{0}; iEve < (int)mcReader->getNEvents(0); ++iEve) {
      const auto& head = mcReader->getMCEventHeader(0, iEve);
      hVtxEffDenZ->Fill(head.GetZ());
    }

    if (mcReader) {
      delete mcReader;
      mcReader = nullptr;
    }
    trkTree = nullptr;
    f->Close();
    delete f;

    if (chdir(base.string().c_str()) != 0) {
      perror("chdir back failed");
    }
  }

  auto c = new TCanvas();
  c->Divide(3, 3);
  c->cd(1);
  hVtxX->Draw();
  fitGaus(hVtxX);
  gPad->SetLogy();
  c->cd(2);
  hVtxY->Draw();
  fitGaus(hVtxY);
  gPad->SetLogy();
  c->cd(3);
  hVtxZ->Draw();
  fitGaus(hVtxZ);
  gPad->SetLogy();
  c->cd(4);
  hVtxResX->Draw();
  c->cd(5);
  hVtxResY->Draw();
  c->cd(6);
  hVtxResZ->Draw();
  c->cd(7);
  hVtxPullX->Draw();
  fitGaus(hVtxPullX);
  gPad->SetLogy();
  c->cd(8);
  hVtxPullY->Draw();
  fitGaus(hVtxPullY);
  gPad->SetLogy();
  c->cd(9);
  hVtxPullZ->Draw();
  fitGaus(hVtxPullZ);
  gPad->SetLogy();
  c->Draw();

  hVtxEffNumZ->Divide(hVtxEffNumZ, hVtxEffDenZ, 1., 1., "B");
  c = new TCanvas();
  c->Divide(2, 2);
  c->cd(1);
  hVtxEffNumZ->Draw();
  c->cd(2);
  hVtxPurity->Draw();
  c->cd(3);
  hVtxPurityVsNContrib->Draw("colz");
  c->cd(4);
  hVtxNContribVsNPrim->Draw("colz");
  c->Draw();
}

void fitGaus(TH1* h)
{
  if (!h) {
    return;
  }

  // fit
  TF1* f = new TF1("fG", "gaus");
  h->Fit(f, "QMS"); // quiet fit

  // draw histogram
  h->Draw();

  // get parameters
  double mean = f->GetParameter(1);
  double sigma = f->GetParameter(2);
  double chi2 = f->GetChisquare();
  double ndf = f->GetNDF();

  // add text box
  auto t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.04);
  t->DrawLatex(0.15, 0.85, Form("mean = %.4f", mean));
  t->DrawLatex(0.15, 0.80, Form("sigma = %.4f", sigma));
  t->DrawLatex(0.15, 0.75, Form("#chi^{2}/ndf = %.2f / %.0f", chi2, ndf));
}
