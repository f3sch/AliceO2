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

// AOD benchmark tool

#include "boost/program_options.hpp"
#include "Framework/Logger.h"
#include "fmt/format.h"
#include "O2Version.h"

#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TBranch.h"
#include "TTreeCache.h"
#include "TStopwatch.h"
#include "TFileCacheRead.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "ROOT/InternalTreeUtils.hxx"

#include "aodReadSpeed.h"

#include <vector>
#include <string>
#include <filesystem>
#include <memory>
#include <algorithm>
#include <functional>
#include <regex>

namespace bpo = boost::program_options;
namespace fs = std::filesystem;
using CRef = std::reference_wrapper<const Result>;

namespace
{
Results gResultMap;
bpo::variables_map gVM;
std::unordered_map<std::string, std::vector<CRef>> gRefMap;
std::unordered_map<std::string, CRef> gRefSummaryMap;
TreeMap gTreeMap;
RegexHandler gTreeRegex;
RegexHandler gBranchRegex;
ProgramOptions gOpts;
} // namespace

bool initOptionsAndParse(bpo::options_description& options, int argc, char* argv[])
{
  bpo::options_description general_options("General options");
  general_options.add_options()(
    "input,i", bpo::value<std::string>(&gOpts.input)->default_value(""), "Input files (separated by semicolon file0;file1)")(
    "input-regex", bpo::value<std::string>(&gOpts.inputReg)->default_value(""), "Input files (separated by semicolon file0;file1)")(
    "output,o", bpo::value<std::string>(&gOpts.output)->default_value("aod-readspeed.root"), "Output file name")(
    "output-prefix", bpo::value<std::string>(&gOpts.outputSuffix)->default_value(""), "Output prefix file name")(
    "repeat,r", bpo::value<unsigned int>(&gOpts.repeat)->default_value(10), "Number of repeatitions")(
    "batch,p", bpo::bool_switch(&gOpts.batch)->default_value(false), "Enable batch processing")(
    "verbose,v", bpo::bool_switch(&gOpts.verbose)->default_value(false), "Enable verbose statisics printing")(
    "debug,d", bpo::bool_switch(&gOpts.debug)->default_value(false), "Enable debug printing")(
    "help,h", "Produce help message.");
  options.add(general_options);

  bpo::options_description result_options("Result options");
  result_options.add_options()(
    "ratio-default", bpo::value<std::string>(&gOpts.ratioDefault)->default_value(""), "Produce ratio plot with this key as baseline")(
    "disable-pdf", bpo::bool_switch(&gOpts.disablePDF)->default_value(false), "Disable pdf production")(
    "disable-root", bpo::bool_switch(&gOpts.disableRoot)->default_value(false), "Disable root production")(
    "disable-graphs", bpo::bool_switch(&gOpts.disableGraphs)->default_value(false), "Disable graphs")(
    "disable-system-info", bpo::bool_switch(&gOpts.disableSystemInfo)->default_value(false), "Print system info; only in verbose output")(
    "disable-summary", bpo::bool_switch(&gOpts.disableSummary)->default_value(false), "Disable summary")(
    "disable-summary-plot", bpo::bool_switch(&gOpts.disableSummaryPlot)->default_value(false), "Disable summary plot")(
    "disable-aod-summary", bpo::bool_switch(&gOpts.disableAODSummary)->default_value(false), "Disable overall AO2D summary")(
    "disable-result-supression", bpo::bool_switch(&gOpts.disableResultSupression)->default_value(false), "Disable inconclusive result supression")(
    "disable-empty-results-truncation", bpo::bool_switch(&gOpts.disableEmptyResultTruncation)->default_value(false), "Disable empty results truncation");
  options.add(result_options);

  bpo::options_description file_options("File options");
  file_options.add_options()(
    "disable-file-prefetch", bpo::bool_switch(&gOpts.disableFilePrefetch)->default_value(false), "Disable file prefetching")(
    "file-readahead-size", bpo::value<Int_t>(&gOpts.fileReadAheadSize)->default_value(-1), "Set File readahead size");
  options.add(file_options);

  bpo::options_description branch_options("Branch options");
  branch_options.add_options()(
    "disable-all-branches", bpo::bool_switch()->default_value(false), "report on all branches")(
    "branch-regex-inclusive", bpo::value<std::string>(&gOpts.branchIncReg)->default_value(""), "Analyse only trees matching this regex (C++ regex)")(
    "branch-regex-exclusive", bpo::value<std::string>(&gOpts.branchExReg)->default_value(""), "Analyse only trees not matching this regex (C++ regex)");
  options.add(branch_options);

  bpo::options_description tree_options("Tree options");
  tree_options.add_options()(
    "tree-regex-inclusive", bpo::value<std::string>(&gOpts.treeIncReg)->default_value(""), "Analyse only trees matching this regex (C++ regex)")(
    "tree-regex-exclusive", bpo::value<std::string>(&gOpts.treeExReg)->default_value(""), "Analyse only trees not matching this regex (C++ regex)")(
    "tree-cache-size", bpo::value<Int_t>(&gOpts.treeCacheSize)->default_value(-1), "tree cache size");
  options.add(tree_options);

  // Try parsing
  try {
    bpo::store(parse_command_line(argc, argv, options), gVM);

    // help
    if (gVM.count("help") != 0u) {
      LOG(info) << options;
      LOGP(info, "GitVersion: {}   Revision: {}", o2::fullVersion(), o2::gitRevision());
      LOG(info) << o2::getBuildInfo();
      return false;
    }

    bpo::notify(gVM);
  } catch (const bpo::error& e) {
    LOG(error) << e.what();
    LOG(error) << "Error parsing command line arguments; Available options:";
    LOG(error) << options;
    return false;
  }

  // Root specific stuff
  gErrorIgnoreLevel = kError;
  if (gOpts.verbose) {
    gErrorIgnoreLevel = kInfo;
  }
  if (gOpts.debug) {
    gErrorIgnoreLevel = kPrint;
  }

  return true;
}

void readTree(TFile* file, TTree* tree, TStopwatch& sw, Result& res)
{

  // Pick active branches
  tree->SetBranchStatus("*", false);
  std::vector<TBranch*> branches;
  for (const auto& branchName : ROOT::Internal::TreeUtils::GetTopLevelBranchNames(*tree)) {
    auto b = tree->GetBranch(branchName.c_str());
    if (gBranchRegex.matchesAny(branchName)) {
      b->SetStatus(true);
      branches.push_back(b);
    }
  }

  if (!gOpts.disableSystemInfo) {
    gSystem->GetProcInfo(&res.fBefore);
  }

  ULong64_t bytesRead = 0;
  const ULong64_t fileStartBytes = file->GetBytesRead();

  sw.Start(true);

  for (Long64_t iEntry{0}; tree->LoadTree(iEntry) >= 0; ++iEntry) {
    for (auto b : branches) {
      bytesRead += b->GetEntry(iEntry);
    }
  }

  sw.Stop();

  if (!gOpts.disableSystemInfo) {
    gSystem->GetProcInfo(&res.fAfter);
  }

  const ULong64_t fileBytesRead = file->GetBytesRead() - fileStartBytes;

  res.fUncompressedBytesReadTot = bytesRead;
  res.fCompressedBytesReadTot = fileBytesRead;
  res.fRealTime = sw.RealTime();
  res.fCpuTime = sw.CpuTime();
}

void evalThrougputST(TFile* file, TTree* tree)
{
  LOG_IF(info, gOpts.debug) << "EvalThrougputST for " << tree->GetName() << " in " << file->GetName();

  const int comp = file->GetCompressionSettings();
  const std::string name = std::to_string(comp) + "." + tree->GetName();

  TStopwatch sw;

  Result res;
  res.fAlgorithm = file->GetCompressionAlgorithm(),
  res.fLevel = file->GetCompressionLevel(),
  res.fName = name,
  res.fCounter = 1,
  res.fRepeat = gOpts.repeat,

  // read actual tree
    readTree(file, tree, sw, res);

  if (gResultMap.find(name) != gResultMap.end()) {
    gResultMap[name] += res;
  } else {
    gResultMap[name] = res;
  }
}

void produceResults()
{
  LOGP(info, "Starting Benchmark Loop");
  TStopwatch sw;
  for (const auto& [treeName, mapping] : gTreeMap) {
    sw.Start(true);
    for (auto i{gOpts.repeat}; i--;) {
      for (const auto& [fileName, dfs] : mapping.fFileMap) {
        LOG_IF(info, gOpts.debug) << "  `-> " << fileName;
        for (const auto& df : dfs) {
          LOG_IF(info, gOpts.debug) << "     `-> " << df;
          // opening the file again for every data frame avoids caching and resets counters
          std::unique_ptr<TFile> file(TFile::Open(fileName.c_str(), "READ"));
          if (!file || file->IsZombie()) {
            LOGP(error, "File '{}' cannot be opened or is zombie!", fileName);
            continue;
          }

          if (auto fc = file->GetCacheRead(); gOpts.disableFilePrefetch && fc != nullptr) {
            fc->SetEnablePrefetching(false);
          }

          if (gOpts.fileReadAheadSize >= 0) {
            file->SetReadaheadSize(gOpts.fileReadAheadSize);
          }

          auto treePath = df + "/";
          treePath += treeName;
          std::unique_ptr<TTree> tree(file->Get<TTree>(treePath.c_str()));
          if (tree == nullptr) {
            LOGP(error, "Failed retrieving {} from {} in {}", treeName, df, fileName);
            continue;
          }
          tree->SetCacheSize(gOpts.treeCacheSize);

          evalThrougputST(file.get(), tree.get());

          // Close out all caches
          if (auto cache = file->GetCacheRead(); cache != nullptr) {
            cache->Close();
            LOG_IF(info, gOpts.debug) << "Dropped file cache for " << fileName;
          }
        }
      }
    }
    sw.Stop();
    LOGP(info, "  `-> Finished Tree {} in {:.2f} s", treeName, sw.RealTime());
  }
}

void makeTreeMap(const std::string& fileName)
{
  std::unique_ptr<TFile> file(TFile::Open(fileName.data(), "READ"));
  if (!file || file->IsZombie()) {
    LOGP(error, "File '{}' cannot be opened or is zombie!", fileName);
    return;
  }
  for (auto keyDF : *file->GetListOfKeys()) {
    const std::string df = (dynamic_cast<TKey*>(keyDF))->ReadObj()->GetName();
    if (!df.starts_with("DF_")) {
      continue;
    }
    LOG_IF(info, gOpts.debug) << "Entering in " << fileName << " -> " << df;
    auto dir = file->GetDirectory(df.c_str());
    for (auto key : *dir->GetListOfKeys()) {
      TObject* keyTree = (dynamic_cast<TKey*>(key))->ReadObj();
      if (!keyTree->InheritsFrom(TTree::Class())) {
        continue;
      }
      const std::string tree = keyTree->GetName();

      LOG_IF(info, gOpts.debug) << "Tree " << tree;
      // Check if regex match
      if (!gTreeRegex.matchesAny(tree)) {
        LOG_IF(info, gOpts.debug) << "   not accepted";
        continue;
      }
      LOG_IF(info, gOpts.debug) << "   accepted";

      // Add to Map
      if (gTreeMap.find(tree) != gTreeMap.end()) {
        gTreeMap[tree].addEntry(fileName, df);
      } else {
        gTreeMap[tree] = {.fTreeName = tree};
      }
    }
  }
  LOG_IF(info, gOpts.verbose) << "Completed read of file";
}

void reportResults()
{
  LOGP(info, "{:~^100}", " Reporting START ");
  for (const auto& [name, vres] : gRefMap) {
    LOGP(info, "{:*^80}", fmt::format(" REPORT {} ", name));
    for (const auto& res : vres) {
      res.get().print(!gOpts.disableSystemInfo, gOpts.debug);
    }
  }
  LOGP(info, "{:~^100}", " Reporting END ");
}

void graphResults()
{
  LOGP(info, "Graphing Results");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kOcean);

  std::unique_ptr<TFile> file;
  if (!gOpts.disableRoot) {
    file.reset(TFile::Open(gOpts.output.c_str(), "RECREATE"));
  }

  // Summary plots
  auto hsUnTot = new THStack("hs_untot", "Uncompressed Total;;MB");
  auto hsCompTot = new THStack("hs_comptot", "Compressed Total;;MB");
  auto hsUnThrough = new THStack("hs_unthrough", "Uncompressed Throughput (higher is better);;MB/s");
  auto hsCompThrough = new THStack("hs_compthrough", "Compressed Throughput;;MB/s");
  auto hsCpuEff = new THStack("hs_cpueff", "CPU Efficiency (higher is better);;Efficiency * 100");
  auto hsRatio = new THStack("hs_Ratio", "Compression Ratio (higher is better);;Ratio");
  auto leg = new TLegend(0.6, 0.6, 0.9, 0.9);

  Int_t color = 1;

  TVirtualPad* pad;
  for (const auto& [name, vres] : gRefMap) {
    LOGP(info, "  - {} : {}", name, vres.front().get().getName());
    auto hUnTot = new TH1F(Form("untot%s", name.c_str()), Form("Uncompressed Total;;MB"), vres.size(), 0, vres.size());
    auto hCompTot = new TH1F(Form("comptot%s", name.c_str()), Form("Compressed Total;;MB"), vres.size(), 0, vres.size());
    auto hUnThrough = new TH1F(Form("unthroughput%s", name.c_str()), Form("Uncompressed Throughput (higher is better);;Throughput MB/s"), vres.size(), 0, vres.size());
    auto hCompThrough = new TH1F(Form("compthroughput%s", name.c_str()), Form("Compressed Throughput;;Throughput MB/s"), vres.size(), 0, vres.size());
    auto hCpuEff = new TH1F(Form("cpuEff%s", name.c_str()), Form("CPU Efficiency (higher is better);;Efficiency  * 100"), vres.size(), 0, vres.size());
    auto hRatio = new TH1F(Form("ratio%s", name.c_str()), Form("Compression Ratio (higher is better);;Ratio"), vres.size(), 0, vres.size());

    for (size_t i{0}; i < vres.size(); ++i) {
      const auto& res = vres[i].get();
      auto branchName = res.fName.substr(4);

      hUnTot->GetXaxis()->SetBinLabel(i + 1, branchName.c_str());
      hUnTot->SetBinContent(i + 1, res.fUncompressedBytesReadTot / MB);

      hCompTot->GetXaxis()->SetBinLabel(i + 1, branchName.c_str());
      hCompTot->SetBinContent(i + 1, res.fCompressedBytesReadTot / MB);

      hUnThrough->GetXaxis()->SetBinLabel(i + 1, branchName.c_str());
      hUnThrough->SetBinContent(i + 1, res.fThroughputUncomp);

      hCompThrough->GetXaxis()->SetBinLabel(i + 1, branchName.c_str());
      hCompThrough->SetBinContent(i + 1, res.fThroughputComp);

      hCpuEff->GetXaxis()->SetBinLabel(i + 1, branchName.c_str());
      hCpuEff->SetBinContent(i + 1, res.fCpuEff);

      hRatio->GetXaxis()->SetBinLabel(i + 1, branchName.c_str());
      hRatio->SetBinContent(i + 1, res.fCompRatio);
    }

    hUnTot->SetFillColor(color);
    hUnTot->SetLineColor(color);
    hCompTot->SetFillColor(color);
    hCompTot->SetLineColor(color);
    hUnThrough->SetFillColor(color);
    hUnThrough->SetLineColor(color);
    hCompThrough->SetFillColor(color);
    hCompThrough->SetLineColor(color);
    hCpuEff->SetFillColor(color);
    hCpuEff->SetLineColor(color);
    hRatio->SetFillColor(color);
    hRatio->SetLineColor(color);

    auto c = new TCanvas(name.c_str());
    c->Divide(3, 2);
    pad = c->cd(1);
    hUnThrough->Draw("hist");
    pad->SetGrid();
    pad = c->cd(2);
    hCompThrough->Draw("hist");
    pad->SetGrid();
    pad = c->cd(3);
    hCpuEff->Draw("hist");
    pad->SetGrid();
    pad = c->cd(4);
    hUnTot->Draw("hist");
    pad->SetGrid();
    pad->SetLogy();
    pad = c->cd(5);
    hCompTot->Draw("hist");
    pad->SetGrid();
    pad->SetLogy();
    pad = c->cd(6);
    hRatio->Draw("hist");
    pad->SetGrid();
    c->Write();

    ++color;

    if (!gOpts.disablePDF) {
      c->SaveAs(Form("%s%s.pdf", gOpts.outputSuffix.c_str(), vres.front().get().getName().c_str()));
    }

    if (gRefMap.size() < 2) {
      continue;
    }

    if (!gOpts.disableSummaryPlot) {
      hsUnTot->Add(hUnTot);
      hsCompTot->Add(hCompTot);
      hsUnThrough->Add(hUnThrough);
      hsCompThrough->Add(hCompThrough);
      hsCpuEff->Add(hCpuEff);
      hsRatio->Add(hRatio);
    }

    leg->AddEntry(hUnTot, vres.front().get().getName().c_str(), "f");
  }

  if (gRefMap.size() < 2) {
    LOGP(info, " - Skipping comparison plots");
    return;
  }

  if (!gOpts.disableSummaryPlot) {
    auto csum = new TCanvas("summary");
    csum->Divide(3, 2);
    pad = csum->cd(1);
    hsUnThrough->Draw("nostackb");
    pad->SetGrid();
    pad = csum->cd(2);
    hsCompThrough->Draw("nostackb");
    pad->SetGrid();
    pad = csum->cd(3);
    hsCpuEff->Draw("nostackb");
    pad->SetGrid();
    pad = csum->cd(4);
    hsUnTot->Draw("nostackb");
    pad->SetGrid();
    pad->SetLogy();
    leg->Draw();
    pad = csum->cd(5);
    hsCompTot->Draw("nostackb");
    pad->SetGrid();
    pad->SetLogy();
    pad = csum->cd(6);
    hsRatio->Draw("nostackb");
    pad->SetGrid();
    csum->Write();

    if (!gOpts.disablePDF) {
      csum->SaveAs(Form("%saod-readspeed-summary.pdf", gOpts.outputSuffix.c_str()));
    }
  }

  if (!gOpts.ratioDefault.empty()) {
    if (gRefMap.find(gOpts.ratioDefault) == gRefMap.end()) {
      LOGP(error, "Did not find {} key for default ratio", gOpts.ratioDefault);
    } else {
      const auto& def = gRefMap[gOpts.ratioDefault];
      Int_t color = 1;
      auto hsUnTotRatio = new THStack("hs_untot_ratio", Form("Uncompressed Total (Basline %s);;Ratio", gOpts.ratioDefault.c_str()));
      auto hsCompTotRatio = new THStack("hs_comptot_ratio", Form("Compressed Total (Basline %s);;Ratio", gOpts.ratioDefault.c_str()));
      auto hsUnThroughRatio = new THStack("hs_unthrough_ratio", Form("Uncompressed Throughput (higher is better) (Basline %s);;Ratio", gOpts.ratioDefault.c_str()));
      auto hsCompThroughRatio = new THStack("hs_compthrough_ratio", Form("Compressed Throughput (higher is better) (Basline %s);;Ratio", gOpts.ratioDefault.c_str()));
      auto hsCpuEffRatio = new THStack("hs_cpueff_ratio", Form("CPU Efficiency (higher is better) (Basline %s);;Ratio", gOpts.ratioDefault.c_str()));
      auto hsRatioRatio = new THStack("hs_Ratio_ratio", Form("Compression Ratio (higher is better) (Basline %s);;Ratio", gOpts.ratioDefault.c_str()));
      auto legRatio = new TLegend(0.1, 0.1, 0.4, 0.4);
      for (const auto& [name, vres] : gRefMap) {
        if (name == gOpts.ratioDefault) {
          continue;
        }

        auto hUnTot = new TH1F(Form("untotrat%s", name.c_str()), Form("Uncompressed Total;;Ratio"), vres.size(), 0, vres.size());
        auto hCompTot = new TH1F(Form("comptotrat%s", name.c_str()), Form("Compressed Total;;Ratio"), vres.size(), 0, vres.size());
        auto hUnThrough = new TH1F(Form("unthroughputrat%s", name.c_str()), Form("Uncompressed Throughput (higher is better);;Ratio"), vres.size(), 0, vres.size());
        auto hCompThrough = new TH1F(Form("compthroughputrat%s", name.c_str()), Form("Compressed Throughput (higher is better);;Ratio"), vres.size(), 0, vres.size());
        auto hCpuEff = new TH1F(Form("cpuEffrat%s", name.c_str()), Form("CPU Efficiency (higher is better);;Ratio"), vres.size(), 0, vres.size());
        auto hRatio = new TH1F(Form("ratiorat%s", name.c_str()), Form("Compression Ratio (higher is better);;Ratio"), vres.size(), 0, vres.size());

        for (size_t i{0}; i < vres.size(); ++i) {
          const auto& defRes = def[i].get();
          const auto& res = vres[i].get();
          auto branchName = res.fName.substr(4);

          hUnTot->GetXaxis()->SetBinLabel(i + 1, branchName.c_str());
          hUnTot->SetBinContent(i + 1, static_cast<Double_t>(res.fUncompressedBytesReadTot) / static_cast<Double_t>(defRes.fUncompressedBytesReadTot));

          hCompTot->GetXaxis()->SetBinLabel(i + 1, branchName.c_str());
          hCompTot->SetBinContent(i + 1, static_cast<Double_t>(res.fCompressedBytesReadTot) / static_cast<Double_t>(defRes.fCompressedBytesReadTot));

          hUnThrough->GetXaxis()->SetBinLabel(i + 1, branchName.c_str());
          hUnThrough->SetBinContent(i + 1, static_cast<Double_t>(res.fThroughputUncomp) / static_cast<Double_t>(defRes.fThroughputUncomp));

          hCompThrough->GetXaxis()->SetBinLabel(i + 1, branchName.c_str());
          hCompThrough->SetBinContent(i + 1, static_cast<Double_t>(res.fThroughputComp) / static_cast<Double_t>(defRes.fThroughputUncomp));

          hCpuEff->GetXaxis()->SetBinLabel(i + 1, branchName.c_str());
          hCpuEff->SetBinContent(i + 1, static_cast<Double_t>(res.fCpuEff) / static_cast<Double_t>(defRes.fCpuEff));

          hRatio->GetXaxis()->SetBinLabel(i + 1, branchName.c_str());
          hRatio->SetBinContent(i + 1, static_cast<Double_t>(res.fCompRatio) / static_cast<Double_t>(defRes.fCompRatio));
        }

        hUnTot->SetFillColor(color);
        hUnTot->SetLineColor(color);
        hCompTot->SetFillColor(color);
        hCompTot->SetLineColor(color);
        hUnThrough->SetFillColor(color);
        hUnThrough->SetLineColor(color);
        hCompThrough->SetFillColor(color);
        hCompThrough->SetLineColor(color);
        hCpuEff->SetFillColor(color);
        hCpuEff->SetLineColor(color);
        hRatio->SetFillColor(color);
        hRatio->SetLineColor(color);

        ++color;

        hsUnTotRatio->Add(hUnTot);
        hsCompTotRatio->Add(hCompTot);
        hsUnThroughRatio->Add(hUnThrough);
        hsCompThroughRatio->Add(hCompThrough);
        hsCpuEffRatio->Add(hCpuEff);
        hsRatioRatio->Add(hRatio);

        legRatio->AddEntry(hUnTot, vres.front().get().getName().c_str(), "f");
      }

      auto cratio = new TCanvas("ratio");
      cratio->Divide(3, 2);
      pad = cratio->cd(1);
      hsUnThroughRatio->Draw("nostackb");
      pad->SetGrid();
      pad = cratio->cd(2);
      hsCompThroughRatio->Draw("nostackb");
      pad->SetGrid();
      pad = cratio->cd(3);
      hsCpuEffRatio->Draw("nostackb");
      pad->SetGrid();
      pad = cratio->cd(4);
      hsUnTotRatio->Draw("nostackb");
      pad->SetGrid();
      legRatio->Draw();
      pad = cratio->cd(5);
      hsCompTotRatio->Draw("nostackb");
      pad->SetGrid();
      pad = cratio->cd(6);
      hsRatioRatio->Draw("nostackb");
      pad->SetGrid();
      cratio->Write();

      if (!gOpts.disablePDF) {
        cratio->SaveAs(Form("%saod-readspeed-ratio.pdf", gOpts.outputSuffix.c_str()));
      }
      LOGP(info, " - Ratio plot");
    }
  }
}

void cleanResults()
{
  LOGP(info, "Cleaning results");

  unsigned truncated{0};
  if (!gOpts.disableEmptyResultTruncation) { // skip empty reads
    for (auto it = gResultMap.begin(); it != gResultMap.end();) {
      if (it->second.fUncompressedBytesReadTot == 0) {
        ++truncated;
        it = gResultMap.erase(it);
      } else {
        ++it;
      }
    }
  }
  LOGP(info, "    `-> Truncated {} Results", truncated);

  unsigned suppressed{0};
  if (!gOpts.disableResultSupression) {
    for (auto it = gResultMap.begin(); it != gResultMap.end();) {
      if (!it->second.fAccountable) {
        ++suppressed;
        it = gResultMap.erase(it);
      } else {
        ++it;
      }
    }
  }
  LOGP(info, "    `-> Suppressed {} Results", suppressed);

  LOGP(info, "Keeping {} Results", gResultMap.size());
}

void makeReferences()
{
  LOGP(info, "Making references");
  for (const auto& [name, res] : gResultMap) {
    const auto algoName = name.substr(0, 3);
    if (gRefMap.find(algoName) != gRefMap.end()) {
      gRefMap[algoName].emplace_back(res);
    } else {
      gRefMap[algoName] = {res};
    }
  }
  // Sort references using string compare
  for (auto& [algoName, vres] : gRefMap) {
    std::sort(vres.begin(), vres.end(), [](auto a, auto b) {
      return a.get().fName < b.get().fName;
    });
  }
}

void doAODSummary()
{
  LOGP(info, "AO2D summary for:");
  for (auto& [algoName, vres] : gRefMap) {
    LOGP(info, " - {}", algoName);
    Result res;
    for (const auto& r : vres) {
      res += r;
    }
    res.fCounter /= vres.size();
    res.fName = algoName + ".AO2D";
    res.setAlgorithm();
    res.setLevel();
    res.normalize();
    gResultMap[algoName] = res;
    vres.insert(vres.begin(), std::cref(gResultMap[algoName]));
  }
}

int main(int argc, char* argv[])
{
  const char* desc{
    "AO2D Read Speed tester"
    "\nThis program allows to benchmark AO2D files for commonly used observables:"
    "\n   - Un-/Compressed throughput"
    "\n   - CPU Effiency"
    "\nThe observalbles are created for every compression settings of different trees to compare compression configuration of ROOT files"};
  bpo::options_description options(desc);
  if (!initOptionsAndParse(options, argc, argv)) {
    return 1;
  }

  // Input Files
  std::vector<std::string> fileNames;
  if (!gOpts.input.empty()) {
    fileNames = tokenize(gOpts.input);
  } else if (!gOpts.inputReg.empty()) {
    RegexHandler matcher;
    matcher.init(gOpts.inputReg, "");
    for (const auto& entry : fs::directory_iterator(fs::current_path())) {
      if (fs::is_regular_file(entry.path()) && matcher.matchesAny(entry.path().filename())) {
        fileNames.emplace_back(entry.path().filename());
      }
    }
  }

  // Quick existence check
  fileNames.erase(std::remove_if(fileNames.begin(), fileNames.end(), [](const auto& file) {
                    const bool exists = fs::exists(file);
                    LOG_IF(info, gOpts.verbose) << "Checking for " << file << " -> " << exists;
                    LOG_IF(error, !exists) << "Removing " << file << " from fileNames does not exists!";
                    return !exists;
                  }),
                  fileNames.end());
  LOGP(info, "FileList: ");
  for (size_t i{0}; i < fileNames.size(); ++i) {
    LOGP(info, "   {}. {}", i, fileNames[i]);
  }
  if (fileNames.empty()) {
    LOGP(error, "No viable input files provided!");
  }

  if (!gOpts.treeIncReg.empty() || !gOpts.treeExReg.empty()) {
    LOGP(info, "Building tree regex");
    gTreeRegex.init(gOpts.treeIncReg, gOpts.treeExReg);
  }

  if (!gOpts.branchIncReg.empty() || !gOpts.branchExReg.empty()) {
    LOGP(info, "Building branch regex");
    gBranchRegex.init(gOpts.branchIncReg, gOpts.branchExReg);
  }

  // Make the mapping for tree to file and dfs
  std::for_each(fileNames.cbegin(), fileNames.cend(), makeTreeMap);
  finalizeTreeMap(gTreeMap);
  if (gOpts.verbose) {
    printTreeMap(gTreeMap);
  }
  produceResults();
  for (auto& [_, res] : gResultMap) {
    res.normalize();
  }
  cleanResults();
  makeReferences();

  if (!gOpts.disableSummary) {
    LOG_IF(info, gOpts.debug) << "Reporting results";
    if (!gOpts.disableAODSummary) {
      doAODSummary();
    }
    reportResults();
  }

  if (!gOpts.disableGraphs) {
    gROOT->SetBatch(gOpts.batch);
    graphResults();
  }

  return 0;
}
