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

#ifndef AODREADSPEED_H_
#define AODREADSPEED_H_

#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <numeric>
#include <regex>
#include <cmath>

#include "Framework/Logger.h"

#include "RZip.h"
#include "TSystem.h"

constexpr Double_t KB = 1024.;
constexpr Double_t MB = KB * KB;

struct ProgramOptions {
  // General
  std::string input;
  std::string inputReg;
  std::string output;
  std::string outputSuffix;
  unsigned int repeat;
  bool batch;
  bool verbose;
  bool debug;

  // Results
  std::string ratioDefault;
  bool disablePDF;
  bool disableRoot;
  bool disableGraphs;
  bool disableSystemInfo;
  bool disableSummary;
  bool disableSummaryPlot;
  bool disableAODSummary;
  bool disableResultSupression;
  bool disableEmptyResultTruncation;

  // File
  bool disableFilePrefetch;
  Int_t fileReadAheadSize;

  // Branch
  bool disableAllBranches;
  std::string branchIncReg;
  std::string branchExReg;

  // Tree
  std::string treeIncReg;
  std::string treeExReg;
  Int_t treeCacheSize;
};

// Remove leading whitespace
inline std::string ltrimSpace(std::string src)
{
  return src.erase(0, src.find_first_not_of(' '));
}

// Remove trailing whitespace
inline std::string rtrimSpace(std::string src)
{
  return src.erase(src.find_last_not_of(' ') + 1);
}

// Remove leading/trailing whitespace
inline std::string trimSpace(std::string const& src)
{
  return ltrimSpace(rtrimSpace(src));
}

// Split a given string on a delim character, return vector of tokens
// If trim is true, then also remove leading/trailing whitespace of each token.
inline std::vector<std::string> tokenize(const std::string& src, char delim = ',', bool trim = true)
{
  std::stringstream ss(src);
  std::string token;
  std::vector<std::string> tokens;

  while (std::getline(ss, token, delim)) {
    token = (trim ? trimSpace(token) : token);
    if (!token.empty()) {
      tokens.push_back(std::move(token));
    }
  }

  return tokens;
}

// Simple regex handler
struct RegexHandler {
  void init(const std::string& regsInc, const std::string& regsEx)
  {
    if (!regsInc.empty()) {
      for (const auto& sreg : tokenize(regsInc)) {
        regInc.emplace_back(sreg);
      }
    }

    if (!regsEx.empty()) {
      for (const auto& sreg : tokenize(regsEx)) {
        regEx.emplace_back(sreg);
      }
    }
  }

  bool matchesAny(const std::string name)
  {
    bool isInc{false}, isEx{false};

    for (const auto& reg : regInc) {
      if (std::regex_match(name, reg)) {
        isInc = true;
        break;
      }
    }

    for (const auto& reg : regEx) {
      if (std::regex_match(name, reg)) {
        isEx = true;
        break;
      }
    }

    if (regInc.empty() && !regEx.empty()) {
      return !isEx;
    } else if (!regInc.empty() && regEx.empty()) {
      return isInc;
    } else if (regInc.empty() && regEx.empty()) {
      return true;
    } else {
      return isInc && !isEx;
    }
  }

  std::vector<std::regex> regInc; // inclusive
  std::vector<std::regex> regEx;  // inclusive
};

struct Result {
  Double_t fRealTime{0.};                 // [s]
  Double_t fRealTimeAvg{0.};              // [s]
  Double_t fCpuTime{0.};                  // [s]
  Double_t fCpuTimeAvg{0.};               // [s]
  Double_t fCpuEff{0.};                   // [eff]
  bool fCpuEffEstimated{false};           // too short
  ULong64_t fUncompressedBytesReadTot{0}; // [B]
  Double_t fUncompressedBytesReadAvg;     // [MB]
  Double_t fThroughputUncomp;             // [MB/s]
  ULong64_t fCompressedBytesReadTot{0};   // [B]
  Double_t fCompressedBytesReadAvg;       // [MB]
  Double_t fThroughputComp;               // [MB/s]
  Double_t fCompRatio{0.};                // [ratio]
  bool fCompRatioEstimated{false};        // no compression
  bool fAccountable{true};
  unsigned int fCounter{0};
  unsigned int fRepeat{0};
  float fCpuSys{0.};
  float fCpuUser{0.};
  float fMemResident{0.};
  float fMemVirtual{0.};
  int fAlgorithm;
  int fLevel;
  ProcInfo_t fBefore;
  ProcInfo_t fAfter;
  std::string fName;
  bool fNormalized{false};

  Result& operator+=(const Result& other)
  {
    this->fRealTime += other.fRealTime;
    this->fCpuTime += other.fCpuTime;
    this->fUncompressedBytesReadTot += other.fUncompressedBytesReadTot;
    this->fCompressedBytesReadTot += other.fCompressedBytesReadTot;
    this->fCounter += other.fCounter;

    this->fCpuSys = std::max(this->fAfter.fCpuSys, fBefore.fCpuSys);
    this->fCpuUser = std::max(this->fAfter.fCpuUser, fBefore.fCpuUser);
    this->fMemResident = std::max(this->fAfter.fMemResident, fBefore.fMemResident);
    this->fMemVirtual = std::max(this->fAfter.fMemVirtual, fBefore.fMemVirtual);
    return *this;
  }

  void normalize()
  {
    if (fNormalized) {
      return;
    }

    fRealTimeAvg = fRealTime / static_cast<Double_t>(fCounter);
    fCpuTimeAvg = fCpuTime / static_cast<Double_t>(fCounter);
    fUncompressedBytesReadAvg = static_cast<Double_t>(fUncompressedBytesReadTot) / static_cast<Double_t>(fCounter) / MB;
    fCompressedBytesReadAvg = static_cast<Double_t>(fCompressedBytesReadTot) / static_cast<Double_t>(fCounter) / MB;

    fThroughputUncomp = fUncompressedBytesReadAvg / fRealTimeAvg;
    fThroughputComp = fCompressedBytesReadAvg / fRealTimeAvg;

    fCpuEff = (fCpuTimeAvg / fRealTimeAvg) * 100.;
    if (fCpuTime == 0. || fCpuTimeAvg > fRealTimeAvg) {
      fCpuEffEstimated = true;
      fCpuEff = 100.;
    }

    fCompRatio = static_cast<Double_t>(fUncompressedBytesReadTot) / static_cast<Double_t>(fCompressedBytesReadTot);
    if (fCompressedBytesReadTot == fUncompressedBytesReadTot) {
      fCompRatioEstimated = true;
      fCompRatio = 0;
    } else if (fCompressedBytesReadTot > fUncompressedBytesReadTot) {
      fAccountable = false;
    }

    fNormalized = true;
  }

  void print(bool systemInfo = true, bool debug = false) const
  {
    auto getPrecision = [](auto value, unsigned int minimum = 1) -> unsigned int {
      if (value == static_cast<decltype(value)>(0)) {
        return 1;
      }
      const float abs = std::abs(value);
      auto log = std::log10(abs);
      unsigned int val = static_cast<int>(-std::floor(log)) + 1;
      if (val > 10) { // max precision
        val = 5;
      } else if (val == 1) {
        ++val;
      }
      return std::max(val, minimum);
    };

    LOGP(info, "{:*>30}: {} -> {} DFs with {} repeats", "", fName, fCounter, fRepeat);
    LOGP(info, "Real time: {:.{}f} s/DF; total: {:.{}f} s", fRealTimeAvg, getPrecision(fRealTimeAvg), fRealTime, getPrecision(fRealTime));
    LOGP(info, "CPU time: {:.{}f} s/DF; total: {:.{}f} s", fCpuTimeAvg, getPrecision(fCpuTimeAvg), fCpuTime, getPrecision(fCpuTime));
    LOGP(info, "Uncompressed Bytes read: {}\t({} MB)", fUncompressedBytesReadTot, static_cast<ULong64_t>(fUncompressedBytesReadTot / MB));
    LOGP(info, "  Compressed Bytes read: {}\t({} MB)\tRatio={}", fCompressedBytesReadTot, static_cast<ULong64_t>(fCompressedBytesReadTot / MB),
         (fCompRatioEstimated) ? "no compression" : fmt::format("{:.{}f}", fCompRatio, getPrecision(fCompRatio, 3)));
    LOGP(info, "Uncompressed throughput: {:.{}f} MB/s", fThroughputUncomp, getPrecision(fThroughputUncomp));
    LOGP(info, "  Compressed throughput: {:.{}f} MB/s", fThroughputComp, getPrecision(fThroughputComp));
    LOGP(info, "CPU efficiency: {:.2f}% {}", fCpuEff, (fCpuEffEstimated) ? "*estimated*" : "");
    if (systemInfo) {
      LOGP(info, "CPU Sys={:.{}f}; CPU User={:.{}f}", fCpuSys, getPrecision(fCpuSys), fCpuUser, getPrecision(fCpuUser));
      LOGP(info, "MemResident={:.{}f} MB; MemVirutal={:.{}f} MB", fMemResident / MB, getPrecision(fMemResident / MB), fMemVirtual, getPrecision(fMemVirtual / MB));
    }
    LOGP(info, "Algorithm: {} with Level: {}", getAlgorithm(), getLevel(false));
    LOG_IF(info, !fAccountable) << "--> will not be further accounted for!";
  }

  std::string getAlgorithm() const
  {
    switch (fAlgorithm) {
      case ROOT::RCompressionSetting::EAlgorithm::kInherit:
        return "inherited";
      case ROOT::RCompressionSetting::EAlgorithm::kUseGlobal:
        return "global";
      case ROOT::RCompressionSetting::EAlgorithm::kZLIB:
        return "ZLIB";
      case ROOT::RCompressionSetting::EAlgorithm::kLZMA:
        return "LZMA";
      case ROOT::RCompressionSetting::EAlgorithm::kOldCompressionAlgo:
        return "Old";
      case ROOT::RCompressionSetting::EAlgorithm::kLZ4:
        return "LZ4";
      case ROOT::RCompressionSetting::EAlgorithm::kZSTD:
        return "ZSTD";
      case ROOT::RCompressionSetting::EAlgorithm::kUndefined:
      default:
        return "undefined";
    }
  }

  std::string getLevel(bool symbName = true) const
  {
    if (symbName) {
      switch (fLevel) {
        case ROOT::RCompressionSetting::ELevel::kInherit:
          return "inherited";
        case ROOT::RCompressionSetting::ELevel::kUncompressed:
          return "uncompressed";
        case ROOT::RCompressionSetting::ELevel::kUseMin:
          return "defaultZLIB";
        case ROOT::RCompressionSetting::ELevel::kDefaultLZ4:
          return "defaultZLZ4";
        case ROOT::RCompressionSetting::ELevel::kDefaultZSTD:
          return "defaultZSTD";
        case ROOT::RCompressionSetting::ELevel::kDefaultOld:
          return "defaultOld";
        case ROOT::RCompressionSetting::ELevel::kDefaultLZMA:
          return "defaultLZMA";
        default:
          return "undefined";
      }
    } else {
      return std::to_string(fLevel);
    }
  }

  std::string getName() const
  {
    return getAlgorithm() + "-" + getLevel(false);
  }

  void setAlgorithm()
  {
    fAlgorithm = std::stoi(fName.substr(0, 1));
  }

  void setLevel()
  {
    fLevel = std::stoi(fName.substr(2, 1));
  }
};

using Results = std::unordered_map<std::string, Result>;

struct TreeMapping {
  std::string fTreeName;
  std::unordered_map<std::string, std::vector<std::string>> fFileMap; // Maps file -> DF

  void addEntry(const std::string& file, const std::string& df)
  {
    if (fFileMap.find(file) != fFileMap.end()) {
      fFileMap[file].emplace_back(df);
    } else {
      fFileMap[file] = {df};
    }
  }

  void print() const
  {
    LOGP(info, "TreeMapping for {}", fTreeName);
    for (const auto& [file, dfs] : fFileMap) {
      LOGP(info, "   File: {}", file);
      for (const auto& df : dfs) {
        LOGP(info, "      - {}", df);
      }
    }
  }

  void finalize()
  {
    for (auto& [file, dfs] : fFileMap) {
      std::sort(dfs.begin(), dfs.end(), [](auto& a, auto& b) {
        // Order in timestamp
        const ULong64_t at = std::stoull(a.substr(3)), bt = std::stoull(b.substr(3));
        return a < b;
      });
    }
  }
};

using TreeMap = std::unordered_map<std::string, TreeMapping>;

inline void finalizeTreeMap(TreeMap& map)
{
  for (auto& [_, mapping] : map) {
    mapping.finalize();
  }
}

inline void printTreeMap(const TreeMap& map)
{
  for (const auto& [_, mapping] : map) {
    mapping.print();
  }
}

#endif // AODREADSPEED_H_
