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

#include <unordered_map>
#include <numeric>
#include <algorithm>
#include <getopt.h>

#include "TSystem.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TRegexp.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TList.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TObjString.h"
#include "TGrid.h"
#include "TMap.h"
#include "TLeaf.h"

#include "aodMerger.h"

// AOD reduction tool
//   Designed for the 2022 pp data with specific selections:
//   - Remove all TPC only tracks, keep TPC-only V0 tracks
//   - Remove all V0s which refer to any removed track
//   - Remove all cascade which refer to any removed V0
//   - Remove all ambiguous track entries which point to a track with collision
//   - Adjust all indices
//   adjusted also for 2023

// Collection of possible error codes.
// Some of these are transient, e.g., overwritten by others.
// If any is encountered check the logs!
enum ErrorCodes : int {
  SUCCESS = 0,                   // Succesful exit
  ORDERED_KEYS = 5,              // Encounter unoreded keys
  MISSING_TRKEXTRA = 6,          // Branch O2trackextra* was not present in file
  MISSING_TRKIU = 7,             // Branch O2trackiu* was not present in file
  MISSING_V0 = 8,                // Branch O2v0* was not present in file
  NO_VLA_SUPPORT = 9,            // Branch with VLA detected but is not supported by thinner
  BRANCH_DETECTION = 10,         // Sanity check for branches failed
  EMPTY_FILE = 11,               // Input/Output file empty after thinning
  MISSING_TRACKEDV0 = 12,        // Branch O2trackedv0* not present in file
  MISSING_TRACKEDCASC = 13,      // Branch O2trackedcascade* not present in file
  MISSING_TRACKED3BODY = 14,     // Branch O2tracked3body* not present in file
  SANITY_TRACKED = 15,           // Late sanity check if O2tracked* is present
  NOT_EXPECTED_REDUCTION = 16,   // Reduced a tree which should not be reduced
  MISSING_CASCADES = 17,         // Branch O2cascades not preset in file
  MISSING_DECAY3BODYS = 18,      // Branch O2decay3body not preset in file
  MISSING_TRACKQA = 19,          // Branch O2trackqa not preset in file
  PROCESSING_MORE_THAN_ONE = 20, // Processing more than once actively check branch
};

int main(int argc, char* argv[])
{
  std::string inputFileName("AO2D.root");
  std::string outputFileName("AO2D_thinned.root");
  int exitCode = SUCCESS; // 0: success, !=0: failure
  bool bOverwrite = false;

  int option_index = 1;

  const char* const short_opts = "i:o:KOh";
  static struct option long_options[] = {
    {"input", required_argument, nullptr, 'i'},
    {"output", required_argument, nullptr, 'o'},
    {"overwrite", no_argument, nullptr, 'O'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}};

  while (true) {
    const auto opt = getopt_long(argc, argv, short_opts, long_options, &option_index);
    if (opt == -1) {
      break; // use defaults
    }
    switch (opt) {
      case 'i':
        inputFileName = optarg;
        break;
      case 'o':
        outputFileName = optarg;
        break;
      case 'O':
        bOverwrite = true;
        printf("Overwriting existing output file if existing\n");
        break;
      case 'h':
      case '?':
      default:
        printf("AO2D thinning tool. Options: \n");
        printf("  --input/-i <inputfile.root>     Contains input file path to the file to be thinned. Default: %s\n", inputFileName.c_str());
        printf("  --output/-o <outputfile.root>   Target output ROOT file. Default: %s\n", outputFileName.c_str());
        printf("\n");
        printf("  Optional Arguments:\n");
        printf("  --overwrite/-O                  Overwrite existing output file\n");
        return -1;
    }
  }

  printf("AOD reduction started with:\n");
  printf("  Input file: %s\n", inputFileName.c_str());
  printf("  Ouput file name: %s\n", outputFileName.c_str());

  TStopwatch clock;
  clock.Start(kTRUE);

  auto outputFile = TFile::Open(outputFileName.c_str(), (bOverwrite) ? "RECREATE" : "CREATE", "", 505);
  if (outputFile == nullptr) {
    printf("Error: File %s exists or cannot be created!\n", outputFileName.c_str());
    return 1;
  }
  TDirectory* outputDir = nullptr;

  if (inputFileName.find("alien:") == 0) {
    printf("Connecting to AliEn...");
    TGrid::Connect("alien:");
  }

  auto inputFile = TFile::Open(inputFileName.c_str());
  if (!inputFile) {
    printf("Error: Could not open input file %s.\n", inputFileName.c_str());
    return 1;
  }

  TList* keyList = inputFile->GetListOfKeys();
  keyList->Sort();

  for (auto key1 : *keyList) {
    // Keep metaData
    if (((TObjString*)key1)->GetString().EqualTo("metaData")) {
      auto metaData = (TMap*)inputFile->Get("metaData");
      outputFile->cd();
      metaData->Write("metaData", TObject::kSingleKey);
    }

    // Keep parentFiles
    if (((TObjString*)key1)->GetString().EqualTo("parentFiles")) {
      auto parentFiles = (TMap*)inputFile->Get("parentFiles");
      outputFile->cd();
      parentFiles->Write("parentFiles", TObject::kSingleKey);
    }

    // Skip everything else, except 'DF_*'
    if (!((TObjString*)key1)->GetString().BeginsWith("DF_")) {
      continue;
    }

    auto dfName = ((TObjString*)key1)->GetString().Data();

    printf("  Processing folder %s\n", dfName);
    auto folder = (TDirectoryFile*)inputFile->Get(dfName);
    auto treeList = folder->GetListOfKeys();
    treeList->Sort();

    // purging keys from duplicates
    for (auto i = 0; i < treeList->GetEntries(); ++i) {
      TKey* ki = (TKey*)treeList->At(i);
      for (int j = i + 1; j < treeList->GetEntries(); ++j) {
        TKey* kj = (TKey*)treeList->At(j);
        if (std::strcmp(ki->GetName(), kj->GetName()) == 0 && std::strcmp(ki->GetTitle(), kj->GetTitle()) == 0) {
          if (ki->GetCycle() < kj->GetCycle()) {
            printf("    *** FATAL *** we had ordered the keys, first cycle should be higher, please check");
            exitCode = ORDERED_KEYS;
          } else {
            // key is a duplicate, let's remove it
            treeList->Remove(kj);
            j--;
          }
        } else {
          // we changed key, since they are sorted, we won't have the same anymore
          break;
        }
      }
    }

    // Scan versions e.g. 001 or 002 ...
    TString v0Name{"O2v0*"};
    TRegexp v0Re(v0Name, kTRUE);
    TString cascName{"O2cascade*"};
    TRegexp cascRe(cascName, kTRUE);
    TString decay3bodyName{"O2decay3body*"};
    TRegexp decay3bodyRe(decay3bodyName, kTRUE);
    TString trkExtraName{"O2trackextra*"};
    TRegexp trkExtraRe(trkExtraName, kTRUE);
    bool hasTrackedV0{false};
    TString trackedV0Name{"O2trackedv0*"};
    TRegexp trackedV0Re(trackedV0Name, kTRUE);
    bool hasTrackedCasc{false};
    TString trackedCascName{"O2trackedcascade*"};
    TRegexp trackedCascRe(trackedCascName, kTRUE);
    bool hasTracked3Body{false};
    TString tracked3BodyName{"O2tracked3body*"};
    TRegexp tracked3BodyRe(tracked3BodyName, kTRUE);
    bool hasTrackQA{false};
    TString trackQAName{"O2trackqa*"};
    TRegexp trackQARe(trackQAName, kTRUE);
    for (TObject* obj : *treeList) {
      TString st = obj->GetName();
      if (st.Index(v0Re) != kNPOS) {
        v0Name = st;
      } else if (st.Index(cascRe) != kNPOS) {
        cascName = st;
      } else if (st.Index(decay3bodyRe) != kNPOS) {
        decay3bodyName = st;
      } else if (st.Index(trkExtraRe) != kNPOS) {
        trkExtraName = st;
      } else if (st.Index(trackedV0Re) != kNPOS) {
        hasTrackedV0 = true;
        trackedV0Name = st;
      } else if (st.Index(trackedCascRe) != kNPOS) {
        hasTrackedCasc = true;
        trackedCascName = st;
      } else if (st.Index(tracked3BodyRe) != kNPOS) {
        hasTracked3Body = true;
        tracked3BodyName = st;
      } else if (st.Index(trackQARe) != kNPOS) {
        hasTrackQA = true;
        trackQAName = st;
      }
    }

    // Certain order needed in order to populate vectors of skipped entries
    auto v0Entry = (TObject*)treeList->FindObject(v0Name);
    treeList->Remove(v0Entry);
    treeList->AddFirst(v0Entry);

    // Prepare maps for track skimming
    auto trackExtraTree = (TTree*)inputFile->Get(Form("%s/%s", dfName, trkExtraName.Data()));
    if (trackExtraTree == nullptr) {
      printf("%s table not found\n", trkExtraName.Data());
      exitCode = MISSING_TRKEXTRA;
      break;
    }
    auto track_iu = (TTree*)inputFile->Get(Form("%s/%s", dfName, "O2track_iu"));
    if (track_iu == nullptr) {
      printf("O2track_iu table not found\n");
      exitCode = MISSING_TRKIU;
      break;
    }
    auto v0s = (TTree*)inputFile->Get(Form("%s/%s", dfName, v0Name.Data()));
    if (v0s == nullptr) {
      printf("%s table not found\n", v0Name.Data());
      exitCode = MISSING_V0;
      break;
    }
    auto cascades = (TTree*)inputFile->Get(Form("%s/%s", dfName, cascName.Data()));
    if (cascades == nullptr) {
      printf("%s table not found\n", cascName.Data());
      exitCode = MISSING_CASCADES;
      break;
    }
    auto decay3bodys = (TTree*)inputFile->Get(Form("%s/%s", dfName, decay3bodyName.Data()));
    if (decay3bodys == nullptr) {
      printf("%s table not found\n", decay3bodyName.Data());
      exitCode = MISSING_DECAY3BODYS;
      break;
    }

    // Strangeness Tracking Trees
    TTree *trackedV0s{nullptr}, *trackedCasc{nullptr}, *tracked3Body{nullptr};
    if (hasTrackedV0 && (trackedV0s = (TTree*)inputFile->Get(Form("%s/%s", dfName, trackedV0Name.Data()))) == nullptr) {
      exitCode = MISSING_TRACKEDV0;
      break;
    }
    if (hasTrackedCasc && (trackedCasc = (TTree*)inputFile->Get(Form("%s/%s", dfName, trackedCascName.Data()))) == nullptr) {
      exitCode = MISSING_TRACKEDCASC;
      break;
    }
    if (hasTracked3Body && (tracked3Body = (TTree*)inputFile->Get(Form("%s/%s", dfName, tracked3BodyName.Data()))) == nullptr) {
      exitCode = MISSING_TRACKED3BODY;
      break;
    }
    // TrackQA
    TTree* trackQA;
    if (hasTrackQA && (trackQA = (TTree*)inputFile->Get(Form("%s/%s", dfName, trackQAName.Data()))) == nullptr) {
      exitCode = MISSING_TRACKQA;
      break;
    }

    std::vector<int> acceptedTracks(trackExtraTree->GetEntriesFast(), -1); // New track Idx for accepted Tracks
    std::vector<bool> hasCollision(trackExtraTree->GetEntriesFast(), false);

    // We need to loop over the V0s once and flag the prong indices
    int trackIdxPos = 0, trackIdxNeg = 0;
    std::vector<bool> keepV0s(trackExtraTree->GetEntriesFast(), false);
    v0s->SetBranchAddress("fIndexTracks_Pos", &trackIdxPos);
    v0s->SetBranchAddress("fIndexTracks_Neg", &trackIdxNeg);
    for (int i{0}; i < v0s->GetEntriesFast(); ++i) {
      v0s->GetEntry(i);
      keepV0s[trackIdxPos] = true;
      keepV0s[trackIdxNeg] = true;
    }

    int trackIdxCasc = -1;
    std::vector<bool> keepCascades(trackExtraTree->GetEntriesFast(), false);
    cascades->SetBranchAddress("fIndexTracks", &trackIdxCasc);
    for (int i{0}; i < cascades->GetEntriesFast(); ++i) {
      cascades->GetEntry(i);
      keepCascades[trackIdxCasc] = true;
    }

    int trackIdx0{-1}, trackIdx1{-1}, trackIdx2{-1};
    std::vector<bool> keep3Bodys(trackExtraTree->GetEntriesFast(), false);
    decay3bodys->SetBranchAddress("fIndexTracks_0", &trackIdx0);
    decay3bodys->SetBranchAddress("fIndexTracks_1", &trackIdx1);
    decay3bodys->SetBranchAddress("fIndexTracks_2", &trackIdx2);
    for (int i{0}; i < decay3bodys->GetEntriesFast(); ++i) {
      decay3bodys->GetEntry(i);
      keep3Bodys[trackIdx0] = true;
      keep3Bodys[trackIdx1] = true;
      keep3Bodys[trackIdx2] = true;
    }

    /// Strangeness Tracked
    // V0
    std::vector<bool> keepTrackedV0s;
    if (hasTrackedV0) {
      keepTrackedV0s = std::vector<bool>(trackExtraTree->GetEntriesFast(), false);
      trackedV0s->SetBranchAddress("fIndexTracks", &trackIdx0);
      trackedV0s->SetBranchAddress("fIndexTracks_ITS", &trackIdx1);
      for (int i{0}; i < trackedV0s->GetEntriesFast(); ++i) {
        trackedV0s->GetEntry(i);
        keepTrackedV0s[trackIdx0] = true;
        keepTrackedV0s[trackIdx1] = true;
      }
    }
    // Cascades
    std::vector<bool> keepTrackedCascs;
    if (hasTrackedCasc) {
      keepTrackedCascs = std::vector<bool>(trackExtraTree->GetEntriesFast(), false);
      trackedCasc->SetBranchAddress("fIndexTracks", &trackIdx0);
      trackedCasc->SetBranchAddress("fIndexTracks_ITS", &trackIdx1);
      for (int i{0}; i < trackedCasc->GetEntriesFast(); ++i) {
        trackedCasc->GetEntry(i);
        keepTrackedCascs[trackIdx0] = true;
        keepTrackedCascs[trackIdx1] = true;
      }
    }
    // 3Bodys
    std::vector<bool> keepTracked3Bodys;
    if (hasTracked3Body) {
      keepTracked3Bodys = std::vector<bool>(trackExtraTree->GetEntriesFast(), false);
      tracked3Body->SetBranchAddress("fIndexTracks", &trackIdx0);
      tracked3Body->SetBranchAddress("fIndexTracks_ITS", &trackIdx1);
      for (int i{0}; i < tracked3Body->GetEntriesFast(); ++i) {
        tracked3Body->GetEntry(i);
        keepTracked3Bodys[trackIdx0] = true;
        keepTracked3Bodys[trackIdx1] = true;
      }
    }

    // TrackQA
    std::vector<bool> keepTrackQA;
    if (hasTrackQA) {
      keepTrackQA = std::vector<bool>(trackExtraTree->GetEntriesFast(), false);
      trackQA->SetBranchAddress("fIndexTracks", &trackIdx0);
      for (int i{0}; i < trackQA->GetEntriesFast(); ++i) {
        trackQA->GetEntry(i);
        keepTrackQA[trackIdx0] = true;
      }
    }

    uint8_t tpcNClsFindable = 0;
    bool bTPClsFindable = false;
    uint8_t ITSClusterMap = 0;
    UInt_t ITSClusterSizes = 0;
    bool bITSClusterMap = false;
    bool bITSClusterSizes = false;
    uint8_t TRDPattern = 0;
    bool bTRDPattern = false;
    float_t TOFChi2 = 0;
    bool bTOFChi2 = false;

    // Test if track properties exist
    TBranch* br = nullptr;
    TIter next(trackExtraTree->GetListOfBranches());
    while ((br = (TBranch*)next())) {
      TString brName = br->GetName();
      if (brName == "fTPCNClsFindable") {
        trackExtraTree->SetBranchAddress("fTPCNClsFindable", &tpcNClsFindable);
        bTPClsFindable = true;
      } else if (brName == "fITSClusterMap") {
        trackExtraTree->SetBranchAddress("fITSClusterMap", &ITSClusterMap);
        bITSClusterMap = true;
      } else if (brName == "fITSClusterSizes") {
        trackExtraTree->SetBranchAddress("fITSClusterSizes", &ITSClusterSizes);
        bITSClusterSizes = true;
      } else if (brName == "fTRDPattern") {
        trackExtraTree->SetBranchAddress("fTRDPattern", &TRDPattern);
        bTRDPattern = true;
      } else if (brName == "fTOFChi2") {
        trackExtraTree->SetBranchAddress("fTOFChi2", &TOFChi2);
        bTOFChi2 = true;
      }
    }

    // Sanity-Check
    // If any (%ITSClusterMap or %ITSClusterSizes) of these are not found, continuation is not possible, hence fataling
    if (!bTPClsFindable || !bTRDPattern || !bTOFChi2 ||
        (!bITSClusterMap && !bITSClusterSizes)) {
      printf("    *** FATAL *** Branch detection failed in %s for trackextra.[(fITSClusterMap=%d,fITSClusterSizes=%d),fTPCNClsFindable=%d,fTRDPattern=%d,fTOFChi2=%d]\n", dfName, bITSClusterMap, bITSClusterSizes, bTPClsFindable, bTRDPattern, bTOFChi2);
      exitCode = BRANCH_DETECTION;
      break;
    }

    int fIndexCollisions = 0;
    track_iu->SetBranchAddress("fIndexCollisions", &fIndexCollisions);

    // loop over all tracks and reindex them
    auto entries = trackExtraTree->GetEntriesFast();
    int counter = 0;
    for (int i = 0; i < entries; i++) {
      trackExtraTree->GetEntry(i);
      track_iu->GetEntry(i);

      // Flag collisions
      hasCollision[i] = (fIndexCollisions >= 0);

      // Remove TPC only tracks, if they are not assoc. to a V0/Cascade/3Body
      if (tpcNClsFindable > 0 && TRDPattern == 0 && TOFChi2 < -1. &&
          (!bITSClusterMap || ITSClusterMap == 0) &&
          (!bITSClusterSizes || ITSClusterSizes == 0) &&
          !keepV0s[i] && !keepCascades[i] && !keep3Bodys[i] &&
          (!hasTrackQA || !keepTrackQA[i]) &&
          (!hasTrackedV0 || !keepTrackedV0s[i]) &&
          (!hasTrackedCasc || !keepTrackedCascs[i]) &&
          (!hasTracked3Body || !keepTracked3Bodys[i])) {
        counter++;
      } else {
        acceptedTracks[i] = i - counter;
      }
    }

    for (auto key2 : *treeList) {
      TString treeName = ((TObjString*)key2)->GetString().Data();

      // Connect trees but do not copy entries (using the clone function)
      // NOTE Basket size etc. are copied in CloneTree()
      if (outputDir == nullptr) {
        outputDir = outputFile->mkdir(dfName);
        printf("Writing to output folder %s\n", dfName);
      }
      outputDir->cd();

      auto inputTree = (TTree*)inputFile->Get(Form("%s/%s", dfName, treeName.Data()));
      printf("    Processing tree %s with %lld entries with total size %lld\n", treeName.Data(), inputTree->GetEntriesFast(), inputTree->GetTotBytes());
      auto outputTree = inputTree->CloneTree(0);
      outputTree->SetAutoFlush(0);

      std::vector<int*> indexList;
      std::vector<char*> vlaPointers;
      std::vector<int*> indexPointers;
      TObjArray* branches = inputTree->GetListOfBranches();
      for (int i = 0; i < branches->GetEntriesFast(); ++i) {
        TBranch* br = (TBranch*)branches->UncheckedAt(i);
        TString branchName(br->GetName());
        TString tableName(getTableName(branchName, treeName.Data()));
        // register index of track index ONLY
        if (!tableName.EqualTo("O2track")) {
          continue;
        }
        // detect VLA
        if (((TLeaf*)br->GetListOfLeaves()->First())->GetLeafCount() != nullptr) {
          printf("  *** FATAL ***: VLA detection is not supported\n");
          exitCode = NO_VLA_SUPPORT;
        } else if (branchName.BeginsWith("fIndexSlice")) {
          int* buffer = new int[2];
          memset(buffer, 0, 2 * sizeof(buffer[0]));
          vlaPointers.push_back(reinterpret_cast<char*>(buffer));
          inputTree->SetBranchAddress(br->GetName(), buffer);
          outputTree->SetBranchAddress(br->GetName(), buffer);

          indexList.push_back(buffer);
          indexList.push_back(buffer + 1);
        } else if (branchName.BeginsWith("fIndex") && !branchName.EndsWith("_size")) {
          int* buffer = new int;
          *buffer = 0;
          indexPointers.push_back(buffer);

          inputTree->SetBranchAddress(br->GetName(), buffer);
          outputTree->SetBranchAddress(br->GetName(), buffer);

          indexList.push_back(buffer);
        }
      }

      const bool processingTracked = treeName.BeginsWith("O2tracked");
      const bool processingTrackQA = treeName.BeginsWith("O2trackqa");
      const bool processingTracks = treeName.BeginsWith("O2track") && !processingTracked && !processingTrackQA; // matches any of the track tables and not tracked{v0s,cascase,3body} or trackqa;
      const bool processingCascades = treeName.BeginsWith("O2cascade");
      const bool processing3Body = treeName.BeginsWith("O2decay3body");
      const bool processingV0s = treeName.BeginsWith("O2v0");
      const bool processingAmbiguousTracks = treeName.BeginsWith("O2ambiguoustrack");

      const std::array<bool, 7> processingCheck{processingTracked, processingTrackQA, processingTracks, processingCascades, processing3Body, processingV0s, processingAmbiguousTracks};
      if (std::accumulate(processingCheck.begin(), processingCheck.end(), 0) > 1) {
        exitCode = PROCESSING_MORE_THAN_ONE;
        printf("  Processing more than branch (Mask=%d/%d/%d/%d/%d/%d/%d)\n", processingTracked, processingTrackQA, processingTracks, processingCascades, processing3Body, processingV0s, processingAmbiguousTracks);
        break;
      }

      int indexTrack0{-1}, indexTrack1{-1}, indexTrack2{-1};
      if (processingV0s) {
        inputTree->SetBranchAddress("fIndexTracks_Pos", &indexTrack0);
        outputTree->SetBranchAddress("fIndexTracks_Pos", &indexTrack0);
        inputTree->SetBranchAddress("fIndexTracks_Neg", &indexTrack1);
        outputTree->SetBranchAddress("fIndexTracks_Neg", &indexTrack1);
      } else if (processingCascades || processingTrackQA) {
        inputTree->SetBranchAddress("fIndexTracks", &indexTrack0);
        outputTree->SetBranchAddress("fIndexTracks", &indexTrack0);
      } else if (processing3Body) {
        inputTree->SetBranchAddress("fIndexTracks_0", &indexTrack0);
        outputTree->SetBranchAddress("fIndexTracks_0", &indexTrack0);
        inputTree->SetBranchAddress("fIndexTracks_1", &indexTrack1);
        outputTree->SetBranchAddress("fIndexTracks_1", &indexTrack1);
        inputTree->SetBranchAddress("fIndexTracks_2", &indexTrack2);
        outputTree->SetBranchAddress("fIndexTracks_2", &indexTrack2);
      } else if (processingTracked) {
        inputTree->SetBranchAddress("fIndexTracks", &indexTrack0);
        outputTree->SetBranchAddress("fIndexTracks", &indexTrack0);
        inputTree->SetBranchAddress("fIndexTracks_ITS", &indexTrack1);
        outputTree->SetBranchAddress("fIndexTracks_ITS", &indexTrack1);
      }

      // Here we loop over all tracks
      auto entries = inputTree->GetEntriesFast();
      for (int i = 0; i < entries; i++) {
        inputTree->GetEntry(i);
        bool fillThisEntry = true; // Assume by default tracks are accepted
        // Special case for Tracks, TracksExtra, TracksCov
        if (processingTracks) {
          if (acceptedTracks[i] < 0) {
            fillThisEntry = false;
          }
        } else if (processingTrackQA) {
          // Do nothing, just fill
        } else {
          // Other table than Tracks* --> reassign indices to Tracks
          for (const auto& idx : indexList) {
            int oldTrackIndex = *idx;

            // if negative, the index is unassigned.
            if (oldTrackIndex >= 0) {
              if (acceptedTracks[oldTrackIndex] < 0) {
                fillThisEntry = false;
              } else {
                *idx = acceptedTracks[oldTrackIndex];
              }
            }
          }
        }

        // Reassign prong track(s)
        if (processingV0s || processingTracked) { // Reassign two Prongs or reassign the add. ITS Prong + track
          indexTrack0 = acceptedTracks[indexTrack0];
          indexTrack1 = acceptedTracks[indexTrack1];
        } else if (processing3Body) { // Reassign three Prongs
          indexTrack0 = acceptedTracks[indexTrack0];
          indexTrack1 = acceptedTracks[indexTrack1];
          indexTrack2 = acceptedTracks[indexTrack2];
        } else if (processingCascades || processingTrackQA) { // Reassign the add. prong or the TPC ref. track
          indexTrack0 = acceptedTracks[indexTrack0];
        }

        // Keep only tracks which have no collision, see O2-3601
        if (processingAmbiguousTracks) {
          if (hasCollision[i]) {
            fillThisEntry = false;
          }
        }

        if (fillThisEntry) {
          outputTree->Fill();
        }
      }

      if (entries != outputTree->GetEntriesFast()) {
        printf("      Reduced from %lld to %lld entries\n", entries, outputTree->GetEntriesFast());
        // sanity check by hardcoding the trees for which we expect a reduction
        const TString tableName{removeVersionSuffix(outputTree->GetName())};
        static const std::array<TString, 4> checkNames{"O2track_iu", "O2trackextra", "O2trackcov_iu", "O2ambiguoustrack"};
        std::vector<bool> checks(checkNames.size(), false);
        for (size_t i{0}; i < checkNames.size(); ++i) {
          if (tableName.EqualTo(checkNames[i])) {
            checks[i] = true;
          }
        }
        if (std::none_of(checks.begin(), checks.end(), [](bool b) { return b; })) {
          exitCode = NOT_EXPECTED_REDUCTION;
          printf("       -> Reduction is not expected for this tree!\n");
          break;
        }
      }

      delete inputTree;

      for (auto& buffer : indexPointers) {
        delete buffer;
      }
      for (auto& buffer : vlaPointers) {
        delete[] buffer;
      }

      outputDir->cd();
      outputTree->Write();
      delete outputTree;
    }
    if (exitCode > 0) {
      break;
    }

    outputDir = nullptr;
  }
  inputFile->Close();

  outputFile->Write();
  outputFile->Close();

  // in case of failure, remove the incomplete file
  if (exitCode != 0) {
    printf("Removing incomplete output file %s.\n", outputFile->GetName());
    gSystem->Unlink(outputFile->GetName());
    return exitCode; // skip output below
  }

  clock.Stop();

  // Report savings
  auto sBefore = inputFile->GetSize();
  auto sAfter = outputFile->GetSize();
  if (sBefore <= 0 || sAfter <= 0) {
    printf("Warning: Empty input or output file after thinning!\n");
    exitCode = EMPTY_FILE;
  }
  auto spaceSaving = (1 - ((double)sAfter) / ((double)sBefore)) * 100;
  printf("Stats: After=%lld / Before=%lld Bytes ---> Saving %.1f%% diskspace!\n", sAfter, sBefore, spaceSaving);
  printf("Timing: CPU=%.2f (s);   Real=%.2f (s)\n", clock.CpuTime(), clock.RealTime());
  printf("End of AOD thinning.\n");

  return exitCode;
}
