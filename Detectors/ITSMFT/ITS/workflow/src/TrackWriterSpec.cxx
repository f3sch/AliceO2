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

/// @file   TrackWriterSpec.cxx

#include <vector>

#include "ITSWorkflow/TrackWriterSpec.h"
#include "DPLUtils/MakeRootTreeWriterSpec.h"
#include "DataFormatsITS/TrackITS.h"
#include "ITStracking/Definitions.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"

namespace o2::its
{

using namespace o2::framework;
template <typename T>
using BranchDefinition = MakeRootTreeWriterSpec::BranchDefinition<T>;
using LabelsType = std::vector<o2::MCCompLabel>;

DataProcessorSpec getTrackWriterSpec(bool useMC)
{
  auto logger = [](std::vector<o2::its::TrackITS> const& tracks) {
    LOG(info) << "ITSTrackWriter pulled " << tracks.size() << " tracks";
  };
  return MakeRootTreeWriterSpec("its-track-writer",
                                "o2trac_its.root",
                                MakeRootTreeWriterSpec::TreeAttributes{.name = "o2sim", .title = "Tree with ITS tracks"},
                                BranchDefinition<std::vector<o2::its::TrackITS>>{InputSpec{"tracks", "ITS", "TRACKS", 0}, "ITSTrack", logger},
                                BranchDefinition<std::vector<int>>{InputSpec{"trackClIdx", "ITS", "TRACKCLSID", 0}, "ITSTrackClusIdx"},
                                BranchDefinition<std::vector<Vertex>>{InputSpec{"vertices", "ITS", "VERTICES", 0}, "ITSVertices"},
                                BranchDefinition<LabelsType>{InputSpec{"labels", "ITS", "TRACKSMCTR", 0}, "ITSTrackMCTruth", (useMC ? 1 : 0), ""},
                                BranchDefinition<LabelsType>{InputSpec{"labelsVertices", "ITS", "VERTICESMCTR", 0}, "ITSVertexMCTruth", (useMC ? 1 : 0), ""},
                                BranchDefinition<std::vector<float>>{InputSpec{"purityVertices", "ITS", "VERTICESMCPUR", 0}, "ITSVertexMCPurity", (useMC ? 1 : 0), ""})();
}

} // namespace o2::its
