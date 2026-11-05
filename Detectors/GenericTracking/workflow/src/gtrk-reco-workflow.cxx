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

/// \file gtrk-reco-workflow.cxx
/// \brief The ITS's and the MFT's tracking: O2's clusters in, GenericTracking's
///        tracks and vertices out (TrackerSpec.h).
///
/// `--detectors its,mft`, a chain per detector, each detector's devices and
/// messages its own (gtrk-<det>-..., GTRK/<DET>...): O2's cluster reader
/// (unless `--clusters-from-upstream`), the converter into our clusters, the
/// measurement producer and the tracker. With `--vertices-from-upstream` the
/// tracker takes the native vertex seeds (GTRK/<DET>VERTICES) from the
/// workflow before it in the pipe.

#include "CommonUtils/ConfigurableParam.h"
#include "DetectorsRaw/HBFUtilsInitializer.h"
#include "Framework/CallbacksPolicy.h"
#include "Framework/ConfigParamSpec.h"

using namespace o2::framework;

// DPL looks both of these up by symbol before runDataProcessing.h is included,
// so they have to keep external linkage.
// NOLINTNEXTLINE(misc-use-internal-linkage)
void customize(std::vector<o2::framework::CallbacksPolicy>& policies) { o2::raw::HBFUtilsInitializer::addNewTimeSliceCallback(policies); }

// NOLINTNEXTLINE(misc-use-internal-linkage)
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  workflowOptions.push_back(ConfigParamSpec{"disable-mc", VariantType::Bool, false, {"do not propagate MC labels"}});
  workflowOptions.push_back(ConfigParamSpec{"configKeyValues", VariantType::String, "", {"semicolon separated key=value strings"}});
  workflowOptions.push_back(ConfigParamSpec{"detectors", VariantType::String, "", {"comma-separated: its, mft, a chain each"}});
  workflowOptions.push_back(
    ConfigParamSpec{"clusters-from-upstream", VariantType::Bool, false, {"the clusters from the workflow before this one in the pipe, not O2's readers"}});
  workflowOptions.push_back(
    ConfigParamSpec{"vertices-from-upstream", VariantType::Bool, false, {"the native vertex seeds (GTRK/<DET>VERTICES) from the workflow before this one in the pipe"}});
  o2::raw::HBFUtilsInitializer::addConfigOption(workflowOptions, o2::raw::HBFUtilsInitializer::HBFUSrc);
}

#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Workflow/ClusterConverterSpec.h"
#include "GenericTracking/Workflow/MeasurementProducerSpec.h"
#include "GenericTracking/Workflow/TrackerSpec.h"

#include "Framework/WorkflowSpec.h"
#include "ITSMFTWorkflow/ClusterReaderSpec.h"

#include <sstream>
#include <stdexcept>
#include <string>

#include "Framework/runDataProcessing.h" // IWYU pragma: keep

WorkflowSpec defineDataProcessing(const ConfigContext& context)
{
  o2::conf::ConfigurableParam::updateFromString(context.options().get<std::string>("configKeyValues"));

  namespace gw = o2::gtrk::workflow;
  namespace det = o2::gtrk::detectors;
  const bool useMC = !context.options().get<bool>("disable-mc");
  const bool upstream = context.options().get<bool>("clusters-from-upstream");
  const bool withVertices = context.options().get<bool>("vertices-from-upstream");
  WorkflowSpec specs;
  const auto chain = [&]<o2::gtrk::DetectorLayout Det>(const DataProcessorSpec& reader) {
    if (!upstream) {
      specs.push_back(reader);
    }
    specs.push_back(gw::getClusterConverterSpec<Det>(useMC));
    specs.push_back(gw::getMeasurementProducerSpec<Det>());
    specs.push_back(gw::getTrackerSpec<Det>(withVertices));
  };
  std::istringstream detectors{context.options().get<std::string>("detectors")};
  for (std::string detector; std::getline(detectors, detector, ',');) {
    if (detector == "its") {
      chain.operator()<det::ITS>(o2::itsmft::getITSClusterReaderSpec(useMC, false, true, false));
    } else if (detector == "mft") {
      chain.operator()<det::MFT>(o2::itsmft::getMFTClusterReaderSpec(useMC, false, true, false));
    } else {
      throw std::invalid_argument("--detectors takes its and mft, not '" + detector + "' (ALICE 3: o2-gtrk-alice3-reco-workflow)");
    }
  }

  const o2::raw::HBFUtilsInitializer hbfIni(context, specs);

  return specs;
}
