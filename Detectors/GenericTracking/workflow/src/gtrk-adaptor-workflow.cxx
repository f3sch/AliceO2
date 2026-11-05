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

/// \file gtrk-adaptor-workflow.cxx
/// \brief The O2-format adaptors, piped after o2-gtrk-reco-workflow: per
///        detector the tracking device's products as O2's tracks. Nothing
///        else is spawned; O2's writers are piped after it when files are wanted.

#include "CommonUtils/ConfigurableParam.h"
#include "Framework/ConfigParamSpec.h"

using namespace o2::framework;

// DPL looks it up by symbol before runDataProcessing.h is included, so it has
// to keep external linkage.
// NOLINTNEXTLINE(misc-use-internal-linkage)
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  workflowOptions.push_back(ConfigParamSpec{"disable-mc", VariantType::Bool, false, {"no MC labels"}});
  workflowOptions.push_back(ConfigParamSpec{"configKeyValues", VariantType::String, "", {"semicolon separated key=value strings"}});
  workflowOptions.push_back(ConfigParamSpec{"detectors", VariantType::String, "", {"comma-separated: its, mft, an adaptor each"}});
}

#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Workflow/AdaptorSpec.h"

#include "Framework/WorkflowSpec.h"

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
  WorkflowSpec specs;
  std::istringstream detectors{context.options().get<std::string>("detectors")};
  for (std::string detector; std::getline(detectors, detector, ',');) {
    if (detector == "its") {
      specs.push_back(gw::getAdaptorSpec<det::ITS>(useMC));
    } else if (detector == "mft") {
      specs.push_back(gw::getAdaptorSpec<det::MFT>(useMC));
    } else {
      throw std::invalid_argument("--detectors takes its and mft, not '" + detector + "'");
    }
  }
  return specs;
}
