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

/// \file gtrk-its3-reco-workflow.cxx
/// \brief ITS3's tracking: O2's ITS3 clusters in (the ITS's cluster reader,
///        unless `--clusters-from-upstream`), the converter into our
///        clusters, the measurement producer and the tracker, its tracks and
///        vertices GTRK/IT3...; no O2 formats yet.

#include "GenericTracking/Workflow/RecoWorkflow.h"
#include "GenericTracking/WorkflowITSMFT/ITSMFTOptions.h"

#include "CommonUtils/ConfigurableParam.h"
#include "DataFormatsITSMFT/DPLAlpideParamInitializer.h"
#include "Framework/CallbacksPolicy.h"

using namespace o2::framework;

// DPL looks these up by symbol before runDataProcessing.h is included, so
// they have to keep external linkage.
// NOLINTNEXTLINE(misc-use-internal-linkage)
void customize(std::vector<o2::framework::CallbacksPolicy>& policies) { o2::raw::HBFUtilsInitializer::addNewTimeSliceCallback(policies); }

// NOLINTNEXTLINE(misc-use-internal-linkage)
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  o2::gtrk::workflow::addRecoOptions(workflowOptions);
  o2::itsmft::DPLAlpideParamInitializer::addConfigOption(workflowOptions);
}

#include "GenericTracking/Detectors/ITS3.h"

#include "Framework/runDataProcessing.h" // IWYU pragma: keep
#include "ITSMFTWorkflow/ClusterReaderSpec.h"

#include <string>

WorkflowSpec defineDataProcessing(const ConfigContext& context)
{
  o2::conf::ConfigurableParam::updateFromString(context.options().get<std::string>("configKeyValues"));

  namespace gw = o2::gtrk::workflow;
  using Det = o2::gtrk::detectors::ITS3;
  const bool useMC = !context.options().get<bool>("disable-mc");
  auto specs = gw::trackingChain<Det>(context, {o2::itsmft::getITSClusterReaderSpec(useMC, gw::staggeredInput<Det>(context), true, false)});

  const o2::raw::HBFUtilsInitializer hbfIni(context, specs);

  return specs;
}
