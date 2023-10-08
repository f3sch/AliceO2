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

/// \author felix.schlepper@cern.ch

#include "GlobalTrackingStudy/TPCV0Study.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "CommonUtils/ConfigurableParam.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/CallbacksPolicy.h"
#include "GlobalTrackingWorkflowHelpers/InputHelper.h"
#include "DetectorsRaw/HBFUtilsInitializer.h"

using namespace o2::framework;
using GID = o2::dataformats::GlobalTrackID;
using DetID = o2::detectors::DetID;

// ------------------------------------------------------------------
void customize(std::vector<o2::framework::CallbacksPolicy>& policies)
{
  o2::raw::HBFUtilsInitializer::addNewTimeSliceCallback(policies);
}

// we need to add workflow options before including Framework/runDataProcessing
void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  // option allowing to set parameters
  std::vector<o2::framework::ConfigParamSpec> options{
    {"disable-mc", o2::framework::VariantType::Bool, false, {"disable MC propagation"}},
    {"track-sources", VariantType::String, std::string{GID::ALL}, {"comma-separated list of track sources to use"}},
    {"disable-root-input", VariantType::Bool, false, {"disable root-files input reader"}},
    {"configKeyValues", VariantType::String, "", {"Semicolon separated key=value strings ..."}}};
  o2::raw::HBFUtilsInitializer::addConfigOption(options);
  std::swap(workflowOptions, options);
}

// ------------------------------------------------------------------

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(ConfigContext const& configcontext)
{
  WorkflowSpec specs;

  GID::mask_t allowedSourcesTrc = GID::getSourcesMask("ITS,ITS-TPC,TPC");

  // Update the (declared) parameters if changed from the command line
  o2::conf::ConfigurableParam::updateFromString(configcontext.options().get<std::string>("configKeyValues"));
  auto useMC = !configcontext.options().get<bool>("disable-mc");
  GID::mask_t src = allowedSourcesTrc & GID::getSourcesMask(configcontext.options().get<std::string>("track-sources"));
  src |= GID::getSourcesMask("TPC"); // guarantee at least TPConly tracks
  o2::globaltracking::InputHelper::addInputSpecs(configcontext, specs, src, src, src, useMC);
  o2::globaltracking::InputHelper::addInputSpecsPVertex(configcontext, specs, useMC);
  o2::globaltracking::InputHelper::addInputSpecsSVertex(configcontext, specs);
  specs.emplace_back(o2::trackstudy::getTPCV0StudySpec(src, useMC));

  // configure dpl timer to inject correct firstTForbit: start from the 1st orbit of TF containing 1st sampled orbit
  o2::raw::HBFUtilsInitializer hbfIni(configcontext, specs);
  // write the configuration used for the studies workflow
  o2::conf::ConfigurableParam::writeINI("o2_tpc_v0_study.ini");

  return specs;
}