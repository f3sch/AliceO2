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

/// \file gtrk-trk-reco-workflow.cxx
/// \brief ALICE 3's tracking on the TRK's and FT3's clusters: GenericTracking's
///        tracks and vertices out (TrackerSpec.h).
///
/// O2's TRK and FT3 cluster readers (unless `--clusters-from-upstream`), the
/// converter of their per-layer streams, the measurement producer and the
/// tracker of the TRKFT3 layout, the topology dictionaries as CCDB
/// conditions. With `--vertices-from-upstream` the tracker takes the native
/// vertex seeds (GTRK/TRKVERTICES) from the workflow before it in the pipe.
/// With `--tof-layers` the layout is TRKFT3TOF (Detectors/TRKFT3TOF.h): the
/// IOTOF's clusters too (O2's TF3 cluster reader, unless
/// `--clusters-from-upstream`), found, owned and fitted by the tracking as
/// the TRK's, the material LUT `TRK/Param/MatLUTTOF`. The products are the
/// same: GTRK/TRKTRACKS, whose cluster references then hold TOF rows.
/// No adaptor: ALICE 3 has no O2 track format yet.

#include "GenericTracking/Workflow/RecoWorkflow.h"

#include "CommonUtils/ConfigurableParam.h"
#include "Framework/CallbacksPolicy.h"

using namespace o2::framework;

// DPL looks both of these up by symbol before runDataProcessing.h is included,
// so they have to keep external linkage.
// NOLINTNEXTLINE(misc-use-internal-linkage)
void customize(std::vector<o2::framework::CallbacksPolicy>& policies) { o2::raw::HBFUtilsInitializer::addNewTimeSliceCallback(policies); }

// NOLINTNEXTLINE(misc-use-internal-linkage)
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  o2::gtrk::workflow::addRecoOptions(workflowOptions);
  workflowOptions.push_back(ConfigParamSpec{"tof-layers", VariantType::Bool, false, {"track with the IOTOF's layers (the TRKFT3TOF layout)"}});
}

#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Detectors/TRKFT3TOF.h"

#include "Framework/runDataProcessing.h" // IWYU pragma: keep
#include "IOTOFWorkflow/ClusterReaderSpec.h"
#include "TRKWorkflow/ClusterReaderSpec.h"

#include <string>

WorkflowSpec defineDataProcessing(const ConfigContext& context)
{
  o2::conf::ConfigurableParam::updateFromString(context.options().get<std::string>("configKeyValues"));

  const bool useMC = !context.options().get<bool>("disable-mc");
  auto specs = context.options().get<bool>("tof-layers")
                 ? o2::gtrk::workflow::trackingChain<o2::gtrk::detectors::TRKFT3TOF>(context, {o2::trk::getTRKClusterReaderSpec(useMC), o2::trk::getFT3ClusterReaderSpec(useMC), o2::iotof::getIOTOFClusterReaderSpec(useMC)})
                 : o2::gtrk::workflow::trackingChain<o2::gtrk::detectors::TRKFT3>(context, {o2::trk::getTRKClusterReaderSpec(useMC), o2::trk::getFT3ClusterReaderSpec(useMC)});

  const o2::raw::HBFUtilsInitializer hbfIni(context, specs);

  return specs;
}
