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

/// \file gtrk-study-workflow.cxx
/// \brief O2's clusters from files, with the study's sink.
///
/// `--detectors its,mft,trk`, a chain per detector, each detector's devices
/// and messages its own (gtrk-<det>-..., GTRK/<DET>...), each
/// device pipelined on its own (`--pipeline gtrk-its-measurement-producer:2`)
/// and configured on its own (`--gtrk-mft-measurement-sink "--qa-output
/// ..."`). With `--clusters-from-upstream` the clusters come from an O2
/// workflow before it in the pipe (`o2-its-reco-workflow --disable-tracking`),
/// without it from O2's cluster readers.
///
/// `its` or `mft` runs on O2's clusters: O2's own cluster reader, the
/// converter into our clusters, the producer and the sink.
/// Run it
/// in the timeframe's directory of an O2DPG production. With `--mc-vertices`
/// the MC vertices go to its sink too, for the iterations that
/// GTRKReco<DET>.iterationVertexSource seeds from them; with
/// GTRKReco<DET>.vertexer the sink finds its own, and takes
/// one or the other.
///
/// `trk` runs on ALICE 3's TRK and FT3 clusters: O2's TRK and FT3 cluster readers, the converter of
/// their per-layer streams,
/// the producer and the sink of the TRKFT3 layout, in the timeframe's
/// directory of MC/ALICE3.sh, with the TRK's and FT3's topology dictionaries
/// as CCDB conditions.
///
/// With `--compare-tracks` each chain also runs the tracking device
/// (gtrk-<det>-tracker, workflow/TrackerSpec.h) on the same measurements and
/// seeds, and its sink compares the device's tracks and vertices with its
/// own tracker's, every timeframe.

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
  workflowOptions.push_back(ConfigParamSpec{"detectors", VariantType::String, "", {"comma-separated: its, mft, trk (ALICE 3's barrel), a chain each"}});
  workflowOptions.push_back(
    ConfigParamSpec{"clusters-from-upstream", VariantType::Bool, false, {"the clusters from the workflow before this one in the pipe, not O2's readers"}});
  workflowOptions.push_back(
    ConfigParamSpec{"mc-vertices", VariantType::Bool, false, {"the MC vertices as the vertex list (any detector with an MC collision context)"}});
  workflowOptions.push_back(ConfigParamSpec{
    "compare-tracks", VariantType::Bool, false, {"also run the tracking device; each sink compares its tracks and vertices with its own tracker's"}});
  workflowOptions.push_back(ConfigParamSpec{
    "tracker-digests", VariantType::Bool, false, {"also run the tracking device; each sink logs its tracks' and vertices' digests, as the bench prints a device's"}});
  // Use the HBFUtils instance updated by --configKeyValues. This also makes
  // the timer advance TimingInfo::firstTForbit by nHBFPerTF for every entry.
  o2::raw::HBFUtilsInitializer::addConfigOption(workflowOptions, o2::raw::HBFUtilsInitializer::HBFUSrc);
}

#include "GenericTracking/Workflow/ClusterConverterSpec.h"
#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "GenericTracking/Detectors/TRKFT3.h"
#include "GenericTracking/Study/MCVertexReaderSpec.h"
#include "GenericTracking/Workflow/MeasurementProducerSpec.h"
#include "GenericTracking/Workflow/TrackerSpec.h"
#include "GenericTracking/Study/MeasurementSinkSpec.h"

#include "Framework/WorkflowSpec.h"
#include "ITSMFTWorkflow/ClusterReaderSpec.h"
#include "TRKWorkflow/ClusterReaderSpec.h"

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
  const bool mcVertices = context.options().get<bool>("mc-vertices");
  const bool upstream = context.options().get<bool>("clusters-from-upstream");
  const bool compareTracks = context.options().get<bool>("compare-tracks");
  const bool trackerDigests = context.options().get<bool>("tracker-digests");
  if (compareTracks && trackerDigests) {
    throw std::invalid_argument("--compare-tracks and --tracker-digests: one or the other");
  }
  namespace gs = o2::gtrk::study;
  const auto products = compareTracks ? gs::TrackerProducts::kCompare : (trackerDigests ? gs::TrackerProducts::kDigests : gs::TrackerProducts::kNone);
  WorkflowSpec specs;
  const auto chain = [&]<o2::gtrk::DetectorLayout Det>(const WorkflowSpec& readers) {
    if (!upstream) {
      specs.insert(specs.end(), readers.begin(), readers.end());
    }
    specs.push_back(gw::getClusterConverterSpec<Det>(useMC));
    specs.push_back(gw::getMeasurementProducerSpec<Det>());
    specs.push_back(o2::gtrk::study::getMeasurementSinkSpec<Det>(useMC, mcVertices, products));
    if (products != gs::TrackerProducts::kNone) {
      specs.push_back(gw::getTrackerSpec<Det>(mcVertices));
    }
    if (mcVertices) {
      specs.push_back(o2::gtrk::study::getMCVertexReaderSpec<Det>());
    }
  };
  std::istringstream detectors{context.options().get<std::string>("detectors")};
  for (std::string detector; std::getline(detectors, detector, ',');) {
    if (detector == "its") {
      chain.operator()<det::ITS>({o2::itsmft::getITSClusterReaderSpec(useMC, false, true, false)});
    } else if (detector == "mft") {
      chain.operator()<det::MFT>({o2::itsmft::getMFTClusterReaderSpec(useMC, false, true, false)});
    } else if (detector == "trk") {
      chain.operator()<det::TRKFT3>({o2::trk::getTRKClusterReaderSpec(useMC), o2::trk::getFT3ClusterReaderSpec(useMC)});
    } else {
      throw std::invalid_argument("--detectors takes its, mft and trk, not '" + detector + "'");
    }
  }

  const o2::raw::HBFUtilsInitializer hbfIni(context, specs);

  return specs;
}
