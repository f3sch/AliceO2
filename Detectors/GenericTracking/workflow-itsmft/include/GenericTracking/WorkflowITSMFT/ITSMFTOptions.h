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

#ifndef GTRK_WORKFLOW_ITSMFTOPTIONS_H_
#define GTRK_WORKFLOW_ITSMFTOPTIONS_H_

#include "GenericTracking/Detectors/ITS.h"
#include "GenericTracking/Detectors/MFT.h"
#include "DataFormatsITSMFT/DPLAlpideParamInitializer.h"
#include "Framework/ConfigContext.h"
#include "Framework/ConfigParamSpec.h"

#include <concepts>
#include <vector>

namespace o2::gtrk::workflow
{

/// The ITS's and MFT's reco workflow options beside every reco workflow's
/// (Workflow/RecoWorkflow.h, `addRecoOptions`): their staggered readout, the
/// adaptor's (`--disable-adaptor`, `--disable-root-output`), and the two
/// halves of a tracking inside o2-gpu-reco-workflow: `--tracking-in-gpu-reco`
/// before it in the pipe, `--adaptor-only` after it.
inline void addITSMFTRecoOptions(std::vector<o2::framework::ConfigParamSpec>& options)
{
  using o2::framework::ConfigParamSpec;
  using o2::framework::VariantType;
  o2::itsmft::DPLAlpideParamInitializer::addConfigOption(options);
  options.push_back(ConfigParamSpec{"disable-adaptor", VariantType::Bool, false, {"only the native tracks and vertices (GTRK/<DET>...), not O2's formats nor their files"}});
  options.push_back(ConfigParamSpec{"disable-root-output", VariantType::Bool, false, {"O2's formats without O2's writers"}});
  options.push_back(ConfigParamSpec{"tracking-in-gpu-reco",
                                    VariantType::Bool,
                                    false,
                                    {"the tracking runs inside o2-gpu-reco-workflow (--output-type gtrk-<det>-tracks), on the TPC's GPU: only the converter and "
                                     "the producer, before it in the pipe; the adaptor after it (--adaptor-only)"}});
  options.push_back(ConfigParamSpec{"adaptor-only", VariantType::Bool, false, {"only the adaptor (and O2's writers), after o2-gpu-reco-workflow in the pipe: its tracks, the measurements and clusters from upstream"}});
}

template <DetectorLayout Det>
[[nodiscard]] bool staggeredInput(const o2::framework::ConfigContext& context)
{
  if constexpr (std::same_as<Det, detectors::ITS>) {
    return o2::itsmft::DPLAlpideParamInitializer::isITSStaggeringEnabled(context);
  } else if constexpr (std::same_as<Det, detectors::MFT>) {
    return o2::itsmft::DPLAlpideParamInitializer::isMFTStaggeringEnabled(context);
  } else {
    return false;
  }
}

} // namespace o2::gtrk::workflow

#endif // GTRK_WORKFLOW_ITSMFTOPTIONS_H_
