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

/// \file TruthEvent.h
/// \brief What a truth study reads for one timeframe, the same from the bench
///        and from the sink.

#ifndef GTRK_STUDY_TRUTHEVENT_H_
#define GTRK_STUDY_TRUTHEVENT_H_

#include "GenericTracking/Study/Alice/Hits.h"
#include "GenericTracking/Study/ParticleChains.h"
#include "GenericTracking/Study/TrackComparison.h"
#include "GenericTracking/Tracking/Parallel.h"
#include "GenericTracking/Tracking/Propagator.h"
#include "GenericTracking/Tracking/TimeFrame.h"
#include "GenericTracking/Tracking/TrackingSurface.h"

#include "DetectorsBase/MatLayerCylSet.h"

#include <cstdint>
#include <span>

namespace o2::gtrk::study
{

/// One timeframe as a truth study sees it. It owns nothing: the host keeps
/// every part alive for the call it hands the event to, and sets every
/// pointer but `hits`.
template <DetectorLayout Det>
struct TruthEvent {
  std::uint32_t timeframe = 0U; ///< the host's count from 0
  const TimeFrame<Det>* frame = nullptr;
  const TrackingGeometry<Det>* geometry = nullptr;
  const o2::base::MatLayerCylSet* material = nullptr;
  Propagator propagator;          ///< the tracking's, with its field
  std::span<const TruthRow> rows; ///< one per row of the frame
  const ParticleTable* particles = nullptr;
  const HitIndex* hits = nullptr; ///< none without the MC hits
  const ParticleChains* chains = nullptr;
  o2::gtrk::detail::TaskArena* arena = nullptr; ///< the studies' threads
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_TRUTHEVENT_H_
