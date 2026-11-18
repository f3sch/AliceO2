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

/// \file DeviceDecayFitter.h
/// \brief The decay vertex fit on a device: a batch of decays (DecayFitBatch.h),
///        one thread each, vendor-neutral.
///
/// A fit uploads the batch's columns, has every decay's workspace size
/// computed on the device (the fitter's own sizing), and fits the decays in
/// launches of as many as the arena's free memory holds, each decay's
/// workspace a slice of it at its offset in the launch. A decay whose
/// workspace alone exceeds the free memory is not fitted: its result says
/// `kTooLarge`. A launch's slices are reused only after it finished.

#ifndef GTRK_TRACKING_GPU_DEVICEDECAYFITTER_H_
#define GTRK_TRACKING_GPU_DEVICEDECAYFITTER_H_

#include "GenericTracking/Tracking/DecayFitBatch.h"
#include "GenericTracking/Tracking/Definitions.h"

#include <cstddef>
#include <span>

namespace o2::gtrk
{

template <DetectorLayout Det, typename Material>
class DeviceDecayFitter
{
 public:
  /// The wall times of a fit's parts [s], each ended by the stream's
  /// synchronisation, and its launches.
  struct Timing {
    double uploadSeconds = 0.;
    double sizeSeconds = 0.;
    double fitSeconds = 0.;
    double downloadSeconds = 0.;
    std::size_t launches = 0U;
    std::size_t refused = 0U;   ///< decays whose workspace did not fit
    std::size_t peakBytes = 0U; ///< the arena's peak
  };

  DeviceDecayFitter() = default;
  virtual ~DeviceDecayFitter() = default;
  DeviceDecayFitter(const DeviceDecayFitter&) = delete;
  DeviceDecayFitter& operator=(const DeviceDecayFitter&) = delete;
  DeviceDecayFitter(DeviceDecayFitter&&) = delete;
  DeviceDecayFitter& operator=(DeviceDecayFitter&&) = delete;

  /// \a material uploaded for the fits that follow; false if the arena cannot
  /// hold it.
  [[nodiscard]] virtual bool configure(const Material& material) = 0;
  /// \a batch's decays under \a variant into \a results (one per decay, host
  /// memory); false if the arena cannot hold the batch's inputs and results.
  [[nodiscard]] virtual bool fit(DecayFitVariant variant, const DecayFitBatch& batch, std::span<DecayVertexFitResult> results, Timing& timing) = 0;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GPU_DEVICEDECAYFITTER_H_
