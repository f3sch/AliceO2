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

/// \file DeviceCompileCheck.cuh
/// \brief Every per-element function of the tracking, in a kernel, for every
///        detector: the device half of "one algorithm, two backends".
///
/// Nothing here runs; the translation units
/// that instantiate it, one per detector (DeviceCompileCheck<Detector>.cu, so
/// that the detectors compile in parallel), only have to compile, with nvcc as
/// written and with hipcc after hipify-perl. A stage that stops compiling for
/// a device fails the `gtrk-device-compile-check` target, whichever backend is
/// being worked on.
/// The kernels are the generic launcher shape of the device backend: one
/// thread per element,
/// and a stage's emitter where it has one.

#ifndef GTRK_TEST_DEVICE_DEVICECOMPILECHECK_CUH_
#define GTRK_TEST_DEVICE_DEVICECOMPILECHECK_CUH_

#include "GenericTracking/Tracking/Material.h"
#include "GenericTracking/Tracking/NeighbourStage.h"
#include "GenericTracking/Tracking/Follow.h"
#include "GenericTracking/Tracking/RoadCount.h"
#include "GenericTracking/Tracking/RoadFit.h"
#include "GenericTracking/Tracking/Selection.h"
#include "GenericTracking/Tracking/TrackRefit.h"
#include "GenericTracking/Tracking/Tracklet.h"
#include "GenericTracking/Tracking/TrackletSort.h"
#include "GenericTracking/Tracking/TrackletStage.h"
#include "GenericTracking/Tracking/TripletStage.h"

#include "DetectorsBase/MatLayerCylSet.h"
#include "UniformMaterial.h"

#include <cstdint>

namespace o2::gtrk::device_check
{

/// Counts what a generator emits, which is all a kernel has to do with it.
struct CountingEmitter {
  std::uint32_t emitted = 0U;
  template <typename Record>
  GPUdi() void emit(const Record& /*record*/)
  {
    ++emitted;
  }
};

template <typename Body>
GPUg() void forEachElement(Body body, std::uint32_t elements)
{
  const std::uint32_t element = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (element < elements) {
    body(element);
  }
}

template <typename Generate>
GPUg() void generateElements(Generate generate, std::uint32_t elements, std::uint32_t* counts)
{
  const std::uint32_t element = (blockIdx.x * blockDim.x) + threadIdx.x;
  if (element < elements) {
    CountingEmitter emitter;
    generate(element, emitter);
    counts[element] = emitter.emitted;
  }
}

/// Taking a kernel's address instantiates its device code.
template <typename Kernel>
void use(Kernel* kernel)
{
  static_cast<void>(kernel);
}

template <DetectorLayout Det, typename Material>
void instantiateWithMaterial()
{
  constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  use(&generateElements<RoadFitStage<Det, Material, false, false, false, false>>);
  use(&forEachElement<TrackRefitElement<Det, Material>>);
  use(&generateElements<RoadFitStage<Det, Material, true, false, false, false>>);
  use(&generateElements<RoadFitStage<Det, Material, false, false, false, true>>);
  use(&forEachElement<FollowStage<Det, Material>>);
  use(&forEachElement<FollowStage<Det, Material, true>>);
  use(&forEachElement<FollowStage<Det, Material, true, false, false, false, true>>);
  use(&forEachElement<follow::KeepWitnesses<kMaxClusters>>);
  use(&forEachElement<follow::Refit<Det, Material>>);
  use(&forEachElement<follow::Bid<kMaxClusters>>);
  use(&forEachElement<follow::Resolve<kMaxClusters>>);
  use(&forEachElement<follow::Mark<kMaxClusters>>);
}

template <DetectorLayout Det>
void instantiate()
{
  constexpr std::size_t kMaxClusters = maxChainLayers<Det>;
  use(&generateElements<TrackletStage<Det>>);
  use(&generateElements<TrackletStageFor<Det, true>>);
  use(&generateElements<TrackletStageFor<Det, true, true>>);
  use(&generateElements<TrackletStageFor<Det, false, false, true>>);
  use(&generateElements<TrackletStageFor<Det, true, false, true>>);
  use(&generateElements<TrackletStageFor<Det, true, true, true>>);
  use(&generateElements<TripletStage<Det>>);
  use(&generateElements<TripletStage<Det, true>>);
  use(&generateElements<NeighbourStage<Det>>);
  use(&forEachElement<selection::BidReset<kMaxClusters>>);
  use(&forEachElement<selection::Bid<kMaxClusters>>);
  use(&forEachElement<selection::Decide<kMaxClusters>>);
  use(&forEachElement<selection::Prune<kMaxClusters>>);
  instantiateWithMaterial<Det, UniformMaterial>();
  instantiateWithMaterial<Det, o2::base::MatLayerCylSet>();
}

} // namespace o2::gtrk::device_check
#endif // GTRK_TEST_DEVICE_DEVICECOMPILECHECK_CUH_
