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

/// \file DecayFitBatch.h
/// \brief Many decays' inputs to the decay vertex fit packed into columns, as
///        a device takes them, and their fit one decay at a time.
///
/// A batch holds every decay's daughters, tracked mothers and their rows one
/// after the other in flat arrays; a decay refers to its share by first index
/// and count. `fitBatchDecay` rebuilds one decay's spans over those arrays
/// and fits it: the same code on the host (`fitDecayBatch`) and in a device
/// kernel, so the two differ only in the arithmetic of their backends.

#ifndef GTRK_TRACKING_DECAYFITBATCH_H_
#define GTRK_TRACKING_DECAYFITBATCH_H_

#include "GenericTracking/Tracking/DecayVertexFit.h"

#include "GPUCommonDef.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace o2::gtrk
{

/// One decay of a batch: its share of the batch's daughters and tracked
/// mothers, its constraints, its tree and its vertices' start (`starts` of
/// them; none: kNoStart).
struct DecayFitInput {
  std::array<Vec3D, decay_fit_detail::kMaxConstrainedVertices> start{};
  std::uint32_t starts = 0U;
  std::uint32_t firstDaughter = 0U;
  std::uint32_t daughters = 0U;
  std::uint32_t firstMother = 0U;
  std::uint32_t mothers = 0U;
  DecayConstraints constraints;
  DecayTopology topology;
};

/// A tracked mother of a batch: `DecayMother` with its rows as a share of
/// the batch's row arrays, each array's own share kept, so a mother whose
/// arrays differ in length reaches the fit's input check as it was given.
struct DecayMotherRows {
  int vertex = 0;
  Propagator propagator;
  std::uint32_t firstRow = 0U;
  std::uint32_t rows = 0U; ///< planes
  std::uint32_t firstPosition = 0U;
  std::uint32_t positions = 0U;
  std::uint32_t firstCovariance = 0U;
  std::uint32_t covariances = 0U;
};

/// A batch's arrays where the fit reads them: host or device memory.
struct DecayFitBatchView {
  const DecayFitInput* decays = nullptr;
  const DecayDaughter* daughters = nullptr;
  const DecayMotherRows* mothers = nullptr;
  const TrackingSurface* planes = nullptr;
  const Vec2D* positions = nullptr;
  const Vec3D* covariance = nullptr;
  std::uint32_t count = 0U; ///< decays
};

static_assert(std::is_trivially_copyable_v<DecayFitInput> && std::is_trivially_copyable_v<DecayDaughter> && std::is_trivially_copyable_v<DecayMotherRows> && std::is_trivially_copyable_v<TrackingSurface>, "a batch's columns upload as bytes");
static_assert(std::is_trivially_copyable_v<DecayVertexFitResult>, "results download as bytes");

namespace decay_fit_detail
{

/// \a input's tracked mothers into \a mothers, their spans over \a batch's
/// rows.
GPUhdi() void batchMothers(const DecayFitBatchView& batch, const DecayFitInput& input, std::array<DecayMother, kMaxConstrainedVertices>& mothers)
{
  for (std::uint32_t a = 0; a < input.mothers; ++a) {
    const auto& rows = batch.mothers[input.firstMother + a];
    mothers[a] = {.vertex = rows.vertex, .planes = {batch.planes + rows.firstRow, rows.rows}, .positions = {batch.positions + rows.firstPosition, rows.positions}, .covariance = {batch.covariance + rows.firstCovariance, rows.covariances}, .propagator = rows.propagator};
  }
}

} // namespace decay_fit_detail

/// A batch on the host, filled decay by decay.
class DecayFitBatch
{
 public:
  /// Takes one decay's inputs (copied), its vertices' \a start among them,
  /// and returns its index in the batch. More tracked mothers or starts than
  /// the fit has vertices is refused here, as the fit could not hold them
  /// (std::invalid_argument).
  std::uint32_t add(std::span<const DecayDaughter> daughters, const DecayConstraints& constraints, const DecayTopology& topology, std::span<const DecayMother> mothers, std::span<const Vec3D> start)
  {
    if (mothers.size() > static_cast<std::size_t>(decay_fit_detail::kMaxConstrainedVertices) || start.size() > static_cast<std::size_t>(decay_fit_detail::kMaxConstrainedVertices)) {
      throw std::invalid_argument("decay fit batch: more tracked mothers or starts than vertices");
    }
    DecayFitInput input{.starts = static_cast<std::uint32_t>(start.size()), .firstDaughter = static_cast<std::uint32_t>(mDaughters.size()), .daughters = static_cast<std::uint32_t>(daughters.size()), .firstMother = static_cast<std::uint32_t>(mMothers.size()), .mothers = static_cast<std::uint32_t>(mothers.size()), .constraints = constraints, .topology = topology};
    std::ranges::copy(start, input.start.begin());
    mDaughters.insert(mDaughters.end(), daughters.begin(), daughters.end());
    for (const auto& mother : mothers) {
      mMothers.push_back({.vertex = mother.vertex, .propagator = mother.propagator, .firstRow = static_cast<std::uint32_t>(mPlanes.size()), .rows = static_cast<std::uint32_t>(mother.planes.size()), .firstPosition = static_cast<std::uint32_t>(mPositions.size()), .positions = static_cast<std::uint32_t>(mother.positions.size()), .firstCovariance = static_cast<std::uint32_t>(mCovariance.size()), .covariances = static_cast<std::uint32_t>(mother.covariance.size())});
      mPlanes.insert(mPlanes.end(), mother.planes.begin(), mother.planes.end());
      mPositions.insert(mPositions.end(), mother.positions.begin(), mother.positions.end());
      mCovariance.insert(mCovariance.end(), mother.covariance.begin(), mother.covariance.end());
    }
    mDecays.push_back(input);
    return static_cast<std::uint32_t>(mDecays.size() - 1U);
  }

  /// Takes \a other's decays after this batch's, in their order.
  void append(const DecayFitBatch& other)
  {
    const auto view = other.view();
    for (const auto& input : other.mDecays) {
      std::array<DecayMother, decay_fit_detail::kMaxConstrainedVertices> mothers{};
      decay_fit_detail::batchMothers(view, input, mothers);
      add({other.mDaughters.data() + input.firstDaughter, input.daughters}, input.constraints, input.topology, {mothers.data(), input.mothers}, {input.start.data(), input.starts});
    }
  }

  [[nodiscard]] std::size_t size() const noexcept { return mDecays.size(); }
  [[nodiscard]] bool empty() const noexcept { return mDecays.empty(); }
  void clear() noexcept
  {
    mDecays.clear();
    mDaughters.clear();
    mMothers.clear();
    mPlanes.clear();
    mPositions.clear();
    mCovariance.clear();
  }

  [[nodiscard]] DecayFitBatchView view() const noexcept
  {
    return {.decays = mDecays.data(), .daughters = mDaughters.data(), .mothers = mMothers.data(), .planes = mPlanes.data(), .positions = mPositions.data(), .covariance = mCovariance.data(), .count = static_cast<std::uint32_t>(mDecays.size())};
  }

  /// The arrays, for an upload.
  [[nodiscard]] std::span<const DecayFitInput> decays() const noexcept { return mDecays; }
  [[nodiscard]] std::span<const DecayDaughter> daughters() const noexcept { return mDaughters; }
  [[nodiscard]] std::span<const DecayMotherRows> mothers() const noexcept { return mMothers; }
  [[nodiscard]] std::span<const TrackingSurface> planes() const noexcept { return mPlanes; }
  [[nodiscard]] std::span<const Vec2D> positions() const noexcept { return mPositions; }
  [[nodiscard]] std::span<const Vec3D> covariance() const noexcept { return mCovariance; }

 private:
  std::vector<DecayFitInput> mDecays;
  std::vector<DecayDaughter> mDaughters;
  std::vector<DecayMotherRows> mMothers;
  std::vector<TrackingSurface> mPlanes;
  std::vector<Vec2D> mPositions;
  std::vector<Vec3D> mCovariance;
};

namespace decay_fit_detail
{

/// Decay \a index of \a batch's fitter, its tracked mothers' spans over the
/// batch's rows in \a mothers.
template <DetectorLayout Det, DecayFitVariant Variant, typename Material>
[[nodiscard]] GPUd() DecayVertexFitter<Det, Variant, Material> batchFitter(const DecayFitBatchView& batch, std::uint32_t index, const Material& material, std::array<DecayMother, kMaxConstrainedVertices>& mothers)
{
  const auto& input = batch.decays[index];
  batchMothers(batch, input, mothers);
  return DecayVertexFitter<Det, Variant, Material>({batch.daughters + input.firstDaughter, input.daughters}, material, input.constraints, input.topology, {mothers.data(), input.mothers});
}

} // namespace decay_fit_detail

/// The workspace [bytes] decay \a index of \a batch takes (0 for inputs the
/// fit refuses).
template <DetectorLayout Det, DecayFitVariant Variant, typename Material>
[[nodiscard]] GPUd() std::size_t batchDecayWorkspace(const DecayFitBatchView& batch, std::uint32_t index, const Material& material)
{
  std::array<DecayMother, decay_fit_detail::kMaxConstrainedVertices> mothers{};
  auto fitter = decay_fit_detail::batchFitter<Det, Variant>(batch, index, material, mothers);
  return fitter.workspace();
}

/// Decay \a index of \a batch fitted into \a result, in \a workspace.
template <DetectorLayout Det, DecayFitVariant Variant, typename Material>
GPUd() void fitBatchDecay(const DecayFitBatchView& batch, std::uint32_t index, const Material& material, std::span<std::byte> workspace, DecayVertexFitResult& result)
{
  std::array<DecayMother, decay_fit_detail::kMaxConstrainedVertices> mothers{};
  auto fitter = decay_fit_detail::batchFitter<Det, Variant>(batch, index, material, mothers);
  const auto& input = batch.decays[index];
  fitter.fit({input.start.data(), input.starts}, workspace, result);
}

/// Every decay of \a batch fitted on the host into \a results (one per
/// decay), one workspace reused, grown to the largest decay's.
template <DetectorLayout Det, DecayFitVariant Variant, typename Material>
void fitDecayBatch(const DecayFitBatchView& batch, const Material& material, std::span<DecayVertexFitResult> results)
{
  std::vector<std::max_align_t> storage;
  for (std::uint32_t i = 0; i < batch.count; ++i) {
    std::array<DecayMother, decay_fit_detail::kMaxConstrainedVertices> mothers{};
    auto fitter = decay_fit_detail::batchFitter<Det, Variant>(batch, i, material, mothers);
    const std::size_t bytes = fitter.workspace();
    const std::size_t words = (bytes + sizeof(std::max_align_t) - 1U) / sizeof(std::max_align_t);
    if (storage.size() < words) {
      storage.resize(words);
    }
    const auto& input = batch.decays[i];
    fitter.fit({input.start.data(), input.starts}, std::as_writable_bytes(std::span{storage}).first(bytes), results[i]);
  }
}

} // namespace o2::gtrk

#endif // GTRK_TRACKING_DECAYFITBATCH_H_
