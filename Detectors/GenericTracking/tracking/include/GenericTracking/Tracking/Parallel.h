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

/// \file Parallel.h
/// \brief Host parallelism without TBB in any header.
///
/// A parallel loop and a parallel sort, both in a task arena, implemented with
/// TBB in Parallel.cxx. VertexStage.h includes this header and is compiled
/// into the CUDA and HIP units, so TBB here would reach nvcc and hipcc.

#ifndef GTRK_TRACKING_PARALLEL_H_
#define GTRK_TRACKING_PARALLEL_H_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <type_traits>

namespace o2::gtrk
{

/// The consecutive source elements one host task of a stage takes: small
/// enough that the threads share a stage's work evenly, large enough that a
/// task's overhead stays small beside its work.
inline constexpr std::uint32_t kHostChunkSize = 256U;

} // namespace o2::gtrk

namespace o2::gtrk::detail
{

using ParallelForBody = void (*)(void*, size_t, size_t);
using ParallelSortLess = bool (*)(void*, std::uint32_t, std::uint32_t);

/// A TBB task arena kept behind the tracking library boundary.
class TaskArena
{
 public:
  explicit TaskArena(std::int32_t maxConcurrency = 0);
  ~TaskArena();
  TaskArena(const TaskArena&) = delete;
  TaskArena& operator=(const TaskArena&) = delete;
  TaskArena(TaskArena&&) = delete;
  TaskArena& operator=(TaskArena&&) = delete;

  [[nodiscard]] std::int32_t maxConcurrency() const noexcept;

 private:
  struct Impl;
  std::unique_ptr<Impl> mImpl;

  friend void parallelFor(TaskArena& arena, size_t begin, size_t end, size_t grainSize, void* context, ParallelForBody body);
  friend void parallelSortIndices(TaskArena& arena, size_t count, void* context, ParallelSortLess less, uint32_t* first);
  friend void parallelSortKeys(TaskArena& arena, std::span<std::uint64_t> keys);
};

/// The type-erased body of `parallelSort`.
void parallelSortIndices(TaskArena& arena, size_t count, void* context, ParallelSortLess less, uint32_t* first);

/// Sort \a keys ascending on \a arena's threads: the comparison of plain
/// integers, where `parallelSort` calls a comparator through a pointer.
void parallelSortKeys(TaskArena& arena, std::span<std::uint64_t> keys);

void parallelFor(TaskArena& arena, size_t begin, size_t end, size_t grainSize, void* context, ParallelForBody body);
[[nodiscard]] std::int32_t currentArenaThreadIndex() noexcept;

/// Sort \a indices in place by \a less on \a arena's threads (a sort outside
/// it would take every core TBB has, not the tracker's threads). The sort is
/// unstable, so \a less has to be a strict total order -- no two distinct
/// indices equivalent -- for the result to be a serial sort's; every caller
/// ends its comparison in an index tie-break.
template <typename Less>
void parallelSort(TaskArena& arena, std::span<std::uint32_t> indices, Less less)
{
  parallelSortIndices(arena, indices.size(), static_cast<void*>(&less), [](void* context, std::uint32_t lhs, std::uint32_t rhs) { return (*static_cast<Less*>(context))(lhs, rhs); }, indices.data());
}

/// Run the loop inside an arena, which limits its concurrency independently
/// of the process-wide TBB scheduler setting; \a grainSize > 0. The bounds
/// \a function gets are of \a end's type.
template <typename Index, typename Function>
void parallelFor(TaskArena& arena, std::type_identity_t<Index> begin, Index end, std::type_identity_t<Index> grainSize, Function function)
{
  parallelFor(arena, begin, end, grainSize, static_cast<void*>(&function), [](void* context, std::size_t first, std::size_t last) { (*static_cast<Function*>(context))(static_cast<Index>(first), static_cast<Index>(last)); });
}

} // namespace o2::gtrk::detail

#endif // GTRK_TRACKING_PARALLEL_H_
