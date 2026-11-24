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

/// \file DeviceStorage.h
/// \brief Typed, non-owning descriptors over DeviceArena storage.

#ifndef GTRK_TRACKING_GPU_DEVICESTORAGE_H_
#define GTRK_TRACKING_GPU_DEVICESTORAGE_H_

#include "GenericTracking/Base/Columns.h"
#include "GenericTracking/Tracking/GPU/DeviceArena.h"
#include "GenericTracking/Tracking/GPU/DeviceRuntime.h"

#include "GPUCommonDef.h"

#include <memlayout/wrapper.h>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace o2::gtrk
{

namespace device_storage_detail
{
template <typename T>
[[nodiscard]] DeviceAllocation allocateArray(DeviceArena& arena, std::size_t count, ArenaEnd end)
{
  static_assert(std::is_trivially_copyable_v<T>);
  if (count > std::numeric_limits<std::size_t>::max() / sizeof(T)) {
    return arena.allocate(std::numeric_limits<std::size_t>::max(), kArenaAlignment, end);
  }
  return arena.allocate(count * sizeof(T), kArenaAlignment, end);
}

/// \a copy(to, from) for each column of \a from and the same column of \a to,
/// both column views of one record type.
template <typename From, typename To, typename Copy>
void forEachColumnPair(const From& from, const To& to, const Copy& copy)
{
  from.apply([&](const auto&... source) { to.apply([&](auto... target) { (copy(target, source), ...); }); });
}

/// The arena of a descriptor, \a what naming it in the error; with \a tag,
/// also checked live (std::logic_error for storage of a rewound scope).
[[nodiscard]] inline DeviceArena& requireAllocated(DeviceArena* arena, const char* what, const ArenaTag* tag = nullptr)
{
  if (arena == nullptr) {
    throw std::logic_error(std::string{what} + " has not been allocated");
  }
  if (tag != nullptr) {
    arena->requireLive(*tag);
  }
  return *arena;
}
} // namespace device_storage_detail

template <typename T>
class DeviceSpan
{
 public:
  static_assert(std::is_trivially_copyable_v<T>);

  [[nodiscard]] bool allocate(DeviceArena& arena, std::size_t count, ArenaEnd end = ArenaEnd::kBase)
  {
    const DeviceAllocation allocation = device_storage_detail::allocateArray<T>(arena, count, end);
    if (!allocation) {
      return false;
    }
    mArena = &arena;
    mData = static_cast<T*>(allocation.data);
    mSize = count;
    mTag = allocation.tag;
    return true;
  }

  [[nodiscard]] T* data()
  {
    requireLive();
    return mData;
  }
  [[nodiscard]] const T* data() const
  {
    requireLive();
    return mData;
  }
  [[nodiscard]] std::size_t size() const
  {
    requireLive();
    return mSize;
  }
  [[nodiscard]] DeviceArena& arena() const { return device_storage_detail::requireAllocated(mArena, "device span"); }

 private:
  void requireLive() const { static_cast<void>(device_storage_detail::requireAllocated(mArena, "device span", &mTag)); }

  DeviceArena* mArena = nullptr;
  T* mData = nullptr;
  std::size_t mSize = 0;
  ArenaTag mTag;
};

template <typename T>
class DeviceObject
{
 public:
  static_assert(std::is_trivially_copyable_v<T>);

  [[nodiscard]] bool allocate(DeviceArena& arena, ArenaEnd end = ArenaEnd::kBase) { return mStorage.allocate(arena, 1U, end); }
  [[nodiscard]] T* data() { return mStorage.data(); }
  [[nodiscard]] const T* data() const { return mStorage.data(); }
  [[nodiscard]] DeviceArena& arena() const { return mStorage.arena(); }

 private:
  DeviceSpan<T> mStorage;
};

template <template <template <class> class> class S>
class DeviceColumns
{
 public:
  using View = memlayout::Wrapper<S, memlayout::pointer, memlayout::Layout::soa>;
  using ConstView = memlayout::Wrapper<S, memlayout::const_pointer, memlayout::Layout::soa>;
  using Pointers = S<memlayout::pointer>;

  /// Allocate all columns or none. Successful allocations inherit the parent
  /// scope's lifetime; a failure rewinds the temporary scope.
  [[nodiscard]] bool allocate(DeviceArena& arena, std::size_t rows, ArenaEnd end = ArenaEnd::kBase)
  {
    auto transaction = arena.scope(end);
    Pointers pointers{};
    bool complete = true;
    pointers.apply([&](auto&... pointer) {
      const auto allocateOne = [&](auto& out) {
        if (!complete) {
          return;
        }
        using Element = std::remove_pointer_t<std::remove_reference_t<decltype(out)>>;
        const DeviceAllocation allocation = device_storage_detail::allocateArray<Element>(arena, rows, end);
        if (!allocation) {
          complete = false;
          return;
        }
        out = static_cast<Element*>(allocation.data);
      };
      (allocateOne(pointer), ...);
    });
    if (!complete) {
      return false;
    }
    transaction.commit();
    mArena = &arena;
    mPointers = pointers;
    mRows = rows;
    mTag = arena.currentTag(end);
    return true;
  }

  [[nodiscard]] std::size_t size() const
  {
    requireLive();
    return mRows;
  }
  [[nodiscard]] View view()
  {
    requireLive();
    return View{mPointers};
  }
  [[nodiscard]] ConstView view() const
  {
    requireLive();
    return ConstView{mPointers};
  }
  [[nodiscard]] DeviceArena& arena() const { return device_storage_detail::requireAllocated(mArena, "device columns"); }

 private:
  void requireLive() const { static_cast<void>(device_storage_detail::requireAllocated(mArena, "device columns", &mTag)); }

  DeviceArena* mArena = nullptr;
  Pointers mPointers{};
  std::size_t mRows = 0;
  ArenaTag mTag;
};

template <typename T>
void upload(const T* source, DeviceSpan<T>& destination, std::size_t count, DeviceStream& stream)
{
  if (count > destination.size()) {
    throw std::out_of_range("device span upload exceeds its destination");
  }
  if (count == 0U) {
    return;
  }
  destination.arena().runtime().copyToDevice(destination.data(), source, count * sizeof(T), stream);
}

template <typename T>
void download(const DeviceSpan<T>& source, T* destination, std::size_t count, DeviceStream& stream)
{
  if (count > source.size()) {
    throw std::out_of_range("device span download exceeds its source");
  }
  if (count == 0U) {
    return;
  }
  source.arena().runtime().copyToHost(destination, source.data(), count * sizeof(T), stream);
}

/// The same into \a destination, sized to \a count.
template <typename T>
void download(const DeviceSpan<T>& source, std::vector<T>& destination, std::size_t count, DeviceStream& stream)
{
  destination.resize(count);
  download(source, destination.data(), count, stream);
}

template <typename T>
void uploadObject(const T& source, DeviceObject<T>& destination, DeviceStream& stream)
{
  destination.arena().runtime().copyToDevice(destination.data(), &source, sizeof(T), stream);
}

/// Upload \a rows rows of any host column view -- a DPL message's columns
/// included, which are not a BoundedColumns. The record type cannot be deduced
/// through the nested ConstView, so callers name it: `upload<MeasurementT>(...)`.
template <template <template <class> class> class S>
void upload(typename DeviceColumns<S>::ConstView source, std::size_t rows, DeviceColumns<S>& destination, DeviceStream& stream)
{
  if (rows > destination.size()) {
    throw std::out_of_range("device column upload exceeds its destination");
  }
  if (rows == 0U) {
    return;
  }
  auto& runtime = destination.arena().runtime();
  device_storage_detail::forEachColumnPair(source, destination.view(), [&](auto* to, const auto* from) { runtime.copyToDevice(to, from, rows * sizeof(*to), stream); });
}

/// Download \a rows rows of \a source into \a destination from row \a
/// destinationOffset, which must hold them.
template <template <template <class> class> class S>
void download(const DeviceColumns<S>& source, BoundedColumns<S>& destination, std::size_t destinationOffset, std::size_t rows, DeviceStream& stream)
{
  if (rows > source.size() || destinationOffset > destination.size() || rows > destination.size() - destinationOffset) {
    throw std::out_of_range("device column download exceeds its source or its destination");
  }
  if (rows == 0U) {
    return;
  }
  auto& runtime = source.arena().runtime();
  device_storage_detail::forEachColumnPair(source.view(), destination.view(), [&](auto* to, const auto* from) { runtime.copyToHost(to + destinationOffset, from, rows * sizeof(*to), stream); });
}

template <template <template <class> class> class S>
void download(const DeviceColumns<S>& source, BoundedColumns<S>& destination, std::size_t rows, DeviceStream& stream)
{
  if (rows > source.size()) {
    throw std::out_of_range("device column download exceeds its source");
  }
  destination.resize(rows);
  download(source, destination, 0U, rows, stream);
}

namespace device_storage_detail
{
/// The copies of \a rows rows of \a source, from row \a sourceOffset, into
/// \a destination at row \a destinationOffset, a copy per column, appended to
/// \a copies.
template <template <template <class> class> class S>
void columnCopies(const DeviceColumns<S>& source, std::size_t sourceOffset, DeviceColumns<S>& destination, std::size_t destinationOffset, std::size_t rows, std::vector<DeviceCopy>& copies)
{
  if (sourceOffset > source.size() || rows > source.size() - sourceOffset || destinationOffset > destination.size() || rows > destination.size() - destinationOffset) {
    throw std::out_of_range("device column copy exceeds its source or its destination");
  }
  if (rows == 0U) {
    return;
  }
  forEachColumnPair(source.view(), destination.view(), [&](auto* to, const auto* from) { copies.push_back({.to = to + destinationOffset, .from = from + sourceOffset, .bytes = rows * sizeof(*to)}); });
}
} // namespace device_storage_detail

/// Copy \a rows rows of \a source, from row \a sourceOffset, into
/// \a destination at row \a destinationOffset, on the device: all columns
/// by one kernel, queued on \a stream.
template <template <template <class> class> class S>
void copyColumns(const DeviceColumns<S>& source, std::size_t sourceOffset, DeviceColumns<S>& destination, std::size_t destinationOffset, std::size_t rows, DeviceStream& stream)
{
  std::vector<DeviceCopy> copies;
  device_storage_detail::columnCopies(source, sourceOffset, destination, destinationOffset, rows, copies);
  if (!copies.empty()) {
    destination.arena().runtime().copyEachOnDevice(copies, stream);
  }
}

/// Records of a column set kept as chunks, each exactly sized, in the order
/// they came, and gathered into one column set when they are all there.
/// The chunks live in the arena scope that was
/// current on their end when they were appended; gathering them after it
/// closed throws, as any use of stale device storage does.
template <template <template <class> class> class S>
class DeviceColumnChunks
{
 public:
  /// Copy \a rows rows of \a source, from row \a sourceOffset, into a new
  /// chunk on \a end of \a arena. False, and nothing kept, if the arena refuses.
  [[nodiscard]] bool append(DeviceArena& arena, ArenaEnd end, const DeviceColumns<S>& source, std::size_t sourceOffset, std::size_t rows, DeviceStream& stream)
  {
    if (rows == 0U) {
      return true;
    }
    DeviceColumns<S> chunk;
    if (!chunk.allocate(arena, rows, end)) {
      return false;
    }
    copyColumns(source, sourceOffset, chunk, 0U, rows, stream);
    mChunks.push_back(chunk);
    mRows += rows;
    return true;
  }

  /// \a columns, all their rows, as the next chunk, where they were allocated.
  void adopt(const DeviceColumns<S>& columns)
  {
    mChunks.push_back(columns);
    mRows += columns.size();
  }

  /// Copy every chunk, in the order appended, into \a destination from row
  /// \a destinationOffset: all chunks' columns by one kernel.
  void gather(DeviceColumns<S>& destination, std::size_t destinationOffset, DeviceStream& stream) const
  {
    std::vector<DeviceCopy> copies;
    std::size_t at = destinationOffset;
    for (const auto& chunk : mChunks) {
      device_storage_detail::columnCopies(chunk, 0U, destination, at, chunk.size(), copies);
      at += chunk.size();
    }
    if (!copies.empty()) {
      destination.arena().runtime().copyEachOnDevice(copies, stream);
    }
  }

  /// The inverse of `gather` for the first \a chunks chunks: each copied
  /// back from \a source, the chunks laid out in order from row
  /// \a sourceOffset, all chunks' columns by one kernel.
  void scatter(const DeviceColumns<S>& source, std::size_t sourceOffset, std::size_t chunks, DeviceStream& stream) const
  {
    std::vector<DeviceCopy> copies;
    std::size_t at = sourceOffset;
    for (std::size_t k = 0; k < chunks && k < mChunks.size(); ++k) {
      auto chunk = mChunks[k];
      device_storage_detail::columnCopies(source, at, chunk, 0U, chunk.size(), copies);
      at += chunk.size();
    }
    if (!copies.empty()) {
      source.arena().runtime().copyEachOnDevice(copies, stream);
    }
  }

  /// Download every chunk, in the order appended, into \a destination from
  /// row \a destinationOffset, which must hold them: no gather on the device.
  void download(BoundedColumns<S>& destination, std::size_t destinationOffset, DeviceStream& stream) const
  {
    std::size_t at = destinationOffset;
    for (const auto& chunk : mChunks) {
      gtrk::download(chunk, destination, at, chunk.size(), stream);
      at += chunk.size();
    }
  }

  /// Forget the chunks; their storage goes with its scope.
  void clear() noexcept
  {
    mChunks.clear();
    mRows = 0U;
  }

  [[nodiscard]] std::size_t rows() const noexcept { return mRows; }
  [[nodiscard]] std::span<const DeviceColumns<S>> chunks() const noexcept { return mChunks; }

 private:
  std::vector<DeviceColumns<S>> mChunks;
  std::size_t mRows = 0U;
};

} // namespace o2::gtrk

#endif // GTRK_TRACKING_GPU_DEVICESTORAGE_H_
