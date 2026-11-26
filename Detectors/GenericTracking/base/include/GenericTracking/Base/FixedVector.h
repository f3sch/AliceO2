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

/// \file FixedVector.h
/// \brief At most a compile-time number of elements in place, a runtime count
///        of them used: storage bounded by a capacity, without allocation, on
///        the host and on the device.

#ifndef GTRK_BASE_FIXEDVECTOR_H_
#define GTRK_BASE_FIXEDVECTOR_H_

#include "GPUCommonDef.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>

namespace o2::gtrk
{

/// Elements a resize adds are value-initialised, as std::vector's. What
/// would go beyond the capacity is not written: the vector keeps what fits
/// and records the overflow (`overflowed`), which a consumer of the vector
/// as an input refuses.
template <typename T, std::size_t Capacity>
struct FixedVector {
  std::array<T, Capacity> items{};
  std::uint32_t count = 0;
  bool overflow = false;

  FixedVector() = default;
  /// The elements of \a values (host code's brace initialisation).
  FixedVector(std::initializer_list<T> values)
  {
    for (const auto& value : values) {
      push_back(value);
    }
  }

  [[nodiscard]] GPUhdi() static constexpr std::uint32_t capacity() { return Capacity; }
  [[nodiscard]] GPUhdi() std::uint32_t size() const { return count; }
  [[nodiscard]] GPUhdi() bool empty() const { return count == 0; }
  [[nodiscard]] GPUhdi() T* data() { return items.data(); }
  [[nodiscard]] GPUhdi() const T* data() const { return items.data(); }
  [[nodiscard]] GPUhdi() T* begin() { return items.data(); }
  [[nodiscard]] GPUhdi() T* end() { return items.data() + count; }
  [[nodiscard]] GPUhdi() const T* begin() const { return items.data(); }
  [[nodiscard]] GPUhdi() const T* end() const { return items.data() + count; }
  [[nodiscard]] GPUhdi() T& operator[](std::uint32_t i) { return items[i]; }
  [[nodiscard]] GPUhdi() const T& operator[](std::uint32_t i) const { return items[i]; }
  [[nodiscard]] GPUhdi() T& front() { return items[0]; }
  [[nodiscard]] GPUhdi() const T& front() const { return items[0]; }
  [[nodiscard]] GPUhdi() T& back() { return items[count - 1]; }
  [[nodiscard]] GPUhdi() const T& back() const { return items[count - 1]; }
  /// Whether a write went beyond the capacity.
  [[nodiscard]] GPUhdi() bool overflowed() const { return overflow; }
  GPUhdi() void clear()
  {
    count = 0;
    overflow = false;
  }
  GPUhdi() void push_back(const T& value) // NOLINT(readability-identifier-naming): the std containers' name
  {
    if (count < Capacity) {
      items[count++] = value;
    } else {
      overflow = true;
    }
  }
  GPUhdi() void pop_back() { --count; } // NOLINT(readability-identifier-naming): the std containers' name
  GPUhdi() void resize(std::uint32_t size)
  {
    if (size > Capacity) {
      overflow = true;
      size = Capacity;
    }
    for (std::uint32_t i = count; i < size; ++i) {
      items[i] = T{};
    }
    count = size;
  }
};

} // namespace o2::gtrk

#endif // GTRK_BASE_FIXEDVECTOR_H_
