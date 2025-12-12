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

#include <vector>
#include <array>
#include <algorithm>
#include <limits>
#include <cstdint>

namespace o2::its::utils
{

using Bracket = std::pair<int32_t, int32_t>;
constexpr Bracket InvalidBracket{std::numeric_limits<int32_t>::min(), std::numeric_limits<int32_t>::min()};

template <size_t N>
Bracket computeSmallestBracket(const std::array<Bracket, N>& intervals)
{
  // Collect valid intervals
  std::vector<Bracket> validIntervals;
  for (const auto& inter : intervals) {
    if (inter != InvalidBracket) {
      validIntervals.push_back(inter);
    }
  }

  if (validIntervals.empty()) {
    return InvalidBracket;
  }

  // Check for gaps: sort by start, verify each overlaps/touches the next
  std::sort(validIntervals.begin(), validIntervals.end());
  for (size_t i = 1; i < validIntervals.size(); ++i) {
    if (validIntervals[i].first > validIntervals[i - 1].second) {
      // Gap found: intervals are disconnected
      return InvalidBracket;
    }
  }

  // Collect unique endpoints
  std::array<int32_t, 2 * N> ends;
  size_t count = 0;
  for (const auto& inter : validIntervals) {
    ends[count++] = inter.first;
    ends[count++] = inter.second;
  }

  std::sort(ends.begin(), ends.begin() + count);
  auto last = std::unique(ends.begin(), ends.begin() + count);
  count = last - ends.begin();

  int32_t bestLength = std::numeric_limits<int32_t>::max();
  Bracket bestBracket{InvalidBracket};

  for (size_t i = 0; i < count; ++i) {
    const int32_t L = ends[i];
    for (size_t j = i + 1; j < count; ++j) {
      const int32_t R = ends[j];
      int32_t length = R - L;
      if (length >= bestLength) {
        break;
      }

      bool overlaps = true;
      for (const auto& inter : validIntervals) {
        if (L >= inter.second || R <= inter.first) {
          overlaps = false;
          break;
        }
      }

      if (overlaps) {
        bestLength = length;
        bestBracket = {L, R};
        break;
      }
    }
  }

  return bestBracket;
}

} // namespace o2::its::utils
