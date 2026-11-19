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

/// \file Statistics.h
/// \brief The summaries of a sample the studies, tools and tests share: the
///        median, a quantile, robust widths and shares. Host only.

#ifndef GTRK_BASE_STATISTICS_H_
#define GTRK_BASE_STATISTICS_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <ranges>
#include <utility>
#include <vector>

namespace o2::gtrk
{

/// The median of \a values, the upper of the two middle ones for an even
/// count, 0 for none. Reorders \a values (nth_element); pass a copy to keep
/// their order.
template <std::ranges::random_access_range R>
[[nodiscard]] std::ranges::range_value_t<R> median(R&& values)
{
  using T = std::ranges::range_value_t<R>;
  if (std::ranges::empty(values)) {
    return T{};
  }
  const auto middle = std::ranges::begin(values) + static_cast<std::ptrdiff_t>(std::ranges::size(values) / 2U);
  std::ranges::nth_element(std::forward<R>(values), middle);
  return *middle;
}

/// The element at rank floor(\a q (n - 1)) of \a values (n of them, not
/// empty: the caller says what an empty set gives), \a q in [0, 1]: 0 the
/// smallest, 1 the largest. Reorders \a values (nth_element).
template <std::ranges::random_access_range R>
[[nodiscard]] std::ranges::range_value_t<R> quantile(R&& values, double q)
{
  const auto at = std::ranges::begin(values) + static_cast<std::ptrdiff_t>(q * static_cast<double>(std::ranges::size(values) - 1U));
  std::ranges::nth_element(std::forward<R>(values), at);
  return *at;
}

/// Half the range of the central 68.27 % of \a values (not empty), between
/// the quantiles a Gaussian has at -1 and +1 sigma: its sigma, blind to a
/// tail. Reorders \a values.
template <std::ranges::random_access_range R>
[[nodiscard]] std::ranges::range_value_t<R> halfWidth68(R&& values)
{
  constexpr double kBelowOneSigma = 0.15865525393145707; // Phi(-1)
  const auto low = quantile(values, kBelowOneSigma);
  return static_cast<std::ranges::range_value_t<R>>(0.5 * (quantile(std::forward<R>(values), 1. - kBelowOneSigma) - low));
}

/// 1.4826 times the median absolute deviation, the factor making it a
/// Gaussian's sigma: a width blind to a tail; 0 for no values. Reorders
/// \a values.
template <std::ranges::random_access_range R>
[[nodiscard]] std::ranges::range_value_t<R> robustSigma(R&& values)
{
  using T = std::ranges::range_value_t<R>;
  const T centre = median(values);
  std::vector<T> deviations;
  deviations.reserve(std::ranges::size(values));
  for (const T value : std::forward<R>(values)) {
    deviations.push_back(std::abs(value - centre));
  }
  return static_cast<T>(1.4826 * median(deviations));
}

/// \a part as a fraction of \a whole, and as a percentage; 0 for none.
[[nodiscard]] inline double fraction(double part, double whole) { return whole == 0. ? 0. : part / whole; }
[[nodiscard]] inline double percent(std::uint64_t part, std::uint64_t whole) { return whole == 0U ? 0. : 100. * static_cast<double>(part) / static_cast<double>(whole); }

} // namespace o2::gtrk

#endif // GTRK_BASE_STATISTICS_H_
