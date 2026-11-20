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

/// \file OptionText.h
/// \brief Reading the studies' option values: lists and numbers.

#ifndef GTRK_STUDY_OPTIONTEXT_H_
#define GTRK_STUDY_OPTIONTEXT_H_

#include <cctype>
#include <cerrno>
#include <charconv>
#include <cstdlib>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <vector>

namespace o2::gtrk::study
{

/// \a text split at \a separator, empty pieces left out.
[[nodiscard]] inline std::vector<std::string_view> splitList(std::string_view text, char separator)
{
  std::vector<std::string_view> pieces;
  while (!text.empty()) {
    const auto at = text.find(separator);
    const auto piece = text.substr(0, at);
    if (!piece.empty()) {
      pieces.push_back(piece);
    }
    if (at == std::string_view::npos) {
      break;
    }
    text.remove_prefix(at + 1U);
  }
  return pieces;
}

/// The number \a text holds, all of it; none for anything else (an empty
/// text, trailing characters, a sign on an unsigned type, out of range).
template <typename T>
[[nodiscard]] std::optional<T> parseNumber(std::string_view text)
{
  T result{};
  if constexpr (std::is_floating_point_v<T>) {
    // Floating-point std::from_chars is unavailable in libc++ below macOS 26.
    const std::string copy{text};
    char* end = nullptr;
    errno = 0;
    result = static_cast<T>(std::strtod(copy.c_str(), &end));
    // strtod skips leading white space, from_chars does not.
    if (copy.empty() || std::isspace(static_cast<unsigned char>(copy.front())) != 0 || errno != 0 || end != copy.c_str() + copy.size()) {
      return std::nullopt;
    }
  } else {
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), result);
    if (error != std::errc{} || end != text.data() + text.size()) {
      return std::nullopt;
    }
  }
  return result;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_OPTIONTEXT_H_
