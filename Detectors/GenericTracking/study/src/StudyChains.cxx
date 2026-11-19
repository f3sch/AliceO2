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

#include "GenericTracking/Study/StudyChains.h"

#include "GenericTracking/Study/TruthStudies.h"

#include <charconv>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>

namespace o2::gtrk::study
{
namespace
{

/// \a text split at \a separator, empty pieces left out.
std::vector<std::string_view> split(std::string_view text, char separator)
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

template <typename T>
T number(std::string_view key, std::string_view value)
{
  T result{};
  const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
  if (error != std::errc{} || end != value.data() + value.size()) {
    throw std::invalid_argument(std::format("truth selection: {}={} is not a number", key, value));
  }
  return result;
}

} // namespace

ChainSelection parseChainSelection(std::string_view text)
{
  ChainSelection selection;
  for (const auto item : split(text, ';')) {
    const auto equals = item.find('=');
    const auto key = item.substr(0, equals);
    const auto value = equals == std::string_view::npos ? std::string_view{} : item.substr(equals + 1U);
    if (key == "findable" && equals == std::string_view::npos) {
      selection.findable = true;
      continue;
    }
    if (value.empty()) {
      throw std::invalid_argument(std::format("truth selection: {} has no value", item));
    }
    if (key == "origins") {
      selection.origins = 0U;
      for (const auto origin : split(value, '+')) {
        const auto bit = origin == "primary" ? ParticleOrigin::kPhysicalPrimary : origin == "decay"    ? ParticleOrigin::kDecay
                                                                                : origin == "material" ? ParticleOrigin::kMaterial
                                                                                                       : ParticleOrigin::kCount;
        if (bit == ParticleOrigin::kCount) {
          throw std::invalid_argument(std::format("truth selection: origin {}: primary, decay or material", origin));
        }
        selection.origins |= 1U << static_cast<unsigned>(bit);
      }
    } else if (key == "pdg") {
      for (const auto code : split(value, '+')) {
        const auto pdg = number<std::int64_t>(key, code);
        if (pdg == 0 || std::abs(pdg) > std::numeric_limits<std::int32_t>::max()) {
          throw std::invalid_argument(std::format("truth selection: pdg {} is no PDG code", code));
        }
        selection.pdgs.push_back(static_cast<std::int32_t>(std::abs(pdg)));
      }
    } else if (key == "minPt") {
      selection.minPtGeV = number<double>(key, value);
      if (!std::isfinite(selection.minPtGeV) || selection.minPtGeV < 0.) {
        throw std::invalid_argument(std::format("truth selection: minPt={} is no pT", value));
      }
    } else if (key == "minRows") {
      selection.minRows = number<std::uint32_t>(key, value);
    } else if (key == "layers") {
      selection.layers = number<std::uint64_t>(key, value);
    } else {
      throw std::invalid_argument(std::format("truth selection: {}: origins, pdg, minPt, minRows, layers or findable", item));
    }
  }
  return selection;
}

TruthStudySpec parseTruthStudy(std::string_view text)
{
  // The path is all after the second colon, colons included.
  const auto first = text.find(':');
  const auto second = first == std::string_view::npos ? std::string_view::npos : text.find(':', first + 1U);
  if (second == std::string_view::npos || first == 0U || second == first + 1U || second + 1U == text.size()) {
    throw std::invalid_argument(std::format("truth study {}: <module>:<source>:<path>", text));
  }
  const auto source = text.substr(first + 1U, second - first - 1U);
  TruthStudySpec spec{.module = std::string{text.substr(0, first)}, .path = std::string{text.substr(second + 1U)}};
  if (source == "truth") {
    spec.source = ChainSource::kTruth;
  } else if (source == "tracks") {
    spec.source = ChainSource::kTracks;
  } else {
    throw std::invalid_argument(std::format("truth study {}: the sources are truth and tracks", text));
  }
  return spec;
}

std::string describe(const ChainCounts& counts, std::string_view label)
{
  return std::format("{}: {} particles selected, {} chains; left out: {} untimed, {} first leg inward, {} short, {} without truth", label, counts.selected, counts.chains, counts.untimed, counts.inward, counts.tooShort, counts.noTruth);
}

} // namespace o2::gtrk::study
