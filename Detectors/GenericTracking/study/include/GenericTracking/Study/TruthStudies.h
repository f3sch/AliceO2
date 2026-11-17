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

/// \file TruthStudies.h
/// \brief The truth studies a host runs: each a module, a source of chains
///        and an output, as `<module>:<source>:<path>` names them.
///
/// Modules: `gsf` (GSFStudy.h), `species` (SpeciesFitStudy.h), `clusters`
/// (ClusterPullStudy.h), `decaycensus`
/// (DecayCensusStudy.h). Sources: `truth`, the particles' own rows,
/// and `tracks`, the tracker's good tracks (StudyChains.h). The path's
/// extension picks the format (RecordWriter.h). The bench and the sink both
/// hold one of these for the run and hand it each timeframe.

#ifndef GTRK_STUDY_TRUTHSTUDIES_H_
#define GTRK_STUDY_TRUTHSTUDIES_H_

#include "GenericTracking/Study/ClusterPullStudy.h"
#include "GenericTracking/Study/DecayCensusStudy.h"
#include "GenericTracking/Study/GSFStudy.h"
#include "GenericTracking/Study/RecordWriter.h"
#include "GenericTracking/Study/SpeciesFitStudy.h"
#include "GenericTracking/Study/StudyChains.h"
#include "GenericTracking/Study/TruthEvent.h"

#include <array>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace o2::gtrk::study
{

/// A study as `<module>:<source>:<path>` names it.
struct TruthStudySpec {
  std::string module;
  ChainSource source = ChainSource::kTruth;
  std::string path;
};

/// \a text as a study; throws std::invalid_argument if it names no module,
/// source or path.
[[nodiscard]] TruthStudySpec parseTruthStudy(std::string_view text);

template <DetectorLayout Det>
class TruthStudies
{
 public:
  /// The studies \a specs name, every one's chains taken by \a selection.
  /// Throws std::invalid_argument for an unknown module, an output of another
  /// format or one named twice, before any output is opened.
  TruthStudies(std::span<const TruthStudySpec> specs, ChainSelection selection) : mSelection(std::move(selection))
  {
    std::vector<std::string_view> paths;
    for (const auto& spec : specs) {
      if (spec.module != "gsf" && spec.module != "species" && spec.module != "decaycensus" && spec.module != "clusters") {
        throw std::invalid_argument("truth study " + spec.module + ": the modules are gsf, species, decaycensus and clusters");
      }
      if (!RecordWriter::writes(spec.path) || std::ranges::find(paths, spec.path) != paths.end()) {
        throw std::invalid_argument("truth study " + spec.module + ": its output " + spec.path + " is not a new .csv or .root path");
      }
      paths.push_back(spec.path);
    }
    for (const auto& spec : specs) {
      auto& study = mStudies.emplace_back(std::make_unique<Study>());
      study->source = spec.source;
      study->writer = std::make_unique<RecordWriter>(spec.path);
      if (spec.module == "gsf") {
        study->module.template emplace<GSFStudy<Det>>(*study->writer);
      } else if (spec.module == "species") {
        study->module.template emplace<SpeciesFitStudy<Det>>(*study->writer);
      } else if (spec.module == "clusters") {
        study->module.template emplace<ClusterPullStudy<Det>>(*study->writer);
      } else {
        study->module.template emplace<DecayCensusStudy<Det>>(*study->writer);
      }
    }
  }

  [[nodiscard]] bool empty() const { return mStudies.empty(); }
  /// Whether a study takes chains from \a source.
  [[nodiscard]] bool takes(ChainSource source) const
  {
    return std::ranges::any_of(mStudies, [&](const auto& study) { return study->source == source; });
  }

  /// The truth-chain studies on \a event.
  void run(const TruthEvent<Det>& event)
  {
    if (!takes(ChainSource::kTruth)) {
      return;
    }
    const auto chains = truthChains(event, mSelection, mCounts[0]);
    dispatch(event, chains, ChainSource::kTruth);
  }

  /// The track-chain studies on \a event's \a tracks.
  void run(const TruthEvent<Det>& event, const typename TrackRefit<Det>::Tracks& tracks)
  {
    if (!takes(ChainSource::kTracks)) {
      return;
    }
    const auto chains = trackChains(event, tracks, mSelection, mCounts[1]);
    dispatch(event, chains, ChainSource::kTracks);
  }

  /// What every source took and every study did, in log lines.
  [[nodiscard]] std::vector<std::string> summary() const
  {
    std::vector<std::string> lines;
    if (takes(ChainSource::kTruth)) {
      lines.push_back(describe(mCounts[0], "truth chains"));
    }
    if (takes(ChainSource::kTracks)) {
      lines.push_back(describe(mCounts[1], "track chains"));
    }
    for (const auto& study : mStudies) {
      std::visit([&](const auto& module) {
        if constexpr (!std::same_as<std::decay_t<decltype(module)>, std::monostate>) {
          for (auto& line : module.summary()) {
            lines.push_back(std::move(line));
          }
        }
      },
                 study->module);
    }
    return lines;
  }

 private:
  struct Study {
    ChainSource source = ChainSource::kTruth;
    std::unique_ptr<RecordWriter> writer;
    std::variant<std::monostate, GSFStudy<Det>, SpeciesFitStudy<Det>, DecayCensusStudy<Det>, ClusterPullStudy<Det>> module;
  };

  void dispatch(const TruthEvent<Det>& event, std::span<const Chain> chains, ChainSource source)
  {
    for (auto& study : mStudies) {
      if (study->source != source) {
        continue;
      }
      std::visit([&](auto& module) {
        if constexpr (!std::same_as<std::decay_t<decltype(module)>, std::monostate>) {
          module.process(event, chains, source);
        }
      },
                 study->module);
    }
  }

  ChainSelection mSelection;
  std::vector<std::unique_ptr<Study>> mStudies;
  std::array<ChainCounts, 2> mCounts{}; ///< [truth, tracks]
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_TRUTHSTUDIES_H_
