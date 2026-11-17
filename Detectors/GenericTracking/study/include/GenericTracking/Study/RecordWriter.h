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

/// \file RecordWriter.h
/// \brief A study's records into a CSV or a ROOT tree, from one list of
///        columns.
///
/// A record is a flat struct of scalars that names its columns once:
///
///     struct Record {
///       static constexpr std::string_view kName = "gsf";
///       std::uint64_t particle = 0U;
///       double pTrue = 0.;
///       template <typename Visit>
///       void columns(const Visit& visit) const { visit("particle", particle); visit("pTrue", pTrue); }
///     };
///
/// `columns` hands the visitor the record's own members, in their own types:
/// the tree keeps their addresses until the entry is filled. A `.csv` path
/// gets a header from the names and a line per record, integers (bytes and
/// flags included) as decimals and floating point at round-trip precision,
/// `nan`, `inf` and `-inf` spelled out; a `.root` path a tree named `kName`,
/// a branch per column. A writer takes one record type, which `begin` names
/// before any record (the CSV's header, written even if no record follows),
/// and lives for a run.

#ifndef GTRK_STUDY_RECORDWRITER_H_
#define GTRK_STUDY_RECORDWRITER_H_

#include "CommonUtils/TreeStreamRedirector.h"

#include <cmath>
#include <concepts>
#include <cstdint>
#include <format>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

namespace o2::gtrk::study
{

class RecordWriter
{
 public:
  /// A writer into \a path, a `.csv` or a `.root` file, created anew.
  /// Throws std::invalid_argument for another extension, std::runtime_error
  /// if the file cannot be written.
  explicit RecordWriter(const std::string& path);
  ~RecordWriter();
  RecordWriter(const RecordWriter&) = delete;
  RecordWriter& operator=(const RecordWriter&) = delete;
  RecordWriter(RecordWriter&&) = delete;
  RecordWriter& operator=(RecordWriter&&) = delete;

  /// Whether \a path names a format a writer writes.
  [[nodiscard]] static bool writes(std::string_view path) { return path.ends_with(".csv") || path.ends_with(".root"); }

  [[nodiscard]] const std::string& path() const { return mPath; }

  /// The writer's record type, once, before any record. Throws
  /// std::logic_error for a second call.
  template <typename Record>
  void begin();

  /// \a record, of the type `begin` named.
  template <typename Record>
  void write(const Record& record);

 private:
  /// \a value as the CSV writes it.
  template <typename T>
  static void appendValue(std::string& line, const T& value);

  std::string mPath;
  std::unique_ptr<std::ofstream> mCsv;
  std::string_view mRecord; ///< the record type's `kName`, once `begin` named it
  std::unique_ptr<o2::utils::TreeStreamRedirector> mTree;
  std::string mLine;
};

template <typename T>
void RecordWriter::appendValue(std::string& line, const T& value)
{
  if constexpr (std::same_as<T, bool>) {
    line += value ? '1' : '0';
  } else if constexpr (std::is_integral_v<T>) {
    std::format_to(std::back_inserter(line), "{}", static_cast<std::conditional_t<std::is_signed_v<T>, std::int64_t, std::uint64_t>>(value));
  } else {
    static_assert(std::is_floating_point_v<T>, "a record column is a scalar");
    if (std::isnan(value)) {
      line += "nan";
    } else {
      std::format_to(std::back_inserter(line), "{}", value);
    }
  }
}

template <typename Record>
void RecordWriter::begin()
{
  if (!mRecord.empty()) {
    throw std::logic_error("RecordWriter: " + mPath + " already takes " + std::string{mRecord});
  }
  mRecord = Record::kName;
  if (mCsv) {
    mLine.clear();
    Record{}.columns([&](std::string_view name, const auto& /*value*/) {
      mLine += mLine.empty() ? "" : ",";
      mLine += name;
    });
    *mCsv << mLine << '\n';
  }
}

template <typename Record>
void RecordWriter::write(const Record& record)
{
  if (mRecord != Record::kName) {
    throw std::logic_error("RecordWriter: " + mPath + " takes " + std::string{mRecord} + ", not " + std::string{Record::kName});
  }
  if (mCsv) {
    mLine.clear();
    bool first = true;
    record.columns([&](std::string_view /*name*/, const auto& value) {
      if (!first) {
        mLine += ',';
      }
      first = false;
      appendValue(mLine, value);
    });
    *mCsv << mLine << '\n';
    return;
  }
  auto& tree = (*mTree) << Record::kName.data();
  record.columns([&](std::string_view name, const auto& value) {
    const std::string label = std::string{name} + "=";
    tree << label.c_str() << value;
  });
  tree << "\n";
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_RECORDWRITER_H_
