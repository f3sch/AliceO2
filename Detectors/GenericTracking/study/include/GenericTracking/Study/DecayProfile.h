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

/// \file DecayProfile.h
/// \brief A decay vertex's position along a line from its likelihood profile
///        and a prior over cells of the line: a conversion's point, which two
///        nearly parallel daughters measure along the flight only through
///        their curvature (a shift d turns them by rho d), weakly and far from
///        linearly.
///
/// At every cell the fit with the vertex fixed at the cell's middle along the
/// line (`LinePosition`); each cell weighted by its prior probability and the
/// fit's likelihood exp(-chi2 / 2). The result: the fit at the weights' mean;
/// its covariance the fixed fit's plus the weights' spread along the line
/// (each cell counted as uniform over its length) carried to every parameter
/// by its derivative by the position, J sigma^2 J^T; its ndf without the
/// position, which the profile estimated.

#ifndef GTRK_STUDY_DECAYPROFILE_H_
#define GTRK_STUDY_DECAYPROFILE_H_

#include "GenericTracking/Study/DecayConstraints.h"
#include "GenericTracking/Study/DecayVertexFit.h"

#include <array>
#include <cmath>
#include <limits>
#include <span>

namespace o2::gtrk::study
{

/// A stretch of a line where the vertex may be: its middle along the line,
/// half its length, and its prior probability.
struct LineCell {
  double s = 0.;
  double halfLength = 0.;
  double probability = 0.;
};

/// The profile's result: the fit at the mean, the mean and spread along the
/// line.
struct LineProfile {
  DecayVertexFitResult fit;
  double mean = 0.;
  double sigma = 0.;
};

/// \a constraints' fit (\a fitUnder) profiled along the line from \a origin
/// along the unit \a direction over \a cells.
template <typename FitUnder>
[[nodiscard]] LineProfile profileAlongLine(const FitUnder& fitUnder, const DecayConstraints& constraints, const std::array<double, 3>& origin, const std::array<double, 3>& direction, std::span<const LineCell> cells)
{
  LineProfile result;
  struct Point {
    double chi2;
    const LineCell* cell;
  };
  std::vector<Point> points;
  double best = std::numeric_limits<double>::infinity();
  for (const auto& cell : cells) {
    auto fixed = constraints;
    fixed.position = LinePosition{.origin = origin, .direction = direction, .s = cell.s};
    const auto atCell = fitUnder(fixed);
    if (atCell.status == DecayFitStatus::kOk) {
      points.push_back({atCell.chi2, &cell});
      best = std::min(best, atCell.chi2);
    }
  }
  double weight = 0.;
  double mean = 0.;
  double square = 0.;
  for (const auto& point : points) {
    const double w = point.cell->probability * std::exp(-0.5 * (point.chi2 - best));
    weight += w;
    mean += w * point.cell->s;
    square += w * ((point.cell->s * point.cell->s) + (point.cell->halfLength * point.cell->halfLength / 3.));
  }
  if (!(weight > 0.)) {
    result.fit.status = DecayFitStatus::kNoProfile;
    return result;
  }
  result.mean = mean / weight;
  result.sigma = std::sqrt(std::max((square / weight) - (result.mean * result.mean), 0.));
  auto fixed = constraints;
  fixed.position = LinePosition{.origin = origin, .direction = direction, .s = result.mean};
  result.fit = fitUnder(fixed);
  if (result.fit.status == DecayFitStatus::kOk) {
    const auto& j = result.fit.positionSensitivity;
    const std::size_t n = j.size();
    const double variance = result.sigma * result.sigma;
    for (std::size_t a = 0; a < n; ++a) {
      for (std::size_t b = 0; b < n; ++b) {
        result.fit.covariance[(a * n) + b] += variance * j[a] * j[b];
      }
    }
    --result.fit.ndf;
  }
  return result;
}

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYPROFILE_H_
