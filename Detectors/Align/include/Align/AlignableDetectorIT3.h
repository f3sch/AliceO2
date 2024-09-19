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

/// @file   AlignableDetectorIT3.h
/// @brief  IT3 detector wrapper

#ifndef ALIGNABLEDETECTORIT3_H
#define ALIGNABLEDETECTORIT3_H

#include "Align/AlignableDetector.h"
#include "Align/utils.h"
#include "ReconstructionDataFormats/BaseCluster.h"

namespace o2
{
namespace its3
{
class TopologyDictionary;
}

namespace align
{

class Controller;

class AlignableDetectorIT3 : public AlignableDetector
{
 public:
  //
  using ClusterD = o2::BaseCluster<double>;
  AlignableDetectorIT3() = default;
  AlignableDetectorIT3(AlignableDetectorIT3&&) = delete;
  AlignableDetectorIT3& operator=(AlignableDetectorIT3&&) = delete;
  AlignableDetectorIT3(Controller* ctr);
  ~AlignableDetectorIT3() override = default;

  void defineVolumes() final;
  int processPoints(GIndex gid, int npntCut, bool inv) final;
  bool prepareDetectorData() final;
  void SetAddErrorLr(int ilr, double sigY, double sigZ);
  void SetSkipLr(int ilr);
  //
  void updatePointByTrackInfo(AlignmentPoint* pnt, const trackParam_t* t) const final {}
  void setUseErrorParam(int v = 0) final;
  //
  void setIT3Dictionary(const o2::its3::TopologyDictionary* d) { mITSDict = d; }
  //
  void Print(const Option_t* opt = "") const override;
  //
 protected:
  //
  std::vector<ClusterD> mITSClustersArray;
  const o2::its3::TopologyDictionary* mITSDict{nullptr}; // cluster patterns dictionary
  //
  ClassDefOverride(AlignableDetectorIT3, 1);
};
} // namespace align
} // namespace o2
#endif
