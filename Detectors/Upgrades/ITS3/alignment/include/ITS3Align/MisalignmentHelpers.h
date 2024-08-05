// Copyright 2020-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file MisalignmentHelpers.h
/// \brief Misalign helpers
/// \author felix.schlepper@cern.ch

#ifndef ITS3_MISALIGNMENT_HELPERS_H_
#define ITS3_MISALIGNMENT_HELPERS_H_

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "ITSBase/GeometryTGeo.h"
#include "ITSMFTSimulation/Hit.h"
#include "ITS3Base/SpecsV2.h"
#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"

#include "Framework/Logger.h"
#include "fmt/format.h"
#include "boost/property_tree/ptree.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "Math/Transform3D.h"
#include "Math/Translation3D.h"
#include "Math/Rotation3D.h"
#include "Math/EulerAngles.h"
#include "Math/PositionVector3D.h"

#include <string>
#include <filesystem>
#include <memory>
#include <algorithm>
#include <optional>
#include <vector>
#endif

namespace o2::its3::align
{

namespace fs = std::filesystem;
namespace pt = boost::property_tree;
namespace its3c = o2::its3::constants;
using Point3D = ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>, ROOT::Math::DefaultCoordinateSystemTag>;
using Trans3D = ROOT::Math::Translation3DF;
using Rot3D = ROOT::Math::Rotation3D;
using Euler3D = ROOT::Math::EulerAngles;
using Trafo3D = ROOT::Math::Transform3DF;

#define DECLARE_SENSOR(id)       \
  float Sensor##id##Dx = 0.f;    \
  float Sensor##id##Dy = 0.f;    \
  float Sensor##id##Dz = 0.f;    \
  float Sensor##id##Phi = 0.f;   \
  float Sensor##id##Theta = 0.f; \
  float Sensor##id##Psi = 0.f;

struct MisAlignGlobalParams : public o2::conf::ConfigurableParamHelper<MisAlignGlobalParams> {
  DECLARE_SENSOR(0)
  DECLARE_SENSOR(1)
  DECLARE_SENSOR(2)
  DECLARE_SENSOR(3)
  DECLARE_SENSOR(4)
  DECLARE_SENSOR(5)

  O2ParamDef(MisAlignGlobalParams, "MisAlignGlobalParams");
};
O2ParamImpl(MisAlignGlobalParams);

void backup_file(const fs::path& src, const fs::path& dest)
{
  if (fs::exists(dest)) {
    LOGP(info, "Previous orignal file found, using this as src");
  } else {
    if (!fs::exists(src)) {
      LOGP(fatal, "File {} does not exist", src.c_str());
    }
    LOGP(info, "Trying to backup file to {}", dest.c_str());
    try {
      fs::rename(src, dest);
    } catch (const fs::filesystem_error& err) {
      LOGP(fatal, "Cannot create backup file for Hit-File: {}", err.what());
    }
  }
}

std::string append_stem(const std::string& filename, const std::string& add)
{
  fs::path filepath{filename};
  auto stem = filepath.stem().string();
  auto extension = filepath.extension().string();
  return stem + add + extension;
}

// Utility function to split a string by a delimiter
std::vector<std::string> split(const std::string& s, char delimiter = '/')
{
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)) {
    if (!token.empty()) {
      tokens.push_back(token);
    }
  }
  return tokens;
}

inline void navigate(const std::string& path)
{
  if (!gGeoManager->cd(path.c_str())) {
    LOGP(fatal, "Cannot navigate to {}", path);
  }
}

std::string composePathSensor(int sensor)
{
  const int layerID{sensor / 2};
  const int sensorID{sensor % 2};
  return fmt::format("/cave/barrel_1/ITSV_2/ITSUWrapVol0_1/ITS3Layer{}_0/ITS3CarbonForm{}_{}",
                     layerID, layerID, sensorID);
}

void ApplyGlobalMatrixVolume(const std::string& path, const TGeoHMatrix& globalMatrix)
{
  gGeoManager->CdTop();
  TGeoHMatrix* pgMatrix{nullptr};
  TGeoHMatrix gAccMatrix;
  std::string curPath{};
  for (const auto& comp : split(path)) {
    curPath += "/" + comp;
    navigate(curPath);
    pgMatrix = gGeoManager->GetCurrentMatrix();
    gAccMatrix.Multiply(pgMatrix);
  }
  navigate(path);
  auto node = gGeoManager->GetCurrentNode();
  if (node == nullptr) {
    LOGP(fatal, "Nullptr for node at {}", path);
  }
  auto motherVol = node->GetMotherVolume();
  if (motherVol == nullptr) {
    LOGP(fatal, "Nullptr for motherVol at {}", path);
  }
  // Compute the inverse of the accumulated global transformation matrix
  auto gAccMatrix1 = gAccMatrix.Inverse();
  // Compute the relative transformation matrix for the volume
  auto relativeMatrix = globalMatrix;
  relativeMatrix.MultiplyLeft(gAccMatrix1);

  auto nodemat = dynamic_cast<TGeoNodeMatrix*>(node);
  nodemat->SetMatrix(new TGeoHMatrix(globalMatrix));

  // Force the container volume of the object to update itself
  motherVol->Voxelize("");
}

} // namespace o2::its3::align

#endif
