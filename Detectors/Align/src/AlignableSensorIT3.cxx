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

/// @file   AlignableSensorIT3.h
/// @brief  IT3 sensor

#include "ITSBase/GeometryTGeo.h"
#include "Align/AlignableSensorIT3.h"
#include "Align/utils.h"
#include "Framework/Logger.h"
#include "Align/AlignmentPoint.h"
#include "Align/AlignableDetector.h"
#include "ITS3Base/SpecsV2.h"

ClassImp(o2::align::AlignableSensorIT3);

using namespace o2::align::utils;
using namespace TMath;

namespace o2::align
{

//_________________________________________________________
AlignableSensorIT3::AlignableSensorIT3(const char* name, int vid, int iid, Controller* ctr)
  : AlignableSensor(name, vid, iid, ctr){}

//____________________________________________
void AlignableSensorIT3::prepareMatrixT2L()
{
  // extract geometry T2L matrix
  TGeoHMatrix t2l;
  const auto& l2g = getMatrixL2GIdeal();
  double locA[3] = {-100., 0., 0.}, locB[3] = {100., 0., 0.}, gloA[3], gloB[3];

  if(o2::its3::constants::detID::isDetITS3(getVolID())){
    int iLayer = o2::its3::constants::detID::getDetID2Layer(getVolID());
    // We need to calcualte the line tangent at the mid-point in the geometry
    const auto radius = o2::its3::constants::radii[iLayer];
    const auto phi1 = o2::its3::constants::tile::width / radius;
    const auto phi2 = o2::its3::constants::pixelarray::width / radius + phi1;
    const auto phi3 = (phi2 - phi1) / 2.; // mid-point in phi
    const auto x = radius * std::cos(phi3);
    const auto y = radius * std::sin(phi3);
    // For the tangent we make the parametric line equation y = m * x - c
    const auto m = x / y;
    const auto c = y - m * x;
    // Now we can given any x calulate points along this line, we pick points far away,
    // the calculation of the normal should work then below.
    locA[1] = m * locA[0] + c;
    locB[1] = m * locB[0] + c;
  }

  l2g.LocalToMaster(locA, gloA);
  l2g.LocalToMaster(locB, gloB);
  double dx = gloB[0] - gloA[0], dy = gloB[1] - gloA[1];
  double t = (gloB[0] * dx + gloB[1] * dy) / (dx * dx + dy * dy);
  double xp = gloB[0] - dx * t, yp = gloB[1] - dy * t;
  mX = std::sqrt(xp * xp + yp * yp);
  float alp = std::atan2(yp, xp);
  o2::math_utils::bringTo02Pi(alp);
  mAlp = alp;
  t2l.RotateZ(mAlp * RadToDeg()); // rotate in direction of normal to the sensor plane
  const TGeoHMatrix l2gi = l2g.Inverse();
  t2l.MultiplyLeft(&l2gi);
  setMatrixT2L(t2l);
}

} // namespace o2
