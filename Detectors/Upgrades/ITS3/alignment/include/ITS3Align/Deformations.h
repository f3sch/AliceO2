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

#ifndef ITS3_DEFORMATIONS_H_
#define ITS3_DEFORMATIONS_H_

#include "ITS3Align/MisalignmentParameters.h"
#include "MathUtils/LegendrePols.h"
#include "GPUCommonMath.h"

#include <filesystem>

namespace o2::its3::align
{

class Deformations
{
 public:
  // init deformations from the parameter file
  void init(const std::filesystem::path&);

  std::tuple<double, double, double> getDeformation(unsigned int id, double radius, double phi, double z, double u, double v, double fac = 1.0) const
  {
    // Calculate f_def(phi,z) = ((r+dr)*cos(phi), (r+dr)*sin(phi), z)^T + (dx, dy, dz)^T
    const double dr = mLegendre[id](u, v);
    const double drr = radius + dr * fac;
    double sn, cs;
    o2::gpu::GPUCommonMath::SinCosd(phi, sn, cs);
    const auto& global = mParams.getGlobal(id);
    return {drr * cs + global.getX() * fac, drr * sn + global.getY() * fac, z + global.getZ() * fac};
  }

  double getDeformation(unsigned int id, double u, double v)
  {
    return mLegendre[id](u, v);
  }

  const o2::math_utils::Legendre2DPolynominal& getLegendre(unsigned int id) { return mLegendre[id]; }

  unsigned int getOrder(unsigned int id) { return mLegendre[id].NOrder(); }

 private:
  MisalignmentParameters mParams;

  // Legendre polynominals to model deformations in radius; parameterized by normalized phi (u) and z (v) coordinates
  std::array<o2::math_utils::Legendre2DPolynominal, 6> mLegendre;
};

} // namespace o2::its3::align

#endif
