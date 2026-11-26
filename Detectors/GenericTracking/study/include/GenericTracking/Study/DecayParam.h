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

/// \file DecayParam.h
/// \brief The decay finder's selections: GTRKDecay.

#ifndef GTRK_STUDY_DECAYPARAM_H_
#define GTRK_STUDY_DECAYPARAM_H_

#include "CommonUtils/ConfigurableParam.h"
#include "CommonUtils/ConfigurableParamHelper.h"

namespace o2::gtrk::study
{

/// The Lambda candidates a cascade takes (`decayfinder`, CascadeFinder.h):
/// a fitted Lambda or anti-Lambda whose mass lies within `lambdaNSigma` of
/// the Lambda's in its own sigma, whose legs' displacement S^2 against the
/// collision's point are both at least `lambdaMinS2` (a leg without one, in a
/// collision without a point, passes only 0), and whose vertex lies before
/// both legs' first measurements (always: a decay is no later than its
/// daughters' first hits). The rows of a cascade's own mother found from
/// its fit (MotherSearch.h): each within the gate chiSquareQuantile(2,
/// `motherNSigma`). A charm three-prong (D+, Ds+, Lambda_c+: CharmFinder.h)
/// is fitted only where its tracks' mass from their first measurements'
/// momenta under its species lies within `charmPreMassWindow` [GeV] of the
/// hypothesis's mass.
// O2ParamDef deletes the copy constructor and supplies the rest of the boilerplate.
// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions)
struct DecayParam : public o2::conf::ConfigurableParamHelper<DecayParam> {
  double lambdaNSigma = 3.;
  double lambdaMinS2 = 0.;
  double motherNSigma = 3.;
  double charmPreMassWindow = 0.3;

  O2ParamDef(DecayParam, "GTRKDecay")
};

} // namespace o2::gtrk::study

#endif // GTRK_STUDY_DECAYPARAM_H_
