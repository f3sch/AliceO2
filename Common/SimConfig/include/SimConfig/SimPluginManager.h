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

#ifndef SIMPLUGINMANAGER_H_
#define SIMPLUGINMANAGER_H_

#include "CommonUtils/PluginManager.h"

namespace o2::conf
{
class SimPluginManager : public o2::utils::PluginManagerBase<SimPluginManager>
{
  friend class o2::utils::PluginManagerBase<PluginManagerBase>;
};

} // namespace o2::conf

#endif // SIMPLUGINMANAGER_H_
