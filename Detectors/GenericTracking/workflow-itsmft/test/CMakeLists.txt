# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

function(gtrk_workflow_itsmft_test name source)
  o2_add_test(${name} SOURCES ${source} COMPONENT_NAME generic-tracking
              PUBLIC_LINK_LIBRARIES O2::GTrkWorkflowITSMFT O2::GTrkTracking ${ARGN}
              LABELS generic-tracking TARGETVARNAME target)
  if(target)
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../tracking/test)
  endif()
endfunction()

gtrk_workflow_itsmft_test(alice-layouts testAliceLayouts.cxx O2::ITSBase O2::MFTBase O2::MathUtils)
gtrk_workflow_itsmft_test(alice-clusters testAliceClusters.cxx O2::DataFormatsITSMFT O2::CommonConstants O2::CommonUtils)
gtrk_workflow_itsmft_test(disc-chain-pulls testDiscChainPulls.cxx O2::MFTBase O2::MathUtils)
