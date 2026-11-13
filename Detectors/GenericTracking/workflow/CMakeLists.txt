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

# The detector-independent part of the workflows: the parameters, the
# placement's frame, the clusters' common steps, the conditions' declaration,
# and the devices every detector instantiates in its own library
# (MeasurementProducer.h, TrackerDevice.h, RecoWorkflow.h): the ITS's and
# MFT's in workflow-itsmft, ALICE 3's in workflow-trk.
o2_add_library(GTrkWorkflow
               SOURCES src/Placement.cxx src/RecoParam.cxx
               PUBLIC_LINK_LIBRARIES O2::GTrkTracking O2::Framework O2::CommonUtils O2::DataFormatsCalibration O2::DetectorsBase
                                     O2::DetectorsRaw O2::MathUtils ROOT::Core
               PRIVATE_LINK_LIBRARIES O2::CommonConstants O2::ReconstructionDataFormats
               TARGETVARNAME targetName)

# o2_add_library installs the headers of include/<library name> only.
target_include_directories(${targetName} PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(DIRECTORY include/GenericTracking DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# GTRKReco, a ConfigurableParam, which needs a dictionary for
# --configKeyValues to reach its fields by name, and the cluster-shape band,
# a CCDB object.
o2_target_root_dictionary(GTrkWorkflow
                          HEADERS include/GenericTracking/Workflow/RecoParam.h
                                  include/GenericTracking/Workflow/ClusterShapeBand.h
                          LINKDEF src/GTrkWorkflowLinkDef.h)

add_subdirectory(test)
