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

# ITS3's workflow: its chips as our sensors, O2's ITS3 clusters as ours, the
# converter, the measurement producer and the tracker (its unit instantiating
# the tracker), with the ITS's conditions (WorkflowITSMFT/ConditionsITSMFT.h).
# It reads O2's ITS3 clusters as workflow-itsmft reads the ITS's
# (ITSMFTInput.h, a header of that library's own).
set(trackerUnits src/TrackerITS3.cxx)
o2_add_library(GTrkWorkflowITS3
               SOURCES src/ITS3.cxx src/Clusters.cxx ${trackerUnits}
               PUBLIC_LINK_LIBRARIES O2::GTrkWorkflowITSMFT O2::ITSBase
               PRIVATE_LINK_LIBRARIES O2::ITS3Base O2::ITS3Reconstruction O2::DataFormatsITSMFT O2::SimulationDataFormat O2::MathUtils
               TARGETVARNAME targetName)

# The tracker's unit instantiates the header tracker (GTRK_TRACKER_OPTIMIZATION).
set_source_files_properties(${trackerUnits} PROPERTIES COMPILE_OPTIONS "${GTRK_TRACKER_OPTIMIZATION}")

target_include_directories(${targetName} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../workflow-itsmft/src)

# o2_add_library installs the headers of include/<library name> only.
target_include_directories(${targetName} PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(DIRECTORY include/GenericTracking DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

o2_add_executable(its3-reco-workflow
                  SOURCES src/gtrk-its3-reco-workflow.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkWorkflowITS3 O2::DetectorsRaw O2::ITSMFTWorkflow)
