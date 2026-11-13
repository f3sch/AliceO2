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

# ALICE 3's TRK and FT3 workflow: their sensors, clusters and conditions, the
# converter, the measurement producer and the tracker (its unit instantiating
# the tracker).
set(trackerUnits src/TrackerTRKFT3.cxx)
o2_add_library(GTrkWorkflowTRK
               SOURCES src/TRK.cxx ${trackerUnits}
               PUBLIC_LINK_LIBRARIES O2::GTrkWorkflow O2::DataFormatsTRKFT3 O2::TRKReconstruction O2::TRKBase
               PRIVATE_LINK_LIBRARIES O2::FT3Base O2::Field O2::SimulationDataFormat O2::CommonConstants ROOT::Geom
               TARGETVARNAME targetName)

# The tracker's unit instantiates the header tracker (GTRK_TRACKER_OPTIMIZATION).
set_source_files_properties(${trackerUnits} PROPERTIES COMPILE_OPTIONS "${GTRK_TRACKER_OPTIMIZATION}")

# o2_add_library installs the headers of include/<library name> only.
target_include_directories(${targetName} PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(DIRECTORY include/GenericTracking DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

o2_add_executable(trk-reco-workflow
                  SOURCES src/gtrk-trk-reco-workflow.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkWorkflowTRK O2::DetectorsRaw O2::TRKWorkflow)

add_subdirectory(test)
