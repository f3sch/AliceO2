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

# The ITS's and the MFT's workflows: their chips as our sensors, O2's clusters
# as ours, their conditions and readout frames, the converter, the measurement
# producer and the tracker of each (the tracker's unit instantiating the
# tracker), the ITS's frame selection, the adaptor into O2's formats, and
# their tracking inside the GPU reconstruction.
set(trackerUnits src/TrackerITS.cxx src/TrackerMFT.cxx)
o2_add_library(GTrkWorkflowITSMFT
               SOURCES src/Clusters.cxx src/Placement.cxx src/ITS.cxx src/MFT.cxx src/ITSROFSelection.cxx src/AdaptorITS.cxx
                       src/AdaptorMFT.cxx src/GPURecoTracking.cxx ${trackerUnits}
               PUBLIC_LINK_LIBRARIES O2::GTrkWorkflow O2::DataFormatsITSMFT O2::ITSBase O2::MFTBase
               PRIVATE_LINK_LIBRARIES O2::ITStracking O2::DataFormatsTRD O2::ITSMFTBase O2::ITSMFTTracking O2::DataFormatsITS
                                      O2::DataFormatsMFT O2::SimulationDataFormat O2::CommonConstants O2::ReconstructionDataFormats
                                      ROOT::Geom ROOT::RIO ROOT::Tree
               TARGETVARNAME targetName)

# The tracker's units instantiate the header tracker (GTRK_TRACKER_OPTIMIZATION).
set_source_files_properties(${trackerUnits} PROPERTIES COMPILE_OPTIONS "${GTRK_TRACKER_OPTIMIZATION}")

# o2_add_library installs the headers of include/<library name> only.
target_include_directories(${targetName} PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(DIRECTORY include/GenericTracking DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# The tracking workflows: O2's clusters in, GenericTracking's tracks and
# vertices out, and O2's formats and files.
o2_add_executable(its-reco-workflow
                  SOURCES src/gtrk-its-reco-workflow.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkWorkflowITSMFT O2::DetectorsRaw O2::ITSMFTWorkflow O2::ITSWorkflow O2::GlobalTrackingWorkflowWriters)
o2_add_executable(mft-reco-workflow
                  SOURCES src/gtrk-mft-reco-workflow.cxx
                  COMPONENT_NAME gtrk
                  PUBLIC_LINK_LIBRARIES O2::GTrkWorkflowITSMFT O2::DetectorsRaw O2::ITSMFTWorkflow O2::ITSMFTCAWriter)

add_subdirectory(test)
