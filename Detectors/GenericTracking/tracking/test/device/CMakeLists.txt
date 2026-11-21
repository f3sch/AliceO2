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

# The device tests on each backend, labelled gpu because they need one.

if(NOT CUDA_ENABLED AND NOT HIP_ENABLED)
  return()
endif()

set(fixtures "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/.." "${CMAKE_CURRENT_SOURCE_DIR}/../../../study/include")

function(gtrk_device_test name source backend)
  string(TOLOWER ${backend} lower)
  o2_add_test(${name}-${lower} SOURCES ${source} COMPONENT_NAME generic-tracking NAME generic-tracking-${name}-${lower}
              PUBLIC_LINK_LIBRARIES O2::GTrkTracking${backend} O2::GTrkDeviceTest${backend} O2::GTrkTracking O2::GPUTracking
              LABELS generic-tracking gpu ${ARGN} TARGETVARNAME target)
  if(target)
    target_include_directories(${target} PRIVATE ${fixtures})
    target_compile_definitions(${target} PRIVATE GTRK_TEST_${backend})
  endif()
endfunction()

foreach(backend IN ITEMS CUDA HIP)
  if(NOT ${backend}_ENABLED)
    continue()
  endif()
  string(TOLOWER ${backend} lower)

  gtrk_device_library(GTrkDeviceTest${backend} ${backend} DeviceTestKernels.cu)
  set(target ${GTrkDeviceTest${backend}_TARGET})
  target_include_directories(${target} PRIVATE ${fixtures})
  target_link_libraries(${target} PRIVATE O2::GTrkTracking)
  gtrk_device_target(${target} ${backend} SHARED "_ZN2o24gtrk*testKernels${backend}E*")

  gtrk_device_test(device-runtime testDeviceRuntime.cxx ${backend})
  gtrk_device_test(gpu-stage testGpuStage.cxx ${backend})
  gtrk_device_test(flat-stage testFlatStage.cxx ${backend})
  gtrk_device_test(device-tracker testDeviceTracker.cxx ${backend})
  gtrk_device_test(slice-glue testSliceGlue.cxx ${backend})
  gtrk_device_test(device-storage testDeviceStorage.cxx ${backend})
  gtrk_device_test(road-glue testRoadGlue.cxx ${backend})
  gtrk_device_test(selection-glue testSelectionGlue.cxx ${backend})
endforeach()
