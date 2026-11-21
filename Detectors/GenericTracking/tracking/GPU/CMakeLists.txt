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

# The device tracker for CUDA and, hipified from the same sources, for HIP:
# one shared library per backend, each unit whole (no relocatable device
# code: the tracking's device code is all in its headers), which O2's GPU
# backend links and asks for the runtime
# (GPUReconstruction::GetGenericTrackingRuntime).

if(NOT CUDA_ENABLED AND NOT HIP_ENABLED)
  return()
endif()

# gtrk_device_library(<name> <CUDA|HIP> <source.cu>...): an O2 library of
# the sources as written for CUDA, hipified for HIP; its target in
# <name>_TARGET.
function(gtrk_device_library name backend)
  if(backend STREQUAL "HIP")
    o2_add_hipified_library(${name} SOURCES ${ARGN} TARGETVARNAME target)
  else()
    o2_add_library(${name} SOURCES ${ARGN} TARGETVARNAME target)
  endif()
  set(${name}_TARGET ${target} PARENT_SCOPE)
endfunction()

# gtrk_device_target(<target> <CUDA|HIP> <SHARED|OBJECT> [<exported symbol pattern>...]):
# a library of device code for one backend. The CUDA runtime is the shared
# one, so every library of a process talks to the same runtime state;
# cross-execution-space calls are errors. A shared library exports only its
# entries.
set(GTRK_DEVICE_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/cuda" CACHE INTERNAL "")
function(gtrk_device_target target backend type)
  # The device headers next to the CUDA sources: a hipified .hip lives in the
  # build tree, so a quoted include does not find them beside it.
  target_include_directories(${target} PRIVATE "${CMAKE_SOURCE_DIR}" "${GTRK_DEVICE_HEADERS}")
  string(TOLOWER "${backend}" backend_namespace)
  target_compile_definitions(${target} PRIVATE GTRK_BACKEND_NAMESPACE=${backend_namespace} GTRK_BACKEND_${backend})
  set_target_gpu_arch("${backend}" ${target})
  # The architecture's parameters (GPUCA_WARP_SIZE): GPUDefParametersWrapper.h.
  target_link_libraries(${target} PRIVATE O2::GPUTracking)
  set_target_properties(${target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
  if(type STREQUAL "SHARED")
    set(version_script "${CMAKE_CURRENT_BINARY_DIR}/${target}.exports")
    list(JOIN ARGN "; " exported)
    if(exported)
      file(GENERATE OUTPUT "${version_script}" CONTENT "{ global: ${exported}; local: *; };\n")
    else()
      file(GENERATE OUTPUT "${version_script}" CONTENT "{ local: *; };\n")
    endif()
    target_link_options(${target} PRIVATE "LINKER:--version-script=${version_script}")
    set_property(TARGET ${target} APPEND PROPERTY LINK_DEPENDS "${version_script}")
  endif()
  if(backend STREQUAL "CUDA")
    set_target_properties(${target} PROPERTIES CUDA_RUNTIME_LIBRARY Shared)
    target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Werror=cross-execution-space-call>)
  else()
    # FindO2GPU compiles HIP with -ffast-math outside the deterministic mode,
    # and that includes -ffinite-math-only: comparisons with an infinite
    # start value or a NaN guard become undefined. Only the finite-math
    # assumption is taken out (CUDA's -use_fast_math makes none).
    target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:HIP>:-fno-finite-math-only>)
    # hipCUB, CUB's interface on rocPRIM, for the scans.
    target_link_libraries(${target} PRIVATE hip::hipcub)
    if(type STREQUAL "SHARED")
      target_link_libraries(${target} PRIVATE hip::host)
    endif()
  endif()
endfunction()

foreach(backend IN ITEMS CUDA HIP)
  if(NOT ${backend}_ENABLED)
    continue()
  endif()
  string(TOLOWER ${backend} lower)
  string(LENGTH ${lower} lower_length)

  # One device tracker per layout and unit, compiled in parallel.
  gtrk_device_library(GTrkTracking${backend} ${backend} cuda/DeviceRuntime.cu cuda/DeviceTrackerITS.cu cuda/DeviceTrackerMFT.cu cuda/DeviceTrackerTRKFT3.cu cuda/DeviceTrackerTRKFT3TOF.cu cuda/DeviceDecayFitterTRKFT3.cu)
  set(target ${GTrkTracking${backend}_TARGET})
  target_link_libraries(${target} PUBLIC O2::GTrkTracking)
  # Register caps of the CUDA units. Fewer registers per thread leave more
  # warps resident to hide the memory waits of the fits and the follow; each
  # cap is where its unit runs fastest. nvcc applies one cap to every
  # architecture a unit is compiled for. HIP has none: there the follow
  # already runs at the register file's limit, and a lower one only spills.
  if(backend STREQUAL "CUDA" AND NOT GPUCA_DETERMINISTIC_MODE GREATER_EQUAL ${GPUCA_DETERMINISTIC_MODE_MAP_NO_FAST_MATH})
    set_source_files_properties(cuda/DeviceTrackerITS.cu cuda/DeviceTrackerTRKFT3.cu cuda/DeviceTrackerTRKFT3TOF.cu PROPERTIES COMPILE_OPTIONS -maxrregcount=96)
    set_source_files_properties(cuda/DeviceTrackerMFT.cu PROPERTIES COMPILE_OPTIONS -maxrregcount=128)
  endif()
  gtrk_device_target(${target} ${backend} SHARED "_ZN2o24gtrk3gpu${lower_length}${lower}11makeRuntimeEi")
endforeach()
