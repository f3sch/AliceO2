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
# Every unit that instantiates the header tracker -- the tracking device's, the
# study library's, the sink's and the benches' -- compiles at -O3 in
# RelWithDebInfo, where the tracker runs faster than at O2's -O2.
set(GTRK_TRACKER_OPTIMIZATION "$<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANGUAGE:CXX>>:-O3>")

# The targets created in `dir` and its subdirectories.
function(gtrk_targets dir out)
  get_property(targets DIRECTORY ${dir} PROPERTY BUILDSYSTEM_TARGETS)
  get_property(subdirs DIRECTORY ${dir} PROPERTY SUBDIRECTORIES)
  foreach(subdir IN LISTS subdirs)
    gtrk_targets(${subdir} more)
    list(APPEND targets ${more})
  endforeach()
  set(${out} ${targets} PARENT_SCOPE)
endfunction()

# The include directories of the targets outside GenericTracking whose headers
# `target` sees: those of its links and of their links in turn, crossing the
# GenericTracking targets on the way without taking theirs. Read from each
# target's own INTERFACE_INCLUDE_DIRECTORIES, not through
# $<TARGET_PROPERTY:...>: evaluating a link's transitive includes inside the
# target's own marks the targets it passes as seen, and CMake then drops their
# directories from the target's regular includes.
function(gtrk_foreign_includes target root out)
  get_target_property(queue ${target} LINK_LIBRARIES)
  set(seen)
  set(includes)
  while(queue)
    list(POP_FRONT queue link)
    # O2's conditional link, $<$<BOOL:${ENABLE_UPGRADES}>:O2::ITS3Base>: the
    # condition is already a constant here.
    if(link MATCHES "^\\$<\\$<BOOL:([^>]*)>:([^<>$]+)>$")
      if(NOT CMAKE_MATCH_1)
        continue()
      endif()
      set(link ${CMAKE_MATCH_2})
    endif()
    if(NOT TARGET ${link})
      continue() # a plain library, a $<LINK_ONLY:...> or a ::@ directory marker: no headers
    endif()
    get_target_property(aliased ${link} ALIASED_TARGET)
    if(aliased)
      set(link ${aliased})
    endif()
    get_target_property(imported ${link} IMPORTED)
    if(link IN_LIST seen OR imported)
      continue() # an imported target's headers are system headers already
    endif()
    list(APPEND seen ${link})
    get_target_property(dir ${link} SOURCE_DIR)
    if(NOT (dir STREQUAL root OR dir MATCHES "^${root}/"))
      get_target_property(dirs ${link} INTERFACE_INCLUDE_DIRECTORIES)
      foreach(include IN LISTS dirs)
        if(include AND NOT include MATCHES "^\\$<INSTALL_INTERFACE:" AND NOT include MATCHES "/Detectors/GenericTracking/")
          list(APPEND includes "${include}")
        endif()
      endforeach()
    endif()
    get_target_property(links ${link} INTERFACE_LINK_LIBRARIES)
    if(links)
      list(APPEND queue ${links})
    endif()
  endwhile()
  list(REMOVE_DUPLICATES includes)
  set(${out} ${includes} PARENT_SCOPE)
endfunction()

# Every GenericTracking target gets the warnings for its C++ units. The headers
# of the rest of O2 that it includes become system headers for it (O2 and the
# compilers then pass them as -isystem), so that the warnings and the errors
# are those of GenericTracking's code only: O2's own headers are not ours to
# fix. Called at the end of the configuration, when every target that
# GenericTracking links exists.
function(gtrk_strict_warnings root)
  # GenericTracking's own C++ compiles with a strict warning set, as errors.
  # Silent narrowing between the floating-point payloads and the integer indices
  # is the main bug class here, hence -Wconversion and -Wsign-conversion.
  # -Warray-bounds and -Wnull-dereference stay warnings: GCC derives them from
  # its optimiser's view after inlining, which cannot see a caller's invariants
  # (an index below a table's size, a range that is never empty) and changes
  # with the compiler version.
  set(GTRK_WARNINGS -Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor -Woverloaded-virtual -Wold-style-cast -Wcast-align
                    -Wconversion -Wsign-conversion -Wextra-semi -Wimplicit-fallthrough -Wnull-dereference
                    -Wmisleading-indentation -Wformat=2 -Wundef -Werror -Wno-error=array-bounds -Wno-error=null-dereference)
  # GCC-only spellings: Clang does not know them, and under -Werror an unknown
  # warning option is itself an error. Clang knows -Wdouble-promotion and
  # -Wfloat-equal but reads them more broadly (every float-to-double
  # initialisation, every defaulted operator== on a struct with float members);
  # GCC's reading is the one meant.
  set(GTRK_GCC_WARNINGS -Wdouble-promotion -Wfloat-equal -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wuseless-cast)
  gtrk_targets(${root} targets)
  foreach(target IN LISTS targets)
    get_target_property(type ${target} TYPE)
    if(type STREQUAL "UTILITY" OR type STREQUAL "INTERFACE_LIBRARY")
      continue()
    endif()
    target_compile_options(${target} PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:${GTRK_WARNINGS}>"
                                             "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:GNU>>:${GTRK_GCC_WARNINGS}>")
    gtrk_foreign_includes(${target} ${root} includes)
    if(includes)
      target_include_directories(${target} SYSTEM PRIVATE ${includes})
    endif()
    # ROOT's generated dictionary code is not ours either.
    get_target_property(sources ${target} SOURCES)
    foreach(source IN LISTS sources)
      if(source MATCHES "/G__[^/]*\\.cxx$")
        set_source_files_properties(${source} TARGET_DIRECTORY ${target} PROPERTIES COMPILE_OPTIONS -w)
      endif()
    endforeach()
  endforeach()
endfunction()
cmake_language(EVAL CODE "cmake_language(DEFER DIRECTORY [[${CMAKE_SOURCE_DIR}]] CALL gtrk_strict_warnings [[${CMAKE_CURRENT_SOURCE_DIR}]])")

add_subdirectory(base)
add_subdirectory(tracking)
add_subdirectory(workflow)
add_subdirectory(workflow-itsmft)
if(ENABLE_UPGRADES)
  add_subdirectory(workflow-trk)
  add_subdirectory(workflow-its3)
  add_subdirectory(tools)
  add_subdirectory(study)
endif()
