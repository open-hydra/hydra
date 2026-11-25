# Checks out components of hydra at their releases without building anything: what figaro runs before it builds a
# component on its own (figaro build ATLAS, GROOT, MOSE, ...).
#   cmake -S cmake/deps -B build/.deps -DHYDRA_DEPS="ATLAS;ORION" [-DSEMVER_UPDATE=ALL]
# The releases come from deps.lock of the source root, else the newest in the ranges of cmake/Releases.cmake, and
# deps.lock is written as in the build of an application.
cmake_minimum_required(VERSION 3.23)
project(hydra-deps NONE)

if(NOT HYDRA_DEPS)
  message(FATAL_ERROR "hydra-deps: give the components to check out, -DHYDRA_DEPS=<NAME>;<NAME>...")
endif()
include(${CMAKE_CURRENT_SOURCE_DIR}/../Releases.cmake)
hydra_releases(${HYDRA_DEPS})
foreach(_name IN LISTS HYDRA_DEPS)
  message(STATUS "hydra-deps: ${_name} in ${HYDRA_${_name}_DIR}")
endforeach()
