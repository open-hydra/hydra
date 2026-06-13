# CMake project file for 'hydra'

cmake_minimum_required(VERSION 3.23)

project (hydra LANGUAGES Fortran CXX C)

set(CMAKE_VERBOSE_MAKEFILE FALSE)
# search path for cmake modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
#OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
OPTION(USE_OPENMP "Use OpenMP for parallelization" ON)

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING. You should  review this file and make sure the flags 
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
# Locate and set parallelization libraries. There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

# Set directories
include(GNUInstallDirs)
set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

set(MASTER "hydra")
set ( ORION_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/ORION/")
set ( FINER_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/third_party/FiNeR/")
set ( OSLO_PATH  "${CMAKE_CURRENT_SOURCE_DIR}/lib/OSlo/")
set ( FLINT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/FLINT/")
IF(EXE_NAME STREQUAL "MF")
    add_subdirectory(lib/ORION)
    add_subdirectory(lib/third_party/FiNeR/)
    add_subdirectory(src/lib/FUSS/)
    add_subdirectory(lib/FLINT/)
    add_subdirectory(src/lib/MOSE/)
    add_subdirectory(src/app)
ELSEIF(EXE_NAME STREQUAL "AF")
    add_subdirectory(lib/ORION)
    add_subdirectory(lib/third_party/FiNeR/)
    add_subdirectory(src/lib/FUSS/)
    add_subdirectory(lib/FLINT/)
    add_subdirectory(src/lib/ARES/)
    add_subdirectory(src/app)
ELSE()
    message(FATAL_ERROR "EXE_NAME not set correctly. Choose from: MF, AF")
ENDIF()
