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

# ifx 2024.1 compiles the ghost metrics of MOSE and ARES wrongly at -O2 and above without IPO: the check MOSE, ARES and
# FUSS make when built alone, which does not run inside hydra. It only reads the flags of the build type.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 2024.1
   AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 2024.2)
  string(TOUPPER "${CMAKE_BUILD_TYPE}" _ifx_bt)
  separate_arguments(_ifx_flags UNIX_COMMAND "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${_ifx_bt}}")
  set(_ifx_level "no -O option (-O2)")
  set(_ifx_high TRUE)
  set(_ifx_ipo FALSE)
  if(CMAKE_INTERPROCEDURAL_OPTIMIZATION OR CMAKE_INTERPROCEDURAL_OPTIMIZATION_${_ifx_bt})
    set(_ifx_ipo TRUE)
  endif()
  foreach(_ifx_f IN LISTS _ifx_flags)
    if(_ifx_f MATCHES "^[-/]O[01]$")
      set(_ifx_level "${_ifx_f}")
      set(_ifx_high FALSE)
    elseif(_ifx_f MATCHES "^[-/]O([23s]|fast)?$")
      set(_ifx_level "${_ifx_f}")
      set(_ifx_high TRUE)
    elseif(_ifx_f MATCHES "^[-/]Q?ipo$")
      set(_ifx_ipo TRUE)
    endif()
  endforeach()
  if(_ifx_high AND NOT _ifx_ipo)
    message(FATAL_ERROR "ifx ${CMAKE_Fortran_COMPILER_VERSION} with ${_ifx_level} and without IPO compiles a loop of the "
                        "ghost metrics of MOSE and ARES wrongly (an error of the compiler: a minimal program in standard "
                        "Fortran gives wrong values at -O2 and -O3 and the right ones at -O0, -O1 and with -ipo); a solver "
                        "built so stops at set-up on the check of its ghost metrics or on a negative volume, or runs with "
                        "wrong ghost metrics. Build with -O0 (figaro build hydra-<SUFFIX> --type=DEBUG), with gfortran "
                        "(--compilers=gnu), or with IPO (-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON).")
  endif()
endif()
# Locate and set parallelization libraries. There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

# Set directories
include(GNUInstallDirs)
set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

set(MASTER "hydra")

# ORION, OSLO and FLINT by release range (cmake/SemVerDeps.cmake), checked out in lib/ORION, lib/OSLO and lib/FLINT.
# hydra is the top of its build tree: these requirements select the releases, before any solver or library checks
# its own. deps.lock (tracked) holds the releases a build takes; -DSEMVER_UPDATE=ALL selects the newest ones in the
# ranges again, -DFETCHCONTENT_SOURCE_DIR_<NAME>=<path> builds a local checkout instead. FLINT sets FLINT_SRC, the
# sources MOSE and ARES compile, in this directory
include(${CMAKE_SOURCE_DIR}/cmake/SemVerDeps.cmake)
# Without the submodules lib/sundials of OSLO and lib/cantera of FLINT, which they clone only when they build them
semver_dependency(OSLO GIT_REPOSITORY https://github.com/MarcoGrossi92/OSLO.git VERSION 1.1 GIT_SUBMODULES ""
                  SOURCE_DIR lib/OSLO)
semver_dependency(ORION GIT_REPOSITORY https://github.com/MarcoGrossi92/ORION.git VERSION 2.1 SOURCE_DIR lib/ORION)
semver_dependency(FLINT GIT_REPOSITORY https://github.com/MarcoGrossi92/FLINT.git VERSION 2.1 GIT_SUBMODULES ""
                  SOURCE_DIR lib/FLINT)
add_subdirectory(lib/third_party/FiNeR/)
FetchContent_MakeAvailable(OSLO ORION FLINT)

IF(EXE_NAME STREQUAL "MF")
    add_subdirectory(src/lib/FUSS/)
    add_subdirectory(src/lib/MOSE/)
    add_subdirectory(src/app)
ELSEIF(EXE_NAME STREQUAL "AF")
    add_subdirectory(src/lib/FUSS/)
    add_subdirectory(src/lib/ARES/)
    add_subdirectory(src/app)
ELSEIF(EXE_NAME STREQUAL "MI2")
    # MOSE (Eulerian gas) + ICE (Eulerian condensed) + IGLOO (Lagrangian).
    #
    # MPI is refused at CONFIGURE time, deliberately and permanently -- not deferred work.
    # MOSE decomposes over blocks (is_local_block) while IGLOO requires the FULL gas field
    # replicated identically on every rank (see src/lib/IGLOO/src/lib/obj_IGLOO.f90 setup_static),
    # so a block-partitioned MOSExIGLOO would silently integrate particles through holes.
    # ICE partitions its own blocks under USE_MPI, but nothing makes its decomposition
    # match MOSE's, and IGLOO's requirement above stands regardless. Refusing to configure
    # beats linking a wrong binary.
    IF(USE_MPI)
        message(FATAL_ERROR
            "hydra-MI2 does not support MPI. MOSE partitions over blocks but IGLOO requires the "
            "full gas field on every rank. "
            "Rebuild without --use-mpi.")
    ENDIF()

    add_subdirectory(src/lib/MOSE/)
    add_subdirectory(src/lib/ICE/)
    add_subdirectory(src/lib/IGLOO/)
    add_subdirectory(src/app)
ELSE()
    message(FATAL_ERROR "EXE_NAME not set correctly. Choose from: MF, AF, MI2")
ENDIF()
