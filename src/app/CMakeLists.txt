function(add_fortran_executable_from_dirs TARGET_NAME)
    # Usage:
    #   add_fortran_executable_from_dirs(name SRCDIRS dir1 dir2 ... LIBS lib1 lib2 ... [EXCLUDE file1 file2 ...])
    set(options)
    set(oneValueArgs)
    set(multiValueArgs SRCDIRS LIBS EXCLUDE)
    cmake_parse_arguments(AFE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT AFE_SRCDIRS)
        message(FATAL_ERROR "add_fortran_executable_from_dirs: no SRCDIRS specified")
    endif()

    # Collect sources from all provided directories
    set(ALL_SRC_FILES "")
    foreach(DIR ${AFE_SRCDIRS})
        file(GLOB DIR_SOURCES CONFIGURE_DEPENDS "${DIR}/*.f90")
        list(APPEND ALL_SRC_FILES ${DIR_SOURCES})
        # A few helpers have to be C: the cycle counter calls perf_event_open,
        # which has no Fortran binding.
        file(GLOB DIR_C_SOURCES CONFIGURE_DEPENDS "${DIR}/*.c")
        list(APPEND ALL_SRC_FILES ${DIR_C_SOURCES})
    endforeach()

    # Leave out the files named in EXCLUDE (file names, matched at the end of the path)
    foreach(NAME ${AFE_EXCLUDE})
        list(FILTER ALL_SRC_FILES EXCLUDE REGEX "/${NAME}$")
    endforeach()

    if(NOT ALL_SRC_FILES)
        message(FATAL_ERROR "No Fortran source files found in directories: ${AFE_SRCDIRS}")
    endif()

    # Create the executable
    add_executable(${TARGET_NAME} ${ALL_SRC_FILES})

    # All sources here are Fortran, but TecIO is a C++ archive: without this
    # CMake infers CXX for the link step and the Fortran main is never seen
    # ("undefined reference to `main'").
    set_target_properties(${TARGET_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

    # Link libraries
    target_link_libraries(${TARGET_NAME} PRIVATE ${AFE_LIBS})

    # Fortran module directory handling
    target_include_directories(${TARGET_NAME} PRIVATE
        $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    )

    # Optional: OpenMP support
    if(USE_OPENMP)
        target_link_libraries(${TARGET_NAME} PRIVATE OpenMP::OpenMP_Fortran)
    endif()

    # Optional: TECIO / pthread support
    if(USE_TECIO)
        target_link_libraries(${TARGET_NAME} PRIVATE pthread)
    endif()
endfunction()


# -------------------------------------------------------------------------------------

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# hydra-MF
if(EXE_NAME STREQUAL "MF")
    add_fortran_executable_from_dirs(
        hydra_MF
        SRCDIRS
            ${CMAKE_SOURCE_DIR}/src/app/common/base
            ${CMAKE_SOURCE_DIR}/src/app/common/config
            ${CMAKE_SOURCE_DIR}/src/app/MF
        LIBS FUSSL MOSEL FiNeR)
ENDIF()

# hydra-AF
if(EXE_NAME STREQUAL "AF")
    add_fortran_executable_from_dirs(
        hydra_AF
        SRCDIRS
            ${CMAKE_SOURCE_DIR}/src/app/common/base
            ${CMAKE_SOURCE_DIR}/src/app/common/config
            ${CMAKE_SOURCE_DIR}/src/app/AF
        LIBS FUSSL ARESL FiNeR)
ENDIF()

# hydra-MI2
# ORION is listed explicitly: ICEL links it PRIVATE, so type(orion_data) -- which the MI2 coupler
# declares -- does not arrive transitively.
# The shared timers are left out: they use FUSS's MPI module, and only MF and AF call them.
if(EXE_NAME STREQUAL "MI2")
    add_fortran_executable_from_dirs(
        hydra_MI2
        SRCDIRS
            ${CMAKE_SOURCE_DIR}/src/app/common/base
            ${CMAKE_SOURCE_DIR}/src/app/common/config
            ${CMAKE_SOURCE_DIR}/src/app/MI2
        LIBS MOSEL ICEL IGLOOL ORION FiNeR
        EXCLUDE Mod_Timers.f90 hydra_cycles.c)

    # The commit of each repository the executable is built from, written again at every build
    set(MI2_STAMP ${CMAKE_BINARY_DIR}/generated/MI2_Build_Stamp.f90)
    set(MI2_STAMP_CMD ${CMAKE_COMMAND} -DSRC=${CMAKE_SOURCE_DIR} -DOUT=${MI2_STAMP}
                      -P ${CMAKE_SOURCE_DIR}/src/app/MI2/MI2_Build_Stamp.cmake)
    execute_process(COMMAND ${MI2_STAMP_CMD})
    add_custom_target(hydra_MI2_stamp COMMAND ${MI2_STAMP_CMD} BYPRODUCTS ${MI2_STAMP} VERBATIM)
    target_sources(hydra_MI2 PRIVATE ${MI2_STAMP})
    add_dependencies(hydra_MI2 hydra_MI2_stamp)
ENDIF()
