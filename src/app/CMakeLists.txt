function(add_fortran_executable_from_dirs TARGET_NAME)
    # Usage:
    #   add_fortran_executable_from_dirs(name SRCDIRS dir1 dir2 ... LIBS lib1 lib2 ...)
    set(options)
    set(oneValueArgs)
    set(multiValueArgs SRCDIRS LIBS)
    cmake_parse_arguments(AFE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT AFE_SRCDIRS)
        message(FATAL_ERROR "add_fortran_executable_from_dirs: no SRCDIRS specified")
    endif()

    # Collect sources from all provided directories
    set(ALL_SRC_FILES "")
    foreach(DIR ${AFE_SRCDIRS})
        file(GLOB DIR_SOURCES CONFIGURE_DEPENDS "${DIR}/*.f90")
        list(APPEND ALL_SRC_FILES ${DIR_SOURCES})
    endforeach()

    if(NOT ALL_SRC_FILES)
        message(FATAL_ERROR "No Fortran source files found in directories: ${AFE_SRCDIRS}")
    endif()

    # Create the executable
    add_executable(${TARGET_NAME} ${ALL_SRC_FILES})

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
            ${CMAKE_SOURCE_DIR}/src/common
            ${CMAKE_SOURCE_DIR}/src/MF
        LIBS FUSSL MOSEL FiNeR)
ENDIF()

# hydra-MI2
if(EXE_NAME STREQUAL "MI2")
    add_fortran_executable_from_dirs(
        hydra_MI2
        SRCDIRS
            ${CMAKE_SOURCE_DIR}/src/common
            ${CMAKE_SOURCE_DIR}/src/MI2
        LIBS MOSEL ICEL IGLOOL FiNeR)
ENDIF()

# hydra-MQ2D
if(EXE_NAME STREQUAL "MQ2D")
    add_fortran_executable_from_dirs(
        hydra_MQ2D
        SRCDIRS
            ${CMAKE_SOURCE_DIR}/src/common
            ${CMAKE_SOURCE_DIR}/src/MQ2D
        LIBS MOSKAL Q2DL FiNeR)
ENDIF()

# hydra-MM
if(EXE_NAME STREQUAL "MM")
    add_fortran_executable_from_dirs(
        hydra_MM
        SRCDIRS
            ${CMAKE_SOURCE_DIR}/src/common
            ${CMAKE_SOURCE_DIR}/src/MM
        LIBS MOSKAL MOSEL FiNeR)
ENDIF()