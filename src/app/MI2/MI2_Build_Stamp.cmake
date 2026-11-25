# Writes OUT, the Fortran module hydraMI2_Build_Stamp: the commit of each repository hydra-MI2 is built from, after
# its release tag when it is at one. Script mode: cmake -DSRC=<hydra root> -DOUT=<file.f90> [-DDIR_<NAME>=<dir>...]
# -P MI2_Build_Stamp.cmake (the folders the build used, for MOSE ICE IGLOO ORION OSLO FLINT; default the folders of
# cmake/Releases.cmake)
find_program(GIT_EXE git)
set(repos "hydra=." "MOSE=src/lib/MOSE" "ICE=src/lib/ICE" "IGLOO=src/lib/IGLOO"
          "ORION=lib/ORION" "OSLO=lib/OSLO" "FLINT=lib/FLINT" "FiNeR=lib/third_party/FiNeR")
set(lines "")
foreach(entry ${repos})
  string(REPLACE "=" ";" pair "${entry}")
  list(GET pair 0 name)
  list(GET pair 1 rel)
  set(stamp "unknown")
  set(dir "${SRC}/${rel}")
  if(DIR_${name})
    set(dir "${DIR_${name}}")
  endif()
  if(GIT_EXE AND EXISTS "${dir}/.git")
    execute_process(COMMAND ${GIT_EXE} -C "${dir}" rev-parse --short=7 HEAD
                    OUTPUT_VARIABLE sha RESULT_VARIABLE rc OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(rc EQUAL 0)
      set(stamp "${sha}")
      # a component at a release: its tag first
      if(NOT name STREQUAL "hydra")
        execute_process(COMMAND ${GIT_EXE} -C "${dir}" describe --tags --exact-match HEAD
                        OUTPUT_VARIABLE tag RESULT_VARIABLE rc OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
        if(rc EQUAL 0)
          set(stamp "${tag} ${sha}")
        endif()
      endif()
      # modified tracked files only; a submodule's own state is its own line
      execute_process(COMMAND ${GIT_EXE} --no-optional-locks -C "${dir}" status --porcelain
                              --untracked-files=no --ignore-submodules=all
                      OUTPUT_VARIABLE dirty RESULT_VARIABLE rc OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
      if(rc EQUAL 0 AND NOT dirty STREQUAL "")
        string(APPEND stamp "-dirty")
      endif()
    endif()
  endif()
  string(APPEND lines "    write(*,'(A)') '# [MI2] built from ${name} ${stamp}'\n")
endforeach()
get_filename_component(dir "${OUT}" DIRECTORY)
file(MAKE_DIRECTORY "${dir}")
file(WRITE "${OUT}.tmp"
  "module hydraMI2_Build_Stamp\n  implicit none\n  private\n  public :: Print_Build_Stamp\ncontains\n"
  "  !> One line per repository: its commit, and -dirty when tracked files are modified.\n"
  "  subroutine Print_Build_Stamp()\n${lines}  end subroutine Print_Build_Stamp\nend module hydraMI2_Build_Stamp\n")
file(COPY_FILE "${OUT}.tmp" "${OUT}" ONLY_IF_DIFFERENT)
file(REMOVE "${OUT}.tmp")
