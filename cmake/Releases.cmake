# Releases.cmake: the in-house components of hydra and the releases it requires, selected by cmake/SemVerDeps.cmake
#
#   include(Releases)
#   hydra_releases(<NAME>...)       select and check out each component (no build)
#
# A requirement is a caret range on the vX.Y.Z tags of the component's repository: 1.1 accepts 1.1.0 and every later
# 1.x, not 2.0.0; 0.2 accepts 0.2.x only (a 0.y release may break 0.y-1). The first configure of a clone takes the
# newest release in the range on the remote and writes it to deps.lock (in the source root, local: git-ignored); the
# next ones take the release deps.lock holds while it satisfies the range. Each component
# is checked out in its folder below (git-ignored), out of the build tree, at a detached tag.
#   -DSEMVER_UPDATE=ALL (or =<NAME>)            select the newest releases in the ranges again
#   -DFETCHCONTENT_SOURCE_DIR_<NAME>=<folder>   build a checkout as it is (a branch, local changes), e.g. src/lib/MOSE
# FiNeR is not here: it is the submodule lib/third_party/FiNeR, pinned at one release.

# <NAME>  <range>  <repository>  <folder>  [NOSUB: none of its submodules]
set(HYDRA_RELEASES
  "OSLO   1.1.1  https://github.com/MarcoGrossi92/OSLO.git   lib/OSLO        NOSUB"
  "ORION  2.2    https://github.com/MarcoGrossi92/ORION.git  lib/ORION"
  "FLINT  2.1    https://github.com/MarcoGrossi92/FLINT.git  lib/FLINT       NOSUB"
  "MOSE   0.2.1  https://github.com/open-hydra/MOSE.git      src/lib/MOSE    NOSUB"
  "FUSS   1.0    https://github.com/open-hydra/FUSS.git      src/lib/FUSS    NOSUB"
  "ARES   1.1    https://github.com/open-hydra/ARES.git      src/lib/ARES    NOSUB"
  "ICE    1.1    https://github.com/open-hydra/ICE.git       src/lib/ICE     NOSUB"
  "IGLOO  1.1.1  https://github.com/open-hydra/IGLOO.git     src/lib/IGLOO   NOSUB"
  "GROOT  1.0    https://github.com/open-hydra/GROOT.git     src/lib/GROOT   NOSUB"
  "ATLAS  0.2    https://github.com/open-hydra/ATLAS.git     utils/ATLAS     NOSUB")
# The components of each application, besides the libraries OSLO, ORION and FLINT
set(HYDRA_COMPONENTS_MF  FUSS MOSE)
set(HYDRA_COMPONENTS_AF  FUSS ARES)
set(HYDRA_COMPONENTS_MI2 MOSE ICE IGLOO)

get_filename_component(HYDRA_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
# One lock file for every configure of the tree: an application's build and cmake/deps alike
set(SEMVER_LOCK_FILE "${HYDRA_ROOT}/deps.lock" CACHE FILEPATH "Releases of the dependencies selected")
include(${CMAKE_CURRENT_LIST_DIR}/SemVerDeps.cmake)

function(hydra_releases)
  foreach(_name IN LISTS ARGN)
    set(_found FALSE)
    foreach(_entry IN LISTS HYDRA_RELEASES)
      string(REGEX REPLACE "[ \t]+" ";" _f "${_entry}")
      list(GET _f 0 _n)
      if(NOT _n STREQUAL _name)
        continue()
      endif()
      set(_found TRUE)
      list(GET _f 1 _range)
      list(GET _f 2 _url)
      list(GET _f 3 _dir)
      set(_sub "")
      if("NOSUB" IN_LIST _f)
        set(_sub GIT_SUBMODULES "")
      endif()
      # A checkout given by a relative path is taken from the source root
      if(FETCHCONTENT_SOURCE_DIR_${_name})
        get_filename_component(_local "${FETCHCONTENT_SOURCE_DIR_${_name}}" ABSOLUTE BASE_DIR "${HYDRA_ROOT}")
        set(FETCHCONTENT_SOURCE_DIR_${_name} "${_local}" CACHE PATH "Local checkout of ${_name}" FORCE)
      endif()
      semver_dependency(${_name} GIT_REPOSITORY ${_url} VERSION ${_range} ${_sub} SOURCE_DIR "${HYDRA_ROOT}/${_dir}")
      # Where the component is: the checkout given, or its folder
      if(FETCHCONTENT_SOURCE_DIR_${_name})
        set(HYDRA_${_name}_DIR "${FETCHCONTENT_SOURCE_DIR_${_name}}" PARENT_SCOPE)
      else()
        set(HYDRA_${_name}_DIR "${HYDRA_ROOT}/${_dir}" PARENT_SCOPE)
      endif()
    endforeach()
    if(NOT _found)
      message(FATAL_ERROR "hydra_releases: no component ${_name} in cmake/Releases.cmake")
    endif()
  endforeach()
endfunction()
