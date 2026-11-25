#!/bin/bash

# building/submods.sh - Builds a selected PROG module

source "$SCRIPT_DIR/utils/common.sh"

# Preliminaries
PROG="$1"
if [[ $PROG == "ATLAS" || $PROG == "GRIB" ]]; then
  where=utils
elif [[ $PROG == "MOSE" || $PROG == "FUSS" || $PROG == "GROOT" || $PROG == "ARES" ]]; then
  where=src/lib
else
  where=lib
fi

# A component of cmake/Releases.cmake: checked out at its release (deps.lock), with ORION, which all of them read
if grep -q "^ *\"$PROG " "$HYDRADIR/cmake/Releases.cmake"; then
  source "$SCRIPT_DIR/modules/deps.sh" ORION "$PROG"
fi
cd $HYDRADIR/$where/$PROG || log_error "PROG $PROG not found."

# The solvers are handed hydra's FiNeR: make sure its dependencies are there. ATLAS builds with its own
# pinned FiNeR, the one its CI tests, so hydra's FiNeR can move ahead of ATLAS's pointer.
if [[ $PROG != "ATLAS" ]]; then
  source "$SCRIPT_DIR/utils/finer_deps.sh"
  fill_finer_deps "$HYDRADIR"
fi

# Execute the install script with 'build' argument
log_info "Building $PROG"
log_info "Passing jurisdiction to install.sh in $HYDRADIR/$where/$PROG"
echo
# hydra's libraries, once its configure has checked their releases out in lib/: the module builds on the same ones;
# before that, it selects its own
include() { [[ -e "$HYDRADIR/$2/.git" ]] && echo "--include-$1=$HYDRADIR/$2"; true; }
if [[ $PROG == "ATLAS" ]]; then
  # ATLAS requires ORION 2.1 or later, which deps.lock selects; it builds with its own pinned FiNeR
  ./install.sh build $(include orion lib/ORION) "${@:2}"
elif [[ $PROG == "MOSE" || $PROG == "ARES" ]]; then
  ./install.sh build $(include orion lib/ORION) \
                    $(include finer lib/third_party/FiNeR) \
                    $(include oslo lib/OSLO) \
                    $(include flint lib/FLINT) \
                  "${@:2}"
else
  # FUSS and GROOT read neither OSLO nor FLINT, and their install.sh refuse those options
  ./install.sh build $(include orion lib/ORION) \
                    $(include finer lib/third_party/FiNeR) \
                  "${@:2}"
fi
echo
log_success "$PROG build completed"
