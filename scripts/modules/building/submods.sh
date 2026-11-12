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
if [[ $PROG == "ATLAS" ]]; then
  # ATLAS selects its own ORION release (2.1 or later, cloned at its first configure) and refuses an older checkout,
  # which hydra's lib/ORION still is
  ./install.sh build "${@:2}"
else
  ./install.sh build --include-orion=$HYDRADIR/lib/ORION \
                    --include-finer=$HYDRADIR/lib/third_party/FiNeR \
                    --include-oslo=$HYDRADIR/lib/OSLO \
                    --include-flint=$HYDRADIR/lib/FLINT \
                  "${@:2}"
fi
echo
log_success "$PROG build completed"
