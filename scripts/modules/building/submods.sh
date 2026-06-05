#!/bin/bash

# building/submods.sh - Builds a selected PROG module

source "$SCRIPT_DIR/utils/common.sh"

# Preliminaries
PROG="$1"
if [[ $PROG == "ATLAS" || $PROG == "GRIB" ]]; then
  where=utils
elif [[ $PROG == "MOSE" || $PROG == "FUSS" || $PROG == "GROOT" || $PROG == "PESTE" ]]; then
  where=src/lib
else
  where=lib
fi

cd $HYDRADIR/$where/$PROG || log_error "PROG $PROG not found."

# Execute the install script with 'build' argument
log_info "Building $PROG"
log_info "Passing jurisdiction to install.sh in $HYDRADIR/$where/$PROG"
echo
./install.sh build --include-orion=$HYDRADIR/lib/ORION \
                   --include-finer=$HYDRADIR/lib/third_party/FiNeR \
                   --include-oslo=$HYDRADIR/lib/OSLO \
                   --include-flint=$HYDRADIR/lib/FLINT \
                "${@:2}"
echo
log_success "$PROG build completed"
