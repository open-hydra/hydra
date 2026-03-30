#!/bin/bash

# compile/submods.sh - Compiles a selected PROG module

# Preliminaries
PROG="$1"

if [[ $PROG == "ATLAS" || $PROG == "GRIB" ]]; then
  where=utils
else
  where=lib
fi

cd "$HYDRADIR/$where/$PROG" || log_error "PROG $PROG not found"
log_info "Compiling $PROG"

# Edit CMakePresets.json if requested
PRESET_PATH="$HYDRADIR/$where/$PROG/CMakePresets.json"
if [[ ! -f "$PRESET_PATH" ]]; then
  log_error "CMakePresets.json not found in $HYDRADIR/$where/$PROG"
  exit 1
fi
if [[ $edit == true ]]; then
  log_info "Editing CMakePresets.json before compilation"
  ${EDITOR:-vim} "$PRESET_PATH"
fi

# Compile the PROG
log_info "Passing jurisdiction to install.sh in $HYDRADIR/$where/$PROG"
echo
./install.sh compile
echo
log_success "$PROG compilation completed"
