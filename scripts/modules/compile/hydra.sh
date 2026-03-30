#!/usr/bin/env bash

# compile/hydra.sh - Compilation handler for hydra

# Preliminaries
suffix="${1#hydra-}"
BUILD_DIR=$HYDRADIR/build/${suffix}
if [[ ! -d "$BUILD_DIR" ]]; then
    log_error "Build directory not found. Run 'figaro build hydra' first."
    exit 1
fi
log_info "Compiling hydra-${suffix}"

# Edit CMakePresets.json if requested
PRESET_PATH="$HYDRADIR/CMakePresets.json"
if [[ ! -f "$PRESET_PATH" ]]; then
  log_error "CMakePresets.json not found in hydra-${BUILD_DIR}. Run 'figaro build hydra-${suffix}' first."
  exit 1
fi
if [[ $edit == true ]]; then
  log_info "Editing CMakePresets.json before compilation"
  ${EDITOR:-vim} "$PRESET_PATH"
fi

# Compile hydra*
log_info "Compiling hydra-${suffix} using CMakePresets"
cmake --preset $suffix || exit 1
cmake --build "$BUILD_DIR" || exit 1
