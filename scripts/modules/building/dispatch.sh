#!/usr/bin/env bash

# building/dispatch.sh - Dispatches building tasks to appropriate handlers

case "$1" in
    hydra*)
        source "$SCRIPT_DIR/modules/building/hydra.sh" "$@" ;;
    ATLAS|GRIB|MOSE|FUSS|GROOT)
        source "$SCRIPT_DIR/modules/building/submods.sh" "$@" ;;
    *)
        log_error "Unknown build target: $@"
        exit 1
        ;;
esac