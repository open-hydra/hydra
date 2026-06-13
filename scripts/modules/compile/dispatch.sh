#!/usr/bin/env bash

# compile/dispatch.sh - Dispatches compilation tasks to appropriate handlers

source "$SCRIPT_DIR/utils/common.sh"

# Argument parsing with getopt
TEMP=$(getopt -o e --long edit, -- "$@") || exit 1
eval set -- "$TEMP"

# Default flags
edit=false

while [ : ]; do
  case "$1" in
    -e|--edit) edit=true; shift;;
    --) shift; break;;
    *) echo "Unknown option: $1"; exit 1;;
  esac
done

case "$@" in
    hydra*)
        source "$SCRIPT_DIR/modules/compile/hydra.sh" "$@" ;;
    ATLAS|GRIB|MOSE|FUSS|GROOT|PESTE|ARES)
        source "$SCRIPT_DIR/modules/compile/submods.sh" "$@" ;;
    *)
        log_error "Unknown compile target: $@"
        exit 1
        ;;
esac