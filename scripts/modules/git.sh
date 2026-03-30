#!/bin/bash

# git.sh - management of git submodules

source "$SCRIPT_DIR/utils/common.sh"
cd $HYDRADIR

# Argument parsing with getopt
TEMP=$(getopt -o rf --long remote,force -- "$@") || exit 1
eval set -- "$TEMP"

# Default flags
remote=false
force=false
while [ : ]; do
  case "$1" in
    -r|--remote) remote=true; shift;;
    -f|--force) force=true; shift;;
    --) shift; break;;
    *) echo "Unknown option: $1"; exit 1;;
  esac
done


if [ $# -eq 1 ]; then
    arg=""
else
    if [[ "$2" == "ATLAS" || "$2" == "GRIB" ]]; then
        arg="utils/$2"
    else
        arg="lib/$2"
    fi
fi

case "$1" in
    update)
        if [ "$remote" = true ]; then
            git submodule update --init --remote "$arg"
        else
            git submodule update --init "$arg"
        fi
        ;;

    delete)
        if [ "$force" = true ]; then
            git submodule deinit -f -all "$arg"
        else
            git submodule deinit -all "$arg"
        fi
        ;;

    *)
        log_error "Unknown command for git submodules: $@"
        exit 1
        ;;
esac