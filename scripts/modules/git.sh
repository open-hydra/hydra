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
    elif [[ "$2" == "MOSE" || "$2" == "FUSS" || "$2" == "GROOT" || "$2" == "ARES" || "$2" == "ICE" || "$2" == "IGLOO" ]]; then
        arg="src/lib/$2"
    else
        arg="lib/$2"
    fi
fi

# ICE is skipped by a plain update (update = none in .gitmodules): asked by name, it is checked out.
# Without a name the update covers every submodule.
target=()
[[ -n "$arg" ]] && target=("$arg")
[[ "$arg" == "src/lib/ICE" ]] && target=(--checkout "$arg")

case "$1" in
    update)
        if [ "$remote" = true ]; then
            git submodule update --init --remote ${target[@]+"${target[@]}"}
        else
            git submodule update --init ${target[@]+"${target[@]}"}
        fi
        ;;

    delete)
        if [ "$force" = true ]; then
            if [ "$arg" = "" ]; then
                git submodule deinit -f --all
            else
                git submodule deinit -f "$arg"
            fi
        else
            if [ "$arg" = "" ]; then
                git submodule deinit --all
            else
                git submodule deinit "$arg"
            fi
        fi
        ;;

    *)
        log_error "Unknown command for git submodules: $@"
        exit 1
        ;;
esac