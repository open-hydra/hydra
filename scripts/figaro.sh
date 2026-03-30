#!/bin/bash
# figaro.sh - Entry point for all tasks
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/utils/config.sh"
source "$SCRIPT_DIR/utils/common.sh"
source "$SCRIPT_DIR/utils/print_usage.sh"

[[ $# -lt 1 ]] && print_usage
COMMAND=$1
shift

case $COMMAND in
  # Sanity checks
  doctor)         source "$SCRIPT_DIR/doctor.sh" ;;

  # Compiling management tasks
  compile)        source "$SCRIPT_DIR/modules/compile/dispatch.sh" "$@" ;;

  # Building management tasks
  build)          source "$SCRIPT_DIR/modules/building/dispatch.sh" "$@" ;;

  # git submodules tasks
  submodule)         source "$SCRIPT_DIR/modules/git.sh" "$@";;

  # Simulations tasks
  clean)          source "$SCRIPT_DIR/modules/simulations/clean.sh" ;;
  ATLAS)          source "$SCRIPT_DIR/modules/simulations/ATLAS.sh" "$@" ;;
  kill)           source "$SCRIPT_DIR/modules/simulations/kill.sh" ;;
  hydra*|MOSE|IGLOO|ICE|FUSS|MOSKA|Q2D)
                  source "$SCRIPT_DIR/modules/simulations/solve.sh" "$COMMAND" "$@" ;;

  *)
    log_error "Unknown command: $COMMAND"
    print_usage
    exit 1
    ;;
esac
