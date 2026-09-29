#!/bin/bash

# solve.sh - logic for running supported solvers (MOSE, FUSS, GROOT, ARES, hydra)

source "$SCRIPT_DIR/utils/common.sh"

# Argument parsing with getopt
TEMP=$(getopt -o bnp: --long background,nosleep,parallel:, -- "$@") || exit 1
eval set -- "$TEMP"

# Default flags
BG=0
NTHREADS=1

while true; do
  case "$1" in
    -b|--background) BG=1; shift;;
    -n|--nosleep) NOSLEEP=1; shift;;
    -p|--parallel) NTHREADS=$2; shift 2;;
    --) shift; break;;
    *) echo "Unknown option: $1"; exit 1;;
  esac
done

run_solver() {
  local solver=$1
  local action=$2
  shift 2

  require_command ulimit
  ulimit -s unlimited

  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=${NTHREADS:-1}

  local solver_exec
  case "$solver" in
    MOSE|FUSS|GROOT|ARES)
      local_var="${solver}LOCAL"
      master_var="${solver}MASTER"

      local="${!local_var}"
      master="${!master_var}"

      [[ "$master" -nt "$local" ]] && cp "$master" "$local"
      [[ "$action" == "solve" ]] && run_generic_solver "$local"
      ;;

    hydra*)
      ls INPUT/*phase.txt > filelist.txt 2>/dev/null
      # CMake writes the master as bin/hydra_<SUFFIX> (underscore); the case-local copy keeps the
      # figaro name bin/hydra-<SUFFIX> (dash). The lookup used the dash name for both and never found
      # the master, so `figaro hydra-MI2 solve` could not run (fixed 2026-09-16).
      local="./bin/$solver"
      master="$HYDRADIR/bin/${solver/-/_}"
      if [[ ! -f "$master" && ! -f "$local" ]]; then
        log_error "No executable for $solver: build it first (figaro build $solver ...) -- expected $master"
        exit 1
      fi
      [[ "$master" -nt "$local" ]] && mkdir -p ./bin && cp "$master" "$local"
      [[ "$action" == "solve" ]] && run_generic_solver "$local"
      sleep 3
      rm -f filelist.txt
      ;;

    *)
      log_error "Unknown solver: $solver"
      exit 1
      ;;
  esac
}

run_generic_solver() {
  local solver_exec=$1
  mkdir -p OUTPUT
  if [[ ${BG:-0} -eq 0 ]]; then
    "$solver_exec"
  else
    "$solver_exec" 2>errors_file >logfile &
    echo $! > .ID
  fi
}

run_solver "$@"
