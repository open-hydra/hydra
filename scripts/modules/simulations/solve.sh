#!/bin/bash

# solve.sh - logic for running supported solvers (MOSE, FUSS, GROOT, PESTE, hydra)

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
    MOSE|FUSS|GROOT|PESTE)
      local_var="${solver}LOCAL"
      master_var="${solver}MASTER"

      local="${!local_var}"
      master="${!master_var}"

      [[ "$master" -nt "$local" ]] && cp "$master" "$local"
      [[ "$action" == "solve" ]] && run_generic_solver "$local"
      ;;

    hydra*)
      ls INPUT/*phase.txt > filelist.txt 2>/dev/null
      local="./bin/$solver"
      [[ "$HYDRADIR/bin/$solver" -nt "$local" ]] && cp "$HYDRADIR/bin/$solver" "$local"
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
