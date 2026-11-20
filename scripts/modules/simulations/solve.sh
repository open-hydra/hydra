#!/bin/bash

# solve.sh - logic for running supported solvers (MOSE, FUSS, GROOT, ARES, hydra)

source "$SCRIPT_DIR/utils/common.sh"

# Argument parsing with getopt
TEMP=$(getopt -o bnp:m:c: --long background,nosleep,parallel:,mpi:,cpus:, -- "$@") || exit 1
eval set -- "$TEMP"

# Default flags
BG=0
NTHREADS=1
NRANKS=0
CPUS=

while true; do
  case "$1" in
    -b|--background) BG=1; shift;;
    -n|--nosleep) NOSLEEP=1; shift;;
    -p|--parallel) NTHREADS=$2; shift 2;;
    -m|--mpi) NRANKS=$2; shift 2;;
    -c|--cpus) CPUS=$2; shift 2;;
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
      # FORT_BUFFERED is left to the caller: =true makes an ifx build write its files a buffer at a time, and a SIGKILL,
      # SIGBUS or SIGHUP then loses what the buffer holds (an end, an error stop, SIGTERM/INT/SEGV/FPE/ABRT flush it)
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
  local launch=("$solver_exec")
  # MPI: launch through mpirun; .ID then holds mpirun's PID, which forwards a kill to every rank
  if [[ ${NRANKS:-0} -gt 0 ]]; then
    if ! command -v mpirun > /dev/null; then
      log_error "mpirun not found: load the MPI environment first (e.g. source /opt/intel/oneapi/mpi/<version>/env/vars.sh)"
      exit 1
    fi
    case ${I_MPI_PIN:-on} in
      off|0|no|disable|false)
        # Pinning off: OMP_PLACES=cores alone would still bind the threads of every rank to the
        # same first cores of the machine (sprop2, 2026-10-08), so no thread binding either
        export OMP_PROC_BIND=${OMP_PROC_BIND:-false} ;;
      *)
        # Give each rank its own OMP_NUM_THREADS physical cores (Intel MPI; ignored by other MPI
        # libraries): without PIN_CELL=core the domain counts hyperthreads, i.e. half the cores
        export I_MPI_PIN_DOMAIN=${I_MPI_PIN_DOMAIN:-omp}
        export I_MPI_PIN_CELL=${I_MPI_PIN_CELL:-core}
        # and one thread per core inside it, instead of letting threads share a core
        export OMP_PLACES=${OMP_PLACES:-cores}
        export OMP_PROC_BIND=${OMP_PROC_BIND:-close} ;;
    esac
    launch=(mpirun -np "$NRANKS" "$solver_exec")
  fi
  # --cpus: run inside that set of logical CPUs. Intel MPI keeps its pinning domains inside the set
  # (checked 2026-10-08), so two cases on one machine get disjoint cores; without it both start from core 0
  if [[ -n ${CPUS:-} ]]; then
    if ! taskset -c "$CPUS" true 2> /dev/null; then
      log_error "--cpus $CPUS: not a valid CPU list for taskset on $(hostname)"
      exit 1
    fi
    launch=(taskset -c "$CPUS" "${launch[@]}")
  fi
  mkdir -p OUTPUT
  if [[ ${BG:-0} -eq 0 ]]; then
    "${launch[@]}"
  else
    "${launch[@]}" 2>errors_file >logfile &
    echo $! > .ID
  fi
}

run_solver "$@"
