#!/bin/bash
# run.sh [-x exe] [-s scratch] [-t threads] [case ...]
# Runs the hydra-MI2 test cases of cases/ and judges each one. A case with a check.sh is judged by
# `check.sh <run-dir> <exe>`, which runs its own decks; any other is run in <scratch>/<case> by run_case.sh (the log is
# run.log) and judged by `check.py <run-dir>`. Without cases: all, in the order below (a channel case compares against
# its base, which runs before it); fpe, which builds a trapping binary, runs only when named. The scratch directory
# must lie outside test/: each run directory is emptied first.
# Exit 1 when a case fails.
set -u
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exe=$here/../../bin/hydra_MI2; scratch=${MI2_SCRATCH:-/tmp/mi2-tests-$USER}; nt=1
while getopts x:s:t: o; do
  case $o in x) exe=$OPTARG;; s) scratch=$OPTARG;; t) nt=$OPTARG;; *) exit 2;; esac
done
shift $((OPTIND - 1))
cases=("$@")
[ ${#cases[@]} -gt 0 ] || cases=(
  igloo-evap igloo-evap-sa igloo-two-materials igloo-burn igloo-al igloo-al-release
  ice-exchange ice-burn ice-burnout ice-burnout-steady ice-al ice-igloo-burn
  ch-gas-ice ch-gas-igloo ch-steady refusals
  twins burnout routing ini-reread floor-guard multigrid chimera ice-bc mass-flux-inlet restart)
exe=$(cd "$(dirname "$exe")" && pwd)/$(basename "$exe")
mkdir -p "$scratch"; scratch=$(cd "$scratch" && pwd)
case $scratch/ in "$(cd "$here/.." && pwd)"/*)
  echo "run.sh: the scratch directory $scratch is inside test/: its run directories would replace the cases"; exit 2;;
esac
# a private copy of the executable, one per invocation (ctest runs several at once): a rebuild does not touch it
copy=$scratch/hydra_MI2.$$; cp -p "$exe" "$copy"; exe=$copy; trap 'rm -f "$copy"' EXIT

failed=()
for c in "${cases[@]}"; do
  src=$here/cases/$c; run=$scratch/$c
  [ -d "$src" ] || { echo "== $c: no such case"; failed+=("$c"); continue; }
  if [ -f "$src/check.sh" ]; then
    OMP_NUM_THREADS=$nt bash "$src/check.sh" "$run" "$exe" > "$scratch/$c.check" 2>&1; rc=$?
  else
    rm -rf "$run"
    bash "$here/run_case.sh" "$src" "$run" "$nt" "$run.log" "$exe" > /dev/null 2>&1; rc=$?
    cp "$run.log" "$run/run.log" 2> /dev/null
    if [ $rc -ne 0 ]; then
      echo "== $c: hydra_MI2 exit $rc (see $run.log)"; tail -3 "$run.log"; failed+=("$c"); continue
    fi
    python3 -B "$src/check.py" "$run" > "$scratch/$c.check" 2>&1; rc=$?
  fi
  grep -v '^PASS' "$scratch/$c.check"
  [ $rc -eq 0 ] || failed+=("$c")
done
echo "== mi2 tests: ${#cases[@]} cases, ${#failed[@]} failed${failed[*]:+: ${failed[*]}} (each case's lines: $scratch/<case>.check)"
[ ${#failed[@]} -eq 0 ]
