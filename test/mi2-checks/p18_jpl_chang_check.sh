#!/bin/bash
# P18: the JPL nozzle against Chang, Hourng & Chien (1996), Fig. 11 (JPL_chang.dat; zone MPL, 10 micron radius, is the
# case the runs reproduce) as a gate. The three JPL cases as shipped (JPL-Lagrangian-20micron: IGLOO,
# JPL-Eulerian-20micron: ICE, JPL-3solver-20micron: both, 0.15 per phase), each with two exit probes of the gas Mach
# number added to a scratch copy of its input.ini (jpl_chang_gate.py --inject), run at full length at the same time,
# P18_THREADS threads each (default 4); then jpl_chang_gate.py --live judges the rows R1-R6 of its docstring against
# the values of its MEASURED block (one-thread runs of the three cases). Byte identity of these runs is the business of
# the one-thread references, not of this gate: at 4 threads every judged row stays within 1.6e-7 of one thread.
# RED on: MOSE alone (R3, R4, R5 and R6 of the ICE case: axis exit M +1.08, margins 0.06 and 0.35); particles of
# dp 2e-6 in the ICE case (R3, R4, R5, R6: margins 0.06 and 0.38, axis mean -0.19); 0.17647 per phase, 27 % more
# loading (R3, R4 axis 1.20, R5 axis mean -0.116 against -0.027, R6 loading 0.5455); NoDrag in IGLOO's phase of the
# three-solver case (R3, R4 axis 1.83, R5, R6 SUM rho_p 1.37). Ranz-Marshall in ICE stays green. Twice as long: green.
# The three runs take 520 to 550 s at a load of about 20 on 96 cores.
#   usage: p18_jpl_chang_check.sh [scratch-root] [exe]
set -u
export LC_ALL=C
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P18}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
NT=${P18_THREADS:-4}
RUN=$H/test/mi2-checks/run_case.sh
GATE=$H/test/mi2-checks/jpl_chang_gate.py
CASES="JPL-Lagrangian-20micron JPL-Eulerian-20micron JPL-3solver-20micron"
mkdir -p "$T/src"
echo "  exe $EXE ($(md5sum < "$EXE" | cut -c1-8)); $NT threads per run, the three runs at once"
t0=$(date +%s); pids=""
for c in $CASES; do
  s=$T/src/$c; rm -rf "$s" "$T/$c" "$T/$c.log"; mkdir -p "$s"
  cp -rp "$H/test/$c/INPUT" "$H/test/$c/MESH" "$H/test/$c/input.ini" "$s/"
  python3 -B "$GATE" --inject "$s/input.ini" || { echo "== 0 passed, 1 failed"; exit 1; }
  "$RUN" "$s" "$T/$c" "$NT" "$T/$c.log" "$EXE" > /dev/null 2>&1 &
  pids="$pids $!"
done
for p in $pids; do wait "$p"; done
echo "  wall of the three runs $(( $(date +%s) - t0 )) s"
python3 -B "$GATE" --live --require-probes IGLOO="$T/JPL-Lagrangian-20micron" ICE="$T/JPL-Eulerian-20micron" \
  both="$T/JPL-3solver-20micron"
