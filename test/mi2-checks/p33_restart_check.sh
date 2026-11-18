#!/bin/bash
# P33: a MOSE + ICE restart inside hydra-MI2 ([MOSE-Parameters] and [ICE-Parameters] newrun = false) goes on as the
# run in one go, every value equal (1 thread). ICE's field file does not hold what its last coupling step hands the
# gas (the source of its last step's stages) nor what its next step starts from (the relaxation times that bound its
# steps, a smoke family's birth rate): MI2 ends a run with ICE by writing them, OUTPUT/<ICE prefix>coupling.tec, and an
# ICE restart reads them back. Three cases on short schedules, four legs, part 1 and its restart for as many coupling steps
# again against the run in one go:
#  jpl     : test/JPL-Eulerian-20micron, MOSE and ICE steady, coupling steps of 100 and 50 iterations (2 + 2 against 4):
#            ICE's source alone;
#  smoke   : test/steady-ice-smoke on the same schedule: a burning family whose smoke ICE births as a rate (srate) and
#            takes out of what the gas receives (birth); the smoke family's steps bounded by its relaxation time;
#  removal : test/steady-ice-burnout with MOSE steady and ICE time-accurate (first order, cfl 0.4, K-burn x100) in
#            alternate cycles, MOSE-iter 100 0, ICE-iter 0 50 (4 cycles + 4 against 8): the removals ICE stores and the
#            time they were made over (goneRate), ICE's source averaged over its block.
#  odd     : the jpl leg with sol-diter 300 (MOSE) and 150 (ICE), which divide neither part's iterations: the end of
#            part 1 is no output iteration, and its final fields are still written and restarted from.
# p33_restart.py judges each restart against its run in one go: every value of gas-field.tec, part-field.tec and
# part-coupling.tec equal, the output numbers included (a variable below 1e-20 in both runs compared as zero: w on a
# 2-D mesh, which ICE zeroes at every set-up); gas-diagnostic.tec and part-residuals.tec equal; each residual history
# equal row for row; part 1's fields moved, the restarted counters those of the whole run. And the end of each run is
# its end-of-run output: no history holds a row twice in a row, each solver's output number is its iterations //
# sol-diter. Every run: GAME OVER, no error line, exit 0.
# Two refusals, each judged on its message and on the run not ending:
#  missing : the jpl restart without part-coupling.tec stops and names the file;
#  igloo   : test/evap-box-ta (MOSE + IGLOO) restarted with [MOSE-Parameters] newrun = false stops: IGLOO has no
#            restart, the restarted run would begin without IGLOO's source.
#   usage: p33_restart_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P33}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
RUN=$H/test/mi2-checks/run_case.sh
PY=$H/test/mi2-checks/p33_restart.py
M=100; K=50
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
echo "  exe $EXE ($(md5sum < "$EXE" | cut -c1-8))"

sched() {  # sched <input.ini> <N-coupling> <MOSE-iter> <ICE-iter>
  sed -i -e "s/^N-coupling *=.*/N-coupling    = $2/" -e "s/^MOSE-iter *=.*/MOSE-iter     = $3/" \
         -e "s/^ICE-iter *=.*/ICE-iter      = $4/" "$1"; }
mk() {  # mk <case> <name> <N-coupling> <MOSE-iter> <ICE-iter> <MOSE sol-diter> <ICE sol-diter> [ini sed] [models sed]
  local d=$T/$2; rm -rf "$d"; mkdir -p "$d"; cp -rp "$1/INPUT" "$1/MESH" "$1/input.ini" "$d/"
  sched "$d/input.ini" "$3" "$4" "$5"
  [ -n "${8:-}" ] && sed -i "$8" "$d/input.ini"
  [ -n "${9:-}" ] && sed -i "$9" "$d/INPUT/part-models.txt"
  awk -v m="$6" -v k="$7" '/^\[/ {s=$0} /^sol-diter/ && s=="[MOSE-IO]" {$0="sol-diter     = " m}
                           /^sol-diter/ && s=="[ICE-IO]" {$0="sol-diter   = " k} {print}' "$d/input.ini" > "$d/ini" &&
    mv "$d/ini" "$d/input.ini"
  echo "$d"; }
restart() {  # restart <part-1 run dir> <dir> <N-coupling> <MOSE-iter> <ICE-iter> [file to remove]: newrun = false, run
  rm -rf "$2"; cp -rp "$1" "$2"; sched "$2/input.ini" "$3" "$4" "$5"; [ -n "${6:-}" ] && rm -f "$2/$6"
  sed -i -e '/^\[MOSE-Parameters\]/a newrun = false' -e '/^\[ICE-Parameters\]/a newrun = false' "$2/input.ini"
  ( cd "$2" && ulimit -s unlimited && export KMP_STACKSIZE=100M OMP_NUM_THREADS=1 &&
    echo "# exe=$EXE restart threads=1 start=$(date -Is)" > ../$(basename "$2").log &&
    if [ -x /usr/bin/time ]; then /usr/bin/time -f "# wall=%e s  maxrss=%M kB" "$EXE" >> ../$(basename "$2").log 2>&1
    else "$EXE" >> ../$(basename "$2").log 2>&1; fi )
  local rc=$?; echo "# exit=$rc end=$(date -Is)" >> "$(dirname "$2")/$(basename "$2").log"; return $rc; }
clean() {  # clean <name> <log> <exit>
  local bad; bad=$(grep -c -e '\[ERROR\]' -e 'STOP' -e 'p<0 detected' -e 'NaN' -e 'ABORT:' "$2")
  if [ "$3" = 0 ] && grep -q 'GAME OVER' "$2" && [ "$bad" = 0 ]; then ok "$1: clean run (GAME OVER, no error line, exit 0, $(grep -o 'wall=[0-9.]* s' "$2"))"
  else ko "$1: run not clean (exit=$3, GAME OVER $(grep -c 'GAME OVER' "$2"), error lines $bad: $(grep -m2 -e '\[ERROR\]' -e 'STOP' -e 'p<0 detected' -e 'NaN' "$2" | tr '\n' ' '))"; fi; }
refused() {  # refused <name> <log> <exit> <message>: the run stops on the message, before GAME OVER
  if [ "$3" != 0 ] && grep -qF -- "$4" "$2" && ! grep -q 'GAME OVER' "$2"; then ok "$1: refused (exit $3): $(grep -m1 -F -- "$4" "$2" | cut -c1-150)"
  else ko "$1: not refused as expected (exit=$3, GAME OVER $(grep -c 'GAME OVER' "$2"), message found $(grep -cF -- "$4" "$2"))"; fi; }
judge() {  # judge <mode> <args...>: the .py's lines counted; a crash counts as a failure
  local out rc; out=$(python3 -B "$PY" "$@" 2>&1); rc=$?; echo "$out"
  pass=$((pass + $(echo "$out" | grep -c '^PASS'))); fail=$((fail + $(echo "$out" | grep -c '^FAIL')))
  if [ $rc != 0 ] && [ "$(echo "$out" | grep -c '^FAIL')" = 0 ]; then ko "$1: the checks did not run (exit $rc)"; fi; }

leg() {  # leg <tag> <case> <part-1 N-coupling> <MOSE-iter> <ICE-iter> <one-go N-coupling> <MOSE-iter> <ICE-iter>
        #     [ini sed] [models sed]: part 1, the run in one go, part 1 restarted on part 1's schedule; judged
  local d t=$1 c=$2 sm=${SDM:-$M} sk=${SDK:-$K}
  d=$(mk "$c" "$t-part1" "$3" "$4" "$5" $sm $sk "${9:-}" "${10:-}"); "$RUN" "$d" "$d/run" 1 "$T/$t-part1.log" "$EXE" > /dev/null 2>&1
  clean "$t-part1" "$T/$t-part1.log" $?
  d=$(mk "$c" "$t-one-go" "$6" "$7" "$8" $sm $sk "${9:-}" "${10:-}"); "$RUN" "$d" "$d/run" 1 "$T/$t-one-go.log" "$EXE" > /dev/null 2>&1
  clean "$t-one-go" "$T/$t-one-go.log" $?
  restart "$T/$t-part1/run" "$T/$t-part2" "$3" "$4" "$5"; clean "$t-part2" "$T/$t-part2.log" $?
  if grep -q 'Initial conditions file *OUTPUT/gas-field.tec' "$T/$t-part2.log" && \
     grep -q 'ICE restarts with what its last coupling step handed the gas, from OUTPUT/part-coupling.tec' "$T/$t-part2.log"
  then ok "$t-part2: MOSE reads OUTPUT/gas-field.tec, MI2 reads ICE's hand-off from OUTPUT/part-coupling.tec"
  else ko "$t-part2: the restart does not name OUTPUT/gas-field.tec and OUTPUT/part-coupling.tec as what it reads"; fi
  judge exact "$T/$t-one-go/run/OUTPUT" "$T/$t-part2/OUTPUT" "$T/$t-part1/run/OUTPUT" $sm $sk; }

FO='s/^time-scheme          = RK3$/time-scheme          = euler/; s/^space-reconstruction = MUSCL$/space-reconstruction = first-order/; /^flux-limiter         = vanleer$/d'
TA='s/^time-accurate        = false$/time-accurate        = true/; s/^cfl                  = 0.2$/cfl                  = 0.4/'
leg jpl     "$H/test/JPL-Eulerian-20micron" 2 $M $K 4 $M $K
SDM=300 SDK=150 leg odd "$H/test/JPL-Eulerian-20micron" 2 $M $K 4 $M $K
leg smoke   "$H/test/steady-ice-smoke"      2 $M $K 4 $M $K
leg removal "$H/test/steady-ice-burnout"    "1 1 1 1" "$M 0 $M 0" "0 $K 0 $K" \
            "1 1 1 1 1 1 1 1" "$M 0 $M 0 $M 0 $M 0" "0 $K 0 $K 0 $K 0 $K" "$TA; $FO" 's/ K-burn=1.39e-06 / K-burn=1.39e-04 /'
# --- missing: the restart without the hand-off ---
restart "$T/jpl-part1/run" "$T/missing" 2 $M $K OUTPUT/part-coupling.tec
refused missing "$T/missing.log" $? 'ICE restarts ([ICE-Parameters] newrun = false) and OUTPUT/part-coupling.tec is missing'
# --- igloo: a restart with IGLOO is refused ---
d=$(mk "$H/test/evap-box-ta" igloo-part1 1 50 0 50 1); "$RUN" "$d" "$d/run" 1 "$T/igloo-part1.log" "$EXE" > /dev/null 2>&1
clean igloo-part1 "$T/igloo-part1.log" $?
restart "$T/igloo-part1/run" "$T/igloo" 1 50 0
refused igloo "$T/igloo.log" $? 'newrun = false with IGLOO: IGLOO has no restart'
echo "== $pass passed, $fail failed"; [ $fail = 0 ]
