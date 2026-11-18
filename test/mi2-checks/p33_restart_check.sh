#!/bin/bash
# P33: a MOSE + ICE restart inside hydra-MI2 ([MOSE-Parameters] and [ICE-Parameters] newrun = false) goes on as the
# run in one go. test/JPL-Eulerian-20micron in a short schedule (coupling steps of 100 MOSE and 50 ICE iterations,
# both steady, 1 thread): a case where both solvers step in every coupling step, each reads back its own field file,
# and the coupler holds no store of its own (no evaporation, no burning). The time-accurate MOSE + ICE boxes advance
# one solver per cycle (a cycle advancing both is refused); the steady boxes of p31 burn, and hand removals and smoke.
#  one-go : N-coupling 2 1 2, MOSE-iter 100 0 100, ICE-iter 50 50 50, in one run;
#  part1  : N-coupling 2 (100, 50);
#  part2  : part 1's directory with newrun = false in both solvers, N-coupling 1 2, MOSE-iter 0 100, ICE-iter 50 50.
# The coupling step with MOSE-iter 0 primes the restart: the coupler's ICE source is in neither field file (the
# first MOSE block after a restart would run without it), so ICE recomputes it from its restored state before MOSE
# needs it. p33_restart.py then asserts the files of part 2 against those of the run in one go: every value of
# gas-field.tec and part-field.tec equal, each solver's output number one higher (MI2's final postprocess writes an
# output more than the run made, hydra_MI2.f90 TODO = 2); gas-diagnostic.tec and part-residuals.tec equal; each
# residual history the one-go file with the row of part 1's last iteration once more (the same postprocess writes
# it twice); part 1's fields moved, the restarted counters those of the whole run. sol-diter is the iterations of
# one coupling step, so that part 1 writes its last state (the final postprocess writes a field only on a multiple of
# sol-diter). Every run: GAME OVER, no error line, exit 0.
#  plain  : the restart without the priming step (one-go N-coupling 4, part 2 N-coupling 2) must differ from its run
#           in one go: an XFAIL line with the size; equal fields fail the gate (the coupler then keeps its source and
#           the priming step can go).
#   usage: p33_restart_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P33}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
S=$H/test/JPL-Eulerian-20micron
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
mk() {  # mk <name> <N-coupling> <MOSE-iter> <ICE-iter> -> a copy of the case on that schedule, fields written per step
  local d=$T/$1; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
  sched "$d/input.ini" "$2" "$3" "$4"
  awk -v m=$M -v k=$K '/^\[/ {s=$0} /^sol-diter/ && s=="[MOSE-IO]" {$0="sol-diter     = " m}
                       /^sol-diter/ && s=="[ICE-IO]" {$0="sol-diter   = " k} {print}' "$d/input.ini" > "$d/ini" &&
    mv "$d/ini" "$d/input.ini"
  echo "$d"; }
restart() {  # restart <part-1 run dir> <dir> <N-coupling> <MOSE-iter> <ICE-iter>: part 1's files, newrun = false, run
  rm -rf "$2"; cp -rp "$1" "$2"; sched "$2/input.ini" "$3" "$4" "$5"
  sed -i -e '/^\[MOSE-Parameters\]/a newrun = false' -e '/^\[ICE-Parameters\]/a newrun = false' "$2/input.ini"
  ( cd "$2" && ulimit -s unlimited && export KMP_STACKSIZE=100M OMP_NUM_THREADS=1 &&
    echo "# exe=$EXE restart threads=1 start=$(date -Is)" > ../$(basename "$2").log &&
    /usr/bin/time -f "# wall=%e s  maxrss=%M kB" "$EXE" >> ../$(basename "$2").log 2>&1 )
  local rc=$?; echo "# exit=$rc end=$(date -Is)" >> "$(dirname "$2")/$(basename "$2").log"; return $rc; }
clean() {  # clean <name> <log> <exit>
  local bad; bad=$(grep -c -e '\[ERROR\]' -e 'STOP' -e 'p<0 detected' -e 'NaN' -e 'ABORT:' "$2")
  if [ "$3" = 0 ] && grep -q 'GAME OVER' "$2" && [ "$bad" = 0 ]; then ok "$1: clean run (GAME OVER, no error line, exit 0, $(grep -o 'wall=[0-9.]* s' "$2"))"
  else ko "$1: run not clean (exit=$3, GAME OVER $(grep -c 'GAME OVER' "$2"), error lines $bad: $(grep -m2 -e '\[ERROR\]' -e 'STOP' -e 'p<0 detected' -e 'NaN' "$2" | tr '\n' ' '))"; fi; }
judge() {  # judge <mode> <args...>: the .py's lines counted; a crash counts as a failure
  local out rc; out=$(python3 -B "$PY" "$@" 2>&1); rc=$?; echo "$out"
  pass=$((pass + $(echo "$out" | grep -c '^PASS'))); fail=$((fail + $(echo "$out" | grep -c '^FAIL')))
  if [ $rc != 0 ] && [ "$(echo "$out" | grep -c '^FAIL')" = 0 ]; then ko "$1: the checks did not run (exit $rc)"; fi; }

# --- part 1, shared by both restarts ---
d=$(mk part1 2 $M $K); "$RUN" "$d" "$d/run" 1 "$T/part1.log" "$EXE" > /dev/null 2>&1; clean part1 "$T/part1.log" $?
# --- primed: the run in one go and the restart ---
d=$(mk one-go "2 1 2" "$M 0 $M" "$K $K $K"); "$RUN" "$d" "$d/run" 1 "$T/one-go.log" "$EXE" > /dev/null 2>&1
clean one-go "$T/one-go.log" $?
restart "$T/part1/run" "$T/part2" "1 2" "0 $M" "$K $K"; clean part2 "$T/part2.log" $?
if grep -q 'Initial conditions file *OUTPUT/gas-field.tec' "$T/part2.log"; then ok "part2: MOSE reads its restart file OUTPUT/gas-field.tec"
else ko "part2: MOSE does not name OUTPUT/gas-field.tec as its initial condition"; fi
judge exact "$T/one-go/run/OUTPUT" "$T/part2/OUTPUT" "$T/part1/run/OUTPUT"
# --- plain: no priming step ---
d=$(mk one-go-plain 4 $M $K); "$RUN" "$d" "$d/run" 1 "$T/one-go-plain.log" "$EXE" > /dev/null 2>&1
clean one-go-plain "$T/one-go-plain.log" $?
restart "$T/part1/run" "$T/part2-plain" 2 $M $K; clean part2-plain "$T/part2-plain.log" $?
judge plain "$T/one-go-plain/run/OUTPUT" "$T/part2-plain/OUTPUT"
echo "== $pass passed, $fail failed"; [ $fail = 0 ]
