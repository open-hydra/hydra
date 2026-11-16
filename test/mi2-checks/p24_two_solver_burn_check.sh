#!/bin/bash
# two-solver burning gate (test/ice-igloo-burn-box): one run with two burning phases on two models, an Eulerian phase
# solved by ICE (combustion = Beckstead) and a Lagrangian phase solved by IGLOO (combustion = direct), each with its
# own material line, its own products file and its own record, which gives a species of the gas that only that phase
# makes. check.py compares what each solver hands the gas with its own closed form, the gas's energy rise of every
# block with the sum, the coupler's closed-domain books, the mass each of the two species gained with what its phase
# burnt, and the warming of the gas with the first law of the rigid box.
#   usage: p24_two_solver_burn_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P24}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/ice-igloo-burn-box
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
echo "  exe: $EXE ($(md5sum < "$EXE" | cut -c1-8))"

for f in partE-phase.txt partL-phase.txt partE-models.txt partL-models.txt partE-products.txt partL-products.txt; do
  [ -s "$S/INPUT/$f" ] || ko "INPUT/$f is missing"
done
d=$T/run; log=$T/run.log; rm -rf "$d"
"$RUN" "$S" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
echo "  exit=$rc  $(grep -m1 '^# wall' "$log")"
[ $rc = 0 ] && ok "run exits 0" || ko "run exit=$rc: $(grep -m1 -e 'ERROR' -e 'Wrong combustion' "$log" | cut -c1-160)"
grep -q "GAME OVER" "$log" && ok "GAME OVER printed" || ko "no GAME OVER"
bad=$(grep -c -e '\[ERROR\]' -e ' FAIL' -e 'on-physical state' -e 'Integration failed at' -e 'p<0 detected' \
              -e 'ABORT:' -e 'stuck in cell' -e 'no net progress' -e 'hit outer maxIter' -e 'floor guard' "$log")
[ "$bad" = 0 ] && ok "log clean (no error, failure, parcel-abort or floor-guard line)" || ko "log has $bad error/abort line(s)"
n=$(grep -c -- "- A \[combustion\] :" "$log")
[ "$n" = 2 ] && ok "the coupler reports a combustion record for each of the two phases" || ko "$n combustion records reported, not 2"
grep -q "combustion = Beckstead" "$log" && grep -q "direct d2 burn law" "$log" \
  && ok "ICE reports Beckstead and IGLOO reports direct" || ko "the two solvers do not report the two models"
out=$(python3 -B "$S/check.py" "$d" "$log" "$H/test/ice-box" 2>&1); echo "$out" | grep -v '^== check.py'
pass=$((pass + $(echo "$out" | grep -c '^PASS'))); fail=$((fail + $(echo "$out" | grep -c '^FAIL')))
echo "$out" | grep -q '^== check.py' || ko "check.py did not finish: $(echo "$out" | tail -1)"

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
