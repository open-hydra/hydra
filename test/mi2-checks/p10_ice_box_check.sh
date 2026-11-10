#!/bin/bash
# ice-box gate (test/ice-box): the ICE source of cycle 1 -- one forward-Euler step evaluated at the
# initial state of a uniform x-periodic box -- is frozen through cycle 2 (200 coupling steps of one
# time-accurate MOSE iteration). The gas must gain exactly the momentum and energy the coupling
# report says ICE hands it, and both must equal the closed forms of test/ice-box/README.md:
# check.py compares tau, the momentum and energy gains (total and per step), the non-vacuity of
# both effects and the uniformity of the final gas; this script adds the log hygiene and the
# absence of the mass-channel warning.
#  refuse : the same case with one cycle that advances MOSE and ICE together (MOSE-iter 200, ICE-iter 1,
#           N-coupling 1) while both are time-accurate: the run must stop with "advances MOSE and ICE
#           together" and a non-zero exit, while the case's own schedule (ICE alone, then MOSE alone on
#           the frozen source) is accepted.
#   usage: p10_ice_box_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P10}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/ice-box
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }

d=$T/run; log=$T/run.log; rm -rf "$d"
echo "  exe: $EXE ($(md5sum < "$EXE" | cut -c1-8))"
"$RUN" "$S" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
echo "  exit=$rc  $(grep -m1 '^# wall' "$log")"
[ $rc = 0 ] && ok "run exits 0" || ko "run exit=$rc"
grep -q "GAME OVER" "$log" && ok "GAME OVER printed" || ko "no GAME OVER"
bad=$(grep -c -e '\[ERROR\]' -e ' FAIL' -e 'on-physical state' -e 'Integration failed at' -e 'p<0 detected' -e 'ABORT:' "$log")
if [ "$bad" = 0 ]; then ok "log clean"; else
  ko "log has $bad error line(s): $(grep -m1 -e '\[ERROR\]' -e ' FAIL' -e 'on-physical state' -e 'Integration failed at' -e 'p<0 detected' -e 'ABORT:' "$log")"; fi
grep -q "NO LONGER ZERO" "$log" && ko "the ICE mass-channel warning fired" || ok "ICE mass-channel warning absent"
grep -q "advances MOSE and ICE together" "$log" && ko "the case's own schedule was refused" \
                                               || ok "the case's own schedule (ICE alone, then MOSE alone) is accepted"

out=$(python3 -B "$S/check.py" "$d" "$log" 2>&1); echo "$out" | grep -v '^== check.py'
pass=$((pass + $(echo "$out" | grep -c '^PASS'))); fail=$((fail + $(echo "$out" | grep -c '^FAIL')))

# --- refuse: one cycle advancing MOSE and ICE together while time-accurate ---
d=$T/refuse; log=$T/refuse.log; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
sed -i -e 's/^N-coupling *=.*/N-coupling    = 1/' -e 's/^MOSE-iter *=.*/MOSE-iter     = 200/' \
       -e 's/^ICE-iter *=.*/ICE-iter      = 1/' "$d/input.ini"
grep -q '^MOSE-iter *= *200$' "$d/input.ini" && grep -q '^ICE-iter *= *1$' "$d/input.ini" \
  || ko "could not set the forbidden schedule in $d/input.ini"
"$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
n=$(grep -c "advances MOSE and ICE together" "$log")
echo "  refuse: exit=$rc, refusal lines=$n"
if [ $rc != 0 ] && [ "$n" -ge 1 ]; then ok "a cycle advancing MOSE and ICE together while time-accurate is refused (exit=$rc)"
else ko "the forbidden schedule was not refused (exit=$rc, refusal lines=$n)"; fi

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
