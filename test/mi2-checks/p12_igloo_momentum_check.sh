#!/bin/bash
# Momentum and energy bookkeeping of IGLOO's source coupling (IGLOOxMOSE with IGLOO-iter 1) on evap-box-mom:
# non-evaporating drops decelerate in a box periodic in x and y, so the gas momentum changes only through the
# source. From the coupling report of every MOSE block that had IGLOO sources in force, check.py verifies
#   px1 - px0 = Fx (t1 - t0) to 1e-8 (y, z to 1e-8 of the x term), E1 - E0 = Q (t1 - t0) to 1e-5 (MOSE's
#   per-stage energy round trip is the floor), Fx = SUM mdot (u_inj - u_end) and Q = SUM mdot [cp (T0 - T_end)
#   + (u_inj^2 - v_end^2)/2] against the streams and IGLOO's own snapshot end states, non-vacuity bounds.
# A binary that multiplies sourceMom by the cell volume in the source branch fails the x check by exactly
# V_cell (1.25e-7 here) while the reported Fx is unchanged; a binary without the report fails at once.
#   usage: p12_igloo_momentum_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P12}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/evap-box-mom
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
echo "  exe $EXE ($(md5sum < "$EXE" | cut -c1-8))"

# --- the case files: faces 1-4 periodic (201) and 5-6 symmetry (300) for the gas and the dispersed phase ---
codes() { awk -v f="$2" 'NF==6 && $5==f {print $6}' "$S/INPUT/$1-bc.txt" | sort -u | tr '\n' ' '; }
bcok=1
for p in gas part; do
  for f in 1 2 3 4; do [ "$(codes $p $f)" = "201 " ] || { bcok=0; echo "  $p-bc.txt face $f: codes $(codes $p $f)"; }; done
  for f in 5 6;     do [ "$(codes $p $f)" = "300 " ] || { bcok=0; echo "  $p-bc.txt face $f: codes $(codes $p $f)"; }; done
done
[ $bcok = 1 ] && ok "case: faces 1-4 periodic and 5-6 symmetry in gas-bc.txt and part-bc.txt" || ko "case: unexpected BC codes"
grep -q "^evaporation" "$S/input.ini" && ko "case: an evaporation key is set" || ok "case: no evaporation key (drops keep their mass)"

# --- the run ---
d=$T/run; rm -rf "$d"
"$RUN" "$S" "$d" 1 "$T/run.log" "$EXE" > /dev/null 2>&1; rc=$?
echo "  run: exit=$rc  $(grep -o 'wall=[0-9.]* s' "$T/run.log")"
[ $rc = 0 ] && ok "run exits 0" || ko "run exit=$rc"
grep -q "GAME OVER" "$T/run.log" && ok "run reaches GAME OVER" || ko "no GAME OVER in the log"
bad=$(grep -c -e '\[ERROR\]' -e ' FAIL' -e 'Non-physical state' -e 'Integration failed at' -e 'NaN or p<0 detected' \
              -e 'ABORT:' -e 'stuck in cell' -e 'no net progress' -e 'hit outer maxIter' "$T/run.log")
[ "$bad" = 0 ] && ok "log clean (no error, failure or parcel-abort line)" || ko "log has $bad error/abort line(s)"
n=$(grep -c "Stop condition : time-end" "$T/run.log")
[ "$n" -ge 1 ] && ok "every IGLOO sweep ends at time-end ($n sweeps)" || ko "no sweep ended at time-end"

# --- the bookkeeping ---
out=$(python3 -B "$S/check.py" "$d" --log "$T/run.log" 2>&1); crc=$?
echo "$out" | sed 's/^/  /'
pass=$((pass + $(echo "$out" | grep -c '^PASS'))); fail=$((fail + $(echo "$out" | grep -c '^FAIL')))
[ $crc = 0 ] || [ "$(echo "$out" | grep -c '^FAIL')" -gt 0 ] || ko "check.py exit=$crc without a FAIL line"

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
