#!/bin/bash
# P17: two injected IGLOO materials. test/evap-box-2mat (evap-box-ta's closed box, time accurate): water drops evaporate
# into H2O and graphite particles burn to CO2 with the gas's O2, both materials on the same four DB streams; time-end
# stops the water drops part-way, so the two lose different amounts. p17_two_materials.py checks, per material, the
# species the gas gained (the closed-domain balance) against the mass IGLOO's parcels of that material lost (its
# outloc records of each sweep, over the MOSE blocks that carried them), per block IGLOO's records against the
# coupler's report, the coupler's own audit, the energy bookkeeping, and that the case discriminates (M7). RED: the
# material index swapped in IGLOOxMOSE (sourceMass(nm+1-imat), a scratch build) fails M4 and M5 on H2O, CO2 and O2
# (5.2e-1) and M6 (9.6e-6). One core, about 10 s.
#   usage: p17_two_materials_check.sh [scratch-root] [exe]
set -u
export LC_ALL=C
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P17}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
RUN=$H/test/mi2-checks/run_case.sh
PY=$H/test/mi2-checks/p17_two_materials.py
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
echo "  exe $EXE ($(md5sum < "$EXE" | cut -c1-8))"

c=evap-box-2mat; d=$T/$c; log=$T/$c.log; rm -rf "$d"
"$RUN" "$H/test/$c" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
bad=$(grep -c -e '\[ERROR\]' -e ' FAIL' -e 'on-physical state' -e 'Integration failed at' -e 'p<0 detected' -e 'ABORT:' \
              -e 'LK-interface:' "$log")
if [ $rc = 0 ] && grep -q 'GAME OVER' "$log" && [ "$bad" = 0 ]; then
  ok "$c: clean run (exit 0, GAME OVER, no failure signature; $(grep -o 'wall=[0-9.]* s' "$log"))"
else
  ko "$c: run not clean (exit=$rc, GAME OVER $(grep -c 'GAME OVER' "$log"), failure signatures $bad: $(grep -a -m1 -e '\[ERROR\]' -e 'ABORT:' "$log" | cut -c1-160))"
fi
out=$(python3 -B "$PY" "$d" "$log" "$H/test/evap-box-mom" 2>&1)
echo "$out" | grep -v '^== p17'
pass=$((pass + $(echo "$out" | grep -c '^PASS'))); fail=$((fail + $(echo "$out" | grep -c '^FAIL')))
echo "$out" | grep -q '^== p17' || ko "$c: p17_two_materials.py did not finish: $(echo "$out" | tail -1)"

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
