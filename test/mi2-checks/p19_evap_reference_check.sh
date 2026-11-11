#!/bin/bash
# P19: the evaporation reference check (test/evap-box/check_reference.py, REFERENCE.md) as a gate. The four
# evaporating boxes as shipped, at full length and one thread: evap-box (steady), evap-box-ord2 (gas-order 2),
# evap-box-ta (time accurate) and evap-box-sa (RANS); then check_reference.py on each run with its log. It integrates
# each stream's droplet ODE again (CEM, Langmuir-Knudsen interface, burnout at min(1e-6 m0, 1e-20 kg)) in the gas
# IGLOO sampled and compares the exit rows, every scatter sample, the plateau, the sample clock, the source field and
# the Eulerian integrals (rows E S Q C G), and the run log (row L: no LK-interface: warning of IGLOO).
# Per box: the run clean (exit 0, GAME OVER, none of the failure signatures), and the reference check with no failed
# row (its count and the first failed row quoted otherwise). RED on IGLOO bdb4b91, whose unconverged Langmuir-Knudsen
# iteration burnt the drops out 4.3-4.9 us early at 297.8-299.1 K against the reference's 329.7-334.7 K: the E rows
# of every stream fail (8, 4, 4, 4 rows). One core: the four runs about a minute, the reference check one to two
# minutes per box on an idle core (19 to 20 minutes in all at a load of 85 to 120 on 96 cores).
#   usage: p19_evap_reference_check.sh [scratch-root] [exe]
set -u
export LC_ALL=C
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P19}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
RUN=$H/test/mi2-checks/run_case.sh
CR=$H/test/evap-box/check_reference.py
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
echo "  exe $EXE ($(md5sum < "$EXE" | cut -c1-8))"
t0=$(date +%s)

for c in evap-box evap-box-ord2 evap-box-ta evap-box-sa; do
  d=$T/$c; rm -rf "$d"
  "$RUN" "$H/test/$c" "$d" 1 "$T/$c.log" "$EXE" > /dev/null 2>&1; rc=$?
  bad=$(grep -c -e '\[ERROR\]' -e ' FAIL' -e 'on-physical state' -e 'Integration failed at' -e 'p<0 detected' -e 'ABORT:' \
                -e 'LK-interface:' "$T/$c.log")
  if [ $rc = 0 ] && grep -q 'GAME OVER' "$T/$c.log" && [ "$bad" = 0 ]; then
    ok "$c: clean run (exit 0, GAME OVER, no failure signature; $(grep -o 'wall=[0-9.]* s' "$T/$c.log"))"
  else
    ko "$c: run not clean (exit=$rc, GAME OVER $(grep -c 'GAME OVER' "$T/$c.log"), failure signatures $bad: $(grep -o -e '\[ERROR\]' -e 'Integration failed at' -e 'p<0 detected' -e 'on-physical state' -e 'ABORT:' -e 'LK-interface:' "$T/$c.log" | sort -u | tr '\n' ' '))"
  fi
  python3 -B "$CR" "$d" --log "$T/$c.log" > "$T/$c.reference.out" 2>&1; crc=$?
  nf=$(grep -c '^ *\[FAIL\]' "$T/$c.reference.out")
  last=$(grep -E '^(PASS|FAIL): [0-9]+ check' "$T/$c.reference.out" | tail -1)
  if [ $crc = 0 ] && [ "$nf" = 0 ] && [ "${last#PASS}" != "$last" ]; then
    ok "$c: reference check, no failed row (E S Q C G L; $T/$c.reference.out)"
  elif [ "$nf" -gt 0 ]; then
    ko "$c: reference check, $nf failed row(s); first: $(grep -m1 '^ *\[FAIL\]' "$T/$c.reference.out" | sed 's/^ *//' | cut -c1-160)"
  else
    ko "$c: reference check did not run (exit $crc): $(tail -1 "$T/$c.reference.out" | cut -c1-160)"
  fi
done

echo "  wall $(( $(date +%s) - t0 )) s"
echo "== $pass passed, $fail failed"; [ $fail = 0 ]
