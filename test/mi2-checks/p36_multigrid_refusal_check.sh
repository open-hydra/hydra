#!/bin/bash
# Grid levels: hydra-MI2 couples MOSE's finest grid only. Each leg runs test/evap-box with one key added, at 1 thread:
#  mose    : [MOSE-Multigrid] levels = 2, level2-iter = 300, with ATLAS's coarse boundary file of the case
#            (data/evap-box-gas-bc2.txt, BCB with [ATLAS-Parameters] MG-levels = 2): the run stops after MOSE's set-up
#            with the coupler's refusal, before any MOSE iteration; MOSE's set-up table shows the 2 levels it read.
#  hydra   : [HYDRA-Multigrid] levels = 2: the run stops with its refusal.
#  control : [MOSE-Multigrid] levels = 1 (one coupling step): the run ends cleanly.
# RED on a build without the refusal: the mose leg runs to GAME OVER, MOSE stepping its coarse grid while the
# coupler hands IGLOO the finest one (which then never moves when the coarse level spans the run).
#   usage: p36_multigrid_refusal_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P36}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/evap-box
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
echo "  exe $EXE ($(md5sum < "$EXE" | cut -c1-8))"
# variant <name> <text appended to input.ini>: a copy of evap-box with the text added
variant() { local d=$T/$1-src; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
            printf '\n%b' "$2" >> "$d/input.ini"; echo "$d"; }

# --- mose ---
d=$(variant mose '[MOSE-Multigrid]\nlevels      = 2\nlevel2-iter = 300\n')
cp -p "$H/test/mi2-checks/data/evap-box-gas-bc2.txt" "$d/INPUT/gas-bc2.txt"
rm -rf "$T/mose"; "$RUN" "$d" "$T/mose" 1 "$T/mose.log" "$EXE" > /dev/null 2>&1; rc=$?
grep -q 'Multigrid levels *2' "$T/mose.log" && ok "mose: MOSE's set-up read the 2 levels and the coarse boundary file" \
  || ko "mose: MOSE's set-up shows no 2 levels ($T/mose.log)"
grep -q '\[ERROR\] MI2: \[MOSE-Multigrid\] levels = 2: ' "$T/mose.log" && ok "mose: the coupler refuses MOSE's grid levels" \
  || ko "mose: no refusal of [MOSE-Multigrid] levels = 2"
! grep -q -e 'GAME OVER' -e 'MOSE Grid Level' "$T/mose.log" && [ $rc != 0 ] \
  && ok "mose: the run stops before any MOSE iteration (exit $rc)" \
  || ko "mose: the run went on: exit $rc, $(grep -c 'MOSE Grid Level' "$T/mose.log") MOSE level lines, $(grep -c 'GAME OVER' "$T/mose.log") GAME OVER"

# --- hydra ---
d=$(variant hydra '[HYDRA-Multigrid]\nlevels = 2\n')
rm -rf "$T/hydra"; "$RUN" "$d" "$T/hydra" 1 "$T/hydra.log" "$EXE" > /dev/null 2>&1; rc=$?
grep -q 'hydra-MI2 does not support multigrid' "$T/hydra.log" && [ $rc != 0 ] && ok "hydra: [HYDRA-Multigrid] levels = 2 is refused" \
  || ko "hydra: no refusal of [HYDRA-Multigrid] levels = 2 (exit $rc)"

# --- control ---
d=$(variant control '[MOSE-Multigrid]\nlevels = 1\n')
sed -i 's/^N-coupling    = 5$/N-coupling    = 1/' "$d/input.ini"
grep -q '^N-coupling    = 1$' "$d/input.ini" || ko "control: N-coupling could not be cut"
rm -rf "$T/control"; "$RUN" "$d" "$T/control" 1 "$T/control.log" "$EXE" > /dev/null 2>&1; rc=$?
[ $rc = 0 ] && grep -q 'GAME OVER' "$T/control.log" && ! grep -q -e '\[ERROR' -e 'p<0 detected' "$T/control.log" \
  && ok "control: [MOSE-Multigrid] levels = 1 runs to its end" || ko "control: exit $rc ($T/control.log)"

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
