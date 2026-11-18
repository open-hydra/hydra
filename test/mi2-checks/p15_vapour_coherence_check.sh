#!/bin/bash
# P15: vapour coherence of the two condensed solvers. ATLAS GPB writes the vapour data of a dispersed phase as columns
# of <phase>-properties.dat (Psat, Lv, Mv, Cpv); ICE and IGLOO each read them by name, with their INI keys as the
# fallback. Static leg (the checkout's sources, p15_vapour_coherence.py static): the names both classifiers accept,
# the grammar they print, the shared reader routines the same code in both, every column of the tables ATLAS wrote
# for test/ read by both; the names both accept and ATLAS never writes printed as the expected gap. Run leg: the
# table of test/evap-box-ta (water, its four vapour columns) read by IGLOO (evap-box-ta, one coupling step) and by ICE
# (the same INPUT/ with the phase claimed by ICE, MOSE + ICE, one ICE step); each solver's setup lines and the
# coupler's datum lines against the table's columns (p15_vapour_coherence.py run). RED: the Lv column renamed in
# ICE's copy of the table and latent-heat given to ICE (a column one solver ignores, its fallback taken): R0, R1,
# R2 fail. One core, about 10 s.
#   usage: p15_vapour_coherence_check.sh [scratch-root] [exe]
set -u
export LC_ALL=C
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P15}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
RUN=$H/test/mi2-checks/run_case.sh
PY=$H/test/mi2-checks/p15_vapour_coherence.py
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
take() { pass=$((pass + $(echo "$1" | grep -c '^PASS'))); fail=$((fail + $(echo "$1" | grep -c '^FAIL'))); echo "$1" | grep -v '^== p15'
         echo "$1" | grep -q '^== p15' || ko "$2: p15_vapour_coherence.py did not finish: $(echo "$1" | tail -1)"; }
echo "  exe $EXE ($(md5sum < "$EXE" | cut -c1-8))"

echo "--- static: the two readers' sources ---"
take "$(python3 -B "$PY" static "$H" 2>&1)" static

echo "--- run: evap-box-ta's table through IGLOO and through ICE ---"
src=$H/test/evap-box-ta
vi=$T/case-igloo; vc=$T/case-ice; rm -rf "$vi" "$vc"; mkdir -p "$vi" "$vc"
cp -rp "$src/INPUT" "$src/MESH" "$vi/"; cp -rp "$src/INPUT" "$src/MESH" "$vc/"
sed 's/^N-coupling    = 20$/N-coupling    = 1/' "$src/input.ini" > "$vi/input.ini"
sed -i '1s/ modeling=lagrangian$/ modeling=eulerian/' "$vc/INPUT/part-phase.txt"
cat > "$vc/input.ini" << 'EOF'
; p15: evap-box-ta's INPUT/ with its water phase claimed by ICE; one ICE step at the initial state
[HYDRA-Parameters]
N-coupling    = 1
MOSE-iter     = 0
ICE-iter      = 1
relaxation    = false
time-accurate = true

[HYDRA-IO]
ini-diter = 1000000000

[MOSE-Parameters]
res-threshold = 1e-30

[MOSE-Physics]
equations = euler

[MOSE-Numerics]
cfl                   = 0.6
time-scheme           = RK3
time-accurate         = true
integration-variables = cons
space-reconstruction  = MUSCL
flux-limiter          = vanleer
riemann-solver        = HLLC

[MOSE-IO]
shell-diter   = 50
sol-diter     = 200
sol-format    = tecplot ascii
sol-variables = thermo

[ICE-Parameters]
phase         = part
res-threshold = 0

[ICE-Numerics]
time-scheme          = euler
cfl                  = 0.4
time-accurate        = true
space-reconstruction = first-order

; the models and constants of evap-box-ta's [IGLOO-*] sections; Lv, Mv and the vapour cp from the table
[ICE-Physics]
drag                  = Stokes
heat-transfer         = Ranz-Marshall
emissivity            = 0
evaporation           = CEM
evaporation-interface = LK
boiling-temperature   = 373.15
lewis-number          = 1.0
vapour-mass-fraction  = 0.0

[ICE-Family1]
closure = MK

[ICE-IO]
shell-diter = 1
sol-diter   = 1
sol-format  = tecplot ascii
ic-format   = tecplot ascii
EOF
if grep -q '^N-coupling    = 1$' "$vi/input.ini" && head -1 "$vc/INPUT/part-phase.txt" | grep -q 'modeling=eulerian$'; then
  "$RUN" "$vi" "$T/igloo" 1 "$T/igloo.log" "$EXE" > /dev/null 2>&1; ri=$?
  "$RUN" "$vc" "$T/ice" 1 "$T/ice.log" "$EXE" > /dev/null 2>&1; rc=$?
  echo "  IGLOO run exit=$ri ($(grep -o 'wall=[0-9.]* s' "$T/igloo.log")), ICE run exit=$rc ($(grep -o 'wall=[0-9.]* s' "$T/ice.log"))"
  take "$(python3 -B "$PY" run "$T/igloo" "$T/igloo.log" "$T/ice" "$T/ice.log" 2>&1)" run
else
  ko "run: the two variants of evap-box-ta could not be made"
fi

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
