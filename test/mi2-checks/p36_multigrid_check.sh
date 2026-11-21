#!/bin/bash
# Grid levels under the coupler: hydra-MI2 couples MOSE's finest grid, and of MOSE's levels only its multigrid cycle.
# Each leg runs test/evap-box with keys added, at 1 thread (ATLAS's coarse boundary file of the case is
# data/evap-box-gas-bc2.txt: BCB with [ATLAS-Parameters] MG-levels = 2):
#  nocycle : [MOSE-Multigrid] levels = 2 without a cycle: the run stops after MOSE's set-up with the coupler's refusal,
#            before any MOSE iteration; MOSE's set-up table shows the 2 levels it read.
#  noiter  : levels = 2, cycle = v, no level1-iter: the run stops with the refusal that names the key.
#  hydra   : [HYDRA-Multigrid] levels = 2: the run stops with its refusal.
#  control : [MOSE-Multigrid] levels = 1 (one coupling step): the run ends cleanly.
#  cycle   : levels = 2, cycle = v, level1-iter set: the run ends cleanly on the finest grid with as many gas
#            iterations as the stock case, which is run beside it.
# RED on a build that accepts every levels > 1: nocycle and noiter run on to GAME OVER, and the cycle leg makes other
# gas iterations than the stock case (MOSE starts on its coarse grid); on a build that refuses every levels > 1 the
# cycle leg stops.
#   usage: p36_multigrid_check.sh [scratch-root] [exe]
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
stopped() { ! grep -q -e 'GAME OVER' -e 'MOSE Grid Level' "$1" && [ "$2" != 0 ]; }
iters() { sed -n 's/^ *Iterations  *\([0-9][0-9]*\) *$/\1/p' "$1" | tail -1; }

# --- nocycle ---
d=$(variant nocycle '[MOSE-Multigrid]\nlevels      = 2\nlevel2-iter = 300\n')
cp -p "$H/test/mi2-checks/data/evap-box-gas-bc2.txt" "$d/INPUT/gas-bc2.txt"
rm -rf "$T/nocycle"; "$RUN" "$d" "$T/nocycle" 1 "$T/nocycle.log" "$EXE" > /dev/null 2>&1; rc=$?
grep -q 'Multigrid levels *2' "$T/nocycle.log" && ok "nocycle: MOSE's set-up read the 2 levels and the coarse boundary file" \
  || ko "nocycle: MOSE's set-up shows no 2 levels ($T/nocycle.log)"
grep -q '\[ERROR\] MI2: \[MOSE-Multigrid\] levels = 2 without a multigrid cycle' "$T/nocycle.log" \
  && ok "nocycle: the coupler refuses a start on MOSE's coarse grids" || ko "nocycle: no refusal of levels = 2 without a cycle"
stopped "$T/nocycle.log" $rc && ok "nocycle: the run stops before any MOSE iteration (exit $rc)" \
  || ko "nocycle: the run went on: exit $rc, $(grep -c 'MOSE Grid Level' "$T/nocycle.log") MOSE level lines, $(grep -c 'GAME OVER' "$T/nocycle.log") GAME OVER"

# --- noiter ---
d=$(variant noiter '[MOSE-Multigrid]\nlevels = 2\ncycle  = v\n')
cp -p "$H/test/mi2-checks/data/evap-box-gas-bc2.txt" "$d/INPUT/gas-bc2.txt"
rm -rf "$T/noiter"; "$RUN" "$d" "$T/noiter" 1 "$T/noiter.log" "$EXE" > /dev/null 2>&1; rc=$?
grep -q '\[ERROR\] MI2: \[MOSE-Multigrid\] level1-iter is not set' "$T/noiter.log" && stopped "$T/noiter.log" $rc \
  && ok "noiter: the cycle without level1-iter is refused (exit $rc)" \
  || ko "noiter: no refusal of the cycle without level1-iter (exit $rc, $(grep -c 'GAME OVER' "$T/noiter.log") GAME OVER)"

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

# --- cycle ---
rm -rf "$T/stock"; "$RUN" "$S" "$T/stock" 1 "$T/stock.log" "$EXE" > /dev/null 2>&1
d=$(variant cycle '[MOSE-Multigrid]\nlevels      = 2\ncycle       = v\nlevel1-iter = 1000000000\n')
cp -p "$H/test/mi2-checks/data/evap-box-gas-bc2.txt" "$d/INPUT/gas-bc2.txt"
rm -rf "$T/cycle"; "$RUN" "$d" "$T/cycle" 1 "$T/cycle.log" "$EXE" > /dev/null 2>&1; rc=$?
[ $rc = 0 ] && grep -q 'GAME OVER' "$T/cycle.log" && ! grep -q -e '\[ERROR' -e 'p<0 detected' "$T/cycle.log" \
  && grep -q "MI2: MOSE's multigrid cycle on 2 levels" "$T/cycle.log" \
  && ok "cycle: levels = 2 with cycle = v runs to its end, the coupler saying what it couples" \
  || ko "cycle: exit $rc ($T/cycle.log)"
n0=$(iters "$T/stock.log"); n1=$(iters "$T/cycle.log"); n2=$(grep -c 'MOSE Grid Level *[2-9]' "$T/cycle.log")
[ "${n0:-0}" -gt 0 ] && [ "$n1" = "$n0" ] && [ "$n2" = 0 ] \
  && ok "cycle: MOSE steps its finest grid only, with the stock case's gas iterations ($n1)" \
  || ko "cycle: ${n1:-no} gas iterations against the stock case's ${n0:-none}, $n2 lines of a coarse level"

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
