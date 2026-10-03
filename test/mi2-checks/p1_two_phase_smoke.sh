#!/bin/bash
# Two-phase smoke test: MOSE + ICE(partE) + IGLOO(partL) on the JPL nozzle, short schedule.
# partL-* = JPL-Lagrangian part-* (bc f3=300 f4=301); partE-* = JPL-Eulerian part-* (bc f3=300 f4=300).
# Asserts: routing table shows both prefixes, both solvers run, OUTPUT carries partE-field* and partL-*.
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P1}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
LAG=$H/test/JPL-Lagrangian-20micron
EUL=$H/test/JPL-Eulerian-20micron
d=$T/${CASE:-jpl-two-phase}; rm -rf "$d"; mkdir -p "$d/INPUT" "$d/OUTPUT"
cp -rp "$LAG/MESH" "$d/"; cp -p "$LAG/JPL_chang.dat" "$d/"
cp -p "$LAG"/INPUT/gas-* "$d/INPUT/"
for f in phase.txt properties.dat bc.txt; do cp -p "$LAG/INPUT/part-$f" "$d/INPUT/partL-$f"; done
for f in phase.txt properties.dat bc.txt ic.tec; do cp -p "$EUL/INPUT/part-$f" "$d/INPUT/partE-$f"; done
# INI: Lagrangian INI, keys renamed, short schedule, ICE enabled, ICE sections appended from the Eulerian INI
# schedule overridable from the environment (P2 uses NC="2 2" MI="100 100" II="2 2" IC="50 0")
NC=${NC:-2}; MI=${MI:-100}; II=${II:-2}; IC=${IC:-50}
sed -e "s/^N-coupling *=.*/N-coupling    = $NC/" -e "s/^MOSE-iter *=.*/MOSE-iter     = $MI/" \
    -e "s/^IGLOO-iter *=.*/IGLOO-iter    = $II\nICE-iter      = $IC/" \
    -e 's/^phase = part .*/phase = partL/' "$LAG/input.ini" > "$d/input.ini"
sed -n '/^\[ICE-Parameters\]/,$p' "$EUL/input.ini" | sed -e 's/^phase = part .*/phase = partE/' \
    -e 's/^sol-diter *=.*/sol-diter   = 50   ; smoke: ICE writes its field only on this cadence (review S3 #11)/' >> "$d/input.ini"
grep -n "^phase\|^ICE-iter\|^IGLOO-iter\|^\[ICE-\|^\[IGLOO-General" "$d/input.ini"
( cd "$d" && ls INPUT/*phase.txt > filelist.txt && ulimit -s unlimited && OMP_NUM_THREADS=${NT:-2} KMP_STACKSIZE=100M "$EXE" > run.log 2>&1 ); rc=$?
echo "exit=$rc"
grep -n "phase routing" -A 4 "$d/run.log"
grep -n "ERROR\|WARNING\|DEPRECATED" "$d/run.log" | head -5
grep -n "run ICE\|run IGLOO\|ICE   blocks\|IGLOO geometry" "$d/run.log"
echo "--- OUTPUT:"; ls "$d/OUTPUT" | sed 's/-sweep.*//' | sort -u | tr '\n' ' '; echo
[ $rc = 0 ] && ls "$d/OUTPUT" | grep -q '^partE-field' && ls "$d/OUTPUT" | grep -q '^partL-' && echo "TWO-PHASE SMOKE: PASS" || echo "TWO-PHASE SMOKE: FAIL"
