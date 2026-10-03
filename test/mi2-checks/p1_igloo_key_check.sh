#!/bin/bash
# IGLOO standalone phase-key check on tests/standard/drag-stokes:
#   unnamed  : INPUT/{phase.txt,bc.txt,properties.dat} as shipped (prefix '')
#   keyed    : files renamed INPUT/A-*, [IGLOO-General] phase = A
#   badkey   : phase = nosuch -> must fail
# 1 thread; unnamed vs keyed OUTPUT must be byte-identical (except the A- prefix on file names).
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P1}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
S=$H/src/lib/IGLOO/tests/standard/drag-stokes
EXE=$H/src/lib/IGLOO/bin/IGLOO
for v in unnamed keyed badkey; do
  d=$T/igloo-$v; rm -rf "$d"; mkdir -p "$d/OUTPUT"
  cp -rpL "$S/INPUT" "$S/input.ini" "$d/"
  case $v in
    keyed)  for f in phase.txt bc.txt properties.dat; do mv "$d/INPUT/$f" "$d/INPUT/A-$f"; done
            sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = A/' "$d/input.ini" ;;
    badkey) sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = nosuch/' "$d/input.ini" ;;
  esac
  ( cd "$d" && ulimit -s unlimited && OMP_NUM_THREADS=1 "$EXE" > run.log 2>&1 ); echo "$v exit=$?"
done
grep -n "Error reading phase\|phase" "$T/igloo-badkey/run.log" | tail -2
ls "$T/igloo-unnamed/OUTPUT" | head -5; ls "$T/igloo-keyed/OUTPUT" | head -5
ok=1
for f in "$T"/igloo-unnamed/OUTPUT/*; do
  b=$(basename "$f"); g="$T/igloo-keyed/OUTPUT/A-$b"
  [ -f "$g" ] || { echo "missing $g"; ok=0; continue; }
  cmp -s "$f" "$g" || { echo "DIFF $b"; ok=0; }
done
[ $ok = 1 ] && echo "IGLOO keyed vs unnamed: OUTPUT byte-identical (modulo A- prefix)" || echo "IGLOO keyed vs unnamed: DIFFER"
