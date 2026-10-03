#!/bin/bash
# ICE phase-key check: JPL-Eulerian short schedule, legacy (the key removed) vs keyed ([ICE-Parameters] phase = part,
# as shipped), 1 thread -> OUTPUT must be byte-identical.
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P1}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
S=$H/test/JPL-Eulerian-20micron
EXE=${2:-$H/bin/hydra_MI2}
for v in legacy keyed; do
  d=$T/jplE-$v
  rm -rf "$d"; mkdir -p "$d"
  cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
  sed -i -e 's/^N-coupling *=.*/N-coupling    = 2/' -e 's/^MOSE-iter *=.*/MOSE-iter     = 100/' \
         -e 's/^ICE-iter *=.*/ICE-iter      = 50/' "$d/input.ini"
  if [ "$v" = legacy ]; then
    sed -i '/^phase *= *part/d' "$d/input.ini"
  fi
  echo "$v: phase keys in the INI: $(grep -c '^phase *=' "$d/input.ini")"
  ( cd "$d" && mkdir -p OUTPUT && ls INPUT/*phase.txt > filelist.txt && ulimit -s unlimited \
    && OMP_NUM_THREADS=1 KMP_STACKSIZE=100M "$EXE" > run.log 2>&1 ); echo "$v exit=$?"
  grep -n "phase routing" -A 3 "$d/run.log" | sed -n 2,4p
  grep -n "DEPRECATED" "$d/run.log" | head -2
done
if diff -rq "$T/jplE-legacy/OUTPUT" "$T/jplE-keyed/OUTPUT"; then echo "ICE keyed vs legacy: OUTPUT byte-identical"; else echo "ICE keyed vs legacy: DIFFER"; fi
