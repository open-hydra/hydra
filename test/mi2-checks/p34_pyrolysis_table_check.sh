#!/bin/bash
# BCB re-run alone on an imported case with a pyrolysing-ablator wall (BC 506). BCB reads the pyrolysis-gas tables
# of the phase, <phase>pyrolysis.dat (written by GPB), from the case root first; after an import they are in INPUT/.
#  regen   : regen_case.sh <case> <scratch> BCB writes the bc.txt of ATLAS's own reference;
#  figaro  : scripts/modules/simulations/ATLAS.sh BCB (what `figaro ATLAS BCB` runs) does the same, imports it, and
#            leaves no table in the case root.
# Both for the unnamed phase (pyrolysis.dat) and for a named one (gas-pyrolysis.dat). The cases are built from the
# fixtures of ATLAS's tests IG-carbon-phenolic and DP-gsi-manifold. Without INPUT/*pyrolysis.dat copied back, BCB
# stops: "[ERROR] key pyrolysis-gas of section [...]: <phase>pyrolysis.dat is neither in the case directory nor in ...".
# Needs a built ATLAS that knows BC 506 (ATLASDIR, default the checked-out utils/ATLAS).
#   usage: p34_pyrolysis_table_check.sh [scratch-root]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P34}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
ATLASDIR=${ATLASDIR:-$H/utils/ATLAS}; export ATLASDIR
REGEN=$H/test/mi2-checks/regen_case.sh
FIGARO_ATLAS=$H/scripts/modules/simulations/ATLAS.sh
F=$ATLASDIR/test/BCB
pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
if [ ! -x "$ATLASDIR/bin/BCB" ]; then echo "SKIP  no built BCB at $ATLASDIR/bin: nothing was checked"; exit 0; fi
for d in IG-carbon-phenolic DP-gsi-manifold; do
  [ -f "$F/$d/input.ini" ] || { echo "FAIL  ATLAS's fixture $F/$d is missing"; exit 2; }
done
rm -rf "$T"; mkdir -p "$T"

# an imported case: input.ini, MESH/, and in INPUT/ what an import of GPB + BCB leaves there for BCB to read again
make_case() {   # <fixture> <case-dir>
  local s=$F/$1 c=$2 f
  mkdir -p "$c/INPUT" "$c/MESH"
  cp -L "$s/input.ini" "$c/"; cp -L "$s/mesh.tec" "$c/MESH/"
  for f in "$s"/*phase.txt "$s"/*pyrolysis.dat; do cp -L "$f" "$c/INPUT/"; done
}
# every bc file of ATLAS's reference, byte for byte
same_bc() {     # <fixture> <dir with the bc files> ; prints the first file that differs
  local r=$F/$1/reference f n=0
  for f in "$r"/*bc.txt; do
    n=$((n+1)); cmp -s "$f" "$2/$(basename "$f")" || { basename "$f"; return 1; }
  done
  [ $n -gt 0 ]
}

for fx in IG-carbon-phenolic DP-gsi-manifold; do
  tab=$(cd "$F/$fx" && ls *pyrolysis.dat)
  # --- regen_case.sh
  make_case $fx "$T/$fx-regen"
  "$REGEN" "$T/$fx-regen" "$T/$fx-regen.out" BCB > "$T/$fx-regen.log" 2>&1; rc=$?
  if [ $rc = 0 ] && bad=$(same_bc $fx "$T/$fx-regen.out/fromATLAStoSolver"); then
    ok "$fx ($tab): regen_case.sh BCB writes the bc.txt of ATLAS's reference"
  else
    ko "$fx ($tab): regen_case.sh BCB: exit $rc, ${bad:+$bad differs, }$(grep -m1 '\[ERROR\]' "$T/$fx-regen.out/atlas.log" 2>/dev/null)"
  fi
  # --- figaro's script, answering y to its import question
  c=$T/$fx-figaro; make_case $fx "$c"
  ( cd "$c" && SCRIPT_DIR=$H/scripts bash "$FIGARO_ATLAS" BCB <<< y ) > "$T/$fx-figaro.log" 2>&1
  if bad=$(same_bc $fx "$c/INPUT"); then
    ok "$fx ($tab): figaro's ATLAS.sh BCB imports the bc.txt of ATLAS's reference"
  else
    ko "$fx ($tab): figaro's ATLAS.sh BCB: ${bad:-no bc.txt} not ATLAS's reference in INPUT/, $(grep -m1 '\[ERROR\]' "$T/$fx-figaro.log")"
  fi
  left=$(cd "$c" && ls *pyrolysis.dat *phase.txt 2>/dev/null | tr '\n' ' ')
  [ -z "$left" ] && ok "$fx: the import leaves no table in the case root" || ko "$fx: left in the case root after the import: $left"
  cmp -s "$F/$fx/$tab" "$c/INPUT/$tab" && ok "$fx: INPUT/$tab untouched" || ko "$fx: INPUT/$tab changed or gone"
done
echo "== $pass passed, $fail failed"; [ $fail = 0 ]
