#!/bin/bash
# aluminium smoke gate: liquid aluminium (AL(L) of ATLAS's database, with its tables) burns with steam,
# 2 AL(L) + 3 H2O => AL2O3(L) + 3 H2, and the liquid alumina AL2O3(L) is the smoke, born in an ICE family (the tail of
# the combustion record ATLAS GPB writes for a material of an eulerian phase). Four cases, a pair per solver:
#   IGLOO  test/igloo-al-box (Beckstead, the oxide cap, xi-cap 0.2, its residue left; the rest of the oxide the smoke)
#          and test/igloo-al-box-direct (direct at T-burn down to d-stop); the smoke in the ICE phase partE
#   ICE    test/ice-al-box (Beckstead down to d-stop) and test/ice-al-box-direct; the smoke in family 2
# and a variant made here: igloo-al-box with residue = off, its cap handed at burnout born as smoke through the
# release record of the oxide. Each against p30_al_smoke.py: the IGLOO cases close mass, species, momentum and energy
# over the gas, the smoke ICE holds and the residue; the ICE cases close one ICE step (aluminium lost = the gas's share
# + the smoke, in mass, momentum and energy), the smoke's birth state and number, the gas's species over cycle 2.
#   usage: p30_al_smoke_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P30}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
RUN=$H/test/mi2-checks/run_case.sh
PY=$H/test/mi2-checks/p30_al_smoke.py
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
take() { pass=$((pass + $(echo "$1" | grep -c '^PASS'))); fail=$((fail + $(echo "$1" | grep -c '^FAIL'))); echo "$1" | grep -v '^== p30'
         echo "$1" | grep -q '^== p30' || ko "$2: p30_al_smoke.py did not finish: $(echo "$1" | tail -1)"; }
echo "  exe: $EXE ($(md5sum < "$EXE" | cut -c1-8))"

one() {  # one <mode> <name> <case dir>
  local d=$T/$2 log=$T/$2.log
  rm -rf "$d"
  "$RUN" "$3" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; local rc=$?
  [ $rc = 0 ] && grep -q "GAME OVER" "$log" && ok "$2: run exits 0 with GAME OVER" \
                                           || ko "$2: run exit=$rc: $(grep -a -m1 'ERROR' "$log" | cut -c1-200)"
  local o; o=$(python3 -B "$PY" "$1" "$d" "$log" "$H/test/evap-box-mom" 2>&1); take "$o" "$2"
}
for c in igloo-al-box igloo-al-box-direct; do echo "--- $c ---"; one igloo $c "$H/test/$c"; done
for c in ice-al-box ice-al-box-direct; do echo "--- $c ---"; one ice $c "$H/test/$c"; done

echo "--- release: igloo-al-box with residue = off ---"
v=$T/release-case; rm -rf "$v"; mkdir -p "$v"; cp -rp "$H/test/igloo-al-box/INPUT" "$H/test/igloo-al-box/MESH" "$H/test/igloo-al-box/input.ini" "$v/"
sed -i 's/ residue=on$/ residue=off/' "$v/INPUT/partL-models.txt"; sed -i 's/^residue    = on$/residue    = off/' "$v/input.ini"
if grep -q ' residue=off$' "$v/INPUT/partL-models.txt" && grep -q '^residue    = off$' "$v/input.ini"; then
  one igloo release "$v"
  grep -q "IGLOO material AL(L) hands its cap at burnout to the smoke AL2O3(L) (ICE family 1)" "$T/release.log" \
    && ok "release: the coupler reports the cap's smoke" || ko "release: the coupler does not report the cap's smoke"
else ko "release: the variant could not be made"; fi

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
