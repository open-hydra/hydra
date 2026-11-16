#!/bin/bash
# two-solver burning gate (test/ice-igloo-burn-box): one run with two burning phases on two models, an Eulerian phase
# solved by ICE (combustion = Beckstead) and a Lagrangian phase solved by IGLOO (combustion = direct), each with its
# own material line, its own products file and its own record, which gives a species of the gas that only that phase
# makes. check.py compares what each solver hands the gas with its own closed form, the gas's energy rise of every
# block with the sum, the coupler's closed-domain books, the mass each of the two species gained with what its phase
# burnt, and the warming of the gas with the first law of the rigid box.
#  default: the material A is in both solvers and burns on two models; both phase files say so and nothing is
#           printed. The same case with ICE's model given by [ICE-Physics] combustion instead of its material line
#           writes the same bytes, and the coupler prints one warning that names the material and the two models.
#   usage: p24_two_solver_burn_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P24}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/ice-igloo-burn-box
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
echo "  exe: $EXE ($(md5sum < "$EXE" | cut -c1-8))"

for f in partE-phase.txt partL-phase.txt partE-models.txt partL-models.txt partE-products.txt partL-products.txt; do
  [ -s "$S/INPUT/$f" ] || ko "INPUT/$f is missing"
done
d=$T/run; log=$T/run.log; rm -rf "$d"
"$RUN" "$S" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
echo "  exit=$rc  $(grep -m1 '^# wall' "$log")"
[ $rc = 0 ] && ok "run exits 0" || ko "run exit=$rc: $(grep -m1 -e 'ERROR' -e 'Wrong combustion' "$log" | cut -c1-160)"
grep -q "GAME OVER" "$log" && ok "GAME OVER printed" || ko "no GAME OVER"
bad=$(grep -c -e '\[ERROR\]' -e ' FAIL' -e 'on-physical state' -e 'Integration failed at' -e 'p<0 detected' \
              -e 'ABORT:' -e 'stuck in cell' -e 'no net progress' -e 'hit outer maxIter' -e 'floor guard' "$log")
[ "$bad" = 0 ] && ok "log clean (no error, failure, parcel-abort or floor-guard line)" || ko "log has $bad error/abort line(s)"
n=$(grep -c -- "- A \[combustion\] :" "$log")
[ "$n" = 2 ] && ok "the coupler reports a combustion record for each of the two phases" || ko "$n combustion records reported, not 2"
grep -q "combustion = Beckstead" "$log" && grep -q "direct d2 burn law" "$log" \
  && ok "ICE reports Beckstead and IGLOO reports direct" || ko "the two solvers do not report the two models"
out=$(python3 -B "$S/check.py" "$d" "$log" "$H/test/ice-box" 2>&1); echo "$out" | grep -v '^== check.py'
pass=$((pass + $(echo "$out" | grep -c '^PASS'))); fail=$((fail + $(echo "$out" | grep -c '^FAIL')))
echo "$out" | grep -q '^== check.py' || ko "check.py did not finish: $(echo "$out" | tail -1)"

# --- default: one material name in the two solvers, the model of one of them from its section ---
W='burns as Beckstead in ICE and as direct in IGLOO'
n=$(grep -c "$W" "$log")
[ "$n" = 0 ] && ok "default: both phase files say the model: no warning" || ko "default: the case as it is prints $n warning(s)"
v=$T/default; vlog=$T/default.log; rm -rf "$v"; mkdir -p "$v"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$v/"
sed -i 's/^A 1 combustion=Beckstead$/A 1/' "$v/INPUT/partE-phase.txt"
sed -i 's/^\[ICE-Physics\]$/[ICE-Physics]\ncombustion = Beckstead/' "$v/input.ini"
if grep -q '^A 1$' "$v/INPUT/partE-phase.txt" && grep -q '^A combustion Beckstead K-burn=' "$v/INPUT/partE-models.txt" \
   && grep -q '^combustion = Beckstead$' "$v/input.ini"; then
  "$RUN" "$v" "$v/run" 1 "$vlog" "$EXE" > /dev/null 2>&1; rc=$?
  n=$(grep -c "$W, and one phase file does not say its model" "$vlog"); k=$(ls "$d/OUTPUT" | wc -l)
  [ $rc = 0 ] && [ "$n" = 1 ] && ok "default: ICE's model from [ICE-Physics]: one warning that names the material and the two models" \
                             || ko "default: exit=$rc, $n warning(s)"
  if [ "$k" -gt 0 ] && diff -rq "$d/OUTPUT" "$v/run/OUTPUT" > /dev/null; then ok "default: the same $k OUTPUT files"
  else ko "default: the OUTPUT differs: $(diff -rq "$d/OUTPUT" "$v/run/OUTPUT" 2>&1 | head -1)"; fi
else ko "default: could not move ICE's model to [ICE-Physics]"; fi

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
