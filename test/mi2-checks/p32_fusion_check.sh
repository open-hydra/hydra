#!/bin/bash
# fusion gate (test/ice-al-freeze-box): the smoke of burning aluminium, liquid alumina, cools below its melting point
# (2327 K, where AL2O3(a) of the database ends) and freezes; ICE's solidification model on the smoke family carries
# the heat of fusion to the gas (ATLAS GPB stops a smoke whose table crosses its melting point without it).
#  freeze : the case as it is: aluminium burns direct at 2800 K in a gas held at 1750 K, the smoke born above T-melt,
#           ending frozen at the gas's temperature; check.py: the books (B1, B3, B4), the crossing (F0), the frozen
#           smoke (F1), and the books with the database's solid AL2O3(a) (F2: the gas received the heat of fusion, to
#           the constant cp-solid's error, 2.5e-2 of it at 1750 K).
#  liquid : the same with the smoke's line of INPUT/part-phase.txt without its solidification tokens (the smoke on the
#           liquid's table alone, as GPB wrote it before the melting point was read): the books still close, and F2b
#           must fail by more than half the heat of fusion (0.89 measured): the gate sees the heat of fusion.
#   usage: p32_fusion_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P32}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/ice-al-freeze-box
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
echo "  exe: $EXE ($(md5sum < "$EXE" | cut -c1-8))"

grep -q '^AL2O3(L) 1 combustion=none solidification=on T-melt=2327 ' "$S/INPUT/part-phase.txt" \
  && ok "the smoke AL2O3(L), family 2, with ICE's solidification model, T-melt 2327 K" || ko "the smoke has no solidification model"

echo "--- freeze: the case as it is ---"
d=$T/freeze; log=$T/freeze.log; rm -rf "$d"
"$RUN" "$S" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
[ $rc = 0 ] && grep -q "GAME OVER" "$log" && ok "freeze: run exits 0 with GAME OVER" || ko "freeze: run exit=$rc"
o=$(python3 -B "$S/check.py" "$d" "$log" "$H/test/ice-box" 2>&1)
pass=$((pass + $(echo "$o" | grep -c '^PASS'))); fail=$((fail + $(echo "$o" | grep -c '^FAIL')))
echo "$o" | grep -v '^== check.py'
echo "$o" | grep -q '^== check.py' || ko "freeze: check.py did not finish: $(echo "$o" | tail -1)"

echo "--- liquid: the smoke on the liquid's table alone ---"
v=$T/liquid-case; rm -rf "$v" "$T/liquid"; mkdir -p "$v"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$v/"
sed -i 's/^\(AL2O3(L) 1 combustion=none\) .*$/\1/' "$v/INPUT/part-phase.txt"
if grep -q '^AL2O3(L) 1 combustion=none$' "$v/INPUT/part-phase.txt"; then
  "$RUN" "$v" "$T/liquid" 1 "$T/liquid.log" "$EXE" > /dev/null 2>&1; rc=$?
  [ $rc = 0 ] && grep -q "GAME OVER" "$T/liquid.log" && ok "liquid: run exits 0 with GAME OVER" || ko "liquid: run exit=$rc"
  o=$(python3 -B "$S/check.py" "$T/liquid" "$T/liquid.log" "$H/test/ice-box" 2>&1)
  echo "$o" | grep -q '^PASS  B3 energy' && ok "liquid: the books close without the heat of fusion (B3)" \
                                         || ko "liquid: B3 $(echo "$o" | grep 'B3' | cut -c1-120)"
  r=$(echo "$o" | sed -n 's/^FAIL  F2b .* = \([0-9.e+-]*\) (tol.*/\1/p')
  if [ -n "$r" ] && python3 -c "import sys; sys.exit(0 if float('$r') > 0.5 else 1)"; then
    ok "liquid: F2b fails by $r of M_s h-fus (> 0.5): the gas misses the heat of fusion"
  else ko "liquid: F2b did not fail by more than half the heat of fusion (${r:-absent})"; fi
else ko "liquid: the variant could not be made"; fi

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
