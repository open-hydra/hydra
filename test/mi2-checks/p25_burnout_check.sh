#!/bin/bash
# burn-out gate (test/ice-burnout-box, test/igloo-burnout-box): a burning material without a residue burns out in ICE
# and in IGLOO, one burning section in the two cases: what is left of a drop at d-burnout is handed to the gas whole.
#  books  : test/ice-burnout-box. ICE advances alone through the whole burn; the mass, the momentum and the energy its
#           sources hand, summed with the step, are those the cloud held at the start (check.py), the removal fires
#           by the size rule at the first step that leaves the drops at or below d-burnout, the field ends empty.
#  long   : the same case with a step twenty times longer (dt-max 2.0e-4, the relaxation-time ceiling of the step
#           off): too long for the end of the burn. The removal fires at the first step that leaves less than it
#           burnt, well above d-burnout; the books close, and no density is floored.
#  frozen : the same case with d-burnout above the drops and one ICE step: the removal fires in it, its source is the
#           burn alone and the rest of the cloud is stored; frozen through 5 MOSE iterations, the gas gains the frozen
#           burn over their time and the store once (the inter-phase books of a removal).
#  frozen2: the same with ICE-iter 2: the second step finds the cloud empty, and the gas still gains the first one's
#           store.
#  cold   : the two boxes with T-ign above any temperature their drops reach and d-burnout above the drops: no drop
#           burns, none is removed, no mass reaches the gas.
#  twin   : test/igloo-burnout-box against its closed forms (test/igloo-burn-box/check.py: everything a stream
#           carries in ends in the gas), and the energy per kg handed, ICE against IGLOO.
#  refuse : d-stop beside d-burnout on the material line stops each solver before the first step.
#   usage: p25_burnout_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P25}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
RUN=$H/test/mi2-checks/run_case.sh
I=$H/test/ice-burnout-box; L=$H/test/igloo-burnout-box
TOL_E=1e-9        # energy per kg handed, ICE against IGLOO (round-off of the two sums: 1e-15 measured)
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
num() { echo "$1" | grep -m1 "$2" | sed -e "s/.*$2 *//" | awk '{print $1}'; }
close() { python3 -c "import sys; a, b, t = map(float, sys.argv[1:4]); sys.exit(0 if abs(a - b) <= t * abs(b) else 1)" "$1" "$2" "$3"; }
rel() { python3 -c "import sys; a, b = map(float, sys.argv[1:3]); print('%.2e' % (abs(a - b) / abs(b)))" "$1" "$2"; }
take() { pass=$((pass + $(echo "$1" | grep -c '^PASS'))); fail=$((fail + $(echo "$1" | grep -c '^FAIL'))); echo "$1" | grep -v '^== check.py'
         echo "$1" | grep -q '^== check.py' || ko "$2: check.py did not finish: $(echo "$1" | tail -1)"; }
echo "  exe: $EXE ($(md5sum < "$EXE" | cut -c1-8))"

a=$(grep '^A ' "$I/INPUT/part-phase.txt"); b=$(grep '^A ' "$L/INPUT/part-phase.txt")
[ -n "$a" ] && [ "$a" = "$b" ] && ok "one material line in the two phase files ($a)" || ko "the material lines differ: '$a' / '$b'"
a=$(grep '^A ' "$I/INPUT/part-models.txt"); b=$(grep '^A ' "$L/INPUT/part-models.txt")
[ -n "$a" ] && [ "$a" = "$b" ] && echo "$a" | grep -q " d-burnout=" && ! echo "$a" | grep -q " d-stop=" \
  && ok "one row in the two models files, with d-burnout and no d-stop ($a)" || ko "the rows of the models files differ: '$a' / '$b'"

echo "--- books: test/ice-burnout-box ---"
d=$T/books; log=$T/books.log; rm -rf "$d"
"$RUN" "$I" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
echo "  exit=$rc  $(grep -m1 '^# wall' "$log")"
[ $rc = 0 ] && grep -q "GAME OVER" "$log" && ok "books: run exits 0 with GAME OVER" || ko "books: run exit=$rc"
grep -q "no residue: the drops of a cell are removed at d-burnout" "$log" && ok "books: ICE reports the removal at d-burnout" \
                                                                         || ko "books: ICE does not report the removal"
oi=$(python3 -B "$I/check.py" "$d" "$log" "$H/test/ice-box" 2>&1); take "$oi" books

echo "--- long: a step twenty times longer ---"
v=$T/long-case; d=$T/long; log=$T/long.log; rm -rf "$v" "$d"; mkdir -p "$v"; cp -rp "$I/INPUT" "$I/MESH" "$I/input.ini" "$v/"
sed -i -e 's/^dt-max  *= 1.0e-5$/dt-max               = 2.0e-4\ntau-factor           = 0/' \
       -e 's/^N-coupling    = 600 1$/N-coupling    = 60 1/' "$v/input.ini"
if grep -q '^dt-max               = 2.0e-4$' "$v/input.ini" && grep -q '^tau-factor           = 0$' "$v/input.ini" \
   && grep -q '^N-coupling    = 60 1$' "$v/input.ini"; then
  "$RUN" "$v" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  [ $rc = 0 ] && grep -q "GAME OVER" "$log" && ok "long: run exits 0 with GAME OVER" || ko "long: run exit=$rc"
  o=$(python3 -B "$I/check.py" "$d" "$log" "$H/test/ice-box" 2>&1); take "$o" long
  echo "$o" | grep -q 'R4 the removal fires by the rule of the long step' && ok "long: the removal fires by the rule of the long step" \
                                                                         || ko "long: the removal did not fire by the rule of the long step"
else ko "long: the variant could not be made"; fi

echo "--- frozen: the removal's source applied to the gas ---"
v=$T/frozen-case; d=$T/frozen; log=$T/frozen.log; rm -rf "$v" "$d"; mkdir -p "$v"; cp -rp "$I/INPUT" "$I/MESH" "$I/input.ini" "$v/"
sed -i -e 's/^dt-max  *= 1.0e-5$/dt-max               = 1.0e-4/' -e 's/^N-coupling    = 600 1$/N-coupling    = 1 5/' \
       -e 's/^d-burnout  = 3.0e-6$/d-burnout  = 4.0e-5/' -e 's/^sol-diter     = 200$/sol-diter     = 5/' "$v/input.ini"
sed -i 's/ d-burnout=3e-06$/ d-burnout=4e-05/' "$v/INPUT/part-models.txt"
if grep -q '^dt-max               = 1.0e-4$' "$v/input.ini" && grep -q '^N-coupling    = 1 5$' "$v/input.ini" \
   && grep -q '^d-burnout  = 4.0e-5$' "$v/input.ini" && grep -q '^sol-diter     = 5$' "$v/input.ini" \
   && grep -q ' d-burnout=4e-05$' "$v/INPUT/part-models.txt"; then
  "$RUN" "$v" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  [ $rc = 0 ] && grep -q "GAME OVER" "$log" && ok "frozen: run exits 0 with GAME OVER" || ko "frozen: run exit=$rc"
  o=$(python3 -B "$I/check.py" "$d" "$log" "$H/test/ice-box" frozen 2>&1); take "$o" frozen
else ko "frozen: the variant could not be made"; fi

echo "--- frozen2: two ICE steps in the coupling step ---"
v=$T/frozen2-case; d=$T/frozen2; log=$T/frozen2.log; rm -rf "$v" "$d"; mkdir -p "$v"; cp -rp "$T/frozen-case/INPUT" "$T/frozen-case/MESH" "$T/frozen-case/input.ini" "$v/"
sed -i 's/^ICE-iter      = 1 0$/ICE-iter      = 2 0/' "$v/input.ini"
if grep -q '^ICE-iter      = 2 0$' "$v/input.ini" && grep -q ' d-burnout=4e-05$' "$v/INPUT/part-models.txt"; then
  "$RUN" "$v" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  [ $rc = 0 ] && grep -q "GAME OVER" "$log" && ok "frozen2: run exits 0 with GAME OVER" || ko "frozen2: run exit=$rc"
  o=$(python3 -B "$I/check.py" "$d" "$log" "$H/test/ice-box" frozen2 2>&1); take "$o" frozen2
else ko "frozen2: the variant could not be made"; fi

echo "--- cold: drops at or below d-burnout that do not burn ---"
v=$T/cold-case; d=$T/cold; log=$T/cold.log; rm -rf "$v" "$d"; mkdir -p "$v"; cp -rp "$I/INPUT" "$I/MESH" "$I/input.ini" "$v/"
sed -i -e 's/^d-burnout  = 3.0e-6$/d-burnout  = 4.0e-5/' -e 's/^T-ign      = 290$/T-ign      = 500/' "$v/input.ini"
sed -i -e 's/ T-ign=290 / T-ign=500 /' -e 's/ d-burnout=3e-06$/ d-burnout=4e-05/' "$v/INPUT/part-models.txt"
if grep -q '^d-burnout  = 4.0e-5$' "$v/input.ini" && grep -q '^T-ign      = 500$' "$v/input.ini" \
   && grep -q ' T-ign=500 .* d-burnout=4e-05$' "$v/INPUT/part-models.txt"; then
  "$RUN" "$v" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  [ $rc = 0 ] && grep -q "GAME OVER" "$log" && ok "cold: run exits 0 with GAME OVER" || ko "cold: run exit=$rc"
  o=$(python3 -B "$I/check.py" "$d" "$log" "$H/test/ice-box" cold 2>&1); take "$o" cold
else ko "cold: the variant could not be made"; fi
v=$T/cold-igloo-case; d=$T/cold-igloo; log=$T/cold-igloo.log; rm -rf "$v" "$d"; mkdir -p "$v"; cp -rp "$L/INPUT" "$L/MESH" "$L/input.ini" "$v/"
sed -i -e 's/ T-ign=290 / T-ign=500 /' -e 's/ d-burnout=3e-06$/ d-burnout=4e-05/' "$v/INPUT/part-models.txt"
if grep -q ' T-ign=500 .* d-burnout=4e-05$' "$v/INPUT/part-models.txt"; then
  "$RUN" "$v" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  [ $rc = 0 ] && grep -q "GAME OVER" "$log" && ok "cold: igloo-burnout-box run exits 0 with GAME OVER" || ko "cold: igloo-burnout-box run exit=$rc"
  # round-off of IGLOO's in - out sums: 2.5e-21 kg/s measured; a removal hands 4.0e-06 kg/s (an IGLOO that removes parcels that do not burn)
  m=$(grep -a '\[MI2-report\] igloo-to-gas ' "$log" | awk '{ n++; v = $NF; sub(/^-/, "", v); if (v + 0 > big) big = v + 0 } END { printf "%d %.3e %d", n, big, (n > 0 && big <= 1e-15) }')
  [ "${m##* }" = 1 ] && ok "cold: IGLOO hands the gas no mass (${m% *}: blocks, largest kg/s; tol 1e-15)" \
                     || ko "cold: IGLOO hands the gas mass (${m% *}: blocks, largest kg/s; tol 1e-15)"
else ko "cold: the IGLOO variant could not be made"; fi

echo "--- twin: test/igloo-burnout-box ---"
d=$T/twin; log=$T/twin.log; rm -rf "$d"
"$RUN" "$L" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
echo "  exit=$rc  $(grep -m1 '^# wall' "$log")"
[ $rc = 0 ] && grep -q "GAME OVER" "$log" && ok "twin: run exits 0 with GAME OVER" || ko "twin: run exit=$rc"
bad=$(grep -c -e '\[ERROR\]' -e ' FAIL' -e 'on-physical state' -e 'Integration failed at' -e 'p<0 detected' \
              -e 'ABORT:' -e 'stuck in cell' -e 'no net progress' -e 'hit outer maxIter' -e 'floor guard' "$log")
[ "$bad" = 0 ] && ok "twin: log clean (no error, failure, parcel-abort or floor-guard line)" || ko "twin: log has $bad error/abort line(s)"
ol=$(python3 -B "$H/test/igloo-burn-box/check.py" "$d" "$log" "$H/test/evap-box-mom" 2>&1); take "$ol" twin
ei=$(num "$oi" '(SUM Q dt)\/(SUM M dt)'); el=$(num "$ol" '(Q - residue)\/M')
if [ -n "$ei" ] && [ -n "$el" ] && close "$ei" "$el" $TOL_E; then ok "twin: energy per kg handed, ICE $ei against IGLOO $el J/kg: $(rel "$ei" "$el") (tol $TOL_E)"
else ko "twin: energy per kg handed, ICE '$ei' against IGLOO '$el' J/kg (tol $TOL_E)"; fi

echo "--- refuse: d-stop beside d-burnout ---"
for S in "$I" "$L"; do
  c=$(basename "$S"); v=$T/refuse-$c; log=$T/refuse-$c.log; rm -rf "$v"; mkdir -p "$v"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$v/"
  sed -i 's/ d-burnout=3e-06$/ d-stop=1e-06 d-burnout=3e-06/' "$v/INPUT/part-models.txt"
  grep -q ' d-stop=1e-06 d-burnout=3e-06$' "$v/INPUT/part-models.txt" || ko "refuse: the row of $c could not be edited"
  "$RUN" "$v" "$v/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  n=$(grep -c "d-stop beside d-burnout" "$log")
  if [ $rc != 0 ] && [ "$n" -ge 1 ] && ! grep -q "GAME OVER" "$log"; then ok "refuse: $c stops on d-stop beside d-burnout (exit=$rc)"
  else ko "refuse: $c did not stop on d-stop beside d-burnout (exit=$rc, lines=$n)"; fi
done

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
