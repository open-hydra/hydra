#!/bin/bash
# smoke gate (test/ice-smoke-box): a burning ICE material whose reaction makes a condensed product, the smoke, born in
# an ICE family of the MK closure (the tail of its combustion record: the smoke material, kg per kg burnt, its h298,
# ICE's phase).
#  books  : the case as it is: the cloud burns out, ICE alone; the gas's share the coupler hands, the smoke ICE holds
#           at the end and the cloud's start close mass, momentum and energy (check.py B1-B5).
#  birth  : one ICE step: the smoke against its birth (mass nu_s M_b dt, the drop's velocity and temperature, each
#           drop at d-smoke).
#  gas    : one ICE step, its sources frozen through 5 MOSE iterations: the gas gains the record's species without the
#           smoke, (1 - nu_s) M_b t, and the reported momentum and energy.
#  refuse : the smoke is not a material of ICE; it is born in another phase; the burning row has no d-smoke; an
#           IGLOO material's record makes a smoke (C-IGLOO9 not built): each stops before the first step.
#   usage: p28_smoke_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P28}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/ice-smoke-box
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
take() { pass=$((pass + $(echo "$1" | grep -c '^PASS'))); fail=$((fail + $(echo "$1" | grep -c '^FAIL'))); echo "$1" | grep -v '^== check.py'
         echo "$1" | grep -q '^== check.py' || ko "$2: check.py did not finish: $(echo "$1" | tail -1)"; }
echo "  exe: $EXE ($(md5sum < "$EXE" | cut -c1-8))"

grep -q '^S 1$' "$S/INPUT/part-phase.txt" && grep -q '^\[ICE-Family2\]$' "$S/input.ini" \
  && ok "the smoke S: a material of ICE's phase with one population, family 2" || ko "the case has no smoke family"
grep -q '^material A combustion 1 formation [^ ]* S 5.0000000000e-01 [^ ]* part-$' "$S/INPUT/part-products.txt" \
  && ok "the record of A makes 0.5 kg of S per kg burnt, born in part-" || ko "the record of A has no smoke tail"
grep -q ' d-smoke=3e-06$' "$S/INPUT/part-models.txt" && ok "the row of A gives d-smoke 3e-06 m" || ko "the row of A has no d-smoke"

echo "--- books: the burn to the end ---"
d=$T/books; log=$T/books.log; rm -rf "$d"
"$RUN" "$S" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
echo "  exit=$rc  $(grep -m1 '^# wall' "$log")"
[ $rc = 0 ] && grep -q "GAME OVER" "$log" && ok "books: run exits 0 with GAME OVER" || ko "books: run exit=$rc: $(grep -a -m1 'ERROR' "$log" | cut -c1-160)"
grep -q "ICE material A burns into the smoke S (family 2)" "$log" && ok "books: the coupler reports the smoke and its family" \
                                                                 || ko "books: the coupler does not report the smoke"
o=$(python3 -B "$S/check.py" "$d" "$log" "$H/test/ice-box" 2>&1); take "$o" books

variant() {  # variant <name> <sed script for input.ini>...: a copy of the case in $T/<name>-case, run into $T/<name>
  local name=$1; shift
  local v=$T/$name-case; rm -rf "$v" "$T/$name"; mkdir -p "$v"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$v/"
  for e in "$@"; do sed -i "$e" "$v/input.ini"; done
}
echo "--- birth: one ICE step ---"
variant birth 's/^N-coupling    = 600 1$/N-coupling    = 1 1/'
if grep -q '^N-coupling    = 1 1$' "$T/birth-case/input.ini"; then
  "$RUN" "$T/birth-case" "$T/birth" 1 "$T/birth.log" "$EXE" > /dev/null 2>&1; rc=$?
  [ $rc = 0 ] && grep -q "GAME OVER" "$T/birth.log" && ok "birth: run exits 0 with GAME OVER" || ko "birth: run exit=$rc"
  o=$(python3 -B "$S/check.py" "$T/birth" "$T/birth.log" "$H/test/ice-box" birth 2>&1); take "$o" birth
else ko "birth: the variant could not be made"; fi

echo "--- gas: one ICE step, frozen through 5 MOSE iterations ---"
variant gas 's/^N-coupling    = 600 1$/N-coupling    = 1 5/' 's/^sol-diter     = 200$/sol-diter     = 5/'
if grep -q '^N-coupling    = 1 5$' "$T/gas-case/input.ini" && grep -q '^sol-diter     = 5$' "$T/gas-case/input.ini"; then
  "$RUN" "$T/gas-case" "$T/gas" 1 "$T/gas.log" "$EXE" > /dev/null 2>&1; rc=$?
  [ $rc = 0 ] && grep -q "GAME OVER" "$T/gas.log" && ok "gas: run exits 0 with GAME OVER" || ko "gas: run exit=$rc"
  o=$(python3 -B "$S/check.py" "$T/gas" "$T/gas.log" "$H/test/ice-box" gas 2>&1); take "$o" gas
else ko "gas: the variant could not be made"; fi

echo "--- refuse ---"
stop() {  # stop <name> <message> <case dir> <file of the case> <sed script>
  local d=$T/refuse-$1 log=$T/refuse-$1.log
  rm -rf "$d"; mkdir -p "$d"; cp -rp "$3/INPUT" "$3/MESH" "$3/input.ini" "$d/"
  sed -i "$5" "$d/$4"
  if cmp -s "$d/$4" "$3/$4"; then ko "refuse $1: the edit changed nothing"; return; fi
  "$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  if [ $rc != 0 ] && grep -a -q -F "$2" "$log" && ! grep -q "GAME OVER" "$log"; then ok "refuse $1: stops ($2)"
  else ko "refuse $1: exit=$rc, message \"$2\" $(grep -a -q -F "$2" "$log" && echo printed || echo absent)"; fi
}
stop not-ice   "burns into the smoke Q, which is not a material of ICE's phase" "$S" INPUT/part-products.txt \
               's/^\(material A combustion 1 formation [^ ]*\) S /\1 Q /'
stop phase     "is born in the phase drop-, and ICE runs the phase part-"       "$S" INPUT/part-products.txt \
               's/^\(material A combustion .*\) part-$/\1 drop-/'
stop no-dsmoke "gives no d-smoke, the diameter each smoke drop is born with"    "$S" INPUT/part-models.txt ' s/ d-smoke=3e-06$//'
L=$H/test/igloo-burnout-box
stop igloo     "IGLOO does not yet hand the coupler what each parcel burns"     "$L" INPUT/part-products.txt \
               's/^\(material A combustion 1 formation [^ ]*\)$/\1 S 5.0000000000e-01 -1.0000000000e+06 part-/; s/^  N2 1.0000000000e+00 1.0000000000e+00$/  N2 5.0000000000e-01 5.0000000000e-01/'

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
