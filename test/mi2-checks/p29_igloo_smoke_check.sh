#!/bin/bash
# IGLOO smoke gate (test/igloo-smoke-box): a burning IGLOO material whose reaction makes a condensed product, the smoke,
# born in an ICE family (the tail of its combustion record: the smoke material, kg per kg burnt, its h298, ICE's phase).
#  books  : the case as it is: the smoke of every MOSE block joins ICE after the block; the gas, the smoke ICE holds and
#           the streams close mass, momentum and energy; the smoke's temperature and number at birth (check.py K0-K7).
#  rates  : MOSE steady, ICE time-accurate (RK3, three steps): the smoke is handed as rates, a source of its family
#           over ICE's own steps (DECISION 21); what ICE holds = the rates times ICE's time (check.py R0-R5).
#  refuse : ICE not run; relaxation = true; the burning row has no d-smoke; the smoke is not a material of ICE; it is
#           born in another phase; the Eulerian coupling (IGLOO-iter = 2): each stops before the first step. (A smoke
#           family of another closure is refused too, but no case reaches it: ICE reads ATLAS ICB's initial condition,
#           which holds the six variables of MK, and stops on any other closure first.)
#   usage: p29_igloo_smoke_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P29}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/igloo-smoke-box
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
take() { pass=$((pass + $(echo "$1" | grep -c '^PASS'))); fail=$((fail + $(echo "$1" | grep -c '^FAIL'))); echo "$1" | grep -v '^== check.py'
         echo "$1" | grep -q '^== check.py' || ko "$2: check.py did not finish: $(echo "$1" | tail -1)"; }
echo "  exe: $EXE ($(md5sum < "$EXE" | cut -c1-8))"

grep -q '^AL2O3(a) 1$' "$S/INPUT/partE-phase.txt" && grep -q '^\[ICE-Family1\]$' "$S/input.ini" \
  && ok "the smoke AL2O3(a): a material of ICE's phase partE with one population" || ko "the case has no smoke family"
grep -q '^material A combustion 1 formation [^ ]* AL2O3(a) 5.00000[0-9]*e-01 [^ ]* partE-$' "$S/INPUT/partL-products.txt" \
  && ok "the record of A makes 0.5 kg of AL2O3(a) per kg burnt, born in partE-" || ko "the record of A has no smoke tail"
grep -q ' d-smoke=3e-06$' "$S/INPUT/partL-models.txt" && ok "the row of A gives d-smoke 3e-06 m" || ko "the row of A has no d-smoke"

echo "--- books ---"
d=$T/books; log=$T/books.log; rm -rf "$d"
"$RUN" "$S" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
echo "  exit=$rc  $(grep -m1 '^# wall' "$log")"
[ $rc = 0 ] && grep -q "GAME OVER" "$log" && ok "books: run exits 0 with GAME OVER" || ko "books: run exit=$rc: $(grep -a -m1 'ERROR' "$log" | cut -c1-200)"
o=$(python3 -B "$S/check.py" "$d" "$log" "$H/test/evap-box-mom" 2>&1); take "$o" books

echo "--- rates ---"
d=$T/rates; c=$T/rates-case; log=$T/rates.log; rm -rf "$d" "$c"; mkdir -p "$c"
cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$c/"
sed -i -e 's/^time-accurate         = true$/time-accurate         = false/' -e 's/^ICE-iter      = 0 1$/ICE-iter      = 0 3/' \
       -e 's/^time-scheme          = euler$/time-scheme          = RK3/' "$c/input.ini"
n=$(diff "$S/input.ini" "$c/input.ini" | grep -c '^>')
[ "$n" = 3 ] && ok "rates: the deck differs from the case in 3 lines (MOSE steady; ICE RK3, 3 steps)" \
  || ko "rates: the deck edit changed $n lines, not 3"
"$RUN" "$c" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
echo "  exit=$rc  $(grep -m1 '^# wall' "$log")"
[ $rc = 0 ] && grep -q "GAME OVER" "$log" && ok "rates: run exits 0 with GAME OVER" || ko "rates: run exit=$rc: $(grep -a -m1 'ERROR' "$log" | cut -c1-200)"
o=$(python3 -B "$S/check.py" "$d" "$log" "$H/test/evap-box-mom" 2>&1); take "$o" rates

echo "--- refuse ---"
stop() {  # stop <name> <message> <file of the case> <sed script>
  local d=$T/refuse-$1 log=$T/refuse-$1.log
  rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
  sed -i "$4" "$d/$3"
  if cmp -s "$d/$3" "$S/$3"; then ko "refuse $1: the edit changed nothing"; return; fi
  "$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  if [ $rc != 0 ] && grep -a -q -F "$2" "$log" && ! grep -q "GAME OVER" "$log"; then ok "refuse $1: stops ($2)"
  else ko "refuse $1: exit=$rc, message \"$2\" $(grep -a -q -F "$2" "$log" && echo printed || echo absent)"; fi
}
stop no-ice     "and the smoke has no family: ICE is not run"                     input.ini '/^ICE-iter /d'
stop relax      "relaxation = true blends two coupling steps of the gas's share"   input.ini 's/^relaxation    = false$/relaxation    = true/'
stop no-dsmoke  "gives no d-smoke, the diameter each smoke drop is born with"     INPUT/partL-models.txt 's/ d-smoke=3e-06$//'
stop not-ice    "burns into the smoke Q, which is not a material of ICE's phase"  INPUT/partL-products.txt \
                's/^\(material A combustion 1 formation [^ ]*\) AL2O3(a) /\1 Q /'
stop phase      "is born in the phase drop-, and ICE runs the phase partE-"       INPUT/partL-products.txt \
                's/^\(material A combustion .*\) partE-$/\1 drop-/'
stop eulerian   "a cycle has IGLOO-iter = 2"                                       input.ini 's/^IGLOO-iter    = 1 0$/IGLOO-iter    = 2 0/'

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
