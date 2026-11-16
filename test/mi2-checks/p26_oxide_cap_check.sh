#!/bin/bash
# oxide-cap gate (test/igloo-cap-box): an IGLOO material that burns with an oxide cap. IGLOO hands the coupler, beside
# the net mass of the material, the oxide its caps gained and the cap oxide handed to the gas; the coupler splits the
# gas species between the three records of INPUT/part-products.txt (combustion, cap, release).
#  on     : the case as it is (residue = on: the cap stays, a residue of the oxide). check.py: the net mass, each gas
#           species, the energy and the residues against closed forms; the coupler's books.
#  off    : the same case with residue = off: at burnout the cap is handed to the gas through the release record, so
#           everything a stream carries in ends in the gas, as if no cap had formed: the closed forms of
#           test/igloo-burn-box/check.py.
#  refuse : the products file without the cap record; residue = off without the release record; a cap record whose
#           oxide ratio is not the one IGLOO's cap law holds: each stops before the first step and names the record.
#   usage: p26_oxide_cap_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P26}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/igloo-cap-box
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
take() { pass=$((pass + $(echo "$1" | grep -c '^PASS'))); fail=$((fail + $(echo "$1" | grep -c '^FAIL'))); echo "$1" | grep -v '^== check.py'
         echo "$1" | grep -q '^== check.py' || ko "$2: check.py did not finish: $(echo "$1" | tail -1)"; }
clean() { grep -c -e '\[ERROR\]' -e ' FAIL' -e 'on-physical state' -e 'Integration failed at' -e 'p<0 detected' \
                  -e 'ABORT:' -e 'stuck in cell' -e 'no net progress' -e 'hit outer maxIter' -e 'floor guard' "$1"; }
echo "  exe: $EXE ($(md5sum < "$EXE" | cut -c1-8))"

grep -q '^A 1 combustion=Beckstead cap=constant$' "$S/INPUT/part-phase.txt" && grep -q '^OX 0$' "$S/INPUT/part-phase.txt" \
  && ok "the phase file: A burns with a cap, OX is a material of 0 groups" || ko "the phase file does not hold the capped material and its oxide"
grep -q '^A cap constant xi-cap=0.2 product=OX residue=on$' "$S/INPUT/part-models.txt" && grep -q '^A combustion Beckstead ' "$S/INPUT/part-models.txt" \
  && ok "the models file: the combustion row of A and its cap row (xi-cap 0.2, product OX, residue on)" || ko "the models file does not hold the two rows of A"
n=$(grep -c -e '^material A combustion ' -e '^material A cap ' -e '^material OX release ' "$S/INPUT/part-products.txt")
[ "$n" = 3 ] && ok "the products file holds the three records (combustion, cap, release)" || ko "the products file holds $n of the three records"

echo "--- on: the residue stays ---"
d=$T/on; log=$T/on.log; rm -rf "$d"
"$RUN" "$S" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
echo "  exit=$rc  $(grep -m1 '^# wall' "$log")"
[ $rc = 0 ] && grep -q "GAME OVER" "$log" && ok "on: run exits 0 with GAME OVER" || ko "on: run exit=$rc: $(grep -m1 'ERROR' "$log" | cut -c1-160)"
bad=$(clean "$log"); [ "$bad" = 0 ] && ok "on: log clean (no error, failure, parcel-abort or floor-guard line)" || ko "on: log has $bad error/abort line(s)"
grep -q "its oxide cap takes OX" "$log" && grep -q "the cap stays condensed at burnout (residue = on)" "$log" \
  && ok "on: the coupler reports the cap record and that the cap stays" || ko "on: the coupler does not report the cap"
o=$(python3 -B "$S/check.py" "$d" "$log" "$H/test/evap-box-mom" 2>&1); take "$o" on

echo "--- off: the cap is handed to the gas ---"
v=$T/off-case; d=$T/off; log=$T/off.log; rm -rf "$v" "$d"; mkdir -p "$v"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$v/"
sed -i 's/^\(A cap constant .*\) residue=on$/\1 residue=off/' "$v/INPUT/part-models.txt"
sed -i 's/^residue    = on$/residue    = off/' "$v/input.ini"
if grep -q ' residue=off$' "$v/INPUT/part-models.txt" && grep -q '^residue    = off$' "$v/input.ini"; then
  "$RUN" "$v" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  [ $rc = 0 ] && grep -q "GAME OVER" "$log" && ok "off: run exits 0 with GAME OVER" || ko "off: run exit=$rc: $(grep -m1 'ERROR' "$log" | cut -c1-160)"
  bad=$(clean "$log"); [ "$bad" = 0 ] && ok "off: log clean" || ko "off: log has $bad error/abort line(s)"
  grep -q "the cap is handed to the gas at burnout (residue = off)" "$log" && ok "off: the coupler reports the release record" \
                                                                           || ko "off: the coupler does not report the release record"
  o=$(python3 -B "$H/test/igloo-burn-box/check.py" "$d" "$log" "$H/test/evap-box-mom" 2>&1); take "$o" off
else ko "off: the variant could not be made"; fi

echo "--- refuse ---"
refuse() {  # refuse <name> <message> : the variant in $T/refuse-<name> must stop and print <message>
  local v=$T/refuse-$1 log=$T/refuse-$1.log
  "$RUN" "$v" "$v/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  n=$(grep -c "$2" "$log")
  if [ $rc != 0 ] && [ "$n" -ge 1 ] && ! grep -q "GAME OVER" "$log"; then ok "refuse: $1 stops and names it (exit=$rc)"
  else ko "refuse: $1 was not refused (exit=$rc, lines=$n)"; fi
}
v=$T/refuse-no-cap-record; rm -rf "$v"; mkdir -p "$v"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$v/"
sed -i -e '/^material A cap /,+1d' "$v/INPUT/part-products.txt"
grep -q '^material A cap ' "$v/INPUT/part-products.txt" && ko "could not remove the cap record" || refuse no-cap-record "has no cap record for it"
v=$T/refuse-no-release-record; rm -rf "$v"; mkdir -p "$v"; cp -rp "$T/off-case/INPUT" "$S/MESH" "$T/off-case/input.ini" "$v/"
sed -i -e '/^material OX release /,+2d' "$v/INPUT/part-products.txt"; echo "material OX none 0 none 0.0000000000e+00" >> "$v/INPUT/part-products.txt"
grep -q '^material OX release ' "$v/INPUT/part-products.txt" && ko "could not remove the release record" || refuse no-release-record "has no release record for OX"
v=$T/refuse-oxide-ratio; rm -rf "$v"; mkdir -p "$v"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$v/"
sed -i -e 's/^\(material A cap 1 formation [^ ]*\) OX 1.8894415674e+00 /\1 OX 1.5000000000e+00 /' \
       -e 's/^  O2 -8.8944156736e-01 /  O2 -5.0000000000e-01 /' "$v/INPUT/part-products.txt"
grep -q ' OX 1.5000000000e+00 ' "$v/INPUT/part-products.txt" && refuse oxide-ratio "the cap law of IGLOO holds" || ko "could not change the oxide ratio"

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
