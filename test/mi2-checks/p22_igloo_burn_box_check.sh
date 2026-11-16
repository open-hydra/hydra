#!/bin/bash
# igloo-burn-box gate: four cases of IGLOO drops burning in a closed box. test/igloo-burn-box (combustion =
# Beckstead): test/evap-box-mom with a gas of one species, the burnt mass handed as that species (a-n2.yaml), its heat
# of reaction in the enthalpy h0 of the material; test/igloo-burn-box-rx: graphite burning to CO2 in air (c-o2.yaml);
# test/igloo-burn-box-beckstead and test/igloo-burn-box-direct: the twins of the two ICE burning cases (the same
# burning section, down to a residue d-stop), one per model. Every stream burns out, or down to its residue, inside
# the box, so the gas must gain what the streams carry in less what the residue keeps: test/igloo-burn-box/check.py
# compares the reported mass, force and energy rates with the closed forms of the README.md files, the gas energy
# rise of every block with the reported rate, the coupler's closed-domain books and the warming of the gas (a gas of
# several species: against the heat of reaction, and the species; a gas of one species: against the first law of the
# rigid box).
#  default: test/igloo-burn-box with its model given by [IGLOO-Models] combustion instead of the material line must
#           write the same bytes (its row of the models file names the model the section names).
#  line   : the same case with a constant of its burning model left on the material line must stop and name the
#           models file.
#  heat   : the coupler prints, for each burning material, the record's heat from the material at the temperature
#           it burns at and q-comb beside it; it warns when the two differ by more than 5 % for Beckstead (the row's
#           q-comb at 930000: 4.55 %, no warning; at 1030000: 5.71 %, one), never for direct (2000000: none).
#  column : test/igloo-burn-box-direct with its latent heat as the Lv column of the property table (2.26e6 J/kg at
#           the burning temperature, 1000 J/kg more per kelvin) and no Lv key must write the same bytes; with the key
#           beside the column, and with neither, it must stop and say so.
#  lv     : test/igloo-burn-box with IGLOO's phase change on ([IGLOO-Models] evaporation) and a latent heat listed for
#           its burning material ([IGLOO-Properties] Lv) must write the same bytes: the datum offset of a burning
#           material holds no latent heat (the pin subtracted it).
#  refuse : the same case whose products file has no combustion record for the material must stop before the first
#           step and name the record.
#   usage: p22_igloo_burn_box_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P22}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/igloo-burn-box
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
echo "  exe: $EXE ($(md5sum < "$EXE" | cut -c1-8))"

for c in igloo-burn-box igloo-burn-box-rx igloo-burn-box-beckstead igloo-burn-box-direct; do
  d=$T/$c; log=$T/$c.log; rm -rf "$d"
  "$RUN" "$H/test/$c" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  echo "  $c: exit=$rc  $(grep -m1 '^# wall' "$log")"
  [ $rc = 0 ] && ok "$c: run exits 0" || ko "$c: run exit=$rc"
  grep -q "GAME OVER" "$log" && ok "$c: GAME OVER printed" || ko "$c: no GAME OVER"
  bad=$(grep -c -e '\[ERROR\]' -e ' FAIL' -e 'on-physical state' -e 'Integration failed at' -e 'p<0 detected' \
                -e 'ABORT:' -e 'stuck in cell' -e 'no net progress' -e 'hit outer maxIter' -e 'floor guard' "$log")
  [ "$bad" = 0 ] && ok "$c: log clean (no error, failure, parcel-abort or floor-guard line)" || ko "$c: log has $bad error/abort line(s)"
  grep -q -- "\[combustion\] :" "$log" && ok "$c: the coupler reports the combustion record" \
                                        || ko "$c: no combustion record reported"
  out=$(python3 -B "$S/check.py" "$d" "$log" "$H/test/evap-box-mom" 2>&1); echo "$out" | grep -v '^== check.py'
  pass=$((pass + $(echo "$out" | grep -c '^PASS'))); fail=$((fail + $(echo "$out" | grep -c '^FAIL')))
  echo "$out" | grep -q '^== check.py' || ko "$c: check.py did not finish: $(echo "$out" | tail -1)"
done

# --- lv: phase change on and a latent heat listed for the burning material ---
d=$T/lv; log=$T/lv.log; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
sed -i 's/^heat = Ranz-Marshall$/heat = Ranz-Marshall\nevaporation = CEM\n[IGLOO-Properties]\nTboil = 373.15\nLe = 1.0\nYinf = 0.0\npsat = 3.0e3\nMv = 18.015\nLv = 2.26e6\ncpv = 1900/' "$d/input.ini"
if grep -q '^Lv = 2.26e6$' "$d/input.ini" && grep -q '^evaporation = CEM$' "$d/input.ini"; then
  "$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  n=$(ls "$T/igloo-burn-box/OUTPUT" | wc -l)
  grep -q 'DISABLING evaporation for this material' "$log" && ok "lv: IGLOO has phase change on and drops the evaporation of the burning material" \
                                                            || ko "lv: IGLOO did not see the phase change"
  if [ $rc = 0 ] && [ "$n" -gt 0 ] && diff -rq "$T/igloo-burn-box/OUTPUT" "$d/run/OUTPUT" > /dev/null; then
    ok "lv: a latent heat listed for the burning material changes none of the $n OUTPUT files"
  else ko "lv: the run with a latent heat listed differs (exit=$rc): $(diff -rq "$T/igloo-burn-box/OUTPUT" "$d/run/OUTPUT" 2>&1 | head -1)"; fi
else ko "lv: the variant could not be made"; fi

# --- default: the model given by [IGLOO-Models], no token on the material line ---
d=$T/default; log=$T/default.log; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
sed -i 's/^A 1 combustion=Beckstead$/A 1/' "$d/INPUT/part-phase.txt"
sed -i 's/^heat = Ranz-Marshall$/heat = Ranz-Marshall\ncombustion = Beckstead/' "$d/input.ini"
if grep -q '^A 1$' "$d/INPUT/part-phase.txt" && grep -q '^A combustion Beckstead K-burn=' "$d/INPUT/part-models.txt" \
   && grep -q '^combustion = Beckstead$' "$d/input.ini"; then
  "$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  n=$(ls "$T/igloo-burn-box/OUTPUT" | wc -l)
  if [ $rc = 0 ] && [ "$n" -gt 0 ] && diff -rq "$T/igloo-burn-box/OUTPUT" "$d/run/OUTPUT" > /dev/null; then
    ok "default: the model given by [IGLOO-Models] writes the same $n OUTPUT files"
  else ko "default: the model given by [IGLOO-Models] differs (exit=$rc): $(diff -rq "$T/igloo-burn-box/OUTPUT" "$d/run/OUTPUT" 2>&1 | head -1)"; fi
else ko "default: could not move the model to [IGLOO-Models]"; fi

# --- refuse: no combustion record for the burning material ---
d=$T/refuse; log=$T/refuse.log; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
sed -i -e 's/^material A combustion 1 formation .*/material A none 0 none 0.0000000000e+00/' -e '/^  N2 /d' "$d/INPUT/part-products.txt"
grep -q '^material A none 0 none' "$d/INPUT/part-products.txt" || ko "could not remove the combustion record"
"$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
n=$(grep -c "has no combustion record for it" "$log")
echo "  refuse: exit=$rc, refusal lines=$n"
if [ $rc != 0 ] && [ "$n" -ge 1 ] && ! grep -q "GAME OVER" "$log"; then ok "a burning material without a combustion record is refused (exit=$rc)"
else ko "the missing record was not refused (exit=$rc, refusal lines=$n)"; fi

# --- line: a constant of the burning model left on the material line ---
d=$T/line; log=$T/line.log; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
sed -i 's/^A 1 combustion=Beckstead$/A 1 combustion=Beckstead K-burn=1.39e-06/' "$d/INPUT/part-phase.txt"
"$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
if grep -q ' K-burn=' "$d/INPUT/part-phase.txt" && [ $rc != 0 ] && ! grep -q "GAME OVER" "$log" \
   && grep -q "the constants of the burning models are rows of the models file" "$log"; then
  ok "line: a constant of the burning model on the material line is refused (exit=$rc)"
else ko "line: a constant on the material line was not refused (exit=$rc)"; fi

# --- column: the latent heat of combustion = direct from the Lv column of the property table ---
C=$H/test/igloo-burn-box-direct
column() {  # column <case dir>: an Lv column, 2.26e6 J/kg at 300 K and 1000 J/kg more per kelvin, appended to the table
  awk '/^VARIABLES/ { print $0 ", \"Lv\""; next }
       /^ *[0-9.+-]+([eE][0-9+-]+)? / && NF > 2 { printf "%s %.10e\n", $0, 2.26e6 + 1000.0*($1 - 300.0); next }
       { print }' "$1/INPUT/part-properties.dat" > "$1/INPUT/part-properties.new" && mv "$1/INPUT/part-properties.new" "$1/INPUT/part-properties.dat"
}
d=$T/column; log=$T/column.log; rm -rf "$d"; mkdir -p "$d"; cp -rp "$C/INPUT" "$C/MESH" "$C/input.ini" "$d/"
column "$d"; sed -i '/^Lv = 2.26e6$/d' "$d/input.ini"
if grep -q '^VARIABLES.*"Lv"$' "$d/INPUT/part-properties.dat" && ! grep -q '^Lv = ' "$d/input.ini" \
   && [ "$(awk '$1 == "300.0" || $1 == "300" { print $NF }' "$d/INPUT/part-properties.dat")" = 2.2600000000e+06 ]; then
  "$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  n=$(ls "$T/igloo-burn-box-direct/OUTPUT" | wc -l)
  grep -q "combustion=direct: Lv( 3.00000E+02 K) =  2.26000E+06 J/kg" "$log" && ok "column: IGLOO reads the Lv column at the burning temperature" \
                                                                           || ko "column: IGLOO does not report Lv(300 K) = 2.26e6 J/kg"
  if [ $rc = 0 ] && [ "$n" -gt 0 ] && diff -rq "$T/igloo-burn-box-direct/OUTPUT" "$d/run/OUTPUT" > /dev/null; then
    ok "column: the latent heat from the column writes the same $n OUTPUT files"
  else ko "column: the run with the Lv column differs (exit=$rc): $(diff -rq "$T/igloo-burn-box-direct/OUTPUT" "$d/run/OUTPUT" 2>&1 | head -1)"; fi
else ko "column: the variant could not be made"; fi
d=$T/column-both; log=$T/column-both.log; rm -rf "$d"; mkdir -p "$d"; cp -rp "$C/INPUT" "$C/MESH" "$C/input.ini" "$d/"; column "$d"
"$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
if [ $rc != 0 ] && grep -q "Lv given both in \[IGLOO-Properties\] and as a properties.dat column" "$log" && ! grep -q "GAME OVER" "$log"; then
  ok "column: Lv beside the Lv column is refused"
else ko "column: Lv beside the Lv column was not refused (exit=$rc)"; fi
d=$T/column-none; log=$T/column-none.log; rm -rf "$d"; mkdir -p "$d"; cp -rp "$C/INPUT" "$C/MESH" "$C/input.ini" "$d/"
sed -i -e '/^\[IGLOO-Properties\]$/d' -e '/^Lv = 2.26e6$/d' "$d/input.ini"
"$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
if ! grep -q '^Lv = ' "$d/input.ini" && [ $rc != 0 ] && grep -q "combustion=direct requires its latent heat" "$log" && ! grep -q "GAME OVER" "$log"; then
  ok "column: neither the column nor the key is refused"
else ko "column: a direct material without a latent heat was not refused (exit=$rc)"; fi

# --- heat: q-comb beside the heat the gas receives from the material where it burns ---
heat() {  # heat <tag> <case> <q-comb of the row> <its printed form> <warnings expected>
  local tag=$1 c=$2 q=$3 d=$T/heat-$1 log=$T/heat-$1.log n w rc
  rm -rf "$d"; mkdir -p "$d"; cp -rp "$H/test/$c/INPUT" "$H/test/$c/MESH" "$H/test/$c/input.ini" "$d/"
  sed -i "s/ q-comb=1000000\( \|\$\)/ q-comb=$q\1/" "$d/INPUT/part-models.txt"
  if grep -q " q-comb=$q\( \|\$\)" "$d/INPUT/part-models.txt"; then
    "$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
    n=$(grep -a -c 'the same from the material at .* K, where it burns' "$log"); w=$(grep -a -c '\[WARNING\] q-comb = ' "$log")
    if [ $rc = 0 ] && [ "$n" = 1 ] && grep -a -q "q-comb of the solver  *:  *$4\$" "$log" && [ "$w" = "$5" ]; then
      ok "heat $tag: $c with q-comb = $q: the heat from the material where it burns is printed beside q-comb, $w warning(s)"
    else ko "heat $tag: $c with q-comb = $q: exit=$rc, $n heat line(s), $w warning(s) where $5 expected"; fi
  else ko "heat $tag: the row of $c could not be edited"; fi
}
for c in igloo-burn-box-beckstead igloo-burn-box-direct; do
  n=$(grep -a -c 'the same from the material at .* K, where it burns' "$T/$c.log"); w=$(grep -a -c '\[WARNING\] q-comb = ' "$T/$c.log")
  [ "$n" = 1 ] && [ "$w" = 0 ] && ok "heat: $c as shipped prints the heat where it burns and no warning" \
                               || ko "heat: $c as shipped: $n heat line(s), $w warning(s)"
done
heat near   igloo-burn-box-beckstead   930000  9.3000000E+05 0
heat far    igloo-burn-box-beckstead   1030000 1.0300000E+06 1
heat direct igloo-burn-box-direct 2000000 2.0000000E+06 0
grep -a -q 'combustion = direct: q-comb sets the rate' "$T/heat-direct.log" && ! grep -a -q 'combustion = direct: q-comb sets the rate' "$T/heat-near.log" \
  && ok "heat direct: the coupler says what q-comb is for direct, and only for direct" || ko "heat direct: the line on q-comb of direct is missing, or printed for Beckstead"

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
