#!/bin/bash
# ice-burn-box gate (test/ice-burn-box, test/ice-burn-box-beckstead): test/ice-box with its material burning, once per
# model (combustion = direct, combustion = Beckstead), the model on the material line of the phase file, its constants
# in the material's row of the models file. The gas must
# gain the mass, the momentum and the energy of the closed forms of test/ice-burn-box/README.md: check.py compares
# the reported rates and the gas's gains with them, the coupler's closed-domain books, the non-vacuity of the
# mass and of the heat of reaction, and the warming of the gas with the first law of the rigid box.
#  default: the first case with its model and constants given by [ICE-Physics] instead of the material line and the
#           models file must write the same bytes.
#  refuse : the first case whose products file has no combustion record for the material must stop before the first
#           step and name the record; so must a record with the former word "explicit", the former key d-res, a
#           constant left on the material line, the oxide cap (its word on the line, its row), a row of an unknown
#           kind, of another model than the line's, with a key the model does not read.
#  heat   : the coupler prints, for each burning material, the record's heat from the material at the temperature
#           it burns at and q-comb beside it; it warns when the two differ by more than 5 % for Beckstead (the row's
#           q-comb at 930000: 4.55 %, no warning; at 1030000: 5.71 %, one), never for direct (2000000: none).
#  lv     : the first case with its latent heat as the Lv column of the property table (2.26e6 J/kg at the burning
#           temperature, 1000 J/kg more per kelvin) and no latent-heat key must write the same bytes; with the key
#           beside the column, and with neither, it must stop and say so.
#   usage: p21_ice_burn_box_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P21}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
BAD=(-e '\[ERROR\]' -e ' FAIL' -e 'on-physical state' -e 'Integration failed at' -e 'p<0 detected' -e 'ABORT:')
echo "  exe: $EXE ($(md5sum < "$EXE" | cut -c1-8))"

for C in ice-burn-box ice-burn-box-beckstead; do
  S=$H/test/$C
  word=direct; [ $C = ice-burn-box-beckstead ] && word=Beckstead
  echo "--- $C (combustion = $word) ---"
  diff -rq -x part-phase.txt -x part-models.txt -x part-products.txt -x part-properties.dat "$H/test/ice-box/INPUT" "$S/INPUT" > /dev/null \
    && diff -rq "$H/test/ice-box/MESH" "$S/MESH" > /dev/null \
    && ok "$C: INPUT (but the phase, models, products and property files) and MESH are those of test/ice-box" || ko "$C: INPUT or MESH differ from test/ice-box"
  grep -q "^A 1 combustion=$word\$" "$S/INPUT/part-phase.txt" && ok "$C: the material line names combustion=$word and nothing else" \
                                                              || ko "$C: the material line is not 'A 1 combustion=$word'"
  [ "$(grep -c "^A combustion $word .*q-comb=" "$S/INPUT/part-models.txt")" = 1 ] && [ "$(grep -v -e '^#' -e '^$' -c "$S/INPUT/part-models.txt")" = 1 ] \
    && ok "$C: the models file holds one row, the constants of combustion=$word" || ko "$C: the models file does not hold the one row of A"
  d=$T/$C; log=$T/$C.log; rm -rf "$d"
  "$RUN" "$S" "$d" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  echo "  exit=$rc  $(grep -m1 '^# wall' "$log")"
  [ $rc = 0 ] && ok "$C: run exits 0" || ko "$C: run exit=$rc"
  grep -q "GAME OVER" "$log" && ok "$C: GAME OVER printed" || ko "$C: no GAME OVER"
  bad=$(grep -c "${BAD[@]}" "$log")
  if [ "$bad" = 0 ]; then ok "$C: log clean"; else ko "$C: log has $bad error line(s): $(grep -m1 "${BAD[@]}" "$log")"; fi
  grep -q "combustion = $word" "$log" && ok "$C: ICE reports the burning model" || ko "$C: ICE does not report combustion = $word"
  grep -q -- "- A \[combustion\] :" "$log" && ok "$C: the coupler reports the combustion record of A" \
                                           || ko "$C: no combustion record reported for A"
  out=$(python3 -B "$H/test/ice-burn-box/check.py" "$d" "$log" "$H/test/ice-box" 2>&1); echo "$out" | grep -v '^== check.py'
  pass=$((pass + $(echo "$out" | grep -c '^PASS'))); fail=$((fail + $(echo "$out" | grep -c '^FAIL')))
  echo "$out" | grep -q '^== check.py' || ko "$C: check.py did not finish: $(echo "$out" | tail -1)"
done

S=$H/test/ice-burn-box
# --- default: the model given by [ICE-Physics], no token on the material line; its constants stay in their row, whose
#     model word "-" stands for the model the material burns on ---
d=$T/default; log=$T/default.log; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
sed -i 's/^A 1 combustion=direct$/A 1/' "$d/INPUT/part-phase.txt"; sed -i 's/^A combustion direct /A combustion - /' "$d/INPUT/part-models.txt"
sed -i 's/^\[ICE-Physics\]$/[ICE-Physics]\ncombustion = direct/' "$d/input.ini"
if grep -q '^A 1$' "$d/INPUT/part-phase.txt" && grep -q '^combustion = direct$' "$d/input.ini" && grep -q '^A combustion - ' "$d/INPUT/part-models.txt"; then
  "$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  n=$(ls "$T/ice-burn-box/OUTPUT" | wc -l)
  if [ $rc = 0 ] && [ "$n" -gt 0 ] && diff -rq "$T/ice-burn-box/OUTPUT" "$d/run/OUTPUT" > /dev/null; then
    ok "default: the model given by [ICE-Physics] writes the same $n OUTPUT files"
  else ko "default: the model given by [ICE-Physics] differs (exit=$rc)"; fi
else ko "default: could not move the model to [ICE-Physics]"; fi

# --- refuse: no combustion record for the burning material ---
d=$T/refuse; log=$T/refuse.log; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
sed -i -e 's/^material A combustion 1 formation .*/material A none 0 none 0.0000000000e+00/' -e '/^  N2 /d' "$d/INPUT/part-products.txt"
grep -q '^material A none 0 none' "$d/INPUT/part-products.txt" || ko "could not remove the combustion record"
"$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
n=$(grep -c "has no combustion record for it" "$log")
echo "  refuse: exit=$rc, refusal lines=$n"
if [ $rc != 0 ] && [ "$n" -ge 1 ] && ! grep -q "GAME OVER" "$log"; then ok "a burning material without a combustion record is refused (exit=$rc)"
else ko "the missing record was not refused (exit=$rc, refusal lines=$n)"; fi

# --- refuse: the former word of the record, the former key of the residue, the reader of the models file ---
stop() {  # stop <name> <message> <file of the case> <sed script>
  local d=$T/$1 log=$T/$1.log
  rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
  sed -i "$4" "$d/$3"
  if cmp -s "$d/$3" "$S/$3"; then ko "$1: the edit changed nothing"; return; fi
  "$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  if [ $rc != 0 ] && grep -q -F "$2" "$log" && ! grep -q "GAME OVER" "$log"; then ok "$1: refused ($2)"
  else ko "$1: exit=$rc, message \"$2\" $(grep -q -F "$2" "$log" && echo printed || echo absent)"; fi
}
stop old-word "a combustion record that is not \"formation\"" INPUT/part-products.txt 's/ combustion 1 formation / combustion 1 explicit /'
stop old-key  "d-res is now d-stop"                           INPUT/part-models.txt   's/ d-stop=1e-06$/ d-res=1e-06/'
stop old-line "the constants of the burning models are rows of INPUT/part-models.txt" \
                                                              INPUT/part-phase.txt    's/^A 1 combustion=direct$/A 1 combustion=direct q-comb=1000000/'
stop cap-word "names the oxide cap of a burning drop"         INPUT/part-phase.txt    's/^A 1 combustion=direct$/A 1 combustion=direct cap=constant/'
stop cap-row  "a cap row gives the oxide cap of a burning drop" INPUT/part-models.txt '$a A cap constant xi-cap=0.2 product=OX'
stop kind     'unknown kind "breakup"'                        INPUT/part-models.txt   '$a A breakup TAB Comega=8'
stop model    "the material burns with combustion=direct"     INPUT/part-models.txt   's/^A combustion direct /A combustion Beckstead /'
stop key      "K-burn is not a constant of combustion=direct" INPUT/part-models.txt   's/^A combustion direct /A combustion direct K-burn=1.39e-06 /'
stop section  "[ICE-Physics] q-comb: a constant of the burning model, in the material's row of the models file" \
                                                              input.ini               's/^\[ICE-Physics\]$/[ICE-Physics]\nq-comb = 1.0e6/'
stop no-qcomb "requires q-comb, the heat of reaction per kg burnt, in its row of the models file" \
                                                              INPUT/part-models.txt   's/ q-comb=1000000 / /'

# --- lv: the latent heat from the Lv column of the property table ---
column() {  # column <case dir>: an Lv column, 2.26e6 J/kg at 300 K and 1000 J/kg more per kelvin, appended to the table
  awk '/^VARIABLES/ { print $0 ", \"Lv\""; next }
       /^ *[0-9.+-]+([eE][0-9+-]+)? / && NF > 2 { printf "%s %.10e\n", $0, 2.26e6 + 1000.0*($1 - 300.0); next }
       { print }' "$1/INPUT/part-properties.dat" > "$1/INPUT/part-properties.new" && mv "$1/INPUT/part-properties.new" "$1/INPUT/part-properties.dat"
}
d=$T/lv; log=$T/lv.log; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
column "$d"; sed -i '/^latent-heat   = 2.26e6$/d' "$d/input.ini"
if grep -q '^VARIABLES.*"Lv"$' "$d/INPUT/part-properties.dat" && ! grep -q '^latent-heat' "$d/input.ini" \
   && [ "$(awk '$1 == "300.0" || $1 == "300" { print $NF }' "$d/INPUT/part-properties.dat")" = 2.2600000000e+06 ]; then
  "$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
  n=$(ls "$T/ice-burn-box/OUTPUT" | wc -l)
  grep -q "latent heat of combustion = direct from the table (Lv column)" "$log" && ok "lv: ICE reads the Lv column for the burning material" \
                                                                               || ko "lv: ICE does not report the Lv column"
  if [ $rc = 0 ] && [ "$n" -gt 0 ] && diff -rq "$T/ice-burn-box/OUTPUT" "$d/run/OUTPUT" > /dev/null; then
    ok "lv: the latent heat from the column, read at the burning temperature, writes the same $n OUTPUT files"
  else ko "lv: the run with the Lv column differs (exit=$rc): $(diff -rq "$T/ice-burn-box/OUTPUT" "$d/run/OUTPUT" 2>&1 | head -1)"; fi
else ko "lv: the variant could not be made"; fi
d=$T/lv-both; log=$T/lv-both.log; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"; column "$d"
"$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
if [ $rc != 0 ] && grep -q "is set and the table gives the Lv column" "$log" && ! grep -q "GAME OVER" "$log"; then
  ok "lv: latent-heat beside the Lv column is refused (the table owns the property)"
else ko "lv: latent-heat beside the Lv column was not refused (exit=$rc)"; fi
d=$T/lv-none; log=$T/lv-none.log; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
sed -i '/^latent-heat   = 2.26e6$/d' "$d/input.ini"
"$RUN" "$d" "$d/run" 1 "$log" "$EXE" > /dev/null 2>&1; rc=$?
if ! grep -q '^latent-heat' "$d/input.ini" && [ $rc != 0 ] && grep -q "combustion = direct requires its latent heat" "$log" && ! grep -q "GAME OVER" "$log"; then
  ok "lv: neither the column nor the key is refused"
else ko "lv: a direct material without a latent heat was not refused (exit=$rc)"; fi

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
for c in ice-burn-box-beckstead ice-burn-box; do
  n=$(grep -a -c 'the same from the material at .* K, where it burns' "$T/$c.log"); w=$(grep -a -c '\[WARNING\] q-comb = ' "$T/$c.log")
  [ "$n" = 1 ] && [ "$w" = 0 ] && ok "heat: $c as shipped prints the heat where it burns and no warning" \
                               || ko "heat: $c as shipped: $n heat line(s), $w warning(s)"
done
heat near   ice-burn-box-beckstead   930000  9.3000000E+05 0
heat far    ice-burn-box-beckstead   1030000 1.0300000E+06 1
heat direct ice-burn-box 2000000 2.0000000E+06 0
grep -a -q 'combustion = direct: q-comb sets the rate' "$T/heat-direct.log" && ! grep -a -q 'combustion = direct: q-comb sets the rate' "$T/heat-near.log" \
  && ok "heat direct: the coupler says what q-comb is for direct, and only for direct" || ko "heat direct: the line on q-comb of direct is missing, or printed for Beckstead"
# q-comb = 0 written: allowed; direct says its rate ignores the heat of reaction, Beckstead (beta-part 0.3) warns
heat zero      ice-burn-box           0 0.0000000E+00 0
heat zero-beck ice-burn-box-beckstead 0 0.0000000E+00 1
grep -a -q 'q-comb = 0: the rate is pure vaporisation' "$T/heat-zero.log" && ! grep -a -q 'q-comb = 0: the rate is pure vaporisation' "$T/heat-direct.log" \
  && ok "heat zero: direct with q-comb = 0 says its rate ignores the heat of reaction, and only then" \
  || ko "heat zero: the line on q-comb = 0 is missing, or printed with q-comb > 0"

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
