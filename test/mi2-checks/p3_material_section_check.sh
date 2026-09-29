#!/bin/bash
# P3 check: IGLOO per-material keys come from the material line of the phase file (key=value after
# "<name> <groups>", written by ATLAS GPB from [GPB-Phase*]); no INI section carries them any more.
#  inert   : evap-box-ta as shipped -> OUTPUT byte-identical to the reference run (REF), log never mentions
#            [IGLOO-Material
#  live    : token alpha-e=0.5 on line 2 of INPUT/part-phase.txt -> the GAS FIELD must change.
#            NOT the audited evaporated mass: on this case every stream evaporates completely inside the box,
#            so the audited total equals the injected mass whatever the rate -- alpha-e moves WHERE the vapour
#            is released, not how much. gas-field.tec (and part-source*/euler*) are the rate-sensitive artefacts.
#  refuse  : an unknown key on the material line stops the run naming the key
#  ini     : a stale [IGLOO-Material1] alpha-e = 0.5 in the INI is ignored (gas field == inert) -- the section is gone
#   usage: p3_material_section_check.sh [scratch-root]   (REF=<OUTPUT dir of the reference run> to enable the byte check)
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P3}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
S=$H/test/evap-box-ta
RUN=$H/test/mi2-checks/run_case.sh
REF=${REF:-$S/OUTPUT}
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
mk() { local d=$T/$1; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"; echo "$d"; }
inj() { grep -A3 "closed-domain mass balance" "$1" | grep "^   1  " | awk '{print $2}'; }

d=$(mk inert); "$RUN" "$d" "$d/run" 1 "$d/run.log" > /dev/null 2>&1; rc=$?
grep -q "^\[IGLOO-Material" "$d/input.ini" && ko "shipped INI still carries an [IGLOO-Material*] section" || ok "shipped INI has no [IGLOO-Material*] section"
[ $rc = 0 ] && ok "shipped case runs (exit 0)" || ko "shipped case exit=$rc"
grep -q "IGLOO-Material" "$d/run.log" && ko "log mentions IGLOO-Material" || ok "log never mentions IGLOO-Material"
if [ -d "$REF" ]; then
  diff -rq "$REF" "$d/run/OUTPUT" > /dev/null && ok "OUTPUT byte-identical to the reference run" || ko "OUTPUT differs from the reference run ($REF)"
else
  echo "SKIP  no reference run at $REF: the byte comparison was not made"
fi
m1=$(inj "$d/run.log"); echo "  evaporated (alpha-e default): $m1"

d=$(mk live); sed -i "2s/.*/A 1 alpha-e=${ALPHA:-0.5} interface=LK/" "$d/INPUT/part-phase.txt"
"$RUN" "$d" "$d/run" 1 "$d/run.log" > /dev/null 2>&1; rc=$?
m2=$(inj "$d/run.log"); echo "  evaporated (alpha-e = ${ALPHA:-0.5}): $m2  (expected EQUAL: complete evaporation in the box)"
if [ $rc = 0 ] && ! cmp -s "$T/inert/run/OUTPUT/gas-field.tec" "$d/run/OUTPUT/gas-field.tec"; then
  ok "material-line token alpha-e is LIVE (gas field changed)"; else ko "alpha-e token had no effect on the gas field (exit=$rc)"; fi
printed=$(grep -o ">> \[material 1\].*alpha-e= *[0-9.E+-]*" "$d/run.log" | head -1 | awk '{print $NF}')
python3 -c "import sys; sys.exit(0 if abs(float('$printed') - float('${ALPHA:-0.5}')) < 1e-12 else 1)" 2>/dev/null && grep -q ">> \[material 1\].* intf=1 " "$d/run.log" && ok "resolved-settings line reports the token values" || ko "resolved-settings line missing/wrong"

d=$(mk refuse); sed -i '2s/.*/A 1 no-such-key=1/' "$d/INPUT/part-phase.txt"
"$RUN" "$d" "$d/run" 1 "$d/run.log" > /dev/null 2>&1; rc=$?
[ $rc != 0 ] && grep -q 'unknown key "no-such-key"' "$d/run.log" && ok "unknown material-line key refused, key named (exit=$rc)" || ko "unknown key not refused (exit=$rc)"

d=$(mk ini); printf '\n[IGLOO-Material1]\nalpha-e = 0.5\n' >> "$d/input.ini"
"$RUN" "$d" "$d/run" 1 "$d/run.log" > /dev/null 2>&1; rc=$?
[ $rc = 0 ] && cmp -s "$T/inert/run/OUTPUT/gas-field.tec" "$d/run/OUTPUT/gas-field.tec" && ok "a stale [IGLOO-Material1] section is ignored (gas field unchanged)" || ko "stale [IGLOO-Material1] section changed the run (exit=$rc)"

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
