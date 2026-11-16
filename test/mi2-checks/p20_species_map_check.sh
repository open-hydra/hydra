#!/bin/bash
# Species-map check (MI2_Stoichiometry.f90 Load_Species_Map): a material that evaporates gets its gas species from
# INPUT/<phase>products.txt, and a file that does not say them stops the run before the first step.
#  control : evap-box-ta as shipped runs, and its log names the file and the process
#  each row below: evap-box-ta with ONE edit must stop with exit /= 0 and the message
#   usage: p20_species_map_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P20}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/evap-box-ta
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
mk() { local d=$T/$1; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"; echo "$d"; }
F=INPUT/part-products.txt

d=$(mk control)
"$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
[ $rc = 0 ] && grep -q "IGLOO material -> gas species .*$F" "$d/run.log" && grep -q -- "- A \[evaporation\] :" "$d/run.log" \
  && ok "control: the shipped case runs and reports the file and the process" || ko "control: exit=$rc, or the report lines are missing"

# refuse <name> <message> <sed script on the products file | a shell command run in the case directory>
refuse() {
  local d; d=$(mk "$1")
  ( cd "$d" && eval "$3" ) || { ko "$1: the edit failed"; return; }
  if cmp -s "$d/$F" "$S/$F" && cmp -s "$d/input.ini" "$S/input.ini"; then ko "$1: the edit changed nothing"; return; fi
  "$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
  if [ $rc != 0 ] && grep -q -F "$2" "$d/run.log" && ! grep -q "GAME OVER" "$d/run.log"; then ok "$1: refused ($2)"
  else ko "$1: exit=$rc, message \"$2\" $(grep -q -F "$2" "$d/run.log" && echo printed || echo absent)"; fi
}
refuse missing-file   "is missing: run ATLAS GPB"                    "rm $F"
refuse other-material "names the material B, which IGLOO does not"   "sed -i 's/^material A /material B /' $F"
refuse other-species  "names the gas species CO2, which the gas"     "sed -i 's/^  H2O /  CO2 /' $F"
refuse other-gas      "was written for another gas phase"            "sed -i 's/^gas gas- 2\$/gas gas- 3/' $F"
refuse no-record      "has no evaporation record for it"             "sed -i -e 's/^material A evaporation 1 latent .*/material A none 0 none 0.0000000000e+00/' -e '/^  H2O /d' $F"
refuse row-sum        "does not sum to 1"                            "sed -i 's/^  H2O 1.0000000000e+00 /  H2O 9.0000000000e-01 /' $F"
refuse energy-word    "an evaporation record that is not"            "sed -i 's/ evaporation 1 latent / evaporation 1 formation /' $F"
refuse truncated      "fewer terms than its last record announces"   "sed -i '/^  H2O /d' $F"
refuse retired-key    "products-A is no longer read"                 "sed -i 's/^\[MI2-Coupling\]\$/[MI2-Coupling]\nproducts-A = H2O:1.0/' input.ini"
refuse euler-cycle    "a cycle has IGLOO-iter = 2"                   "sed -i 's/^IGLOO-iter    = 1\$/IGLOO-iter    = 2/' input.ini"
refuse old-word       "an evaporation record that is not"            "sed -i 's/ evaporation 1 latent / evaporation 1 explicit /' $F"

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
