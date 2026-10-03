#!/bin/bash
# Routing checks for hydra-MI2 Assign_Phases. Each case is a scratch copy of test/evap-box-ta
# with one permutation applied; asserts exit code and a log fragment.
#   usage: p1_routing_checks.sh <scratch-root> [exe]
set -u
ROOT=${1:?usage: p1_routing_checks.sh <scratch-root> [exe]}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
SRC=$H/test/evap-box-ta
pass=0; fail=0

mk() {  # mk <name>  -> fresh scratch case dir, fast schedule
  local d=$ROOT/$1
  rm -rf "$d"; mkdir -p "$d"
  cp -rp "$SRC/INPUT" "$SRC/MESH" "$SRC/input.ini" "$d/"
  sed -i -e 's/^N-coupling *=.*/N-coupling    = 1/' -e 's/^MOSE-iter *=.*/MOSE-iter     = 1/' -e '/^phase *= *part/d' "$d/input.ini"   # start every scenario WITHOUT a solver phase key
  sed -i '1s/ modeling=.*//' "$d/INPUT/part-phase.txt"                                                                                    # ... and WITHOUT the ATLAS line-1 token
  echo "$d"
}
run() {  # run <dir> <expected-exit(0|nz)> <grep-fragment> <label>
  local d=$1 want=$2 frag=$3 label=$4
  ( cd "$d" && rm -rf OUTPUT && mkdir -p OUTPUT && ls INPUT/*phase.txt > filelist.txt \
    && ulimit -s unlimited && OMP_NUM_THREADS=1 KMP_STACKSIZE=100M "$EXE" > run.log 2>&1 )
  local rc=$?
  local ok=1
  if [ "$want" = 0 ] && [ $rc -ne 0 ]; then ok=0; fi
  if [ "$want" = nz ] && [ $rc -eq 0 ]; then ok=0; fi
  grep -q -- "$frag" "$d/run.log" || ok=0
  if [ $ok = 1 ]; then pass=$((pass+1)); echo "PASS  $label  (exit=$rc)"
  else fail=$((fail+1)); echo "FAIL  $label  (exit=$rc, wanted $want, fragment '$frag')"; grep -n "ERROR\|DEPRECATED\|routing" "$d/run.log" | head -5; fi
}

# 1 legacy: no key, one untagged condensed phase, one solver -> runs, deprecation
d=$(mk legacy);                 run "$d" 0  "DEPRECATED: \[IGLOO-General\] phase not set" "legacy fallback"
# 2 key ok
d=$(mk key-ok);                 sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = part/' "$d/input.ini"
                                run "$d" 0  "IGLOO  prefix \"part-\"  (from its phase key)" "key names the phase"
# 3 key names nothing
d=$(mk key-wrong);              sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = nosuch/' "$d/input.ini"
                                run "$d" nz "names no condensed phase file" "key names no file"
# 4 ATLAS token agrees with the solver
d=$(mk token-ok);               sed -i '1s/.*/condensed-dispersed modeling=lagrangian/' "$d/INPUT/part-phase.txt"
                                sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = part/' "$d/input.ini"
                                run "$d" 0  "modeling=lagrangian)  --> IGLOO" "token agrees"
# 5 ATLAS token contradicts the solver
d=$(mk token-mismatch);         sed -i '1s/.*/condensed-dispersed modeling=eulerian/' "$d/INPUT/part-phase.txt"
                                sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = part/' "$d/input.ini"
                                run "$d" nz "is tagged \"modeling=eulerian\" by ATLAS" "token contradicts"
# 6 token contradicts, no key -> no candidate for IGLOO
d=$(mk token-nokey);            sed -i '1s/.*/condensed-dispersed modeling=eulerian/' "$d/INPUT/part-phase.txt"
                                run "$d" nz "0 condensed phase file(s) could be meant" "token excludes the only file"
# 7 both solvers, one untagged phase, no keys -> double-injection refused
d=$(mk both-one-phase);         sed -i 's/^IGLOO-iter *=.*/IGLOO-iter    = 1\nICE-iter      = 1/' "$d/input.ini"
                                run "$d" nz "ICE and IGLOO both resolve to INPUT/part-phase.txt" "both solvers one phase"
# 8 both solvers, both keys name the same phase -> refused
d=$(mk both-same-key);          sed -i 's/^IGLOO-iter *=.*/IGLOO-iter    = 1\nICE-iter      = 1/' "$d/input.ini"
                                sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = part/' "$d/input.ini"; printf '\n[ICE-Parameters]\nphase = part\n' >> "$d/input.ini"
                                run "$d" nz "ICE and IGLOO both resolve to INPUT/part-phase.txt" "both keys same phase"
# 9 dash in the phase name
d=$(mk dash-name);              for f in "$d"/INPUT/part-*; do mv "$f" "${f/part-/part-E-}"; done
                                run "$d" nz "contains \"-\"" "dash in name"
# 10 substring names (part / partL)
d=$(mk substring);              for f in "$d"/INPUT/part-*; do cp -p "$f" "${f/part-/partL-}"; done
                                sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = partL/' "$d/input.ini"
                                run "$d" nz "is a substring of" "substring names"
# 11 unnamed phase next to named ones
d=$(mk unnamed);                for f in "$d"/INPUT/part-*; do mv "$f" "${f/part-/}"; done
                                run "$d" nz "is an UNNAMED phase but the case has" "unnamed + named"
# 12 key ok but a required file is missing
d=$(mk missing-bc);             sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = part/' "$d/input.ini"; rm -f "$d/INPUT/part-bc.txt"
                                run "$d" nz "INPUT/part-bc.txt does not exist" "missing bc file"
# 13 IGLOO enabled but no condensed phase at all
d=$(mk no-condensed);           rm -f "$d"/INPUT/part-*
                                run "$d" nz "0 condensed phase file(s) could be meant" "no condensed phase"
# 14 condensed phase present, its solver disabled -> UNUSED, runs (MOSE only)
d=$(mk unused);                 sed -i 's/^IGLOO-iter *=.*//' "$d/input.ini"
                                run "$d" 0  "UNUSED (no enabled solver names it)" "phase unused"

# 15 phase name longer than the solver's prefix buffer (MOSE/ICE clen=16 -> 15 chars max)
d=$(mk long-name);              for f in "$d"/INPUT/gas-*; do mv "$f" "${f/gas-/averyveryverylonggasname-}"; done
                                run "$d" nz "does not fit MOSE's prefix buffer" "name overflows prefix buffer"

# 16 ATLAS type words (P3.4): liquid-dispersed and solid-dispersed route to the condensed solver
d=$(mk liquid-dispersed);       sed -i '1s/.*/liquid-dispersed modeling=lagrangian/' "$d/INPUT/part-phase.txt"
                                sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = part/' "$d/input.ini"
                                run "$d" 0  "(liquid-dispersed, modeling=lagrangian)  --> IGLOO" "liquid-dispersed routes to IGLOO"
d=$(mk solid-dispersed);        sed -i '1s/.*/solid-dispersed/' "$d/INPUT/part-phase.txt"
                                sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = part/' "$d/input.ini"
                                run "$d" 0  "(solid-dispersed, no modeling token)  --> IGLOO" "solid-dispersed routes to IGLOO"
# 17 solid-bulk is refused (hydra-MI2 has no solid solver)
d=$(mk solid-bulk);             sed -i '1s/.*/solid-bulk/' "$d/INPUT/part-phase.txt"
                                sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = part/' "$d/input.ini"
                                run "$d" nz "is a SOLID BULK phase" "solid-bulk refused"
# 18 material-line tokens (P3.3): IGLOO applies key=value after "<name> <groups>"
d=$(mk material-tokens);        sed -i '2s/.*/A 1 alpha-e=0.5 interface=LK/' "$d/INPUT/part-phase.txt"
                                sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = part/' "$d/input.ini"
                                run "$d" 0  "\[material 1\] evap=2 liq=0 intf=1 boil=0 comb=0 solid=0 alpha-e= 5.000E-01" "material-line tokens applied"
# 19 an unknown material-line key stops the run
d=$(mk material-badkey);        sed -i '2s/.*/A 1 no-such-key=1/' "$d/INPUT/part-phase.txt"
                                sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = part/' "$d/input.ini"
                                run "$d" nz "unknown key \"no-such-key\"" "unknown material key refused"

echo "== $pass passed, $fail failed"
[ $fail = 0 ]
