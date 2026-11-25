#!/bin/bash
# check.sh <run-dir> <exe>
# The phase routing of hydra-MI2 (MI2_Mod_Import): which solver takes which INPUT/<phase>phase.txt. Each deck is
# igloo-evap cut to one MOSE iteration, without its solver phase key and without ATLAS's modeling token on line 1 of
# part-phase.txt, then one edit; each must exit as expected and print the expected line.
source "$(dirname "${BASH_SOURCE[0]}")/../../checklib.sh"

mk() {
  local d; d=$(variant "$1" "$cases/igloo-evap")
  sed -i -e 's/^N-coupling *=.*/N-coupling    = 1/' -e 's/^MOSE-iter *=.*/MOSE-iter     = 1/' \
         -e '/^phase *= *part/d' "$d/input.ini"
  sed -i '1s/ modeling=.*//' "$d/INPUT/part-phase.txt"
  echo "$d"
}
key() { sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = part/' "$1/input.ini"; }
both() { sed -i 's/^IGLOO-iter *=.*/IGLOO-iter    = 1\nICE-iter      = 1/' "$1/input.ini"; }
# leg <dir> <0|nz> <fragment> <label>
leg() {
  runin "$1"; local rc=$?
  if { [ "$2" = 0 ] && [ $rc -eq 0 ]; } || { [ "$2" = nz ] && [ $rc -ne 0 ]; }; then
    grep -qF -- "$3" "$1/run.log" && { ok "$4 (exit $rc)"; return; }
  fi
  ko "$4: exit $rc, wanted $2 and \"$3\" (see $1/run.log)"
}

d=$(mk legacy);          leg "$d" 0  "DEPRECATED: [IGLOO-General] phase not set" "no key, one condensed phase: taken, with a warning"
d=$(mk key-ok);          key "$d"; leg "$d" 0 'IGLOO  prefix "part-"  (from its phase key)' "the key names the phase"
d=$(mk key-wrong);       sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = nosuch/' "$d/input.ini"
                         leg "$d" nz "names no condensed phase file" "a key that names no file"
d=$(mk token-ok);        sed -i '1s/.*/condensed-dispersed modeling=lagrangian/' "$d/INPUT/part-phase.txt"; key "$d"
                         leg "$d" 0  "modeling=lagrangian)  --> IGLOO" "ATLAS's token agrees with the solver"
d=$(mk token-mismatch);  sed -i '1s/.*/condensed-dispersed modeling=eulerian/' "$d/INPUT/part-phase.txt"; key "$d"
                         leg "$d" nz 'is tagged "modeling=eulerian" by ATLAS' "ATLAS's token contradicts the key"
d=$(mk token-nokey);     sed -i '1s/.*/condensed-dispersed modeling=eulerian/' "$d/INPUT/part-phase.txt"
                         leg "$d" nz "0 condensed phase file(s) could be meant" "the token excludes the only file"
d=$(mk both-one-phase);  both "$d"
                         leg "$d" nz "ICE and IGLOO both resolve to INPUT/part-phase.txt" "two solvers, one phase, no keys"
d=$(mk both-same-key);   both "$d"; key "$d"; printf '\n[ICE-Parameters]\nphase = part\n' >> "$d/input.ini"
                         leg "$d" nz "ICE and IGLOO both resolve to INPUT/part-phase.txt" "two keys name one phase"
d=$(mk dash-name);       for f in "$d"/INPUT/part-*; do mv "$f" "${f/part-/part-E-}"; done
                         leg "$d" nz 'contains "-"' "a dash in a phase name"
d=$(mk substring);       for f in "$d"/INPUT/part-*; do cp -p "$f" "${f/part-/partL-}"; done
                         sed -i 's/^\[IGLOO-General\]/[IGLOO-General]\nphase = partL/' "$d/input.ini"
                         leg "$d" nz "is a substring of" "one phase name inside another"
d=$(mk unnamed);         for f in "$d"/INPUT/part-*; do mv "$f" "${f/part-/}"; done
                         leg "$d" nz "is an UNNAMED phase but the case has" "an unnamed phase beside named ones"
d=$(mk missing-bc);      key "$d"; rm -f "$d/INPUT/part-bc.txt"
                         leg "$d" nz "INPUT/part-bc.txt does not exist" "a file of the phase missing"
d=$(mk no-condensed);    rm -f "$d"/INPUT/part-*
                         leg "$d" nz "0 condensed phase file(s) could be meant" "IGLOO on, no condensed phase"
d=$(mk unused);          sed -i 's/^IGLOO-iter *=.*//' "$d/input.ini"
                         leg "$d" 0  "UNUSED (no enabled solver names it)" "a phase no solver takes: MOSE alone"
d=$(mk long-name);       for f in "$d"/INPUT/gas-*; do mv "$f" "${f/gas-/averyveryverylonggasname-}"; done
                         leg "$d" nz "does not fit MOSE's prefix buffer" "a name longer than MOSE's prefix"
d=$(mk liquid);          sed -i '1s/.*/liquid-dispersed modeling=lagrangian/' "$d/INPUT/part-phase.txt"; key "$d"
                         leg "$d" 0  "(liquid-dispersed, modeling=lagrangian)  --> IGLOO" "liquid-dispersed to IGLOO"
d=$(mk solid);           sed -i '1s/.*/solid-dispersed/' "$d/INPUT/part-phase.txt"; key "$d"
                         leg "$d" 0  "(solid-dispersed, no modeling token)  --> IGLOO" "solid-dispersed to IGLOO"
d=$(mk solid-bulk);      sed -i '1s/.*/solid-bulk/' "$d/INPUT/part-phase.txt"; key "$d"
                         leg "$d" nz "is a SOLID BULK phase" "solid-bulk refused"
d=$(mk material-badkey); sed -i '2s/.*/A 1 no-such-key=1/' "$d/INPUT/part-phase.txt"; key "$d"
                         leg "$d" nz 'unknown key "no-such-key"' "an unknown key on the material line"
d=$(mk schedule);        sed -i -e 's/^MOSE-iter *=.*/MOSE-iter     = 1 1/' -e 's/^IGLOO-iter *=.*/IGLOO-iter    = 1 1/' "$d/input.ini"
                         leg "$d" nz "N-coupling has 1 values but MOSE-iter has 2" "a schedule list shorter than MOSE-iter"
finish
