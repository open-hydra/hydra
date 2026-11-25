#!/bin/bash
# check.sh <run-dir> <exe>
# The [MI2-Channels] lines hydra-MI2 refuses at setup. Each deck is a time-accurate case of test/MI2/cases with one
# channel line (igloo-evap's for relaxation = true: igloo-al's generated mass refuses it first), and must stop before
# the first step with the expected message.
source "$(dirname "${BASH_SOURCE[0]}")/../../checklib.sh"

# refuse <name> <base case> <channel line> <expected message> [sed edit of input.ini]
refuse() {
  local d; d=$(variant "$1" "$cases/$2")
  printf '\n[MI2-Channels]\nbad = %s\n' "$3" >> "$d/input.ini"
  [ -n "${5:-}" ] && sed -i "$5" "$d/input.ini"
  RUN_TIMEOUT=300 runin "$d"; local rc=$?
  if [ $rc -ne 0 ] && grep -qF -- "$4" "$d/run.log" && ! grep -q 'Iter =' "$d/run.log"; then ok "$1: refused at setup ($4)"
  else ko "$1: exit $rc; expected \"$4\" before any step (see $d/run.log)"; fi
}

refuse syntax        igloo-al   'gas:H2O'                                   'does not read "<from> <to> <law>'
refuse unknown-law   igloo-al   'gas:H2O IGLOO:AL2O3(L) nolaw rate=1'       'no law is called nolaw'
refuse no-constant   igloo-al   'gas:H2O IGLOO:AL2O3(L) constant'           'the law constant needs rate='
refuse bad-d         igloo-al   'gas:H2O IGLOO:AL2O3(L) constant rate=1 birth-d=-1' 'is not a diameter > 0'
refuse no-species    igloo-al   'gas:XX IGLOO:AL2O3(L) constant rate=1'     'XX is not a species of the gas'
refuse no-solver     igloo-al   'gas:H2O MOON:AL2O3(L) constant rate=1'     'names no solver'
refuse to-gas        igloo-al   'IGLOO:AL2O3(L) gas:H2O constant rate=1'    'a channel ends in ICE or IGLOO'
refuse same-solver   igloo-al   'IGLOO:AL(L) IGLOO:AL2O3(L) constant rate=1' 'a channel joins two solvers'
refuse igloo-donor   igloo-al   'IGLOO:AL(L) ICE:AL2O3(L) relaxation k=1'   "needs the donor's mass in a cell"
refuse ice-unsteady  igloo-al   'gas:H2O ICE:AL2O3(L) constant rate=1 birth-d=3e-6' 'a channel with an ICE end needs a steady'
refuse relaxation    igloo-evap 'gas:H2O IGLOO:A constant rate=1' \
       'relaxation = true blends two sweeps' 's/^relaxation *= *false/relaxation = true/'
refuse solidifies    ice-al     'gas:H2O ICE:AL2O3(L) constant rate=1'      'solidifies: an absorbed mass'

finish
