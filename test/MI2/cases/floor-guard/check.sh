#!/bin/bash
# check.sh <run-dir> <exe>
# The consumed-species floor guard (MI2_Coupling): it compares rates [kg/s], never a rate with a mass.
#   consume: igloo-evap with its record set to H2O 21, N2 -20 (20 kg of N2 consumed per kg of water evaporated; the
#            row still sums to 1). The N2 sink, ~1e-6 kg/s per cell against ~1e-7 kg of N2 per cell, empties no cell
#            within a step: no clamp, and every species row of the mass balance closes to 1e-11. A guard that compares
#            kg/s with kg clamps it and the N2 balance is off by O(1).
#   control: the row H2O -1, N2 2 on a box without vapour (the gas at the floor): the guard fires, its warning once.
source "$(dirname "${BASH_SOURCE[0]}")/../../checklib.sh"

# row <dir> <nu H2O> <nu N2>: the one-term record of A as two terms
row() {
  sed -i -e 's/^material A evaporation 1 /material A evaporation 2 /' \
         -e "s/^  H2O 1.0000000000e+00 1.0000000000e+00\$/  H2O $2 0.0\n  N2 $3 0.0/" "$1/INPUT/part-products.txt"
  grep -q "^  N2 $3 0.0$" "$1/INPUT/part-products.txt" || ko "could not write the record of A in $1"
}
worst() {  # the largest rel.err of the closed-domain mass balance's species rows
  awk '/closed-domain mass balance/{p=1; next} p && NF == 4 && $1 ~ /^[0-9]+$/ {if ($4 + 0 > m) m = $4 + 0}
       p && /ENERGY/ {p=0} END {printf "%.3e", m}' "$1"
}

d=$(variant consume "$cases/igloo-evap"); row "$d" 2.1000000000e+01 -2.0000000000e+01
runin "$d"; rc=$?; n=$(grep -c 'floor guard' "$d/run.log"); e=$(worst "$d/run.log")
echo "      consume: exit $rc, $n floor-guard lines, worst species rel.err $e"
expect "consume: the run ends" ended "$d/run.log"
expect "consume: no clamp of a sink far below the cell-emptying rate" [ "$n" = 0 ]
expect "consume: every species row of the mass balance closes (<= 1e-11)" python3 -c "import sys; sys.exit(float('$e') > 1e-11)"

d=$(variant control "$cases/igloo-evap"); row "$d" -1.0000000000e+00 2.0000000000e+00
runin "$d"; rc=$?; n=$(grep -c 'floor guard' "$d/run.log")
echo "      control: exit $rc, $n floor-guard lines"
expect "control: the guard fires on an empty species" [ "$n" -ge 1 ]
expect "control: its warning is printed once" [ "$n" = 1 ]
finish
