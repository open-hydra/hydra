#!/bin/bash
# check.sh <run-dir> <exe>
# The burning twins: one burning section [GPB-Phase2] in an ICE case and in an IGLOO case. Each case meets its own
# closed forms (its base's check.py), and the two solvers agree on the energy the gas receives per kg burnt and on
# the constant of the burning law: ICE's from the rate of one drop in its frozen source, IGLOO's from the diameters
# of its parcels at the end of a run cut inside the burn (the law variant: time-end 2 ms, the drops entering at
# 1e-6 m/s, so that the Nusselt number of direct is 2 to the tolerance).
#   direct   : ice-burn, and igloo/ (igloo-burn with ice-burn's burning section and gas);
#   Beckstead: the two with the section's burning keys set to Beckstead down to a 1 um residue;
#   burnout  : ice-burnout, and igloo/ on Beckstead removed at d-burnout 3 um: the energy per kg handed.
source "$(dirname "${BASH_SOURCE[0]}")/../../checklib.sh"
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TOL_E=1e-9        # energy per kg, ICE against IGLOO: the residue's end state is read with seven figures
TOL_BECK=2e-6     # K-burn X-eff, ICE against IGLOO: seven printed figures of the parcel's diameter
TOL_DIRECT=3e-4   # K of direct: Ranz-Marshall at 1e-6 m/s is 2.3e-4 above Nu 2 at injection

source "$here/beckstead.sh"
num() { echo "$1" | grep -m1 -F -- "$2" | awk -v k="$2" '{print substr($0, index($0, k) + length(k))}' | awk '{print $1}'; }
close() { python3 -c "import sys; a, b, t = map(float, sys.argv[1:4]); sys.exit(not abs(a - b) <= t * abs(b))" "$1" "$2" "$3"; }
rel() { python3 -c "import sys; a, b = map(float, sys.argv[1:3]); print('%.2e' % (abs(a - b) / abs(b)))" "$1" "$2"; }
# judged <label> <dir> <check.py> [args]: run, then the base's check; its output on stdout of the caller's $out
judged() {
  local label=$1 d=$2 py=$3; shift 3; runin "$d"; local rc=$?
  out=$(python3 -B "$py" "$d" "$d/run.log" "$@" 2>&1)
  if [ $rc = 0 ] && echo "$out" | grep -q '^== .*[0-9]* passed, 0 failed'; then ok "$label: runs and meets its closed forms"
  else ko "$label: exit $rc or its check fails (see $d/run.log)"; echo "$out" | grep '^FAIL' | head -3; fi
}
law() {  # law <dir>: the IGLOO run cut inside the burn
  sed -i -e 's/^time-end    = 0.05$/time-end    = 2.0e-3/' -e 's/^up    = 1.0$/up    = 1.0e-6/' "$1/input.ini"
  grep -q '^time-end    = 2.0e-3$' "$1/input.ini" && grep -q '^up    = 1.0e-6$' "$1/input.ini" || ko "the law variant of $1 could not be made"
}
pair() {  # pair <word> <tol of the law constant>: the ICE twin in $run/<word>-ice, the IGLOO one in $run/<word>-igloo
  local w=$1 oi ol ow
  a=$(grep '^A ' "$run/$w-ice/INPUT/part-models.txt"); b=$(grep '^A ' "$run/$w-igloo/INPUT/part-models.txt")
  expect "$w: one row in the two models files ($a)" [ -n "$a" -a "$a" = "$b" ]
  judged "$w: ICE" "$run/$w-ice" "$cases/ice-burn/check.py"; oi=$out
  judged "$w: IGLOO" "$run/$w-igloo" "$cases/igloo-burn/check.py"; ol=$out
  cp -rp "$run/$w-igloo" "$run/$w-law"; law "$run/$w-law"
  judged "$w: IGLOO cut inside the burn" "$run/$w-law" "$cases/igloo-burn/check.py" law; ow=$out
  ei=$(num "$oi" 'M + beta q'); el=$(num "$ol" '(Q - residue)/M')
  ki=$(num "$oi" 'twin: law constant'); kl=$(num "$ow" 'twin: law constant')
  expect "$w: energy per kg burnt, ICE ${ei:-?} against IGLOO ${el:-?} J/kg (tol $TOL_E)" close "${ei:-nan}" "${el:-nan}" $TOL_E
  expect "$w: constant of the law, ICE ${ki:-?} against IGLOO ${kl:-?} (tol $2)" close "${ki:-nan}" "${kl:-nan}" "$2"
}

variant direct-ice "$cases/ice-burn" > /dev/null; variant direct-igloo "$here/igloo" > /dev/null
pair direct $TOL_DIRECT
d=$(variant Beckstead-ice "$cases/ice-burn"); beckstead "$d" 'd-stop     = 1.0e-6' 'd-stop=1e-06'
d=$(variant Beckstead-igloo "$here/igloo"); beckstead "$d" 'd-stop     = 1.0e-6' 'd-stop=1e-06'
pair Beckstead $TOL_BECK

d=$(variant burnout-ice "$cases/ice-burnout")
judged "burnout: ICE" "$d" "$cases/ice-burnout/check.py"; oi=$out
d=$(variant burnout-igloo "$here/igloo"); beckstead "$d" 'd-burnout  = 3.0e-6' 'd-burnout=3e-06'
a=$(grep '^A ' "$run/burnout-ice/INPUT/part-models.txt"); b=$(grep '^A ' "$d/INPUT/part-models.txt")
expect "burnout: one row in the two models files ($a)" [ -n "$a" -a "$a" = "$b" ]
judged "burnout: IGLOO" "$d" "$cases/igloo-burn/check.py"; ol=$out
ei=$(num "$oi" '(SUM Q dt)/(SUM M dt)'); el=$(num "$ol" '(Q - residue)/M')
expect "burnout: energy per kg handed, ICE ${ei:-?} against IGLOO ${el:-?} J/kg (tol $TOL_E)" close "${ei:-nan}" "${el:-nan}" $TOL_E
finish
