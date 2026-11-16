#!/bin/bash
# burning twins gate: the same burning section [GPB-Phase2] in an ICE case and in an IGLOO case, once per model.
#   direct     test/ice-burn-box            and  test/igloo-burn-box-direct
#   Beckstead  test/ice-burn-box-beckstead  and  test/igloo-burn-box-beckstead
# For each pair: the two material lines of INPUT/part-phase.txt and the two rows of INPUT/part-models.txt are the
# same; the energy the gas receives per kg
# burnt is the same in the two solvers; the constant of the burning law is the same, ICE's from the rate of one drop
# in its frozen source, IGLOO's from the diameter of its parcels at the end of a run cut inside the burn (the law
# variant: time-end 2 ms, the drops entering at 1e-6 m/s so that the Nusselt number of direct is 2 within the
# tolerance check.py prints). Each case against its own closed forms: p21 (ICE), p22 (IGLOO).
#   usage: p23_burn_twins_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P23}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
RUN=$H/test/mi2-checks/run_case.sh
TOL_E=1e-9        # energy per kg burnt, ICE against IGLOO: the residue's end state is read with seven figures
TOL_BECK=2e-6     # K-burn X-eff, ICE against IGLOO: seven printed figures of the parcel's diameter
TOL_DIRECT=3e-4   # K of direct, ICE against IGLOO: Ranz-Marshall at 1e-6 m/s is 2.3e-4 above Nu 2 at injection
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
num() { echo "$1" | grep -m1 "$2" | sed -e "s/.*$2 *//" | awk '{print $1}'; }
close() { python3 -c "import sys; a, b, t = map(float, sys.argv[1:4]); sys.exit(0 if abs(a - b) <= t * abs(b) else 1)" "$1" "$2" "$3"; }
rel() { python3 -c "import sys; a, b = map(float, sys.argv[1:3]); print('%.2e' % (abs(a - b) / abs(b)))" "$1" "$2"; }
echo "  exe: $EXE ($(md5sum < "$EXE" | cut -c1-8))"

for word in direct Beckstead; do
  if [ $word = direct ]; then I=ice-burn-box; L=igloo-burn-box-direct; tol=$TOL_DIRECT
  else I=ice-burn-box-beckstead; L=igloo-burn-box-beckstead; tol=$TOL_BECK; fi
  echo "--- combustion = $word: $I and $L ---"
  a=$(grep '^A ' "$H/test/$I/INPUT/part-phase.txt"); b=$(grep '^A ' "$H/test/$L/INPUT/part-phase.txt")
  [ -n "$a" ] && [ "$a" = "$b" ] && echo "$a" | grep -q "combustion=$word\$" \
    && ok "$word: one material line in the two phase files ($a)" || ko "$word: the material lines differ: '$a' / '$b'"
  a=$(grep '^A ' "$H/test/$I/INPUT/part-models.txt"); b=$(grep '^A ' "$H/test/$L/INPUT/part-models.txt")
  [ -n "$a" ] && [ "$a" = "$b" ] && echo "$a" | grep -q "^A combustion $word " \
    && ok "$word: one row in the two models files ($a)" || ko "$word: the rows of the models files differ: '$a' / '$b'"

  d=$T/$I; rm -rf "$d"; "$RUN" "$H/test/$I" "$d" 1 "$T/$I.log" "$EXE" > /dev/null 2>&1; rc=$?
  oi=$(python3 -B "$H/test/ice-burn-box/check.py" "$d" "$T/$I.log" "$H/test/ice-box" 2>&1)
  [ $rc = 0 ] && echo "$oi" | grep -q '^== check.py: [0-9]* passed, 0 failed' && ok "$word: ICE case runs and meets its closed forms" \
                                                                           || ko "$word: ICE case exit=$rc or its check fails"
  d=$T/$L; rm -rf "$d"; "$RUN" "$H/test/$L" "$d" 1 "$T/$L.log" "$EXE" > /dev/null 2>&1; rc=$?
  ol=$(python3 -B "$H/test/igloo-burn-box/check.py" "$d" "$T/$L.log" "$H/test/evap-box-mom" 2>&1)
  [ $rc = 0 ] && echo "$ol" | grep -q '^== check.py: [0-9]* passed, 0 failed' && ok "$word: IGLOO case runs and meets its closed forms" \
                                                                           || ko "$word: IGLOO case exit=$rc or its check fails"
  v=$T/$L-law-case; d=$T/$L-law; rm -rf "$v" "$d"; mkdir -p "$v"; cp -rp "$H/test/$L/INPUT" "$H/test/$L/MESH" "$H/test/$L/input.ini" "$v/"
  sed -i -e 's/^time-end    = 0.05$/time-end    = 2.0e-3/' -e 's/^up    = 1.0$/up    = 1.0e-6/' "$v/input.ini"
  if grep -q '^time-end    = 2.0e-3$' "$v/input.ini" && grep -q '^up    = 1.0e-6$' "$v/input.ini"; then
    "$RUN" "$v" "$d" 1 "$T/$L-law.log" "$EXE" > /dev/null 2>&1; rc=$?
    ow=$(python3 -B "$H/test/igloo-burn-box/check.py" "$d" "$T/$L-law.log" "$H/test/evap-box-mom" law 2>&1)
    [ $rc = 0 ] && echo "$ow" | grep -q '^== check.py: [0-9]* passed, 0 failed' && ok "$word: IGLOO parcels follow the law at time-end" \
                                                                             || ko "$word: IGLOO law run exit=$rc or its check fails"
  else ko "$word: the law variant could not be made"; ow=; fi

  ei=$(num "$oi" 'M + beta q'); el=$(num "$ol" '(Q - residue)\/M')
  ki=$(num "$oi" 'twin: law constant'); kl=$(num "$ow" 'twin: law constant')
  if [ -n "$ei" ] && [ -n "$el" ] && close "$ei" "$el" $TOL_E; then ok "$word: energy per kg burnt, ICE $ei against IGLOO $el J/kg: $(rel "$ei" "$el") (tol $TOL_E)"
  else ko "$word: energy per kg burnt, ICE '$ei' against IGLOO '$el' J/kg (tol $TOL_E)"; fi
  if [ -n "$ki" ] && [ -n "$kl" ] && close "$ki" "$kl" $tol; then ok "$word: constant of the law, ICE $ki against IGLOO $kl: $(rel "$ki" "$kl") (tol $tol)"
  else ko "$word: constant of the law, ICE '$ki' against IGLOO '$kl' (tol $tol)"; fi
done

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
