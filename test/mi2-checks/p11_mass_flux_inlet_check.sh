#!/bin/bash
# P11: the gas inflow mass flow the coupler hands IGLOO at a gas inlet (MOSExIGLOO, [MI2-Coupling] report = true).
# MOSE's bc%mdot is the face mass flow it computes [kg/s] for 401/402/405-407 but the PRESCRIBED MASS FLUX g
# [kg/(s m2)] read from gas-bc.txt for 403/404/420; IGLOO injects krho/(1-krho) x |mdotGas| per face cell and
# needs kg/s, so the coupler multiplies g by IGLOO's face-cell area for 403/404/420 and passes the others as
# they are. JPL-Lagrangian-20micron (MOSE + IGLOO, gas inlet 401 on face 1, dispersed 401, krho = 0.30) in a
# short schedule (N-coupling 2, MOSE-iter 100, 1 thread), three runs:
#  inlet403 : the 40 gas inlet records 401 -> 403 with g = 219.39 kg/(s m2) (the converged rho*u of the 401 case),
#             part-bc.txt untouched -> every ' [MI2-report] igloo-inlet-mdot <sum> <n>' line has n = 40 and
#             sum = g * SUM(cell%area) to 1e-8 (the areas recomputed from MESH/mesh.tec as IGLOO does), and the
#             inlet-column condensed flux of IGLOO's Eulerian projection is krho/(1-krho) x the gas rho*u of the
#             same run's field within 5 % (and krho/(1-krho) x g within 5 %);
#  one403   : the same with N-coupling 1 -> the flux check on a binary whose over-injection kills the gas in the
#             second MOSE block (no final write, no euler file) still yields a number;
#  ctrl401  : the case as shipped -> n = 40, MOSE's positive computed inflow within [0.5, 2] x g * SUM(area), and
#             the same 5 % flux check.
# Every run: GAME OVER, none of the failure signatures, exit 0, parcels injected in every sweep.
# Without the area factor the report holds g x n (1/mean(area) = 1.1e6 too large) or is absent, and the second
# MOSE block dies with 'NaN or p<0 detected' ('NaN, infinity or p<0 detected' from a later MOSE). The numeric checks are in p11_mass_flux_inlet.py.
#   usage: p11_mass_flux_inlet_check.sh [scratch-root] [exe]
set -u
export LC_ALL=C
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P11}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
S=$H/test/JPL-Lagrangian-20micron
EXE=${2:-$H/bin/hydra_MI2}
RUN=$H/test/mi2-checks/run_case.sh
PY=$H/test/mi2-checks/p11_mass_flux_inlet.py
G=219.39               # prescribed mass flux of the 403 records [kg/(s m2)]
GTXT=0.219390E+03      # the same number as written in gas-bc.txt
MOSE_ITER=100
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
echo "  exe $EXE ($(md5sum < "$EXE" | cut -c1-8))"
python3 -c "import sys; sys.exit(0 if abs(float('$GTXT') - $G) <= 1e-12 * $G else 1)" || { ko "GTXT $GTXT is not G = $G"; echo "== $pass passed, $fail failed"; exit 1; }

mk() {  # mk <name> <N-coupling> -> scratch copy of the case with the short schedule and the coupling report
  local d=$T/$1; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
  sed -i -e "s/^N-coupling *=.*/N-coupling    = $2/" -e "s/^MOSE-iter *=.*/MOSE-iter     = $MOSE_ITER/" \
         -e "s/^sol-diter *=.*/sol-diter     = $MOSE_ITER/" -e "s/^shell-diter *=.*/shell-diter   = 1/" "$d/input.ini"
  if grep -q '^\[MI2-Coupling\]' "$d/input.ini"; then sed -i '/^\[MI2-Coupling\]/a report = true' "$d/input.ini"
  else printf '\n[MI2-Coupling]\nreport = true\n' >> "$d/input.ini"; fi
  echo "$d"; }
to403() {  # to403 <gas-bc.txt in> <out>: face-1 records 401 -> 403, payload field 2 (p0) -> g, in place
  awk -v gt="$GTXT," 'BEGIN{pl=0}
    { if (pl==1) { n=index($0,$2); $0=substr($0,1,n-1) gt substr($0,n+length($2)); pl=0; print; next }
      if (NF==6 && $1 ~ /^[0-9]+$/ && $6 ~ /^[0-9]+$/ && $5==1 && $6==401) { $6=403; pl=1
        printf "%8d%8d%8d%8d%8d%8d\n",$1,$2,$3,$4,$5,$6; next }
      print }' "$1" > "$2"; }
hygiene() {  # hygiene <name> <log> <exit> <sweeps>
  local bad; bad=$(grep -c -e '\[ERROR\]' -e ' FAIL' -e 'on-physical state' -e 'Integration failed at' -e 'p<0 detected' -e 'ABORT:' "$2")
  if [ "$3" = 0 ] && grep -q 'GAME OVER' "$2" && [ "$bad" = 0 ]; then ok "$1: clean run (GAME OVER, no failure signature, exit 0)"
  else ko "$1: run not clean (exit=$3, GAME OVER $(grep -c 'GAME OVER' "$2"), failure signatures $bad: $(grep -o -e 'Integration failed at' -e 'p<0 detected' -e 'on-physical state' -e '\[ERROR\]' -e 'ABORT:' "$2" | sort -u | tr '\n' ' '))"; fi
  local np; np=$(grep 'number of particles' "$2" | awk '{print $NF}' | tr '\n' ' ')
  if [ "$(grep -c 'number of particles' "$2")" = "$4" ] && ! echo " $np" | grep -q ' 0 '; then ok "$1: parcels injected in every sweep ($np)"
  else ko "$1: parcels per sweep: $np (expected $4 sweeps, all > 0)"; fi
  echo "  $1: $(grep '^# wall' "$2" | head -1)"; }
check() {  # check <name> <dir> <mode>: the numeric checks; a crash of the script counts as a failure
  local out rc; out=$(python3 -B "$PY" "$2/run" "$2/run.log" "$3" "$G" 2>&1); rc=$?
  echo "$out" | sed "s/^\(PASS\|FAIL\)  /\1  $1: /"
  pass=$((pass + $(echo "$out" | grep -c '^PASS'))); fail=$((fail + $(echo "$out" | grep -c '^FAIL')))
  if [ $rc != 0 ] && [ "$(echo "$out" | grep -c '^FAIL')" = 0 ]; then ko "$1: numeric checks did not run (exit $rc)"; fi; }

# --- inlet403: gas 403 with g, two coupling steps ---
d=$(mk inlet403 2); to403 "$S/INPUT/gas-bc.txt" "$d/INPUT/gas-bc.txt"
n401=$(awk 'NF==6 && $1 ~ /^[0-9]+$/ && $6==401' "$S/INPUT/gas-bc.txt" | wc -l); n403=$(awk 'NF==6 && $1 ~ /^[0-9]+$/ && $6==403' "$d/INPUT/gas-bc.txt" | wc -l)
nd=$(diff "$S/INPUT/gas-bc.txt" "$d/INPUT/gas-bc.txt" | grep -c '^[<>]')
[ "$n403" -gt 0 ] && [ "$n403" = "$n401" ] && [ "$nd" = $((4*n403)) ] && ok "inlet403: $n403 gas inlet records 401 -> 403 with g = $G (2 lines each changed)" || ko "inlet403: rewrite of gas-bc.txt (401 records $n401, 403 records $n403, diff lines $nd)"
sed -n 1,2p "$d/INPUT/gas-bc.txt"
"$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
hygiene inlet403 "$d/run.log" $rc 2; check inlet403 "$d" 403

# --- one403: the same, one coupling step ---
d=$(mk one403 1); cp -p "$T/inlet403/INPUT/gas-bc.txt" "$d/INPUT/gas-bc.txt"
"$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
hygiene one403 "$d/run.log" $rc 1; check one403 "$d" 403

# --- ctrl401: as shipped ---
d=$(mk ctrl401 2)
"$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
hygiene ctrl401 "$d/run.log" $rc 2; check ctrl401 "$d" 401

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
