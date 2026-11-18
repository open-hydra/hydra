#!/bin/bash
# Solidification coherence of ICE and IGLOO on one stream, both solvers standalone (they meet in hydra-MI2):
#  leg A : IGLOO's solid-box stream, a liquid at 2400 K into a 600 K gas -- cooling, nucleation, the freezing
#          plateau, the solid;
#  leg B : a solid at 1500 K into a 3000 K gas, h-fus = 4e5 -- heating to T-melt, the melting plateau, the liquid;
#  leg C : leg A with a liquid heat capacity that varies, c_l = 50 + 0.5 T (a Cp column, the Enthalpy its integral) in
#          both solvers' property tables; the liquid, the event cells, the exit temperature and the total heat against
#          the closed form of the table besides.
# p8_solid_coherence.py writes both cases on one x-spacing (nx = 240 over 0.15 m), this script runs IGLOO once and
# ICE twice per leg (first order, MUSCL), and the .py compares them: (i) the event cells +-1, (ii) T = T-melt on
# the plateau in both, (iii) the heat per unit mass flow in the plateau cells to 1e-6, (iv) the total heat within
# the first-order front error; leg C also (v) the liquid's temperature upstream of nucleation (ICE first order
# against its backward-Euler march of the table and IGLOO; MUSCL against IGLOO) and IGLOO against the closed form.
# ICE must be a build that knows solidification=on (ICE 494668b or later) and a varying cp under it (fbecd78), IGLOO
# one that melts a heated solid (leg B) and freezes on a varying cp (leg C). IGLOO runs on 5 threads in legs A and
# B, on one in the others.
#   usage: p8_solid_coherence_check.sh [scratch-root] [ice-exe] [igloo-exe] [legs]   (default: every leg)
set -u
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
G=$H/test/mi2-checks
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P8}
ICE=${2:-${ICE_BIN:-$H/src/lib/ICE/bin/ICE}}
IGLOO=${3:-${IGLOO_BIN:-$H/src/lib/IGLOO/bin/IGLOO}}
LEGS=${4:-A B C}
mkdir -p "$T"; pass=0; fail=0
echo "  ICE   $ICE ($(md5sum < "$ICE" | cut -c1-8))"
echo "  IGLOO $IGLOO ($(md5sum < "$IGLOO" | cut -c1-8))"
for leg in $LEGS; do
  d=$T/leg$leg; rm -rf "$d"
  if ! python3 "$G/p8_solid_coherence.py" make "$d" $leg; then
    echo "FAIL  leg $leg: case generation"; fail=$((fail+1)); continue
  fi
  case $leg in A|B) nt=5;; *) nt=1;; esac
  (cd "$d/iglo" && ulimit -s unlimited && OMP_NUM_THREADS=$nt KMP_STACKSIZE=100M timeout 1800 "$IGLOO" \
     > run_out.txt 2> run_err.txt); echo "  leg $leg IGLOO exit=$?"
  for o in o1 o2; do
    (cd "$d/ice-$o" && ulimit -s unlimited && OMP_NUM_THREADS=1 KMP_STACKSIZE=100M OMP_STACKSIZE=100M \
       timeout 1800 "$ICE" > log 2> err); echo "  leg $leg ICE $o exit=$?"
  done
  out=$(python3 "$G/p8_solid_coherence.py" check "$d" $leg); echo "$out"
  pass=$((pass + $(echo "$out" | grep -c '^PASS'))); fail=$((fail + $(echo "$out" | grep -c '^FAIL')))
done
echo "== $pass passed, $fail failed"; [ $fail = 0 ]
