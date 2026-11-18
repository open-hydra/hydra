#!/bin/bash
# Coherence of ICE and IGLOO on one stream, both solvers standalone (they meet in hydra-MI2): solidification (legs
# A-C), the drop's internal heating and boiling (legs D, E):
#  leg A : IGLOO's solid-box stream, a liquid at 2400 K into a 600 K gas -- cooling, nucleation, the freezing
#          plateau, the solid;
#  leg B : a solid at 1500 K into a 3000 K gas, h-fus = 4e5 -- heating to T-melt, the melting plateau, the liquid;
#  leg C : leg A with a liquid heat capacity that varies, c_l = 50 + 0.5 T (a Cp column, the Enthalpy its integral) in
#          both solvers' property tables; the liquid, the event cells, the exit temperature and the total heat against
#          the closed form of the table besides;
#  leg D : the drop's internal heating, liquid-conduction = P2T: IGLOO's temp-relax-p2t drop (11.89 um, 2950 kg/m3,
#          c 1250, k-liq = k_g, no mass transfer) from 300 K in a 600 K gas at 10 m/s, against T_s = (5 T + T_g)/6
#          and T - T_g = (T_0 - T_g) exp(-x/L), L = u 1.2 tau_ITC;
#  leg E : boiling = ZGR through the boiling limit: a 20 um water drop 4 K above its boiling temperature at the gas
#          pressure (T_B(rho R T_g) = T_g = 600 K, IGLOO's boiling-zgr-switch gas) at 1 m/s over 0.075 m: the flash to
#          T_B, the switch to CEM, the wet bulb; the flash identity, the switch cell, the wet bulb's closed form and the
#          mass at the exit.
# p8_solid_coherence.py writes both cases on one x-spacing (dx = 0.625 mm), this script runs IGLOO once and
# ICE twice per leg (first order, MUSCL), and the .py compares them: (i) the event cells +-1, (ii) T = T-melt on
# the plateau in both, (iii) the heat per unit mass flow in the plateau cells to 1e-6, (iv) the total heat within
# the first-order front error; leg C also (v) the liquid's temperature upstream of nucleation (ICE first order
# against its backward-Euler march of the table and IGLOO; MUSCL against IGLOO) and IGLOO against the closed form;
# legs D and E read ICE's cells from one probe per cell (T, rho, n and the surface temperature Ts). ICE must be a
# build that knows solidification=on (ICE 494668b or later), a varying cp under it, P2T and ZGR (fbecd78), IGLOO one
# that melts a heated solid (leg B), freezes on a varying cp (leg C) and has P2T and ZGR (legs D, E). IGLOO runs on
# 5 threads in legs A and B, on one in the others.
#   usage: p8_solid_coherence_check.sh [scratch-root] [ice-exe] [igloo-exe] [legs]   (default: every leg)
set -u
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
G=$H/test/mi2-checks
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P8}
ICE=${2:-${ICE_BIN:-$H/src/lib/ICE/bin/ICE}}
IGLOO=${3:-${IGLOO_BIN:-$H/src/lib/IGLOO/bin/IGLOO}}
LEGS=${4:-A B C D E}
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
