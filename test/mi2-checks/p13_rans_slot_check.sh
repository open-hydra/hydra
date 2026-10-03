#!/bin/bash
# RANS-slot check: evap-box-sa (MOSE Navier-Stokes + Spalart-Allmaras, IGLOO droplets, closed box), one thread.
#  With a turbulence model the gas carries the SA variable after the energy slot (nprim > np), and the
#  coupler's residual array must leave that slot alone. evap-box-sa/check.py reads the run log and the last
#  gas field:
#   R1  closed-domain mass balance of species 1 closes to <= 1e-12, with >= 5e-9 kg injected
#   R2  energy bookkeeping rel.err <= 1e-7, with >= 1e-2 J handed to the gas
#   R3  the SA variable finite and within 5 % of its initial value everywhere; the eddy viscosity ratio in [0.1, 0.4]
#   R4  log hygiene; the run is Navier-Stokes with the SA closure; the field carries mi_tilde
#   usage: p13_rans_slot_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P13}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/evap-box-sa
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"

d=$T/sa; rm -rf "$d"
"$RUN" "$S" "$d" 1 "$T/sa.log" "$EXE" > /dev/null 2>&1; rc=$?
echo "  evap-box-sa on $EXE: exit=$rc  $(grep -o 'wall=[0-9.]* s' "$T/sa.log")  log $T/sa.log"
python3 -B "$S/check.py" "$T/sa.log" "$d/OUTPUT" "$S/INPUT"
