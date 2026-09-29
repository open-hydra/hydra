#!/bin/bash
# Floor-guard check (MI2_Coupling.f90 IGLOOxMOSE): the consumed-species guard compares RATES (kg/s).
#  consume : evap-box-ta with [MI2-Coupling] products-A = H2O:1.0, N2:-20.0 (20 kg of N2 consumed per kg of
#            water evaporated). The N2 sink is ~1e-6 kg/s per cell against ~1e-7 kg of N2 per cell: a guard
#            that compares kg/s with kg clamps it (the run then reports the clamp and the N2 mass balance is
#            off by O(1)); the rate guard leaves it alone (no clamp, every species row closes to <= 1e-12).
#  control : products-A = H2O:-1.0 on a box with no water vapour (yH2O = 0, gas at the 1e-20 floor): the
#            guard MUST fire, and the warning is printed once (the duplicate end-of-run print is gone).
#   usage: p6_floor_guard_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P6}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/evap-box-ta
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
mk() { local d=$T/$1; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"; echo "$d"; }
worst_relerr() {  # largest rel.err over the species rows of the closed-domain mass balance
  awk '/closed-domain mass balance/{p=1;next} p&&/^   [0-9]+ +[0-9.E+-]+ +[0-9.E+-]+ +[0-9.E+-]+/{if($4+0>m)m=$4+0} p&&/^ *$/{if(seen)p=0} p{seen=1} END{printf "%.3e", m}' "$1"
}

d=$(mk consume); sed -i 's/^products-A = H2O:1.0$/products-A = H2O:1.0, N2:-20.0/' "$d/input.ini"
grep -q "N2:-20.0" "$d/input.ini" || ko "could not set the consuming stoichiometry"
"$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
n=$(grep -c "floor guard" "$d/run.log"); e=$(worst_relerr "$d/run.log")
echo "  consume: exit=$rc, 'floor guard' lines=$n, worst species rel.err=$e"
[ $rc = 0 ] && ok "consuming case runs" || ko "consuming case exit=$rc"
[ "$n" = 0 ] && ok "no clamp on a sink far below the cell-emptying rate" || ko "the guard clamped a legitimate sink ($n warning line(s))"
python3 -c "import sys; sys.exit(0 if float('$e') <= 1e-12 else 1)" && ok "every species row of the mass balance closes (<= 1e-12)" || ko "mass balance off: worst rel.err $e"

d=$(mk control); sed -i 's/^products-A = H2O:1.0$/products-A = H2O:-1.0/' "$d/input.ini"
"$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
n=$(grep -c "floor guard" "$d/run.log")
echo "  control: exit=$rc, 'floor guard' lines=$n"
[ $rc = 0 ] && [ "$n" -ge 1 ] && ok "the guard still fires on an empty species" || ko "guard did not fire on an empty species (exit=$rc, lines=$n)"
[ "$n" = 1 ] && ok "the clamp warning is printed once" || ko "clamp warning printed $n times"

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
