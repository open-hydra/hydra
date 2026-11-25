#!/bin/bash
# check.sh <run-dir> <exe>
# ice-burnout-steady: MOSE and ICE steady, the drops burning out in the nozzle. In pseudo-time every transfer between
# families and solvers is handed as a rate, so the converged field does not depend on the local time steps: the deck
# runs as given (ICE cfl 0.2) and with ICE cfl 0.1, and every cell-centred variable of every OUTPUT field agrees within
# 1e-12 of its peak (ten times the floor: the deck run 40 against 50 cycles differs by 3e-14). The coupler reports the
# burn-out handed to the gas.
source "$(dirname "${BASH_SOURCE[0]}")/../../checklib.sh"
src=$cases/ice-burnout-steady
RUN=$cases/../run_case.sh
claim() { if [ "$1" = 0 ]; then ok "$2"; else ko "$2"; fi; }

rm -rf "$run"; mkdir -p "$run/variant"
cp -rL "$src/INPUT" "$src/MESH" "$src/input.ini" "$run/variant/"
sed -i 's/^cfl                  = 0.2$/cfl                  = 0.1/' "$run/variant/input.ini"
grep -q '^cfl                  = 0.1$' "$run/variant/input.ini"; claim $? "the variant halves ICE's cfl"
"$RUN" "$src" "$run/a" 1 "$run/a.log" "$exe" > /dev/null 2>&1 &
"$RUN" "$run/variant" "$run/b" 1 "$run/b.log" "$exe" > /dev/null 2>&1 &
wait
for r in a b; do grep -q 'GAME OVER' "$run/$r.log"; claim $? "run $r ends (GAME OVER)"; done
r=$(grep -a 'MI2-report\] ice-to-gas-mass' "$run/a.log" | tail -1 | awk '{print $3}')
python3 -c "import sys; sys.exit(0 if float('${r:-0}') > 0 else 1)"; claim $? "the coupler hands the gas a burn of ${r:-none} kg/s"
out=$(python3 -B "$src/step_cmp.py" "$run/a" "$run/b" 1e-12); rc=$?
echo "$out" | sed 's/^/      /'
claim $rc "converged fields independent of ICE's step: $(echo "$out" | tail -1) (tol 1e-12)"
finish
