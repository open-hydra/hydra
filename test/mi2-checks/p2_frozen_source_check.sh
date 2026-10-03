#!/bin/bash
# P2 check: frozen-source notice cadence, plus inertness when nothing is frozen.
#  frozen : evap-box-ta, 2 cycles: MOSE-iter 50 50, N-coupling 3 3, IGLOO-iter 1 0, [MOSE-IO] shell-diter 100
#           -> cycle 2 keeps IGLOO's cycle-1 sources; notices expected at cycle-2 steps 1 and 3 only.
#              Cumulative MOSE iters at step start are 101 / 151 / 201, NOT 150 / 200 / 250: the very first
#              coupling step on the quiescent IC converges at MOSE iteration 1 (zero residual -> TODO = 3 ->
#              the remaining 49 iterations of that step are skipped), so cycle 1 = 1 + 50 + 50.
#              step 1 -> first of the cycle; step 2 -> 151/100 == 101/100, silent; step 3 -> 201/100 > 1.
#  inert  : evap-box-ta as shipped -> no notice, OUTPUT byte-identical to the case's reference OUTPUT/
#   usage: p2_frozen_source_check.sh [scratch-root] [exe]   (REF=<OUTPUT dir of the reference run> for the byte check)
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P2}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/evap-box-ta
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"
pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }

# --- frozen ---
d=$T/frozen; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
sed -i -e 's/^N-coupling *=.*/N-coupling    = 3 3/' -e 's/^MOSE-iter *=.*/MOSE-iter     = 50 50/' \
       -e 's/^IGLOO-iter *=.*/IGLOO-iter    = 1 0/' -e 's/^shell-diter *=.*/shell-diter   = 100/' "$d/input.ini"
"$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
n=$(grep -c "IGLOO source FROZEN" "$d/run.log")
steps=$(grep -o "\[cycle [0-9]*, step [0-9]*, MOSE iter [0-9]*\]" "$d/run.log" | tr '\n' ' ')
echo "  frozen notices: $n  at: $steps  (exit=$rc)"
[ $rc = 0 ] && [ "$n" = 2 ] && [ "$steps" = "[cycle 2, step 1, MOSE iter 101] [cycle 2, step 3, MOSE iter 201] " ] && ok "frozen-source notice at first step + shell-diter cadence" || ko "frozen-source notice cadence"
grep -q "produced in cycle 1, coupling step 3" "$d/run.log" && ok "notice names the producing (cycle, step)" || ko "notice provenance"
# the frozen sources must still be applied: the mass audit of cycle 2 must keep growing
inj=$(grep -A3 "closed-domain mass balance" "$d/run.log" | grep "^   1  " | awk '{print $2}')
echo "  injected total with a frozen cycle: $inj"

# --- inert ---
d=$T/inert; rm -rf "$d"
"$RUN" "$S" "$d" 1 "$T/inert.log" "$EXE" > /dev/null 2>&1; rc=$?
if grep -q "FROZEN" "$T/inert.log"; then ko "no notice when nothing is frozen"; else ok "no notice when nothing is frozen"; fi
REF=${REF:-$S/OUTPUT}   # the case's maintained reference run (1 thread)
if [ -d "$REF" ]; then
  diff -rq "$REF" "$d/OUTPUT" > /dev/null && ok "OUTPUT byte-identical to the reference run" || ko "OUTPUT differs from the reference run ($REF)"
else
  echo "SKIP  no reference run at $REF: the byte comparison was not made"
fi
echo "== $pass passed, $fail failed"; [ $fail = 0 ]
