#!/bin/bash
# Hot re-read of input.ini ([HYDRA-IO] ini-diter), on evap-box-ta with a re-read at every second coupling step:
#  reread : the run ends (exit 0, GAME OVER) without an [ERROR] line -- the keys of the solvers MI2 does not run
#           (FUSS-iter, the multigrid levels) are not required of an MI2 input;
#  inert  : re-reading an unchanged file changes nothing, OUTPUT byte-identical to the run of the case as shipped.
#   usage: p9_ini_reread_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P9}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/evap-box-ta
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }

"$RUN" "$S" "$T/shipped" 1 "$T/shipped.log" "$EXE" > /dev/null 2>&1

d=$T/case; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
grep -q '^ini-diter' "$d/input.ini" || { echo "FAIL  the case has no ini-diter line to edit"; exit 2; }
sed -i 's/^ini-diter *=.*/ini-diter = 2/' "$d/input.ini"
"$RUN" "$d" "$T/reread" 1 "$T/reread.log" "$EXE" > /dev/null 2>&1; rc=$?
[ $rc = 0 ] && grep -q "GAME OVER" "$T/reread.log" && ok "the run with the re-read ends (exit 0, GAME OVER)" || ko "the run with the re-read: exit=$rc"
if grep -q '\[ERROR\]' "$T/reread.log"; then ko "an [ERROR] line: $(grep -m1 '\[ERROR\]' "$T/reread.log")"; else ok "no [ERROR] line"; fi
diff -rq "$T/shipped/OUTPUT" "$T/reread/OUTPUT" > /dev/null && ok "OUTPUT byte-identical to the run without re-read" || ko "OUTPUT differs from the run without re-read"
echo "== $pass passed, $fail failed"; [ $fail = 0 ]
