#!/bin/bash
# Overset-mesh warning check (MI2_Report.f90 Report_Start): the coupler makes no correction for a chimera overlap,
# and its sums over the domain count it once per block that covers it; a run on such a mesh must say so.
#  chimera : MOSE's two-block overset case (its test/checks/connections/chimera-face, case plane: two blocks of
#            8 x 6 cells overlapping by 4 cells in x, face 2 of block 1 and face 1 of block 2 of kind 102), run by
#            hydra-MI2 with MOSE alone for one iteration: the run ends and prints the warning once, with the number
#            of chimera boundary cells MOSE's own set-up table gives (12) on 2 blocks.
#  control : the same generator's one-block case (wedgeref, no 102 face): the run ends without the warning.
# RED on a build without the warning: the chimera run ends with no warning line.
#   usage: p35_chimera_warning_check.sh [scratch-root] [exe]      (MOSE_DIR: the MOSE tree, default src/lib/MOSE)
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P35}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
M=${MOSE_DIR:-$H/src/lib/MOSE}/test
GEN=$M/checks/connections/chimera-face/make_case.py
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
[ -f "$GEN" ] && [ -f "$M/common/Air-viscous/phase.txt" ] || { echo "FAIL  no MOSE test tree under $M"; exit 1; }

mk() {  # mk <case>: the generator's case as a hydra-MI2 case with MOSE alone
  local d=$T/$1; rm -rf "$d"; mkdir -p "$d/INPUT"
  ( cd "$d" && python3 -B "$GEN" case "$1" ) || return 1
  cp -p "$M"/common/Air-viscous/{phase.txt,thermo.dat,transport.dat} "$d/INPUT/"
  cat > "$d/input.ini" <<'EOF'
[HYDRA-Parameters]
N-coupling    = 1
MOSE-iter     = 1
time-accurate = false

[MOSE-Parameters]
iter-threshold = 1
newrun = true

[MOSE-Physics]
equations = navier-stokes

[MOSE-Numerics]
cfl = 0.5
time-accurate = false
time-scheme = euler
integration-variables = cons
space-reconstruction = MUSCL
flux-limiter = vanleer
riemann-solver = HLLC

[MOSE-IO]
sol-format = tecplot ascii
EOF
  echo "$d"
}
ended() { grep -q "GAME OVER" "$1" && ! grep -q -i "error" "$1"; }
WARN="\[WARNING\] MI2: an overset mesh"

d=$(mk plane) || ko "could not write the chimera case"
"$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1
n=$(grep -c "$WARN" "$d/run.log")
mose=$(awk '$1=="Chimera"{print $2; exit}' "$d/run.log")
said=$(sed -n "s/.*$WARN: \([0-9]*\) chimera (102) boundary cells on \([0-9]*\) block.*/\1 \2/p" "$d/run.log")
echo "  chimera: warning lines=$n, MOSE's table: ${mose:-none} chimera cells, the warning: '${said:-none}' (cells blocks)"
ended "$d/run.log" && ok "the chimera case runs to its end" || ko "the chimera case did not end cleanly (see $d/run.log)"
[ "$n" = 1 ] && ok "the warning is printed once on an overset mesh" || ko "the warning is printed $n times on an overset mesh"
[ "${mose:-0}" = 12 ] && [ "$said" = "12 2" ] && ok "it counts the 12 chimera cells of MOSE's table on 2 blocks" \
  || ko "the warning's counts '${said:-none}' are not 12 cells (MOSE's table: ${mose:-none}) on 2 blocks"

d=$(mk wedgeref) || ko "could not write the control case"
"$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1
n=$(grep -c "$WARN" "$d/run.log")
echo "  control: warning lines=$n"
ended "$d/run.log" && ok "the control case runs to its end" || ko "the control case did not end cleanly (see $d/run.log)"
[ "$n" = 0 ] && ok "no warning on a mesh without a chimera face" || ko "the warning is printed on a mesh without a chimera face"

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
