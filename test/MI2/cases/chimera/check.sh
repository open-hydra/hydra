#!/bin/bash
# check.sh <run-dir> <exe>
# The overset-mesh warning (MI2_Report Report_Start): the coupler makes no correction for a chimera overlap, and its
# sums over the domain count it once per block that covers it; a run on such a mesh says so.
#   chimera: MOSE's two-block overset case (its test/checks/connections/chimera-face, case plane: two blocks of 8 x 6
#            cells overlapping by 4 cells in x, face 2 of block 1 and face 1 of block 2 of kind 102), run by hydra-MI2
#            with MOSE alone for one iteration: the run ends and prints the warning once, with the 12 chimera boundary
#            cells of MOSE's own set-up table on 2 blocks;
#   control: the generator's one-block case (wedgeref, no 102 face): no warning.
# MOSE_DIR: the MOSE tree (default src/lib/MOSE).
source "$(dirname "${BASH_SOURCE[0]}")/../../checklib.sh"
M=${MOSE_DIR:-$cases/../../../src/lib/MOSE}/test
GEN=$M/checks/connections/chimera-face/make_case.py
[ -f "$GEN" ] && [ -f "$M/common/Air-viscous/phase.txt" ] || { ko "no MOSE test tree under $M"; finish; }

mk() {  # mk <case>: the generator's case as a hydra-MI2 case with MOSE alone
  local d=$run/$1; rm -rf "$d"; mkdir -p "$d/INPUT"
  ( cd "$d" && python3 -B "$GEN" case "$1" ) > /dev/null || return 1
  cp -p "$M"/common/Air-viscous/{phase.txt,thermo.dat,transport.dat} "$d/INPUT/"
  cat > "$d/input.ini" <<'INI'
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
INI
  echo "$d"
}
WARN='\[WARNING\] MI2: an overset mesh'

d=$(mk plane) || ko "could not write the chimera case"
runin "$d"
n=$(grep -c "$WARN" "$d/run.log")
mose=$(awk '$1 == "Chimera" {print $2; exit}' "$d/run.log")
said=$(sed -n "s/.*$WARN: \([0-9]*\) chimera (102) boundary cells on \([0-9]*\) block.*/\1 \2/p" "$d/run.log")
expect "chimera: the run ends" ended "$d/run.log"
expect "chimera: the warning is printed once ($n)" [ "$n" = 1 ]
expect "chimera: it counts the 12 chimera cells of MOSE's table (${mose:-none}) on 2 blocks ('${said:-none}')" \
  [ "${mose:-0}" = 12 -a "$said" = "12 2" ]

d=$(mk wedgeref) || ko "could not write the control case"
runin "$d"
expect "control: the run ends" ended "$d/run.log"
expect "control: no warning without a chimera face" [ "$(grep -c "$WARN" "$d/run.log")" = 0 ]
finish
