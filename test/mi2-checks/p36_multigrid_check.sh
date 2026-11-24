#!/bin/bash
# Grid levels under the coupler: MOSE's start on its coarse grids and its multigrid cycle, with IGLOO and ICE. Each leg
# runs a stock case with keys added, at 1 thread, with ATLAS's coarse boundary files (BCB with [ATLAS-Parameters]
# MG-levels = 2): data/evap-box-gas-bc2.txt, data/JPL-gas-bc2.txt (the gas of both JPL nozzles), data/ice-box-part-bc2.txt.
#  coarse    : evap-box, IGLOO's solve first (MOSE-iter = 0), then [MOSE-Multigrid] levels = 2, level2-iter = 300, with
#              [MI2-Coupling] report = true: the run ends cleanly with the gas iterations of the schedule, MOSE steps its
#              coarse level, then its finest; the gas the coupler reports (the finest grid's) moves over the coupling
#              step stepped on the coarse level; the totals of IGLOO's sources restricted to the coarse level equal the
#              finest grid's to round-off.
#  coarse-ice: JPL-Eulerian-20micron (a curved mesh: a coarse cell's volume is not its sub-cells'), 2 levels, the run
#              ends on the coarse level: ICE's per-volume sources restricted with the totals of the finest grid; the
#              field written at the end (MOSE names it after the level the run ends on: gas-field-level2) is the finest
#              grid filled from the coarse level, not the initial one.
#  coarse-igloo: JPL-Lagrangian-20micron, 2 levels, IGLOO's Eulerian coupling: its fields and sources restricted with
#              the totals of the finest grid.
#  ta, relax : a coarse start in a time-accurate run (evap-box-ta), and with relaxation = true: each stops after MOSE's
#              set-up with its refusal, before any MOSE iteration.
#  noiter    : levels = 2, cycle = v, no level1-iter: the run stops with the refusal that names the key.
#  hydra     : [HYDRA-Multigrid] levels = 2: the run stops with its refusal.
#  icelevels : [ICE-Multigrid] levels = 2 (ice-box): the run stops with its refusal, after ICE's set-up.
#  control   : [MOSE-Multigrid] levels = 1 (one coupling step): the run ends cleanly.
#  cycle     : levels = 2, cycle = v, level1-iter set: the run ends cleanly on the finest grid with as many gas
#              iterations as the stock case, which is run beside it.
# These legs: about 40 s at 1 thread on a Xeon W-2265 with the scratch on a local disk.
# Long leg, with P36_LONG=1 (about 52 minutes at 1 thread on a Xeon W-2145, the four runs one after the other): the coarse
# level's own answer on JPL-Lagrangian-20micron, IGLOO's Eulerian sources frozen after one solve on the initial gas
# (before any gas iteration). C: 2 levels, then 12000 iterations on the coarse level and 30000 on the finest; A: no
# levels, 30000 iterations; C0 and A0: the same without IGLOO's solve, the gas alone. The exit Mach number of the
# coarse level (each coarse cell of the last column against its four sub-cells): with the sources, within twice the
# difference of the two levels of the gas alone; a coarse level without the sources (C0 against A, which a coupler
# adding nothing on the coarse level would give) is out of it. After the switch the finest grid ends on A's field
# (max |d| over the variable's scale below 1e-8).
# RED on a build that refuses every coarse start: coarse, coarse-ice, coarse-igloo and the refusals of ta and relax.
# Needs python3 with numpy.
#   usage: p36_multigrid_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P36}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
D=$H/test/mi2-checks/data
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
echo "  exe $EXE ($(md5sum < "$EXE" | cut -c1-8))"
# variant <name> <case> <text appended to input.ini>: a copy of the case with the text added
variant() { local d=$T/$1-src; rm -rf "$d"; mkdir -p "$d"; cp -rp "$H/test/$2/INPUT" "$H/test/$2/MESH" "$H/test/$2/input.ini" "$d/"
            [ -f "$H/test/$2/JPL_chang.dat" ] && cp -p "$H/test/$2/JPL_chang.dat" "$d/"
            printf '\n%b' "$3" >> "$d/input.ini"; echo "$d"; }
# schedule <dir> <N-coupling> <MOSE-iter> <key> <values>: the case's [HYDRA-Parameters] lists replaced
schedule() { sed -i -e "s/^N-coupling *=.*/N-coupling    = $2/" -e "s/^MOSE-iter *=.*/MOSE-iter     = $3/" \
                    -e "s/^$4 *=.*/$4 = $5/" "$1/input.ini"; }
run() { rm -rf "$T/$1"; "$RUN" "$2" "$T/$1" 1 "$T/$1.log" "$EXE" > /dev/null 2>&1; }
clean() { [ "$2" = 0 ] && grep -q 'GAME OVER' "$1" && ! grep -q -e '\[ERROR' -e 'p<0 detected' "$1"; }
stopped() { ! grep -q -e 'GAME OVER' -e 'MOSE Grid Level' "$1" && [ "$2" != 0 ]; }
iters() { sed -n 's/^ *Iterations  *\([0-9][0-9]*\) *$/\1/p' "$1" | tail -1; }
# restricted <log>: the report lines of the restriction, '<n> <worst>': n lines, the largest |fine - coarse| over the
# sum of magnitudes (0 where all three are 0; 1 where only the total is 0)
restricted() { awk '/\[MI2-report\] restricted block/ { n++; f=$(NF-2); c=$(NF-1); a=$NF; d=f-c; if (d<0) d=-d;
                    r = (a > 0) ? d/a : ((d > 0) ? 1 : 0); if (r > w) w = r }
                    END { printf "%d %.3e\n", n, w }' "$1"; }
coarse_lines() { grep -c 'MOSE Grid Level 2' "$1"; }
L2='[MOSE-Multigrid]\nlevels      = 2\nlevel2-iter = 300\nlevel1-iter = 1000000000\n'

# --- coarse -------------------------------------------------------------------------------------------------------
d=$(variant coarse evap-box "$L2[MI2-Coupling]\nreport = true\n")
schedule "$d" '1 4' '0 200' IGLOO-iter '1 1'
cp -p "$D/evap-box-gas-bc2.txt" "$d/INPUT/gas-bc2.txt"
run coarse "$d"; rc=$?; L=$T/coarse.log
clean "$L" $rc && grep -q "MI2: MOSE starts on its coarse grids (levels = 2)" "$L" \
  && ok "coarse: a start on the coarse grid runs to its end, the coupler saying what it couples" || ko "coarse: exit $rc ($L)"
n=$(iters "$L"); n2=$(coarse_lines "$L"); n1=$(grep -c 'MOSE Grid Level 1' "$L")
[ "$n" = 800 ] && [ "$n2" -gt 0 ] && [ "$n1" -gt 0 ] \
  && ok "coarse: the 800 gas iterations of the schedule, on the coarse level ($n2 lines), then on the finest ($n1)" \
  || ko "coarse: ${n:-no} gas iterations of the schedule's 800, $n2 lines of the coarse level, $n1 of the finest"
# the first coupling step of cycle 2 is stepped on the coarse level only (MOSE switches levels after its report, not
# before: a coarse level with no source stands still and its residual reaches 0, which switches early)
s1=$(awk '/HYDRA-MI2: CYCLE 2 /{c=1} c && /MOSE Grid Level 1/{print "fine"; exit}
          c && /writing prolongated solution/{print "switch"; exit} c && /\[MI2-report\] cycle 2 step 1 /{print "report"; exit}' "$L")
e=$(awk '/\[MI2-report\] cycle 2 step 1 /{s=1} s && /\[MI2-report\] gas-energy/{print $3, $4; exit}' "$L")
[ "$s1" = report ] && python3 -c "import sys; a, b = map(float, sys.argv[1:3]); sys.exit(0 if abs(b - a) > 1e-6*abs(a) else 1)" ${e:-0 0} \
  && ok "coarse: the gas handed to IGLOO moves over a coupling step stepped on the coarse level only (energy $e)" \
  || ko "coarse: over the first step of cycle 2 (seen first: ${s1:-nothing}) the gas handed to IGLOO moves by no more than 1e-6 (energy ${e:-no report})"
read -r m w <<< "$(restricted "$L")"
[ "$m" -gt 0 ] && python3 -c "import sys; sys.exit(0 if float('$w') <= 1e-12 else 1)" \
  && ok "coarse: $m restricted totals of IGLOO's sources equal the finest grid's (worst $w of the magnitudes)" \
  || ko "coarse: $m restricted totals, worst $w of the magnitudes (bound 1e-12)"

# --- coarse-ice -----------------------------------------------------------------------------------------------------
d=$(variant coarse-ice JPL-Eulerian-20micron \
    '[MOSE-Multigrid]\nlevels      = 2\nlevel2-iter = 1000000\nlevel1-iter = 1000000000\n[MI2-Coupling]\nreport = true\n')
schedule "$d" '1 1' '10 10' ICE-iter '10 0'
cp -p "$D/JPL-gas-bc2.txt" "$d/INPUT/gas-bc2.txt"
run coarse-ice "$d"; rc=$?; L=$T/coarse-ice.log
clean "$L" $rc && [ "$(iters "$L")" = 20 ] && grep 'MOSE ' "$L" | tail -1 | grep -q 'MOSE Grid Level 2 | Iter = *20 ' \
  && ok "coarse-ice: the run ends cleanly on the coarse level after its 20 gas iterations" \
  || ko "coarse-ice: exit $rc, $(iters "$L") gas iterations, last MOSE line: $(grep 'MOSE ' "$L" | tail -1)"
read -r m w <<< "$(restricted "$L")"
[ "$m" -gt 0 ] && grep -q 'restricted block 1 ice-source' "$L" && python3 -c "import sys; sys.exit(0 if float('$w') <= 1e-12 else 1)" \
  && ok "coarse-ice: $m restricted totals of ICE's sources equal the finest grid's on the curved mesh (worst $w)" \
  || ko "coarse-ice: $m restricted totals, worst $w of the magnitudes (bound 1e-12)"
# the end-of-run field, named by MOSE after the level the run ends on, holds the finest grid filled from it
r=$(python3 -c "
import sys; sys.path.insert(0, '$H/test/mi2-checks'); import numpy as np; from tec_cc import read_cc
a = read_cc('$d/INPUT/gas-ic.tec'); b = read_cc('$T/coarse-ice/OUTPUT/gas-field-level2.tec')
assert a['p'].shape == b['p'].shape
print('%.3e' % max(float(np.abs(a[v] - b[v]).max() / (a[v].max() - a[v].min() or 1.0)) for v in ('u', 'p')))" 2>/dev/null)
python3 -c "import sys; sys.exit(0 if float('${r:-0}') > 1e-3 else 1)" \
  && ok "coarse-ice: the field written at the end is the finest grid filled from the coarse level (max |d|/range $r against the initial one)" \
  || ko "coarse-ice: the field written at the end is not the finest grid moved from the initial one (max |d|/range ${r:-unread} of u and p)"

# --- coarse-igloo ---------------------------------------------------------------------------------------------------
d=$(variant coarse-igloo JPL-Lagrangian-20micron \
    '[MOSE-Multigrid]\nlevels      = 2\nlevel2-iter = 1000000\nlevel1-iter = 1000000000\n[MI2-Coupling]\nreport = true\n')
schedule "$d" '1 1' '10 10' IGLOO-iter '2 0'
cp -p "$D/JPL-gas-bc2.txt" "$d/INPUT/gas-bc2.txt"
run coarse-igloo "$d"; rc=$?; L=$T/coarse-igloo.log
read -r m w <<< "$(restricted "$L")"
clean "$L" $rc && [ "$m" -gt 0 ] && grep -q 'restricted block 1 igloo-euler(1)-mass' "$L" \
  && python3 -c "import sys; sys.exit(0 if float('$w') <= 1e-12 else 1)" \
  && ok "coarse-igloo: $m restricted totals of IGLOO's Eulerian fields and sources equal the finest grid's (worst $w)" \
  || ko "coarse-igloo: exit $rc, $m restricted totals, worst $w of the magnitudes (bound 1e-12)"

# --- ta, relax --------------------------------------------------------------------------------------------------------
d=$(variant ta evap-box-ta "$L2")
cp -p "$D/evap-box-gas-bc2.txt" "$d/INPUT/gas-bc2.txt"
run ta "$d"; rc=$?
grep -q 'starts MOSE on its coarse grids, and the run is time-accurate' "$T/ta.log" && stopped "$T/ta.log" $rc \
  && ok "ta: a coarse start in a time-accurate run is refused (exit $rc)" \
  || ko "ta: no refusal of a coarse start in a time-accurate run (exit $rc, $(coarse_lines "$T/ta.log") coarse lines)"
d=$(variant relax evap-box "$L2")
sed -i 's/^relaxation *=.*/relaxation    = true/' "$d/input.ini"
cp -p "$D/evap-box-gas-bc2.txt" "$d/INPUT/gas-bc2.txt"
run relax "$d"; rc=$?
grep -q 'starts MOSE on its coarse grids, and relaxation = true' "$T/relax.log" && stopped "$T/relax.log" $rc \
  && ok "relax: a coarse start with relaxation = true is refused (exit $rc)" \
  || ko "relax: no refusal of a coarse start with relaxation (exit $rc, $(coarse_lines "$T/relax.log") coarse lines)"

# --- noiter -----------------------------------------------------------------------------------------------------------
d=$(variant noiter evap-box '[MOSE-Multigrid]\nlevels = 2\ncycle  = v\n')
cp -p "$D/evap-box-gas-bc2.txt" "$d/INPUT/gas-bc2.txt"
run noiter "$d"; rc=$?
grep -q '\[ERROR\] MI2: \[MOSE-Multigrid\] level1-iter is not set' "$T/noiter.log" && stopped "$T/noiter.log" $rc \
  && ok "noiter: the cycle without level1-iter is refused (exit $rc)" \
  || ko "noiter: no refusal of the cycle without level1-iter (exit $rc, $(grep -c 'GAME OVER' "$T/noiter.log") GAME OVER)"

# --- hydra ------------------------------------------------------------------------------------------------------------
d=$(variant hydra evap-box '[HYDRA-Multigrid]\nlevels = 2\n')
run hydra "$d"; rc=$?
grep -q 'hydra-MI2 does not support multigrid' "$T/hydra.log" && [ $rc != 0 ] && ok "hydra: [HYDRA-Multigrid] levels = 2 is refused" \
  || ko "hydra: no refusal of [HYDRA-Multigrid] levels = 2 (exit $rc)"

# --- icelevels --------------------------------------------------------------------------------------------------------
d=$(variant icelevels ice-box '[ICE-Multigrid]\nlevels = 2\n')
cp -p "$D/ice-box-part-bc2.txt" "$d/INPUT/part-bc2.txt"
run icelevels "$d"; rc=$?
grep -q '\[ERROR\] MI2: \[ICE-Multigrid\] levels = 2' "$T/icelevels.log" && ! grep -q 'GAME OVER' "$T/icelevels.log" && [ $rc != 0 ] \
  && ok "icelevels: [ICE-Multigrid] levels = 2 is refused (exit $rc)" || ko "icelevels: no refusal of ICE's levels (exit $rc)"

# --- control ----------------------------------------------------------------------------------------------------------
d=$(variant control evap-box '[MOSE-Multigrid]\nlevels = 1\n')
sed -i 's/^N-coupling    = 5$/N-coupling    = 1/' "$d/input.ini"
grep -q '^N-coupling    = 1$' "$d/input.ini" || ko "control: N-coupling could not be cut"
run control "$d"; rc=$?
clean "$T/control.log" $rc && ok "control: [MOSE-Multigrid] levels = 1 runs to its end" || ko "control: exit $rc ($T/control.log)"

# --- cycle ------------------------------------------------------------------------------------------------------------
run stock "$H/test/evap-box"
d=$(variant cycle evap-box '[MOSE-Multigrid]\nlevels      = 2\ncycle       = v\nlevel1-iter = 1000000000\n')
cp -p "$D/evap-box-gas-bc2.txt" "$d/INPUT/gas-bc2.txt"
run cycle "$d"; rc=$?
clean "$T/cycle.log" $rc && grep -q "MI2: MOSE's multigrid cycle on 2 levels" "$T/cycle.log" \
  && ok "cycle: levels = 2 with cycle = v runs to its end, the coupler saying what it couples" \
  || ko "cycle: exit $rc ($T/cycle.log)"
n0=$(iters "$T/stock.log"); n1=$(iters "$T/cycle.log"); n2=$(grep -c 'MOSE Grid Level *[2-9]' "$T/cycle.log")
[ "${n0:-0}" -gt 0 ] && [ "$n1" = "$n0" ] && [ "$n2" = 0 ] \
  && ok "cycle: MOSE steps its finest grid only, with the stock case's gas iterations ($n1)" \
  || ko "cycle: ${n1:-no} gas iterations against the stock case's ${n0:-none}, $n2 lines of a coarse level"

# --- long: the coarse level's own answer --------------------------------------------------------------------------
if [ "${P36_LONG:-0}" = 1 ]; then
  L='[MOSE-Multigrid]\nlevels      = 2\nlevel2-iter = 12000\nlevel1-iter = 1000000000\n'
  for v in A C A0 C0; do
    case $v in
      A)  d=$(variant g3-$v JPL-Lagrangian-20micron '');   k='0 30000'; ig='2 0';;
      C)  d=$(variant g3-$v JPL-Lagrangian-20micron "$L"); k='0 42000'; ig='2 0';;
      A0) d=$(variant g3-$v JPL-Lagrangian-20micron '');   k='0 30000'; ig='0 0';;
      C0) d=$(variant g3-$v JPL-Lagrangian-20micron "$L"); k='0 42000'; ig='0 0';;
    esac
    case $v in C|C0) cp -p "$D/JPL-gas-bc2.txt" "$d/INPUT/gas-bc2.txt";; esac
    schedule "$d" '1 1' "$k" IGLOO-iter "$ig"
    run g3-$v "$d"; rc=$?
    clean "$T/g3-$v.log" $rc && ok "g3-$v: runs to its end ($(iters "$T/g3-$v.log") gas iterations)" || ko "g3-$v: exit $rc"
  done
  P=$H/test/mi2-checks/p36_coarse_phase.py
  ex() { python3 "$P" "$1" "$2" 2>/dev/null | sed -n 's/^exit column: max |M_coarse - M_fine| \([^ ]*\).*/\1/p'; }
  ds=$(ex "$T/g3-C/OUTPUT/gas-field-level2.tec" "$T/g3-A/OUTPUT/gas-field.tec")
  d0=$(ex "$T/g3-C0/OUTPUT/gas-field-level2.tec" "$T/g3-A0/OUTPUT/gas-field.tec")
  dm=$(ex "$T/g3-C0/OUTPUT/gas-field-level2.tec" "$T/g3-A/OUTPUT/gas-field.tec")
  python3 -c "import sys; s, z, m = map(float, sys.argv[1:4]); sys.exit(0 if 0 < s <= 2*z < m else 1)" "${ds:-1}" "${d0:-0}" "${dm:-0}" \
    && ok "g3: exit Mach of the coarse level with IGLOO's restricted sources against the finest grid's: $ds, within twice the two levels' difference of the gas alone ($d0); without the sources: $dm" \
    || ko "g3: exit Mach, coarse against fine: $ds with the sources, $d0 for the gas alone, $dm without the sources (bound twice the gas alone's, which the coarse level without sources must exceed)"
  s=$(python3 "$P" --same "$T/g3-A/OUTPUT/gas-field.tec" "$T/g3-C/OUTPUT/gas-field.tec" 2>/dev/null | sed -n 's/^largest max|d|\/scale \(.*\)/\1/p')
  python3 -c "import sys; sys.exit(0 if float(sys.argv[1]) <= 1e-8 else 1)" "${s:-1}" \
    && ok "g3: after the switch the finest grid ends on the field of the run without levels (max |d|/scale $s)" \
    || ko "g3: the finest grid after the switch ends ${s:-unread} of scale away from the run without levels (bound 1e-8)"
fi

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
