#!/bin/bash
# check.sh <run-dir> <exe>
# Grid levels under the coupler: MOSE's start on its coarse grids and its multigrid cycle, with IGLOO and ICE. Each leg
# is a case of test/MI2 with keys added, at 1 thread, with ATLAS's coarse boundary files (BCB with [ATLAS-Parameters]
# MG-levels = 2) of this folder: <case>-gas-bc2.txt, <case>-part-bc2.txt. The coarse legs run with [MI2-Coupling]
# report = true, whose lines they judge: 'restricted' (at each rebuild of the coarse level's sources, per block and
# row of what a solver hands the gas: the finest grid's total, the coarse level's, the sum of magnitudes; totals only)
# and 'residual' (at each rebuild, per row of the residual, what the per-stage routine of IGLOO, ICE or the law
# channels subtracts on a copy of the finest grid and of the coarse level, the floor guard off). Every judged quantity
# must have a non-zero magnitude.
#  coarse    : validation/evap-box, IGLOO's solve first (MOSE-iter = 0), then levels = 2, level2-iter = 300: the run
#              ends with the gas iterations of the schedule, on the coarse level, then the finest; the gas handed to
#              IGLOO (the finest grid's) moves over the coupling step stepped on the coarse level; IGLOO's restricted
#              totals and residual rows equal the finest grid's to 1e-12; its sources are restricted again at each
#              coupling step that starts on the coarse level: ceil(level2-iter / MOSE-iter) times.
#  coarse-ice: validation/JPL-Eulerian-20micron (a curved mesh: a coarse cell's volume is not its sub-cells'), the run
#              ending on the coarse level: ICE's totals and rows as above; the field written at the end
#              (gas-field-level2) is the finest grid filled from the coarse level, not the initial one.
#  coarse-igloo: validation/JPL-Lagrangian-20micron, IGLOO's Eulerian coupling: its restricted fields as above; the
#              inflow handed to IGLOO after the coarse block positive; the Eulerian rows non-zero (they depend on the
#              level's gas state: not compared).
#  coarse-gen: ch-steady (igloo-al in pseudo-time: IGLOO's aluminium generates alumina in ICE, three law channels),
#              2 coupling steps: IGLOO's rows, the generated mass withheld from the gas included, and the channels'
#              rows equal the finest grid's to 1e-12 (ICE's drops are born at the gas's state: its rows are zero here).
#  coarse-removal: ice-burnout-steady with K-burn x100 and ICE time-accurate (first order, cfl 0.4), MOSE steady:
#              ICE alone, then MOSE alone on the coarse level: ICE's totals and rows, its removals included (their
#              rate in the report non-zero), equal the finest grid's to 1e-12.
#  ta, relax, restart: a coarse start in a time-accurate run (igloo-evap), with relaxation = true, and with
#              newrun = false on the field of the stock evap-box run: each stops after MOSE's set-up with its refusal.
#  noiter, hydra, icelevels: levels = 2, cycle = v without level1-iter; [HYDRA-Multigrid] levels = 2;
#              [ICE-Multigrid] levels = 2 (ice-exchange): each refused.
#  control   : [MOSE-Multigrid] levels = 1 (one coupling step): the run ends cleanly.
#  cycle     : levels = 2, cycle = v, level1-iter set: the finest grid only, with the stock case's gas iterations.
# Not seen: the refusals after input.ini is read again during a run; the per-cell parent map (totals only); the
# Eulerian branch's volume; the floor guard on the coarse level.
# Long leg, MULTIGRID_LONG=1 (about 52 minutes): the coarse level's own answer on JPL-Lagrangian-20micron, IGLOO's
# Eulerian sources frozen after one solve on the initial gas. coarse-src: 2 levels, then 12000 iterations on the
# coarse level and 30000 on the finest; fine-src: no levels, 30000 iterations; coarse-gas and fine-gas: the same
# without IGLOO, the gas alone. coarse-answer: the exit Mach number of the coarse level with the sources within twice
# the two levels' difference of the gas alone, and the coarse level without the sources outside it;
# fine-after-switch: the finest grid ends on fine-src's field (max |d| over the variable's scale below 1e-8).
source "$(dirname "${BASH_SOURCE[0]}")/../../checklib.sh"
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); M2=$(cd "$here/../.." && pwd); V=$M2/validation
export OMP_NUM_THREADS=1

# mg <name> <base dir> <text appended to input.ini> [coarse file -> INPUT name]...: a variant
mg() {
  local d; d=$(variant "$1" "$2"); printf '\n%b' "$3" >> "$d/input.ini"; shift 3
  while [ $# -ge 2 ]; do cp -p "$here/$1" "$d/INPUT/$2"; shift 2; done
  echo "$d"
}
# schedule <dir> <N-coupling> <MOSE-iter> <key> <values>
schedule() { sed -i -e "s/^N-coupling *=.*/N-coupling    = $2/" -e "s/^MOSE-iter *=.*/MOSE-iter     = $3/" \
                    -e "s/^$4 *=.*/$4 = $5/" "$1/input.ini"; }
clean() { [ "$2" = 0 ] && ended "$1" && ! grep -q 'p<0 detected' "$1"; }
stopped() { ! grep -q -e 'GAME OVER' -e 'MOSE Grid Level' "$1" && [ "$2" != 0 ]; }
iters() { sed -n 's/^ *Iterations  *\([0-9][0-9]*\) *$/\1/p' "$1" | tail -1; }
# rows <log> <pattern>: the report lines matching the pattern, '<n> <nonzero> <worst>': n lines, those of a non-zero
# magnitude, the largest |fine - coarse| over the sum of magnitudes (0 where all three are 0; 1 where only that is)
rows() { awk -v p="$2" '$0 ~ p { n++; f=$(NF-2); c=$(NF-1); a=$NF; d=f-c; if (d<0) d=-d; if (a > 0) z++
                        r = (a > 0) ? d/a : ((d > 0) ? 1 : 0); if (r > w) w = r } END { printf "%d %d %.3e\n", n, z, w }' "$1"; }
restricted() { rows "$1" '[[]MI2-report[]] restricted block'; }
residual() { rows "$1" "[[]MI2-report[]] residual level [0-9]+ $2 row"; }
# magnitude <log> <quantity>: the sum of magnitudes of the first restricted line of block 1 for it (0 if none)
magnitude() { awk -v q="$2" '/\[MI2-report\] restricted block 1 / && $5 == q { print $NF; f=1; exit } END { if (!f) print 0 }' "$1"; }
le() { python3 -c "import sys; sys.exit(0 if float(sys.argv[1]) <= float(sys.argv[2]) else 1)" "$1" "$2"; }
gt0() { python3 -c "import sys; sys.exit(0 if float(sys.argv[1]) > 0 else 1)" "$1"; }
# same <label> <log> <routine>: the routine's rows on the coarse level equal the finest grid's
same() {
  local m z w; read -r m z w <<< "$(residual "$2" "$3")"
  expect "$1: what $3's per-stage routine subtracts on the coarse level equals what it subtracts on the finest grid, $m rows ($z non-zero, worst $w; bound 1e-12)" \
    eval "[ $m -gt 0 ] && [ $z -gt 0 ] && le $w 1e-12"
}
# totals <label> <log> <quantity of a non-zero magnitude>: the restricted totals equal the finest grid's
totals() {
  local m z w a; read -r m z w <<< "$(restricted "$2")"; a=$(magnitude "$2" "$3")
  expect "$1: $m restricted totals equal the finest grid's (worst $w; $z non-zero, $3's magnitude $a; bound 1e-12)" \
    eval "[ $m -gt 0 ] && [ $z -gt 0 ] && gt0 $a && le $w 1e-12"
}
coarse_lines() { grep -c 'MOSE Grid Level 2' "$1"; }
L2='[MOSE-Multigrid]\nlevels      = 2\nlevel2-iter = 300\nlevel1-iter = 1000000000\n'
LONG='[MOSE-Multigrid]\nlevels      = 2\nlevel2-iter = 1000000\nlevel1-iter = 1000000000\n'
R='[MI2-Coupling]\nreport = true\n'
EB=evap-gas-bc2.txt

# --- coarse ---------------------------------------------------------------------------------------------------------
d=$(mg coarse "$V/evap-box" "$L2$R" $EB gas-bc2.txt); mi=200; schedule "$d" '1 4' "0 $mi" IGLOO-iter '1 1'
runin "$d"; rc=$?; L=$d/run.log
expect "coarse: a start on the coarse grid runs to its end, the coupler saying what it couples" \
  eval "clean $L $rc && grep -q 'MI2: MOSE starts on its coarse grids (levels = 2)' $L"
n=$(iters "$L"); n2=$(coarse_lines "$L"); n1=$(grep -c 'MOSE Grid Level 1' "$L")
expect "coarse: the 800 gas iterations of the schedule, on the coarse level ($n2 lines), then on the finest ($n1)" \
  eval "[ '$n' = 800 ] && [ $n2 -gt 0 ] && [ $n1 -gt 0 ]"
# the first coupling step of cycle 2 is stepped on the coarse level only (MOSE switches levels after its report)
s1=$(awk '/HYDRA-MI2: CYCLE 2 /{c=1} c && /MOSE Grid Level 1/{print "fine"; exit}
          c && /writing prolongated solution/{print "switch"; exit} c && /\[MI2-report\] cycle 2 step 1 /{print "report"; exit}' "$L")
e=$(awk '/\[MI2-report\] cycle 2 step 1 /{s=1} s && /\[MI2-report\] gas-energy/{print $3, $4; exit}' "$L")
expect "coarse: the gas handed to IGLOO moves over a coupling step stepped on the coarse level only (energy ${e:-none})" \
  eval "[ '$s1' = report ] && python3 -c \"import sys; a, b = map(float, sys.argv[1:3]); sys.exit(0 if abs(b - a) > 1e-6*abs(a) else 1)\" ${e:-0 0}"
totals coarse "$L" igloo-cell-energy
same coarse "$L" igloo-source
nb=$(( (300 + mi - 1) / mi )); nr=$(grep -c '\[MI2-report\] restricted IGLOO to level 2' "$L")
expect "coarse: IGLOO's sources restricted again at each of the $nb coupling steps that start on the coarse level ($nr)" [ "$nr" = "$nb" ]

# --- coarse-ice -----------------------------------------------------------------------------------------------------
d=$(mg coarse-ice "$V/JPL-Eulerian-20micron" "$LONG$R" JPL-gas-bc2.txt gas-bc2.txt)
schedule "$d" '1 1' '10 10' ICE-iter '10 0'
runin "$d"; rc=$?; L=$d/run.log
expect "coarse-ice: the run ends cleanly on the coarse level after its 20 gas iterations" \
  eval "clean $L $rc && [ '$(iters "$L")' = 20 ] && grep 'MOSE ' $L | tail -1 | grep -q 'MOSE Grid Level 2 | Iter = *20 '"
totals coarse-ice "$L" 'ice(1)-mom(1)'
same coarse-ice "$L" ice
r=$(python3 -c "
import sys; sys.path.insert(0, '$M2'); import numpy as np; from tec_cc import read_cc
a = read_cc('$d/INPUT/gas-ic.tec'); b = read_cc('$d/OUTPUT/gas-field-level2.tec')
assert a['p'].shape == b['p'].shape
print('%.3e' % max(float(np.abs(a[v] - b[v]).max() / (a[v].max() - a[v].min() or 1.0)) for v in ('u', 'p')))" 2>/dev/null)
expect "coarse-ice: the field written at the end is the finest grid filled from the coarse level (max |d|/range ${r:-unread} against the initial one)" \
  python3 -c "import sys; sys.exit(0 if float('${r:-0}') > 1e-3 else 1)"

# --- coarse-igloo ---------------------------------------------------------------------------------------------------
d=$(mg coarse-igloo "$V/JPL-Lagrangian-20micron" "$LONG$R" JPL-gas-bc2.txt gas-bc2.txt)
schedule "$d" '1 1' '10 10' IGLOO-iter '2 0'
runin "$d"; rc=$?; L=$d/run.log
expect "coarse-igloo: the run ends" clean "$L" $rc
totals coarse-igloo "$L" 'igloo-euler(1)-mass'
q=$(awk '/\[MI2-report\] igloo-inlet-mdot/ { n++; if ($3 > 0) p++ } END { printf "%d %d", n, p }' "$L"); read -r n p <<< "$q"
expect "coarse-igloo: the inflow handed to IGLOO after the coarse block is positive ($p of $n lines)" \
  eval "[ $n -gt 0 ] && [ $p = $n ]"
read -r m z w <<< "$(residual "$L" igloo-euler)"
expect "coarse-igloo: IGLOO's Eulerian per-stage routine subtracts on the coarse level ($m rows, $z non-zero; against the finest grid $w, not judged)" \
  eval "[ $m -gt 0 ] && [ $z -gt 0 ]"

# --- coarse-gen -----------------------------------------------------------------------------------------------------
d=$(mg coarse-gen "$cases/ch-steady" "$LONG$R" igloo-al-gas-bc2.txt gas-bc2.txt)
schedule "$d" 2 10 ICE-iter 10
runin "$d"; rc=$?; L=$d/run.log
expect "coarse-gen: a coarse start with IGLOO's generated mass and the law channels runs to its end" clean "$L" $rc
totals coarse-gen "$L" igloo-cell-energy
same coarse-gen "$L" igloo-source
same coarse-gen "$L" channels

# --- coarse-removal -------------------------------------------------------------------------------------------------
d=$(mg coarse-removal "$cases/ice-burnout-steady" "$LONG$R" ice-burnout-steady-gas-bc2.txt gas-bc2.txt)
schedule "$d" '1 1' '0 10' ICE-iter '20 0'
sed -i -e 's/^time-scheme          = RK3$/time-scheme          = euler/' \
       -e 's/^space-reconstruction = MUSCL$/space-reconstruction = first-order/' -e '/^flux-limiter         = vanleer$/d' \
       -e 's/^time-accurate        = false$/time-accurate        = true/' -e 's/^cfl                  = 0.2$/cfl                  = 0.4/' \
       "$d/input.ini"
sed -i 's/ K-burn=1.39e-06 / K-burn=1.39e-04 /' "$d/INPUT/part-models.txt"
runin "$d"; rc=$?; L=$d/run.log
rm_=$(grep -a '\[MI2-report\] ice-removed-rate' "$L" | tail -1 | awk '{print $3}')
expect "coarse-removal: ICE alone, then MOSE alone on the coarse level: 10 gas iterations, the removals' rate in force ${rm_:-none} kg/s" \
  eval "clean $L $rc && [ '$(iters "$L")' = 10 ] && gt0 ${rm_:-0}"
totals coarse-removal "$L" 'ice(1)-rho(1)'
same coarse-removal "$L" ice

# --- ta, relax, restart ---------------------------------------------------------------------------------------------
d=$(mg ta "$cases/igloo-evap" "$L2" $EB gas-bc2.txt); runin "$d"; rc=$?
expect "ta: a coarse start in a time-accurate run is refused (exit $rc)" \
  eval "grep -q 'starts MOSE on its coarse grids, and the run is time-accurate' $d/run.log && stopped $d/run.log $rc"
d=$(mg relax "$V/evap-box" "$L2" $EB gas-bc2.txt); sed -i 's/^relaxation *=.*/relaxation    = true/' "$d/input.ini"
runin "$d"; rc=$?
expect "relax: a coarse start with relaxation = true is refused (exit $rc)" \
  eval "grep -q 'starts MOSE on its coarse grids, and relaxation = true' $d/run.log && stopped $d/run.log $rc"

# --- noiter, hydra, icelevels ---------------------------------------------------------------------------------------
d=$(mg noiter "$V/evap-box" '[MOSE-Multigrid]\nlevels = 2\ncycle  = v\n' $EB gas-bc2.txt); runin "$d"; rc=$?
expect "noiter: the cycle without level1-iter is refused (exit $rc)" \
  eval "grep -q '\[ERROR\] MI2: \[MOSE-Multigrid\] level1-iter is not set' $d/run.log && stopped $d/run.log $rc"
d=$(mg hydra "$V/evap-box" '[HYDRA-Multigrid]\nlevels = 2\n'); runin "$d"; rc=$?
expect "hydra: [HYDRA-Multigrid] levels = 2 is refused (exit $rc)" \
  eval "[ $rc != 0 ] && grep -q 'hydra-MI2 does not support multigrid' $d/run.log"
d=$(mg icelevels "$cases/ice-exchange" '[ICE-Multigrid]\nlevels = 2\n' ice-exchange-part-bc2.txt part-bc2.txt)
runin "$d"; rc=$?
expect "icelevels: [ICE-Multigrid] levels = 2 is refused (exit $rc)" \
  eval "[ $rc != 0 ] && grep -q '\[ERROR\] MI2: \[ICE-Multigrid\] levels = 2' $d/run.log && ! grep -q 'GAME OVER' $d/run.log"

# --- control, cycle -------------------------------------------------------------------------------------------------
d=$(mg control "$V/evap-box" '[MOSE-Multigrid]\nlevels = 1\n'); sed -i 's/^N-coupling *=.*/N-coupling    = 1/' "$d/input.ini"
runin "$d"; rc=$?
expect "control: [MOSE-Multigrid] levels = 1 runs to its end (exit $rc)" clean "$d/run.log" $rc
s=$(variant stock "$V/evap-box"); runin "$s"
d=$(mg cycle "$V/evap-box" '[MOSE-Multigrid]\nlevels      = 2\ncycle       = v\nlevel1-iter = 1000000000\n' $EB gas-bc2.txt)
runin "$d"; rc=$?
expect "cycle: levels = 2 with cycle = v runs to its end, the coupler saying what it couples" \
  eval "clean $d/run.log $rc && grep -q \"MI2: MOSE's multigrid cycle on 2 levels\" $d/run.log"
n0=$(iters "$s/run.log"); n1=$(iters "$d/run.log"); n2=$(grep -c 'MOSE Grid Level *[2-9]' "$d/run.log")
expect "cycle: the finest grid only, the stock case's gas iterations (${n1:-none} vs ${n0:-none})" \
  eval "[ '${n0:-0}' -gt 0 ] && [ '$n1' = '$n0' ] && [ $n2 = 0 ]"

# --- restart: a coarse start with newrun = false, on the stock run's field -------------------------------------------
d=$run/restart; rm -rf "$d"; mkdir -p "$d"; cp -rp "$s/INPUT" "$s/MESH" "$s/input.ini" "$s/OUTPUT" "$d/"
printf '\n%b' "$L2" >> "$d/input.ini"; sed -i '/^\[MOSE-Parameters\]/a newrun = false' "$d/input.ini"
cp -p "$here/$EB" "$d/INPUT/gas-bc2.txt"
( cd "$d" && ls INPUT/*phase.txt > filelist.txt && ulimit -s unlimited && KMP_STACKSIZE=100M "$exe" > run.log 2>&1 ); rc=$?
expect "restart: a coarse start that restarts MOSE is refused (exit $rc)" \
  eval "grep -q 'starts MOSE on its coarse grids, and newrun = false' $d/run.log && stopped $d/run.log $rc"

# --- long: the coarse level's own answer ----------------------------------------------------------------------------
if [ "${MULTIGRID_LONG:-0}" = 1 ]; then
  LL='[MOSE-Multigrid]\nlevels      = 2\nlevel2-iter = 12000\nlevel1-iter = 1000000000\n'
  for v in fine-src coarse-src fine-gas coarse-gas; do
    case $v in
      fine-src)   d=$(mg $v "$V/JPL-Lagrangian-20micron" '');  k='0 30000'; ig='2 0';;
      coarse-src) d=$(mg $v "$V/JPL-Lagrangian-20micron" "$LL" JPL-gas-bc2.txt gas-bc2.txt); k='0 42000'; ig='2 0';;
      fine-gas)   d=$(mg $v "$V/JPL-Lagrangian-20micron" '');  k='0 30000'; ig='0 0';;
      coarse-gas) d=$(mg $v "$V/JPL-Lagrangian-20micron" "$LL" JPL-gas-bc2.txt gas-bc2.txt); k='0 42000'; ig='0 0';;
    esac
    schedule "$d" '1 1' "$k" IGLOO-iter "$ig"
    runin "$d"; rc=$?
    expect "$v: runs to its end ($(iters "$d/run.log") gas iterations)" clean "$d/run.log" $rc
  done
  P=$here/coarse_phase.py
  ex() { python3 "$P" "$1" "$2" 2>/dev/null | sed -n 's/^exit column: max |M_coarse - M_fine| \([^ ]*\).*/\1/p'; }
  ds=$(ex "$run/coarse-src/OUTPUT/gas-field-level2.tec" "$run/fine-src/OUTPUT/gas-field.tec")
  d0=$(ex "$run/coarse-gas/OUTPUT/gas-field-level2.tec" "$run/fine-gas/OUTPUT/gas-field.tec")
  dm=$(ex "$run/coarse-gas/OUTPUT/gas-field-level2.tec" "$run/fine-src/OUTPUT/gas-field.tec")
  expect "coarse-answer: exit Mach of the coarse level with IGLOO's restricted sources against the finest grid's ${ds:-?}, within twice the two levels' difference of the gas alone (${d0:-?}); without the sources ${dm:-?}" \
    python3 -c "import sys; s, z, m = map(float, sys.argv[1:4]); sys.exit(0 if 0 < s <= 2*z < m else 1)" "${ds:-1}" "${d0:-0}" "${dm:-0}"
  sw=$(python3 "$P" --same "$run/fine-src/OUTPUT/gas-field.tec" "$run/coarse-src/OUTPUT/gas-field.tec" 2>/dev/null | sed -n 's/^largest max|d|\/scale \(.*\)/\1/p')
  expect "fine-after-switch: the finest grid ends on the field of the run without levels (max |d|/scale ${sw:-unread}; bound 1e-8)" \
    python3 -c "import sys; sys.exit(0 if float(sys.argv[1]) <= 1e-8 else 1)" "${sw:-1}"
fi
finish
