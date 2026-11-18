#!/bin/bash
# P14: the hot JPL nozzle, solidification on in both condensed solvers, ICE against IGLOO on one gas:
# JPL-3solver-hot-2micron (T0 2500 K, 2 um alumina that freezes completely: --design full) and
# JPL-3solver-hot-20micron (T0 2400 K, 20 um, T-nuc 2227 K, drops that leave on the plateau: --design partial).
# Per case: G1 the inputs (p14_hot_nozzle.py g1), a scratch copy with p18's two exit probes of the gas Mach number
# (p14_hot_nozzle.py inject), run_case.sh at P14_THREADS threads (default 8), G2 the run witnesses (p14_hot_nozzle.py
# g2: hygiene, solidification active, W1 MOSE residual, W1b exit-probe range, W2 ICE's last two fields, W3 IGLOO's
# last sweeps), the comparison (JPL-3solver-hot-2micron/compare_freezing.py) and its refutation (refute_compare.py:
# every check fails on a copy doctored against it, and only that one).
# Measured 2026-11-18 (hydra 2afbfc8, ifx, 8 threads, 50 x 500/2/250 as the cases ship): each run 15-20 min; RED today
# on the gas witnesses of both cases, nothing else: 2 um W1 4.9e-8 (still falling; 1.5e-8, a floor, in blocks of 250),
# W1b 5.6e-9; 20 um W1 1.7e-7, W1b 6.4e-7 (the axis probe jumps by ~1e-6 at every other coupling step, in blocks of 250
# too). The comparison passes both (9 and 7 checks), the refutation too (12 and 9). The cold JPL-3solver-20micron with
# these output settings: W1 3.7e-9, W1b 1.1e-9, W2, W3 PASS. RED on: the cold case directory (G1 (a), (c), (d)); a
# T-nuc of 1900 K in IGLOO's phase only (G1 (a)); T-nuc 1961.6 K in ICE's only, 2 um (T2: event 18-22 cells past ICE's
# front, against N_front 4; T4, T6, C1 too); h-fus 1.635e6 in IGLOO's only, 2 um (T4: 24-27 cells, C1 -34 %, T6 67 K);
# T-nuc 2277 K in ICE's only, 20 um (T2: 23-25 cells); h-fus 1.635e6 in IGLOO's only, 20 um (C2p 0.081 against 0.009,
# C1 -35 %). The cold reference judged with P14_RUN: 'solidification active' and the comparison fail, W1 and W3 pass.
#   usage: p14_hot_nozzle_check.sh [scratch-root] [exe]
#   P14_VARIANT=2micron|20micron (default both, one after the other); P14_CASE=<dir> judges a copy of the case (a
#   mutated input) instead of test/JPL-3solver-hot-<variant>; P14_RUN=<dir> judges a stored run (input.ini, INPUT/,
#   OUTPUT/ and run.log, or <dir>.log; one variant; no solver run; without OUTPUT/ the run checks print SKIP).
set -u
export LC_ALL=C
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P14}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
NT=${P14_THREADS:-8}
M=$H/test/mi2-checks
CMP=$H/test/JPL-3solver-hot-2micron
VARIANTS=${P14_VARIANT:-2micron 20micron}
if [ -n "${P14_RUN:-}" ] && [ -z "${P14_VARIANT:-}" ]; then
  echo "P14_RUN judges one run: set P14_VARIANT=2micron or 20micron"; echo "== 0 passed, 1 failed"; exit 1
fi
mkdir -p "$T"
P=0; F=0
tally() {  # tally <file>: count its PASS / FAIL lines
  P=$(( P + $(grep -c '^PASS' "$1") )); F=$(( F + $(grep -c '^FAIL' "$1") ))
}
for v in $VARIANTS; do
  case $v in
    2micron)  d=2um;  design=full;;
    20micron) d=20um; design=partial;;
    *) echo "unknown variant $v"; F=$(( F + 1 )); continue;;
  esac
  case_dir=${P14_CASE:-$H/test/JPL-3solver-hot-$v}
  echo "== $v (${case_dir})"
  python3 -B "$M/p14_hot_nozzle.py" g1 "$case_dir" "$d" > "$T/$v.g1" 2>&1; cat "$T/$v.g1"; tally "$T/$v.g1"
  if [ -n "${P14_RUN:-}" ]; then
    run=$P14_RUN; log=$run/run.log; [ -f "$run.log" ] && log=$run.log
  else
    s=$T/src/$v; run=$T/$v; log=$T/$v.log
    rm -rf "$s" "$run" "$log"; mkdir -p "$s"
    cp -rp "$case_dir/INPUT" "$case_dir/MESH" "$case_dir/input.ini" "$s/"
    python3 -B "$M/p14_hot_nozzle.py" inject "$s/input.ini" || { echo "FAIL  probes not injected"; F=$(( F + 1 )); continue; }
    echo "  exe $EXE ($(md5sum < "$EXE" | cut -c1-8)), $NT threads"
    t0=$(date +%s)
    "$M/run_case.sh" "$s" "$run" "$NT" "$log" "$EXE" > /dev/null 2>&1
    echo "  run $(( $(date +%s) - t0 )) s"
  fi
  if [ ! -d "$run/OUTPUT" ]; then
    echo "SKIP  no OUTPUT/ in $run: G2, the comparison and its refutation are not judged"; continue
  fi
  python3 -B "$M/p14_hot_nozzle.py" g2 "$run" "$log" "$d" > "$T/$v.g2" 2>&1; cat "$T/$v.g2"; tally "$T/$v.g2"
  python3 -B "$CMP/compare_freezing.py" "$run" --design "$design" --log "$log" > "$T/$v.cmp" 2>&1; cat "$T/$v.cmp"
  tally "$T/$v.cmp"
  python3 -B "$CMP/refute_compare.py" "$run" "$T/refute-$v" --design "$design" > "$T/$v.ref" 2>&1; cat "$T/$v.ref"
  tally "$T/$v.ref"
  if [ -z "${P14_RUN:-}" ]; then   # keep the last three ICE fields of the scratch run (each ~3 MB, 100 of them)
    ls "$run"/OUTPUT/partE-field[0-9]*.tec 2>/dev/null | sed 's/.*partE-field\([0-9]*\)\.tec/\1 &/' | sort -n | head -n -3 \
      | cut -d' ' -f2 | xargs -r rm -f
  fi
done
echo "== $P passed, $F failed"
[ "$F" -eq 0 ]
