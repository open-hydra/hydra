#!/bin/bash
# P14: the hot JPL nozzle, solidification on in both condensed solvers, ICE against IGLOO on one gas. The gate's leg:
# JPL-3solver-hot-2micron (T0 2500 K, 2 um alumina that freezes completely, --design full). Per case: G1 the inputs
# (p14_hot_nozzle.py g1), a scratch copy with p18's two exit probes of the gas Mach number (p14_hot_nozzle.py inject),
# run_case.sh at P14_THREADS threads (default 1: the witnesses' floors are one-thread numbers), G2 the run witnesses
# (p14_hot_nozzle.py g2: hygiene, solidification active, W1 MOSE residual, W1b exit-probe range, W2 ICE's last two
# fields, W3 IGLOO's last sweeps), the comparison (JPL-3solver-hot-2micron/compare_freezing.py) and its refutation
# (refute_compare.py: every check fails on a copy doctored against it, and only that one).
# W1 and W1b are printed for the 2 um case with the cold case's values and their bounds, not judged: the user's decision
# of 2026-11-19 ("print the steadiness figures for 2 um without judging them, and accept the ICE-against-IGLOO
# comparison on a gas that still moves slightly"). At one thread the 2 um gas still moves (W1 4.9e-8 of the first
# residual, still falling, against 1e-8; W1b 2.0e-8 against 1e-7) while every comparison number is the same to its
# printed digits at one and at eight threads. W2, W3, the hygiene and the comparison are judged.
# JPL-3solver-hot-20micron (T0 2400 K, 20 um, T-nuc 2227 K, --design partial) is no leg: P14_VARIANT=20micron judges
# it on demand with W1 and W1b as the plan writes them, and it is RED today (one thread: W1 7.8e-8, W1b 7.1e-7; the
# gas exit Mach on the axis switches between two states ~6.5e-7 apart at irregular coupling steps).
# Every python step's verdict is counted; a step that exits non-zero without a FAIL line, or with a traceback, or with
# no summary line, is one more FAIL; a witness that cannot be evaluated is a FAIL line of its own (p14_hot_nozzle.py).
# Measured 2026-11-19 (hydra 74a90a4, ifx, one thread, 50 x 500/2/250 as the case ships): the run 49-51 min on a
# pcprop; G1 6, G2 5 passed (W1 4.9e-8, W1b 2.0e-8 printed), compare 9, refute 13 passed. RED: the cold case directory
# (G1 (a), (c), (d)); T-nuc 1900 K in IGLOO's phase only (G1 (a)); T-nuc 1961.6 K in ICE's phase only (T2: events 18-22
# cells past ICE's front against N_front 4); h-fus 1.635e6 in IGLOO's phase only (T4 24-27 cells against 5, C1 -34 %);
# the cold reference judged with P14_RUN ('solidification active', T0, T2, T3, T4, T5 fail); a run without
# gas-residual-history.dat and a run read while it was written (each fails, rc 1; before: passed with rc 0).
#   usage: p14_hot_nozzle_check.sh [scratch-root] [exe]
#   P14_VARIANT=2micron|20micron (default 2micron); P14_CASE=<dir> judges a copy of the case (a mutated input) instead
#   of test/JPL-3solver-hot-<variant>; P14_RUN=<dir> judges a stored run (input.ini, INPUT/, OUTPUT/ and run.log, or
#   <dir>.log; one variant; no solver run; a directory without OUTPUT/ fails).
set -u
export LC_ALL=C
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P14}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
NT=${P14_THREADS:-1}
M=$H/test/mi2-checks
CMP=$H/test/JPL-3solver-hot-2micron
VARIANTS=${P14_VARIANT:-2micron}
if [ -n "${P14_RUN:-}" ] && [ "$(echo "$VARIANTS" | wc -w)" -ne 1 ]; then
  echo "FAIL  P14_RUN judges one run: set P14_VARIANT=2micron or 20micron"; echo "== 0 passed, 1 failed"; exit 1
fi
mkdir -p "$T"
P=0; F=0
step() {  # step <name> <file> <python args...>: run one python step, show it, count its verdict and its exit
  local name=$1 f=$2; shift 2
  python3 -B "$@" > "$f" 2>&1; local rc=$?
  cat "$f"
  P=$(( P + $(grep -c '^PASS' "$f") )); F=$(( F + $(grep -c '^FAIL' "$f") ))
  if [ "$rc" -ne 0 ] && { ! grep -q '^FAIL' "$f" || grep -q '^Traceback' "$f"; }; then
    echo "FAIL  $name exited $rc without a verdict for everything it checks"; F=$(( F + 1 ))
  elif [ "$rc" -eq 0 ] && ! grep -q '^== ' "$f"; then
    echo "FAIL  $name printed no summary line"; F=$(( F + 1 ))
  fi
}
for v in $VARIANTS; do
  case $v in
    2micron)  d=2um;  design=full;;
    20micron) d=20um; design=partial;;
    *) echo "FAIL  unknown variant $v"; F=$(( F + 1 )); continue;;
  esac
  case_dir=${P14_CASE:-$H/test/JPL-3solver-hot-$v}
  echo "== $v (${case_dir})"
  step G1 "$T/$v.g1" "$M/p14_hot_nozzle.py" g1 "$case_dir" "$d"
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
    echo "FAIL  no OUTPUT/ in $run: G2, the comparison and its refutation cannot be judged"; F=$(( F + 1 )); continue
  fi
  step G2 "$T/$v.g2" "$M/p14_hot_nozzle.py" g2 "$run" "$log" "$d"
  step compare_freezing "$T/$v.cmp" "$CMP/compare_freezing.py" "$run" --design "$design" --log "$log"
  step refute_compare "$T/$v.ref" "$CMP/refute_compare.py" "$run" "$T/refute-$v" --design "$design"
  if [ -z "${P14_RUN:-}" ]; then   # keep the last three ICE fields of the scratch run (each ~3 MB, 100 of them)
    ls "$run"/OUTPUT/partE-field[0-9]*.tec 2>/dev/null | sed 's/.*partE-field\([0-9]*\)\.tec/\1 &/' | sort -n | head -n -3 \
      | cut -d' ' -f2 | xargs -r rm -f
  fi
done
echo "== $P passed, $F failed"
[ "$F" -eq 0 ]
