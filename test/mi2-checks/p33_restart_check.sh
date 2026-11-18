#!/bin/bash
# P33: a MOSE + ICE restart inside hydra-MI2 ([MOSE-Parameters] and [ICE-Parameters] newrun = false) goes on as the
# run in one go, every value equal (1 thread). ICE's field file does not hold what its last coupling step hands the
# gas (the source of its last step's stages) nor what its next step starts from (the relaxation times that bound its
# steps, a smoke family's birth rate): MI2 writes them, OUTPUT/<ICE prefix>coupling.tec, at the end of every coupling
# step that leaves both solvers' fields on disk and at the end of the run, and an ICE restart reads them back.
# Restarts judged against the run in one go (p33_restart.py):
#  jpl     : test/JPL-Eulerian-20micron, MOSE and ICE steady, coupling steps of 100 and 50 iterations (2 + 2 against 4):
#            ICE's source alone;
#  odd     : the jpl leg with sol-diter 300 (MOSE) and 150 (ICE), which divide neither part's iterations: the end of
#            part 1 is no output iteration, and its final fields are still written and restarted from;
#  smoke   : test/steady-ice-smoke on the jpl schedule: a burning family whose smoke ICE births as a rate (srate) and
#            takes out of what the gas receives (birth); the smoke family's steps bounded by its relaxation time;
#  removal : test/steady-ice-burnout with MOSE steady and ICE time-accurate (first order, cfl 0.4, K-burn x100) in
#            alternate cycles, MOSE-iter 100 0, ICE-iter 0 50 (4 cycles + 4 against 8): the removals ICE stores and the
#            time they were made over (goneRate), ICE's source averaged over its block;
#  kill    : the jpl case run for 6 coupling steps with outputs at the step ends only, killed (kill -9 of its PID) in
#            the MOSE block of step 3 or later: restarted for the steps left, it equals the run of 6 in one go, the
#            residual histories included (each solver drops the rows the killed run wrote past its restart iteration).
# The end of each run is its end-of-run output: no history holds an iteration twice in a row, each solver's output
# number is its iterations // sol-diter. Every run: GAME OVER, no error line, exit 0.
# Restarts whose hand-off cannot be used (missing, out of step with the fields, counters not resumed), each judged on
# its message and on how the run ends: the first MOSE block would run without ICE's source, so the run stops when a
# cycle in which ICE iterates is time-accurate (the gas steps in physical time, or ICE does) or when ICE iterates in no
# cycle (the hand-off would be the gas's only source of ICE); otherwise one warning and a completed run:
#  midstep : the kill leg with MOSE outputs in the middle of a step (sol-diter 50), killed after a MOSE output in
#            the middle of step 3 or later (250, 350, ...): the fields are not of one coupling step; steady: a warning,
#            the hand-off unused, the run completes;
#  kind    : steady-ice-burnout run steady (2 cycles), restarted with ICE time-accurate: ICE starts its counters anew;
#            ICE time-accurate: stops;
#  stale   : a new run (newrun = true) in an OUTPUT/ holding jpl part 1's files removes their hand-off; killed after
#            ICE's output of iteration 100 (the stale file's), its restart never reads it: steady, a warning, completes;
#  missing : the jpl restart without part-coupling.tec: steady, a warning naming the file, the run completes;
#  missing-ta, outofstep-ta: the removal leg's part 1 (ICE time-accurate) restarted without its hand-off, or with the
#            hand-off's MOSE iteration changed: stops, names the way out; the way out followed runs clean
#            (missing-ta-way-out);
#  missing-gas-ta: test/ice-burn-box (MOSE and ICE time-accurate) without its hand-off: stops on the gas's clock;
#  zero-missing: the jpl restart with ICE-iter 0 and no hand-off: stops (ICE iterates in no cycle);
#  layout  : the jpl restart with a variable of part-coupling.tec renamed stops ("does not match");
#  zero    : the jpl restart with ICE-iter 0: ICE writes no end-of-run output (its field and history as part 1 left
#            them);
# A restart with IGLOO is a new IGLOO start (IGLOO keeps no parcels and no sources from one run to the next): it stops
# when a cycle in which IGLOO iterates is time-accurate (MOSE's time-accurate key on, MOSE iterating in the schedule;
# or ICE iterating in that cycle with its own key on), one warning and a completed run otherwise, nothing when no cycle
# runs IGLOO:
#  igloo       : test/evap-box-ta (MOSE time-accurate + IGLOO) restarted: stops before any iteration, names the way out
#                (the field as the initial condition, newrun = true), the fields as part 1 left them; the way out
#                followed (each field copied as the message says, newrun = true) runs clean (igloo-way-out);
#  igloo-steady: test/evap-box (MOSE steady + IGLOO) restarted: one warning, the run completes;
#  igloo-off   : test/evap-box-ta with IGLOO-iter 0: no stop, no warning;
#  mixed       : test/igloo-al-box with MOSE steady and ICE time-accurate (MOSE + IGLOO, then ICE alone): IGLOO in
#                MOSE's steady cycle, one warning, the run completes (ICE's hand-off read); with IGLOO moved to ICE's
#                cycle (mixed-ice) it stops on ICE's key, and its way out (both fields) runs clean;
#  3solver     : test/JPL-3solver-20micron (steady), 1 coupling step restarted for 1 more: ICE's hand-off read, IGLOO's
#                warning, the run completes.
#   usage: p33_restart_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P33}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
RUN=$H/test/mi2-checks/run_case.sh
PY=$H/test/mi2-checks/p33_restart.py
M=100; K=50
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
echo "  exe $EXE ($(md5sum < "$EXE" | cut -c1-8))"
IGW='newrun = false with IGLOO: IGLOO keeps no parcels and no sources'
IGS='newrun = false with IGLOO, and IGLOO iterates in cycle'
READ='ICE restarts with what its last coupling step handed the gas, from OUTPUT/part-coupling.tec'
HOW="ICE's hand-off is not used: "
HOS='ICE restarts ([ICE-Parameters] newrun = false) and its hand-off cannot be used: '

sched() {  # sched <input.ini> <N-coupling> <MOSE-iter> <ICE-iter>
  sed -i -e "s/^N-coupling *=.*/N-coupling    = $2/" -e "s/^MOSE-iter *=.*/MOSE-iter     = $3/" \
         -e "s/^ICE-iter *=.*/ICE-iter      = $4/" "$1"; }
mk() {  # mk <case> <name> <N-coupling> <MOSE-iter> <ICE-iter> <MOSE sol-diter> <ICE sol-diter> [ini sed] [models sed]
  local d=$T/$2; rm -rf "$d"; mkdir -p "$d"; cp -rp "$1/INPUT" "$1/MESH" "$1/input.ini" "$d/"
  sched "$d/input.ini" "$3" "$4" "$5"
  [ -n "${8:-}" ] && sed -i "$8" "$d/input.ini"
  [ -n "${9:-}" ] && sed -i "$9" "$d/INPUT/part-models.txt"
  awk -v m="$6" -v k="$7" '/^\[/ {s=$0} /^sol-diter/ && s=="[MOSE-IO]" {$0="sol-diter     = " m}
                           /^sol-diter/ && s=="[ICE-IO]" {$0="sol-diter   = " k} {print}' "$d/input.ini" > "$d/ini" &&
    mv "$d/ini" "$d/input.ini"
  echo "$d"; }
go() {  # go <dir> <log>: run in <dir> (OUTPUT/ kept), 1 thread; exit status
  ( cd "$1" && mkdir -p OUTPUT && ls INPUT/*phase.txt > filelist.txt && ulimit -s unlimited &&
    export KMP_STACKSIZE=100M OMP_NUM_THREADS=1 && echo "# exe=$EXE threads=1 start=$(date -Is)" > "$2" &&
    if [ -x /usr/bin/time ]; then /usr/bin/time -f "# wall=%e s  maxrss=%M kB" "$EXE" >> "$2" 2>&1
    else "$EXE" >> "$2" 2>&1; fi )
  local rc=$?; echo "# exit=$rc end=$(date -Is)" >> "$2"; return $rc; }
restart() {  # restart <part-1 run dir> <dir> <N-coupling> <MOSE-iter> <ICE-iter> [file to remove]: newrun = false, run
  rm -rf "$2"; cp -rp "$1" "$2"; sched "$2/input.ini" "$3" "$4" "$5"; [ -n "${6:-}" ] && rm -f "$2/$6"
  sed -i -e '/^\[MOSE-Parameters\]/a newrun = false' -e '/^\[ICE-Parameters\]/a newrun = false' "$2/input.ini"
  go "$2" "$T/$(basename "$2").log"; }
bad_lines() { grep -c -e '\[ERROR\]' -e 'STOP' -e 'p<0 detected' -e 'NaN' -e 'ABORT:' "$1"; }
clean() {  # clean <name> <log> <exit>
  local bad; bad=$(bad_lines "$2")
  if [ "$3" = 0 ] && grep -q 'GAME OVER' "$2" && [ "$bad" = 0 ]; then ok "$1: clean run (GAME OVER, no error line, exit 0, $(grep -o 'wall=[0-9.]* s' "$2"))"
  else ko "$1: run not clean (exit=$3, GAME OVER $(grep -c 'GAME OVER' "$2"), error lines $bad: $(grep -m2 -e '\[ERROR\]' -e 'STOP' -e 'p<0 detected' -e 'NaN' "$2" | tr '\n' ' '))"; fi; }
refused() {  # refused <name> <log> <exit> <message>: the run stops on the message, before GAME OVER
  if [ "$3" != 0 ] && grep -qF -- "$4" "$2" && ! grep -q 'GAME OVER' "$2"; then ok "$1: refused (exit $3): $(grep -m1 -F -- "$4" "$2" | cut -c1-150)"
  else ko "$1: not refused as expected (exit=$3, GAME OVER $(grep -c 'GAME OVER' "$2"), message found $(grep -cF -- "$4" "$2"))"; fi; }
warned() {  # warned <name> <log> <exit> <warning> [text that must be absent]: one warning, the run completes
  local n; n=$(grep -cF -- "$4" "$2")
  if [ "$3" = 0 ] && [ "$n" = 1 ] && grep -q 'GAME OVER' "$2" && [ "$(bad_lines "$2")" = 0 ] && \
     { [ -z "${5:-}" ] || ! grep -qF -- "$5" "$2"; }; then ok "$1: one warning, the run completes: $(grep -m1 -F -- "$4" "$2" | cut -c1-140)"
  else ko "$1: expected one warning and a completed run (exit=$3, warnings $n, GAME OVER $(grep -c 'GAME OVER' "$2"), error lines $(bad_lines "$2")${5:+, \"$5\" $(grep -cF -- "$5" "$2")})"; fi; }
judge() {  # judge <mode> <args...>: the .py's lines counted; a crash counts as a failure
  local out rc; out=$(python3 -B "$PY" "$@" 2>&1); rc=$?; echo "$out"
  pass=$((pass + $(echo "$out" | grep -c '^PASS'))); fail=$((fail + $(echo "$out" | grep -c '^FAIL')))
  if [ $rc != 0 ] && [ "$(echo "$out" | grep -c '^FAIL')" = 0 ]; then ko "$1: the checks did not run (exit $rc)"; fi; }
auxv() {  # auxv <file> <key>: the value of a DATASETAUXDATA record (last one), empty if absent
  [ -f "$1" ] && tail -c 4000 "$1" | grep -a -o "DATASETAUXDATA $2 = \"[^\"]*\"" | tail -1 | sed 's/.*= "//; s/"$//'; }
inum() { awk -v x="${1:-}" 'BEGIN { if (x == "") print -1; else printf "%d\n", x + 0 }'; }
lastrow() { [ -s "$1" ] && tail -n 1 "$1" | awk 'END { print (NR ? $1 + 0 : -1) }' || echo -1; }
alive() { [ -e "/proc/$1" ] && ! grep -q '^State:.*zombie' "/proc/$1/status" 2>/dev/null; }
until_or_end() {  # until_or_end <pid> <seconds> <test command...>: poll every 0.05 s; 1 if the run ended or time ran out
  local pid=$1 t=$2 n=0; shift 2
  until "$@"; do alive "$pid" || return 1; sleep 0.05; n=$((n+1)); [ $n -ge $((t*20)) ] && return 1; done; return 0; }
aux_ge() { [ "$(inum "$(auxv "$1" "$2")")" -ge "$3" ]; }
aux_eq() { [ "$(inum "$(auxv "$1" "$2")")" -eq "$3" ]; }
aux_mid() { local v; v=$(inum "$(auxv "$1" "$2")"); [ "$v" -ge "$3" ] && [ $((v % $4)) -eq $(($4 / 2)) ]; }
row_ge() { [ "$(lastrow "$1")" -ge "$2" ]; }
way_out() {  # way_out <name> <stopped run dir>: follow the stop's way out in a copy (each 'copy <field> to <ic>' of its
             # message, newrun = false lines removed, OUTPUT/ emptied) and expect a clean new run
  local d=$T/$1 log=$T/$1.log pairs n=0 a b
  rm -rf "$d"; cp -rp "$2" "$d"; rm -rf "$d/OUTPUT"; mkdir -p "$d/OUTPUT"
  pairs=$(grep -o 'OUTPUT/[^ ]*field[^ ]* to INPUT/[^ ]*ic\.[a-z]*' "$T/$(basename "$2").log" | sed 's/ to / /')
  while read -r a b; do [ -n "$a" ] && cp -p "$2/$a" "$d/$b" && n=$((n+1)); done <<< "$pairs"
  sed -i '/^newrun = false$/d' "$d/input.ini"
  go "$d" "$log"; local rc=$?
  if [ $n -ge 1 ]; then clean "$1 ($n field(s) as initial condition, newrun = true)" "$log" $rc
  else ko "$1: the stop names no field to copy"; fi; }
killed_run() {  # killed_run <dir> <log> <test command...>: start the run, kill -9 its PID once the test holds
  ( cd "$1" && mkdir -p OUTPUT && ls INPUT/*phase.txt > filelist.txt && ulimit -s unlimited &&
    export KMP_STACKSIZE=100M OMP_NUM_THREADS=1 && exec "$EXE" > "$2" 2>&1 ) &
  local pid=$! d=$1; shift 2
  if until_or_end $pid 600 "$@"; then kill -9 $pid; wait $pid 2>/dev/null; return 0; fi
  wait $pid 2>/dev/null; return 1; }

leg() {  # leg <tag> <case> <part-1 N-coupling> <MOSE-iter> <ICE-iter> <one-go N-coupling> <MOSE-iter> <ICE-iter>
        #     [ini sed] [models sed]: part 1, the run in one go, part 1 restarted on part 1's schedule; judged
  local d t=$1 c=$2 sm=${SDM:-$M} sk=${SDK:-$K}
  d=$(mk "$c" "$t-part1" "$3" "$4" "$5" $sm $sk "${9:-}" "${10:-}"); "$RUN" "$d" "$d/run" 1 "$T/$t-part1.log" "$EXE" > /dev/null 2>&1
  clean "$t-part1" "$T/$t-part1.log" $?
  d=$(mk "$c" "$t-one-go" "$6" "$7" "$8" $sm $sk "${9:-}" "${10:-}"); "$RUN" "$d" "$d/run" 1 "$T/$t-one-go.log" "$EXE" > /dev/null 2>&1
  clean "$t-one-go" "$T/$t-one-go.log" $?
  restart "$T/$t-part1/run" "$T/$t-part2" "$3" "$4" "$5"; clean "$t-part2" "$T/$t-part2.log" $?
  if grep -q 'Initial conditions file *OUTPUT/gas-field.tec' "$T/$t-part2.log" && grep -qF "$READ" "$T/$t-part2.log"
  then ok "$t-part2: MOSE reads OUTPUT/gas-field.tec, MI2 reads ICE's hand-off from OUTPUT/part-coupling.tec"
  else ko "$t-part2: the restart does not name OUTPUT/gas-field.tec and OUTPUT/part-coupling.tec as what it reads"; fi
  judge exact "$T/$t-one-go/run/OUTPUT" "$T/$t-part2/OUTPUT" "$T/$t-part1/run/OUTPUT" $sm $sk; }

JPL=$H/test/JPL-Eulerian-20micron
FO='s/^time-scheme          = RK3$/time-scheme          = euler/; s/^space-reconstruction = MUSCL$/space-reconstruction = first-order/; /^flux-limiter         = vanleer$/d'
TA='s/^time-accurate        = false$/time-accurate        = true/; s/^cfl                  = 0.2$/cfl                  = 0.4/'
leg jpl     "$JPL" 2 $M $K 4 $M $K
SDM=300 SDK=150 leg odd "$JPL" 2 $M $K 4 $M $K
leg smoke   "$H/test/steady-ice-smoke"      2 $M $K 4 $M $K
leg removal "$H/test/steady-ice-burnout"    "1 1 1 1" "$M 0 $M 0" "0 $K 0 $K" \
            "1 1 1 1 1 1 1 1" "$M 0 $M 0 $M 0 $M 0" "0 $K 0 $K 0 $K 0 $K" "$TA; $FO" 's/ K-burn=1.39e-06 / K-burn=1.39e-04 /'

# --- kill: a run killed with kill -9 in the MOSE block of step 3 or later restarts as the run in one go ---
d=$(mk "$JPL" kill-one-go 6 $M $K $M $K); "$RUN" "$d" "$d/run" 1 "$T/kill-one-go.log" "$EXE" > /dev/null 2>&1
clean kill-one-go "$T/kill-one-go.log" $?
d=$(mk "$JPL" kill-part1 6 $M $K $M $K); mkdir -p "$d/run"; cp -rp "$d/INPUT" "$d/MESH" "$d/input.ini" "$d/run/"
o=$d/run/OUTPUT
# step 2's end is past once MOSE's history holds iteration 210: the fields and the hand-off of step 2 are on disk
if killed_run "$d/run" "$T/kill-part1.log" row_ge "$o/gas-residual-history.dat" 210; then
  im=$(inum "$(auxv "$o/gas-field.tec" MOSE_ITERATION)"); ic=$(inum "$(auxv "$o/part-field.tec" ICE_ITERATION)")
  ok "kill-part1: killed (kill -9) with MOSE's history at iteration $(lastrow "$o/gas-residual-history.dat"); fields on disk of MOSE iteration $im, ICE iteration $ic"
  cp -rp "$d/run" "$T/kill-part1-killed"
  k=$((im / M))
  if [ $im = $((k * M)) ] && [ $ic = $((k * K)) ] && [ $k -ge 2 ] && [ $k -lt 6 ]; then
    restart "$d/run" "$T/kill-part2" $((6 - k)) $M $K; clean kill-part2 "$T/kill-part2.log" $?
    grep -qF "$READ" "$T/kill-part2.log" && ok "kill-part2: MI2 reads ICE's hand-off of the killed run's step $k" \
                                         || ko "kill-part2: the hand-off of step $k was not read: $(grep -m1 -e WARNING -e ERROR "$T/kill-part2.log" | cut -c1-160)"
    judge exact "$T/kill-one-go/run/OUTPUT" "$T/kill-part2/OUTPUT" "$T/kill-part1-killed/OUTPUT" $M $K $im $ic
  else ko "kill-part1: the fields on disk are not of one step end ($im, $ic): no restart point"; fi
else ko "kill-part1: the run ended before MOSE's history reached iteration 210 (nothing killed)"; fi

# --- midstep: killed after MOSE's output in the middle of step 3: not one coupling step ---
d=$(mk "$JPL" midstep-part1 6 $M $K $K $K); mkdir -p "$d/run"; cp -rp "$d/INPUT" "$d/MESH" "$d/input.ini" "$d/run/"
o=$d/run/OUTPUT
if killed_run "$d/run" "$T/midstep-part1.log" aux_mid "$o/gas-field.tec" MOSE_ITERATION 250 $M; then
  ok "midstep-part1: killed after MOSE's output of iteration $(auxv "$o/gas-field.tec" MOSE_ITERATION) (ICE's field of iteration $(inum "$(auxv "$o/part-field.tec" ICE_ITERATION)"))"
  restart "$d/run" "$T/midstep" 3 $M $K
  warned midstep "$T/midstep.log" $? "${HOW}OUTPUT/part-coupling.tec is of ICE iteration" "$READ"
  grep -qF 'the run that wrote them stopped between two of its outputs. The first MOSE block runs without ICE' "$T/midstep.log" \
    && ok "midstep: the warning says why (stopped between two outputs) and what follows (no ICE source in the first MOSE block)" \
    || ko "midstep: the warning does not say why and what follows"
else ko "midstep-part1: the run ended before MOSE's output of iteration 250 (nothing killed)"; fi

# --- kind: a steady run restarted with ICE time-accurate ---
d=$(mk "$H/test/steady-ice-burnout" kind-part1 "1 1" "$M 0" "0 $K" $M $K "$FO" 's/ K-burn=1.39e-06 / K-burn=1.39e-04 /')
"$RUN" "$d" "$d/run" 1 "$T/kind-part1.log" "$EXE" > /dev/null 2>&1; clean kind-part1 "$T/kind-part1.log" $?
rm -rf "$T/kind"; cp -rp "$d/run" "$T/kind-src"; sed -i "$TA" "$T/kind-src/input.ini"
restart "$T/kind-src" "$T/kind" "1 1" "$M 0" "0 $K"; rc=$?; rm -rf "$T/kind-src"
refused kind "$T/kind.log" $rc "${HOS}ICE restarts without the counters of the run it continues (a run of the other kind, steady or time-accurate, or a field without them). ICE iterates in cycle 2, a time-accurate cycle: ICE iterates in it in physical time"

# --- stale: a new run removes the hand-off an earlier run left; killed after ICE's output of iteration 100 ---
d=$(mk "$JPL" stale-part1 3 200 200 200 $M); mkdir -p "$d/run"; cp -rp "$d/INPUT" "$d/MESH" "$d/input.ini" "$d/run/"
mkdir -p "$d/run/OUTPUT"; cp -p "$T/jpl-part1/run/OUTPUT/part-coupling.tec" "$d/run/OUTPUT/"; o=$d/run/OUTPUT
echo "  stale: jpl part 1's hand-off (ICE $(auxv "$o/part-coupling.tec" CYCLE), MOSE $(auxv "$o/part-coupling.tec" MI2_MOSE_ITERATION)) planted"
if killed_run "$d/run" "$T/stale-part1.log" aux_eq "$o/part-field.tec" ICE_ITERATION 100; then
  if grep -qF 'ICE starts anew: OUTPUT/part-coupling.tec of an earlier run removed' "$T/stale-part1.log" && [ ! -e "$o/part-coupling.tec" ]
  then ok "stale-part1: the new run removed the earlier hand-off; killed with its fields of MOSE $(auxv "$o/gas-field.tec" MOSE_ITERATION), ICE $(auxv "$o/part-field.tec" ICE_ITERATION)"
  else ko "stale-part1: the earlier hand-off was not removed ($(ls "$o" | grep -c coupling) file(s) left)"; fi
  restart "$d/run" "$T/stale" 1 200 200
  warned stale "$T/stale.log" $? "${HOW}OUTPUT/part-coupling.tec is missing" "$READ"
else ko "stale-part1: the run ended before ICE's output of iteration 100 (nothing killed)"; fi

# --- missing, layout: the jpl restart without its hand-off, with a renamed variable ---
restart "$T/jpl-part1/run" "$T/missing" 2 $M $K OUTPUT/part-coupling.tec
warned missing "$T/missing.log" $? "${HOW}OUTPUT/part-coupling.tec is missing" "$READ"
restart "$T/jpl-part1/run" "$T/zero-missing" 2 $M 0 OUTPUT/part-coupling.tec
refused zero-missing "$T/zero-missing.log" $? "${HOS}OUTPUT/part-coupling.tec is missing (MI2 writes it with the fields at the end of each coupling step that leaves both solvers' fields on disk and at the end of a run with ICE; an older binary did not). ICE iterates in no cycle of this run"

# --- missing-ta, outofstep-ta, missing-gas-ta: a hand-off that cannot be used in a time-accurate schedule stops ---
restart "$T/removal-part1/run" "$T/missing-ta" "1 1 1 1" "$M 0 $M 0" "0 $K 0 $K" OUTPUT/part-coupling.tec
refused missing-ta "$T/missing-ta.log" $? "${HOS}OUTPUT/part-coupling.tec is missing"
grep -qF 'ICE iterates in cycle 2, a time-accurate cycle: ICE iterates in it in physical time ([ICE-Numerics] time-accurate = true). Without the hand-off the first MOSE block would run without ICE' "$T/missing-ta.log" && \
  grep -qF 'start a new run from them: copy OUTPUT/gas-field.tec to INPUT/gas-ic.tec and OUTPUT/part-field.tec to INPUT/part-ic.tec' "$T/missing-ta.log" \
  && ok "missing-ta: the stop names the time-accurate cycle (ICE's key) and the way out (both fields)" \
  || ko "missing-ta: the stop does not name the cycle, its clock or the way out"
way_out missing-ta-way-out "$T/missing-ta"
rm -rf "$T/oos-src"; cp -rp "$T/removal-part1/run" "$T/oos-src"
sed -i 's/\(DATASETAUXDATA MI2_MOSE_ITERATION = "\)[^"]*"/\1-7"/' "$T/oos-src/OUTPUT/part-coupling.tec"
restart "$T/oos-src" "$T/outofstep-ta" "1 1 1 1" "$M 0 $M 0" "0 $K 0 $K"; rc=$?; rm -rf "$T/oos-src"
refused outofstep-ta "$T/outofstep-ta.log" $rc "and MOSE iteration -7, the fields restored of ICE iteration"
d=$(mk "$H/test/ice-burn-box" gas-ta-part1 "1 2" "0 1" "1 0" 1 1); "$RUN" "$d" "$d/run" 1 "$T/gas-ta-part1.log" "$EXE" > /dev/null 2>&1
clean gas-ta-part1 "$T/gas-ta-part1.log" $?
restart "$T/gas-ta-part1/run" "$T/missing-gas-ta" "1 2" "0 1" "1 0" OUTPUT/part-coupling.tec
refused missing-gas-ta "$T/missing-gas-ta.log" $? "OUTPUT/part-coupling.tec is missing (MI2 writes it with the fields at the end of each coupling step that leaves both solvers' fields on disk and at the end of a run with ICE; an older binary did not). ICE iterates in cycle 1, a time-accurate cycle: the gas ICE hands its sources to steps in physical time ([MOSE-Numerics] time-accurate = true)"
rm -rf "$T/layout-src"; cp -rp "$T/jpl-part1/run" "$T/layout-src"; sed -i '1s/"S1_p1"/"X1_p1"/; 1s/ S1_p1 / X1_p1 /' "$T/layout-src/OUTPUT/part-coupling.tec"
restart "$T/layout-src" "$T/layout" 2 $M $K; rc=$?; rm -rf "$T/layout-src"
refused layout "$T/layout.log" $rc 'does not match the ICE set-up of this run'

# --- zero: ICE restarted with no iteration writes no end-of-run output ---
restart "$T/jpl-part1/run" "$T/zero" 2 $M 0; clean zero "$T/zero.log" $?
p=$T/jpl-part1/run/OUTPUT; z=$T/zero/OUTPUT
if cmp -s "$p/part-field.tec" "$z/part-field.tec" && cmp -s "$p/part-residual-history.dat" "$z/part-residual-history.dat"; then
  ok "zero: ICE made no iteration: part-field.tec and part-residual-history.dat as part 1 left them ($(grep -c . "$z/part-residual-history.dat") rows, output number $(auxv "$z/part-field.tec" ICE_OUTPUT_NUMBER))"
else ko "zero: ICE made no iteration and rewrote its files (rows $(grep -c . "$p/part-residual-history.dat") -> $(grep -c . "$z/part-residual-history.dat"), output number $(auxv "$p/part-field.tec" ICE_OUTPUT_NUMBER) -> $(auxv "$z/part-field.tec" ICE_OUTPUT_NUMBER))"; fi

# --- igloo: a restart with IGLOO is a new IGLOO start: it stops in a time-accurate cycle, warns in a steady one ---
d=$(mk "$H/test/evap-box-ta" igloo-part1 1 50 0 50 1); "$RUN" "$d" "$d/run" 1 "$T/igloo-part1.log" "$EXE" > /dev/null 2>&1
clean igloo-part1 "$T/igloo-part1.log" $?
restart "$T/igloo-part1/run" "$T/igloo" 1 50 0
refused igloo "$T/igloo.log" $? "$IGS 1, a time-accurate cycle: the gas IGLOO hands its sources to steps in physical time"
p=$T/igloo-part1/run/OUTPUT; z=$T/igloo/OUTPUT
if grep -qF 'start a new run from them: copy OUTPUT/gas-field.tec to INPUT/gas-ic.tec' "$T/igloo.log" && \
   grep -qF 'set newrun = true in [MOSE-Parameters]' "$T/igloo.log" && cmp -s "$p/gas-field.tec" "$z/gas-field.tec" && \
   cmp -s "$p/gas-residual-history.dat" "$z/gas-residual-history.dat"
then ok "igloo: the stop names the way out (OUTPUT/gas-field.tec to INPUT/gas-ic.tec, newrun = true); gas-field.tec and gas-residual-history.dat as part 1 left them"
else ko "igloo: the way out is not named, or the stopped run changed the fields (field equal: $(cmp -s "$p/gas-field.tec" "$z/gas-field.tec" && echo yes || echo no))"; fi
way_out igloo-way-out "$T/igloo"
d=$(mk "$H/test/evap-box" igloo-steady-part1 1 50 0 50 1); "$RUN" "$d" "$d/run" 1 "$T/igloo-steady-part1.log" "$EXE" > /dev/null 2>&1
clean igloo-steady-part1 "$T/igloo-steady-part1.log" $?
restart "$T/igloo-steady-part1/run" "$T/igloo-steady" 1 50 0
warned igloo-steady "$T/igloo-steady.log" $? "$IGW" "$IGS"
d=$(mk "$H/test/evap-box-ta" igloo-off-part1 1 50 0 50 1 's/^IGLOO-iter *=.*/IGLOO-iter    = 0/')
"$RUN" "$d" "$d/run" 1 "$T/igloo-off-part1.log" "$EXE" > /dev/null 2>&1; clean igloo-off-part1 "$T/igloo-off-part1.log" $?
restart "$T/igloo-off-part1/run" "$T/igloo-off" 1 50 0; rc=$?; clean igloo-off "$T/igloo-off.log" $rc
grep -qF 'newrun = false with IGLOO' "$T/igloo-off.log" && ko "igloo-off: IGLOO-iter = 0 everywhere and an IGLOO stop or warning printed" \
                                                     || ok "igloo-off: IGLOO-iter = 0 everywhere: no IGLOO stop, no IGLOO warning"

# --- mixed: MOSE steady with IGLOO, then ICE time-accurate alone; IGLOO moved to ICE's cycle stops on ICE's key ---
MS='/^\[MOSE-Numerics\]/,/^\[/ s/^time-accurate *= *true/time-accurate         = false/'
d=$(mk "$H/test/igloo-al-box" mixed-part1 "1 1" "50 0" "0 1" 50 1 "$MS"); "$RUN" "$d" "$d/run" 1 "$T/mixed-part1.log" "$EXE" > /dev/null 2>&1
clean mixed-part1 "$T/mixed-part1.log" $?
restart "$T/mixed-part1/run" "$T/mixed" "1 1" "50 0" "0 1"; rc=$?
warned mixed "$T/mixed.log" $rc "$IGW" "$IGS"
grep -qF 'ICE restarts with what its last coupling step handed the gas, from OUTPUT/partE-coupling.tec' "$T/mixed.log" \
  && ok "mixed: MI2 reads ICE's hand-off (OUTPUT/partE-coupling.tec)" \
  || ko "mixed: ICE's hand-off not read: $(grep -m1 -e WARNING -e ERROR "$T/mixed.log" | cut -c1-160)"
rm -rf "$T/mixed-ice-src"; cp -rp "$T/mixed-part1/run" "$T/mixed-ice-src"
sed -i 's/^IGLOO-iter *=.*/IGLOO-iter    = 0 1/' "$T/mixed-ice-src/input.ini"
restart "$T/mixed-ice-src" "$T/mixed-ice" "1 1" "50 0" "0 1"; rc=$?; rm -rf "$T/mixed-ice-src"
refused mixed-ice "$T/mixed-ice.log" $rc "$IGS 2, a time-accurate cycle: ICE iterates in it in physical time ([ICE-Numerics] time-accurate = true)"
way_out mixed-ice-way-out "$T/mixed-ice"

# --- 3solver: MOSE + ICE + IGLOO, 1 coupling step restarted for 1 more ---
d=$(mk "$H/test/JPL-3solver-20micron" 3solver-part1 1 $M $K $M $K); "$RUN" "$d" "$d/run" 1 "$T/3solver-part1.log" "$EXE" > /dev/null 2>&1
clean 3solver-part1 "$T/3solver-part1.log" $?
restart "$T/3solver-part1/run" "$T/3solver" 1 $M $K
warned 3solver "$T/3solver.log" $? "$IGW" "$IGS"
grep -qF 'ICE restarts with what its last coupling step handed the gas, from OUTPUT/partE-coupling.tec' "$T/3solver.log" \
  && ok "3solver: MI2 reads ICE's hand-off (OUTPUT/partE-coupling.tec) with IGLOO in the case" \
  || ko "3solver: ICE's hand-off not read: $(grep -m1 -e WARNING -e ERROR "$T/3solver.log" | cut -c1-160)"
echo "== $pass passed, $fail failed"; [ $fail = 0 ]
