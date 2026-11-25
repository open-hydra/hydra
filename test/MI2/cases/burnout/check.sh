#!/bin/bash
# check.sh <run-dir> <exe>
# The variants of ice-burnout (a material without a residue burns out; what is left of a drop at d-burnout is handed to
# the gas whole), each judged by ice-burnout's check.py in the mode named:
#   long   : a step twenty times longer (dt-max 2e-4, the relaxation-time ceiling off), too long for the end of the
#            burn: the removal fires at the first step that leaves less than it burnt, the books close;
#   frozen : d-burnout above the drops and one ICE step: the removal fires in it, its source is the burn alone and the
#            rest of the cloud is stored; frozen through 5 MOSE iterations, the gas gains the frozen burn and the
#            store once (the inter-phase books of a removal);
#   frozen2: the same with ICE-iter 2: the second step finds the cloud empty; the gas takes the block's rate (the two
#            steps' sources weighted by their dt) and the first one's store;
#   cold   : T-ign above any temperature the drops reach and d-burnout above them: no drop burns, none is removed; the
#            same on the IGLOO twin (twins/igloo on Beckstead removed at d-burnout): IGLOO hands the gas no mass;
#   refuse : d-stop beside d-burnout on the models row stops ICE and IGLOO before the first step.
source "$(dirname "${BASH_SOURCE[0]}")/../../checklib.sh"
source "$cases/twins/beckstead.sh"
I=$cases/ice-burnout; py=$I/check.py

# edit <dir> <file> <sed>...: each edit must change the file
edit() {
  local d=$1 f=$2; shift 2
  for e in "$@"; do
    cp "$d/$f" "$d/$f.0"; sed -i "$e" "$d/$f"
    cmp -s "$d/$f" "$d/$f.0" && ko "$(basename "$d"): the edit '$e' of $f changed nothing"; rm -f "$d/$f.0"
  done
}
# mode <name> <dir> [mode]: run, then ice-burnout's check in the mode; its PASS/FAIL lines counted
mode() {
  runin "$2"; local rc=$? out
  expect "$1: the run ends (exit $rc)" ended "$2/run.log"
  out=$(python3 -B "$py" "$2" "$2/run.log" "$cases/ice-exchange" ${3:-} 2>&1)
  echo "$out" | grep -E '^(PASS|FAIL)' | sed "s/^\(PASS\|FAIL\) */\1  $1: /"
  npass=$((npass + $(echo "$out" | grep -c '^PASS'))); nfail=$((nfail + $(echo "$out" | grep -c '^FAIL')))
  echo "$out" | grep -q '^== ' || ko "$1: the check did not finish: $(echo "$out" | tail -1)"
  last=$out
}

d=$(variant long "$I")
edit "$d" input.ini 's/^dt-max  *= 1.0e-5$/dt-max               = 2.0e-4\ntau-factor           = 0/' \
                    's/^N-coupling    = 600 1$/N-coupling    = 60 1/'
mode long "$d"
expect "long: the removal fires by the rule of the long step" grep -q 'R4 the removal fires by the rule of the long step' <<< "$last"

d=$(variant frozen "$I")
edit "$d" input.ini 's/^dt-max  *= 1.0e-5$/dt-max               = 1.0e-4/' 's/^N-coupling    = 600 1$/N-coupling    = 1 5/' \
                    's/^d-burnout  = 3.0e-6$/d-burnout  = 4.0e-5/' 's/^sol-diter     = 200$/sol-diter     = 5/'
edit "$d" INPUT/part-models.txt 's/ d-burnout=3e-06$/ d-burnout=4e-05/'
mode frozen "$d" frozen

d=$(variant frozen2 "$run/frozen"); edit "$d" input.ini 's/^ICE-iter      = 1 0$/ICE-iter      = 2 0/'
mode frozen2 "$d" frozen2

d=$(variant cold "$I")
edit "$d" input.ini 's/^d-burnout  = 3.0e-6$/d-burnout  = 4.0e-5/' 's/^T-ign      = 290$/T-ign      = 500/'
edit "$d" INPUT/part-models.txt 's/ T-ign=290 / T-ign=500 /' 's/ d-burnout=3e-06$/ d-burnout=4e-05/'
mode cold "$d" cold

d=$(variant cold-igloo "$cases/twins/igloo"); beckstead "$d" 'd-burnout  = 4.0e-5' 'd-burnout=4e-05'
edit "$d" INPUT/part-models.txt 's/ T-ign=290 / T-ign=500 /'
runin "$d"; rc=$?
expect "cold-igloo: the run ends (exit $rc)" ended "$d/run.log"
# round-off of IGLOO's in - out sums: 2.5e-21 kg/s measured; a removal of parcels that do not burn hands 4e-6 kg/s
m=$(grep -a '\[MI2-report\] igloo-to-gas ' "$d/run.log" | awk '{n++; v = $NF; sub(/^-/, "", v); if (v + 0 > big) big = v + 0}
                                                             END {printf "%d %.3e %d", n, big, (n > 0 && big <= 1e-15)}')
expect "cold-igloo: IGLOO hands the gas no mass (${m% *}: blocks, largest kg/s; tol 1e-15)" [ "${m##* }" = 1 ]

for c in ice igloo; do
  if [ $c = ice ]; then d=$(variant refuse-ice "$I")
  else d=$(variant refuse-igloo "$cases/twins/igloo"); beckstead "$d" 'd-burnout  = 3.0e-6' 'd-burnout=3e-06'; fi
  edit "$d" INPUT/part-models.txt 's/ d-burnout=3e-06$/ d-stop=1e-06 d-burnout=3e-06/'
  runin "$d"; rc=$?
  expect "refuse-$c: d-stop beside d-burnout stops the run (exit $rc)" \
    eval "[ $rc != 0 ] && grep -q 'd-stop beside d-burnout' '$d/run.log' && ! grep -q 'GAME OVER' '$d/run.log'"
done
finish
