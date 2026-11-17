#!/bin/bash
# Steady gate: in a run where MOSE or ICE is steady, every transfer between families and solvers is handed as a rate,
# so the converged field does not depend on the local time steps. Three cases on test/JPL-Eulerian-20micron's nozzle,
# coarsened to 38 x 8 cells, MOSE and ICE steady on JPL's schedule (MOSE-iter : ICE-iter = 2 : 1, 40 cycles) and ICE on
# JPL's scheme (MUSCL, vanleer, RK3; cfl 0.4, steady-ice-burnout 0.2: its README). Each leg runs a deck twice, as given and with one local time step changed (ICE
# tau-factor 1 vs 0.5 where the smoke family's step is relaxation-limited, ICE cfl 0.2 vs 0.1 where the drops' is not);
# p31_steady_cmp.py takes the largest difference over the cells of each variable, relative to its peak (a velocity
# component relative to its phase's peak speed). Each case also runs at first order (forward Euler, first order):
#  births   test/steady-ice-smoke    ICE's smoke births (srate)
#  removal  test/steady-ice-burnout  the steady burn-out: f of what leaves a cell burns at its exit, f rising from 0 at
#           d-burnout to 1 at the band's lower edge (ICE's `burnout_band` = [0.90, 1], x d-burnout) (cfl)
#  igloo    test/steady-igloo-smoke  IGLOO's smoke born in ICE (ext)
#  mixed    test/steady-ice-burnout with MOSE steady and ICE time-accurate (first order, cfl 0.4) in alternate cycles
#           (the removal at once, handed as a rate): MOSE cfl 0.5 vs 0.25, the gas field only
#  burn-fo  test/steady-ice-burnout at first order, K-burn x100 (the drops burn out within the first cells), ICE cfl
#           0.2 vs 0.1
#  d0-fo    the same without d-burnout (the removal at m_burnout, d -> 0): one run, its books
#  agree    removal-fo's steady run against mixed's time-accurate one (after both): the metal at the outlet within
#           1e-6 of the hand-off, judged (the band's drops leave in part: 3e-8 at [0.90, 1]);
#           the mass and the energy handed to the gas against the time-accurate cycle's mean and range (its blocks
#           over the second half), stated (a STATE line), not judged: the time-accurate run keeps the hard switch, the
#           steady form the smooth sink, two models on this deck
#  agree-d0k10, agree-w200: steady against time-accurate where the two models coincide, each judged inside its own
#           time-accurate run's range (the second half of 120 pairs), mass and energy, and the outlet as agree's: the
#           d0-fo deck at K-burn x10 (no d-burnout: the removal at m_burnout), and the shipped deck at first order with
#           ICE's wall taking no drops (face-4 records type 200 in INPUT/part-bc.txt)
#  agree-d0 d0-fo's steady run against a 120-pair time-accurate twin (after d0-fo): the outlet judged, the hand-off
#           stated (a STATE line): there the switch fires by over-burn in the first cell, the steady sink takes the band
#           at m_burnout's drop, and the energy differs by about the range's half-width
#  band     ICE's `burnout_band` the same in the log (removal-fo's), Mod_BC_Fluxes.f90 and the texts that name it
# Tolerances: ten times each leg's floor rounded up to a power of ten; floor: the same deck run 40 against 50 cycles
# (30 against 31 cycle pairs, mixed), ifx: births 1.4e-14, births-fo 5.5e-14, removal 3.0e-14,
# removal-fo 1.2e-14, igloo 1.2e-9, igloo-fo 1.0e-9, mixed 1.2e-2. burn-fo: the pair itself (1.1e-14), its
# 40-against-50 difference being the coupled field still moving. mixed keeps 1e-1, 8 times its floor and not ten:
# the floor rose from 4.6e-3 when the switch began to keep an emptied cell's velocity, and the rule's 1 judges nothing.
# For mixed 1e-1 is a sanity bound: the leg is the extra one, step independence is judged on the strict legs.
# Balance (p31_balance.py) of the first-order legs' converged fields, which a step-independent wrong model does not
# pass: births-fo the smoke leaving against nu_s/(1 - nu_s) ice-to-gas-mass (5.6e-11, the print precision of nu and
# M; tol 1e-9), igloo-fo the smoke leaving against igloo-smoke-rate (3.6e-15; 1e-13), removal-fo, burn-fo and d0-fo
# the burn-out's books, ice-to-gas-mass = SUM burn vol + SUM f out of the flagged cells (3.4e-15, 1e-13; 4.0e-7 and
# 4.1e-7, the field one block after the report, 1e-5). agree: the outlet judged; the steady hand-off against the
# time-accurate cycle stated, not judged (the time-accurate run keeps the hard switch).
# Non-vacuity: each knob changes one line of its deck and the residual history of the run; the smoke families hold
# more than 1e-4 kg/m3 at their peak;
# drops below d-burnout fill cells of the removal decks, and the coupler hands the gas a burn.
#   usage: p31_steady_check.sh [scratch-root] [exe] [legs]   (legs: a space-separated subset, default every leg)
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P31}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
LEGS=${3:-births births-fo removal removal-fo igloo igloo-fo mixed burn-fo d0-fo agree agree-d0 agree-d0k10 agree-w200 band}
RUN=$H/test/mi2-checks/run_case.sh
CMP=$H/test/mi2-checks/p31_steady_cmp.py
FLD=$H/test/mi2-checks/p31_field.py
BAL=$H/test/mi2-checks/p31_balance.py
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
echo "  exe: $EXE ($(md5sum < "$EXE" | cut -c1-8))"
TAU='s/^cfl                  = 0.4$/cfl                  = 0.4\ntau-factor           = 0.5/'
MOSECFL='s/^cfl                   = 0.5$/cfl                   = 0.25/'
ICECFL='s/^cfl                  = 0.2$/cfl                  = 0.1/'
FO='s/^time-scheme          = RK3$/time-scheme          = euler/; s/^space-reconstruction = MUSCL$/space-reconstruction = first-order/; /^flux-limiter         = vanleer$/d'
ONES=$(printf '1 %.0s' $(seq 60)); MOSE=$(printf '500 0 %.0s' $(seq 30)); ICE=$(printf '0 250 %.0s' $(seq 30))
MIX="s/^N-coupling    = 40$/N-coupling    = ${ONES% }/; s/^MOSE-iter     = 500$/MOSE-iter     = ${MOSE% }/; s/^ICE-iter      = 250$/ICE-iter      = ${ICE% }/; s/^time-accurate        = false$/time-accurate        = true/; s/^cfl                  = 0.2$/cfl                  = 0.4/"
K100='s/ K-burn=1.39e-06 / K-burn=1.39e-04 /'
NODB='s/ d-burnout=1.7e-05$//'
K10='s/ K-burn=1.39e-06 / K-burn=1.39e-05 /'
# the agreement legs' time-accurate side: as MIX, 120 pairs (at 30 the K-burn x10 deck's run has not settled: its
# second half lies 1 % below its mean at 120)
ONES120=$(printf '1 %.0s' $(seq 240)); MOSE120=$(printf '500 0 %.0s' $(seq 120)); ICE120=$(printf '0 250 %.0s' $(seq 120))
MIX120="s/^N-coupling    = 40$/N-coupling    = ${ONES120% }/; s/^MOSE-iter     = 500$/MOSE-iter     = ${MOSE120% }/; s/^ICE-iter      = 250$/ICE-iter      = ${ICE120% }/; s/^time-accurate        = false$/time-accurate        = true/; s/^cfl                  = 0.2$/cfl                  = 0.4/"

# deck <name> <case> [sed script] [models sed script]: the case's deck, or a copy changed by the scripts (input.ini,
# INPUT/part-models.txt), in $T/<name>-deck
deck() {
  local d=$T/$1-deck c=$H/test/$2
  rm -rf "$d"; mkdir -p "$d"; cp -rp "$c/INPUT" "$c/MESH" "$c/input.ini" "$d/"
  [ -n "${3:-}" ] && sed -i "$3" "$d/input.ini"
  [ -n "${4:-}" ] && sed -i "$4" "$d/INPUT/part-models.txt"
  echo "$d"
}
# nodeposit <deck>: ICE's wall takes no drops: the deck's face-4 records (the nozzle wall, type 300) set to type 200,
# which ICE's MK flux crosses with nothing; the gas's records untouched
nodeposit() {
  python3 -B - "$1/INPUT/part-bc.txt" <<'EOF'
import sys
p = sys.argv[1]
L = open(p).read().split('\n')
n = 0
for k, l in enumerate(L):
    w = l.split()
    if len(w) == 6 and w[4] == '4' and w[5] == '300':
        L[k] = l[:-3] + '200'
        n += 1
open(p, 'w').write('\n'.join(L))
print(n)
EOF
}
# twin <leg> <steady deck> <time-accurate deck>: the two runs of an agreement leg, at once, one thread each
twin() {
  rm -rf "$T/$1-st" "$T/$1-ta"
  OMP_NUM_THREADS=1 "$RUN" "$2" "$T/$1-st" 1 "$T/$1-st.log" "$EXE" > /dev/null 2>&1 &
  OMP_NUM_THREADS=1 "$RUN" "$3" "$T/$1-ta" 1 "$T/$1-ta.log" "$EXE" > /dev/null 2>&1 &
  wait
  for r in st ta; do
    grep -q "GAME OVER" "$T/$1-$r.log" && ok "$1: the $r run ends with GAME OVER" \
      || ko "$1: the $r run: $(grep -a -m1 'ERROR\|Error in block\|NaN' "$T/$1-$r.log" | cut -c1-200)"
  done
}
# agreein <leg> <steady run> <time-accurate log>: judged inside the time-accurate range, and the outlet
agreein() {
  local o; o=$(python3 -B "$BAL" agree-in "$2" "$2.log" "$3" 2>&1); local rc=$?
  echo "$o" | sed 's/^/      /'
  [ $rc = 0 ] && ok "$1: $(echo "$o" | tail -1 | cut -c9-)" || ko "$1: $(echo "$o" | tail -1)"
}
# pair <leg> <deck> <sed script>: the deck as given and its variant, at once, one thread each
pair() {
  local leg=$1 c=$2 v=$T/$1-variant-case
  rm -rf "$T/$leg-base" "$T/$leg-variant" "$v"; mkdir -p "$v"
  cp -rp "$c/INPUT" "$c/MESH" "$c/input.ini" "$v/"
  sed -i "$3" "$v/input.ini"
  local n; n=$(diff "$c/input.ini" "$v/input.ini" | grep -c '^>')
  [ "$n" = 1 ] && ok "$leg: the variant differs from its deck in one line ($(diff "$c/input.ini" "$v/input.ini" | grep '^>' | cut -c3-))" \
               || ko "$leg: the variant differs in $n lines"
  OMP_NUM_THREADS=1 "$RUN" "$c" "$T/$leg-base" 1 "$T/$leg-base.log" "$EXE" > /dev/null 2>&1 &
  OMP_NUM_THREADS=1 "$RUN" "$v" "$T/$leg-variant" 1 "$T/$leg-variant.log" "$EXE" > /dev/null 2>&1 &
  wait
  for r in base variant; do
    grep -q "GAME OVER" "$T/$leg-$r.log" && ok "$leg: the $r run ends with GAME OVER" \
      || ko "$leg: the $r run: $(grep -a -m1 'ERROR\|Error in block\|NaN' "$T/$leg-$r.log" | cut -c1-200)"
  done
  local f d=0
  for f in "$T/$leg-base"/OUTPUT/*-residual-history.dat; do
    cmp -s "$f" "$T/$leg-variant/OUTPUT/$(basename "$f")" || d=$((d+1))
  done
  [ $d -gt 0 ] && ok "$leg: the knob changes the steps ($d residual histories differ)" \
               || ko "$leg: the knob changes no step (every residual history is the same)"
}
judge() {   # judge <leg> <tol> [files]: only when both runs ended (a stopped run has no converged field)
  local leg=$1 tol=$2; shift 2
  if ! grep -q "GAME OVER" "$T/$leg-base.log" || ! grep -q "GAME OVER" "$T/$leg-variant.log"; then
    ko "$leg: no converged field to compare (a run did not end with GAME OVER)"; return
  fi
  local o; o=$(python3 -B "$CMP" "$T/$leg-base" "$T/$leg-variant" "$tol" "$@" 2>&1); local rc=$?
  echo "$o" | sed 's/^/      /'
  [ $rc = 0 ] && ok "$leg: the converged field does not depend on the step ($(echo "$o" | tail -1 | cut -d' ' -f2) <= $tol)" \
              || ko "$leg: the converged field depends on the step ($(echo "$o" | tail -1) > $tol)"
}
full() {    # full <leg> <file> <var>: the family is not empty
  local p; p=$(python3 -B "$FLD" peak "$T/$1-base/OUTPUT/$2" "$3" 2>/dev/null)
  python3 -c "import sys; sys.exit(0 if float('${p:-0}') > 1e-4 else 1)" && ok "$1: the smoke family holds $p kg/m3 at its peak" \
    || ko "$1: the smoke family is empty (${p:-no field})"
}
balance() { # balance <leg> <p31_balance.py arguments>: the converged field's balance
  local leg=$1; shift
  local o; o=$(python3 -B "$BAL" "$@" 2>&1); local rc=$?
  echo "$o" | sed 's/^/      /'
  [ $rc = 0 ] && ok "$leg: $(echo "$o" | tail -1 | cut -d, -f1) balances ($(echo "$o" | tail -1 | grep -o 'relative [^ ]*'))" \
              || ko "$leg: the balance fails ($(echo "$o" | tail -1))"
}
burns() {   # burns <leg>: drops below d-burnout fill cells, and the coupler hands the gas a burn
  local s; s=$(python3 -B "$FLD" small "$T/$1-base/OUTPUT/part-field.tec" 1 1000 1.7e-5 2>/dev/null)
  [ "${s%% *}" -gt 0 ] 2>/dev/null && ok "$1: $(echo $s | cut -d' ' -f1) cells hold drops at or below d-burnout (smallest $(echo $s | cut -d' ' -f2) m)" \
    || ko "$1: no cell holds drops at or below d-burnout (${s:-no field})"
  local r; r=$(grep -a 'ice-to-gas-mass' "$T/$1-base.log" | tail -1 | awk '{print $3}')
  python3 -c "import sys; sys.exit(0 if float('${r:-0}') > 0 else 1)" && ok "$1: the coupler hands the gas $r kg/s" \
    || ko "$1: no burn handed to the gas (${r:-none})"
}

for leg in $LEGS; do
  case $leg in
    births)     echo "--- births: ICE's smoke births, ICE tau-factor 1 vs 0.5 ---"
                pair births "$H/test/steady-ice-smoke" "$TAU"; judge births 1e-12; full births part-field.tec rho_p2 ;;
    births-fo)  echo "--- births-fo: the same at first order ---"
                pair births-fo "$(deck births-fo steady-ice-smoke "$FO")" "$TAU"; judge births-fo 1e-12; full births-fo part-field.tec rho_p2
                balance births-fo smoke "$T/births-fo-base" "$T/births-fo-base.log" part-field.tec 2 1e-9 ;;
    removal)    echo "--- removal: the steady burn-out (f of what leaves a cell in the band burns at its exit), ICE cfl 0.2 vs 0.1 ---"
                pair removal "$H/test/steady-ice-burnout" "$ICECFL"; judge removal 1e-12; burns removal ;;
    removal-fo) echo "--- removal-fo: the same at first order ---"
                pair removal-fo "$(deck removal-fo steady-ice-burnout "$FO")" "$ICECFL"; judge removal-fo 1e-12; burns removal-fo
                balance removal-fo burnout "$T/removal-fo-base" "$T/removal-fo-base.log" 1e-13 ;;
    igloo)      echo "--- igloo: IGLOO's smoke born in ICE, ICE tau-factor 1 vs 0.5 ---"
                pair igloo "$H/test/steady-igloo-smoke" "$TAU"; judge igloo 1e-7; full igloo partE-field.tec rho_p1 ;;
    igloo-fo)   echo "--- igloo-fo: the same at first order ---"
                pair igloo-fo "$(deck igloo-fo steady-igloo-smoke "$FO")" "$TAU"; judge igloo-fo 1e-8; full igloo-fo partE-field.tec rho_p1
                balance igloo-fo igloo "$T/igloo-fo-base" "$T/igloo-fo-base.log" partE-field.tec 1 1e-13 ;;
    mixed)      echo "--- mixed: MOSE steady, ICE time-accurate in alternate cycles, MOSE cfl 0.5 vs 0.25 ---"
                pair mixed "$(deck mixed steady-ice-burnout "$MIX; $FO")" "$MOSECFL"; judge mixed 1e-1 gas-field.tec
                r=$(grep -a 'ice-removed-rate' "$T/mixed-base.log" | tail -1 | awk '{print $3}')
                python3 -c "import sys; sys.exit(0 if float('${r:-0}') > 0 else 1)" && ok "mixed: the coupler hands a removal of $r kg/s" \
                  || ko "mixed: no removal reported (${r:-none})" ;;
    burn-fo)    echo "--- burn-fo: K-burn x100 at first order, the drops burning out in the first cells, ICE cfl 0.2 vs 0.1 ---"
                pair burn-fo "$(deck burn-fo steady-ice-burnout "$FO" "$K100")" "$ICECFL"; judge burn-fo 1e-12; burns burn-fo
                balance burn-fo burnout "$T/burn-fo-base" "$T/burn-fo-base.log" 1e-5 ;;
    d0-fo)      echo "--- d0-fo: the same without d-burnout, the removal at m_burnout (d -> 0): one run ---"
                c=$(deck d0-fo steady-ice-burnout "$FO" "$K100; $NODB")
                n=$(grep -c 'd-burnout' "$c/INPUT/part-models.txt")
                [ "$n" = 0 ] && ok "d0-fo: the deck has no d-burnout" || ko "d0-fo: d-burnout still in the models file"
                rm -rf "$T/d0-fo-base"
                OMP_NUM_THREADS=1 "$RUN" "$c" "$T/d0-fo-base" 1 "$T/d0-fo-base.log" "$EXE" > /dev/null 2>&1
                grep -q "GAME OVER" "$T/d0-fo-base.log" && ok "d0-fo: the run ends with GAME OVER" \
                  || ko "d0-fo: the run: $(grep -a -m1 'ERROR\|Error in block\|NaN' "$T/d0-fo-base.log" | cut -c1-200)"
                grep -q 'negative densities floored' "$T/d0-fo-base.log" && ko "d0-fo: densities floored ($(grep -a 'negative densities floored' "$T/d0-fo-base.log" | tail -1 | cut -c1-120))" \
                  || ok "d0-fo: no density floored"
                balance d0-fo burnout "$T/d0-fo-base" "$T/d0-fo-base.log" 1e-5
                f=$(python3 -B "$BAL" burnout "$T/d0-fo-base" "$T/d0-fo-base.log" 1 2>/dev/null | head -1 | awk '{print $3}' | tr -d ';')
                [ "${f:-0}" -gt 0 ] 2>/dev/null && ok "d0-fo: $f cells burn their drops out at m_burnout" || ko "d0-fo: no cell flagged (${f:-none})" ;;
    agree)      echo "--- agree: removal-fo's steady run against mixed's time-accurate one (the outlet judged, the hand-off stated) ---"
                if [ -f "$T/removal-fo-base.log" ] && [ -f "$T/mixed-base.log" ]; then
                  ao=$(python3 -B "$BAL" agree "$T/removal-fo-base" "$T/removal-fo-base.log" "$T/mixed-base.log" 2>&1); arc=$?
                  echo "$ao" | sed 's/^/      /'
                  [ $arc = 0 ] && ok "agree: the metal leaving the steady run's outlet is within 1e-6 of its hand-off ($(echo "$ao" | grep -o 'outflow records [^ ]* kg/s'))" \
                               || ko "agree: $(echo "$ao" | tail -1)"
                  [ $arc = 0 ] && echo "STATE agree: $(echo "$ao" | tail -1)"
                else
                  ko "agree: run removal-fo and mixed first"
                fi ;;
    agree-d0)   echo "--- agree-d0: d0-fo's steady run against its time-accurate twin (the outlet judged, the hand-off stated) ---"
                if [ -f "$T/d0-fo-base.log" ]; then
                  c=$(deck agree-d0-ta steady-ice-burnout "$MIX120; $FO" "$K100; $NODB")
                  rm -rf "$T/agree-d0-ta"
                  OMP_NUM_THREADS=1 "$RUN" "$c" "$T/agree-d0-ta" 1 "$T/agree-d0-ta.log" "$EXE" > /dev/null 2>&1
                  grep -q "GAME OVER" "$T/agree-d0-ta.log" && ok "agree-d0: the ta run ends with GAME OVER" \
                    || ko "agree-d0: the ta run: $(grep -a -m1 'ERROR\|Error in block\|NaN' "$T/agree-d0-ta.log" | cut -c1-200)"
                  ao=$(python3 -B "$BAL" agree "$T/d0-fo-base" "$T/d0-fo-base.log" "$T/agree-d0-ta.log" 2>&1); arc=$?
                  echo "$ao" | sed 's/^/      /'
                  [ $arc = 0 ] && ok "agree-d0: the metal leaving the steady run's outlet is within 1e-6 of its hand-off" || ko "agree-d0: $(echo "$ao" | tail -1)"
                  [ $arc = 0 ] && echo "STATE agree-d0: $(echo "$ao" | tail -1)"
                else
                  ko "agree-d0: run d0-fo first"
                fi ;;
    agree-d0k10) echo "--- agree-d0k10: the same at K-burn x10 (the burn-out spread over cells) ---"
                twin agree-d0k10 "$(deck agree-d0k10-st steady-ice-burnout "$FO" "$K10; $NODB")" \
                  "$(deck agree-d0k10-ta steady-ice-burnout "$MIX120; $FO" "$K10; $NODB")"
                agreein agree-d0k10 "$T/agree-d0k10-st" "$T/agree-d0k10-ta.log" ;;
    agree-w200) echo "--- agree-w200: the shipped deck at first order, ICE's wall taking no drops (the models coincide but for the band) ---"
                cs=$(deck agree-w200-st steady-ice-burnout "$FO"); ct=$(deck agree-w200-ta steady-ice-burnout "$MIX120; $FO")
                ns=$(nodeposit "$cs"); nt=$(nodeposit "$ct")
                [ "$ns" -gt 0 ] && [ "$ns" = "$nt" ] && ok "agree-w200: $ns wall records of ICE set to 200 in each deck" \
                  || ko "agree-w200: wall records set: $ns and $nt"
                twin agree-w200 "$cs" "$ct"
                agreein agree-w200 "$T/agree-w200-st" "$T/agree-w200-ta.log" ;;
    band)       echo "--- band: ICE's burnout_band, one value in the log, the source and the texts that name it ---"
                if [ -f "$T/removal-fo-base.log" ]; then
                  bo=$(python3 -B "$BAL" band "$T/removal-fo-base.log" "$H/src/lib/ICE/src/lib/numerics/fluxes/bc/Mod_BC_Fluxes.f90" \
                    "$H/test/steady-ice-burnout/README.md" "$H/src/lib/ICE/docs/theory/physics.md" "$H/test/mi2-checks/p31_steady_check.sh" 2>&1)
                  brc=$?
                  echo "$bo" | sed 's/^/      /'
                  [ $brc = 0 ] && ok "band: $(echo "$bo" | tail -1)" || ko "band: $(echo "$bo" | tail -1)"
                else
                  ko "band: run removal-fo first"
                fi ;;
    *) ko "unknown leg $leg" ;;
  esac
done

echo "== p31: $pass passed, $fail failed"
[ $fail = 0 ]
