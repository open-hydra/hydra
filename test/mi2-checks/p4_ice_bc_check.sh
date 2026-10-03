#!/bin/bash
# P4 checks on the ICE BC reader (JPL-Eulerian, short schedule, 1 thread):
#  301 : face-4 records 300 -> 301 (ATLAS dispersed wall) -> run OK and OUTPUT byte-identical to the 300 run
#        (ICE's type 3 already absorbs toward-moving particles: 300 and 301 coincide in ICE)
#  403 : face-1 records 401 -> 403 with gp = 94.105 kg/m2/s (= krho/(1-krho)*mean(rho_g u_g)|inlet of the CONVERGED
#        401 run) -> the inlet-column condensed flux must equal gp within 1 %. The 401 run is checked against ITS OWN
#        prescription, krho/(1-krho)*mean(rho_g u_g) from the same run's gas field (a short run is not converged to
#        94.1, so the two inlets are not compared with each other).
#  403z: same with gp = 0 -> the inlet column drains (mass flux < 50 % of the 401 run). Before P4, 402/403 mapped
#        to ICE type 0 = NO boundary condition, i.e. exactly this "open face" behaviour.
#  grp : INPUT/<prefix>phase.txt declaring 2 populations while [ICE-Family*] has 1 -> WARNING printed; shipped -> none
#   usage: p4_ice_bc_check.sh [scratch-root] [exe]
set -u
export LC_ALL=C
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P4}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
S=$H/test/JPL-Eulerian-20micron
EXE=${2:-$H/bin/hydra_MI2}
RUN=$H/test/mi2-checks/run_case.sh
TEC=$H/test/mi2-checks/tec_cc.py
GP=94.105
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
mk() { local d=$T/$1; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"
       sed -i -e 's/^N-coupling *=.*/N-coupling    = 2/' -e 's/^MOSE-iter *=.*/MOSE-iter     = 100/' \
              -e 's/^ICE-iter *=.*/ICE-iter      = 50/' -e 's/^sol-diter *=.*/sol-diter   = 50/' "$d/input.ini"; echo "$d"; }
flux() { python3 - "$1" "$H/test/mi2-checks" <<'EOF'
import sys; sys.path.insert(0, sys.argv[2])
from tec_cc import read_cc
p = read_cc(sys.argv[1]); print('%.6e' % (p['rho_p1'][0,:,0]*p['u_p1'][0,:,0]).mean())
EOF
}
gasflux() { python3 - "$1" "$H/test/mi2-checks" <<'EOF'
import sys; sys.path.insert(0, sys.argv[2])
from tec_cc import read_cc
g = read_cc(sys.argv[1]); print('%.6e' % (0.3/0.7*(g['rho(1)'][0,:,0]*g['u'][0,:,0]).mean()))
EOF
}

# --- reference (401, face-4 300) ---
d=$(mk ref401); "$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
[ $rc = 0 ] && ok "reference 401 run (exit 0)" || ko "reference run exit=$rc"
f401=$(flux "$d/run/OUTPUT/part-field.tec"); e401=$(gasflux "$d/run/OUTPUT/gas-field.tec")
echo "  401 inlet-column condensed flux: $f401   prescribed krho/(1-krho)*rho_g*u_g: $e401"
python3 -c "import sys; a,b=float('$f401'),float('$e401'); sys.exit(0 if abs(a-b)/b < 0.02 else 1)" && ok "401 inlet reproduces its krho prescription within 2 %" || ko "401 inlet flux off its prescription"
grep -q "WARNING.*phase.txt declares" "$d/run.log" && ko "spurious group-count warning on the shipped case" || ok "no group-count warning on the shipped case"

# --- 301 ---
d=$(mk wall301); awk '{ if ($5==4 && $6==300) $6=301; if (NF==6) printf "%8d%8d%8d%8d%8d%8d\n",$1,$2,$3,$4,$5,$6; else print }' "$S/INPUT/part-bc.txt" > "$d/INPUT/part-bc.txt"
echo "  301 records: $(awk '$6==301' "$d/INPUT/part-bc.txt" | wc -l)"
"$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
[ $rc = 0 ] && ok "301 (dispersed wall) accepted by ICE" || { ko "301 run exit=$rc"; grep -n "unsupported\|ERROR" "$d/run.log" | head -2; }
diff -rq "$T/ref401/run/OUTPUT" "$d/run/OUTPUT" > /dev/null && ok "301 run byte-identical to the 300 run" || ko "301 run differs from the 300 run"

# --- 403 equivalent gp ---
d=$(mk inlet403); awk -v gp=$GP -v tp=555.56 'BEGIN{pl=0} { if (pl==1) { $1=sprintf("%.5E",gp); $5=sprintf("%.5E",tp); printf "  %13s  %13s %15s %15s  %13s  %13s  %13s %s  %13s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9; pl=0; next }
  if ($5==1 && $6==401) { $6=403; pl=1 } if (NF==6) printf "%8d%8d%8d%8d%8d%8d\n",$1,$2,$3,$4,$5,$6; else print }' "$S/INPUT/part-bc.txt" > "$d/INPUT/part-bc.txt"
sed -n 1,2p "$d/INPUT/part-bc.txt"
"$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
[ $rc = 0 ] && ok "403 inlet accepted by ICE" || { ko "403 run exit=$rc"; grep -n "unsupported\|ERROR" "$d/run.log" | head -2; }
f403=$(flux "$d/run/OUTPUT/part-field.tec"); echo "  403 (gp=$GP) inlet-column condensed flux: $f403"
python3 -c "import sys; a,b=float('$f403'),$GP; sys.exit(0 if abs(a-b)/b < 0.01 else 1)" && ok "403 inlet reproduces gp within 1 % (was an open face before P4)" || ko "403 inlet flux off gp"

# --- 403 with gp = 0 (open face) ---
d=$(mk inlet403zero); sed 's/^    9.41050E+01 /    0.00000E+00 /' "$T/inlet403/INPUT/part-bc.txt" > "$d/INPUT/part-bc.txt"
echo "  gp=0 payload lines: $(grep -c '^    0.00000E+00 ' "$d/INPUT/part-bc.txt")"
"$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
f403z=$(flux "$d/run/OUTPUT/part-field.tec"); echo "  403 (gp=0) inlet-column condensed flux: $f403z"
python3 -c "import sys; a,b=float('$f401'),float('$f403z'); sys.exit(0 if b < 0.5*a else 1)" && ok "gp = 0 drains the inlet column (the pre-P4 'no BC' behaviour is not what 403 does now)" || ko "gp = 0 did not drain the inlet column"

# --- group-count warning ---
d=$(mk groups); sed -i '2s/ 1$/ 2/' "$d/INPUT/part-phase.txt"; sed -n 1,2p "$d/INPUT/part-phase.txt"
"$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
grep -q "phase.txt declares 2 population(s) over 1 material(s), but \[ICE-Family\*\] defines" "$d/run.log" && ok "group-count mismatch warned" || ko "group-count mismatch not warned"

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
