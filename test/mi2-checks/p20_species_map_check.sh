#!/bin/bash
# Species-map check (MI2_Stoichiometry.f90 Load_Species_Map): a material that evaporates gets its gas species from
# INPUT/<phase>products.txt, and a file that does not say them stops the run before the first step.
#  control : evap-box-ta as shipped runs, and its log names the file and the process
#  each row below: evap-box-ta with ONE edit must stop with exit /= 0 and the message
#  lv      : evap-box-ta gives the latent heat as the Lv column of the property table: no message. The same case with
#            the column taken out and a constant in [IGLOO-Properties] runs, and the coupler prints one warning with
#            the size of what the constant misses at the boiling temperature (against the two tables, here).
#   usage: p20_species_map_check.sh [scratch-root] [exe]
set -u
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P20}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=${2:-$H/bin/hydra_MI2}
S=$H/test/evap-box-ta
RUN=$H/test/mi2-checks/run_case.sh
mkdir -p "$T"; pass=0; fail=0
ok() { pass=$((pass+1)); echo "PASS  $1"; }
ko() { fail=$((fail+1)); echo "FAIL  $1"; }
mk() { local d=$T/$1; rm -rf "$d"; mkdir -p "$d"; cp -rp "$S/INPUT" "$S/MESH" "$S/input.ini" "$d/"; echo "$d"; }
F=INPUT/part-products.txt

d=$(mk control)
"$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
[ $rc = 0 ] && grep -q "IGLOO material -> gas species .*$F" "$d/run.log" && grep -q -- "- A \[evaporation\] :" "$d/run.log" \
  && ok "control: the shipped case runs and reports the file and the process" || ko "control: exit=$rc, or the report lines are missing"

# refuse <name> <message> <sed script on the products file | a shell command run in the case directory>
refuse() {
  local d; d=$(mk "$1")
  ( cd "$d" && eval "$3" ) || { ko "$1: the edit failed"; return; }
  if cmp -s "$d/$F" "$S/$F" && cmp -s "$d/input.ini" "$S/input.ini"; then ko "$1: the edit changed nothing"; return; fi
  "$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
  if [ $rc != 0 ] && grep -q -F "$2" "$d/run.log" && ! grep -q "GAME OVER" "$d/run.log"; then ok "$1: refused ($2)"
  else ko "$1: exit=$rc, message \"$2\" $(grep -q -F "$2" "$d/run.log" && echo printed || echo absent)"; fi
}
refuse missing-file   "is missing: run ATLAS GPB"                    "rm $F"
refuse other-material "names the material B, which IGLOO does not"   "sed -i 's/^material A /material B /' $F"
refuse other-species  "names the gas species CO2, which the gas"     "sed -i 's/^  H2O /  CO2 /' $F"
refuse other-gas      "was written for another gas phase"            "sed -i 's/^gas gas- 2\$/gas gas- 3/' $F"
refuse no-record      "has no evaporation record for it"             "sed -i -e 's/^material A evaporation 1 latent .*/material A none 0 none 0.0000000000e+00/' -e '/^  H2O /d' $F"
refuse row-sum        "does not sum to 1"                            "sed -i 's/^  H2O 1.0000000000e+00 /  H2O 9.0000000000e-01 /' $F"
refuse energy-word    "an evaporation record that is not"            "sed -i 's/ evaporation 1 latent / evaporation 1 formation /' $F"
refuse truncated      "fewer terms than its last record announces"   "sed -i '/^  H2O /d' $F"
refuse retired-key    "products-A is no longer read"                 "sed -i 's/^\[MI2-Coupling\]\$/[MI2-Coupling]\nproducts-A = H2O:1.0/' input.ini"
refuse euler-cycle    "a cycle has IGLOO-iter = 2"                   "sed -i 's/^IGLOO-iter    = 1\$/IGLOO-iter    = 2/' input.ini"
refuse old-word       "an evaporation record that is not"            "sed -i 's/ evaporation 1 latent / evaporation 1 explicit /' $F"

# --- lv: the latent heat as a column (no message) and as a constant (one warning, its size checked) ---
W='the latent heat is a constant'
n=$(grep -c "$W" "$T/control/run.log")
[ "$n" = 0 ] && ok "lv: the shipped case (Lv column) prints no warning" || ko "lv: the shipped case prints $n warning(s)"
d=$(mk lv-constant)
python3 - "$d" <<'EOF'
import re, sys
d = sys.argv[1]
p = d + '/INPUT/part-properties.dat'
out = []
for line in open(p):
    if line.startswith('VARIABLES'):
        names = re.findall(r'"([^"]*)"', line)
        k = names.index('Lv')
        out.append('VARIABLES = ' + ', '.join('"%s"' % v for i, v in enumerate(names) if i != k) + '\n')
    else:
        f = line.split()
        try:
            [float(x) for x in f]
            out.append(' '.join(v for i, v in enumerate(f) if i != k) + ' \n' if len(f) == len(names) else line)
        except ValueError:
            out.append(line)
open(p, 'w').write(''.join(out))
t = open(d + '/input.ini').read()
assert t.count('[IGLOO-Properties]\n') == 1
open(d + '/input.ini', 'w').write(t.replace('[IGLOO-Properties]\n', '[IGLOO-Properties]\nLv    = 2.4417e6\n'))
EOF
if ! grep -q '"Lv"' "$d/INPUT/part-properties.dat" && grep -q '^Lv    = 2.4417e6$' "$d/input.ini"; then
  "$RUN" "$d" "$d/run" 1 "$d/run.log" "$EXE" > /dev/null 2>&1; rc=$?
  n=$(grep -c "$W" "$d/run.log")
  [ $rc = 0 ] && grep -q "GAME OVER" "$d/run.log" && [ "$n" = 1 ] && ok "lv: a constant latent heat runs and the coupler prints one warning"                                                                || ko "lv: exit=$rc, $n warning(s)"
  out=$(python3 - "$S" "$d/run.log" <<'EOF'
import re, sys
case, log = sys.argv[1], sys.argv[2]
def zones(path):
    out, cur, names = {}, None, []
    for line in open(path):
        if line.startswith('VARIABLES'):
            names = re.findall(r'"([^"]*)"', line)
        elif line.startswith('ZONE'):
            cur = out.setdefault(re.search(r'T="([^"]*)"', line).group(1), {})
        elif cur is not None:
            try:
                v = [float(x) for x in line.split()]
            except ValueError:
                continue
            if v:
                cur[v[0]] = v
    return out, names
def at(rows, T, col):
    lo = float(int(T))
    return rows[lo][col] if T == lo else rows[lo][col] + (rows[lo + 1.0][col] - rows[lo][col]) * (T - lo)
m = re.search(r'the latent heat is a constant: at\s+(\S+) K it differs by\s+(\S+) J/kg \(\s*(\S+) % of Lv\)', open(log, errors='replace').read())
if not m:
    print('FAIL  lv: the warning does not read "at <T> K it differs by <x> J/kg (<p> % of Lv)"'); sys.exit()
Tb, off = float(m.group(1)), float(m.group(2))
pt, pn = zones(case + '/INPUT/part-properties.dat'); gt, gn = zones(case + '/INPUT/gas-thermo.dat')
hc = [i for i, v in enumerate(pn) if v.startswith('Enthalpy')][0]
an = (at(pt['A'], Tb, hc) - at(pt['A'], 298.15, hc)) - (at(gt['H2O'], Tb, 2) - at(gt['H2O'], 298.15, 2))
e = abs(off - an) / abs(an)
print('%s  lv: the size printed, %.4e J/kg at %.2f K, against the two tables, %.6e J/kg: rel.err %.1e (tol 1e-4: five printed figures)'
      % ('PASS' if e <= 1e-4 and abs(an) > 1e4 else 'FAIL', off, Tb, an, e))
EOF
)
  echo "$out"; pass=$((pass + $(echo "$out" | grep -c '^PASS'))); fail=$((fail + $(echo "$out" | grep -c '^FAIL')))
else ko "lv: the variant could not be made"; fi

echo "== $pass passed, $fail failed"; [ $fail = 0 ]
