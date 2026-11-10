#!/bin/bash
# regen_case.sh <case-dir> <scratch-dir> <tools: GPB|BCB|ICB ...> [--install]
# Runs the ATLAS tools on a scratch copy of the case (MESH/mesh.tec copied to the case root; without GPB among the
# tools also INPUT/*phase.txt, *thermo.dat, *transport.dat and *properties.dat, which BCB and ICB then read), then
# diffs fromATLAStoSolver/ against the case's INPUT/. The ATLAS that wrote the files (commit, deps.lock, the ORION it
# was configured with, the md5 of BCB and ICB, the Fortran compiler of its build, the Python that ran GPB) is printed
# and kept as atlas-provenance.txt in the scratch copy. With --install the regenerated files are copied into INPUT/,
# the previous INPUT/ kept as INPUT.pre-regen-<date>-<tag>/ (checked against INPUT/ before anything is installed)
# with the provenance beside it.
# ATLASDIR selects the ATLAS tree (default: the checked-out utils/ATLAS); ATLAS_PYTHON bypasses conda.
# Byte identity with a stored INPUT/ needs ATLAS built with the compiler that wrote it (gfortran for the cases of this
# repository): an ifx build of ICB moves the *-ic.tec values by about 5e-9 relative.
set -u
usage="usage: regen_case.sh <case-dir> <scratch-dir> <tools: GPB|BCB|ICB ...> [--install]"
src=$(cd "${1:?$usage}" && pwd); dst=${2:?$usage}; shift 2
install=0; tools=()
for a in "$@"; do [ "$a" = --install ] && install=1 || tools+=("$a"); done
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
ATLASDIR=${ATLASDIR:-$H/utils/ATLAS}
export ATLASDIR
[ -n "${ATLAS_PYTHON:-}" ] && export ATLAS_PYTHON

rm -rf "$dst"; mkdir -p "$dst"
cp -rpL "$src/INPUT" "$dst/"; cp -p "$src/input.ini" "$dst/"
cd "$dst" || exit 2
# BCB and ICB read a phase's tables from the case root when its phase.txt is there, else from fromATLAStoSolver/
# (ATLAS phase_dir): after GPB they must read GPB's tables, so the stored ones are copied only without GPB
gpb=0; for t in "${tools[@]}"; do [ "$t" = GPB ] && gpb=1; done
[ $gpb = 1 ] || cp INPUT/*phase.txt INPUT/*thermo.dat INPUT/*transport.dat INPUT/*properties.dat . 2>/dev/null
cp "$src/MESH/mesh.tec" . 2>/dev/null || cp "$src/MESH/mesh.p3d" . 2>/dev/null
"$ATLASDIR/ATLAS.sh" "${tools[@]}" > atlas.log 2>&1; rc=$?
if [ $rc -ne 0 ]; then echo "ATLAS exit $rc -- see $dst/atlas.log"; tail -5 atlas.log; exit $rc; fi
{ echo "ATLASDIR $ATLASDIR"
  m=$(git -C "$ATLASDIR" status --porcelain --untracked-files=no 2>/dev/null)
  echo "commit $(git -C "$ATLASDIR" rev-parse HEAD 2>/dev/null || echo '?')${m:+ (modified)}"
  if [ -s "$ATLASDIR/deps.lock" ]; then sed 's/^/deps.lock /' "$ATLASDIR/deps.lock"
  else echo "deps.lock absent or empty (a local ORION is not locked)"; fi
  rev() { local m; m=$(git -C "$1" status --porcelain --untracked-files=no 2>/dev/null)
          echo "commit $(git -C "$1" rev-parse HEAD 2>/dev/null || echo '?') ($(git -C "$1" describe --tags --always \
            2>/dev/null || echo '?')${m:+, modified})"; }
  o=$(sed -n 's/^FETCHCONTENT_SOURCE_DIR_ORION:[A-Z]*=//p' "$ATLASDIR/build/CMakeCache.txt" 2>/dev/null)
  [ -n "$o" ] && echo "ORION local $o $(rev "$o")"
  [ -d "$ATLASDIR/lib/ORION" ] && echo "ORION lib/ORION $(rev "$ATLASDIR/lib/ORION")"
  for b in BCB ICB; do [ -f "$ATLASDIR/bin/$b" ] && echo "bin/$b $(md5sum < "$ATLASDIR/bin/$b" | cut -c1-32)"; done
  fc=$(sed -n 's/^CMAKE_Fortran_COMPILER:[A-Z]*=//p' "$ATLASDIR/build/CMakeCache.txt" 2>/dev/null)
  echo "Fortran ${fc:-?} ($("${fc:-false}" --version 2>/dev/null | head -1))"
  # the Python GPB ran: ATLAS_PYTHON, else python3 of ATLAS's conda environment (ATLAS.sh)
  py=${ATLAS_PYTHON:-}
  if [ -z "$py" ] && command -v conda > /dev/null 2>&1; then
    py=$(conda env list 2>/dev/null | awk -v e="${ATLAS_CONDA_ENV:-ct-env}" '$1 == e {print $NF}')/bin/python3
  fi
  if [ -n "$py" ] && "$py" -c 1 2> /dev/null; then
    echo "python $("$py" -c 'import sys
v = []
for n in ("cantera", "numpy", "CoolProp"):
    try:
        v.append("%s %s" % (n, __import__(n).__version__))
    except Exception:
        v.append(n + " absent")
print(sys.executable, sys.version.split()[0] + ";", ", ".join(v))')"
  else echo "python ? (GPB's interpreter not found)"; fi
  echo "tools ${tools[*]}"; } > atlas-provenance.txt
echo "== $(basename "$src"): ${tools[*]} (ATLAS $ATLASDIR)"
sed 's/^/   /' atlas-provenance.txt
for f in fromATLAStoSolver/*; do
  b=$(basename "$f")
  if [ ! -e "INPUT/$b" ]; then echo "   $b: NEW (not in INPUT)"
  elif cmp -s "$f" "INPUT/$b"; then echo "   $b: identical"
  else echo "   $b: DIFFERS ($(diff "$f" "INPUT/$b" | grep -c '^[<>]') diff lines)"; diff "$f" "INPUT/$b" | head -4 | sed 's/^/      /'; fi
done
if [ $install = 1 ]; then
  bak=$(mktemp -d "$src/INPUT.pre-regen-$(date +%Y%m%d-%H%M%S)-XXXX") \
    || { echo "   no backup directory: nothing installed"; exit 2; }
  if ! cp -rp "$src/INPUT/." "$bak/" || ! diff -rq "$src/INPUT" "$bak" > /dev/null; then
    echo "   $bak is not a copy of INPUT/: nothing installed"; exit 2
  fi
  cp -p fromATLAStoSolver/* "$src/INPUT/"
  cp -p atlas-provenance.txt "$bak.atlas-provenance.txt"
  echo "   installed into $src/INPUT (backup $(basename "$bak"), provenance $(basename "$bak").atlas-provenance.txt)"
fi
