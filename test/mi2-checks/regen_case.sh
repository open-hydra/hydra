#!/bin/bash
# regen_case.sh <case-dir> <scratch-dir> <tools: GPB|BCB|ICB ...> [--install]
# Runs the ATLAS tools on a scratch copy of the case (INPUT/*phase.txt, *thermo.dat, *transport.dat,
# *properties.dat and MESH/mesh.tec copied to the case root as figaro's ATLAS.sh does), then diffs
# fromATLAStoSolver/ against the case's INPUT/. With --install the regenerated files are copied into
# INPUT/ (a backup of the previous INPUT/ is left as INPUT.pre-regen/ unless one already exists).
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
cp INPUT/*phase.txt INPUT/*thermo.dat INPUT/*transport.dat INPUT/*properties.dat . 2>/dev/null
cp "$src/MESH/mesh.tec" . 2>/dev/null || cp "$src/MESH/mesh.p3d" . 2>/dev/null
"$ATLASDIR/ATLAS.sh" "${tools[@]}" > atlas.log 2>&1; rc=$?
if [ $rc -ne 0 ]; then echo "ATLAS exit $rc -- see $dst/atlas.log"; tail -5 atlas.log; exit $rc; fi
echo "== $(basename "$src"): ${tools[*]} (ATLAS $ATLASDIR)"
for f in fromATLAStoSolver/*; do
  b=$(basename "$f")
  if [ ! -e "INPUT/$b" ]; then echo "   $b: NEW (not in INPUT)"
  elif cmp -s "$f" "INPUT/$b"; then echo "   $b: identical"
  else echo "   $b: DIFFERS ($(diff "$f" "INPUT/$b" | grep -c '^[<>]') diff lines)"; diff "$f" "INPUT/$b" | head -4 | sed 's/^/      /'; fi
done
if [ $install = 1 ]; then
  [ -d "$src/INPUT.pre-regen" ] || cp -rp "$src/INPUT" "$src/INPUT.pre-regen"
  cp -p fromATLAStoSolver/* "$src/INPUT/"
  echo "   installed into $src/INPUT (backup INPUT.pre-regen)"
fi
