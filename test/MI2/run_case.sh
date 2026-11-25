#!/bin/bash
# run_case.sh <src-case-dir> <scratch-dir> <nthreads> <logfile> [exe]
# Copies INPUT/, MESH/, input.ini (and JPL_chang.dat) into <scratch-dir>, runs hydra_MI2 there.
# The log header records the executable's mtime and the IGLOO / ICE commits currently checked out
# (read from their git ref files -- what the tree holds, not necessarily what the binary was linked
# against: rebuild before trusting the stamp).
set -u
usage="usage: run_case.sh <src-case-dir> <scratch-dir> <nthreads> <logfile> [exe]"
src=${1:?$usage}; dst=${2:?$usage}; nt=${3:?$usage}; log=${4:?$usage}
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
exe=${5:-$H/bin/hydra_MI2}

head_of() {  # head_of <checkout>  -> short sha of HEAD (follows a gitdir file and a symbolic ref; '?' if unreadable)
  local g=$1/.git ref
  [ -f "$g" ] && g=$(cd "$1" 2>/dev/null && cd "$(sed -n 's/^gitdir: //p' .git)" 2>/dev/null && pwd)
  [ -f "$g/HEAD" ] || { echo '?'; return; }
  ref=$(cat "$g/HEAD")
  case "$ref" in
    ref:*) ref=${ref#ref: }; [ -f "$g/$ref" ] && cut -c1-7 "$g/$ref" || echo '?';;
    *)     echo "${ref:0:7}";;
  esac
}
igloo=$(head_of "$H/src/lib/IGLOO")
ice=$(head_of "$H/src/lib/ICE")

mkdir -p "$dst"
cp -rpL "$src/INPUT" "$src/input.ini" "$dst/"
[ -d "$src/MESH" ] && cp -rpL "$src/MESH" "$dst/"
[ -f "$src/JPL_chang.dat" ] && cp -p "$src/JPL_chang.dat" "$dst/"
cd "$dst" || exit 2
rm -rf OUTPUT; mkdir -p OUTPUT
ls INPUT/*phase.txt > filelist.txt
ulimit -s unlimited
export KMP_STACKSIZE=100M
export OMP_NUM_THREADS=$nt
echo "# exe=$exe ($(stat -c %y "$exe")) IGLOO-tree=$igloo ICE-tree=$ice threads=$nt start=$(date -Is)" > "$log"
if [ -x /usr/bin/time ]; then
  /usr/bin/time -f "# wall=%e s  maxrss=%M kB" "$exe" >> "$log" 2>&1
else
  "$exe" >> "$log" 2>&1
fi
rc=$?
echo "# exit=$rc end=$(date -Is)" >> "$log"
exit $rc
