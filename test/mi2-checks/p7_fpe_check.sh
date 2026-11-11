#!/bin/bash
# p7_fpe_check.sh [scratch-root]
# FPE cleanliness of the coupling paths: builds hydra_MI2 with -fpe0 -traceback (RELEASE flags otherwise) in a
# scratch build tree and runs test/evap-box-ord2 (gas-order = 2, 1 thread). PASS = GAME OVER; a trap prints the
# traceback. RED on the pre-H5 source (2026-09-18): MOSExIGLOO evaluated Gas_State_From_Prim on MOSE's never-written
# edge/corner ghosts -> f_Rtot 0/0 ("floating invalid" at CYCLE 1); GREEN after F5 (interior-only copy).
# PASS also needs no LK-interface: line (IGLOO's warning that a Langmuir-Knudsen fixed point was not bracketed).
# Every build tree links the same bin/hydra_MI2 (plain set() of CMAKE_RUNTIME_OUTPUT_DIRECTORY in src/app/CMakeLists.txt):
# the production binary is saved first and restored on every exit path. ~1.5 min.
set -u
H=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
T=${1:-${MI2_SCRATCH:-/tmp/mi2-checks-${USER:-$(id -un)}}/P7}
mkdir -p "$T"
cp -p "$H/bin/hydra_MI2" "$T/hydra_MI2.saved" || exit 2
restore() {
  cp "$T/hydra_MI2.saved" "$H/bin/hydra_MI2"
  cmp -s "$H/bin/hydra_MI2" "$T/hydra_MI2.saved" || echo "!! bin/hydra_MI2 NOT restored -- copy $T/hydra_MI2.saved back by hand"
}
trap restore EXIT

cmake -S "$H" -B "$T/build" -DEXE_NAME=MI2 \
  -DCMAKE_Fortran_COMPILER=${IFX:-ifx} \
  -DCMAKE_BUILD_TYPE=RELEASE -DUSE_OPENMP=true -DUSE_CANTERA=false -DUSE_SUNDIALS=false -DUSE_TECIO=false -DUSE_MPI=false \
  -DCMAKE_Fortran_FLAGS="-xHost -cpp -extend-source -fpe0 -traceback -g" > "$T/configure.log" 2>&1 || { echo "== p7 FAIL: configure ($T/configure.log)"; exit 2; }
cmake --build "$T/build" -j 8 > "$T/build.log" 2>&1 || { echo "== p7 FAIL: build ($T/build.log)"; exit 2; }
mv "$H/bin/hydra_MI2" "$T/hydra_MI2.fpe0"
cp "$T/hydra_MI2.saved" "$H/bin/hydra_MI2"

"$H/test/mi2-checks/run_case.sh" "$H/test/evap-box-ord2" "$T/run" 1 "$T/run.log" "$T/hydra_MI2.fpe0"
rc=$?
if [ $rc -eq 0 ] && grep -q "GAME OVER" "$T/run.log" && ! grep -q 'LK-interface:' "$T/run.log"; then
  echo "== p7 PASS: -fpe0 build of evap-box-ord2 reaches GAME OVER ($(grep '# wall' "$T/run.log"))"
  exit 0
fi
echo "== p7 FAIL: exit $rc"
grep -m1 'LK-interface:' "$T/run.log"
grep -n -A12 "forrtl" "$T/run.log" | head -20
exit 1
