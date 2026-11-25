#!/bin/bash
# check.sh <run-dir> <exe>
# FPE cleanliness of the coupling paths: igloo-evap with gas-order = 2 (IGLOO samples the gas bilinearly, through
# MOSE's ghost cells) run by a hydra_MI2 that traps floating-point exceptions (-fpe0 -traceback). PASS: the run ends,
# and no LK-interface: line (IGLOO's warning that a Langmuir-Knudsen fixed point was not bracketed). A trap prints the
# traceback.
# The trapping binary: FPE_EXE, else built here in <run-dir>/build (RELEASE flags with -fpe0 -traceback -g, ifx or
# IFX; about 1.5 min). Every build tree links bin/hydra_MI2 (CMAKE_RUNTIME_OUTPUT_DIRECTORY), so the binary there is
# saved first and restored on every exit path. <exe> is not used.
source "$(dirname "${BASH_SOURCE[0]}")/../../checklib.sh"
H=$(cd "$cases/../../.." && pwd)

if [ -n "${FPE_EXE:-}" ]; then exe=$FPE_EXE
else
  saved=$run/hydra_MI2.saved
  if [ -e "$H/bin/hydra_MI2" ]; then cp -p "$H/bin/hydra_MI2" "$saved"; fi
  restore() {
    if [ -e "$saved" ]; then cp -p "$saved" "$H/bin/hydra_MI2"
      cmp -s "$H/bin/hydra_MI2" "$saved" || echo "!! bin/hydra_MI2 NOT restored: copy $saved back by hand"
    else rm -f "$H/bin/hydra_MI2"; fi
  }
  trap restore EXIT
  cmake -S "$H" -B "$run/build" -DEXE_NAME=MI2 -DCMAKE_Fortran_COMPILER="${IFX:-ifx}" -DCMAKE_BUILD_TYPE=RELEASE \
    -DUSE_OPENMP=true -DUSE_CANTERA=false -DUSE_SUNDIALS=false -DUSE_TECIO=false -DUSE_MPI=false \
    -DCMAKE_Fortran_FLAGS="-xHost -cpp -extend-source -fpe0 -traceback -g" > "$run/configure.log" 2>&1 \
    || { ko "configure ($run/configure.log)"; finish; }
  cmake --build "$run/build" -j 8 > "$run/build.log" 2>&1 || { ko "build ($run/build.log)"; finish; }
  mv "$H/bin/hydra_MI2" "$run/hydra_MI2.fpe0"; restore; exe=$run/hydra_MI2.fpe0
  ok "the -fpe0 build"
fi

d=$(variant ord2 "$cases/igloo-evap"); sed -i 's/^gas-order   = 1$/gas-order   = 2/' "$d/input.ini"
expect "the deck has gas-order = 2" grep -q '^gas-order   = 2$' "$d/input.ini"
runin "$d"; rc=$?
expect "the trapping binary runs igloo-evap at gas-order 2 to its end (exit $rc)" ended "$d/run.log"
grep -n -A12 forrtl "$d/run.log" | head -20
expect "no LK-interface: line" eval "! grep -q 'LK-interface:' '$d/run.log'"
finish
