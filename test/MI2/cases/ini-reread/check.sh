#!/bin/bash
# check.sh <run-dir> <exe>
# The hot re-read of input.ini ([HYDRA-IO] ini-diter): igloo-evap re-reading its file at every second coupling step.
#   the run ends without an [ERROR]: the keys of the solvers MI2 does not run (FUSS, the multigrid levels) are not
#   required of an MI2 input;
#   re-reading an unchanged file changes nothing: OUTPUT byte-identical to the run without the re-read.
source "$(dirname "${BASH_SOURCE[0]}")/../../checklib.sh"

a=$(variant once "$cases/igloo-evap"); runin "$a"
b=$(variant reread "$cases/igloo-evap"); sed -i 's/^ini-diter *=.*/ini-diter = 2/' "$b/input.ini"
grep -q '^ini-diter = 2$' "$b/input.ini" || ko "the base has no ini-diter line to edit"
runin "$b"; rc=$?
expect "the run with the re-read ends (exit $rc)" ended "$b/run.log"
expect "OUTPUT byte-identical to the run without the re-read" diff -rq "$a/OUTPUT" "$b/OUTPUT"
finish
