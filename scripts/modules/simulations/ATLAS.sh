#!/bin/bash

# ATLAS.sh - logic for running ATLAS

source "$SCRIPT_DIR/utils/common.sh"

function ATLAS_ops () {
  [[ $# == 0 ]] && print_usage
  if [[ $1 == import ]]; then
    import_from_ATLAS
    exit 1
  fi
  cp INPUT/*phase.txt . 2>/dev/null || :
  cp INPUT/*thermo.dat . 2>/dev/null || :
  cp INPUT/*transport.dat . 2>/dev/null || :
  cp INPUT/*properties.dat . 2>/dev/null || :
  [[ ! -f "mesh.tec" ]] && cp MESH/mesh.tec . 2>/dev/null || :
  [[ ! -f "mesh.p3d" ]] && cp MESH/mesh.p3d . 2>/dev/null || :
  for j in $@; do
    [[ $j == ATLAS ]] && continue
    $ATLASDIR/ATLAS.sh $j
  done
  ask_user "type y to import"
  echo $reply
  if [[ $reply == "y" ]]; then 
    import_from_ATLAS
    log_success "imported"
  fi
  exit
}

function import_from_ATLAS {
  cp -r fromATLAStoSolver/* INPUT/. 2>/dev/null || :
  rm -rf fromATLAStoSolver mesh.tec mesh.p3d
  rm -f *phase.txt
  rm -f *thermo.dat
  rm -f *transport.dat
  rm -f *properties.dat
}

ATLAS_ops $@
