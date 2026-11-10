#!/bin/bash

# ATLAS.sh - logic for running ATLAS

source "$SCRIPT_DIR/utils/common.sh"

function ATLAS_ops () {
  [[ $# == 0 ]] && print_usage
  if [[ $1 == import ]]; then
    import_from_ATLAS
    exit 1
  fi
  # BCB and ICB read a phase's tables from the case root when its phase.txt is there, else from fromATLAStoSolver/:
  # after GPB they must read GPB's, so the stored tables are copied only when GPB is not asked for
  local gpb=0 j
  for j in "$@"; do [[ $j == GPB ]] && gpb=1; done
  if [[ $gpb == 0 ]]; then
    cp INPUT/*phase.txt . 2>/dev/null || :
    cp INPUT/*thermo.dat . 2>/dev/null || :
    cp INPUT/*transport.dat . 2>/dev/null || :
    cp INPUT/*properties.dat . 2>/dev/null || :
  else
    # tables already in the case root (an earlier run, or written by hand) are kept: say that BCB and ICB read them
    local f stale=()
    for f in *phase.txt *thermo.dat *transport.dat *properties.dat; do
      if [[ -f $f ]]; then stale+=("$f"); fi
    done
    if (( ${#stale[@]} )); then
      log_warn "the case root holds ${stale[*]}: for a phase whose phase.txt is there, BCB and ICB read these, not the tables GPB writes now"
    fi
  fi
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
