#!/bin/bash

# deps.sh - the components of hydra at their releases (cmake/Releases.cmake, deps.lock)
#   figaro deps [--update] <NAME>...   check each component out at its release of deps.lock, or the newest in its
#                                      range; --update selects the newest releases in the ranges again
# Without a name: every component but ICE, which is private (name it to check it out).

source "$SCRIPT_DIR/utils/common.sh"
cd "$HYDRADIR"

update=false
names=()
for a in "$@"; do
  case "$a" in
    -u|--update) update=true ;;
    -*) log_error "Unknown option for deps: $a"; exit 1 ;;
    *) names+=("$a") ;;
  esac
done
if [[ ${#names[@]} -eq 0 ]]; then
  names=($(sed -n 's/^ *"\([A-Z][A-Za-z0-9]*\) .*/\1/p' cmake/Releases.cmake | grep -vx ICE))
fi

args=(-DHYDRA_DEPS="$(IFS=';'; echo "${names[*]}")")
[[ $update == true ]] && args+=(-DSEMVER_UPDATE="$(IFS=';'; echo "${names[*]}")")
log_info "Checking out ${names[*]} at their releases"
cmake -S cmake/deps -B build/.deps "${args[@]}" > build/.deps.log 2>&1 \
  || { tail -20 build/.deps.log; log_error "deps failed (see build/.deps.log)"; exit 1; }
grep -E '^-- +[A-Za-z0-9]+ v[0-9][^,]*, ' build/.deps.log | sed 's/^-- */  /'
log_success "deps done"
