#!/bin/bash
# kill.sh — Kill background solver jobs safely

source "$SCRIPT_DIR/utils/common.sh"
log_info "Killing simulation"

if [[ ! -f ".ID" ]]; then
    log_error "No running solver process found"
    exit 1
fi

read PID < .ID && kill $PID
if [[ -f ".nosleep/KID" ]]; then
    read PID < .nosleep/KID && kill $PID
fi
rm -rf .ID .nosleep

log_success "Simulation killed"