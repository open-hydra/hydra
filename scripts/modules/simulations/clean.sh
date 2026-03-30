#!/bin/bash

# clean.sh - Cleans simulation folder

source "$SCRIPT_DIR/utils/common.sh"

log_info "Cleaning simulation folder"

rm -rf error* .ID
rm -rf OUTPUT bin

log_success "Clean completed"
