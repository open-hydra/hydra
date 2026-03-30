#!/usr/bin/env bash
set -euo pipefail

# source logging functions
source "./scripts/utils/common.sh"

echo "Bootstrapping hydra"

# -------------------------
# Sanity checks
# -------------------------
source "./scripts/doctor.sh"
echo "If you see warnings, hope they will be fixed while building!"
echo

# -------------------------
# Create HYDRADIR and build directory
# -------------------------
HYDRADIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$HYDRADIR/build"

# -------------------------
# Shell integration
# -------------------------
log_info "Setting up shell integration"
rm -f $HYDRADIR/.setvars.sh
echo 'export HYDRADIR='$HYDRADIR >> .setvars.sh
echo 'source '$HYDRADIR'/scripts/figaro_completion.sh' >> .setvars.sh
if [[ $SHELL == *"zsh"* ]]; then
    echo 'figaro () { '$HYDRADIR'/scripts/figaro.sh $@; }' >> .setvars.sh
    RCFILE=$HOME/.zshrc
    log_success "Detected zsh shell"
elif [[ $SHELL == *"bash"* ]]; then
    echo 'function figaro () { '$HYDRADIR'/scripts/figaro.sh $@; }' >> .setvars.sh
    RCFILE=$HOME/.bashrc
    log_success "Detected bash shell"
else
    log_error "Unsupported shell, skipping shell integration"
    exit 1
fi
echo 'export -f figaro' >> .setvars.sh
grep -v "hydra" $RCFILE > tmpfile && mv tmpfile $RCFILE
echo 'source '$HYDRADIR'/.setvars.sh' >> $RCFILE

# -------------------------
# End
# -------------------------
echo
log_success "Bootstrap complete"
echo "Restart your shell or run: source $RCFILE"
