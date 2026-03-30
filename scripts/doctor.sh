#!/usr/bin/env bash
set -euo pipefail

# source logging functions
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/utils/common.sh"

require_cmd() {
    command -v "$1" >/dev/null 2>&1
}

require_file() {
    [[ -f "$1" ]]
}

require_dir() {
    [[ -d "$1" ]]
}

check_version() {
    local have="$1"
    local need="$2"
    [[ "$(printf '%s\n' "$need" "$have" | sort -V | head -n1)" == "$need" ]]
}

compiler_supports_openmp() {
    local cc="$1"
    echo 'int main(){return 0;}' > /tmp/omp_test.c
    if "$cc" -fopenmp /tmp/omp_test.c -o /tmp/omp_test 2>/dev/null; then
        rm -f /tmp/omp_test.c /tmp/omp_test
        return 0
    fi
    rm -f /tmp/omp_test.c
    return 1
}


# -------------------------
# Header
# -------------------------
echo
log_info "Running sanity checks"

# -------------------------
# Environment checks
# -------------------------
echo
echo "Environment checks"
if [[ -z "${HYDRADIR:-}" ]]; then
    log_warn "HYDRADIR is not set"
    echo "      Fix: run ./bootstrap.sh or export HYDRADIR=/path/to/hydra"
else
    log_success "HYDRADIR is set"
    require_dir "$HYDRADIR" || log_error "HYDRADIR does not exist"
fi

# -------------------------
# Project structure
# -------------------------
echo
echo "Project structure checks"
require_file "$HYDRADIR/CMakeLists.txt" \
    && log_success "CMakeLists.txt found" \
    || log_error "Missing CMakeLists.txt"

require_dir "$HYDRADIR/build" \
    && log_success "build/ directory exists" \
    || log_warn "build/ directory missing (will be created on build)"

# -------------------------
# Tools
# -------------------------
echo
echo "Toolchain checks"
if require_cmd cmake; then
    cmake_ver=$(cmake --version | head -n1 | awk '{print $3}')
    if check_version "$cmake_ver" "3.23"; then
        log_success "cmake $cmake_ver"
    else
        log_error "cmake $cmake_ver (>= 3.23 required)"
        exit 1
    fi
else
    log_error "cmake not found"
    exit 1
fi

require_cmd git && log_success "git found" || log_error "git not found" || exit 1
require_cmd bash && log_success "bash found" || log_error "bash not found" || exit 1

if [[ -n "${CONDA_EXE:-}" ]]; then
    log_success "Conda installed"
else
    log_error "Conda not found"
fi


# -------------------------
# Compilers
# -------------------------
echo
echo "Compiler checks"

if require_cmd gfortran; then
    log_success "GNU Fortran (gfortran)"
else
    log_warn "gfortran not found"
fi

if require_cmd gcc && require_cmd g++; then
    log_success "GNU C/C++"
else
    log_warn "GNU C/C++ incomplete"
fi

if require_cmd ifx || require_cmd ifort; then
    log_success "Intel Fortran"
else
    log_warn "Intel Fortran not found"
fi

# -------------------------
# Optional features
# -------------------------
echo
echo "Optional features"
# MPI
require_cmd mpicc && log_success "MPI available" || log_warn "MPI not found"
# OpenMP
if require_cmd gcc; then
    if compiler_supports_openmp gcc; then
        log_success "GNU OpenMP supported"
    else
        log_error "GNU OpenMP not supported"
    fi
else
    log_warn "gcc not found (OpenMP not checked)"
fi
if require_cmd icx; then
    compiler_supports_openmp icx && log_success "Intel OpenMP supported"
elif require_cmd icc; then
    compiler_supports_openmp icc && log_success "Intel OpenMP supported"
else
    log_warn "Intel compilers not found"
fi
# macOS note
if [[ "$(uname)" == "Darwin" ]]; then
    if ! ls /usr/local/lib/libomp.* >/dev/null 2>&1 && \
       ! ls /opt/homebrew/lib/libomp.* >/dev/null 2>&1; then
        log_warn "libomp not found (macOS)"
        echo "      Fix: brew install libomp"
    fi
fi

# -------------------------
# Presets sanity
# -------------------------
echo
echo "CMake presets"

if require_file "$HYDRADIR/CMakePresets.json"; then
    if cmake --list-presets >/dev/null 2>&1; then
        log_success "CMakePresets.json valid"
    else
        log_warn "CMakePresets.json exists but is invalid"
        echo "      Fix: run figaro build hydra-<solver>"
    fi
else
    log_info "No CMakePresets.json (will be generated on build)"
fi

# -------------------------
# End
# -------------------------
echo
echo "Doctor finished"
