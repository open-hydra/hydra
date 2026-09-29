#!/usr/bin/env bash
set -euo pipefail

# building/hydra.sh - Handles building tasks for hydra solvers (Option A)

# -------------------------
# Defaults
# -------------------------
COMPILERS=""
USE_OPENMP="false"
USE_MPI="false"
USE_TECIO="false"
USE_SUNDIALS="false"
USE_CANTERA="false"
BUILD_TYPE="Release"
SUFFIX=""

# -------------------------
# Parse solver name + options
# -------------------------
SUFFIX="${1#hydra-}"
shift || true

while [[ $# -gt 0 ]]; do
    case "$1" in
        --compilers=*)
            [[ "$1" =~ ^--compilers=(intel|gnu)$ ]] \
                || { error "Invalid --compilers (intel|gnu)"; exit 1; }
            COMPILERS="${1#*=}"
            ;;
        --use-openmp)    USE_OPENMP="true" ;;
        --use-mpi)       USE_MPI="true" ;;
        --use-tecio)     USE_TECIO="true" ;;
        --use-sundials)  USE_SUNDIALS="true" ;;
        --use-cantera)   USE_CANTERA="true" ;;
        --type=*)        BUILD_TYPE="${1#*=}" ;;
        *)
            error "Unknown option '$1'"
            exit 1
            ;;
    esac
    shift
done

[[ -n "$SUFFIX" ]] || { error "Missing hydra suffix"; exit 1; }

# -------------------------
# Compiler selection
# -------------------------
if [[ $COMPILERS == "intel" ]]; then
    export FC="ifx" CC="icx" CXX="icpx"
elif [[ $COMPILERS == "gnu" ]]; then
    export FC="gfortran" CC="gcc" CXX="g++"
fi

# -------------------------
# Build directory
# -------------------------
cd "$HYDRADIR"

BUILD_DIR="$HYDRADIR/build/$SUFFIX"

log_info "Building hydra-$SUFFIX"
log_info "Build dir: $BUILD_DIR"

rm -rf "$BUILD_DIR"

git submodule update --init lib/ORION
git submodule update --init lib/OSLO
git submodule update --init lib/FLINT
git submodule update --init --recursive lib/third_party/FiNeR

git submodule update --init src/lib/MOSE
git submodule update --init src/lib/FUSS
git submodule update --init src/lib/GROOT
git submodule update --init src/lib/ARES

# hydra-MI2 only: its two condensed-phase solvers. ICE is a private repository, skipped by a
# recursive clone (update = none in .gitmodules): it is checked out here, on request.
if [[ "$SUFFIX" == "MI2" ]]; then
    git submodule update --init --checkout src/lib/ICE
    git submodule update --init src/lib/IGLOO
fi

# -------------------------
# Configure + build
# -------------------------
cmake -S "$HYDRADIR" -B "$BUILD_DIR" \
    -DEXE_NAME="$SUFFIX" \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DUSE_OPENMP="$USE_OPENMP" \
    -DUSE_MPI="$USE_MPI" \
    -DUSE_TECIO="$USE_TECIO" \
    -DUSE_SUNDIALS="$USE_SUNDIALS" \
    -DUSE_CANTERA="$USE_CANTERA"

cmake --build "$BUILD_DIR"

# -------------------------
# Write / update CMakePresets.json (SOURCE ROOT)
# -------------------------
# cached <build-dir> <variable> <fallback>: the value of a variable in the CMake cache of a build tree
cached() {
    local value
    value=$(grep -m1 "^$2:" "$1/CMakeCache.txt" | cut -d= -f2- || true)
    echo "${value:-$3}"
}

write_presets() {
    local presets_file="$HYDRADIR/CMakePresets.json"

    log_info "Regenerating CMakePresets.json"

    {
        echo "{"
        echo '  "version": 3,'
        echo '  "configurePresets": ['

        local first=true
        for build_dir in "$HYDRADIR"/build/*; do
            [[ -d "$build_dir" ]] || continue
            [[ -f "$build_dir/CMakeCache.txt" ]] || continue

            local solver
            solver="$(basename "$build_dir")"

            local FC CC CXX
            FC=$(grep '^CMAKE_Fortran_COMPILER:FILEPATH=' "$build_dir/CMakeCache.txt" | cut -d= -f2-)
            CC=$(grep '^CMAKE_C_COMPILER:FILEPATH=' "$build_dir/CMakeCache.txt" | cut -d= -f2-)
            CXX=$(grep '^CMAKE_CXX_COMPILER:FILEPATH=' "$build_dir/CMakeCache.txt" | cut -d= -f2-)

            # Each preset carries the options its own build tree was configured with
            local btype openmp mpi tecio sundials cantera
            btype=$(cached "$build_dir" CMAKE_BUILD_TYPE "$BUILD_TYPE")
            openmp=$(cached "$build_dir" USE_OPENMP "$USE_OPENMP")
            mpi=$(cached "$build_dir" USE_MPI "$USE_MPI")
            tecio=$(cached "$build_dir" USE_TECIO "$USE_TECIO")
            sundials=$(cached "$build_dir" USE_SUNDIALS "$USE_SUNDIALS")
            cantera=$(cached "$build_dir" USE_CANTERA "$USE_CANTERA")

            [[ "$first" == true ]] && first=false || echo ","

            cat <<EOF
    {
      "name": "$solver",
      "binaryDir": "\${sourceDir}/build/$solver",
      "cacheVariables": {
        "CMAKE_BUILD_TYPE": "$btype",
        "CMAKE_Fortran_COMPILER": "$FC",
        "CMAKE_C_COMPILER": "$CC",
        "CMAKE_CXX_COMPILER": "$CXX",
        "EXE_NAME": "$solver",
        "USE_OPENMP": "$openmp",
        "USE_MPI": "$mpi",
        "USE_TECIO": "$tecio",
        "USE_SUNDIALS": "$sundials",
        "USE_CANTERA": "$cantera"
      }
    }
EOF
        done

        echo
        echo "  ]"
        echo "}"
    } > "$presets_file"

    log_info "CMakePresets.json written successfully"
}

# Usage in build script
log_info "Writing CMakePresets.json"
write_presets
