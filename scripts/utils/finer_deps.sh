#!/bin/bash
# finer_deps.sh [hydra-dir] - makes lib/third_party/FiNeR buildable: FiNeR no longer tracks its five dependencies,
# yet its CMakeLists adds them. Fills them from the last FiNeR commit that pinned them (v2.0.6, cfc9194), then
# returns to the pin. Sourced by building/hydra.sh and building/submods.sh; run it alone before any other build
# that is handed hydra's FiNeR (e.g. ICE's install.sh --include-finer=...).

fill_finer_deps() {
    local root="${1:-$HYDRADIR}"
    local finer="$root/lib/third_party/FiNeR"
    [[ -f "$finer/CMakeLists.txt" ]] || git -C "$root" submodule update --init --recursive lib/third_party/FiNeR
    if [[ ! -f "$finer/src/third_party/PENF/CMakeLists.txt" ]]; then
        ( cd "$finer" && pin=$(git rev-parse HEAD) \
          && git checkout -q cfc9194 \
          && git submodule update --init --recursive \
          && git checkout -q "$pin" )
    fi
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    set -euo pipefail
    fill_finer_deps "${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
fi
