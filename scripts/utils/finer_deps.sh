#!/bin/bash
# finer_deps.sh [hydra-dir] - makes lib/third_party/FiNeR buildable: FiNeR tracks none of its five dependencies, yet
# its CMakeLists adds them. Puts each one at the pinned commit below (cloned, or moved there if it is a checkout at
# another commit). Sourced by building/hydra.sh and building/submods.sh; run it alone before any other build that is
# handed hydra's FiNeR (e.g. ICE's install.sh --include-finer=...). The same five commits are pinned by IGLOO, ICE and
# ATLAS for their own FiNeR.

FINER_DEPS=(
    "PENF       34e10af852f81822bac0df7c5f73936e0655f54f"
    "StringiFor 77b75b1d7f7012d2984135b231008f751ab32dda"
    "FACE       95226cd6135b9daeb8ff2c748d2065deee230358"
    "FLAP       7f5ec8e1f26be5e2a263fdd2118501b340c0d461"
    "BeFoR64    40d19c5f6127088126606b90448e99e02d8b641c"
)

fill_finer_deps() {
    local root="${1:-$HYDRADIR}"
    local finer="$root/lib/third_party/FiNeR"
    local entry dep sha dir url
    [[ -f "$finer/CMakeLists.txt" ]] || git -C "$root" submodule update --init lib/third_party/FiNeR || return 1
    for entry in "${FINER_DEPS[@]}"; do
        read -r dep sha <<< "$entry"
        dir="$finer/src/third_party/$dep"
        url="https://github.com/szaghi/$dep"
        if [[ ! -e "$dir" ]] || [[ -d "$dir" && -z "$(ls -A "$dir")" ]]; then
            git clone -q "$url" "$dir" || return 1
        elif ! git -C "$dir" rev-parse --git-dir > /dev/null 2>&1; then
            echo "finer_deps: $dir exists but is not a git checkout; remove it and rerun" >&2; return 1
        fi
        if [[ "$(git -C "$dir" rev-parse HEAD 2>/dev/null)" != "$sha" ]]; then
            git -C "$dir" cat-file -e "$sha^{commit}" 2>/dev/null || git -C "$dir" fetch -q "$url" "$sha" || return 1
            git -C "$dir" checkout -q "$sha" || return 1
        fi
    done
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    set -euo pipefail
    fill_finer_deps "${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
fi
