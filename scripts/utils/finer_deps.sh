#!/bin/bash
# finer_deps.sh [hydra-dir] - makes lib/third_party/FiNeR buildable: FiNeR tracks none of its five dependencies, yet
# its CMakeLists adds them. Puts each one at the pinned commit below (cloned, or moved there if it is a checkout at
# another commit). Sourced by building/hydra.sh and building/submods.sh; run it alone before any other build that is
# handed hydra's FiNeR (e.g. ICE's install.sh --include-finer=...). The five releases MOSE, ARES and FUSS pin for their
# own FiNeR v2.3.0, and the BACKUP set of ATLAS and ICE.

FINER_DEPS=(
    "PENF       4cb9d326c08f109ed5126d83c27e8e70c98b7eaa"   # v2.0.1
    "StringiFor b12f62bab32dde9fed99b697ef6eff4c120d77ec"   # v1.8.0
    "FACE       da7b7407b51e56257a09970099997098e3a52d99"   # v1.1.17
    "FLAP       a3741e61f360addf4772833d078619cb8b72e0e1"   # v2.5.3
    "BeFoR64    c90ce9507a24e166ebaeed53d67bf26512bb9263"   # v1.2.4
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
