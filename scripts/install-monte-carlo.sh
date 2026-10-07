#!/usr/bin/env bash
set -euo pipefail

if [[ ${1:-} == --help || ${1:-} == -h ]]; then
  printf 'Usage: %s [PREFIX]\nBuild native Release rs-product-monte-carlo and rs-product-test and install both (default: ~/.local). Tests remain disabled.\nMonte Carlo runtime dimensions: --n1 256 --k1 224 --n2 256 --k2 254. Strong N,R remain aligned powers of two; weak R=2 shortening (e.g. --n2 175 --k2 173) and full R=4 (--n2 256 --k2 252) are supported. Small codes need explicit flip bounds <=8*n1*n2.\nUse rs-product-test --help for codeword generation and corruption fixtures.\nRequires CMake, a C++20 compiler, and Git/network for pinned nlohmann/json v3.12.0 headers. No Python runtime or sudo.\n' "$0"
  printf 'Optional --postprocessing (default off) enables one final validated repair stage. Counters report repaired stall events, detected strong miscorrection columns (or an existence lower bound of one), and corrected strong columns. Enabled saved-flip runs use RSFLIP02; replay verifies supplemental counters as well as the legacy 22 metrics. See rs-product-monte-carlo --help.\n'
  exit 0
fi
if (( $# > 1 )) || [[ ${1:-} == -* ]]; then
  printf 'Usage: %s [PREFIX]\n' "$0" >&2
  exit 2
fi

prefix=${1:-"$HOME/.local"}
if [[ -z $prefix ]]; then
  printf 'Installation prefix must not be empty.\n' >&2
  exit 2
fi
# Resolve relative prefixes before changing to the repository root.
[[ $prefix == /* ]] || prefix="$PWD/$prefix"
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd -- "$root"
cmake --preset experimental
cmake --build --preset experimental --parallel
cmake --install "$root/build/experimental-preset" \
  --prefix "$prefix" --component product-monte-carlo
cmake --install "$root/build/experimental-preset" \
  --prefix "$prefix" --component product-test
printf 'Installed rs-product-monte-carlo and rs-product-test under %s. Ensure its bin directory is on PATH.\n' "$prefix"
