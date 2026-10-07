#!/usr/bin/env bash
set -euo pipefail

if [[ ${1:-} == --help || ${1:-} == -h ]]; then
  printf 'Usage: %s [PREFIX]\nBuild native Release lch-rs and install it (default: ~/.local).\nRequires CMake and a C++20 compiler; dependency downloads may require network access. No sudo is invoked.\n' "$0"
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
cmake --preset cli -DGF256_BUILD_MONTE_CARLO=OFF
cmake --build --preset cli --parallel
cmake --install "$root/build/cli-preset" \
  --prefix "$prefix" --component lch-rs
printf 'Installed lch-rs under %s. Ensure its bin directory is on PATH.\n' "$prefix"
