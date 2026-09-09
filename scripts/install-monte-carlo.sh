#!/usr/bin/env bash
set -euo pipefail

if [[ ${1:-} == --help || ${1:-} == -h ]]; then
  printf 'Usage: %s [PREFIX]\nBuild native Release rs-product-monte-carlo and install it (default: ~/.local).\nRequires CMake, a C++20 compiler, OpenSSL 3 development files, and Git/network for pinned jsoncons headers. No Python runtime or sudo.\n' "$0"
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
printf 'Installed rs-product-monte-carlo under %s. Ensure its bin directory is on PATH.\n' "$prefix"
