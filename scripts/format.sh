#!/usr/bin/env bash
set -euo pipefail

case "${1:---check}" in
  --check) flags=(--dry-run --Werror) ;;
  --write) flags=(-i) ;;
  *) printf 'Usage: bash scripts/format.sh [--check|--write]\n' >&2; exit 2 ;;
esac
if [[ $# -gt 1 ]]; then
  printf 'Expected at most one argument\n' >&2
  exit 2
fi

root="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
formatter="${CLANG_FORMAT:-clang-format-18}"
version="$("$formatter" --version)"
if [[ ! "$version" =~ (^|[[:space:]])version[[:space:]]18\.1\.3([[:space:]]|$) ]]; then
  printf 'Expected clang-format 18.1.3, got: %s\n' "$version" >&2
  exit 2
fi

# Git's index defines ownership: never traverse submodules or untracked assets.
mapfile -d '' -t tracked < <(git -C "$root" ls-files -z)
files=()
for file in "${tracked[@]}"; do
  case "$file" in
    third_party/*|agentic/*|assets/*|reference/*|lessons/*) continue ;;
  esac
  case "$file" in
    *.c|*.cc|*.cpp|*.cxx|*.h|*.hh|*.hpp|*.hxx|*.inc|*.inl)
      [[ ! -f "$root/$file" || -L "$root/$file" ]] || files+=("$root/$file") ;;
  esac
done
if [[ ${#files[@]} -eq 0 ]]; then
  printf 'No tracked owned C/C++ files found\n' >&2
  exit 2
fi
printf '%s: %s (%d tracked owned C/C++ files)\n' "${1:---check}" "$version" "${#files[@]}"
"$formatter" --style="file:$root/.clang-format" "${flags[@]}" "${files[@]}"
