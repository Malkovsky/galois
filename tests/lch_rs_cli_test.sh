#!/usr/bin/env bash
set -euo pipefail

LCH_RS="${1:?usage: lch_rs_cli_test.sh <lch-rs>}"
WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cd "$WORKDIR"

expect_exit() {
  local want=$1
  shift
  set +e
  "$@"
  local got=$?
  set -e
  if [[ "$got" -ne "$want" ]]; then
    echo "expected exit $want, got $got: $*" >&2
    exit 1
  fi
}

assert_contains() {
  local path=$1
  local text=$2
  if ! grep -q -- "$text" "$path"; then
    echo "expected $path to contain: $text" >&2
    exit 1
  fi
}

assert_empty() {
  local path=$1
  if [[ -s "$path" ]]; then
    echo "expected $path to be empty" >&2
    exit 1
  fi
}

write_bytes() {
  local path=$1
  local size=$2
  python3 - "$path" "$size" <<'PY'
import sys
path, size = sys.argv[1], int(sys.argv[2])
data = bytes((i * 31 + 7) % 256 for i in range(size))
with open(path, "wb") as fh:
    fh.write(data)
PY
}

flip_payload() {
  local path=$1
  python3 - "$path" <<'PY'
import sys
from pathlib import Path
path = Path(sys.argv[1])
data = bytearray(path.read_bytes())
name_len = data[12] | (data[13] << 8)
offset = 64 + name_len
if offset >= len(data):
    raise SystemExit(f"no payload in {path}")
data[offset] ^= 0xFF
path.write_bytes(data)
PY
}

roundtrip() {
  local name=$1
  local size=$2
  local k=4
  local r=2
  local c=32
  write_bytes "$name" "$size"
  cp "$name" "$name.orig"
  expect_exit 0 "$LCH_RS" encode "$name" -k "$k" -r "$r" --chunk-size "$c" --jobs 1
  expect_exit 0 "$LCH_RS" verify --jobs 1 "$name".lch.6.*
  expect_exit 0 "$LCH_RS" decode -o restored --jobs 1 "$name".lch.6.*
  cmp -s "$name.orig" restored
  rm -f "$name".lch.6.* restored "$name" "$name.orig"
}

roundtrip empty 0
roundtrip one 1
roundtrip almost_chunk 31
roundtrip one_stripe 128
roundtrip one_stripe_plus 129

small_name=small_padding
write_bytes "$small_name" 1205
cp "$small_name" "$small_name.orig"
expect_exit 0 "$LCH_RS" encode "$small_name" -k 3 -r 2 --jobs 1
expected_share_size=$((64 + ${#small_name} + 402 + 8))
for share in small_padding.lch.5.*; do
  actual_share_size=$(stat -c %s "$share")
  if [[ "$actual_share_size" -ne "$expected_share_size" ]]; then
    echo "expected $share to use $expected_share_size bytes, got $actual_share_size" >&2
    exit 1
  fi
done
expect_exit 0 "$LCH_RS" decode -o small_padding.restored --jobs 1 small_padding.lch.5.*
cmp -s small_padding.orig small_padding.restored
rm -f small_padding.lch.5.* small_padding.restored small_padding small_padding.orig

write_bytes fewmeg 2097152
cp fewmeg fewmeg.orig
expect_exit 0 "$LCH_RS" encode fewmeg -k 4 -r 2 --chunk-size 4096 --jobs 2
expect_exit 0 "$LCH_RS" verify --jobs 2 fewmeg.lch.6.*
expect_exit 0 "$LCH_RS" decode -o fewmeg.restored --jobs 2 fewmeg.lch.6.*
cmp -s fewmeg.orig fewmeg.restored
rm -f fewmeg.lch.6.* fewmeg.restored fewmeg fewmeg.orig

write_bytes damaged 200
cp damaged damaged.orig
expect_exit 0 "$LCH_RS" encode damaged -k 4 -r 2 --chunk-size 32 --jobs 1
flip_payload damaged.lch.6.0
expect_exit 1 "$LCH_RS" verify --jobs 1 damaged.lch.6.*
expect_exit 0 "$LCH_RS" decode -o damaged.restored --jobs 1 damaged.lch.6.*
cmp -s damaged.orig damaged.restored
rm -f damaged.lch.6.* damaged.restored damaged damaged.orig

write_bytes truncated 128
cp truncated truncated.orig
expect_exit 0 "$LCH_RS" encode truncated -k 4 -r 2 --chunk-size 32 --jobs 1
truncate -s 80 truncated.lch.6.0
expect_exit 0 "$LCH_RS" decode -o truncated.restored --jobs 1 truncated.lch.6.*
cmp -s truncated.orig truncated.restored
rm -f truncated.lch.6.* truncated.restored truncated truncated.orig

write_bytes missing 200
cp missing missing.orig
expect_exit 0 "$LCH_RS" encode missing -k 2 -r 2 --chunk-size 32 --jobs 1
rm missing missing.lch.4.0 missing.lch.4.1
expect_exit 0 "$LCH_RS" decode -o missing.restored --jobs 1 missing.lch.4.*
cmp -s missing.orig missing.restored
rm -f missing.lch.4.* missing.restored missing.orig

write_bytes too_damaged 256
expect_exit 0 "$LCH_RS" encode too_damaged -k 4 -r 2 --chunk-size 32 --jobs 1
flip_payload too_damaged.lch.6.0
flip_payload too_damaged.lch.6.1
flip_payload too_damaged.lch.6.2
expect_exit 2 "$LCH_RS" decode -o too_damaged.restored --jobs 1 too_damaged.lch.6.*
if [[ -e too_damaged.restored ]]; then
  echo "corrupt destination was kept" >&2
  exit 1
fi
rm -f too_damaged.lch.6.* too_damaged

write_bytes exists_dest 64
expect_exit 0 "$LCH_RS" encode exists_dest -k 4 -r 2 --chunk-size 32 --jobs 1
echo already > exists_dest.out
expect_exit 3 "$LCH_RS" decode -o exists_dest.out exists_dest.lch.6.*
rm -f exists_dest.lch.6.* exists_dest exists_dest.out

expect_exit 3 "$LCH_RS" encode missing_kr.bin
write_bytes invalid_kr 32
expect_exit 3 "$LCH_RS" encode invalid_kr -k 200 -r 200 --chunk-size 32
rm -f invalid_kr

write_bytes conflict 128
cp conflict conflict.orig
write_bytes other2 128
expect_exit 0 "$LCH_RS" encode conflict -k 2 -r 2 --chunk-size 32 --jobs 1
expect_exit 0 "$LCH_RS" encode other2 -k 4 -r 1 --chunk-size 32 --jobs 1
rm conflict.lch.4.0
expect_exit 0 "$LCH_RS" decode -o conflict.restored --jobs 1 conflict.lch.4.* other2.lch.5.*
cmp -s conflict.orig conflict.restored
rm -f conflict.lch.4.* other2.lch.5.* other2 conflict conflict.orig conflict.restored

write_bytes progress 4096
expect_exit 0 "$LCH_RS" encode progress -k 4 -r 2 --chunk-size 32 \
  --jobs 1 --progress=auto 2>progress.auto.log
assert_empty progress.auto.log
expect_exit 0 "$LCH_RS" verify --jobs 1 --progress=always progress.lch.6.* \
  >progress.verify.out 2>progress.verify.log
assert_contains progress.verify.log "Verifying"
assert_contains progress.verify.log "%"
expect_exit 0 "$LCH_RS" decode -o progress.restored --jobs 1 \
  --progress=always progress.lch.6.* 2>progress.decode.log
assert_contains progress.decode.log "Decoding"
assert_contains progress.decode.log "%"
cmp -s progress progress.restored
expect_exit 0 "$LCH_RS" encode progress -k 4 -r 2 --chunk-size 32 \
  --jobs 1 --force --progress=never -v 2>progress.verbose.log
assert_contains progress.verbose.log "encode k=4 r=2"
assert_contains progress.verbose.log "backend="
expect_exit 3 "$LCH_RS" verify --progress=sometimes progress.lch.6.* \
  >/dev/null 2>/dev/null
rm -f progress.lch.6.* progress progress.restored progress.*.log progress.verify.out

write_bytes piped 200
cp piped piped.orig
expect_exit 0 "$LCH_RS" encode -k 4 -r 2 --chunk-size 32 --name piped \
  --jobs 1 --progress=always < piped 2>piped.progress.log
assert_contains piped.progress.log "Encoding"
rm piped
expect_exit 0 "$LCH_RS" decode -o - --jobs 1 --progress=never \
  piped.lch.6.* > piped.restored
cmp -s piped.orig piped.restored
rm -f piped.lch.6.* piped.restored piped.orig piped.progress.log

if command -v gzip >/dev/null 2>&1; then
  write_bytes gzip_pipe 65537
  cp gzip_pipe gzip_pipe.orig
  gzip -c gzip_pipe | "$LCH_RS" encode -k 4 -r 2 --chunk-size 1024 \
    --name gzip_pipe.gz --jobs 2 --progress=never
  rm gzip_pipe
  "$LCH_RS" decode -o - --jobs 1 --progress=never gzip_pipe.gz.lch.6.* |
    gzip -dc >gzip_pipe.restored
  cmp -s gzip_pipe.orig gzip_pipe.restored
  rm -f gzip_pipe.gz.lch.6.* gzip_pipe.orig gzip_pipe.restored
fi
