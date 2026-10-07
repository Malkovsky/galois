#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>

namespace mc {

/** @brief Compute SHA-256 over bytes for private Monte Carlo file identities.
 * @param bytes Input bytes, including any embedded zero bytes.
 * @return The 32 digest bytes in standard big-endian order.
 */
inline std::array<uint8_t, 32> Sha256(std::span<const uint8_t> bytes) {
  constexpr std::array<uint32_t, 64> k = {
      0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1,
      0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
      0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786,
      0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
      0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147,
      0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
      0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
      0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
      0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a,
      0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
      0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
  std::array<uint32_t, 8> state = {0x6a09e667, 0xbb67ae85, 0x3c6ef372,
                                   0xa54ff53a, 0x510e527f, 0x9b05688c,
                                   0x1f83d9ab, 0x5be0cd19};
  const auto compress = [&](const uint8_t* block) {
    std::array<uint32_t, 64> w{};
    for (size_t i = 0; i < 16; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        w[i] = (w[i] << 8) | block[4 * i + j];
      }
    }
    for (size_t i = 16; i < 64; ++i) {
      const auto x = w[i - 15], y = w[i - 2];
      w[i] = w[i - 16] + (std::rotr(x, 7) ^ std::rotr(x, 18) ^ (x >> 3)) +
             w[i - 7] + (std::rotr(y, 17) ^ std::rotr(y, 19) ^ (y >> 10));
    }
    auto [a, b, c, d, e, f, g, h] = state;
    for (size_t i = 0; i < 64; ++i) {
      const uint32_t t1 =
          h + (std::rotr(e, 6) ^ std::rotr(e, 11) ^ std::rotr(e, 25)) +
          ((e & f) ^ (~e & g)) + k[i] + w[i];
      const uint32_t t2 =
          (std::rotr(a, 2) ^ std::rotr(a, 13) ^ std::rotr(a, 22)) +
          ((a & b) ^ (a & c) ^ (b & c));
      h = g;
      g = f;
      f = e;
      e = d + t1;
      d = c;
      c = b;
      b = a;
      a = t1 + t2;
    }
    const std::array<uint32_t, 8> working = {a, b, c, d, e, f, g, h};
    for (size_t i = 0; i < state.size(); ++i) {
      state[i] += working[i];
    }
  };
  // Convert before multiplying so the bit count is modulo 2^64, not size_t.
  const uint64_t bits = static_cast<uint64_t>(bytes.size()) * uint64_t{8};
  while (bytes.size() >= 64) {
    compress(bytes.data());
    bytes = bytes.subspan(64);
  }
  std::array<uint8_t, 128> tail{};
  std::copy(bytes.begin(), bytes.end(), tail.begin());
  tail[bytes.size()] = 0x80;
  const size_t padded = bytes.size() < 56 ? 64 : 128;
  for (size_t i = 0; i < 8; ++i) {
    tail[padded - 1 - i] = static_cast<uint8_t>(bits >> (8 * i));
  }
  compress(tail.data());
  if (padded == 128) {
    compress(tail.data() + 64);
  }
  std::array<uint8_t, 32> digest{};
  for (size_t i = 0; i < digest.size(); ++i) {
    digest[i] = static_cast<uint8_t>(state[i / 4] >> (8 * (3 - i % 4)));
  }
  return digest;
}

}  // namespace mc
