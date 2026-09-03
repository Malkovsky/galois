#pragma once

// Reproduce GF-Complete's Autotools feature definitions from the exact ISA
// macros active for each directly compiled translation unit.
#if defined(__SSE2__) && !defined(INTEL_SSE2)
#define INTEL_SSE2
#endif

#if defined(__SSE3__) && !defined(INTEL_SSE3)
#define INTEL_SSE3
#endif

#if defined(__SSSE3__) && !defined(INTEL_SSSE3)
#define INTEL_SSSE3
#endif

#if (defined(__SSE4_1__) || defined(__SSE4_2__)) && !defined(INTEL_SSE4)
#define INTEL_SSE4
#endif

#if defined(__PCLMUL__) && defined(INTEL_SSE4) && !defined(INTEL_SSE4_PCLMUL)
#define INTEL_SSE4_PCLMUL
#endif
