#pragma once

#include <cstdint>

// Private trial ABI retained for legacy differential tests and RSFLIP01 replay.
extern "C" {
/**
 * @brief Runs a trial with the production or known-error ideal decoder.
 * @param ideal Zero selects RS correction; one selects byte-error peeling.
 * @param supplemental Three final-repair counters, published only on success.
 * @return Zero on success, nonzero on invalid input or a failed trial.
 * @details Older entry points remain compatible with saved records and
 * reference tests.
 */
int product_trial_decoder(uint64_t seed,
                          uint64_t batch,
                          uint64_t trial,
                          uint64_t k,
                          uint64_t passes,
                          int anchors,
                          int binary,
                          uint64_t* output,
                          int sampler,
                          uint32_t* positions,
                          int random,
                          uint8_t* residual,
                          uint64_t n1,
                          uint64_t k1,
                          uint64_t n2,
                          uint64_t k2,
                          int postprocessing,
                          int ideal,
                          uint64_t* supplemental);
// Supplemental counters are separate from the unchanged legacy 22 metrics.
int product_trial_postprocessing(uint64_t seed,
                                 uint64_t batch,
                                 uint64_t trial,
                                 uint64_t k,
                                 uint64_t passes,
                                 int anchors,
                                 int binary,
                                 uint64_t* output,
                                 int sampler,
                                 uint32_t* positions,
                                 int random,
                                 uint8_t* residual,
                                 uint64_t n1,
                                 uint64_t k1,
                                 uint64_t n2,
                                 uint64_t k2,
                                 int postprocessing,
                                 uint64_t* supplemental);
int product_trial_dimensions(uint64_t seed,
                             uint64_t batch,
                             uint64_t trial,
                             uint64_t k,
                             uint64_t passes,
                             int anchors,
                             int binary,
                             uint64_t* output,
                             int sampler,
                             uint32_t* positions,
                             int random,
                             uint8_t* residual,
                             uint64_t n1,
                             uint64_t k1,
                             uint64_t n2,
                             uint64_t k2);
int product_interrupt_install();
int product_interrupted();
uint64_t product_batch_k(uint64_t seed,
                         uint64_t batch,
                         uint64_t lo,
                         uint64_t hi);
int product_trial_reference(uint64_t seed,
                            uint64_t batch,
                            uint64_t trial,
                            uint64_t k,
                            uint64_t passes,
                            int anchors,
                            int binary,
                            uint64_t* output,
                            int sampler,
                            uint32_t* positions,
                            int random,
                            uint8_t* residual);
}
