/* Layout-2 genotype payload encoding. SPDX-License-Identifier: MIT */
#ifndef BGEN_ENCODE_H
#define BGEN_ENCODE_H

#include "bgen/export.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum bgen_encode_type {
    BGEN_ENCODE_FLOAT32 = 1,
    BGEN_ENCODE_FLOAT64 = 2,
    BGEN_ENCODE_INT8 = 3,
    BGEN_ENCODE_INT16 = 4,
    BGEN_ENCODE_INT32 = 5,
    BGEN_ENCODE_INT64 = 6
};

/** One variant. All buffers are contiguous, native-endian, and borrowed.
 * data is a row-major [nsamples, width] array of the indicated type.
 * Float input requires ploidy and missing arrays of length nsamples.
 * Unphased active width is C(P+K-1,P); phased active width is P*K.
 * Only active entries of nonmissing samples are read; padding is ignored.
 * Integer input is hard calls: -1 missing allele, -2 trailing ploidy fill.
 * Its ploidy and missing arrays are optional. Explicit ploidy must agree
 * with trailing fill. Without missing, any -1 marks the whole call missing.
 * A nonmissing call may not contain -1. Even missing hard calls must have
 * valid allele indexes and trailing fill. Byte lengths bound every read.
 * Inputs must not overlap outputs or be mutated during either API call.
 */
struct bgen_layout2_variant {
    void const *data;
    size_t data_bytes;
    size_t width;
    uint8_t const *ploidy;
    size_t ploidy_bytes;
    uint8_t const *missing;
    size_t missing_bytes;
    uint32_t nsamples;
    uint16_t nalleles;
    uint8_t phased;
    uint8_t bits;
    enum bgen_encode_type type;
};

/** Errors are per-call, not global. SIZE_MAX denotes no variant/sample. */
struct bgen_encode_error {
    size_t variant;
    size_t sample;
    char message[192];
};

/** Measure payload bytes, validating metadata, buffer bounds and hard calls.
 * Probability numeric validation occurs in encode. Each block must fit a
 * uint32 BGEN length field; total size must not exceed max_output_bytes.
 * Returns 0 on success, -1 on error (including allocation failure).
 */
BGEN_EXPORT int bgen_layout2_size(
    struct bgen_layout2_variant const *variants, size_t nvariants,
    size_t max_output_bytes, size_t *size, struct bgen_encode_error *error);

/** Encode concatenated uncompressed layout-2 payloads (no length prefixes,
 * variant identifiers, compression, or file header). offsets has capacity
 * at least nvariants+1, starts at zero and ends at the total byte count.
 * Capacity is also the allocation limit. On error output may be partially
 * written and must not be used. No allocations escape this function.
 * Float simplexes must be finite, in [0,1], sum to 1 within absolute 1e-6.
 * Normalized largest-remainder quantization breaks ties by original index.
 * Missing calls have zero stored values. The last state is implicit.
 */
BGEN_EXPORT int bgen_layout2_encode(
    struct bgen_layout2_variant const *variants, size_t nvariants,
    uint8_t *output, size_t capacity, uint64_t *offsets, size_t offsets_capacity,
    struct bgen_encode_error *error);

/** Specialized compatibility kernels return these error codes (and context).
 * Other failures return -1. Output contents are unspecified after failure.
 */
enum bgen_encode_status {
    BGEN_ENCODE_INVALID_PLOIDY = -2,
    BGEN_ENCODE_INVALID_ALLELE = -3,
    BGEN_ENCODE_MIXED_PLOIDY = -4,
    BGEN_ENCODE_STRING_LENGTH = -5,
    BGEN_ENCODE_BUFFER_TOO_SMALL = -6
};

/** Row-major int8 [nvariants, nsamples, 2], bool-byte [nvariants].
 * Specialized biallelic, 8-bit legacy policy: P=1/2 only, slot 0 -2 is
 * rejected, and zero samples report Pmin/Pmax=2. Half-missing is missing.
 * Output rows have stride >= 10+3*nsamples; lengths gives actual row sizes.
 * Explicit byte counts bound inputs and outputs. No hidden allocations.
 */
BGEN_EXPORT int bgen_layout2_hardcalls_legacy(
    int8_t const *genotypes, size_t genotypes_bytes,
    uint8_t const *phased, size_t phased_bytes, size_t nvariants, uint32_t nsamples,
    uint8_t *output, size_t output_bytes, size_t row_stride,
    uint32_t *lengths, size_t lengths_capacity, struct bgen_encode_error *error);

/** Padded byte strings: row ends at first NUL or width, without decoding.
 * Three identifier strings have uint16 lengths, alleles uint32 lengths.
 */
struct bgen_encode_strings {
    uint8_t const *data;
    size_t bytes;
    size_t width;
};

struct bgen_layout2_fixed_batch {
    struct bgen_encode_strings varid, rsid, chrom, allele1, allele2;
    int32_t const *position;
    size_t position_bytes;
    int8_t const *genotypes;
    size_t genotypes_bytes;
    uint8_t const *phased;
    size_t phased_bytes;
    size_t nvariants;
    uint32_t nsamples;
    uint8_t uniform_ploidy;
    size_t total_string_length;
};

/** Exact fixed full-variant row width, with RFC1950 stored-block compression.
 * Ploidy must be 1/2. Returns -1 on arithmetic/BGEN uint32 framing overflow.
 */
BGEN_EXPORT int bgen_layout2_fixed_size(
    uint32_t nsamples, uint8_t uniform_ploidy, size_t total_string_length,
    size_t *size, struct bgen_encode_error *error);

/** Write concatenated fixed-stride complete variant records (no file header).
 * Includes identifying fields, C/D lengths and deterministic zlib stored
 * blocks of at most 65535 bytes each. Same hard-call policy as legacy above.
 * Actual sum of the five string lengths must equal total_string_length in
 * every variant. One worst-case genotype scratch buffer is reused. Inputs
 * and output must not alias. Negative positions are rejected.
 */
BGEN_EXPORT int bgen_layout2_fixed_encode(
    struct bgen_layout2_fixed_batch const *batch, uint8_t *output,
    size_t output_bytes, struct bgen_encode_error *error);

#ifdef __cplusplus
}
#endif
#endif
