/** Open, close, and query bgen file.
 * @file bgen/file.h
 */
#ifndef BGEN_FILE_H
#define BGEN_FILE_H

#include "bgen/export.h"
#include "bgen/s3stream.h"
#include <stdbool.h>
#include <stdint.h>

/** Bgen file handler.
 * @struct bgen_file
 */
struct bgen_file;

/** Return the handler's underlying stream, for callers that read variant
 * records themselves (e.g. resolving a `.bgi` offset the same way sequential
 * scanning does).
 *
 * @param bgen_file Bgen file handler.
 * @return Stream handle. Owned by `bgen_file`; do not close it directly.
 */
BGEN_EXPORT stream_handle* bgen_file_stream(struct bgen_file const* bgen_file);

/** Open bgen file and return a handler.
 *
 * Remember to call @ref bgen_file_close to close the file and release
 * resources after the interaction has finished.
 *
 * @param filepath File path to the bgen file.
 * @return Bgen file handler. Return `NULL` on failure.
 */
BGEN_EXPORT struct bgen_file* bgen_file_open(char const* filepath);
/** Close bgen file handler.
 *
 * @param bgen_file Bgen file handler.
 */
BGEN_EXPORT void bgen_file_close(struct bgen_file const* bgen_file);
/** Get the file's layout version.
 *
 * @param bgen_file Bgen file handler.
 * @return `1` for Layout 1 or `2` for Layout 2.
 */
BGEN_EXPORT unsigned bgen_file_layout(struct bgen_file const* bgen_file);
/** Get the file's genotype block compression method.
 *
 * @param bgen_file Bgen file handler.
 * @return `0` for no compression, `1` for zlib, `2` for zstd.
 */
BGEN_EXPORT unsigned bgen_file_compression(struct bgen_file const* bgen_file);
/** Get the number of samples.
 *
 * @param bgen_file Bgen file handler.
 * @return Number of samples.
 */
BGEN_EXPORT uint32_t bgen_file_nsamples(struct bgen_file const* bgen_file);
/** Get the number of variants.
 *
 * @param bgen_file Bgen file handler.
 * @return Number of variants.
 */
BGEN_EXPORT uint32_t bgen_file_nvariants(struct bgen_file const* bgen_file);
/** Check if the file contain sample identifications.
 *
 * @param bgen_file Bgen file handler.
 * @return `true` if bgen file contains the sample ids; `false` otherwise.
 */
BGEN_EXPORT bool bgen_file_contain_samples(struct bgen_file const* bgen_file);
/** Return all sample identifications.
 *
 * @param bgen_file Bgen file handler.
 * @return Sample identifications. Return `NULL` on failure.
 */
BGEN_EXPORT struct bgen_samples* bgen_file_read_samples(struct bgen_file* bgen_file);
/** Open a variant for genotype queries.
 *
 * @param bgen_file Bgen file handler.
 * @param genotype_offset Genotype offset obtained from @ref bgen_variant.genotype_offset.
 * @return Variant genotype handler. Return `NULL` on failure.
 */
BGEN_EXPORT struct bgen_genotype* bgen_file_open_genotype(struct bgen_file* bgen_file,
                                                          uint64_t          genotype_offset);

/** Read probabilities for multiple variants in one call.
 *
 * Reads @p n_offsets variants into the caller-allocated array @p out, which must
 * point to @p n_offsets × @ref bgen_file_nsamples × @p ncombs contiguous doubles.
 * All variants must have the same @p ncombs; if any variant differs, the function
 * returns 1. Pass the value from @ref bgen_genotype_ncombs of a representative
 * variant (or detect it with one @ref bgen_file_open_genotype call first).
 * Supports the same layouts and compression methods as single-record reads.
 *
 * If @p offsets are sorted in ascending order, the function will skip redundant
 * seeks for consecutive variants where the file pointer is already positioned.
 *
 * This function is NOT thread-safe when sharing a single @p bgen_file handle.
 * For multithreaded use, open one @p bgen_file per thread.
 *
 * @param bgen_file   Bgen file handler.
 * @param offsets     Array of genotype offsets (from metafile).
 * @param n_offsets   Number of offsets.
 * @param ncombs      Expected genotype combinations per sample.
 * @param out         Pre-allocated output array [n_offsets × nsamples × ncombs].
 * @return 0 on success; 1 on error.
 */
BGEN_EXPORT int bgen_file_read_genotypes_batch(struct bgen_file* bgen_file,
                                               uint64_t const*   offsets,
                                               uint32_t          n_offsets,
                                               uint32_t          ncombs,
                                               double*           out);

/**
 * Read the number of genotype combinations for multiple variants (header-only,
 * no probability decode).
 * Compressed layout-2 blocks are still decompressed.
 *
 * @param ncombs_out  Output array of length n_offsets; filled with gt.ncombs.
 */
BGEN_EXPORT int bgen_file_read_ncombs_batch(struct bgen_file* bgen_file,
                                            uint64_t const*   offsets,
                                            uint32_t          n_offsets,
                                            uint32_t*         ncombs_out);

/**
 * Read probabilities for multiple variants into a padded rectangular buffer.
 *
 * Output shape: (n_offsets, nsamples, max_ncombs) — row-major.
 * Unused combination slots (when gt.ncombs < max_ncombs) are NaN-filled.
 * Returns an error if any variant's ncombs > max_ncombs.
 */
BGEN_EXPORT int bgen_file_read_genotypes_batch_padded(struct bgen_file* bgen_file,
                                                      uint64_t const*   offsets,
                                                      uint32_t          n_offsets,
                                                      uint32_t          max_ncombs,
                                                      double*           out);

/** Caller-allocated batch metadata. Each pointer is optional (NULL to omit).
 *
 * Variant arrays have n_offsets elements; sample arrays have
 * n_offsets * bgen_file_nsamples elements, in variant-major order.
 * probability_count is the structural probability width of a sample, even
 * when missing: ploidy * nalleles if phased, otherwise the number of unordered
 * genotypes. ncombs is the maximum-ploidy width used by single-record reads.
 * Missing probabilities and unused padding are both NaN; use missing and
 * probability_count rather than inspecting NaNs to distinguish them.
 */
struct bgen_batch_metadata
{
    bool* phased;                 /**< Variant array. */
    uint8_t* ploidy;               /**< Sample array. */
    bool* missing;                /**< Sample array. */
    uint32_t* ncombs;              /**< Variant array. */
    uint32_t* probability_count;   /**< Sample array. */
    uint16_t* nalleles;            /**< Variant array. */
    uint8_t* nbits;                /**< Variant array. */
};

/** Read genotype metadata without decoding or allocating probability arrays.
 * Supports every layout/compression accepted by single-record reads.
 * Compressed layout 2 still requires decompression of the genotype block.
 * Layout 1 reads the triplets to determine original missingness.
 * Validates framing, headers and payload lengths, not probability sums.
 * Ordering and duplicate offsets are preserved. An empty request succeeds.
 * Returns 0 on success, 1 on failure; outputs may be partially filled on error.
 * Caller must provide enough storage for every non-NULL metadata pointer.
 * These functions are not thread-safe on a shared file handle.
 */
BGEN_EXPORT int bgen_file_read_genotype_metadata(
    struct bgen_file* file, uint64_t const* offsets, uint32_t n_offsets,
    struct bgen_batch_metadata* metadata);

/** Decode a structured padded batch using the single-record decoding core.
 * out has n_offsets * nsamples * max_ncombs elements (variant-major).
 * Both extra variant-width and sample-width slots are NaN. max_ncombs may be
 * zero for zero-width records. Returns 1 if any variant's width exceeds it.
 * Explicit width avoids a pre-scan: one block read/decompression per offset.
 * NULL metadata omits metadata; a NULL out is allowed only for zero elements.
 * All requested dimensions and byte counts are checked for size overflow.
 * Other lifetime, ordering, error and thread-safety rules match metadata reads.
 */
BGEN_EXPORT int bgen_file_read_genotypes_batch64(
    struct bgen_file* file, uint64_t const* offsets, uint32_t n_offsets,
    uint32_t max_ncombs, double* out, struct bgen_batch_metadata* metadata);

/** Float32 counterpart of bgen_file_read_genotypes_batch64. */
BGEN_EXPORT int bgen_file_read_genotypes_batch32(
    struct bgen_file* file, uint64_t const* offsets, uint32_t n_offsets,
    uint32_t max_ncombs, float* out, struct bgen_batch_metadata* metadata);

/** Probability-model assumption for lightweight dimension discovery.
 * AUTO reads phase; other modes are assumptions, not phase validation.
 * MIXED reserves both GP and HP widths for every selected variant.
 */
enum bgen_probability_model {
    BGEN_MODEL_AUTO = 0,
    BGEN_MODEL_UNPHASED = 1,
    BGEN_MODEL_PHASED = 2,
    BGEN_MODEL_MIXED = 3
};

/** Maxima for discovery, without per-sample output allocation.
 * gp_width is zero if no unphased model is selected. hp_width is at least
 * one when phased is selected, even at zero ploidy, to retain phase presence.
 * decode_width = max(gp_width, hp_width). Empty requests return all zeros.
 */
struct bgen_dimensions {
    uint32_t gp_width;
    uint32_t hp_width;
    uint32_t decode_width;
};

/** Discover dimensions using actual sample ploidy or declared max ploidy.
 * actual_ploidy=true scans sample ploidies; false uses the declared header.
 * Layout 1 uses fixed diploid/biallelic dimensions without payload access.
 * Raw/zlib layout 2 read only necessary prefixes; declared+explicit model
 * need only the first eight uncompressed bytes. AUTO must also reach phase
 * after the sample-ploidy bytes. Zstd uses full validated block loading.
 * Explicit models also apply to layout 1 (PHASED reserves width 4, MIXED both).
 * This is NOT full integrity validation: skipped data/phase assumptions must
 * be validated during later full reads. In particular, trailing truncation
 * may go undetected by prefix strategies. Width overflow is an error.
 * Outputs are written only on success. Caller owns output; no arrays returned.
 * File-handle thread-safety rules match the batch APIs.
 */
BGEN_EXPORT int bgen_file_read_dimensions(
    struct bgen_file* file, uint64_t const* offsets, uint32_t n_offsets,
    bool actual_ploidy, enum bgen_probability_model model,
    struct bgen_dimensions* out);

/** Dimension discovery with a caller's ploidy capacity assumption.
 * max_ploidy in [0,63] overrides the ploidy used to calculate capacities,
 * without validating sample ploidies against it. Pass -1 for no override.
 * Raw/zlib AUTO reads phase but need not scan sample ploidy; explicit models
 * need only the header. Zstd still loads/validates a full compact block.
 * Later full decode must validate any user-supplied ploidy bound.
 */
BGEN_EXPORT int bgen_file_read_dimensions_with_ploidy(
    struct bgen_file* file, uint64_t const* offsets, uint32_t n_offsets,
    bool actual_ploidy, enum bgen_probability_model model, int max_ploidy,
    struct bgen_dimensions* out);

/** Single-pass compact-capacity batch, without a width pre-scan.
 * Capacity bounds actual sample structural widths, not declared max ploidy.
 * Full headers/payload lengths are validated before each variant is written.
 * Missing samples count towards bounds. Effective ncombs in metadata is the
 * actual maximum sample width (zero at zero samples), unlike default batches.
 * No declared-width probability buffer is allocated. Earlier variants may be
 * filled on error; a failing variant is not written. Capacity is never a
 * truncation limit. All other semantics match batch64/32.
 */
BGEN_EXPORT int bgen_file_read_genotypes_batch_compact64(
    struct bgen_file* file, uint64_t const* offsets, uint32_t n_offsets,
    uint32_t capacity, double* out, struct bgen_batch_metadata* metadata);
BGEN_EXPORT int bgen_file_read_genotypes_batch_compact32(
    struct bgen_file* file, uint64_t const* offsets, uint32_t n_offsets,
    uint32_t capacity, float* out, struct bgen_batch_metadata* metadata);

#endif
