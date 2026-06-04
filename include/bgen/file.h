/** Open, close, and query bgen file.
 * @file bgen/file.h
 */
#ifndef BGEN_FILE_H
#define BGEN_FILE_H

#include "bgen/export.h"
#include <stdbool.h>
#include <stdint.h>

/** Bgen file handler.
 * @struct bgen_file
 */
struct bgen_file;

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

/** Read probabilities for multiple variants in one call (Plan C: batch API).
 *
 * Reads @p n_offsets variants into the caller-allocated array @p out, which must
 * point to @p n_offsets × @ref bgen_file_nsamples × @p ncombs contiguous doubles.
 * All variants must have the same @p ncombs; if any variant differs, the function
 * returns 1. Pass the value from @ref bgen_genotype_ncombs of a representative
 * variant (or detect it with one @ref bgen_file_open_genotype call first).
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

#endif
