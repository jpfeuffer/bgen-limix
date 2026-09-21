#ifndef BGEN_VARIANT_H
#define BGEN_VARIANT_H

#include "bgen/export.h"
#include <inttypes.h>

struct bgen_file;

/** Variant metadata.
 * @struct bgen_variant
 */
struct bgen_variant
{
    uint64_t                   genotype_offset; /**< Genotype offset (bgen file). */
    struct bgen_string const*  id;              /**< Variant identification. */
    struct bgen_string const*  rsid;            /**< RSID. */
    struct bgen_string const*  chrom;           /**< Chromossome name. */
    uint32_t                   position;        /**< Base-pair position. */
    uint16_t                   nalleles;        /**< Number of alleles. */
    struct bgen_string const** allele_ids;      /**< Allele ids. */
};

/** Seek to and parse the first variant.
 *
 * @param bgen_file Bgen file handler.
 * @param error Set to non-zero on failure; left at `0` on success, including
 *   the end-of-file case (`NULL` return with no data left to read).
 * @return Variant metadata, owned by the caller (@ref bgen_variant_destroy).
 *   Return `NULL` at end of file or on failure; check `error` to tell them apart.
 */
BGEN_EXPORT struct bgen_variant* bgen_variant_begin(struct bgen_file* bgen_file, int* error);
/** Parse the variant that follows the stream's current position.
 *
 * @param bgen_file Bgen file handler.
 * @param error Set to non-zero on failure; left at `0` on success, including
 *   the end-of-file case (`NULL` return with no data left to read).
 * @return Variant metadata, owned by the caller (@ref bgen_variant_destroy).
 *   Return `NULL` at end of file or on failure; check `error` to tell them apart.
 */
BGEN_EXPORT struct bgen_variant* bgen_variant_next(struct bgen_file* bgen_file, int* error);
/** Release a variant returned by @ref bgen_variant_begin or @ref bgen_variant_next.
 *
 * @param variant Variant metadata.
 */
BGEN_EXPORT void bgen_variant_destroy(struct bgen_variant const* variant);

#endif
