/** Create and query a metafile.
 * @file bgen/metafile.h
 *
 * A bgen metafile (v05) is defined as follows:
 *
 * [ char[13] : signature "bgen index 05"       ],  \
 * [ uint32_t : number of variants              ],   | Header block
 * [ uint32_t : number of partitions            ],   |
 * [ uint64_t : metadata block size             ],   |
 * [ uint8_t  : all_biallelic flag (0=no,1=yes) ],   /
 * [                                                \
 *   [                                              |
 *     uint64_t : partition offset (this file)      | Offsets block
 *   ], ...                                         |
 * ]                                                /
 * [                                                \
 *   uint64_t        : genotype offset (bgen file)  |
 *   uint16_t, str,  : variant id                   |
 *   uint16_t, str,  : variant rsid                 |
 *   uint16_t, str,  : variant chrom                |
 *   uint32_t,       : genetic position             | Metadata block
 *   uint16_t,       : number of alleles            |
 *   [                                              |
 *     uint32_t, str : allele id                    |
 *   ], ...                                         |
 * ], ...                                           /
 *
 * Version 05 adds the all_biallelic byte at offset 29.
 * Version 04 files are read with all_biallelic == 2 (unknown).
 */
#ifndef BGEN_METAFILE_H_PRIVATE
#define BGEN_METAFILE_H_PRIVATE

#include <inttypes.h>
#include <stdio.h>

#define BGEN_METAFILE_SIGNATURE_V04 "bgen index 04"
#define BGEN_METAFILE_SIGNATURE_V05 "bgen index 05"
#define BGEN_METAFILE_SIGNATURE     BGEN_METAFILE_SIGNATURE_V05  /* written by default */
#define BGEN_METAFILE_HEADER_SIZE   (13 + 4 + 4 + 8 + 1)        /* 30 bytes (v05) */
#define BGEN_METAFILE_HEADER_SIZE_V04 (13 + 4 + 4 + 8)          /* 29 bytes (v04) */

struct bgen_metafile
{
    char*     filepath;
    FILE*     stream;
    uint32_t  nvariants;
    uint32_t  npartitions;
    uint64_t  metadata_block_size;
    uint64_t* partition_offset; /**< Array of partition offsets */
    /** 1 = all variants are biallelic, 0 = at least one multiallelic, 2 = unknown (v04) */
    uint8_t   all_biallelic;
};

uint32_t bgen_metafile_partition_size(uint32_t nvariants, uint32_t npartitions);

#endif
