#ifndef BGEN_FILE_H_PRIVATE
#define BGEN_FILE_H_PRIVATE

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

struct bgen_file;

FILE*       bgen_file_stream(struct bgen_file const* bgen_file);
char const* bgen_file_filepath(struct bgen_file const* bgen_file);
unsigned    bgen_file_layout(struct bgen_file const* bgen_file);
unsigned    bgen_file_compression(struct bgen_file const* bgen_file);
int         bgen_file_seek_variants_start(struct bgen_file* bgen_file);

/**
 * Read and decompress the compressed genotype block at the current stream position
 * into the file's reusable scratch buffer (Plan A + B).
 *
 * Reads:  4 bytes compressed_length, 4 bytes uncompressed_length, then the data.
 * Result: *chunk points to file->chunk_scratch (owned by bgen_file, do NOT free).
 *         *size   is the uncompressed byte count.
 * Returns 0 on success, 1 on error.
 */
int bgen_file_decompress_block(struct bgen_file* bgen_file, char** chunk, size_t* size);

/**
 * Return pointer to the file's ploidy scratch buffer (nsamples bytes, Plan A).
 * Valid until the next call to bgen_file_decompress_block or bgen_file_close.
 */
uint8_t* bgen_file_ploidy_scratch(struct bgen_file* bgen_file);

#endif
