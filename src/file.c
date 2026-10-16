#include "file.h"
#include "bgen/file.h"
#include "bgen/genotype.h"
#include "bgen/s3stream.h"
#include "bstring.h"
#include "bmath.h"
#include "free.h"
#include "genotype.h"
#include "layout1.h"
#include "layout2.h"
#include "mem.h"
#include "report.h"
#include "samples.h"
#include "strdup.h"
#include "unzlib.h"
#include "unzstd.h"
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

struct bgen_file
{
    char*             filepath;
    stream_handle*  stream;
    uint32_t          nvariants;
    uint32_t nsamples;
    unsigned compression;
    unsigned layout;
    bool     contain_sample;
    int64_t  samples_start;
    int64_t  variants_start;
    /* Plan A: reusable decompression scratch buffers */
    char*    chunk_scratch;       /**< reusable decompressed chunk buffer */
    size_t   chunk_scratch_size;  /**< current allocation size of chunk_scratch */
    uint8_t* ploidy_scratch;      /**< nsamples bytes, pre-allocated at open */
    /* Plan E: reusable probability scratch for padded batch (scatter with NaN-padding) */
    double*  prob_scratch;        /**< temp probability buffer for scatter-with-padding */
    size_t   prob_scratch_size;   /**< current allocation in doubles */
    /* Plan B: persistent decompressor contexts */
    bgen_zlib_ctx* zlib_ctx;
    bgen_zstd_ctx* zstd_ctx;
};

static struct bgen_file* bgen_file_create(char const* filepath);
static int               bgen_file_read_header(struct bgen_file* bgen);

struct bgen_file* bgen_file_open(char const* filepath)
{
    struct bgen_file* bgen = bgen_file_create(filepath);
    if (bgen == NULL)
        return NULL;

    bgen->variants_start = 0;
    if (bgen_stream_fread(bgen->stream, &bgen->variants_start, 4) != 1) {
        bgen_perror_eof(stream_handle_eof(bgen->stream), "could not read the `variants_start` field");
        goto err;
    }

    bgen->variants_start += 4;

    if (bgen_file_read_header(bgen)) {
        bgen_error("could not read bgen header");
        goto err;
    }

    /* Pre-allocate ploidy scratch once nsamples is known (Plan A). */
    if (bgen->nsamples > 0) {
        bgen->ploidy_scratch = malloc(bgen->nsamples * sizeof(uint8_t));
        if (!bgen->ploidy_scratch) {
            bgen_error("could not allocate ploidy scratch buffer");
            goto err;
        }
    }

    /* if they actually exist */
    if ((bgen->samples_start = stream_handle_tell(bgen->stream)) < 0) {
        bgen_perror("could not ftell");
        goto err;
    }

    return bgen;

err:
    bgen_file_close(bgen);
    return NULL;
}

void bgen_file_close(struct bgen_file const* bgen)
{
    if (bgen->stream != NULL)
        stream_handle_close(bgen->stream);
    bgen_free(bgen->filepath);
    bgen_free(bgen->chunk_scratch);
    bgen_free(bgen->ploidy_scratch);
    bgen_free(bgen->prob_scratch);
    if (bgen->zlib_ctx) bgen_zlib_ctx_destroy(bgen->zlib_ctx);
    if (bgen->zstd_ctx) bgen_zstd_ctx_destroy(bgen->zstd_ctx);
    bgen_free(bgen);
}

uint32_t bgen_file_nsamples(struct bgen_file const* bgen) { return bgen->nsamples; }

uint32_t bgen_file_nvariants(struct bgen_file const* bgen) { return bgen->nvariants; }

bool bgen_file_contain_samples(struct bgen_file const* bgen) { return bgen->contain_sample; }

struct bgen_samples* bgen_file_read_samples(struct bgen_file* bgen)
{
    char* block = NULL;

    if (stream_handle_seek(bgen->stream, bgen->samples_start, SEEK_SET)) {
        bgen_perror("could not fseek to `samples_start`");
        return NULL;
    }

    if (!bgen->contain_sample) {
        bgen_warning("file does not contain sample ids");
        return NULL;
    }

    struct bgen_samples* samples = bgen_samples_create(bgen->nsamples);

    uint32_t block_size = 0;
    if (bgen_stream_fread(bgen->stream, &block_size, sizeof(block_size)) != 1) {
        bgen_perror_eof(stream_handle_eof(bgen->stream), "could not read block size");
        goto err;
    }

    block = malloc(block_size - sizeof(block_size));
    if (bgen_stream_fread(bgen->stream, block, block_size - sizeof(block_size)) != 1) {
        bgen_perror_eof(stream_handle_eof(bgen->stream), "could not read samples block");
        goto err;
    }

    char const* block_ptr = block;
    uint32_t    nsamples = 0;
    bgen_memfread(&nsamples, &block_ptr, sizeof(nsamples));

    if (nsamples != bgen->nsamples) {
        bgen_error("number of samples mismatch (corrupted file?)");
        goto err;
    }

    for (uint32_t i = 0; i < bgen->nsamples; ++i) {
        struct bgen_string const* sample_id = bgen_string_memfread(&block_ptr, 2);
        bgen_samples_set(samples, i, sample_id);
    }

    if ((bgen->variants_start = stream_handle_tell(bgen->stream)) < 0) {
        bgen_error("could not ftell `variants_start`");
        goto err;
    }

    bgen_free(block);
    return samples;

err:
    bgen_samples_destroy(samples);
    bgen_free(block);
    return NULL;
}

static int open_genotype(struct bgen_file* bgen, uint64_t genotype_offset,
                         struct bgen_genotype* genotype)
{
    genotype->layout = bgen->layout;
    genotype->offset = genotype_offset;

    if (genotype_offset > INT64_MAX || genotype_offset < (uint64_t)bgen->variants_start) {
        bgen_error("invalid genotype offset");
        goto err;
    }

    if (stream_handle_tell(bgen->stream) != (int64_t)genotype_offset &&
        stream_handle_seek(bgen_file_stream(bgen), (int64_t)genotype_offset, SEEK_SET)) {
        bgen_perror("could not fseek a variant");
        goto err;
    }

    if (bgen_file_layout(bgen) == 1) {
        if (bgen_layout1_read_header(bgen, genotype))
            goto err;
    } else if (bgen_file_layout(bgen) == 2) {
        if (bgen_layout2_read_header(bgen, genotype))
            goto err;
    } else {
        bgen_error("unrecognized layout type %d", bgen_file_layout(bgen));
        goto err;
    }

    return 0;
err:
    return 1;
}

struct bgen_genotype* bgen_file_open_genotype(struct bgen_file* bgen, uint64_t genotype_offset)
{
    if (!bgen) {
        bgen_error("file handle is NULL");
        return NULL;
    }
    struct bgen_genotype* genotype = bgen_genotype_create();
    if (!genotype) {
        bgen_error("could not allocate genotype");
        return NULL;
    }
    if (!open_genotype(bgen, genotype_offset, genotype))
        return genotype;
    bgen_genotype_close(genotype);
    return NULL;
}

stream_handle* bgen_file_stream(struct bgen_file const* bgen_file) { return bgen_file->stream; }

/* See file.h: an fread(buf, size, 1, stream)-compatible wrapper, so call
 * sites' `!= 1` / `< 1` checks stay correct without inspecting byte counts. */
size_t bgen_stream_fread(stream_handle* stream, void* buf, size_t size)
{
    return (size_t)(stream_handle_read(stream, buf, size) == (int64_t)size);
}

char const* bgen_file_filepath(struct bgen_file const* bgen_file)
{
    return bgen_file->filepath;
}

unsigned bgen_file_layout(struct bgen_file const* bgen_file) { return bgen_file->layout; }

unsigned bgen_file_compression(struct bgen_file const* bgen_file)
{
    return bgen_file->compression;
}

int bgen_file_seek_variants_start(struct bgen_file* bgen_file)
{
    if (stream_handle_seek(bgen_file->stream, bgen_file->variants_start, SEEK_SET)) {
        bgen_perror("could not jump to variants start");
        return 1;
    }
    return 0;
}

int bgen_file_decompress_block(struct bgen_file* f, char** chunk, size_t* size)
{
    char* compressed = NULL;

    uint32_t compressed_length = 0;
    if (bgen_stream_fread(f->stream, &compressed_length, 4) < 1) {
        bgen_perror_eof(stream_handle_eof(f->stream), "could not read compressed length");
        goto err;
    }
    if (compressed_length < 4) {
        bgen_error("compressed length too small (corrupted file?)");
        goto err;
    }
    compressed_length -= 4;

    uint32_t uncompressed_length = 0;
    if (bgen_stream_fread(f->stream, &uncompressed_length, 4) < 1) {
        bgen_perror_eof(stream_handle_eof(f->stream), "could not read uncompressed length");
        goto err;
    }

    /* Grow scratch buffer if needed (Plan A: reuse instead of malloc each call). */
    if (uncompressed_length > f->chunk_scratch_size) {
        char* p = realloc(f->chunk_scratch, uncompressed_length);
        if (!p) { bgen_error("could not realloc chunk scratch"); goto err; }
        f->chunk_scratch      = p;
        f->chunk_scratch_size = uncompressed_length;
    }

    compressed = malloc(compressed_length);
    if (!compressed) { bgen_error("could not malloc compressed chunk"); goto err; }
    if (bgen_stream_fread(f->stream, compressed, compressed_length) < 1) {
        bgen_perror_eof(stream_handle_eof(f->stream), "could not read compressed data");
        goto err;
    }

    size_t out_size = uncompressed_length;

    if (f->compression == 1) {
        /* Lazily create context (Plan B). */
        if (!f->zlib_ctx) f->zlib_ctx = bgen_zlib_ctx_create();
        if (bgen_unzlib_reuse(f->zlib_ctx, compressed, compressed_length,
                              &f->chunk_scratch, &out_size))
            goto err;
    } else if (f->compression == 2) {
        if (!f->zstd_ctx) f->zstd_ctx = bgen_zstd_ctx_create();
        void* scratch = f->chunk_scratch;
        if (bgen_unzstd_reuse(f->zstd_ctx, compressed, compressed_length, &scratch, &out_size))
            goto err;
    } else {
        bgen_error("unrecognized compression method %u", f->compression);
        goto err;
    }

    bgen_free(compressed);
    *chunk = f->chunk_scratch;
    *size  = uncompressed_length;
    return 0;

err:
    bgen_free(compressed);
    return 1;
}

uint8_t* bgen_file_ploidy_scratch(struct bgen_file* f)
{
    return f->ploidy_scratch;
}

static int copy_metadata(struct bgen_genotype const* gt, uint32_t index,
                          struct bgen_batch_metadata* metadata)
{
    if (!metadata)
        return 0;
    if (metadata->phased) metadata->phased[index] = gt->phased != 0;
    if (metadata->ncombs) metadata->ncombs[index] = gt->ncombs;
    if (metadata->nalleles) metadata->nalleles[index] = gt->nalleles;
    if (metadata->nbits) metadata->nbits[index] = bgen_genotype_nbits(gt);
    if (!metadata->ploidy && !metadata->missing && !metadata->probability_count)
        return 0;
    uint32_t counts[64];
    if (metadata->probability_count) {
        for (unsigned p = gt->min_ploidy; p <= gt->max_ploidy; ++p) {
            if (gt->layout == 1)
                counts[p] = 3;
            else if (gt->phased)
                counts[p] = p * gt->nalleles;
            else {
                unsigned count;
                if (bgen_checked_ncombs(p, gt->nalleles, &count))
                    return 1;
                counts[p] = count;
            }
        }
    }
    size_t start = (size_t)index * gt->nsamples;
    for (uint32_t sample = 0; sample < gt->nsamples; ++sample) {
        uint8_t value = gt->ploidy_missingness[sample];
        uint8_t ploidy = value & 127;
        if (metadata->ploidy) metadata->ploidy[start + sample] = ploidy;
        if (metadata->missing) metadata->missing[start + sample] = value >> 7;
        if (metadata->probability_count) metadata->probability_count[start + sample] = counts[ploidy];
    }
    return 0;
}

static int read_batch(struct bgen_file* f, uint64_t const* offsets, uint32_t n_offsets,
                      uint32_t width, void* out, unsigned precision, bool exact,
                      struct bgen_batch_metadata* metadata)
{
    if (!f || (n_offsets && !offsets)) {
        bgen_error("invalid batch file handle or offsets");
        return 1;
    }
    size_t item_size = precision == 32 ? sizeof(float) : sizeof(double);
    uint64_t rows = (uint64_t)n_offsets * f->nsamples;
    if (rows > SIZE_MAX / sizeof(uint32_t) ||
        (precision && width && rows > SIZE_MAX / item_size / width)) {
        bgen_error("batch dimensions overflow");
        return 1;
    }
    if (precision && rows && width && !out) {
        bgen_error("batch probability output is NULL");
        return 1;
    }
    bool padding_filled = false;
    for (uint32_t i = 0; i < n_offsets; ++i) {
        struct bgen_genotype gt;
        bgen_genotype_init(&gt);
        if (open_genotype(f, offsets[i], &gt))
            goto err;
        if (precision && (gt.ncombs > width || (exact && gt.ncombs != width))) {
            bgen_error("batch max_ncomb=%u does not accommodate ncombs=%u at offset %" PRIu64,
                       width, gt.ncombs, offsets[i]);
            goto err;
        }
        if (copy_metadata(&gt, i, metadata))
            goto err;
        if (precision && f->nsamples && width) {
            size_t slice_start = (size_t)i * f->nsamples * width;
            char* slice = (char*)out + slice_start * item_size;
            void* destination = slice;
            if (gt.ncombs != width && gt.ncombs) {
                size_t bytes = (size_t)f->nsamples * gt.ncombs * item_size;
                size_t needed = bytes / sizeof(double) + (bytes % sizeof(double) != 0);
                if (needed > SIZE_MAX / sizeof(double)) {
                    bgen_error("batch scratch size overflow");
                    goto err;
                }
                if (needed > f->prob_scratch_size) {
                    double* p = realloc(f->prob_scratch, needed * sizeof(double));
                    if (!p) {
                        bgen_error("could not allocate padded batch scratch");
                        goto err;
                    }
                    f->prob_scratch = p;
                    f->prob_scratch_size = needed;
                }
                destination = f->prob_scratch;
            }
            int error = precision == 32 ? bgen_genotype_read32(&gt, destination) :
                                          bgen_genotype_read64(&gt, destination);
            if (error)
                goto err;
            if (gt.ncombs != width) {
                if (!padding_filled) {
                    /* Fill the untouched suffix once, keeping padding vectorizable. */
                    size_t elements = (size_t)(n_offsets - i) * f->nsamples * width;
                    if (precision == 32) {
                        for (size_t c = 0; c < elements; ++c)
                            ((float*)slice)[c] = NAN;
                    } else {
                        for (size_t c = 0; c < elements; ++c)
                            ((double*)slice)[c] = NAN;
                    }
                    padding_filled = true;
                }
                if (gt.ncombs) {
                    for (uint32_t s = 0; s < f->nsamples; ++s)
                        memcpy(slice + (size_t)s * width * item_size,
                               (char*)destination + (size_t)s * gt.ncombs * item_size,
                               (size_t)gt.ncombs * item_size);
                }
            }
        }
        bgen_genotype_release_buffers(&gt);
        continue;
err:
        bgen_genotype_release_buffers(&gt);
        return 1;
    }
    return 0;
}

int bgen_file_read_genotypes_batch64(struct bgen_file* f, uint64_t const* offsets,
                                     uint32_t n_offsets, uint32_t width, double* out,
                                     struct bgen_batch_metadata* metadata)
{
    return read_batch(f, offsets, n_offsets, width, out, 64, false, metadata);
}

int bgen_file_read_genotypes_batch32(struct bgen_file* f, uint64_t const* offsets,
                                     uint32_t n_offsets, uint32_t width, float* out,
                                     struct bgen_batch_metadata* metadata)
{
    return read_batch(f, offsets, n_offsets, width, out, 32, false, metadata);
}

int bgen_file_read_genotype_metadata(struct bgen_file* f, uint64_t const* offsets,
                                     uint32_t n_offsets, struct bgen_batch_metadata* metadata)
{
    return read_batch(f, offsets, n_offsets, 0, NULL, 0, false, metadata);
}

int bgen_file_read_genotypes_batch(struct bgen_file* f, uint64_t const* offsets,
                                   uint32_t n_offsets, uint32_t ncombs, double* out)
{
    return read_batch(f, offsets, n_offsets, ncombs, out, 64, true, NULL);
}

int bgen_file_read_ncombs_batch(struct bgen_file* f, uint64_t const* offsets,
                                uint32_t n_offsets, uint32_t* ncombs_out)
{
    if (n_offsets && !ncombs_out) {
        bgen_error("batch ncombs output is NULL");
        return 1;
    }
    struct bgen_batch_metadata metadata = {0};
    metadata.ncombs = ncombs_out;
    return bgen_file_read_genotype_metadata(f, offsets, n_offsets, &metadata);
}

int bgen_file_read_genotypes_batch_padded(struct bgen_file* f, uint64_t const* offsets,
                                          uint32_t n_offsets, uint32_t width, double* out)
{
    return bgen_file_read_genotypes_batch64(f, offsets, n_offsets, width, out, NULL);
}

static struct bgen_file* bgen_file_create(char const* filepath)
{
    struct bgen_file* bgen = malloc(sizeof(struct bgen_file));
    bgen->filepath = bgen_strdup(filepath);
    bgen->stream = NULL;
    bgen->nvariants = 0;
    bgen->nsamples = 0;
    bgen->compression = 0;
    bgen->layout = 0;
    bgen->contain_sample = 0;
    bgen->samples_start = 0;
    bgen->variants_start = 0;
    bgen->chunk_scratch = NULL;
    bgen->chunk_scratch_size = 0;
    bgen->ploidy_scratch = NULL;
    bgen->prob_scratch = NULL;
    bgen->prob_scratch_size = 0;
    bgen->zlib_ctx = NULL;
    bgen->zstd_ctx = NULL;

    if (!(bgen->stream = stream_handle_open(bgen->filepath))) {
        bgen_perror("could not open file %s", bgen->filepath);
        bgen_file_close(bgen);
        return NULL;
    }

    return bgen;
}

/*
 * Read the header block defined as follows:
 *
 *   header length: 4 bytes
 *   number of variants: 4 bytes
 *   number of samples: 4 bytes
 *   magic number: 4 bytes
 *   unused space: header length minus 20 bytes
 *   bgen flags: 4 bytes
 */
static int bgen_file_read_header(struct bgen_file* bgen)
{
    uint32_t header_length = 0;
    uint32_t magic_number = 0;
    uint32_t flags = 0;

    if (bgen_stream_fread(bgen->stream, &header_length, sizeof(header_length)) != 1) {
        bgen_perror_eof(stream_handle_eof(bgen->stream), "could not read header length");
        return 1;
    }
    if (header_length < 20 || (uint64_t)header_length + 4 > (uint64_t)bgen->variants_start) {
        bgen_error("invalid bgen header length or variant offset");
        return 1;
    }

    if (bgen_stream_fread(bgen->stream, &bgen->nvariants, sizeof(bgen->nvariants)) != 1) {
        bgen_perror_eof(stream_handle_eof(bgen->stream), "could not read number of variants");
        return 1;
    }

    if (bgen_stream_fread(bgen->stream, &bgen->nsamples, sizeof(bgen->nsamples)) != 1) {
        bgen_perror_eof(stream_handle_eof(bgen->stream), "could not read number of samples");
        return 1;
    }

    if (bgen_stream_fread(bgen->stream, &magic_number, sizeof(magic_number)) != 1) {
        bgen_perror_eof(stream_handle_eof(bgen->stream), "could not read magic number");
        return 1;
    }

    if (magic_number != 1852139362 && magic_number != 0) {
        bgen_error("magic number mismatch");
        return 1;
    }

    if (stream_handle_seek(bgen->stream, header_length - 20, SEEK_CUR)) {
        bgen_perror("fseek error while reading bgen file");
        return 1;
    }

    if (bgen_stream_fread(bgen->stream, &flags, sizeof(flags)) != 1) {
        bgen_perror_eof(stream_handle_eof(bgen->stream), "could not read bgen flags");
        return 1;
    }

    bgen->compression = flags & 3;
    bgen->layout = (flags & (15 << 2)) >> 2;
    bgen->contain_sample = ((flags & ((uint32_t)1 << 31)) >> 31) == 1;
    if ((flags & ~UINT32_C(0x8000003f)) != 0 ||
        (bgen->layout != 1 && bgen->layout != 2) || bgen->compression > 2 ||
        (bgen->layout == 1 && bgen->compression > 1)) {
        bgen_error("unsupported bgen flags: %" PRIu32, flags);
        return 1;
    }

    return 0;
}
