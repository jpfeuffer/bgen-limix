#include "file.h"
#include "bgen/file.h"
#include "bgen/genotype.h"
#include "bstring.h"
#include "free.h"
#include "genotype.h"
#include "io.h"
#include "layout1.h"
#include "layout2.h"
#include "mem.h"
#include "report.h"
#include "samples.h"
#include "stream.h"
#include "strdup.h"
#include "unzlib.h"
#include "unzstd.h"
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

struct bgen_file
{
    char*    filepath;
    FILE*    stream;
    uint32_t nvariants;
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
    if (fread(&bgen->variants_start, 4, 1, bgen->stream) != 1) {
        bgen_perror_eof(bgen->stream, "could not read the `variants_start` field");
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
    if ((bgen->samples_start = bgen_ftell(bgen->stream)) < 0) {
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
    if (bgen->stream != NULL && fclose(bgen->stream))
        bgen_perror("could not close %s file", bgen->filepath);
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

    if (bgen_fseek(bgen->stream, bgen->samples_start, SEEK_SET)) {
        bgen_perror("could not fseek to `samples_start`");
        return NULL;
    }

    if (!bgen->contain_sample) {
        bgen_warning("file does not contain sample ids");
        return NULL;
    }

    struct bgen_samples* samples = bgen_samples_create(bgen->nsamples);

    uint32_t block_size = 0;
    if (fread(&block_size, sizeof(block_size), 1, bgen->stream) != 1) {
        bgen_perror_eof(bgen->stream, "could not read block size");
        goto err;
    }

    block = malloc(block_size - sizeof(block_size));
    if (fread(block, block_size - sizeof(block_size), 1, bgen->stream) != 1) {
        bgen_perror_eof(bgen->stream, "could not read samples block");
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

    if ((bgen->variants_start = bgen_ftell(bgen->stream)) < 0) {
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

struct bgen_genotype* bgen_file_open_genotype(struct bgen_file* bgen, uint64_t genotype_offset)
{
    struct bgen_genotype* genotype = bgen_genotype_create();
    genotype->layout = bgen->layout;
    genotype->offset = genotype_offset;

    if (genotype_offset > INT64_MAX) {
        bgen_error("variant offset overflow");
        goto err;
    }

    if (bgen_fseek(bgen_file_stream(bgen), (int64_t)genotype_offset, SEEK_SET)) {
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

    return genotype;
err:
    bgen_genotype_close(genotype);
    return NULL;
}

FILE* bgen_file_stream(struct bgen_file const* bgen_file) { return bgen_file->stream; }

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
    if (bgen_fseek(bgen_file->stream, bgen_file->variants_start, SEEK_SET)) {
        bgen_perror("could not jump to variants start");
        return 1;
    }
    return 0;
}

int bgen_file_decompress_block(struct bgen_file* f, char** chunk, size_t* size)
{
    char* compressed = NULL;

    size_t compressed_length = 0;
    if (fread(&compressed_length, 4, 1, f->stream) < 1) {
        bgen_perror_eof(f->stream, "could not read compressed length");
        goto err;
    }
    if (compressed_length < 4) {
        bgen_error("compressed length too small (corrupted file?)");
        goto err;
    }
    compressed_length -= 4;

    size_t uncompressed_length = 0;
    if (fread(&uncompressed_length, 4, 1, f->stream) < 1) {
        bgen_perror_eof(f->stream, "could not read uncompressed length");
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
    if (fread(compressed, compressed_length, 1, f->stream) < 1) {
        bgen_perror_eof(f->stream, "could not read compressed data");
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

int bgen_file_read_genotypes_batch(struct bgen_file* f, uint64_t const* offsets,
                                   uint32_t n_offsets, uint32_t ncombs, double* out)
{
    if (n_offsets == 0) return 0;
    if (f->layout != 2) {
        bgen_error("bgen_file_read_genotypes_batch only supports layout 2");
        return 1;
    }
    if (f->compression == 0) {
        bgen_error("bgen_file_read_genotypes_batch requires compressed layout 2");
        return 1;
    }

    uint32_t nsamples = f->nsamples;

    for (uint32_t i = 0; i < n_offsets; ++i) {
        /* Skip seek if the file pointer is already at the right position.
         * This happens naturally when iterating over pre-sorted offsets where
         * the previous read left the pointer at the next metadata block
         * and the caller already seeked here (or we got lucky on consecutive
         * single-variant calls). */
        int64_t cur_pos = bgen_ftell(f->stream);
        if (cur_pos < 0 || cur_pos != (int64_t)offsets[i]) {
            if (bgen_fseek(f->stream, (int64_t)offsets[i], SEEK_SET)) {
                bgen_perror("batch: could not fseek to offset %" PRIu64, offsets[i]);
                return 1;
            }
        }

        /* Use stack-allocated genotype — no heap allocation per variant. */
        struct bgen_genotype gt;
        bgen_genotype_init(&gt);
        gt.layout   = f->layout;
        gt.nsamples = nsamples;

        if (bgen_layout2_read_header(f, &gt)) {
            bgen_error("batch: failed to read header at offset %" PRIu64, offsets[i]);
            return 1;
        }

        if (gt.ncombs != ncombs) {
            bgen_error("batch: ncombs mismatch at offset %" PRIu64
                       " (expected %u, got %u)",
                       offsets[i], ncombs, gt.ncombs);
            return 1;
        }

        bgen_layout2_read_genotype64(&gt, out + (uint64_t)i * nsamples * ncombs);
        /* gt.chunk and gt.ploidy_missingness are owned by f (scratch buffers).
         * gt is on the stack — nothing to free. */
    }
    return 0;
}

/* ---- Plan E: ncombs-scan batch (header-only, no probability decode) ---------- */

int bgen_file_read_ncombs_batch(struct bgen_file* f, uint64_t const* offsets,
                                uint32_t n_offsets, uint32_t* ncombs_out)
{
    if (n_offsets == 0) return 0;
    if (f->layout != 2) {
        bgen_error("bgen_file_read_ncombs_batch only supports layout 2");
        return 1;
    }
    if (f->compression == 0) {
        bgen_error("bgen_file_read_ncombs_batch requires compressed layout 2");
        return 1;
    }

    uint32_t nsamples = f->nsamples;

    for (uint32_t i = 0; i < n_offsets; ++i) {
        int64_t cur_pos = bgen_ftell(f->stream);
        if (cur_pos < 0 || cur_pos != (int64_t)offsets[i]) {
            if (bgen_fseek(f->stream, (int64_t)offsets[i], SEEK_SET)) {
                bgen_perror("ncombs_batch: could not fseek to offset %" PRIu64, offsets[i]);
                return 1;
            }
        }

        struct bgen_genotype gt;
        bgen_genotype_init(&gt);
        gt.layout   = f->layout;
        gt.nsamples = nsamples;

        if (bgen_layout2_read_header(f, &gt)) {
            bgen_error("ncombs_batch: failed to read header at offset %" PRIu64, offsets[i]);
            return 1;
        }

        ncombs_out[i] = gt.ncombs;
        /* Chunk owned by f->chunk_scratch — nothing to free. */
    }
    return 0;
}

/* ---- Plan E: padded batch (NaN-fills unused ncombs slots per variant) -------- */

int bgen_file_read_genotypes_batch_padded(struct bgen_file* f, uint64_t const* offsets,
                                          uint32_t n_offsets, uint32_t max_ncombs,
                                          double* out)
{
    if (n_offsets == 0) return 0;
    if (max_ncombs == 0) {
        bgen_error("bgen_file_read_genotypes_batch_padded: max_ncombs cannot be zero");
        return 1;
    }
    if (f->layout != 2) {
        bgen_error("bgen_file_read_genotypes_batch_padded only supports layout 2");
        return 1;
    }
    if (f->compression == 0) {
        bgen_error("bgen_file_read_genotypes_batch_padded requires compressed layout 2");
        return 1;
    }

    uint32_t nsamples = f->nsamples;

    /* Pre-fill the entire output buffer with NaN so that unused padding slots
     * are already correct without per-variant bookkeeping. */
    uint64_t total = (uint64_t)n_offsets * nsamples * max_ncombs;
    for (uint64_t k = 0; k < total; ++k) out[k] = NAN;

    for (uint32_t i = 0; i < n_offsets; ++i) {
        int64_t cur_pos = bgen_ftell(f->stream);
        if (cur_pos < 0 || cur_pos != (int64_t)offsets[i]) {
            if (bgen_fseek(f->stream, (int64_t)offsets[i], SEEK_SET)) {
                bgen_perror("padded_batch: could not fseek to offset %" PRIu64, offsets[i]);
                return 1;
            }
        }

        struct bgen_genotype gt;
        bgen_genotype_init(&gt);
        gt.layout   = f->layout;
        gt.nsamples = nsamples;

        if (bgen_layout2_read_header(f, &gt)) {
            bgen_error("padded_batch: failed to read header at offset %" PRIu64, offsets[i]);
            return 1;
        }

        if (gt.ncombs > max_ncombs) {
            bgen_error("padded_batch: variant at offset %" PRIu64
                       " has ncombs=%u which exceeds max_ncombs=%u",
                       offsets[i], gt.ncombs, max_ncombs);
            return 1;
        }

        double* slice = out + (uint64_t)i * nsamples * max_ncombs;

        if (gt.ncombs == max_ncombs) {
            /* No padding needed — decode directly into the output slice. */
            bgen_layout2_read_genotype64(&gt, slice);
        } else {
            /* ncombs < max_ncombs: decode into prob_scratch, then scatter each
             * sample's probabilities into the wider padded row. */
            uint64_t scratch_needed = (uint64_t)nsamples * gt.ncombs;
            if (scratch_needed > f->prob_scratch_size) {
                double* p = realloc(f->prob_scratch, scratch_needed * sizeof(double));
                if (!p) {
                    bgen_error("padded_batch: could not realloc prob_scratch");
                    return 1;
                }
                f->prob_scratch      = p;
                f->prob_scratch_size = scratch_needed;
            }

            bgen_layout2_read_genotype64(&gt, f->prob_scratch);

            /* Scatter: copy each sample's ncombs probabilities into the
             * appropriate slot; the remaining (max_ncombs - ncombs) slots
             * were already set to NaN by the pre-fill above. */
            uint32_t nc = gt.ncombs;
            for (uint32_t s = 0; s < nsamples; ++s) {
                memcpy(slice + (uint64_t)s * max_ncombs,
                       f->prob_scratch + (uint64_t)s * nc,
                       nc * sizeof(double));
            }
        }
    }
    return 0;
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

    if (!(bgen->stream = bgen_stream_open(bgen->filepath))) {
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

    if (fread(&header_length, sizeof(header_length), 1, bgen->stream) != 1) {
        bgen_perror_eof(bgen->stream, "could not read header length");
        return 1;
    }

    if (fread(&bgen->nvariants, sizeof(bgen->nvariants), 1, bgen->stream) != 1) {
        bgen_perror_eof(bgen->stream, "could not read number of variants");
        return 1;
    }

    if (fread(&bgen->nsamples, sizeof(bgen->nsamples), 1, bgen->stream) != 1) {
        bgen_perror_eof(bgen->stream, "could not read number of samples");
        return 1;
    }

    if (fread(&magic_number, sizeof(magic_number), 1, bgen->stream) != 1) {
        bgen_perror_eof(bgen->stream, "could not read magic number");
        return 1;
    }

    if (magic_number != 1852139362)
        bgen_warning("magic number mismatch");

    if (bgen_fseek(bgen->stream, header_length - 20, SEEK_CUR)) {
        bgen_perror("fseek error while reading bgen file");
        return 1;
    }

    if (fread(&flags, sizeof(flags), 1, bgen->stream) != 1) {
        bgen_perror_eof(bgen->stream, "could not read bgen flags");
        return 1;
    }

    bgen->compression = flags & 3;
    bgen->layout = (flags & (15 << 2)) >> 2;
    bgen->contain_sample = ((flags & ((uint32_t)1 << 31)) >> 31) == 1;

    return 0;
}
