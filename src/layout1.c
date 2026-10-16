#include "layout1.h"
#include "bgen/bgen.h"
#include "file.h"
#include "free.h"
#include "genotype.h"
#include "mem.h"
#include "report.h"
#include "unzlib.h"
#include <math.h>
#include <stdlib.h>

static void  read_unphased64(struct bgen_genotype* vg, double* probs);
static void  read_unphased32(struct bgen_genotype* vg, float* probs);
static char* decompress(struct bgen_file* bgen_file, size_t size);

int bgen_layout1_read_header(struct bgen_file* bgen_file, struct bgen_genotype* genotype)
{
    char* chunk = NULL;
    uint32_t nsamples = bgen_file_nsamples(bgen_file);
    if ((uint64_t)nsamples * 6 > SIZE_MAX) {
        bgen_error("layout 1 genotype size overflow");
        return 1;
    }
    size_t size = (size_t)nsamples * 6;

    if (bgen_file_compression(bgen_file) > 0) {
        chunk = decompress(bgen_file, size);
    } else {
        chunk = malloc(size ? size : 1);
        if (!chunk) {
            bgen_error("could not allocate layout 1 genotype");
            return 1;
        }

        stream_handle* stream = bgen_file_stream(bgen_file);
        if (bgen_stream_fread(stream, chunk, size) < 1) {
            bgen_perror_eof(stream_handle_eof(stream), "could not read layout 1 genotype");
            bgen_free(chunk);
            return 1;
        }
    }
    if (!chunk)
        return 1;

    genotype->chunk = chunk;
    genotype->owns_chunk = true;
    genotype->ploidy_missingness = malloc(nsamples ? nsamples : 1);
    genotype->owns_ploidy = true;
    if (!genotype->ploidy_missingness) {
        bgen_error("could not allocate layout 1 ploidy");
        return 1;
    }
    for (uint32_t j = 0; j < nsamples; ++j) {
        bool missing = true;
        for (size_t i = 0; i < 6; ++i)
            missing = missing && chunk[(size_t)j * 6 + i] == 0;
        genotype->ploidy_missingness[j] = 2 | (missing ? 128 : 0);
    }
    genotype->nsamples = nsamples;
    genotype->nalleles = 2;
    genotype->ncombs = 3;
    genotype->min_ploidy = 2;
    genotype->max_ploidy = 2;
    genotype->phased = false;
    genotype->chunk_ptr = chunk;

    return 0;
}

void bgen_layout1_read_genotype64(struct bgen_genotype* genotype, double* probs)
{
    read_unphased64(genotype, probs);
}

void bgen_layout1_read_genotype32(struct bgen_genotype* genotype, float* probs)
{
    read_unphased32(genotype, probs);
}

#define MAKE_READ_UNPHASED(BITS, FPTYPE)                                                      \
    static void read_unphased##BITS(struct bgen_genotype* genotype, FPTYPE* probs)            \
    {                                                                                         \
        uint16_t ui_prob = 0;                                                                 \
        FPTYPE   denom = 32768;                                                               \
                                                                                              \
        char const* restrict chunk = genotype->chunk_ptr;                                     \
                                                                                              \
        for (uint32_t j = 0; j < genotype->nsamples; ++j) {                                   \
            unsigned int ui_prob_sum = 0;                                                     \
                                                                                              \
            for (size_t i = 0; i < 3; ++i) {                                                  \
                bgen_memfread(&ui_prob, &chunk, 2);                                           \
                probs[(size_t)j * 3 + i] = ui_prob / denom;                                   \
                ui_prob_sum += ui_prob;                                                       \
            }                                                                                 \
                                                                                              \
            if (ui_prob_sum == 0) {                                                           \
                for (size_t i = 0; i < 3; ++i)                                                \
                    probs[(size_t)j * 3 + i] = NAN;                                           \
            }                                                                                 \
        }                                                                                     \
    }

MAKE_READ_UNPHASED(64, double)
MAKE_READ_UNPHASED(32, float)

static char* decompress(struct bgen_file* bgen_file, size_t size)
{
    uint32_t compressed_length = 0;
    stream_handle* stream = bgen_file_stream(bgen_file);

    if (bgen_stream_fread(stream, &compressed_length, sizeof(compressed_length)) < 1) {
        bgen_perror_eof(stream_handle_eof(stream), "could not read chunk size");
        return NULL;
    }

    char* compressed_chunk = malloc(compressed_length ? compressed_length : 1);
    if (!compressed_chunk) {
        bgen_error("could not allocate compressed layout 1 genotype");
        return NULL;
    }

    if (bgen_stream_fread(stream, compressed_chunk, compressed_length) < 1) {
        bgen_perror_eof(stream_handle_eof(stream), "could not read compressed chunk");
        bgen_free(compressed_chunk);
        return NULL;
    }

    if (bgen_file_compression(bgen_file) != 1) {
        bgen_error("compression flag should be 1; not %u", bgen_file_compression(bgen_file));
        bgen_free(compressed_chunk);
        return NULL;
    }

    char* chunk = malloc(size ? size : 1);
    if (!chunk) {
        bgen_error("could not allocate layout 1 genotype");
        bgen_free(compressed_chunk);
        return NULL;
    }
    if (bgen_unzlib(compressed_chunk, compressed_length, &chunk, &size)) {
        bgen_free(chunk);
        chunk = NULL;
    }

    bgen_free(compressed_chunk);

    return chunk;
}
