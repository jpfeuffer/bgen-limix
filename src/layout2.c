#include "layout2.h"
#include "bgen/file.h"
#include "bmath.h"
#include "file.h"
#include "free.h"
#include "genotype.h"
#include "mem.h"
#include "unzlib.h"
#include "unzstd.h"
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#define BIT(var, bit) ((var & (1 << bit)) != 0)

static void  read_phased_genotype64(struct bgen_genotype* genotype, double* probs);
static void  read_phased_genotype32(struct bgen_genotype* genotype, float* probs);
static void  read_unphased_genotype64(struct bgen_genotype* genotype, double* probs);
static void  read_unphased_genotype32(struct bgen_genotype* genotype, float* probs);

static inline uint8_t read_ploidy(uint8_t ploidy_miss) { return ploidy_miss & 127; }

static inline int read_missingness(uint8_t ploidy_miss) { return ploidy_miss >> 7; }

static inline int get_bit(char const* mem, uint64_t bit_idx)
{
    uint64_t bytes = bit_idx / 8;
    return BIT(*(mem + bytes), bit_idx % 8);
}

static inline void set_array_nan64(double* p, size_t n)
{
    for (size_t i = 0; i < n; ++i)
        p[i] = NAN;
}

static inline void set_array_nan32(float* p, size_t n)
{
    for (size_t i = 0; i < n; ++i)
        p[i] = NAN;
}

static int read_header(struct bgen_file* bgen_file, struct bgen_genotype* genotype,
                       bool compact)
{
    uint32_t nsamples = 0;

    char const* chunk_ptr = NULL;
    char*       chunk     = NULL;
    bool        chunk_is_scratch = false;
    size_t size = 0;

    if (bgen_file_compression(bgen_file) > 0) {
        /* Plan A+B: decompress into file's scratch buffer, reusing zlib/zstd context. */
        if (bgen_file_decompress_block(bgen_file, &chunk, &size))
            goto err;
        chunk_is_scratch = true;
    } else {
        uint32_t length = 0;
        if (bgen_stream_fread(bgen_file_stream(bgen_file), &length, sizeof(length)) < 1) {
            bgen_perror_eof(stream_handle_eof(bgen_file_stream(bgen_file)),
                            "could not read genotype block length");
            goto err;
        }
        size = length;
        chunk = malloc(size ? size : 1);
        if (!chunk) { bgen_error("could not malloc uncompressed chunk"); goto err; }

        if (bgen_stream_fread(bgen_file_stream(bgen_file), chunk, size) < 1) {
            bgen_perror_eof(stream_handle_eof(bgen_file_stream(bgen_file)), "could not read chunk");
            goto err;
        }

    }

    if (size < 10) {
        bgen_error("layout 2 genotype header is truncated");
        goto err;
    }
    chunk_ptr = chunk;
    bgen_memfread(&nsamples, &chunk_ptr, sizeof(nsamples));
    if (nsamples != bgen_file_nsamples(bgen_file) || nsamples > size - 10) {
        bgen_error("layout 2 sample count mismatch or truncated ploidy");
        goto err;
    }
    uint16_t nalleles = 0;
    bgen_memfread(&nalleles, &chunk_ptr, 2);
    uint8_t min_ploidy = 0;
    bgen_memfread(&min_ploidy, &chunk_ptr, 1);
    uint8_t max_ploidy = 0;
    bgen_memfread(&max_ploidy, &chunk_ptr, 1);
    genotype->min_ploidy = min_ploidy;
    genotype->max_ploidy = max_ploidy;

    uint8_t* plo_miss = bgen_file_ploidy_scratch(bgen_file);
    uint8_t const* encoded_ploidy = (uint8_t const*)chunk_ptr;
    chunk_ptr += nsamples;

    uint8_t phased = 0;
    bgen_memfread(&phased, &chunk_ptr, 1);
    uint8_t nbits = 0;
    bgen_memfread(&nbits, &chunk_ptr, 1);

    genotype->nsamples            = nsamples;
    genotype->nalleles            = nalleles;
    genotype->phased              = phased;
    genotype->nbits               = nbits;
    genotype->ploidy_missingness  = plo_miss;
    genotype->owns_ploidy         = false;   /* owned by bgen_file scratch */

    if (genotype->nalleles == 0) {
        bgen_error("`nalleles` cannot be zero");
        goto err;
    }
    if (phased > 1 || nbits < 1 || nbits > 32 || min_ploidy > max_ploidy ||
        max_ploidy > 63) {
        bgen_error("invalid layout 2 phase, bit depth or ploidy range");
        goto err;
    }

    if (compact && (min_ploidy != max_ploidy || !nsamples)) {
        unsigned actual_min = 63, actual_max = 0;
        for (uint32_t i = 0; i < nsamples; ++i) {
            unsigned p = read_ploidy(encoded_ploidy[i]);
            if (p < min_ploidy || p > max_ploidy) {
                bgen_error("sample ploidy outside genotype range");
                goto err;
            }
            if (p < actual_min) actual_min = p;
            if (p > actual_max) actual_max = p;
        }
        min_ploidy = (uint8_t)(nsamples ? actual_min : 0);
        max_ploidy = (uint8_t)actual_max;
        genotype->min_ploidy = min_ploidy;
        genotype->max_ploidy = max_ploidy;
    }
    if (phased)
        genotype->ncombs = (unsigned)nalleles * (unsigned)genotype->max_ploidy;
    else if (bgen_checked_ncombs(max_ploidy, nalleles, &genotype->ncombs))
        goto err;
    if (compact && !nsamples)
        genotype->ncombs = 0;

    uint64_t probability_bits = 0;
    size_t payload_size = size - (size_t)(chunk_ptr - chunk);
    if (min_ploidy == max_ploidy) {
        /* Keep common fixed-ploidy validation vectorizable while copying. */
        unsigned invalid_ploidy = 0;
        for (uint32_t i = 0; i < nsamples; ++i) {
            uint8_t value = encoded_ploidy[i];
            plo_miss[i] = value;
            invalid_ploidy |= read_ploidy(value) != max_ploidy;
        }
        if (invalid_ploidy) {
            bgen_error("sample ploidy outside genotype range");
            goto err;
        }
        uint64_t bits_per_sample = (uint64_t)nbits * (phased ?
            (unsigned)max_ploidy * (unsigned)(nalleles - 1) :
            (genotype->ncombs ? genotype->ncombs - 1 : 0));
        if (bits_per_sample && nsamples > (uint64_t)payload_size * 8 / bits_per_sample) {
            bgen_error("layout 2 probability data is truncated");
            goto err;
        }
        probability_bits = (uint64_t)nsamples * bits_per_sample;
    } else {
        unsigned counts[64];
        for (unsigned ploidy = min_ploidy; ploidy <= max_ploidy; ++ploidy) {
            if (phased)
                counts[ploidy] = ploidy * (unsigned)(nalleles - 1);
            else {
                if (bgen_checked_ncombs(ploidy, nalleles, &counts[ploidy]))
                    goto err;
                --counts[ploidy];
            }
        }
        for (uint32_t i = 0; i < nsamples; ++i) {
            uint8_t value = encoded_ploidy[i];
            plo_miss[i] = value;
            unsigned ploidy = read_ploidy(value);
            if (ploidy < min_ploidy || ploidy > max_ploidy) {
                bgen_error("sample ploidy outside genotype range");
                goto err;
            }
            probability_bits += (uint64_t)nbits * counts[ploidy];
            if ((probability_bits + 7) / 8 > payload_size) {
                bgen_error("layout 2 probability data is truncated");
                goto err;
            }
        }
    }
    if ((probability_bits + 7) / 8 != payload_size) {
        bgen_error("layout 2 probability length mismatch");
        goto err;
    }

    genotype->chunk     = chunk;
    genotype->chunk_ptr = chunk_ptr;
    genotype->owns_chunk = !chunk_is_scratch;  /* false = file owns scratch; true = malloc'd */

    return 0;

err:
    if (!chunk_is_scratch) bgen_free(chunk);
    genotype->chunk              = NULL;
    genotype->ploidy_missingness = NULL;
    return 1;
}

int bgen_layout2_read_header(struct bgen_file* file, struct bgen_genotype* genotype)
{
    return read_header(file, genotype, false);
}

int bgen_layout2_read_header_compact(struct bgen_file* file, struct bgen_genotype* genotype)
{
    return read_header(file, genotype, true);
}

void bgen_layout2_read_genotype64(struct bgen_genotype* genotype, double* probs)
{
    if (genotype->ncombs == 0)
        return;
    if (genotype->phased)
        read_phased_genotype64(genotype, probs);
    else
        read_unphased_genotype64(genotype, probs);
}

void bgen_layout2_read_genotype32(struct bgen_genotype* genotype, float* probs)
{
    if (genotype->ncombs == 0)
        return;
    if (genotype->phased)
        read_phased_genotype32(genotype, probs);
    else
        read_unphased_genotype32(genotype, probs);
}

/* ── Plan D: fast byte-aligned read for nbits divisible by 8 ────────────────
 * When nbits is a multiple of 8, each stored value occupies exactly nbits/8
 * contiguous bytes in the chunk.  sample_start always stays byte-aligned
 * because the per-sample advance (nbits × ncombs-per-sample) is also a
 * multiple of 8 bits.  We replace the inner bit-extraction loop with a single
 * memcpy, which compilers auto-vectorise on modern targets.
 * ────────────────────────────────────────────────────────────────────────── */
#define READ_UINT_BYTES(chunk_ptr, bit_offset, nbytes, out_u64)              \
    do {                                                                     \
        uint64_t _v = 0;                                                     \
        memcpy(&_v, (chunk_ptr) + (bit_offset) / 8, (nbytes));              \
        (out_u64) = _v;                                                      \
    } while (0)

#define MAKE_READ_PHASED_GENOTYPE(BITS, FPTYPE)                                               \
    static void read_phased_genotype##BITS(struct bgen_genotype* genotype, FPTYPE* probs)     \
    {                                                                                         \
        unsigned nbits = genotype->nbits;                                                     \
        unsigned nalleles = genotype->nalleles;                                               \
        uint8_t  max_ploidy = genotype->max_ploidy;                                           \
        uint64_t denom = ((uint64_t)1 << nbits) - 1;                                         \
        int      aligned = (nbits % 8 == 0);                                                  \
        unsigned nbytes  = nbits / 8;                                                         \
                                                                                              \
        uint64_t sample_start = 0;                                                            \
        for (uint32_t j = 0; j < genotype->nsamples; ++j) {                                   \
            unsigned ploidy = read_ploidy(genotype->ploidy_missingness[j]);                   \
            FPTYPE*  pend = probs + max_ploidy * nalleles;                                    \
                                                                                              \
            if (read_missingness(genotype->ploidy_missingness[j]) != 0) {                     \
                set_array_nan##BITS(probs, (size_t)(pend - probs));                           \
                probs = pend;                                                                 \
                sample_start += nbits * ploidy * (nalleles - 1);                              \
                continue;                                                                     \
            }                                                                                 \
                                                                                              \
            uint64_t haplo_start = 0;                                                         \
            for (uint8_t i = 0; i < ploidy; ++i) {                                            \
                uint64_t uip_sum = 0;                                                         \
                uint64_t allele_start = 0;                                                    \
                for (uint16_t ii = 0; ii < nalleles - 1; ++ii) {                              \
                    uint64_t ui_prob = 0;                                                     \
                    uint64_t offset = sample_start + haplo_start + allele_start;              \
                    if (aligned) {                                                             \
                        READ_UINT_BYTES(genotype->chunk_ptr, offset, nbytes, ui_prob);        \
                    } else {                                                                   \
                        for (uint8_t bi = 0; bi < nbits; ++bi) {                              \
                            if (get_bit(genotype->chunk_ptr, bi + offset))                    \
                                ui_prob |= ((uint64_t)1 << bi);                               \
                        }                                                                     \
                    }                                                                         \
                    *probs = (FPTYPE)ui_prob / (FPTYPE)denom;                                 \
                    ++probs;                                                                  \
                    uip_sum += ui_prob;                                                       \
                    allele_start += nbits;                                                    \
                }                                                                             \
                *probs = (FPTYPE)(denom - uip_sum) / (FPTYPE)denom;                           \
                ++probs;                                                                      \
                haplo_start += nbits * (nalleles - 1);                                        \
            }                                                                                 \
            sample_start += nbits * ploidy * (nalleles - 1);                                  \
            set_array_nan##BITS(probs, (size_t)(pend - probs));                               \
            probs = pend;                                                                     \
        }                                                                                     \
    }

MAKE_READ_PHASED_GENOTYPE(64, double)
MAKE_READ_PHASED_GENOTYPE(32, float)

#define MAKE_UNPHASED_GENOTYPE(BITS, FPTYPE)                                                  \
    static void read_unphased_genotype##BITS(struct bgen_genotype* genotype, FPTYPE* probs)   \
    {                                                                                         \
        uint8_t  nbits = genotype->nbits;                                                     \
        uint16_t nalleles = genotype->nalleles;                                               \
        uint64_t denom = ((uint64_t)1 << nbits) - 1;                                         \
        int      aligned = (nbits % 8 == 0);                                                  \
        unsigned nbytes  = nbits / 8;                                                         \
        unsigned max_ncombs = genotype->ncombs;                                               \
        unsigned combinations[64];                                                           \
        for (unsigned ploidy = genotype->min_ploidy; ploidy <= genotype->max_ploidy; ++ploidy) \
            bgen_checked_ncombs(ploidy, nalleles, &combinations[ploidy]);                     \
                                                                                              \
        uint64_t sample_start = 0;                                                            \
        for (uint32_t j = 0; j < genotype->nsamples; ++j) {                                   \
            FPTYPE*  pend = probs + max_ncombs;                                               \
            uint8_t  ploidy = read_ploidy(genotype->ploidy_missingness[j]);                   \
            unsigned ncombs = combinations[ploidy];                                           \
                                                                                              \
            if (read_missingness(genotype->ploidy_missingness[j]) != 0) {                     \
                set_array_nan##BITS(probs, (size_t)(pend - probs));                           \
                probs = pend;                                                                 \
                sample_start += (uint64_t)nbits * (ncombs - 1);                               \
                continue;                                                                     \
            }                                                                                 \
                                                                                              \
            uint64_t uip_sum = 0;                                                             \
            uint64_t geno_start = 0;                                                          \
            for (unsigned i = 0; i < ncombs - 1; ++i) {                                      \
                uint64_t ui_prob = 0;                                                         \
                uint64_t offset = sample_start + geno_start;                                  \
                if (aligned) {                                                                 \
                    READ_UINT_BYTES(genotype->chunk_ptr, offset, nbytes, ui_prob);             \
                } else {                                                                       \
                    for (uint8_t bi = 0; bi < (uint8_t)nbits; ++bi) {                         \
                        if (get_bit(genotype->chunk_ptr, bi + offset))                        \
                            ui_prob |= ((uint64_t)1 << bi);                                   \
                    }                                                                         \
                }                                                                             \
                *probs = (FPTYPE)ui_prob / (FPTYPE)denom;                                     \
                ++probs;                                                                      \
                uip_sum += ui_prob;                                                          \
                geno_start += nbits;                                                          \
            }                                                                                 \
            *probs = (FPTYPE)(denom - uip_sum) / (FPTYPE)denom;                               \
            ++probs;                                                                          \
            sample_start += (uint64_t)nbits * (ncombs - 1);                                   \
            set_array_nan##BITS(probs, (size_t)(pend - probs));                               \
            probs = pend;                                                                     \
        }                                                                                     \
    }

MAKE_UNPHASED_GENOTYPE(64, double)
MAKE_UNPHASED_GENOTYPE(32, float)
