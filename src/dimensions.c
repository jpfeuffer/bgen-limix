#include "bgen/file.h"
#include "bgen/genotype.h"
#include "bmath.h"
#include "file.h"
#include "genotype.h"
#include "layout2.h"
#include "report.h"
#include "zlib-ng.h"
#include <stdio.h>
#include <string.h>

struct prefix_reader {
    struct bgen_file* file;
    bool compressed;
    zng_stream* zlib;
    uint32_t remaining;
    unsigned char input[8192];
};

static unsigned uint_le(unsigned char const* p, unsigned bytes)
{
    unsigned value = 0;
    for (unsigned i = 0; i < bytes; ++i)
        value |= (unsigned)p[i] << (8 * i);
    return value;
}

static int read_prefix(struct prefix_reader* r, void* out, uint32_t count)
{
    if (!r->compressed) {
        if (count > r->remaining ||
            !bgen_stream_fread(bgen_file_stream(r->file), out, count)) {
            bgen_error("could not read genotype dimension prefix");
            return 1;
        }
        r->remaining -= count;
        return 0;
    }
    zng_stream* z = r->zlib;
    z->next_out = out;
    z->avail_out = count;
    while (z->avail_out) {
        if (!z->avail_in) {
            size_t n = r->remaining < sizeof(r->input) ? r->remaining : sizeof(r->input);
            if (!n || !bgen_stream_fread(bgen_file_stream(r->file), r->input, n)) {
                bgen_error("truncated zlib genotype dimension prefix");
                return 1;
            }
            r->remaining -= (uint32_t)n;
            z->next_in = r->input;
            z->avail_in = (uint32_t)n;
        }
        uint32_t before_in = z->avail_in, before_out = z->avail_out;
        int error = zng_inflate(z, Z_NO_FLUSH);
        if ((error != Z_OK && error != Z_STREAM_END) ||
            (error == Z_STREAM_END && z->avail_out) ||
            (before_in == z->avail_in && before_out == z->avail_out)) {
            bgen_error("could not inflate genotype dimension prefix");
            return 1;
        }
    }
    return 0;
}

static int discard_ploidy(struct prefix_reader* r, uint32_t count)
{
    if (!r->compressed) {
        if (count > r->remaining ||
            stream_handle_seek(bgen_file_stream(r->file), count, SEEK_CUR)) {
            bgen_error("could not seek past genotype ploidy");
            return 1;
        }
        r->remaining -= count;
        return 0;
    }
    unsigned char discard[4096];
    while (count) {
        uint32_t n = count < sizeof(discard) ? count : (uint32_t)sizeof(discard);
        if (read_prefix(r, discard, n))
            return 1;
        count -= n;
    }
    return 0;
}

static int add_dimensions(struct bgen_dimensions* dimensions, unsigned k, unsigned p,
                          bool phased, enum bgen_probability_model model)
{
    if (model == BGEN_MODEL_UNPHASED || model == BGEN_MODEL_MIXED ||
        (model == BGEN_MODEL_AUTO && !phased)) {
        unsigned count;
        if (bgen_checked_ncombs(p, k, &count))
            return 1;
        if (count > dimensions->gp_width)
            dimensions->gp_width = count;
    }
    if (model == BGEN_MODEL_PHASED || model == BGEN_MODEL_MIXED ||
        (model == BGEN_MODEL_AUTO && phased)) {
        unsigned count = p * k;
        if (!count) count = 1; /* Keep phase presence visible at zero ploidy. */
        if (count > dimensions->hp_width)
            dimensions->hp_width = count;
    }
    return 0;
}

int bgen_file_read_dimensions_with_ploidy(struct bgen_file* file, uint64_t const* offsets,
                                          uint32_t n_offsets, bool actual,
                                          enum bgen_probability_model model, int override,
                                          struct bgen_dimensions* out)
{
    if (!file || !out || (n_offsets && !offsets) ||
        model < BGEN_MODEL_AUTO || model > BGEN_MODEL_MIXED || override < -1 || override > 63) {
        bgen_error("invalid genotype dimension arguments");
        return 1;
    }
    struct bgen_dimensions result = {0};
    zng_stream zlib = {0};
    bool zlib_initialized = false;
    unsigned layout = bgen_file_layout(file), compression = bgen_file_compression(file);
    uint32_t nsamples = bgen_file_nsamples(file);
    if (override >= 0) actual = false;
    for (uint32_t i = 0; i < n_offsets; ++i) {
        if (bgen_file_seek_genotype(file, offsets[i]))
            goto err;
        if (layout == 1) {
            if (add_dimensions(&result, 2, override >= 0 ? (unsigned)override : 2, false, model))
                goto err;
            continue;
        }
        if (compression == 2) {
            struct bgen_genotype gt;
            bgen_genotype_init(&gt);
            if ((actual || override >= 0) ? bgen_layout2_read_header_compact(file, &gt) :
                         bgen_layout2_read_header(file, &gt)) {
                bgen_genotype_release_buffers(&gt);
                goto err;
            }
            unsigned maximum = override >= 0 ? (unsigned)override : gt.max_ploidy;
            int error = add_dimensions(&result, gt.nalleles, maximum, gt.phased != 0, model);
            bgen_genotype_release_buffers(&gt);
            if (error) goto err;
            continue;
        }
        struct prefix_reader reader = {.file = file, .compressed = compression == 1,
                                        .zlib = &zlib};
        unsigned char lengths[8], header[8], trailer[2];
        unsigned length_bytes = reader.compressed ? 8 : 4;
        if (!bgen_stream_fread(bgen_file_stream(file), lengths, length_bytes)) {
            bgen_error("could not read genotype dimension block length");
            goto err;
        }
        uint32_t size = uint_le(lengths, 4);
        uint64_t prefix_size = (uint64_t)nsamples + 10;
        if (reader.compressed) {
            if (size < 4 || uint_le(lengths + 4, 4) < prefix_size) {
                bgen_error("invalid compressed genotype dimension block length");
                goto err;
            }
            reader.remaining = size - 4;
            if (!zlib_initialized) {
                if (zng_inflateInit(&zlib) != Z_OK) {
                    bgen_error("could not initialize dimension inflater");
                    goto err;
                }
                zlib_initialized = true;
            } else if (zng_inflateReset(&zlib) != Z_OK) {
                bgen_error("could not reset dimension inflater");
                goto err;
            }
            zlib.avail_in = 0;
        } else {
            if (size < prefix_size) {
                bgen_error("invalid raw genotype dimension block length");
                goto err;
            }
            reader.remaining = size;
        }
        if (read_prefix(&reader, header, 8))
            goto err;
        unsigned k = uint_le(header + 4, 2), minimum = header[6], maximum = header[7];
        if (uint_le(header, 4) != nsamples || !k || minimum > maximum || maximum > 63) {
            bgen_error("invalid genotype dimension sample, allele or ploidy header");
            goto err;
        }
        if (actual) {
            unsigned char* ploidies = bgen_file_ploidy_scratch(file);
            if (read_prefix(&reader, ploidies, nsamples))
                goto err;
            unsigned actual_maximum = 0;
            for (uint32_t s = 0; s < nsamples; ++s) {
                unsigned p = ploidies[s] & 127;
                if (p < minimum || p > maximum) {
                    bgen_error("sample ploidy outside dimension header range");
                    goto err;
                }
                if (p > actual_maximum) actual_maximum = p;
            }
            maximum = actual_maximum;
        } else if (model == BGEN_MODEL_AUTO && discard_ploidy(&reader, nsamples))
            goto err;
        bool phased = false;
        if (model == BGEN_MODEL_AUTO) {
            if (read_prefix(&reader, trailer, 2))
                goto err;
            if (trailer[0] > 1 || trailer[1] < 1 || trailer[1] > 32) {
                bgen_error("invalid genotype dimension phase or bit depth");
                goto err;
            }
            phased = trailer[0] != 0;
        }
        if (add_dimensions(&result, k, override >= 0 ? (unsigned)override : maximum, phased, model))
            goto err;
    }
    result.decode_width = result.gp_width > result.hp_width ? result.gp_width : result.hp_width;
    if (zlib_initialized) zng_inflateEnd(&zlib);
    *out = result;
    return 0;
err:
    if (zlib_initialized) zng_inflateEnd(&zlib);
    return 1;
}

int bgen_file_read_dimensions(struct bgen_file* file, uint64_t const* offsets,
                              uint32_t n_offsets, bool actual,
                              enum bgen_probability_model model,
                              struct bgen_dimensions* out)
{
    return bgen_file_read_dimensions_with_ploidy(file, offsets, n_offsets, actual,
                                                model, -1, out);
}
