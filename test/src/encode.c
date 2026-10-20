#include "bgen/bgen.h"
#undef NDEBUG
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static void sweep(void)
{
    double const gp[3] = {0.5, 0.25, 0.25};
    double const hp[4] = {0.5, 0.5, 0.25, 0.75};
    float const gp32[3] = {0.5f, 0.25f, 0.25f};
    float const hp32[4] = {0.5f, 0.5f, 0.25f, 0.75f};
    int64_t gt[63] = {0};
    uint8_t p = 2, missing = 0;
    uint8_t out[512];
    uint64_t offsets[2];
    size_t size = 0;
    struct bgen_encode_error error;
    struct bgen_layout2_variant v = {0};
    v.nsamples = 1;
    v.nalleles = 2;
    v.ploidy = &p;
    v.ploidy_bytes = 1;
    v.missing = &missing;
    v.missing_bytes = 1;
    for (unsigned bits = 1; bits <= 32; ++bits) {
        v.bits = (uint8_t)bits;
        for (unsigned phase = 0; phase < 2; ++phase) {
            v.phased = (uint8_t)phase;
            for (unsigned mode = 0; mode < 3; ++mode) {
                p = 2;
                if (mode == 0) {
                    v.type = BGEN_ENCODE_FLOAT64;
                    v.data = phase ? hp : gp;
                    v.width = phase ? 4 : 3;
                    v.data_bytes = v.width * sizeof(double);
                } else if (mode == 1) {
                    v.type = BGEN_ENCODE_FLOAT32;
                    v.data = phase ? hp32 : gp32;
                    v.width = phase ? 4 : 3;
                    v.data_bytes = v.width * sizeof(float);
                } else {
                    v.type = BGEN_ENCODE_INT64;
                    v.data = gt;
                    v.width = 63;
                    v.data_bytes = sizeof(gt);
                    p = 63;
                }
                for (missing = 0; missing < 2; ++missing) {
                    assert(bgen_layout2_size(&v, 1, sizeof(out), &size, &error) == 0);
                    memset(out, 0xcd, sizeof(out));
                    assert(bgen_layout2_encode(&v, 1, out, size, offsets, 2, &error) == 0);
                    assert(offsets[1] == size);
                    assert(out[size] == 0xcd);
                    if (missing) {
                        for (size_t j = 11; j < size; ++j) assert(out[j] == 0);
                    }
                }
            }
        }
    }
}

int main(void)
{
    double probabilities[3] = {0.5, 0.25, 0.25};
    uint8_t p = 2, missing = 0;
    struct bgen_layout2_variant v = {
        probabilities, sizeof(probabilities), 3, &p, 1, &missing, 1,
        1, 2, 0, 8, BGEN_ENCODE_FLOAT64
    };
    struct bgen_encode_error error;
    size_t size = 0;
    uint8_t out[13];
    uint64_t offsets[2];
    uint8_t const expected[13] = {1, 0, 0, 0, 2, 0, 2, 2, 2, 0, 8, 127, 64};
    assert(bgen_layout2_size(&v, 1, sizeof(out), &size, &error) == 0);
    assert(size == sizeof(out));
    assert(bgen_layout2_encode(&v, 1, out, sizeof(out), offsets, 2, &error) == 0);
    assert(memcmp(out, expected, sizeof(out)) == 0);
    assert(offsets[0] == 0 && offsets[1] == sizeof(out));
    assert(bgen_layout2_encode(&v, 1, out, sizeof(out) - 1, offsets, 2, &error) == -1);
    assert(bgen_layout2_encode(&v, 1, out, sizeof(out), offsets, 1, &error) == -1);
    assert(bgen_layout2_encode(&v, 1, NULL, sizeof(out), offsets, 2, &error) == -1);
    assert(bgen_layout2_size(NULL, 1, SIZE_MAX, &size, &error) == -1);
    assert(bgen_layout2_size(&v, 1, SIZE_MAX, NULL, &error) == -1);
    v.data_bytes--;
    assert(bgen_layout2_size(&v, 1, SIZE_MAX, &size, &error) == -1);
    v.data_bytes++;
    v.width = SIZE_MAX;
    assert(bgen_layout2_size(&v, 1, SIZE_MAX, &size, &error) == -1);
    v.width = 3;
    v.ploidy_bytes = 0;
    assert(bgen_layout2_size(&v, 1, SIZE_MAX, &size, &error) == -1);
    v.ploidy_bytes = 1;
    v.missing_bytes = 0;
    assert(bgen_layout2_size(&v, 1, SIZE_MAX, &size, &error) == -1);
    v.missing_bytes = 1;
    missing = 2;
    assert(bgen_layout2_size(&v, 1, SIZE_MAX, &size, &error) == -1);
    assert(error.variant == 0 && error.sample == 0);
    missing = 0;
    v.phased = 2;
    assert(bgen_layout2_size(&v, 1, SIZE_MAX, &size, &error) == -1);
    v.phased = 0;
    v.type = (enum bgen_encode_type)0;
    assert(bgen_layout2_size(&v, 1, SIZE_MAX, &size, &error) == -1);
    assert(error.variant == 0);
    v.type = BGEN_ENCODE_FLOAT64;
    probabilities[0] = -1;
    assert(bgen_layout2_encode(&v, 1, out, sizeof(out), offsets, 2, &error) == -1);
    assert(error.variant == 0 && error.sample == 0);
    /* No data is read by measurement, even for this deliberately huge row. */
    p = 63;
    v.nalleles = 8;
    v.bits = 32;
    v.width = UINT32_MAX;
    v.data_bytes = SIZE_MAX;
    assert(bgen_layout2_size(&v, 1, SIZE_MAX, &size, &error) == -1);
    assert(strstr(error.message, "uint32 length") != NULL);
    assert(bgen_layout2_encode(NULL, 0, NULL, 0, offsets, 1, &error) == 0);
    assert(offsets[0] == 0);
    sweep();
    puts("layout-2 encoder checks passed");
    return 0;
}
