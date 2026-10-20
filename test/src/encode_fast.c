#include "bgen/bgen.h"
#undef NDEBUG
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint32_t le32(uint8_t const *p)
{
    return (uint32_t)p[0] | (uint32_t)p[1] << 8 |
        (uint32_t)p[2] << 16 | (uint32_t)p[3] << 24;
}

static void legacy_checks(void)
{
    int8_t gt[4] = {0, 1, -1, -2};
    uint8_t phase = 0, out[16];
    uint32_t length;
    struct bgen_encode_error e;
    uint8_t const expected[] = {2, 0, 0, 0, 2, 0, 1, 2, 2, 129, 0, 8, 0, 255, 0};
    memset(out, 0xCC, sizeof(out));
    assert(bgen_layout2_hardcalls_legacy(gt, 4, &phase, 1, 1, 2, out, 16, 16,
                                       &length, 1, &e) == 0);
    assert(length == 15 && memcmp(out, expected, 15) == 0 && out[15] == 0xCC);
    assert(bgen_layout2_hardcalls_legacy(gt, 3, &phase, 1, 1, 2, out, 16, 16,
                                       &length, 1, &e) == -1);
    assert(bgen_layout2_hardcalls_legacy(gt, 4, &phase, 0, 1, 2, out, 16, 16,
                                       &length, 1, &e) == -1);
    assert(bgen_layout2_hardcalls_legacy(gt, 4, &phase, 1, 1, 2, out, 16, 15,
                                       &length, 1, &e) == BGEN_ENCODE_BUFFER_TOO_SMALL);
    assert(bgen_layout2_hardcalls_legacy(gt, 4, &phase, 1, 1, 2, out, 15, 16,
                                       &length, 1, &e) == BGEN_ENCODE_BUFFER_TOO_SMALL);
    assert(bgen_layout2_hardcalls_legacy(gt, 4, &phase, 1, 1, 2, out, 16, 16,
                                       &length, 0, &e) == -1);
    assert(bgen_layout2_hardcalls_legacy(gt, 4, &phase, 1, 1, 2, NULL, 16, 16,
                                       &length, 1, &e) == -1);
    for (int invalid = -128; invalid < 128; ++invalid) {
        if (invalid >= -2 && invalid <= 1) continue;
        gt[2] = (int8_t)invalid;
        assert(bgen_layout2_hardcalls_legacy(gt, 4, &phase, 1, 1, 2, out, 16, 16,
                                           &length, 1, &e) == BGEN_ENCODE_INVALID_ALLELE);
        assert(e.variant == 0 && e.sample == 1);
    }
    gt[2] = -2;
    assert(bgen_layout2_hardcalls_legacy(gt, 4, &phase, 1, 1, 2, out, 16, 16,
                                       &length, 1, &e) == BGEN_ENCODE_INVALID_PLOIDY);
    gt[2] = -1;
    phase = 2;
    assert(bgen_layout2_hardcalls_legacy(gt, 4, &phase, 1, 1, 2, out, 16, 16,
                                       &length, 1, &e) == -1);
    assert(bgen_layout2_hardcalls_legacy(NULL, 0, NULL, 0, 0, 0, NULL, 0, 10,
                                       NULL, 0, &e) == 0);
    assert(bgen_layout2_hardcalls_legacy(NULL, 0, &phase, 1, SIZE_MAX, UINT32_MAX,
                                       NULL, 0, SIZE_MAX, NULL, 0, &e) == -1);
}

static void check_stored(uint8_t const *encoded, uint8_t const *raw, size_t bytes)
{
    assert(encoded[0] == 0x78 && encoded[1] == 1);
    size_t offset = 0, at = 2;
    uint32_t a = 1, b = 0;
    for (;;) {
        uint8_t final = encoded[at++];
        size_t n = (size_t)encoded[at] | (size_t)encoded[at + 1] << 8;
        unsigned inverse = (unsigned)encoded[at + 2] | (unsigned)encoded[at + 3] << 8;
        assert((n ^ inverse) == 65535 && n <= bytes - offset);
        at += 4;
        assert(memcmp(encoded + at, raw + offset, n) == 0);
        for (size_t j = 0; j < n; ++j) {
            a = (a + raw[offset + j]) % 65521;
            b = (b + a) % 65521;
        }
        offset += n;
        at += n;
        if (final) break;
        assert(n == 65535);
    }
    assert(offset == bytes);
    uint32_t checksum = (uint32_t)encoded[at] << 24 | (uint32_t)encoded[at + 1] << 16 |
        (uint32_t)encoded[at + 2] << 8 | encoded[at + 3];
    assert(checksum == (b << 16 | a));
}

static void fixed_checks(uint32_t ns, uint8_t ploidy, uint8_t phased)
{
    struct bgen_encode_error e;
    size_t row;
    assert(bgen_layout2_fixed_size(ns, ploidy, 5, &row, &e) == 0);
    int8_t *gt = malloc(ns ? (size_t)ns * 2 : 1);
    uint8_t *out = malloc(row + 1);
    size_t stride = 10 + (size_t)ns * 3;
    uint8_t *raw = malloc(stride);
    assert(gt && out && raw);
    for (size_t s = 0; s < ns; ++s) {
        gt[s * 2] = s % 5 == 0 ? -1 : (int8_t)(s % 2);
        gt[s * 2 + 1] = ploidy == 1 ? -2 : (int8_t)((s + 1) % 2);
    }
    uint8_t str = '.';
    int32_t pos = INT32_MAX;
    struct bgen_encode_strings string = {&str, 1, 1};
    struct bgen_layout2_fixed_batch batch = {
        string, string, string, string, string, &pos, 4, gt, (size_t)ns * 2,
        &phased, 1, 1, ns, ploidy, 5
    };
    memset(out, 0xCC, row + 1);
    assert(bgen_layout2_fixed_encode(&batch, out, row, &e) == 0);
    assert(out[row] == 0xCC);
    uint32_t raw_len;
    assert(bgen_layout2_hardcalls_legacy(gt, (size_t)ns * 2, &phased, 1, 1, ns,
                                       raw, stride, stride, &raw_len, 1, &e) == 0);
    assert(le32(out + 29) == raw_len);
    assert(le32(out + 25) == row - 29);
    check_stored(out + 33, raw, raw_len);
    assert(bgen_layout2_fixed_encode(&batch, out, row - 1, &e) == BGEN_ENCODE_BUFFER_TOO_SMALL);
    batch.varid.bytes = 0;
    assert(bgen_layout2_fixed_encode(&batch, out, row, &e) == -1);
    batch.varid.bytes = 1;
    batch.total_string_length = 4;
    assert(bgen_layout2_fixed_encode(&batch, out, row, &e) == BGEN_ENCODE_STRING_LENGTH);
    batch.total_string_length = 5;
    pos = -1;
    assert(bgen_layout2_fixed_encode(&batch, out, row, &e) == -1);
    pos = 1;
    if (ns) {
        gt[0] = 0;
        gt[1] = ploidy == 1 ? 0 : -2;
        assert(bgen_layout2_fixed_encode(&batch, out, row, &e) == BGEN_ENCODE_MIXED_PLOIDY);
        gt[0] = -2;
        assert(bgen_layout2_fixed_encode(&batch, out, row, &e) == BGEN_ENCODE_INVALID_PLOIDY);
        gt[0] = 3;
        assert(bgen_layout2_fixed_encode(&batch, out, row, &e) == BGEN_ENCODE_INVALID_ALLELE);
    }
    free(gt);
    free(out);
    free(raw);
}

int main(void)
{
    legacy_checks();
    uint32_t const samples[] = {0, 1, 21841, 21842, 32763, 43688};
    for (size_t j = 0; j < sizeof(samples) / sizeof(*samples); ++j) {
        for (uint8_t p = 1; p <= 2; ++p) {
            for (uint8_t phase = 0; phase < 2; ++phase) fixed_checks(samples[j], p, phase);
        }
    }
    struct bgen_encode_error e;
    size_t size;
    assert(bgen_layout2_fixed_size(UINT32_MAX, 2, 0, &size, &e) == -1);
    assert(bgen_layout2_fixed_size(1, 2, SIZE_MAX, &size, &e) == -1);
    assert(bgen_layout2_fixed_size(1, 0, 0, &size, &e) == -1);
    assert(bgen_layout2_fixed_encode(NULL, NULL, 0, &e) == -1);
    puts("specialized layout-2 encoding checks passed");
    return 0;
}
