#include "bgen/bgen.h"
#include "cass.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static void write_uint(FILE* file, uint32_t value, unsigned width)
{
    for (unsigned i = 0; i < width; ++i)
        cass_cond(fputc((int)((value >> (8 * i)) & 255), file) != EOF);
}

static void write_string(FILE* file, char const* value, unsigned width)
{
    size_t size = strlen(value);
    write_uint(file, (uint32_t)size, width);
    cass_cond(fwrite(value, 1, size, file) == size);
}

static FILE* write_header(char const* path, unsigned layout)
{
    FILE* file = fopen(path, "wb");
    if (!file) {
        cass_cond(file != NULL);
        return NULL;
    }
    write_uint(file, 20, 4);
    write_uint(file, 20, 4);
    write_uint(file, 2, 4);
    write_uint(file, 2, 4);
    cass_cond(fwrite("bgen", 1, 4, file) == 4);
    write_uint(file, layout << 2, 4);
    return file;
}

static uint64_t write_variant(FILE* file, unsigned layout, unsigned nalleles)
{
    if (layout == 1)
        write_uint(file, 2, 4);
    write_string(file, "snp", 2);
    write_string(file, "rs", 2);
    write_string(file, "chr1", 2);
    write_uint(file, 42, 4);
    if (layout == 2)
        write_uint(file, nalleles, 2);
    for (unsigned i = 0; i < nalleles; ++i)
        write_string(file, i == 0 ? "A" : "T", 4);
    return (uint64_t)ftell(file);
}

static void check_index(struct bgen_file* file, uint64_t const* offsets)
{
    struct bgen_metafile* mf =
        bgen_metafile_create(file, "native_layouts.tmp/index", 1, 0);
    cass_cond(mf != NULL);
    if (!mf)
        return;
    struct bgen_partition const* part = bgen_metafile_read_partition(mf, 0);
    cass_cond(part != NULL);
    if (part) {
        cass_cond(bgen_partition_nvariants(part) == 2);
        for (uint32_t i = 0; i < 2; ++i) {
            struct bgen_variant const* variant = bgen_partition_get_variant(part, i);
            cass_cond(variant->genotype_offset == offsets[i]);
            cass_cond(variant->position == 42);
            cass_cond(bgen_string_equal(BGEN_STRING("snp"), *variant->id));
        }
        bgen_partition_destroy(part);
    }
    cass_cond(bgen_metafile_close(mf) == 0);
}

static void test_layout1(void)
{
    char const* path = "native_layouts.tmp/layout1.bgen";
    FILE* output = write_header(path, 1);
    if (!output)
        return;
    uint64_t offsets[2];
    for (unsigned i = 0; i < 2; ++i) {
        offsets[i] = write_variant(output, 1, 2);
        write_uint(output, 6553, 2);
        write_uint(output, 10922, 2);
        write_uint(output, 15292, 2);
        for (unsigned j = 0; j < 3; ++j)
            write_uint(output, 0, 2);
    }
    cass_cond(fclose(output) == 0);
    struct bgen_file* file = bgen_file_open(path);
    cass_cond(file != NULL);
    if (!file)
        return;
    check_index(file, offsets);
    struct bgen_genotype* first = bgen_file_open_genotype(file, offsets[0]);
    struct bgen_genotype* second = bgen_file_open_genotype(file, offsets[1]);
    cass_cond(first != NULL && second != NULL);
    bgen_file_close(file);
    if (second)
        bgen_genotype_close(second);
    if (!first)
        return;
    double probabilities[6];
    float probabilities32[6];
    cass_cond(bgen_genotype_read64(first, probabilities) == 0);
    cass_cond(bgen_genotype_read32(first, probabilities32) == 0);
    cass_cond(!bgen_genotype_phased(first));
    cass_cond(bgen_genotype_nbits(first) == 16);
    cass_cond(bgen_genotype_min_ploidy(first) == 2);
    cass_cond(bgen_genotype_max_ploidy(first) == 2);
    for (uint32_t i = 0; i < 2; ++i) {
        cass_cond(bgen_genotype_ploidy(first, i) == 2);
        cass_cond(bgen_genotype_missing(first, i) == (i == 1));
    }
    cass_cond(probabilities[0] == 6553.0 / 32768);
    cass_cond(probabilities[1] == 10922.0 / 32768);
    cass_cond(probabilities[2] == 15292.0 / 32768);
    cass_cond(probabilities[0] + probabilities[1] + probabilities[2] == 32767.0 / 32768);
    for (unsigned i = 0; i < 3; ++i) {
        cass_cond(probabilities32[i] == (float)probabilities[i]);
        cass_cond(isnan(probabilities[i + 3]));
        cass_cond(isnan(probabilities32[i + 3]));
    }
    bgen_genotype_close(first);
}

static void test_layout2(unsigned bits, unsigned phased, unsigned zero_ploidy)
{
    char const* path = "native_layouts.tmp/layout2.bgen";
    FILE* output = write_header(path, 2);
    if (!output)
        return;
    uint64_t offsets[2];
    uint32_t denominator = (uint32_t)((UINT64_C(1) << bits) - 1);
    unsigned count = phased ? 3 : 6;
    uint32_t value = denominator / count;
    unsigned stored = zero_ploidy ? 0 : (phased ? 4 : 5);
    unsigned payload_size = (stored * bits + 7) / 8;
    for (unsigned record = 0; record < 2; ++record) {
        offsets[record] = write_variant(output, 2, 3);
        write_uint(output, 12 + payload_size, 4);
        write_uint(output, 2, 4);
        write_uint(output, 3, 2);
        write_uint(output, 0, 1);
        write_uint(output, zero_ploidy ? 0 : 2, 1);
        write_uint(output, zero_ploidy ? 128 : 0, 1);
        write_uint(output, zero_ploidy ? 0 : 2, 1);
        write_uint(output, phased, 1);
        write_uint(output, bits, 1);
        uint8_t packed[32] = {0};
        for (unsigned i = 0; i < stored; ++i)
            for (unsigned bit = 0; bit < bits; ++bit)
                if ((value >> bit) & 1)
                    packed[(i * bits + bit) / 8] |= (uint8_t)(1U << ((i * bits + bit) % 8));
        cass_cond(fwrite(packed, 1, payload_size, output) == payload_size);
    }
    cass_cond(fclose(output) == 0);
    struct bgen_file* file = bgen_file_open(path);
    cass_cond(file != NULL);
    if (!file)
        return;
    check_index(file, offsets);
    for (unsigned record = 0; record < 2; ++record) {
        struct bgen_genotype* gt = bgen_file_open_genotype(file, offsets[record]);
        cass_cond(gt != NULL);
        if (!gt)
            break;
        unsigned ncombs = zero_ploidy ? 0 : 6;
        cass_cond(bgen_genotype_ncombs(gt) == ncombs);
        cass_cond(bgen_genotype_phased(gt) == (phased != 0));
        cass_cond(bgen_genotype_ploidy(gt, 0) == 0);
        cass_cond(bgen_genotype_ploidy(gt, 1) == (zero_ploidy ? 0 : 2));
        cass_cond(bgen_genotype_missing(gt, 0) == (zero_ploidy != 0));
        double probabilities[12];
        float probabilities32[12];
        cass_cond(bgen_genotype_read64(gt, zero_ploidy ? NULL : probabilities) == 0);
        cass_cond(bgen_genotype_read32(gt, zero_ploidy ? NULL : probabilities32) == 0);
        if (!zero_ploidy) {
            for (unsigned i = 0; i < 6; ++i) {
                double expected = (i % count == count - 1) ?
                    (double)(denominator - (count - 1) * value) / denominator :
                    (double)value / denominator;
                cass_cond(probabilities[6 + i] == expected);
                cass_cond(fabs((double)probabilities32[6 + i] - expected) < 1e-7);
            }
            cass_cond(phased ? isnan(probabilities[0]) : probabilities[0] == 1);
        }
        bgen_genotype_close(gt);
    }
    bgen_file_close(file);
}

int main(void)
{
    test_layout1();
    for (unsigned bits = 1; bits <= 32; ++bits)
        for (unsigned phased = 0; phased < 2; ++phased)
            test_layout2(bits, phased, 0);
    test_layout2(8, 1, 1);
    return cass_status();
}
