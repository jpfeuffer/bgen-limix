#include "bgen/bgen.h"
#include "cass.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
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

static void check_batch(struct bgen_file* file, uint64_t const* offsets, unsigned width)
{
    uint64_t selected[3] = {offsets[1], offsets[0], offsets[1]};
    bool phased[3], missing[6];
    uint8_t ploidy[6], bits[3];
    uint16_t alleles[3];
    uint32_t counts[6], widths[3];
    struct bgen_batch_metadata metadata = {
        phased, ploidy, missing, widths, counts, alleles, bits
    };
    unsigned padded = width + 2;
    double probabilities[48], exact[36];
    float probabilities32[48];
    cass_cond(bgen_file_read_genotypes_batch64(file, selected, 3, padded,
                                              probabilities, &metadata) == 0);
    cass_cond(bgen_file_read_genotypes_batch32(file, selected, 3, padded,
                                              probabilities32, NULL) == 0);
    cass_cond(bgen_file_read_genotypes_batch(file, selected, 3, width,
                                            width ? exact : NULL) == 0);
    double compact[48];
    float compact32[48];
    cass_cond(bgen_file_read_genotypes_batch_compact64(file, selected, 3, padded,
                                                       compact, NULL) == 0);
    cass_cond(bgen_file_read_genotypes_batch_compact32(file, selected, 3, padded,
                                                       compact32, NULL) == 0);
    for (unsigned c = 0; c < 6 * padded; ++c) {
        cass_cond(isnan(probabilities[c]) ? isnan(compact[c]) : probabilities[c] == compact[c]);
        cass_cond(isnan(probabilities32[c]) ? isnan(compact32[c]) :
            probabilities32[c] == compact32[c]);
    }
    struct bgen_dimensions dimensions = {0};
    for (unsigned actual = 0; actual < 2; ++actual) {
        cass_cond(bgen_file_read_dimensions(file, selected, 3, actual != 0,
                                             BGEN_MODEL_AUTO, &dimensions) == 0);
        cass_cond(dimensions.gp_width == (phased[0] ? 0 : width));
        cass_cond(dimensions.hp_width == (phased[0] ? (width ? width : 1) : 0));
        cass_cond(dimensions.decode_width == (width ? width : 1));
    }
    cass_cond(bgen_file_read_dimensions(file, NULL, 0, true,
                                         BGEN_MODEL_AUTO, &dimensions) == 0);
    cass_cond(dimensions.gp_width == 0 && dimensions.hp_width == 0 &&
              dimensions.decode_width == 0);
    dimensions.gp_width = 42;
    cass_cond(bgen_file_read_dimensions(NULL, selected, 1, true,
                                         BGEN_MODEL_AUTO, &dimensions) != 0);
    cass_cond(dimensions.gp_width == 42);
    cass_cond(bgen_file_read_dimensions(file, NULL, 1, true,
                                         BGEN_MODEL_AUTO, &dimensions) != 0);
    cass_cond(bgen_file_read_dimensions(file, selected, 1, true,
                                         BGEN_MODEL_AUTO, NULL) != 0);
    cass_cond(bgen_file_read_dimensions(file, selected, 1, true,
                                         (enum bgen_probability_model)99, &dimensions) != 0);
    for (unsigned i = 0; i < 3; ++i) {
        struct bgen_genotype* gt = bgen_file_open_genotype(file, selected[i]);
        cass_cond(gt != NULL);
        if (!gt)
            return;
        double single[12];
        float single32[12];
        cass_cond(bgen_genotype_read64(gt, width ? single : NULL) == 0);
        cass_cond(bgen_genotype_read32(gt, width ? single32 : NULL) == 0);
        cass_cond(widths[i] == width);
        cass_cond(phased[i] == bgen_genotype_phased(gt));
        cass_cond(bits[i] == bgen_genotype_nbits(gt));
        cass_cond(alleles[i] == (bgen_file_layout(file) == 1 ? 2 : 3));
        for (unsigned s = 0; s < 2; ++s) {
            cass_cond(ploidy[i * 2 + s] == bgen_genotype_ploidy(gt, s));
            cass_cond(missing[i * 2 + s] == bgen_genotype_missing(gt, s));
            unsigned p = ploidy[i * 2 + s];
            unsigned expected = phased[i] ? p * alleles[i] :
                (alleles[i] == 2 ? p + 1 : (p + 1) * (p + 2) / 2);
            cass_cond(counts[i * 2 + s] == expected);
            for (unsigned c = 0; c < padded; ++c) {
                size_t index = (i * 2 + s) * padded + c;
                if (c >= width) {
                    cass_cond(isnan(probabilities[index]));
                    cass_cond(isnan(probabilities32[index]));
                } else {
                    double want = single[s * width + c];
                    float want32 = single32[s * width + c];
                    cass_cond(isnan(want) ? isnan(probabilities[index]) :
                        probabilities[index] == want);
                    cass_cond(isnan(want32) ? isnan(probabilities32[index]) :
                        probabilities32[index] == want32);
                    cass_cond(isnan(want) ? isnan(exact[(i * 2 + s) * width + c]) :
                        exact[(i * 2 + s) * width + c] == want);
                }
            }
        }
        bgen_genotype_close(gt);
    }
    memset(widths, 0xff, sizeof(widths));
    cass_cond(bgen_file_read_genotype_metadata(file, selected, 3, &metadata) == 0);
    for (unsigned i = 0; i < 3; ++i)
        cass_cond(widths[i] == width);
    cass_cond(bgen_file_read_genotype_metadata(file, selected, 3, NULL) == 0);
    cass_cond(bgen_file_read_ncombs_batch(file, selected, 3, widths) == 0);
    cass_cond(bgen_file_read_genotypes_batch_padded(file, selected, 3, padded,
                                                   probabilities) == 0);
    cass_cond(bgen_file_read_genotypes_batch64(file, NULL, 0, UINT32_MAX, NULL, NULL) == 0);
    cass_cond(bgen_file_read_genotype_metadata(file, NULL, 0, NULL) == 0);
    cass_cond(bgen_file_read_genotypes_batch64(NULL, selected, 1, width, NULL, NULL) != 0);
    cass_cond(bgen_file_read_genotype_metadata(file, NULL, 1, NULL) != 0);
    cass_cond(bgen_file_read_ncombs_batch(file, selected, 1, NULL) != 0);
    cass_cond(bgen_file_read_genotypes_batch64(file, selected, UINT32_MAX,
                                              UINT32_MAX, probabilities, NULL) != 0);
    if (width) {
        probabilities[0] = 42;
        cass_cond(bgen_file_read_genotypes_batch64(file, selected, 1, width - 1,
                                                  probabilities, NULL) != 0);
        cass_cond(probabilities[0] == 42);
        cass_cond(bgen_file_read_genotypes_batch_compact64(file, selected, 1, width - 1,
                                                            probabilities, NULL) != 0);
        cass_cond(probabilities[0] == 42);
        cass_cond(bgen_file_read_genotypes_batch64(file, selected, 1, width, NULL, NULL) != 0);
    }
    uint64_t invalid[2] = {offsets[0], UINT64_MAX};
    cass_cond(bgen_file_read_genotype_metadata(file, invalid, 2, &metadata) != 0);
    cass_cond(bgen_file_read_genotypes_batch32(file, invalid, 2, padded,
                                              probabilities32, NULL) != 0);
    cass_cond(bgen_file_read_genotypes_batch64(file, selected, 1, width,
                                              width ? probabilities : NULL, &metadata) == 0);
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
    check_batch(file, offsets, 3);
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
    check_batch(file, offsets, zero_ploidy ? 0 : 6);
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

static void test_compact_declared_overflow(void)
{
    char const* path = "native_layouts.tmp/compact.bgen";
    FILE* output = write_header(path, 2);
    if (!output) return;
    uint64_t offsets[2];
    for (unsigned record = 0; record < 2; ++record) {
        offsets[record] = write_variant(output, 2, 16);
        unsigned payload = record ? 30 : 0;
        write_uint(output, 12 + payload, 4);
        write_uint(output, 2, 4);
        write_uint(output, 16, 2);
        write_uint(output, 0, 1);
        write_uint(output, 63, 1);
        write_uint(output, record, 1);
        write_uint(output, record | 128, 1);
        write_uint(output, 0, 1);
        write_uint(output, 8, 1);
        for (unsigned c = 0; c < payload; ++c)
            write_uint(output, 0, 1);
    }
    cass_cond(fclose(output) == 0);
    struct bgen_file* file = bgen_file_open(path);
    cass_cond(file != NULL);
    if (!file) return;
    struct bgen_dimensions dimensions = {0};
    cass_cond(bgen_file_read_dimensions(file, offsets, 2, true,
                                         BGEN_MODEL_AUTO, &dimensions) == 0);
    cass_cond(dimensions.gp_width == 16 && dimensions.decode_width == 16);
    dimensions.gp_width = 42;
    cass_cond(bgen_file_read_dimensions(file, offsets, 2, false,
                                         BGEN_MODEL_AUTO, &dimensions) != 0);
    cass_cond(dimensions.gp_width == 42);
    cass_cond(bgen_file_read_dimensions_with_ploidy(file, offsets, 2, false,
                                                     BGEN_MODEL_AUTO, 0, &dimensions) == 0);
    cass_cond(dimensions.gp_width == 1);
    double probs[4] = {42, 42, 42, 42};
    float probs32[4] = {42, 42, 42, 42};
    uint32_t widths[2] = {42, 42};
    struct bgen_batch_metadata metadata = {0};
    metadata.ncombs = widths;
    cass_cond(bgen_file_read_genotypes_batch_compact64(file, offsets, 2, 1,
                                                       probs, &metadata) != 0);
    cass_cond(probs[0] == 1 && isnan(probs[1]));
    cass_cond(probs[2] == 42 && probs[3] == 42);
    cass_cond(widths[0] == 1 && widths[1] == 42);
    cass_cond(bgen_file_read_genotypes_batch_compact32(file, offsets, 2, 1,
                                                       probs32, NULL) != 0);
    cass_cond(probs32[0] == 1 && isnan(probs32[1]));
    cass_cond(probs32[2] == 42 && probs32[3] == 42);
    cass_cond(bgen_file_read_genotypes_batch_compact64(file, offsets, 1, 1, probs, NULL) == 0);
    cass_cond(bgen_file_read_dimensions_with_ploidy(file, offsets, 1, true,
                                                     BGEN_MODEL_AUTO, 64, &dimensions) != 0);
    bgen_file_close(file);
}

static void check_compressed_dimensions(char const* path)
{
    struct bgen_file* file = bgen_file_open(path);
    cass_cond(file != NULL);
    if (!file) return;
    struct bgen_metafile* index =
        bgen_metafile_create(file, "native_layouts.tmp/dimensions.index", 1, 0);
    cass_cond(index != NULL);
    if (!index) { bgen_file_close(file); return; }
    struct bgen_partition const* part = bgen_metafile_read_partition(index, 0);
    cass_cond(part != NULL);
    if (!part) { bgen_metafile_close(index); bgen_file_close(file); return; }
    uint32_t n = bgen_partition_nvariants(part);
    uint64_t* offsets = malloc((size_t)n * sizeof(uint64_t));
    cass_cond(offsets != NULL);
    if (!offsets) goto cleanup;
    uint32_t declared_gp = 0, declared_hp = 0;
    for (uint32_t i = 0; i < n; ++i) {
        offsets[i] = bgen_partition_get_variant(part, i)->genotype_offset;
        struct bgen_genotype* gt = bgen_file_open_genotype(file, offsets[i]);
        cass_cond(gt != NULL);
        if (!gt) continue;
        uint32_t width = bgen_genotype_ncombs(gt);
        if (bgen_genotype_phased(gt)) {
            if (!width) width = 1;
            if (width > declared_hp) declared_hp = width;
        } else if (width > declared_gp) declared_gp = width;
        bgen_genotype_close(gt);
    }
    struct bgen_dimensions declared, actual;
    cass_cond(bgen_file_read_dimensions(file, offsets, n, false,
                                         BGEN_MODEL_AUTO, &declared) == 0);
    cass_cond(declared.gp_width == declared_gp && declared.hp_width == declared_hp);
    cass_cond(bgen_file_read_dimensions(file, offsets, n, true,
                                         BGEN_MODEL_AUTO, &actual) == 0);
    cass_cond(actual.gp_width <= declared.gp_width && actual.hp_width <= declared.hp_width);
    for (int model = BGEN_MODEL_AUTO; model <= BGEN_MODEL_MIXED; ++model) {
        struct bgen_dimensions assumed;
        cass_cond(bgen_file_read_dimensions_with_ploidy(file, offsets, n, false,
                   (enum bgen_probability_model)model, 0, &assumed) == 0);
        cass_cond(assumed.gp_width <= 1 && assumed.hp_width <= 1);
        cass_cond(assumed.decode_width == (n ? 1 : 0));
    }
    free(offsets);
cleanup:
    bgen_partition_destroy(part);
    bgen_metafile_close(index);
    bgen_file_close(file);
}

int main(void)
{
    test_layout1();
    for (unsigned bits = 1; bits <= 32; ++bits)
        for (unsigned phased = 0; phased < 2; ++phased)
            test_layout2(bits, phased, 0);
    test_layout2(8, 1, 1);
    test_compact_declared_overflow();
    check_compressed_dimensions(TEST_DATADIR "example.1bits.bgen");
    check_compressed_dimensions(TEST_DATADIR "example.v11.bgen");
    check_compressed_dimensions(TEST_DATADIR "haplotypes.bgen");
    check_compressed_dimensions(TEST_DATADIR "complex.23bits.bgen");
    return cass_status();
}
