/* S3 integration test: open a BGEN file from an S3 URI and verify basic metadata.
 *
 * The file path is taken from the env var BGEN_TEST_S3_URI, which should be
 * set to an "s3://bucket/key" pointing to example.14bits.bgen.
 * If the env var is not set, the test is skipped (exits 0).
 *
 * Also requires the matching metafile to be passed as BGEN_TEST_S3_METAFILE
 * (a local path to the .metafile, which is not streamed).
 */
#include "bgen/bgen.h"
#include "cass.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    const char* uri = getenv("BGEN_TEST_S3_URI");
    if (!uri) {
        printf("BGEN_TEST_S3_URI not set — skipping S3 test.\n");
        return 0;
    }

    const char* metafile_path = getenv("BGEN_TEST_S3_METAFILE");
    if (!metafile_path) {
        printf("BGEN_TEST_S3_METAFILE not set — skipping S3 test.\n");
        return 0;
    }

    printf("Testing S3 open: %s\n", uri);

    /* ── Open bgen file from S3 ─── */
    struct bgen_file* bgen = bgen_file_open(uri);
    cass_cond(bgen != NULL);
    if (bgen == NULL)
        return cass_status();

    /* example.14bits.bgen has 500 samples and 199 variants */
    cass_equal_int(bgen_file_nsamples(bgen), 500);
    cass_equal_int(bgen_file_nvariants(bgen), 199);
    cass_cond(bgen_file_contain_samples(bgen));

    struct bgen_samples* samples = bgen_file_read_samples(bgen);
    cass_cond(samples != NULL);
    if (samples) {
        cass_cond(bgen_string_equal(BGEN_STRING("sample_001"),
                                    *bgen_samples_get(samples, 0)));
        cass_cond(bgen_string_equal(BGEN_STRING("sample_500"),
                                    *bgen_samples_get(samples, 499)));
        bgen_samples_destroy(samples);
    }

    /* ── Open metafile locally (index) ─── */
    struct bgen_metafile* metafile = bgen_metafile_open(metafile_path);
    cass_cond(metafile != NULL);
    if (metafile == NULL) {
        bgen_file_close(bgen);
        return cass_status();
    }

    cass_equal_int(bgen_metafile_nvariants(metafile), 199);

    /* Read first partition and open a genotype from S3 */
    struct bgen_partition const* part = bgen_metafile_read_partition(metafile, 0);
    cass_cond(part != NULL);
    if (part) {
        cass_cond(bgen_partition_nvariants(part) == 67);
        cass_cond(bgen_string_equal(BGEN_STRING("SNPID_2"),
                                    *bgen_partition_get_variant(part, 0)->id));

        /* Open and read probabilities for the first variant — exercises S3 seek+read */
        struct bgen_variant const* v = bgen_partition_get_variant(part, 0);
        struct bgen_genotype*      gt = bgen_file_open_genotype(bgen, v->genotype_offset);
        cass_cond(gt != NULL);
        if (gt) {
            cass_equal_int(bgen_genotype_ncombs(gt), 3);
            cass_equal_int(bgen_genotype_max_ploidy(gt), 2);

            uint32_t ncombs = bgen_genotype_ncombs(gt);
            double*  probs = calloc(500 * ncombs, sizeof(*probs));
            cass_cond(bgen_genotype_read(gt, probs) == 0);
            free(probs);
            bgen_genotype_close(gt);
        }
        bgen_partition_destroy(part);
    }

    bgen_metafile_close(metafile);
    bgen_file_close(bgen);

    if (cass_status() == 0)
        printf("All S3 checks passed.\n");

    return cass_status();
}
