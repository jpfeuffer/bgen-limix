#ifndef BGEN_GENOTYPE_H_PRIVATE
#define BGEN_GENOTYPE_H_PRIVATE

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

struct bgen_genotype
{
    unsigned    layout;
    uint32_t    nsamples;
    uint16_t    nalleles;
    uint8_t     phased;
    uint8_t     nbits;
    uint8_t*    ploidy_missingness;
    unsigned    ncombs;
    uint8_t     min_ploidy;
    uint8_t     max_ploidy;
    char*       chunk;
    char const* chunk_ptr;
    uint64_t    offset;
    /* ownership flags: false when buffer points into bgen_file scratch (Plan A) */
    bool        owns_chunk;   /**< if false, chunk is owned by bgen_file; do not free */
    bool        owns_ploidy;  /**< if false, ploidy_missingness is owned by bgen_file */
};

static inline struct bgen_genotype* bgen_genotype_create(void)
{
    struct bgen_genotype* genotype = malloc(sizeof(struct bgen_genotype));
    genotype->layout = 0;
    genotype->nsamples = 0;
    genotype->nalleles = 0;
    genotype->phased = 0;
    genotype->nbits = 0;
    genotype->ploidy_missingness = NULL;
    genotype->ncombs = 0;
    genotype->min_ploidy = 0;
    genotype->max_ploidy = 0;
    genotype->chunk = NULL;
    genotype->chunk_ptr = NULL;
    genotype->offset = 0;
    genotype->owns_chunk  = true;
    genotype->owns_ploidy = true;
    return genotype;
}

/** Initialise a stack-allocated genotype (no heap allocation, no free needed). */
static inline void bgen_genotype_init(struct bgen_genotype* genotype)
{
    genotype->layout = 0;
    genotype->nsamples = 0;
    genotype->nalleles = 0;
    genotype->phased = 0;
    genotype->nbits = 0;
    genotype->ploidy_missingness = NULL;
    genotype->ncombs = 0;
    genotype->min_ploidy = 0;
    genotype->max_ploidy = 0;
    genotype->chunk = NULL;
    genotype->chunk_ptr = NULL;
    genotype->offset = 0;
    genotype->owns_chunk  = false;
    genotype->owns_ploidy = false;
}

#endif
