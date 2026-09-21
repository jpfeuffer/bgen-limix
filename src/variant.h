#ifndef BGEN_VARIANT_H_PRIVATE
#define BGEN_VARIANT_H_PRIVATE

#include <inttypes.h>

struct bgen_file;
struct bgen_variant;

// bgen_variant_begin/next/destroy are public (bgen/variant.h); only the
// end-of-iteration sentinel and internal constructors stay private.
struct bgen_variant* bgen_variant_create(void);
void bgen_variant_create_alleles(struct bgen_variant* variant, uint16_t nalleles);
struct bgen_variant* bgen_variant_end(struct bgen_file const* bgen_file);

#endif
