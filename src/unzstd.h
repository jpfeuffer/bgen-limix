#ifndef BGEN_UNZSTD_H
#define BGEN_UNZSTD_H

#include <stddef.h>

/** Opaque persistent zstd decompressor context (Plan B). */
struct bgen_zstd_ctx_s;
typedef struct bgen_zstd_ctx_s bgen_zstd_ctx;

bgen_zstd_ctx* bgen_zstd_ctx_create(void);
void           bgen_zstd_ctx_destroy(bgen_zstd_ctx* ctx);

int bgen_unzstd(char const* src, size_t src_size, void** dst, size_t* dst_size);
int bgen_unzstd_reuse(bgen_zstd_ctx* ctx, char const* src, size_t src_size, void** dst,
                      size_t* dst_size);

#endif
