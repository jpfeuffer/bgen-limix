#ifndef BGEN_UNZLIB_H
#define BGEN_UNZLIB_H

#include <stddef.h>

/** Opaque persistent zlib-ng decompressor context (Plan B). */
struct bgen_zlib_ctx_s;
typedef struct bgen_zlib_ctx_s bgen_zlib_ctx;

bgen_zlib_ctx* bgen_zlib_ctx_create(void);
void           bgen_zlib_ctx_destroy(bgen_zlib_ctx* ctx);

int bgen_unzlib(char const* src, size_t src_size, char** dst, size_t* dst_size);
int bgen_unzlib_reuse(bgen_zlib_ctx* ctx, char const* src, size_t src_size, char** dst,
                      size_t* dst_size);
int bgen_unzlib_chunked(char const* src, size_t src_size, char** dst, size_t* dst_size);

#endif
