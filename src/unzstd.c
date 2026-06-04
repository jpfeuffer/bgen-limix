#include "unzstd.h"
#include "report.h"
#include "zstd.h"
#include <stdlib.h>

struct bgen_zstd_ctx_s
{
    ZSTD_DCtx* dctx;
};

bgen_zstd_ctx* bgen_zstd_ctx_create(void)
{
    bgen_zstd_ctx* ctx = malloc(sizeof(bgen_zstd_ctx));
    if (!ctx) return NULL;
    ctx->dctx = ZSTD_createDCtx();
    if (!ctx->dctx) {
        free(ctx);
        return NULL;
    }
    return ctx;
}

void bgen_zstd_ctx_destroy(bgen_zstd_ctx* ctx)
{
    if (!ctx) return;
    ZSTD_freeDCtx(ctx->dctx);
    free(ctx);
}

/* Reuse an existing context — avoids ZSTD_createDCtx/freeDCtx on every call. */
int bgen_unzstd_reuse(bgen_zstd_ctx* ctx, char const* src, size_t src_size, void** dst,
                      size_t* dst_size)
{
    if (!ctx) return bgen_unzstd(src, src_size, dst, dst_size);

    size_t dSize = ZSTD_decompressDCtx(ctx->dctx, *dst, *dst_size, src, src_size);
    if (ZSTD_isError(dSize)) {
        bgen_error("zstd decoding (%s)", ZSTD_getErrorName(dSize));
        return 1;
    }
    return 0;
}

int bgen_unzstd(char const* src, size_t src_size, void** dst, size_t* dst_size)
{
    size_t dSize = ZSTD_decompress(*dst, *dst_size, src, src_size);

    if (ZSTD_isError(dSize)) {
        bgen_error("zstd decoding (%s)", ZSTD_getErrorName(dSize));
        return 1;
    }

    return 0;
}
