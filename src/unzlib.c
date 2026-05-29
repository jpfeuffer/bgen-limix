#include "unzlib.h"
#include "report.h"
#include "zlib-ng.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

int bgen_unzlib(char const* src, size_t src_size, char** dst, size_t* dst_size)
{
    zng_stream strm = {
        .zalloc = Z_NULL,
        .zfree = Z_NULL,
        .opaque = Z_NULL,
        .avail_in = 0,
        .next_in = (unsigned char const*)src,
    };

    int e = zng_inflateInit(&strm);

    if (e != Z_OK) {
        bgen_error("zlib-ng failed to init (%s)", zng_zError(e));
        goto err;
    }

    if (src_size > UINT_MAX) {
        bgen_error("zlib-ng src_size overflow");
        goto err;
    }
    strm.avail_in = (unsigned)src_size;

    if (*dst_size > UINT_MAX) {
        bgen_error("zlib-ng *dst_size overflow");
        goto err;
    }
    strm.avail_out = (unsigned)*dst_size;
    strm.next_out = (unsigned char*)*dst;

    e = zng_inflate(&strm, Z_FINISH);
    if (e != Z_STREAM_END) {
        bgen_error("zlib-ng failed to inflate (%s)", zng_zError(e));
        goto err;
    }

    if (zng_inflateEnd(&strm) != Z_OK) {
        bgen_error("zlib-ng failed to inflateEnd (%s)", zng_zError(e));
        return 1;
    }
    return 0;

err:
    zng_inflateEnd(&strm);
    return 1;
}

int bgen_unzlib_chunked(char const* src, size_t src_size, char** dst, size_t* dst_size)
{
    if (*dst_size > UINT_MAX) {
        bgen_error("zlib-ng *dst_size overflow");
        return 1;
    }

    unsigned       unused = (unsigned)*dst_size;
    unsigned char* cdst = (unsigned char*)*dst;

    zng_stream strm = {
        .zalloc = Z_NULL, .zfree = Z_NULL, .opaque = Z_NULL, .avail_in = 0, .next_in = Z_NULL};

    int ret = zng_inflateInit(&strm);

    if (ret != Z_OK) {
        bgen_error("zlib-ng failed to uncompress (%s)", zng_zError(ret));
        goto err;
    }

    strm.avail_in = (unsigned)src_size;
    strm.next_in = (unsigned char const*)src;

    while (1) {
        strm.avail_out = unused;
        strm.next_out = cdst;

        ret = zng_inflate(&strm, Z_NO_FLUSH);

        if (ret == Z_NEED_DICT) {
            ret = Z_DATA_ERROR;
            goto err;
        }

        if (ret == Z_DATA_ERROR || ret == Z_MEM_ERROR)
            goto err;

        unsigned just_wrote = unused - strm.avail_out;

        cdst += just_wrote;
        unused -= just_wrote;

        if (ret == Z_STREAM_END) {
            *dst_size -= unused;
            *dst = (char*)realloc(*dst, *dst_size);
            break;
        }

        if (strm.avail_out == 0) {
            if (unused > 0) {
                bgen_error("zlib-ng failed to uncompress (unknown error)");
                goto err;
            }

            unused = (unsigned)*dst_size;
            *dst_size += unused;
            *dst = (char*)realloc(*dst, *dst_size);
            cdst = (unsigned char*)*dst + unused;
        }
    }

    zng_inflateEnd(&strm);
    return 0;

err:
    zng_inflateEnd(&strm);
    return 1;
}
