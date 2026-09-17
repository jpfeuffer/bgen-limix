/* bgen's stream API, implemented on top of the reusable s3stream module in
 * src/s3stream.  Kept as a thin adapter so s3stream can be lifted into a
 * standalone library without touching bgen. */

#include "stream.h"

#include "bgen/s3stream.h"

#include <stdio.h>

extern "C" {

int bgen_stream_is_s3(char const* path) { return s3stream_is_s3_uri(path); }

int bgen_stream_is_http(char const* path) { return s3stream_is_http_url(path); }

FILE* bgen_stream_open(char const* path)
{
    FILE* file = s3stream_open(path);
    if (!file && path && s3stream_is_remote(path))
        fprintf(stderr, "bgen: %s\n", s3stream_last_error());
    return file;
}

} /* extern "C" */

