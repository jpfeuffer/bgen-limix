/* bgen's stream API, implemented on top of the reusable s3stream module in
 * src/s3stream.  Kept as a thin adapter so s3stream can be lifted into a
 * standalone library without touching bgen. */

#include "stream.h"

#include "s3stream/s3stream.h"

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

#ifdef BGEN_S3_SUPPORT

#include "s3stream/s3stream_internal.h"

/* Test hook: lets the C test suite observe the credential chain without
 * reaching into s3stream's internals. */
extern "C" int bgen_test_resolve_credentials(char* access, size_t access_len,
                                             char* secret, size_t secret_len,
                                             char* token, size_t token_len)
{
    s3stream::Credentials creds;
    if (!s3stream::ResolveCredentials(&creds))
        return 0;
    snprintf(access, access_len, "%s", creds.access_key.c_str());
    snprintf(secret, secret_len, "%s", creds.secret_key.c_str());
    snprintf(token, token_len, "%s", creds.session_token.c_str());
    return 1;
}

#endif /* BGEN_S3_SUPPORT */
