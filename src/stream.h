#ifndef BGEN_STREAM_H
#define BGEN_STREAM_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Open a file for reading. Handles both local paths and S3 URIs (s3://bucket/key).
 *
 * For local files, this is equivalent to fopen(path, "rb").
 * For S3 URIs, this returns a FILE* backed by HTTP Range requests,
 * supporting fread, fseek, ftell, and feof transparently.
 *
 * S3 credentials are read from environment variables:
 *   - AWS_ACCESS_KEY_ID (required for S3)
 *   - AWS_SECRET_ACCESS_KEY (required for S3)
 *   - AWS_DEFAULT_REGION (optional, defaults to "us-east-1")
 *   - AWS_ENDPOINT_URL (optional, for S3-compatible services like MinIO)
 *
 * The returned FILE* must be closed with fclose() as usual.
 *
 * @param path Local file path or S3 URI.
 * @return FILE* on success, NULL on failure.
 */
FILE* bgen_stream_open(char const* path);

/**
 * Check if a path is an S3 URI.
 * @return non-zero if path starts with "s3://".
 */
int bgen_stream_is_s3(char const* path);

#ifdef __cplusplus
}
#endif

#endif
