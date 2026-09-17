/* Test-only hook: lets the C test suite observe the credential chain
 * through the same path production code uses.
 *
 * Lives in bgen_s3 rather than bgen/stream.cpp so the symbols it touches
 * (s3stream::Credentials, s3stream::ResolveCredentials) never cross a DLL
 * boundary: on a shared-library build, an unexported internal symbol of one
 * target cannot be resolved from another target's translation unit. */

#include "s3stream_internal.h"

#include "bgen/s3_export.h"

#include <cstdio>

#ifdef S3STREAM_ENABLE

extern "C" S3_EXPORT int s3stream_test_resolve_credentials(
    char* access, size_t access_len, char* secret, size_t secret_len, char* token,
    size_t token_len)
{
  s3stream::Credentials creds;
  if (!s3stream::ResolveCredentials(&creds)) {
    return 0;
  }
  std::snprintf(access, access_len, "%s", creds.access_key.c_str());
  std::snprintf(secret, secret_len, "%s", creds.secret_key.c_str());
  std::snprintf(token, token_len, "%s", creds.session_token.c_str());
  return 1;
}

#endif  // S3STREAM_ENABLE
