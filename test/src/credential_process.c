/* Unit tests for the credential_process provider (POSIX only: the helper
 * commands are printf/false). */
#if !defined(_WIN32) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE 200809L
#endif
#include "cass.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
int main(void) { return 0; }
#else
#include <unistd.h>

extern int s3stream_test_resolve_credentials(char* access, size_t access_len,
                                             char* secret, size_t secret_len,
                                             char* token,  size_t token_len);

static void write_config(const char* path)
{
    FILE* f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", path); exit(1); }
    fputs("[profile good]\n", f);
    fputs("credential_process = printf '{\"Version\":1,\"AccessKeyId\":\"AKIAPROC\","
          "\"SecretAccessKey\":\"procsecret\",\"SessionToken\":\"proctok\","
          "\"Expiration\":\"2099-01-01T00:00:00Z\"}'\n", f);
    fputs("\n[profile failing]\n", f);
    fputs("credential_process = false\n", f);
    fputs("\n[profile badversion]\n", f);
    fputs("credential_process = printf '{\"Version\":2,\"AccessKeyId\":\"A\","
          "\"SecretAccessKey\":\"B\"}'\n", f);
    fputs("\n[profile garbage]\n", f);
    fputs("credential_process = printf 'not json'\n", f);
    fputs("\n[profile both]\n", f);
    fputs("aws_access_key_id = AKIASTATIC\n", f);
    fputs("aws_secret_access_key = staticsecret\n", f);
    fputs("credential_process = printf '{\"Version\":1,\"AccessKeyId\":\"AKIAPROC2\","
          "\"SecretAccessKey\":\"procsecret2\"}'\n", f);
    fclose(f);
}

int main(void)
{
    unsetenv("AWS_ACCESS_KEY_ID");
    unsetenv("AWS_SECRET_ACCESS_KEY");
    unsetenv("AWS_SESSION_TOKEN");
    unsetenv("AWS_CONTAINER_CREDENTIALS_RELATIVE_URI");
    unsetenv("AWS_CONTAINER_CREDENTIALS_FULL_URI");
    setenv("AWS_EC2_METADATA_DISABLED", "true", 1);
    setenv("AWS_SHARED_CREDENTIALS_FILE", "/nonexistent/bgen_test_credentials", 1);

    char path[] = "/tmp/bgen_test_cred_process_XXXXXX";
    int fd = mkstemp(path);
    if (fd < 0) { fprintf(stderr, "cannot create temp file\n"); return 1; }
    close(fd);
    write_config(path);
    setenv("AWS_CONFIG_FILE", path, 1);

    char access[256], secret[256], token[256];

    setenv("AWS_PROFILE", "good", 1);
    int found = s3stream_test_resolve_credentials(
        access, sizeof access, secret, sizeof secret, token, sizeof token);
    cass_cond(found == 1);
    cass_cond(strcmp(access, "AKIAPROC") == 0);
    cass_cond(strcmp(secret, "procsecret") == 0);
    cass_cond(strcmp(token, "proctok") == 0);

    setenv("AWS_PROFILE", "failing", 1);
    found = s3stream_test_resolve_credentials(
        access, sizeof access, secret, sizeof secret, token, sizeof token);
    cass_cond(found == 0);

    setenv("AWS_PROFILE", "badversion", 1);
    found = s3stream_test_resolve_credentials(
        access, sizeof access, secret, sizeof secret, token, sizeof token);
    cass_cond(found == 0);

    setenv("AWS_PROFILE", "garbage", 1);
    found = s3stream_test_resolve_credentials(
        access, sizeof access, secret, sizeof secret, token, sizeof token);
    cass_cond(found == 0);

    /* credential_process outranks static keys in the same config profile. */
    setenv("AWS_PROFILE", "both", 1);
    found = s3stream_test_resolve_credentials(
        access, sizeof access, secret, sizeof secret, token, sizeof token);
    cass_cond(found == 1);
    cass_cond(strcmp(access, "AKIAPROC2") == 0);
    cass_cond(strcmp(token, "") == 0);

    unlink(path);
    return cass_status();
}
#endif
