/* Unit tests for the shared-credentials-file provider.
 *
 * Uses bgen_test_resolve_credentials() (compiled into libbgen when
 * BGEN_ENABLE_S3=ON) to exercise the credential chain without a real
 * S3 connection.
 *
 * env vars manipulated here are restored before the test returns.
 */
#include "cass.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#  include <windows.h>
#  include <io.h>
#  include <fcntl.h>
/* Portability shims for POSIX env/file helpers on Windows. */
static int bgen_setenv(const char* name, const char* value)
{
    return _putenv_s(name, value) == 0 ? 0 : -1;
}
static int bgen_unsetenv(const char* name)
{
    return _putenv_s(name, "") == 0 ? 0 : -1;
}
static int bgen_mkstemp(char* buf, size_t bufsz)
{
    char dir[MAX_PATH];
    char path[MAX_PATH];
    if (!GetTempPathA(MAX_PATH, dir)) return -1;
    if (!GetTempFileNameA(dir, "bgn", 0, path)) return -1;
    if (strlen(path) + 1 > bufsz) return -1;
    memcpy(buf, path, strlen(path) + 1);
    return _open(path, _O_RDWR | _O_BINARY | _O_CREAT | _O_TRUNC,
                 _S_IREAD | _S_IWRITE);
}
#  define setenv(n, v, o) bgen_setenv((n), (v))
#  define unsetenv(n)     bgen_unsetenv((n))
#  define close(fd)       _close((fd))
#  define unlink(p)       _unlink((p))
#else
#  include <unistd.h>
#endif

/* Forward declaration — symbol is provided by libbgen when S3 is enabled. */
extern int bgen_test_resolve_credentials(char* access, size_t access_len,
                                         char* secret, size_t secret_len,
                                         char* token,  size_t token_len);

static void write_creds_file(const char* path)
{
    FILE* f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", path); exit(1); }
    fputs("[default]\n", f);
    fputs("aws_access_key_id     = AKIATEST0000000001\n", f);
    fputs("aws_secret_access_key = secret0000000000000000000000000001\n", f);
    fputs("aws_session_token     = SESSTOKEN01\n", f);
    fputs("\n", f);
    fputs("[myprofile]\n", f);
    fputs("aws_access_key_id     = AKIATEST0000000002\n", f);
    fputs("aws_secret_access_key = secret0000000000000000000000000002\n", f);
    /* no session_token in myprofile — intentional */
    fclose(f);
}

static void save_env(const char* name, char* buf, size_t len)
{
    const char* v = getenv(name);
    if (v) { strncpy(buf, v, len - 1); buf[len - 1] = '\0'; }
    else   { buf[0] = '\0'; }
}

int main(void)
{
    /* ── Save original env state ────────────────────────────────────────── */
    char orig_access[256], orig_secret[256], orig_profile[256];
    char orig_creds_file[1024], orig_meta_disabled[32];
    save_env("AWS_ACCESS_KEY_ID",           orig_access,        sizeof orig_access);
    save_env("AWS_SECRET_ACCESS_KEY",       orig_secret,        sizeof orig_secret);
    save_env("AWS_PROFILE",                 orig_profile,       sizeof orig_profile);
    save_env("AWS_SHARED_CREDENTIALS_FILE", orig_creds_file,    sizeof orig_creds_file);
    save_env("AWS_EC2_METADATA_DISABLED",   orig_meta_disabled, sizeof orig_meta_disabled);

    /* Clear env-var provider so file provider is reached. */
    unsetenv("AWS_ACCESS_KEY_ID");
    unsetenv("AWS_SECRET_ACCESS_KEY");
    unsetenv("AWS_PROFILE");
    /* Disable remote providers so tests are instant and deterministic. */
    setenv("AWS_EC2_METADATA_DISABLED", "true", 1);
    unsetenv("AWS_CONTAINER_CREDENTIALS_RELATIVE_URI");
    unsetenv("AWS_CONTAINER_CREDENTIALS_FULL_URI");

    /* ── Write a temp credentials file ─────────────────────────────────── */
#ifdef _WIN32
    char tmppath[MAX_PATH];
    int  fd = bgen_mkstemp(tmppath, sizeof tmppath);
#else
    char tmppath[] = "/tmp/bgen_test_creds_XXXXXX";
    int  fd = mkstemp(tmppath);
#endif
    if (fd < 0) { fprintf(stderr, "cannot create temp file\n"); return 1; }
    close(fd);
    write_creds_file(tmppath);
    setenv("AWS_SHARED_CREDENTIALS_FILE", tmppath, 1);

    char access[256], secret[256], token[256];

    /* ── Test 1: default profile ────────────────────────────────────────── */
    int found = bgen_test_resolve_credentials(
        access, sizeof access, secret, sizeof secret, token, sizeof token);
    cass_cond(found == 1);
    cass_cond(strcmp(access, "AKIATEST0000000001") == 0);
    cass_cond(strcmp(secret, "secret0000000000000000000000000001") == 0);
    cass_cond(strcmp(token,  "SESSTOKEN01") == 0);

    /* ── Test 2: named profile via AWS_PROFILE ──────────────────────────── */
    setenv("AWS_PROFILE", "myprofile", 1);
    found = bgen_test_resolve_credentials(
        access, sizeof access, secret, sizeof secret, token, sizeof token);
    cass_cond(found == 1);
    cass_cond(strcmp(access, "AKIATEST0000000002") == 0);
    cass_cond(strcmp(secret, "secret0000000000000000000000000002") == 0);
    cass_cond(strcmp(token,  "") == 0); /* no session token in myprofile */
    unsetenv("AWS_PROFILE");

    /* ── Test 3: non-existent profile → no credentials ──────────────────── */
    setenv("AWS_PROFILE", "doesnotexist", 1);
    found = bgen_test_resolve_credentials(
        access, sizeof access, secret, sizeof secret, token, sizeof token);
    cass_cond(found == 0);
    unsetenv("AWS_PROFILE");

    /* ── Test 4: env-var provider takes precedence over file ────────────── */
    setenv("AWS_ACCESS_KEY_ID",     "ENVKEY",    1);
    setenv("AWS_SECRET_ACCESS_KEY", "ENVSECRET", 1);
    found = bgen_test_resolve_credentials(
        access, sizeof access, secret, sizeof secret, token, sizeof token);
    cass_cond(found == 1);
    cass_cond(strcmp(access, "ENVKEY")    == 0);
    cass_cond(strcmp(secret, "ENVSECRET") == 0);
    unsetenv("AWS_ACCESS_KEY_ID");
    unsetenv("AWS_SECRET_ACCESS_KEY");

    /* ── Cleanup ─────────────────────────────────────────────────────────── */
    unlink(tmppath);
    unsetenv("AWS_SHARED_CREDENTIALS_FILE");
    unsetenv("AWS_EC2_METADATA_DISABLED");

    if (orig_access[0])        setenv("AWS_ACCESS_KEY_ID",           orig_access,        1);
    if (orig_secret[0])        setenv("AWS_SECRET_ACCESS_KEY",       orig_secret,        1);
    if (orig_profile[0])       setenv("AWS_PROFILE",                 orig_profile,       1);
    if (orig_creds_file[0])    setenv("AWS_SHARED_CREDENTIALS_FILE", orig_creds_file,    1);
    if (orig_meta_disabled[0]) setenv("AWS_EC2_METADATA_DISABLED",   orig_meta_disabled, 1);

    return cass_status();
}
