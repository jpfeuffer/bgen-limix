#include "stream.h"

#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifdef BGEN_S3_SUPPORT

#include <algorithm>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <curl/curl.h>

/* --------------------------------------------------------------------------
 * Credentials
 * -------------------------------------------------------------------------- */

struct S3Credentials
{
    std::string access_key;
    std::string secret_key;
    std::string session_token; /* empty for long-term credentials */
};

/* Simple string collector used by all internal curl calls. */
static size_t curl_collect_cb(char* ptr, size_t size, size_t nmemb, void* userdata)
{
    auto* buf = static_cast<std::string*>(userdata);
    buf->append(ptr, size * nmemb);
    return size * nmemb;
}

/* Extract a JSON string value for a given key from a flat JSON object.
 * This avoids pulling in a JSON library for a simple two-field response. */
static std::string json_get(const std::string& json, const char* key)
{
    /* Look for "key" : "value" */
    std::string needle = std::string("\"") + key + "\"";
    auto pos = json.find(needle);
    if (pos == std::string::npos)
        return {};
    pos = json.find('"', pos + needle.size());
    if (pos == std::string::npos)
        return {};
    ++pos;
    auto end = json.find('"', pos);
    if (end == std::string::npos)
        return {};
    return json.substr(pos, end - pos);
}

/* Fetch credentials from the EC2 Instance Metadata Service (IMDSv2).
 * Returns std::nullopt if the IMDS is not available (e.g. not on EC2).
 * Respects AWS_EC2_METADATA_DISABLED=true|1 so tests don't wait on non-EC2. */
static std::optional<S3Credentials> imds_credentials()
{
    char const* disabled = getenv("AWS_EC2_METADATA_DISABLED");
    if (disabled && (disabled[0] == '1' || disabled[0] == 't' || disabled[0] == 'T'))
        return std::nullopt;

    /* Step 1: obtain an IMDSv2 session token (PUT with TTL header). */
    std::string token;
    {
        CURL* c = curl_easy_init();
        if (!c)
            return std::nullopt;

        struct curl_slist* hdrs =
            curl_slist_append(nullptr, "X-aws-ec2-metadata-token-ttl-seconds: 21600");
        curl_easy_setopt(c, CURLOPT_URL,
                         "http://169.254.169.254/latest/api/token");
        curl_easy_setopt(c, CURLOPT_CUSTOMREQUEST, "PUT");
        curl_easy_setopt(c, CURLOPT_HTTPHEADER, hdrs);
        curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, curl_collect_cb);
        curl_easy_setopt(c, CURLOPT_WRITEDATA, &token);
        curl_easy_setopt(c, CURLOPT_TIMEOUT, 2L); /* fail fast when not on EC2 */
        curl_easy_setopt(c, CURLOPT_FAILONERROR, 1L);
        CURLcode res = curl_easy_perform(c);
        curl_slist_free_all(hdrs);
        curl_easy_cleanup(c);
        if (res != CURLE_OK || token.empty())
            return std::nullopt;
    }

    /* Step 2: discover the attached IAM role name. */
    std::string role;
    {
        CURL* c = curl_easy_init();
        if (!c)
            return std::nullopt;

        std::string token_hdr = "X-aws-ec2-metadata-token: " + token;
        struct curl_slist* hdrs = curl_slist_append(nullptr, token_hdr.c_str());
        curl_easy_setopt(c, CURLOPT_URL,
                         "http://169.254.169.254/latest/meta-data/iam/security-credentials/");
        curl_easy_setopt(c, CURLOPT_HTTPHEADER, hdrs);
        curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, curl_collect_cb);
        curl_easy_setopt(c, CURLOPT_WRITEDATA, &role);
        curl_easy_setopt(c, CURLOPT_TIMEOUT, 2L);
        curl_easy_setopt(c, CURLOPT_FAILONERROR, 1L);
        CURLcode res = curl_easy_perform(c);
        curl_slist_free_all(hdrs);
        curl_easy_cleanup(c);
        if (res != CURLE_OK || role.empty())
            return std::nullopt;
        /* Strip trailing newline */
        while (!role.empty() && (role.back() == '\n' || role.back() == '\r'))
            role.pop_back();
    }

    /* Step 3: retrieve the temporary credentials JSON for that role. */
    std::string creds_json;
    {
        CURL* c = curl_easy_init();
        if (!c)
            return std::nullopt;

        std::string token_hdr = "X-aws-ec2-metadata-token: " + token;
        struct curl_slist* hdrs = curl_slist_append(nullptr, token_hdr.c_str());
        std::string url =
            "http://169.254.169.254/latest/meta-data/iam/security-credentials/" + role;
        curl_easy_setopt(c, CURLOPT_URL, url.c_str());
        curl_easy_setopt(c, CURLOPT_HTTPHEADER, hdrs);
        curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, curl_collect_cb);
        curl_easy_setopt(c, CURLOPT_WRITEDATA, &creds_json);
        curl_easy_setopt(c, CURLOPT_TIMEOUT, 2L);
        curl_easy_setopt(c, CURLOPT_FAILONERROR, 1L);
        CURLcode res = curl_easy_perform(c);
        curl_slist_free_all(hdrs);
        curl_easy_cleanup(c);
        if (res != CURLE_OK || creds_json.empty())
            return std::nullopt;
    }

    S3Credentials creds;
    creds.access_key   = json_get(creds_json, "AccessKeyId");
    creds.secret_key   = json_get(creds_json, "SecretAccessKey");
    creds.session_token = json_get(creds_json, "Token");

    if (creds.access_key.empty() || creds.secret_key.empty())
        return std::nullopt;

    return creds;
}

/* Parse ~/.aws/credentials (or AWS_SHARED_CREDENTIALS_FILE).
 * Profile is chosen by AWS_PROFILE (default: "default"). */
static std::optional<S3Credentials> shared_file_credentials()
{
    char const* path_env = getenv("AWS_SHARED_CREDENTIALS_FILE");
    std::string path;
    if (path_env && *path_env) {
        path = path_env;
    } else {
        char const* home = getenv("HOME");
#ifdef _WIN32
        if (!home) home = getenv("USERPROFILE");
#endif
        if (!home || !*home)
            return std::nullopt;
        path = std::string(home) + "/.aws/credentials";
    }

    FILE* f = fopen(path.c_str(), "r");
    if (!f)
        return std::nullopt;

    char const* profile_env = getenv("AWS_PROFILE");
    std::string section_target =
        std::string("[") + (profile_env && *profile_env ? profile_env : "default") + "]";

    bool          in_section = false;
    S3Credentials creds;
    char          line[1024];

    while (fgets(line, static_cast<int>(sizeof(line)), f)) {
        int len = static_cast<int>(strlen(line));
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r' ||
                           line[len-1] == ' '  || line[len-1] == '\t'))
            line[--len] = '\0';
        if (len == 0) continue;

        if (line[0] == '[') { in_section = (section_target == line); continue; }
        if (!in_section || line[0] == '#' || line[0] == ';') continue;

        char* eq = strchr(line, '=');
        if (!eq) continue;
        *eq = '\0';
        char* k = line;
        char* v = eq + 1;

        while (*k == ' ' || *k == '\t') ++k;
        int klen = static_cast<int>(strlen(k));
        while (klen > 0 && (k[klen-1] == ' ' || k[klen-1] == '\t')) k[--klen] = '\0';
        while (*v == ' ' || *v == '\t') ++v;
        int vlen = static_cast<int>(strlen(v));
        while (vlen > 0 && (v[vlen-1] == ' ' || v[vlen-1] == '\t')) v[--vlen] = '\0';

        if      (strcmp(k, "aws_access_key_id")    == 0) creds.access_key    = v;
        else if (strcmp(k, "aws_secret_access_key") == 0) creds.secret_key    = v;
        else if (strcmp(k, "aws_session_token")     == 0) creds.session_token = v;
    }
    fclose(f);

    if (creds.access_key.empty() || creds.secret_key.empty())
        return std::nullopt;
    return creds;
}

/* Fetch credentials from the ECS container metadata endpoint.
 * Supports AWS_CONTAINER_CREDENTIALS_RELATIVE_URI and _FULL_URI. */
static std::optional<S3Credentials> ecs_credentials()
{
    char const* rel  = getenv("AWS_CONTAINER_CREDENTIALS_RELATIVE_URI");
    char const* full = getenv("AWS_CONTAINER_CREDENTIALS_FULL_URI");
    std::string url;
    if      (rel  && *rel)  url = std::string("http://169.254.170.2") + rel;
    else if (full && *full) url = full;
    else                    return std::nullopt;

    std::string body;
    CURL* c = curl_easy_init();
    if (!c) return std::nullopt;
    curl_easy_setopt(c, CURLOPT_URL,           url.c_str());
    curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, curl_collect_cb);
    curl_easy_setopt(c, CURLOPT_WRITEDATA,     &body);
    curl_easy_setopt(c, CURLOPT_TIMEOUT,       2L);
    curl_easy_setopt(c, CURLOPT_FAILONERROR,   1L);
    CURLcode res = curl_easy_perform(c);
    curl_easy_cleanup(c);
    if (res != CURLE_OK || body.empty()) return std::nullopt;

    S3Credentials creds;
    creds.access_key    = json_get(body, "AccessKeyId");
    creds.secret_key    = json_get(body, "SecretAccessKey");
    creds.session_token = json_get(body, "Token");
    if (creds.access_key.empty() || creds.secret_key.empty())
        return std::nullopt;
    return creds;
}

/* Resolve credentials — full AWS chain:
 *   1. Environment variables (AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY)
 *   2. Shared credentials file (~/.aws/credentials or AWS_SHARED_CREDENTIALS_FILE)
 *   3. ECS container metadata (AWS_CONTAINER_CREDENTIALS_RELATIVE/FULL_URI)
 *   4. EC2 Instance Metadata Service v2 (IMDSv2)
 *   5. None → unauthenticated
 */
static std::optional<S3Credentials> resolve_credentials()
{
    /* 1. Environment variables */
    char const* access = getenv("AWS_ACCESS_KEY_ID");
    char const* secret = getenv("AWS_SECRET_ACCESS_KEY");
    if (access && *access && secret && *secret) {
        S3Credentials c;
        c.access_key = access;
        c.secret_key = secret;
        char const* token = getenv("AWS_SESSION_TOKEN");
        if (token && *token) c.session_token = token;
        return c;
    }

    /* 2. Shared credentials file */
    auto file_creds = shared_file_credentials();
    if (file_creds) return file_creds;

    /* 3. ECS task role */
    auto ecs = ecs_credentials();
    if (ecs) return ecs;

    /* 4. EC2 IMDS */
    auto imds = imds_credentials();
    if (imds) return imds;

    return std::nullopt;
}

/* --------------------------------------------------------------------------
 * S3 cookie: holds state for a FILE* backed by S3 Range requests.
 * -------------------------------------------------------------------------- */

static constexpr int64_t CHUNK_SIZE = 4 * 1024 * 1024; /* 4 MiB read-ahead */

struct S3Cookie
{
    std::string url;      /* full URL for range requests */
    CURL*       curl;     /* persistent handle */
    struct curl_slist* extra_headers; /* e.g. x-amz-security-token */

    int64_t file_size;
    int64_t position;
    bool    eof_flag;

    std::vector<char> buffer;
    int64_t           buf_start;
    int64_t           buf_valid;
};

/* --------------------------------------------------------------------------
 * libcurl range read helper
 * -------------------------------------------------------------------------- */

static size_t curl_write_cb(char* ptr, size_t size, size_t nmemb, void* userdata)
{
    auto* ctx = static_cast<std::vector<char>*>(userdata);
    ctx->insert(ctx->end(), ptr, ptr + size * nmemb);
    return size * nmemb;
}

static int s3_range_read(S3Cookie* cookie, int64_t range_start, int64_t range_end,
                         std::vector<char>& dest)
{
    dest.clear();
    char range_hdr[128];
    snprintf(range_hdr, sizeof(range_hdr), "%" PRId64 "-%" PRId64, range_start, range_end);
    curl_easy_setopt(cookie->curl, CURLOPT_RANGE, range_hdr);
    curl_easy_setopt(cookie->curl, CURLOPT_WRITEDATA, &dest);

    CURLcode res = curl_easy_perform(cookie->curl);
    if (res != CURLE_OK) {
        fprintf(stderr, "bgen: S3 range read failed: %s\n", curl_easy_strerror(res));
        return -1;
    }
    return 0;
}

#ifndef _WIN32
static int s3_fill_buffer(S3Cookie* cookie)
{
    int64_t start = cookie->position;
    int64_t end   = std::min(start + CHUNK_SIZE - 1, cookie->file_size - 1);
    if (start > end) {
        cookie->buf_valid = 0;
        return 0;
    }
    std::vector<char> data;
    if (s3_range_read(cookie, start, end, data) != 0)
        return -1;
    cookie->buffer    = std::move(data);
    cookie->buf_start = start;
    cookie->buf_valid = static_cast<int64_t>(cookie->buffer.size());
    return 0;
}
#endif /* !_WIN32 */

/* --------------------------------------------------------------------------
 * funopen (macOS/BSD) callbacks
 * -------------------------------------------------------------------------- */

#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)

static int s3_read_fn(void* cookie_ptr, char* buf, int nbytes)
{
    auto* c = static_cast<S3Cookie*>(cookie_ptr);
    if (nbytes <= 0) return 0;
    if (c->position >= c->file_size) { c->eof_flag = true; return 0; }

    int total_read = 0;
    while (total_read < nbytes && c->position < c->file_size) {
        if (c->buf_valid > 0 && c->position >= c->buf_start &&
            c->position < c->buf_start + c->buf_valid) {
            int64_t off     = c->position - c->buf_start;
            int     to_copy = static_cast<int>(
                std::min(static_cast<int64_t>(nbytes - total_read), c->buf_valid - off));
            memcpy(buf + total_read, c->buffer.data() + off, static_cast<size_t>(to_copy));
            c->position  += to_copy;
            total_read   += to_copy;
        } else {
            if (s3_fill_buffer(c) != 0) return -1;
            if (c->buf_valid == 0) { c->eof_flag = true; break; }
        }
    }
    if (c->position >= c->file_size) c->eof_flag = true;
    return total_read;
}

static fpos_t s3_seek_fn(void* cookie_ptr, fpos_t offset, int whence)
{
    auto* c = static_cast<S3Cookie*>(cookie_ptr);
    int64_t new_pos;
    switch (whence) {
    case SEEK_SET: new_pos = offset; break;
    case SEEK_CUR: new_pos = c->position + offset; break;
    case SEEK_END: new_pos = c->file_size + offset; break;
    default: return -1;
    }
    if (new_pos < 0) return -1;
    c->position = new_pos;
    c->eof_flag = false;
    return static_cast<fpos_t>(c->position);
}

static int s3_close_fn(void* cookie_ptr)
{
    auto* c = static_cast<S3Cookie*>(cookie_ptr);
    curl_easy_cleanup(c->curl);
    curl_slist_free_all(c->extra_headers);
    delete c;
    return 0;
}

static FILE* s3_create_file(S3Cookie* cookie)
{
    return funopen(cookie, s3_read_fn, nullptr, s3_seek_fn, s3_close_fn);
}

/* --------------------------------------------------------------------------
 * fopencookie (Linux/glibc) callbacks
 * -------------------------------------------------------------------------- */

#elif defined(__linux__) && defined(__GLIBC__)

static ssize_t s3_cookie_read(void* cookie_ptr, char* buf, size_t size)
{
    auto* c = static_cast<S3Cookie*>(cookie_ptr);
    if (size == 0) return 0;
    if (c->position >= c->file_size) { c->eof_flag = true; return 0; }

    ssize_t total_read = 0;
    while (static_cast<size_t>(total_read) < size && c->position < c->file_size) {
        if (c->buf_valid > 0 && c->position >= c->buf_start &&
            c->position < c->buf_start + c->buf_valid) {
            int64_t off     = c->position - c->buf_start;
            ssize_t to_copy = static_cast<ssize_t>(
                std::min(static_cast<int64_t>(size - static_cast<size_t>(total_read)),
                         c->buf_valid - off));
            memcpy(buf + total_read, c->buffer.data() + off, static_cast<size_t>(to_copy));
            c->position += to_copy;
            total_read  += to_copy;
        } else {
            if (s3_fill_buffer(c) != 0) return -1;
            if (c->buf_valid == 0) { c->eof_flag = true; break; }
        }
    }
    if (c->position >= c->file_size) c->eof_flag = true;
    return total_read;
}

static int s3_cookie_seek(void* cookie_ptr, off64_t* offset, int whence)
{
    auto* c = static_cast<S3Cookie*>(cookie_ptr);
    int64_t new_pos;
    switch (whence) {
    case SEEK_SET: new_pos = *offset; break;
    case SEEK_CUR: new_pos = c->position + *offset; break;
    case SEEK_END: new_pos = c->file_size + *offset; break;
    default: return -1;
    }
    if (new_pos < 0) return -1;
    c->position = new_pos;
    c->eof_flag = false;
    *offset = new_pos;
    return 0;
}

static int s3_cookie_close(void* cookie_ptr)
{
    auto* c = static_cast<S3Cookie*>(cookie_ptr);
    curl_easy_cleanup(c->curl);
    curl_slist_free_all(c->extra_headers);
    delete c;
    return 0;
}

static FILE* s3_create_file(S3Cookie* cookie)
{
    cookie_io_functions_t funcs = {};
    funcs.read  = s3_cookie_read;
    funcs.seek  = s3_cookie_seek;
    funcs.close = s3_cookie_close;
    return fopencookie(cookie, "rb", funcs);
}

#elif defined(_WIN32)

#include <windows.h>
#include <io.h>
#include <fcntl.h>

/* Windows has neither funopen nor fopencookie.  Download the whole S3 object
 * to a temporary file and return a seekable FILE* backed by it.  The temp
 * file is opened with FILE_FLAG_DELETE_ON_CLOSE so it is removed when the
 * returned FILE* is closed. */
static FILE* s3_create_file(S3Cookie* cookie)
{
    char tmpdir[MAX_PATH];
    char tmppath[MAX_PATH];
    if (!GetTempPathA(MAX_PATH, tmpdir) ||
        !GetTempFileNameA(tmpdir, "bgn", 0, tmppath)) {
        fprintf(stderr, "bgen: failed to obtain temp path for S3 stream\n");
        return nullptr;
    }

    FILE* tmp = fopen(tmppath, "w+b");
    if (!tmp) {
        fprintf(stderr, "bgen: failed to open temp file for S3 stream\n");
        _unlink(tmppath);
        return nullptr;
    }

    /* Download in CHUNK_SIZE pieces via range requests. */
    int64_t pos = 0;
    while (pos < cookie->file_size) {
        int64_t end = std::min(pos + CHUNK_SIZE - 1, cookie->file_size - 1);
        std::vector<char> data;
        if (s3_range_read(cookie, pos, end, data) != 0) {
            fclose(tmp);
            _unlink(tmppath);
            return nullptr; /* caller owns cookie cleanup on failure */
        }
        if (!data.empty())
            fwrite(data.data(), 1, data.size(), tmp);
        if (data.empty()) break; /* server returned 0 bytes: treat as EOF */
        pos += static_cast<int64_t>(data.size());
    }

    /* Cookie is fully consumed; clean it up now. */
    curl_easy_cleanup(cookie->curl);
    curl_slist_free_all(cookie->extra_headers);
    delete cookie;

    fflush(tmp);
    fclose(tmp);

    /* Reopen with FILE_FLAG_DELETE_ON_CLOSE so the file disappears on close. */
    HANDLE h = CreateFileA(tmppath, GENERIC_READ,
                           FILE_SHARE_READ | FILE_SHARE_DELETE,
                           nullptr, OPEN_EXISTING,
                           FILE_FLAG_DELETE_ON_CLOSE, nullptr);
    if (h == INVALID_HANDLE_VALUE) {
        _unlink(tmppath);
        return nullptr;
    }
    int fd = _open_osfhandle(reinterpret_cast<intptr_t>(h),
                             _O_RDONLY | _O_BINARY);
    if (fd < 0) {
        CloseHandle(h);
        _unlink(tmppath);
        return nullptr;
    }
    return _fdopen(fd, "rb");
}

#else
#error "S3 streaming requires funopen (macOS/BSD), fopencookie (Linux/glibc), or _WIN32"
#endif

/* --------------------------------------------------------------------------
 * S3 URI parsing
 * -------------------------------------------------------------------------- */

static int parse_s3_uri(char const* uri, std::string& bucket, std::string& key)
{
    std::string_view sv(uri);
    if (!sv.starts_with("s3://")) return -1;
    sv.remove_prefix(5);
    auto slash = sv.find('/');
    if (slash == std::string_view::npos || slash == 0) return -1;
    bucket = std::string(sv.substr(0, slash));
    key    = std::string(sv.substr(slash + 1));
    return key.empty() ? -1 : 0;
}

/* --------------------------------------------------------------------------
 * HEAD request via curl: returns Content-Length, or -1 on failure.
 * Uses SigV4 when creds is non-null, unauthenticated otherwise.
 * extra_headers is shared with the cookie (must outlive this call).
 * -------------------------------------------------------------------------- */

static int64_t s3_head_object(const std::string& url, const std::string& region,
                              const std::optional<S3Credentials>& creds,
                              struct curl_slist* extra_headers)
{
    CURL* c = curl_easy_init();
    if (!c) return -1;

    std::string discard;
    curl_easy_setopt(c, CURLOPT_URL, url.c_str());
    curl_easy_setopt(c, CURLOPT_NOBODY, 1L); /* HEAD */
    curl_easy_setopt(c, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(c, CURLOPT_FAILONERROR, 1L);
    curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, curl_collect_cb);
    curl_easy_setopt(c, CURLOPT_WRITEDATA, &discard);

    if (creds) {
        std::string sigv4 = "aws:amz:" + region + ":s3";
        curl_easy_setopt(c, CURLOPT_AWS_SIGV4, sigv4.c_str());
        std::string userpwd = creds->access_key + ":" + creds->secret_key;
        curl_easy_setopt(c, CURLOPT_USERPWD, userpwd.c_str());
    }

    if (extra_headers)
        curl_easy_setopt(c, CURLOPT_HTTPHEADER, extra_headers);

    CURLcode res = curl_easy_perform(c);
    if (res != CURLE_OK) {
        fprintf(stderr, "bgen: S3 HEAD failed for %s: %s\n", url.c_str(),
                curl_easy_strerror(res));
        curl_easy_cleanup(c);
        return -1;
    }

    curl_off_t content_length = -1;
    curl_easy_getinfo(c, CURLINFO_CONTENT_LENGTH_DOWNLOAD_T, &content_length);
    curl_easy_cleanup(c);
    return static_cast<int64_t>(content_length);
}

/* --------------------------------------------------------------------------
 * Main open function
 * -------------------------------------------------------------------------- */

static FILE* open_s3_stream(char const* path)
{
    std::string bucket, key;
    if (parse_s3_uri(path, bucket, key) != 0) {
        fprintf(stderr, "bgen: invalid S3 URI: %s\n", path);
        return nullptr;
    }

    /* ── Region ─────────────────────────────────────────────────────────── */
    char const* region_env = getenv("AWS_DEFAULT_REGION");
    std::string region = region_env ? region_env : "us-east-1";

    /* ── Endpoint / URL ──────────────────────────────────────────────────── */
    char const* endpoint_env = getenv("AWS_ENDPOINT_URL");
    std::string url;
    if (endpoint_env) {
        std::string ep(endpoint_env);
        if (!ep.empty() && ep.back() == '/') ep.pop_back();
        url = ep + "/" + bucket + "/" + key; /* path-style */
    } else {
        /* virtual-hosted style: https://bucket.s3.region.amazonaws.com/key */
        url = "https://" + bucket + ".s3." + region + ".amazonaws.com/" + key;
    }

    /* ── Credentials (chain: env → IMDS → unauthenticated) ──────────────── */
    bool no_sign = false;
    {
        char const* v = getenv("AWS_NO_SIGN_REQUEST");
        if (v && (v[0] == '1' || v[0] == 't' || v[0] == 'T' || v[0] == 'y' || v[0] == 'Y'))
            no_sign = true;
    }

    std::optional<S3Credentials> creds;
    if (!no_sign) {
        creds = resolve_credentials();
        if (!creds)
            fprintf(stderr,
                    "bgen: no AWS credentials found (env, IMDS); "
                    "trying unauthenticated (set AWS_NO_SIGN_REQUEST=1 to suppress this)\n");
    }

    /* ── Build the security-token header list (used by both HEAD and GETs) ─ */
    struct curl_slist* extra_headers = nullptr;
    if (creds && !creds->session_token.empty()) {
        std::string hdr = "x-amz-security-token: " + creds->session_token;
        extra_headers = curl_slist_append(nullptr, hdr.c_str());
    }

    /* ── HEAD to get content-length ──────────────────────────────────────── */
    int64_t file_size = s3_head_object(url, region, creds, extra_headers);
    if (file_size <= 0) {
        curl_slist_free_all(extra_headers);
        return nullptr;
    }

    /* ── Build persistent curl handle for range GETs ─────────────────────── */
    CURL* ch = curl_easy_init();
    if (!ch) {
        curl_slist_free_all(extra_headers);
        return nullptr;
    }
    curl_easy_setopt(ch, CURLOPT_URL, url.c_str());
    curl_easy_setopt(ch, CURLOPT_WRITEFUNCTION, curl_write_cb);
    curl_easy_setopt(ch, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(ch, CURLOPT_FAILONERROR, 1L);

    if (creds) {
        std::string sigv4 = "aws:amz:" + region + ":s3";
        curl_easy_setopt(ch, CURLOPT_AWS_SIGV4, sigv4.c_str());
        /* USERPWD must stay alive; store in cookie below */
    }
    if (extra_headers)
        curl_easy_setopt(ch, CURLOPT_HTTPHEADER, extra_headers);

    /* ── Assemble cookie ─────────────────────────────────────────────────── */
    auto* cookie = new S3Cookie();
    cookie->url            = url;
    cookie->curl           = ch;
    cookie->extra_headers  = extra_headers;
    cookie->file_size      = file_size;
    cookie->position       = 0;
    cookie->eof_flag       = false;
    cookie->buf_start      = 0;
    cookie->buf_valid      = 0;

    /* USERPWD string must outlive the curl handle — store it in the cookie's url
     * field is wrong, so keep it as a separate member via a small lambda trick
     * by setting it directly now (the string is owned by creds which is local,
     * so we must persist it via curl_easy_setopt on the already-stored handle). */
    if (creds) {
        /* We need the userpwd string to survive beyond this function scope.
         * Curl copies the string internally on setopt, so this is safe. */
        std::string userpwd = creds->access_key + ":" + creds->secret_key;
        curl_easy_setopt(ch, CURLOPT_USERPWD, userpwd.c_str());
    }

    FILE* fp = s3_create_file(cookie);
    if (!fp) {
        fprintf(stderr, "bgen: failed to create S3 stream\n");
        curl_easy_cleanup(ch);
        curl_slist_free_all(extra_headers);
        delete cookie;
        return nullptr;
    }

    return fp;
}

/* --------------------------------------------------------------------------
 * Test hook: exposes resolve_credentials() to C unit tests without requiring
 * access to internal headers.  Compiled whenever BGEN_S3_SUPPORT is enabled.
 * -------------------------------------------------------------------------- */
extern "C" int bgen_test_resolve_credentials(char* access, size_t access_len,
                                             char* secret, size_t secret_len,
                                             char* token,  size_t token_len)
{
    auto creds = resolve_credentials();
    if (!creds) return 0;
    snprintf(access, access_len, "%s", creds->access_key.c_str());
    snprintf(secret, secret_len, "%s", creds->secret_key.c_str());
    snprintf(token,  token_len,  "%s", creds->session_token.c_str());
    return 1;
}

#endif /* BGEN_S3_SUPPORT */

/* --------------------------------------------------------------------------
 * Public C API
 * -------------------------------------------------------------------------- */

extern "C" {

int bgen_stream_is_s3(char const* path)
{
    return (path && strncmp(path, "s3://", 5) == 0) ? 1 : 0;
}

FILE* bgen_stream_open(char const* path)
{
    if (!path)
        return nullptr;

    if (bgen_stream_is_s3(path)) {
#ifdef BGEN_S3_SUPPORT
        return open_s3_stream(path);
#else
        fprintf(stderr, "bgen: S3 support not enabled (rebuild with -DBGEN_ENABLE_S3=ON)\n");
        return nullptr;
#endif
    }

    return fopen(path, "rb");
}

} /* extern "C" */
