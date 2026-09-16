#ifndef S3STREAM_HPP
#define S3STREAM_HPP

/* C++ RAII wrapper over the s3stream handle API (bgen/s3stream.h).
 *
 * Same lazy, Range-request-backed access on every platform, including
 * Windows -- see the "Handle API" section of s3stream.h for how this
 * compares to the FILE*-based s3stream_open(). Header-only: no separate
 * compiled sources, no ABI to keep stable across this header's changes.
 */

#include "bgen/s3stream.h"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>

namespace s3stream {

class Stream {
public:
  explicit Stream(const std::string& path) : handle_(open_or_throw(path, nullptr)) {}

  Stream(const std::string& path, const s3stream_credentials& creds)
      : handle_(open_or_throw(path, &creds)) {}

  ~Stream() { close(); }

  Stream(const Stream&) = delete;
  Stream& operator=(const Stream&) = delete;

  Stream(Stream&& other) noexcept : handle_(std::exchange(other.handle_, nullptr)) {}

  Stream& operator=(Stream&& other) noexcept {
    if (this != &other) {
      close();
      handle_ = std::exchange(other.handle_, nullptr);
    }
    return *this;
  }

  /* Reads up to n bytes into buf. Returns the number of bytes read (0 at
   * EOF), or -1 on error. */
  int64_t read(void* buf, size_t n) { return s3stream_handle_read(handle_, buf, n); }

  /* whence is SEEK_SET, SEEK_CUR or SEEK_END. Returns 0 on success, -1 on
   * error. */
  int seek(int64_t offset, int whence) { return s3stream_handle_seek(handle_, offset, whence); }

  int64_t tell() const { return s3stream_handle_tell(handle_); }

  bool eof() const { return s3stream_handle_eof(handle_) != 0; }

private:
  static s3stream_handle* open_or_throw(const std::string& path,
                                        const s3stream_credentials* creds) {
    s3stream_handle* handle = creds ? s3stream_handle_open_with_credentials(path.c_str(), creds)
                                     : s3stream_handle_open(path.c_str());
    if (!handle) {
      throw std::runtime_error(s3stream_last_error());
    }
    return handle;
  }

  void close() {
    if (handle_) {
      s3stream_handle_close(handle_);
      handle_ = nullptr;
    }
  }

  s3stream_handle* handle_;
};

}  // namespace s3stream

#endif /* S3STREAM_HPP */
