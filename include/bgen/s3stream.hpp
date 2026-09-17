#ifndef S3STREAM_HPP
#define S3STREAM_HPP

/* C++ RAII wrapper over the s3stream handle API (bgen/s3stream.h).
 *
 * Same lazy, Range-request-backed access on every platform, including
 * Windows -- see the "Handle API" section of s3stream.h for how this
 * compares to the FILE*-based s3stream_open(). Header-only: no separate
 * compiled sources, no ABI to keep stable across this header's changes.
 *
 * Three layers, from thinnest to most convenient:
 *   Stream     -- RAII handle with read/seek/tell/eof
 *   StreamBuf  -- std::streambuf, for composing with stream adaptors
 *   IStream    -- std::istream, for getline/operator>> on a remote path
 */

#include "bgen/s3stream.h"

#include <cstdint>
#include <cstdio>
#include <istream>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <utility>
#include <vector>

namespace s3stream {

/* Default read-ahead, in bytes. One range request per buffer refill, so
 * larger values trade memory for fewer round trips. */
const std::size_t kDefaultBufferSize = 1 << 16;

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
  int64_t read(void* buf, size_t n) { return stream_handle_read(handle_, buf, n); }

  /* whence is SEEK_SET, SEEK_CUR or SEEK_END. Returns 0 on success, -1 on
   * error. */
  int seek(int64_t offset, int whence) { return stream_handle_seek(handle_, offset, whence); }

  int64_t tell() const { return stream_handle_tell(handle_); }

  bool eof() const { return stream_handle_eof(handle_) != 0; }

private:
  static stream_handle* open_or_throw(const std::string& path,
                                        const s3stream_credentials* creds) {
    stream_handle* handle = creds ? stream_handle_open_with_credentials(path.c_str(), creds)
                                     : stream_handle_open(path.c_str());
    if (!handle) {
      throw std::runtime_error(s3stream_last_error());
    }
    return handle;
  }

  void close() {
    if (handle_) {
      stream_handle_close(handle_);
      handle_ = nullptr;
    }
  }

  stream_handle* handle_;
};

/* std::streambuf over a Stream, so a remote object can be read through the
 * ordinary iostreams machinery -- std::getline, operator>>, and adaptors
 * that wrap a streambuf or an istream (Boost.Iostreams filtering streams,
 * for instance, which is how a gzipped remote file gets decompressed on the
 * fly).
 *
 * Only the get area is implemented; this is a read-only stream. */
class StreamBuf : public std::streambuf {
public:
  explicit StreamBuf(const std::string& path, std::size_t buffer_size = kDefaultBufferSize)
      : stream_(path), buffer_(buffer_size ? buffer_size : kDefaultBufferSize) {
    reset_get_area();
  }

  StreamBuf(const std::string& path, const s3stream_credentials& creds,
            std::size_t buffer_size = kDefaultBufferSize)
      : stream_(path, creds), buffer_(buffer_size ? buffer_size : kDefaultBufferSize) {
    reset_get_area();
  }

  explicit StreamBuf(Stream&& stream, std::size_t buffer_size = kDefaultBufferSize)
      : stream_(std::move(stream)), buffer_(buffer_size ? buffer_size : kDefaultBufferSize) {
    reset_get_area();
  }

  StreamBuf(const StreamBuf&) = delete;
  StreamBuf& operator=(const StreamBuf&) = delete;

protected:
  int_type underflow() override {
    if (gptr() < egptr()) {
      return traits_type::to_int_type(*gptr());
    }

    const int64_t nread = stream_.read(&buffer_[0], buffer_.size());
    if (nread <= 0) {  /* 0 is EOF, -1 an error the istream reports as failbit */
      return traits_type::eof();
    }

    setg(&buffer_[0], &buffer_[0], &buffer_[0] + nread);
    return traits_type::to_int_type(*gptr());
  }

  pos_type seekoff(off_type off, std::ios_base::seekdir dir,
                   std::ios_base::openmode which = std::ios_base::in) override {
    if (!(which & std::ios_base::in)) {
      return pos_type(off_type(-1));
    }

    int whence = SEEK_SET;
    if (dir == std::ios_base::end) {
      whence = SEEK_END;
    } else if (dir == std::ios_base::cur) {
      off += logical_pos();
    }

    if (stream_.seek(off, whence) != 0) {
      return pos_type(off_type(-1));
    }

    reset_get_area();  /* read-ahead is stale after a seek */
    return pos_type(stream_.tell());
  }

  pos_type seekpos(pos_type pos, std::ios_base::openmode which = std::ios_base::in) override {
    return seekoff(off_type(pos), std::ios_base::beg, which);
  }

  std::streamsize showmanyc() override { return egptr() - gptr(); }

private:
  /* The handle sits ahead of the logical stream position by whatever is
   * still unread in the buffer. */
  off_type logical_pos() const {
    return off_type(stream_.tell()) - off_type(egptr() - gptr());
  }

  void reset_get_area() { setg(&buffer_[0], &buffer_[0], &buffer_[0]); }

  Stream stream_;
  std::vector<char> buffer_;
};

/* std::istream over a remote or local path, owning its StreamBuf:
 *
 *   s3stream::IStream in("s3://bucket/phenotypes.tsv");
 *   std::string line;
 *   while (std::getline(in, line)) { ... }
 *
 * Throws std::runtime_error if the path cannot be opened. */
class IStream : public std::istream {
public:
  explicit IStream(const std::string& path, std::size_t buffer_size = kDefaultBufferSize)
      : std::istream(nullptr), buf_(path, buffer_size) {
    rdbuf(&buf_);
  }

  IStream(const std::string& path, const s3stream_credentials& creds,
          std::size_t buffer_size = kDefaultBufferSize)
      : std::istream(nullptr), buf_(path, creds, buffer_size) {
    rdbuf(&buf_);
  }

  explicit IStream(Stream&& stream, std::size_t buffer_size = kDefaultBufferSize)
      : std::istream(nullptr), buf_(std::move(stream), buffer_size) {
    rdbuf(&buf_);
  }

private:
  StreamBuf buf_;
};

}  // namespace s3stream

#endif /* S3STREAM_HPP */
