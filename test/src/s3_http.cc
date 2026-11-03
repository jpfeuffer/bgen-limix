#include "bgen/s3stream.h"
#include "cass.h"

#include <arpa/inet.h>
#include <poll.h>
#include <signal.h>
#include <sys/socket.h>
#include <unistd.h>

#include <atomic>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

class HttpServer {
 public:
  std::atomic<int> connections{0};
  std::atomic<int> metadata_requests{0};
  std::atomic<int> generation{1};
  std::atomic<bool> metadata_failure{false};
  std::atomic<bool> metadata_expired{false};

  HttpServer() {
    listener_ = socket(AF_INET, SOCK_STREAM, 0);
    sockaddr_in address{};
    address.sin_family = AF_INET;
    address.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
    if ((listener_ < 0) ||
        (bind(listener_, reinterpret_cast<sockaddr*>(&address), sizeof(address)) < 0) ||
        (listen(listener_, 16) < 0)) {
      throw std::runtime_error("cannot start loopback HTTP server");
    }
    socklen_t address_size = sizeof(address);
    getsockname(listener_, reinterpret_cast<sockaddr*>(&address), &address_size);
    url = "http://127.0.0.1:" + std::to_string(ntohs(address.sin_port));
    worker_ = std::thread([this] { Run(); });
  }

  ~HttpServer() {
    stopping_ = true;
    worker_.join();
    close(listener_);
  }

  bool SawAuthorization(const std::string& value) {
    std::lock_guard<std::mutex> lock(requests_mutex_);
    for (const auto& request : requests_) {
      if (request.find(value) != std::string::npos) {
        return true;
      }
    }
    return false;
  }

  std::string url;

 private:
  int listener_;
  std::atomic<bool> stopping_{false};
  std::thread worker_;
  std::mutex requests_mutex_;
  std::vector<std::string> requests_;

  std::string Response(const std::string& request) {
    {
      std::lock_guard<std::mutex> lock(requests_mutex_);
      requests_.push_back(request);
    }
    const bool head = request.compare(0, 5, "HEAD ") == 0;
    int status = 200;
    std::string body;
    std::string headers;
    size_t length = 16384;
    if (request.find(" /metadata") != std::string::npos) {
      ++metadata_requests;
      if (metadata_failure) {
        status = 404;
      } else {
        body = "{\"AccessKeyId\":\"CACHE" + std::to_string(generation.load()) +
               "\",\"SecretAccessKey\":\"secret\",\"Token\":\"token" +
               std::to_string(generation.load()) + "\",\"Expiration\":\"" +
               (metadata_expired ? "2000-01-01" : "2099-01-01") +
               "T00:00:00Z\"}";
      }
      length = body.size();
    } else if (request.find("/missing") != std::string::npos) {
      status = 404;
      length = 0;
    } else if (request.find("/denied") != std::string::npos) {
      status = 403;
      length = 0;
    } else if ((request.find("/anonymous") != std::string::npos) &&
               ((request.find("Authorization:") != std::string::npos) ||
                (request.find("x-amz-security-token:") != std::string::npos))) {
      status = 403;
      length = 0;
    } else if (!head && (request.find("/refresh") != std::string::npos) &&
               (request.find("Credential=CACHE1/") != std::string::npos)) {
      status = 403;
      body = "<Error><Code>ExpiredToken</Code></Error>";
      length = body.size();
    } else {
      headers = "ETag: \"fixture\"\r\n";
      const size_t range = request.find("\r\nRange: bytes=");
      if (!head && (range != std::string::npos)) {
        size_t first = 0;
        size_t last = 0;
        sscanf(request.c_str() + range, "\r\nRange: bytes=%zu-%zu", &first, &last);
        if ((first >= length) || (last < first) || (last >= length)) {
          status = 416;
          length = 0;
        } else {
          status = 206;
          headers += "Content-Range: bytes " + std::to_string(first) + "-" +
                     std::to_string(last) + "/" + std::to_string(length) + "\r\n";
          length = last - first + 1;
          body.assign(length, 'x');
        }
      }
    }
    return "HTTP/1.1 " + std::to_string(status) + " Test\r\n" + headers +
           "Content-Length: " + std::to_string(length) + "\r\n\r\n" +
           (head ? std::string() : body);
  }

  void Run() {
    std::vector<int> clients;
    std::vector<std::string> pending;
    while (!stopping_) {
      std::vector<pollfd> watches{{listener_, POLLIN, 0}};
      for (int client : clients) {
        watches.push_back({client, POLLIN, 0});
      }
      if (poll(watches.data(), watches.size(), 100) <= 0) {
        continue;
      }
      for (size_t client_index = clients.size(); client_index > 0;) {
        --client_index;
        if (!watches[client_index + 1].revents) {
          continue;
        }
        char buffer[4096];
        const ssize_t count = recv(clients[client_index], buffer, sizeof(buffer), 0);
        if (count <= 0) {
          close(clients[client_index]);
          clients.erase(clients.begin() + client_index);
          pending.erase(pending.begin() + client_index);
          continue;
        }
        pending[client_index].append(buffer, static_cast<size_t>(count));
        const size_t end = pending[client_index].find("\r\n\r\n");
        if (end == std::string::npos) {
          continue;
        }
        const std::string response = Response(pending[client_index].substr(0, end + 4));
        pending[client_index].erase(0, end + 4);
        size_t sent = 0;
        while (sent < response.size()) {
          const ssize_t count_sent = send(clients[client_index], response.data() + sent,
                                          response.size() - sent, 0);
          if (count_sent <= 0) {
            break;
          }
          sent += static_cast<size_t>(count_sent);
        }
      }
      if (watches[0].revents & POLLIN) {
        const int client = accept(listener_, nullptr, nullptr);
        if (client >= 0) {
          clients.push_back(client);
          pending.emplace_back();
          ++connections;
        }
      }
    }
    for (int client : clients) {
      close(client);
    }
  }
};

int main() {
  signal(SIGPIPE, SIG_IGN);
  setenv("NO_PROXY", "127.0.0.1", 1);
  setenv("no_proxy", "127.0.0.1", 1);
  cass_cond(s3stream_init() == 0);
  errno = 0;
  cass_cond(s3stream_open(nullptr) == nullptr);
  cass_cond(errno == EINVAL);
  cass_cond(stream_handle_open("s3://invalid") == nullptr);
  cass_cond(errno == EINVAL);
  cass_cond(stream_handle_open("/nonexistent/bgen-fixture") == nullptr);
  cass_cond(errno == ENOENT);

  HttpServer server;
  s3stream_credentials failure_creds{};
  failure_creds.endpoint_url = server.url.c_str();
  failure_creds.region = "us-east-1";
  failure_creds.no_sign_request = 1;
  failure_creds.force_path_style = 1;
  for (const auto& failure : std::vector<std::pair<std::string, int>>{
           {"missing", ENOENT}, {"denied", EACCES}}) {
    const std::string url = server.url + "/" + failure.first;
    const std::string uri = "s3://bucket/" + failure.first;
    errno = 0;
    cass_cond(s3stream_open(url.c_str()) == nullptr);
    cass_equal_int(errno, failure.second);
    cass_cond(std::string(s3stream_last_error()).find("HTTP ") != std::string::npos);
    errno = 0;
    cass_cond(stream_handle_open(url.c_str()) == nullptr);
    cass_equal_int(errno, failure.second);
    cass_cond(std::string(s3stream_last_error()).find("HTTP ") != std::string::npos);
    errno = 0;
    cass_cond(s3stream_open_with_credentials(uri.c_str(), &failure_creds) == nullptr);
    cass_equal_int(errno, failure.second);
    cass_cond(std::string(s3stream_last_error()).find("HTTP ") != std::string::npos);
    errno = 0;
    cass_cond(stream_handle_open_with_credentials(uri.c_str(), &failure_creds) == nullptr);
    cass_equal_int(errno, failure.second);
    cass_cond(std::string(s3stream_last_error()).find("HTTP ") != std::string::npos);
  }

  const int before = server.connections;
  stream_handle* handle = stream_handle_open((server.url + "/object").c_str());
  cass_cond(handle != nullptr);
  if (handle) {
    cass_cond(stream_handle_set_chunk_size(handle, 4096) == 0);
    char buffer[4096];
    for (int chunk = 0; chunk < 4; ++chunk) {
      cass_cond(stream_handle_read(handle, buffer, sizeof(buffer)) == sizeof(buffer));
      cass_cond(std::string(buffer, sizeof(buffer)) == std::string(sizeof(buffer), 'x'));
    }
    cass_cond(stream_handle_read(handle, buffer, 1) == 0);
    cass_cond(stream_handle_seek(handle, 0, SEEK_SET) == 0);
    cass_cond(stream_handle_read(handle, buffer, 1) == 1);
    cass_cond(server.connections == before + 1);
    stream_handle_close(handle);
  }

  unsetenv("AWS_ACCESS_KEY_ID");
  unsetenv("AWS_SECRET_ACCESS_KEY");
  unsetenv("AWS_SESSION_TOKEN");
  unsetenv("AWS_CONTAINER_CREDENTIALS_RELATIVE_URI");
  unsetenv("AWS_NO_SIGN_REQUEST");
  setenv("AWS_EC2_METADATA_DISABLED", "true", 1);
  setenv("AWS_SHARED_CREDENTIALS_FILE", "/nonexistent/bgen-credentials", 1);
  setenv("AWS_CONFIG_FILE", "/nonexistent/bgen-config", 1);
  setenv("AWS_ENDPOINT_URL_S3", server.url.c_str(), 1);
  setenv("AWS_REGION", "us-east-1", 1);
  setenv("AWS_CONTAINER_CREDENTIALS_FULL_URI", (server.url + "/metadata").c_str(), 1);

  const int metadata_before = server.metadata_requests;
  for (int open_index = 0; open_index < 3; ++open_index) {
    handle = stream_handle_open("s3://bucket/object");
    cass_cond(handle != nullptr);
    stream_handle_close(handle);
  }
  cass_cond(server.metadata_requests == metadata_before + 1);
  cass_cond(server.SawAuthorization("Credential=CACHE1/"));

  handle = stream_handle_open("s3://bucket/refresh");
  cass_cond(handle != nullptr);
  if (handle) {
    server.generation = 2;
    char buffer;
    cass_cond(stream_handle_read(handle, &buffer, 1) == 1);
    cass_cond(buffer == 'x');
    cass_cond(server.metadata_requests == metadata_before + 2);
    cass_cond(server.SawAuthorization("Credential=CACHE2/"));
    stream_handle_close(handle);
  }

  s3stream_credentials explicit_creds{};
  explicit_creds.access_key_id = "EXPLICIT";
  explicit_creds.secret_access_key = "explicit-secret";
  handle = stream_handle_open_with_credentials("s3://bucket/object", &explicit_creds);
  cass_cond(handle != nullptr);
  stream_handle_close(handle);
  cass_cond(server.SawAuthorization("Credential=EXPLICIT/"));
  cass_cond(server.metadata_requests == metadata_before + 2);

  handle = stream_handle_open("s3://bucket/anonymous");
  cass_cond(handle != nullptr);
  if (handle) {
    char buffer;
    cass_cond(stream_handle_read(handle, &buffer, 1) == 1);
    stream_handle_close(handle);
  }

  setenv("AWS_ACCESS_KEY_ID", "ENVKEY", 1);
  setenv("AWS_SECRET_ACCESS_KEY", "env-secret", 1);
  handle = stream_handle_open("s3://bucket/object");
  cass_cond(handle != nullptr);
  stream_handle_close(handle);
  cass_cond(server.SawAuthorization("Credential=ENVKEY/"));
  cass_cond(server.metadata_requests == metadata_before + 2);
  unsetenv("AWS_ACCESS_KEY_ID");
  unsetenv("AWS_SECRET_ACCESS_KEY");

  server.metadata_failure = true;
  setenv("AWS_PROFILE", "negative-cache", 1);
  for (int open_index = 0; open_index < 2; ++open_index) {
    handle = stream_handle_open("s3://bucket/object");
    cass_cond(handle != nullptr);
    stream_handle_close(handle);
  }
  cass_cond(server.metadata_requests == metadata_before + 3);

  server.metadata_failure = false;
  server.metadata_expired = true;
  setenv("AWS_PROFILE", "expired-cache", 1);
  handle = stream_handle_open("s3://bucket/object");
  cass_cond(handle != nullptr);
  stream_handle_close(handle);
  const int expired_before = server.metadata_requests;
  server.metadata_expired = false;
  server.generation = 3;
  handle = stream_handle_open("s3://bucket/object");
  cass_cond(handle != nullptr);
  stream_handle_close(handle);
  cass_cond(server.metadata_requests == expired_before + 1);
  cass_cond(server.SawAuthorization("Credential=CACHE3/"));

  setenv("AWS_PROFILE", "concurrent-cache", 1);
  const int concurrent_before = server.metadata_requests;
  std::atomic<int> opened{0};
  std::vector<std::thread> readers;
  for (int reader_index = 0; reader_index < 4; ++reader_index) {
    readers.emplace_back([&opened] {
      stream_handle* reader = stream_handle_open("s3://bucket/object");
      if (reader) {
        ++opened;
        stream_handle_close(reader);
      }
    });
  }
  for (auto& reader : readers) {
    reader.join();
  }
  cass_cond(opened == 4);
  cass_cond(server.metadata_requests == concurrent_before + 1);
  return cass_status();
}