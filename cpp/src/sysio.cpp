#include "sysio.hpp"

#include <atomic>
#include <cerrno>
#include <csignal>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <stdexcept>

#include "core.hpp"
#include "seqio.hpp"

#ifdef _WIN32
#include <io.h>
#include <process.h>
#define COAL_UNLINK _unlink
#define COAL_WRITE _write
#define COAL_EXIT _exit
#else
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#define COAL_UNLINK ::unlink
#define COAL_WRITE ::write
#define COAL_EXIT _exit
#endif

namespace fs = std::filesystem;

namespace coal {
namespace {

// Fixed storage, so the signal handler reads it without allocating.
constexpr int SLOTS = 64;
constexpr size_t PATH_MAX_LEN = 4096;
char g_paths[SLOTS][PATH_MAX_LEN];
std::atomic<bool> g_used[SLOTS];

extern "C" void on_signal(int sig) {
    for (int i = 0; i < SLOTS; ++i)
        if (g_used[i].load()) COAL_UNLINK(g_paths[i]);
    if (sig == SIGINT) {
        static const char msg[] = "coal: interrupted\n";
        (void)!COAL_WRITE(2, msg, (unsigned)(sizeof msg - 1));
    }
    COAL_EXIT(128 + sig);
}

std::string unique_name() {
    static std::atomic<unsigned> counter{0};
    std::random_device rd;
    static const char* hex = "0123456789abcdef";
    std::string s;
    unsigned v = rd() ^ (counter++ * 2654435761u);
    for (int i = 0; i < 8; ++i, v >>= 4) s += hex[v & 15];
    return s;
}

}  // namespace

void install_signal_handlers() {
    std::signal(SIGTERM, on_signal);
    std::signal(SIGINT, on_signal);
#ifdef SIGHUP
    std::signal(SIGHUP, on_signal);
#endif
#ifdef SIGPIPE
    std::signal(SIGPIPE, SIG_IGN);  // a closed stdout pipe makes writes fail with EPIPE: exit 141 (DEC-3)
#endif
}

void register_temp(const std::string& path) {
    if (path.empty() || path.size() >= PATH_MAX_LEN) return;
    for (int i = 0; i < SLOTS; ++i)
        if (!g_used[i].load()) {
            std::memcpy(g_paths[i], path.c_str(), path.size() + 1);
            g_used[i].store(true);
            return;
        }
}

void unregister_temp(const std::string& path) {
    for (int i = 0; i < SLOTS; ++i)
        if (g_used[i].load() && path == g_paths[i]) g_used[i].store(false);
}

std::string os_error(int err, const std::string& path) {
    return "[Errno " + std::to_string(err) + "] " + std::strerror(err) + ": " + py_repr(path);
}

std::vector<std::string> spool_inputs(const std::vector<std::string>& args,
                                      std::vector<std::pair<std::string, std::string>>& spooled) {
    std::vector<std::string> out;
    for (const std::string& a : args) {
        std::error_code ec;
        if (fs::is_regular_file(a, ec)) {
            out.push_back(a);
            continue;
        }
        if (fs::is_directory(a, ec)) throw InputError(os_error(EISDIR, a), "file_error");
        std::ifstream src(a, std::ios::binary);
        if (!src) throw InputError(os_error(errno ? errno : ENOENT, a), "file_error");
        std::string base = fs::path(a).filename().string();
        for (char& c : base)
            if (!((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '.' ||
                  c == '_' || c == '-'))
                c = '_';
        if (base.size() > 40) base = base.substr(base.size() - 40);
        std::string tmp = (fs::temp_directory_path() / (".coal_in_" + unique_name() + "_" + base)).string();
        register_temp(tmp);
        spooled.push_back({tmp, a});
        std::ofstream dst(tmp, std::ios::binary | std::ios::trunc);
        if (!dst) throw InputError(os_error(errno ? errno : EACCES, tmp), "file_error");
        std::vector<char> buf(1 << 20);
        while (src) {
            src.read(buf.data(), (std::streamsize)buf.size());
            std::streamsize n = src.gcount();
            if (n > 0) dst.write(buf.data(), n);
        }
        dst.flush();
        if (!dst) throw InputError("could not copy " + a + " to a temporary file (disk full?)", "file_error");
        out.push_back(tmp);
    }
    return out;
}

void remove_spooled(const std::vector<std::pair<std::string, std::string>>& spooled) {
    for (auto& s : spooled) {
        std::error_code ec;
        fs::remove(s.first, ec);
        unregister_temp(s.first);
    }
}

// Python's os.path.realpath (non-strict): follow the final component's symlinks, resolve the folder.
static std::string real_path(const std::string& p) {
    fs::path cur = fs::absolute(fs::path(p));
    for (int i = 0; i < 40; ++i) {
        std::error_code ec;
        if (!fs::is_symlink(cur, ec)) break;
        fs::path t = fs::read_symlink(cur, ec);
        if (ec) break;
        cur = t.is_absolute() ? t : cur.parent_path() / t;
    }
    std::error_code ec;
    fs::path parent = fs::weakly_canonical(cur.parent_path(), ec);
    if (ec) parent = cur.parent_path().lexically_normal();
    return (parent / cur.filename()).string();
}

OutTarget resolve_output(const std::string& out) {
#ifndef _WIN32
    struct stat lst, st;
    if (::lstat(out.c_str(), &lst) == 0 && S_ISLNK(lst.st_mode) && ::stat(out.c_str(), &st) != 0 && errno == ELOOP)
        throw CoalError(out + ": symbolic link loop; cannot write the output there", "file_error");
#endif
    OutTarget t;
    t.real = real_path(out);
    std::error_code ec;
    t.existing = fs::is_regular_file(t.real, ec);
    bool exists = fs::exists(t.real, ec);
    t.direct = exists && (!t.existing || fs::hard_link_count(t.real, ec) > 1);
    if (t.direct) return t;
    // the temp file goes next to the real target: its folder must be usable
    std::string folder = fs::absolute(fs::path(out)).lexically_normal().parent_path().string();
    fs::path dir = fs::path(t.real).parent_path();
#ifndef _WIN32
    if (::stat(dir.c_str(), &st) != 0 || !S_ISDIR(st.st_mode)) {
        int e = errno;
        if (::stat(dir.c_str(), &st) == 0) e = ENOTDIR;
        std::string why = e == ENOENT ? "folder " + folder + " does not exist"
                          : e == ELOOP ? "folder " + folder + " is a symbolic link loop"
                          : e == ENOTDIR ? folder + " is not a folder"
                                         : std::string(std::strerror(e));
        throw CoalError("cannot write " + out + ": " + why, "file_error");
    }
#else
    if (!fs::is_directory(dir, ec)) {
        std::string why = fs::exists(dir, ec) ? folder + " is not a folder" : "folder " + folder + " does not exist";
        throw CoalError("cannot write " + out + ": " + why, "file_error");
    }
#endif
    return t;
}

void copy_mode(const std::string& real, const std::string& tmp) {
    std::error_code ec;
    fs::file_status s = fs::status(real, ec);
    if (!ec) fs::permissions(tmp, s.permissions(), fs::perm_options::replace, ec);
}

void write_stdout(const char* data, size_t n) {
    if (n == 0) return;
    errno = 0;
    size_t w = std::fwrite(data, 1, n, stdout);
    if (w == n) return;
    int e = errno;
#ifdef _WIN32
    if (e == EINVAL || e == EPIPE) throw BrokenPipe();
#else
    if (e == EPIPE) throw BrokenPipe();
#endif
    throw std::runtime_error(std::string("[Errno ") + std::to_string(e) + "] " + std::strerror(e));
}

void flush_stdout() {
    errno = 0;
    if (std::fflush(stdout) == 0) return;
    int e = errno;
#ifdef _WIN32
    if (e == EINVAL || e == EPIPE) throw BrokenPipe();
#else
    if (e == EPIPE) throw BrokenPipe();
#endif
    throw std::runtime_error(std::string("[Errno ") + std::to_string(e) + "] " + std::strerror(e));
}

}  // namespace coal
