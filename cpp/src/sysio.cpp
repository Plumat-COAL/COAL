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
#include <fcntl.h>
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
    std::signal(SIGPIPE, SIG_IGN);  // a closed stdout pipe makes writes fail with EPIPE: exit 141
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

// `-`: standard input copied, as bytes, to a temp file (pipeline.py _spool_stdin). The same mechanism as a pipe or
// /dev/stdin, so CRLF and BOM handling are those of file input; also on Windows (no /dev/stdin; binary mode).
static std::string spool_stdin(std::vector<std::pair<std::string, std::string>>& spooled) {
    const std::string closed = std::string("cannot read ") + STDIN_NAME + ": it is closed";
#ifdef _WIN32
    if (_fileno(stdin) < 0 || _get_osfhandle(_fileno(stdin)) == -1) throw InputError(closed, "file_error");
    _setmode(_fileno(stdin), _O_BINARY);
#else
    if (::fcntl(0, F_GETFD) == -1) throw InputError(closed, "file_error");
#endif
    std::string tmp = (fs::temp_directory_path() / (".coal_in_" + unique_name() + "_stdin")).string();
    register_temp(tmp);
    spooled.push_back({tmp, STDIN_NAME});
    std::ofstream dst(tmp, std::ios::binary | std::ios::trunc);
    if (!dst) throw InputError(os_error(errno ? errno : EACCES, tmp), "file_error");
    std::vector<char> buf(1 << 20);
    for (;;) {
        errno = 0;
        size_t n = std::fread(buf.data(), 1, buf.size(), stdin);
        if (n > 0) dst.write(buf.data(), (std::streamsize)n);
        if (n < buf.size()) {
            if (std::ferror(stdin)) {
                int e = errno ? errno : EIO;
                throw InputError(std::string("cannot read ") + STDIN_NAME + ": " +
                                     (e == EBADF ? std::string("it is closed") : std::string(std::strerror(e))),
                                 "file_error");
            }
            break;
        }
    }
    dst.flush();
    if (!dst) throw InputError(std::string("could not copy ") + STDIN_NAME + " to a temporary file (disk full?)",
                               "file_error");
    return tmp;
}

int stdin_count(const std::vector<std::string>& args) {
    int n = 0;
    for (const std::string& a : args) n += a == STDIN_ARG;
    return n;
}

std::vector<std::string> spool_inputs(const std::vector<std::string>& args,
                                      std::vector<std::pair<std::string, std::string>>& spooled) {
    std::vector<std::string> out;
    for (const std::string& a : args) {
        if (a == STDIN_ARG) {
            out.push_back(spool_stdin(spooled));
            continue;
        }
        std::error_code ec;
        if (fs::is_regular_file(a, ec)) {
            out.push_back(a);
            continue;
        }
        // "cannot open '<path>': <reason>"
        if (fs::is_directory(a, ec)) throw InputError(open_error(a), "file_error");
        std::ifstream src(a, std::ios::binary);
        if (!src) throw InputError(open_error(a), "file_error");
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
    // -o naming a folder is refused here, before any work
    if (fs::is_directory(t.real, ec))
        throw CoalError("cannot write " + out + ": it is a folder, not a file (give a file name with -o)", "file_error");
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

// As pipeline.py _output_error: ENOENT in a folder that exists is not "does not exist", permission denied says
// what to do, no -o names the temp folder.
std::string output_error(const std::string& out, int err) {
    std::error_code ec;
    std::string folder = out.empty() ? fs::temp_directory_path(ec).string()
                                     : fs::absolute(fs::path(out)).lexically_normal().parent_path().string();
    std::string why;
    if (err == ENOENT && !fs::is_directory(folder, ec)) why = "folder " + folder + " does not exist";
#ifdef ELOOP
    else if (err == ELOOP) why = "folder " + folder + " is a symbolic link loop";
#endif
    else if (err == ENOTDIR) why = folder + " is not a folder";
    else if (err == EACCES || err == EPERM) why = "permission denied (choose a folder you can write to)";
    else why = std::strerror(err);
    return "cannot write " + (out.empty() ? "a temporary file in " + folder : out) + ": " + why;
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
