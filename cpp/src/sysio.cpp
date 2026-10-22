#include "sysio.hpp"

#include <atomic>
#include <cerrno>
#include <csignal>
#include <cstdio>
#include <cstdlib>
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

// The folder for temporary files, found as Python's tempfile.gettempdir() finds it: the first of TMPDIR, TEMP,
// TMP, then the system's usual folders, then the current folder, in which a file can be created. Cached. So an
// unusable TMPDIR (a missing folder, a file) falls back as the Python engine does instead of failing the run.
std::string temp_dir() {
    static std::string cached;
    if (!cached.empty()) return cached;
    std::vector<std::string> cands;
    for (const char* v : {"TMPDIR", "TEMP", "TMP"}) {
        const char* e = std::getenv(v);
        if (e && *e) cands.push_back(e);
    }
#ifdef _WIN32
    for (const char* d : {"C:\\TEMP", "C:\\TMP", "\\TEMP", "\\TMP"}) cands.push_back(d);
#else
    for (const char* d : {"/tmp", "/var/tmp", "/usr/tmp"}) cands.push_back(d);
#endif
    std::error_code ec;
    cands.push_back(fs::current_path(ec).string());
    for (const std::string& c : cands) {
        if (c.empty() || !fs::is_directory(c, ec)) continue;
        std::string probe = (fs::path(c) / (".coal_probe_" + unique_name())).string();
        bool ok;
        {
            std::ofstream f(probe, std::ios::binary | std::ios::trunc);
            ok = (bool)f;
            if (ok) {
                f << "blat";  // a folder on a full disk is not usable either (Python writes a probe too)
                f.flush();
                ok = (bool)f;
            }
        }
        fs::remove(probe, ec);
        if (!ok) continue;
        cached = fs::absolute(c, ec).string();
        if (cached.size() > 1 && (cached.back() == '/' || cached.back() == '\\')) cached.pop_back();
        return cached;
    }
    throw CoalError("no usable temporary folder: set TMPDIR to a folder you can write to", "file_error");
}

// "No space left on device" -> "no space left on device" (pipeline.py _lower_reason)
std::string lower_reason(int err) {
    std::string s = std::strerror(err);
    if (!s.empty() && s[0] >= 'A' && s[0] <= 'Z') s[0] = char(s[0] + 32);
    return s;
}

// The error for an input that cannot be copied to its temporary file (pipeline.py _copy_error)
std::string copy_error(const std::string& name, int err) {
    return "cannot copy " + name + " to a temporary file in " + temp_dir() + ": " + lower_reason(err);
}

std::string os_error(int err, const std::string& path) {
    return "[Errno " + std::to_string(err) + "] " + std::strerror(err) + ": " + py_repr(path);
}

// A character or block device as an input (`coal /dev/zero n.fa`, `coal - n.fa </dev/urandom`) would be copied to a
// temp file without end, filling the disk: refused before anything is read (pipeline.py _refuse_device). A terminal
// on standard input is read until EOF as before, and the null device is still an empty input.
#ifndef _WIN32  // Windows has no device nodes a path could name; the function would be unused (MSVC C4505 under /WX)
static void refuse_device(int fd, const std::string& name) {
    struct stat st;
    if (fd < 0 || ::fstat(fd, &st) != 0) return;
    if (!S_ISCHR(st.st_mode) && !S_ISBLK(st.st_mode)) return;
    if (::isatty(fd)) return;
    struct stat null_st;
    if (S_ISCHR(st.st_mode) && ::stat("/dev/null", &null_st) == 0 && null_st.st_rdev == st.st_rdev) return;
    throw InputError("cannot read " + name + ": it is a device, not a file", "file_error");
}
#endif

#ifndef _WIN32
// Reads a file descriptor with the small part of std::istream's surface the spool loop uses (read, gcount,
// bool), and closes it at the end of its scope.
class FdSource {
public:
    explicit FdSource(int fd) : fd_(fd) {}
    ~FdSource() { ::close(fd_); }
    FdSource(const FdSource&) = delete;
    FdSource& operator=(const FdSource&) = delete;
    void read(char* buf, std::streamsize n) {
        ssize_t got = ::read(fd_, buf, (size_t)n);
        if (got < 0) { err_ = errno; got = 0; eof_ = true; }
        else if (got == 0) eof_ = true;
        gcount_ = got;
    }
    std::streamsize gcount() const { return gcount_; }
    explicit operator bool() const { return !eof_; }
    int error() const { return err_; }
private:
    int fd_;
    std::streamsize gcount_ = 0;
    bool eof_ = false;
    int err_ = 0;
};
#endif

// `-`: standard input copied, as bytes, to a temp file (pipeline.py _spool_stdin). The same mechanism as a pipe or
// /dev/stdin, so CRLF and BOM handling are those of file input; also on Windows (no /dev/stdin; binary mode).
static std::string spool_stdin(std::vector<std::pair<std::string, std::string>>& spooled) {
    const std::string closed = std::string("cannot read ") + STDIN_NAME + ": it is closed";
#ifdef _WIN32
    if (_fileno(stdin) < 0 || _get_osfhandle(_fileno(stdin)) == -1) throw InputError(closed, "file_error");
    _setmode(_fileno(stdin), _O_BINARY);
#else
    if (::fcntl(0, F_GETFD) == -1) throw InputError(closed, "file_error");
    refuse_device(0, STDIN_NAME);
#endif
    std::string tmp = (fs::path(temp_dir()) / (".coal_in_" + unique_name() + "_stdin")).string();
    register_temp(tmp);
    spooled.push_back({tmp, STDIN_NAME});
    errno = 0;
    std::ofstream dst(tmp, std::ios::binary | std::ios::trunc);
    if (!dst) throw InputError(copy_error(STDIN_NAME, errno ? errno : EACCES), "file_error");
    std::vector<char> buf(1 << 20);
    for (;;) {
        errno = 0;
        size_t n = std::fread(buf.data(), 1, buf.size(), stdin);
        if (n > 0) {
            errno = 0;
            dst.write(buf.data(), (std::streamsize)n);
            if (!dst) throw InputError(copy_error(STDIN_NAME, errno ? errno : ENOSPC), "file_error");
        }
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
    if (!dst) throw InputError(copy_error(STDIN_NAME, errno ? errno : ENOSPC), "file_error");
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
        // One open only: a FIFO's writer is consumed by the first reader (pipeline.py checks the open fd too), so
        // the device check runs on the descriptor the copy will read from, never on a second open of the path.
#ifndef _WIN32
        int fd = ::open(a.c_str(), O_RDONLY | O_NOCTTY);
        if (fd < 0) throw InputError(open_error(a), "file_error");
        try {
            refuse_device(fd, a);
        } catch (...) {
            ::close(fd);
            throw;
        }
        FdSource src(fd);
#else
        std::ifstream src(a, std::ios::binary);
        if (!src) throw InputError(open_error(a), "file_error");
#endif
        std::string base = fs::path(a).filename().string();
        for (char& c : base)
            if (!((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '.' ||
                  c == '_' || c == '-'))
                c = '_';
        if (base.size() > 40) base = base.substr(base.size() - 40);
        std::string tmp = (fs::path(temp_dir()) / (".coal_in_" + unique_name() + "_" + base)).string();
        register_temp(tmp);
        spooled.push_back({tmp, a});
        errno = 0;
        std::ofstream dst(tmp, std::ios::binary | std::ios::trunc);
        if (!dst) throw InputError(copy_error(a, errno ? errno : EACCES), "file_error");
        std::vector<char> buf(1 << 20);
        while (src) {
            src.read(buf.data(), (std::streamsize)buf.size());
            std::streamsize n = src.gcount();
            if (n > 0) {
                errno = 0;
                dst.write(buf.data(), n);
                if (!dst) throw InputError(copy_error(a, errno ? errno : ENOSPC), "file_error");
            }
        }
#ifndef _WIN32
        if (src.error()) throw InputError("cannot read '" + a + "': " + lower_reason(src.error()), "file_error");
#endif
        dst.flush();
        if (!dst) throw InputError(copy_error(a, errno ? errno : ENOSPC), "file_error");
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
std::string real_path(const std::string& p) {
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
    std::error_code ec;
    // An existing non-regular target (a FIFO, /dev/stdout, /dev/fd/1, a device) is written through as it is, before
    // any path resolution: the real path of /dev/fd/1 on a pipe is "pipe:[N]", which no folder holds
    if (fs::exists(out, ec) && !fs::is_directory(out, ec) && !fs::is_regular_file(out, ec)) {
        t.real = out;
        t.direct = true;
        return t;
    }
    t.real = real_path(out);
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
    std::string folder = out.empty() ? temp_dir()
                                     : fs::absolute(fs::path(out)).lexically_normal().parent_path().string();
    std::string why;
    if (err == ENOENT && !fs::is_directory(folder, ec)) why = "folder " + folder + " does not exist";
#ifdef ELOOP
    else if (err == ELOOP) why = "folder " + folder + " is a symbolic link loop";
#endif
    else if (err == ENOTDIR) why = folder + " is not a folder";
    else if (err == EACCES || err == EPERM) why = "permission denied (choose a folder you can write to)";
    else why = lower_reason(err);
    return "cannot write " + (out.empty() ? "a temporary file in " + folder : out) + ": " + why;
}

bool same_file(const std::string& a, const std::string& b) {
    std::error_code ec;
    bool eq = fs::equivalent(a, b, ec);
    if (!ec) return eq;
    return real_path(a) == real_path(b);
}

std::string report_path_conflict(const std::string& report, const std::vector<std::string>& inputs,
                                 const std::string& out) {
    if (report.empty()) return std::string();
    for (const std::string& p : inputs)
        if (p != STDIN_ARG && same_file(report, p)) return "the report path " + report + " is an input file";
    if (!out.empty() && same_file(report, out)) return "the report path " + report + " is the output file";
    return std::string();
}

int write_json_file(const std::string& path, const std::string& text) {
    std::string real = real_path(path);
    std::error_code ec;
    bool exists = fs::exists(real, ec);
    bool regular = exists && fs::is_regular_file(real, ec);
    if (exists && (!regular || fs::hard_link_count(real, ec) > 1)) {
        // a FIFO, device or hard-linked file: written in place, not replaced by a regular file
        errno = 0;
        std::ofstream f(path, std::ios::binary | std::ios::trunc);
        if (!f) return errno ? errno : EIO;
        f << text;
        f.flush();
        return f ? 0 : (errno ? errno : EIO);
    }
    std::string tmp = real + ".tmp" + unique_name();
    int err = 0;
    {
        errno = 0;
        std::ofstream f(tmp, std::ios::binary | std::ios::trunc);
        if (!f) return errno ? errno : EIO;
        f << text;
        f.flush();
        if (!f) err = errno ? errno : EIO;
    }
    if (!err) {
        if (regular) copy_mode(real, tmp);  // keep the existing file's mode
        fs::rename(tmp, real, ec);
        if (ec) err = ec.value();
    }
    if (err) {
        std::error_code rm;
        fs::remove(tmp, rm);
    }
    return err;
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
