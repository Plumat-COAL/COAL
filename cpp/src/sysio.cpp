#include "sysio.hpp"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <thread>

#include "core.hpp"
#include "seqio.hpp"

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#include <process.h>
#include <sys/stat.h>
#include <windows.h>
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
    // an ignored SIGHUP (nohup, setsid) stays ignored, as the Python engine does
    if (std::signal(SIGHUP, SIG_IGN) != SIG_IGN) std::signal(SIGHUP, on_signal);
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

std::string report_folder_missing(const std::string& report) {
    std::error_code ec;
    if (report.empty() || fs::exists(report, ec)) return std::string();
    std::string folder = fs::absolute(fs::path(report)).lexically_normal().parent_path().string();
    if (fs::is_directory(folder, ec)) return std::string();
    if (fs::exists(folder, ec)) return "the report path " + report + ": " + folder + " is not a folder";
    return "the report path " + report + ": folder " + folder + " does not exist";
}

// Test hook: COAL_TEST_HOOKS=1 and COAL_TEST_PAUSE_IN_REPORT=1 stop here while the report's temp file exists (so a
// test can send SIGTERM at that point), as report.py _test_pause_in_report
static void test_pause_in_report() {
    const char* h = std::getenv("COAL_TEST_HOOKS");
    const char* v = std::getenv("COAL_TEST_PAUSE_IN_REPORT");
    if (!h || std::string(h) != "1" || !v || std::string(v) != "1") return;
    std::cerr << "coal: test: paused in report" << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(3600));
}

int write_json_file(const std::string& path, const std::string& text) {
    return write_json_file(path, [&](std::ostream& f) { f << text; });
}

int write_json_file(const std::string& path, const std::function<void(std::ostream&)>& writer) {
    std::error_code ec;
    // An existing non-regular target is written through as it is, before any path resolution (as -o does): the
    // real path of /dev/fd/1 on a pipe is "pipe:[N]", which no folder holds
    if (fs::exists(path, ec) && !fs::is_directory(path, ec) && !fs::is_regular_file(path, ec)) {
        errno = 0;
        std::ofstream f(path, std::ios::binary | std::ios::trunc);
        if (!f) return errno ? errno : EIO;
        writer(f);
        f.flush();
        return f ? 0 : (errno ? errno : EIO);
    }
    std::string real = real_path(path);
    bool exists = fs::exists(real, ec);
    bool regular = exists && fs::is_regular_file(real, ec);
    if (exists && (!regular || fs::hard_link_count(real, ec) > 1)) {
        // a FIFO, device or hard-linked file: written in place, not replaced by a regular file
        errno = 0;
        std::ofstream f(path, std::ios::binary | std::ios::trunc);
        if (!f) return errno ? errno : EIO;
        writer(f);
        f.flush();
        return f ? 0 : (errno ? errno : EIO);
    }
    std::string tmp = real + ".tmp" + unique_name();
    register_temp(tmp);  // SIGTERM while the report is written leaves nothing behind
    int err = 0;
    {
        errno = 0;
        std::ofstream f(tmp, std::ios::binary | std::ios::trunc);
        if (!f) {
            unregister_temp(tmp);
            return errno ? errno : EIO;
        }
        writer(f);
        f.flush();
        if (!f) err = errno ? errno : EIO;
        test_pause_in_report();
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
    unregister_temp(tmp);
    return err;
}

// ------------------------------------------------------------------------------------------------ the output lock
//
// `<real output>.coal-lock`, one line: "<pid> <host> <start> <token>\n", the format of the Python engine
// (src/coal/pipeline.py, where each rule below is explained; docs/pipeline.md "One run per folder"):
//   pid    the process id (all that COAL up to 0.5.0 wrote);
//   host   the machine's name: a lock from another machine counts as held and is never replaced automatically;
//   start  when that process started (Linux: field 22 of /proc/<pid>/stat; "-" elsewhere): a process id that is
//          alive but started at another time is another process (the id was reused);
//   token  random, so that a process recognises exactly its own lock.
// Missing fields (an older lock) are unknown and judged by the pid alone. The lock is created with its content in
// one step (a hard link to a finished temporary file), so it is never seen empty, and a stale lock is replaced by
// one process at a time (`<lock>.takeover`, created exclusively), which removes it only if it still holds the
// content that was judged stale. So runs of either engine that start together onto a stale lock have one holder.

static bool pid_alive(long long pid) {
    if (pid <= 0) return false;
#ifdef _WIN32
    // (no /proc on Windows: a handle to the process and its exit code. Written to mirror pipeline.py's
    // _pid_alive; compiled by CI, its behaviour on Windows is not verified by a test of this lock.)
    HANDLE h = OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION, FALSE, (DWORD)pid);
    if (!h) return false;
    DWORD code = 0;
    bool alive = !GetExitCodeProcess(h, &code) || code == STILL_ACTIVE;
    CloseHandle(h);
    return alive;
#else
    if (::kill((pid_t)pid, 0) == 0) return true;
    return errno == EPERM;
#endif
}

namespace {

long long own_pid() {
#ifdef _WIN32
    return (long long)_getpid();
#else
    return (long long)::getpid();
#endif
}

// The machine's name as pipeline.py's _hostname() writes it: printable ASCII as it is, other ASCII (blanks,
// control characters) as '_', every other character as '?'; "-" when there is none.
std::string host_name() {
    std::string out;
#ifdef _WIN32
    wchar_t buf[256];
    DWORD n = 256;
    if (GetComputerNameExW(ComputerNamePhysicalDnsHostname, buf, &n))   // what Python's socket.gethostname() asks
        for (DWORD k = 0; k < n; ++k)
            out += buf[k] >= 128 ? '?' : (buf[k] > 32 && buf[k] < 127 ? (char)buf[k] : '_');
#else
    char buf[256];
    if (::gethostname(buf, sizeof buf) == 0) {
        buf[sizeof buf - 1] = '\0';
        for (const char* q = buf; *q; ++q) {
            const unsigned char ch = (unsigned char)*q;
            if (ch >= 0x80) {
                if (ch >= 0xC0) out += '?';   // one '?' per character: UTF-8 continuation bytes are skipped
            } else {
                out += (ch > 32 && ch < 127) ? (char)ch : '_';
            }
        }
    }
#endif
    return out.empty() ? "-" : out;
}

// The start time of a process as the system counts it, or "-" (unknown: not Linux, or no such process). The same
// text as pipeline.py's _proc_start: the 20th field after the last ')' of /proc/<pid>/stat.
std::string proc_start(long long pid) {
#ifdef __linux__
    std::ifstream in("/proc/" + std::to_string(pid) + "/stat", std::ios::binary);
    if (!in) return "-";
    std::string text(4096, '\0');
    in.read(&text[0], (std::streamsize)text.size());
    text.resize((size_t)in.gcount());
    const size_t close = text.rfind(')');
    if (close == std::string::npos) return "-";
    std::istringstream fields(text.substr(close + 1));
    std::string f;
    for (int k = 0; k <= 19; ++k)
        if (!(fields >> f)) return "-";
    return f;
#else
    (void)pid;
    return "-";
#endif
}

struct LockInfo {
    bool exists = false;     // there is a file at that path
    std::string raw;         // its bytes (at most 4096)
    bool has_pid = false;    // the first field is a number
    long long pid = -1;
    std::string host = "-", start = "-";
};

LockInfo read_lock(const std::string& path) {
    LockInfo li;
    errno = 0;
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        std::error_code ec;
        li.exists = fs::exists(fs::symlink_status(path, ec));   // there, but unreadable: not a lock we can judge
        return li;
    }
    li.exists = true;
    li.raw.assign(4096, '\0');
    in.read(&li.raw[0], (std::streamsize)li.raw.size());
    li.raw.resize((size_t)in.gcount());
    std::istringstream fields(li.raw);
    std::string f0, f1, f2;
    if (!(fields >> f0)) return li;
    size_t k = (f0[0] == '-' || f0[0] == '+') ? 1 : 0;
    if (k == f0.size() || f0.size() > 18) return li;
    for (size_t q = k; q < f0.size(); ++q)
        if (f0[q] < '0' || f0[q] > '9') return li;
    li.pid = std::atoll(f0.c_str());
    li.has_pid = true;
    if (fields >> f1) li.host = f1;
    if (fields >> f2) li.start = f2;
    return li;
}

enum class Holder { nobody, other_host, process };

// Who holds a lock with these fields (pipeline.py _lock_holder).
Holder lock_holder(const LockInfo& li) {
    if (li.host != "-" && li.host != host_name()) return Holder::other_host;
    if (!pid_alive(li.pid)) return Holder::nobody;
    if (li.start != "-") {
        const std::string now = proc_start(li.pid);
        if (now != "-" && now != li.start) return Holder::nobody;   // the id is in use again, by another process
    }
    if (li.pid == own_pid()) return Holder::nobody;   // our id, but this process holds no lock yet: a leftover
    return Holder::process;
}

bool write_all(int fd, const std::string& data) {
    size_t done = 0;
    while (done < data.size()) {
        auto n = COAL_WRITE(fd, data.data() + done, (unsigned)(data.size() - done));
        if (n <= 0) return false;
        done += (size_t)n;
    }
    return true;
}

int open_excl(const std::string& path) {
#ifdef _WIN32
    return ::_open(path.c_str(), _O_WRONLY | _O_CREAT | _O_EXCL | _O_BINARY, _S_IREAD | _S_IWRITE);
#else
    return ::open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0644);
#endif
}

void close_fd(int fd) {
#ifdef _WIN32
    ::_close(fd);
#else
    ::close(fd);
#endif
}

// Create the file `path` holding `raw` in one step: a hard link to a finished temporary file, so that the lock is
// never seen empty; where links are not possible, exclusively and then written. 1: created; 0: it exists;
// < 0: -errno of what went wrong (a folder that cannot be written).
int create_lock(const std::string& path, const std::string& raw) {
    const std::string tmp = path + "." + std::to_string(own_pid()) + "." + unique_name() + ".tmp";
    register_temp(tmp);
    int fd = open_excl(tmp);
    if (fd < 0 && errno != EEXIST) {
        const int e = errno;
        unregister_temp(tmp);
        return -e;
    }
    if (fd >= 0) {
        const bool written = write_all(fd, raw);
        close_fd(fd);
        int linked = -1;   // 1 linked, 0 the lock exists, -1 no hard links here (or the write failed)
        if (written) {
#ifdef _WIN32
            if (CreateHardLinkA(path.c_str(), tmp.c_str(), nullptr)) linked = 1;
            else if (GetLastError() == ERROR_ALREADY_EXISTS || GetLastError() == ERROR_FILE_EXISTS) linked = 0;
#else
            if (::link(tmp.c_str(), path.c_str()) == 0) linked = 1;
            else if (errno == EEXIST) linked = 0;
#endif
        }
        COAL_UNLINK(tmp.c_str());
        unregister_temp(tmp);
        if (linked >= 0) return linked;
    } else {
        unregister_temp(tmp);
    }
    fd = open_excl(path);
    if (fd < 0) return errno == EEXIST ? 0 : -errno;
    const bool written = write_all(fd, raw);
    const int e = errno;
    close_fd(fd);
    if (!written) {
        COAL_UNLINK(path.c_str());
        return -(e ? e : EIO);
    }
    return 1;
}

// Remove the lock at `path` if it still holds `stale` (the content that was judged stale). One process at a time
// does this for a path: it first creates `<lock>.takeover` exclusively (pipeline.py _replace_stale_lock). Returns
// 0 or the errno of a guard that cannot be created.
int replace_stale_lock(const std::string& path, const std::string& stale) {
    const std::string guard = path + ".takeover";
    const std::string mine = std::to_string(own_pid()) + " " + host_name() + " " + proc_start(own_pid()) + "\n";
    const int made = create_lock(guard, mine);
    if (made < 0) return -made;
    if (made == 0) {
        // another run is replacing the lock right now (it takes microseconds). A guard whose process is gone, or
        // one that has been there for 30 s, is a leftover: moved aside under a unique name (one mover wins)
        const LockInfo g = read_lock(guard);
        std::error_code ec;
        bool old = false;
        const auto mtime = fs::last_write_time(guard, ec);
        if (!ec) old = fs::file_time_type::clock::now() - mtime > std::chrono::seconds(30);
        if (g.exists && (old || (g.has_pid && lock_holder(g) == Holder::nobody))) {
            const std::string aside = guard + "." + std::to_string(own_pid()) + "." + unique_name();
            fs::rename(guard, aside, ec);
            if (!ec) fs::remove(aside, ec);
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        return 0;
    }
    register_temp(guard);
    const LockInfo now = read_lock(path);
    std::error_code ec;
    if (now.exists && now.raw == stale) fs::remove(path, ec);   // still the stale one: nobody else replaced it
    fs::remove(guard, ec);
    unregister_temp(guard);
    return 0;
}

std::string g_lock_raw;   // what this process wrote into the lock it holds (one -o per run)

}  // namespace

std::string take_output_lock(const std::string& real, const std::string& name) {
    const std::string path = real + ".coal-lock";
    const std::string raw = std::to_string(own_pid()) + " " + host_name() + " " + proc_start(own_pid()) + " " +
                            unique_name() + unique_name() + "\n";
    // Test hook (tests/test_output_lock.py): start at this time (milliseconds since the epoch), so that several
    // runs reach the lock in the same instant.
    if (const char* at = std::getenv("COAL_TEST_LOCK_START_AT")) {
        const auto when = std::chrono::system_clock::time_point(std::chrono::milliseconds(std::atoll(at)));
        if (when - std::chrono::system_clock::now() < std::chrono::seconds(60))
            while (std::chrono::system_clock::now() < when) {}
    }
    double waited = 0.0;
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(30);
    for (;;) {
        const int made = create_lock(path, raw);
        if (made < 0) throw CoalError(output_error(name, -made), "file_error");
        if (made == 1) {
            const LockInfo got = read_lock(path);
            if (!(got.exists && got.raw == raw)) continue;   // replaced under us by a COAL that predates this
            register_temp(path);
            g_lock_raw = raw;
            // Test hook (tests/test_output_lock.py): hold the lock until this file exists, then check that it is
            // still ours, so that a test can count the runs that hold one output at the same moment.
            if (const char* until = std::getenv("COAL_TEST_HOLD_LOCK_UNTIL")) {
                std::error_code ec;
                for (int k = 0; k < 12000 && !fs::exists(until, ec); ++k)
                    std::this_thread::sleep_for(std::chrono::milliseconds(5));
                if (read_lock(path).raw != raw)
                    throw CoalError("test hook: the lock " + path + " was replaced under its holder", "file_error");
            }
            return path;
        }
        const LockInfo cur = read_lock(path);
        if (!cur.exists) continue;   // gone in the meantime: try again
        if (!cur.has_pid) {
            // Empty (or not a lock): a COAL up to 0.5.0 created it an instant ago and has not written its id yet,
            // so it is not stale; only one that stays empty for 2 s is
            if (waited < 2.0) {
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
                waited += 0.05;
                continue;
            }
        } else {
            const Holder h = lock_holder(cur);
            if (h == Holder::other_host)
                throw CoalError("another COAL run (process " + std::to_string(cur.pid) + " on the machine " + cur.host +
                                    ") is writing " + name + ", or was when it last ran: a process on another machine "
                                    "cannot be checked from here. Wait for it to finish, or write to another file (if "
                                    "that run is gone, remove " + path + ")",
                                "file_error");
            if (h == Holder::process)
                throw CoalError("another COAL run (process " + std::to_string(cur.pid) + ") is writing " + name +
                                    ": wait for it to finish, or write to another file (if that process is gone, remove " +
                                    path + ")",
                                "file_error");
        }
        if (std::chrono::steady_clock::now() > deadline)
            throw CoalError("cannot write " + name + ": " + path + " is in the way and cannot be removed", "file_error");
        // Test hook: wait between judging the lock stale and replacing it (the window of the race)
        if (const char* ms = std::getenv("COAL_TEST_LOCK_BEFORE_TAKEOVER_MS"))
            std::this_thread::sleep_for(std::chrono::milliseconds(std::min(std::atoll(ms), 60000LL)));
        const int err = replace_stale_lock(path, cur.raw);   // a lock whose process is gone (or one empty for 2 s)
        if (err) throw CoalError(output_error(name, err), "file_error");
    }
}

// Removes the lock this process took. A file at that path that is not ours any more (someone removed the lock by
// hand and another run took it) is left alone.
void release_output_lock(const std::string& path) {
    if (path.empty()) return;
    const LockInfo cur = read_lock(path);
    unregister_temp(path);
    if (!cur.exists || (!g_lock_raw.empty() && cur.raw != g_lock_raw && cur.has_pid)) return;
    std::error_code ec;
    fs::remove(path, ec);
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

StdoutBuf::StdoutBuf() : buf_(1 << 16) { setp(buf_.data(), buf_.data() + buf_.size()); }

StdoutBuf::~StdoutBuf() = default;

bool StdoutBuf::send(const char* data, size_t n) {
    if (broken_ || !error_.empty()) return false;
    try {
        write_stdout(data, n);
        return true;
    } catch (const BrokenPipe&) {
        broken_ = true;
    } catch (const std::exception& e) {
        error_ = e.what();
    }
    return false;
}

bool StdoutBuf::flush_buffer() {
    size_t n = (size_t)(pptr() - pbase());
    setp(buf_.data(), buf_.data() + buf_.size());
    return n == 0 || send(buf_.data(), n);
}

int StdoutBuf::overflow(int ch) {
    if (!flush_buffer()) return traits_type::eof();
    if (ch != traits_type::eof()) {
        *pptr() = (char)ch;
        pbump(1);
    }
    return traits_type::not_eof(ch);
}

std::streamsize StdoutBuf::xsputn(const char* s, std::streamsize n) {
    if (n <= 0) return 0;
    if ((size_t)n <= (size_t)(epptr() - pptr())) {
        std::memcpy(pptr(), s, (size_t)n);
        pbump((int)n);
        return n;
    }
    if (!flush_buffer()) return 0;
    if ((size_t)n >= buf_.size()) return send(s, (size_t)n) ? n : 0;  // large pieces go straight out
    std::memcpy(pptr(), s, (size_t)n);
    pbump((int)n);
    return n;
}

int StdoutBuf::sync() { return flush_buffer() ? 0 : -1; }

void StdoutBuf::finish() {
    if (flush_buffer()) {
        try {
            flush_stdout();
        } catch (const BrokenPipe&) {
            broken_ = true;
        } catch (const std::exception& e) {
            error_ = e.what();
        }
    }
    if (broken_) throw BrokenPipe();
    if (!error_.empty()) throw std::runtime_error(error_);
}

}  // namespace coal
