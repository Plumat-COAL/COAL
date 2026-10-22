// Operating-system side of coal-cpp, as the reference engine:
//   * temp files are registered so a signal handler can remove them: SIGTERM -> exit 143,
//     SIGHUP -> 129, Ctrl-C -> "coal: interrupted", 130; a closed stdout pipe -> exit 141, silent;
//   * non-regular inputs (`-` = standard input, pipes, /dev/stdin, FIFOs) are copied to temp files first;
//   * -o: symlinks written through, FIFOs / devices / hard-linked files written in place, an existing file keeps
//     its permissions, clear errors for symlink loops and missing folders.
#pragma once
#include <string>
#include <utility>
#include <vector>

namespace coal {

// Thrown when stdout is a closed pipe: the caller exits 141 without a message.
struct BrokenPipe {};

// Signal handlers (SIGTERM, SIGHUP, SIGINT; SIGPIPE ignored so writes fail with EPIPE instead).
void install_signal_handlers();
// Temp files removed by the signal handlers (and by unregister_temp's caller as usual).
void register_temp(const std::string& path);
void unregister_temp(const std::string& path);

// "[Errno N] <text>: '<path>'" (Python's OSError message) for the current errno.
std::string os_error(int err, const std::string& path);

// `-` as an input path means standard input (pipeline.py STDIN); messages name it STDIN_NAME. At most one `-` per
// command line (STDIN_TWICE, a usage error, as cli.py).
constexpr const char* STDIN_ARG = "-";
constexpr const char* STDIN_NAME = "standard input";
constexpr const char* STDIN_TWICE =
    "'-' (standard input) can be given only once: as the protein alignment or as one DNA file";
int stdin_count(const std::vector<std::string>& args);

// Inputs that are not regular files (and `-`, standard input) are copied to temp files. Returns the paths to read
// and fills `spooled` with (temp path, name for messages) pairs (removed by remove_spooled).
std::vector<std::string> spool_inputs(const std::vector<std::string>& args,
                                      std::vector<std::pair<std::string, std::string>>& spooled);
void remove_spooled(const std::vector<std::pair<std::string, std::string>>& spooled);

struct OutTarget {
    std::string real;     // the file finally written (symlinks resolved)
    bool existing = false;  // a regular file exists there
    bool direct = false;    // write in place (FIFO, device, a file with other hard links)
};
// Checks -o PATH before any work: throws CoalError for a symlink loop or an unusable folder.
OutTarget resolve_output(const std::string& out);
// "cannot write <out>: <reason in plain words>" (pipeline.py _output_error); out "" = stdout's temp file.
std::string output_error(const std::string& out, int err);

// Python's os.path.realpath (non-strict): the final component's symlinks followed, the folder resolved.
std::string real_path(const std::string& p);
// Whether two paths name the same file (by identity when both exist, else by resolved path).
bool same_file(const std::string& a, const std::string& b);
// Why --report REPORT cannot be written (it names an input or the -o output; a usage error), or "" (report.py
// report_path_conflict).
std::string report_path_conflict(const std::string& report, const std::vector<std::string>& inputs,
                                 const std::string& out);
// Writes a JSON text as the alignment is written: temp file + rename through symlinks, in place for a FIFO, device
// or hard-linked file (report.py write_json). Returns 0 or the errno.
int write_json_file(const std::string& path, const std::string& text);

// Gives `tmp` the permissions `real` has (existing file) before it replaces it.
void copy_mode(const std::string& real, const std::string& tmp);

// fwrite to stdout; throws BrokenPipe on EPIPE (and Windows EINVAL), std::runtime_error on other errors.
void write_stdout(const char* data, size_t n);
// fflush(stdout) with the same error handling.
void flush_stdout();

}  // namespace coal
