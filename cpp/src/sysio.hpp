// Operating-system side of coal-cpp, as the reference engine (overnight/integration, engine-fixes):
//   * temp files are registered so a signal handler can remove them (ENG-f, R2-E8): SIGTERM -> exit 143,
//     SIGHUP -> 129, Ctrl-C -> "coal: interrupted", 130; a closed stdout pipe -> exit 141, silent (DEC-3, R2-D1b);
//   * non-regular inputs (pipes, /dev/stdin, FIFOs) are copied to temp files first (ENG-b, R2-E7);
//   * -o: symlinks written through, FIFOs / devices / hard-linked files written in place, an existing file keeps
//     its permissions, clear errors for symlink loops and missing folders (ENG-a, R2-E5, R3-E2).
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

// Inputs that are not regular files are copied to temp files. Returns the paths to read and fills
// `spooled` with (temp path, original argument) pairs (removed by remove_spooled).
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
// Gives `tmp` the permissions `real` has (existing file) before it replaces it.
void copy_mode(const std::string& real, const std::string& tmp);

// fwrite to stdout; throws BrokenPipe on EPIPE (and Windows EINVAL), std::runtime_error on other errors.
void write_stdout(const char* data, size_t n);
// fflush(stdout) with the same error handling.
void flush_stdout();

}  // namespace coal
