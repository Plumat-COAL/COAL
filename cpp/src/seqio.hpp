// Readers matching src/coal/seqio.py: FASTA / CLUSTAL protein alignments and an
// on-disk index of DNA FASTA files (IDs + byte offsets only).
#pragma once
#include <cstdint>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace coal {

// code: "input_format" (unreadable/malformed input) or "file_error" (cannot open), as in the Python engine.
struct InputError : std::runtime_error {
    std::string code = "input_format";
    explicit InputError(const std::string& msg, std::string c = "input_format")
        : std::runtime_error(msg), code(std::move(c)) {}
};

// Python's repr() of an ID string ('abc', or "it's").
std::string py_repr(const std::string& s);

// "cannot open '<path>': no such file" (or: permission denied, it is a folder, or strerror), as
// open_error_message in src/coal/seqio.py, for an input that cannot be opened.
std::string open_error(const std::string& path);

// The reader messages (src/coal/seqio.py NON_ASCII, NON_ASCII_ID, EMPTY_HEADER, BEFORE_HEADER_*).
extern const char* const NON_ASCII;
extern const char* const NON_ASCII_ID;
extern const char* const EMPTY_HEADER;
extern const char* const BEFORE_HEADER_PROTEIN;
extern const char* const BEFORE_HEADER_DNA;
std::string duplicate_message(const std::string& where, const std::string& rid);
// Python str.split() whitespace (ASCII).
bool py_space(unsigned char c);

// Reads lines with Python "universal newlines" (\n, \r\n or \r end a line).
class LineReader {
public:
    explicit LineReader(const std::string& path);
    bool next(std::string& line);  // line without its terminator
    long long lineno() const { return lineno_; }

private:
    bool fill();
    std::string path_;
    std::ifstream in_;
    std::vector<char> buf_;
    size_t pos_ = 0, len_ = 0;
    long long lineno_ = 0;
    bool eof_ = false;
};

// Streaming protein alignment reader: FASTA (one record at a time) or CLUSTAL
// (read whole, as the Python version does).
class ProteinReader {
public:
    explicit ProteinReader(const std::string& path);
    bool next(std::string& id, std::string& row);  // raw row (not yet prepared)
    bool is_clustal() const { return clustal_; }

private:
    std::string path_;
    bool clustal_ = false;
    std::vector<std::pair<std::string, std::string>> clustal_recs_;
    size_t ci_ = 0;
    std::unique_ptr<LineReader> lr_;
    // The '>' line that ended the previous record, kept raw: like iter_fasta in the
    // reference, its header checks (empty, duplicate) run only when the next record is
    // asked for, after the caller has seen the previous one.
    std::string pending_line_;
    long long pending_lineno_ = 0;
    bool have_pending_ = false;
    std::unordered_map<std::string, int> seen_;
};

bool is_clustal(const std::string& path);
// PAL2NAL's '#' mark line of a CLUSTAL alignment (one character per column; see seqio.cpp).
std::string clustal_mask(const std::string& path);

class DnaIndex {
public:
    explicit DnaIndex(const std::vector<std::string>& paths);
    size_t size() const { return start_.size(); }
    const std::vector<std::string>& ids() const { return ids_; }
    const std::unordered_map<std::string, size_t>& pos() const { return pos_; }
    void forget_ids();
    // Sequence i: '#' lines skipped, every non-letter removed (as PAL2NAL).
    std::string get(size_t i);
    // Number of records in each input file.
    std::vector<long long> counts_per_file() const {
        std::vector<long long> c(paths_.size(), 0);
        for (auto f : fno_) ++c[f];
        return c;
    }
    // Length of record i before cleaning (bytes between header and next record).
    long long raw_size(size_t i) const { return end_[i] - start_[i]; }

private:
    void scan(int fno, const std::string& path);
    std::vector<std::string> paths_;
    std::vector<std::string> ids_;
    std::unordered_map<std::string, size_t> pos_;
    std::vector<uint16_t> fno_;
    std::vector<long long> start_, end_;
    std::vector<std::unique_ptr<std::ifstream>> fh_;
    std::string raw_;
};

}  // namespace coal
