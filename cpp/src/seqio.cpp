#include "seqio.hpp"

#include <algorithm>

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>

namespace coal {

// Plain-words reader messages, as src/coal/seqio.py.
const char* const NON_ASCII =
    "non-ASCII character: sequence files must be plain text (is this a binary or word-processor file?)";
const char* const NON_ASCII_ID =
    "non-ASCII character in the sequence ID: use only ASCII letters, digits and punctuation in IDs";
const char* const CONTROL_CHAR =
    "control character: sequence files must be plain text (is this a binary or word-processor file?)";

// NUL and the other C0 control bytes, except tab, LF, VT, FF and CR (seqio.py _CONTROL)
static inline bool control_byte(unsigned char c) { return c < 0x20 && !(c >= 0x09 && c <= 0x0d); }
const char* const EMPTY_HEADER =
    "empty FASTA header: a '>' line without a sequence ID (every record needs a name after '>')";
#define COAL_BEFORE "sequence data before the first '>' header: not a FASTA file? (each FASTA record starts with a '>name' line"
const char* const BEFORE_HEADER_PROTEIN = COAL_BEFORE "; a CLUSTAL file must start with 'CLUSTAL')";
const char* const BEFORE_HEADER_DNA = COAL_BEFORE ")";
#undef COAL_BEFORE

std::string short_id(const std::string& rid) {
    return rid.size() <= 200 ? rid : rid.substr(0, 100) + "..." + rid.substr(rid.size() - 20);
}

std::string duplicate_message(const std::string& where, const std::string& rid) {
    // (raised before the IDs are known to the message layer: shortened here)
    return where + ": duplicate sequence ID " + py_repr(short_id(rid)) + ": every sequence needs its own ID";
}

std::string open_error(const std::string& path) {
    std::error_code ec;
    std::string why;
    if (std::filesystem::is_directory(path, ec)) {
        why = "it is a folder, not a file";
    } else {
        errno = 0;
        std::FILE* f = std::fopen(path.c_str(), "rb");
        int e = errno;
        if (f) std::fclose(f);
        if (!f && e == ENOENT) why = "no such file";
        else if (!f && (e == EACCES || e == EPERM)) why = "permission denied";
        else why = std::strerror(e ? e : EIO);
    }
    return "cannot open '" + path + "': " + why;
}

bool py_space(unsigned char c) {
    return c == ' ' || (c >= 9 && c <= 13) || (c >= 0x1c && c <= 0x1f);
}

static bool bytes_space(unsigned char c) { return c == ' ' || (c >= 9 && c <= 13); }

std::string py_repr(const std::string& s) {
    bool sq = s.find('\'') != std::string::npos, dq = s.find('"') != std::string::npos;
    char q = (sq && !dq) ? '"' : '\'';
    std::string out(1, q);
    for (unsigned char c : s) {
        if (c == '\\') out += "\\\\";
        else if (c == (unsigned char)q) out += std::string("\\") + char(q);
        else if (c == '\t') out += "\\t";
        else if (c == '\n') out += "\\n";
        else if (c == '\r') out += "\\r";
        else if (c < 32 || c == 127) {
            static const char* hex = "0123456789abcdef";
            out += std::string("\\x") + hex[c >> 4] + hex[c & 15];
        } else out += char(c);
    }
    out += q;
    return out;
}

// Length of the Python str.split() whitespace character (ASCII or UTF-8 encoded) at s[i], or 0.
static size_t uni_space(const std::string& s, size_t i) {
    unsigned char c = (unsigned char)s[i];
    if (c < 128) return py_space(c) ? 1 : 0;
    auto at = [&](size_t k) -> unsigned char { return i + k < s.size() ? (unsigned char)s[i + k] : 0; };
    if (c == 0xC2 && (at(1) == 0x85 || at(1) == 0xA0)) return 2;                    // U+0085, U+00A0
    if (c == 0xE1 && at(1) == 0x9A && at(2) == 0x80) return 3;                      // U+1680
    if (c == 0xE2 && at(1) == 0x80 && (at(2) <= 0x8A || at(2) == 0xA8 || at(2) == 0xA9 || at(2) == 0xAF) &&
        at(2) >= 0x80)
        return 3;                                                                    // U+2000-200A, 2028, 2029, 202F
    if (c == 0xE2 && at(1) == 0x81 && at(2) == 0x9F) return 3;                      // U+205F
    if (c == 0xE3 && at(1) == 0x80 && at(2) == 0x80) return 3;                      // U+3000
    return 0;
}

// Text inputs must be ASCII (a UTF-8 BOM is skipped): "<file>:<line>: non-ASCII character", as the reference
// engine. A header's description (after the ID) may hold any text: the line becomes ">ID"; a non-ASCII ID is
// "<file>:<line>: non-ASCII character in the sequence ID".
static void check_ascii(std::string& s, const std::string& path, long long lineno) {
    bool ascii = true;
    for (unsigned char c : s) {
        if (control_byte(c)) throw InputError(path + ":" + std::to_string(lineno) + ": " + CONTROL_CHAR);
        if (c > 127) ascii = false;
    }
    if (ascii) return;
    if (!s.empty() && s[0] == '>') {
        size_t i = 1, n = s.size(), w;
        while (i < n && (w = uni_space(s, i)) > 0) i += w;
        size_t j = i;
        while (j < n && uni_space(s, j) == 0) ++j;
        bool id_ascii = j > i;
        for (size_t k = i; k < j && id_ascii; ++k) id_ascii = (unsigned char)s[k] < 128;
        if (id_ascii) {
            s = ">" + s.substr(i, j - i);
            return;
        }
        throw InputError(path + ":" + std::to_string(lineno) + ": " + NON_ASCII_ID);
    }
    throw InputError(path + ":" + std::to_string(lineno) + ": " + NON_ASCII);
}

// ------------------------------------------------------------------ LineReader

LineReader::LineReader(const std::string& path) : path_(path), in_(path, std::ios::binary), buf_(1 << 16) {
    if (!in_) throw InputError(open_error(path), "file_error");
    if (fill() && len_ >= 3 && (unsigned char)buf_[0] == 0xEF && (unsigned char)buf_[1] == 0xBB &&
        (unsigned char)buf_[2] == 0xBF)
        pos_ = 3;  // UTF-8 byte-order mark
}

bool LineReader::fill() {
    if (eof_) return false;
    in_.read(buf_.data(), (std::streamsize)buf_.size());
    len_ = (size_t)in_.gcount();
    pos_ = 0;
    if (len_ == 0) eof_ = true;
    return len_ > 0;
}

bool LineReader::next(std::string& line) {
    line.clear();
    bool got = false;
    for (;;) {
        if (pos_ >= len_ && !fill()) {
            if (got) {
                ++lineno_;
                check_ascii(line, path_, lineno_);
            }
            return got;
        }
        got = true;
        char* b = buf_.data() + pos_;
        size_t n = len_ - pos_;
        size_t k = 0;
        while (k < n && b[k] != '\n' && b[k] != '\r') ++k;
        line.append(b, k);
        pos_ += k;
        if (k < n) {
            char t = b[k];
            ++pos_;
            if (t == '\r') {  // \r\n counts as one terminator
                if (pos_ >= len_) fill();
                if (pos_ < len_ && buf_[pos_] == '\n') ++pos_;
            }
            ++lineno_;
            check_ascii(line, path_, lineno_);
            return true;
        }
    }
}

// ------------------------------------------------------------------ protein

static bool blank(const std::string& s);

// The first non-blank line starts with CLUSTAL.
bool is_clustal(const std::string& path) {
    LineReader lr(path);
    std::string line;
    while (lr.next(line))
        if (!blank(line)) return line.compare(0, 7, "CLUSTAL") == 0;
    return false;
}

static std::vector<std::string> split_ws(const std::string& s, size_t from = 0) {
    std::vector<std::string> parts;
    size_t i = from;
    while (i < s.size()) {
        while (i < s.size() && py_space((unsigned char)s[i])) ++i;
        size_t j = i;
        while (j < s.size() && !py_space((unsigned char)s[j])) ++j;
        if (j > i) parts.push_back(s.substr(i, j - i));
        i = j;
    }
    return parts;
}

static bool blank(const std::string& s) {
    for (unsigned char c : s)
        if (!py_space(c)) return false;
    return true;
}

ProteinReader::ProteinReader(const std::string& path) : path_(path) {
    clustal_ = is_clustal(path);
    if (clustal_) {
        LineReader lr(path);
        std::string line;
        while (lr.next(line) && blank(line)) {
        }
        std::vector<std::string> order;
        std::unordered_map<std::string, std::string> seqs;
        while (lr.next(line)) {
            // '#' mark lines and repeated CLUSTAL header lines are skipped too
            if (blank(line) || line[0] == ' ' || line[0] == '\t' || line[0] == '#' || line.compare(0, 7, "CLUSTAL") == 0)
                continue;
            auto parts = split_ws(line);
            if (parts.size() < 2)
                throw InputError(path + ":" + std::to_string(lr.lineno()) + ": malformed CLUSTAL line " + py_repr(line) +
                                 ": each line of a block must hold a sequence name, then (after spaces) its aligned "
                                 "sequence");
            auto it = seqs.find(parts[0]);
            if (it == seqs.end()) {
                order.push_back(parts[0]);
                seqs[parts[0]] = parts[1];
            } else {
                it->second += parts[1];
            }
        }
        for (auto& n : order) clustal_recs_.push_back({n, seqs[n]});
    } else {
        lr_ = std::make_unique<LineReader>(path);
    }
}

bool ProteinReader::next(std::string& id, std::string& row) {
    if (clustal_) {
        if (ci_ >= clustal_recs_.size()) return false;
        id = clustal_recs_[ci_].first;
        row = clustal_recs_[ci_].second;
        ++ci_;
        return true;
    }
    std::string line;
    row.clear();
    bool in_rec = have_pending_;
    if (in_rec) id = pending_id_;
    while (lr_->next(line)) {
        if (!line.empty() && line[0] == '>') {
            auto parts = split_ws(line, 1);
            if (parts.empty())
                throw InputError(path_ + ":" + std::to_string(lr_->lineno()) + ": " + EMPTY_HEADER);
            const std::string& rid = parts[0];
            if (seen_.count(rid)) throw InputError(duplicate_message(path_, rid));
            seen_[rid] = 1;
            if (in_rec) {
                pending_id_ = rid;
                have_pending_ = true;
                return true;
            }
            id = rid;
            in_rec = true;
            have_pending_ = false;
        } else if (!blank(line)) {
            if (!in_rec)
                throw InputError(path_ + ":" + std::to_string(lr_->lineno()) + ": " + BEFORE_HEADER_PROTEIN);
            bool plain = true;
            for (unsigned char c : line)
                if (py_space(c)) {
                    plain = false;
                    break;
                }
            if (plain) row += line;
            else
                for (unsigned char c : line)
                    if (!py_space(c)) row += char(c);
        }
    }
    have_pending_ = false;
    return in_rec;
}

// As clustal_mask in the reference: PAL2NAL's '#' mark line, one character per alignment column. After
// each sequence line, the next line that is not a sequence line is padded with spaces to that sequence line's
// length and the part under its sequence is appended. A block whose last sequence line ends the file adds nothing.
std::string clustal_mask(const std::string& path) {
    std::vector<std::string> lines;
    {
        LineReader lr(path);
        std::string line;
        while (lr.next(line)) lines.push_back(line);
    }
    while (!lines.empty() && lines.back().empty()) lines.pop_back();
    std::string mask;
    bool get = false;
    size_t tmplen = 0, idspc = 0, sublen = 0;
    for (std::string line : lines) {
        unsigned char c0 = line.empty() ? 0 : (unsigned char)line[0];
        bool seq_line = line.compare(0, 7, "CLUSTAL") != 0 && !line.empty() && c0 != '#' && !py_space(c0);
        if (seq_line) {
            size_t e = line.size();
            while (e > 0 && py_space((unsigned char)line[e - 1])) --e;
            line.resize(e);
            tmplen = line.size();
            size_t i = 0;
            while (i < line.size() && !py_space((unsigned char)line[i])) ++i;
            if (i > 0 && i < line.size()) {  // re.match(r"\S+\s+")
                while (i < line.size() && py_space((unsigned char)line[i])) ++i;
                idspc = i;
            }
            auto parts = split_ws(line);
            sublen = parts.size() > 1 ? parts[1].size() : 0;
            get = true;
        } else if (get) {
            if (line.size() < tmplen) line.append(tmplen - line.size(), ' ');
            if (idspc < line.size()) mask += line.substr(idspc, sublen);
            get = false;
        }
    }
    return mask;
}

// ------------------------------------------------------------------ DNA index

DnaIndex::DnaIndex(const std::vector<std::string>& paths) : paths_(paths) {
    for (size_t f = 0; f < paths_.size(); ++f) scan((int)f, paths_[f]);
    fh_.resize(paths_.size());
}

void DnaIndex::scan(int fno, const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw InputError(open_error(path), "file_error");
    long long offset = 0;
    {  // UTF-8 byte-order mark
        char bom[3] = {};
        in.read(bom, 3);
        if (in.gcount() == 3 && (unsigned char)bom[0] == 0xEF && (unsigned char)bom[1] == 0xBB &&
            (unsigned char)bom[2] == 0xBF)
            offset = 3;
        in.clear();
        in.seekg(offset);
    }
    // Classic Mac files: a lone \r also ends a line (decided from the first 64 kB, as the Python engine)
    bool cr_mode = false;
    {
        std::string c(1 << 16, '\0');
        in.read(&c[0], (std::streamsize)c.size());
        c.resize((size_t)in.gcount());
        if (!c.empty() && c.back() == '\r') {
            char x;
            if (in.get(x)) c.push_back(x);
        }
        size_t r = 0, rn = 0;
        for (size_t i = 0; i < c.size(); ++i)
            if (c[i] == '\r') {
                ++r;
                if (i + 1 < c.size() && c[i + 1] == '\n') ++rn;
            }
        cr_mode = r != rn + (!c.empty() && c.back() == '\r' ? 1 : 0);
        in.clear();
        in.seekg(offset);
    }
    std::string joined;
    for (size_t k = 0; k < paths_.size(); ++k) joined += (k ? ", " : "") + paths_[k];
    std::vector<char> buf(1 << 20);
    long long start = -1, lineno = 0, line_off = offset;
    std::string line;  // header lines only
    bool line_start = true, header = false, hash = false, data_seen = false, pending_cr = false;
    auto finish_line = [&]() {
        ++lineno;
        if (header) {
            size_t i = 1;
            while (i < line.size() && bytes_space((unsigned char)line[i])) ++i;
            size_t j = i;
            while (j < line.size() && !bytes_space((unsigned char)line[j])) ++j;
            if (start >= 0) {
                fno_.push_back((uint16_t)fno);
                start_.push_back(start);
                end_.push_back(line_off);
            }
            for (unsigned char ch : line)  // (the sequence data are checked when they are read: get())
                if (control_byte(ch)) throw InputError(path + ":" + std::to_string(lineno) + ": " + CONTROL_CHAR);
            if (j == i) throw InputError(path + ":" + std::to_string(lineno) + ": " + EMPTY_HEADER);
            std::string rid = line.substr(i, j - i);
            for (unsigned char ch : rid)
                if (ch > 127) throw InputError(path + ":" + std::to_string(lineno) + ": " + NON_ASCII_ID);
            if (pos_.count(rid)) throw InputError(duplicate_message(joined, rid));
            pos_[rid] = ids_.size();
            ids_.push_back(rid);
            start = offset;  // byte after the header line
        } else if (start < 0 && data_seen && !hash) {
            throw InputError(path + ":" + std::to_string(lineno) + ": " + BEFORE_HEADER_DNA);
        }
        line_start = true;
    };
    for (;;) {
        in.read(buf.data(), (std::streamsize)buf.size());
        size_t n = (size_t)in.gcount();
        if (n == 0) break;
        for (size_t k = 0; k < n; ++k) {
            char c = buf[k];
            if (pending_cr) {
                pending_cr = false;
                if (c == '\n') {  // \r\n: one terminator
                    ++offset;
                    finish_line();
                    continue;
                }
                finish_line();
            }
            if (!line_start && !header && start >= 0) {
                // inside sequence data: only the end of the line matters
                size_t q = k;
                if (cr_mode) {
                    while (q < n && buf[q] != '\n' && buf[q] != '\r') ++q;
                } else {
                    const char* nl = (const char*)std::memchr(buf.data() + k, '\n', n - k);
                    q = nl ? (size_t)(nl - buf.data()) : n;
                }
                offset += (long long)(q - k);
                k = q;
                if (k >= n) break;
                c = buf[k];
            }
            if (line_start) {
                line_start = false;
                line.clear();
                header = c == '>';
                hash = c == '#';
                data_seen = false;
                line_off = offset;
            }
            ++offset;
            if (header) line.push_back(c);
            else if (!bytes_space((unsigned char)c)) data_seen = true;
            if (c == '\n') finish_line();
            else if (cr_mode && c == '\r') pending_cr = true;
        }
    }
    if (pending_cr || !line_start) finish_line();
    if (start >= 0) {
        fno_.push_back((uint16_t)fno);
        start_.push_back(start);
        end_.push_back(offset);
    }
}

void DnaIndex::forget_ids() {
    std::vector<std::string>().swap(ids_);
    std::unordered_map<std::string, size_t>().swap(pos_);
}

// The 1-based line number of a byte offset (for an error message only: reads the file up to there).
long long DnaIndex::line_of(int f, long long offset) const {
    std::ifstream in(paths_[f], std::ios::binary);
    std::vector<char> buf(1 << 20);
    long long nl = 0, cr = 0, left = offset;
    while (left > 0 && in) {
        in.read(buf.data(), (std::streamsize)std::min<long long>((long long)buf.size(), left));
        std::streamsize got = in.gcount();
        if (got <= 0) break;
        for (std::streamsize k = 0; k < got; ++k) {
            if (buf[(size_t)k] == '\n') ++nl;
            else if (buf[(size_t)k] == '\r') ++cr;
        }
        left -= got;
    }
    return (nl ? nl : cr) + 1;  // old Mac line ends when there is no \n
}

std::string DnaIndex::get(size_t i) {
    int f = fno_[i];
    if (!fh_[f]) fh_[f] = std::make_unique<std::ifstream>(paths_[f], std::ios::binary);
    std::ifstream& in = *fh_[f];
    // Never fall back to the previous record's bytes: a file that cannot be (re)opened, or that
    // shrank or changed after it was indexed, is an error, not another gene's DNA.
    if (!in.is_open()) throw InputError(open_error(paths_[f]), "file_error");
    long long n = end_[i] - start_[i];
    raw_.resize((size_t)n);
    in.clear();
    in.seekg(start_[i]);
    if (n > 0) in.read(&raw_[0], (std::streamsize)n);  // (gcount() is only meaningful after a read)
    if (!in || (n > 0 && (long long)in.gcount() != n)) {
        throw InputError(paths_[f] + ": could not read the DNA of record " + std::to_string(i + 1) +
                             " (the file changed after it was read, or cannot be read)",
                         "file_error");
    }
    for (long long k = 0; k < n; ++k)
        if (control_byte((unsigned char)raw_[(size_t)k]))
            throw InputError(paths_[f] + ":" + std::to_string(line_of(f, start_[i] + k)) + ": " + CONTROL_CHAR);
    std::string out;
    out.reserve((size_t)n);
    bool at_line_start = true, skip = false;
    for (long long k = 0; k < n; ++k) {
        char c = raw_[(size_t)k];
        if (at_line_start) {
            skip = c == '#';
            at_line_start = false;
        }
        if (c == '\n' || c == '\r') {
            at_line_start = true;
            continue;
        }
        if (!skip && ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'))) out.push_back(c);
    }
    return out;
}

}  // namespace coal
