#include "seqio.hpp"

#include <cstring>

namespace coal {

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

static std::string no_such_file(const std::string& path) {
    return "[Errno 2] No such file or directory: " + py_repr(path);
}

static void check_ascii(const std::string& s, long long offset) {
    for (size_t i = 0; i < s.size(); ++i)
        if ((unsigned char)s[i] > 127) {
            static const char* hex = "0123456789abcdef";
            unsigned char c = (unsigned char)s[i];
            throw InputError(std::string("'ascii' codec can't decode byte 0x") + hex[c >> 4] + hex[c & 15] +
                             " in position " + std::to_string(offset + (long long)i) +
                             ": ordinal not in range(128)");
        }
}

// ------------------------------------------------------------------ LineReader

LineReader::LineReader(const std::string& path) : path_(path), in_(path, std::ios::binary), buf_(1 << 16) {
    if (!in_) throw InputError(no_such_file(path));
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
                check_ascii(line, offset_);
                ++lineno_;
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
            check_ascii(line, offset_);
            ++lineno_;
            return true;
        }
    }
}

// ------------------------------------------------------------------ protein

bool is_clustal(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw InputError(no_such_file(path));
    char head[7] = {};
    in.read(head, 7);
    return in.gcount() == 7 && std::memcmp(head, "CLUSTAL", 7) == 0;
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
        lr.next(line);
        std::vector<std::string> order;
        std::unordered_map<std::string, std::string> seqs;
        while (lr.next(line)) {
            if (blank(line) || line[0] == ' ' || line[0] == '\t') continue;
            auto parts = split_ws(line);
            if (parts.size() < 2) throw InputError(path + ": malformed CLUSTAL line: " + py_repr(line));
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
                throw InputError(path_ + ":" + std::to_string(lr_->lineno()) + ": empty FASTA header");
            const std::string& rid = parts[0];
            if (seen_.count(rid)) throw InputError(path_ + ": duplicate sequence ID " + py_repr(rid));
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
                throw InputError(path_ + ":" + std::to_string(lr_->lineno()) +
                                 ": sequence data before first '>' header");
            for (unsigned char c : line)
                if (!py_space(c)) row += char(c);
        }
    }
    have_pending_ = false;
    return in_rec;
}

// ------------------------------------------------------------------ DNA index

DnaIndex::DnaIndex(const std::vector<std::string>& paths) : paths_(paths) {
    for (size_t f = 0; f < paths_.size(); ++f) scan((int)f, paths_[f]);
    fh_.resize(paths_.size());
}

void DnaIndex::scan(int fno, const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw InputError(no_such_file(path));
    std::vector<char> buf(1 << 20);
    long long offset = 0, start = -1;
    std::string line;  // current line (only the first bytes are kept for long lines)
    bool line_start = true, header = false, data_seen = false;
    long long line_off = 0;
    auto finish_line = [&]() {
        if (header) {
            // parse the ID from the header line
            size_t i = 1;
            while (i < line.size() && bytes_space((unsigned char)line[i])) ++i;
            size_t j = i;
            while (j < line.size() && !bytes_space((unsigned char)line[j])) ++j;
            if (j == i) throw InputError(path + ": empty FASTA header");
            std::string rid = line.substr(i, j - i);
            check_ascii(rid, line_off + (long long)i);
            if (start >= 0) {
                fno_.push_back((uint16_t)fno);
                start_.push_back(start);
                end_.push_back(line_off);
            }
            std::string joined;
            for (size_t k = 0; k < paths_.size(); ++k) joined += (k ? ", " : "") + paths_[k];
            if (pos_.count(rid)) throw InputError(joined + ": duplicate sequence ID " + py_repr(rid));
            pos_[rid] = ids_.size();
            ids_.push_back(rid);
            start = offset;  // byte after the header line
        } else if (start < 0 && data_seen) {
            throw InputError(path + ": sequence data before first '>' header");
        }
    };
    for (;;) {
        in.read(buf.data(), (std::streamsize)buf.size());
        size_t n = (size_t)in.gcount();
        if (n == 0) break;
        for (size_t k = 0; k < n; ++k) {
            char c = buf[k];
            if (line_start) {
                line_start = false;
                line.clear();
                header = c == '>';
                data_seen = false;
                line_off = offset;
            }
            ++offset;
            if (header) line.push_back(c);
            else if (!bytes_space((unsigned char)c)) data_seen = true;
            if (c == '\n') {
                finish_line();
                line_start = true;
            }
        }
    }
    if (!line_start) finish_line();
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

std::string DnaIndex::get(size_t i) {
    int f = fno_[i];
    if (!fh_[f]) fh_[f] = std::make_unique<std::ifstream>(paths_[f], std::ios::binary);
    std::ifstream& in = *fh_[f];
    long long n = end_[i] - start_[i];
    raw_.resize((size_t)n);
    in.clear();
    in.seekg(start_[i]);
    in.read(&raw_[0], (std::streamsize)n);
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
