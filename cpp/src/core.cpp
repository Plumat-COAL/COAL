#include "core.hpp"

#include <algorithm>
#include <set>

namespace coal {
namespace {

enum Mode { NORMAL = 0, START = 1, ANCHOR_START = 2 };
constexpr int ANCHOR = 10;
constexpr int MIN_CHECKABLE = 10;            // reference R3-E1
constexpr int AMBIGUOUS_MAX_RESIDUES = 3;    // reference FZ-6

inline bool is_digit(char c) { return c >= '0' && c <= '9'; }
inline bool as_x(char c) { return c == 'B' || c == 'Z' || c == 'J' || c == 'O' || c == 'U'; }
inline bool valid_char(char c) {
    if (c == '-' || c == '*' || is_digit(c)) return true;
    if (c < 'A' || c > 'Z') return false;
    return true;  // every upper-case letter: 20 standard + X + BZJOU
}

// One element of the pattern searched in the DNA: a residue (3 nt that must
// encode it, unless `any`) or a run of n nucleotides (frameshift digit or wildcard).
struct Token {
    bool residue;
    char aa;          // residue letter ('X' for B/Z/J/O/U)
    int mode;
    int n;            // length in nt
    int col;          // 1-based alignment column
    bool any;         // 'X': anything matches
    uint64_t allowed; // codons that encode it
};

uint64_t allowed_for(const GeneticCode& code, char aa, int mode) {
    if (aa == 'M' && mode == START) return code.starts;
    // CPP-22 (reference FIN-E1): the fallback's first Met also accepts the selected code's start and Met codons
    // (was ATG/CTG/GTG only, PAL2NAL's rule; deliberate difference from PAL2NAL v14)
    if (aa == 'M' && mode == ANCHOR_START) return anchor_start_mask() | code.starts | code.codons[(unsigned char)'M'];
    return code.codons[(unsigned char)aa];
}

// PAL2NAL's "first Met": an M before any non-X residue gets start_mode.
std::vector<Token> tokens(const std::string& row, size_t from, size_t to, int start_mode, int col0,
                          const GeneticCode& code) {
    std::vector<Token> toks;
    toks.reserve(to - from);
    bool seen_specific = false;
    for (size_t i = from; i < to; ++i) {
        char ch = row[i];
        int col = col0 + int(i - from);
        if (ch == '-') continue;
        if (is_digit(ch)) {
            toks.push_back({false, 0, NORMAL, ch - '0', col, false, 0});
            continue;
        }
        char aa = as_x(ch) ? 'X' : ch;
        int mode = (aa == 'M' && !seen_specific) ? start_mode : NORMAL;
        if (aa != 'X') seen_specific = true;
        toks.push_back({true, aa, mode, 3, col, aa == 'X', allowed_for(code, aa, mode)});
    }
    return toks;
}

long long total_length(const std::vector<Token>& toks) {
    long long n = 0;
    for (const Token& t : toks) n += t.n;
    return n;
}

// DNA upper-cased with U->T, as IUPAC codes (255 = not an IUPAC letter).
struct Target {
    std::vector<uint8_t> code;
    const uint64_t* triplet = iupac_triplet_table();

    // Does the triplet at pos encode a residue with these allowed codons?
    inline bool ok(long long pos, uint64_t allowed) const {
        const size_t i = (size_t)pos;  // CPP-7: explicit (size_t is 32-bit on 32-bit Windows: MSVC C4244)
        uint8_t a = code[i], b = code[i + 1], c = code[i + 2];
        if ((a | b | c) & 0xF0) return false;  // 255 = not IUPAC (valid codes are 0..14)
        return (triplet[a * 225 + b * 15 + c] & ~allowed) == 0;
    }
};

inline bool residue_ok(const Target& t, long long pos, const Token& tok) { return t.ok(pos, tok.allowed); }

bool matches_at(const Target& t, long long s, const std::vector<Token>& toks) {
    long long pos = s;
    for (const Token& tok : toks) {
        if (tok.residue && !tok.any && !residue_ok(t, pos, tok)) return false;
        pos += tok.n;
    }
    return pos <= (long long)t.code.size();
}

// Leftmost start where every token matches consecutively, or -1.
long long find(const Target& t, const std::vector<Token>& toks) {
    long long last = (long long)t.code.size() - total_length(toks);
    if (last < 0) return -1;
    // Filter on the first checkable residue (fixed offset from the start).
    long long off = 0;
    const Token* first = nullptr;
    for (const Token& tok : toks) {
        if (tok.residue && !tok.any) {
            first = &tok;
            break;
        }
        off += tok.n;
    }
    for (long long s = 0; s <= last; ++s) {
        if (first && !residue_ok(t, s + off, *first)) continue;
        if (matches_at(t, s, toks)) return s;
    }
    return -1;
}

// Chunks of 10 non-gap characters; a final chunk shorter than 10 characters
// is merged into the previous one (PAL2NAL rule). Returns [from, to) ranges.
std::vector<std::pair<size_t, size_t>> anchors(const std::string& row) {
    std::vector<std::pair<size_t, size_t>> pre;
    size_t start = 0;
    int n = 0;
    for (size_t i = 0; i < row.size(); ++i) {
        if (row[i] != '-') ++n;
        if (n == ANCHOR || i == row.size() - 1) {
            pre.push_back({start, i + 1});
            start = i + 1;
            n = 0;
        }
    }
    if (!pre.empty() && pre.back().second - pre.back().first < (size_t)ANCHOR) {
        if (pre.size() >= 2) {
            size_t end = pre.back().second;
            pre.pop_back();
            pre.back().second = end;
        } else {
            pre.clear();
        }
    }
    return pre;
}

// Python repr() of a one-character string, for the invalid-character message.
std::string py_repr_char(unsigned char c) {
    if (c == '\'') return "\"'\"";
    if (c == '\\') return "'\\\\'";
    if (c == '\t') return "'\\t'";
    if (c == '\n') return "'\\n'";
    if (c == '\r') return "'\\r'";
    if (c < 32 || c == 127) {
        static const char* hex = "0123456789abcdef";
        return std::string("'\\x") + hex[c >> 4] + hex[c & 15] + "'";
    }
    return std::string("'") + char(c) + "'";
}

}  // namespace

std::string prepare_row(const std::string& aln) {
    std::string r = aln;
    bool tfastx = false;
    for (char& c : r) {
        if (c >= 'a' && c <= 'z') c = char(c - 32);
        else if (c == '\\' || c == '/') tfastx = true;
    }
    if (tfastx) {
        // CPP-14 (reference FZ-3): tfastx/tfasty frameshift notation, read as PAL2NAL v14 reads it: '\\' = 1 base
        // (digit 1); '/' + gaps + residue = a gap, the gaps, then 2 bases (digit 2). Any other '/' stays (invalid).
        std::string t;
        t.reserve(r.size());
        for (size_t i = 0; i < r.size();) {
            char c = r[i];
            if (c == '\\') {
                t += '1';
                ++i;
            } else if (c == '/') {
                size_t j = i + 1;
                while (j < r.size() && r[j] == '-') ++j;
                if (j < r.size() && ((r[j] >= 'A' && r[j] <= 'Z') || r[j] == '*')) {
                    t += '-';
                    t.append(r, i + 1, j - i - 1);
                    t += '2';
                    i = j + 1;
                } else {
                    t += '/';
                    ++i;
                }
            } else {
                t += c;
                ++i;
            }
        }
        r.swap(t);
    }
    for (char& c : r)
        if (c == '.') c = '-';
    return r;
}

// CPP-15 (reference FZ-6): a protein of at most 3 residues that fits in several places of the DNA gets a notice;
// the leftmost is used, as PAL2NAL. stdout is unchanged.
static void note_ambiguous(Threaded& res, const std::string& seq_id, const Target& target,
                           const std::vector<Token>& toks) {
    long long nres = 0;
    for (const Token& t : toks) nres += t.residue;
    if (nres > AMBIGUOUS_MAX_RESIDUES) return;
    long long last = (long long)target.code.size() - total_length(toks), n = 0;
    for (long long st = 0; st <= last; ++st) n += matches_at(target, st, toks);
    if (n > 1)
        res.notices.push_back(seq_id + ": placement is ambiguous: " + std::to_string(n) +
                              " possible positions in the DNA for this " + std::to_string(nres) +
                              "-residue protein; used the leftmost (as PAL2NAL)");
}

Threaded thread_one(const std::string& seq_id, const std::string& row, const std::string& raw_dna,
                    const GeneticCode& code, bool strict) {
    std::set<unsigned char> bad;
    for (char c : row)
        if (!valid_char(c)) bad.insert((unsigned char)c);
    if (!bad.empty()) {
        std::string list = "[";
        bool firstc = true;
        for (unsigned char c : bad) {
            if (!firstc) list += ", ";
            list += py_repr_char(c);
            firstc = false;
        }
        throw CoalError(seq_id + ": invalid protein alignment characters " + list + "]", "invalid_chars", seq_id);
    }

    Threaded res;
    for (size_t i = 0; i < row.size(); ++i)
        if (as_x(row[i])) {
            // CPP-23 (reference FIN-C14): was "U unknown AA type. Taken as 'X'" (PAL2NAL's wording)
            const char* name = row[i] == 'U' ? "selenocysteine" : row[i] == 'O' ? "pyrrolysine"
                               : row[i] == 'B' ? "D or N" : row[i] == 'Z' ? "E or Q" : "I or L";
            res.warnings.push_back("pepAlnPos " + std::to_string(i + 1) + ": " + row[i] + " (" + name +
                                   ") is not checked; any codon is accepted");
            res.mismatch_cols.push_back(int(i + 1));  // CPP-11 (DEC-4f): -nomismatch removes it, as PAL2NAL
        }

    bool any_residue = false;
    for (char c : row)
        if (c != '-') {
            any_residue = true;
            break;
        }
    if (!any_residue) throw CoalError(seq_id + ": protein sequence is empty", "empty_protein", seq_id);

    Target target;
    target.code.resize(raw_dna.size());
    {
        // upper-case, U -> T, then IUPAC code: one table lookup per base
        uint8_t map[256];
        const uint8_t* iu = iupac_table();
        for (int c = 0; c < 256; ++c) {
            int u = (c >= 'a' && c <= 'z') ? c - 32 : c;
            if (u == 'U') u = 'T';
            map[c] = iu[u];
        }
        const unsigned char* src = (const unsigned char*)raw_dna.data();
        uint8_t* dst = target.code.data();
        for (size_t i = 0; i < raw_dna.size(); ++i) dst[i] = map[src[i]];
    }

    // 1. exact fit
    std::vector<Token> toks = tokens(row, 0, row.size(), START, 1, code);
    long long len = total_length(toks);
    long long s = find(target, toks);
    if (s >= 0) {
        res.dna = raw_dna.substr((size_t)s, (size_t)len);
        note_ambiguous(res, seq_id, target, toks);
        return res;
    }

    // 2. anchor fallback
    auto ranges = anchors(row);
    if (ranges.empty()) throw CoalError(seq_id + ": protein and DNA are inconsistent (no match found)", "no_match", seq_id);
    std::vector<Token> pattern;
    int col = 1;
    bool found_any = false;
    for (size_t k = 0; k < ranges.size(); ++k) {
        auto [a, b] = ranges[k];
        std::vector<Token> at = tokens(row, a, b, k == 0 ? ANCHOR_START : NORMAL, col, code);
        bool found = find(target, at) >= 0;
        found_any = found_any || found;
        for (Token& t : at) {
            if (!found) {  // wildcard of the same length
                t.residue = false;
                t.any = false;
            }
            pattern.push_back(t);
        }
        col += int(b - a);
    }
    s = find(target, pattern);
    if (s < 0) throw CoalError(seq_id + ": protein and DNA are inconsistent (no match found)", "no_match", seq_id);
    res.dna = raw_dna.substr((size_t)s, (size_t)len);

    // mismatch report: first residue Met -> ATG/CTG/GTG, others -> genetic code
    long long pos = 0, nres = 0, nchk = 0, nbad_chk = 0;  // CPP-15: checkable = A/C/G/T(U) bases only
    bool first = true;
    for (const Token& t : toks) {
        if (!t.residue) {
            pos += t.n;
            continue;
        }
        int mode = (first && t.aa == 'M') ? ANCHOR_START : NORMAL;
        ++nres;
        bool checkable = true;
        for (long long q = 0; q < 3; ++q) {
            char ch = (size_t)(pos + q) < res.dna.size() ? res.dna[(size_t)(pos + q)] : 'A';
            if (ch >= 'a' && ch <= 'z') ch = char(ch - 32);
            if (ch != 'A' && ch != 'C' && ch != 'G' && ch != 'T' && ch != 'U') checkable = false;
        }
        nchk += checkable;
        bool ok;
        if (t.any) {
            ok = true;
        } else {
            ok = target.ok(s + pos, allowed_for(code, t.aa, mode));
        }
        if (!ok) {
            res.warnings.push_back("pepAlnPos " + std::to_string(t.col) + ": " + row[t.col - 1] +
                                   " does not correspond to " + res.dna.substr((size_t)pos, 3));
            res.mismatch_cols.push_back(t.col);
            nbad_chk += checkable;
        }
        first = false;
        pos += 3;
    }
    // CPP-15 (reference DEC-1 with ENG-c, R2-E3, R3-E1, FZ-4): probably-wrong data. By default the row is written as
    // PAL2NAL writes it, with a notice; --strict makes the first two rules errors.
    const long long namb = nres - nchk;
    const std::string N = std::to_string(nres), K = std::to_string(nchk), B = std::to_string(nbad_chk);
    if (nchk && 2 * nbad_chk > nchk) {
        if (strict)
            throw CoalError(seq_id + ": protein and DNA are inconsistent (" + B + " of " + K +
                                " checkable codons do not match the protein)",
                            "unrelated", seq_id);
        res.notices.push_back(seq_id + ": protein and DNA are probably unrelated: " + B + " of " + K +
                              " checkable codons do not match the protein. The row was written as PAL2NAL writes "
                              "it, but it is probably the wrong DNA or the wrong genetic code (--strict makes this an "
                              "error)");
        return res;
    }
    if (!found_any && (long long)raw_dna.size() > len && nchk < MIN_CHECKABLE) {
        const std::string few = namb ? "only " + K + " of " + N + " codons have no ambiguous bases"
                                     : "the gene has only " + N + " codons";
        const std::string what = namb ? "the gene" : "it";
        if (strict)
            throw CoalError(seq_id + ": protein and DNA are inconsistent (no 10-residue segment of the protein was "
                                "found in the DNA and " + few + ": too few to tell where " + what +
                                " starts in the longer DNA)",
                            "position_guess", seq_id);
        res.notices.push_back(seq_id + ": position of the gene in the DNA is a guess: no 10-residue segment of the "
                              "protein was found in the DNA and " + few + ", too few to tell where " + what +
                              " starts in the longer DNA. The codons were taken from the start of the DNA as PAL2NAL "
                              "does, so the row may be shifted (--strict makes this an error)");
        return res;
    }
    if (!found_any)
        res.notices.push_back(seq_id + ": no 10-residue segment of the protein was found in the DNA; codons were "
                              "assigned by position (" + B + " of " + K + " checkable codons do not match" +
                              (namb ? "; " + std::to_string(namb) + " contain N or other ambiguous bases" : "") + ")");
    return res;
}

void update_column_widths(std::vector<int>& widths, const std::string& row) {
    for (size_t i = 0; i < row.size(); ++i)
        if (is_digit(row[i])) {
            int d = row[i] - '0';
            int w = (std::max(d - 1, 0) / 3 + 1) * 3;
            if (w > widths[i]) widths[i] = w;
        }
}

std::string expand(const std::string& row, const std::string& dna, const std::vector<int>& widths) {
    std::string out;
    long long total = 0;
    for (int w : widths) total += w;
    out.reserve((size_t)total);
    size_t pos = 0;
    bool plain = std::all_of(widths.begin(), widths.end(), [](int w) { return w == 3; }) &&
                 std::none_of(row.begin(), row.end(), is_digit);
    if (plain) {
        // every column 3 nt: copy runs of residues, '---' per gap
        size_t i = 0, n = row.size();
        while (i < n) {
            if (row[i] == '-') {
                out.append(3, '-');
                ++i;
                continue;
            }
            size_t j = i;
            while (j < n && row[j] != '-') ++j;
            out.append(dna, pos, 3 * (j - i));
            pos += 3 * (j - i);
            i = j;
        }
        return out;
    }
    for (size_t i = 0; i < row.size(); ++i) {
        int w = widths[i];
        char ch = row[i];
        if (ch == '-') {
            out.append((size_t)w, '-');
            continue;
        }
        int d = is_digit(ch) ? ch - '0' : 3;
        size_t take = std::min((size_t)d, dna.size() > pos ? dna.size() - pos : 0);
        out.append(dna, pos, take);
        out.append((size_t)(w - (int)take), '-');
        pos += d;
    }
    return out;
}

void mark_nogap_triplets(const std::string& e, std::vector<char>& bad, const GeneticCode& code) {
    size_t n = e.size() / 3;
    if (bad.size() < n) bad.resize(n, 0);
    const uint8_t* iu = iupac_table();
    const uint64_t* trip = iupac_triplet_table();
    const uint64_t stops = code.codons[(unsigned char)'*'];
    auto norm = [&](char ch) -> uint8_t {
        int u = (unsigned char)ch;
        if (u >= 'a' && u <= 'z') u -= 32;
        if (u == 'U') u = 'T';
        return iu[u];
    };
    for (size_t k = 0; k < n; ++k) {
        const char* t = e.data() + 3 * k;
        if (t[0] == '-' || t[1] == '-' || t[2] == '-') {
            bad[k] = 1;
            continue;
        }
        uint8_t a = norm(t[0]), b = norm(t[1]), c = norm(t[2]);
        if ((a | b | c) & 0xF0) continue;  // not an IUPAC triplet
        if (stops && (trip[a * 225 + b * 15 + c] & ~stops) == 0) bad[k] = 1;
    }
}

std::vector<std::pair<long long, long long>> keep_slices(const std::vector<int>& widths,
                                                         const std::vector<char>& drop_cols,
                                                         const std::vector<char>& bad_triplets) {
    std::vector<std::pair<long long, long long>> sl;
    long long pos = 0;
    for (size_t i = 0; i < widths.size(); ++i) {
        int w = widths[i];
        if (!(i < drop_cols.size() && drop_cols[i])) {
            for (long long k = pos; k < pos + w; k += 3) {
                size_t t = (size_t)(k / 3);
                if (t < bad_triplets.size() && bad_triplets[t]) continue;
                if (!sl.empty() && sl.back().second == k) sl.back().second = k + 3;
                else sl.push_back({k, k + 3});
            }
        }
        pos += w;
    }
    return sl;
}

}  // namespace coal
