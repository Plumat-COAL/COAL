// Mirrors src/coal/report.py. See report.hpp.
#include "report.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <cctype>
#include <set>

#include "core.hpp"
#include "tables.hpp"

namespace coal {

static const char* VERSION = "0.5.0rc3";
static const char* SCHEMA = "coal.report/1";
static const char* REPO_URL = "https://github.com/plumat-coal/COAL";
// The hint re-tests at most 50 sequences and drops a code as soon as it cannot fit as many as the best so far;
// break locations are found for up to 200 failed genes (report.py HINT_SAMPLE, BREAK_SAMPLE)
static const int HINT_SAMPLE = 50;
static const int BREAK_SAMPLE = 200;

static std::string code_name(int t) {  // CODE_NAMES (names inside a sentence)
    switch (t) {
        case 1: return "standard";
        case 2: return "vertebrate mitochondrial";
        case 3: return "yeast mitochondrial";
        case 4: return "mold, protozoan and coelenterate mitochondrial and Mycoplasma/Spiroplasma";
        case 5: return "invertebrate mitochondrial";
        case 6: return "ciliate, dasycladacean and Hexamita nuclear";
        case 9: return "echinoderm and flatworm mitochondrial";
        case 10: return "euplotid nuclear";
        case 11: return "bacterial, archaeal and plant plastid";
        case 12: return "alternative yeast nuclear";
        case 13: return "ascidian mitochondrial";
        case 14: return "alternative flatworm mitochondrial";
        case 15: return "Blepharisma nuclear";
        case 16: return "chlorophycean mitochondrial";
        case 21: return "trematode mitochondrial";
        case 22: return "Scenedesmus obliquus mitochondrial";
        case 23: return "Thraustochytrium mitochondrial";
    }
    return "";
}

static std::string app_code_name(int t) {  // APP_CODE_NAMES
    switch (t) {
        case 1: return "Standard (nuclear genes of most organisms)";
        case 2: return "Vertebrate mitochondrial";
        case 3: return "Yeast mitochondrial";
        case 4: return "Mold, protozoan and coelenterate mitochondrial; Mycoplasma";
        case 5: return "Invertebrate mitochondrial";
        case 6: return "Ciliate, dasycladacean and Hexamita nuclear";
        case 9: return "Echinoderm and flatworm mitochondrial";
        case 10: return "Euplotid nuclear";
        case 11: return "Bacterial, archaeal and plant plastid";
        case 12: return "Alternative yeast nuclear";
        case 13: return "Ascidian mitochondrial";
        case 14: return "Alternative flatworm mitochondrial";
        case 15: return "Blepharisma nuclear";
        case 16: return "Chlorophycean mitochondrial";
        case 21: return "Trematode mitochondrial";
        case 22: return "Scenedesmus obliquus mitochondrial";
        case 23: return "Thraustochytrium mitochondrial";
    }
    return "";
}

static const char* as_x_name(char c) {
    switch (c) {
        case 'U': return "selenocysteine";
        case 'O': return "pyrrolysine";
        case 'B': return "D or N";
        case 'Z': return "E or Q";
        case 'J': return "I or L";
    }
    return nullptr;
}

// ------------------------------------------------------------------ small helpers

std::string fmt_n(long long x) {
    std::string s = std::to_string(x < 0 ? -x : x), out;
    int c = 0;
    for (size_t i = s.size(); i-- > 0;) {
        out += s[i];
        if (++c % 3 == 0 && i > 0) out += ',';
    }
    if (x < 0) out += '-';
    std::reverse(out.begin(), out.end());
    return out;
}

static std::string seqs(long long n) { return fmt_n(n) + " sequence" + (n == 1 ? "" : "s"); }

static std::string shlex_quote(const std::string& s) {
    if (s.empty()) return "''";
    bool safe = true;
    for (unsigned char c : s)
        if (c > 127 || !(isalnum(c) || c == '_' || c == '@' || c == '%' || c == '+' || c == '=' || c == ':' ||
                         c == ',' || c == '.' || c == '/' || c == '-')) {
            safe = false;
            break;
        }
    if (safe) return s;
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') out += "'\"'\"'";
        else out += c;
    }
    return out + "'";
}

std::string shlex_join(const std::vector<std::string>& args) {
    std::string out;
    for (size_t i = 0; i < args.size(); ++i) {
        if (i) out += ' ';
        out += shlex_quote(args[i]);
    }
    return out;
}

static std::string replace_all(std::string s, const std::string& a, const std::string& b) {
    if (a.empty()) return s;
    size_t p = 0;
    while ((p = s.find(a, p)) != std::string::npos) {
        s.replace(p, a.size(), b);
        p += b.size();
    }
    return s;
}

static std::string basename_of(const std::string& p) {  // os.path.basename
    size_t i = p.find_last_of("/\\");
    return i == std::string::npos ? p : p.substr(i + 1);
}

static std::string strip(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && isspace((unsigned char)s[a])) ++a;
    while (b > a && isspace((unsigned char)s[b - 1])) --b;
    return s.substr(a, b - a);
}

static std::string strip_version(const std::string& x) {  // re.sub(r"\.\d+$", "", x)
    size_t dot = x.rfind('.');
    if (dot == std::string::npos || dot + 1 >= x.size()) return x;
    for (size_t i = dot + 1; i < x.size(); ++i)
        if (!isdigit((unsigned char)x[i])) return x;
    return x.substr(0, dot);
}

static bool is_ambiguous_codon(const std::string& codon) {
    for (char ch : codon) {
        char c = (char)toupper((unsigned char)ch);
        if (c == 'U') c = 'T';
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return true;
    }
    return false;
}

static Json* field(Json& j, const std::string& k) {
    for (auto& p : j.o)
        if (p.first == k) return &p.second;
    return nullptr;
}

static std::string sget(const Json& j, const char* k) {
    const Json* v = j.get(k);
    return v && v->type == Json::String ? v->s : std::string();
}

static const std::vector<Json>& msgs_of(const Json& e) {
    static const std::vector<Json> none;
    const Json* m = e.get("messages");
    return m && m->type == Json::Array ? m->a : none;
}

static Json str_or_null(const std::string& s, bool present) { return present ? Json(s) : Json(); }

// digits at s[i..]: value, advance i; false if none
static bool read_num(const std::string& s, size_t& i, long long& v) {
    size_t j = i;
    v = 0;
    while (j < s.size() && isdigit((unsigned char)s[j])) v = v * 10 + (s[j++] - '0');
    if (j == i) return false;
    i = j;
    return true;
}

static bool lit(const std::string& s, size_t& i, const std::string& t) {
    if (s.compare(i, t.size(), t) != 0) return false;
    i += t.size();
    return true;
}

static std::string token(const std::string& s, size_t& i) {  // \S+
    size_t j = i;
    while (j < s.size() && !isspace((unsigned char)s[j])) ++j;
    std::string t = s.substr(i, j - i);
    i = j;
    return t;
}

// ^pepAlnPos (\d+): (\S+) does not correspond to (\S+)
static bool parse_mismatch(const std::string& w, long long& col, std::string& aa, std::string& codon) {
    size_t i = 0;
    if (!lit(w, i, "pepAlnPos ") || !read_num(w, i, col) || !lit(w, i, ": ")) return false;
    aa = token(w, i);
    if (aa.empty() || !lit(w, i, " does not correspond to ")) return false;
    codon = token(w, i);
    return !codon.empty();
}

// ^pepAlnPos (\d+): (\S+) (?:unknown AA type|\([^)]*\) is not checked)
static bool parse_unknown(const std::string& w, long long& col, std::string& aa) {
    size_t i = 0;
    if (!lit(w, i, "pepAlnPos ") || !read_num(w, i, col) || !lit(w, i, ": ")) return false;
    aa = token(w, i);
    if (aa.empty() || !lit(w, i, " ")) return false;
    if (w.compare(i, 15, "unknown AA type") == 0) return true;
    if (i < w.size() && w[i] == '(') {
        size_t j = w.find(')', i);
        return j != std::string::npos && w.compare(j, 16, ") is not checked") == 0;
    }
    return false;
}

static Json warning_message(const std::string& w, bool removed = false) {
    long long col;
    std::string aa, codon;
    Json m = Json::object();
    if (parse_mismatch(w, col, aa, codon)) {
        m.set("code", "mismatch");
        m.set("column", col);
        m.set("text", "Column " + std::to_string(col) + ": amino acid " + aa + " does not match codon " + codon + ".");
        m.set("residue", aa);
        m.set("codon", codon);
        if (removed) m.set("removed", true);
        return m;
    }
    if (parse_unknown(w, col, aa)) {
        const char* what = aa.size() == 1 ? as_x_name(aa[0]) : nullptr;
        m.set("code", "unknown_aa");
        m.set("column", col);
        m.set("residue", aa);
        m.set("text", "Column " + std::to_string(col) + ": '" + aa + "'" + (what ? std::string(" (") + what + ")" : "") +
                          " is not one of the 20 standard amino acids, so COAL did not check its codon (any codon "
                          "accepted).");
        return m;
    }
    m.set("code", "warning");
    m.set("text", "COAL warning: " + strip(w));
    return m;
}

static Json notice_message(const std::string& pid, const std::string& n) {
    std::string rest = n.compare(0, pid.size() + 2, pid + ": ") == 0 ? n.substr(pid.size() + 2) : n;
    Json m = Json::object();
    size_t i = 0;
    long long a = 0, b = 0;
    // ^protein and DNA are probably unrelated: (\d+) of (\d+) checkable codons do not match
    if (lit(rest, i, "protein and DNA are probably unrelated: ") && read_num(rest, i, a) && lit(rest, i, " of ") &&
        read_num(rest, i, b) && lit(rest, i, " checkable codons do not match")) {
        m.set("code", "suspect_unrelated");
        m.set("text", std::to_string(a) + " of the " + std::to_string(b) +
                          " codons COAL could check do not match this protein, so this is almost certainly the "
                          "wrong DNA for it (or the wrong genetic code). COAL still wrote its row, codon by "
                          "position, as PAL2NAL does: do not trust it. Run COAL with --strict (in the app: 'Refuse "
                          "probably-wrong genes') to refuse such genes.");
        return m;
    }
    i = 0;
    const std::string guess = "position of the gene in the DNA is a guess: ";
    if (lit(rest, i, guess)) {
        for (size_t p = i; p <= rest.size(); ++p) {  // lazy .*? then the first alternative that matches
            size_t q = p;
            long long nchk = 0, nres = 0, shrt = 0;
            if (lit(rest, q, "only ") && read_num(rest, q, nchk) && lit(rest, q, " of ") && read_num(rest, q, nres) &&
                lit(rest, q, " codons have no ambiguous bases")) {
                m.set("code", "suspect_position");
                m.set("text", "Only " + std::to_string(nchk) + " of this gene's " + std::to_string(nres) +
                                  " codons have no unknown bases (N), and no 10-residue stretch of the protein was "
                                  "found in the DNA, so COAL cannot tell where the gene starts in the longer DNA. "
                                  "Its row was taken from the start of the DNA (as PAL2NAL does) and may be "
                                  "shifted: do not trust it. Run COAL with --strict (in the app: 'Refuse "
                                  "probably-wrong genes') to refuse such genes.");
                return m;
            }
            q = p;
            if (lit(rest, q, "the gene has only ") && read_num(rest, q, shrt) && lit(rest, q, " codons")) {
                m.set("code", "suspect_position");
                m.set("text", "The gene has only " + std::to_string(shrt) +
                                  " codons, and no 10-residue stretch of the protein was found in the DNA, so COAL "
                                  "cannot tell where the gene starts in the longer DNA. Its row was taken from the "
                                  "start of the DNA (as PAL2NAL does) and may be shifted: do not trust it. Run COAL "
                                  "with --strict (in the app: 'Refuse probably-wrong genes') to refuse such genes.");
                return m;
            }
        }
    }
    i = 0;
    long long nbad = 0, nres2 = 0, namb = -1;
    if (lit(rest, i, "no 10-residue segment of the protein was found in the DNA; codons were assigned by position (") &&
        read_num(rest, i, nbad) && lit(rest, i, " of ") && read_num(rest, i, nres2) && lit(rest, i, " ")) {
        lit(rest, i, "checkable codons ");
        bool ok = lit(rest, i, "do not match");
        size_t save = i;
        if (ok && !(lit(rest, i, "; ") && read_num(rest, i, namb) && lit(rest, i, " contain N or other ambiguous bases"))) {
            i = save;
            namb = -1;
        }
        if (ok && lit(rest, i, ")")) {
            std::string amb = namb >= 0 ? "; " + std::to_string(namb) + " contain N or other ambiguous bases" : "";
            std::string ask = namb >= 0 && nbad == 0  // Python: `namb and int(nbad) == 0` (namb = the group's text)
                                  ? "Every codon COAL could check matches, so the placement looks right."
                                  : "Is this the right DNA for this protein, and the right genetic code?";
            m.set("code", "no_anchor");
            m.set("text", "No 10-residue stretch of the protein was found in the DNA, so codons were assigned by "
                          "position (" + std::to_string(nbad) + " of " + std::to_string(nres2) +
                              " checkable codons do not match" + amb + "). " + ask);
            return m;
        }
    }
    m.set("code", "warning");
    m.set("text", "COAL warning: " + strip(rest));
    return m;
}

static Json failure_message(const std::string& pid, const std::string& code, const std::string& message) {
    std::string rest = message.compare(0, pid.size() + 2, pid + ": ") == 0 ? message.substr(pid.size() + 2) : message;
    Json m = Json::object();
    m.set("code", code);
    m.set("text", rest);
    return m;
}

static bool is_suspect(const std::string& c) { return c == "suspect_unrelated" || c == "suspect_position"; }

static std::string gene_kind(const Json& e) {
    std::string st = sget(e, "status");
    const auto& msgs = msgs_of(e);
    if (st == "failed") {
        for (auto& m : msgs) {
            std::string c = sget(m, "code");
            if (is_suspect(c) || c == "suspect_pairing") return "flagged";
        }
        return "excluded";
    }
    if (st == "warning") {
        bool any = false;
        for (auto& m : msgs)
            if (sget(m, "code") == "mismatch") {
                any = true;
                if (!is_ambiguous_codon(sget(m, "codon"))) return "mismatch";
            }
        return any ? "unchecked" : "other";
    }
    return st;
}

static char kind_char(const std::string& k) {
    if (k == "ok") return 'o';
    if (k == "mismatch") return 'm';
    if (k == "unchecked") return 'u';
    if (k == "other") return 'x';
    if (k == "flagged") return 'f';
    if (k == "excluded") return 'e';
    return 'n';
}

static std::vector<const Json*> real_mismatches(const Json& e) {
    std::vector<const Json*> out;
    for (auto& m : msgs_of(e))
        if (sget(m, "code") == "mismatch" && !is_ambiguous_codon(sget(m, "codon"))) out.push_back(&m);
    return out;
}

static bool has_code(const Json& e, const char* code) {
    for (auto& m : msgs_of(e))
        if (sget(m, "code") == code) return true;
    return false;
}

static bool fixable(const Json& e) {
    std::string k = gene_kind(e);
    if (k == "mismatch") return !has_code(e, "possible_frameshift");
    if (k == "flagged") return has_code(e, "suspect_unrelated") || !real_mismatches(e).empty();
    if (k == "not_checked") return has_code(e, "not_written") ? !real_mismatches(e).empty() : true;
    return k == "excluded";
}

static long long code_distance(int a, int b) {  // the first part of the app's _code_distance (codons that differ)
    const GeneticCode* ca = get_code(a);
    const GeneticCode* cb = get_code(b);
    long long d = 0;
    for (int c = 0; c < 64; ++c) d += ca->aa[c] != cb->aa[c];
    return d;
}

static void annotate_mismatch_pattern(const std::string& row, Json& e) {
    std::map<long long, long long> index;
    long long k = 0;
    for (size_t c = 0; c < row.size(); ++c)
        if (row[c] != '-' && row[c] != '.') index[(long long)c + 1] = k++;
    std::set<long long> idxs;
    bool all_m = true;
    Json* msgs = field(e, "messages");
    if (!msgs) return;
    std::vector<size_t> mm;
    for (size_t i = 0; i < msgs->a.size(); ++i) {
        const Json& m = msgs->a[i];
        if (sget(m, "code") != "mismatch" || is_ambiguous_codon(sget(m, "codon"))) continue;
        mm.push_back(i);
        const Json* col = m.get("column");
        if (col && col->type == Json::Int && index.count(col->i)) idxs.insert(index[col->i]);
        std::string r = sget(m, "residue");
        if (!(r.size() == 1 && toupper((unsigned char)r[0]) == 'M')) all_m = false;
    }
    std::vector<long long> idx(idxs.begin(), idxs.end());
    if (idx.empty()) return;
    if (idx.size() == 1 && idx[0] == 0 && all_m) {
        for (size_t i : mm) msgs->a[i].set("start_codon", true);
        return;
    }
    std::vector<std::vector<long long>> runs;
    // push_back, not `cur = {x}`: GCC 14's -Warray-bounds misreads the initializer_list temporary of the assignment
    // as an out-of-bounds memcpy (a false positive that broke the manylinux wheel build under -Werror)
    std::vector<long long> cur;
    cur.push_back(idx[0]);
    for (size_t i = 1; i < idx.size(); ++i) {
        if (idx[i] - cur.back() <= 2) cur.push_back(idx[i]);
        else {
            runs.push_back(cur);
            cur.clear();
            cur.push_back(idx[i]);
        }
    }
    runs.push_back(cur);
    size_t bi = 0;
    for (size_t i = 1; i < runs.size(); ++i)
        if (runs[i].size() > runs[bi].size()) bi = i;
    const auto& best = runs[bi];
    long long span = best.back() - best.front() + 1;
    bool to_end = best.back() >= k - 2;
    long long nb = (long long)best.size();
    if (nb >= 3 && (double)nb / (double)span >= 0.6 && 2 * nb >= (long long)idx.size() && (to_end || nb >= 4)) {
        long long c0 = 0, c1 = 0;
        for (auto& [c, v] : index) {
            if (v == best.front()) c0 = c;
            if (v == best.back()) c1 = c;
        }
        std::string where = to_end ? "to the end of the gene" : "up to column " + std::to_string(c1);
        Json m = Json::object();
        m.set("code", "possible_frameshift");
        m.set("column", c0);
        m.set("end_column", c1);
        m.set("text", "The mismatches form one run, from column " + std::to_string(c0) + " " + where +
                          ": this looks like a missing or extra base (a frameshift) or a small insertion/deletion "
                          "in the DNA near column " + std::to_string(c0) +
                          ", not a genetic-code problem. Check the DNA there (a sequencing or assembly error?).");
        msgs->a.push_back(m);
    }
}

// ------------------------------------------------------------------ where a failed gene breaks (app viewer rule)

static bool aa_ok(char aa, const std::string& dna, long long pos, const GeneticCode& code, bool first) {
    if (pos < 0 || pos + 3 > (long long)dna.size()) return false;
    if (aa == 'X' || aa == 'B' || aa == 'Z' || aa == 'J' || aa == 'O' || aa == 'U') return true;
    const uint8_t* iu = iupac_table();
    uint8_t a = iu[(unsigned char)dna[(size_t)pos]], b = iu[(unsigned char)dna[(size_t)pos + 1]],
            c = iu[(unsigned char)dna[(size_t)pos + 2]];
    if ((a | b | c) & 0xF0) return false;
    uint64_t mask = iupac_triplet_table()[a * 225 + b * 15 + c];
    if (!mask) return false;
    uint64_t ok = code.codons[(unsigned char)aa];
    if (first && aa == 'M') ok |= code.starts;
    return (mask & ~ok) == 0;
}

// As report.py break_event: the earliest credible start, a shift must hold for the next 20
// residues and may include a lost / extra codon (+-4, +-5), "at" vs "near" column, "the DNA ends" only after a
// credible walk.
static bool break_event(const std::string& row, const std::string& raw_dna, int table, Json& ev) {
    const GeneticCode& code = *get_code(table);
    std::string dna;
    for (char ch : raw_dna) {
        if (!isalpha((unsigned char)ch)) continue;
        char c = (char)toupper((unsigned char)ch);
        dna += c == 'U' ? 'T' : c;
    }
    std::vector<std::pair<long long, char>> res;
    for (size_t i = 0; i < row.size(); ++i)
        if (isalpha((unsigned char)row[i])) res.push_back({(long long)i + 1, (char)toupper((unsigned char)row[i])});
    if (res.empty() || dna.size() < 3) return false;
    const long long nres = (long long)res.size(), nd = (long long)dna.size();
    long long want = std::min<long long>(nres, 12), allowed = std::max<long long>(0, want / 4);
    const long long need = std::max(std::min<long long>(3, nres), want - allowed);
    long long best_off = 0, best_score = -1;
    for (long long off = 0; off < nd - 2; ++off) {
        long long score = 0;
        for (long long k = 0; k < want; ++k) {
            if (aa_ok(res[(size_t)k].second, dna, off + 3 * k, code, k == 0)) ++score;
            else if (k + 1 - score > allowed) break;
        }
        if (score >= need) {  // the earliest credible start
            best_off = off;
            best_score = score;
            break;
        }
        if (score > best_score) {
            best_off = off;
            best_score = score;
        }
    }
    auto event = [&](const char* c, long long col, const std::string& text) {
        ev = Json::object();
        ev.set("code", c);
        ev.set("column", col);
        ev.set("text", text);
    };
    if (best_score < need) {
        event("inconsistent", res[0].first, "The DNA does not spell the start of this protein in any reading frame.");
        return true;
    }
    auto hits = [&](long long start, long long k0, long long n) {
        long long h = 0;
        for (long long j = 0; j < n && k0 + j < nres; ++j)
            h += aa_ok(res[(size_t)(k0 + j)].second, dna, start + 3 * j, code, false);
        return h;
    };
    long long pos = best_off, k = 0, matched = 0, last_good = -1;
    std::map<long long, std::string> flags;
    while (k < nres) {
        long long col = res[(size_t)k].first;
        char aa = res[(size_t)k].second;
        if (pos + 3 > nd) {
            if (k && 5 * matched >= 4 * k) {
                event("missing_dna", col, "The DNA ends before the protein does (at column " + std::to_string(col) + ").");
                return true;
            }
            return false;
        }
        if (aa_ok(aa, dna, pos, code, k == 0)) {
            flags[col] = "";
            pos += 3;
            ++k;
            ++matched;
            last_good = k - 1;
            continue;
        }
        long long ahead = 0;
        for (long long j = 1; j <= 5 && k + j < nres; ++j)
            ahead += aa_ok(res[(size_t)(k + j)].second, dna, pos + 3 * j, code, false);
        long long window = std::min<long long>(5, nres - k - 1);
        if (window && ahead >= window - 1) {
            flags[col] = "mismatch";
            pos += 3;
            ++k;
            continue;
        }
        int shifted = 0;
        bool have = false;
        const long long n = std::min<long long>(10, nres - k - 1), n2 = std::min<long long>(20, nres - k - 1);
        for (int s : {-1, 1, -2, 2, -4, 4, -5, 5}) {
            if (n <= 0) break;
            if (hits(pos + 3 + s, k + 1, n) >= std::max<long long>(1, (long long)(0.8 * (double)n)) &&
                (double)hits(pos + 3 + s, k + 1, n2) >= 0.8 * (double)n2) {
                shifted = s;
                have = true;
                break;
            }
        }
        if (have) {
            const char* what = shifted == -1   ? "one base is missing"
                               : shifted == -2 ? "two bases are missing"
                               : shifted == 1  ? "there is one extra base"
                               : shifted == 2  ? "there are two extra bases"
                               : shifted == -4 ? "four bases are missing"
                               : shifted == -5 ? "five bases are missing"
                               : shifted == 4  ? "there are four extra bases"
                                               : "there are five extra bases";
            bool exact = last_good == k - 1 && aa_ok(res[(size_t)(k + 1)].second, dna, pos + 3 + shifted, code, false);
            ev = Json::object();
            ev.set("code", "frameshift");
            ev.set("column", col);
            ev.set("shift", (long long)shifted);
            ev.set("text", std::string("Frameshift ") + (exact ? "at" : "near") + " column " + std::to_string(col) +
                               ": " + what + " in the DNA, so from here on the codons are read out of step.");
            return true;
        }
        flags[col] = "mismatch";
        pos += 3;
        ++k;
    }
    long long bad_after = 0, first = -1;
    for (long long i = nres - 1; i >= 0; --i) {
        auto it = flags.find(res[(size_t)i].first);
        std::string f = it == flags.end() ? "missing" : it->second;
        bad_after += (f == "mismatch" || f == "missing");
        long long nn = nres - i;
        if (f == "mismatch" && nn >= 10 && 2 * bad_after > nn) first = i;
    }
    if (first >= 0) {
        long long col = res[(size_t)first].first;
        event("inconsistent", col, "From column " + std::to_string(col) + " on, the DNA no longer spells this protein.");
        return true;
    }
    return false;
}

// ------------------------------------------------------------------ genetic-code search (report.py search_codes)

// One gene fits a code: it threads, with no definite (A/C/G/T) mismatch and no suspect notice
static bool gene_fits(const HintGene& g, const GeneticCode& gc) {
    Threaded it;
    try {
        it = thread_one(g.id, g.row, g.dna, gc);
    } catch (const CoalError&) {
        return false;
    }
    long long col;
    std::string aa, codon;
    for (auto& w : it.warnings)
        if (parse_mismatch(w, col, aa, codon) && !is_ambiguous_codon(w.substr(w.rfind(' ') + 1))) return false;
    for (auto& x : it.notices)
        if (is_suspect(sget(notice_message(g.id, x), "code"))) return false;
    return true;
}

// The hint's search, shared by the end-of-run hint and the early notice; a code is dropped
// as soon as it can no longer fit as many genes as the best so far (same result as testing all). tied empty = none.
// memo: (sequence, code) -> fits, shared by the early notice and the end-of-run hint
static CodeSearch search_codes(const std::vector<long long>& tested, const std::vector<long long>& problems,
                               const std::map<long long, HintGene>& data, int table, FitMemo& memo) {
    const long long nt = (long long)tested.size();
    auto bad_under = [&](int c, long long at_least, bool& ok) {
        std::set<long long> bad;
        ok = true;
        const GeneticCode* gc = get_code(c);
        for (long long k : tested) {
            auto key = std::make_pair(k, c);
            auto hit = memo.find(key);
            bool fits = hit != memo.end() ? hit->second : (memo[key] = gene_fits(data.at(k), *gc));
            if (!fits) {
                bad.insert(k);
                if ((long long)bad.size() > nt - at_least) {
                    ok = false;
                    break;
                }
            }
        }
        return bad;
    };
    CodeSearch out;
    bool ok_now;
    std::set<long long> bad_now = bad_under(table, 0, ok_now);
    out.fits_now = nt - (long long)bad_now.size();
    std::vector<std::pair<long long, int>> results;
    long long best = out.fits_now + 1;
    for (int c : supported_tables()) {
        if (c == table) continue;
        bool ok;
        std::set<long long> bad = bad_under(c, best, ok);
        if (!ok) continue;
        long long fit = nt - (long long)bad.size(), fixed = 0;
        for (long long k : problems) fixed += bad_now.count(k) && !bad.count(k);
        if (fixed && fit > out.fits_now) {
            results.push_back({fit, c});
            best = std::max(best, fit);
        }
    }
    if (results.empty()) return out;
    for (auto& r : results) out.top = std::max(out.top, r.first);
    for (auto& r : results)
        if (r.first == out.top) out.tied.push_back(r.second);
    std::stable_sort(out.tied.begin(), out.tied.end(), [&](int a, int b) {
        return std::make_pair(code_distance(table, a), a) < std::make_pair(code_distance(table, b), b);
    });
    // The standard code 1 first when it is among the tied codes (as report.py)
    out.standard = out.tied.size() > 1 && std::find(out.tied.begin(), out.tied.end(), 1) != out.tied.end();
    if (out.standard) {
        out.tied.erase(std::find(out.tied.begin(), out.tied.end(), 1));
        out.tied.insert(out.tied.begin(), 1);
    }
    return out;
}

static std::string do_(long long k) {
    return k == 0 ? std::string("none do") : k == 1 ? std::string("only 1 does") : std::to_string(k) + " do";
}

// As report.py early_code_notice: the suggestion said before the per-codon warnings, from the
// first HINT_SAMPLE genes. Empty: no code does better.
std::string early_code_notice(const RowsFn& rows, const DnaFn& dna, int table, long long total, FitMemo& memo) {
    std::map<long long, HintGene> data;
    std::vector<long long> tested, problems;
    const GeneticCode* gc = get_code(table);
    rows([&](long long k, const std::string& pid, const std::string& row) {
        if (k >= HINT_SAMPLE) return false;
        data[k] = HintGene{pid, row, dna(k)};
        tested.push_back(k);
        bool fits = memo[std::make_pair(k, table)] = gene_fits(data[k], *gc);
        if (!fits) problems.push_back(k);
        return k + 1 < HINT_SAMPLE;
    });
    if (problems.empty()) return "";
    CodeSearch f = search_codes(tested, problems, data, table, memo);
    if (f.tied.empty()) return "";
    const int c = f.tied[0];
    const long long nt = (long long)tested.size();
    std::string who = nt == total ? (f.top == total ? "all " + std::to_string(total) + " sequences match"
                                                    : std::to_string(f.top) + " of " + std::to_string(total) +
                                                          " sequences match")
                                  : std::to_string(f.top) + " of the first " + std::to_string(nt) + " sequences match";
    return "wrong genetic code? Under genetic code " + std::to_string(c) + " (" + app_code_name(c) + "), " + who +
           " their protein; under code " + std::to_string(table) + ", " + do_(f.fits_now) + ". Re-run with -codontable " +
           std::to_string(c) + ".";
}

// ------------------------------------------------------------------ RunReport

RunReport::RunReport(const std::vector<std::string>& argv, const std::string& pep, const std::vector<std::string>& nucs,
                     const std::string& out, const std::string& fmt, int table, bool nogap, bool nomismatch,
                     bool keep_going, bool keep_entries_, bool strict, bool blockonly)
    : keep_entries(keep_entries_), t0_(std::chrono::steady_clock::now()), argv_(argv), pep_(pep), out_(out), fmt_(fmt),
      nucs_(nucs), table_(table), nogap_(nogap), nomismatch_(nomismatch), keep_going_(keep_going), strict_(strict),
      blockonly_(blockonly) {}

void RunReport::set_inputs(std::vector<std::string>&& ids) {
    ids_ = std::move(ids);
    has_ids_ = true;
}

void RunReport::set_dna_inputs(const std::vector<long long>& dna_counts, const std::vector<std::string>& dna_file_ids) {
    dna_counts_ = dna_counts;
    has_dna_counts_ = true;
    dna_file_ids_ = dna_file_ids;
    has_dna_file_ids_ = true;
}

void RunReport::set_pairing(const std::string& mode, const std::vector<std::string>& per_seq,
                            const std::vector<size_t>& pairing, const std::vector<std::string>& notes) {
    pairing_ = mode;
    seq_pairing_ = per_seq;
    pair_idx_ = pairing;
    has_dna_ids_ = has_dna_file_ids_;
    for (const std::string& n : notes) {
        size_t i = 0;
        long long k = 0, tot = 0;
        if (lit(n, i, "only ") && read_num(n, i, k) && lit(n, i, " of ") && read_num(n, i, tot) &&
            lit(n, i, " protein IDs are DNA IDs; ")) {
            std::string rest = n.substr(i);
            partial_ = true;
            std::string eg;
            size_t p = rest.find("protein '");
            while (p != std::string::npos) {  // re.search(r"protein '([^']*)' gets DNA '([^']*)'")
                size_t a = p + 9, b = rest.find('\'', a);
                if (b != std::string::npos && rest.compare(b, 12, "' gets DNA '") == 0) {
                    size_t c = b + 12, d = rest.find('\'', c);
                    if (d != std::string::npos) {
                        eg = " (for example protein '" + rest.substr(a, b - a) + "' got DNA '" + rest.substr(c, d - c) + "')";
                        break;
                    }
                }
                p = rest.find("protein '", p + 1);
            }
            Json h = Json::object();
            h.set("code", "partial_id_match");
            h.set("text", "Only " + std::to_string(k) + " of " + std::to_string(tot) +
                              " protein IDs were found in the DNA file, so COAL paired the sequences by their order "
                              "in the files instead" + eg + ". Some proteins are probably paired with the wrong DNA: "
                              "make the names match (or put both files in the same order) and run again.");
            pair_hints_.push_back(h);
        } else if (mode == "normalised" && !n.empty()) {
            std::string cap = n;
            cap[0] = (char)toupper((unsigned char)cap[0]);
            Json h = Json::object();
            h.set("code", "normalised_ids");
            h.set("text", "Some protein IDs match their DNA only after tidying them (a '.1'-style version suffix, "
                          "upper/lower case or database fields such as 'sp|'), so COAL paired those sequences by "
                          "their tidied IDs. PAL2NAL pairs by file order here, which would have given some proteins "
                          "the DNA of another gene. " + cap + ".");
            pair_hints_.push_back(h);
        }
    }
}

void RunReport::set_kind(long long k, char ch) {
    if (k < (long long)kinds_.size()) kinds_[(size_t)k] = ch;
    else kinds_ += ch;
}

static void add_sorted(std::vector<long long>& v, long long x) {
    auto it = std::lower_bound(v.begin(), v.end(), x);
    if (it == v.end() || *it != x) v.insert(it, x);
}

const std::string* RunReport::add_sequence(long long k, const std::string& pid, const std::vector<std::string>& warnings,
                                    const std::vector<std::string>& notices, const std::string& fail_code,
                                    const std::string& fail_msg) {
    while ((long long)kinds_.size() < k) kinds_ += 'n';
    std::vector<Json> msgs;
    std::string status;
    if (!fail_code.empty()) {
        status = "failed";
        msgs.push_back(failure_message(pid, fail_code, fail_msg));
    } else {
        long long col;
        std::string aa, codon;
        if (nomismatch_) {
            for (auto& w : warnings)
                if (parse_mismatch(w, col, aa, codon)) {
                    msgs.push_back(warning_message(w, true));
                    add_sorted(cols_, col);
                    if (!is_ambiguous_codon(codon)) add_sorted(real_cols_, col);
                }
            for (auto& w : warnings)
                if (!parse_mismatch(w, col, aa, codon)) msgs.push_back(warning_message(w));
        } else {
            for (auto& w : warnings) msgs.push_back(warning_message(w));
        }
        for (auto& n : notices) msgs.push_back(notice_message(pid, n));
        bool has_d = has_dna_ids_ && k < (long long)pair_idx_.size();
        std::string d = has_d ? dna_id((size_t)k) : "";
        bool mispaired = pairing_ == "order" && partial_ && has_d && d != pid;
        if (mispaired) {
            if (!paired_set_built_) {
                for (size_t i = 0; i < pair_idx_.size(); ++i) paired_set_.insert(dna_id(i));
                paired_set_built_ = true;
            }
            bool in_dna = paired_set_.count(pid) != 0;
            Json m = Json::object();
            m.set("code", "suspect_pairing");
            m.set("dna_id", d);
            m.set("text", "Paired with DNA '" + d + "' only by its position in the files" +
                              (in_dna ? " (this protein's own name is also in the DNA file)" : "") +
                              ": check that this is the right DNA.");
            msgs.push_back(m);
        }
        bool sus = false, mm = false;
        for (auto& m : msgs) {
            sus = sus || is_suspect(sget(m, "code"));
            mm = mm || sget(m, "code") == "mismatch";
        }
        if (sus || (mispaired && mm)) {
            status = "failed";
            std::stable_sort(msgs.begin(), msgs.end(), [](const Json& a, const Json& b) {
                auto key = [](const Json& m) {
                    std::string c = sget(m, "code");
                    return std::make_pair(!is_suspect(c), c != "suspect_pairing");
                };
                return key(a) < key(b);
            });
        } else {
            status = msgs.empty() ? "ok" : "warning";
        }
    }
    if (status == "ok") {
        set_kind(k, 'o');
        return nullptr;
    }
    Json e = Json::object();
    e.set("id", pid);
    e.set("index", k);
    e.set("dna_id", has_dna_ids_ && k < (long long)pair_idx_.size() ? Json(dna_id((size_t)k)) : Json());
    e.set("pairing", !seq_pairing_.empty() ? Json(seq_pairing_[(size_t)k]) : str_or_null(pairing_, !pairing_.empty()));
    e.set("status", status);
    Json arr = Json::array();
    for (auto& m : msgs) arr.push(m);
    e.set("messages", arr);
    store(k, e);
    set_kind(k, kind_char(gene_kind(e)));
    return &entries_[k].json;
}

void RunReport::store(long long k, const Json& e) {
    Kept& kp = entries_[k];
    kp.json = json_dump(e);
    bool col = false;
    for (auto& m : msgs_of(e)) col = col || (m.get("column") && m.get("column")->type == Json::Int);
    kp.failed_nocol = sget(e, "status") == "failed" && !col;
    kp.real_mm = !real_mismatches(e).empty();
    kp.unknown_aa = has_code(e, "unknown_aa");
}

Json RunReport::entry(long long k) const { return json_parse(entries_.at(k).json); }

void RunReport::fail(const std::string& code, const std::string& message, bool has_seq_id) {
    error_ = Json::object();
    error_.set("code", code);
    error_.set("text", message);
    error_.set("detail", message);
    has_error_ = true;
    if (!has_seq_id) {
        Json e = Json::object();
        e.set("code", code);
        e.set("text", message);
        errors_.push_back(e);
    }
}

void RunReport::set_output(const std::string& path, const std::string& sha, long long n) {
    output_ = Json::object();
    output_.set("path", str_or_null(path, !path.empty()));
    output_.set("format", fmt_);
    output_.set("sha256", str_or_null(sha, !sha.empty()));
    output_.set("sequences", n);
    has_output_ = true;
}

Json RunReport::checkpoint_state() const {
    Json s = Json::object();
    s.set("kinds", kinds_);
    Json c = Json::array(), r = Json::array();
    for (long long x : cols_) c.push(x);
    for (long long x : real_cols_) r.push(x);
    s.set("cols", c);
    s.set("real_cols", r);
    return s;
}

void RunReport::restore(const Json& state, const std::vector<Json>& nonok, long long completed) {
    const Json* k = state.get("kinds");
    kinds_ = k && k->type == Json::String ? k->s : std::string((size_t)completed, 'o');
    if ((long long)kinds_.size() > completed) kinds_.resize((size_t)completed);
    cols_.clear();
    real_cols_.clear();
    if (const Json* c = state.get("cols"))
        for (auto& x : c->a) add_sorted(cols_, x.as_int());
    if (const Json* c = state.get("real_cols"))
        for (auto& x : c->a) add_sorted(real_cols_, x.as_int());
    entries_.clear();
    for (auto& e : nonok) {
        long long idx = e.get("index") ? e.get("index")->as_int() : -1;
        if (idx >= 0 && idx < completed) store(idx, e);
    }
    resumed_from_ = completed;
}

long long RunReport::total() const { return has_ids_ ? (long long)ids_.size() : (long long)kinds_.size(); }

long long RunReport::first_fail() const {
    size_t p = kinds_.find('e');
    return p == std::string::npos ? -1 : (long long)p;
}

Json RunReport::entry_view(long long k, int aborted) const {
    bool ab = aborted < 0 ? !has_output_ : aborted != 0;
    Json e;
    auto it = entries_.find(k);
    if (it != entries_.end()) {
        e = json_parse(it->second.json);
    } else {
        std::string st = k < (long long)kinds_.size() && kinds_[(size_t)k] == 'o' ? "ok" : "not_checked";
        e = Json::object();
        e.set("id", has_ids_ && k < (long long)ids_.size() ? Json(ids_[(size_t)k]) : Json());
        e.set("index", k);
        bool paired = !pairing_.empty();
        e.set("dna_id", paired && has_dna_ids_ && k < (long long)pair_idx_.size() ? Json(dna_id((size_t)k)) : Json());
        e.set("pairing", paired ? (!seq_pairing_.empty() ? Json(seq_pairing_[(size_t)k]) : Json(pairing_)) : Json());
        e.set("status", st);
        e.set("messages", Json::array());
    }
    if (!ab || (k < (long long)kinds_.size() && kinds_[(size_t)k] == 'e')) return e;
    long long ff = first_fail();
    Json msgs = Json::array();
    if (ff >= 0) {
        std::string fid = sget(entry(ff), "id");
        Json m = Json::object();
        m.set("code", k < ff ? "not_written" : "not_reached");
        m.set("failed_id", fid);
        msgs.push(m);
        if (k < ff)
            for (auto& x : msgs_of(e))
                if (sget(x, "code") != "suspect_pairing") msgs.push(x);
    }
    e.set("status", "not_checked");
    e.set("messages", msgs);
    return e;
}

std::map<std::string, long long> RunReport::counts() const {
    // By what COAL found for each sequence it processed (its kind: failed for f/e, warning for m/u/x, ok for o), plus
    // not_checked for the ones it never reached (report.py counts): in a run that stopped at an error these are the
    // numbers of the "coal: summary:" lines, not the entry views' statuses (which mark the sequences checked before
    // the error not_checked / not_written too, as no output holds them)
    std::map<std::string, long long> c{{"ok", 0}, {"warning", 0}, {"failed", 0}, {"not_checked", 0}};
    for (char ch : kinds_) {
        switch (ch) {
            case 'o': c["ok"]++; break;
            case 'm': case 'u': case 'x': c["warning"]++; break;
            case 'f': case 'e': c["failed"]++; break;
            default: c["not_checked"]++;
        }
    }
    c["not_checked"] += std::max<long long>(0, total() - (long long)kinds_.size());
    return c;
}

void RunReport::finish(const RowsFn& rows, const DnaFn& dna, bool aborted) {
    finished = true;
    try {
        finish_impl(rows, dna, aborted);
    } catch (const std::exception& e) {
        // extras must never hide the real result. COAL_DEBUG=1 says what failed (report.py)
        if (std::getenv("COAL_DEBUG"))
            std::fprintf(stderr, "coal: debug: the report's finish step failed (the run goes on): %s\n", e.what());
    } catch (...) {
        if (std::getenv("COAL_DEBUG"))
            std::fprintf(stderr, "coal: debug: the report's finish step failed (the run goes on)\n");
    }
}

void RunReport::finish_impl(const RowsFn& rows, const DnaFn& dna, bool aborted) {
    const long long n = total();
    std::vector<long long> failed, marks;
    for (auto& [k, kp] : entries_) {
        if (kp.failed_nocol) failed.push_back(k);
        if (kp.real_mm) marks.push_back(k);
    }
    std::set<long long> fset(failed.begin(), failed.begin() + (std::ptrdiff_t)std::min<size_t>(failed.size(), BREAK_SAMPLE));
    std::set<long long> mset(marks.begin(), marks.end());
    std::set<long long> need = fset;
    need.insert(mset.begin(), mset.end());
    if (!need.empty()) {
        rows([&](long long k, const std::string&, const std::string& row) {
            if (need.count(k)) {
                Json e = entry(k);
                if (std::find(failed.begin(), failed.end(), k) != failed.end()) {  // (all failed, as Python)
                    Json ev;
                    if (break_event(row, dna(k), table_, ev)) {
                        Json* ms = field(e, "messages");
                        ms->a.insert(ms->a.begin(), ev);
                    }
                }
                if (mset.count(k)) annotate_mismatch_pattern(row, e);
                store(k, e);
            }
            return k + 1 < n;
        });
    }
    if (!want_hint) return;  // -q without --report/--json never shows it
    // One entry view at a time (keeping all n of them took ~1.1 KB per sequence: 114 MB at 100,000, and keeping the
    // problem ones ~2 KB per mismatching sequence); a clean sequence of a run that was not aborted is never fixable,
    // so its view is not even built. Of the problems' views only the start-codon summary is kept: whether every
    // definite mismatch is a first-Met one, and the first such codon (the text below).
    std::vector<long long> problems;
    bool any_start = false, all_start = true;
    std::string first_codon;
    for (long long k = 0; k < n; ++k) {
        if (!aborted && k < (long long)kinds_.size() && kinds_[(size_t)k] == 'o' && !entries_.count(k)) continue;
        Json v = entry_view(k, aborted ? 1 : 0);
        if (fixable(v)) {
            problems.push_back(k);
            for (auto* m : real_mismatches(v)) {
                if (!any_start) first_codon = sget(*m, "codon");
                any_start = true;
                const Json* sc = m->get("start_codon");
                all_start = all_start && sc && sc->type == Json::Bool && sc->b;
            }
        }
    }
    all_start = all_start && any_start;
    if (problems.empty()) return;
    std::set<long long> pset(problems.begin(), problems.end());
    std::vector<long long> tested = problems;
    for (long long k = 0; k < n && (long long)tested.size() < HINT_SAMPLE; ++k)
        if (!pset.count(k)) tested.push_back(k);
    if ((long long)tested.size() > HINT_SAMPLE) tested.resize(HINT_SAMPLE);
    if (tested.empty()) return;
    std::set<long long> want(tested.begin(), tested.end());
    std::map<long long, HintGene> data;
    rows([&](long long k, const std::string& pid, const std::string& row) {
        if (want.count(k)) data[k] = HintGene{pid, row, dna(k)};
        return data.size() < want.size();
    });
    const long long nt = (long long)tested.size();
    CodeSearch found = search_codes(tested, problems, data, table_, fit_memo);
    if (found.tied.empty()) return;
    const std::vector<int>& tied = found.tied;
    const bool standard = found.standard;
    const long long top = found.top, fits_now = found.fits_now;
    int c = tied[0];
    std::string name = app_code_name(c);
    std::string text;
    std::string T = std::to_string(table_), C = std::to_string(c);
    if (nt == n) {
        std::string who = top == n ? "all " + std::to_string(n) + " sequences match"
                                   : std::to_string(top) + " of " + std::to_string(n) + " sequences match";
        text = "Under genetic code " + C + " (" + name + "), " + who + " their protein; under code " + T + ", " +
               do_(fits_now) + ".";
    } else {
        std::string which = (long long)problems.size() >= nt
                                ? "a sample of the " + std::to_string(problems.size()) + " with problems"
                                : "all " + std::to_string(problems.size()) + " with problems, then others";
        text = "COAL tested " + std::to_string(nt) + " of the " + std::to_string(n) + " sequences (" + which +
               "): under genetic code " + C + " (" + name + "), " + std::to_string(top) +
               " of them match their protein; under code " + T + ", " + do_(fits_now) + ".";
    }
    if (all_start) {
        size_t np = pset.size();
        text = "The only mismatches are at the first codon of " + std::to_string(np) + " sequence" + (np != 1 ? "s" : "") +
               " (for example " + first_codon +
               " where the protein starts with M): an alternative start codon, which code " + C + " accepts. " + text;
    }
    if (tied.size() > 1) {
        std::string others;
        if (tied.size() > 2) {
            std::string mid;
            for (size_t i = 1; i + 1 < tied.size(); ++i) mid += (i > 1 ? ", " : "") + std::to_string(tied[i]);
            others = mid + " and " + std::to_string(tied.back());
        } else {
            others = std::to_string(tied[1]);
        }
        bool closer = code_distance(table_, c) < code_distance(table_, tied[1]);
        text += std::string(" Code") + (tied.size() > 2 ? "s" : "") + " " + others + " fit" + (tied.size() > 2 ? "" : "s") + " " +
                (standard ? std::string("as well, but code 1 is the standard code, used by most nuclear genes.")
                 : closer ? "as well, but code " + C + " differs least from code " + T + "."
                          : std::string("these sequences too: choose the code of the organism and organelle the genes come from."));
    }
    // The summary line names the option; the report keeps the app's text (report.py)
    gc_hint_cli_ = text + " Re-run with -codontable " + C + ".";
    text += " Re-run with code " + C + "?";
    gc_hint_ = Json::object();
    gc_hint_.set("code", "genetic_code");
    gc_hint_.set("text", text);
    gc_hint_.set("suggested_codontable", c);
    gc_hint_.set("checked", nt);
    gc_hint_.set("fit", top);
    gc_hint_.set("fit_current", fits_now);
    has_gc_hint_ = true;
}

std::vector<Json> RunReport::id_hints() const {
    if (partial_ || !(pairing_ == "order" || pairing_.empty()) || !has_ids_ || !has_dna_file_ids_) return {};
    std::set<std::string> dset(dna_file_ids_.begin(), dna_file_ids_.end());
    std::vector<std::string> missing;
    for (auto& p : ids_)
        if (!dset.count(p)) missing.push_back(p);
    if (missing.empty()) return {};
    auto c = counts();
    long long tot = (long long)ids_.size();
    bool clean = exit_code == 0 && tot && c["ok"] == tot;
    std::string verdict = pairing_.empty() ? "COAL stopped before pairing them."
                          : clean ? "Sequences were paired by their order in the files instead, and every pair "
                                    "matched codon for codon, so the pairing looks right; renaming makes it safe."
                                  : "Sequences were therefore paired by their order in the files, which is probably "
                                    "wrong for the sequences that failed or have warnings.";
    std::map<std::string, std::string> dnorm;
    for (auto& d : dna_file_ids_) dnorm.emplace(strip_version(d), d);
    std::vector<std::string> fix;
    for (auto& p : missing)
        if (dnorm.count(strip_version(p))) fix.push_back(p);
    Json h = Json::object();
    h.set("code", "id_mismatch");
    if (!fix.empty() && fix.size() == missing.size()) {
        const std::string& ex = fix[0];
        h.set("text", "The protein IDs don't exactly match the DNA IDs: they differ only by a version suffix (for "
                      "example protein '" + ex + "' vs DNA '" + dnorm[strip_version(ex)] + "'). " + verdict +
                          " Rename the sequences so the IDs match exactly, for example by removing the '.1'-style "
                          "version suffixes from the protein alignment ('" + ex + "' -> '" + strip_version(ex) +
                          "'), then run again.");
    } else {
        h.set("text", std::to_string(missing.size()) + " protein ID(s) were not found in the DNA file (for example '" +
                          missing[0] + "'). " + verdict +
                          " Check that both files describe the same genes, in the same order.");
    }
    return {h};
}

std::vector<Json> RunReport::hints() const {
    std::vector<Json> h = pair_hints_;
    for (auto& x : id_hints()) h.push_back(x);
    if (has_gc_hint_) h.push_back(gc_hint_);
    return h;
}

bool RunReport::completed() const {
    return (exit_code == 0 || (exit_code == 2 && keep_going_)) && errors_.empty() && has_output_;
}

std::string RunReport::command() const {
    std::vector<std::string> parts{"coal", pep_};
    parts.insert(parts.end(), nucs_.begin(), nucs_.end());
    if (fmt_ != "fasta") {
        parts.push_back("-output");
        parts.push_back(fmt_);
    }
    if (table_ != 1) {
        parts.push_back("-codontable");
        parts.push_back(std::to_string(table_));
    }
    if (nogap_) parts.push_back("-nogap");
    if (nomismatch_) parts.push_back("-nomismatch");
    if (blockonly_) parts.push_back("-blockonly");
    if (strict_) parts.push_back("--strict");
    if (keep_going_) parts.push_back("--keep-going");
    if (!out_.empty()) {
        parts.push_back("-o");
        parts.push_back(out_);
    }
    return shlex_join(parts);
}

std::vector<std::string> RunReport::summary_lines() const {
    std::map<char, long long> kc;
    for (char ch : kinds_) kc[ch]++;
    const long long n = total();
    std::vector<std::string> lines;
    long long k = nomismatch_ ? 0 : kc['m'];
    if (k) {
        std::string fate = nomismatch_ ? "removed by -nomismatch" : "kept in the alignment";
        std::string s = fmt_n(k) + " of " + fmt_n(n) + " sequences " + (k == 1 ? "has" : "have") +
                        " codons that do not match their amino acid (" + fate + ").";
        if (has_gc_hint_) s += " " + gc_hint_cli_;
        lines.push_back(s);
    }
    k = nomismatch_ ? 0 : kc['u'];
    if (k)
        lines.push_back(fmt_n(k) + " of " + fmt_n(n) + " sequences " + (k == 1 ? "has" : "have") +
                        " codons with unknown or ambiguous bases (e.g. N) that could not be checked: missing "
                        "sequence data, not a genetic-code problem.");
    std::vector<long long> fl, ex;
    for (auto& [i, kp] : entries_) {
        // gene_kind(entry) is the kind recorded in kinds_ (set_kind at add_sequence; the end-of-run step's
        // break_event and mismatch-pattern messages never change a kind)
        char g = i < (long long)kinds_.size() ? kinds_[(size_t)i] : kind_char(gene_kind(json_parse(kp.json)));
        if (g == 'f') fl.push_back(i);
        if (g == 'e') ex.push_back(i);
    }
    if (!fl.empty()) {
        const Json first = entry(fl[0]);
        bool one = fl.size() == 1;
        lines.push_back(seqs((long long)fl.size()) + " of " + fmt_n(n) + " " + (one ? "was" : "were") +
                        " written to the output but " + (one ? "is" : "are") + " probably wrong (first: " +
                        sget(first, "id") + ": " + sget(msgs_of(first)[0], "text") + ") Remove " + (one ? "it" : "them") +
                        " before analysis, or use --strict to refuse such genes.");
    }
    auto c = counts();
    if (!ex.empty()) {
        const Json first = entry(ex[0]);
        const auto& ms = msgs_of(first);
        std::string why = ms.empty() ? "" : sget(ms.back(), "text"), where;
        for (auto& m : ms) {
            std::string cd = sget(m, "code");
            if (cd != "frameshift" && cd != "missing_dna" && !(cd == "inconsistent" && m.get("column"))) {
                why = sget(m, "text");
                break;
            }
        }
        for (auto& m : ms)
            if (m.get("column")) {
                where = sget(m, "text");
                break;
            }
        std::string gc = has_gc_hint_ && (nomismatch_ || !kc['m']) ? " " + gc_hint_cli_ : "";
        const std::string fid = sget(first, "id");
        if (keep_going_) {
            std::string left = errors_.empty() ? std::string(" and ") + (ex.size() == 1 ? "was" : "were") + " left out" : "";
            lines.push_back(seqs((long long)ex.size()) + " of " + fmt_n(n) + " could not be converted" + left + " (first: " +
                            fid + ": " + why + "." + (where.empty() ? "" : " " + where) + ")" + gc);
        } else {
            // Where the gene breaks without repeating the error, then what to do (report.py)
            if (!where.empty()) lines.push_back(fid + ": " + where + gc);
            else if (!gc.empty()) lines.push_back(gc.substr(1));
            long long left = c["not_checked"];
            std::string s = "The run stopped at " + fid;
            s += left ? ": " + seqs(left) + " " + (left == 1 ? "was" : "were") + " not checked and no output was written."
                      : std::string(", so no output was written.");
            if (n > 1) s += " Add --keep-going to write the other sequences and leave this one out.";
            lines.push_back(s);
        }
    }
    // A run-level error is not repeated as "The run failed: ..." (report.py)
    if (pairing_ == "normalised") {
        long long kk = 0, j = -1;
        for (size_t i = 0; i < seq_pairing_.size(); ++i)
            if (seq_pairing_[i] == "normalised") {
                ++kk;
                if (j < 0) j = (long long)i;
            }
        std::string eg = j >= 0 && has_ids_ ? ", e.g. protein " + ids_[(size_t)j] + " = DNA " + dna_id((size_t)j) : "";
        lines.push_back(fmt_n(kk) + " of " + fmt_n(n) + " IDs matched only after normalising (version suffix, case or "
                        "database fields)" + eg + "; COAL paired by those IDs where PAL2NAL's file order would mispair them.");
    }
    long long unk = 0;
    for (auto& [i, kp] : entries_) unk += kp.unknown_aa;
    if (unk)
        lines.push_back(fmt_n(unk) + " of " + fmt_n(n) + " sequences contain non-standard amino acids (U, O, B, Z or J) "
                        "whose codons were not checked.");
    return lines;
}

std::string RunReport::methods_text(const Json& rep) const {
    if (!completed()) return "";
    const Json* s = rep.get("summary");
    long long tot = s && s->get("total") ? s->get("total")->as_int() : 0;
    long long ok = s && s->get("ok") ? s->get("ok")->as_int() : 0;
    std::string T = std::to_string(table_);
    std::vector<std::string> parts{
        "Codon alignments were built with COAL v" + std::string(VERSION) + " (" + REPO_URL +
        "), which places each coding DNA sequence onto the protein multiple sequence alignment codon by codon (the "
        "method of PAL2NAL; Suyama et al. 2006), using the " + code_name(table_) + " genetic code (NCBI translation "
        "table " + T + ")."};
    // gene_kind of every sequence's entry view: the run completed (output written), so a view is its entry, whose
    // kind is the one recorded in kinds_ (report.py counts them over the report's "sequences")
    std::map<std::string, size_t> kinds;
    for (char ch : kinds_) {
        const char* name = ch == 'o' ? "ok" : ch == 'm' ? "mismatch" : ch == 'u' ? "unchecked" : ch == 'x' ? "other"
                           : ch == 'f' ? "flagged" : ch == 'e' ? "excluded" : "not_checked";
        kinds[name] += 1;
    }
    auto were = [](size_t k) { return std::string(k == 1 ? "was" : "were"); };
    std::vector<std::string> clauses;
    if (kinds["mismatch"])
        clauses.push_back(std::to_string(kinds["mismatch"]) + " contained codons that did not match the aligned "
                          "amino acid" + (nomismatch_ ? "" : " (kept in the alignment)"));
    if (kinds["unchecked"])
        clauses.push_back(std::to_string(kinds["unchecked"]) + " contained codons with unknown or ambiguous bases "
                          "(e.g. N) that could not be checked" + (nomismatch_ ? "" : " (kept in the alignment)"));
    if (kinds["other"]) {
        std::set<std::string> codes;
        for (auto& [k, kp] : entries_)
            if (k < (long long)kinds_.size() && kinds_[(size_t)k] == 'x') {
                const Json e = json_parse(kp.json);
                for (auto& m : msgs_of(e)) codes.insert(sget(m, "code"));
            }
        std::vector<std::string> why;
        if (codes.count("unknown_aa"))
            why.push_back("contained non-standard amino acids (such as U, B or Z) whose codons were not checked");
        if (codes.count("suspect_pairing")) why.push_back("were paired with their DNA only by their order in the files");
        bool other = false;
        for (auto& x : codes) other = other || (x != "unknown_aa" && x != "suspect_pairing");
        if (other) why.push_back("had other warnings (see the report)");
        std::string w;
        if (why.size() == 1) w = why[0];
        else {
            w = "had warnings: they ";
            for (size_t i = 0; i < why.size(); ++i) w += (i ? " or " : "") + why[i];
        }
        clauses.push_back(std::to_string(kinds["other"]) + " " + w);
    }
    if (kinds["flagged"]) {
        size_t k = kinds["flagged"];
        clauses.push_back(std::to_string(k) + " " + were(k) + " written to the output but flagged as probably wrong "
                          "(most likely the wrong DNA, a shifted position or the wrong genetic code); remove " +
                          (k == 1 ? "it" : "them") + " before analysis");
    }
    if (kinds["excluded"]) {
        size_t k = kinds["excluded"];
        clauses.push_back(std::to_string(k) + " could not be aligned and " + were(k) + " excluded");
    }
    std::string cnt = "Of " + std::to_string(tot) + " sequences, " + std::to_string(ok) + " matched their protein exactly";
    if (clauses.size() > 1) {
        cnt += "; ";
        for (size_t i = 0; i + 1 < clauses.size(); ++i) cnt += (i ? "; " : "") + clauses[i];
        cnt += "; and " + clauses.back();
    } else if (clauses.size() == 1) {
        cnt += (kinds["flagged"] ? "; " : ", and ") + clauses[0];
    }
    parts.push_back(cnt + ".");
    long long removed = s && s->get("columns_removed") ? s->get("columns_removed")->as_int() : 0;
    if (nomismatch_ && (kinds["mismatch"] || kinds["unchecked"]))
        parts.push_back(removed ? "The " + std::to_string(removed) + " codon columns holding such codons were removed (-nomismatch)."
                                : std::string("The alignment columns holding such codons were removed (-nomismatch)."));
    if (nogap_) parts.push_back("Codon columns containing gaps or stop codons were removed.");
    std::string cmd = sget(rep, "command");
    std::vector<std::string> paths;
    if (const Json* in = rep.get("inputs"))
        for (auto& x : in->a) paths.push_back(sget(x, "path"));
    if (const Json* o = rep.get("output"))
        if (o->type == Json::Object) paths.push_back(sget(*o, "path"));
    for (auto& p : paths) {
        if (p.empty()) continue;
        // As report.py: replace the quoted full path first (Windows paths are always quoted)
        std::string q = shlex_quote(p);
        if (q != p) cmd = replace_all(cmd, q, shlex_quote(basename_of(p)));
        cmd = replace_all(cmd, p, basename_of(p));
    }
    if (!cmd.empty()) parts.push_back("Command: " + cmd);
    std::string out;
    for (size_t i = 0; i < parts.size(); ++i) out += (i ? " " : "") + parts[i];
    if (tot && 2 * ok < tot)  // a note to fix the warnings first (report.py)
        out = "Fix the warnings first: only " + std::to_string(ok) + " of " + std::to_string(tot) +
              " sequences matched their protein exactly, so this alignment is probably not ready to use (see the hints "
              "and the messages for each sequence). Methods text for this run as it stands: " + out;
    return out;
}

void RunReport::write_json(std::ostream& out) const {
    auto c = counts();
    long long n = total();
    Json inputs = Json::array();
    auto hash = [&](size_t i, const std::string& p) {
        // "path" is the argument as given (`-` for standard input); the hash is of the data read (the temp copy)
        bool have = hash_paths_.size() == nucs_.size() + 1;
        if (!have && p == "-") return Json();  // never a file named '-' in the current folder (report.py)
        std::string h = sha256_file(have ? hash_paths_[i] : p);
        return h.empty() ? Json() : Json(h);
    };
    Json pin = Json::object();
    pin.set("role", "protein");
    pin.set("path", pep_);
    pin.set("sha256", hash(0, pep_));
    pin.set("sequences", has_ids_ ? Json((long long)ids_.size()) : Json());
    inputs.push(pin);
    for (size_t i = 0; i < nucs_.size(); ++i) {
        Json x = Json::object();
        x.set("role", "dna");
        x.set("path", nucs_[i]);
        x.set("sha256", hash(i + 1, nucs_[i]));
        x.set("sequences", has_dna_counts_ && i < dna_counts_.size() ? Json(dna_counts_[i]) : Json());
        inputs.push(x);
    }
    std::vector<std::string> p2n{"pal2nal.pl", pep_};
    p2n.insert(p2n.end(), nucs_.begin(), nucs_.end());
    p2n.push_back("-output");
    p2n.push_back(fmt_);
    if (table_ != 1) {
        p2n.push_back("-codontable");
        p2n.push_back(std::to_string(table_));
    }
    if (nogap_) p2n.push_back("-nogap");
    if (nomismatch_) p2n.push_back("-nomismatch");
    if (blockonly_) p2n.push_back("-blockonly");
    Json summary = Json::object();
    summary.set("total", n);
    summary.set("ok", c["ok"]);
    summary.set("warning", c["warning"]);
    summary.set("failed", c["failed"]);
    summary.set("pairing", str_or_null(pairing_, !pairing_.empty()));
    if (c["not_checked"]) summary.set("not_checked", c["not_checked"]);
    if (nomismatch_ && has_output_) {
        summary.set("columns_removed", (long long)cols_.size());
        long long amb = 0;
        for (long long x : cols_) amb += !std::binary_search(real_cols_.begin(), real_cols_.end(), x);
        summary.set("columns_removed_ambiguous", amb);
    }
    Json d = Json::object();
    d.set("schema", SCHEMA);
    d.set("coal_version", VERSION);
    d.set("source", "coal");
    d.set("engine", "coal-cpp");  // the same string as the --resume checkpoint's "engine" (main.cpp)
    d.set("command", command());
    std::vector<std::string> cl{"coal"};
    cl.insert(cl.end(), argv_.begin(), argv_.end());
    d.set("command_line", shlex_join(cl));
    d.set("pal2nal_equivalent", shlex_join(p2n));
    Json opts = Json::object();
    opts.set("output_format", fmt_);
    opts.set("codontable", table_);
    opts.set("nogap", nogap_);
    opts.set("nomismatch", nomismatch_);
    opts.set("keep_going", keep_going_);
    opts.set("strict", strict_);
    opts.set("blockonly", blockonly_);
    d.set("options", opts);
    d.set("inputs", inputs);
    d.set("output", has_output_ ? output_ : Json());
    d.set("summary", summary);
    Json h = Json::array();
    for (auto& x : hints()) h.push(x);
    d.set("hints", h);
    Json errs = Json::array();
    for (auto& x : errors_) errs.push(x);
    d.set("errors", errs);
    // json_dump(d, 1) with the "sequences" array written one entry at a time, in its place
    auto item = [&](const char* key, const Json& v, bool first = false) {
        std::string s = first ? "{\n " : ",\n ";
        json_dump_to(s, Json(key), 1, 1);
        s += ": ";
        json_dump_to(s, v, 1, 1);
        out << s;
    };
    bool first = true;
    for (auto& [key, v] : d.o) {
        item(key.c_str(), v, first);
        first = false;
    }
    out << ",\n \"sequences\": ";
    if (n == 0) {
        out << "[]";
    } else {
        out << '[';
        for (long long k = 0; k < n; ++k) {
            std::string s = k ? ",\n  " : "\n  ";
            json_dump_to(s, entry_view(k), 1, 2);
            out << s;
        }
        out << "\n ]";
    }
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_).count();
    item("runtime_seconds", Json((double)(long long)(secs * 1000.0 + 0.5) / 1000.0));
    item("exit_code", Json(exit_code));
    Json st = Json::array();
    for (auto& l : summary_lines()) st.push(l);
    item("summary_text", st);
    item("methods_text", Json(methods_text(d)));
    if (has_error_) item("error", error_);
    if (resumed_from_ >= 0) item("resumed_from", Json(resumed_from_));
    out << "\n}";
}

}  // namespace coal
