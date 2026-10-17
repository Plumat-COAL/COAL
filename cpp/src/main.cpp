// coal-cpp: command-line front end and streaming, multithreaded pipeline.
// Mirrors src/coal/cli.py + src/coal/pipeline.py: same options, same output,
// same warnings, same exit codes.
//
// Memory does not grow with the input beyond the per-sequence index (IDs and
// byte offsets): pass 1 reads protein IDs and the DNA index; pass 2 reads one
// batch (~2M nt) at a time, threads it on N worker threads while the next batch
// is read, and writes results in input order. Output goes to a temporary file
// that replaces the destination only when everything succeeded.
#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdlib>
#include <functional>
#include <chrono>
#include <climits>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <unordered_map>
#include <string>
#include <thread>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include "args.hpp"
#include "core.hpp"
#include "help_text.hpp"
#include "pool.hpp"
#include "json.hpp"
#include "report.hpp"
#include "seqio.hpp"
#include "sysio.hpp"
#include "tables.hpp"

#ifdef _WIN32
#define COAL_ISATTY(fd) _isatty(fd)
#else
#include <fcntl.h>
#include <unistd.h>
#define COAL_ISATTY(fd) isatty(fd)
#endif

namespace fs = std::filesystem;
using namespace coal;

static const char* VERSION = "0.4.0";
static long long env_ll(const char* name, long long dflt) {
    const char* v = std::getenv(name);
    if (!v || !*v) return dflt;
    try {
        return std::stoll(v);
    } catch (...) {
        return dflt;
    }
}
static double env_double(const char* name, double dflt) {
    const char* v = std::getenv(name);
    if (!v || !*v) return dflt;
    try {
        return std::stod(v);
    } catch (...) {
        return dflt;
    }
}
// Same batch boundaries as pipeline.py (COAL_BATCH_NT overrides, for tests/tuning).
static const long long BATCH_NT = env_ll("COAL_BATCH_NT", 2000000);
static const int WRAP = 60;
static const char* const MISMATCH_WARNING = " does not correspond to ";
static const long long BLOCK_CACHE = 64LL << 20;  // bytes of rows held while writing blocks

struct Options {
    std::string pep;
    std::vector<std::string> nucs;
    int table = 1;
    bool nogap = false, nomismatch = false, quiet = false, blockonly = false, strict = false;
    std::string fmt = "fasta";
    std::string out;
    int jobs = 0;
    bool keep_going = false, progress = false, resume = false, json = false;
    bool resumable = true;  // CPP-28 (reference MRG-4): false when an input came through a pipe (spooled)
    int max_warnings = 20;  // CHANGED(IMP-E-C): was 0 (no limit); cli.py MAX_WARNINGS
    std::string report;
    std::vector<std::string> argv;  // raw arguments (for the report's "command")
};

// ------------------------------------------------------------------ temp files

static std::string unique_suffix() {
    static std::atomic<unsigned> counter{0};
    std::random_device rd;
    return std::to_string(rd()) + "_" + std::to_string(counter++);
}


// ------------------------------------------------------------------ formatting

static void wrap_into(std::string& out, const std::string& seq) {
    if (seq.size() <= (size_t)WRAP) {
        out += seq;
        out += '\n';
        return;
    }
    for (size_t i = 0; i < seq.size(); i += WRAP) {
        out.append(seq, i, std::min((size_t)WRAP, seq.size() - i));
        out += '\n';
    }
}

static void record_into(std::string& out, const std::string& fmt, const std::string& id, const std::string& seq) {
    if (fmt == "fasta") out += '>';
    out += id;
    out += '\n';
    wrap_into(out, seq);
}

static std::string paml_header(long long n, long long len) {
    char b[64];
    std::snprintf(b, sizeof b, " %3lld %6lld\n", n, len);  // CPP-20: as PAL2NAL (was %4lld%7lld: 1,000+ sequences)
    return b;
}

static std::string apply_slices(const std::string& s, const std::vector<std::pair<long long, long long>>& sl) {
    std::string out;
    for (auto& [a, b] : sl) out.append(s, (size_t)a, (size_t)(b - a));
    return out;
}

// ------------------------------------------------------------------ ID pairing
// Same rules as pair_ids() in src/coal/core.py.

static std::string strip_version(const std::string& x) {
    size_t dot = x.rfind('.');
    if (dot == std::string::npos || dot == 0 || dot + 1 >= x.size()) return x;
    for (size_t i = dot + 1; i < x.size(); ++i)
        if (x[i] < '0' || x[i] > '9') return x;
    return x.substr(0, dot);
}

static std::string lower(std::string s) {
    for (char& c : s)
        if (c >= 'A' && c <= 'Z') c = char(c + 32);
    return s;
}

static std::vector<std::string> id_keys(const std::string& rid) {
    static const std::set<std::string> DB = {"sp", "tr", "gi", "ref", "gb", "emb", "dbj", "pdb", "lcl", "gnl",
                                             "pir", "prf", "tpg", "tpe", "tpd", "bbs", "bbm", "pat"};
    std::vector<std::string> keys{lower(strip_version(rid))};
    std::vector<std::string> fields;
    size_t a = 0;
    while (a <= rid.size()) {
        size_t b = rid.find('|', a);
        if (b == std::string::npos) b = rid.size();
        if (b > a) fields.push_back(rid.substr(a, b - a));
        a = b + 1;
    }
    if (fields.size() > 1)
        for (auto& f : fields)
            if (!DB.count(lower(f))) keys.push_back(lower(strip_version(f)));
    return keys;
}

static std::string records(long long n, const std::string& what) {
    return std::to_string(n) + " " + what + (n == 1 ? "" : "s");
}

// CHANGED(ENG-UX-6): the count error with what to do; an example record when the IDs show which records are extra
// (every protein ID is a DNA ID) or which proteins have no DNA (exactly that many are missing). As core.py.
static std::string count_mismatch_message(const std::vector<std::string>& prot, const std::vector<std::string>& dna_ids,
                                          const std::unordered_map<std::string, size_t>& pos) {
    const long long p = (long long)prot.size(), d = (long long)dna_ids.size();
    std::string head = "number of sequences differs (protein: " + std::to_string(p) + ", DNA: " + std::to_string(d) +
                       "): each protein needs exactly one DNA record";
    if (d > p) {
        std::string eg;
        bool all = true;
        for (auto& x : prot) all = all && pos.count(x);
        if (all) {
            std::set<std::string> pset(prot.begin(), prot.end());
            for (auto& x : dna_ids)
                if (!pset.count(x)) {
                    eg = " (e.g. '" + x + "', which is not in the protein alignment)";
                    break;
                }
        }
        return head + "; remove the " + records(d - p, "extra DNA record") + eg + " or check you used the right files";
    }
    std::vector<std::string> missing;
    for (auto& x : prot)
        if (!pos.count(x)) missing.push_back(x);
    std::string eg = !missing.empty() && (long long)missing.size() == p - d ? " (e.g. '" + missing[0] + "')" : "";
    return head + "; " + std::to_string(p - d) + " protein" + (p - d == 1 ? " has" : "s have") + " no DNA record" + eg +
           ": add the missing DNA, remove those proteins from the alignment, or check you used the right files";
}

// CPP-25 (reference MRG-1, core.py pair_ids): returns "id", "normalised" or "order"; fills idx, per_seq (only for
// "normalised") and the notices.
static std::string pair_ids(const std::vector<std::string>& prot, const std::vector<std::string>& dna_ids,
                            const std::unordered_map<std::string, size_t>& pos, std::vector<size_t>& idx,
                            std::vector<std::string>& per_seq, std::vector<std::string>& notes, bool strict) {
    if (prot.size() != dna_ids.size())  // CHANGED(ENG-UX-6): + what to do (core.py count_mismatch_message)
        throw CoalError(count_mismatch_message(prot, dna_ids, pos), "count_mismatch");
    if (prot.empty()) throw CoalError("no sequences in input", "no_sequences");
    idx.assign(prot.size(), 0);
    per_seq.clear();
    long long hits = 0;
    for (size_t k = 0; k < prot.size(); ++k) {
        auto it = pos.find(prot[k]);
        if (it != pos.end()) {
            idx[k] = it->second;
            ++hits;
        }
    }
    if (hits == (long long)prot.size()) return "id";
    auto first_diff = [&](std::string& p, std::string& d) {
        for (size_t k = 0; k < prot.size(); ++k)
            if (prot[k] != dna_ids[k]) {
                p = prot[k];
                d = dna_ids[k];
                return;
            }
    };
    const std::string H = std::to_string(hits), T = std::to_string(prot.size());
    if (hits && strict) {
        std::string p, d;
        first_diff(p, d);
        throw CoalError("only " + H + " of " + T + " protein IDs are DNA IDs, so sequences would be paired by order "
                            "(e.g. protein '" + p + "' with DNA '" + d + "'); make the IDs match, or run without "
                            "--strict to pair by order as PAL2NAL does",
                        "partial_id_match");
    }
    // normalised_pairing: one candidate per protein (an exact ID wins), one-to-one
    std::unordered_map<std::string, std::vector<size_t>> by_key;
    for (size_t i = 0; i < dna_ids.size(); ++i)
        for (auto& key : id_keys(dna_ids[i])) {
            auto& v = by_key[key];
            if (v.empty() || v.back() != i) v.push_back(i);
        }
    std::vector<size_t> nidx(prot.size());
    std::vector<std::string> modes;
    std::vector<char> used(dna_ids.size(), 0);
    bool ok = true;
    for (size_t k = 0; k < prot.size() && ok; ++k) {
        size_t j;
        auto it = pos.find(prot[k]);
        if (it != pos.end()) {
            j = it->second;
            modes.push_back("id");
        } else {
            std::set<size_t> cands;
            for (auto& key : id_keys(prot[k])) {
                auto f = by_key.find(key);
                if (f != by_key.end()) cands.insert(f->second.begin(), f->second.end());
            }
            if (cands.size() != 1) {
                ok = false;
                break;
            }
            j = *cands.begin();
            modes.push_back("normalised");
        }
        if (used[j]) {
            ok = false;
            break;
        }
        used[j] = 1;
        nidx[k] = j;
    }
    bool moved_any = false;
    for (size_t k = 0; ok && k < prot.size(); ++k) moved_any = moved_any || nidx[k] != k;
    if (ok && moved_any) {
        long long nn = 0, nmoved = 0;
        size_t first = prot.size();
        for (size_t k = 0; k < prot.size(); ++k) {
            nn += modes[k] == "normalised";
            if (nidx[k] != k) {
                ++nmoved;
                if (first == prot.size()) first = k;
            }
        }
        notes.push_back(std::to_string(nn) + " of " + T + " protein IDs match a DNA ID only after normalising (version "
                        "suffix, case or database fields), so sequences were paired by those IDs, not by file order "
                        "as PAL2NAL pairs them: by order, " + std::to_string(nmoved) + " of " + T +
                        " proteins would get the DNA of another gene (e.g. by file order, protein '" + prot[first] +
                        "' would get DNA '" + dna_ids[first] + "'; COAL paired it with '" + dna_ids[nidx[first]] +
                        "')");  // CHANGED(ENG-UX-1): the example read backwards
        idx = nidx;
        per_seq = modes;
        return "normalised";
    }
    if (hits) {
        std::string p, d;
        first_diff(p, d);
        notes.push_back("only " + H + " of " + T + " protein IDs are DNA IDs; pairing sequences by order instead (e.g. "
                        "protein '" + p + "' gets DNA '" + d + "') (--strict makes this an error)");
    }
    for (size_t k = 0; k < prot.size(); ++k) idx[k] = k;
    return "order";
}

// ------------------------------------------------------------------ pipeline

struct Item {
    std::string id, row, dna;  // input (row prepared)
    std::string out;           // finished record, or the expanded row when staged
    std::vector<std::string> warnings;
    std::vector<int> mcols;
    std::vector<std::string> notices;  // CPP-15
    bool failed = false;
    std::string error, error_code;
};

struct Ctx {
    const GeneticCode* code;
    std::vector<int> widths;
    bool staged;
    std::string fmt;
    bool strict;  // CPP-15
};

static void work(Item& it, const Ctx& c) {
    try {
        Threaded t = thread_one(it.id, it.row, it.dna, *c.code, c.strict);
        std::string exp = expand(it.row, t.dna, c.widths);
        it.warnings = std::move(t.warnings);
        it.mcols = std::move(t.mismatch_cols);
        it.notices = std::move(t.notices);
        if (c.staged) {
            it.out = std::move(exp);
        } else {
            it.out.clear();
            record_into(it.out, c.fmt, it.id, exp);
        }
    } catch (const CoalError& e) {
        it.failed = true;
        it.error = e.what();
        it.error_code = e.code;
    }
    std::string().swap(it.dna);
}

static void read_next(std::istream& in, long long count, std::string& dst) {
    dst.resize((size_t)count);
    if (count) in.read(&dst[0], (std::streamsize)count);
}

static void read_rows(std::istream& in, long long offset, long long count, std::string& dst) {
    dst.resize((size_t)count);
    if (count == 0) return;
    in.clear();
    in.seekg(offset);
    in.read(&dst[0], (std::streamsize)count);
}

// CHANGED(ENG-UX-20): the residue-line widths of -output codon come from the last protein row of the input, also
// when --keep-going left that sequence out (as pipeline.py codon_aa_widths(last_row)); was the last row written, so
// the C++ and Python outputs differed (found by TestCppMatchesPython, seed 31337).
static std::string g_last_row;

// CPP-24 (reference FIN-E2: bounded memory): clustal / codon output for many sequences. Every line of a block has a
// fixed length, so each line's place in the output is known in advance: the rows are read a chunk of consecutive
// sequences at a time (sequential reads) and each block's lines for that chunk are written at their offset. Memory
// stays bounded by BLOCK_CHUNK_BYTES, not by the number of sequences, and no block of every sequence is held.
static const long long BLOCK_SEQ_LIMIT = env_ll("COAL_BLOCK_MEMORY", 16000000);  // as the reference's BLOCK_MEMORY
static const long long BLOCK_CHUNK_BYTES = env_ll("COAL_BLOCK_CHUNK", 32LL << 20);  // tests: small chunks

static void write_blocks_chunked(std::ostream& out, const std::string& fmt, std::fstream& grid, long long length,
                                 std::fstream* prot, long long ncol, const std::vector<std::string>& ids,
                                 const std::string& mask) {
    const bool codon = fmt == "codon";
    const long long n = (long long)ids.size();
    size_t w0 = 10;
    for (auto& i : ids) w0 = std::max(w0, i.size());
    const long long W = (long long)w0 + 4;
    // residue-line tokens (codon): as in write_blocks
    std::vector<long long> toks;
    if (codon && n > 0) {
        const std::string& last = g_last_row;  // ENG-UX-20: the last input row (pipeline.py last_row)
        for (long long i = 0; i < (long long)last.size(); ++i) {
            toks.push_back(i);
            char ch = last[(size_t)i];
            if (ch >= '0' && ch <= '9')
                for (int e = 0; e < std::max(ch - '0' - 1, 0) / 3; ++e) toks.push_back(-1);
        }
    }
    const long long ntok = (long long)toks.size();
    const long long nb = codon ? std::max((ntok + 19) / 20, (length / 3 + 19) / 20) : (length + 59) / 60;
    auto trip_chars = [](long long len) { return len > 0 ? len + (len + 2) / 3 - 1 : 0; };  // "AAA CCC G"
    const size_t nbz = (size_t)nb;
    std::vector<long long> seg0(nbz), seg1(nbz), ta(nbz), te(nbz), per(nbz), off(nbz);
    long long pos = fmt == "clustal" ? 39 : 0;  // "CLUSTAL W multiple sequence alignment\n\n"
    for (long long b = 0; b < nb; ++b) {
        seg0[(size_t)b] = std::min(length, b * 60);
        seg1[(size_t)b] = std::min(length, (b + 1) * 60);
        long long sl = seg1[(size_t)b] - seg0[(size_t)b];
        ta[(size_t)b] = std::min(ntok, b * 20);
        te[(size_t)b] = std::min(ntok, (b + 1) * 20);
        long long m = te[(size_t)b] - ta[(size_t)b];
        per[(size_t)b] = codon ? (W + (m ? 4 * m - 3 : 0) + 1) + (W + trip_chars(sl) + 1) : W + sl + 1;
        off[(size_t)b] = pos;
        long long mline = mask.empty() ? 0 : W + (codon ? trip_chars(sl) : sl) + 1;
        pos += n * per[(size_t)b] + mline + 1;
    }
    if (fmt == "clustal") {
        out.seekp(0);
        out << "CLUSTAL W multiple sequence alignment\n\n";
    }
    const long long chunk = std::max(1LL, BLOCK_CHUNK_BYTES / std::max(1LL, length + ncol + 2 * W + 64));
    std::string rows, prows, text;
    const std::string pad((size_t)W, ' ');
    for (long long k0 = 0; k0 < n; k0 += chunk) {
        long long k1 = std::min(n, k0 + chunk);
        read_rows(grid, k0 * length, (k1 - k0) * length, rows);
        if (codon) read_rows(*prot, k0 * ncol, (k1 - k0) * ncol, prows);
        for (long long b = 0; b < nb; ++b) {
            text.clear();
            text.reserve((size_t)((k1 - k0) * per[(size_t)b]));
            for (long long k = k0; k < k1; ++k) {
                const char* row = rows.data() + (size_t)((k - k0) * length);
                if (codon) {
                    const char* prow = prows.data() + (size_t)((k - k0) * ncol);
                    text += pad;
                    for (long long t = ta[(size_t)b]; t < te[(size_t)b]; ++t) {
                        if (t > ta[(size_t)b]) text += "   ";
                        long long c = toks[(size_t)t];
                        text += c >= 0 ? prow[(size_t)c] : '-';
                    }
                    text += '\n';
                    text += ids[(size_t)k];
                    text.append((size_t)(W - (long long)ids[(size_t)k].size()), ' ');
                    for (long long i = seg0[(size_t)b]; i < seg1[(size_t)b]; i += 3) {
                        if (i > seg0[(size_t)b]) text += ' ';
                        text.append(row + (size_t)i, (size_t)std::min(3LL, seg1[(size_t)b] - i));
                    }
                    text += '\n';
                } else {
                    text += ids[(size_t)k];
                    text.append((size_t)(W - (long long)ids[(size_t)k].size()), ' ');
                    text.append(row + (size_t)seg0[(size_t)b], (size_t)(seg1[(size_t)b] - seg0[(size_t)b]));
                    text += '\n';
                }
            }
            out.seekp(off[(size_t)b] + k0 * per[(size_t)b]);
            out.write(text.data(), (std::streamsize)text.size());
        }
    }
    for (long long b = 0; b < nb; ++b) {  // mark line and blank line closing each block
        text.clear();
        if (!mask.empty()) {
            text += pad;
            if (codon) {
                for (long long i = seg0[(size_t)b]; i < seg1[(size_t)b]; i += 3) {
                    if (i > seg0[(size_t)b]) text += ' ';
                    text.append(mask, (size_t)i, (size_t)std::min(3LL, seg1[(size_t)b] - i));
                }
            } else {
                text.append(mask, (size_t)seg0[(size_t)b], (size_t)(seg1[(size_t)b] - seg0[(size_t)b]));
            }
            text += '\n';
        }
        text += '\n';
        out.seekp(off[(size_t)b] + n * per[(size_t)b]);
        out.write(text.data(), (std::streamsize)text.size());
    }
    out.seekp(pos);
}

// CLUSTAL / codon output from fixed-stride row files, a few blocks at a time.
// mask: the '#' mark line (length characters) printed under each block, or empty (CPP-13).
static void write_blocks(std::ostream& out, const std::string& fmt, std::fstream& grid, long long length,
                         std::fstream* prot, long long ncol, const std::vector<std::string>& ids,
                         const std::string& mask) {
    long long n = (long long)ids.size();
    // CPP-24: the chunked writer unless one row is very long (then this one, a few blocks of every row at a time)
    if (length > 0 && (n * (60 + (fmt == "codon" ? 20 : 0)) > BLOCK_SEQ_LIMIT || (length + ncol) * 4 <= BLOCK_CHUNK_BYTES)) {
        write_blocks_chunked(out, fmt, grid, length, prot, ncol, ids, mask);
        return;
    }
    size_t width = 10;
    for (auto& i : ids) width = std::max(width, i.size());
    width += 4;
    std::vector<std::string> names;
    names.reserve(ids.size());
    for (auto& i : ids) names.push_back(i + std::string(width - i.size(), ' '));
    std::string chunk, pchunk, text;

    if (fmt == "clustal") {
        out << "CLUSTAL W multiple sequence alignment\n\n";
        long long nb = (length + 59) / 60;
        long long per = std::max(1LL, BLOCK_CACHE / std::max(1LL, n * 60));
        std::vector<std::string> rows((size_t)n);
        for (long long sb = 0; sb < nb; sb += per) {
            long long s0 = sb * 60, s1 = std::min(length, (sb + per) * 60);
            if (s0 == 0 && s1 == length) {  // whole rows: read sequentially
                grid.seekg(0);
                for (long long k = 0; k < n; ++k) read_next(grid, length, rows[(size_t)k]);
            } else {
                for (long long k = 0; k < n; ++k) read_rows(grid, k * length + s0, s1 - s0, rows[(size_t)k]);
            }
            for (long long b = sb; b < std::min(nb, sb + per); ++b) {
                long long a = b * 60 - s0, e = std::min(s1 - s0, a + 60);
                text.clear();
                for (long long k = 0; k < n; ++k) {
                    text += names[(size_t)k];
                    text.append(rows[(size_t)k], (size_t)a, (size_t)(e - a));
                    text += '\n';
                }
                if (!mask.empty()) {
                    text.append(width, ' ');
                    text.append(mask, (size_t)(b * 60), 60);
                    text += '\n';
                }
                text += '\n';
                out << text;
            }
        }
        return;
    }
    // codon: residue line and codon line cut into blocks of 20 separately (as PAL2NAL).
    // Residue-line tokens: one per column, plus extra '-' tokens after a column where the
    // LAST sequence has a frameshift digit (see residue_line_blocks in pipeline.py).
    long long ntrip = length / 3;
    const std::string& last = g_last_row;  // ENG-UX-20
    std::vector<long long> toks;  // column index, or -1 for an extra '-'
    for (long long i = 0; i < (long long)last.size(); ++i) {
        toks.push_back(i);
        char ch = last[(size_t)i];
        if (ch >= '0' && ch <= '9') {
            int extra = std::max(ch - '0' - 1, 0) / 3;
            for (int e = 0; e < extra; ++e) toks.push_back(-1);
        }
    }
    long long nres = ((long long)toks.size() + 19) / 20;
    long long nb = std::max(nres, (ntrip + 19) / 20);
    long long per = std::max(1LL, BLOCK_CACHE / std::max(1LL, n * 80));
    std::vector<std::string> rows((size_t)n), prows((size_t)n);
    std::string pad(width, ' ');
    for (long long sb = 0; sb < nb; sb += per) {
        long long sb_end = std::min(nb, sb + per);
        long long d0 = std::min(length, sb * 60), d1 = std::min(length, sb_end * 60);
        // columns needed by the residue tokens of these blocks
        long long p0 = ncol, p1 = 0;
        for (long long t = std::min((long long)toks.size(), sb * 20); t < std::min((long long)toks.size(), sb_end * 20); ++t)
            if (toks[(size_t)t] >= 0) {
                p0 = std::min(p0, toks[(size_t)t]);
                p1 = std::max(p1, toks[(size_t)t] + 1);
            }
        if (p1 < p0) p0 = p1 = 0;
        bool whole = d0 == 0 && d1 == length && p0 == 0 && p1 == ncol;
        if (whole) {  // whole rows: read sequentially
            grid.seekg(0);
            prot->seekg(0);
        }
        for (long long k = 0; k < n; ++k) {
            if (whole) {
                read_next(grid, length, rows[(size_t)k]);
                read_next(*prot, ncol, prows[(size_t)k]);
            } else {
                read_rows(grid, k * length + d0, d1 - d0, rows[(size_t)k]);
                read_rows(*prot, k * ncol + p0, p1 - p0, prows[(size_t)k]);
            }
        }
        for (long long b = sb; b < sb_end; ++b) {
            long long da = std::min(length, b * 60) - d0, de = std::min(length, (b + 1) * 60) - d0;
            long long ta = std::min((long long)toks.size(), b * 20), te = std::min((long long)toks.size(), (b + 1) * 20);
            text.clear();
            for (long long k = 0; k < n; ++k) {
                text += pad;
                for (long long t = ta; t < te; ++t) {
                    if (t > ta) text += "   ";
                    long long c = toks[(size_t)t];
                    text += c >= 0 ? prows[(size_t)k][(size_t)(c - p0)] : '-';
                }
                text += '\n';
                text += names[(size_t)k];
                for (long long i = da; i < de; i += 3) {
                    if (i > da) text += ' ';
                    text.append(rows[(size_t)k], (size_t)i, (size_t)std::min(3LL, de - i));
                }
                text += '\n';
            }
            if (!mask.empty()) {
                text += pad;
                long long m0 = std::min(length, b * 60), m1 = std::min(length, (b + 1) * 60);
                for (long long i = m0; i < m1; i += 3) {
                    if (i > m0) text += ' ';
                    text.append(mask, (size_t)i, (size_t)std::min(3LL, m1 - i));
                }
                text += '\n';
            }
            text += '\n';
            out << text;
        }
    }
}

// ------------------------------------------------------------------ durability + checkpoints

static void fsync_path(const std::string& path) {
#ifdef _WIN32
    int fd = _open(path.c_str(), _O_WRONLY | _O_APPEND | _O_BINARY);
    if (fd >= 0) {
        _commit(fd);
        _close(fd);
    }
#else
    int fd = ::open(path.c_str(), O_RDONLY);
    if (fd >= 0) {
        ::fsync(fd);
        ::close(fd);
    }
#endif
}

// Flush to disk and return the size in bytes.
static long long durable_size(std::fstream& f, const std::string& path) {
    f.flush();
    fsync_path(path);
    std::error_code ec;
    auto n = fs::file_size(path, ec);
    return ec ? 0 : (long long)n;
}

static const char* CHANGED =
    "The input files or options changed since the interrupted run, so it can't be resumed safely. "
    "Start a new run without --resume.";

// Same design as src/coal/checkpoint.py (engine "coal-cpp"; each engine resumes only its own runs).
struct Checkpoint {
    std::string partial, stage, prot, entries, path, options_hash;
    std::vector<std::string> inputs;
    int count = 0;
    std::string why = "no checkpoint";  // why load_for_resume() returned false

    Checkpoint(const std::string& out, const std::string& canon, const std::vector<std::string>& in) : inputs(in) {
        std::string base = fs::absolute(fs::path(out)).string();
        partial = base + ".coal-partial";
        stage = base + ".coal-stage";
        prot = base + ".coal-prot";
        entries = base + ".coal-entries.jsonl";
        path = base + ".coal-checkpoint.json";
        std::string c = canon + "|out=" + base;
        for (auto& p : inputs) c += "|in=" + fs::absolute(fs::path(p)).string();
        options_hash = sha256_hex(c);
    }

    Json signature() const {
        Json a = Json::array();
        for (auto& p : inputs) {
            Json x = Json::object();
            std::error_code ec;
            x.set("path", fs::absolute(fs::path(p)).string());
            x.set("size", (long long)fs::file_size(p, ec));
            x.set("mtime", (long long)fs::last_write_time(p, ec).time_since_epoch().count());
            a.push(x);
        }
        return a;
    }

    // true + state if this run can resume; false if there is nothing to resume.
    // Throws CoalError(resume_mismatch) if inputs or options changed.
    bool load_for_resume(Json& st, const std::vector<std::string>& needed) {
        std::ifstream in(path, std::ios::binary);
        if (!in) {
            why = "no checkpoint";
            return false;
        }
        std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        try {
            st = json_parse(text);
        } catch (...) {
            why = "an unreadable checkpoint";
            return false;
        }
        if (st.type != Json::Object) {
            why = "an unreadable checkpoint";
            return false;
        }
        auto str = [&](const char* k) {
            const Json* v = st.get(k);
            return v ? v->as_str() : std::string();
        };
        const Json* ins = st.get("inputs");
        if (str("engine") != "coal-cpp" || str("coal_version") != VERSION || str("options_hash") != options_hash ||
            !ins || json_dump(*ins) != json_dump(signature()))
            throw CoalError(CHANGED, "resume_mismatch");
        why = "a checkpoint whose partial output is missing or incomplete";
        static const char* keys[] = {"sink_bytes", "entries_bytes", "prot_bytes"};
        for (size_t i = 0; i < needed.size(); ++i) {
            if (needed[i].empty()) continue;
            std::error_code ec;
            auto n = fs::file_size(needed[i], ec);
            const Json* v = st.get(keys[i]);
            if (ec || !v || (long long)n < v->as_int()) return false;
        }
        return true;
    }

    void save(Json st) {
        ++count;
        st.set("engine", "coal-cpp");
        st.set("coal_version", VERSION);
        st.set("options_hash", options_hash);
        st.set("inputs", signature());
        std::string tmp = path + ".tmp";
        {
            std::ofstream f(tmp, std::ios::binary | std::ios::trunc);
            f << json_dump(st);
            f.flush();
        }
        fsync_path(tmp);
        std::error_code ec;
        fs::rename(tmp, path, ec);
    }

    void remove(bool keep_partial = false) {
        std::error_code ec;
        for (auto* p : {&path, &entries, &stage, &prot}) fs::remove(*p, ec);
        fs::remove(stage + ".grid", ec);
        fs::remove(path + ".tmp", ec);
        if (!keep_partial) fs::remove(partial, ec);
    }
};

// CHANGED(IMP-E-I): test hooks act only with COAL_TEST_HOOKS=1 as well (pipeline.py _test_hooks)
static bool test_hooks() {
    const char* v = std::getenv("COAL_TEST_HOOKS");
    return v && std::string(v) == "1";
}

static void test_pause(int count) {
    if (!test_hooks()) return;
    const char* v = std::getenv("COAL_TEST_PAUSE_AFTER_CHECKPOINT");
    if (v && std::to_string(count) == v) {
        std::cerr << "coal: test: paused after checkpoint " << count << std::endl;
        std::this_thread::sleep_for(std::chrono::seconds(3600));
    }
}

// The -o argument ("" = stdout), for output_error() when a file of the run cannot be created (ENG-UX-12).
static std::string g_out_arg;

// A file written by the pipeline: fresh, or (resume) truncated to a byte length and appended to.
struct OutFile {
    std::string path;
    std::vector<char> buf;  // declared before f: f's destructor flushes into it
    std::fstream f;
    OutFile() : buf(1 << 20) {}
    void open(const std::string& p, long long resume_bytes = -1) {
        path = p;
        f.rdbuf()->pubsetbuf(buf.data(), (std::streamsize)buf.size());
        errno = 0;
        if (resume_bytes >= 0) {
            fs::resize_file(p, (uintmax_t)resume_bytes);
            f.open(p, std::ios::in | std::ios::out | std::ios::binary);
            f.seekp(0, std::ios::end);
        } else {
            f.open(p, std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);
        }
        // CHANGED(ENG-UX-12): was "[Errno 2] No such file or directory: '<COAL's hidden file>'" whatever the reason;
        // now the user's -o in plain words, as pipeline.py _output_error
        if (!f) throw CoalError(output_error(g_out_arg, errno ? errno : ENOENT), "file_error");
    }
};

// ------------------------------------------------------------------ progress

static const long long PROGRESS_MIN_SEQS = env_ll("COAL_PROGRESS_MIN_SEQS", 10000);
static const double PROGRESS_EVERY = env_double("COAL_PROGRESS_EVERY", 2.0);
static const double CHECKPOINT_EVERY = env_double("COAL_CHECKPOINT_EVERY", 2.0);

struct Progress {
    using clock = std::chrono::steady_clock;
    long long total;
    clock::time_point t0, last;
    bool tty, printed = false;
    explicit Progress(long long t) : total(t), t0(clock::now()), last(t0), tty(COAL_ISATTY(2) != 0) {}
    void update(long long done, bool final = false) {
        auto now = clock::now();
        if (!final && std::chrono::duration<double>(now - last).count() < PROGRESS_EVERY) return;
        last = now;
        long long pct = total ? 100 * done / total : 100;
        double el = std::chrono::duration<double>(now - t0).count();
        char b[64];
        if (final) {
            std::snprintf(b, sizeof b, "done in %.1f s", el);
        } else if (done) {
            double left = el / (double)done * (double)(total - done);
            if (left < 120) std::snprintf(b, sizeof b, "about %.0f s left", left);
            else std::snprintf(b, sizeof b, "about %.0f min left", left / 60);
        } else {
            std::snprintf(b, sizeof b, "estimating time left");
        }
        std::string msg = "coal: progress: " + std::to_string(pct) + "% (" + fmt_n(done) + " of " + fmt_n(total) +
                          " sequences), " + b;
        if (tty) {
            if (msg.size() < 78) msg += std::string(78 - msg.size(), ' ');
            std::cerr << '\r' << msg << (final ? "\n" : "");
        } else {
            std::cerr << msg << '\n';
        }
        std::cerr.flush();
        printed = true;
    }
};

// ------------------------------------------------------------------ run

struct RunState {
    std::vector<size_t> pairing;
    std::vector<std::string> dna_ids_paired;
};

using WarnFn = std::function<void(const std::string&)>;

static long long run(const Options& o, RunReport& rec, RunState& rs, const WarnFn& warn, const WarnFn& note) {
    const GeneticCode* code = get_code(o.table);
    if (o.fmt == "codon" && (o.nogap || o.nomismatch || o.blockonly))
        throw CoalError("-output codon cannot be combined with -nogap, -nomismatch or -blockonly (as in PAL2NAL)");
    const bool block = o.fmt == "clustal" || o.fmt == "codon";
    // CPP-13 (reference DEC-4): PAL2NAL's '#' block marks of a CLUSTAL input: -blockonly keeps only those
    // columns; clustal / codon output show the marks under each block when -blockonly is not given
    const std::string mask = is_clustal(o.pep) ? clustal_mask(o.pep) : std::string();
    const bool has_marks = mask.find('#') != std::string::npos;
    const bool block_cols = o.blockonly && has_marks;
    const bool two_pass = o.nogap || o.nomismatch || block_cols;
    // keep_going + paml: the header needs the final count, so stage the rows
    const bool staged = two_pass || block || (o.keep_going && o.fmt == "paml");

    // pass 1: protein IDs, row length, frameshift column widths
    std::vector<std::string> ids_read;
    long long ncol = -1;
    std::vector<int> widths;
    {
        ProteinReader pr(o.pep);
        std::string id, row, first_row;
        while (pr.next(id, row)) {
            row = prepare_row(row);
            if (ncol < 0) {
                ncol = (long long)row.size();
                widths.assign((size_t)ncol, 3);
                first_row = row;
            } else if ((long long)row.size() != ncol) {
                // CHANGED(ENG-UX-3): files swapped (the first file is DNA) / protein not aligned told apart (core.py)
                throw CoalError(row_length_message(o.pep, ncol, (long long)row.size(), id, first_row + row),
                                "row_lengths");
            }
            update_column_widths(widths, row);
            ids_read.push_back(id);
        }
    }
    if (ids_read.empty())  // CHANGED(ENG-UX-6): was "number of sequences differs (protein: 0, DNA: N)"
        throw CoalError(o.pep + ": no sequences found: the file is empty or has no '>' header lines (expected a "
                                "protein alignment in FASTA or CLUSTAL format)",
                        "no_sequences");
    rec.set_inputs(std::move(ids_read));
    const std::vector<std::string>& ids = rec.ids();  // CPP-33: one copy, owned by the report
    DnaIndex dna(o.nucs);
    if (dna.size() == 0) {  // CHANGED(ENG-UX-6): was "number of sequences differs (protein: N, DNA: 0)"
        std::string joined;
        for (size_t i = 0; i < o.nucs.size(); ++i) joined += (i ? ", " : "") + o.nucs[i];
        const bool many = o.nucs.size() > 1;
        throw CoalError(joined + ": no sequences found: the file" + (many ? "s are" : " is") + " empty or ha" +
                            (many ? "ve" : "s") + " no '>' header lines (expected the coding DNA in FASTA format)",
                        "no_sequences");
    }
    rec.set_dna_inputs(dna.counts_per_file(), dna.ids());
    // CPP-25 (reference MRG-1): exact IDs; else (--strict: a partial exact match is an error) normalised IDs only
    // where file order would give a protein another gene's DNA, with a notice; else file order (ENG-d notice when
    // some IDs match, silent when none do)
    std::vector<std::string> per_seq, pair_notes;
    std::string mode = pair_ids(ids, dna.ids(), dna.pos(), rs.pairing, per_seq, pair_notes, o.strict);
    const std::vector<size_t>& pairing = rs.pairing;
    for (auto& n : pair_notes) warn(n);
    rec.set_pairing(mode, per_seq, pairing, pair_notes);
    // CHANGED(ENG-UX-18): sequences paired by file order whose IDs differ: a failure names the DNA record used
    // (pipeline.py _paired)
    std::string order_why;
    std::vector<std::string> order_dna;
    if (mode == "order") {
        long long hits = 0;
        for (auto& p : ids) hits += dna.pos().count(p) ? 1 : 0;
        order_why = hits ? "because only " + std::to_string(hits) + " of " + std::to_string(ids.size()) +
                               " protein IDs are DNA IDs"
                         : std::string("because no protein ID matches a DNA ID");
        order_dna.resize(ids.size());
        for (size_t k = 0; k < ids.size(); ++k)
            if (dna.ids()[pairing[k]] != ids[k]) order_dna[k] = dna.ids()[pairing[k]];
    }
    auto paired = [&](long long k, const std::string& ecode, const std::string& msg) {
        if (!order_why.empty() && ecode == "inconsistent" && !order_dna[(size_t)k].empty())
            return msg + " (paired by file order with DNA '" + order_dna[(size_t)k] + "', " + order_why + ")";
        return msg;
    };
    dna.forget_ids();
    const long long nseq = (long long)ids.size();
    long long wtotal = 0;
    for (int w : widths) wtotal += w;

    Ctx ctx{code, widths, staged, o.fmt, o.strict};
    // CPP-27: the report's end-of-run extras read the protein rows and the paired DNA again
    RowsFn rows_fn = [&](const std::function<bool(long long, const std::string&, const std::string&)>& fn) {
        ProteinReader rr(o.pep);
        std::string id, row;
        for (long long k = 0; rr.next(id, row); ++k)
            if (!fn(k, id, prepare_row(row))) break;
    };
    DnaFn dna_fn = [&](long long k) { return dna.get(pairing[(size_t)k]); };

    // ---- files: deterministic names next to -o (resumable), temp files for stdout
    std::unique_ptr<Checkpoint> ck;
    Json state;
    bool resuming = false;
    std::string out_tmp, stage_path, prot_path;
    // CPP-16 (reference ENG-a, R2-E5, R3-E2): -o PATH is checked before any work. Symlinks are written through;
    // FIFOs, devices and files with other hard links are written in place (from a temp file, no checkpoint).
    OutTarget target;
    g_out_arg = o.out;
    if (!o.out.empty()) target = resolve_output(o.out);
    // CHANGED(ENG-UX-12): -o naming an input file overwrote that input; now refused before any work (pipeline.py)
    if (!o.out.empty() && target.existing) {
        std::vector<std::string> ins{o.pep};
        ins.insert(ins.end(), o.nucs.begin(), o.nucs.end());
        for (auto& p : ins) {
            std::error_code eq;
            if (fs::equivalent(p, target.real, eq) && !eq)
                throw CoalError("cannot write " + o.out + ": it is one of the input files; choose another name for the "
                                "output",
                                "file_error");
        }
    }
    const bool direct = !o.out.empty() && target.direct;
    // CPP-28 (reference MRG-4): a checkpoint only for an output that is (or will be) a regular file - hard-linked
    // files included (written in place at the end) - and inputs that are regular files
    std::error_code ex_ec;
    const bool want_ck = !o.out.empty() && o.resumable && (target.existing || !fs::exists(target.real, ex_ec));
    bool ck_registered = false;
    if (want_ck) {
        std::string canon = "fmt=" + o.fmt + "|table=" + std::to_string(o.table) + "|nogap=" + std::to_string(o.nogap) +
                            "|nomismatch=" + std::to_string(o.nomismatch) + "|blockonly=" + std::to_string(o.blockonly) + "|strict=" + std::to_string(o.strict) + "|keep_going=" +
                            std::to_string(o.keep_going) + "|staged=" + std::to_string(staged);
        std::vector<std::string> inputs{o.pep};
        inputs.insert(inputs.end(), o.nucs.begin(), o.nucs.end());
        ck = std::make_unique<Checkpoint>(target.real, canon, inputs);
        out_tmp = ck->partial;
        if (staged) stage_path = ck->stage;
        if (o.fmt == "codon") prot_path = ck->prot;
        if (o.resume) {
            resuming = ck->load_for_resume(state, {staged ? stage_path : out_tmp, ck->entries, prot_path});
            if (resuming)
                note("resuming " + o.out + ": " + fmt_n(state.get("completed")->as_int()) + " of " + fmt_n(nseq) +
                     " sequences were already done" +
                     (staged ? " (this format/option needs a second pass, which is redone in full)" : ""));
            else
                note((ck->why == "no checkpoint" ? std::string("no checkpoint to resume")
                                                  : "found " + ck->why + ", so nothing can be resumed") +
                     " for " + o.out + "; starting a normal run");
        }
        if (!resuming) {
            ck->remove(true);  // a stale checkpoint must never describe this run's files
            // CPP-28 (reference MRG-4 / ENG-f): until a checkpoint describes them, a signal removes these files
            for (auto* q : {&ck->partial, &ck->stage, &ck->prot, &ck->entries}) register_temp(*q);
            register_temp(ck->stage + ".grid");
            ck_registered = true;
        }
    } else {
        if (o.resume && !o.out.empty())
            note(o.out + " cannot be resumed (" +
                 (o.resumable ? "the output is not a regular file" : "an input is not a regular file") +
                 "); starting a normal run");
        // the output's temp file goes next to the real target when it is renamed onto it (same file system)
        fs::path d = (!o.out.empty() && !direct) ? fs::path(target.real).parent_path() : fs::temp_directory_path();
        out_tmp = (d / (".coal_" + unique_suffix() + ".tmp")).string();
        if (staged) stage_path = (d / (".coal_rows_" + unique_suffix() + ".tmp")).string();
        if (o.fmt == "codon") prot_path = (d / (".coal_prot_" + unique_suffix() + ".tmp")).string();
        // CPP-18: removed by the signal handlers too (ENG-f, R2-E8)
        for (auto* q : {&out_tmp, &stage_path, &prot_path}) register_temp(*q);
        if (!stage_path.empty()) register_temp(stage_path + ".grid");
    }
    const std::string sink_path = staged ? stage_path : out_tmp;
    bool finished = false, failed = false;
    struct Cleanup {
        std::function<void()> fn;
        ~Cleanup() { fn(); }
    } cleanup{[&]() {
        if (ck) {
            std::error_code ec;
            // otherwise (crash, kill, disk full ...) keep partial + checkpoint for --resume, if one was written
            if (finished || failed || !fs::exists(ck->path, ec)) ck->remove();
            if (ck_registered) {
                for (auto* q : {&ck->partial, &ck->stage, &ck->prot, &ck->entries}) unregister_temp(*q);
                unregister_temp(ck->stage + ".grid");
            }
        } else {
            std::error_code ec;
            for (auto* p : {&out_tmp, &stage_path, &prot_path})
                if (!p->empty()) {
                    fs::remove(*p, ec);
                    unregister_temp(*p);
                }
            if (!stage_path.empty()) {
                fs::remove(stage_path + ".grid", ec);
                unregister_temp(stage_path + ".grid");
            }
        }
    }};

    try {
        auto get = [&](const char* k) -> long long {
            const Json* v = state.get(k);
            return v ? v->as_int() : 0;
        };
        long long start = resuming ? get("completed") : 0;
        OutFile sink, entries, prot;
        sink.open(sink_path, resuming ? get("sink_bytes") : -1);
        if (ck) entries.open(ck->entries, resuming ? get("entries_bytes") : -1);
        if (!prot_path.empty()) prot.open(prot_path, resuming ? get("prot_bytes") : -1);
        std::vector<char> drop((o.nomismatch || block_cols) ? (size_t)ncol : 0, 0), bad;
        long long skipped = 0;
        std::set<long long> skip_idx;       // sequences skipped by --keep-going (resume)
        if (resuming) {
            skipped = get("skipped");
            if (const Json* v = state.get("drop"))
                for (auto& x : v->a)
                    if (x.as_int() >= 0 && x.as_int() < (long long)drop.size()) drop[(size_t)x.as_int()] = 1;
            if (const Json* v = state.get("bad"))
                for (auto& x : v->a) {
                    long long t = x.as_int();
                    if (t >= 0) {
                        if ((size_t)t >= bad.size()) bad.resize((size_t)t + 1, 0);
                        bad[(size_t)t] = 1;
                    }
                }
            std::vector<Json> nonok;
            {
                std::ifstream in(ck->entries, std::ios::binary);
                std::string line;
                while (std::getline(in, line))
                    if (!line.empty()) nonok.push_back(json_parse(line));
            }
            const Json* rst = state.get("rec");
            rec.restore(rst ? *rst : Json::object(), nonok, start);
            if (const Json* v = state.get("skip_idx"))
                for (auto& x : v->a) skip_idx.insert(x.as_int());
        } else if (o.fmt == "paml" && !staged) {
            sink.f << paml_header(nseq, wtotal);
        }

        // CHANGED(IMP-E-B): the genetic-code suggestion before the per-codon warnings (pipeline.py early_code_notice)
        if (!o.quiet && !o.nomismatch) {
            std::string hint = early_code_notice(rows_fn, dna_fn, o.table, nseq, rec.fit_memo);
            if (!hint.empty()) note(hint);
        }
        // pass 2: batches of ~BATCH_NT nucleotides; thread batch k while reading batch k+1
        ProteinReader pr(o.pep);
        {
            std::string id, row;
            for (long long k = 0; k < start && pr.next(id, row); ++k) {
            }
        }
        size_t k = (size_t)start;
        auto read_batch = [&](std::vector<Item>& b) {
            b.clear();
            long long nt = 0;
            std::string id, row;
            while (pr.next(id, row)) {
                Item it;
                it.id = id;
                it.row = prepare_row(row);
                it.dna = dna.get(pairing[k++]);
                nt += (long long)it.dna.size();
                b.push_back(std::move(it));
                if (nt >= BATCH_NT) break;
            }
        };
        std::unique_ptr<Progress> prog;
        if (o.progress && nseq > PROGRESS_MIN_SEQS) prog = std::make_unique<Progress>(nseq);
        auto last_ck = std::chrono::steady_clock::now();
        long long done = start;
        // thread each batch while the next one is read; results written in input order (pool.hpp)
        // CPP-13: with -blockonly, as PAL2NAL, only warnings about kept ('#') columns are printed
        auto kept_col = [&](const std::string& w) {
            if (!block_cols || w.compare(0, 10, "pepAlnPos ") != 0) return true;
            size_t i = 10;
            long long c = 0;
            while (i < w.size() && w[i] >= '0' && w[i] <= '9') c = c * 10 + (w[i++] - '0');
            if (i == 10 || i >= w.size() || w[i] != ':') return true;
            return c >= 1 && c <= (long long)mask.size() && mask[(size_t)(c - 1)] == '#';
        };
        auto shown = [&](const std::vector<std::string>& ws) {  // the warnings -blockonly keeps (as on stderr)
            std::vector<std::string> v;
            for (auto& w : ws)
                if (kept_col(w)) v.push_back(w);
            return v;
        };
        auto emit = [&](std::vector<Item>& cur) {
            for (size_t i = 0; i < cur.size(); ++i)
                if (cur[i].failed) cur[i].error = paired(done + (long long)i, cur[i].error_code, cur[i].error);
            if (!o.keep_going)
                for (size_t i = 0; i < cur.size(); ++i)
                    if (cur[i].failed) {
                        // the sequences before it in its batch were checked (not written): report only
                        for (size_t j = 0; j < i; ++j)
                            rec.add_sequence(done + (long long)j, cur[j].id, shown(cur[j].warnings), cur[j].notices);
                        rec.add_sequence(done + (long long)i, cur[i].id, {}, {}, cur[i].error_code, cur[i].error);
                        throw CoalError(cur[i].error, cur[i].error_code, cur[i].id);
                    }
            for (Item& it : cur) {
                g_last_row = it.row;  // ENG-UX-20
                // CPP-26 (reference MRG-3): per sequence its warnings, its notices, then "skipped ..." (--keep-going)
                for (auto& w : it.warnings)  // CPP-11 (ENG-i2): no mismatch warnings with -nomismatch, as PAL2NAL
                    if (!(o.nomismatch && w.find(MISMATCH_WARNING) != std::string::npos) && kept_col(w))
                        warn(it.id + " " + w);
                for (auto& n : it.notices) warn(n);
                if (it.failed) warn("skipped " + it.error);
                const Json* e = it.failed ? rec.add_sequence(done, it.id, {}, {}, it.error_code, it.error)
                                          : rec.add_sequence(done, it.id, shown(it.warnings), it.notices);
                if (ck && e) entries.f << json_dump(*e) << '\n';
                if (it.failed) {
                    ++skipped;
                    skip_idx.insert(done);
                    ++done;
                    continue;
                }
                ++done;
                if (staged) {
                    sink.f.write(it.out.data(), (std::streamsize)it.out.size());
                    if (!prot_path.empty()) prot.f.write(it.row.data(), (std::streamsize)it.row.size());
                    if (o.nogap) mark_nogap_triplets(it.out, bad, *code);
                    if (o.nomismatch)
                        for (int c : it.mcols) drop[(size_t)(c - 1)] = 1;
                } else {
                    sink.f.write(it.out.data(), (std::streamsize)it.out.size());
                }
            }
            std::cerr.flush();
            if (prog) prog->update(done);
            if (ck && ck->count == 0 && done < nseq && test_hooks() && std::getenv("COAL_TEST_PAUSE_BEFORE_CHECKPOINT")) {
                sink.f.flush();  // tests (CPP-28): partial files exist, no checkpoint describes them yet
                std::cerr << "coal: test: paused before any checkpoint" << std::endl;
                std::this_thread::sleep_for(std::chrono::seconds(3600));
            }
            if (ck && done < nseq &&
                std::chrono::duration<double>(std::chrono::steady_clock::now() - last_ck).count() >= CHECKPOINT_EVERY) {
                Json st = Json::object();
                st.set("completed", done);
                st.set("nseq", nseq);
                st.set("sink_bytes", durable_size(sink.f, sink.path));
                st.set("entries_bytes", durable_size(entries.f, entries.path));
                if (!prot_path.empty()) st.set("prot_bytes", durable_size(prot.f, prot.path));
                st.set("skipped", skipped);
                Json dl = Json::array(), bl = Json::array();
                for (size_t i = 0; i < drop.size(); ++i)
                    if (drop[i]) dl.push((long long)i);
                for (size_t i = 0; i < bad.size(); ++i)
                    if (bad[i]) bl.push((long long)i);
                st.set("drop", dl);
                st.set("bad", bl);
                Json sk = Json::array();
                for (long long x : skip_idx) sk.push(x);
                st.set("skip_idx", sk);
                st.set("rec", rec.checkpoint_state());
                ck->save(st);
                if (ck_registered) {  // CPP-28: a checkpoint describes them now: kept after a signal
                    for (auto* q : {&ck->partial, &ck->stage, &ck->prot, &ck->entries}) unregister_temp(*q);
                    unregister_temp(ck->stage + ".grid");
                    ck_registered = false;
                }
                last_ck = std::chrono::steady_clock::now();
                test_pause(ck->count);
            }
        };
        run_ordered<Item>(read_batch, [&](Item& it) { work(it, ctx); }, emit, o.jobs);
        if (prog && prog->printed) prog->update(done, true);
        if (o.keep_going && skipped == nseq)
            throw CoalError("none of the " + std::to_string(nseq) + " sequences could be converted: check that the DNA "
                            "file holds the coding sequences of these proteins and that the genetic code (-codontable) "
                            "is right",
                            "all_failed");  // CHANGED(ENG-UX-17): + what to check
        rec.finish(rows_fn, dna_fn);  // CPP-27: break locations, mismatch patterns, genetic-code hint
        const long long nout = nseq - skipped;
        // CPP-33: the final pass names the written sequences: the report's IDs, or (--keep-going skipped some)
        // a list without the skipped ones
        std::vector<std::string> kept_ids;
        if (staged && !skip_idx.empty()) {
            kept_ids.reserve((size_t)nout);
            for (long long q = 0; q < nseq; ++q)
                if (!skip_idx.count(q)) kept_ids.push_back(ids[(size_t)q]);
        }
        const std::vector<std::string>& out_ids = skip_idx.empty() ? ids : kept_ids;

        OutFile out;
        if (staged) {
            out.open(out_tmp);
            std::vector<std::pair<long long, long long>> sl;
            if (block_cols)  // CPP-13: -blockonly drops the columns not marked '#'
                for (long long i = 0; i < ncol; ++i)
                    if (i >= (long long)mask.size() || mask[(size_t)i] != '#') drop[(size_t)i] = 1;
            if (two_pass) sl = keep_slices(widths, drop, bad);
            else if (wtotal > 0) sl.push_back({0, wtotal});
            long long length = 0;
            for (auto& [a, b] : sl) length += b - a;
            bool identity = sl.empty() ? wtotal == 0 : (sl.size() == 1 && sl[0].first == 0 && sl[0].second == wtotal);
            sink.f.flush();
            std::string row;
            if (!block) {
                if (o.fmt == "paml") out.f << paml_header(nout, length);
                std::string text;
                sink.f.seekg(0);
                for (long long i = 0; i < nout; ++i) {
                    read_next(sink.f, wtotal, row);
                    text.clear();
                    record_into(text, o.fmt, out_ids[(size_t)i], identity ? row : apply_slices(row, sl));
                    out.f << text;
                }
            } else {
                OutFile grid;
                std::fstream* g = &sink.f;
                if (!identity) {
                    std::string gpath = stage_path + ".grid";  // removed by the cleanup on any exit
                    grid.open(gpath);
                    sink.f.seekg(0);
                    for (long long i = 0; i < nout; ++i) {
                        read_next(sink.f, wtotal, row);
                        std::string s2 = apply_slices(row, sl);
                        grid.f.write(s2.data(), (std::streamsize)s2.size());
                    }
                    grid.f.flush();
                    g = &grid.f;
                }
                if (!prot_path.empty()) prot.f.flush();
                std::string mask_exp;  // the mark line under each block, one character per output nucleotide
                if (has_marks && !o.blockonly) {
                    std::string full;
                    for (long long i = 0; i < ncol; ++i)
                        full.append((size_t)widths[(size_t)i], i < (long long)mask.size() && mask[(size_t)i] == '#' ? '#' : ' ');
                    mask_exp = identity ? full : apply_slices(full, sl);
                }
                write_blocks(out.f, o.fmt, *g, length, prot_path.empty() ? nullptr : &prot.f, ncol, out_ids, mask_exp);
                if (!identity) {
                    grid.f.close();
                    std::error_code ec;
                    fs::remove(grid.path, ec);
                }
            }
            out.f.flush();
            if (!out.f) throw CoalError("could not write the output (disk full?)");
            out.f.close();
        } else {
            sink.f.flush();
            if (!sink.f) throw CoalError("could not write the output (disk full?)");
        }
        sink.f.close();
        entries.f.close();
        prot.f.close();
        // hashing a big output takes time: only when a report was asked for
        rec.set_output(o.out, rec.keep_entries ? sha256_file(out_tmp) : std::string(), nout);
        if (!o.out.empty() && !direct) {
            std::error_code ec;
            if (target.existing) copy_mode(target.real, out_tmp);  // CPP-16 (R2-E5): keep its permissions
            fs::rename(out_tmp, target.real, ec);
            if (ec) throw CoalError(output_error(o.out, ec.value()), "file_error");  // CHANGED(ENG-UX-12)
        } else if (direct) {  // CPP-16: FIFO / device / hard-linked file: written in place
            std::ifstream in(out_tmp, std::ios::binary);
            errno = 0;
            std::ofstream dst(o.out, std::ios::binary | std::ios::trunc);
            if (!dst) throw CoalError(output_error(o.out, errno ? errno : EIO), "file_error");  // CHANGED(ENG-UX-12)
            dst << in.rdbuf();
            dst.flush();
            if (!dst) throw CoalError("could not write " + o.out + " (disk full?)", "file_error");
        } else {
#ifdef _WIN32
            _setmode(_fileno(stdout), _O_BINARY);
#endif
            std::ifstream in(out_tmp, std::ios::binary);
            std::vector<char> buf(1 << 20);
            for (;;) {
                in.read(buf.data(), (std::streamsize)buf.size());
                std::streamsize n = in.gcount();
                if (n <= 0) break;
                write_stdout(buf.data(), (size_t)n);  // CPP-18: closed pipe -> BrokenPipe (exit 141)
            }
            flush_stdout();
        }
        finished = true;
        return skipped;
    } catch (const CoalError&) {
        failed = true;  // a real error: nothing to resume
        if (!rec.finished) rec.finish(rows_fn, dna_fn, true);
        throw;
    }
}


// ------------------------------------------------------------------ CLI

// -------------------------------------------------------------- command line (CPP-4)
// Same options, help and error messages as `python -m coal` on the reviewed branch (argparse); see
// args.hpp. Differences (documented in overnight/core/CHANGES.md): -output also accepts clustal, paml
// and codon; -strict is accepted but not implemented yet; the C++-only extras match exact names only.

static std::string usage_text() {
    std::string h = COAL_HELP;
    return h.substr(0, h.find("\n\n") + 1);
}

static int usage_error(const std::string& msg) {
    std::cerr << usage_text() << "coal: error: " << msg << '\n';
    return 2;
}

// argparse-normalised decimal -> int (saturating; the values are thread counts / warning limits)
static int to_int(const std::string& norm) {
    bool neg = !norm.empty() && norm[0] == '-';
    std::string d = neg ? norm.substr(1) : norm;
    if (d.size() > 9) return neg ? INT_MIN : INT_MAX;
    int v = std::stoi(d);
    return neg ? -v : v;
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    Options o;
    for (int i = 1; i < argc; ++i) o.argv.push_back(argv[i]);

    using coal::args::Kind;
    using coal::args::Spec;
    coal::args::Parser ap("coal", usage_text(), COAL_HELP, std::string("coal ") + VERSION + " (C++)");
    std::vector<std::string> tables;
    for (int t : supported_tables()) tables.push_back(std::to_string(t));
    ap.add({Kind::Help, {"-h", "--help"}, "", false, {}, true});
    const int a_pep = ap.add({Kind::Positional, {}, "pep_aln", false, {}, true});
    const int a_nuc = ap.add({Kind::PositionalPlus, {}, "nuc_fasta", false, {}, true});
    const int a_table = ap.add({Kind::Value, {"-codontable", "--codontable"}, "", true, tables, true});
    const int a_nogap = ap.add({Kind::Flag, {"-nogap", "--nogap"}, "", false, {}, true});
    const int a_nomis = ap.add({Kind::Flag, {"-nomismatch", "--nomismatch"}, "", false, {}, true});
    const int a_fmt =
        ap.add({Kind::Value, {"-output", "--output-format"}, "", false, {"fasta", "clustal", "paml", "codon"}, true});
    const int a_block = ap.add({Kind::Flag, {"-blockonly", "--blockonly"}, "", false, {}, true});  // CPP-13
    const int a_out = ap.add({Kind::Value, {"-o", "--out"}, "", false, {}, true});
    const int a_jobs = ap.add({Kind::Value, {"-j", "--jobs"}, "", true, {}, true});
    const int a_quiet = ap.add({Kind::Flag, {"-q", "--quiet"}, "", false, {}, true});
    const int a_strict = ap.add({Kind::Flag, {"-strict", "--strict"}, "", false, {}, true});
    // CPP-29 (reference MRG-3/MRG-4 cli.py): the extras are ordinary options of the merged parser, registered in its
    // order (before --version) and abbreviated like every other option (`--j` is now ambiguous: --jobs, --json)
    const int a_keep = ap.add({Kind::Flag, {"--keep-going", "-keep-going"}, "", false, {}, true});
    const int a_maxw = ap.add({Kind::Value, {"--max-warnings"}, "", true, {}, true});
    const int a_prog = ap.add({Kind::Flag, {"--progress"}, "", false, {}, true});
    const int a_resume = ap.add({Kind::Flag, {"--resume"}, "", false, {}, true});
    const int a_report = ap.add({Kind::Value, {"--report"}, "", false, {}, true});
    const int a_json = ap.add({Kind::Flag, {"--json"}, "", false, {}, true});
    ap.add({Kind::Version, {"--version"}, "", false, {}, true});

    coal::args::Result args = ap.parse(o.argv);
    if (args.exit) {
        std::cout << args.out << std::flush;
        std::cerr << args.err << std::flush;
        return args.code;
    }
    o.pep = args.get(a_pep);
    o.nucs = args.values[(size_t)a_nuc];
    if (args.given(a_table)) o.table = to_int(args.get(a_table));
    o.nogap = args.given(a_nogap);
    o.nomismatch = args.given(a_nomis);
    o.blockonly = args.given(a_block);
    if (args.given(a_fmt)) o.fmt = args.get(a_fmt);
    if (args.given(a_out)) o.out = args.get(a_out);
    if (args.given(a_jobs)) o.jobs = to_int(args.get(a_jobs));
    o.quiet = args.given(a_quiet);
    o.keep_going = args.given(a_keep);
    if (args.given(a_maxw)) o.max_warnings = to_int(args.get(a_maxw));
    o.progress = args.given(a_prog);
    o.resume = args.given(a_resume);
    if (args.given(a_report)) o.report = args.get(a_report);
    o.json = args.given(a_json);
    o.strict = args.given(a_strict);  // CPP-15 (reference DEC-1)
    if (o.json && o.out.empty())
        return usage_error("--json requires -o (the alignment and the report cannot both go to stdout)");
    if (o.resume && o.out.empty())
        return usage_error("--resume requires -o (only runs that write to a file can be resumed)");
    if (o.max_warnings < 0) return usage_error("--max-warnings must be 0 (no limit) or more");  // CPP-29 (MRG-3)

    RunReport rec(o.argv, o.pep, o.nucs, o.out, o.fmt, o.table, o.nogap, o.nomismatch, o.keep_going,
                  !o.report.empty() || o.json, o.strict, o.blockonly);
    rec.want_hint = !o.quiet || !o.report.empty() || o.json;  // CHANGED(IMP-E-D)
    install_signal_handlers();  // CPP-18 (reference ENG-f, R2-E8, DEC-3)
    // CPP-17 (reference ENG-b, R2-E7): pipes, /dev/stdin and FIFOs are copied to temp files first; messages name
    // the argument, not the temp copy
    std::vector<std::pair<std::string, std::string>> spooled;
    struct SpoolCleanup {
        std::vector<std::pair<std::string, std::string>>& s;
        ~SpoolCleanup() { remove_spooled(s); }
    } spool_cleanup{spooled};
    auto restore_names = [&](std::string msg) {
        for (auto& sp : spooled)
            for (size_t at; (at = msg.find(sp.first)) != std::string::npos;) msg.replace(at, sp.first.size(), sp.second);
        return msg;
    };
    long long shown = 0, hidden = 0;
    WarnFn warn = [&](const std::string& w) {
        if (o.quiet) return;
        if (o.max_warnings > 0 && w.find("pepAlnPos") != std::string::npos) {
            if (shown >= o.max_warnings) {
                ++hidden;
                return;
            }
            ++shown;
        }
        std::cerr << "coal: warning: " << w << '\n';
    };
    WarnFn note = [&](const std::string& t) {
        // CHANGED(ENG-UX-10): was "coal: note:"; -q: no notices (as the Python CLI)
        if (!o.quiet) std::cerr << "coal: notice: " << t << std::endl;
    };
    RunState rs;
    int rc = 0;
    Options ro = o;  // the same options, reading the (possibly spooled) input copies
    try {
        std::vector<std::string> in_args{o.pep};
        in_args.insert(in_args.end(), o.nucs.begin(), o.nucs.end());
        std::vector<std::string> paths = spool_inputs(in_args, spooled);
        ro.pep = paths[0];
        ro.nucs.assign(paths.begin() + 1, paths.end());
        ro.resumable = spooled.empty();
        rec.set_hash_paths(paths);
        long long skipped = run(ro, rec, rs, warn, note);
        if (skipped) rc = 2;  // partial success
    } catch (const BrokenPipe&) {
        return 141;  // CPP-18 (reference DEC-3): a closed stdout pipe ends the run silently; temp files are gone
    } catch (const CoalError& e) {
        std::cerr << "coal: error: " << restore_names(e.what()) << '\n';
        rc = 1;
        rec.fail(e.code, restore_names(e.what()), !e.seq_id.empty());
    } catch (const InputError& e) {
        std::cerr << "coal: error: " << restore_names(e.what()) << '\n';
        rc = 1;
        rec.fail(e.code, restore_names(e.what()), false);
    } catch (const std::exception& e) {
        std::cerr << "coal: error: " << restore_names(e.what()) << '\n';
        rc = 1;
        rec.fail("file_error", restore_names(e.what()), false);
    }
    if (hidden) {
        std::string where = (!o.report.empty() || o.json) ? "see the report"
                                                           : "see them all with --max-warnings 0 or --report FILE.json";
        std::cerr << "coal: warning: ... and " << fmt_n(hidden) << " more warnings (" << where << ")\n";
    }
    rec.exit_code = rc;
    if (!o.quiet)
        for (auto& l : rec.summary_lines()) std::cerr << "coal: summary: " << l << '\n';
    std::cerr.flush();
    if (!o.report.empty()) {
        std::string tmp = o.report + ".tmp" + unique_suffix();
        bool ok;
        int err = 0;
        {
            errno = 0;
            std::ofstream f(tmp, std::ios::binary | std::ios::trunc);
            f << json_dump(rec.to_dict(), 1) << '\n';
            f.flush();
            ok = (bool)f;
            if (!ok) err = errno ? errno : EIO;
        }
        std::error_code ec;
        if (ok) fs::rename(tmp, o.report, ec);
        if (!ok || ec) {
            if (ec) err = ec.value();
            std::error_code rm;
            fs::remove(tmp, rm);
            // CHANGED(ENG-UX-12): was "could not write the report: <path>" (no reason); now as cli.py
            std::string why = fs::is_directory(o.report, rm) ? "it is a folder, not a file" : std::strerror(err);
            std::cerr << "coal: error: could not write the report " << o.report << ": " << why << '\n';
            if (!rc) rc = 1;
        }
    }
    if (o.json) std::cout << json_dump(rec.to_dict(), 1) << '\n' << std::flush;
    return rc;
}
