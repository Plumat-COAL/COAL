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
// COAL_CPP_ONLY_LINES=0 (CPP-18, for comparing with the reference engine in tests): leave out the stderr lines only
// coal-cpp prints ("coal: summary:" lines, the no-ID-match pairing warning).
static const bool CPP_ONLY_LINES = env_ll("COAL_CPP_ONLY_LINES", 1) != 0;
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
    int max_warnings = 0;
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

// Returns mode "id", "normalised" or "order"; fills idx (and per_seq for "normalised").
static std::string pair_ids(const std::vector<std::string>& prot, const std::vector<std::string>& dna_ids,
                            const std::unordered_map<std::string, size_t>& pos, std::vector<size_t>& idx,
                            std::vector<std::string>& per_seq) {
    idx.assign(prot.size(), 0);
    per_seq.clear();
    bool all = true;
    for (auto& p : prot)
        if (!pos.count(p)) {
            all = false;
            break;
        }
    if (all) {
        for (size_t k = 0; k < prot.size(); ++k) idx[k] = pos.at(prot[k]);
        return "id";
    }
    std::unordered_map<std::string, std::vector<size_t>> by_key;
    for (size_t i = 0; i < dna_ids.size(); ++i)
        for (auto& key : id_keys(dna_ids[i])) {
            auto& v = by_key[key];
            if (v.empty() || v.back() != i) v.push_back(i);
        }
    std::vector<char> used(dna_ids.size(), 0);
    bool ok = true;
    for (size_t k = 0; k < prot.size() && ok; ++k) {
        size_t j;
        std::string m;
        auto it = pos.find(prot[k]);
        if (it != pos.end()) {
            j = it->second;
            m = "id";
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
            m = "normalised";
        }
        if (used[j]) {
            ok = false;
            break;
        }
        used[j] = 1;
        idx[k] = j;
        per_seq.push_back(m);
    }
    if (ok) return "normalised";
    for (size_t k = 0; k < prot.size(); ++k) idx[k] = k;
    per_seq.clear();
    return "order";
}

static std::string order_pairing_warning(const std::vector<std::string>& prot, const std::vector<std::string>& dna) {
    std::string pairs;
    for (size_t i = 0; i < std::min<size_t>(2, prot.size()); ++i) {
        if (i) pairs += ", ";
        pairs += prot[i] + (i == 0 ? " <-> DNA " : " <-> ") + dna[i];
    }
    return "IDs don't match, so sequences were paired by file order (protein " + pairs +
           (prot.size() > 2 ? ", ..." : "") + "). Check this is what you want.";
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

// CLUSTAL / codon output from fixed-stride row files, a few blocks at a time.
// mask: the '#' mark line (length characters) printed under each block, or empty (CPP-13).
static void write_blocks(std::ostream& out, const std::string& fmt, std::fstream& grid, long long length,
                         std::fstream* prot, long long ncol, const std::vector<std::string>& ids,
                         const std::string& mask) {
    long long n = (long long)ids.size();
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
    std::string last;
    if (n > 0) read_rows(*prot, (n - 1) * ncol, ncol, last);
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

static void test_pause(int count) {
    const char* v = std::getenv("COAL_TEST_PAUSE_AFTER_CHECKPOINT");
    if (v && std::to_string(count) == v) {
        std::cerr << "coal: test: paused after checkpoint " << count << std::endl;
        std::this_thread::sleep_for(std::chrono::seconds(3600));
    }
}

// A file written by the pipeline: fresh, or (resume) truncated to a byte length and appended to.
struct OutFile {
    std::string path;
    std::vector<char> buf;  // declared before f: f's destructor flushes into it
    std::fstream f;
    OutFile() : buf(1 << 20) {}
    void open(const std::string& p, long long resume_bytes = -1) {
        path = p;
        f.rdbuf()->pubsetbuf(buf.data(), (std::streamsize)buf.size());
        if (resume_bytes >= 0) {
            fs::resize_file(p, (uintmax_t)resume_bytes);
            f.open(p, std::ios::in | std::ios::out | std::ios::binary);
            f.seekp(0, std::ios::end);
        } else {
            f.open(p, std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);
        }
        if (!f) throw InputError("[Errno 2] No such file or directory: " + py_repr(p), "file_error");
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
    std::vector<std::string> ids;
    long long ncol = -1;
    std::vector<int> widths;
    {
        ProteinReader pr(o.pep);
        std::string id, row;
        while (pr.next(id, row)) {
            row = prepare_row(row);
            if (ncol < 0) {
                ncol = (long long)row.size();
                widths.assign((size_t)ncol, 3);
            } else if ((long long)row.size() != ncol) {
                throw CoalError("protein alignment rows have different lengths (" + std::to_string(ncol) + " and " +
                                    std::to_string(row.size()) + ", at " + id + ")",
                                "row_lengths");
            }
            update_column_widths(widths, row);
            ids.push_back(id);
        }
    }
    DnaIndex dna(o.nucs);
    rec.set_inputs((long long)ids.size(), dna.counts_per_file());
    if (ids.size() != dna.size())
        throw CoalError("number of sequences differs (protein: " + std::to_string(ids.size()) +
                            ", DNA: " + std::to_string(dna.size()) + ")",
                        "count_mismatch");
    if (ids.empty()) throw CoalError("no sequences in input", "no_sequences");
    std::vector<std::string> per_seq;
    std::string mode = pair_ids(ids, dna.ids(), dna.pos(), rs.pairing, per_seq);
    const std::vector<size_t>& pairing = rs.pairing;
    if (mode == "order") {
        // CPP-15 (reference ENG-d, DEC-1): some but not all protein IDs are DNA IDs: a notice (an error with
        // --strict). No ID matching at all keeps coal-cpp's own pairing warning (see LOG, needs reconciling).
        long long hits = 0;
        for (auto& p : ids) hits += dna.pos().count(p) ? 1 : 0;
        if (hits) {
            std::string p, d;
            for (size_t k = 0; k < ids.size(); ++k)
                if (ids[k] != dna.ids()[k]) {
                    p = ids[k];
                    d = dna.ids()[k];
                    break;
                }
            const std::string H = std::to_string(hits), T = std::to_string(ids.size());
            if (o.strict)
                throw CoalError("only " + H + " of " + T + " protein IDs are DNA IDs, so sequences would be paired by "
                                    "order (e.g. protein '" + p + "' with DNA '" + d + "'); make the IDs match, or "
                                    "run without --strict to pair by order as PAL2NAL does",
                                "partial_ids");
            warn("only " + H + " of " + T + " protein IDs are DNA IDs; pairing sequences by order instead (e.g. "
                 "protein '" + p + "' gets DNA '" + d + "') (--strict makes this an error)");
        } else {
            if (CPP_ONLY_LINES) warn(order_pairing_warning(ids, dna.ids()));
        }
    }
    std::string ex_prot, ex_dna;
    for (size_t k = 0; k < per_seq.size(); ++k)
        if (per_seq[k] == "normalised") {
            ex_prot = ids[k];
            ex_dna = dna.ids()[pairing[k]];
            break;
        }
    if (rec.keep_entries) {
        rs.dna_ids_paired.reserve(ids.size());
        for (size_t k = 0; k < ids.size(); ++k) rs.dna_ids_paired.push_back(dna.ids()[pairing[k]]);
    }
    rec.set_pairing(mode, rec.keep_entries ? &rs.dna_ids_paired : nullptr, &rs.pairing, per_seq, ex_prot, ex_dna);
    dna.forget_ids();
    const long long nseq = (long long)ids.size();
    long long wtotal = 0;
    for (int w : widths) wtotal += w;

    Ctx ctx{code, widths, staged, o.fmt, o.strict};

    // ---- files: deterministic names next to -o (resumable), temp files for stdout
    std::unique_ptr<Checkpoint> ck;
    Json state;
    bool resuming = false;
    std::string out_tmp, stage_path, prot_path;
    // CPP-16 (reference ENG-a, R2-E5, R3-E2): -o PATH is checked before any work. Symlinks are written through;
    // FIFOs, devices and files with other hard links are written in place (from a temp file, no checkpoint).
    OutTarget target;
    if (!o.out.empty()) target = resolve_output(o.out);
    const bool direct = !o.out.empty() && target.direct;
    if (!o.out.empty() && !direct) {
        std::string canon = "fmt=" + o.fmt + "|table=" + std::to_string(o.table) + "|nogap=" + std::to_string(o.nogap) +
                            "|nomismatch=" + std::to_string(o.nomismatch) + "|blockonly=" + std::to_string(o.blockonly) + "|keep_going=" +
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
        if (!resuming) ck->remove(true);  // a stale checkpoint must never describe this run's files
    } else {
        fs::path d = fs::temp_directory_path();
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
            if (finished || failed) ck->remove();
            // otherwise (crash, kill, disk full ...) keep partial + checkpoint for --resume
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
        std::vector<std::string> out_ids;  // staged formats: IDs of the written sequences, in order
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
            std::vector<Entry> nonok;
            {
                std::ifstream in(ck->entries, std::ios::binary);
                std::string line;
                while (std::getline(in, line))
                    if (!line.empty()) nonok.push_back(Entry::from_json(json_parse(line)));
            }
            const Json* rst = state.get("rec");
            rec.restore(rst ? *rst : Json::object(), nonok, start, ids);
            if (staged) {
                std::vector<char> is_failed((size_t)start, 0);
                for (auto& e : nonok)
                    if (e.status == "failed" && e.index >= 0 && e.index < start) is_failed[(size_t)e.index] = 1;
                for (long long k = 0; k < start; ++k)
                    if (!is_failed[(size_t)k]) out_ids.push_back(ids[(size_t)k]);
            }
        } else if (o.fmt == "paml" && !staged) {
            sink.f << paml_header(nseq, wtotal);
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
        auto emit = [&](std::vector<Item>& cur) {
            if (!o.keep_going)
                for (size_t i = 0; i < cur.size(); ++i)
                    if (cur[i].failed) {
                        rec.add_sequence(cur[i].id, {}, cur[i].error_code, cur[i].error, done + (long long)i);
                        throw CoalError(cur[i].error, cur[i].error_code, cur[i].id);
                    }
            // CPP-15: as the reference, a batch's per-codon warnings first, then its notices
            for (Item& it : cur) {
                for (auto& w : it.warnings)  // CPP-11 (ENG-i2): no mismatch warnings with -nomismatch, as PAL2NAL
                    if (!(o.nomismatch && w.find(MISMATCH_WARNING) != std::string::npos) && kept_col(w))
                        warn(it.id + " " + w);
                if (it.failed) warn("skipped " + it.error);
            }
            for (Item& it : cur)
                for (auto& n : it.notices) warn(n);
            for (Item& it : cur) {
                const Entry* e = rec.add_sequence(it.id, it.warnings, it.failed ? it.error_code : "",
                                                  it.failed ? it.error : "", done);
                if (ck && e && e->status != "ok") entries.f << json_dump(e->to_json()) << '\n';
                ++done;
                if (it.failed) {
                    ++skipped;
                    continue;
                }
                if (staged) {
                    sink.f.write(it.out.data(), (std::streamsize)it.out.size());
                    if (!prot_path.empty()) prot.f.write(it.row.data(), (std::streamsize)it.row.size());
                    if (o.nogap) mark_nogap_triplets(it.out, bad, *code);
                    if (o.nomismatch)
                        for (int c : it.mcols) drop[(size_t)(c - 1)] = 1;
                    out_ids.push_back(it.id);
                } else {
                    sink.f.write(it.out.data(), (std::streamsize)it.out.size());
                }
            }
            std::cerr.flush();
            if (prog) prog->update(done);
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
                st.set("rec", rec.checkpoint_state());
                ck->save(st);
                last_ck = std::chrono::steady_clock::now();
                test_pause(ck->count);
            }
        };
        run_ordered<Item>(read_batch, [&](Item& it) { work(it, ctx); }, emit, o.jobs);
        if (prog && prog->printed) prog->update(done, true);
        if (o.keep_going && skipped == nseq)
            throw CoalError("none of the " + std::to_string(nseq) + " sequences could be converted", "all_failed");
        const long long nout = nseq - skipped;

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
            if (ec) throw CoalError("could not write " + o.out + ": " + ec.message(), "file_error");
        } else if (direct) {  // CPP-16: FIFO / device / hard-linked file: written in place
            std::ifstream in(out_tmp, std::ios::binary);
            std::ofstream dst(o.out, std::ios::binary | std::ios::trunc);
            if (!dst) throw CoalError("cannot write " + o.out + ": " + std::strerror(errno), "file_error");
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
        throw;
    }
}

// Genetic-code hint: re-test the problem sequences with every supported code.
static void compute_hints(const Options& o, RunReport& rec, const RunState& rs) {
    if (!rec.needs_hint() || rs.pairing.empty()) return;
    std::vector<long long> want = rec.problems;
    std::sort(want.begin(), want.end());
    std::vector<std::pair<std::string, std::string>> rows;  // (id, row) for wanted positions in order
    std::vector<long long> pos;
    {
        ProteinReader pr(o.pep);
        std::string id, row;
        long long k = 0;
        size_t w = 0;
        while (w < want.size() && pr.next(id, row)) {
            if (k == want[w]) {
                rows.push_back({id, prepare_row(row)});
                pos.push_back(k);
                ++w;
            }
            ++k;
        }
    }
    DnaIndex dna(o.nucs);
    std::vector<std::string> seqs;
    for (long long k : pos) seqs.push_back(dna.get(rs.pairing[(size_t)k]));
    std::vector<std::pair<int, int>> fits;
    for (int t : supported_tables()) {
        const GeneticCode* g = get_code(t);
        int n = 0;
        for (size_t i = 0; i < rows.size(); ++i) {
            try {
                if (thread_one(rows[i].first, rows[i].second, seqs[i], *g).mismatch_cols.empty()) ++n;
            } catch (const CoalError&) {
            }
        }
        fits.push_back({t, n});
    }
    rec.apply_fits(fits);
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
    ap.add({Kind::Version, {"--version"}, "", false, {}, true});
    // C++-only options: exact names (no abbreviations), so they never make a Python command line ambiguous
    const int a_keep = ap.add({Kind::Flag, {"-keep-going", "--keep-going"}, "", false, {}, false});
    const int a_maxw = ap.add({Kind::Value, {"--max-warnings"}, "", true, {}, false});
    const int a_prog = ap.add({Kind::Flag, {"--progress"}, "", false, {}, false});
    const int a_resume = ap.add({Kind::Flag, {"--resume"}, "", false, {}, false});
    const int a_report = ap.add({Kind::Value, {"--report"}, "", false, {}, false});
    const int a_json = ap.add({Kind::Flag, {"--json"}, "", false, {}, false});

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

    RunReport rec(o.argv, o.pep, o.nucs, o.out, o.fmt, o.table, o.nogap, o.nomismatch, o.keep_going,
                  !o.report.empty() || o.json);
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
        if (!o.quiet) std::cerr << "coal: note: " << t << std::endl;  // -q: no notices (as the Python CLI)
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
        rec.set_hash_paths(paths);
        long long skipped = run(ro, rec, rs, warn, note);
        if (skipped) rc = 2;  // partial success
    } catch (const BrokenPipe&) {
        return 141;  // CPP-18 (reference DEC-3): a closed stdout pipe ends the run silently; temp files are gone
    } catch (const CoalError& e) {
        std::cerr << "coal: error: " << restore_names(e.what()) << '\n';
        rc = 1;
        rec.fail(e.code, restore_names(e.what()));
    } catch (const InputError& e) {
        std::cerr << "coal: error: " << restore_names(e.what()) << '\n';
        rc = 1;
        rec.fail(e.code, restore_names(e.what()));
    } catch (const std::exception& e) {
        std::cerr << "coal: error: " << restore_names(e.what()) << '\n';
        rc = 1;
        rec.fail("file_error", restore_names(e.what()));
    }
    if (hidden) {
        std::string where = (!o.report.empty() || o.json) ? "see the report" : "see them all with --report FILE.json";
        std::cerr << "coal: warning: ... and " << fmt_n(hidden) << " more warnings (" << where << ")\n";
    }
    rec.exit_code = rc;
    try {
        compute_hints(ro, rec, rs);
    } catch (...) {
        // a hint is optional; never let it hide the real result
    }
    if (!o.quiet)
        for (auto& l : rec.summary_lines())
            if (CPP_ONLY_LINES) std::cerr << "coal: summary: " << l << '\n';
    std::cerr.flush();
    if (!o.report.empty()) {
        std::string tmp = o.report + ".tmp" + unique_suffix();
        bool ok;
        {
            std::ofstream f(tmp, std::ios::binary | std::ios::trunc);
            f << json_dump(rec.to_json(), 1) << '\n';
            f.flush();
            ok = (bool)f;
        }
        std::error_code ec;
        if (ok) fs::rename(tmp, o.report, ec);
        if (!ok || ec) {
            fs::remove(tmp, ec);
            std::cerr << "coal: error: could not write the report: " << o.report << '\n';
            if (!rc) rc = 1;
        }
    }
    if (o.json) std::cout << json_dump(rec.to_json(), 1) << '\n' << std::flush;
    return rc;
}
