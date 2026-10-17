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

#include "core.hpp"
#include "json.hpp"
#include "report.hpp"
#include "seqio.hpp"
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

static const char* VERSION = "0.3.0";
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
static const long long BLOCK_CACHE = 64LL << 20;  // bytes of rows held while writing blocks

struct Options {
    std::string pep;
    std::vector<std::string> nucs;
    int table = 1;
    bool nogap = false, nomismatch = false, quiet = false;
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

struct TempFile {
    std::string path;
    std::fstream f;
    std::vector<char> buf;
    explicit TempFile(const std::string& dir, const std::string& prefix) : buf(1 << 20) {
        fs::path d = dir.empty() ? fs::temp_directory_path() : fs::path(dir);
        path = (d / (prefix + unique_suffix() + ".tmp")).string();
        f.rdbuf()->pubsetbuf(buf.data(), (std::streamsize)buf.size());
        f.open(path, std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);
        if (!f) throw InputError("[Errno 2] No such file or directory: " + py_repr(path));
    }
    void close() {
        if (f.is_open()) f.close();
    }
    ~TempFile() {
        close();
        std::error_code ec;
        if (!path.empty()) fs::remove(path, ec);
    }
};

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
    std::snprintf(b, sizeof b, "%4lld%7lld\n", n, len);
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
    bool failed = false;
    std::string error, error_code;
};

struct Ctx {
    const GeneticCode* code;
    std::vector<int> widths;
    bool staged;
    std::string fmt;
};

static void work(Item& it, const Ctx& c) {
    try {
        Threaded t = thread_one(it.id, it.row, it.dna, *c.code);
        std::string exp = expand(it.row, t.dna, c.widths);
        it.warnings = std::move(t.warnings);
        it.mcols = std::move(t.mismatch_cols);
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

static void process(std::vector<Item>& batch, const Ctx& c, int nthreads) {
    size_t n = batch.size();
    int k = (int)std::min<size_t>((size_t)std::max(1, nthreads), n);
    if (k <= 1) {
        for (Item& it : batch) work(it, c);
        return;
    }
    std::atomic<size_t> next{0};
    std::vector<std::thread> pool;
    for (int t = 0; t < k; ++t)
        pool.emplace_back([&]() {
            size_t i;
            while ((i = next++) < n) work(batch[i], c);
        });
    for (auto& th : pool) th.join();
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
static void write_blocks(std::ostream& out, const std::string& fmt, std::fstream& grid, long long length,
                         std::fstream* prot, long long ncol, const std::vector<std::string>& ids) {
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
                text += '\n';
                out << text;
            }
        }
        return;
    }
    // codon: residue line and codon line cut into blocks of 20 separately (as PAL2NAL)
    long long ntrip = length / 3;
    long long nb = std::max((ncol + 19) / 20, (ntrip + 19) / 20);
    long long per = std::max(1LL, BLOCK_CACHE / std::max(1LL, n * 80));
    std::vector<std::string> rows((size_t)n), prows((size_t)n);
    std::string pad(width, ' ');
    for (long long sb = 0; sb < nb; sb += per) {
        long long d0 = std::min(length, sb * 60), d1 = std::min(length, (sb + per) * 60);
        long long p0 = std::min(ncol, sb * 20), p1 = std::min(ncol, (sb + per) * 20);
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
        for (long long b = sb; b < std::min(nb, sb + per); ++b) {
            long long da = std::min(length, b * 60) - d0, de = std::min(length, (b + 1) * 60) - d0;
            long long pa = std::min(ncol, b * 20) - p0, pe = std::min(ncol, (b + 1) * 20) - p0;
            text.clear();
            for (long long k = 0; k < n; ++k) {
                text += pad;
                for (long long i = pa; i < pe; ++i) {
                    if (i > pa) text += "   ";
                    text += prows[(size_t)k][(size_t)i];
                }
                text += '\n';
                text += names[(size_t)k];
                for (long long i = da; i < de; i += 3) {
                    if (i > da) text += ' ';
                    text.append(rows[(size_t)k], (size_t)i, (size_t)std::min(3LL, de - i));
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
        if (!in) return false;
        std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        try {
            st = json_parse(text);
        } catch (...) {
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
    std::fstream f;
    std::vector<char> buf;
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
    if (o.fmt == "codon" && (o.nogap || o.nomismatch))
        throw CoalError("\"-output codon\" cannot be combined with -nogap or -nomismatch (as in PAL2NAL)");
    const bool block = o.fmt == "clustal" || o.fmt == "codon";
    const bool two_pass = o.nogap || o.nomismatch;
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
    if (mode == "order") warn(order_pairing_warning(ids, dna.ids()));
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

    Ctx ctx{code, widths, staged, o.fmt};
    int nthreads = o.jobs > 0 ? o.jobs : (int)std::max(1u, std::thread::hardware_concurrency());

    // ---- files: deterministic names next to -o (resumable), temp files for stdout
    std::unique_ptr<Checkpoint> ck;
    Json state;
    bool resuming = false;
    std::string out_tmp, stage_path, prot_path;
    if (!o.out.empty()) {
        std::string canon = "fmt=" + o.fmt + "|table=" + std::to_string(o.table) + "|nogap=" + std::to_string(o.nogap) +
                            "|nomismatch=" + std::to_string(o.nomismatch) + "|keep_going=" +
                            std::to_string(o.keep_going) + "|staged=" + std::to_string(staged);
        std::vector<std::string> inputs{o.pep};
        inputs.insert(inputs.end(), o.nucs.begin(), o.nucs.end());
        ck = std::make_unique<Checkpoint>(o.out, canon, inputs);
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
                note("no checkpoint to resume for " + o.out + "; starting a normal run");
        }
        if (!resuming) ck->remove(true);  // a stale checkpoint must never describe this run's files
    } else {
        fs::path d = fs::temp_directory_path();
        out_tmp = (d / (".coal_" + unique_suffix() + ".tmp")).string();
        if (staged) stage_path = (d / (".coal_rows_" + unique_suffix() + ".tmp")).string();
        if (o.fmt == "codon") prot_path = (d / (".coal_prot_" + unique_suffix() + ".tmp")).string();
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
                if (!p->empty()) fs::remove(*p, ec);
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
        std::vector<char> drop(o.nomismatch ? (size_t)ncol : 0, 0), bad;
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
        std::vector<Item> cur, nxt;
        read_batch(cur);
        while (!cur.empty()) {
            std::thread worker([&]() { process(cur, ctx, nthreads); });
            try {
                read_batch(nxt);
            } catch (...) {
                worker.join();
                throw;
            }
            worker.join();
            if (!o.keep_going)
                for (size_t i = 0; i < cur.size(); ++i)
                    if (cur[i].failed) {
                        rec.add_sequence(cur[i].id, {}, cur[i].error_code, cur[i].error, done + (long long)i);
                        throw CoalError(cur[i].error, cur[i].error_code, cur[i].id);
                    }
            for (Item& it : cur) {
                for (auto& w : it.warnings) warn(it.id + " " + w);
                if (it.failed) warn("skipped " + it.error);
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
                    if (o.nogap) mark_nogap_triplets(it.out, bad);
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
            std::swap(cur, nxt);
        }
        if (prog && prog->printed) prog->update(done, true);
        if (o.keep_going && skipped == nseq)
            throw CoalError("none of the " + std::to_string(nseq) + " sequences could be converted", "all_failed");
        const long long nout = nseq - skipped;

        OutFile out;
        if (staged) {
            out.open(out_tmp);
            std::vector<std::pair<long long, long long>> sl;
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
                    std::string gpath = ck ? ck->stage + ".grid" : stage_path + ".grid";
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
                write_blocks(out.f, o.fmt, *g, length, prot_path.empty() ? nullptr : &prot.f, ncol, out_ids);
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
        if (!o.out.empty()) {
            std::error_code ec;
            fs::rename(out_tmp, o.out, ec);
            if (ec) throw CoalError("could not write " + o.out + ": " + ec.message(), "file_error");
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
                std::fwrite(buf.data(), 1, (size_t)n, stdout);
            }
            std::fflush(stdout);
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

static const char* USAGE =
    "usage: coal-cpp [-h] [-codontable N] [-nogap] [-nomismatch] [-output {fasta,clustal,paml,codon}]\n"
    "                [-o OUT] [-j N] [-q] [--keep-going] [--max-warnings N] [--progress] [--resume]\n"
    "                [--report FILE.json] [--json] [--version] pep_aln nuc_fasta [nuc_fasta ...]\n";

static const char* HELP =
    "\nConvert a protein alignment + unaligned coding DNA into a codon alignment (C++ engine of COAL;\n"
    "same options and output as `python -m coal`).\n\n"
    "  -codontable N       NCBI genetic code (default 1)\n"
    "  -nogap              remove codon columns containing gaps or stop codons\n"
    "  -nomismatch         remove alignment columns where a codon does not match its residue\n"
    "  -output FORMAT      fasta (default), clustal, paml or codon\n"
    "  -o, --out FILE      output file (default: stdout)\n"
    "  -j, --jobs N        worker threads (default 0 = all hardware threads)\n"
    "  -q, --quiet         do not print warnings or the summary\n"
    "  --keep-going        skip sequences that cannot be converted; exit code 2 if any were skipped\n"
    "  --max-warnings N    print at most N per-codon warning lines (default 0 = no limit)\n"
    "  --progress          percent done and time left on stderr (inputs over 10,000 sequences)\n"
    "  --resume            continue an interrupted -o run from its checkpoint (same inputs/options);\n"
    "                      with -nogap, -nomismatch, clustal or codon output only the first pass\n"
    "                      resumes and the final pass is redone\n"
    "  --report FILE.json  write a JSON run report (schema coal.report/1), also when the run fails\n"
    "  --json              print the run report to stdout (requires -o)\n";

static int usage_error(const std::string& msg) {
    std::cerr << USAGE << "coal-cpp: error: " << msg << '\n';
    return 2;
}

static bool parse_int(const std::string& s, int& v) {
    if (s.empty()) return false;
    size_t i = (s[0] == '-' || s[0] == '+') ? 1 : 0;
    if (i == s.size()) return false;
    for (size_t k = i; k < s.size(); ++k)
        if (s[k] < '0' || s[k] > '9') return false;
    try {
        v = std::stoi(s);
    } catch (...) {
        return false;
    }
    return true;
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    Options o;
    std::vector<std::string> pos;
    for (int i = 1; i < argc; ++i) o.argv.push_back(argv[i]);
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a.size() > 1 && a[0] == '-' && !(a[1] >= '0' && a[1] <= '9')) {
            std::string name = a, val;
            bool has_val = false;
            size_t eq = a.find('=');
            if (eq != std::string::npos) {
                name = a.substr(0, eq);
                val = a.substr(eq + 1);
                has_val = true;
            }
            auto take = [&](std::string& dst) -> bool {
                if (has_val) {
                    dst = val;
                    return true;
                }
                if (i + 1 >= argc) return false;
                dst = argv[++i];
                return true;
            };
            if (name == "-h" || name == "--help") {
                std::cout << USAGE << HELP;
                return 0;
            } else if (name == "--version") {
                std::cout << "coal " << VERSION << " (C++)\n";
                return 0;
            } else if (name == "-codontable" || name == "--codontable") {
                std::string v;
                int t;
                if (!take(v)) return usage_error("argument -codontable/--codontable: expected one argument");
                if (!parse_int(v, t)) return usage_error("argument -codontable/--codontable: invalid int value: '" + v + "'");
                if (!get_code(t)) return usage_error("argument -codontable/--codontable: invalid choice: " + v);
                o.table = t;
            } else if (name == "-nogap" || name == "--nogap") {
                o.nogap = true;
            } else if (name == "-nomismatch" || name == "--nomismatch") {
                o.nomismatch = true;
            } else if (name == "-output" || name == "--output-format") {
                std::string v;
                if (!take(v)) return usage_error("argument -output/--output-format: expected one argument");
                if (v != "fasta" && v != "clustal" && v != "paml" && v != "codon")
                    return usage_error("argument -output/--output-format: invalid choice: '" + v + "'");
                o.fmt = v;
            } else if (name == "-o" || name == "--out") {
                if (!take(o.out)) return usage_error("argument -o/--out: expected one argument");
            } else if (name == "-j" || name == "--jobs") {
                std::string v;
                if (!take(v) || !parse_int(v, o.jobs)) return usage_error("argument -j/--jobs: invalid int value");
            } else if (name == "-q" || name == "--quiet") {
                o.quiet = true;
            } else if (name == "--keep-going" || name == "-keep-going") {
                o.keep_going = true;
            } else if (name == "--max-warnings") {
                std::string v;
                if (!take(v) || !parse_int(v, o.max_warnings))
                    return usage_error("argument --max-warnings: invalid int value");
            } else if (name == "--progress") {
                o.progress = true;
            } else if (name == "--resume") {
                o.resume = true;
            } else if (name == "--report") {
                if (!take(o.report)) return usage_error("argument --report: expected one argument");
            } else if (name == "--json") {
                o.json = true;
            } else {
                return usage_error("unrecognized arguments: " + a);
            }
        } else {
            pos.push_back(a);
        }
    }
    if (pos.size() < 2) return usage_error("the following arguments are required: pep_aln, nuc_fasta");
    if (o.json && o.out.empty())
        return usage_error("--json requires -o (the alignment and the report cannot both go to stdout)");
    if (o.resume && o.out.empty())
        return usage_error("--resume requires -o (only runs that write to a file can be resumed)");
    o.pep = pos[0];
    o.nucs.assign(pos.begin() + 1, pos.end());

    RunReport rec(o.argv, o.pep, o.nucs, o.out, o.fmt, o.table, o.nogap, o.nomismatch, o.keep_going,
                  !o.report.empty() || o.json);
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
    WarnFn note = [](const std::string& t) { std::cerr << "coal: note: " << t << std::endl; };
    RunState rs;
    int rc = 0;
    try {
        long long skipped = run(o, rec, rs, warn, note);
        if (skipped) rc = 2;  // partial success
    } catch (const CoalError& e) {
        std::cerr << "coal: error: " << e.what() << '\n';
        rc = 1;
        rec.fail(e.code, e.what());
    } catch (const InputError& e) {
        std::cerr << "coal: error: " << e.what() << '\n';
        rc = 1;
        rec.fail(e.code, e.what());
    } catch (const std::exception& e) {
        std::cerr << "coal: error: " << e.what() << '\n';
        rc = 1;
        rec.fail("file_error", e.what());
    }
    if (hidden) {
        std::string where = (!o.report.empty() || o.json) ? "see the report" : "see them all with --report FILE.json";
        std::cerr << "coal: warning: ... and " << fmt_n(hidden) << " more warnings (" << where << ")\n";
    }
    rec.exit_code = rc;
    try {
        compute_hints(o, rec, rs);
    } catch (...) {
        // a hint is optional; never let it hide the real result
    }
    if (!o.quiet)
        for (auto& l : rec.summary_lines()) std::cerr << "coal: summary: " << l << '\n';
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
