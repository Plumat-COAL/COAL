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
#include "seqio.hpp"
#include "tables.hpp"

namespace fs = std::filesystem;
using namespace coal;

static const char* VERSION = "0.3.0";
static const long long BATCH_NT = 2000000;  // same batch boundaries as pipeline.py
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
    std::string error;
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

static int run(const Options& o) {
    const GeneticCode* code = get_code(o.table);
    if (o.fmt == "codon" && (o.nogap || o.nomismatch))
        throw CoalError("\"-output codon\" cannot be combined with -nogap or -nomismatch (as in PAL2NAL)");
    const bool block = o.fmt == "clustal" || o.fmt == "codon";
    const bool two_pass = o.nogap || o.nomismatch;
    const bool staged = two_pass || block;

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
                                std::to_string(row.size()) + ", at " + id + ")");
            }
            update_column_widths(widths, row);
            ids.push_back(id);
        }
    }
    DnaIndex dna(o.nucs);
    if (ids.size() != dna.size())
        throw CoalError("number of sequences differs (protein: " + std::to_string(ids.size()) +
                        ", DNA: " + std::to_string(dna.size()) + ")");
    if (ids.empty()) throw CoalError("no sequences in input");
    std::vector<size_t> pairing;
    std::vector<std::string> per_seq;
    std::string mode = pair_ids(ids, dna.ids(), dna.pos(), pairing, per_seq);
    if (mode == "order" && !o.quiet) std::cerr << "coal: warning: " << order_pairing_warning(ids, dna.ids()) << '\n';
    dna.forget_ids();
    const long long nseq = (long long)ids.size();
    if (!block) std::vector<std::string>().swap(ids);  // IDs come from pass 2 again
    long long wtotal = 0;
    for (int w : widths) wtotal += w;

    Ctx ctx{code, widths, staged, o.fmt};
    int nthreads = o.jobs > 0 ? o.jobs : (int)std::max(1u, std::thread::hardware_concurrency());

    std::string out_dir;
    if (!o.out.empty()) {
        fs::path p = fs::absolute(fs::path(o.out)).parent_path();
        out_dir = p.string();
    }
    TempFile tmp(out_dir, ".coal_");
    std::unique_ptr<TempFile> stage, prot;
    if (staged) stage = std::make_unique<TempFile>("", ".coal_rows_");
    if (o.fmt == "codon") prot = std::make_unique<TempFile>("", ".coal_prot_");
    std::vector<char> drop(o.nomismatch ? (size_t)ncol : 0, 0), bad;
    std::vector<std::string> out_ids;  // staged FASTA/PAML: IDs in output order

    if (o.fmt == "paml" && !staged) tmp.f << paml_header(nseq, wtotal);

    // pass 2: batches of ~BATCH_NT nucleotides; thread batch k while reading batch k+1
    ProteinReader pr(o.pep);
    size_t k = 0;
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
        for (Item& it : cur)
            if (it.failed) throw CoalError(it.error);
        for (Item& it : cur) {
            if (!o.quiet)
                for (auto& w : it.warnings) std::cerr << "coal: warning: " << it.id << ' ' << w << '\n';
            if (staged) {
                stage->f.write(it.out.data(), (std::streamsize)it.out.size());
                if (prot) prot->f.write(it.row.data(), (std::streamsize)it.row.size());
                if (o.nogap) mark_nogap_triplets(it.out, bad);
                if (o.nomismatch)
                    for (int c : it.mcols) drop[(size_t)(c - 1)] = 1;
                if (!block) out_ids.push_back(it.id);
            } else {
                tmp.f.write(it.out.data(), (std::streamsize)it.out.size());
            }
        }
        std::swap(cur, nxt);
    }
    std::cerr.flush();

    if (staged) {
        std::vector<std::pair<long long, long long>> sl;
        if (two_pass) sl = keep_slices(widths, drop, bad);
        else if (wtotal > 0) sl.push_back({0, wtotal});
        long long length = 0;
        for (auto& [a, b] : sl) length += b - a;
        bool identity = sl.empty() ? wtotal == 0 : (sl.size() == 1 && sl[0].first == 0 && sl[0].second == wtotal);
        stage->f.flush();
        std::string row;
        if (!block) {
            if (o.fmt == "paml") tmp.f << paml_header(nseq, length);
            std::string text;
            stage->f.seekg(0);
            for (long long i = 0; i < nseq; ++i) {
                read_next(stage->f, wtotal, row);
                text.clear();
                record_into(text, o.fmt, out_ids[(size_t)i], identity ? row : apply_slices(row, sl));
                tmp.f << text;
            }
        } else {
            std::unique_ptr<TempFile> grid;
            std::fstream* g = &stage->f;
            if (!identity) {
                grid = std::make_unique<TempFile>("", ".coal_grid_");
                stage->f.seekg(0);
                for (long long i = 0; i < nseq; ++i) {
                    read_next(stage->f, wtotal, row);
                    std::string s = apply_slices(row, sl);
                    grid->f.write(s.data(), (std::streamsize)s.size());
                }
                grid->f.flush();
                g = &grid->f;
            }
            if (prot) prot->f.flush();
            write_blocks(tmp.f, o.fmt, *g, length, prot ? &prot->f : nullptr, ncol, ids);
        }
    }
    tmp.f.flush();
    if (!tmp.f) throw CoalError("could not write the output (disk full?)");
    tmp.close();
    if (!o.out.empty()) {
        std::error_code ec;
        fs::rename(tmp.path, o.out, ec);
        if (ec) throw CoalError("could not write " + o.out + ": " + ec.message());
        tmp.path.clear();
    } else {
#ifdef _WIN32
        _setmode(_fileno(stdout), _O_BINARY);
#endif
        std::ifstream in(tmp.path, std::ios::binary);
        std::vector<char> buf(1 << 20);
        for (;;) {
            in.read(buf.data(), (std::streamsize)buf.size());
            std::streamsize n = in.gcount();
            if (n <= 0) break;
            std::fwrite(buf.data(), 1, (size_t)n, stdout);
        }
        std::fflush(stdout);
    }
    return 0;
}

// ------------------------------------------------------------------ CLI

static const char* USAGE =
    "usage: coal-cpp [-h] [-codontable N] [-nogap] [-nomismatch] [-output {fasta,clustal,paml,codon}]\n"
    "                [-o OUT] [-j N] [-q] [--version] pep_aln nuc_fasta [nuc_fasta ...]\n";

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
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a.size() > 1 && a[0] == '-' && !(a.size() > 1 && a[1] >= '0' && a[1] <= '9')) {
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
                std::cout << USAGE;
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
            } else {
                return usage_error("unrecognized arguments: " + a);
            }
        } else {
            pos.push_back(a);
        }
    }
    if (pos.size() < 2) return usage_error("the following arguments are required: pep_aln, nuc_fasta");
    o.pep = pos[0];
    o.nucs.assign(pos.begin() + 1, pos.end());
    try {
        return run(o);
    } catch (const CoalError& e) {
        std::cerr << "coal: error: " << e.what() << '\n';
    } catch (const InputError& e) {
        std::cerr << "coal: error: " << e.what() << '\n';
    } catch (const std::exception& e) {
        std::cerr << "coal: error: " << e.what() << '\n';
    }
    return 1;
}
