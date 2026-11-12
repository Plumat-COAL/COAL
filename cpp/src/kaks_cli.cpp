// `coal kaks` (KaKs_Calculator-compatible) and `coal dnds` (tidy TSV): the command lines around kaks.hpp.
// Input: AXT (KaKs_Calculator's format: a name line, the two aligned sequences, a blank line) or aligned FASTA
// (every pair of the alignment, or consecutive pairs). Pairs are computed on a thread pool (pool.hpp) and written
// in input order. A malformed pair never stops the run: it gets status "error" (dnds) or is left out with a
// message on stderr (kaks, as KaKs_Calculator leaves invalid pairs out). docs/kaks.md has the details.
#include "kaks_cli.hpp"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include "kaks.hpp"
#include "pool.hpp"
#include "seqio.hpp"
#include "sysio.hpp"
#include "tables.hpp"

namespace coal {
namespace {

using kaks::Method;
using kaks::Result;
using kaks::Status;

const char* KAKS_HELP =
    "usage: coal kaks -i IN -o OUT [-m METHOD ...] [-c TABLE] [-j N] [--pairs all|consecutive]\n"
    "                 [--digits N] [--yn-flavour kaks|paml] [--yn-weighting 0|1]\n"
    "\n"
    "Ka/Ks for pairs of aligned coding sequences, written in KaKs_Calculator 2.0's tab-separated format\n"
    "(a drop-in for `KaKs_Calculator -i IN -o OUT -m METHOD`). Native engine; see docs/kaks.md.\n"
    "\n"
    "  -i IN          AXT (name line, two aligned sequences, blank line) or aligned FASTA\n"
    "  -o OUT         output file ('-' or none: standard output)\n"
    "  -m METHOD      NG, LWL, MLWL, LPB, MLPB, YN (default), MYN; repeat for several. KaKs_Calculator's\n"
    "                 GY, MS, MA, NONE and gamma (G*) methods are not implemented natively and are refused\n"
    "  -c TABLE       NCBI genetic code (default 1), the tables COAL's --codontable accepts\n"
    "  -j, --jobs N   threads (default: all hardware threads)\n"
    "  --pairs P      FASTA input: all (every pair, default) or consecutive (records 1+2, 3+4, ...)\n"
    "  --digits N     significant digits of numbers (default 6, as KaKs_Calculator)\n"
    "  --yn-flavour F YN00 iteration as KaKs_Calculator (kaks, default) or PAML yn00 (paml)\n"
    "  --yn-weighting W  1: weight pathways (default); 0: equal weights (PAML yn00 weighting = 0)\n"
    "  -q, --quiet    no summary on stderr\n";

const char* DNDS_HELP =
    "usage: coal dnds IN [-o OUT] [-m METHOD ...] [-c TABLE] [-j N] [--pairs all|consecutive]\n"
    "                 [--format auto|axt|fasta] [--digits N] [--yn-flavour kaks|paml] [--yn-weighting 0|1]\n"
    "\n"
    "dN/dS for pairs of aligned coding sequences as a tidy TSV, one row per pair and method:\n"
    "seq1 seq2 method dN dS dNdS SE_dN SE_dS N_sites S_sites Nd Sd codons_compared status note\n"
    "status: ok, warning (internal stop codons skipped), saturated (dN or dS NA), no_data, error.\n"
    "\n"
    "  IN             aligned FASTA or AXT file (or -i IN)\n"
    "  -o OUT         output file (default: standard output)\n"
    "  -m METHOD      YN (default), NG, LWL, MLWL, LPB, MLPB, MYN; repeat for several\n"
    "  -c, --codontable TABLE  NCBI genetic code (default 1)\n"
    "  -j, --jobs N   threads (default: all hardware threads)\n"
    "  --pairs P      FASTA input: all (every pair, default) or consecutive (records 1+2, 3+4, ...)\n"
    "  --format F     input format (default auto: '>' first = FASTA, else AXT)\n"
    "  --digits N     significant digits of numbers (default 8)\n"
    "  --yn-flavour F YN00 iteration as KaKs_Calculator (kaks, default) or PAML yn00 (paml)\n"
    "  --yn-weighting W  1: weight pathways (default); 0: equal weights (PAML yn00 weighting = 0)\n"
    "  -q, --quiet    no summary on stderr\n";

struct Cli {
    bool compat = false;
    std::string in, out;
    std::vector<Method> methods;
    int table = 1;
    int jobs = 0;
    int digits = -1;
    bool all_pairs = true;
    std::string format = "auto";
    kaks::Options opt;
    bool quiet = false;
};

int usage(const Cli& c, const std::string& msg) {
    std::cerr << (c.compat ? "usage: coal kaks -i IN -o OUT [-m METHOD ...] [-c TABLE]\n"
                           : "usage: coal dnds IN [-o OUT] [-m METHOD ...] [-c TABLE]\n")
              << "coal " << (c.compat ? "kaks" : "dnds") << ": error: " << msg << "\n";
    return 2;
}

bool to_int(const std::string& s, int& v) {
    if (s.empty()) return false;
    char* end = nullptr;
    long x = std::strtol(s.c_str(), &end, 10);
    if (*end != '\0' || x < -1000000 || x > 1000000) return false;
    v = (int)x;
    return true;
}

std::string upper(std::string s) {
    for (char& ch : s) ch = (char)std::toupper((unsigned char)ch);
    return s;
}

// ------------------------------------------------------------------------------------------------ input
struct Pair {
    std::string name, id1, id2;
    std::string own1, own2;                       // AXT: the rows themselves
    const std::string* r1 = nullptr;              // FASTA: rows in the record table
    const std::string* r2 = nullptr;
    std::string read_error;                       // malformed record (status error)
    std::vector<Result> res;
    const std::string& row1() const { return r1 ? *r1 : own1; }
    const std::string& row2() const { return r2 ? *r2 : own2; }
};

std::string strip_ws(const std::string& s) {
    std::string o;
    o.reserve(s.size());
    for (char ch : s)
        if (!py_space((unsigned char)ch)) o += ch;
    return o;
}

// AXT as KaKs_Calculator reads it: a name line, then sequence lines up to a blank line; with exactly two lines
// they are the two sequences, with three or more the concatenation is split in half (wrapped sequences). One
// line is an error for that record (KaKs_Calculator would split it in half as well).
class AxtReader {
public:
    explicit AxtReader(const std::string& path) : lr_(path) {}
    bool next(Pair& p) {
        std::string line;
        std::string name;
        bool have = false;
        while (lr_.next(line)) {
            std::string t = strip_ws(line);
            if (t.empty() || (!have && line[0] == '#')) continue;
            name = trim(line);
            have = true;
            break;
        }
        if (!have) return false;
        std::vector<std::string> rows;
        while (lr_.next(line)) {
            std::string t = strip_ws(line);
            if (t.empty()) break;
            rows.push_back(t);
        }
        p = Pair();
        p.name = name;
        const size_t amp = name.find('&');
        if (amp != std::string::npos) {
            p.id1 = name.substr(0, amp);
            p.id2 = name.substr(amp + 1);
        } else {
            p.id1 = name;
            p.id2 = "NA";
        }
        if (rows.size() == 2) {
            p.own1 = rows[0];
            p.own2 = rows[1];
        } else if (rows.size() == 1) {
            // a name and one sequence line: the second sequence is missing (a truncated file, or a blank line
            // lost between two records). Never cut the one line in half and compare the halves.
            p.read_error = "record has one sequence line: the name line must be followed by two aligned sequences";
        } else {
            std::string all;
            for (auto& r : rows) all += r;
            if (rows.empty() || all.size() % 2 != 0)
                p.read_error = "record does not hold two sequences of equal length";
            p.own1 = all.substr(0, all.size() / 2);
            p.own2 = all.substr(all.size() / 2);
        }
        return true;
    }

private:
    static std::string trim(const std::string& s) {
        size_t a = 0, b = s.size();
        while (a < b && py_space((unsigned char)s[a])) ++a;
        while (b > a && py_space((unsigned char)s[b - 1])) --b;
        return s.substr(a, b - a);
    }
    LineReader lr_;
};

struct FastaRec {
    std::string id, seq;
};

std::vector<FastaRec> read_fasta(const std::string& path) {
    LineReader lr(path);
    std::vector<FastaRec> recs;
    std::string line;
    while (lr.next(line)) {
        if (!line.empty() && line[0] == '>') {
            FastaRec r;
            size_t a = 1;
            while (a < line.size() && py_space((unsigned char)line[a])) ++a;
            size_t b = a;
            while (b < line.size() && !py_space((unsigned char)line[b])) ++b;
            r.id = line.substr(a, b - a);
            if (r.id.empty()) throw InputError("line " + std::to_string(lr.lineno()) + ": empty FASTA header");
            recs.push_back(std::move(r));
        } else {
            std::string t = strip_ws(line);
            if (t.empty()) continue;
            if (recs.empty())
                throw InputError("line " + std::to_string(lr.lineno()) + ": sequence before the first '>' header");
            recs.back().seq += t;
        }
    }
    return recs;
}

std::string sniff(const std::string& path) {
    LineReader lr(path);
    std::string line;
    while (lr.next(line)) {
        std::string t = strip_ws(line);
        if (t.empty()) continue;
        return t[0] == '>' ? "fasta" : "axt";
    }
    return "axt";
}

// ------------------------------------------------------------------------------------------------ output
std::string num(double v, int digits) {
    if (std::isnan(v)) return "NA";
    if (v == 0) return "0";
    if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*g", digits, v);
    return buf;
}

const double SMALL = 1e-6;

// KaKs_Calculator 2.0's header (KaKs.cpp), verbatim.
const char* COMPAT_HEADER =
    "Sequence\tMethod\tKa\tKs\tKa/Ks\tP-Value(Fisher)\tLength\tS-Sites\tN-Sites\tFold-Sites(0:2:4)\t"
    "Substitutions\tS-Substitutions\tN-Substitutions\tFold-S-Substitutions(0:2:4)\tFold-N-Substitutions(0:2:4)\t"
    "Divergence-Time\tSubstitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)\tGC(1:2:3)\tML-Score\tAICc\t"
    "Akaike-Weight\tModel\n";

// One row as KaKs_Calculator's parseOutput writes it: values below 1e-6 are NA (so a Ka of 0 is NA there).
std::string compat_row(const std::string& name, const Result& r, int dg) {
    auto pos = [&](double v) { return (std::isnan(v) || v < SMALL) ? std::string("NA") : num(v, dg); };
    std::string o;
    auto add = [&](const std::string& s, char end = '\t') {
        o += s;
        o += end;
    };
    add(name);
    add(kaks::method_name(r.method));
    add(pos(r.ka));
    add(pos(r.ks));
    add((std::isnan(r.ks) || r.ks < SMALL || std::isnan(r.ka)) ? "NA" : num(r.ka / r.ks, dg));
    add(num(r.p_fisher, dg));
    add(std::to_string(r.length));
    add(pos(r.S));
    add(pos(r.N));
    if (!r.has_fold || (r.fold_sites[0] < SMALL && r.fold_sites[1] < SMALL && r.fold_sites[2] < SMALL))
        add("NA");
    else
        add(num(r.fold_sites[0], dg) + ":" + num(r.fold_sites[1], dg) + ":" + num(r.fold_sites[2], dg));
    add(num(r.substitutions, dg));
    add(r.Sd > SMALL ? num(r.Sd, dg) : "NA");
    add(r.Nd > SMALL ? num(r.Nd, dg) : "NA");
    if (r.has_fold && (r.fold_s[0] != 0 || r.fold_s[1] != 0 || r.fold_s[2] != 0))
        add(num(r.fold_s[0], dg) + ":" + num(r.fold_s[1], dg) + ":" + num(r.fold_s[2], dg));
    else
        add("NA");
    if (r.has_fold && (r.fold_n[0] != 0 || r.fold_n[1] != 0 || r.fold_n[2] != 0))
        add(num(r.fold_n[0], dg) + ":" + num(r.fold_n[1], dg) + ":" + num(r.fold_n[2], dg));
    else
        add("NA");
    add(pos(r.divergence));
    {
        std::string s;
        for (int i = 0; i < 6; ++i) s += (i ? ":" : "") + num(r.rates[i], dg);
        add(s);
    }
    if (std::isnan(r.gc[0]))
        add("NA");
    else
        add(num(r.gc[0], dg) + "(" + num(r.gc[1], dg) + ":" + num(r.gc[2], dg) + ":" + num(r.gc[3], dg) + ")");
    add("NA");  // ML-Score
    add("NA");  // AICc
    add("NA");  // Akaike-Weight
    add("NA", '\n');  // Model
    return o;
}

const char* TIDY_HEADER =
    "seq1\tseq2\tmethod\tdN\tdS\tdNdS\tSE_dN\tSE_dS\tN_sites\tS_sites\tNd\tSd\tcodons_compared\tstatus\tnote\n";

std::string clean_note(const std::string& s) {
    std::string o = s;
    for (char& ch : o)
        if (ch == '\t' || ch == '\n' || ch == '\r') ch = ' ';
    return o.empty() ? "" : o;
}

std::string tidy_row(const Pair& p, const Result& r, int dg) {
    std::string o;
    auto add = [&](const std::string& s, char end = '\t') {
        o += s;
        o += end;
    };
    add(p.id1);
    add(p.id2);
    add(kaks::method_name(r.method));
    add(num(r.ka, dg));
    add(num(r.ks, dg));
    add(num(r.ka_ks, dg));
    add(num(r.se_ka, dg));
    add(num(r.se_ks, dg));
    add(num(r.N, dg));
    add(num(r.S, dg));
    add(num(r.Nd, dg));
    add(num(r.Sd, dg));
    add(std::to_string(r.codons_compared));
    add(kaks::status_name(r.status));
    add(clean_note(r.note), '\n');
    return o;
}

int parse(const std::vector<std::string>& a, Cli& c) {
    std::vector<std::string> pos;
    for (size_t i = 0; i < a.size(); ++i) {
        std::string s = a[i];
        std::string val;
        bool has_eq = false;
        if (s.size() > 2 && s[0] == '-' && s[1] == '-') {
            size_t eq = s.find('=');
            if (eq != std::string::npos) {
                val = s.substr(eq + 1);
                s = s.substr(0, eq);
                has_eq = true;
            }
        }
        // KaKs_Calculator reads its own options in any case (-I, -O, -M, -C, -D)
        const std::string u =
            c.compat && s.size() == 2 && s[0] == '-' && std::strchr("iomcdIOMCD", s[1]) ? upper(s) : s;
        auto need = [&](std::string& out) -> bool {
            if (has_eq) {
                out = val;
                return true;
            }
            if (i + 1 >= a.size()) return false;
            out = a[++i];
            return true;
        };
        std::string v;
        if (u == "-h" || u == "-H" || u == "--help") {
            std::cout << (c.compat ? KAKS_HELP : DNDS_HELP);
            return 0;
        } else if (u == "-I" || u == "-i" || u == "--in") {
            if (!need(v)) return usage(c, s + ": expected one argument");
            c.in = v;
        } else if (u == "-O" || u == "-o" || u == "--out") {
            if (!need(v)) return usage(c, s + ": expected one argument");
            c.out = v;
        } else if (u == "-M" || u == "-m" || u == "--method") {
            if (!need(v)) return usage(c, s + ": expected one argument");
            Method m;
            std::string why;
            if (!kaks::parse_method(v, m, why)) return usage(c, why);
            bool dup = false;
            for (Method x : c.methods) dup = dup || x == m;
            if (!dup) c.methods.push_back(m);
        } else if (u == "-C" || u == "-c" || u == "--codontable" || u == "-codontable") {
            if (!need(v) || !to_int(v, c.table)) return usage(c, s + ": expected an integer");
            if (!get_code(c.table)) {
                std::string list;
                for (int t : supported_tables()) list += (list.empty() ? "" : ", ") + std::to_string(t);
                return usage(c, "genetic code " + v + " is not supported (choose from " + list + ")");
            }
        } else if (u == "-j" || u == "--jobs") {
            if (!need(v) || !to_int(v, c.jobs) || c.jobs < 0) return usage(c, s + ": expected a number >= 0");
        } else if (u == "--digits") {
            if (!need(v) || !to_int(v, c.digits) || c.digits < 1 || c.digits > 17)
                return usage(c, "--digits: expected 1 to 17");
        } else if (u == "--pairs") {
            if (!need(v) || (v != "all" && v != "consecutive")) return usage(c, "--pairs: expected all or consecutive");
            c.all_pairs = v == "all";
        } else if (u == "--format") {
            if (!need(v) || (v != "auto" && v != "axt" && v != "fasta"))
                return usage(c, "--format: expected auto, axt or fasta");
            c.format = v;
        } else if (u == "--yn-flavour" || u == "--yn-flavor") {
            if (!need(v) || (v != "kaks" && v != "paml")) return usage(c, s + ": expected kaks or paml");
            c.opt.yn_flavour = v == "kaks" ? kaks::Flavour::kaks_calculator : kaks::Flavour::paml;
        } else if (u == "--yn-weighting") {
            if (!need(v) || (v != "0" && v != "1")) return usage(c, "--yn-weighting: expected 0 or 1");
            c.opt.yn_weighting = v == "1";
        } else if (u == "-q" || u == "--quiet") {
            c.quiet = true;
        } else if (u == "-D" && c.compat) {
            return usage(c, "-d (model-selection details) belongs to KaKs_Calculator's MS/MA methods, which COAL "
                            "does not implement natively");
        } else if (s == "-") {
            pos.push_back(s);
        } else if (!s.empty() && s[0] == '-') {
            return usage(c, "unrecognized argument: " + s);
        } else {
            pos.push_back(s);
        }
    }
    if (!pos.empty()) {
        if (c.compat) return usage(c, "unrecognized argument: " + pos[0] + " (the input is given with -i)");
        if (pos.size() > 1 || !c.in.empty()) return usage(c, "one input file expected");
        c.in = pos[0];
    }
    if (c.in.empty()) return usage(c, c.compat ? "-i IN is required" : "the input file is required");
    if (c.methods.empty()) c.methods.push_back(Method::YN);
    if (c.digits < 0) c.digits = c.compat ? 6 : 8;
    if (c.out == "-") c.out.clear();
    if (c.in == "-") return usage(c, "reading standard input ('-') is not supported: give a file");
    // The output is opened (and emptied) before the input is read to the end, so writing onto the input would
    // destroy it: refused, whatever the spelling of the two paths (./x, a symbolic or hard link).
    if (!c.out.empty() && same_file(c.in, c.out))
        return usage(c, "the output file " + c.out + " is the input file: choose another name for -o");
    return -1;
}

}  // namespace

const char* kaks_compat_header() { return COMPAT_HEADER; }
std::string kaks_compat_row(const std::string& name, const kaks::Result& r, int digits) {
    return compat_row(name, r, digits);
}
std::string kaks_num(double v, int digits) { return num(v, digits); }

int kaks_main(const std::string& sub, const std::vector<std::string>& argv) {
    Cli c;
    c.compat = sub == "kaks";
    {
        int rc = parse(argv, c);
        if (rc >= 0) return rc;
    }
    const GeneticCode* code = get_code(c.table);
    // methods in KaKs_Calculator's output order (it ignores the order of -m)
    std::vector<Method> methods;
    for (Method m : kaks::output_order())
        for (Method x : c.methods)
            if (x == m) methods.push_back(m);
    const char* prog = c.compat ? "coal kaks" : "coal dnds";
    try {
        std::string fmt = c.format == "auto" ? sniff(c.in) : c.format;
        std::unique_ptr<AxtReader> axt;
        std::vector<FastaRec> recs;
        size_t pi = 0, pj = 1;  // next FASTA pair
        if (fmt == "axt") {
            axt.reset(new AxtReader(c.in));
        } else {
            recs = read_fasta(c.in);
            if (recs.size() < 2) throw InputError("the FASTA input needs at least two sequences");
            if (!c.all_pairs && recs.size() % 2 != 0)
                throw InputError("--pairs consecutive needs an even number of sequences (" +
                                 std::to_string(recs.size()) + " found)");
        }
        std::ofstream fout;
        if (!c.out.empty()) {
            // (checked in parse() too, before anything is read; opening the output truncates it)
            fout.open(c.out, std::ios::binary);
            if (!fout) throw InputError("cannot write '" + c.out + "'", "file_error");
        }
#ifdef _WIN32
        if (c.out.empty()) _setmode(_fileno(stdout), _O_BINARY);  // \n line ends, as KaKs_Calculator's files
#endif
        std::ostream& os = c.out.empty() ? std::cout : fout;
        os << (c.compat ? COMPAT_HEADER : TIDY_HEADER);
        long n_pairs = 0, n_status[5] = {0, 0, 0, 0, 0};
        std::string first_error;
        auto read_batch = [&](std::vector<Pair>& b) {
            size_t nt = 0;
            while (b.size() < 256 && nt < (4u << 20)) {
                Pair p;
                if (axt) {
                    if (!axt->next(p)) break;
                } else {
                    if (c.all_pairs) {
                        if (pi >= recs.size() - 1) break;
                        p.r1 = &recs[pi].seq;
                        p.r2 = &recs[pj].seq;
                        p.id1 = recs[pi].id;
                        p.id2 = recs[pj].id;
                        if (++pj >= recs.size()) {
                            ++pi;
                            pj = pi + 1;
                        }
                    } else {
                        if (pi >= recs.size()) break;
                        p.r1 = &recs[pi].seq;
                        p.r2 = &recs[pi + 1].seq;
                        p.id1 = recs[pi].id;
                        p.id2 = recs[pi + 1].id;
                        pi += 2;
                    }
                    p.name = p.id1 + "&" + p.id2;
                }
                nt += p.row1().size() + p.row2().size();
                b.push_back(std::move(p));
            }
        };
        auto work = [&](Pair& p) {
            if (!p.read_error.empty()) {
                Result r;
                r.status = Status::error;
                r.note = p.read_error;
                for (Method m : methods) {
                    r.method = m;
                    p.res.push_back(r);
                }
                return;
            }
            kaks::Prepared prep = kaks::prepare(p.row1(), p.row2(), *code);
            for (Method m : methods) p.res.push_back(kaks::compute(prep, m, *code, c.opt));
        };
        auto emit = [&](std::vector<Pair>& b) {
            std::string buf;
            for (Pair& p : b) {
                ++n_pairs;
                for (const Result& r : p.res) {
                    ++n_status[(int)r.status];
                    if (r.status == Status::error) {
                        if (first_error.empty()) first_error = p.name + ": " + r.note;
                        if (c.compat) continue;  // KaKs_Calculator leaves invalid pairs out
                    }
                    buf += c.compat ? compat_row(p.name, r, c.digits) : tidy_row(p, r, c.digits);
                }
            }
            os << buf;
            if (!os) throw InputError("write error on " + (c.out.empty() ? std::string("stdout") : c.out), "file_error");
        };
        run_ordered<Pair>(read_batch, work, emit, c.jobs);
        os.flush();
        if (!os) throw InputError("write error", "file_error");
        if (!c.quiet) {
            std::cerr << prog << ": " << n_pairs << " pair(s), " << methods.size() << " method(s): " << n_status[0]
                      << " ok, " << n_status[1] << " warning, " << n_status[2] << " saturated, " << n_status[3]
                      << " no_data, " << n_status[4] << " error\n";
            if (n_status[4] > 0)
                std::cerr << prog << ": warning: " << n_status[4] << " result(s) could not be computed"
                          << (c.compat ? " and were left out" : " (status error)") << "; first: " << first_error
                          << "\n";
        }
        return 0;
    } catch (const InputError& e) {
        std::cerr << prog << ": error: " << e.what() << "\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << prog << ": error: " << e.what() << "\n";
        return 1;
    }
}

}  // namespace coal
