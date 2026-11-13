// `coal-cpp translate` and `coal-cpp codon-dnds`: the two native steps of `coal pipeline` (docs/pipeline.md).
//
//   translate   CDS FASTA -> protein FASTA (terminal stop dropped) + a per-gene table of what was found
//   codon-dnds  a stream of groups (pairs or clusters: member IDs, their CDS, optionally their rows of a protein
//               alignment or a DIAMOND hit) -> protein alignment (NW for pairs without rows) -> codon alignment ->
//               Ka/Ks for the pairs asked for -> one TSV row per pair and method.
//
// codon-dnds reads the groups from a file or standard input and holds one batch of them at a time, so memory is
// bounded by the batch (and by the largest single group), not by the input. Groups are prepared and pairs computed
// on threads; output is written in input order whatever -j is. A gene, pair or group that cannot be used never
// stops the run: its pairs get status "error" with the reason, and --failures lists the genes and groups.
#include "dnds_cli.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

#include "align.hpp"
#include "core.hpp"
#include "dnds.hpp"
#include "kaks.hpp"
#include "kaks_cli.hpp"
#include "pool.hpp"
#include "seqio.hpp"
#include "shortopts.hpp"
#include "sysio.hpp"
#include "tables.hpp"

namespace coal {
namespace {

using kaks::Method;
using kaks::Result;
using kaks::Status;

const char* TRANSLATE_HELP =
    "usage: coal translate CDS.fasta [CDS.fasta ...] [-o PROTEINS.fasta] [--genes GENES.tsv] [-c TABLE] [-q]\n"
    "\n"
    "Translate coding sequences from their first base (a step of `coal pipeline`, docs/pipeline.md).\n"
    "A stop codon at the end is left out; an internal stop is written as '*'; a codon with an ambiguous base is X\n"
    "unless every reading gives the same residue; bases after the last complete codon are ignored. Genes that\n"
    "cannot be translated (empty, a character that is no nucleotide code, a repeated ID, an ID with a\n"
    "non-ASCII character) are left out of the protein file and listed as failed in GENES.tsv.\n"
    "\n"
    "  -o FILE        protein FASTA, one line per sequence (default: standard output)\n"
    "  --genes FILE   per-gene table: id nt_length aa_length terminal_stop internal_stops ambiguous_codons\n"
    "                 trailing_nt status code note   (status: ok, warning, failed)\n"
    "  -c, --codontable TABLE  NCBI genetic code (default 1)\n"
    "  -q, --quiet    no summary on stderr\n";

const char* CODON_DNDS_HELP =
    "usage: coal codon-dnds [-i GROUPS] -o PAIRS.tsv [-m METHOD ...] [-c TABLE] [-j N] [--group-column NAME]\n"
    "                 [--codon-out FILE] [--codon-dir DIR] [--kaks-out FILE] [--failures FILE] [--digits N]\n"
    "                 [--matrix NAME] [--gap-open X] [--gap-extend X] [--end-gaps free|penalized] [--max-cells N]\n"
    "                 [--flanks gap|align] [--yn-flavour kaks|paml] [--yn-weighting 0|1]\n"
    "                 [--comeron-flavour coal|gestimator] [--li-flavour coal|seqinr]\n"
    "                 [--terminal-stops drop|compare] [-q]\n"
    "\n"
    "Protein alignment -> codon alignment -> Ka/Ks for groups of coding sequences (the engine of\n"
    "`coal pipeline`, docs/pipeline.md). GROUPS (a file, or standard input when -i is - or absent) is\n"
    "tab-separated, one group after the other:\n"
    "  G <group id> <all|rep|list> [nocodon]  start of a group; which pairs: all, member 1 against the others,\n"
    "                                  P lines; nocodon: do not write this group's codon alignment\n"
    "  M <id> <CDS> [<protein row>]    a member; rows (when given, for every member) are the protein alignment\n"
    "  F <id> <reason>                 a member already known to be unusable (its pairs get status error)\n"
    "  H <qstart> <qend> <sstart> <send> <qseq_gapped> <sseq_gapped>   2 members: use this DIAMOND hit\n"
    "  P <i> <j>                       a pair of member numbers (from 0), with `list`\n"
    "  X <code> <reason>               the group is already known to have failed (its pairs get status error)\n"
    "  E                               end of the group\n"
    "Two members without rows are aligned here (global, BLOSUM62, needle defaults).\n"
    "\n"
    "  -o FILE          the pair table: [NAME] seq1 seq2 method dN dS dNdS SE_dN SE_dS N_sites S_sites Nd Sd\n"
    "                   codons_compared prot_aln_length prot_identity status note\n"
    "  --group-column NAME  write the group id as a first column called NAME\n"
    "  -m METHOD        YN (default), NG, LWL, MLWL, LPB, MLPB, MYN, Comeron (alias C95), Li; repeat for\n"
    "                   several (Comeron and Li: genetic code 1 only)\n"
    "  --comeron-flavour F, --li-flavour F  as in `coal dnds` (coal, the default; gestimator; seqinr)\n"
    "  --terminal-stops S  with the flavours gestimator / seqinr: drop (default) ends a pair before the\n"
    "                   terminal stop codons; compare takes the two stop codons as the last compared column\n"
    "                   when both genes end in one (as orthologr 0.4.3 does)\n"
    "  -c, --codontable TABLE  NCBI genetic code (default 1)\n"
    "  -j, --jobs N     threads (default: all hardware threads)\n"
    "  --codon-out FILE the codon alignments, one after the other (pairs: records 1+2, 3+4, ...)\n"
    "  --codon-dir DIR  one codon alignment file per group in DIR\n"
    "  --kaks-out FILE  the same results in KaKs_Calculator's format (pairs named seq1&seq2)\n"
    "  --failures FILE  genes and groups that could not be used: group seq code message\n"
    "  --digits N       significant digits (default 8)\n"
    "  -q, --quiet      no summary on stderr\n";

int usage(const char* prog, const std::string& msg) {
    std::cerr << "coal " << prog << ": error: " << msg << " (see coal " << prog << " --help)\n";
    return 2;
}

bool to_int(const std::string& s, long long& v) {
    if (s.empty()) return false;
    char* end = nullptr;
    v = std::strtoll(s.c_str(), &end, 10);
    return *end == '\0';
}

bool to_double(const std::string& s, double& v) {
    if (s.empty()) return false;
    char* end = nullptr;
    v = std::strtod(s.c_str(), &end);
    return *end == '\0' && std::isfinite(v);
}

std::string clean(const std::string& s) {
    std::string o = s;
    for (char& ch : o)
        if (ch == '\t' || ch == '\n' || ch == '\r') ch = ' ';
    return o;
}

struct ArgIter {
    std::vector<std::string> a;
    std::string value_opts, flag_opts;   // one-letter options with / without a value (-j4, -qj4: shortopts.hpp)
    size_t i = 0;
    std::string opt, val;
    bool has_eq = false;
    ArgIter(const std::vector<std::string>& args, const char* with_value, const char* flags)
        : a(args), value_opts(with_value), flag_opts(flags) {}
    bool next() {
        if (i >= a.size()) return false;
        expand_short_option(a, i, value_opts, flag_opts, {"-codontable"});
        opt = a[i++];
        has_eq = false;
        if (opt.size() > 2 && opt[0] == '-' && opt[1] == '-') {
            size_t eq = opt.find('=');
            if (eq != std::string::npos) {
                val = opt.substr(eq + 1);
                opt = opt.substr(0, eq);
                has_eq = true;
            }
        }
        return true;
    }
    bool value(std::string& out) {
        if (has_eq) {
            out = val;
            return true;
        }
        if (i >= a.size()) return false;
        out = a[i++];
        return true;
    }
};

void binary_stdio(bool in, bool out) {
#ifdef _WIN32
    if (in) _setmode(_fileno(stdin), _O_BINARY);
    if (out) _setmode(_fileno(stdout), _O_BINARY);
#else
    (void)in;
    (void)out;
#endif
}

// A sequence ID with non-ASCII bytes, for the genes table: unchanged when it is valid UTF-8 (so that it can be
// found again in the input), otherwise every byte that is not part of a valid UTF-8 sequence as \xhh (two lower-case
// hexadecimal digits). The tables are therefore always valid UTF-8, two different IDs never look the same, and the
// text is what Python's "backslashreplace" error handler gives for the same bytes (the pipeline reads the IDs of
// the pairs and clusters files with it, so both name a gene the same way). "Valid" is the Unicode definition:
// no overlong forms, no surrogates (ED A0..BF), nothing above U+10FFFF.
std::string printable_id(const std::string& id) {
    static const char* hex = "0123456789abcdef";
    std::string o;
    o.reserve(id.size());
    for (size_t i = 0; i < id.size();) {
        const unsigned char c = (unsigned char)id[i];
        size_t n = c < 0x80 ? 1 : (c >= 0xC2 && c <= 0xDF) ? 2 : (c >= 0xE0 && c <= 0xEF) ? 3
                   : (c >= 0xF0 && c <= 0xF4) ? 4 : 0;
        if (n > 1) {
            if (i + n > id.size()) n = 0;
            for (size_t k = 1; k < n; ++k)
                if (((unsigned char)id[i + k] & 0xC0) != 0x80) n = 0;
        }
        if (n > 1) {
            const unsigned char d = (unsigned char)id[i + 1];
            if ((c == 0xE0 && d < 0xA0) || (c == 0xED && d > 0x9F) || (c == 0xF0 && d < 0x90) || (c == 0xF4 && d > 0x8F))
                n = 0;
        }
        if (n == 0) {
            o += "\\x";
            o += hex[c >> 4];
            o += hex[c & 15];
            ++i;
        } else {
            o.append(id, i, n);
            i += n;
        }
    }
    return o;
}

// ------------------------------------------------------------------------------------------------ translate
int translate_main_(const std::vector<std::string>& argv) {
    const char* prog = "translate";
    std::vector<std::string> ins;
    std::string out, genes;
    int table = 1;
    bool quiet = false;
    ArgIter it(argv, "oc", "q");
    while (it.next()) {
        const std::string& s = it.opt;
        std::string v;
        if (s == "-h" || s == "--help") {
            std::cout << TRANSLATE_HELP;
            return 0;
        } else if (s == "-o" || s == "--out") {
            if (!it.value(out)) return usage(prog, s + ": expected one argument");
        } else if (s == "--genes") {
            if (!it.value(genes)) return usage(prog, s + ": expected one argument");
        } else if (s == "-c" || s == "--codontable" || s == "-codontable") {
            long long t;
            if (!it.value(v) || !to_int(v, t) || !get_code((int)t))
                return usage(prog, "genetic code " + v + " is not supported");
            table = (int)t;
        } else if (s == "-q" || s == "--quiet") {
            quiet = true;
        } else if (s.size() > 1 && s[0] == '-') {
            return usage(prog, "unrecognized argument: " + s);
        } else {
            ins.push_back(s);
        }
    }
    if (ins.empty()) return usage(prog, "the CDS FASTA file is required");
    if (out == "-") out.clear();
    for (const std::string& in : ins)  // an output opened onto an input would empty it before it is read
        for (const std::string* o : {&out, &genes})
            if (!o->empty() && same_file(in, *o))
                return usage(prog, "the output file " + *o + " is the input file " + in + ": choose another name");
    const GeneticCode* code = get_code(table);
    try {
        std::ofstream fout, fgenes;
        if (!out.empty()) {
            fout.open(out, std::ios::binary);
            if (!fout) throw InputError("cannot write '" + out + "'", "file_error");
        } else {
            binary_stdio(false, true);
        }
        if (!genes.empty()) {
            fgenes.open(genes, std::ios::binary);
            if (!fgenes) throw InputError("cannot write '" + genes + "'", "file_error");
            fgenes << "id\tnt_length\taa_length\tterminal_stop\tinternal_stops\tambiguous_codons\ttrailing_nt\t"
                      "status\tcode\tnote\n";
        }
        std::ostream& os = out.empty() ? std::cout : fout;
        std::unordered_set<std::string> seen;
        long n = 0, n_ok = 0, n_warn = 0, n_fail = 0;
        std::string id, seq, line, obuf, gbuf;
        long long id_line = 0;
        bool have = false;
        auto flush_rec = [&]() {
            if (!have) return;
            ++n;
            dnds::Translation t;
            std::string status = "ok", codestr, note;
            bool ascii_id = true;
            for (unsigned char ch : id) ascii_id = ascii_id && ch < 128;
            if (id.empty()) {
                t.error_code = "empty_id";
                t.error = "the FASTA header has no ID";
            } else if (!ascii_id) {
                // one gene's failure, not the end of the run: the other genes of the file are translated
                t.error_code = "id_not_ascii";
                t.error = "non-ASCII character in the sequence ID (line " + std::to_string(id_line) +
                          "): use only ASCII letters, digits and punctuation in IDs";
                id = printable_id(id);
            } else if (!seen.insert(id).second) {
                t.error_code = "duplicate_id";
                t.error = "the ID occurs more than once in the file (the first record is used)";
            } else {
                t = dnds::translate_cds(seq, *code);
            }
            if (!t.error.empty()) {
                status = "failed";
                codestr = t.error_code;
                note = t.error;
                ++n_fail;
            } else {
                note = dnds::translation_note(t);
                if (!note.empty()) {
                    status = "warning";
                    codestr = t.trailing_nt ? (t.internal_stops ? "cds_length;internal_stop" : "cds_length")
                                            : "internal_stop";
                    ++n_warn;
                } else {
                    ++n_ok;
                }
                obuf += '>';
                obuf += id;
                obuf += '\n';
                obuf += t.protein;
                obuf += '\n';
            }
            if (fgenes.is_open()) {
                gbuf += id.empty() ? "(no ID)" : id;
                gbuf += '\t' + std::to_string(t.nt_length) + '\t' + std::to_string(t.protein.size()) + '\t' +
                        (t.terminal_stop ? "1" : "0") + '\t' + std::to_string(t.internal_stops) + '\t' +
                        std::to_string(t.ambiguous_codons) + '\t' + std::to_string(t.trailing_nt) + '\t' + status +
                        '\t' + codestr + '\t' + clean(note) + '\n';
            }
            if (obuf.size() > (1u << 20)) {
                os << obuf;
                obuf.clear();
            }
            if (gbuf.size() > (1u << 20)) {
                fgenes << gbuf;
                gbuf.clear();
            }
        };
        for (const std::string& in : ins) {
            LineReader lr(in, true);
            have = false;
            const long before = n;
            while (lr.next(line)) {
                if (!line.empty() && line[0] == '>') {
                    flush_rec();
                    size_t a = 1;
                    while (a < line.size() && py_space((unsigned char)line[a])) ++a;
                    size_t b = a;
                    while (b < line.size() && !py_space((unsigned char)line[b])) ++b;
                    id = line.substr(a, b - a);
                    id_line = lr.lineno();
                    seq.clear();
                    have = true;
                } else if (have) {
                    seq += line;
                } else {
                    for (unsigned char ch : line)
                        if (!py_space(ch))
                            throw InputError(in + ": line " + std::to_string(lr.lineno()) +
                                             ": sequence before the first '>' header (is this a FASTA file?)");
                }
            }
            flush_rec();
            if (n == before) throw InputError(in + ": no FASTA record found");
        }
        os << obuf;
        os.flush();
        if (!os) throw InputError("write error on " + (out.empty() ? std::string("stdout") : out), "file_error");
        if (fgenes.is_open()) {
            fgenes << gbuf;
            fgenes.flush();
            if (!fgenes) throw InputError("write error on " + genes, "file_error");
        }
        if (!quiet)
            std::cerr << "coal translate: " << n << " gene(s): " << n_ok << " ok, " << n_warn << " warning, "
                      << n_fail << " failed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "coal translate: error: " << e.what() << "\n";
        return 1;
    }
}

// ------------------------------------------------------------------------------------------------ codon-dnds
struct Cli {
    std::string in, out, group_column, codon_out, codon_dir, kaks_out, failures;
    dnds::Options opt;
    int jobs = 0;
    int digits = 8;
    bool quiet = false;
};

struct Item {
    dnds::Group g;
    dnds::GroupResult r;
    std::vector<std::pair<int, int>> pairs;
    std::string fail_code, fail;  // an X line: the group failed before it got here
    bool no_codon = false;        // G ... nocodon: leave this group out of --codon-out / --codon-dir
};

struct Task {
    Item* item = nullptr;
    int i = 0, j = 0;
    dnds::PairResult res;
};

std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> f;
    size_t a = 0;
    for (;;) {
        size_t b = line.find('\t', a);
        if (b == std::string::npos) {
            f.push_back(line.substr(a));
            break;
        }
        f.push_back(line.substr(a, b - a));
        a = b + 1;
    }
    return f;
}

class GroupReader {
public:
    explicit GroupReader(std::istream& in, std::string name) : in_(in), name_(std::move(name)) {}
    // false at the end of the input. Throws InputError for a malformed stream.
    bool next(Item& it) {
        std::string line;
        bool open = false;
        while (std::getline(in_, line)) {
            ++lineno_;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            std::vector<std::string> f = split_tabs(line);
            const std::string& k = f[0];
            if (k == "G") {
                if (open) bad("G before the previous group's E");
                if (f.size() < 3) bad("G needs a group id and all|rep|list");
                it = Item();
                it.g.id = f[1];
                if (f[2] == "all")
                    it.g.pairs = dnds::Pairs::all;
                else if (f[2] == "rep")
                    it.g.pairs = dnds::Pairs::representative;
                else if (f[2] == "list")
                    it.g.pairs = dnds::Pairs::list;
                else
                    bad("G: the third field must be all, rep or list");
                it.no_codon = f.size() > 3 && f[3] == "nocodon";
                open = true;
            } else if (!open) {
                bad("a line before the first G");
            } else if (k == "M") {
                if (f.size() < 3) bad("M needs an id and a CDS");
                dnds::Member m;
                m.id = f[1];
                m.cds = std::move(f[2]);
                if (f.size() > 3) m.row = std::move(f[3]);
                it.g.members.push_back(std::move(m));
            } else if (k == "F") {
                if (f.size() < 3) bad("F needs an id and a reason");
                dnds::Member m;
                m.id = f[1];
                m.failed = f[2].empty() ? "not usable" : f[2];
                it.g.members.push_back(std::move(m));
            } else if (k == "H") {
                long long v[4];
                if (f.size() < 7 || !to_int(f[1], v[0]) || !to_int(f[2], v[1]) || !to_int(f[3], v[2]) ||
                    !to_int(f[4], v[3]))
                    bad("H needs qstart qend sstart send qseq_gapped sseq_gapped");
                it.g.has_hit = true;
                it.g.hit.qstart = v[0];
                it.g.hit.qend = v[1];
                it.g.hit.sstart = v[2];
                it.g.hit.send = v[3];
                it.g.hit.qseq_gapped = std::move(f[5]);
                it.g.hit.sseq_gapped = std::move(f[6]);
            } else if (k == "P") {
                long long a, b;
                if (f.size() < 3 || !to_int(f[1], a) || !to_int(f[2], b) || a < 0 || b < 0 || a > 2000000000 ||
                    b > 2000000000)
                    bad("P needs two member numbers");
                it.g.pair_list.emplace_back((int)a, (int)b);
            } else if (k == "X") {
                if (f.size() < 3) bad("X needs a code and a reason");
                it.fail_code = f[1].empty() ? "failed_before" : f[1];
                it.fail = f[2].empty() ? "not usable" : f[2];
            } else if (k == "E") {
                return true;
            } else {
                bad("unknown line type '" + k.substr(0, 20) + "'");
            }
        }
        if (in_.bad()) throw InputError(name_ + ": read error", "file_error");
        if (open) bad("the input ends inside a group (no E line)");
        return false;
    }

private:
    [[noreturn]] void bad(const std::string& why) {
        throw InputError(name_ + ": line " + std::to_string(lineno_) + ": " + why);
    }
    std::istream& in_;
    std::string name_;
    long long lineno_ = 0;
};

std::string safe_name(const std::string& id) {
    std::string o;
    for (unsigned char ch : id) {
        const bool ok = (ch >= '0' && ch <= '9') || (ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') ||
                        ch == '.' || ch == '_' || ch == '-';
        o += ok ? (char)ch : '_';
    }
    if (o.empty() || o == "." || o == "..") o = "_" + o;
    if (o.size() > 150) o.resize(150);
    return o;
}

int codon_dnds_main_(const std::vector<std::string>& argv) {
    const char* prog = "codon-dnds";
    Cli c;
    std::vector<Method> methods;
    ArgIter it(argv, "iomcj", "q");
    while (it.next()) {
        const std::string& s = it.opt;
        std::string v;
        long long n;
        double d;
        if (s == "-h" || s == "--help") {
            std::cout << CODON_DNDS_HELP;
            return 0;
        } else if (s == "-i" || s == "--in") {
            if (!it.value(c.in)) return usage(prog, s + ": expected one argument");
        } else if (s == "-o" || s == "--out") {
            if (!it.value(c.out)) return usage(prog, s + ": expected one argument");
        } else if (s == "--group-column") {
            if (!it.value(c.group_column)) return usage(prog, s + ": expected one argument");
        } else if (s == "--codon-out") {
            if (!it.value(c.codon_out)) return usage(prog, s + ": expected one argument");
        } else if (s == "--codon-dir") {
            if (!it.value(c.codon_dir)) return usage(prog, s + ": expected one argument");
        } else if (s == "--kaks-out") {
            if (!it.value(c.kaks_out)) return usage(prog, s + ": expected one argument");
        } else if (s == "--failures") {
            if (!it.value(c.failures)) return usage(prog, s + ": expected one argument");
        } else if (s == "-m" || s == "--method") {
            if (!it.value(v)) return usage(prog, s + ": expected one argument");
            Method m;
            std::string why;
            if (!kaks::parse_method(v, m, why)) return usage(prog, why);
            methods.push_back(m);
        } else if (s == "-c" || s == "--codontable" || s == "-codontable") {
            if (!it.value(v) || !to_int(v, n) || !get_code((int)n))
                return usage(prog, "genetic code " + v + " is not supported");
            c.opt.table = (int)n;
        } else if (s == "-j" || s == "--jobs") {
            if (!it.value(v) || !to_int(v, n) || n < 0 || n > 4096) return usage(prog, s + ": expected a number >= 0");
            c.jobs = (int)n;
        } else if (s == "--digits") {
            if (!it.value(v) || !to_int(v, n) || n < 1 || n > 17) return usage(prog, "--digits: expected 1 to 17");
            c.digits = (int)n;
        } else if (s == "--matrix") {
            if (!it.value(c.opt.align.matrix)) return usage(prog, s + ": expected one argument");
        } else if (s == "--gap-open") {
            if (!it.value(v) || !to_double(v, d)) return usage(prog, s + ": expected a number");
            c.opt.align.gap_open = d;
        } else if (s == "--gap-extend") {
            if (!it.value(v) || !to_double(v, d)) return usage(prog, s + ": expected a number");
            c.opt.align.gap_extend = d;
        } else if (s == "--end-gaps") {
            if (!it.value(v) || (v != "free" && v != "penalized" && v != "penalised"))
                return usage(prog, "--end-gaps: expected free or penalized");
            c.opt.align.free_end_gaps = v == "free";
        } else if (s == "--max-cells") {
            if (!it.value(v) || !to_int(v, n) || n < 0) return usage(prog, s + ": expected a number >= 0");
            c.opt.align.max_cells = (uint64_t)n;
        } else if (s == "--flanks") {
            if (!it.value(v) || (v != "gap" && v != "align")) return usage(prog, "--flanks: expected gap or align");
            c.opt.flanks = v == "gap" ? Flanks::Gap : Flanks::Align;
        } else if (s == "--yn-flavour" || s == "--yn-flavor") {
            if (!it.value(v) || (v != "kaks" && v != "paml")) return usage(prog, s + ": expected kaks or paml");
            c.opt.kaks.yn_flavour = v == "kaks" ? kaks::Flavour::kaks_calculator : kaks::Flavour::paml;
        } else if (s == "--yn-weighting") {
            if (!it.value(v) || (v != "0" && v != "1")) return usage(prog, "--yn-weighting: expected 0 or 1");
            c.opt.kaks.yn_weighting = v == "1";
        } else if (s == "--comeron-flavour" || s == "--comeron-flavor") {
            if (!it.value(v) || (v != "coal" && v != "gestimator"))
                return usage(prog, s + ": expected coal or gestimator");
            c.opt.kaks.comeron_flavour = v == "coal" ? kaks::ComeronFlavour::coal : kaks::ComeronFlavour::gestimator;
        } else if (s == "--li-flavour" || s == "--li-flavor") {
            if (!it.value(v) || (v != "coal" && v != "seqinr")) return usage(prog, s + ": expected coal or seqinr");
            c.opt.kaks.li_flavour = v == "coal" ? kaks::LiFlavour::coal : kaks::LiFlavour::seqinr;
        } else if (s == "--terminal-stops") {
            if (!it.value(v) || (v != "drop" && v != "compare")) return usage(prog, s + ": expected drop or compare");
            c.opt.compare_terminal_stops = v == "compare";
        } else if (s == "-q" || s == "--quiet") {
            c.quiet = true;
        } else {
            return usage(prog, "unrecognized argument: " + s);
        }
    }
    if (c.out.empty()) return usage(prog, "-o PAIRS.tsv is required");
    if (!methods.empty()) c.opt.methods = methods;
    c.opt.methods = dnds::ordered_methods(c.opt.methods);
    for (Method m : c.opt.methods)
        if ((m == Method::Comeron || m == Method::Li) && c.opt.table != 1)
            return usage(prog, std::string("method ") + kaks::method_name(m) +
                                   " is implemented for the standard genetic code only (-c 1), as its reference");
    try {
        check_align_params(c.opt.align);
    } catch (const CoalError& e) {
        return usage(prog, e.what());
    }
    const std::vector<Method>& ms = c.opt.methods;
    const bool from_stdin = c.in.empty() || c.in == "-";
    if (!from_stdin)
        for (const std::string* o : {&c.out, &c.codon_out, &c.kaks_out, &c.failures})
            if (!o->empty() && same_file(c.in, *o))
                return usage(prog, "the output file " + *o + " is the input file: choose another name");
    try {
        std::ifstream fin;
        if (from_stdin) {
            binary_stdio(true, false);
        } else {
            fin.open(c.in, std::ios::binary);
            if (!fin) throw InputError(open_error(c.in), "file_error");
        }
        GroupReader reader(from_stdin ? std::cin : fin, from_stdin ? "<stdin>" : c.in);
        auto open_out = [](std::ofstream& f, const std::string& path) {
            if (path.empty()) return;
            f.open(path, std::ios::binary);
            if (!f) throw InputError("cannot write '" + path + "'", "file_error");
        };
        std::ofstream fout, fcodon, fkaks, ffail;
        open_out(fout, c.out);
        open_out(fcodon, c.codon_out);
        open_out(fkaks, c.kaks_out);
        open_out(ffail, c.failures);
        if (!c.group_column.empty()) fout << clean(c.group_column) << '\t';
        fout << "seq1\tseq2\tmethod\tdN\tdS\tdNdS\tSE_dN\tSE_dS\tN_sites\tS_sites\tNd\tSd\tcodons_compared\t"
                "prot_aln_length\tprot_identity\tstatus\tnote\n";
        if (fkaks.is_open()) fkaks << kaks_compat_header();
        if (ffail.is_open()) ffail << "group\tseq\tcode\tmessage\n";
        std::unordered_map<std::string, int> used_names;  // --codon-dir file names

        long long n_groups = 0, n_pairs = 0, n_status[kaks::STATUS_COUNT] = {}, n_member_fail = 0, n_group_fail = 0;
        const int nthreads = resolve_jobs(c.jobs);
        const size_t BATCH_GROUPS = 1024, BATCH_NT = 16u << 20, BATCH_PAIRS = 32768, CHUNK = 8192;

        auto prepare = [&](Item& item) {
            if (!item.fail.empty()) {
                item.r.members.resize(item.g.members.size());
                item.r.protein_rows.assign(item.g.members.size(), "");
                item.r.error = item.fail;
                item.r.error_code = item.fail_code;
            } else {
                item.r = dnds::prepare_group(item.g, c.opt);
            }
            item.pairs = dnds::group_pairs(item.g);
        };
        auto compute = [&](Task& t) { t.res = dnds::compute_pair(t.item->r, t.i, t.j, c.opt); };

        auto write_tasks = [&](std::vector<Task>& tasks) {
            std::string buf, kbuf;
            for (Task& t : tasks) {
                const dnds::Group& g = t.item->g;
                const std::string& id1 = g.members[(size_t)t.i].id;
                const std::string& id2 = g.members[(size_t)t.j].id;
                ++n_pairs;
                auto lead = [&](Method m) {
                    if (!c.group_column.empty()) {
                        buf += g.id;
                        buf += '\t';
                    }
                    buf += id1;
                    buf += '\t';
                    buf += id2;
                    buf += '\t';
                    buf += kaks::method_name(m);
                    buf += '\t';
                };
                if (!t.res.error.empty()) {
                    for (Method m : ms) {
                        lead(m);
                        buf += "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0\tNA\tNA\terror\t";
                        buf += clean(t.res.error);
                        buf += '\n';
                        ++n_status[(int)Status::error];
                    }
                    continue;
                }
                char ident[32];
                std::snprintf(ident, sizeof ident, "%.2f",
                              t.res.prot_aln_length
                                  ? 100.0 * (double)t.res.prot_identities / (double)t.res.prot_aln_length
                                  : 0.0);
                for (const Result& r : t.res.results) {
                    Status st = r.status;
                    std::string note = r.note;
                    if (!t.res.note.empty()) {
                        note += (note.empty() ? "" : "; ") + t.res.note;
                        if (st == Status::ok) st = Status::warning;
                    }
                    ++n_status[(int)st];
                    lead(r.method);
                    const int dg = c.digits;
                    buf += kaks_num(r.ka, dg) + '\t' + kaks_num(r.ks, dg) + '\t' + kaks_num(r.ka_ks, dg) + '\t' +
                           kaks_num(r.se_ka, dg) + '\t' + kaks_num(r.se_ks, dg) + '\t' + kaks_num(r.N, dg) + '\t' +
                           kaks_num(r.S, dg) + '\t' + kaks_num(r.Nd, dg) + '\t' + kaks_num(r.Sd, dg) + '\t' +
                           std::to_string(r.codons_compared) + '\t' + std::to_string(t.res.prot_aln_length) + '\t' +
                           ident + '\t' + kaks::status_name(st) + '\t' + clean(note) + '\n';
                    if (fkaks.is_open() && r.status != Status::error) kbuf += kaks_compat_row(id1 + "&" + id2, r, 6);
                }
            }
            fout << buf;
            if (!fout) throw InputError("write error on " + c.out, "file_error");
            if (fkaks.is_open()) {
                fkaks << kbuf;
                if (!fkaks) throw InputError("write error on " + c.kaks_out, "file_error");
            }
        };

        using Clock = std::chrono::steady_clock;
        double t_align = 0, t_kaks = 0;  // wall seconds: translate + align + thread codons; Ka/Ks
        auto since = [](Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); };
        std::vector<Item> batch;
        Item pending;
        bool have_pending = false, eof = false;
        while (!eof || have_pending) {
            batch.clear();
            size_t nt = 0, np = 0;
            while (batch.size() < BATCH_GROUPS && nt < BATCH_NT && np < BATCH_PAIRS) {
                if (!have_pending) {
                    if (eof || !reader.next(pending)) {
                        eof = true;
                        break;
                    }
                    have_pending = true;
                }
                size_t gnt = 0;
                for (auto& m : pending.g.members) gnt += m.cds.size();
                const size_t k = pending.g.members.size();
                const size_t gp = pending.g.pairs == dnds::Pairs::all ? k * (k ? k - 1 : 0) / 2
                                  : pending.g.pairs == dnds::Pairs::representative ? (k ? k - 1 : 0)
                                                                                    : pending.g.pair_list.size();
                nt += gnt;
                np += gp;
                batch.push_back(std::move(pending));
                have_pending = false;
            }
            if (batch.empty()) break;
            {
                const auto t0 = Clock::now();
                parallel_batch(batch, prepare, nthreads);
                t_align += since(t0);
            }
            // pairs, in order, a chunk at a time (a large group's pairs are not all held at once)
            std::vector<Task> tasks;
            tasks.reserve(CHUNK);
            for (Item& item : batch) {
                ++n_groups;
                for (auto& ij : item.pairs) {
                    Task t;
                    t.item = &item;
                    t.i = ij.first;
                    t.j = ij.second;
                    tasks.push_back(std::move(t));
                    if (tasks.size() >= CHUNK) {
                        const auto t0 = Clock::now();
                        parallel_batch(tasks, compute, nthreads);
                        t_kaks += since(t0);
                        write_tasks(tasks);
                        tasks.clear();
                    }
                }
            }
            if (!tasks.empty()) {
                const auto t0 = Clock::now();
                parallel_batch(tasks, compute, nthreads);
                t_kaks += since(t0);
                write_tasks(tasks);
                tasks.clear();
            }
            // codon alignments and failures
            for (Item& item : batch) {
                const dnds::Group& g = item.g;
                size_t n_rows = 0;
                for (auto& m : item.r.members) n_rows += m.codon_row.empty() ? 0 : 1;
                if (!item.r.error.empty()) {
                    ++n_group_fail;
                    if (ffail.is_open() && item.fail.empty())
                        ffail << g.id << "\t\t" << item.r.error_code << '\t' << clean(item.r.error) << '\n';
                }
                for (size_t k = 0; k < g.members.size(); ++k) {
                    const dnds::MemberResult& m = item.r.members[k];
                    if (m.error.empty() || m.error_code == "failed_before") continue;
                    ++n_member_fail;
                    if (ffail.is_open())
                        ffail << g.id << '\t' << g.members[k].id << '\t' << m.error_code << '\t' << clean(m.error)
                              << '\n';
                }
                const bool whole = !item.no_codon && item.r.error.empty() && n_rows >= 2 &&
                                   (g.members.size() != 2 || n_rows == 2);
                if (fcodon.is_open() && whole) {
                    std::string buf;
                    for (size_t k = 0; k < g.members.size(); ++k) {
                        if (item.r.members[k].codon_row.empty()) continue;
                        buf += '>' + g.members[k].id + '\n' + item.r.members[k].codon_row + '\n';
                    }
                    fcodon << buf;
                    if (!fcodon) throw InputError("write error on " + c.codon_out, "file_error");
                }
                if (!c.codon_dir.empty() && whole) {
                    std::string name = safe_name(g.id);
                    int& seen = used_names[name];
                    if (seen++ > 0) name += "~" + std::to_string(seen);
                    const std::string path = c.codon_dir + "/" + name + ".fasta";
                    std::ofstream f(path, std::ios::binary);
                    if (!f) throw InputError("cannot write '" + path + "'", "file_error");
                    for (size_t k = 0; k < g.members.size(); ++k) {
                        if (item.r.members[k].codon_row.empty()) continue;
                        f << '>' << g.members[k].id << '\n' << item.r.members[k].codon_row << '\n';
                    }
                    f.flush();
                    if (!f) throw InputError("write error on " + path, "file_error");
                }
            }
        }
        for (std::ofstream* f : {&fout, &fcodon, &fkaks, &ffail}) {
            if (!f->is_open()) continue;
            f->flush();
            if (!*f) throw InputError("write error", "file_error");
        }
        if (!c.quiet)
            std::cerr << "coal codon-dnds: " << n_groups << " group(s), " << n_pairs << " pair(s), " << ms.size()
                      << " method(s): " << n_status[0] << " ok, " << n_status[1] << " warning, " << n_status[2]
                      << " not_converged, " << n_status[3] << " saturated, " << n_status[4] << " no_data, "
                      << n_status[5] << " error; " << n_member_fail
                      << " gene(s) and " << n_group_fail << " group(s) failed\n"
                      << "coal codon-dnds: timing: align+codon " << t_align << " s, kaks " << t_kaks << " s, "
                      << nthreads << " thread(s)\n";
        return 0;
    } catch (const std::bad_alloc&) {
        std::cerr << "coal codon-dnds: error: out of memory (use fewer threads with -j, or smaller groups)\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "coal codon-dnds: error: " << e.what() << "\n";
        return 1;
    }
}

}  // namespace

int translate_main(const std::vector<std::string>& argv) { return translate_main_(argv); }
int codon_dnds_main(const std::vector<std::string>& argv) { return codon_dnds_main_(argv); }

}  // namespace coal
