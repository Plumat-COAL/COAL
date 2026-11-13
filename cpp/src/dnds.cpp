// The one-step dN/dS core (dnds.hpp): translate -> align -> thread codons -> Ka/Ks for pairs and groups.
#include "dnds.hpp"

#include <cctype>
#include <exception>

#include "core.hpp"

namespace coal {
namespace dnds {

Translation translate_cds(const std::string& raw, const GeneticCode& code) {
    Translation t;
    std::string& s = t.cds;
    s.reserve(raw.size());
    const uint8_t* iu = iupac_table();
    for (unsigned char ch : raw) {
        if (ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || ch == '\f' || ch == '\v') continue;
        unsigned char u = (unsigned char)std::toupper(ch);
        if (u == 'U') u = 'T';
        if (iu[u] == 255) {
            t.error_code = "cds_invalid_character";
            t.error = "the CDS holds a character that is no nucleotide code: ";
            if (ch >= 33 && ch < 127)
                t.error += std::string("'") + (char)ch + "'";
            else
                t.error += "byte " + std::to_string((int)ch);
            t.error += " at position " + std::to_string(s.size() + 1);
            s.clear();
            return t;
        }
        s += (char)ch;  // kept as written (case, U): the codon alignment copies the DNA verbatim
    }
    t.nt_length = (long)s.size();
    if (s.empty()) {
        t.error_code = "cds_empty";
        t.error = "the CDS is empty";
        return t;
    }
    if (s.size() < 3) {
        t.error_code = "cds_too_short";
        t.error = "the CDS is shorter than one codon (" + std::to_string(s.size()) + " nt)";
        s.clear();
        return t;
    }
    t.trailing_nt = (int)(s.size() % 3);
    s.resize(s.size() - (size_t)t.trailing_nt);
    const size_t n = s.size() / 3;
    const uint64_t* masks = iupac_triplet_table();
    t.protein.reserve(n);
    for (size_t k = 0; k < n; ++k) {
        auto norm = [](char ch) {
            const unsigned char u = (unsigned char)std::toupper((unsigned char)ch);
            return u == 'U' ? (unsigned char)'T' : u;
        };
        const unsigned char a = norm(s[3 * k]), b = norm(s[3 * k + 1]), c = norm(s[3 * k + 2]);
        const int ia = base_index(a), ib = base_index(b), ic = base_index(c);
        char aa;
        if (ia >= 0 && ib >= 0 && ic >= 0) {
            aa = code.aa[ia * 16 + ib * 4 + ic];
        } else {
            ++t.ambiguous_codons;
            uint64_t m = masks[(int)iu[a] * 225 + (int)iu[b] * 15 + (int)iu[c]];
            aa = 0;
            for (int idx = 0; idx < 64 && m; ++idx, m >>= 1) {
                if (!(m & 1)) continue;
                if (aa == 0)
                    aa = code.aa[idx];
                else if (aa != code.aa[idx]) {
                    aa = 'X';
                    break;
                }
            }
            if (aa == 0 || aa == '*') aa = 'X';  // an ambiguous codon is never read as a stop
        }
        if (aa == '*') {
            if (k + 1 == n) {
                t.terminal_stop = true;
                continue;
            }
            ++t.internal_stops;
        }
        t.protein += aa;
    }
    if (t.protein.empty()) {
        t.error_code = "cds_only_stop";
        t.error = "the CDS holds no codon besides a stop codon";
        s.clear();
    }
    return t;
}

std::string translation_note(const Translation& t) {
    std::string o;
    if (t.trailing_nt)
        o += "CDS length " + std::to_string(t.nt_length) + " is not a multiple of 3 (last " +
             std::to_string(t.trailing_nt) + " nt ignored)";
    if (t.internal_stops) {
        if (!o.empty()) o += "; ";
        o += std::to_string(t.internal_stops) + " internal stop codon" + (t.internal_stops == 1 ? "" : "s") +
             " in the CDS";
    }
    return o;
}

std::vector<kaks::Method> ordered_methods(const std::vector<kaks::Method>& methods) {
    std::vector<kaks::Method> out;
    for (kaks::Method m : kaks::output_order())
        for (kaks::Method x : methods)
            if (x == m) {
                out.push_back(m);
                break;
            }
    if (out.empty()) out.push_back(kaks::Method::YN);
    return out;
}

std::vector<std::pair<int, int>> group_pairs(const Group& g) {
    std::vector<std::pair<int, int>> out;
    const int n = (int)g.members.size();
    if (g.pairs == Pairs::list) {
        for (auto& p : g.pair_list)
            if (p.first >= 0 && p.second >= 0 && p.first < n && p.second < n && p.first != p.second) out.push_back(p);
    } else if (g.pairs == Pairs::representative) {
        for (int j = 1; j < n; ++j) out.emplace_back(0, j);
    } else {
        out.reserve((size_t)n * (size_t)(n > 0 ? n - 1 : 0) / 2);
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) out.emplace_back(i, j);
    }
    return out;
}

namespace {

std::string strip_row(const std::string& row) {
    std::string o;
    o.reserve(row.size());
    for (unsigned char ch : row)
        if (!(ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r')) o += (char)ch;
    return o;
}

}  // namespace

GroupResult prepare_group(const Group& g, const Options& opt) {
    GroupResult r;
    const size_t n = g.members.size();
    r.members.resize(n);
    r.protein_rows.assign(n, "");
    const GeneticCode* code = get_code(opt.table);
    if (!code) {
        r.error_code = "genetic_code";
        r.error = "genetic code " + std::to_string(opt.table) + " is not supported";
        return r;
    }
    // 1. translate
    std::vector<size_t> ok;
    for (size_t k = 0; k < n; ++k) {
        MemberResult& m = r.members[k];
        if (!g.members[k].failed.empty()) {
            m.error = g.members[k].failed;
            m.error_code = "failed_before";
            continue;
        }
        m.tr = translate_cds(g.members[k].cds, *code);
        if (!m.tr.error.empty()) {
            m.error = m.tr.error;
            m.error_code = m.tr.error_code;
            continue;
        }
        ok.push_back(k);
    }
    // 2. the protein alignment
    bool any_row = false, all_rows = true;
    for (size_t k : ok) {
        const bool has = !g.members[k].row.empty();
        any_row = any_row || has;
        all_rows = all_rows && has;
    }
    try {
        if (ok.size() < 2) {
            for (size_t k : ok) r.protein_rows[k] = r.members[k].tr.protein;  // nothing to align
        } else if (any_row) {
            if (!all_rows) {
                r.error_code = "alignment_incomplete";
                r.error = "only some members have a row of the protein alignment";
                return r;
            }
            size_t len = 0;
            for (size_t k : ok) {
                r.protein_rows[k] = prepare_row(strip_row(g.members[k].row));
                if (k == ok[0])
                    len = r.protein_rows[k].size();
                else if (r.protein_rows[k].size() != len) {
                    r.error_code = "alignment_ragged";
                    r.error = "the rows of the protein alignment differ in length (" + g.members[ok[0]].id + ": " +
                              std::to_string(len) + ", " + g.members[k].id + ": " +
                              std::to_string(r.protein_rows[k].size()) + ")";
                    return r;
                }
            }
            r.aligner = "given";
        } else if (n == 2) {
            const std::string& a = r.members[0].tr.protein;
            const std::string& b = r.members[1].tr.protein;
            PairAlignment pa;
            if (g.has_hit) {
                pa = expand_local_hit(a, b, g.hit, opt.flanks, opt.align);
                r.aligner = "local_hit";
            } else {
                pa = align_pair(a, b, opt.align);
                r.aligner = "nw";
            }
            r.protein_rows[0] = std::move(pa.a);
            r.protein_rows[1] = std::move(pa.b);
        } else {
            r.error_code = "alignment_missing";
            r.error = "a group of " + std::to_string(n) + " members needs a protein alignment (rows)";
            return r;
        }
    } catch (const CoalError& e) {
        r.error_code = e.code;
        r.error = std::string("protein alignment failed: ") + e.what();
        return r;
    } catch (const std::bad_alloc&) {
        r.error_code = "out_of_memory";
        r.error = "protein alignment failed: out of memory";
        return r;
    } catch (const std::exception& e) {
        r.error_code = "align_error";
        r.error = std::string("protein alignment failed: ") + e.what();
        return r;
    }
    // 3. codon rows
    for (size_t k : ok) {
        MemberResult& m = r.members[k];
        const std::string& row = r.protein_rows[k];
        try {
            std::vector<int> widths(row.size(), 3);
            Threaded t = thread_one(g.members[k].id, row, m.tr.cds, *code);
            m.codon_row = expand(row, t.dna, widths);
        } catch (const CoalError& e) {
            m.error = std::string("codon alignment failed: ") + e.what();
            m.error_code = e.code;
        } catch (const std::exception& e) {
            m.error = std::string("codon alignment failed: ") + e.what();
            m.error_code = "codon_error";
        }
    }
    return r;
}

PairResult compute_pair(const GroupResult& r, int i, int j, const Options& opt) {
    PairResult p;
    p.i = i;
    p.j = j;
    const int n = (int)r.members.size();
    if (i < 0 || j < 0 || i >= n || j >= n || i == j) {
        p.error = "pair indices out of range";
        return p;
    }
    if (!r.error.empty()) {
        p.error = r.error;
        return p;
    }
    const MemberResult& a = r.members[(size_t)i];
    const MemberResult& b = r.members[(size_t)j];
    if (!a.error.empty() || !b.error.empty()) {
        if (!a.error.empty()) p.error = "seq1: " + a.error;
        if (!b.error.empty()) p.error += (p.error.empty() ? "" : "; ") + ("seq2: " + b.error);
        return p;
    }
    const GeneticCode* code = get_code(opt.table);
    if (!code) {
        p.error = "genetic code " + std::to_string(opt.table) + " is not supported";
        return p;
    }
    const std::string& pa = r.protein_rows[(size_t)i];
    const std::string& pb = r.protein_rows[(size_t)j];
    const size_t len = pa.size() < pb.size() ? pa.size() : pb.size();
    for (size_t c = 0; c < len; ++c) {
        const char x = pa[c], y = pb[c];
        if (x == '-' && y == '-') continue;
        ++p.prot_aln_length;
        if (x == y) ++p.prot_identities;
    }
    {
        // which gene has the issue (the Ka/Ks policy counts internal-stop columns, but not whose they are, and it
        // does not see a stop that the protein alignment put against a gap)
        const std::string na = translation_note(a.tr), nb = translation_note(b.tr);
        if (!na.empty()) p.note = "seq1: " + na;
        if (!nb.empty()) p.note += (p.note.empty() ? "" : "; ") + ("seq2: " + nb);
    }
    kaks::Prepared prep = kaks::prepare(a.codon_row, b.codon_row, *code);
    // Options::compare_terminal_stops: the reference flavours of Comeron and Li (ComeronFlavour::gestimator,
    // LiFlavour::seqinr) compare stop codons, and a pipeline that translates the stop codon to '*', aligns it
    // and threads it like a residue (orthologr 0.4.3) has the two terminal stop codons as the last column of
    // the codon alignment. The codon rows here end before a terminal stop, so that column is put back when
    // both genes have a terminal stop. (One gene only: there is no column to restore; see docs/kaks.md.)
    const bool stops = opt.compare_terminal_stops && a.tr.terminal_stop && b.tr.terminal_stop &&
                       a.tr.cds.size() >= 3 && b.tr.cds.size() >= 3;
    kaks::Prepared prep_stops;
    bool have_stops = false;
    for (kaks::Method m : ordered_methods(opt.methods)) {
        const bool reference =
            (m == kaks::Method::Comeron && opt.kaks.comeron_flavour == kaks::ComeronFlavour::gestimator) ||
            (m == kaks::Method::Li && opt.kaks.li_flavour == kaks::LiFlavour::seqinr);
        if (reference && stops) {
            if (!have_stops) {
                prep_stops = kaks::prepare(a.codon_row + a.tr.cds.substr(a.tr.cds.size() - 3),
                                           b.codon_row + b.tr.cds.substr(b.tr.cds.size() - 3), *code);
                have_stops = true;
            }
            p.results.push_back(kaks::compute(prep_stops, m, *code, opt.kaks));
        } else {
            p.results.push_back(kaks::compute(prep, m, *code, opt.kaks));
        }
    }
    return p;
}

GroupResult run_group(const Group& g, const Options& opt) {
    GroupResult r = prepare_group(g, opt);
    for (auto& ij : group_pairs(g)) r.pairs.push_back(compute_pair(r, ij.first, ij.second, opt));
    return r;
}

}  // namespace dnds
}  // namespace coal
