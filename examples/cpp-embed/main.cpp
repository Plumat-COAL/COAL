// Codon-aligns two in-memory sequences with the COAL library (coal::core): the per-sequence steps that coal-cpp
// runs for every row of a protein alignment. Prints the codon alignment as FASTA on stdout, warnings on stderr.
//
//   1. prepare_row: upper-case each aligned protein row, '.' -> '-';
//   2. update_column_widths: a column is 3 nt, wider for frameshift digits (PAL2NAL's convention);
//   3. thread_one: find the protein in its coding DNA (any frame; PAL2NAL's fallback with mismatch warnings);
//      throws coal::CoalError when a sequence cannot be placed;
//   4. expand: the codon-aligned DNA row.
#include <coal/core.hpp>
#include <coal/tables.hpp>
#include <coal/version.hpp>

#include <iostream>
#include <string>
#include <vector>

int main() {
    // (id, aligned protein row, unaligned coding DNA); mouse has a 5' UTR and a stop codon, which are left out
    struct Seq {
        std::string id, aln, dna;
    };
    const std::vector<Seq> seqs = {
        {"human", "MKV-AL", "ATGAAAGTTGCTCTG"},
        {"mouse", "MKVSAL", "ccgATGAAGGTCTCAGCACTATAA"},
    };

    const coal::GeneticCode* code = coal::get_code(1);  // NCBI table 1; nullptr for an unsupported number
    if (!code) return 1;

    // Pass 1: prepared rows and column widths (every row of an alignment has the same length)
    std::vector<std::string> rows;
    std::vector<int> widths;
    for (const Seq& s : seqs) {
        rows.push_back(coal::prepare_row(s.aln));
        if (widths.empty()) widths.assign(rows.back().size(), 3);
        coal::update_column_widths(widths, rows.back());
    }

    // Pass 2: thread each DNA onto its row
    std::vector<std::string> warnings;
    try {
        for (size_t k = 0; k < seqs.size(); ++k) {
            coal::Threaded t = coal::thread_one(seqs[k].id, rows[k], seqs[k].dna, *code);
            std::cout << '>' << seqs[k].id << '\n' << coal::expand(rows[k], t.dna, widths) << '\n';
            for (const std::string& w : t.warnings) warnings.push_back(seqs[k].id + " " + w);
            for (const std::string& n : t.notices) warnings.push_back(n);
        }
    } catch (const coal::CoalError& e) {  // e.code: the stable reason ("inconsistent", ...), e.seq_id: which one
        std::cerr << "error (" << e.code << "): " << e.what() << '\n';
        return 1;
    }
    for (const std::string& w : warnings) std::cerr << "warning: " << w << '\n';
    std::cerr << "COAL " << COAL_VERSION_STRING << '\n';
    return 0;
}
