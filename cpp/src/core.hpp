// Core conversion: protein alignment row + coding DNA -> codon-aligned row.
// Same algorithm and behaviour as src/coal/core.py (see its docstring):
//   1. locate the whole protein in the DNA (any frame, leftmost fit);
//   2. anchor fallback (10-residue chunks) with mismatch warnings;
//   3. thread codons into the alignment columns.
#pragma once
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "tables.hpp"

namespace coal {

// An error reported to the user. `code` is the stable reason used in the run
// report (same codes as the Python engine); `seq_id` is set when one sequence failed.
struct CoalError : std::runtime_error {
    std::string code = "error";
    std::string seq_id;
    explicit CoalError(const std::string& msg, std::string c = "error", std::string sid = "")
        : std::runtime_error(msg), code(std::move(c)), seq_id(std::move(sid)) {}
};

struct Threaded {
    std::string dna;                    // matched DNA, original case
    std::vector<std::string> warnings;  // without the leading sequence ID
    std::vector<int> mismatch_cols;     // 1-based alignment columns
    std::vector<std::string> notices;   // probably-wrong data (with the sequence ID), as the reference
};

// Upper-case the row and turn '.' into '-' (as prepare_row in core.py).
std::string prepare_row(const std::string& aln);

// As core.py: at least 20 letters and >= 90% of them A/C/G/T/U/N; at least 20 letters and
// >= 10% of them amino acids that are no nucleotide code (E F I J L O P Q Z).
bool looks_like_dna(const std::string& text);
bool looks_like_protein(const std::string& text);
// The error for protein rows of different lengths (core.py row_length_message).
std::string row_length_message(const std::string& label, long long ncol, long long n, const std::string& pid,
                               const std::string& rows);

// Throws CoalError (message starts with seq_id) when the sequence cannot be threaded. strict (--strict):
// probably-wrong data is an error instead of a notice.
Threaded thread_one(const std::string& seq_id, const std::string& row, const std::string& raw_dna,
                    const GeneticCode& code, bool strict = false);

// Widen columns for frameshift digits: width 3*ceil(N/3) (min 3).
void update_column_widths(std::vector<int>& widths, const std::string& row);

// Codon-aligned row (before any column removal).
std::string expand(const std::string& row, const std::string& dna, const std::vector<int>& widths);

// Mark (bad[k] = 1) the 3-nt output units -nogap removes: a gap, or a stop codon of the selected
// genetic code (any case, U read as T; an IUPAC triplet counts if every reading is a stop).
// (PAL2NAL v14 removes only the standard code's TAA/TAG/TGA, upper case: a documented difference.)
void mark_nogap_triplets(const std::string& expanded, std::vector<char>& bad, const GeneticCode& code);

// Character ranges [a, b) of the expanded row that survive -nomismatch / -nogap.
std::vector<std::pair<long long, long long>> keep_slices(const std::vector<int>& widths,
                                                         const std::vector<char>& drop_cols,
                                                         const std::vector<char>& bad_triplets);

}  // namespace coal
