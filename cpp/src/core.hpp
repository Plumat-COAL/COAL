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

struct CoalError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct Threaded {
    std::string dna;                    // matched DNA, original case
    std::vector<std::string> warnings;  // without the leading sequence ID
    std::vector<int> mismatch_cols;     // 1-based alignment columns
};

// Upper-case the row and turn '.' into '-' (as prepare_row in core.py).
std::string prepare_row(const std::string& aln);

// Throws CoalError (message starts with seq_id) when the sequence cannot be threaded.
Threaded thread_one(const std::string& seq_id, const std::string& row, const std::string& raw_dna,
                    const GeneticCode& code);

// Widen columns for frameshift digits: width 3*ceil(N/3) (min 3).
void update_column_widths(std::vector<int>& widths, const std::string& row);

// Codon-aligned row (before any column removal).
std::string expand(const std::string& row, const std::string& dna, const std::vector<int>& widths);

// Mark (bad[k] = 1) the 3-nt output units -nogap removes: a gap, or a stop codon
// by PAL2NAL's case-sensitive rule ([TU]A[AGR] or [TU]GA).
void mark_nogap_triplets(const std::string& expanded, std::vector<char>& bad);

// Character ranges [a, b) of the expanded row that survive -nomismatch / -nogap.
std::vector<std::pair<long long, long long>> keep_slices(const std::vector<int>& widths,
                                                         const std::vector<char>& drop_cols,
                                                         const std::vector<char>& bad_triplets);

}  // namespace coal
