// Vendored from COAL cpp/src/align.hpp by scripts/sync_r_sources.py. Do not edit here: edit cpp/src and sync.
// Licence: GPL-3.0-only, as COAL (https://github.com/Plumat-COAL/COAL); copyright and third-party notices
// (KaKs_Calculator 2.0 and PAML for kaks.cpp / kaks.hpp, libsequence for kaks_comeron.cpp, seqinr for
// kaks_li.cpp): inst/COPYRIGHTS of this package.
// Pairwise protein alignment for the codon step: global Needleman-Wunsch / Gotoh with affine gaps, and the
// reconstruction of a full-length pairwise alignment from a local (DIAMOND) hit. docs/alignment.md has the
// defaults, the validation against Biopython / EMBOSS needle and the benchmark.
//
// Conventions (the same as EMBOSS needle and Biopython's PairwiseAligner):
//   * a gap of length k costs gap_open + (k - 1) * gap_extend (gap_open is the cost of the FIRST gap position);
//   * scores are given as positive costs here; gap_open >= gap_extend >= 0;
//   * free_end_gaps: gaps before the first or after the last residue of either sequence cost nothing (EMBOSS
//     needle's default, -endweight off). Off: they cost the same as internal gaps (pwalign type = "global").
//   * Output is always full length: every residue of both sequences is in the alignment, overhangs are gapped,
//     never cut off (orthologr truncates them).
//
// Residues: letters A-Z (either case) and '*' (stop). Letters missing from the matrix (J, O, U) score as X.
// Anything else ('-', '.', digits, spaces) is an error (coal::CoalError, code "align_input").
//
// Like the rest of coal_core: no global mutable state, no output, errors are exceptions; every function can be
// called from several threads at once.
#pragma once
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace coal {

// A substitution matrix over "ARNDCQEGHILKMFPSTWYVBZX*" (NCBI order).
struct SubstitutionMatrix {
    const char* name;
    int8_t score[24][24];
};

// "BLOSUM62" (default), "BLOSUM45", "BLOSUM80", "PAM250" (case-insensitive); nullptr if unknown.
const SubstitutionMatrix* substitution_matrix(const std::string& name);
// The names substitution_matrix() knows, in that order.
const std::vector<std::string>& substitution_matrix_names();

struct AlignParams {
    std::string matrix = "BLOSUM62";
    double gap_open = 10.0;    // cost of a gap of length 1
    double gap_extend = 0.5;   // cost of each further gap position
    bool free_end_gaps = true;
    // Refuse a pair whose dynamic-programming table would have more cells than this (0 = no limit). The traceback
    // takes half a byte per cell: the default (4e9 cells, ~2 GB) admits titin against titin (35k x 35k ~ 1.2e9).
    uint64_t max_cells = 4000000000ULL;

    // EMBOSS needle's defaults: BLOSUM62, open 10, extend 0.5, end gaps free. COAL's default.
    static AlignParams needle() { return AlignParams(); }
    // pwalign::pairwiseAlignment(type = "global", substitutionMatrix = "BLOSUM62") with its default gapOpening = 10,
    // gapExtension = 4. pwalign charges gapOpening + k * gapExtension for a gap of length k, i.e. open 14, extend 4
    // in the convention above; end gaps are penalised.
    static AlignParams pwalign() {
        AlignParams p;
        p.gap_open = 14.0;
        p.gap_extend = 4.0;
        p.free_end_gaps = false;
        return p;
    }
};

struct PairAlignment {
    std::string a, b;          // the two gapped rows, same length; residues upper-cased
    double score = 0;          // under the AlignParams used
    size_t length = 0;         // alignment columns
    size_t identities = 0;     // columns with the same residue in both rows
    size_t similarities = 0;   // columns with a positive substitution score (identities included)
    size_t gaps = 0;           // columns with a gap in one row
    std::string error;         // set (and the rest empty) when align_pairs could not align this pair
};

// Checks the parameters (known matrix, gap_open >= gap_extend >= 0, finite, multiples of 1/1000);
// throws CoalError (code "align_params") otherwise.
void check_align_params(const AlignParams& p);

// Global alignment of a against b. Throws CoalError for bad residues ("align_input"), bad parameters
// ("align_params") or a pair above max_cells ("align_too_large").
PairAlignment align_pair(const std::string& a, const std::string& b, const AlignParams& p = AlignParams());

// Score of two given gapped rows under p (same rules as align_pair; columns that are gaps in both rows are
// ignored). Throws CoalError if the rows differ in length or hold bad characters.
double score_alignment(const std::string& row_a, const std::string& row_b, const AlignParams& p = AlignParams());

// Many pairs on `jobs` threads (<= 0: one per hardware thread). Results in input order. A pair that cannot be
// aligned gets its `error` set; the others are unaffected.
std::vector<PairAlignment> align_pairs(const std::vector<std::pair<std::string, std::string>>& pairs,
                                       const AlignParams& p = AlignParams(), int jobs = 0);

// ---------------------------------------------------------------- local hits (DIAMOND) -> full length

// One local alignment as DIAMOND reports it (--outfmt 6 ... qstart qend sstart send qseq_gapped sseq_gapped):
// 1-based inclusive coordinates on the full query / subject, and the aligned parts with '-' for gaps.
struct LocalHit {
    long long qstart = 0, qend = 0, sstart = 0, send = 0;
    std::string qseq_gapped, sseq_gapped;
};

enum class Flanks {
    Gap,    // unaligned ends against gaps, staggered: no residue of one flank is paired with the other (no
            // realignment at all)
    Align,  // unaligned ends aligned with align_pair (global; leading flanks right-anchored, trailing left-anchored
            // through free end gaps)
};

// The full-length pairwise alignment of query and subject that contains the hit's aligned block unchanged.
// Only the gap pattern of qseq_gapped / sseq_gapped is used; the residues come from query / subject, and the
// gapped strings must agree with them at the given coordinates (case-insensitive; 'X' in the hit, e.g. masking,
// matches any residue). Throws CoalError (code "local_hit") when the hit does not fit the sequences.
// score: of the result under p.
PairAlignment expand_local_hit(const std::string& query, const std::string& subject, const LocalHit& hit,
                               Flanks flanks = Flanks::Gap, const AlignParams& p = AlignParams());

}  // namespace coal
