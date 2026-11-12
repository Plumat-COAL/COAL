// Vendored from COAL cpp/src/dnds.hpp by scripts/sync_r_sources.py. Do not edit here: edit cpp/src and sync.
// Licence: GPL-3.0-only, as COAL (https://github.com/Plumat-COAL/COAL); copyright and third-party notices
// (KaKs_Calculator 2.0 and PAML for kaks.cpp / kaks.hpp): inst/COPYRIGHTS of this package.
// The one-step dN/dS core: coding sequences (CDS) -> translation -> protein alignment -> codon alignment ->
// Ka/Ks, for a pair of genes or for a group (cluster) of genes. It ties together tables.hpp (genetic codes),
// align.hpp (pairwise alignment, DIAMOND-hit expansion), core.hpp (codon threading) and kaks.hpp (Ka/Ks).
// `coal pipeline` (src/coal/dnds_pipeline.py, docs/pipeline.md) drives it through `coal-cpp codon-dnds`.
//
//   #include <coal/dnds.hpp>
//   coal::dnds::Options opt;                                       // YN00, genetic code 1, BLOSUM62 needle defaults
//   coal::dnds::Group g;
//   g.members = {{"geneA", cdsA, ""}, {"geneB", cdsB, ""}};        // two members, no rows: aligned here (NW)
//   coal::dnds::GroupResult r = coal::dnds::run_group(g, opt);     // never throws for bad data
//   // r.pairs[0].results[0].ka / .ks / .ka_ks, r.pairs[0].error, r.codon_rows, r.members[i].error
//
// Like the rest of coal_core: no global mutable state, no output, thread-safe. Bad data never throws: a gene that
// cannot be used has Member::error set, a pair that cannot be computed PairResult::error, a group GroupResult::error;
// the others are unaffected.
#pragma once
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "align.hpp"
#include "kaks.hpp"
#include "tables.hpp"

namespace coal {
namespace dnds {

// A CDS translated codon by codon from its first base.
struct Translation {
    std::string cds;            // the CDS as used: whitespace and an incomplete last codon removed, letters as given
    std::string protein;        // one letter per codon; a stop codon at the end is left out; an internal stop is
                                // '*'; a codon with an ambiguous base is X unless every reading gives one residue
    long nt_length = 0;         // letters in the input
    int trailing_nt = 0;        // nt_length % 3: bases after the last complete codon (ignored)
    bool terminal_stop = false; // the last codon is a stop (left out of `protein`)
    long internal_stops = 0;    // stop codons before the last codon
    long ambiguous_codons = 0;  // codons with a letter other than A, C, G, T/U
    std::string error;          // not usable: why (empty sequence, a character that is no nucleotide code, ...)
    std::string error_code;     // cds_empty | cds_invalid_character | cds_too_short | cds_only_stop
};

// Whitespace in `cds` is ignored. Never throws.
Translation translate_cds(const std::string& cds, const GeneticCode& code);
// The issues of a usable translation as one note ("" when there is none): length not a multiple of 3, internal
// stop codons. (Ambiguous codons and a terminal stop are normal and not mentioned.)
std::string translation_note(const Translation& t);

struct Member {
    std::string id;
    std::string cds;     // unaligned coding DNA
    std::string row;     // this member's row of the group's protein alignment ('-' gaps), or "" (see Group)
    std::string failed;  // set by the caller: the gene is already known to be unusable (why); cds/row are ignored
};

enum class Pairs {
    all,             // every pair i < j
    representative,  // member 0 against every other member
    list,            // Group::pair_list
};

struct Group {
    std::string id;
    std::vector<Member> members;
    Pairs pairs = Pairs::all;
    std::vector<std::pair<int, int>> pair_list;  // member indices (Pairs::list)
    // Two members without rows: use this local hit (member 0 = query, member 1 = subject) instead of aligning.
    bool has_hit = false;
    LocalHit hit;
};

struct Options {
    std::vector<kaks::Method> methods = {kaks::Method::YN};
    int table = 1;                   // NCBI genetic code
    AlignParams align;               // for two-member groups without rows, and the flanks of a local hit
    Flanks flanks = Flanks::Align;   // how a local hit's unaligned ends are completed
    kaks::Options kaks;
};

struct MemberResult {
    Translation tr;
    std::string codon_row;   // the member's row of the codon alignment ("" when it failed)
    std::string error;       // why this member is not in the codon alignment
    std::string error_code;  // the translation's code, "failed_before" (Member::failed), or the codon engine's
};

struct PairResult {
    int i = 0, j = 0;                    // member indices
    size_t prot_aln_length = 0;          // protein alignment columns with a residue in at least one of the two
    size_t prot_identities = 0;          // of them, columns with the same residue in both
    std::vector<kaks::Result> results;   // one per Options::methods, in kaks::output_order(); empty when `error`
    std::string error;                   // the pair could not be computed (a member failed, alignment failed, ...)
    std::string note;                    // CDS issues of the two genes (added to every result's note by the CLI)
};

struct GroupResult {
    std::vector<MemberResult> members;
    std::string aligner;                 // "given" | "nw" | "local_hit" | "" (nothing to align)
    std::string error;                   // the whole group failed (rows of different lengths, no alignment, ...)
    std::string error_code;
    std::vector<std::string> protein_rows;  // the protein alignment used (one row per member; "" for failed ones)
    std::vector<PairResult> pairs;       // filled by run_group / compute_pair
};

// The pairs Group::pairs asks for, in output order.
std::vector<std::pair<int, int>> group_pairs(const Group& g);

// Steps 1-3: translate, align (when no rows are given), thread the codons. No pairs yet.
GroupResult prepare_group(const Group& g, const Options& opt);
// Step 4 for one pair of a prepared group. Thread-safe on a const GroupResult.
PairResult compute_pair(const GroupResult& r, int i, int j, const Options& opt);
// All four steps (pairs computed one after the other; for threads call compute_pair yourself).
GroupResult run_group(const Group& g, const Options& opt);

// The requested methods without duplicates, in KaKs_Calculator's output order.
std::vector<kaks::Method> ordered_methods(const std::vector<kaks::Method>& methods);

}  // namespace dnds
}  // namespace coal
