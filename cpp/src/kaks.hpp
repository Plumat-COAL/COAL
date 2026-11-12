// SPDX-License-Identifier: GPL-3.0-only
//
// Ka/Ks (dN/dS) for pairs of aligned coding sequences: NG86, LWL85, MLWL85, LPB93, MLPB93, YN00 and MYN, in
// process (no external binaries). The implementation (kaks.cpp) is derived from KaKs_Calculator 2.0
// (Copyright (C) 2009, Beijing Institute of Genomics of CAS; Da-Peng Wang, Yu-Bin Zhang, Zhang Zhang,
// Kristian K Ullrich; GPL v3) and from PAML's yn00.c (Copyright, 1998, Ziheng Yang; GPL v3); kaks.cpp's header
// and the NOTICE file say what was taken and what was changed. Results match KaKs_Calculator 2.0 (and PAML's
// yn00 for YN00, with Options::yn_flavour = Flavour::paml). Methods, the gap / stop / ambiguity policy, the
// validation and every deliberate difference from the reference tools: docs/kaks.md.
//
// Thread-safe: no global or static mutable state, and no call into a C library function that has one (the
// log-gamma of Fisher's test is computed in kaks.cpp, not with lgamma(), which writes `signgam`); compute() may
// be called from many threads at once (coal-cpp kaks/dnds does, from its thread pool). Never throws for bad data: problems are reported in Result::status / Result::note, so one
// bad pair never stops a run. (parse_method() is the only function that reports through its return value.)
//
//   #include <coal/kaks.hpp>
//   const coal::GeneticCode* code = coal::get_code(1);
//   coal::kaks::Result r = coal::kaks::compute(row1, row2, coal::kaks::Method::YN, *code);
//   if (r.status != coal::kaks::Status::error) use(r.ka, r.ks, r.ka_ks);   // NaN = not available
#pragma once
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "tables.hpp"

namespace coal {
namespace kaks {

// KaKs_Calculator's method names: NG, LWL, MLWL, LPB, MLPB, YN, MYN.
enum class Method { NG, LWL, MLWL, LPB, MLPB, YN, MYN };

// The KaKs_Calculator name ("NG", "YN", ...).
const char* method_name(Method m);
// Case-insensitive KaKs_Calculator name -> Method. False for a name that is not a method, or one this engine does
// not implement natively (GY, MS, MA, NONE, the gamma variants G*); `why` then says which.
bool parse_method(const std::string& name, Method& out, std::string& why);
// The order KaKs_Calculator writes its methods in for each pair (whatever the order of -m on its command line).
const std::vector<Method>& output_order();

// ok: estimated. warning: estimated, but see note (internal stop codons were skipped).
// saturated: Ka and/or Ks could not be corrected for multiple hits (log of a non-positive number): that value is
// NaN. no_data: nothing to compare (no codon left after the policy, or no synonymous / nonsynonymous sites).
// error: the pair itself is malformed (rows of different lengths, length not a multiple of 3).
enum class Status { ok, warning, saturated, no_data, error };
const char* status_name(Status s);

// Which reference the YN00 iteration follows where KaKs_Calculator and PAML's yn00 differ (docs/kaks.md):
// the start value of t, the kappa used when none can be estimated (2 vs 4.6) and the cap on kappa (2 vs 999).
enum class Flavour { kaks_calculator, paml };

struct Options {
    Flavour yn_flavour = Flavour::kaks_calculator;
    // YN00 pathway weighting (iteration with P(t)). KaKs_Calculator always weights; PAML yn00's weighting = 1.
    // false = equal weighting of pathways in one round (PAML yn00's default control file, weighting = 0).
    bool yn_weighting = true;
};

constexpr double NA = std::numeric_limits<double>::quiet_NaN();

struct Result {
    Method method = Method::YN;
    Status status = Status::ok;
    std::string note;                // why the status is not ok ("" when ok)

    double ka = NA, ks = NA, ka_ks = NA;  // dN, dS, dN/dS (NaN = not available)
    double se_ka = NA, se_ks = NA;        // standard errors where the method defines one (NG86, YN00, MYN)
    double S = NA, N = NA;                // synonymous / nonsynonymous sites
    double Sd = NA, Nd = NA;              // synonymous / nonsynonymous differences
    double substitutions = 0;             // nucleotide differences in the compared codons
    double p_fisher = NA;                 // Fisher's exact test on Sd, Nd, S-Sd, N-Nd (KaKs_Calculator's version)
    double divergence = NA;               // (S*Ks + N*Ka) / (S + N), as KaKs_Calculator's Divergence-Time
    double kappa = NA;                    // ts/tv ratio used (MYN: kappa from TN93, kappa_tc / kappa_ag in rates)
    double rates[6] = {1, 1, 1, 1, 1, 1}; // KaKs_Calculator's Substitution-Rate-Ratio
    bool has_fold = false;                // LWL / LPB family: the next three are set
    double fold_sites[3] = {0, 0, 0};     // L0, L2, L4
    double fold_s[3] = {0, 0, 0};         // transitional differences per class (Si)
    double fold_n[3] = {0, 0, 0};         // transversional differences per class (Vi)
    double gc[4] = {NA, NA, NA, NA};      // GC content: all, codon positions 1, 2, 3 (as KaKs_Calculator)

    long codons_compared = 0;             // codon columns used
    long length = 0;                      // nucleotides compared (3 * codons_compared), KaKs_Calculator's Length
    long internal_stops = 0;              // columns skipped because one sequence has a stop codon there
    long gap_codons = 0;                  // columns skipped for a gap in either sequence
    long ambiguous_codons = 0;            // columns skipped for a non-ACGT/U letter in either sequence
    bool terminal_stop_dropped = false;   // a final stop codon was left out
    int iterations = 0;                   // YN00 / MYN rounds
    bool converged = true;                // YN00 / MYN: the iteration met its tolerance
};

// A pair after the policy (docs/kaks.md "Policy"): the codon columns to compare, as codon indices 0..63
// (NCBI order, T C A G), and the counts of what was left out. Prepare once, compute many methods.
struct Prepared {
    std::vector<uint8_t> c1, c2;
    Result base;  // status/note/counts/gc filled in; estimates NaN
};

// Rows: aligned coding DNA of equal length (A C G T U, any case; '-' or '.' gaps; other letters ambiguous).
Prepared prepare(const std::string& row1, const std::string& row2, const GeneticCode& code);
Result compute(const Prepared& p, Method m, const GeneticCode& code, const Options& opt = {});
Result compute(const std::string& row1, const std::string& row2, Method m, const GeneticCode& code,
               const Options& opt = {});

}  // namespace kaks
}  // namespace coal
