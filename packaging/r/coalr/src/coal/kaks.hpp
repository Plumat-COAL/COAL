// Vendored from COAL cpp/src/kaks.hpp by scripts/sync_r_sources.py. Do not edit here: edit cpp/src and sync.
// Licence: GPL-3.0-only, as COAL (https://github.com/Plumat-COAL/COAL); copyright and third-party notices
// (KaKs_Calculator 2.0 and PAML for kaks.cpp / kaks.hpp, libsequence for kaks_comeron.cpp, seqinr for
// kaks_li.cpp): inst/COPYRIGHTS of this package.
// SPDX-License-Identifier: GPL-3.0-only
//
// Ka/Ks (dN/dS) for pairs of aligned coding sequences: NG86, LWL85, MLWL85, LPB93, MLPB93, YN00, MYN, Comeron
// (1995) and Li (1993), in process (no external binaries). Comeron (kaks_comeron.cpp) is derived from libsequence
// (Copyright (C) 2003-2009 Kevin Thornton; GPL v3 or later) as the R package orthologr compiles it in (gestimator;
// modified by Sarah Scharfenberg and Hajk-Georg Drost, 2014); Li (kaks_li.cpp) is derived from kaks.c of the R
// package seqinr (Manolo Gouy's code; seqinr: D. Charif, J. R. Lobry and contributors; GPL >= 2). The headers of
// those two files and the NOTICE file say what was taken and what was changed.
// The implementation of the other methods (kaks.cpp) is derived from KaKs_Calculator 2.0
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

// KaKs_Calculator's method names: NG, LWL, MLWL, LPB, MLPB, YN, MYN. Two methods that KaKs_Calculator does not
// have, under the names the R package orthologr uses for them: Comeron (Comeron 1995, as libsequence's gestimator
// computes it; alias C95) and Li (Li 1993, as seqinr::kaks computes it). Both are for the standard genetic code
// only (as their references): with another code the result has status `error`.
enum class Method { NG, LWL, MLWL, LPB, MLPB, YN, MYN, Comeron, Li };

// The KaKs_Calculator name ("NG", "YN", ...), or "Comeron" / "Li".
const char* method_name(Method m);
// Case-insensitive name -> Method ("C95" is Comeron). False for a name that is not a method, or one this engine
// does not implement natively (GY, MS, MA, NONE, the gamma variants G*); `why` then says which.
bool parse_method(const std::string& name, Method& out, std::string& why);
// The order KaKs_Calculator writes its methods in for each pair (whatever the order of -m on its command line).
const std::vector<Method>& output_order();

// ok: estimated. warning: estimated, but see note (internal stop codons were skipped).
// not_converged: the YN00 / MYN iteration ended at its round limit without settling; there is no estimate (Ka, Ks,
// their standard errors, Sd and Nd are NaN; S, N and kappa are set), because the values of the last round are
// not reproducible between platforms (docs/kaks.md "Platform reproducibility"). They are in last_ka / last_ks.
// saturated: Ka and/or Ks could not be corrected for multiple hits (log of a non-positive number): that value is
// NaN. no_data: nothing to compare (no codon left after the policy, or no synonymous / nonsynonymous sites).
// error: the pair itself is malformed (rows of different lengths, length not a multiple of 3).
// The order is the order of severity (a later one replaces an earlier one), with one exception: an iteration
// that did not converge is not_converged even if its last round was saturated, since that too differs by platform.
enum class Status { ok, warning, not_converged, saturated, no_data, error };
constexpr int STATUS_COUNT = 6;
const char* status_name(Status s);

// Which reference the YN00 iteration follows where KaKs_Calculator and PAML's yn00 differ (docs/kaks.md):
// the start value of t, the kappa used when none can be estimated (2 vs 4.6) and the cap on kappa (2 vs 999).
enum class Flavour { kaks_calculator, paml };

// Comeron: `coal` applies COAL's policy (prepare(): terminal stop dropped, internal stops skipped) and gives NA
// with status `saturated` where a logarithm of the correction has an argument <= 0. `gestimator` reproduces the
// numbers of orthologr 0.4.3's gestimator (libsequence): stop codons are compared like any other codon, U is not
// read as T (the codon is skipped), and an undefined logarithm is replaced by 0, so a finite number comes out
// (status `warning`); only where gestimator itself prints 999 is the value NA. docs/kaks.md lists the differences.
enum class ComeronFlavour { coal, gestimator };
// Comeron, codons that differ at two positions: the four branches of the two pathways. `comeron`: codon1 -
// int1 - codon2 and codon1 - int2 - codon2, as Comeron (1995) defines them; the estimate does not depend on
// which sequence is first. `libsequence`: the branches libsequence's TwoSubs counts (codon1 - int1, int1 - int2,
// codon1 - codon2, codon2 - int2; see kaks_comeron.cpp), gestimator's numbers, which change when the two
// sequences are swapped. `by_flavour`: comeron with ComeronFlavour::coal, libsequence with
// ComeronFlavour::gestimator. docs/kaks.md "Two-position codon differences".
enum class ComeronTwoSubs { by_flavour, comeron, libsequence };
// the order a run uses (by_flavour resolved)
inline ComeronTwoSubs resolved_twosubs(ComeronTwoSubs t, ComeronFlavour f) {
    if (t != ComeronTwoSubs::by_flavour) return t;
    return f == ComeronFlavour::gestimator ? ComeronTwoSubs::libsequence : ComeronTwoSubs::comeron;
}
inline const char* twosubs_name(ComeronTwoSubs t) {
    return t == ComeronTwoSubs::libsequence ? "libsequence" : t == ComeronTwoSubs::comeron ? "comeron" : "by_flavour";
}
// Comeron: a site class (L0, L2S, L2V, L4) with fewer than this many sites cannot make dN or dS NA. Its
// undefined logarithm is set to 0, as gestimator does for every class, and the row says so (status warning,
// note "Ka: a class with fewer than 3 sites has an undefined correction ..."). 0: every class with sites can
// (COAL 0.5.x before this option). Never changes a value that is not NA otherwise. docs/kaks.md.
constexpr double COMERON_MIN_SITES_DEFAULT = 3.0;
// Li: `coal` applies COAL's policy; `seqinr` compares stop codons as seqinr::kaks does and does not read U as T.
// In both, seqinr's sentinel 9.999999 ("saturated") is NA.
enum class LiFlavour { coal, seqinr };

struct Options {
    Flavour yn_flavour = Flavour::kaks_calculator;
    // YN00 pathway weighting (iteration with P(t)). KaKs_Calculator always weights; PAML yn00's weighting = 1.
    // false = equal weighting of pathways in one round (PAML yn00's default control file, weighting = 0).
    bool yn_weighting = true;
    // Testing only: compute the transition-probability matrix P(t) of YN / MYN in every round, also for a pair in
    // which no codon differs at two or three positions, where nothing reads it (kaks.cpp needs_pathway_weights).
    // The results are the same bits either way; the tests compare the two (`--yn-always-pmat`).
    bool yn_always_pmat = false;
    ComeronFlavour comeron_flavour = ComeronFlavour::coal;
    ComeronTwoSubs comeron_twosubs = ComeronTwoSubs::by_flavour;
    double comeron_min_sites = COMERON_MIN_SITES_DEFAULT;
    LiFlavour li_flavour = LiFlavour::coal;
};

constexpr double NA = std::numeric_limits<double>::quiet_NaN();

struct Result {
    Method method = Method::YN;
    Status status = Status::ok;
    std::string note;                // why the status is not ok ("" when ok)

    double ka = NA, ks = NA, ka_ks = NA;  // dN, dS, dN/dS (NaN = not available)
    double se_ka = NA, se_ks = NA;        // standard errors where the method defines one (NG86, YN00, MYN, Li)
    double S = NA, N = NA;                // synonymous / nonsynonymous sites
    double Sd = NA, Nd = NA;              // synonymous / nonsynonymous differences
    double substitutions = 0;             // nucleotide differences in the compared codons
    double p_fisher = NA;                 // Fisher's exact test on Sd, Nd, S-Sd, N-Nd (KaKs_Calculator's version)
    double divergence = NA;               // (S*Ks + N*Ka) / (S + N), as KaKs_Calculator's Divergence-Time
    double kappa = NA;                    // ts/tv ratio used (MYN: kappa from TN93, kappa_tc / kappa_ag in rates)
    double rates[6] = {1, 1, 1, 1, 1, 1}; // KaKs_Calculator's Substitution-Rate-Ratio
    bool has_fold = false;                // LWL / LPB family and Li: the next three are set
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
    double last_ka = NA, last_ks = NA;    // not converged: Ka and Ks of the last round (NaN: saturated there or
                                          // no such sites). Not estimates, and not the same on every platform
};

// A pair after the policy (docs/kaks.md "Policy"): the codon columns to compare, as codon indices 0..63
// (NCBI order, T C A G), and the counts of what was left out. Prepare once, compute many methods.
struct Prepared {
    std::vector<uint8_t> c1, c2;
    Result base;  // status/note/counts/gc filled in; estimates NaN

    // The pair as the reference programs of Comeron and Li read it (ComeronFlavour::gestimator,
    // LiFlavour::seqinr): every codon column in which both codons consist of A, C, G, T (any case; U is not T
    // here), stop codons included. raw_skipped: the other columns.
    std::vector<uint8_t> r1, r2;
    long raw_skipped = 0;
    // Kimura's (1980) two-parameter counts over the whole rows as libsequence's Kimura80 takes them (every
    // position where neither row has '-'): gestimator uses that distance in one ad-hoc rule (kaks_comeron.cpp).
    long k80_ts = 0, k80_tv = 0, k80_sites = 0;
};

// Rows: aligned coding DNA of equal length (A C G T U, any case; '-' or '.' gaps; other letters ambiguous).
Prepared prepare(const std::string& row1, const std::string& row2, const GeneticCode& code);
Result compute(const Prepared& p, Method m, const GeneticCode& code, const Options& opt = {});
Result compute(const std::string& row1, const std::string& row2, Method m, const GeneticCode& code,
               const Options& opt = {});

// Internal: the two methods that live in their own files (kaks_comeron.cpp, kaks_li.cpp). Use compute().
namespace detail {
Result run_comeron(const Prepared& p, const GeneticCode& code, const Options& opt, Result r);
Result run_li(const Prepared& p, const GeneticCode& code, const Options& opt, Result r);
}  // namespace detail

}  // namespace kaks
}  // namespace coal
