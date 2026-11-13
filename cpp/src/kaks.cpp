// SPDX-License-Identifier: GPL-3.0-only
//
// Ka/Ks engine (kaks.hpp): NG86, LWL85, MLWL85, LPB93, MLPB93, YN00 and MYN.
//
// This file is derived from two programs that are distributed under the GNU General Public License, version 3,
// and it is distributed under the same licence (COAL as a whole is GPL-3.0-only; see LICENSE and NOTICE):
//
//   KaKs_Calculator 2.0 (the kakscalculator2 2.0.1 sources, https://github.com/kullrich/kakscalculator2;
//   original project https://sourceforge.net/projects/kakscalculator2/): base.cpp (Fisher's exact test, the
//   derived output values), NG86.cpp, LWL85.cpp, LPB93.cpp, YN00.cpp and MYN.cpp. (KaKs.cpp, the seventh upstream
//   file named in NOTICE, is not used here: its header line and column order are reproduced in kaks_cli.cpp.)
//     Copyright (C) 2009, Beijing Institute of Genomics of CAS. All rights reserved.
//     Version 2.0: Da-Peng Wang, Yu-Bin Zhang. Version 1.0: Zhang Zhang.
//     Modified versions 2.0.1 and 2.0.2: Kristian K Ullrich.
//
//   PAML 4.10 (https://github.com/abacus-gene/paml), yn00.c: "Pairwise estimation of dS and dN by the method of
//   Yang & Nielsen (2000 Mol. Biol. Evol. 17:32-43)". Copyright, 1998, Ziheng Yang. KaKs_Calculator's YN00.cpp is
//   itself adapted from yn00.c; the functions here follow yn00.c's GetKappa, GetFreqs, CountSites, GetPMatCodon,
//   CountDiffs, DistanceF84 and DistanceYN00 (and eigenQREV of PAML's tools.c) in structure, constants and order
//   of operations, so that the numbers are the same.
//
// Modified for COAL by Chloe Newton and Keanu Murati, 2026. What was changed:
//   - rewritten as free functions on codon indices 0..63 (NCBI order T C A G) for any genetic code of COAL's
//     tables; no classes, no global or static state (the references keep the pair, kappa, t and omega in
//     members and globals), so pairs can be computed on several threads;
//   - one gap / ambiguity / stop-codon policy applied to the pair before every method (prepare());
//   - an undefined multiple-hit correction is reported as NA with the status `saturated` where the references
//     substitute a capped value (JC69 with P = 0.75(n-1)/n, t = 99, or K = 0);
//   - Fisher's exact test: log-factorials by log_gamma() below (the references use an approximation; the C
//     library's lgamma is not thread-safe), tables as probable as the observed one within 1e-7 count, and no
//     test when a cell would be negative;
//   - YN00 and MYN: standard errors are returned; with sites of one kind only (S = 0 or N = 0) the estimate
//     that has sites is still made; base frequencies at the sites are averaged over the sequences that have
//     such sites; PAML's variant of the edge cases is an option (Options::yn_flavour, yn_weighting);
//   - MYN: a class of sites with weight 0 is left out of the mean kappa (the reference adds 0 x NaN), and TN93
//     with no T-C pair or no A-G pair is the F84 distance (the reference divides by zero);
//   - the symmetric eigen solver is the EISPACK tred2/tql2 algorithm on transposed storage.
// docs/kaks.md lists every place where this engine gives something else than the references, and why.
//
// References:
//   NG86  Nei M, Gojobori T (1986) Mol Biol Evol 3:418-426
//   LWL85 Li W-H, Wu C-I, Luo C-C (1985) Mol Biol Evol 2:150-174
//   LPB93 Li W-H (1993) J Mol Evol 36:96-99; Pamilo P, Bianchi NO (1993) Mol Biol Evol 10:271-281
//   MLWL / MLPB  Tzeng Y-H, Pan R, Li W-H (2004) Mol Biol Evol 21:2290-2298
//   YN00  Yang Z, Nielsen R (2000) Mol Biol Evol 17:32-43
//   MYN   Zhang Z, Li J, Yu J (2006) BMC Evol Biol 6:44
#include "kaks.hpp"

// No floating-point contraction in this file: a * b + c is computed as written, never fused into one FMA
// instruction. With contraction the last digits of every estimate depend on whether the compiler and the target
// have FMA (Clang contracts by default; arm64 always has the instruction, x86-64 only with -march). Said here as
// well as in CMakeLists.txt so that it also holds where this file is compiled by something else (R's build of
// coalr, an embedding project). docs/kaks.md "Platform reproducibility" says what is and is not guaranteed.
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#elif defined(_MSC_VER)
#pragma fp_contract(off)
#endif

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <initializer_list>

namespace coal {
namespace kaks {

namespace {

constexpr double SMALL = 1e-6;  // KaKs_Calculator's SMALLVALUE

inline bool is_stop(const GeneticCode& g, int c) { return g.aa[c] == '*'; }
inline int base_at(int c, int pos) { return pos == 0 ? c / 16 : pos == 1 ? (c / 4) % 4 : c % 4; }
inline int with_base(int c, int pos, int b) {
    static const int by[3] = {16, 4, 1};
    return c + (b - base_at(c, pos)) * by[pos];
}
// T=0 C=1 A=2 G=3: T<->C sums to 1, A<->G to 5
inline bool transition(int a, int b) { return a + b == 1 || a + b == 5; }

// ---------------------------------------------------------------------------------------------------------
// Fisher's exact test as KaKs_Calculator computes it on the (non-integer) table [[Sd, S-Sd], [Nd, N-Nd]]:
// two-tailed sum of the tables with the same margins, stepping the first cell by 1 from 0, that are no more
// probable than the observed one, skipping tables with a cell <= 1e-6; log-factorials via lgamma. A table counts
// as "no more probable" within a relative 1e-7, so the observed table itself is never lost to rounding.
//
// log Gamma(x) for x > 0, computed here and not with the C library's lgamma(): that function writes the global
// `signgam` on glibc, macOS, the BSDs and MinGW, which is a data race when pairs are computed on several threads,
// and its thread-safe variant lgamma_r is not available everywhere (MSVC has none). Stirling's series above 16,
// reached with the recurrence Gamma(x + 1) = x Gamma(x). Checked against glibc's lgamma on 2e7 values of
// 1e-7 < x < 1e7: the difference is at most 1.7e-14 of max(1, |lgamma|), a few units in the last place. Pure: no
// global or static state.
double log_gamma(double x) {
    double prod = 1;
    while (x < 16) {
        prod *= x;
        x += 1;
    }
    const double xi = 1 / x, x2 = xi * xi;
    const double series =
        xi * (1 / 12.0 -
              x2 * (1 / 360.0 -
                    x2 * (1 / 1260.0 - x2 * (1 / 1680.0 - x2 * (1 / 1188.0 - x2 * (691 / 360360.0 - x2 / 156.0))))));
    return (x - 0.5) * std::log(x) - x + 0.91893853320467274178 + series - std::log(prod);
}

double lfact(double n) { return n > 0 ? log_gamma(n + 1) : 0.0; }

double fisher(double sd, double nd, double s, double n) {
    double m[4] = {sd, s, nd, n};
    const double R0 = m[0] + m[2], R1 = m[1] + m[3], C0 = m[0] + m[1], C1 = m[2] + m[3];
    const double sum = R0 + R1;
    const double num = lfact(R0) + lfact(R1) + lfact(C0) + lfact(C1);
    const double fs = lfact(sum);
    double den = fs;
    for (double v : m) den += lfact(v);
    const double cur = std::exp(num - den);
    double total = 0;
    for (long i = 0; (double)i < R0 + SMALL; ++i) {
        const double a = (double)i, b = C0 - a, c = R0 - a, d = R1 - C0 + a;
        if (a > SMALL && b > SMALL && c > SMALL && d > SMALL) {
            const double p = std::exp(num - (fs + lfact(a) + lfact(b) + lfact(c) + lfact(d)));
            if (p <= cur * (1 + 1e-7)) total += p;  // ties with the observed table count (as R fisher.test)
        }
    }
    return total;
}

// Fill the common output fields and the KaKs_Calculator-style derived values.
void finish(Result& r) {
    for (double* x : {&r.ka, &r.ks, &r.se_ka, &r.se_ks, &r.divergence})
        if (*x == 0) *x = 0.0;  // no negative zeros (-0.75 * log(1))
    if (std::isfinite(r.ka) && std::isfinite(r.ks) && r.ks >= SMALL) r.ka_ks = r.ka / r.ks;
    // (no test when a cell would be negative: LPB93's S = Sd/Ks is below Sd once Ks > 1; KaKs_Calculator then
    // feeds the negative cell to its log-factorial, which returns 1 for it)
    // A cell that is zero up to rounding (weighted counts that add up to the number of sites: 2 against
    // 2.0000000000000004) is zero, not negative: otherwise the test is made or not depending on the last bit.
    const double restS = r.S - r.Sd, restN = r.N - r.Nd, tie = 1e-9;
    if (r.Sd >= SMALL && r.Nd >= SMALL && r.S >= SMALL && r.N >= SMALL && restS > -tie && restN > -tie)
        r.p_fisher = fisher(r.Sd, r.Nd, std::max(restS, 0.0), std::max(restN, 0.0));
}

void mark_saturated(Result& r, const char* which) {
    if (r.status != Status::no_data && r.status != Status::error) r.status = Status::saturated;
    std::string msg = std::string(which) + " saturated (multiple-hit correction undefined)";
    r.note = r.note.empty() ? msg : r.note + "; " + msg;
}

void add_note(Result& r, const std::string& msg) { r.note = r.note.empty() ? msg : r.note + "; " + msg; }

// The YN00 / MYN iteration (pathway weights from the current estimates, estimates from the weighted counts)
// ended at its round limit without settling. This happens with very short or nearly saturated pairs, where the
// estimates jump between rounds. The numbers of the last round are then not a fixed point of anything, and they
// are not reproducible either: the jumps amplify a difference in the 16th digit (another C library's exp or log,
// a fused multiply-add) into the third digit, and can move a value across the saturation limit. Measured for the
// pair AAA / ATG with YN: Ka is 3.169, 3.174 or "saturated" depending on the platform. So there is no estimate:
// Ka, Ks, their standard errors and the weighted difference counts Sd and Nd are NA, the status is
// `not_converged` (whatever the last round looked like, saturated or not; only no_data and error, which are
// decided before the iteration, stay), and the last round's values are kept in Result::last_ka / last_ks and, with
// six digits, in the note. S, N and kappa do not depend on the iteration and are kept.
void not_converged(Result& r, const char* what, int rounds, bool hasN, double dN, bool satN, bool hasS, double dS,
                   bool satS) {
    if (r.status == Status::ok || r.status == Status::warning) r.status = Status::not_converged;
    auto last = [](bool has, double v, bool sat) -> std::string {
        if (!has) return "NA";
        if (sat || !std::isfinite(v)) return "saturated";
        char buf[40];
        std::snprintf(buf, sizeof buf, "%.6g", v);
        return buf;
    };
    r.last_ka = (hasN && !satN && std::isfinite(dN)) ? dN : NA;
    r.last_ks = (hasS && !satS && std::isfinite(dS)) ? dS : NA;
    r.ka = r.ks = r.se_ka = r.se_ks = r.divergence = r.Sd = r.Nd = NA;
    add_note(r, std::string(what) + " iteration did not converge in " + std::to_string(rounds) +
                    " rounds: no estimate (last round: Ka " + last(hasN, dN, satN) + ", Ks " + last(hasS, dS, satS) +
                    "; such values can differ between platforms)");
}

// Ka or Ks cannot be estimated because there is no site of that kind (e.g. a gene of Met and Trp codons only).
void no_sites(Result& r, const char* msg) {
    if (r.status != Status::error) r.status = Status::no_data;
    add_note(r, msg);
}

// ---------------------------------------------------------------------------------------------------------
// NG86 (KaKs_Calculator's NG): sites from positions 1 and 3 (changes to stop codons there are left out),
// pathways through stop codons excluded, Jukes-Cantor correction; S + N scaled to the compared length.
Result run_ng(const Prepared& p, const GeneticCode& g, Result r) {
    double S = 0, N = 0, Sd = 0, Nd = 0, snp = 0;
    auto sites = [&](int c) {
        int stop = 0;
        double syn = 0;
        for (int pos = 0; pos < 3; pos += 2)
            for (int b = 0; b < 4; ++b) {
                if (b == base_at(c, pos)) continue;
                int t = with_base(c, pos, b);
                if (is_stop(g, t))
                    ++stop;
                else if (g.aa[t] == g.aa[c])
                    syn += 1;
            }
        S += syn / 3.0;
        N += 3 - stop / 3.0 - syn / 3.0;
    };
    const size_t n = p.c1.size();
    for (size_t h = 0; h < n; ++h) {
        const int a = p.c1[h], b = p.c2[h];
        sites(a);
        sites(b);
        if (a == b) continue;
        int diff[3], num = 0;
        for (int i = 0; i < 3; ++i)
            if (base_at(a, i) != base_at(b, i)) diff[num++] = i;
        snp += num;
        double sdt = 0, ndt = 0;
        int stop = 0, path = num == 1 ? 1 : num == 2 ? 2 : 6;
        auto step = [&](int x, int y) {
            if (g.aa[x] == g.aa[y])
                sdt += 1;
            else
                ndt += 1;
        };
        if (num == 1) {
            step(a, b);
        } else if (num == 2) {
            for (int i = 0; i < 2; ++i) {
                int t1 = with_base(a, diff[i], base_at(b, diff[i]));
                if (is_stop(g, t1)) {
                    ++stop;
                    continue;
                }
                step(a, t1);
                step(t1, b);
            }
        } else {
            static const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
            for (const auto& o : perm) {
                int t1 = with_base(a, diff[o[0]], base_at(b, diff[o[0]]));
                int t2 = with_base(t1, diff[o[1]], base_at(b, diff[o[1]]));
                if (is_stop(g, t1) || is_stop(g, t2)) {
                    ++stop;
                    continue;
                }
                step(a, t1);
                step(t1, t2);
                step(t2, b);
            }
        }
        if (path == stop) {  // every pathway passes a stop codon
            if (num == 2) {
                Sd += 0.5;
                Nd += 1.5;
            } else {
                Sd += 1.0;
                Nd += 2.0;
            }
        } else {
            Sd += sdt / (path - stop);
            Nd += ndt / (path - stop);
        }
    }
    S /= 2;
    N /= 2;
    const double y = (double)(3 * n) / (S + N);
    S *= y;
    N *= y;
    r.S = S;
    r.N = N;
    r.Sd = Sd;
    r.Nd = Nd;
    r.substitutions = snp;
    auto jc = [&](double d_, double sites_, double& se, const char* which) {
        if (!(sites_ > 0)) return NA;
        double pp = d_ / sites_;
        double x = 1 - 4 * pp / 3;
        if (!(x > 0)) {
            mark_saturated(r, which);
            return NA;
        }
        se = std::sqrt(pp * (1 - pp) / sites_) / x;
        return -0.75 * std::log(x);
    };
    if (!(S > 0)) no_sites(r, "no synonymous sites");
    if (!(N > 0)) no_sites(r, "no nonsynonymous sites");
    r.ks = jc(Sd, S, r.se_ks, "Ks");
    r.ka = jc(Nd, N, r.se_ka, "Ka");
    if (std::isfinite(r.ks) && std::isfinite(r.ka)) r.divergence = (S * r.ks + N * r.ka) / (S + N);
    return r;
}

// ---------------------------------------------------------------------------------------------------------
// LWL85 family (KaKs_Calculator's LWL, MLWL, LPB, MLPB). Sites are classed 0-, 2- or 4-fold degenerate
// (number of synonymous, non-stop changes: 0, 1-2, 3); differences are split into transitions (Si) and
// transversions (Vi) per class, averaged over the pathways that avoid stop codons; Kimura 2-parameter per class.
int fold_class(const GeneticCode& g, int c, int pos) {  // 0, 1 (2-fold), 2 (4-fold)
    int k = 0;
    for (int b = 0; b < 4; ++b) {
        if (b == base_at(c, pos)) continue;
        int t = with_base(c, pos, b);
        if (!is_stop(g, t) && g.aa[t] == g.aa[c]) ++k;
    }
    return k == 0 ? 0 : k == 3 ? 2 : 1;
}

// codon indices of the special cases KaKs_Calculator (following Li's kakstools.py, and Tzeng et al. 2004 for
// the M variants) hard-codes; they are standard-code codons and are applied with every genetic code, as there.
constexpr int cod(int a, int b, int c) { return a * 16 + b * 4 + c; }
constexpr int T_ = 0, C_ = 1, A_ = 2, G_ = 3;
constexpr int CGA = cod(C_, G_, A_), AGA = cod(A_, G_, A_), CGG = cod(C_, G_, G_), AGG = cod(A_, G_, G_);
constexpr int ATA = cod(A_, T_, A_), ATG = cod(A_, T_, G_), ATC = cod(A_, T_, C_), ATT = cod(A_, T_, T_);

struct LwlAcc {
    double L[3] = {0, 0, 0}, Si[3] = {0, 0, 0}, Vi[3] = {0, 0, 0};
    double snp = 0;
};

// One substitution x -> y at position pos, adding half a difference to the class of each end codon.
void lwl_step(const GeneticCode& g, bool modified, int x, int y, int pos, double* si, double* vi) {
    const int c1 = fold_class(g, x, pos), c2 = fold_class(g, y, pos);
    int syn = -1;
    if (!modified) {
        if (pos == 0 && ((x == CGA && y == AGA) || (x == AGA && y == CGA) || (x == CGG && y == AGG) ||
                         (x == AGG && y == CGG))) {
            si[c1] += 0.5;  // as KaKs_Calculator's LWL85: half a transition, half a transversion
            vi[c2] += 0.5;
            return;
        }
    } else {
        if (pos == 2 && ((x == ATA && y == ATG) || (x == ATG && y == ATA))) syn = 0;
        if (pos == 2 && ((x == ATA && (y == ATC || y == ATT)) || ((x == ATC || x == ATT) && y == ATA))) syn = 1;
        if (pos == 0 && ((x == CGA && y == AGA) || (x == AGA && y == CGA) || (x == CGG && y == AGG) ||
                         (x == AGG && y == CGG)))
            syn = 1;
    }
    if (syn < 0) syn = transition(base_at(x, pos), base_at(y, pos)) ? 1 : 0;
    double* t = syn ? si : vi;
    t[c1] += 0.5;
    t[c2] += 0.5;
}

LwlAcc lwl_count(const Prepared& p, const GeneticCode& g, bool modified) {
    LwlAcc acc;
    const size_t n = p.c1.size();
    for (size_t h = 0; h < n; ++h) {
        const int a = p.c1[h], b = p.c2[h];
        for (int i = 0; i < 3; ++i) {
            acc.L[fold_class(g, a, i)] += 1;
            acc.L[fold_class(g, b, i)] += 1;
        }
        if (a == b) continue;
        int diff[3], num = 0;
        for (int i = 0; i < 3; ++i)
            if (base_at(a, i) != base_at(b, i)) diff[num++] = i;
        acc.snp += num;
        double si[3] = {0, 0, 0}, vi[3] = {0, 0, 0};
        int stop = 0, paths = num == 1 ? 1 : num == 2 ? 2 : 6;
        if (num == 1) {
            lwl_step(g, modified, a, b, diff[0], si, vi);
        } else if (num == 2) {
            for (int i = 0; i < 2; ++i) {
                int j = 1 - i;
                int t1 = with_base(a, diff[i], base_at(b, diff[i]));
                if (is_stop(g, t1)) {
                    ++stop;
                    continue;
                }
                lwl_step(g, modified, a, t1, diff[i], si, vi);
                lwl_step(g, modified, t1, b, diff[j], si, vi);
            }
        } else {
            static const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
            for (const auto& o : perm) {
                int t1 = with_base(a, diff[o[0]], base_at(b, diff[o[0]]));
                int t2 = with_base(t1, diff[o[1]], base_at(b, diff[o[1]]));
                if (is_stop(g, t1) || is_stop(g, t2)) {
                    ++stop;
                    continue;
                }
                lwl_step(g, modified, a, t1, diff[o[0]], si, vi);
                lwl_step(g, modified, t1, t2, diff[o[1]], si, vi);
                lwl_step(g, modified, t2, b, diff[o[2]], si, vi);
            }
        }
        if (paths - stop > 0)
            for (int k = 0; k < 3; ++k) {
                acc.Si[k] += si[k] / (paths - stop);
                acc.Vi[k] += vi[k] / (paths - stop);
            }
    }
    for (double& l : acc.L) l /= 2;
    return acc;
}

Result run_lwl(const Prepared& p, const GeneticCode& g, Method m, Result r) {
    const bool modified = (m == Method::MLWL || m == Method::MLPB);
    LwlAcc acc = lwl_count(p, g, modified);
    double A[3] = {0, 0, 0}, B[3] = {0, 0, 0}, K[3] = {0, 0, 0};
    bool sat[3] = {false, false, false};
    double ts = 0, tv = 0;
    for (int k = 0; k < 3; ++k) {
        ts += acc.Si[k];
        tv += acc.Vi[k];
        if (!(acc.L[k] > 0)) continue;  // no site of this class: contributes nothing
        const double P = acc.Si[k] / acc.L[k], Q = acc.Vi[k] / acc.L[k];
        const double x = 1 - 2 * P - Q, y = 1 - 2 * Q;
        if (x > 0 && y > 0) {
            const double la = std::log(1 / x), lb = std::log(1 / y);
            B[k] = lb >= 0 ? 0.5 * lb : 0;
            const double av = 0.5 * la - 0.25 * lb;
            A[k] = av >= 0 ? av : 0;  // KaKs_Calculator truncates a negative A at 0
            K[k] = A[k] + B[k];
        } else {
            sat[k] = true;
        }
    }
    const double* L = acc.L;
    double kappa;
    if (m == Method::MLWL) {
        kappa = (ts < SMALL || tv < SMALL) ? 1.0 : 2.0 * ts / tv;
    } else {
        kappa = tv > 0 ? 2 * ts / tv : 2.0;
    }
    r.kappa = kappa;
    r.rates[0] = r.rates[1] = kappa;
    r.has_fold = true;
    for (int k = 0; k < 3; ++k) {
        r.fold_sites[k] = L[k];
        r.fold_s[k] = acc.Si[k];
        r.fold_n[k] = acc.Vi[k];
    }
    r.substitutions = acc.snp;
    const bool ks_sat = sat[1] || sat[2], ka_sat = sat[0] || sat[1];
    double Sd = L[1] * A[1] + L[2] * K[2];
    double Nd = L[0] * K[0] + L[1] * B[1];
    double S, N, ks, ka;
    if (m == Method::LWL || m == Method::MLWL) {
        if (m == Method::MLWL && kappa > 2.0) {
            S = (kappa - 1) * L[1] / (kappa + 1) + L[2];
            N = L[0] + 2 * L[1] / (kappa + 1);
        } else if (m == Method::MLWL && kappa > 0.5) {
            S = (kappa - 0.5) * L[1] / (kappa + 1.5) + L[2];
            N = L[0] + 2 * L[1] / (kappa + 1.5);
        } else {
            S = L[1] / 3 + L[2];
            N = L[0] + 2 * L[1] / 3;
        }
        ks = S > 0 ? Sd / S : NA;
        ka = N > 0 ? Nd / N : NA;
        r.S = S;
        r.N = N;
        if (std::isfinite(ks) && std::isfinite(ka) && !ks_sat && !ka_sat) r.divergence = (S * ks + N * ka) / (S + N);
    } else {  // LPB, MLPB
        ks = (L[1] + L[2]) > 0 ? B[2] + (L[1] * A[1] + L[2] * A[2]) / (L[1] + L[2]) : NA;
        ka = (L[0] + L[1]) > 0 ? A[0] + (L[0] * B[0] + L[1] * B[1]) / (L[0] + L[1]) : NA;
        r.S = (std::isfinite(ks) && ks > 0) ? Sd / ks : NA;  // KaKs_Calculator: S = Sd/Ks, N = Nd/Ka
        r.N = (std::isfinite(ka) && ka > 0) ? Nd / ka : NA;
        const double lt = L[0] + L[1] + L[2];
        if (lt > 0 && !sat[0] && !sat[1] && !sat[2]) r.divergence = (L[0] * K[0] + L[1] * K[1] + L[2] * K[2]) / lt;
    }
    r.Sd = Sd;
    r.Nd = Nd;
    if (ks_sat) {
        mark_saturated(r, "Ks");
        ks = NA;
    }
    if (ka_sat) {
        mark_saturated(r, "Ka");
        ka = NA;
    }
    r.ks = ks;
    r.ka = ka;
    if (!(L[1] + L[2] > 0)) no_sites(r, "no synonymous sites");
    if (!(L[0] + L[1] > 0)) no_sites(r, "no nonsynonymous sites");
    return r;
}

// ---------------------------------------------------------------------------------------------------------
// Symmetric eigen solver (Householder tridiagonalisation + implicit QL; the classical EISPACK tred2/tql2
// algorithm). V: the symmetric n x n matrix, replaced by the eigenvectors, one per ROW (V[k*n + i] = component i
// of eigenvector k); d: eigenvalues. The algorithm is written for eigenvectors in columns, v(i, j); it runs on the
// transposed storage so that its inner loops (over the first index) walk memory contiguously.
void eigen_sym(std::vector<double>& V, int n, std::vector<double>& d) {
    std::vector<double> e((size_t)n, 0.0);
    d.assign((size_t)n, 0.0);
    double* const Z = V.data();
    auto v = [Z, n](int i, int j) -> double& { return Z[(size_t)j * n + i]; };
    for (int j = 0; j < n; ++j) d[j] = v(n - 1, j);
    for (int i = n - 1; i > 0; --i) {
        double scale = 0, h = 0;
        for (int k = 0; k < i; ++k) scale += std::fabs(d[k]);
        if (scale == 0) {
            e[i] = d[i - 1];
            for (int j = 0; j < i; ++j) {
                d[j] = v(i - 1, j);
                v(i, j) = 0;
                v(j, i) = 0;
            }
        } else {
            for (int k = 0; k < i; ++k) {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            double f = d[i - 1], g = std::sqrt(h);
            if (f > 0) g = -g;
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;
            for (int j = 0; j < i; ++j) e[j] = 0;
            for (int j = 0; j < i; ++j) {
                f = d[j];
                v(j, i) = f;
                g = e[j] + v(j, j) * f;
                for (int k = j + 1; k <= i - 1; ++k) {
                    g += v(k, j) * d[k];
                    e[k] += v(k, j) * f;
                }
                e[j] = g;
            }
            f = 0;
            for (int j = 0; j < i; ++j) {
                e[j] /= h;
                f += e[j] * d[j];
            }
            double hh = f / (h + h);
            for (int j = 0; j < i; ++j) e[j] -= hh * d[j];
            for (int j = 0; j < i; ++j) {
                f = d[j];
                g = e[j];
                for (int k = j; k <= i - 1; ++k) v(k, j) -= (f * e[k] + g * d[k]);
                d[j] = v(i - 1, j);
                v(i, j) = 0;
            }
        }
        d[i] = h;
    }
    for (int i = 0; i < n - 1; ++i) {
        v(n - 1, i) = v(i, i);
        v(i, i) = 1;
        double h = d[i + 1];
        if (h != 0) {
            for (int k = 0; k <= i; ++k) d[k] = v(k, i + 1) / h;
            for (int j = 0; j <= i; ++j) {
                double g = 0;
                for (int k = 0; k <= i; ++k) g += v(k, i + 1) * v(k, j);
                for (int k = 0; k <= i; ++k) v(k, j) -= g * d[k];
            }
        }
        for (int k = 0; k <= i; ++k) v(k, i + 1) = 0;
    }
    for (int j = 0; j < n; ++j) {
        d[j] = v(n - 1, j);
        v(n - 1, j) = 0;
    }
    v(n - 1, n - 1) = 1;
    e[0] = 0;
    // tql2
    for (int i = 1; i < n; ++i) e[i - 1] = e[i];
    e[n - 1] = 0;
    double f = 0, tst1 = 0;
    const double eps = std::ldexp(1.0, -52);
    for (int l = 0; l < n; ++l) {
        tst1 = std::max(tst1, std::fabs(d[l]) + std::fabs(e[l]));
        int m = l;
        while (m < n) {
            if (std::fabs(e[m]) <= eps * tst1) break;
            ++m;
        }
        if (m == n) m = n - 1;
        if (m > l) {
            int iter = 0;
            do {
                ++iter;
                double g = d[l];
                double p = (d[l + 1] - g) / (2 * e[l]);
                double r = std::hypot(p, 1.0);
                if (p < 0) r = -r;
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                double dl1 = d[l + 1];
                double h = g - d[l];
                for (int i = l + 2; i < n; ++i) d[i] -= h;
                f += h;
                p = d[m];
                double c = 1, c2 = c, c3 = c, el1 = e[l + 1], s = 0, s2 = 0;
                for (int i = m - 1; i >= l; --i) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = std::hypot(p, e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    for (int k = 0; k < n; ++k) {
                        h = v(k, i + 1);
                        v(k, i + 1) = s * v(k, i) + c * h;
                        v(k, i) = c * v(k, i) - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
            } while (std::fabs(e[l]) > eps * tst1 && iter < 100);
        }
        d[l] += f;
        e[l] = 0;
    }
}

// Transition probabilities P(t) = exp(Q t / mr) between codons one substitution apart, for the pathway
// weighting of YN00 / MYN. Q: rate 1, x kappa (transition; MYN: kappa_tc or kappa_ag), x omega (nonsynonymous),
// x pi_j; stop codons excluded; mr = mean rate. Codons with pi = 0 keep P = identity (as the references).
struct PMat {
    double p[64][64];
};

void pmat_codon(const GeneticCode& g, const double pi[64], double kappa_tc, double kappa_ag, double omega,
                double t, PMat& out) {
    int idx[64], ns = 0, back[64];
    for (int c = 0; c < 64; ++c) {
        back[c] = -1;
        if (pi[c] > 0) {
            back[c] = ns;
            idx[ns++] = c;
        }
    }
    // rate scaling factor over all codons
    double mr = 0;
    std::vector<double> A((size_t)ns * ns, 0.0);
    double s[64][64];  // symmetric exchangeabilities
    for (int i = 0; i < 64; ++i)
        for (int j = 0; j < 64; ++j) s[i][j] = 0;
    for (int i = 0; i < 64; ++i) {
        if (is_stop(g, i)) continue;
        for (int j = 0; j < i; ++j) {
            if (is_stop(g, j)) continue;
            int nd = 0, pos = 0;
            for (int k = 0; k < 3; ++k)
                if (base_at(i, k) != base_at(j, k)) {
                    ++nd;
                    pos = k;
                }
            if (nd != 1) continue;
            const int x = base_at(i, pos), y = base_at(j, pos);
            double q = 1;
            if (x + y == 1)
                q *= kappa_tc;
            else if (x + y == 5)
                q *= kappa_ag;
            if (g.aa[i] != g.aa[j]) q *= omega;
            s[i][j] = s[j][i] = q;
        }
    }
    for (int i = 0; i < 64; ++i) {
        double row = 0;
        for (int j = 0; j < 64; ++j) row += s[i][j] * pi[j];
        mr += pi[i] * row;
    }
    for (int a = 0; a < ns; ++a) {
        const int i = idx[a];
        double row = 0;
        for (int j = 0; j < 64; ++j) row += s[i][j] * pi[j];
        A[(size_t)a * ns + a] = -row;
        for (int b = 0; b < a; ++b) {
            const int j = idx[b];
            const double v = s[i][j] * std::sqrt(pi[i] * pi[j]);
            A[(size_t)a * ns + b] = A[(size_t)b * ns + a] = v;
        }
    }
    std::vector<double> d;
    eigen_sym(A, ns, d);  // A: eigenvector k in row k
    // U[a][k] = R[a][k] * exp(t d_k / (2 mr)), so that sum_k R_ak R_bk exp(t d_k / mr) = U[a] . U[b]
    std::vector<double> U((size_t)ns * ns);
    for (int k = 0; k < ns; ++k) {
        const double h = std::exp(t * d[k] / (2 * mr));
        const double* rk = &A[(size_t)k * ns];
        for (int a = 0; a < ns; ++a) U[(size_t)a * ns + k] = rk[a] * h;
    }
    for (int i = 0; i < 64; ++i)
        for (int j = 0; j < 64; ++j) out.p[i][j] = (i == j && back[i] < 0) ? 1.0 : 0.0;
    // only the entries the pathway weighting reads: one substitution apart (and the diagonal, for completeness)
    for (int a = 0; a < ns; ++a) {
        const int i = idx[a];
        const double* ra = &U[(size_t)a * ns];
        for (int pos = 0; pos < 3; ++pos)
            for (int bb = 0; bb < 4; ++bb) {
                const int j = with_base(i, pos, bb);
                if (j == i && pos > 0) continue;
                const int b = back[j];
                if (b < 0) continue;
                const double* rb = &U[(size_t)b * ns];
                double sum = 0;
                for (int k = 0; k < ns; ++k) sum += ra[k] * rb[k];
                double v = sum * std::sqrt(pi[j] / pi[i]);
                out.p[i][j] = v < 0 ? 0 : v;
            }
    }
}

// Distance from the proportions of transitional (P) and transversional (Q) differences under F84, falling back
// to K80, then JC69 (Yang & Nielsen 2000 / PAML yn00's DistanceF84). Returns failure flags; saturated = the
// distance is not defined (P + Q > 1, or JC69 fails), and t is then the references' cap.
struct F84Out {
    double t = 0, se = 0, kappa = -1;
    bool saturated = false;
    int fails = 0;
};

F84Out distance_f84(double n, double P, double Q, const double pi4[4], double maxkappa_cap, double maxkappa_set) {
    F84Out o;
    const double Qsmall = std::min(1e-10, 0.1 / n), maxt = 99;
    const double Y = pi4[0] + pi4[1], R = pi4[2] + pi4[3], tc = pi4[0] * pi4[1], ag = pi4[2] * pi4[3];
    if (P + Q > 1) {
        o.t = maxt;
        o.kappa = 1;
        o.fails = 3;
        o.saturated = true;
        return o;
    }
    if (P < -1e-10 || Q < -1e-10 || std::fabs(Y + R - 1) > 1e-8) {
        o.fails = 3;
        o.t = NA;
        return o;
    }
    bool f84 = false, k80 = false, jc = false;
    double a = 0, b = 0, A = 0, B = 0, C = 0;
    if (Q < Qsmall) {
        f84 = k80 = true;
    } else if (Y <= 0 || R <= 0 || (tc <= 0 && ag <= 0)) {
        f84 = true;
    } else {
        A = tc / Y + ag / R;
        B = tc + ag;
        C = Y * R;
        a = (2 * B + 2 * (tc * R / Y + ag * Y / R) * (1 - Q / (2 * C)) - P) / (2 * A);
        b = 1 - Q / (2 * C);
        if (a <= 0 || b <= 0) f84 = true;
    }
    if (!f84) {
        a = -0.5 * std::log(a);
        b = -0.5 * std::log(b);
        if (b <= 0) {
            f84 = true;
        } else {
            const double kF84 = a / b - 1;
            o.t = 4 * b * (tc * (1 + kF84 / Y) + ag * (1 + kF84 / R) + C);
            o.kappa = (B + (tc / Y + ag / R) * kF84) / B;
            const double aa = A * C / (A * C - C * P / 2 - (A - B) * Q / 2);
            const double bb = A * (A - B) / (A * C - C * P / 2 - (A - B) * Q / 2) - (A - B - C) / (C - Q / 2);
            o.se = std::sqrt((aa * aa * P + bb * bb * Q - (aa * P + bb * Q) * (aa * P + bb * Q)) / n);
        }
    }
    if (f84 && !k80) {
        a = 1 - 2 * P - Q;
        b = 1 - 2 * Q;
        if (a <= 0 || b <= 0) {
            k80 = true;
        } else {
            a = -std::log(a);
            b = -std::log(b);
            if (b <= 0) {
                k80 = true;
            } else {
                o.kappa = (0.5 * a - 0.25 * b) / (0.25 * b);
                o.t = 0.5 * a + 0.25 * b;
            }
            const double aa = 1 / (1 - 2 * P - Q), bb = (aa + 1 / (1 - 2 * Q)) / 2;
            o.se = std::sqrt((aa * aa * P + bb * bb * Q - (aa * P + bb * Q) * (aa * P + bb * Q)) / n);
        }
    }
    if (k80) {
        P += Q;
        if (P >= 0.75) {
            jc = true;
            o.saturated = true;
            P = 0.75 * (n - 1.0) / n;
        }
        o.t = -0.75 * std::log(1 - P * 4 / 3.0);
        if (o.t > maxt) o.t = maxt;
        o.se = std::sqrt(9 * P * (1 - P) / n) / (3 - 4 * P);
    }
    if (o.kappa > maxkappa_cap) o.kappa = maxkappa_set;
    o.fails = (int)f84 + (int)k80 + (int)jc;
    return o;
}

// Codon frequencies from the position-specific base frequencies of both sequences (F3x4), stops removed.
void f3x4(const Prepared& p, const GeneticCode& g, double pi[64]) {
    double f[12] = {0};
    const size_t n = p.c1.size();
    for (size_t h = 0; h < n; ++h)
        for (int pos = 0; pos < 3; ++pos) {
            f[pos * 4 + base_at(p.c1[h], pos)] += 1;
            f[pos * 4 + base_at(p.c2[h], pos)] += 1;
        }
    for (double& x : f) x /= (double)(2 * n);
    double fstop = 0;
    for (int c = 0; c < 64; ++c) {
        pi[c] = f[base_at(c, 0)] * f[4 + base_at(c, 1)] * f[8 + base_at(c, 2)];
        if (is_stop(g, c)) {
            fstop += pi[c];
            pi[c] = 0;
        }
    }
    for (int c = 0; c < 64; ++c) pi[c] /= (1 - fstop);
}

// Base-pair counts at non-degenerate sites (F[0]) and at 4-fold degenerate third positions of codons coding
// the same amino acid (F[1]), symmetrised, for estimating kappa (YN00 GetKappa).
void kappa_tables(const Prepared& p, const GeneticCode& g, double F[2][16]) {
    for (int k = 0; k < 2; ++k)
        for (int j = 0; j < 16; ++j) F[k][j] = 0;
    const size_t n = p.c1.size();
    for (size_t h = 0; h < n; ++h) {
        const int c[2] = {p.c1[h], p.c2[h]};
        for (int pos = 0; pos < 3; ++pos) {
            int nondeg = 0;
            for (int k = 0; k < 2; ++k) {
                int b;
                for (b = 0; b < 4; ++b) {
                    if (b == base_at(c[k], pos)) continue;
                    if (g.aa[with_base(c[k], pos, b)] == g.aa[c[k]]) break;
                }
                if (b == 4) ++nondeg;
            }
            if (nondeg == 2) {
                const int x = base_at(c[0], pos), y = base_at(c[1], pos);
                F[0][x * 4 + y] += 0.5;
                F[0][y * 4 + x] += 0.5;
            }
        }
        int four = 0;
        for (int k = 0; k < 2; ++k) {
            int b;
            for (b = 0; b < 4; ++b)
                if (b != base_at(c[k], 2) && g.aa[with_base(c[k], 2, b)] != g.aa[c[k]]) break;
            if (g.aa[c[0]] == g.aa[c[1]] && b == 4) ++four;
        }
        if (four == 2) {
            const int x = base_at(c[0], 2), y = base_at(c[1], 2);
            F[1][x * 4 + y] += 0.5;
            F[1][y * 4 + x] += 0.5;
        }
    }
}

// Synonymous / nonsynonymous sites of one sequence weighted by kappa(s) and pi (YN00 CountSites), and the
// base frequencies at those sites; S + N scaled to the sequence length.
void count_sites(const std::vector<uint8_t>& z, const GeneticCode& g, const double pi[64], double ktc, double kag,
                 double& St, double& Nt, double fbS[4], double fbN[4]) {
    St = Nt = 0;
    for (int k = 0; k < 4; ++k) fbS[k] = fbN[k] = 0;
    for (uint8_t c0 : z) {
        double S = 0, N = 0;
        for (int j = 0; j < 3; ++j) {
            const int bj = base_at(c0, j);
            for (int k = 0; k < 4; ++k) {
                if (k == bj) continue;
                const int c1 = with_base(c0, j, k);
                if (is_stop(g, c1)) continue;
                double r = pi[c1];
                if (k + bj == 1)
                    r *= ktc;
                else if (k + bj == 5)
                    r *= kag;
                if (g.aa[c0] == g.aa[c1]) {
                    S += r;
                    fbS[bj] += r;
                } else {
                    N += r;
                    fbN[bj] += r;
                }
            }
        }
        St += S;
        Nt += N;
    }
    const double r = (double)(3 * z.size()) / (St + Nt);
    St *= r;
    Nt *= r;
    double sS = 0, sN = 0;
    for (int k = 0; k < 4; ++k) {
        sS += fbS[k];
        sN += fbN[k];
    }
    for (int k = 0; k < 4; ++k) {
        fbS[k] /= sS;
        fbN[k] /= sN;
    }
}

// Sites of the pair: the mean of the two sequences (YN00 DistanceYN00). The base frequencies at the synonymous
// (nonsynonymous) sites are the mean over the sequences that have such sites: a sequence with none has no
// frequencies (0/0), and the references average that NaN in, which loses the estimate whenever the F84 or TN93
// formula then needs the frequencies.
void pair_sites(const Prepared& p, const GeneticCode& g, const double pi[64], double ktc, double kag, double& S,
                double& N, double fbS[4], double fbN[4]) {
    S = N = 0;
    int withS = 0, withN = 0;
    for (int j = 0; j < 4; ++j) fbS[j] = fbN[j] = 0;
    for (int s = 0; s < 2; ++s) {
        double St, Nt, a[4], b[4];
        count_sites(s == 0 ? p.c1 : p.c2, g, pi, ktc, kag, St, Nt, a, b);
        S += St / 2;
        N += Nt / 2;
        if (St > 0) {
            ++withS;
            for (int j = 0; j < 4; ++j) fbS[j] += a[j] / 2;
        }
        if (Nt > 0) {
            ++withN;
            for (int j = 0; j < 4; ++j) fbN[j] += b[j] / 2;
        }
    }
    for (int j = 0; j < 4; ++j) {
        if (withS == 1) fbS[j] *= 2;
        if (withN == 1) fbN[j] *= 2;
    }
}

// Synonymous / nonsynonymous differences split into transitions (ts1: A<->G, ts2: T<->C) and transversions,
// pathways weighted by P (nullptr: equal weights), pathways through stop codons excluded (YN00 CountDiffs).
struct Diffs {
    double sts1 = 0, sts2 = 0, stv = 0, nts1 = 0, nts2 = 0, ntv = 0, snp = 0;
};

Diffs count_diffs(const Prepared& p, const GeneticCode& g, const PMat* P, bool myn) {
    Diffs D;
    const size_t n = p.c1.size();
    for (size_t h = 0; h < n; ++h) {
        const int c0 = p.c1[h], c1 = p.c2[h];
        if (c0 == c1) continue;
        int dmark[3], nd = 0;
        for (int k = 0; k < 3; ++k)
            if (base_at(c0, k) != base_at(c1, k)) dmark[nd++] = k;
        D.snp += nd;
        double s1 = 0, s2 = 0, sv = 0, n1 = 0, n2 = 0, nv = 0;
        auto classify = [](int x, int y, bool syn, double w, double& a1, double& a2, double& av, double& b1,
                           double& b2, double& bv) {
            const int sum = x + y;
            if (syn) {
                if (sum == 5)
                    a1 += w;
                else if (sum == 1)
                    a2 += w;
                else
                    av += w;
            } else {
                if (sum == 5)
                    b1 += w;
                else if (sum == 1)
                    b2 += w;
                else
                    bv += w;
            }
        };
        if (nd == 1) {
            const int pos = dmark[0];
            classify(base_at(c0, pos), base_at(c1, pos), g.aa[c0] == g.aa[c1], 1.0, s1, s2, sv, n1, n2, nv);
        } else {
            const int npath = nd == 2 ? 2 : 6;
            double pp[6], ps[6][6];
            int nstop = 0;
            for (int k = 0; k < npath; ++k) {
                int step[3];
                if (nd == 2) {
                    step[0] = dmark[k];
                    step[1] = dmark[1 - k];
                } else {
                    step[0] = k / 2;
                    step[1] = k % 2;
                    if (step[0] <= step[1]) step[1]++;
                    step[2] = 3 - step[0] - step[1];
                }
                for (double& x : ps[k]) x = 0;
                pp[k] = 1;
                int cur = c0;
                for (int i = 0; i < nd; ++i) {
                    const int pos = step[i];
                    const int nxt = with_base(cur, pos, base_at(c1, pos));
                    if (P) pp[k] *= P->p[cur][nxt];
                    if (is_stop(g, nxt)) {
                        ++nstop;
                        pp[k] = 0;
                        break;
                    }
                    classify(base_at(c0, pos), base_at(c1, pos), g.aa[cur] == g.aa[nxt], 1.0, ps[k][0], ps[k][1],
                             ps[k][2], ps[k][3], ps[k][4], ps[k][5]);
                    cur = nxt;
                }
            }
            if (nstop == npath) {  // every pathway passes a stop codon
                if (myn) {
                    n1 = 0.25;
                    n2 = 0.25;
                } else {
                    n1 = 0.5;  // (counted as nonsynonymous transitions, as the references)
                }
                nv = nd == 2 ? 1.5 : 2.5;
            } else {
                double sump = 0;
                for (int k = 0; k < npath; ++k) sump += pp[k];
                if (sump > 1e-20)
                    for (int k = 0; k < npath; ++k) {
                        const double w = pp[k] / sump;
                        s1 += ps[k][0] * w;
                        s2 += ps[k][1] * w;
                        sv += ps[k][2] * w;
                        n1 += ps[k][3] * w;
                        n2 += ps[k][4] * w;
                        nv += ps[k][5] * w;
                    }
            }
        }
        D.sts1 += s1;
        D.sts2 += s2;
        D.stv += sv;
        D.nts1 += n1;
        D.nts2 += n2;
        D.ntv += nv;
    }
    return D;
}

Result run_yn(const Prepared& p, const GeneticCode& g, const Options& opt, Result r) {
    const bool kc = opt.yn_flavour == Flavour::kaks_calculator;
    const double cap = kc ? 99 : 999, capset = kc ? 2 : 999;
    double pi[64];
    f3x4(p, g, pi);
    // kappa from non-degenerate and 4-fold sites, weighted by their numbers
    double F[2][16];
    kappa_tables(p, g, F);
    double ka[2] = {0, 0}, wk[2] = {0, 0};
    double t = 0.4;  // both references start the iteration at t = 0.4, omega = 1
    for (int k = 0; k < 2; ++k) {
        double S = 0;
        for (double x : F[k]) S += x;
        if (S <= 0) continue;
        for (double& x : F[k]) x /= S;
        const double T = (F[k][0 * 4 + 1] + F[k][2 * 4 + 3]) * 2;
        const double V = 1 - T - (F[k][0] + F[k][5] + F[k][10] + F[k][15]);
        double pi4[4];
        for (int j = 0; j < 4; ++j) pi4[j] = F[k][j * 4] + F[k][j * 4 + 1] + F[k][j * 4 + 2] + F[k][j * 4 + 3];
        F84Out o = distance_f84(S, T, V, pi4, cap, capset);
        ka[k] = o.kappa;
        wk[k] = o.kappa > 0 ? S : 0;
    }
    const double kappa = (wk[0] + wk[1] == 0) ? (kc ? 2.0 : 4.6) : (ka[0] * wk[0] + ka[1] * wk[1]) / (wk[0] + wk[1]);
    r.kappa = kappa;
    r.rates[0] = r.rates[1] = kappa;

    double S = 0, N = 0, fbS[4], fbN[4];
    pair_sites(p, g, pi, kappa, kappa, S, N, fbS, fbN);
    r.S = S;
    r.N = N;
    // No site of one kind (S = 0: e.g. Met/Trp codons only; N = 0: e.g. GCT against GCA, where F3x4 gives weight
    // to synonymous neighbours only): that estimate is NA, and the other one is still estimated, from its own
    // sites alone. No site of either kind: nothing to estimate.
    const bool hasS = S > 0, hasN = N > 0;
    if (!hasS) no_sites(r, "no synonymous sites");
    if (!hasN) no_sites(r, "no nonsynonymous sites");
    if (!hasS && !hasN) return r;
    double omega = 1, dS = 0.1, dN = 0.1, w0 = 0, dS0 = 0, dN0 = 0;
    if (!hasS) dS = 0;
    if (!hasN) dN = 0;
    const double accu = 5e-4, minomega = 1e-5, maxomega = 99;
    const int nround = opt.yn_weighting ? 10 : 1;
    if (t == 0) t = 0.5;
    if (opt.yn_weighting) {
        if (t < 0.001 || t > 5) t = 0.5;
        if (omega < 0.01 || omega > 5) omega = 0.5;
    }
    PMat* P = nullptr;
    std::vector<PMat> store(opt.yn_weighting ? 1 : 0);
    if (opt.yn_weighting) P = &store[0];
    F84Out os, on;
    Diffs D;
    int ir;
    for (ir = 0; ir < nround; ++ir) {
        if (P) pmat_codon(g, pi, kappa, kappa, omega, t, *P);
        D = count_diffs(p, g, P, false);
        const double Sdts = D.sts1 + D.sts2, Sdtv = D.stv, Ndts = D.nts1 + D.nts2, Ndtv = D.ntv;
        if (hasS) os = distance_f84(S, Sdts / S, Sdtv / S, fbS, cap, capset);
        if (hasN) on = distance_f84(N, Ndts / N, Ndtv / N, fbN, cap, capset);
        if (hasS && (os.fails != 3 || os.saturated)) dS = os.t;
        if (hasN && (on.fails != 3 || on.saturated)) dN = on.t;
        // With one kind of site missing, omega has no meaning; the references' own arithmetic then leaves it at
        // minomega (max(minomega, NaN)), and it has no effect: every change between codons of non-zero frequency
        // is of the one kind, so omega cancels in the scaled rate matrix. t comes from the kind that exists
        // (the references add NaN here and lose the pair when a codon differs at two or three positions).
        if (!hasS || !hasN)
            omega = minomega;
        else
            omega = dS < 1e-9 ? maxomega : std::max(minomega, dN / dS);
        t = dS * 3 * S / (S + N) + dN * 3 * N / (S + N);
        if (std::fabs(dS - dS0) < accu && std::fabs(dN - dN0) < accu && std::fabs(omega - w0) < accu) break;
        dS0 = dS;
        dN0 = dN;
        w0 = omega;
    }
    r.iterations = std::min(ir + 1, nround);
    r.converged = !(opt.yn_weighting && ir == nround);
    r.Sd = D.sts1 + D.sts2 + D.stv;
    r.Nd = D.nts1 + D.nts2 + D.ntv;
    r.substitutions = D.snp;
    if (!r.converged) {
        not_converged(r, "YN00", nround, hasN, dN, hasN && on.saturated, hasS, dS, hasS && os.saturated);
        return r;
    }
    r.ks = hasS ? dS : NA;
    r.ka = hasN ? dN : NA;
    r.se_ks = hasS ? os.se : NA;
    r.se_ka = hasN ? on.se : NA;
    if (hasS && os.saturated) {
        mark_saturated(r, "Ks");
        r.ks = r.se_ks = NA;
    }
    if (hasN && on.saturated) {
        mark_saturated(r, "Ka");
        r.ka = r.se_ka = NA;
    }
    if (std::isfinite(r.ks) && std::isfinite(r.ka)) r.divergence = (S * r.ks + N * r.ka) / (S + N);
    return r;
}

// ---------------------------------------------------------------------------------------------------------
// MYN (Zhang, Li & Yu 2006): YN00 with separate transition rates between purines (kappa_ag) and between
// pyrimidines (kappa_tc), and Tamura-Nei (TN93) distances; falls back to F84/K80/JC69 where TN93 fails.
void kappa_tn93(double n, double P1, double P2, double Q, const double pi4[4], double& ktc, double& kag) {
    const double Qsmall = std::min(1e-10, 0.1 / n);
    ktc = kag = -1;
    const double Y = pi4[0] + pi4[1], R = pi4[2] + pi4[3], tc = pi4[0] * pi4[1], ag = pi4[2] * pi4[3];
    if ((P1 + P2 + Q) > 1 || (P1 + P2) < -1e-10 || Q < -1e-10 || std::fabs(Y + R - 1) > 1e-8) return;
    bool fail = false;
    if (Q < Qsmall) {
        fail = true;
    } else if (Y <= 0 || R <= 0 || (tc <= 0 && ag <= 0)) {
        fail = true;
    } else {
        double a1 = 1 - R * P1 / (2 * ag) - Q / (2 * R);
        double a2 = 1 - Y * P2 / (2 * tc) - Q / (2 * Y);
        double b = 1 - Q / (2 * Y * R);
        if (a1 < 0 || a2 < 0 || b < 0) {
            fail = true;
        } else {
            a1 = std::log(a1);
            a2 = std::log(a2);
            b = std::log(b);
            kag = (Y * b - a1) / (-R * b);
            ktc = (R * b - a2) / (-Y * b);
        }
    }
    if (fail) ktc = kag = 2;
}

F84Out distance_tn93(double n, double P1, double P2, double Q, const double pi4[4]) {
    const double Y = pi4[0] + pi4[1], R = pi4[2] + pi4[3], tc = pi4[0] * pi4[1], ag = pi4[2] * pi4[3];
    bool fail = false;
    F84Out o;
    // tc <= 0 or ag <= 0 (the sites have no T and C together, or no A and G together): that kind of transition
    // cannot occur, its term 2*tc*log(a2)/Y (2*ag*log(a1)/R) is zero, and TN93 reduces exactly to F84, which the
    // fallback below computes. The reference tests only "tc <= 0 && ag <= 0", then divides P by zero: 0/0 = NaN
    // when there is no such difference (and its Ka or Ks is nan), -inf < 0 = "TN93 fails" when there is one.
    if (P1 + P2 + Q > 1 || std::fabs(Y + R - 1) > 1e-10 || Y <= 0 || R <= 0 || tc <= 0 || ag <= 0) {
        fail = true;
    } else {
        const double C = Y * R;
        double a1 = 1 - R * P1 / (2 * ag) - Q / (2 * R);
        double a2 = 1 - Y * P2 / (2 * tc) - Q / (2 * Y);
        double b = 1 - Q / (2 * C);
        if (a1 < 0 || a2 < 0 || b < 0) {
            fail = true;
        } else {
            const double la1 = std::log(a1), la2 = std::log(a2), lb = std::log(b);
            o.t = (-2 * ag * la1 / R) + (-2 * tc * la2 / Y) + (-2 * (C - ag * Y / R - tc * R / Y) * lb);
            const double cc1 = 2 * ag * R / (2 * ag * R - R * R * P1 - ag * Q);
            const double cc2 = 2 * tc * Y / (2 * tc * Y - Y * Y * P2 - tc * Q);
            double cc3 = 2 * ag * ag / (R * (2 * ag * R - R * R * P1 - ag * Q));
            cc3 += 2 * tc * tc / (Y * (2 * tc * Y - Y * Y * P2 - tc * Q));
            cc3 += (R * R * (Y * Y - 2 * tc) + Y * Y * (R * R - 2 * ag)) / (2 * R * R * Y * Y - R * Y * Q);
            const double var = (cc1 * cc1 * P1 + cc2 * cc2 * P2 + cc3 * cc3 * Q -
                                (cc1 * P1 + cc2 * P2 + cc3 * Q) * (cc1 * P1 + cc2 * P2 + cc3 * Q)) /
                               n;
            o.se = var >= 0 ? std::sqrt(var) : NA;  // (KaKs_Calculator keeps the variance; not printed there)
            if (!std::isfinite(o.t)) o.saturated = true;
        }
    }
    if (fail) o = distance_f84(n, P1 + P2, Q, pi4, 99, 2);
    return o;
}

Result run_myn(const Prepared& p, const GeneticCode& g, Result r) {
    double pi[64];
    f3x4(p, g, pi);
    double F[2][16];
    kappa_tables(p, g, F);
    double ktcK[2] = {0, 0}, kagK[2] = {0, 0}, kK[2] = {0, 0}, wk[2] = {0, 0};
    for (int k = 0; k < 2; ++k) {
        double S = 0;
        for (double x : F[k]) S += x;
        if (S <= 0) continue;
        for (double& x : F[k]) x /= S;
        const double T1 = 2 * F[k][2 * 4 + 3], T2 = 2 * F[k][0 * 4 + 1];
        const double V = 1 - T1 - T2 - (F[k][0] + F[k][5] + F[k][10] + F[k][15]);
        double pi4[4];
        for (int j = 0; j < 4; ++j) pi4[j] = F[k][j * 4] + F[k][j * 4 + 1] + F[k][j * 4 + 2] + F[k][j * 4 + 3];
        kappa_tn93(S, T1, T2, V, pi4, ktcK[k], kagK[k]);
        wk[k] = (ktcK[k] > 0 && kagK[k] > 0) ? S : 0;
        kK[k] = 2 * (ktcK[k] * pi4[0] * pi4[1] + kagK[k] * pi4[2] * pi4[3]) / ((pi4[0] + pi4[1]) * (pi4[2] + pi4[3]));
    }
    double ktc, kag, kappa;
    if (wk[0] + wk[1] == 0) {
        ktc = kag = kappa = 2;
    } else {
        const double w = wk[0] + wk[1];
        auto mix = [&](const double* v) {
            double s = 0;
            for (int k = 0; k < 2; ++k)
                if (wk[k] != 0) s += v[k] * wk[k];
            return s / w;
        };
        ktc = mix(ktcK);
        kag = mix(kagK);
        kappa = mix(kK);
    }
    r.kappa = kappa;
    r.rates[0] = ktc;
    r.rates[1] = kag;
    double S = 0, N = 0, fbS[4], fbN[4];
    pair_sites(p, g, pi, ktc, kag, S, N, fbS, fbN);
    r.S = S;
    r.N = N;
    const bool hasS = S > 0, hasN = N > 0;  // as in run_yn: the estimate that has sites is still made
    if (!hasS) no_sites(r, "no synonymous sites");
    if (!hasN) no_sites(r, "no nonsynonymous sites");
    if (!hasS && !hasN) return r;
    double t = 0.09, omega = 0.5, dS = 0, dN = 0, dS0 = 0, dN0 = 0, w0 = 0;
    const double accu = 5e-8, minomega = 1e-5, maxomega = 99;
    const int nround = 100;
    PMat P;
    F84Out os, on;
    Diffs D;
    int ir;
    for (ir = 0; ir < nround; ++ir) {
        pmat_codon(g, pi, ktc, kag, omega, t, P);
        D = count_diffs(p, g, &P, true);
        double Sd = D.sts1 + D.sts2 + D.stv, Nd = D.nts1 + D.nts2 + D.ntv;
        double s1 = D.sts1, s2 = D.sts2, sv = D.stv, n1 = D.nts1, n2 = D.nts2, nv = D.ntv;
        if (hasS && Sd > S) {
            s1 *= S / Sd;
            s2 *= S / Sd;
            sv *= S / Sd;
        }
        if (hasN && Nd > N) {
            n1 *= N / Nd;
            n2 *= N / Nd;
            nv *= N / Nd;
        }
        if (hasS) os = distance_tn93(S, s1 / S, s2 / S, sv / S, fbS);
        if (hasN) on = distance_tn93(N, n1 / N, n2 / N, nv / N, fbN);
        dS = hasS ? os.t : 0;
        dN = hasN ? on.t : 0;
        if (!hasS || !hasN)
            omega = minomega;  // see run_yn
        else
            omega = dS < 1e-9 ? maxomega : std::max(minomega, dN / dS);
        t = dS * 3 * S / (S + N) + dN * 3 * N / (S + N);
        if (std::fabs(dS - dS0) < accu && std::fabs(dN - dN0) < accu && std::fabs(omega - w0) < accu) break;
        dS0 = dS;
        dN0 = dN;
        w0 = omega;
    }
    r.iterations = std::min(ir + 1, nround);
    r.converged = ir < nround;
    r.Sd = D.sts1 + D.sts2 + D.stv;
    r.Nd = D.nts1 + D.nts2 + D.ntv;
    r.substitutions = D.snp;
    if (!r.converged) {
        not_converged(r, "MYN", nround, hasN, dN, hasN && on.saturated, hasS, dS, hasS && os.saturated);
        return r;
    }
    r.ks = hasS ? dS : NA;
    r.ka = hasN ? dN : NA;
    r.se_ks = hasS ? os.se : NA;
    r.se_ka = hasN ? on.se : NA;
    if (hasS && (os.saturated || !std::isfinite(dS))) {
        mark_saturated(r, "Ks");
        r.ks = r.se_ks = NA;
    }
    if (hasN && (on.saturated || !std::isfinite(dN))) {
        mark_saturated(r, "Ka");
        r.ka = r.se_ka = NA;
    }
    if (std::isfinite(r.ks) && std::isfinite(r.ka)) r.divergence = (S * r.ks + N * r.ka) / (S + N);
    return r;
}

}  // namespace

const char* method_name(Method m) {
    switch (m) {
        case Method::NG: return "NG";
        case Method::LWL: return "LWL";
        case Method::MLWL: return "MLWL";
        case Method::LPB: return "LPB";
        case Method::MLPB: return "MLPB";
        case Method::YN: return "YN";
        case Method::MYN: return "MYN";
    }
    return "?";
}

const std::vector<Method>& output_order() {
    static const std::vector<Method> o = {Method::NG,   Method::LWL, Method::MLWL, Method::LPB,
                                          Method::MLPB, Method::YN,  Method::MYN};
    return o;
}

bool parse_method(const std::string& name, Method& out, std::string& why) {
    std::string u;
    for (char ch : name) u += (char)std::toupper((unsigned char)ch);
    for (Method m : output_order())
        if (u == method_name(m)) {
            out = m;
            return true;
        }
    static const char* const known[] = {"GY",   "MS",    "MA",    "NONE", "GNG", "GLWL",
                                        "GLPB", "GMLWL", "GMLPB", "GYN",  "GMYN", "ALL"};
    for (const char* k : known)
        if (u == k) {
            why = "method " + u + " is a KaKs_Calculator method that COAL does not implement natively "
                  "(available: NG, LWL, MLWL, LPB, MLPB, YN, MYN)";
            return false;
        }
    why = "unknown method '" + name + "' (available: NG, LWL, MLWL, LPB, MLPB, YN, MYN)";
    return false;
}

const char* status_name(Status s) {
    switch (s) {
        case Status::ok: return "ok";
        case Status::warning: return "warning";
        case Status::not_converged: return "not_converged";
        case Status::saturated: return "saturated";
        case Status::no_data: return "no_data";
        case Status::error: return "error";
    }
    return "?";
}

Prepared prepare(const std::string& row1, const std::string& row2, const GeneticCode& code) {
    Prepared p;
    Result& r = p.base;
    // GC content as KaKs_Calculator reports it: over both rows as given (gaps count in the length)
    {
        const std::string* rows[2] = {&row1, &row2};
        double gcp[3] = {0, 0, 0};
        size_t total = row1.size() + row2.size();
        for (const std::string* s : rows)
            for (size_t i = 0; i < s->size(); ++i) {
                const char c = (char)std::toupper((unsigned char)(*s)[i]);
                if (c == 'G' || c == 'C') gcp[i % 3] += 1;
            }
        if (total > 0 && total % 3 == 0 && row1.size() % 3 == 0) {
            r.gc[0] = (gcp[0] + gcp[1] + gcp[2]) / (double)total;
            for (int k = 0; k < 3; ++k) r.gc[k + 1] = gcp[k] / (double)(total / 3);
        }
    }
    if (row1.size() != row2.size()) {
        r.status = Status::error;
        r.note = "aligned sequences differ in length (" + std::to_string(row1.size()) + " and " +
                 std::to_string(row2.size()) + ")";
        return p;
    }
    if (row1.size() % 3 != 0) {
        r.status = Status::error;
        r.note = "alignment length " + std::to_string(row1.size()) + " is not a multiple of 3";
        return p;
    }
    const size_t nc = row1.size() / 3;
    // per column: -1 gap, -2 ambiguous, else codon index
    auto codon_at = [](const std::string& s, size_t col) {
        int idx = 0;
        bool gap = false, amb = false;
        for (int k = 0; k < 3; ++k) {
            char c = (char)std::toupper((unsigned char)s[col * 3 + k]);
            if (c == 'U') c = 'T';
            const int b = base_index((unsigned char)c);
            if (c == '-' || c == '.')
                gap = true;
            else if (b < 0)
                amb = true;
            else
                idx = idx * 4 + b;
        }
        return gap ? -1 : amb ? -2 : idx;
    };
    // terminal stop: the last codon of a sequence that is not a gap column, if it is a stop codon
    auto last_codon_col = [&](const std::string& s) -> long {
        for (size_t c = nc; c-- > 0;) {
            bool allgap = true;
            for (int k = 0; k < 3; ++k) {
                const char ch = s[c * 3 + k];
                if (ch != '-' && ch != '.') allgap = false;
            }
            if (!allgap) return (long)c;
        }
        return -1;
    };
    long term[2] = {-1, -1};
    for (int s = 0; s < 2; ++s) {
        const std::string& row = s == 0 ? row1 : row2;
        long c = last_codon_col(row);
        if (c >= 0) {
            int ci = codon_at(row, (size_t)c);
            if (ci >= 0 && code.aa[ci] == '*') term[s] = c;
        }
    }
    p.c1.reserve(nc);
    p.c2.reserve(nc);
    for (size_t c = 0; c < nc; ++c) {
        const int a = codon_at(row1, c), b = codon_at(row2, c);
        if ((long)c == term[0] || (long)c == term[1]) {
            r.terminal_stop_dropped = true;
            continue;
        }
        if (a == -1 || b == -1) {
            ++r.gap_codons;
            continue;
        }
        if (a == -2 || b == -2) {
            ++r.ambiguous_codons;
            continue;
        }
        if (code.aa[a] == '*' || code.aa[b] == '*') {
            ++r.internal_stops;
            continue;
        }
        p.c1.push_back((uint8_t)a);
        p.c2.push_back((uint8_t)b);
    }
    r.codons_compared = (long)p.c1.size();
    r.length = 3 * r.codons_compared;
    if (r.internal_stops > 0) {
        r.status = Status::warning;
        r.note = std::to_string(r.internal_stops) + " codon column(s) with an internal stop codon skipped";
    }
    if (r.codons_compared == 0) {
        r.status = Status::no_data;
        add_note(r, "no codon column left to compare (gaps, ambiguous bases or stop codons)");
    }
    return p;
}

Result compute(const Prepared& p, Method m, const GeneticCode& code, const Options& opt) {
    Result r = p.base;
    r.method = m;
    if (r.status == Status::error || r.status == Status::no_data) return r;
    switch (m) {
        case Method::NG: r = run_ng(p, code, r); break;
        case Method::LWL:
        case Method::MLWL:
        case Method::LPB:
        case Method::MLPB: r = run_lwl(p, code, m, r); break;
        case Method::YN: r = run_yn(p, code, opt, r); break;
        case Method::MYN: r = run_myn(p, code, r); break;
    }
    finish(r);
    return r;
}

Result compute(const std::string& row1, const std::string& row2, Method m, const GeneticCode& code,
               const Options& opt) {
    return compute(prepare(row1, row2, code), m, code, opt);
}

}  // namespace kaks
}  // namespace coal
