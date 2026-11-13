// SPDX-License-Identifier: GPL-3.0-only
//
// Ka/Ks engine (kaks.hpp): the method of Comeron (1995), Method::Comeron.
//
// This file is derived from libsequence as the R package orthologr compiles it into its function gestimator()
// (orthologr 0.4.3, https://github.com/drostlab/orthologr, src/: Comeron95.h, Sites.h, RedundancyCom95.h,
// SingleSub.h, TwoSubs.h, ThreeSubs.h, GranthamWeights.h, Grantham.h, Kimura80.h, Comparisons.h, Sequence.h and
// gestimator.cpp). Each of those files carries this notice:
//
//   Copyright (C) 2003-2009 Kevin Thornton, krthornt[]@[]uci.edu
//   This file is part of libsequence.
//   libsequence is free software: you can redistribute it and/or modify it under the terms of the GNU General
//   Public License as published by the Free Software Foundation, either version 3 of the License, or (at your
//   option) any later version.
//   [...]
//   Modified by Sarah Scharfenberg and Hajk-Georg Drost 2014 to work in orthologr without using external
//   libraries from libsequence.
//   All changes are also free under the terms of GNU General Public License version 2 of the License, or any
//   later version.
//
// This file is distributed under the GNU General Public License, version 3 (COAL as a whole is GPL-3.0-only; see
// LICENSE and NOTICE), which both statements allow. It is a modified version of that code, not an independent
// implementation: the functions follow the upstream ones in structure, constants and order of floating-point
// operations, so that the numbers are the same.
//
// Modified for COAL by Chloe Newton and Keanu Murati, 2026. What was changed:
//   - rewritten as free functions on codon indices 0..63 (COAL's order T C A G) instead of classes on strings;
//     the degeneracy tables are built once, into an immutable function-local static (no mutable global state),
//     so pairs can be computed on several threads; no Rcpp, no output;
//   - the Grantham distances are a 21 x 21 table with the upstream numbers (upstream: a 60 x 60 array indexed by
//     the letter, of which only these 441 entries are set);
//   - ComeronFlavour::coal (the default): COAL's policy is applied to the pair first (prepare(): a terminal stop
//     codon is dropped, columns with an internal stop are skipped, U is read as T), and where a logarithm of the
//     correction has an argument <= 0 the distance concerned is NA with status `saturated` (upstream replaces
//     that logarithm by 0 and prints the finite number that results); upstream's rule "Ba = 0 if Kimura's
//     distance over the whole alignment is below 1" is not applied;
//   - ComeronFlavour::gestimator: upstream's behaviour throughout (stop codons compared as codons, U not read
//     as T, undefined logarithms replaced by 0, the Kimura rule); upstream's 999 ("cannot be calculated") is NA;
//   - only maxHits = 3 (gestimator's default, the one orthologr uses) is implemented.
// One property of the upstream code is kept on purpose in both flavours, because the task of this method is to
// give gestimator's numbers: for codons that differ at two positions, TwoSubs calls its Calculate() with the
// arguments (codon1, intermediate1, codon2, intermediate2) where the parameters are (codon1, intermediate1,
// intermediate2, codon2), so the four branches counted are codon1-int1, int1-int2, codon1-codon2 and
// codon2-int2 instead of codon1-int1, int1-codon2, codon1-int2, int2-codon2, and SingleSub, given two codons that
// differ at two positions, classifies the change at the last of them. The current libsequence (1.9.8,
// src/TwoSubs.cc) makes the same call. docs/kaks.md ("Comeron") describes what follows from it.
//
// Reference: Comeron JM (1995) A method for estimating the numbers of synonymous and nonsynonymous substitutions
// per site. J Mol Evol 41:1152-1159.

// No floating-point contraction in this file (as in kaks.cpp, which says why): a * b + c is computed as written,
// never fused into one FMA instruction, so the estimates do not depend on whether the compiler and the target have
// FMA. Said here as well as in CMakeLists.txt so that it also holds where this file is compiled by something else
// (R's build of coalr). docs/kaks.md "Platform reproducibility".
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#elif defined(_MSC_VER)
#pragma fp_contract(off)
#endif

#include <cfloat>
#include <cmath>
#include <string>

#include "kaks.hpp"

namespace coal {
namespace kaks {

namespace {

inline int base_at(int c, int pos) { return pos == 0 ? c / 16 : pos == 1 ? (c / 4) % 4 : c % 4; }
inline int with_base(int c, int pos, int b) {
    static const int by[3] = {16, 4, 1};
    return c + (b - base_at(c, pos)) * by[pos];
}
// T=0 C=1 A=2 G=3. libsequence's TsTv(): "Ts" when the sum of its numbers (A 0, T 1, G 2, C 3) is even, which
// holds for A-G, T-C and for two equal bases; "Tv" otherwise.
inline bool is_ts(int a, int b) { return a == b || a + b == 1 || a + b == 5; }

// Grantham's (1974) distances as libsequence's Grantham.h has them (rows and columns in the order of AA; '*' is
// a stop codon, 999 from every amino acid).
const char AA[] = "ACDEFGHIKLMNPQRSTVWY*";
const double GRANTHAM[21][21] = {
    {0., 195.0037439, 125.7485902, 106.7125514, 112.5814896, 59.99535652, 85.49373612, 93.66243338, 105.5709011, 96.18707223, 84.43958222, 110.418026, 26.87433834, 90.76116997, 111.2645916, 99.13249708, 58.0233276, 64.43025059, 201.5127248, 111.5882137, 999},  // A
    {195.0037439, 0., 153.5743572, 169.4374593, 204.3884874, 158.2317729, 173.8426819, 197.2497963, 201.5793361, 197.4289714, 195.5543318, 138.7664683, 168.5906208, 154.2308431, 179.5348377, 111.6724957, 148.9294552, 191.232501, 127.4031665, 193.7247005, 999},  // C
    {125.7485902, 153.5743572, 0., 44.60463645, 176.5261931, 93.76755083, 81.2379796, 168.0809621, 101.4568985, 171.7572022, 160.0267338, 23.005355, 107.9096634, 61.28447445, 95.81352768, 65.47081516, 85.07613668, 152.013397, 175.3689846, 159.9374066, 999},  // D
    {106.7125514, 169.4374593, 44.60463645, 0., 140.2120632, 97.85190048, 40.79415368, 134.1621451, 56.86795522, 138.3430673, 126.0849231, 40.85895231, 93.73085383, 29.27370624, 54.01884479, 79.78756921, 65.50131254, 121.3285063, 160.4382217, 122.7768651, 999},  // E
    {112.5814896, 204.3884874, 176.5261931, 140.2120632, 0., 153.125519, 99.99505513, 21.27701376, 102.141346, 21.82392074, 28.52784356, 158.1121886, 113.7271541, 115.5869365, 97.03272116, 154.8079023, 102.8587191, 49.93518986, 142.6771181, 21.60966489, 999},  // F
    {59.99535652, 158.2317729, 93.76755083, 97.85190048, 153.125519, 0., 97.53340215, 135.4188097, 126.4579761, 137.6916434, 126.9448727, 79.34585612, 41.62910033, 87.1663871, 125.1291926, 55.26736276, 59.15645198, 108.7873884, 198.7740455, 146.9262927, 999},  // G
    {85.49373612, 173.8426819, 81.2379796, 40.79415368, 99.99505513, 97.53340215, 0., 94.33752851, 32.40288701, 98.69301471, 86.34392904, 68.35485568, 76.27712277, 24.08398867, 28.81922224, 88.93552755, 46.75286971, 83.892947, 147.2364297, 83.32844466, 999},  // H
    {93.66243338, 197.2497963, 168.0809621, 134.1621451, 21.27701376, 135.4188097, 94.33752851, 0., 101.6123327, 4.855121198, 10.12098632, 148.9140957, 95.3760884, 108.566243, 97.58656354, 141.7979605, 89.2829021, 29.60907063, 149.826083, 33.04769862, 999},  // I
    {105.5709011, 201.5793361, 101.4568985, 56.86795522, 102.141346, 126.4579761, 32.40288701, 101.6123327, 0., 106.3355126, 94.4900972, 93.88268262, 102.7139942, 53.22834015, 26.00394626, 120.5325805, 77.74105567, 96.9973994, 157.923034, 84.78644482, 999},  // K
    {96.18707223, 197.4289714, 171.7572022, 138.3430673, 21.82392074, 137.6916434, 98.69301471, 4.855121198, 106.3355126, 0., 14.30316651, 152.3309797, 97.77619789, 112.4414772, 101.8802928, 144.0792346, 92.35427468, 31.78478998, 150.0095229, 35.6773121, 999},  // L
    {84.43958222, 195.5543318, 160.0267338, 126.0849231, 28.52784356, 126.9448727, 86.34392904, 10.12098632, 94.4900972, 14.30316651, 0., 141.1526466, 86.59479895, 100.8991435, 91.63827667, 134.863884, 81.03969279, 21.52179938, 152.3965508, 35.22269546, 999},  // M
    {110.418026, 138.7664683, 23.005355, 40.85895231, 158.1121886, 79.34585612, 68.35485568, 148.9140957, 93.88268262, 152.3309797, 141.1526466, 0., 90.15739727, 45.75219661, 85.11393761, 46.23977508, 64.77428784, 132.8777227, 159.7688325, 142.2293696, 999},  // N
    {26.87433834, 168.5906208, 107.9096634, 93.73085383, 113.7271541, 41.62910033, 76.27712277, 95.3760884, 102.7139942, 97.77619789, 86.59479895, 90.15739727, 0., 75.1358782, 102.7148709, 73.35256887, 37.56362476, 67.78659846, 182.7653305, 109.6155723, 999},  // P
    {90.76116997, 154.2308431, 61.28447445, 29.27370624, 115.5869365, 87.1663871, 24.08398867, 108.566243, 53.22834015, 112.4414772, 100.8991435, 45.75219661, 75.1358782, 0., 42.81393193, 68.19787424, 41.10476252, 96.32571878, 141.5619982, 98.78325151, 999},  // Q
    {111.2645916, 179.5348377, 95.81352768, 54.01884479, 97.03272116, 125.1291926, 28.81922224, 97.58656354, 26.00394626, 101.8802928, 91.63827667, 85.11393761, 102.7148709, 42.81393193, 0., 109.213153, 70.97102635, 95.79707649, 132.5874819, 77.10764416, 999},  // R
    {99.13249708, 111.6724957, 65.47081516, 79.78756921, 154.8079023, 55.26736276, 88.93552755, 141.7979605, 120.5325805, 144.0792346, 134.863884, 46.23977508, 73.35256887, 68.19787424, 109.213153, 0., 57.74902766, 123.0339581, 157.855071, 143.1073844, 999},  // S
    {58.0233276, 148.9294552, 85.07613668, 65.50131254, 102.8587191, 59.15645198, 46.75286971, 89.2829021, 77.74105567, 92.35427468, 81.03969279, 64.77428784, 37.56362476, 41.10476252, 70.97102635, 57.74902766, 0., 69.49621703, 150.752701, 92.24764162, 999},  // T
    {64.43025059, 191.232501, 152.013397, 121.3285063, 49.93518986, 108.7873884, 83.892947, 29.60907063, 96.9973994, 31.78478998, 21.52179938, 132.8777227, 67.78659846, 96.32571878, 95.79707649, 123.0339581, 69.49621703, 0., 162.8552264, 54.66278181, 999},  // V
    {201.5127248, 127.4031665, 175.3689846, 160.4382217, 142.6771181, 198.7740455, 147.2364297, 149.826083, 157.923034, 150.0095229, 152.3965508, 159.7688325, 182.7653305, 141.5619982, 132.5874819, 157.855071, 150.752701, 162.8552264, 0., 128.9732994, 999},  // W
    {111.5882137, 193.7247005, 159.9374066, 122.7768651, 21.60966489, 146.9262927, 83.32844466, 33.04769862, 84.78644482, 35.6773121, 35.22269546, 142.2293696, 109.6155723, 98.78325151, 77.10764416, 143.1073844, 92.24764162, 54.66278181, 128.9732994, 0., 999},  // Y
    {999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 0.},  // *
};

// RedundancyCom95: for every codon, the degeneracy of its first and third position, and the L values.
struct Tables {
    double firstNon[64], first2S[64], first2V[64];
    double thirdNon[64], third2S[64], third2V[64], thirdFour[64];
    double l0[64], l2S[64], l2V[64], l4[64];
    int aa[64];  // index into AA
};

Tables build_tables(const GeneticCode& g) {
    Tables t{};
    for (int c = 0; c < 64; ++c) {
        int k = 0;
        while (AA[k] && AA[k] != g.aa[c]) ++k;
        t.aa[c] = k;  // the standard code has the 20 amino acids and '*' only
    }
    for (int pos = 0; pos <= 2; pos += 2) {
        for (int c = 0; c < 64; ++c) {
            int numTsS = 0, numTvS = 0, numPossTs = 0, numPossTv = 0;
            const int state = base_at(c, pos);
            for (int l = 0; l <= 3; ++l) {
                if (l == state) continue;
                const int m = with_base(c, pos, l);
                const bool ts = is_ts(state, l);
                if (g.aa[c] != '*' && g.aa[m] != '*') {
                    const bool S = g.aa[c] == g.aa[m];
                    if (ts)
                        ++numPossTs;
                    else
                        ++numPossTv;
                    if (S && ts)
                        ++numTsS;
                    else if (S)
                        ++numTvS;
                }
            }
            if (pos == 0) {
                double FirstN, First2S, First2V;
                if ((numPossTs + numPossTv) == 0) {  // stop codons
                    First2S = 0;
                    First2V = 0;
                    FirstN = 0;
                } else if ((numTsS + numTvS) == 0) {  // non-degenerate sites
                    First2S = 0;
                    First2V = 0;
                    FirstN = 1.0;
                } else if (numPossTs != 0 && numPossTv != 0) {  // special cases
                    if (double(numTsS) / double(numPossTs) != 1.0 && double(numTvS) / double(numPossTv) != 1.0) {
                        // fractional redundancy
                        First2S = double(numTsS) / double(numPossTs + numPossTv);
                        First2V = double(numTvS) / double(numPossTs + numPossTv);
                        FirstN = 1.0 - First2S - First2V;
                    } else {  // odd types of degeneracy
                        First2S = double(numTsS) / double(numPossTs);
                        First2V = double(numTvS) / double(numPossTv);
                        FirstN = 1.0 - First2S - First2V;
                    }
                } else {
                    if (numPossTs > 0)
                        First2S = double(numTsS) / double(numPossTs + numPossTv);
                    else
                        First2S = 0.0;
                    if (numPossTv > 0)
                        First2V = double(numTvS) / double(numPossTv + numPossTs);
                    else
                        First2V = 0.0;
                    FirstN = 1.0 - First2S - First2V;
                }
                t.firstNon[c] = FirstN;
                t.first2S[c] = First2S;
                t.first2V[c] = First2V;
            } else {
                double ThirdN = 0.0, Third2S = 0.0, Third2V = 0.0, Third4 = 0.0;
                if (numPossTs + numPossTv == 0) {  // stop codons
                } else if (numTsS + numTvS == 0) {  // nondegenerate
                    ThirdN = 1.0;
                } else if (numTsS + numTvS == 3) {  // 4fold degenerate
                    Third4 = 1.0;
                } else if (numTsS == 0 || numTvS == 0) {  // non-degeneracy for one type
                    if (numPossTs > 0)
                        Third2S = double(numTsS) / double(numPossTs);
                    else
                        Third2S = 0.0;
                    if (numPossTv > 0)
                        Third2V = double(numTvS) / double(numPossTv);
                    else
                        Third2V = 0.0;
                    ThirdN = 1.0 - Third2S - Third2V;
                    Third4 = 0.0;
                } else if (numTsS > 0 && numTvS > 0 &&
                           (double(numTsS) / double(numPossTs) != 1.0 ||
                            double(numTvS) / double(numPossTv) != 1.0)) {  // fractional degeneracy
                    Third2S = 1.0 / 3.0;
                    Third2V = 1.0 / 3.0;
                    ThirdN = 1.0 / 3.0;
                    Third4 = 0.0;
                }  // (no other case sets anything: all four stay 0, as upstream)
                t.thirdNon[c] = ThirdN;
                t.third2S[c] = Third2S;
                t.third2V[c] = Third2V;
                t.thirdFour[c] = Third4;
            }
        }
    }
    for (int c = 0; c < 64; ++c) {
        if (g.aa[c] != '*') {
            t.l0[c] = 1. + t.firstNon[c] + t.thirdNon[c];
            t.l2S[c] = t.first2S[c] + t.third2S[c];
            t.l2V[c] = t.first2V[c] + t.third2V[c];
            t.l4[c] = t.thirdFour[c];
        } else {
            t.l0[c] = t.l2S[c] = t.l2V[c] = t.l4[c] = 0.0;
        }
    }
    return t;
}

const Tables& tables() {
    // built once (C++11 guarantees a thread-safe initialisation) and never written again
    static const Tables t = build_tables(*get_code(1));
    return t;
}

// the eight counts: transitions (P) and transversions (Q) at non-degenerate (0), transitional-degenerate (2S),
// transversional-degenerate (2V) and fourfold-degenerate (4) sites
enum { P0, P2S, P2V, P4, Q0, Q2S, Q2V, Q4, NCOUNT };
struct Counts {
    double v[NCOUNT] = {0, 0, 0, 0, 0, 0, 0, 0};
};

// SingleSub::Calculate for one of the two codons: the class the change at `pos` is counted in
inline int single_class(const Tables& t, int codon, int pos, bool ts) {
    switch (pos) {
        case 0:
            if (ts) {
                if (t.firstNon[codon] < 1.) return t.first2S[codon] > 0. ? P2S : P2V;
                return P0;
            }
            if (t.firstNon[codon] < 1.) return t.first2V[codon] > 0. ? Q2V : Q2S;
            return Q0;
        case 1: return ts ? P0 : Q0;
        default:
            if (ts) {
                if (t.thirdNon[codon] < 1.) {
                    if (t.thirdFour[codon] < 1.) return t.third2S[codon] > 0. ? P2S : P2V;
                    return P4;
                }
                return P0;
            }
            if (t.thirdNon[codon] < 1.) {
                if (t.thirdFour[codon] < 1.) return t.third2V[codon] > 0. ? Q2V : Q2S;
                return Q4;
            }
            return Q0;
    }
}

// SingleSub: the change between two codons, classified in each codon and averaged. As upstream, the position
// is the LAST one at which the codons differ (position 0 when they do not differ), whatever the number of
// differences.
Counts single_sub(const Tables& t, int c1, int c2) {
    int pos = 0;
    for (int k = 0; k <= 2; ++k)
        if (base_at(c1, k) != base_at(c2, k)) pos = k;
    const bool ts = is_ts(base_at(c1, pos), base_at(c2, pos));
    double vi[NCOUNT] = {0, 0, 0, 0, 0, 0, 0, 0}, vj[NCOUNT] = {0, 0, 0, 0, 0, 0, 0, 0};
    vi[single_class(t, c1, pos, ts)] += 1.;
    vj[single_class(t, c2, pos, ts)] += 1.;
    Counts out;
    for (int k = 0; k < NCOUNT; ++k) out.v[k] = (vi[k] + vj[k]) / 2.0;
    return out;
}

inline double gdist(const Tables& t, int c1, int c2) { return GRANTHAM[t.aa[c1]][t.aa[c2]]; }

// TwoSubs with GranthamWeights2
Counts two_subs(const Tables& t, int c1, int c2) {
    int pos[2] = {0, 0};
    for (int i = 0, j = 0; i <= 2; ++i)
        if (base_at(c1, i) != base_at(c2, i)) pos[j++] = i;
    const int int0 = with_base(c1, pos[0], base_at(c2, pos[0]));
    const int int1 = with_base(c1, pos[1], base_at(c2, pos[1]));
    double len_path_1 = 0.0, len_path_2 = 0.0;
    len_path_1 += gdist(t, c1, int0);
    len_path_1 += gdist(t, int0, c2);
    len_path_2 += gdist(t, c1, int1);
    len_path_2 += gdist(t, int1, c2);
    double w[2] = {0., 0.}, w_tot = 0.;
    if (std::fabs(len_path_1 - 0.) <= DBL_EPSILON && std::fabs(len_path_2 - 0.) <= DBL_EPSILON) {
        w[0] = w[1] = 0.5;
        w_tot = 1.;
    } else {
        w[0] = 1. - len_path_1 / (len_path_1 + len_path_2);
        w[1] = 1. - len_path_2 / (len_path_1 + len_path_2);
        w_tot = w[0] + w[1];
    }
    w[0] /= w_tot;
    w[1] /= w_tot;
    // upstream: Calculate(sitesObj, codon1, intermediates[0], codon2, intermediates[1], ...) into the parameters
    // (codon1, int_1, int_2, codon2): see the header of this file
    const Counts b1 = single_sub(t, c1, int0);
    const Counts b2 = single_sub(t, int0, int1);
    const Counts b3 = single_sub(t, c1, c2);
    const Counts b4 = single_sub(t, c2, int1);
    Counts out;
    for (int k = 0; k < NCOUNT; ++k) out.v[k] = (b1.v[k] + b2.v[k]) * w[0] + (b3.v[k] + b4.v[k]) * w[1];
    return out;
}

// ThreeSubs with GranthamWeights3
Counts three_subs(const Tables& t, int c1, int c2) {
    auto mix = [&](bool p0, bool p1, bool p2) {  // codon1 with the listed positions taken from codon2
        int c = c1;
        if (p0) c = with_base(c, 0, base_at(c2, 0));
        if (p1) c = with_base(c, 1, base_at(c2, 1));
        if (p2) c = with_base(c, 2, base_at(c2, 2));
        return c;
    };
    const int im[9] = {mix(1, 0, 0), mix(1, 1, 0), mix(1, 0, 1), mix(0, 1, 0), mix(1, 1, 0),
                       mix(0, 1, 1), mix(0, 0, 1), mix(1, 0, 1), mix(0, 1, 1)};
    // the six paths: codon1 - first - second - codon2
    static const int first[6] = {0, 0, 3, 3, 6, 6}, second[6] = {1, 2, 4, 5, 7, 8};
    double w[6], w_tot;
    for (int k = 0; k < 6; ++k) {
        double len = 0.;
        len += gdist(t, c1, im[first[k]]);
        len += gdist(t, im[first[k]], im[second[k]]);
        len += gdist(t, im[second[k]], c2);
        w[k] = 1. / len;
    }
    w_tot = w[0] + w[1] + w[2] + w[3] + w[4] + w[5];
    for (int k = 0; k < 6; ++k) w[k] /= w_tot;
    // the 15 branches, in upstream's numbering
    const Counts b[15] = {single_sub(t, c1, im[0]),    single_sub(t, im[0], im[1]), single_sub(t, im[1], c2),
                          single_sub(t, im[0], im[2]), single_sub(t, im[2], c2),    single_sub(t, c1, im[3]),
                          single_sub(t, im[3], im[4]), single_sub(t, im[4], c2),    single_sub(t, im[3], im[5]),
                          single_sub(t, im[5], c2),    single_sub(t, c1, im[6]),    single_sub(t, im[6], im[7]),
                          single_sub(t, im[7], c2),    single_sub(t, im[6], im[8]), single_sub(t, im[8], c2)};
    Counts out;
    for (int k = 0; k < NCOUNT; ++k)
        out.v[k] = (b[0].v[k] + b[1].v[k] + b[2].v[k]) * w[0] + (b[0].v[k] + b[3].v[k] + b[4].v[k]) * w[1] +
                   (b[5].v[k] + b[6].v[k] + b[7].v[k]) * w[2] + (b[5].v[k] + b[8].v[k] + b[9].v[k]) * w[3] +
                   (b[10].v[k] + b[11].v[k] + b[12].v[k]) * w[4] + (b[10].v[k] + b[13].v[k] + b[14].v[k]) * w[5];
    return out;
}

inline int num_diffs(int c1, int c2) {
    return (base_at(c1, 0) != base_at(c2, 0)) + (base_at(c1, 1) != base_at(c2, 1)) +
           (base_at(c1, 2) != base_at(c2, 2));
}

void add_note(Result& r, const std::string& msg) { r.note = r.note.empty() ? msg : r.note + "; " + msg; }

}  // namespace

namespace detail {

Result run_comeron(const Prepared& p, const GeneticCode& code, const Options& opt, Result r) {
    (void)code;  // the standard code (compute() refuses the others for this method)
    const bool legacy = opt.comeron_flavour == ComeronFlavour::gestimator;
    const std::vector<uint8_t>& s1 = legacy ? p.r1 : p.c1;
    const std::vector<uint8_t>& s2 = legacy ? p.r2 : p.c2;
    const Tables& t = tables();

    // Sites::count_sites (maxhits = 3: every codon column counts)
    double L0 = 0., L2S = 0., L2V = 0.0, L4 = 0.0;
    for (size_t i = 0; i < s1.size(); ++i) {
        const int a = s1[i], b = s2[i];
        L0 += (t.l0[a] + t.l0[b]) / 2.0;
        L2S += (t.l2S[a] + t.l2S[b]) / 2.0;
        L2V += (t.l2V[a] + t.l2V[b]) / 2.0;
        L4 += (t.l4[a] + t.l4[b]) / 2.0;
    }

    // Comeron95::diverge
    double c[NCOUNT] = {0, 0, 0, 0, 0, 0, 0, 0};
    double subs = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        const int a = s1[i], b = s2[i];
        if (a == b) continue;
        const int ndiff = num_diffs(a, b);
        subs += ndiff;
        const Counts k = ndiff == 1 ? single_sub(t, a, b) : ndiff == 2 ? two_subs(t, a, b) : three_subs(t, a, b);
        for (int j = 0; j < NCOUNT; ++j) c[j] += k.v[j];
    }
    for (int j = 0; j < NCOUNT; ++j)
        if (!std::isfinite(c[j])) c[j] = 0.0;
    const double p0 = c[P0], p2S = c[P2S], p2V = c[P2V], p4 = c[P4];
    const double q0 = c[Q0], q2S = c[Q2S], q2V = c[Q2V], q4 = c[Q4];
    r.substitutions = subs;

    // Comeron95::omega
    bool saturated_s = false, saturated_a = false;  // a logarithm with an argument <= 0
    // upstream: log(x), and 0 when the result is not finite (NaN argument from 0/0: no site of the class; an
    // argument <= 0: more differences than the correction allows)
    auto checked_log = [](double x, bool& flag) {
        if (x != x) return 0.0;
        if (x <= 0.0) {
            flag = true;
            return 0.0;
        }
        const double l = std::log(x);
        if (!std::isfinite(l)) return 0.0;
        return l;
    };
    double log1, log2;
    double Qs = (q2V + q4) / (L2V + L4);
    if (!std::isfinite(Qs)) Qs = 0.0;
    double Bs = (-0.5) * std::log(1.0 - (2.0 * Qs));
    if (1.0 - (2.0 * Qs) <= 0.0) saturated_s = true;  // Bs is then inf or NaN: upstream prints 999 for dS
    double Qa = (q0 + q2S) / (L0 + L2S);
    if (!std::isfinite(Qa)) Qa = 0.0;
    double Ba = (-0.5) * std::log(1.0 - (2.0 * Qa));
    bool ba_kimura_rule = false;
    if (!std::isfinite(Ba)) {
        saturated_a = true;  // 1 - 2 Qa <= 0
        if (legacy) {
            // upstream's ad-hoc rule: "if sites aren't saturated in general (i.e. Kimura's distance < 1.0) [...]
            // Ba should equal 0.0". Kimura80 on the whole rows:
            const double P = double(p.k80_ts) / double(p.k80_sites);
            const double Q = double(p.k80_tv) / double(p.k80_sites);
            double divergence;
            if (std::fabs(1.0 - 2.0 * P - Q) > DBL_EPSILON)
                divergence = -1.0 * 0.5 * std::log((1.0 - 2.0 * P - Q) * std::pow((1 - 2.0 * Q), 0.5));
            else
                divergence = 0.;
            if (divergence <= 0.0 - DBL_EPSILON) divergence = 0.0;
            const double K = std::isfinite(divergence) ? divergence : 999.0;
            if (K < 1.0) {
                Ba = 0.0;
                ba_kimura_rule = true;
            }
        }
    }
    const double P2S_site = p2S / L2S;
    const double P2V_site = p2V / L2V;
    const double P0_site = p0 / L0;
    const double Q0_site = q0 / L0;
    const double P4_site = p4 / L4;
    const double Q4_site = q4 / L4;

    // (a class without sites has weight 0 in As / Aa: an undefined logarithm there is no saturation)
    bool none = false;
    log1 = checked_log(1.0 - (2.0 * P2S_site) - Qa, L2S > 0 ? saturated_s : none);
    log2 = checked_log(1.0 - (2.0 * Qa), L2S > 0 ? saturated_s : none);
    const double A2S = (-0.5) * log1 + (0.25) * log2;
    log1 = checked_log(1.0 - (2.0 * P4_site) - Q4_site, L4 > 0 ? saturated_s : none);
    log2 = checked_log(1.0 - (2.0 * Q4_site), L4 > 0 ? saturated_s : none);
    const double A4 = (-0.5) * log1 + (0.25) * log2;
    double As = (L2S * A2S + L4 * A4) / (L2S + L4);

    log1 = checked_log(1.0 - (2.0 * P2V_site) - Qs, L2V > 0 ? saturated_a : none);
    log2 = checked_log(1.0 - (2.0 * Qs), L2V > 0 ? saturated_a : none);
    const double A2V = (-0.5) * log1 + (0.25) * log2;
    log1 = checked_log(1.0 - (2.0 * P0_site) - Q0_site, L0 > 0 ? saturated_a : none);
    log2 = checked_log(1.0 - (2.0 * Q0_site), L0 > 0 ? saturated_a : none);
    const double A0 = (-0.5) * log1 + (0.25) * log2;
    double Aa = (L2V * A2V + L0 * A0) / (L2V + L0);

    if (As <= 0.0) As = 0.0;
    if (Bs <= 0.0) Bs = 0.0;
    if (Aa <= 0.0) Aa = 0.0;
    if (Ba <= 0.0) Ba = 0.0;
    const double dS = As + Bs;
    const double dN = Aa + Ba;

    // what upstream prints as 999 ("cannot be calculated"): not finite
    const bool s_undefined = !std::isfinite(dS), a_undefined = !std::isfinite(dN);
    const bool no_s_sites = !(L2S + L4 > 0), no_a_sites = !(L2V + L0 > 0);
    auto saturated = [&](const char* which) {
        if (r.status != Status::no_data) r.status = Status::saturated;
        add_note(r, std::string(which) + " saturated (multiple-hit correction undefined)");
    };
    auto no_sites = [&](const char* msg) {
        r.status = Status::no_data;
        add_note(r, msg);
    };
    if (legacy) {
        // gestimator's numbers. NA only where it prints 999; where it replaced an undefined logarithm by 0 the
        // number is given, with a warning.
        r.ka = a_undefined ? NA : dN;
        r.ks = s_undefined ? NA : dS;
        if (a_undefined) {
            if (no_a_sites)
                no_sites("no nonsynonymous sites");
            else
                saturated("Ka");
        }
        if (s_undefined) {
            if (no_s_sites)
                no_sites("no synonymous sites");
            else
                saturated("Ks");
        }
        if (!a_undefined && saturated_a) {
            if (r.status == Status::ok) r.status = Status::warning;
            add_note(r, ba_kimura_rule
                            ? "Ka: correction undefined, gestimator's value (Ba set to 0 by its Kimura rule, "
                              "undefined logarithms set to 0)"
                            : "Ka: correction undefined, gestimator's value (undefined logarithms set to 0)");
        }
        if (!s_undefined && saturated_s) {
            if (r.status == Status::ok) r.status = Status::warning;
            add_note(r, "Ks: correction undefined, gestimator's value (undefined logarithms set to 0)");
        }
        // gestimator's ratio: 999 when dN or dS is 999 or dS == 0
        if (!a_undefined && !s_undefined && dS != 0.) r.ka_ks = dN / dS;
    } else {
        if (no_a_sites) {
            no_sites("no nonsynonymous sites");
        } else if (saturated_a || a_undefined) {
            saturated("Ka");
        } else {
            r.ka = dN;
        }
        if (no_s_sites) {
            no_sites("no synonymous sites");
        } else if (saturated_s || s_undefined) {
            saturated("Ks");
        } else {
            r.ks = dS;
        }
    }
    return r;
}

}  // namespace detail
}  // namespace kaks
}  // namespace coal
