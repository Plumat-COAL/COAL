// Vendored from COAL cpp/src/align.cpp by scripts/sync_r_sources.py. Do not edit here: edit cpp/src and sync.
// Licence: GPL-3.0-only, as COAL (https://github.com/Plumat-COAL/COAL); copyright and third-party notices
// (KaKs_Calculator 2.0 and PAML for kaks.cpp / kaks.hpp, libsequence for kaks_comeron.cpp, seqinr for
// kaks_li.cpp): inst/COPYRIGHTS of this package.
// Pairwise protein alignment (align.hpp): Gotoh global alignment with a row-wise, vectorisable kernel, and the
// expansion of local (DIAMOND) hits to full length.
//
// Kernel. With H = best score ending at (i, j), X = ending in a gap in b (a[i] against '-', vertical),
// Y = ending in a gap in a (horizontal), M = H[i-1][j-1] + s(a[i], b[j]) and D = max(M, X):
//     X[i][j] = max(H[i-1][j] - open, X[i-1][j] - extend)       depends on the previous row only
//     Y[i][j] = max(D[i][j-1] - open, Y[i][j-1] - extend)       (H[i][j-1] - open would add nothing:
//                                                              Y - open <= Y - extend when open >= extend)
//     H[i][j] = max(D[i][j], Y[i][j])
// Y unrolls to a max-plus prefix scan: Y[i][j] = q[j-1] - open with q[k] = max_{l <= k} (D[l] - (k - l) * extend)
// (D[i][0] := H[i][0]). So a row is three loops: M, X, D (vectorised by the compiler), the scan (the only
// loop-carried dependency; eight lanes at a time with GCC / Clang vector extensions), and Y, H and the traceback
// bits (vectorised). On x86 with GCC / Clang the row loops also exist as an AVX2 build, used when the CPU has AVX2. The query profile (one row of scores per residue
// type, along b) keeps the inner loop free of table lookups. Traceback: 4 bits per cell (H from M / X / Y, X and Y
// open or extend), packed two cells per byte. Scores are integers: the parameters are scaled by the smallest
// factor that makes the gap costs whole numbers (2 for EMBOSS's extend 0.5).
#include "align.hpp"

#include <algorithm>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstring>
#include <memory>

#include "core.hpp"
#include "pool.hpp"

// x86 with GCC / Clang: the row loops are compiled twice, for the baseline (SSE2) and for AVX2, and the AVX2 copy is
// used when the CPU has it (checked once). Elsewhere (MSVC, ARM) only the baseline, which the compiler vectorises
// for its target. COAL_ALIGN_SCALAR: baseline only, everywhere (for benchmarks and tests of that path).
#if (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__)) && \
    !defined(COAL_ALIGN_SCALAR)
#define COAL_ALIGN_AVX2 1
#define COAL_TARGET_AVX2 __attribute__((target("avx2")))
#else
#define COAL_ALIGN_AVX2 0
#endif
#if defined(__GNUC__) || defined(__clang__)
#define COAL_INLINE inline __attribute__((always_inline))
#else
#define COAL_INLINE inline
#endif

namespace coal {
namespace {

const char ALPHABET[] = "ARNDCQEGHILKMFPSTWYVBZX*";
constexpr int X_CODE = 22;

// Source: the NCBI matrices (ftp.ncbi.nih.gov/blast/matrices), as distributed with Biopython
// (Bio/Align/substitution_matrices/data); the same values as EMBOSS's EBLOSUM62.
const SubstitutionMatrix MATRICES[] = {
    {"BLOSUM62", {
        {  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
        { -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
        { -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
        { -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
        {  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
        { -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
        { -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
        {  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
        { -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
        { -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
        { -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
        { -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
        { -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
        { -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
        { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
        {  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
        {  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
        { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
        { -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
        {  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
        { -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
        { -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
        {  0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
        { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1},
    }},
    {"BLOSUM45", {
        {  5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -2, -2,  0, -1, -1,  0, -5},
        { -2,  7,  0, -1, -3,  1,  0, -2,  0, -3, -2,  3, -1, -2, -2, -1, -1, -2, -1, -2, -1,  0, -1, -5},
        { -1,  0,  6,  2, -2,  0,  0,  0,  1, -2, -3,  0, -2, -2, -2,  1,  0, -4, -2, -3,  4,  0, -1, -5},
        { -2, -1,  2,  7, -3,  0,  2, -1,  0, -4, -3,  0, -3, -4, -1,  0, -1, -4, -2, -3,  5,  1, -1, -5},
        { -1, -3, -2, -3, 12, -3, -3, -3, -3, -3, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -2, -3, -2, -5},
        { -1,  1,  0,  0, -3,  6,  2, -2,  1, -2, -2,  1,  0, -4, -1,  0, -1, -2, -1, -3,  0,  4, -1, -5},
        { -1,  0,  0,  2, -3,  2,  6, -2,  0, -3, -2,  1, -2, -3,  0,  0, -1, -3, -2, -3,  1,  4, -1, -5},
        {  0, -2,  0, -1, -3, -2, -2,  7, -2, -4, -3, -2, -2, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -5},
        { -2,  0,  1,  0, -3,  1,  0, -2, 10, -3, -2, -1,  0, -2, -2, -1, -2, -3,  2, -3,  0,  0, -1, -5},
        { -1, -3, -2, -4, -3, -2, -3, -4, -3,  5,  2, -3,  2,  0, -2, -2, -1, -2,  0,  3, -3, -3, -1, -5},
        { -1, -2, -3, -3, -2, -2, -2, -3, -2,  2,  5, -3,  2,  1, -3, -3, -1, -2,  0,  1, -3, -2, -1, -5},
        { -1,  3,  0,  0, -3,  1,  1, -2, -1, -3, -3,  5, -1, -3, -1, -1, -1, -2, -1, -2,  0,  1, -1, -5},
        { -1, -1, -2, -3, -2,  0, -2, -2,  0,  2,  2, -1,  6,  0, -2, -2, -1, -2,  0,  1, -2, -1, -1, -5},
        { -2, -2, -2, -4, -2, -4, -3, -3, -2,  0,  1, -3,  0,  8, -3, -2, -1,  1,  3,  0, -3, -3, -1, -5},
        { -1, -2, -2, -1, -4, -1,  0, -2, -2, -2, -3, -1, -2, -3,  9, -1, -1, -3, -3, -3, -2, -1, -1, -5},
        {  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -3, -1, -2, -2, -1,  4,  2, -4, -2, -1,  0,  0,  0, -5},
        {  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1,  2,  5, -3, -1,  0,  0, -1,  0, -5},
        { -2, -2, -4, -4, -5, -2, -3, -2, -3, -2, -2, -2, -2,  1, -3, -4, -3, 15,  3, -3, -4, -2, -2, -5},
        { -2, -1, -2, -2, -3, -1, -2, -3,  2,  0,  0, -1,  0,  3, -3, -2, -1,  3,  8, -1, -2, -2, -1, -5},
        {  0, -2, -3, -3, -1, -3, -3, -3, -3,  3,  1, -2,  1,  0, -3, -1,  0, -3, -1,  5, -3, -3, -1, -5},
        { -1, -1,  4,  5, -2,  0,  1, -1,  0, -3, -3,  0, -2, -3, -2,  0,  0, -4, -2, -3,  4,  2, -1, -5},
        { -1,  0,  0,  1, -3,  4,  4, -2,  0, -3, -2,  1, -1, -3, -1,  0, -1, -2, -2, -3,  2,  4, -1, -5},
        {  0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -2, -1, -1, -1, -1, -1, -5},
        { -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1},
    }},
    {"BLOSUM80", {
        {  7, -3, -3, -3, -1, -2, -2,  0, -3, -3, -3, -1, -2, -4, -1,  2,  0, -5, -4, -1, -3, -2, -1, -8},
        { -3,  9, -1, -3, -6,  1, -1, -4,  0, -5, -4,  3, -3, -5, -3, -2, -2, -5, -4, -4, -2,  0, -2, -8},
        { -3, -1,  9,  2, -5,  0, -1, -1,  1, -6, -6,  0, -4, -6, -4,  1,  0, -7, -4, -5,  5, -1, -2, -8},
        { -3, -3,  2, 10, -7, -1,  2, -3, -2, -7, -7, -2, -6, -6, -3, -1, -2, -8, -6, -6,  6,  1, -3, -8},
        { -1, -6, -5, -7, 13, -5, -7, -6, -7, -2, -3, -6, -3, -4, -6, -2, -2, -5, -5, -2, -6, -7, -4, -8},
        { -2,  1,  0, -1, -5,  9,  3, -4,  1, -5, -4,  2, -1, -5, -3, -1, -1, -4, -3, -4, -1,  5, -2, -8},
        { -2, -1, -1,  2, -7,  3,  8, -4,  0, -6, -6,  1, -4, -6, -2, -1, -2, -6, -5, -4,  1,  6, -2, -8},
        {  0, -4, -1, -3, -6, -4, -4,  9, -4, -7, -7, -3, -5, -6, -5, -1, -3, -6, -6, -6, -2, -4, -3, -8},
        { -3,  0,  1, -2, -7,  1,  0, -4, 12, -6, -5, -1, -4, -2, -4, -2, -3, -4,  3, -5, -1,  0, -2, -8},
        { -3, -5, -6, -7, -2, -5, -6, -7, -6,  7,  2, -5,  2, -1, -5, -4, -2, -5, -3,  4, -6, -6, -2, -8},
        { -3, -4, -6, -7, -3, -4, -6, -7, -5,  2,  6, -4,  3,  0, -5, -4, -3, -4, -2,  1, -7, -5, -2, -8},
        { -1,  3,  0, -2, -6,  2,  1, -3, -1, -5, -4,  8, -3, -5, -2, -1, -1, -6, -4, -4, -1,  1, -2, -8},
        { -2, -3, -4, -6, -3, -1, -4, -5, -4,  2,  3, -3,  9,  0, -4, -3, -1, -3, -3,  1, -5, -3, -2, -8},
        { -4, -5, -6, -6, -4, -5, -6, -6, -2, -1,  0, -5,  0, 10, -6, -4, -4,  0,  4, -2, -6, -6, -3, -8},
        { -1, -3, -4, -3, -6, -3, -2, -5, -4, -5, -5, -2, -4, -6, 12, -2, -3, -7, -6, -4, -4, -2, -3, -8},
        {  2, -2,  1, -1, -2, -1, -1, -1, -2, -4, -4, -1, -3, -4, -2,  7,  2, -6, -3, -3,  0, -1, -1, -8},
        {  0, -2,  0, -2, -2, -1, -2, -3, -3, -2, -3, -1, -1, -4, -3,  2,  8, -5, -3,  0, -1, -2, -1, -8},
        { -5, -5, -7, -8, -5, -4, -6, -6, -4, -5, -4, -6, -3,  0, -7, -6, -5, 16,  3, -5, -8, -5, -5, -8},
        { -4, -4, -4, -6, -5, -3, -5, -6,  3, -3, -2, -4, -3,  4, -6, -3, -3,  3, 11, -3, -5, -4, -3, -8},
        { -1, -4, -5, -6, -2, -4, -4, -6, -5,  4,  1, -4,  1, -2, -4, -3,  0, -5, -3,  7, -6, -4, -2, -8},
        { -3, -2,  5,  6, -6, -1,  1, -2, -1, -6, -7, -1, -5, -6, -4,  0, -1, -8, -5, -6,  6,  0, -3, -8},
        { -2,  0, -1,  1, -7,  5,  6, -4,  0, -6, -5,  1, -3, -6, -2, -1, -2, -5, -4, -4,  0,  6, -1, -8},
        { -1, -2, -2, -3, -4, -2, -2, -3, -2, -2, -2, -2, -2, -3, -3, -1, -1, -5, -3, -2, -3, -1, -2, -8},
        { -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1},
    }},
    {"PAM250", {
        {  2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0,  0,  0,  0, -8},
        { -2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0, -4,  0,  0, -1,  2, -4, -2, -1,  0, -1, -8},
        {  0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2, -3,  0,  1,  0, -4, -2, -2,  2,  1,  0, -8},
        {  0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3, -6, -1,  0,  0, -7, -4, -2,  3,  3, -1, -8},
        { -2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3,  0, -2, -8,  0, -2, -4, -5, -3, -8},
        {  0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1, -5,  0, -1, -1, -5, -4, -2,  1,  3, -1, -8},
        {  0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2, -5, -1,  0,  0, -7, -4, -2,  3,  3, -1, -8},
        {  1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3, -5,  0,  1,  0, -7, -5, -1,  0,  0, -1, -8},
        { -1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2, -2,  0, -1, -1, -3,  0, -2,  1,  2, -1, -8},
        { -1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,  1, -2, -1,  0, -5, -1,  4, -2, -2, -1, -8},
        { -2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,  2, -3, -3, -2, -2, -1,  2, -3, -3, -1, -8},
        { -1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0, -5, -1,  0,  0, -3, -4, -2,  1,  0, -1, -8},
        { -1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,  0, -2, -2, -1, -4, -2,  2, -2, -2, -1, -8},
        { -3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9, -5, -3, -3,  0,  7, -1, -4, -5, -2, -8},
        {  1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6,  1,  0, -6, -5, -1, -1,  0, -1, -8},
        {  1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2,  1, -2, -3, -1,  0,  0,  0, -8},
        {  1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3, -5, -3,  0,  0, -1,  0, -8},
        { -6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17,  0, -6, -5, -6, -4, -8},
        { -3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10, -2, -3, -4, -2, -8},
        {  0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4, -2, -2, -1, -8},
        {  0, -1,  2,  3, -4,  1,  3,  0,  1, -2, -3,  1, -2, -4, -1,  0,  0, -5, -3, -2,  3,  2, -1, -8},
        {  0,  0,  1,  3, -5,  3,  3,  0,  2, -2, -3,  0, -2, -5,  0,  0, -1, -6, -4, -2,  2,  3, -1, -8},
        {  0, -1,  0, -1, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1,  0,  0, -4, -2, -1, -1, -1, -1, -8},
        { -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1},
    }},
};

// residue byte -> 0..23, 255 = not a residue
struct Codes {
    uint8_t c[256];
    Codes() {
        std::memset(c, 255, sizeof c);
        for (int ch = 'A'; ch <= 'Z'; ++ch) {
            const char* p = std::strchr(ALPHABET, ch);
            c[ch] = c[std::tolower(ch)] = (uint8_t)(p ? p - ALPHABET : X_CODE);
        }
        c['*'] = 23;
    }
};
const Codes& codes() {
    static const Codes k;
    return k;
}

std::string upper_name(std::string s) {
    for (char& ch : s) ch = (char)std::toupper((unsigned char)ch);
    return s;
}

struct Scaled {
    const SubstitutionMatrix* m = nullptr;
    int scale = 1;
    int open = 0, ext = 0;
    int sub[24][24];
    int max_abs = 0;
};

Scaled scale_params(const AlignParams& p) {
    check_align_params(p);
    Scaled s;
    s.m = substitution_matrix(p.matrix);
    static const int SCALES[] = {1, 2, 4, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000};
    for (int k : SCALES) {
        double o = p.gap_open * k, e = p.gap_extend * k;
        if (std::fabs(o - std::round(o)) < 1e-6 && std::fabs(e - std::round(e)) < 1e-6) {
            s.scale = k;
            s.open = (int)std::lround(o);
            s.ext = (int)std::lround(e);
            break;
        }
    }
    s.max_abs = s.open;
    for (int i = 0; i < 24; ++i)
        for (int j = 0; j < 24; ++j) {
            s.sub[i][j] = s.m->score[i][j] * s.scale;
            s.max_abs = std::max(s.max_abs, std::abs(s.sub[i][j]));
        }
    return s;
}

std::vector<uint8_t> encode(const std::string& seq, const char* which) {
    const Codes& k = codes();
    std::vector<uint8_t> out(seq.size());
    for (size_t i = 0; i < seq.size(); ++i) {
        uint8_t c = k.c[(unsigned char)seq[i]];
        if (c == 255) {
            unsigned char ch = (unsigned char)seq[i];
            std::string shown = std::isprint(ch) ? std::string("'") + (char)ch + "'" : "byte " + std::to_string(ch);
            throw CoalError(std::string(which) + " sequence: " + shown + " at position " + std::to_string(i + 1) +
                                " is not an amino acid letter or '*' (give ungapped protein sequences)",
                            "align_input");
        }
        out[i] = c;
    }
    return out;
}

std::string upper(const std::string& s) {
    std::string r(s);
    for (char& ch : r) ch = (char)std::toupper((unsigned char)ch);
    return r;
}

constexpr int NEG = INT_MIN / 4;

// One row of the kernel, j = 1..m, in three loops; the first and third have no loop-carried dependency.
// 1. X, M, D and the X / D traceback bits.
COAL_INLINE void row_pass1(const int* __restrict hp, const int* __restrict xp, const int* __restrict prof,
                           int* __restrict xc, int* __restrict dc, int* __restrict t, int m, int open, int ext) {
    for (int j = 1; j <= m; ++j) {
        int xo = hp[j] - open, xe = xp[j] - ext;
        int x = xe > xo ? xe : xo;
        int mm = hp[j - 1] + prof[j];
        int d = x > mm ? x : mm;
        xc[j] = x;
        dc[j] = d;
        t[j] = (xe > xo) | ((x > mm) << 1);
    }
}

// 2. q[j] = max_{k <= j} (D[k] - (j - k) * ext) for j = 0..m-1 (D[0] = H[i][0]), so that Y[j] = q[j-1] - open:
//    Y[j] = max(D[j-1] - open, Y[j-1] - ext) unrolled. The only loop-carried dependency of the kernel.
COAL_INLINE void row_pass2_scalar(const int* __restrict d, int* __restrict q, int m, int ext) {
    int run = NEG;
    for (int j = 0; j < m; ++j) {
        run = run - ext > d[j] ? run - ext : d[j];
        q[j] = run;
    }
}

#if COAL_ALIGN_AVX2
// The same eight lanes at a time with AVX2: a log-step max-plus scan in the register, then the carry from the
// blocks before (2.5x faster than the scalar loop; on SSE2 the emulated shuffles make it slower, so baseline
// builds use the scalar loop).
typedef int v8i __attribute__((vector_size(32)));
#if defined(__clang__)
#define COAL_SHIFT(v, f, k) __builtin_shufflevector(v, f, SHIFT##k)
#define SHIFT1 8, 0, 1, 2, 3, 4, 5, 6
#define SHIFT2 8, 9, 0, 1, 2, 3, 4, 5
#define SHIFT4 8, 9, 10, 11, 0, 1, 2, 3
#define COAL_LAST(v) __builtin_shufflevector(v, v, 7, 7, 7, 7, 7, 7, 7, 7)
#else
#define COAL_SHIFT(v, f, k) __builtin_shuffle(v, f, shift##k)
#define COAL_LAST(v) __builtin_shuffle(v, last7)
#endif
#define COAL_VMAX(a, b) (((a) & ((a) > (b))) | ((b) & ~((a) > (b))))
COAL_INLINE void row_pass2_avx2(const int* __restrict d, int* __restrict q, int m, int ext) {
    const v8i fill = {NEG, NEG, NEG, NEG, NEG, NEG, NEG, NEG}, zero = {0, 0, 0, 0, 0, 0, 0, 0};
    const v8i lane = {1, 2, 3, 4, 5, 6, 7, 8};
#if !defined(__clang__)
    const v8i shift1 = {8, 0, 1, 2, 3, 4, 5, 6}, shift2 = {8, 9, 0, 1, 2, 3, 4, 5}, shift4 = {8, 9, 10, 11, 0, 1, 2, 3};
    const v8i last7 = {7, 7, 7, 7, 7, 7, 7, 7};
#endif
    const v8i e1 = zero + ext, e2 = zero + 2 * ext, e4 = zero + 4 * ext, ecarry = lane * ext;
    v8i carry = fill;  // q of the block before, in every lane
    int j = 0;
    for (; j + 8 <= m; j += 8) {
        v8i v, w;
        std::memcpy(&v, d + j, sizeof v);
        w = COAL_SHIFT(v, fill, 1) - e1;
        v = COAL_VMAX(v, w);
        w = COAL_SHIFT(v, fill, 2) - e2;
        v = COAL_VMAX(v, w);
        w = COAL_SHIFT(v, fill, 4) - e4;
        v = COAL_VMAX(v, w);
        w = carry - ecarry;
        v = COAL_VMAX(v, w);
        std::memcpy(q + j, &v, sizeof v);
        carry = COAL_LAST(v);
    }
    int run = carry[0];
    for (; j < m; ++j) {
        run = run - ext > d[j] ? run - ext : d[j];
        q[j] = run;
    }
}
#undef COAL_VMAX
#undef COAL_SHIFT
#undef COAL_LAST
#endif

// 3. Y, H and the Y / H traceback bits; packs the row's traceback two cells per byte: cell j <= h (h = ceil(m/2))
//    in the low nibble of byte j-1, cell j > h in the high nibble of byte j-1-h (contiguous halves: no shuffles).
COAL_INLINE void row_pass3(const int* __restrict dc, const int* __restrict q, int* __restrict hc, int* __restrict t,
                           uint8_t* __restrict out, int m, int open) {
    for (int j = 1; j <= m; ++j) {
        int y = q[j - 1] - open;
        int yo = dc[j - 1] - open;  // dc[0] = H[i][0]
        int d = dc[j];
        int hb = y > d;
        hc[j] = hb ? y : d;
        t[j] = t[j] | ((y > yo) << 2) | (hb << 3);
    }
    const int h = (m + 1) / 2, lo = m - h;
    for (int k = 0; k < lo; ++k) out[k] = (uint8_t)(t[k + 1] | (t[k + 1 + h] << 4));
    if (h > lo) out[h - 1] = (uint8_t)t[h];
}

// The forward pass over all rows: fills the traceback (n rows of `stride` bytes), leaves row n in hp and the last
// column in lastcol.
struct Dp {
    const uint8_t* a;
    size_t n, m;
    int* const* prof;
    int open, ext;
    bool free_ends;
    uint8_t* tb;
    size_t stride;
    int *hp, *xp, *hc, *xc, *dc, *q, *t, *lastcol;
};

template <bool AVX2>
COAL_INLINE void forward_body(Dp& p) {
    const int mi = (int)p.m, o = p.open, e = p.ext;
    for (size_t i = 1; i <= p.n; ++i) {
        const int h0 = p.free_ends ? 0 : -(o + (int)(i - 1) * e);
        p.hc[0] = h0;
        p.xc[0] = NEG;
        p.dc[0] = h0;
        row_pass1(p.hp, p.xp, p.prof[p.a[i - 1]], p.xc, p.dc, p.t, mi, o, e);
#if COAL_ALIGN_AVX2
        if (AVX2) row_pass2_avx2(p.dc, p.q, mi, e);
        else row_pass2_scalar(p.dc, p.q, mi, e);
#else
        row_pass2_scalar(p.dc, p.q, mi, e);
#endif
        row_pass3(p.dc, p.q, p.hc, p.t, p.tb + (i - 1) * p.stride, mi, o);
        p.lastcol[i] = p.hc[p.m];
        std::swap(p.hp, p.hc);
        std::swap(p.xp, p.xc);
    }
}

void forward_base(Dp& p) { forward_body<false>(p); }
#if COAL_ALIGN_AVX2
COAL_TARGET_AVX2 void forward_avx2(Dp& p) { forward_body<true>(p); }
bool cpu_has_avx2() {
    static const bool yes = [] {
        __builtin_cpu_init();
        return __builtin_cpu_supports("avx2") != 0;
    }();
    return yes;
}
#endif

void forward(Dp& p) {
#if COAL_ALIGN_AVX2
    if (cpu_has_avx2()) return forward_avx2(p);
#endif
    forward_base(p);
}

struct Raw {
    std::string a, b;
    long long score = 0;  // scaled
};

Raw kernel(const std::string& sa, const std::string& sb, const std::vector<uint8_t>& a, const std::vector<uint8_t>& b,
           const Scaled& s, bool free_ends, uint64_t max_cells) {
    const size_t n = a.size(), m = b.size();
    const int o = s.open, e = s.ext;
    Raw r;
    auto gap_cost = [&](size_t k) -> long long { return k ? (long long)o + (long long)(k - 1) * e : 0; };
    if (n == 0 || m == 0) {
        r.a = upper(sa) + std::string(m, '-');
        r.b = std::string(n, '-') + upper(sb);
        r.score = free_ends ? 0 : -gap_cost(n + m);
        return r;
    }
    if (max_cells && (uint64_t)n * (uint64_t)m > max_cells)
        throw CoalError("pair too large to align: " + std::to_string(n) + " x " + std::to_string(m) +
                            " residues is more than max_cells = " + std::to_string(max_cells) + " (" +
                            std::to_string(max_cells / 2 / 1000000) + " MB of traceback)",
                        "align_too_large");
    // int32 headroom: every score is within (n + m + 2) * max_abs of 0
    if ((double)(n + m + 2) * (double)s.max_abs * 1.0 > (double)(INT_MAX / 8))
        throw CoalError("pair too long for 32-bit scores at this scale (" + std::to_string(n) + " + " +
                            std::to_string(m) + " residues)",
                        "align_too_large");

    // query profile: prof[c][j] = s(c, b[j-1]) for the residue types c that occur in a
    bool used[24] = {};
    for (uint8_t c : a) used[c] = true;
    std::vector<int> prof_store;
    int* prof[24] = {};
    {
        size_t k = 0;
        for (int c = 0; c < 24; ++c) k += used[c];
        prof_store.assign(k * (m + 1), 0);
        size_t at = 0;
        for (int c = 0; c < 24; ++c) {
            if (!used[c]) continue;
            prof[c] = prof_store.data() + at * (m + 1);
            for (size_t j = 1; j <= m; ++j) prof[c][j] = s.sub[c][b[j - 1]];
            ++at;
        }
    }
    const size_t stride = (m + 1) / 2;
    std::unique_ptr<uint8_t[]> tb(new uint8_t[n * stride]);
    std::vector<int> hp(m + 1), xp(m + 1), hc(m + 1), xc(m + 1), dc(m + 1), q(m + 1), t(m + 1);
    std::vector<int> lastcol(n + 1);
    hp[0] = 0;
    for (size_t j = 1; j <= m; ++j) {
        hp[j] = free_ends ? 0 : -(int)gap_cost(j);
        xp[j] = NEG;
    }
    lastcol[0] = hp[m];
    Dp dp{a.data(), n, m, prof, o, e, free_ends, tb.get(), stride, hp.data(), xp.data(), hc.data(), xc.data(),
          dc.data(), q.data(), t.data(), lastcol.data()};
    forward(dp);
    const int* last = dp.hp;  // row n
    size_t bi = n, bj = m;
    long long best = last[m];
    if (free_ends) {
        for (size_t j = 0; j < m; ++j)
            if (last[j] > best) best = last[j], bi = n, bj = j;
        for (size_t i = 0; i < n; ++i)
            if (lastcol[i] > best) best = lastcol[i], bi = i, bj = m;
    }
    r.score = best;
    std::string ra, rb;
    ra.reserve(n + m);
    rb.reserve(n + m);
    // trailing overhang (free end gaps), written reversed
    for (size_t j = m; j > bj; --j) ra.push_back('-'), rb.push_back(sb[j - 1]);
    for (size_t i = n; i > bi; --i) ra.push_back(sa[i - 1]), rb.push_back('-');
    size_t i = bi, j = bj;
    int state = 0;  // 0 = H, 1 = X, 2 = Y
    auto nib = [&](size_t ii, size_t jj) -> unsigned {
        const uint8_t* rowp = tb.get() + (ii - 1) * stride;
        return jj <= stride ? (rowp[jj - 1] & 15u) : (unsigned)(rowp[jj - 1 - stride] >> 4);
    };
    while (i > 0 && j > 0) {
        unsigned c = nib(i, j);
        if (state == 0) {
            if (c & 8u) state = 2;
            else if (c & 2u) state = 1;
            else {
                ra.push_back(sa[i - 1]);
                rb.push_back(sb[j - 1]);
                --i, --j;
                continue;
            }
        }
        if (state == 1) {
            ra.push_back(sa[i - 1]);
            rb.push_back('-');
            state = (c & 1u) ? 1 : 0;
            --i;
        } else {
            ra.push_back('-');
            rb.push_back(sb[j - 1]);
            state = (c & 4u) ? 2 : 0;
            --j;
        }
    }
    for (; i > 0; --i) ra.push_back(sa[i - 1]), rb.push_back('-');
    for (; j > 0; --j) ra.push_back('-'), rb.push_back(sb[j - 1]);
    std::reverse(ra.begin(), ra.end());
    std::reverse(rb.begin(), rb.end());
    r.a = upper(ra);
    r.b = upper(rb);
    return r;
}

void fill_stats(PairAlignment& out, const Scaled& s) {
    const Codes& k = codes();
    out.length = out.a.size();
    out.identities = out.similarities = out.gaps = 0;
    for (size_t c = 0; c < out.a.size(); ++c) {
        char x = out.a[c], y = out.b[c];
        if (x == '-' && y == '-') continue;
        if (x == '-' || y == '-') {
            ++out.gaps;
            continue;
        }
        if (x == y) ++out.identities;
        if (s.m->score[k.c[(unsigned char)x]][k.c[(unsigned char)y]] > 0) ++out.similarities;
    }
}

long long score_rows_scaled(const std::string& ra, const std::string& rb, const Scaled& s, bool free_ends) {
    if (ra.size() != rb.size())
        throw CoalError("alignment rows differ in length (" + std::to_string(ra.size()) + " and " +
                            std::to_string(rb.size()) + ")",
                        "align_input");
    const Codes& k = codes();
    const size_t L = ra.size();
    auto check = [&](char ch, const char* which, size_t col) {
        if (ch != '-' && k.c[(unsigned char)ch] == 255)
            throw CoalError(std::string(which) + " row: '" + ch + "' at column " + std::to_string(col + 1) +
                                " is not an amino acid letter, '*' or '-'",
                            "align_input");
    };
    size_t firstA = L, lastA = 0, firstB = L, lastB = 0;
    for (size_t c = 0; c < L; ++c) {
        check(ra[c], "first", c);
        check(rb[c], "second", c);
        if (ra[c] != '-') firstA = std::min(firstA, c), lastA = c;
        if (rb[c] != '-') firstB = std::min(firstB, c), lastB = c;
    }
    long long total = 0;
    int prev = 0;  // 0 = pair, 1 = gap in b (a residue), 2 = gap in a
    for (size_t c = 0; c < L; ++c) {
        char x = ra[c], y = rb[c];
        if (x == '-' && y == '-') continue;
        if (x != '-' && y != '-') {
            total += s.sub[k.c[(unsigned char)x]][k.c[(unsigned char)y]];
            prev = 0;
            continue;
        }
        int st = (y == '-') ? 1 : 2;
        // the gapped row is b for st 1, a for st 2
        bool end = (st == 1) ? (firstB == L || c < firstB || c > lastB) : (firstA == L || c < firstA || c > lastA);
        if (!(end && free_ends)) total -= (prev == st) ? s.ext : s.open;
        prev = st;
    }
    return total;
}

}  // namespace

const SubstitutionMatrix* substitution_matrix(const std::string& name) {
    std::string u = upper_name(name);
    for (const auto& m : MATRICES)
        if (u == m.name) return &m;
    return nullptr;
}

const std::vector<std::string>& substitution_matrix_names() {
    static const std::vector<std::string> names = [] {
        std::vector<std::string> v;
        for (const auto& m : MATRICES) v.emplace_back(m.name);
        return v;
    }();
    return names;
}

void check_align_params(const AlignParams& p) {
    if (!substitution_matrix(p.matrix)) {
        std::string known;
        for (const auto& n : substitution_matrix_names()) known += (known.empty() ? "" : ", ") + n;
        throw CoalError("unknown substitution matrix '" + p.matrix + "' (known: " + known + ")", "align_params");
    }
    if (!std::isfinite(p.gap_open) || !std::isfinite(p.gap_extend) || p.gap_extend < 0 || p.gap_open < p.gap_extend ||
        p.gap_open > 1000)
        throw CoalError("gap costs must satisfy 1000 >= gap_open >= gap_extend >= 0 (given open " +
                            std::to_string(p.gap_open) + ", extend " + std::to_string(p.gap_extend) + ")",
                        "align_params");
    for (double v : {p.gap_open, p.gap_extend}) {
        double k = v * 1000;
        if (std::fabs(k - std::round(k)) > 1e-6)
            throw CoalError("gap costs must be multiples of 0.001 (given " + std::to_string(v) + ")", "align_params");
    }
}

PairAlignment align_pair(const std::string& a, const std::string& b, const AlignParams& p) {
    Scaled s = scale_params(p);
    std::vector<uint8_t> ea = encode(a, "first"), eb = encode(b, "second");
    Raw r = kernel(a, b, ea, eb, s, p.free_end_gaps, p.max_cells);
    PairAlignment out;
    out.a = std::move(r.a);
    out.b = std::move(r.b);
    out.score = (double)r.score / s.scale;
    fill_stats(out, s);
    return out;
}

double score_alignment(const std::string& row_a, const std::string& row_b, const AlignParams& p) {
    Scaled s = scale_params(p);
    return (double)score_rows_scaled(row_a, row_b, s, p.free_end_gaps) / s.scale;
}

std::vector<PairAlignment> align_pairs(const std::vector<std::pair<std::string, std::string>>& pairs,
                                       const AlignParams& p, int jobs) {
    check_align_params(p);
    std::vector<PairAlignment> out(pairs.size());
    std::vector<size_t> idx(pairs.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    // longest pairs first, so one giant pair does not start last and leave the other threads idle
    std::stable_sort(idx.begin(), idx.end(), [&](size_t x, size_t y) {
        return pairs[x].first.size() * pairs[x].second.size() > pairs[y].first.size() * pairs[y].second.size();
    });
    auto work = [&](size_t& i) {
        try {
            out[i] = align_pair(pairs[i].first, pairs[i].second, p);
        } catch (const std::exception& ex) {
            out[i] = PairAlignment();
            out[i].error = ex.what();
        } catch (...) {
            out[i] = PairAlignment();
            out[i].error = "unknown error";
        }
    };
    parallel_batch(idx, work, resolve_jobs(jobs));
    return out;
}

PairAlignment expand_local_hit(const std::string& query, const std::string& subject, const LocalHit& h,
                               Flanks flanks, const AlignParams& p) {
    Scaled s = scale_params(p);
    auto fail = [](const std::string& msg) { throw CoalError("local hit: " + msg, "local_hit"); };
    const long long nq = (long long)query.size(), ns = (long long)subject.size();
    if (h.qstart < 1 || h.qend < h.qstart || h.qend > nq)
        fail("query coordinates " + std::to_string(h.qstart) + ".." + std::to_string(h.qend) +
             " do not fit a query of " + std::to_string(nq) + " residues");
    if (h.sstart < 1 || h.send < h.sstart || h.send > ns)
        fail("subject coordinates " + std::to_string(h.sstart) + ".." + std::to_string(h.send) +
             " do not fit a subject of " + std::to_string(ns) + " residues");
    if (h.qseq_gapped.size() != h.sseq_gapped.size())
        fail("qseq_gapped and sseq_gapped differ in length (" + std::to_string(h.qseq_gapped.size()) + " and " +
             std::to_string(h.sseq_gapped.size()) + ")");
    auto residues = [](const std::string& g) {
        long long k = 0;
        for (char c : g) k += c != '-';
        return k;
    };
    if (residues(h.qseq_gapped) != h.qend - h.qstart + 1)
        fail("qseq_gapped has " + std::to_string(residues(h.qseq_gapped)) + " residues, qstart..qend spans " +
             std::to_string(h.qend - h.qstart + 1));
    if (residues(h.sseq_gapped) != h.send - h.sstart + 1)
        fail("sseq_gapped has " + std::to_string(residues(h.sseq_gapped)) + " residues, sstart..send spans " +
             std::to_string(h.send - h.sstart + 1));
    // core block: gap pattern from the hit, residues from the full sequences
    std::string ca, cb;
    ca.reserve(h.qseq_gapped.size());
    cb.reserve(h.qseq_gapped.size());
    size_t qi = (size_t)h.qstart - 1, si = (size_t)h.sstart - 1;
    auto same = [](char hit, char full) {
        int x = std::toupper((unsigned char)hit), y = std::toupper((unsigned char)full);
        return x == y || x == 'X';
    };
    for (size_t c = 0; c < h.qseq_gapped.size(); ++c) {
        char x = h.qseq_gapped[c], y = h.sseq_gapped[c];
        if (x == '-' && y == '-') continue;
        if (x != '-') {
            if (!same(x, query[qi]))
                fail(std::string("qseq_gapped has '") + x + "' where the query has '" + query[qi] + "' (position " +
                     std::to_string(qi + 1) + ")");
            ca.push_back(query[qi++]);
        } else {
            ca.push_back('-');
        }
        if (y != '-') {
            if (!same(y, subject[si]))
                fail(std::string("sseq_gapped has '") + y + "' where the subject has '" + subject[si] +
                     "' (position " + std::to_string(si + 1) + ")");
            cb.push_back(subject[si++]);
        } else {
            cb.push_back('-');
        }
    }
    const std::string ql = query.substr(0, (size_t)h.qstart - 1), sl = subject.substr(0, (size_t)h.sstart - 1);
    const std::string qr = query.substr((size_t)h.qend), sr = subject.substr((size_t)h.send);
    PairAlignment out;
    if (flanks == Flanks::Gap) {
        out.a = ql + std::string(sl.size(), '-') + ca + qr + std::string(sr.size(), '-');
        out.b = std::string(ql.size(), '-') + sl + cb + std::string(qr.size(), '-') + sr;
    } else {
        PairAlignment L = align_pair(ql, sl, p), R = align_pair(qr, sr, p);
        out.a = L.a + ca + R.a;
        out.b = L.b + cb + R.b;
    }
    out.a = upper(out.a);
    out.b = upper(out.b);
    // checks the residues of the flanks too
    out.score = (double)score_rows_scaled(out.a, out.b, s, p.free_end_gaps) / s.scale;
    fill_stats(out, s);
    return out;
}

}  // namespace coal
