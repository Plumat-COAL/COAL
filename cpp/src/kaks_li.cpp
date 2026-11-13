// SPDX-License-Identifier: GPL-3.0-only
//
// Ka/Ks engine (kaks.hpp): the method of Li (1993) as the R package seqinr computes it, Method::Li.
//
// This file is derived from src/kaks.c of the R package seqinr, version 4.2-44
// (https://cran.r-project.org/package=seqinr, https://github.com/lbbe-software/seqinr): the functions kaks(),
// num(), fastlwl(), catsite(), transf(), titv1(), titv2(), titv3() and prefastlwl(). kaks.c carries no copyright
// notice of its own; its comments name M. Gouy ("ajoute par M. Gouy", "MG") and J. R. Lobry ("JRL", "JLO") as
// authors of parts of it, and seqinr's documentation (?kaks, author field "D. Charif, J.R. Lobry") says: "The C
// programm implementing this method was provided by Manolo Gouy. More info is needed here to trace back the
// original C source so as to credit correct source." The package as a whole:
//   seqinr: Biological Sequences Retrieval and Analysis. Authors: Delphine Charif [aut], Olivier Clerc [ctb],
//   Carolin Frank [ctb], Jean R. Lobry [aut, cph], Anamaria Necsulea [ctb], Leonor Palmeira [ctb], Simon Penel
//   [cre], Guy Perriere [ctb].
//   License: GPL (>= 2)   (the DESCRIPTION file of seqinr 4.2-44)
// This file is distributed under the GNU General Public License, version 3 (COAL as a whole is GPL-3.0-only; see
// LICENSE and NOTICE), which "GPL (>= 2)" allows. It is a modified version of that code, not an independent
// implementation: the functions follow the upstream ones in structure, constants (the amino-acid similarity
// matrix, the weights 0.333333) and order of floating-point operations, so that the numbers are the same.
//
// Modified for COAL by Chloe Newton and Keanu Murati, 2026. What was changed:
//   - one pair at a time on COAL's codon indices, no R API; the 64 x 64 tables of prefastlwl() are built once,
//     into an immutable function-local static, so pairs can be computed on several threads;
//   - LiFlavour::coal (the default): COAL's policy is applied to the pair first (prepare(): a terminal stop
//     codon is dropped, columns with an internal stop are skipped, U is read as T); LiFlavour::seqinr compares
//     the pair as seqinr::kaks does (stop codons are compared as codons; a codon with any letter other than
//     A, C, G, T is removed with its column);
//   - seqinr's sentinel 9.999999 for Ka, Ks and their variances ("saturated", also returned when a class of
//     sites is empty) is NA here, with status `saturated` or `no_data`;
//   - standard errors (the square roots of seqinr's vka and vks) are returned; a negative Ka or Ks (which the
//     formulas give for some short pairs, and seqinr returns) is returned with status `warning`;
//   - gaps are removed per pair (seqinr removes every column in which any sequence of the alignment has a gap;
//     for two sequences that is the same).
//
// References: Li W-H (1993) Unbiased estimation of the rates of synonymous and nonsynonymous substitution. J Mol
// Evol 36:96-99. Charif D, Lobry JR (2007) SeqinR 1.0-2, in: Structural approaches to sequence evolution,
// Springer, pp 207-232.
#include <cmath>
#include <cstring>
#include <string>

#include "kaks.hpp"

namespace coal {
namespace kaks {

namespace {

// seqinr's amino-acid similarity matrix (kaks(): mat[19][19]); rl[i][j] below is built from it as upstream
const double MAT[19][19] = {
    {.382, .382, .343, .382, .382, .382, .382, .128, .040, .128, .040, .128, .040, .128, .040, .128, .343, .128, .040},
    {.382, .382, .128, .343, .343, .343, .343, .128, .040, .128, .040, .128, .040, .128, .040, .128, .128, .040, .040},
    {.343, .128, .343, .382, .382, .382, .343, .128, .040, .128, .128, .343, .128, .343, .128, .343, .343, .128, .040},
    {.382, .343, .382, .343, .343, .343, .343, .343, .040, .343, .343, .382, .343, .382, .343, .382, .382, .382, .343},
    {.382, .343, .382, .343, .382, .382, .382, .343, .040, .343, .128, .343, .128, .128, .128, .343, .343, .128, .040},
    {.382, .343, .382, .343, .382, .382, .382, .343, .040, .343, .128, .343, .128, .128, .040, .128, .128, .128, .040},
    {.382, .343, .343, .343, .382, .382, .382, .343, .040, .343, .128, .343, .128, .128, .128, .128, .343, .128, .040},
    {.128, .128, .128, .343, .343, .343, .343, .343, .040, .343, .128, .343, .128, .343, .128, .343, .343, .128, .040},
    {.040, .040, .040, .040, .040, .040, .040, .040, .040, .382, .382, .382, .343, .343, .343, .128, .128, .343, .128},
    {.128, .128, .128, .343, .343, .343, .343, .343, .382, .040, .040, .128, .128, .040, .128, .040, .040, .040, .040},
    {.040, .040, .128, .343, .128, .128, .128, .128, .382, .040, .343, .343, .343, .343, .128, .128, .128, .128, .128},
    {.128, .128, .343, .382, .343, .343, .343, .343, .382, .128, .343, .343, .343, .343, .343, .128, .128, .343, .343},
    {.040, .040, .128, .343, .128, .128, .128, .128, .343, .128, .343, .343, .343, .382, .343, .343, .343, .343, .343},
    {.128, .128, .343, .382, .128, .128, .128, .343, .343, .040, .343, .343, .382, .343, .382, .128, .128, .343, .343},
    {.040, .040, .128, .343, .128, .040, .128, .128, .343, .128, .128, .343, .343, .382, .382, .343, .382, .382, .343},
    {.128, .128, .343, .382, .343, .128, .128, .343, .128, .040, .128, .128, .343, .128, .343, .343, .343, .382, .382},
    {.343, .128, .343, .382, .343, .128, .343, .343, .128, .040, .128, .128, .343, .128, .382, .343, .382, .343, .128},
    {.128, .040, .128, .382, .128, .128, .128, .128, .343, .040, .128, .343, .343, .343, .382, .382, .343, .343, .343},
    {.040, .040, .040, .343, .040, .040, .040, .040, .128, .040, .128, .343, .343, .343, .343, .382, .128, .343, .382}};

// seqinr's number of each codon's amino acid (prefastlwl(): aa[]; codons in the order A C G T, 0 = stop), for the
// standard code (upstream's code_mt = 0)
const int AA_OF[64] = {17, 16, 17, 16, 13, 13, 13, 13, 18, 14, 18, 14, 7,  7,  5,  7,  15, 4,  15, 4,  9,  9,
                       9,  9,  18, 18, 18, 18, 6,  6,  6,  6,  19, 20, 19, 20, 11, 11, 11, 11, 12, 12, 12, 12,
                       8,  8,  8,  8,  0,  3,  0,  3,  14, 14, 14, 14, 0,  10, 2,  10, 6,  1,  6,  1};

int num(const char* cod) {
    int n1 = 0, n2 = 0, n3 = 0;
    if (cod[0] == 'C') n1 = 1;
    if (cod[1] == 'C') n2 = 1;
    if (cod[2] == 'C') n3 = 1;
    if (cod[0] == 'G') n1 = 2;
    if (cod[1] == 'G') n2 = 2;
    if (cod[2] == 'G') n3 = 2;
    if (cod[0] == 'T') n1 = 3;
    if (cod[1] == 'T') n2 = 3;
    if (cod[2] == 'T') n3 = 3;
    return 16 * n1 + 4 * n2 + n3;
}

// 0: site i (1..3) of the codon is non-degenerate; 1: 2-fold degenerate; 2: 4-fold degenerate
int catsite(char c1, char c2, char c3, int i) {
    if (i == 3) {
        if ((c1 == 'A') && (c2 == 'T') && (c3 == 'G')) return 0;
        if ((c1 == 'T') && (c2 == 'G') && (c3 == 'A')) return 0;
        if ((c1 == 'T') && (c2 == 'G') && (c3 == 'G')) return 0;
        if (c2 == 'C') return 2;
        if ((c1 == 'C') && (c2 == 'T')) return 2;
        if ((c1 == 'G') && (c2 == 'T')) return 2;
        if ((c1 == 'G') && (c2 == 'G')) return 2;
        if ((c1 == 'C') && (c2 == 'G')) return 2;
        return 1;
    } else if (i == 1) {
        if ((c1 == 'C') && (c2 == 'T') && (c3 == 'A')) return 1;
        if ((c1 == 'C') && (c2 == 'T') && (c3 == 'G')) return 1;
        if ((c1 == 'T') && (c2 == 'T') && (c3 == 'A')) return 1;
        if ((c1 == 'T') && (c2 == 'T') && (c3 == 'G')) return 1;
        if ((c1 == 'A') && (c2 == 'G') && (c3 == 'A')) return 1;
        if ((c1 == 'A') && (c2 == 'G') && (c3 == 'G')) return 1;
        if ((c1 == 'C') && (c2 == 'G') && (c3 == 'A')) return 1;
        if ((c1 == 'C') && (c2 == 'G') && (c3 == 'G')) return 1;
        return 0;
    }
    return 0;
}

// 'i' transition, 'v' transversion ('S' for two equal bases; upstream prints a message then)
char transf(char nt1, char nt2) {
    if (nt1 == nt2) return 'S';
    if ((nt1 == 'A' && nt2 == 'G') || (nt1 == 'G' && nt2 == 'A') || (nt1 == 'C' && nt2 == 'T') ||
        (nt1 == 'T' && nt2 == 'C'))
        return 'i';
    return 'v';
}

inline bool is3(const char* c, char a, char b, char d) { return c[0] == a && c[1] == b && c[2] == d; }

void titv1(const char* cod1, const char* cod2, double poids, double* ti, double* tv, double* l) {
    const char ci1 = cod1[0], ci2 = cod1[1], ci3 = cod1[2];
    const char cj1 = cod2[0], cj2 = cod2[1], cj3 = cod2[2];
    for (int i = 0; i <= 2; i++)
        if (cod1[i] != cod2[i]) {
            l[catsite(ci1, ci2, ci3, i + 1)] += 0.5 * poids;
            l[catsite(cj1, cj2, cj3, i + 1)] += 0.5 * poids;
            const char a = cod1[i], b = cod2[i];
            if (transf(a, b) == 'i') {
                ti[catsite(ci1, ci2, ci3, i + 1)] += 0.5 * poids;
                ti[catsite(cj1, cj2, cj3, i + 1)] += 0.5 * poids;
            } else {
                tv[catsite(ci1, ci2, ci3, i + 1)] += 0.5 * poids;
                tv[catsite(cj1, cj2, cj3, i + 1)] += 0.5 * poids;
            }
            if (((ci2 == 'T') && (cj2 == 'T')) || ((ci2 == 'G') && (cj2 == 'G'))) {
                if (i == 0) {
                    // transitions at a 2-fold site that are nonsynonymous: from the 2-fold transitions to the
                    // 2-fold transversions
                    if ((is3(cod1, 'C', 'G', 'A') && is3(cod2, 'T', 'G', 'A')) ||
                        (is3(cod1, 'C', 'G', 'G') && is3(cod2, 'T', 'G', 'G')) ||
                        (is3(cod1, 'A', 'G', 'G') && is3(cod2, 'G', 'G', 'G')) ||
                        (is3(cod1, 'A', 'G', 'A') && is3(cod2, 'G', 'G', 'A')) ||
                        (is3(cod1, 'T', 'G', 'A') && is3(cod2, 'C', 'G', 'A')) ||
                        (is3(cod1, 'T', 'G', 'G') && is3(cod2, 'C', 'G', 'G')) ||
                        (is3(cod1, 'G', 'G', 'G') && is3(cod2, 'A', 'G', 'G')) ||
                        (is3(cod1, 'G', 'G', 'A') && is3(cod2, 'A', 'G', 'A'))) {
                        ti[1] -= 0.5 * poids;
                        tv[1] += 0.5 * poids;
                    }
                    // synonymous transversions at 2-fold sites: from the transversions to the transitions
                    if ((is3(cod1, 'C', 'G', 'A') && is3(cod2, 'A', 'G', 'A')) ||
                        (is3(cod1, 'A', 'G', 'A') && is3(cod2, 'C', 'G', 'A')) ||
                        (is3(cod1, 'C', 'G', 'G') && is3(cod2, 'A', 'G', 'G')) ||
                        (is3(cod1, 'A', 'G', 'G') && is3(cod2, 'C', 'G', 'G'))) {
                        tv[1] -= poids;
                        ti[1] += poids;
                    }
                }
                if (i == 2) {
                    if ((is3(cod1, 'A', 'T', 'A') && is3(cod2, 'A', 'T', 'T')) ||
                        (is3(cod1, 'A', 'T', 'T') && is3(cod2, 'A', 'T', 'A')) ||
                        (is3(cod1, 'A', 'T', 'A') && is3(cod2, 'A', 'T', 'C')) ||
                        (is3(cod1, 'A', 'T', 'C') && is3(cod2, 'A', 'T', 'A'))) {
                        tv[1] -= poids;
                        ti[1] += poids;
                    }
                    if ((is3(cod1, 'A', 'T', 'A') && is3(cod2, 'A', 'T', 'G')) ||
                        (is3(cod1, 'A', 'T', 'G') && is3(cod2, 'A', 'T', 'A'))) {
                        ti[1] -= 0.5 * poids;
                        tv[1] += 0.5 * poids;
                    }
                }
            }
        }
}

void titv2(const char* cod1, const char* cod2, double* ti, double* tv, double* l, const double rl[21][21],
           const int* pos) {
    char codint1[3], codint2[3];
    int i, j, n = 0;
    std::memcpy(codint1, cod1, 3);
    std::memcpy(codint2, cod1, 3);
    for (i = 0; i < 2; i++) {
        if (cod1[i] != cod2[i]) {
            codint1[i] = cod2[i];
            break;
        }
    }
    for (j = i + 1; j <= 2; j++) {
        if (cod1[j] != cod2[j]) {
            codint2[j] = cod2[j];
            break;
        }
    }
    const int aa1 = AA_OF[num(cod1)], aa2 = AA_OF[num(cod2)];
    const int aaint1 = AA_OF[num(codint1)], aaint2 = AA_OF[num(codint2)];
    const double l1 = rl[aa1][aaint1] * rl[aaint1][aa2];
    const double l2 = rl[aa1][aaint2] * rl[aaint2][aa2];
    const double p1 = (l1 + l2) ? l1 / (l1 + l2) : 0.;
    const double p2 = (l1 + l2) ? 1. - p1 : 0.;
    for (i = 0; i < 3; i++)
        if (pos[i] == 0) n = i + 1;
    l[catsite(cod1[0], cod1[1], cod1[2], n)] += 0.333333;
    l[catsite(cod2[0], cod2[1], cod2[2], n)] += 0.333333;
    l[catsite(codint1[0], codint1[1], codint1[2], n)] += 0.333333 * p1;
    l[catsite(codint2[0], codint2[1], codint2[2], n)] += 0.333333 * p2;
    titv1(cod1, codint1, p1, ti, tv, l);
    titv1(cod2, codint1, p1, ti, tv, l);
    titv1(cod1, codint2, p2, ti, tv, l);
    titv1(cod2, codint2, p2, ti, tv, l);
}

void titv3(const char* cod1, const char* cod2, double* ti, double* tv, double* l, const double rl[21][21]) {
    char codint1[6][3], codint2[6][3];
    double like[6], p[6], somli;
    int ii;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            if (j != i) {
                if ((i == 0) || ((i == 1) && (j == 0))) {
                    ii = 3 * i + j - 1;
                } else {
                    ii = 3 * i + j - 2;
                }
                std::memcpy(codint1[ii], cod1, 3);
                codint1[ii][i] = cod2[i];
                std::memcpy(codint2[ii], codint1[ii], 3);
                codint2[ii][j] = cod2[j];
                const int aaa = AA_OF[num(cod1)], aab = AA_OF[num(codint1[ii])];
                const int aac = AA_OF[num(codint2[ii])], aad = AA_OF[num(cod2)];
                const double rlab = rl[aaa][aab], rlbc = rl[aab][aac], rlcd = rl[aac][aad];
                like[ii] = rlab * rlbc * rlcd;
            }
    }
    somli = 0;
    for (int i = 0; i < 6; i++) somli += like[i];
    for (int i = 0; i < 6; i++) {
        p[i] = like[i] / somli;
        titv1(cod1, codint1[i], p[i], ti, tv, l);
        titv1(codint1[i], codint2[i], p[i], ti, tv, l);
        titv1(codint2[i], cod2, p[i], ti, tv, l);
    }
}

// prefastlwl(): for every pair of codons (seqinr's numbering), the sites l[0..2] and the transitional and
// transversional differences ti[0..2], tv[0..2] per class (0: non-degenerate, 1: 2-fold, 2: 4-fold)
struct Tables {
    double tl[3][64][64], tti[3][64][64], ttv[3][64][64];
};

void build_tables(Tables& T) {
    double rl[21][21];
    for (int i = 0; i < 21; ++i)
        for (int j = 0; j < 21; ++j) rl[i][j] = 0;
    for (int i = 2; i < 21; i++)
        for (int j = 1; j < i; j++) rl[i][j] = MAT[j - 1][i - 2];
    for (int i = 1; i <= 20; i++) {
        rl[i][i] = 1.0;
        for (int j = i + 1; j <= 20; j++) rl[i][j] = rl[j][i];
    }
    // "ajoute par M. Gouy": a stop codon (amino acid 0) gets the smallest similarity of the table
    double minrl = rl[1][1];
    for (int i = 1; i <= 20; i++)
        for (int j = i + 1; j <= 20; j++)
            if (rl[i][j] < minrl) minrl = rl[i][j];
    for (int i = 0; i <= 20; i++) rl[0][i] = rl[i][0] = minrl;

    static const char B[] = "ACGT";
    for (int i = 0; i < 64; i++) {
        for (int j = i; j < 64; j++) {
            double l[3] = {0, 0, 0}, ti[3] = {0, 0, 0}, tv[3] = {0, 0, 0};
            const char cod1[3] = {B[i / 16], B[(i / 4) % 4], B[i % 4]};
            const char cod2[3] = {B[j / 16], B[(j / 4) % 4], B[j % 4]};
            int nbdiff = 0, pos[3] = {0, 0, 0};
            for (int k = 0; k < 3; ++k)
                if (cod1[k] != cod2[k]) {
                    nbdiff++;
                    pos[k] = 1;
                }
            if (nbdiff != 2)
                for (int jj = 0; jj < 3; jj++)
                    if (pos[jj] == 0) {
                        l[catsite(cod1[0], cod1[1], cod1[2], jj + 1)] += 0.5;
                        l[catsite(cod2[0], cod2[1], cod2[2], jj + 1)] += 0.5;
                    }
            if (nbdiff == 1) titv1(cod1, cod2, 1.0, ti, tv, l);
            if (nbdiff == 2) titv2(cod1, cod2, ti, tv, l, rl, pos);
            if (nbdiff == 3) titv3(cod1, cod2, ti, tv, l, rl);
            for (int k = 0; k < 3; ++k) {
                T.tl[k][i][j] = T.tl[k][j][i] = l[k];
                T.tti[k][i][j] = T.tti[k][j][i] = ti[k];
                T.ttv[k][i][j] = T.ttv[k][j][i] = tv[k];
            }
        }
    }
}

struct TablesHolder {
    Tables t;
    TablesHolder() { build_tables(t); }
};

const Tables& tables() {
    // built once (C++11 guarantees a thread-safe initialisation) and never written again
    static const TablesHolder h;
    return h.t;
}

// COAL's codon index (T C A G) -> seqinr's (A C G T)
inline int to_seqinr(int c) {
    static const int m[4] = {3, 1, 0, 2};
    return 16 * m[c / 16] + 4 * m[(c / 4) % 4] + m[c % 4];
}

void add_note(Result& r, const std::string& msg) { r.note = r.note.empty() ? msg : r.note + "; " + msg; }

}  // namespace

namespace detail {

Result run_li(const Prepared& pr, const GeneticCode& code, const Options& opt, Result r) {
    (void)code;  // the standard code (compute() refuses the others for this method)
    const bool raw = opt.li_flavour == LiFlavour::seqinr;
    const std::vector<uint8_t>& s1 = raw ? pr.r1 : pr.c1;
    const std::vector<uint8_t>& s2 = raw ? pr.r2 : pr.c2;
    const Tables& T = tables();

    // fastlwl()
    double l[3], a[3], b[3], p[3], q[3], ti[3], tv[3], cc[3], aaa[3], bb[3], va[3], vb[3];
    l[0] = l[1] = l[2] = 0;
    ti[0] = ti[1] = ti[2] = tv[0] = tv[1] = tv[2] = 0;
    double subs = 0;
    for (size_t ii = 0; ii < s1.size(); ii++) {
        const int num1 = to_seqinr(s1[ii]), num2 = to_seqinr(s2[ii]);
        l[0] += T.tl[0][num1][num2];
        l[1] += T.tl[1][num1][num2];
        l[2] += T.tl[2][num1][num2];
        ti[0] += T.tti[0][num1][num2];
        ti[1] += T.tti[1][num1][num2];
        ti[2] += T.tti[2][num1][num2];
        tv[0] += T.ttv[0][num1][num2];
        tv[1] += T.ttv[1][num1][num2];
        tv[2] += T.ttv[2][num1][num2];
        const int x = s1[ii] ^ s2[ii];
        subs += ((x & 0x30) != 0) + ((x & 0x0c) != 0) + ((x & 0x03) != 0);
    }
    for (int ii = 0; ii < 3; ii++) {
        p[ii] = ti[ii] / l[ii];
        q[ii] = tv[ii] / l[ii];
        aaa[ii] = 1 / (1 - 2 * p[ii] - q[ii]);
        bb[ii] = 1 / (1 - 2 * q[ii]);
        cc[ii] = (aaa[ii] + bb[ii]) / 2;
        if (bb[ii] <= 0 || !std::isfinite(bb[ii])) {
            b[ii] = 10.0;
        } else {
            b[ii] = 0.5 * (double)std::log(bb[ii]);
        }
        if ((aaa[ii] <= 0) || (bb[ii] <= 0) || !std::isfinite(aaa[ii]) || !std::isfinite(bb[ii])) {
            a[ii] = 10.0;
        } else {
            a[ii] = 0.5 * (double)std::log(aaa[ii]) - 0.25 * std::log(bb[ii]);
        }
        va[ii] = (aaa[ii] * aaa[ii] * p[ii] + cc[ii] * cc[ii] * q[ii] -
                  (aaa[ii] * p[ii] + cc[ii] * q[ii]) * (aaa[ii] * p[ii] + cc[ii] * q[ii])) /
                 l[ii];
        vb[ii] = bb[ii] * bb[ii] * q[ii] * (1 - q[ii]) / l[ii];
    }
    r.has_fold = true;
    for (int k = 0; k < 3; ++k) {
        r.fold_sites[k] = l[k];
        r.fold_s[k] = ti[k];
        r.fold_n[k] = tv[k];
    }
    r.substitutions = subs;

    auto undefined = [&](const char* which, bool empty_class, const char* empty_msg) {
        if (empty_class) {
            r.status = Status::no_data;
            add_note(r, empty_msg);
        } else {
            if (r.status != Status::no_data) r.status = Status::saturated;
            add_note(r, std::string(which) + " saturated (multiple-hit correction undefined)");
        }
    };
    // seqinr: ks = vks = 9.999999 unless a[1], a[2] and b[2] are all below 10 (10 marks an undefined logarithm)
    if ((a[1] < 10) && (a[2] < 10) && (b[2] < 10)) {
        r.ks = (l[1] * a[1] + l[2] * a[2]) / (l[2] + l[1]) + b[2];
        const double vks = (l[1] * l[1] * va[1] + l[2] * l[2] * va[2]) / ((l[1] + l[2]) * (l[1] + l[2])) + vb[2] -
                           bb[2] * q[2] * (2 * aaa[2] * p[2] - cc[2] * (1 - q[2])) / (l[1] + l[2]);
        r.se_ks = std::sqrt(vks);
    } else {
        undefined("Ks", !(l[1] > 0) || !(l[2] > 0),
                  "Li: no 2-fold or no 4-fold degenerate sites (Ks needs both; seqinr returns 9.999999)");
    }
    if ((a[0] < 10) && (b[0] < 10) && (b[1] < 10)) {
        r.ka = a[0] + (l[0] * b[0] + l[1] * b[1]) / (l[0] + l[1]);
        const double vka = (l[0] * l[0] * vb[0] + l[1] * l[1] * vb[1]) / ((l[1] + l[0]) * (l[1] + l[0])) + va[0] -
                           bb[0] * q[0] * (2 * aaa[0] * p[0] - cc[0] * (1 - q[0])) / (l[1] + l[0]);
        r.se_ka = std::sqrt(vka);
    } else {
        undefined("Ka", !(l[0] > 0) || !(l[1] > 0),
                  "Li: no non-degenerate or no 2-fold degenerate sites (Ka needs both; seqinr returns 9.999999)");
    }
    // The formulas can give a negative distance (few differences of one kind in a short pair): seqinr returns
    // it, and so does COAL, with a warning (orthologr turns a negative value into NA).
    for (const char* which : {"Ka", "Ks"}) {
        const double v = which[1] == 'a' ? r.ka : r.ks;
        if (v < 0) {
            if (r.status == Status::ok) r.status = Status::warning;
            add_note(r, std::string(which) + " is negative (Li's correction for these counts; not a distance)");
        }
    }
    return r;
}

}  // namespace detail
}  // namespace kaks
}  // namespace coal
