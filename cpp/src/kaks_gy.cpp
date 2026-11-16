// SPDX-License-Identifier: GPL-3.0-only
//
// Method GY: maximum-likelihood dN / dS of a pair under the codon model of Goldman & Yang (1994), as PAML's
// codeml estimates it for a pair (runmode = -2, CodonFreq = 2) and as KaKs_Calculator's GY (GY-HKY) does.
//
//   Q(i, j) = pi_j                      i, j one nucleotide apart, transversion, synonymous
//           = kappa pi_j                transition, synonymous
//           = omega pi_j                transversion, nonsynonymous
//           = kappa omega pi_j          transition, nonsynonymous
//           = 0                         otherwise (two or three positions; a stop codon of the genetic code)
//
// scaled to one substitution per codon at stationarity, pi = F3x4 from the two sequences (stop codons removed and
// renormalised: kaks.cpp f3x4, the frequencies YN uses). The log-likelihood of the compared codon columns is
// sum_h log(pi_i P_ij(t)); it is maximised over the branch length t, kappa and omega. dN and dS are those of Yang &
// Nielsen (2000), which codeml prints: with rhoS (rhoN) the proportion of the rate that is synonymous
// (nonsynonymous) and rhoS1 (rhoN1) the same at omega = 1, S = 3 L rhoS1, N = 3 L rhoN1, dS = t rhoS / (3 rhoS1),
// dN = t rhoN / (3 rhoN1), so dN / dS = omega.
//
// The algorithm is COAL's own (no code of codeml or KaKs_Calculator is used; the model, the bounds' roles and the
// definitions are theirs): Q is symmetrised with pi and decomposed once per (kappa, omega) (kaks.cpp's symmetric
// eigen solver), which makes the likelihood and its two derivatives in t cheap; t is profiled out by a safeguarded
// Newton iteration and (log kappa, log omega) are found by a projected BFGS iteration with the exact gradient,
// from fixed starting points (one for a close pair, four for a diverged one). Standard errors come from the
// observed information (central differences of that gradient) and the delta method. Nothing is random and
// nothing depends on the thread: the same pair gives the same numbers in any run. docs/kaks.md, "GY: maximum
// likelihood", has the statuses, the bounds, the validation against codeml and the differences from codeml and
// KaKs_Calculator.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "kaks.hpp"

namespace coal {
namespace kaks {
namespace {

inline bool is_stop(const GeneticCode& g, int c) { return g.aa[c] == '*'; }
inline int base_at(int c, int pos) { return pos == 0 ? c / 16 : pos == 1 ? (c / 4) % 4 : c % 4; }

void add_note(Result& r, const std::string& msg) { r.note = r.note.empty() ? msg : r.note + "; " + msg; }

struct GyModel {
    const GeneticCode& g;
    double pi[64];
    int ns = 0, idx[64], back[64];
    // neighbours (one substitution apart, both sense codons with pi > 0), by compact index
    int nnb[64], nb[64][9];
    bool nb_ts[64][9], nb_non[64][9];
    // data: distinct unordered codon pairs (compact indices) and their counts
    std::vector<int> pa, pb;
    std::vector<double> pn;
    double lconst = 0;  // sum n * log sqrt(pi_a pi_b)
    long ncodons = 0;
    // state for the current (kappa, omega)
    std::vector<double> A, d, e, C, mu, E;
    std::vector<double> T, W, U, Vt;   // scratch of gradient()
    double sq[64][9];              // the exchangeabilities of the neighbours at the current (kappa, omega)
    // the symmetrised, unscaled rate matrix B in sparse form, for series(): off-diagonal entries per neighbour,
    // minus the diagonal, and the parts of the diagonal that carry kappa (transitions) and omega (nonsynonymous)
    double bo[64][9], bdiag[64], bdiag_ts[64], bdiag_non[64], lambda = 0;
    double mr = 1, mr_ts = 0, rhoS = 0, rhoN = 0, rhoS1 = 0, rhoN1 = 0;
    long neval = 0;

    explicit GyModel(const GeneticCode& code) : g(code) {}

    void setup(const Prepared& p) {
        // F3x4 from both sequences, stop codons removed (the same arithmetic as kaks.cpp f3x4)
        double f[12] = {0};
        const size_t n = p.c1.size();
        ncodons = (long)n;
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
        ns = 0;
        for (int c = 0; c < 64; ++c) {
            back[c] = -1;
            if (pi[c] > 0) {
                back[c] = ns;
                idx[ns++] = c;
            }
        }
        for (int a = 0; a < ns; ++a) {
            const int i = idx[a];
            nnb[a] = 0;
            for (int pos = 0; pos < 3; ++pos)
                for (int b4 = 0; b4 < 4; ++b4) {
                    const int x = base_at(i, pos);
                    if (b4 == x) continue;
                    const int j = pos == 0 ? i + (b4 - x) * 16 : pos == 1 ? i + (b4 - x) * 4 : i + (b4 - x);
                    if (back[j] < 0) continue;
                    nb[a][nnb[a]] = back[j];
                    nb_ts[a][nnb[a]] = (x + b4 == 1 || x + b4 == 5);
                    nb_non[a][nnb[a]] = g.aa[i] != g.aa[j];
                    ++nnb[a];
                }
        }
        // patterns: unordered pairs (the model is reversible), in a fixed order
        std::vector<double> cnt((size_t)ns * ns, 0.0);
        for (size_t h = 0; h < n; ++h) {
            int a = back[p.c1[h]], b = back[p.c2[h]];
            if (a < b) std::swap(a, b);
            cnt[(size_t)a * ns + b] += 1;
        }
        lconst = 0;
        for (int a = 0; a < ns; ++a)
            for (int b = 0; b <= a; ++b) {
                const double c = cnt[(size_t)a * ns + b];
                if (c == 0) continue;
                pa.push_back(a);
                pb.push_back(b);
                pn.push_back(c);
                lconst += c * 0.5 * std::log(pi[idx[a]] * pi[idx[b]]);
            }
        A.resize((size_t)ns * ns);
        d.resize(ns);
        e.resize(ns);
        mu.resize(ns);
        E.resize(ns);
        C.resize(pa.size() * (size_t)ns);
        T.resize((size_t)ns * ns);
        W.resize((size_t)ns * ns);
        U.resize((size_t)ns * ns);
        Vt.resize((size_t)ns * ns);
    }

    // The synonymous and nonsynonymous parts of the mean rate (before scaling): syn with its kappa factors, non0
    // the nonsynonymous part without omega.
    void rates(double kappa, double& syn, double& non0) const {
        syn = non0 = 0;
        for (int a = 0; a < ns; ++a) {
            double s = 0, v = 0;
            for (int k = 0; k < nnb[a]; ++k) {
                const double q = (nb_ts[a][k] ? kappa : 1.0) * pi[idx[nb[a][k]]];
                if (nb_non[a][k])
                    v += q;
                else
                    s += q;
            }
            syn += pi[idx[a]] * s;
            non0 += pi[idx[a]] * v;
        }
    }

    // Decompose for (kappa, omega): afterwards lnl(t) is cheap.
    void build(double kappa, double omega) {
        ++neval;
        double syn, non0;
        rates(kappa, syn, non0);
        mr = syn + omega * non0;
        rhoS = syn / mr;
        rhoN = omega * non0 / mr;
        rhoS1 = syn / (syn + non0);
        rhoN1 = non0 / (syn + non0);
        std::fill(A.begin(), A.end(), 0.0);
        mr_ts = 0;
        lambda = 0;
        for (int a = 0; a < ns; ++a) {
            double row = 0, row_ts = 0, row_non = 0;
            for (int k = 0; k < nnb[a]; ++k) {
                const int b = nb[a][k];
                const double s = (nb_ts[a][k] ? kappa : 1.0) * (nb_non[a][k] ? omega : 1.0);
                sq[a][k] = s;
                bo[a][k] = s * std::sqrt(pi[idx[a]] * pi[idx[b]]);
                if (nb_ts[a][k]) mr_ts += pi[idx[a]] * s * pi[idx[b]];
                if (nb_ts[a][k]) row_ts += s * pi[idx[b]];
                if (nb_non[a][k]) row_non += s * pi[idx[b]];
                row += s * pi[idx[b]];
                if (b < a) A[(size_t)a * ns + b] = A[(size_t)b * ns + a] = s * std::sqrt(pi[idx[a]] * pi[idx[b]]);
            }
            A[(size_t)a * ns + a] = -row;
            bdiag[a] = row;
            bdiag_ts[a] = row_ts;
            bdiag_non[a] = row_non;
            lambda = std::max(lambda, row);
        }
        detail::eigen_sym_rows(A.data(), ns, d.data(), e.data());  // eigenvector k in row k
        for (int k = 0; k < ns; ++k) mu[k] = d[k] / mr;
        const size_t np = pa.size();
        for (size_t h = 0; h < np; ++h) {
            double* c = &C[h * ns];
            const int a = pa[h], b = pb[h];
            for (int k = 0; k < ns; ++k) c[k] = A[(size_t)k * ns + a] * A[(size_t)k * ns + b];
        }
    }

    // Row a of M = exp(B tau) (tau = t / mr) without cancellation, for the entries the eigen sum cannot give:
    // sum_k R_ak R_bk e^(d_k tau) carries an absolute error of about 1e-16, which is the whole value of a
    // probability below that (codons two or three substitutions apart at a small t: t^2, t^3) and a relative
    // error of 1e-8 at 1e-8. Uniformisation: B = lambda (R - I) with R = I + B / lambda >= 0 entry by entry, so
    // exp(B tau) = e^(-x) sum_k x^k / k! R^k (x = lambda tau) is a sum of non-negative terms, exact to rounding in
    // every entry however small. Only the row is needed: v_k = R v_(k-1) from the unit vector, with the sparse R.
    // M: the row; M1, M2: its first and second derivative in t; Gk, Gw: its derivatives in log kappa and log
    // omega at fixed tau (tangent vectors of the same recursion; lambda is a constant of the method).
    // False when x is too large for the series (the caller keeps the eigen value).
    struct Row {
        double M[64], M1[64], M2[64], Gk[64], Gw[64];
    };
    bool series(int a, double t, Row& r, bool with_grad) const {
        const double tau = t / mr, x = lambda * tau;
        if (!(x < 600) || !(lambda > 0)) return false;
        const int K = (int)(x + 12 * std::sqrt(x) + 30);
        double v[64], w[64], vk[64], wk[64], vw[64], ww[64];
        for (int b = 0; b < ns; ++b) {
            v[b] = b == a ? 1.0 : 0.0;
            vk[b] = vw[b] = 0;
            r.M[b] = r.M1[b] = r.M2[b] = r.Gk[b] = r.Gw[b] = 0;
        }
        double p = std::exp(-x), p1 = 0, p2 = 0;   // p_k, p_(k-1), p_(k-2): e^-x x^k / k!
        const double il = 1 / lambda;
        for (int k = 0;; ++k) {
            const double c1 = (p1 - p) * lambda / mr, c2 = (p2 - 2 * p1 + p) * lambda * lambda / (mr * mr);
            for (int b = 0; b < ns; ++b) {
                r.M[b] += p * v[b];
                r.M1[b] += c1 * v[b];
                r.M2[b] += c2 * v[b];
            }
            if (with_grad)
                for (int b = 0; b < ns; ++b) {
                    r.Gk[b] += p * vk[b];
                    r.Gw[b] += p * vw[b];
                }
            if (k >= K) break;
            // next power: v <- R v, and the tangents dv <- dR v + R dv
            for (int b = 0; b < ns; ++b) {
                double s = 0, sk = 0, sw = 0, tk = 0, tw = 0;
                for (int j = 0; j < nnb[b]; ++j) {
                    const int c = nb[b][j];
                    const double q = bo[b][j];
                    s += q * v[c];
                    if (with_grad) {
                        tk += q * vk[c];
                        tw += q * vw[c];
                        if (nb_ts[b][j]) sk += q * v[c];
                        if (nb_non[b][j]) sw += q * v[c];
                    }
                }
                const double keep = 1 - bdiag[b] * il;
                w[b] = keep * v[b] + s * il;
                if (with_grad) {
                    wk[b] = keep * vk[b] + (tk + sk - bdiag_ts[b] * v[b]) * il;
                    ww[b] = keep * vw[b] + (tw + sw - bdiag_non[b] * v[b]) * il;
                }
            }
            for (int b = 0; b < ns; ++b) {
                v[b] = w[b];
                if (with_grad) {
                    vk[b] = wk[b];
                    vw[b] = ww[b];
                }
            }
            p2 = p1;
            p1 = p;
            p *= x / (k + 1);
        }
        return true;
    }
    // below this the eigen sum is replaced by series() (its relative error there is 1e-8 and more)
    static constexpr double TINY_PROB = 1e-8;

    // log-likelihood at t with the current decomposition; g1, g2: first and second derivative in t
    double lnl(double t, double* g1 = nullptr, double* g2 = nullptr) {
        for (int k = 0; k < ns; ++k) E[k] = std::exp(mu[k] * t);
        double f = lconst, d1 = 0, d2 = 0;
        const size_t np = pa.size();
        Row row;
        int row_of = -1;   // the source codon `row` holds (patterns come sorted by pa)
        bool row_ok = false;
        for (size_t h = 0; h < np; ++h) {
            const double* c = &C[h * ns];
            double s0 = 0, s1 = 0, s2 = 0;
            for (int k = 0; k < ns; ++k) {
                const double x = c[k] * E[k];
                s0 += x;
                s1 += x * mu[k];
                s2 += x * mu[k] * mu[k];
            }
            if (!(s0 > TINY_PROB)) {
                if (row_of != pa[h]) {
                    row_ok = series(pa[h], t, row, false);
                    row_of = pa[h];
                }
                if (row_ok) {
                    s0 = row.M[pb[h]];
                    s1 = row.M1[pb[h]];
                    s2 = row.M2[pb[h]];
                }
            }
            // a probability lost in rounding (codons two or three substitutions apart at a t far below the data's)
            if (!(s0 > 1e-300)) {
                s0 = 1e-300;
                s1 = s2 = 0;
            }
            f += pn[h] * std::log(s0);
            const double q = s1 / s0;
            d1 += pn[h] * q;
            d2 += pn[h] * (s2 / s0 - q * q);
        }
        if (g1) *g1 = d1;
        if (g2) *g2 = d2;
        return f;
    }

    // The gradient of lnl at t with respect to log kappa (g[0]) and log omega (g[1]), for the current
    // decomposition; slope_t: d lnl / dt there (0 at an interior maximum in t).
    // With B the symmetrised, unscaled rate matrix (B = V' D V, eigenvector k in row k of A), tau = t / mr and
    // M = exp(B tau): dM = sum_kl v_k v_l' G_kl F_kl, G = V dB V', F_kl = (e^(d_k tau) - e^(d_l tau)) / (d_k - d_l)
    // (tau e^(d_k tau) on the diagonal). So d lnl = sum_ab dB_ab Y_ab with Y = V' ((V w V') o F) V, where w has
    // n_p / M_p at each pattern's (a, b); and the scaling adds -t (d log mr) slope_t.
    void gradient(double t, double slope_t, double* grad) {
        const size_t n = (size_t)ns, np = pa.size();
        const double tau = t / mr;
        for (int k = 0; k < ns; ++k) E[k] = std::exp(mu[k] * t);
        std::fill(T.begin(), T.end(), 0.0);
        Row row;
        int row_of = -1;
        bool row_ok = false;
        double gk_tiny = 0, gw_tiny = 0;
        for (size_t h = 0; h < np; ++h) {
            const double* c = &C[h * n];
            double s0 = 0;
            for (int k = 0; k < ns; ++k) s0 += c[k] * E[k];
            const int a = pa[h], b = pb[h];
            if (!(s0 > TINY_PROB)) {
                // a probability the eigen sum cannot give (lnl): its derivative from the series as well
                if (row_of != a) {
                    row_ok = series(a, t, row, true);
                    row_of = a;
                }
                if (row_ok) {
                    if (row.M[b] > 1e-300) {
                        gk_tiny += pn[h] * row.Gk[b] / row.M[b];
                        gw_tiny += pn[h] * row.Gw[b] / row.M[b];
                    }
                    continue;
                }
            }
            if (!(s0 > 1e-300)) continue;
            const double w = pn[h] / s0;
            for (size_t k = 0; k < n; ++k) T[k * n + b] += w * A[k * n + a];
        }
        // W = T V' (as sums of rows of V', so that no sum waits on the one before), times F (symmetric)
        for (size_t l = 0; l < n; ++l)
            for (size_t b = 0; b < n; ++b) Vt[b * n + l] = A[l * n + b];
        std::fill(W.begin(), W.end(), 0.0);
        for (size_t k = 0; k < n; ++k) {
            double* wk = &W[k * n];
            const double* tk = &T[k * n];
            for (size_t b = 0; b < n; ++b) {
                const double x = tk[b];
                if (x == 0) continue;   // (columns b that no pattern has)
                const double* vb = &Vt[b * n];
                for (size_t l = 0; l < n; ++l) wk[l] += x * vb[l];
            }
        }
        for (size_t k = 0; k < n; ++k) {
            W[k * n + k] *= tau * E[k];
            for (size_t l = 0; l < k; ++l) {
                const double x = (d[k] - d[l]) * tau;
                const double F = tau * E[l] * (std::fabs(x) < 1e-12 ? 1.0 : std::expm1(x) / x);
                W[k * n + l] *= F;
                W[l * n + k] *= F;
            }
        }
        // U = W V (U[k][b] = sum_l W_kl V[l][b])
        std::fill(U.begin(), U.end(), 0.0);
        for (size_t k = 0; k < n; ++k)
            for (size_t l = 0; l < n; ++l) {
                const double x = W[k * n + l];
                const double* vl = &A[l * n];
                double* uk = &U[k * n];
                for (size_t b = 0; b < n; ++b) uk[b] += x * vl[b];
            }
        auto Y = [&](int a, int b) {
            double sum = 0;
            for (size_t k = 0; k < n; ++k) sum += A[k * n + a] * U[k * n + b];
            return sum;
        };
        double gk = 0, gw = 0;
        for (int a = 0; a < ns; ++a) {
            const double yaa = Y(a, a);
            for (int k = 0; k < nnb[a]; ++k) {
                if (!nb_ts[a][k] && !nb_non[a][k]) continue;
                const int b = nb[a][k];
                const double v = sq[a][k] * (std::sqrt(pi[idx[a]] * pi[idx[b]]) * Y(a, b) - pi[idx[b]] * yaa);
                if (nb_ts[a][k]) gk += v;
                if (nb_non[a][k]) gw += v;
            }
        }
        grad[0] = gk + gk_tiny - t * (mr_ts / mr) * slope_t;
        grad[1] = gw + gw_tiny - t * rhoN * slope_t;
    }
};

constexpr double T_MIN = 1e-6, T_MAX = 50;
constexpr double PAR_MIN = 1e-6, PAR_MAX = 999;   // kappa and omega
// the stopping rule tries the bound itself from this near (in log units: below 1e-3, above 99)
const double BOUND_ZONE_LOW = std::log(1e-3 / 1e-6), BOUND_ZONE_HIGH = std::log(999.0 / 99.0);
constexpr double SE_PAR_MIN = 1e-4, SE_PAR_MAX = 99;   // between these, kappa or omega counts as estimated (SE)

// max over t in [T_MIN, T_MAX] of lnl(t), by Newton steps in u = log t inside a bracket; t: start and result.
double profile_t(GyModel& m, double& t, double* slope = nullptr) {
    const double ulo = std::log(T_MIN), uhi = std::log(T_MAX);
    double u = std::log(std::min(std::max(t, T_MIN), T_MAX));
    double lo = ulo, hi = uhi;   // bracket of the maximum, narrowed by the sign of the slope
    double g1, g2;
    double f = m.lnl(std::exp(u), &g1, &g2);
    for (int it = 0; it < 200; ++it) {
        const double tt = std::exp(u);
        const double du1 = tt * g1, du2 = tt * tt * g2 + tt * g1;   // derivatives in u
        if (du1 > 0)
            lo = u;
        else
            hi = u;
        if ((u <= ulo && du1 <= 0) || (u >= uhi && du1 >= 0)) break;   // at a bound, slope pointing out
        double step = du2 < 0 ? -du1 / du2 : (du1 > 0 ? 1.0 : -1.0);
        if (step > 2) step = 2;
        if (step < -2) step = -2;
        double un = u + step;
        // outside the bracket: the bound itself where the bracket still ends at it, else the middle
        if (un >= hi) un = hi >= uhi ? uhi : 0.5 * (lo + hi);
        if (un <= lo) un = lo <= ulo ? ulo : 0.5 * (lo + hi);
        if (hi - lo < 1e-13) break;
        // ascent, or halve towards u (the slope's sign says the maximum of a unimodal profile lies that way)
        double fn = 0, n1 = 0, n2 = 0;
        int half = 0;
        for (;; ++half) {
            fn = m.lnl(std::exp(un), &n1, &n2);
            if (fn >= f - 1e-13 * std::fabs(f) || half >= 60) break;
            un = 0.5 * (un + u);
        }
        const double moved = std::fabs(un - u);
        u = un;
        f = fn;
        g1 = n1;
        g2 = n2;
        if (moved < 1e-13 || std::fabs(std::exp(u) * g1) < 1e-10) break;
    }
    t = std::exp(u);
    if (u <= ulo) t = T_MIN;
    if (u >= uhi) t = T_MAX;
    if (slope) *slope = g1;
    return f;
}

struct Fit {
    double lnl = -1e300, t = NA, kappa = NA, omega = NA;
    bool converged = false;
    bool at_limit = false;   // the iteration used all its steps
    bool crawled = false;    // it took a step that the second run of optimise() would lengthen
    // kappa and omega where the iteration first stopped by itself, before the ends of a bound's zone were tried
    // (the result itself when they were not): what later starts are seeded from, so that none begins on a bound
    double seed_kappa = NA, seed_omega = NA;
    int iterations = 0;
    // optimise(): where the first run stopped (its seed_kappa, seed_omega and t) when the second (lengthened) run
    // replaced it
    bool has_first_run_stop = false;
    double first_run_stop_kappa = NA, first_run_stop_omega = NA, first_run_stop_t = NA;
};

// -lnl at y = (log kappa, log omega), profiled over t (tt: start and result), and, with g, its gradient. The model
// is left decomposed at y (an evaluation changes the model's state, so it is not a const operation).
struct ProfileEval {
    GyModel& m;
    int nfree;   // 1: omega is fixed, and its gradient is reported as 0
    double operator()(const double* y, double& tt, double* g, double* dt = nullptr) {
        m.build(std::exp(y[0]), std::exp(y[1]));
        double slope = 0;
        const double f = -profile_t(m, tt, &slope);
        if (dt) *dt = slope;   // d lnl / dt where the profile stopped
        if (g) {
            double gl[2];
            m.gradient(tt, slope, gl);
            g[0] = -gl[0];
            g[1] = nfree > 1 ? -gl[1] : 0.0;
        }
        return f;
    }
};

// The stopping rule of optimise_run near a bound. A small gradient in log kappa and log omega means a maximum only
// inside the bounds. Where the data want kappa or omega at 0 the likelihood is linear in the parameter itself, so
// its gradient in the logarithm is the parameter times a constant: it falls below the tolerance at about 1e-5, a
// factor 10 short of the bound and up to 1e-5 in lnL below the value there, at a point that depends on the path.
// And it fades in the other direction as well: an iteration that is near or on the bound sees no gradient towards
// a maximum that lies further inside (a start on the bound "converged" there, 6.6e-3 below a maximum at kappa
// 0.0075: PR #97). So when the iteration stops with a parameter in the zone of a bound (below 1e-3, above 99), the
// two ends of the zone are tried: the bound itself, taken if the likelihood there is not lower, and the inner edge
// (1e-3 or 99), where the logarithmic gradient is alive again, taken if the likelihood there is higher. The
// iteration then goes on from the better one: true is returned and x, f, t, g (and H, reset at the inner edge) are
// those of that point. Twice per parameter at most (probes). A stop outside the zones is final (false).
// Before the first probe of a run, the point where the iteration stopped by itself goes to fit.seed_* (seeded).
bool probe_zone_ends(ProfileEval& eval, double x[2], double& f, double& t, double g[2], double H[2][2],
                     int probes[2], Fit& fit, bool& seeded) {
    const double lo = std::log(PAR_MIN), hi = std::log(PAR_MAX);
    const int nfree = eval.nfree;
    bool moved = false;
    for (int i = 0; i < nfree && !moved; ++i) {
        const bool low = x[i] < lo + BOUND_ZONE_LOW, high = x[i] > hi - BOUND_ZONE_HIGH;
        if (probes[i] >= 2 || !(low || high)) continue;
        if (probes[0] + probes[1] == 0) {   // where the iteration stopped by itself: the seed of later starts
            fit.seed_kappa = std::exp(x[0]);
            fit.seed_omega = std::exp(x[1]);
            seeded = true;
        }
        ++probes[i];
        const double bound = low ? lo : hi, edge = low ? lo + BOUND_ZONE_LOW : hi - BOUND_ZONE_HIGH;
        double fb = f, tb = t, gb[2] = {g[0], g[1]}, xb[2] = {x[0], x[1]};
        bool take = false;
        if (x[i] != bound) {
            double xs[2] = {x[0], x[1]}, gs[2], ts = t;
            xs[i] = bound;
            const double fs = eval(xs, ts, gs);
            if (fs <= fb) {
                take = true;
                fb = fs;
                tb = ts;
                xb[i] = bound;
                gb[0] = gs[0];
                gb[1] = gs[1];
            }
        }
        {
            double xs[2] = {x[0], x[1]}, gs[2], ts = t;
            xs[i] = edge;
            const double fs = eval(xs, ts, gs);
            if (fs < fb - 1e-9) {
                take = true;
                fb = fs;
                tb = ts;
                xb[i] = edge;
                gb[0] = gs[0];
                gb[1] = gs[1];
                H[0][0] = H[1][1] = 1;   // (what was learnt near the bound says nothing here)
                H[0][1] = H[1][0] = 0;
            }
        }
        if (take) {
            x[0] = xb[0];
            x[1] = xb[1];
            f = fb;
            t = tb;
            g[0] = gb[0];
            g[1] = gb[1];
            moved = true;
        }
    }
    return moved;
}

// A step of optimise_run that the second run of optimise() lengthens (see there and at the call): from x along
// dir (largest component dmax, directional derivative `slope`), accepted with factor alpha at xn (likelihood fn,
// branch length tn). The factor is doubled while the likelihood keeps improving, up to a step of 2 log units;
// alpha, xn, fn, tn and sn (the slope in t) are those of the last accepted trial, and the model is left decomposed
// at xn.
void lengthen_step(ProfileEval& eval, const double x[2], const double dir[2], double dmax, double slope,
                   double& alpha, double xn[2], double& fn, double& tn, double& sn) {
    const double lo = std::log(PAR_MIN), hi = std::log(PAR_MAX);
    double a2 = 2, xa[2];
    bool stale = false;   // the model's state is that of a rejected trial
    // (up to 2 log units, the limit of every later step, also for the first step: its limit of 0.5 is
    // there for the unscaled gradient of an ordinary pair, hundreds per log unit, and every doubling here
    // has to improve the likelihood to be kept)
    while (a2 * dmax <= 2.0) {
        for (int i = 0; i < 2; ++i) xa[i] = std::min(hi, std::max(lo, x[i] + a2 * dir[i]));
        double ta = tn, sa = 0;
        const double ft = eval(xa, ta, nullptr, &sa);
        if (!(ft < fn + 1e-4 * (a2 - alpha) * slope)) {
            stale = true;
            break;
        }
        stale = false;
        alpha = a2;
        fn = ft;
        tn = ta;
        sn = sa;
        xn[0] = xa[0];
        xn[1] = xa[1];
        a2 *= 2;
    }
    if (stale) fn = eval(xn, tn, nullptr, &sn);
}

// Maximise over (log kappa, log omega) (omega fixed when fix_omega) the likelihood profiled over t: a projected
// BFGS iteration with the analytic gradient (at the maximum in t the gradient of the profile is the partial
// gradient; at a bound of t the term in the slope accounts for the scaling of the rate matrix).
// lengthen: see optimise() below and the comment at the line search.
Fit optimise_run(GyModel& m, double kappa0, double omega0, double t0, bool fix_omega, bool lengthen) {
    const double lo = std::log(PAR_MIN), hi = std::log(PAR_MAX);
    const int nfree = fix_omega ? 1 : 2;
    double x[2] = {std::log(kappa0), std::log(omega0)};
    double t = t0;
    ProfileEval eval{m, nfree};
    double g[2];
    double f = eval(x, t, g);
    double H[2][2] = {{1, 0}, {0, 1}};
    bool act_prev[2] = {false, false};
    Fit fit;
    const int maxit = 100;
    int it = 0;
    bool restarted = false;
    int small_steps = 0;   // successive steps that changed the likelihood by less than 1e-9
    int at_tmax = 0;       // successive iterates with t at its upper bound
    bool curvature = false;   // a BFGS update has been made (H holds curvature of this surface)
    int probes[2] = {0, 0};   // times the ends of a bound's zone were tried for this parameter (stopping rule)
    bool seeded = false;      // fit.seed_* hold the point of the first stop
    for (; it < maxit; ++it) {
        bool act[2];
        double pg[2];
        double gmax = 0;
        for (int i = 0; i < 2; ++i) {
            act[i] = i >= nfree || (x[i] <= lo && g[i] > 0) || (x[i] >= hi && g[i] < 0);
            pg[i] = act[i] ? 0 : g[i];
            gmax = std::max(gmax, std::fabs(pg[i]));
        }
        // (a gradient of 1e-5 per unit of log kappa or log omega is, at a maximum inside the bounds and with the
        // curvatures of real pairs, a distance of 1e-6 or less from the maximum and 1e-11 in lnL; the rounding of
        // the likelihood sits just below)
        if (gmax < 1e-5 || (small_steps >= 2 && gmax < 1e-3)) {
            // (not so near a bound: probe_zone_ends)
            if (probe_zone_ends(eval, x, f, t, g, H, probes, fit, seeded)) {
                small_steps = 0;
                continue;
            }
            fit.converged = true;
            break;
        }
        // On the ridge of saturation (t at its bound) the likelihood keeps creeping up along omega * t = constant
        // without a maximum inside the bounds: the pair is saturated, and more steps change nothing in that
        // (stop when two successive steps there gained less than 1e-9; the caller reads t, not `converged`).
        at_tmax = t >= T_MAX * (1 - 1e-3) ? at_tmax + 1 : 0;
        if (at_tmax >= 4 && small_steps >= 2) break;
        if (act[0] != act_prev[0] || act[1] != act_prev[1]) {
            H[0][0] = H[1][1] = 1;
            H[0][1] = H[1][0] = 0;
        }
        act_prev[0] = act[0];
        act_prev[1] = act[1];
        double dir[2];
        for (int i = 0; i < 2; ++i) dir[i] = act[i] ? 0 : -(H[i][0] * pg[0] + H[i][1] * pg[1]);
        double slope = dir[0] * pg[0] + dir[1] * pg[1];
        if (!(slope < 0)) {
            H[0][0] = H[1][1] = 1;
            H[0][1] = H[1][0] = 0;
            dir[0] = -pg[0];
            dir[1] = -pg[1];
            slope = -(pg[0] * pg[0] + pg[1] * pg[1]);
        }
        const double dmax = std::max(std::fabs(dir[0]), std::fabs(dir[1]));
        const double cap = it == 0 ? 0.5 : 2.0;   // log units per step
        if (dmax > cap) {
            dir[0] *= cap / dmax;
            dir[1] *= cap / dmax;
            slope *= cap / dmax;
        }
        double alpha = 1, xn[2], fn = f, tn = t, sn = 0, gn[2];
        bool ok = false;
        for (int ls = 0; ls < 20; ++ls) {
            for (int i = 0; i < 2; ++i) xn[i] = std::min(hi, std::max(lo, x[i] + alpha * dir[i]));
            tn = t;
            fn = eval(xn, tn, nullptr, &sn);
            if (fn <= f + 1e-4 * alpha * slope) {
                ok = true;
                break;
            }
            alpha *= 0.5;
        }
        // (Only in the second run of optimise().) A step of the gradient's own size (no curvature seen since the
        // start: H is the first identity) on a nearly flat surface: a gradient of 1e-4 moves 1e-4 per step (in
        // log kappa, log omega), and with next to no curvature the BFGS update never gets a scale, so the
        // iteration would crawl to its limit (two random codons, every site different: 100 steps, 3e-7 gained,
        // `not_converged`). Lengthen such a step while the likelihood keeps improving, up to the step limit: the
        // iteration then reaches the curvature or the bounds, where it can converge. Steps with a larger
        // gradient, or after the first BFGS update, are not touched.
        const bool crawling = ok && alpha == 1 && dmax < 1e-2 && !curvature;
        if (crawling) fit.crawled = true;
        if (lengthen && crawling) lengthen_step(eval, x, dir, dmax, slope, alpha, xn, fn, tn, sn);
        if (!ok) {
            // no descent along this direction: once more from the plain gradient, then accept the point
            if (!restarted && (H[0][0] != 1 || H[1][1] != 1 || H[0][1] != 0)) {
                restarted = true;
                H[0][0] = H[1][1] = 1;
                H[0][1] = H[1][0] = 0;
                continue;
            }
            fit.converged = gmax < 1e-2;
            break;
        }
        restarted = false;
        {   // the gradient at the accepted point (the model's state and sn are those of the last evaluation: xn)
            double gl[2];
            // (sn is the slope where the profile stopped; at a bound tn is the bound itself, which the exponential
            // of its logarithm need not be to the last digit: there the slope is taken at tn)
            if (tn <= T_MIN || tn >= T_MAX) m.lnl(tn, &sn);
            m.gradient(tn, sn, gl);
            gn[0] = -gl[0];
            gn[1] = nfree > 1 ? -gl[1] : 0.0;
        }
        const double s[2] = {xn[0] - x[0], xn[1] - x[1]};
        const double y[2] = {gn[0] - g[0], gn[1] - g[1]};
        const double sy = s[0] * y[0] + s[1] * y[1];
        if (sy > 1e-12 * (std::fabs(s[0]) + std::fabs(s[1]))) {
            // the first update after a start from the identity: scale it to the curvature just seen, so that the
            // direction no step has explored yet does not get a step of the gradient's size
            if (H[0][0] == 1 && H[1][1] == 1 && H[0][1] == 0 && H[1][0] == 0) {
                const double yy = y[0] * y[0] + y[1] * y[1];
                if (yy > 0) H[0][0] = H[1][1] = sy / yy;
            }
            curvature = true;
            // BFGS update of the inverse Hessian
            const double Hy[2] = {H[0][0] * y[0] + H[0][1] * y[1], H[1][0] * y[0] + H[1][1] * y[1]};
            const double yHy = y[0] * Hy[0] + y[1] * Hy[1];
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j)
                    H[i][j] += (1 + yHy / sy) * s[i] * s[j] / sy - (Hy[i] * s[j] + s[i] * Hy[j]) / sy;
        }
        small_steps = f - fn < 1e-9 ? small_steps + 1 : 0;
        x[0] = xn[0];
        x[1] = xn[1];
        f = fn;
        t = tn;
        g[0] = gn[0];
        g[1] = gn[1];
    }
    fit.iterations = it;
    // the state of the model is that of the last evaluation, not of x: rebuild at x
    fit.lnl = -eval(x, t, nullptr);
    fit.t = t;
    fit.kappa = std::exp(x[0]);
    fit.omega = std::exp(x[1]);
    if (x[0] <= lo) fit.kappa = PAR_MIN;
    if (x[0] >= hi) fit.kappa = PAR_MAX;
    if (x[1] <= lo) fit.omega = PAR_MIN;
    if (x[1] >= hi) fit.omega = PAR_MAX;
    if (!seeded) {
        fit.seed_kappa = fit.kappa;
        fit.seed_omega = fit.omega;
    }
    fit.at_limit = it >= maxit;
    return fit;
}

// Two log-likelihoods within LNL_TIE of each other are the same maximum as far as the optimiser can tell (its
// stopping rule leaves about 1e-11, a start from elsewhere up to 1e-8): the earlier or the converged one is kept.
constexpr double LNL_TIE = 1e-7;
// An iteration that did not converge displaces the best usable start (run_gy) only if its likelihood is higher by
// more than this: below it the two are, for every use of the result, the same likelihood.
constexpr double LNL_MATERIAL = 1e-4;

// One start. An iteration that ends at its limit of steps without converging, away from the bound of t, and that
// took steps of the gradient's size without curvature (`crawled`: over a nearly flat surface; two random codons
// with every site different: 100 steps, 3e-7 gained) is run once more with such steps lengthened. That run is
// taken if it converges and its likelihood is not lower by more than LNL_TIE (1e-7); a second run that converges
// lower than that found another, lower maximum and is dropped. An iteration that never took such a step would be
// repeated step for step and is not run again. Every start that converges, or ends on the ridge of saturation,
// is exactly what it was without this second run.
Fit optimise(GyModel& m, double kappa0, double omega0, double t0, bool fix_omega) {
    Fit fit = optimise_run(m, kappa0, omega0, t0, fix_omega, false);
    if (!fit.converged && fit.at_limit && fit.crawled && fit.t < T_MAX * (1 - 1e-3)) {
        Fit again = optimise_run(m, kappa0, omega0, t0, fix_omega, true);
        const int its = fit.iterations + again.iterations;
        if (again.converged && again.lnl >= fit.lnl - LNL_TIE) {
            again.has_first_run_stop = true;
            again.first_run_stop_kappa = fit.seed_kappa;
            again.first_run_stop_omega = fit.seed_omega;
            again.first_run_stop_t = fit.t;
            fit = again;
        }
        fit.iterations = its;
    }
    return fit;
}

// SE of dN and dS at the maximum `fit` into r.se_ka / r.se_ks: the inverse of the observed information (minus
// the Hessian of the log-likelihood, by central differences of the exact gradient) of the parameters that are
// estimated inside their bounds, and the delta method. use[0..2]: t, log kappa, log omega take part; a parameter
// at a bound is held fixed there (the caller then gives no SE to the distance that the bound makes 0). Left NaN
// when the information is not positive definite (a flat maximum: no standard error to give).
void standard_errors(GyModel& m, const Fit& fit, const bool use[3], Result& r, bool want_s, bool want_n) {
    const double x0[3] = {fit.t, std::log(fit.kappa), std::log(fit.omega)};
    const double h[3] = {1e-4 * fit.t, 1e-4, 1e-4};
    int id[3], n = 0;
    for (int i = 0; i < 3; ++i)
        if (use[i]) id[n++] = i;
    auto grad3 = [&](const double* x, double* g) {
        m.build(std::exp(x[1]), std::exp(x[2]));
        double slope = 0, gl[2];
        m.lnl(x[0], &slope);
        m.gradient(x[0], slope, gl);
        g[0] = slope;
        g[1] = gl[0];
        g[2] = gl[1];
    };
    double I[3][3];   // minus the Hessian, of the parameters in id[]
    for (int b = 0; b < n; ++b) {
        const int j = id[b];
        double xp[3] = {x0[0], x0[1], x0[2]}, xm[3] = {x0[0], x0[1], x0[2]}, gp[3], gm[3];
        xp[j] += h[j];
        xm[j] -= h[j];
        grad3(xp, gp);
        grad3(xm, gm);
        for (int a = 0; a < n; ++a) I[a][b] = -(gp[id[a]] - gm[id[a]]) / (2 * h[j]);
    }
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < a; ++b) I[a][b] = I[b][a] = 0.5 * (I[a][b] + I[b][a]);
    // Cholesky I = L L' (fails unless positive definite), then V = I^-1 through L^-1
    double L[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int a = 0; a < n; ++a)
        for (int b = 0; b <= a; ++b) {
            double sum = I[a][b];
            for (int k = 0; k < b; ++k) sum -= L[a][k] * L[b][k];
            if (a == b) {
                if (!(I[a][a] > 0) || !(sum > 1e-10 * I[a][a])) return;
                L[a][a] = std::sqrt(sum);
            } else {
                L[a][b] = sum / L[b][b];
            }
        }
    double Li[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int a = 0; a < n; ++a) {
        Li[a][a] = 1 / L[a][a];
        for (int b = 0; b < a; ++b) {
            double sum = 0;
            for (int k = b; k < a; ++k) sum -= L[a][k] * Li[k][b];
            Li[a][b] = sum / L[a][a];
        }
    }
    double V[3][3];
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b) {
            double sum = 0;
            for (int k = 0; k < n; ++k) sum += Li[k][a] * Li[k][b];
            V[a][b] = sum;
        }
    // dS = t (syn + non0) / (3 (syn + omega non0)), dN = omega dS: the Jacobian in the same parameters
    auto dsdn = [&](const double* x, double& dS, double& dN) {
        double syn, non0;
        m.rates(std::exp(x[1]), syn, non0);
        const double w = std::exp(x[2]);
        dS = x[0] * (syn + non0) / (3 * (syn + w * non0));
        dN = w * dS;
    };
    double JS[3], JN[3];
    for (int b = 0; b < n; ++b) {
        const int j = id[b];
        double xp[3] = {x0[0], x0[1], x0[2]}, xm[3] = {x0[0], x0[1], x0[2]}, sp, np, sm, nm;
        xp[j] += h[j];
        xm[j] -= h[j];
        dsdn(xp, sp, np);
        dsdn(xm, sm, nm);
        JS[b] = (sp - sm) / (2 * h[j]);
        JN[b] = (np - nm) / (2 * h[j]);
    }
    double vS = 0, vN = 0;
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b) {
            vS += JS[a] * V[a][b] * JS[b];
            vN += JN[a] * V[a][b] * JN[b];
        }
    if (want_s && vS > 0 && std::isfinite(vS)) r.se_ks = std::sqrt(vS);
    if (want_n && vN > 0 && std::isfinite(vN)) r.se_ka = std::sqrt(vN);
}

std::string fmt6(double v) {
    char buf[40];
    std::snprintf(buf, sizeof buf, "%.6g", v);
    return buf;
}

// omega of a start on the ridge of saturation (t at T_MAX) from a point (omega, t): omega * t, about dN, kept
// (within [PAR_MIN, 1])
double ridge_omega(double omega, double t) {
    return std::min(std::max(omega * t / T_MAX, PAR_MIN), 1.0);
}

// Which of the starts of run_gy is reported (the indices of Choice are into `fits`, the starts in the order in
// which they were run; -1: none).
//
// A start is usable if it converged or ended on the ridge (t at its bound: saturated whatever else). Among
// the usable starts the highest likelihood wins, the earlier or the converged one within LNL_TIE (1e-7). A
// start that is not usable (it ran out of steps somewhere else) never displaces a usable one by itself: if
// its likelihood is higher than that of the best usable start (which may itself be an unconverged ridge
// start), the iteration is continued from where it stopped (go_on; run_gy does that twice at most, and not
// again once a continuation ends lower), and only if it still has not converged and is higher than the best
// usable start by more than LNL_MATERIAL (1e-4) is the pair `not_converged`: the usable start is then known not
// to be the highest point. (No pair of 25,200 random ones gets there: tests/cpp/kaks_gy_choose_start.cpp does,
// with made-up starts.) Without a usable start the best point is reported, as `not_converged`, after the same
// continuations.
struct Choice {
    int best = -1;     // the best usable start
    int cand = -1;     // the best start that is not usable
    int report = -1;   // the one of the two whose estimate is the result
    bool go_on = false;   // continue the iteration of `cand` before reporting
};

bool usable_start(const Fit& f) {
    return f.converged || f.t >= T_MAX * (1 - 1e-3);
}

Choice choose_start(const std::vector<Fit>& fits) {
    Choice c;
    for (int i = 0; i < (int)fits.size(); ++i) {
        const Fit& f = fits[i];
        if (usable_start(f)) {
            if (c.best < 0 || f.lnl > fits[c.best].lnl + LNL_TIE ||
                (!fits[c.best].converged && f.converged && f.lnl >= fits[c.best].lnl - LNL_TIE))
                c.best = i;
        } else if (c.cand < 0 || f.lnl > fits[c.cand].lnl) {
            c.cand = i;
        }
    }
    c.go_on = c.cand >= 0 && (c.best < 0 || fits[c.cand].lnl > fits[c.best].lnl + LNL_TIE);
    const bool displaced = c.cand >= 0 && (c.best < 0 || fits[c.cand].lnl > fits[c.best].lnl + LNL_MATERIAL);
    c.report = displaced ? c.cand : c.best;
    return c;
}

// The coarse scan of run_gy (see there): the likelihood, maximised over t from t0, at each point of a grid in kappa and
// omega (kappa alone when omega is fixed), highest first (ties in the order of the grid). The grid avoids the
// bounds (an iteration started on a bound sees no gradient in the logarithm and stalls) and spans their zones'
// inner edges, 1e-3 and 99.
constexpr double SCAN_VALUES[] = {1e-3, 1e-2, 0.1, 1, 10, 99};
constexpr size_t SCAN_SMALL_PAIR = 100;   // codons compared: every such pair is scanned
constexpr int SCAN_STARTS = 2;            // starts from the highest points of the scan

struct ScanPoint {
    double lnl, kappa, omega, t;
};

std::vector<ScanPoint> scan_grid(GyModel& m, double t0, bool fix_omega) {
    std::vector<ScanPoint> pts;
    for (const double k : SCAN_VALUES)
        for (const double w : SCAN_VALUES) {
            const double omega = fix_omega ? 1.0 : w;
            m.build(k, omega);
            double t = t0;
            const double l = profile_t(m, t);
            pts.push_back({l, k, omega, t});
            if (fix_omega) break;
        }
    std::stable_sort(pts.begin(), pts.end(), [](const ScanPoint& a, const ScanPoint& b) { return a.lnl > b.lnl; });
    return pts;
}

// The continuation that choose_start asks for (go_on): the iteration of the best start that is not usable is run
// again from where it stopped (cont(from): the start's result), twice at most. Before each continuation the
// starts that are not usable are dropped (only the best of them is continued, and the continuation takes its
// place). A continuation that ends lower than the point it started from, without becoming usable, is replaced by
// that point, and there is no further one (it would repeat the same run). Returns the final choice; the
// continuations' iterations are added to `iterations`.
template <class Cont>
Choice continue_starts(std::vector<Fit>& fits, Choice ch, Cont&& cont, int& iterations) {
    for (int c = 0; c < 2 && ch.go_on; ++c) {
        const Fit from = fits[(size_t)ch.cand];
        fits.erase(std::remove_if(fits.begin(), fits.end(), [](const Fit& f) { return !usable_start(f); }),
                   fits.end());
        const Fit next = cont(from);
        iterations += next.iterations;
        const bool lower = !usable_start(next) && next.lnl < from.lnl;
        fits.push_back(lower ? from : next);
        ch = choose_start(fits);
        if (lower) break;
    }
    return ch;
}

// Saturation: an estimate above GY_MAX_DISTANCE substitutions per site is NA (t at its bound of 50 per codon
// always has one: the mean of dS and dN over the sites is t / 3). The other estimate is kept, as for the other
// methods, unless the optimiser did not settle on the ridge or the likelihood is flat in t, where nothing is an
// estimate (the caller then sets both satS and satN, where the pair has such sites). Sets the status, the notes
// and the NAs of a pair with satS or satN; dS and dN are the estimates.
void apply_saturation(Result& r, bool satS, bool satN, double dS, double dN, bool flat, bool t_at_max) {
    if (r.status != Status::no_data && r.status != Status::error) r.status = Status::saturated;
    char val[32];
    auto why = [&](double v) {
        if (flat)
            return std::string("the likelihood is flat in the branch length up to its upper bound, t = 50: the "
                               "data put no upper limit on the distance");
        if (!(v > GY_MAX_DISTANCE)) return std::string("the branch length is at its upper bound, t = 50");
        std::snprintf(val, sizeof val, "%.3g", v);
        return std::string("the maximum-likelihood estimate, ") + val + ", is above 10 substitutions per site" +
               (t_at_max ? "; the branch length is at its upper bound, t = 50" : "");
    };
    if (satS) {
        add_note(r, "Ks saturated (" + why(dS) + ")");
        r.ks = NA;
    }
    if (satN) {
        add_note(r, "Ka saturated (" + why(dN) + ")");
        r.ka = NA;
    }
    r.omega = NA;
}

}  // namespace

namespace detail {

// Test hooks (tests/cpp/kaks_gy_choose_start.cpp, which declares them; no header has them). They are compiled into
// every build of the library, not only with COAL_TEST_HOOKS: that option is for coal-cpp's environment-driven
// hooks and applies to the coal-cpp target alone, and no CI job builds a hooked library, so behind it the test
// would run nowhere. They are pure functions of their arguments, read nothing from the environment and are not
// reached from the engine.
//
// gy_choose_start: choose_start for n starts given by their likelihood, branch length and convergence. out: best,
// cand, report, go_on.
void gy_choose_start(const double* lnl, const double* t, const int* converged, int n, int out[4]);
void gy_choose_start(const double* lnl, const double* t, const int* converged, int n, int out[4]) {
    std::vector<Fit> fits((size_t)n);
    for (int i = 0; i < n; ++i) {
        fits[(size_t)i].lnl = lnl[i];
        fits[(size_t)i].t = t[i];
        fits[(size_t)i].converged = converged[i] != 0;
    }
    const Choice c = choose_start(fits);
    out[0] = c.best;
    out[1] = c.cand;
    out[2] = c.report;
    out[3] = c.go_on ? 1 : 0;
}

// gy_continue_starts: choose_start and continue_starts for n starts as above, the continuations scripted: the
// k-th continuation returns start k of the script (lnl, t, converged), whatever it starts from (past the end of
// the script: an unconverged start at -1e300). out: continuations made, the start reported (an index into the
// final list of starts), the number of starts in that list, go_on of the final choice; lnl_out: the likelihood
// of the start reported (-1e300: none).
void gy_continue_starts(const double* lnl, const double* t, const int* converged, int n, const double* s_lnl,
                        const double* s_t, const int* s_converged, int ns, int out[4], double* lnl_out);
void gy_continue_starts(const double* lnl, const double* t, const int* converged, int n, const double* s_lnl,
                        const double* s_t, const int* s_converged, int ns, int out[4], double* lnl_out) {
    std::vector<Fit> fits((size_t)n);
    for (int i = 0; i < n; ++i) {
        fits[(size_t)i].lnl = lnl[i];
        fits[(size_t)i].t = t[i];
        fits[(size_t)i].converged = converged[i] != 0;
    }
    int calls = 0, iterations = 0;
    const Choice c = continue_starts(
        fits, choose_start(fits),
        [&](const Fit&) {
            Fit f;
            if (calls < ns) {
                f.lnl = s_lnl[calls];
                f.t = s_t[calls];
                f.converged = s_converged[calls] != 0;
            }
            f.iterations = 1;
            ++calls;
            return f;
        },
        iterations);
    out[0] = calls;
    out[1] = c.report;
    out[2] = (int)fits.size();
    out[3] = c.go_on ? 1 : 0;
    *lnl_out = c.report >= 0 ? fits[(size_t)c.report].lnl : -1e300;
}

Result run_gy(const Prepared& p, const GeneticCode& g, const Options&, Result r) {
    const size_t n = p.c1.size();
    long snp = 0;
    for (size_t h = 0; h < n; ++h) {
        int k = 0;
        for (int pos = 0; pos < 3; ++pos) k += base_at(p.c1[h], pos) != base_at(p.c2[h], pos);
        snp += k;
    }
    r.substitutions = (double)snp;
    GyModel m(g);
    m.setup(p);
    const double L3 = 3.0 * (double)n;
    double syn1, non1;
    m.rates(1.0, syn1, non1);
    const bool hasS = syn1 > 0, hasN = non1 > 0;
    if (!hasS) {
        if (r.status != Status::error) r.status = Status::no_data;
        add_note(r, "no synonymous sites");
    }
    if (!hasN) {
        if (r.status != Status::error) r.status = Status::no_data;
        add_note(r, "no nonsynonymous sites");
    }
    if (!hasS && !hasN) return r;   // one codon, all its neighbours stop codons or of frequency 0
    if (snp == 0) {
        // Identical sequences: the likelihood is largest at t = 0, where it is sum n_i log pi_i whatever kappa
        // and omega are (lconst: with a = b, log sqrt(pi_a pi_b) is log pi_a; it was printed doubled before
        // 0.6.1). dN = dS = 0; kappa, omega and the site counts (which depend on kappa) have no estimate.
        r.lnl = m.lconst;
        r.t = 0;
        r.ks = hasS ? 0.0 : NA;
        r.ka = hasN ? 0.0 : NA;
        r.Sd = r.Nd = 0;
        r.divergence = 0;
        r.iterations = 0;
        return r;
    }
    // start of t: Jukes-Cantor distance per codon from the proportion of differing nucleotides
    const double pd = (double)snp / L3;
    double t0 = pd < 0.7 ? -3 * 0.75 * std::log(1 - 4 * pd / 3) : 3.0;
    t0 = std::min(std::max(t0, 1e-4), 5.0);
    const bool fix_omega = !hasS || !hasN;   // one kind of change only: omega cancels in the scaled matrix
    // Starting points. The likelihood of a diverged pair can have two maxima: one at a moderate t and one on the
    // ridge of saturation (t at or near its bound with a small omega); which one is higher depends on the pair,
    // and which one an iteration reaches depends on its start. Measured on 4,388 S. cerevisiae / C. glabrata
    // pairs: starts disagree (by more than 1e-4 in lnL) for 22 pairs, all with t above 2 at both maxima. So a pair
    // whose first maximum has t below 1 substitution per codon is done; any other pair is also started from
    // another (kappa, omega), from a small omega at t = 10, and from the ridge itself.
    //
    // Which start is reported: choose_start.
    std::vector<Fit> fits;
    int iterations = 0;
    const auto offer = [&](const Fit& f) {
        iterations += f.iterations;
        fits.push_back(f);
    };
    const Fit first = optimise(m, 2.0, fix_omega ? 1.0 : 0.3, t0, fix_omega);
    offer(first);
    if (!first.converged || first.t >= 1.0) {
        // the ridge start: t at its bound, omega such that omega * t (about dN) is that of the first maximum
        // (from where the first iteration stopped by itself, not from a bound it was then moved to: at a bound
        // the gradient in the logarithm is 0 in both directions and a start would stall)
        const double w_ridge = ridge_omega(first.seed_omega, first.t);
        const double start[3][3] = {{6.0, 1.5, t0}, {2.0, 0.02, 10.0}, {first.seed_kappa, w_ridge, T_MAX}};
        for (const double* st : start) offer(optimise(m, st[0], fix_omega ? 1.0 : st[1], st[2], fix_omega));
        // a second ridge start from where the first run of the first start stopped, when its second run replaced
        // it (PR #101): from the second run's stop the ridge start can end lower on the ridge. Not made where it
        // would be the ridge start above once more.
        if (first.has_first_run_stop) {
            const double w2 = ridge_omega(first.first_run_stop_omega, first.first_run_stop_t);
            if (first.first_run_stop_kappa != first.seed_kappa || w2 != w_ridge)
                offer(optimise(m, first.first_run_stop_kappa, fix_omega ? 1.0 : w2, T_MAX, fix_omega));
        }
    }
    // The coarse scan. On a pair of few codons the likelihood is flat and has several maxima, often with kappa or
    // omega at a bound, which none of the starts above may reach; and a start that ends on the ridge of saturation
    // stops on a plateau (at t = 50 the likelihood barely depends on kappa and omega), where it cannot see a higher
    // maximum at a moderate t or with another kappa (GTGGACAAC / TCCGGGAAG under code 11: every start ended on
    // the ridge, at lnL -17.3590 or below; the maximum is -16.5761 at t 7.96, kappa 1e-6). So a pair of
    // SCAN_SMALL_PAIR (100) codons or fewer, and a pair whose best start ended on the ridge (or that has no usable
    // start), has its likelihood evaluated on a grid of 36 points in (kappa, omega), each maximised over t
    // (scan_grid), and is started once more from the two highest points: a small pair always, any other pair only
    // from a point higher than its best start (there the extra starts would only repeat what the starts found).
    // The new starts are offered to choose_start like the others, so the likelihood reported can only be higher
    // (by more than LNL_TIE; within it a converged start can still replace an unconverged one on the ridge).
    {
        const Choice c0 = choose_start(fits);
        const bool small = n <= SCAN_SMALL_PAIR;
        const bool ridge = c0.best < 0 || fits[(size_t)c0.best].t >= T_MAX * (1 - 1e-3);
        if (small || ridge) {
            const double bar = c0.best >= 0 ? fits[(size_t)c0.best].lnl : -1e300;
            const std::vector<ScanPoint> pts = scan_grid(m, t0, fix_omega);
            for (int i = 0; i < SCAN_STARTS && i < (int)pts.size(); ++i) {
                if (!small && !(pts[(size_t)i].lnl > bar + LNL_TIE)) break;
                offer(optimise(m, pts[(size_t)i].kappa, pts[(size_t)i].omega, pts[(size_t)i].t, fix_omega));
            }
        }
    }
    const Choice ch = continue_starts(
        fits, choose_start(fits),
        [&](const Fit& from) { return optimise(m, from.kappa, fix_omega ? 1.0 : from.omega, from.t, fix_omega); },
        iterations);
    Fit best = fits[ch.report];
    best.iterations = iterations;
    r.iterations = best.iterations;
    m.build(best.kappa, best.omega);
    // A likelihood that is flat in t from the estimate to the bound (within 1e-6: one codon, one kind of change,
    // every site different) has its maximum at the bound as much as anywhere: the optimiser stops where the slope
    // falls below its tolerance, which is not an estimate of t. Such a pair is saturated, with t at the bound.
    bool flat = false;
    if (best.converged && best.t < T_MAX * (1 - 1e-3)) {
        const double at_bound = m.lnl(T_MAX);
        if (at_bound >= best.lnl - 1e-6) {
            flat = true;
            best.t = T_MAX;
            best.lnl = at_bound;
        }
    }
    r.kappa = best.kappa;
    r.rates[0] = r.rates[1] = best.kappa;
    r.S = L3 * m.rhoS1;
    r.N = L3 * m.rhoN1;
    // The branch length at (or within 0.1% of) its upper bound: saturated, whether or not the optimiser settled
    // on the ridge there (that it often does not is a property of the ridge, not a different result).
    const bool t_at_max = best.t >= T_MAX * (1 - 1e-3);
    if (!best.converged && !t_at_max) {
        if (r.status == Status::ok || r.status == Status::warning) r.status = Status::not_converged;
        r.converged = false;
        add_note(r, "GY: the likelihood optimisation did not converge: no estimate (last point: t " +
                        fmt6(best.t) + ", kappa " + fmt6(best.kappa) + ", omega " + fmt6(best.omega) +
                        "; such values can differ between platforms)");
        return r;
    }
    r.converged = true;
    r.lnl = best.lnl;
    r.t = best.t;
    const double dS = hasS ? best.t * m.rhoS / (3 * m.rhoS1) : NA;
    const double dN = hasN ? best.t * m.rhoN / (3 * m.rhoN1) : NA;
    r.ks = dS;
    r.ka = dN;
    // (dN / dS is computed from dN and dS as for every method, and equals omega; omega at a bound, 1e-6 or 999,
    // says that the data support no nonsynonymous, or no synonymous, change: docs/kaks.md)
    if (!fix_omega) r.omega = best.omega;
    const bool unsettled = !best.converged || flat;   // (then t is at its bound: see above)
    const bool satS = hasS && (unsettled || dS > GY_MAX_DISTANCE);
    const bool satN = hasN && (unsettled || dN > GY_MAX_DISTANCE);
    if (satS || satN) apply_saturation(r, satS, satN, dS, dN, flat, t_at_max);
    if (std::isfinite(r.ks) && std::isfinite(r.ka)) {
        r.divergence = (r.S * r.ks + r.N * r.ka) / (r.S + r.N);   // = t / 3
        // differences as KaKs_Calculator's GY gives them: the expected numbers S dS and N dN, scaled to the
        // observed number of nucleotide differences
        const double es = r.S * r.ks, en = r.N * r.ka;
        if (es + en > 0) {
            r.Sd = es * (double)snp / (es + en);
            r.Nd = en * (double)snp / (es + en);
        }
    } else if (!satS && !satN) {
        // one kind of site only: every difference is of that kind
        if (hasS) r.Sd = (double)snp;
        if (hasN) r.Nd = (double)snp;
    }
    // Standard errors of dN and dS, as codeml's getSE = 1 gives them (there with a Hessian by differences of the
    // likelihood): from the observed information of the parameters estimated inside their bounds. A kappa or an
    // omega near a bound (below 1e-4 or above 99; the bounds are 1e-6 and 999) is the bound in all but name: the
    // likelihood has no curvature to speak of in it, and it is held fixed; with omega there, the distance that the
    // bound makes 0 (dN at the lower bound, dS at the upper one) has no standard error. None either for a
    // saturated pair, or with one kind of site only.
    if (!fix_omega && !satS && !satN && std::isfinite(r.ks) && std::isfinite(r.ka) && best.t > T_MIN) {
        const bool k_in = best.kappa > SE_PAR_MIN && best.kappa < SE_PAR_MAX;
        const bool w_low = best.omega <= SE_PAR_MIN, w_high = best.omega >= SE_PAR_MAX;
        const bool use[3] = {true, k_in, !w_low && !w_high};
        standard_errors(m, best, use, r, !w_high, !w_low);
    }
    return r;
}

}  // namespace detail
}  // namespace kaks
}  // namespace coal
