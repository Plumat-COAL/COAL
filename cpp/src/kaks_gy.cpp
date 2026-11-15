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
// from fixed starting points (one for a close pair, four for a diverged one). Nothing is random and nothing depends
// on the thread: the same pair gives the same numbers in any run. docs/kaks.md, "GY: maximum likelihood", has the statuses, the bounds, the validation against
// codeml and the differences from codeml and KaKs_Calculator.
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
        for (int a = 0; a < ns; ++a) {
            double row = 0;
            for (int k = 0; k < nnb[a]; ++k) {
                const int b = nb[a][k];
                const double s = (nb_ts[a][k] ? kappa : 1.0) * (nb_non[a][k] ? omega : 1.0);
                sq[a][k] = s;
                if (nb_ts[a][k]) mr_ts += pi[idx[a]] * s * pi[idx[b]];
                row += s * pi[idx[b]];
                if (b < a) A[(size_t)a * ns + b] = A[(size_t)b * ns + a] = s * std::sqrt(pi[idx[a]] * pi[idx[b]]);
            }
            A[(size_t)a * ns + a] = -row;
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

    // log-likelihood at t with the current decomposition; g1, g2: first and second derivative in t
    double lnl(double t, double* g1 = nullptr, double* g2 = nullptr) {
        for (int k = 0; k < ns; ++k) E[k] = std::exp(mu[k] * t);
        double f = lconst, d1 = 0, d2 = 0;
        const size_t np = pa.size();
        for (size_t h = 0; h < np; ++h) {
            const double* c = &C[h * ns];
            double s0 = 0, s1 = 0, s2 = 0;
            for (int k = 0; k < ns; ++k) {
                const double x = c[k] * E[k];
                s0 += x;
                s1 += x * mu[k];
                s2 += x * mu[k] * mu[k];
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
        for (size_t h = 0; h < np; ++h) {
            const double* c = &C[h * n];
            double s0 = 0;
            for (int k = 0; k < ns; ++k) s0 += c[k] * E[k];
            if (!(s0 > 1e-300)) continue;
            const double w = pn[h] / s0;
            const int a = pa[h], b = pb[h];
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
        grad[0] = gk - t * (mr_ts / mr) * slope_t;
        grad[1] = gw - t * rhoN * slope_t;
    }
};

constexpr double T_MIN = 1e-6, T_MAX = 50;
constexpr double PAR_MIN = 1e-6, PAR_MAX = 999;   // kappa and omega
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
    int iterations = 0;
};

// Maximise over (log kappa, log omega) (omega fixed when fix_omega) the likelihood profiled over t: a projected
// BFGS iteration with the analytic gradient (at the maximum in t the gradient of the profile is the partial
// gradient; at a bound of t the term in the slope accounts for the scaling of the rate matrix).
Fit optimise(GyModel& m, double kappa0, double omega0, double t0, bool fix_omega) {
    const double lo = std::log(PAR_MIN), hi = std::log(PAR_MAX);
    const int nfree = fix_omega ? 1 : 2;
    double x[2] = {std::log(kappa0), std::log(omega0)};
    double t = t0;
    // -lnl at y, profiled over t (tt: start and result), and its gradient
    auto eval = [&](const double* y, double& tt, double* g) -> double {
        m.build(std::exp(y[0]), std::exp(y[1]));
        double slope = 0;
        const double f = -profile_t(m, tt, &slope);
        if (g) {
            double gl[2];
            m.gradient(tt, slope, gl);
            g[0] = -gl[0];
            g[1] = nfree > 1 ? -gl[1] : 0.0;
        }
        return f;
    };
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
    for (; it < maxit; ++it) {
        bool act[2];
        double pg[2];
        double gmax = 0;
        for (int i = 0; i < 2; ++i) {
            act[i] = i >= nfree || (x[i] <= lo && g[i] > 0) || (x[i] >= hi && g[i] < 0);
            pg[i] = act[i] ? 0 : g[i];
            gmax = std::max(gmax, std::fabs(pg[i]));
        }
        // (a gradient of 1e-5 per unit of log kappa or log omega is, with the curvatures of real pairs, a distance
        // of 1e-6 or less from the maximum and 1e-11 in lnL; the rounding of the likelihood sits just below)
        if (gmax < 1e-5 || (small_steps >= 2 && gmax < 1e-3)) {
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
        double alpha = 1, xn[2], fn = f, tn = t, gn[2];
        bool ok = false;
        for (int ls = 0; ls < 20; ++ls) {
            for (int i = 0; i < 2; ++i) xn[i] = std::min(hi, std::max(lo, x[i] + alpha * dir[i]));
            tn = t;
            fn = eval(xn, tn, nullptr);
            if (fn <= f + 1e-4 * alpha * slope) {
                ok = true;
                break;
            }
            alpha *= 0.5;
        }
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
        {   // the gradient at the accepted point (the model's state is that of the last evaluation: xn)
            double gl[2];
            double sl = 0;
            m.lnl(tn, &sl);
            m.gradient(tn, sl, gl);
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

}  // namespace

namespace detail {

Result run_gy(const Prepared& p, const GeneticCode& g, const Options& opt, Result r) {
    (void)opt;
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
        // and omega are. dN = dS = 0; kappa, omega and the site counts (which depend on kappa) have no estimate.
        r.lnl = 2 * m.lconst;
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
    // another (kappa, omega), from a small omega at t = 10, and from the ridge itself, and the highest maximum
    // wins (the earlier one on a tie within the noise of the optimisation, 1e-7).
    Fit best = optimise(m, 2.0, fix_omega ? 1.0 : 0.3, t0, fix_omega);
    if (!best.converged || best.t >= 1.0) {
        // the ridge start: t at its bound, omega such that omega * t (about dN) is that of the first maximum
        const double w_ridge = std::min(std::max(best.omega * best.t / T_MAX, PAR_MIN), 1.0);
        const double start[3][3] = {{6.0, 1.5, t0}, {2.0, 0.02, 10.0}, {best.kappa, w_ridge, T_MAX}};
        for (const double* st : start) {
            Fit next = optimise(m, st[0], fix_omega ? 1.0 : st[1], st[2], fix_omega);
            const int its = best.iterations + next.iterations;
            if (next.lnl > best.lnl + 1e-7 || (!best.converged && next.converged && next.lnl >= best.lnl - 1e-7))
                best = next;
            best.iterations = its;
        }
    }
    r.iterations = best.iterations;
    m.build(best.kappa, best.omega);
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
    // Saturation: an estimate above GY_MAX_DISTANCE substitutions per site is NA (t at its bound of 50 per codon
    // always has one: the mean of dS and dN over the sites is t / 3). The other estimate is kept, as for the other
    // methods, unless the optimiser did not settle on the ridge, where nothing is an estimate.
    const bool unsettled = !best.converged;   // (then t is at its bound: see above)
    const bool satS = hasS && (unsettled || dS > GY_MAX_DISTANCE);
    const bool satN = hasN && (unsettled || dN > GY_MAX_DISTANCE);
    if (satS || satN) {
        if (r.status != Status::no_data && r.status != Status::error) r.status = Status::saturated;
        char val[32];
        auto why = [&](double v) {
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
    // likelihood has no curvature to speak of in it, and it is held fixed; with omega there, the distance that the bound makes 0 (dN at the lower bound, dS at the upper one)
    // has no standard error. None either for a saturated pair, or with one kind of site only.
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
