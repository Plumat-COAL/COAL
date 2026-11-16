// GY: which of a pair's starts is reported (choose_start in cpp/src/kaks_gy.cpp), with made-up starts, and the
// continuation of a start that is not usable (continue_starts), with scripted continuations. The rule that an
// unconverged start displaces a usable one only when it is higher by more than 1e-4 in lnL, and the continuation
// itself, are reached by no pair of 25,200 random ones, so no test of the command can hold them in place; this
// one does, through the hooks the library has for them. Built and run by tests/test_kaks_gy.py against the build's libcoal_core. Exit 0 when
// every check holds; otherwise the failed checks on stderr and exit 1.
#include <cstddef>
#include <cstdio>
#include <vector>

namespace coal {
namespace kaks {
namespace detail {
// n starts by likelihood, branch length and convergence; out: the best usable start, the best start that is not
// usable, the one that is reported (indices; -1: none), and whether the iteration of the second is to be continued
void gy_choose_start(const double* lnl, const double* t, const int* converged, int n, int out[4]);
// the same starts, then the continuations, the k-th of which returns start k of the script; out: continuations
// made, the start reported (index into the final list), the size of that list, go_on at the end; lnl_out: the
// likelihood reported
void gy_continue_starts(const double* lnl, const double* t, const int* converged, int n, const double* s_lnl,
                        const double* s_t, const int* s_converged, int ns, int out[4], double* lnl_out);
// whether the coarse scan is made for a pair of n codons whose result has likelihood lnl, branch length t and
// convergence `converged`, with the plateau of saturation at `plateau`
bool gy_scan_wanted(size_t n, double lnl, double t, int converged, double plateau);
}  // namespace detail
}  // namespace kaks
}  // namespace coal

namespace {

struct Start {
    double lnl, t;
    int converged;
};

int bad = 0;

void expect(const char* what, const std::vector<Start>& s, int best, int cand, int report, int go_on) {
    std::vector<double> lnl, t;
    std::vector<int> conv;
    for (const Start& x : s) {
        lnl.push_back(x.lnl);
        t.push_back(x.t);
        conv.push_back(x.converged);
    }
    int out[4] = {-9, -9, -9, -9};
    coal::kaks::detail::gy_choose_start(lnl.data(), t.data(), conv.data(), (int)s.size(), out);
    if (out[0] != best || out[1] != cand || out[2] != report || out[3] != go_on) {
        std::fprintf(stderr, "FAILED: %s: best %d cand %d report %d go_on %d, expected %d %d %d %d\n", what, out[0],
                     out[1], out[2], out[3], best, cand, report, go_on);
        ++bad;
    }
}

void expect_cont(const char* what, const std::vector<Start>& s, const std::vector<Start>& script, int calls,
                 int report, int size, int go_on, double lnl) {
    std::vector<double> l, t, sl, st;
    std::vector<int> c, sc;
    for (const Start& x : s) {
        l.push_back(x.lnl);
        t.push_back(x.t);
        c.push_back(x.converged);
    }
    for (const Start& x : script) {
        sl.push_back(x.lnl);
        st.push_back(x.t);
        sc.push_back(x.converged);
    }
    int out[4] = {-9, -9, -9, -9};
    double got = 0;
    coal::kaks::detail::gy_continue_starts(l.data(), t.data(), c.data(), (int)s.size(), sl.data(), st.data(),
                                           sc.data(), (int)script.size(), out, &got);
    if (out[0] != calls || out[1] != report || out[2] != size || out[3] != go_on || got != lnl) {
        std::fprintf(stderr,
                     "FAILED: %s: continuations %d report %d starts %d go_on %d lnl %g, expected %d %d %d %d %g\n",
                     what, out[0], out[1], out[2], out[3], got, calls, report, size, go_on, lnl);
        ++bad;
    }
}

}  // namespace

int main() {
    const double T = 50;   // the upper bound of the branch length: a start that ends there is usable
    // one converged start
    expect("one converged start", {{-10, 0.5, 1}}, 0, -1, 0, 0);
    // the highest converged start, wherever it comes
    expect("highest of three", {{-10, 2, 1}, {-9, 3, 1}, {-9.5, 4, 1}}, 1, -1, 1, 0);
    // within 1e-7 the earlier one stays; beyond it the higher one wins
    expect("tie: the earlier start", {{-10, 2, 1}, {-10 + 5e-8, 3, 1}}, 0, -1, 0, 0);
    expect("just above the tie", {{-10, 2, 1}, {-10 + 2e-7, 3, 1}}, 1, -1, 1, 0);
    // an unconverged start on the ridge (t at its bound, or within 0.1% of it) is usable ...
    expect("unconverged ridge start", {{-10, T, 0}}, 0, -1, 0, 0);
    expect("ridge start within 0.1% of the bound", {{-10, T * (1 - 5e-4), 0}}, 0, -1, 0, 0);
    expect("not on the ridge at 1% below the bound", {{-10, T * 0.99, 0}}, -1, 0, 0, 1);
    // ... and wins by its likelihood ...
    expect("ridge start higher than a converged one", {{-10, 2, 1}, {-9, T, 0}}, 1, -1, 1, 0);
    // ... but within the tie a converged start replaces it, also from 1e-7 below, and not the other way round
    expect("tie: converged replaces the ridge start", {{-10, T, 0}, {-10 - 5e-8, 2, 1}}, 1, -1, 1, 0);
    expect("tie: the ridge start does not replace", {{-10, 2, 1}, {-10 + 5e-8, T, 0}}, 0, -1, 0, 0);
    expect("converged start lower than the tie", {{-10, T, 0}, {-10 - 2e-7, 2, 1}}, 0, -1, 0, 0);
    // a start that is not usable and lower: nothing to do
    expect("lower unconverged start", {{-10, 2, 1}, {-11, 3, 0}}, 0, 1, 0, 0);
    expect("unconverged start within the tie", {{-10, 2, 1}, {-10 + 5e-8, 3, 0}}, 0, 1, 0, 0);
    // higher by less than 1e-4: continued, but the usable start is the one reported
    expect("unconverged, 5e-5 higher", {{-10, 2, 1}, {-10 + 5e-5, 3, 0}}, 0, 1, 0, 1);
    expect("unconverged, 5e-5 higher than a ridge start", {{-10, T, 0}, {-10 + 5e-5, 3, 0}}, 0, 1, 0, 1);
    // higher by more than 1e-4: continued, and it is reported (not_converged) if that is how it stays
    expect("unconverged, 2e-4 higher", {{-10, 2, 1}, {-10 + 2e-4, 3, 0}}, 0, 1, 1, 1);
    expect("unconverged, 2e-4 higher than a ridge start", {{-10, T, 0}, {-10 + 2e-4, 3, 0}}, 0, 1, 1, 1);
    expect("unconverged first, 2e-4 higher", {{-10 + 2e-4, 3, 0}, {-10, 2, 1}}, 1, 0, 0, 1);
    // of the starts that are not usable the highest counts (the earlier one of two equal ones)
    expect("highest unconverged start", {{-10, 2, 1}, {-9.5, 3, 0}, {-9, 4, 0}, {-9.2, 5, 0}}, 0, 2, 2, 1);
    expect("equal unconverged starts", {{-9, 3, 0}, {-9, 4, 0}}, -1, 0, 0, 1);
    // no usable start: the best point is reported, after its continuation
    expect("no usable start", {{-10, 2, 0}, {-9, 3, 0}}, -1, 1, 1, 1);

    // The continuation (continue_starts): calls made, start reported, starts kept, go_on at the end, lnl reported.
    // nothing to continue: no call, and the starts are left as they are
    expect_cont("nothing to continue", {{-10, 2, 1}, {-11, 3, 0}}, {{0, 3, 1}}, 0, 0, 2, 0, -10);
    // a continuation that converges higher is reported
    expect_cont("converges higher", {{-10, 2, 1}, {-9, 3, 0}}, {{-8.5, 3, 1}}, 1, 1, 2, 0, -8.5);
    // still not usable and higher each time: two continuations, no third; the last one is reported (higher than
    // the usable start by more than 1e-4)
    expect_cont("two continuations at most", {{-10, 2, 1}, {-9, 3, 0}}, {{-8.9, 3, 0}, {-8.8, 3, 0}, {-8.7, 3, 0}},
                2, 1, 2, 1, -8.8);
    // a continuation that ends lower without becoming usable: the point it started from is kept, and no second
    // continuation is made
    expect_cont("ends lower: keep its start, stop", {{-10, 2, 1}, {-9, 3, 0}}, {{-9.5, 3, 0}, {-8, 3, 1}}, 1, 1, 2,
                1, -9);
    // lower than its start but usable (on the ridge): it is kept
    expect_cont("lower but usable", {{-10, 2, 1}, {-9, 3, 0}}, {{-9.5, T, 0}}, 1, 1, 2, 0, -9.5);
    // the other starts that are not usable are dropped before the continuation: here the one at -9.5 would
    // otherwise be continued next
    expect_cont("the other unusable starts are dropped", {{-10, 2, 1}, {-9, 3, 0}, {-9.5, 4, 0}},
                {{-9.6, 3, 1}, {-1, 3, 1}}, 1, 1, 2, 0, -9.6);
    // no usable start: the best point is continued, and its continuation is reported
    expect_cont("no usable start", {{-10, 2, 0}, {-9, 3, 0}}, {{-8, 3, 1}}, 1, 0, 1, 0, -8);

    // Which pairs are scanned (scan_wanted): 40 codons or fewer always; a larger one if its result did not converge
    // off the ridge, or is within 0.05 lnL per codon of the plateau of saturation.
    struct Gate {
        const char* what;
        size_t n;
        double lnl, t;
        int converged;
        bool want;
    };
    const double plateau = -1000;
    const Gate gates[] = {
        {"40 codons, far from the plateau", 40, plateau + 40 * 0.5, 2, 1, true},
        {"41 codons, far from the plateau", 41, plateau + 41 * 0.5, 2, 1, false},
        {"0.049 per codon above the plateau", 100, plateau + 100 * 0.049, 2, 1, true},
        {"0.051 per codon above the plateau", 100, plateau + 100 * 0.051, 2, 1, false},
        {"on the ridge, 0.049 above", 100, plateau + 100 * 0.049, T, 0, true},
        {"on the ridge, 0.051 above", 100, plateau + 100 * 0.051, T, 0, false},
        {"not converged off the ridge", 100, plateau + 100 * 0.5, 2, 0, true},
        {"converged on the ridge, far", 100, plateau + 100 * 0.5, T, 1, false},
    };
    for (const Gate& g : gates)
        if (coal::kaks::detail::gy_scan_wanted(g.n, g.lnl, g.t, g.converged, plateau) != g.want) {
            std::fprintf(stderr, "FAILED: scan_wanted: %s: expected %d\n", g.what, (int)g.want);
            ++bad;
        }
    std::printf("%s\n", bad ? "failed" : "ok");
    return bad ? 1 : 0;
}
