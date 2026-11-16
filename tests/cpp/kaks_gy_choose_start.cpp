// GY: which of a pair's starts is reported (choose_start in cpp/src/kaks_gy.cpp), with made-up starts. The rule
// that an unconverged start displaces a usable one only when it is higher by more than 1e-4 in lnL is reached by
// no pair of 25,200 random ones, so no test of the command can hold it in place; this one does, through the hook
// the library has for it. Built and run by tests/test_kaks_gy.py against the build's libcoal_core. Exit 0 when
// every check holds; otherwise the failed checks on stderr and exit 1.
#include <cstdio>
#include <vector>

namespace coal {
namespace kaks {
namespace detail {
// n starts by likelihood, branch length and convergence; out: the best usable start, the best start that is not
// usable, the one that is reported (indices; -1: none), and whether the iteration of the second is to be continued
void gy_choose_start(const double* lnl, const double* t, const int* converged, int n, int out[4]);
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
    std::printf("%s\n", bad ? "failed" : "ok");
    return bad ? 1 : 0;
}
