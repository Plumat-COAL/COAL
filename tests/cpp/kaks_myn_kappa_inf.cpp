// MYN with kappa_AG estimated as infinite, through the C++ API (coal/kaks.hpp): status saturated, a note that names
// kappa, Ka / Ks / sites NaN, the rate ratio infinite, the substitution count kept. Built and run by
// tests/test_kaks_myn_kappa_inf.py against the build's libcoal_core; argv[1], argv[2] are the two rows of a pair of
// tests/kaks/myn_kappa_inf.axt. Exit 0 when every check holds; otherwise the failed check on stderr and exit 1.
#include <coal/kaks.hpp>

#include <cmath>
#include <cstdio>
#include <string>

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: %s ROW1 ROW2\n", argv[0]);
        return 2;
    }
    const coal::GeneticCode* code = coal::get_code(1);
    if (!code) return 2;
    int bad = 0;
    auto check = [&](bool ok, const char* what) {
        if (!ok) {
            std::fprintf(stderr, "FAILED: %s\n", what);
            ++bad;
        }
    };
    const coal::kaks::Result r = coal::kaks::compute(argv[1], argv[2], coal::kaks::Method::MYN, *code);
    check(r.status == coal::kaks::Status::saturated, "status saturated");
    check(r.note.find("kappa_AG estimated as infinite") != std::string::npos, "note names kappa_AG");
    check(std::isnan(r.ka) && std::isnan(r.ks) && std::isnan(r.ka_ks), "Ka, Ks, Ka/Ks NA");
    check(std::isnan(r.S) && std::isnan(r.N), "no site counts");
    check(std::isinf(r.rates[1]) && std::isfinite(r.rates[0]), "rate ratio kappa_AG infinite, kappa_TC finite");
    check(r.substitutions > 0, "substitutions kept");
    check(r.codons_compared == 100, "codons compared");
    // the other methods are unaffected: YN has an estimate
    const coal::kaks::Result y = coal::kaks::compute(argv[1], argv[2], coal::kaks::Method::YN, *code);
    check(y.status == coal::kaks::Status::ok && std::isfinite(y.ka) && std::isfinite(y.ks), "YN estimates");
    // the same through prepare() + compute(Prepared)
    const coal::kaks::Prepared p = coal::kaks::prepare(argv[1], argv[2], *code);
    const coal::kaks::Result r2 = coal::kaks::compute(p, coal::kaks::Method::MYN, *code);
    check(r2.status == r.status && r2.note == r.note && r2.substitutions == r.substitutions, "prepare() path");
    std::printf("%s\n", bad ? "failed" : "ok");
    return bad ? 1 : 0;
}
