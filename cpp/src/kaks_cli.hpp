// `coal-cpp kaks ...` (KaKs_Calculator-compatible output) and `coal-cpp dnds ...` (tidy TSV); docs/kaks.md.
#pragma once
#include <string>
#include <vector>

#include "kaks.hpp"

namespace coal {

// sub: "kaks" or "dnds"; argv: the arguments after the subcommand. Returns the exit code
// (0 done, even if some pairs have a non-ok status; 1 input/output error; 2 command-line error).
int kaks_main(const std::string& sub, const std::vector<std::string>& argv);

// For the other subcommands that write the same formats (dnds_cli.cpp):
// KaKs_Calculator 2.0's header line and one of its rows (values below 1e-6 are NA), and a number as the TSVs print it.
const char* kaks_compat_header();
std::string kaks_compat_row(const std::string& name, const kaks::Result& r, int digits);
std::string kaks_num(double v, int digits);

// --min-codons N / --max-ds X of `coal dnds` and `coal codon-dnds` (docs/kaks.md "Filtering before
// interpretation"): flags for estimates that rest on little data or on a saturated dS. They never drop a row and
// never change a number: a result that has an estimate (status ok or warning) and fewer than N compared codons, or
// a dS above X, gets the note "fewer than N codons compared" / "dS above X" appended, and its status becomes
// warning when it was ok. Results without an estimate (saturated, not_converged, no_data, error) are left as they
// are. Both are off by default.
struct PairFlags {
    long min_codons = 0;       // 0: off
    double max_ds = -1;        // < 0: off
    long n_short = 0, n_high_ds = 0;   // results flagged (counted by apply; call it from one thread)
    bool on() const { return min_codons > 0 || max_ds >= 0; }
    // Reads the value of --min-codons / --max-ds; false with the reason in `why`.
    bool parse(const std::string& option, const std::string& value, std::string& why);
    void apply(kaks::Status& status, std::string& note, long codons_compared, double ks);
    // "3 result(s) flagged: 2 with fewer than 50 codons compared, 1 with dS above 2" ("" when off).
    std::string summary() const;
};

}  // namespace coal
