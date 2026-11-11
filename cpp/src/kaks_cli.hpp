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

}  // namespace coal
