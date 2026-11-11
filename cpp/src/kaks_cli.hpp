// `coal-cpp kaks ...` (KaKs_Calculator-compatible output) and `coal-cpp dnds ...` (tidy TSV); docs/kaks.md.
#pragma once
#include <string>
#include <vector>

namespace coal {

// sub: "kaks" or "dnds"; argv: the arguments after the subcommand. Returns the exit code
// (0 done, even if some pairs have a non-ok status; 1 input/output error; 2 command-line error).
int kaks_main(const std::string& sub, const std::vector<std::string>& argv);

}  // namespace coal
