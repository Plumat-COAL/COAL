// `coal-cpp translate ...` and `coal-cpp codon-dnds ...`: the native steps of `coal pipeline` (docs/pipeline.md).
#pragma once
#include <string>
#include <vector>

namespace coal {

// argv: the arguments after the subcommand. Exit codes: 0 done (genes / pairs / groups that could not be used are
// reported in the output files, they do not change the exit code); 1 input/output error; 2 command-line error.
int translate_main(const std::vector<std::string>& argv);
int codon_dnds_main(const std::vector<std::string>& argv);

}  // namespace coal
