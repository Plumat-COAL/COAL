// Attached values of one-letter options for the small command lines (coal kaks, dnds, translate, codon-dnds and
// coal-cpp-align), which read their arguments by hand: `-j4`, `-j=4`, `-mYN`, `-ofile` and chained flags (`-qj4`)
// mean what `-j 4`, `-m YN`, `-o file` and `-q -j 4` mean, as in the main command line (args.cpp, which follows
// Python's argparse there). Header only: coal-cpp-align is a program of its own.
#pragma once

#include <string>
#include <vector>

namespace coal {

// If `tok` is one dash followed by known one-letter options, the last of which may carry an attached value
// (everything after a value option's letter, less one leading '='), return true and its parts ("-q", "-j", "4").
// value_opts / flag_opts: the letters that take a value / take none. `exact`: single-dash words that are options
// themselves (-codontable) and are left alone. Anything else (a plain `-j`, `--long`, a negative number, an
// unknown letter, a flag followed by '=' or '-') returns false: the caller reports it as it always did.
inline bool split_short_option(const std::string& tok, const std::string& value_opts, const std::string& flag_opts,
                               const std::vector<std::string>& exact, std::vector<std::string>& parts) {
    parts.clear();
    if (tok.size() <= 2 || tok[0] != '-' || tok[1] == '-') return false;
    for (const std::string& e : exact)
        if (tok == e) return false;
    for (size_t k = 1; k < tok.size();) {
        const char ch = tok[k];
        if (value_opts.find(ch) != std::string::npos) {
            parts.push_back(std::string("-") + ch);
            std::string rest = tok.substr(k + 1);
            if (!rest.empty()) parts.push_back(rest[0] == '=' ? rest.substr(1) : rest);
            return true;
        }
        const bool last = k + 1 == tok.size();
        if (flag_opts.find(ch) == std::string::npos || (!last && (tok[k + 1] == '=' || tok[k + 1] == '-'))) {
            parts.clear();
            return false;
        }
        parts.push_back(std::string("-") + ch);
        ++k;
    }
    return true;
}

// The same on an argument vector, at position i (an option position: the caller has consumed everything before
// it): the token is replaced by its parts. Values of earlier options are never looked at, so `-o -q1` keeps its
// file name.
inline void expand_short_option(std::vector<std::string>& args, size_t i, const std::string& value_opts,
                                const std::string& flag_opts, const std::vector<std::string>& exact = {}) {
    std::vector<std::string> parts;
    if (i < args.size() && split_short_option(args[i], value_opts, flag_opts, exact, parts)) {
        args.erase(args.begin() + (std::ptrdiff_t)i);
        args.insert(args.begin() + (std::ptrdiff_t)i, parts.begin(), parts.end());
    }
}

}  // namespace coal
