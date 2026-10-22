// Command-line parsing that behaves like Python 3.13's argparse, so coal-cpp accepts the same
// command lines as `python -m coal` and prints the same error messages. Covers what COAL's
// parser uses: flags, one-value options (string or int, optional choices), one positional with one
// value and one with one-or-more; long-option abbreviations, single-dash options with attached values,
// `=` values, chained short flags (-qq), negative numbers as values, `--`.
#pragma once

#include <string>
#include <vector>

namespace coal {
namespace args {

enum class Kind { Flag, Value, Help, Version, Positional, PositionalPlus };

struct Spec {
    Kind kind;
    std::vector<std::string> flags;    // option strings in registration order (empty for positionals)
    std::string name;                  // positional name (pep_aln); unused for options
    bool is_int = false;               // type=int
    std::vector<std::string> choices;  // allowed values (ints as decimal strings)
    bool abbrev = true;                // false: exact option names only (C++ extras)
};

struct Result {
    bool exit = false;  // help/version or error: print `out` / `err` and exit with `code`
    int code = 0;
    std::string out, err;
    // per spec (same index as add()): flags -> "1" when given; values -> last value (ints normalised);
    // positionals -> all values
    std::vector<std::vector<std::string>> values;
    bool given(int i) const { return !values[(size_t)i].empty(); }
    const std::string& get(int i) const { return values[(size_t)i].back(); }
};

class Parser {
public:
    Parser(std::string prog, std::string usage, std::string help, std::string version)
        : prog_(std::move(prog)), usage_(std::move(usage)), help_(std::move(help)), version_(std::move(version)) {}
    int add(Spec s) {
        specs_.push_back(std::move(s));
        return (int)specs_.size() - 1;
    }
    Result parse(const std::vector<std::string>& argv) const;
    // "usage: ...\n<prog>: error: <msg>\n" (what argparse prints before exiting with 2)
    std::string error_text(const std::string& msg) const { return usage_ + prog_ + ": error: " + msg + "\n"; }

private:
    std::string prog_, usage_, help_, version_;
    std::vector<Spec> specs_;
    friend struct ParseRun;
};

// Python int(): surrounding whitespace, sign, '_' between digits; returns false if not an int.
// `norm` gets the canonical decimal form ("007" -> "7", "1_0" -> "10", "-0" -> "0").
bool py_int(const std::string& s, std::string& norm);
// Python repr() of a str (quote choice and escapes), as argparse uses in messages.
std::string py_repr(const std::string& s);

}  // namespace args
}  // namespace coal
