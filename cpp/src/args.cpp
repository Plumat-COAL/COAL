// argparse-compatible command-line parsing (CPP-4). The structure follows the documented behaviour of
// Python 3.10's argparse (parse_known_args): first every argument is classified as an option ('O'),
// a value ('A') or the '--' marker ('-'); then options and runs of positionals are consumed left to
// right. Positional runs are matched greedily, like argparse's nargs patterns "-*A-*" (one value) and
// "-*A[A-]*" (one or more). Written from the behaviour, not copied.
#include "args.hpp"

#include <climits>
#include <cstdio>

namespace coal {
namespace args {

namespace {

struct ArgError {
    std::string msg;
};
struct Done {
    int code;
    std::string out;
};

bool is_space(unsigned char c) { return c == ' ' || (c >= 9 && c <= 13) || (c >= 0x1c && c <= 0x1f); }
bool is_digit(char c) { return c >= '0' && c <= '9'; }
bool starts_with(const std::string& s, const std::string& p) { return s.compare(0, p.size(), p) == 0; }

// argparse's negative-number test: ^-\d+$|^-\d*\.\d+$  ('$' also matches before a final newline)
bool negative_number(std::string a) {
    if (!a.empty() && a.back() == '\n') a.pop_back();
    if (a.size() < 2 || a[0] != '-') return false;
    size_t i = 1, d1 = 0, d2 = 0;
    while (i < a.size() && is_digit(a[i])) ++i, ++d1;
    if (i == a.size()) return d1 > 0;
    if (a[i] != '.') return false;
    ++i;
    while (i < a.size() && is_digit(a[i])) ++i, ++d2;
    return i == a.size() && d2 > 0;
}

struct Opt {
    int spec = -1;  // -1: looks like an option but is unknown
    std::string opt;
    bool has_exp = false;
    std::string exp;
};

// one element of a positional nargs pattern: '-' (the '--' marker), 'A' (a value), 'B' (either)
struct El {
    char cls;
    bool star;
    int group;  // >= 0: last element of this positional's group
};

bool match_els(const std::vector<El>& els, size_t ei, const std::string& s, size_t si, std::vector<size_t>& ends) {
    if (ei == els.size()) return true;
    const El& e = els[ei];
    auto ok = [&](char c) { return e.cls == 'B' ? (c == 'A' || c == '-') : c == e.cls; };
    if (!e.star) {
        if (si < s.size() && ok(s[si])) {
            if (e.group >= 0) ends[(size_t)e.group] = si + 1;
            return match_els(els, ei + 1, s, si + 1, ends);
        }
        return false;
    }
    size_t k = si;
    while (k < s.size() && ok(s[k])) ++k;
    for (size_t j = k + 1; j-- > si;) {  // greedy first, then give back (regex backtracking)
        if (e.group >= 0) ends[(size_t)e.group] = j;
        if (match_els(els, ei + 1, s, j, ends)) return true;
    }
    return false;
}

}  // namespace

bool py_int(const std::string& s, std::string& norm) {
    size_t b = 0, e = s.size();
    while (b < e && is_space((unsigned char)s[b])) ++b;
    while (e > b && is_space((unsigned char)s[e - 1])) --e;
    if (b == e) return false;
    bool neg = false;
    if (s[b] == '+' || s[b] == '-') {
        neg = s[b] == '-';
        ++b;
    }
    if (b == e) return false;
    std::string digits;
    bool prev_digit = false;
    for (size_t i = b; i < e; ++i) {
        char c = s[i];
        if (is_digit(c)) {
            digits += c;
            prev_digit = true;
        } else if (c == '_' && prev_digit && i + 1 < e && is_digit(s[i + 1])) {
            prev_digit = false;
        } else {
            return false;
        }
    }
    size_t z = digits.find_first_not_of('0');
    digits = z == std::string::npos ? "0" : digits.substr(z);
    norm = (neg && digits != "0" ? "-" : "") + digits;
    return true;
}

std::string py_repr(const std::string& s) {
    char q = (s.find('\'') != std::string::npos && s.find('"') == std::string::npos) ? '"' : '\'';
    std::string r(1, q);
    for (char ch : s) {
        unsigned char c = (unsigned char)ch;
        if (c == '\\') {
            r += "\\\\";
        } else if (ch == q) {
            r += '\\';
            r += ch;
        } else if (c == '\n') {
            r += "\\n";
        } else if (c == '\r') {
            r += "\\r";
        } else if (c == '\t') {
            r += "\\t";
        } else if (c < 0x20 || c == 0x7f) {
            char buf[8];
            std::snprintf(buf, sizeof buf, "\\x%02x", (unsigned)c);
            r += buf;
        } else {
            r += ch;
        }
    }
    r += q;
    return r;
}

struct ParseRun {
    const Parser& p;
    const std::vector<std::string>& argv;
    Result& res;
    std::vector<std::pair<std::string, int>> order;  // option strings in registration order
    std::string pat;
    std::vector<Opt> opts;
    std::vector<int> pos_left;
    std::vector<bool> seen;
    std::vector<std::string> extras;

    ParseRun(const Parser& parser, const std::vector<std::string>& a, Result& r) : p(parser), argv(a), res(r) {
        for (size_t i = 0; i < p.specs_.size(); ++i) {
            const Spec& s = p.specs_[i];
            if (s.kind == Kind::Positional || s.kind == Kind::PositionalPlus)
                pos_left.push_back((int)i);
            else
                for (auto& f : s.flags) order.emplace_back(f, (int)i);
        }
        seen.assign(p.specs_.size(), false);
    }

    const Spec& spec(int i) const { return p.specs_[(size_t)i]; }

    std::string name(int i) const {
        const Spec& s = spec(i);
        if (s.flags.empty()) return s.name;
        std::string n;
        for (auto& f : s.flags) n += (n.empty() ? "" : "/") + f;
        return n;
    }

    [[noreturn]] void fail(int i, const std::string& msg) const { throw ArgError{"argument " + name(i) + ": " + msg}; }

    int exact(const std::string& a) const {
        for (auto& o : order)
            if (o.first == a) return o.second;
        return -1;
    }

    std::vector<Opt> tuples(const std::string& a) const {
        std::vector<Opt> t;
        if (a[1] == '-') {  // --long: unique prefixes, optional =value
            std::string prefix = a, exp;
            bool has = false;
            size_t eq = a.find('=');
            if (eq != std::string::npos) {
                prefix = a.substr(0, eq);
                exp = a.substr(eq + 1);
                has = true;
            }
            for (auto& o : order)
                if (spec(o.second).abbrev && starts_with(o.first, prefix)) t.push_back({o.second, o.first, has, exp});
        } else {  // -x: a two-character option with the rest attached, or a prefix of a single-dash option
            std::string sp = a.substr(0, 2), sexp = a.substr(2);
            for (auto& o : order) {
                if (!spec(o.second).abbrev) continue;
                if (o.first == sp)
                    t.push_back({o.second, o.first, true, sexp});
                else if (starts_with(o.first, a))
                    t.push_back({o.second, o.first, false, ""});
            }
        }
        return t;
    }

    bool classify(const std::string& a, Opt& o) const {
        if (a.empty() || a[0] != '-') return false;
        int i = exact(a);
        if (i >= 0) {
            o = {i, a, false, ""};
            return true;
        }
        if (a.size() == 1) return false;
        size_t eq = a.find('=');
        if (eq != std::string::npos) {
            std::string n = a.substr(0, eq);
            int j = exact(n);
            if (j >= 0) {
                o = {j, n, true, a.substr(eq + 1)};
                return true;
            }
        }
        std::vector<Opt> t = tuples(a);
        if (t.size() > 1) {
            std::string m;
            for (auto& x : t) m += (m.empty() ? "" : ", ") + x.opt;
            throw ArgError{"ambiguous option: " + a + " could match " + m};
        }
        if (t.size() == 1) {
            o = t[0];
            return true;
        }
        if (negative_number(a)) return false;  // no option of ours looks like a negative number
        if (a.find(' ') != std::string::npos) return false;
        o = {-1, a, false, ""};
        return true;
    }

    std::string convert(int i, const std::string& v) const {
        const Spec& s = spec(i);
        std::string out = v;
        if (s.is_int && !py_int(v, out)) fail(i, "invalid int value: " + py_repr(v));
        if (!s.choices.empty()) {
            bool ok = false;
            for (auto& c : s.choices) ok = ok || c == out;
            if (!ok) {
                std::string list;
                for (auto& c : s.choices) list += (list.empty() ? "" : ", ") + (s.is_int ? c : py_repr(c));
                fail(i, "invalid choice: " + (s.is_int ? out : py_repr(out)) + " (choose from " + list + ")");
            }
        }
        return out;
    }

    void take(int i, std::vector<std::string> vals) {
        seen[(size_t)i] = true;
        const Spec& s = spec(i);
        switch (s.kind) {
            case Kind::Help:
                throw Done{0, p.help_};
            case Kind::Version:
                throw Done{0, p.version_ + "\n"};
            case Kind::Flag:
                res.values[(size_t)i] = {"1"};
                break;
            case Kind::Value:
                res.values[(size_t)i] = {convert(i, vals.at(0))};
                break;
            case Kind::Positional:
            case Kind::PositionalPlus:
                for (size_t k = 0; k < vals.size(); ++k)
                    if (vals[k] == "--") {  // argparse drops one '--' from a positional's values
                        vals.erase(vals.begin() + (long)k);
                        break;
                    }
                res.values[(size_t)i] = vals;
                break;
        }
    }

    std::vector<size_t> match_partial(const std::string& s) const {
        for (size_t n = pos_left.size(); n > 0; --n) {
            std::vector<El> els;
            for (size_t g = 0; g < n; ++g) {
                bool plus = spec(pos_left[g]).kind == Kind::PositionalPlus;
                els.push_back({'-', true, -1});
                els.push_back({'A', false, -1});
                els.push_back({plus ? 'B' : '-', true, (int)g});
            }
            std::vector<size_t> ends(n, 0);
            if (match_els(els, 0, s, 0, ends)) {
                std::vector<size_t> counts;
                size_t prev = 0;
                for (size_t e : ends) {
                    counts.push_back(e - prev);
                    prev = e;
                }
                return counts;
            }
        }
        return {};
    }

    size_t consume_positionals(size_t start) {
        std::vector<size_t> counts = match_partial(pat.substr(start));
        for (size_t k = 0; k < counts.size(); ++k) {
            std::vector<std::string> vals(argv.begin() + (long)start, argv.begin() + (long)(start + counts[k]));
            start += counts[k];
            take(pos_left[k], vals);
        }
        pos_left.erase(pos_left.begin(), pos_left.begin() + (long)counts.size());
        return start;
    }

    size_t consume_optional(size_t start) {
        Opt o = opts[start];
        std::vector<std::pair<int, std::vector<std::string>>> acts;
        size_t stop = start + 1;
        for (;;) {
            if (o.spec < 0) {
                extras.push_back(argv[start]);
                return start + 1;
            }
            int nargs = spec(o.spec).kind == Kind::Value ? 1 : 0;
            if (o.has_exp) {
                if (nargs == 0 && o.opt[1] != '-') {  // -qq: chained single-character flags
                    acts.push_back({o.spec, {}});
                    if (o.exp.empty()) fail(o.spec, "ignored explicit argument " + py_repr(o.exp));
                    std::string next = std::string("-") + o.exp[0];
                    int j = exact(next);
                    if (j < 0) fail(o.spec, "ignored explicit argument " + py_repr(o.exp));
                    std::string rest = o.exp.substr(1);
                    o = {j, next, !rest.empty(), rest};
                } else if (nargs == 1) {
                    stop = start + 1;
                    acts.push_back({o.spec, {o.exp}});
                    break;
                } else {
                    fail(o.spec, "ignored explicit argument " + py_repr(o.exp));
                }
            } else {
                size_t s0 = start + 1, cnt = 0;
                if (nargs == 1) {
                    if (s0 < pat.size() && pat[s0] == 'A')
                        cnt = 1;
                    else
                        fail(o.spec, "expected one argument");
                }
                stop = s0 + cnt;
                acts.push_back({o.spec, std::vector<std::string>(argv.begin() + (long)s0, argv.begin() + (long)stop)});
                break;
            }
        }
        for (auto& a : acts) take(a.first, a.second);
        return stop;
    }

    void run() {
        size_t n = argv.size();
        pat.assign(n, 'A');
        opts.assign(n, Opt());
        bool rest = false;
        for (size_t i = 0; i < n; ++i) {
            if (rest) continue;
            if (argv[i] == "--") {
                pat[i] = '-';
                rest = true;
            } else if (classify(argv[i], opts[i])) {
                pat[i] = 'O';
            }
        }
        size_t last = pat.rfind('O');
        size_t start = 0;
        while (last != std::string::npos && start <= last) {
            size_t next = pat.find('O', start);
            if (start != next) {
                size_t end = consume_positionals(start);
                if (end > start) {
                    start = end;
                    continue;
                }
                start = end;
            }
            if (pat[start] != 'O') {
                extras.insert(extras.end(), argv.begin() + (long)start, argv.begin() + (long)next);
                start = next;
            }
            start = consume_optional(start);
        }
        size_t stop = consume_positionals(start);
        extras.insert(extras.end(), argv.begin() + (long)stop, argv.end());

        std::string missing;
        for (size_t i = 0; i < p.specs_.size(); ++i) {
            Kind k = p.specs_[i].kind;
            if ((k == Kind::Positional || k == Kind::PositionalPlus) && !seen[i])
                missing += (missing.empty() ? "" : ", ") + p.specs_[i].name;
        }
        if (!missing.empty()) throw ArgError{"the following arguments are required: " + missing};
        if (!extras.empty()) {
            std::string e;
            for (auto& x : extras) e += (e.empty() ? "" : " ") + x;
            throw ArgError{"unrecognized arguments: " + e};
        }
    }
};

Result Parser::parse(const std::vector<std::string>& argv) const {
    Result res;
    res.values.resize(specs_.size());
    try {
        ParseRun(*this, argv, res).run();
    } catch (const ArgError& e) {
        res.exit = true;
        res.code = 2;
        res.err = error_text(e.msg);
    } catch (const Done& d) {
        res.exit = true;
        res.code = d.code;
        res.out = d.out;
    }
    return res;
}

}  // namespace args
}  // namespace coal
