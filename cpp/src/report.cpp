#include "report.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>

namespace coal {

const char* const TABLE_NAMES[24] = {
    nullptr, "standard", "vertebrate mitochondrial", "yeast mitochondrial",
    "mold, protozoan and coelenterate mitochondrial", "invertebrate mitochondrial",
    "ciliate, dasycladacean and Hexamita nuclear", nullptr, nullptr, "echinoderm and flatworm mitochondrial",
    "euplotid nuclear", "bacterial, archaeal and plant plastid", "alternative yeast nuclear",
    "ascidian mitochondrial", "alternative flatworm mitochondrial", "Blepharisma nuclear",
    "chlorophycean mitochondrial", nullptr, nullptr, nullptr, nullptr, "trematode mitochondrial",
    "Scenedesmus obliquus mitochondrial", "Thraustochytrium mitochondrial"};

static const char* VERSION = "0.4.0";
static const char* REPO_URL = "https://github.com/chloe-newt/COAL";
static const int HINT_SAMPLE = 200;
static const double HINT_THRESHOLD = 0.20;

std::string fmt_n(long long x) {
    std::string s = std::to_string(x < 0 ? -x : x), out;
    int c = 0;
    for (size_t i = s.size(); i-- > 0;) {
        out += s[i];
        if (++c % 3 == 0 && i > 0) out += ',';
    }
    if (x < 0) out += '-';
    std::reverse(out.begin(), out.end());
    return out;
}

static std::string seqs(long long n) { return fmt_n(n) + " sequence" + (n == 1 ? "" : "s"); }

static std::string shlex_quote(const std::string& s) {
    if (s.empty()) return "''";
    bool safe = true;
    for (unsigned char c : s)
        if (!(isalnum(c) || c == '_' || c == '@' || c == '%' || c == '+' || c == '=' || c == ':' || c == ',' ||
              c == '.' || c == '/' || c == '-') || c > 127) {
            safe = false;
            break;
        }
    if (safe) return s;
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') out += "'\"'\"'";
        else out += c;
    }
    return out + "'";
}

std::string shlex_join(const std::vector<std::string>& args) {
    std::string out;
    for (size_t i = 0; i < args.size(); ++i) {
        if (i) out += ' ';
        out += shlex_quote(args[i]);
    }
    return out;
}

static std::string replace_all(std::string s, const std::string& a, const std::string& b) {
    size_t p = 0;
    while ((p = s.find(a, p)) != std::string::npos) {
        s.replace(p, a.size(), b);
        p += b.size();
    }
    return s;
}

std::string plain_reason(const std::string& code, const std::string& message) {
    std::string detail = message;
    if ((code == "invalid_chars" || code == "no_match")) {
        size_t p = message.find(": ");
        if (p != std::string::npos) detail = message.substr(p + 2);
    }
    if (code == "count_mismatch" || code == "row_lengths") {
        // last "(...)" group at the end of the message
        size_t end = message.find_last_not_of(" \t\n\r");
        if (end != std::string::npos && message[end] == ')') {
            size_t open = message.rfind('(', end);
            if (open != std::string::npos) {
                std::string inner = message.substr(open + 1, end - open - 1);
                if (inner.find('(') == std::string::npos && inner.find(')') == std::string::npos) detail = inner;
            }
        }
    }
    if (code == "invalid_chars") {
        size_t p = message.rfind("characters ");
        detail = p == std::string::npos ? message : message.substr(p + 11);
    }
    std::string tpl = "{detail}";
    if (code == "no_match")
        tpl = "The protein could not be found in its DNA sequence in any reading frame, even allowing for a few "
              "wrong codons. Check that this is the right DNA for this protein and the genetic code (-codontable).";
    else if (code == "invalid_chars")
        tpl = "The protein alignment contains characters that are not amino acids ({detail}).";
    else if (code == "empty_protein")
        tpl = "This protein row contains only gaps.";
    else if (code == "count_mismatch")
        tpl = "The protein alignment and the DNA file(s) contain different numbers of sequences ({detail}); every "
              "protein needs exactly one DNA sequence.";
    else if (code == "no_sequences")
        tpl = "No sequences were found in the input files.";
    else if (code == "row_lengths")
        tpl = "The protein alignment rows have different lengths ({detail}), so it is not an alignment.";
    else if (code == "input_format")
        tpl = "An input file could not be read: {detail}.";
    else if (code == "file_error")
        tpl = "A file could not be opened or written: {detail}.";
    else if (code == "all_failed")
        tpl = "None of the sequences could be converted, so no output was written.";
    return replace_all(tpl, "{detail}", detail);
}

// ---------------------------------------------------------------- Entry

Json Entry::to_json() const {
    Json j = Json::object();
    j.set("id", id);
    j.set("index", index);
    j.set("dna_id", has_dna_id ? Json(dna_id) : Json());
    j.set("pairing", pairing.empty() ? Json() : Json(pairing));
    j.set("status", status);
    Json ms = Json::array();
    for (auto& m : messages) {
        Json x = Json::object();
        x.set("code", m.code);
        if (m.column) x.set("column", m.column);
        x.set("text", m.text);
        ms.push(x);
    }
    j.set("messages", ms);
    return j;
}

Entry Entry::from_json(const Json& j) {
    Entry e;
    if (auto v = j.get("id")) e.id = v->as_str();
    if (auto v = j.get("index")) e.index = v->as_int();
    if (auto v = j.get("dna_id"); v && v->type == Json::String) {
        e.dna_id = v->s;
        e.has_dna_id = true;
    }
    if (auto v = j.get("pairing")) e.pairing = v->as_str();
    if (auto v = j.get("status")) e.status = v->as_str();
    if (auto v = j.get("messages"))
        for (auto& m : v->a) {
            Message x;
            if (auto c = m.get("code")) x.code = c->as_str();
            if (auto c = m.get("column")) x.column = (int)c->as_int();
            if (auto c = m.get("text")) x.text = c->as_str();
            e.messages.push_back(x);
        }
    return e;
}

// ---------------------------------------------------------------- RunReport

RunReport::RunReport(const std::vector<std::string>& argv, const std::string& pep,
                     const std::vector<std::string>& nucs, const std::string& out, const std::string& fmt, int table,
                     bool nogap, bool nomismatch, bool keep_going, bool keep_entries_)
    : keep_entries(keep_entries_), t0_(std::chrono::steady_clock::now()), pep_(pep), out_(out), fmt_(fmt),
      nucs_(nucs), table_(table), nogap_(nogap), nomismatch_(nomismatch), keep_going_(keep_going) {
    std::vector<std::string> c{"coal"};
    c.insert(c.end(), argv.begin(), argv.end());
    command_ = shlex_join(c);
    std::vector<std::string> p{"pal2nal.pl", pep};
    p.insert(p.end(), nucs.begin(), nucs.end());
    p.push_back("-output");
    p.push_back(fmt);
    if (table != 1) {
        p.push_back("-codontable");
        p.push_back(std::to_string(table));
    }
    if (nogap) p.push_back("-nogap");
    if (nomismatch) p.push_back("-nomismatch");
    pal2nal_ = shlex_join(p);
}

void RunReport::set_inputs(long long n_protein_, const std::vector<long long>& dna_counts) {
    n_protein = n_protein_;
    dna_counts_ = dna_counts;
    has_dna_counts_ = true;
}

void RunReport::set_pairing(const std::string& mode, const std::vector<std::string>* dna_ids,
                            const std::vector<size_t>*, const std::vector<std::string>& per_seq,
                            const std::string& ex_prot, const std::string& ex_dna) {
    pairing_ = mode;
    pairing_set = true;
    dna_ids_ = dna_ids;
    seq_pairing_ = per_seq;
    ex_prot_ = ex_prot;
    ex_dna_ = ex_dna;
}

static bool parse_warning(const std::string& w, const char* tail, int& col, char& aa, std::string& rest) {
    // "pepAlnPos <col>: <aa><tail><rest>"
    const std::string head = "pepAlnPos ";
    if (w.compare(0, head.size(), head) != 0) return false;
    size_t p = head.size(), q = p;
    while (q < w.size() && isdigit((unsigned char)w[q])) ++q;
    if (q == p || w.compare(q, 2, ": ") != 0 || q + 3 > w.size()) return false;
    std::string t = tail;
    if (w.compare(q + 3, t.size(), t) != 0) return false;
    col = std::stoi(w.substr(p, q - p));
    aa = w[q + 2];
    rest = w.substr(q + 3 + t.size());
    return true;
}

const Entry* RunReport::add_sequence(const std::string& pid, const std::vector<std::string>& warnings,
                                     const std::string& fail_code, const std::string& fail_msg, long long k) {
    bool failure = !fail_code.empty();
    if (warnings.empty() && !failure && !keep_entries) {
        ++n_ok;
        return nullptr;
    }
    static Entry scratch;
    Entry e;
    e.id = pid;
    e.index = k;
    bool mismatch = false, unknown = false;
    for (auto& w : warnings) {
        int col;
        char aa;
        std::string rest;
        if (parse_warning(w, " does not correspond to ", col, aa, rest) && !rest.empty() &&
            rest.find(' ') == std::string::npos) {
            mismatch = true;
            e.messages.push_back({"mismatch", "Column " + std::to_string(col) + ": " + aa + " does not match codon " +
                                                  rest, col});
            continue;
        }
        if (parse_warning(w, " unknown AA type", col, aa, rest)) {
            unknown = true;
            e.messages.push_back({"unknown_aa", "Column " + std::to_string(col) + ": " + aa +
                                                    " is not a standard amino acid; treated as X (any codon)", col});
            continue;
        }
        e.messages.push_back({"warning", w, 0});
    }
    if (failure) {
        e.messages.push_back({fail_code, plain_reason(fail_code, fail_msg), 0});
        e.status = "failed";
    } else {
        e.status = e.messages.empty() ? "ok" : "warning";
    }
    if ((failure || mismatch) && problems.size() < (size_t)HINT_SAMPLE) problems.push_back(k);
    if (e.status == "ok") ++n_ok;
    else if (e.status == "warning") ++n_warning;
    else ++n_failed;
    n_mismatch += mismatch;
    n_unknown += unknown;
    if (dna_ids_ && k >= 0 && (size_t)k < dna_ids_->size()) {
        e.dna_id = (*dna_ids_)[(size_t)k];
        e.has_dna_id = true;
    }
    e.pairing = !seq_pairing_.empty() && (size_t)k < seq_pairing_.size() ? seq_pairing_[(size_t)k] : pairing_;
    if (failure) failed_entries_.push_back(e);
    if (keep_entries) {
        entries_.push_back(std::move(e));
        return &entries_.back();
    }
    scratch = std::move(e);
    return &scratch;
}

void RunReport::fail(const std::string& code, const std::string& message) {
    has_error = true;
    error_code = code;
    error_text = plain_reason(code, message);
    error_detail = message;
}

void RunReport::set_output(const std::string& path, const std::string& sha, long long n) {
    has_output_ = true;
    out_path_ = path;
    out_sha_ = sha;
    out_n_ = n;
}

bool RunReport::needs_hint() const {
    long long done = n_ok + n_warning + n_failed;
    long long mism = n_mismatch + n_failed;
    if (!done || (double)mism / (double)done <= HINT_THRESHOLD || problems.empty() || !pairing_set) return false;
    return true;
}

void RunReport::apply_fits(const std::vector<std::pair<int, int>>& fits) {
    if (fits.empty()) return;
    int n = (int)problems.size();
    int cur = 0;
    for (auto& f : fits)
        if (f.first == table_) cur = f.second;
    // best: most fits, then lowest table number
    std::pair<int, int> best = fits[0];
    for (auto& f : fits)
        if (f.second > best.second || (f.second == best.second && f.first < best.first)) best = f;
    int t = best.first, k = best.second;
    if (t == table_ || (double)k < 0.8 * n || (double)(k - cur) < 0.2 * n) return;
    std::string name = TABLE_NAMES[t];
    std::string quant = k == n ? "They all fit" : fmt_n(k) + " of the " + fmt_n(n) + " checked fit";
    Hint h;
    h.table = t;
    h.checked = n;
    h.fit = k;
    h.fit_current = cur;
    h.text = "Most sequences fit genetic code " + std::to_string(t) + " (" + name + "). Try -codontable " +
             std::to_string(t) + ".";
    h.summary = quant + " genetic code " + std::to_string(t) + " (" + name + "): try -codontable " +
                std::to_string(t) + ".";
    hints.push_back(h);
}

std::vector<std::string> RunReport::summary_lines() const {
    long long done = n_ok + n_warning + n_failed;
    long long total = n_protein >= 0 ? n_protein : done;
    std::vector<std::string> lines;
    const Hint* hint = hints.empty() ? nullptr : &hints[0];
    if (n_mismatch) {
        std::string s = fmt_n(n_mismatch) + " of " + fmt_n(total) + " sequences " + (n_mismatch == 1 ? "has" : "have") +
                        " codon mismatches.";
        if (hint) s += " " + hint->summary;
        lines.push_back(s);
    }
    const auto& failed = failed_entries_;
    if (!failed.empty()) {
        const Entry& first = failed[0];
        std::string left_out;
        if (keep_going_ && !has_error) left_out = std::string(" and ") + (failed.size() == 1 ? "was" : "were") + " left out";
        std::string s = seqs((long long)failed.size()) + " of " + fmt_n(total) + " could not be converted" + left_out +
                        " (first: " + first.id + ": " + first.messages.back().text + ")";
        if (hint && !n_mismatch) s += " " + hint->summary;
        lines.push_back(s);
    }
    if (has_error && failed.empty()) lines.push_back("The run failed: " + error_text);
    if (pairing_ == "order" && done) {
        lines.push_back("Protein and DNA IDs don't match, so sequences were paired by file order. "
                        "Check that the files list the same genes in the same order.");
    } else if (pairing_ == "normalised" && done) {
        long long k = std::count(seq_pairing_.begin(), seq_pairing_.end(), std::string("normalised"));
        std::string ex = ex_prot_.empty() ? "" : ", e.g. protein " + ex_prot_ + " = DNA " + ex_dna_;
        lines.push_back(fmt_n(k) + " of " + fmt_n(total) +
                        " IDs matched only after normalising (version suffix, case or pipe fields)" + ex + ".");
    }
    if (n_unknown)
        lines.push_back(fmt_n(n_unknown) + " of " + fmt_n(total) +
                        " sequences contain non-standard amino acids (treated as X).");
    if (has_error && n_protein > 0 && done < n_protein && !failed.empty()) {
        long long left = n_protein - done;
        lines.push_back("The run stopped there: " + seqs(left) + " " + (left == 1 ? "was" : "were") +
                        " not processed and no output was written.");
    }
    return lines;
}

std::string RunReport::methods_text() const {
    long long n = has_output_ ? out_n_ : n_ok + n_warning + n_failed;
    int t = table_;
    std::string code = t == 1 ? "standard genetic code" : std::string(TABLE_NAMES[t]) + " genetic code";
    std::vector<std::string> opts;
    if (nogap_) opts.push_back("removing codon columns containing gaps or stop codons (-nogap)");
    if (nomismatch_) opts.push_back("removing columns with codons that do not match their amino acid (-nomismatch)");
    std::string opt_txt;
    if (!opts.empty()) opt_txt = ", " + (opts.size() == 2 ? opts[0] + " and " + opts[1] : opts[0]);
    std::string skipped;
    if (n_failed && keep_going_ && !has_error)
        skipped = " (" + seqs(n_failed) + " that could not be matched to " + (n_failed == 1 ? "its" : "their") +
                  " DNA " + (n_failed == 1 ? "was" : "were") + " excluded)";
    std::string verified = "outputs are verified byte-identical to PAL2NAL v14 on the COAL test suite";
    if (t == 10 || t == 11 || t == 13)
        verified += " (COAL follows the current NCBI definition of table " + std::to_string(t) +
                    ", a documented difference)";
    std::string fname = fmt_ == "fasta" ? "FASTA" : fmt_ == "clustal" ? "CLUSTAL" : fmt_ == "paml" ? "PAML"
                                                                                   : "codon-table (PAL2NAL '-output codon')";
    return std::string("Codon alignments were generated with COAL v") + VERSION + " (" + REPO_URL + ") from " +
           fmt_n(n) + " aligned protein sequence" + (n == 1 ? "" : "s") + " and their coding sequences" + skipped +
           ", using the " + code + " (NCBI table " + std::to_string(t) + ") and " + fname + " output" + opt_txt +
           "; " + verified + ".";
}

Json RunReport::to_json() const {
    long long done = n_ok + n_warning + n_failed;
    long long total = n_protein >= 0 ? n_protein : done;
    Json d = Json::object();
    d.set("schema", "coal.report/1");
    d.set("coal_version", VERSION);
    d.set("command", command_);
    d.set("pal2nal_equivalent", pal2nal_);
    Json opts = Json::object();
    opts.set("output_format", fmt_);
    opts.set("codontable", table_);
    opts.set("nogap", nogap_);
    opts.set("nomismatch", nomismatch_);
    opts.set("keep_going", keep_going_);
    d.set("options", opts);
    Json inputs = Json::array();
    auto sha_or_null = [](const std::string& p) {
        std::string h = sha256_file(p);
        return h.empty() ? Json() : Json(h);
    };
    Json pin = Json::object();
    pin.set("role", "protein");
    pin.set("path", pep_);
    pin.set("sha256", sha_or_null(hash_paths_.size() == nucs_.size() + 1 ? hash_paths_[0] : pep_));
    pin.set("sequences", n_protein >= 0 ? Json(n_protein) : Json());
    inputs.push(pin);
    for (size_t i = 0; i < nucs_.size(); ++i) {
        Json x = Json::object();
        x.set("role", "dna");
        x.set("path", nucs_[i]);
        x.set("sha256", sha_or_null(hash_paths_.size() == nucs_.size() + 1 ? hash_paths_[i + 1] : nucs_[i]));
        x.set("sequences", has_dna_counts_ ? Json(dna_counts_[i]) : Json());
        inputs.push(x);
    }
    d.set("inputs", inputs);
    Json out = Json::object();
    if (has_output_) {
        out.set("path", out_path_.empty() ? Json() : Json(out_path_));
        out.set("format", fmt_);
        out.set("sha256", out_sha_.empty() ? Json() : Json(out_sha_));
        out.set("sequences", out_n_);
    } else {
        out.set("path", out_.empty() ? Json() : Json(out_));
        out.set("format", fmt_);
        out.set("sha256", Json());
        out.set("sequences", 0);
    }
    d.set("output", out);
    Json sum = Json::object();
    sum.set("total", total);
    sum.set("ok", n_ok);
    sum.set("warning", n_warning);
    sum.set("failed", n_failed);
    sum.set("pairing", pairing_.empty() ? Json() : Json(pairing_));
    sum.set("not_processed", total ? std::max(0LL, total - done) : 0LL);
    d.set("summary", sum);
    Json hs = Json::array();
    for (auto& h : hints) {
        Json x = Json::object();
        x.set("code", "genetic_code");
        x.set("text", h.text);
        x.set("suggested_codontable", h.table);
        x.set("checked", h.checked);
        x.set("fit", h.fit);
        x.set("fit_current", h.fit_current);
        hs.push(x);
    }
    d.set("hints", hs);
    Json sq = Json::array();
    for (auto& e : entries_) sq.push(e.to_json());
    d.set("sequences", sq);
    double rt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_).count();
    d.set("runtime_seconds", (double)((long long)(rt * 1000 + 0.5)) / 1000.0);
    d.set("exit_code", exit_code);
    Json st = Json::array();
    for (auto& l : summary_lines()) st.push(l);
    d.set("summary_text", st);
    d.set("methods_text", methods_text());
    if (has_error) {
        Json e = Json::object();
        e.set("code", error_code);
        e.set("text", error_text);
        e.set("detail", error_detail);
        d.set("error", e);
    }
    if (resumed_from >= 0) d.set("resumed_from", resumed_from);
    return d;
}

Json RunReport::checkpoint_state() const {
    Json n = Json::object();
    n.set("ok", n_ok);
    n.set("warning", n_warning);
    n.set("failed", n_failed);
    n.set("mismatch", n_mismatch);
    n.set("unknown_aa", n_unknown);
    Json p = Json::array();
    for (long long k : problems) p.push(k);
    Json s = Json::object();
    s.set("n", n);
    s.set("problems", p);
    return s;
}

void RunReport::restore(const Json& state, const std::vector<Entry>& nonok, long long completed,
                        const std::vector<std::string>& prot_ids) {
    if (const Json* n = state.get("n")) {
        n_ok = n->get("ok") ? n->get("ok")->as_int() : 0;
        n_warning = n->get("warning") ? n->get("warning")->as_int() : 0;
        n_failed = n->get("failed") ? n->get("failed")->as_int() : 0;
        n_mismatch = n->get("mismatch") ? n->get("mismatch")->as_int() : 0;
        n_unknown = n->get("unknown_aa") ? n->get("unknown_aa")->as_int() : 0;
    }
    problems.clear();
    if (const Json* p = state.get("problems"))
        for (auto& v : p->a) problems.push_back(v.as_int());
    failed_entries_.clear();
    for (auto& e : nonok)
        if (e.status == "failed") failed_entries_.push_back(e);
    resumed_from = completed;
    if (keep_entries) {
        entries_.clear();
        size_t j = 0;
        for (long long k = 0; k < completed; ++k) {
            if (j < nonok.size() && nonok[j].index == k) {
                entries_.push_back(nonok[j++]);
                continue;
            }
            Entry e;
            e.id = (size_t)k < prot_ids.size() ? prot_ids[(size_t)k] : "";
            e.index = k;
            if (dna_ids_ && (size_t)k < dna_ids_->size()) {
                e.dna_id = (*dna_ids_)[(size_t)k];
                e.has_dna_id = true;
            }
            e.pairing = !seq_pairing_.empty() ? seq_pairing_[(size_t)k] : pairing_;
            e.status = "ok";
            entries_.push_back(e);
        }
    }
}

}  // namespace coal
