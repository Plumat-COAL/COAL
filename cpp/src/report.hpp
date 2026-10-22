// Run report (schema coal.report/1), the "coal: summary:" lines, the genetic-code hint and the methods text.
// Mirrors src/coal/report.py: the app's message codes and texts, statuses (ok / warning / failed / not_checked),
// where a failed gene breaks (the app viewer's rule), start-codon and frameshift-like mismatch patterns, the hint
// rule and wording, the methods text.
#pragma once
#include <chrono>
#include <functional>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "json.hpp"

namespace coal {

std::string fmt_n(long long x);  // 1234567 -> "1,234,567"
std::string shlex_join(const std::vector<std::string>& args);

// Rows of the protein alignment for the end-of-run extras: call fn(k, id, prepared row) for every sequence in input
// order until it returns false.
using RowsFn = std::function<void(const std::function<bool(long long, const std::string&, const std::string&)>&)>;
using DnaFn = std::function<std::string(long long)>;  // DNA paired with sequence k

// The genetic-code search shared by the end-of-run hint and the early notice
struct HintGene {
    std::string id, row, dna;
};
struct CodeSearch {
    std::vector<int> tied;  // best first (empty: no code does better)
    bool standard = false;  // code 1 is among the tied codes (and first)
    long long top = 0, fits_now = 0;
};
using FitMemo = std::map<std::pair<long long, int>, bool>;  // (sequence, code) -> fits
// The suggestion said before the per-codon warnings ("" = none): the first genes tested as the end-of-run hint
std::string early_code_notice(const RowsFn& rows, const DnaFn& dna, int table, long long total, FitMemo& memo);

class RunReport {
public:
    RunReport(const std::vector<std::string>& argv, const std::string& pep, const std::vector<std::string>& nucs,
              const std::string& out, const std::string& fmt, int table, bool nogap, bool nomismatch,
              bool keep_going, bool keep_entries, bool strict, bool blockonly);

    void set_hash_paths(std::vector<std::string> p) { hash_paths_ = std::move(p); }
    // The report owns the protein IDs (the pipeline reads them through ids()), one copy per run
    void set_inputs(std::vector<std::string>&& ids);
    const std::vector<std::string>& ids() const { return ids_; }
    const std::vector<std::string>& dna_file_ids() const { return dna_file_ids_; }
    void set_dna_inputs(const std::vector<long long>& dna_counts, const std::vector<std::string>& dna_file_ids);
    // mode id / order / normalised; per_seq: per-protein "id"/"normalised" (empty unless normalised);
    // pairing[k] = the DNA record (index into dna_file_ids) of protein k
    void set_pairing(const std::string& mode, const std::vector<std::string>& per_seq,
                     const std::vector<size_t>& pairing, const std::vector<std::string>& notes);
    // warnings: thread_one's ("pepAlnPos ..."); notices: with the ID; failure: code + message (code empty = none).
    // Returns the entry as one compact JSON line (json_dump(entry): the entries file's format), or nullptr for an ok
    // sequence.
    const std::string* add_sequence(long long k, const std::string& pid, const std::vector<std::string>& warnings,
                                    const std::vector<std::string>& notices, const std::string& fail_code = "",
                                    const std::string& fail_msg = "");
    void fail(const std::string& code, const std::string& message, bool has_seq_id);
    void set_output(const std::string& path, const std::string& sha, long long n);
    // End of run, inputs still open: break locations, mismatch patterns, genetic-code hint.
    void finish(const RowsFn& rows, const DnaFn& dna, bool aborted = false);
    std::vector<std::string> summary_lines() const;
    // The report (json.dump(report.to_dict(), indent=1) of report.py, without the final newline), streamed: one
    // sequence entry at a time, so a run of a million sequences does not hold a million entries at once.
    void write_json(std::ostream& out) const;
    // --resume
    Json checkpoint_state() const;
    void restore(const Json& state, const std::vector<Json>& nonok, long long completed);

    bool keep_entries;
    bool want_hint = true;  // false under -q without --report/--json (report.py want_hint)
    FitMemo fit_memo;  // shared with early_code_notice
    bool finished = false;
    int exit_code = 0;

private:
    long long total() const;
    std::string command() const;
    Json entry_view(long long k, int aborted = -1) const;  // -1: aborted = no output
    long long first_fail() const;
    // Statuses per kinds_ (views = false), or as the entry views of a run that stopped show them (views = true, the
    // report's summary: every sequence not_checked but the failed one)
    std::map<std::string, long long> counts(bool views = false) const;
    // A non-ok sequence's entry (parsed from its compact line) and the store that keeps it
    Json entry(long long k) const;
    void store(long long k, const Json& e);
    std::vector<Json> id_hints() const;
    std::vector<Json> hints() const;
    bool completed() const;
    std::string methods_text(const Json& rep) const;
    void finish_impl(const RowsFn& rows, const DnaFn& dna, bool aborted);
    void set_kind(long long k, char ch);

    std::chrono::steady_clock::time_point t0_;
    std::vector<std::string> argv_;
    std::string pep_, out_, fmt_;
    std::vector<std::string> nucs_, hash_paths_;
    int table_;
    bool nogap_, nomismatch_, keep_going_, strict_, blockonly_;
    bool has_ids_ = false;
    std::vector<std::string> ids_;
    bool has_dna_counts_ = false;
    std::vector<long long> dna_counts_;
    bool has_dna_file_ids_ = false;
    std::vector<std::string> dna_file_ids_;
    std::string pairing_;  // empty = not paired
    std::vector<std::string> seq_pairing_;
    bool has_dna_ids_ = false;
    std::vector<size_t> pair_idx_;  // protein k's DNA ID is dna_file_ids_[pair_idx_[k]]
    const std::string& dna_id(size_t k) const { return dna_file_ids_[pair_idx_[k]]; }
    mutable std::set<std::string> paired_set_;  // built on first use (suspect_pairing text)
    mutable bool paired_set_built_ = false;
    std::vector<Json> pair_hints_;
    bool partial_ = false;
    std::string kinds_;
    // Non-ok sequences only, as report.py keeps them (ok sequences are counted in kinds_), each as its compact
    // JSON line rather than a Json tree (about a tenth of the memory per sequence), plus what the summary lines and
    // the end-of-run step read from every entry without parsing it
    struct Kept {
        std::string json;
        bool failed_nocol = false;  // status failed with no column message (break_event candidates)
        bool real_mm = false;       // a definite (A/C/G/T) mismatch
        bool unknown_aa = false;    // an unknown_aa message
    };
    std::map<long long, Kept> entries_;
    std::vector<long long> cols_, real_cols_;  // sorted unique
    bool has_error_ = false;
    Json error_;
    std::vector<Json> errors_;
    bool has_output_ = false;
    Json output_;
    bool has_gc_hint_ = false;
    Json gc_hint_;
    std::string gc_hint_cli_;  // the hint as the command line says it
    long long resumed_from_ = -1;
};

}  // namespace coal
