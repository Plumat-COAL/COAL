// Run report (schema coal.report/1), end-of-run summary lines, genetic-code hint
// and methods text. Same content as src/coal/report.py.
#pragma once
#include <chrono>
#include <string>
#include <vector>

#include "json.hpp"

namespace coal {

struct Message {
    std::string code, text;
    int column = 0;  // 0 = none
};

struct Entry {
    std::string id;
    long long index = 0;
    std::string dna_id;
    bool has_dna_id = false;
    std::string pairing;
    std::string status;  // ok | warning | failed
    std::vector<Message> messages;
    Json to_json() const;
    static Entry from_json(const Json& j);
};

struct Hint {
    int table = 0, checked = 0, fit = 0, fit_current = 0;
    std::string text, summary;
};

std::string plain_reason(const std::string& code, const std::string& message);
std::string fmt_n(long long x);  // 1234567 -> "1,234,567"
std::string shlex_join(const std::vector<std::string>& args);
extern const char* const TABLE_NAMES[];  // index by table number (nullptr if unsupported)

class RunReport {
public:
    RunReport(const std::vector<std::string>& argv, const std::string& pep, const std::vector<std::string>& nucs,
              const std::string& out, const std::string& fmt, int table, bool nogap, bool nomismatch,
              bool keep_going, bool keep_entries);

    void set_inputs(long long n_protein, const std::vector<long long>& dna_counts);
    void set_pairing(const std::string& mode, const std::vector<std::string>* dna_ids,
                     const std::vector<size_t>* dna_index, const std::vector<std::string>& per_seq,
                     const std::string& ex_prot, const std::string& ex_dna);
    // Returns the entry when one was built (not for ok sequences without keep_entries).
    const Entry* add_sequence(const std::string& pid, const std::vector<std::string>& warnings,
                              const std::string& fail_code, const std::string& fail_msg, long long k);
    void fail(const std::string& code, const std::string& message);
    void set_output(const std::string& path, const std::string& sha, long long n);
    // Genetic-code hint: needs_hint() says whether to re-test the problem sequences;
    // apply_fits() takes (table, number of them that thread with no mismatch) for every table.
    bool needs_hint() const;
    void apply_fits(const std::vector<std::pair<int, int>>& fits);
    std::vector<std::string> summary_lines() const;
    std::string methods_text() const;
    Json to_json() const;
    // --resume support
    Json checkpoint_state() const;
    void restore(const Json& state, const std::vector<Entry>& nonok, long long completed,
                 const std::vector<std::string>& prot_ids);

    bool keep_entries;
    int exit_code = 0;
    std::vector<long long> problems;  // input positions of mismatching / failed sequences (<= 200)
    long long n_ok = 0, n_warning = 0, n_failed = 0, n_mismatch = 0, n_unknown = 0;
    long long n_protein = -1;         // -1 = unknown
    long long resumed_from = -1;
    bool has_error = false;
    bool pairing_set = false;
    std::string error_code, error_text, error_detail;
    std::vector<Hint> hints;

private:
    std::chrono::steady_clock::time_point t0_;
    std::string pep_, out_, fmt_, command_, pal2nal_;
    std::vector<std::string> nucs_;
    int table_;
    bool nogap_, nomismatch_, keep_going_;
    std::vector<long long> dna_counts_;
    bool has_dna_counts_ = false;
    std::string pairing_;
    const std::vector<std::string>* dna_ids_ = nullptr;
    std::vector<std::string> seq_pairing_;
    std::string ex_prot_, ex_dna_;
    std::vector<Entry> entries_, failed_entries_;
    bool has_output_ = false;
    std::string out_path_, out_sha_;
    long long out_n_ = 0;
};

}  // namespace coal
