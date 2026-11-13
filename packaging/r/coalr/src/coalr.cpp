// coalr: R bindings to COAL's C++ core (src/coal, vendored from cpp/src by scripts/sync_r_sources.py).
//
// Rules kept by every function here:
//   * R's API is used on the calling (main) thread only. Worker threads see std::string copies of the inputs and
//     write into pre-sized C++ vectors; the R objects of the result are built after the workers have been joined.
//   * No worker lets an exception escape: a pair that cannot be computed becomes a row with status "error".
//   * The main thread polls for a user interrupt (Ctrl-C / Esc) about every 100 ms while workers run, stops them
//     after the item each is working on, and only then signals the interrupt.
//   * Nothing is written to files, stdout or stderr, and no process is started.

// No floating-point contraction here either (as in coal/kaks.cpp, which says why): the quartiles of the cluster
// summary are interpolated with a * b + c, and R's Makeconf decides the flags of this file.
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#elif defined(_MSC_VER)
#pragma fp_contract(off)
#endif

#include <Rcpp.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "align.hpp"
#include "coal_version.h"
#include "core.hpp"
#include "dnds.hpp"
#include "kaks.hpp"
#include "tables.hpp"

using namespace coal;

namespace {

// ------------------------------------------------------------------------------------------------ threads

void check_interrupt_fn(void*) { R_CheckUserInterrupt(); }
// True when the user asked to stop. Main thread only. (R_ToplevelExec keeps R from jumping out of our frames.)
bool user_interrupt() { return R_ToplevelExec(check_interrupt_fn, nullptr) == FALSE; }

int resolve_threads(int threads) {
    if (threads > 0) return threads;
    const unsigned hc = std::thread::hardware_concurrency();
    return hc > 0 ? (int)hc : 1;
}

// fn(i) for every i in [0, n), on up to `threads` threads. fn must not use R and must not throw.
// Throws Rcpp's interrupt exception (after every worker has stopped) when the user interrupts.
template <class F>
void run_tasks(size_t n, int threads, F&& fn) {
    if (n == 0) return;
    using Clock = std::chrono::steady_clock;
    const size_t k = std::min<size_t>((size_t)std::max(1, threads), n);
    if (k <= 1) {
        auto last = Clock::now();
        for (size_t i = 0; i < n; ++i) {
            fn(i);
            const auto now = Clock::now();
            if (now - last > std::chrono::milliseconds(100)) {
                last = now;
                if (user_interrupt()) throw Rcpp::internal::InterruptedException();
            }
        }
        return;
    }
    std::atomic<size_t> next{0};
    std::atomic<bool> stop{false};
    std::mutex m;
    std::condition_variable cv;
    size_t finished = 0;
    auto body = [&]() {
        size_t i;
        while (!stop.load(std::memory_order_relaxed) && (i = next.fetch_add(1)) < n) fn(i);
        {
            std::lock_guard<std::mutex> g(m);
            ++finished;
        }
        cv.notify_one();
    };
    std::vector<std::thread> workers;
    workers.reserve(k);
    try {
        for (size_t t = 0; t < k; ++t) workers.emplace_back(body);
    } catch (...) {
        // could not start a thread: the ones that run finish the work; none is left behind
    }
    if (workers.empty()) {  // no thread at all: do the work here
        body();
        return;
    }
    bool interrupted = false;
    {
        std::unique_lock<std::mutex> lk(m);
        while (finished < workers.size()) {
            cv.wait_for(lk, std::chrono::milliseconds(100));
            if (finished >= workers.size()) break;
            if (!interrupted) {
                lk.unlock();
                const bool intr = user_interrupt();
                lk.lock();
                if (intr) {
                    interrupted = true;
                    stop.store(true);
                }
            }
        }
    }
    for (auto& th : workers) th.join();
    if (interrupted) throw Rcpp::internal::InterruptedException();
}

// ------------------------------------------------------------------------------------------------ helpers

inline bool is_na(SEXP x, R_xlen_t i) { return STRING_ELT(x, i) == NA_STRING; }
inline std::string str_at(SEXP x, R_xlen_t i) {
    SEXP s = STRING_ELT(x, i);
    return std::string(CHAR(s), (size_t)LENGTH(s));
}
inline double na_real(double v) { return std::isnan(v) ? NA_REAL : v; }

Rcpp::CharacterVector to_r(const std::vector<std::string>& v) {
    Rcpp::CharacterVector out(v.size());
    for (size_t i = 0; i < v.size(); ++i) SET_STRING_ELT(out, (R_xlen_t)i, Rf_mkCharLen(v[i].data(), (int)v[i].size()));
    return out;
}
Rcpp::NumericVector to_r(const std::vector<double>& v) {
    Rcpp::NumericVector out(v.size());
    for (size_t i = 0; i < v.size(); ++i) out[(R_xlen_t)i] = na_real(v[i]);
    return out;
}
Rcpp::IntegerVector to_r(const std::vector<int>& v) { return Rcpp::IntegerVector(v.begin(), v.end()); }

const GeneticCode* need_code(int table) {
    const GeneticCode* code = get_code(table);
    if (!code) {
        std::string ok;
        for (int t : supported_tables()) ok += (ok.empty() ? "" : ", ") + std::to_string(t);
        Rcpp::stop("genetic code %d is not supported (NCBI tables: %s)", table, ok.c_str());
    }
    return code;
}

std::vector<kaks::Method> need_methods(const Rcpp::CharacterVector& methods) {
    std::vector<kaks::Method> ms;
    for (R_xlen_t i = 0; i < methods.size(); ++i) {
        if (is_na(methods, i)) Rcpp::stop("a method is NA");
        kaks::Method m;
        std::string why;
        if (!kaks::parse_method(str_at(methods, i), m, why)) Rcpp::stop("%s", why.c_str());
        ms.push_back(m);
    }
    if (ms.empty()) Rcpp::stop("no Ka/Ks method given");
    return dnds::ordered_methods(ms);
}

AlignParams need_align(const Rcpp::List& a) {
    AlignParams p;
    p.matrix = Rcpp::as<std::string>(a["matrix"]);
    p.gap_open = Rcpp::as<double>(a["gap_open"]);
    p.gap_extend = Rcpp::as<double>(a["gap_extend"]);
    p.free_end_gaps = Rcpp::as<bool>(a["free_end_gaps"]);
    const double mc = Rcpp::as<double>(a["max_cells"]);
    if (!(mc >= 0) || !std::isfinite(mc)) Rcpp::stop("max_cells must be a number >= 0 (0 = no limit)");
    p.max_cells = (uint64_t)mc;
    try {
        check_align_params(p);
    } catch (const std::exception& e) {
        Rcpp::stop("%s", e.what());
    }
    return p;
}

kaks::Options need_kaks(const std::string& flavour, bool weighting, const std::string& comeron_flavour,
                        const std::string& li_flavour) {
    kaks::Options o;
    if (comeron_flavour == "gestimator")
        o.comeron_flavour = kaks::ComeronFlavour::gestimator;
    else if (comeron_flavour != "coal")
        Rcpp::stop("comeron_flavour must be \"coal\" or \"gestimator\"");
    if (li_flavour == "seqinr")
        o.li_flavour = kaks::LiFlavour::seqinr;
    else if (li_flavour != "coal")
        Rcpp::stop("li_flavour must be \"coal\" or \"seqinr\"");
    if (flavour == "paml")
        o.yn_flavour = kaks::Flavour::paml;
    else if (flavour != "kaks")
        Rcpp::stop("yn_flavour must be \"kaks\" or \"paml\"");
    o.yn_weighting = weighting;
    return o;
}

// The columns of the pair table (one entry per pair and method).
struct PairColumns {
    std::vector<int> unit, m1, m2, method, status, codons, plen;
    std::vector<double> dn, ds, w, se_dn, se_ds, N, S, Nd, Sd, ident;
    std::vector<std::string> note;
    void resize(size_t n) {
        for (auto* v : {&unit, &m1, &m2, &method, &status, &codons, &plen}) v->resize(n, NA_INTEGER);
        for (auto* v : {&dn, &ds, &w, &se_dn, &se_ds, &N, &S, &Nd, &Sd, &ident}) v->resize(n, kaks::NA);
        note.resize(n);
    }
    size_t size() const { return unit.size(); }
};

struct Item {
    dnds::Group g;
    dnds::GroupResult r;
    std::vector<std::pair<int, int>> pairs;
    size_t unit = 0;    // index of the unit
    size_t mstart = 0;  // position of member 0 in `members`
};

struct Task {
    Item* item;
    int i, j;
    size_t row;  // first of the pair's rows
};

}  // namespace

// ================================================================================================ dN/dS

// Units are groups of coding sequences. `members` holds, unit after unit, 0-based indices into `cds` (NA: no such
// record); `unit_start[u] .. unit_start[u + 1]` are unit u's members. `rows` (or NULL) holds each member's row of
// the unit's protein alignment ("" or NA: none). `unit_pairs`: 0 = every pair, 1 = member 0 against the others,
// 2 = the pairs pl_i/pl_j[pl_start[u] .. pl_start[u + 1]] (member numbers within the unit, from 0).
// [[Rcpp::export(rng = false)]]
Rcpp::List coalr_dnds(SEXP cds, SEXP member_ids, Rcpp::IntegerVector members, SEXP rows,
                      Rcpp::IntegerVector unit_start, Rcpp::IntegerVector unit_pairs, Rcpp::IntegerVector pl_start,
                      Rcpp::IntegerVector pl_i, Rcpp::IntegerVector pl_j, Rcpp::CharacterVector methods, int table,
                      int threads, Rcpp::List align, std::string yn_flavour, bool yn_weighting,
                      std::string comeron_flavour, std::string li_flavour, bool compare_terminal_stops,
                      bool keep_alignments) {
    if (TYPEOF(cds) != STRSXP || TYPEOF(member_ids) != STRSXP) Rcpp::stop("internal: cds and ids must be character");
    const bool have_rows = !Rf_isNull(rows);
    if (have_rows && (TYPEOF(rows) != STRSXP || XLENGTH(rows) != members.size()))
        Rcpp::stop("internal: rows must be a character vector parallel to members");
    if (XLENGTH(member_ids) != members.size()) Rcpp::stop("internal: ids must be parallel to members");
    const size_t n_units = (size_t)unit_pairs.size();
    if ((size_t)unit_start.size() != n_units + 1 || (size_t)pl_start.size() != n_units + 1)
        Rcpp::stop("internal: unit offsets have the wrong length");
    const R_xlen_t n_cds = XLENGTH(cds);

    dnds::Options opt;
    opt.methods = need_methods(methods);
    opt.table = table;
    need_code(table);
    opt.align = need_align(align);
    opt.kaks = need_kaks(yn_flavour, yn_weighting, comeron_flavour, li_flavour);
    opt.compare_terminal_stops = compare_terminal_stops;
    const std::vector<kaks::Method> ms = opt.methods;
    const size_t nm = ms.size();
    const int nthreads = resolve_threads(threads);

    PairColumns col;
    std::vector<std::string> unit_aligner(n_units), unit_error(n_units), unit_error_code(n_units);
    std::vector<int> fail_member;
    std::vector<std::string> fail_code, fail_message;
    std::vector<std::string> codon_rows, protein_rows;
    if (keep_alignments) {
        codon_rows.resize((size_t)members.size());
        protein_rows.resize((size_t)members.size());
    }
    using Clock = std::chrono::steady_clock;
    double t_prepare = 0, t_kaks = 0;
    auto since = [](Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); };

    const size_t BATCH_GROUPS = 1024, BATCH_NT = 16u << 20, BATCH_PAIRS = 32768;
    size_t u = 0;
    std::vector<Item> batch;
    std::vector<Task> tasks;
    while (u < n_units) {
        // ---- build one batch on this thread (the only place the inputs are read)
        batch.clear();
        size_t nt = 0, np = 0;
        while (u < n_units && batch.size() < BATCH_GROUPS && nt < BATCH_NT && np < BATCH_PAIRS) {
            Item it;
            it.unit = u;
            const int a = unit_start[(R_xlen_t)u], b = unit_start[(R_xlen_t)u + 1];
            if (a < 0 || b < a || b > members.size()) Rcpp::stop("internal: bad unit offsets");
            it.mstart = (size_t)a;
            it.g.members.resize((size_t)(b - a));
            for (int k = a; k < b; ++k) {
                dnds::Member& mem = it.g.members[(size_t)(k - a)];
                if (!is_na(member_ids, k)) mem.id = str_at(member_ids, k);
                const int idx = members[k];
                if (idx == NA_INTEGER || idx < 0 || idx >= n_cds) {
                    mem.failed = "no CDS record with this ID";
                } else if (is_na(cds, idx)) {
                    mem.failed = "the CDS is NA";
                } else {
                    mem.cds = str_at(cds, idx);
                    nt += mem.cds.size();
                    if (have_rows && !is_na(rows, k)) mem.row = str_at(rows, k);
                }
            }
            const int mode = unit_pairs[(R_xlen_t)u];
            if (mode == 1) {
                it.g.pairs = dnds::Pairs::representative;
            } else if (mode == 2) {
                it.g.pairs = dnds::Pairs::list;
                const int pa = pl_start[(R_xlen_t)u], pb = pl_start[(R_xlen_t)u + 1];
                if (pa < 0 || pb < pa || pb > pl_i.size() || pb > pl_j.size())
                    Rcpp::stop("internal: bad pair-list offsets");
                it.g.pair_list.reserve((size_t)(pb - pa));
                for (int k = pa; k < pb; ++k) it.g.pair_list.emplace_back(pl_i[k], pl_j[k]);
            }
            it.pairs = dnds::group_pairs(it.g);
            np += it.pairs.size();
            batch.push_back(std::move(it));
            ++u;
        }
        // ---- translate, align, thread the codons
        {
            const auto t0 = Clock::now();
            run_tasks(batch.size(), nthreads, [&](size_t k) {
                Item& it = batch[k];
                try {
                    it.r = dnds::prepare_group(it.g, opt);
                } catch (const std::bad_alloc&) {
                    it.r = dnds::GroupResult();
                    it.r.error = "out of memory";
                    it.r.error_code = "out_of_memory";
                } catch (const std::exception& e) {
                    it.r = dnds::GroupResult();
                    it.r.error = e.what();
                    it.r.error_code = "error";
                } catch (...) {
                    it.r = dnds::GroupResult();
                    it.r.error = "unknown error";
                    it.r.error_code = "error";
                }
                if (it.r.members.size() != it.g.members.size()) {
                    it.r.members.resize(it.g.members.size());
                    it.r.protein_rows.assign(it.g.members.size(), "");
                }
            });
            t_prepare += since(t0);
        }
        // ---- Ka/Ks for every pair of the batch
        tasks.clear();
        size_t row = col.size();
        for (Item& it : batch)
            for (auto& ij : it.pairs) {
                tasks.push_back(Task{&it, ij.first, ij.second, row});
                row += nm;
            }
        col.resize(row);
        {
            const auto t0 = Clock::now();
            run_tasks(tasks.size(), nthreads, [&](size_t k) {
                const Task& t = tasks[k];
                dnds::PairResult p;
                try {
                    p = dnds::compute_pair(t.item->r, t.i, t.j, opt);
                    if (p.error.empty() && p.results.size() != nm) p.error = "internal: wrong number of results";
                } catch (const std::bad_alloc&) {
                    p.error = "out of memory";
                } catch (const std::exception& e) {
                    p.error = e.what();
                } catch (...) {
                    p.error = "unknown error";
                }
                for (size_t q = 0; q < nm; ++q) {
                    const size_t r = t.row + q;
                    col.unit[r] = (int)t.item->unit;
                    col.m1[r] = (int)(t.item->mstart + (size_t)t.i);
                    col.m2[r] = (int)(t.item->mstart + (size_t)t.j);
                    col.method[r] = (int)ms[q];
                    if (!p.error.empty()) {
                        col.status[r] = (int)kaks::Status::error;
                        col.codons[r] = 0;
                        col.note[r] = p.error;
                        continue;
                    }
                    const kaks::Result& res = p.results[q];
                    kaks::Status st = res.status;
                    std::string note = res.note;
                    if (!p.note.empty()) {
                        note += (note.empty() ? "" : "; ") + p.note;
                        if (st == kaks::Status::ok) st = kaks::Status::warning;
                    }
                    col.method[r] = (int)res.method;
                    col.status[r] = (int)st;
                    col.note[r] = std::move(note);
                    col.dn[r] = res.ka;
                    col.ds[r] = res.ks;
                    col.w[r] = res.ka_ks;
                    col.se_dn[r] = res.se_ka;
                    col.se_ds[r] = res.se_ks;
                    col.N[r] = res.N;
                    col.S[r] = res.S;
                    col.Nd[r] = res.Nd;
                    col.Sd[r] = res.Sd;
                    col.codons[r] = (int)res.codons_compared;
                    col.plen[r] = (int)p.prot_aln_length;
                    // two decimals, rounded as the command line's "%.2f" rounds them
                    char ident[32];
                    std::snprintf(ident, sizeof ident, "%.2f",
                                  p.prot_aln_length
                                      ? 100.0 * (double)p.prot_identities / (double)p.prot_aln_length
                                      : 0.0);
                    col.ident[r] = std::strtod(ident, nullptr);
                }
            });
            t_kaks += since(t0);
        }
        // ---- what failed, and the alignments when asked for
        for (Item& it : batch) {
            unit_aligner[it.unit] = it.r.aligner;
            unit_error[it.unit] = it.r.error;
            unit_error_code[it.unit] = it.r.error_code;
            for (size_t k = 0; k < it.g.members.size(); ++k) {
                const dnds::MemberResult& m = it.r.members[k];
                if (!m.error.empty()) {
                    fail_member.push_back((int)(it.mstart + k));
                    fail_code.push_back(m.error_code);
                    fail_message.push_back(m.error);
                }
                if (keep_alignments) {
                    codon_rows[it.mstart + k] = m.codon_row;
                    if (k < it.r.protein_rows.size()) protein_rows[it.mstart + k] = it.r.protein_rows[k];
                }
            }
        }
        batch.clear();
        Rcpp::checkUserInterrupt();  // no worker is running here
    }

    std::vector<std::string> method_names, status_names;
    for (kaks::Method m : kaks::output_order()) method_names.push_back(kaks::method_name(m));
    for (int s = 0; s < kaks::STATUS_COUNT; ++s)   // the codes index status_names
        status_names.push_back(kaks::status_name((kaks::Status)s));
    // Method's enumerators are in output order (NG, LWL, MLWL, LPB, MLPB, YN, MYN, Comeron, Li): the codes index
    // method_names.
    std::vector<std::string> enum_names(method_names.size());
    for (kaks::Method m : kaks::output_order()) enum_names[(size_t)(int)m] = kaks::method_name(m);

    Rcpp::List pairs = Rcpp::List::create(
        Rcpp::_["unit"] = to_r(col.unit), Rcpp::_["m1"] = to_r(col.m1), Rcpp::_["m2"] = to_r(col.m2),
        Rcpp::_["method"] = to_r(col.method), Rcpp::_["dN"] = to_r(col.dn), Rcpp::_["dS"] = to_r(col.ds),
        Rcpp::_["dNdS"] = to_r(col.w), Rcpp::_["SE_dN"] = to_r(col.se_dn), Rcpp::_["SE_dS"] = to_r(col.se_ds),
        Rcpp::_["N_sites"] = to_r(col.N), Rcpp::_["S_sites"] = to_r(col.S), Rcpp::_["Nd"] = to_r(col.Nd),
        Rcpp::_["Sd"] = to_r(col.Sd), Rcpp::_["codons_compared"] = to_r(col.codons),
        Rcpp::_["prot_aln_length"] = to_r(col.plen), Rcpp::_["prot_identity"] = to_r(col.ident),
        Rcpp::_["status"] = to_r(col.status), Rcpp::_["note"] = to_r(col.note));
    Rcpp::List out = Rcpp::List::create(
        Rcpp::_["pairs"] = pairs, Rcpp::_["method_names"] = to_r(enum_names),
        Rcpp::_["status_names"] = to_r(status_names), Rcpp::_["unit_aligner"] = to_r(unit_aligner),
        Rcpp::_["unit_error"] = to_r(unit_error), Rcpp::_["unit_error_code"] = to_r(unit_error_code),
        Rcpp::_["fail_member"] = to_r(fail_member), Rcpp::_["fail_code"] = to_r(fail_code),
        Rcpp::_["fail_message"] = to_r(fail_message),
        Rcpp::_["seconds"] = Rcpp::NumericVector::create(Rcpp::_["align_codon"] = t_prepare, Rcpp::_["kaks"] = t_kaks),
        Rcpp::_["threads"] = nthreads);
    if (keep_alignments) {
        out["codon_rows"] = to_r(codon_rows);
        out["protein_rows"] = to_r(protein_rows);
    }
    return out;
}

// ================================================================================================ summary

namespace {

// Linear interpolation between order statistics (R's type 7); v sorted, not empty.
double quantile7(const std::vector<double>& v, double q) {
    const double h = (double)(v.size() - 1) * q;
    const size_t lo = (size_t)std::floor(h);
    const size_t hi = std::min(lo + 1, v.size() - 1);
    return v[lo] + (v[hi] - v[lo]) * (h - (double)lo);
}

}  // namespace

// Per cluster and method: counts by status and median / mean / quartiles / min / max of dN, dS, dN/dS over the
// pairs with status ok or warning (the ones with an estimate). `cluster` (0-based, non-decreasing) and the other vectors are the pair table's
// columns; n_methods is the number of method codes. One row per (cluster, method) that has pairs, in order.
// [[Rcpp::export(rng = false)]]
Rcpp::List coalr_summary(Rcpp::IntegerVector cluster, Rcpp::IntegerVector method, Rcpp::IntegerVector status,
                         Rcpp::NumericVector dn, Rcpp::NumericVector ds, Rcpp::NumericVector w, int n_methods) {
    const R_xlen_t n = cluster.size();
    std::vector<int> o_cluster, o_method, o_n, o_ok, o_nc, o_sat, o_nodata, o_failed, o_nw;
    std::vector<std::vector<double>> stats(18);
    struct Acc {
        int n = 0, st[kaks::STATUS_COUNT] = {};
        std::vector<double> v[3];
    };
    R_xlen_t i = 0;
    while (i < n) {
        const int c = cluster[i];
        std::vector<Acc> acc((size_t)n_methods);
        R_xlen_t j = i;
        for (; j < n && cluster[j] == c; ++j) {
            const int m = method[j], s = status[j];
            if (m < 0 || m >= n_methods || s < 0 || s >= kaks::STATUS_COUNT) Rcpp::stop("internal: bad method or status code");
            Acc& a = acc[(size_t)m];
            ++a.n;
            ++a.st[s];
            if (s == (int)kaks::Status::ok || s == (int)kaks::Status::warning) {
                if (!ISNAN(dn[j])) a.v[0].push_back(dn[j]);
                if (!ISNAN(ds[j])) a.v[1].push_back(ds[j]);
                if (!ISNAN(w[j])) a.v[2].push_back(w[j]);
            }
        }
        for (int m = 0; m < n_methods; ++m) {
            Acc& a = acc[(size_t)m];
            if (!a.n) continue;
            o_cluster.push_back(c);
            o_method.push_back(m);
            o_n.push_back(a.n);
            o_ok.push_back(a.st[(int)kaks::Status::ok] + a.st[(int)kaks::Status::warning]);
            o_nc.push_back(a.st[(int)kaks::Status::not_converged]);
            o_sat.push_back(a.st[(int)kaks::Status::saturated]);
            o_nodata.push_back(a.st[(int)kaks::Status::no_data]);
            o_failed.push_back(a.st[(int)kaks::Status::error]);
            o_nw.push_back((int)a.v[2].size());
            for (int k = 0; k < 3; ++k) {
                std::vector<double>& v = a.v[k];
                double st6[6] = {kaks::NA, kaks::NA, kaks::NA, kaks::NA, kaks::NA, kaks::NA};
                if (!v.empty()) {
                    std::sort(v.begin(), v.end());
                    long double sum = 0;
                    for (double x : v) sum += x;
                    st6[0] = quantile7(v, 0.5);
                    st6[1] = (double)(sum / (long double)v.size());
                    st6[2] = quantile7(v, 0.25);
                    st6[3] = quantile7(v, 0.75);
                    st6[4] = v.front();
                    st6[5] = v.back();
                }
                for (int q = 0; q < 6; ++q) stats[(size_t)(k * 6 + q)].push_back(st6[q]);
            }
        }
        i = j;
        if ((o_cluster.size() & 4095) == 0) Rcpp::checkUserInterrupt();
    }
    Rcpp::List out = Rcpp::List::create(
        Rcpp::_["cluster"] = to_r(o_cluster), Rcpp::_["method"] = to_r(o_method), Rcpp::_["n_pairs"] = to_r(o_n),
        Rcpp::_["n_ok"] = to_r(o_ok), Rcpp::_["n_not_converged"] = to_r(o_nc), Rcpp::_["n_saturated"] = to_r(o_sat),
        Rcpp::_["n_no_data"] = to_r(o_nodata), Rcpp::_["n_failed"] = to_r(o_failed));
    const char* vars[3] = {"dN", "dS", "dNdS"};
    const char* names[6] = {"median", "mean", "q25", "q75", "min", "max"};
    for (int k = 0; k < 3; ++k)
        for (int q = 0; q < 6; ++q)
            out[std::string(vars[k]) + "_" + names[q]] = to_r(stats[(size_t)(k * 6 + q)]);
    out["n_dNdS"] = to_r(o_nw);
    return out;
}

// ================================================================================================ Ka/Ks

// Ka/Ks for pairs of codon-aligned sequences. One row per pair and method.
// [[Rcpp::export(rng = false)]]
Rcpp::List coalr_kaks(SEXP a, SEXP b, Rcpp::CharacterVector methods, int table, int threads, std::string yn_flavour,
                      bool yn_weighting, std::string comeron_flavour, std::string li_flavour) {
    if (TYPEOF(a) != STRSXP || TYPEOF(b) != STRSXP || XLENGTH(a) != XLENGTH(b))
        Rcpp::stop("internal: the two alignments must be character vectors of one length");
    const std::vector<kaks::Method> ms = need_methods(methods);
    const GeneticCode* code = need_code(table);
    const kaks::Options ko = need_kaks(yn_flavour, yn_weighting, comeron_flavour, li_flavour);
    const size_t n = (size_t)XLENGTH(a), nm = ms.size();
    const int nthreads = resolve_threads(threads);
    const size_t total = n * nm;

    std::vector<int> pair(total), method(total), status(total, (int)kaks::Status::error), codons(total, 0),
        length(total, 0), stops(total, 0), gaps(total, 0), amb(total, 0);
    std::vector<double> dn(total, kaks::NA), ds(total, kaks::NA), w(total, kaks::NA), se_dn(total, kaks::NA),
        se_ds(total, kaks::NA), N(total, kaks::NA), S(total, kaks::NA), Nd(total, kaks::NA), Sd(total, kaks::NA),
        pf(total, kaks::NA), kappa(total, kaks::NA), div(total, kaks::NA), subs(total, kaks::NA), gc(total, kaks::NA);
    std::vector<std::string> note(total);

    const size_t BATCH = 4096, BATCH_NT = 32u << 20;
    std::vector<std::pair<std::string, std::string>> rows;
    std::vector<char> missing;
    size_t i0 = 0;
    while (i0 < n) {
        rows.clear();
        missing.clear();
        size_t nt = 0;
        while (i0 + rows.size() < n && rows.size() < BATCH && nt < BATCH_NT) {
            const R_xlen_t i = (R_xlen_t)(i0 + rows.size());
            if (is_na(a, i) || is_na(b, i)) {
                rows.emplace_back();
                missing.push_back(1);
            } else {
                rows.emplace_back(str_at(a, i), str_at(b, i));
                missing.push_back(0);
                nt += rows.back().first.size() * 2;
            }
        }
        run_tasks(rows.size(), nthreads, [&](size_t k) {
            const size_t base = (i0 + k) * nm;
            for (size_t q = 0; q < nm; ++q) {
                pair[base + q] = (int)(i0 + k);
                method[base + q] = (int)ms[q];
            }
            if (missing[k]) {
                for (size_t q = 0; q < nm; ++q) note[base + q] = "a sequence is NA";
                return;
            }
            try {
                const kaks::Prepared prep = kaks::prepare(rows[k].first, rows[k].second, *code);
                for (size_t q = 0; q < nm; ++q) {
                    const kaks::Result r = kaks::compute(prep, ms[q], *code, ko);
                    const size_t x = base + q;
                    status[x] = (int)r.status;
                    note[x] = r.note;
                    dn[x] = r.ka;
                    ds[x] = r.ks;
                    w[x] = r.ka_ks;
                    se_dn[x] = r.se_ka;
                    se_ds[x] = r.se_ks;
                    N[x] = r.N;
                    S[x] = r.S;
                    Nd[x] = r.Nd;
                    Sd[x] = r.Sd;
                    pf[x] = r.p_fisher;
                    kappa[x] = r.kappa;
                    div[x] = r.divergence;
                    subs[x] = r.substitutions;
                    gc[x] = r.gc[0];
                    codons[x] = (int)r.codons_compared;
                    length[x] = (int)r.length;
                    stops[x] = (int)r.internal_stops;
                    gaps[x] = (int)r.gap_codons;
                    amb[x] = (int)r.ambiguous_codons;
                }
            } catch (const std::bad_alloc&) {
                for (size_t q = 0; q < nm; ++q) note[base + q] = "out of memory";
            } catch (const std::exception& e) {
                for (size_t q = 0; q < nm; ++q) note[base + q] = e.what();
            } catch (...) {
                for (size_t q = 0; q < nm; ++q) note[base + q] = "unknown error";
            }
        });
        i0 += rows.size();
        Rcpp::checkUserInterrupt();
    }
    std::vector<std::string> enum_names(kaks::output_order().size()), status_names;
    for (kaks::Method m : kaks::output_order()) enum_names[(size_t)(int)m] = kaks::method_name(m);
    for (int s = 0; s < kaks::STATUS_COUNT; ++s)   // the codes index status_names
        status_names.push_back(kaks::status_name((kaks::Status)s));
    return Rcpp::List::create(
        Rcpp::_["pair"] = to_r(pair), Rcpp::_["method"] = to_r(method), Rcpp::_["dN"] = to_r(dn),
        Rcpp::_["dS"] = to_r(ds), Rcpp::_["dNdS"] = to_r(w), Rcpp::_["SE_dN"] = to_r(se_dn),
        Rcpp::_["SE_dS"] = to_r(se_ds), Rcpp::_["N_sites"] = to_r(N), Rcpp::_["S_sites"] = to_r(S),
        Rcpp::_["Nd"] = to_r(Nd), Rcpp::_["Sd"] = to_r(Sd), Rcpp::_["codons_compared"] = to_r(codons),
        Rcpp::_["status"] = to_r(status), Rcpp::_["note"] = to_r(note), Rcpp::_["p_fisher"] = to_r(pf),
        Rcpp::_["kappa"] = to_r(kappa), Rcpp::_["divergence"] = to_r(div), Rcpp::_["substitutions"] = to_r(subs),
        Rcpp::_["length"] = to_r(length), Rcpp::_["gc"] = to_r(gc), Rcpp::_["internal_stops"] = to_r(stops),
        Rcpp::_["gap_codons"] = to_r(gaps), Rcpp::_["ambiguous_codons"] = to_r(amb),
        Rcpp::_["method_names"] = to_r(enum_names), Rcpp::_["status_names"] = to_r(status_names));
}

// ================================================================================================ alignment

// Global protein alignment of a[i] against b[i].
// [[Rcpp::export(rng = false)]]
Rcpp::List coalr_align(SEXP a, SEXP b, Rcpp::List align, int threads) {
    if (TYPEOF(a) != STRSXP || TYPEOF(b) != STRSXP || XLENGTH(a) != XLENGTH(b))
        Rcpp::stop("internal: a and b must be character vectors of one length");
    const AlignParams p = need_align(align);
    const size_t n = (size_t)XLENGTH(a);
    const int nthreads = resolve_threads(threads);
    std::vector<std::string> ra(n), rb(n), err(n);
    std::vector<double> score(n, kaks::NA);
    std::vector<int> length(n, NA_INTEGER), ident(n, NA_INTEGER), simil(n, NA_INTEGER), gaps(n, NA_INTEGER);
    std::vector<char> failed(n, 0);

    const size_t BATCH = 1024, BATCH_AA = 8u << 20;
    std::vector<std::pair<std::string, std::string>> in;
    std::vector<char> missing;
    size_t i0 = 0;
    while (i0 < n) {
        in.clear();
        missing.clear();
        size_t aa = 0;
        while (i0 + in.size() < n && in.size() < BATCH && aa < BATCH_AA) {
            const R_xlen_t i = (R_xlen_t)(i0 + in.size());
            if (is_na(a, i) || is_na(b, i)) {
                in.emplace_back();
                missing.push_back(1);
            } else {
                in.emplace_back(str_at(a, i), str_at(b, i));
                missing.push_back(0);
                aa += in.back().first.size() + in.back().second.size();
            }
        }
        run_tasks(in.size(), nthreads, [&](size_t k) {
            const size_t x = i0 + k;
            if (missing[k]) {
                failed[x] = 1;
                err[x] = "a sequence is NA";
                return;
            }
            try {
                PairAlignment r = align_pair(in[k].first, in[k].second, p);
                ra[x] = std::move(r.a);
                rb[x] = std::move(r.b);
                score[x] = r.score;
                length[x] = (int)r.length;
                ident[x] = (int)r.identities;
                simil[x] = (int)r.similarities;
                gaps[x] = (int)r.gaps;
            } catch (const std::bad_alloc&) {
                failed[x] = 1;
                err[x] = "out of memory";
            } catch (const std::exception& e) {
                failed[x] = 1;
                err[x] = e.what();
            } catch (...) {
                failed[x] = 1;
                err[x] = "unknown error";
            }
        });
        i0 += in.size();
        Rcpp::checkUserInterrupt();
    }
    Rcpp::CharacterVector oa = to_r(ra), ob = to_r(rb), oe = to_r(err);
    for (size_t i = 0; i < n; ++i) {
        if (failed[i]) {
            SET_STRING_ELT(oa, (R_xlen_t)i, NA_STRING);
            SET_STRING_ELT(ob, (R_xlen_t)i, NA_STRING);
        } else {
            SET_STRING_ELT(oe, (R_xlen_t)i, NA_STRING);
        }
    }
    return Rcpp::List::create(Rcpp::_["aligned_a"] = oa, Rcpp::_["aligned_b"] = ob, Rcpp::_["score"] = to_r(score),
                              Rcpp::_["length"] = to_r(length), Rcpp::_["identities"] = to_r(ident),
                              Rcpp::_["similarities"] = to_r(simil), Rcpp::_["gaps"] = to_r(gaps),
                              Rcpp::_["error"] = oe);
}

// ================================================================================================ codon alignment

// Protein alignment rows + their coding DNA -> codon alignment rows (what PAL2NAL does). One row per sequence;
// NA where a sequence could not be threaded (its message is in `error`).
// [[Rcpp::export(rng = false)]]
Rcpp::List coalr_codon_align(SEXP rows, SEXP cds, SEXP ids, int table, bool strict, bool nogap, bool nomismatch) {
    if (TYPEOF(rows) != STRSXP || TYPEOF(cds) != STRSXP || TYPEOF(ids) != STRSXP || XLENGTH(rows) != XLENGTH(cds) ||
        XLENGTH(rows) != XLENGTH(ids))
        Rcpp::stop("internal: rows, cds and ids must be character vectors of one length");
    const GeneticCode* code = need_code(table);
    const size_t n = (size_t)XLENGTH(rows);
    std::vector<std::string> prow(n), out(n), err(n), warn(n);
    std::vector<char> failed(n, 0);
    std::vector<std::vector<int>> mcols(n);
    std::vector<int> widths;
    long long ncol = -1;
    size_t first = 0;
    for (size_t i = 0; i < n; ++i) {
        if (is_na(rows, (R_xlen_t)i) || is_na(cds, (R_xlen_t)i)) {
            failed[i] = 1;
            err[i] = "the protein row or the CDS is NA";
            continue;
        }
        std::string raw = str_at(rows, (R_xlen_t)i), stripped;
        stripped.reserve(raw.size());
        for (unsigned char ch : raw)
            if (!(ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r')) stripped += (char)ch;
        prow[i] = prepare_row(stripped);
        if (ncol < 0) {
            ncol = (long long)prow[i].size();
            widths.assign((size_t)ncol, 3);
            first = i;
        } else if ((long long)prow[i].size() != ncol) {
            Rcpp::stop("the rows of the protein alignment differ in length (%s: %d, %s: %d)",
                       CHAR(STRING_ELT(ids, (R_xlen_t)first)), (int)ncol, CHAR(STRING_ELT(ids, (R_xlen_t)i)),
                       (int)prow[i].size());
        }
        update_column_widths(widths, prow[i]);
    }
    long long total = 0;
    for (int w : widths) total += w;
    std::vector<char> drop((size_t)(ncol > 0 ? ncol : 0), 0), bad((size_t)(total / 3), 0);
    for (size_t i = 0; i < n; ++i) {
        if (failed[i]) continue;
        const std::string id = is_na(ids, (R_xlen_t)i) ? std::string("sequence ") + std::to_string(i + 1)
                                                        : str_at(ids, (R_xlen_t)i);
        std::string dna;
        {
            const std::string raw = str_at(cds, (R_xlen_t)i);
            dna.reserve(raw.size());
            for (unsigned char ch : raw)
                if ((ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z')) dna += (char)ch;  // as PAL2NAL: letters only
        }
        try {
            Threaded t = thread_one(id, prow[i], dna, *code, strict);
            out[i] = expand(prow[i], t.dna, widths);
            for (auto& s : t.warnings) warn[i] += (warn[i].empty() ? "" : "; ") + s;
            for (auto& s : t.notices) warn[i] += (warn[i].empty() ? "" : "; ") + s;
            mcols[i] = std::move(t.mismatch_cols);
            if (nogap) mark_nogap_triplets(out[i], bad, *code);
            if (nomismatch)
                for (int c : mcols[i])
                    if (c >= 1 && (size_t)c <= drop.size()) drop[(size_t)(c - 1)] = 1;
        } catch (const std::bad_alloc&) {
            failed[i] = 1;
            err[i] = "out of memory";
        } catch (const std::exception& e) {
            failed[i] = 1;
            err[i] = e.what();
        }
        if ((i & 255) == 255) Rcpp::checkUserInterrupt();
    }
    if (nogap || nomismatch) {
        const auto sl = keep_slices(widths, drop, bad);
        for (size_t i = 0; i < n; ++i) {
            if (failed[i]) continue;
            std::string kept;
            for (auto& ab : sl)
                if ((size_t)ab.second <= out[i].size()) kept.append(out[i], (size_t)ab.first, (size_t)(ab.second - ab.first));
            out[i] = std::move(kept);
        }
    }
    Rcpp::CharacterVector o = to_r(out), e = to_r(err), w = to_r(warn);
    Rcpp::List mc(n);
    for (size_t i = 0; i < n; ++i) {
        if (failed[i])
            SET_STRING_ELT(o, (R_xlen_t)i, NA_STRING);
        else
            SET_STRING_ELT(e, (R_xlen_t)i, NA_STRING);
        if (warn[i].empty()) SET_STRING_ELT(w, (R_xlen_t)i, NA_STRING);
        mc[(R_xlen_t)i] = Rcpp::IntegerVector(mcols[i].begin(), mcols[i].end());
    }
    return Rcpp::List::create(Rcpp::_["aligned"] = o, Rcpp::_["error"] = e, Rcpp::_["warning"] = w,
                              Rcpp::_["mismatch_columns"] = mc);
}

// ================================================================================================ translation

// [[Rcpp::export(rng = false)]]
Rcpp::List coalr_translate(SEXP cds, int table) {
    if (TYPEOF(cds) != STRSXP) Rcpp::stop("internal: cds must be a character vector");
    const GeneticCode* code = need_code(table);
    const size_t n = (size_t)XLENGTH(cds);
    Rcpp::CharacterVector protein(n), status(n), codestr(n), note(n);
    Rcpp::IntegerVector nt(n), aa(n), term(n), stops(n), amb(n), trail(n);
    for (size_t k = 0; k < n; ++k) {
        const R_xlen_t i = (R_xlen_t)k;
        if (is_na(cds, i)) {
            protein[i] = NA_STRING;
            status[i] = "failed";
            codestr[i] = "cds_na";
            note[i] = "the CDS is NA";
            nt[i] = aa[i] = term[i] = stops[i] = amb[i] = trail[i] = NA_INTEGER;
            continue;
        }
        const dnds::Translation t = dnds::translate_cds(str_at(cds, i), *code);
        nt[i] = (int)t.nt_length;
        term[i] = t.terminal_stop ? 1 : 0;
        stops[i] = (int)t.internal_stops;
        amb[i] = (int)t.ambiguous_codons;
        trail[i] = t.trailing_nt;
        if (!t.error.empty()) {
            protein[i] = NA_STRING;
            aa[i] = 0;
            status[i] = "failed";
            codestr[i] = t.error_code;
            note[i] = t.error;
            continue;
        }
        SET_STRING_ELT(protein, i, Rf_mkCharLen(t.protein.data(), (int)t.protein.size()));
        aa[i] = (int)t.protein.size();
        const std::string nn = dnds::translation_note(t);
        if (nn.empty()) {
            status[i] = "ok";
            codestr[i] = "";
            note[i] = "";
        } else {
            status[i] = "warning";
            codestr[i] = t.trailing_nt ? (t.internal_stops ? "cds_length;internal_stop" : "cds_length") : "internal_stop";
            note[i] = nn;
        }
        if ((k & 1023) == 1023) Rcpp::checkUserInterrupt();
    }
    return Rcpp::List::create(Rcpp::_["protein"] = protein, Rcpp::_["nt_length"] = nt, Rcpp::_["aa_length"] = aa,
                              Rcpp::_["terminal_stop"] = term, Rcpp::_["internal_stops"] = stops,
                              Rcpp::_["ambiguous_codons"] = amb, Rcpp::_["trailing_nt"] = trail,
                              Rcpp::_["status"] = status, Rcpp::_["code"] = codestr, Rcpp::_["note"] = note);
}

// ================================================================================================ info

// [[Rcpp::export(rng = false)]]
Rcpp::List coalr_info() {
    std::vector<std::string> methods;
    for (kaks::Method m : kaks::output_order()) methods.push_back(kaks::method_name(m));
    return Rcpp::List::create(Rcpp::_["coal_version"] = std::string(COAL_VERSION_STRING),
                              Rcpp::_["methods"] = to_r(methods),
                              Rcpp::_["genetic_codes"] = Rcpp::IntegerVector(supported_tables().begin(),
                                                                             supported_tables().end()),
                              Rcpp::_["matrices"] = to_r(substitution_matrix_names()),
                              Rcpp::_["hardware_threads"] = (int)std::thread::hardware_concurrency());
}
