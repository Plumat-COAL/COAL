// coal-cpp-dump: test helper for tests/test_cpp_*.py. Exposes the C++ readers and genetic code tables
// as JSON on stdout so the tests can compare them with the reviewed Python engine.
//
//   coal-cpp-dump protein FILE            records of a protein alignment (FASTA or CLUSTAL), raw rows
//   coal-cpp-dump dna FILE [FILE ...]     DNA records: IDs and cleaned sequences
//   coal-cpp-dump tables                  every supported genetic code: amino acids + start flags
//   coal-cpp-dump pool N JOBS BATCH [FAIL_AT]
//                                         the ordered thread pool (pool.hpp) over items 0..N-1 with a
//                                         stand-in per-item function; prints "i value" per item, then a
//                                         JSON line with the pool's stats (or {"error": ...}, exit 1)
// Errors are printed as {"error": "<message>"} (the message the coal-cpp CLI would show).
#include <atomic>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "json.hpp"
#include "pool.hpp"
#include "seqio.hpp"
#include "tables.hpp"

using namespace coal;

static int dump_protein(const std::string& path) {
    ProteinReader pr(path);
    Json recs = Json::array();
    std::string id, row;
    while (pr.next(id, row)) {
        Json r = Json::array();
        r.push(id);
        r.push(row);
        recs.push(r);
    }
    Json out = Json::object();
    out.set("records", recs);
    if (is_clustal(path)) out.set("mask", clustal_mask(path));  // CPP-13
    std::cout << json_dump(out) << '\n';
    return 0;
}

static int dump_dna(const std::vector<std::string>& paths) {
    DnaIndex dna(paths);
    Json recs = Json::array();
    for (size_t i = 0; i < dna.size(); ++i) {
        Json r = Json::array();
        r.push(dna.ids()[i]);
        r.push(dna.get(i));
        recs.push(r);
    }
    Json out = Json::object();
    out.set("records", recs);
    std::cout << json_dump(out) << '\n';
    return 0;
}

static int dump_tables() {
    Json out = Json::object();
    for (int t : supported_tables()) {
        const GeneticCode* g = get_code(t);
        std::string aa(g->aa, 64), starts(64, '-');
        for (int c = 0; c < 64; ++c)
            if (g->starts >> c & 1) starts[(size_t)c] = 'M';
        Json x = Json::object();
        x.set("aa", aa);
        x.set("starts", starts);
        out.set(std::to_string(t), x);
    }
    std::cout << json_dump(out) << '\n';
    return 0;
}


// Stand-in per-sequence function: an LCG run for a variable number of steps, so items finish out of
// order; the first item of every batch is much slower than the rest.
static uint64_t stand_in(uint64_t i, bool slow) {
    uint64_t x = i, steps = (i % 50 + 1) * 200 + (slow ? 4000000 : 0);
    for (uint64_t k = 0; k < steps; ++k) x = x * 6364136223846793005ULL + 1442695040888963407ULL;
    return x;
}

struct PoolItem {
    long long index = 0;
    bool first = false;
    uint64_t value = 0;
    long long finished = -1;  // completion order
};

static int dump_pool(long long n, int jobs, long long batch, long long fail_at) {
    long long next = 0;
    std::atomic<long long> finished{0};
    bool reordered = false;
    long long expect = 0;
    auto read_batch = [&](std::vector<PoolItem>& b) {
        b.clear();
        while (next < n && (long long)b.size() < batch) {
            PoolItem it;
            it.index = next++;
            it.first = b.empty();
            b.push_back(it);
        }
    };
    auto work = [&](PoolItem& it) {
        if (it.index == fail_at) throw std::runtime_error("stand-in failure at item " + std::to_string(it.index));
        it.value = stand_in((uint64_t)it.index, it.first);
        it.finished = finished++;
    };
    long long base = 0;
    auto emit = [&](std::vector<PoolItem>& b) {
        for (size_t k = 0; k < b.size(); ++k) {
            if (b[k].index != expect++) throw std::runtime_error("out of order");
            if (b[k].finished != base + (long long)k) reordered = true;
            std::cout << b[k].index << ' ' << b[k].value << '\n';
        }
        base += (long long)b.size();
    };
    PoolStats st = run_ordered<PoolItem>(read_batch, work, emit, jobs);
    Json o = Json::object();
    o.set("batches", (long long)st.batches);
    o.set("items", (long long)st.items);
    o.set("max_live", (long long)st.max_live);
    o.set("max_threads", (long long)st.max_threads);
    o.set("threads_default", (long long)resolve_jobs(0));
    o.set("reordered", reordered);
    std::cout << json_dump(o) << '\n';
    return 0;
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    if (argc < 2) {
        std::cerr << "usage: coal-cpp-dump protein FILE | dna FILE... | tables | pool N JOBS BATCH [FAIL_AT]\n";
        return 2;
    }
    std::string cmd = argv[1];
    try {
        if (cmd == "protein" && argc == 3) return dump_protein(argv[2]);
        if (cmd == "dna" && argc >= 3) return dump_dna(std::vector<std::string>(argv + 2, argv + argc));
        if (cmd == "tables") return dump_tables();
        if (cmd == "pool" && (argc == 5 || argc == 6))
            return dump_pool(std::stoll(argv[2]), std::stoi(argv[3]), std::stoll(argv[4]), argc == 6 ? std::stoll(argv[5]) : -1);
    } catch (const std::exception& e) {
        Json o = Json::object();
        o.set("error", std::string(e.what()));
        std::cout << json_dump(o) << '\n';
        return 1;
    }
    std::cerr << "coal-cpp-dump: bad arguments\n";
    return 2;
}
