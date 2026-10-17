// coal-cpp-dump: test helper for tests/test_cpp_*.py. Exposes the C++ readers and genetic code tables
// as JSON on stdout so the tests can compare them with the reviewed Python engine.
//
//   coal-cpp-dump protein FILE            records of a protein alignment (FASTA or CLUSTAL), raw rows
//   coal-cpp-dump dna FILE [FILE ...]     DNA records: IDs and cleaned sequences
//   coal-cpp-dump tables                  every supported genetic code: amino acids + start flags
// Errors are printed as {"error": "<message>"} (the message the coal-cpp CLI would show).
#include <iostream>
#include <string>
#include <vector>

#include "json.hpp"
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

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    if (argc < 2) {
        std::cerr << "usage: coal-cpp-dump protein FILE | dna FILE... | tables\n";
        return 2;
    }
    std::string cmd = argv[1];
    try {
        if (cmd == "protein" && argc == 3) return dump_protein(argv[2]);
        if (cmd == "dna" && argc >= 3) return dump_dna(std::vector<std::string>(argv + 2, argv + argc));
        if (cmd == "tables") return dump_tables();
    } catch (const std::exception& e) {
        Json o = Json::object();
        o.set("error", std::string(e.what()));
        std::cout << json_dump(o) << '\n';
        return 1;
    }
    std::cerr << "coal-cpp-dump: bad arguments\n";
    return 2;
}
