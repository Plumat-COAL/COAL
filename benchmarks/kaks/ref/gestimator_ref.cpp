// SPDX-License-Identifier: GPL-3.0-only
// Reference values for Method::Comeron: orthologr 0.4.3's own copies of the libsequence headers (Comeron95.h and
// what it includes), compiled UNMODIFIED from an orthologr source tree and called exactly as orthologr's
// gestimator() calls them (maxHits = 3, one RedundancyCom95 object), but printed with 17 significant digits
// (gestimator() prints 4). Nothing of orthologr or libsequence is copied into COAL's repository for this: the
// headers are found with -I.
//
//   g++ -O2 -I benchmarks/kaks/ref/stub -I /path/to/orthologr/src benchmarks/kaks/ref/gestimator_ref.cpp \
//       -o gestimator_ref
//   gestimator_ref pairs.axt > ref.tsv
//
// Input: AXT (name line, two aligned rows, blank line). Output, one line per pair:
//   name dN dS ratio P0 P2S P2V P4 Q0 Q2S Q2V Q4 L0 L2S L2V L4     (999 = gestimator's "cannot be calculated")
// A pair that gestimator 0.4.3 cannot take (rows of unequal length, shorter than one codon, or not a multiple
// of 3: it reads beyond the rows or aborts) is written as "name ERR".
#include <Rcpp.h>

#include <cstdio>
#include <fstream>
#include <utility>

#include "Comeron95.h"

int main(int argc, char** argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: gestimator_ref pairs.axt\n");
        return 2;
    }
    std::ifstream in(argv[1]);
    if (!in) {
        std::fprintf(stderr, "cannot read %s\n", argv[1]);
        return 1;
    }
    RedundancyCom95 code_degeneracy;
    std::string line;
    std::vector<std::string> rec;
    auto flush = [&]() {
        if (rec.empty()) return;
        if (rec.size() != 3) {
            std::printf("%s\tERR\n", rec[0].c_str());
            rec.clear();
            return;
        }
        std::pair<std::string, std::string> a("1", rec[1]), b("2", rec[2]);
        const size_t n = rec[1].size();
        if (n != rec[2].size() || n < 3 || n % 3 != 0) {
            std::printf("%s\tERR\n", rec[0].c_str());
        } else {
            Comeron95 C(&a, &b, 3, &code_degeneracy);
            std::printf("%s\t%.17g\t%.17g\t%.17g", rec[0].c_str(), C.dn(), C.ds(), C.ratio());
            const double v[12] = {C.P0(), C.P2S(), C.P2V(), C.P4(), C.Q0(),  C.Q2S(),
                                  C.Q2V(), C.Q4(), C.L0(),  C.L2S(), C.L2V(), C.L4()};
            for (double x : v) std::printf("\t%.17g", x);
            std::printf("\n");
        }
        rec.clear();
    };
    while (std::getline(in, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
        if (line.empty()) {
            flush();
            continue;
        }
        rec.push_back(line);
    }
    flush();
    return 0;
}
