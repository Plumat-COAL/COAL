// NCBI genetic codes supported by PAL2NAL v14, and IUPAC codon matching.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace coal {

// Codon index in NCBI order: first base T,C,A,G; then second; then third.
inline int base_index(unsigned char c) {
    switch (c) {
        case 'T': return 0;
        case 'C': return 1;
        case 'A': return 2;
        case 'G': return 3;
        default: return -1;
    }
}

struct GeneticCode {
    int number = 0;
    char aa[64] = {};          // amino acid for each codon ('*' = stop)
    uint64_t starts = 0;       // start codons (bit per codon index)
    uint64_t codons[128] = {}; // amino acid letter -> its codons
};

// nullptr if the table number is not supported.
const GeneticCode* get_code(int number);
const std::vector<int>& supported_tables();

// IUPAC alphabet used for DNA after upper-casing and U->T: 15 letters.
// iupac_code(c) is 0..14, or 255 for any other byte.
uint8_t iupac_code(unsigned char c);
// All concrete codons an IUPAC triplet can stand for (bit per codon index).
// idx = a*225 + b*15 + c with a, b, c from iupac_code.
uint64_t iupac_triplet_mask(int idx);

// Whole lookup tables, for hot loops (avoid a call per base).
const uint8_t* iupac_table();          // 256 entries
const uint64_t* iupac_triplet_table(); // 15*15*15 entries

// The first-Met codons PAL2NAL accepts in fallback (anchor) mode.
uint64_t anchor_start_mask();

}  // namespace coal
