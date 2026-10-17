#include "tables.hpp"

#include <array>

namespace coal {
namespace {

struct Raw {
    int number;
    const char* aas;
    const char* starts;
};

// Source: NCBI "The Genetic Codes" (same data as src/coal/tables.py).
const Raw RAW[] = {
    {1, "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "---M------**--*----M---------------M----------------------------"},
    {2, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        "----------**--------------------MMMM----------**---M------------"},
    // ENG-g: starts ATA ATG GTG (current NCBI)
    {3, "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "----------**----------------------MM---------------M------------"},
    {4, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "--MM------**-------M------------MMMM---------------M------------"},
    {5, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        "---M------**--------------------MMMM---------------M------------"},
    {6, "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "--------------*--------------------M----------------------------"},
    {9, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "----------**-----------------------M---------------M------------"},
    {10, "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "----------**-----------------------M----------------------------"},
    {11, "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "---M------**--*----M------------MMMM---------------M------------"},
    {12, "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "----------**--*----M---------------M----------------------------"},
    {13, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
         "---M------**----------------------MM---------------M------------"},
    {14, "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
         "-----------*-----------------------M----------------------------"},
    {15, "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "----------*---*--------------------M----------------------------"},
    {16, "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "----------*---*--------------------M----------------------------"},
    {21, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
         "----------**-----------------------M---------------M------------"},
    {22, "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "------*---*---*--------------------M----------------------------"},
    {23, "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "--*-------**--*-----------------M--M---------------M------------"},
};
constexpr int NRAW = sizeof(RAW) / sizeof(RAW[0]);

const char* IUPAC_LETTERS = "ACGTRYSWKMBDHVN";
const char* IUPAC_BASES[15] = {"A", "C", "G", "T", "AG", "CT", "CG", "AT", "GT", "AC",
                               "CGT", "AGT", "ACT", "ACG", "ACGT"};

struct Tables {
    std::array<GeneticCode, NRAW> codes;
    std::vector<int> numbers;
    std::array<uint8_t, 256> iupac;
    std::array<uint64_t, 15 * 15 * 15> triplet;
    uint64_t anchor_start = 0;

    Tables() {
        for (int i = 0; i < NRAW; ++i) {
            GeneticCode& g = codes[i];
            g.number = RAW[i].number;
            numbers.push_back(g.number);
            for (int c = 0; c < 64; ++c) {
                g.aa[c] = RAW[i].aas[c];
                g.codons[(unsigned char)g.aa[c]] |= uint64_t(1) << c;
                if (RAW[i].starts[c] == 'M') g.starts |= uint64_t(1) << c;
            }
        }
        iupac.fill(255);
        for (int i = 0; i < 15; ++i) iupac[(unsigned char)IUPAC_LETTERS[i]] = (uint8_t)i;
        for (int a = 0; a < 15; ++a)
            for (int b = 0; b < 15; ++b)
                for (int c = 0; c < 15; ++c) {
                    uint64_t m = 0;
                    for (const char* x = IUPAC_BASES[a]; *x; ++x)
                        for (const char* y = IUPAC_BASES[b]; *y; ++y)
                            for (const char* z = IUPAC_BASES[c]; *z; ++z)
                                m |= uint64_t(1) << (base_index(*x) * 16 + base_index(*y) * 4 + base_index(*z));
                    triplet[a * 225 + b * 15 + c] = m;
                }
        for (const char* s : {"ATG", "CTG", "GTG"})
            anchor_start |= uint64_t(1) << (base_index(s[0]) * 16 + base_index(s[1]) * 4 + base_index(s[2]));
    }
};

const Tables& tables() {
    static const Tables t;
    return t;
}

}  // namespace

const GeneticCode* get_code(int number) {
    for (const GeneticCode& g : tables().codes)
        if (g.number == number) return &g;
    return nullptr;
}

const std::vector<int>& supported_tables() { return tables().numbers; }

uint8_t iupac_code(unsigned char c) { return tables().iupac[c]; }

uint64_t iupac_triplet_mask(int idx) { return tables().triplet[idx]; }

uint64_t anchor_start_mask() { return tables().anchor_start; }

const uint8_t* iupac_table() { return tables().iupac.data(); }

const uint64_t* iupac_triplet_table() { return tables().triplet.data(); }

}  // namespace coal
