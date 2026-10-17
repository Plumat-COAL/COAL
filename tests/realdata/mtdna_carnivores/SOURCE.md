# mtdna_carnivores

- **Source:** HyPhy (github.com/veg/hyphy), commit `1fa8fd0469328b8d6df2c7273077170765beee5f`
  - `tests/hbltests/data/mtDNA.fas` -> `mtDNA.fas.gz` (gzip of the unchanged file; 677,593 bytes, sha256 `3b0d6a2cf3292d82fc1989fb5dba4b161170baa36776aaa9a1b014e59442f44f`)
    <https://github.com/veg/hyphy/blob/1fa8fd0469328b8d6df2c7273077170765beee5f/tests/hbltests/data/mtDNA.fas>
- **Licence:** MIT (copyright the HyPhy core developers; full text: ../LICENSE-hyphy.txt). The sequences themselves are public GenBank/INSDC records.
- **Organism / gene:** Concatenated mitochondrial protein-coding genes of 62 carnivores
- **Size:** 62 × 3,623 codons
- **Genetic code (-codontable):** 2 (vertebrate mitochondrial)
- **Kind in datasets.json:** roundtrip
- **Notes:** Long genes; an internal stop column (TAA, codon 227: a gene boundary) in every row plus a few more; 154 ambiguous codons (N, IUPAC). Newick tree appended (ignored).
