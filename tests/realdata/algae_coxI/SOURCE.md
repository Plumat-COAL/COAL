# algae_coxI

- **Source:** HyPhy (github.com/veg/hyphy), commit `1fa8fd0469328b8d6df2c7273077170765beee5f`
  - `tests/hbltests/libv3/data/algae_coxI_aligned.fas` -> `algae_coxI_aligned.fas.gz` (gzip of the unchanged file; 13,314 bytes, sha256 `a7156861385f9500437d22ed65c95bd43ef888a5d74d92c10f42f7a2af78e75f`)
    <https://github.com/veg/hyphy/blob/1fa8fd0469328b8d6df2c7273077170765beee5f/tests/hbltests/libv3/data/algae_coxI_aligned.fas>
- **Licence:** MIT (copyright the HyPhy core developers; full text: ../LICENSE-hyphy.txt). The sequences themselves are public GenBank/INSDC records.
- **Organism / gene:** Mitochondrial COX1 of 19 algae (GenBank KT7799xx)
- **Size:** 19 × 214 codons
- **Genetic code (-codontable):** 4 (mold/protozoan mitochondrial; HyPhy's algae-mtDNA.wbf uses "Mold-Protozoan-mtDNA")
- **Kind in datasets.json:** roundtrip
- **Notes:** 64 TGA (= Trp in code 4) codons, so code 1 fails. Newick tree appended (ignored by the loader).
