# hiv_rt

- **Source:** HyPhy (github.com/veg/hyphy), commit `1fa8fd0469328b8d6df2c7273077170765beee5f`
  - `tests/data/HIV_RT.nex` -> `HIV_RT.nex.gz` (gzip of the unchanged file; 497,214 bytes, sha256 `7364b7170cbfb44737635d36f2e130c7c7e3a739abbda3a81cc4495de8624dea`)
    <https://github.com/veg/hyphy/blob/1fa8fd0469328b8d6df2c7273077170765beee5f/tests/data/HIV_RT.nex>
- **Licence:** MIT (copyright the HyPhy core developers; full text: ../LICENSE-hyphy.txt). The sequences themselves are public GenBank/INSDC records.
- **Organism / gene:** HIV-1 reverse transcriptase, 476 sequences
- **Size:** 476 × 335 codons
- **Genetic code (-codontable):** 1 (standard)
- **Kind in datasets.json:** roundtrip
- **Notes:** 3,110 IUPAC-ambiguous codons (population sequencing). PAL2NAL v14 prints 113 spurious 'does not correspond' warnings (YTA/YTG = L, ATM/ATW = I, MGA/MGG = R), COAL none: the documented IUPAC difference; output identical.
