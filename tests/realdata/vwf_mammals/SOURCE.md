# vwf_mammals

- **Source:** HyPhy (github.com/veg/hyphy), commit `1fa8fd0469328b8d6df2c7273077170765beee5f`
  - `tests/data/vwf.nex` -> `vwf.nex.gz` (gzip of the unchanged file; 77,459 bytes, sha256 `fcb454a8a08315bc925d4050fed07c89600c9e025ffb819621072660ec324cd3`)
    <https://github.com/veg/hyphy/blob/1fa8fd0469328b8d6df2c7273077170765beee5f/tests/data/vwf.nex>
- **Licence:** MIT (copyright the HyPhy core developers; full text: ../LICENSE-hyphy.txt). The sequences themselves are public GenBank/INSDC records.
- **Organism / gene:** von Willebrand factor, partial coding sequences of 62 mammals
- **Size:** 62 × 392 codons
- **Genetic code (-codontable):** 1 (standard)
- **Kind in datasets.json:** partial
- **Notes:** Not strictly codon-aligned at fragment ends: 15 partly gapped codons (e.g. '--G'). kind=partial: those become a protein gap and their bases stay in the DNA (a fragment translated in its own frame); COAL must skip them like a UTR.
