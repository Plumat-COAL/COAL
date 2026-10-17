# ciliates

- **Source:** HyPhy (github.com/veg/hyphy), commit `1fa8fd0469328b8d6df2c7273077170765beee5f`
  - `tests/hbltests/libv3/data/ciliates.fas` -> `ciliates.fas.gz` (gzip of the unchanged file; 9,552 bytes, sha256 `f868d491bf805b8335c18019d7f6ad8f3c99c139e618629b097a14512d79f90e`)
    <https://github.com/veg/hyphy/blob/1fa8fd0469328b8d6df2c7273077170765beee5f/tests/hbltests/libv3/data/ciliates.fas>
- **Licence:** MIT (copyright the HyPhy core developers; full text: ../LICENSE-hyphy.txt). The sequences themselves are public GenBank/INSDC records.
- **Organism / gene:** One ortholog group (OG5_149820) from 16 ciliate transcriptome contigs
- **Size:** 16 × 153 codons
- **Genetic code (-codontable):** 6 (ciliate nuclear; HyPhy's tests/hbltests/libv3/ciliate-code.wbf uses "Ciliate-Nuclear")
- **Kind in datasets.json:** roundtrip
- **Notes:** FASTA with a Newick tree appended after the last record (HyPhy convention; the loader ignores lines starting with '('). The alignment contains no TAA/TAG codon, so code 6 and code 1 give the same protein here: this set does not distinguish them.
