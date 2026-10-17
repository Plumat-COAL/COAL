# egfr

- **Source:** Biopython (github.com/biopython/biopython), commit `4eda92397bdedff82f04ac218a862ae25ebdc4a8`
  - `Tests/codonalign/egfr_pro.aln` -> `egfr_pro.aln.gz` (gzip of the unchanged file; 16,569 bytes, sha256 `1b9225c60ac50f00f9b1dc1ce23cdfc05fed71f6a3147120191dfbdac9e25d8b`)
    <https://github.com/biopython/biopython/blob/4eda92397bdedff82f04ac218a862ae25ebdc4a8/Tests/codonalign/egfr_pro.aln>
  - `Tests/codonalign/egfr_nucl.fa` -> `egfr_nucl.fa.gz` (gzip of the unchanged file; 29,173 bytes, sha256 `40de3eefa413a9ac90e84516fb7aeb1023007f8f399f9e92e25e45c14c836ea0`)
    <https://github.com/biopython/biopython/blob/4eda92397bdedff82f04ac218a862ae25ebdc4a8/Tests/codonalign/egfr_nucl.fa>
  - `Tests/codonalign/egfr_id` -> `egfr_id.gz` (gzip of the unchanged file; 345 bytes, sha256 `aa9fb01b7f6536b6c044f3d0a0080c806a37fb099c8968c82e17211f03cb0352`)
    <https://github.com/biopython/biopython/blob/4eda92397bdedff82f04ac218a862ae25ebdc4a8/Tests/codonalign/egfr_id>
- **Licence:** Biopython License Agreement / BSD 3-Clause (dual; full text: ../LICENSE-biopython.rst). The sequences themselves are public GenBank/RefSeq records.
- **Organism / gene:** EGFR: 6 proteins (CLUSTAL, NP_/gb IDs) and their mRNAs with UTRs (NM_/gb IDs, full NCBI descriptions; another order)
- **Size:** 6 × 1,482 alignment columns
- **Genetic code (-codontable):** 1 (standard)
- **Kind in datasets.json:** pair
- **Notes:** A real protein-alignment + mRNA pair, not a codon alignment. egfr_id maps protein ID -> mRNA ID; the test renames the mRNAs with it. One codon in one gene does not match its residue (1 warning, the same in PAL2NAL). Without the map, pairing is by file order and the run fails.
