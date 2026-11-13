# dnds: two species, 12 ortholog pairs, for `coal pipeline`, `coal dnds` and `coal kaks`

Coding sequences of two "species" with known orthologs (**synthetic, not real genes**; `examples/make_dnds_example.py`
makes them from a fixed seed with a codon-substitution simulator, described at the top of that script).

| File | Content |
|---|---|
| `species_a_cds.fasta`, `species_b_cds.fasta` | 13 CDS each: the 12 orthologous genes (`spA_g01` ... / `spB_g01` ..., in another order in the second file) and one gene without an ortholog |
| `orthologs.tsv` | the true pairs, and what each was simulated with: length in codons, omega (dN/dS), dS, and the size of an in-frame insertion or deletion in species B |
| `codon.fasta` | the 12 pairwise codon alignments as COAL writes them (`--keep-alignments`), records 1+2, 3+4, ... |

The genes differ in length (45 to 480 codons), in the simulated dN/dS (0.04 to 1.2) and in dS (0.08 to 0.9), so one
run shows a range of values. Source and licence: generated for COAL, part of COAL (GPL-3.0-only); no sequence is
taken from a database.

## From coding sequences to dN/dS in one step

Needs the platform wheel (it has `coal-cpp`) and DIAMOND (`coal pipeline install diamond`):

```
coal pipeline pairs --query examples/dnds/species_a_cds.fasta --subject examples/dnds/species_b_cds.fasta -o out
```

`out/dnds.tsv` then has the 12 reciprocal best hits (YN, the default method; columns shortened here):

| seq1 | seq2 | dN | dS | dNdS | codons_compared | prot_identity | status |
|---|---|---|---|---|---|---|---|
| spA_g01 | spB_g01 | 0.0102 | 0.399 | 0.0255 | 320 | 97.81 | ok |
| spA_g02 | spB_g02 | 0.0478 | 0.723 | 0.0662 | 207 | 88.57 | ok |
| spA_g03 | spB_g03 | 0.0289 | 0.230 | 0.126 | 480 | 93.75 | ok |
| spA_g04 | spB_g04 | 0.0828 | 0.523 | 0.158 | 150 | 82.89 | ok |
| spA_g05 | spB_g05 | 0.204 | 1.003 | 0.204 | 258 | 66.41 | ok |
| spA_g06 | spB_g06 | 0.118 | 0.277 | 0.427 | 384 | 77.18 | ok |
| spA_g07 | spB_g07 | 0.104 | 0.272 | 0.380 | 180 | 81.11 | ok |
| spA_g08 | spB_g08 | 0.305 | 0.457 | 0.667 | 300 | 54.28 | ok |
| spA_g09 | spB_g09 | 0.152 | 0.274 | 0.553 | 120 | 70.83 | ok |
| spA_g10 | spB_g10 | 0.364 | 0.427 | 0.853 | 240 | 46.67 | ok |
| spA_g11 | spB_g11 | 0.102 | 0.0676 | 1.50 | 200 | 81.00 | ok |
| spA_g12 | spB_g12 | 0.133 | 0.160 | 0.833 | 45 | 71.11 | ok |

and `out/unpaired.tsv` lists `spA_orphan` and `spB_orphan` (`no_hit`). The estimates scatter around the simulated
values, as estimates from one pair of sequences do, and most for the short genes: `g12` (45 codons, simulated
omega 0.3) comes out at 0.83, and `g11` (simulated 1.2, but only dS 0.08, so few synonymous changes) at 1.5.
`--min-codons 50` flags `g12` (status `warning`, note `fewer than 50 codons compared`); see "Filtering before
interpretation" in `docs/kaks.md`.

Without DIAMOND, with the pairs given:

```
coal pipeline pairs --query examples/dnds/species_a_cds.fasta --subject examples/dnds/species_b_cds.fasta \
     --pairs examples/dnds/orthologs.tsv -o out --keep-alignments
```

(`out/codon_alignments.fasta` is then `examples/dnds/codon.fasta`.)

## From a codon alignment

```
coal dnds examples/dnds/codon.fasta --pairs consecutive                       # the tidy table, on standard output
coal kaks -i examples/dnds/codon.fasta -o out.kaks --pairs consecutive -m YN  # KaKs_Calculator's format
```

Both give the numbers of the table above (`coal kaks` with 6 significant digits).
