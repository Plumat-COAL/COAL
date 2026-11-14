# Tutorial 2: two genomes in, a dN/dS table out (`coal pipeline pairs`)

**Goal:** go from the coding sequences of two bacterial genomes to a table of dN, dS and dN/dS for every pair of
orthologs, and learn to read it: the statuses, saturation, what to filter, and how to reproduce the run.

| | |
|---|---|
| Data | *Escherichia coli* K-12 MG1655 (4,318 CDS) and *Salmonella enterica* serovar Typhimurium LT2 (4,554 CDS), NCBI RefSeq |
| Needs | the platform wheel (`coal-cpp`), DIAMOND 2.1.0 or later (`coal pipeline install diamond`), Python 3 for the summary script |
| Runtime | about 16 s with one thread (`-j 1`): DIAMOND 6 s, alignment + Ka/Ks 10 s; more threads are faster |
| Tested with | COAL 0.6.0.dev0 (commit 009bdd5), DIAMOND 2.2.8, Linux x86_64 |

*E. coli* and *Salmonella* split more than 100 million years ago, and this pair has been a textbook case for
synonymous and nonsynonymous divergence since the 1980s: their genes are close enough to align over their whole
length, and their synonymous sites are far enough apart to show what saturation looks like.

## 0. Set up

```
export T=/path/to/COAL/docs/tutorials        # see the README of this folder
cd coal-tutorials                             # the scratch folder of tutorial 1, or a new one
bash $T/data/get_enterobacteria.sh
```

```
ok                 coal-tutorial-data/ecoli_K12_MG1655.cds.fna.gz
ok                 coal-tutorial-data/styphimurium_LT2.cds.fna.gz
```

(GCF_000005845.2 and GCF_000006945.2, NCBI's `cds_from_genomic.fna.gz`; checksums in the script.) COAL reads
gzip-compressed FASTA as it is.

## 1. Run it

```
coal pipeline pairs \
    --query   coal-tutorial-data/ecoli_K12_MG1655.cds.fna.gz \
    --subject coal-tutorial-data/styphimurium_LT2.cds.fna.gz \
    -o ecoli_vs_salmonella -j 1
```

```
coal pipeline: 2 gzip-compressed input file(s) decompressed into .../ecoli_vs_salmonella/work/inputs
coal pipeline: translated 4,318 query and 4,554 subject records (0 not usable)
coal pipeline: warning: the input looks like several transcripts per gene (query: 4,318 records of 4,308 genes, 9 genes with several records (e.g. b0149: lcl|NC_000913.3_cds_NP_414691.1_149, lcl|NC_000913.3_cds_YP_010051172.1_150), read from the header field [locus_tag=). The pipeline treats each record as a gene: such a gene is compared once per transcript, and reciprocal best hits are split between its transcripts. Use --longest-isoform to keep the longest CDS of each gene (with --gene-map TSV if the genes are not in the headers), or give one sequence per gene
coal pipeline: 3,142 pairs (rbh)
coal pipeline: summary: 8,872 records (0 not usable); 3,142 pairs: 2,973 ok, 3 warning, 11 not_converged, 155 saturated, 0 no_data, 0 error
coal pipeline: summary: decompress_inputs 0.0 s, translate_query 0.0 s, translate_subject 0.1 s, diamond_search 6.0 s, orthologs 0.0 s, align_codon_kaks 9.7 s
coal pipeline: summary: results in .../ecoli_vs_salmonella/dnds.tsv
coal pipeline: summary: what to cite for this run: coal cite .../ecoli_vs_salmonella/report.json
```

Exit status 0. What happened, step by step: every CDS was translated (standard code, `-c 1`; bacterial code 11
has the same codons), DIAMOND searched the proteins both ways, the 3,142 **reciprocal best hits** became the
pairs, each pair's proteins were aligned (global Needleman-Wunsch, BLOSUM62), the CDS were threaded onto that
alignment, and Ka/Ks was estimated with Yang & Nielsen's (2000) method, YN00.

**The warning about transcripts** is worth reading rather than skipping. *E. coli* has no splicing, so these are
not transcripts: they are 9 genes annotated with two or three alternative start sites (*mrcB*'s PBP1B alpha and
gamma, *infB*'s three IF2 forms, ...). Here they cost nothing (the pairs are the same with and without the fix),
but in a eukaryote with five isoforms per gene they split the reciprocal best hits; `--longest-isoform` keeps
one record per gene (step 6).

## 2. What is in the folder

```
ls ecoli_vs_salmonella
```

```
dnds.tsv  failed.tsv  genes.tsv  orthologs.tsv  report.json  unpaired.tsv  work
```

| File | Rows here | What it is |
|---|---|---|
| `dnds.tsv` | 3,142 | the result: one row per pair (and per method with several `-m`) |
| `orthologs.tsv` | 3,142 | the pairs, with DIAMOND's `pident length evalue bitscore` |
| `unpaired.tsv` | 2,588 | genes without a partner: `no_hit` (884 *E. coli*, 1,081 *Salmonella*) or `not_reciprocal` (292, 331) |
| `genes.tsv` | 8,872 | one row per CDS record: length, stops, status |
| `failed.tsv` | 0 | anything that could not be used (only the header here) |
| `report.json` | | everything needed to describe and repeat the run |
| `work/` | | proteins, DIAMOND tables, logs of every external command, the `--resume` state |

About 70% of the genes of each genome have a reciprocal best hit. The `no_hit` genes include what makes the
two species different (prophages, genomic islands, horizontally acquired genes); `not_reciprocal` genes are usually members of
paralogous families, where the best hit of A is B but B's best hit is A's paralog.

## 3. The table

```
head -3 ecoli_vs_salmonella/dnds.tsv | cut -f1-6,13-16
```

```
seq1	seq2	method	dN	dS	dNdS	codons_compared	prot_aln_length	prot_identity	status
lcl|NC_000913.3_cds_NP_414543.1_2	lcl|NC_003197.2_cds_NP_459007.1_2	YN	0.031469968	1.9097809	0.016478313	820	820	94.51	ok
lcl|NC_000913.3_cds_NP_414544.1_3	lcl|NC_003197.2_cds_NP_459008.1_3	YN	0.035680572	2.0166797	0.017692731	309	310	92.90	ok
```

The first row is *thrA* (aspartokinase / homoserine dehydrogenase): the proteins are 94.5% identical over 820
columns, 820 codons were compared, dN = 0.031 nonsynonymous substitutions per nonsynonymous site, dS = 1.91
synonymous substitutions per synonymous site, so dN/dS = 0.016: strong purifying selection, as for most genes.
The other columns (`SE_dN`, `SE_dS`, `N_sites`, `S_sites`, `Nd`, `Sd`, `note`) are described in
[pipeline.md](../pipeline.md#dndstsv).

## 4. The statuses

Count them:

```
cut -f16 ecoli_vs_salmonella/dnds.tsv | sort | uniq -c
```

```
     11 not_converged
   2973 ok
    155 saturated
      1 status
      3 warning
```

| Status | Here | What it means | What to do |
|---|---|---|---|
| `ok` | 2,973 | the method gave dN and dS | still filter (step 5) |
| `warning` | 3 | estimated, but something in the note needs a look | read the note |
| `saturated` | 155 | dS (here always dS) could not be corrected for multiple hits: `NA` | no estimate; the pair is too divergent for this method |
| `not_converged` | 11 | YN00's iteration did not settle: no estimate | no estimate; the note's "last round" values are *not* estimates |
| `no_data`, `error` | 0 | nothing to compare / the pair could not be computed | see the note and `failed.tsv` |

**The three warnings are selenoproteins**, a nice example of why a warning should be read:

```
awk -F'\t' '$16=="warning"' ecoli_vs_salmonella/dnds.tsv | cut -f1,17
```

```
lcl|NC_000913.3_cds_NP_415991.1_1460	1 codon column(s) with an internal stop codon skipped; seq1: 1 internal stop codon in the CDS; seq2: 1 internal stop codon in the CDS
lcl|NC_000913.3_cds_NP_418330.1_3826	1 codon column(s) with an internal stop codon skipped; seq1: 1 internal stop codon in the CDS; seq2: 1 internal stop codon in the CDS
lcl|NC_000913.3_cds_NP_418503.1_3999	1 codon column(s) with an internal stop codon skipped; seq1: 1 internal stop codon in the CDS; seq2: 1 internal stop codon in the CDS
```

These are *fdnG*, *fdoG* and *fdhF*, the three formate dehydrogenases, whose headers say
`[transl_except=(pos:586..588,aa:Sec)]`: the internal TGA is a selenocysteine. One column was skipped; the
estimate is fine. In `genes.tsv` two other *E. coli* records have internal stops (*insN*, *phnE1*); their
headers say `[pseudo=true]`, and neither found a partner (`no_hit`).

**Saturated and not-converged pairs are the divergent ones.** The median protein identity of the saturated pairs
is 66%, against 90% for the `ok` pairs. When dS is near 2, a small extra divergence tips the correction over the edge; step 5 shows
that many `ok` pairs are close to it too.

## 5. Filter, then interpret

A status of `ok` only says that the method produced numbers. [kaks.md, "Filtering before
interpretation"](../kaks.md#filtering-before-interpretation) gives four checks and starting thresholds. The
summary script of this folder applies them in order and says what each removes:

```
python3 $T/scripts/dnds_summary.py ecoli_vs_salmonella/dnds.tsv \
    --cds coal-tutorial-data/ecoli_K12_MG1655.cds.fna.gz coal-tutorial-data/styphimurium_LT2.cds.fna.gz \
    --write-kept kept.tsv
```

```
ecoli_vs_salmonella/dnds.tsv: 3142 rows of method YN
status: ok 2973, saturated 155, not_converged 11, warning 3
filters (each applied to what the previous ones kept):
  has an estimate (status ok or warning)           2976  (-166)
  codons_compared >= 50                            2953  (-23)
  codons_compared / prot_aln_length >= 0.5         2947  (-6)
  prot_identity >= 40                              2916  (-31)
  dS <= 2                                          2468  (-448)
  dS >= 0.01                                       2468  (-0)
quartiles over the 2468 rows kept (25% / median / 75%):
  dN    0.0254 / 0.0460 / 0.0797
  dS    0.8539 / 1.2228 / 1.4946
  dNdS  0.0258 / 0.0408 / 0.0668
  dN/dS > 1: 3 of 2468
highest dN/dS among the rows kept:
  yqgC       dN 0.2079  dS 0.0524  dN/dS 3.9651  codons   54  protein YqgC
  yqgD       dN 0.3557  dS 0.2598  dN/dS 1.3690  codons   81  DUF2684 domain-containing protein YqgD
  yaaY       dN 0.2527  dS 0.1933  dN/dS 1.3078  codons   71  DUF2575 domain-containing protein YaaY
  rpsT       dN 0.0104  dS 0.0149  dN/dS 0.7009  codons   87  30S ribosomal subunit protein S20
  rpmD       dN 0.0152  dS 0.0229  dN/dS 0.6654  codons   59  50S ribosomal subunit protein L30
  ...
lowest dN/dS among the rows kept:
  ybaB       dN 0.0000  dS 0.6997  dN/dS 0.0000  codons  109  putative nucleoid-associated protein YbaB
  infA       dN 0.0000  dS 0.1284  dN/dS 0.0000  codons   72  translation initiation factor IF-1
  rpmF       dN 0.0000  dS 0.0585  dN/dS 0.0000  codons   57  50S ribosomal subunit protein L32
  ...
wrote 2468 rows to kept.tsv
```

How to read it:

- **The biggest filter is dS ≤ 2** (448 pairs). The median dS of what is left is 1.22: on average every
  synonymous site of these genes has changed about once since the split. That is the upper end of what a
  pairwise method can measure, so treat dS (and the ratio) of individual genes with care; the *distribution*
  is robust, single values less so. A lower cut (`--max-ds 1`) is defensible too; state which you used.
- **dN/dS is low almost everywhere**: median 0.04, 75% of genes below 0.07. Most proteins are under strong
  purifying selection, which is the expected picture for bacterial core genes.
- **The "top hits" are a lesson, not a finding.** The three pairs above 1 are short (54 to 81 codons), poorly
  characterised genes. *yqgC* has dS 0.05: about two synonymous differences, against 10 nonsynonymous ones. A
  pairwise whole-gene ratio from 54 codons is noise; testing for positive selection needs site or branch models
  (codeml, HyPhy) on more than two sequences. Raise `--min-codons` to 100 and they drop out.
- **The ribosomal proteins have a tenth of the usual dS**, a classic observation in exactly these two species
  (Sharp & Li 1987): highly expressed genes use a restricted set of codons, which constrains synonymous change.
  Check it with the script's `--only-ids`:

  ```
  zcat coal-tutorial-data/ecoli_K12_MG1655.cds.fna.gz | grep -E '\[protein=(30S|50S) ribosomal' \
      | sed 's/^>//; s/ .*//' > ribo.ids
  python3 $T/scripts/dnds_summary.py ecoli_vs_salmonella/dnds.tsv --only-ids ribo.ids --top 0
  ```

  ```
  ecoli_vs_salmonella/dnds.tsv: 61 rows of method YN (only the IDs in ribo.ids)
  ...
  quartiles over the 57 rows kept (25% / median / 75%):
    dN    0.0000 / 0.0051 / 0.0121
    dS    0.0581 / 0.0833 / 0.1322
    dNdS  0.0000 / 0.0448 / 0.1192
  ```

  Median dS 0.08 against 1.22 for all genes; dN/dS about the same. This is why dS is a poor molecular clock
  across genes of different expression, and why comparing dN/dS between gene classes needs care.

The same filters are one line in R (from kaks.md) or awk:

```
awk -F'\t' 'NR==1 || ($16=="ok" || $16=="warning") && $13>=50 && $13>=0.5*$14 && $15>=40 && $5!="NA" && $5<=2 && $5>=0.01' \
    ecoli_vs_salmonella/dnds.tsv > kept_awk.tsv
```

### Or let COAL flag the rows

`--min-codons`, `--max-ds` and `--min-identity` mark rows instead of removing them: status `warning`, a note,
numbers unchanged.

```
coal pipeline pairs --query coal-tutorial-data/ecoli_K12_MG1655.cds.fna.gz \
    --subject coal-tutorial-data/styphimurium_LT2.cds.fna.gz -o ecsa_filtered -j 1 \
    --longest-isoform --min-codons 50 --max-ds 2 --min-identity 40 --keep-alignments
```

```
coal pipeline: --longest-isoform: query: 10 of 4,318 records dropped, 4,308 genes (from header field [locus_tag=); subject: 0 of 4,554 records dropped, 4,554 genes (from header field [locus_tag=); 10 in all, listed in genes.tsv with status dropped
coal pipeline: 3,142 pairs (rbh)
coal pipeline: summary: 8,872 records (0 not usable, 10 dropped by --longest-isoform: 8,862 genes); 3,142 pairs: 2,465 ok, 511 warning, 11 not_converged, 155 saturated, 0 no_data, 0 error
...
coal pipeline: summary: flagged (status warning, numbers unchanged): 23 pair(s) with fewer than 50 codons compared, 472 with dS above 2, 38 with protein identity below 40%
```

`report.json` has the same counts (`counts.flagged`), `genes.tsv` lists the 10 dropped records with the record
kept for each gene, and `codon_alignments.fasta` has the 3,142 codon alignments (records 1+2, 3+4, ...), ready for
`coal dnds --pairs consecutive` or another program.

## 6. Reproduce a run

`report.json` records what is needed to describe and repeat the run:

```
python3 - <<'EOF'
import json
r = json.load(open("ecoli_vs_salmonella/report.json"))
print(r["status"], r["exit_code"], r["coal_version"], r["started"])
print("coal " + " ".join(r["command"]))
for i in r["inputs"]:
    print(i["option"], i["sha256"])
for name, t in r["tools"].items():
    print(name, t["version"])
EOF
```

```
completed 0 0.6.0.dev0 2026-11-14T14:57:54Z
coal pipeline pairs --query coal-tutorial-data/ecoli_K12_MG1655.cds.fna.gz --subject coal-tutorial-data/styphimurium_LT2.cds.fna.gz -o ecoli_vs_salmonella -j 1
--query 8b7b48b7b905313082141db4b77807b7aafebca7a5c40bbd42832baacd1cfd9b
--subject 111a4660298f1d1edd99e5e7171538dd805f36606aab1ac9aa0bc0ac986ee0e5
coal-cpp 0.6.0.dev0
diamond 2.2.8
```

The input checksums are those of the download script, so anyone can check they have the same files. `parameters`
has every option with the value used (defaults included), `tools.*.commands` the exact DIAMOND and engine
command lines, `platform` and `threads` the machine. Running the same command again gives the same bytes:

```
coal pipeline pairs --query coal-tutorial-data/ecoli_K12_MG1655.cds.fna.gz \
    --subject coal-tutorial-data/styphimurium_LT2.cds.fna.gz -o ecsa_rerun -j 1 -q
sha256sum ecoli_vs_salmonella/dnds.tsv ecsa_rerun/dnds.tsv
```

```
7ed75371df946874b4270a0c17e08e496ba75f04ba3ef5d91678a11092e1b8e8  ecoli_vs_salmonella/dnds.tsv
7ed75371df946874b4270a0c17e08e496ba75f04ba3ef5d91678a11092e1b8e8  ecsa_rerun/dnds.tsv
```

Rows never depend on `-j` either. With the same COAL and DIAMOND versions and inputs you should get this exact
checksum; another DIAMOND version can find slightly different hits.

**Change the method without searching again:** `--pairs ecoli_vs_salmonella/orthologs.tsv` skips DIAMOND and
uses those pairs; `--resume` in the *same* output folder reuses the translation and the search when the inputs
and search options are unchanged (the report's `stages` then say `cached: true`). Adding a method changes
neither:

```
coal pipeline pairs --query coal-tutorial-data/ecoli_K12_MG1655.cds.fna.gz \
    --subject coal-tutorial-data/styphimurium_LT2.cds.fna.gz -o ecsa_rerun -j 1 --resume -m YN -m NG
```

```
coal pipeline: translate_query: reused from the earlier run (--resume)
coal pipeline: translate_subject: reused from the earlier run (--resume)
...
coal pipeline: diamond_search: reused from the earlier run (--resume)
coal pipeline: 3,142 pairs (rbh)
coal pipeline: summary: 8,872 records (0 not usable); 3,142 pairs: NG 2,924 ok, 3 warning, 0 not_converged, 215 saturated, 0 no_data, 0 error; YN 2,973 ok, 3 warning, 11 not_converged, 155 saturated, 0 no_data, 0 error
coal pipeline: summary: decompress_inputs 0.0 s, translate_query 0.0 s (reused), translate_subject 0.0 s (reused), diamond_search 0.0 s (reused), orthologs 0.0 s, align_codon_kaks 7.0 s
```

(`--resume` reuses `work/` of the same folder: a copied or renamed folder starts again.) Which method to pick,
and why NG saturates more often here: [tutorial 5](05-choosing-a-method.md).

**What to cite** is in the report too:

```
coal cite ecoli_vs_salmonella/report.json
```

```
[coal] Newton C, Murati K (2026). COAL: codon alignments from protein alignments, version 0.6.0.dev0. https://github.com/plumat-coal/COAL
[suyama2006] Suyama M, Torrents D, Bork P (2006). PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments. Nucleic Acids Research 34:W609-W612. https://doi.org/10.1093/nar/gkl315
[yang2000] Yang Z, Nielsen R (2000). Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. Molecular Biology and Evolution 17:32-43. https://doi.org/10.1093/oxfordjournals.molbev.a026236
[wang2010] Wang D, Zhang Y, Zhang Z, Zhu J, Yu J (2010). KaKs_Calculator 2.0: a toolkit incorporating gamma-series methods and sliding window strategies. Genomics, Proteomics & Bioinformatics 8:77-80. https://doi.org/10.1016/S1672-0229(10)60008-3
[yang2007] Yang Z (2007). PAML 4: phylogenetic analysis by maximum likelihood. Molecular Biology and Evolution 24:1586-1591. https://doi.org/10.1093/molbev/msm088
[buchfink2021] Buchfink B, Reuter K, Drost H-G (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods 18:366-368. https://doi.org/10.1038/s41592-021-01101-x
```

## A methods paragraph

Built from the report, for this run:

> Coding sequences of *E. coli* K-12 MG1655 (GCF_000005845.2) and *S. enterica* Typhimurium LT2
> (GCF_000006945.2) were compared with `coal pipeline pairs` (COAL 0.6.0.dev0). Orthologs were reciprocal best
> hits of DIAMOND 2.2.8 blastp (e-value 1e-5; Buchfink et al. 2021); each pair's proteins were aligned globally
> (BLOSUM62, gap open 10, extend 0.5), the CDS threaded onto the alignment (Suyama et al. 2006), and dN and dS
> estimated with YN00 (Yang & Nielsen 2000) as in KaKs_Calculator 2.0 (Wang et al. 2010). Of 3,142 pairs, we kept
> those with an estimate, at least 50 compared codons covering half the alignment, protein identity ≥ 40% and
> 0.01 ≤ dS ≤ 2 (2,468 pairs).

## Clean up

Tutorials 4 and 5 use `ecoli_vs_salmonella/` and `ecsa_filtered/`; once you are done with them:

```
rm -rf ecoli_vs_salmonella ecsa_filtered ecsa_rerun kept.tsv kept_awk.tsv ribo.ids
```

## Reference

Sharp PM, Li WH (1987). The rate of synonymous substitution in enterobacterial genes is inversely related to
codon usage bias. *Mol Biol Evol* 4:222-230.
