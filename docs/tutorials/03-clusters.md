# Tutorial 3: gene families across several genomes (`coal pipeline clusters`)

**Goal:** cluster the genes of four genomes into families, get dN/dS for every pair within each family and a
per-family summary, and learn what decides the result: how tight the clusters are, which aligner builds the
family alignments, and what happens to very large families.

| | |
|---|---|
| Data | four *Chlamydia* genomes (3,738 CDS): *C. trachomatis* D/UW-3/CX and L2/434/Bu (two strains of one species), *C. muridarum* Nigg (its sister species, from mouse) and *C. pneumoniae* CWL029 (distant), NCBI RefSeq |
| Needs | the platform wheel (`coal-cpp`), DIAMOND 2.1.0 or later, Clustal Omega (`coal pipeline install diamond clustalo`); MAFFT for one step |
| Runtime | about 85 s per run with one thread (DeepClust 10 s, 814 Clustal Omega alignments 58 s, Ka/Ks the rest); the MAFFT run 6.3 min |
| Tested with | COAL 0.6.0.dev0 (commit 009bdd5), DIAMOND 2.2.8, Clustal Omega 1.2.4, MAFFT 7.526, Linux x86_64 |

Four genomes at three distances, on purpose: two strains that differ at about 1% of synonymous sites, a sister
species at about 1.3 substitutions per synonymous site, and an outgroup beyond what pairwise dS can measure.
Reading the result means telling these apart.

## 0. Set up

```
export T=/path/to/COAL/docs/tutorials
cd coal-tutorials
bash $T/data/get_chlamydia.sh
D=coal-tutorial-data
CDS="$D/ctrachomatis_D.cds.fna.gz $D/ctrachomatis_L2.cds.fna.gz $D/cmuridarum_Nigg.cds.fna.gz $D/cpneumoniae_CWL029.cds.fna.gz"
```

```
ok                 coal-tutorial-data/ctrachomatis_D.cds.fna.gz
ok                 coal-tutorial-data/ctrachomatis_L2.cds.fna.gz
ok                 coal-tutorial-data/cmuridarum_Nigg.cds.fna.gz
ok                 coal-tutorial-data/cpneumoniae_CWL029.cds.fna.gz
```

(GCF_000008725.1, GCF_000068585.1, GCF_000006685.1, GCF_000008745.1; checksums in the script.) NCBI's record
IDs (`lcl|NC_000117.1_cds_NP_219506.1_5`) start with the molecule's accession, so they are unique across the four
files. With files whose IDs repeat across species (`gene1`, `g1.t1` in each), COAL refuses to start and
`--prefix-ids` puts each file's name in front of its IDs.

## 1. Run it

```
coal pipeline clusters --cds $CDS -o chlamydia -j 1
```

```
coal pipeline: 4 gzip-compressed input file(s) decompressed into .../chlamydia/work/inputs
coal pipeline: translated 3,738 records (0 not usable)
coal pipeline: 1,074 clusters, 3,738 sequences
coal pipeline: warning: 3,395 of the 6,873 pairs (49%) are saturated or not_converged: dS (sometimes dN too) has no estimate there: these sequences are too divergent for synonymous changes to be counted (dN and the protein identity are still given). Compare closer species for dS; docs/pipeline.md "Filtering before interpretation"
coal pipeline: summary: 3,738 records (0 not usable); 6,873 pairs: 3,386 ok, 34 warning, 105 not_converged, 3,290 saturated, 58 no_data, 0 error
coal pipeline: summary: 1,074 clusters: 853 with every pair, 221 without a pair, 0 above --max-cluster-size (0 skipped)
coal pipeline: summary: decompress_inputs 0.0 s, translate_proteins 0.1 s, deepclust 10.1 s, align_codon_kaks 72.6 s, cluster_summary 0.1 s
coal pipeline: summary: results in .../chlamydia/dnds.tsv
coal pipeline: summary: what to cite for this run: coal cite .../chlamydia/report.json
```

What happened: all CDS were translated, DIAMOND DeepClust clustered the proteins with its defaults (member
coverage 80%, no identity threshold), each cluster of 3 or more was aligned by Clustal Omega, clusters of 2 by
COAL's pairwise aligner, the CDS were threaded onto those alignments, and Ka/Ks (YN00) was computed for every
pair of rows in every cluster.

**Half the pairs have no dS, and the warning says why.** That is not a failed run: it is what default clusters
over these four genomes give. Steps 3 and 4 show where those pairs come from and how to get fewer of them.

## 2. The per-cluster summary

```
ls chlamydia
cut -f3 chlamydia/cluster_summary.tsv | tail -n +2 | sort -n | uniq -c
```

```
cluster_summary.tsv  clusters.tsv  dnds.tsv  failed.tsv  genes.tsv  report.json  work
    221 1
     39 2
     75 3
    685 4
      4 5
      9 6
      5 7
     26 8
      5 12
      2 13
      1 15
      1 16
      1 40
```

685 of the 1,074 clusters have exactly 4 members: one gene per genome, the shape of a family of orthologs. The 221
singletons are genes found in one genome only (or too different to cluster). Clusters of 8, 12 or 16 are
families with paralogs, and one cluster has 40 members.

The first cluster, column by column (abridged):

```
head -2 chlamydia/cluster_summary.tsv | cut -f3-13,20,26,33
```

```
n_members	n_used	method	pair_selection	aligner	n_pairs	n_ok	n_not_converged	n_saturated	n_no_data	n_failed	dS_median	dNdS_median	status
4	4	YN	all	clustalo	6	3	0	3	0	0	1.4807088	0.011209224	ok
```

How to read it: 4 genes, all usable, all 6 pairs computed from one Clustal Omega alignment. 3 pairs have
estimates, 3 are saturated (the three pairs with *C. pneumoniae*). The medians are over the 3 pairs with an
estimate, which here are one pair at dS 0 (the two *C. trachomatis* strains, identical at synonymous sites, so
dN/dS is `NA` for it) and two at dS about 1.5. **A cluster median mixes pairs at very different distances**; it is
a quick screen, not a result. `status` is `ok` when at least one pair has an estimate; the full rules and the
other columns are in [pipeline.md](../pipeline.md#cluster_summarytsv).

## 3. Split the pairs by species before reading them

Every row of `dnds.tsv` has `cluster_id`, `seq1` and `seq2`. The summary script groups the rows by the files the
two genes come from (`--by-species`):

```
python3 $T/scripts/dnds_summary.py chlamydia/dnds.tsv --cds $CDS --by-species --top 5
```

```
chlamydia/dnds.tsv: 6873 rows of method YN
status: ok 3386, saturated 3290, not_converged 105, no_data 58, warning 34
by species pair: rows, with an estimate, saturated + not_converged, median dS, median dN/dS (rows with an estimate)
  cmuridarum_Nigg x cmuridarum_Nigg: 94, 24, 68, 1.8161, 0.4946
  cmuridarum_Nigg x cpneumoniae_CWL029: 1032, 165, 865, 2.4003, 0.1371
  cmuridarum_Nigg x ctrachomatis_D: 1027, 864, 160, 1.3204, 0.0531
  cmuridarum_Nigg x ctrachomatis_L2: 1020, 870, 146, 1.3121, 0.0517
  cpneumoniae_CWL029 x cpneumoniae_CWL029: 443, 180, 230, 1.8904, 0.4278
  cpneumoniae_CWL029 x ctrachomatis_D: 1029, 195, 832, 2.3080, 0.1594
  cpneumoniae_CWL029 x ctrachomatis_L2: 1033, 187, 844, 2.2758, 0.1549
  ctrachomatis_D x ctrachomatis_D: 86, 16, 66, 2.1551, 0.4326
  ctrachomatis_D x ctrachomatis_L2: 1029, 902, 124, 0.0098, 0.2886
  ctrachomatis_L2 x ctrachomatis_L2: 80, 17, 60, 2.2559, 0.3944
filters (each applied to what the previous ones kept):
  has an estimate (status ok or warning)           3420  (-3453)
  codons_compared >= 50                            3399  (-21)
  codons_compared / prot_aln_length >= 0.5         3275  (-124)
  prot_identity >= 40                              2924  (-351)
  dS <= 2                                          2445  (-479)
  dS >= 0.01                                       1987  (-458)
...
```

This table is the key to the run. Reading it line by line:

- ***C. trachomatis* D x L2** (1,029 pairs, nearly all estimated): median dS 0.0098. The two strains differ at
  one synonymous site in a hundred, so most pairs rest on a handful of synonymous differences. Their median
  dN/dS, 0.29, is five times that of the species comparison below. This is the well-known **time dependence of
  dN/dS** between close bacterial genomes (Rocha et al. 2006): slightly deleterious nonsynonymous changes have
  not yet been removed by selection. It is a property of the comparison, not a sign of positive selection, and
  the `dS >= 0.01` filter removes about half of these pairs for that reason (all 458 pairs it removes are D x L2
  pairs).
- ***C. muridarum* x *C. trachomatis*** (about 1,020 pairs each): median dS 1.3, median dN/dS 0.05. This is the
  informative comparison: far enough apart for selection to have acted, close enough to measure (84% of the
  pairs have an estimate).
- ***C. pneumoniae* x anything**: more than 80% of the pairs are saturated and the rest have median dS 2.3. Pairwise dS
  cannot be measured at this distance. dN still can (it is in the table for saturated pairs too), and so can
  protein identity; for dN/dS at this depth, use codon models on a tree (codeml, HyPhy) with the codon alignments
  of `--keep-alignments` as input.
- **Within one genome** (`x` the same species): paralog pairs, mostly ancient duplications, mostly saturated.
  Their dN/dS values (0.4 to 0.5) come from the few that are not, and mean little.

To keep only one species pair, filter on the accessions at the start of the IDs, for example *C. trachomatis* D
(`NC_000117`) against *C. muridarum* (`NC_002620`):

```
awk -F'\t' 'NR==1 || ($2 ~ /NC_000117/ && $3 ~ /NC_002620/) || ($2 ~ /NC_002620/ && $3 ~ /NC_000117/)' \
    chlamydia/dnds.tsv > ct_vs_cm.tsv
python3 $T/scripts/dnds_summary.py ct_vs_cm.tsv --cds $CDS --top 3
```

```
ct_vs_cm.tsv: 1026 rows of method YN
status: ok 860, saturated 145, not_converged 14, warning 4, no_data 3
...
quartiles over the 711 rows kept (25% / median / 75%):
  dN    0.0366 / 0.0629 / 0.1012
  dS    0.8945 / 1.2427 / 1.4405
  dNdS  0.0331 / 0.0538 / 0.0918
  dN/dS > 1: 0 of 711
highest dN/dS among the rows kept:
  CT_081/TC_RS04930 dN 0.3712  dS 0.4436  dN/dS 0.8367  codons   97  hypothetical protein
  CT_849.1/TC_RS01200 dN 0.0844  dS 0.1949  dN/dS 0.4332  codons   62  hypothetical protein
  omcA/TC_RS03685 dN 0.0821  dS 0.2021  dN/dS 0.4063  codons   88  lipoprotein OmcA
...
```

711 ortholog pairs pass the filters; none has dN/dS above 1, and the highest are short, unnamed or
surface-exposed proteins (OmcA, an outer-membrane lipoprotein), the usual candidates for faster protein
evolution. The names in the list are `gene/locus tag` of the two genomes when they differ.

## 4. Tighter clusters

DeepClust's defaults put distant paralogs in one cluster ([pipeline.md](../pipeline.md#cluster-mode-input-order-and-what-a-default-cluster-is)).
For dN/dS within families of orthologs or recent duplicates, ask for more identity and coverage of both
sequences:

```
coal pipeline clusters --cds $CDS -o chlamydia_tight --approx-id 60 --mutual-cover 80 -j 1
```

```
coal pipeline: 1,526 clusters, 3,738 sequences
coal pipeline: warning: 1,325 of the 4,074 pairs (33%) are saturated or not_converged: ...
coal pipeline: summary: 3,738 records (0 not usable); 4,074 pairs: 2,731 ok, 18 warning, 27 not_converged, 1,298 saturated, 0 no_data, 0 error
coal pipeline: summary: 1,526 clusters: 879 with every pair, 647 without a pair, 0 above --max-cluster-size (0 skipped)
coal pipeline: summary: decompress_inputs 0.0 s, translate_proteins 0.1 s, deepclust 10.3 s, align_codon_kaks 65.5 s, cluster_summary 0.0 s
```

```
python3 $T/scripts/dnds_summary.py chlamydia_tight/dnds.tsv --cds $CDS --by-species --top 0
```

```
by species pair: rows, with an estimate, saturated + not_converged, median dS, median dN/dS (rows with an estimate)
  cmuridarum_Nigg x cmuridarum_Nigg: 2, 2, 0, 1.0058, 0.4891
  cmuridarum_Nigg x cpneumoniae_CWL029: 528, 88, 440, 2.3951, 0.0629
  cmuridarum_Nigg x ctrachomatis_D: 806, 797, 9, 1.2964, 0.0495
  cmuridarum_Nigg x ctrachomatis_L2: 814, 807, 7, 1.2970, 0.0490
  cpneumoniae_CWL029 x cpneumoniae_CWL029: 3, 3, 0, 0.2659, 0.5265
  cpneumoniae_CWL029 x ctrachomatis_D: 527, 94, 433, 2.2652, 0.0636
  cpneumoniae_CWL029 x ctrachomatis_L2: 528, 94, 434, 2.2587, 0.0591
  ctrachomatis_D x ctrachomatis_D: 1, 1, 0, 2.2906, 0.1422
  ctrachomatis_D x ctrachomatis_L2: 864, 862, 2, 0.0092, 0.2731
  ctrachomatis_L2 x ctrachomatis_L2: 1, 1, 0, 2.2040, 0.1503
```

The within-genome paralog pairs are gone (from 703 to 7), about half the *C. pneumoniae* genes (515 of the 647
singletons) now form clusters of their own, the *C. muridarum* x *C. trachomatis* comparison lost 20% of its pairs but nearly
all its saturated ones (9 instead of 160), and the `no_data` pairs are gone (below). What you trade: fewer genes
in families. Which settings are right depends on the question; look at `n_saturated` and the dS distribution,
as here, and say in the methods what you used (`report.json` records it).

## 5. Unusual clusters: `no_data`, and a family of 40

**`no_data` pairs share no codon.** The 58 `no_data` pairs of the default run come from a few clusters. One of
them:

```
grep NP_219847.1_344 chlamydia/cluster_summary.tsv | cut -f3-13
```

```
12	12	YN	all	clustalo	66	9	0	41	16	0
```

Its representative is a 678-residue fused alpha/beta subunit of a 2-oxoacid dehydrogenase E1 component; its
members include the separate alpha (340 residues) and beta (328 residues) subunits of pyruvate dehydrogenase
from each genome. Each half matches its part of the fused protein well enough to join the cluster, but an alpha
and a beta subunit have no column in common in the alignment, so their pairs have `codons_compared` 0 and status
`no_data`. Such domain-fusion clusters are one more reason to tighten the clustering (`--mutual-cover` needs
both sequences covered) or to read `n_no_data`.

**Large clusters.** The number of pairs grows with the square of the cluster size. The 40-member cluster here is
the polymorphic membrane protein (Pmp) autotransporter family (19 genes in *C. pneumoniae*, 7 in each of the
others): 780 pairs, 672 of them saturated. Clusters above `--max-cluster-size` (default 200) are not aligned as
one MSA; `--large-clusters` decides which pairs are computed (`representative`, the default: n − 1 pairs;
`sample`; `skip`). With a cap of 20, reusing the clusters of the first run (`--clusters` skips DeepClust):

```
coal pipeline clusters --cds $CDS --clusters chlamydia/clusters.tsv -o chlamydia_cap --max-cluster-size 20 -j 1
awk -F'\t' '$3==40' chlamydia_cap/cluster_summary.tsv | cut -f3-13,33,34
```

```
40	40	YN	representative	nw_pairwise	39	1	0	38	0	0	ok	40 members: above --max-cluster-size 20
```

The representative against each of the other 39 members, aligned pairwise; the note says why. With
`--large-clusters sample --pair-sample K --seed N` the draw is reproducible.

## 6. Picking the aligner

Clusters of 3 or more are aligned by Clustal Omega by default; `--aligner mafft` uses MAFFT (`mafft --auto`,
which picks its accurate L-INS-i mode for small families). Same clusters, other aligner:

```
coal pipeline clusters --cds $CDS --clusters chlamydia/clusters.tsv -o chlamydia_mafft --aligner mafft -j 1
```

```
coal pipeline: summary: 3,738 records (0 not usable); 6,873 pairs: 3,377 ok, 33 warning, 82 not_converged, 3,358 saturated, 23 no_data, 0 error
coal pipeline: summary: decompress_inputs 0.0 s, translate_proteins 0.1 s, clusters 0.0 s, align_codon_kaks 377.3 s, cluster_summary 0.0 s
```

MAFFT took 359 s for the 814 alignments, Clustal Omega 58 s. Compare the two `dnds.tsv` files pair by pair,
by protein identity:

```
python3 - <<'EOF'
import csv
def load(d):
    return {(r["seq1"], r["seq2"]): r for r in csv.DictReader(open(d + "/dnds.tsv"), delimiter="\t")}
a, b = load("chlamydia"), load("chlamydia_mafft")
est = lambda r: r["status"] in ("ok", "warning")
both = [k for k in a if est(a[k]) and est(b[k])]
print(len(both), "pairs with an estimate in both")
for lo, hi in ((90, 101), (70, 90), (50, 70), (0, 50)):
    ks = [k for k in both if lo <= float(a[k]["prot_identity"]) < hi]
    diff = sum(1 for k in ks if a[k]["dS"] != b[k]["dS"])
    big = sum(1 for k in ks if abs(float(a[k]["dS"]) - float(b[k]["dS"])) > 0.1 * float(a[k]["dS"]))
    print(f"identity {lo}-{hi}: {len(ks)} pairs, dS not identical {diff}, more than 10% apart {big}")
EOF
```

```
3262 pairs with an estimate in both
identity 90-101: 1488 pairs, dS not identical 49, more than 10% apart 14
identity 70-90: 1092 pairs, dS not identical 276, more than 10% apart 75
identity 50-70: 242 pairs, dS not identical 161, more than 10% apart 61
identity 0-50: 440 pairs, dS not identical 412, more than 10% apart 223
```

In a table:

| Protein identity | Pairs | dS not identical | dS differs by more than 10% |
|---|---|---|---|
| 90% or more | 1,488 | 49 | 14 |
| 70-90% | 1,092 | 276 | 75 |
| 50-70% | 242 | 161 | 61 |
| below 50% | 440 | 412 | 223 |

Between close sequences the aligners agree (an alignment with few gaps has one answer). Below 50% identity they
almost never do, and a quarter of the pairs move by more than 10%: there the alignment, not the method, is the
largest uncertainty. That is one more reason to filter on identity and to name the aligner in the methods
(`report.json`: `parameters.aligner`, the tool's version and command template). `--keep-protein-alignments` and
`--keep-alignments` write each cluster's alignments to look at. For timing, `-j N` runs N alignments at a time;
[alignment.md](../alignment.md) compares the aligners on other data.

## Clean up

Tutorial 4 uses `chlamydia/`; once you are done with it:

```
rm -rf chlamydia chlamydia_tight chlamydia_cap chlamydia_mafft ct_vs_cm.tsv
```

## What to cite

`coal cite chlamydia/report.json` lists COAL, PAL2NAL, YN00 / KaKs_Calculator 2.0 / PAML, DIAMOND (Buchfink et
al. 2021), DeepClust (Buchfink et al. 2026) and Clustal Omega (Sievers et al. 2011); with `--aligner mafft`,
MAFFT (Katoh & Standley 2013) instead.

## Reference

Rocha EPC, Smith JM, Hurst LD, Holden MTG, Cooper JE, Smith NH, Feil EJ (2006). Comparisons of dN/dS are time
dependent for closely related bacterial genomes. *J Theor Biol* 239:226-235.
