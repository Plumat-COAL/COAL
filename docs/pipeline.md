# One-step dN/dS: `coal pipeline`

`coal pipeline` takes coding sequences (CDS, FASTA) and writes a table of dN, dS and dN/dS. It runs every step in
between: translation, the DIAMOND search or clustering, the protein alignment, the codon alignment and Ka/Ks. It
has two modes:

```
coal pipeline pairs    --query q_cds.fa --subject s_cds.fa -o out_dir      # orthologs between two species
coal pipeline clusters --cds all_cds.fa [--clusters clusters.tsv] -o out_dir   # every pair within each cluster
coal pipeline install diamond clustalo                                      # the external tools, into COAL's cache (verified)
```

| Step | Pair mode | Cluster mode | Done by |
|---|---|---|---|
| 1. Translate | both files | all files | `coal-cpp translate` |
| 2. Find the groups | DIAMOND blastp both ways, reciprocal best hits (`--ortho bh`: best hits); or `--pairs TSV` | DIAMOND DeepClust; or `--clusters TSV` | DIAMOND |
| 3. Align the proteins | global NW, BLOSUM62 (`--aligner clustalo`, `mafft`, `diamond`) | Clustal Omega per cluster of 3 or more (`--aligner mafft`); NW for 2 | `coal-cpp codon-dnds`; Clustal Omega / MAFFT |
| 4. Codon alignment | per pair | per cluster | `coal-cpp codon-dnds` |
| 5. Ka/Ks | per pair | every pair of the cluster | `coal-cpp codon-dnds` |
| 6. Tables | `dnds.tsv` | `dnds.tsv`, `cluster_summary.tsv` | Python |

The Python module `coal.dnds_pipeline` only orchestrates. It needs `coal-cpp` (the platform wheels have it;
`coal --cpp-path`), and nothing outside the engine package: it works without `coal_app`. DIAMOND, Clustal Omega
and MAFFT are separate programs. They are looked for on `PATH`, in `$COAL_DIAMOND` / `$COAL_CLUSTALO` /
`$COAL_MAFFT`, and in COAL's cache; `--install-tools` (or `coal pipeline install ...`) installs a missing one into
the cache. They are not bundled (licences: [alignment.md](alignment.md#licences-why-clustal-omega-is-not-bundled)).

## Pair mode

```
coal pipeline pairs --query thaliana_cds.fa --subject lyrata_cds.fa -o out -j 8
coal pipeline pairs --query q.fa --subject s.fa -o out -m YN -m NG --kaks-format --keep-alignments
coal pipeline pairs --query q.fa --subject s.fa -o out --pairs my_orthologs.tsv        # no DIAMOND
```

- **Orthologs.** `diamond blastp --max-target-seqs 1 --evalue 1e-5`, query against subject and subject against
  query. A pair is kept when each is the other's best hit (`--ortho rbh`, the default). `--ortho bh` keeps every
  query's best hit, so a subject can appear in several pairs. `--evalue`, `--sensitivity` (DIAMOND's modes, from
  `faster` to `ultra-sensitive`) and `--diamond-args '...'` are passed on. Pairs are written in the order of the
  query file.
- **`--pairs TSV`** takes the pairs from a file (query ID, subject ID; a header line `query subject` is skipped)
  and DIAMOND is not run. `out/orthologs.tsv` of an earlier run can be given.
- **Alignment (`--aligner`).**
  - `nw` (default): COAL's global aligner, BLOSUM62, gap open 10, extend 0.5, end gaps free, full length. Change
    it with `--matrix`, `--gap-open`, `--gap-extend`, `--end-gaps`, or `--preset pwalign`.
  - `diamond`: DIAMOND's own alignment of the hit, completed to full length (`--flanks align`, the default, aligns
    the ends DIAMOND left out; `gap` puts them against gaps). Nothing is realigned. The search then runs with
    `--masking 0`, which can change which pairs are found.
  - `clustalo`, `mafft`: one process per pair, `-j` at a time. Slow, for comparison only.
- A gene without a partner is listed in `out/unpaired.tsv` (`no_hit` or `not_reciprocal`).

## Cluster mode

```
coal pipeline clusters --cds all_cds.fa -o out -j 8                          # DeepClust, then every pair
coal pipeline clusters --cds a.fa b.fa -o out --approx-id 50 --member-cover 80 --realign
coal pipeline clusters --cds all_cds.fa --clusters deepclust.tsv -o out --max-cluster-size 100 \
     --large-clusters sample --pair-sample 500 --seed 1
```

- **Clusters.** Without `--clusters`, `diamond deepclust` runs on the translated proteins with DIAMOND's defaults
  (member coverage 80%); `--approx-id`, `--member-cover`, `--mutual-cover`, `--cluster-evalue` and
  `--deepclust-args '...'` are passed on, and the result is kept as `out/clusters.tsv`. `--realign` also runs
  `diamond realign` (members against their representatives, the columns orthologr's `deepclust_realign()` asks
  for) into `out/realign.tsv`.
- **`--clusters TSV`** is DIAMOND's two-column table: representative, member, tab-separated, one line per member.
  A header line (`centroid member`, `representative member`, ...) and `#` lines are skipped. The cluster ID is the
  representative's ID. A member named in two clusters stays in its first one.
- **Alignment.** Clusters of 3 or more get one Clustal Omega run each (`--aligner mafft` for MAFFT), `-j`
  processes at a time, each with `--msa-threads` threads (default 1). Two-member clusters are aligned by COAL's
  global aligner. The codon MSA is threaded from the protein MSA, and Ka/Ks is computed for every pair of its rows.
- **Large clusters.** The number of pairs grows with the square of the cluster size (200 members: 19,900 pairs).
  A cluster with more than `--max-cluster-size` usable members (default 200) is handled as `--large-clusters` says:
  - `representative` (default): the representative against every other member, n − 1 pairs;
  - `sample`: `--pair-sample K` pairs (default 1,000) drawn without replacement. The draw depends only on `--seed`
    and the cluster ID, so it is the same in every run, for every `-j`, and whatever else is in the input;
  - `skip`: no pairs.

  The pairs of such a cluster are aligned pairwise by the global aligner, without an MSA, so the cost stays
  proportional to the number of pairs. `pair_selection` and `aligner` in the summary say what was done.
- **Singletons** and clusters with fewer than two usable members are listed in the summary with `n_pairs` 0.
- A cluster whose MSA run fails is not lost: its pairs are aligned pairwise instead, the summary's `aligner` is
  `nw_pairwise` and its note gives the reason.

## Options common to both modes

| Option | Meaning | Default |
|---|---|---|
| `-o DIR` | output folder | required |
| `-m METHOD` | Ka/Ks method, repeatable: YN, NG, LWL, MLWL, LPB, MLPB, MYN ([kaks.md](kaks.md)) | YN |
| `-c N` | NCBI genetic code, used for translation, the codon alignment and Ka/Ks | 1 |
| `-j N` | threads of DIAMOND and of the engine, and MSA processes at a time | one per CPU |
| `--kaks-format` | also write `dnds.kaks` in KaKs_Calculator's 22-column format | off |
| `--keep-alignments` | also write the codon alignments | off |
| `--resume` | reuse the translation, search and clustering results in `DIR/work` when the inputs (path, size, modification time), the options and the DIAMOND version are the same | off |
| `--install-tools` | install a missing DIAMOND / Clustal Omega / MAFFT into COAL's cache | off |
| `--diamond PATH` | the DIAMOND to run | found automatically |
| `--report FILE` | the JSON report | `DIR/report.json` |
| `--msa-timeout S` | give up on one Clustal Omega / MAFFT run after S seconds | none |
| `--yn-flavour`, `--yn-weighting`, `--digits`, `--max-cells` | as in `coal dnds` and `coal-cpp-align` | |

## Output

Everything is in the output folder. Tables are tab-separated with a header line, `NA` for a missing value, `\n`
line ends. Rows are in input order (pairs: the query file's order; clusters: the clusters file's order, pairs
1-2, 1-3, ..., 2-3, ...), whatever `-j` is.

| File | Mode | Content |
|---|---|---|
| `dnds.tsv` | both | one row per pair and method |
| `cluster_summary.tsv` | clusters | one row per cluster and method |
| `genes.tsv` | both | one row per CDS record: what translation found |
| `orthologs.tsv` | pairs | the pairs: `query subject` and, from DIAMOND, `pident length evalue bitscore` |
| `unpaired.tsv` | pairs | `set id reason`: usable genes without a partner (`no_hit`, `not_reciprocal`) |
| `clusters.tsv`, `realign.tsv` | clusters | DeepClust's clusters (when it ran); `diamond realign`'s table (`--realign`) |
| `failed.tsv` | both | everything that could not be used: `level cluster_id seq1 seq2 code message` |
| `report.json` | both | the run report (below) |
| `dnds.kaks` | both | `--kaks-format`: the same results as KaKs_Calculator writes them, pairs named `seq1&seq2` |
| `codon_alignments.fasta` | pairs | `--keep-alignments`: the codon alignments, records 1+2, 3+4, ... (`coal dnds --pairs consecutive` reads it) |
| `codon_alignments/<cluster>.fasta` | clusters | `--keep-alignments`: one codon MSA per cluster (not for clusters aligned pairwise); characters other than letters, digits, `.`, `_`, `-` in the cluster ID become `_` |
| `work/` | both | proteins, DIAMOND tables, logs of every external command (`work/logs/`), the `--resume` state |

### `dnds.tsv`

| Column | Meaning |
|---|---|
| `cluster_id` | cluster mode only (first column): the cluster's representative |
| `seq1`, `seq2` | the two genes. Pair mode: query, subject |
| `method` | NG, LWL, MLWL, LPB, MLPB, YN or MYN; with several methods a pair's rows come in this order |
| `dN`, `dS`, `dNdS` | Ka, Ks, Ka/Ks. `NA` when not available; `dNdS` is `NA` when dS < 1e-6 |
| `SE_dN`, `SE_dS` | standard errors (YN, MYN, NG); `NA` for the other methods |
| `N_sites`, `S_sites` | nonsynonymous and synonymous sites |
| `Nd`, `Sd` | nonsynonymous and synonymous differences |
| `codons_compared` | codon columns used (no gap, ambiguous base or stop in either sequence) |
| `prot_aln_length` | columns of the protein alignment with a residue in at least one of the two sequences |
| `prot_identity` | percent of those columns with the same residue in both (two decimals) |
| `status` | `ok`; `warning` (estimated, see note); `saturated` (dN or dS could not be corrected: `NA`); `no_data` (nothing to compare); `error` (the pair could not be computed: every value `NA`) |
| `note` | why the status is not `ok`: internal stop codons, a CDS length that is not a multiple of 3 (with `seq1:` / `seq2:`), the saturated estimate, or the error |

The numeric columns, statuses and notes are those of `coal dnds` ([kaks.md](kaks.md#output)); the pipeline adds
the two protein-alignment columns and what it knows about each CDS.

### `cluster_summary.tsv`

| Column | Meaning |
|---|---|
| `cluster_id`, `representative` | the representative's ID (both columns) |
| `n_members` | members in the clusters file |
| `n_used` | members with a usable CDS (the others are in `failed.tsv`) |
| `method` | one row per method; `NA` on the single row of a cluster without pairs |
| `pair_selection` | `all`, `representative`, `sample:K` (K pairs drawn), `skipped`, `singleton`, `none` (several members, fewer than two usable) |
| `aligner` | `clustalo` or `mafft` (MSA), `nw` (two members), `nw_pairwise` (above the size cap, or the MSA failed), `none` |
| `n_pairs` | rows in `dnds.tsv` for this cluster and method |
| `n_ok` | pairs with status `ok` or `warning` (both have estimates) |
| `n_saturated`, `n_no_data`, `n_failed` | pairs with status `saturated`, `no_data`, `error` |
| `dN_median`, `dN_mean`, `dN_q25`, `dN_q75`, `dN_min`, `dN_max` | over the `n_ok` pairs |
| `dS_...`, `dNdS_...` | the same six for dS and dN/dS |
| `n_dNdS` | how many dN/dS values the `dNdS_` columns are over (dN/dS is `NA` when dS is 0) |
| `status` | `ok`; `failed` (pairs, but none with an estimate); `skipped` |
| `note` | members left out, the size cap, an MSA that failed |

Quartiles are interpolated linearly between order statistics (R's type 7, numpy's default).

### `genes.tsv`

`set` (`query`, `subject` or `cds`), `id`, `nt_length`, `aa_length`, `terminal_stop` (1 when the last codon is a
stop; it is left out of the protein), `internal_stops`, `ambiguous_codons`, `trailing_nt` (bases after the last
complete codon), `status`, `code`, `note`.

| What is found | Status, code | What happens |
|---|---|---|
| a stop codon at the end | `ok` | dropped from the protein; not counted in Ka/Ks |
| lower case, RNA (`U`) | `ok` | read as upper-case DNA; the codon alignment keeps the letters as written |
| a codon with an ambiguous base (N, IUPAC) | `ok` | translated to the residue every reading gives, else X; that codon column is skipped in Ka/Ks |
| length not a multiple of 3 | `warning`, `cds_length` | the last 1 or 2 bases are ignored; the gene's pairs get status `warning` |
| an internal stop codon | `warning`, `internal_stop` | written as `*` and aligned; those codon columns are skipped; the gene's pairs get status `warning` |
| empty, shorter than one codon, only a stop codon | `failed`, `cds_empty` / `cds_too_short` / `cds_only_stop` | the gene is left out |
| a character that is no nucleotide code | `failed`, `cds_invalid_character` | the gene is left out |
| an ID that already occurred | `failed`, `duplicate_id` | the first record is used, this one is left out |

### `report.json`

Schema `coal.pipeline/1`: [pipeline.schema.json](pipeline.schema.json). It is written also when the run fails.

- `mode`, `command`, `out_dir`, `parameters` (every option as used), `tools` (path and version of coal-cpp,
  DIAMOND, Clustal Omega, MAFFT).
- `stages`: name, wall `seconds`, `cached` (reused by `--resume`). The stage `align_codon_kaks` also has
  `align_codon_seconds` and `kaks_seconds` (inside the engine) and, when an MSA program ran, `msa_seconds`.
- `counts`: `genes` (total, ok, warning, failed), `pairs` (total and one count per status), `orthologs` and
  `unpaired` (pair mode), `clusters` (total, singletons, all_pairs, above_max_size, skipped, msa_fallback,
  members_dropped).
- `failures`: the first 1,000 entries of `failed.tsv` (`failures_total`, `failures_truncated`), each with `level`
  (`gene`, `pair`, `cluster`), `code` and `message`.
- `exit_code`, and `error` when the run could not be done.

### Failures and exit codes

One bad gene, pair or cluster never stops a run.

- A gene that cannot be used is left out and listed (`genes.tsv`, `failed.tsv`). In pair mode a pair that names
  it (`--pairs`) gets a row with status `error`; in cluster mode it is dropped from its cluster.
- A pair that cannot be computed (alignment too large for `--max-cells`, the codon engine refused a sequence)
  gets a row with status `error` and the reason.
- `saturated`, `no_data` and `warning` are results, not failures.

| Exit code | Meaning |
|---|---|
| 0 | finished, nothing failed |
| 2 | finished, but `failed.tsv` has entries (as `coal --keep-going` and `coal-cpp-align`); also a command-line error |
| 1 | the run could not be done: an input cannot be read, a tool is missing or failed, the folder cannot be written |
| 130 | interrupted |

## Defaults, and why

- **YN00** is the Ka/Ks default, as in the rest of COAL ([kaks.md](kaks.md)).
- **Reciprocal best hits at e-value 1e-5** follow orthologr's `dNdS()` defaults.
- **NW with needle's defaults for pairs.** It is exact and takes about 0.4 s of CPU per 1,000 pairs. The DIAMOND
  route saves almost nothing and changes dN/dS (r = 0.988 against NW here, below); Clustal Omega and MAFFT per
  pair took about 80 and 500 times as long as NW here ([alignment.md](alignment.md#recommended-defaults)).
- **Clustal Omega for clusters**, as the master plan asks; MAFFT is the option.
- **`--max-cluster-size 200`** bounds one cluster at 19,900 pairs and its MSA at 200 rows. The default above it is
  `representative` because it needs no MSA and no random draw, and its cost is linear in the cluster size.
- **Full length.** Overhanging ends of the longer protein stay in the alignment as gaps. The final stop codon is
  dropped. (orthologr truncates the overhang and scores the stop-codon column.)

## How it runs

- **No process and no temporary file per pair.** The engine step is one `coal-cpp codon-dnds` process. Python
  writes the groups (IDs, CDS, and for clusters the protein MSA rows) to its standard input; the engine translates,
  aligns the pairs that have no rows, threads the codons and computes Ka/Ks on `-j` threads, and writes the tables.
  Cluster MSAs are one Clustal Omega process per cluster, `-j` at a time.
- **Bounded memory.** Sequences are fetched from the FASTA files by byte offset when they are needed. Python
  holds the IDs, the pair or cluster lists and one batch; the engine holds one batch of groups (at most 1,024
  groups or 16 Mb of CDS) and writes a large group's pairs 8,192 at a time. A cluster's own MSA, and one
  alignment table per thread (cells / 2 bytes; a 23,000-residue pair needs about 270 MB), are held whole.
- **Order.** Output order never depends on `-j`: the same bytes for `-j 1` and `-j 2` in every test here.
- **Subprocesses** are started without a shell (argument lists), so paths with spaces and Windows work the same.
- **`--resume`** works per stage (translation, search, clustering, realign). The engine step itself is always
  redone; it is the fast part in pair mode. In cluster mode the MSAs are redone too (see Limits).

The native steps can be run by hand:

```
coal translate cds.fa -o proteins.fa --genes genes.tsv          # = coal-cpp translate
coal codon-dnds -i groups.txt -o dnds.tsv -m YN -j 8            # = coal-cpp codon-dnds (--help: the input format)
```

### From Python

```python
from coal import dnds_pipeline as P
res = P.run_pairs("q_cds.fa", "s_cds.fa", "out", method=["YN"], jobs=8)           # the options, with underscores
res = P.run_clusters(["all_cds.fa"], "out", clusters="clusters.tsv", max_cluster_size=200)
res.exit_code; res.report["counts"]; res.outputs["dnds"]                          # PipelineResult
P.codon_dnds([("p1", [("a", cds_a), ("b", cds_b)])], "pairs.tsv")                 # groups held in Python
P.find_diamond(); P.install_diamond(); P.read_clusters(path); P.sample_pairs(n, k, seed, cluster_id)
```

`run_pairs` / `run_clusters` do not raise for a failed run; they return exit code 1 and the report. Like
`coal.align` and `coal.msa`, this module is not yet covered by the 1.0 promise of [python-api.md](python-api.md).

### From C++

`#include <coal/dnds.hpp>` (in `coal::core`): the same steps for sequences held in memory, without files.

```cpp
coal::dnds::Options opt;                       // methods {YN}, table 1, AlignParams (needle), flanks, kaks options
coal::dnds::Group g;                           // id, members {id, cds, row}, pairs (all / representative / list)
g.members = {{"a", cds_a, "", ""}, {"b", cds_b, "", ""}};
coal::dnds::GroupResult r = coal::dnds::run_group(g, opt);       // never throws for bad data
// r.pairs[k]: i, j, results (one kaks::Result per method), prot_aln_length, prot_identities, error, note
// r.members[k]: tr (Translation), codon_row, error; r.protein_rows; r.aligner; r.error
coal::dnds::Translation t = coal::dnds::translate_cds(cds, *coal::get_code(1));
// for threads: prepare_group(g, opt) once, then compute_pair(r, i, j, opt) from any thread
```

## Validation

All of this is from 2026-11-11, branch `feat/dnds-pipeline`, on a 2-core cloud VM, with DIAMOND 2.2.8, Clustal
Omega 1.2.4, MAFFT 7.526 and KaKs_Calculator 2.0.1. The scripts are in `benchmarks/pipeline/`.

**Pair mode, orthologr's 1,000-gene example** (`ortho_thal_cds_1000.fasta` against `ortho_lyra_cds_1000.fasta`,
defaults; `validate_pairs.py`):

- 2,000 genes translated (12 with a length that is not a multiple of 3), 999 reciprocal best hits, 2 genes
  `not_reciprocal`. 986 pairs `ok`, 12 `warning` (the 12 CDS lengths), 1 `saturated` (AT1G04625.1, as in
  [kaks.md](kaks.md)).
- **By hand.** The same pairs were translated by an independent Python translation, aligned with `coal-cpp-align`,
  threaded with one `coal-cpp` run per pair and scored with `coal-cpp dnds` and `coal-cpp kaks`. All 999 codon
  alignments are identical to the pipeline's, all 999 rows have the same text in every numeric column, and
  `dnds.kaks` is byte-identical.
- **KaKs_Calculator 2.0** `-m YN` on the same codon alignments (AXT): Ka, Ks and Ka/Ks agree in all 999 pairs
  within its 6 printed digits (largest relative difference 4.9e-6); the one exception is the saturated pair, where
  it prints a capped Ks of 3.68 and COAL gives `NA` (a documented difference, kaks.md).
- `-j 1` and `-j 2` give byte-identical tables.

**Edge cases** (`tests/tools/make_pipeline_fixture.py --big`, 22 query and 20 subject genes): the run finishes
with exit code 2 and reports each one.

| Case | Result |
|---|---|
| a 23,000-codon gene pair | `ok`, 22,982 codons compared, about 1.3 s for the alignment; equal to the by-hand route and to KaKs_Calculator |
| an internal stop codon | `warning`; the note names the gene and the skipped column |
| a CDS that is not a multiple of 3 | `warning`; the note gives the length |
| ambiguous bases (N, R, Y) | `ok`; the codons that hold them are skipped |
| identical sequences | `ok`, dN = dS = 0, dN/dS `NA`, identity 100 |
| lower case and `U` | `ok` |
| genes without a hit | `unpaired.tsv`, `no_hit` |
| a bad character, an empty record, 2 bases, a repeated ID | `failed` in `genes.tsv` and `failed.tsv`; the first record of the repeated ID is used |

**Cluster mode** (`validate_clusters.py`). Input: the two 1,000-gene files plus synthetic families of 6, 12, 2, 1,
30 and 300 members, 2,351 CDS; DeepClust with its defaults gave 844 clusters (27 singletons, 712 pairs, 104 with 3
to 30 members, the 300-member family above the size cap).

- **Per pair.** `coal-cpp dnds` on each of the 816 kept codon MSAs gives 2,933 rows; all are equal to `dnds.tsv`.
  The 299 representative pairs of the 300-member family were redone by hand (`coal-cpp-align`, `coal-cpp`,
  `coal-cpp dnds`): all equal. With `--large-clusters sample --pair-sample 500`: 300 of the 500 redone, all equal.
- **Summary.** All 844 rows recomputed from `dnds.tsv` with pandas and numpy (counts by status; median, mean,
  quartiles, minimum and maximum: 14,706 values): all equal. The 27 singletons are listed with `n_pairs` 0.
- **Reproducible.** The sample is the same for `-j 1` and `-j 2` and across runs (byte-identical `dnds.tsv` and
  summary), and different for another seed. DeepClust's clusters and the whole table are also identical for
  `-j 1` and `-j 2`.
- Of the 3,232 pairs, 576 are `saturated` and 8 `no_data`: DeepClust's default clusters hold distant paralogs
  (dS above 2). This is what the per-cluster `n_saturated` column is for.

**Tests.** `tests/test_dnds_pipeline.py` (33 tests) repeats these checks on small synthetic fixtures
(`tests/fixtures/pipeline/`): the by-hand comparison, the summary recomputed with the `statistics` module, the
sample's reproducibility, every edge case above, failures and exit codes, the report against its schema,
`--resume`, and a failing MSA program. Parts that need DIAMOND or Clustal Omega are skipped without them.

## Timings

Wall seconds on the 2-core VM (best of 3), orthologr's 1,000-gene example, 999 pairs, YN:

| Stage | `-j 1` | `-j 2` |
|---|---:|---:|
| translate (2 × 1,000 CDS) | 0.04 | 0.05 |
| DIAMOND: 2 × makedb + blastp | 2.18 | 1.40 |
| reciprocal best hits | 0.003 | 0.004 |
| engine step: alignment + codon alignment | 0.44 | 0.24 |
| engine step: Ka/Ks (YN00) | 1.50 | 0.75 |
| **whole run** (with Python start-up) | **4.2** | **2.5** |

- All seven methods instead of YN: the Ka/Ks part takes 1.8 s (`-j 2`).
- `--aligner diamond`: 0.04 s for the alignment step instead of 0.24 s. It finds 1,000 pairs, not 999, because
  of `--masking 0`, and dN/dS correlates with the default's at r = 0.988 (619 of 998 common pairs identical).
- `--aligner clustalo`: 19 s for the alignment (r = 0.986). `--aligner mafft`: 123 s (r = 0.995).
- **19,980 pairs** (the same set 20 times, `--pairs`): 38.3 s with `-j 1`, 20.5 s with `-j 2`; peak memory 41 MB
  for Python and 41 MB for coal-cpp.
- **Cluster mode** on the 2,351-CDS set (844 clusters, 3,232 pairs), `-j 2`: 16.4 s. DeepClust 5.0 s, Clustal
  Omega 7.1 s (104 clusters), Ka/Ks 3.9 s, alignment + codon 0.2 s, summary 0.03 s. With `-j 1`: 28.6 s.

For scale, orthologr 0.4.3's `dNdS()` spent 377 s in `pairwise_aln` alone for the same 1,000 genes in the 29–30
October tests. That was another machine, so it is not a like-for-like benchmark; the 3-way comparison on one
machine is still to be done.

## Limits

- **`--resume` does not cover the engine step or the cluster MSAs.** A cluster run that is interrupted during the
  MSAs repeats them.
- **One MSA per cluster is held in memory**, and Clustal Omega's time grows faster than linearly with the cluster
  size. Lower `--max-cluster-size`, or use `--aligner mafft`, for inputs with many clusters of hundreds of members.
- **Large clusters are aligned pairwise**, so their values come from pairwise alignments and those of smaller
  clusters from an MSA. The summary says which (`aligner`).
- **Pair alignment memory** is cells / 2 bytes per thread (`--max-cells`, default 4·10⁹ cells): a titin-sized pair
  needs about 600 MB. Above the limit that pair gets status `error`.
- **`--aligner diamond` uses one HSP** and DIAMOND's heuristic alignment ([alignment.md](alignment.md)); a hit that
  does not fit the translated sequences is that pair's error.
- **Best hits are DIAMOND's first target per query** (`--max-target-seqs 1`). Ties are broken by DIAMOND.
- **Protein input is not accepted**: the pipeline starts from CDS, so that protein and CDS always match.
- **Python holds every ID** (and a byte offset per record), so memory grows with the number of sequences, not
  with their length.
- **Installers verify what they download** ([alignment.md](alignment.md#clusters-clustal-omega-and-mafft-coalmsa),
  "Verification"). The DIAMOND 2.2.8 release archives (Linux x86_64, macOS, Windows x64) and micromamba are pinned
  by SHA-256 in COAL's source; DIAMOND publishes no checksums, so those hashes are the ones of the files as
  downloaded on 2026-11-12. The official Clustal Omega binary has no known hash and is refused unless
  `coal pipeline install clustalo` gets `--sha256 HEX` or `--allow-unverified`; `--install-tools` during a run
  never installs an unverified file. Elsewhere Bioconda is used.
- Thread scaling beyond 2 cores and whole-proteome wall times are to be measured on the 16-thread desktop.
