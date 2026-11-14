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

## One record is one gene

**The pipeline treats every FASTA record as a gene.** It does not know about transcripts. CDS files from genome
annotations (Ensembl, NCBI, TAIR, Phytozome, FlyBase, a MAKER or Trinity assembly) usually hold every transcript of every gene, and
with such a file a gene with five isoforms is five "genes": it is compared once per transcript, the reciprocal
best hits are split between the isoforms of a gene (most transcripts end up in `unpaired.tsv` as
`not_reciprocal`), and the same gene pair appears several times with different values. On an Ensembl human / mouse
subset (36,827 and 9,931 transcripts) an independent review found 3,967 pairs for only 2,854 gene pairs, 685 gene
pairs appearing 2 to 18 times, and for 202 of them a spread of more than 1.5-fold in dN/dS between the isoform
pairs. **Give one sequence per gene**, or let COAL choose one:

```
coal pipeline pairs --query q_cds.fa --subject s_cds.fa -o out --longest-isoform
coal pipeline pairs --query q_cds.fa --subject s_cds.fa -o out --longest-isoform --gene-map genes.tsv
```

- **A warning by default.** When the records of a file look like several transcripts per gene, the run prints
  one warning with the counts (`query: 36,827 records of 3,734 genes, ...`), an example and where the genes were
  read from, and goes on as before: nothing else changes without the options below.
- **Where the gene of a record comes from.** (1) `--gene-map TSV`: record ID and gene ID, tab-separated (a header
  line and `#` lines are skipped; one map serves all input files). The map wins for every record it names. A
  record it does not name falls back to (2) and (3), read for its file as without a map, so a map for one
  species, or for part of a file, does not switch the detection off for the rest (`--isoform-rule map` does: then
  a record the map does not name is a gene of its own). A map that names no record of any input file is an error
  (exit 1; the message shows IDs of the map and of the files), and so is a map that gives one record two
  different genes (the records are listed); the same line twice is fine. (2) A field of the FASTA header that at
  least 90% of a file's records have, stable gene identifiers before gene symbols: `[locus_tag=...]` or
  `locus_tag=...`, `GeneID:...` (also inside `[db_xref=...]`), `gene:...` (Ensembl), `parent=FBgn...` (FlyBase),
  `locus=...` (Phytozome), then `[gene=...]` or `gene=...` (NCBI `cds_from_genomic`, GFF-derived files) and
  `gene_symbol:...`. (3) The ID itself, when IDs share a stem and differ in a transcript suffix (the table
  below). When a file's annotation follows none of these conventions, give a gene map.
- **Header styles.** What `coal pipeline` recognises, with two records of one gene each (one small test file per
  row: `tests/fixtures/isoforms/styles/`):

  | Annotation | Two records of one gene | Read from |
  |---|---|---|
  | NCBI, with locus tags | `lcl\|NC_1_cds_NP_1.1_1`, `lcl\|NC_1_cds_NP_2.1_2` with `[locus_tag=...` | header field [locus_tag= |
  | NCBI, with GeneIDs | `lcl\|NC_1_cds_NP_1.2_1`, `lcl\|NC_1_cds_XP_2.1_2` with `GeneID:...` | header field GeneID: |
  | Ensembl | `ENST00000380152.8`, `ENST00000544455.6` with `gene:...` | header field gene: |
  | FlyBase | `FBpp0070000`, `FBpp0070001` with `parent=FBgn...` | header field parent=FBgn |
  | Phytozome | `Potri.001G000100.1`, `Potri.001G000100.2` with `locus=...` | header field locus= |
  | NCBI, gene symbols only | `lcl\|NC_1_cds_NP_1.1_1`, `lcl\|NC_1_cds_XP_2.1_2` with `[gene=...` | header field [gene= |
  | GFF-derived (gffread, WormBase) | `mRNA1`, `mRNA2` with `gene=...` | header field gene= |
  | TransDecoder on Trinity | `TRINITY_DN1000_c0_g1_i1.p1`, `TRINITY_DN1000_c0_g1_i2.p1` | ID suffix _iN.pN (TransDecoder on Trinity) |
  | Trinity | `TRINITY_DN1000_c0_g1_i1`, `TRINITY_DN1000_c0_g1_i2` | ID suffix _iN (Trinity) |
  | AUGUSTUS, BRAKER | `g1.t1`, `g1.t2` | ID suffix .tN |
  | MAKER | `maker-scaf1-gene-0.5-mRNA-1`, `maker-scaf1-gene-0.5-mRNA-2` | ID suffix -mRNA-N (MAKER) |
  | MAKER (-RA, -RB) | `Gene0001-RA`, `Gene0001-RB` | ID suffix -RA (MAKER) |
  | maize (B73 v4) | `Zm00001d027230_T001`, `Zm00001d027230_T002` | ID suffix _TNNN (maize) |
  | Phytozome protein IDs | `Potri.001G000100.1.p`, `Potri.001G000100.2.p` | ID suffix .N.p (Phytozome) |
  | TAIR, Ensembl Plants, Phytozome | `AT1G01020.1`, `AT1G01020.2` | ID suffix .N |
  | RAP-DB (rice) | `Os01t0100200-01`, `Os01t0100200-02` | ID suffix -NN (RAP-DB) |

  Not recognised: UniProt headers (`GN=`; protein files anyway), IDs without a common stem (`rna-XM_...`), and
  **RefSeq transcript accessions** (`NM_000059.4`, `XM_011535000.3`): the accession holds no gene, so without a
  gene field in the header several transcripts of one gene cannot be seen; a warning says so. Use NCBI's
  `cds_from_genomic` file (its headers name the gene) or a gene map.
- **These are heuristics, and this is how they fail.**
  - *An ID suffix is a guess.* An ID that is alone with its stem is its own gene (the `.8` of
    `ENST00000380152.8` is a version, not an isoform), and a suffix rule is not applied when a stem has more than
    four records on average (`scaffold7.1` ... `scaffold7.340` numbers genes). But a numbering of paralogs
    (**`HSP70.1 / HSP70.2`**) or of a few genes per scaffold (**`scaffold7.1`, `scaffold7.2`, `scaffold7.3`**)
    looks exactly like transcripts, and `--longest-isoform` would merge different genes and keep one. Two
    guards: the summary line names the rule and the counts (`query: 6 of 26 records dropped, 20 genes (from ID
    suffix .N)`), and **when a suffix rule alone (no header field, no gene map) would drop more than half of a
    file's records, the run stops** (exit 1, nothing dropped) and asks for `--gene-map`, or for
    `--isoform-rule suffix` if the suffixes really are transcripts. Below that share a family such as
    `HSP70.1 / HSP70.2` among otherwise plain IDs is merged without a refusal: look at the `dropped` rows of
    `genes.tsv` when your IDs end in numbers and you have not checked what the numbers mean.
  - *A gene symbol is not a locus.* In a file of several strains or genomes every `[gene=dnaA]` is another gene,
    and one genome can hold two copies of `tuf`. Locus tags and GeneIDs are therefore preferred to symbols, and
    when only symbols are there, records that share a symbol but sit on different contigs (`lcl|CONTIG_cds_...`)
    or carry different locus tags or GeneIDs are **not merged**: each stays a gene of its own and a warning
    names the symbols and an example. Two copies of a gene on one contig, with nothing but the symbol in their
    headers, cannot be told from two transcripts: give a gene map.
  - *The 90% rule.* A field that fewer than 90% of a file's records have is not used (a gene map is the remedy).
- **`--longest-isoform`** keeps, for every gene, the record with the longest CDS (in nucleotides; of equally long
  ones the ID that sorts first) among the records that can be used, in both modes, before DIAMOND runs. The other
  records are listed in `genes.tsv` with status `dropped`, code `not_longest_isoform` and a note naming the gene
  and the record kept; they are not failures (the exit code is not affected). A `--pairs` or `--clusters` file
  that names a dropped record gets that pair's `error` row, or that member left out, with the reason.
  `report.json` records the choice (`isoforms`: the option, the gene map, `isoform_rule`, and per input set the
  records, genes, genes with several records, where the genes were read from, how many records the map names, how
  many were dropped and how many of those on the evidence of an ID suffix), and
  `counts.genes.dropped` the number. With `--longest-isoform` and no gene information in a file, every record
  of it is kept and a warning says so.
- The longest CDS is a convention, not the biologically "main" transcript. For a curated choice (Ensembl
  canonical, MANE Select, APPRIS) filter the FASTA file yourself and run without the option.

## Genes of another genetic code

**One run has one genetic code** (`-c N`, default 1, the standard code): it is used for the translation, the
codon alignment and Ka/Ks of every gene of the run. A CDS file of a whole genome usually also holds the
**mitochondrial** genes (13 in a vertebrate) and, in plants and algae, the **plastid** genes, and those are read
with another table: vertebrate mitochondria use code 2 (TGA is Trp, AGA / AGG are stops, ATA is Met), invertebrate
mitochondria code 5, yeast mitochondria code 3, mold and protozoan mitochondria and Mycoplasma code 4. Translated
with the standard code such a gene has internal stop codons: its row in `genes.tsv` has status `warning` and
code `internal_stop`, the codon columns with a stop are skipped in Ka/Ks, and what is left is estimated from
codons that mean something else (a pair is often `saturated`). Plastid and bacterial genes use code 11, whose
codons and stops are those of the standard code (only start codons differ), so they need nothing special.

- **The hint.** When at least two genes of an input set have two or more internal stop codons under the chosen
  code and none under another code COAL supports, the run prints one warning per set: how many of the genes with
  internal stops are clean under which code, what the choice rests on, three examples, and the recipe below.
  `report.json` has the same as `genetic_code_hints` (`clean_under`, `candidates`, `evidence`). A gene with a
  single internal TGA does not count: a selenoprotein gene would look "clean" under every code where TGA is Trp.
  The hint is a hint: it does not change a row, and a set whose every gene is flagged simply needs another `-c`.
- **How the code is chosen** (it used to be "the code that makes the most genes stop-free", which named code 14
  for yeast mitochondria and for frameshifted pseudogenes and code 4 for human mitochondria). Only genes that
  look complete count: not marked as pseudogenes (`[pseudo=true]`, a `*pseudogene` biotype), and with a length
  that is a multiple of 3 or ends in the `T` / `TA` of an incomplete stop codon. When every header names its
  molecule (NCBI `lcl|NC_..._cds_`, Ensembl `chromosome:GRCh38:MT:`), only genes on a molecule of which at least
  half the records have internal stops count (the mitochondrial chromosome, not nuclear genes here and there).
  Of the codes under
  which the genes are clean, the one wins whose extra stop codons or unused codons the genes avoid where the
  set's other genes use them: AGA / AGG (stops in code 2) for vertebrate mitochondria, CGA / CGC (unused in code
  3) for yeast mitochondria; then the one that reassigns no stop codon the genes do not need (code 14 also reads
  TAA). When several codes are equally good, all are named and none is recommended: codes 4, 5, 9, 10, 13 and 21
  all read TGA as Trp (or Cys) and differ only in codons whose meaning a DNA sequence does not show; start codons
  that only some of them allow narrow the list (Drosophila mitochondria: 4 or 5). A code that wins alone on its
  codons is named alone only when its start codons fit the genes: *C. elegans* mitochondrial genes avoid CGA /
  CGC as yeast mitochondrial genes do, but most start with ATT or TTG, which code 3 does not allow, so the hint
  names 4 and 5 (their code is 5), not 3; when nothing but such a code would be left, there is no hint. A code
  that reads TAA or TAG as an amino acid is suggested only when at least half of the set has internal stops (the
  whole set's code, e.g. ciliate nuclear genes). Sets with more than 5,000 such genes are judged on 5,000 spread evenly over the
  file (not the first ones). On real genomes: *S. cerevisiae*, *C. glabrata*, *K. lactis*, *A. gossypii* with
  their mitochondrial genes: code 3; human and mouse Ensembl: code 2 (12 genes each); *D. melanogaster*: 4 or 5;
  *E. coli*, *Y. pestis* (pseudogenes), *A. lyrata* and Ensembl mouse all-transcript files: no hint.
- **The whole set in another code.** When at least half of a set's records have internal stops, the hint is
  about the whole set, not about organelle genes in it: "most genes (376 of the 524) have internal stop codons
  under genetic code 11, ... The whole set looks like code 4 ... Run again with -c 4". The molecule rule above
  is not applied then (it left out *Mycoplasma genitalium*'s one chromosome, and the hint was silent), and when
  one molecule has 150 records or more (a genome, not an organelle) only the codes of whole genomes are
  candidates (1, 4, 6, 10, 11, 12, 15); in such a genome-size set (or one with 50 or more such genes), of codes
  that tie, those whose start codons do not fit at least a tenth of the genes drop out (bacteria start 10-15% of
  their genes with GTG / TTG); a small set keeps the quarter (*Ciona* mitochondrial genes alone: 4, 5 or 13; a
  tenth dropped its code 13 on two odd start codons). On real sets: *M. genitalium*
  with `-c 1` or `-c 11`: code 4; *Tetrahymena thermophila* and *Paramecium tetraurelia* with `-c 1` (also with
  their records shuffled): code 6; human Ensembl with `-c 2`: codes 1, 11 or 12 (nothing in the sequences tells
  them apart); a file of yeast mitochondrial genes alone with `-c 1`: code 3. `report.json` marks such a hint
  `whole_set: true`.
- **`pipeline clusters` judges each `--cds` file on its own** (each is one species): one hint per file, naming
  it (`cds sp_two.fa: most genes ...`), and `file` in `report.json`. A *Mycoplasma* among other bacteria gets
  code 4, a ciliate among animals code 6, and yeast and human mitochondrial genes get codes 3 and 2, each in its
  own file (judged over all files together, these got no hint or codes that fit neither).
- **The two-run recipe.** Take the organelle genes out of both CDS files and run them on their own:

  ```
  coal pipeline pairs --query q_nuclear.fa --subject s_nuclear.fa -o out_nuclear
  coal pipeline pairs --query q_mito.fa    --subject s_mito.fa    -o out_mito -c 2
  ```

  The two `dnds.tsv` have the same columns and can be concatenated (the second without its header line). Ensembl
  headers name the chromosome (`chromosome:GRCh38:MT:...`), NCBI's `cds_from_genomic` IDs the molecule
  (`lcl|NC_012920.1_cds_...`), so the genes are one `grep` or `seqkit grep` away.
- **`Comeron` and `Li` are for the standard code only** (`-c 1`): the two reference programs they are checked
  against have no other, and COAL refuses the combination (exit 2) instead of guessing. The other seven methods
  take every table (`coal --help` lists them).
- There is no per-gene code option (`--code-map`); [development.md](development.md#a-genetic-code-per-gene-design-note)
  has the design and why it was not built in this release.

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
  query file. `--diamond-args` is one string of further `diamond blastp` options (`--diamond-args '--block-size 4
  --tmpdir /scratch'`); it is split like a command line and never goes through a shell. The options COAL sets
  itself are refused in it (exit 2, before anything runs): `--query` / `-q`, `--db` / `-d`, `--in`, `--out` /
  `-o`, `--outfmt` / `-f`, `--threads` / `-p` (use `-j`), `--max-target-seqs` / `-k` (1: the best hit is the
  rule), `--evalue` / `-e` (use `--evalue`), `--quiet`, a sensitivity such as `--very-sensitive` (use
  `--sensitivity`), and with `--aligner diamond` also `--masking`. Before, `--diamond-args '--out x.tsv'` sent
  the hits elsewhere and the run ended with an unrelated error.
- **`--pairs TSV`** takes the pairs from a file (query ID, subject ID; a header line `query subject` is skipped)
  and DIAMOND is not run. `out/orthologs.tsv` of an earlier run can be given, also with the same `-o out`: the
  file is then read and left exactly as it is (not written again).
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

- **Clusters.** Without `--clusters`, `diamond deepclust` runs on the translated proteins, sorted by ID so that
  the result does not depend on the order of the input
  ([what the defaults do](#cluster-mode-input-order-and-what-a-default-cluster-is)), with DIAMOND's defaults
  (member coverage 80%); `--approx-id`, `--member-cover`, `--mutual-cover`, `--cluster-evalue` and
  `--deepclust-args '...'` are passed on, and the result is kept as `out/clusters.tsv`. `--realign` also runs
  `diamond realign` (members against their representatives, the columns orthologr's `deepclust_realign()` asks
  for) into `out/realign.tsv`: `member representative approx_pident evalue bitscore qcovhsp scovhsp`, with a
  header line. `--deepclust-args '...'` and `--realign-args '...'` are each one string of further options for
  that DIAMOND command (`--deepclust-args '--memory-limit 12G'`), split like a command line, no shell. The
  options COAL sets itself are refused (exit 2): in `--deepclust-args` `--db` / `-d`, `--out` / `-o`,
  `--threads` / `-p`, `--quiet`, and `--approx-id`, `--member-cover`, `--mutual-cover`, `--evalue` when COAL's
  own option for it is given as well; in `--realign-args` `--db`, `--clusters`, `--out`, `--outfmt` / `-f` (the
  columns of `realign.tsv` are fixed), `--threads`, `--quiet`.
- **Several `--cds` files are one set of genes**, so an ID must name one record across all of them. IDs that
  occur in two files (BRAKER / AUGUSTUS `g1.t1`, `gene1` in every species) are refused before anything runs, with
  both files and examples named: only the first file's record of such an ID could be used, and with every ID
  repeated a whole species would be left out. `--prefix-ids` puts each file's name (without folder and FASTA /
  gzip suffix) and `|` before its IDs instead: `S_cer.fa`'s `g1.t1` is `S_cer|g1.t1` in every output, the report's
  `inputs` give each file's `id_prefix`, a `--clusters` file names the prefixed IDs, and a `--gene-map` may name
  either (a gene ID is then a gene of its species: `S_cer|GENE1` and `C_gla|GENE1` are two genes). Two files with
  the same name in different folders cannot be told apart this way and are refused with `--prefix-ids`.
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
| `-m METHOD` | Ka/Ks method, repeatable: YN, NG, LWL, MLWL, LPB, MLPB, MYN, Comeron (alias C95), Li ([kaks.md](kaks.md)); names are read in any case; Comeron and Li need `-c 1` | YN |
| `-c N` | NCBI genetic code, used for translation, the codon alignment and Ka/Ks | 1 |
| `-j N` | threads of DIAMOND and of the engine, and MSA processes at a time. The default is one per CPU the process may use (its affinity: `taskset`, a SLURM or cgroup allocation), not one per CPU of the machine; a larger `-j` is cut to that number for every program, with a warning, and `report.json` `threads` has what was asked and used | one per usable CPU |
| `--kaks-format` | also write `dnds.kaks` in KaKs_Calculator's 22-column format | off |
| `--keep-alignments` | also write the codon alignments | off |
| `--keep-protein-alignments` | also write the protein alignments the codon alignments are made from | off |
| `--resume` | reuse the translation, search and clustering results in `DIR/work` when the inputs (path, size, modification time), the options and the DIAMOND version are the same | off |
| `--install-tools` | install a missing DIAMOND / Clustal Omega / MAFFT into COAL's cache | off |
| `--diamond PATH` | the DIAMOND to run | found automatically |
| `--report FILE` | the JSON report | `DIR/report.json` |
| `--hash-inputs auto\|always\|never` | SHA-256 of the input files in the report: files up to 2 GiB, every file, or none | `auto` |
| `--msa-timeout S` | give up on one Clustal Omega / MAFFT run after S seconds | none |
| `--longest-isoform`, `--gene-map TSV` | keep the longest CDS per gene; the genes from a map where it names the record, else from the headers and IDs ([One record is one gene](#one-record-is-one-gene)) | off |
| `--isoform-rule auto\|suffix\|map` | `auto`: header fields, else ID suffixes, and `--longest-isoform` stops when a suffix rule alone would drop more than half of a file's records; `suffix`: the ID suffixes are transcripts, whatever their share; `map`: only `--gene-map`, nothing is read from headers or IDs | `auto` |
| `--min-codons N`, `--max-ds X`, `--min-identity P` | flag estimates from fewer than N compared codons / with dS above X / from a pair whose protein identity is below P percent: status `warning` and a note; no row dropped, no number changed ([Filtering before interpretation](#filtering-before-interpretation)) | off |
| `--yn-flavour`, `--yn-weighting`, `--digits`, `--max-cells` | as in `coal dnds` and `coal-cpp-align` | |
| `--comeron-flavour coal\|gestimator`, `--li-flavour coal\|seqinr` | how `Comeron` and `Li` read a pair: COAL's policy, or as their reference programs (stop codons compared; [kaks.md](kaks.md), "Comeron and Li") | `coal` |
| `--terminal-stops drop\|compare` | only with the reference flavours: `compare` takes the two terminal stop codons of a pair as the last compared codon column when both genes end in one, as orthologr 0.4.3's pipeline does; `drop` ends the pair before them | `drop` |

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
| `protein_alignments.fasta` | pairs | `--keep-protein-alignments`: the pairwise protein alignments (COAL's aligner, the MSA tool's, or DIAMOND's completed to full length, as `--aligner` says), records 1+2, 3+4, ...: the same records in the same order as `codon_alignments.fasta`, and each codon row is three times its protein row, column for column |
| `protein_alignments/<cluster>.fasta` | clusters | `--keep-protein-alignments`: the protein MSA of each cluster as it was used (the rows of Clustal Omega / MAFFT for the members that could be used; the engine's pairwise alignment for a cluster of two), under the same file names as `codon_alignments/`; not for clusters aligned pairwise |
| `work/` | both | proteins, DIAMOND tables, logs of every external command (`work/logs/`), the `--resume` state |

**A file is complete or absent.** Every table and alignment file is written as `NAME.tmp` and renamed when it is
complete: the engine's files (`dnds.tsv`, `dnds.kaks`, `codon_alignments.fasta`, `protein_alignments.fasta`) when
the Ka/Ks stage has ended with success, each per-cluster alignment file when it is written, `clusters.tsv` when
DIAMOND has ended, and the tables COAL writes itself (`genes.tsv`, `orthologs.tsv`, `unpaired.tsv`,
`cluster_summary.tsv`, `report.json`) with their last line; `failed.tsv` is `failed.tsv.tmp` while the run goes on
and gets its name when the run ends, also when it ends with an error. So after a stage that failed (a full disk,
an engine error) or was interrupted there is **no cut-off table or alignment under a final name**: the file is not
there, the report lists only the files that are, and the `.tmp` files are removed (what a `kill -9` leaves is
removed by the next run into the folder, which knows the names from `work/outputs.json`). The temporary files are
created anew (`O_EXCL`, and `O_NOFOLLOW` where the system has it), and the rename replaces whatever has the final
name: **a symbolic link in the output folder is replaced, never written through** (an output name that is a link
to another file leaves that file untouched). `report.json` shows the engine's command line with the `.tmp` names
and `--new-files`, as it ran.

**One run per folder, and one run's files only.** A run takes the lock `<out_dir>.coal-lock` (next to the real
folder: a symbolic link to the folder and the folder are one output; COAL's convention for `-o` files) and a
second run into the same folder is refused with exit 1 while the first one is alive. The lock is one line,
`<process id> <host name> <process start> <token>`. The lock of a process that is gone (or whose id now belongs
to a process that started at another time; Linux) is replaced, by exactly one of the runs that find it. A lock
written on another machine cannot be checked and is never replaced: the message names the machine and the file
to remove if that run is gone. The host name is the same text in both engines: the computer's DNS host name
(Windows: `GetComputerNameExW(ComputerNamePhysicalDnsHostname)`; elsewhere `gethostname`), printable ASCII as it
is, any other ASCII character as `_`, every character outside ASCII as `?`, at most 255 characters, `-` when there
is none; two names are the same machine when they are equal (on Windows without regard to case). While a stale
lock is replaced, the file `<lock>.takeover` (the same four fields) exists for an instant; the run that created it
removes it. Should one stay behind (a run killed in that instant), the next run clears it when its process is gone
or when it is older than 15 s. When a run starts, the files
of the table above that an earlier COAL run wrote into the folder are removed (a `dnds.kaks` or
`codon_alignments.fasta` or `protein_alignments.fasta` of a run with other options, the `.fasta` files of a
`codon_alignments/` or `protein_alignments/` that a cluster run with `--keep-alignments` /
`--keep-protein-alignments` wrote (`work/outputs.json` records both folders), a half-written `dnds.tsv.tmp`), so the folder never mixes two runs. COAL
removes only what a run is known to have written: each run lists the files it may write in `work/outputs.json`
before it writes the first one (for a folder of a version without that list, the outputs named in its
`report.json`). A file with one of those names that no run wrote (your own `clusters.tsv`, your own files in
`codon_alignments/` or `protein_alignments/`) stays, files with other names are not touched, and `work/` is kept for `--resume`.

**Inputs are never removed or written.** A file given as `--query`, `--subject`, `--pairs`, `--cds`, `--clusters`
or `--diamond` is recognised under any name (a link, another spelling of the path). If the run would write a file
of that name (`--pairs out/dnds.tsv -o out`, `--report` onto an input, an input inside `work/`,
`codon_alignments/` or `protein_alignments/`) it ends with exit 1 before anything in the folder changes; the one exception is
`--pairs out/orthologs.tsv`, which is read and kept as it is. `--keep-alignments` in cluster mode is refused in
the same way when `codon_alignments/` exists, is not empty and was not written by a COAL run, and
`--keep-protein-alignments` when `protein_alignments/` is such a folder.

The alignment files end before a terminal stop codon (also with `--terminal-stops compare`, which changes what
Ka/Ks compares and not the files), and hold the members that could be used: a gene that failed is in
`failed.tsv` and in neither file.

### `dnds.tsv`

| Column | Meaning |
|---|---|
| `cluster_id` | cluster mode only (first column): the cluster's representative |
| `seq1`, `seq2` | the two genes. Pair mode: query, subject |
| `method` | NG, LWL, MLWL, LPB, MLPB, YN, MYN, Comeron or Li; with several methods a pair's rows come in this order |
| `dN`, `dS`, `dNdS` | Ka, Ks, Ka/Ks. `NA` when not available; `dNdS` is `NA` when dS < 1e-6 |
| `SE_dN`, `SE_dS` | standard errors (YN, MYN, NG, Li); `NA` for the other methods |
| `N_sites`, `S_sites` | nonsynonymous and synonymous sites |
| `Nd`, `Sd` | nonsynonymous and synonymous differences |
| `codons_compared` | codon columns used (no gap, ambiguous base or stop in either sequence) |
| `prot_aln_length` | columns of the protein alignment with a residue in at least one of the two sequences |
| `prot_identity` | percent of those columns with the same residue in both (two decimals) |
| `status` | `ok`; `warning` (estimated, see note); `not_converged` (YN / MYN: the iteration did not converge, so there is no estimate: dN, dS, their SEs, Nd and Sd are `NA`; [kaks.md](kaks.md#platform-reproducibility)); `saturated` (dN or dS could not be corrected: `NA`); `no_data` (nothing to compare); `error` (the pair could not be computed: every value `NA`) |
| `note` | why the status is not `ok`: internal stop codons, a CDS length that is not a multiple of 3 (with `seq1:` / `seq2:`), the saturated estimate, the iteration that did not converge (with the values of its last round, which are not estimates and can differ between platforms), or the error |

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
| `n_not_converged`, `n_saturated`, `n_no_data`, `n_failed` | pairs with status `not_converged`, `saturated`, `no_data`, `error`. None of them is in `n_ok` or in the statistics |
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
| not the longest record of its gene, with `--longest-isoform` | `dropped`, `not_longest_isoform` | the record is left out before the search; the note names the gene and the record kept ([One record is one gene](#one-record-is-one-gene)) |
| empty, shorter than one codon, only a stop codon | `failed`, `cds_empty` / `cds_too_short` / `cds_only_stop` | the gene is left out |
| a character that is no nucleotide code | `failed`, `cds_invalid_character` | the gene is left out |
| an ID that already occurred | `failed`, `duplicate_id` | the first record is used, this one is left out |

### `report.json`

Schema `coal.pipeline/1`: [pipeline.schema.json](pipeline.schema.json). It is written also when the run fails.

- `mode`, `command`, `out_dir`, `parameters` (every option as used), `tools` (path and version of coal-cpp,
  DIAMOND, Clustal Omega, MAFFT).
- **Provenance, for a methods section.** `started` and `finished` (UTC, ISO 8601: `2026-11-13T09:15:46Z`);
  `platform` (`system`, `release`, `machine`, `python`, `cpus`); `coal_version`; `threads` (`requested`: `-j` as
  given, 0 for the default; `usable_cpus`: the CPUs this process may run on; `engine`: the threads coal-cpp was
  given; `tools`: what DIAMOND and the aligners were given; neither is ever more than `usable_cpus`); `inputs`: every
  file the run read, with `option`, absolute `path`, `bytes`, `modified` and `sha256` (of the file as given: a
  `.gz` input is hashed compressed and has `"compression": "gzip"`, and also `sha256_uncompressed`, the hash of
  its content: the same FASTA compressed twice has two `sha256` and one `sha256_uncompressed`). Hashing reads each input once more; with
  the default `--hash-inputs auto` a file above 2 GiB is not hashed (`"sha256": null`, and `sha256_skipped` says
  why); `--hash-inputs always` hashes every file, `never` none. `tools.<name>.commands`: the exact command
  lines, as argv arrays, of every DIAMOND and coal-cpp run (the same lines as in `work/logs/`; kept in the
  `--resume` state, so a reused stage still shows how its files were made), and for Clustal Omega / MAFFT,
  which run once per cluster, `command_template` (one run's command with `IN.fasta` / `OUT.fasta`). `cite`: the
  citation keys for what this run used; `coal cite DIR/report.json` prints the references
  ([kaks.md](kaks.md#what-to-cite)). A `clusters` run that used DeepClust lists DIAMOND
  (`buchfink2021`) and the DeepClust paper (`buchfink2026`: Buchfink et al. 2026, Nature Methods 23:724-727,
  <https://doi.org/10.1038/s41592-026-03030-z>; reports written before it appeared name its preprint,
  `buchfink2023`, which `coal cite` still resolves). `warnings`: the warning lines of the run, when there were any. None of
  this changes a table: `dnds.tsv` and the others are byte-identical whatever `--hash-inputs` is.
- `stages`: name, wall `seconds`, `cached` (reused by `--resume`). The stage `align_codon_kaks` also has
  `align_codon_seconds` and `kaks_seconds` (inside the engine) and, when an MSA program ran, `msa_seconds`.
- `counts`: `genes` (total, ok, warning, failed), `pairs` (total and one count per status), `flagged`
  (with `--min-codons` / `--max-ds`: `min_codons`, `max_ds`, the pairs each flag marked), `orthologs` and
  `unpaired` (pair mode), `clusters` (total, singletons, all_pairs, above_max_size, skipped, msa_fallback,
  members_dropped).
- `failures`: the first 1,000 entries of `failed.tsv` (`failures_total`, `failures_truncated`), each with `level`
  (`gene`, `pair`, `cluster`), `code` and `message`.
- `exit_code`, `status` (`completed`, `completed_with_failures`, `failed`, `interrupted`), and `error` when the
  run could not be done or was interrupted.

### Failures and exit codes

One bad gene, pair or cluster never stops a run.

- A gene that cannot be used is left out and listed (`genes.tsv`, `failed.tsv`). In pair mode a pair that names
  it (`--pairs`) gets a row with status `error`; in cluster mode it is dropped from its cluster.
- A pair that cannot be computed (alignment too large for `--max-cells`, the codon engine refused a sequence)
  gets a row with status `error` and the reason.
- `saturated`, `no_data`, `not_converged` and `warning` are results, not failures (the exit code stays 0). Only
  `ok` and `warning` have estimates.

| Exit code | Meaning |
|---|---|
| 0 | finished, nothing failed |
| 2 | finished, but `failed.tsv` has entries (as `coal --keep-going` and `coal-cpp-align`); also a command-line error |
| 1 | the run could not be done: an input cannot be read, a tool is missing or failed, the folder cannot be written (a full disk included) |
| 130 | interrupted by Ctrl-C (SIGINT); 143 for SIGTERM, 129 for SIGHUP (128 + the signal's number) |

**Interrupts.** When Ctrl-C or a termination signal reaches the `coal pipeline` process (also when only that
process gets it and not its children, as from `kill`, a scheduler or a notebook), the programs it started are
ended: the `coal-cpp codon-dnds` engine is killed and not left to finish as an orphan, and so are DIAMOND and the
running Clustal Omega / MAFFT processes. `dnds.tsv.tmp` is removed, `report.json` is written with
`"status": "interrupted"` and the lock is released; `dnds.tsv` is absent, so an interrupted folder cannot be taken
for a finished one. `--resume` then reuses the finished stages. A process killed with SIGKILL cannot do any of
this: the engine then ends when its input pipe closes, and the next run replaces the stale lock.

**A full disk, or any other error of the operating system, ends the run with one line.** The message names the
file and the system's reason (`cannot write out/work/state.json.tmp: No space left on device. The disk (or quota)
of the output folder is full: ...`), the exit code is 1 and the lock is released. `report.json` then says
`"status": "failed"` with that error; it is written through a temporary file, so it is complete or absent, never
cut short. When the whole report does not fit, a short one is tried (`schema`, `exit_code`, `status`, `error`,
`"report_truncated": true`). A run that finished its tables but could not write `failed.tsv` or its report does
not end with exit 0. A tool that was ended by a signal is reported as that: `diamond blastp failed (killed by
signal 9 (SIGKILL; out of memory?))`.

**Warnings of the run** (the lines `coal pipeline: warning: ...` on stderr) are also in `report.json` as
`warnings`: the same file given as `--query` and `--subject` (every gene's best hit is then itself: dN = dS = 0;
with `--pairs`, only the lines that pair a gene with itself are counted), a cached MAFFT or Clustal Omega that is
not the pinned version, a `-j` above the number of CPUs, an old DIAMOND, a `--pairs` / `--clusters` file whose IDs are not those of the
sequences, input that looks like several transcripts per gene, and the counts of the `--min-codons` / `--max-ds`
flags.

**IDs that match nothing.** When 90% or more of the IDs of a `--pairs` or `--clusters` file are not records of
the sequence files, one line says so with an example from each side (`12 of the 12 query IDs of --pairs p.tsv are
not in the query sequence file(s): e.g. 'AT1G01010' in p.tsv, while the sequences are named like 'AT1G01010.1'`)
and, when that is what it looks like, that the IDs differ by a version or transcript suffix. Every pair still gets
its row with status `error` (pair mode) and its line in `failed.tsv`, and the exit code is 2 as before.

**Inputs that cannot work are reported before the long steps.** When a CDS file has no usable gene at all, the
run stops after translation, before DIAMOND or an aligner is started, and says why: "does not hold coding DNA ...
Is this a protein file?" when most sequences have letters that are no nucleotide codes, otherwise the count of each
reason (`cds_empty`, `cds_too_short`, ...). `genes.tsv` and `failed.tsv` list every gene. A sequence ID with a
non-ASCII character is that gene's failure (`id_not_ascii`), not the end of the run; the next paragraph says how
the ID is shown. Non-ASCII text elsewhere in a sequence file (in the sequence lines) still stops the run: such a
file is not a FASTA file of nucleotides.

**Text encoding: the same bytes on every system.** Sequence IDs are bytes that the pipeline never re-encodes.
Every file it writes (`dnds.tsv`, `genes.tsv`, `failed.tsv`, `orthologs.tsv`, `cluster_summary.tsv`, `report.json`,
the alignments) is UTF-8 with `\n` line ends, on Windows too: neither the system's default encoding (cp1252 on
Windows, ASCII in the C locale) nor its line ends are ever used, for reading or for writing, and `coal-cpp` writes
its tables and its standard output in binary mode. An ID, in a sequence file or in a `--pairs` / `--clusters`
file (which may have `\r\n` line ends), that is valid UTF-8 is written as it is (`gène_αβ` stays `gène_αβ`). Bytes
that are not valid UTF-8 (an ID saved as Latin-1, say) are written as `\xhh`, two lower-case hexadecimal digits
per byte: `gène` in Latin-1 is shown as `g\xe8ne`. Nothing is dropped or replaced by `?`, two different IDs never
become the same text, and every output file is valid UTF-8 whatever the input was. The same rule is used by
`coal-cpp translate` (the genes table) and by the Python layer (Python's `backslashreplace` error handler), so a
gene has one name in all tables. Messages on the console follow the console: what it cannot show is escaped.
Every input (sequence files, `--pairs`, `--clusters`, `--gene-map`, `coal --batch` job files) may start with a
UTF-8 byte order mark (Windows Notepad writes one), which is skipped, and may have `\n`, `\r\n` or classic Mac
`\r` line ends, in every reader alike (`coal-cpp` and the Python layer). Read the tables as UTF-8 (`encoding="utf-8"` in Python on Windows; `fileEncoding = "UTF-8"` in R).

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

## Filtering before interpretation

`dnds.tsv` holds every pair that could be computed, and status `ok` means that the method produced numbers, not
that they are informative. [kaks.md](kaks.md#filtering-before-interpretation) explains what to check and why; in
short, for a pipeline table:

- **Short alignments**: `codons_compared` below about 50. A handful of substitutions decides the ratio.
- **Low coverage or identity**: `codons_compared` well below `prot_aln_length` (under half), or `prot_identity`
  below 30 to 40%. Reciprocal best hits are not always orthologs over their whole length, and cluster members
  need not be alignable at all.
- **Saturation**: dS above about 2 (1 to be strict). `saturated` marks only the pairs without any value; a pair
  with dS 4 and status `ok` is no better. When more than 30% of a run's pairs (of at least 20) are `saturated` or
  `not_converged`, the run says so in one warning that points here (*S. cerevisiae* / *C. glabrata*: 86%). Default DeepClust clusters hold many such pairs
  ([above](#cluster-mode-input-order-and-what-a-default-cluster-is)).
- **Tiny dS**: dN/dS is a ratio, and with dS near 0 (very close sequences) it is unstable and often large for no
  biological reason; `dNdS` is `NA` only when dS is below 1e-6.

These thresholds are conventions, with the usual caveat: report the ones you use, and check that the result does
not hinge on them. `--min-codons N`, `--max-ds X` and `--min-identity P` flag such rows (status `warning`, notes
`fewer than N codons compared`, `dS above X` and `protein identity below P%`) **without dropping a row or
changing a number**; the counts are in the summary on stderr and in `report.json` (`counts.flagged`). They are
off by default. `--min-identity` compares the column `prot_identity` (identical residues as a percentage of the
protein alignment's columns, gap columns included) with P.

```
coal pipeline pairs --query q_cds.fa --subject s_cds.fa -o out --min-codons 50 --max-ds 2 --min-identity 30
```

**A cluster can hold pairs that are not related.** Clustering puts a sequence with a representative it matches;
two members of one cluster need not match each other, and with DeepClust's defaults many do not. Every pair of a
cluster is still aligned and estimated, and an alignment of two unrelated proteins gives numbers like any other:
in a 3,000-member cluster of the second review, 306 of 1,000 sampled pairs had status `ok` at a protein identity
below 15% (dN around 1, dS around 1.7). Status `ok` says that the method ran. Use `--min-identity` (30 to 40 is
a common floor), or filter on `prot_identity` yourself, before reading such rows.

## How much the aligner matters

The protein alignment decides which codons are compared, so `--aligner` moves dN/dS. On 150 reciprocal best hits
of *Arabidopsis thaliana* and *A. lyrata* (the first 150 genes of orthologr's example files, method YN, DIAMOND
2.2.8; second review, 2026-11-13), against the default `nw`:

| `--aligner` | pairs with exactly the same dN, dS and dN/dS as `nw` | dN/dS differs by more than 0.01 | by more than 0.05 | largest difference in dN/dS |
|---|---|---|---|---|
| `mafft` | 120 of 150 | 20 | 5 | 0.18 |
| `diamond` (DIAMOND's alignment, completed to full length) | 90 of 150 | 33 | 8 | 0.43 |
| `clustalo` | 75 of 150 | 53 | 20 | 0.30 |

The pairs themselves do not depend on the aligner, and `--ortho bh` gave the same 150 pairs and values as `rbh`
on this set. Half to four fifths of the pairs are untouched (alignments without a gap have one answer); the rest
are the pairs with gaps, where the programs place them differently. None of the four is "the" alignment: **say
in your methods which aligner was used** (`report.json` has it: `parameters.aligner`, the tool's version and
command line, and `cite` lists its paper), and when a conclusion rests on single pairs, look at them with a second
aligner (`--keep-protein-alignments` writes the alignments that were used).

## What the search and clustering defaults do

The numbers in this section come from an independent review of the pipeline (November 2026, DIAMOND 2.2.8); they
describe the examples named and are not guarantees for other data.

### Pair mode: the best hit, and DIAMOND's sensitivity

- **The best hit is DIAMOND's first target.** Each search runs with `--max-target-seqs 1`, and COAL takes the one
  target DIAMOND reports for a query. DIAMOND orders a query's targets by e-value and then by score, so this is
  the hit with the lowest e-value, not always the one with the highest bit score: a shorter, cleaner hit can have
  a slightly lower bit score and a lower e-value than a longer one. On yeast (*S. cerevisiae* against
  *C. glabrata*, 4,387 reciprocal best hits) `--max-target-seqs 1` gave exactly the first row of a 25-target
  search for every query; 15 queries had a later target with a higher bit score, and a reciprocal-best-hit set
  built with a pure bit-score rule (as orthologr's COAL fork chooses) differed from COAL's in 20 and 17 pairs
  (pairs only in one set and only in the other), about 0.4% of the pairs. On orthologr's 1,000-gene
  *Arabidopsis* example the two rules gave the same pairs. To use another rule, make the pairs yourself and give
  them with `--pairs`.
- **DIAMOND's default sensitivity** is used unless `--sensitivity` says otherwise. It is slower than `--fast`
  and finds more of the short proteins. `--sensitivity fast` is a reasonable choice for closely related
  species, where almost every ortholog is a strong hit.
- **A worked comparison with orthologr's COAL fork** (`ortho_thal_cds_1000.fasta` against
  `ortho_lyra_cds_1000.fasta`): COAL finds 999 reciprocal best hits; the fork's DIAMOND route finds 995, or 983
  with `delete_corrupt_cds = TRUE`. The whole difference is explained:

  | Pairs | In COAL, not in the fork | Why |
  |---|---|---|
  | 2 | AT1G02074.1, AT1G10747.1 | the fork searches with `--fast`, which gives these 48 to 83 residue proteins no hit or another one |
  | 2 | AT1G04150.1, AT1G08230.2 | the fork drops a CDS that has an `N`; COAL keeps it (the codon with the `N` is skipped in Ka/Ks) |
  | 12 | (only with `delete_corrupt_cds = TRUE`) | the fork deletes a CDS whose length is not a multiple of 3; COAL ignores the last 1 or 2 bases and gives the pair status `warning` |

  The best-hit rule changes nothing in this example. COAL's choices keep more genes; the rows they add are marked
  (`warning`, and the note says why), so they can be filtered out
  ([Filtering before interpretation](#filtering-before-interpretation)).

### Cluster mode: input order, and what a default cluster is

- **The clusters do not depend on the order of the input.** DIAMOND DeepClust is order-dependent: which sequence
  becomes a representative, and with that which sequences end up together, changes with the order of the
  sequences it is given. COAL therefore gives DeepClust the proteins **sorted by ID** (the bytes of the ID;
  `work/proteins.sorted.faa`), so the same set of sequences gives the same `clusters.tsv`, `dnds.tsv` and
  `cluster_summary.tsv` whatever the order of the `--cds` files and of the records in them (tested by shuffling:
  `tests/test_pipeline_order.py`). Up to this change the proteins were clustered in input order: on
  orthologr's two 1,000-gene files, giving the files in the other order changed 838 clusters into 836 and moved
  about 45 clusters' members. Renaming sequences still changes the order, and so can the clusters; a clustering
  you want to keep is best saved (`clusters.tsv`) and given back with `--clusters`. `genes.tsv` keeps the input
  order.
- **Two things still differ from orthologr's DeepClust route** on the same input: it drops the two genes with
  an `N` and leaves the terminal `*` on the proteins it clusters (three clusters differ for that reason in the
  example above).
- **DeepClust's defaults make loose clusters.** Without `--approx-id` DeepClust has no identity threshold beyond
  what its e-value and the 80% member coverage imply, so a default cluster holds distant paralogs as well as
  close ones. On the 2,000 *Arabidopsis* CDS of the example the default clusters gave 2,427 pairs, of which 574
  (24%) were `saturated`, 63 `not_converged` and 8 `no_data`: 27% of the pairs have no estimate, because the
  sequences are too far apart for dS to be estimated. That is the expected result of comparing distant family
  members, not a fault of the run. For dN/dS within families of recent duplicates or orthologs, tighten the
  clustering: `--approx-id 50` to `--approx-id 70` (percent identity to the representative) and
  `--mutual-cover 70` to `--mutual-cover 80` (both sequences covered, instead of only the member) are values
  to start from; look at `n_saturated` in `cluster_summary.tsv` and at the dS distribution, and tighten until
  most pairs have dS below about 2 ([Filtering before interpretation](#filtering-before-interpretation)).

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
| YN / MYN: the iteration does not converge (short or nearly saturated pairs; [kaks.md](kaks.md), difference 12) | `not_converged`; dN and dS are `NA`, the note has the last round. Counted in `n_not_converged`, and in the report under `counts.pairs.not_converged` |
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
- **The lock is a file, not an operating-system lock.** A lock left by a run on another machine that was killed
  stays until you remove it (the message says which file); COAL cannot tell it from a run that is still going.
  The process start time that tells a reused process id from the holder is recorded on Linux only. `coal-cpp`
  follows the same rules for its own `-o` lock (the same four fields, the host check, the start time, the
  one-at-a-time replacement of a stale lock), so Python and C++ runs honour each other's locks and runs of
  either engine that start together onto a stale lock have one holder (`tests/test_output_lock.py`: 200 trials
  of six simultaneous `coal-cpp -o` starts, and mixed Python and C++ starts). On Windows the C++ code checks the
  holder through a process handle and creates the lock through a hard link, like the Python code; it is compiled
  by CI, but these lock tests have only been run on Linux.
- **A file the run writes replaces a file of that name that was already in the folder**, as with any `-o`: a
  `dnds.tsv` or `genes.tsv` of your own in the output folder is overwritten (unless it is one of this run's
  inputs, which is refused). Only the removal at the start is limited to what an earlier run wrote.
- **SIGKILL (and a power cut) leave an interrupted folder without `"status": "interrupted"`**: `report.json` is
  then absent (it is removed when a run starts and written when it ends), and so is `dnds.tsv`.
- **One MSA per cluster is held in memory**, and Clustal Omega's time grows faster than linearly with the cluster
  size. Lower `--max-cluster-size`, or use `--aligner mafft`, for inputs with many clusters of hundreds of members.
- **Large clusters are aligned pairwise**, so their values come from pairwise alignments and those of smaller
  clusters from an MSA. The summary says which (`aligner`).
- **Pair alignment memory** is cells / 2 bytes per thread up to 256 MiB; a larger pair is aligned with
  checkpoints in about that much (the same alignment, its forward pass run twice): a titin-sized pair needs about
  130 MB instead of 630 MB. Above `--max-cells` (default 4·10⁹ cells) that pair gets status `error`.
- **`--aligner diamond` uses one HSP** and DIAMOND's heuristic alignment ([alignment.md](alignment.md)); a hit that
  does not fit the translated sequences is that pair's error.
- **Best hits are DIAMOND's first target per query** (`--max-target-seqs 1`). Ties are broken by DIAMOND
  ([what the defaults do](#pair-mode-the-best-hit-and-diamonds-sensitivity)).
- **DIAMOND 2.1.0 or later.** 2.1.0 is the first release with `diamond deepclust`; every test and benchmark here
  ran with 2.2.8, the version `coal pipeline install diamond` installs. The version is checked when a run starts:
  with an older DIAMOND (Ubuntu 22.04's package is 2.0.14) cluster mode without `--clusters` stops and says how to
  get a newer one, and pair mode goes on with a warning, because no version older than 2.1.0 has been tested.
- **`-j` above the number of CPUs** is passed to the engine as it is, but DIAMOND and the aligner processes get
  at most one thread or process per CPU of the machine (`-j 5000` made DIAMOND reserve memory for 5,000 threads
  and be killed); a warning says so.
- **Compressed input.** A gzip-compressed `--query`, `--subject`, `--cds`, `--pairs` or `--clusters` file
  (recognised by its first two bytes, whatever its name) is decompressed into `DIR/work/inputs/` and read from
  there, so the folder needs room for the plain copies; they are kept for `--resume` and reused while the
  compressed file is unchanged. Messages of the translation step name the copy. Other compression formats
  (bzip2, xz, zstd) are not read: decompress them first.
- **Protein input is not accepted**: the pipeline starts from CDS, so that protein and CDS always match. A protein
  file is recognised only when it leaves no usable gene (see "Inputs that cannot work"); a file that mixes CDS and
  protein records runs, with the protein records listed as failed genes.
- **Python holds every ID** (and a byte offset per record), so memory grows with the number of sequences, not
  with their length. The check for several transcripts per gene reads every FASTA header once more and holds the
  header lines while it looks at them.
- **Tools in COAL's cache are used only with their checksum record.** A DIAMOND, Clustal Omega or MAFFT that an
  earlier COAL (0.5.0 or before) installed into the cache has no `<binary>.coal-manifest.json` and is not used:
  the run says so once and asks for `coal pipeline install <tool> --force` (or `--install-tools`). The record
  lies next to the binary, so it guards against old, incomplete or altered installs, not against someone who can
  write to the cache folder ([alignment.md](alignment.md), "What is in the cache is checked").
- **Installers verify what they download** ([alignment.md](alignment.md#clusters-clustal-omega-and-mafft-coalmsa),
  "Verification"). The DIAMOND 2.2.8 release archives (Linux x86_64, macOS, Windows x64) and micromamba are pinned
  by SHA-256 in COAL's source; DIAMOND publishes no checksums, so those hashes are the ones of the files as
  downloaded on 2026-11-12. Clustal Omega and MAFFT come from Bioconda (conda / mamba / micromamba; Linux and
  macOS); `--install-tools` during a run never installs an unverified file. DIAMOND's macOS archive is one
  universal binary (an x86_64 slice for macOS 10.15 and later, an arm64 slice for macOS 11 and later that carries a code signature; only
  system libraries), so the same pinned file serves Intel and Apple-silicon Macs without Rosetta (read from the
  2.2.8 archive, 2026-11-13; COAL's CI installs and requires it on the arm64 macOS runner).
- **Windows: no automatic Clustal Omega or MAFFT install.** The pre-compiled Clustal Omega binaries of the
  Clustal site are gone (the address serves an unrelated site since 2026-11) and Bioconda has no Windows packages,
  so `--install-tools` and `coal pipeline install clustalo` cannot install one there and say what works: run COAL
  under WSL; or put a `clustalo.exe` you trust on `PATH` (or set `COAL_CLUSTALO`); or `--aligner mafft` with a
  MAFFT you installed yourself (`COAL_MAFFT`; not tested on Windows). Pair mode and DIAMOND are not affected
  (DIAMOND's Windows release is pinned). [alignment.md](alignment.md#clusters-clustal-omega-and-mafft-coalmsa).
- Thread scaling beyond 2 cores and whole-proteome wall times are to be measured on the 16-thread desktop.
