# Changelog

## Unreleased

**In short.** This release makes COAL more than a PAL2NAL replacement, and changes its licence.

- **New scope.** Next to codon alignments, COAL now estimates Ka/Ks (dN/dS) itself (`coal kaks` in
  KaKs_Calculator's format, `coal dnds` as a tidy table; nine methods), runs the whole way from coding sequences to
  a dN/dS table in one step (`coal pipeline pairs` / `clusters`, with DIAMOND and Clustal Omega / MAFFT), aligns
  proteins (`coal-cpp-align`) and has an R package (`coalr`).
- **Licence: GPL-3.0-only from this release** (0.5.0 and earlier stay MIT), because the Ka/Ks engine is derived
  from GPL-3 code. Running COAL puts no condition on your data or scripts; see the first entry under "Changed".
- **For users of 0.5.0, the existing commands behave as before** except for what "Changed" lists: the `-o` lock
  (two runs onto one output file are refused), tools in COAL's cache need their checksum record (re-install them
  once: `coal pipeline install <tool> --force`), `pal2nal` without arguments is a usage error, and `coal --help`
  ends with the list of the other commands.
- **What each install gives:** the platform wheels have everything; the source package has the Python engine
  only (codon alignment), and `coal kaks`, `coal dnds` and `coal pipeline` say so there (docs/install.md).
- The version of the development tree is `0.6.0.dev0`, so that no GPL build reports 0.5.0.

Much of what follows was found by the stress test of 2026-10-23, the full verification of 2026-10-24, two
independent reviews (2026-11-12 and 13), the first CI runs on Windows and macOS, and a review of the project as a
new user would meet it (2026-11-13). Fixes to features that were never in a release are part of that feature's
entry under "Added", not under "Fixed".

### Added

- **Ka/Ks (dN/dS) engine: `coal kaks` and `coal dnds`.** In-process C++ (no external program runs), **derived
  from the source code of KaKs_Calculator 2.0 and of PAML's yn00** (both GPL version 3; attribution in
  `cpp/src/kaks.cpp`, `NOTICE` and docs/kaks.md; this is why COAL is GPL-3.0-only from this release). Methods:
  YN00 (default), NG86, LWL85, MLWL, LPB93, MLPB and MYN, for AXT or aligned FASTA (every pair of an MSA, or
  consecutive pairs), any of COAL's genetic codes, on all threads. `coal kaks -i in.axt -o out.kaks -m YN` takes
  KaKs_Calculator 2.0's `-i -o -m -c` and writes its format (22 columns, NA rules) for these seven methods; it is
  not a drop-in for every KaKs_Calculator call: GY94 ML, MS/MA and the gamma variants are refused, and the default
  method is YN where KaKs_Calculator's is MA. `coal dnds` writes a tidy TSV with standard errors and a per-pair
  status (ok, warning for internal stop codons, saturated, no_data, error), so one bad pair never stops a run.
  NG is KaKs_Calculator's NG, whose site counting is not the textbook Nei-Gojobori one (docs/kaks.md).
  Validation on 999 orthologr pairs and 600 simulated pairs (docs/kaks.md): against a KaKs_Calculator 2.0.1 built
  to print 15 digits, Ka and Ks differ by at most 1.4e-13 (YN, MYN) or not at all (NG, LWL family) wherever both
  give a number; against PAML 4.10.10 yn00 built to print 12 digits, run one pair at a time on the columns
  without gaps, ambiguous bases or stop codons and compared with `coal dnds --yn-flavour paml`, dN and dS differ
  by at most 5e-12, which is yn00's printing precision. Every row that differs from the stock KaKs_Calculator's
  text is classified (docs/kaks.md lists the deliberate differences). C++ API `<coal/kaks.hpp>` in `coal::core`;
  the library has no global mutable state and can be called from many threads.

  Settled before the release (states that no release ever had are not listed as fixes):

  - **`--yn-flavour paml` agrees with yn00 at the K80 boundary 2P + Q = 1** (deep test of 2026-11-14, D01).
    When the transition and transversion proportions of the 4-fold (or non-degenerate) sites satisfy 2P + Q = 1
    exactly, yn00 finds K80 undefined while COAL, forming Q in KaKs_Calculator's order, found a residue of 1e-17
    and a kappa of 60 there. The paml flavour now forms Q as yn00's `GetKappa` does. The S. cerevisiae /
    C. glabrata pair of the finding: before S 194.22, N 462.78, dN 0.5024, dS 2.1598 (`ok`); after S 128.0,
    N 529.0, dN 0.4123, as yn00 (dS: yn00 prints 3.6389, COAL `saturated`, as for every pair whose Ks
    correction is undefined). All pairs with `--yn-flavour paml`: 10 of 4,388 Scer / Cgla rows change, 8 of
    4,208 Scer / Klac, 1 of 3,685 E. coli K-12 / Sakai, 0 of 2,862 Salmonella / Shigella; the 10 Scer / Cgla
    pairs disagreed with yn00 in S and N before and agree to yn00's printed digits after (weighting 1 and 0). The
    300-pair samples against yn00 (Scer / Cgla, Scer / Klac, Salmonella / Shigella, K-12 / Sakai, human / mouse,
    A. thaliana / lyrata): every value agrees to yn00's printed digits where both give one (Scer / Cgla before:
    1 pair of 300 differed). The default flavour (KaKs_Calculator's numbers) is unchanged: 0 rows differ on all
    four sets.
  - **Exit status 2 when no pair has an estimate because of the input.** `coal dnds` / `coal kaks` end with 2
    (not 0) when every result is `error` or `no_data`, with one line naming the usual cause: the input looks like
    protein sequences, or the sequences of the pairs differ in length (not aligned), and the way to a codon
    alignment. When only some pairs differ in length the exit status is 0 and a warning says how many. A run with
    any estimate, `saturated` or `not_converged` row exits 0 (docs/kaks.md "Exit status").
  - **`--min-codons N`, `--max-ds X`** (`coal dnds`, `coal codon-dnds`, `coal pipeline`; `min_codons`, `max_ds`
    in coalr): flags for estimates that rest on few codons or a saturated dS. They never drop a row and never
    change a number: the status becomes `warning` and the note says `fewer than N codons compared` / `dS above X`.
    Off by default. docs/kaks.md and docs/pipeline.md have a new section "Filtering before interpretation"
    (short alignments, coverage and identity, dS saturation, very small dS, starting thresholds).
  - **gzip-compressed input** (recognised by its first bytes): the `coal` command decompresses it for
    `coal kaks`, `coal dnds` and `coal translate`; `coal-cpp` itself says "is gzip-compressed: decompress it first
    (gunzip -k FILE)" instead of "control character". Write errors name the file and the cause
    (`write error on 'out.tsv': No space left on device`).
  - **YN and MYN are faster on close sequences, with the same bits.** The transition-probability matrix is not
    computed for a pair in which no codon differs at two or three positions (nothing reads it then): checked
    byte for byte at 17 digits on the validation sets with GCC and Clang; about half the CPU time on a set of
    closely related pairs, 10 to 20% (YN) on the 999 orthologr pairs (docs/kaks.md "Benchmark").
  - **Attribution corrected.** The documentation said the Ka/Ks engine was an independent implementation with no
    copied code. It is derived from KaKs_Calculator 2.0 and PAML's yn00 (GPL v3). The sources, `NOTICE`, the R
    package and docs/kaks.md now say so, with the upstream copyright notices and the list of changes; COAL's licence
    changed accordingly (first entry). `coal kaks` is described as KaKs_Calculator's format for seven of its methods
    with the default YN, not as a drop-in; `-m NG` is documented as KaKs_Calculator's site convention, which is not
    the textbook NG86 one.
  - **Ka/Ks: YN and MYN keep the estimate that has sites.** With no nonsynonymous sites (`GCT` against `GCA`) or no
    synonymous sites (Met/Trp only) both Ka and Ks were NA; the one that has sites is now estimated, with the
    numbers KaKs_Calculator's YN and PAML give (docs/kaks.md, difference 8).
  - **Ka/Ks: MYN gives the defined value where KaKs_Calculator prints `nan`** (a kappa class of weight 0 added as
    0 × NaN; TN93 dividing 0 by 0 when the sites have no T-C or no A-G pair; site base frequencies of a sequence
    without such sites). An identical short pair was reported as "Ka saturated"; it is Ka = 0 now (differences 9
    and 10). Checked on 4,000 random short pairs against the reference and against the reference with those two
    lines corrected.
  - **Ka/Ks: no data race.** Fisher's test took log-factorials from `lgamma()`, which writes the global `signgam`;
    the engine has its own log-gamma now (ThreadSanitizer: one race before, none after; results unchanged).
  - **`coal kaks` / `coal dnds` (and `coal translate`, `coal codon-dnds`) refuse an output file that is the input**
    (exit 2) instead of emptying it, and **an AXT record with one sequence line is that record's error**, no longer
    cut in half and compared with itself.
  - **Ka/Ks: what the comparison with KaKs_Calculator does and does not show** (docs/kaks.md, "Differences",
    rewritten; nothing changes on real genes: the output for the 999 orthologr and 600 simulated pairs is
    byte-identical). The statements "no value that both tools print differs by more than 1e-9" and "the corrected
    reference gives COAL's numbers" were made from 4,000 random pairs in which no codon differs at more than one
    position, and are now scoped to them. On 3,000 short pairs with codons drawn independently YN and MYN print
    different numbers in 33 fields (23 since the next entry but one: where the iteration does not converge COAL
    now prints none), each now explained and tested (differences 11 to 13): KaKs_Calculator drops
    the differences of codons that differ at several positions once a NaN is in its iteration and prints a Ka
    that is too small (COAL is right: it gives yn00's unweighted value, or the corrected reference's); an
    iteration that does not converge is not an estimate in either program; a saturated kind of site changes the
    path of the iteration within its tolerance. KaKs_Calculator's MYN also depends on the order of the pairs in a
    file (49 of 3,000 rows), so the comparison now runs it one pair per process.
    **Not said before: the MYN fix above ("gives the defined value where KaKs_Calculator prints `nan`") also
    changed MYN numbers that COAL used to print** on such short pairs (13 Ka and 3 Ks among the 3,000; for
    example Ka 0.0781 became 0.1673, where YN gives 0.1672): COAL had followed the reference into the lost codons.
  - **Ka/Ks: a YN or MYN iteration that does not converge gives no estimate: status `not_converged`, NA.** It was
    `ok` with a note, then `warning`, with the values of the last round printed, as KaKs_Calculator and PAML do.
    Those values are not reproducible: for the one-codon pair `AAA` / `ATG` the same source printed Ka 3.1688 on
    Linux, 3.1736 on Windows and macOS, and "saturated" with another compiler and C library (CI of 2026-11-13).
    The wandering iteration amplifies a last-bit difference of `exp` / `log`, or a fused multiply-add, into the
    third digit. Now `coal dnds`, `coal kaks`, the pipeline and coalr give `NA` for Ka, Ks, Ka/Ks, the standard
    errors and the weighted difference counts (in `coal kaks` also Fisher's p and Divergence-Time), the status is
    `not_converged` (new; between `warning` and `saturated`, and also where the last round was saturated, since
    that differs by platform too), and the note keeps the last round for inspection. The site counts are kept.
    `cluster_summary.tsv` has a new column `n_not_converged` (after `n_ok`, which does not count these pairs), the
    pipeline report has `counts.pairs.not_converged` (docs/pipeline.schema.json), the summaries on stderr name the
    count, and the C++ `Result` has `last_ka` / `last_ks`. Among the real-gene validation sets this changes one
    row (a simulated pair, MYN: Ka 0.58117 and Ks 1.43957 became NA); the 999 orthologr pairs' output is
    byte-identical. Among 3,000 random short pairs: 123 YN and 87 MYN rows. docs/kaks.md, difference 12.
  - **Ka/Ks: the same digits on every compiler.** The engine is compiled without floating-point contraction
    (a pragma in `kaks.cpp`, `-ffp-contract=off` / `/fp:precise` in CMake), so a build for a CPU with FMA (every
    arm64 Mac) no longer differs from an x86-64 one in the last digits. Eight builds (GCC 13 and 16, Clang 18,
    `-O2` / `-O3`, with and without `-march=x86-64-v3`) are byte-identical on about 67,000 rows. What remains is
    the C library's `exp` and `log`: with another one, no row of the 999 + 600 real-gene pairs changes at the
    default precision, and the last one to three of 15 digits change in 8 to 12% of them. No measurable cost in
    speed. docs/kaks.md, "Platform reproducibility", says what is guaranteed and what is not.
  - **Ka/Ks: Fisher's test with a cell that is zero up to rounding.** When the weighted differences add up to the
    number of sites (Sd = S = 2, or 2.0000000000000004), the test was made or left out by the last bit (p = 0 or
    NA, by platform). The cell now counts as zero (p = 0, as KaKs_Calculator prints). No row of the real-gene sets
    changes.
  - **`coal kaks -j1`, `coal dnds -j1`: attached values of one-letter options are read**, as the main command
    does (`-j4`, `-j=4`, `-mYN`, `-ofile`, chained flags `-qj4`), in `coal kaks`, `coal dnds`, `coal translate`,
    `coal codon-dnds` and `coal-cpp-align` (which also takes `--jobs=4` now). They were "unrecognized argument"
    while `-j 1` worked. The value of an option is never taken apart (`-o -q1` still writes a file called `-q1`).
- **Ka/Ks: the methods `Comeron` and `Li`**, the two dN/dS methods of the R package orthologr that do not need
  KaKs_Calculator (`Comeron`, alias `C95`: Comeron 1995 as libsequence's gestimator computes it, orthologr's
  default; `Li`: Li 1993 as `seqinr::kaks()` computes it). In `coal kaks` (22 columns, what the method defines),
  `coal dnds` (Li with standard errors), `coal pipeline` (`-m Comeron`), `<coal/kaks.hpp>` (`Method::Comeron`,
  `Method::Li`) and coalr (`coal_kaks()`, `coal_dnds_pairs()`, `coal_dnds_clusters()`); standard genetic code
  only. Both are ports: `cpp/src/kaks_comeron.cpp` is derived from libsequence (Copyright (C) 2003-2009 Kevin
  Thornton, GPL version 3 or later; orthologr's copy, modified by Sarah Scharfenberg and Hajk-Georg Drost under
  GPL version 2 or later) and `cpp/src/kaks_li.cpp` from `kaks.c` of seqinr 4.2-44 (`GPL (>= 2)`); the file
  headers, `NOTICE`, coalr's `inst/COPYRIGHTS` and `Authors@R` say so. Two flavours each: COAL's policy by
  default (terminal stop dropped, internal stops skipped, `NA` with status `saturated` where the correction is
  undefined), and the reference's own reading of the pair (`--comeron-flavour gestimator`, `--li-flavour
  seqinr`; in the pipeline and coalr with `--terminal-stops compare` for orthologr 0.4.3, which compares the two
  stop codons). Validation (docs/kaks.md, "Comeron and Li"): on 14,153 pairs (999 + 984 real, 600 simulated,
  570 fixtures, 11,000 random short ones) every number that gestimator's code or `seqinr::kaks` returns is the
  same double in the reference flavours, and every 999 / 9.999999 is `NA`; with orthologr 0.4.3's own
  alignments of its 1,000-gene example coalr gives orthologr's table for all 984 pairs (at the 4 digits
  gestimator prints). From the coding sequences alone it does not: orthologr 0.4.3 aligns proteins without a
  substitution matrix, COAL with BLOSUM62, and 540 of the 984 protein alignments differ. `Li` is not `LPB`
  (other site classes and pathway weights; the numbers differ). libsequence's count for codons that differ at
  two positions (it passes its arguments in another order than its function takes them) is kept in both
  flavours and documented.
  With the merge of the review fixes of 2026-11-13: `kaks_comeron.cpp` and `kaks_li.cpp` are compiled without
  floating-point contraction like `kaks.cpp` (pragma and build flags; the vendored copies in coalr too), so the
  two methods print the same bytes whether or not the target has FMA: GCC 13 and Clang 18, default and
  `-march=x86-64-v3`, with and without the build-system flag, `-O0`: 0 of 169,828 rows differ on the 14,153
  validation pairs, where 3,045 of 28,306 reference-flavour rows differed at 17 digits with contraction
  (docs/kaks.md, "Platform reproducibility"). The status `not_converged` (YN and MYN only) is in every status
  list, count and summary that Comeron and Li rows go through; they never have it.
- **One step from coding sequences to dN/dS: `coal pipeline` (docs/pipeline.md).**
  `coal pipeline pairs --query q_cds.fa --subject s_cds.fa -o out` translates both files, finds orthologs with
  DIAMOND (reciprocal best hits, or best hits; or takes a pairs file), aligns each pair's proteins (COAL's global
  aligner by default; Clustal Omega, MAFFT or DIAMOND's own alignment as options), builds the codon alignment and
  computes Ka/Ks (YN00 by default, several methods allowed). `coal pipeline clusters --cds all_cds.fa -o out`
  clusters with DIAMOND DeepClust (or takes its two-column clusters file), aligns each cluster with Clustal Omega
  (two members: the global aligner), and computes Ka/Ks for every pair of a cluster, with a per-cluster summary
  (counts by status; median, mean, quartiles, minimum and maximum of dN, dS and dN/dS). Clusters above
  `--max-cluster-size` (200) are computed against the representative, sampled reproducibly (`--pair-sample`,
  `--seed`) or skipped. Output: `dnds.tsv`, `cluster_summary.tsv`, `genes.tsv`, `failed.tsv` and a JSON report
  (`docs/pipeline.schema.json`); optionally KaKs_Calculator's format and the codon alignments. One bad gene, pair or
  cluster never stops a run (exit 2 when something failed). No process or temporary file per pair: the new native
  steps `coal-cpp translate` and `coal-cpp codon-dnds` (C++ API `<coal/dnds.hpp>`) do translation, pair alignment,
  codon threading and Ka/Ks in one threaded process fed through a pipe; output order does not depend on `-j`.
  `coal pipeline install diamond` installs DIAMOND into COAL's cache, like Clustal Omega (pinned by SHA-256). On orthologr's
  1,000-gene example (999 pairs): 2.5 s on 2 cores, every value equal to the steps run one by one and to
  KaKs_Calculator 2.0 on the same codon alignments; 19,980 pairs in 20.5 s with 41 MB.

  Settled before the release (states that no release ever had are not listed as fixes):

  - **Every FASTA record is a gene, and the run says so when that looks wrong.** Input with several transcripts
    per gene (found from `gene:`, `[gene=`, `locus_tag=`, `gene_symbol:` header fields, or from shared ID stems
    such as `AT1G01010.1` / `.2`) gets one warning with counts. `--longest-isoform` keeps the longest CDS per gene
    before the search (`--gene-map TSV` when the headers do not name the gene); dropped records are in `genes.tsv`
    with status `dropped`, and `report.json` records the choice (docs/pipeline.md "One record is one gene").
    Second review (2026-11-13): stable gene identifiers are read before gene symbols (`locus_tag`, `GeneID:`,
    Ensembl `gene:`, FlyBase `parent=FBgn`, Phytozome `locus=`, then `[gene=` / `gene_symbol:`), and records that
    share a symbol but are clearly different loci (other contig, locus tag or GeneID: several strains in one
    file, two copies of `tuf`) are not merged, with a warning. More ID styles: TransDecoder `_i1.p1`, MAKER
    `-RA` / `-mRNA-1`, maize `_T001`, RAP-DB `-01`, Phytozome proteins `.1.p`; RefSeq `NM_` / `XM_` accessions get
    a warning that the ID holds no gene. `--gene-map` wins for the records it names and the others fall back to
    the headers and IDs (it used to switch detection off for every file); a map that names no record of the
    input, or gives a record two genes, is an error (exit 1). **`--longest-isoform` never drops more than half
    of a file's records on the evidence of ID suffixes alone**: it stops with exit 1 and asks for `--gene-map`
    or the new `--isoform-rule suffix` (`map`: use the map only). docs/pipeline.md has the table of header
    styles (one test file each) and how the heuristics fail (`HSP70.1 / HSP70.2`, `scaffold7.1`). The progress
    and summary lines count "records" until `--longest-isoform` has made genes of them, and a repeated ID is
    counted as not usable in both.
  - **A file is complete or absent** (second review, 2026-11-13). Every table and alignment file of the pipeline
    is written as `NAME.tmp`, created anew, and renamed when it is complete; before, only `dnds.tsv` was, and a
    stage that failed on a full disk left a cut-off `codon_alignments.fasta`, `protein_alignments.fasta` and
    `dnds.kaks` under their final names. After a failed or interrupted stage none of them is there, the report
    lists only the files that exist, and the `.tmp` files are removed (the next run removes what a killed run
    left). A symbolic link in the output folder under an output's name is replaced, never written through
    (`dnds.kaks -> another file` used to overwrite that file). `coal-cpp codon-dnds` has `--new-files` for this
    and writes each per-cluster file through `NAME.fasta.tmp`. docs/pipeline.md "Output".
  - **A hint for genes of another genetic code** (second review, 2026-11-13). Mitochondrial genes in a nuclear
    CDS set got `internal_stop` warnings and estimates from mistranslated codons, and nothing said why. When at
    least two genes of a set have two or more internal stop codons under `-c` and none under another table COAL
    supports, the run prints one warning per set (how many genes, which table, three examples, the two-run
    recipe); `report.json` has it as `genetic_code_hints`. docs/pipeline.md "Genes of another genetic code":
    mitochondrial and plastid genes, the recipe, and that `Comeron` / `Li` are standard-code only. A per-gene
    code option (`--code-map`) was evaluated and not built: the engine has one code per run; the design is in
    docs/development.md.
  - **Smaller things from the second review** (2026-11-13). `--diamond-args`, `--deepclust-args` and
    `--realign-args` refuse the DIAMOND options COAL sets itself (`--out`, `--db`, `--query`, `--outfmt`,
    `--threads`, ...; exit 2 before anything runs: `--diamond-args '--out x.tsv'` used to send the hits elsewhere
    and end with an unrelated error), all three are documented (`--realign-args` was not), and `realign.tsv`
    always has its header line. The hint about IDs of `--pairs` / `--clusters` that match no record comes when
    half of one side is unknown or no pair can be made (it needed 90% of one side). Messages name the file the
    user gave, not its decompressed copy (`coal pipeline`, `coal kaks`, `coal dnds`). The default `-j` is one
    thread per CPU the process may use (its affinity: `taskset`, SLURM, cgroups), and `coal-cpp` is never given
    more threads than that either (`report.json` `threads`: `requested`, `usable_cpus`, `engine`, `tools`).
    DIAMOND's version is read as the first dotted number of `diamond version` (`2.1` is 2.1.0, not "too old"),
    and a version that cannot be read gives one warning. A gzip input has `sha256_uncompressed` next to `sha256`
    in the report.
  - **`--min-identity P`**, the third flag next to `--min-codons` and `--max-ds` (`coal pipeline`,
    `coal codon-dnds`; `min_identity` in coalr): a pair whose protein identity (the column `prot_identity`) is
    below P percent gets status `warning` and the note `protein identity below P%`; no row is dropped, no number
    changes. docs/pipeline.md says why: a cluster can hold pairs that are not related at all, and they get
    status `ok`.
  - **Documentation from the second review:** docs/pipeline.md "How much the aligner matters" (150 *Arabidopsis*
    pairs: 120, 90 and 75 of them have the same values as the default aligner with MAFFT, DIAMOND's alignment and
    Clustal Omega; up to 0.43 in dN/dS for single pairs); docs/kaks.md says what `P-Value(Fisher)` tests and
    that `SE_dN` / `SE_dS` are standard errors of dN and of dS (NG, YN, MYN, Li), with none for the ratio; README
    example 5 (the orthologr fork) passes `sensitivity_mode = "default"`, so that it finds the 12 pairs of
    example 3 (orthologr's default, DIAMOND `--fast`, finds 8 on the example files).
  - **Cluster mode does not depend on the order of the input.** DIAMOND DeepClust is order-dependent, so it is
    given the proteins sorted by ID: the same sequences give the same `clusters.tsv` whatever the order of the
    `--cds` files and records. (Development builds before 2026-11-13 clustered in input order; their clusters can
    differ from this release's.) docs/pipeline.md "What the search and clustering defaults do" documents the
    best-hit rule, DIAMOND's sensitivity, a worked comparison with orthologr's fork, and why default DeepClust
    clusters hold many saturated pairs.
  - **Provenance in `report.json`:** start and end time (UTC), platform, Python and COAL versions, the thread
    counts as resolved, every input file with size and SHA-256 (`--hash-inputs auto|always|never`; above 2 GiB the
    hash is skipped with the reason), the exact DIAMOND and coal-cpp command lines, the aligner's command
    template, the citation keys (`cite`), the run's warnings. The tables are unchanged by it.
  - **Robustness:** gzip-compressed inputs are decompressed into `work/inputs/`; a full disk (any operating-system
    error) ends the run with one line, exit 1 and a report that is complete or absent; a `--pairs` / `--clusters`
    file whose IDs match (almost) nothing gets one line with example IDs of both sides; DIAMOND older than 2.1.0
    stops cluster mode with the reason (pair mode warns; tested with 2.2.8); `-j` above the CPU count is clamped
    for DIAMOND and the aligners; a tool ended by a signal is reported as "killed by signal 9 (SIGKILL; out of
    memory?)"; the same file as `--query` and `--subject` is a warning.
  - **`coal pipeline ... --keep-protein-alignments`** writes the protein alignments the codon alignments are made
    from: `protein_alignments.fasta` in pair mode (records 1+2, 3+4, ...), one file per cluster in
    `protein_alignments/` in cluster mode, with the same records, order and file names as `--keep-alignments`
    writes, so a codon row is three times its protein row, column for column. Until now the protein alignments
    could not be obtained from the pipeline (review of 2026-11-12). The folder is recorded in `work/outputs.json`
    (`protein_alignments_dir`), so the next run removes that run's files and nothing else; a `protein_alignments/`
    that no COAL run wrote is never written into. The engine step has `coal codon-dnds --protein-out FILE` and
    `--protein-dir DIR`; the report's `outputs` has `protein_alignments` (docs/pipeline.md, the schema).
  - **`coal pipeline`: one run per output folder** (`<out_dir>.coal-lock`; a second run is refused), **an earlier
    run's files are removed when a run starts** (no `dnds.kaks` or `codon_alignments.fasta` of another run next to
    a new `dnds.tsv`), **an interrupt ends the engine and is reported** (SIGINT/SIGTERM to the Python process kill
    `coal-cpp codon-dnds` and the aligners; `report.json` gets `"status": "interrupted"`, a new report field, and
    the exit code is 130 or 128 + signal), **a protein file is explained before DIAMOND runs**, and **a non-ASCII
    sequence ID is that gene's failure** (`id_not_ascii`), not the end of the run (docs/pipeline.md).
  - **`coal pipeline` no longer deletes its inputs or files it did not write.** The clean-up at the start of a run
    (previous entry) removed every file with an output's name, also when it was this run's input or the user's own:
    `--pairs out/orthologs.tsv -o out` failed with "cannot open" after deleting the pairs and the earlier
    `dnds.tsv`, and `clusters --clusters proj/clusters.tsv -o proj` deleted the clusters file and every
    `proj/codon_alignments/*.fasta`. Now a run removes only what an earlier run is known to have written (the
    list each run keeps in `work/outputs.json`, or the outputs named in an older `report.json`), never one of its
    inputs, and it never writes onto an input: `--pairs out/orthologs.tsv` is read and left as it is, any other
    input under a name the run writes is refused before anything changes, and so is `--keep-alignments` into a
    `codon_alignments/` folder that holds files COAL did not write (docs/pipeline.md).
  - **`coal pipeline`: one rule for text, on every system** (docs/pipeline.md, "Text encoding"). Every table is
    read and written as UTF-8 with `\n` line ends, never with the system's default encoding or line ends; sequence
    IDs are bytes that are not re-encoded. Bytes that are not valid UTF-8 are shown as `\xhh` (an ID in Latin-1:
    `g\xe8ne`) by `coal-cpp translate` and by the Python layer alike: before, `coal-cpp` wrote `?` for each such
    byte and Python U+FFFD, so two different IDs could look the same. `coal-cpp-align` writes its alignment to
    standard output in binary mode (`\n` on Windows, as its `-o` file). The Windows CI failure itself was the
    test reading `failed.tsv` with the default encoding (cp1252); the tests now run the pipeline under an ASCII
    locale, with a cp1252 console and with Python's `-X warn_default_encoding` as an error, and compare every
    output byte with the UTF-8 run.
  - **Deep test of 2026-11-14 (D01-D08).** **The same IDs in two `--cds` files no longer drop a species
    without a word** (D05). `clusters --cds a.fa b.fa` with IDs repeated across the files (BRAKER / AUGUSTUS
    `g1.t1` in every species) kept the first file's record of each ID and left out the other's: 300 + 300
    Salmonella / Shigella CDS named `gene0`..`gene299` in both gave 289 clusters and 12 pairs instead of about
    300 cross-species pairs, exit 2, and 300 failures saying "the ID occurs more than once in the file" (each
    occurred once per file). Now such a run stops before anything runs (exit 1) and names both files, the count
    and examples; the new `--prefix-ids` puts each file's name and `|` before its IDs (`q_sameid|gene0`; the
    report's `inputs` give each file's `id_prefix`), which gives 600 usable records, 291 clusters and 352 pairs
    on the same files. `coal-cpp translate` with several files says "the ID occurs in two input files, A and B"
    for such a record (within one file the wording is as before).
    **The genetic-code hint no longer names a wrong code** (D02). It named the code that made the most flagged
    genes stop-free, and codes that also read TAA / TAG "clean" frameshifted pseudogenes: it fired in 6 of the 9
    real pair datasets of the deep test and was wrong in each (S. cerevisiae mitochondrial genes: 14, should be
    3; C. glabrata / K. lactis mitochondria: 4, should be 3; human and mouse Ensembl: 4, should be 2; E. coli K-12,
    A. lyrata NCBI and Ensembl mouse all-transcript files: 14 for pseudogenes and NMD transcripts). Now
    pseudogenes (`[pseudo=true]`, `*pseudogene` biotype) and frameshift-length CDS are not evidence, genes
    scattered over nuclear chromosomes are not either when the headers name the molecule, the code whose extra
    stop or unused codons the genes avoid wins (AGA / AGG for 2, CGA / CGC for 3) before the one with the fewest
    unneeded stop reassignments, and when the sequences cannot tell the codes apart all candidates are named and
    none is recommended (report: `clean_under` null, new `candidates` and `evidence`). After, on the same 9
    datasets: Scer / Cgla and Scer / Klac: code 3 for both sets (17, 8 and 7 genes); human / mouse whole Ensembl
    sets: code 2 (12 genes each); no hint for Salmonella / Shigella, K-12 / Sakai, the A. thaliana / lyrata and
    human / mouse subsets, A. thaliana / lyrata whole genomes, the human / mouse all-transcript files. Also
    A. gossypii: 3 (was 4); Y. pestis (78 genes, almost all `[pseudo=true]`): none (was 14); D. melanogaster
    Ensembl: "4 or 5" (was 14; 5 is right). A synthetic set of four nuclear genes given two internal TGA each
    now gets the candidates 3, 4, 5, 9, 10, 13 and 21 instead of "code 4".
    **A byte order mark or classic Mac line ends no longer lose genes** (D03). `coal-cpp translate` skipped a
    UTF-8 byte order mark and split lone `\r` line ends, the Python-side readers did not: with a BOM the first
    query gene was "not in the query file" (300 Salmonella CDS: 1 pair error, exit 2), with `\r` line ends every
    pair was "the CDS is empty" (300 errors), the first line of a `--gene-map` was silently ignored (4 records
    dropped instead of 5, exit 0), and the first member of a `--clusters` file was dropped with an invisible U+FEFF
    in its cluster ID (exit 2). Now the FASTA index, the header reader of `--longest-isoform`, the `--pairs`,
    `--clusters` and `--gene-map` readers and `coal --batch` job files read such files as coal-cpp does; on the
    deep test's files every output is byte-identical to the run on the plain file (exit 0).
    **Messages** (D04, D06, D07, D08; no number changes). A `coal --batch` job line that starts with `coal` (as
    the README's wording suggested) is refused before any job runs, naming the fix, instead of 400 errors "cannot
    open 'coal'"; the README says "without the program name". A pipeline run with no pair says so ("no pairs:
    DIAMOND found no reciprocal best hit ..." / "none of the N clusters has two usable members"; exit status
    unchanged), a `--pairs` file whose columns look swapped and one with repeated pairs get one warning each (the
    pairs are still computed as given), and `$COAL_CLUSTALO` / `$COAL_MAFFT` / `$COAL_DIAMOND` set to something
    that is not a file is named as the cause instead of "not found (PATH, $COAL_CLUSTALO, COAL's cache)".
    `coal dnd x.fa` (and `coal-cpp dnd x.fa`) says "did you mean `coal dnds`?". `coal pipeline install mafft`
    installs MAFFT 7.526 (pinned, as Clustal Omega 1.2.4 is; it took whatever Bioconda had) and prints "MAFFT
    7.526 installed". When more than 30% of the pairs are `saturated` or `not_converged` the run says that the
    sequences are too divergent for dS and points to "Filtering before interpretation" (S. cerevisiae /
    C. glabrata: 3,775 of 4,388, 86%).
  - **Retest of 2026-11-14 (R01-R04).** **The genetic-code hint names the code of a whole set read with the
    wrong `-c`** (R01). *Mycoplasma genitalium* (524 CDS, code 4) with `-c 1` or `-c 11` got no hint at all
    although 376 genes have internal TGA: the rule that keeps only molecules with at least half of their records
    as evidence dropped its one chromosome (253 of 524). *Tetrahymena* with its records shuffled got "81 of the
    26,927 genes ... mitochondrial or plastid genes in a nuclear set? ... the rest with -c 1 and these with -c 6"
    (the first 5,000 flagged genes looked at, spread over scaffolds), i.e. advice to keep 26,846 genes under the
    wrong code. Now, when at least half of a set's records have internal stops, the hint is about the whole set
    ("most genes (376 of the 524) have internal stop codons under genetic code 11, and 253 of the 264 with two
    or more have none under NCBI genetic code 4 ... The whole set looks like code 4 ... Run again with -c 4"),
    the molecule rule does not apply, a molecule with 150 records or more is a genome (only the codes of whole
    genomes are candidates: 1, 4, 6, 10, 11, 12, 15), and of codes that tie those whose start codons do not fit
    a tenth of the genes drop out (Mycoplasma starts 34 of 264 genes with GTG / TTG / TTA, which code 10 does
    not allow). A set with more than 5,000 flagged genes is judged on 5,000 spread evenly over the file, not the
    first ones. `report.json`: new `records` and `whole_set` in `genetic_code_hints`. The organelle wording and
    the two-run recipe are kept for sets in which fewer than half of the records have internal stops. Before
    (0153e3b) and after, the warning of `coal pipeline pairs` for the query set, on every real CDS set of the
    deep test, the retest and the D02 evaluation (43 runs); 33 are unchanged (listed below the table):

    | Set | `-c` | Before | After |
    |---|---|---|---|
    | *M. genitalium* RefSeq (524) | 1 | none | whole set: code 4 (376 of 524 with internal stops) |
    | *M. genitalium* RefSeq | 11 | none | whole set: code 4 |
    | *T. thermophila* RefSeq (26,996) | 1 | "4,896 of the 26,927 genes": code 6, organelle wording | whole set: code 6 (4,966 of 5,000 looked at) |
    | *T. thermophila*, records shuffled | 1 | "81 of the 26,927 genes": code 6, organelle wording | whole set: code 6 (4,971 of 5,000) |
    | *P. tetraurelia* RefSeq (39,642) | 1 | "4,850 of the 39,534 genes": code 6, organelle wording | whole set: code 6 (4,995 of 5,000) |
    | human Ensembl (19,746) | 2 | "4,988 of the 19,311 genes": codes 1, 11, 12, "mitochondrial genes in a nuclear set?" | whole set: codes 1, 11 or 12 |
    | *S. cerevisiae* mitochondrial genes alone (19) | 1 | code 3, organelle wording | whole set: code 3 |
    | *C. glabrata* mitochondrial genes alone (11) | 1 | code 3, organelle wording | whole set: code 3 |
    | *K. lactis* mitochondrial genes alone (9) | 1 | codes 3, 4, 5, 13, organelle wording | whole set: codes 3, 4, 5 or 13 |

    Unchanged: *S. cerevisiae*, *C. glabrata*, *K. lactis*, *A. gossypii* whole genomes (`-c 1`: code 3 for
    their mitochondrial genes); human and mouse Ensembl (`-c 1`: code 2, 12 genes each); *D. melanogaster*
    Ensembl (twice) and RefSeq (`-c 1`: codes 4 or 5, 11 genes); no hint for *Salmonella* (twice), *Shigella*
    (twice), *Y. pestis* (twice), *E. coli* K-12 and Sakai (`-c 11`), *A. thaliana* (whole, 1,000, TAIR10
    RefSeq), *A. lyrata* (whole, 1,000), the human giant-gene and exon sets, the human and mouse subsets and
    all-isoform files; none with the right code (*M. genitalium* `-c 4`, *T. thermophila* and *P. tetraurelia*
    `-c 6`, *S. cerevisiae* mitochondrial genes `-c 3`), and none where the code reads fewer codons as stops than
    the genes' own, so that no gene has an internal stop to explain (*S. cerevisiae* genome `-c 3`,
    *D. melanogaster* `-c 5`). No run names a wrong single code.
- **Protein alignment for the codon step (docs/alignment.md).** `<coal/align.hpp>` and the `coal-cpp-align` tool: a
  global Needleman-Wunsch/Gotoh aligner (BLOSUM62/45/80, PAM250, affine gaps; EMBOSS needle defaults, a pwalign
  preset), full length (overhangs gapped, never cut off), AVX2 row loops picked at run time, multithreaded; optimal
  scores equal Biopython's and EMBOSS needle's. The expansion of DIAMOND `qseq_gapped`/`sseq_gapped` hits to
  full-length alignments without realigning. Python: `coal.align` (the same, with a pure-Python fallback) and
  `coal.msa` (Clustal Omega / MAFFT for clusters: find, install into COAL's cache, run, batch; stop codons handled).
  Clustal Omega (GPL-2.0-or-later) is installed on demand, not bundled (bundling would be licence-compatible now
  that COAL is GPL-3.0-only; it is not done). Benchmark on orthologr's 1,000-gene set:
  `benchmarks/results/align_layer_linux_20261111/`. `coal::core` now links `Threads::Threads`.

  Settled before the release (states that no release ever had are not listed as fixes):

  - **The platform wheels ship `coal-cpp-align` next to `coal-cpp`** (`coal/bin/` inside the package), so
    `coal.align.align_pairs` runs the native aligner after `pip install coal-align` on Linux, macOS and Windows
    instead of the pure-Python kernel (about 1,000 times slower). `packaging/wheels/build_platform_wheel.py` takes
    `--align-binary`; `build_linux_binary.sh` and the release workflow build and check both binaries;
    `packaging/wheels/check_wheel_binaries.py` runs in the environment the wheel was installed into (both lookups
    return the wheel's binaries, one alignment, one `coal pipeline pairs` run with an empty `PATH`).
    `coal-cpp-align --version` is new. Verified on Linux x86-64 with the repository's build script outside the
    container; the macOS and Windows steps of the workflow are changed the same way and have not run yet.
  - **An interrupt can no longer leave a Clustal Omega / MAFFT process behind.** `coal.msa` started an aligner
    and registered it a moment later; an interrupt in between left it running while the run reported that nothing
    was. Starting and registering now happen under the lock that `kill_running()` takes, and after
    `kill_running()` the interrupted batch starts no further aligner (its clusters get the error `interrupted`);
    the next batch runs normally.
  - **Very long pairs no longer take 600 MB or more** (deep test of 2026-11-14). The traceback is held at most
    256 MiB at a time (`AlignParams::traceback_bytes`); a larger pair is aligned with checkpoints, recomputing
    one block of rows at a time with the same code, so its alignment is identical. Human / mouse Ensembl pipeline
    (17,118 pairs with titin): peak 711 MB before, 227 MB after, 49.9 s → 53.2 s on 2 cores; `dnds.tsv` and
    the codon alignments byte-identical, and identical to the deep test's run. Forced on every pair (test hook,
    blocks of 1 and more rows): Salmonella / Shigella and S. cerevisiae / C. glabrata `dnds.tsv`, codon and
    protein alignments byte-identical.
- **R package `coalr` (packaging/r/coalr, docs/r.md).** COAL's C++ core compiled into an R package with Rcpp:
  `coal_dnds_pairs()` (CDS of ortholog pairs to dN, dS, dN/dS: the columns of `dnds.tsv`), `coal_dnds_clusters()`
  (every pair of each cluster and the per-cluster summary of `cluster_summary.tsv`; protein alignments from given
  rows, an aligner function such as `coal_msa_aligner("clustalo")`, or the built-in pairwise aligner),
  `coal_kaks()`, `coal_codon_align()` (PAL2NAL's job), `coal_align_pairs()` and `coal_translate()`. Everything runs
  inside the R process: no program is started, no file is written, worker threads never touch R, a user interrupt
  stops the run within about 0.1 s, and bad data is a row with a status, not an error. `coal_pipeline_pairs()` /
  `coal_pipeline_clusters()` run the whole `coal pipeline` through reticulate when Python with `coal-align` is
  there. The package holds a copy of `cpp/src` (`scripts/sync_r_sources.py`; `tests/test_r_bindings.py` checks it
  is current), so `remotes::install_github("Plumat-COAL/COAL", subdir = "packaging/r/coalr")` works. On
  orthologr's 1,000-gene example its 999 pairs equal `coal pipeline pairs` in every column and take 1.5 s on one
  thread. `R CMD check --as-cran`: no error or warning. Built and tested on Linux only so far.

  Settled before the release (states that no release ever had are not listed as fixes):

  - Version 0.6.0 (0.6.0.9000 in the development tree); `min_codons` / `max_ds` in `coal_dnds_pairs()`,
    `coal_dnds_clusters()` and `coal_kaks()`; its messages name COAL 0.6 and the release wheel (not PyPI).
  - **coalr: a path-like string that names no file is the error "file not found"**, not a sequence (docs/r.md).
- **`coal --batch --batch-parallel N`: N jobs at a time.** Worker processes each import COAL once and take jobs
  from a queue (`0` = one per CPU core; the default is still one job after the other), every job's output
  byte-identical to running it alone; a job's messages are printed together when it ends, so the lines of two
  jobs never interleave, and the batch report keeps job order. Ctrl-C and SIGTERM end every worker and remove
  half-written outputs; a worker killed from outside fails its job and the batch goes on. coal-cpp refuses the
  option like the other batch options. 64 jobs of 2,000 × 300 (16-core WSL 2): 4.15 s one after the other, 1.20 s
  with 4 workers, 0.71 s with 16; tiny jobs gain nothing (docs/usage.md).

  Settled before the release (states that no release ever had are not listed as fixes):

  - **`--batch-parallel`: a worker killed while it sends a result no longer hangs the batch.** Jobs and results
    went through two queues shared by all workers; a worker killed in the middle of writing a result died holding
    the queue's lock, with half a message in the pipe, and the batch waited for ever. Every worker now has a pipe
    of its own, the batch process hands out the jobs and watches the workers themselves, and a killed worker fails
    its job with the documented one line (`error: the job's worker process died (out of memory?)`) while the
    batch goes on; the dead job's lock file is removed. Nothing else about `--batch` changes.
- **What to cite for what you ran.** `coal cite` prints a table (codon alignment: PAL2NAL; each Ka/Ks method: its
  estimator and the program COAL's code is derived from, KaKs_Calculator 2.0, PAML, libsequence or seqinr; the
  pipeline: DIAMOND, DeepClust, Clustal Omega, MAFFT) and the references; `coal cite YN` those of a method;
  `coal cite out/report.json` those of one pipeline run. The same table is in README.md, docs/kaks.md,
  `CITATION.cff` (references) and on the website. Every cited work has its DOI (second review, 2026-11-13: each
  reference was checked against the publisher's or the institution's page): `coal cite` prints it as a
  `https://doi.org/` link, `CITATION.cff` has `doi:` fields, the documents and the website link it. COAL itself
  still has no DOI. **DIAMOND DeepClust is cited by its paper**, Buchfink, Barbé, Ashkenazy, Reuter, Kennedy &
  Drost (2026), Nature Methods 23:724-727 (key `buchfink2026`); the key of its preprint, `buchfink2023`, is no
  longer listed but still resolves (`coal cite buchfink2023`, and a report that names it), followed by the paper.
  The seqinr reference has the book's full title.
  On the website's Cite page the table wraps and, on a phone, stacks the two cells of each row (it showed only
  its first column at 390 px).
- **A documented Python API for Ka/Ks and the pipeline.** docs/python-api.md has a new section "dN/dS and the
  pipeline": `coal.dnds_pipeline.run_pairs` / `run_clusters` (signatures, `PipelineResult`, the exit codes, and
  that a failed run returns exit code 1 with `report["error"]` instead of raising; `TypeError` / `ValueError` for
  a wrong call), `codon_dnds`, the column constants `PAIR_COLUMNS` / `SUMMARY_COLUMNS` / `FAILED_COLUMNS`,
  `METHODS`, `SCHEMA`, `find_diamond` / `install_diamond`, `coal.citations`, and the entry points of `coal.align`
  and `coal.msa`, with an example on `examples/dnds` that the tests run. The package ships `py.typed` and the
  public functions have accurate type hints. The versioning promise now covers the column names and order of
  `dnds.tsv`, `cluster_summary.tsv` and `failed.tsv` and the report schema id `coal.pipeline/1`. `import coal`
  exports what it did (no facade); its docstring no longer calls the rest of the package internal.
- **coalr: a vignette, and one name for the methods.** `vignette("from-two-cds-files-to-dnds")` goes from two
  CDS files to a filtered dN/dS table on an example shipped with the package (in process, in seconds, no
  network; built by `R CMD check`). `method` and `methods` are accepted by every function (a `method =` given
  to `coal_pipeline_pairs()` used to be ignored), `min_identity` is the R form of `--min-identity`, and `?coalr`
  and docs/r.md have one table of the command line's options and the R arguments, with the thread defaults.
- **`coal --help` lists every command, and `coal help [COMMAND]` works.** The help of `coal` (and of `coal-cpp`)
  ends with one line per command: `kaks`, `dnds`, `pipeline`, `translate`, `codon-dnds`, `cite`, `help`,
  `--cpp-path`. `coal kaks --help`, `coal dnds --help`, `coal translate --help` and `coal codon-dnds --help` print
  their help in an install without `coal-cpp` too. The "needs coal-cpp" message of the Ka/Ks commands and the
  pipeline says what the install is (the Python engine: codon alignment only) and names the platform wheel.
- **`examples/dnds`:** two small synthetic species with 12 known ortholog pairs (a range of dN/dS, all `ok`), a
  ready codon alignment, and the true pairs, for `coal pipeline`, `coal dnds` and `coal kaks`
  (`examples/make_dnds_example.py`; the older examples give only saturated pairs in `coal dnds`). README,
  docs/usage.md and the website use it.
- **GUI: a "Colour residues by class" setting.** Off, the alignment view and its downloadable page show the
  letters only (the mismatch outlines, frameshift hatching and the key's mark entries stay): for colour-blind
  readers who prefer no tints, and for printing. On by default; saved with the other settings.

### Changed

- **The app and its MCP server say that they cover codon alignment only.** The app's README, the GUI's About
  section and the MCP server's instructions now state that Ka/Ks (dN/dS) and the pipeline are in the `coal` command
  line (`coal dnds`, `coal kaks`, `coal pipeline`), so that neither a person nor an AI agent looks for them in
  the app. No tool changed. docs/mcp.md has the proposal for the MCP tools `dnds_pairs` and `pipeline_pairs`
  (planned, not built).

- **Licence: GPL-3.0-only from the next release.** COAL (the Python and C++ engines, the `pal2nal` drop-in, the app
  and everything else in this repository) is now under the GNU General Public License, version 3 only (`LICENSE`,
  `NOTICE`; package metadata, `CITATION.cff`, the Bioconda recipe, the website and the Windows app say so). **0.5.0
  and every earlier release remain under the MIT License**: nothing already published changes. Why: the Ka/Ks
  (dN/dS) engine that the next release adds is derived from KaKs_Calculator 2.0 and from PAML's yn00, which are
  both GPL version 3 (without an "or later" clause), so COAL as a whole has to be distributed under GPL-3. All
  contributors agreed to the change. What it means for users: running COAL (the commands, the app, a pipeline step)
  does not affect the licence of your data, results, scripts or workflows. A program that includes COAL's code or
  links its library (the C++ `coal::core`, `import coal` in Python) and is distributed must be under a
  GPL-3-compatible licence; programs that only start `coal` / `coal-cpp` as a separate process are not affected.
  The third-party notices already in the repository (the website's font, the real-data test sets) are unchanged.
- **The development tree has the version `0.6.0.dev0`** (`coal --version`, `coal-cpp --version`, reports), so a
  build of it is never taken for the MIT release 0.5.0; coalr is 0.6.0.9000. The release sets 0.6.0
  (`scripts/release_version.py`). The run report's schema accepts such a version.
- **`coal --help` is longer:** after the options it lists the other commands (see "Added"); a script that
  compares the help text byte for byte sees a difference. The options and their help lines are unchanged.
- **Two runs writing the same `-o` file at the same time are refused** instead of sharing one `<out>.coal-partial`
  file, which could end a run with exit 0 and an empty or another run's output. A run with `-o` keeps
  `<out>.coal-lock` while it runs; a lock left by a killed run is replaced once its process is gone (a lock
  file that is still empty is one being created by another run, not a stale one: the second iteration of the
  stress test caught two runs reading each other's empty lock and both going ahead).

  Also part of the lock as released:

  - **The output lock (`<out>.coal-lock`; `coal -o` and `coal pipeline`) holds under simultaneous starts, through
    links and across machines** (the Python side). Several runs that started together and found a stale or empty
    lock could all take it; a stale lock is now replaced by one run at a time (`<out>.coal-lock.takeover`, held
    for an instant) and only if it still is the lock that was judged stale. The lock is next to the real output,
    so a symbolic link to a pipeline folder no longer has a lock of its own. The lock line is now
    `<pid> <host> <process start> <token>` (the pid stays first, which is all `coal-cpp` and earlier versions
    read): a lock from another machine counts as held and the message names the machine, a process id that is
    alive but started at another time (Linux) is not the holder, a second call into the same output from the same
    process is refused, and releasing never removes a lock that is not the run's own (docs/pipeline.md).
  - **A refused `--resume` (or a signal) no longer leaves `<out>.coal-lock` behind, both engines.** The lock was taken
    before the checkpoint set-up and only released by the block after it, so an error raised in between (the inputs
    changed since the checkpoint, a checkpoint of the other engine) or a signal there left the lock, and the next run
    reported a stale lock of a dead process. Found by the new `--batch-parallel` SIGTERM test on a macOS CI runner;
    a test refuses a resume in both engines and checks the lock is gone.
  - **`coal-cpp -o`: the output lock of the C++ engine follows the Python engine's rules.** `<out>.coal-lock`
    holds `<pid> <host> <start> <token>` in both engines; it is created with its content in one step (never seen
    empty), a lock written on another machine is held and never replaced, a process id that was reused (another
    start time; Linux) is not the holder, and a stale lock is replaced by one run at a time (`<lock>.takeover`),
    which removes it only if it is still the one that was judged stale. Before, several `coal-cpp -o` runs that
    found a stale lock at the same moment could each remove "the stale lock", the later ones removing the lock an
    earlier one had just taken, and both wrote the output; and a lock from another machine counted as stale.
    Tests: the race step by step (the unguarded code fails it), 200 trials of six simultaneous starts onto a stale
    lock (one holder each time), Python and C++ starts mixed, each engine refused by the other's lock. The
    Windows code path (process handle, `CreateHardLink`) is compiled by CI and not exercised by these tests there.
    **On Windows** the first CI run of these tests found the guard `<lock>.takeover` left behind after simultaneous
    starts (one holder each time, as required): Windows does not remove a file that another process has open, the
    runs that lost the race were reading the guard at that instant, and the removal was tried once. Both engines
    now try every removal of a lock, guard or temporary file again for about 2 s, remove the guard on every path,
    and clear a guard that is left all the same when its process is gone or it is older than 15 s (before: 30 s,
    the time a run waits). The host name in the lock has one documented source and spelling in both engines and is
    compared without regard to case on Windows (`docs/pipeline.md`). Tests on every system: a hook makes removals
    fail as they do on Windows.
- **Installers: what is in COAL's cache is checked before it is run.** Each installed binary gets a
  `<binary>.coal-manifest.json` with its SHA-256, and a cache binary is used only while it has that hash.
  **Tools installed into the cache by COAL 0.5.0 or earlier have no manifest and are no longer used
  automatically** (that includes a Clustal Omega fetched without verification): COAL says so once and
  `coal pipeline install <tool> --force` installs them again. Tools on `PATH` or given by `COAL_CLUSTALO`,
  `COAL_MAFFT`, `COAL_DIAMOND`, `--diamond` are unaffected. Also: downloads go to a new randomly named file
  instead of the fixed `<name>.part` (two installs collided there, and a link at that name was written
  through), the hash is taken from the file on disk, one install per tool runs at a time and the others wait, an
  archive that was cut short (possible with `--allow-unverified`) is a clear error instead of an `EOFError`,
  and a zip member name with `:`, an absolute path, `..` or a reserved device name is refused
  (docs/alignment.md).

  Also:

  - **Installers verify before they run anything.** Downloads are kept only if their SHA-256 is the pinned one:
    micromamba is one fixed release (2.10.0-0) with pinned hashes, not `latest`; the DIAMOND 2.2.8 archives are
    pinned. Clustal Omega is installed from Bioconda only; see the entry on the Clustal site below
    (docs/alignment.md).
  - **The download of Clustal Omega from the Clustal site is removed.** That site's address now serves an
    unrelated "people search" site and its three pre-compiled binaries (Linux, macOS, Windows) answer "not found"
    (checked from a browser on 2026-11-12). There is nothing to pin, and anything that appears there later would not
    be from the Clustal authors, so COAL no longer holds an address of that site: `install_clustalo(method=
    "official")` is refused with this explanation, and `sha256=` / `allow_unverified=True` (`--sha256`,
    `--allow-unverified`) cannot bring it back (with Clustal Omega they are an error). Linux and macOS install it
    from Bioconda (conda / mamba / the pinned micromamba) as before. **On Windows there is now no automatic Clustal
    Omega install**: the message names what works (WSL; a `clustalo.exe` of your own on `PATH` or in
    `COAL_CLUSTALO`; `--aligner mafft` with a MAFFT you installed, for which COAL has no verified Windows download
    either). Bundling a self-built `clustalo.exe` is the realistic way to support Windows and is not done yet
    (docs/alignment.md, "Licences").
- **`pal2nal` without arguments is a usage error:** the usage text goes to standard error with exit 1 (as the
  drop-in's help says and as PAL2NAL's own usage goes to stderr), not to standard output with exit 0, so a pipeline
  step that forgot its files no longer gets the usage text as its "alignment". `pal2nal -h` is unchanged (stdout, 0).
- **`-q` without `--report`/`--json` keeps no per-sequence entries** and writes none to the `--resume` entries
  file (a million genes of unrelated DNA wrote 55 GB of it), so memory and disk stay flat on data where every
  codon mismatches (the summary lines, the hint and the report, which read them, are not produced under `-q`).
- **The genetic-code hint no longer re-tests genes of more than 100,000 residues** (28 s of a 30 s run on one
  1,000,000-residue gene with one wrong codon; `-q` already skipped it).
- **Messages name at most three input files** ("a.fa, b.fa, c.fa and 997 more files") instead of every path.
- **coal-cpp: the exact-fit search is linear in the DNA** (a string search over the translated frames) instead of
  trying every start codon by codon: a 200,000-residue low-complexity gene behind a near-copy of itself took 138 s
  (an hour at 1,000,000 residues), now well under a second, as in the Python engine. Same answer everywhere
  (the differential tests and the engine fuzz).
- **GUI: the alignment view's residue-class colours are readable for colour-blind users.** The seven classes keep
  their Okabe-Ito hue families but now sit on an even lightness ladder, so they stay apart under protanopia,
  deuteranopia and tritanopia and in greyscale (the pastel tints of equal lightness merged there; verification
  2026-10-24, V-17); the text stays at 5:1 or better on every tint, in the light and the dark theme. A test holds
  the simulated distances. The letters, outlines and legend carried the meaning before and still do.

### Fixed

In commands that 0.5.0 has (both engines unless said):

- **A helper process killed before it had started no longer ends the run with exit 141 and no message.** With
  the spawn and forkserver start methods (macOS, Windows, Python 3.14 on Linux) the first child of a run with
  `-j` is Python's resource tracker; if it was killed at once, the pool's first write to it failed with a broken
  pipe, which COAL took for a closed standard output. It is now the usual
  `coal: error: a worker process died (out of memory?)` (exit 1). (This was so before the review; `--batch-parallel`
  reports a worker that cannot be started the same way.)
- **Python engine: a worker killed while it sends a result no longer hangs the run.** The pool's management
  thread then waits inside a read that never completes and never reports the dead worker, so `coal -j N` waited
  for ever (seen as an intermittent CI failure on macOS with Python 3.10; reproduced on Linux with the spawn
  start method, about 1 run in 16). The main process now also watches the workers itself while it waits and ends
  with "a worker process died" (exit 1). Ending the broken pool no longer closes the call queue under the
  management thread, which on Python 3.10 / 3.11 printed a `ValueError: Queue ... is closed` traceback before
  that message (about 1 run in 5 there).
- **Thousands of DNA files:** at most about as many files as the open-file limit allows stay open (least recently
  used closed first, the limit raised towards its hard limit when needed) and their read buffers together take
  about 64 MB at most. 1,000 files failed with "too many open files" in the Python engine with workers, 60,000 open
  files took 8 GB (129 KB each on Python 3.14), and more than 65,535 files overflowed a 16-bit file number: a
  traceback in the Python engine, a silent wrap in coal-cpp.
- **A lone `\r` ends a line anywhere in a DNA file** (universal newlines, as the protein reader): the style used to
  be decided from the first 64 kB, so a file whose line ends change part-way lost every record after the first lone
  `\r` beyond that and failed with a misleading count error. PAL2NAL v14 reads such files.
- **`nohup coal ...` survives the terminal closing:** an already ignored SIGHUP stays ignored (the handler used to
  replace nohup's setting; the run died with exit 129).
- **Out of memory in the main process is one `coal: error: out of memory ...` line** (exit 1, report code
  `out_of_memory`) instead of a traceback (a worker dying of it was already handled).
- **Memory no longer grows with the number of per-codon warnings, in either engine** (finding 9 of the stress
  test, both halves). The run report keeps, per sequence, one byte for its kind, one for a few flags and the
  codon of its first definite mismatch (plus the first flagged and the first excluded gene whole, for the
  summary lines, and where up to 200 failed genes break), instead of every per-sequence entry with every message
  until the end of the run; the entries go to the `--resume` entries file next to the output as the run goes
  (or to a temporary file when a report is wanted without `-o`; nothing with `-q` and no report), and `--report`
  / `--json` read them back in order. A resumed run rebuilds its state from that file one line at a time (the
  checkpoint check no longer loads the whole file either). The Python engine's worker processes hand back each
  sequence's mismatch columns and codons instead of its warning strings (about 2 KB against 35 KB for a
  400-codon gene of unrelated DNA, per gene of the batch in flight in every worker and in the main process), and
  the main process formats a warning only when it prints it; the strings, the report, stderr, the exit status
  and `--resume` are unchanged byte for byte (`tests/golden_memory`: recorded from the previous `main`, both
  engines). Measured on 400-codon genes of unrelated DNA, every codon a warning (peak memory as the summed PSS
  over the process tree; 16 threads): 20,000 genes in default mode, coal-cpp 1,172 MB before and 103 MB after,
  the Python engine 2,822 MB and 45 MB (`-j 1`), 3,226 MB and 338 MB (`-j 16`, 7 workers); 100,000 genes,
  coal-cpp 5,654 MB and 118 MB, the Python engine 13,761 MB and 69 MB (`-j 1`), 14,647 MB and 885 MB (`-j 16`,
  16 workers: 545 MB of that is what the same input in clean form costs with 17 processes, the rest the batches
  in flight). The default mode is also faster on such
  data (100,000 genes, `-j 1`: 213 s to 158 s in Python, 185 s to 161 s in coal-cpp) since no entry store grows;
  the one cost is the main process formatting the strings the workers formatted in parallel before, which only
  shows when stderr takes most of a run. A run with `-o` on such data writes the entries file next to the
  output while it runs (about 55 KB per 400-codon gene of unrelated DNA, removed at the end).
- **App: `python -m coal_app` no longer dies on a non-UTF-8 console** (Gradio's MCP banner emoji on a cp1252
  Windows console or pipe).
- **GUI: a run that stops before pairing says why in the banner** (the engine's own sentence, e.g. "control character:
  sequence files must be plain text (is this a binary or word-processor file?)", file names only); it used to say
  only "COAL stopped: no alignment was written" with the reason inside the collapsed raw-messages box.
- **GUI: the genetic-code tip on a large input says it tested a sample** ("50 of the 50 genes COAL tested (a sample
  of the 100,000)") instead of "50 of 50 genes".
- **GUI accessibility:** the "Methods paragraph" box is announced by its label (its copy button, which sat inside the
  label and took the label for itself, is now a separate "Copy methods paragraph" button), and Gradio's hidden
  duplicate tab bar no longer holds focusable buttons (axe `aria-hidden-focus`).

### Internal

- **Bioconda recipe (packaging/bioconda), second review.** The test line for `pal2nal` without arguments expected
  exit 0 and the usage on stdout and would have failed the first build: it now tests what the command does (exit
  1, usage on stderr, nothing on stdout). `python` is unpinned in `host:` / `run:` and Pythons below 3.10 are
  skipped by a selector (the form for a compiled, non-noarch package); `linux-aarch64` and `osx-arm64` are
  additional platforms; the package is the command line and the Python package only (`build.sh` no longer
  installs the static library, the headers and the CMake package); the two test commands with tab characters are
  quoted YAML strings. Rendered with `conda-render`, `build.sh` and all 23 test commands run against a virtual
  environment; a real `conda build` and `bioconda-utils lint` are still to do (packaging/bioconda/README.md).
  `tests/test_bioconda_recipe.py`.

- **Package metadata follows PEP 639:** `license = "GPL-3.0-only"` (an SPDX expression) with
  `license-files = ["LICENSE", "NOTICE"]` in both `pyproject.toml` files, no licence classifier; the build needs
  setuptools 77.0.3 or later (Python 3.9+, below the floor of 3.10; the Bioconda recipe asks for it too). The
  wheels' metadata says `License-Expression: GPL-3.0-only` and carries both files under `licenses/`; the
  deprecation warnings of every build are gone.

Documentation, tests, continuous integration and packaging:

- **Continuous integration.** A workflow of its own for the R package (`.github/workflows/r-package.yml`): `R
  CMD check --as-cran --no-manual` of `packaging/r/coalr` on Linux (R release and the previous release), macOS
  and Windows, and the check that the vendored C++ sources are `cpp/src`; it runs only when the package, the C++
  core or the sync script change, so the main matrix is not slowed. The nightly run has a ThreadSanitizer job
  (`tests/tools/tsan_threads.py`: `coal-cpp kaks`, `dnds`, `codon-dnds -j 4` with all nine methods and both
  flavours, and `coal-cpp-align -j 4`, each also compared with `-j 1`; 28 commands, no report locally). The C++
  job requires `tests.test_kaks_comeron_li`, `tests.test_short_options` and `tests.test_output_lock` (a missing
  binary fails), checks the manylinux wheel with both binaries in a fresh environment, and installs and
  requires DIAMOND on macOS too: DIAMOND's macOS archive is a universal binary with a native arm64 slice.
- **`coal-cpp`: test hooks only in test builds.** The `COAL_TEST_*` environment variables with which the tests
  pause a run or line up several runs are compiled into `coal-cpp` only with the CMake option
  `-DCOAL_TEST_HOOKS=ON` (default `OFF`): release builds, wheels and packages do not contain them
  (`docs/development.md` "Test hooks").
- **Line ends in a Windows checkout** (`.gitattributes`): the R package (an R source package has LF line ends;
  its vendored C++ sources are compared byte for byte), the Ka/Ks goldens, the real-data cases, the stored
  benchmark results and shell scripts are LF in every checkout, as the fixtures and goldens already were.
  `scripts/sync_r_sources.py --check` reads `cpp/src` with either line end and names a copy that differs only by
  its line ends as that. (With `core.autocrlf`, as on the Windows CI machines, it called all 13 vendored files
  "out of date".)
- **Docs: Linux benchmark numbers.** `docs/benchmarks.md` has the 2026-10-24 run on the same desktop as the Windows
  table but under WSL 2 (up to 1,000,000 sequences, PAL2NAL v14 against both engines, raw data in
  `benchmarks/results/linux_wsl_4159e10/`): COAL's times and memory hold, PAL2NAL's Perl runs about 3.6× faster there
  than on Windows, so the speed-ups are 123× (coal-cpp) and 29× (Python) at 1M on that setup.
- **Docs:** the desktop README's privacy line says none of *your data* leaves the computer and names the one helper
  script Gradio's page requests from a CDN (not needed; the page works with it blocked); the installer name is the
  version-free `COAL-Setup-windows-x64.exe`; `docs/usage.md` says a `kill -9` of a run writing to standard output can
  leave `.coal_<random>.tmp` in the temp folder; `tests/README.md` says which of the "17 local codon alignments" are
  in the repository (12) and which are local-only (5); SKILL.md's "after the engine update" and "not yet verified on
  Windows" sections are current (both engines verified on Windows 2026-10-24); `demo/mcp_demo.md` regenerated
  (`codon_align` keeps going past a broken gene and returns `completed: true` with `exit_code: 2`).
- **Tests:** `test_docs_explainer` finds Playwright's Chromium on every platform (it looked for a Linux-only path and
  was always skipped elsewhere); `test_cpp_vs_reference` copies instead of symlinking where symlinks need a privilege
  (a stock Windows account); `test_resume_with_workers` ignores a venv launcher's same-argv child on Windows; new
  tests for the GUI's download buttons (the one gap found by mutation testing: 12 of 15 planted bugs were caught),
  the stop reason, the sampled tip and the methods label. The PAL2NAL fuzz harness no longer reports a difference
  for `-codontable 10` genes that both programs place by position 1 or 2 nt off the generator's frame (the
  neutralised TGA then straddles two codons of the alignment grid; seeds 20261024 and 20261025 each had one such
  case, with COAL's output equal to PAL2NAL's on the original input). `test_cpp_vs_reference` is skipped on
  Windows with the reason: its shim runs coal-cpp as a child process there, so the reference tests that interrupt
  a run watch the wrong process (before, the test errored on the symlink privilege before getting that far).
- `tests/tools/stress_*.py`: the stress-test generator, harness and helpers; `tests/test_stress_findings.py`.
- **Tests added with the review of 2026-11-13:** the pair-table header is the same in `dnds_cli.cpp`,
  `kaks_cli.cpp`, `dnds_pipeline.py` and coalr (`tests/test_pair_table.py`); a `coal-cpp` built without test hooks
  makes the reference-shim resume tests skip instead of fail; five fixed sleeps in the lock, MSA and pipeline
  tests are event or polling waits (test hook `COAL_TEST_LOCK_BEFORE_TAKEOVER_UNTIL`); the YN / MYN shortcut is
  compared with the full computation on every platform (`--yn-always-pmat`, a testing option).
- **Bioconda recipe** (`packaging/bioconda`): written for 0.6.0 with the version and sha256 marked as
  placeholders for release time; both binaries; GPL-3.0-only with `LICENSE` and `NOTICE`; the new scope in the
  summary; DIAMOND, Clustal Omega and MAFFT under `run_constrained`; tests of `coal dnds`, `coal kaks`,
  `coal pipeline --help` and `coal cite`.
- **Website and README** present the full scope (codon alignment, Ka/Ks, the pipeline, coalr, orthologr) with
  one runnable example each; the website's new sections are marked `data-since="0.6.0"` and stay hidden on the
  deployed site until that release. No document tells users to `pip install coal-align` as if it were on PyPI.

## 0.5.0

- **Run report: a new top-level `engine` field** (`coal.report/1`, next to `source`): `"coal-python"` or
  `"coal-cpp"`, the same strings as the `--resume` checkpoint's `engine`, so a pipeline or the app can tell which
  engine wrote a report (also when `coal` handed a large input to `coal-cpp`). The per-job reports of a
  `coal --batch` report (`coal.batch/1`) carry it too. Apart from `engine` and `runtime_seconds` the two engines'
  reports stay identical (#18).
- **Changed: `summary.not_checked` of a run that stopped at an error now counts only the sequences the run never
  reached**, the number the `coal: summary:` line "N sequences were not checked" gives; the sequences checked
  before the error are counted as what COAL found (`ok`, `warning`, `failed`), so `ok + warning + failed +
  not_checked == total`. It used to count every sequence but the failed one (e.g. `examples/one_broken_gene`
  without `--keep-going`: `ok 0, failed 1, not_checked 19`; now `ok 13, failed 1, not_checked 6`). The per-sequence
  entries are unchanged (`not_checked` with `not_written` before the error, `not_reached` after it). Both engines;
  documented in `docs/report-schema.md` (#18).
- **`pip install coal-align` brings the C++ engine** (issue #16): each release now has one wheel per platform with
  `coal-cpp` inside, at `coal/bin/coal-cpp` in the installed package: Linux x86_64 and aarch64 (`manylinux_2_28`,
  glibc 2.28+, `libstdc++` linked statically), macOS 11+ (one `universal2` wheel for Apple silicon and Intel) and
  Windows x64 (static CRT). pip picks it by itself from PyPI; from a GitHub Release, install the wheel file for your
  platform. The pure `py3-none-any` wheel and the sdist stay as the fallback everywhere else. `coal` finds the bundled
  binary with no settings: the lookup order is now `COAL_CPP`, then `coal/bin` inside the package, then next to
  Python, then `PATH` (docs/cpp-engine.md "Where coal looks for coal-cpp"). New `coal --cpp-path` prints the
  `coal-cpp` that `coal` would use (exit 1 with the reason when there is none or it is another version), for
  pipelines that call the binary directly. The wheels are made by `packaging/wheels/build_platform_wheel.py` (no
  cibuildwheel, no compiled extension) in `release.yml`'s new `cpp-wheels` job and checked there (auditwheel, lipo, a
  fresh-venv install whose output equals the Python engine's); `tests/test_wheels.py` does the same with a local build.
- **Nightly fuzz workflow** (`.github/workflows/nightly.yml`, 03:17 UTC, or by hand with a seed): the fuzz against
  the real PAL2NAL v14 with a fresh seed each night, once per engine (`COAL_FUZZ_ENGINE=cpp` in
  `tests/fuzz_vs_pal2nal.py` runs `coal-cpp`), the new `tests/tools/engine_fuzz.py` (coal-cpp vs the Python engine
  on malformed inputs and every option; Linux and macOS) and the C++ tests plus the generator under
  AddressSanitizer and UndefinedBehaviorSanitizer (`COAL_SKIP_MEMORY_TESTS=1` skips the peak-memory test there).
  Differing cases are uploaded as artifacts; CONTRIBUTING.md "Nightly checks" says how to replay them.
- **coal-cpp: no read past the end of a short gene's DNA** when the expected position of an anchor (after a
  frameshift digit or a mismatch) lies beyond the sequence end; found by the first sanitizer run. Output and exit
  codes are unchanged (the out-of-bounds read could never have produced a match).
- **MCP server on mcp 1.x (the Windows app's) is quiet and names its version:** `python -m coal_app --mcp` /
  `COAL.exe -m coal_app --mcp` with the `mcp` 1.x that Gradio's extra installs reported the mcp library's version
  ("1.10.1") as `serverInfo.version` and logged about 37 KB of INFO lines to standard error per session (into
  Claude Desktop's MCP log). It now reports COAL's version, as on mcp 2.x, and logs warnings and errors only. The
  message for a missing `mcp` package no longer claims that version 2 is required (1.x works too).
- **Packaging:** the source distributions of `coal-align` and `coal-align-app` no longer ship `tests/` (the tests
  need the repository: helpers, fixtures, examples, coal-cpp), and the `coal-align-app` wheel and sdist now carry
  the `LICENSE` (MIT at the time; the app had none).
- **App: a settings file that is not text no longer stops COAL from starting.** `settings.json` in another encoding
  (or binary garbage) gave a traceback at start-up; it now counts as a damaged file: the defaults are used and the
  Settings tab says so, as for a half-written file.
- **App: an unusable results folder is one plain sentence, not a traceback.** With `COAL_RUNS_DIR` pointing where no
  folder can be made (a path under a missing directory, an existing file, no permission), the GUI shows "Could not
  create the results folder <folder>: <reason>" in the banner, the MCP tool `codon_align` raises a tool error with
  that sentence, and `coal-gui` / `python -m coal_app` print it and exit 1 instead of failing in `os.makedirs`.

- **Device inputs are refused** (robustness review): `coal /dev/zero nuc.fa` or `coal - nuc.fa </dev/zero` used to
  copy the device to a temporary file without end, filling the disk silently. Both engines now stop before reading
  with `coal: error: cannot read /dev/zero: it is a device, not a file` (exit 1; `standard input` for `-`). A
  terminal on standard input and `/dev/null` behave as before; pipes, FIFOs and `<(...)` are unchanged.
- **`--report PATH` no longer replaces what is at PATH** (robustness review): it used to rename a temp file over
  anything there, so `--report nuc.fa` destroyed an input, `--report out.fa` the output, and a symlink, FIFO or
  device node (`/dev/null` became a regular file!) was replaced. Both engines now give the report path the rules of
  `-o`: an input or the output file is a usage error before anything runs (`coal: error: the report path r.json is
  an input file` / `is the output file`, exit 2; also for a job line's `--report` and the batch report), symlinks
  are written through, and FIFOs, devices and hard-linked files are written in place. A report into a folder that
  does not exist is a usage error too (`the report path results/r.json: folder /data/results does not exist`, exit
  2) instead of a full run that ends by losing its report, and the report's temp file is removed on SIGTERM /
  SIGHUP / Ctrl-C while it is written, in both engines.
- **A killed worker no longer hangs the Python engine** (robustness review): when the out-of-memory killer (or
  `kill -9`) took one worker of a `-j` run, the main process could block forever in the pool's shutdown. It now
  ends the pool without waiting and stops with `coal: error: a worker process died (out of memory?)` (exit 1, report
  code `worker_died`, no partial output, no temp files). Workers also let SIGTERM end them (the inherited handler
  used to be swallowed inside the work) and **die with the main process** (`PR_SET_PDEATHSIG` on Linux, a parent
  watchdog everywhere), so a `kill -9` of `coal` no longer leaves workers holding its stdout open and a pipeline
  reading it waiting for EOF forever.
- **Python 3.14 on Linux: every multi-worker run failed** with `a worker process died (out of memory?)` after
  half a second. 3.14 made `forkserver` the default start method there, so the workers' parent is the fork server,
  not `coal`, and the watchdog above took that for a dead parent. The workers now wait for EOF on a pipe whose only
  write end `coal` holds, whatever the start method; that also fixes **Windows, where a killed `coal` left its worker
  processes running** (there `getppid()` keeps returning the dead parent's pid). `COAL_MIN_NT_PER_JOB` (default 4,000,000)
  can now be set in the environment, as `COAL_BATCH_NT`, so tests can force worker processes on a small input; CI
  runs the suite on 3.14 (Linux).
- **Worker processes end as soon as the main process fails** (Python engine): after an error on the writing side
  of a `-j` run (a full disk, a file-size limit) the pool was only ended when the cyclic garbage collector got to it
  (the suspended result generator and an exception held in a local kept it alive), so terminated workers and pool
  threads could outlive the error by seconds. The generator is now closed on the way out and the pool ended at once.
- **Checkpoints are checked before `--resume` trusts them** (robustness review): an edited or corrupted
  `<out>.coal-checkpoint.json` (`completed` beyond the input or negative, a string or float where an integer belongs,
  a list at the top level, a missing key, a deleted entries file, a partial file shorter than recorded, a truncated
  or deeply nested file, a folder at the checkpoint's path) made both engines write a truncated output with exit 0,
  or crash with a traceback / raw C++ text. Both engines now check every field's type and range and the partial
  files' sizes, treat any problem as the existing "unreadable checkpoint" notice and start a normal run; a
  checkpoint that cannot be written (`cannot write the checkpoint ...: is a directory`) is a notice and the run
  goes on without one; and `coal-cpp`'s JSON parser has a nesting limit (64) instead of overflowing the stack.
  The `drop` and `bad` column lists are bounded by the run's protein and codon columns too: a crafted
  `"bad": [4000000000000]` made `coal-cpp` fail with a raw `std::bad_alloc`; it is an unreadable checkpoint now, in
  both engines.
- **`--resume` of a `--keep-going` run in fasta or paml format** (Python engine): a checkpoint written after a
  skipped sequence was rejected as "unreadable" (the skipped sequence's index was recorded for clustal / codon
  only, and the checkpoint check requires one index per skipped sequence), so the run started over. The index is
  now recorded in every format, as `coal-cpp` does, and the resumed output is byte-identical to an uninterrupted run.
- **A checkpoint written by the other engine** (robustness review): `coal --resume` now picks the engine that wrote
  the checkpoint (so a `coal-cpp` installed or removed between the run and the resume does not get in the way),
  and when an engine is asked to resume the other one's run it says so instead of blaming the inputs: `out.fa was
  interrupted by the Python engine (coal-python), which must resume it: run again with COAL_ENGINE=python, or start
  a new run without --resume`.
- **`-o /dev/fd/1` and `-o /dev/stdout` work when stdout is a pipe** (robustness review): both engines resolved
  the path to `pipe:[N]` and failed with `cannot write /dev/fd/1: No such file or directory`. An existing target
  that is not a regular file is now written through before any path resolution.
- **`--report /dev/fd/1` and `--report /dev/stdout` work when stdout is a pipe**, in both engines, for the same
  reason (the report failed with `could not write the report /dev/fd/1: No such file or directory` after the
  whole run).
- **`--json | head` is a silent exit 141** (robustness review) in both engines (the Python engine printed a
  traceback and exited 1, `coal-cpp` exited 0), also for `--batch ... --json`. `coal-cpp` now streams the `--json`
  report to standard output one sequence at a time instead of building it whole in memory first (the bytes are
  the same).
- **`coal-cpp` falls back when `TMPDIR` is unusable** (robustness review): a `TMPDIR` naming a missing folder or a
  file made it fail with `filesystem error: temp_directory_path`; it now searches `TMPDIR`, `TEMP`, `TMP`, `/tmp`,
  `/var/tmp`, `/usr/tmp` and the current folder as Python's `tempfile` does, so both engines behave the same.
- **Disk-full errors name the file** (robustness review): the Python engine printed `[Errno 28] No space left on
  device` (and could hang at exit after terminating its workers), `coal-cpp` `could not write the output (disk
  full?)` (and, with the rows going through a full temp folder, exit 0 with an empty output). Both now say
  `cannot write out.fa: no space left on device`, `cannot write a temporary file in /tmp: ...` or `cannot copy
  standard input to a temporary file in /tmp: ...`, exit 1, and leave no temp file behind.
- **Python API checks its arguments** (robustness review): `coal.codon_align` raises `CoalError` with
  `code="duplicate_id"` for a repeated protein or DNA ID (it used to pair them silently; the command line refuses
  the file), and both `codon_align` and `coal.align_files` raise `TypeError` with a message naming the argument for
  a dict, generator, bare string, `None`, a non-`str` ID or a pair of the wrong length, instead of a `KeyError` or
  `AttributeError` from inside the engine. `align_files(pep, "nuc.fa")` with one DNA path as a string or
  `os.PathLike` now means that one file (the string used to be iterated character by character). `nuc_paths` is
  one path or a `list` / `tuple` of paths: a generator (or any other iterator) of paths, which used to work, is a
  `TypeError` now; pass `list(...)`.
- **Smaller robustness fixes** (both engines, identical messages): a NUL or other control byte (except tab, LF,
  VT, FF, CR) anywhere in a protein or DNA file is `<file>:<line>: control character: sequence files must be plain
  text (is this a binary or word-processor file?)` (exit 1, code `input_format`) instead of a confusing row-length
  error, a mangled ID or silently dropped bytes; IDs longer than 200 characters are shortened on stderr to their
  first 100 and last 20 characters (the output and the report keep them whole); `-o -` means standard output; the
  Python engine reports a DNA file that shrank between indexing and reading as `could not read the DNA of record N
  (the file changed after it was read, or cannot be read)`, as `coal-cpp` does, instead of reading another
  record's bytes; an empty DNA file among several is documented as accepted.
- **Mismatching long genes no longer take quadratic time.** When a protein does not fit its DNA exactly, both engines
  locate PAL2NAL's 10-residue chunks in the DNA; every chunk was a scan from the start of the DNA (coal-cpp; the
  Python engine scanned for every chunk that was not where the previous ones put it), so a 35,000-codon gene with one
  wrong codon cost about 10 s per gene with coal-cpp (and the genetic-code hint repeated that under 17 other codes).
  The DNA is now translated once per frame and its 10-residue windows indexed, so a chunk of ordinary residues is a
  lookup; chunks with a frameshift digit, an X or the first Met keep the old scan. Same chunks found, same output and
  messages: 10 genes × 35,000 codons with one wrong codon each, 104 s → 3.0 s (coal-cpp), 27 s → 16 s (Python);
  exactly fitting genes are unchanged. The Python engine's fallback is also about 1.5× faster per short gene.
- **coal-cpp uses far less memory on runs with many warnings or with `--report`.** It kept a JSON tree for every
  sequence with a warning (about 4 KB each: 1 GB at 300,000 mismatching sequences) and, for `--report`/`--json`, built
  every sequence's entry before writing (3 KB each: 3.2 GB at 1,000,000 sequences). It now keeps each non-ok entry as
  its compact JSON line (about 0.6 KB; the Python engine's dict entries take 1.4 KB) and streams the report one
  sequence at a time:
  300,000 mismatching sequences 994 → 266 MB; 1,000,000 clean sequences with `--report` 3.2 GB → 204 MB. The
  Python engine streams its report the same way and no longer builds a view of every sequence for the genetic-code
  hint (1,000,000 sequences with `--report`, one process: 635 → 367 MB). The report JSON is unchanged.
- **Docs:** `docs/usage.md` gives the measured memory per sequence of both engines (about 350 bytes with the Python
  engine after the change above, 190 with coal-cpp; the old figures, 300 bytes and 122/139/145 MB, matched coal-cpp
  only), where each engine puts its temporary files and how much disk the clustal/codon formats need, that
  `-q --progress` still prints progress, and that the checkpoint cadence is best-effort under I/O contention.

- **Less memory with many workers:** COAL's own processes now limit numpy's BLAS thread pool to one thread
  (`OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, unless already set). COAL uses no BLAS, but every
  worker that imported numpy started a pool sized to the machine. On a 20-thread Windows machine each process that
  imported numpy committed about 630 MB (9-11 GB at the default `-j`); the BLAS pool is the likely cause, to be
  confirmed by re-measuring. Programs that use COAL as a library are not affected.
- **Faster start-up** for runs on small files (one COAL call per gene family in a pipeline): about 40 ms less per call
  (104 -> 64 ms on Linux), by not importing `dataclasses`, `subprocess`, `hashlib`, `json` and `shlex` unless needed.

- **`-` reads standard input** (`zcat nuc.fa.gz | coal pep.aln -`), for the protein alignment or one DNA file, in
  `coal`, `coal-cpp` and the `pal2nal` drop-in, also on Windows (no `/dev/stdin` needed). The data is read as bytes
  into a temporary file, as for other pipes, so the output is byte-identical to giving the file's path; messages name
  it `standard input`; a second `-` is a usage error. The run report records `"path": "-"` with the sha256 of the data
  read (the report of a run reading `/dev/stdin` or another pipe now hashes the data too, in both engines, instead of
  the Python engine's hash of an empty re-read).
- **C++ library for embedding** (e.g. in DIAMOND): `cmake --install` puts `coal-cpp`, `libcoal_core.a`, the headers
  (`<coal/core.hpp>`, `<coal/tables.hpp>`, `<coal/seqio.hpp>`, `<coal/version.hpp>`) and a CMake package in place;
  `find_package(coal CONFIG)` or `add_subdirectory` / FetchContent give `coal::core`. As a subproject only the library
  is built (`COAL_BUILD_CLI`, `COAL_BUILD_TOOLS`, `COAL_INSTALL` to add more). Example: `examples/cpp-embed/`; docs:
  docs/cpp-engine.md "Embedding the library".
- **orthologr stand-in takes every PAL2NAL option:** `packaging/orthologr/pal2nal.pl` sent every command line to
  `coal-cpp` when installed, which rejects PAL2NAL-only options such as `-nostderr` and prints `coal`'s help for
  `-h`. It now runs `coal-cpp` only for the options it understands (`-output`, `-codontable`, `-nogap`,
  `-nomismatch`, `-blockonly` and the files) and COAL's `pal2nal` drop-in for anything else; `COAL_PAL2NAL` still
  overrides.
- **orthologr:** `packaging/orthologr/pal2nal.pl` is a stand-in for the `pal2nal.pl` bundled in the orthologr R
  package; copied over it, orthologr's `codon_aln()` / `dNdS()` run COAL (coal-cpp when installed) with their
  unchanged command line and get the same output. Notes on what to expect and on a proper `tool = "coal"` option
  for orthologr: `packaging/orthologr/README.md`.

- **Python API: every error is a `CoalError`** (issue #18). `coal.codon_align` and `coal.align_files` with an
  unsupported `table` now raise `CoalError` with `code="unknown_table"` (the message lists the supported NCBI codes)
  instead of a bare `ValueError`.
- **Python API: `coal.seqio.InputFormatError` is a `CoalError`** with `code="input_format"` (issue #18), so one
  `except CoalError` catches a file that is not FASTA or CLUSTAL too; it stays a `ValueError` for code that catches
  that. The `coal` command's messages, exit code and report are unchanged.
- **`pal2nal | head` no longer prints a traceback** (issue #18): the drop-in's usage text into a closed pipe is a
  silent exit 141, as for `coal` and for the alignment output.
- **Batch mode (`coal --batch JOBS`, #15):** thousands of small conversions in one process, one `coal` command line
  per line of JOBS (`-` = stdin), the command line's options as the defaults of every job, `--batch-out-dir DIR` for
  jobs without `-o`, a `coal.batch/1` report with `--report`/`--json`, per-job stderr prefixes
  (`coal: [k/n PEP] warning: ...`) and a `coal: batch:` summary line. A failing job does not stop the batch (exit 1
  at the end; 2 if jobs only left sequences out; mistakes in JOBS are reported with their line number before
  anything runs). Every output is byte-identical to the single run; about 40× less wall time on tiny jobs. Python
  only: `coal-cpp --batch` is a clear error and batches are never dispatched to it. `--resume` and cross-job
  parallelism are not supported in a batch yet.
- **The app no longer waits forever for a stopped engine:** when an engine run hits the app's time limit (or the
  demo pipeline's), it is sent SIGTERM and, after the grace period, SIGKILL; the kill now takes the engine's whole
  process group, so its worker processes die with it. Before, a run that outlived the grace period on a slow machine
  was killed alone, its orphaned workers kept its output pipes open and the app waited for them indefinitely (the
  CI app suite hung until the job timeout). Windows already killed the tree. Test commands that start the engine
  now give it no stdin and always carry a timeout.
- **`coal-cpp` mirrors Python 3.13's argparse** (it mirrored 3.10's): `-qx` is "unrecognized arguments: -x" rather
  than an ignored explicit argument, `-hx` prints help, `-cod=2` is accepted, invalid choices are quoted
  (`invalid choice: '7' (choose from '1', '2', ...)`), and only the `--` marker is dropped from the positionals (a
  second `--` is a file name). `tests/test_cpp_cli.py` compares the two command lines live on an interpreter
  whose argparse behaves that way (probed, not keyed on the version) and skips elsewhere with the differences named.
- **Documented for 1.0:** `docs/python-api.md` (the public Python API: `coal.codon_align`, the new re-export
  `coal.align_files` = the `coal` command as a function, `coal.CoalError` with its `code`/`seq_id`/`index`, and the
  semantic-versioning promise for command-line flags, output bytes, exit codes, the report schema and these
  functions), `docs/report-schema.md` with a JSON Schema (`docs/report.schema.json`, draft 2020-12) for the run report
  `coal.report/1` (every field, status, pairing mode, message and hint code; additive changes only within `/1`),
  `CONTRIBUTING.md`, GitHub issue forms and a pull-request template. New tests run the documented examples
  (`tests/test_python_api.py`) and validate freshly generated reports from both engines against the schema
  (`tests/test_report_schema.py`); `tests/report_example.json` is now a real `examples/one_broken_gene --keep-going`
  report. No change to COAL's behaviour or output.
- **orthologr:** `packaging/orthologr/pal2nal.pl` is a stand-in for the `pal2nal.pl` bundled in the orthologr R
  package; copied over it, orthologr's `codon_aln()` / `dNdS()` run COAL (coal-cpp when installed) with their
  unchanged command line and get the same output. Notes on what to expect and on a proper `tool = "coal"` option
  for orthologr: `packaging/orthologr/README.md`.
- **Repository tidied for outside users:** a short README with a section on using COAL in pipelines; the detailed
  documentation is in `docs/`. Hackathon logs and pitch material are removed (kept on the `hackathon-archive` tag),
  review-round markers are gone from comments, and engine test modules are named after what they test. No change to
  COAL's behaviour or output.
- **`coal-cpp` reports the same error as `coal` when a protein file has two defects.** A FASTA record is now
  handed back before the `>` line that ends it is checked, so a short row is reported before a later duplicate or
  empty header; and a CLUSTAL file's rows are parsed before its `#` mark line is read, so a malformed CLUSTAL line is
  reported before a later non-ASCII character. Found by differential fuzzing; `tests/test_cpp_error_precedence.py`.
- **Package renamed to `coal-align`** (`coal` on PyPI is another project). The `coal` command, the `coal` Python
  module and `pal2nal` are unchanged; release files are now `coal_align-<version>-py3-none-any.whl` and
  `coal_align-<version>.tar.gz` (`coal_src.tar.gz` stays).
- **The app is its own package, `coal-align-app`:** the web GUI, the MCP server, the viewer, the Windows app and the
  Hugging Face Space moved to `apps/coal-app/` (package `coal_app` unchanged; `coal-gui` and `python -m coal_app` work
  as before). `coal-align` is now the engine alone (`coal`, `pal2nal`, the examples) with no dependencies and no app
  code, so pipelines (OrthologR, DIAMOND workflows) install only that. The `gui` and `mcp` extras are on the app:
  `pip install "coal-align @ .../coal_src.tar.gz" "coal-align-app[gui] @ .../coal_app_src.tar.gz"`, or from a checkout
  `pip install -e . -e "apps/coal-app[gui]"`. Both packages always have the same version, and the app depends on that
  exact `coal-align`. Releases add `coal_align_app-<version>-py3-none-any.whl`, `coal_align_app-<version>.tar.gz`
  and the version-free `coal_app_src.tar.gz`. The Windows launchers are now `apps\coal-app\desktop\run_from_source.bat`
  and `apps\coal-app\desktop\setup_mcp.bat`. The app's tests (`apps/coal-app/tests/`) are grouped by feature; the
  checks of the docs, the website and CI moved to `tests/` (`test_docs_*.py`, `test_repo_ci.py`).

- **Bioconda recipe prepared** (`packaging/bioconda/`): `meta.yaml`, `build.sh` (coal-cpp via CMake into `bin/`, then the
  Python package), a fixture pair and a test that `coal` finds the installed `coal-cpp` and hands a 1 MB input to it.
  Not submitted yet; `packaging/bioconda/README.md` explains the sha256 step and the pull request to `bioconda-recipes`.

## 0.4.2

- **COAL in Claude Desktop, from the installer:** the Windows installer's last page has **Connect COAL to Claude
  Desktop** (ticked). Claude Desktop then starts COAL itself (`COAL.exe -m coal_app --mcp`) and can run it on your
  files; uninstalling removes the entry. `codon_align` can save the alignment where you ask (`save_as`), e.g. next to
  the input files.
- **C++ engine included:** the Windows app ships `coal-cpp.exe`, and `coal` hands inputs of 1 MB or more to it when it
  is found and is the same version (same output, faster on large inputs); `COAL_ENGINE=python|cpp` chooses.
- **Quieter, clearer engine output:** a wrong genetic code is reported once, before at most 20 per-codon warnings
  (`--max-warnings`, 0 = all; the `pal2nal` drop-in still prints them all); `--help` is grouped and short.
- **MCP results tell the truth:** `codon_align` returns `completed` (COAL finished), `all_ok` / `needs_attention`
  (whether every sequence matched its protein), `counts`, `suggested_codontable` and `next_step`, and a short summary
  of COAL's messages instead of raw stderr. The old `ok` field (which meant "finished") is gone. A missing file, an
  unknown genetic code or an unsupported format is now a tool error (`isError`), and the arguments have descriptions
  and allowed values. The server reports its version.
- **Windows, from source:** `desktop\run_from_source.bat` (the window) and `desktop\setup_mcp.bat` (the MCP server
  for Claude, printing ready-to-paste registrations) use separate environments, pick Python 3.12 when available,
  check the Python version and repair half-finished installs.
- **Windows app:** the MCP endpoint is off by default (it was on a random port no client could be registered with);
  `COAL.exe --mcp` turns it on.

## 0.4.1

The first release with a working Windows installer (the 0.4.0 build failed, so 0.4.0 has no downloads).

- **Clearer app:** a short header, section headings (Files, Genetic code, Output, Results), the genetic code first,
  a "Start alignment" button, and no empty results area before a run. When the genes fit another genetic code, the
  result says "Wrong genetic code?" with a one-click re-run instead of a page of warnings. Long lists of genes with the
  same problem are folded. One settings place (Gradio's own settings link is hidden). Amber accent colour.
- **Large inputs:** the alignment view is an honest preview at scale (problems first, overview strips, centre on any
  column), and a Settings tab remembers your defaults.
- **Better messages:** errors say what to do next (e.g. "Add --keep-going to write the other sequences"), swapped or
  unaligned input files are recognised, missing files are named plainly, and frameshift locations are only given as
  "at column N" when the evidence confirms it.
- **Windows app fixes:** the packaged app no longer crashes on start-up (UTF-8 output), and file paths with any
  characters work.
- **Examples** set their own genetic code, so one example can't leave the wrong code behind for the next.
- **Tested on real data:** 17 published codon alignments are rebuilt exactly; 12 ship as tests.

## 0.4.0

COAL turns a protein alignment plus the unaligned coding DNA into a codon alignment, like PAL2NAL v14, with the same
output apart from a short list of documented, deliberate differences (see the README).

- **Windows app:** an installer (or portable zip) with a local GUI: a plain-language health check for every gene, an
  alignment view that marks where a gene breaks, one-click re-run with a suggested genetic code, and a ready-to-paste
  methods sentence. Your data stays on your computer.
- **All four PAL2NAL output formats** (`-output fasta|clustal|paml|codon`), `-blockonly`, and a drop-in `pal2nal`
  command that accepts an unchanged PAL2NAL command line.
- **Run report** (`--report FILE.json` / `--json`) with each gene's status, the genetic-code suggestion and the methods
  sentence; `--keep-going` writes every other gene when one can't be aligned (exit 2); `--strict` refuses
  probably-wrong data; `--max-warnings`, `--progress` and `--resume` for very large runs.
- **Safer pairing:** IDs that match only after removing a version suffix, case or `sp|…|` fields are paired by ID,
  and COAL says so.
- **AI agents:** an MCP server with two tools, `codon_align` and `explain_report`.
- **C++ engine** (`coal-cpp`, built from source with CMake): same options and byte-identical output, 3-6x faster
  than the Python engine in our measurements.
- **Tested:** 285 stored checks against PAL2NAL v14's own output and 20,000 random comparisons with the real PAL2NAL
  v14, with 0 unexplained differences.

## 0.3.0

- **Large inputs:** a streaming pipeline with constant memory, automatic multi-process `-j` (only when each worker
  gets enough work) and an optional numpy fast path (`".[fast]"`); 500,000 sequences: 2.5 GB -> 143 MB of memory.
- All four PAL2NAL output formats and `-blockonly`, `--strict`, the drop-in `pal2nal` command, the GUI, MCP tools
  and the Windows app (see README.md).

## 0.2.0

- NCBI genetic codes 1-23 (`-codontable`), extra sequence around the gene (UTRs, stop codon) skipped, PAL2NAL's
  frameshift digits, `-nogap` and `-nomismatch`, and PAL2NAL's handling of codons that don't match.
- About twice as fast; on Windows 29x faster than PAL2NAL on 10,000 sequences and 119x on 10,000-codon genes,
  with identical output.

## 0.1.0

- First version: a pure-Python, dependency-free `coal` command that builds a codon alignment from a protein
  alignment and coding DNA, checked against stored PAL2NAL v14 outputs.
