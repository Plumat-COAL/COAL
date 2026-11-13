# Protein alignment for the codon step

COAL turns a protein alignment plus the coding DNA into a codon alignment. This page covers the step before it,
which makes the protein alignment:

| Input | Aligner | Where |
|---|---|---|
| ortholog **pairs** (e.g. DIAMOND reciprocal best hits) | global Needleman–Wunsch/Gotoh, BLOSUM62, affine gaps, full length | C++ `coal::align_pair` / `align_pairs` (`<coal/align.hpp>`), the `coal-cpp-align` tool, Python `coal.align` |
| **clusters** of 3 or more (e.g. DIAMOND DeepClust) | Clustal Omega (default) or MAFFT, run as external programs | Python `coal.msa` |
| DIAMOND hits **without realigning** | the hit's own gapped alignment, padded to full length | C++ `coal::expand_local_hit`, `coal-cpp-align --diamond`, Python `coal.align.expand_local_hit` |

These parts are new and not yet covered by the 1.0 promise of [python-api.md](python-api.md) (`import coal` does
not export them).

## Pairwise alignment

### Scoring and defaults

- A gap of length *k* costs `gap_open + (k - 1) * gap_extend`. This is the EMBOSS `needle` and Biopython
  `PairwiseAligner` convention: `gap_open` is the cost of the first gap position.
- **Default: BLOSUM62, open 10, extend 0.5, end gaps free.** These are EMBOSS needle's defaults. "End gaps free"
  means that gaps before the first or after the last residue of either sequence cost nothing (needle with
  `-endweight` off).
- **`pwalign` preset: BLOSUM62, open 14, extend 4, end gaps penalised.** This is what
  `pwalign::pairwiseAlignment(type = "global", substitutionMatrix = "BLOSUM62")` scores with its default
  `gapOpening = 10, gapExtension = 4`. pwalign charges `gapOpening + k * gapExtension`, which is open 14, extend 4 in
  the convention above, and its global mode penalises end gaps. Use this preset to match the orthologr fork when it
  switches pwalign to BLOSUM62.
- Matrices: BLOSUM62, BLOSUM45, BLOSUM80 and PAM250. These are the NCBI matrices, the same values as Biopython's
  files and EMBOSS's `EBLOSUM62`.
- Gap costs: 1000 ≥ open ≥ extend ≥ 0, in multiples of 0.001. Scores are computed exactly as integers, scaled by the
  smallest factor that makes the gap costs whole numbers.
- **Full length.** Every residue of both sequences is in the alignment. Overhangs are gapped, never cut off.
  orthologr cuts off the overhanging ends of the longer protein.
- Residues: letters in either case, and `*`. Letters missing from the matrix (J, O, U) score as X. Any other
  character (`-`, `.`, digits) is an error for that pair, so give ungapped sequences.

Why free end gaps are the default: with penalised end gaps, the overhang of a longer protein (an extra exon, an
N-terminal extension) is paid for as one long gap, and the aligner tends to pull terminal residues into wrong
pairings to save it. The two settings gave dN/dS with r = 0.992 on the 1,000-gene set (see the benchmark).

### Algorithm

This is Gotoh's three-state recurrence written so that most of the work runs without a dependency between
neighbouring cells:

- **Vectorised row passes.** X (gap in b) and M/D depend only on the previous row, so they run as one vectorised
  pass.
- **Max-plus scan for Y.** Y (gap in a) unrolls to a max-plus prefix scan of D: `Y[j] = max_k(D[k] - (j-1-k)·extend)
  - open`. With AVX2 this runs 8 lanes at a time as a log-step scan; elsewhere it is a 1-operation-per-cell loop.
- **AVX2 dispatch.** With GCC or Clang on x86 the row loops are compiled twice, for the baseline (SSE2) and for
  AVX2, and the CPU picks one at run time. MSVC and ARM builds get the baseline, which the compiler vectorises for
  its own target.
- **Traceback** takes 4 bits per cell, so memory is `n·m/2` bytes. Titin against titin (35k × 35k) needs ~600 MB.
  Pairs above `max_cells` (default 4·10⁹ cells, ~2 GB) are refused for that pair only.
- **Threads.** `align_pairs` runs on the repository's thread pool (`pool.hpp`) and starts the biggest pairs first.

**parasail** (BSD-style licence, compatible with COAL's GPL-3.0-only) was measured, not vendored. Its 16-bit `sg_trace_scan`
kernel is about 1.5–2× faster per pair. Using it would mean dropping a large generated-code dependency into a
standard-library-only project, plus 16-bit overflow fallbacks and integer-only gap costs. COAL's kernel aligns
about 1,500 pairs a second on one thread, so alignment is no longer the slow step.

### C++ API (`#include <coal/align.hpp>`, namespace `coal`, linked with `coal::core`)

```cpp
struct AlignParams { std::string matrix = "BLOSUM62"; double gap_open = 10.0, gap_extend = 0.5;
                     bool free_end_gaps = true; uint64_t max_cells = 4000000000ULL;
                     static AlignParams needle(); static AlignParams pwalign(); };
struct PairAlignment { std::string a, b; double score; size_t length, identities, similarities, gaps;
                       std::string error; };
PairAlignment align_pair(const std::string& a, const std::string& b, const AlignParams& p = {});
std::vector<PairAlignment> align_pairs(const std::vector<std::pair<std::string, std::string>>& pairs,
                                       const AlignParams& p = {}, int jobs = 0);   // errors per pair in .error
double score_alignment(const std::string& row_a, const std::string& row_b, const AlignParams& p = {});
void check_align_params(const AlignParams& p);
const SubstitutionMatrix* substitution_matrix(const std::string& name);
struct LocalHit { long long qstart, qend, sstart, send; std::string qseq_gapped, sseq_gapped; };
enum class Flanks { Gap, Align };
PairAlignment expand_local_hit(const std::string& query, const std::string& subject, const LocalHit& hit,
                               Flanks flanks = Flanks::Gap, const AlignParams& p = {});
```

Errors are `coal::CoalError` with these codes: `align_input` (a bad residue), `align_params`, `align_too_large`
and `local_hit`. As with the rest of `coal_core` there is no global state and no output, and the functions are
thread-safe. `coal::core` now links `Threads::Threads`, and the CMake package finds it.

### `coal-cpp-align`

```
coal-cpp-align [options] A.fasta [B.fasta]
  pairs: n-th of A with n-th of B; one file: records 1+2, 3+4, ...; --pairs TSV: ID pairs
  --diamond TSV [--flanks gap|align]     expand DIAMOND hits instead of aligning
  -m BLOSUM62|BLOSUM45|BLOSUM80|PAM250  --gap-open 10  --gap-extend 0.5  --end-gaps free|penalized
  --preset needle|pwalign  --max-cells N  -j N  -o out.fasta  --scores scores.tsv
```

- **Output:** the two gapped rows of every pair, one after the other, in input order.
- **`--scores`:** a TSV with the columns `pair id_a id_b len_a len_b length score identities similarities gaps
  status`.
- **Exit status:** 0 when every pair aligned, 2 when some pairs failed (each is named on stderr and in the TSV, and
  the others are written), 1 for unusable input or options.
- The tool is built and installed with `coal-cpp` (`COAL_BUILD_CLI`). The platform wheels do not ship it yet (see
  Known gaps).

### Python (`coal.align`)

```python
from coal import align as A
A.AlignParams(matrix="BLOSUM62", gap_open=10.0, gap_extend=0.5, free_end_gaps=True)   # .needle(), .pwalign()
A.align_pairs(pairs, params=A.AlignParams(), jobs=0, binary=None, python_fallback=True)  # -> [PairAlignment]
A.align_pair(a, b, params)              # pure Python, the same kernel and tie-breaking (slow)
A.score_alignment(row_a, row_b, params)
A.expand_local_hit(query, subject, hit, flanks="gap", params)
A.read_diamond_hits(path, first_only=True)   # DIAMOND_COLUMNS first
A.find_align_binary()                  # $COAL_CPP_ALIGN, coal/bin, next to Python / coal-cpp, PATH
```

`pairs` are `(seq_a, seq_b)` or `(id_a, seq_a, id_b, seq_b)` tuples. `align_pairs` runs `coal-cpp-align` when it
finds it. Otherwise it uses the pure-Python kernel, which gives byte-identical results at roughly 1/1000 of the
speed.

### Validation

`tests/test_align.py` checks the following:

- **Stored optima.** The optimal score equals Biopython 1.88 `PairwiseAligner` and EMBOSS needle 6.6.0 exactly. The
  check runs on 48 fixed pairs × 7 scoring settings (`tests/fixtures/align/reference_scores.json`, made by
  `tests/tools/make_align_fixture.py`), in CI.
  - The pairs cover related, unrelated and extended sequences, 1–4 residues, lower case, B/Z/X and `*`.
  - The settings are the needle default, needle with penalised ends, the pwalign preset, BLOSUM45 12/2,
    BLOSUM80 9/1, PAM250 11/1 and DIAMOND's 11/1.
- **Live Biopython comparison.** When Biopython is installed, 120 new random pairs × 4 settings are compared. During
  development: 2 × 300 random pairs × 8 settings against Biopython and 150 × 2 against needle, with 0 score
  differences.
- **Alignments.** Every alignment returned scores exactly the optimum under an independent rescoring. Where it
  differs from needle's, it is equally scored: needle chose the same alignment in 143/150 pairs with free end gaps
  and 137/150 with penalised end gaps.
- **Full length.** Each row with its gaps removed is the input sequence, and no column is gaps in both rows.
- **Python and C++ agree.** The pure-Python kernel and the C++ tool return byte-identical alignments and scores.
  `-j 1/2/0` give identical output, as do GCC, Clang and the baseline/AVX2 builds.
- **Edge cases.** Errors are isolated per pair, and `--max-cells`, the pairing modes and the options are covered.
  The tests also run under ASan/UBSan (the nightly sanitizer job).
- **needle quirk.** With free end gaps, needle reports a score of 0 when one sequence is a single residue (6 of the
  fixture's comparisons). Biopython and COAL report the optimum, for example 5 for K against a sequence containing
  K.

## The DIAMOND route (no realignment)

DIAMOND 2.2 can write its alignment:
`--outfmt 6 qseqid sseqid qstart qend sstart send qseq_gapped sseq_gapped`. `expand_local_hit` turns one hit into a
full-length pairwise alignment that COAL's codon step accepts:

- **The block.** Only the gap pattern of `qseq_gapped`/`sseq_gapped` is used, and the residues come from the full
  sequences. The hit must agree with them, except that `X` in the hit matches anything, because DIAMOND prints
  masked residues as X. The first line of each query–subject pair is used, which is DIAMOND's best HSP.
- **Flanks.** With `flanks="gap"` (the default, no realignment at all), the unaligned ends face gaps, staggered: the
  query's left flank, then the subject's left flank, then the block, then the right flanks. With `flanks="align"`
  only the flanks are aligned globally, which is cheap.

**Local-vs-global caveats:**

1. **The hit can cover only part of a protein**, such as one domain, or a gene fusion or split. With
   `flanks="gap"` the residues outside it never pair, so those codons drop out of dN/dS. With `flanks="align"` they
   are aligned, but independently of the block.
2. **Only one HSP is used.** A second HSP of the same pair, such as the other half around a large indel, is ignored.
3. **Different scoring.** DIAMOND scores BLOSUM62 with open 11, extend 1. The default NW uses 10 and 0.5.
4. **Masking.** DIAMOND masks low-complexity regions by default, and inside a masked stretch its alignment is
   arbitrary. Pass `--masking 0` when its alignments are going to be used.
5. **Heuristic alignment.** DIAMOND's banded, X-drop extension is not guaranteed to be the optimal alignment even
   where it covers the whole protein.

## Clusters: Clustal Omega and MAFFT (`coal.msa`)

```python
from coal import msa
msa.find_clustalo()          # $COAL_CLUSTALO, `clustalo` on PATH, else COAL's cache
msa.install_clustalo(method="auto")   # -> path; "conda" | "micromamba" (Bioconda; Linux and macOS)
msa.run_msa(records, tool="clustalo", threads=1, exe=None, extra_args=(), timeout=None)  # -> [(id, row)], input order
msa.align_clusters(clusters, tool="clustalo", jobs=0, threads=1, pairs_with="nw")         # -> {name: MsaResult}
msa.find_mafft(); msa.install_mafft(); msa.cache_dir(); msa.tool_version(exe)
```

**Install.** `install_clustalo()` works like orthologr's `diamond_install()`. It installs into COAL's cache folder:
`$COAL_CACHE_DIR`, or `~/.cache/coal`, `~/Library/Caches/coal`, `%LOCALAPPDATA%\coal\Cache`.

- On Linux and macOS it uses the conda, mamba or micromamba already installed, or else micromamba fetched from its
  GitHub releases. It creates `envs/clustalo` with `clustalo=1.2.4` from Bioconda.
- That is the only route. The pre-compiled binaries of the Clustal site, which earlier versions of this code
  could fetch (`method="official"`), no longer exist: checked from a browser on 2026-11-12, the site's address
  serves an unrelated "people search" site and the three binary paths (1.2.4 for Linux, 1.2.3 for macOS, 1.2.2
  for Windows) answer with "not found" pages. There is nothing to pin, and whatever appears there later would not
  come from the Clustal authors, so COAL holds no address of that site. `method="official"` is refused with this
  explanation, and `sha256=` / `allow_unverified=True` (`--sha256`, `--allow-unverified`) do not bring it back:
  with Clustal Omega they are an error.
- **Windows: there is no automatic Clustal Omega install.** Bioconda has no Windows package and the only
  pre-compiled Windows build was the one of that site. `install_clustalo()`, `coal pipeline install clustalo` and
  `--install-tools` say so and name what works:
  1. run COAL under WSL, where `coal pipeline install clustalo` installs it from Bioconda;
  2. use a `clustalo.exe` you obtained and trust: put it on `PATH` or set `COAL_CLUSTALO`;
  3. use MAFFT (`--aligner mafft`) with a MAFFT you installed yourself (`COAL_MAFFT`). COAL has no verified
     MAFFT download for Windows either (Bioconda has no Windows package; no hash is pinned for the Windows
     package of mafft.cbrc.jp, so it is not fetched), and MAFFT on Windows is not tested with COAL.

  Pair mode needs neither tool (COAL's own aligner), and clusters of two members do not either.

**Verification: nothing that was downloaded runs before its SHA-256 was checked.**

| Download | How it is verified |
|---|---|
| micromamba | one fixed release (`msa.MICROMAMBA_VERSION`, 2.10.0-0; not "latest") over HTTPS, kept only if its SHA-256 is the one pinned in `src/coal/msa.py` (four assets; each equals the `.sha256` file of that release) |
| Clustal Omega, MAFFT, DIAMOND through conda / mamba / micromamba | conda's own: HTTPS, each package checked against the channel index |
| DIAMOND release archive (`coal pipeline install diamond`) | HTTPS, kept and unpacked only if its SHA-256 is the one pinned in `src/coal/dnds_pipeline.py` for DIAMOND 2.2.8 (Linux, macOS, Windows). DIAMOND publishes no checksums: the pinned values are those of the files downloaded on 2026-11-12. Another version needs `sha256=` or `allow_unverified=True` |
| Clustal Omega from the Clustal site ("official") | **removed.** The files no longer exist and the address serves an unrelated site (2026-11-12), so COAL never downloads from it, with or without a hash |

A hash that does not match never falls back to "run it anyway".

**What is in the cache is checked before it is used.** When an install has been verified, COAL writes
`<binary>.coal-manifest.json` next to the binary (for DIAMOND on Windows also its DLLs) with the SHA-256 of the
installed files, where they came from and how they were verified. `find_clustalo()`, `find_mafft()`,
`find_diamond()`, the installers and the pipeline use a binary from the cache only while it still has that hash.
A cache binary **without a manifest is not used**: that is every tool installed by COAL 0.5.0 or earlier (whose
Clustal Omega download was not verified at all) and anything that was not installed by COAL. COAL says so once,
names the file, and `coal pipeline install <tool> --force` (or `install_clustalo(force=True)` ...) installs it
again. Tools on `PATH` or named by `COAL_CLUSTALO` / `COAL_MAFFT` / `COAL_DIAMOND` / `--diamond` are yours and are
not checked. The manifest lies in the same folder as the binary, so this detects old, incomplete or altered
installs; it is not a defence against someone who can write to your cache folder. For a conda-installed tool the
manifest covers the program that is started (`clustalo`, `mafft`, `diamond`), not the rest of the environment.

Downloads go to a new file with a random name in the destination folder (never to a fixed `<name>.part`), the
hash is taken from that file as it is on disk, and binaries are put in place by renaming. One install of a tool
runs at a time (`<cache>/tools/<tool>.install.coal-lock`): a second one waits and then finds the tool installed.
A zip archive (DIAMOND for Windows) is unpacked flat, and a member whose name has a `:` (drive letter, alternate
data stream), is absolute, has a `..` component or is a reserved device name is refused with the whole archive.

- Tested here: the micromamba route end to end with the pinned release (Clustal Omega 1.2.4 installed in 30 s),
  and the pinned DIAMOND 2.2.8 archive for Linux.

`install_mafft()` installs MAFFT from Bioconda the same way (Linux and macOS). On Windows it installs nothing
(see above): install MAFFT's own Windows package yourself and set `COAL_MAFFT`.

**Running.**

- **Clustal Omega** runs as `clustalo -i in -o out --outfmt=fa --seqtype=Protein --output-order=input-order
  --threads=N`.
- **MAFFT** runs as `mafft --auto --amino --anysymbol --thread N`.
- **Sequence names.** Sequences are passed under positional names, so any ID works, duplicates included. The result
  is checked: the same residues, rows of equal length, and input order.
- **Stops.** Clustal Omega would align `*` as a residue, so a trailing `*` is removed and comes back as a last column
  of its own, and an internal `*` is aligned as X and then put back.
- **Long sequences.** Neither tool has a length limit. Memory grows with length, and Clustal Omega's HMM stage is
  the slow part.

**Batching** (`align_clusters`):

- One aligner process per cluster, `jobs` at a time (default: one per CPU), each with `threads` threads, biggest
  clusters first.
- Many small clusters are fastest with `jobs = CPUs, threads = 1`. Clustal Omega's threading only helps large
  clusters, so use few jobs × many threads for a few big ones.
- One failing cluster sets that cluster's `.error`, and the others are unaffected.
- Single sequences are returned as they are.
- With `pairs_with="nw"` (the default), 2-sequence clusters go to the pairwise aligner instead of an MSA process.
  This is about 60× faster for pairs (benchmark below) and gives results closer to the NW reference.

### Licences: why Clustal Omega is not bundled

COAL is GPL-3.0-only (from the release after 0.5.0; `LICENSE`). **Clustal Omega is GPL-2.0-or-later** (Bioconda's
metadata and the notice in its source files; its `COPYING` is the GPL v2), linked to argtable2 (LGPL-2.0). DIAMOND
is GPL-3.0-or-later. **MAFFT is BSD-3-Clause.**

**Bundling Clustal Omega would now be licence-compatible.** "Version 2 or later" code can be distributed under
version 3, so a COAL wheel or app with the Clustal Omega binary inside would be GPL-3.0-only as a whole, with no
second licence to declare (COAL 0.5.0 and earlier are MIT: there such a wheel would have been
`MIT AND GPL-2.0-or-later`). COAL only *runs* Clustal Omega as a separate program and exchanges files with it, so
nothing else changes.

It is still **not bundled**, for practical reasons:

- distributing the binary means shipping its licence text and offering the exact corresponding source (Clustal
  Omega, argtable2, and the statically linked libstdc++/libgomp) and keeping that offer for three years;
- no pre-compiled build by its authors can be downloaded any more, for any system (see "Install" above), so we
  would have to build and maintain the binaries ourselves, with the same obligations.

**Windows is where this matters.** With the Clustal site gone there is no Clustal Omega for Windows that COAL
could fetch and verify, so bundling a `clustalo.exe` that we build from the Clustal Omega source (with the
licence text and the source offer above) is now the realistic way to make cluster mode work on Windows without
WSL. That is not done here.

Instead:

- auto-install into the user's cache (above), so that Clustal Omega is available after one call;
- in the Bioconda recipe for COAL, list `clustalo` as a run dependency, which "always there" can rely on and which
  is not distributing it inside COAL;
- in the R package, call `install_clustalo()` from R in the same way as `diamond_install()`.

MAFFT could also be bundled under its BSD licence. It is a set of shell scripts and binaries (a Cygwin bundle on
Windows) and is slow per small cluster (below), so it stays optional. This is a technical reading, not legal advice.

## Benchmark

Data: orthologr's 1,000-gene *A. thaliana* / *A. lyrata* example (`inst/seqs/ortho_{thal,lyra}_cds_1000.fasta`).

- **Pairs.** The CDS were translated, and DIAMOND 2.2.8 blastp in both directions gave 999 reciprocal best hits:
  313.8 M DP cells, median length 384.
- **Coverage.** The DIAMOND hit covers ≥ 90% of the shorter protein in 964 pairs and < 50% in 6.
- **Machine.** A 2-core cloud VM shared with other jobs, so CPU times are the steadier numbers and thread scaling
  could not be measured here.
- Full table and notes: [benchmarks/results/align_layer_linux_20261111/results.md](../benchmarks/results/align_layer_linux_20261111/results.md).
  Reproduce with `python benchmarks/align_bench.py --cds-a ... --cds-b ... --diamond diamond --clustalo clustalo
  --mafft mafft --out bench/`.

**Time to align the 999 pairs** (1 thread unless stated):

| Method | CPU (s) | Pairs/s | GCUPS |
|---|---:|---:|---:|
| `coal-cpp-align` NW, needle defaults (GCC 13, AVX2 picked at run time) | 0.43 | ~1,550 | 0.74 |
| `coal-cpp-align` NW, Clang 18 build (AVX2) | 0.32 | ~2,100 | 0.99 |
| `coal-cpp-align` NW, pwalign preset | 0.47 | ~1,470 | 0.66 |
| `coal-cpp-align` NW, baseline x86-64 (SSE2: MSVC, no AVX2) | 0.66 | ~1,000 | 0.48 |
| `coal-cpp-align` NW, no SIMD at all | 1.44 | ~460 | 0.22 |
| Biopython PairwiseAligner (C), same scoring | 3.97 | ~170 | 0.08 |
| parasail 2.6.1 `sg_trace_scan_16` (16-bit, integer gaps) | 0.21 | ~3,100 | 1.48 |
| Clustal Omega 1.2.4, one process per pair (2 jobs) | 30.9 | 26 | – |
| MAFFT 7.526, one process per pair (2 jobs) | 205.6 | 6 | – |
| DIAMOND hit expanded, flanks = gap / align | 0.02 / 0.03 | ~25,000 | – |

For scale, orthologr 0.4.3 spent 377 s on `pairwise_aln` for 1,000 pairs, or 24 s with BLOSUM62 (29–30 Oct tests,
on another machine). DIAMOND blastp for the RBH search took 2.0 s both ways, and 2.2 s with the gapped columns: the
columns cost about nothing.

**dN/dS** (NG86 from Biopython's `codonalign`, used only to compare the alignments; codon alignments by `coal-cpp`;
codons with a gap or stop dropped), against the needle-default NW:

| Protein alignment | r dN | r dS | r dN/dS | pairs \|Δω\| > 0.05 |
|---|---:|---:|---:|---:|
| NW pwalign preset (open 14 / extend 4, ends penalised) | 0.954 | 0.943 | 0.992 | 33 |
| NW open 11 / extend 1 (DIAMOND's gaps) | 0.991 | 0.993 | 0.998 | 13 |
| Clustal Omega, pair as a cluster | 0.977 | 0.987 | 0.988 | 77 |
| MAFFT, pair as a cluster | 0.990 | 0.993 | 0.996 | 23 |
| DIAMOND hit, flanks = gap | 0.884 | 0.891 | 0.936 | 99 |
| DIAMOND hit, flanks = align | 0.932 | 0.934 | 0.972 | 81 |
| DIAMOND `--masking 0`, flanks = gap | | | 0.983 | 53 |
| DIAMOND `--masking 0`, flanks = align | | | 0.988 | 45 |

r dN/dS for DIAMOND split by hit coverage:

| DIAMOND route | coverage ≥ 0.9 (964 pairs) | coverage < 0.9 (34 pairs) |
|---|---:|---:|
| flanks = gap | 0.971 | 0.600 |
| flanks = align | 0.973 | 0.950 |

**What the benchmark shows:**

- **Skipping realignment does change dN/dS.** It changes it more than any choice among global aligners.
- **Gap costs explain little of it.** NW with DIAMOND's own gap costs stays at r = 0.998.
- **DIAMOND's block itself is the cause.** In 409 of the 964 well-covered pairs it differs from the optimal
  alignment of the same region, 219 of them because of masking.
- **Partial hits make it worse** when the flanks are gapped (r = 0.60 below 90% coverage).
- **The time saved is negligible.** NW takes about 0.4 s of CPU per 1,000 pairs, against seconds for the DIAMOND
  search itself.

## Recommended defaults

- **Pairs: `coal-cpp-align` / `coal::align_pairs`, BLOSUM62, open 10, extend 0.5, end gaps free, full length.** It
  is exact and fast enough that alignment is no longer the bottleneck.
  - Use `--preset pwalign` only to reproduce a pwalign-based orthologr run.
- **Clusters of 3 or more: Clustal Omega**, via `align_clusters(jobs = CPUs, threads = 1)`, with pairs by NW
  (`pairs_with="nw"`). Use MAFFT (`tool="mafft"`) as the alternative for very large clusters, where its per-call
  overhead does not matter.
- **The DIAMOND route: optional, not the default.** Offer it as a "no realignment" mode only with
  `--masking 0 --flanks align`, and say that results differ (r ≈ 0.99 against NW on this set).

## Known gaps

- The platform wheels ship `coal-cpp` but not yet `coal-cpp-align` (`packaging/wheels/build_platform_wheel.py`
  copies one binary). Without it, `coal.align.align_pairs` falls back to pure Python.
- There is no alignment-only subcommand in `coal`/`coal-cpp`. `coal pipeline` ([pipeline.md](pipeline.md)) runs
  the aligner as one of its steps; alone, it is the `coal-cpp-align` tool and the library and Python APIs.
- **No Clustal Omega (or MAFFT) install on Windows** (see "Install"): cluster mode with clusters of three or more
  needs WSL or a tool the user provides. On Linux and macOS the install is Bioconda's, through conda / mamba /
  micromamba.
- The AVX2 path is x86 only. ARM (Apple Silicon, Linux aarch64) gets the compiler-vectorised baseline and the scalar
  scan, and was not measured.
- No linear-space (Hirschberg) fallback: memory is `n·m/2` bytes per pair in flight, times `-j`.
- Thread scaling and the headline numbers should be rerun on the 16-thread PC.
