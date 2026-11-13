# COAL

**Codon alignments from a protein alignment and its coding DNA, Ka/Ks (dN/dS) for pairs, and one step from coding
sequences to a dN/dS table.** COAL is a fast, open-source alternative to
[PAL2NAL](https://www.bork.embl.de/pal2nal/): it maps unaligned coding DNA onto an existing protein multiple sequence
alignment to produce a codon alignment, ready for dN/dS (PAML codeml, HyPhy) or DNA phylogenies. From 0.6 it also
estimates Ka/Ks itself and runs the whole route from two CDS files to dN/dS.

| What | Command | Needs |
|---|---|---|
| Codon alignment from a protein alignment + coding DNA | `coal pep.aln nuc.fasta` (or `pal2nal ...`) | any install |
| Ka/Ks (dN/dS) for pairs of aligned coding sequences | `coal dnds` (tidy table), `coal kaks` (KaKs_Calculator's format) | the platform wheel (`coal-cpp`) |
| Coding sequences in, dN/dS table out, in one step | `coal pipeline pairs` (orthologs of two species), `coal pipeline clusters` (gene families) | the platform wheel, DIAMOND (and Clustal Omega for clusters) |
| The same from R | the `coalr` package; orthologr's `dNdS(engine = "coal")` | R ([docs/r.md](docs/r.md)) |

- **Compatible:** COAL's output is identical to PAL2NAL v14's except for a short list of documented, deliberate
  differences ([docs/differences.md](docs/differences.md)). Every stored PAL2NAL v14 reference check passes, and
  random comparisons with the real PAL2NAL v14 find 0 unexplained differences
  ([tests/README.md](tests/README.md#test-numbers)).
- **Fast and scalable:** with the C++ engine (`coal-cpp`, one thread), about **300× faster** than PAL2NAL v14 on
  1 million sequences (12.5 s against 62 min, and 172 MB against 18.6 GB of memory) and over **1,000× faster** on
  long genes (1,000 sequences × 10,000 codons: 0.25 s against 280 s), measured on Windows; on Linux, where Perl
  runs faster, the same desktop gave about 120× at 1 million sequences and 410× on 100 × 10,000 codons
  ([docs/benchmarks.md](docs/benchmarks.md)). PAL2NAL v14 cannot process proteins longer than about 22,000
  residues; COAL can. The Python engine gives byte-identical output, more slowly
  ([docs/benchmarks.md](docs/benchmarks.md)).
- **Pipeline-friendly:** non-zero exit codes on errors (PAL2NAL always exits 0), `--strict` to refuse
  probably-wrong genes, a JSON run report, and a `pal2nal` command that accepts an unchanged PAL2NAL v14 command
  line.
- **Tells you what went wrong:** warnings name the gene and position, with genetic-code suggestions.
- **Ka/Ks without a second program:** `coal kaks` reads and writes KaKs_Calculator 2.0's formats (its methods NG,
  LWL, MLWL, LPB, MLPB, YN, MYN), `coal dnds` writes a tidy table with standard errors and a status per pair, and
  both have Comeron (1995) and Li (1993) as the R package orthologr computes them
  ([docs/kaks.md](docs/kaks.md): what was validated against what).
- **CDS to dN/dS in one step:** `coal pipeline` translates, finds orthologs (DIAMOND reciprocal best hits) or
  clusters (DIAMOND DeepClust), aligns the proteins, threads the codons and estimates Ka/Ks, with a report that
  records versions, command lines and input checksums ([docs/pipeline.md](docs/pipeline.md)).

## Install

Needs Python 3.10+. The package is called `coal-align` (the name `coal` on PyPI belongs to another project); the
command and the Python module are both `coal`. **Publication on PyPI is pending**, so install from a release. There
are two routes, and they do not give the same thing:

| Route | What you get |
|---|---|
| **The platform wheel** (Linux x86_64 / aarch64 with glibc 2.28+, macOS 11+ on Apple silicon and Intel, Windows x64): download `coal_align-<version>-py3-none-<platform>.whl` from the [Releases page](https://github.com/plumat-coal/COAL/releases/latest), then `python -m pip install` that file | **Everything**: codon alignment with both engines, Ka/Ks (`coal kaks`, `coal dnds`) and `coal pipeline`. The C++ engine `coal-cpp` is inside the wheel |
| **The source package** (the command below), the `py3-none-any` wheel, a checkout without a build | **The Python engine only: codon alignment** (`coal`, `pal2nal`, `import coal`). `coal kaks`, `coal dnds` and `coal pipeline` say that they need `coal-cpp` and how to get it |

```
python -m pip install "coal-align @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_src.tar.gz"
```

`coal --cpp-path` prints where `coal-cpp` is (exit 1 when this install has none); `coal` hands large inputs to it
with no settings. `coal pipeline` also runs **DIAMOND** (2.1.0 or later) and, for clusters, **Clustal Omega** or
MAFFT; they are separate programs, found on `PATH` or installed into COAL's cache, checksum-verified, with

```
coal pipeline install diamond clustalo
```

Once it is on PyPI, `python -m pip install coal-align` will pick the platform wheel by itself.

From source, with the C++ engine (needs a C++17 compiler and CMake 3.16+; optional but recommended for large inputs
and for many small runs):

```
git clone https://github.com/plumat-coal/COAL && cd COAL
python -m pip install -e .
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --config Release
export COAL_CPP=$PWD/build/coal-cpp       # coal hands large inputs to coal-cpp when it finds it
```

The engine (`coal-align`) has no dependencies and no app code: it is all a pipeline needs. The web GUI, the MCP
server for AI agents and the Windows app are a second package in this repository, `coal-align-app`
([apps/coal-app/](apps/coal-app/)), installed next to the engine:

```
python -m pip install "coal-align @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_src.tar.gz" \
    "coal-align-app[gui] @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_app_src.tar.gz"
```

The Windows app, the web GUI and the other install options: [docs/install.md](docs/install.md).

## Quick start

```
coal pep.aln nuc.fasta > codon_aln.fasta                  # FASTA protein alignment (or CLUSTAL) + coding DNA
coal pep.aln nuc.fasta -codontable 2 -nogap -o out.fasta  # PAL2NAL options work
coal pep.aln nuc.fasta -output paml -o out.paml           # fasta (default), clustal, paml, codon
coal pep.aln nuc.fasta -o out.fasta --report run.json --keep-going   # per-gene report; skip genes that fail
```

Example inputs are in `examples/` (`clean`, `mito`, `messy_ids`, `one_broken_gene`, and `dnds` for Ka/Ks and the
pipeline). Every option: `coal --help` (it lists every command; `coal help COMMAND`) and
[docs/usage.md](docs/usage.md). One runnable example of each thing COAL does, from a checkout:

```
# 1. codon alignment: a protein alignment + the coding DNA of its genes
coal examples/clean/protein.aln.fasta examples/clean/dna.fasta -o codon.fasta

# 2. Ka/Ks (dN/dS) for pairs of a codon alignment (needs coal-cpp: the platform wheel)
coal dnds examples/dnds/codon.fasta --pairs consecutive -o dnds.tsv           # tidy table
coal kaks -i examples/dnds/codon.fasta -o out.kaks --pairs consecutive -m YN   # KaKs_Calculator's format

# 3. one step: the coding sequences of two species in, dN/dS of their orthologs out (needs DIAMOND too)
coal pipeline pairs --query examples/dnds/species_a_cds.fasta --subject examples/dnds/species_b_cds.fasta -o out
#    -> out/dnds.tsv (12 pairs), out/orthologs.tsv, out/genes.tsv, out/report.json; `coal cite out/report.json`
```

```r
# 4. from R: the coalr package, in process (no Python, no DIAMOND; the pairs are given)
library(coalr)
q <- coal_read_fasta("examples/dnds/species_a_cds.fasta")
s <- coal_read_fasta("examples/dnds/species_b_cds.fasta")
coal_dnds_pairs(q, s, pairs = read.delim("examples/dnds/orthologs.tsv")[, 1:2])

# 5. inside orthologr (the fork with COAL built in; packaging/orthologr/README.md)
orthologr::dNdS("examples/dnds/species_a_cds.fasta", "examples/dnds/species_b_cds.fasta",
                dnds_est.method = "YN", engine = "coal")
```

**In the pipeline every FASTA record is a gene.** A CDS file with several transcripts per gene (Ensembl, NCBI,
TAIR) gives several rows per gene and splits the reciprocal best hits; the run warns about it, and
`--longest-isoform` (with `--gene-map TSV` if the headers do not name the gene) keeps the longest CDS of each gene
([docs/pipeline.md](docs/pipeline.md#one-record-is-one-gene)). Before reading the values, see
[Filtering before interpretation](docs/kaks.md#filtering-before-interpretation): short alignments, saturated dS and
very small dS give numbers that mean little.

## Using COAL in a pipeline

- **Replacing PAL2NAL:** `pal2nal pep.aln nuc.fasta [PAL2NAL options]` accepts an unchanged PAL2NAL v14 command line
  (clustal output by default, as PAL2NAL). To call COAL's own options, use `coal`, whose default output is FASTA.
- **orthologr:** its `dNdS()` runs a bundled `pal2nal.pl`; `packaging/orthologr/` has a stand-in for that file,
  notes on what COAL changes there, and the fork whose `dNdS(engine = "coal")` hands all pairs to COAL at once.
- **R:** `remotes::install_github("Plumat-COAL/COAL", subdir = "packaging/r/coalr")` installs `coalr`: dN/dS for
  ortholog pairs or clusters, Ka/Ks, codon alignment and protein alignment inside the R process
  ([docs/r.md](docs/r.md)).
- **Exit codes:** 0 success; 1 bad input or a gene that cannot be converted; 2 usage error (`coal`), or some
  sequences left out with `--keep-going`; 141 closed output pipe; 130 / 143 / 129 interrupted. Messages go to stderr,
  each line starting `coal: error:`, `coal: warning:`, `coal: notice:` or `coal: summary:`.
- **Fail rather than continue:** like PAL2NAL, COAL writes genes with many mismatching codons (with a notice).
  `--strict` makes those an error, for pipelines that would rather stop.
- **Machine-readable results:** `--report run.json` (schema `coal.report/1`) gives every sequence's status, the
  pairing and a genetic-code suggestion, also when the run fails ([docs/usage.md](docs/usage.md#run-report-and-long-runs)).
- **Parallel jobs:** with many COAL processes running at once (one per gene family), pass `-j 1` so each uses one
  core. Without it, COAL uses all cores for large inputs.
- **Thousands of small alignments:** `coal --batch jobs.txt` runs one job per line of `jobs.txt` (a `coal` command
  line each, every one with `-o`) in a single process, so Python starts once instead of once per gene: about 40×
  less wall time on tiny jobs, with byte-identical outputs; `--batch-parallel N` runs N jobs at a time
  ([docs/usage.md](docs/usage.md#batch-mode)).
- **Reading from a pipe:** pass `-` as the protein alignment or one DNA file (`zcat nuc.fa.gz | coal pep.aln -`);
  works the same on Windows, with `coal-cpp` and with the `pal2nal` drop-in. Only one `-` per command line.
- **From Python** ([docs/python-api.md](docs/python-api.md): `codon_align`, `align_files`, `CoalError`, and what 1.0 will
  keep stable):

  ```python
  import coal
  rows, warnings = coal.codon_align(
      [("a", "MK-V"), ("b", "MKLV")],                   # (id, aligned protein row)
      [("a", "ATGAAAGTT"), ("b", "ATGAAACTGGTT")],      # (id, coding DNA)
      table=1)
  # rows == [("a", "ATGAAA---GTT"), ("b", "ATGAAACTGGTT")]; raises coal.CoalError on bad input
  ```

## Documentation

| Topic | Where |
|---|---|
| Every option, input and output rules, large inputs, run report | [docs/usage.md](docs/usage.md) |
| Using COAL from Python (`codon_align`, `align_files`, `CoalError`) and the 1.0 versioning promise | [docs/python-api.md](docs/python-api.md) |
| The run report's fields (`--report` / `--json`, schema `coal.report/1`) and its JSON Schema | [docs/report-schema.md](docs/report-schema.md), [docs/report.schema.json](docs/report.schema.json) |
| How and why COAL's output differs from PAL2NAL v14 | [docs/differences.md](docs/differences.md) |
| The C++ engine | [docs/cpp-engine.md](docs/cpp-engine.md) |
| Ka/Ks (dN/dS): `coal kaks` (KaKs_Calculator's format, seven of its methods; Comeron 1995 and Li 1993 as the R package orthologr computes them), `coal dnds`, policy, validation, where the code comes from | [docs/kaks.md](docs/kaks.md) |
| One step from CDS to dN/dS: `coal pipeline pairs` (DIAMOND reciprocal best hits) and `coal pipeline clusters` (DeepClust, per-cluster summary), output tables, validation, timings | [docs/pipeline.md](docs/pipeline.md), [docs/pipeline.schema.json](docs/pipeline.schema.json) |
| R: the `coalr` package (in-process `coal_dnds_pairs()`, `coal_dnds_clusters()`, `coal_kaks()`, `coal_codon_align()`; the whole pipeline through reticulate), install, examples, when to use which | [docs/r.md](docs/r.md), [packaging/r/coalr](packaging/r/coalr) |
| Protein alignment before the codon step: pairwise (BLOSUM62 NW), Clustal Omega / MAFFT for clusters, DIAMOND hits | [docs/alignment.md](docs/alignment.md) |
| Speed and memory vs PAL2NAL | [docs/benchmarks.md](docs/benchmarks.md) |
| Windows app, web GUI, other install options | [docs/install.md](docs/install.md), [apps/coal-app/desktop/README.md](apps/coal-app/desktop/README.md) |
| Running COAL from an AI agent (MCP) | [docs/mcp.md](docs/mcp.md) |
| Tests, the PAL2NAL reference, releasing | [docs/development.md](docs/development.md), [tests/README.md](tests/README.md) |
| Contributing: setup, the two-engine rule, adding a deliberate difference | [CONTRIBUTING.md](CONTRIBUTING.md) |
| Changes between versions | [CHANGELOG.md](CHANGELOG.md) |

## Citation and licence

If you use COAL, please cite it (`CITATION.cff`, with the version `coal --version` prints) and PAL2NAL:

Suyama M, Torrents D, Bork P (2006) PAL2NAL: robust conversion of protein sequence alignments into the
corresponding codon alignments. *Nucleic Acids Res.* 34:W609–W612. <https://doi.org/10.1093/nar/gkl315>

COAL's Ka/Ks numbers come from published estimators, in code derived from KaKs_Calculator 2.0, PAML, libsequence
and seqinr, and `coal pipeline` runs DIAMOND and an aligner. Please cite what your numbers come from:

| What you ran | Cite |
|---|---|
| Anything | COAL |
| A codon alignment (`coal PEP NUC`, `pal2nal`, every Ka/Ks route of the pipeline) | Suyama et al. 2006 (PAL2NAL) |
| Ka/Ks, method NG | Nei & Gojobori 1986; Wang et al. 2010 (KaKs_Calculator 2.0) |
| Ka/Ks, method LWL | Li, Wu & Luo 1985; Wang et al. 2010 (KaKs_Calculator 2.0) |
| Ka/Ks, method MLWL | Li, Wu & Luo 1985; Tzeng, Pan & Li 2004; Wang et al. 2010 (KaKs_Calculator 2.0) |
| Ka/Ks, method LPB | Li 1993; Pamilo & Bianchi 1993; Wang et al. 2010 (KaKs_Calculator 2.0) |
| Ka/Ks, method MLPB | Li 1993; Pamilo & Bianchi 1993; Tzeng, Pan & Li 2004; Wang et al. 2010 (KaKs_Calculator 2.0) |
| Ka/Ks, method YN (the default) | Yang & Nielsen 2000; Wang et al. 2010 (KaKs_Calculator 2.0); Yang 2007 (PAML 4) |
| Ka/Ks, method MYN | Zhang, Li & Yu 2006; Yang & Nielsen 2000; Wang et al. 2010 (KaKs_Calculator 2.0); Yang 2007 (PAML 4) |
| Ka/Ks, method Comeron | Comeron 1995; Thornton 2003 (libsequence) |
| Ka/Ks, method Li | Li 1993; Charif & Lobry 2007 (seqinr) |
| `coal pipeline pairs` with DIAMOND (orthologs) | Buchfink et al. 2021 (DIAMOND) |
| `coal pipeline clusters` with DIAMOND DeepClust | Buchfink et al. 2021 (DIAMOND); Buchfink et al. 2026 (DIAMOND DeepClust) |
| Protein alignments by Clustal Omega (`--aligner clustalo`, the cluster default) | Sievers et al. 2011 (Clustal Omega) |
| Protein alignments by MAFFT (`--aligner mafft`) | Katoh & Standley 2013 (MAFFT) |
| Through the R package orthologr | Drost et al. 2015 (orthologr) |

`coal cite` prints this table and the full references; `coal cite out/report.json` prints the references of one
pipeline run (its report has them as `cite`). The references are in
[docs/kaks.md](docs/kaks.md#what-to-cite).

COAL is free software under the GNU General Public License, version 3 only (SPDX `GPL-3.0-only`; the text is in
`LICENSE`, the copyright notice in `NOTICE`). This applies from the release after 0.5.0: **0.5.0 and earlier were
released under the MIT License and stay MIT.** In practice: running COAL puts no condition on your data, results,
scripts or pipelines; if you copy COAL's code into your own program or link its library (`coal::core`,
`import coal`), what you distribute has to be under a GPL-3-compatible licence. The reason for the change is in
[CHANGELOG.md](CHANGELOG.md). COAL contains no PAL2NAL code; PAL2NAL (GPL-2.0) is used only as an external reference
for testing and is not redistributed here.

COAL was started by Chloe Newton and Keanu Murati at the Reinvent the Wheel AI × Science Tool-building Hackathon,
London, October 2026.
