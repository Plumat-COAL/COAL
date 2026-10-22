# COAL

**Codon alignments from a protein alignment and its coding DNA.** COAL is a fast, open-source alternative to
[PAL2NAL](https://www.bork.embl.de/pal2nal/): it maps unaligned coding DNA onto an existing protein multiple sequence
alignment to produce a codon alignment, ready for dN/dS (PAML codeml, HyPhy) or DNA phylogenies.

- **Compatible:** COAL's output is identical to PAL2NAL v14's except for a short list of documented, deliberate
  differences ([docs/differences.md](docs/differences.md)). Every stored PAL2NAL v14 reference check passes, and
  random comparisons with the real PAL2NAL v14 find 0 unexplained differences
  ([tests/README.md](tests/README.md#test-numbers)).
- **Fast and scalable:** with the C++ engine (`coal-cpp`, one thread), about **300× faster** than PAL2NAL v14 on
  1 million sequences (12.5 s against 62 min, and 172 MB against 18.6 GB of memory) and over **1,000× faster** on
  long genes (1,000 sequences × 10,000 codons: 0.25 s against 280 s). PAL2NAL v14 cannot process proteins longer
  than about 22,000 residues; COAL can. The Python engine gives byte-identical output, more slowly
  ([docs/benchmarks.md](docs/benchmarks.md)).
- **Pipeline-friendly:** non-zero exit codes on errors (PAL2NAL always exits 0), `--strict` to refuse
  probably-wrong genes, a JSON run report, and a `pal2nal` command that accepts an unchanged PAL2NAL v14 command
  line.
- **Tells you what went wrong:** warnings name the gene and position, with genetic-code suggestions.

## Install

Needs Python 3.10+. The package is called `coal-align` (the name `coal` on PyPI belongs to another project); the
command and the Python module are both `coal`. Until it is on PyPI, install from a release or from source:

```
python -m pip install "coal-align @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_src.tar.gz"
```

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

Example inputs are in `examples/` (`clean`, `mito`, `messy_ids`, `one_broken_gene`). Every option:
`coal --help` and [docs/usage.md](docs/usage.md).

## Using COAL in a pipeline

- **Replacing PAL2NAL:** `pal2nal pep.aln nuc.fasta [PAL2NAL options]` accepts an unchanged PAL2NAL v14 command line
  (clustal output by default, as PAL2NAL). To call COAL's own options, use `coal`, whose default output is FASTA.
- **Exit codes:** 0 success; 1 bad input or a gene that cannot be converted; 2 usage error (`coal`), or some
  sequences left out with `--keep-going`; 141 closed output pipe; 130 / 143 / 129 interrupted. Messages go to stderr,
  each line starting `coal: error:`, `coal: warning:`, `coal: notice:` or `coal: summary:`.
- **Fail rather than continue:** like PAL2NAL, COAL writes genes with many mismatching codons (with a notice).
  `--strict` makes those an error, for pipelines that would rather stop.
- **Machine-readable results:** `--report run.json` (schema `coal.report/1`) gives every sequence's status, the
  pairing and a genetic-code suggestion, also when the run fails ([docs/usage.md](docs/usage.md#run-report-and-long-runs)).
- **Parallel jobs:** with many COAL processes running at once (one per gene family), pass `-j 1` so each uses one
  core. Without it, COAL uses all cores for large inputs.
- **Reading from a pipe:** pass `/dev/stdin` as an input path.
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
| Speed and memory vs PAL2NAL | [docs/benchmarks.md](docs/benchmarks.md) |
| Windows app, web GUI, other install options | [docs/install.md](docs/install.md), [apps/coal-app/desktop/README.md](apps/coal-app/desktop/README.md) |
| Running COAL from an AI agent (MCP) | [docs/mcp.md](docs/mcp.md) |
| Tests, the PAL2NAL reference, releasing | [docs/development.md](docs/development.md), [tests/README.md](tests/README.md) |
| Contributing: setup, the two-engine rule, adding a deliberate difference | [CONTRIBUTING.md](CONTRIBUTING.md) |
| Changes between versions | [CHANGELOG.md](CHANGELOG.md) |

## Citation and licence

If you use COAL, please cite it (`CITATION.cff`) and PAL2NAL:

Suyama M, Torrents D, Bork P (2006) PAL2NAL: robust conversion of protein sequence alignments into the
corresponding codon alignments. *Nucleic Acids Res.* 34:W609–W612.

COAL is released under the MIT License (`LICENSE`). It contains no PAL2NAL code; PAL2NAL (GPL-2.0) is used only as an
external reference for testing and is not redistributed here.

COAL was started by Chloe Newton and Keanu Murati at the Reinvent the Wheel AI × Science Tool-building Hackathon,
London, October 2026.
