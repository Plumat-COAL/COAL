# coal-align-app: the COAL web GUI and MCP server

The app for [COAL](https://github.com/plumat-coal/COAL), fast codon alignments (a PAL2NAL alternative), Ka/Ks
(dN/dS) for pairs, and a one-step pipeline from coding sequences to dN/dS; free software under GPL-3.0-only:

- a web GUI on this computer (`coal-gui`, or `python -m coal_app`) with four tabs:
  - **Align**: a codon alignment from a protein alignment and the coding DNA, with a per-gene health check, an
    alignment view that marks where a gene breaks, genetic-code suggestions and a ready-to-paste methods sentence;
  - **Ka/Ks**: dN, dS and dN/dS for pairs of aligned coding sequences (an aligned codon FASTA, every pair or
    consecutive pairs, or KaKs_Calculator's AXT; uploaded or pasted), with the methods YN (the default), NG, LWL,
    MLWL, LPB, MLPB, MYN, Comeron and Li. The table shows each result's status (`ok`, `warning`, `not_converged`,
    `saturated`, `no_data`, `error`) with what it means, and downloads as TSV (`coal dnds`,
    [docs/kaks.md](../../docs/kaks.md));
  - **dN/dS pipeline**: two FASTA files of coding sequences, one per species: translation, ortholog pairs by
    DIAMOND reciprocal best hits, pairwise alignment, codon alignment and Ka/Ks in one run, with its stages and
    progress lines shown as it goes and a Stop button; `dnds.tsv`, `report.json` and a zip of all result files to
    download (`coal pipeline pairs`, [docs/pipeline.md](../../docs/pipeline.md)). It runs only when DIAMOND is
    found (`$COAL_DIAMOND`, `PATH` or COAL's cache); otherwise the tab says how to install it:
    `coal pipeline install diamond`;
  - **Settings**: the defaults of the Align tab.
- an MCP server for AI agents (`python -m coal_app --mcp`), with the tools `codon_align`, `explain_report`,
  `dnds_pairs` (Ka/Ks on a codon alignment or AXT file) and `pipeline_pairs` (two CDS files to dN/dS per ortholog
  pair; DIAMOND, or a table of given pairs);
- the Windows app (`desktop/`) and the Hugging Face Space (`space/`), built from this package.

Ka/Ks and the pipeline need COAL's C++ engine, `coal-cpp`, which the platform wheels and the Windows app include
(`coal --cpp-path` shows it); an install with the Python engine only makes codon alignments, and the two tabs and
tools say so. The pipeline's cluster mode (`coal pipeline clusters`) and KaKs_Calculator's output format
(`coal kaks`) are in the `coal` command line only.

**Limits.** The web GUI reads only uploaded files, pasted text and the bundled examples in the Ka/Ks and pipeline
tabs (no typed paths). A Ka/Ks input can be up to 100 MB and 500,000 results (pairs × methods); a pipeline input
up to 1,000 MB per file (one run at a time). The MCP tools take Ka/Ks inputs of the same size and pipeline inputs
up to 500 MB per file and one hour. The shared Hugging Face Space takes Ka/Ks inputs up to 5 MB and 20,000 results
and does not run the pipeline. Larger runs belong on the command line. Results are written to a run folder that
is deleted after 24 hours; the pipeline's intermediate files (its `work` folder) are removed when a run ends.

The engine is a separate package, `coal-align` (the repository root), which this one depends on at exactly the same
version. Pipelines that only need the `coal` / `pal2nal` commands install `coal-align` alone.

## Install

Neither package is on PyPI yet. From the latest release:

```
pip install "coal-align @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_src.tar.gz" \
    "coal-align-app[gui] @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_app_src.tar.gz"
```

`coal_src.tar.gz` is the Python engine only (codon alignments); for the Ka/Ks and pipeline tabs install the
platform wheel of `coal-align` instead (`coal_align-<version>-py3-none-<platform>.whl` from the same release, which
includes `coal-cpp`), and for the pipeline DIAMOND: `coal pipeline install diamond`.

From a checkout, at the repository root: `pip install -e . -e "apps/coal-app[gui]"` (then build `coal-cpp`,
[docs/cpp-engine.md](../../docs/cpp-engine.md), and set `COAL_CPP` for Ka/Ks and the pipeline). Use `[mcp]` instead of `[gui]`
for the stand-alone MCP server, in a separate environment (the two extras pin conflicting `mcp` versions).

## Tests

From the repository root: `python -m unittest discover -s apps/coal-app/tests`, and in `apps/coal-app/desktop/`:
`python -m unittest test_wrapper test_smoke`. More: [docs/install.md](../../docs/install.md),
[docs/mcp.md](../../docs/mcp.md), [docs/development.md](../../docs/development.md).
