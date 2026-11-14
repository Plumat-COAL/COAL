# AI agents (MCP)

COAL's stdio MCP server gives an agent four tools:

| Tool | What it does | Needs |
|---|---|---|
| `codon_align` | protein alignment + coding DNA -> codon alignment, with a report; `save_as` puts the alignment where you want it | the engine (any install) |
| `explain_report` | that report in plain English | the engine (any install) |
| `dnds_pairs` | Ka/Ks (dN, dS and dN/dS) for pairs of an aligned codon FASTA or a KaKs_Calculator AXT file (`coal dnds`) | `coal-cpp` |
| `pipeline_pairs` | the coding sequences of two species -> ortholog pairs (DIAMOND reciprocal best hits, or a table of given pairs) -> pairwise alignment, codon alignment and dN/dS per pair (`coal pipeline pairs`) | `coal-cpp`, and DIAMOND unless the pairs are given |

`coal-cpp` is COAL's C++ engine: the platform wheels of `coal-align` and the Windows app include it (`coal --cpp-path`
shows it); an install with the Python engine only makes codon alignments, and `dnds_pairs` / `pipeline_pairs` then
answer with a tool error that says so. DIAMOND is not included anywhere: `coal pipeline install diamond` puts it in
COAL's cache (COAL looks at `$COAL_DIAMOND`, then `PATH`, then that cache).

**Windows app:** keep **Connect COAL to Claude Desktop** ticked on the installer's last page (or run
`COAL.exe --connect-claude` later), then quit and reopen Claude Desktop. This adds a `coal` entry that starts
`COAL.exe -m coal_app --mcp` to `claude_desktop_config.json`, keeping everything else in it (the old file is kept as
`claude_desktop_config.json.bak-coal`); uninstalling COAL removes the entry again. The app includes `coal-cpp`; for
`pipeline_pairs` install DIAMOND once with `"C:\path\to\COAL.exe" -m coal pipeline install diamond` (the tool's error
message gives the exact command).

**From source:** the MCP server is part of the app package (`apps/coal-app`, `coal-align-app`); install it with
the engine in its **own** virtual environment (not the one with the `gui` extra), from the repository root:

```
python -m venv .venv-mcp
.venv-mcp/bin/python -m pip install -e . -e "apps/coal-app[mcp]"   # Windows: apps\coal-app\desktop\setup_mcp.bat
claude mcp add coal -- "$PWD/.venv-mcp/bin/python" -m coal_app --mcp
```

That is the Python engine: for `dnds_pairs` and `pipeline_pairs` also build `coal-cpp`
([cpp-engine.md](cpp-engine.md)) and put it on `PATH` or set `COAL_CPP` to it in the server's environment.

Without a checkout: `pip install "coal-align @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_src.tar.gz" "coal-align-app[mcp] @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_app_src.tar.gz"`.
`coal_src.tar.gz` is the Python engine only; for Ka/Ks and the pipeline install the platform wheel of `coal-align`
from the same release instead (`coal_align-<version>-py3-none-<platform>.whl`, which includes `coal-cpp`).

For Claude Desktop, add this to `claude_desktop_config.json` (Settings > Developer > Edit Config), with the full path
of that environment's Python (on Windows, `apps\coal-app\desktop\setup_mcp.bat` prints it ready to paste,
backslashes escaped):

```json
{"mcpServers": {"coal": {"command": "/full/path/to/COAL/.venv-mcp/bin/python", "args": ["-m", "coal_app", "--mcp"]}}}
```

(add `"env": {"COAL_CPP": "/full/path/to/coal-cpp"}` next to `"args"` if `coal-cpp` is not on that `PATH`).

The web GUI (`python -m coal_app`, `coal-gui`) also serves the same four tools at
`http://127.0.0.1:7860/gradio_api/mcp/sse`; the Windows app does not by default (`COAL.exe --mcp`). Agent
instructions: `.claude/skills/coal/SKILL.md`; a stage demo: `demo/mcp_demo.md`.

## How the tools answer

Give the tools **absolute paths** (a relative path resolves against the server's working directory, which the error
names). Every tool writes into a fresh run folder of its own (in `coal_app_runs` in the system's temporary
folder, or `$COAL_RUNS_DIR`; run folders older than 24 hours are deleted when the server starts) unless you name a
place: `save_as` (`codon_align`, `dnds_pairs`: a file path, or an existing folder to save in under COAL's file
name) or `out_dir` (`pipeline_pairs`: a new or empty folder). COAL never overwrites a file.

A result says what happened truthfully: `completed` only means the run finished and wrote its output; whether the
data can be used is a separate field (`all_ok` for `codon_align`, `all_estimated` for `dnds_pairs`, `exit_code` and
the status counts for `pipeline_pairs`), and `next_step` says in one sentence what to do now. Something the caller
can fix before the run starts (a missing file, an unknown genetic code or method, too large an input, a missing
`coal-cpp` or DIAMOND, a `save_as` / `out_dir` that is not allowed) is a **tool error** (`isError: true`, one plain
sentence), not a result. A pair or gene that fails during a run is never an error: it is a row with a status, and
the run goes on.

Every Ka/Ks row has a status: `ok`; `warning` (estimated, but read the note: internal stop codons skipped, or a
filter you set flags it); `not_converged`, `saturated`, `no_data` or `error` (no estimate, or only part of one; dN or
dS is `NA`). Read the status of every row before using its numbers; the result's `statuses` lists what each status
that occurs means ([kaks.md](kaks.md#output)).

## `codon_align`

Arguments: `protein_alignment_path` and `dna_fasta_path` (required: a protein alignment, FASTA or CLUSTAL, and the
unaligned coding DNA, paired by ID); `output_format` (`fasta`, the default, `clustal`, `paml`, `codon`);
`codontable` (an NCBI genetic code: 1-6, 9-16, 21-23; default 1); `nogap`, `nomismatch`, `strict`, `blockonly`
(as on the command line, [usage.md](usage.md)); `save_as`.

```json
{"protein_alignment_path": "/data/genes.aln.fasta", "dna_fasta_path": "/data/genes.cds.fasta", "codontable": 1,
 "save_as": "/data/"}
```

Returns `completed` (COAL finished), `all_ok` / `needs_attention` (whether every sequence matched its protein and no
hint says the data looks wrong), `counts` (total / ok / warning / failed / not_checked), `suggested_codontable` (a
better genetic code, or null: re-run with it), `next_step`, `output_path` and `report_path` (a `coal.report/1`
JSON file, [report-schema.md](report-schema.md)), the first ten `problems` (`problems_total`), `hints`,
`methods_text` (a ready-to-paste methods sentence), `coal_command` and `pal2nal_command` (the closest PAL2NAL v14
command line). Use the alignment only if `all_ok` is true; otherwise call `explain_report`.

## `explain_report`

Argument: `report_path` (a `coal.report/1` JSON file: `codon_align`'s `report_path`, or `coal --report`). Returns
plain English: the counts, how proteins and DNA were paired, one line per problem sequence with its likely cause
(wrong genetic code, mismatched IDs, frameshift), fixes and the methods sentence. A report that cannot be read is a
tool error. (The pipeline's `report.json` is another schema, `coal.pipeline/1`; `pipeline_pairs` already returns
its summary.)

```json
{"report_path": "/tmp/coal_app_runs/run_x1y2z3/genes.report.json"}
```

## `dnds_pairs`

dN, dS and dN/dS for pairs of aligned coding sequences: a codon alignment you already have (for example
`codon_align`'s `output_path`) or an AXT file. For a protein alignment plus DNA, call `codon_align` first; for two
files of unaligned coding sequences, use `pipeline_pairs`.

| Argument | Default | |
|---|---|---|
| `codon_alignment_path` | (required) | an aligned codon FASTA or a KaKs_Calculator AXT file (told apart by the first character, `>` = FASTA) |
| `pairs` | `all` | FASTA: `all` = every pair (n(n-1)/2), `consecutive` = records 1+2, 3+4, ...; ignored for AXT, which is read pair by pair |
| `methods` | `["YN"]` | any of `YN`, `NG`, `LWL`, `MLWL`, `LPB`, `MLPB`, `MYN`, `Comeron`, `Li` (one row per pair and method; Comeron and Li: genetic code 1 only) |
| `codontable` | 1 | NCBI genetic code, as for `codon_align` |
| `min_codons` | 0 (off) | flag (status `warning`; nothing is dropped) estimates from fewer compared codons |
| `max_ds` | 0 (off) | flag (status `warning`) estimates with dS above this (saturation) |
| `save_as` | | a file path or an existing folder for a copy of the table |

```json
{"codon_alignment_path": "/home/me/COAL/examples/dnds/codon.fasta", "pairs": "consecutive", "methods": ["YN"]}
```

On that example (12 ortholog pairs of two synthetic species) the result is, shortened:

```json
{"completed": true, "all_estimated": true, "exit_code": 0, "error": "",
 "counts": {"total": 12, "ok": 12, "warning": 0, "not_converged": 0, "saturated": 0, "no_data": 0, "error": 0},
 "pairs": 12, "methods": ["YN"], "input_format": "fasta",
 "output_path": "/tmp/coal_app_runs/run_x1y2z3/codon.dnds.tsv",
 "rows": [{"seq1": "spA_g01", "seq2": "spB_g01", "method": "YN", "dN": "0.010176579", "dS": "0.39934548",
           "dNdS": "0.025483147", "SE_dN": "0.0038534006", "SE_dS": "0.053627934", "codons_compared": "320",
           "status": "ok", "note": ""}, "..."],
 "rows_total": 12, "flagged": [], "flagged_total": 0,
 "statuses": [{"status": "ok", "count": 12, "meaning": "Estimated. Gaps, ambiguous bases and ..."}],
 "cite": ["coal", "yang2000", "wang2010", "yang2007"], "references": "[coal] Newton C, Murati K (2026). ...",
 "next_step": "Every result is estimated; the table is at output_path. Cite the references when you publish.",
 "coal_command": "coal dnds /home/me/COAL/examples/dnds/codon.fasta -o ... -m YN -c 1 --pairs consecutive",
 "limits": {"max_bytes": 100000000, "max_results": 500000}}
```

`output_path` is the whole table (TSV: `seq1 seq2 method dN dS dNdS SE_dN SE_dS N_sites S_sites Nd Sd
codons_compared status note`, as `coal dnds` writes it, [kaks.md](kaks.md#output)); `rows` holds the first 20 rows
and `flagged` the first ten whose status is not `ok`, with their note. `all_estimated` is true only when every row
is `ok`. When the run finishes but no pair has an estimate (`exit_code` 2; the input is probably protein sequences,
or rows of different lengths), `error` says why and `next_step` points to `codon_align`. `cite` and `references` are
what to cite for this run (`coal cite`); `coal_command` is the same run on the command line, with full paths.

Tool errors: a missing file; an empty input or a FASTA with one sequence; more than 100 MB, or more than 500,000
results (pairs × methods; use `consecutive`, fewer methods or `coal dnds`); an unknown method or genetic code;
Comeron or Li with a code other than 1; `pairs` other than `all` / `consecutive`; a negative `min_codons` /
`max_ds`; a `save_as` that is relative, exists already or is in a folder that does not exist; no `coal-cpp`.

## `pipeline_pairs`

The pipeline's pair mode in one call: two FASTA files of unaligned coding sequences, one per species (may be
gzipped) -> translation -> ortholog pairs by DIAMOND (or a table of given pairs) -> pairwise protein alignment,
codon alignment and Ka/Ks for every pair. The same run as `coal pipeline pairs` ([pipeline.md](pipeline.md#pair-mode)).

| Argument | Default | |
|---|---|---|
| `query_cds_path`, `subject_cds_path` | (required) | the coding sequences of the two species, up to 500 MB each |
| `methods` | `["YN"]` | Ka/Ks methods, as for `dnds_pairs` |
| `codontable` | 1 | NCBI genetic code for translation, the codon alignment and Ka/Ks |
| `longest_isoform` | false | keep one record per gene, its longest CDS (files with several transcripts per gene) |
| `ortho` | `rbh` | `rbh` = reciprocal best hits, `bh` = every query's best hit |
| `pairs_path` | | a two-column TSV of ortholog pairs (query ID, subject ID): those pairs are used and DIAMOND is not needed |
| `min_codons`, `max_ds` | 0 (off) | flags, as for `dnds_pairs` |
| `keep_alignments` | false | also write the codon alignments (`codon_alignments.fasta`) |
| `out_dir` | | a new or empty folder for the results (made if its parent exists); default: a run folder of COAL's |
| `threads` | 0 | 0 = one per CPU |

```json
{"query_cds_path": "/home/me/COAL/examples/dnds/species_a_cds.fasta",
 "subject_cds_path": "/home/me/COAL/examples/dnds/species_b_cds.fasta", "out_dir": "/home/me/results/ab"}
```

The result, shortened:

```json
{"completed": true, "status": "completed", "exit_code": 0, "error": "",
 "counts": {"genes": {"total": 26, "ok": 26, "warning": 0, "failed": 0},
            "unpaired": {"no_hit": 2, "not_reciprocal": 0}, "orthologs": 12, "pairs": {"total": 12, "ok": 12, "...": 0}},
 "pair_counts": {"total": 12, "ok": 12, "warning": 0, "not_converged": 0, "saturated": 0, "no_data": 0, "error": 0},
 "rows": [{"seq1": "spA_g01", "seq2": "spB_g01", "method": "YN", "dN": "0.010176579", "dS": "0.39934548",
           "dNdS": "0.025483147", "codons_compared": "320", "prot_identity": "97.81", "status": "ok", "...": ""}],
 "rows_total": 12, "flagged": [], "flagged_total": 0, "failures": [], "failures_total": 0,
 "out_dir": "/home/me/results/ab",
 "outputs": {"dnds": "/home/me/results/ab/dnds.tsv", "failed": "/home/me/results/ab/failed.tsv",
             "genes": "/home/me/results/ab/genes.tsv", "orthologs": "/home/me/results/ab/orthologs.tsv",
             "unpaired": "/home/me/results/ab/unpaired.tsv"},
 "report_path": "/home/me/results/ab/report.json", "zip_path": "/home/me/results/ab/coal_dnds_results.zip",
 "diamond": "/home/me/.cache/coal/tools/diamond/diamond",
 "cite": ["coal", "suyama2006", "yang2000", "wang2010", "yang2007", "buchfink2021"], "references": "...",
 "next_step": "Every pair is estimated; the table is outputs['dnds'] (and everything in zip_path). Cite references.",
 "coal_command": "coal pipeline pairs --query ... --subject ... -o /home/me/results/ab -m YN -c 1 --ortho rbh --diamond ...",
 "progress": ["coal pipeline: translated 13 query and 13 subject records (0 not usable)", "coal pipeline: 12 pairs (rbh)", "..."],
 "limits": {"max_bytes_per_file": 500000000, "timeout_seconds": 3600}}
```

`status` is `completed`, `failed` or `interrupted` (stopped at the one-hour limit), and `exit_code` is the
pipeline's: 0 nothing failed, 2 some genes or pairs could not be used (see `failures`, the first ten entries of
`failed.tsv`), 1 the run failed (`error` says why, e.g. DIAMOND failed). `counts` are the report's counts of genes,
orthologs and pairs; `pair_counts` the Ka/Ks statuses of `dnds.tsv`. `outputs` names every result file written
(the files are described in [pipeline.md](pipeline.md#output)), `report_path` is `report.json` (schema
`coal.pipeline/1`) and `zip_path` a zip of all of them. `progress` holds the pipeline's last 30 progress lines; the
result comes when the run ends. The work folder (intermediate files) is removed when a run ends.

Tool errors: a missing input or `pairs_path` file; an empty input or one over 500 MB; an unknown method or genetic
code; Comeron or Li with a code other than 1; `ortho` other than `rbh` / `bh`; a negative `threads`, `min_codons` or
`max_ds`; an `out_dir` that is relative, not empty, a file, or in a folder that does not exist; no `coal-cpp`; no
DIAMOND (unless `pairs_path` is given: the message names `coal pipeline install diamond`).

## Limits and the shared server

| | `dnds_pairs` | `pipeline_pairs` |
|---|---|---|
| on this computer | 100 MB, 500,000 results (pairs × methods), one hour | 500 MB per file, one hour |
| shared server (`COAL_HOSTED=1`, the Hugging Face Space) | 5 MB, 20,000 results, only files in the bundled examples, its run folder or Gradio's uploads; no `save_as` | refused |

The limits are checked before the engine starts; `limits` in each result says which apply. Larger runs belong on
the command line (`coal dnds`, `coal pipeline pairs`). The public Space does not serve MCP
at all (the tools take paths on the machine that runs them, so agents run COAL locally), and it has the Python
engine only, so its Ka/Ks needs a local install; on a shared server `codon_align` and `explain_report` read only
those same folders and `codon_align` has no `save_as`.

## Scope: what stays on the command line

The pipeline's cluster mode (`coal pipeline clusters`), KaKs_Calculator's output format (`coal kaks`), the choice
of aligner and inputs over the limits above are not MCP tools: an agent should run the `coal` command line
([kaks.md](kaks.md), [pipeline.md](pipeline.md)) if it has a shell, or tell the person the command. The server's
instructions say so, and tell an agent to call `codon_align` first and then `dnds_pairs` on its `output_path` when
it is asked for dN/dS from a protein alignment and DNA.
