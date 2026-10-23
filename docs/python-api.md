# Python API

The engine package `coal-align` (`pip install coal-align`, no dependencies) can be used from Python as well as from
the command line. Its public API is what `import coal` exports:

| Name | What | Stable for 1.0 |
|---|---|---|
| `coal.codon_align(proteins, dnas, ...)` | in memory: lists of `(id, sequence)` in, codon-aligned rows out | yes |
| `coal.align_files(pep_path, nuc_paths, ...)` | files in, a file (or stdout) out, exactly what the `coal` command does | yes (the arguments listed under [Stable arguments](#stable-arguments)) |
| `coal.CoalError` | the exception both raise for bad input, with `.code` and `.seq_id` | yes |
| `coal.__version__` | the version string, e.g. `"0.4.2"` | yes |

Everything else in the package (`coal.core`, `coal.pipeline`, `coal.report`, `coal.seqio`, ...) is internal and may
change in any release. The `coal.report.RunReport` class is documented below because `align_files` takes one to
produce the run report; only the parts named there are covered by the promise.

The examples on this page were run with COAL 0.5.0rc4; `tests/test_python_api.py` runs them on every change.

## `codon_align`

```python
coal.codon_align(proteins, dnas, table=1, nogap=False, nomismatch=False, strict=False)
    -> (rows, warnings)
```

| Argument | Type | Meaning |
|---|---|---|
| `proteins` | `list[tuple[str, str]]` | `(id, aligned protein row)`, in alignment order. A row holds amino-acid letters (any case), `-` or `.` for a gap, `*` for a stop, the frameshift digits `1`-`9` (that many bases) and PAL2NAL's tfastx marks `\` (1 base) and `/` (2 bases). `B`, `Z`, `J`, `O` and `U` are accepted as `X` (any codon), with a warning. Every row must have the same length. |
| `dnas` | `list[tuple[str, str]]` | `(id, coding DNA)`, the bare sequence: letters only (the file readers strip spaces, digits and line breaks; `codon_align` does not). Upper or lower case, `U` read as `T`. Extra sequence before and after the gene (UTRs, the stop codon) is allowed and dropped. |
| `table` | `int` | NCBI genetic code: one of 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23 (`coal.tables.SUPPORTED`). Any other number raises `ValueError`. |
| `nogap` | `bool` | PAL2NAL `-nogap`: remove every 3-nt column that holds a gap or a stop codon of `table` in any sequence. |
| `nomismatch` | `bool` | PAL2NAL `-nomismatch`: remove every column where a codon does not match its amino acid (and the B/Z/J/O/U columns), and drop the mismatch warnings. |
| `strict` | `bool` | `--strict`: raise `CoalError` instead of writing a probably-wrong row (see [Errors](#errors)). |

**Pairing.** Each protein needs exactly one DNA record (`len(proteins) == len(dnas)`, else `CoalError`
`count_mismatch`). The DNA record for each protein is chosen in this order (the same rule as the command line):

1. **by ID** (`"id"`): every protein ID is a DNA ID. Order and extra records in `dnas` do not matter; the output follows
   `proteins`.
2. **by normalised ID** (`"normalised"`): otherwise, IDs are compared after removing a trailing `.<digits>` version
   suffix, lower-casing, and for `sp|P12345|NAME_HUMAN`-style headers trying each `|` field except the database tags.
   Used only when that gives every protein exactly one DNA record, one-to-one, *and* that pairing differs from file
   order (otherwise PAL2NAL's order pairing gives the same result and COAL keeps PAL2NAL's exact output and messages).
   A warning says so.
3. **by order** (`"order"`): otherwise the k-th protein gets the k-th DNA record, as PAL2NAL does. If some (not all)
   protein IDs are DNA IDs a warning says so; if none are, COAL is silent, like PAL2NAL. With `strict=True` a partial
   exact-ID match is `CoalError` `partial_id_match` instead.

**Return value.** `rows` is `list[tuple[str, str]]`: `(protein_id, codon_aligned_dna)`, one per protein in
`proteins` order. Each residue becomes its codon (the DNA's own letters and case, so synonymous differences are kept),
each gap `---`, each frameshift digit N its N bases padded with `-` to the column width (3, or 6 for digits 4-6,
9 for 7-9); every row has the same length. With `nogap` / `nomismatch` the removed columns are missing from every row.

`warnings` is `list[str]`, in this order: pairing warnings, then per-sequence warnings (in sequence order), then
per-sequence notices. The formats, which are fixed:

| Kind | Format |
|---|---|
| mismatching codon (dropped with `nomismatch`) | `<id> pepAlnPos <column>: <aa> does not correspond to <codon>` |
| B/Z/J/O/U residue | `<id> pepAlnPos <column>: <letter> (<name>) is not checked; any codon is accepted` |
| probably-wrong gene (notice) | `<id>: protein and DNA are probably unrelated: ...`, `<id>: position of the gene in the DNA is a guess: ...`, `<id>: no 10-residue segment of the protein was found in the DNA; ...`, `<id>: placement is ambiguous: ...` |
| pairing | `only <k> of <n> protein IDs are DNA IDs; pairing sequences by order instead ...`, `<k> of <n> protein IDs match a DNA ID only after normalising ...` |

`<column>` is the 1-based column of the protein alignment. These are the same lines the `coal` command prints after
`coal: warning:`.

### Example

```python
>>> from coal import codon_align
>>> proteins = [("g1", "MK-VL"), ("g2", "MKQVL")]
>>> dnas = [("g2", "ATGAAACAGGTTCTG"), ("g1", "ATGAAAGTTCTGTAA")]   # any order; g1 has a stop codon after the gene
>>> rows, warnings = codon_align(proteins, dnas)
>>> rows
[('g1', 'ATGAAA---GTTCTG'), ('g2', 'ATGAAACAGGTTCTG')]
>>> warnings
[]
```

A gene with one wrong codon (PAL2NAL's fallback: the row is written, the codon is reported):

```python
>>> prot = "MFWDKESRSPHESAPQYARKIWEMAAAVAP"
>>> dna = ("ATGTTCTGGGACAAAGAAAGCAGAAGCCCAGATGAAAGCGCACCACAATACGCAAGAAAAATATGGGAAATGGCAGCAGCAGTAGCACCA")
>>> rows, warnings = codon_align([("x", prot)], [("x", dna)])    # column 11 is H, but the codon is GAT (D)
>>> warnings
['x pepAlnPos 11: H does not correspond to GAT']
>>> rows, warnings = codon_align([("x", prot)], [("x", dna)], nomismatch=True)
>>> len(rows[0][1]), warnings                                     # the column is gone, and so is the warning
(87, [])
```

### Errors

`codon_align` raises `coal.CoalError` (a subclass of `Exception`) for bad input; `str(e)` is the message the `coal`
command prints after `coal: error:`. Its attributes:

- `e.code` (`str`): a stable machine-readable reason, the same string as `error.code` in the run report.
- `e.seq_id` (`str | None`): the protein ID when one sequence failed, else `None` (a problem with the whole input).
- `e.index` (`int | None`): set by `align_files` to the failed sequence's 0-based position in the input; `None` from
  `codon_align`.

| `code` | `seq_id` | When |
|---|---|---|
| `row_lengths` | `None` | protein rows of different lengths (the message says whether the first file looks like DNA) |
| `count_mismatch` | `None` | `len(proteins) != len(dnas)` |
| `no_sequences` | `None` | empty input |
| `partial_id_match` | `None` | `strict=True` and only some protein IDs are DNA IDs |
| `invalid_protein` | set | a character that is not an amino acid, gap, `*` or frameshift digit (e.g. `@`) |
| `empty` | set | a row that is all gaps |
| `inconsistent` | set | the DNA does not spell the protein (no match found, even by 10-residue anchors); with `strict=True` also a gene whose checkable codons mostly mismatch or whose position in a longer DNA is a guess |
| `unknown_table` | `None` | `table` is not a supported NCBI genetic-code number (the message lists them); the `coal` command rejects it earlier, as a usage error |
| `duplicate_id` | `None` | two proteins (or two DNA records) have the same ID (`the protein alignment: duplicate sequence ID 'a': every sequence needs its own ID`); the `coal` command refuses the file with the same wording (as `input_format`) |

Arguments of the wrong type raise `TypeError` with a message naming the argument (`proteins must be a list of
(id, sequence) pairs, not dict`, `proteins[0]: the ID must be a str, not int`, `dnas[1] ('g2'): the sequence must be
a str, not bytes`): `proteins` and `dnas` must be lists or tuples of 2-item `(str, str)` tuples or lists (a dict, a
generator or a bare string is refused, not iterated), never a `KeyError` or `AttributeError` from inside the engine.

`align_files` can also raise `file_error` (an input cannot be opened or the output cannot be written),
`input_format` (a file that is not FASTA or CLUSTAL, or not ASCII text: the exception is `coal.seqio.InputFormatError`,
a subclass of both `CoalError` and `ValueError`), `all_failed` (`keep_going=True` and every sequence failed),
`resume_mismatch` (`resume=True` with changed inputs or options), `worker_died` (a worker process of a `jobs > 1`
run was killed, most likely by the out-of-memory killer) and, with no `code` argument, `error` (an unknown
`out_format`, or `-output codon` with `nogap`/`nomismatch`/`blockonly`).

```python
>>> from coal import CoalError
>>> try:
...     codon_align([("a", "MK")], [("a", "CCCCCC")])
... except CoalError as e:
...     print(e.code, e.seq_id, e)
inconsistent a a: protein and DNA are inconsistent (no match found)
```

## `align_files`

```python
coal.align_files(pep_path, nuc_paths, out_path=None, table=1, nogap=False, nomismatch=False, jobs=0,
                 warn=None, notice=None, note=None, strict=False, out_format="fasta", blockonly=False,
                 rec=None, keep_going=False, progress=False, resume=False) -> int
```

The `coal` command as a function (`coal.pipeline.run`): reads a protein alignment (FASTA or CLUSTAL) and one or more
coding-DNA FASTA files, streams the codon alignment to `out_path` (or to stdout when `out_path` is `None`) in
`out_format`, and returns the number of sequences left out by `keep_going` (0 otherwise). The output is written to a
temporary file first and moved into place only when the run succeeded, so an error never leaves a partial file. Large
inputs are threaded in worker processes (`jobs`): call it from under `if __name__ == "__main__":` in a script, as
with `multiprocessing`. Output bytes, warnings and errors are identical to the command line's.

### Stable arguments

| Argument | Type | Meaning |
|---|---|---|
| `pep_path` | path | the protein alignment (FASTA or CLUSTAL); a pipe or `/dev/stdin` is copied to a temp file first |
| `nuc_paths` | a path, or a `list` / `tuple` of paths | the coding DNA (FASTA), one or more files; a single `str` or `os.PathLike` is one file (not iterated as characters). Anything else (a generator or other iterator, a set, a dict, `None`), or a non-path element, is a `TypeError` naming the argument (`nuc_paths must be a path or a list of paths, not generator`, `nuc_paths[1] must be a path (str or os.PathLike), not int`) before any file is opened; a generator used to be accepted and is not any more: pass `list(...)` |
| `out_path` | path or `None` | the output file; `None` writes to `sys.stdout` |
| `table`, `nogap`, `nomismatch`, `strict` | as for `codon_align` | |
| `out_format` | `str` | `"fasta"` (default), `"clustal"`, `"paml"` or `"codon"` (PAL2NAL's four formats, same bytes) |
| `blockonly` | `bool` | PAL2NAL `-blockonly`: keep only the columns marked `#` in a CLUSTAL input |
| `keep_going` | `bool` | `--keep-going`: leave out sequences that cannot be converted (or, with `strict`, are refused) and write the rest; the return value counts them. If every sequence fails, `CoalError` `all_failed` |
| `jobs` | `int` | worker processes: `0` = up to all cores (small inputs always run in one process), `N` = at most N, `1` = in-process |
| `warn` | callable or `None` | receives each per-codon warning string (`<id> pepAlnPos ...`, the formats above) |
| `notice` | callable or `None` | receives the pairing warnings, the probably-wrong-gene notices and, with `keep_going`, `skipped <id>: <reason>` for each sequence left out |
| `note` | callable or `None` | receives informational lines (the early "wrong genetic code?" suggestion, `--resume` messages). The command prints them after `coal: notice:` |
| `rec` | `coal.report.RunReport` or `None` | collects the run report (below) |

Available but **not** covered by the 1.0 promise: `progress` (percent done on stderr for inputs over 10,000
sequences), `resume` (continue an interrupted run from its checkpoint; needs `out_path`), `min_nt_per_job`,
`fmt` (an alias of `out_format`) and `_resumable`.

### Example

Run `examples/one_broken_gene` (19 good genes and one frameshift), leaving the broken gene out, and read the report:

```python
from coal import align_files, CoalError
from coal.report import RunReport

pep, dna, out = "protein.aln.fasta", ["dna.fasta"], "codon.fasta"
messages = []
rec = RunReport([], pep, dna, out, "fasta", 1, False, False, keep_going=True)
skipped = align_files(pep, dna, out_path=out, keep_going=True, warn=messages.append, notice=messages.append, rec=rec)
rec.exit_code = 2 if skipped else 0
report = rec.to_dict()           # the same dict `coal --report` writes (docs/report-schema.md)
```

gives

```python
>>> skipped
1
>>> messages
['skipped OG0001_sp14: protein and DNA are inconsistent (no match found)']
>>> report["summary"]
{'total': 20, 'ok': 19, 'warning': 0, 'failed': 1, 'pairing': 'id'}
>>> report["sequences"][13]["messages"][0]["code"]
'frameshift'
>>> rec.summary_lines()
['1 sequence of 20 could not be converted and was left out (first: OG0001_sp14: protein and DNA are inconsistent (no match found). Frameshift at column 102: one base is missing in the DNA, so from here on the codons are read out of step.)']
```

Without `keep_going` the same call raises, and nothing is written:

```python
>>> try:
...     align_files(pep, dna, out_path="codon.fasta")
... except CoalError as e:
...     print(e.code, e.seq_id, e.index)
inconsistent OG0001_sp14 13
```

### The run report from Python

`coal.report.RunReport` collects what the `--report` / `--json` options write. The stable parts:

- `RunReport(argv, pep, nucs, out, fmt, table, nogap, nomismatch, keep_going=False, strict=False, blockonly=False)`:
  `argv` is the list of command-line arguments shown in the report's `command_line` (pass `[]`), the rest are the
  same values you give `align_files`.
- `rec.exit_code`: set it before reading the report (0; 1 after a `CoalError`; 2 when `keep_going` skipped
  sequences), as the command line does.
- `rec.fail(code, message, seq_id)`: after catching a `CoalError e`, `rec.fail(e.code, str(e), e.seq_id)` records it,
  so the report says why the run failed (the command line does this).
- `rec.to_dict()`: the report, a JSON-ready dict with schema `coal.report/1` ([docs/report-schema.md](report-schema.md)).
- `rec.write(path)`: `to_dict()` as JSON, written atomically, as `--report` does.
- `rec.summary_lines()`: the plain-language `coal: summary:` lines (a list of strings, empty for a clean run).

## Versioning promise

From 1.0 on COAL follows [semantic versioning](https://semver.org/). Within a major version:

- **Command line:** every documented option of `coal` and `pal2nal` keeps its name and meaning. Options can be
  added; none are removed or renamed.
- **Output bytes:** for the same input and options the alignment written by `coal`, `pal2nal` and `coal-cpp` is
  byte-identical across releases, in all four formats. A change in output bytes is a major release, or a documented
  fix of an unexplained difference from PAL2NAL v14 (listed in the changelog as such).
- **Exit codes:** 0 success; 1 bad input or a gene that cannot be converted; 2 usage error of `coal`, or sequences
  left out with `--keep-going`; 141 closed output pipe; 130 / 143 / 129 interrupted.
- **Report schema:** `coal.report/1` only gains fields; a removed or renamed field, or a changed meaning, bumps the
  schema id to `coal.report/2` (policy in [docs/report-schema.md](report-schema.md)).
- **Python:** `codon_align`, `align_files` (the stable arguments above), `CoalError` and its `code`/`seq_id`/`index`
  attributes, the `code` values listed above, the warning formats, and the parts of `RunReport` named above. New
  keyword arguments and new `code` values can be added in minor releases.
- **Messages on stderr** (warning and notice wording, the summary lines) and the text fields of the report are
  **not** frozen: parse the codes, not the text.

Before 1.0 (0.x) the same rules are followed in spirit, but a minor release may still change any of them; the
changelog says so when it does.
