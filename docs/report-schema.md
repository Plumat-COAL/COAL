# Run report (`coal.report/1`)

`coal --report FILE.json` writes a JSON run report, also when the run fails; `--json` prints the same report to
stdout (needs `-o`). `coal-cpp` writes the identical report (same fields, same values; `tests/test_report_schema.py`
checks both engines on every change), and the COAL app reads it. The machine-readable schema is
[`docs/report.schema.json`](report.schema.json) (JSON Schema draft 2020-12); this page says what each field means.

The report is written from what the engine itself did, not parsed from its messages, so a pipeline can read it
instead of stderr: **read the `code` fields, not the `text` fields**. Every `text` is a sentence for a person and
may be reworded in any release.

## Stability policy

- The schema id is `coal.report/1`. Within it, **fields are only added**: a new top-level field, a new key in an
  object, a new message or hint `code`, a new value of an open string field. Readers must ignore what they do not
  know.
- **Removing or renaming a field, changing a field's type or meaning, or removing a documented `code` or enum
  value bumps the id to `coal.report/2`.** A reader should check `schema` first.
- The sets marked *closed* below (statuses, pairing modes, output formats, hint codes, `exit_code`) do not change
  within `/1`. Message codes and error codes are *open*: new ones may appear; the ones listed here keep their meaning.
- The example at the end of this page is the report of `examples/one_broken_gene` run with `--keep-going`, as the
  engine writes it today; it is also `tests/report_example.json`.

## Top-level fields

All of these are present in every report, in this order, except `error` and `resumed_from`.

| Field | Type | Meaning |
|---|---|---|
| `schema` | string | always `"coal.report/1"` |
| `coal_version` | string | the engine's version, as `coal --version` prints it (e.g. `"0.4.2"`) |
| `source` | string | always `"coal"`: the engine wrote the report (the COAL app writes `coal_app ...` here when it has to build a report from an engine without `--report`; such a report is not covered by this page). There is no separate engine field: `coal` and `coal-cpp` write the same report |
| `command` | string | a command that reproduces the alignment: `coal`, the input paths, the options that change the output (`-output`, `-codontable`, `-nogap`, `-nomismatch`, `-blockonly`, `--strict`, `--keep-going`) and `-o`; shell-quoted. `--report`, `-j`, `-q` and the like are left out |
| `command_line` | string | the command line as given, shell-quoted, starting with `coal` |
| `pal2nal_equivalent` | string | the PAL2NAL v14 command line with the same meaning: `pal2nal.pl PEP NUC -output FMT [-codontable N] [-nogap] [-nomismatch] [-blockonly]` |
| `options` | object | the options in effect: `output_format` (`fasta`, `clustal`, `paml`, `codon`; closed), `codontable` (integer, one of the 17 supported NCBI codes), and the booleans `nogap`, `nomismatch`, `keep_going`, `strict`, `blockonly` |
| `inputs` | array of objects | the protein alignment first (`role` = `"protein"`), then each DNA file (`"dna"`) in command-line order. Each: `path` (string, as given), `sha256` (64 hex characters, or `null` if the file could not be read), `sequences` (integer, or `null` if the run failed before the file was read, e.g. the DNA files when the protein alignment is malformed) |
| `output` | object or `null` | `null` when no alignment was written (the run failed). Else `path` (string, or `null` when the alignment went to stdout), `format` (as `options.output_format`), `sha256` of the written file (string; `null` only if it could not be hashed), `sequences` (integer: rows written, i.e. the total minus the sequences `--keep-going` left out) |
| `summary` | object | counts, see below |
| `hints` | array of objects | run-level advice, see [Hints](#hints); `[]` when there is none |
| `errors` | array of objects | run-level errors (`{"code", "text"}`): the ones that are not about one sequence (e.g. `count_mismatch`, `file_error`). A failure of one sequence is in that sequence's `messages` and in `error`, not here. Empty in a successful run |
| `sequences` | array of objects | one entry per protein in alignment order, see [Sequence entries](#sequence-entries) |
| `runtime_seconds` | number | wall-clock time of the run, 3 decimals |
| `exit_code` | integer | the process exit status as the report is written: `0`, `1` (error) or `2` (`--keep-going` left sequences out); closed. The signal exits (130, 141, 143, 129) write no report |
| `summary_text` | array of strings | the plain-language `coal: summary:` lines, as printed on stderr; `[]` for a clean run |
| `methods_text` | string | a ready-to-paste methods sentence for a paper (version, method, genetic code, what happened to how many sequences, the command); `""` when the run did not complete (no output, or an error other than `--keep-going` skips). A run where fewer than half the sequences matched exactly starts with `Fix the warnings first: ...` |
| `error` | object | **only when the run failed** (`exit_code` 1): `code` (the error code, see below), `text` (the message printed after `coal: error:`), `detail` (currently the same as `text`) |
| `resumed_from` | integer | **only when `--resume` continued a run**: the number of sequences the interrupted run had already finished |

### `summary`

| Field | Type | Meaning |
|---|---|---|
| `total` | integer | proteins in the alignment (`len(sequences)`) |
| `ok`, `warning`, `failed` | integer | sequences with that `status` |
| `not_checked` | integer | **only when non-zero**: sequences with status `not_checked`. In a run that stopped at an error this counts every sequence except the failed one, i.e. also the ones that had been checked before it (they are `not_written`); the `summary_text` line "N sequences were not checked" counts only the ones never reached |
| `pairing` | string or `null` | how proteins were paired with DNA records: `"id"` (every protein ID is a DNA ID), `"normalised"` (some IDs matched only after removing a version suffix, case or database fields, where file order would have mispaired them), `"order"` (by position in the files, as PAL2NAL), or `null` when the run failed before pairing (count mismatch, unreadable file). Closed |
| `columns_removed` | integer | **only with `-nomismatch` and an output written**: alignment columns removed because a codon did not match its amino acid (or held a B/Z/J/O/U residue) |
| `columns_removed_ambiguous` | integer | **with `columns_removed`**: how many of those columns were removed only because of codons with N or other ambiguous bases (no definite mismatch) |

## Sequence entries

One object per protein, in alignment order:

| Field | Type | Meaning |
|---|---|---|
| `id` | string or `null` | the protein ID (the output row's name); `null` only when the protein file could not be read at all |
| `index` | integer | 0-based position in the alignment |
| `dna_id` | string or `null` | the ID of the DNA record it was paired with; `null` before pairing (and for a sequence that was never paired) |
| `pairing` | string or `null` | this sequence's pairing: `"id"`, `"normalised"`, `"order"` or `null` (as `summary.pairing`; in a `normalised` run each entry says whether its own ID matched exactly or only after normalising) |
| `status` | string | `ok`, `warning`, `failed` or `not_checked` (closed), see below |
| `messages` | array of objects | what happened, see [Messages](#messages); `[]` for `ok` |

Statuses:

- `ok`: the DNA spells the protein exactly; the row was written.
- `warning`: the row was written and is usable, but something is worth knowing: a mismatching codon (`mismatch`),
  a non-standard amino acid (`unknown_aa`), another notice (`warning`, `no_anchor`), a pairing by file order
  with differing IDs but matching codons (`suspect_pairing` alone).
- `failed`: either **not written** (a `CoalError` for this sequence: `inconsistent`, `invalid_protein`, `empty`;
  the run stopped there, or `--keep-going` left it out) or **written but probably wrong** (`suspect_unrelated`,
  `suspect_position`, or `suspect_pairing` together with a `mismatch`): the first message's code tells which. A
  written-but-flagged row should be removed before analysis; `--strict` refuses such genes instead. A sequence
  refused by `--strict` carries the `inconsistent` failure message.
- `not_checked`: the run stopped at an error before this sequence was finished with. In a stopped run every
  sequence but the failed one is `not_checked`: the ones before it carry `not_written` plus their own messages (they
  were checked but no output exists), the ones after it carry `not_reached`; when the run failed before pairing, all
  are `not_checked` with no messages.

## Messages

Every message has a `code`. All but `not_written` / `not_reached` have a `text`. The codes the engine emits (open
set; new codes may be added):

| `code` | Extra fields | Status | Meaning |
|---|---|---|---|
| `mismatch` | `column` (1-based alignment column), `residue`, `codon`; `removed: true` with `-nomismatch`; `start_codon: true` when the sequence's only mismatches are at its first Met | warning | the codon at that column does not encode the amino acid under the genetic code used. A codon with N or other ambiguous bases is also a `mismatch`; such entries count as "could not be checked", not as definite mismatches |
| `unknown_aa` | `column`, `residue` | warning | `B`, `Z`, `J`, `O` or `U` in the protein: its codon was not checked (any codon accepted) |
| `warning` | - | warning | another per-sequence notice, verbatim after `COAL warning: ` (today: a protein of at most 3 residues that fits in several places of the DNA; the leftmost was used, as PAL2NAL) |
| `suspect_unrelated` | - | failed (written) | more than half of the checkable codons mismatch: almost certainly the wrong DNA or genetic code |
| `suspect_position` | - | failed (written) | no 10-residue stretch of the protein was found in the DNA and fewer than 10 codons could be checked, in a DNA longer than the gene: the row may be shifted |
| `no_anchor` | - | warning | no 10-residue stretch was found but at most half of the checkable codons mismatch: codons assigned by position (the text says how many mismatch and how many are ambiguous) |
| `suspect_pairing` | `dna_id` | warning, or failed with a `mismatch` | paired with this DNA record only by file order, in a run where only some protein IDs are DNA IDs |
| `possible_frameshift` | `column`, `end_column` | with `mismatch` | the definite mismatches form one dense run from `column` to `end_column`: looks like an indel in the DNA, not a genetic-code problem |
| `frameshift` | `column`, `shift` (-1, 1, -2, 2, -4, 4, -5, 5: bases missing (<0) or extra (>0)) | failed (not written) | where a failed gene breaks: from this column the codons read out of step |
| `missing_dna` | `column` | failed (not written) | the DNA ends before the protein does |
| `inconsistent` | `column` (the break location) or none (the error itself) | failed (not written) | with `column`: the DNA does not spell the start of the protein in any frame, or from this column on no longer spells it. Without: the engine's error (`protein and DNA are inconsistent (...)`) |
| `invalid_protein`, `empty` | - | failed (not written) | the engine's error for this sequence (an invalid character; a row of only gaps) |
| `not_written` | `failed_id` | not_checked | checked before the run stopped at sequence `failed_id`; no output exists |
| `not_reached` | `failed_id` | not_checked | never checked: the run stopped at `failed_id` before it |

Order within `messages`: for a failed gene the break location (`frameshift`, `missing_dna`, `inconsistent` with
`column`) comes first, then the error; for a flagged gene `suspect_*` first, then `suspect_pairing`, then the
per-codon messages in column order; `possible_frameshift` is appended last.

## Hints

Run-level advice, each `{"code", "text", ...}`, in this order (closed set of codes):

| `code` | Extra fields | When |
|---|---|---|
| `partial_id_match` | - | only some protein IDs are DNA IDs, so sequences were paired by order (`summary.pairing` = `order`) |
| `normalised_ids` | - | `summary.pairing` = `normalised`: which IDs matched only after tidying, and what file order would have done |
| `id_mismatch` | - | pairing by order (or none) with protein IDs that are not DNA IDs; says whether they differ only by a version suffix and whether the order pairing looked right (every pair matched) |
| `genetic_code` | `suggested_codontable` (integer), `checked` (sequences tested), `fit` (how many match under the suggested code), `fit_current` (how many match under the code used) | another genetic code makes more sequences match their protein and fixes at least one with problems. At most 50 sequences are tested, problem genes first; `text` ends `Re-run with code N?` |

## Error codes

`error.code`, `errors[].code` and the failure message of a sequence use the same codes as `CoalError.code` in the
Python API ([docs/python-api.md](python-api.md#errors)): per sequence `inconsistent`, `invalid_protein`, `empty`;
for the run `count_mismatch`, `row_lengths`, `no_sequences`, `partial_id_match` (`--strict`), `all_failed`
(`--keep-going` and every sequence failed), `file_error` (an input cannot be opened, the output cannot be written),
`input_format` (not FASTA / CLUSTAL, undecodable bytes), `resume_mismatch` (`--resume` with changed inputs or
options) and `error` (an invalid option combination). The code `stopped` is reserved for the COAL app, which adds it
to `errors` when it finds a run was interrupted; the engine never writes it (an interrupted run writes no report).
`unknown_table` (Python API only: an unsupported `table`) never reaches a report written by the command line, which
rejects the number as a usage error before the run starts. `worker_died` is written by the Python engine only (a
worker process of a `-j` run was killed, most likely by the out-of-memory killer; `coal-cpp` uses threads, so it has
no such error).

## Validating and regenerating

```
python -m unittest tests.test_report_schema          # both engines' reports against docs/report.schema.json
pip install jsonschema && python -c "import json, jsonschema; jsonschema.validate(json.load(open('run.json')), json.load(open('docs/report.schema.json')))"
```

To refresh the example below and `tests/report_example.json` after a change to the report:

```
cd examples/one_broken_gene && coal protein.aln.fasta dna.fasta --keep-going -o codon.fasta --report report.json
```

(then set `runtime_seconds` to `0.044` so the two copies stay equal, and paste the file here).

## Complete example

`examples/one_broken_gene` (20 genes, one with a deleted base) with `--keep-going`: 19 rows written, the broken
gene left out, exit status 2.

```json
{
 "schema": "coal.report/1",
 "coal_version": "0.4.2",
 "source": "coal",
 "command": "coal protein.aln.fasta dna.fasta --keep-going -o codon.fasta",
 "command_line": "coal protein.aln.fasta dna.fasta --keep-going -o codon.fasta --report report.json",
 "pal2nal_equivalent": "pal2nal.pl protein.aln.fasta dna.fasta -output fasta",
 "options": {
  "output_format": "fasta",
  "codontable": 1,
  "nogap": false,
  "nomismatch": false,
  "keep_going": true,
  "strict": false,
  "blockonly": false
 },
 "inputs": [
  {
   "role": "protein",
   "path": "protein.aln.fasta",
   "sha256": "f3fd705c35aad60a8e761f622067f209f89d98c8e4a9c235a8d2ffb9de0aebfb",
   "sequences": 20
  },
  {
   "role": "dna",
   "path": "dna.fasta",
   "sha256": "19bd9e5bd96d60e0fa3435984b46ed7b8c31c6d7a71fcc0065476ae6ac8032bc",
   "sequences": 20
  }
 ],
 "output": {
  "path": "codon.fasta",
  "format": "fasta",
  "sha256": "a0093d9c754e40e961eb39bcc00660fa1cbff3ee3b3bd62927256635fe59514c",
  "sequences": 19
 },
 "summary": {
  "total": 20,
  "ok": 19,
  "warning": 0,
  "failed": 1,
  "pairing": "id"
 },
 "hints": [],
 "errors": [],
 "sequences": [
  {
   "id": "OG0001_sp01",
   "index": 0,
   "dna_id": "OG0001_sp01",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp02",
   "index": 1,
   "dna_id": "OG0001_sp02",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp03",
   "index": 2,
   "dna_id": "OG0001_sp03",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp04",
   "index": 3,
   "dna_id": "OG0001_sp04",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp05",
   "index": 4,
   "dna_id": "OG0001_sp05",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp06",
   "index": 5,
   "dna_id": "OG0001_sp06",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp07",
   "index": 6,
   "dna_id": "OG0001_sp07",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp08",
   "index": 7,
   "dna_id": "OG0001_sp08",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp09",
   "index": 8,
   "dna_id": "OG0001_sp09",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp10",
   "index": 9,
   "dna_id": "OG0001_sp10",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp11",
   "index": 10,
   "dna_id": "OG0001_sp11",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp12",
   "index": 11,
   "dna_id": "OG0001_sp12",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp13",
   "index": 12,
   "dna_id": "OG0001_sp13",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp14",
   "index": 13,
   "dna_id": "OG0001_sp14",
   "pairing": "id",
   "status": "failed",
   "messages": [
    {
     "code": "frameshift",
     "column": 102,
     "shift": -1,
     "text": "Frameshift at column 102: one base is missing in the DNA, so from here on the codons are read out of step."
    },
    {
     "code": "inconsistent",
     "text": "protein and DNA are inconsistent (no match found)"
    }
   ]
  },
  {
   "id": "OG0001_sp15",
   "index": 14,
   "dna_id": "OG0001_sp15",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp16",
   "index": 15,
   "dna_id": "OG0001_sp16",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp17",
   "index": 16,
   "dna_id": "OG0001_sp17",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp18",
   "index": 17,
   "dna_id": "OG0001_sp18",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp19",
   "index": 18,
   "dna_id": "OG0001_sp19",
   "pairing": "id",
   "status": "ok",
   "messages": []
  },
  {
   "id": "OG0001_sp20",
   "index": 19,
   "dna_id": "OG0001_sp20",
   "pairing": "id",
   "status": "ok",
   "messages": []
  }
 ],
 "runtime_seconds": 0.044,
 "exit_code": 2,
 "summary_text": [
  "1 sequence of 20 could not be converted and was left out (first: OG0001_sp14: protein and DNA are inconsistent (no match found). Frameshift at column 102: one base is missing in the DNA, so from here on the codons are read out of step.)"
 ],
 "methods_text": "Codon alignments were built with COAL v0.4.2 (https://github.com/plumat-coal/COAL), which places each coding DNA sequence onto the protein multiple sequence alignment codon by codon (the method of PAL2NAL; Suyama et al. 2006), using the standard genetic code (NCBI translation table 1). Of 20 sequences, 19 matched their protein exactly, and 1 could not be aligned and was excluded. Command: coal protein.aln.fasta dna.fasta --keep-going -o codon.fasta"
}
```
