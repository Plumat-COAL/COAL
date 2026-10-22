# Using COAL

- **Input:** a protein alignment (FASTA or CLUSTAL) and the unaligned coding DNA (FASTA, one or more files).
  Non-letters in the DNA (gaps, digits, spaces) are ignored, as in PAL2NAL.
- **Output:** a codon alignment, using the protein-alignment IDs, in any of PAL2NAL v14's formats, as PAL2NAL writes
  them (byte for byte, apart from the documented differences): `-output fasta` (COAL's default; 60 columns per line),
  `clustal` (PAL2NAL's default), `paml` (codeml's input) or `codon` (amino acids above their codons; not combinable
  with `-nogap`/`-nomismatch`/`-blockonly`, as in PAL2NAL). The drop-in `pal2nal` command defaults to clustal, like
  PAL2NAL.
- **Pairing:** sequences are matched by ID if every protein ID is in the DNA file(s), otherwise by order, as in
  PAL2NAL; except that where the order would give a protein the DNA of another gene whose ID matches it after
  normalising (version suffix such as `.1`, upper/lower case, `sp|…|` database fields; one-to-one only), COAL pairs
  by those IDs and says so (a deliberate difference, below).
- **Extra sequence around the gene** (UTRs, stop codon, poly-A) is skipped: COAL finds the leftmost place in the
  DNA, in any reading frame, where the whole protein fits.
- **Genetic codes:** `-codontable N` for the 17 NCBI codes PAL2NAL supports (1–6, 9–16, 21–23). The first Met may use
  any start codon of that code (also when the gene doesn't fit exactly and COAL needs its fallback placement, unlike
  PAL2NAL v14). When the genes fit another code better, COAL says so once, before the per-codon warnings.
- **Frameshifts:** a digit N in the protein alignment means N bases (PAL2NAL convention). That column is widened
  to the next multiple of 3 for every sequence. tfastx/tfasty notation is read as PAL2NAL reads it: `\` = 1 base,
  `/` followed by (gaps and) a residue = a gap, then 2 bases.
- **Mismatches:** if the protein doesn't fit exactly, COAL uses PAL2NAL's fallback: it lines up 10-residue
  chunks and prints a warning for each codon that doesn't match its residue (at most 20 by default,
  `--max-warnings N`, 0 = all; the `pal2nal` drop-in prints them all, like PAL2NAL). `-nomismatch` removes those
  columns (and, as in PAL2NAL, then prints no mismatch warnings). If even that fails, it's an error.
- **Probably-wrong data:** like PAL2NAL, COAL does not refuse a gene for having many mismatches, so its output is
  the same as PAL2NAL's. It adds a notice on stderr naming the gene when (1) more than half of a gene's checkable
  codons (made only of A/C/G/T; codons with N or other ambiguity codes don't count) don't match: probably the wrong
  DNA or the wrong genetic code; (2) no 10-residue chunk matched, fewer than 10 codons are checkable and the DNA is
  longer than the gene: the start can't be told (e.g. a masked gene with a 5' UTR), so the row may be shifted; or
  (3) only some protein IDs are DNA IDs, so sequences are paired by order. **`--strict`** (or `-strict`) turns these
  three cases into errors (exit 1, nothing written), for pipelines that would rather stop than continue.
- **`-nogap`** removes codon columns that contain a gap or a stop codon of the selected genetic code (any case).
- **`-blockonly`** keeps only the columns marked `#` under the blocks of a CLUSTAL protein alignment (PAL2NAL's
  convention, e.g. Gblocks-selected positions); without `-blockonly` those marks are shown under each block of
  clustal and codon output, as PAL2NAL does. No effect if the input has no `#` marks.
- **Messages and exit codes:** the exit code is non-zero on any error (PAL2NAL always exits 0). Warnings and notices
  go to stderr; `-q` (like PAL2NAL's `-nostderr`) silences all of them, including the notices about probably-wrong
  data, so check the exit code or use `--strict` when running quietly. Errors are always printed. stdout is the same
  with or without `-q`. Every stderr line starts with `coal: error:` (the run or a gene failed; the message says
  which file, gene or position and what to do), `coal: warning:` (data that may be wrong: mismatches,
  probably-wrong genes, pairing), `coal: notice:` (information, e.g. about `--resume`) or `coal: summary:` (the
  end-of-run summary; it does not repeat the error).
- **Closed output pipe** (e.g. `coal pep.aln nuc.fa | head`): COAL stops silently with exit status **141**
  (128 + SIGPIPE, the Unix convention), also on Windows, so a pipeline can tell "the reader stopped early" from a
  real error (1). Other exit codes: 130 Ctrl-C, 143 SIGTERM, 129 SIGHUP; a usage error is 2 for `coal` and 1 for the
  drop-in `pal2nal`, which reads its command line as PAL2NAL v14 does (a dangling `-output` or `-codontable` is
  ignored, `-codontable 2x` means code 2).

**Not supported:** PAL2NAL's `-html` web-server mode (a clear error) and Gblocks-format protein input.

## Large inputs

COAL streams the data. It reads the protein alignment and indexes the DNA file(s) (IDs and file positions only),
then handles one sequence at a time, fetching its DNA from disk and writing its output row, so memory doesn't grow
with the size of the sequences. By default the Python engine uses all CPU cores for large inputs and one process for
small ones (each extra process needs about 4 million nucleotides of DNA to pay for its start-up); `coal-cpp` uses
threads. `-j N` caps the number of workers and `-j 1` forces one. Output is always in the input order. Output is
written to a temporary file and only moved into place if the whole run succeeds, so an error never leaves a partial
file.

- `-nogap`/`-nomismatch`/`-blockonly` need every row before deciding which columns to remove, so the rows are first
  written to a temporary file (disk space about the size of the output). The clustal and codon formats interleave
  all sequences block by block, so they use temporary files the same way (about twice the output size on disk) and
  hold at most 16 MB of rows in memory at a time (measured with 500,000 sequences of 150 codons, `-j 1`: peak memory
  122 MB for fasta, 139 MB for clustal, 145 MB for codon). `-j` gives the same output in every format.
- CLUSTAL input is read into memory (FASTA input is streamed).
- Memory still holds the sequence IDs and their file positions, about 300 bytes per sequence (≈150 MB per
  500,000 sequences).

## Run report and long runs

- `--report FILE.json` writes a run report (schema `coal.report/1`, the one the GUI reads), also when the run fails:
  every sequence's status (`ok`, `warning`, `failed` = written but probably wrong, or left out; `not_checked`), its
  messages, the pairing, a genetic-code suggestion and a ready-to-paste methods sentence. `--json` prints it to
  stdout (needs `-o`). Every field: [report-schema.md](report-schema.md) (JSON Schema: `docs/report.schema.json`). At the end of a run with problems, `coal: summary:` lines on stderr sum them up (nothing about
  mismatches with `-nomismatch`, which like PAL2NAL prints none; `-q` silences them).
- `--keep-going` skips sequences that cannot be converted (or, with `--strict`, are refused) instead of stopping:
  the others are written in order, each skipped one is named on stderr, exit status 2. A frameshift digit of a
  skipped row still widens its column for the others.
- `--max-warnings N` prints at most N per-codon warnings (default 20, 0 = all), then a count (notices and errors are
  never capped).
- `--progress` shows percent done and time left on stderr for inputs over 10,000 sequences.
- `--resume` continues an interrupted `-o` run from `<out>.coal-checkpoint.json` (written at most every 2 s next to
  the real output file); the output is byte-identical to an uninterrupted run. If no checkpoint was written yet, an
  interrupted run leaves nothing behind. Not for piped inputs or an output that is not a regular file (COAL says so
  and runs normally).
