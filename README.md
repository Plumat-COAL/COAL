# COAL

**Codon alignments from a protein alignment and its coding DNA.** COAL is a fast, open-source alternative to
[PAL2NAL](https://www.bork.embl.de/pal2nal/): it maps unaligned coding DNA onto an existing protein multiple sequence
alignment to produce a codon alignment, ready for dN/dS (PAML codeml, HyPhy) or DNA phylogenies.

Built by Chloe Newton and Keanu Murati at the Reinvent the Wheel AI × Science Tool-building Hackathon, London,
October 2026.

## Download and install

### Option 1: download the app from the website (recommended)

**https://plumat-coal.github.io/COAL/** has the Windows installer and step-by-step install instructions for
Windows, macOS and Linux. No terminal is needed on Windows.

- **Windows:** click **Download Windows installer** on the website
  ([direct link to the latest installer](https://github.com/plumat-coal/COAL/releases/latest/download/COAL-Setup-windows-x64.exe);
  all files, including a portable zip and `SHA256SUMS.txt`, are on the
  [Releases page](https://github.com/plumat-coal/COAL/releases/latest)). The app is not code-signed yet, so Windows
  shows **"Windows protected your PC"**: click **More info**, then **Run anyway**. The
  [install page](https://plumat-coal.github.io/COAL/install.html) walks through every screen, including what to do
  if Smart App Control blocks it. User guide: `desktop/README.md`.
- **macOS and Linux:** follow the [install page](https://plumat-coal.github.io/COAL/install.html#macos) (three
  commands, also below).

### Option 2: macOS and Linux (or Windows with Python)

Needs Python 3.10+:

```
python3 -m venv ~/coal-env
~/coal-env/bin/python -m pip install "coal[gui] @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_src.tar.gz"
~/coal-env/bin/coal-gui          # the app, in your browser at http://127.0.0.1:7860
~/coal-env/bin/coal --help       # the command line
```

(Install from the release link above, not `pip install coal`: that name on PyPI belongs to another project.)

### Option 3: from source

```
git clone https://github.com/plumat-coal/COAL && cd COAL
pip install -e .            # command line;  .[gui] for the app, .[mcp] for AI agents, .[fast] for numpy
```

## Why COAL

- **Fast:** up to **~119× faster** than PAL2NAL on the benchmarked datasets (see [Benchmark results](#benchmark-results)),
  with a Python engine and an optional C++ engine, `coal-cpp`, that gives byte-identical output.
- **Checked against PAL2NAL itself:** COAL's output is identical to PAL2NAL v14's except for a short list of documented, deliberate differences (newer NCBI genetic-code tables for codes 3, 10, 11 and 13; `-nogap` removing the selected code's stop codons; the selected code's start codons for the first Met; pairing by IDs that match only after normalising where PAL2NAL v14's file order would give a protein another gene's DNA; and inputs where PAL2NAL v14 writes malformed output or fails). Tested: all 285 stored reference checks pass (238 byte-identical, 4 of them with deliberately different warnings; 24 where both tools fail; 23 deliberate output or exit-status differences: 6 deliberate errors, 3 different outputs, 14 `-nogap` checks with a non-standard code), and in 20,000 random comparisons with the real PAL2NAL v14 (5,000 cases x 4 output formats) 13,826 were byte-identical, 4,681 both failed and 1,493 were deliberate differences, each checked to be exactly the documented change; 0 unexplained differences.
- **Tells you what went wrong:** a per-gene health check, warnings that name the gene and position, genetic-code
  suggestions and a ready-to-paste methods sentence.
- **Drop-in:** the `pal2nal` command accepts an unchanged PAL2NAL v14 command line.

## Quick start

```
coal pep.aln nuc.fasta > codon_aln.fasta            # or: python -m coal ...
coal pep.aln nuc.fasta -codontable 2 -nogap         # PAL2NAL-style options work too
coal big_pep.fasta big_nuc.fasta -o out.fasta       # large data: uses all cores automatically
coal pep.aln nuc.fasta -o out.fasta --report run.json --keep-going    # per-gene report; skip genes that fail
pal2nal pep.aln nuc.fasta -output paml -nogap       # an unchanged PAL2NAL v14 command line
```

Example inputs to try are in `examples/` (`clean`, `mito`, `messy_ids`, `one_broken_gene`).

## Ways to use COAL

- **Windows app:** installer from the [website](https://plumat-coal.github.io/COAL/); a local GUI, nothing is uploaded.
- **Web GUI on any OS:** `coal-gui` (or `python -m coal_app`) runs on this computer at http://127.0.0.1:7860, with a
  per-gene health check, an alignment view that marks where a gene breaks, genetic-code suggestions and a
  ready-to-paste methods sentence.
- **Command line:** `coal` / `python -m coal`.
- **Drop-in for PAL2NAL pipelines:** `pal2nal pep.aln nuc.fasta [PAL2NAL options]` (or `python -m coal_compat`)
  accepts an unchanged PAL2NAL v14 command line, with PAL2NAL's default clustal output, all four formats,
  `-blockonly` and `-nostderr` (= `-q`); `-html` is a clear error.
- **AI agents (MCP):** `pip install -e .[mcp]`, then `python -m coal_app --mcp` (two tools: `codon_align`,
  `explain_report`); the GUI also serves MCP. Agent instructions: `.claude/skills/coal/SKILL.md`.
- **C++ engine for the largest jobs:** `coal-cpp`, built with CMake, same options and output
  ([below](#c-engine-coal-cpp)).
- **Two-species preset:** `demo/diamond_pipeline.py` (DIAMOND reciprocal best hits, then COAL on every pair).
- **Hosted demo:** a Hugging Face Space is prepared but not deployed yet: `space/README.md`.

## How COAL works

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
  PAL2NAL v14).
- **Frameshifts:** a digit N in the protein alignment means N bases (PAL2NAL convention). That column is widened
  to the next multiple of 3 for every sequence. tfastx/tfasty notation is read as PAL2NAL reads it: `\` = 1 base,
  `/` followed by (gaps and) a residue = a gap, then 2 bases.
- **Mismatches:** if the protein doesn't fit exactly, COAL uses PAL2NAL's fallback: it lines up 10-residue
  chunks and prints a warning for each codon that doesn't match its residue. `-nomismatch` removes those
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

### Large inputs

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

### Run report and long runs

- `--report FILE.json` writes a run report (schema `coal.report/1`, the one the GUI reads), also when the run fails:
  every sequence's status (`ok`, `warning`, `failed` = written but probably wrong, or left out; `not_checked`), its
  messages, the pairing, a genetic-code suggestion and a ready-to-paste methods sentence. `--json` prints it to
  stdout (needs `-o`). At the end of a run with problems, `coal: summary:` lines on stderr sum them up (nothing about
  mismatches with `-nomismatch`, which like PAL2NAL prints none; `-q` silences them).
- `--keep-going` skips sequences that cannot be converted (or, with `--strict`, are refused) instead of stopping:
  the others are written in order, each skipped one is named on stderr, exit status 2. A frameshift digit of a
  skipped row still widens its column for the others.
- `--max-warnings N` prints at most N per-codon warnings, then a count (notices and errors are never capped).
- `--progress` shows percent done and time left on stderr for inputs over 10,000 sequences.
- `--resume` continues an interrupted `-o` run from `<out>.coal-checkpoint.json` (written at most every 2 s next to
  the real output file); the output is byte-identical to an uninterrupted run. If no checkpoint was written yet, an
  interrupted run leaves nothing behind. Not for piped inputs or an output that is not a regular file (COAL says so
  and runs normally).

## Deliberate differences from PAL2NAL

Each difference has a reference point other than PAL2NAL v14 itself:
**NCBI** (the current genetic-code tables: codes 3, 10, 11 and 13; the selected code's start codons for the first Met
in the fallback; `-nogap` removing the selected code's stop codons), **the sequence formats' definitions** (`.` as a
gap, IUPAC ambiguity codes, `U` = selenocysteine, RNA read as DNA, FASTA with a BOM, a blank first line or `> id`),
**PAL2NAL's own consistency** (its fallback treats the first Met, RNA and B/Z/J/O differently from its exact-fit
path; rows of different lengths; duplicate IDs) or **policy**, a choice of COAL's and not a PAL2NAL error (malformed
input is an error instead of output, `coal` writes FASTA by default, COAL's own wording of warnings). Some inputs make
PAL2NAL v14 fail with an ERROR (a BOM, a blank first line), others make it print a warning (codes 10, 11, 13):
not every difference is a silent one. Details and evidence in `tests/README.md`.

- An invalid character in the protein alignment (e.g. `@`) is an error. PAL2NAL treats it as X but its output can
  lose a codon silently (fixture `invalid_chars`: `ATGAAAGTTZZZ`, the final CTG lost). PAL2NAL's tfastx/tfasty
  frameshift marks `\` (1 base) and `/` (before a residue: 2 bases) are not invalid: COAL reads them exactly as
  PAL2NAL does (e.g. `MK\VL` → `ATGAAAG--GTTCTG`).
- Genetic codes follow the current NCBI tables. In code 10, TGA = Cys (PAL2NAL treats it as a stop). In codes 11
  and 13, COAL also accepts NCBI's alternative start codons (ATT/ATC/ATA, and ATA, respectively). Code 3 (yeast
  mitochondrial) also accepts GTG as a start codon (PAL2NAL v14's table 3 accepts only ATA/ATG, measured).
- `.` is treated as a gap. In PAL2NAL it misaligns the output row.
- Ambiguous DNA bases (IUPAC codes, including plain N) are accepted when every base they could stand for gives the
  right amino acid, e.g. `NTG` as the first Met in table 4 (all of ATG/CTG/GTG/TTG are table-4 starts). PAL2NAL
  accepts only some fixed combinations, and otherwise falls back to another, frame-shifted placement.
- Probably-wrong data (see above) gives the same output as PAL2NAL plus a notice on stderr (PAL2NAL prints only its
  per-codon warnings, or nothing for a partial ID match). When no ID matches COAL stays silent, like PAL2NAL.
  With `--strict` (COAL's own option) these cases are errors.
- `U` (selenocysteine) in the protein is read as any codon, like B/Z/J/O, with a warning. PAL2NAL's search pattern
  gives `U` no codon at all, so its row is misaligned or the gene fails. Repro: protein `M-UVE`, DNA `ATGTCAGTTGAG`:
  COAL `ATG---TCAGTTGAG`; PAL2NAL fails ("inconsistency").
- With `-nomismatch`, the warnings for B/Z/J/O/U (`U (selenocysteine) is not checked; any codon is accepted`;
  PAL2NAL words it `unknown AA type. Taken as 'X'`) are still printed; PAL2NAL prints nothing there. Mismatch
  warnings are not printed, as in PAL2NAL. The output is the same: both remove those columns.
- B/Z/J/O in a gene that needs the fallback: PAL2NAL also checks their codons against the wrong pattern (B against
  its start-codon list) and prints spurious `does not correspond` warnings; COAL accepts any codon there (stderr only;
  with `-nomismatch` both remove those columns).
- The first Met in the fallback (a gene that doesn't fit exactly, e.g. one N or one mismatching codon): COAL
  accepts the selected code's start and Met codons, as it does when the gene fits exactly (NCBI tables).
  PAL2NAL v14 accepts only ATG/CTG/GTG there, whatever `-codontable` says, so a code-2 gene starting with ATT or
  ATA gets a spurious mismatch at column 1, and `-nomismatch` deletes column 1 for every sequence. Repro: a code-2
  gene starting `ATT` with one `N` further on, `-codontable 2`: PAL2NAL "M does not correspond to ATT"; COAL no
  warning (the same as without the `N`).
- RNA input (U instead of T): COAL reads U as T everywhere. PAL2NAL does too, except in its fallback (used when the
  protein doesn't fit exactly), where a first Met on AUG/CUG/GUG is not accepted: PAL2NAL then flags that codon
  (and `-nomismatch` removes the column) or places the gene elsewhere in the DNA, even in another reading frame.
- Protein rows of different lengths are an error (`protein alignment rows have different lengths`); PAL2NAL silently
  pads or truncates.
- Input that PAL2NAL misreads is read correctly: a protein FASTA starting with a blank line (PAL2NAL takes it for
  CLUSTAL and fails), a UTF-8 BOM (PAL2NAL fails), `> id` with a space after `>` (PAL2NAL finds no protein
  sequences: "number of input seqs differ (aa: 0 ...)", no output).
- Input that PAL2NAL silently mangles is a clear error: duplicate protein or DNA IDs (PAL2NAL concatenates or drops
  records), an empty DNA `>` header, a second DNA file whose first record has no header (PAL2NAL appends it to the
  previous file's last record), and `_` as a stop in the protein (PAL2NAL's alias; COAL: invalid character).
- IDs that match only after normalising (e.g. protein `AT1G01010.1`, DNA `AT1G01010`) in another order: PAL2NAL pairs
  by file order, so proteins get the DNA of other genes; COAL pairs by the normalised IDs (with a notice), because a
  codon alignment built from another gene's DNA is wrong by definition. Where the file order already pairs them
  correctly, COAL's pairing, output and messages are PAL2NAL's. `--strict` still refuses a partial exact-ID match.
  Example: `examples/messy_ids`.
- `coal`'s own default output format is FASTA (PAL2NAL's is clustal); the drop-in `pal2nal` command uses clustal.
- `-nogap` removes the stop codons **of the selected genetic code**, in upper or lower case. PAL2NAL always removes
  the standard code's TAA/TAG/TGA, upper case only, whatever `-codontable` says: with `-codontable 2` (vertebrate
  mitochondria) it deletes TGA columns, which code for Trp there, and keeps AGA/AGG stops; lower-case `taa` is never
  removed. Repro: protein `MW*`, DNA `ATGTGAAGA`, `-codontable 2 -nogap`: PAL2NAL `ATGAGA`, COAL `ATGTGA`.

## C++ engine (coal-cpp)

`coal-cpp` is a C++17 implementation of the same engine (standard library only, MIT). It takes the same options as
`coal` and writes the same output, warnings, report and exit codes, so it can replace `coal` in any command line.
It is a separate executable, built from source; the Python package and the Windows app use the Python engine.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
build/coal-cpp pep.aln nuc.fasta -o codon_aln.fasta       # Windows (Visual Studio): build\Release\coal-cpp.exe
```

- **Same behaviour as the Python engine:** all 285 stored PAL2NAL v14 reference checks pass with `COAL_IMPL=cpp`;
  3,000 random inputs (all formats and options, report JSON included) give identical results in both engines for
  two seeds; and the random comparison with the real PAL2NAL v14 (`tests/fuzz_vs_pal2nal.py`, 1,500 cases per seed,
  seeds 20261017, 777 and 4242) finds 0 unexplained differences (these larger runs were done locally). CI builds it
  with warnings as errors on Linux, Windows and macOS and runs the reference checks and 1,000 random C++-vs-Python
  inputs (one seed) on every push.
- **Speed:** 3-6x faster than the Python engine, with a third to a half of its peak memory (measured on a 2-core
  Linux VM, commit `4084d1f`; `benchmarks/results/cpp_vs_python_linux_4084d1f.md`):

| Sequences × codons | Python | C++ | C++ speed-up | Peak memory Python / C++ |
|---|---:|---:|---:|---:|
| 10,000 × 500 | 0.651 s | 0.109 s | 6.0× | 57 / 20 MB |
| 100 × 10,000 | 0.175 s | 0.028 s | 6.3× | 27 / 14 MB |
| 100,000 × 500 | 2.945 s | 0.964 s | 3.1× | 125 / 37 MB |

  (Default `-j` for both engines; outputs byte-identical. `--impl both` in `benchmarks/run_bench.py` reproduces it.)
- **Tests:** `COAL_IMPL=cpp COAL_CPP_BIN=build/coal-cpp python -m unittest tests.test_differential
  tests.test_cpp_vs_python tests.test_cpp_cli` (the C++ tests skip when the binary isn't built).

## Tests

```
python -m unittest discover -s tests -v
```

Compares COAL with stored PAL2NAL v14 outputs (`tests/golden/`): 39 input cases × all option combinations in all
four output formats = 285 checks (the test prints these counts). In the 6 `invalid_chars` checks COAL deliberately
fails (see above). Exit status in the other 279: PAL2NAL always exits 0, so COAL must exit 0 exactly when PAL2NAL
produced output. stdout is byte-identical in 238 checks, and the per-codon warnings are identical in 234 of them
(48 with `-nomismatch`); in the other 4 the warnings differ deliberately (NCBI tables 10, 11 and 13); in 24 both
tools fail. 23 checks are deliberate output or exit-status differences: 6 deliberate errors (`invalid_chars`, where
COAL refuses), 3 with different output and 14 `-nogap` checks with a non-standard genetic code (COAL's output must
equal PAL2NAL's output without `-nogap` with that code's stop columns removed). 238 + 24 + 23 = 285. See
`tests/README.md` for exactly what is compared and for the list of deliberate differences.

**Random differential test against the real PAL2NAL** (needs Perl; PAL2NAL is fetched, never committed):

```
python scripts/fetch_pal2nal.py
PAL2NAL=PAL2NAL/pal2nal.v14/pal2nal.pl python tests/fuzz_vs_pal2nal.py      # FUZZ_N=3000 cases by default
```

Each random case (FASTA or CLUSTAL input, sometimes with `#` marks; clean, mismatching, unrelated, masked or
wrong-code DNA; UTRs, stops, frameshift digits, RNA, lower case; same, shuffled, different or partly matching IDs;
random `-codontable`, `-nogap`, `-nomismatch`, `-blockonly`) runs through PAL2NAL and the drop-in `pal2nal` with the
same command line in all four formats, and stdout must be byte-identical, or both must fail, except for the
deliberate differences listed in the script (each is in the list above). A difference counts as deliberate only if
it is exactly the documented one: the script also runs PAL2NAL on a copy of the case with that feature neutralised
(e.g. `U` → `B`, RNA → DNA, an alternative start codon → ATG) and COAL's output must equal that, apart from the
documented change. CI runs it (3,000 cases) on every push (job `pal2nal-fuzz`). Result here, seed 20261017, 5,000
cases: 20,000 comparisons, 13,826 byte-identical, 4,681 where both fail, 1,493 deliberate differences (each
checked), 0 unexplained (per-format counts in `tests/README.md`; seeds 4242 and 777 also give 0 unexplained).

App, website and docs tests: `python -m unittest discover -s tests_app` (some need the `gui` extra).

## Benchmark results

COAL vs PAL2NAL v14 on synthetic datasets (see `benchmarks/README.md`). Each run is one
command-line process, timed from start to finish (including Perl/Python start-up).
Times are the median of 3 runs; COAL ran as a single process (`-j 1`), with the Python engine at an earlier
commit (`58a1219`, see below). **On these benchmark datasets every COAL output was byte-identical to PAL2NAL v14's.**

| Sequences × codons | Total codons | PAL2NAL | COAL | Speedup | Peak memory PAL2NAL / COAL |
|---|---:|---:|---:|---:|---:|
| 10 × 500 | 5 k | 0.051 s | 0.061 s | 0.8× | 13 / 15 MB |
| 100 × 100 | 10 k | 0.077 s | 0.062 s | 1.3× | 11 / 15 MB |
| 100 × 500 | 50 k | 0.265 s | 0.072 s | 3.7× | 15 / 15 MB |
| 100 × 2,000 | 200 k | 1.36 s | 0.093 s | 14.6× | 23 / 17 MB |
| 1,000 × 500 | 500 k | 2.39 s | 0.139 s | 17.3× | 35 / 19 MB |
| 100 × 10,000 | 1 M | 24.5 s | 0.207 s | **118.6×** | 55 / 23 MB |
| 10,000 × 500 | 5 M | 24.7 s | 0.864 s | **28.6×** | 219 / 55 MB |

What the numbers show:

- **More sequences:** COAL is 17–29× faster from 1,000 sequences upward.
- **Longer genes:** PAL2NAL's time grows much faster than the gene length: 5× more codons
  (2,000 → 10,000) took 18× longer. COAL's time grows in proportion to the input, so the gap
  widens with gene length.
- **Tiny inputs:** both tools take about 0.05 s, which is mostly start-up time.

Machine: Windows 11, Intel 64-bit CPU (Family 6 Model 183), Python 3.13.9, Strawberry Perl 5.42.3, COAL commit
`58a1219`. Raw data: `benchmarks/results/bench_windows_20261016_2356.{csv,txt}`.

**Limitations:** synthetic data only (random ancestral gene, about 10% of codons changed, about 5% deleted; no UTRs or
frameshifts in the benchmark data, although COAL handles both); one machine; 3 repeats. Speedups depend on the
machine and dataset.

### Scaling (v0.3)

Synthetic data, 500 codons per gene, on the same Windows 11 machine as above (20 logical CPUs). One run per
configuration. All outputs were byte-identical.

| Sequences | Input size | 1 process | 4 processes | All cores (20) | Peak memory, 1 process |
|---:|---:|---:|---:|---:|---:|
| 100,000 | 194 MB | 7.1 s | 2.4 s | 2.5 s | 55 MB |
| 500,000 | 974 MB | 25.2 s | 9.8 s | 10.0 s | 131 MB |

- **Memory** stays small: v0.2 needed about 2.5 GB for the 500,000-sequence input (measured on Linux).
- **Throughput** is about 100 MB of input per second with 4 processes. PAL2NAL was not run at these sizes; from its
  10,000-sequence time (24.7 s), and assuming it grows in proportion, it would need roughly 4 and 21 minutes.
- **More than 4 processes don't help yet.** The main process reads all the DNA, passes it to the workers and writes
  all the output, and that is now the bottleneck.

Raw data: `benchmarks/results/bench_windows_20261017_0059.csv` (500,000) and the 100,000 row of
`bench_windows_20261017_0055.csv`. That file's 500,000 row is not used: its 1-process column was accidentally
timed with automatic `-j`. Memory is measured for the main process only.

On the 2-core Linux development machine (v0.2 → v0.3, 1 process): 100,000 sequences 6.3 s → 3.6 s,
515 → 87 MB; 500,000 sequences 32.0 s → 17.0 s, 2,523 → 143 MB.

Reproduce with `python benchmarks/run_bench.py --no-pal2nal --sizes 100000x500 500000x500 --coal-jobs 4 0`.

### Reproduce

```
python scripts/fetch_pal2nal.py            # once; PAL2NAL is not redistributed
pip install -e . psutil                    # psutil only for the memory columns
python benchmarks/run_bench.py --repeats 3 --perl /path/to/perl
```

The script generates the same datasets every time (fixed random seed), times both tools,
checks their outputs match, and writes `benchmarks/results/<name>.csv` plus a `.txt` file
with the machine and software versions.

## Reference implementation (PAL2NAL)

COAL is an independent implementation. PAL2NAL v14 is used only as a black-box reference for differential testing
and benchmarking. PAL2NAL is distributed under the GNU General Public License v2.0 (GPL-2.0). COAL contains no
PAL2NAL code, so **PAL2NAL is not included or redistributed in this repository**, and its code and test data must
not be copied here.

To fetch it locally (into `PAL2NAL/`, which is git-ignored):

```
python scripts/fetch_pal2nal.py
```

If the official site can't be reached, the script falls back to the copy of `pal2nal.pl` in the public orthologr
R package on GitHub, accepted only if its SHA-256 equals that of the official script. If both downloads are blocked,
put `pal2nal.v14.tar.gz` (downloaded from the official site in a browser) in the repo root and run
`python scripts/fetch_pal2nal.py --archive pal2nal.v14.tar.gz`. Run it (requires Perl 5):

```
perl PAL2NAL/pal2nal.v14/pal2nal.pl pep.aln nuc.fasta -output fasta > out.fa 2> out.err
```

Reference: Suyama M, Torrents D, Bork P (2006) PAL2NAL: robust conversion of protein sequence alignments into the
corresponding codon alignments. *Nucleic Acids Res.* 34:W609–W612.

## Releasing COAL

1. Set the same version in `src/coal/__init__.py`, `pyproject.toml` and the `id="version"` span in
   `site/index.html` (a test fails if they differ), optionally add a `## <version>` section to `CHANGELOG.md` for
   the release notes, and commit.
2. Tag that commit: `git tag vX.Y.Z && git push origin vX.Y.Z`.
3. `.github/workflows/release.yml` checks the tag against the source version, builds and checks the wheel, the sdist
   and the Windows app (`desktop/build.ps1`), writes `SHA256SUMS.txt` and publishes the GitHub Release, with
   version-free copies of the downloads. **Wait until it has published, then merge to `main`**, which redeploys the
   website (https://plumat-coal.github.io/COAL/, from `site/`, via `.github/workflows/pages.yml`). The site links to
   `releases/latest/download/<version-free name>`, so it offers the new release at once.

The Windows app is not code-signed yet (the release notes say so, with the SmartScreen and Smart App Control steps).
Details and the one-time setup checklist: `site/README.md`. Building the Windows app: `desktop/BUILD.md`.

## Requirements

- Python 3.10+ (the Windows app needs nothing else)
- A C++17 compiler and CMake 3.16+, only for `coal-cpp`
- Perl 5, only for running the PAL2NAL reference

## Citation and licence

If you use COAL, please cite it (`CITATION.cff`) and PAL2NAL (above).

COAL is released under the MIT License (see `LICENSE`). This does not cover PAL2NAL, which is a separate program
licensed under GPL-2.0 and used only as an external reference for testing.
