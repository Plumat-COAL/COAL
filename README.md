# COAL

A fast, codon-aware alternative to [PAL2NAL](https://www.bork.embl.de/pal2nal/):
it maps unaligned coding DNA onto an existing protein multiple sequence alignment
to produce a codon alignment.

Built at the Reinvent the Wheel AI × Science Tool-building Hackathon, London, October 2026.

> **Status:** v0.3 (pure Python, no required dependencies; numpy optional). Up to **~119× faster** than PAL2NAL
> with byte-identical output on the benchmarked datasets — see [Benchmark results](#benchmark-results).

## Install and use

```
pip install -e .
coal pep.aln nuc.fasta > codon_aln.fasta            # or: python -m coal ...
coal pep.aln nuc.fasta -codontable 2 -nogap         # PAL2NAL-style options work too
coal big_pep.fasta big_nuc.fasta -o out.fasta        # large data: uses all cores automatically
pip install -e .[fast]                              # optional numpy speed-up for large inputs
```

- **Input:** a protein alignment (FASTA or CLUSTAL) and the unaligned coding DNA (FASTA, one or more files).
  Non-letters in the DNA (gaps, digits, spaces) are ignored, as in PAL2NAL.
- **Output:** a codon alignment in FASTA format (60 columns per line), using the protein-alignment IDs.
- **Pairing:** sequences are matched by ID if every protein ID is in the DNA file(s), otherwise by order.
- **Extra sequence around the gene** (UTRs, stop codon, poly-A) is skipped: COAL finds the leftmost place in the
  DNA, in any reading frame, where the whole protein fits.
- **Genetic codes:** `-codontable N` for the 17 NCBI codes PAL2NAL supports
  (1–6, 9–16, 21–23). The first Met may use any start codon of that code.
- **Frameshifts:** a digit N in the protein alignment means N bases (PAL2NAL convention). That column is widened
  to the next multiple of 3 for every sequence.
- **Mismatches:** if the protein doesn't fit exactly, COAL uses PAL2NAL's fallback: it lines up 10-residue
  chunks and prints a warning for each codon that doesn't match its residue. `-nomismatch` removes those
  columns (and, as in PAL2NAL, then prints no mismatch warnings). If even that fails, it's an error. It is also
  an error if more than half of a gene's checkable codons (made only of A/C/G/T) don't match: that is not the
  right DNA for the protein. Codons with N or other ambiguity codes don't count, so masked genes still align,
  except that a gene with no matching 10-residue chunk and fewer than 10 checkable codons is an error when its DNA
  is longer than the gene (the start can't be told, e.g. with a 5' UTR).
  <!-- CHANGED(R3-E1): new exception. See overnight/integration/CHANGES.md -->
- **`-nogap`** removes codon columns that contain a gap or a stop codon.
- **Exit code** is non-zero on any error (PAL2NAL always exits 0). Warnings go to stderr; `-q` silences the
  per-codon warnings, but notices about probably-wrong data (no 10-residue segment found, only some IDs
  matching) are always printed.
<!-- CHANGED(R2-D6, R2-E4): was "-q silences them" and no mention of the more-than-half rule (ENG-c, as amended by R2-E3). See overnight/integration/CHANGES.md -->

**Large inputs:** COAL streams the data. It reads the protein alignment and indexes the DNA file(s) (IDs and file
positions only), then handles one sequence at a time, fetching its DNA from disk and writing its output row. So
memory doesn't grow with the size of the sequences. By default COAL uses all CPU cores for large inputs and one process
for small ones (each extra process needs about 4 million nucleotides of DNA to pay for its start-up). `-j N`
caps the number of processes and `-j 1` forces a single process. Output is always in the input order. Output is written to a temporary file and only moved into place
if the whole run succeeds, so an error never leaves a partial file. Notes:

- `-nogap`/`-nomismatch` need every row before deciding which columns to remove, so the rows are first written to
  a temporary file (disk space about the size of the output).
- CLUSTAL input is read into memory (FASTA input is streamed).
- Memory still holds the sequence IDs and their file positions, about 300 bytes per sequence (≈150 MB per
  500,000 sequences).

**Not yet supported:** output formats other than FASTA (`clustal`, `paml`, `codon`), and `-blockonly`/Gblocks input.

**Deliberate differences from PAL2NAL** (details in `tests/README.md`):

- An invalid character in the protein alignment is an error. PAL2NAL silently drops a codon.
- Genetic codes follow the current NCBI tables. In code 10, TGA = Cys (PAL2NAL treats it as a stop). In codes 11
  and 13, COAL also accepts NCBI's alternative start codons (ATT/ATC/ATA, and ATA, respectively).
- `.` is treated as a gap. In PAL2NAL it misaligns the output row.
- Ambiguous DNA bases (IUPAC codes) are accepted when every base they could stand for gives the right amino acid.
  PAL2NAL accepts only some fixed combinations.
- More than half of a gene's checkable codons not matching is an error (`protein and DNA are inconsistent`).
  PAL2NAL's fallback method (which COAL otherwise follows) accepts such input and writes unrelated DNA under
  the protein. Codons containing N/IUPAC codes are not counted.
- When some but not all protein IDs are DNA IDs, COAL pairs by order like PAL2NAL but prints a notice naming an
  example pair (PAL2NAL is silent); when no ID matches it stays silent, like PAL2NAL.
- Genetic code 3 (yeast mitochondrial) also accepts GTG as a start codon (current NCBI table; whether PAL2NAL v14
  does is unverified, so this is a possible difference).
- With `-nomismatch`, `unknown AA type. Taken as 'X'` warnings (for B/Z/J/O/U) are still printed; PAL2NAL prints
  nothing there. Mismatch warnings are not printed, as in PAL2NAL.
<!-- CHANGED(R2-D6): ENG-c, ENG-d, ENG-g, ENG-i2 added. See overnight/integration/CHANGES.md -->

## Ways to use COAL

<!-- CHANGED(R2-D6): new section. See overnight/integration/CHANGES.md -->
- **Command line:** `coal` / `python -m coal` (above).
- **Drop-in for PAL2NAL pipelines:** `pal2nal pep.aln nuc.fasta -output fasta ...` (or `python -m coal_compat`)
  accepts PAL2NAL's command line and runs COAL (see `src/coal_compat/cli.py` for what is supported).
- **Web GUI for bench scientists:** `pip install -e .[gui]`, then `python -m coal_app` (or `coal-gui`): runs on
  this computer at http://127.0.0.1:7860, with a per-gene health check, an alignment view that marks where a
  gene breaks, genetic-code suggestions and a ready-to-paste methods sentence.
- **AI agents (MCP):** `pip install -e .[mcp]`, then `python -m coal_app --mcp` (two tools: `codon_align`,
  `explain_report`); the GUI also serves MCP. Agent instructions: `.claude/skills/coal/SKILL.md`.
- **Windows app:** a desktop wrapper with an installer, no terminal needed: `desktop/README.md` (users) and
  `desktop/BUILD.md` (building it).
- **Hosted demo:** a Hugging Face Space is prepared but not deployed: `space/README.md`, `space/DEPLOY.md`.
- **Two-species preset:** `demo/diamond_pipeline.py` (DIAMOND reciprocal best hits, then COAL on every pair).

## Tests

```
python -m unittest discover -s tests -v
```

<!-- CHANGED(ENG-i4): was "165 checks. COAL's output must be byte-identical, and its mismatch warnings must equal PAL2NAL's" (rc files were not read, warnings not compared with -nomismatch, 4 checks compared only the exit code). See overnight/engine-fixes/CHANGES.md -->
Compares COAL with stored PAL2NAL v14 outputs (`tests/golden/`): 38 input cases × all FASTA option
combinations = 165 checks (the test prints these counts). In the 3 `invalid_chars` checks COAL deliberately fails
(see below). Exit status in the other 162: PAL2NAL always exits 0, so COAL must exit 0 exactly when PAL2NAL produced
output.<!-- CHANGED(R3-D5): was "in all 165" checks; the harness checks exit status in 162. See overnight/integration/CHANGES.md --> stdout is byte-identical in 147 checks, and the per-codon warnings
are identical in 141 of them (47 with `-nomismatch`); in 12 both tools fail. 12 checks are deliberate differences:
the other 6 of the 147 (same output, different warnings), 3 with different output and 3 where COAL deliberately
fails. See `tests/README.md` for exactly what is compared and for the list of deliberate differences.

## Benchmark results

COAL vs PAL2NAL v14 on synthetic datasets (see `benchmarks/README.md`). Each run is one
command-line process, timed from start to finish (including Perl/Python start-up).
Times are the median of 3 runs. **Every COAL output was byte-identical to PAL2NAL's.**

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

**Limitations:** synthetic data only (random ancestral gene, about 10% of codons changed,
about 5% deleted; no UTRs or frameshifts, which COAL v0.1 doesn't support); one machine;
3 repeats. Speedups depend on the machine and dataset.

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
  <!-- CHANGED(R3-D4): was 27 s, 4.5 and 22 minutes; the table and bench_windows_20261016_2356.csv say 24.7 s.
  See overnight/integration/CHANGES.md -->
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

COAL is an independent implementation. PAL2NAL v14 is used only as a
black-box reference for differential testing and benchmarking. PAL2NAL is
distributed under the GNU General Public License v2.0 (GPL-2.0). COAL contains
no PAL2NAL code, so **PAL2NAL is not included or redistributed in this
repository**, and its code and test data must not be copied here.

To fetch it locally (into `PAL2NAL/`, which is git-ignored):

```
python scripts/fetch_pal2nal.py
```

If the download is blocked, put `pal2nal.v14.tar.gz` (downloaded from the
official site in a browser) in the repo root and run:

```
python scripts/fetch_pal2nal.py --archive pal2nal.v14.tar.gz
```

Run it (requires Perl 5):

```
perl PAL2NAL/pal2nal.v14/pal2nal.pl pep.aln nuc.fasta -output fasta > out.fa 2> out.err
```

Reference: Suyama M, Torrents D, Bork P (2006) PAL2NAL: robust conversion of
protein sequence alignments into the corresponding codon alignments.
*Nucleic Acids Res.* 34:W609–W612.

## Releasing COAL

<!-- WEB-5: new section. See overnight/integration/CHANGES.md -->
Tag, and everything else is automatic:

1. Set the same version in `src/coal/__init__.py`, `pyproject.toml` and the `id="version"` span in
   `site/index.html` (a test fails if they differ) (optionally add a `## <version>` section to
   `CHANGELOG.md` for the release notes; without one, GitHub generates notes from merged pull requests), commit,
   merge to `main`.
2. `git tag v0.3.0 && git push origin v0.3.0`
3. `.github/workflows/release.yml` checks the tag against the source version, builds and checks the wheel, the sdist
   and the Windows app (`desktop/build.ps1`: tests, PyInstaller, frozen-engine check, Inno Setup installer), writes
   `SHA256SUMS.txt` and publishes the GitHub Release, with version-free copies of the downloads.
4. The website (https://plumat-coal.github.io/COAL/, from `site/`) links to
   `releases/latest/download/<version-free name>`, so it offers the new release at once; only the version shown
   on the page is edited by hand (step 1).<!-- CHANGED(WEB-8): was "nothing to edit"; the page now names the
   current version. See overnight/integration/CHANGES.md -->
   The site itself is redeployed by `.github/workflows/pages.yml` when `site/` changes on `main`.

The Windows app is not code-signed yet (the release notes say so, with the SmartScreen and Smart App Control
steps). Details, file names and the one-time setup checklist (enable Pages, first tag, custom domain, Zenodo DOI):
`site/README.md`.

## Requirements

- Python 3.10+
- Perl 5 (only for running the PAL2NAL reference)

## Licence

COAL is released under the MIT License (see `LICENSE`). This does not cover PAL2NAL, which is a separate program licensed under GPL-2.0 and used only as an external reference for testing.
