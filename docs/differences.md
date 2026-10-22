# Deliberate differences from PAL2NAL

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
  accepts only some fixed combinations, and otherwise falls back to another, frame-shifted placement, or fails.
  Measured with PAL2NAL v14 on a 40-residue gene with one such codon: `CTN`, `CGN`, `AGY`, `TCN`, `GTN`, `GTK`,
  `ATY`, `AAR`, `GAY`, `CCN` are accepted by both tools; `YTA` and `YTG` (both readings Leu), `MGA` and `MGG` (both
  Arg) and `ATH` (all three Ile) are accepted by COAL but make PAL2NAL v14 fail the gene ("inconsistency", exit 0,
  no output). Real coding sequences with IUPAC codes (e.g. from Ensembl) can hit this.
- Probably-wrong data (see [Using COAL](usage.md)) gives the same output as PAL2NAL plus a notice on stderr (PAL2NAL prints only its
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
