# COAL

A fast, codon-aware alternative to [PAL2NAL](https://www.bork.embl.de/pal2nal/):
it maps unaligned coding DNA onto an existing protein multiple sequence alignment
to produce a codon alignment.

Built at the Reinvent the Wheel AI × Science Tool-building Hackathon, London, October 2026.

> **Status:** v0.1 MVP (pure Python, no dependencies).

## Install and use

```
pip install -e .
coal pep.aln nuc.fasta > codon_aln.fasta      # or: python -m coal ...
```

Input: protein alignment (FASTA or CLUSTAL) and unaligned coding DNA (FASTA).
Output: codon alignment in FASTA (60 columns per line), using protein-alignment IDs.
Sequences are paired by ID when all protein IDs exist in the DNA file, otherwise by order.
Exit code is non-zero on any error (unlike PAL2NAL, which exits 0).

**v0.1 limitations:** genetic code 1 only; FASTA output only; DNA must be exactly the
coding sequence (an optional trailing stop codon is trimmed); no UTRs, frameshift markers,
or mismatch tolerance — these give a clear error.

## Tests

```
python -m unittest discover -s tests -v
```

Compares COAL with stored PAL2NAL v14 outputs (`tests/golden/`); see `tests/README.md`.

## Reference implementation (PAL2NAL)

COAL is an independent implementation. PAL2NAL v14 is used only as a
black-box reference for differential testing and benchmarking. Its distribution
has no licence file, so **it is not included in this repository** and its code
and test data must not be copied here.

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

## Requirements

- Python 3.10+
- Perl 5 (only for running the PAL2NAL reference)

## Licence

COAL is released under the MIT License (see `LICENSE`). This does not cover PAL2NAL.
