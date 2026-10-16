"""Command-line interface: coal PEP_ALN NUC_FASTA [-o OUT]"""
from __future__ import annotations

import argparse
import sys

from . import __version__
from .core import CoalError, codon_align
from .seqio import InputFormatError, read_alignment, read_fasta

WRAP = 60


def write_fasta(records, fh, wrap: int = WRAP) -> None:
    for rid, seq in records:
        fh.write(f">{rid}\n")
        for i in range(0, len(seq), wrap):
            fh.write(seq[i:i + wrap] + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="coal",
        description="Convert a protein alignment + unaligned coding DNA into a codon alignment.",
    )
    ap.add_argument("pep_aln", help="protein alignment (FASTA or CLUSTAL)")
    ap.add_argument("nuc_fasta", help="coding DNA sequences (FASTA)")
    ap.add_argument("-o", "--output", help="output file (default: stdout)")
    ap.add_argument("--version", action="version", version=f"coal {__version__}")
    args = ap.parse_args(argv)

    try:
        result = codon_align(read_alignment(args.pep_aln), read_fasta(args.nuc_fasta))
    except (CoalError, InputFormatError, UnicodeDecodeError, OSError) as e:
        print(f"coal: error: {e}", file=sys.stderr)
        return 1

    if args.output:
        with open(args.output, "w", newline="\n") as fh:
            write_fasta(result, fh)
    else:
        sys.stdout.reconfigure(newline="\n")
        write_fasta(result, sys.stdout)
    return 0
