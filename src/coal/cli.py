"""Command-line interface.

    coal PEP_ALN NUC_FASTA [NUC_FASTA ...] [-codontable N] [-nogap] [-nomismatch] [-o OUT]

PAL2NAL-style single-dash options are accepted, so COAL can be dropped into
existing pipelines (only '-output fasta' is implemented so far).
"""
from __future__ import annotations

import argparse
import sys

from . import __version__
from .core import CoalError, codon_align
from .seqio import InputFormatError, read_alignment, read_dna
from .tables import SUPPORTED

WRAP = 60


def write_fasta(records, fh, wrap: int = WRAP) -> None:
    for rid, seq in records:
        fh.write(f">{rid}\n")
        if not seq:
            fh.write("\n")  # PAL2NAL prints an empty line for an empty sequence
        for i in range(0, len(seq), wrap):
            fh.write(seq[i:i + wrap] + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="coal",
        description="Convert a protein alignment + unaligned coding DNA into a codon alignment.",
    )
    ap.add_argument("pep_aln", help="protein alignment (FASTA or CLUSTAL)")
    ap.add_argument("nuc_fasta", nargs="+", help="coding DNA sequences (FASTA; one or more files)")
    ap.add_argument("-codontable", "--codontable", type=int, default=1, choices=SUPPORTED,
                    metavar="N", help=f"NCBI genetic code (default 1); one of {', '.join(map(str, SUPPORTED))}")
    ap.add_argument("-nogap", "--nogap", action="store_true",
                    help="remove codon columns containing gaps or stop codons")
    ap.add_argument("-nomismatch", "--nomismatch", action="store_true",
                    help="remove alignment columns where a codon does not match its residue")
    ap.add_argument("-output", "--output-format", default="fasta", choices=["fasta"],
                    help="output format (only fasta so far)")
    ap.add_argument("-o", "--out", help="output file (default: stdout)")
    ap.add_argument("-q", "--quiet", action="store_true", help="do not print warnings")
    ap.add_argument("--version", action="version", version=f"coal {__version__}")
    args = ap.parse_args(argv)

    try:
        result, warnings = codon_align(
            read_alignment(args.pep_aln), read_dna(args.nuc_fasta),
            table=args.codontable, nogap=args.nogap, nomismatch=args.nomismatch,
        )
    except (CoalError, InputFormatError, UnicodeDecodeError, OSError) as e:
        print(f"coal: error: {e}", file=sys.stderr)
        return 1

    if not args.quiet:
        for w in warnings:
            print(f"coal: warning: {w}", file=sys.stderr)
    if args.out:
        with open(args.out, "w", newline="\n") as fh:
            write_fasta(result, fh)
    else:
        if hasattr(sys.stdout, "reconfigure"):
            sys.stdout.reconfigure(newline="\n")
        write_fasta(result, sys.stdout)
    return 0
