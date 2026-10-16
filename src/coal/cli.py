"""Command-line interface.

    coal PEP_ALN NUC_FASTA [NUC_FASTA ...] [-codontable N] [-nogap] [-nomismatch] [-o OUT]

PAL2NAL-style single-dash options are accepted, so COAL can be dropped into
existing pipelines (only '-output fasta' is implemented so far).
"""
from __future__ import annotations

import argparse
import sys

from . import __version__
from .core import CoalError
from .pipeline import run
from .seqio import InputFormatError
from .tables import SUPPORTED


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
    ap.add_argument("-j", "--jobs", type=int, default=0, metavar="N",
                    help="maximum worker processes (default 0 = automatic: all CPU cores for large "
                         "inputs, one process for small ones; 1 = never use extra processes)")
    ap.add_argument("-q", "--quiet", action="store_true", help="do not print warnings")
    ap.add_argument("--version", action="version", version=f"coal {__version__}")
    args = ap.parse_args(argv)

    warn = None if args.quiet else (lambda w: print(f"coal: warning: {w}", file=sys.stderr))
    try:
        run(args.pep_aln, args.nuc_fasta, out_path=args.out, table=args.codontable,
            nogap=args.nogap, nomismatch=args.nomismatch, jobs=args.jobs, warn=warn)
    except (CoalError, InputFormatError, UnicodeDecodeError, OSError) as e:
        print(f"coal: error: {e}", file=sys.stderr)
        return 1
    return 0
