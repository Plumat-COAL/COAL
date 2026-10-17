"""Command-line interface.

    coal PEP_ALN NUC_FASTA [NUC_FASTA ...] [-codontable N] [-nogap] [-nomismatch] [-o OUT]

PAL2NAL-style single-dash options are accepted, so COAL can be dropped into
existing pipelines. Output formats: fasta (default), clustal, paml, codon.
"""
from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .core import CoalError
from .pipeline import FORMATS, run
from .report import RunReport
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
    ap.add_argument("-output", "--output-format", default="fasta", choices=list(FORMATS),
                    help="output format: fasta (default), clustal, paml or codon")
    ap.add_argument("-o", "--out", help="output file (default: stdout)")
    ap.add_argument("-j", "--jobs", type=int, default=0, metavar="N",
                    help="maximum worker processes (default 0 = automatic: all CPU cores for large "
                         "inputs, one process for small ones; 1 = never use extra processes)")
    ap.add_argument("-q", "--quiet", action="store_true", help="do not print warnings")
    ap.add_argument("--report", metavar="FILE.json",
                    help="write a JSON run report (schema coal.report/1), also when the run fails")
    ap.add_argument("--json", action="store_true",
                    help="print the run report to stdout (requires -o, so stdout is free)")
    ap.add_argument("--version", action="version", version=f"coal {__version__}")
    argv = sys.argv[1:] if argv is None else list(argv)
    args = ap.parse_args(argv)
    if args.json and not args.out:
        ap.error("--json requires -o (the alignment and the report cannot both go to stdout)")

    # Always collected: it also drives the end-of-run "coal: summary:" lines.
    rec = RunReport(argv, args.pep_aln, args.nuc_fasta, args.out, args.output_format, args.codontable,
                    args.nogap, args.nomismatch, keep_entries=bool(args.report or args.json))
    warn = None if args.quiet else (lambda w: print(f"coal: warning: {w}", file=sys.stderr))
    rc = 0
    try:
        run(args.pep_aln, args.nuc_fasta, out_path=args.out, table=args.codontable,
            nogap=args.nogap, nomismatch=args.nomismatch, jobs=args.jobs, warn=warn,
            fmt=args.output_format, rec=rec)
    except (CoalError, InputFormatError, UnicodeDecodeError, OSError) as e:
        print(f"coal: error: {e}", file=sys.stderr)
        rc = 1
        if rec is not None:
            code = (e.code if isinstance(e, CoalError) else
                    "input_format" if isinstance(e, (InputFormatError, UnicodeDecodeError)) else "file_error")
            rec.fail(code, str(e))
    if rec is not None:
        rec.exit_code = rc
        try:
            rec.compute_hints()
        except (CoalError, InputFormatError, UnicodeDecodeError, OSError):
            pass  # a hint is optional; never let it hide the real result
        if not args.quiet:
            for line in rec.summary_lines():
                print(f"coal: summary: {line}", file=sys.stderr)
        if args.report:
            try:
                rec.write(args.report)
            except OSError as e:
                print(f"coal: error: could not write the report: {e}", file=sys.stderr)
                rc = rc or 1
        if args.json:
            print(json.dumps(rec.to_dict(), indent=1))
    return rc
