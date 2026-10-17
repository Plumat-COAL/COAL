"""Command-line interface.

    coal PEP_ALN NUC_FASTA [NUC_FASTA ...] [-codontable N] [-nogap] [-nomismatch] [-o OUT]

PAL2NAL-style single-dash options are accepted, so COAL can be dropped into
existing pipelines; -output clustal|paml|fasta|codon and -blockonly as in PAL2NAL v14.
"""
from __future__ import annotations

import argparse
import json
import os
import signal
import sys

from . import __version__
from .core import CoalError
from .pipeline import FORMATS, run
from .report import RunReport
from .seqio import InputFormatError
from .tables import SUPPORTED


def _on_sigterm(signum, frame):
    raise SystemExit(128 + signum)


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
                    help="remove codon columns containing gaps or stop codons (stop codons of the -codontable "
                         "code, any case; PAL2NAL uses the standard code's TAA/TAG/TGA, upper case only)")
    # CHANGED(FZ-5): help said "stop codons" (the standard code's, upper case, as PAL2NAL); now the selected code's
    ap.add_argument("-nomismatch", "--nomismatch", action="store_true",
                    help="remove alignment columns where a codon does not match its residue")
    # CHANGED(DEC-4): was choices=["fasta"]; now PAL2NAL v14's four formats, byte-identical to PAL2NAL. The default
    # stays fasta (deliberate difference; the drop-in `pal2nal` command defaults to clustal like PAL2NAL)
    ap.add_argument("-output", "--output-format", default="fasta", choices=list(FORMATS),
                    help="output format, as PAL2NAL v14 writes it (default fasta; PAL2NAL's default is clustal)")
    ap.add_argument("-blockonly", "--blockonly", action="store_true",
                    help="keep only the columns marked '#' in the CLUSTAL input's conservation lines (PAL2NAL's "
                         "-blockonly; no effect if the input has no '#' marks)")
    ap.add_argument("-o", "--out", help="output file (default: stdout)")
    ap.add_argument("-j", "--jobs", type=int, default=0, metavar="N",
                    help="maximum worker processes (default 0 = automatic: all CPU cores for large "
                         "inputs, one process for small ones; 1 = never use extra processes)")
    # CHANGED(R2-E4): help was "do not print warnings" (notices are printed anyway)
    # CHANGED(DEC-2): was "notices about probably-wrong data are still printed"; now -q silences them too
    ap.add_argument("-q", "--quiet", action="store_true",
                    help="do not print warnings or notices (errors are always printed)")
    # CHANGED(DEC-1): new. Was: these three rules were always errors (ENG-c, R3-E1) or a warning (ENG-d). Now
    # PAL2NAL's behaviour by default (output + notice) and errors only with --strict. See overnight/integration/CHANGES.md
    # CHANGED(REV-1): help said "(PAL2NAL never refuses)"; PAL2NAL does refuse some inputs (FIN-C17), it writes these
    ap.add_argument("-strict", "--strict", action="store_true",
                    help="refuse probably-wrong data instead of writing it as PAL2NAL does (PAL2NAL writes these): "
                         "exit 1 with no output if (1) more than half of a sequence's checkable codons (A/C/G/T "
                         "only) do not match its protein; (2) no 10-residue segment of a protein is found in its "
                         "DNA, the DNA is longer than the gene and fewer than 10 codons are checkable (where the gene "
                         "starts can't be told); (3) only some protein IDs are DNA IDs (pairing by order would "
                         "follow). Without --strict these are notices on stderr")
    # CHANGED(MRG-3, MRG-4): the teammate's options (new). See overnight/integration/CHANGES.md
    ap.add_argument("--keep-going", "-keep-going", action="store_true",
                    help="skip sequences that cannot be converted (or, with --strict, are refused) instead of "
                         "stopping; the others are written in input order, skipped ones are named on stderr and in "
                         "the report; exit code 2 if any was skipped")
    ap.add_argument("--max-warnings", type=int, default=0, metavar="N",
                    help="print at most N per-codon warning lines, then a count of the rest (default 0 = no limit, "
                         "as PAL2NAL); notices and errors are always printed")
    ap.add_argument("--progress", action="store_true",
                    help="show percent done and time left on stderr (inputs over 10,000 sequences)")
    ap.add_argument("--resume", action="store_true",
                    help="continue an interrupted run (needs -o, same inputs and options) from its checkpoint "
                         "<out>.coal-checkpoint.json; without a checkpoint it runs normally. Output is identical to "
                         "an uninterrupted run. With -nogap, -nomismatch, -blockonly, clustal or codon output only "
                         "the first pass resumes; the final pass is redone. Not for piped inputs or an output that "
                         "is not a regular file")
    ap.add_argument("--report", metavar="FILE.json",
                    help="write a JSON run report (schema coal.report/1), also when the run fails")
    ap.add_argument("--json", action="store_true",
                    help="print the run report to stdout (requires -o, so stdout is free)")
    ap.add_argument("--version", action="version", version=f"coal {__version__}")
    argv = sys.argv[1:] if argv is None else list(argv)
    args = ap.parse_args(argv)
    if args.json and not args.out:
        ap.error("--json requires -o (the alignment and the report cannot both go to stdout)")
    if args.resume and not args.out:
        ap.error("--resume requires -o (only runs that write to a file can be resumed)")
    if args.max_warnings < 0:
        ap.error("--max-warnings must be 0 (no limit) or more")

    # CHANGED(MRG-2): always collected (it also drives the end-of-run "coal: summary:" lines); per-sequence entries
    # of ok sequences only when a report is asked for
    rec = RunReport(argv, args.pep_aln, args.nuc_fasta, args.out, args.output_format, args.codontable,
                    args.nogap, args.nomismatch, keep_going=args.keep_going,
                    keep_entries=bool(args.report or args.json), strict=args.strict, blockonly=args.blockonly)
    hidden = [0]
    shown = [0]

    def _warn(w: str):
        # CHANGED(MRG-3): --max-warnings caps the per-codon lines only (teammate's option)
        if args.max_warnings and "pepAlnPos" in w:
            if shown[0] >= args.max_warnings:
                hidden[0] += 1
                return
            shown[0] += 1
        print(f"coal: warning: {w}", file=sys.stderr)
    warn = None if args.quiet else _warn
    # CHANGED(ENG-c): new; notices are printed even with -q. See overnight/engine-fixes/CHANGES.md
    # CHANGED(DEC-2): was printed even with -q; now -q silences notices too (GNU / clig.dev "quiet", and PAL2NAL's
    # -nostderr). Errors are still printed. See overnight/integration/CHANGES.md
    notice = None if args.quiet else (lambda w: print(f"coal: warning: {w}", file=sys.stderr))
    # CHANGED(ENG-UX-10): was "coal: note:"; now "coal: notice:" (the three prefixes: error, warning, notice)
    note = None if args.quiet else (lambda t: print(f"coal: notice: {t}", file=sys.stderr))
    # CHANGED(ENG-f): was default signal handling (SIGTERM left .coal_*.tmp files; Ctrl-C printed a
    # traceback; a closed stdout pipe gave "[Errno 32] Broken pipe", exit 1); now SIGTERM -> SystemExit
    # (cleanup runs, exit 143), Ctrl-C -> "coal: interrupted", exit 130, broken pipe -> silent exit 141.
    # See overnight/engine-fixes/CHANGES.md
    signal.signal(signal.SIGTERM, _on_sigterm)
    if hasattr(signal, "SIGHUP"):  # CHANGED(R2-E8): terminal closed -> same cleanup, exit 129
        signal.signal(signal.SIGHUP, _on_sigterm)
    rc = 0
    try:
        skipped = run(args.pep_aln, args.nuc_fasta, out_path=args.out, table=args.codontable,
                      nogap=args.nogap, nomismatch=args.nomismatch, jobs=args.jobs, warn=warn, notice=notice,
                      strict=args.strict, out_format=args.output_format, blockonly=args.blockonly,
                      rec=rec, keep_going=args.keep_going, progress=args.progress, resume=args.resume, note=note)
        if skipped:
            rc = 2  # CHANGED(MRG-3): --keep-going skipped some sequences (partial success)
    except KeyboardInterrupt:
        print("coal: interrupted", file=sys.stderr)
        return 130
    except BrokenPipeError:
        try:  # stop Python from reporting the pipe again when it flushes stdout at exit
            os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        except (OSError, ValueError):
            pass
        return 141
    except (CoalError, InputFormatError, UnicodeDecodeError, OSError) as e:
        print(f"coal: error: {e}", file=sys.stderr)
        rc = 1
        code = (e.code if isinstance(e, CoalError) else
                "input_format" if isinstance(e, (InputFormatError, UnicodeDecodeError)) else "file_error")
        rec.fail(code, str(e), getattr(e, "seq_id", None))
    if hidden[0]:
        where = "see the report" if (args.report or args.json) else "see them all with --report FILE.json"
        print(f"coal: warning: ... and {hidden[0]:,} more warnings ({where})", file=sys.stderr)
    rec.exit_code = rc
    if not args.quiet:
        for line in rec.summary_lines():
            print(f"coal: summary: {line}", file=sys.stderr)
    if args.report:
        try:
            rec.write(args.report)
        except OSError as e:
            # CHANGED(ENG-UX-12): a folder gave "Is a directory"; now said in plain words (C++ said only the path)
            why = "it is a folder, not a file" if os.path.isdir(args.report) else (e.strerror or e)
            print(f"coal: error: could not write the report {args.report}: {why}", file=sys.stderr)
            rc = rc or 1
    if args.json:
        print(json.dumps(rec.to_dict(), indent=1))
    return rc
