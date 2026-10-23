"""Command-line interface.

    coal PEP_ALN NUC_FASTA [NUC_FASTA ...] [-codontable N] [-nogap] [-nomismatch] [-o OUT]
    coal --batch JOBS [defaults for every job]          (coal/batch.py: many small runs in one process)

`-` as one input path (the protein alignment or one DNA file) reads standard input.

PAL2NAL-style single-dash options are accepted, so COAL can be dropped into
existing pipelines; -output clustal|paml|fasta|codon and -blockonly as in PAL2NAL v14.
"""
from __future__ import annotations

import argparse
import os
import signal
import sys

from . import __version__
from .core import CoalError, limit_blas_threads
from .pipeline import FORMATS, STDIN_TWICE, run, stdin_count
from .report import IdShortener, RunReport, report_folder_missing, report_path_conflict
from .tables import SUPPORTED


MAX_WARNINGS = 20  # default --max-warnings
# The same words in coal-cpp (main.cpp, std::bad_alloc)
OUT_OF_MEMORY = ("out of memory: this run needs more memory than is free on this machine; run with -j 1 (fewer worker processes), split the input, or use a machine with more memory")
BATCH_DESTS = ("batch", "batch_out_dir")  # the Batch group's options (coal only; coal-cpp rejects them)


def _on_sigterm(signum, frame):
    raise SystemExit(128 + signum)


def closed_stdout_exit() -> int:
    """After a BrokenPipeError on stdout (`coal ... | head`): the exit status 141 (128 + SIGPIPE, docs/usage.md
    "Closed output pipe"), with stdout pointed at the null device first so Python does not report the pipe again when
    it flushes stdout at exit. Also used by the `pal2nal` drop-in."""
    try:
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
    except (OSError, ValueError):
        pass
    return 141


def _ranges(nums) -> str:
    """1, 2, 3, 5 -> '1-3, 5'."""
    out, run = [], [nums[0]]
    for n in nums[1:]:
        if n == run[-1] + 1:
            run.append(n)
        else:
            out.append(run)
            run = [n]
    out.append(run)
    return ", ".join(f"{r[0]}-{r[-1]}" if len(r) > 2 else ", ".join(map(str, r)) for r in out)


def build_parser(batch: bool = False) -> argparse.ArgumentParser:
    """The `coal` command line.
    Five groups (Core / PAL2NAL options / Large runs / Reporting / Batch) with one-line help; the details are in the
    README. coal-cpp's cpp/src/help_text.hpp is generated from this parser (tests/test_cpp_cli.py --write-help).
    batch=True: the parser of a `coal --batch` command line, whose inputs come from the JOBS file (the positionals
    are optional and refused in cli.main); everything else is the same."""
    ap = argparse.ArgumentParser(
        prog="coal",
        description="Convert a protein alignment + unaligned coding DNA into a codon alignment.",
        epilog="Every option in detail: README.md (https://github.com/plumat-coal/COAL).",
    )
    core = ap.add_argument_group("Core")
    pal = ap.add_argument_group("PAL2NAL options")
    large = ap.add_argument_group("Large runs")
    rep = ap.add_argument_group("Reporting")
    # registration order = the usage line's order (unchanged)
    core.add_argument("pep_aln", nargs="?" if batch else None,
                      help="protein alignment (FASTA or CLUSTAL; - reads standard input)")
    core.add_argument("nuc_fasta", nargs="*" if batch else "+",
                      help="coding DNA (FASTA; one or more files; - reads standard input)")
    core.add_argument("-codontable", "--codontable", type=int, default=1, choices=SUPPORTED, metavar="N",
                      help=f"NCBI genetic code (default 1): {_ranges(list(SUPPORTED))}")
    pal.add_argument("-nogap", "--nogap", action="store_true",
                     help="remove codon columns with a gap or a stop codon (of the -codontable code)")
    pal.add_argument("-nomismatch", "--nomismatch", action="store_true",
                     help="remove columns where a codon does not match its amino acid")
    # the {fasta,...} choices stay in the help: the app reads the formats from it (coal_app.core.detect_features)
    pal.add_argument("-output", "--output-format", default="fasta", choices=list(FORMATS),
                     help="output format (default fasta; PAL2NAL's default is clustal)")
    pal.add_argument("-blockonly", "--blockonly", action="store_true",
                     help="keep only the columns marked '#' in a CLUSTAL input")
    core.add_argument("-o", "--out", metavar="FILE", help="output file (default: standard output)")
    large.add_argument("-j", "--jobs", type=int, default=0, metavar="N",
                       help="worker processes (default 0 = automatic; 1 = one process)")
    core.add_argument("-q", "--quiet", action="store_true",
                      help="do not print warnings or notices (errors are always printed)")
    core.add_argument("-strict", "--strict", action="store_true",
                      help="refuse probably-wrong genes (exit 1) instead of writing them as PAL2NAL does: more than "
                           "half of the checkable codons do not match, no 10-residue segment places the gene, or "
                           "only some protein IDs are DNA IDs")
    large.add_argument("--keep-going", "-keep-going", action="store_true",
                       help="leave out sequences that cannot be converted, write the rest; exit 2 if any was left out")
    # Capped by default: a wrong genetic code would print hundreds of lines (834 on coxI). The `pal2nal` drop-in
    # passes 0 (PAL2NAL prints them all).
    rep.add_argument("--max-warnings", type=int, default=MAX_WARNINGS, metavar="N",
                     help=f"print at most N per-codon warnings (default {MAX_WARNINGS}; 0 = all)")
    large.add_argument("--progress", action="store_true",
                       help="percent done and time left on stderr (inputs over 10,000 sequences)")
    large.add_argument("--resume", action="store_true",
                       help="continue an interrupted run from its checkpoint (needs -o and the same inputs and "
                            "options)")
    rep.add_argument("--report", metavar="FILE.json", help="write a JSON run report, also when the run fails")
    rep.add_argument("--json", action="store_true", help="print the run report to stdout (needs -o)")
    # The Batch group comes after every other group, so the help of a single run is unchanged but for these lines.
    # The help strings say "not in coal-cpp" because coal-cpp prints this same help and refuses the two options.
    bat = ap.add_argument_group("Batch")
    bat.add_argument("--batch", metavar="JOBS",
                     help="run many small conversions in one process: one coal command line per line of JOBS "
                          "(- = stdin), each with -o; the options given here are the defaults of every job "
                          "(not in coal-cpp)")
    bat.add_argument("--batch-out-dir", metavar="DIR",
                     help="with --batch: a job without -o writes DIR/<its protein file's name>.<fasta|aln|paml|codon>")
    ap.add_argument("--version", action="version", version=f"coal {__version__}")
    return ap


def wants_batch(argv) -> bool:
    """True if the command line has `--batch` (exactly, or `--batch=JOBS`) before any `--`. `--bat` and shorter are
    ambiguous with --batch-out-dir, so argparse refuses them either way."""
    for a in argv:
        if a == "--":
            return False
        if a == "--batch" or a.startswith("--batch="):
            return True
    return False


def _maybe_cpp(ap, args, argv) -> int | None:
    """Run a large input with coal-cpp (coal/dispatch.py); its exit code, or None to go on in Python."""
    from . import dispatch
    note = None if args.quiet else (lambda t: print(f"coal: notice: {t}", file=sys.stderr))
    try:
        path = dispatch.choose(args, argv, {a.dest: a.default for a in ap._actions}, note)
    except dispatch.EngineError as e:
        print(f"coal: error: {e}", file=sys.stderr)
        return 1
    if path is None:
        return None
    rc = dispatch.run_cpp(path, argv)
    if isinstance(rc, int):
        return rc
    if (os.environ.get("COAL_ENGINE") or "").strip().lower() == "cpp":
        print(f"coal: error: COAL_ENGINE=cpp: {rc}", file=sys.stderr)
        return 1
    if note:
        note(f"{rc}; using the Python engine")
    return None


def install_signal_handlers() -> None:
    """SIGTERM -> SystemExit (cleanup runs, so no .coal_*.tmp files are left; exit 143), Ctrl-C -> "coal: interrupted",
    exit 130, broken pipe -> silent exit 141. SIGHUP (terminal closed) -> same cleanup, exit 129."""
    signal.signal(signal.SIGTERM, _on_sigterm)
    # An ignored SIGHUP (nohup, setsid) stays ignored: installing the handler regardless made `nohup coal ...` die
    # with exit 129 when the terminal closed
    if hasattr(signal, "SIGHUP") and signal.getsignal(signal.SIGHUP) != signal.SIG_IGN:
        signal.signal(signal.SIGHUP, _on_sigterm)


def run_job(args, argv, prefix: str = "", want_report: bool = False) -> tuple[int, RunReport]:
    """One conversion from parsed arguments: the body of a single `coal` run, also used for every job of a batch.
    Returns (exit code, the run report); 130 and 141 mean interrupted / closed pipe (the caller stops).
    Every stderr line is `coal: <prefix><kind>: ...` (prefix "" for a single run, "[k/n PEP] " in a batch).
    want_report: keep the per-sequence entries although args has no --report/--json (the batch report wants them)."""
    reporting = bool(args.report or args.json or want_report)
    # Always collected (it also drives the end-of-run "coal: summary:" lines); per-sequence entries
    # of ok sequences only when a report is asked for
    rec = RunReport(argv, args.pep_aln, args.nuc_fasta, args.out, args.output_format, args.codontable,
                    args.nogap, args.nomismatch, keep_going=args.keep_going, keep_entries=reporting,
                    strict=args.strict, blockonly=args.blockonly,
                    hint=not args.quiet or reporting)  # -q never shows it
    shorten = IdShortener(rec)  # IDs over 200 characters are shortened on stderr (the report keeps them whole)

    def say(text: str):
        print(f"coal: {prefix}{shorten(text)}", file=sys.stderr)
    hidden = [0]
    shown = [0]

    def _warn(w: str):
        # --max-warnings caps the per-codon lines only
        if args.max_warnings and "pepAlnPos" in w:
            if shown[0] >= args.max_warnings:
                hidden[0] += 1
                return
            shown[0] += 1
        say(f"warning: {w}")
    warn = None if args.quiet else _warn
    # -q silences notices too (GNU / clig.dev "quiet", and PAL2NAL's -nostderr). Errors are still printed.
    notice = None if args.quiet else (lambda w: say(f"warning: {w}"))
    # the three prefixes: error, warning, notice
    note = None if args.quiet else (lambda t: say(f"notice: {t}"))
    rc = 0
    try:
        skipped = run(args.pep_aln, args.nuc_fasta, out_path=args.out, table=args.codontable,
                      nogap=args.nogap, nomismatch=args.nomismatch, jobs=args.jobs, warn=warn, notice=notice,
                      strict=args.strict, out_format=args.output_format, blockonly=args.blockonly,
                      rec=rec, keep_going=args.keep_going, progress=args.progress, resume=args.resume, note=note)
        if skipped:
            rc = 2  # --keep-going skipped some sequences (partial success)
    except KeyboardInterrupt:
        say("interrupted")
        return 130, rec
    except BrokenPipeError:
        return closed_stdout_exit(), rec
    except (CoalError, UnicodeDecodeError, OSError) as e:  # InputFormatError is a CoalError (code input_format)
        say(f"error: {e}")
        rc = 1
        code = (e.code if isinstance(e, CoalError) else
                "input_format" if isinstance(e, UnicodeDecodeError) else "file_error")
        rec.fail(code, str(e), getattr(e, "seq_id", None))
    except MemoryError:
        # The main process (a worker that dies of it is "a worker process died"): one line, exit 1, no traceback
        say(f"error: {OUT_OF_MEMORY}")
        rc = 1
        rec.fail("out_of_memory", OUT_OF_MEMORY, None)
    if hidden[0]:
        where = "see the report" if reporting else "see them all with --max-warnings 0 or --report FILE.json"
        say(f"warning: ... and {hidden[0]:,} more warnings ({where})")
    rec.exit_code = rc
    if not args.quiet:
        for line in rec.summary_lines():
            say(f"summary: {line}")
    if args.report:
        try:
            rec.write(args.report)
        except OSError as e:
            why = "it is a folder, not a file" if os.path.isdir(args.report) else (e.strerror or e)
            say(f"error: could not write the report {args.report}: {why}")
            rc = rc or 1
    if args.json:
        try:
            rec.write_json(sys.stdout)
            print(flush=True)
        except BrokenPipeError:  # `--json | head`: silent exit 141, as for the alignment
            return closed_stdout_exit(), rec
    return rc, rec


def cpp_path_main() -> int:
    """`coal --cpp-path`: print the path of the coal-cpp that `coal` would hand large inputs to (the lookup order of
    coal/dispatch.py: COAL_CPP, the binary bundled in a platform wheel, next to Python, PATH) and exit 0; exit 1 with
    the reason on stderr when there is none or it is not this version. For pipelines and install checks:
    `cpp=$(coal --cpp-path) && "$cpp" pep.aln nuc.fasta`. Handled before argparse so coal-cpp's own help and parser,
    which mirror `coal --help`, stay unchanged."""
    from . import dispatch
    path = dispatch.find_cpp()
    if path is None:
        print(f"coal: coal-cpp not found (looked at {dispatch.LOOKUP_ORDER})", file=sys.stderr)
        return 1
    ok, why = dispatch._version_ok(path)
    if not ok:
        print(f"coal: {why}", file=sys.stderr)
        return 1
    print(path)
    return 0


def main(argv=None) -> int:
    command_line = argv is None  # only the command line may hand the run to coal-cpp
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv == ["--cpp-path"]:
        return cpp_path_main()
    batch = wants_batch(argv)
    ap = build_parser(batch=batch)
    args = ap.parse_args(argv)
    if batch:  # never handed to coal-cpp: the batch driver runs every job in this process
        from .batch import main as batch_main
        return batch_main(ap, args, argv)
    if args.batch_out_dir:
        ap.error("--batch-out-dir needs --batch")
    if args.out == "-":
        args.out = None  # `-o -` is standard output (the default), as for most Unix tools
    if args.json and not args.out:
        ap.error("--json requires -o (the alignment and the report cannot both go to stdout)")
    if args.resume and not args.out:
        ap.error("--resume requires -o (only runs that write to a file can be resumed)")
    if args.max_warnings < 0:
        ap.error("--max-warnings must be 0 (no limit) or more")
    if stdin_count([args.pep_aln, *args.nuc_fasta]) > 1:
        ap.error(STDIN_TWICE)
    why = report_path_conflict(args.report, [args.pep_aln, *args.nuc_fasta], args.out) or report_folder_missing(args.report)
    if why:
        ap.error(why)
    if command_line:
        limit_blas_threads()  # this process and its workers are COAL's own
        rc = _maybe_cpp(ap, args, argv)
        if rc is not None:
            return rc
    install_signal_handlers()
    return run_job(args, argv)[0]
