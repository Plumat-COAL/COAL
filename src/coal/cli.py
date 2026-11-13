"""Command-line interface.

    coal PEP_ALN NUC_FASTA [NUC_FASTA ...] [-codontable N] [-nogap] [-nomismatch] [-o OUT]
    coal --batch JOBS [defaults for every job]          (coal/batch.py: many small runs in one process)
    coal kaks -i IN.axt -o OUT.kaks [-m YN ...]         (Ka/Ks, KaKs_Calculator format; run by coal-cpp)
    coal dnds IN [-o OUT.tsv] [-m YN ...]               (dN/dS, tidy TSV; run by coal-cpp)
    coal pipeline pairs|clusters ... -o DIR             (CDS in, dN/dS table out: coal/dnds_pipeline.py)
    coal translate / coal codon-dnds                    (the pipeline's native steps; run by coal-cpp)

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


# The epilog of `coal --help` (and of `coal-cpp --help`, which prints the same text): everything else `coal` does.
# One line per command; tests/test_cli_commands.py checks that every command main() dispatches is listed.
COMMANDS = (
    ("kaks", "Ka/Ks for pairs of aligned coding sequences, in KaKs_Calculator's format"),
    ("dnds", "dN/dS for pairs of aligned coding sequences, as a tidy table"),
    ("pipeline", "coding sequences in, dN/dS table out (pairs | clusters | install)"),
    ("translate", "coding sequences to proteins (a step of the pipeline)"),
    ("codon-dnds", "protein alignment, codon alignment and Ka/Ks for groups (the pipeline's engine)"),
    ("cite", "what to cite for what you ran"),
    ("help", "this help; `coal help COMMAND` is `coal COMMAND --help`"),
)
COMMANDS_EPILOG = (
    "Other commands (`coal COMMAND --help`):\n"
    + "".join(f"  coal {name:<11} {text}\n" for name, text in COMMANDS)
    + "  coal --cpp-path  where the C++ engine coal-cpp is (Ka/Ks and the pipeline need it)\n"
    "pipeline, cite and help are commands of `coal` (Python), not of coal-cpp.\n"
    "\n"
    "Every option in detail: README.md (https://github.com/plumat-coal/COAL).")
# What an install without coal-cpp is, and how to get one with it (the same words in `coal pipeline`)
NEEDS_CPP = ("this install has the Python engine only, which makes codon alignments (`coal PEP_ALN NUC_FASTA`); "
             "Ka/Ks (`coal kaks`, `coal dnds`) and `coal pipeline` need the C++ engine coal-cpp. It is inside the "
             "platform wheels (Linux, macOS, Windows): download coal_align-<version>-py3-none-<platform>.whl from "
             "https://github.com/plumat-coal/COAL/releases/latest and `python -m pip install` that file (the "
             "source package coal_src.tar.gz and the py3-none-any wheel do not have it); or build it from a "
             "checkout (docs/cpp-engine.md) and set COAL_CPP. Details: docs/install.md")
MAX_WARNINGS = 20  # default --max-warnings
# The same words in coal-cpp (main.cpp, std::bad_alloc)
OUT_OF_MEMORY = ("out of memory: this run needs more memory than is free on this machine; run with -j 1 (fewer worker processes), split the input, or use a machine with more memory")
BATCH_DESTS = ("batch", "batch_out_dir", "batch_parallel")  # the Batch group's options (coal only; coal-cpp rejects them)


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
        epilog=COMMANDS_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,   # (the epilog's list keeps its lines)
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
    bat.add_argument("--batch-parallel", metavar="N", type=int, default=None,  # None: one after the other
                     help="with --batch: N jobs at a time in worker processes (0 = one per CPU core; default: "
                          "one after the other), each job's messages printed together when it ends")
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


def run_job(args, argv, prefix: str = "", want_report: bool = False, close: bool = True) -> tuple[int, RunReport]:
    """One conversion from parsed arguments: the body of a single `coal` run, also used for every job of a batch.
    Returns (exit code, the run report); 130 and 141 mean interrupted / closed pipe (the caller stops).
    Every stderr line is `coal: <prefix><kind>: ...` (prefix "" for a single run, "[k/n PEP] " in a batch).
    want_report: keep the per-sequence entries although args has no --report/--json (the batch report wants them).
    close: close the report's entries file on the way out (close=False: the caller reads the report first)."""
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
    try:
        return _run_job(args, rec, say, warn, notice, note, reporting, hidden)
    finally:
        if close:
            rec.close()  # the entries file (also on SIGTERM / Ctrl-C between the run and the report)


def _run_job(args, rec, say, warn, notice, note, reporting, hidden) -> tuple[int, RunReport]:
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


def kaks_main(argv) -> int:
    """`coal kaks ...` / `coal dnds ...` (and `coal translate`, `coal codon-dnds`, the native steps of
    `coal pipeline`): the Ka/Ks engine is native code only (coal-cpp, docs/kaks.md), so the
    command line is handed to the coal-cpp of this version (the lookup of `coal --cpp-path`)."""
    from . import dispatch
    path = dispatch.find_cpp()
    ok, why = dispatch._version_ok(path) if path is not None else (False, "")
    wants_help = any(a in ("-h", "--help") or (a == "-H" and argv[0] == "kaks") for a in argv[1:])
    if wants_help and not ok:
        # the help does not need the engine: the copy of coal-cpp's text that ships with the package
        # (coal/_cpp_help.py, compared with the binary's by tests/test_cli_commands.py)
        from ._cpp_help import HELP
        try:
            sys.stdout.write(HELP[argv[0]])
            sys.stdout.flush()
        except BrokenPipeError:
            return closed_stdout_exit()
        print(f"coal {argv[0]}: note: this command cannot run here: "
              + (why if path is not None else f"coal-cpp was not found (looked at {dispatch.LOOKUP_ORDER}); "
                 + NEEDS_CPP), file=sys.stderr)
        return 0
    if path is None:
        print(f"coal {argv[0]}: needs coal-cpp, which was not found (looked at {dispatch.LOOKUP_ORDER}): "
              + NEEDS_CPP, file=sys.stderr)
        return 1
    if not ok:
        print(f"coal {argv[0]}: {why}", file=sys.stderr)
        return 1
    try:
        argv, tmp = _gunzip_inputs(argv)
    except (OSError, EOFError, ValueError) as e:
        print(f"coal {argv[0]}: error: {e}", file=sys.stderr)
        return 1
    if tmp is not None:   # coal-cpp as a child, so that the decompressed copies can be removed when it ends
        return _run_cpp_child(path, argv, tmp)
    rc = dispatch.run_cpp(path, argv)
    if isinstance(rc, str):
        print(f"coal {argv[0]}: {rc}", file=sys.stderr)
        return 1
    return rc


GZIP_MAGIC = b"\x1f\x8b"
# options of kaks / dnds / translate / codon-dnds whose value is a file the command writes (never decompressed)
_CPP_OUTPUT_OPTIONS = ("-o", "-O", "--out", "--genes", "--failures", "--codon-out", "--codon-dir", "--protein-out",
                       "--protein-dir", "--kaks-out")


def is_gzip(path) -> bool:
    """True when the file starts with the gzip magic bytes (whatever its name is)."""
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == GZIP_MAGIC
    except OSError:
        return False


def gunzip_to(src: str, dst: str) -> str:
    """Decompress the gzip file src into dst (streamed) and return the SHA-256 of the decompressed bytes (hex).
    ValueError with a readable reason when src is not a complete gzip file; OSError when dst cannot be written."""
    import gzip
    import hashlib
    import zlib
    h = hashlib.sha256()
    try:
        with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
            for block in iter(lambda: fin.read(1 << 20), b""):
                h.update(block)
                fout.write(block)
        return h.hexdigest()
    except (EOFError, zlib.error, gzip.BadGzipFile) as e:
        try:
            os.remove(dst)
        except OSError:
            pass
        raise ValueError(f"cannot decompress {src!r}: it starts like a gzip file but is damaged or cut short "
                         f"({e})") from None


_GUNZIPPED: dict[str, str] = {}   # decompressed copy -> the file as the user named it (for messages)


def _gunzip_inputs(argv) -> tuple[list[str], str | None]:
    """`coal kaks|dnds|translate|codon-dnds ...`: every input file that is gzip-compressed (by its first bytes) is
    decompressed into a temporary folder and the argument replaced by that copy, because coal-cpp reads plain text
    only. Returns (the command line to run, the temporary folder or None). Arguments after an output option
    (-o, --genes, ...) are never touched."""
    import tempfile
    out, tmp = [argv[0]], None

    def plain(path: str) -> str:
        nonlocal tmp
        if tmp is None:
            tmp = tempfile.mkdtemp(prefix="coal_gunzip_")
        d = os.path.join(tmp, str(len(os.listdir(tmp))))
        os.mkdir(d)
        name = os.path.basename(path)
        dst = os.path.join(d, name[:-3] if name.lower().endswith(".gz") and len(name) > 3 else name + ".txt")
        gunzip_to(path, dst)
        _GUNZIPPED[dst] = path
        return dst
    try:
        for k, arg in enumerate(argv[1:], 1):
            prev = argv[k - 1]
            if prev in _CPP_OUTPUT_OPTIONS or arg == "-":
                out.append(arg)
            elif arg.startswith("--") and "=" in arg:
                opt, val = arg.split("=", 1)
                out.append(f"{opt}={plain(val)}" if opt not in _CPP_OUTPUT_OPTIONS and is_gzip(val) else arg)
            elif arg[:2] in ("-i", "-I") and len(arg) > 2 and is_gzip(arg[2:].lstrip("=")):
                out.append(arg[:2] + plain(arg[2:].lstrip("=")))
            elif not arg.startswith("-") and is_gzip(arg):
                out.append(plain(arg))
            else:
                out.append(arg)
    except BaseException:
        if tmp is not None:
            import shutil
            shutil.rmtree(tmp, ignore_errors=True)
        raise
    return out, tmp


def _run_cpp_child(path: str, argv, tmp: str) -> int:
    """Run coal-cpp as a child process (the same standard streams), remove the folder tmp when it has ended, and
    return its exit status (128 + signal number when a signal ended it)."""
    import shutil
    import subprocess
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.flush()
        except (OSError, ValueError, AttributeError):
            pass
    proc = None
    old = None
    # coal-cpp was given decompressed copies in a temporary folder: its messages name those. They are passed on
    # line by line with the copies' paths replaced by the files the user named.
    names = {os.fsencode(copy): os.fsencode(orig) for copy, orig in _GUNZIPPED.items()
             if tmp and copy.startswith(tmp)}
    relay = None
    try:
        try:
            proc = subprocess.Popen([path, *argv], stderr=subprocess.PIPE if names else None)
        except OSError as e:
            print(f"coal {argv[0]}: could not start {path} ({e.strerror or e})", file=sys.stderr)
            return 1

        def forward(signum, frame):
            proc.terminate()
        try:
            old = signal.signal(signal.SIGTERM, forward)
        except (ValueError, OSError):
            old = None
        if names:
            import threading

            def pump():
                out = getattr(sys.stderr, "buffer", None)
                for line in proc.stderr:
                    for copy, orig in names.items():
                        line = line.replace(copy, orig)
                    try:
                        if out is not None:
                            out.write(line)
                            out.flush()
                        else:
                            sys.stderr.write(line.decode("utf-8", "backslashreplace"))
                    except (OSError, ValueError):
                        pass
            relay = threading.Thread(target=pump, daemon=True)
            relay.start()
        while True:
            try:
                rc = proc.wait()
                break
            except KeyboardInterrupt:   # Ctrl-C reaches coal-cpp too; it cleans up and exits 130
                continue
        if relay is not None:
            relay.join(timeout=30)
        return rc if rc >= 0 else 128 - rc
    finally:
        if old is not None:
            try:
                signal.signal(signal.SIGTERM, old)
            except (ValueError, OSError):
                pass
        shutil.rmtree(tmp, ignore_errors=True)


def main(argv=None) -> int:
    command_line = argv is None  # only the command line may hand the run to coal-cpp
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv == ["--cpp-path"]:
        return cpp_path_main()
    if argv and argv[0] == "help":   # `coal help` = `coal --help`; `coal help pipeline pairs` = `... --help`
        known = [name for name, _ in COMMANDS if name != "help"]
        if len(argv) > 1 and argv[1] not in known:
            print(f"coal help: {argv[1]!r} is not a command of coal: " + ", ".join(known)
                  + " (`coal help` alone: the codon-alignment command and the list)", file=sys.stderr)
            return 2
        argv = [*argv[1:], "--help"]
    if argv and argv[0] in ("kaks", "dnds", "translate", "codon-dnds"):
        return kaks_main(argv)
    if argv and argv[0] == "cite":   # what to cite for what was run (coal/citations.py)
        from .citations import main as cite_main
        try:
            rc = cite_main(argv[1:])
            sys.stdout.flush()
            return rc
        except BrokenPipeError:   # `coal cite | head`
            return closed_stdout_exit()
    if argv and argv[0] == "pipeline":  # CDS -> orthologs / clusters -> alignment -> codon alignment -> Ka/Ks
        from .dnds_pipeline import main as pipeline_main
        return pipeline_main(argv[1:])
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
