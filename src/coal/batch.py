"""Batch mode: `coal --batch JOBS` runs many small conversions in one process (GitHub issue #15).

A pipeline with thousands of gene families pays Python's start-up (~65 ms) once instead of once per job. Each job is
a COAL command line without the program name, one per line of JOBS; the options of the `coal --batch` command line
are the defaults of every job and a job line may override them. Every job runs in this process through the same
code path as a single run (cli.run_job), so its output is byte-identical to running it alone. Jobs run one after the
other; `-j` keeps its per-job meaning (cross-job parallelism is a follow-up). Batches are never handed to coal-cpp
(dispatch.py): coal-cpp has no batch mode, and the Python engine is where start-up is worth saving.

Exit status: 0 every job ok; 2 some job left sequences out (--keep-going) and none failed; 1 some job failed
(the batch goes on after a failed job); 2 before anything runs for a mistake in JOBS (the message names the line);
130 / 143 / 141 interrupted, as a single run.
"""
from __future__ import annotations

import argparse
import json
import os
import shlex
import sys

from . import __version__

SCHEMA = "coal.batch/1"
# --batch-out-dir: the output file's extension per -output format
EXTENSION = {"fasta": "fasta", "clustal": "aln", "paml": "paml", "codon": "codon"}
# options a job line may not carry: batch-level ones, and --resume (not supported in a batch yet)
NOT_IN_A_JOB = {"batch": "--batch", "batch_out_dir": "--batch-out-dir", "resume": "--resume",
                "json": "--json (print the batch report with --json on the coal --batch command line)"}


class JobsError(Exception):
    """A mistake in the JOBS file (usage error, exit 2): `line` is its 1-based line number, 0 = the file itself."""

    def __init__(self, line: int, message: str):
        super().__init__(message)
        self.line = line


class Job:
    def __init__(self, index: int, line_no: int, line: str, tokens: list[str]):
        self.index, self.line_no, self.line, self.tokens = index, line_no, line, tokens
        self.args: argparse.Namespace | None = None
        self.argv: list[str] = []  # the single-run command line this job equals (defaults, then the job line)
        self.output = ""


def read_jobs(path: str) -> list[tuple[int, str]]:
    """(line number, line) for every non-empty, non-comment line of JOBS ('-' = stdin)."""
    if path == "-":
        text = sys.stdin.read()
    else:
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
    out = []
    for n, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if line and not line.startswith("#"):
            out.append((n, line))
    return out


def _job_parser() -> argparse.ArgumentParser:
    """The single-run parser, raising JobsError instead of printing and exiting (also for -h / --version)."""
    from .cli import build_parser
    jp = build_parser()

    def error(message):
        raise JobsError(0, message)

    def exit_(status=0, message=None):
        raise JobsError(0, "-h and --version are not allowed on a job line")
    jp.error = error
    jp.exit = exit_
    jp._print_message = lambda *a, **k: None  # -h / --version print nothing before exit_ raises
    return jp


def _default_tokens(ap: argparse.ArgumentParser, argv: list[str]) -> list[str]:
    """The `coal --batch` command line without its batch-level options (--batch, --batch-out-dir, --report, --json):
    what every job's command line starts with."""
    drop = {"batch", "batch_out_dir", "report", "json"}
    out, skip = [], False
    for i, tok in enumerate(argv):
        if skip:
            skip = False
            continue
        if tok == "--":
            out.extend(argv[i:])
            break
        if tok.startswith("--"):
            name = tok.split("=", 1)[0]
            matches = [o for o in ap._option_string_actions if o.startswith(name)]
            act = ap._option_string_actions.get(name) or (ap._option_string_actions[matches[0]]
                                                            if len(matches) == 1 else None)
            if act is not None and act.dest in drop:
                skip = act.nargs != 0 and "=" not in tok
                continue
        out.append(tok)
    return out


def plan(ap: argparse.ArgumentParser, args: argparse.Namespace, argv: list[str], lines) -> list[Job]:
    """Parse and check every job before anything runs. Raises JobsError (with the line number)."""
    if args.pep_aln is not None or args.nuc_fasta:
        raise JobsError(0, "--batch takes the input files from the JOBS file; give none on the command line")
    if args.resume:
        raise JobsError(0, "--resume is not supported with --batch yet")
    if args.max_warnings < 0:
        raise JobsError(0, "--max-warnings must be 0 (no limit) or more")
    if not lines:
        raise JobsError(0, f"no jobs in {args.batch} (one coal command line per line; # starts a comment)")
    jp = _job_parser()
    defaults = _default_tokens(ap, argv)
    jobs = []
    for k, (n, line) in enumerate(lines, 1):
        job = Job(k, n, line, [])
        try:
            job.tokens = shlex.split(line, posix=True)
        except ValueError as e:
            raise JobsError(n, f"{e} (quotes and backslashes follow POSIX shell rules)") from None
        # the job's namespace starts as the defaults of the main command line; the job line's options come last
        ns = argparse.Namespace(**vars(args))
        for dest in ("pep_aln", "nuc_fasta", *NOT_IN_A_JOB, "report"):
            delattr(ns, dest)
        try:
            jp.parse_args(job.tokens, namespace=ns)
        except JobsError as e:
            raise JobsError(n, str(e)) from None
        for dest, name in NOT_IN_A_JOB.items():
            if getattr(ns, dest):
                raise JobsError(n, f"{name} is not allowed on a job line")
        if "-" in (ns.pep_aln, *ns.nuc_fasta):
            raise JobsError(n, "'-' (standard input) cannot be a job's input: standard input can be read once, and "
                               "`--batch -` may be reading the job list from it; give a file path")
        if not ns.out:
            if not args.batch_out_dir:
                raise JobsError(n, "the job has no -o (every job of a batch writes to a file; or give "
                                   "--batch-out-dir DIR on the coal --batch command line)")
            stem = os.path.splitext(os.path.basename(ns.pep_aln))[0] or os.path.basename(ns.pep_aln)
            ns.out = os.path.join(args.batch_out_dir, f"{stem}.{EXTENSION[ns.output_format]}")
            job.argv = [*defaults, *job.tokens, "-o", ns.out]
        else:
            job.argv = [*defaults, *job.tokens]
        job.args, job.output = ns, ns.out
        jobs.append(job)
    seen: dict[str, Job] = {}
    for job in jobs:
        key = os.path.normcase(os.path.abspath(job.output))
        if key in seen:
            raise JobsError(job.line_no, f"writes {job.output}, the output of the job on line "
                                         f"{seen[key].line_no} (give each job its own -o)")
        seen[key] = job
    return jobs


def write_json(path: str, data: dict) -> None:
    """Write atomically (temp file + rename); never leaves the temp file behind."""
    tmp = f"{path}.tmp{os.getpid()}"
    try:
        with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(data, fh, indent=1)
            fh.write("\n")
        os.replace(tmp, path)
    except OSError:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def main(ap: argparse.ArgumentParser, args: argparse.Namespace, argv: list[str]) -> int:
    """Run a `coal --batch` command line (parsed by cli.main with build_parser(batch=True))."""
    from .cli import install_signal_handlers, run_job
    from .core import limit_blas_threads
    try:
        lines = read_jobs(args.batch)
    except KeyboardInterrupt:  # Ctrl-C while `--batch -` waits for stdin
        print("coal: interrupted", file=sys.stderr)
        return 130
    except OSError as e:
        ap.error(f"cannot read the JOBS file {args.batch}: {e.strerror or e}")
    except UnicodeDecodeError:
        ap.error(f"cannot read the JOBS file {args.batch}: it is not UTF-8 text")
    try:
        jobs = plan(ap, args, argv, lines)
    except JobsError as e:
        ap.error(f"{'standard input' if args.batch == '-' else args.batch} line {e.line}: {e}" if e.line else str(e))
    if args.batch_out_dir:
        try:
            os.makedirs(args.batch_out_dir, exist_ok=True)
        except OSError as e:
            ap.error(f"cannot create --batch-out-dir {args.batch_out_dir}: {e.strerror or e}")
    reporting = bool(args.report or args.json)
    limit_blas_threads()  # this process and its workers are COAL's own
    install_signal_handlers()
    n = len(jobs)
    done, ok, failed, skipped = [], 0, 0, 0
    for job in jobs:
        prefix = f"[{job.index}/{n} {job.args.pep_aln}] "
        if args.progress:
            print(f"coal: progress: job {job.index}/{n} {job.args.pep_aln}", file=sys.stderr, flush=True)
        code, rec = run_job(job.args, job.argv, prefix=prefix, want_report=reporting)
        if code in (130, 141):  # interrupted / closed pipe: stop here, as a single run does
            return code
        entry = {"index": job.index, "line": job.line, "exit_code": code, "output": job.output}
        if reporting:
            entry["report"] = rec.to_dict()
        done.append(entry)
        if code == 1:
            failed += 1
        else:
            ok += 1
        if code == 2:
            left = (rec.output or {}).get("sequences")
            total = len(rec.ids) if rec.ids is not None else None
            skipped += (total - left) if (left is not None and total is not None) else 0
    partial = sum(1 for e in done if e["exit_code"] == 2)
    rc = 1 if failed else (2 if partial else 0)
    if not args.quiet:
        extra = f", {partial} with sequences left out" if partial else ""
        print(f"coal: batch: {n} jobs, {ok} ok{extra}, {failed} failed", file=sys.stderr)
    if reporting:
        report = {"schema": SCHEMA, "version": __version__,
                  "command_line": " ".join(shlex.quote(a) for a in ["coal", *argv]), "jobs": done,
                  "summary": {"jobs": n, "ok": ok, "failed": failed, "skipped_sequences": skipped}}
        if args.report:
            try:
                write_json(args.report, report)
            except OSError as e:
                why = "it is a folder, not a file" if os.path.isdir(args.report) else (e.strerror or e)
                print(f"coal: error: could not write the report {args.report}: {why}", file=sys.stderr)
                rc = 1
        if args.json:
            print(json.dumps(report, indent=1))
    return rc
