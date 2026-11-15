#!/usr/bin/env python3
"""Run the tutorials of docs/tutorials/ as tests: their commands, in order, on the real NCBI data, and compare.

The tutorials show real commands and their real outputs. This script reads each tutorial as it is written (no
markers needed), runs its shell blocks in a scratch folder, and checks two things:

1. **Outputs.** Every line of every output block shown after a command block must appear in what that block printed
   (stdout and stderr together), after a few normalisations that only remove what cannot be the same on another
   machine or day: run times (``6.0 s``), timestamps, the scratch folder's path (the docs write ``.../``), COAL's
   own version string, and ``ok (already here)`` from a download script that found its file in the cache. ``...``
   in a documented line matches anything; a line that is only ``...`` is skipped. ``exit N`` is checked against
   the block's exit status. A block whose commands fail, where the documented output shows no ``exit N`` other
   than 0, fails too.
2. **Facts in the prose.** A short list (``FACTS`` below) of numbers the prose states (rows of each output table,
   alignment lengths, cluster sizes, ...): each takes its number *from the tutorial's text* (so the check follows the
   docs when they change) and recomputes it from the files the run wrote.

What is not run, and why, is listed per tutorial in ``SKIPS`` and printed in the report: the clean-up sections
(later tutorials use the earlier ones' results), ``export T=...`` (T is set here), PAL2NAL itself when
``pal2nal.pl`` is not available, tutorial 3's MAFFT run (6 minutes; ``--full`` runs it) and tutorial 4 (R).

Usage (from a checkout, with coal-cpp found by coal (COAL_CPP), DIAMOND 2.2.8 and Clustal Omega 1.2.4):

    python scripts/run_tutorials.py                       # tutorials 01 02 03 05, scratch folder in /tmp
    python scripts/run_tutorials.py -t 01 -t 02 --data-cache DIR --jobs 4 --report report.json

``--data-cache DIR`` holds the downloaded ``*.cds.fna.gz`` (filled by the tutorials' own download scripts, which
check every file's SHA-256); ``--jobs N`` runs tutorial 3's and 5's commands with ``-j N`` instead of ``-j 1``
(their outputs do not depend on it; tutorials 1 and 2 print their command lines, so they always run as written).
Exit status: 0 when everything matched, 1 on any difference (each one is listed), 2 on a usage or setup error.
Standard library only.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
import time
from dataclasses import dataclass, field

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TUTORIALS_DIR = os.path.join(ROOT, "docs", "tutorials")
TUTORIALS = {
    "01": "01-codon-alignment.md",
    "02": "02-pairs-dnds.md",
    "03": "03-clusters.md",
    "04": "04-r-coalr.md",
    "05": "05-choosing-a-method.md",
}
RUNNABLE = ["01", "02", "03", "05"]  # 04 is R (coalr, orthologr): not run here
DEPENDS = {"05": ["01", "02"]}        # tutorial 5 reads t01/ and ecsa_filtered/
SCRATCH = "coal-tutorials"             # the scratch folder the tutorials create and cd into
DATA = "coal-tutorial-data"            # the download scripts' default destination, inside SCRATCH
JOBS_SUBSTITUTION = {"03", "05"}       # tutorials whose commands may run with -j N (their outputs do not show it)

# Commands that start a shell block. A bare ``` block whose first line starts with one of these (or with
# VAR=value, or with "$ ") is a command block; anything else is output. Output lines that start like a command
# ("coal pipeline: ...", "coal: error: ...") have a colon right after the command name(s) and are excluded.
COMMANDS = {
    "awk", "bash", "cat", "cd", "clustalo", "cmp", "coal", "cut", "diff", "echo", "export", "for", "git", "grep",
    "gunzip", "gzcat", "head", "ls", "mafft", "mkdir", "pal2nal", "pal2nal.pl", "python", "python3", "R", "rm",
    "Rscript", "sed", "sh", "sha256sum", "shasum", "sort", "tail", "wc", "while", "zcat",
}
SHELL_LANGS = {"bash", "sh", "shell", "console", "shell-session"}
OUTPUT_LANGS = {"text", "txt", "output", "tsv", "plain"}

KNOWN_DRIFT = os.path.join(ROOT, "tests", "fixtures", "tutorials", "known_drift.json")

# What is not run. Each rule: (tutorial, a regex searched in the block's text, the reason, condition key or None).
# Condition keys: "no_pal2nal" (pal2nal.pl not found), "not_full" (unless --full).
SKIPS = [
    ("01", r"pal2nal\.pl", "PAL2NAL v14 (pal2nal.pl) not found: give --pal2nal or set PAL2NAL", "no_pal2nal"),
    ("03", r"--aligner mafft", "the MAFFT run takes about 6 minutes (run it with --full)", "not_full"),
    ("03", r"chlamydia_mafft", "needs the MAFFT run (run it with --full)", "not_full"),
]


@dataclass
class Block:
    """One fenced block of a tutorial."""
    tutorial: str
    index: int          # 0-based, in order of the file
    line: int           # 1-based line of the opening fence
    heading: str        # the last ## heading above it
    lang: str
    text: str           # dedented content
    kind: str = ""      # "command", "output", "prompt" ($ cmd + output in one block), "other" (r, ...)
    outputs: list = field(default_factory=list)  # output blocks shown for a command block
    skip: str = ""      # why it is not run

    @property
    def where(self) -> str:
        return f"{TUTORIALS.get(self.tutorial, self.tutorial)}:{self.line} ({self.heading})"


FENCE = re.compile(r"^( *)(```+)\s*([\w+-]*)\s*$")


def parse_blocks(markdown: str, tutorial: str = "") -> list:
    """The fenced blocks of a Markdown text, with the ## heading each is under. Fences may be indented (inside a
    list item); the content is dedented by the fence's indent."""
    blocks, heading = [], ""
    lines = markdown.splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("## "):
            heading = ln[3:].strip()
        m = FENCE.match(ln)
        if not m:
            i += 1
            continue
        indent, ticks, lang = len(m.group(1)), m.group(2), m.group(3).lower()
        start, body = i + 1, []
        i += 1
        while i < len(lines) and not re.match(rf"^ {{0,{indent + 3}}}{ticks}\s*$", lines[i]):
            body.append(lines[i][indent:] if lines[i][:indent].strip() == "" else lines[i])
            i += 1
        blocks.append(Block(tutorial, len(blocks), start, heading, lang, "\n".join(body)))
        i += 1
    return blocks


def _first_line(text: str) -> str:
    for ln in text.splitlines():
        if ln.strip():
            return ln.strip()
    return ""


def classify(block: Block) -> str:
    if block.lang in SHELL_LANGS:
        return "prompt" if _first_line(block.text).startswith("$ ") else "command"
    if block.lang in OUTPUT_LANGS:
        return "output"
    if block.lang:
        return "other"  # r, python, json, ...: not run
    first = _first_line(block.text)
    if not first:
        return "output"
    if first.startswith("$ "):
        return "prompt"
    if re.match(r"^[A-Za-z_][A-Za-z0-9_]*=", first):
        return "command"
    if re.match(r"^coal( [\w-]+)*:", first) or first.startswith("["):
        return "output"
    return "command" if first.split()[0] in COMMANDS else "output"


def plan(blocks: list) -> list:
    """Classify the blocks and attach each output block to the command block before it (output blocks shown one
    after another all belong to that command block). Returns the blocks to run (command and prompt blocks)."""
    runnable, last = [], None
    for b in blocks:
        b.kind = classify(b)
        if b.kind in ("command", "prompt"):
            runnable.append(b)
            last = b
        elif b.kind == "output":
            if last is not None:
                last.outputs.append(b)
        else:
            last = None  # an R or other block in between: what follows is not the shell block's output
    return runnable


def split_prompt(text: str):
    """A "$ cmd" block: the commands (without "$ ") and the lines shown as their output."""
    cmds, out = [], []
    for ln in text.splitlines():
        (cmds.append(ln[2:]) if ln.startswith("$ ") else out.append(ln))
    return "\n".join(cmds), "\n".join(out)


def script_of(block: Block, jobs: int | None) -> str:
    """The shell text that is run for a command block."""
    text = split_prompt(block.text)[0] if block.kind == "prompt" else block.text
    # T is set by the runner; the tutorials' `export T=/path/to/COAL/docs/tutorials` is a placeholder
    text = "\n".join(("# (set by the runner) " + ln) if re.match(r"^\s*export T=", ln) else ln
                     for ln in text.splitlines())
    if jobs and block.tutorial in JOBS_SUBSTITUTION:
        text = re.sub(r"(?<=\s)-j 1(?=\s|$)", f"-j {jobs}", text)
    return text + "\n"


def documented_lines(block: Block) -> list:
    """The documented output lines of a command block, in order."""
    lines = []
    if block.kind == "prompt":
        lines += split_prompt(block.text)[1].splitlines()
    for o in block.outputs:
        lines += o.text.splitlines()
    return [ln for ln in lines if ln.strip()]


# --- normalisation and matching -------------------------------------------------------------------------------

def normalise(line: str, coal_version: str = "", doc_version: str = "", scratch: str = "") -> str:
    s = line.replace("\r", "")
    if scratch:
        s = s.replace(scratch + "/", ".../").replace(scratch, "...")
    # the longer version first: a release (0.6.0) checked against a page run with its development version
    # (0.6.0.dev0) would otherwise turn the page's "0.6.0.dev0" into "<COAL_VERSION>.dev0"
    for v in sorted({coal_version, doc_version} - {""}, key=len, reverse=True):
        s = s.replace(v, "<COAL_VERSION>")
    s = re.sub(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(\.\d+)?(Z|[+-]\d\d:?\d\d)?", "<TIME>", s)
    s = re.sub(r"\b\d+\.\d+ s\b", "<T> s", s)                  # stage timings: "diamond_search 6.0 s"
    s = re.sub(r"^ok \(already here\)\s+", "ok ", s.strip())   # the download scripts, with a cached file
    return re.sub(r"\s+", " ", s).strip()


def line_pattern(doc_line: str) -> re.Pattern:
    """A documented line as a regex: "..." matches anything (paths shortened as ".../x", elided text)."""
    parts = doc_line.split("...")
    return re.compile("^" + ".*".join(re.escape(p) for p in parts) + "$")


def compare(doc: list, actual: list, status: int, **norm) -> list:
    """Documented lines (after normalisation) that are not in the actual output. Each actual line is used at most
    once; the order is not checked (stdout and stderr of a pipe interleave differently from run to run).
    "exit N" is satisfied by an "exit N" line or by the block's exit status. Returns the missing lines."""
    pool = [normalise(a, **norm) for a in actual]
    pool = [a for a in pool if a]
    missing = []
    for d in doc:
        nd = normalise(d, **norm)
        if not nd or nd == "...":
            continue
        m = re.fullmatch(r"exit (\d+)", nd)
        pat = line_pattern(nd)
        hit = next((k for k, a in enumerate(pool) if a == nd or pat.match(a)), None)
        if hit is not None:
            pool.pop(hit)
        elif m and int(m.group(1)) == status:
            continue
        elif len(nd.split()) > 1 and _consume_words(nd.split(), pool):
            continue  # `ls` on a terminal prints columns, into a pipe one name per line
        else:
            missing.append(d)
    return missing


def documented_exits(block: Block) -> set:
    """The exit statuses the documented output shows ("exit N")."""
    return {int(m.group(1)) for ln in documented_lines(block) for m in [re.fullmatch(r"\s*exit (\d+)\s*", ln)] if m}


def _consume_words(words: list, pool: list) -> bool:
    """When every word of a documented line is a line of its own in the output (and looks like a file name), use
    those lines up and say so."""
    if not all(re.fullmatch(r"[\w.+-]+/?", w) for w in words):
        return False
    idx = []
    for w in words:
        k = next((k for k, a in enumerate(pool) if a == w and k not in idx), None)
        if k is None:
            return False
        idx.append(k)
    for k in sorted(idx, reverse=True):
        pool.pop(k)
    return True


def expected_exit(block: Block) -> set:
    """The exit statuses the block itself may end with: 0 when it prints them (`echo "exit $?"`) or shows none,
    else the documented ones (a command shown with "exit 1" under it, as tutorial 1's --strict run)."""
    if "$?" in block.text:
        return {0}
    return documented_exits(block) or {0}


def failures_allowed(block: Block) -> bool:
    """Commands of the block may fail when the documented output shows a non-zero exit status."""
    return any(c != 0 for c in documented_exits(block))


# --- facts stated in the prose ----------------------------------------------------------------------------------
# Each fact: (tutorial, regex over the tutorial's text with one group that captures the stated number, a shell
# command run in the scratch folder after the tutorial that prints the actual number, a short description).
# Numbers are compared after removing thousands separators. A regex that no longer matches the text is reported
# too (the prose changed: update the fact).
FACTS = [
    ("01", r"Each sequence is now ([\d,]+) columns long",
     "awk '/^>/{if (n) exit; n=1; next} {L += length($0)} END {print L}' t01/glpK.codon.fasta",
     "columns of the codon alignment"),
    ("01", r"`-nogap` leaves ([\d,]+) columns",
     "awk '/^>/{if (n) exit; n=1; next} {L += length($0)} END {print L}' t01/glpK.nogap.fasta",
     "columns left by -nogap"),
    ("01", r"`wrong_code\.fasta` is\s+byte-identical to `glpK\.codon\.fasta` \(()",
     "cmp -s t01/wrong_code.fasta t01/glpK.codon.fasta && echo", "wrong_code.fasta == glpK.codon.fasta"),
    ("02", r"\| `dnds\.tsv` \| ([\d,]+) \|", "tail -n +2 ecoli_vs_salmonella/dnds.tsv | wc -l", "rows of dnds.tsv"),
    ("02", r"\| `orthologs\.tsv` \| ([\d,]+) \|", "tail -n +2 ecoli_vs_salmonella/orthologs.tsv | wc -l",
     "rows of orthologs.tsv"),
    ("02", r"\| `unpaired\.tsv` \| ([\d,]+) \|", "tail -n +2 ecoli_vs_salmonella/unpaired.tsv | wc -l",
     "rows of unpaired.tsv"),
    ("02", r"\| `genes\.tsv` \| ([\d,]+) \|", "tail -n +2 ecoli_vs_salmonella/genes.tsv | wc -l", "rows of genes.tsv"),
    ("02", r"\| `failed\.tsv` \| ([\d,]+) \|", "tail -n +2 ecoli_vs_salmonella/failed.tsv | wc -l",
     "rows of failed.tsv"),
    ("02", r"`no_hit` \(([\d,]+) \*E\. coli\*",
     "awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h[\"reason\"]==\"no_hit\" && $h[\"set\"]==\"query\"' "
     "ecoli_vs_salmonella/unpaired.tsv | wc -l", "no_hit E. coli genes"),
    ("02", r"`codon_alignments\.fasta` has the ([\d,]+) codon alignments",
     "echo $(( $(grep -c '^>' ecsa_filtered/codon_alignments.fasta) / 2 ))", "codon alignments kept"),
    ("02", r"`genes\.tsv` lists the ([\d,]+) dropped records",
     "awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h[\"status\"]==\"dropped\"' ecsa_filtered/genes.tsv | wc -l",
     "records dropped by --longest-isoform"),
    ("03", r"([\d,]+) of the [\d,]+ clusters have exactly 4 members",
     "awk -F'\\t' 'NR>1 && $3==4' chlamydia/cluster_summary.tsv | wc -l", "clusters of 4"),
    ("03", r"[\d,]+ of the ([\d,]+) clusters have exactly 4 members",
     "tail -n +2 chlamydia/cluster_summary.tsv | wc -l", "clusters"),
    ("03", r"\(Pmp\) autotransporter family .*?: ([\d,]+) pairs",
     "awk -F'\\t' 'NR>1 && $3==40 {print $8}' chlamydia/cluster_summary.tsv", "pairs of the 40-member cluster"),
    ("03", r"\(Pmp\) autotransporter family .*?: [\d,]+ pairs, ([\d,]+) of them saturated",
     "awk -F'\\t' 'NR>1 && $3==40 {print $11}' chlamydia/cluster_summary.tsv", "saturated pairs of that cluster"),
]


def load_known_drift(path: str = KNOWN_DRIFT) -> list:
    if not path or not os.path.isfile(path):
        return []
    with open(path, encoding="utf-8") as fh:
        return json.load(fh).get("entries", [])


def split_known(tutorial: str, missing: list, known: list, used: set) -> tuple:
    """(new, known) differences: a documented line is known drift when an entry of `known` for this tutorial is a
    prefix of it (whitespace collapsed). Adds the indices of the entries used to `used`."""
    new, old = [], []
    for m in missing:
        nm = re.sub(r"\s+", " ", m).strip()
        k = next((i for i, e in enumerate(known) if e["tutorial"] == tutorial
                  and nm.startswith(re.sub(r"\s+", " ", e["line"]).strip())), None)
        if k is None:
            new.append(m)
        else:
            used.add(k)
            old.append(m)
    return new, old


def check_facts(tutorial: str, doc_text: str, scratch: str) -> list:
    """[(description, stated, actual, ok)] for the FACTS of one tutorial."""
    out = []
    flat = re.sub(r"\s+", " ", doc_text)
    for t, rx, cmd, what in FACTS:
        if t != tutorial:
            continue
        m = re.search(rx, flat, re.S)
        if not m:
            out.append((what, "(not found in the text)", "", False))
            continue
        stated = m.group(1).replace(",", "")
        r = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], cwd=scratch, capture_output=True, text=True)
        actual = r.stdout.strip() if r.returncode == 0 else f"(failed: {r.stderr.strip()[:200]})"
        out.append((what, stated or "(identical)", actual or "(identical)", r.returncode == 0 and actual == stated))
    return out


# --- running ----------------------------------------------------------------------------------------------------

def doc_coal_version(doc_text: str) -> str:
    m = re.search(r"Tested with \| COAL (\S+)", doc_text)
    return m.group(1) if m else ""


def installed_coal_version() -> str:
    r = subprocess.run([sys.executable, "-c", "import coal; print(coal.__version__)"], capture_output=True, text=True)
    return r.stdout.strip()


def apply_skips(tutorial: str, blocks: list, conditions: set):
    for b in blocks:
        if b.heading.lower().startswith("clean up"):
            b.skip = "clean-up step (later tutorials use these results)"
            continue
        for t, rx, why, cond in SKIPS:
            if t == tutorial and re.search(rx, b.text) and (cond is None or cond in conditions):
                b.skip = why
                break


def run_tutorial(tutorial: str, blocks: list, root: str, out_dir: str, env: dict, jobs: int | None,
                 timeout: float) -> dict:
    """Run the blocks of one tutorial in one bash session (so cd and variables carry over), from `root` (the folder
    that holds SCRATCH) or from SCRATCH itself for a tutorial that never changes into it. Writes out_dir/blk_N.sh and reads back out_N (output), st_N (status), err_N (commands
    that failed). Returns {index: (status, output_lines, failed_commands)}."""
    os.makedirs(out_dir, exist_ok=True)
    # Tutorials 1 to 3 start with `cd coal-tutorials` (1 creates it); the others are run from inside it
    start = os.path.join(root, SCRATCH)
    if any(re.search(rf"\bcd {SCRATCH}\b", b.text) for b in blocks) or not os.path.isdir(start):
        start = root
    driver = ["set -o pipefail", f"cd {shq(start)}"]
    for b in blocks:
        if b.skip:
            continue
        n = b.index
        with open(os.path.join(out_dir, f"blk_{n}.sh"), "w", encoding="utf-8") as fh:
            fh.write(script_of(b, jobs))
        o = shq(out_dir)
        driver += [
            f"trap 'echo \"$? $BASH_COMMAND\" >> {o}/err_{n}' ERR",
            f"source {o}/blk_{n}.sh > {o}/out_{n} 2>&1 < /dev/null",
            f"echo $? > {o}/st_{n}",
            "trap - ERR",
        ]
    with open(os.path.join(out_dir, "driver.sh"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(driver) + "\n")
    t0 = time.time()
    try:
        subprocess.run(["bash", os.path.join(out_dir, "driver.sh")], env=env, timeout=timeout, check=False)
    except subprocess.TimeoutExpired:
        pass  # the blocks that did not finish have no st_N: reported as such
    results = {}
    for b in blocks:
        if b.skip:
            continue
        n = b.index
        st = _read(os.path.join(out_dir, f"st_{n}"))
        out = _read(os.path.join(out_dir, f"out_{n}"))
        errs = [e for e in _read(os.path.join(out_dir, f"err_{n}")).splitlines() if e and " source " not in e
                and not e.split(" ", 1)[-1].startswith("source ")]
        results[n] = (int(st) if st.strip() else None, out.splitlines(), errs)
    results["seconds"] = time.time() - t0
    return results


def shq(s: str) -> str:
    return "'" + s.replace("'", "'\\''") + "'"


def _read(p: str) -> str:
    try:
        with open(p, encoding="utf-8", errors="replace") as fh:
            return fh.read()
    except FileNotFoundError:
        return ""


def find_pal2nal(arg: str | None) -> str | None:
    for c in (arg, os.environ.get("PAL2NAL"), shutil.which("pal2nal.pl"),
              os.path.join(ROOT, "PAL2NAL", "pal2nal.v14", "pal2nal.pl")):
        if c and os.path.isfile(c):
            return os.path.abspath(c)
    return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-t", "--tutorial", action="append", choices=RUNNABLE,
                    help="tutorial to run (repeatable; default: all of %s; 05 also runs 01 and 02)" % " ".join(RUNNABLE))
    ap.add_argument("--workdir", help="folder to run in (default: a new temporary folder, removed unless --keep)")
    ap.add_argument("--data-cache", help="folder with the downloaded *.cds.fna.gz (copied into the scratch folder; "
                                         "the download scripts still check every SHA-256)")
    ap.add_argument("--jobs", type=int, help="run tutorials 3 and 5 with -j N instead of -j 1")
    ap.add_argument("--pal2nal", help="pal2nal.pl of PAL2NAL v14 (default: $PAL2NAL, PATH, PAL2NAL/ of the checkout)")
    ap.add_argument("--full", action="store_true", help="also run the slow steps (tutorial 3's MAFFT run)")
    ap.add_argument("--keep", action="store_true", help="keep the temporary folder")
    ap.add_argument("--timeout", type=float, default=1800, help="seconds per tutorial (default 1800)")
    ap.add_argument("--report", help="write the results as JSON here")
    ap.add_argument("--summary", help="append a Markdown summary here (e.g. $GITHUB_STEP_SUMMARY)")
    ap.add_argument("--known-drift", default=KNOWN_DRIFT,
                    help="JSON of documented lines known to differ, reported but not failing (default: %(default)s; "
                         "'' for none)")
    ap.add_argument("--strict-known", action="store_true",
                    help="fail on known drift too, and on known-drift entries that no longer match anything")
    ap.add_argument("--list", action="store_true", help="only print what would be run and checked")
    a = ap.parse_args(argv)

    wanted = a.tutorial or RUNNABLE
    order = [t for t in RUNNABLE if t in wanted or any(t in DEPENDS.get(w, []) for w in wanted)]
    pal2nal = find_pal2nal(a.pal2nal)
    conditions = set()
    if not pal2nal:
        conditions.add("no_pal2nal")
    if not a.full:
        conditions.add("not_full")

    texts, plans = {}, {}
    for t in order:
        with open(os.path.join(TUTORIALS_DIR, TUTORIALS[t]), encoding="utf-8") as fh:
            texts[t] = fh.read()
        blocks = plan(parse_blocks(texts[t], t))
        apply_skips(t, blocks, conditions)
        plans[t] = blocks

    if a.list:
        for t in order:
            for b in plans[t]:
                state = f"SKIP ({b.skip})" if b.skip else f"run, {len(documented_lines(b))} output line(s) checked"
                print(f"{b.where}: {state}")
                print(textwrap.indent(script_of(b, a.jobs).rstrip(), "    | "))
            for _, rx, _, what in (f for f in FACTS if f[0] == t):
                print(f"{TUTORIALS[t]}: fact: {what}")
        return 0

    for tool in ("bash", "zcat", "awk", "sha256sum"):
        if not shutil.which(tool):
            print(f"run_tutorials: error: {tool} not found", file=sys.stderr)
            return 2

    workdir = a.workdir or tempfile.mkdtemp(prefix="coal-tutorials-run-")
    os.makedirs(workdir, exist_ok=True)
    root = os.path.join(workdir, "run")
    if os.path.exists(root):
        shutil.rmtree(root)
    scratch = os.path.join(root, SCRATCH)
    os.makedirs(os.path.join(scratch, DATA))
    if a.data_cache:
        for f in sorted(os.listdir(a.data_cache)):
            if f.endswith(".fna.gz"):
                shutil.copy2(os.path.join(a.data_cache, f), os.path.join(scratch, DATA, f))
    bindir = os.path.join(workdir, "bin")
    os.makedirs(bindir, exist_ok=True)
    env = dict(os.environ)
    env["T"] = TUTORIALS_DIR
    env["LC_ALL"] = "C"  # sort / uniq -c in the same order as the docs
    if pal2nal:
        link = os.path.join(bindir, "pal2nal.pl")
        if not os.path.exists(link):
            os.symlink(pal2nal, link)
    env["PATH"] = bindir + os.pathsep + os.path.dirname(os.path.abspath(sys.executable)) + os.pathsep + env["PATH"]
    coal_version = installed_coal_version()

    known, used = load_known_drift(a.known_drift), set()
    report, failures, t_start = {"tutorials": {}}, 0, time.time()
    for t in order:
        res = run_tutorial(t, plans[t], root, os.path.join(workdir, "logs", t), env, a.jobs, a.timeout)
        norm = dict(coal_version=coal_version, doc_version=doc_coal_version(texts[t]), scratch=scratch)
        rows = []
        for b in plans[t]:
            entry = {"where": b.where, "skip": b.skip}
            if not b.skip:
                st, out, errs = res[b.index]
                doc = documented_lines(b)
                missing = compare(doc, out, st if st is not None else -1, **norm) if st is not None else doc
                missing, known_here = split_known(t, missing, known, used) if not a.strict_known else (missing, [])
                problems = []
                if st is None:
                    problems.append("did not finish (timeout, or an earlier block ended the shell)")
                if st is not None and st not in expected_exit(b):
                    problems.append(f"exit status {st} (expected {sorted(expected_exit(b))})"
                                    + ("; failed: " + "; ".join(errs[:3]) if errs else ""))
                elif errs and not failures_allowed(b):
                    problems.append("a command failed (the documented output shows no error): " + "; ".join(errs[:3]))
                if missing:
                    problems.append(f"{len(missing)} of {len(doc)} documented line(s) not in the output")
                entry.update(status=st, checked=len(doc), missing=missing, known=known_here, problems=problems,
                             output=out[-80:])
                failures += bool(problems)
            rows.append(entry)
        facts = check_facts(t, texts[t], scratch)
        failures += sum(1 for f in facts if not f[3])
        report["tutorials"][t] = {"file": TUTORIALS[t], "seconds": round(res["seconds"], 1), "blocks": rows,
                                  "facts": [dict(what=w, stated=s, actual=x, ok=o) for w, s, x, o in facts]}
    stale = [e for i, e in enumerate(known) if i not in used and e["tutorial"] in order]
    report["stale_known_drift"] = stale
    if a.strict_known:
        failures += len(stale)
    report["seconds"] = round(time.time() - t_start, 1)
    report["coal_version"] = coal_version
    report["failures"] = failures

    lines = _render(report)
    print("\n".join(lines))
    if os.environ.get("GITHUB_ACTIONS") == "true":
        _annotate(report, plans)
        _annotate_stale(report)
    if a.report:
        with open(a.report, "w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=1)
    if a.summary:
        with open(a.summary, "a", encoding="utf-8") as fh:
            fh.write("\n".join(_render(report, markdown=True)) + "\n")
    if not a.workdir and not a.keep:
        shutil.rmtree(workdir, ignore_errors=True)
    else:
        print(f"run_tutorials: outputs and logs in {workdir}")
    return 1 if failures else 0


def _annotate(report: dict, plans: dict):
    """One GitHub Actions error annotation per difference, at the tutorial's line (readable from the checks API)."""
    esc = lambda x: x.replace("%", "%25").replace("\r", "").replace("\n", "%0A")
    for t, r in report["tutorials"].items():
        path = f"docs/tutorials/{r['file']}"
        for b, entry in zip(plans[t], r["blocks"]):
            if entry.get("problems"):
                body = "; ".join(entry["problems"]) + "".join("\nnot printed: " + m for m in entry["missing"][:5])
                print(f"::error file={path},line={b.line},title=tutorial drift::{esc(body)}")
        for b, entry in zip(plans[t], r["blocks"]):
            for m in entry.get("known") or []:
                print(f"::warning file={path},line={b.line},title=known tutorial drift::{esc(m[:300])}")
        for f in r["facts"]:
            if not f["ok"]:
                print(f"::error file={path},title=tutorial drift::{esc(f['what'] + ': the text says ' + f['stated'] + ', the run gives ' + f['actual'])}")


def _annotate_stale(report: dict):
    for e in report.get("stale_known_drift", []):
        print(f"::warning file=tests/fixtures/tutorials/known_drift.json,title=stale known drift::"
              f"tutorial {e['tutorial']} no longer shows '{e['line'][:120]}': remove the entry")


def _render(report: dict, markdown: bool = False) -> list:
    out = ["## Tutorials as tests" if markdown else "== tutorials as tests =="]
    for t, r in report["tutorials"].items():
        blocks = r["blocks"]
        ran = [b for b in blocks if not b["skip"]]
        bad = [b for b in ran if b["problems"]]
        badf = [f for f in r["facts"] if not f["ok"]]
        head = (f"{r['file']}: {len(ran)} block(s) run ({sum(b['checked'] for b in ran)} output lines checked), "
                f"{len(blocks) - len(ran)} skipped, {len(r['facts'])} fact(s); {r['seconds']} s: "
                + ("OK" if not bad and not badf else f"{len(bad)} block(s) and {len(badf)} fact(s) differ")
                + (f" ({nk} known-drift line(s), not failing)" if (nk := sum(len(b.get('known') or []) for b in ran))
                   else ""))
        out.append(("### " if markdown else "") + head)
        for b in blocks:
            if b["skip"]:
                out.append(f"- skipped {b['where']}: {b['skip']}" if markdown else f"  skipped {b['where']}: {b['skip']}")
        for b in bad:
            out.append(f"  DIFFERS {b['where']}: " + "; ".join(b["problems"]))
            for m in b["missing"][:15]:
                out.append(f"      documented, not printed: {m}")
            if len(b["missing"]) > 15:
                out.append(f"      ... and {len(b['missing']) - 15} more")
            if b["missing"] or b["status"] not in (0, None):
                out.append("      printed (last lines):")
                out += [f"      | {x}" for x in b["output"][-25:]]
        for f in r["facts"]:
            if not f["ok"]:
                out.append(f"  DIFFERS fact '{f['what']}': the text says {f['stated']}, the run gives {f['actual']}")
        for b in ran:
            for m in b.get("known") or []:
                out.append(f"  known drift (not failing) {b['where']}: {m[:160]}")
    for e in report.get("stale_known_drift", []):
        out.append(f"  stale known-drift entry (the tutorial no longer shows it; remove it from "
                   f"tests/fixtures/tutorials/known_drift.json): {e['tutorial']}: {e['line'][:120]}")
    out.append(f"total {report['seconds']} s, COAL {report['coal_version']}: "
               + ("no differences" if not report["failures"] else f"{report['failures']} difference(s)")
            + (f"; {nk} known-drift line(s)" if (nk := sum(len(b.get("known") or []) for r in report["tutorials"].values()
                                                         for b in r["blocks"])) else ""))
    if markdown:
        out = [ln if not ln.startswith("  ") else "- " + ln.strip().replace("|", "\\|") for ln in out]
    return out


if __name__ == "__main__":
    sys.exit(main())
