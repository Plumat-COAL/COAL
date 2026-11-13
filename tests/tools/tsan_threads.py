#!/usr/bin/env python3
"""The threaded commands of a ThreadSanitizer build: `coal-cpp kaks / dnds / codon-dnds -j 4` and `coal-cpp-align -j 4`.

    cmake -S . -B build-tsan -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCOAL_WARNINGS_AS_ERRORS=ON \
        -DCMAKE_CXX_FLAGS="-fsanitize=thread -g -O1" -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=thread"
    cmake --build build-tsan -j 2
    python tests/tools/tsan_threads.py --build build-tsan [--seed N] [--jobs 4] [--keep DIR]

nightly.yml's `tsan` job runs this. Every command runs with -j JOBS and with -j 1; a run fails when

  * ThreadSanitizer reports anything (the text "ThreadSanitizer" on stderr, or its exit code 66), or
  * the command fails, or its output with -j JOBS is not, byte for byte, its output with -j 1.

What runs (the inputs are large enough for every thread to get work, and partly random: --seed):

  kaks, dnds    all nine methods (NG, LWL, MLWL, LPB, MLPB, YN, MYN, Comeron, Li), the default flavours and the
                reference ones (--yn-flavour paml, --comeron-flavour gestimator, --li-flavour seqinr), on the
                engine's fixtures (tests/kaks/*.axt) and on random short pairs (benchmarks/kaks/random_short_pairs.py)
  codon-dnds    the pipeline's engine step (protein alignment, codon alignment, Ka/Ks) on groups made from
                tests/fixtures/pipeline, all nine methods, with every optional output
  align         coal-cpp-align: global alignments of random protein pairs, both presets, with --scores

Prints one line per command and a summary line starting with "tsan threads:". Exit 0 when nothing failed.
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
METHODS = ("NG", "LWL", "MLWL", "LPB", "MLPB", "YN", "MYN", "Comeron", "Li")
REFERENCE = ("--yn-flavour", "paml", "--comeron-flavour", "gestimator", "--li-flavour", "seqinr")
AA = "ACDEFGHIKLMNPQRSTVWY"


def find(build: Path, name: str) -> Path:
    for cand in (build / name, build / (name + ".exe"), build / "Release" / (name + ".exe")):
        if cand.is_file():
            return cand
    sys.exit(f"tsan_threads: no {name} in {build}")


def read_fasta(path: Path) -> dict:
    out, name = {}, None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = ""
        elif name is not None:
            out[name] += line.strip()
    return out


def make_inputs(tmp: Path, seed: int) -> dict:
    """The input files, under tmp. Returns {"axt": [...], "groups": path, "proteins": path}."""
    spec = importlib.util.spec_from_file_location("rsp", ROOT / "benchmarks" / "kaks" / "random_short_pairs.py")
    rsp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rsp)
    axt = [ROOT / "tests" / "kaks" / n for n in ("comeron_li_cases.axt", "edge_cases.axt", "sim_small.axt")]
    for name, pairs in (("short_single", rsp.make_pairs(1500, seed)), ("short_multi", rsp.make_pairs_multi(1500, seed))):
        p = tmp / f"{name}.axt"
        with open(p, "w", encoding="utf-8", newline="\n") as fh:
            for k, (_, a, b) in enumerate(pairs):
                fh.write(f"{name}{k}\n{a}\n{b}\n\n")
        axt.append(p)
    # many copies of the simulated genes (300 codons each): long pairs, so that a batch is real work for 4 threads
    sim = (ROOT / "tests" / "kaks" / "sim_small.axt").read_text(encoding="utf-8").split("\n\n")
    with open(tmp / "sim_x8.axt", "w", encoding="utf-8", newline="\n") as fh:
        for rep in range(8):
            for rec in sim:
                lines = rec.strip().splitlines()
                if len(lines) == 3:
                    fh.write(f"{lines[0]}_r{rep}\n{lines[1]}\n{lines[2]}\n\n")
    axt.append(tmp / "sim_x8.axt")

    # groups for codon-dnds: the gene families (aligned pairwise or, for two members, by the engine) and the
    # ortholog pairs of the pipeline fixtures, repeated so that there are a few hundred groups
    fx = ROOT / "tests" / "fixtures" / "pipeline"
    q, s, fam = (read_fasta(fx / n) for n in ("query_cds.fasta", "subject_cds.fasta", "family_cds.fasta"))
    rng = random.Random(seed)
    lines = []
    qs, ss = list(q.items()), list(s.items())
    for rep in range(12):
        for k in range(min(len(qs), len(ss))):
            (a, ca), (b, cb) = qs[k], ss[(k + rep) % len(ss)]
            lines += [f"G\tp{rep}_{k}\tall", f"M\t{a}\t{ca}", f"M\t{b}\t{cb}", "E"]
        members = rng.sample(sorted(fam), 4)
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                lines += [f"G\tf{rep}_{i}_{j}\tall", f"M\t{members[i]}\t{fam[members[i]]}",
                          f"M\t{members[j]}\t{fam[members[j]]}", "E"]
    groups = tmp / "groups.tsv"
    groups.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")

    proteins = tmp / "proteins.fasta"
    with open(proteins, "w", encoding="utf-8", newline="\n") as fh:
        for k in range(300):
            a = "M" + "".join(rng.choice(AA) for _ in range(rng.randint(60, 400)))
            b = "".join(c if rng.random() > 0.2 else rng.choice(AA) for c in a)
            cut = rng.randint(0, 20)
            fh.write(f">a{k}\n{a}\n>b{k}\n{b[cut:]}\n")
    return {"axt": axt, "groups": groups, "proteins": proteins}


class Runner:
    def __init__(self, jobs: int, keep: Path | None):
        self.jobs, self.keep = jobs, keep
        self.n = self.bad = 0
        self.env = dict(os.environ)
        # exitcode=66: TSan's own, told apart from the program's; halt_on_error=0: every report of a run is printed
        self.env.setdefault("TSAN_OPTIONS", "exitcode=66:halt_on_error=0:second_deadlock_stack=1")

    def one(self, label: str, cmd: list, outputs: list[Path]) -> None:
        """cmd has the placeholder "@J" where the number of threads goes; outputs: the files it writes (stdout is
        compared as well). Runs it with self.jobs threads and with one, and compares."""
        self.n += 1
        got = {}
        for j in (self.jobs, 1):
            for p in outputs:
                if p.is_dir():
                    shutil.rmtree(p)
                    p.mkdir()
                elif p.exists():
                    p.unlink()
            r = subprocess.run([str(c) if c != "@J" else str(j) for c in cmd], capture_output=True, env=self.env,
                               stdin=subprocess.DEVNULL, timeout=3600)
            err = r.stderr.decode("utf-8", "replace")
            files = {}
            for p in outputs:
                for f in ([p] if p.is_file() else sorted(p.glob("*")) if p.is_dir() else []):
                    files[str(f)] = f.read_bytes()
            got[j] = (r.returncode, r.stdout, files)
            why = None
            if "ThreadSanitizer" in err or r.returncode == 66:
                why = "ThreadSanitizer report"
            elif r.returncode != 0:
                why = f"exit {r.returncode}"
            elif not r.stdout and not files:
                why = "no output"
            if why:
                self.fail(label, j, why, err, cmd)
                return
        if got[self.jobs] != got[1]:
            self.fail(label, self.jobs, f"output with -j {self.jobs} differs from -j 1", "", cmd)
            return
        print(f"ok    {label}", flush=True)

    def fail(self, label, j, why, err, cmd) -> None:
        self.bad += 1
        line = " ".join(str(c) if c != "@J" else str(j) for c in cmd)
        print(f"FAIL  {label} (-j {j}): {why}\n      {line}", flush=True)
        if err:
            print("\n".join("      " + ln for ln in err.splitlines()[:60]), flush=True)
        if self.keep is not None:
            self.keep.mkdir(parents=True, exist_ok=True)
            with open(self.keep / f"tsan{self.bad}.txt", "w", encoding="utf-8") as fh:
                fh.write(f"{label} (-j {j}): {why}\n{line}\n\n{err}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--build", default="build-tsan", help="the build folder with coal-cpp and coal-cpp-align")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--keep", help="write the report of each failed command into this folder")
    a = ap.parse_args(argv)
    build = Path(a.build)
    cpp, align = find(build, "coal-cpp"), find(build, "coal-cpp-align")
    run = Runner(a.jobs, Path(a.keep) if a.keep else None)
    methods = [x for m in METHODS for x in ("-m", m)]
    with tempfile.TemporaryDirectory(prefix="coal_tsan_") as t:
        tmp = Path(t)
        inp = make_inputs(tmp, a.seed)
        for axt in inp["axt"]:
            for name, flav in (("default", ()), ("reference flavours", REFERENCE)):
                run.one(f"kaks {axt.name}, {name}", [cpp, "kaks", "-i", axt, "-j", "@J", "-q", *methods, *flav], [])
                run.one(f"dnds {axt.name}, {name}",
                        [cpp, "dnds", axt, "--format", "axt", "-j", "@J", "-q", "--digits", "15", *methods, *flav], [])
        out = [tmp / n for n in ("pairs.tsv", "codon.fa", "protein.fa", "out.kaks", "failures.tsv", "codon_dir")]
        out[-1].mkdir()
        for name, flav in (("default", ()), ("reference flavours, terminal stops compared",
                                             (*REFERENCE, "--terminal-stops", "compare"))):
            run.one(f"codon-dnds, {name}",
                    [cpp, "codon-dnds", "-i", inp["groups"], "-o", out[0], "--codon-out", out[1], "--protein-out",
                     out[2], "--kaks-out", out[3], "--failures", out[4], "--codon-dir", out[5], "--group-column",
                     "group", "-j", "@J", "-q", *methods, *flav], out)
        scores = tmp / "scores.tsv"
        for preset in ("needle", "pwalign"):
            run.one(f"align, preset {preset}",
                    [align, inp["proteins"], "--preset", preset, "-j", "@J", "--scores", scores], [scores])
    print(f"tsan threads: seed {a.seed}, {run.n} commands with -j {a.jobs} and -j 1, {run.bad} failed", flush=True)
    return 1 if run.bad else 0


if __name__ == "__main__":
    sys.exit(main())
