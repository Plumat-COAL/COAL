"""`coal pipeline` never writes through a symbolic link in its output folder (AUDIT-99, L1).

A link put at any name the pipeline writes, in DIR or DIR/work (outputs.json.tmp, state.json.tmp, report.json.tmp,
the .tmp of every table, the per-cluster .fasta.tmp, the logs, the decompressed and --prefix-ids copies of the
inputs, the translated proteins, DIAMOND's hits and databases, the sorted proteins for DeepClust) pointed at a
file outside the folder: the next run (and the same run with --resume) wrote into that file, among them the
JSON record, the logs and the proteins. Each test runs once, replaces every file of the folder with a link to a
"victim" outside it and plants a link at every .tmp and log name, runs again and with --resume, and checks that
no victim changed and that the names COAL wrote are files again, not links. POSIX (symbolic links); needs
build/coal-cpp (COAL_CPP_BIN), and DIAMOND for the DIAMOND routes.
"""
from __future__ import annotations

import gzip
import io
import os
import shutil
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from test_dnds_pipeline import (CLUSTALO, FAM, PLAIN, PipelineCase, Q, S, needs_cpp, needs_diamond,  # noqa: E402
                                read_fasta)

from coal import dnds_pipeline as P  # noqa: E402

# names a run writes that the first run may have removed again or written only through a temporary name
EXTRA = ("work/outputs.json.tmp", "work/state.json.tmp", "work/inputs/sources.json.tmp", "work/logs/translate.log",
         "work/logs/diamond.log", "work/logs/codon-dnds.log", "work/subject.dmnd", "work/query.dmnd",
         "realign.tsv.body")
VICTIM = b"keep me\n"


def snapshot(folder: Path) -> dict:
    """Every name under `folder` (links not followed) and the bytes of its files."""
    got = {}
    for dp, dn, fn in os.walk(folder):
        for n in dn + fn:
            p = Path(dp, n)
            got[str(p.relative_to(folder))] = (os.readlink(p) if p.is_symlink() else
                                               p.read_bytes() if p.is_file() else None)
    return got


def run(fn, *args, **kw):
    with redirect_stderr(io.StringIO()) as err:
        res = fn(*args, **kw)
    return res, err.getvalue()


class _Planting(PipelineCase):
    def plant(self, out: Path) -> dict[str, Path]:
        """Every file of `out` (but the record itself, which is read: a link there makes the folder another's) and
        every .tmp and log name becomes a link to a file of its own outside `out`."""
        names = {str(p.relative_to(out)) for p in out.rglob("*") if p.is_file() and not p.is_symlink()}
        names |= {n + ".tmp" for n in names if not n.endswith(".tmp")} | set(EXTRA)
        names.discard("work/outputs.json")
        victims = self.tmp / f"victims-{out.name}"
        victims.mkdir(exist_ok=True)
        planted = {}
        for n in sorted(names):
            p = out / n
            if not p.parent.is_dir():
                continue
            v = victims / n.replace("/", "__")
            v.write_bytes(VICTIM)
            if os.path.lexists(p):
                p.unlink()
            p.symlink_to(v)
            planted[n] = v
        return planted

    def check(self, out: Path, planted: dict[str, Path], what: str) -> None:
        written = sorted(n for n, v in planted.items() if v.read_bytes() != VICTIM)
        self.assertEqual(written, [], f"{what}: written through a link")
        links = sorted(str(p.relative_to(out)) for p in out.rglob("*") if p.is_symlink() and p.exists())
        # what is still a link is a name this run did not write (and so did not read either)
        self.assertFalse([n for n in links if n in ("work/outputs.json", "work/state.json", "report.json",
                                                    "dnds.tsv", "failed.tsv", "genes.tsv")], (what, links))

    def twice(self, fn, *args, ok=(0, 2), **kw) -> Path:
        """fn(*args, **kw) (the output folder is its last argument), then twice more over planted links."""
        out = Path(args[-1])
        res, err = run(fn, *args, **kw)
        self.assertIn(res.exit_code, ok, err)
        for resume in (False, True):
            planted = self.plant(out)
            self.assertGreater(len(planted), 20)
            res, err = run(fn, *args, **kw, resume=resume)
            self.assertIn(res.exit_code, ok, (resume, err))
            self.check(out, planted, f"resume={resume}")
        return out


@unittest.skipUnless(os.name == "posix", "symbolic links (POSIX)")
@needs_cpp
class Planted(_Planting):
    def test_pairs_from_given_pairs_with_a_compressed_input(self):
        qgz = self.tmp / "q.fa.gz"
        qgz.write_bytes(gzip.compress(Q.read_bytes()))
        pf = self.tmp / "pairs.tsv"
        pf.write_text("".join(f"{a}\t{b}\n" for a, b in PLAIN[:4]), encoding="utf-8")
        out = self.twice(P.run_pairs, qgz, S, self.tmp / "o", pairs=pf, kaks_format=True, keep_alignments=True,
                         keep_protein_alignments=True)
        self.assertFalse((out / "work" / "inputs").is_symlink())
        self.assertTrue(any(p.is_file() and not p.is_symlink() for p in (out / "work" / "inputs").iterdir()))

    def test_clusters_with_both_alignment_folders(self):
        ids = list(read_fasta(FAM))      # clusters of two (aligned by NW: no Clustal Omega needed)
        tsv = self.tmp / "clusters.tsv"
        tsv.write_text("".join(f"{ids[k]}\t{ids[k]}\n{ids[k]}\t{ids[k + 1]}\n" for k in range(0, 8, 2)),
                       encoding="utf-8")
        out = self.twice(P.run_clusters, FAM, self.tmp / "c", clusters=tsv, keep_alignments=True,
                         keep_protein_alignments=True, kaks_format=True)
        for d in ("codon_alignments", "protein_alignments"):
            files = list((out / d).iterdir())
            self.assertTrue(files and all(p.is_file() and not p.is_symlink() for p in files), d)


@unittest.skipUnless(os.name == "posix", "symbolic links (POSIX)")
@needs_cpp
class FolderLinks(_Planting):
    """A folder COAL writes into that is a symbolic link (AUDIT-102, M1): os.makedirs(exist_ok=True) took it for a
    folder and O_NOFOLLOW guards only the last part of a path, so work/logs or work/inputs linked elsewhere had the
    logs, the decompressed copies and sources.json written there. DIR/work, work/logs, work/inputs and an
    alignment folder that is a link are refused (exit 1) before anything changes, with --overwrite too."""

    def inputs(self):
        qgz = self.tmp / "q.fa.gz"
        qgz.write_bytes(gzip.compress(Q.read_bytes()))
        pf = self.tmp / "pairs.tsv"
        pf.write_text("".join(f"{a}\t{b}\n" for a, b in PLAIN[:4]), encoding="utf-8")
        return qgz, pf

    def elsewhere(self, name: str) -> Path:
        d = self.tmp / name
        d.mkdir()
        for f in ("translate.log", "codon-dnds.log", "0_q.fa", "sources.json", "outputs.json", "state.json"):
            (d / f).write_bytes(VICTIM)
        return d

    def test_folder_links_in_work_are_refused(self):
        qgz, pf = self.inputs()
        out = self.tmp / "o"
        res, err = run(P.run_pairs, qgz, S, out, pairs=pf)
        self.assertIn(res.exit_code, (0, 2), err)
        for sub in ("work/logs", "work/inputs", "work"):
            target = self.elsewhere(sub.replace("/", "_"))
            real = out / (sub + ".real")
            (out / sub).rename(real)
            (out / sub).symlink_to(target, target_is_directory=True)
            for kw in ({}, {"overwrite": True}, {"resume": True}):
                before = snapshot(out)
                res, err = run(P.run_pairs, qgz, S, out, pairs=pf, **kw)
                self.assertEqual(res.exit_code, 1, (sub, kw, err))
                self.assertEqual(snapshot(out), before, (sub, kw))        # refused before anything changed
                self.assertIn(f"{out / sub} is a symbolic link (to {target})", err)
                self.assertIn("give -o a folder on that file system", err)
                self.assertEqual(sorted(p.name for p in target.iterdir()),
                                 ["0_q.fa", "codon-dnds.log", "outputs.json", "sources.json", "state.json",
                                  "translate.log"], (sub, kw))
                self.assertTrue(all(p.read_bytes() == VICTIM for p in target.iterdir()), (sub, kw))
            (out / sub).unlink()
            real.rename(out / sub)
            res, err = run(P.run_pairs, qgz, S, out, pairs=pf)
            self.assertIn(res.exit_code, (0, 2), (sub, err))

    def test_a_folder_link_made_while_the_run_starts_is_refused(self):
        """make_folder() checks again after os.makedirs: a link put at work/logs after the checks."""
        qgz, pf = self.inputs()
        out = self.tmp / "o"
        target = self.elsewhere("late")
        real = os.makedirs

        def late(path, *a, **k):     # (the link appears while the folder is made)
            if os.path.basename(os.fspath(path)) == "logs" and not os.path.lexists(path):
                os.symlink(target, path)
            return real(path, *a, **k)
        with mock.patch.object(os, "makedirs", late):
            res, err = run(P.run_pairs, qgz, S, out, pairs=pf)
        self.assertEqual(res.exit_code, 1, err)
        self.assertIn("logs is a symbolic link", err)
        self.assertTrue(all(p.read_bytes() == VICTIM for p in target.iterdir()))


@unittest.skipUnless(os.name == "posix", "symbolic links (POSIX)")
@needs_cpp
class LinksMadeAfterTheRemoval(_Planting):
    """AUDIT-102, L1: a name COAL removed before a program wrote it could get a link again before the program
    opened it (DIAMOND opens -o after its search). The programs now write into a new private folder in work/ whose
    files are renamed into place; gunzip_to creates its copy with O_EXCL | O_NOFOLLOW. Simulated: every
    remove_file() of the run puts a link to a victim at the name right after removing it."""

    def run_with_late_links(self, fn, *args, **kw):
        victims = self.tmp / "late-victims"
        victims.mkdir(exist_ok=True)
        planted = {}
        real = P.remove_file

        def remove_then_link(path):
            real(path)
            if os.path.isdir(os.path.dirname(path)):
                v = victims / f"{len(planted)}_{os.path.basename(path)}"
                v.write_bytes(VICTIM)
                os.symlink(v, path)
                planted[path] = v
        with mock.patch.object(P, "remove_file", remove_then_link):
            res, err = run(fn, *args, **kw)
        self.assertTrue(planted)
        self.assertEqual([p for p, v in planted.items() if v.read_bytes() != VICTIM], [], err)
        return res, err, planted

    def test_programs_never_write_through_a_late_link(self):
        pf = self.tmp / "pairs.tsv"
        pf.write_text("".join(f"{a}\t{b}\n" for a, b in PLAIN[:4]), encoding="utf-8")
        out = self.tmp / "o"
        res, err = run(P.run_pairs, Q, S, out, pairs=pf)
        self.assertIn(res.exit_code, (0, 2), err)
        res, err, planted = self.run_with_late_links(P.run_pairs, Q, S, out, pairs=pf)
        self.assertIn(res.exit_code, (0, 2), err)
        for name in ("query.faa", "query.genes.tsv", "subject.faa", "subject.genes.tsv"):
            p = out / "work" / name
            self.assertTrue(str(p) in planted and p.is_file() and not p.is_symlink(), name)   # replaced, not followed
        self.assertFalse([n for n in os.listdir(out / "work") if ".staging" in n])          # staging folders gone

    def test_a_late_link_at_the_decompressed_copy_stops_the_run(self):
        qgz = self.tmp / "q.fa.gz"
        qgz.write_bytes(gzip.compress(Q.read_bytes()))
        pf = self.tmp / "pairs.tsv"
        pf.write_text("".join(f"{a}\t{b}\n" for a, b in PLAIN[:4]), encoding="utf-8")
        res, err, planted = self.run_with_late_links(P.run_pairs, qgz, S, self.tmp / "o", pairs=pf)
        self.assertEqual(res.exit_code, 1, err)          # gunzip_to's O_EXCL: "File exists", never followed
        self.assertIn("File exists", err)


@unittest.skipUnless(os.name == "posix", "symbolic links (POSIX)")
@needs_cpp
class ResumeNeverReusesALink(_Planting):
    """--resume reuses a stage's output, or the decompressed copy of an input, only when it is a file of the
    folder's own, not a link (AUDIT-102, S1: the other tests also link state.json, so nothing was reused there)."""

    def test_a_linked_stage_output_or_copy_is_made_again(self):
        qgz = self.tmp / "q.fa.gz"
        qgz.write_bytes(gzip.compress(Q.read_bytes()))
        pf = self.tmp / "pairs.tsv"
        pf.write_text("".join(f"{a}\t{b}\n" for a, b in PLAIN[:4]), encoding="utf-8")
        out = self.tmp / "o"
        res, err = run(P.run_pairs, qgz, S, out, pairs=pf)
        self.assertIn(res.exit_code, (0, 2), err)
        copy = next(p for p in (out / "work" / "inputs").iterdir() if p.name.startswith("0_"))
        links = {}
        for p in (out / "work" / "subject.faa", copy):   # a link to a file with the same bytes (and times)
            v = self.tmp / ("victim-" + p.name)
            shutil.copy2(p, v)
            p.unlink()
            p.symlink_to(v)
            links[p] = v
        stat_before = {p: v.stat().st_mtime_ns for p, v in links.items()}
        res, err = run(P.run_pairs, qgz, S, out, pairs=pf, resume=True)
        self.assertIn(res.exit_code, (0, 2), err)
        stages = {s["name"]: s for s in res.report["stages"]}
        self.assertFalse(stages["translate_subject"]["cached"])     # its output was a link: made again
        for p, v in links.items():
            self.assertTrue(p.is_file() and not p.is_symlink(), p.name)
            self.assertEqual(v.stat().st_mtime_ns, stat_before[p], p.name)


def masked_report(path: Path) -> dict:
    """report.json without its times and durations, input modification times and resume signatures."""
    import json

    def mask(o, key=None):
        if isinstance(o, dict):
            return {k: ("<t>" if k in ("started", "finished", "seconds", "modified", "sig") or k.endswith("_seconds")
                        else mask(v, k)) for k, v in o.items()}
        if isinstance(o, list):
            return [mask(v) for v in o]
        return o
    return mask(json.loads(path.read_text(encoding="utf-8")))


@unittest.skipUnless(os.name == "posix", "symbolic links and a POSIX shell stand-in for DIAMOND")
@needs_cpp
class StagingFolders(PipelineCase):
    """The folder an external program writes into is work/NAME.staging, made anew (AUDIT-102 re-audit, M2 and M3):
    a killed run's folder (it can hold DIAMOND's databases) is removed by the next run, whatever its kind, and the
    commands report.json records name the same paths in every run."""

    def stand_in(self) -> Path:
        """A DIAMOND that writes one cluster per sequence, and stalls inside `deepclust` (after writing part of its
        output and a 'database') while the file STALL exists."""
        fake = self.tmp / "diamond"
        stall = self.tmp / "STALL"
        fake.write_text(
            "#!/bin/sh\n"
            'if [ "$1" = version ]; then echo "diamond version 2.2.8"; exit 0; fi\n'
            'db=; out=\n'
            'while [ $# -gt 0 ]; do case "$1" in --db) db=$2; shift;; --out) out=$2; shift;; esac; shift; done\n'
            'echo partial > "$out"; echo big > "$(dirname "$out")/tmp.dmnd"\n'
            f'while [ -e "{stall}" ]; do sleep 0.1; done\n'
            "grep '^>' \"$db\" | sed 's/^>//; s/[[:space:]].*//' | awk '{print $1 \"\\t\" $1}' > \"$out\"\n",
            encoding="utf-8")
        fake.chmod(0o755)
        return fake

    def cli(self, out: Path, *extra, wait=True):
        import subprocess
        from subprocess_env import python_env
        args = [sys.executable, "-m", "coal", "pipeline", "clusters", "--cds", str(FAM), "-o", str(out), "-q",
                "-j", "1", "--diamond", str(self.tmp / "diamond"), *extra]
        if wait:
            return subprocess.run(args, env=python_env(), stdin=subprocess.DEVNULL, capture_output=True, text=True,
                                  timeout=300, check=False)
        return subprocess.Popen(args, env=python_env(), stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                                stderr=subprocess.DEVNULL, start_new_session=True)

    def staging(self, out: Path) -> list[str]:
        return sorted(n for n in os.listdir(out / "work") if ".staging" in n) if (out / "work").is_dir() else []

    def test_a_killed_runs_staging_folder_is_removed_by_the_next_run(self):
        import signal
        import time
        self.stand_in()
        out = self.tmp / "o"
        for again in (["--resume"], []):
            (self.tmp / "STALL").write_text("x", encoding="utf-8")
            proc = self.cli(out, wait=False)
            t0 = time.monotonic()
            while not (out / "work" / "deepclust.staging" / "tmp.dmnd").exists() and time.monotonic() - t0 < 120:
                time.sleep(0.05)
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait(timeout=60)
            (self.tmp / "STALL").unlink()
            self.assertEqual(self.staging(out), ["deepclust.staging"])          # what a kill -9 leaves
            # an old-style folder of a development build, and a link with a staging name, too
            (out / "work" / "translate_proteins.staging-ab12cd34").mkdir()
            elsewhere = self.tmp / f"elsewhere{len(again)}"
            elsewhere.mkdir()
            (elsewhere / "keep.txt").write_text("keep me\n", encoding="utf-8")
            (out / "work" / "realign.staging").symlink_to(elsewhere, target_is_directory=True)
            Path(str(out) + ".coal-lock").unlink(missing_ok=True)
            r = self.cli(out, *again)
            self.assertIn(r.returncode, (0, 2), r.stderr)
            self.assertEqual(self.staging(out), [], again)
            self.assertEqual((elsewhere / "keep.txt").read_text(encoding="utf-8"), "keep me\n")
            self.assertTrue((out / "clusters.tsv").is_file())

    def test_the_recorded_commands_are_the_same_in_every_run(self):
        self.stand_in()
        reports = []
        for k in range(2):
            r = self.cli(self.tmp / "o")
            self.assertIn(r.returncode, (0, 2), r.stderr)
            reports.append(masked_report(self.tmp / "o" / "report.json"))
        self.assertEqual(reports[0], reports[1])
        cmds = [c for st in reports[0]["stages"] for c in st.get("commands", ())]
        dc = next(c for c in cmds if c[1:2] == ["deepclust"])
        self.assertEqual(dc[dc.index("--out") + 1], str(self.tmp / "o" / "work" / "deepclust.staging" / "clusters.tsv"))
        tr = next(c for c in cmds if c[1:2] == ["translate"])
        self.assertEqual(tr[tr.index("-o") + 1],
                         str(self.tmp / "o" / "work" / "translate_proteins.staging" / "proteins.faa"))


@unittest.skipUnless(os.name == "posix", "symbolic links (POSIX)")
@needs_diamond
class PlantedDiamond(_Planting):
    """The routes that run DIAMOND: makedb's databases, blastp's hits, DeepClust's sorted input, realign."""

    def test_pairs_found_by_diamond(self):
        self.twice(P.run_pairs, Q, S, self.tmp / "o", kaks_format=True)

    @unittest.skipUnless(CLUSTALO, "needs Clustal Omega (DeepClust's clusters have more than two members)")
    def test_clusters_by_deepclust_with_prefixed_ids(self):
        self.twice(P.run_clusters, [FAM, S], self.tmp / "c", prefix_ids=True, realign=True, keep_alignments=True)


if __name__ == "__main__":
    unittest.main()
