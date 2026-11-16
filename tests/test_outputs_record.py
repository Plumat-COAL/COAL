"""coal.outputs_record: whose files are in a `coal pipeline` output folder (work/outputs.json and its rules), tested
on plain folders, without coal-cpp or any tool. Each case is a behaviour the reviews and audits of #87, #95, #96
and #98 established:

  * a .fasta of the user's in codon_alignments/ or protein_alignments/ is kept by the next run's clean-up, and a
    run that writes that folder again is refused;
  * a run killed while its engine wrote (no list in the record) leaves files that are all its own: the next run
    removes them, and the same run with --resume is not refused;
  * a work/ that no COAL run made is refused unless --overwrite (an empty one, one with only an empty logs/, one
    with a cut-off record that begins as COAL writes it are COAL's); a `work` that is not a folder always;
  * --report onto a file that is not a COAL report is refused unless --overwrite; one that is a folder, a named
    pipe or a device always;
  * a tampered record or report cannot make the clean-up remove anything outside the folder or the files COAL
    names (links are removed, never followed);
  * a record of 0.6.0 (no per-folder lists) clears the engine's lower-case .fasta files as before;
  * .fasta in any other case (X.FASTA) is never COAL's;
  * a folder that is not there when the engine ends (it could not be made) gets an empty list: with none, a folder
    of that name that the user made later had every .fasta in it removed by the next run.

The last case also runs through the pipeline with os.makedirs made to fail (class AlignmentFolderNotMade; needs
build/coal-cpp, COAL_CPP_BIN).
"""
from __future__ import annotations

import errno
import io
import json
import os
import re
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import __version__  # noqa: E402
from coal import dnds_pipeline as P  # noqa: E402
from coal import outputs_record as R  # noqa: E402
from coal.outputs_record import OutputsRecord, OutputsRefused  # noqa: E402
from test_dnds_pipeline import FAM, PipelineCase, needs_cpp, read_fasta  # noqa: E402

SCHEMA = "coal.pipeline/1"
CODON, PROT = "codon_alignments", "protein_alignments"
# what a cluster run with --keep-alignments plans to write (dnds_pipeline._Run._planned_outputs)
PLANNED = ["failed.tsv", "genes.tsv", "dnds.tsv", "report.json", "cluster_summary.tsv"]
PLANNED += [n + ".tmp" for n in PLANNED]


def snapshot(folder: Path) -> dict:
    return {str(p.relative_to(folder)): (p.read_bytes() if p.is_file() and not p.is_symlink() else None)
            for p in sorted(folder.rglob("*"))}


class Folder(unittest.TestCase):
    """An output folder in a temporary directory, and the steps of a run on it as _Run makes them."""

    def symlink(self, target, link) -> None:
        try:
            os.symlink(target, link)
        except (OSError, NotImplementedError) as e:     # (Windows without the privilege)
            self.skipTest(f"cannot make a symbolic link: {e}")

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="coal-owned-"))
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.out = self.tmp / "out"
        self.out.mkdir()

    def put(self, rel: str, text: str = "x\n") -> Path:
        p = self.out / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
        return p

    def record_on_disk(self) -> dict:
        return json.loads((self.out / "work" / "outputs.json").read_text(encoding="utf-8"))

    def write_record(self, rec: dict) -> None:
        (self.out / "work").mkdir(exist_ok=True)
        (self.out / "work" / "outputs.json").write_text(json.dumps(rec), encoding="utf-8")

    def earlier(self, writes=(), report=None, report_is_default=False) -> OutputsRecord:
        owned = OutputsRecord(str(self.out), SCHEMA)
        owned.read_earlier(report, report_is_default, writes)
        return owned

    def start(self, writes=(), planned=PLANNED, keep=(), is_input=lambda p: None, **check) -> OutputsRecord:
        """A run's start: what the earlier run wrote, the checks, its clean-up, the new record (written)."""
        owned = self.earlier(writes, check.get("report"), check.get("report_is_default", False))
        owned.check_foreign(planned, **check)
        owned.clear(keep, is_input)
        (self.out / "work" / "logs").mkdir(parents=True, exist_ok=True)
        owned.new_record([n for n in planned if not (n == "orthologs.tsv" and check.get("pairs_in_place"))], writes)
        owned.write()
        return owned

    def engine(self, owned: OutputsRecord, files: dict, killed: bool = False) -> None:
        """The engine stage: the list dropped, the folders made and written, the list made when it ends."""
        if owned.engine_starts():
            owned.write()
        for folder, names in files.items():
            (self.out / folder).mkdir(exist_ok=True)
            for n in names:
                (self.out / folder / n).write_text(">a\nATG\n", encoding="utf-8")
        if not killed and owned.list_alignment_files():
            owned.write()

    def cluster_run(self, writes=(CODON,), files=None, killed=False, **kw) -> OutputsRecord:
        owned = self.start(writes, **kw)
        files = files if files is not None else {w: ["c1.fasta", "c2.fasta"] for w in writes}
        self.engine(owned, files, killed)
        for n in PLANNED:
            if not n.endswith(".tmp"):
                self.put(n)
        return owned

    def names(self, folder: str) -> list[str]:
        d = self.out / folder
        return sorted(os.listdir(d)) if d.is_dir() else []


class TheRecord(Folder):
    def test_what_a_run_records_and_when(self):
        owned = self.start(writes=(CODON,))
        rec = self.record_on_disk()
        self.assertEqual(rec, {"coal_version": __version__, "files": sorted(PLANNED), "read_in_place": [],
                               "codon_alignments_dir": True, "protein_alignments_dir": False,
                               "codon_alignments_files": []})
        self.assertTrue(owned.engine_starts())
        owned.write()
        self.assertNotIn("codon_alignments_files", self.record_on_disk())
        self.assertFalse(owned.engine_starts())          # nothing more to drop
        (self.out / CODON).mkdir()
        for n in ("b.fasta", "a.fasta", "c.fasta.tmp", "notes.txt"):
            (self.out / CODON / n).write_text("x\n", encoding="utf-8")
        self.assertTrue(owned.list_alignment_files())
        owned.write()
        self.assertEqual(self.record_on_disk()["codon_alignments_files"], ["a.fasta", "b.fasta"])
        self.assertEqual(sorted(os.listdir(self.out / "work")), ["logs", "outputs.json"])   # no .tmp left

    def test_a_pairs_run_has_no_folder_lists(self):
        owned = self.start(planned=["dnds.tsv", "orthologs.tsv"])
        self.assertNotIn("codon_alignments_files", owned.record)
        self.assertFalse(owned.engine_starts())
        self.assertFalse(owned.list_alignment_files())

    def test_a_folder_that_is_not_there_when_the_engine_ends_gets_an_empty_list(self):
        """os.makedirs of protein_alignments/ failed (EACCES, ENOSPC) after engine_starts() dropped both lists: the
        record kept protein_alignments_dir: true with no list, and a protein_alignments/ the user made later had
        every .fasta in it removed by the next run."""
        owned = self.start(writes=(CODON, PROT))
        self.engine(owned, {CODON: ["c1.fasta"]})          # protein_alignments/ was never made
        rec = self.record_on_disk()
        self.assertEqual((rec["codon_alignments_files"], rec["protein_alignments_files"]), (["c1.fasta"], []))
        # the same for a folder that is a link (never COAL's) or a file
        self.put(PROT, "a file\n")
        self.assertTrue(owned.list_alignment_files())
        self.assertEqual(owned.record["protein_alignments_files"], [])
        (self.out / PROT).unlink()
        # the user's protein_alignments/ made afterwards: kept by a run that does not write it, refused by one that does
        self.put(f"{PROT}/mine.fasta")
        with self.assertRaisesRegex(OutputsRefused, "holds mine.fasta, which no COAL run wrote there"):
            self.earlier(writes=(PROT,))
        self.start()
        self.assertEqual(self.names(PROT), ["mine.fasta"])
        self.assertEqual(self.names(CODON), [])

    def test_the_record_cannot_be_read_back_as_a_list_when_the_folder_cannot_be_read(self):
        owned = self.start(writes=(CODON,))
        self.engine(owned, {CODON: ["c1.fasta"]}, killed=True)
        with mock.patch.object(R.os, "listdir", side_effect=PermissionError(errno.EACCES, "Permission denied")), \
                self.assertRaises(OSError):
            owned.list_alignment_files()
        # (the caller does not write it: on disk the record has no list, as while the engine wrote)
        self.assertNotIn("codon_alignments_files", self.record_on_disk())

    def test_read_in_place_is_carried_forward(self):
        """#87: an output read in place (--pairs DIR/orthologs.tsv) is kept in read_in_place while it is there and
        no run writes it."""
        self.write_record({"coal_version": "0.6.0", "files": ["dnds.tsv"], "read_in_place": ["orthologs.tsv"]})
        self.put("orthologs.tsv", "mine\n")
        self.put("dnds.tsv")
        owned = self.start(planned=["dnds.tsv", "genes.tsv"])
        self.assertEqual(owned.record["read_in_place"], ["orthologs.tsv"])
        self.assertEqual((self.out / "orthologs.tsv").read_text(encoding="utf-8"), "mine\n")   # never removed
        owned = self.start(planned=["dnds.tsv", "orthologs.tsv"])     # a run that writes it: in files, not read
        self.assertEqual(owned.record["read_in_place"], [])
        self.assertIn("orthologs.tsv", owned.record["files"])
        self.write_record({"coal_version": "0.6.0", "files": [], "read_in_place": ["orthologs.tsv"]})
        (self.out / "orthologs.tsv").unlink()
        self.assertEqual(self.start(planned=["dnds.tsv"]).record["read_in_place"], [])   # gone: dropped


class UsersFilesInTheAlignmentFolders(Folder):
    def test_a_users_fasta_is_kept(self):
        self.cluster_run(writes=(CODON, PROT))
        self.put(f"{CODON}/mine.fasta", "mine\n")
        self.put(f"{PROT}/mine.fasta", "mine\n")
        self.start()                           # a run that writes neither folder
        self.assertEqual(self.names(CODON), ["mine.fasta"])
        self.assertEqual(self.names(PROT), ["mine.fasta"])
        self.assertEqual((self.out / CODON / "mine.fasta").read_text(encoding="utf-8"), "mine\n")
        self.start()                           # and the one after it (no list now: the folder is not COAL's)
        self.assertEqual(self.names(CODON), ["mine.fasta"])

    def test_a_run_that_writes_the_folder_again_is_refused(self):
        self.cluster_run()
        self.put(f"{CODON}/mine.fasta")
        before = snapshot(self.out)
        with self.assertRaises(OutputsRefused) as cm:
            self.earlier(writes=(CODON,))
        msg = str(cm.exception)
        self.assertIn(f"{self.out / CODON} holds mine.fasta, which no COAL run wrote there: --keep-alignments", msg)
        self.assertIn("Nothing was changed", msg)
        self.assertIn("--overwrite does not apply to a folder", msg)
        self.assertEqual(snapshot(self.out), before)
        for k in range(5):
            self.put(f"{CODON}/m{k}.fasta")
        with self.assertRaisesRegex(OutputsRefused, r"holds m0.fasta, m1.fasta, m2.fasta and 3 more"):
            self.earlier(writes=(CODON,))

    def test_a_folder_no_coal_run_wrote_is_refused_when_not_empty(self):
        self.put(f"{PROT}/data.txt")
        with self.assertRaisesRegex(OutputsRefused, "exists and was not written by a COAL run: "
                                                    "--keep-protein-alignments writes one file per cluster"):
            self.earlier(writes=(PROT,))
        (self.out / PROT / "data.txt").unlink()
        self.earlier(writes=(PROT,))           # an empty folder: nothing to mix with
        self.earlier()                         # and a run that does not write it never looks

    def test_normal_reruns_go_on(self):
        for writes in ((CODON,), (CODON,), (CODON, PROT), (), (PROT,)):
            self.cluster_run(writes=writes)
            self.assertEqual(self.names(CODON), ["c1.fasta", "c2.fasta"] if CODON in writes else [])
            self.assertEqual(self.names(PROT), ["c1.fasta", "c2.fasta"] if PROT in writes else [])


class KilledRuns(Folder):
    def test_a_run_killed_while_the_engine_writes_is_resumed_and_cleared(self):
        """#98, M1: kill -9 while the engine writes; the same run again goes on, a run without the folder clears
        it (also the .fasta.tmp it was writing)."""
        self.cluster_run(files={CODON: ["c1.fasta", "c2.fasta", "c3.fasta.tmp"]}, killed=True)
        rec = self.record_on_disk()
        self.assertTrue(rec["codon_alignments_dir"])
        self.assertNotIn("codon_alignments_files", rec)
        owned = self.earlier(writes=(CODON,))              # not refused
        self.assertEqual(owned.foreign_alignment_files(CODON), [])
        self.assertEqual(owned.dirs, {CODON})
        self.cluster_run(files={CODON: ["c1.fasta"]}, killed=True)   # killed again
        self.start()
        self.assertFalse((self.out / CODON).exists())

    def test_a_run_killed_before_its_engine_leaves_an_empty_list(self):
        """#98 (audit of #96, L1): a run that ends before its engine started makes no .fasta there its own."""
        self.cluster_run()
        self.start(writes=(CODON,))                        # ... and ends here
        self.put(f"{CODON}/later.fasta")
        self.assertEqual(self.earlier().foreign_alignment_files(CODON), ["later.fasta"])
        self.start()
        self.assertEqual(self.names(CODON), ["later.fasta"])

    def test_leftovers_of_a_killed_run_are_removed(self):
        self.write_record({"coal_version": "0.6.0", "files": sorted(PLANNED), "read_in_place": []})
        for n in ("dnds.tsv.tmp", "genes.tsv", "report.json.tmp", "mine.tsv", "dnds.kaks"):
            self.put(n)
        self.start()
        self.assertEqual(sorted(p.name for p in self.out.iterdir()), ["dnds.kaks", "mine.tsv", "work"])


class WorkFolder(Folder):
    def check(self, **kw):
        self.earlier().check_foreign(PLANNED, **kw)

    def test_a_users_work_is_refused_unless_overwrite(self):
        self.put("work/notes.txt", "mine\n")
        with self.assertRaises(OutputsRefused) as cm:
            self.check()
        self.assertIn("work exists and no COAL run made it (it has no work/outputs.json)", str(cm.exception))
        self.assertIn("add --overwrite to write into it (files there with COAL's own names are replaced; nothing "
                      "else is touched)", str(cm.exception))
        self.check(overwrite=True)
        self.start(overwrite=True)
        self.assertEqual((self.out / "work" / "notes.txt").read_text(encoding="utf-8"), "mine\n")
        self.check()                                       # now COAL's: work/outputs.json

    def test_work_folders_that_are_coals(self):
        (self.out / "work").mkdir()
        self.check()                                       # empty
        (self.out / "work" / "logs").mkdir()
        self.check()                                       # only an empty logs/ (#98: killed as it started)
        self.put("work/outputs.json", '{"coal_version": "0.6')
        self.check()                                       # cut short, but begins as COAL writes it
        (self.out / "work" / "outputs.json").unlink()
        self.put("work/report.json")                       # (a folder's own COAL report.json: below)
        with self.assertRaises(OutputsRefused):
            self.check()
        (self.out / "work" / "report.json").unlink()
        self.put("report.json", json.dumps({"schema": SCHEMA, "outputs": {"dnds": "dnds.tsv"}}))
        self.put("work/state.json")
        self.check()                                       # a version without work/outputs.json

    def test_work_folders_that_are_not(self):
        for k, (rel, text) in enumerate((("work/outputs.json", '{"my": "settings"}\n'),
                                         ("work/outputs.json", "[1, 2]\n"),
                                         ("work/outputs.json", '{"coal_version" : 1'),
                                         ("work/logs/mine.log", "x\n"))):
            with self.subTest(rel=rel, text=text):
                shutil.rmtree(self.out / "work", ignore_errors=True)
                self.put(rel, text)
                with self.assertRaisesRegex(OutputsRefused, "work exists and no COAL run made it"):
                    self.check()
        self.assertTrue(R.starts_like_coal_record.__doc__.startswith("True when `path` is a regular file"))
        self.assertFalse(R.starts_like_coal_record(str(self.out / "work")))     # a folder

    def test_a_work_that_is_not_a_folder_is_always_refused(self):
        self.put("work", "a file\n")
        for kw in ({}, {"overwrite": True}):
            with self.assertRaisesRegex(OutputsRefused, "work exists and is not a folder"):
                self.check(**kw)
        (self.out / "work").unlink()
        self.symlink(self.tmp / "nowhere", self.out / "work")                    # a broken link
        with self.assertRaisesRegex(OutputsRefused, "work exists and is not a folder"):
            self.check(overwrite=True)


class ReportFile(Folder):
    def test_a_report_that_is_not_coals_is_refused_unless_overwrite(self):
        for where in (self.tmp / "notes.json", self.out / "mine.json"):
            where.write_text('{"my": "data"}\n', encoding="utf-8")
            with self.assertRaisesRegex(OutputsRefused,
                                        re.escape(f"--report {where} already exists and is not a COAL report")):
                self.earlier(report=str(where)).check_foreign(PLANNED, report=str(where))
            self.earlier().check_foreign(PLANNED, report=str(where), overwrite=True)
        rp = self.tmp / "coal.json"
        for rep in ({"schema": SCHEMA}, {"schema": "coal.report/1", "exit_code": 0},
                    {"schema": SCHEMA, "report_truncated": True}):
            rp.write_text(json.dumps(rep), encoding="utf-8")
            self.earlier().check_foreign(PLANNED, report=str(rp))
        for text in ("", "[]", '{"schema": "other/1"}', '{"schema": "coal.pipeline/1 x"}', "not json"):
            rp.write_text(text, encoding="utf-8")
            self.assertFalse(R.is_coal_report(str(rp)), text)

    def test_the_default_report_spelled_out_is_the_default(self):
        self.put("report.json", '{"my": "data"}\n')
        with self.assertRaisesRegex(OutputsRefused, "already holds report.json"):
            self.earlier(report=str(self.out / "report.json"), report_is_default=True).check_foreign(
                PLANNED, report=str(self.out / "report.json"), report_is_default=True)

    def test_a_report_that_is_not_a_regular_file_is_always_refused(self):
        d = self.tmp / "reports"
        d.mkdir()
        with self.assertRaisesRegex(OutputsRefused, re.escape(f"--report {d} is a folder: ") + ".*"
                                    + re.escape(str(d / "report.json"))):
            R.check_report_path(str(d))
        targets = [Path(os.devnull)] if os.name == "posix" else []   # (Windows: NUL is not stat-ed as a device)
        if hasattr(os, "mkfifo"):
            fifo = self.tmp / "pipe.json"
            os.mkfifo(fifo)
            link = self.tmp / "link.json"
            self.symlink(fifo, link)
            targets += [fifo, link]
            # a DIR/report.json that is a named pipe is not read (it would block) and so lists nothing
            os.mkfifo(self.out / "report.json")
            self.assertEqual(self.earlier().known, set())
        for t in targets:
            with self.assertRaisesRegex(OutputsRefused, "is not a regular file \\(a device or a named pipe\\)"):
                R.check_report_path(str(t))
        for fine in (None, "", str(self.tmp / "new.json"), str(self.tmp / "missing" / "r.json")):
            R.check_report_path(fine)
        (self.tmp / "r.json").write_text("{}", encoding="utf-8")
        R.check_report_path(str(self.tmp / "r.json"))

    def test_a_report_elsewhere_counts_when_it_is_of_this_folder(self):
        rp = self.tmp / "latest.json"
        rp.write_text(json.dumps({"schema": SCHEMA, "out_dir": str(self.out), "mode": "pairs",
                                  "outputs": {"dnds": "dnds.tsv", "genes": "genes.tsv"}}), encoding="utf-8")
        self.put("dnds.tsv")
        owned = self.earlier(report=str(rp))
        self.assertEqual(owned.known, {"dnds.tsv", "genes.tsv", "failed.tsv", "dnds.tsv.tmp", "genes.tsv.tmp",
                                       "failed.tsv.tmp"})
        self.assertFalse(owned.work_is_coal)               # (a report elsewhere does not make work/ COAL's)
        rp.write_text(json.dumps({"schema": SCHEMA, "out_dir": str(self.tmp / "other"),
                                  "outputs": {"dnds": "dnds.tsv"}}), encoding="utf-8")
        self.assertEqual(self.earlier(report=str(rp)).known, set())
        with self.assertRaisesRegex(OutputsRefused, "already holds dnds.tsv"):
            self.earlier(report=str(rp)).check_foreign(PLANNED, report=str(rp))


class ForeignOutputs(Folder):
    def test_files_no_run_recorded_are_refused_unless_overwrite(self):
        for n in ("genes.tsv", "dnds.tsv", "report.json"):
            self.put(n, "mine\n")
        self.put("dnds.tsv.tmp")                           # temporary names: COAL's own, not checked
        with self.assertRaises(OutputsRefused) as cm:
            self.earlier().check_foreign(PLANNED)
        self.assertIn(f"{self.out} already holds genes.tsv, dnds.tsv and report.json, which this run would replace "
                      "and which no COAL run wrote there (work/outputs.json does not list them)", str(cm.exception))
        self.earlier().check_foreign(PLANNED, overwrite=True)
        self.symlink(self.tmp / "nowhere", self.out / "failed.tsv")     # a link counts as such a file
        with self.assertRaisesRegex(OutputsRefused, "already holds failed.tsv, genes.tsv"):
            self.earlier().check_foreign(PLANNED)

    def test_recorded_and_read_in_place_files_go_on(self):
        self.put("orthologs.tsv", "mine\n")
        self.earlier().check_foreign(["orthologs.tsv"], pairs_in_place=True)    # --pairs DIR/orthologs.tsv
        self.write_record({"coal_version": "0.6.0", "files": ["dnds.tsv"], "read_in_place": ["orthologs.tsv"]})
        self.put("dnds.tsv")
        self.earlier().check_foreign(["dnds.tsv", "orthologs.tsv"])
        self.put("clusters.tsv")                           # an output this run does not write: no reason to refuse
        self.earlier().check_foreign(["dnds.tsv"])


class TamperedRecords(Folder):
    def outside(self) -> dict:
        """Files outside the output folder, and a link inside it to each."""
        victims = {}
        for n in ("victim.txt", "victim.fasta"):
            p = self.tmp / n
            p.write_text("keep me\n", encoding="utf-8")
            victims[n] = p
        return victims

    def test_names_in_the_record_cannot_escape_the_folder(self):
        victims = self.outside()
        evil = ["../victim.txt", str(victims["victim.txt"]), "../victim.fasta", "work/outputs.json", "..", "",
                "work", "inputs/x", 5, None]
        self.write_record({"coal_version": "0.6.0", "files": evil, "read_in_place": evil,
                           "codon_alignments_dir": True, "codon_alignments_files": evil})
        self.put(f"{CODON}/c1.fasta")
        self.put("work/inputs/x")
        owned = self.earlier()
        self.assertEqual(owned.known, set())
        owned.clear()
        owned.new_record(PLANNED)
        self.assertEqual(owned.record["read_in_place"], [])
        for p in victims.values():
            self.assertEqual(p.read_text(encoding="utf-8"), "keep me\n")
        self.assertTrue((self.out / "work" / "outputs.json").is_file())
        self.assertTrue((self.out / "work" / "inputs" / "x").is_file())
        self.assertEqual(self.names(CODON), ["c1.fasta"])  # not in its (tampered) list: the user's

    def test_a_tampered_report_cannot_name_other_files(self):
        victims = self.outside()
        self.put("report.json", json.dumps({"schema": SCHEMA, "mode": "clusters", "outputs": {
            "a": "../victim.txt", "b": str(victims["victim.fasta"]), "c": "work/state.json", "d": ["dnds.tsv"],
            "dnds": "dnds.tsv"}}))
        self.put("dnds.tsv")
        self.put("work/state.json")
        owned = self.earlier()
        self.assertEqual(owned.known, {"dnds.tsv", "dnds.tsv.tmp", "report.json", "report.json.tmp", "failed.tsv",
                                       "failed.tsv.tmp"})
        owned.clear()
        self.assertEqual(sorted(p.name for p in self.out.iterdir()), ["work"])
        self.assertTrue((self.out / "work" / "state.json").is_file())
        for p in victims.values():
            self.assertEqual(p.read_text(encoding="utf-8"), "keep me\n")

    @unittest.skipUnless(hasattr(os, "symlink"), "symbolic links")
    def test_links_are_removed_never_followed(self):
        victims = self.outside()
        self.symlink(victims["victim.txt"], self.out / "dnds.tsv")
        elsewhere = self.tmp / "elsewhere"
        elsewhere.mkdir()
        (elsewhere / "theirs.fasta").write_text("keep me\n", encoding="utf-8")
        self.symlink(elsewhere, self.out / CODON)          # an alignment folder that is a link: never COAL's
        self.write_record({"coal_version": "0.6.0", "files": ["dnds.tsv"], "read_in_place": [],
                           "codon_alignments_dir": True})
        owned = self.earlier()
        self.assertEqual(owned.dirs, set())
        owned.clear()
        self.assertFalse(os.path.lexists(self.out / "dnds.tsv"))        # the link went
        self.assertEqual(victims["victim.txt"].read_text(encoding="utf-8"), "keep me\n")
        self.assertEqual((elsewhere / "theirs.fasta").read_text(encoding="utf-8"), "keep me\n")
        with self.assertRaisesRegex(OutputsRefused, "exists and was not written by a COAL run"):
            self.earlier(writes=(CODON,))

    def test_inputs_and_kept_files_are_not_removed(self):
        self.write_record({"coal_version": "0.6.0", "files": ["clusters.tsv", "dnds.tsv"], "read_in_place": []})
        self.put("clusters.tsv")
        self.put("dnds.tsv")
        self.start(keep={"clusters.tsv"})                  # --resume keeps clusters.tsv
        self.assertTrue((self.out / "clusters.tsv").is_file())
        self.write_record({"coal_version": "0.6.0", "files": ["clusters.tsv", "dnds.tsv"], "read_in_place": []})
        self.put("dnds.tsv")
        owned = self.earlier()
        owned.clear(is_input=lambda p: os.path.basename(p) == "dnds.tsv")   # one of this run's inputs
        self.assertTrue((self.out / "dnds.tsv").is_file())
        self.assertFalse((self.out / "clusters.tsv").exists())


class OlderRecords(Folder):
    def test_a_record_of_0_6_0_clears_the_engines_fasta(self):
        """0.6.0's work/outputs.json has no per-folder lists: every .fasta with the engine's lower-case extension is
        the earlier run's (as in 0.6.0), any other file there stays."""
        self.write_record({"coal_version": "0.6.0", "files": sorted(PLANNED), "read_in_place": [],
                           "codon_alignments_dir": True, "protein_alignments_dir": False})
        for n in ("c1.fasta", "c2.fasta", "c3.fasta.tmp", "Mine.FASTA", "notes.txt"):
            self.put(f"{CODON}/{n}")
        owned = self.earlier()
        self.assertEqual(owned.foreign_alignment_files(CODON), ["Mine.FASTA"])
        owned.clear()
        self.assertEqual(self.names(CODON), ["Mine.FASTA", "notes.txt"])

    def test_a_report_without_a_record_names_the_folder(self):
        """work/ deleted (or a version before work/outputs.json): the folder's report.json lists the outputs, and a
        cluster report with codon_alignments names that folder."""
        self.put("report.json", json.dumps({"schema": SCHEMA, "mode": "clusters", "outputs": {
            "dnds": "dnds.tsv", "codon_alignments": "codon_alignments"}}))
        self.put(f"{CODON}/c1.fasta")
        self.put(f"{PROT}/p1.fasta")                       # not named by the report: not COAL's
        owned = self.earlier()
        self.assertTrue(owned.work_is_coal)
        self.assertEqual(owned.dirs, {CODON})
        owned.clear()
        self.assertFalse((self.out / CODON).exists())
        self.assertEqual(self.names(PROT), ["p1.fasta"])
        self.put("report.json", json.dumps({"schema": "other/1", "outputs": {"dnds": "dnds.tsv"}}))
        self.assertEqual(self.earlier().known, set())

    def test_a_development_record_without_read_in_place_keeps_the_reports_list(self):
        self.write_record({"coal_version": "0.6.0.dev", "files": ["dnds.tsv"]})
        self.put("report.json", json.dumps({"schema": SCHEMA, "outputs": {"orthologs": "orthologs.tsv"}}))
        self.put("orthologs.tsv")
        owned = self.earlier()
        self.assertEqual((owned.known, owned.listed), ({"dnds.tsv"}, {"orthologs.tsv"}))
        owned.check_foreign(["orthologs.tsv"])             # the report says it is COAL's


class CaseOfFasta(Folder):
    def test_fasta_in_another_case_is_the_users(self):
        """#98 (audit of #96, L5): on a case-insensitive file system the engine's x.fasta would replace X.FASTA."""
        self.cluster_run()
        for n in ("Mine.FASTA", "other.Fasta"):
            self.put(f"{CODON}/{n}")
        self.assertEqual(self.earlier().foreign_alignment_files(CODON), ["Mine.FASTA", "other.Fasta"])
        with self.assertRaisesRegex(OutputsRefused, "holds Mine.FASTA, other.Fasta, which no COAL run"):
            self.earlier(writes=(CODON,))
        rec = self.record_on_disk()
        del rec["codon_alignments_files"]
        self.write_record(rec)                             # without a list: the same
        with self.assertRaisesRegex(OutputsRefused, "holds Mine.FASTA"):
            self.earlier(writes=(CODON,))
        self.start()
        self.assertEqual(self.names(CODON), ["Mine.FASTA", "other.Fasta"])
        # what the engine writes is listed in its own lower case only
        owned = self.start(writes=(PROT,))
        self.engine(owned, {PROT: ["a.fasta", "B.FASTA"]})
        self.assertEqual(owned.record["protein_alignments_files"], ["a.fasta"])


@needs_cpp
class AlignmentFolderNotMade(PipelineCase):
    """The same through the pipeline: os.makedirs of protein_alignments/ fails (EACCES) in the engine step, after
    engine_starts() dropped both lists."""

    def clusters(self) -> Path:
        ids = list(read_fasta(FAM))
        tsv = self.tmp / "clusters.tsv"
        tsv.write_text("".join(f"{ids[k]}\t{ids[k]}\n{ids[k]}\t{ids[k + 1]}\n" for k in range(0, 8, 2)),
                       encoding="utf-8")
        return tsv

    def test_the_record_lists_nothing_for_a_folder_that_was_not_made(self):
        tsv = self.clusters()
        out = self.tmp / "c"
        real = os.makedirs

        def makedirs(name, *a, **k):
            if os.path.basename(os.fspath(name)) == PROT:
                raise PermissionError(errno.EACCES, "Permission denied", os.fspath(name))
            return real(name, *a, **k)
        with mock.patch.object(P.os, "makedirs", makedirs), redirect_stderr(io.StringIO()) as err:
            res = P.run_clusters(FAM, out, clusters=tsv, keep_alignments=True, keep_protein_alignments=True)
        self.assertEqual(res.exit_code, 1, err.getvalue())
        self.assertIn("Permission denied", err.getvalue())
        rec = json.loads((out / "work" / "outputs.json").read_text(encoding="utf-8"))
        self.assertTrue(rec["protein_alignments_dir"])
        self.assertEqual(rec["protein_alignments_files"], [])
        self.assertFalse((out / PROT).exists())
        # a protein_alignments/ the user makes afterwards is not the run's
        (out / PROT).mkdir()
        (out / PROT / "mine.fasta").write_text(">mine\nMK\n", encoding="utf-8")
        with redirect_stderr(io.StringIO()) as err:
            res = P.run_clusters(FAM, out, clusters=tsv)
        self.assertEqual(res.exit_code, 2, err.getvalue())     # (one record of the fixture is not usable)
        self.assertEqual(sorted(os.listdir(out / PROT)), ["mine.fasta"])


class NoLinkIsWrittenThrough(Folder):
    """L1 of AUDIT-99: work/outputs.json.tmp was opened with open(..., "w"), which follows a symbolic link: a link
    planted there made the next run write its record into the link's target, outside the folder. Every file COAL
    writes in the folder is now created anew (open_new), a log is appended to only when it is no link
    (open_append), and what a program writes is removed before it runs (remove_file)."""

    def victim(self, name="victim.txt") -> Path:
        v = self.tmp / name
        v.write_text("keep me\n", encoding="utf-8")
        return v

    def test_the_record_is_not_written_through_a_link(self):
        v = self.victim()
        owned = self.start()
        self.symlink(v, self.out / "work" / "outputs.json.tmp")
        owned.write()
        self.assertEqual(v.read_text(encoding="utf-8"), "keep me\n")
        self.assertFalse(os.path.islink(self.out / "work" / "outputs.json"))
        self.assertEqual(self.record_on_disk(), owned.record)
        os.link(v, self.out / "work" / "outputs.json.tmp")   # a hard link: removed, not written into
        owned.write()
        self.assertEqual(v.read_text(encoding="utf-8"), "keep me\n")

    def test_open_new(self):
        v = self.victim()
        p = self.out / "t.tmp"
        for make in (lambda: self.symlink(v, p), lambda: self.symlink(self.tmp / "nowhere", p),
                     lambda: os.link(v, p), lambda: p.write_text("old\n", encoding="utf-8")):
            make()
            with R.open_new(str(p)) as fh:
                fh.write("new\n")
            self.assertFalse(os.path.islink(p))
            self.assertEqual(p.read_bytes(), b"new\n")
            self.assertEqual(v.read_text(encoding="utf-8"), "keep me\n")
            p.unlink()
        with R.open_new(str(p), binary=True) as fh:
            fh.write(b"\x00\xff")
        self.assertEqual(p.read_bytes(), b"\x00\xff")

    def test_a_link_made_after_the_removal_is_not_followed(self):
        """O_EXCL (and O_NOFOLLOW): a link that appears between the removal and the open makes the open fail."""
        v = self.victim()
        p = self.out / "t.tmp"
        self.symlink(v, p)
        with mock.patch.object(R, "remove_file", lambda path: None), self.assertRaises(FileExistsError):
            R.open_new(str(p))
        self.assertEqual(v.read_text(encoding="utf-8"), "keep me\n")

    def test_open_append(self):
        v = self.victim()
        log = self.out / "work" / "logs" / "tool.log"
        log.parent.mkdir(parents=True)
        for k in range(2):                                 # a log of the folder's own grows
            with R.open_append(str(log)) as fh:
                fh.write(b"line %d\n" % k)
        self.assertEqual(log.read_bytes(), b"line 0\nline 1\n")
        log.unlink()
        for make in (lambda: self.symlink(v, log), lambda: os.link(v, log)):
            make()
            with R.open_append(str(log)) as fh:
                fh.write(b"mine\n")
            self.assertEqual(log.read_bytes(), b"mine\n")
            self.assertEqual(v.read_text(encoding="utf-8"), "keep me\n")
            log.unlink()

    def test_remove_file(self):
        v = self.victim()
        for make in (lambda p: self.symlink(v, p), lambda p: self.symlink(self.tmp / "nowhere", p),
                     lambda p: p.write_text("x", encoding="utf-8")):
            p = self.out / "x.tsv"
            make(p)
            R.remove_file(str(p))
            self.assertFalse(os.path.lexists(p))
        R.remove_file(str(self.out / "missing"))
        (self.out / "d").mkdir()
        R.remove_file(str(self.out / "d"))                 # a folder stays
        self.assertTrue((self.out / "d").is_dir())
        self.assertEqual(v.read_text(encoding="utf-8"), "keep me\n")


if __name__ == "__main__":
    unittest.main()
