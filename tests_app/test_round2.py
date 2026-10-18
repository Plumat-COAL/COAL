"""Round 2 review fixes (overnight/integration/CHANGES.md, R2-A*, R2-D*, R2-E*)."""
import os
import shutil
import tempfile
import unittest

from coal_app import core, tools

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FX = os.path.join(ROOT, "tests", "fixtures")


class RunsDirCase(unittest.TestCase):
    def setUp(self):
        old = core.RUNS_DIR
        core.RUNS_DIR = tempfile.mkdtemp(prefix="coal_test_runs_")
        self.addCleanup(shutil.rmtree, core.RUNS_DIR, True)
        self.addCleanup(setattr, core, "RUNS_DIR", old)

    def tmp(self):
        d = tempfile.mkdtemp(prefix="coal_test_r2_")
        self.addCleanup(shutil.rmtree, d, True)
        return d


class TestR2A1Nomismatch(RunsDirCase):  # R2-A1 / R2-D5 / R2-E1
    def test_mito_code1_nomismatch_is_not_all_ok(self):
        p, d = core.example_paths("mito")
        r = tools.codon_align(p, d, codontable=1, nomismatch=True)
        self.assertTrue(r["completed"])
        s = r["summary"]
        self.assertEqual(s["ok"], 0, s)
        self.assertEqual(s["warning"], 24, s)
        self.assertEqual(s["columns_removed"], 46, s)
        self.assertIn("46 codon columns were removed", r["summary_text"])
        self.assertTrue(any(h.get("suggested_codontable") == 2 for h in r["hints"]), r["hints"])
        self.assertNotIn("24 matched", r["methods_text"])
        self.assertIn("46 codon columns", r["methods_text"])
        for prob in r["problems"]:
            self.assertIn("removed", prob["explanation"])
        text = tools.explain_report(r["report_path"])
        self.assertNotIn("Nothing went wrong", text)
        self.assertNotEqual(core.banner_head(r["completed"], s), core.banner_head(True, {"ok": 24, "total": 24}))

    def test_probe_leaves_no_files(self):
        p, d = core.example_paths("mito")
        out = self.tmp()
        core.run_coal(p, d, codontable=1, nomismatch=True, out_dir=out)
        self.assertEqual(sorted(os.listdir(out)), ["protein.codon.fasta", "protein.report.json"])

    def test_output_unchanged_by_probe(self):
        p, d = core.example_paths("mito")
        a = core.run_coal(p, d, codontable=1, nomismatch=True, out_dir=self.tmp(), probe_codes=False)
        with open(a.output_path, encoding="utf-8") as fh:
            got = fh.read()
        import subprocess, sys
        ref = subprocess.run([sys.executable, "-m", "coal", p, d, "-nomismatch"], capture_output=True,
                             text=True, timeout=120).stdout
        self.assertEqual(got, ref)

    def test_table3_nomismatch_not_nothing_went_wrong(self):
        p, d = os.path.join(FX, "table3", "pep.fasta"), os.path.join(FX, "table3", "nuc.fasta")
        r = core.run_coal(p, d, codontable=1, nomismatch=True, out_dir=self.tmp())
        self.assertNotIn("Nothing went wrong", core.explain(r.report))
        self.assertGreater(r.report["summary"]["columns_removed"], 0)

    def test_long_mismatch_one_warning(self):  # R2-E1 repro
        p, d = os.path.join(FX, "long_mismatch", "pep.fasta"), os.path.join(FX, "long_mismatch", "nuc.fasta")
        r = core.run_coal(p, d, nomismatch=True, out_dir=self.tmp(), probe_codes=False)
        s = r.report["summary"]
        self.assertEqual((s["ok"], s["warning"], s["columns_removed"]), (1, 1, 1))
        bad = [x for x in r.report["sequences"] if x["status"] == "warning"]
        self.assertEqual(bad[0]["id"], "s2")
        self.assertEqual(bad[0]["messages"][0]["column"], 13)

    def test_clean_nomismatch_is_ok(self):
        p, d = core.example_paths("clean")
        r = core.run_coal(p, d, nomismatch=True, out_dir=self.tmp())
        s = r.report["summary"]
        self.assertEqual((s["ok"], s["warning"], s.get("columns_removed")), (8, 0, 0))
        self.assertIn("Nothing went wrong", core.explain(r.report))


class TestR2A2Notices(RunsDirCase):  # R2-A2 / R2-E2
    def write(self, pep, nuc):
        d = self.tmp()
        p, n = os.path.join(d, "p.fa"), os.path.join(d, "n.fa")
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(pep)
        with open(n, "w", encoding="utf-8") as fh:
            fh.write(nuc)
        return p, n

    def test_parse_no_anchor_notice_strips_colon(self):
        line = ("coal: warning: s0: no 10-residue segment of the protein was found in the DNA; codons were "
                "assigned by position (7 of 25 do not match)")
        sid, msg = core._parse_warning(line)
        self.assertEqual(sid, "s0")
        self.assertEqual(msg["code"], "no_anchor")

    def test_table3_notice_marks_gene(self):
        p, d = os.path.join(FX, "table3", "pep.fasta"), os.path.join(FX, "table3", "nuc.fasta")
        r = core.run_coal(p, d, codontable=1, out_dir=self.tmp(), probe_codes=False)
        s0 = next(x for x in r.report["sequences"] if x["id"] == "s0")
        self.assertEqual(s0["status"], "warning")
        self.assertTrue(any(m["code"] == "no_anchor" for m in s0["messages"]))
        self.assertIn("assigned by position", core.sequence_sentence(s0))
        self.assertNotIn("Nothing went wrong", core.explain(r.report))

    def test_partial_id_match_is_a_prominent_hint(self):
        p, n = self.write(">geneA\nMKVAL\n>geneB\nMKVAL\n", ">geneB\nATGAAAGTTGCTCTG\n>geneX\nATGAAGGTCGCACTA\n")
        r = core.run_coal(p, n, out_dir=self.tmp(), probe_codes=False)
        self.assertIn("1 of 2 protein IDs", r.stderr)
        codes = [h["code"] for h in r.hints]
        self.assertIn("partial_id_match", codes)
        h = next(h for h in r.hints if h["code"] == "partial_id_match")
        self.assertIn("geneA", h["text"])
        self.assertIn("geneB", h["text"])
        self.assertNotIn("looks right", " ".join(x["text"] for x in r.hints))
        text = core.explain(r.report)
        self.assertNotIn("Nothing went wrong", text)
        self.assertIn("1 of 2", text)
        self.assertIn("problems", core.banner_head(r.ok, r.report["summary"], r.hints))

    def test_unknown_warning_line_is_never_dropped(self):
        rep = core.build_fallback_report(*self.write(">a\nMK\n", ">a\nATGAAA\n"), "out.fa", "fasta", 1, False,
                                         False, False, 0, "coal: warning: something new happened here\n", 0.1)
        self.assertTrue(any("something new happened here" in h["text"] for h in rep["hints"]), rep["hints"])
        self.assertNotIn("Nothing went wrong", core.explain(rep))


class TestR2A6Stopped(RunsDirCase):  # R2-A6
    def report(self, code, stderr):
        p, d = core.example_paths("clean")
        return core.build_fallback_report(p, d, os.path.join(self.tmp(), "out.fa"), "fasta", 1, False, False, False,
                                          code, stderr, 0.1)

    def test_killed_or_interrupted_is_not_nothing_went_wrong(self):
        for code, err in ((143, ""), (-9, ""), (130, "coal: interrupted\n"), (1, "coal: interrupted\n"),
                          (137, "")):
            rep = self.report(code, err)
            text = core.explain(rep)
            self.assertNotIn("Nothing went wrong", text, code)
            self.assertIn("stopped before it finished", text, code)
            self.assertEqual(rep["summary"]["ok"], 0, code)
            self.assertIsNone(rep["output"])

    def test_timeout_is_not_nothing_went_wrong(self):
        rep = self.report(1, "coal: error: COAL did not finish within 5 seconds\n")
        self.assertNotIn("Nothing went wrong", core.explain(rep))

    def test_success_still_says_nothing_went_wrong(self):
        self.assertIn("Nothing went wrong", core.explain(self.report(0, "")))


class TestR2A7UsageError(RunsDirCase):  # R2-A7
    def test_exit_2_usage_is_not_success(self):
        p, d = core.example_paths("clean")
        rep = core.build_fallback_report(p, d, os.path.join(self.tmp(), "o.fa"), "fasta", 1, False, False, False,
                                         2, "usage: coal [-h] ...\ncoal: error: the following arguments are "
                                            "required: nuc_fasta\n", 0.1)
        self.assertEqual(rep["summary"]["ok"], 0)
        self.assertNotIn("Nothing went wrong", core.explain(rep))
        self.assertEqual(core.methods_text(rep), "")
        rr = core.RunResult(2, p, "", rep, "usage: coal", "", "")
        self.assertFalse(rr.ok)

    def test_file_names_starting_with_dash(self):
        d = self.tmp()
        p0, d0 = core.example_paths("clean")
        shutil.copy(p0, os.path.join(d, "-p.fa"))
        shutil.copy(d0, os.path.join(d, "-d.fa"))
        old = os.getcwd()
        os.chdir(d)
        self.addCleanup(os.chdir, old)
        r = tools.codon_align("-p.fa", "-d.fa")
        self.assertTrue(r["completed"], r)
        self.assertEqual(r["summary"]["ok"], 8, r)


class TestR2A8Prune(RunsDirCase):  # R2-A8
    def old(self, path):
        t = __import__("time").time() - 48 * 3600
        os.utime(path, (t, t))

    def test_only_coal_run_folders_are_pruned(self):
        mine = core.new_run_dir()
        user = os.path.join(core.RUNS_DIR, "thesis_data")
        os.makedirs(user)
        with open(os.path.join(user, "chapter1.fa"), "w", encoding="utf-8") as fh:
            fh.write(">x\nA\n")
        fake = os.path.join(core.RUNS_DIR, "run_not_ours")
        os.makedirs(fake)
        for p in (mine, user, fake):
            self.old(p)
        core.prune_runs()
        self.assertFalse(os.path.exists(mine))
        self.assertTrue(os.path.isfile(os.path.join(user, "chapter1.fa")))
        self.assertTrue(os.path.isdir(fake))

    def test_fresh_run_folder_kept(self):
        mine = core.new_run_dir()
        core.prune_runs()
        self.assertTrue(os.path.isdir(mine))


class TestR2A11PastedPaths(unittest.TestCase):  # R2-A11
    def setUp(self):
        d = tempfile.mkdtemp(prefix="coal test r2 ")
        self.addCleanup(shutil.rmtree, d, True)
        self.f = os.path.join(d, "my protein.fa")
        with open(self.f, "w", encoding="utf-8") as fh:
            fh.write(">a\nMK\n")

    def test_file_url_with_percent_escapes(self):
        from pathlib import Path
        url = Path(self.f).as_uri()
        self.assertIn("%20", url)
        self.assertEqual(core.resolve_input(None, url)[0], os.path.abspath(self.f))
        self.assertEqual(core.resolve_input(None, "  '" + url + "'  ")[0], os.path.abspath(self.f))

    def test_invisible_marks_and_quotes(self):
        for junk in ("\u202a", "\u200e", "\u200f", "\ufeff", " \u202a"):
            got = core.resolve_input(None, junk + '"' + self.f + '"' + junk + " \n")[0]
            self.assertEqual(got, os.path.abspath(self.f), repr(junk))

    def test_bad_path_still_rejected(self):
        with self.assertRaises(ValueError):
            core.resolve_input(None, "file:///no/such%20file.fa")


def _procs_mentioning(text):
    """PIDs of live processes whose command line contains `text` (Linux /proc)."""
    found = []
    for pid in os.listdir("/proc"):
        if not pid.isdigit() or int(pid) == os.getpid():
            continue
        try:
            with open(f"/proc/{pid}/cmdline", "rb") as fh:
                cmd = fh.read().decode("utf-8", "replace")
            with open(f"/proc/{pid}/stat", "rb") as fh:
                state = fh.read().decode().rsplit(")", 1)[1].split()[0]
        except OSError:
            continue
        if text in cmd and state != "Z":
            found.append(int(pid))
    return found


def _wait_gone(text, limit=10.0):
    import time
    end = time.monotonic() + limit
    left = _procs_mentioning(text)
    while left and time.monotonic() < end:
        time.sleep(0.2)
        left = _procs_mentioning(text)
    return left


@unittest.skipUnless(os.path.isdir("/proc/self") and os.name != "nt", "R2-A4: needs Linux /proc and POSIX signals")
class TestR2A4NoOrphans(unittest.TestCase):  # R2-A4
    @classmethod
    def setUpClass(cls):
        import subprocess, sys
        cls.data = tempfile.mkdtemp(prefix="coal_test_orphans_")
        r = subprocess.run([sys.executable, os.path.join(ROOT, "benchmarks", "gen_data.py"), "40000", "600",
                            cls.data], capture_output=True, text=True, timeout=300)
        assert r.returncode == 0, r.stderr
        cls.pep, cls.nuc = os.path.join(cls.data, "pep.fasta"), os.path.join(cls.data, "nuc.fasta")

    @classmethod
    def tearDownClass(cls):
        for pid in _procs_mentioning(cls.data):
            try:
                os.kill(pid, 9)
            except OSError:
                pass
        shutil.rmtree(cls.data, True)

    def test_app_timeout_leaves_no_processes_or_temp_files(self):
        out_dir = tempfile.mkdtemp(prefix="coal_test_orphan_out_")
        self.addCleanup(shutil.rmtree, out_dir, True)
        code, err, _ = core._run_engine(self.pep, self.nuc, os.path.join(out_dir, "o.fasta"), "fasta", 1,
                                        False, True, False, None, timeout=0.8)
        self.assertNotEqual(code, 0)
        self.assertIn("did not finish within", err)
        self.assertEqual(_wait_gone(self.data), [])
        self.assertEqual([n for n in os.listdir(out_dir) if n.startswith(".coal_")], [])

    def _compat(self, sig):
        import signal, subprocess, sys, time
        p = subprocess.Popen([sys.executable, "-m", "coal_compat", self.pep, self.nuc, "-output", "fasta",
                              "-nomismatch"], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        try:
            time.sleep(0.8)  # let the engine start its work (a full run takes ~2.5 s)
            p.send_signal(sig)
            err = p.communicate(timeout=60)[1].decode()
        finally:
            if p.poll() is None:
                p.kill()
                p.communicate(timeout=60)
        return p.returncode, err

    def test_compat_ctrl_c(self):
        import signal
        rc, err = self._compat(signal.SIGINT)
        self.assertEqual(rc, 130, err)
        self.assertNotIn("Traceback", err)
        self.assertIn("coal: interrupted", err)
        self.assertEqual(_wait_gone(self.data), [])

    def test_compat_sigterm_stops_the_engine(self):
        import signal
        rc, err = self._compat(signal.SIGTERM)
        self.assertEqual(rc, 128 + signal.SIGTERM, err)
        self.assertNotIn("Traceback", err)
        self.assertEqual(_wait_gone(self.data), [])


class TestR2A12ViewerNoFalseFailures(RunsDirCase):  # R2-A12
    def test_mito_code1_diagnose_calls_no_break(self):
        from coal_app import viewer
        p, d = core.example_paths("mito")
        dna = dict(viewer.read_alignment(d))
        for pid, row in viewer.read_alignment(p):
            ev = viewer.diagnose(row, dna[pid], 1)["events"]
            self.assertEqual(ev, [], pid)

    def test_standalone_viewer_runs_the_engine(self):
        import subprocess, sys
        p, d = core.example_paths("mito")
        out = os.path.join(self.tmp(), "v.html")
        r = subprocess.run([sys.executable, "-m", "coal_app.viewer", p, "--dna", d, "--codontable", "1", "-o", out],
                           capture_output=True, text=True, timeout=120,
                           env=dict(os.environ, COAL_RUNS_DIR=core.RUNS_DIR))
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(out, encoding="utf-8") as fh:
            page = fh.read()
        self.assertNotIn("Failed", page)
        self.assertNotIn("no longer spells", page)
        self.assertNotIn("COAL wrote no codon alignment", page)

    def test_engine_aligned_gene_is_never_failed(self):
        from coal_app import viewer
        p, d = core.example_paths("mito")
        rep = {"sequences": [{"id": pid, "status": "warning", "messages": []}
                             for pid, _ in viewer.read_alignment(p)]}
        rows, _, _ = viewer.build_rows(viewer.read_alignment(p), None, rep, viewer.read_alignment(d), 1)
        self.assertNotIn("failed", {r["status"] for r in rows})

    def test_real_frameshift_and_wrong_dna_still_found(self):
        from coal_app import viewer
        p, d = core.example_paths("one_broken_gene")
        prot, dna = dict(viewer.read_alignment(p)), dict(viewer.read_alignment(d))
        ev = viewer.diagnose(prot["OG0001_sp14"], dna["OG0001_sp14"])["events"]
        self.assertEqual(ev[0]["code"], "frameshift")
        import random
        rng = random.Random(1)
        row = prot["OG0001_sp01"]
        junk = dna["OG0001_sp01"][:60] + "".join(rng.choice("ACGT") for _ in range(len(dna["OG0001_sp01"]) - 60))
        ev = viewer.diagnose(row, junk)["events"]
        self.assertTrue(ev and ev[0]["code"] == "inconsistent", ev)


if __name__ == "__main__":
    unittest.main()
