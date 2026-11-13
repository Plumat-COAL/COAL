"""The example inputs in examples/: present, small, LF everywhere, and the engine gives the documented result on each."""
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX = os.path.join(ROOT, "examples")
NAMES = ["clean", "mito", "messy_ids", "one_broken_gene"]


def coal(*args):
    return subprocess.run([sys.executable, "-m", "coal", *args], capture_output=True, text=True, timeout=120,
                          stdin=subprocess.DEVNULL)


class TestExamples(unittest.TestCase):
    def test_files_present_and_small(self):
        for n in NAMES:
            d = os.path.join(EX, n)
            for f in ("protein.aln.fasta", "dna.fasta", "README.md"):
                self.assertTrue(os.path.isfile(os.path.join(d, f)), f"{n}/{f}")
            size = sum(os.path.getsize(os.path.join(d, f)) for f in os.listdir(d))
            self.assertLess(size, 50 * 1024, n)

    def test_clean_and_mito_code2_pass(self):
        with tempfile.TemporaryDirectory() as t:
            for n, extra in (("clean", []), ("mito", ["-codontable", "2"])):
                d = os.path.join(EX, n)
                r = coal(os.path.join(d, "protein.aln.fasta"), os.path.join(d, "dna.fasta"),
                         "-o", os.path.join(t, "out.fa"), *extra)
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertNotIn("warning", r.stderr)

    def test_mito_needs_code2(self):
        d = os.path.join(EX, "mito")
        with tempfile.TemporaryDirectory() as t:
            r = coal(os.path.join(d, "protein.aln.fasta"), os.path.join(d, "dna.fasta"),
                     "-o", os.path.join(t, "out.fa"))
        self.assertIn("does not correspond", r.stderr)

    def test_broken_gene_fails(self):
        d = os.path.join(EX, "one_broken_gene")
        with tempfile.TemporaryDirectory() as t:
            r = coal(os.path.join(d, "protein.aln.fasta"), os.path.join(d, "dna.fasta"),
                     "-o", os.path.join(t, "out.fa"))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("OG0001_sp14", r.stderr)


DNDS = os.path.join(EX, "dnds")


def cpp_bin():
    p = os.environ.get("COAL_CPP_BIN", os.path.join(ROOT, "build", "coal-cpp"))
    p = p if os.path.isabs(p) else os.path.join(ROOT, p)
    for cand in (p, p + ".exe", os.path.join(os.path.dirname(p), "Release", os.path.basename(p) + ".exe")):
        if os.path.isfile(cand) and (os.name != "nt" or cand.lower().endswith(".exe")):
            return cand
    return None


class TestDndsExample(unittest.TestCase):
    """examples/dnds: two synthetic species with 12 known ortholog pairs (examples/make_dnds_example.py), for
    `coal pipeline`, `coal dnds` and `coal kaks`."""
    FILES = ("README.md", "species_a_cds.fasta", "species_b_cds.fasta", "orthologs.tsv", "codon.fasta")

    def need_cpp(self):
        cpp = cpp_bin()
        if cpp is None:
            if os.environ.get("COAL_REQUIRE_CPP"):
                self.fail("build/coal-cpp (COAL_CPP_BIN) is not built (COAL_REQUIRE_CPP is set)")
            self.skipTest("build/coal-cpp (COAL_CPP_BIN) is not built")
        return cpp

    def coal(self, *args, cpp):
        env = {**os.environ, "COAL_CPP": cpp, "PYTHONPATH": os.path.join(ROOT, "src")}
        return subprocess.run([sys.executable, "-m", "coal", *args], capture_output=True, text=True, timeout=300,
                              stdin=subprocess.DEVNULL, env=env, encoding="utf-8")

    def test_files_present_small_and_packaged(self):
        self.assertEqual(sorted(os.listdir(DNDS)), sorted(self.FILES))
        self.assertLess(sum(os.path.getsize(os.path.join(DNDS, f)) for f in self.FILES), 60 * 1024)
        with open(os.path.join(ROOT, "pyproject.toml"), encoding="utf-8") as fh:
            proj = fh.read()
        for f in self.FILES:
            self.assertIn(f'"examples/dnds/{f}"', proj, "pyproject.toml does not install it under share/coal/examples")
        with open(os.path.join(DNDS, "README.md"), encoding="utf-8") as fh:
            readme = fh.read()
        self.assertIn("synthetic, not real genes", readme)

    def test_generator_reproduces_the_inputs(self):
        with tempfile.TemporaryDirectory() as t:
            script = os.path.join(t, "make_dnds_example.py")
            shutil.copyfile(os.path.join(EX, "make_dnds_example.py"), script)
            r = subprocess.run([sys.executable, script], capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr)
            for f in ("species_a_cds.fasta", "species_b_cds.fasta", "orthologs.tsv"):
                with open(os.path.join(t, "dnds", f), "rb") as a, open(os.path.join(DNDS, f), "rb") as b:
                    self.assertEqual(a.read(), b.read(), f"{f}: run examples/make_dnds_example.py")

    def readme_table(self):
        with open(os.path.join(DNDS, "README.md"), encoding="utf-8") as fh:
            rows = [ln.strip("|\n").split("|") for ln in fh if ln.startswith("| spA_g")]
        return [[c.strip() for c in r] for r in rows]

    def check_rows(self, rows):
        """12 pairs, all with an estimate, a range of dN/dS, and the numbers the README shows."""
        self.assertEqual(len(rows), 12)
        self.assertEqual({r["status"] for r in rows}, {"ok"})
        w = [float(r["dNdS"]) for r in rows]
        self.assertLess(min(w), 0.05)
        self.assertGreater(max(w), 1.0)
        self.assertTrue(all(0.05 < float(r["dS"]) < 1.1 for r in rows))
        for r, shown in zip(rows, self.readme_table()):
            self.assertEqual((r["seq1"], r["seq2"]), (shown[0], shown[1]))
            for col, text in zip(("dN", "dS", "dNdS"), shown[2:5]):
                self.assertAlmostEqual(float(r[col]) / float(text), 1.0, delta=0.006, msg=(r["seq1"], col))
            if "codons_compared" in r:
                self.assertEqual(r["codons_compared"], shown[5])

    def test_coal_dnds_and_kaks_on_the_codon_alignment(self):
        cpp = self.need_cpp()
        r = self.coal("dnds", os.path.join(DNDS, "codon.fasta"), "--pairs", "consecutive", cpp=cpp)
        self.assertEqual(r.returncode, 0, r.stderr)
        lines = r.stdout.splitlines()
        head = lines[0].split("\t")
        self.check_rows([dict(zip(head, ln.split("\t"))) for ln in lines[1:]])
        with tempfile.TemporaryDirectory() as t:
            out = os.path.join(t, "o.kaks")
            r = self.coal("kaks", "-i", os.path.join(DNDS, "codon.fasta"), "-o", out, "--pairs", "consecutive",
                          "-m", "YN", cpp=cpp)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(out, encoding="utf-8") as fh:
                kaks = [ln.split("\t") for ln in fh.read().splitlines()[1:]]
        self.assertEqual([k[0] for k in kaks], [ln.split("\t")[0] + "&" + ln.split("\t")[1] for ln in lines[1:]])
        for k, ln in zip(kaks, lines[1:]):
            self.assertAlmostEqual(float(k[4]), float(ln.split("\t")[5]), delta=1e-5 * max(1.0, float(k[4])))

    def test_pipeline_with_the_given_pairs_writes_codon_fasta(self):
        cpp = self.need_cpp()
        with tempfile.TemporaryDirectory() as t:
            out = os.path.join(t, "out")
            r = self.coal("pipeline", "pairs", "--query", os.path.join(DNDS, "species_a_cds.fasta"), "--subject",
                          os.path.join(DNDS, "species_b_cds.fasta"), "--pairs", os.path.join(DNDS, "orthologs.tsv"),
                          "-o", out, "--keep-alignments", "-q", cpp=cpp)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(os.path.join(out, "codon_alignments.fasta"), "rb") as a, \
                    open(os.path.join(DNDS, "codon.fasta"), "rb") as b:
                self.assertEqual(a.read(), b.read(), "examples/dnds/codon.fasta is not what the pipeline writes")
            with open(os.path.join(out, "dnds.tsv"), encoding="utf-8") as fh:
                lines = fh.read().splitlines()
            head = lines[0].split("\t")
            self.check_rows([dict(zip(head, ln.split("\t"))) for ln in lines[1:]])

    def test_pipeline_finds_the_orthologs_with_diamond(self):
        cpp = self.need_cpp()
        sys.path.insert(0, os.path.join(ROOT, "src"))
        from coal import dnds_pipeline as P
        if P.find_diamond() is None:
            if os.environ.get("COAL_REQUIRE_DIAMOND"):
                self.fail("DIAMOND was not found (COAL_REQUIRE_DIAMOND is set)")
            self.skipTest("DIAMOND was not found")
        with tempfile.TemporaryDirectory() as t:
            out = os.path.join(t, "out")
            # the README's first command, with a flag
            r = self.coal("pipeline", "pairs", "--query", os.path.join(DNDS, "species_a_cds.fasta"), "--subject",
                          os.path.join(DNDS, "species_b_cds.fasta"), "-o", out, "--min-codons", "50", cpp=cpp)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertNotIn("transcripts per gene", r.stderr)
            with open(os.path.join(out, "dnds.tsv"), encoding="utf-8") as fh:
                lines = fh.read().splitlines()
            head = lines[0].split("\t")
            rows = [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]
            with open(os.path.join(DNDS, "orthologs.tsv"), encoding="utf-8") as fh:
                truth = [tuple(ln.split("\t")[:2]) for ln in fh.read().splitlines()[1:]]
            self.assertEqual([(r_["seq1"], r_["seq2"]) for r_ in rows], truth)
            self.assertEqual([(r_["status"], r_["note"]) for r_ in rows if r_["status"] != "ok"],
                             [("warning", "fewer than 50 codons compared")])
            for r_ in rows:
                r_["status"] = "ok"
            self.check_rows(rows)
            with open(os.path.join(out, "unpaired.tsv"), encoding="utf-8") as fh:
                self.assertEqual(fh.read().splitlines()[1:], ["query\tspA_orphan\tno_hit", "subject\tspB_orphan\tno_hit"])


class TestExamplesAreLF(unittest.TestCase):
    def test_gitattributes_forces_lf_for_examples(self):
        with open(os.path.join(ROOT, ".gitattributes"), encoding="utf-8") as fh:
            rules = [ln.split() for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
        self.assertIn(["examples/**", "text", "eol=lf"], rules)

    def test_git_reports_lf_for_example_files(self):
        if not os.path.isdir(os.path.join(ROOT, ".git")) and not os.path.isfile(os.path.join(ROOT, ".git")):
            self.skipTest("not a git checkout")
        if shutil.which("git") is None:
            self.skipTest("git not installed")
        files = sorted(glob.glob(os.path.join(ROOT, "examples", "*", "*.fasta")))
        self.assertTrue(files)
        rel = [os.path.relpath(f, ROOT).replace(os.sep, "/") for f in files]
        out = subprocess.run(["git", "check-attr", "eol", "--"] + rel, cwd=ROOT, capture_output=True,
                             text=True, timeout=60).stdout
        for line in out.splitlines():
            self.assertTrue(line.endswith(": eol: lf"), line)

    def test_example_files_have_no_carriage_returns(self):
        files = glob.glob(os.path.join(ROOT, "examples", "*", "*.fasta"))
        self.assertTrue(files)
        for f in files:
            with open(f, "rb") as fh:
                self.assertNotIn(b"\r", fh.read(), f)


if __name__ == "__main__":
    unittest.main()
