"""The fixes of the new-user audit of 2026-11-14 (U01-U21, COAL side): messages a first-time user meets.

  * Units (nothing needed): the coal-cpp-missing text of a development version and of a release (U01), the banner
    of the development docs (U01), merged genetic codes 7 / 8 (U21), orthologr's option names (U21), --matrix
    checked by the parser (U08), the --pairs CSV message (U11), Windows paths (U17), the interrupt line (U13).
  * With build/coal-cpp (COAL_CPP_BIN): --pairs / --clusters that name no record stop with exit 1, a header line
    is skipped (U07); the first unusable record on stderr (U12); a genome as --query (U04); the notice when -o
    holds an earlier run (U16); `coal translate`'s genetic-code hint (U10); "kept in the alignment" only when an
    output was written, in both engines (U20); `coal dnds` / `coal kaks` on an input without pairs exit 2 (fuzz
    finding K05).
"""
from __future__ import annotations

import io
import os
import random
import re
import subprocess
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
from coal import cli  # noqa: E402
from coal import dnds_pipeline as P  # noqa: E402
from test_dnds_pipeline import CLUSTALO, CPP, FAM, FAM_CLUSTERS, PLAIN, Q, S, needs_cpp  # noqa: E402

BANNER_FILES = [ROOT / "README.md", ROOT / "docs" / "install.md", ROOT / "docs" / "tutorials" / "README.md"]
BANNER_MARK = "dev-version banner"
EXAMPLES = ROOT / "examples" / "clean"


def coal(*args, env=None, cwd=None):
    e = {**os.environ, "PYTHONPATH": str(ROOT / "src"), **({"COAL_CPP": str(CPP)} if CPP else {}), **(env or {})}
    return subprocess.run([sys.executable, "-m", "coal", *map(str, args)], capture_output=True, text=True,
                          timeout=300, stdin=subprocess.DEVNULL, env=e, encoding="utf-8", cwd=cwd)


def pipeline_main(argv) -> tuple[int, str]:
    """P.main(argv) with stderr captured: (exit status, stderr); argparse's SystemExit is the status."""
    err = io.StringIO()
    with redirect_stderr(err):
        try:
            rc = P.main(argv)
        except SystemExit as e:
            rc = e.code
    return rc, err.getvalue()


class TestInstallTexts(unittest.TestCase):
    def test_needs_cpp_of_a_development_version_does_not_send_users_to_the_release(self):
        dev = cli.needs_cpp_text("0.6.0.dev0")
        for must in ("0.6.0.dev0", "not in any release yet", "0.5.0", "cmake -S . -B build", "COAL_CPP",
                     "docs/install.md", "coal translate", "coal pipeline"):
            self.assertIn(must, dev)
        self.assertNotIn("releases/latest", dev)
        rel = cli.needs_cpp_text("0.6.0")
        for must in ("coal_align-<version>-py3-none-<platform>.whl", "https://github.com/plumat-coal/COAL/releases/latest",
                     "docs/install.md"):
            self.assertIn(must, rel)
        self.assertEqual(cli.NEEDS_CPP, cli.needs_cpp_text(__version__))

    def test_the_development_banner_is_there_while_the_version_is_one(self):
        """README, install.md and the tutorials' README say at the top that they describe the development
        version and how to install it, until a release: then the banner must go (this test fails)."""
        for p in BANNER_FILES:
            text = p.read_text(encoding="utf-8")
            if ".dev" in __version__:
                self.assertIn(BANNER_MARK, text, p)
                block = text[text.index(BANNER_MARK):][:1500]
                self.assertIn(__version__, block, p)
                for must in ("0.5.0", "git clone https://github.com/plumat-coal/COAL", "python -m pip install -e .",
                             "cmake -S . -B build -DCMAKE_BUILD_TYPE=Release", "export COAL_CPP=$PWD/build/coal-cpp"):
                    self.assertIn(must, block, p)
                self.assertLess(text.index(BANNER_MARK), 1200, f"{p}: the banner belongs at the top")
            else:
                self.assertNotIn(BANNER_MARK, text, f"{p}: {__version__} is a release: remove the banner")

    def test_install_from_source_builds_coal_cpp_and_says_how_to_keep_it(self):
        text = (ROOT / "docs" / "install.md").read_text(encoding="utf-8")
        sec = text[text.index("## Option 3: from source"):]
        sec = sec[:sec.index("\n## ", 5)]
        for must in ("cmake -S . -B build", "cmake --build build", "export COAL_CPP", "~/.bashrc",
                     'cp build/coal-cpp "$(dirname "$(command -v python)")/"'):
            self.assertIn(must, sec)


class TestOptionMessages(unittest.TestCase):
    def test_merged_genetic_codes(self):
        self.assertIn("NCBI merged code 7 into 4", cli.merged_code_error(["-codontable", "7"], ("-codontable",)))
        self.assertIn("use --codontable 1", cli.merged_code_error(["--codontable=8"], ("--codontable",)))
        self.assertIsNone(cli.merged_code_error(["-codontable", "11"], ("-codontable",)))
        self.assertIsNone(cli.merged_code_error(["x7.fa", "7"], ("-c",)))
        r = coal(EXAMPLES / "protein.aln.fasta", EXAMPLES / "dna.fasta", "-codontable", "7")
        self.assertEqual(r.returncode, 2)
        self.assertIn("use -codontable 4", r.stderr)
        rc, err = pipeline_main(["pairs", "--query", str(Q), "--subject", str(S), "-o", "x", "-c", "8"])
        self.assertEqual(rc, 2)
        self.assertIn("coal pipeline pairs: error: -c 8: NCBI genetic code 8 no longer exists", err)

    def test_orthologr_option_names_and_typos(self):
        rc, err = pipeline_main(["pairs", "--query", str(Q), "--subject", str(S), "-o", "x",
                                 "--sensitivity-mode", "default"])
        self.assertEqual(rc, 2)
        self.assertIn("unrecognized arguments: --sensitivity-mode default (--sensitivity-mode is orthologr's "
                      "sensitivity_mode: here it is --sensitivity)", err)
        rc, err = pipeline_main(["pairs", "--query", str(Q), "--subject", str(S), "-o", "x", "--sensitivty", "fast"])
        self.assertIn("(did you mean --sensitivity?)", err)
        rc, err = pipeline_main(["pairs", "--query", str(Q), "--subject", str(S), "-o", "x", "--zzz"])
        self.assertIn("unrecognized arguments: --zzz\n", err)

    def test_codontable_spelled_as_in_pal2nal(self):
        a = P.build_parser().parse_args(["pairs", "--query", "q", "--subject", "s", "-o", "o", "-codontable", "11"])
        self.assertEqual(a.codontable, 11)

    def test_matrix_is_checked_by_the_parser(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "out"
            rc, err = pipeline_main(["pairs", "--query", str(Q), "--subject", str(S), "-o", str(out),
                                     "--matrix", "BLOSUM26"])
            self.assertEqual(rc, 2)
            self.assertIn("argument --matrix: unknown matrix 'BLOSUM26' (choose from BLOSUM62", err)
            self.assertNotIn("codon-dnds", err)
            self.assertFalse(out.exists())   # nothing ran
            with self.assertRaisesRegex(ValueError, "unknown matrix 'BLOSUM26'"):
                P.run_pairs(Q, S, out, matrix="BLOSUM26")
        self.assertEqual(P._matrix_name("blosum45"), "blosum45")

    def test_pairs_file_in_csv(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "pairs.csv"
            p.write_text("query,subject\ng01,h01\n", encoding="utf-8")
            with self.assertRaises(P.PipelineError) as cm:
                P._read_pairs_file(str(p))
            msg = str(cm.exception)
            self.assertIn("line 1: expected two tab-separated columns (query ID, subject ID), found 1 column "
                          '"query,subject"', msg)
            self.assertIn("is this a CSV file? Convert it with tr ',' '\\t' < pairs.csv > pairs.tsv", msg)

    @unittest.skipIf(os.name == "nt", "a Windows path is a path there")
    def test_windows_path_on_linux(self):
        msg = P._cannot_open("C:\\Users\\me\\a.fasta")
        self.assertEqual(msg, "cannot open 'C:\\Users\\me\\a.fasta': no such file. This is a Windows path; under WSL "
                              "(Windows Subsystem for Linux) the same file is /mnt/c/Users/me/a.fasta")
        self.assertEqual(P._cannot_open("/x/y.fa"), "cannot open '/x/y.fa': no such file")
        # normalised as seqio / coal-cpp (#88): 'C:\\' is /mnt/c, doubled or mixed separators collapse (follow-up audit)
        for path, wsl in (("C:\\", "/mnt/c"), ("C:\\\\data\\\\x.fa", "/mnt/c/data/x.fa"),
                          ("D:/a//b\\", "/mnt/d/a/b")):
            self.assertEqual(P._cannot_open(path), f"cannot open '{path}': no such file. This is a Windows path; under "
                                                   f"WSL (Windows Subsystem for Linux) the same file is {wsl}")

    def test_interrupt_says_how_to_go_on(self):
        with mock.patch.object(P, "_execute", side_effect=KeyboardInterrupt):
            rc, err = pipeline_main(["pairs", "--query", str(Q), "--subject", str(S), "-o", "OUTDIR"])
        self.assertEqual(rc, 130)
        self.assertIn("coal pipeline: interrupted; the tables in OUTDIR are incomplete. Run the same command with "
                      "--resume", err)


@needs_cpp
class TestPipelineMessages(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(self._td.name)
        self._env = mock.patch.dict(os.environ, {"COAL_CPP": str(CPP)})
        self._env.start()

    def tearDown(self):
        self._env.stop()
        self._td.cleanup()

    def run_main(self, *argv) -> tuple[int, str]:
        return pipeline_main([str(a) for a in argv])

    def pairs(self, name, lines) -> Path:
        p = self.tmp / name
        p.write_text("".join(lines), encoding="utf-8")
        return p

    def test_pairs_that_name_no_record_stop_the_run(self):
        pf = self.pairs("wrong.tsv", [f"{k}\t{k * 3}\n" for k in range(1, 6)])
        rc, err = self.run_main("pairs", "--query", Q, "--subject", S, "-o", self.tmp / "o", "--pairs", pf)
        self.assertEqual(rc, 1)
        self.assertIn("error: none of the 5 pairs of --pairs", err)
        self.assertIn("can be made, so nothing was computed. 5 of the 5 query IDs of --pairs", err)
        self.assertIn("e.g. '1' in wrong.tsv, while the sequences are named like", err)
        self.assertIn("--pairs needs two tab-separated columns", err)
        self.assertFalse((self.tmp / "o" / "dnds.tsv").exists())

    def test_swapped_columns_stop_the_run(self):
        pf = self.pairs("swapped.tsv", [f"{b}\t{a}\n" for a, b in PLAIN])
        rc, err = self.run_main("pairs", "--query", Q, "--subject", S, "-o", self.tmp / "o", "--pairs", pf)
        self.assertEqual(rc, 1)
        self.assertIn("so nothing was computed: the columns look swapped", err)

    def test_a_header_line_is_not_a_pair(self):
        plain = self.pairs("plain.tsv", [f"{a}\t{b}\n" for a, b in PLAIN])
        head = self.pairs("head.tsv", ["gene_a\tgene_b\textra\n"] + [f"{a}\t{b}\n" for a, b in PLAIN])
        rc0, err = self.run_main("pairs", "--query", Q, "--subject", S, "-o", self.tmp / "a", "--pairs", plain, "-j", 1)
        self.assertIn(rc0, (0, 2), err)   # (2: the fixture has records that are not usable)
        self.assertNotIn("column header", err)
        rc, err = self.run_main("pairs", "--query", Q, "--subject", S, "-o", self.tmp / "b", "--pairs", head, "-j", 1)
        self.assertEqual(rc, rc0, err)
        self.assertIn("notice: line 1 of --pairs", err)
        self.assertIn("(gene_a gene_b) names no record of --query or --subject: it is read as a column header", err)
        self.assertEqual((self.tmp / "a" / "dnds.tsv").read_bytes(), (self.tmp / "b" / "dnds.tsv").read_bytes())
        # one pair that is no record (a gene missing on one side) still runs, as before: a row with its error
        some = self.pairs("some.tsv", [f"{a}\t{b}\n" for a, b in PLAIN] + ["nope\th01\n"])
        rc, err = self.run_main("pairs", "--query", Q, "--subject", S, "-o", self.tmp / "c", "--pairs", some, "-q")
        self.assertEqual(rc, 2)
        self.assertIn("nope", (self.tmp / "c" / "dnds.tsv").read_text(encoding="utf-8"))

    def test_clusters_that_name_no_record_stop_the_run(self):
        cf = self.pairs("clusters.tsv", ["rep1\tmem1\n", "rep1\tmem2\n"])
        rc, err = self.run_main("clusters", "--cds", FAM, "-o", self.tmp / "o", "--clusters", cf)
        self.assertEqual(rc, 1)
        self.assertIn("none of the 3 IDs of --clusters", err)
        self.assertIn("is a record of --cds, so nothing was computed", err)
        self.assertIn("e.g. 'rep1' in clusters.tsv", err)
        if not CLUSTALO:   # (the clusters route aligns each cluster with Clustal Omega: none on Windows / macOS CI)
            self.skipTest("the run with clusters that name records needs Clustal Omega")
        rc, err = self.run_main("clusters", "--cds", FAM, "-o", self.tmp / "ok", "--clusters", FAM_CLUSTERS, "-q")
        self.assertIn(rc, (0, 2), err)

    def test_the_first_unusable_record_is_named(self):
        text = Q.read_text(encoding="utf-8")
        first = re.search(r"^>(\S+)", text, re.M).group(1)
        bad = self.tmp / "q.fasta"
        bad.write_text(text.replace(">" + first, ">gène01", 1), encoding="utf-8")
        pf = self.pairs("p.tsv", [f"{a}\t{b}\n" for a, b in PLAIN])
        rc, err = self.run_main("pairs", "--query", bad, "--subject", S, "-o", self.tmp / "o", "--pairs", pf)
        self.assertIn("warning: ", err)
        line = next(x for x in err.splitlines() if "not usable (first:" in x)
        self.assertRegex(line, r"warning: \d+ records? not usable \(first: [^:]+: ")
        self.assertIn("gène01: non-ASCII character in the sequence ID", err)
        self.assertIn("failed.tsv", line)

    def test_a_genome_as_query(self):
        rnd = random.Random(7)
        genome = self.tmp / "genome.fna"
        seq = "".join(rnd.choice("ACGT") for _ in range(150_000))
        genome.write_text(">NC_1 a chromosome\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n",
                          encoding="utf-8")
        pf = self.pairs("p.tsv", ["NC_1\th01\n"])
        rc, err = self.run_main("pairs", "--query", genome, "--subject", S, "-o", self.tmp / "o", "--pairs", pf)
        self.assertIn("warning: --query ", err)
        self.assertIn(": 1 record of 150,000 nt with ", err)
        self.assertIn("this looks like a genome, chromosomes or contigs, not coding sequences", err)
        self.assertIn("_cds_from_genomic.fna", err)
        pf = self.pairs("p2.tsv", [f"{a}\t{b}\n" for a, b in PLAIN])
        rc, err = self.run_main("pairs", "--query", Q, "--subject", S, "-o", self.tmp / "o2", "--pairs", pf)
        self.assertNotIn("looks like a genome", err)

    def test_notice_when_the_output_folder_holds_an_earlier_run(self):
        pf = self.pairs("p.tsv", [f"{a}\t{b}\n" for a, b in PLAIN])
        out = self.tmp / "o"
        args = ("pairs", "--query", Q, "--subject", S, "-o", out, "--pairs", pf)
        rc, err = self.run_main(*args)
        self.assertNotIn("holds the results of an earlier run", err)
        first = (out / "dnds.tsv").read_bytes()
        rc, err = self.run_main(*args)
        self.assertIn(f"notice: {out} holds the results of an earlier run (completed", err)
        self.assertIn("Use another -o to keep them", err)
        self.assertEqual((out / "dnds.tsv").read_bytes(), first)
        rc, err = self.run_main(*args, "--resume")
        self.assertNotIn("holds the results of an earlier run", err)
        (self.tmp / "mine").mkdir()
        (self.tmp / "mine" / "notes.txt").write_text("x", encoding="utf-8")
        rc, err = self.run_main("pairs", "--query", Q, "--subject", S, "-o", self.tmp / "mine", "--pairs", pf)
        self.assertNotIn("holds the results of an earlier run", err)   # a folder of the user's own: nothing of COAL's is replaced
        self.assertTrue((self.tmp / "mine" / "notes.txt").is_file())


@needs_cpp
class TestTranslateHint(unittest.TestCase):
    def cds(self, path: Path, n: int = 4) -> None:
        """n genes whose TGG codons (Trp) are TGA: internal stops under code 1, clean under code 4."""
        from coal.tables import get_code
        recs = []
        rnd = random.Random(3)
        sense = sorted(c for c, aa in get_code(1).code.items() if aa not in "*W")
        for k in range(n):
            codons = ["ATG"] + [rnd.choice(sense) if j % 9 else "TGA" for j in range(1, 120)] + ["TAA"]
            recs.append(f">gene{k}\n{''.join(codons)}\n")
        path.write_text("".join(recs), encoding="ascii")

    def test_wrong_code_gets_a_hint(self):
        with tempfile.TemporaryDirectory() as td:
            f = Path(td) / "cds.fasta"
            self.cds(f)
            r = coal("translate", f, "-o", Path(td) / "c1.faa")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("coal translate: 4 gene(s): 0 ok, 4 warning, 0 failed", r.stderr)
            self.assertIn("coal translate: warning: 4 of the 4 genes have internal stop codons under genetic code 1 "
                          "(first: gene0, 13 internal stops), and 4 of the 4 with two or more have none under", r.stderr)
            self.assertIn("--genes FILE lists every gene", r.stderr)
            r = coal("translate", f, "-o", Path(td) / "c1.faa", "--genes", Path(td) / "g.tsv")
            self.assertIn(f"Every gene is in {Path(td) / 'g.tsv'} (column internal_stops)", r.stderr)
            for quiet in (["-c", "4"], ["-q"]):
                r = coal("translate", f, "-o", Path(td) / "c4.faa", *quiet)
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertNotIn("coal translate: warning", r.stderr)
            r = coal("translate", Q, "-o", Path(td) / "q.faa")   # one gene with stops among many: no line
            self.assertNotIn("coal translate: warning", r.stderr)


@needs_cpp
class TestNoPair(unittest.TestCase):
    def test_coal_dnds_and_kaks_on_an_input_without_pairs_exit_2(self):
        """Fuzz finding K05 through `coal` (which hands the command to coal-cpp), also for a gzip input."""
        import gzip
        with tempfile.TemporaryDirectory() as t:
            plain, gz = Path(t) / "empty.axt", Path(t) / "empty.axt.gz"
            plain.write_text("\n\n", encoding="utf-8")
            gz.write_bytes(gzip.compress(b"# a comment\n"))
            for f in (plain, gz):
                r = coal("dnds", f, "-q")
                self.assertEqual(r.returncode, 2, r.stderr)
                self.assertIn("coal dnds: error: the input holds no pair of sequences", r.stderr)
                r = coal("kaks", "-i", f, "-o", Path(t) / "o.kaks")
                self.assertEqual(r.returncode, 2, r.stderr)


class TestKeptInTheAlignment(unittest.TestCase):
    def inputs(self, td: Path) -> tuple[Path, Path]:
        """examples/clean with one mismatching codon in gene01 and a gene02 whose DNA cannot be placed."""
        recs = [r.split("\n", 1) for r in (EXAMPLES / "dna.fasta").read_text(encoding="utf-8").split(">")[1:]]
        d = {h.split()[0]: s.replace("\n", "") for h, s in recs}
        ids = list(d)
        g = d[ids[0]]
        d[ids[0]] = g[:27] + ("TGG" if g[27:30] != "TGG" else "GCC") + g[30:]
        d[ids[1]] = "".join(random.Random(1).choice("ACGT") for _ in range(90))
        dna = td / "dna.fasta"
        dna.write_text("".join(f">{k}\n{v}\n" for k, v in d.items()), encoding="ascii")
        return EXAMPLES / "protein.aln.fasta", dna

    def test_no_kept_when_nothing_was_written(self):
        engines = ["python"] + (["cpp"] if CPP else [])
        with tempfile.TemporaryDirectory() as t:
            pep, dna = self.inputs(Path(t))
            for eng in engines:
                r = coal(pep, dna, "-o", Path(t) / "out.fasta", env={"COAL_ENGINE": eng})
                self.assertEqual(r.returncode, 1, (eng, r.stderr))
                self.assertIn("coal: summary: 1 of 8 sequences has codons that do not match their amino acid.\n",
                              r.stderr, eng)
                self.assertIn("no output was written", r.stderr)
                r = coal(pep, dna, "-o", Path(t) / "out.fasta", "--keep-going", env={"COAL_ENGINE": eng})
                self.assertEqual(r.returncode, 2, (eng, r.stderr))
                self.assertIn("do not match their amino acid (kept in the alignment).", r.stderr, eng)


if __name__ == "__main__":
    unittest.main()
