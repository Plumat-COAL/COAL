"""Pairwise protein alignment for the codon step: cpp/src/align.{hpp,cpp}, coal-cpp-align and src/coal/align.py
(docs/alignment.md).

  * Optimal scores equal Biopython's PairwiseAligner and EMBOSS needle under 7 scoring settings
    (tests/fixtures/align/reference_scores.json, made by tests/tools/make_align_fixture.py; the live comparison
    with Biopython runs too when it is installed). The alignments returned score exactly that (independent
    rescoring) and hold every residue of both sequences (full length, nothing cut off).
  * The pure-Python kernel and coal-cpp-align give byte-identical alignments; -j 1, 2 and 0 give identical output.
  * Edge cases: empty and one-residue sequences, lower case, '*', B/Z/X/J/O/U; bad residues, bad options, pairs
    over --max-cells (that pair fails, exit 2, the others are written); --pairs and pairing by order.
  * The DIAMOND route: expand_local_hit in C++ (--diamond) and Python agree; the hit's block is kept unchanged;
    flanks gapped or aligned; inconsistent hits fail with a clear message; DIAMOND's masking ('X') is accepted.
    With a real `diamond` on PATH (or COAL_TEST_DIAMOND), its qseq_gapped/sseq_gapped output is expanded.
The C++ checks skip if build/coal-cpp-align (next to COAL_CPP_BIN) is not built.
"""
from __future__ import annotations

import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from test_cpp_cli import ROOT, cpp_binary  # noqa: E402

from coal import align as A  # noqa: E402

FIXTURE = ROOT / "tests" / "fixtures" / "align" / "reference_scores.json"
AA = "ARNDCQEGHILKMFPSTWYV"
TIMEOUT = 120


def align_binary() -> Path | None:
    env = os.environ.get("COAL_CPP_ALIGN")
    if env:
        return Path(env) if Path(env).is_file() else None
    cpp = cpp_binary()
    base = cpp.parent if cpp else ROOT / "build"
    for cand in (base / "coal-cpp-align", base / "coal-cpp-align.exe", base / "Release" / "coal-cpp-align.exe"):
        if cand.is_file() and (os.name != "nt" or cand.suffix.lower() == ".exe"):
            return cand
    return None


def no_binary(test):
    """Skip a test that needs coal-cpp-align, or fail it when COAL_REQUIRE_ALIGN is set (CI)."""
    if os.environ.get("COAL_REQUIRE_ALIGN"):
        test.fail("coal-cpp-align not built (COAL_REQUIRE_ALIGN is set)")
    test.skipTest("coal-cpp-align not built")


def mutate(s, rng, sub=0.3, indel=0.05):
    out = []
    for c in s:
        r = rng.random()
        if r < indel / 2:
            continue
        if r < indel:
            out.append("".join(rng.choice(AA) for _ in range(rng.randint(1, 6))))
        out.append(rng.choice(AA) if rng.random() < sub else c)
    return "".join(out)


def random_pairs(seed, n=80, maxlen=120):
    rng = random.Random(seed)
    pairs = []
    for k in range(n):
        a = "".join(rng.choice(AA) for _ in range(rng.randint(0 if k % 9 == 0 else 1, maxlen)))
        b = mutate(a, rng) if k % 3 else "".join(rng.choice(AA) for _ in range(rng.randint(0, maxlen)))
        if k % 4 == 0:
            b = "".join(rng.choice(AA) for _ in range(rng.randint(0, 20))) + b
        pairs.append((a, b))
    return pairs


def params_of(cfg) -> A.AlignParams:
    return A.AlignParams(cfg["matrix"], cfg["gap_open"], cfg["gap_extend"], cfg["free_end_gaps"])


def run_cpp(binary, args, **kw):
    return subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True, timeout=TIMEOUT, encoding="utf-8",
                          **kw)


def read_rows(path):
    rows = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith(">"):
                rows.append(line.rstrip("\n"))
    return [(rows[k], rows[k + 1]) for k in range(0, len(rows), 2)]


def write_pairs(path, pairs):
    with open(path, "w", newline="\n", encoding="utf-8") as fh:
        for i, (a, b) in enumerate(pairs):
            fh.write(f">a{i}\n{a}\n>b{i}\n{b}\n")


class FullLength:
    def check_rows(self, a, b, ra, rb):
        self.assertEqual(len(ra), len(rb))
        self.assertEqual(ra.replace("-", ""), a.upper(), "first sequence changed or cut")
        self.assertEqual(rb.replace("-", ""), b.upper(), "second sequence changed or cut")
        self.assertFalse(any(x == "-" and y == "-" for x, y in zip(ra, rb)), "an all-gap column")


class ReferenceScores(unittest.TestCase, FullLength):
    """Against the stored Biopython / EMBOSS needle optima (no tool needed)."""

    @classmethod
    def setUpClass(cls):
        with open(FIXTURE, encoding="utf-8") as fh:
            cls.data = json.load(fh)
        cls.binary = align_binary()

    def test_fixture_covers_the_defaults_and_presets(self):
        cfgs = {(c["matrix"], c["gap_open"], c["gap_extend"], c["free_end_gaps"]) for c in self.data["configs"]}
        p, w = A.AlignParams(), A.AlignParams.pwalign()
        self.assertIn((p.matrix, p.gap_open, p.gap_extend, p.free_end_gaps), cfgs)
        self.assertIn((w.matrix, w.gap_open, w.gap_extend, w.free_end_gaps), cfgs)
        self.assertIsNotNone(self.data["configs"][0]["needle"], "the default setting is checked against needle")

    def test_needle_agrees_with_biopython_except_one_residue_sequences(self):
        # EMBOSS needle reports 0 for a one-residue sequence with free end gaps (its end-gap rule treats the only
        # residue as an end): a documented needle quirk, docs/alignment.md. Everywhere else they agree.
        for cfg in self.data["configs"]:
            if cfg["needle"] is None:
                continue
            for (a, b), x, y in zip(self.data["pairs"], cfg["biopython"], cfg["needle"]):
                if y is None:
                    continue
                if min(len(a), len(b)) == 1 and cfg["free_end_gaps"]:
                    continue
                self.assertAlmostEqual(x, y, places=6, msg=(cfg["matrix"], a, b))

    def test_python_kernel(self):
        pairs = self.data["pairs"]
        for cfg in self.data["configs"]:
            p = params_of(cfg)
            for k in range(0, len(pairs), 3):  # pure Python is slow: every third pair
                a, b = pairs[k]
                r = A.align_pair(a, b, p)
                self.assertEqual(r.score, cfg["biopython"][k], (cfg, k))
                self.check_rows(a, b, r.a, r.b)
                self.assertEqual(A.score_alignment(r.a, r.b, p), r.score)

    def test_cpp(self):
        if self.binary is None:
            no_binary(self)
        pairs = self.data["pairs"]
        with tempfile.TemporaryDirectory() as tmp:
            inp, out, sc = (os.path.join(tmp, n) for n in ("in.fa", "out.fa", "s.tsv"))
            write_pairs(inp, pairs)
            for cfg in self.data["configs"]:
                p = params_of(cfg)
                r = run_cpp(self.binary, [inp, "-o", out, "--scores", sc, *p.args()])
                self.assertEqual(r.returncode, 0, r.stderr)
                rows = read_rows(out)
                with open(sc, encoding="utf-8") as fh:
                    scores = [float(line.split("\t")[6]) for line in fh.read().splitlines()[1:]]
                for k, ((a, b), (ra, rb)) in enumerate(zip(pairs, rows)):
                    self.assertEqual(scores[k], cfg["biopython"][k], (cfg, k))
                    self.check_rows(a, b, ra, rb)
                    self.assertEqual(A.score_alignment(ra, rb, p), scores[k])


class LiveBiopython(unittest.TestCase, FullLength):
    """New random pairs against Biopython itself, when it is installed."""

    def test_random_pairs(self):
        try:
            from Bio.Align import PairwiseAligner, substitution_matrices
        except ImportError:
            self.skipTest("Biopython not installed")
        binary = align_binary()
        if binary is None:
            no_binary(self)
        pairs = [(a, b) for a, b in random_pairs(7, n=120, maxlen=300) if a and b]
        for mat, o, e, free in (("BLOSUM62", 10, 0.5, True), ("BLOSUM62", 14, 4, False), ("PAM250", 3, 0, True),
                                ("BLOSUM80", 5, 5, False)):
            p = A.AlignParams(mat, o, e, free)
            al = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load(mat),
                                 open_gap_score=-o, extend_gap_score=-e)
            if free:
                al.end_gap_score = 0
            got = A.align_pairs(pairs, p, jobs=2, binary=str(binary))
            for (a, b), r in zip(pairs, got):
                self.assertIsNone(r.error)
                self.assertEqual(r.score, al.score(a, b), (mat, o, e, free, a, b))
                self.check_rows(a, b, r.a, r.b)


class CppTool(unittest.TestCase, FullLength):
    @classmethod
    def setUpClass(cls):
        cls.binary = align_binary()

    def setUp(self):
        if self.binary is None:
            no_binary(self)
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def path(self, name):
        return os.path.join(self.tmp, name)

    def test_python_and_cpp_identical(self):
        pairs = random_pairs(11, n=60, maxlen=90)
        for p in (A.AlignParams(), A.AlignParams.pwalign(), A.AlignParams("BLOSUM45", 7.25, 0.75, True)):
            cpp = A.align_pairs(pairs, p, binary=str(self.binary))
            for (a, b), c in zip(pairs, cpp):
                py = A.align_pair(a, b, p)
                self.assertEqual((c.a, c.b, c.score, c.identities, c.similarities, c.gaps),
                                 (py.a, py.b, py.score, py.identities, py.similarities, py.gaps))

    def test_threads_do_not_change_output(self):
        write_pairs(self.path("in.fa"), random_pairs(3, n=150, maxlen=200))
        outs = []
        for j in (1, 2, 0):
            r = run_cpp(self.binary, [self.path("in.fa"), "-j", j, "--scores", self.path("s.tsv")])
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(self.path("s.tsv"), encoding="utf-8") as fh:
                outs.append(r.stdout + fh.read())
        self.assertEqual(outs[0], outs[1])
        self.assertEqual(outs[0], outs[2])

    def test_stdout_and_files_have_the_same_bytes_on_every_system(self):
        # Windows: standard output in text mode would write \r\n where the -o file has \n
        write_pairs(self.path("in.fa"), random_pairs(5, n=20, maxlen=60))
        r = subprocess.run([str(self.binary), self.path("in.fa"), "--scores", self.path("s.tsv")], capture_output=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        r2 = subprocess.run([str(self.binary), self.path("in.fa"), "-o", self.path("o.fa")], capture_output=True)
        self.assertEqual(r2.returncode, 0, r2.stderr)
        with open(self.path("o.fa"), "rb") as fh:
            self.assertEqual(r.stdout, fh.read())
        self.assertGreater(r.stdout.count(b"\n"), 39)
        with open(self.path("s.tsv"), "rb") as fh:
            self.assertNotIn(b"\r", r.stdout + fh.read())

    def test_edge_residues(self):
        pairs = [("", ""), ("", "MKV"), ("M", ""), ("m", "M"), ("mkvJOU*", "MKVXXX*"), ("BZX", "DEA"),
                 ("W" * 50, "W"), ("A", "C" * 40)]
        write_pairs(self.path("in.fa"), pairs)
        for p in (A.AlignParams(), A.AlignParams.pwalign()):
            r = run_cpp(self.binary, [self.path("in.fa"), "--scores", self.path("s.tsv"), *p.args()])
            self.assertEqual(r.returncode, 0, r.stderr)
            rows = [x for x in r.stdout.splitlines() if not x.startswith(">")]
            with open(self.path("s.tsv"), encoding="utf-8") as fh:
                scores = [float(x.split("\t")[6]) for x in fh.read().splitlines()[1:]]
            for k, (a, b) in enumerate(pairs):
                ra, rb = rows[2 * k], rows[2 * k + 1]
                self.check_rows(a, b, ra, rb)
                py = A.align_pair(a, b, p)
                self.assertEqual((ra, rb, scores[k]), (py.a, py.b, py.score))
        # J, O, U score as X; '*' against '*' is +1 in BLOSUM62
        self.assertEqual(A.score_alignment("JOU*", "XXX*"), A.score_alignment("XXX*", "XXX*"))
        self.assertEqual(A.score_alignment("*", "*", A.AlignParams.pwalign()), 1)
        # empty against non-empty: free end gaps cost nothing; penalised, one gap of length 3
        self.assertEqual(A.align_pair("", "MKV").score, 0)
        self.assertEqual(A.align_pair("", "MKV", A.AlignParams.pwalign()).score, -(14 + 2 * 4))

    def test_bad_residue_fails_that_pair_only(self):
        write_pairs(self.path("in.fa"), [("MKV", "MKV"), ("MK-V", "MKV"), ("MK1", "MKV"), ("MKV", "MKVL")])
        r = run_cpp(self.binary, [self.path("in.fa"), "--scores", self.path("s.tsv")])
        self.assertEqual(r.returncode, 2)
        self.assertIn("pair 2 (a1, b1) skipped: first sequence: '-' at position 3 is not an amino acid letter",
                      r.stderr)
        self.assertIn("pair 3 (a2, b2) skipped", r.stderr)
        self.assertIn("2 of 4 pairs failed", r.stderr)
        self.assertEqual(r.stdout.count(">"), 4)
        with open(self.path("s.tsv"), encoding="utf-8") as fh:
            status = [x.split("\t")[10] for x in fh.read().splitlines()[1:]]
        self.assertEqual(status[0], "ok")
        self.assertTrue(status[1].startswith("first sequence"))
        self.assertEqual(status[3], "ok")

    def test_max_cells(self):
        write_pairs(self.path("in.fa"), [("M" * 100, "M" * 100), ("MK", "MK")])
        r = run_cpp(self.binary, [self.path("in.fa"), "--max-cells", "5000"])
        self.assertEqual(r.returncode, 2)
        self.assertIn("pair too large to align: 100 x 100 residues is more than max_cells = 5000", r.stderr)
        self.assertEqual(r.stdout, ">a1\nMK\n>b1\nMK\n")

    def test_options(self):
        write_pairs(self.path("in.fa"), [("MKV", "MKV")])
        for args, msg in (([ "--matrix", "BLOSUM99"], "unknown substitution matrix 'BLOSUM99'"),
                          (["--gap-open", "1", "--gap-extend", "2"], "gap_open >= gap_extend"),
                          (["--gap-open", "x"], "not a number"),
                          (["--end-gaps", "maybe"], "free or penalized"),
                          (["--bogus"], "unknown option --bogus"),
                          (["--gap-open", "10.0001"], "multiples of 0.001")):
            r = run_cpp(self.binary, [self.path("in.fa"), *args])
            self.assertEqual(r.returncode, 1, args)
            self.assertIn(msg, r.stderr, args)
        r = run_cpp(self.binary, ["--help"])
        self.assertEqual(r.returncode, 0)
        self.assertIn("--preset needle|pwalign", r.stdout)
        # --preset then an override
        r1 = run_cpp(self.binary, [self.path("in.fa"), "--preset", "pwalign", "--scores", self.path("a.tsv")])
        r2 = run_cpp(self.binary, [self.path("in.fa"), "--gap-open", "14", "--gap-extend", "4", "--end-gaps",
                                   "penalized", "--scores", self.path("b.tsv")])
        self.assertEqual(r1.returncode, 0)
        with open(self.path("a.tsv"), encoding="utf-8") as f1, open(self.path("b.tsv"), encoding="utf-8") as f2:
            self.assertEqual(f1.read(), f2.read())

    def test_pairing_modes(self):
        with open(self.path("A.fa"), "w", encoding="utf-8") as fh:
            fh.write(">x1 desc\nMKVL\nAAG\n>x2\nWWWW\n")
        with open(self.path("B.fa"), "w", encoding="utf-8") as fh:
            fh.write(">y1\r\nMKVLAG\r\n>y2\r\nWW\r\n")
        r = run_cpp(self.binary, [self.path("A.fa"), self.path("B.fa")])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(r.stdout.splitlines()[0::2], [">x1", ">y1", ">x2", ">y2"])
        with open(self.path("p.tsv"), "w", encoding="utf-8") as fh:
            fh.write("# a comment\nx2\ty1\nx1\ty2\n")
        r = run_cpp(self.binary, [self.path("A.fa"), self.path("B.fa"), "--pairs", self.path("p.tsv")])
        self.assertEqual(r.stdout.splitlines()[0::2], [">x2", ">y1", ">x1", ">y2"])
        r = run_cpp(self.binary, [self.path("A.fa"), "--pairs", self.path("p.tsv")])
        self.assertEqual(r.returncode, 1)
        self.assertIn("ID 'y1' is not in", r.stderr)
        with open(self.path("C.fa"), "w", encoding="utf-8") as fh:
            fh.write(">c1\nMK\n")
        r = run_cpp(self.binary, [self.path("A.fa"), self.path("C.fa")])
        self.assertEqual(r.returncode, 1)
        self.assertIn("pairing by order needs the same number", r.stderr)
        r = run_cpp(self.binary, [self.path("C.fa")])
        self.assertIn("odd number of records", r.stderr)


def local_hit_from(a, b, rng):
    """A DIAMOND-like hit: a block of the global alignment of a and b (residue coordinates 1-based)."""
    g = A.align_pair(a, b, A.AlignParams(gap_open=11, gap_extend=1))
    cols = [k for k, (x, y) in enumerate(zip(g.a, g.b)) if x != "-" and y != "-"]
    c0 = rng.choice(cols[: max(1, len(cols) // 3)])
    c1 = rng.choice(cols[-max(1, len(cols) // 3):])
    if c1 < c0:
        c0, c1 = c1, c0
    qs = len(g.a[:c0].replace("-", "")) + 1
    ss = len(g.b[:c0].replace("-", "")) + 1
    qg, sg = g.a[c0:c1 + 1], g.b[c0:c1 + 1]
    return A.LocalHit(qs, qs + len(qg.replace("-", "")) - 1, ss, ss + len(sg.replace("-", "")) - 1, qg, sg)


class LocalHits(unittest.TestCase, FullLength):
    def test_block_kept_and_full_length(self):
        rng = random.Random(5)
        for k in range(60):
            a = "".join(rng.choice(AA) for _ in range(rng.randint(30, 150)))
            b = "".join(rng.choice(AA) for _ in range(rng.randint(0, 30))) + mutate(a, rng) + "MKV"
            h = local_hit_from(a, b, rng)
            for fl in ("gap", "align"):
                r = A.expand_local_hit(a, b, h, flanks=fl)
                self.check_rows(a, b, r.a, r.b)
                core_a, core_b = h.qseq_gapped, h.sseq_gapped
                self.assertIn(core_a, r.a)
                i = r.a.index(core_a) if fl == "gap" else None
                if fl == "gap":
                    self.assertEqual(r.b[i:i + len(core_b)], core_b)
                    # staggered flanks: no residue pair outside the block
                    left = len(a[:h.qstart - 1]) + len(b[:h.sstart - 1])
                    self.assertEqual(i, left)
                self.assertEqual(r.score, A.score_alignment(r.a, r.b))

    def test_errors_and_masking(self):
        a, b = "MKVLAAGIW", "MKVLSAGIW"
        ok = A.LocalHit(2, 8, 2, 8, "KVLAAGI", "KVLSAGI")
        r = A.expand_local_hit(a, b, ok)
        # staggered flanks: the unaligned M and W of each sequence face gaps
        self.assertEqual((r.a, r.b), ("M-KVLAAGIW-", "-MKVLSAGI-W"))
        r = A.expand_local_hit(a, b, ok, flanks="align")
        self.assertEqual((r.a, r.b), ("MKVLAAGIW", "MKVLSAGIW"))
        masked = A.LocalHit(2, 8, 2, 8, "KXXXXGI", "KVLSAGI")
        self.assertEqual(A.expand_local_hit(a, b, masked).a, "M-KVLAAGIW-")
        for hit, msg in ((A.LocalHit(2, 8, 2, 8, "KVLWAGI", "KVLSAGI"), "qseq_gapped has 'W' where the query has 'A'"),
                         (A.LocalHit(0, 8, 2, 8, "MKVLAAGI", "KVLSAGI"), "query coordinates 0..8"),
                         (A.LocalHit(2, 8, 2, 12, "KVLAAGI", "KVLSAGI"), "subject coordinates 2..12"),
                         (A.LocalHit(2, 8, 2, 8, "KVLAAGI", "KVLSAG"), "differ in length"),
                         (A.LocalHit(2, 9, 2, 8, "KVLAAGI", "KVLSAGI"), "qseq_gapped has 7 residues")):
            with self.assertRaises(A.AlignError) as cm:
                A.expand_local_hit(a, b, hit)
            self.assertIn(msg, str(cm.exception))
            self.assertEqual(cm.exception.code, "local_hit")

    def test_cpp_matches_python(self):
        binary = align_binary()
        if binary is None:
            no_binary(self)
        rng = random.Random(9)
        with tempfile.TemporaryDirectory() as tmp:
            qa, sb, tsv = (os.path.join(tmp, n) for n in ("q.fa", "s.fa", "hits.tsv"))
            hits, seqs = [], []
            with open(qa, "w", encoding="utf-8") as fq, open(sb, "w", encoding="utf-8") as fs, open(tsv, "w", encoding="utf-8") as ft:
                for k in range(40):
                    a = "".join(rng.choice(AA) for _ in range(rng.randint(30, 150)))
                    b = "".join(rng.choice(AA) for _ in range(rng.randint(0, 30))) + mutate(a, rng)
                    h = local_hit_from(a, b, rng)
                    fq.write(f">q{k}\n{a}\n")
                    fs.write(f">s{k}\n{b}\n")
                    ft.write(f"q{k}\ts{k}\t{h.qstart}\t{h.qend}\t{h.sstart}\t{h.send}\t{h.qseq_gapped}\t"
                             f"{h.sseq_gapped}\t99.9\n")
                    if k % 5 == 0:  # a second HSP of the same pair: ignored
                        ft.write(f"q{k}\ts{k}\t1\t1\t1\t1\tX\tX\t1.0\n")
                    hits.append(h)
                    seqs.append((a, b))
                ft.write("q1\ts2\t1\t200\t1\t3\tMKV\tMKV\t1\n")  # another pair: bad coordinates
                ft.write("q2\ts1\t1\t1\t1\t1\tX\tX\t1\n")
            for fl in ("gap", "align"):
                r = run_cpp(binary, [qa, sb, "--diamond", tsv, "--flanks", fl, "--scores", os.path.join(tmp, "s")])
                self.assertEqual(r.returncode, 2, r.stderr)
                self.assertIn("local hit: query coordinates 1..200 do not fit a query of", r.stderr)
                rows = [x for x in r.stdout.splitlines() if not x.startswith(">")]
                for k, ((a, b), h) in enumerate(zip(seqs, hits)):
                    py = A.expand_local_hit(a, b, h, flanks=fl)
                    self.assertEqual((rows[2 * k], rows[2 * k + 1]), (py.a, py.b))
                parsed = list(A.read_diamond_hits(tsv))
                self.assertEqual(len(parsed), 42)
                self.assertEqual(parsed[0].extra, ["99.9"])

    def test_real_diamond(self):
        diamond = os.environ.get("COAL_TEST_DIAMOND") or shutil.which("diamond")
        binary = align_binary()
        if not diamond or binary is None:
            if os.environ.get("COAL_REQUIRE_DIAMOND"):   # CI's C++ job, after it installed DIAMOND
                self.fail("diamond or coal-cpp-align not available (COAL_REQUIRE_DIAMOND is set)")
            self.skipTest("diamond or coal-cpp-align not available")
        rng = random.Random(1)
        with tempfile.TemporaryDirectory() as tmp:
            q, s, db, out = (os.path.join(tmp, n) for n in ("q.fa", "s.fa", "db", "hits.tsv"))
            qs, ss = {}, {}
            for k in range(30):
                a = "".join(rng.choice(AA) for _ in range(rng.randint(120, 300)))
                qs[f"q{k}"] = "M" + a
                ss[f"s{k}"] = "MSTE" + mutate(a, rng, sub=0.15, indel=0.02) + "GGKK"
            for path, recs in ((q, qs), (s, ss)):
                with open(path, "w", encoding="utf-8") as fh:
                    fh.writelines(f">{i}\n{x}\n" for i, x in recs.items())
            subprocess.run([diamond, "makedb", "--in", s, "-d", db, "--quiet"], check=True, timeout=TIMEOUT)
            subprocess.run([diamond, "blastp", "-q", q, "-d", db, "-o", out, "--outfmt", "6", *A.DIAMOND_COLUMNS,
                            "--max-target-seqs", "1", "--quiet", "--threads", "1"], check=True, timeout=TIMEOUT)
            hits = list(A.read_diamond_hits(out))
            self.assertGreater(len(hits), 20)
            for h in hits:
                r = A.expand_local_hit(qs[h.qseqid], ss[h.sseqid], h)
                self.check_rows(qs[h.qseqid], ss[h.sseqid], r.a, r.b)
            r = run_cpp(binary, [q, s, "--diamond", out])
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(r.stdout.count(">"), 2 * len(hits))


class ParamsAndHelpers(unittest.TestCase):
    def test_presets(self):
        self.assertEqual(A.AlignParams.needle(), A.AlignParams("BLOSUM62", 10, 0.5, True))
        self.assertEqual(A.AlignParams.pwalign(), A.AlignParams("BLOSUM62", 14, 4, False))
        with self.assertRaises(A.AlignError):
            A.AlignParams("BLOSUM99").check()
        with self.assertRaises(A.AlignError):
            A.align_pair("MK-V", "MKV")

    def test_matrices_match_cpp_source(self):
        src = (ROOT / "cpp" / "src" / "align.cpp").read_text(encoding="utf-8")
        for name, mat in A.MATRICES.items():
            block = src[src.index(f'{{"{name}", {{'):]
            nums = [int(x) for x in block[block.index("{", 3):block.index("}},")].replace("{", " ").replace(
                "}", " ").replace(",", " ").split()]
            self.assertEqual(nums, [v for row in mat for v in row], name)

    def test_python_fallback_without_binary(self):
        old = os.environ.get("COAL_CPP_ALIGN")
        os.environ["COAL_CPP_ALIGN"] = os.path.join(tempfile.gettempdir(), "no-such-coal-cpp-align")
        try:
            self.assertIsNone(A.find_align_binary())
            r = A.align_pairs([("id1", "MKV", "id2", "MKIV")])
            self.assertEqual((r[0].id_a, r[0].a, r[0].b), ("id1", "MKV-", "MKIV"))  # end gaps free
            with self.assertRaises(RuntimeError):
                A.align_pairs([("MKV", "MKV")], python_fallback=False)
        finally:
            if old is None:
                del os.environ["COAL_CPP_ALIGN"]
            else:
                os.environ["COAL_CPP_ALIGN"] = old


if __name__ == "__main__":
    unittest.main()
