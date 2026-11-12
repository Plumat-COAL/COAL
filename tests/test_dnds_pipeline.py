"""`coal pipeline` (src/coal/dnds_pipeline.py) and its native steps `coal-cpp translate` / `coal-cpp codon-dnds`
(cpp/src/dnds.cpp, dnds_cli.cpp); docs/pipeline.md.

Inputs: tests/fixtures/pipeline (synthetic, made by tests/tools/make_pipeline_fixture.py; no real sequence data):
ordinary ortholog pairs plus an internal stop, a CDS that is not a multiple of 3, ambiguous bases, identical
sequences, RNA / lower case, genes without a hit, a bad character, an empty record, a repeated ID; and four gene
families (6, 12, 2 and 1 members) with a clusters file in DIAMOND's two-column format.

  * Units (nothing needed): the pair sample, the clusters reader, the quartiles, the DIAMOND installer with a faked
    download (nothing is fetched).
  * With build/coal-cpp (COAL_CPP_BIN): translate; the codon-dnds stream (pairs, rows, hits, failed members, bad
    input); pair mode with a given pairs file and cluster mode with a given clusters file, against the same steps
    done one by one (coal.align -> coal-cpp -> coal-cpp dnds) and against statistics recomputed here; failures,
    exit codes, the report (docs/pipeline.schema.json), --resume, the same bytes for every -j.
  * With DIAMOND too (PATH, COAL_DIAMOND or COAL's cache): reciprocal best hits, best hits, the DIAMOND-alignment
    route, DeepClust. With Clustal Omega / MAFFT too: cluster MSAs against `coal-cpp dnds` on the kept codon MSA.
Each part is skipped when its tool is missing. No network.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
import os
import stat
import statistics
import subprocess
import sys
import tarfile
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import align as A  # noqa: E402
from coal import dnds_pipeline as P  # noqa: E402
from coal import msa  # noqa: E402
from test_report_schema import validate  # noqa: E402

FIX = ROOT / "tests" / "fixtures" / "pipeline"
Q, S = FIX / "query_cds.fasta", FIX / "subject_cds.fasta"
FAM, FAM_CLUSTERS = FIX / "family_cds.fasta", FIX / "family_clusters.tsv"
SCHEMA = json.loads((ROOT / "docs" / "pipeline.schema.json").read_text())
PLAIN = [(f"g{k:02d}", f"h{k:02d}") for k in range(1, 13)]
EDGE = [("stop_q", "stop_s"), ("frame_q", "frame_s"), ("amb_q", "amb_s"), ("same_q", "same_s"), ("rna_q", "rna_s")]


def cpp_bin() -> Path | None:
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / (p.name + ".exe")):
        if cand.is_file() and (os.name != "nt" or cand.suffix.lower() == ".exe"):
            return cand
    return None


CPP = cpp_bin()
DIAMOND = P.find_diamond()
CLUSTALO = msa.find_clustalo()
MAFFT = msa.find_mafft()
needs_cpp = unittest.skipUnless(CPP, "build/coal-cpp (COAL_CPP_BIN) is not built")
needs_diamond = unittest.skipUnless(CPP and DIAMOND, "needs coal-cpp and DIAMOND")
needs_clustalo = unittest.skipUnless(CPP and CLUSTALO, "needs coal-cpp and Clustal Omega")


def read_fasta(path) -> dict[str, str]:
    """{id: sequence}; of a repeated ID the first record."""
    recs, cur = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = (line[1:].split() or [""])[0]
            cur = None if name in recs else recs.setdefault(name, [])
        elif cur is not None:
            cur.append(line.strip())
    return {k: "".join(v) for k, v in recs.items()}


def read_tsv(path) -> list[dict]:
    lines = Path(path).read_text().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]


def run_cpp(*args, stdin: bytes | None = None):
    return subprocess.run([str(CPP), *map(str, args)], input=stdin, capture_output=True)


def stream(groups) -> bytes:
    """The codon-dnds input for [(gid, mode, lines)], lines given as tuples of fields."""
    out = []
    for gid, mode, lines in groups:
        out.append(f"G\t{gid}\t{mode}\n")
        out += ["\t".join(map(str, ln)) + "\n" for ln in lines]
        out.append("E\n")
    return "".join(out).encode()


def by_hand(pairs, q: dict, s: dict, methods=("YN",)) -> dict:
    """The steps one by one, for a few pairs: translate in Python (coal.codon-free: the standard code), align with
    coal.align (coal-cpp-align, or the byte-identical pure-Python kernel), thread the codons with `coal-cpp PEP NUC`,
    then `coal-cpp dnds --pairs consecutive`. -> {(q, s, method): row dict}."""
    bases = "TCAG"
    aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    code = {a + b + c: aas[16 * i + 4 * j + k] for i, a in enumerate(bases) for j, b in enumerate(bases)
            for k, c in enumerate(bases)}

    def tr(cds: str) -> str:
        cds = cds.upper().replace("U", "T")
        cds = cds[:len(cds) - len(cds) % 3]
        p = "".join(code.get(cds[i:i + 3], "X") for i in range(0, len(cds), 3))
        return p[:-1] if p.endswith("*") else p
    alns = A.align_pairs([(tr(q[x]), tr(s[y])) for x, y in pairs])
    with tempfile.TemporaryDirectory() as td:
        stacked = os.path.join(td, "stacked.fasta")
        with open(stacked, "w") as fs:
            for (x, y), al in zip(pairs, alns):
                pep, nuc = os.path.join(td, "pep.fa"), os.path.join(td, "nuc.fa")
                Path(pep).write_text(f">a\n{al.a}\n>b\n{al.b}\n")
                Path(nuc).write_text(f">a\n{q[x]}\n>b\n{s[y]}\n")
                r = subprocess.run([str(CPP), pep, nuc, "-q"], capture_output=True, text=True, check=True)
                r1, r2 = ("".join(p.splitlines()[1:]) for p in r.stdout.split(">")[1:])
                fs.write(f">{x}\n{r1}\n>{y}\n{r2}\n")
        out = os.path.join(td, "o.tsv")
        subprocess.run([str(CPP), "dnds", stacked, "--pairs", "consecutive", "-o", out, "-q",
                        *[a for m in methods for a in ("-m", m)]], check=True)
        return {(r["seq1"], r["seq2"], r["method"]): r for r in read_tsv(out)}


VALUES = ("dN", "dS", "dNdS", "SE_dN", "SE_dS", "N_sites", "S_sites", "Nd", "Sd", "codons_compared")


class PipelineCase(unittest.TestCase):
    """Runs with COAL_CPP pointing at the built binary and a temporary folder."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = Path(self._td.name)
        env = {"COAL_CPP": str(CPP)} if CPP else {}
        self._env = mock.patch.dict(os.environ, env)
        self._env.start()

    def tearDown(self):
        self._env.stop()
        self._td.cleanup()

    def pairs_file(self, pairs) -> Path:
        p = self.tmp / "pairs.tsv"
        p.write_text("".join(f"{a}\t{b}\n" for a, b in pairs))
        return p

    def check_report(self, res: P.PipelineResult):
        validate(res.report, SCHEMA)
        on_disk = json.loads(Path(res.outputs["report"]).read_text())
        self.assertEqual(on_disk, json.loads(json.dumps(res.report)))
        self.assertEqual(res.report["exit_code"], res.exit_code)


# ---------------------------------------------------------------------------------------------- units
class Units(unittest.TestCase):
    def test_sample_pairs(self):
        for n in (2, 3, 7, 30):
            total = n * (n - 1) // 2
            every = [(i, j) for i in range(n) for j in range(i + 1, n)]
            self.assertEqual(P.sample_pairs(n, total, 1, "c"), every)
            self.assertEqual(P.sample_pairs(n, total + 5, 1, "c"), every)
            # the index -> pair mapping: sampling all but one leaves exactly one pair out
            if total > 1:
                got = P.sample_pairs(n, total - 1, 3, "c")
                self.assertEqual(len(set(got)), total - 1)
                self.assertTrue(set(got) < set(every))
                self.assertEqual(got, sorted(got))
        a = P.sample_pairs(300, 500, 7, "famF")
        self.assertEqual(a, P.sample_pairs(300, 500, 7, "famF"))
        self.assertNotEqual(a, P.sample_pairs(300, 500, 8, "famF"))
        self.assertNotEqual(a, P.sample_pairs(300, 500, 7, "other"))
        self.assertEqual(len(set(a)), 500)
        self.assertTrue(all(0 <= i < j < 300 for i, j in a))
        # stable across Python versions and platforms (random.Random seeded with a string)
        self.assertEqual(P.sample_pairs(10, 3, 1, "x"), P.sample_pairs(10, 3, "1", "x"))

    def test_read_clusters(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "c.tsv"
            p.write_text("centroid\tmember\n# a comment\nr1\tr1\nr1\ta\nr2\tb\n\nr1\tc\nr2\ta\nr3\tr3\r\n")
            clusters, dup = P.read_clusters(str(p))
            self.assertEqual(clusters, {"r1": ["r1", "a", "c"], "r2": ["r2", "b"], "r3": ["r3"]})
            self.assertEqual(dup, 1)
            p.write_text("only_one_column\n")
            with self.assertRaises(P.PipelineError):
                P.read_clusters(str(p))
        with self.assertRaises(P.PipelineError):
            P.read_clusters(os.path.join(td, "missing.tsv"))

    def test_quartiles_are_type_7(self):
        for v in ([1.0], [1.0, 2.0], [0.3, 0.1, 0.9, 0.4, 0.5], [float(x * x) for x in range(11)]):
            s = sorted(v)
            if len(s) > 1:
                q = statistics.quantiles(s, n=4, method="inclusive")
                self.assertAlmostEqual(P._quantile(s, 0.25), q[0])
                self.assertAlmostEqual(P._quantile(s, 0.75), q[2])
            self.assertAlmostEqual(P._quantile(s, 0.5), statistics.median(s))
        self.assertEqual(P._stats([]), ["NA"] * 6)
        self.assertEqual(P._stats([0.0, 0.5]), ["0.25", "0.25", "0.125", "0.375", "0", "0.5"])

    def test_fasta_index(self):
        idx = P.FastaIndex([str(Q), str(S)])
        seqs = read_fasta(Q)
        self.assertEqual(idx.get("g02").decode(), seqs["g02"])
        self.assertEqual(idx.get("g01").decode(), seqs["g01"])     # a repeated ID: the first record
        self.assertEqual(idx.get("empty_q"), b"")
        self.assertIn("h05", idx)
        self.assertNotIn("nope", idx)
        idx.close()

    def test_find_diamond_env(self):
        with tempfile.TemporaryDirectory() as td, mock.patch.dict(os.environ, {"COAL_DIAMOND": td + "/none"}):
            self.assertIsNone(P.find_diamond())

    @unittest.skipIf(os.name == "nt", "the stand-in tool is a POSIX shell script")
    def test_install_diamond_with_a_faked_download(self):
        with tempfile.TemporaryDirectory() as td, \
                mock.patch.dict(os.environ, {"COAL_CACHE_DIR": td}), \
                mock.patch.object(msa, "_platform", return_value=("Linux", "x86_64")):
            script = b"#!/bin/sh\necho 'diamond version 2.2.8'\n"

            def archive(body: bytes) -> bytes:
                buf = io.BytesIO()
                with tarfile.open(fileobj=buf, mode="w:gz") as tf:
                    info = tarfile.TarInfo("diamond")
                    info.size = len(body)
                    info.mode = 0o755
                    tf.addfile(info, io.BytesIO(body))
                return buf.getvalue()
            good = archive(script)
            asset = ("2.2.8", "diamond-linux64.tar.gz")

            def serving(data, calls):
                def fake(url, dest, quiet):
                    calls.append(url)
                    with open(dest + ".part", "wb") as fh:
                        fh.write(data)
                    return hashlib.sha256(data).hexdigest()
                return mock.patch.object(msa, "_download", side_effect=fake)
            calls = []
            # the real pinned hash: a file that is not the release asset is deleted, never unpacked or run
            self.assertRegex(P._DIAMOND_SHA256[asset], r"^[0-9a-f]{64}$")
            with serving(good, calls):
                with self.assertRaises(P.PipelineError) as cm:
                    P.install_diamond(method="official", quiet=True)
            self.assertIn("failed its SHA-256 check and was deleted", str(cm.exception))
            self.assertIn("/v2.2.8/diamond-linux64.tar.gz", calls[0])
            self.assertTrue(calls[0].startswith("https://github.com/bbuchfink/diamond/releases/download/"))
            self.assertEqual(os.listdir(os.path.join(td, "tools", "diamond")), [])
            # with the hash of the served file pinned, the install goes through
            calls.clear()
            with serving(good, calls) as dl, \
                    mock.patch.dict(P._DIAMOND_SHA256, {asset: hashlib.sha256(good).hexdigest()}):
                exe = P.install_diamond(method="official", quiet=True)
                self.assertTrue(os.access(exe, os.X_OK))
                self.assertEqual(P._diamond_version(exe), "2.2.8")
                self.assertEqual(P.install_diamond(method="official", quiet=True), exe)   # the cached copy
                self.assertEqual(dl.call_count, 1)
            with mock.patch.dict(os.environ, {"PATH": ""}):
                os.environ.pop("COAL_DIAMOND", None)
                self.assertEqual(P.find_diamond(), exe)
            # a version COAL has no hash for: refused before downloading; fine with the caller's hash
            calls.clear()
            with serving(good, calls):
                with self.assertRaises(P.PipelineError) as cm:
                    P.install_diamond(method="official", version="2.1.0", quiet=True, force=True)
                self.assertIn("refusing to download and run DIAMOND 2.1.0", str(cm.exception))
                self.assertEqual(calls, [])
                self.assertEqual(P.install_diamond(method="official", version="2.1.0", quiet=True, force=True,
                                                   sha256=hashlib.sha256(good).hexdigest()), exe)
            # a binary that does not run is refused
            bad = archive(b"#!/bin/sh\nexit 3\n")
            with serving(bad, calls), mock.patch.dict(P._DIAMOND_SHA256, {asset: hashlib.sha256(bad).hexdigest()}):
                with self.assertRaises(P.PipelineError):
                    P.install_diamond(method="official", quiet=True, force=True)

    def test_install_command_takes_a_hash_or_explicit_permission(self):
        seen = {}

        def fake_clustalo(**kw):
            seen.update(kw)
            return "/cache/clustalo"
        with mock.patch.object(msa, "install_clustalo", fake_clustalo), redirect_stdout(io.StringIO()):
            self.assertEqual(P.main(["install", "clustalo", "--sha256", "A" * 64, "-q"]), 0)
            self.assertEqual((seen["sha256"], seen["allow_unverified"]), ("A" * 64, False))
            self.assertEqual(P.main(["install", "clustalo", "--allow-unverified", "-q"]), 0)
            self.assertEqual((seen["sha256"], seen["allow_unverified"]), (None, True))
        for argv in (["install", "clustalo", "--sha256", "abc"],
                     ["install", "clustalo", "diamond", "--sha256", "a" * 64]):
            with self.assertRaises(SystemExit) as cm, redirect_stderr(io.StringIO()):
                P.main(argv)
            self.assertEqual(cm.exception.code, 2, argv)

    def test_command_line_errors(self):
        for argv in ([], ["pairs", "-o", "x"], ["clusters", "--cds", "a", "-o", "x", "--max-cluster-size", "0"],
                     ["pairs", "--query", "a", "--subject", "b", "-o", "x", "-m", "GY"]):
            with self.assertRaises(SystemExit) as cm, redirect_stderr(io.StringIO()):
                raise SystemExit(P.main(argv))
            self.assertEqual(cm.exception.code, 2, argv)

    def test_no_coal_cpp_is_exit_1(self):
        with tempfile.TemporaryDirectory() as td, mock.patch.dict(os.environ, {"COAL_CPP": td + "/none"}), \
                redirect_stderr(io.StringIO()) as err:
            res = P.run_pairs(Q, S, os.path.join(td, "out"))
        self.assertEqual(res.exit_code, 1)
        self.assertIn("coal-cpp", err.getvalue())


# ---------------------------------------------------------------------------------------------- native steps
@needs_cpp
class Translate(PipelineCase):
    def test_translate(self):
        out, genes = self.tmp / "p.faa", self.tmp / "g.tsv"
        r = run_cpp("translate", Q, "-o", out, "--genes", genes)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn(b"21 gene(s): 17 ok, 1 warning, 3 failed", r.stderr)
        prot = read_fasta(out)
        g = {}
        for row in read_tsv(genes):
            g.setdefault(row["id"], []).append(row)
        self.assertEqual(len(prot), 18)
        cds = read_fasta(Q)
        self.assertEqual(len(prot["g02"]), len(cds["g02"]) // 3 - 1)          # the terminal stop is left out
        self.assertTrue(prot["g02"].startswith("M") and "*" not in prot["g02"])
        self.assertEqual(g["g02"][0]["terminal_stop"], "1")
        self.assertEqual(prot["stop_q"].count("*"), 1)
        self.assertEqual((g["stop_q"][0]["status"], g["stop_q"][0]["code"], g["stop_q"][0]["internal_stops"]),
                         ("warning", "internal_stop", "1"))
        self.assertEqual(g["amb_q"][0]["status"], "ok")
        self.assertEqual(g["amb_q"][0]["ambiguous_codons"], "2")
        self.assertEqual(prot["amb_q"][10], "X")                              # NNN
        self.assertEqual(g["rna_q"][0]["status"], "ok")                       # lower case, U
        self.assertEqual([(x["status"], x["code"]) for x in g["g01"]], [("ok", ""), ("failed", "duplicate_id")])
        self.assertEqual((g["bad_q"][0]["status"], g["bad_q"][0]["code"]), ("failed", "cds_invalid_character"))
        self.assertEqual(g["empty_q"][0]["code"], "cds_empty")
        for bad in ("bad_q", "empty_q"):
            self.assertNotIn(bad, prot)
        r = run_cpp("translate", S, "--genes", self.tmp / "s.tsv")
        s = {row["id"]: row for row in read_tsv(self.tmp / "s.tsv")}
        self.assertEqual((s["frame_s"]["status"], s["frame_s"]["trailing_nt"], s["frame_s"]["code"]),
                         ("warning", "1", "cds_length"))
        self.assertEqual(s["short_s"]["code"], "cds_too_short")
        self.assertIn(b">h01\n", r.stdout)

    def test_translate_other_code_and_errors(self):
        f = self.tmp / "c.fa"
        f.write_text(">a\nATGTGAAAATAA\n")                                    # TGA: stop (1), Trp (4)
        self.assertIn(b"M*K", run_cpp("translate", f, "-q").stdout)
        self.assertIn(b"MWK", run_cpp("translate", f, "-q", "-c", "4").stdout)
        self.assertEqual(run_cpp("translate", self.tmp / "missing.fa").returncode, 1)
        f.write_text("ATGAAA\n")
        r = run_cpp("translate", f)
        self.assertEqual(r.returncode, 1)
        self.assertIn(b"before the first '>'", r.stderr)
        self.assertEqual(run_cpp("translate").returncode, 2)
        self.assertEqual(run_cpp("translate", f, "-c", "7").returncode, 2)
        # through `coal`
        r = subprocess.run([sys.executable, "-m", "coal", "translate", str(Q), "-q"], capture_output=True,
                           env={**os.environ, "PYTHONPATH": str(ROOT / "src")})
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn(b">g02\nM", r.stdout)


@needs_cpp
class CodonDnds(PipelineCase):
    def setUp(self):
        super().setUp()
        self.q, self.s = read_fasta(Q), read_fasta(S)

    def run_stream(self, data: bytes, *args):
        out = self.tmp / "o.tsv"
        r = run_cpp("codon-dnds", "-o", out, *args, stdin=data)
        return r, (read_tsv(out) if r.returncode == 0 else [])

    def test_pairs_equal_the_steps_one_by_one(self):
        pairs = PLAIN[1:4] + EDGE
        data = stream([(f"{x}&{y}", "all", [("M", x, self.q[x]), ("M", y, self.s[y])]) for x, y in pairs])
        r, rows = self.run_stream(data, "-m", "YN", "-m", "NG", "--codon-out", self.tmp / "c.fa")
        self.assertEqual(r.returncode, 0, r.stderr)
        hand = by_hand(pairs, self.q, self.s, ("YN", "NG"))
        self.assertEqual(len(rows), 2 * len(pairs))
        self.assertEqual([r_["method"] for r_ in rows[:2]], ["NG", "YN"])     # KaKs_Calculator's order
        for row in rows:
            h = hand[(row["seq1"], row["seq2"], row["method"])]
            self.assertEqual([row[c] for c in VALUES], [h[c] for c in VALUES], row["seq1"])
        st = {(r_["seq1"], r_["method"]): r_ for r_ in rows}
        self.assertEqual(st[("g02", "YN")]["status"], "ok")
        self.assertEqual(st[("stop_q", "YN")]["status"], "warning")
        self.assertIn("seq1: 1 internal stop codon", st[("stop_q", "YN")]["note"])
        self.assertEqual(st[("frame_q", "YN")]["status"], "warning")
        self.assertIn("seq2: CDS length 301 is not a multiple of 3", st[("frame_q", "YN")]["note"])
        self.assertEqual((st[("same_q", "YN")]["dN"], st[("same_q", "YN")]["dS"], st[("same_q", "YN")]["dNdS"],
                          st[("same_q", "YN")]["prot_identity"]), ("0", "0", "NA", "100.00"))
        self.assertEqual(st[("amb_q", "YN")]["status"], "ok")
        self.assertEqual(st[("rna_q", "YN")]["status"], "ok")
        kept = (self.tmp / "c.fa").read_text().splitlines()
        self.assertEqual(len(kept), 4 * len(pairs))
        k = kept.index(">rna_q")
        self.assertTrue(kept[k + 1].startswith("aug"))                        # the DNA is copied as written
        self.assertEqual(len(kept[k + 1]), len(kept[k + 3]))

    def test_failed_members_groups_and_pair_lists(self):
        q, s = self.q, self.s
        fam = read_fasta(FAM)
        three = [("M", n, fam[n]) for n in ("famA_001", "famA_002", "famA_003")]
        data = stream([
            ("p1", "all", [("M", "g02", q["g02"]), ("F", "h02", "no CDS record")]),
            ("p2", "all", [("M", "g03", q["g03"]), ("M", "bad", "ATGAAA!AA")]),
            ("p3", "all", [("X", "msa_failed", "clustalo failed: boom"), ("M", "g04", q["g04"]),
                           ("M", "h04", s["h04"])]),
            ("c1", "all", three),                                             # three members, no rows
            ("c2", "rep", [("M", "g05", q["g05"]), ("M", "h05", s["h05"])]),
            ("c3", "list", [("M", "g06", q["g06"]), ("M", "h06", s["h06"]), ("P", 1, 0), ("P", 0, 1), ("P", 0, 5)]),
        ])
        r, rows = self.run_stream(data, "--group-column", "cluster_id", "--failures", self.tmp / "f.tsv",
                                  "--kaks-out", self.tmp / "k.kaks", "--codon-out", self.tmp / "c.fa")
        self.assertEqual(r.returncode, 0, r.stderr)
        got = [(x["cluster_id"], x["seq1"], x["seq2"], x["status"]) for x in rows]
        self.assertEqual(got, [("p1", "g02", "h02", "error"), ("p2", "g03", "bad", "error"),
                               ("p3", "g04", "h04", "error"),
                               ("c1", "famA_001", "famA_002", "error"), ("c1", "famA_001", "famA_003", "error"),
                               ("c1", "famA_002", "famA_003", "error"),
                               ("c2", "g05", "h05", "ok"), ("c3", "h06", "g06", "ok"), ("c3", "g06", "h06", "ok")])
        self.assertEqual(rows[0]["note"], "seq2: no CDS record")
        self.assertIn("'!' at position 7", rows[1]["note"])
        self.assertEqual(rows[2]["note"], "clustalo failed: boom")
        self.assertIn("needs a protein alignment", rows[3]["note"])
        self.assertEqual(rows[0]["dN"], "NA")
        self.assertEqual(rows[7]["dN"], rows[8]["dN"])                        # the same pair, both orders
        fails = read_tsv(self.tmp / "f.tsv")
        self.assertEqual([(f["group"], f["seq"], f["code"]) for f in fails],
                         [("p2", "bad", "cds_invalid_character"), ("c1", "", "alignment_missing")])
        kaks = (self.tmp / "k.kaks").read_text().splitlines()
        self.assertTrue(kaks[0].startswith("Sequence\tMethod\tKa\tKs\tKa/Ks\tP-Value(Fisher)"))
        self.assertEqual([ln.split("\t")[0] for ln in kaks[1:]], ["g05&h05", "h06&g06", "g06&h06"])
        self.assertEqual((self.tmp / "c.fa").read_text().count(">"), 4)       # only the two whole pairs

    def test_rows_and_hits(self):
        # rows given (an MSA): equal to `coal-cpp dnds` on the codon MSA written for the same group
        fam = read_fasta(FAM)
        names = ["famC_001", "famC_002"]
        rows_in = []
        r0 = run_cpp("translate", FAM, "-q")
        prot = {}
        cur = None
        for line in r0.stdout.decode().splitlines():
            if line.startswith(">"):
                cur = line[1:]
            else:
                prot[cur] = line
        al = A.align_pairs([(prot[names[0]], prot[names[1]])])[0]
        for n, row in zip(names, (al.a, al.b)):
            rows_in.append(("M", n, fam[n], row))
        d = self.tmp / "codon"
        d.mkdir()
        r, rows = self.run_stream(stream([("fam C/1", "all", rows_in)]), "--codon-dir", d)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual([p.name for p in d.iterdir()], ["fam_C_1.fasta"])   # a safe file name
        ref = run_cpp("dnds", d / "fam_C_1.fasta", "-q").stdout.decode().splitlines()[1].split("\t")
        self.assertEqual([rows[0][c] for c in VALUES], ref[3:13])
        # no rows: the same alignment is made here
        r, rows2 = self.run_stream(stream([("c", "all", [x[:3] for x in rows_in])]))
        self.assertEqual([rows2[0][c] for c in VALUES], ref[3:13])
        # ragged rows
        bad = [rows_in[0], rows_in[1][:3] + (rows_in[1][3] + "-",)]
        r, rows3 = self.run_stream(stream([("c", "all", bad)]))
        self.assertEqual(rows3[0]["status"], "error")
        self.assertIn("differ in length", rows3[0]["note"])
        # a local hit covering the whole of both proteins = the full alignment
        a, b = prot[names[0]], prot[names[1]]
        hit = ("H", 1, len(a), 1, len(b), al.a, al.b)
        r, rows4 = self.run_stream(stream([("c", "all", [x[:3] for x in rows_in] + [hit])]))
        self.assertEqual([rows4[0][c] for c in VALUES], ref[3:13])
        # a hit that does not fit the sequences is that pair's error
        bad_hit = ("H", 1, 5, 1, 5, "WWWWW", "WWWWW")
        r, rows5 = self.run_stream(stream([("c", "all", [x[:3] for x in rows_in] + [bad_hit])]))
        self.assertEqual(rows5[0]["status"], "error")

    def test_threads_give_the_same_bytes_and_bad_input(self):
        pairs = PLAIN + EDGE
        data = stream([(f"{x}&{y}", "all", [("M", x, self.q[x]), ("M", y, self.s[y])]) for x, y in pairs] * 3)
        outs = []
        for j in (1, 2, 5):
            out = self.tmp / f"o{j}.tsv"
            self.assertEqual(run_cpp("codon-dnds", "-o", out, "-j", j, "-q", stdin=data).returncode, 0)
            outs.append(out.read_bytes())
        self.assertEqual(outs[0], outs[1])
        self.assertEqual(outs[0], outs[2])
        f = self.tmp / "in.txt"
        f.write_bytes(data)
        self.assertEqual(run_cpp("codon-dnds", "-i", f, "-o", self.tmp / "o.tsv", "-q").returncode, 0)
        self.assertEqual((self.tmp / "o.tsv").read_bytes(), outs[0])
        for broken, why in ((b"M\ta\tATG\n", b"before the first G"), (b"G\tg\tall\nM\ta\tATGAAA\n", b"no E line"),
                            (b"G\tg\tsome\nE\n", b"all, rep or list"), (b"G\tg\tall\nZ\tq\nE\n", b"unknown line")):
            r = run_cpp("codon-dnds", "-o", self.tmp / "x.tsv", stdin=broken)
            self.assertEqual(r.returncode, 1, broken)
            self.assertIn(why, r.stderr)
        self.assertEqual(run_cpp("codon-dnds", stdin=b"").returncode, 2)                      # no -o
        self.assertEqual(run_cpp("codon-dnds", "-o", self.tmp / "x", "-m", "GY", stdin=b"").returncode, 2)
        self.assertEqual(run_cpp("codon-dnds", "-o", self.tmp / "x", "--gap-open", "1", "--gap-extend", "2",
                                 stdin=b"").returncode, 2)

    def test_python_codon_dnds(self):
        out = self.tmp / "o.tsv"
        n = P.codon_dnds([("a", [("g02", self.q["g02"]), ("h02", self.s["h02"])]),
                          ("b", [("g03", self.q["g03"]), ("h03", self.s["h03"])], "rep")], str(out), method=["NG"],
                         group_column="grp")
        self.assertEqual(n, 2)
        rows = read_tsv(out)
        self.assertEqual([(r["grp"], r["seq1"], r["method"], r["status"]) for r in rows],
                         [("a", "g02", "NG", "ok"), ("b", "g03", "NG", "ok")])


# ---------------------------------------------------------------------------------------------- pair mode
@needs_cpp
class PairsGiven(PipelineCase):
    """Pair mode with --pairs (no DIAMOND needed)."""

    def test_run_matches_the_steps_one_by_one(self):
        pairs = PLAIN + EDGE
        res = P.run_pairs(Q, S, self.tmp / "out", pairs=self.pairs_file(pairs), method=["YN", "MYN"], jobs=2,
                          kaks_format=True, keep_alignments=True, quiet=True)
        self.assertEqual(res.exit_code, 2)            # the query file holds genes that cannot be used
        self.check_report(res)
        rows = read_tsv(res.outputs["dnds"])
        self.assertEqual(list(rows[0]), list(P.PAIR_COLUMNS))
        self.assertEqual([(r["seq1"], r["seq2"]) for r in rows[::2]], pairs)
        hand = by_hand(pairs, read_fasta(Q), read_fasta(S), ("YN", "MYN"))
        for r in rows:
            h = hand[(r["seq1"], r["seq2"], r["method"])]
            self.assertEqual([r[c] for c in VALUES], [h[c] for c in VALUES], r["seq1"])
        c = res.report["counts"]
        self.assertEqual(c["genes"], {"total": 40, "ok": 34, "warning": 2, "failed": 4})
        self.assertEqual(c["pairs"], {"total": 17, "ok": 15, "warning": 2, "saturated": 0, "no_data": 0, "error": 0})
        failed = read_tsv(res.outputs["failed"])
        self.assertEqual(list(failed[0]), list(P.FAILED_COLUMNS))
        self.assertEqual(sorted((f["level"], f["seq1"], f["code"]) for f in failed),
                         [("gene", "bad_q", "cds_invalid_character"), ("gene", "empty_q", "cds_empty"),
                          ("gene", "g01", "duplicate_id"), ("gene", "short_s", "cds_too_short")])
        self.assertEqual(res.report["failures_total"], 4)
        genes = read_tsv(res.outputs["genes"])
        self.assertEqual({g["set"] for g in genes}, {"query", "subject"})
        self.assertEqual(len(genes), 40)
        kept = Path(res.outputs["codon_alignments"]).read_text().splitlines()
        self.assertEqual([ln[1:] for ln in kept[0::4]], [a for a, _ in pairs])
        kaks = Path(res.outputs["kaks"]).read_text().splitlines()
        self.assertEqual(len(kaks), 1 + 2 * len(pairs))
        self.assertEqual([s["name"] for s in res.report["stages"]],
                         ["translate_query", "translate_subject", "orthologs", "align_codon_kaks"])

    def test_bad_pairs_are_rows_not_crashes(self):
        pairs = [("g02", "h02"), ("bad_q", "h03"), ("g04", "ghost"), ("g05", "short_s"), ("g06", "h06")]
        res = P.run_pairs(Q, S, self.tmp / "out", pairs=self.pairs_file(pairs), quiet=True)
        self.assertEqual(res.exit_code, 2)
        rows = read_tsv(res.outputs["dnds"])
        self.assertEqual([r["status"] for r in rows], ["ok", "error", "error", "error", "ok"])
        self.assertIn("seq1: the CDS holds a character that is no nucleotide code", rows[1]["note"])
        self.assertIn("seq2: no CDS record with this ID in the subject file", rows[2]["note"])
        self.assertIn("shorter than one codon", rows[3]["note"])
        self.assertEqual(res.report["counts"]["pairs"]["error"], 3)
        pf = [f for f in read_tsv(res.outputs["failed"]) if f["level"] == "pair"]
        self.assertEqual([(f["seq1"], f["seq2"]) for f in pf], pairs[1:4])
        self.check_report(res)

    def test_exit_0_jobs_and_resume(self):
        q, s = self.tmp / "q.fa", self.tmp / "s.fa"
        allq, alls = read_fasta(Q), read_fasta(S)
        q.write_text("".join(f">{a}\n{allq[a]}\n" for a, _ in PLAIN))
        s.write_text("".join(f">{b}\n{alls[b]}\n" for _, b in PLAIN))
        pf = self.pairs_file(PLAIN)
        outs = []
        for j in (1, 3):
            res = P.run_pairs(q, s, self.tmp / f"o{j}", pairs=pf, jobs=j, quiet=True)
            self.assertEqual(res.exit_code, 0)
            self.assertEqual(res.report["failures"], [])
            outs.append(Path(res.outputs["dnds"]).read_bytes())
        self.assertEqual(outs[0], outs[1])
        again = P.run_pairs(q, s, self.tmp / "o1", pairs=pf, jobs=1, quiet=True, resume=True)
        self.assertEqual([st["cached"] for st in again.report["stages"]], [True, True, False, False])
        self.assertEqual(Path(again.outputs["dnds"]).read_bytes(), outs[0])
        # a changed input is translated again
        q.write_text(q.read_text() + f">extra\n{allq['same_q']}\n")
        again = P.run_pairs(q, s, self.tmp / "o1", pairs=pf, jobs=1, quiet=True, resume=True)
        self.assertEqual([st["cached"] for st in again.report["stages"]][:2], [False, True])
        self.assertEqual(again.report["counts"]["genes"]["total"], 25)
        # the pwalign preset reaches the aligner
        other = P.run_pairs(q, s, self.tmp / "o_pw", pairs=pf, preset="pwalign", quiet=True)
        self.assertEqual(other.report["parameters"]["gap_open"], 14.0)
        self.assertNotEqual(Path(other.outputs["dnds"]).read_bytes(), outs[0])

    def test_runs_that_cannot_start(self):
        with redirect_stderr(io.StringIO()) as err:
            res = P.run_pairs(self.tmp / "missing.fa", S, self.tmp / "out", pairs=self.pairs_file(PLAIN))
        self.assertEqual(res.exit_code, 1)
        self.assertIn("missing.fa", err.getvalue())
        self.assertEqual(json.loads((self.tmp / "out" / "report.json").read_text())["exit_code"], 1)
        blocker = self.tmp / "file"
        blocker.write_text("x")
        with redirect_stderr(io.StringIO()):
            res = P.run_pairs(Q, S, blocker / "out", pairs=self.pairs_file(PLAIN))
        self.assertEqual(res.exit_code, 1)
        with self.assertRaises(TypeError):
            P.run_pairs(Q, S, self.tmp / "o", no_such_option=1)
        with self.assertRaises(ValueError):
            P.run_pairs(Q, S, self.tmp / "o", method=["GY"])

    def test_command_line(self):
        r = subprocess.run([sys.executable, "-m", "coal", "pipeline", "pairs", "--query", str(Q), "--subject",
                            str(S), "--pairs", str(self.pairs_file(PLAIN[1:3])), "-o", str(self.tmp / "out"),
                            "-m", "ng"], capture_output=True, text=True,
                           env={**os.environ, "PYTHONPATH": str(ROOT / "src")})
        self.assertEqual(r.returncode, 2, r.stderr)   # finished; the unusable genes of the query file are failures
        self.assertEqual(r.stdout, "")
        self.assertIn("coal pipeline: summary: 40 genes (4 not usable); 2 pairs: 2 ok", r.stderr)
        rows = read_tsv(self.tmp / "out" / "dnds.tsv")
        self.assertEqual([x["method"] for x in rows], ["NG", "NG"])


@needs_diamond
class PairsDiamond(PipelineCase):
    def test_reciprocal_best_hits(self):
        res = P.run_pairs(Q, S, self.tmp / "out", jobs=2, quiet=True)
        self.assertEqual(res.exit_code, 2)
        self.check_report(res)
        rows = read_tsv(res.outputs["dnds"])
        self.assertEqual([(r["seq1"], r["seq2"]) for r in rows], PLAIN + EDGE)   # query-file order
        self.assertTrue(all(r["status"] in ("ok", "warning") for r in rows))
        self.assertEqual(read_tsv(res.outputs["unpaired"]),
                         [{"set": "query", "id": "nohit_q", "reason": "no_hit"},
                          {"set": "subject", "id": "nohit_s", "reason": "no_hit"}])
        orth = read_tsv(res.outputs["orthologs"])
        self.assertEqual(list(orth[0]), ["query", "subject", "pident", "length", "evalue", "bitscore"])
        self.assertEqual(res.report["counts"]["orthologs"], 17)
        self.assertEqual(res.report["tools"]["diamond"]["version"], P._diamond_version(DIAMOND))
        # the same pairs given as a file give the same table
        given = P.run_pairs(Q, S, self.tmp / "given", pairs=res.outputs["orthologs"], quiet=True)
        self.assertEqual(Path(given.outputs["dnds"]).read_bytes(), Path(res.outputs["dnds"]).read_bytes())
        # --resume reuses the search
        again = P.run_pairs(Q, S, self.tmp / "out", jobs=2, quiet=True, resume=True)
        self.assertEqual({s["name"]: s["cached"] for s in again.report["stages"]}["diamond_search"], True)
        self.assertEqual(Path(again.outputs["dnds"]).read_bytes(), Path(given.outputs["dnds"]).read_bytes())

    def test_best_hits_and_the_diamond_alignment(self):
        bh = P.run_pairs(Q, S, self.tmp / "bh", ortho="bh", quiet=True)
        self.assertEqual(len(read_tsv(bh.outputs["dnds"])), 17)
        self.assertEqual([u["set"] for u in read_tsv(bh.outputs["unpaired"])], ["query"])
        dm = P.run_pairs(Q, S, self.tmp / "dm", aligner="diamond", quiet=True)
        rows = read_tsv(dm.outputs["dnds"])
        self.assertEqual(len(rows), 17)
        self.assertTrue(all(r["status"] in ("ok", "warning") for r in rows))
        with redirect_stderr(io.StringIO()) as err:
            bad = P.run_pairs(Q, S, self.tmp / "x", aligner="diamond", pairs=bh.outputs["orthologs"])
        self.assertEqual(bad.exit_code, 1)
        self.assertIn("cannot be combined with --pairs", err.getvalue())

    def test_missing_diamond_is_explained(self):
        with mock.patch.object(P, "find_diamond", return_value=None), redirect_stderr(io.StringIO()) as err:
            res = P.run_pairs(Q, S, self.tmp / "out")
        self.assertEqual(res.exit_code, 1)
        self.assertIn("--install-tools", err.getvalue())


@needs_clustalo
class PairsMsaTool(PipelineCase):
    def test_clustalo_per_pair(self):
        pairs = PLAIN[:4] + [("bad_q", "h05")]
        res = P.run_pairs(Q, S, self.tmp / "out", pairs=self.pairs_file(pairs), aligner="clustalo", jobs=2, quiet=True)
        rows = read_tsv(res.outputs["dnds"])
        self.assertEqual([r["status"] for r in rows], ["ok"] * 4 + ["error"])
        self.assertIn("msa_seconds", res.report["stages"][-1])
        self.assertIn("clustalo", res.report["tools"])


# ---------------------------------------------------------------------------------------------- cluster mode
def recompute_summary(rows: list[dict]) -> dict:
    """{(cluster, method): expected summary cells} from dnds.tsv rows, with the statistics module."""
    out = {}
    groups: dict = {}
    for r in rows:
        groups.setdefault((r["cluster_id"], r["method"]), []).append(r)
    for key, g in groups.items():
        ok = [r for r in g if r["status"] in ("ok", "warning")]
        e = {"n_pairs": len(g), "n_ok": len(ok), "n_saturated": sum(r["status"] == "saturated" for r in g),
             "n_no_data": sum(r["status"] == "no_data" for r in g), "n_failed": sum(r["status"] == "error" for r in g)}
        for c in ("dN", "dS", "dNdS"):
            v = sorted(float(r[c]) for r in ok if r[c] != "NA")
            if c == "dNdS":
                e["n_dNdS"] = len(v)
            if not v:
                continue
            q = statistics.quantiles(v, n=4, method="inclusive") if len(v) > 1 else [v[0]] * 3
            e.update({f"{c}_median": statistics.median(v), f"{c}_mean": math.fsum(v) / len(v), f"{c}_q25": q[0],
                      f"{c}_q75": q[2], f"{c}_min": v[0], f"{c}_max": v[-1]})
        out[key] = e
    return out


@needs_cpp
class ClustersGiven(PipelineCase):
    """Cluster mode with --clusters and a size cap of 2, so that no MSA program is needed."""

    def run_it(self, name="out", **kw):
        kw.setdefault("max_cluster_size", 2)
        return P.run_clusters(FAM, self.tmp / name, clusters=FAM_CLUSTERS, quiet=True, **kw)

    def check_summary(self, res):
        rows = read_tsv(res.outputs["dnds"])
        summ = read_tsv(res.outputs["cluster_summary"])
        self.assertEqual(list(summ[0]), list(P.SUMMARY_COLUMNS))
        want = recompute_summary(rows)
        seen = set()
        for s in summ:
            e = want.get((s["cluster_id"], s["method"]))
            if e is None:
                self.assertEqual((s["n_pairs"], s["dN_median"], s["method"]), ("0", "NA", "NA"))
                continue
            seen.add((s["cluster_id"], s["method"]))
            for k, v in e.items():
                if isinstance(v, int):
                    self.assertEqual(int(s[k]), v, (s["cluster_id"], k))
                else:
                    self.assertAlmostEqual(float(s[k]), v, delta=1e-7 * max(1.0, abs(v)), msg=(s["cluster_id"], k))
        self.assertEqual(seen, set(want))
        return rows, summ

    def test_representative_and_singletons(self):
        res = self.run_it(method=["YN", "NG"], jobs=2)
        self.assertEqual(res.exit_code, 2)            # famB holds two members that cannot be used
        self.check_report(res)
        rows, summ = self.check_summary(res)
        s = {(x["cluster_id"], x["method"]): x for x in summ}
        a = s[("famA_001", "YN")]
        self.assertEqual((a["n_members"], a["n_used"], a["pair_selection"], a["aligner"], a["n_pairs"], a["status"]),
                         ("6", "6", "representative", "nw_pairwise", "5", "ok"))
        self.assertIn("above --max-cluster-size 2", a["note"])
        b = s[("famB_001", "NG")]
        self.assertEqual((b["n_members"], b["n_used"], b["n_pairs"]), ("13", "11", "10"))
        self.assertIn("2 member(s) left out", b["note"])
        c = s[("famC_001", "YN")]
        self.assertEqual((c["pair_selection"], c["aligner"], c["n_pairs"]), ("all", "nw", "1"))
        d = s[("famD_001", "NA")]
        self.assertEqual((d["n_members"], d["pair_selection"], d["aligner"], d["n_pairs"], d["dS_mean"], d["status"]),
                         ("1", "singleton", "none", "0", "NA", "ok"))
        self.assertEqual([x["cluster_id"] for x in summ[::2]][:3], ["famA_001", "famB_001", "famC_001"])
        self.assertEqual(list(rows[0])[0], "cluster_id")
        self.assertTrue(all(r["seq1"] == r["cluster_id"] for r in rows))      # representative against the others
        fam = read_fasta(FAM)
        hand = by_hand([(r["seq1"], r["seq2"]) for r in rows[::2]], fam, fam, ("YN", "NG"))
        for r in rows:
            h = hand[(r["seq1"], r["seq2"], r["method"])]
            self.assertEqual([r[k] for k in VALUES], [h[k] for k in VALUES])
        failed = read_tsv(res.outputs["failed"])
        self.assertEqual(sorted((f["level"], f["cluster_id"], f["seq1"], f["code"]) for f in failed),
                         [("gene", "", "famB_007", "cds_invalid_character"),
                          ("gene", "famB_001", "famB_007", "gene_failed"),
                          ("gene", "famB_001", "famB_ghost", "no_cds")])
        self.assertEqual(res.report["counts"]["clusters"],
                         {"total": 4, "singletons": 1, "all_pairs": 1, "above_max_size": 2, "skipped": 0,
                          "msa_fallback": 0, "members_dropped": 2})

    def test_sample_is_reproducible(self):
        a = self.run_it("a", large_clusters="sample", pair_sample=7, seed=5, jobs=1)
        b = self.run_it("b", large_clusters="sample", pair_sample=7, seed=5, jobs=3)
        c = self.run_it("c", large_clusters="sample", pair_sample=7, seed=6, jobs=1)
        self.assertEqual(Path(a.outputs["dnds"]).read_bytes(), Path(b.outputs["dnds"]).read_bytes())
        self.assertEqual(Path(a.outputs["cluster_summary"]).read_bytes(),
                         Path(b.outputs["cluster_summary"]).read_bytes())
        self.assertNotEqual(Path(a.outputs["dnds"]).read_bytes(), Path(c.outputs["dnds"]).read_bytes())
        rows, summ = self.check_summary(a)
        s = {x["cluster_id"]: x for x in summ}
        self.assertEqual((s["famA_001"]["pair_selection"], s["famA_001"]["n_pairs"]), ("sample:7", "7"))
        self.assertEqual((s["famB_001"]["pair_selection"], s["famB_001"]["n_pairs"]), ("sample:7", "7"))
        pairs_a = [(r["seq1"], r["seq2"]) for r in rows if r["cluster_id"] == "famA_001"]
        self.assertEqual(len(set(pairs_a)), 7)
        # more pairs asked for than there are: all of them
        big = self.run_it("big", large_clusters="sample", pair_sample=1000)
        s = {x["cluster_id"]: x for x in read_tsv(big.outputs["cluster_summary"])}
        self.assertEqual((s["famA_001"]["pair_selection"], s["famA_001"]["n_pairs"]), ("sample:15", "15"))

    def test_skip(self):
        res = self.run_it(large_clusters="skip")
        rows, summ = self.check_summary(res)
        self.assertEqual({r["cluster_id"] for r in rows}, {"famC_001"})
        s = {x["cluster_id"]: x for x in summ}
        self.assertEqual((s["famA_001"]["pair_selection"], s["famA_001"]["n_pairs"], s["famA_001"]["status"]),
                         ("skipped", "0", "skipped"))
        self.assertEqual(res.report["counts"]["clusters"]["skipped"], 2)

    def test_several_cds_files_and_kept_alignments(self):
        fam = read_fasta(FAM)
        one, two = self.tmp / "one.fa", self.tmp / "two.fa"
        names = [n for n in fam if not n.startswith("famB")]
        one.write_text("".join(f">{n}\n{fam[n]}\n" for n in names[:4]))
        two.write_text("".join(f">{n}\n{fam[n]}\n" for n in names[4:]))
        cl = self.tmp / "cl.tsv"
        cl.write_text("".join(ln + "\n" for ln in FAM_CLUSTERS.read_text().splitlines() if not ln.startswith("famB")))
        res = P.run_clusters([one, two], self.tmp / "out", clusters=cl, max_cluster_size=2, keep_alignments=True,
                             kaks_format=True, quiet=True)
        self.assertEqual(res.exit_code, 0)
        self.assertEqual(res.report["counts"]["genes"]["total"], 9)
        self.assertEqual(sorted(os.listdir(res.outputs["codon_alignments"])), ["famC_001.fasta"])
        self.assertEqual(len(Path(res.outputs["kaks"]).read_text().splitlines()), 1 + 5 + 1)

    @unittest.skipIf(os.name == "nt", "the stand-in tool is a POSIX shell script")
    def test_a_failing_msa_falls_back_to_pairwise(self):
        tool = self.tmp / "clustalo"
        tool.write_text("#!/bin/sh\necho 'cannot align' >&2\nexit 1\n")
        tool.chmod(tool.stat().st_mode | stat.S_IXUSR)
        with mock.patch.dict(os.environ, {"COAL_CLUSTALO": str(tool)}):
            res = self.run_it(max_cluster_size=8)
        self.assertEqual(res.exit_code, 2)
        rows, summ = self.check_summary(res)
        s = {x["cluster_id"]: x for x in summ}
        self.assertEqual((s["famA_001"]["pair_selection"], s["famA_001"]["aligner"], s["famA_001"]["n_pairs"],
                          s["famA_001"]["n_ok"]), ("all", "nw_pairwise", "15", "15"))
        self.assertIn("clustalo failed", s["famA_001"]["note"])
        self.assertEqual(res.report["counts"]["clusters"]["msa_fallback"], 1)
        self.assertIn(("cluster", "famA_001", "msa_failed"),
                      [(f["level"], f["cluster_id"], f["code"]) for f in read_tsv(res.outputs["failed"])])

    def test_missing_msa_tool_is_explained(self):
        with mock.patch.object(msa, "find_clustalo", return_value=None), \
                mock.patch.dict(os.environ, {"COAL_CLUSTALO": ""}), redirect_stderr(io.StringIO()) as err:
            res = self.run_it(max_cluster_size=8)
        self.assertEqual(res.exit_code, 1)
        self.assertIn("Clustal Omega was not found", err.getvalue())


@needs_clustalo
class ClustersMsa(PipelineCase):
    def check(self, tool):
        res = P.run_clusters(FAM, self.tmp / tool, clusters=FAM_CLUSTERS, max_cluster_size=8, aligner=tool, jobs=2,
                             keep_alignments=True, quiet=True)
        self.assertEqual(res.exit_code, 2)
        rows = read_tsv(res.outputs["dnds"])
        summ = {x["cluster_id"]: x for x in read_tsv(res.outputs["cluster_summary"])}
        self.assertEqual((summ["famA_001"]["pair_selection"], summ["famA_001"]["aligner"],
                          summ["famA_001"]["n_pairs"]), ("all", tool, "15"))
        self.assertEqual((summ["famB_001"]["pair_selection"], summ["famB_001"]["aligner"]),
                         ("representative", "nw_pairwise"))
        self.assertEqual(sorted(os.listdir(res.outputs["codon_alignments"])), ["famA_001.fasta", "famC_001.fasta"])
        # every pair of the cluster equals `coal-cpp dnds` on the kept codon MSA
        ref = run_cpp("dnds", Path(res.outputs["codon_alignments"]) / "famA_001.fasta", "-q").stdout.decode()
        ref_rows = [ln.split("\t") for ln in ref.splitlines()[1:]]
        mine = [r for r in rows if r["cluster_id"] == "famA_001"]
        self.assertEqual(len(mine), 15)
        for r, h in zip(mine, ref_rows):
            self.assertEqual([r["seq1"], r["seq2"], r["method"]] + [r[c] for c in VALUES], h[:13])
        self.assertIn(tool, res.report["tools"])
        return res

    def test_clustalo(self):
        a = self.check("clustalo")
        b = P.run_clusters(FAM, self.tmp / "j1", clusters=FAM_CLUSTERS, max_cluster_size=8, jobs=1, quiet=True)
        self.assertEqual(Path(a.outputs["dnds"]).read_bytes(), Path(b.outputs["dnds"]).read_bytes())

    @unittest.skipUnless(MAFFT, "needs MAFFT")
    def test_mafft(self):
        self.check("mafft")


@needs_diamond
class ClustersDeepClust(PipelineCase):
    def test_deepclust_finds_the_families(self):
        res = P.run_clusters(FAM, self.tmp / "out", max_cluster_size=2, realign=True, jobs=2, quiet=True)
        self.assertIn(res.exit_code, (0, 2))
        self.check_report(res)
        clusters, _ = P.read_clusters(res.outputs["clusters"])
        fams = sorted(sorted({m.split("_")[0] for m in mem}) for mem in clusters.values())
        self.assertEqual(fams, [["famA"], ["famB"], ["famC"], ["famD"]])
        self.assertEqual(sorted(len(m) for m in clusters.values()), [1, 2, 6, 11])
        summ = read_tsv(res.outputs["cluster_summary"])
        self.assertEqual(len(summ), 4)
        self.assertEqual(sum(int(s["n_pairs"]) for s in summ), 5 + 10 + 1)
        realign = Path(res.outputs["realign"]).read_text().splitlines()
        self.assertEqual(realign[0].split("\t")[:3], ["member", "representative", "approx_pident"])
        self.assertGreaterEqual(len(realign), 1 + 16)
        again = P.run_clusters(FAM, self.tmp / "out", max_cluster_size=2, realign=True, jobs=1, quiet=True,
                               resume=True)
        st = {s["name"]: s["cached"] for s in again.report["stages"]}
        self.assertTrue(st["deepclust"] and st["realign"] and st["translate_proteins"])


if __name__ == "__main__":
    unittest.main()
