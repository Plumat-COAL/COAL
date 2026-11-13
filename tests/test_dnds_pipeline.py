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
import signal
import stat
import statistics
import subprocess
import sys
import tarfile
import tempfile
import threading
import time
import unittest
import zipfile
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
SCHEMA = json.loads((ROOT / "docs" / "pipeline.schema.json").read_text(encoding="utf-8"))
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


def needs(have, why: str, *require_env: str):
    """Skip the test or class without `have`; fail it instead when one of the environment variables is set (CI's
    C++ job sets COAL_REQUIRE_CPP, and COAL_REQUIRE_DIAMOND after installing DIAMOND: there a missing tool must
    not turn the tests of the default route into skips)."""
    if have:
        return lambda obj: obj
    required = [e for e in require_env if os.environ.get(e)]
    if not required:
        return unittest.skip(why)

    def deco(obj):
        def fail(self, *a, **k):
            self.fail(f"{why} ({' / '.join(required)} is set)")
        if isinstance(obj, type):
            obj.setUp = fail
            return obj
        return fail
    return deco


needs_cpp = needs(CPP, "build/coal-cpp (COAL_CPP_BIN) is not built", "COAL_REQUIRE_CPP")
# A missing coal-cpp fails where COAL_REQUIRE_CPP is set; a missing DIAMOND only where COAL_REQUIRE_DIAMOND is
# (CI's macOS job builds coal-cpp but installs no DIAMOND: there these tests skip).
needs_diamond = needs(CPP and DIAMOND, "needs coal-cpp and DIAMOND",
                      *(("COAL_REQUIRE_DIAMOND",) if CPP else ("COAL_REQUIRE_CPP", "COAL_REQUIRE_DIAMOND")))
needs_clustalo = unittest.skipUnless(CPP and CLUSTALO, "needs coal-cpp and Clustal Omega")


def read_fasta(path) -> dict[str, str]:
    """{id: sequence}; of a repeated ID the first record."""
    recs, cur = {}, None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = (line[1:].split() or [""])[0]
            cur = None if name in recs else recs.setdefault(name, [])
        elif cur is not None:
            cur.append(line.strip())
    return {k: "".join(v) for k, v in recs.items()}


def read_tsv(path) -> list[dict]:
    lines = Path(path).read_text(encoding="utf-8").splitlines()
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
        with open(stacked, "w", encoding="utf-8") as fs:
            for (x, y), al in zip(pairs, alns):
                pep, nuc = os.path.join(td, "pep.fa"), os.path.join(td, "nuc.fa")
                Path(pep).write_text(f">a\n{al.a}\n>b\n{al.b}\n", encoding="utf-8")
                Path(nuc).write_text(f">a\n{q[x]}\n>b\n{s[y]}\n", encoding="utf-8")
                r = subprocess.run([str(CPP), pep, nuc, "-q"], capture_output=True, text=True, check=True, encoding="utf-8")
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
        p.write_text("".join(f"{a}\t{b}\n" for a, b in pairs), encoding="utf-8")
        return p

    def check_report(self, res: P.PipelineResult):
        validate(res.report, SCHEMA)
        on_disk = json.loads(Path(res.outputs["report"]).read_text(encoding="utf-8"))
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
            p.write_text("centroid\tmember\n# a comment\nr1\tr1\nr1\ta\nr2\tb\n\nr1\tc\nr2\ta\nr3\tr3\r\n", encoding="utf-8")
            clusters, dup = P.read_clusters(str(p))
            self.assertEqual(clusters, {"r1": ["r1", "a", "c"], "r2": ["r2", "b"], "r3": ["r3"]})
            self.assertEqual(dup, 1)
            p.write_text("only_one_column\n", encoding="utf-8")
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
                    with open(dest + ".fetched", "wb") as fh:
                        fh.write(data)
                    return dest + ".fetched", hashlib.sha256(data).hexdigest()
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
            # the installed binary has a checksum record; one that changed since, or has none, is not used
            manifest = json.loads(Path(exe + msa.MANIFEST_SUFFIX).read_text(encoding="utf-8"))
            self.assertEqual(manifest["files"], {"diamond": hashlib.sha256(script).hexdigest()})
            self.assertEqual(manifest["verified"], "SHA-256 given by the caller")
            with mock.patch.dict(os.environ, {"PATH": ""}), redirect_stderr(io.StringIO()) as note:
                os.environ.pop("COAL_DIAMOND", None)
                with open(exe, "ab") as fh:
                    fh.write(b"echo planted\n")
                self.assertIsNone(P.find_diamond())
                self.assertIn("is not the file COAL installed", note.getvalue())
                os.remove(exe + msa.MANIFEST_SUFFIX)
                self.assertIsNone(P.find_diamond())
            # a download that was cut short (possible only when nobody verifies it): a clear error, nothing kept
            with serving(good[:len(good) // 2], calls), self.assertRaises(P.PipelineError) as cm, \
                    redirect_stderr(io.StringIO()):
                P.install_diamond(method="official", version="9.9.9", quiet=True, force=True, allow_unverified=True)
            self.assertIn("is not a complete archive", str(cm.exception))
            self.assertIn("Nothing was installed", str(cm.exception))
            self.assertEqual([n for n in os.listdir(os.path.join(td, "tools", "diamond")) if "diamond-" in n], [])
            # a binary that does not run is refused
            bad = archive(b"#!/bin/sh\nexit 3\n")
            with serving(bad, calls), mock.patch.dict(P._DIAMOND_SHA256, {asset: hashlib.sha256(bad).hexdigest()}):
                with self.assertRaises(P.PipelineError):
                    P.install_diamond(method="official", quiet=True, force=True)

    def test_zip_member_names_cannot_leave_the_tool_folder(self):
        """Finding 5d: a name with ':' (a drive, an alternate data stream), an absolute path, `..` or a reserved
        device name is refused, on every system, and nothing of such an archive is installed."""
        for good in ("diamond.exe", "bin/libgcc_s_seh-1.dll", "x\\\\diamond.exe".replace("\\\\", "\\"), "a/b/zlib1.dll"):
            self.assertEqual(P._safe_member_name(good, "a.zip"), good.replace("\\", "/").rsplit("/", 1)[-1])
        for bad in ("C:evil.exe", "diamond.exe:stream", "x/C:/evil.exe", "/abs/diamond.exe", "\\\\srv\\share\\d.exe",
                    "../diamond.exe", "a/../../b.dll", "NUL", "con.txt", "COM1.exe", "diamond.exe.", "dia mond.exe",
                    ".hidden", "a/..", "d\u00e9.exe"):
            with self.subTest(name=bad), self.assertRaises(msa.MsaError):
                P._safe_member_name(bad, "a.zip")
        with tempfile.TemporaryDirectory() as td, mock.patch.dict(os.environ, {"COAL_CACHE_DIR": td}), \
                mock.patch.object(msa, "_platform", return_value=("Windows", "AMD64")):
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w") as zf:
                zf.writestr("diamond.exe", b"MZ")
                zf.writestr("diamond.exe:evil", b"stream")
            data = buf.getvalue()

            def fake(url, dest, quiet):
                with open(dest + ".fetched", "wb") as fh:
                    fh.write(data)
                return dest + ".fetched", hashlib.sha256(data).hexdigest()
            with mock.patch.object(msa, "_download", side_effect=fake), self.assertRaises(msa.MsaError) as cm:
                P._install_diamond_official("9.9.9", True, hashlib.sha256(data).hexdigest())
            self.assertIn("unsafe name", str(cm.exception))
            self.assertEqual(os.listdir(os.path.join(td, "tools", "diamond")), [])

    def test_concurrent_installs_wait_and_all_succeed(self):
        """Finding 5a: two installs wrote the same `<archive>.part`. Now one installs and the others wait for it."""
        if os.name == "nt":
            self.skipTest("POSIX shell stand-in for DIAMOND")
        script = b"#!/bin/sh\necho 'diamond version 2.2.8'\n"
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:gz") as tf:
            info = tarfile.TarInfo("diamond")
            info.size, info.mode = len(script), 0o755
            tf.addfile(info, io.BytesIO(script))
        data = buf.getvalue()
        calls = []

        def urlopen(url, timeout=None):
            calls.append(url)
            time.sleep(0.3)
            return io.BytesIO(data)
        with tempfile.TemporaryDirectory() as td, mock.patch.dict(os.environ, {"COAL_CACHE_DIR": td}), \
                mock.patch.object(msa, "_platform", return_value=("Linux", "x86_64")), \
                mock.patch.object(msa.urllib.request, "urlopen", urlopen), \
                mock.patch.dict(P._DIAMOND_SHA256, {("2.2.8", "diamond-linux64.tar.gz"): hashlib.sha256(data).hexdigest()}):
            got, errors = [], []

            def one():
                try:
                    got.append(P.install_diamond(method="official", quiet=True))
                except Exception as e:   # noqa: BLE001
                    errors.append(repr(e))
            threads = [threading.Thread(target=one) for _ in range(5)]
            for t in threads:
                t.start()
            for t in threads:
                t.join(60)
            self.assertEqual(errors, [])
            self.assertEqual(len(set(got)), 1)
            self.assertEqual(len(got), 5)
            self.assertEqual(len(calls), 1, "the installs that waited found the tool installed")
            self.assertEqual(sorted(os.listdir(os.path.join(td, "tools", "diamond"))),
                             ["diamond", "diamond" + msa.MANIFEST_SUFFIX])
            self.assertEqual(sorted(os.listdir(os.path.join(td, "tools"))), ["diamond"])
            # --force installs again, also while the installed binary is "running" (it is replaced, not rewritten)
            self.assertEqual(P.install_diamond(method="official", quiet=True, force=True), got[0])
            self.assertEqual(len(calls), 2)

    def test_install_command_takes_a_hash_or_explicit_permission(self):
        seen = {}

        def fake_diamond(**kw):
            seen.update(kw)
            return "/cache/diamond"
        with mock.patch.object(P, "install_diamond", fake_diamond), redirect_stdout(io.StringIO()):
            self.assertEqual(P.main(["install", "diamond", "--sha256", "A" * 64, "-q"]), 0)
            self.assertEqual((seen["sha256"], seen["allow_unverified"]), ("A" * 64, False))
            self.assertEqual(P.main(["install", "diamond", "--allow-unverified", "-q"]), 0)
            self.assertEqual((seen["sha256"], seen["allow_unverified"]), (None, True))
        # Clustal Omega: Bioconda only. The options of the former download of the Clustal site (gone) are refused
        # with the explanation, before anything is installed or fetched
        for extra in (["--sha256", "a" * 64], ["--allow-unverified"]):
            with mock.patch.object(msa, "_conda_install", side_effect=AssertionError("installed")), \
                    mock.patch.object(msa.urllib.request, "urlopen", side_effect=AssertionError("fetched")), \
                    redirect_stderr(io.StringIO()) as err, redirect_stdout(io.StringIO()):
                self.assertEqual(P.main(["install", "clustalo", *extra, "-q"]), 1)
            self.assertIn("do not apply to Clustal Omega", err.getvalue())
            self.assertIn("binaries are gone", err.getvalue())
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
        f.write_text(">a\nATGTGAAAATAA\n", encoding="utf-8")                                    # TGA: stop (1), Trp (4)
        self.assertIn(b"M*K", run_cpp("translate", f, "-q").stdout)
        self.assertIn(b"MWK", run_cpp("translate", f, "-q", "-c", "4").stdout)
        self.assertEqual(run_cpp("translate", self.tmp / "missing.fa").returncode, 1)
        f.write_text("ATGAAA\n", encoding="utf-8")
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
        kept = (self.tmp / "c.fa").read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(kept), 4 * len(pairs))
        k = kept.index(">rna_q")
        self.assertTrue(kept[k + 1].startswith("aug"))                        # the DNA is copied as written
        self.assertEqual(len(kept[k + 1]), len(kept[k + 3]))

    def test_iteration_that_does_not_converge_is_not_an_estimate(self):
        """AAA against ATG, YN: the iteration does not converge (tests/test_kaks.py). The pair has status
        not_converged and no dN, dS; the KaKs_Calculator-format file has NA there as well. Nothing platform
        dependent is asserted (the last round, kept in the note, differs between platforms)."""
        x, y = PLAIN[1]
        # (protein rows given: the pairwise aligner would not put K against M)
        data = stream([("nc", "all", [("M", "a", "AAA", "K"), ("M", "b", "ATG", "M")]),
                       ("ok", "all", [("M", x, self.q[x]), ("M", y, self.s[y])])])
        r, rows = self.run_stream(data, "-m", "YN", "-m", "NG", "--kaks-out", self.tmp / "k.kaks")
        self.assertEqual(r.returncode, 0, r.stderr)
        st = {(x["seq1"], x["method"]): x for x in rows}
        nc = st[("a", "YN")]
        self.assertEqual(nc["status"], "not_converged")
        self.assertEqual([nc[c] for c in ("dN", "dS", "dNdS", "SE_dN", "SE_dS", "Nd", "Sd")], ["NA"] * 7)
        self.assertIn("did not converge", nc["note"])
        self.assertNotEqual(st[("a", "NG")]["status"], "not_converged")     # no iteration there
        self.assertEqual(st[(x, "YN")]["status"], "ok")
        k = {(f[0], f[1]): f for f in (line.split("\t") for line in
                                       (self.tmp / "k.kaks").read_text(encoding="utf-8").splitlines()[1:])}
        self.assertEqual(k[("a&b", "YN")][2:6], ["NA"] * 4)
        self.assertNotEqual(k[(f"{x}&{y}", "YN")][3], "NA")

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
        kaks = (self.tmp / "k.kaks").read_text(encoding="utf-8").splitlines()
        self.assertTrue(kaks[0].startswith("Sequence\tMethod\tKa\tKs\tKa/Ks\tP-Value(Fisher)"))
        self.assertEqual([ln.split("\t")[0] for ln in kaks[1:]], ["g05&h05", "h06&g06", "g06&h06"])
        self.assertEqual((self.tmp / "c.fa").read_text(encoding="utf-8").count(">"), 4)       # only the two whole pairs

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
        self.assertEqual(c["pairs"], {"total": 17, "ok": 15, "warning": 2, "not_converged": 0, "saturated": 0,
                                      "no_data": 0, "error": 0})
        failed = read_tsv(res.outputs["failed"])
        self.assertEqual(list(failed[0]), list(P.FAILED_COLUMNS))
        self.assertEqual(sorted((f["level"], f["seq1"], f["code"]) for f in failed),
                         [("gene", "bad_q", "cds_invalid_character"), ("gene", "empty_q", "cds_empty"),
                          ("gene", "g01", "duplicate_id"), ("gene", "short_s", "cds_too_short")])
        self.assertEqual(res.report["failures_total"], 4)
        genes = read_tsv(res.outputs["genes"])
        self.assertEqual({g["set"] for g in genes}, {"query", "subject"})
        self.assertEqual(len(genes), 40)
        kept = Path(res.outputs["codon_alignments"]).read_text(encoding="utf-8").splitlines()
        self.assertEqual([ln[1:] for ln in kept[0::4]], [a for a, _ in pairs])
        kaks = Path(res.outputs["kaks"]).read_text(encoding="utf-8").splitlines()
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
        q.write_text("".join(f">{a}\n{allq[a]}\n" for a, _ in PLAIN), encoding="utf-8")
        s.write_text("".join(f">{b}\n{alls[b]}\n" for _, b in PLAIN), encoding="utf-8")
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
        q.write_text(q.read_text(encoding="utf-8") + f">extra\n{allq['same_q']}\n", encoding="utf-8")
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
        self.assertEqual(json.loads((self.tmp / "out" / "report.json").read_text(encoding="utf-8"))["exit_code"], 1)
        blocker = self.tmp / "file"
        blocker.write_text("x", encoding="utf-8")
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
                           env={**os.environ, "PYTHONPATH": str(ROOT / "src")}, encoding="utf-8")
        self.assertEqual(r.returncode, 2, r.stderr)   # finished; the unusable genes of the query file are failures
        self.assertEqual(r.stdout, "")
        self.assertIn("coal pipeline: summary: 40 genes (4 not usable); 2 pairs: 2 ok", r.stderr)
        rows = read_tsv(self.tmp / "out" / "dnds.tsv")
        self.assertEqual([x["method"] for x in rows], ["NG", "NG"])


@needs_cpp
class Robustness(PipelineCase):
    """Review findings of 2026-11-12 (e): stale outputs, one run per folder, interrupts, inputs that cannot work."""

    def small(self):
        q, s = self.tmp / "q.fa", self.tmp / "s.fa"
        allq, alls = read_fasta(Q), read_fasta(S)
        q.write_text("".join(f">{a}\n{allq[a]}\n" for a, _ in PLAIN), encoding="utf-8")
        s.write_text("".join(f">{b}\n{alls[b]}\n" for _, b in PLAIN), encoding="utf-8")
        return q, s, self.pairs_file(PLAIN)

    def test_outputs_of_an_earlier_run_are_removed(self):
        q, s, pf = self.small()
        out = self.tmp / "out"
        first = P.run_pairs(q, s, out, pairs=pf, quiet=True, kaks_format=True, keep_alignments=True)
        self.assertEqual(first.exit_code, 0)
        for name in ("dnds.kaks", "codon_alignments.fasta", "dnds.tsv"):
            self.assertTrue((out / name).is_file(), name)
        (out / "dnds.tsv.tmp").write_text("half a table of a killed run\n", encoding="utf-8")
        # files with a run's names that no run wrote (the user's): they stay
        (out / "cluster_summary.tsv").write_text("the user's own table\n", encoding="utf-8")
        (out / "codon_alignments").mkdir()
        (out / "codon_alignments" / "mine.fasta").write_text(">x\nATG\n", encoding="utf-8")
        (out / "my_notes.txt").write_text("the user's own file\n", encoding="utf-8")
        second = P.run_pairs(q, s, out, pairs=pf, quiet=True)   # the same folder, without the optional outputs
        self.assertEqual(second.exit_code, 0)
        left = sorted(p.name for p in out.iterdir())
        self.assertEqual(left, ["cluster_summary.tsv", "codon_alignments", "dnds.tsv", "failed.tsv", "genes.tsv",
                                "my_notes.txt", "orthologs.tsv", "report.json", "work"])
        self.assertEqual(sorted(p.name for p in (out / "codon_alignments").iterdir()), ["mine.fasta"])
        self.assertEqual((out / "cluster_summary.tsv").read_text(encoding="utf-8"), "the user's own table\n")
        self.assertEqual(sorted(second.report["outputs"]), ["dnds", "failed", "genes", "orthologs"])
        self.assertEqual((out / "dnds.tsv").read_bytes(), Path(first.outputs["dnds"]).read_bytes())
        # --resume reuses work/ and still leaves no optional output of the earlier run
        third = P.run_pairs(q, s, out, pairs=pf, quiet=True, kaks_format=True)
        fourth = P.run_pairs(q, s, out, pairs=pf, quiet=True, resume=True)
        self.assertTrue(Path(third.outputs["kaks"]).name == "dnds.kaks" and not (out / "dnds.kaks").exists())
        self.assertEqual([st["cached"] for st in fourth.report["stages"]][:2], [True, True])

    def test_outputs_of_a_run_without_the_record_are_known_from_its_report(self):
        """A folder written before work/outputs.json existed: the outputs its report.json lists are that run's."""
        q, s, pf = self.small()
        out = self.tmp / "out"
        self.assertEqual(P.run_pairs(q, s, out, pairs=pf, quiet=True, kaks_format=True).exit_code, 0)
        (out / "work" / "outputs.json").unlink()
        (out / "unpaired.tsv").write_text("not listed in the report: the user's\n", encoding="utf-8")
        self.assertEqual(P.run_pairs(q, s, out, pairs=pf, quiet=True).exit_code, 0)
        self.assertFalse((out / "dnds.kaks").exists())
        self.assertTrue((out / "unpaired.tsv").is_file())
        # a report that cannot be read names nothing: nothing is removed
        (out / "work" / "outputs.json").unlink()
        (out / "dnds.kaks").write_text("whose?\n", encoding="utf-8")
        (out / "report.json").write_text("{not json", encoding="utf-8")
        self.assertEqual(P.run_pairs(q, s, out, pairs=pf, quiet=True).exit_code, 0)
        self.assertEqual((out / "dnds.kaks").read_text(encoding="utf-8"), "whose?\n")

    def test_the_pairs_of_an_earlier_run_can_be_the_input(self):
        """Independent review 2026-11-13, finding 1: `--pairs out/orthologs.tsv -o out` deleted its own input
        (and the earlier dnds.tsv) and failed with "cannot open"."""
        q, s, pf = self.small()
        out = self.tmp / "out"
        first = P.run_pairs(q, s, out, pairs=pf, quiet=True)
        self.assertEqual(first.exit_code, 0)
        orth = out / "orthologs.tsv"
        orth_bytes, dnds_bytes = orth.read_bytes(), (out / "dnds.tsv").read_bytes()
        r = subprocess.run([sys.executable, "-m", "coal", "pipeline", "pairs", "--query", str(q), "--subject", str(s),
                            "--pairs", str(orth), "-o", str(out), "-q"], capture_output=True, text=True,
                           encoding="utf-8", env={**os.environ, "PYTHONPATH": str(ROOT / "src")})
        self.assertEqual((r.returncode, r.stderr), (0, ""))
        self.assertEqual(orth.read_bytes(), orth_bytes, "the input is left exactly as it was")
        self.assertEqual((out / "dnds.tsv").read_bytes(), dnds_bytes)
        rep = json.loads((out / "report.json").read_text(encoding="utf-8"))
        self.assertEqual(rep["outputs"]["orthologs"], "orthologs.tsv")
        # the user's own pairs file under that name (no header, a third column) is not rewritten either
        mine = "".join(f"{a}\t{b}\tnote\n" for a, b in PLAIN[:5])
        orth.write_text(mine, encoding="utf-8")
        third = P.run_pairs(q, s, out, pairs=orth, quiet=True)
        self.assertEqual((third.exit_code, third.report["counts"]["orthologs"]), (0, 5))
        self.assertEqual(orth.read_text(encoding="utf-8"), mine)
        # by another name (a link) it is the same file
        if hasattr(os, "symlink") and os.name != "nt":
            link = self.tmp / "link.tsv"
            link.symlink_to(orth)
            self.assertEqual(P.run_pairs(q, s, out, pairs=link, quiet=True).exit_code, 0)
            self.assertEqual(orth.read_text(encoding="utf-8"), mine)
        # and a later run that searches or is given other pairs writes orthologs.tsv again
        self.assertEqual(P.run_pairs(q, s, out, pairs=pf, quiet=True).exit_code, 0)
        self.assertEqual(orth.read_bytes(), orth_bytes)

    def test_a_run_never_removes_or_writes_its_inputs(self):
        """An input under a name the run writes is refused before anything in the folder changes."""
        q, s, pf = self.small()
        out = self.tmp / "out"
        self.assertEqual(P.run_pairs(q, s, out, pairs=pf, quiet=True, kaks_format=True).exit_code, 0)

        def snapshot():
            return {str(p.relative_to(out)): p.read_bytes() for p in sorted(out.rglob("*")) if p.is_file()}
        before = snapshot()
        cases = [dict(query=q, subject=s, pairs=out / "dnds.tsv"),                 # seq1/seq2 are its first columns
                 dict(query=q, subject=s, pairs=pf, report=q),                     # --report onto the query
                 dict(query=out / "genes.tsv", subject=s, pairs=pf),
                 dict(query=q, subject=s, pairs=out / "work" / "state.json")]
        if os.name != "nt":
            (self.tmp / "alias.tsv").symlink_to(out / "dnds.tsv")
            cases.append(dict(query=q, subject=s, pairs=self.tmp / "alias.tsv"))
        for kw in cases:
            with self.subTest(**{k: str(v) for k, v in kw.items()}):
                qq, ss = kw.pop("query"), kw.pop("subject")
                with redirect_stderr(io.StringIO()) as err:
                    res = P.run_pairs(qq, ss, out, quiet=True, **kw)
                self.assertEqual(res.exit_code, 1)
                self.assertRegex(err.getvalue(), r"error: --(pairs|query) \S+ is (a file|inside \S+, a folder) this "
                                                 r"run writes")
                self.assertEqual(snapshot(), before, "nothing in the folder changed")
                self.assertFalse(Path(str(out) + ".coal-lock").exists())
                self.assertEqual(q.read_bytes()[:1], b">")

    def test_cluster_run_into_the_users_project_folder(self):
        """Finding 1, second command: `clusters --cds proj/cds.fasta --clusters proj/clusters.tsv -o proj` deleted
        proj/clusters.tsv and every proj/codon_alignments/*.fasta the user had."""
        proj = self.tmp / "proj"
        (proj / "codon_alignments").mkdir(parents=True)
        cds, clu = proj / "family_cds.fasta", proj / "clusters.tsv"
        cds.write_bytes(FAM.read_bytes())
        clu.write_bytes(FAM_CLUSTERS.read_bytes())
        mine = proj / "codon_alignments" / "famA_001.fasta"
        mine.write_text(">my own alignment\nATG\n", encoding="utf-8")
        (proj / "orthologs.tsv").write_text("mine too\n", encoding="utf-8")

        def run(**kw):
            return P.run_clusters(cds, proj, clusters=clu, quiet=True, max_cluster_size=2, **kw)
        res = run()
        self.assertEqual(res.exit_code, 2)
        self.assertEqual(clu.read_bytes(), FAM_CLUSTERS.read_bytes())
        self.assertEqual(cds.read_bytes(), FAM.read_bytes())
        self.assertEqual(mine.read_text(encoding="utf-8"), ">my own alignment\nATG\n")
        self.assertEqual((proj / "orthologs.tsv").read_text(encoding="utf-8"), "mine too\n")
        self.assertTrue((proj / "dnds.tsv").is_file())
        # --keep-alignments would write one file per cluster into the user's folder: refused, nothing changed
        dnds = (proj / "dnds.tsv").read_bytes()
        with redirect_stderr(io.StringIO()) as err:
            res = run(keep_alignments=True)
        self.assertEqual(res.exit_code, 1)
        self.assertIn("was not written by a COAL run", err.getvalue())
        self.assertEqual(mine.read_text(encoding="utf-8"), ">my own alignment\nATG\n")
        self.assertEqual((proj / "dnds.tsv").read_bytes(), dnds)
        # in a folder of its own COAL writes them, knows them as its own, and removes them with the next run
        (proj / "codon_alignments").rename(proj / "my_alignments")
        res = run(keep_alignments=True)
        self.assertEqual(res.exit_code, 2)
        written = sorted(p.name for p in (proj / "codon_alignments").iterdir())
        self.assertTrue(written and all(n.endswith(".fasta") for n in written), written)
        (proj / "codon_alignments" / "notes.txt").write_text("kept\n", encoding="utf-8")
        self.assertEqual(run().exit_code, 2)
        self.assertEqual(sorted(p.name for p in (proj / "codon_alignments").iterdir()), ["notes.txt"])
        self.assertEqual(clu.read_bytes(), FAM_CLUSTERS.read_bytes())
        self.assertTrue((proj / "my_alignments" / "famA_001.fasta").is_file())

    @unittest.skipIf(os.name == "nt", "uses a POSIX sleep as the other process")
    def test_one_run_per_output_folder(self):
        q, s, pf = self.small()
        out = self.tmp / "out"
        lock = Path(str(out) + ".coal-lock")
        other = subprocess.Popen(["sleep", "60"])   # a live process that "owns" the folder
        try:
            out.mkdir()
            (out / "report.json").write_text("the other run's report\n", encoding="utf-8")
            lock.write_text(f"{other.pid}\n", encoding="utf-8")
            with redirect_stderr(io.StringIO()) as err:
                res = P.run_pairs(q, s, out, pairs=pf, quiet=True)
            self.assertEqual(res.exit_code, 1)
            self.assertIn(f"another COAL run (process {other.pid}) is writing", err.getvalue())
            self.assertIn("use another output folder", err.getvalue())
            self.assertEqual(res.report["status"], "failed")
            # nothing of the folder was touched, and the other run's lock is still there
            self.assertEqual(sorted(p.name for p in out.iterdir()), ["report.json"])
            self.assertEqual((out / "report.json").read_text(encoding="utf-8"), "the other run's report\n")
            self.assertEqual(lock.read_text(encoding="utf-8"), f"{other.pid}\n")
        finally:
            other.kill()
            other.wait()
        # the lock of a process that is gone is replaced, and a finished run leaves none
        res = P.run_pairs(q, s, out, pairs=pf, quiet=True)
        self.assertEqual(res.exit_code, 0)
        self.assertFalse(lock.exists())
        self.assertEqual(res.report["status"], "completed")
        with redirect_stderr(io.StringIO()):
            bad = P.run_pairs(self.tmp / "missing.fa", s, out, pairs=pf, quiet=True)
        self.assertEqual((bad.exit_code, bad.report["status"]), (1, "failed"))
        self.assertFalse(lock.exists(), "a failed run releases the folder too")

    @unittest.skipIf(os.name == "nt", "POSIX signals and a shell stand-in for the engine")
    def test_interrupt_ends_the_engine_and_is_reported(self):
        """SIGINT or SIGTERM sent to the Python process only (not its process group): the engine it started is
        killed (no orphan that goes on writing), report.json says `interrupted`, the lock is released."""
        q, s, pf = self.small()
        for signum, code, name in ((signal.SIGINT, 130, "SIGINT"), (signal.SIGTERM, 143, "SIGTERM")):
            with self.subTest(signal=name):
                out = self.tmp / f"out_{name}"
                pidfile = self.tmp / f"engine_{name}.pid"
                # the real coal-cpp for everything but the engine step, which stands still (as on a huge input)
                engine = self.tmp / f"coal-cpp-{name}"
                engine.write_text(f'#!/bin/sh\ncase "$1" in\n  codon-dnds) echo $$ > "{pidfile}"; '
                                  f'echo partial > "{out}/dnds.tsv.tmp"; exec sleep 120;;\n'
                                  f'  *) exec "{CPP}" "$@";;\nesac\n', encoding="utf-8")
                engine.chmod(0o755)
                proc = subprocess.Popen([sys.executable, "-m", "coal", "pipeline", "pairs", "--query", str(q),
                                         "--subject", str(s), "--pairs", str(pf), "-o", str(out)],
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                        env={**os.environ, "PYTHONPATH": str(ROOT / "src"), "COAL_CPP": str(engine)},
                                        encoding="utf-8", start_new_session=True)   # only the test can signal the child
                try:
                    for _ in range(600):
                        if pidfile.is_file() and pidfile.read_text(encoding="utf-8").strip() and (out / "dnds.tsv.tmp").exists():
                            break
                        time.sleep(0.05)
                    else:
                        self.fail("the engine step did not start")
                    engine_pid = int(pidfile.read_text(encoding="utf-8"))
                    os.kill(engine_pid, 0)            # it is running
                    time.sleep(0.2)
                    os.kill(proc.pid, signum)         # the Python process only
                    _, err = proc.communicate(timeout=60)
                finally:
                    if proc.poll() is None:
                        proc.kill()
                        proc.wait()
                self.assertEqual(proc.returncode, code, err)
                self.assertIn("interrupted", err)
                for _ in range(100):                  # the engine is gone (reaped by the pipeline before it exits)
                    try:
                        os.kill(engine_pid, 0)
                    except ProcessLookupError:
                        break
                    time.sleep(0.05)
                else:
                    os.kill(engine_pid, signal.SIGKILL)
                    self.fail("the engine process was left running")
                report = json.loads((out / "report.json").read_text(encoding="utf-8"))
                validate(report, SCHEMA)
                self.assertEqual((report["status"], report["exit_code"]), ("interrupted", code))
                self.assertIn(f"interrupted by {name}", report["error"])
                self.assertFalse((out / "dnds.tsv.tmp").exists())
                self.assertFalse((out / "dnds.tsv").exists())
                self.assertFalse(Path(str(out) + ".coal-lock").exists())

    def test_interrupted_msa_batch_ends_its_aligners(self):
        """coal.msa: an interrupt while a batch of clusters is aligned ends the running aligner processes and
        drops the clusters not started yet (the executor used to align all of them first)."""
        if os.name == "nt":
            self.skipTest("POSIX shell stand-in")
        started = self.tmp / "started"
        started.mkdir()
        tool = self.tmp / "clustalo"
        tool.write_text(f'#!/bin/sh\ntouch "{started}/$$"\nexec sleep 120\n', encoding="utf-8")
        tool.chmod(0o755)
        clusters = {k: [("a", "MKV"), ("b", "MKL"), ("c", "MRV")] for k in range(40)}

        def interrupt():
            while len(list(started.iterdir())) < 2:
                time.sleep(0.02)
            os.kill(os.getpid(), signal.SIGINT)   # a real Ctrl-C: wakes the main thread where it waits
        threading.Thread(target=interrupt, daemon=True).start()
        t0 = time.time()
        with self.assertRaises(KeyboardInterrupt):
            msa.align_clusters(clusters, tool="clustalo", jobs=2, exe=str(tool), pairs_with="msa")
            time.sleep(30)   # (the interrupt arrives here at the latest)
        self.assertLess(time.time() - t0, 20)
        pids = [int(p.name) for p in started.iterdir()]
        self.assertLess(len(pids), 10, "the clusters not started yet were dropped")
        time.sleep(0.3)
        for pid in pids:
            with self.assertRaises(ProcessLookupError, msg=f"aligner {pid} was left running"):
                os.kill(pid, 0)
                os.kill(pid, signal.SIGKILL)

    @unittest.skipIf(os.name == "nt", "the stand-in DIAMOND is a POSIX shell script")
    def test_a_file_without_coding_dna_is_explained_before_diamond(self):
        q, s, _ = self.small()
        marker = self.tmp / "diamond_was_started"
        fake = self.tmp / "diamond"
        fake.write_text(f'#!/bin/sh\ntouch "{marker}"\necho "diamond version 2.2.8"\n', encoding="utf-8")
        fake.chmod(0o755)
        prot = self.tmp / "proteins.fa"
        r = run_cpp("translate", q, "-o", prot, "-q")
        self.assertEqual(r.returncode, 0, r.stderr)
        with redirect_stderr(io.StringIO()) as err:
            res = P.run_pairs(prot, s, self.tmp / "o_prot", diamond=str(fake), quiet=True)
        self.assertEqual(res.exit_code, 1)
        msg = err.getvalue()
        self.assertIn("proteins.fa does not hold coding DNA: 12 of its 12 sequences", msg)
        self.assertIn("Is this a protein file?", msg)
        self.assertNotIn("Error detecting input file format", msg)
        self.assertFalse(marker.exists(), "DIAMOND was started")
        self.assertEqual([st["name"] for st in res.report["stages"]], ["translate_query", "translate_subject"])
        self.assertEqual(res.report["counts"]["genes"]["failed"], 12)
        self.assertEqual(len(read_tsv(self.tmp / "o_prot" / "failed.tsv")), 12)   # every gene with its reason
        # all genes unusable for other reasons: the subject file this time; and cluster mode
        empty = self.tmp / "empty.fa"
        empty.write_text(">a\n\n>b\nAT\n>c\n\n", encoding="utf-8")
        with redirect_stderr(io.StringIO()) as err:
            res = P.run_pairs(q, empty, self.tmp / "o_empty", diamond=str(fake), quiet=True)
        self.assertEqual(res.exit_code, 1)
        self.assertIn("none of the 3 sequences of", err.getvalue())
        self.assertIn("2 cds_empty, 1 cds_too_short", err.getvalue())
        with redirect_stderr(io.StringIO()) as err:
            res = P.run_clusters(prot, self.tmp / "o_cl", diamond=str(fake), quiet=True)
        self.assertEqual(res.exit_code, 1)
        self.assertIn("does not hold coding DNA", err.getvalue())
        self.assertFalse(marker.exists(), "DIAMOND was started")

    def test_a_non_ascii_id_is_that_genes_failure(self):
        q, s, pf = self.small()
        text = q.read_text(encoding="utf-8")
        for label, raw, shown in (("utf8", "gène_αβ".encode("utf-8"), "gène_αβ"),
                                  ("latin1", "gène".encode("latin-1"), "g\\xe8ne")):
            with self.subTest(encoding=label):
                bad = self.tmp / f"q_{label}.fa"
                bad.write_bytes(text.encode().replace(b">g03\n", b">" + raw + b"\n"))
                r = run_cpp("translate", bad, "--genes", self.tmp / "genes.tsv", "-o", self.tmp / "p.fa")
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertIn(b"12 gene(s): 11 ok, 0 warning, 1 failed", r.stderr)
                res = P.run_pairs(bad, s, self.tmp / f"o_{label}", pairs=pf, quiet=True)
                self.assertEqual(res.exit_code, 2)               # finished, with failures
                self.assertEqual(res.report["status"], "completed_with_failures")
                failed = read_tsv(self.tmp / f"o_{label}" / "failed.tsv")
                gene = [f for f in failed if f["level"] == "gene"]
                self.assertEqual([(f["seq1"], f["code"]) for f in gene], [(shown, "id_not_ascii")])
                self.assertIn("non-ASCII character in the sequence ID (line 5)", gene[0]["message"])
                rows = read_tsv(res.outputs["dnds"])
                self.assertEqual(sum(r["status"] == "ok" for r in rows), 11)   # the other genes are computed
                self.assertEqual([r["status"] for r in rows if r["seq1"] == "g03"], ["error"])


@needs_cpp
class Encoding(PipelineCase):
    """IDs are bytes that are never re-encoded, and every table is UTF-8 with \\n line ends whatever the system's
    default encoding and line ends are (cp1252 and \\r\\n on Windows). The Windows condition is made here on every
    system: the pipeline runs in its own interpreter with UTF-8 mode off and the C locale (default encoding ASCII
    on POSIX, the ANSI code page on Windows), and with -X warn_default_encoding as an error, which fails any open()
    or text-mode pipe that does not name its encoding."""

    LATIN1 = "gène".encode("latin-1")            # not valid UTF-8: shown as g\xe8ne
    UTF8 = "gène_αβ".encode("utf-8")             # valid UTF-8: shown as it is
    SHOWN = {"g03": "gène_αβ", "g05": "g\\xe8ne", "g07": "g\\xe2ne"}

    small = Robustness.small

    def inputs(self):
        q, s, _ = self.small()
        raw = q.read_text(encoding="utf-8").encode()   # (\n line ends, whatever write_text made of them)
        for old, new in ((b"g03", self.UTF8), (b"g05", self.LATIN1), (b"g07", "gâne".encode("latin-1"))):
            self.assertIn(b">" + old + b"\n", raw)
            raw = raw.replace(b">" + old + b"\n", b">" + new + b"\n")
        bad = self.tmp / "q_ids.fa"
        bad.write_bytes(raw)
        # the pairs file as a Windows editor saves it (\r\n), naming the genes by the bytes of the FASTA file, and
        # one ID (Latin-1) that is in no file
        pf = self.tmp / "pairs_ids.tsv"
        names = {b"g03": self.UTF8, b"g05": self.LATIN1, b"g07": "gâne".encode("latin-1")}
        lines = [names.get(a.encode(), a.encode()) + b"\t" + b.encode() for a, b in PLAIN]
        lines.append("señor".encode("latin-1") + b"\th01")
        pf.write_bytes(b"query\tsubject\r\n" + b"\r\n".join(lines) + b"\r\n")
        return bad, s, pf

    def run_cli(self, out, bad, s, pf, env, *python_opts):
        return subprocess.run([sys.executable, *python_opts, "-m", "coal", "pipeline", "pairs", "--query", str(bad),
                               "--subject", str(s), "--pairs", str(pf), "-o", str(out), "--kaks-format",
                               "--keep-alignments"], capture_output=True, env=env)

    def tables(self, out: Path) -> dict[str, bytes]:
        return {p.name: p.read_bytes() for p in sorted(out.iterdir()) if p.is_file() and p.name != "report.json"}

    def test_the_same_bytes_whatever_the_default_encoding(self):
        bad, s, pf = self.inputs()
        base = {k: v for k, v in os.environ.items()
                if k not in ("PYTHONUTF8", "PYTHONIOENCODING", "LC_ALL", "LC_CTYPE", "LANG", "PYTHONWARNINGS",
                             "PYTHONWARNDEFAULTENCODING")}
        base.update(PYTHONPATH=str(ROOT / "src"), COAL_CPP=str(CPP))
        legacy = dict(base, PYTHONUTF8="0", PYTHONCOERCECLOCALE="0", LC_ALL="C", LANG="C")
        runs = {
            "utf8": (dict(base, PYTHONUTF8="1"), ()),
            "legacy": (legacy, ()),                     # ASCII (POSIX) or the ANSI code page (Windows) by default
            "legacy_strict": (legacy, ("-X", "warn_default_encoding", "-W", "error::EncodingWarning")),
            "cp1252_stdio": (dict(legacy, PYTHONIOENCODING="cp1252"), ()),   # a Windows console's streams
        }
        got = {}
        for name, (env, opts) in runs.items():
            with self.subTest(run=name):
                out = self.tmp / f"o_{name}"
                r = self.run_cli(out, bad, s, pf, env, *opts)
                self.assertEqual(r.returncode, 2, r.stderr.decode("utf-8", "backslashreplace"))
                self.assertNotIn(b"Traceback", r.stderr)
                self.assertNotIn(b"EncodingWarning", r.stderr)
                got[name] = self.tables(out)
                report = json.loads((out / "report.json").read_bytes().decode("utf-8"))
                self.assertEqual(report["status"], "completed_with_failures")
                seq1 = {f.get("seq1") for f in report["failures"]}
                self.assertLessEqual(set(self.SHOWN.values()) | {"se\\xf1or"}, seq1)
        ref = got["utf8"]
        self.assertLessEqual({"dnds.tsv", "genes.tsv", "failed.tsv", "orthologs.tsv", "dnds.kaks",
                              "codon_alignments.fasta"}, set(ref))
        for name, files in got.items():
            self.assertEqual(sorted(files), sorted(ref), name)
            for fn, data in files.items():
                self.assertEqual(data, ref[fn], f"{fn} of the {name} run differs from the UTF-8 run")
        for fn, data in ref.items():
            with self.subTest(file=fn):
                self.assertNotIn(b"\r", data)              # \n line ends on every system
                data.decode("utf-8")                       # valid UTF-8, always (strict)
                self.assertNotIn("\ufffd".encode(), data)  # and nothing replaced by a look-alike

        # the IDs: valid UTF-8 unchanged, other bytes as \xhh; two Latin-1 IDs that differ stay two IDs
        failed = [ln.split(b"\t") for ln in ref["failed.tsv"].split(b"\n")[1:] if ln]
        genes = {f[2]: f[4] for f in failed if f[0] == b"gene"}
        for shown in self.SHOWN.values():
            self.assertEqual(genes.get(shown.encode("utf-8")), b"id_not_ascii", shown)
        self.assertIn(self.UTF8, ref["genes.tsv"])
        self.assertIn(b"g\\xe8ne\t", ref["genes.tsv"])
        self.assertIn(b"g\\xe2ne\t", ref["genes.tsv"])
        self.assertIn(b"se\\xf1or\th01\n", ref["orthologs.tsv"])   # the pairs file's ID, by the same rule
        pair_rows = [f for f in failed if f[0] == b"pair"]
        self.assertEqual(sorted(f[2] for f in pair_rows),
                         sorted(x.encode("utf-8") for x in (*self.SHOWN.values(), "se\\xf1or")))
        rows = [ln.split(b"\t") for ln in ref["dnds.tsv"].split(b"\n")[1:] if ln]
        self.assertEqual(sum(f[-2] == b"ok" for f in rows), 9)      # the nine other genes are computed

    def test_messages_with_such_ids_reach_any_console(self):
        # every gene unusable, the first with a non-ASCII ID: the message names it; an ASCII or cp1252 stderr
        # shows what it cannot encode as an escape and the run ends with its own error, not a UnicodeEncodeError
        _, s, _ = self.small()
        q = self.tmp / "q_bad.fa"
        q.write_bytes(b">" + self.UTF8 + b"\nAT\n>b\n\n")
        pf = self.tmp / "p.tsv"
        pf.write_bytes(self.UTF8 + b"\th01\n")
        base = {k: v for k, v in os.environ.items() if k not in ("PYTHONUTF8", "PYTHONIOENCODING", "LC_ALL", "LANG")}
        base.update(PYTHONPATH=str(ROOT / "src"), COAL_CPP=str(CPP), PYTHONUTF8="0", PYTHONCOERCECLOCALE="0",
                    LC_ALL="C", LANG="C")
        for name, env in (("default", base), ("cp1252", dict(base, PYTHONIOENCODING="cp1252")),
                          ("ascii", dict(base, PYTHONIOENCODING="ascii"))):
            with self.subTest(stderr=name):
                r = self.run_cli(self.tmp / f"o_{name}", q, s, pf, env)
                self.assertEqual(r.returncode, 1, r.stderr)
                self.assertNotIn(b"Traceback", r.stderr)
                self.assertIn(b"none of the 2 sequences of", r.stderr)
                self.assertIn(b"the first: g", r.stderr)

    def test_no_text_file_is_opened_with_the_default_encoding(self):
        # the modules of the alignment layer, the pipeline and the dispatcher name the encoding of every text file
        # and text-mode pipe they open (what -X warn_default_encoding checks at run time, checked here in the source)
        import ast
        for mod in ("dnds_pipeline.py", "align.py", "msa.py", "dispatch.py"):
            tree = ast.parse((ROOT / "src" / "coal" / mod).read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                kw = {k.arg: k.value for k in node.keywords}
                name = getattr(node.func, "attr", getattr(node.func, "id", ""))
                where = f"src/coal/{mod}:{node.lineno}"
                if name == "open" and isinstance(node.func, ast.Name):
                    mode = node.args[1].value if len(node.args) > 1 else "r"
                    if "b" in mode:
                        continue
                    self.assertIn("encoding", kw, where)
                    self.assertEqual(kw["encoding"].value, "utf-8", where)
                    self.assertEqual(getattr(kw.get("errors"), "value", None), "backslashreplace", where)
                    if "w" in mode or "a" in mode:
                        self.assertEqual(getattr(kw.get("newline"), "value", None), "\n", where)
                elif name in ("read_text", "write_text"):
                    self.assertIn("encoding", kw, where)
                elif getattr(kw.get("text"), "value", False) or "universal_newlines" in kw:
                    self.assertEqual(getattr(kw.get("encoding"), "value", None), "utf-8", where)
                    self.assertEqual(getattr(kw.get("errors"), "value", None), "backslashreplace", where)
                elif name == "decode":
                    self.assertEqual([getattr(a, "value", None) for a in node.args] if node.args else ["utf-8", "strict"],
                                     ["utf-8", "backslashreplace"] if node.args else ["utf-8", "strict"], where)


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
_STANDARD = dict(zip((a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"),
                     "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


@needs_cpp
class ProteinAlignments(PipelineCase):
    """--keep-protein-alignments: the protein alignments the codon alignments were threaded on, in the layout of
    --keep-alignments (pairs: protein_alignments.fasta; clusters: protein_alignments/<cluster>.fasta)."""

    def inputs(self):
        q, s = self.tmp / "q.fa", self.tmp / "s.fa"
        allq, alls = read_fasta(Q), read_fasta(S)
        q.write_text("".join(f">{a}\n{allq[a]}\n" for a, _ in PLAIN), encoding="utf-8")
        s.write_text("".join(f">{b}\n{alls[b]}\n" for _, b in PLAIN), encoding="utf-8")
        return q, s, self.pairs_file(PLAIN)

    def assert_threaded(self, protein: Path, codon: Path) -> int:
        """Record k of the protein file is record k of the codon file, column for column: a gap against `---`,
        a residue against the codon that codes for it. Returns the number of records."""
        prot = Path(protein).read_text(encoding="utf-8").splitlines()
        cod = Path(codon).read_text(encoding="utf-8").splitlines()
        self.assertEqual(prot[0::2], cod[0::2])
        self.assertTrue(prot and all(h.startswith(">") for h in prot[0::2]))
        for name, p, c in zip(prot[0::2], prot[1::2], cod[1::2]):
            self.assertEqual(len(c), 3 * len(p), name)
            for k, aa in enumerate(p):
                codon_ = c[3 * k:3 * k + 3]
                self.assertEqual(aa == "-", codon_ == "---", (name, k))
                if aa != "-":
                    self.assertEqual(_STANDARD.get(codon_.upper().replace("U", "T"), "X"), aa, (name, k))
        return len(prot) // 2

    def test_pairs(self):
        q, s, pf = self.inputs()
        out = self.tmp / "out"
        res = P.run_pairs(q, s, out, pairs=pf, quiet=True, keep_alignments=True, keep_protein_alignments=True)
        self.assertEqual(res.exit_code, 0)
        self.check_report(res)
        self.assertEqual(res.outputs["protein_alignments"], str(out / "protein_alignments.fasta"))
        self.assertEqual(res.report["outputs"]["protein_alignments"], "protein_alignments.fasta")
        self.assertIs(res.report["parameters"]["keep_protein_alignments"], True)
        n = self.assert_threaded(out / "protein_alignments.fasta", out / "codon_alignments.fasta")
        self.assertEqual(n, 2 * len(PLAIN))
        names = (out / "protein_alignments.fasta").read_text(encoding="utf-8").splitlines()[0::2]
        self.assertEqual([x[1:] for x in names], [g for pair in PLAIN for g in pair])
        # the rows are the alignment the table's prot_aln_length counts
        rows = read_tsv(res.outputs["dnds"])
        lens = [len(x) for x in (out / "protein_alignments.fasta").read_text(encoding="utf-8").splitlines()[1::4]]
        self.assertEqual([int(r["prot_aln_length"]) for r in rows], lens)
        rec = json.loads((out / "work" / "outputs.json").read_text(encoding="utf-8"))
        self.assertIn("protein_alignments.fasta", rec["files"])
        self.assertIs(rec["protein_alignments_dir"], False)
        # alone (no codon alignments), from the command line; then a run without it removes the file
        first = (out / "protein_alignments.fasta").read_bytes()
        with redirect_stderr(io.StringIO()):
            rc = P.main(["pairs", "--query", str(q), "--subject", str(s), "--pairs", str(pf), "-o", str(out), "-q",
                         "--keep-protein-alignments"])
        self.assertEqual(rc, 0)
        self.assertEqual((out / "protein_alignments.fasta").read_bytes(), first)
        self.assertFalse((out / "codon_alignments.fasta").exists())
        res = P.run_pairs(q, s, out, pairs=pf, quiet=True)
        self.assertEqual(res.exit_code, 0)
        self.assertFalse((out / "protein_alignments.fasta").exists())
        self.assertNotIn("protein_alignments", res.report["outputs"])

    def test_pairs_with_compared_terminal_stops(self):
        """--terminal-stops compare (orthologr's reading) changes what Ka/Ks compares, not the alignment files:
        they end before the terminal stop codons either way."""
        q, s = self.tmp / "q.fa", self.tmp / "s.fa"
        q.write_text(">a\nATGGCTAAAGGTCTGTTTGATCCGGAATGGTAA\n", encoding="utf-8")
        s.write_text(">b\nATGGCAAAGGGTTTGTTCCCGGAGTGGTGA\n", encoding="utf-8")
        for stops, last in (("drop", "W"), ("compare", "W")):
            out = self.tmp / stops
            res = P.run_pairs(q, s, out, pairs=self.pairs_file([("a", "b")]), quiet=True, method=["Comeron"],
                              comeron_flavour="gestimator", terminal_stops=stops, keep_alignments=True,
                              keep_protein_alignments=True)
            self.assertEqual(res.exit_code, 0, stops)
            self.assertEqual(self.assert_threaded(out / "protein_alignments.fasta", out / "codon_alignments.fasta"), 2)
            rows = (out / "protein_alignments.fasta").read_text(encoding="utf-8").splitlines()
            self.assertEqual((rows[1][-1], rows[3][-1]), (last, last), stops)
            self.assertIn("-", rows[3])

    def test_clusters(self):
        """Clusters aligned as a whole get one file each (here: clusters of two members, aligned by the engine);
        the folder is recorded in work/outputs.json, so the next run removes its files; a folder of that name
        that no COAL run wrote is not written into."""
        out = self.tmp / "out"
        clu = self.tmp / "clusters.tsv"
        ids = list(read_fasta(FAM))
        clu.write_text(f"{ids[0]}\t{ids[0]}\n{ids[0]}\t{ids[1]}\n{ids[2]}\t{ids[2]}\n{ids[2]}\t{ids[3]}\n",
                       encoding="utf-8")

        def run(**kw):
            return P.run_clusters(FAM, out, clusters=clu, quiet=True, **kw)
        res = run(keep_alignments=True, keep_protein_alignments=True)
        self.assertIn(res.exit_code, (0, 2))
        self.check_report(res)
        self.assertEqual(res.outputs["protein_alignments"], str(out / "protein_alignments"))
        files = sorted(os.listdir(out / "protein_alignments"))
        self.assertEqual(files, sorted(os.listdir(out / "codon_alignments")))
        self.assertEqual(len(files), 2)
        for f in files:
            self.assertEqual(self.assert_threaded(out / "protein_alignments" / f, out / "codon_alignments" / f), 2)
        rec = json.loads((out / "work" / "outputs.json").read_text(encoding="utf-8"))
        self.assertIs(rec["protein_alignments_dir"], True)
        (out / "protein_alignments" / "notes.txt").write_text("kept\n", encoding="utf-8")
        res = run()
        self.assertIn(res.exit_code, (0, 2))
        self.assertEqual(os.listdir(out / "protein_alignments"), ["notes.txt"])    # COAL's files gone, the user's kept
        self.assertFalse((out / "codon_alignments").exists())
        dnds = (out / "dnds.tsv").read_bytes()
        with redirect_stderr(io.StringIO()) as err:
            res = run(keep_protein_alignments=True)                                 # not COAL's folder any more
        self.assertEqual(res.exit_code, 1)
        self.assertIn("protein_alignments exists and was not written by a COAL run: --keep-protein-alignments",
                      err.getvalue())
        self.assertEqual((out / "dnds.tsv").read_bytes(), dnds)
        self.assertEqual(os.listdir(out / "protein_alignments"), ["notes.txt"])
        # an input inside the folder the run writes is refused before anything changes
        (out / "protein_alignments" / "cds.fasta").write_bytes(FAM.read_bytes())
        with redirect_stderr(io.StringIO()) as err:
            res = P.run_clusters(out / "protein_alignments" / "cds.fasta", out, clusters=clu, quiet=True)
        self.assertEqual(res.exit_code, 1)
        self.assertIn("a folder this run writes", err.getvalue())

    def test_engine_options(self):
        """`coal-cpp codon-dnds --protein-out / --protein-dir`: the same groups as --codon-out / --codon-dir; a
        group marked nocodon and a group that failed have no protein alignment either."""
        groups = (b"G\tg1\tall\nM\ta\tATGGCTAAAGGTCTGTTTGATCCGGAATGG\nM\tb\tATGGCAAAGGGTTTGTTCCCGGAGTGG\nE\n"
                  b"G\tg2\tall\tnocodon\nM\tc\tATGGCTAAAGGTCTG\nM\td\tATGGCAAAGGGTTTG\nE\n"
                  b"G\tg3\tall\nM\te\tATGGCTAAAGGTCTG\tMAKGL\nM\tf\tATGGCAAAGGGT\tMAKG-\n"
                  b"M\tg\tATGGCAAAGGGTTTG\tMAKGL\nE\n")
        cd, pd = self.tmp / "cd", self.tmp / "pd"
        cd.mkdir()
        pd.mkdir()
        r = run_cpp("codon-dnds", "-o", self.tmp / "o.tsv", "--codon-out", self.tmp / "c.fa", "--protein-out",
                    self.tmp / "p.fa", "--codon-dir", cd, "--protein-dir", pd, "-q", stdin=groups)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.assert_threaded(self.tmp / "p.fa", self.tmp / "c.fa"), 5)
        self.assertEqual(sorted(os.listdir(pd)), ["g1.fasta", "g3.fasta"])
        self.assertEqual(sorted(os.listdir(cd)), ["g1.fasta", "g3.fasta"])
        self.assertEqual((pd / "g3.fasta").read_text(encoding="utf-8"), ">e\nMAKGL\n>f\nMAKG-\n>g\nMAKGL\n")
        self.assertEqual(self.assert_threaded(pd / "g1.fasta", cd / "g1.fasta"), 2)
        r = run_cpp("codon-dnds", "-o", self.tmp / "o2.tsv", "--protein-out", self.tmp / "p2.fa", "-q", stdin=groups)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual((self.tmp / "p2.fa").read_bytes(), (self.tmp / "p.fa").read_bytes())


def recompute_summary(rows: list[dict]) -> dict:
    """{(cluster, method): expected summary cells} from dnds.tsv rows, with the statistics module."""
    out = {}
    groups: dict = {}
    for r in rows:
        groups.setdefault((r["cluster_id"], r["method"]), []).append(r)
    for key, g in groups.items():
        ok = [r for r in g if r["status"] in ("ok", "warning")]
        e = {"n_pairs": len(g), "n_ok": len(ok), "n_saturated": sum(r["status"] == "saturated" for r in g),
             "n_not_converged": sum(r["status"] == "not_converged" for r in g),
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


class SummaryCounts(unittest.TestCase):
    """cluster_summary.tsv and the report's pair counts with a pair whose iteration did not converge: it has no
    estimate, so it is not in n_ok and not in the statistics; it has a column and a count of its own."""

    HEAD = ("cluster_id",) + tuple(P.PAIR_COLUMNS)

    def row(self, cid, a, b, status, dn="NA", ds="NA"):
        w = "NA" if "NA" in (dn, ds) else repr(float(dn) / float(ds))
        return [cid, a, b, "YN", dn, ds, w, "NA", "NA", "2", "1", "NA", "NA", "1", "1", "0.00", status, ""]

    def test_not_converged_is_counted_apart(self):
        with tempfile.TemporaryDirectory() as td:
            d = os.path.join(td, "dnds.tsv")
            rows = [self.row("c1", "a", "b", "ok", "0.1", "0.5"), self.row("c1", "a", "c", "not_converged"),
                    self.row("c1", "b", "c", "warning", "0.2", "0.5"), self.row("c1", "a", "d", "saturated", "0.3"),
                    self.row("c2", "x", "y", "not_converged")]
            with open(d, "w", encoding="utf-8", newline="\n") as fh:
                fh.write("\t".join(self.HEAD) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
            plans = [{"cluster_id": c, "representative": r, "n_members": n, "n_used": n, "pair_selection": "all",
                      "aligner": "pairwise", "note": ""} for c, r, n in (("c1", "a", 4), ("c2", "x", 2), ("c3", "z", 1))]
            out = os.path.join(td, "s.tsv")
            P.summarise(d, plans, out)
            summ = read_tsv(out)
            self.assertEqual(list(summ[0]), list(P.SUMMARY_COLUMNS))
            self.assertIn("n_not_converged", P.SUMMARY_COLUMNS)
            c1, c2, c3 = summ
            self.assertEqual([c1[k] for k in ("n_pairs", "n_ok", "n_not_converged", "n_saturated", "n_no_data",
                                              "n_failed", "n_dNdS", "status")], ["4", "2", "1", "1", "0", "0", "2", "ok"])
            self.assertAlmostEqual(float(c1["dN_mean"]), 0.15)
            # a cluster whose only pair did not converge has no pair with an estimate
            self.assertEqual([c2[k] for k in ("n_pairs", "n_ok", "n_not_converged", "dN_median", "status")],
                             ["1", "0", "1", "NA", "failed"])
            self.assertEqual([c3[k] for k in ("n_pairs", "n_ok", "n_not_converged", "n_saturated", "status")],
                             ["0", "0", "0", "0", "ok"])
            self.assertTrue(all(len(x) == len(P.SUMMARY_COLUMNS) for x in summ))

            class Run:
                counts: dict = {}
                def fail(self, *a, **k):
                    raise AssertionError("no pair failed")
            run = Run()
            P._count_pairs(run, d, True)
            self.assertEqual(run.counts["pairs"], {"total": 5, "ok": 1, "warning": 1, "not_converged": 2,
                                                   "saturated": 1, "no_data": 0, "error": 0})


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
        one.write_text("".join(f">{n}\n{fam[n]}\n" for n in names[:4]), encoding="utf-8")
        two.write_text("".join(f">{n}\n{fam[n]}\n" for n in names[4:]), encoding="utf-8")
        cl = self.tmp / "cl.tsv"
        cl.write_text("".join(ln + "\n" for ln in FAM_CLUSTERS.read_text(encoding="utf-8").splitlines() if not ln.startswith("famB")), encoding="utf-8")
        res = P.run_clusters([one, two], self.tmp / "out", clusters=cl, max_cluster_size=2, keep_alignments=True,
                             kaks_format=True, quiet=True)
        self.assertEqual(res.exit_code, 0)
        self.assertEqual(res.report["counts"]["genes"]["total"], 9)
        self.assertEqual(sorted(os.listdir(res.outputs["codon_alignments"])), ["famC_001.fasta"])
        self.assertEqual(len(Path(res.outputs["kaks"]).read_text(encoding="utf-8").splitlines()), 1 + 5 + 1)

    @unittest.skipIf(os.name == "nt", "the stand-in tool is a POSIX shell script")
    def test_a_failing_msa_falls_back_to_pairwise(self):
        tool = self.tmp / "clustalo"
        tool.write_text("#!/bin/sh\necho 'cannot align' >&2\nexit 1\n", encoding="utf-8")
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
                             keep_alignments=True, keep_protein_alignments=True, quiet=True)
        self.assertEqual(res.exit_code, 2)
        # the MSA tool's protein rows, kept under the same names as the codon MSAs and threaded column by column
        self.assertEqual(sorted(os.listdir(res.outputs["protein_alignments"])), ["famA_001.fasta", "famC_001.fasta"])
        n = ProteinAlignments.assert_threaded(self, Path(res.outputs["protein_alignments"]) / "famA_001.fasta",
                                              Path(res.outputs["codon_alignments"]) / "famA_001.fasta")
        self.assertEqual(n, 6)
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
        realign = Path(res.outputs["realign"]).read_text(encoding="utf-8").splitlines()
        self.assertEqual(realign[0].split("\t")[:3], ["member", "representative", "approx_pident"])
        self.assertGreaterEqual(len(realign), 1 + 16)
        again = P.run_clusters(FAM, self.tmp / "out", max_cluster_size=2, realign=True, jobs=1, quiet=True,
                               resume=True)
        st = {s["name"]: s["cached"] for s in again.report["stages"]}
        self.assertTrue(st["deepclust"] and st["realign"] and st["translate_proteins"])


if __name__ == "__main__":
    unittest.main()
