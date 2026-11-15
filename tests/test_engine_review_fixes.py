"""Engine fixes of the review of 2026-11-15 (branch fix/engine-review), each checked in both engines where both
have the code path:

- a gzip-compressed input of the codon alignment (`coal PEP DNA`): "X is gzip-compressed: decompress it first",
  not "control character" / "sequence data before the first '>' header"; coal-cpp says it byte for byte the same
  (stderr, exit code, the report's error code), for the protein alignment, any DNA file and standard input;
- the WSL hint of a missing Windows path: one '/' between the parts and none at the end ('C:\\' is /mnt/c,
  'C:\\\\double\\\\x.fa' is /mnt/c/double/x.fa), the same in both engines;
- `coal-cpp codon-dnds` (the native step of `coal pipeline`) accepts the testing flags --kaks-baseline and
  --kaks-isa of `dnds` / `kaks`, without showing them in its help, and --kaks-baseline gives the same bytes;
- `coal-cpp kaks`, MYN with kappa_TC (or kappa_AG) estimated as infinite: the Substitution-Rate-Ratio field is NA
  (KaKs_Calculator prints `inf`; COAL's output has no nan / inf), the rest of the row and the dnds table unchanged;
- MLWL at kappa exactly 2 or 0.5 (2 ts / tv, ts and tv sums of thirds and halves of a difference): the site formula
  is chosen in exact arithmetic, so permuting the codon columns no longer moves a pair into the other formula
  (property fuzz P6, seed 4242: FIXCI-B2; KaKs_Calculator has the same order dependence);
- MLWL dS above 10 without a rounding residue (property fuzz P10, seed 4243: FIXCI-B3): a genuine value of the
  method, 5 times LWL's correction at kappa = 1 when there is no 4-fold site; KNOWN M01 in test_kaks_properties.py.

The Python halves run in-process; the C++ halves skip if coal-cpp is not built (COAL_CPP_BIN, default
./build/coal-cpp).
"""
from __future__ import annotations

import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal.seqio import GZIP_INPUT, InputFormatError, read_alignment, read_dna, windows_path_hint  # noqa: E402
from test_engine_ux import D, P, make, run_py  # noqa: E402
from test_differential import cpp_binary  # noqa: E402

GZ_P, GZ_D = gzip.compress(P.encode(), mtime=0), gzip.compress(D.encode(), mtime=0)
# name -> (files, arguments, the error line both engines print)
GZIP_CASES = {
    "protein": ({"p.fa.gz": GZ_P, "d.fa": D}, ["p.fa.gz", "d.fa"], f"coal: error: p.fa.gz {GZIP_INPUT}\n"),
    "dna": ({"p.fa": P, "d.fa.gz": GZ_D}, ["p.fa", "d.fa.gz"], f"coal: error: d.fa.gz {GZIP_INPUT}\n"),
    "second_dna": ({"p.fa": P, "d.fa": D, "e.fa.gz": GZ_D}, ["p.fa", "d.fa", "e.fa.gz"],
                   f"coal: error: e.fa.gz {GZIP_INPUT}\n"),
    "both": ({"p.fa.gz": GZ_P, "d.fa.gz": GZ_D}, ["p.fa.gz", "d.fa.gz"], f"coal: error: p.fa.gz {GZIP_INPUT}\n"),
    # found by its first bytes, whatever its name
    "no_suffix": ({"p.fa": P, "dna": GZ_D}, ["p.fa", "dna"], f"coal: error: dna {GZIP_INPUT}\n"),
    "clustal_output": ({"p.fa": P, "d.fa.gz": GZ_D}, ["p.fa", "d.fa.gz", "-output", "clustal"],
                       f"coal: error: d.fa.gz {GZIP_INPUT}\n"),
    # the errors before reading are unchanged: a missing file is named first
    "missing_first": ({"d.fa.gz": GZ_D}, ["nope.fa", "d.fa.gz"], "coal: error: cannot open 'nope.fa': no such file\n"),
}
# Missing Windows paths (U17): argument -> its WSL form in the hint
WSL = {
    "C:\\": "/mnt/c",
    "C:/": "/mnt/c",
    "C:\\\\double\\\\x.fa": "/mnt/c/double/x.fa",
    "D:/data//x y/n.fa": "/mnt/d/data/x y/n.fa",
    "e:\\data\\": "/mnt/e/data",
    "C:\\Users\\me\\p.fa": "/mnt/c/Users/me/p.fa",   # (unchanged)
    "C:\\/a\\/b": "/mnt/c/a/b",
}


def hint_line(arg: str, wsl: str) -> str:
    return (f"coal: error: cannot open '{arg}': no such file. This is a Windows path; under WSL (Windows Subsystem "
            f"for Linux) the same file is {wsl}\n")


class TestGzipInputPython(unittest.TestCase):
    def setUp(self):
        self.d = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def test_messages(self):
        for name, (files, args, line) in GZIP_CASES.items():
            with self.subTest(case=name):
                d = self.d / name
                d.mkdir()
                make(d, files)
                rc, out, err = run_py(d, args)
                self.assertEqual((rc, out, err), (1, "", line))
                self.assertNotIn("control character", err)

    def test_message_wording(self):
        self.assertTrue(GZIP_INPUT.startswith("is gzip-compressed: decompress it first"))

    def test_report_error_code(self):
        make(self.d, {"p.fa": P, "d.fa.gz": GZ_D})
        rc, _, _ = run_py(self.d, ["p.fa", "d.fa.gz", "-o", "o.fa", "--report", "r.json"])
        self.assertEqual(rc, 1)
        r = json.loads((self.d / "r.json").read_text())
        self.assertEqual(r["error"]["code"], "input_format")
        self.assertEqual(r["error"]["text"], f"d.fa.gz {GZIP_INPUT}")

    def test_readers(self):
        make(self.d, {"p.fa.gz": GZ_P, "d.fa.gz": GZ_D, "p.fa": P, "d.fa": D})
        for fn, path in ((read_alignment, "p.fa.gz"), (lambda p: read_dna([p]), "d.fa.gz")):
            with self.assertRaises(InputFormatError) as cm:
                fn(str(self.d / path))
            self.assertEqual(str(cm.exception), f"{self.d / path} {GZIP_INPUT}")
        # plain files are read as before
        self.assertEqual(read_alignment(str(self.d / "p.fa")), [("a", "MKV-L"), ("b", "MKVAL")])
        # a file that only starts with one of the two magic bytes is not gzip: the reader's usual message
        (self.d / "x.fa").write_bytes(b"\x1f>a\nMKV\n")
        with self.assertRaises(InputFormatError) as cm:
            read_alignment(str(self.d / "x.fa"))
        self.assertIn("control character", str(cm.exception))

    def test_standard_input(self):
        make(self.d, {"p.fa": P})
        env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), COAL_ENGINE="python")
        r = subprocess.run([sys.executable, "-m", "coal", "p.fa", "-"], cwd=self.d, input=GZ_D, capture_output=True,
                           env=env, timeout=120)
        err = r.stderr.decode().replace("\r\n", "\n")   # (Windows: CRLF)
        self.assertEqual((r.returncode, err), (1, f"coal: error: standard input {GZIP_INPUT}\n"))


@unittest.skipIf(os.name == "nt", "the WSL hint is given on Linux and macOS only")
class TestWslHintPython(unittest.TestCase):
    def test_hint(self):
        for arg, wsl in WSL.items():
            with self.subTest(path=arg):
                self.assertTrue(windows_path_hint(arg).endswith(f" the same file is {wsl}"), windows_path_hint(arg))

    def test_command(self):
        with tempfile.TemporaryDirectory() as tmp:
            make(Path(tmp), {"d.fa": D})
            for arg, wsl in WSL.items():
                with self.subTest(path=arg):
                    rc, _, err = run_py(Path(tmp), [arg, "d.fa"])
                    self.assertEqual((rc, err), (1, hint_line(arg, wsl)))


@unittest.skipUnless(cpp_binary(), "coal-cpp not built (COAL_CPP_BIN)")
class TestCppSaysTheSame(unittest.TestCase):
    def setUp(self):
        self.bin = str(cpp_binary())
        self.d = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def cpp(self, d: Path, args, stdin: bytes | None = None):
        r = subprocess.run([self.bin, *args], cwd=d, input=stdin, capture_output=True, timeout=120)
        return r.returncode, r.stdout.decode(), r.stderr.decode("utf-8", "replace").replace("\r\n", "\n")

    def test_gzip_messages(self):
        for name, (files, args, line) in GZIP_CASES.items():
            with self.subTest(case=name):
                d = self.d / name
                d.mkdir()
                make(d, files)
                py = run_py(d, args)
                self.assertEqual(self.cpp(d, args), py)
                self.assertEqual(py[2], line)

    def test_gzip_report(self):
        make(self.d, {"p.fa.gz": GZ_P, "d.fa": D})
        args = ["p.fa.gz", "d.fa", "-o", "o.fa", "--report", "r.json"]
        run_py(self.d, args)
        py = json.loads((self.d / "r.json").read_text())["error"]
        self.cpp(self.d, args)
        self.assertEqual(json.loads((self.d / "r.json").read_text())["error"], py)

    def test_gzip_standard_input(self):
        make(self.d, {"p.fa": P})
        self.assertEqual(self.cpp(self.d, ["p.fa", "-"], stdin=GZ_D),
                         (1, "", f"coal: error: standard input {GZIP_INPUT}\n"))
        self.assertEqual(self.cpp(self.d, ["-", "p.fa"], stdin=GZ_P),
                         (1, "", f"coal: error: standard input {GZIP_INPUT}\n"))

    def test_gzip_through_the_python_command(self):
        # `coal` handing a large or COAL_ENGINE=cpp run to coal-cpp: the same message
        make(self.d, {"p.fa": P, "d.fa.gz": GZ_D})
        env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), COAL_ENGINE="cpp", COAL_CPP=self.bin)
        r = subprocess.run([sys.executable, "-m", "coal", "p.fa", "d.fa.gz"], cwd=self.d, capture_output=True,
                           env=env, timeout=120)
        err = r.stderr.decode().replace("\r\n", "\n")   # (Windows: CRLF)
        self.assertEqual((r.returncode, err), (1, f"coal: error: d.fa.gz {GZIP_INPUT}\n"))

    def test_kaks_commands_keep_their_message(self):
        # coal-cpp kaks / dnds refuse a gzip input with their own words (they name the `coal` command, which
        # decompresses it): unchanged
        (self.d / "x.axt.gz").write_bytes(gzip.compress(b"p\nATG\nATG\n", mtime=0))
        r = subprocess.run([self.bin, "dnds", "x.axt.gz", "-o", "o.tsv"], cwd=self.d, capture_output=True)
        self.assertIn(b"'x.axt.gz' is gzip-compressed: decompress it first (gunzip -k 'x.axt.gz'). coal-cpp reads "
                      b"plain text only", r.stderr)

    @unittest.skipIf(os.name == "nt", "the WSL hint is given on Linux and macOS only")
    def test_wsl_hint(self):
        make(self.d, {"d.fa": D})
        for arg, wsl in WSL.items():
            with self.subTest(path=arg):
                self.assertEqual(self.cpp(self.d, [arg, "d.fa"]), (1, "", hint_line(arg, wsl)))
                self.assertEqual(self.cpp(self.d, [arg, "d.fa"]), run_py(self.d, [arg, "d.fa"]))


@unittest.skipUnless(cpp_binary(), "coal-cpp not built (COAL_CPP_BIN)")
class TestCodonDndsIsaFlags(unittest.TestCase):
    """`coal pipeline`'s native step takes the testing flags of `dnds` / `kaks` (AUDIT-85 L2)."""

    def setUp(self):
        self.bin = str(cpp_binary())
        self.d = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def run_cpp(self, *args, stdin: bytes | None = None):
        return subprocess.run([self.bin, *map(str, args)], input=stdin, capture_output=True, timeout=600)

    def test_isa_is_printed_as_dnds_prints_it(self):
        for args in ((), ("--kaks-baseline",)):
            with self.subTest(args=args):
                a = self.run_cpp("codon-dnds", "--kaks-isa", *args)
                b = self.run_cpp("dnds", "--kaks-isa", *args)
                self.assertEqual((a.returncode, a.stdout, a.stderr), (0, b.stdout, b""))
                self.assertIn(a.stdout, (b"avx2\n", b"baseline\n"))
        self.assertEqual(self.run_cpp("codon-dnds", "--kaks-isa", "--kaks-baseline").stdout, b"baseline\n")

    def test_hidden_from_the_help(self):
        text = self.run_cpp("codon-dnds", "--help").stdout
        self.assertNotIn(b"--kaks-baseline", text)
        self.assertNotIn(b"--kaks-isa", text)

    def test_baseline_gives_the_same_bytes(self):
        from test_dnds_pipeline import PLAIN, Q, S, read_fasta, stream
        q, s = read_fasta(Q), read_fasta(S)
        data = stream([(f"{x}&{y}", "all", [("M", x, q[x]), ("M", y, s[y])]) for x, y in PLAIN])
        outs = []
        for extra in ((), ("--kaks-baseline",)):
            o = self.d / f"o{len(outs)}.tsv"
            r = self.run_cpp("codon-dnds", "-o", o, "-m", "YN", "-m", "MYN", "--digits", "17", "-j", "2", *extra,
                             stdin=data)
            self.assertEqual(r.returncode, 0, r.stderr)
            outs.append(o.read_bytes())
        self.assertEqual(outs[0], outs[1])
        self.assertGreater(outs[0].count(b"\tYN\t"), 5)


@unittest.skipUnless(cpp_binary(), "coal-cpp not built (COAL_CPP_BIN)")
class TestMynInfiniteKappaRateRatio(unittest.TestCase):
    """KaKs_Calculator prints `inf:0.994949:1:1:1:1` for this pair (kappa_TC infinite); COAL prints NA."""
    AXT = ROOT / "tests" / "kaks" / "myn_kappa_tc_inf.axt"   # (the property fuzz's pair, seed 8585)

    def test_rate_ratio_is_na(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "o.kaks"
            r = subprocess.run([str(cpp_binary()), "kaks", "-i", str(self.AXT), "-m", "MYN", "-o", str(out)],
                               capture_output=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr)
            head, row = out.read_text().splitlines()
            cols = dict(zip(head.split("\t"), row.split("\t")))
            self.assertEqual(cols["Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)"], "NA")
            self.assertNotIn("inf", row.lower())
            self.assertNotIn("nan", row.lower())
            # the rest of the row as before the fix
            self.assertEqual(row, "c1230\tMYN\tNA\tNA\tNA\tNA\t147\tNA\tNA\tNA\t97\tNA\tNA\tNA\tNA\tNA\tNA\t"
                                  "0.377551(0.561224:0.438776:0.132653)\tNA\tNA\tNA\tNA")
            r = subprocess.run([str(cpp_binary()), "dnds", str(self.AXT), "-m", "MYN", "-q"], capture_output=True,
                               timeout=120)
            self.assertEqual(r.stdout.decode().splitlines()[1].split("\t")[13:],
                             ["saturated", "MYN: kappa_TC estimated as infinite (TN93 transition correction undefined "
                                           "in a class of sites): no site counts, Ka and Ks NA"])



@unittest.skipUnless(cpp_binary(), "coal-cpp not built (COAL_CPP_BIN)")
class TestMlwlKappaBoundary(unittest.TestCase):
    """tests/kaks/mlwl_kappa_boundary.axt: pairs whose MLWL kappa = 2 ts / tv is exactly 0.5 or 2. The doubles ts
    and tv are summed in column order and came out as 0.50000000000000011 or 2.0000000000000004 for some orders of
    the same codon pairs, which chose another site formula (N and S a whole site or more apart). c1359_code13: the
    property fuzz's pair (seed 4242; ambiguous base and gap codon, code 13), kappa 0.5 (N / S were 13.67 / 4.33 or
    15 / 3 by the order); mlwl_order_code1: a random pair, code 1, kappa 2 (N / S were 79 / 26 or 83.5 / 21.5 by the
    order; KaKs_Calculator prints 83.5 / 21.5 for the order of the file); mlwl_kappa2_code1: code 1, kappa 2, the
    same in every order before too (13.67 / 4.33, as KaKs_Calculator prints it), now 13.29 / 4.71, the formula of
    kappa = 2 exactly, as in every other case."""
    AXT = ROOT / "tests" / "kaks" / "mlwl_kappa_boundary.axt"
    EXPECT = {"c1359_code13": (13, ("13.6666666667", "4.33333333333")),
              "mlwl_order_code1": (1, ("79", "26")),
              "mlwl_kappa2_code1": (1, ("13.2857142857", "4.71428571429"))}

    def test_every_order_of_the_codon_columns_gives_the_same_sites(self):
        import random
        recs = {}
        for block in self.AXT.read_text(encoding="utf-8").split("\n\n"):
            lines = block.split()
            if lines:
                recs[lines[0]] = lines[1], lines[2]
        self.assertEqual(sorted(recs), sorted(self.EXPECT))
        rng = random.Random(4242)
        with tempfile.TemporaryDirectory() as tmp:
            for name, (code, sites) in self.EXPECT.items():
                with self.subTest(pair=name):
                    a, b = recs[name]
                    ca, cb = [a[i:i + 3] for i in range(0, len(a), 3)], [b[i:i + 3] for i in range(0, len(b), 3)]
                    text = []
                    for k in range(40):
                        order = list(range(len(ca)))
                        if k:
                            rng.shuffle(order)
                        text.append(f"o{k}\n{''.join(ca[i] for i in order)}\n{''.join(cb[i] for i in order)}\n")
                    path = Path(tmp) / f"{name}.axt"
                    path.write_text("\n".join(text), encoding="utf-8")
                    r = subprocess.run([str(cpp_binary()), "dnds", str(path), "-c", str(code), "-m", "MLWL", "-q",
                                        "--digits", "12"], capture_output=True, timeout=120)
                    self.assertEqual(r.returncode, 0, r.stderr)
                    rows = [ln.split("\t") for ln in r.stdout.decode().splitlines()[1:]]
                    self.assertEqual(len(rows), 40)
                    self.assertEqual({(x[8], x[9]) for x in rows}, {sites})
                    self.assertEqual(len({(x[3], x[4]) for x in rows}), 1)   # dN, dS (12 digits)



@unittest.skipUnless(cpp_binary(), "coal-cpp not built (COAL_CPP_BIN)")
class TestMlwlLargeDsIsGenuine(unittest.TestCase):
    """FIXCI-B3 (seed 4243, c1370): 513 Gln codon pairs CAG / CAA (256 of them differ, all transitions at the third
    position, 2-fold sites only) and a final TAG. 1 - 2P = 1 - 512/513 = 1/513 exactly: no rounding residue. LWL:
    S = L2 / 3 = 171, dS = 3 A2 = 1.5 ln 513; MLWL (kappa = 1: no transversion): S = L2 / 5 = 102.6, dS = 5 A2 =
    2.5 ln 513 = 15.6007, as KaKs_Calculator 2.0 prints. The same with code 1 and code 14."""

    def test_values(self):
        import math
        a = ["CAG"] * 215 + ["CAG"] * 210 + ["CAA"] * 47 + ["CAA"] * 41 + ["TAG"]
        b = ["CAA"] * 215 + ["CAG"] * 210 + ["CAA"] * 47 + ["CAG"] * 41 + ["TAG"]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "c1370.axt"
            path.write_text(f"c1370\n{''.join(a)}\n{''.join(b)}\n", encoding="utf-8")
            for code in (1, 14):
                with self.subTest(code=code):
                    r = subprocess.run([str(cpp_binary()), "dnds", str(path), "-c", str(code), "-m", "LWL", "-m",
                                        "MLWL", "--digits", "17", "-q"], capture_output=True, timeout=120)
                    self.assertEqual(r.returncode, 0, r.stderr)
                    rows = {x.split("\t")[2]: x.split("\t") for x in r.stdout.decode().splitlines()[1:]}
                    for m, k, s in (("LWL", 1.5, "171"), ("MLWL", 2.5, "102.59999999999999")):
                        row = rows[m]
                        self.assertEqual((row[3], row[9], row[13]), ("0", s, "ok"), m)
                        self.assertAlmostEqual(float(row[4]) / (k * math.log(513)), 1, delta=1e-13)

    def test_the_property_fuzz_excludes_it_and_nothing_larger(self):
        import test_kaks_properties as tkp
        case = tkp.Case("c1370", "CAG" * 514, "CAA" * 514, 14)
        f = tkp.Finding("P10", case, "MLWL", "P10: dS = 15.600689612926853 with status ok: a rounding residue")
        self.assertEqual(tkp.known(f), "M01")
        for other in (tkp.Finding("P10", case, "MLWL", "P10: dS = 60.1 with status ok: a rounding residue"),
                      tkp.Finding("P10", case, "LWL", "P10: dS = 15.6 with status ok: a rounding residue"),
                      tkp.Finding("P10", case, "MLWL", "P10: dN = 15.6 with status ok: a rounding residue")):
            self.assertEqual(tkp.known(other), "", other.detail)


if __name__ == "__main__":
    unittest.main()
