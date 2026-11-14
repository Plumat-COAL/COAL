"""MYN when kappa_AG (or kappa_TC) is estimated as infinite: status `saturated`, not `no_data`.

MYN estimates kappa_TC and kappa_AG from TN93 at two classes of sites. When the A-G term of a class is exactly at
its limit (1 - R*P1/(2*ag) - Q/(2*R) == 0: P1 = 3/25, Q = 12/25 with equal A and G frequencies in the first pair
below), its logarithm is -inf and kappa_AG is +inf. Until 2026-11-14 the site counts then came out as NaN and the
row was `no_data` with the note "no synonymous sites; no nonsynonymous sites", although the pair has 100 codons of
sites, and `coal kaks` on such pairs alone exited 2 suggesting that the input was not a codon alignment. The
multiple-hit correction of that class is undefined, which the status policy calls `saturated` (docs/kaks.md,
"Policy" and "Differences", 15). KaKs_Calculator 2.0 prints `nan` for every number of such a row and the same rate
ratio (`1.23581:inf:1:1:1:1` for the first pair); COAL prints NA. The simulation study of 2026-11-14 had 79 such
MYN rows among 240,000 pairs (dS 0.5 to 3); in real genes 3 of about 16,000 MYN rows (one Salmonella / Shigella
pair, two S. cerevisiae / C. glabrata pairs).

The pairs are tests/kaks/myn_kappa_inf.axt. tests/cpp/kaks_myn_kappa_inf.cpp checks the same through the C++ API
(compiled here against the build's library when a C++ compiler and libcoal_core are found). Needs build/coal-cpp
(COAL_CPP_BIN).
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests"))
from test_dnds_pipeline import CPP, needs_cpp, run_cpp  # noqa: E402

AXT = ROOT / "tests" / "kaks" / "myn_kappa_inf.axt"
NAMES = ["w0.2_dS1.5_r19", "w0.2_dS3_r1", "w1_dS1.5_r18"]
# What KaKs_Calculator 2.0 (bioconda kakscalculator2) prints for the three pairs with -m MYN: every number nan,
# the sites -nan, and this rate ratio; Substitutions is its count of the iteration that lost the codons with
# several differences (docs/kaks.md, "Differences", 11), so it is not compared.
KC_RATES = {"w0.2_dS1.5_r19": "1.23581:inf:1:1:1:1", "w0.2_dS3_r1": "2.15044:inf:1:1:1:1",
            "w1_dS1.5_r18": "1.99961:inf:1:1:1:1"}
SUBSTITUTIONS = {"w0.2_dS1.5_r19": "121", "w0.2_dS3_r1": "155", "w1_dS1.5_r18": "185"}
NOTE = ("MYN: kappa_AG estimated as infinite (TN93 transition correction undefined in a class of sites): no site "
        "counts, Ka and Ks NA")


def tsv(text: str) -> list[dict]:
    lines = text.splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]


def one_pair(name: str, d: Path) -> Path:
    recs = AXT.read_text(encoding="utf-8").split("\n\n")
    for r in recs:
        lines = [x for x in r.splitlines() if x and not x.startswith("#")]
        if lines and lines[0].split("&")[0] == name:
            p = d / f"{name}.axt"
            p.write_text("\n".join(lines) + "\n\n", encoding="utf-8")
            return p
    raise KeyError(name)


@needs_cpp
class TestMynKappaInfinite(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="myn_kinf_"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_dnds_rows_are_saturated_with_a_note_naming_kappa(self):
        r = run_cpp("dnds", AXT, "-m", "MYN", "-m", "YN", "--digits", "8")
        self.assertEqual(r.returncode, 0, r.stderr)
        rows = {(x["seq1"], x["method"]): x for x in tsv(r.stdout.decode())}
        for n in NAMES:
            m = rows[(n, "MYN")]
            self.assertEqual(m["status"], "saturated", n)
            self.assertEqual(m["note"], NOTE, n)
            self.assertNotIn("no synonymous sites", m["note"])
            for c in ("dN", "dS", "dNdS", "SE_dN", "SE_dS", "N_sites", "S_sites", "Nd", "Sd"):
                self.assertEqual(m[c], "NA", f"{n} {c}")
            self.assertEqual(m["codons_compared"], "100")
            # YN on the same pairs is not affected
            y = rows[(n, "YN")]
            self.assertEqual(y["status"], "ok", n)
            self.assertNotEqual(y["dN"], "NA")
            self.assertNotEqual(y["dS"], "NA")
        self.assertIn(b"3 saturated, 0 no_data", r.stderr)

    def test_kaks_format_exit_status_and_kaks_calculator_columns(self):
        out = self.tmp / "o.kaks"
        r = run_cpp("kaks", "-i", AXT, "-m", "MYN", "-o", out)
        self.assertEqual(r.returncode, 0, r.stderr)        # an estimate-less but explained result: not exit 2
        self.assertNotIn(b"codon alignment", r.stderr)
        self.assertIn(b"0 ok, 0 warning, 0 not_converged, 3 saturated, 0 no_data, 0 error", r.stderr)
        rows = {x["Sequence"].split("&")[0]: x for x in tsv(out.read_text(encoding="utf-8"))}
        for n in NAMES:
            x = rows[n]
            self.assertEqual(x["Method"], "MYN")
            for c in ("Ka", "Ks", "Ka/Ks", "P-Value(Fisher)", "S-Sites", "N-Sites", "S-Substitutions",
                      "N-Substitutions", "Divergence-Time"):
                self.assertEqual(x[c], "NA", f"{n} {c}")
            self.assertEqual(x["Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)"], KC_RATES[n], n)
            self.assertEqual(x["Substitutions"], SUBSTITUTIONS[n], n)
            self.assertEqual(x["Length"], "300")

    def test_a_single_such_pair_exits_0_without_the_codon_alignment_hint(self):
        for n in NAMES:
            p = one_pair(n, self.tmp)
            for sub in (["kaks", "-i", p, "-m", "MYN", "-o", self.tmp / "o.kaks"],
                        ["dnds", p, "-m", "MYN", "-o", self.tmp / "o.tsv"]):
                r = run_cpp(*sub)
                self.assertEqual(r.returncode, 0, f"{n} {sub[0]}: {r.stderr!r}")
                self.assertNotIn(b"none of the", r.stderr)
                self.assertNotIn(b"codon alignment", r.stderr)
                r = run_cpp(*sub, "-q")
                self.assertEqual((r.returncode, r.stderr), (0, b""), f"{n} {sub[0]} -q")

    def test_python_command_reports_it_the_same_way(self):
        env = dict(os.environ, COAL_CPP=str(CPP), PYTHONPATH=str(ROOT / "src"))
        a, b = self.tmp / "py.kaks", self.tmp / "cpp.kaks"
        r = subprocess.run([sys.executable, "-m", "coal", "kaks", "-i", str(AXT), "-m", "MYN", "-o", str(a)],
                           capture_output=True, env=env)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn(b"3 saturated, 0 no_data", r.stderr)
        self.assertNotIn(b"codon alignment", r.stderr)
        self.assertEqual(run_cpp("kaks", "-i", AXT, "-m", "MYN", "-o", b).returncode, 0)
        self.assertEqual(a.read_bytes(), b.read_bytes())

    def test_swapped_rows_and_threads_give_the_same_rows(self):
        recs = [r for r in AXT.read_text(encoding="utf-8").split("\n\n") if r.strip()]
        sw = []
        for r in recs:
            lines = [x for x in r.splitlines() if x and not x.startswith("#")]
            sw.append(f"{lines[0]}\n{lines[2]}\n{lines[1]}\n")
        p = self.tmp / "swapped.axt"
        p.write_text("\n".join(sw) + "\n", encoding="utf-8")
        base = run_cpp("dnds", AXT, "-m", "MYN", "-q", "-j", "1").stdout
        self.assertEqual(run_cpp("dnds", p, "-m", "MYN", "-q", "-j", "1").stdout, base)
        self.assertEqual(run_cpp("dnds", AXT, "-m", "MYN", "-q", "-j", "3").stdout, base)


def _api_build():
    """(compiler, include dir, library) for tests/cpp/kaks_myn_kappa_inf.cpp, or None."""
    if CPP is None:
        return None
    build = Path(CPP).resolve().parent
    inc = build / "coal_include"
    lib = next((p for p in (build / "libcoal_core.a",) if p.exists()), None)
    cxx = os.environ.get("CXX") or shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")
    if not (inc / "coal" / "kaks.hpp").exists() or lib is None or cxx is None:
        return None
    return cxx, inc, lib


@unittest.skipUnless(_api_build(), "needs a C++ compiler and the build's coal_include/ and libcoal_core.a")
class TestMynKappaInfiniteCppApi(unittest.TestCase):
    def test_cpp_api(self):
        cxx, inc, lib = _api_build()
        with tempfile.TemporaryDirectory() as td:
            exe = Path(td) / "kaks_myn_kappa_inf"
            src = ROOT / "tests" / "cpp" / "kaks_myn_kappa_inf.cpp"
            r = subprocess.run([cxx, "-std=c++17", "-O1", "-Wall", "-Wextra", "-I", str(inc), str(src), str(lib),
                                "-o", str(exe), "-pthread"], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            for r_ in [r for r in AXT.read_text(encoding="utf-8").split("\n\n") if r.strip()]:
                lines = [x for x in r_.splitlines() if x and not x.startswith("#")]
                p = subprocess.run([str(exe), lines[1], lines[2]], capture_output=True, text=True)
                self.assertEqual((p.returncode, p.stdout.strip()), (0, "ok"), f"{lines[0]}: {p.stderr}")


if __name__ == "__main__":
    unittest.main()
