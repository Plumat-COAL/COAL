"""pwalign's copies of the substitution matrices (BLOSUM62-PWALIGN, BLOSUM45-PWALIGN, BLOSUM80-PWALIGN) and
--preset pwalign with any matrix: the protein alignments of pwalign::pairwiseAlignment, alignment and score.

tests/fixtures/pwalign_matrices: 104 real protein pairs (make_pairs.py) with B, Z, X, J, U, O and '*' in some of them,
and pwalign 1.6.0's alignments of them with BLOSUM62, BLOSUM45, BLOSUM80 and PAM250 (make_expected.R). pwalign refuses
U and O (and J with PAM250, which has no J row): those pairs' expected alignments are pwalign's with the letters
read as X, which is what COAL does ("as_x=" in the name line).
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import align as A  # noqa: E402
from coal import dnds_pipeline as P  # noqa: E402
from test_align import align_binary, no_binary  # noqa: E402
from test_dnds_pipeline import Q, S, PipelineCase, needs_cpp  # noqa: E402
from test_integration1_findings import _fixture_pair_ids  # noqa: E402

FX = ROOT / "tests" / "fixtures" / "pwalign_matrices"
MATRICES = ("BLOSUM62", "BLOSUM45", "BLOSUM80", "PAM250")
TIMEOUT = 300


def fasta(path: Path) -> list[tuple[str, str]]:
    out: list[list[str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            out.append([line[1:], ""])
        elif out:
            out[-1][1] += line.strip()
    return [(n, s) for n, s in out]


def fixture(matrix: str):
    """[(a, b, expected_a, expected_b, score, letters pwalign read as X)] for one matrix."""
    pairs, exp = fasta(FX / "pairs.fasta"), fasta(FX / f"expected_{matrix}.fasta")
    out = []
    for k in range(0, len(pairs), 2):
        name = exp[k][0]
        score = float(re.search(r" score=(\S+)", name).group(1))
        as_x = re.search(r" as_x=(\S+)", name)
        out.append((pairs[k][1], pairs[k + 1][1], exp[k][1], exp[k + 1][1], score, as_x.group(1) if as_x else ""))
    return out


def as_pwalign(row: str, as_x: str) -> str:
    return row.translate(str.maketrans({c: "X" for c in as_x})) if as_x else row


def finish(*argv: str):
    a = P.build_parser().parse_args(["pairs", "--query", "q", "--subject", "s", "-o", "o", *argv])
    P._finish_args(a)
    return a


class Matrices(unittest.TestCase):
    def test_alphabet_and_shape(self):
        self.assertEqual(A.ALPHABET, "ARNDCQEGHILKMFPSTWYVBZX*J")
        for name, m in A.MATRICES.items():
            self.assertEqual(len(m), 25, name)
            self.assertTrue(all(len(row) == 25 for row in m), name)
            self.assertTrue(all(m[i][j] == m[j][i] for i in range(25) for j in range(25)), f"{name} not symmetric")
        for name in ("BLOSUM62-PWALIGN", "BLOSUM45-PWALIGN", "BLOSUM80-PWALIGN"):
            self.assertIn(name, A.MATRICES)
        self.assertNotIn("PAM250-PWALIGN", A.MATRICES)   # pwalign's PAM250 is NCBI's

    def test_ncbi_matrices_read_j_o_u_as_x(self):
        x, j = A.ALPHABET.index("X"), A.ALPHABET.index("J")
        for name in ("BLOSUM62", "BLOSUM45", "BLOSUM80", "PAM250"):
            m = A.MATRICES[name]
            self.assertEqual(m[j], [m[x][k] if k != j else m[x][x] for k in range(25)], name)
        for a in "OU":
            self.assertEqual(A._CODE[a], x)
            self.assertEqual(A._CODE[a.lower()], x)
        self.assertEqual(A._CODE["J"], j)

    def test_pwalign_copies_differ_where_pwalign_differs(self):
        al = A.ALPHABET
        for base, letters in (("BLOSUM62", "BZXJ"), ("BLOSUM45", "BZXJ")):
            pw, ncbi = A.MATRICES[base + "-PWALIGN"], A.MATRICES[base]
            diff = {(al[i], al[j]) for i in range(25) for j in range(25) if pw[i][j] != ncbi[i][j]}
            self.assertTrue(diff, base)
            self.assertTrue(all(set(d) & set(letters) for d in diff), (base, diff))
        for name in ("BLOSUM62-PWALIGN", "BLOSUM45-PWALIGN", "BLOSUM80-PWALIGN"):
            m, x = A.MATRICES[name], al.index("X")
            self.assertTrue(all(m[x][k] == -1 for k in range(25) if al[k] != "*"), name)
            self.assertGreater(m[al.index("J")][al.index("I")], 0, name)   # J is I or L, not X
        # pwalign's BLOSUM80 is the matrix on another scale
        pw, ncbi = A.MATRICES["BLOSUM80-PWALIGN"], A.MATRICES["BLOSUM80"]
        self.assertEqual((pw[0][0], pw[4][4], pw[0][23]), (5, 9, -6))
        self.assertEqual((ncbi[0][0], ncbi[4][4], ncbi[0][23]), (7, 13, -8))

    def test_pwalign_matrix_names(self):
        self.assertEqual([A.pwalign_matrix(n) for n in ("BLOSUM62", "blosum45", "Blosum80", "pam250")],
                         ["BLOSUM62-PWALIGN", "BLOSUM45-PWALIGN", "BLOSUM80-PWALIGN", "PAM250"])
        self.assertEqual(A.pwalign_matrix("blosum45-pwalign"), "BLOSUM45-PWALIGN")
        self.assertEqual(A.pwalign_matrix("BLOSUM99"), "BLOSUM99")
        for name in A.MATRICES:
            self.assertIn(A.pwalign_matrix(name), A.MATRICES)
        self.assertEqual(A.AlignParams.pwalign(), A.AlignParams("BLOSUM62-PWALIGN", 14, 4, False, "pwalign"))
        self.assertEqual(A.AlignParams.pwalign("BLOSUM80"), A.AlignParams("BLOSUM80-PWALIGN", 14, 4, False, "pwalign"))
        self.assertEqual(A.AlignParams.pwalign("PAM250").matrix, "PAM250")

    def test_matrices_match_cpp_and_vendored_sources(self):
        for src in (ROOT / "cpp" / "src" / "align.cpp", ROOT / "packaging" / "r" / "coalr" / "src" / "coal" / "align.cpp"):
            text = src.read_text(encoding="utf-8")
            names = re.findall(r'^    \{"([A-Z0-9-]+)", \{$', text, re.M)
            self.assertEqual(names, list(A.MATRICES), src)
            for name, mat in A.MATRICES.items():
                block = text[text.index(f'{{"{name}", {{'):]
                nums = [int(x) for x in block[block.index("{", 3):block.index("}},")].replace("{", " ").replace(
                    "}", " ").replace(",", " ").split()]
                self.assertEqual(nums, [v for row in mat for v in row], (src, name))


class Preset(unittest.TestCase):
    def test_pipeline_preset_takes_pwalign_copy_of_the_matrix(self):
        self.assertEqual(finish("--preset", "pwalign").matrix, "BLOSUM62-PWALIGN")
        for given, want in (("BLOSUM62", "BLOSUM62-PWALIGN"), ("blosum45", "BLOSUM45-PWALIGN"),
                            ("BLOSUM80", "BLOSUM80-PWALIGN"), ("PAM250", "PAM250"),
                            ("BLOSUM45-PWALIGN", "BLOSUM45-PWALIGN")):
            self.assertEqual(finish("--preset", "pwalign", "--matrix", given).matrix, want, given)
        # without the preset (or with needle) a name is the matrix of that name
        for given in ("BLOSUM62", "BLOSUM45", "BLOSUM80", "PAM250", "BLOSUM80-PWALIGN"):
            self.assertEqual(finish("--matrix", given).matrix, given)
            self.assertEqual(finish("--preset", "needle", "--matrix", given).matrix, given)
        self.assertEqual(finish().matrix, "BLOSUM62")


class Identity(unittest.TestCase):
    """COAL's --preset pwalign against pwalign 1.6.0, every matrix: the same rows and the same score."""

    def test_fixture(self):
        pairs = fasta(FX / "pairs.fasta")
        self.assertGreaterEqual(len(pairs) // 2, 100)
        text = "".join(s for _, s in pairs)
        for c in "BZXJUO*":
            self.assertIn(c, text)
        for m in MATRICES:
            rows = fixture(m)
            self.assertEqual(len(rows), len(pairs) // 2)
            # the pairs pwalign could align as given: every pair without U and O (and J for PAM250)
            native = [r for r in rows if not r[5]]
            self.assertGreaterEqual(len(native), 70, m)
            for a, b, *_rest, as_x in rows:
                foreign = "UOJ" if m == "PAM250" else "UO"
                self.assertEqual(bool(as_x), any(c in a + b for c in foreign), (m, a))

    def test_matrices_change_the_alignments(self):
        # the fixture only tests the copies if NCBI's matrices give other alignments or scores for some pairs
        for m in ("BLOSUM62", "BLOSUM45", "BLOSUM80"):
            ncbi = A.AlignParams(m, 14, 4, False, "pwalign")
            differ = 0
            for a, b, ea, eb, score, as_x in fixture(m)[::4]:
                r = A.align_pair(a, b, ncbi)
                differ += (as_pwalign(r.a, as_x), as_pwalign(r.b, as_x), r.score) != (ea, eb, score)
            self.assertGreater(differ, 0, m)

    def test_python_kernel(self):
        for m in MATRICES:
            p = A.AlignParams.pwalign(m)
            for a, b, ea, eb, score, as_x in fixture(m)[::4]:   # pure Python is slow: every fourth pair
                r = A.align_pair(a, b, p)
                self.assertEqual((as_pwalign(r.a, as_x), as_pwalign(r.b, as_x)), (ea, eb), (m, a))
                self.assertEqual(r.score, score, (m, a))
                self.assertEqual(A.score_alignment(r.a, r.b, p), score)

    def test_cpp_align_binary(self):
        exe = align_binary()
        if exe is None:
            no_binary(self)
        for m in MATRICES:
            rows = fixture(m)
            scores = Path(self.tmp_scores())
            out = subprocess.run([str(exe), "--preset", "pwalign", "-m", m, "-j", "2", "--scores", str(scores),
                                  str(FX / "pairs.fasta")], capture_output=True, text=True, timeout=TIMEOUT)
            self.assertEqual(out.returncode, 0, out.stderr)
            got = [s for _, s in fasta_text(out.stdout)]
            got_scores = [float(line.split("\t")[6]) for line in scores.read_text(encoding="utf-8").splitlines()[1:]]
            self.assertEqual(len(got), 2 * len(rows))
            for k, (a, b, ea, eb, score, as_x) in enumerate(rows):
                self.assertEqual((as_pwalign(got[2 * k], as_x), as_pwalign(got[2 * k + 1], as_x)), (ea, eb), (m, a))
                self.assertEqual(got_scores[k], score, (m, a))
        # the matrix names as given: explicit *-PWALIGN names, and -m after --preset pwalign
        for args, want in ((["-m", "BLOSUM45-PWALIGN", "--preset", "pwalign"], "BLOSUM62"),
                           (["--preset", "pwalign", "-m", "BLOSUM45-PWALIGN"], "BLOSUM45"),
                           (["--preset", "pwalign", "-mblosum45"], "BLOSUM45")):
            out = subprocess.run([str(exe), *args, str(FX / "pairs.fasta")], capture_output=True, text=True,
                                 timeout=TIMEOUT)
            self.assertEqual(out.returncode, 0, out.stderr)
            got = [s for _, s in fasta_text(out.stdout)]
            for k, (a, b, ea, eb, _score, as_x) in enumerate(fixture(want)[:20]):
                self.assertEqual((as_pwalign(got[2 * k], as_x), as_pwalign(got[2 * k + 1], as_x)), (ea, eb), args)
        bad = subprocess.run([str(exe), "-m", "PAM250-PWALIGN", str(FX / "pairs.fasta")], capture_output=True,
                             text=True, timeout=TIMEOUT)
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn("BLOSUM80-PWALIGN", bad.stderr)   # the known names

    def tmp_scores(self) -> str:
        import tempfile
        fd, path = tempfile.mkstemp(suffix=".tsv")
        os.close(fd)
        self.addCleanup(os.remove, path)
        return path


def fasta_text(text: str) -> list[tuple[str, str]]:
    out: list[list[str]] = []
    for line in text.splitlines():
        if line.startswith(">"):
            out.append([line[1:], ""])
        elif out:
            out[-1][1] += line.strip()
    return [(n, s) for n, s in out]


@needs_cpp
class Pipeline(PipelineCase):
    def test_preset_pwalign_with_blosum45_is_recorded(self):
        pairs = str(self.pairs_file(_fixture_pair_ids()))
        res = P.run_pairs(Q, S, self.tmp / "pw45", pairs=pairs, preset="pwalign", matrix="BLOSUM45", jobs=1,
                          quiet=True)
        self.assertIn(res.exit_code, (0, 2))
        self.check_report(res)
        self.assertEqual(res.report["parameters"]["matrix"], "BLOSUM45-PWALIGN")
        cmd = [c for c in res.report["tools"]["coal-cpp"]["commands"] if "codon-dnds" in c][0]
        self.assertEqual(cmd[cmd.index("--matrix") + 1], "BLOSUM45-PWALIGN")
        plain = P.run_pairs(Q, S, self.tmp / "b45", pairs=pairs, matrix="BLOSUM45", jobs=1, quiet=True)
        self.assertEqual(plain.report["parameters"]["matrix"], "BLOSUM45")


if __name__ == "__main__":
    unittest.main()
