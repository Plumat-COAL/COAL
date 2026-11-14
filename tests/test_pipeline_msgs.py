"""Pipeline messages and one ordering rule, found by the tutorial writers and testers (2026-11-14):

  * CLUSTER-ORDER: the pairs of a cluster aligned by Clustal Omega / MAFFT are computed in input order (the gene of
    the earlier --cds file, then the earlier record, first), as the pairwise-aligned pairs are; the rows still list
    the clusters file's order. Engine side: the `O` line of `coal-cpp codon-dnds` (dnds.hpp Group::rank).
  * --resume does not depend on where the output folder is: a renamed or moved folder reuses its stages.
  * The warning for several CDS records per gene names alternative start sites (bacterial annotations) as well as
    isoforms.

Synthetic inputs only; coal-cpp (COAL_CPP_BIN). The MSA tool is a stand-in script (POSIX) that returns its input,
so no Clustal Omega is needed.
"""
from __future__ import annotations

import csv
import io
import os
import random
import shutil
import stat
import subprocess
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from coal import dnds_pipeline as P  # noqa: E402
from test_dnds_pipeline import CPP, PipelineCase, needs_cpp  # noqa: E402

TIMEOUT = 300
SENSE = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG" if a + b + c not in ("TAA", "TAG", "TGA")]
# libsequence's two-position branches: with them Comeron depends on which gene is first, so the order is visible
LIB = {"comeron_twosubs": "libsequence"}
METHODS = ["YN", "LWL", "Comeron"]


def family(rng: random.Random, n: int, length: int = 150) -> list[str]:
    """n related CDS of the same length (substitutions only): every alignment of them is gapless, so the MSA a
    stand-in tool returns (its input) is the alignment any aligner would make."""
    base = ["ATG"] + [rng.choice(SENSE) for _ in range(length)]
    out = []
    for _ in range(n):
        s = list(base)
        for i in range(1, len(s)):
            if rng.random() < 0.12:
                s[i] = rng.choice(SENSE)
        out.append("".join(s) + "TAA")
    return out


CODE = dict(zip([a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"],
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def translate(cds: str) -> str:
    return "".join(CODE[cds[i:i + 3]] for i in range(0, len(cds), 3))


def rows(path) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def engine(text: str, tmp: Path, *extra: str):
    out = tmp / "engine.tsv"
    r = subprocess.run([str(CPP), "codon-dnds", "-o", str(out), "-m", "YN", "-m", "Comeron", "-q",
                        "--comeron-twosubs", "libsequence", *extra],
                       input=text, capture_output=True, text=True, timeout=TIMEOUT)
    return r, (out.read_text(encoding="utf-8").splitlines()[1:] if r.returncode in (0, 2) else [])


@needs_cpp
class EngineRank(PipelineCase):
    """`O <rank> ...`: every pair of a group is computed with the member of the lower rank first, listed as given."""

    def test_rank_line(self):
        a, b, c = family(random.Random(3), 3)
        # rows of a gapless MSA: the translations
        row = {k: translate(s[:-3]) for k, s in zip("abc", (a, b, c))}

        def group(order, rank=None):
            g = "G\tg\tall\n" + "".join(f"M\t{k}\t{dict(a=a, b=b, c=c)[k]}\t{row[k]}\n" for k in order)
            return g + (("O\t" + "\t".join(map(str, rank)) + "\n") if rank is not None else "") + "E\n"
        # members listed c, a, b; input order a < b < c
        _, plain = engine(group("cab"), self.tmp)
        _, ranked = engine(group("cab", [2, 0, 1]), self.tmp)
        _, ordered = engine(group("abc"), self.tmp)
        self.assertEqual(len(ranked), 6)
        # listed as given: c-a, c-b, a-b
        self.assertEqual([r.split("\t")[:3] for r in ranked], [r.split("\t")[:3] for r in plain])
        self.assertEqual([tuple(r.split("\t")[:2]) for r in ranked][::2], [("c", "a"), ("c", "b"), ("a", "b")])
        # computed in input order: the values of the group given in input order (pairs a-c, b-c, a-b)
        by = {(r.split("\t")[0], r.split("\t")[1], r.split("\t")[2]): r.split("\t")[3:] for r in ordered}
        for r in ranked:
            f = r.split("\t")
            x, y = sorted(f[:2])
            self.assertEqual(f[3:], by[(x, y, f[2])], f[:3])
        # and the test can tell: Comeron differs without the ranks for a swapped pair
        self.assertNotEqual(ranked, plain)
        # equal to member order when the ranks are in member order
        _, same = engine(group("abc", [0, 1, 2]), self.tmp)
        self.assertEqual(same, ordered)

    def test_bad_rank_lines(self):
        a, b, c = family(random.Random(4), 3)
        g = f"G\tg\tall\nM\ta\t{a}\nM\tb\t{b}\nM\tc\t{c}\n"
        bad, _ = engine(g + "O\t1\t0\nE\n", self.tmp)
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn("O gives 2 ranks for 3 members", bad.stderr)
        bad, _ = engine(g + "O\t1\tx\t0\nE\n", self.tmp)
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn("O needs one input rank", bad.stderr)


@needs_cpp
@unittest.skipIf(os.name == "nt", "the stand-in MSA tool is a POSIX shell script")
class MsaClusterOrder(PipelineCase):
    """A cluster aligned by the MSA tool gives the pairs route's values for each pair (input order), and lists its
    rows in the clusters file's order (the representative first)."""

    def setUp(self):
        super().setUp()
        tool = self.tmp / "clustalo"
        # returns its input (-i IN -o OUT ...): the sequences have one length, so that is their gapless MSA
        tool.write_text('#!/bin/sh\nif [ "$1" = "--version" ]; then echo 1.2.4; exit 0; fi\n'
                        'while [ $# -gt 0 ]; do case "$1" in -i) i="$2"; shift;; -o) o="$2"; shift;; esac; shift; '
                        'done\ncp "$i" "$o"\n', encoding="utf-8")
        tool.chmod(tool.stat().st_mode | stat.S_IXUSR)
        p = mock.patch.dict(os.environ, {"COAL_CLUSTALO": str(tool)})
        p.start()
        self.addCleanup(p.stop)

    def files(self):
        rng = random.Random(11)
        fams = [family(rng, 4) for _ in range(4)]
        f1, f2 = self.tmp / "sp1.fa", self.tmp / "sp2.fa"
        # sp1 (first --cds file): y, z, w of each family; sp2: x, the representative in the clusters file
        f1.write_text("".join(f">{n}{k}\n{s}\n" for k, fam in enumerate(fams) for n, s in zip("yzw", fam[1:])),
                      encoding="utf-8")
        f2.write_text("".join(f">x{k}\n{fam[0]}\n" for k, fam in enumerate(fams)), encoding="utf-8")
        cl = self.tmp / "clusters.tsv"
        # member order: x (last in the input), w, y, z: most pairs are listed later gene first
        cl.write_text("".join(f"x{k}\t{m}{k}\n" for k in range(len(fams)) for m in "xwyz"), encoding="utf-8")
        return f1, f2, cl

    def test_msa_pairs_in_input_order(self):
        f1, f2, cl = self.files()
        msa = P.run_clusters([f1, f2], self.tmp / "msa", clusters=cl, method=METHODS, jobs=1, quiet=True, **LIB)
        pw = P.run_clusters([f1, f2], self.tmp / "pw", clusters=cl, method=METHODS, jobs=1, quiet=True,
                            max_cluster_size=2, large_clusters="sample", pair_sample=100, **LIB)
        self.assertEqual(msa.exit_code, 0)
        self.assertEqual(pw.exit_code, 0)
        self.assertEqual({r["aligner"] for r in rows(msa.outputs["cluster_summary"])}, {"clustalo"})
        self.assertEqual({r["aligner"] for r in rows(pw.outputs["cluster_summary"])}, {"nw_pairwise"})
        m, p = rows(msa.outputs["dnds"]), rows(pw.outputs["dnds"])
        # the clusters file's order, the representative first: x-w, x-y, x-z, w-y, w-z, y-z
        self.assertEqual([(r["seq1"][0], r["seq2"][0]) for r in m if r["method"] == "YN"][:6],
                         [("x", "w"), ("x", "y"), ("x", "z"), ("w", "y"), ("w", "z"), ("y", "z")])
        key = lambda r: (frozenset((r["seq1"], r["seq2"])), r["method"])  # noqa: E731
        cols = ("dN", "dS", "dNdS", "SE_dN", "SE_dS", "N_sites", "S_sites", "codons_compared", "status")
        pv = {key(r): tuple(r[c] for c in cols) for r in p}
        self.assertEqual(len(m), len(p))
        self.assertEqual({key(r): tuple(r[c] for c in cols) for r in m}, pv)
        # the other order of the files gives other Comeron values (so the test sees the order)
        rev = P.run_clusters([f2, f1], self.tmp / "rev", clusters=cl, method=METHODS, jobs=1, quiet=True, **LIB)
        rv = {key(r): tuple(r[c] for c in cols) for r in rows(rev.outputs["dnds"])}
        self.assertTrue(any(rv[k] != pv[k] for k in pv if k[1] == "Comeron"))
        self.assertEqual({k: v for k, v in rv.items() if k[1] == "YN"}, {k: v for k, v in pv.items() if k[1] == "YN"})


def run_quiet_err(fn, *args, **kw):
    err = io.StringIO()
    with redirect_stderr(err):
        res = fn(*args, **kw)
    return res, err.getvalue()


@needs_cpp
class AlternativeStarts(PipelineCase):
    """NCBI's bacterial CDS files have genes with two or three CDS records of one locus tag that differ in their start
    (E. coli mrcB, infB): the warning names them as well as isoforms."""

    def files(self):
        rng = random.Random(5)
        a, b = family(rng, 2)
        c, d = family(rng, 2)
        q = self.tmp / "ecoli.fna"
        # mrcB: the long form and a form that starts 20 codons later (same stop)
        q.write_text(f">lcl|NC_1_cds_NP_1.1_1 [gene=mrcB] [locus_tag=b0149]\n{a}\n"
                     f">lcl|NC_1_cds_YP_2.1_2 [gene=mrcB] [locus_tag=b0149]\nATG{a[63:]}\n"
                     f">lcl|NC_1_cds_NP_3.1_3 [gene=other] [locus_tag=b0150]\n{c}\n", encoding="utf-8")
        s = self.tmp / "salm.fna"
        s.write_text(f">lcl|NC_2_cds_NP_9.1_1 [gene=mrcB] [locus_tag=STM0211]\n{b}\n"
                     f">lcl|NC_2_cds_NP_8.1_2 [gene=other] [locus_tag=STM0212]\n{d}\n", encoding="utf-8")
        return q, s

    def test_pairs_warning(self):
        q, s = self.files()
        pf = self.pairs_file([("lcl|NC_1_cds_NP_1.1_1", "lcl|NC_2_cds_NP_9.1_1"),
                              ("lcl|NC_1_cds_NP_3.1_3", "lcl|NC_2_cds_NP_8.1_2")])
        res, err = run_quiet_err(P.run_pairs, q, s, self.tmp / "o", pairs=pf, jobs=1)
        self.assertEqual(res.exit_code, 0, err)
        w = [ln for ln in err.splitlines() if "warning:" in ln]
        self.assertEqual(len(w), 1, err)
        self.assertIn("the input has several CDS records per gene (isoforms or alternative start sites; query: 3 "
                      "records of 2 genes, 1 gene with several records (e.g. b0149: ", w[0])
        self.assertIn("such a gene is compared once per record, and its reciprocal best hits are split between its "
                      "records. Use --longest-isoform", w[0])
        self.assertNotIn("transcript", w[0])
        # --longest-isoform keeps the long form
        res, err = run_quiet_err(P.run_pairs, q, s, self.tmp / "o2", pairs=pf, jobs=1, longest_isoform=True)
        self.assertEqual(res.exit_code, 0, err)
        genes = {r["id"]: r["status"] for r in rows(res.outputs["genes"])}
        self.assertEqual(genes["lcl|NC_1_cds_YP_2.1_2"], "dropped")
        self.assertEqual(genes["lcl|NC_1_cds_NP_1.1_1"], "ok")

    def test_clusters_warning_has_no_best_hits(self):
        q, s = self.files()
        cl = self.tmp / "cl.tsv"
        cl.write_text("lcl|NC_1_cds_NP_1.1_1\tlcl|NC_1_cds_NP_1.1_1\nlcl|NC_1_cds_NP_1.1_1\tlcl|NC_2_cds_NP_9.1_1\n",
                      encoding="utf-8")
        res, err = run_quiet_err(P.run_clusters, [q, s], self.tmp / "c", clusters=cl, jobs=1)
        w = [ln for ln in err.splitlines() if "several CDS records" in ln]
        self.assertEqual(len(w), 1, err)
        self.assertIn("such a gene is compared once per record. Use --longest-isoform", w[0])


if __name__ == "__main__":
    unittest.main()
