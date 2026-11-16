"""`coal pipeline` never writes through a symbolic link in its output folder (AUDIT-99, L1).

A link put at any name the pipeline writes, in DIR or DIR/work (outputs.json.tmp, state.json.tmp, report.json.tmp,
the .tmp of every table, the per-cluster .fasta.tmp, the logs, the decompressed and --prefix-ids copies of the
inputs, the translated proteins, DIAMOND's hits and databases, the sorted proteins for DeepClust) pointed at a
file outside the folder: the next run (and the same run with --resume) wrote into that file, among them the
JSON record, the logs and the proteins. Each test runs once, replaces every file of the folder with a link to a
"victim" outside it and plants a link at every .tmp and log name, runs again and with --resume, and checks that
no victim changed and that the names COAL wrote are files again, not links. POSIX (symbolic links); needs
build/coal-cpp (COAL_CPP_BIN), and DIAMOND for the DIAMOND routes.
"""
from __future__ import annotations

import gzip
import io
import os
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
from test_dnds_pipeline import (CLUSTALO, FAM, PLAIN, PipelineCase, Q, S, needs_cpp, needs_diamond,  # noqa: E402
                                read_fasta)

from coal import dnds_pipeline as P  # noqa: E402

# names a run writes that the first run may have removed again or written only through a temporary name
EXTRA = ("work/outputs.json.tmp", "work/state.json.tmp", "work/inputs/sources.json.tmp", "work/logs/translate.log",
         "work/logs/diamond.log", "work/logs/codon-dnds.log", "work/subject.dmnd", "work/query.dmnd",
         "realign.tsv.body")
VICTIM = b"keep me\n"


def run(fn, *args, **kw):
    with redirect_stderr(io.StringIO()) as err:
        res = fn(*args, **kw)
    return res, err.getvalue()


class _Planting(PipelineCase):
    def plant(self, out: Path) -> dict[str, Path]:
        """Every file of `out` (but the record itself, which is read: a link there makes the folder another's) and
        every .tmp and log name becomes a link to a file of its own outside `out`."""
        names = {str(p.relative_to(out)) for p in out.rglob("*") if p.is_file() and not p.is_symlink()}
        names |= {n + ".tmp" for n in names if not n.endswith(".tmp")} | set(EXTRA)
        names.discard("work/outputs.json")
        victims = self.tmp / f"victims-{out.name}"
        victims.mkdir(exist_ok=True)
        planted = {}
        for n in sorted(names):
            p = out / n
            if not p.parent.is_dir():
                continue
            v = victims / n.replace("/", "__")
            v.write_bytes(VICTIM)
            if os.path.lexists(p):
                p.unlink()
            p.symlink_to(v)
            planted[n] = v
        return planted

    def check(self, out: Path, planted: dict[str, Path], what: str) -> None:
        written = sorted(n for n, v in planted.items() if v.read_bytes() != VICTIM)
        self.assertEqual(written, [], f"{what}: written through a link")
        links = sorted(str(p.relative_to(out)) for p in out.rglob("*") if p.is_symlink() and p.exists())
        # what is still a link is a name this run did not write (and so did not read either)
        self.assertFalse([n for n in links if n in ("work/outputs.json", "work/state.json", "report.json",
                                                    "dnds.tsv", "failed.tsv", "genes.tsv")], (what, links))

    def twice(self, fn, *args, ok=(0, 2), **kw) -> Path:
        """fn(*args, **kw) (the output folder is its last argument), then twice more over planted links."""
        out = Path(args[-1])
        res, err = run(fn, *args, **kw)
        self.assertIn(res.exit_code, ok, err)
        for resume in (False, True):
            planted = self.plant(out)
            self.assertGreater(len(planted), 20)
            res, err = run(fn, *args, **kw, resume=resume)
            self.assertIn(res.exit_code, ok, (resume, err))
            self.check(out, planted, f"resume={resume}")
        return out


@unittest.skipUnless(os.name == "posix", "symbolic links (POSIX)")
@needs_cpp
class Planted(_Planting):
    def test_pairs_from_given_pairs_with_a_compressed_input(self):
        qgz = self.tmp / "q.fa.gz"
        qgz.write_bytes(gzip.compress(Q.read_bytes()))
        pf = self.tmp / "pairs.tsv"
        pf.write_text("".join(f"{a}\t{b}\n" for a, b in PLAIN[:4]), encoding="utf-8")
        out = self.twice(P.run_pairs, qgz, S, self.tmp / "o", pairs=pf, kaks_format=True, keep_alignments=True,
                         keep_protein_alignments=True)
        self.assertFalse((out / "work" / "inputs").is_symlink())
        self.assertTrue(any(p.is_file() and not p.is_symlink() for p in (out / "work" / "inputs").iterdir()))

    def test_clusters_with_both_alignment_folders(self):
        ids = list(read_fasta(FAM))      # clusters of two (aligned by NW: no Clustal Omega needed)
        tsv = self.tmp / "clusters.tsv"
        tsv.write_text("".join(f"{ids[k]}\t{ids[k]}\n{ids[k]}\t{ids[k + 1]}\n" for k in range(0, 8, 2)),
                       encoding="utf-8")
        out = self.twice(P.run_clusters, FAM, self.tmp / "c", clusters=tsv, keep_alignments=True,
                         keep_protein_alignments=True, kaks_format=True)
        for d in ("codon_alignments", "protein_alignments"):
            files = list((out / d).iterdir())
            self.assertTrue(files and all(p.is_file() and not p.is_symlink() for p in files), d)


@unittest.skipUnless(os.name == "posix", "symbolic links (POSIX)")
@needs_diamond
class PlantedDiamond(_Planting):
    """The routes that run DIAMOND: makedb's databases, blastp's hits, DeepClust's sorted input, realign."""

    def test_pairs_found_by_diamond(self):
        self.twice(P.run_pairs, Q, S, self.tmp / "o", kaks_format=True)

    @unittest.skipUnless(CLUSTALO, "needs Clustal Omega (DeepClust's clusters have more than two members)")
    def test_clusters_by_deepclust_with_prefixed_ids(self):
        self.twice(P.run_clusters, [FAM, S], self.tmp / "c", prefix_ids=True, realign=True, keep_alignments=True)


if __name__ == "__main__":
    unittest.main()
