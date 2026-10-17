#!/usr/bin/env python3
"""Build the simulated test sets (deterministic: same seed -> same files).

    python validation/build_sim.py accuracy            # edge-case families with ground truth (~minutes)
    python validation/build_sim.py speed --max-n 100000
    python validation/build_sim.py speed               # full grid incl. 1,000,000 sequences (large, GBs)
    python validation/build_sim.py families --count 2000

Output under validation/data/sim/ (git-ignored), plus a manifest.jsonl per tier listing every case.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import zlib
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cases import make_case, sha256, write_case  # noqa: E402
from simulate import Family, simulate_clades, simulate_family, translate  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "sim"

# name: (families, n range, L range, sim kwargs, features, note)
ACCURACY = {
    "baseline":        (30, (4, 60), (60, 600), {}, set(), ""),
    "divergent":       (20, (4, 40), (100, 500), {"height": 2.5, "omega": 0.6}, set(), "near-saturated synonymous sites"),
    "indel_heavy":     (20, (4, 60), (100, 600), {"indel_rate": 0.15, "indel_mean": 3.0}, set(), "very gappy alignments"),
    "large_family":    (3, (2000, 2000), (400, 400), {}, set(), ""),
    "long_genes":      (3, (8, 8), (10000, 34350), {}, set(), "titin-scale"),
    "short_genes":     (30, (3, 20), (20, 40), {}, set(), ""),
    "pairs":           (20, (2, 2), (50, 800), {}, set(), ""),
    "utr_stop":        (30, (4, 40), (60, 500), {}, {"utr", "stop3"}, ""),
    "protein_stop":    (15, (4, 30), (60, 400), {}, {"protein_stop"}, ""),
    "internal_stop":   (15, (4, 30), (60, 400), {}, {"internal_stop"}, ""),
    "alt_start_t1":    (8, (4, 30), (60, 400), {}, {"alt_start"}, ""),
    "alt_start_t11":   (8, (4, 30), (60, 400), {"table": 11}, {"alt_start"}, ""),
    "lowercase":       (10, (4, 30), (60, 400), {}, {"lower"}, ""),
    "mixedcase_utr":   (10, (4, 30), (60, 400), {}, {"mixed", "utr"}, ""),
    "rna":             (10, (4, 30), (60, 400), {}, {"rna"}, ""),
    "iupac":           (15, (4, 30), (100, 400), {}, {"iupac"}, "PAL2NAL accepts only some ambiguity codes"),
    "frameshift":      (15, (4, 30), (60, 400), {}, {"frameshift"}, "truth follows the PAL2NAL digit convention"),
    "mismatch":        (15, (4, 30), (200, 600), {}, {"mismatch"}, "expected_warnings.tsv lists every wrong codon"),
    "shuffled_ids":    (10, (4, 40), (60, 400), {}, {"shuffle"}, ""),
    "pair_by_position": (10, (4, 40), (60, 400), {}, {"rename"}, ""),
    "split_nuc":       (10, (4, 40), (60, 400), {}, {"split", "shuffle"}, ""),
    "headers_wrap":    (10, (4, 40), (60, 400), {}, {"headers", "wrap"}, ""),
    "crlf":            (5, (4, 40), (60, 400), {}, {"crlf", "utr"}, "Windows line endings"),
    "kitchen_sink":    (10, (4, 40), (150, 500), {}, {"utr", "stop3", "mixed", "iupac", "frameshift", "shuffle", "wrap"}, ""),
}
for t in (2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23):
    ACCURACY[f"table{t}"] = (5, (4, 40), (60, 400), {"table": t}, {"utr", "stop3"},
                             "PAL2NAL treats TGA as stop in table 10" if t == 10 else "")


def seed_for(*parts) -> int:
    return zlib.crc32("/".join(map(str, parts)).encode())


def log_uniform_int(rng, lo, hi):
    return int(round(np.exp(rng.uniform(np.log(lo), np.log(hi))))) if hi > lo else lo


def build_accuracy(only=None):
    out = DATA / "accuracy"
    manifest = []
    for name, (nfam, nr, lr, kw, feats, note) in ACCURACY.items():
        if only and name not in only:
            continue
        t0 = time.time()
        for f in range(nfam):
            rng = np.random.default_rng(seed_for("accuracy", name, f))
            n = log_uniform_int(rng, *nr)
            L = log_uniform_int(rng, *lr)
            fam = simulate_family(n, L, rng, id_prefix=f"{name}_f{f + 1}_s", **kw)
            case = make_case(fam, rng, set(feats))
            d = out / name / f"fam{f + 1:03d}"
            rec = write_case(d, case, set(feats), rng, {"tier": "accuracy", "scenario": name, "note": note,
                                                        "dir": str(d.relative_to(DATA.parent))})
            manifest.append(rec)
        print(f"{name}: {nfam} families in {time.time() - t0:.1f}s", flush=True)
    _write_manifest(out, manifest)


# ---------------------------------------------------------------- speed tier (streamed writer)
def write_family_fast(d: Path, fam: Family, table: int, meta: dict) -> dict:
    """Write pep/nuc/truth for big families straight from the state matrix (no edits)."""
    d.mkdir(parents=True, exist_ok=True)
    st = fam.states
    cod = np.array([s.encode() for s in st] + [b"---"], dtype="S3")
    aa = np.array([translate(s, table).encode() for s in st] + [b"-"], dtype="S1")
    gap = len(st)
    cds_nt = 0
    with open(d / "pep.fasta", "wb") as pf, open(d / "nuc.fasta", "wb") as nf, open(d / "truth.fasta", "wb") as tf:
        for i, sid in enumerate(fam.ids):
            r = fam.mat[i].astype(np.int64)
            r[r < 0] = gap
            h = b">" + sid.encode() + b"\n"
            pf.write(h + aa[r].tobytes() + b"\n")
            ug = r[r != gap]
            dna = cod[ug].tobytes()
            cds_nt += len(dna)
            nf.write(h + dna + b"\n")
            tf.write(h + cod[r].tobytes() + b"\n")
    rec = {**meta, "pep": "pep.fasta", "nuc": ["nuc.fasta"], "truth": "truth.fasta",
           "opts": [] if table == 1 else ["-codontable", str(table)], "features": [],
           "n_seqs": len(fam.ids), "aln_cols": fam.ncol, "cds_nt": cds_nt, "expected_warnings": 0,
           "bytes_in": (d / "pep.fasta").stat().st_size + (d / "nuc.fasta").stat().st_size}
    (d / "meta.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec


SPEED_GRID = [  # (n_seqs, n_codons)
    (1_000, 400), (10_000, 400), (100_000, 400), (1_000_000, 400),
    (1_000, 100), (1_000, 1_000), (1_000, 10_000), (100, 35_000), (10, 35_000),
]


def build_speed(max_n: int):
    out = DATA / "speed"
    manifest = []
    for n, L in SPEED_GRID:
        if n > max_n:
            continue
        t0 = time.time()
        rng = np.random.default_rng(seed_for("speed", n, L))
        if n <= 2000:
            fam = simulate_family(n, L, rng, height=0.6)
        else:
            fam = simulate_clades(n, L, rng, n_clades=512, height=0.6)
        d = out / f"n{n}_L{L}"
        rec = write_family_fast(d, fam, 1, {"tier": "speed", "scenario": f"n{n}_L{L}", "note": "",
                                            "dir": str(d.relative_to(DATA.parent))})
        rec["sha256_nuc"] = sha256(d / "nuc.fasta")
        manifest.append(rec)
        print(f"n={n} L={L}: {rec['bytes_in'] / 1e6:.0f} MB input, {time.time() - t0:.1f}s", flush=True)
    _write_manifest(out, manifest)


def build_families(count: int):
    """Many realistic gene families (sizes/lengths roughly like real orthogroups) for per-family throughput."""
    out = DATA / "families"
    manifest = []
    t0 = time.time()
    for f in range(count):
        rng = np.random.default_rng(seed_for("families", f))
        n = min(500, max(2, int(rng.lognormal(np.log(25), 1.0))))
        L = min(5000, max(30, int(rng.lognormal(np.log(400), 0.6))))
        fam = simulate_family(n, L, rng, height=float(rng.uniform(0.2, 1.5)), id_prefix=f"fam{f + 1}_s")
        d = out / f"fam{f + 1:05d}"
        manifest.append(write_family_fast(d, fam, 1, {"tier": "families", "scenario": "families", "note": "",
                                                      "dir": str(d.relative_to(DATA.parent))}))
        if (f + 1) % 500 == 0:
            print(f"{f + 1} families, {time.time() - t0:.0f}s", flush=True)
    _write_manifest(out, manifest)


def _write_manifest(out: Path, manifest: list[dict]):
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "manifest.jsonl", "w") as fh:
        for r in manifest:
            fh.write(json.dumps(r) + "\n")
    print(f"wrote {out / 'manifest.jsonl'} ({len(manifest)} cases)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="tier", required=True)
    a = sub.add_parser("accuracy")
    a.add_argument("--only", nargs="*")
    s = sub.add_parser("speed")
    s.add_argument("--max-n", type=int, default=1_000_000)
    f = sub.add_parser("families")
    f.add_argument("--count", type=int, default=2000)
    args = ap.parse_args()
    if args.tier == "accuracy":
        build_accuracy(args.only)
    elif args.tier == "speed":
        build_speed(args.max_n)
    else:
        build_families(args.count)


if __name__ == "__main__":
    main()
