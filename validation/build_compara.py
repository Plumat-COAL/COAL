#!/usr/bin/env python3
"""Build real-data test cases from Ensembl Compara gene trees (independent reference answer).

Ensembl publishes, for every protein gene tree, the protein alignment (aa.fasta.gz) and the CDS
alignment Ensembl derived from it (cds.fasta.gz); alignments are separated by '//' lines.
For each tree we give COAL the protein alignment + the ungapped CDS, and Ensembl's CDS alignment
is the expected answer ("truth").

    python validation/build_compara.py --download            # fetch release files (~3.3 GB) into validation/data/raw/
    python validation/build_compara.py                       # stratified sample (default 400 trees per size bin)
    python validation/build_compara.py --per-bin 2000 --special 500
    python validation/build_compara.py --all                 # every tree (large: tens of GB)

Size bins (sequences per tree): 2-10, 11-50, 51-200, 201-1000, >1000. On top of the random sample,
up to --special trees of each kind are kept: selenocysteine (U), X residues, CDS length not
3 x residues (partial codons), other unusual characters. A stats file for EVERY tree
(stats.tsv) records size and flags, so the sample can be compared with the whole release.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import random
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / "data" / "raw" / "compara"
OUT = HERE / "data" / "real" / "compara"
BASE = "https://ftp.ensembl.org/pub/release-{rel}/emf/ensembl-compara/homologies/"
FILES = ["Compara.{rel}.protein_default.aa.fasta.gz", "Compara.{rel}.protein_default.cds.fasta.gz"]
BINS = [(2, 10), (11, 50), (51, 200), (201, 1000), (1001, 10**9)]
STD_AA = set("ACDEFGHIKLMNPQRSTVWY-")


def download(rel: int, raw: Path):
    raw.mkdir(parents=True, exist_ok=True)
    base = BASE.format(rel=rel)
    md5 = urllib.request.urlopen(base + "MD5SUM").read().decode()
    sums = {l.split()[1]: l.split()[0] for l in md5.splitlines() if len(l.split()) == 2}
    for f in FILES:
        name = f.format(rel=rel)
        dest = raw / name
        if not dest.exists():
            print(f"downloading {name} ...", flush=True)
            if shutil.which("curl"):
                subprocess.run(["curl", "-sS", "-C", "-", "-o", str(dest) + ".part", base + name], check=True)
                Path(str(dest) + ".part").rename(dest)
            else:
                urllib.request.urlretrieve(base + name, dest)
        h = hashlib.md5()
        with open(dest, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 22), b""):
                h.update(b)
        ok = sums.get(name) == h.hexdigest()
        print(f"{name}: md5 {'OK' if ok else 'MISMATCH'}", flush=True)
        if not ok:
            sys.exit("checksum mismatch; delete the file and retry")


def _open(path: Path):
    for tool in ("pigz", "gzip"):
        if shutil.which(tool):
            p = subprocess.Popen([tool, "-dc", str(path)], stdout=subprocess.PIPE, bufsize=1 << 20)
            return p.stdout
    import gzip
    return gzip.open(path, "rb")


def blocks(path: Path):
    """Yield each alignment as [(id, seq)]."""
    recs: list[list] = []
    for raw in _open(path):
        line = raw.rstrip(b"\r\n")
        if line == b"//":
            if recs:
                yield [(i, b"".join(s).decode()) for i, s in recs]
            recs = []
        elif line.startswith(b">"):
            recs.append((line[1:].split()[0].decode(), []))
        elif line and recs:
            recs[-1][1].append(line.strip())
    if recs:
        yield [(i, b"".join(s).decode()) for i, s in recs]


def flags_for(aa: list[tuple[str, str]], cds: list[tuple[str, str]]) -> set[str]:
    f = set()
    chars = set("".join(s for _, s in aa).upper())
    if "U" in chars:
        f.add("selenocysteine")
    if "X" in chars:
        f.add("x_residue")
    if "*" in chars:
        f.add("stop_in_protein")
    if chars - STD_AA - {"U", "X", "*"}:
        f.add("other_chars")
    for (_, a), (_, c) in zip(aa, cds):
        nres = sum(1 for ch in a if ch != "-")
        nnt = sum(1 for ch in c if ch != "-")
        if nnt != 3 * nres:
            f.add("cds_len_not_3x")
            break
    for _, c in cds:
        if set(c.upper()) - set("ACGTN-"):
            f.add("ambiguous_bases")
            break
    return f


def bin_of(n: int) -> int:
    for k, (lo, hi) in enumerate(BINS):
        if lo <= n <= hi:
            return k
    return 0


def write_tree(d: Path, idx: int, aa, cds, flags, rel: int) -> dict:
    d.mkdir(parents=True, exist_ok=True)
    with open(d / "pep.fasta", "w", newline="\n") as fh:
        for i, s in aa:
            fh.write(f">{i}\n")
            for k in range(0, len(s), 60):
                fh.write(s[k:k + 60] + "\n")
    with open(d / "nuc.fasta", "w", newline="\n") as fh:
        for i, s in cds:
            fh.write(f">{i}\n{s.replace('-', '')}\n")
    with open(d / "truth.fasta", "w", newline="\n") as fh:
        for i, s in cds:
            fh.write(f">{i}\n{s}\n")
    rec = {"tier": "compara", "scenario": f"bin{bin_of(len(aa))}_{BINS[bin_of(len(aa))][0]}-"
                                          f"{min(BINS[bin_of(len(aa))][1], 99999)}",
           "note": f"Ensembl Compara {rel} tree #{idx} (first member {aa[0][0]})",
           "dir": str(d.relative_to(HERE / "data")), "pep": "pep.fasta", "nuc": ["nuc.fasta"],
           "truth": "truth.fasta", "opts": [], "features": sorted(flags), "n_seqs": len(aa),
           "aln_cols": len(aa[0][1]), "cds_nt": sum(len(s) - s.count("-") for _, s in cds),
           "expected_warnings": 0}
    (d / "meta.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--release", type=int, default=116)
    ap.add_argument("--raw", type=Path, default=RAW)
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--per-bin", type=int, default=400)
    ap.add_argument("--special", type=int, default=200, help="max extra trees per unusual feature")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--seed", type=int, default=116)
    a = ap.parse_args()
    if a.download:
        download(a.release, a.raw)
        return
    aa_path, cds_path = (a.raw / f.format(rel=a.release) for f in FILES)
    for p in (aa_path, cds_path):
        if not p.exists():
            sys.exit(f"{p} missing; run with --download first")
    rng = random.Random(a.seed)
    OUT.mkdir(parents=True, exist_ok=True)
    stats_path = OUT / "stats.tsv"
    # pass 1: stats for every tree (sizes + flags); cached, so re-sampling is quick
    if not stats_path.exists():
        bad_pairs = 0
        with open(str(stats_path) + ".part", "w") as st:
            st.write("tree\tfirst_id\tn_seqs\taln_cols\tcds_nt\tflags\n")
            for idx, (aa, cds) in enumerate(zip(blocks(aa_path), blocks(cds_path)), 1):
                if [i for i, _ in aa] != [i for i, _ in cds]:
                    bad_pairs += 1
                    st.write(f"{idx}\t{aa[0][0]}\t{len(aa)}\t{len(aa[0][1])}\t0\tid_order_mismatch\n")
                    continue
                fl = flags_for(aa, cds)
                st.write(f"{idx}\t{aa[0][0]}\t{len(aa)}\t{len(aa[0][1])}\t"
                         f"{sum(len(s) - s.count('-') for _, s in cds)}\t{','.join(sorted(fl))}\n")
                if idx % 2000 == 0:
                    print(f"pass 1: {idx} trees", flush=True)
        Path(str(stats_path) + ".part").rename(stats_path)
        print(f"pass 1 done; aa/cds ID-order mismatches: {bad_pairs}")
    rows = [l.rstrip("\n").split("\t") for l in open(stats_path).readlines()[1:]]
    trees = [(int(r[0]), int(r[2]), set(filter(None, r[5].split(",")))) for r in rows]
    trees = [t for t in trees if "id_order_mismatch" not in t[2]]
    if a.all:
        chosen = {t[0] for t in trees}
    else:
        chosen = set()
        for b in range(len(BINS)):
            pool = [t[0] for t in trees if bin_of(t[1]) == b]
            chosen |= set(rng.sample(pool, min(a.per_bin, len(pool))))
        kinds = sorted({f for t in trees for f in t[2]})
        for f in kinds:
            pool = [t[0] for t in trees if f in t[2]]
            chosen |= set(rng.sample(pool, min(a.special, len(pool))))
    counts = {f"{lo}-{hi}": sum(1 for t in trees if bin_of(t[1]) == k) for k, (lo, hi) in enumerate(BINS)}
    feats = {}
    for t in trees:
        for f in t[2]:
            feats[f] = feats.get(f, 0) + 1
    print(f"{len(trees)} trees; per size bin: {counts}; unusual: {feats}; writing {len(chosen)}", flush=True)
    # pass 2: write the chosen trees
    sub = "all" if a.all else "sample"
    manifest = []
    for idx, (aa, cds) in enumerate(zip(blocks(aa_path), blocks(cds_path)), 1):
        if idx in chosen:
            fl = flags_for(aa, cds)
            manifest.append(write_tree(OUT / sub / f"t{idx:06d}", idx, aa, cds, fl, a.release))
    with open(OUT / "manifest.jsonl", "w") as fh:
        for r in manifest:
            fh.write(json.dumps(r) + "\n")
    print(f"wrote {len(manifest)} cases to {OUT}")


if __name__ == "__main__":
    main()
