#!/usr/bin/env python3
"""Build real tree-of-life test cases from NCBI RefSeq genomes (no reference answer; scored by
the reference-free checks and against PAL2NAL).

For a spread of RefSeq reference genomes (bacteria, archaea, fungi, protists, plants, animals,
including organisms with alternative genetic codes), download protein.faa + cds_from_genomic.fna,
pair each protein with its CDS by protein_id, group proteins into families by gene symbol (longest
isoform per genome), align each family's proteins with FAMSA (pyfamsa), and write
pep.fasta (alignment) + nuc.fasta (the raw annotated CDS, stop codon included).

These are the real-world quirks COAL will meet at scale: stop codons, partial CDS, selenocysteine
(transl_except), programmed frameshifts / ribosomal slippage (protein length != CDS/3), IUPAC bases,
extreme GC, and genetic codes 4, 6, 11, 12. Each sequence's quirks are recorded in meta.json.
The genetic code of each genome is inferred from the data (the table under which its CDSs
translate to its proteins), and families are built per genetic code.

    python validation/build_ncbi.py                       # ~300 genomes, 3,000 families
    python validation/build_ncbi.py --bacteria 1000 --families 20000
"""
from __future__ import annotations

import argparse
import gzip
import json
import random
import re
import subprocess
import sys
import time
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from simulate import translate  # noqa: E402

RAW = HERE / "data" / "raw" / "ncbi"
OUT = HERE / "data" / "real" / "ncbi"
SUMMARY = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{g}/assembly_summary.txt"
GROUPS = {"bacteria": 250, "archaea": 30, "fungi": 10, "protozoa": 8, "plant": 5, "invertebrate": 6,
          "vertebrate_mammalian": 3, "vertebrate_other": 4}
# Always included when present as reference genomes (alternative codes, extreme GC, model organisms)
PINNED = ["Mycoplasmoides genitalium", "Mycoplasmoides pneumoniae", "Mesomycoplasma hyopneumoniae",
          "Tetrahymena thermophila", "Paramecium tetraurelia", "Candida albicans", "Escherichia coli",
          "Bacillus subtilis", "Streptomyces coelicolor", "Plasmodium falciparum", "Saccharomyces cerevisiae",
          "Arabidopsis thaliana", "Drosophila melanogaster", "Homo sapiens", "Methanocaldococcus jannaschii",
          "Danio rerio", "Caenorhabditis elegans"]
CANDIDATE_TABLES = [11, 1, 4, 6, 12, 25]
HDR = re.compile(r"\[(\w+)=([^\]]*)\]")


def fetch(url: str, dest: Path, tries=3):
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    for k in range(tries):
        try:
            tmp = dest.with_suffix(dest.suffix + ".part")
            with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as fh:
                while b := r.read(1 << 20):
                    fh.write(b)
            tmp.rename(dest)
            return dest
        except Exception as e:  # noqa: BLE001
            if k == tries - 1:
                print(f"  failed {url}: {e}", flush=True)
                return None
            time.sleep(2)


def choose_genomes(rng, counts) -> list[dict]:
    chosen = []
    for g, k in counts.items():
        p = fetch(SUMMARY.format(g=g), RAW / f"assembly_summary_{g}.txt")
        rows = []
        for line in open(p, encoding="utf-8", errors="replace"):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 28 or f[4] != "reference genome" or f[19] == "na":
                continue
            if g in ("bacteria", "archaea") and f[11] not in ("Complete Genome", "Chromosome"):
                continue
            rows.append({"group": g, "acc": f[0], "organism": f[7], "ftp": f[19].replace("ftp://", "https://"),
                         "gc": float(f[27]) if f[27] not in ("", "na") else None})
        pinned = [r for r in rows if any(r["organism"].startswith(p) for p in PINNED)]
        seen_sp = {r["organism"].split()[0:2].__str__() for r in pinned}
        pool = [r for r in rows if r not in pinned and str(r["organism"].split()[0:2]) not in seen_sp]
        # spread over GC content: sort by GC and take evenly spaced picks from a shuffled-within-stratum list
        pool.sort(key=lambda r: (r["gc"] or 50, rng.random()))
        need = max(0, k - len(pinned))
        step = len(pool) / need if need else 0
        picks = [pool[int(i * step)] for i in range(need)] if need and pool else []
        chosen += pinned[: max(k, len(pinned))] + picks
        print(f"{g}: {len(rows)} reference genomes, chose {len(pinned)} pinned + {len(picks)}", flush=True)
    return chosen


def download_genome(gen: dict) -> dict | None:
    base = gen["ftp"].rstrip("/")
    stem = base.split("/")[-1]
    d = RAW / "genomes" / stem
    prot = fetch(f"{base}/{stem}_protein.faa.gz", d / f"{stem}_protein.faa.gz")
    cds = fetch(f"{base}/{stem}_cds_from_genomic.fna.gz", d / f"{stem}_cds_from_genomic.fna.gz")
    if not prot or not cds:
        return None
    return {**gen, "prot": prot, "cds": cds, "stem": stem}


def read_gz_fasta(path: Path):
    hdr, seq = None, []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(seq)
                hdr, seq = line[1:].rstrip("\n"), []
            else:
                seq.append(line.strip())
    if hdr is not None:
        yield hdr, "".join(seq)


_DROP = re.compile(r"^(low quality protein:|putative|probable|predicted:?)\s+|\s+(isoform|transcript variant)\b.*$|,? partial$|-like\b")


def _norm_product(p: str) -> str:
    p = p.strip().lower()
    for _ in range(3):
        p = _DROP.sub("", p).strip()
    if not p or "hypothetical" in p or "uncharacterized" in p or "unknown" in p or "domain-containing" in p:
        return ""
    return "PRODUCT:" + p


def load_genome(gen: dict, max_check=300):
    prots = {h.split()[0]: s for h, s in read_gz_fasta(gen["prot"])}
    recs = []
    for h, s in read_gz_fasta(gen["cds"]):
        tags = dict(HDR.findall(h))
        pid = tags.get("protein_id")
        if not pid or pid not in prots or tags.get("pseudo") == "true":
            continue
        sym = (tags.get("gene") or "").upper()
        if not sym or sym.startswith(("LOC", "ORF")) or re.fullmatch(r"[A-Z]+_?\d{3,}", sym):
            sym = ""  # locus-tag-like or missing: fall back to the product name below
        prod = _norm_product(tags.get("protein", ""))
        quirks = [k for k in ("partial", "transl_except", "exception") if k in tags]
        recs.append({"pid": pid, "sym": sym or prod, "cds": s, "prot": prots[pid], "quirks": quirks,
                     "quirk_detail": {k: tags[k] for k in quirks}})
    # infer genetic code: which table makes the most CDSs translate to their protein
    sample = [r for r in recs if not r["quirks"]][:max_check]
    prok = gen["group"] in ("bacteria", "archaea")
    best, best_score = (11 if prok else 1), -1
    for t in ([11, 4, 25] if prok else [1, 6, 12, 4, 10, 15, 26]):
        ok = 0
        for r in sample:
            c = r["cds"].upper()
            aa = "".join(translate(c[i:i + 3], t) for i in range(0, len(c) - len(c) % 3, 3))
            ok += aa[1:len(r["prot"])] == r["prot"][1:]
        if ok > best_score:
            best, best_score = t, ok
    gen["table"] = best
    gen["table_fit"] = round(best_score / max(1, len(sample)), 3)
    # longest protein per family key
    by = {}
    for r in recs:
        if not r["sym"]:
            continue
        if r["sym"] not in by or len(r["prot"]) > len(by[r["sym"]]["prot"]):
            by[r["sym"]] = r
    return by


def align(prots: list[tuple[str, str]], threads: int) -> list[tuple[str, str]]:
    from pyfamsa import Aligner, Sequence
    al = Aligner(threads=threads, guide_tree="upgma" if len(prots) < 1000 else "sl")
    std = set("ACDEFGHIKLMNPQRSTVWY")
    msa = al.align([Sequence(i.encode(), "".join(c if c in std else "X" for c in s.upper()).encode())
                    for i, s in prots])
    out = {s.id.decode(): s.sequence.decode() for s in msa}
    rows = []
    for i, s in prots:
        it = iter(s)  # put the original residues back (U, *, B, Z, ... were aligned as X)
        rows.append((i, "".join(next(it) if ch != "-" else "-" for ch in out[i])))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    for g, k in GROUPS.items():
        ap.add_argument(f"--{g.replace('_', '-')}", type=int, default=k, dest=g)
    ap.add_argument("--families", type=int, default=3000)
    ap.add_argument("--min-members", type=int, default=4)
    ap.add_argument("--threads", type=int, default=0, help="FAMSA threads (0 = all)")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    gens = choose_genomes(rng, {g: getattr(a, g) for g in GROUPS})
    with ThreadPoolExecutor(8) as ex:
        gens = [g for g in ex.map(download_genome, gens) if g]
    print(f"downloaded {len(gens)} genomes", flush=True)
    fams: dict[tuple, list] = defaultdict(list)
    genome_info = []
    for k, g in enumerate(gens):
        by = load_genome(g)
        genome_info.append({k2: g[k2] for k2 in ("group", "acc", "organism", "gc", "table", "table_fit")})
        tag = re.sub(r"\W+", "_", "_".join(g["organism"].split()[:2]))
        for sym, r in by.items():
            fams[(sym, g["table"])].append({**r, "genome": g["acc"], "id": f"{tag}|{r['pid']}"})
        if (k + 1) % 25 == 0:
            print(f"  loaded {k + 1}/{len(gens)} genomes", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "genomes.json").write_text(json.dumps(genome_info, indent=1))
    print("genetic codes inferred:", Counter(g["table"] for g in genome_info), flush=True)
    eligible = [k for k, v in fams.items() if len(v) >= a.min_members]
    # stratify by family size so big families are represented
    bins = defaultdict(list)
    for k in eligible:
        n = len(fams[k])
        bins[0 if n <= 10 else 1 if n <= 50 else 2].append(k)
    per = a.families // 3
    chosen = []
    for b in sorted(bins):
        rng.shuffle(bins[b])
        chosen += bins[b][:per]
    rest = [k for b in sorted(bins) for k in bins[b][per:]]
    rng.shuffle(rest)
    chosen += rest[: a.families - len(chosen)]
    # always keep every family on a non-standard code
    chosen += [k for k in eligible if k[1] not in (1, 11) and k not in set(chosen)][:500]
    print(f"{len(eligible)} families with >= {a.min_members} members; writing {len(chosen)}", flush=True)
    manifest = []
    t0 = time.time()
    for i, (sym, table) in enumerate(chosen):
        members = fams[(sym, table)]
        aln = align([(m["id"], m["prot"]) for m in members], a.threads)
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", sym)[:80]
        d = OUT / "families" / f"t{table}" / f"{safe}"
        d.mkdir(parents=True, exist_ok=True)
        with open(d / "pep.fasta", "w", newline="\n") as fh:
            for sid, s in aln:
                fh.write(f">{sid}\n")
                for k in range(0, len(s), 60):
                    fh.write(s[k:k + 60] + "\n")
        with open(d / "nuc.fasta", "w", newline="\n") as fh:
            for m in members:
                fh.write(f">{m['id']}\n{m['cds']}\n")
        quirks = Counter(q for m in members for q in m["quirks"])
        feats = sorted(quirks) + (["selenocysteine"] if any("U" in m["prot"] for m in members) else []) + \
            (["cds_len_not_3x_protein"] if any(len(m["cds"]) not in (3 * len(m["prot"]), 3 * len(m["prot"]) + 3)
                                                for m in members) else [])
        rec = {"tier": "ncbi", "scenario": f"table{table}_{'small' if len(members) <= 10 else 'medium' if len(members) <= 50 else 'large'}",
               "note": f"gene symbol {sym}", "dir": str(d.relative_to(HERE / "data")), "pep": "pep.fasta",
               "nuc": ["nuc.fasta"], "truth": None, "opts": [] if table == 1 else ["-codontable", str(table)],
               "features": sorted(set(feats)), "n_seqs": len(members), "aln_cols": len(aln[0][1]),
               "cds_nt": sum(len(m["cds"]) for m in members), "expected_warnings": 0,
               "genomes": sorted({m["genome"] for m in members}),
               "quirky_members": {m["id"]: m["quirk_detail"] for m in members if m["quirks"]}}
        (d / "meta.json").write_text(json.dumps(rec, indent=1) + "\n")
        manifest.append(rec)
        if (i + 1) % 250 == 0:
            print(f"  aligned {i + 1}/{len(chosen)} families ({time.time() - t0:.0f}s)", flush=True)
    with open(OUT / "manifest.jsonl", "w") as fh:
        for r in manifest:
            fh.write(json.dumps(r) + "\n")
    print(f"wrote {len(manifest)} families to {OUT}")


if __name__ == "__main__":
    main()
