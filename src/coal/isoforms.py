"""Several transcripts of one gene in a CDS file: finding them, and keeping one per gene (docs/pipeline.md
"One record is one gene").

`coal pipeline` treats every FASTA record as a gene. Genome annotations (Ensembl, NCBI, TAIR, Phytozome, Trinity
assemblies) list every transcript, so a gene with five isoforms is five "genes": it is compared five times, and
reciprocal best hits are split between isoforms. This module reads the gene of each record from

  * a gene map (--gene-map TSV: record ID, gene ID), or
  * a field of the FASTA header: `gene:ENSG...` (Ensembl), `[gene=...]` / `gene=...` (NCBI, GFF-derived files),
    `[locus_tag=...]` / `locus_tag=...`, `gene_symbol:...`; the first of these that at least 90% of the records of
    a file have, or
  * the ID itself, when IDs share a stem and differ in a transcript suffix: `AT1G01010.1` / `AT1G01010.2` (TAIR,
    Ensembl Plants, Phytozome), `g12.t1` / `g12.t2` (AUGUSTUS, BRAKER), `..._g1_i1` / `..._g1_i2` (Trinity).

Detection is conservative: a suffix rule applies only when some stems are shared and a stem has four records or
fewer on average (`scaffold12.1 ... scaffold12.340` is a numbering of genes, not of isoforms), and an ID that is
alone with its stem (the version of `ENST00000380152.8`, `NM_000059.4`) is its own gene either way.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

__all__ = ["read_headers", "detect", "read_gene_map", "choose_longest", "Detection", "HEADER_FIELDS"]

# (name shown in messages, pattern on the header text after the ID), in order of preference
HEADER_FIELDS: tuple[tuple[str, re.Pattern], ...] = (
    ("gene:", re.compile(r"(?:^|\s)gene:(\S+)")),
    ("[gene=", re.compile(r"\[gene=([^\]\s][^\]]*)\]")),
    ("gene=", re.compile(r"(?:^|[\s;|])gene=([^\s;|\]]+)")),
    ("[locus_tag=", re.compile(r"\[locus_tag=([^\]\s][^\]]*)\]")),
    ("locus_tag=", re.compile(r"(?:^|[\s;|])locus_tag=([^\s;|\]]+)")),
    ("gene_symbol:", re.compile(r"(?:^|\s)gene_symbol:(\S+)")),
)
# (name, pattern on the ID -> the stem), tried in this order
SUFFIX_RULES: tuple[tuple[str, re.Pattern], ...] = (
    ("ID suffix _iN (Trinity)", re.compile(r"^(.+_g\d+)_i\d+$")),
    ("ID suffix .tN", re.compile(r"^(.+)\.t\d+$")),
    ("ID suffix .N", re.compile(r"^(.+)\.\d+$")),
)
FIELD_SHARE = 0.9          # a header field counts when this share of a file's records has it
MAX_MEAN_GROUP = 4.0       # a suffix rule counts when a stem has at most this many records on average


@dataclass
class Detection:
    source: str | None = None          # "gene map", "header field gene:", "ID suffix .N", or None: nothing found
    gene_of: dict = field(default_factory=dict)    # record ID -> gene (only records whose gene is known)
    records: int = 0
    genes: int = 0                     # distinct genes among the records (a record without a gene is its own)
    multi: int = 0                     # genes with more than one record
    extra: int = 0                     # records beyond the first of their gene (what --longest-isoform drops)
    example: tuple | None = None       # (gene, [record IDs]) of one gene with several records


def read_headers(paths) -> list[tuple[str, str]]:
    """[(record ID, the rest of its header line)] of FASTA files, in order (bytes that are not UTF-8 as \\xhh, as
    the pipeline names IDs)."""
    out = []
    for path in paths:
        with open(path, "rb") as fh:
            for line in fh:
                if line[:1] == b">":
                    text = line[1:].decode("utf-8", "backslashreplace").strip()
                    words = text.split(None, 1)
                    rid = words[0] if words else ""
                    rest = words[1] if len(words) > 1 else ""
                    out.append((rid, rest))
    return out


def _summarise(det: Detection, ids) -> Detection:
    groups: dict[str, list[str]] = {}
    seen = set()
    for rid in ids:
        if rid in seen:   # a repeated ID: its first record stands
            continue
        seen.add(rid)
        groups.setdefault(det.gene_of.get(rid, "\0" + rid), []).append(rid)
    det.records = len(seen)
    det.genes = len(groups)
    many = [(g, m) for g, m in groups.items() if len(m) > 1]
    det.multi = len(many)
    det.extra = sum(len(m) - 1 for _, m in many)
    det.example = (many[0][0], many[0][1][:4]) if many else None
    return det


def detect(headers) -> Detection:
    """The gene of each record from the headers [(id, rest)] of one set of files: a header field, else a shared
    ID stem, else nothing (source None)."""
    ids = [rid for rid, _ in headers]
    n = len(headers)
    if n == 0:
        return Detection()
    for name, pat in HEADER_FIELDS:
        gene_of = {}
        for rid, rest in headers:
            m = pat.search(rest)
            if m and rid not in gene_of:
                gene_of[rid] = m.group(1).strip()
        if len(gene_of) >= FIELD_SHARE * n:
            return _summarise(Detection(source=f"header field {name}", gene_of=gene_of), ids)
    for name, pat in SUFFIX_RULES:
        stems: dict[str, list[str]] = {}
        for rid in dict.fromkeys(ids):
            m = pat.match(rid)
            if m:
                stems.setdefault(m.group(1), []).append(rid)
        shared = {s: m for s, m in stems.items() if len(m) > 1}
        if not shared:
            continue
        with_suffix = sum(len(m) for m in stems.values())
        if with_suffix / len(stems) > MAX_MEAN_GROUP:
            continue   # a numbering of genes (scaffold7.1 ... scaffold7.340), not of transcripts
        gene_of = {rid: s for s, m in shared.items() for rid in m}
        return _summarise(Detection(source=name, gene_of=gene_of), ids)
    return _summarise(Detection(), ids)


def read_gene_map(path: str, shown: str | None = None) -> dict[str, str]:
    """--gene-map: record ID and gene ID in the first two columns (tab-separated, or blank-separated when a line
    has no tab); `#` lines, empty lines and a header line (id / record / transcript ... gene) are skipped. A record
    named twice keeps its first gene. ValueError for a line with one column."""
    shown = shown or path
    out: dict[str, str] = {}
    with open(path, encoding="utf-8", errors="backslashreplace") as fh:
        for n, line in enumerate(fh):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t") if "\t" in line else line.split()
            if len(f) < 2 or not f[0].strip() or not f[1].strip():
                raise ValueError(f"{shown}: line {n + 1}: expected two columns (record ID, gene ID)")
            a, b = f[0].strip(), f[1].strip()
            if not out and a.lower() in ("id", "record", "record_id", "transcript", "transcript_id", "seq", "seqid",
                                         "protein", "protein_id", "cds") and "gene" in b.lower():
                continue
            out.setdefault(a, b)
    return out


def from_gene_map(gene_map: dict, ids) -> Detection:
    ids = list(ids)
    return _summarise(Detection(source="gene map", gene_of={r: gene_map[r] for r in ids if r in gene_map}), ids)


def choose_longest(det: Detection, lengths: dict) -> dict:
    """{record ID: (gene, the record kept)} of the records to drop: of each gene with several usable records
    (those in `lengths`: ID -> CDS length in nucleotides) the longest stays, and of equally long ones the one whose
    ID sorts first (bytes)."""
    groups: dict[str, list[str]] = {}
    for rid in lengths:
        g = det.gene_of.get(rid)
        if g is not None:
            groups.setdefault(g, []).append(rid)
    dropped = {}
    for g, members in groups.items():
        if len(members) < 2:
            continue
        keep = min(members, key=lambda r: (-lengths[r], r.encode("utf-8", "backslashreplace")))
        for r in members:
            if r != keep:
                dropped[r] = (g, keep)
    return dropped
