"""Several transcripts of one gene in a CDS file: finding them, and keeping one per gene (docs/pipeline.md
"One record is one gene").

`coal pipeline` treats every FASTA record as a gene. Genome annotations (Ensembl, NCBI, TAIR, Phytozome, FlyBase,
MAKER, Trinity assemblies) list every transcript, so a gene with five isoforms is five "genes": it is compared five
times, and reciprocal best hits are split between isoforms. This module reads the gene of each record from

  * a gene map (--gene-map TSV: record ID, gene ID): it wins for every record it names; the records it does not
    name fall back to the next two, or
  * a field of the FASTA header, stable gene identifiers before gene symbols: `[locus_tag=...]` / `locus_tag=...`,
    `GeneID:...` (also as `[db_xref=GeneID:...]`), `gene:...` (Ensembl), `parent=FBgn...` (FlyBase), `locus=...`
    (Phytozome), then `[gene=...]` / `gene=...` (NCBI, GFF-derived files) and `gene_symbol:...`; the first of
    these that at least 90% of the records of a file have, or
  * the ID itself, when IDs share a stem and differ in a transcript suffix (SUFFIX_RULES): `..._g1_i1.p1`
    (TransDecoder on Trinity), `..._g1_i1` (Trinity), `g12.t1` (AUGUSTUS, BRAKER), `-mRNA-1` and `-RA` / `-RB`
    (MAKER), `_T001` / `_P001` (maize), `.1.p` (Phytozome proteins), `AT1G01010.1` (TAIR, Ensembl Plants,
    Phytozome), `Os01t0100100-01` (RAP-DB).

Detection is conservative, and a heuristic all the same:

  * a suffix rule applies only when some stems are shared and a stem has four records or fewer on average
    (`scaffold12.1 ... scaffold12.340` is a numbering of genes, not of isoforms), and an ID that is alone with its
    stem (the version of `ENST00000380152.8`, `NM_000059.4`) is its own gene either way. It cannot tell a
    numbering of paralogs (`HSP70.1`, `HSP70.2`) or of a few genes per scaffold (`scaffold7.1 ... .3`) from
    transcripts: the pipeline refuses --longest-isoform when a suffix rule alone would drop more than half of a
    set's records (unless --isoform-rule suffix says the suffixes are transcripts);
  * a gene symbol (`[gene=dnaA]`) names different loci in a file of several strains or genomes, and two copies of
    a gene in one genome (`tuf`): records that share a symbol but sit on different contigs, or carry different
    locus tags or GeneIDs, are not merged (Detection.conflicts counts the symbols; the caller warns);
  * RefSeq transcript accessions (NM_, XM_, NR_, XR_) hold no gene: without a gene field in the header nothing
    can be found (Detection.refseq_ids says so).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

__all__ = ["read_headers", "detect", "read_gene_map", "from_gene_map", "choose_longest", "Detection",
           "HEADER_FIELDS", "SUFFIX_RULES", "STYLES"]

# (name shown in messages, pattern on the header text after the ID, True for a stable identifier / False for a
# gene symbol), in order of preference
HEADER_FIELDS: tuple[tuple[str, re.Pattern, bool], ...] = (
    ("[locus_tag=", re.compile(r"\[locus_tag=([^\]\s][^\]]*)\]"), True),
    ("locus_tag=", re.compile(r"(?:^|[\s;|])locus_tag=([^\s;|\]]+)"), True),
    ("GeneID:", re.compile(r"(?:^|[\s\[=,;|])GeneID:(\d+)"), True),
    ("gene:", re.compile(r"(?:^|\s)gene:(\S+)"), True),
    ("parent=FBgn", re.compile(r"(?:^|[\s;])parent=(FBgn\d+)"), True),
    ("locus=", re.compile(r"(?:^|[\s;|])locus=([^\s;|\]]+)"), True),
    ("[gene=", re.compile(r"\[gene=([^\]\s][^\]]*)\]"), False),
    ("gene=", re.compile(r"(?:^|[\s;|])gene=([^\s;|\]]+)"), False),
    ("gene_symbol:", re.compile(r"(?:^|\s)gene_symbol:(\S+)"), False),
)
# (name, pattern on the ID -> the stem), tried in this order
SUFFIX_RULES: tuple[tuple[str, re.Pattern], ...] = (
    ("ID suffix _iN.pN (TransDecoder on Trinity)", re.compile(r"^(.+_g\d+)_i\d+\.p\d+$")),
    ("ID suffix _iN (Trinity)", re.compile(r"^(.+_g\d+)_i\d+$")),
    ("ID suffix .tN", re.compile(r"^(.+)\.t\d+$")),
    ("ID suffix -mRNA-N (MAKER)", re.compile(r"^(.+)-mRNA-\d+$")),
    ("ID suffix -RA (MAKER)", re.compile(r"^(.+)-R[A-Z]{1,2}$")),
    ("ID suffix _TNNN (maize)", re.compile(r"^(.+)_[TP]\d{3,}$")),
    ("ID suffix .N.p (Phytozome)", re.compile(r"^(.+)\.\d+\.p$")),
    ("ID suffix .N", re.compile(r"^(.+)\.\d+$")),
    ("ID suffix -NN (RAP-DB)", re.compile(r"^(.+)-\d{2}$")),
)
# The table of docs/pipeline.md "Header styles" (tests/test_pipeline_isoforms.py runs every example):
# (annotation style, two records of one gene as header lines, the source detect() names)
STYLES: tuple[tuple[str, tuple[str, str], str], ...] = (
    ("NCBI, with locus tags", ("lcl|NC_1_cds_NP_1.1_1 [gene=thrA] [locus_tag=b0002] [protein=x]",
                               "lcl|NC_1_cds_NP_2.1_2 [gene=thrA] [locus_tag=b0002]"), "header field [locus_tag="),
    ("NCBI, with GeneIDs", ("lcl|NC_1_cds_NP_1.2_1 [gene=OR4F5] [db_xref=CCDS:CCDS1.2,GeneID:79501]",
                            "lcl|NC_1_cds_XP_2.1_2 [gene=OR4F5] [db_xref=GeneID:79501]"), "header field GeneID:"),
    ("Ensembl", ("ENST00000380152.8 cds chromosome:GRCh38:13:1:2:1 gene:ENSG00000139618.17 gene_symbol:BRCA2",
                 "ENST00000544455.6 cds chromosome:GRCh38:13:1:2:1 gene:ENSG00000139618.17 gene_symbol:BRCA2"),
     "header field gene:"),
    ("FlyBase", ("FBpp0070000 type=CDS; name=Nep3-PA; parent=FBgn0031081,FBtr0070000; species=Dmel;",
                 "FBpp0070001 type=CDS; name=Nep3-PB; parent=FBgn0031081,FBtr0070001; species=Dmel;"),
     "header field parent=FBgn"),
    ("Phytozome", ("Potri.001G000100.1 pacid=1 locus=Potri.001G000100 ID=Potri.001G000100.1.v4.1",
                   "Potri.001G000100.2 pacid=2 locus=Potri.001G000100 ID=Potri.001G000100.2.v4.1"),
     "header field locus="),
    ("NCBI, gene symbols only", ("lcl|NC_1_cds_NP_1.1_1 [gene=LOC1] [protein=x]",
                                 "lcl|NC_1_cds_XP_2.1_2 [gene=LOC1] [protein=x isoform X2]"), "header field [gene="),
    ("GFF-derived (gffread, WormBase)", ("mRNA1 gene=g1 CDS=1-300", "mRNA2 gene=g1 CDS=1-200"), "header field gene="),
    ("TransDecoder on Trinity", ("TRINITY_DN1000_c0_g1_i1.p1 ORF type:complete", "TRINITY_DN1000_c0_g1_i2.p1 x"),
     "ID suffix _iN.pN (TransDecoder on Trinity)"),
    ("Trinity", ("TRINITY_DN1000_c0_g1_i1 len=300", "TRINITY_DN1000_c0_g1_i2 len=300"), "ID suffix _iN (Trinity)"),
    ("AUGUSTUS, BRAKER", ("g1.t1", "g1.t2"), "ID suffix .tN"),
    ("MAKER", ("maker-scaf1-gene-0.5-mRNA-1", "maker-scaf1-gene-0.5-mRNA-2"), "ID suffix -mRNA-N (MAKER)"),
    ("MAKER (-RA, -RB)", ("Gene0001-RA", "Gene0001-RB"), "ID suffix -RA (MAKER)"),
    ("maize (B73 v4)", ("Zm00001d027230_T001", "Zm00001d027230_T002"), "ID suffix _TNNN (maize)"),
    ("Phytozome protein IDs", ("Potri.001G000100.1.p", "Potri.001G000100.2.p"), "ID suffix .N.p (Phytozome)"),
    ("TAIR, Ensembl Plants, Phytozome", ("AT1G01020.1 | Symbols: ARV1", "AT1G01020.2 | Symbols: ARV1"),
     "ID suffix .N"),
    ("RAP-DB (rice)", ("Os01t0100200-01", "Os01t0100200-02"), "ID suffix -NN (RAP-DB)"),
)
FIELD_SHARE = 0.9          # a header field counts when this share of a file's records has it
MAX_MEAN_GROUP = 4.0       # a suffix rule counts when a stem has at most this many records on average
_REFSEQ_RNA = re.compile(r"^(?:[a-z]+\|)?[NX][MR]_\d+(?:\.\d+)?(?:\||$)")
_CONTIG = re.compile(r"^(?:lcl\|)?(.+?)_cds(?:_|$)")      # lcl|NC_000913.3_cds_NP_414542.1_1 -> NC_000913.3
_EVIDENCE = (HEADER_FIELDS[0][1], HEADER_FIELDS[1][1], HEADER_FIELDS[2][1])   # locus tags, GeneIDs


@dataclass
class Detection:
    source: str | None = None          # "gene map", "header field gene:", "ID suffix .N", or None: nothing found
    gene_of: dict = field(default_factory=dict)    # record ID -> gene (only records whose gene is known)
    records: int = 0
    genes: int = 0                     # distinct genes among the records (a record without a gene is its own)
    multi: int = 0                     # genes with more than one record
    extra: int = 0                     # records beyond the first of their gene (what --longest-isoform drops)
    example: tuple | None = None       # (gene, [record IDs]) of one gene with several records
    kind: str | None = None            # "map", "header", "symbol" (a header field that holds gene symbols), "suffix"
    by_suffix: set = field(default_factory=set)   # the records whose gene comes from a suffix rule
    conflicts: int = 0                 # gene symbols shared by records at clearly different loci: not merged
    conflict_example: tuple | None = None   # (symbol, [record IDs], why)
    refseq_ids: int = 0                # records whose ID is a RefSeq transcript accession (when nothing was found)
    in_map: int | None = None          # with a gene map: the records of this set it names
    fallback: str | None = None        # with a gene map: the source of the records it does not name


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


def _split_conflicts(headers, gene_of: dict) -> tuple[int, tuple | None]:
    """Records that share a gene symbol but are clearly different loci (another contig in the record ID, another
    locus tag or GeneID in the header) are taken out of gene_of: (symbols concerned, an example)."""
    groups: dict[str, list[tuple[str, str]]] = {}
    for rid, rest in headers:
        if rid in gene_of:
            groups.setdefault(gene_of[rid], []).append((rid, rest))
    n, example = 0, None
    for sym, members in groups.items():
        if len(members) < 2:
            continue
        why = None
        contigs = {m.group(1) for m in (_CONTIG.match(rid) for rid, _ in members) if m}
        if len(contigs) > 1:
            why = "contigs " + ", ".join(sorted(contigs)[:3])
        else:
            for pat, what in zip(_EVIDENCE, ("locus tags", "locus tags", "GeneIDs")):
                vals = {m.group(1) for m in (pat.search(rest) for _, rest in members) if m}
                if len(vals) > 1:
                    why = f"{what} " + ", ".join(sorted(vals)[:3])
                    break
        if why:
            n += 1
            if example is None:
                example = (sym, [rid for rid, _ in members][:4], why)
            for rid, _ in members:
                del gene_of[rid]
    return n, example


def detect(headers) -> Detection:
    """The gene of each record from the headers [(id, rest)] of one set of files: a header field, else a shared
    ID stem, else nothing (source None)."""
    headers = list(headers)
    ids = [rid for rid, _ in headers]
    n = len(headers)
    if n == 0:
        return Detection()
    for name, pat, stable in HEADER_FIELDS:
        gene_of = {}
        for rid, rest in headers:
            m = pat.search(rest)
            if m and rid not in gene_of:
                gene_of[rid] = ("GeneID:" if name == "GeneID:" else "") + m.group(1).strip()
        if len(gene_of) >= FIELD_SHARE * n:
            det = Detection(source=f"header field {name}", gene_of=gene_of, kind="header" if stable else "symbol")
            if not stable:
                det.conflicts, det.conflict_example = _split_conflicts(headers, gene_of)
            return _summarise(det, ids)
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
        return _summarise(Detection(source=name, gene_of=gene_of, kind="suffix", by_suffix=set(gene_of)), ids)
    det = _summarise(Detection(), ids)
    det.refseq_ids = sum(1 for rid in dict.fromkeys(ids) if _REFSEQ_RNA.match(rid))
    return det


def read_gene_map(path: str, shown: str | None = None) -> dict[str, str]:
    """--gene-map: record ID and gene ID in the first two columns (tab-separated, or blank-separated when a line
    has no tab); `#` lines, empty lines and a header line (id / record / transcript ... gene) are skipped. A record
    named twice with the same gene is fine. ValueError for a line with one column, and for records named with
    two different genes (they are listed)."""
    shown = shown or path
    out: dict[str, str] = {}
    conflicts: list[str] = []
    first = True
    with open(path, encoding="utf-8", errors="backslashreplace") as fh:
        for n, line in enumerate(fh):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t") if "\t" in line else line.split()
            if len(f) < 2 or not f[0].strip() or not f[1].strip():
                raise ValueError(f"{shown}: line {n + 1}: expected two columns (record ID, gene ID)")
            a, b = f[0].strip(), f[1].strip()
            if first and a.lower() in ("id", "record", "record_id", "transcript", "transcript_id", "seq", "seqid",
                                       "protein", "protein_id", "cds") and "gene" in b.lower():
                first = False
                continue
            first = False
            if out.setdefault(a, b) != b:
                conflicts.append(f"{a} (line {n + 1}: {b}; before: {out[a]})")
    if conflicts:
        more = f" and {len(conflicts) - 5:,} more" if len(conflicts) > 5 else ""
        raise ValueError(f"{shown}: {len(conflicts):,} record(s) are given two different genes: "
                         + "; ".join(conflicts[:5]) + more + ". Each record ID must have one gene")
    return out


def from_gene_map(gene_map: dict, headers, detection: bool = True) -> Detection:
    """The genes of one set with a gene map: the map's gene for every record it names; the records it does not
    name get what detect() finds for the set (detection=False: nothing, each is its own gene). headers: [(id,
    rest)], or record IDs alone."""
    headers = [(h, "") if isinstance(h, str) else tuple(h) for h in headers]
    ids = [rid for rid, _ in headers]
    found = detect(headers) if detection else Detection()
    in_map = sum(1 for r in dict.fromkeys(ids) if r in gene_map)
    if in_map == 0:      # the map says nothing about this set
        found = _summarise(found, ids)
        found.in_map = 0
        return found
    gene_of = {r: g for r, g in found.gene_of.items() if r not in gene_map}
    by_suffix = {r for r in found.by_suffix if r in gene_of}
    fallback = found.source if gene_of else None
    gene_of.update((r, gene_map[r]) for r in ids if r in gene_map)
    det = Detection(source="gene map", gene_of=gene_of, kind="map", by_suffix=by_suffix, in_map=in_map,
                    fallback=fallback, conflicts=found.conflicts, conflict_example=found.conflict_example)
    return _summarise(det, ids)


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
