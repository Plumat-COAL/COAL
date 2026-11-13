"""Write tests/fixtures/isoforms: CDS files with several transcripts per gene, in the header styles of Ensembl,
NCBI and TAIR (synthetic: the sequences are the ortholog pairs g01..g06 / h01..h06 of tests/fixtures/pipeline).

    python tests/tools/make_isoform_fixture.py

Per gene: the full CDS, and for some genes shorter "isoforms" (the CDS without an internal block of codons, or
without its last codons), so that the longest record of a gene is known: transcript 1 of gene 3 is NOT the longest
(the full CDS is its transcript 2), and gene 5 has two transcripts of the same length (the ID that sorts first wins).

  ensembl_query_cds.fasta   >ENSTQ00000000101.1 cds chromosome:T1:1:1000:1 gene:ENSGQ00000000001.4 gene_biotype:...
  ncbi_subject_cds.fasta    >lcl|NC_000001.1_cds_XP_000101.1_1 [gene=HOM1] [protein=...] [protein_id=XP_000101.1]
  tair_query_cds.fasta      >QT1G00010.1 ...     (the gene is the ID without its .N suffix; no gene field)
  gene_map.tsv              record ID -> gene for the TAIR-style and the NCBI-style file
  expected.tsv              file, gene, the record --longest-isoform keeps
  styles/NN_*.fasta         one file per header style of coal.isoforms.STYLES (the table of docs/pipeline.md): the two
                            records of the table's example (the first one the shorter transcript) and one record of
                            another gene; styles/expected.tsv: file, style, source, gene, the record kept
  multi_strain_cds.fasta    two strains in one file, the same gene symbols in both, no locus tags: not merged
  scaffold_numbered_cds.fasta   genes numbered per scaffold (s1.1, s1.2, s1.3 ...): a suffix rule would drop most
  paralog_numbered_cds.fasta    HSP70.1 / HSP70.2 among other genes: the documented failure of the .N rule
  refseq_rna.fasta          RefSeq transcript accessions without a gene field
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
FIX = ROOT / "tests" / "fixtures"


def read(path):
    recs, cur = {}, None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            cur = None if name in recs else recs.setdefault(name, [])
        elif cur is not None:
            cur.append(line.strip())
    return {k: "".join(v) for k, v in recs.items()}


def variants(cds: str, k: int) -> list[str]:
    """The transcripts of gene k (1-based), in the order they are written."""
    body, stop = cds[:-3], cds[-3:]
    cut_mid = body[:30] + body[60:] + stop          # an exon of 10 codons skipped
    cut_end = body[:len(body) - 45] + stop          # 15 codons shorter at the 3' end
    alt_end = body[:len(body) - 3] + "GCA" + stop   # same length, another last codon
    return {1: [cds, cut_mid, cut_end], 2: [cds], 3: [cut_mid, cds], 4: [cds, cut_end], 5: [cds, alt_end],
            6: [cds]}[k]


def wrap(seq: str) -> str:
    return "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))


def main():
    out = FIX / "isoforms"
    out.mkdir(exist_ok=True)
    q = read(FIX / "pipeline" / "query_cds.fasta")
    s = read(FIX / "pipeline" / "subject_cds.fasta")
    ens, ncbi, tair, gmap, expected = [], [], [], [], []
    for k in range(1, 7):
        qv, sv = variants(q[f"g{k:02d}"], k), variants(s[f"h{k:02d}"], k)

        def longest(ids, seqs):
            return min(zip(ids, seqs), key=lambda t: (-len(t[1]), t[0].encode()))[0]
        gene = f"ENSGQ{k:011d}"
        ids = [f"ENSTQ{k * 100 + t:011d}.{t}" for t in range(1, len(qv) + 1)]
        for rid, seq in zip(ids, qv):
            ens.append(f">{rid} cds chromosome:T1:{k}:{k * 1000}:{k * 1000 + len(seq)}:1 gene:{gene}.4 "
                       f"gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:QRY{k} "
                       f"description:synthetic test gene {k}\n{wrap(seq)}\n")
        expected.append(("ensembl_query_cds.fasta", f"{gene}.4", longest(ids, qv)))
        ids = [f"lcl|NC_000001.1_cds_XP_{k * 100 + t:06d}.1_{len(ncbi) + t}" for t in range(1, len(sv) + 1)]
        for t, (rid, seq) in enumerate(zip(ids, sv), 1):
            ncbi.append(f">{rid} [gene=HOM{k}] [db_xref=GeneID:{9000 + k}] [protein=synthetic protein {k} isoform "
                        f"X{t}] [protein_id=XP_{k * 100 + t:06d}.1] [location=join({k * 1000}..{k * 1000 + len(seq)})] "
                        f"[gbkey=CDS]\n{wrap(seq)}\n")
            gmap.append((rid, f"HOM{k}"))
        expected.append(("ncbi_subject_cds.fasta", f"HOM{k}", longest(ids, sv)))
        ids = [f"QT1G{k * 10:05d}.{t}" for t in range(1, len(qv) + 1)]
        for rid, seq in zip(ids, qv):
            tair.append(f">{rid} | Symbols: QRY{k} | synthetic test gene {k} | chr1:{k * 1000}-{k * 1000 + len(seq)} "
                        f"FORWARD LENGTH={len(seq)}\n{wrap(seq)}\n")
            gmap.append((rid, f"QT1G{k * 10:05d}"))
        expected.append(("tair_query_cds.fasta", f"QT1G{k * 10:05d}", longest(ids, qv)))
    for name, text in (("ensembl_query_cds.fasta", "".join(ens)), ("ncbi_subject_cds.fasta", "".join(ncbi)),
                       ("tair_query_cds.fasta", "".join(tair)),
                       ("gene_map.tsv", "record\tgene\n" + "".join(f"{a}\t{b}\n" for a, b in gmap)),
                       ("expected.tsv", "file\tgene\tkept\n" + "".join("\t".join(e) + "\n" for e in expected))):
        (out / name).write_text(text, encoding="utf-8", newline="\n")
    print(f"{out}: {len(ens)} + {len(ncbi)} + {len(tair)} records of 6 genes each")
    styles(out, q)


# the record of another gene in each style file, in the order of coal.isoforms.STYLES
OTHER = (
    "lcl|NC_1_cds_NP_3.1_3 [gene=thrB] [locus_tag=b0003] [protein=y]",
    "lcl|NC_1_cds_NP_3.1_3 [gene=OR4F29] [db_xref=GeneID:729759]",
    "ENST00000269305.9 cds chromosome:GRCh38:17:1:2:1 gene:ENSG00000141510.18 gene_symbol:TP53",
    "FBpp0070002 type=CDS; name=CG9570-PA; parent=FBgn0031085,FBtr0070002; species=Dmel;",
    "Potri.001G000200.1 pacid=3 locus=Potri.001G000200 ID=Potri.001G000200.1.v4.1",
    "lcl|NC_1_cds_NP_3.1_3 [gene=LOC2] [protein=y]",
    "mRNA3 gene=g2 CDS=1-9",
    "TRINITY_DN1001_c0_g1_i1.p1 y",
    "TRINITY_DN1001_c0_g1_i1 len=3",
    "g2.t1",
    "maker-scaf1-gene-0.9-mRNA-1",
    "Gene0002-RA",
    "Zm00001d027231_T001",
    "Potri.001G000200.1.p",
    "AT1G01030.1 | x",
    "Os01t0100100-01",
)


def styles(out, q):
    sys.path.insert(0, str(ROOT / "src"))
    from coal import isoforms as I
    assert len(OTHER) == len(I.STYLES)
    d = out / "styles"
    d.mkdir(exist_ok=True)
    for old in d.glob("*.fasta"):
        old.unlink()
    full, other = q["g01"], q["g02"]
    short = full[:len(full) - 48] + full[-3:]      # 15 codons shorter at the 3' end: the first record is NOT kept
    rows = []
    for n, ((style, (h1, h2), source), h3) in enumerate(zip(I.STYLES, OTHER), 1):
        name = f"{n:02d}_" + re.sub(r"[^a-z0-9]+", "_", style.lower()).strip("_") + ".fasta"
        (d / name).write_text(f">{h1}\n{wrap(short)}\n>{h2}\n{wrap(full)}\n>{h3}\n{wrap(other)}\n", encoding="utf-8",
                              newline="\n")
        det = I.detect(I.read_headers([str(d / name)]))
        rows.append((name, style, source, det.gene_of[h1.split()[0]], h2.split()[0]))
    (d / "expected.tsv").write_text("file\tstyle\tsource\tgene\tkept\n" + "".join("\t".join(r) + "\n" for r in rows),
                                    encoding="utf-8", newline="\n")
    genes = [q[f"g{k:02d}"] for k in range(1, 7)]
    files = {
        "multi_strain_cds.fasta": [(f"lcl|STRAIN{st}_cds_WP_{st}{k:03d}.1_{k} [gene={g}] [protein=x]", genes[k - 1])
                                   for st in (1, 2) for k, g in enumerate(("dnaA", "dnaN", "recF"), 1)],
        "scaffold_numbered_cds.fasta": [(f"s{sc}.{k}", genes[(sc - 1) * 3 + k - 1]) for sc in (1, 2) for k in (1, 2, 3)],
        "paralog_numbered_cds.fasta": list(zip(("HSP70.1", "HSP70.2", "ACT1", "ACT2", "TUB1", "RPL3"), genes)),
        "refseq_rna.fasta": [("NM_000001.1 Homo sapiens synthetic gene 1 (SYN1), mRNA", genes[0]),
                             ("XM_000002.3 PREDICTED: Homo sapiens synthetic gene 1 (SYN1), transcript variant X1, "
                              "mRNA", genes[0][:len(genes[0]) - 48] + genes[0][-3:]),
                             ("NM_000003.2 Homo sapiens synthetic gene 2 (SYN2), mRNA", genes[1])],
    }
    for name, recs in files.items():
        (out / name).write_text("".join(f">{h}\n{wrap(seq)}\n" for h, seq in recs), encoding="utf-8", newline="\n")
    print(f"{d}: {len(rows)} header styles; {', '.join(files)}")


if __name__ == "__main__":
    main()
