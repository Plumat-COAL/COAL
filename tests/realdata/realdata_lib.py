"""Helpers for the real-data round-trip tests (tests/test_realdata.py) and the local survey scripts.

A real codon alignment is split into what a user would give COAL:
  * the protein alignment: each codon column translated with the dataset's genetic code (NCBI tables, written
    out here independently of COAL's own tables), a gap codon '---' -> '-';
  * the unaligned coding DNA: the same rows without gaps.
COAL's FASTA output must then reproduce the codon alignment (upper case).

Parsers: FASTA, NEXUS (sequential or interleaved MATRIX, quoted names), PAML/PHYLIP sequential (as in PAML's
examples: name on its own line or followed by 2+ spaces; spaces and position numbers inside the sequence ignored).
Standard library only.
"""
from __future__ import annotations

import gzip
import itertools
import re
from pathlib import Path

# NCBI translation tables, codons in TCAG order (base1 slowest). From the NCBI "The Genetic Codes" page.
_B = "TCAG"
CODONS = ["".join(c) for c in itertools.product(_B, repeat=3)]
NCBI_AAS = {
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    2: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    3: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    5: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    6: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    9: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    10: "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT",
         "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def table(code: int) -> dict[str, str]:
    return dict(zip(CODONS, NCBI_AAS[code]))


def translate_codon(codon: str, code: int) -> str:
    """'---' -> '-'; IUPAC codon -> its amino acid if every reading agrees, else 'X'. Partial gaps -> ValueError."""
    c = codon.upper()
    if c == "---":
        return "-"
    if "-" in c:
        raise ValueError(f"partly gapped codon {codon!r}")
    t = table(code)
    try:
        aas = {t[a + b + d] for a in IUPAC[c[0]] for b in IUPAC[c[1]] for d in IUPAC[c[2]]}
    except KeyError:
        raise ValueError(f"non-IUPAC codon {codon!r}")
    return aas.pop() if len(aas) == 1 else "X"


def _open_text(path: Path) -> str:
    raw = gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes()
    return raw.decode("utf-8")


def parse_fasta(text: str) -> list[tuple[str, str]]:
    recs, name, buf = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        elif name is not None and not line.startswith("("):  # HyPhy appends a Newick tree after the last record
            buf.append(re.sub(r"\s", "", line))
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def parse_nexus(text: str) -> list[tuple[str, str]]:
    text = re.sub(r"\[[^\]]*\]", "", text)  # comments
    m = re.search(r"\bmatrix\b(.*?);", text, re.I | re.S)
    if not m:
        raise ValueError("no MATRIX")
    order, seqs = [], {}
    for line in m.group(1).splitlines():
        line = line.strip()
        if not line:
            continue
        mm = re.match(r"'([^']*)'\s*(.*)$", line) or re.match(r"(\S+)\s*(.*)$", line)
        name, rest = mm.group(1), re.sub(r"\s", "", mm.group(2))
        if name not in seqs:
            order.append(name)
            seqs[name] = []
        seqs[name].append(rest)
    return [(n, "".join(seqs[n])) for n in order]


def parse_paml(text: str) -> list[tuple[str, str]]:
    lines = text.splitlines()
    i = 0
    while not lines[i].strip():
        i += 1
    n, length = (int(x) for x in lines[i].split()[:2])
    i += 1
    recs = []
    for _ in range(n):
        while not lines[i].strip():
            i += 1
        m = re.match(r"^(\S.*?)(?:\s{2,}|\t)(.*)$", lines[i].rstrip())
        name, seq = (m.group(1), m.group(2)) if m else (lines[i].strip(), "")
        seq = re.sub(r"[^A-Za-z\-?.*]", "", seq)
        i += 1
        while len(seq) < length:
            seq += re.sub(r"[^A-Za-z\-?.*]", "", lines[i])
            i += 1
        if len(seq) != length:
            raise ValueError(f"{name}: {len(seq)} != {length}")
        recs.append((name.strip(), seq))
    return recs


def load(path: Path, fmt: str) -> list[tuple[str, str]]:
    text = _open_text(path)
    return {"fasta": parse_fasta, "nexus": parse_nexus, "paml": parse_paml}[fmt](text)


def clean_rows(rows: list[tuple[str, str]], missing_to_n: bool = True) -> list[tuple[str, str]]:
    """Upper case; '?' (NEXUS/PAML missing data) -> N when asked. Nothing else is changed."""
    out = []
    for name, seq in rows:
        s = seq.upper()
        if missing_to_n:
            s = s.replace("?", "N")
        out.append((name, s))
    return out


def survey(rows: list[tuple[str, str]], code: int) -> dict:
    """Facts about a codon alignment: what blocks an exact round trip."""
    lengths = {len(s) for _, s in rows}
    info = {"n": len(rows), "lengths": sorted(lengths), "partial_gap_codons": 0, "bad_codons": 0,
            "internal_stops": 0, "terminal_stops": 0, "ambiguous_codons": 0, "x_codons": 0, "rows_with_partial": []}
    for name, s in rows:
        aas = []
        partial = 0
        for k in range(0, len(s) - len(s) % 3, 3):
            c = s[k:k + 3]
            try:
                aa = translate_codon(c, code)
            except ValueError:
                if "-" in c:
                    partial += 1
                else:
                    info["bad_codons"] += 1
                aa = "?"
            aas.append(aa)
            if aa != "-" and any(b not in "ACGT" for b in c):
                info["ambiguous_codons"] += 1
                if aa == "X":
                    info["x_codons"] += 1
        if partial:
            info["partial_gap_codons"] += partial
            info["rows_with_partial"].append(name)
        prot = "".join(aas).rstrip("-")
        info["terminal_stops"] += prot.endswith("*")
        info["internal_stops"] += prot[:-1].count("*")
    return info


def derive(rows: list[tuple[str, str]], code: int, partial: str = "error") -> tuple[str, str, list]:
    """(protein alignment FASTA, unaligned DNA FASTA, expected codon rows).

    partial="gap": a partly gapped codon (fragment end, e.g. '--G') is a gap in the protein and its bases stay in
    the DNA, as when each fragment is translated in its own frame; the expected row has '---' there and COAL must
    skip the stray bases like a UTR. partial="error": ValueError."""
    pep, nuc, exp = [], [], []
    for name, s in rows:
        if len(s) % 3:
            raise ValueError(f"{name}: length {len(s)} not a multiple of 3")
        p, e = [], []
        for k in range(0, len(s), 3):
            c = s[k:k + 3]
            if partial == "gap" and "-" in c and c != "---":
                p.append("-"); e.append("---")
            else:
                p.append(translate_codon(c, code)); e.append(c)
        pep.append(f">{name}\n{''.join(p)}\n")
        nuc.append(f">{name}\n{s.replace('-', '')}\n")
        exp.append((name.split()[0], "".join(e)))
    return "".join(pep), "".join(nuc), exp


def fasta_rows(text: str) -> list[tuple[str, str]]:
    return [(n.split()[0] if n else n, s.upper()) for n, s in parse_fasta(text)]


def read_clustal(text: str) -> list[tuple[str, str]]:
    order, seqs = [], {}
    for line in text.splitlines()[1:]:
        if not line.strip() or line[0] in " \t":
            continue
        n, s = line.split()[:2]
        if n not in seqs:
            order.append(n)
            seqs[n] = ""
        seqs[n] += s
    return [(n, seqs[n]) for n in order]


def parse_output(fmt: str, text: str) -> list[tuple[str, str]]:
    """(id, codon row) from COAL/PAL2NAL output in any of the four formats, upper case."""
    if fmt == "fasta":
        return fasta_rows(text)
    lines = text.splitlines()
    if fmt == "paml":
        n, length = (int(x) for x in lines[0].split()[:2])
        rows, i = [], 1
        for _ in range(n):
            name = lines[i].strip()
            i += 1
            s = ""
            while len(s) < length:
                s += lines[i].strip()
                i += 1
            rows.append((name, s.upper()))
        return rows
    order, seqs = [], {}
    for line in (lines[1:] if fmt == "clustal" else lines):
        if not line.strip() or line[0] in " \t":
            continue  # blank, amino-acid line (codon format), '#' line
        parts = line.split(None, 1)
        name, rest = parts[0], (parts[1] if len(parts) > 1 else "")
        if name not in seqs:
            order.append(name)
            seqs[name] = []
        seqs[name].append(rest.replace(" ", ""))
    return [(n, "".join(seqs[n]).upper()) for n in order]


def nogap_rule(rows: list[tuple[str, str]], code: int) -> list[tuple[str, str]]:
    """COAL's documented -nogap: drop codon columns with a gap or a stop codon of the selected code (FZ-5)."""
    stops = {c for c, a in table(code).items() if a == "*"}
    ncol = len(rows[0][1]) // 3
    keep = [k for k in range(ncol)
            if all("-" not in s[3 * k:3 * k + 3] and s[3 * k:3 * k + 3] not in stops for _, s in rows)]
    return [(n, "".join(s[3 * k:3 * k + 3] for k in keep)) for n, s in rows]


def back_translation_mismatches(pep_rows: list[tuple[str, str]], out_rows: list[tuple[str, str]], code: int) -> int:
    """Codons in COAL's output that are not the protein alignment's residue (gap <-> '---'; X/B/Z/J/O/U accept any
    codon; an ambiguous codon that translates to X is not counted; the first Met is not checked). No frameshift
    digits expected."""
    bad = 0
    for (_, p), (_, o) in zip(pep_rows, out_rows):
        if len(o) != 3 * len(p):
            return -1
        first = next((i for i, a in enumerate(p) if a != "-"), -1)
        for i, aa in enumerate(p.upper()):
            if i == first and aa == "M":
                continue  # first Met: any start codon of the code (not checked here)
            c = o[3 * i:3 * i + 3]
            if aa == "-":
                bad += c != "---"
                continue
            if aa in "XBZJOU":
                continue
            try:
                t = translate_codon(c, code)
            except ValueError:
                bad += 1
                continue
            bad += t not in (aa, "X")
    return bad
