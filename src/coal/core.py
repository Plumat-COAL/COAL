"""Core conversion: protein alignment + coding DNA -> codon alignment.

MVP scope (v0.1):
  * genetic code 1 only
  * gap characters '-' (each becomes '---')
  * DNA must be exactly the coding sequence (3 nt per residue), optionally
    followed by one stop codon not present in the protein
  * every codon must match its residue (IUPAC-ambiguous codons accepted when
    unambiguous at the amino-acid level; 'X' accepts any codon)
Out of scope, rejected with a clear error: UTRs, frameshift digits, other
genetic codes, mismatch tolerance.
"""
from __future__ import annotations

import re

from .codon import FAST_CODE, codon_matches, translations

GAP = "-"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYXBZJUO*")


_RUNS = re.compile(r"-+|[^-]+")
_CODON = re.compile(r"...", re.S)
_INVALID = re.compile("[^" + re.escape(GAP + "".join(sorted(VALID_AA))) + "]")


class CoalError(Exception):
    pass


def pair_sequences(prot_ids: list[str], dna_ids: list[str]) -> list[int]:
    """Index into the DNA list for each protein.

    Same rule as PAL2NAL v14: counts must be equal; if every protein ID is
    present among the DNA IDs, match by ID, otherwise match by position.
    """
    if len(prot_ids) != len(dna_ids):
        raise CoalError(
            f"number of sequences differs (protein: {len(prot_ids)}, DNA: {len(dna_ids)})"
        )
    if not prot_ids:
        raise CoalError("no sequences in input")
    pos = {d: i for i, d in enumerate(dna_ids)}
    if all(p in pos for p in prot_ids):
        return [pos[p] for p in prot_ids]
    return list(range(len(prot_ids)))


def thread_codons(seq_id: str, aln: str, dna: str) -> str:
    """Replace each residue in the aligned protein with its codon, gaps with '---'."""
    aln_u = aln.upper()
    if _INVALID.search(aln_u):
        bad = {c for c in aln_u if c != GAP and c not in VALID_AA}
        digits = sorted(c for c in bad if c.isdigit())
        if digits:
            raise CoalError(f"{seq_id}: frameshift markers {digits} are not supported yet")
        raise CoalError(f"{seq_id}: invalid protein characters {sorted(bad)}")

    n_res = len(aln_u) - aln_u.count(GAP)
    need = 3 * n_res
    if len(dna) == need + 3 and translations(dna[need:]) == {"*"}:
        dna = dna[:need]  # trailing stop codon not in protein
    if len(dna) != need:
        raise CoalError(
            f"{seq_id}: DNA length {len(dna)} does not match {n_res} residues "
            f"(expected {need} or {need + 3} with stop codon; UTRs not supported yet)"
        )

    # Fast path: translate every codon with one dict lookup and compare the whole
    # protein at once. Only if that fails (mismatch, ambiguous base, 'X', ...)
    # fall back to the per-codon check, which also produces the error message.
    residues = aln_u.replace(GAP, "")
    try:
        ok = "".join(map(FAST_CODE.__getitem__, _CODON.findall(dna))) == residues
    except KeyError:  # a codon with a non-ACGTU character
        ok = False
    if not ok:
        _check_codons_slow(seq_id, aln_u, dna)

    # Build output run by run (a run = consecutive residues or consecutive gaps),
    # copying whole DNA slices instead of one codon at a time.
    out = []
    k = 0
    for m in _RUNS.finditer(aln_u):
        n = m.end() - m.start()
        if aln_u[m.start()] == GAP:
            out.append("---" * n)
        else:
            out.append(dna[k:k + 3 * n])
            k += 3 * n
    return "".join(out)


def _check_codons_slow(seq_id: str, aln_u: str, dna: str) -> None:
    k = 0
    for col, aa in enumerate(aln_u, 1):
        if aa == GAP:
            continue
        codon = dna[k:k + 3]
        if not codon_matches(codon, aa):
            raise CoalError(
                f"{seq_id}: alignment column {col}: codon {codon} does not encode {aa}"
            )
        k += 3


def codon_align(
    proteins: list[tuple[str, str]], dnas: list[tuple[str, str]]
) -> list[tuple[str, str]]:
    """Return [(protein_id, codon_aligned_dna)] in protein-alignment order."""
    lengths = {len(s) for _, s in proteins}
    if len(lengths) > 1:
        raise CoalError(f"protein alignment rows have different lengths: {sorted(lengths)}")
    idx = pair_sequences([p[0] for p in proteins], [d[0] for d in dnas])
    return [(pid, thread_codons(pid, aln, dnas[j][1])) for (pid, aln), j in zip(proteins, idx)]
