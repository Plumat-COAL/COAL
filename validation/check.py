#!/usr/bin/env python3
"""Score a codon-alignment output against its inputs (and the truth, when there is one).

Reference-free checks (work on any data, at any scale):
  shape      every row the same length; the output has exactly the protein-alignment IDs, in order
  translate  each protein column maps to a DNA token of the expected width: a residue -> a codon
             encoding it (first column: any start codon of the table may stand for M; 'X' -> any
             codon; mismatches are allowed only where COAL/PAL2NAL warned), '-' -> all gaps,
             digit d -> d bases then gap padding
  cds        removing the gaps from each row gives a contiguous piece of the input CDS
             (UTRs/stop codons may be trimmed off the ends, nothing else)
Truth comparison (simulated data and Ensembl Compara): rows exactly equal, and per-column agreement.

Use from run.py, or on one output:
    python validation/check.py CASE_DIR OUTPUT.fasta [--stderr ERR.txt]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from simulate import code, translate  # noqa: E402

_NONLETTER = re.compile(r"[^A-Za-z]")
_WARN = re.compile(r"WARNING:?\s+(\S+)\s+pepAlnPos\s+(\d+)", re.I)


def read_fasta(path: Path, letters_only=False) -> list[tuple[str, str]]:
    recs: list[tuple[str, list[str]]] = []
    with open(path, encoding="ascii", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                recs.append((line[1:].split()[0] if line[1:].split() else "", []))
            elif line.strip() and recs:
                if letters_only and line.startswith("#"):
                    continue
                recs[-1][1].append(_NONLETTER.sub("", line) if letters_only else "".join(line.split()))
    return [(i, "".join(s)) for i, s in recs]


def pair_cds(pep_ids: list[str], nuc: list[tuple[str, str]]) -> list[str]:
    d = dict(nuc)
    if all(i in d for i in pep_ids):
        return [d[i] for i in pep_ids]
    if len(nuc) != len(pep_ids):
        raise ValueError("cannot pair: IDs differ and counts differ")
    return [s for _, s in nuc]


def column_widths(pep_rows: list[str]) -> list[int]:
    w = [3] * (len(pep_rows[0]) if pep_rows else 0)
    for r in pep_rows:
        for j, ch in enumerate(r):
            if ch.isdigit():
                w[j] = max(w[j], -(-int(ch) // 3) * 3)
    return w


def _starts(table: int) -> set[str]:
    return set(code(table).starts) | {"ATG"}


def check_case(case_dir: Path, out_path: Path, stderr_text: str = "", meta: dict | None = None,
               max_examples: int = 5) -> dict:
    meta = meta or json.loads((case_dir / "meta.json").read_text())
    opts = meta.get("opts", [])
    table = int(opts[opts.index("-codontable") + 1]) if "-codontable" in opts else 1
    pep = read_fasta(case_dir / meta["pep"])
    nuc = [r for f in meta["nuc"] for r in read_fasta(case_dir / f, letters_only=True)]
    ids = [i for i, _ in pep]
    pep_rows = [s.upper().replace(".", "-") for _, s in pep]
    cds = pair_cds(ids, nuc)
    res = {"rows_total": len(ids), "problems": [], "examples": []}

    def problem(kind, detail):
        res["problems"].append(kind)
        if len(res["examples"]) < max_examples:
            res["examples"].append(f"{kind}: {detail}")

    out = read_fasta(out_path) if out_path.exists() else []
    res["output_rows"] = len(out)
    warned = {(m.group(1), int(m.group(2))) for m in _WARN.finditer(stderr_text)}
    res["warnings_found"] = len(warned)
    res["warnings_expected"] = int(meta.get("expected_warnings", 0))
    if not out:
        problem("no_output", "empty output")
        return _finish(res)
    if [i for i, _ in out] != ids:
        problem("ids", "output IDs/order differ from the protein alignment")
    out_rows = [s for _, s in out]
    widths = column_widths(pep_rows)
    total_w = sum(widths)
    starts = _starts(table)
    ok_rows = 0
    for k, (sid, prow, orow) in enumerate(zip(ids, pep_rows, out_rows)):
        row_ok = True
        if len(orow) != total_w:
            problem("shape", f"{sid}: length {len(orow)} != {total_w}")
            continue
        pos = 0
        first_res = True
        for j, (ch, w) in enumerate(zip(prow, widths)):
            tok = orow[pos:pos + w]
            pos += w
            if ch == "-":
                if tok != "-" * w:
                    row_ok = False
                    problem("gap", f"{sid} col {j + 1}: expected gaps, got {tok}")
                continue
            if ch.isdigit():
                d = int(ch)
                if not (tok[:d].isalpha() and set(tok[d:]) <= {"-"}):
                    row_ok = False
                    problem("frameshift", f"{sid} col {j + 1}: digit {d}, got {tok}")
                continue
            cod, pad = tok[:3], tok[3:]
            if not cod.isalpha() or set(pad) - {"-"}:
                row_ok = False
                problem("translate", f"{sid} col {j + 1}: {ch} -> {tok}")
                continue
            aa = translate(cod, table)
            fine = (aa == ch or ch == "X" or
                    (first_res and ch == "M" and cod.upper().replace("U", "T") in starts) or
                    (sid, j + 1) in warned)
            if not fine:
                row_ok = False
                problem("translate", f"{sid} col {j + 1}: {ch} but codon {cod} = {aa}")
            first_res = False
        degapped = orow.replace("-", "")
        if degapped and degapped not in cds[k]:
            row_ok = False
            problem("cds", f"{sid}: degapped row is not a substring of its CDS")
        ok_rows += row_ok
    res["rows_refree_ok"] = ok_rows

    truth_path = case_dir / meta["truth"] if meta.get("truth") else None
    if truth_path and truth_path.exists():
        truth = dict(read_fasta(truth_path))
        exact = cols_ok = cols_total = 0
        for sid, orow in zip(ids, out_rows):
            t = truth.get(sid, "")
            exact += orow == t
            if orow != t:
                problem("truth", f"{sid}: differs from truth")
            pos = 0
            for w in widths:
                cols_total += 1
                cols_ok += orow[pos:pos + w] == t[pos:pos + w]
                pos += w
        res["rows_truth_exact"] = exact
        res["col_agreement"] = cols_ok / cols_total if cols_total else 1.0
    wfile = case_dir / "expected_warnings.tsv"
    if wfile.exists():
        exp = {(l.split("\t")[0], int(l.split("\t")[1])) for l in wfile.read_text().splitlines()[1:] if l}
        missing, extra = exp - warned, warned - exp
        if missing or extra:
            problem("warnings", f"{len(missing)} expected warnings missing, {len(extra)} unexpected")
    elif warned:
        problem("warnings", f"{len(warned)} unexpected mismatch warnings")
    return _finish(res)


def _finish(res):
    kinds = res.pop("problems")
    res["problem_counts"] = {k: kinds.count(k) for k in sorted(set(kinds))}
    res["pass"] = not kinds
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("case_dir", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--stderr", type=Path)
    a = ap.parse_args()
    err = a.stderr.read_text(errors="replace") if a.stderr else ""
    print(json.dumps(check_case(a.case_dir, a.output, err), indent=1))


if __name__ == "__main__":
    main()
