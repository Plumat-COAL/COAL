"""Alignment viewer: self-contained HTML showing exactly where each gene breaks.

    python -m coal_app.viewer pep.aln codon.fasta [--report r.json] [--dna nuc.fasta]
                              [--codontable N] -o view.html

No Gradio needed. One row per sequence; every codon sits under its amino acid,
tinted by amino-acid property (colour-blind-safe palette; the letter is always
shown). Gaps are light '---'. Mismatched codons get an outline and a '≠' symbol,
frameshift columns an outline and a '↯' symbol. Problems in the health list link
to the exact cell (plain anchor links; a little optional JS also highlights the
whole column).

If COAL wrote no row for a gene (it failed), and the DNA is given, the viewer
places that gene's codons itself and finds where they stop matching: if shifting
the reading frame by 1 or 2 bases makes them match again, that column is marked
as a frameshift.
"""
from __future__ import annotations

import argparse
import html
import json
import os
import re
import sys

from . import gcode

MAX_SEQS = 200
MAX_COLS = 2000
MAX_CELLS = 200_000

# amino-acid property groups -> CSS class (tints of the Okabe-Ito palette)
GROUPS = {
    **dict.fromkeys("AVLIM", "hy"),   # hydrophobic
    **dict.fromkeys("FWY", "ar"),     # aromatic
    **dict.fromkeys("KRH", "po"),     # positive
    **dict.fromkeys("DE", "ne"),      # negative
    **dict.fromkeys("STNQ", "pl"),    # polar
    **dict.fromkeys("GPC", "sp"),     # special
    "*": "st",                        # stop
}
LEGEND = [("hy", "Hydrophobic (A V L I M)"), ("ar", "Aromatic (F W Y)"), ("po", "Positive (K R H)"),
          ("ne", "Negative (D E)"), ("pl", "Polar (S T N Q)"), ("sp", "Special (G P C)"), ("st", "Stop (*)")]


# ---------------------------------------------------------------- reading

def read_alignment(path: str, dna: bool = False, want=None) -> list[tuple[str, str]]:
    """(id, sequence) records from FASTA or CLUSTAL, in file order (dna=True: DNA FASTA rules).
    CHANGED(R3-A1): was plain UTF-8, first line had to start with '>' or 'CLUSTAL' (a BOM or '#' preamble the
    engine accepts was "not FASTA or CLUSTAL"); now the engine's rules, via seqformat.
    See overnight/integration/CHANGES.md
    CHANGED(GUI-10): want: only these IDs (a set), see seqformat.iter_records."""
    from .seqformat import read_records
    return read_records(path, dna=dna, want=want)


def slug(seq_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]", "_", seq_id)


def cell_id(seq_id: str, col: int) -> str:
    return f"coalv-{slug(seq_id)}-c{col}"


def row_id(seq_id: str) -> str:
    return f"coalv-{slug(seq_id)}"


# ---------------------------------------------------------------- preview: which rows and columns are drawn

def select_rows(ids: list[str], report: dict | None, max_seqs: int, prefer: set | None = None) -> tuple[list[str], dict]:
    """CHANGED(GUI-10): new; was always the first max_seqs sequences. Which sequences the preview draws, in the
    order drawn: all of them in file order when they fit; otherwise every failed and then every warning sequence
    (file order, up to the cap, keeping some room) and then OK ones spread evenly through the file for context.
    prefer: problem sequences to put first of all (those with a problem in the columns the preview is centred on).
    Returns (ids, {"problems": n with problems, "problems_shown": k, "context": m, "problems_first": bool})."""
    if len(ids) <= max_seqs:
        return list(ids), {"problems": None, "problems_shown": None, "context": 0, "problems_first": False}
    st = {s.get("id"): s.get("status") for s in (report or {}).get("sequences") or []}
    failed = [i for i in ids if st.get(i) == "failed"]
    warning = [i for i in ids if st.get(i) == "warning"]
    problems = failed + warning
    if prefer:
        problems = [i for i in problems if i in prefer] + [i for i in problems if i not in prefer]
    if not problems:
        others = list(ids)
        picked = [others[k] for k in _even_idx(len(others), max_seqs)] if st else others[:max_seqs]
        return picked, {"problems": 0, "problems_shown": 0, "context": len(picked), "problems_first": False}
    room = min(max(10, max_seqs // 10), max_seqs // 2)  # rows kept for OK sequences, for context
    probset = set(problems)
    others = [i for i in ids if i not in probset]
    k = min(len(problems), max_seqs - min(room, len(others)))
    shown = problems[:k]
    ctx = [others[j] for j in _even_idx(len(others), max_seqs - k)]
    return shown + ctx, {"problems": len(problems), "problems_shown": k, "context": len(ctx), "problems_first": True}


def _even_idx(n: int, k: int) -> list[int]:
    if k >= n:
        return list(range(n))
    if k <= 0:
        return []
    if k == 1:
        return [0]
    return sorted({round(i * (n - 1) / (k - 1)) for i in range(k)})


def preview_columns(nrows: int, ncols: int, max_cols: int = MAX_COLS, max_cells: int = MAX_CELLS) -> int:
    """How many columns the preview draws for nrows sequences (at least 50, at most max_cols, about max_cells
    cells)."""
    return max(0, min(ncols, max_cols, max(50, max_cells // max(1, nrows))))


def column_window(ncols: int, cols: int, centre: int | None = None) -> tuple[int, int]:
    """(first, last) column drawn, 1-based: from column 1, or centred on `centre`, always inside 1..ncols."""
    if ncols <= 0 or cols <= 0:
        return 1, 0
    if centre is None or ncols <= cols:
        return 1, min(ncols, cols)
    first = max(1, min(int(centre) - cols // 2, ncols - cols + 1))
    return first, first + cols - 1


# ---------------------------------------------------------------- placing codons for failed genes

def _ok(aa: str, codon: str, code: dict, starts: frozenset, first: bool) -> bool:
    """Does `codon` code for `aa`? IUPAC ambiguity codes count when every expansion does (the engine's rule)."""
    if len(codon) < 3:
        return False
    # CHANGED(APP-5): was exact codons only and only X accepted; now IUPAC codes (all expansions must fit) and
    # B Z J O U like X, as the engine does. See overnight/app/CHANGES.md
    if aa in "XBZJOU":  # the engine treats B Z J O U as X: any codon is accepted
        return True
    exp = gcode.expand(codon)
    if not exp:
        return False
    return all(code.get(c) == aa or (first and aa == "M" and c in starts) for c in exp)


def diagnose(protein_row: str, dna: str, codontable: int = 1) -> dict:
    """Place a gene's codons under its protein row and find where they stop matching.

    Returns {"cells": {col: (codon, flag)}, "events": [...], "matched": n} where flag is
    "", "mismatch", "frameshift" or "missing", and events are
    {"code": "frameshift", "column": c, "shift": s, "text": ...} or {"code": "inconsistent", "column": c, ...}.
    """
    code, starts = gcode.table(codontable)
    dna = re.sub(r"[^A-Za-z]", "", dna).upper().replace("U", "T")
    res = [(i + 1, a.upper()) for i, a in enumerate(protein_row) if a.isalpha()]
    out: dict = {"cells": {}, "events": [], "matched": 0}
    if not res or len(dna) < 3:
        return out
    # CHANGED(APP-2): was "first 3 residues must match exactly"; now up to 1 in 4 of the first 12 may mismatch. See overnight/app/CHANGES.md
    # where does the gene start? the leftmost offset where most of the first 12 residues match
    # (mismatches allowed: an alternative start codon or one odd early codon must not hide the gene)
    want = min(len(res), 12)
    allowed = max(0, want // 4)  # up to 3 of 12 may mismatch
    best = (0, -1)
    for off in range(0, len(dna) - 2):
        score = 0
        for k in range(want):
            if _ok(res[k][1], dna[off + 3 * k: off + 3 * k + 3], code, starts, k == 0):
                score += 1
            elif k + 1 - score > allowed:
                break
        if score > best[1]:
            best = (off, score)
            if score == want:
                break
    if best[1] < max(min(3, len(res)), want - allowed):
        out["events"].append({"code": "inconsistent", "column": res[0][0],
                              "text": "The DNA does not spell the start of this protein in any reading frame."})
        return out
    pos = best[0]
    reported_break = False
    k = 0
    while k < len(res):
        col, aa = res[k]
        codon = dna[pos:pos + 3]
        if len(codon) < 3:
            for col2, _ in res[k:]:
                out["cells"][col2] = (codon.ljust(3, " ") if col2 == col else "   ", "missing")
            if not reported_break:
                out["events"].append({"code": "missing_dna", "column": col,
                                      "text": f"The DNA ends before the protein does (at column {col})."})
            break
        if _ok(aa, codon, code, starts, k == 0):
            out["cells"][col] = (codon, "")
            out["matched"] += 1
            pos += 3
            k += 1
            continue
        # a mismatch: isolated, or the start of a frameshift?
        ahead = sum(_ok(a, dna[pos + 3 * j: pos + 3 * j + 3], code, starts, False)
                    for j, (_, a) in enumerate(res[k + 1:k + 6], 1))
        window = min(5, len(res) - k - 1)
        if window and ahead >= window - 1:
            out["cells"][col] = (codon, "mismatch")
            pos += 3
            k += 1
            continue
        shifted = None
        for s in (-1, 1, -2, 2):
            n = min(10, len(res) - k - 1)
            if n <= 0:
                break
            hits = sum(_ok(a, dna[pos + 3 + s + 3 * j: pos + 3 + s + 3 * j + 3], code, starts, False)
                       for j, (_, a) in enumerate(res[k + 1:k + 1 + n]))
            if hits >= max(1, int(0.8 * n)):
                shifted = s
                break
        if shifted is not None:
            what = {-1: "one base is missing", -2: "two bases are missing",
                    1: "there is one extra base", 2: "there are two extra bases"}[shifted]
            out["cells"][col] = (codon, "frameshift")
            out["events"].append({"code": "frameshift", "column": col, "shift": shifted,
                                  "text": f"Frameshift at column {col}: {what} in the DNA, so from here on the "
                                          f"codons are read out of step."})
            pos += 3 + shifted
            k += 1
            continue
        # CHANGED(R2-A12): was "the DNA stops matching for good" here (a cluster of TGA/ATA mismatches under the
        # wrong code made engine-aligned genes "Failed"); now no shift signature = an in-frame mismatch, and a
        # break is only called below when most of the rest of the gene mismatches. See overnight/integration/CHANGES.md
        out["cells"][col] = (codon, "mismatch")
        pos += 3
        k += 1
    if not out["events"]:
        _call_break(res, out)
    return out


BREAK_MIN_CODONS = 10  # a break needs at least this many codons after it, more than half of them mismatching


def _call_break(res, out) -> None:
    """Add an 'inconsistent' event at the first column from which more than half of the remaining (at least
    BREAK_MIN_CODONS) codons do not match, if there is one."""
    flags = [out["cells"].get(col, ("", "missing"))[1] for col, _ in res]
    bad_after = 0
    first = None
    for i in range(len(flags) - 1, -1, -1):
        bad_after += flags[i] in ("mismatch", "missing")
        n = len(flags) - i
        if flags[i] == "mismatch" and n >= BREAK_MIN_CODONS and 2 * bad_after > n:
            first = i
    if first is not None:
        col = res[first][0]
        out["events"].append({"code": "inconsistent", "column": col,
                              "text": f"From column {col} on, the DNA no longer spells this protein."})


# ---------------------------------------------------------------- model

def _report_marks(report: dict | None, only: set | None = None) -> dict[str, dict]:
    marks: dict[str, dict] = {}
    for s in (report or {}).get("sequences", []):
        if only is not None and s.get("id") not in only:  # CHANGED(GUI-10): only the sequences drawn
            continue
        m = {"status": s.get("status"), "mismatch": set(), "frameshift": set()}
        for msg in s.get("messages") or []:
            c = msg.get("column")
            if isinstance(c, int):
                if msg.get("code") == "frameshift":
                    m["frameshift"].add(c)
                elif msg.get("code") == "mismatch":
                    m["mismatch"].add(c)
        marks[s.get("id")] = m
    return marks


_DIGIT = re.compile(r"[0-9]")


def _column_widths(rows, ncol: int) -> list[int]:
    """nt per alignment column, as the engine writes them: 3, or 3 * ceil(N / 3) for a frameshift digit N."""
    widths = [3] * ncol
    for row in rows:
        for i, ch in enumerate(row):
            if ch.isdigit():
                widths[i] = max(widths[i], (max(int(ch) - 1, 0) // 3 + 1) * 3)
    return widths


def build_rows(protein: list[tuple[str, str]], codon: list[tuple[str, str]] | None,
               report: dict | None = None, dna: list[tuple[str, str]] | None = None,
               codontable: int | None = None, widths: list[int] | None = None) -> tuple[list[dict], int, list[str]]:
    """Rows for rendering: [{"id", "status", "cells": [(aa, codon, flag)], "note", "events"}], ncols, notes.
    CHANGED(GUI-10): widths: nt per column computed from every sequence (view() passes it when only some are
    built); default: from the given rows."""
    if codontable is None:
        codontable = int(((report or {}).get("options") or {}).get("codontable") or 1)
    code, starts = gcode.table(codontable)
    marks = _report_marks(report, {p for p, _ in protein})
    cod = dict(codon or [])
    dnad = dict(dna or [])
    notes: list[str] = []
    ncols = max((len(s) for _, s in protein), default=0)
    # CHANGED(FIN-A7): was "any row length != 3 x columns", so a column widened by a frameshift digit >= 4 looked like
    # removed columns ("-nogap or -nomismatch" note although neither was used); now the options actually used, or,
    # without a report, the length expected from the column widths (frameshift digits included).
    # See overnight/integration/CHANGES.md
    if widths is None or len(widths) < ncols:
        widths = _column_widths([row for _, row in protein], ncols)
    opts = (report or {}).get("options") if isinstance((report or {}).get("options"), dict) else None
    if opts is not None:
        codon_only = bool(codon) and bool(opts.get("nogap") or opts.get("nomismatch") or opts.get("blockonly"))
    else:
        codon_only = bool(codon) and any(pid in cod and len(cod[pid]) != sum(widths[:len(row)])
                                         for pid, row in protein)
    if codon_only:
        notes.append("Columns were removed from the codon alignment (-nogap or -nomismatch), so columns are "
                     "numbered as in the output and amino acids are translated from the codons. Problems are "
                     "reported in protein-alignment columns, which no longer line up with this view, so no "
                     "codons are marked here. Run again without -nogap / -nomismatch to see where they are.")
        ncols = max((len(cod[p]) // 3 for p, _ in protein if p in cod), default=0)
    rows = []
    for pid, prow in protein:
        m = marks.get(pid, {"status": None, "mismatch": set(), "frameshift": set()})
        r = {"id": pid, "status": m["status"], "cells": [], "note": "", "events": [], "checked": not codon_only}
        if pid in cod:
            crow = cod[pid].upper()
            n = len(crow) // 3 if codon_only else len(prow)
            first_res = True
            bounds = None if codon_only else [0]
            if bounds is not None:  # FIN-A7: cell i holds widths[i] nt (wider for a frameshift digit >= 4)
                for w in widths[:n]:
                    bounds.append(bounds[-1] + w)
            for i in range(n):
                c = crow[3 * i:3 * i + 3] if bounds is None else crow[bounds[i]:bounds[i + 1]]
                aa = gcode.translate(c, codontable) if codon_only else prow[i].upper()
                if codon_only and "-" in c:
                    aa = "-"
                col = i + 1
                if aa in "-.":
                    r["cells"].append(("-", c or "---", "gap"))
                    continue
                flag = ""
                if codon_only:
                    # CHANGED(APP-1): was report columns applied to output columns and aa translated from the
                    # same codon (self-comparison); now nothing is marked in this mode. See overnight/app/CHANGES.md
                    # report columns are protein-alignment columns and aa comes from this very codon:
                    # nothing here can be checked honestly, so nothing is marked
                    pass
                elif col in m["frameshift"]:
                    flag = "frameshift"
                elif col in m["mismatch"] or ("-" not in c and not _ok(aa, c, code, starts, first_res)):
                    flag = "mismatch"
                first_res = False
                r["cells"].append((aa, c, flag))
        else:
            seq_dna = dnad.get(pid)
            if seq_dna is None and dna and len(dna) == len(protein):
                # COAL pairs by order when IDs differ
                seq_dna = dna[[p for p, _ in protein].index(pid)][1]
            if seq_dna is not None:
                d = diagnose(prow, seq_dna, codontable)
                r["events"] = d["events"]
                fs = {e["column"] for e in d["events"] if e["code"] == "frameshift"} | m["frameshift"]
                for i, a in enumerate(prow):
                    col = i + 1
                    if not a.isalpha():
                        r["cells"].append(("-", "---", "gap"))
                        continue
                    c, flag = d["cells"].get(col, ("   ", "missing"))
                    if col in fs:
                        flag = "frameshift"
                    r["cells"].append((a.upper(), c, flag))
                if d["events"]:
                    r["note"] = d["events"][0]["text"]
                if r["status"] is None:
                    # CHANGED(APP-4): was "failed" if events else "ok" (mismatches ignored); now "warning" too
                    flagged = any(f in ("mismatch", "frameshift") for _, _, f in r["cells"])
                    r["status"] = "failed" if d["events"] else ("warning" if flagged else "ok")
            else:
                for a in prow:
                    r["cells"].append(("-" if not a.isalpha() else a.upper(), "   ",
                                       "gap" if not a.isalpha() else "missing"))
                r["note"] = "No codons: COAL could not align this gene (give the DNA file to see where it breaks)."
        if r["status"] is None:
            r["status"] = "warning" if any(f in ("mismatch", "frameshift") for _, _, f in r["cells"]) else "ok"
        rows.append(r)
    return rows, ncols, notes


# ---------------------------------------------------------------- rendering

# CHANGED(GUI-2): was one light palette with no text colour pinned on headings, cells and links, so in Gradio's dark
# theme its light text leaked into the white panel (light grey on white, nearly invisible); now every colour is a
# variable, set for light, for Gradio's dark theme (the ".dark" body class) and, in the downloaded page, for a dark
# system theme; rules are "div.coalv ..." so they win over Gradio's ".prose ..." styles. Residue tints keep the
# Okabe-Ito hues (dark versions: text 7.1-10.2:1); outlines, hatching and icon + word statuses unchanged. Text colours
# are !important: Gradio's own rules for h2/th/td/a/code must not win in either theme.
# See overnight/integration/CHANGES.md
# CHANGED(GUI-10): ov* colours of the overview (graphics: problem marks >= 3:1 on the background, WCAG 1.4.11)
# CHANGED(AMB-1): statuses follow the app (warning ember #D45400 / #FF6A13, was amber; failed crimson #E0284F, was
# ember) so that none looks like the accent; the jump-target ring and preview edge are the accent (deep gold
# #8A5A00 on light, amber on dark; light was #E69F00 / #B07A00). Marks keep their shapes (tall / short / hanging).
_LIGHT = {"bg": "#ffffff", "fg": "#1a1a1a", "muted": "#5f6b75", "line": "#d5dbe0", "gap": "#6e7780",
          "ovok": "#b8b2ab", "ovwarn": "#D45400", "ovfail": "#E0284F", "ovgap": "#8a939b", "ovprob": "#A83E00",
          "ovbox": "#1a1a1a", "ovun": "#d5d0ca", "prevbg": "#fff6dc", "prevline": "#8A5A00",
          "mark": "#000000", "hl": "#fff3b0", "hlring": "#8A5A00", "link": "#0057a0", "notebg": "#eef3f7",
          "noteline": "#0072B2", "hatch": "rgba(0,0,0,0.12)", "hatch2": "rgba(0,0,0,0.06)",
          "hy": "#d4ebf9", "ar": "#f0d6e6", "po": "#c9dcef", "ne": "#f8d5c2", "pl": "#c9eadf", "sp": "#faf4c0",
          "st": "#d9d9d9"}
_DARK = {"bg": "#1c1a19", "fg": "#ede9e4", "muted": "#b3ada6", "line": "#3a3633", "gap": "#9a938c",
         "ovok": "#6a655f", "ovwarn": "#FF6A13", "ovfail": "#E0284F", "ovgap": "#7d7670", "ovprob": "#FF8A45",
         "ovbox": "#ede9e4", "ovun": "#4a4642", "prevbg": "#2e2614", "prevline": "#FFC24B",
         "mark": "#ede9e4", "hl": "#5c4a10", "hlring": "#FFC24B", "link": "#8cc4ff", "notebg": "#262321",
         "noteline": "#56B4E9", "hatch": "rgba(255,255,255,0.18)", "hatch2": "rgba(255,255,255,0.08)",
         "hy": "#1e3b52", "ar": "#4a2a40", "po": "#2f3157", "ne": "#573420", "pl": "#1d4537", "sp": "#48421c",
         "st": "#3d3b39"}


def _vars(d: dict) -> str:
    return "".join(f"--{k}:{v};" for k, v in d.items())


CSS = ("body.coalv-page { margin:0; padding:12px 16px; background:#ffffff; }"
       f"div.coalv {{ {_vars(_LIGHT)} color-scheme:light; }}"
       f".dark div.coalv {{ {_vars(_DARK)} color-scheme:dark; }}"
       f"@media (prefers-color-scheme: dark) {{ body.coalv-page div.coalv {{ {_vars(_DARK)} color-scheme:dark; }}"
       f" body.coalv-page {{ background:{_DARK['bg']}; }} }}") + """
div.coalv { background:var(--bg); color:var(--fg) !important; font-family:system-ui,"Segoe UI",Helvetica,Arial,sans-serif;
  font-size:16px; line-height:1.35; padding:0.2em 0.9em 0.9em; border-radius:6px; }
div.coalv h2 { font-size:1.25em; margin:0.6em 0 0.3em; color:var(--fg) !important; }
div.coalv p, div.coalv li, div.coalv span, div.coalv strong { color:inherit !important; }
div.coalv code { background:none; color:inherit !important; padding:0; font-size:0.95em; }
div.coalv .note { background:var(--notebg); border-left:5px solid var(--noteline); padding:0.4em 0.7em;
  margin:0.4em 0; color:var(--fg) !important; }
div.coalv .legend { display:flex; flex-wrap:wrap; gap:0.4em 1em; margin:0.5em 0; font-size:0.9em; }
div.coalv .legend span { display:inline-flex; align-items:center; gap:0.35em; }
div.coalv .sw { display:inline-block; width:1.1em; height:1.1em; border:1px solid var(--line); border-radius:3px; }
div.coalv .issues { margin:0.3em 0 0.8em; padding-left:1.2em; }
div.coalv .issues li { margin:0.25em 0; }
div.coalv a { color:var(--link) !important; text-decoration:underline; }
div.coalv a:focus-visible, div.coalv .scroll:focus-visible { outline:3px solid var(--hlring); outline-offset:2px; }
div.coalv .scroll { overflow:auto; max-height:75vh; border:1px solid var(--line); border-radius:6px;
  background:var(--bg); }
div.coalv table { border-collapse:separate; border-spacing:0; font-family:ui-monospace,Consolas,Menlo,monospace;
  font-size:15px; margin:0; width:auto; background:var(--bg); color:var(--fg) !important; }
div.coalv thead th { position:sticky; top:0; z-index:2; background:var(--bg); color:var(--muted) !important;
  font-weight:400; font-size:12px; text-align:left; padding:2px 0 2px 2px; border:0;
  border-bottom:1px solid var(--line); height:1.4em; white-space:nowrap; }
div.coalv thead th.flag { color:var(--mark) !important; font-weight:700; }
div.coalv th.name { position:sticky; left:0; z-index:1; background:var(--bg); color:var(--fg) !important; text-align:left;
  font-family:system-ui,"Segoe UI",Helvetica,Arial,sans-serif; font-size:14px; font-weight:600;
  padding:0 0.6em; border:0; border-right:2px solid var(--line); white-space:nowrap; max-width:16em;
  overflow:hidden; text-overflow:ellipsis; }
div.coalv thead th.name { z-index:3; }
div.coalv th.name .st { font-weight:700; margin-right:0.3em; }
div.coalv td { padding:1px 2px; text-align:center; min-width:2.6em; scroll-margin:3em 3em 3em 18em; border:0;
  color:var(--fg) !important; }
div.coalv td b { display:block; font-size:12px; line-height:1.1; font-weight:700; color:inherit !important; }
div.coalv td i { display:block; font-style:normal; letter-spacing:0.02em; color:inherit !important; }
div.coalv td.gap { color:var(--gap) !important; }
div.coalv td.hy{background:var(--hy)} div.coalv td.ar{background:var(--ar)} div.coalv td.po{background:var(--po)}
div.coalv td.ne{background:var(--ne)} div.coalv td.pl{background:var(--pl)} div.coalv td.sp{background:var(--sp)}
div.coalv td.st{background:var(--st)}
div.coalv td.mm { outline:2px solid var(--mark); outline-offset:-2px; }
div.coalv td.mm b::before { content:"\\2260 "; }
div.coalv td.fs { outline:3px double var(--mark); outline-offset:-3px; background-image:repeating-linear-gradient(
  45deg, transparent 0 4px, var(--hatch) 4px 6px); }
div.coalv td.fs b::before { content:"\\21AF "; }
div.coalv td.miss { color:var(--gap) !important; background:repeating-linear-gradient(-45deg, transparent 0 4px,
  var(--hatch2) 4px 6px); }
div.coalv td:target, div.coalv td.hl { background:var(--hl) !important; box-shadow:inset 0 0 0 3px var(--hlring); }
div.coalv tr:target th.name { background:var(--hl); }
div.coalv .gapkey { color:var(--gap) !important; }
div.coalv .rownote { font-weight:400; font-size:12px; color:var(--muted) !important; }
div.coalv .preview { background:var(--prevbg); border-left:5px solid var(--prevline); padding:0.45em 0.8em;
  margin:0.4em 0 0.6em; border-radius:0 6px 6px 0; }
div.coalv .preview p { margin:0.2em 0; color:var(--fg) !important; }
div.coalv .outside { color:var(--muted) !important; }
div.coalv details.hidden-issues { margin:0 0 0.8em; }
div.coalv details.hidden-issues summary { cursor:pointer; color:var(--fg) !important; }
div.coalv .overview { border:1px solid var(--line); border-radius:6px; padding:0.4em 0.8em 0.5em; margin:0.4em 0 0.7em; }
div.coalv .ov-title { margin:0.1em 0 0.3em; }
div.coalv .ov-sum { color:var(--muted) !important; font-size:0.92em; }
div.coalv .ov-row { display:grid; grid-template-columns:6.2em minmax(0, 1fr); column-gap:0.6em; align-items:center;
  margin:0.3em 0; }
div.coalv .ov-lab { grid-row:1; font-size:0.85em; color:var(--muted) !important; }
div.coalv svg.ov { grid-column:2; grid-row:1; display:block; width:100%%; background:transparent; }
div.coalv svg.ov-seq { height:25px; } div.coalv svg.ov-col { height:40px; }
div.coalv .ov-ends { grid-column:2; display:flex; justify-content:space-between; font-size:11px; line-height:1.2;
  color:var(--muted) !important; }
div.coalv .ov-ok { fill:var(--ovok); } div.coalv .ov-warning { fill:var(--ovwarn); }
div.coalv .ov-failed { fill:var(--ovfail); } div.coalv .ov-not_checked { fill:var(--ovun); }
div.coalv .ov-drawn { fill:var(--fg); } div.coalv .ov-gap { fill:var(--ovgap); }
div.coalv .ov-prob { fill:var(--ovprob); }
div.coalv .ov-box { fill:none; stroke:var(--ovbox); stroke-width:2px; }
div.coalv .ov-key { display:flex; flex-wrap:wrap; gap:0.25em 1.1em; margin:0.4em 0 0; font-size:0.82em;
  color:var(--muted) !important; }
div.coalv .ov-key span { display:inline-flex; align-items:center; gap:0.35em; }
div.coalv i.k { display:inline-block; width:0.7em; height:1em; box-sizing:border-box; }
div.coalv i.k-fail { background:var(--ovfail); } div.coalv i.k-warn { background:var(--ovwarn); height:0.75em;
  align-self:flex-end; margin-bottom:0.2em; }
div.coalv i.k-drawn { background:var(--fg); height:0.3em; }
div.coalv i.k-prob { background:var(--ovprob); height:0.5em; align-self:flex-start; margin-top:0.25em; }
div.coalv i.k-gap { background:var(--ovgap); height:0.6em; align-self:flex-end; margin-bottom:0.2em; }
div.coalv i.k-box { border:2px solid var(--ovbox); width:1.1em; }
@media (max-width: 600px) { div.coalv { padding:0.2em 0.5em 0.6em; }
  div.coalv th.name { max-width:7.5em; padding:0 0.3em; font-size:12px; }
  div.coalv .gapkey { color:var(--gap) !important; }
div.coalv .rownote { display:none; } div.coalv td { scroll-margin-left:9em; } }
@media print { div.coalv { %s } }
""" % _vars(_LIGHT)

JS = """
(function(){
  function mark(){
    document.querySelectorAll('div.coalv td.hl').forEach(function(e){e.classList.remove('hl');});
    var id = decodeURIComponent(location.hash.slice(1)); if(!id) return;
    var t = document.getElementById(id); if(!t || t.tagName !== 'TD') return;
    var idx = t.cellIndex, tb = t.closest('table');
    for (var i = 0; i < tb.rows.length; i++) { var c = tb.rows[i].cells[idx]; if (c && c.tagName === 'TD') c.classList.add('hl'); }
  }
  window.addEventListener('hashchange', mark); mark();
})();
"""

_ICON = {"ok": ("✔", "OK"), "warning": ("⚠", "Warning"), "failed": ("✖", "Failed"),
         "not_checked": ("–", "Not checked")}


ISSUES_OPEN = 8  # CHANGED(GUI-10): problems listed open above the grid; more are under a disclosure
HIDDEN_ISSUES_MAX = 200  # CHANGED(R2-A3): problems of sequences beyond the drawn rows, listed from the report


def _issues(rows: list[dict], report: dict | None, window: tuple[int, int], drawn: set | None = None,
            preview: bool = False) -> str:
    """CHANGED(GUI-10): window (first, last) drawn columns, was the count from column 1; problems of sequences not
    drawn are under a closed "not in the preview" disclosure when the view is a preview."""
    c0, c1 = window
    sentences = {}
    hidden = []
    if report:
        from . import core
        sentences = {s.get("id"): core.sequence_sentence(s) for s in report.get("sequences", [])
                     if s.get("status") not in ("ok", None)}
        if drawn is not None:
            hidden = [s for s in report.get("sequences", []) if s.get("id") not in drawn
                      and s.get("status") in ("failed", "warning")]
    items = []
    unchecked = 0
    for r in rows:
        flagged = [(i + 1, f) for i, (_, _, f) in enumerate(r["cells"]) if f in ("mismatch", "frameshift")]
        if r["status"] == "not_checked" and not flagged and not r["events"]:
            unchecked += 1
            continue
        if r["status"] == "ok" and not flagged and not r["events"]:
            continue
        icon, word = _ICON.get(r["status"], ("?", str(r["status"])))
        text = sentences.get(r["id"], "")
        ev = " ".join(e["text"] for e in r["events"])
        links = []
        for col, f in flagged[:20]:
            label = f"column {col:,}" + (" (frameshift)" if f == "frameshift" else "")
            if c0 <= col <= c1:
                links.append(f'<a href="#{cell_id(r["id"], col)}">{label}</a>')
            else:
                links.append(f'{label} <span class="outside">(outside the drawn columns {c0:,}–{c1:,})</span>')
        more = f" and {len(flagged) - 20} more" if len(flagged) > 20 else ""
        jump = f" Jump to: {', '.join(links)}{more}." if links else ""
        items.append(f'<li><strong>{icon} {word}</strong> <a href="#{row_id(r["id"])}"><code>'
                     f'{html.escape(r["id"])}</code></a>: {html.escape(text or ev)}'
                     f'{(" " + html.escape(ev)) if ev and text and ev not in text else ""}{jump}</li>')
    hidden_items = []
    for sq in hidden[:HIDDEN_ISSUES_MAX]:
        icon, word = _ICON.get(sq.get("status"), ("?", str(sq.get("status"))))
        hidden_items.append(f'<li><strong>{icon} {word}</strong> <code>{html.escape(str(sq.get("id")))}</code> (not '
                            f'drawn below): {html.escape(sentences.get(sq.get("id"), ""))}</li>')
    if len(hidden) > HIDDEN_ISSUES_MAX:
        hidden_items.append(f"<li>and {len(hidden) - HIDDEN_ISSUES_MAX:,} more sequences with problems (the report "
                            f"JSON lists every sequence)</li>")
    hidden_html = ""
    if hidden_items and preview:
        hidden_html = (f'<details class="hidden-issues"><summary>{len(hidden):,} more sequence'
                       f'{"s" if len(hidden) != 1 else ""} with problems are not in the preview</summary>'
                       f'<ul class="issues">{"".join(hidden_items)}</ul></details>')
    else:
        items += hidden_items
    if report and drawn is not None:
        unchecked += sum(1 for sq in report.get("sequences", []) if sq.get("id") not in drawn
                         and sq.get("status") == "not_checked")
    if unchecked:
        items.append(f"<li><strong>– Not checked</strong>: {unchecked} other sequence"
                     f"{'s' if unchecked != 1 else ''} (COAL stopped at an earlier error). In the view below their "
                     f"codons match their amino acids.</li>")
    if not items and not hidden_html:
        if any(not r.get("checked", True) for r in rows):  # CHANGED(APP-1): no "every codon matches" claim
            return ('<p class="note">Codons could not be checked against the protein in this view (columns were '
                    'removed).</p>')
        return '<p class="note">✔ Every codon matches its amino acid.</p>'
    if len(items) > ISSUES_OPEN:  # CHANGED(GUI-10): a long list would push the grid off screen; the rest on demand
        rest = items[ISSUES_OPEN:]
        return ('<ul class="issues">' + "".join(items[:ISSUES_OPEN]) + "</ul>"
                f'<details class="hidden-issues"><summary>{len(rest):,} more drawn sequence'
                f'{"s" if len(rest) != 1 else ""} with problems</summary><ul class="issues">{"".join(rest)}</ul>'
                '</details>' + hidden_html)
    return ('<ul class="issues">' + "".join(items) + "</ul>" if items else "") + hidden_html


def preview_sentence(shown: int, total: int, window: tuple[int, int], ncols: int, selection: dict | None) -> str:
    """CHANGED(GUI-10): 'Preview: 100 of 21,800 sequences, problems first, columns 1–1,000 of 12,000.'"""
    sel = selection or {}
    if shown >= total:
        rows = f"all {total:,} sequence{'s' if total != 1 else ''}"
    else:
        rows = f"{shown:,} of {total:,} sequences"
        if sel.get("problems_first"):
            p, k, c = sel.get("problems") or 0, sel.get("problems_shown") or 0, sel.get("context") or 0
            what = f"all {p:,} with problems" if k == p else f"{k:,} of the {p:,} with problems"
            rows += (f", problems first ({what}, then {c:,} OK ones from across the file)" if c
                     else f", problems first ({what})")
        elif sel.get("context"):
            rows += ", spread evenly through the file"
        else:
            rows += " (the first ones in the file)"
    c0, c1 = window
    cols = (f"all {ncols:,} columns" if c0 == 1 and c1 >= ncols else f"columns {c0:,}–{c1:,} of {ncols:,}")
    return f"Preview: {rows}; {cols}."


FULL_VIEW_HELP = ("To see every sequence and column, open the codon alignment file COAL wrote (Download codon "
                  "alignment) in a desktop alignment viewer such as Jalview or AliView (both free, made for large "
                  "alignments: File, then Open).")


def render_fragment(rows: list[dict], ncols: int, notes: list[str] | None = None, report: dict | None = None,
                    title: str = "Codon alignment", max_seqs: int = MAX_SEQS, max_cols: int = MAX_COLS,
                    total_rows: int | None = None, col_start: int = 1, max_cells: int = MAX_CELLS,
                    overview: str = "", selection: dict | None = None) -> str:
    """HTML fragment (a <style> plus a <div class="coalv">) for embedding in another page.
    total_rows: how many sequences the alignment has when `rows` holds only some (see view()).
    CHANGED(GUI-10): col_start (first column drawn), max_cells, the overview block (overview.overview_html) and
    the row selection (select_rows) for the "Preview: ..." line; was "Large alignment: showing the first ..."."""
    notes = list(notes or [])
    shown_rows = rows[:max_seqs]
    total = max(len(rows), total_rows or 0)
    cols = preview_columns(len(shown_rows), max(0, ncols - col_start + 1), max_cols, max_cells)
    c0 = max(1, col_start)
    c1 = c0 + cols - 1
    is_preview = total > len(shown_rows) or c0 > 1 or c1 < ncols
    out = [f"<style>{CSS}</style>", '<div class="coalv">', f"<h2>{html.escape(title)}</h2>"]
    if is_preview:
        out.append(f'<div class="preview" role="note"><p><strong>{html.escape(preview_sentence(len(shown_rows), total, (c0, c1), ncols, selection))}'
                   f'</strong> Drawing all of it would make the page slow, so only this part is drawn; problems '
                   f'outside it are still listed below.</p><p>{html.escape(FULL_VIEW_HELP)}</p></div>')
    out.append(overview)
    for n in notes:
        out.append(f'<p class="note">{html.escape(n)}</p>')
    out.append(_issues(shown_rows, report, (c0, c1), {r["id"] for r in shown_rows}, is_preview))
    out.append('<div class="legend" aria-label="Colour key">')
    for cls, label in LEGEND:
        out.append(f'<span><span class="sw" style="background:var(--{cls})"></span>{label}</span>')
    out.append('<span><b>≠</b> codon does not match its amino acid (outlined)</span>'
               '<span><b>↯</b> frameshift (double outline, hatched)</span>'
               '<span class="gapkey">--- gap</span></div>')  # CHANGED(GUI-2): was an inline colour
    flagged_cols = {c0 + i for r in shown_rows for i, (_, _, f) in enumerate(r["cells"][c0 - 1:c1])
                    if f in ("mismatch", "frameshift")}
    out.append('<div class="scroll" tabindex="0" aria-label="Codon alignment, scrolls sideways">'
               '<table><thead><tr><th class="name" scope="col">Sequence</th>')
    for c in range(c0, c1 + 1):
        if c in flagged_cols:
            out.append(f'<th class="flag" id="coalv-c{c}" title="Column {c}: a problem in at least one sequence">'
                       f'{c}▾</th>')
        elif c == c0 or c % 10 == 0:
            out.append(f"<th>{c}</th>")
        else:
            out.append("<th></th>")
    out.append("</tr></thead><tbody>")
    for r in shown_rows:
        icon, word = _ICON.get(r["status"], ("?", ""))
        note = f'<br><span class="rownote">{html.escape(r["note"])}</span>' if r["note"] else ""
        out.append(f'<tr id="{row_id(r["id"])}"><th class="name" scope="row" title="{html.escape(r["id"])}">'
                   f'<span class="st" aria-label="{word}">{icon}</span>{html.escape(r["id"])}{note}</th>')
        for i, (aa, c, flag) in enumerate(r["cells"][c0 - 1:c1]):
            col = c0 + i
            if flag == "gap":
                out.append(f'<td class="gap"><b>-</b><i>{html.escape(c)}</i></td>')
                continue
            cls = [GROUPS.get(aa, "")]
            attr = ""
            if flag == "mismatch":
                cls.append("mm")
                attr = (f' id="{cell_id(r["id"], col)}" title="Column {col}: {aa} is not coded by '
                        f'{html.escape(c.strip())}"')
            elif flag == "frameshift":
                cls.append("fs")
                attr = f' id="{cell_id(r["id"], col)}" title="Column {col}: frameshift starts here"'
            elif flag == "missing":
                cls = ["miss"]
            out.append(f'<td class="{" ".join(x for x in cls if x)}"{attr}><b>{aa if aa.isalpha() else html.escape(aa)}</b>'
                       f'<i>{(c if c.isalpha() else html.escape(c)) if c.strip() else "&nbsp;"}</i></td>')
        out.append("</tr>")
    out.append("</tbody></table></div></div>")
    return "".join(out)


def render_page(fragment: str, title: str = "COAL alignment view") -> str:
    """A complete, self-contained HTML document (no external scripts or fonts)."""
    return ("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"
            f"<title>{html.escape(title)}</title><meta name=\"color-scheme\" content=\"light dark\">"
            f"</head><body class=\"coalv-page\">{fragment}<script>{JS}</script></body></html>")


def view(protein_path: str, codon_path: str | None = None, report: dict | str | None = None,
         dna_path: str | None = None, codontable: int | None = None, page: bool = True,
         max_seqs: int = MAX_SEQS, max_cols: int = MAX_COLS, max_cells: int = MAX_CELLS,
         centre: int | None = None, info: dict | None = None) -> str:
    """Render a COAL run as HTML. `report` may be a dict, a path to a coal.report/1 file, or None.
    CHANGED(GUI-10): an honest preview at scale: rows chosen by select_rows (problems first, then OK ones from
    across the file), a column window (centre: the column to centre it on; by default column 1, or the first
    problem of a drawn sequence when that is further right), an overview of the whole alignment above the grid
    (overview.py) and max_cells. info (a dict, optional) is filled with what was drawn: preview, ids, first,
    last, ncols, total, shown. See overnight/integration/CHANGES.md"""
    from . import overview as ov
    if isinstance(report, str):
        with open(report, encoding="utf-8") as fh:
            report = json.load(fh)
    protein = read_alignment(protein_path)
    total = len(protein)
    ids = [p for p, _ in protein]
    ncols_all = max((len(r) for _, r in protein), default=0)
    # CHANGED(R2-A3): was rows built for every sequence (22,000 x 600 codons: 29 s, GBs) though only max_seqs
    # are drawn; now only the drawn sequences are built (problems of the others come from the report).
    prefer = None
    if centre is not None and report:
        # centred on a column: sequences with a problem in that window are drawn first
        cw0, cw1 = column_window(ncols_all, preview_columns(min(max_seqs, total), ncols_all, max_cols, max_cells),
                                 centre)
        prefer = {s.get("id") for s in report.get("sequences") or [] if any(
            isinstance(m.get("column"), int) and cw0 <= m["column"] <= cw1 for m in s.get("messages") or [])}
    chosen, selection = select_rows(ids, report, max_seqs, prefer)
    want = set(chosen)
    prot = dict(protein)
    protein_sel = [(p, prot[p]) for p in chosen]
    # column widths from every row (a frameshift digit >= 4 in any sequence widens that column for all of them)
    digits = _DIGIT.search("".join(r for _, r in protein)) is not None
    widths = _column_widths([r for _, r in protein if _DIGIT.search(r)] if digits else [], ncols_all)
    codon = (read_alignment(codon_path, want=want) if codon_path and os.path.isfile(codon_path) else None)
    have = {i for i, _ in codon} if codon is not None else set()
    need = [p for p in chosen if p not in have]  # drawn sequences COAL wrote no codons for
    dna = None
    if dna_path and os.path.isfile(dna_path) and (need or codon is None):
        dna = read_alignment(dna_path, dna=True, want=set(need))  # CHANGED(R3-A1)
        if len(dna) < len(need):
            full = read_alignment(dna_path, dna=True)
            if len(full) == total and not all(p in {i for i, _ in full} for p in need):
                # COAL pairs by order when IDs differ: give each drawn protein the DNA at its own position
                pos = {p: k for k, p in enumerate(ids)}
                dna = [(p, full[pos[p]][1]) for p in need]
            else:
                dna = [(i, r) for i, r in full if i in want]
    rows, ncols, notes = build_rows(protein_sel, codon, report, dna, codontable, widths=widths)
    if codon is None:
        notes.insert(0, "COAL wrote no codon alignment for this run. Codons below were placed by the viewer, "
                        "only to show where each gene breaks." if dna else
                     "COAL wrote no codon alignment for this run, and no DNA file was given.")
    cols = preview_columns(len(rows), ncols, max_cols, max_cells)
    if centre is None and cols < ncols:
        first_problem = min((i + 1 for r in rows for i, (_, _, f) in enumerate(r["cells"])
                             if f in ("mismatch", "frameshift")), default=None)
        if first_problem is not None and first_problem > cols:
            centre = first_problem
    c0, c1 = column_window(ncols, cols, centre)
    is_preview = total > len(rows) or c0 > 1 or c1 < ncols
    codon_only = any(not r.get("checked", True) for r in rows)
    over = ""
    if is_preview:
        over = ov.overview_html(ids, report, [r for _, r in protein], ncols_all, want, (c0, c1),
                                columns=not codon_only)
    if info is not None:
        info.update({"preview": is_preview, "ids": want, "first": c0, "last": c1, "ncols": ncols,
                     "total": total, "shown": len(rows)})
    title = f"Codon alignment: {os.path.basename(protein_path)}"
    frag = render_fragment(rows, ncols, notes, report, title, max_seqs, max_cols, total_rows=total, col_start=c0,
                           max_cells=max_cells, overview=over, selection=selection)
    return render_page(frag, title) if page else frag


def render_run(protein_path: str, output_path: str | None, report: dict, dna_path: str | None, codontable: int,
               out_dir: str, limits: tuple[int, int, int] | None = None, centre: int | None = None,
               info: dict | None = None) -> tuple[str, str]:
    """CHANGED(R2-A3): new; was two full view() calls in ui.run_ui. Renders once: (fragment for the GUI, path of
    the downloadable page written into out_dir).
    CHANGED(GUI-10): limits (sequences, columns, cells: the preview size from the settings), centre (a column to
    centre the preview on) and info (filled by view()); the downloadable page is the same preview."""
    # CHANGED(FIN-A4): was the output read as FASTA/CLUSTAL whatever its format, so PAML and codon output gave
    # "<server path>: not FASTA or CLUSTAL"; now those are drawn from a FASTA copy of the same run (same options),
    # written to a hidden temporary file in out_dir and deleted. See overnight/integration/CHANGES.md
    opts = (report or {}).get("options") or {}
    tmp = None
    if output_path and opts.get("output_format") in ("paml", "codon"):
        from . import core
        tmp = os.path.join(out_dir, ".view_copy.fasta")
        code, _, _ = core._run_engine(protein_path, dna_path, tmp, "fasta", codontable, bool(opts.get("nogap")),
                                      bool(opts.get("nomismatch")), False, None,
                                      blockonly=bool(opts.get("blockonly")))
        output_path = tmp if code == 0 and os.path.isfile(tmp) else None
    ms, mc, mx = limits or (MAX_SEQS, MAX_COLS, MAX_CELLS)
    try:
        frag = view(protein_path, output_path, report, dna_path, codontable, page=False, max_seqs=ms, max_cols=mc,
                    max_cells=mx, centre=centre, info=info)
    finally:
        if tmp and os.path.exists(tmp):
            os.remove(tmp)
    path = os.path.join(out_dir, "alignment_view.html")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(render_page(frag, f"Codon alignment: {os.path.basename(protein_path)}"))
    return frag, path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m coal_app.viewer",
                                 description="Render a COAL codon alignment as a self-contained HTML page.")
    ap.add_argument("protein", help="protein alignment given to COAL (FASTA or CLUSTAL)")
    ap.add_argument("codon", nargs="?", help="codon alignment written by COAL (FASTA); omit if COAL failed")
    ap.add_argument("--report", help="coal.report/1 JSON (from coal --report or the GUI)")
    ap.add_argument("--dna", help="the DNA FASTA given to COAL: lets the viewer show where failed genes break")
    ap.add_argument("--codontable", type=int, help="genetic code (default: from the report, else 1)")
    ap.add_argument("--max-seqs", type=int, default=MAX_SEQS)
    ap.add_argument("--max-cols", type=int, default=MAX_COLS)
    ap.add_argument("-o", "--out", required=True, help="output HTML file")
    a = ap.parse_args(argv)
    report = a.report
    codon = a.codon
    if codon is None and report is None and a.dna:
        # CHANGED(R2-A12): was every gene diagnosed by the viewer alone (it called engine-aligned genes "Failed");
        # now COAL itself is run on the inputs, and the viewer only diagnoses genes COAL rejected.
        from . import core
        try:
            r = core.run_coal(a.protein, a.dna, codontable=a.codontable or 1, probe_codes=False)
            codon, report = r.output_path, r.report
        except (OSError, ValueError) as e:
            print(f"coal viewer: note: could not run COAL ({e}); showing the viewer's own placement",
                  file=sys.stderr)
    try:
        page = view(a.protein, codon, report, a.dna, a.codontable, True, a.max_seqs, a.max_cols)
    except (OSError, ValueError) as e:
        print(f"coal viewer: error: {e}", file=sys.stderr)
        return 1
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(page)
    print(f"wrote {a.out}")
    return 0


def locate_breaks(protein_path: str, dna_path: str, report: dict) -> None:
    """Add the column where each failed gene breaks (frameshift or not) to a report, in place."""
    try:
        protein = read_alignment(protein_path)
        dna = read_alignment(dna_path, dna=True)  # CHANGED(R3-A1)
    except (OSError, ValueError):
        return
    table = int((report.get("options") or {}).get("codontable") or 1)
    prot = dict(protein)
    dnad = dict(dna)
    order = [p for p, _ in protein]
    for s in report.get("sequences", []):
        if s.get("status") != "failed" or any(isinstance(m.get("column"), int) for m in s.get("messages") or []):
            continue
        pid = s.get("id")
        if pid not in prot:
            continue
        d = dnad.get(s.get("dna_id") or pid)
        if d is None and len(dna) == len(protein):
            d = dna[order.index(pid)][1]
        if d is None:
            continue
        ev = diagnose(prot[pid], d, table)["events"]
        if ev:
            s["messages"] = [{k: v for k, v in e.items()} for e in ev[:1]] + list(s.get("messages") or [])


if __name__ == "__main__":
    sys.exit(main())
