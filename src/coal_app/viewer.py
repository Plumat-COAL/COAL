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

def read_alignment(path: str) -> list[tuple[str, str]]:
    """(id, sequence) records from FASTA or CLUSTAL, in file order."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    stripped = text.lstrip()
    recs: list[tuple[str, str]] = []
    if stripped.startswith(">"):
        name, buf = None, []
        for line in text.splitlines():
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                parts = line[1:].split()
                name, buf = (parts[0] if parts else ""), []
            elif name is not None:
                buf.append(line.strip())
        if name is not None:
            recs.append((name, "".join(buf)))
        return recs
    if stripped.upper().startswith("CLUSTAL"):
        seqs: dict[str, list[str]] = {}
        for line in text.splitlines()[1:]:
            if not line.strip() or line[0].isspace():
                continue
            parts = line.split()
            if len(parts) >= 2:
                seqs.setdefault(parts[0], []).append(parts[1])
        return [(k, "".join(v)) for k, v in seqs.items()]
    raise ValueError(f"{path}: not FASTA or CLUSTAL")


def slug(seq_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]", "_", seq_id)


def cell_id(seq_id: str, col: int) -> str:
    return f"coalv-{slug(seq_id)}-c{col}"


def row_id(seq_id: str) -> str:
    return f"coalv-{slug(seq_id)}"


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

def _report_marks(report: dict | None) -> dict[str, dict]:
    marks: dict[str, dict] = {}
    for s in (report or {}).get("sequences", []):
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


def build_rows(protein: list[tuple[str, str]], codon: list[tuple[str, str]] | None,
               report: dict | None = None, dna: list[tuple[str, str]] | None = None,
               codontable: int | None = None) -> tuple[list[dict], int, list[str]]:
    """Rows for rendering: [{"id", "status", "cells": [(aa, codon, flag)], "note", "events"}], ncols, notes."""
    if codontable is None:
        codontable = int(((report or {}).get("options") or {}).get("codontable") or 1)
    code, starts = gcode.table(codontable)
    marks = _report_marks(report)
    cod = dict(codon or [])
    dnad = dict(dna or [])
    notes: list[str] = []
    ncols = max((len(s) for _, s in protein), default=0)
    codon_only = bool(codon) and any(
        pid in cod and len(cod[pid]) != 3 * len(row) for pid, row in protein)
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
            for i in range(n):
                c = crow[3 * i:3 * i + 3]
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

CSS = """
.coalv { --bg:#ffffff; --fg:#1a1a1a; --muted:#5f6b75; --line:#d5dbe0; --gap:#9aa4ad; --mark:#000000;
  --hl:#fff3b0; --hy:#d4ebf9; --ar:#f0d6e6; --po:#c9dcef; --ne:#f8d5c2; --pl:#c9eadf; --sp:#faf4c0; --st:#d9d9d9;
  background:var(--bg); color:var(--fg); font-family:system-ui,"Segoe UI",Helvetica,Arial,sans-serif;
  font-size:16px; line-height:1.35; }
.coalv h2 { font-size:1.25em; margin:0.6em 0 0.3em; }
.coalv .note { background:#eef3f7; border-left:5px solid #0072B2; padding:0.4em 0.7em; margin:0.4em 0; }
.coalv .legend { display:flex; flex-wrap:wrap; gap:0.4em 1em; margin:0.5em 0; font-size:0.9em; }
.coalv .legend span { display:inline-flex; align-items:center; gap:0.35em; }
.coalv .sw { display:inline-block; width:1.1em; height:1.1em; border:1px solid var(--line); border-radius:3px; }
.coalv .issues { margin:0.3em 0 0.8em; padding-left:1.2em; }
.coalv .issues li { margin:0.25em 0; }
.coalv .issues a { color:#0057a0; }
.coalv .scroll { overflow:auto; max-height:75vh; border:1px solid var(--line); border-radius:6px; }
.coalv table { border-collapse:separate; border-spacing:0; font-family:ui-monospace,Consolas,Menlo,monospace;
  font-size:15px; }
.coalv thead th { position:sticky; top:0; z-index:2; background:var(--bg); color:var(--muted);
  font-weight:400; font-size:12px; text-align:left; padding:2px 0 2px 2px; border-bottom:1px solid var(--line);
  height:1.4em; white-space:nowrap; }
.coalv thead th.flag { color:var(--mark); font-weight:700; }
.coalv th.name { position:sticky; left:0; z-index:1; background:var(--bg); text-align:left;
  font-family:system-ui,"Segoe UI",Helvetica,Arial,sans-serif; font-size:14px; font-weight:600;
  padding:0 0.6em; border-right:2px solid var(--line); white-space:nowrap; max-width:16em;
  overflow:hidden; text-overflow:ellipsis; }
.coalv thead th.name { z-index:3; }
.coalv th.name .st { font-weight:700; margin-right:0.3em; }
.coalv td { padding:1px 2px; text-align:center; min-width:2.6em; scroll-margin:3em 3em 3em 18em; }
.coalv td b { display:block; font-size:12px; line-height:1.1; }
.coalv td i { display:block; font-style:normal; letter-spacing:0.02em; }
.coalv td.gap { color:var(--gap); }
.coalv td.hy{background:var(--hy)} .coalv td.ar{background:var(--ar)} .coalv td.po{background:var(--po)}
.coalv td.ne{background:var(--ne)} .coalv td.pl{background:var(--pl)} .coalv td.sp{background:var(--sp)}
.coalv td.st{background:var(--st)}
.coalv td.mm { outline:2px solid var(--mark); outline-offset:-2px; }
.coalv td.mm b::before { content:"\\2260 "; }
.coalv td.fs { outline:3px double var(--mark); outline-offset:-3px; background-image:repeating-linear-gradient(
  45deg, transparent 0 4px, rgba(0,0,0,0.12) 4px 6px); }
.coalv td.fs b::before { content:"\\21AF "; }
.coalv td.miss { color:var(--gap); background:repeating-linear-gradient(-45deg, transparent 0 4px,
  rgba(0,0,0,0.06) 4px 6px); }
.coalv td:target, .coalv td.hl { background:var(--hl) !important; box-shadow:inset 0 0 0 3px #E69F00; }
.coalv tr:target th.name { background:var(--hl); }
.coalv .rownote { font-weight:400; font-size:12px; color:var(--muted); }
@media (max-width: 600px) { .coalv th.name { max-width:7.5em; padding:0 0.3em; font-size:12px; }
  .coalv .rownote { display:none; } .coalv td { scroll-margin-left:9em; } }
"""

JS = """
(function(){
  function mark(){
    document.querySelectorAll('.coalv td.hl').forEach(function(e){e.classList.remove('hl');});
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


HIDDEN_ISSUES_MAX = 200  # CHANGED(R2-A3): problems of sequences beyond the drawn rows, listed from the report


def _issues(rows: list[dict], report: dict | None, shown_cols: int, drawn: set | None = None) -> str:
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
            label = f"column {col}" + (" (frameshift)" if f == "frameshift" else "")
            if col <= shown_cols:
                links.append(f'<a href="#{cell_id(r["id"], col)}">{label}</a>')
            else:
                links.append(f"{label} (beyond the shown columns)")
        more = f" and {len(flagged) - 20} more" if len(flagged) > 20 else ""
        jump = f" Jump to: {', '.join(links)}{more}." if links else ""
        items.append(f'<li><strong>{icon} {word}</strong> <a href="#{row_id(r["id"])}"><code>'
                     f'{html.escape(r["id"])}</code></a>: {html.escape(text or ev)}'
                     f'{(" " + html.escape(ev)) if ev and text and ev not in text else ""}{jump}</li>')
    for sq in hidden[:HIDDEN_ISSUES_MAX]:
        icon, word = _ICON.get(sq.get("status"), ("?", str(sq.get("status"))))
        items.append(f'<li><strong>{icon} {word}</strong> <code>{html.escape(str(sq.get("id")))}</code> (not drawn '
                     f'below): {html.escape(sentences.get(sq.get("id"), ""))}</li>')
    if len(hidden) > HIDDEN_ISSUES_MAX:
        items.append(f"<li>and {len(hidden) - HIDDEN_ISSUES_MAX:,} more sequences with problems (see the report "
                     f"JSON or the table above)</li>")
    if report and drawn is not None:
        unchecked += sum(1 for sq in report.get("sequences", []) if sq.get("id") not in drawn
                         and sq.get("status") == "not_checked")
    if unchecked:
        items.append(f"<li><strong>– Not checked</strong>: {unchecked} other sequence"
                     f"{'s' if unchecked != 1 else ''} (COAL stopped at an earlier error). In the view below their "
                     f"codons match their amino acids.</li>")
    if not items:
        if any(not r.get("checked", True) for r in rows):  # CHANGED(APP-1): no "every codon matches" claim
            return ('<p class="note">Codons could not be checked against the protein in this view (columns were '
                    'removed).</p>')
        return '<p class="note">✔ Every codon matches its amino acid.</p>'
    return '<ul class="issues">' + "".join(items) + "</ul>"


def render_fragment(rows: list[dict], ncols: int, notes: list[str] | None = None, report: dict | None = None,
                    title: str = "Codon alignment", max_seqs: int = MAX_SEQS, max_cols: int = MAX_COLS,
                    total_rows: int | None = None) -> str:
    """HTML fragment (a <style> plus a <div class="coalv">) for embedding in another page.
    total_rows: how many sequences the alignment has when `rows` holds only the first ones (see view())."""
    notes = list(notes or [])
    shown_rows = rows[:max_seqs]
    total = max(len(rows), total_rows or 0)
    cols = min(ncols, max_cols, max(50, MAX_CELLS // max(1, len(shown_rows))))
    if total > len(shown_rows) or ncols > cols:
        notes.append(f"Large alignment: showing the first {len(shown_rows)} of {total} sequences and the first "
                     f"{cols} of {ncols} columns, so your browser stays fast. Problems further on are still listed "
                     f"above the table.")
    out = [f"<style>{CSS}</style>", '<div class="coalv">', f"<h2>{html.escape(title)}</h2>"]
    for n in notes:
        out.append(f'<p class="note">{html.escape(n)}</p>')
    out.append(_issues(rows, report, cols, {r["id"] for r in rows}))
    out.append('<div class="legend" aria-label="Colour key">')
    for cls, label in LEGEND:
        out.append(f'<span><span class="sw" style="background:var(--{cls})"></span>{label}</span>')
    out.append('<span><b>≠</b> codon does not match its amino acid (outlined)</span>'
               '<span><b>↯</b> frameshift (double outline, hatched)</span>'
               '<span style="color:var(--gap)">--- gap</span></div>')
    flagged_cols = {i + 1 for r in shown_rows for i, (_, _, f) in enumerate(r["cells"][:cols])
                    if f in ("mismatch", "frameshift")}
    out.append('<div class="scroll" tabindex="0" aria-label="Codon alignment, scrolls sideways">'
               '<table><thead><tr><th class="name" scope="col">Sequence</th>')
    for c in range(1, cols + 1):
        if c in flagged_cols:
            out.append(f'<th class="flag" id="coalv-c{c}" title="Column {c}: a problem in at least one sequence">'
                       f'{c}▾</th>')
        elif c == 1 or c % 10 == 0:
            out.append(f"<th>{c}</th>")
        else:
            out.append("<th></th>")
    out.append("</tr></thead><tbody>")
    for r in shown_rows:
        icon, word = _ICON.get(r["status"], ("?", ""))
        note = f'<br><span class="rownote">{html.escape(r["note"])}</span>' if r["note"] else ""
        out.append(f'<tr id="{row_id(r["id"])}"><th class="name" scope="row" title="{html.escape(r["id"])}">'
                   f'<span class="st" aria-label="{word}">{icon}</span>{html.escape(r["id"])}{note}</th>')
        for i, (aa, c, flag) in enumerate(r["cells"][:cols]):
            col = i + 1
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
            out.append(f'<td class="{" ".join(x for x in cls if x)}"{attr}><b>{html.escape(aa)}</b>'
                       f'<i>{html.escape(c) if c.strip() else "&nbsp;"}</i></td>')
        out.append("</tr>")
    out.append("</tbody></table></div></div>")
    return "".join(out)


def render_page(fragment: str, title: str = "COAL alignment view") -> str:
    """A complete, self-contained HTML document (no external scripts or fonts)."""
    return ("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"
            f"<title>{html.escape(title)}</title><style>body{{margin:0;padding:12px 16px;background:#fff}}</style>"
            f"</head><body>{fragment}<script>{JS}</script></body></html>")


def view(protein_path: str, codon_path: str | None = None, report: dict | str | None = None,
         dna_path: str | None = None, codontable: int | None = None, page: bool = True,
         max_seqs: int = MAX_SEQS, max_cols: int = MAX_COLS) -> str:
    """Render a COAL run as HTML. `report` may be a dict, a path to a coal.report/1 file, or None."""
    if isinstance(report, str):
        with open(report, encoding="utf-8") as fh:
            report = json.load(fh)
    protein = read_alignment(protein_path)
    codon = read_alignment(codon_path) if codon_path and os.path.isfile(codon_path) else None
    dna = read_alignment(dna_path) if dna_path and os.path.isfile(dna_path) else None
    # CHANGED(R2-A3): was rows built for every sequence (22,000 x 600 codons: 29 s, GBs) though only max_seqs
    # are drawn; now only the drawn sequences are built (problems of the others come from the report).
    total = len(protein)
    if total > max_seqs:
        shown = protein[:max_seqs]
        want = {p for p, _ in shown}
        if codon is not None:
            codon = [(i, r) for i, r in codon if i in want]
        if dna is not None:
            dids = {i for i, _ in dna}
            if all(p in dids for p in want) or len(dna) != total:
                dna = [(i, r) for i, r in dna if i in want]
            else:
                dna = dna[:max_seqs]  # paired by order: keep the same positions
        protein = shown
    rows, ncols, notes = build_rows(protein, codon, report, dna, codontable)
    if codon is None:
        notes.insert(0, "COAL wrote no codon alignment for this run. Codons below were placed by the viewer, "
                        "only to show where each gene breaks." if dna else
                     "COAL wrote no codon alignment for this run, and no DNA file was given.")
    title = f"Codon alignment: {os.path.basename(protein_path)}"
    frag = render_fragment(rows, ncols, notes, report, title, max_seqs, max_cols, total_rows=total)
    return render_page(frag, title) if page else frag


def render_run(protein_path: str, output_path: str | None, report: dict, dna_path: str | None, codontable: int,
               out_dir: str) -> tuple[str, str]:
    """CHANGED(R2-A3): new; was two full view() calls in ui.run_ui. Renders once: (fragment for the GUI, path of
    the downloadable page written into out_dir)."""
    frag = view(protein_path, output_path, report, dna_path, codontable, page=False)
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
        dna = read_alignment(dna_path)
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
