"""HTML pieces of the GUI (header, result banner), kept free of Gradio so they can be tested and rendered on
their own. CHANGED(GUI-1): new file. See overnight/integration/CHANGES.md

Colours come from CSS custom properties defined in ui.CSS (light) and re-defined under Gradio's ``.dark`` body
class (dark), so these fragments follow the app theme. Every status has an icon and a word, never colour alone.
"""
from __future__ import annotations

import html
import re

# The COAL mark (site/assets/coal-mark.svg), one class per tile colour so CSS can swap the two greys in dark mode
# (site/assets/coal-mark-dark.svg uses #8C8782 / #B5AFA8 for them). Decorative: aria-hidden.
_MARK_TILES = {
    "m1": [(440, 0), (544, 0), (168, 112), (272, 112), (504, 112), (712, 112), (0, 224), (208, 224), (336, 224),
           (440, 224), (544, 224), (776, 224), (880, 224), (0, 336), (104, 336), (208, 336), (336, 336),
           (544, 336), (776, 336), (0, 448), (208, 448), (440, 448), (672, 448), (0, 560), (208, 560)],
    "m2": [(336, 0), (376, 112), (608, 112), (104, 224), (672, 224), (440, 336), (672, 336), (104, 448)],
    "m3": [(880, 336), (336, 448), (776, 448), (880, 448), (104, 560), (336, 560), (440, 560), (672, 560),
           (880, 560), (168, 672), (272, 672), (608, 672)],
    "m4": [(544, 448), (544, 560), (776, 560), (376, 672), (504, 672), (712, 672)],
}


def mark_svg(size: int = 40) -> str:
    groups = "".join(
        f'<g class="{cls}">' + "".join(f'<rect x="{x}" y="{y}" width="88" height="88" rx="12.4"/>' for x, y in tiles)
        + "</g>" for cls, tiles in _MARK_TILES.items())
    h = round(size * 760 / 968)
    return (f'<svg class="coal-mark" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 968 760" width="{size}" '
            f'height="{h}" aria-hidden="true" focusable="false">{groups}</svg>')


def header_html(intro: str = "") -> str:
    """Page header: the COAL mark and the name.
    CHANGED(GUI-13): was the mark, "COAL codon alignment", a two-line description and the where-the-data-goes note
    (core.intro_text()). The owner: "just COAL is fine"; the description is gone and the note is in footer_html.
    `intro` is accepted and ignored so older callers keep working."""
    return f'<header class="coal-head">{mark_svg(32)}<h1>COAL</h1></header>'


def footer_html(intro: str) -> str:
    """CHANGED(GUI-13): new. The where-the-data-goes note (core.intro_text()), small, at the foot of the page (as
    Squoosh puts "images never leave your device" below its drop zone, not in its header).
    CHANGED(POL-1): "" when there is no note (the local app has none now; the hosted Space keeps its note)."""
    return f'<p class="coal-foot">{html.escape(intro)}</p>' if intro else ""


# ---------------------------------------------------------------- result banner

_KIND = {  # banner kind: (icon, css class)
    "ok": ("✔", "ok"),
    "warning": ("⚠", "warning"),
    "failed": ("✖", "failed"),
    "info": ("ℹ", "info"),
}
_TILE = {"ok": ("ok", ""), "warning": ("warning", "!"), "failed": ("failed", "✖"), "not_checked": ("unchecked", "–")}
_WORD = {"ok": "OK", "warning": "Warning", "failed": "Failed", "not_checked": "Not checked"}
STRIP_MAX = 400  # at most this many gene tiles in the banner


def _kind_of(head: str) -> str:
    """Banner kind from core.banner_head's text (its icon)."""
    head = head.lstrip("# ").strip()
    if head.startswith("✔"):
        return "ok"
    if head.startswith("✖"):
        return "failed"
    return "warning"


def _plural(n: int, one: str, many: str | None = None) -> str:
    return f"{n:,} {one if n == 1 else (many or one + 's')}"


def headline(summary: dict, kind: str) -> str:
    """'19 of 20 genes aligned' style headline from a report summary (aligned = OK or with warnings)."""
    total = int(summary.get("total") or 0)
    ok = int(summary.get("ok") or 0)
    aligned = ok + int(summary.get("warning") or 0)
    if kind == "failed":
        return "COAL stopped: no alignment was written"
    if total and ok == total:
        return "The gene matches codon for codon" if total == 1 else f"All {total:,} genes match codon for codon"
    return f"{aligned:,} of {_plural(total, 'gene')} aligned"


def _chips(summary: dict) -> str:
    out = []
    for key, word in (("failed", "failed"), ("warning", "with warnings"), ("not_checked", "not checked")):
        n = int(summary.get(key) or 0)
        if n:
            icon = {"failed": "✖", "warning": "⚠", "not_checked": "–"}[key]
            cls = _TILE[key][0]
            out.append(f'<span class="coal-chip coal-{cls}"><span aria-hidden="true">{icon}</span> '
                       f'{n:,} {word}</span>')
    return "".join(out)


def gene_strip(report: dict) -> str:
    """One small tile per gene, in file order (like the tiles of the COAL mark): quiet grey for OK, ember with '!'
    for a warning, crimson with '✖' for a failure (CHANGED(AMB-1): were amber and ember), dashed for not checked.
    Decorative for screen readers (the counts are in the text); each tile's tooltip names the gene."""
    seqs = report.get("sequences") or []
    if not seqs or len(seqs) > STRIP_MAX:
        return ""
    tiles = []
    for s in seqs:
        st = s.get("status")
        cls, glyph = _TILE.get(st, ("unchecked", "?"))
        tip = html.escape(f"{s.get('id', '')}: {_WORD.get(st, str(st))}", quote=True)
        tiles.append(f'<span class="coal-tile coal-{cls}" title="{tip}">{glyph}</span>')
    return f'<div class="coal-strip" aria-hidden="true">{"".join(tiles)}</div>'


def extra_sentence(sentence: str) -> str:
    """CHANGED(GUI-13): the part of core.summary_sentence that the headline and chips do not already say (the
    columns removed by -nomismatch); "" when it only repeats the counts ("20 sequences: 19 OK, ...")."""
    m = re.match(r"\d[\d,]* sequences?: [^.]*\.\s*", sentence or "")
    return (sentence or "")[m.end():].strip() if m else (sentence or "").strip()


def run_banner(head: str, report: dict, sentence: str, files: str = "") -> str:
    """Summary banner after a run. head: core.banner_head(...) (markdown '### ✔ Done ...'); sentence:
    core.summary_sentence(report); files: 'Files used: ...'.
    CHANGED(GUI-13): was a status line ("⚠ Done, with some problems (see the table)") above the headline and the
    summary sentence ("24 sequences: 0 OK, 24 with warnings, 0 failed.") under the tiles, both repeating the
    headline and the chips. Now: icon + headline, chips, tiles, then only what is new (columns removed) and the
    files line."""
    kind = _kind_of(head)
    icon, cls = _KIND[kind]
    s = report.get("summary") or {}
    title = headline(s, kind)
    extra = extra_sentence(sentence)
    return (
        f'<div class="coal-banner coal-{cls}" role="status">'
        f'<p class="coal-headline"><span class="coal-icon" aria-hidden="true">{icon}</span> {html.escape(title)}</p>'
        f'<div class="coal-chips">{_chips(s)}</div>'
        f'{gene_strip(report)}'
        + (f'<p class="coal-detail">{html.escape(extra)}</p>' if extra else "")
        + (f'<p class="coal-files">{html.escape(files)}</p>' if files else "")
        + '</div>')


def _code_name(code, names: dict) -> str:
    name = (names or {}).get(code, "")
    return f"code {code} ({name.split(' (')[0]})" if name else f"code {code}"


def code_fix_banner(fix: dict, report: dict, names: dict | None = None, files: str = "") -> str:
    """CHANGED(AMB-4): new. The banner when the genetic-code hint explains the run (core.code_fix): it reads as
    "the code chosen does not fit these genes; here is the fix", not as a page of warnings. Informational style
    (nothing is broken, the choice of code is), icon + words; no gene tiles, chips only for failed or unchecked
    genes (those are real problems). The re-run button (ui: #coal-hint-rerun, primary) sits under it in the
    same box (ui.CSS #coal-fixbox). names: code -> short name (ui.CODE_SHORT)."""
    s = report.get("summary") or {}
    cur, new = fix.get("current"), fix["code"]
    checked, fit, bad = fix["checked"], fix["fit"], fix["checked"] - fix["fit_current"]
    used = _code_name(cur, names) if cur else "the code used"
    good = f"all {checked:,}" if fit == checked else f"{fit:,} of {checked:,}"
    text = (f"With {used} {bad:,} of {_plural(checked, 'gene')} have mismatched codons; with "
            f"{_code_name(new, names)} {good} match.")
    real = {k: s.get(k) for k in ("failed", "not_checked") if s.get(k)}
    return (
        '<div class="coal-banner coal-info coal-codefix" role="status">'
        '<p class="coal-headline"><span class="coal-icon" aria-hidden="true">ℹ</span> Wrong genetic code?</p>'
        f'<p class="coal-detail">{html.escape(text)}</p>'
        + (f'<div class="coal-chips">{_chips(real)}</div>' if real else "")
        + (f'<p class="coal-small">{html.escape(fix["others"])}</p>' if fix.get("others") else "")
        + (f'<p class="coal-files">{html.escape(files)}</p>' if files else "")
        + '</div>')


def message_banner(kind: str, title: str, text: str = "", text_html: str = "") -> str:
    """A banner without a report: kind 'warning' / 'failed' / 'info'. text is escaped; text_html is not."""
    icon, cls = _KIND.get(kind, _KIND["info"])
    body = text_html or html.escape(text)
    return (f'<div class="coal-banner coal-{cls}" role="status">'
            f'<p class="coal-headline"><span class="coal-icon" aria-hidden="true">{icon}</span> '
            f'{html.escape(title)}</p>' + (f'<p class="coal-detail">{body}</p>' if body else "") + '</div>')


def empty_banner() -> str:
    # CHANGED(GUI-13): was a dashed box "Results appear here: add the two files (or pick an example) and press
    # Make codon alignment." directly under that button; then one quiet line "Results appear here."
    # CHANGED(POL-4): nothing. The whole results area (heading, banner, view, methods) stays hidden until there is
    # something to show (ui.build: #coal-results visible=False, made visible by every result), as Squoosh shows no
    # result panel before an image is opened. The form above already names the next step (the button).
    return ""


def sample_banner(report: dict, sentence: str) -> str:
    """Banner of the demo sample report (core.sample_banner's text, as HTML).
    CHANGED(GUI-13): was "Demo sample, not your data" + a two-sentence explanation + the counts sentence; now the
    status line, the headline, chips, tiles and one short line (the tip under the banner says what to do)."""
    s = report.get("summary") or {}
    kind = "warning" if (s.get("failed") or s.get("warning")) else "ok"
    return (
        '<div class="coal-banner coal-info" role="status">'
        '<p class="coal-status"><span class="coal-icon" aria-hidden="true">ℹ</span> Sample report, not your data</p>'
        f'<p class="coal-headline">{html.escape(headline(s, kind))}</p>'
        f'<div class="coal-chips">{_chips(s)}</div>'
        f'{gene_strip(report)}'
        '<p class="coal-detail">The <em>Mitochondrial (wrong code)</em> example: run with genetic code 1 on purpose.'
        + (f' {html.escape(extra_sentence(sentence))}' if extra_sentence(sentence) else "")
        + '</p></div>')
