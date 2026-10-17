"""HTML pieces of the GUI (header, result banner), kept free of Gradio so they can be tested and rendered on
their own. CHANGED(GUI-1): new file. See overnight/integration/CHANGES.md

Colours come from CSS custom properties defined in ui.CSS (light) and re-defined under Gradio's ``.dark`` body
class (dark), so these fragments follow the app theme. Every status has an icon and a word, never colour alone.
"""
from __future__ import annotations

import html

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


def header_html(intro: str) -> str:
    """Short page header: mark, name, one sentence, and where the data goes (core.intro_text())."""
    return (
        '<header class="coal-head">'
        f'{mark_svg()}'
        '<div><h1>COAL <span class="coal-head-sub">codon alignment</span></h1>'
        '<p>Give COAL a protein alignment and the DNA of the same genes: it lines the DNA up '
        '<abbr title="Three DNA letters that code for one amino acid">codon</abbr> by codon and says, in plain '
        'words, which genes don\'t fit and why.</p>'
        f'<p class="coal-head-note">{html.escape(intro)}</p></div>'
        '</header>')


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
    """One small tile per gene, in file order (like the tiles of the COAL mark): quiet grey for OK, amber with '!'
    for a warning, ember with '✖' for a failure, dashed for not checked. Decorative for screen readers (the counts
    are in the text); each tile's tooltip names the gene."""
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


def run_banner(head: str, report: dict, sentence: str, files: str = "") -> str:
    """Summary banner after a run. head: core.banner_head(...) (markdown '### ✔ Done ...'); sentence:
    core.summary_sentence(report); files: 'Files used: ...'."""
    kind = _kind_of(head)
    icon, cls = _KIND[kind]
    status = head.lstrip("# ").lstrip("✔⚠✖ ").strip()
    s = report.get("summary") or {}
    title = headline(s, kind)
    line = (f'<p class="coal-status"><span class="coal-icon" aria-hidden="true">{icon}</span> '
            f'{html.escape(status)}</p>') if status != title else ""
    lead = "" if line else f'<span class="coal-icon" aria-hidden="true">{icon}</span> '
    return (
        f'<div class="coal-banner coal-{cls}" role="status">{line}'
        f'<p class="coal-headline">{lead}{html.escape(title)}</p>'
        f'<div class="coal-chips">{_chips(s)}</div>'
        f'{gene_strip(report)}'
        f'<p class="coal-detail">{html.escape(sentence)}</p>'
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
    return ('<div class="coal-banner coal-empty" role="status"><p class="coal-detail">Results appear here: add the '
            'two files (or pick an example) and press <strong>Make codon alignment</strong>.</p></div>')


def sample_banner(report: dict, sentence: str) -> str:
    """Banner of the demo sample report (core.sample_banner's text, as HTML)."""
    s = report.get("summary") or {}
    kind = "warning" if (s.get("failed") or s.get("warning")) else "ok"
    return (
        '<div class="coal-banner coal-info" role="status">'
        '<p class="coal-status"><span class="coal-icon" aria-hidden="true">ℹ</span> Demo sample, not your data</p>'
        f'<p class="coal-headline">{html.escape(headline(s, kind))}</p>'
        f'<div class="coal-chips">{_chips(s)}</div>'
        f'{gene_strip(report)}'
        '<p class="coal-detail">The 24 mitochondrial genes of the <em>Mitochondrial genes</em> example, run with '
        'genetic code 1, the wrong code for them: see the tip below to re-run with code 2. '
        f'{html.escape(sentence)}</p></div>')
