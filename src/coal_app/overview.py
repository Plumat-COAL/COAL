"""Overview of a whole alignment, for any size: two thin strips drawn as two small SVGs from aggregated numbers.

CHANGED(GUI-10): new file. See overnight/integration/CHANGES.md and docs/ui-principles.md (Shneiderman's "overview
first, zoom and filter, then details on demand"; Jalview's overview window with its box for the part on screen).

- Sequences strip: the genes in file order, put into at most BINS bins; each bin shows the worst status in it
  (failed: full height, warning: three quarters, OK: a low grey bar), and ticks under it mark the bins that hold a
  sequence drawn in the preview.
- Columns strip: per column bin, the share of gaps (grey bars rising from the bottom, estimated from an even sample of
  at most GAP_SAMPLE sequences) and the number of sequences with a mismatch or frameshift there (ember bars hanging
  from the top, all sequences, from the report); an outlined box marks the columns drawn in the preview.

Never one element per sequence or column: each series is one <path>, so 100,000 sequences cost the same as 10.
Colour is never the only cue: heights and directions differ and the counts are in the text next to it.
"""
from __future__ import annotations

import html

BINS = 400          # most bins per strip
GAP_SAMPLE = 2000   # sequences read for the gap profile
_GAP = bytes(1 if chr(i) in "-." else 0 for i in range(256))  # translate table: 1 for a gap, 0 for anything else
_RANK = {"failed": 3, "warning": 2, "not_checked": 1, "ok": 0}


def _even(n: int, k: int) -> list[int]:
    """k indices spread evenly over range(n), first and last included (k >= n: all)."""
    if k >= n:
        return list(range(n))
    if k <= 1:
        return [0] if n else []
    return sorted({round(i * (n - 1) / (k - 1)) for i in range(k)})


def gap_profile(rows: list[str], ncols: int, sample: int = GAP_SAMPLE) -> tuple[list[float], int]:
    """(share of gaps per column, sequences sampled). Sums columns 255 rows at a time as bytes of one big integer
    (one byte per column), which is ~100 times faster than a loop over characters."""
    if not rows or ncols <= 0:
        return [0.0] * max(ncols, 0), 0
    idx = _even(len(rows), sample)
    counts = [0] * ncols
    acc, k = 0, 0

    def flush(a):
        if a:
            for i, v in enumerate(a.to_bytes(ncols + 1, "little")[:ncols]):
                if v:
                    counts[i] += v

    for i in idx:
        r = rows[i]
        b = r.encode("ascii", "replace")[:ncols].translate(_GAP)
        if len(b) < ncols:  # a short row: the missing end counts as gaps
            b += b"\x01" * (ncols - len(b))
        acc += int.from_bytes(b, "little")
        k += 1
        if k == 255:
            flush(acc)
            acc, k = 0, 0
    flush(acc)
    n = len(idx)
    return [c / n for c in counts], n


def problem_columns(report: dict | None, ncols: int) -> list[int]:
    """Per alignment column: how many sequences have a mismatch or frameshift (or the break of a failed gene) there."""
    counts = [0] * max(ncols, 0)
    for s in (report or {}).get("sequences") or []:
        seen = set()
        for m in s.get("messages") or []:
            c = m.get("column")
            if isinstance(c, int) and 1 <= c <= ncols and c not in seen:
                seen.add(c)
                counts[c - 1] += 1
    return counts


def _bins(n: int, nb: int) -> list[tuple[int, int]]:
    """nb (start, end) half-open ranges covering range(n)."""
    return [(i * n // nb, (i + 1) * n // nb) for i in range(nb)]


def _fmt(x: float) -> str:
    return f"{x:.2f}".rstrip("0").rstrip(".")


def _svg(w: int, h: int, label: str, body: str, cls: str) -> str:
    return (f'<svg class="ov {cls}" viewBox="0 0 {w} {h}" preserveAspectRatio="none" role="img" '
            f'aria-label="{html.escape(label, quote=True)}" focusable="false">{body}</svg>')


def sequences_strip(statuses: list[str | None], drawn: set[int]) -> str:
    """statuses: per sequence in file order; drawn: indexes of the sequences drawn in the preview."""
    n = len(statuses)
    if not n:
        return ""
    nb = min(n, BINS)
    H, T = 20, 5  # bar area, tick area under it
    paths: dict[str, list[str]] = {"failed": [], "warning": [], "not_checked": [], "ok": []}
    ticks = []
    drawn_bins = {i * nb // n for i in drawn}
    for b, (lo, hi) in enumerate(_bins(n, nb)):
        worst = max((statuses[i] for i in range(lo, hi)), key=lambda s: _RANK.get(s, 1), default="ok")
        worst = worst if worst in paths else "not_checked"
        h = {"failed": H, "warning": H * 3 // 4, "not_checked": H // 2, "ok": H // 3}[worst]
        paths[worst].append(f"M{b} {H}v-{h}h1v{h}z")
        if b in drawn_bins:
            ticks.append(f"M{b} {H + 2}h1v{T - 2}h-1z")
    body = "".join(f'<path class="ov-{k}" d="{"".join(v)}"/>' for k, v in paths.items() if v)
    if ticks:
        body += f'<path class="ov-drawn" d="{"".join(ticks)}"/>'
    counts = {k: sum(1 for s in statuses if s == k) for k in ("failed", "warning")}
    label = (f"{n:,} sequences in file order: {counts['failed']:,} failed, {counts['warning']:,} with warnings; "
             f"ticks mark the {len(drawn):,} drawn in the preview")
    return _svg(nb, H + T, label, body, "ov-seq")


def columns_strip(gaps: list[float], problems: list[int], window: tuple[int, int]) -> str:
    """gaps/problems: per column; window: (first, last) column drawn in the preview, 1-based."""
    ncols = len(gaps)
    if not ncols:
        return ""
    nb = min(ncols, BINS)
    H = 36
    top = max(problems) if problems else 0
    gp, pp = [], []
    for b, (lo, hi) in enumerate(_bins(ncols, nb)):
        g = sum(gaps[lo:hi]) / max(1, hi - lo)
        if g > 0.005:
            h = max(1.0, g * (H * 0.6))
            gp.append(f"M{b} {H}v-{_fmt(h)}h1v{_fmt(h)}z")
        p = max(problems[lo:hi], default=0)
        if p:
            # a square-root scale, so one odd sequence is still visible next to a column where hundreds break
            h = max(3.0, (p / top) ** 0.5 * (H * 0.45))
            pp.append(f"M{b} 0v{_fmt(h)}h1v-{_fmt(h)}z")
    body = ""
    if gp:
        body += f'<path class="ov-gap" d="{"".join(gp)}"/>'
    if pp:
        body += f'<path class="ov-prob" d="{"".join(pp)}"/>'
    c0, c1 = window
    x0 = (c0 - 1) * nb / ncols
    x1 = c1 * nb / ncols
    body += (f'<rect class="ov-box" x="{_fmt(x0)}" y="0.5" width="{_fmt(max(x1 - x0, 0.5))}" height="{H - 1}" '
             f'vector-effect="non-scaling-stroke"/>')
    nprob = sum(1 for p in problems if p)
    label = (f"{ncols:,} columns: {nprob:,} with a problem in at least one sequence; the box marks columns "
             f"{c0:,} to {c1:,}, drawn in the preview")
    return _svg(nb, H, label, body, "ov-col")


def overview_html(ids: list[str], report: dict | None, protein_rows: list[str], ncols: int,
                  drawn_ids: set[str], window: tuple[int, int], columns: bool = True) -> str:
    """The overview block (HTML) for the viewer fragment. ids/protein_rows: every sequence, in file order.
    columns=False leaves out the columns strip (when -nogap/-nomismatch removed columns, the drawn columns are
    not protein-alignment columns)."""
    st = {s.get("id"): s.get("status") for s in (report or {}).get("sequences") or []}
    statuses = [st.get(i) for i in ids]
    drawn = {k for k, i in enumerate(ids) if i in drawn_ids}
    gaps, sampled = gap_profile(protein_rows, ncols) if columns else ([], 0)
    probs = problem_columns(report, ncols)
    n = len(ids)
    nfail = sum(1 for s in statuses if s == "failed")
    nwarn = sum(1 for s in statuses if s == "warning")
    nprob_cols = sum(1 for p in probs if p)
    gap_note = (f"gaps estimated from {sampled:,} sequences spread through the file" if sampled < n
                else "share of gaps")
    return (
        '<div class="overview" role="group" aria-label="Overview of the whole alignment">'
        '<p class="ov-title"><strong>Whole alignment</strong> '
        f'<span class="ov-sum">{n:,} sequences: <span class="ov-k ov-k-failed">✖ {nfail:,} failed</span>, '
        f'<span class="ov-k ov-k-warning">⚠ {nwarn:,} with warnings</span>; {ncols:,} columns, '
        f'{nprob_cols:,} with a problem</span></p>'
        f'<div class="ov-row"><span class="ov-lab">Sequences</span>{sequences_strip(statuses, drawn)}'
        f'<span class="ov-ends"><span>1</span><span>{n:,}</span></span></div>'
        + (f'<div class="ov-row"><span class="ov-lab">Columns</span>{columns_strip(gaps, probs, window)}'
           f'<span class="ov-ends"><span>1</span><span>{ncols:,}</span></span></div>' if columns else "") +
        '<p class="ov-key"><span><i class="k k-fail"></i>tall mark: a failed gene</span>'
        '<span><i class="k k-warn"></i>shorter mark: a warning</span>'
        '<span><i class="k k-drawn"></i>tick: drawn below</span>'
        + ('<span><i class="k k-prob"></i>hanging bar: sequences with a problem in that column</span>'
           f'<span><i class="k k-gap"></i>rising bar: {html.escape(gap_note)}</span>'
           '<span><i class="k k-box"></i>box: the columns drawn below</span>' if columns else "")
        + '</p></div>')
