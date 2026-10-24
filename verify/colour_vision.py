"""Human-checklist item: is the alignment view's residue-class palette readable for colour-blind users?
Renders the one-broken-gene example, reads the legend swatch colours and the cell background colours, simulates
protanopia / deuteranopia / tritanopia / achromatopsia with Chromium's own vision-deficiency emulation (screenshots)
and numerically (Machado et al. 2009 matrices, severity 1.0), and reports the smallest pairwise colour distance
between residue classes under each simulation. The view also carries the residue letter in every cell and marks
problems with outlines, so colour is never the only cue; this measures how much the colour adds."""
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, r"D:\coal-stress\verify")
import phase7
phase7.URL = "http://127.0.0.1:7861/"
from phase7 import new_page, run_example
from playwright.sync_api import sync_playwright

OUT = Path(r"D:\coal-stress\verify\p7\colour"); OUT.mkdir(parents=True, exist_ok=True)
# Machado, Oliveira & Fernandes (2009), severity 1.0, applied in linear RGB
SIM = {
    "protanopia": [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
    "deuteranopia": [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]],
    "tritanopia": [[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]],
    "achromatopsia": [[0.2126, 0.7152, 0.0722]] * 3,
}


def lin(c):
    c /= 255
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def to_lab(rgb_lin):
    r, g, b = rgb_lin
    x = (0.4124 * r + 0.3576 * g + 0.1805 * b) / 0.95047
    y = (0.2126 * r + 0.7152 * g + 0.0722 * b)
    z = (0.0193 * r + 0.1192 * g + 0.9505 * b) / 1.08883
    f = lambda t: t ** (1 / 3) if t > 0.008856 else 7.787 * t + 16 / 116
    fx, fy, fz = f(x), f(y), f(z)
    return 116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz)


def simulate(rgb, m):
    v = [lin(float(c)) for c in rgb]
    return [max(0.0, min(1.0, sum(m[i][j] * v[j] for j in range(3)))) for i in range(3)]


def de(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def parse(css):
    nums = [float(x) for x in css.replace("rgba(", "").replace("rgb(", "").rstrip(")").split(",")[:3]]
    return nums


with sync_playwright() as p:
    b = p.chromium.launch(); ctx = b.new_context(viewport={"width": 1366, "height": 900}); pg = new_page(ctx)
    run_example(pg, "One broken gene"); pg.wait_for_timeout(800)
    legend = pg.locator("#coal-view").evaluate("""e => Array.from(e.querySelectorAll('.coalv .legend span, .coalv .leg span, .coalv [class*=legend] span, .coalv [class*=key] span'))
        .map(s => [s.textContent.trim(), getComputedStyle(s).backgroundColor, getComputedStyle(s).borderColor]).filter(x => x[0])""")
    cells = pg.locator("#coal-view").evaluate("""e => { const m = {}; for (const c of e.querySelectorAll('.coalv td, .coalv .cell, .coalv [class*=aa]')) {
        const cls = c.className; const bg = getComputedStyle(c).backgroundColor; if (!m[cls]) m[cls] = bg; } return m; }""")
    view = pg.locator("#coal-view")
    client = ctx.new_cdp_session(pg)
    for kind in ("none", "protanopia", "deuteranopia", "tritanopia", "achromatopsia"):
        client.send("Emulation.setEmulatedVisionDeficiency", {"type": kind})
        pg.wait_for_timeout(300)
        view.screenshot(path=str(OUT / f"view_{kind}.png"))
    client.send("Emulation.setEmulatedVisionDeficiency", {"type": "none"})
    b.close()

(OUT / "legend.json").write_text(json.dumps({"legend": legend, "cells": cells}, indent=1), encoding="utf-8")
print("legend swatches:", [(t, bg) for t, bg, _ in legend][:12])
classes = {}
for t, bg, _ in legend:
    if bg.startswith("rgb") and bg not in ("rgba(0, 0, 0, 0)",):
        classes[t] = parse(bg)
if not classes:
    for cls, bg in cells.items():
        if bg.startswith("rgb") and bg != "rgba(0, 0, 0, 0)":
            classes[cls[:30]] = parse(bg)
names = list(classes)
print(f"{len(names)} coloured classes")
report = {}
for kind in ("none",) + tuple(SIM):
    labs = {n: to_lab(simulate(classes[n], SIM[kind]) if kind != "none" else [lin(float(c)) for c in classes[n]]) for n in names}
    pairs = sorted(((de(labs[a], labs[b]), a, b) for i, a in enumerate(names) for b in names[i + 1:]))
    report[kind] = [(round(d, 1), a, b) for d, a, b in pairs[:4]]
    print(f"{kind:14s} smallest class distances (CIE76 dE):", ", ".join(f"{d:.0f} {a}/{b}" for d, a, b in pairs[:3]))
(OUT / "distances.json").write_text(json.dumps(report, indent=1), encoding="utf-8")
