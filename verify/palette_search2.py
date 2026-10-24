"""Structured palette search: an even lightness ladder (greyscale-safe by construction), Okabe-Ito hue families per
class, every permutation of ladder steps to classes, chroma as high as the gamut and the text contrast allow; keep
the assignment whose worst simulated (protan / deutan / tritan) class distance is largest."""
import itertools
import sys

sys.path.insert(0, r"D:\coal-stress\verify")
from palette_search import HUES, contrast, de, hex2rgb, lab, lab2rgb, lab_of_lin, rgb2hex, scores, simulate, SIM, unlin  # noqa: E402


def in_gamut(L, a, b):
    fy = (L + 16) / 116; fx = fy + a / 500; fz = fy - b / 200
    g = lambda t: t ** 3 if t ** 3 > 0.008856 else (t - 16 / 116) / 7.787
    x, y, z = g(fx) * 0.95047, g(fy), g(fz) * 1.08883
    r = 3.2406 * x - 1.5372 * y - 0.4986 * z; gg = -0.9689 * x + 1.8758 * y + 0.0415 * z; bb = 0.0557 * x - 0.2040 * y + 1.0570 * z
    return all(-0.002 <= v <= 1.002 for v in (r, gg, bb))


def best_colour(name, L, text, min_contrast, cmax):
    a, b = HUES[name]
    for c in range(cmax, -1, -1):
        if name == "st":
            c = 0
        if in_gamut(L, a * c, b * c):
            rgb = lab2rgb(L, a * c, b * c)
            if contrast(rgb, text) >= min_contrast:
                return rgb
        if name == "st":
            return None
    return None


def run(label, ladder, text, cmax):
    names = list(HUES); best = None
    for perm in itertools.permutations(ladder):
        cols = {}
        for n, L in zip(names, perm):
            rgb = best_colour(n, L, hex2rgb(text), 4.5, cmax)
            if rgb is None:
                break
            cols[n] = rgb
        if len(cols) < 7:
            continue
        s = scores(cols); worst = min(s["protan"], s["deutan"], s["tritan"])
        if best is None or worst > best[0]:
            best = (worst, cols, s, perm)
    worst, cols, s, perm = best
    print(label, "worst simulated", round(worst, 1), {k: round(v, 1) for k, v in s.items()})
    for n in names:
        print(f'  "{n}": "{rgb2hex(cols[n])}",  L*={lab(cols[n])[0]:.0f} contrast={contrast(cols[n], hex2rgb(text)):.1f}')
    return cols


for cmax in (28, 34):
    light = run(f"LIGHT cmax={cmax}", [58, 64.3, 70.6, 77, 83.3, 89.6, 96], "#1a1a1a", cmax)
for cmax in (28, 40):
    dark = run(f"DARK cmax={cmax}", [14, 18.7, 23.3, 28, 32.7, 37.3, 42], "#ede9e4", cmax)
