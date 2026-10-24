"""Search a residue-class palette (7 classes: hy ar po ne pl sp st) that keeps the Okabe-Ito hue families but varies
lightness, so the classes stay apart under protanopia / deuteranopia / tritanopia and in greyscale, with the viewer's
text colour at >= 4.5:1 on every tint. Light and dark sets. Prints the best found and its scores."""
import itertools
import math
import random

SIM = {
    "protan": [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
    "deutan": [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]],
    "tritan": [[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]],
}


def hex2rgb(h):
    h = h.lstrip("#"); return [int(h[i:i + 2], 16) for i in (0, 2, 4)]


def rgb2hex(c):
    return "#%02x%02x%02x" % tuple(max(0, min(255, round(x))) for x in c)


def lin(c):
    c /= 255; return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def unlin(v):
    v = max(0.0, min(1.0, v)); return 255 * (12.92 * v if v <= 0.0031308 else 1.055 * v ** (1 / 2.4) - 0.055)


def lum(rgb):
    r, g, b = (lin(float(c)) for c in rgb); return 0.2126 * r + 0.7152 * g + 0.0722 * b


def contrast(a, b):
    la, lb = lum(a), lum(b); hi, lo = max(la, lb), min(la, lb); return (hi + 0.05) / (lo + 0.05)


def lab_of_lin(r, g, b):
    x = (0.4124 * r + 0.3576 * g + 0.1805 * b) / 0.95047; y = 0.2126 * r + 0.7152 * g + 0.0722 * b
    z = (0.0193 * r + 0.1192 * g + 0.9505 * b) / 1.08883
    f = lambda t: t ** (1 / 3) if t > 0.008856 else 7.787 * t + 16 / 116
    fx, fy, fz = f(x), f(y), f(z); return 116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz)


def lab(rgb):
    return lab_of_lin(*(lin(float(c)) for c in rgb))


def lab2rgb(L, a, b):
    fy = (L + 16) / 116; fx = fy + a / 500; fz = fy - b / 200
    g = lambda t: t ** 3 if t ** 3 > 0.008856 else (t - 16 / 116) / 7.787
    x, y, z = g(fx) * 0.95047, g(fy), g(fz) * 1.08883
    r = 3.2406 * x - 1.5372 * y - 0.4986 * z; gg = -0.9689 * x + 1.8758 * y + 0.0415 * z; bb = 0.0557 * x - 0.2040 * y + 1.0570 * z
    return [unlin(r), unlin(gg), unlin(bb)]


def simulate(rgb, m):
    v = [lin(float(c)) for c in rgb]
    return [max(0.0, min(1.0, sum(m[i][j] * v[j] for j in range(3)))) for i in range(3)]


def de(p, q):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(p, q)))


def scores(cols):
    names = list(cols); out = {}
    for kind in ("none", "protan", "deutan", "tritan", "grey"):
        labs = {}
        for n in names:
            if kind == "none":
                labs[n] = lab(cols[n])
            elif kind == "grey":
                L = lab(cols[n])[0]; labs[n] = (L, 0, 0)
            else:
                labs[n] = lab_of_lin(*simulate(cols[n], SIM[kind]))
        out[kind] = min(de(labs[a], labs[b]) for a, b in itertools.combinations(names, 2))
    return out


def search(hues, L_range, text, min_contrast, seed=1, iters=40000):
    """hues: name -> (a, b) direction in Lab (unit), L_range: (lo, hi) for the lightness ladder. Random search over
    per-class L and chroma, keeping text contrast >= min_contrast; maximise the worst simulated distance."""
    rng = random.Random(seed); names = list(hues); best = None
    for _ in range(iters):
        Ls = sorted(rng.uniform(*L_range) for _ in names); rng.shuffle(Ls)
        cols = {}
        ok = True
        for n, L in zip(names, Ls):
            c = rng.uniform(8, 40) if n != "st" else rng.uniform(0, 4)
            a, b = hues[n]; rgb = lab2rgb(L, a * c, b * c)
            if contrast(rgb, text) < min_contrast:
                ok = False; break
            cols[n] = rgb
        if not ok:
            continue
        s = scores(cols); worst = min(s["protan"], s["deutan"], s["tritan"], s["grey"] * 1.5, s["none"])
        if best is None or worst > best[0]:
            best = (worst, cols, s)
    return best


# hue directions in Lab (unit vectors), Okabe-Ito families: hy sky blue, ar reddish purple, po blue, ne vermilion/orange,
# pl bluish green, sp yellow, st neutral grey
HUES = {"hy": (-0.35, -0.94), "ar": (0.8, -0.6), "po": (0.3, -0.95), "ne": (0.75, 0.66), "pl": (-0.9, 0.3), "sp": (0.0, 1.0), "st": (0.0, 0.0)}
light = search(HUES, (62, 95), hex2rgb("#1a1a1a"), 4.5, seed=7, iters=60000)
dark = search(HUES, (18, 48), hex2rgb("#ede9e4"), 4.5, seed=11, iters=60000)
for label, best in (("LIGHT", light), ("DARK", dark)):
    worst, cols, s = best
    print(label, "worst", round(worst, 1), {k: round(v, 1) for k, v in s.items()})
    for n in HUES:
        print(f'  "{n}": "{rgb2hex(cols[n])}",  L*={lab(cols[n])[0]:.0f} contrast={contrast(cols[n], hex2rgb("#1a1a1a" if label == "LIGHT" else "#ede9e4")):.1f}')
cur = {"hy": "#d4ebf9", "ar": "#f0d6e6", "po": "#c9dcef", "ne": "#f8d5c2", "pl": "#c9eadf", "sp": "#faf4c0", "st": "#d9d9d9"}
print("CURRENT light", {k: round(v, 1) for k, v in scores({k: hex2rgb(v) for k, v in cur.items()}).items()})
curd = {"hy": "#1e3b52", "ar": "#4a2a40", "po": "#2f3157", "ne": "#573420", "pl": "#1d4537", "sp": "#48421c", "st": "#3d3b39"}
print("CURRENT dark", {k: round(v, 1) for k, v in scores({k: hex2rgb(v) for k, v in curd.items()}).items()})
