"""Draw the COAL mark: desktop/coal.ico, and with --site the website's mark and favicons. Needs Pillow.

    python desktop/make_icon.py              # writes desktop/coal.ico
    python desktop/make_icon.py --site site  # also writes site/assets/coal-mark.svg, site/favicon.svg and PNGs

CHANGED(WEB-8): was a blue tile with a 3x3 block grid (Okabe-Ito colours); now the COAL logo mark. See
overnight/integration/CHANGES.md

The mark: 51 rounded squares in groups of three (codons), cooling from charcoal at the top to ember at the bottom.
Geometry and colours were measured on the owner's logo PNG (coal-mark-original.png, 1200 x 952): squares 88 px,
corner radius 12.4 px, 104 px pitch inside a group of three and 128 px from one group to the next, rows on a
112 px pitch; the first, second and last rows are shorter. A render of MARK_SVG differs from the PNG by a mean
of 0.13/255 per pixel over the mark.

At 32 px and below the 9 x 7 squares blur into a smudge, so the small sizes use a simplified mark: one codon
block of 3 x 3 squares with the same top-to-bottom cooling (charcoal, slate, ember, amber).
"""
import os
import sys

from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))

COLOURS = {"c": "#2A2A2A", "s": "#4A4744", "e": "#FF6A13", "a": "#FFC24B"}
SQUARE, RADIUS = 88, 12.4
NINE = (0, 104, 208, 336, 440, 544, 672, 776, 880)          # x of each square in a full row
SIX = (168, 272, 376, 504, 608, 712)                         # rows 2 and 7 sit half a group in
ROWS = [  # (x positions, colours), rows 112 px apart
    ((336, 440, 544), "scc"),
    (SIX, "ccscsc"),
    (NINE, "cscccc" "scc"),
    (NINE, "ccccsc" "sce"),
    (NINE, "csceca" "cee"),
    (NINE, "ceceea" "eae"),
    (SIX, "eeaaea"),
]
WIDTH, HEIGHT = 880 + SQUARE, 6 * 112 + SQUARE              # 968 x 760
# Small sizes: one codon block, 3 x 3.
SMALL = ["scc", "cec", "eae"]


def squares():
    for row, (xs, cols) in enumerate(ROWS):
        for x, k in zip(xs, cols):
            yield x, row * 112, k


def mark_svg() -> str:
    groups = {k: [] for k in "csea"}
    for x, y, k in squares():
        groups[k].append(f'<rect x="{x}" y="{y}" width="{SQUARE}" height="{SQUARE}" rx="{RADIUS}"/>')
    body = "".join(f'<g fill="{COLOURS[k]}">{"".join(v)}</g>' for k, v in groups.items())
    return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {WIDTH} {HEIGHT}">{body}</svg>\n'


def favicon_svg() -> str:
    """The 3 x 3 block; in a dark browser theme the charcoal squares get lighter so they stay visible."""
    rects = []
    for r, cols in enumerate(SMALL):
        for c, k in enumerate(cols):
            rects.append(f'<rect class="{k}" x="{c * 11}" y="{r * 11}" width="10" height="10" rx="1.6"/>')
    style = (".c{fill:#2A2A2A}.s{fill:#4A4744}.e{fill:#FF6A13}.a{fill:#FFC24B}"
             "@media (prefers-color-scheme:dark){.c{fill:#8C8782}.s{fill:#B5AFA8}}")
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="-1 -1 34 34"><style>{style}</style>'
            + "".join(rects) + "</svg>\n")


def _rgba(hex_colour):
    h = hex_colour.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4)) + (255,)


def draw(size: int, tile=None) -> Image.Image:
    """The mark at size x size px, drawn 8x larger and downsampled. tile: background colour or None."""
    ss = 8 * size
    im = Image.new("RGBA", (ss, ss), (0, 0, 0, 0))
    d = ImageDraw.Draw(im)
    if tile:
        d.rounded_rectangle([0, 0, ss - 1, ss - 1], radius=ss * 0.18, fill=_rgba(tile))
    if size <= 32:
        pad = ss * (0.16 if tile else 0.03)
        step = (ss - 2 * pad) / 3
        gap = step * 0.12
        for r, cols in enumerate(SMALL):
            for c, k in enumerate(cols):
                x0, y0 = pad + c * step + gap / 2, pad + r * step + gap / 2
                d.rounded_rectangle([x0, y0, x0 + step - gap, y0 + step - gap], radius=step * 0.14,
                                    fill=_rgba(COLOURS[k]))
    else:
        pad = ss * (0.12 if tile else 0.02)
        scale = (ss - 2 * pad) / WIDTH
        oy = (ss - HEIGHT * scale) / 2
        for x, y, k in squares():
            x0, y0 = pad + x * scale, oy + y * scale
            d.rounded_rectangle([x0, y0, x0 + SQUARE * scale, y0 + SQUARE * scale], radius=RADIUS * scale,
                                fill=_rgba(COLOURS[k]))
    return im.resize((size, size), Image.LANCZOS)


# The program icon sits on a light tile so the charcoal squares stay visible on a dark taskbar.
ICO_TILE = "#F4F2EF"


def main(argv):
    sizes = [16, 24, 32, 48, 64, 128, 256]
    big = draw(256, ICO_TILE)
    big.save(os.path.join(HERE, "coal.ico"), sizes=[(n, n) for n in sizes],
             append_images=[draw(n, ICO_TILE) for n in sizes[:-1]])
    print("wrote", os.path.join(HERE, "coal.ico"))
    if "--site" in argv:
        site = argv[argv.index("--site") + 1]
        os.makedirs(os.path.join(site, "assets"), exist_ok=True)
        with open(os.path.join(site, "assets", "coal-mark.svg"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(mark_svg())
        with open(os.path.join(site, "favicon.svg"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(favicon_svg())
        draw(32).save(os.path.join(site, "assets", "favicon-32.png"), optimize=True)
        draw(512).save(os.path.join(site, "assets", "icon-512.png"), optimize=True)
        # iOS fills transparency with black: give the touch icon an opaque light background.
        touch = Image.new("RGBA", (180, 180), _rgba("#FFFFFF"))
        touch.alpha_composite(draw(148), (16, 16))
        touch.convert("RGB").save(os.path.join(site, "assets", "apple-touch-icon.png"), optimize=True)
        print("wrote site icons in", site)


if __name__ == "__main__":
    main(sys.argv[1:])
