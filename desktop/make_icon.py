"""Draw desktop/coal.ico (the program and shortcut icon). Needs Pillow; run once, commit the .ico.

    python desktop/make_icon.py

Design: a blue (Okabe-Ito #0072B2) rounded square with three columns of three blocks,
a codon alignment in miniature; the middle block of the right column is orange (#E69F00),
the one codon that doesn't fit. Readable at 16 px because the blocks merge into stripes.
"""
import os

from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))
BLUE, ORANGE, WHITE = (0, 114, 178, 255), (230, 159, 0, 255), (255, 255, 255, 255)


def draw(size: int) -> Image.Image:
    s = 4 * size  # draw large, then downsample for smooth edges
    im = Image.new("RGBA", (s, s), (0, 0, 0, 0))
    d = ImageDraw.Draw(im)
    d.rounded_rectangle([0, 0, s - 1, s - 1], radius=s // 5, fill=BLUE)
    pad, gap = s * 0.17, s * 0.05
    cell = (s - 2 * pad - 2 * gap) / 3
    for col in range(3):
        for row in range(3):
            x0 = pad + col * (cell + gap)
            y0 = pad + row * (cell + gap)
            colour = ORANGE if (col, row) == (2, 1) else WHITE
            d.rounded_rectangle([x0, y0, x0 + cell, y0 + cell * 0.82], radius=max(1, s // 40), fill=colour)
    return im.resize((size, size), Image.LANCZOS)


def main():
    sizes = [16, 24, 32, 48, 64, 128, 256]
    big = draw(256)
    big.save(os.path.join(HERE, "coal.ico"), sizes=[(n, n) for n in sizes],
             append_images=[draw(n) for n in sizes[:-1]])
    print("wrote", os.path.join(HERE, "coal.ico"))


if __name__ == "__main__":
    main()
