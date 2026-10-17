"""Assemble the files a Hugging Face Space needs into a folder (e.g. a clone of the Space repo).

    python space/build_space.py PATH/TO/SPACE_CLONE

Copies: app.py, README.md, requirements.txt (from space/), src/coal, src/coal_app and examples/.
Does not deploy anything; see space/DEPLOY.md.
"""
from __future__ import annotations

import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IGNORE = shutil.ignore_patterns("__pycache__", "*.pyc", "*.egg-info")


def build(dest: str) -> list[str]:
    os.makedirs(dest, exist_ok=True)
    for f in ("app.py", "README.md", "requirements.txt"):
        shutil.copy2(os.path.join(ROOT, "space", f), os.path.join(dest, f))
    for d in ("src/coal", "src/coal_app", "examples"):
        target = os.path.join(dest, d)
        if os.path.exists(target):
            shutil.rmtree(target)
        shutil.copytree(os.path.join(ROOT, d), target, ignore=IGNORE)
    shutil.copy2(os.path.join(ROOT, "LICENSE"), os.path.join(dest, "LICENSE"))
    return sorted(os.listdir(dest))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    print("wrote:", ", ".join(build(sys.argv[1])))
