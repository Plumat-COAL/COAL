"""Assemble the files a Hugging Face Space needs into a folder (e.g. a clone of the Space repo).

    python apps/coal-app/space/build_space.py PATH/TO/SPACE_CLONE

Copies: app.py, README.md, requirements.txt (from this folder), the engine (src/coal at the repository root), the app
(apps/coal-app/src/coal_app) and the repository's examples/. In the Space folder both packages sit side by side under
src/, so app.py puts that one folder on the path. Does not deploy anything; see DEPLOY.md.
"""
from __future__ import annotations

import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))     # apps/coal-app/space
APP = os.path.dirname(HERE)                           # apps/coal-app
ROOT = os.path.dirname(os.path.dirname(APP))          # the repository root
IGNORE = shutil.ignore_patterns("__pycache__", "*.pyc", "*.egg-info")
# (source, destination in the Space folder)
TREES = ((os.path.join(ROOT, "src", "coal"), "src/coal"),
         (os.path.join(APP, "src", "coal_app"), "src/coal_app"),
         (os.path.join(ROOT, "examples"), "examples"))


def build(dest: str) -> list[str]:
    os.makedirs(dest, exist_ok=True)
    for f in ("app.py", "README.md", "requirements.txt"):
        shutil.copy2(os.path.join(HERE, f), os.path.join(dest, f))
    for source, d in TREES:
        target = os.path.join(dest, d)
        if os.path.exists(target):
            shutil.rmtree(target)
        shutil.copytree(source, target, ignore=IGNORE)
    shutil.copy2(os.path.join(ROOT, "LICENSE"), os.path.join(dest, "LICENSE"))
    return sorted(os.listdir(dest))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    print("wrote:", ", ".join(build(sys.argv[1])))
