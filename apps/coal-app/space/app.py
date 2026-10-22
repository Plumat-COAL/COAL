"""Hugging Face Space entry point (copied to the Space root by space/build_space.py)."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "src")
sys.path.insert(0, SRC)
# the GUI runs the engine as `python -m coal` in a subprocess, which must find it too
os.environ["PYTHONPATH"] = SRC + os.pathsep + os.environ.get("PYTHONPATH", "")

from coal_app.space_app import main  # noqa: E402

main()
