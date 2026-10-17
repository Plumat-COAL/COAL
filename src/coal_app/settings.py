"""User settings of the COAL app: the defaults a user sets once (genetic code, output format, options, preview size).

CHANGED(GUI-9): new file. See overnight/integration/CHANGES.md and docs/ui-principles.md.

Where they live:
- running locally (``python -m coal_app`` or the desktop app): one small JSON file in the user's config folder,
  ``%APPDATA%\\COAL\\settings.json`` on Windows, ``~/Library/Application Support/COAL/settings.json`` on macOS,
  ``$XDG_CONFIG_HOME/coal/settings.json`` (default ``~/.config/coal/settings.json``) elsewhere; ``COAL_SETTINGS_DIR``
  overrides the folder. A file, not the browser's storage, because the desktop app opens a private web view on a
  new port each time, so anything the page keeps in localStorage is gone at the next start.
- on the hosted Space (``COAL_HOSTED=1``): never a file on the server (it would be shared by every visitor);
  the page keeps them in that visitor's browser (``gr.BrowserState``, see ui.py).

Reads are tolerant (a missing, empty, half-written or hand-edited file gives the defaults for whatever is wrong) and
writes are atomic (a temporary file in the same folder, flushed, then ``os.replace``), so a crash or two windows
saving at once can never leave a broken file behind.
"""
from __future__ import annotations

import json
import os
import sys
import tempfile

from . import core

VERSION = 1

# preview size: (sequences drawn, columns drawn, most cells drawn). "large" is the old fixed cap (viewer.MAX_*).
PREVIEW_SIZES: dict[str, tuple[int, int, int]] = {
    "small": (40, 400, 16_000),
    "standard": (100, 1000, 60_000),
    "large": (200, 2000, 200_000),
}
# CHANGED(GUI-13): shorter; the radio's label says what the numbers are ("sequences × columns")
PREVIEW_LABELS = {
    "small": "Small: 40 × 400",
    "standard": "Standard: 100 × 1,000",
    "large": "Large: 200 × 2,000 (slower)",
}

DEFAULTS: dict = {
    "version": VERSION,
    "codontable": 1,
    "output_format": "fasta",
    "nogap": False,
    "nomismatch": False,
    "strict": False,
    "preview": "standard",
}
KEYS = [k for k in DEFAULTS if k != "version"]


def settings_dir() -> str:
    d = os.environ.get("COAL_SETTINGS_DIR")
    if d:
        return d
    if sys.platform == "win32":
        base = os.environ.get("APPDATA") or os.path.join(os.path.expanduser("~"), "AppData", "Roaming")
        return os.path.join(base, "COAL")
    if sys.platform == "darwin":
        return os.path.join(os.path.expanduser("~"), "Library", "Application Support", "COAL")
    base = os.environ.get("XDG_CONFIG_HOME") or os.path.join(os.path.expanduser("~"), ".config")
    return os.path.join(base, "coal")


def settings_path() -> str:
    return os.path.join(settings_dir(), "settings.json")


def stored_in_browser() -> bool:
    """True on the hosted Space: settings live in each visitor's browser, never in a file on the server."""
    return core.hosted()


def clean(data, formats=None) -> dict:
    """A complete, valid settings dict from anything (unknown keys dropped, bad values replaced by defaults).
    formats: the output formats this engine can write (core.detect_features().output_formats)."""
    out = dict(DEFAULTS)
    if not isinstance(data, dict):
        return out
    try:
        code = int(data.get("codontable", DEFAULTS["codontable"]))
        if code in core.GENETIC_CODES:
            out["codontable"] = code
    except (TypeError, ValueError):
        pass
    fmt = data.get("output_format")
    if isinstance(fmt, str) and fmt in (formats or core.OUTPUT_FORMATS):
        out["output_format"] = fmt
    for k in ("nogap", "nomismatch", "strict"):
        if isinstance(data.get(k), bool):
            out[k] = data[k]
    if data.get("preview") in PREVIEW_SIZES:
        out["preview"] = data["preview"]
    return out


def load(formats=None) -> tuple[dict, str]:
    """(settings, note). note is "" or one plain sentence about a file that could not be read."""
    if stored_in_browser():
        return clean(None, formats), ""
    path = settings_path()
    try:
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
    except FileNotFoundError:
        return clean(None, formats), ""
    except OSError as e:
        return clean(None, formats), f"Could not read your saved settings ({e.strerror or e}); using the defaults."
    try:
        data = json.loads(text) if text.strip() else {}
    except ValueError:
        return clean(None, formats), "Your saved settings file was damaged, so the defaults are used until you change one."
    return clean(data, formats), ""


def save(values: dict, formats=None) -> tuple[dict, str]:
    """Write settings atomically. Returns (the cleaned settings, "" or a sentence saying why they were not saved).
    On the hosted Space nothing is written (the page keeps them in the browser)."""
    data = clean(values, formats)
    if stored_in_browser():
        return data, ""
    d = settings_dir()
    tmp = None
    try:
        os.makedirs(d, exist_ok=True)
        fd, tmp = tempfile.mkstemp(prefix=".settings-", suffix=".tmp", dir=d)
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=2, sort_keys=True)
            fh.write("\n")
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, settings_path())
        tmp = None
    except OSError as e:
        return data, f"Could not save your settings ({e.strerror or e}); they apply until you close COAL."
    finally:
        if tmp and os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass
    return data, ""


def reset(formats=None) -> tuple[dict, str]:
    return save(dict(DEFAULTS), formats)


def preview_limits(name) -> tuple[int, int, int]:
    return PREVIEW_SIZES.get(name if name in PREVIEW_SIZES else DEFAULTS["preview"])
