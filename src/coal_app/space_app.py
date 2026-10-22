"""Hugging Face Spaces entry: the COAL GUI + MCP server, listening on 0.0.0.0:7860.

Used by space/app.py. Locally you can try the same thing with ``python -m coal_app.space_app``.
"""
from __future__ import annotations

import os

os.environ.setdefault("GRADIO_ANALYTICS_ENABLED", "False")
# COAL_HOSTED=1: hosted intro text, tools limited to a sandbox.
os.environ["COAL_HOSTED"] = "1"


def main() -> None:
    from . import ui
    port = int(os.environ.get("PORT", os.environ.get("GRADIO_SERVER_PORT", "7860")))
    # mcp=False on the public Space. A remote agent cannot hand the
    # Space its own files (the tools take server paths), so MCP there could only reach the examples; the
    # tools' path sandbox stays on regardless, because Gradio's HTTP API also exposes them.
    ui.launch(host="0.0.0.0", port=port, mcp=False, demo_mode=False)


if __name__ == "__main__":
    main()
